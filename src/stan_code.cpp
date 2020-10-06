// [[Rcpp::depends(BH)]]
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppParallel)]]
// [[Rcpp::depends(StanHeaders)]]
#define STRICT_R_HEADERS
#include <stan/math/fwd/mat.hpp>
#include <stan/math/mix/mat.hpp>
#include <stan/math.hpp>
#include <Rcpp.h>
#include <RcppEigen.h>

#include <Eigen/Dense>
#include <unsupported/Eigen/KroneckerProduct>
#include <string>
#include <memory>
#include <functional>

#include "PolyaGamma.h"
#include "sampler_params.pb.h"
#include "univariate_mixture_state.pb.h"
#include "recordio.h"
#include "utils.h"

//' Simple test with stan/math C++ library
//'
//' Simply computes logN(1|2,3)
//' @export
// [[Rcpp::export]]
void stan_HelloWorld() {
	Rcpp::Rcout << "log normal(1 | 2, 3) = " << stan::math::normal_log(1, 2, 3) << std::endl;
	PolyaGamma pg;
	Rcpp::Rcout << "A PolyaGamma object has been instanciated!" << std::endl;
	Rcpp::Rcout << "Testing for KroneckerProduct..." << std::endl;
	Eigen::MatrixXd A(3,3); A << 1, 2, 3, 4, 5, 6, 7, 8, 9;
	Eigen::MatrixXd B(2,2); B << 4, 3, 2, 1;
	Rcpp::Rcout << "A:\n" << A << std::endl << std::endl;
	Rcpp::Rcout << "B:\n" << B << std::endl << std::endl;
	Eigen::MatrixXd C = Eigen::kroneckerProduct(A,B).eval();
	Rcpp::Rcout << "A x B:\n" << C << std::endl << std::endl;
  return;
}


//' Serialization testing
//' @export
// [[Rcpp::export]]
std::vector<Rcpp::RawVector> fromProto_tostring() {

	Rcpp::Rcout << "Reading params..." << std::endl;
	SamplerParams params;
    params = loadTextProto<SamplerParams>("/home/m_gianella/Documents/R-Packages/SPMIX_input/sampler_params.asciipb");
    Rcpp::Rcout << params.DebugString() << std::endl;

    Rcpp::Rcout << "Creating states..." << std::endl;
    std::vector<UnivariateMixtureState> states;
    UnivariateMixtureAtom* atom;
    for (int i = 0; i < 5; i++) {
        UnivariateMixtureState curr;
        curr.set_num_components(i + 3);
        for (int j = 0; j < i+3; j++) {
            atom = curr.add_atoms();
            atom->set_mean(-1.0 * i);
            atom->set_stdev(0.5 * i);
        }
        states.push_back(curr);
    }

    for (int i = 0; i < states.size(); ++i) {
    	Rcpp::Rcout << "Printing state " << i+1 << std::endl;
    	Rcpp::Rcout << states[i].DebugString() << std::endl << std::endl;
    }

    Rcpp::Rcout << "Serializing messages... ";
    std::vector<std::string> tmp;
    for (int i = 0; i < states.size(); ++i) {
    	std::string s;
		states[i].SerializeToString(&s);
		tmp.push_back(s);
    }

    Rcpp::Rcout << "done!" << std::endl;

    std::vector<Rcpp::RawVector> out;
    for (std::string elem : tmp)
      out.push_back(utils::str2raw(elem));

    return out;
}


//' Translation from serialization testing
//' @export
// [[Rcpp::export]]
void readingStates(std::vector<Rcpp::RawVector> raw_vect){

	Rcpp::Rcout << "Restoring original messages from raw vectors:" << std::endl;
  size_t i = 1;
	for (Rcpp::RawVector elem : raw_vect) {
		Rcpp::Rcout << "Message " << i++ << " is:" << std::endl;
		UnivariateMixtureState curr;
		curr.ParseFromString(utils::raw2str(elem));
		Rcpp::Rcout << curr.DebugString() << std::endl << std::endl;
	}
	return;
}

//' Simple test to see if messages in R can be passed to a C++ function
//' @export
// [[Rcpp::export]]
void messageFromR(Rcpp::S4 params) {
    Rcpp::Rcout << std::boolalpha << "Has pointer: " << params.hasSlot("pointer") << std::endl;
    Rcpp::Rcout << "S4 Class of type Message? " << params.is("Message") << std::endl;
    Rcpp::XPtr<SamplerParams> pt = params.slot("pointer");
    std::string tmp; SamplerParams obj; pt->SerializeToString(&tmp);
    obj.ParseFromString(tmp); Rcpp::Rcout << obj.DebugString() << std::endl;
    //Rcpp::Rcout << (*(pt.get())).DebugString() << std::endl;
    return;
}


/*template<typename T>
class function_test {
	using FunctionType = std::function<T(const Eigen::Matrix<T,Eigen::Dynamic,1> &)>;
	template<typename T> FunctionType<T> function;
public:
	function_test(const FunctionType<T> _function)
};*/


class function {
private:
	Eigen::VectorXd data;
public:
	function(const Eigen::VectorXd & _data): data(_data) {};
	template<typename T>
    T operator()(const Eigen::Matrix<T, Eigen::Dynamic, 1> & data_toadd) const {
		Eigen::Matrix<T, Eigen::Dynamic, 1> data_prom(data);
		Eigen::Matrix<T, Eigen::Dynamic, 1> data_ext(data_prom.size() + data_toadd.size());
		data_ext << data_prom, data_toadd;
		return stan::math::variance(data_ext);
    }
};

//' Test to use stan math library for automatic differentiation
//' @export
// [[Rcpp::export]]
void hessian_test() {
	using StanVector = Eigen::Matrix<stan::math::var, Eigen::Dynamic, 1>;
	Eigen::VectorXd data(4); data << 1., 2., 3., 4.;
	Rcpp::Rcout << "data' = " << data.transpose() << std::endl;
	function fun(data);
	Eigen::VectorXd input(2); input << 2.5, 0.01;
	Rcpp::Rcout << "input' = " << input.transpose() << std::endl;
	Rcpp::Rcout << "output = " << fun(input) << std::endl;

	double eval;
	Eigen::VectorXd grad_x;
	Eigen::MatrixXd hess_x;

	Rcpp::Rcout << "Automatic differentiation..." << std::endl;
	auto start = std::chrono::high_resolution_clock::now();
	stan::math::hessian(fun, input, eval, grad_x, hess_x);
	auto end = std::chrono::high_resolution_clock::now();
	double duration = std::chrono::duration<double>(end - start).count();
	Rcpp::Rcout << "Duration: " << duration << std::endl;
	Rcpp::Rcout << "grad(fun)(x):\n" << grad_x << std::endl;
	Rcpp::Rcout << "hess(fun)(x):\n" << hess_x << std::endl;
	Rcpp::Rcout << std::endl;
	Rcpp::Rcout << "Finite Differences..." << std::endl;
	start = std::chrono::high_resolution_clock::now();
	stan::math::finite_diff_hessian(fun, input, eval, grad_x, hess_x);
	end = std::chrono::high_resolution_clock::now();
	duration = std::chrono::duration<double>(end - start).count();
	Rcpp::Rcout << "Duration: " << duration << std::endl;
	Rcpp::Rcout << "grad(fun)(x):\n" << grad_x << std::endl;
	Rcpp::Rcout << "hess(fun)(x):\n" << hess_x << std::endl;

	return;
}
