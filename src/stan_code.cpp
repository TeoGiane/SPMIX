#define STRICT_R_HEADERS
#include <Rcpp.h>
#include <stan/math.hpp>
#include <Eigen/Dense>
#include <unsupported/Eigen/KroneckerProduct>
#include <string>
#include <memory>

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
    SamplerParams & obj = *pt;
    Rcpp::Rcout << "Ci siamo" << std::endl;
    //Rcpp::Rcout << obj->DebugString() << std::endl;
    return;
}
