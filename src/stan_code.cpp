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
#include "newton_options.pb.h"
#include "recordio.h"
#include "utils.h"
#include "newton_opt.h"
#include "mcmc_utils.h"
#include "sampler_rjmcmc.h"

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
    params = loadTextProto<SamplerParams>("/home/m_gianella/Documents/R-Packages/SPMIX/inst/input_files/sampler_params.asciipb");
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


//' Test to evaluate times and correctness of newton method for optimization
//' @export
// [[Rcpp::export]]
void newton_opt_test(const Rcpp::S4 & state, const std::vector<std::vector<double>> & data,
                     const Eigen::MatrixXd & W, const Rcpp::S4 & params, const Rcpp::S4 & options) {

	// Check S4 class for state
    if (not(state.is("Message") and Rcpp::as<std::string>(state.slot("type")) == "UnivariateState")) {
        throw std::runtime_error("Input 'state' is not of type Message::UnivariateState.");
    }

    // Check S4 class for params
    if (not(params.is("Message") and Rcpp::as<std::string>(params.slot("type")) == "SamplerParams")) {
        throw std::runtime_error("Input 'params' is not of type Message::SamplerParams.");
    }

    // Check S4 class for options
    if (not(options.is("Message") and Rcpp::as<std::string>(options.slot("type")) == "NewtonOptions")) {
        throw std::runtime_error("Input 'options' is not of type Message::NewtonOptions.");
    }

    // Create a deep-copy of the messages with the workaround
    std::string tmp;

    // State copy
    UnivariateState state_cp;
    Rcpp::XPtr<UnivariateState>(Rcpp::as<Rcpp::XPtr<UnivariateState>>(state.slot("pointer")))
    	->SerializeToString(&tmp); state_cp.ParseFromString(tmp);

    // Params copy
    SamplerParams params_cp;
    Rcpp::XPtr<SamplerParams>(Rcpp::as<Rcpp::XPtr<SamplerParams>>(params.slot("pointer")))
    	->SerializeToString(&tmp); params_cp.ParseFromString(tmp);

    // Options copy
    NewtonOptions options_cp;
    Rcpp::XPtr<NewtonOptions>(Rcpp::as<Rcpp::XPtr<NewtonOptions>>(options.slot("pointer")))
    	->SerializeToString(&tmp); options_cp.ParseFromString(tmp);

    utils::spmixLogLikelihood fun(data, W, params_cp, state_cp);
    NewtonOpt solver(fun, options_cp);

    // Initialize and executing Newton Method
    Eigen::VectorXd x0 = solver.init();
    Rcpp::Rcout << "x0: " << x0.transpose() << std::endl;
    Rcpp::Rcout << "Solving..." << std::endl;
    auto start = std::chrono::high_resolution_clock::now();
    solver.solve(x0);
    auto end = std::chrono::high_resolution_clock::now();
    double duration = std::chrono::duration<double>(end - start).count();
    Rcpp::Rcout << "Duration: " << duration << std::endl;
    NewtonState currstate = solver.get_state();
    Rcpp::Rcout << "minimizer: " << currstate.current_minimizer.transpose() << std::endl;
    Rcpp::Rcout << "||grad_f(x)||: " << currstate.current_gradient.norm() << std::endl;
    Rcpp::Rcout << "||hess_f(x)||: " << currstate.current_hessian.norm() << std::endl;

	return;
}


//' Test fot the RJSampler
//' @export
// [[Rcpp::export]]
void RJsampler_test(const std::vector<std::vector<double>> & data,
                    const Eigen::MatrixXd & W, const Rcpp::S4 & params, const Rcpp::S4 & options) {

    // Check S4 class for params
    if (not(params.is("Message") and Rcpp::as<std::string>(params.slot("type")) == "SamplerParams")) {
        throw std::runtime_error("Input 'params' is not of type Message::SamplerParams.");
    }

    // Check S4 class for options
    if (not(options.is("Message") and Rcpp::as<std::string>(options.slot("type")) == "NewtonOptions")) {
        throw std::runtime_error("Input 'options' is not of type Message::NewtonOptions.");
    }

    // Create a deep-copy of the messages with the workaround
    std::string tmp;

    // Params copy
    SamplerParams params_cp;
    Rcpp::XPtr<SamplerParams>(Rcpp::as<Rcpp::XPtr<SamplerParams>>(params.slot("pointer")))
        ->SerializeToString(&tmp); params_cp.ParseFromString(tmp);

    // Options copy
    NewtonOptions options_cp;
    Rcpp::XPtr<NewtonOptions>(Rcpp::as<Rcpp::XPtr<NewtonOptions>>(options.slot("pointer")))
        ->SerializeToString(&tmp); options_cp.ParseFromString(tmp);

    // Various tests
    SpatialMixtureRJSampler sampler(params_cp, data, W, options_cp);
    sampler.init();
    sampler.sampleSigma();
    sampler.betweenModelMove();

    return;
}