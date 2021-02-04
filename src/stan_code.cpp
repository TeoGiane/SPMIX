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

#include <google/protobuf/text_format.h>

#include "PolyaGamma.h"
#include "sampler_params.pb.h"
#include "univariate_mixture_state.pb.h"
#include "optimization_options.pb.h"
#include "recordio.h"
#include "utils.h"
#include "newton_method.h"
#include "gradient_ascent.h"
#include "functors.h"
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


/*//' Serialization testing
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
}*/


/*//' Translation from serialization testing
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
}*/


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
    if (not(options.is("Message") and Rcpp::as<std::string>(options.slot("type")) == "OptimOptions")) {
        throw std::runtime_error("Input 'options' is not of type Message::OptimOptions.");
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
    OptimOptions options_cp;
    Rcpp::XPtr<OptimOptions>(Rcpp::as<Rcpp::XPtr<OptimOptions>>(options.slot("pointer")))
    	->SerializeToString(&tmp); options_cp.ParseFromString(tmp);

    //function::spmixLogLikelihood fun(data, W, params_cp, state_cp);
    function::test_function fun;
    optimization::NewtonMethod<decltype(fun)> solver(fun, options_cp);

    // Initialize and executing Newton Method
    Eigen::VectorXd x0 = fun.init();
    Rcpp::Rcout << "x0: " << x0.transpose() << std::endl;
    Rcpp::Rcout << "Solving..." << std::endl;
    auto start = std::chrono::high_resolution_clock::now();
    solver.solve(x0);
    auto end = std::chrono::high_resolution_clock::now();
    double duration = std::chrono::duration<double>(end - start).count();
    Rcpp::Rcout << "Duration: " << duration << std::endl;
    optimization::OptimState currstate = solver.get_state();
    Rcpp::Rcout << "minimizer: " << currstate.current_maximizer.transpose() << std::endl;
    Rcpp::Rcout << "||grad_f(x)||: " << currstate.current_gradient.norm() << std::endl;
    Rcpp::Rcout << "||hess_f(x)||: " << currstate.current_hessian.norm() << std::endl;

	return;
}

//' Test to evaluate times and correctness of gradient ascent method for optimization on a test function
//' @export
// [[Rcpp::export]]
void grad_ascent_test(const Rcpp::S4 & state, const std::vector<std::vector<double>> & data,
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
    if (not(options.is("Message") and Rcpp::as<std::string>(options.slot("type")) == "OptimOptions")) {
        throw std::runtime_error("Input 'options' is not of type Message::OptimOptions.");
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
    OptimOptions options_cp;
    Rcpp::XPtr<OptimOptions>(Rcpp::as<Rcpp::XPtr<OptimOptions>>(options.slot("pointer")))
        ->SerializeToString(&tmp); options_cp.ParseFromString(tmp);

    // Instanciating functor and solver
    function::test_function fun;
    //function::spmixLogLikelihood fun(data, W, params_cp, state_cp);
    optimization::GradientAscent<decltype(fun)> solver(fun, options_cp);

    // Initialize and executing Gradient Ascent Method
    Eigen::VectorXd x0 = fun.init();
    Rcpp::Rcout << "x0: " << x0.transpose() << std::endl;
    Rcpp::Rcout << "Solving..." << std::endl;
    auto start = std::chrono::high_resolution_clock::now();
    solver.solve(x0);
    auto end = std::chrono::high_resolution_clock::now();
    double duration = std::chrono::duration<double>(end - start).count();
    Rcpp::Rcout << "Duration: " << duration << std::endl;
    optimization::OptimState currstate = solver.get_state();
    Rcpp::Rcout << "minimizer: " << currstate.current_maximizer.transpose() << std::endl;
    Rcpp::Rcout << "||grad_f(x)||: " << currstate.current_gradient_norm << std::endl;

	return;
}

//' Test for loglikelohood with m_tilde contribution
//' @export
// [[Rcpp::export]]
void test_spmixloglikelihood(const std::vector<std::vector<double>> & data,
                             const Eigen::MatrixXd & W, const std::string & params_str,
                             const std::string & state_str) {

    //Omettiamo i controlli che tanto poi si cancella tutto
    // Create a deep-copy of the messages with the workaround
    // std::string tmp;

    // Messages parsing
    SamplerParams params; google::protobuf::TextFormat::ParseFromString(params_str, &params);
    UnivariateState state; google::protobuf::TextFormat::ParseFromString(state_str, &state);

    // Build spmix objects with mtilde
	int numGroups = data.size();
	int numComponents = params.num_components();
    double rho = 0.95;
    Eigen::VectorXd means(numComponents);
    Eigen::VectorXd sqrt_std_devs(numComponents);
    for (int i = 0; i < numComponents; ++i) {
        means(i) = state.atoms()[i].mean();
        sqrt_std_devs(i) = std::sqrt(state.atoms()[i].stdev());
    }
    Rcpp::Rcout << "means: " << means.transpose() << "\n"
    << "sqrt_std_devs: " << sqrt_std_devs.transpose() << std::endl;

    Eigen::MatrixXd postNormalGammaParams(numComponents,4);
    for (int i = 0; i < numComponents; ++i){
    	postNormalGammaParams.row(i) << params.p0_params().mu0(),params.p0_params().a(),
    								    params.p0_params().b(), params.p0_params().lam_();
    }
    Rcpp::Rcout << "postNormalGammaParams:\n" << postNormalGammaParams << std::endl;

    Eigen::MatrixXd transformed_weights = Eigen::MatrixXd::Zero(numGroups, numComponents);
    Eigen::MatrixXd weights(numGroups, numComponents);
    for (int i = 0; i < numGroups; ++i) {
        Eigen::VectorXd tmp(numComponents);
        for (int j = 0; j < numComponents; ++j) {
            tmp(j) = state.groupparams()[i].weights()[j];
            weights(i,j) = state.groupparams()[i].weights()[j];
        }
        transformed_weights.row(i) = utils::Alr(tmp, true);
    }
    Rcpp::Rcout << "weights:\n" << weights << "\n"
    << "transformed_weights:\n" << transformed_weights << std::endl;

    Eigen::MatrixXd Sigma = 0.85 * Eigen::MatrixXd::Identity(state.sigma().rows(), state.sigma().cols());
    /*for (int i = 0; i < state.sigma().rows(); ++i) {
        for (int j = 0; j < state.sigma().cols(); ++j) {
            Sigma(i,j) = state.sigma().data()[i*state.sigma().rows()+j];
        }
    }*/
    Rcpp::Rcout << "Sigma:\n" << Sigma << std::endl << std::endl;

    // Auxiliary variables to instanciate various tests
    int to_drop = 2;
    Eigen::MatrixXd trans_weights = transformed_weights.block(0,0,numGroups,numComponents-1);

    // Instanciating functors
    Rcpp::Rcout << "EVALUATING FUNCTION AND TESTING - EXPANSION:" << std::endl;
    function::spmixLogLikelihood fun_extended(data,W,params,numGroups,numComponents,rho,means,sqrt_std_devs,
    										  trans_weights,Sigma);

    Eigen::VectorXd x(numGroups+2); x << 1,2,3,4,5,6,2.5,1;
    Rcpp::Rcout << "x: " << x.transpose() << std::endl;
    double res = fun_extended(x);
    Rcpp::Rcout << "fun_extended(x): " << res << std::endl;
    res = fun_extended();
    Rcpp::Rcout << "fun_extended(): " << res << std::endl << std::endl;


    Rcpp::Rcout << "EVALUATING FUNCTION AND TESTING - REDUCTION:" << std::endl;
    function::spmixLogLikelihood fun_reduced(data,W,params,numGroups,numComponents-1,rho,
    										 utils::removeElem(means,to_drop), utils::removeElem(sqrt_std_devs,to_drop),
    										 utils::removeColumn(trans_weights,to_drop),
    										 utils::removeRowColumn(Sigma,to_drop),to_drop);
    x << trans_weights.col(to_drop), means(to_drop), sqrt_std_devs(to_drop);
    Rcpp::Rcout << "x: " << x.transpose() << std::endl;
    res = fun_reduced(x);
    Rcpp::Rcout << "fun_reduced(x): " << res << std::endl;
    res = fun_reduced();
    Rcpp::Rcout << "fun_reduced(): " << res << std::endl;

    return;
}

/*
//' Test to evaluate the acceptance rate in case of extension move
//' @export
// [[Rcpp::export]]
void IncreaseMove_test(const std::vector<std::vector<double>> & data,
                       const Eigen::MatrixXd & W, const Rcpp::S4 & params,
                       const Rcpp::S4 & state, const Rcpp::S4 & options) {

    // Check S4 class for params
    if (not(params.is("Message") and Rcpp::as<std::string>(params.slot("type")) == "SamplerParams")) {
        throw std::runtime_error("Input 'params' is not of type Message::SamplerParams.");
    }

    // Check S4 class for state
    if (not(state.is("Message") and Rcpp::as<std::string>(state.slot("type")) == "UnivariateState")) {
        throw std::runtime_error("Input 'state' is not of type Message::UnivariateState.");
    }

    // Check S4 class for options
    if (not(options.is("Message") and Rcpp::as<std::string>(options.slot("type")) == "OptimOptions")) {
        throw std::runtime_error("Input 'options' is not of type Message::OptimOptions.");
    }

    // Create a deep-copy of the messages with the workaround
    std::string tmp;

    // Params copy
    SamplerParams params_cp;
    Rcpp::XPtr<SamplerParams>(Rcpp::as<Rcpp::XPtr<SamplerParams>>(params.slot("pointer")))
        ->SerializeToString(&tmp); params_cp.ParseFromString(tmp);

    // State copy
    UnivariateState state_cp;
    Rcpp::XPtr<UnivariateState>(Rcpp::as<Rcpp::XPtr<UnivariateState>>(state.slot("pointer")))
        ->SerializeToString(&tmp); state_cp.ParseFromString(tmp);

    // Options copy
    OptimOptions options_cp;
    Rcpp::XPtr<OptimOptions>(Rcpp::as<Rcpp::XPtr<OptimOptions>>(options.slot("pointer")))
        ->SerializeToString(&tmp); options_cp.ParseFromString(tmp);

    // Help quantities
    int numGroups = data.size();
    int numComponents = state_cp.num_components();
    std::vector<double> means;
    std::vector<double> stddevs;
    for (int i = 0; i < numComponents; ++i) {
        means.emplace_back(state_cp.atoms()[i].mean());
		stddevs.emplace_back(state_cp.atoms()[i].stdev());
    }
    Eigen::MatrixXd transformed_weights = Eigen::MatrixXd::Zero(numGroups, numComponents);
    Eigen::MatrixXd weights = Eigen::MatrixXd::Zero(numGroups, numComponents);
    for (int i = 0; i < numGroups; ++i) {
        Eigen::VectorXd tmp(numComponents);
        for (int j = 0; j < numComponents; ++j) {
            tmp(j) = state_cp.groupparams()[i].weights()[j];
            weights(i,j) = state_cp.groupparams()[i].weights()[j];
        }
        transformed_weights.row(i) = utils::Alr(tmp, true);
    }
    Eigen::MatrixXd Sigma(state_cp.sigma().rows(), state_cp.sigma().cols());
    for (int i = 0; i < state_cp.sigma().rows(); ++i) {
        for (int j = 0; j < state_cp.sigma().cols(); ++j) {
            Sigma(i,j) = state_cp.sigma().data()[i*state_cp.sigma().rows()+j];
        }
    }
    std::mt19937_64 rng{213513435};
    double alpha;

    // Eliciting the approximated optimal proposal parameters
    function::spmixLogLikelihood loglik_extended(data, W, params_cp, state_cp);
    optimization::GradientAscent<decltype(loglik_extended)> solver(loglik_extended, options_cp);
    Eigen::VectorXd x0 = loglik_extended.init();
    solver.solve(x0);
    Eigen::VectorXd optMean = solver.get_state().current_maximizer;
    Eigen::MatrixXd optCov = -solver.get_state().current_hessian.inverse();

    // Simulating from the approximated optimal posterior
    Eigen::VectorXd x = stan::math::multi_normal_rng(optMean, optCov, rng);

    //Computing Acceptance Rate
    alpha = std::exp( loglik_extended(x)+stan::math::poisson_lpmf(numComponents+1, 1)
                     -loglik_extended()-stan::math::poisson_lpmf(numComponents, 1)
                     -stan::math::multi_normal_lpdf(x, optMean, optCov) );
    alpha = std::min(1., alpha);

    // Print acceptance rate
	Rcpp::Rcout << "alpha: " << alpha << std::endl << std::endl;

	// DEBUG
    Rcpp::Rcout << "State BEFORE expansion:\n"
    << "numComponents: " << numComponents << "\n"
    << "means: " << Eigen::Map<Eigen::VectorXd>(means.data(), means.size()).transpose() << "\n"
    << "stddevs: " << Eigen::Map<Eigen::VectorXd>(stddevs.data(), stddevs.size()).transpose() << "\n"
    << "weights:\n" << weights << "\n"
    << "transformed_weights:\n" << transformed_weights << "\n"
    << "Sigma:\n" << Sigma << std::endl << std::endl;

    // Increase state
	++numComponents;
	means.emplace_back(x(numGroups));
	stddevs.emplace_back(x(numGroups+1)*x(numGroups+1));
	transformed_weights.conservativeResize(numGroups, numComponents);
	transformed_weights.col(numComponents-2) = x.head(numGroups);
	transformed_weights.col(numComponents-1) = Eigen::VectorXd::Zero(numGroups);
	weights.resize(numGroups, numComponents);
	for (int i = 0; i < numGroups; ++i)
		weights.row(i) = utils::InvAlr(static_cast<Eigen::VectorXd>(transformed_weights.row(i)), true);
	Sigma.conservativeResize(numComponents-1, numComponents-1);
	Sigma.col(numComponents-2).head(numComponents-2) = Eigen::VectorXd::Zero(numComponents-1);
	Sigma.row(numComponents-2).head(numComponents-2) = Eigen::VectorXd::Zero(numComponents-1);
	Sigma(numComponents-2,numComponents-2) = Sigma(0,0);

	// DEBUG
    Rcpp::Rcout << "State AFTER expansion:\n"
    << "numComponents: " << numComponents << "\n"
    << "means: " << Eigen::Map<Eigen::VectorXd>(means.data(), means.size()).transpose() << "\n"
    << "stddevs: " << Eigen::Map<Eigen::VectorXd>(stddevs.data(), stddevs.size()).transpose() << "\n"
    << "weights:\n" << weights << "\n"
    << "transformed_weights:\n" << transformed_weights << "\n"
    << "Sigma:\n" << Sigma << std::endl << std::endl;

    return;
}

//' Test to evaluate the acceptance rate in case of reduction move
//' @export
// [[Rcpp::export]]
void ReduceMove_test(const std::vector<std::vector<double>> & data,
                     const Eigen::MatrixXd & W, const Rcpp::S4 & params,
                     const Rcpp::S4 & state, const Rcpp::S4 & options) {

    // Check S4 class for params
    if (not(params.is("Message") and Rcpp::as<std::string>(params.slot("type")) == "SamplerParams")) {
        throw std::runtime_error("Input 'params' is not of type Message::SamplerParams.");
    }

    // Check S4 class for state
    if (not(state.is("Message") and Rcpp::as<std::string>(state.slot("type")) == "UnivariateState")) {
        throw std::runtime_error("Input 'state' is not of type Message::UnivariateState.");
    }

    // Check S4 class for options
    if (not(options.is("Message") and Rcpp::as<std::string>(options.slot("type")) == "OptimOptions")) {
        throw std::runtime_error("Input 'options' is not of type Message::OptimOptions.");
    }

    // Create a deep-copy of the messages with the workaround
    std::string tmp;

    // Params copydata
    SamplerParams params_cp;
    Rcpp::XPtr<SamplerParams>(Rcpp::as<Rcpp::XPtr<SamplerParams>>(params.slot("pointer")))
        ->SerializeToString(&tmp); params_cp.ParseFromString(tmp);

    // State copy
    UnivariateState state_cp;
    Rcpp::XPtr<UnivariateState>(Rcpp::as<Rcpp::XPtr<UnivariateState>>(state.slot("pointer")))
        ->SerializeToString(&tmp); state_cp.ParseFromString(tmp);

    // Options copy
    OptimOptions options_cp;
    Rcpp::XPtr<OptimOptions>(Rcpp::as<Rcpp::XPtr<OptimOptions>>(options.slot("pointer")))
        ->SerializeToString(&tmp); options_cp.ParseFromString(tmp);

    // Help quantities
    int numGroups = data.size();
    int numComponents = state_cp.num_components();
    std::mt19937_64 rng{213513435};
    Eigen::MatrixXd transformed_weights = Eigen::MatrixXd::Zero(numGroups, numComponents);
    Eigen::MatrixXd weights = Eigen::MatrixXd::Zero(numGroups, numComponents);
    for (int i = 0; i < numGroups; ++i) {
        Eigen::VectorXd tmp(numComponents);
        for (int j = 0; j < numComponents; ++j) {
            tmp(j) = state_cp.groupparams()[i].weights()[j];
            weights(i,j) = state_cp.groupparams()[i].weights()[j];
        }
        transformed_weights.row(i) = utils::Alr(tmp, true);
    }
    double rho = state_cp.rho();
    std::vector<double> means;
    std::vector<double> stddevs;
    for (int i = 0; i < numComponents; ++i) {
        means.emplace_back(state_cp.atoms()[i].mean());
        stddevs.emplace_back(state_cp.atoms()[i].stdev());
    }
    Eigen::MatrixXd Sigma(state_cp.sigma().rows(), state_cp.sigma().cols());
    for (int i = 0; i < state_cp.sigma().rows(); ++i) {
        for (int j = 0; j < state_cp.sigma().cols(); ++j) {
            Sigma(i,j) = state_cp.sigma().data()[i*state_cp.sigma().rows()+j];
        }
    }

    // Randomly select the component to drop
    int to_drop = stan::math::categorical_rng(Eigen::VectorXd::Constant(numComponents-1, 1./(numComponents-1)), rng)-1;

    // Build the reduced loglikelihood
    // Eigen::MatrixXd trans_weights = transformed_weights.block(0,0,numGroups,numComponents-1);
    Eigen::MatrixXd trans_weights_reduced = utils::removeColumn(transformed_weights.block(0,0,numGroups,numComponents-1),
                                                                to_drop);
    Eigen::Map<Eigen::VectorXd> trans_weights_vect_reduced(trans_weights_reduced.data(),trans_weights_reduced.size());
    Eigen::Map<Eigen::VectorXd> means_map(means.data(), means.size());
    Eigen::VectorXd sqrt_stddevs(numComponents);
    for (int i = 0; i < numComponents; ++i)
        sqrt_stddevs(i) = std::sqrt(stddevs[i]);


    /*Rcpp::Rcout << "to_drop: " << to_drop << "\n"
    << "trans_weights_reduced:\n" << trans_weights_reduced << "\n"
    << "means_reduced: " << utils::removeElem(means_map, to_drop).transpose() << "\n"
    << "sqrt_stddevs_reduced: " << utils::removeElem(sqrt_stddevs, to_drop).transpose() << "\n"
    << "Sigma_reduced:\n" << utils::removeRowColumn(Sigma,to_drop) << std::endl;*/

    /*function::spmixLogLikelihood loglik_reduced(data, W, params_cp, numGroups, numComponents-1, rho,
                                                utils::removeElem(means_map,to_drop), utils::removeElem(sqrt_stddevs,to_drop),
                                                trans_weights_vect_reduced, utils::removeRowColumn(Sigma,to_drop));

    // Eliciting the approximated optimal proposal parameters
    optimization::GradientAscent<decltype(loglik_reduced)> solver(loglik_reduced, options_cp);
    Eigen::VectorXd x0(numGroups+2);
    x0 << transformed_weights.col(to_drop), means_map(to_drop), sqrt_stddevs(to_drop);
    // Rcpp::Rcout << "x0: " << x0.transpose() << std::endl;
    solver.solve(x0);
    Eigen::VectorXd optMean = solver.get_state().current_maximizer;
    Eigen::MatrixXd optCov = -solver.get_state().current_hessian.inverse();

    double alpha{0.};
    if (solver.get_state().current_iteration < options_cp.max_iter() and !solver.get_state().stagnated) {

        // Compute Acceptance rate
        alpha = std::exp( loglik_reduced()+stan::math::poisson_lpmf(numComponents-1, 1)
                         -loglik_reduced(x0)-stan::math::poisson_lpmf(numComponents, 1)
                         +stan::math::multi_normal_lpdf(x0, optMean, optCov) );
    }

    // Accept of Reject the move
    alpha = std::min(1., alpha);

    // Print acceptance rate
    Rcpp::Rcout << "alpha: " << alpha << std::endl << std::endl;

    // DEBUG
    Rcpp::Rcout << "State BEFORE reduction:\n"
    << "numComponents: " << numComponents << "\n"
    << "means: " << means_map.transpose() << "\n"
    << "sqrt_stddevs: " << sqrt_stddevs.transpose() << "\n"
    << "weights:\n" << weights << "\n"
    << "transformed_weights:\n" << transformed_weights << "\n"
    << "Sigma:\n" << Sigma << std::endl << std::endl;

	// Reduce state
    --numComponents;
    means.erase(means.begin()+to_drop);
    stddevs.erase(stddevs.begin()+to_drop);
    transformed_weights = utils::removeColumn(transformed_weights, to_drop);
    weights.resize(numGroups, numComponents);
    for (int i = 0; i < numGroups; ++i)
        weights.row(i) = utils::InvAlr(static_cast<Eigen::VectorXd>(transformed_weights.row(i)), true);
    Sigma = utils::removeRowColumn(Sigma, to_drop);

	// DEBUG
    Rcpp::Rcout << "State AFTER reduction:\n"
    << "numComponents: " << numComponents << "\n"
    << "means: " << Eigen::Map<Eigen::VectorXd>(means.data(), means.size()).transpose() << "\n"
    << "sqrt_stddevs: " << Eigen::Map<Eigen::VectorXd>(stddevs.data(), stddevs.size()).transpose() << "\n"
    << "weights:\n" << weights << "\n"
    << "transformed_weights:\n" << transformed_weights << "\n"
    << "Sigma:\n" << Sigma << std::endl << std::endl;

    return;
}*/

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
    if (not(options.is("Message") and Rcpp::as<std::string>(options.slot("type")) == "OptimOptions")) {
        throw std::runtime_error("Input 'options' is not of type Message::OptimOptions.");
    }

    // Create a deep-copy of the messages with the workaround
    std::string tmp;

    // Params copy
    SamplerParams params_cp;
    Rcpp::XPtr<SamplerParams>(Rcpp::as<Rcpp::XPtr<SamplerParams>>(params.slot("pointer")))
        ->SerializeToString(&tmp); params_cp.ParseFromString(tmp);

    // Options copy
    OptimOptions options_cp;
    Rcpp::XPtr<OptimOptions>(Rcpp::as<Rcpp::XPtr<OptimOptions>>(options.slot("pointer")))
        ->SerializeToString(&tmp); options_cp.ParseFromString(tmp);

    // Tests
    SpatialMixtureRJSampler sampler(params_cp, data, W, options_cp);
    sampler.init();
    
    //sampler.sampleSigma();
    //sampler.sampleW();
    //sampler.betweenModelMove();
    
    Rcpp::Rcout << "Testing sample()..." << std::endl;
    for (int i = 0; i < 20; ++i)
    	sampler.sample();
    Rcpp::Rcout << "OK HERE!" << std::endl;

    return;
}

//' Test for poisson_lpmf output in stan
//' @export
// [[Rcpp::export]]
void poisson_lpmf(size_t seed) {
    std::mt19937_64 rng{seed};
    int N = 100;
    for (int i = 1; i < N; ++i)
        Rcpp::Rcout << "pmf(" << i-1 << ")/pmf(" << i << "): "
        << std::exp(stan::math::poisson_lpmf(i-1,1)-stan::math::poisson_lpmf(i,1)) << std::endl;
    Rcpp::Rcout << std::endl;

    for (int i = 1; i < N; ++i)
        Rcpp::Rcout << "pmf(" << i-1 << ")/pmf(" << i << "): "
        << std::exp(stan::math::poisson_lpmf(i-1,1))/std::exp(stan::math::poisson_lpmf(i,1)) << std::endl;
    Rcpp::Rcout << std::endl;

    for (int i = 1; i < N; ++i)
        Rcpp::Rcout << "pmf(" << i+1 << ")/pmf(" << i << "): "
        << std::exp(stan::math::poisson_lpmf(i+1,1)-stan::math::poisson_lpmf(i,1)) << std::endl;
    Rcpp::Rcout << std::endl;

    for (int i = 1; i < N; ++i)
        Rcpp::Rcout << "pmf(" << i+1 << ")/pmf(" << i << "): "
        << std::exp(stan::math::poisson_lpmf(i+1,1))/std::exp(stan::math::poisson_lpmf(i,1)) << std::endl;
    Rcpp::Rcout << std::endl;

    return;
}

//' Test to swap two columns of a matrix with eigen
//' @export
// [[Rcpp::export]]
Eigen::MatrixXd swapCols(Eigen::MatrixXd mat, unsigned int index1, unsigned int index2){
    Eigen::MatrixXd out = mat;
    out.col(index1).swap(out.col(index2));
    return out;
}
