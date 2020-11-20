/* This source file contains all utilities that are exported to R through Rcpp. Any other function
 * should be put here in order to preserve coherence in the code
 */

#ifndef SPMIX_EXPORTS
#define SPMIX_EXPORTS

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

#include <deque>
#include <string>
#include <vector>
#include <Eigen/Dense>
#include <progress.hpp>
#include <progress_bar.hpp>
#include <utility>
//#include <google/protobuf/text_format.h>
#include <exception>

#include "utils.h"
#include "functors.h"
#include "sampler.h"
#include "sampler_rjmcmc.h"
#include "univariate_mixture_state.pb.h"

/*TODO: finire documentazione*/

//' Additive log ratio
//'
//' This utility computes the additive log ratio transorm of a given vector. This transformation is
//' defined as: (METTI DEF)
//' @param x Vector of double in \ifelse{html}{\out{R<sup>H</sup>}}{\eqn{R^{H}}}.
//' @param pad_zero Boolean, default to FALSE.
//' @return Vector of double in \ifelse{html}{\out{R<sup>H-1</sup>}}{\eqn{R^{H-1}}} i.e. alr(x).
//' @export
// [[Rcpp::export]]
Eigen::VectorXd alr(Eigen::VectorXd x, bool pad_zero = false) {
	return utils::Alr(x, pad_zero);
}

/*TODO: finire documentazione*/

//' Inverse additive log ratio
//'
//' This utility computes the inverse additive log ratio transorm of a given vector. This transofrmation is
//' defined as: (METTI DEF)
//' @param x Vector of double in \ifelse{html}{\out{R<sup>H-1</sup>}}{\eqn{R^{H-1}}}.
//' @param padded_zero Boolean, default to FALSE.
//' @return Vector of double in \ifelse{html}{\out{R<sup>H</sup>}}{\eqn{R^{H}}} i.e.
//' \ifelse{html}{\out{alr <sup>-1</sup>(x)}}{\eqn{alr^{-1}(x)}}.
//' @export
// [[Rcpp::export]]
Eigen::VectorXd inv_alr(Eigen::VectorXd x, bool padded_zero = false) {
  return utils::InvAlr(x, padded_zero);
}

/* Spatial Sampler execution routine (no RJMCMC,w/ or w/o covariates,data,W and params as strings or proper data types from R)*/
// [[Rcpp::export]]
std::vector<Rcpp::RawVector> runSpatialSampler(int burnin, int niter, int thin, const std::vector<std::vector<double>> & data,
    										   const Eigen::MatrixXd & W, Rcpp::S4 params,
    										   const std::vector<Eigen::MatrixXd> & covariates, bool display_progress) {

	// Deep copy of messages
	std::string tmp;

	// Params copy
    SamplerParams params_cp;
    Rcpp::XPtr<SamplerParams>(Rcpp::as<Rcpp::XPtr<SamplerParams>>(params.slot("pointer")))
    	->SerializeToString(&tmp); params_cp.ParseFromString(tmp);

    // Initializarion
    SpatialMixtureSampler spSampler(params_cp, data, W, covariates);
    spSampler.init();

    // Initialize output container
    std::vector<Rcpp::RawVector> out;

    // Sampling
    auto start = std::chrono::high_resolution_clock::now();
    REprintf("SPMIX Sampler: Burn-in\n");
    Progress p_burn(burnin, display_progress);
    for (int i=0; i < burnin; i++) {
        spSampler.sample();
        p_burn.increment();
    }
    p_burn.cleanup();
    Rcpp::Rcout << std::endl;


    REprintf("SPMIX Sampler: Running\n");
    Progress p_run(niter, display_progress);
    for (int i=0; i < niter; i++) {
        spSampler.sample();
        if ((i+1) % thin == 0) {
            std::string s;
            spSampler.getStateAsProto().SerializeToString(&s);
            out.push_back(utils::str2raw(s));
        }
        p_run.increment();
    }
    Rcpp::Rcout << std::endl;
    auto end = std::chrono::high_resolution_clock::now();

    double duration = std::chrono::duration<double>(end - start).count();
    Rcpp::Rcout << "Duration: " << duration << std::endl;
    return out;
}

/* Spatial Sampler execution routine (RJMCMC,w/ or w/o covariates,data,W and params as strings or proper data types from R)*/
// [[Rcpp::export]]
std::vector<Rcpp::RawVector> runSpatialRJSampler(int burnin, int niter, int thin, const std::vector<std::vector<double>> & data,
    											 const Eigen::MatrixXd & W, Rcpp::S4 params,
    											 const std::vector<Eigen::MatrixXd> & covariates,
    											 const Rcpp::S4 & options, bool display_progress) {

	// Deep copy of messages
	std::string tmp;

	// Params copy
    SamplerParams params_cp;
    Rcpp::XPtr<SamplerParams>(Rcpp::as<Rcpp::XPtr<SamplerParams>>(params.slot("pointer")))
    	->SerializeToString(&tmp); params_cp.ParseFromString(tmp);

    // Options copy
    OptimOptions options_cp;
    Rcpp::XPtr<OptimOptions>(Rcpp::as<Rcpp::XPtr<OptimOptions>>(options.slot("pointer")))
    	->SerializeToString(&tmp); options_cp.ParseFromString(tmp);

    // Initializarion
	SpatialMixtureRJSampler spSampler(params_cp, data, W, options_cp, covariates);
	spSampler.init();

    // Initialize output container
    std::vector<Rcpp::RawVector> out;

    // Sampling
    auto start = std::chrono::high_resolution_clock::now();
    REprintf("SPMIX RJ Sampler: Burn-in\n");
    Progress p_burn(burnin, display_progress);
    for (int i=0; i < burnin; i++) {
        //Rcpp::Rcout << "Burn-in iter: " << i+1 << std::endl;
        spSampler.sample();
        p_burn.increment();
    }
    p_burn.cleanup();
    Rcpp::Rcout << std::endl;


    REprintf("SPMIX RJ Sampler: Running\n");
    Progress p_run(niter, display_progress);
    for (int i=0; i < niter; i++) {
        //Rcpp::Rcout << "Running iter: " << i+1 << std::endl;
        spSampler.sample();
        if ((i+1) % thin == 0) {
            std::string s;
            spSampler.getStateAsProto().SerializeToString(&s);
            out.push_back(utils::str2raw(s));
        }
        p_run.increment();
    }
    Rcpp::Rcout << std::endl;
    auto end = std::chrono::high_resolution_clock::now();

    double duration = std::chrono::duration<double>(end - start).count();
    Rcpp::Rcout << "Duration: " << duration << std::endl;
    Rcpp::Rcout << "Acceptance rate: " << static_cast<double>(spSampler.get_acceptedMoves())/(burnin+niter)*100 << "%" << std::endl;
    return out;
}

//' Just a check for output format in R when reading a matrix
//' @export
// [[Rcpp::export]]
Eigen::MatrixXd readMatrixFromCSV(std::string filename) {
    return utils::readMatrixFromCSV(filename);
}

//' Just a check for output format in R when reading data
//' @export
// [[Rcpp::export]]
std::vector<std::vector<double>> readDataFromCSV(std::string filename) {
    return utils::readDataFromCSV(filename);
}

/*
//' Loglikelihood of a Spatial Mixture model state
//' @export
// [[Rcpp::export]]
double spmixLogLikelihood(const Rcpp::S4 & state, const std::vector<std::vector<double>> & data,
                          const Eigen::MatrixXd & W, const Rcpp::S4 & params) {

    // Check S4 class for state
    if (not(state.is("Message") and Rcpp::as<std::string>(state.slot("type")) == "UnivariateState")) {
        throw std::runtime_error("Input 'state' is not of type Message::UnivariateState.");
    }

    // Check S4 class for params
    if (not(params.is("Message") and Rcpp::as<std::string>(params.slot("type")) == "SamplerParams")) {
        throw std::runtime_error("Input 'params' is not of type Message::SamplerParams.");
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

    function::spmixLogLikelihood fun(data, W, params_cp, state_cp);
    return fun();
}*/

#endif // SPMIX_EXPORTS
