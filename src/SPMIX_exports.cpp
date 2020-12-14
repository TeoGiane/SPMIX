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
#include <exception>

#include "utils.h"
#include "functors.h"
#include "sampler.h"
#include "sampler_rjmcmc.h"
#include "univariate_mixture_state.pb.h"

/*TODO: finire documentazione*/

//' Additive Log Ratio
//'
//' This utility computes the additive log ratio transorm of a given vector. This transformation is
//' defined as: (METTI DEF)
//' @param x Vector of double in \ifelse{html}{\out{R<sup>H</sup>}}{\eqn{R^{H}}}.
//' @return Vector of double in \ifelse{html}{\out{R<sup>H-1</sup>}}{\eqn{R^{H-1}}} i.e. alr(x).
//' @export
// [[Rcpp::export]]
Eigen::VectorXd alr(Eigen::VectorXd x) {
	return utils::Alr(x, false);
}

/*TODO: finire documentazione*/

//' Inverse Additive Log Ratio
//'
//' This utility computes the inverse additive log ratio transorm of a given vector. This transformation is
//' defined as: (METTI DEF)
//' @param x Vector of double in \ifelse{html}{\out{R<sup>H-1</sup>}}{\eqn{R^{H-1}}}.
//' @return Vector of double in \ifelse{html}{\out{R<sup>H</sup>}}{\eqn{R^{H}}} i.e.
//' \ifelse{html}{\out{alr <sup>-1</sup>(x)}}{\eqn{alr^{-1}(x)}}.
//' @export
// [[Rcpp::export]]
Eigen::VectorXd invAlr(Eigen::VectorXd x) {
  return utils::InvAlr(x, false);
}

/* Spatial Sampler execution routine (no RJMCMC, w/ or w/o covariates,data,W and params as strings or proper data types from R)*/
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

/* Spatial Sampler execution routine (RJMCMC, w/ or w/o covariates,data,W and params as strings or proper data types from R)*/
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
        spSampler.sample();
        p_burn.increment();
    }
    p_burn.cleanup();
    Rcpp::Rcout << std::endl;


    REprintf("SPMIX RJ Sampler: Running\n");
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

//' Import Proximity Matrix from File
//'
//' This function simply reads the proximity matrix G of the Spatial Mixture Model from a \code{.csv} file. This file
//' needs to be provided with no columns or row headers and in the visual form of a matrix only composed by 0 or 1.
//' As assumption, the diagonal of this matrix should be 0.
//' @param filename A string identifying the path to a \code{.csv} file from which the matrix will be read.
//' @return The proximity matrix as a usual \code{R matrix} object.
//' @export
// [[Rcpp::export]]
Eigen::MatrixXd readMatrixFromCSV(std::string filename) {
    return utils::readMatrixFromCSV(filename);
}

//' Import Data from File
//'
//' This utility reads the input data for the sampler from a \code{.csv} file. This file needs to be provided with column headers.
//' The file should be organized in two columns for a correct read:
//' \itemize{
//' \item{\strong{group}, an integer (from 0 to I-1) describing the area to which the datum belongs to;}
//' \item{\strong{data}, the actual datum, which will be assigned to area indicated by "group" column.}
//' }
//' @return A list of dimension I, in which the i-th element is a vector containing all data that have been assigned to location i.
//' @export
// [[Rcpp::export]]
std::vector<std::vector<double>> readDataFromCSV(std::string filename) {
    return utils::readDataFromCSV(filename);
}

#endif // SPMIX_EXPORTS
