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
#include <google/protobuf/text_format.h>

#include "utils.h"
#include "functors.h"
#include "sampler.h"
#include "sampler_rjmcmc.h"
#include "univariate_mixture_state.pb.h"


//' Additive Log Ratio
//'
//' \loadmathjax This utility computes the additive log-ratio transform of a given vector. Given a generic vector of the simplex
//' \mjseqn{x \in S^{H}}, the transformation is defined as:
//' \mjsdeqn{ \operatorname{alr}(x)_j = \log \left( \frac{x_j}{x_H} \right) \quad \forall j=1,\dots,H-1 }
//'
//' @param x Vector of double in the simplex \mjseqn{S^H}.
//' @return Vector of double in \mjseqn{\mathbb{R}^{H-1}} i.e. \mjseqn{\operatorname{alr}(x)}.
//' @export
// [[Rcpp::export]]
Eigen::VectorXd Alr(Eigen::VectorXd x) {
	return utils::Alr(x, false);
}

/*
&=& \frac{ \operatorname{e}^{x_j} }{ 1 + \sum_{h=1}^{H-1} \operatorname{e}^{x_h} }
*/


//' Inverse Additive Log Ratio
//'
//' \loadmathjax This utility computes the inverse additive log-ratio transform of a given vector. Given a generic vector 
//' \mjseqn{ x \in \mathbb{R}^{H-1} }, the transformation is defined as:
//' \mjsdeqn{ \begin{eqnarray*}
//'	\operatorname{alr}^{-1}(x)_{j} &=& \textstyle\frac{\operatorname{e}^{x_j}}{\sum _{h} \operatorname{e}^{x_h}} \quad \forall j=1,\dots,H-1 \cr
//' \operatorname{alr}^{-1}(x)_H &=& 1 - \textstyle\sum _{h} \operatorname{alr}^{-1}(x)_h
//' 		 \end{eqnarray*} }
//'
//' @param x Vector of double in \mjseqn{ \mathbb{R}^{H-1} }.
//' @return Vector of double in the simplex \mjseqn{S^H} i.e. \mjseqn{\operatorname{alr}^{-1}(x)}.
//' @export
// [[Rcpp::export]]
Eigen::VectorXd InvAlr(Eigen::VectorXd x) {
  return utils::InvAlr(x, false);
}

/* Spatial Sampler execution routine (no RJMCMC, w/ or w/o covariates,data,W and params as strings or proper data types from R)*/
// [[Rcpp::export]]
std::vector<Rcpp::RawVector> runSpatialSampler(int burnin, int niter, int thin, const std::vector<std::vector<double>> & data,
    										   const Eigen::MatrixXd & W, std::string params_str,
    										   const std::vector<Eigen::MatrixXd> & covariates, bool display_progress) {

	// Messages Parsing
	SamplerParams params; google::protobuf::TextFormat::ParseFromString(params_str, &params);

    // Initializarion
    SpatialMixtureSampler spSampler(params, data, W, covariates);
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
    											 const Eigen::MatrixXd & W, std::string params_str,
    											 const std::vector<Eigen::MatrixXd> & covariates,
    											 const std::string & options_str, bool boundary_detection, bool display_progress) {

	// Messages Parsing
	SamplerParams params; google::protobuf::TextFormat::ParseFromString(params_str, &params);
	OptimOptions options; google::protobuf::TextFormat::ParseFromString(options_str, &options);

	// Initializarion
	SpatialMixtureRJSampler spSampler(params, data, W, options, covariates, boundary_detection);
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
//' \loadmathjax This function simply reads the proximity matrix \mjseqn{G} of the Spatial Mixture Model from a \code{.csv} file.
//' This file should not have columns or row headers and it must be written in the visual form of a matrix 
//' only composed by either \mjseqn{0} or \mjseqn{1}. As assumption, the diagonal of this matrix should be \mjseqn{0}.
//' @param filename A string identifying the path to a \code{.csv} file from which the matrix will be read.
//' @return The proximity matrix as a usual \code{R matrix} object.
//' @export
// [[Rcpp::export]]
Eigen::MatrixXd ReadMatrixFromCSV(std::string filename) {
    return utils::readMatrixFromCSV(filename);
}

//' Import Data from File
//'
//' \loadmathjax This utility reads the input data for the sampler from a \code{.csv} file. This file needs to be provided 
//' with column headers. Moreover, it should be organized in two columns for the correct parsing:
//' \itemize{
//' \item{\strong{group}, an integer (from \mjseqn{0} to \mjseqn{I-1}) describing the area to which the data belongs to;}
//' \item{\strong{data}, the actual data, which will be assigned to area indicated by the "group" column.}
//' }
//' @return A list of dimension \mjseqn{I}, in which the \mjseqn{i}-th element is a vector containing all data that 
//' have been assigned to the \mjseqn{i}-th location.
//' @export
// [[Rcpp::export]]
std::vector<std::vector<double>> ReadDataFromCSV(std::string filename) {
    return utils::readDataFromCSV(filename);
}

#endif // SPMIX_EXPORTS
