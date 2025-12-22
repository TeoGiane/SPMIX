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
// #include <stan/math/fwd/mat.hpp>
// #include <stan/math/mix/mat.hpp>
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
// #include "functors.h"
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
    										   const Eigen::MatrixXd & W, std::string params_filename,
    										   const std::vector<Eigen::MatrixXd> & covariates,
											   bool boundary_detection, bool display_progress, unsigned long seed) {
	
	// Parse Sampler Parameters
	spmix::SamplerParams params;
	std::ifstream params_file(params_filename);
	if (params_file.is_open()) {
		params.ParseFromIstream(&params_file);
		params_file.close();
	} else {
		throw std::runtime_error("Cannot open parameters file: " + params_filename);
	}

	// Initializarion
	SpatialMixtureSampler spSampler(params, data, W, covariates, boundary_detection, seed);
	spSampler.init();

	// Initialize output container
	std::vector<Rcpp::RawVector> out;

	// Sampling
	auto start = std::chrono::high_resolution_clock::now();

	if (burnin > 0) {
		if (display_progress) { REprintf("SPMIX Sampler: Burn-in\n"); }
		Progress p_burn(burnin, display_progress);
		for (int i=0; i < burnin; i++) {
			spSampler.sample();
			p_burn.increment();
		}
		p_burn.cleanup();
		Rcpp::Rcout << std::endl;
	}

	if (niter > 0) {
		if (display_progress) { REprintf("SPMIX Sampler: Running\n"); }
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
	}

	auto end = std::chrono::high_resolution_clock::now();

	double duration = std::chrono::duration<double>(end - start).count();
	if (display_progress) {
		Rcpp::Rcout << "Elapsed Time: " << duration << std::endl << std::endl;
	}
	return out;
}

/* Spatial Sampler execution routine (RJMCMC, w/ or w/o covariates,data,W and params as strings or proper data types from R)*/
// [[Rcpp::export]]
std::vector<Rcpp::RawVector> runSpatialRJSampler(int burnin, int niter, int thin, const std::vector<std::vector<double>> & data,
    											 const Eigen::MatrixXd & W, const std::string & params_filename,
    											 const std::vector<Eigen::MatrixXd> & covariates,
    											 const std::string & options_filename, bool boundary_detection, bool display_progress, unsigned long seed) {
	
	// Parse Sampler Parameters
	spmix::SamplerParams params;
	std::ifstream params_file(params_filename);
	if (params_file.is_open()) {
		params.ParseFromIstream(&params_file);
		params_file.close();
	} else {
		throw std::runtime_error("Cannot open parameters file: " + params_filename);
	}

	// Parse Optimization Options
	spmix::OptimOptions options;
	std::ifstream options_file(options_filename);
	if (options_file.is_open()) {
		options.ParseFromIstream(&options_file);
		options_file.close();
	} else {
		throw std::runtime_error("Cannot open optimization options file: " + options_filename);
	}

	// Initializarion
	SpatialMixtureRJSampler spSampler(params, data, W, options, covariates, boundary_detection, seed);
	spSampler.init();

	// Initialize output container
	std::vector<Rcpp::RawVector> out;

    // Sampling
    auto start = std::chrono::high_resolution_clock::now();

    if (burnin > 0) {
		if (display_progress) { REprintf("SPMIX RJ Sampler: Burn-in\n"); }
		Progress p_burn(burnin, display_progress);
		for (int i=0; i < burnin; i++) {
			spSampler.sample();
			p_burn.increment();
		}
		p_burn.cleanup();
		Rcpp::Rcout << std::endl;
	}

	if (niter > 0) {
		if (display_progress) { REprintf("SPMIX RJ Sampler: Running\n"); }
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
	}

    auto end = std::chrono::high_resolution_clock::now();

	double duration = std::chrono::duration<double>(end - start).count();
	if(display_progress) {
		Rcpp::Rcout << "Elapsed Time: " << duration << std::endl;
		// Rcpp::Rcout << "Acceptance Rate: " << spSampler.computeAcceptanceRate() << std::endl << std::endl;
	}
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

//' Compute the log posterior densities for each group and for each data point
//'
//' \loadmathjax This utility takes as input the deserialized output of the samplers
//' (via \code{\link{DeserializeSPMIXProto}}) and compute the chain of posterior log-likelihood
//' for each data point
//'
//' @param serialized_states A list of `raw` vectore which stores the serialized output of
//' the sampler (either with a fixed or a variable number of components).
//' @param data The data that needs to be fitted by the model. Data are passed as a list of vectors, whose
//' \mjseqn{i}-th element represents the vector of data assigned to the \mjseqn{i}-th location.
//' @param display_progress A bool. If \code{TRUE}, prints the progress of the computation.
//' 
//' @return A \mjseqn{T \times N} matrix, \mjseqn{T} being the number of iterations
//' of the MCMC chain and \mjseqn{N} the total number of data points.
//' Element \mjseqn{t,n} of the matrix is the log-likelihood of the \mjseqn{n}-th
//' data point at \mjseqn{t}-th iteration.
//'
//' @export
// [[Rcpp::export]]
std::vector<Eigen::MatrixXd> ComputePosteriorLPDFs(const std::vector<Rcpp::RawVector> & serialized_states,
												   const std::vector<std::vector<double>> & data, bool display_progress = true) {

	// Specify sizes
	int numIterations = serialized_states.size();
	int numGroups = data.size();

	// Current state buffer
	spmix::UnivariateState curr_state;

	// Prepare output buffer
	std::vector<Eigen::MatrixXd> all_post_lpdfs(numGroups);
	for (int g = 0; g < numGroups; ++g) {
		all_post_lpdfs[g] = Eigen::MatrixXd(numIterations, data[g].size());
	}

	// Populate output buffer
	if(display_progress) { REprintf("Computing log Posterior Densities:\n"); }
	Progress p_compute(numIterations, display_progress);
	for (int i = 0; i < numIterations; i++) {
		curr_state.ParseFromString(utils::raw2str(serialized_states[i]));
		auto post_lpdfs = utils::post_lpdf_from_state(curr_state, data);
        for (int g = 0; g < numGroups; ++g) {
            all_post_lpdfs[g].row(i) = post_lpdfs[g].transpose();
        }
		p_compute.increment();
	}
	return all_post_lpdfs;
}

//' Compute the log predictive densities for each group across a grid of values.
//'
//' Given the serialized MCMC chain from \code{Sampler.DensityEstimation} or \code{Sampler.BoundaryDetection}, this function
//' computes the predictive log density in each area over the points specified by \code{grid}.
//'
//' @param serialized_states A list of `raw` vectors, where each element is a serialized `UnivariateState` protobuf object.
//' Each serialized state represents a sample from the MCMC chain.
//' @param grid An numeric vector representing the grid of values at which to compute the predictive log pdf.
//' @param display_progress (Optional) a bool, if `TRUE`, it display a progress bar during the computation.
//'
//' @return A list of matrices. Each element of the list is associated to a group and contains a matrix of size
//' `number_of_iterations` times `len(grid)`, each element representing the log predictive density for that group
//' at the corresponding iteration and grid point.
//'
//' @export
// [[Rcpp::export]]
std::vector<Eigen::MatrixXd> ComputePredictiveLPDFs(const std::vector<Rcpp::RawVector> & serialized_states,
													const Eigen::VectorXd & grid, bool display_progress = true) {
	
	// Specify sizes
	int numIterations = serialized_states.size();
	int gridSize = grid.size();
	
	// Current state buffer
	spmix::UnivariateState curr_state;

	// Get number of groups
	curr_state.ParseFromString(utils::raw2str(serialized_states[0]));
	int numGroups = curr_state.groupparams_size();

	// Prepare output buffer
	std::vector<Eigen::MatrixXd> all_pred_lpdfs(numGroups, Eigen::MatrixXd(numIterations, gridSize));

	// Populate output buffer
	if(display_progress) { REprintf("Computing log Predictive Densities:\n"); }
	Progress p_compute(numIterations, display_progress);
	for (int i = 0; i < numIterations; i++) {
		curr_state.ParseFromString(utils::raw2str(serialized_states[i]));
		auto pred_lpdfs = utils::pred_lpdf_from_state(curr_state, grid);
        for (int g = 0; g < numGroups; ++g) {
            all_pred_lpdfs[g].row(i) = pred_lpdfs[g].transpose();
        }
		p_compute.increment();
	}
	return all_pred_lpdfs;
}

#endif // SPMIX_EXPORTS
