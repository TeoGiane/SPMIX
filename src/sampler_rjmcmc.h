#ifndef RJMCMC_SAMPLER_HH
#define RJMCMC_SAMPLER_HH

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

#include <algorithm>
#include <Eigen/Dense>

#include "sampler_params.pb.h"
#include "optimization_options.pb.h"
#include "utils.h"
#include "functors.h"
#include "newton_method.h"
#include "sampler_base.h"

class SpatialMixtureRJSampler: public SpatialMixtureSamplerBase {
protected:
	// prior for Sigma --> here is an InvGamma
	double alpha_Sigma;
	double beta_Sigma;

	// Options for Newton Method for Optimization
	OptimOptions options;
public:
	SpatialMixtureRJSampler() = default;
	
	SpatialMixtureRJSampler(const SamplerParams &_params,
							const std::vector<std::vector<double>> &_data,
							const Eigen::MatrixXd &W,
							const OptimOptions &_options);

	SpatialMixtureRJSampler(const SamplerParams &_params,
							const std::vector<std::vector<double>> &_data,
							const Eigen::MatrixXd &W,
							const OptimOptions &_options,
							const std::vector<Eigen::MatrixXd> &X);

	~SpatialMixtureRJSampler() = default;

	void init();

	void sample() override;

	void sampleSigma() override;

	void betweenModelMove();
};

#endif // RJMCMC_SAMPLER_HH
