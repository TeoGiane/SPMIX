#ifndef RJMCMC_SAMPLER_HH
#define RJMCMC_SAMPLER_HH

#include "sampler_base.h"
#include "functors.h"
#include "newton_method.h"
#include "optimization_options.pb.h"

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

	void init();

	void sample() override;

	void sampleSigma() override;

	void betweenModelMove();
};

#endif // RJMCMC_SAMPLER_HH
