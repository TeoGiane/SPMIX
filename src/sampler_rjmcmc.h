#ifndef RJMCMC_SAMPLER_HH
#define RJMCMC_SAMPLER_HH

#include "utils.h"
#include "mcmc_utils.h"
#include "newton_opt.h"
#include "sampler_base.h"

class SpatialMixtureRJSampler: public SpatialMixtureSamplerBase {
protected:
	// prior for Sigma --> here is an InvGamma
	double alpha_Sigma;
	double beta_Sigma;
public:
	SpatialMixtureRJSampler() = default;
	
	SpatialMixtureRJSampler(const SamplerParams &_params,
							const std::vector<std::vector<double>> &_data,
							const Eigen::MatrixXd &W);

	SpatialMixtureRJSampler(const SamplerParams &_params,
							const std::vector<std::vector<double>> &_data,
							const Eigen::MatrixXd &W,
							const std::vector<Eigen::MatrixXd> &X);

	~SpatialMixtureRJSampler() = default;

	void init();

	void sample() override;

	void sampleSigma() override;

	void betweenModelMove();
};

#endif // RJMCMC_SAMPLER_HH
