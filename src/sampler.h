#ifndef SRC_SAMPLER_HPP
#define SRC_SAMPLER_HPP

#include "sampler_base.h"

class SpatialMixtureSampler: public SpatialMixtureSamplerBase {
  protected:

	// prior for Sigma --> here is an InvWishart
	double nu;
	Eigen::MatrixXd V0;

  public:
	SpatialMixtureSampler() {}

	SpatialMixtureSampler(
		const SamplerParams &_params,
		const std::vector<std::vector<double>> &_data,
		const Eigen::MatrixXd &_W);

	SpatialMixtureSampler(
		const SamplerParams &_params,
		const std::vector<std::vector<double>> &_data,
		const Eigen::MatrixXd &_W, const std::vector<Eigen::MatrixXd> &X);

    void init();

    void sample() override;

    /*
     * We use a conjugate Inverse - Wishart prior for Sigma, so the
     * posterior law of Sigma is still Inverse - Wishart with parameters
     * (\nu + I, Psi + \sum_{i=1}^I (tw_i - \mu_i) (tw_i - \mu_i)^T
     * for tw = transformed weights
     * \mu_i = \rho N^{-1} \sum{j \n N(i)} tw_j
     */
    void sampleSigma() override;

    void sample_mtilde();
};

#endif // SRC_SAMPLER_HPP
