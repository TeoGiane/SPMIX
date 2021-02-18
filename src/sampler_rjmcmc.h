#ifndef RJMCMC_SAMPLER_HH
#define RJMCMC_SAMPLER_HH

#include "sampler_base.h"
#include "functors.h"
#include "gradient_ascent.h"
#include "optimization_options.pb.h"

class SpatialMixtureRJSampler: public SpatialMixtureSamplerBase {
  protected:
	// prior for Sigma --> here is an InvGamma
	double alpha_Sigma;
	double beta_Sigma;

	// prior for W --> Beta-Bernoulli prior
	//std::vector<std::vector<int>> neighbors;
	//std::vector<std::vector<double>> p;
	//double alpha_p; double beta_p;
	//double p = 0.5;

	// data range --> used in gradient ascent
	double lowerBound, upperBound;

	// Iteration counter for sample method
	int iter{0};
	//int cutoff{10};
	//int acceptedMoves{0};

	// Options for Optimization Algorithm
	OptimOptions options;

  public:
	SpatialMixtureRJSampler() = default;

	SpatialMixtureRJSampler(const SamplerParams &_params,
							const std::vector<std::vector<double>> &_data,
							const Eigen::MatrixXd &_W,
							const OptimOptions &_options,
							bool _boundary_detection);

	SpatialMixtureRJSampler(const SamplerParams &_params,
							const std::vector<std::vector<double>> &_data,
							const Eigen::MatrixXd &_W,
							const OptimOptions &_options,
							const std::vector<Eigen::MatrixXd> &X,
							bool _boundary_detection);

	void init();

	void sample() override;

	void sampleSigma() override;

	//void sampleW();

	//void sampleP();

	void labelSwitch();

	void betweenModelMove();

	void increaseMove();

	void reduceMove();

	//int get_acceptedMoves() {return acceptedMoves;};
};

#endif // RJMCMC_SAMPLER_HH
