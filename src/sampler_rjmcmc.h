#ifndef RJMCMC_SAMPLER_HH
#define RJMCMC_SAMPLER_HH

#include "cpp_proto/optimization_options.pb.h"
// #include "functors/areal_conditional_posterior_neglpdf.h"
#include "functors/conditional_posterior_neglpdf.h"
#include "optimization/LBFGS.h"
#include "sampler_base.h"

// #include "spmix_neglpdf.h"

// #include "functors.h"
// #include "gradient_ascent.h"

class SpatialMixtureRJSampler: public SpatialMixtureSamplerBase {
  protected:
	// prior for Sigma --> here is an InvGamma
	// double alpha_Sigma;
	// double beta_Sigma;

	// prior for W --> Beta-Bernoulli prior
	//std::vector<std::vector<int>> neighbors;
	//std::vector<std::vector<double>> p;
	//double alpha_p; double beta_p;
	//double p = 0.5;

	// data range --> used in gradient ascent
	double lowerBound, upperBound;

	// Iteration counter for sample method
	int numAccepted = 0;
	int itercounter = 1;
	//int cutoff{10};
	//int acceptedMoves{0};

	// Selected area
	// int selected_area;
	// int subset_size;
	// std::vector<std::vector<double>> subset_data;

	// Options for Optimization Algorithm
	LBFGSpp::LBFGSParam<double> options;
	int jump_every;
	// OptimOptions options;

  public:
	SpatialMixtureRJSampler() = default;

	SpatialMixtureRJSampler(const spmix::SamplerParams &_params,
							const std::vector<std::vector<double>> &_data,
							const Eigen::MatrixXd &_W,
							const spmix::OptimOptions &_options,
							bool _boundary_detection);

	SpatialMixtureRJSampler(const spmix::SamplerParams &_params,
							const std::vector<std::vector<double>> &_data,
							const Eigen::MatrixXd &_W,
							const spmix::OptimOptions &_options,
							const std::vector<Eigen::MatrixXd> &X,
							bool _boundary_detection);

	void init();

	void sample() override;

	// void sampleSigma() override;

	//void sampleW();

	//void sampleP();

	void labelSwitch();

	void betweenModelMove();

	void increaseMove();

	void reduceMove();

	// double computeAcceptanceRate() const { return static_cast<double>(numAccepted) / itercounter; }

	//int get_acceptedMoves() {return acceptedMoves;};
};

#endif // RJMCMC_SAMPLER_HH
