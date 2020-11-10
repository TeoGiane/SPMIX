#ifndef SPMIX_SAMPLER_BASE_HH
#define SPMIX_SAMPLER_BASE_HH

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
#include <random>
#include <vector>
#include <numeric>
#include <stdexcept>
#include <Eigen/Dense>
#include <unsupported/Eigen/KroneckerProduct>

#include "collector.h"
#include "univariate_mixture_state.pb.h"
#include "sampler_params.pb.h"
#include "PolyaGammaHybrid.h"
#include "mcmc_utils.h"
#include "utils.h"

class SpatialMixtureSamplerBase {
  protected:
	// Params
	SamplerParams params;

	// data
	int numGroups;
	std::vector<int> samplesPerGroup;
	std::vector<std::vector<double>> data;
	int numdata;
	Eigen::MatrixXd W_init;

	// mixtures
	int numComponents;
	std::vector<double> means;
	std::vector<double> stddevs;

	Eigen::MatrixXd weights; // one set of weights per location
	Eigen::MatrixXd transformed_weights;
	std::vector<std::vector<int>> cluster_allocs;
	Eigen::MatrixXd mtildes;
	int num_connected_comps;

	// MCAR
	double rho;
	Eigen::MatrixXd Sigma;
	Eigen::MatrixXd W;
	Eigen::MatrixXd SigmaInv;
	Eigen::MatrixXd F;
	std::vector<int> node2comp;
	double mtilde_sigmasq = 1.0;
	std::vector<std::vector<int>> comp2node;
	std::vector<Eigen::MatrixXd> F_by_comp;
	std::vector<Eigen::MatrixXd> G_by_comp;

	// Regression
	bool regression = false;
	int p_size;
	Eigen::VectorXd reg_coeff;
	Eigen::MatrixXd predictors;
	Eigen::VectorXd reg_coeff_mean;
	Eigen::MatrixXd reg_coeff_prec;
	Eigen::VectorXd reg_data;
	Eigen::VectorXd mu;
	Eigen::DiagonalMatrix<double, Eigen::Dynamic, Eigen::Dynamic> V;

	// prior for Sigma --> depends on the derivation

	// prior for Rho
	double alpha;
	double beta;

	// adaptive MCMC for rho
	double sigma_n_rho;
	double rho_sum;
	double rho_sum_sq;
	int iter = 0;

	std::vector<Eigen::VectorXd> pippo;
	Eigen::MatrixXd sigma_star_h;

	// HyperParams for NormalGamma
	double priorMean, priorA, priorB, priorLambda;

	unsigned long seed = 213513435;
	PolyaGammaHybridDouble* pg_rng = nullptr;
	std::mt19937_64 rng{213513435};

	// diagnostic for the MH sampler
	int numAccepted = 0;

  public:
	SpatialMixtureSamplerBase() {}

	SpatialMixtureSamplerBase(
		const SamplerParams &_params,
		const std::vector<std::vector<double>> &_data,
		const Eigen::MatrixXd &W);

	SpatialMixtureSamplerBase(
		const SamplerParams &_params,
		const std::vector<std::vector<double>> &_data,
		const Eigen::MatrixXd &W, const std::vector<Eigen::MatrixXd> &X);

	virtual ~SpatialMixtureSamplerBase() {
		//Rcpp::Rcout << "SamplerBase destructor!" << std::endl;
		delete(pg_rng);
	}

	void init();

	virtual void sample() = 0;
	
	/*
	 * We use a Normal kernel with conjugate Normal - Inverse Gamma
	 * base measure, so the update of the atom is
	 * Law(\mu, \sigma)_h \propto
	 *  P_0(\mu, \sigma) \prod_{(i, j): s_{ij}=h} N(y_{ij} | \mu_h, \sigma_h)
	 * That is a conjugate normal likelihood with normal inverse gamma prior
	 */
	void sampleAtoms();

	void sampleAllocations();

	/*
	 * We use the PolyaGamma trick to sample from the transformed weights
	 */
	void sampleWeights();

	/*
	 * This step requires a Metropolis Hastings step
	 */
	void sampleRho();

	/*
	 * We use a conjugate Inverse - Wishart prior for Sigma, so the
	 * posterior law of Sigma is still Inverse - Wishart with parameters
	 * (\nu + I, Psi + \sum_{i=1}^I (tw_i - \mu_i) (tw_i - \mu_i)^T
	 * for tw = transformed weights
	 * \mu_i = \rho N^{-1} \sum{j \n N(i)} tw_j
	 */
	virtual void sampleSigma() = 0;

	void regress();

	void computeRegressionResiduals();

	void _computeInvSigmaH();

	/*
	 * Sampler the hyperparameters in the base measure P_0
	 */
	void sampleHyperParams();

	//void sample_mtilde();

	void saveState(Collector<UnivariateState>* collector);

	UnivariateState getStateAsProto();

	void printDebugString();

	const int& getNumAccepted() {return numAccepted;}
};

#endif // SPMIX_SAMPLER_BASE_HH