#include "sampler_rjmcmc.h"
#include <algorithm>

SpatialMixtureRJSampler::SpatialMixtureRJSampler(const SamplerParams &_params,
												 const std::vector<std::vector<double>> &_data,
												 const Eigen::MatrixXd &W):
SpatialMixtureSamplerBase(_params, _data, W) {
	if (!_params.sigma_params().is_invgamma()) {
		throw std::runtime_error("Cannot build object of class 'SpatialMixtureRJSampler': expected parameters for an Inverse Gamma distribution.");
	}
}

SpatialMixtureRJSampler::SpatialMixtureRJSampler(const SamplerParams &_params,
												 const std::vector<std::vector<double>> &_data,
												 const Eigen::MatrixXd &W,
												 const std::vector<Eigen::MatrixXd> &X):
SpatialMixtureSamplerBase(_params, _data, W, X) {
	if (!_params.sigma_params().is_invgamma()) {
		throw std::runtime_error("Cannot build object of class 'SpatialMixtureRJSampler': expected parameters for an Inverse Gamma distribution.");
	}
}

void SpatialMixtureRJSampler::init() {
	SpatialMixtureSamplerBase::init();
	alpha_Sigma = params.sigma_params().alpha();
	beta_Sigma = params.sigma_params().beta();
	Rcpp::Rcout << "Init done." << std::endl;
}

void SpatialMixtureRJSampler::sample() {
	/*if (regression) {
    	regress();
    	computeRegressionResiduals();
	}*/
  sampleAtoms();
  sampleAllocations();
  sampleWeights();
  sampleSigma();
  sampleRho();
  //betweenModelMove();
}

void SpatialMixtureRJSampler::sampleSigma() {
	double alpha_n = alpha_Sigma + numGroups*(numComponents - 1);
	double beta_n = beta_Sigma;
	Eigen::MatrixXd F_m_rhoG = F - W_init * rho;

	#pragma omp parallel for collapse(1)
	for (int i = 0; i < numGroups; i++) {
		Eigen::VectorXd wtilde_i = transformed_weights.row(i).head(numComponents - 1);
		Eigen::VectorXd mtilde_i = mtildes.row(node2comp[i]).head(numComponents - 1);
		for (int j = 0; j < numGroups; j++) {
			Eigen::VectorXd wtilde_j = transformed_weights.row(j).head(numComponents - 1);
			Eigen::VectorXd mtilde_j = mtildes.row(node2comp[j]).head(numComponents - 1);
			beta_n += ((wtilde_i - mtilde_i).dot(wtilde_j - mtilde_j)) * F_m_rhoG(i, j);
		}
	}

	double sigma_new = stan::math::inv_gamma_rng(alpha_n, beta_n, rng);
	Sigma = sigma_new * Eigen::MatrixXd::Identity(numComponents - 1, numComponents -1);
	_computeInvSigmaH();
	return;
}

void SpatialMixtureRJSampler::betweenModelMove() {

	// Guess an Increase or a Reduction in the number of components
	bool increase;
	if (numComponents == 1)
		increase = true;
	else
		increase = stan::math::bernoulli_rng(0.5, rng);

	// Reduction case
	if (not increase) {

		Eigen::VectorXd transformed_weights_vect_reduced =
		Eigen::Map<Eigen::VectorXd>(transformed_weights.block(0,0, numGroups, numComponents - 2).data(),
									transformed_weights.block(0,0, numGroups, numComponents - 2).size());
		std::vector<double> sqrt_stddevs;
		for (auto elem : stddevs)
			sqrt_stddevs.emplace_back(std::sqrt(elem));

		utils::spmixLogLikelihood
		loglik_reduced(data, W_init, params, numGroups, numComponents-1, rho,
					   Eigen::Map<Eigen::VectorXd>(means.data(), means.size()).head(numComponents-1),
					   Eigen::Map<Eigen::VectorXd>(sqrt_stddevs.data(), sqrt_stddevs.size()).head(numComponents-1),
					   transformed_weights_vect_reduced, Sigma.block(0,0, numComponents-2, numComponents-2));

		// Compute Acceptance rate
		double alpha = std::exp(loglik_reduced()+stan::math::poisson_lpmf(numComponents-1, 1)-stan::math::poisson_lpmf(numComponents, 1));
		bool accept = stan::math::bernoulli_rng(std::min(1., alpha), rng);

		// Update state to reduce dimension
		if (accept) {
			/* code */
		}
	}

	// if +1 --> newton_opt to simulate form pi_tilde and compute acceptance rates
	// in any case, update the state of the sampler (maybe is required an update() utility? non credo in realtà)
	return;
}
