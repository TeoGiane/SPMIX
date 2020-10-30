#include "sampler_rjmcmc.h"

SpatialMixtureRJSampler::SpatialMixtureRJSampler(const SamplerParams &_params,
												 const std::vector<std::vector<double>> &_data,
												 const Eigen::MatrixXd &W,
												 const OptimOptions &_options):
SpatialMixtureSamplerBase(_params, _data, W), options(_options) {
	if (!_params.sigma_params().has_inv_gamma()) {
		throw std::runtime_error("Cannot build object of class 'SpatialMixtureRJSampler': expected parameters for an Inverse Gamma distribution.");
	}
}

SpatialMixtureRJSampler::SpatialMixtureRJSampler(const SamplerParams &_params,
												 const std::vector<std::vector<double>> &_data,
												 const Eigen::MatrixXd &W,
												 const OptimOptions &_options,
												 const std::vector<Eigen::MatrixXd> &X):
SpatialMixtureSamplerBase(_params, _data, W, X), options(_options) {
	if (!_params.sigma_params().has_inv_gamma()) {
		throw std::runtime_error("Cannot build object of class 'SpatialMixtureRJSampler': expected parameters for an Inverse Gamma distribution.");
	}
}

void SpatialMixtureRJSampler::init() {
	SpatialMixtureSamplerBase::init();
	alpha_Sigma = params.sigma_params().inv_gamma().alpha();
	beta_Sigma = params.sigma_params().inv_gamma().beta();
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

		function::spmixLogLikelihood
		loglik_reduced(data, W_init, params, numGroups, numComponents-1, rho,
					   Eigen::Map<Eigen::VectorXd>(means.data(), means.size()).head(numComponents-1),
					   Eigen::Map<Eigen::VectorXd>(sqrt_stddevs.data(), sqrt_stddevs.size()).head(numComponents-1),
					   transformed_weights_vect_reduced, Sigma.block(0,0, numComponents-2, numComponents-2));

		// Compute Acceptance rate
		double alpha = std::exp(loglik_reduced()+stan::math::poisson_lpmf(numComponents-1, 1)
					   -stan::math::poisson_lpmf(numComponents, 1));
		bool accept = stan::math::bernoulli_rng(std::min(1., alpha), rng);

		// Update state to reduce dimension
		if (accept) {
			--numComponents;
			means.resize(numComponents);
			stddevs.resize(numComponents);
			transformed_weights.conservativeResize(numGroups, numComponents);
			transformed_weights.col(numComponents-1) = Eigen::VectorXd::Zero(numGroups);
			weights.resize(numGroups, numComponents);
			for (int i = 0; i < numGroups; ++i)
				weights.row(i) = utils::InvAlr(static_cast<Eigen::VectorXd>(transformed_weights.row(i)), true);
			mtildes.conservativeResize(num_connected_comps, numComponents);
			Sigma.conservativeResize(numComponents-1, numComponents-1);
			pippo.resize(numComponents - 1);
			sigma_star_h.resize(numGroups, numComponents-1);
			_computeInvSigmaH();
		}
	}
	else { // Increase Case

		// Eliciting the approximated optimal proposal parameters
		function::spmixLogLikelihood loglik_extended(data, W_init, params, getStateAsProto());
		optimization::NewtonMethod<decltype(loglik_extended)> solver(loglik_extended, options);
		Eigen::VectorXd x0 = loglik_extended.init();
		solver.solve(x0);
		Eigen::VectorXd optMean = solver.get_state().current_minimizer;
		Eigen::MatrixXd optCov = solver.get_state().current_hessian.inverse();

		// Simulating from the approximated optimal posterior
		Eigen::VectorXd weightsMean = optMean.head(numGroups);
		Eigen::MatrixXd weightsCov = optCov.topLeftCorner(numGroups, numGroups);
		Eigen::VectorXd weights_new = stan::math::multi_normal_rng(weightsMean, weightsCov, rng);
		double means_new = stan::math::normal_rng(optMean(numGroups), optCov(numGroups, numGroups), rng);
		double var = optMean(numGroups+1)*optMean(numGroups+1)*optMean(numGroups+1)*optMean(numGroups+1);
		double varvar = optCov(numGroups+1,numGroups+1)*optCov(numGroups+1,numGroups+1)*
						optCov(numGroups+1,numGroups+1)*optCov(numGroups+1,numGroups+1);
		double stddevs_new = std::sqrt(1./stan::math::gamma_rng((var*var)/(varvar), var/varvar, rng));

		//Computing Acceptance Rate
		Eigen::VectorXd x(numGroups+2); x << weights_new, means_new, stddevs_new;
		alpha = std::exp( loglik_extended(x)+stan::math::poisson_lpmf(numComponents+1, 1)
				-loglik_extended()-stan::math::poisson_lpmf(numComponents, 1)
				-stan::math::multi_normal_lpdf(weights_new, weightsMean, weightsCov)
				-stan::math::normal_lpdf(means_new,optMean(numGroups), optCov(numGroups, numGroups))
				-stan::math::gamma_lpdf(1./(stddevs_new*stddevs_new), (var*var)/(varvar), var/varvar) );
		bool accept = stan::math::bernoulli_rng(std::min(1., alpha), rng);

		/*Rcpp::Rcout << "numComponents: " << numComponents << "\n"
		<< "means: " << Eigen::Map<Eigen::VectorXd>(means.data(), means.size()).transpose() << "\n"
		<< "stddevs: " << Eigen::Map<Eigen::VectorXd>(stddevs.data(), stddevs.size()).transpose() << "\n"
		<< "transformed_weights:\n" << transformed_weights << "\n"
		<< "weights:\n" << weights << "\n"
		<< "mtildes:\n" << mtildes << "\n"
		<< "Sigma:\n" << Sigma << "\n"
		<< "pippo:\n";
		for (auto elem : pippo) {
			Rcpp::Rcout << "    " << elem.transpose() << std::endl;
		}
		Rcpp::Rcout << "sigma_star_h:\n" << sigma_star_h << std::endl;*/

		// Update state to augment dimension
		if (accept) {
			++numComponents;
			means.emplace_back(means_new);
			stddevs.emplace_back(stddevs_new);
			transformed_weights.conservativeResize(numGroups, numComponents);
			transformed_weights.col(numComponents-2) = weights_new;
			transformed_weights.col(numComponents-1) = Eigen::VectorXd::Zero(numGroups);
			weights.conservativeResize(numGroups, numComponents);
			for (int i = 0; i < numGroups; ++i)
				weights.row(i) = utils::InvAlr(static_cast<Eigen::VectorXd>(transformed_weights.row(i)), true);
			mtildes.conservativeResize(num_connected_comps, numComponents);
			mtildes.col(numComponents-1) = Eigen::VectorXd::Zero(num_connected_comps);
			Sigma.conservativeResize(numComponents-1, numComponents-1);
			Sigma.col(numComponents-2).head(numComponents-2) =
				Sigma.row(numComponents-2).head(numComponents-2) = Eigen::VectorXd::Zero(numComponents-1);
			Sigma(numComponents-2,numComponents-2) = Sigma(0,0);
			pippo.resize(numComponents - 1);
			sigma_star_h.resize(numGroups, numComponents-1);
			_computeInvSigmaH();
		}
	}

	return;
}
