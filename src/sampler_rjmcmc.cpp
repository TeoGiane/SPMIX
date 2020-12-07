#include "sampler_rjmcmc.h"

SpatialMixtureRJSampler::SpatialMixtureRJSampler(const SamplerParams &_params,
												 const std::vector<std::vector<double>> &_data,
												 const Eigen::MatrixXd &W,
												 const OptimOptions &_options):
SpatialMixtureSamplerBase(_params, _data, W), options(_options) {
	if (!_params.sigma_params().has_inv_gamma()) {
		std::string message = "Cannot build object of class 'SpatialMixtureRJSampler': "
							  "expected parameters for an Inverse Gamma distribution.";
		throw std::runtime_error(message);
	}
}

SpatialMixtureRJSampler::SpatialMixtureRJSampler(const SamplerParams &_params,
												 const std::vector<std::vector<double>> &_data,
												 const Eigen::MatrixXd &W,
												 const OptimOptions &_options,
												 const std::vector<Eigen::MatrixXd> &X):
SpatialMixtureSamplerBase(_params, _data, W, X), options(_options) {
	if (!_params.sigma_params().has_inv_gamma()) {
		std::string message = "Cannot build object of class 'SpatialMixtureRJSampler': "
							  "expected parameters for an Inverse Gamma distribution.";
		throw std::runtime_error(message);
	}
}

void SpatialMixtureRJSampler::init() {

	// Base class init
	SpatialMixtureSamplerBase::init();
	
	// DEBUG
	/*while (utils::numComponentsPrior(cutoff,0.3,1.1) > 1e-8)
		cutoff++;*/
	Rcpp::Rcout << "cutoff: " << cutoff << std::endl;
	
	// Setting InvGamma Params
	alpha_Sigma = params.sigma_params().inv_gamma().alpha();
	beta_Sigma = params.sigma_params().inv_gamma().beta();

	// Setting data range
	std::tie(lowerBound, upperBound) = utils::range(data);

	// Confirm
	Rcpp::Rcout << "Init done." << std::endl;
}

void SpatialMixtureRJSampler::sample() {
	/*if (regression) {
    	regress();
    	computeRegressionResiduals();
	}*/
	sampleAtoms();
	sampleWeights();
	sampleSigma();
	sampleRho();
	//labelSwitch();
	betweenModelMove();
	sampleAllocations();
}

void SpatialMixtureRJSampler::sampleSigma() {
	double alpha_n = alpha_Sigma + numGroups*(numComponents-1);
	double beta_n = beta_Sigma;
	Eigen::MatrixXd F_m_rhoG = F - W_init * rho;

	for (int i = 0; i < numGroups; i++) {
		Eigen::VectorXd wtilde_i = transformed_weights.row(i).head(numComponents - 1);
		Eigen::VectorXd mtilde_i = mtildes.row(node2comp[i]).head(numComponents - 1);
		for (int j = 0; j < numGroups; j++) {
			Eigen::VectorXd wtilde_j = transformed_weights.row(j).head(numComponents - 1);
			Eigen::VectorXd mtilde_j = mtildes.row(node2comp[j]).head(numComponents - 1);
			beta_n += ((wtilde_i - mtilde_i).dot(wtilde_j - mtilde_j)) * F_m_rhoG(i, j);
		}
	}

	double sigma_new = stan::math::inv_gamma_rng(alpha_n/2, (1./beta_n)/2, rng);
	Sigma = sigma_new * Eigen::MatrixXd::Identity(numComponents-1, numComponents-1);
	_computeInvSigmaH();
	return;
}

void SpatialMixtureRJSampler::labelSwitch() {

	// Randomly select the component to swap with the last one
	int to_swap = stan::math::categorical_rng(Eigen::VectorXd::Constant(numComponents,1./(numComponents)), rng)-1;

	// Swap columns of weights ad recompute transformed weights
	weights.col(to_swap).swap(weights.col(numComponents-1));
	// #pragma omp parallel for
	for (int i = 0; i < numGroups; ++i)
		transformed_weights.row(i) = utils::Alr(weights.row(i), true);

	// Swap postNormalGammaParams rows
	postNormalGammaParams.row(to_swap).swap(postNormalGammaParams.row(numComponents-1));

	// Swap Atoms
	std::iter_swap(means.begin()+to_swap, means.begin()+(numComponents-1));
	std::iter_swap(stddevs.begin()+to_swap, stddevs.begin()+(numComponents-1));

	return;
}

void SpatialMixtureRJSampler::betweenModelMove() {

	// Guess an Increase or a Reduction in the number of components
	bool increase;
	if (numComponents == 2)
		increase = true;
	else if (numComponents == cutoff)
		increase = false;
	else
		increase = stan::math::bernoulli_rng(0.5, rng);

	// Select the proper move
	if (not increase)
		reduceMove();
	else
		increaseMove();

	return;
};

void SpatialMixtureRJSampler::increaseMove() {

	// Increase Case
	//Rcpp::Rcout << "Increase." << std::endl;

	// Build the extended loglikelihood
	Eigen::MatrixXd trans_weights = transformed_weights.block(0,0,numGroups,numComponents-1);
	Eigen::Map<Eigen::VectorXd> means_map(means.data(), means.size());
	Eigen::VectorXd sqrt_stddevs(numComponents);
	for (int i = 0; i < numComponents; ++i)
		sqrt_stddevs(i) = std::sqrt(stddevs[i]);

	Eigen::MatrixXd priorMatrix(numComponents,4);
	for (int i = 0; i < numComponents; ++i){
		Eigen::VectorXd tmp(4); tmp << priorMean, priorA, priorB, priorLambda;
		priorMatrix.row(i) = tmp;
	}

	function::spmixLogLikelihood loglik_extended(data, W_init, params, numGroups, numComponents, rho,
				   								 means_map, sqrt_stddevs, priorMatrix,
				   								 trans_weights, Sigma);

	// Eliciting the approximated optimal proposal parameters
	optimization::GradientAscent<decltype(loglik_extended)> solver(loglik_extended, options);
	int randComp = stan::math::categorical_rng(Eigen::VectorXd::Constant(numComponents-1, 1./(numComponents-1)), rng)-1;
	//Eigen::VectorXd x0 = loglik_extended.init(randComp);
	Eigen::VectorXd x0(numGroups+2);
	//Rcpp::Rcout << "low_bound: " << lowerBound << ", up_bound: " << upperBound << std::endl;
	//double low_bound = *std::min_element(data[0].begin(),data[0].end());
	//double up_bound = *std::max_element(data[0].begin(),data[0].end());
	x0 << Eigen::VectorXd::Zero(numGroups), stan::math::uniform_rng(lowerBound,upperBound,rng), 1.;
	solver.solve(x0);
	//Rcpp::Rcout << "Ended after " << solver.get_state().current_iteration << " iterations" << std::endl;
	Eigen::VectorXd optMean = solver.get_state().current_minimizer;
	Eigen::MatrixXd optCov = -solver.get_state().current_hessian.inverse();

	double alpha{0.}; Eigen::VectorXd x;
	//Rcpp::Rcout << "Has stagnated? " << std::boolalpha << solver.get_state().stagnated << std::endl;
	//Rcpp::Rcout << "Is negative stdev? " << std::boolalpha << (solver.get_state().current_minimizer(numGroups+1) < 0.) << std::endl;
	if (solver.get_state().current_iteration < options.max_iter() and !solver.get_state().stagnated) {

		// Simulating from the approximated optimal posterior
		x = stan::math::multi_normal_rng(optMean, optCov, rng);

		/*Rcpp::Rcout << "trans_weights: " << transformed_weights << std::endl;
		Rcpp::Rcout << "x added: " << x.transpose() << std::endl;*/
		/*Rcpp::Rcout << "(+)loglik_extended(x):\n" << loglik_extended(x) << "\n"
		<< "(-)loglik_extended():\n" << loglik_extended() << "\n"
		<< "(-)stan::math::multi_normal_lpdf(x,optMean,optCov): " << stan::math::multi_normal_lpdf(x,optMean,optCov) << "\n";*/

		//Computing Acceptance Rate
		alpha = std::exp(loglik_extended(x)-loglik_extended()-stan::math::multi_normal_lpdf(x,optMean,optCov)) *
				//utils::numComponentsPrior(numComponents+1,0.3,1.1)/utils::numComponentsPrior(numComponents,0.3,1.1);
				std::exp(stan::math::poisson_lpmf((numComponents+1-2),1)-stan::math::poisson_lpmf((numComponents-2),1));
	}

	//Rcpp::Rcout << "mean_diff: " << (means_map-Eigen::VectorXd::Constant(means_map.size(),x(numGroups))).transpose() << std::endl;
	// Accept of Reject the move
	bool accept = stan::math::bernoulli_rng(std::min(1., alpha), rng);
	//Rcpp::Rcout << "alpha: " << alpha << ". Accept? " << std::boolalpha << accept << std::endl << std::endl;

	// Update state to augment dimension
	if (accept) {
		++acceptedMoves;
		++numComponents;

		// HERE WE INDUCE A SWITCH BETWEEN LAST AND BEFORE LAST COMPONENT!
		means.resize(numComponents, means[numComponents-2]); means[numComponents-2] = x(numGroups);
		stddevs.resize(numComponents, stddevs[numComponents-2]); stddevs[numComponents-2] = x(numGroups+1)*x(numGroups+1);
		//means.emplace_back(x(numGroups));
		//stddevs.emplace_back(x(numGroups+1)*x(numGroups+1));

		transformed_weights.conservativeResize(numGroups, numComponents);
		transformed_weights.col(numComponents-2) = x.head(numGroups);
		transformed_weights.col(numComponents-1) = Eigen::VectorXd::Zero(numGroups);
		weights.resize(numGroups, numComponents);
		for (int i = 0; i < numGroups; ++i)
			weights.row(i) = utils::InvAlr(static_cast<Eigen::VectorXd>(transformed_weights.row(i)), true);

		// THE ATOMS SWITCH IS FOLLOWED BY THE WEIGHTS ONE!
		/*weights.col(numComponents-2).swap(weights.col(numComponents-1));
		for (int i = 0; i < numGroups; ++i){
			transformed_weights.row(i) = utils::Alr(weights.row(i),true);
		}*/

		mtildes.conservativeResize(num_connected_comps, numComponents);
		mtildes.col(numComponents-1) = Eigen::VectorXd::Zero(num_connected_comps);
		Sigma.conservativeResize(numComponents-1, numComponents-1);
		Sigma.col(numComponents-2).head(numComponents-2) = Eigen::VectorXd::Zero(numComponents-1);
		Sigma.row(numComponents-2).head(numComponents-2) = Eigen::VectorXd::Zero(numComponents-1);
		Sigma(numComponents-2,numComponents-2) = Sigma(0,0);
		pippo.resize(numComponents-1);
		sigma_star_h.resize(numGroups, numComponents-1);
		_computeInvSigmaH();
		//sampleAllocations();
		//Rcpp::Rcout << "State Enlarged" << std::endl;
	}
	return;
}

void SpatialMixtureRJSampler::reduceMove() {

	// Reduction Case
	//Rcpp::Rcout << "Reduce." << std::endl;

	// Randomly select the component to drop
	int to_drop = stan::math::categorical_rng(Eigen::VectorXd::Constant(numComponents-1, 1./(numComponents-1)), rng)-1;

	// Build the reduced loglikelihood
	Eigen::MatrixXd trans_weights = transformed_weights.block(0,0,numGroups,numComponents-1);
	Eigen::Map<Eigen::VectorXd> means_map(means.data(), means.size());
	Eigen::VectorXd sqrt_stddevs(numComponents);
	for (int i = 0; i < numComponents; ++i)
		sqrt_stddevs(i) = std::sqrt(stddevs[i]);

	Eigen::MatrixXd priorMatrix(numComponents,4);
	for (int i = 0; i < numComponents; ++i){
		Eigen::VectorXd tmp(4); tmp << priorMean, priorA, priorB, priorLambda;
		priorMatrix.row(i) = tmp;
	}

	function::spmixLogLikelihood loglik_reduced(data, W_init, params, numGroups, numComponents-1, rho,
				   								utils::removeElem(means_map,to_drop), utils::removeElem(sqrt_stddevs,to_drop),
				   								priorMatrix, utils::removeColumn(trans_weights, to_drop),
				   								utils::removeRowColumn(Sigma,to_drop),to_drop);

	// Eliciting the approximated optimal proposal parameters
	//optimization::GradientAscent<decltype(loglik_reduced)> solver(loglik_reduced, options);
	Eigen::VectorXd x0(numGroups+2); x0 << trans_weights.col(to_drop), means_map(to_drop), sqrt_stddevs(to_drop);
	double fx; Eigen::VectorXd grad_fx; Eigen::MatrixXd hess_fx;
	stan::math::hessian(loglik_reduced, x0, fx, grad_fx, hess_fx);
	Eigen::MatrixXd optCov = -hess_fx.inverse();
	//solver.solve(x0);
	//Rcpp::Rcout << "Ended after " << solver.get_state().current_iteration << " iterations" << std::endl;
	//Eigen::VectorXd optMean = solver.get_state().current_minimizer;
	//Eigen::MatrixXd optCov = -solver.get_state().current_hessian.inverse();

	double alpha{0.};
	//Rcpp::Rcout << "Has stagnated? " << std::boolalpha << solver.get_state().stagnated << std::endl;
	//Rcpp::Rcout << "Is negative stdev? " << std::boolalpha << (solver.get_state().current_minimizer(numGroups+1) < 0.) << std::endl;
	if (Eigen::LDLT<Eigen::MatrixXd>(optCov).isPositive()) { //solver.get_state().current_iteration < options.max_iter() and !solver.get_state().stagnated) {

		/*Rcpp::Rcout << "trans_weights: " << transformed_weights << std::endl;
		Rcpp::Rcout << "x removed: " << x0.transpose() << std::endl;*/
		/*Rcpp::Rcout << "(+)loglik_reduced():\n" << loglik_reduced() << "\n"
		<< "(-)loglik_reduced(x0):\n" << loglik_reduced(x0) << "\n"
		<< "(+)stan::math::multi_normal_lpdf(x0,x0,optCov): " << stan::math::multi_normal_lpdf(x0,x0,optCov) << "\n";*/

		// Compute Acceptance rate
		alpha = std::exp(loglik_reduced()-loglik_reduced(x0)+stan::math::multi_normal_lpdf(x0,x0,optCov)) *
				std::exp(stan::math::poisson_lpmf((numComponents-1-2),1)-stan::math::poisson_lpmf((numComponents-2),1));
				//utils::numComponentsPrior(numComponents-1,0.3,1.1)/utils::numComponentsPrior(numComponents,0.3,1.1);
	} else {
		alpha = 1.;
	}

	// Accept of Reject the move
	bool accept = stan::math::bernoulli_rng(std::min(1., alpha), rng);
	//Rcpp::Rcout << "alpha: " << alpha << ". Accept? " << std::boolalpha << accept << std::endl << std::endl;

	// Update state to reduce dimension
	if (accept) {
		++acceptedMoves;
		--numComponents;
		means.erase(means.begin()+to_drop);
		stddevs.erase(stddevs.begin()+to_drop);
		transformed_weights = utils::removeColumn(transformed_weights, to_drop);
		weights.resize(numGroups, numComponents);
		for (int i = 0; i < numGroups; ++i)
			weights.row(i) = utils::InvAlr(static_cast<Eigen::VectorXd>(transformed_weights.row(i)), true);
		mtildes.conservativeResize(num_connected_comps, numComponents);
		Sigma = utils::removeRowColumn(Sigma, to_drop);
		pippo.resize(numComponents-1);
		sigma_star_h.resize(numGroups, numComponents-1);
		_computeInvSigmaH();
		//sampleAllocations();
		//Rcpp::Rcout << "State Reduced" << std::endl;
	}
	return;
}
