#include "sampler_rjmcmc.h"

SpatialMixtureRJSampler::SpatialMixtureRJSampler(const SamplerParams &_params,
												 const std::vector<std::vector<double>> &_data,
												 const Eigen::MatrixXd &_W,
												 const OptimOptions &_options,
												 bool _boundary_detection):
SpatialMixtureSamplerBase(_params, _data, _W), options(_options) {

	// Setting boundary detection flag
	boundary_detection = _boundary_detection;

	// Control the prior for Sigma
	if (!_params.sigma_params().has_inv_gamma()) {
		std::string message = "Cannot build object of class 'SpatialMixtureRJSampler': "
							  "expected parameters for an Inverse Gamma distribution.";
		throw std::runtime_error(message);
	}
}

SpatialMixtureRJSampler::SpatialMixtureRJSampler(const SamplerParams &_params,
												 const std::vector<std::vector<double>> &_data,
												 const Eigen::MatrixXd &_W,
												 const OptimOptions &_options,
												 const std::vector<Eigen::MatrixXd> &X,
												 bool _boundary_detection):
SpatialMixtureSamplerBase(_params, _data, _W, X), options(_options) {

	// Setting boundary detection flag
	boundary_detection = _boundary_detection;

	// Control the prior for Sigma
	if (!_params.sigma_params().has_inv_gamma()) {
		std::string message = "Cannot build object of class 'SpatialMixtureRJSampler': "
							  "expected parameters for an Inverse Gamma distribution.";
		throw std::runtime_error(message);
	}
}

void SpatialMixtureRJSampler::init() {

	// Switch on boundary detection
	//boundary_detection = true;
	//Rcpp::Rcout << "boundary_detection? " << std::boolalpha << boundary_detection << std::endl;

	// Base class init
	SpatialMixtureSamplerBase::init();

	// Setting InvGamma Params
	alpha_Sigma = params.sigma_params().inv_gamma().alpha();
	beta_Sigma = params.sigma_params().inv_gamma().beta();

	// Setting data range
	std::tie(lowerBound, upperBound) = utils::range(data);

	// Setting sampler for Boundary detection
	//boundary_detection = true;
	//W = W_init;//Eigen::MatrixXd::Zero(numGroups,numGroups);//W_init;
	//Rcpp::Rcout << "W:\n" << W << std::endl;
	/*for (int i = 0; i < numGroups; ++i) {
		//Rcpp::Rcout << "i: " << i << std::endl;
		std::vector<int> tmp;
		std::vector<double> tmp_p;
		for (int j = i+1; j < numGroups; ++j) {
			if (W_init(i,j)){
				tmp.emplace_back(j);
				if (params.graph_params().has_beta())
					tmp_p.emplace_back(stan::math::beta_rng(params.graph_params().beta().a(),
															params.graph_params().beta().b(), rng));
				else
					tmp_p.emplace_back(params.graph_params().fixed());
			}
		}
		neighbors.emplace_back(tmp);
		p.emplace_back(tmp_p);
		Rcpp::Rcout << "neighbors: "
			<< Eigen::Map<Eigen::Matrix<int,Eigen::Dynamic,1>>(neighbors[i].data(),neighbors[i].size()).transpose() << std::endl
			<< "initial probs: " << Eigen::Map<Eigen::VectorXd>(p[i].data(), p[i].size()).transpose() << std::endl;
	}*/

	// Confirm
	Rcpp::Rcout << "Init done." << std::endl;
}

void SpatialMixtureRJSampler::sample() {

	if (regression) {
		// Rcpp::Rcout << "regression, ";
    	regress();
    	computeRegressionResiduals();
	}
	if (boundary_detection) {
		// Rcpp::Rcout << "boundary, ";
		for (int i = 0; i < 2; ++i) {
			sampleP();
			sampleW();
		}
	}
	// Rcpp::Rcout << "atoms, ";
	sampleAtoms();

	if (not boundary_detection){
		// Rcpp::Rcout << "rho, ";
		sampleRho();
	}
	// Rcpp::Rcout << "label, ";
	//labelSwitch();
	if (itercounter % 5 == 0){
		// Rcpp::Rcout << "jump, ";
		// labelSwitch();
		betweenModelMove();
	}
	for (int i = 0; i < 2; ++i) {
		// Rcpp::Rcout << "allocs, ";
		sampleAllocations();
		// Rcpp::Rcout << "weights, ";
		sampleWeights();
	}
	// Rcpp::Rcout << "sigma, ";
	sampleSigma();
	// Rcpp::Rcout << std::endl;
	//Rcpp::Rcout << "W:\n" << W << std::endl;
	++itercounter;
}

void SpatialMixtureRJSampler::sampleSigma() {

	//Rcpp::Rcout << std::endl;

	double alpha_n = alpha_Sigma + numGroups*(numComponents-1);
	double beta_n = beta_Sigma;
	Eigen::MatrixXd F_m_rhoG = F - W * rho;

	// Rcpp::Rcout << "F_m_rhoG:\n" << F_m_rhoG << std::endl;
	// Rcpp::Rcout << "mtildes:\n" << mtildes << std::endl;

	for (int i = 0; i < numGroups; i++) {
		Eigen::VectorXd wtilde_i = transformed_weights.row(i).head(numComponents - 1);
		Eigen::VectorXd mtilde_i = mtildes.row(node2comp[i]).head(numComponents - 1);
		for (int j = 0; j < numGroups; j++) {
			Eigen::VectorXd wtilde_j = transformed_weights.row(j).head(numComponents - 1);
			Eigen::VectorXd mtilde_j = mtildes.row(node2comp[j]).head(numComponents - 1);
			beta_n += ((wtilde_i - mtilde_i).dot(wtilde_j - mtilde_j)) * F_m_rhoG(i, j);

			/*Rcpp::Rcout << "(" << i << "," << j << ")\n"
				<< "wtilde_" << i << ": " << wtilde_i.transpose() << "\n"
				<< "mtilde_" << i << ": " << mtilde_i.transpose() << "\n"
				<< "wtilde_" << j << ": " << wtilde_j.transpose() << "\n"
				<< "mtilde_" << j << ": " << mtilde_j.transpose() << std::endl;*/
		}
	}

	// Rcpp::Rcout << "beta_n: " << beta_n << std::endl;

	double sigma_new = stan::math::inv_gamma_rng(alpha_n/2, beta_n/2, rng);
	Sigma = sigma_new * Eigen::MatrixXd::Identity(numComponents-1, numComponents-1);
	_computeInvSigmaH();
	return;
}

/*void SpatialMixtureRJSampler::sampleW() {

	Rcpp::Rcout << "START:\n"
	<< "W:\n" << W << std::endl << std::endl;

	//Rcpp::Rcout << "transformed_weights:\n" << transformed_weights << std::endl;

	// Initial quantities
	Eigen::MatrixXd W_uppertri = Eigen::MatrixXd::Zero(numGroups,numGroups);
	Eigen::VectorXd logProbas(2);

	for (int i = 0; i < neighbors.size(); ++i) {
		if (neighbors[i].size() > 0) {
			Rcpp::Rcout << "i: " << i;// << std::endl;
			Eigen::VectorXd wtilde_i = transformed_weights.row(i).head(numComponents-1);
			//Rcpp::Rcout << "wtilde_i: " << wtilde_i.transpose() << std::endl;
			Eigen::VectorXd mtilde_i = mtildes.row(node2comp[i]).head(numComponents-1);
			for (int j = 0; j < neighbors[i].size(); ++j) {
				Rcpp::Rcout << " j: " << neighbors[i][j];// << std::endl;
				Eigen::VectorXd wtilde_j = transformed_weights.row(neighbors[i][j]).head(numComponents - 1);
				//Rcpp::Rcout << "wtilde_j: " << wtilde_j.transpose() << std::endl;
				Eigen::VectorXd mtilde_j = mtildes.row(node2comp[neighbors[i][j]]).head(numComponents - 1);

				// Computing probabilities
				double addendum_ij = rho/(2*Sigma(0,0)) * ((wtilde_i - mtilde_i).dot(wtilde_j - mtilde_j));
				logProbas(0) = std::log(1-p[i][j]); logProbas(1) = std::log(p[i][j]) + addendum_ij;
				Eigen::VectorXd probas = logProbas.array().exp(); probas /= probas.sum();
				Rcpp::Rcout << " new_probs: " << probas.transpose() << std::endl;

				// Sampling new edge
				W_uppertri(i,neighbors[i][j]) = stan::math::categorical_rng(probas, rng)-1;
				//Rcpp::Rcout << "W(" << i << "," << neighbors[i][j] << ") = " << W_uppertri(i,neighbors[i][j]) << std::endl;
			}
			//Rcpp::Rcout << std::endl;
		}
	}

	// Computing whole W
	W = W_uppertri + W_uppertri.transpose();
	Rcpp::Rcout << std::endl;
	Rcpp::Rcout << "W:\n" << W << std::endl << "END:\n" << std::endl;

	return;
}*/

/*void SpatialMixtureRJSampler::sampleP() {

	if (params.graph_params().has_beta()) {
		double alpha_p = params.graph_params().beta().a();
		double beta_p = params.graph_params().beta().b();

		for (int i = 0; i < neighbors.size(); ++i) {
			for (int j = 0; j < neighbors[i].size(); ++j) {
				p[i][j] = stan::math::beta_rng(alpha_p+W(i,neighbors[i][j]), beta_p+1-W(i,neighbors[i][j]), rng);
			}
		}
	}
	return;
}*/

void SpatialMixtureRJSampler::labelSwitch() {

	// Randomly select the component to swap with the last one
	int to_swap = stan::math::categorical_rng(Eigen::VectorXd::Constant(numComponents,1./(numComponents)), rng) - 1;

	// Swap columns of weights ad recompute transformed weights
	weights.col(to_swap).swap(weights.col(numComponents-1));
	for (int i = 0; i < numGroups; ++i)
		transformed_weights.row(i) = utils::Alr(weights.row(i), true);

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
	// Rcpp::Rcout << "Increase." << std::endl;

	// Build the extended loglikelihood
	Eigen::MatrixXd trans_weights = transformed_weights.block(0,0,numGroups,numComponents-1);
	Eigen::Map<Eigen::VectorXd> means_map(means.data(), means.size());
	Eigen::VectorXd sqrt_stddevs(numComponents);
	for (int i = 0; i < numComponents; ++i)
		sqrt_stddevs(i) = std::sqrt(stddevs[i]);

	function::spmixLogLikelihood loglik_extended(data, W, params, numGroups, numComponents, rho,
				   								 means_map, sqrt_stddevs, trans_weights, Sigma);

	// Eliciting the approximated optimal proposal parameters
	optimization::GradientAscent<decltype(loglik_extended)> solver(loglik_extended, options);
	Eigen::VectorXd x0(numGroups+2);
	x0 << Eigen::VectorXd::Zero(numGroups), stan::math::uniform_rng(lowerBound,upperBound,rng), 1.;
	solver.solve(x0);
	Eigen::VectorXd optMean = solver.get_state().current_maximizer;
	Eigen::MatrixXd optCov = -solver.get_state().current_hessian.inverse();

	double log_arate = stan::math::negative_infinity(); Eigen::VectorXd x;
	if (solver.get_state().current_iteration < options.max_iter() and
      !solver.get_state().stagnated and
      Eigen::LDLT<Eigen::MatrixXd>(optCov).isPositive()) {

		// Simulating from the approximated optimal posterior
		x = stan::math::multi_normal_rng(optMean, optCov, rng);

		// Computing Log Acceptance Rate
		log_arate = loglik_extended(x) - loglik_extended() +
					stan::math::poisson_lpmf((numComponents+1-2),1) -
					stan::math::poisson_lpmf((numComponents-2),1) -
					stan::math::multi_normal_lpdf(x,optMean,optCov);
	}

	// Accept of Reject the move
	// bool accept = stan::math::bernoulli_rng(std::min(1., alpha), rng);

	// Update state to augment dimension
	if (std::log(stan::math::uniform_rng(0, 1, rng)) < log_arate) {
		++numComponents;
		means.resize(numComponents, means[numComponents-2]); means[numComponents-2] = x(numGroups);
		stddevs.resize(numComponents, stddevs[numComponents-2]); stddevs[numComponents-2] = x(numGroups+1)*x(numGroups+1);
		transformed_weights.conservativeResize(numGroups, numComponents);
		transformed_weights.col(numComponents-2) = x.head(numGroups);
		transformed_weights.col(numComponents-1) = Eigen::VectorXd::Zero(numGroups);
		weights.resize(numGroups, numComponents);
		for (int i = 0; i < numGroups; ++i)
			weights.row(i) = utils::InvAlr(static_cast<Eigen::VectorXd>(transformed_weights.row(i)), true);
		mtildes.conservativeResize(num_connected_comps, numComponents);
		mtildes.col(numComponents-1) = Eigen::VectorXd::Zero(num_connected_comps);
		Sigma.conservativeResize(numComponents-1, numComponents-1);
		Sigma.col(numComponents-2).head(numComponents-2) = Eigen::VectorXd::Zero(numComponents-1);
		Sigma.row(numComponents-2).head(numComponents-2) = Eigen::VectorXd::Zero(numComponents-1);
		Sigma(numComponents-2,numComponents-2) = Sigma(0,0);
		pippo.resize(numComponents-1);
		sigma_star_h.resize(numGroups, numComponents-1);
		_computeInvSigmaH();
	}
	return;
}

void SpatialMixtureRJSampler::reduceMove() {

	// Reduction Case
	// Rcpp::Rcout << "Reduce." << std::endl;

	// Randomly select the component to drop
	int to_drop = stan::math::categorical_rng(Eigen::VectorXd::Constant(numComponents-1, 1./(numComponents-1)), rng)-1;

	// Build the reduced loglikelihood
	Eigen::MatrixXd trans_weights = transformed_weights.block(0,0,numGroups,numComponents-1);
	Eigen::Map<Eigen::VectorXd> means_map(means.data(), means.size());
	Eigen::VectorXd sqrt_stddevs(numComponents);
	for (int i = 0; i < numComponents; ++i)
		sqrt_stddevs(i) = std::sqrt(stddevs[i]);

	function::spmixLogLikelihood loglik_reduced(data, W, params, numGroups, numComponents-1, rho,
				   								utils::removeElem(means_map,to_drop), utils::removeElem(sqrt_stddevs,to_drop),
				   								utils::removeColumn(trans_weights, to_drop),
				   								utils::removeRowColumn(Sigma,to_drop),to_drop);

	// Eliciting the approximated optimal proposal parameters
	Eigen::VectorXd x0(numGroups+2); x0 << trans_weights.col(to_drop), means_map(to_drop), sqrt_stddevs(to_drop);
	double fx; Eigen::VectorXd grad_fx; Eigen::MatrixXd hess_fx;
	stan::math::hessian(loglik_reduced, x0, fx, grad_fx, hess_fx);
	Eigen::MatrixXd optCov = -hess_fx.inverse();

	double log_arate = 0;
	if (Eigen::LDLT<Eigen::MatrixXd>(optCov).isPositive()) {

		// Compute Acceptance rate
		log_arate = loglik_reduced() - loglik_reduced(x0) +
					stan::math::poisson_lpmf((numComponents-1-2),1) -
					stan::math::poisson_lpmf((numComponents-2),1) +
					stan::math::multi_normal_lpdf(x0,x0,optCov);
	}

	// Accept of Reject the move
	// bool accept = stan::math::bernoulli_rng(std::min(1., alpha), rng);

	// Update state to reduce dimension
	if (std::log(stan::math::uniform_rng(0, 1, rng)) < log_arate) {
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
	}
	return;
}
