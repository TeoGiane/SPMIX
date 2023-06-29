#include "sampler.h"

using namespace stan::math;

SpatialMixtureSampler::SpatialMixtureSampler(const SamplerParams &_params,
                                             const std::vector<std::vector<double>> &_data,
                                             const Eigen::MatrixXd &_W,
																						 bool _boundary_detection):
SpatialMixtureSamplerBase(_params, _data, _W) {
	// Setting boundary detection flag
	boundary_detection = _boundary_detection;
	// Control the prior for Sigma
	// if (!_params.sigma().has_inv_wishart_prior()) {
	// 	throw std::runtime_error("Cannot build object of class 'SpatialMixtureSampler': expected parameters for an Inverse Wishart distribution.");
	// }
};

SpatialMixtureSampler::SpatialMixtureSampler(const SamplerParams &_params,
                                             const std::vector<std::vector<double>> &_data,
                                             const Eigen::MatrixXd &_W,
                                             const std::vector<Eigen::MatrixXd> &X,
																						 bool _boundary_detection):
SpatialMixtureSamplerBase(_params, _data, _W, X) {
	// Setting boundary detection flag
	boundary_detection = _boundary_detection;
	// Control the prior for Sigma
	// if (!_params.sigma().has_inv_wishart_prior()) {
	// 	throw std::runtime_error("Cannot build object of class 'SpatialMixtureSampler': expected parameters for an Inverse Wishart distribution.");
	// }
};

void SpatialMixtureSampler::init() {

  // Setting variables for W sampling
	std::cout << "boundary detection: " << std::boolalpha << boundary_detection << std::endl;
  // boundary_detection = false;

	// Base class init
	SpatialMixtureSamplerBase::init();

	// Setting InvWishart Params
	// nu = params.sigma().inv_wishart_prior().nu();
	// if (params.sigma().inv_wishart_prior().identity()){
	//     V0 = Eigen::MatrixXd::Identity(numComponents - 1, numComponents - 1);
	// }
	// else {
	//     V0 = Eigen::MatrixXd::Identity(numComponents - 1, numComponents - 1);
	//     Rcpp::Rcout << "Case not implemented yet, settig V0 to identity" << std::endl;
	// }

	// Confirm
	Rcpp::Rcout << "Init done." << std::endl << std::endl;
}

void SpatialMixtureSampler::sample() {

	if (regression){
		// std::cout << "regression, ";
		regress();
		computeRegressionResiduals();
	}

	if (boundary_detection) {
		// std::cout << "boundary, ";
		sampleP();
		sampleW();
	} else {
		// std::cout << "rho, ";
		sampleRho();
	}

	// std::cout << "atoms, ";
	sampleAtoms();

	for (int i = 0; i < 2; ++i) {
		// std::cout << "allocations, ";
		sampleAllocations();
		// std::cout << "weights, ";
		sampleWeights();
	}

	// std::cout << "sigma";
	sampleSigma();
	// std::cout << "m_tilde";
	// sample_mtilde();
	// std::cout << std::endl;
}

// void SpatialMixtureSampler::sampleSigma() {
// 	Eigen::MatrixXd Vn = V0;
// 	double nu_n = nu + numGroups;
// 	Eigen::MatrixXd F_m_rhoG = F - W * rho;

//   	for (int i = 0; i < numGroups; i++) {
//     	Eigen::VectorXd wtilde_i = transformed_weights.row(i).head(numComponents - 1);
//     	Eigen::VectorXd mtilde_i = mtildes.row(node2comp[i]).head(numComponents - 1);
// 	    for (int j = 0; j < numGroups; j++) {
// 			Eigen::VectorXd wtilde_j = transformed_weights.row(j).head(numComponents - 1);
// 			Eigen::VectorXd mtilde_j = mtildes.row(node2comp[j]).head(numComponents - 1);
// 			Vn += ((wtilde_i - mtilde_i) * (wtilde_j - mtilde_j).transpose()) * F_m_rhoG(i, j);
// 	    }
// 	}
// 	Sigma = inv_wishart_rng(nu_n, Vn, rng);
// 	_computeInvSigmaH();
// }

void SpatialMixtureSampler::sample_mtilde() {

	int H = numComponents;
	Eigen::MatrixXd prec_prior = Eigen::MatrixXd::Identity(numComponents-1, numComponents-1).array()*(1.0 / mtilde_sigmasq);
	Eigen::MatrixXd F_min_rhoG = F - rho * W;

	for (int k = 0; k < num_connected_comps; k++) {

		Eigen::VectorXd currweights(comp2node[k].size() * (numComponents - 1));
		for (int i = 0; i < comp2node[k].size(); i++) {
			currweights.segment(i*(H-1), (H-1)) = transformed_weights.row(comp2node[k][i]).head(H - 1).transpose();
		}

		Eigen::MatrixXd curr_f_min_rhoG = F_by_comp[k] - rho * G_by_comp[k];
		Eigen::MatrixXd curr_prec = Eigen::kroneckerProduct(curr_f_min_rhoG, SigmaInv);
		Eigen::MatrixXd I_star = Eigen::kroneckerProduct(Eigen::VectorXd::Ones(comp2node[k].size()),
														 Eigen::MatrixXd::Identity((H-1), (H-1)));
		Eigen::MatrixXd prec_post = I_star.transpose() * curr_prec * I_star + prec_prior;
		Eigen::VectorXd m_post = prec_post.ldlt().solve(I_star.transpose() * curr_prec * currweights);
		Eigen::VectorXd sampled = stan::math::multi_normal_prec_rng(m_post, prec_post, rng);
		mtildes.row(k).head(H - 1) = sampled;
	}
}
