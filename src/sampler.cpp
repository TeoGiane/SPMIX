#include "sampler.h"

using namespace stan::math;

SpatialMixtureSampler::SpatialMixtureSampler(const SamplerParams &_params,
                                             const std::vector<std::vector<double>> &_data,
                                             const Eigen::MatrixXd &W):
SpatialMixtureSamplerBase(_params, _data, W) {
    if (!_params.sigma_params().has_inv_wishart()) {
        throw std::runtime_error("Cannot build object of class 'SpatialMixtureSampler': expected parameters for an Inverse Wishart distribution.");
    }
};

SpatialMixtureSampler::SpatialMixtureSampler(const SamplerParams &_params,
                                             const std::vector<std::vector<double>> &_data,
                                             const Eigen::MatrixXd &W,
                                             const std::vector<Eigen::MatrixXd> &X):
SpatialMixtureSamplerBase(_params, _data, W, X) {
    if (!_params.sigma_params().has_inv_wishart()) {
        throw std::runtime_error("Cannot build object of class 'SpatialMixtureSampler': expected parameters for an Inverse Wishart distribution.");
    }
};

void SpatialMixtureSampler::init() {
    SpatialMixtureSamplerBase::init();
    nu = params.sigma_params().inv_wishart().nu();
    if (params.sigma_params().inv_wishart().identity()){
        V0 = Eigen::MatrixXd::Identity(numComponents - 1, numComponents - 1);
    }
    else {
        V0 = Eigen::MatrixXd::Identity(numComponents - 1, numComponents - 1);
        Rcpp::Rcout << "Case not implemented yet, settig V0 to identity" << std::endl;
    }
    Rcpp::Rcout << "Init done." << std::endl;
}

/*SpatialMixtureSampler::SpatialMixtureSampler(
    const SamplerParams &_params, const std::vector<std::vector<double>> &_data,
    const Eigen::MatrixXd &W)
    : params(_params), data(_data), W_init(W) {
  numGroups = data.size();
  samplesPerGroup.resize(numGroups);
  for (int i = 0; i < numGroups; i++) {
    samplesPerGroup[i] = data[i].size();
  }
  numdata = std::accumulate(samplesPerGroup.begin(), samplesPerGroup.end(), 0);
}

SpatialMixtureSampler::SpatialMixtureSampler(
    const SamplerParams &_params, const std::vector<std::vector<double>> &_data,
    const Eigen::MatrixXd &W, const std::vector<Eigen::MatrixXd> &X)
    : params(_params), data(_data), W_init(W) {
  numGroups = data.size();
  samplesPerGroup.resize(numGroups);
  for (int i = 0; i < numGroups; i++) {
    samplesPerGroup[i] = data[i].size();
  }
  numdata = std::accumulate(samplesPerGroup.begin(), samplesPerGroup.end(), 0);

  if (X.size() > 0) {
    regression = true;
    p_size = X[0].cols();
    reg_coeff_mean = Eigen::VectorXd::Zero(p_size);
    reg_coeff_prec = Eigen::MatrixXd::Identity(p_size, p_size);
    reg_coeff =
        stan::math::multi_normal_rng(reg_coeff_mean, 10 * reg_coeff_prec, rng);

    predictors.resize(numdata, p_size);
    reg_data.resize(numdata);
    V.resize(numdata);
    mu.resize(numdata);
    int start = 0;
    for (int i = 0; i < numGroups; i++) {
      predictors.block(start, 0, samplesPerGroup[i], p_size) = X[i];
      reg_data.segment(start, samplesPerGroup[i]) =
          Eigen::Map<Eigen::VectorXd>(data[i].data(), samplesPerGroup[i]);
      start += samplesPerGroup[i];
    }
    computeRegressionResiduals();
  }
}*/

/*void SpatialMixtureSampler::init() {
  pg_rng = new PolyaGammaHybridDouble(seed);

  numComponents = params.num_components();
  mtilde_sigmasq = params.mtilde_sigmasq();

  priorMean = params.p0_params().mu0();
  priorA = params.p0_params().a();
  priorB = params.p0_params().b();
  priorLambda = params.p0_params().lam_();

  if (params.sigma_params().has_inv_wishart()) {
    nu = params.sigma_params().inv_wishart().nu();
    if (params.sigma_params().inv_wishart().identity())
      V0 = Eigen::MatrixXd::Identity(numComponents - 1, numComponents - 1);
    else {
      V0 = Eigen::MatrixXd::Identity(numComponents - 1, numComponents - 1);
      Rcpp::Rcout << "Case not implemented yet, settig V0 to identity" << std::endl;
    }
  }
  else {
    std::runtime_error("InvGamma distribution parameters not allowed.");
  }

  alpha = params.rho_params().a();
  beta = params.rho_params().b();

  node2comp = utils::findConnectedComponents(W_init);
  auto it = std::max_element(node2comp.begin(), node2comp.end());
  num_connected_comps = *it + 1;

  comp2node.resize(num_connected_comps);
  for (int i = 0; i < numGroups; i++) {
    comp2node[node2comp[i]].push_back(i);
  }

  // Now proper initialization
  rho = 0.9;
  rho_sum = 0;
  rho_sum_sq = 0;
  Sigma = Eigen::MatrixXd::Identity(numComponents - 1, numComponents - 1);
  means.resize(numComponents);
  stddevs.resize(numComponents);
  weights = Eigen::MatrixXd::Zero(numGroups, numComponents);
  transformed_weights = Eigen::MatrixXd::Zero(numGroups, numComponents);
  cluster_allocs.resize(numGroups);
  for (int i = 0; i < numGroups; i++) {
    cluster_allocs[i].resize(samplesPerGroup[i]);
  }

  for (int h = 0; h < numComponents; h++) {
    means[h] = normal_rng(0.0, 10.0, rng);
    stddevs[h] = uniform_rng(0.5, 2.0, rng);
  }

  for (int i = 0; i < numGroups; i++) {
    weights.row(i) = dirichlet_rng(Eigen::VectorXd::Ones(numComponents), rng);
    transformed_weights.row(i) = utils::Alr(weights.row(i), true);
  }

  for (int i = 0; i < numGroups; i++) {
    for (int j = 0; j < std::min(numComponents, samplesPerGroup[i]); j++)
      cluster_allocs[i][j] = j;

    for (int j = numComponents; j < samplesPerGroup[i]; j++)
      cluster_allocs[i][j] = categorical_rng(weights.row(i), rng) - 1;
  }

  // W = W_init;
  // // normalize W
  // for (int i=0; i < W.rows(); ++i){
  //   W.row(i) *= rho/W.row(i).sum();
  // }

  F = Eigen::MatrixXd::Zero(numGroups, numGroups);
  for (int i = 0; i < numGroups; i++)
    F(i, i) = rho * W_init.row(i).sum() + (1 - rho);

  F_by_comp.resize(num_connected_comps);
  G_by_comp.resize(num_connected_comps);
  for (int k = 0; k < num_connected_comps; k++) {
    Eigen::MatrixXd curr_f =
        Eigen::MatrixXd::Zero(comp2node[k].size(), comp2node[k].size());
    Eigen::MatrixXd curr_g =
        Eigen::MatrixXd::Zero(comp2node[k].size(), comp2node[k].size());
    for (int i = 0; i < comp2node[k].size(); i++) {
      curr_f(i, i) = F(comp2node[k][i], comp2node[k][i]);
      for (int j = 0; j < comp2node[k].size(); j++) {
        curr_g(i, j) = W_init(comp2node[k][i], comp2node[k][j]);
      }
    }
    F_by_comp[k] = curr_f;
    G_by_comp[k] = curr_g;
  }

  // last component is not used!
  mtildes = Eigen::MatrixXd::Zero(num_connected_comps, numComponents);

  pippo.resize(numComponents - 1);
  sigma_star_h = Eigen::MatrixXd::Zero(numGroups, numComponents - 1);

  _computeInvSigmaH();

  Rcpp::Rcout << "Init done" << std::endl;
}*/

void SpatialMixtureSampler::sample() {
  if (regression) {
    regress();
    computeRegressionResiduals();
  }
  sampleAtoms();
  sampleAllocations();
  sampleWeights();
  sampleSigma();
  sampleRho();
  //sample_mtilde();
}

/*void SpatialMixtureSampler::sampleAtoms() {
  std::vector<std::vector<double>> datavec(numComponents);
  for (int h = 0; h < numComponents; h++) datavec[h].reserve(numdata);

  for (int i = 0; i < numGroups; i++) {
    for (int j = 0; j < samplesPerGroup[i]; j++) {
      int comp = cluster_allocs[i][j];
      datavec[comp].push_back(data[i][j]);
    }
  }

#pragma omp parallel for
  for (int h = 0; h < numComponents; h++) {
    std::vector<double> params = utils::normalGammaUpdate(
        datavec[h], priorMean, priorA, priorB, priorLambda);
    double tau = stan::math::gamma_rng(params[1], params[2], rng);
    double sigmaNorm = 1.0 / std::sqrt(tau * params[3]);
    double mu = stan::math::normal_rng(params[0], sigmaNorm, rng);
    means[h] = mu;
    stddevs[h] = 1.0 / std::sqrt(tau);
  }
}

void SpatialMixtureSampler::sampleAllocations() {
	for (int i = 0; i < numGroups; i++) {

		#pragma omp parallel for
		for (int j = 0; j < samplesPerGroup[i]; j++) {
	  		double datum = data[i][j];
	  		Eigen::VectorXd logProbas(numComponents);
	  		for (int h = 0; h < numComponents; h++) {
	    		logProbas(h) = std::log(weights(i, h) + 1e-6) + normal_lpdf(datum, means[h], stddevs[h]);
	  		}

	  		Eigen::VectorXd probas = logProbas.array().exp();
	  		probas /= probas.sum();
	  		cluster_allocs[i][j] = categorical_rng(probas, rng) - 1;
		}
	}
}*/

/*void SpatialMixtureSampler::sampleWeights() {
  for (int i = 0; i < numGroups; i++) {
    std::vector<int> cluster_sizes(numComponents, 0);

    #pragma omp parallel for
    for (int j = 0; j < samplesPerGroup[i]; j++)
      cluster_sizes[cluster_allocs[i][j]] += 1;

    for (int h = 0; h < numComponents - 1; h++) {
  		
   		// we draw omega from a Polya-Gamma distribution
		Eigen::VectorXd weightsForCih =
		  utils::removeElem(transformed_weights.row(i), h);
		double C_ih = stan::math::log_sum_exp(weightsForCih);

		double omega_ih =
		  pg_rng->draw(samplesPerGroup[i], transformed_weights(i, h) - C_ih);

		Eigen::VectorXd mu_i =
		  (W_init.row(i) * transformed_weights).array() * rho +
		  mtildes.row(node2comp[i]).array() * (1 - rho);
		mu_i = mu_i.array() / (W_init.row(i).sum() * rho + 1 - rho);
		mu_i = mu_i.head(numComponents - 1);
		Eigen::VectorXd wtilde =
		  transformed_weights.row(i).head(numComponents - 1);

		double mu_star_ih = mu_i[h] + pippo[h].dot(utils::removeElem(wtilde, h) -
		                                         utils::removeElem(mu_i, h));

		double sigma_hat_ih = 1.0 / (1.0 / sigma_star_h(i, h) + omega_ih);
		int N_ih = cluster_sizes[h];
		double mu_hat_ih = (mu_star_ih / sigma_star_h(i, h) + N_ih -
		                    0.5 * samplesPerGroup[i] + omega_ih * C_ih) *
							(sigma_hat_ih);
		transformed_weights(i, h) = normal_rng(mu_hat_ih, std::sqrt(sigma_hat_ih), rng);
    }

    weights.row(i) = utils::InvAlr(static_cast<Eigen::VectorXd>(transformed_weights.row(i)), true);
  }

  /*#pragma omp parallel for
  for (int i = 0; i < numGroups; i++)
    transformed_weights.row(i) = utils::Alr(weights.row(i), true);
}*/

// We use a MH step with a truncated normal proposal
/*void SpatialMixtureSampler::sampleRho() {
  iter += 1;
  double curr = rho;
  double sigma;
  if (iter < 3) {
    sigma = 0.01;
  } else {
    if (sigma_n_rho == 0)
      sigma = 0.01;
    else {
      if (stan::math::uniform_rng(0.0, 1.0, rng) < 0.05)
        sigma = 0.01;
      else
        sigma = 2.38 * sigma_n_rho;
    }
  }
  double proposed = utils::trunc_normal_rng(curr, sigma, 0.0, 0.9999, rng);

  // compute acceptance ratio
  Eigen::MatrixXd row_prec = F - proposed * W_init;

  Eigen::MatrixXd temp =
      utils::removeColumn(transformed_weights, numComponents - 1);
  Eigen::MatrixXd meanmat = Eigen::MatrixXd::Zero(numGroups, numComponents - 1);
  for (int i = 0; i < numGroups; i++)
    meanmat.row(i) = mtildes.row(node2comp[i]).head(numComponents - 1);

  double num =
      stan::math::beta_lpdf(proposed, alpha, beta) +
      utils::matrix_normal_prec_lpdf(temp, meanmat, row_prec, SigmaInv) +
      utils::trunc_normal_lpdf(proposed, curr, sigma, 0.0, 1);

  row_prec = F - curr * W_init;
  double den =
      stan::math::beta_lpdf(curr, alpha, beta) +
      utils::matrix_normal_prec_lpdf(temp, meanmat, row_prec, SigmaInv) +
      utils::trunc_normal_lpdf(curr, proposed, sigma, 0.0, 1);

  double arate = std::min(1.0, std::exp(num - den));
  if (stan::math::uniform_rng(0.0, 1.0, rng) < arate) {
    rho = proposed;
    numAccepted += 1;

    F = Eigen::MatrixXd::Zero(numGroups, numGroups);
    for (int i = 0; i < numGroups; i++)
      F(i, i) = rho * W_init.row(i).sum() + (1 - rho);

    F_by_comp.resize(num_connected_comps);
    G_by_comp.resize(num_connected_comps);
    for (int k = 0; k < num_connected_comps; k++) {
      Eigen::MatrixXd curr_f =
          Eigen::MatrixXd::Zero(comp2node[k].size(), comp2node[k].size());
      Eigen::MatrixXd curr_g =
          Eigen::MatrixXd::Zero(comp2node[k].size(), comp2node[k].size());
      for (int i = 0; i < comp2node[k].size(); i++) {
        curr_f(i, i) = F(comp2node[k][i], comp2node[k][i]);
        for (int j = 0; j < comp2node[k].size(); j++) {
          curr_g(i, j) = W_init(comp2node[k][i], comp2node[k][j]);
        }
      }
      F_by_comp[k] = curr_f;
      G_by_comp[k] = curr_g;
    }
  }

  // update adaptive MCMC params
  rho_sum += rho;
  rho_sum_sq += rho * rho;
  double rho_mean = rho_sum / iter;
  sigma_n_rho = rho_sum_sq / iter - rho_mean * rho_mean;
}*/

void SpatialMixtureSampler::sampleSigma() {
	Eigen::MatrixXd Vn = V0;
	double nu_n = nu + numGroups;
	Eigen::MatrixXd F_m_rhoG = F - W_init * rho;

	// #pragma omp parallel for collapse(2)
  	for (int i = 0; i < numGroups; i++) {
    	Eigen::VectorXd wtilde_i = transformed_weights.row(i).head(numComponents - 1);
    	Eigen::VectorXd mtilde_i = mtildes.row(node2comp[i]).head(numComponents - 1);
	    for (int j = 0; j < numGroups; j++) {
			Eigen::VectorXd wtilde_j = transformed_weights.row(j).head(numComponents - 1);
			Eigen::VectorXd mtilde_j = mtildes.row(node2comp[j]).head(numComponents - 1);
			Vn += ((wtilde_i - mtilde_i) * (wtilde_j - mtilde_j).transpose()) * F_m_rhoG(i, j);
	    }
	}
	Sigma = inv_wishart_rng(nu_n, Vn, rng);
	_computeInvSigmaH();
}

/*void SpatialMixtureSampler::regress() {
	// Compute mu and v
	int start = 0;
	int s = 0;

	for (int i = 0; i < numGroups; i++) {
		#pragma omp parallel for
    	for (int j = 0; j < samplesPerGroup[i]; j++) {
      		s = cluster_allocs[i][j];
			mu(start + j) = means[s];
			V.diagonal()[start + j] = 1.0 / (stddevs[s] * stddevs[s]);
    	}
		start += samplesPerGroup[i];
	}

	// compute posterior parameters for beta
	Eigen::VectorXd postMean(p_size);
	Eigen::MatrixXd postPrec(p_size, p_size);
	postPrec = predictors.transpose() * V * predictors + reg_coeff_prec;
	postMean = postPrec.ldlt().solve(predictors.transpose() * V * (reg_data - mu));
	reg_coeff = stan::math::multi_normal_prec_rng(postMean, postPrec, rng);
}

void SpatialMixtureSampler::computeRegressionResiduals() {
	Eigen::VectorXd residuals = reg_data - predictors * reg_coeff;
	int start = 0;

	for (int i = 0; i < numGroups; i++) {
	
		#pragma omp parallel for
		for (int j = 0; j < samplesPerGroup[i]; j++) {
			data[i][j] = residuals[start + j];
    	}
    	
    	start += samplesPerGroup[i];
	}
}

void SpatialMixtureSampler::_computeInvSigmaH() {
  SigmaInv = Sigma.llt().solve(
      Eigen::MatrixXd::Identity(numComponents - 1, numComponents - 1));

  Eigen::MatrixXd I =
      Eigen::MatrixXd::Identity(numComponents - 2, numComponents - 2);

  // #pragma omp parallel for
  for (int h = 0; h < numComponents - 1; ++h) {
    pippo[h] = utils::removeColumn(Sigma, h).row(h) *
               utils::removeRowColumn(Sigma, h).llt().solve(I);
  }

  // #pragma omp parallel for
  for (int h = 0; h < numComponents - 1; ++h) {
    double aux = pippo[h].dot(utils::removeRow(Sigma, h).col(h));
    for (int i = 0; i < numGroups; i++) {
      sigma_star_h(i, h) = (Sigma(h, h) - aux) / F(i, i);
    }
  }
}*/

void SpatialMixtureSampler::sample_mtilde() {
  int H = numComponents;
  Eigen::MatrixXd prec_prior =
      Eigen::MatrixXd::Identity(numComponents - 1, numComponents - 1).array() *
      (1.0 / mtilde_sigmasq);

  Eigen::MatrixXd F_min_rhoG = F - rho * W_init;

  for (int k = 1; k < num_connected_comps; k++) {
    Eigen::VectorXd currweights(comp2node[k].size() * (numComponents - 1));
    for (int i = 0; i < comp2node[k].size(); i++) {
      currweights.segment(i * (H - 1), (H - 1)) =
          transformed_weights.row(comp2node[k][i]).head(H - 1).transpose();
    }

    Eigen::MatrixXd curr_f_min_rhoG = F_by_comp[k] - rho * G_by_comp[k];

    Eigen::MatrixXd curr_prec = kroneckerProduct(curr_f_min_rhoG, SigmaInv);

    Eigen::MatrixXd I_star =
        kroneckerProduct(Eigen::VectorXd::Ones(comp2node[k].size()),
                         Eigen::MatrixXd::Identity((H - 1), (H - 1)));

    Eigen::MatrixXd prec_post =
        I_star.transpose() * curr_prec * I_star + prec_prior;
    Eigen::VectorXd m_post =
        prec_post.ldlt().solve(I_star.transpose() * curr_prec * currweights);

    Eigen::VectorXd sampled =
        stan::math::multi_normal_prec_rng(m_post, prec_post, rng);

    mtildes.row(k).head(H - 1) = sampled;
  }
}

/*void SpatialMixtureSampler::saveState(Collector<UnivariateState> *collector) {
  collector->collect(getStateAsProto());
}

UnivariateState SpatialMixtureSampler::getStateAsProto() {
  UnivariateState state;
  state.set_num_components(numComponents);
  for (int i = 0; i < numGroups; i++) {
    UnivariateState::GroupParams *p;
    p = state.add_groupparams();
    Eigen::VectorXd w = weights.row(i);

    *p->mutable_weights() = {w.data(), w.data() + numComponents};
    *p->mutable_cluster_allocs() = {cluster_allocs[i].begin(),
                                    cluster_allocs[i].end()};
  }
  for (int h = 0; h < numComponents; h++) {
    UnivariateMixtureAtom *atom;
    atom = state.add_atoms();
    atom->set_mean(means[h]);
    atom->set_stdev(stddevs[h]);
  }
  state.set_rho(rho);
  state.mutable_sigma()->set_rows(Sigma.rows());
  state.mutable_sigma()->set_cols(Sigma.cols());
  *state.mutable_sigma()->mutable_data() = {Sigma.data(),
                                            Sigma.data() + Sigma.size()};

  if (regression)
    *state.mutable_regression_coefficients() = {reg_coeff.data(),
                                                reg_coeff.data() + p_size};
  return state;
}

void SpatialMixtureSampler::printDebugString() {
  Rcpp::Rcout << "***** Debug String ****" << std::endl;
  Rcpp::Rcout << "numGroups: " << numGroups << ", samplesPerGroup: ";
  for (int n : samplesPerGroup) Rcpp::Rcout << n << ", ";
  Rcpp::Rcout << std::endl;

  // one vector per component
  std::vector<std::vector<std::vector<double>>> datavecs(numGroups);

  for (int i = 0; i < numGroups; i++) {
    datavecs[i].resize(numComponents);
    for (int j = 0; j < samplesPerGroup[i]; j++) {
      int comp = cluster_allocs[i][j];
      datavecs[i][comp].push_back(data[i][j]);
    }
  }

  Rcpp::Rcout << "Sigma: \n" << Sigma << std::endl << std::endl;

  for (int h = 0; h < numComponents; h++) {
    Rcpp::Rcout << "### Component #" << h << std::endl;
    Rcpp::Rcout << "##### Atom: mean=" << means[h] << ", sd=" << stddevs[h]
              << " weights per group: " << weights.col(h).transpose()
              << std::endl;
    Rcpp::Rcout << std::endl;
  }

  if (regression) {
    Rcpp::Rcout << "Regression Coefficients: " << std::endl;
    Rcpp::Rcout << "    " << reg_coeff.transpose() << std::endl;
  }
}*/
