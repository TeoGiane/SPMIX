// Template function definitions

template<typename T>
T spmixLogLikelihood::operator() (const Eigen::Matrix<T, Eigen::Dynamic, 1> & x) const {

	if (x.size() != numGroups + 2)
		throw std::runtime_error("Input vector is not of the proper size.");

	// Setting the dropped index for all cases
	int dropped = (dropped_index == -1) ? (numComponents-1) : (dropped_index);

	// Creating buffers
	int numComponents_ext(numComponents);
	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> transformed_weights_ext;
	Eigen::Matrix<T, Eigen::Dynamic, 1> means_ext;
	Eigen::Matrix<T, Eigen::Dynamic, 1> sqrt_stddevs_ext;
	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> Sigma_ext;
	Eigen::MatrixXd postNormalGammaParams_ext;

	// Splitting input vector
	Eigen::Matrix<T, Eigen::Dynamic, 1> weights_toadd(numGroups);
	weights_toadd << x.head(numGroups);
	T means_toadd(x(numGroups));
	T sqrt_stddevs_toadd(x(numGroups+1));

	// Promoting and extending variables
	++numComponents_ext;
	transformed_weights_ext.resize(numGroups,numComponents_ext-1);
	means_ext.resize(numComponents_ext);
	sqrt_stddevs_ext.resize(numComponents_ext);
	transformed_weights_ext.leftCols(dropped) = transformed_weights.leftCols(dropped);
	transformed_weights_ext.col(dropped) = weights_toadd;
	transformed_weights_ext.rightCols(numComponents-1-dropped) = transformed_weights.rightCols(numComponents-1-dropped);
	means_ext << means.head(dropped), means_toadd, means.tail(numComponents-dropped);
	sqrt_stddevs_ext << sqrt_stddevs.head(dropped), sqrt_stddevs_toadd,sqrt_stddevs.tail(numComponents-dropped);
	if (numComponents > 1) {
		Sigma_ext = Eigen::MatrixXd::Zero(numComponents_ext - 1, numComponents_ext - 1);
		Sigma_ext.block(0, 0, numComponents - 1, numComponents - 1) = Sigma;
		Sigma_ext(numComponents_ext - 2, numComponents_ext - 2) = Sigma(0,0);
	}

	// Computation
	T output{0.};

	// Compute weights matrix
	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> weights(numGroups, numComponents_ext);
	for (int i = 0; i < numGroups; ++i) {
		Eigen::Matrix<T,-1,1> tmp = transformed_weights_ext.row(i);
		weights.row(i) = utils::InvAlr(tmp, false);
	}

	// Computing contribution of data
	for (int i = 0; i < data.size(); ++i) {
	    for (int j = 0; j < data[i].size(); ++j) {
	        std::vector<T> contributions(numComponents_ext);
	        for (int h = 0; h < numComponents_ext; ++h) {
				contributions[h] = log(weights(i,h)) + stan::math::normal_lpdf(data[i][j],
																			   means_ext(h),
																			   sqrt_stddevs_ext(h) * sqrt_stddevs_ext(h));
	        }
	        output += stan::math::log_sum_exp(contributions);
	    }
	}

    // Contributions from kernels
    for (int h = 0; h < numComponents_ext; ++h) {
		T sigma = sqrt_stddevs_ext(h)*sqrt_stddevs_ext(h);
		T sigmasq = sigma*sigma;
		T means_stdev = sigma / std::sqrt(params.p0_params().lam_());
		output += stan::math::inv_gamma_lpdf(sigmasq, params.p0_params().a(), params.p0_params().b()) +
				  stan::math::normal_lpdf(means_ext(h), params.p0_params().mu0(), means_stdev);
    }

	// Contribution from transformed weights
	if (numComponents > 1) {
		Eigen::Matrix<T,Eigen::Dynamic,1> tw_vec = Eigen::Map<Eigen::Matrix<T,-1,1>>(transformed_weights_ext.data(),
																					 transformed_weights_ext.size());
		Eigen::Matrix<T,Eigen::Dynamic,1> weightsMean = Eigen::VectorXd::Zero(tw_vec.size());
		Eigen::MatrixXd I = Eigen::MatrixXd::Identity(numGroups,numGroups);
		Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> F_rhoWInv = (F-rho*W).ldlt().solve(I);
		Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> weightsCov = Eigen::kroneckerProduct(Sigma_ext, F_rhoWInv);
		output += stan::math::multi_normal_lpdf(tw_vec, weightsMean, weightsCov);
	}

	// Contribution from other stuff if needed (rho, m_tilde, H, Sigma)
	return output;
}

template<typename T>
T test_function::operator() (const Eigen::Matrix<T, Eigen::Dynamic, 1> & x) const {
	return -(x(0)-4)*(x(0)-4)*(x(0)+5)*(x(0)+5) - (x(1)-4)*(x(1)-4)*(x(1)+5)*(x(1)+5);
}
