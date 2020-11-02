// Template function definitions

template<typename T>
T spmixLogLikelihood::operator() (const Eigen::Matrix<T, Eigen::Dynamic, 1> & x) const {

	if (x.size() != numGroups + 2)
		throw std::runtime_error("Input vector is not of the proper size.");

	// Creating buffers
	int numComponents_ext(numComponents);
	Eigen::Matrix<T, Eigen::Dynamic, 1> transformed_weights_ext;
	Eigen::Matrix<T, Eigen::Dynamic, 1> means_ext;
	Eigen::Matrix<T, Eigen::Dynamic, 1> sqrt_std_devs_ext;
	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> Sigma_ext;

	// Splitting input vector
	Eigen::Matrix<T, Eigen::Dynamic, 1> weights_toadd(numGroups);
	weights_toadd << x.head(numGroups);
	std::copy(x.data(), x.data()+numGroups, weights_toadd.data());
	T means_toadd(x(numGroups));
	T std_devs_toadd(x(numGroups+1));

	// DEBUG
	/*Rcpp::Rcout << "x: " << x.transpose() << "\n"
	<< "weights_toadd: " << weights_toadd.transpose() << "\n"
	<< "means_toadd: " << means_toadd << "\n"
	<< "std_devs_toadd: " << std_devs_toadd << std::endl << std::endl;

	Rcpp::Rcout << "BEFORE EXPANSION: " << std::endl;
	Rcpp::Rcout << "transformed_weights_vect: " << transformed_weights_vect.transpose() << "\n"
	<< "means: " << means.transpose() << "\n"
	<< "sqrt_std_devs: " << sqrt_std_devs.transpose() << "\n"
	<< "Sigma:\n" << Sigma << std::endl << std::endl;*/

	// Promoting and extending variables
	numComponents_ext++;
	transformed_weights_ext.resize(transformed_weights_vect.size() + weights_toadd.size());
	transformed_weights_ext << transformed_weights_vect, weights_toadd;
	means_ext.resize(means.size() + 1);
	means_ext << means, means_toadd;
	sqrt_std_devs_ext.resize(sqrt_std_devs.size() + 1);
	sqrt_std_devs_ext << sqrt_std_devs, std_devs_toadd;
	Sigma_ext = Eigen::MatrixXd::Zero(numComponents_ext - 1, numComponents_ext - 1);
	Sigma_ext.block(0, 0, numComponents - 1, numComponents - 1) = Sigma;
	Sigma_ext(numComponents_ext - 2, numComponents_ext - 2) = Sigma(0,0);

	// DEBUG
	/*Rcpp::Rcout << "AFTER EXPANSION: " << std::endl;
	Rcpp::Rcout << "transformed_weights_vect: " << transformed_weights_ext.transpose() << "\n"
	<< "means: " << means_ext.transpose() << "\n"
	<< "sqrt_std_devs: " << sqrt_std_devs_ext.transpose() << "\n"
	<< "Sigma:\n" << Sigma_ext << std::endl << std::endl;*/

	// Computation
	T output{0.};

	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> transformed_weights(numGroups, numComponents_ext - 1);
	transformed_weights << transformed_weights_ext;
	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> weights(numGroups, numComponents_ext);
	for (int i = 0; i < weights.rows(); ++i) {
		weights.row(i) = utils::InvAlr(static_cast<Eigen::Matrix<T,-1,1>>(transformed_weights.row(i)), false);
	}

	//DEBUG
	/*Rcpp::Rcout << "transformed_weights:\n" << transformed_weights << "\n"
	<< "weights:\n" << weights << std::endl;*/

	// Computing contribution of data (PARE CORRETTO)
	for (int i = 0; i < data.size(); ++i) {
	    for (int j = 0; j < data[i].size(); ++j) {
	        std::vector<T> contributions(numComponents_ext);
	        for (int h = 0; h < numComponents_ext; ++h) {
				contributions[h] = log(weights(i,h)) + stan::math::normal_lpdf(data[i][j],
																			   means_ext(h),
																			   sqrt_std_devs_ext(h) * sqrt_std_devs_ext(h));
	        }
	        output += stan::math::log_sum_exp(contributions);
	    }
	}

    // Contributions from kernels
    for (int h = 0; h < numComponents_ext; ++h) {
		T sigma = sqrt_std_devs_ext(h)*sqrt_std_devs_ext(h);
		T sigmasq = sigma*sigma;
		T means_stdev = sigma / std::sqrt(params.p0_params().lam_());
		output += stan::math::inv_gamma_lpdf(sigmasq, params.p0_params().a(), params.p0_params().b()) +
                  stan::math::normal_lpdf(means_ext(h), params.p0_params().mu0(), means_stdev);
		//T std_dev = sqrt_std_devs_ext(h)*sqrt_std_devs_ext(h)*sqrt_std_devs_ext(h)*sqrt_std_devs_ext(h);
		//T tau = 1.0/(std_dev * std_dev);
		//T sigmaNorm = std_dev / std::sqrt(params.p0_params().lam_());
		//output += stan::math::gamma_lpdf(tau, params.p0_params().a(), params.p0_params().b()) +
        //          stan::math::normal_lpdf(means_ext[h], params.p0_params().mu0(), sigmaNorm);
    }

    // Contribution from weights
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> F = Eigen::MatrixXd::Zero(numGroups, numGroups);
    for (int i = 0; i < numGroups; i++)
    	F(i, i) = rho * W.row(i).sum() + (1 - rho);
    Eigen::Matrix<T, Eigen::Dynamic, 1> weightsMean = Eigen::VectorXd::Zero(transformed_weights_ext.size());
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> weightsCov = Eigen::KroneckerProduct((F - rho*W),
    	Sigma_ext.inverse()).eval().inverse();
    output += stan::math::multi_normal_lpdf(transformed_weights_ext, weightsMean, weightsCov);

    // Contribution from other stuff if needed (rho, m_tilde, H, Sigma)*/
	return output;
}

template<typename T>
T test_function::operator() (const Eigen::Matrix<T, Eigen::Dynamic, 1> & x) const {
	return -(x(0)-4)*(x(0)-4)*(x(0)+5)*(x(0)+5) - (x(1)-4)*(x(1)-4)*(x(1)+5)*(x(1)+5);
}
