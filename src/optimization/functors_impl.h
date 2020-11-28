// Template function definitions

template<typename T>
T spmixLogLikelihood::operator() (const Eigen::Matrix<T, Eigen::Dynamic, 1> & x) const {

	if (x.size() != numGroups + 2)
		throw std::runtime_error("Input vector is not of the proper size.");

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

	// DEBUG
	/*Rcpp::Rcout << "x: " << x.transpose() << "\n"
	<< "weights_toadd: " << weights_toadd.transpose() << "\n"
	<< "means_toadd: " << means_toadd << "\n"
	<< "sqrt_std_devs_toadd: " << sqrt_stddevs_toadd << std::endl << std::endl;

	Rcpp::Rcout << "BEFORE EXPANSION: " << std::endl;
	Rcpp::Rcout << "transformed_weights:\n" << transformed_weights << "\n"
	<< "means: " << means.transpose() << "\n"
	<< "sqrt_std_devs: " << sqrt_stddevs.transpose() << "\n"
	<< "postNormalGammaParams:\n" << postNormalGammaParams << "\n"
	<< "Sigma:\n" << Sigma << std::endl << std::endl;*/

	// Promoting and extending variables
	++numComponents_ext;
	transformed_weights_ext.resize(numGroups,numComponents_ext-1);
	means_ext.resize(numComponents_ext);
	sqrt_stddevs_ext.resize(numComponents_ext);
	if (dropped_index == -1){
		//Rcpp::Rcout << "dropped_index == -1 case in promotion" << std::endl;
		transformed_weights_ext << transformed_weights, weights_toadd;
		means_ext << means.head(numComponents-1), means_toadd, means.tail(1);
		sqrt_stddevs_ext << sqrt_stddevs.head(numComponents-1), sqrt_stddevs_toadd, sqrt_stddevs.tail(1);

		postNormalGammaParams_ext.resize(numComponents_ext,4);
		postNormalGammaParams_ext.topRows(numComponents-1) = postNormalGammaParams.topRows(numComponents-1);
		postNormalGammaParams_ext.row(numComponents-1) << params.p0_params().mu0(),params.p0_params().a(),
														  params.p0_params().b(),params.p0_params().lam_();
		postNormalGammaParams_ext.bottomRows(1) = postNormalGammaParams.bottomRows(1);
	}
	else {
		//Rcpp::Rcout << "dropped_index != -1 case in promotion" << std::endl;
		transformed_weights_ext.leftCols(dropped_index) = transformed_weights.leftCols(dropped_index);
		transformed_weights_ext.col(dropped_index) = weights_toadd;
		transformed_weights_ext.rightCols(numComponents-1-dropped_index) = transformed_weights.rightCols(numComponents-1-dropped_index);
		means_ext << means.head(dropped_index), means_toadd, means.tail(numComponents-dropped_index);
		sqrt_stddevs_ext << sqrt_stddevs.head(dropped_index), sqrt_stddevs_toadd,sqrt_stddevs.tail(numComponents-dropped_index);
	}

	if (numComponents > 1) {
		Sigma_ext = Eigen::MatrixXd::Zero(numComponents_ext - 1, numComponents_ext - 1);
		Sigma_ext.block(0, 0, numComponents - 1, numComponents - 1) = Sigma;
		Sigma_ext(numComponents_ext - 2, numComponents_ext - 2) = Sigma(0,0);
	}

	// DEBUG
	/*Rcpp::Rcout << "AFTER EXPANSION: " << std::endl;
	Rcpp::Rcout << "transformed_weights:\n" << transformed_weights_ext << "\n"
	<< "means: " << means_ext.transpose() << "\n"
	<< "sqrt_std_devs: " << sqrt_stddevs_ext.transpose() << "\n"
	<< "postNormalGammaParams:\n" << postNormalGammaParams_ext << "\n"
	<< "Sigma:\n" << Sigma_ext << std::endl << std::endl;
	Rcpp::Rcout << "promotion ok" << std::endl;*/

	// Computation
	T output{0.};

	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> weights(numGroups, numComponents_ext);
	for (int i = 0; i < numGroups; ++i) {
		Eigen::Matrix<T,-1,1> tmp = transformed_weights_ext.row(i);
		weights.row(i) = utils::InvAlr(tmp, false);
	}
	//Rcpp::Rcout << "weights:\n" << weights << std::endl;

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
	//Rcpp::Rcout << "data contrib: " << output << std::endl;
	//Rcpp::Rcout << "data ok" << std::endl;

    // Contributions from kernels
    for (int h = 0; h < numComponents_ext; ++h) {
		T sigma = sqrt_stddevs_ext(h)*sqrt_stddevs_ext(h);
		T sigmasq = sigma*sigma;
		T means_stdev;
		if (dropped_index == -1) {
			//Rcpp::Rcout << "dropped_index == -1 case in kernels" << std::endl;
			means_stdev = sigma / std::sqrt(postNormalGammaParams_ext(h,3));
			output += stan::math::inv_gamma_lpdf(sigmasq, postNormalGammaParams_ext(h,1), postNormalGammaParams_ext(h,2)) +
                  stan::math::normal_lpdf(means_ext(h), postNormalGammaParams_ext(h,0), means_stdev);
		}
		else {
			//Rcpp::Rcout << "dropped_index != -1 case in kernels" << std::endl;
			means_stdev = sigma / std::sqrt(postNormalGammaParams(h,3));
			output += stan::math::inv_gamma_lpdf(sigmasq, postNormalGammaParams(h,1), postNormalGammaParams(h,2)) +
            	      stan::math::normal_lpdf(means_ext(h), postNormalGammaParams(h,0), means_stdev);
        }
    }

    //Rcpp::Rcout << "data + kernels contrib: " << output << std::endl;
    //Rcpp::Rcout << "kernels ok" << std::endl;

	// Contribution from weights
	if (numComponents > 1) {
		Eigen::Matrix<T,Eigen::Dynamic,1> tw_vec = Eigen::Map<Eigen::Matrix<T,-1,1>>(transformed_weights_ext.data(),
																					 transformed_weights_ext.size());
		//Rcpp::Rcout << "tw_vec: " << tw_vec.transpose() << std::endl;
		Eigen::Matrix<T,Eigen::Dynamic,1> weightsMean = Eigen::VectorXd::Zero(tw_vec.size());
		Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> F_rhoWInv = (F-rho*W).inverse();
		Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> weightsCov = Eigen::kroneckerProduct(F_rhoWInv,Sigma_ext);
		output += stan::math::multi_normal_lpdf(tw_vec, weightsMean, weightsCov);
	}
	//Rcpp::Rcout << "weights ok" << std::endl;

	//Rcpp::Rcout << "data + kernels + weights contrib: " << output << std::endl;
	// Contribution from other stuff if needed (rho, m_tilde, H, Sigma)
	return output;
}

template<typename T>
T test_function::operator() (const Eigen::Matrix<T, Eigen::Dynamic, 1> & x) const {
	return -(x(0)-4)*(x(0)-4)*(x(0)+5)*(x(0)+5) - (x(1)-4)*(x(1)-4)*(x(1)+5)*(x(1)+5);
}
