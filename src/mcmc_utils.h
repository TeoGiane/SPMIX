#ifndef MCMC_UTILS
#define MCMC_UTILS

#include <algorithm>
#include <cmath>
#include <numeric>
#include <vector>
#include <cmath>
#include <exception>

#include <Eigen/Dense>
#include <unsupported/Eigen/KroneckerProduct>
#include <stan/math/prim/mat.hpp>

#include "univariate_mixture_state.pb.h"
#include "sampler_params.pb.h"
#include "utils.h"

#define STRICT_R_HEADERS
#include <Rcpp.h>

namespace utils {

std::vector<double> normalGammaUpdate(std::vector<double> data, double priorMean,
	double priorA, double priorB, double priorLambda);

double marginalLogLikeNormalGamma(double datum, double mean, double a, double b, double lambda);

class spmixLogLikelihood {
  private:
  	std::vector<std::vector<double>> data;
	Eigen::MatrixXd W;
	SamplerParams params;
	int numGroups;
	int numComponents;
	double rho;
	Eigen::VectorXd means;
	Eigen::VectorXd std_devs;
	Eigen::VectorXd transformed_weights_vect;
	Eigen::MatrixXd Sigma;
  public:
	spmixLogLikelihood(const std::vector<std::vector<double>> & _data, const Eigen::MatrixXd & _W,
					   const SamplerParams & _params, const UnivariateState & _state);
	template<typename T> T operator() (const Eigen::Matrix<T, Eigen::Dynamic, 1> & x) const;
};

// Template function definition
template<typename T>
T spmixLogLikelihood::operator() (const Eigen::Matrix<T, Eigen::Dynamic, 1> & x) const {

	int numComponents_ext(numComponents);
	Eigen::Matrix<T, Eigen::Dynamic, 1> transformed_weights_ext;
	Eigen::Matrix<T, Eigen::Dynamic, 1> means_ext;
	Eigen::Matrix<T, Eigen::Dynamic, 1> std_devs_ext;
	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> Sigma_ext;

	if (x.size() != 0) {

		if (x.size() != numGroups + 2)
			throw std::runtime_error("Input vector is not of size numGroups + 2.");

		// Splitting input vector
		Eigen::Matrix<T, Eigen::Dynamic, 1> weights_toadd(numGroups);
		std::copy(x.data(), x.data()+numGroups, weights_toadd.data());
		T means_toadd(x(numGroups));
		T std_devs_toadd(x(numGroups+1));

		// Promoting and extending variables
		numComponents_ext++;
		transformed_weights_ext.resize(transformed_weights_vect.size() + weights_toadd.size());
		transformed_weights_ext << transformed_weights_vect, weights_toadd;
		means_ext.resize(means.size() + 1);
		means_ext << means, means_toadd;
		std_devs_ext.resize(std_devs.size() + 1);
		std_devs_ext << std_devs, std_devs_toadd;

		// TEMPORARY SECTION NOT CORRECT BUT TO MAKE FUNCTION RUN.
		Sigma_ext = Eigen::MatrixXd::Zero(numComponents_ext - 1, numComponents_ext - 1);
		Sigma_ext.block(0, 0, numComponents - 1, numComponents - 1) = Sigma;
    	Sigma_ext(numComponents_ext - 2, numComponents_ext - 2) = 0.01;
	}
	else {
		transformed_weights_ext.resize(transformed_weights_vect.size());
		transformed_weights_ext << transformed_weights_vect;
		means_ext.resize(means.size());
		means_ext << means;
		std_devs_ext.resize(std_devs.size());
		std_devs_ext << std_devs;
		Sigma_ext.resize(Sigma.rows(), Sigma.cols()); Sigma_ext << Sigma;
	}

	// Computation
	T output{0.};

	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> transformed_weights(numGroups, numComponents_ext - 1);
	transformed_weights << transformed_weights_ext;
	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> weights(numGroups, numComponents_ext);
	for (int i = 0; i < weights.rows(); ++i) {
		weights.row(i) = utils::InvAlr(static_cast<Eigen::Matrix<T,-1,1>>(transformed_weights.row(i)), false);
	}

	// Computing contribution of data (Mettiamo un po' di openmp)
	for (int i = 0; i < data.size(); ++i) {
	    for (int j = 0; j < data[i].size(); ++j) {
	        std::vector<T> contributions(numComponents_ext);
	        for (int h = 0; h < numComponents_ext; ++h) {
	            contributions[h] = log(weights(i,h)) + stan::math::normal_lpdf(data[i][j], means_ext[h], std_devs_ext[h]);
	        }
	        output += stan::math::log_sum_exp(contributions);
	    }
	}

    // Contributions from kernels
    for (int h = 0; h < numComponents_ext; ++h) {
        T tau = 1.0/(std_devs_ext[h]*std_devs_ext[h]);
        T sigmaNorm = std_devs_ext[h] / std::sqrt(params.p0_params().lam_());
        output += stan::math::gamma_lpdf(tau, params.p0_params().a(), params.p0_params().b()) +
                  stan::math::normal_lpdf(means[h], params.p0_params().mu0(), sigmaNorm);
    }

    // Contribution from weights
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> F = Eigen::MatrixXd::Zero(numGroups, numGroups);
    for (int i = 0; i < numGroups; i++)
    	F(i, i) = rho * W.row(i).sum() + (1 - rho);
    Eigen::Matrix<T, Eigen::Dynamic, 1> weightsMean = Eigen::VectorXd::Zero(transformed_weights_ext.size());
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> weightsCov = Eigen::KroneckerProduct((F - rho*W),
    	Sigma_ext.inverse()).eval().inverse();
    output += stan::math::multi_normal_lpdf(transformed_weights_ext, weightsMean, weightsCov);

    // Contribution from other stuff (rho, m_tilde, H, Sigma)
	return output;
}


}; // namespace utils

#endif  // MCMC_UTILS
