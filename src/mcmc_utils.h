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
	Eigen::VectorXd sqrt_std_devs;
	Eigen::VectorXd transformed_weights_vect;
	Eigen::MatrixXd Sigma;
  public:
	spmixLogLikelihood(const std::vector<std::vector<double>> & _data, const Eigen::MatrixXd & _W,
					   const SamplerParams & _params, const UnivariateState & _state);

	spmixLogLikelihood(const std::vector<std::vector<double>> & _data,
					   const Eigen::MatrixXd & _W,
					   const SamplerParams & _params, 
					   int _numGroups,
					   int _numComponents,
					   double _rho,
					   const Eigen::VectorXd & _means,
					   const Eigen::VectorXd & _sqrt_std_devs,
					   const Eigen::VectorXd & _transformed_weights_vect,
					   const Eigen::MatrixXd & _Sigma);

	double operator()() const;

	template<typename T> T operator() (const Eigen::Matrix<T, Eigen::Dynamic, 1> & x) const;

	Eigen::VectorXd init() const;

	int get_numGroups() const {return numGroups;};

	int get_numComponents() const {return numComponents;};

	Eigen::VectorXd get_means() const {return means;};

	Eigen::VectorXd get_sqrt_std_devs() const {return sqrt_std_devs;};

	Eigen::VectorXd get_transformed_weights_vect() const {return transformed_weights_vect;};
};

#include "mcmc_utils_impl.h"

}; // namespace utils

#endif  // MCMC_UTILS
