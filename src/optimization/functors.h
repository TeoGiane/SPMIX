#ifndef FUNCTORS_HPP
#define FUNCTORS_HPP

#define STRICT_R_HEADERS
#include <Rcpp.h>
#include <Eigen/Dense>
#include <unsupported/Eigen/KroneckerProduct>
#include <stan/math/prim/mat.hpp>

#include <exception>

#include <univariate_mixture_state.pb.h>
#include <sampler_params.pb.h>
#include <utils.h>

namespace function {

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
};


class test_function {
  public:
  	double operator()() const {return 0.;};
	template<typename T> T operator() (const Eigen::Matrix<T,Eigen::Dynamic,1> & x) const {
		return -(x(0)-4)*(x(0)-4)*(x(0)+5)*(x(0)+5) - (x(1)-4)*(x(1)-4)*(x(1)+5)*(x(1)+5);
	}
	Eigen::VectorXd init() const {return (Eigen::VectorXd(2) << -3.,-3.).finished();}
};

#include "functors_impl.h"

}; // namespace function

#endif // FUNCTORS_HPP