#ifndef NEWTONTRAITS_HPP
#define NEWTONTRAITS_HPP

// [[Rcpp::depends(BH)]]
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppParallel)]]
// [[Rcpp::depends(StanHeaders)]]
#define STRICT_R_HEADERS
#include <stan/math/fwd/mat.hpp>
#include <stan/math/mix/mat.hpp>
#include <stan/math.hpp>
#include <Rcpp.h>
#include <RcppEigen.h>

#include <Eigen/Dense>
// #include <functional>

#include "functors.h"

namespace optimization {

class OptimizationTraits {
  public:
	using ArgumentType = Eigen::VectorXd;
	using ReturnType = double;
	using GradientType = Eigen::VectorXd;
	using HessianType = Eigen::MatrixXd;
};

class OptimState: public OptimizationTraits {
  public:
	ReturnType current_solution;
	ArgumentType current_minimizer;
	GradientType current_gradient;
	HessianType current_hessian;
	unsigned int current_iteration;
	double current_gradient_norm;
	bool stagnated = false;
	void print() const;
};

} // namespace optimization

#endif /* NEWTONTRAITS_HPP */