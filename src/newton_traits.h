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
#include <functional>

#include "mcmc_utils.h"

class NewtonTraits {
public:
	using ArgumentType = Eigen::VectorXd;
	using ReturnType = double;
	using GradientType = Eigen::VectorXd;
	using HessianType = Eigen::MatrixXd;
	using TargetFunctionType = utils::spmixLogLikelihood;
};

#endif /* NEWTONTRAITS_HPP */