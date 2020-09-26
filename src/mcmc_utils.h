#ifndef MCMC_UTILS
#define MCMC_UTILS

#include <algorithm>
#include <cmath>
#include <numeric>
#include <vector>
#include <cmath>
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

double spmixLogLikelihood(const UnivariateState &state, const std::vector<std::vector<double>> &data,
						  const Eigen::MatrixXd &W, const SamplerParams &params);

}

#endif  // MCMC_UTILS
