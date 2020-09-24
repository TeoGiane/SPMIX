#ifndef MCMC_UTILS
#define MCMC_UTILS

#include <algorithm>
#include <cmath>
#include <numeric>
#include <vector>

#include "univariate_mixture_state.pb.h"
#include "sampler_params.pb.h"

#define STRICT_R_HEADERS
#include <Rcpp.h>

namespace utils {

std::vector<double> normalGammaUpdate(std::vector<double> data, double priorMean,
	double priorA, double priorB, double priorLambda);

double marginalLogLikeNormalGamma(double datum, double mean, double a, double b, double lambda);

void spmixLogLikelihood(const std::vector<std::vector<double>> &data, const UnivariateState &state, const SamplerParams &params);

}

#endif  // MCMC_UTILS
