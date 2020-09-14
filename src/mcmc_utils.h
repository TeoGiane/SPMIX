#ifndef MCMC_UTILS
#define MCMC_UTILS

#include <algorithm>
#include <cmath>
#include <numeric>
#include <vector>

namespace utils {

std::vector<double> normalGammaUpdate(
    std::vector<double> data, double priorMean, double priorA,
    double priorB, double priorLambda);


double marginalLogLikeNormalGamma(
    double datum, double mean, double a, double b, double lambda);

}
#endif  // MCMC_UTILS
