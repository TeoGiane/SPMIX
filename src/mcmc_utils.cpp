#include "mcmc_utils.h"

//#include <Eigen/Dense>
#include <cmath>
#include <stan/math.hpp>

namespace utils {

std::vector<double> normalGammaUpdate(
        std::vector<double> data, double priorMean, double priorA,
        double priorB, double priorLambda) {

    double postMean, postA, postB, postLambda;
    int n = data.size();
    if (n == 0) {
      return std::vector<double>{priorMean, priorA, priorB, priorLambda};
    }

    double sum = std::accumulate(std::begin(data), std::end(data), 0.0);
    double ybar = sum / n;
    postMean = (priorLambda * priorMean + sum) / (priorLambda + n);
    postA = 1.0 * priorA + 1.0 * n / 2;

    double ss = 0.0;
    std::for_each(data.begin(), data.end(), [&ss, &ybar](double x) {
      ss += (x - ybar) * (x - ybar);});

    postB = (
        priorB + 0.5 * ss +
        0.5 * priorLambda / (n + priorLambda) * n *(ybar - priorMean) * (ybar - priorMean));

    postLambda = priorLambda + n;

    return std::vector<double>{postMean, postA, postB, postLambda};
}

double marginalLogLikeNormalGamma(
        double datum, double mean, double a, double b, double lambda) {

    std::vector<double> params = normalGammaUpdate(
        std::vector<double>{datum}, mean, a, b, lambda);

    double out = std::lgamma(params[1]) - std::lgamma(a);
    out += a * std::log(b) - params[1] * std::log(params[2]);
    out += 0.5 * (std::log(lambda) - std::log(params[3]));
    out -= M_PI;
    return out;
}

void spmixLogLikelihood(const std::vector<std::vector<double>> & data, const UnivariateState & state, const SamplerParams & params) {

    // Exporting required data from state
    int num_components{state.num_components()};
    std::vector<double> means; std::vector<double> std_devs;
    for (auto elem : state.atoms()) {
        means.emplace_back(elem.mean());
        std_devs.emplace_back(elem.stdev());
    }

    Eigen::MatrixXd weights(state.groupparams().size(), num_components);
    for (int i = 0; i < state.groupparams().size(); ++i) {
        for (int j = 0; j < num_components; ++j) {
            weights(i,j) = state.groupparams()[i].weights()[j];
        }
    }
    double rho{state.rho()};
    Eigen::MatrixXd Sigma(state.sigma().rows(), state.sigma().cols());
    for (int i = 0; i < state.sigma().rows(); ++i) {
        for (int j = 0; j < state.sigma().cols(); ++j) {
            Sigma(i,j) = state.sigma().data()[i*state.sigma().rows()+j];
        }
    }

    // Initialize output
    double output{0.};

    // Computing contribution of data (Mettiamo un po' di openmp)
    for (int i = 0; i < data.size(); ++i) {
        for (int j = 0; j < data[i].size(); ++j) {
            std::vector<double> contributions(num_components);
            for (int h = 0; h < num_components; ++h) {
                contributions[h] = log(weights(i,h)) + stan::math::normal_lpdf(data[i][j], means[h], std_devs[h]);
            }
            output += stan::math::log_sum_exp(contributions);
        }
    }


    // Contributions from kernels
    for (int h = 0; i < num_components; ++h) {
        double tau = 1.0/(std_devs[h]*std_devs[h]);
        double sigmaNorm = 1.0 / std::sqrt(tau * params.normalgammaparams().lam_())
        output += stan::math::gamma_lpdf(tau, params.normalgammaparams().a(), params.normalgammaparams().b()) +
        stan::math::normal_lpdf(means[h], params.normalgammaparams().mu0(), sigmaNorm);
    }


    // COntribution from weights

    Rcpp::Rcout << "Output: " << output << std::endl;
    return;
}

}
