#include "mcmc_utils.h"

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

spmixLogLikelihood::spmixLogLikelihood(const std::vector<std::vector<double>> & _data, const Eigen::MatrixXd & _W,
                                       const SamplerParams & _params, const UnivariateState & _state):
data(_data), W(_W), params(_params), numGroups(data.size()) {

    //Exporting required info from state
    numComponents = _state.num_components();
    rho = _state.rho();

    means.resize(numComponents);
    //std_devs.resize(numComponents);
    sqrt_std_devs.resize(numComponents);
    for (int i = 0; i < numComponents; ++i) {
        means(i) = _state.atoms()[i].mean();//.emplace_back(elem.mean());

        //std_devs(i) = _state.atoms()[i].stdev();
        sqrt_std_devs(i) = std::sqrt(_state.atoms()[i].stdev());
    }

    Eigen::MatrixXd transformed_weights(_state.groupparams().size(), numComponents-1);
    Eigen::MatrixXd weights(_state.groupparams().size(), numComponents);
    for (int i = 0; i < _state.groupparams().size(); ++i) {
        Eigen::VectorXd tmp(numComponents);
        for (int j = 0; j < numComponents; ++j) {
            tmp(j) = _state.groupparams()[i].weights()[j];
            weights(i,j) = _state.groupparams()[i].weights()[j];
        }
        transformed_weights.row(i) = utils::Alr(tmp, false);
    }
    transformed_weights_vect.resize(transformed_weights.size());
    transformed_weights_vect << transformed_weights;

    Sigma.resize(_state.sigma().rows(), _state.sigma().cols());
    for (int i = 0; i < _state.sigma().rows(); ++i) {
        for (int j = 0; j < _state.sigma().cols(); ++j) {
            Sigma(i,j) = _state.sigma().data()[i*_state.sigma().rows()+j];
        }
    }
}

spmixLogLikelihood::spmixLogLikelihood(const std::vector<std::vector<double>> & _data,
                                       const Eigen::MatrixXd & _W,
                                       const SamplerParams & _params,
                                       int _numGroups,
                                       int _numComponents,
                                       double _rho,
                                       const Eigen::VectorXd & _means,
                                       const Eigen::VectorXd & _sqrt_std_devs,
                                       const Eigen::VectorXd & _transformed_weights_vect,
                                       const Eigen::MatrixXd & _Sigma):
data(_data), W(_W), params(_params), numGroups(_numGroups), numComponents(_numComponents), rho(_rho), means(_means),
sqrt_std_devs(_sqrt_std_devs), transformed_weights_vect(_transformed_weights_vect), Sigma(_Sigma) {};

double spmixLogLikelihood::operator()() const {

    // Computation
    double output{0.};

    Eigen::MatrixXd transformed_weights = Eigen::Map<const Eigen::MatrixXd>(transformed_weights_vect.data(), numGroups, numComponents-1);
    Eigen::MatrixXd weights(numGroups, numComponents);
    for (int i = 0; i < weights.rows(); ++i) {
        weights.row(i) = utils::InvAlr(static_cast<Eigen::VectorXd>(transformed_weights.row(i)), false);
    }

    // Computing contribution of data (Mettiamo un po' di openmp)
    for (int i = 0; i < data.size(); ++i) {
        for (int j = 0; j < data[i].size(); ++j) {
            std::vector<double> contributions(numComponents);
            for (int h = 0; h < numComponents; ++h) {
                //contributions[h] = log(weights(i,h)) + stan::math::normal_lpdf(data[i][j], means_ext[h], std_devs_ext[h]);
                contributions[h] = log(weights(i,h)) + stan::math::normal_lpdf(data[i][j], means[h], sqrt_std_devs[h] * sqrt_std_devs[h]);
            }
            output += stan::math::log_sum_exp(contributions);
        }
    }

    // Contributions from kernels
    for (int h = 0; h < numComponents; ++h) {

        //T tau = 1.0/(std_devs_ext[h]*std_devs_ext[h]);
        double std_dev = sqrt_std_devs[h]*sqrt_std_devs[h];
        double tau = 1.0/(std_dev * std_dev);
        //T sigmaNorm = std_devs_ext[h] / std::sqrt(params.p0_params().lam_());
        double sigmaNorm = std_dev / std::sqrt(params.p0_params().lam_());

        output += stan::math::gamma_lpdf(tau, params.p0_params().a(), params.p0_params().b()) +
                  stan::math::normal_lpdf(means[h], params.p0_params().mu0(), sigmaNorm);
    }

    // Contribution from weights
    Eigen::MatrixXd F = Eigen::MatrixXd::Zero(numGroups, numGroups);
    for (int i = 0; i < numGroups; i++)
        F(i, i) = rho * W.row(i).sum() + (1 - rho);
    Eigen::VectorXd weightsMean = Eigen::VectorXd::Zero(transformed_weights.size()); //TO IMPROVE INDEED
    Eigen::MatrixXd weightsCov = Eigen::KroneckerProduct((F - rho*W), Sigma.inverse()).eval().inverse();
    output += stan::math::multi_normal_lpdf(transformed_weights_vect, weightsMean, weightsCov);

    // Contribution from other stuff, if needed (rho, m_tilde, H, Sigma)
    return output;
};

} // namespace utils
