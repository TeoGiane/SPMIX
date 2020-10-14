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

/*double spmixLogLikelihood(const UnivariateState &state, const std::vector<std::vector<double>> &data,
						  const Eigen::MatrixXd &W, const SamplerParams & params) {

    // Exporting required data from state <- Here we can add the tracking variables stan::math::var
    int num_components{state.num_components()};
    int numGroups(data.size());
    double rho{state.rho()};

    std::vector<double> means; std::vector<double> std_devs;
    for (auto elem : state.atoms()) {
        means.emplace_back(elem.mean());
        std_devs.emplace_back(elem.stdev());
    }

    Eigen::MatrixXd weights(state.groupparams().size(), num_components);
    Eigen::Matrix<double, -1, -1, Eigen::RowMajor> transformed_weights(state.groupparams().size(), num_components-1);
    for (int i = 0; i < state.groupparams().size(); ++i) {
        for (int j = 0; j < num_components; ++j) {
            weights(i,j) = state.groupparams()[i].weights()[j];
        }
        transformed_weights.row(i) = utils::Alr(weights.row(i), false);
    }

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
    for (int h = 0; h < num_components; ++h) {
        double tau = 1.0/(std_devs[h]*std_devs[h]);
        double sigmaNorm = 1.0 / std::sqrt(tau * params.p0_params().lam_());
        output += stan::math::gamma_lpdf(tau, params.p0_params().a(), params.p0_params().b()) +
                  stan::math::normal_lpdf(means[h], params.p0_params().mu0(), sigmaNorm);
    }

    // Contribution from weights
    Eigen::VectorXd tw_vector(numGroups*(num_components-1));
    std::copy(transformed_weights.data(), transformed_weights.data()+transformed_weights.size(), tw_vector.data());
    Eigen::MatrixXd F = Eigen::MatrixXd::Zero(numGroups, numGroups);
    for (int i = 0; i < numGroups; i++)
      F(i, i) = rho * W.row(i).sum() + (1 - rho);
    Eigen::VectorXd weightsMean = Eigen::VectorXd::Zero(tw_vector.size());
    Eigen::MatrixXd weightsCov = Eigen::KroneckerProduct((F - rho*W), Sigma.inverse()).eval().inverse();
    output += stan::math::multi_normal_lpdf(tw_vector, weightsMean, weightsCov);

    // Contribution from other stuff (rho, m_tilde, H, Sigma)
    //Rcpp::Rcout << "Output: " << output << std::endl;
    return output;
}*/

} // namespace utils
