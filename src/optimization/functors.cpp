#include "functors.h"

namespace function {

/* Definitions for spmixLogLikelihood */
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

    Eigen::MatrixXd transformed_weights = Eigen::Map<const Eigen::MatrixXd>(transformed_weights_vect.data(),
                                                                            numGroups, numComponents-1);
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
    	double sigmasq = sqrt_std_devs(h)*sqrt_std_devs(h)*sqrt_std_devs(h)*sqrt_std_devs(h);
		double means_stdev = std::sqrt(sigmasq / params.p0_params().lam_());
		output += stan::math::inv_gamma_lpdf(sigmasq, params.p0_params().a(), params.p0_params().b()) +
                  stan::math::normal_lpdf(means(h), params.p0_params().mu0(), means_stdev);
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

Eigen::VectorXd spmixLogLikelihood::init() const {

	// Generate the initial point for the newton solver.
	Eigen::VectorXd x0(numGroups + 2);
	Eigen::Map<const Eigen::MatrixXd> tw_mat(transformed_weights_vect.data(), numGroups, numComponents-1);
	x0 << tw_mat.rowwise().mean(), means.mean(), sqrt_std_devs.mean();
	return x0;
}


/* Definitions for test_function */
Eigen::VectorXd test_function::init() const {
    return (Eigen::VectorXd(2) << -10., -10.).finished();
}

}; // namespace function
