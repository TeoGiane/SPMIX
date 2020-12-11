#include "functors.h"

namespace function {

/* Definitions for spmixLogLikelihood */
/*spmixLogLikelihood::spmixLogLikelihood(const std::vector<std::vector<double>> & _data, const Eigen::MatrixXd & _W,
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
}*/

spmixLogLikelihood::spmixLogLikelihood(const std::vector<std::vector<double>> & _data,
                                       const Eigen::MatrixXd & _W,
                                       const SamplerParams & _params,
                                       int _numGroups,
                                       int _numComponents,
                                       double _rho,
                                       const Eigen::VectorXd & _means,
                                       const Eigen::VectorXd & _sqrt_stddevs,
                                       const Eigen::MatrixXd & _transformed_weights,
                                       const Eigen::MatrixXd & _Sigma,
                                       int _dropped_index):
data(_data), W(_W), params(_params), numGroups(_numGroups), numComponents(_numComponents), rho(_rho), means(_means),
sqrt_stddevs(_sqrt_stddevs), transformed_weights(_transformed_weights), Sigma(_Sigma), dropped_index(_dropped_index) {

	// Computing F
	F = Eigen::MatrixXd::Zero(numGroups, numGroups);
	for (int i = 0; i < numGroups; ++i)
		F(i, i) = rho*W.row(i).sum()+(1.-rho);
};

double spmixLogLikelihood::operator()() const {

    // Computation
    double output{0.};

    Eigen::MatrixXd weights(numGroups, numComponents);
    for (int i = 0; i < weights.rows(); ++i) {
        weights.row(i) = utils::InvAlr(static_cast<Eigen::VectorXd>(transformed_weights.row(i)), false);
    }

    // Computing contribution of data
    for (int i = 0; i < data.size(); ++i) {
        for (int j = 0; j < data[i].size(); ++j) {
            std::vector<double> contributions(numComponents);
            for (int h = 0; h < numComponents; ++h) {
                contributions[h] = log(weights(i,h)) + stan::math::normal_lpdf(data[i][j],
                															   means[h],
                															   sqrt_stddevs[h]*sqrt_stddevs[h]);
            }
            output += stan::math::log_sum_exp(contributions);
        }
    }

    // Contributions from kernels
    for (int h = 0; h < numComponents; ++h) {
    	double sigmasq = sqrt_stddevs(h)*sqrt_stddevs(h)*sqrt_stddevs(h)*sqrt_stddevs(h);
        double means_stdev = std::sqrt(sigmasq / params.p0_params().lam_());
        output += stan::math::inv_gamma_lpdf(sigmasq, params.p0_params().a(), params.p0_params().b()) +
        		  stan::math::normal_lpdf(means(h), params.p0_params().mu0(), means_stdev);
    }

    // Contribution from weights
    if (numComponents > 1) {
		Eigen::VectorXd tw_vec = Eigen::Map<const Eigen::VectorXd>(transformed_weights.data(), transformed_weights.size());
		Eigen::VectorXd weightsMean = Eigen::VectorXd::Zero(tw_vec.size());
		Eigen::MatrixXd F_rhoWInv = (F-rho*W).inverse();
		Eigen::MatrixXd weightsCov = Eigen::kroneckerProduct(Sigma, F_rhoWInv);
		output += stan::math::multi_normal_lpdf(tw_vec, weightsMean, weightsCov);
    }

    // Contribution from other stuff, if needed (rho, m_tilde, H, Sigma)
    return output;
};

Eigen::VectorXd spmixLogLikelihood::init(const int & randComp) const {

	// Generate the initial point for the newton solver.
	Eigen::VectorXd x0(numGroups + 2);
    x0 << transformed_weights.col(randComp), means(randComp), sqrt_stddevs(randComp);
	//x0 << transformed_weights.rowwise().mean(), means.mean(), sqrt_stddevs.mean();
	return x0;
}


/* Definitions for test_function */
Eigen::VectorXd test_function::init() const {
    return (Eigen::VectorXd(2) << -10., -10.).finished();
}

}; // namespace function
