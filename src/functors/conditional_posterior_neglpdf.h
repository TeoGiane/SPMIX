#ifndef CONDITIONAL_POSTERIOR_NEGLPDF_H
#define CONDITIONAL_POSTERIOR_NEGLPDF_H

#include <stan/math.hpp>
#include <vector>

#include "sampler_params.pb.h"
#include "utils.h"

class conditional_posterior_neglpdf_internal {
  private:

    using DataType = std::vector<std::vector<double>>;
    template <typename T> using Vector = Eigen::Matrix<T, Eigen::Dynamic, 1>;
    template <typename T> using Matrix = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;

    // Members passed at construction time
    DataType data;
    Eigen::MatrixXd transformed_weights;
    Eigen::VectorXd means;
    Eigen::VectorXd log_stddevs;
    Eigen::MatrixXd cov_weights;
    spmix::SamplerParams params;
	// Matrix W;
	// Scalar rho;
	// Scalar sigma;

    // Quantities created when constructor is called
    int numGroups;
    int numComponents;
    std::vector<Eigen::MatrixXd> precomputed_lpdfs;
    double precomputed_atoms_contrib = 0.0;
    double precomputed_weights_contrib = 0.0;
    // Matrix F;
    // Matrix F_rhoWInv;

	template <typename Scalar> Scalar likelihood_contribution(const Vector<Scalar>& x) const;
    template <typename Scalar> Scalar atoms_contribution(const Vector<Scalar>& x) const;
    template <typename Scalar> Scalar weights_contribution(const Vector<Scalar>& x) const;

 public:

	conditional_posterior_neglpdf_internal(const DataType& _data, const Eigen::MatrixXd& _transformed_weights, const Eigen::VectorXd& _means, const Eigen::VectorXd& _log_stddevs, const Eigen::MatrixXd& _cov_weights, const spmix::SamplerParams& _params):
    data(_data), transformed_weights(_transformed_weights), means(_means), log_stddevs(_log_stddevs), cov_weights(_cov_weights), params(_params) {
        
        // Compute sizes
        numGroups = data.size();
        numComponents = means.size();
        
        // Pre-computation of likelihood contributions for existing components
        precomputed_lpdfs.resize(numGroups);
        for (int i = 0; i < numGroups; ++i) {
            precomputed_lpdfs[i].resize(data[i].size(), numComponents);
            for (size_t j = 0; j < data[i].size(); ++j) {
                for (int h = 0; h < numComponents; ++h) {
                    precomputed_lpdfs[i](j, h) = stan::math::normal_lpdf(data[i][j], means(h), exp(log_stddevs(h)));
                }
            }
        }

        // Pre-computation for atoms contribution
        // precomputed_atoms_contrib = 0;
        for (int h = 0; h < numComponents; ++h) {
            double sigma = exp(log_stddevs(h));
            double m_sd = sigma / std::sqrt(params.p0_params().lam_());
            precomputed_atoms_contrib += stan::math::inv_gamma_lpdf(sigma*sigma, params.p0_params().a(), params.p0_params().b()) +
                                         stan::math::normal_lpdf(means(h), params.p0_params().mu0(), m_sd) +
                                         log(2) + 2*log_stddevs(h);
        }

        // Pre-computation for weights contribution
        // precomputed_weights_contrib = 0;
        Eigen::VectorXd mean = Eigen::VectorXd::Zero(numGroups);
        for (int i = 0; i < numComponents - 1; i++) {
            precomputed_weights_contrib += stan::math::multi_normal_lpdf(transformed_weights.col(i), mean, cov_weights);
        }
    };

	template <typename Scalar> Scalar operator() (const Vector<Scalar>& x) const;

    double value() const {
        return this->operator()(Eigen::VectorXd(0));
    };
};

class conditional_posterior_neglpdf {
  private:
    // Type aliases
    using DataType = std::vector<std::vector<double>>;
    template <typename T> using Vector = Eigen::Matrix<T, Eigen::Dynamic, 1>;
    
    // Members
    conditional_posterior_neglpdf_internal lpdf_functor;

  public:
    // Constructor
    conditional_posterior_neglpdf(const DataType& _data, const Eigen::MatrixXd& _transformed_weights, const Eigen::VectorXd& _means, const Eigen::VectorXd& _log_stddevs, const Eigen::MatrixXd& _cov_weights, const spmix::SamplerParams& _params):
    lpdf_functor(_data, _transformed_weights, _means, _log_stddevs, _cov_weights, _params) {};

    double value() const {
        return lpdf_functor.value();
    };

    template <typename Scalar> Scalar operator() (const Vector<Scalar>& x) const;
    template <typename Scalar> Scalar operator() (const Vector<Scalar>& x, Vector<Scalar>& grad) const;
};

#include "conditional_posterior_neglpdf.tpp"

#endif
