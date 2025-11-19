#ifndef CONDITIONAL_POSTERIOR_NEGLPDF_TPP
#define CONDITIONAL_POSTERIOR_NEGLPDF_TPP

#include "conditional_posterior_neglpdf.h"

template <typename Scalar>
Scalar conditional_posterior_neglpdf_internal::likelihood_contribution(const Vector<Scalar>& x) const {
    // Rcpp::Rcout << "likelihood_contribution()" << std::endl;
    // Check if x is passed and set H accordingly
    bool have_x = (x.size() != 0);
    int H = have_x ? (numComponents + 1) : (numComponents);
    // std::cout << "have_x is " << std::boolalpha << have_x << std::endl;
    // std::cout << "H: " << H << std::endl;

    // Compute likelihood contribution
    Scalar output = 0;
    for (int i = 0; i < data.size(); ++i) {
        for (int j = 0; j < data[i].size(); ++j) {
            std::vector<Scalar> contributions(H);
            Vector<Scalar> tw(H), w(H);
            if (have_x) {
                tw << x(i), static_cast<Vector<Scalar>>(transformed_weights.row(i).head(numComponents));
                w = utils::InvAlr(tw, true);
                contributions[0] = log(w(0)) + stan::math::normal_lpdf(data[i][j], x(numGroups), exp(x(numGroups+1)));
                for (int h = 0; h < numComponents; ++h) {
                    contributions[h+1] = log(w(h+1)) + precomputed_lpdfs[i](j, h);
                }
            } else {
                tw << static_cast<Vector<Scalar>>(transformed_weights.row(i).head(numComponents));
                w = utils::InvAlr(tw, true);
                for (int h = 0; h < H; ++h) {
                    contributions[h] = log(w(h)) + precomputed_lpdfs[i](j, h);
                }
            }
            output += stan::math::log_sum_exp(contributions);
        }
    }
    // Return negative log-likelihood
    return -output;
}

template <typename Scalar>
Scalar conditional_posterior_neglpdf_internal::atoms_contribution(const Vector<Scalar>& x) const {
    // Rcpp::Rcout << "atoms_contribution()" << std::endl;
    Scalar output = precomputed_atoms_contrib;
    bool have_x = (x.size() != 0);
    if (have_x) {
        Scalar new_mean = x(numGroups);
        Scalar new_std_dev = exp(x(numGroups+1));
        Scalar m_sd = new_std_dev / std::sqrt(params.p0_params().lam_());
        output += stan::math::inv_gamma_lpdf(new_std_dev*new_std_dev, params.p0_params().a(), params.p0_params().b()) +
                  stan::math::normal_lpdf(new_mean, params.p0_params().mu0(), m_sd) +
                  log(2) + 2*x(numGroups+1);
    }
    return -output;
}

template <typename Scalar>
Scalar conditional_posterior_neglpdf_internal::weights_contribution(const Vector<Scalar>& x) const {
    // Rcpp::Rcout << "weights_contribution()" << std::endl;
    // Check if x is passed
    bool have_x = (x.size() != 0);
    // Compute weights contribution
    Scalar output = precomputed_weights_contrib;
    if(have_x) {
        Vector<Scalar> tw = x.head(numGroups);
        Vector<Scalar> mean = Vector<Scalar>::Zero(numGroups);
        output += stan::math::multi_normal_lpdf(tw, mean, cov_weights);
    }
    return -output;
}

template <typename Scalar>
Scalar conditional_posterior_neglpdf_internal::operator()(const Vector<Scalar>& x) const {
    return likelihood_contribution(x) + atoms_contribution(x) + weights_contribution(x);
}

template <typename Scalar>
Scalar conditional_posterior_neglpdf::operator()(const Vector<Scalar>& x) const {
    return lpdf_functor(x);
}

template <typename Scalar>
Scalar conditional_posterior_neglpdf::operator()(const Vector<Scalar>& x, Vector<Scalar>& grad) const {
    Scalar fx = 0.0;
    stan::math::gradient(lpdf_functor, x, fx, grad);
    return fx;
}

#endif
