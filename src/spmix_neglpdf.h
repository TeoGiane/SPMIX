#ifndef SPMIX_NEGLPDF_H
#define SPMIX_NEGLPDF_H

#include <stan/math.hpp>
#include <vector>

#include "sampler_params.pb.h"


class spmix_neglpdf_internal
{
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
  SamplerParams params;
	// Matrix W;
	// Scalar rho;
	// Scalar sigma;
	int dropped_index;

  // Quantities created when constructor is called
  int numGroups;
	int numComponents;
	// Matrix F;
	// Matrix F_rhoWInv;

	template <typename Scalar> Scalar likelihood_contribution(const Vector<Scalar>& x) const
  {

    // Rcpp::Rcout << "likelihood_contribution()" << std::endl;
    // Check if x is passed
    bool have_x = (x.size() != 0);
    // std::cout << "have_x is " << std::boolalpha << have_x << std::endl;
    int H = have_x ? (numComponents + 1) : (numComponents);
    // std::cout << "H: " << H << std::endl;

    // Compute likelihood contribution
  	Scalar output = 0;
	  for (int i = 0; i < data.size(); ++i)
    {
	    for (int j = 0; j < data[i].size(); ++j)
      {
        std::vector<Scalar> contributions(H);
        Vector<Scalar> tw(H), w(H), m(H), lsd(H);
        if (have_x) {
          // std::cout << "Here" << std::endl;
          tw << x(i), static_cast<Vector<Scalar>>(transformed_weights.row(i));
          // std::cout << "tw: " << tw.transpose() << std::endl;
          m << x(numGroups), means;
          // std::cout << "m: " << m.transpose() << std::endl;
          lsd << x(numGroups+1), log_stddevs;
          // std::cout << "lsd: " << lsd.transpose() << std::endl;
        } else {
          // std::cout << "There" << std::endl;
          tw << static_cast<Vector<Scalar>>(transformed_weights.row(i));
          // std::cout << "tw: " << tw.transpose() << std::endl;
          m << means;
          // std::cout << "m: " << m.transpose() << std::endl;
          lsd << log_stddevs;
          // std::cout << "lsd: " << lsd.transpose() << std::endl;
        }
        w = utils::InvAlr(tw, true);
        // std::cout << "w: " << w.transpose() << std::endl;
        for (int h = 0; h < H; ++h)
        {
          contributions[h] = log(w(h)) + stan::math::normal_lpdf(data[i][j], m(h), exp(lsd(h)));
	      }
	      output += stan::math::log_sum_exp(contributions);
	    }
    }
	  return -output;
  };

  template <typename Scalar> Scalar atoms_contribution(const Vector<Scalar>& x) const
  {
    // Rcpp::Rcout << "atoms_contribution()" << std::endl;
    // Check if x is passed
    bool have_x = (x.size() != 0);
    int H = have_x ? (numComponents + 1) : (numComponents);
    Vector<Scalar> m(H), lsd(H);
    if (have_x) {
      m << means, x(numGroups);
      lsd << log_stddevs, x(numGroups+1);
    } else {
      m << means;
      lsd << log_stddevs;
    }

    // Compute atoms contribution
	  Scalar output = 0;
    for (int h = 0; h < H; ++h) {
		  Scalar sigma = exp(lsd(h)); //sqrt_stddevs_ext(h)*sqrt_stddevs_ext(h);
		  Scalar m_sd = sigma / std::sqrt(params.p0_params().lam_());
		  output += stan::math::inv_gamma_lpdf(sigma*sigma, params.p0_params().a(), params.p0_params().b()) +
                stan::math::normal_lpdf(m(h), params.p0_params().mu0(), m_sd);
    }
	// std::cout << "output: " << output << std::endl;
	return -output;
  };


  template <typename Scalar> Scalar weights_contribution(const Vector<Scalar>& x) const
  {
    // Rcpp::Rcout << "weights_contribution()" << std::endl;
    // Check if x is passed
    bool have_x = (x.size() != 0);
    // std::cout << "have x is " << std::boolalpha << have_x << std::endl;

    // Compute weights contribution
	  Scalar output = 0;
		Vector<Scalar> mean = Vector<Scalar>::Zero(numGroups);
		for (int i = 0; i < numComponents - 1; i++) {
      // std::cout << "i = " << i << std::endl;
			Vector<Scalar> wcol = transformed_weights.col(i);
      // std::cout << "wcol: " << wcol.transpose() << std::endl;
			output += stan::math::multi_normal_lpdf(wcol, mean, cov_weights);
		}
    if(have_x) {
      // std::cout << "here" << std::endl;
      Vector<Scalar> tw = x.head(numGroups);
      // std::cout << "tw: " << tw.transpose() << std::endl;
      output += stan::math::multi_normal_lpdf(tw, mean, cov_weights);
    }
	  return -output;
  };

 public:

	spmix_neglpdf_internal(
    const DataType& _data,
    const Eigen::MatrixXd& _transformed_weights,
    const Eigen::VectorXd& _means,
    const Eigen::VectorXd& _log_stddevs,
    const Eigen::MatrixXd& _cov_weights,
    const SamplerParams& _params,
    int _dropped_index = -1):
      data(_data),
      transformed_weights(_transformed_weights),
      means(_means),
      log_stddevs(_log_stddevs),
      cov_weights(_cov_weights),
      params(_params)
      {
        numGroups = data.size();
        numComponents = means.size();
      };

	template <typename Scalar> Scalar operator() (const Vector<Scalar>& x) const
  {
    return likelihood_contribution(x) + atoms_contribution(x) + weights_contribution(x);
  };

  double value() const
  {
    return this->operator()(Eigen::VectorXd(0));
  };
};

class spmix_neglpdf
{
private:
    // Members
    spmix_neglpdf_internal fun;
    // Type aliases
    using DataType = std::vector<std::vector<double>>;
    template <typename T> using Vector = Eigen::Matrix<T, Eigen::Dynamic, 1>;

public:
    spmix_neglpdf(
    const DataType& _data,
    const Eigen::MatrixXd& _transformed_weights,
    const Eigen::VectorXd& _means,
    const Eigen::VectorXd& _log_stddevs,
    const Eigen::MatrixXd& _cov_weights,
    const SamplerParams& _params,
    int _dropped_index = -1):
      fun(_data, _transformed_weights, _means, _log_stddevs, _cov_weights, _params, _dropped_index) {};

    double value() const
    {
      return fun.value();
    }

    template <typename Scalar> Scalar operator() (const Vector<Scalar>& x) const
    {
      return fun(x);
    }

    template <typename Scalar> Scalar operator() (const Vector<Scalar>& x, Vector<Scalar>& grad) const
    {
        Scalar fx = 0.0;
        stan::math::gradient(fun, x, fx, grad);
        return fx;
    }
};

#endif
