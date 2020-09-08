#include <Rcpp.h>
#include "utils.hpp"

using namespace Rcpp;

/*TODO: finire documentazione*/

//' Additive log ratio
//'
//' This utility computes the additive log ratio transorm of a given vector. This transformation is
//' defined as: (METTI DEF)
//' @param x Vector of double in \ifelse{html}{\out{R<sup>H</sup>}}{\eqn{R^{H}}}.
//' @param pad_zero Boolean, default to FALSE.
//' @return Vector of double in \ifelse{html}{\out{R<sup>H-1</sup>}}{\eqn{R^{H-1}}} i.e. alr(x).
//' @export
// [[Rcpp::export]]
Eigen::VectorXd alr(Eigen::VectorXd x, bool pad_zero = false) {
	return utils::Alr(x, pad_zero);
}

/*TODO: finire documentazione*/

//' Inverse additive log ratio
//'
//' This utility computes the inverse additive log ratio transorm of a given vector. This transofrmation is
//' defined as: (METTI DEF)
//' @param x Vector of double in \ifelse{html}{\out{R<sup>H-1</sup>}}{\eqn{R^{H-1}}}.
//' @param padded_zero Boolean, default to FALSE.
//' @return Vector of double in \ifelse{html}{\out{R<sup>H</sup>}}{\eqn{R^{H}}} i.e.
//' \ifelse{html}{\out{alr <sup>-1</sup>(x)}}{\eqn{alr^{-1}(x)}}.
//' @export
// [[Rcpp::export]]
Eigen::VectorXd inv_alr(Eigen::VectorXd x, bool padded_zero = false) {
  return utils::InvAlr(x, padded_zero);
}
