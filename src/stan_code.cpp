#include <Rcpp.h>
#include <stan/math.hpp>
using namespace Rcpp;


//' Simple test with stan/math C++ library
//'
//' Simply computes logN(1|2,3)
//' @export
// [[Rcpp::export]]
void stan_HelloWorld() {
  Rcout << "log normal(1 | 2, 3) = "
        << stan::math::normal_log(1, 2, 3) << std::endl;
  return;
}
