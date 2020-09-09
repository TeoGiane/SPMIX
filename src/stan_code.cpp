#include <Rcpp.h>
#include <stan/math.hpp>
#include "PolyaGamma.h"
#include <Eigen/Dense>
#include <unsupported/Eigen/KroneckerProduct>
using namespace Rcpp;


//' Simple test with stan/math C++ library
//'
//' Simply computes logN(1|2,3)
//' @export
// [[Rcpp::export]]
void stan_HelloWorld() {
	Rcout << "log normal(1 | 2, 3) = " << stan::math::normal_log(1, 2, 3) << std::endl;
	PolyaGamma pg;
	Rcout << "A PolyaGamma object has been instanciated!" << std::endl;
	Rcout << "Testing for KroneckerProduct..." << std::endl;
	Eigen::MatrixXd A(3,3); A << 1, 2, 3, 4, 5, 6, 7, 8, 9;
	Eigen::MatrixXd B(2,2); B << 4, 3, 2, 1;
	Rcout << "A:\n" << A << std::endl << std::endl;
	Rcout << "B:\n" << B << std::endl << std::endl;
	Eigen::MatrixXd C = Eigen::kroneckerProduct(A,B).eval();
	Rcout << "A x B:\n" << C << std::endl << std::endl;
  return;
}
