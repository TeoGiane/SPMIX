#include "newton_method.h"

namespace optimization
{

void NewtonState::print() const {
	Rcpp::Rcout << "Printing State:\n"
	<< "Iteration: " 	<< current_iteration << "\n"
	<< "Solution: "		<< current_solution << "\n"
	<< "Minimizer: " 	<< current_minimizer.transpose() << "\n"
	<< "Gradient: "  	<< current_gradient.transpose() << "\n"
	<< "Hessian:\n"		<< current_hessian << "\n"
	<< "||Gradient||: " << current_gradient_norm << std::endl;
};

}; // namespace optimization