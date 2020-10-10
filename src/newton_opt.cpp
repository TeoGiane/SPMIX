#include "newton_opt.h"

NewtonOpt::NewtonOpt(const TargetFunctionType & _target_function, const NewtonOptions & _options):
target_function(_target_function), options(_options) {};

void NewtonOpt::solve(const ArgumentType & x0) {

	ArgumentType x_old = x0;
	auto loss_function = [this](const auto & x){return -target_function(x);};

	// Defining aliases
	ReturnType & fx = state.current_solution;
	ArgumentType & x_new = state.current_minimizer;
	GradientType & grad_fx = state.current_gradient;
	HessianType & hess_fx = state.current_hessian;
	state.current_iteration = 0;

	for (int i = 0; i < options.max_iter(); ++i) {
		state.current_iteration++;

		//Printing
		/*Rcpp::Rcout << "iter: " << state.current_iteration << "\n"
		<< "solution: " << state.current_solution << "\n"
		<< "minimizer: " << state.current_minimizer.transpose() << "\n"
		<< "x_old: " << x_old.transpose() << "\n";*/


		stan::math::hessian(loss_function, x_old, fx, grad_fx, hess_fx);
		Eigen::VectorXd h = hess_fx.fullPivLu().solve(grad_fx);
		x_new = x_old - h;

		// Check positivity of last element
		if (x_new(x_new.size()-1) < 0) x_new(x_new.size()-1) = 1e-6;
		state.current_norm_res = h.lpNorm<1>();
		x_old = x_new;

		// Printing
		/*Rcpp::Rcout << "iter: " << state.current_iteration << "\n"
		<< "solution: " << state.current_solution << "\n"
		<< "minimizer: " << state.current_minimizer.transpose() << "\n"
		<< "x_old: " << x_old.transpose() << "\n"
		<< "norm_res: " << state.current_norm_res << "\n\n";*/

		// Convergence check
		if (state.current_norm_res < options.tol())
			break;
	}

	return;
}

NewtonState NewtonOpt::get_state() const {
	return state;
}
