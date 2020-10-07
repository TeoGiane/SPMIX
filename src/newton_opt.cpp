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

	for (int i = 0; i < options.max_iter(); ++i) {
		stan::math::hessian(loss_function, x_old, fx, grad_fx, hess_fx);
		Eigen::VectorXd h = hess_fx.ldlt().solve(grad_fx);
		x_new = x_old - h;
		x_old = x_new;

		// Convergence check
		if (std::fabs(x_new.lpNorm<Eigen::Infinity>()) < options.tol())
			break;
	}

	return;
}

NewtonState NewtonOpt::get_state() const {
	return state;
}
