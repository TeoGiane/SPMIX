/* Template class definitions for NewtonMethod */

template<typename D>
NewtonMethod<D>::NewtonMethod(const D & _target_function, const OptimOptions & _options):
target_function(_target_function), options(_options) {};

template<typename D>
void NewtonMethod<D>::solve(const ArgumentType & x0) {

	ArgumentType x_old = x0;
	auto loss_function = [this](const auto & x){return -target_function(x);};

	// Defining aliases
	state.current_iteration = 0;

	for (int i = 0; i < options.max_iter(); ++i) {

		// Step Iteration
		++state.current_iteration;

		// Initializing buffers
		ReturnType fx;
		ArgumentType x_new;
		GradientType grad_fx;
		HessianType hess_fx;

		// Computing h direction
		stan::math::hessian(loss_function, x_old, fx, grad_fx, hess_fx);
		Eigen::VectorXd h = hess_fx.ldlt().solve(grad_fx);

		// Update state
		state.current_solution = fx;
		state.current_maximizer = x_old;
		state.current_gradient = grad_fx;
		state.current_hessian = hess_fx;
		state.current_gradient_norm = grad_fx.norm();

		// Next step evaluation point
		x_new = x_old - h;
		x_old = x_new;

		// Printing state (FOR DEBUG)
		state.print();
		Rcpp::Rcout << "h: " << h.transpose() << std::endl;
		Rcpp::Rcout << "Next step evaluation point: " << x_old.transpose() << std::endl << std::endl;

		// Convergence check
		if (state.current_gradient.norm() < options.tol())
			break;
	}

	return;
}
