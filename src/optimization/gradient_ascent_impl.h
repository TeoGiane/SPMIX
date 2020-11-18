/* Template class definitions for GradientAscent */

template<typename D>
GradientAscent<D>::GradientAscent(const D & _target_function, const OptimOptions & _options):
target_function(_target_function), options(_options) {};

template<typename D>
void GradientAscent<D>::solve(const ArgumentType & x0) {

	ArgumentType x_old = x0 - Eigen::VectorXd::Constant(x0.size(),0.1);
	ArgumentType x_curr = x0;

	// DEBUG
	/*Rcpp::Rcout << "x_old: " << x_old.transpose() << "\n"
	<< "x_curr: " << x_curr.transpose() << std::endl;*/

	state.current_iteration = 0;

	for (int i = 0; i < options.max_iter(); ++i) {

		// Step Iteration
		++state.current_iteration;

		// Initializing buffers
		double gamma_i; //double gamma_i_old = 1;
		ArgumentType x_new;
		ReturnType fx_curr, fx_old;
		GradientType grad_fx_curr, grad_fx_old;
		
		// Computing gradients and step size gamma
		stan::math::gradient(target_function, x_curr, fx_curr, grad_fx_curr);
		stan::math::gradient(target_function, x_old, fx_old, grad_fx_old);
		gamma_i = std::abs((x_curr - x_old).dot(grad_fx_curr - grad_fx_old)) / (grad_fx_curr - grad_fx_old).squaredNorm();

		if (isnan(gamma_i)) {
			state.stagnated = true;
			break;
		}

		// Update state
		state.current_solution = fx_curr;
		state.current_minimizer = x_curr;
		state.current_gradient = grad_fx_curr;
		state.current_gradient_norm = grad_fx_curr.norm();

		// Next step evaluation points and update
		x_new = x_curr + gamma_i*grad_fx_curr;
		x_old = x_curr;
		x_curr = x_new;

		// DEBUG
		/*state.print();
		Rcpp::Rcout << "gamma_i: " << gamma_i << "\n"
					<< "x_old: " << x_old.transpose() << "\n"
					<< "x_curr: " << x_curr.transpose() << std::endl << std::endl;*/

		// Convergence check
		if (state.current_gradient_norm < options.tol())
			break;

		// Stagnation check
		/*if (gamma_i < 1e-16 and gamma_i_old < 1e-16) {
			state.stagnated = true;
			break;
		}

		// Save gamma_i_old
		gamma_i_old = gamma_i;*/
	}

	stan::math::hessian(target_function, state.current_minimizer,
						state.current_solution, state.current_gradient, state.current_hessian);
	//state.print();
	return;
};
