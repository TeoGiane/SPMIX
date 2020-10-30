#ifndef NEWTON_OPT_H
#define NEWTON_OPT_H

#include "optimization_traits.h"
#include "optimization_options.pb.h"

namespace optimization {

class NewtonState: public OptimizationTraits {
  public:
	ReturnType current_solution;
	ArgumentType current_minimizer;
	GradientType current_gradient;
	HessianType current_hessian;
	unsigned int current_iteration;
	double current_gradient_norm;
	void print() const;
};

template<typename D>
class NewtonMethod: public OptimizationTraits {
  protected:
	function::functorBase<D> target_function;
	OptimOptions options;
	NewtonState state;
  public:
	NewtonMethod(const function::functorBase<D> & _target_function, const OptimOptions & _options);
	void solve(const ArgumentType & x0);
	NewtonState get_state() const {return state;};
};

#include "newton_method_impl.h"

} // namespace optimization

#endif // NEWTON_OPT_H
