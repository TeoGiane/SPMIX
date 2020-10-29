#ifndef NEWTON_OPT_H
#define NEWTON_OPT_H

#include "optimization_traits.h"
#include "newton_options.pb.h"

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


class NewtonMethod: public OptimizationTraits {
  protected:
	TargetFunctionType target_function;
	NewtonOptions options;
	NewtonState state;
  public:
	NewtonMethod(const TargetFunctionType & _target_function, const NewtonOptions & _options);
	void solve(const ArgumentType & x0);
	NewtonState get_state() const {return state;};
};

} // namespace optimization

#endif // NEWTON_OPT_H
