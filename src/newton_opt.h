#ifndef NEWTON_OPT_H
#define NEWTON_OPT_H

#include "newton_traits.h"
#include "newton_options.pb.h"

class NewtonState: public NewtonTraits {
public:
	ReturnType current_solution;
	ArgumentType current_minimizer;
	GradientType current_gradient;
	HessianType current_hessian;
	unsigned int current_iteration;
	double current_norm_res;
};


class NewtonOpt: public NewtonTraits {
public:
	NewtonOpt(const TargetFunctionType & _target_function, const NewtonOptions & _options);
	void solve(const ArgumentType & x0);
	NewtonState get_state() const;
protected:
	TargetFunctionType target_function;
	NewtonOptions options;
	NewtonState state;
};


#endif // NEWTON_OPT_H
