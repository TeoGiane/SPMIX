#ifndef GRADIENT_ASCENT_HH
#define GRADIENT_ASCENT_HH

#include "optimization_traits.h"
#include "optimization_options.pb.h"

namespace optimization {

class GradientState : public OptimizationTraits {
  public:
  	ReturnType current_solution;
	ArgumentType current_minimizer;
	GradientType current_gradient;
	unsigned int current_iteration;
	double current_gradient_norm;
	void print() const;
};

template<typename D>
class GradientAscent : public OptimizationTraits {
  private:
  	function::functorBase<D> target_function;
	OptimOptions options;
	GradientState state;
public:
	GradientAscent(const function::functorBase<D> & _target_function, const OptimOptions & _options);
	void solve(const ArgumentType & x0); // {return;};
	GradientState get_state() const {return state;};
};

#include "gradient_ascent_impl.h"

}; // namespace optimization

#endif // GRADIENT_ASCENT_HH