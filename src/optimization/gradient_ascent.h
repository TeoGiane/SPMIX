#ifndef GRADIENT_ASCENT_HH
#define GRADIENT_ASCENT_HH

#include <cmath>
#include "optimization_traits.h"
#include "optimization_options.pb.h"

namespace optimization {

/*class GradientState : public OptimizationTraits {
  public:
  	ReturnType current_solution;
	ArgumentType current_minimizer;
	GradientType current_gradient;
	unsigned int current_iteration;
	double current_gradient_norm;
	void print() const;
};*/

template<typename D>
class GradientAscent: public OptimizationTraits {
  private:
  	D target_function;
	OptimOptions options;
	OptimState state;
  public:
	GradientAscent(const D & _target_function, const OptimOptions & _options);
	void solve(const ArgumentType & x0);
	OptimState get_state() const {return state;};
};

#include "gradient_ascent_impl.h"

}; // namespace optimization

#endif // GRADIENT_ASCENT_HH