#ifndef GRADIENT_ASCENT_HH
#define GRADIENT_ASCENT_HH

#include <cmath>
#include <memory>
#include <utility>
#include <type_traits>

#include "functors.h"
#include "optimization_traits.h"
#include "optimization_options.pb.h"

namespace optimization {

template<typename D>
class GradientAscent: public OptimizationTraits {
  private:
  	std::unique_ptr<function::functorBase<D>> target_function_ptr;
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