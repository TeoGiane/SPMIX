#ifndef NEWTON_OPT_H
#define NEWTON_OPT_H

#include <cmath>
#include <memory>
#include <utility>
#include <type_traits>

#include "optimization_traits.h"
#include "optimization_options.pb.h"

namespace optimization {

template<typename D>
class NewtonMethod: public OptimizationTraits {
  protected:
  	std::unique_ptr<function::functorBase<D>> target_function_ptr;
	OptimOptions options;
	OptimState state;
  public:
	NewtonMethod(const D & _target_function, const OptimOptions & _options);
	void solve(const ArgumentType & x0);
	OptimState get_state() const {return state;};
};

#include "newton_method_impl.h"

} // namespace optimization

#endif // NEWTON_OPT_H
