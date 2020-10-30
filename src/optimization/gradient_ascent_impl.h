/* Template class definitions for GradientAscent */

template<typename D>
GradientAscent<D>::GradientAscent(const typename function::functorBase<D> & _target_function, const OptimOptions & _options):
target_function(_target_function), options(_options) {};