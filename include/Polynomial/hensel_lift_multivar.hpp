#ifndef HENSEL_LIFT_MULTIVAR_HPP
#define HENSEL_LIFT_MULTIVAR_HPP

#include "polynomial.hpp"
#include "factor_node.hpp"

namespace polynomial {

template <class Coef, class Order>
struct HenselFactor {
    Polynomial<Coef, Order> factor;
    Polynomial<Coef, Order> cofactor;
};

template <class Coef, class Order>
HenselFactor<Coef, Order> henselLiftMultiVarStep(const Polynomial<Coef, Order>& f, Polynomial<Coef, Order>& g,
                                                 Polynomial<Coef, Order>& h, const Polynomial<Coef, Order>& s,
                                                 const Polynomial<Coef, Order>& t, int var, Coef a, int64_t precision);
template <class Coef, class Order>
void henselLiftMultiVar(const Polynomial<Coef, Order>& f, Polynomial<Coef, Order>& g, Polynomial<Coef, Order>& h,
                        int var, const Coef& a);
template <class Coef, class Order>
void henselLiftMultiVarTree(const Polynomial<Coef, Order>& f, FactorNode<Coef, Order>& node, int var, const Coef& a);

}  // namespace polynomial
#endif
