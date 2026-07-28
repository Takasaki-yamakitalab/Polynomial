#ifndef HENSEL_LIFT_PADIC_HPP
#define HENSEL_LIFT_PADIC_HPP

#include "polynomial.hpp"
#include "factor_node.hpp"

namespace polynomial {

template <class Coef, class Order>
void henselPadicStep(const Polynomial<Coef, Order>& f, Polynomial<Coef, Order>& g, Polynomial<Coef, Order>& h,
                     Polynomial<Coef, Order>& s, Polynomial<Coef, Order>& t, const Coef& pk);

template <class Coef, class Order>
void henselPadic(const Polynomial<Coef, Order>& f, Polynomial<Coef, Order>& g, Polynomial<Coef, Order>& h,
                 Polynomial<Coef, Order>& s, Polynomial<Coef, Order>& t, const Coef& p, const Coef& bound);

template <class Coef, class Order>
void henselPadicTree(const Polynomial<Coef, Order>& f, FactorNode<Coef, Order>& node, const Coef& pk);

}

#endif
