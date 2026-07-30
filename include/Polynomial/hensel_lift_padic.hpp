#ifndef HENSEL_LIFT_PADIC_HPP
#define HENSEL_LIFT_PADIC_HPP

#include "polynomial.hpp"
#include "factor_node.hpp"

namespace polynomial {

template <class Integer, class Order>
void henselPadicStep(const Polynomial<Integer, Order>& f, Polynomial<Integer, Order>& g, Polynomial<Integer, Order>& h,
                     Polynomial<Integer, Order>& s, Polynomial<Integer, Order>& t, const Integer& pk, const Integer& p);
template <class Integer, class Order>
void henselPadic(const Polynomial<Integer, Order>& f, Polynomial<Integer, Order>& g, Polynomial<Integer, Order>& h,
                 Polynomial<Integer, Order>& s, Polynomial<Integer, Order>& t, const Integer& p, const Integer& bound);
template <class Integer, class Order>
void henselPadicTree(const Polynomial<Integer, Order>& f, FactorNode<Integer, Order>& node, const Integer& pk,
                     const Integer& p);

template <class Integer, class Order>
void henselPadicStep(const Polynomial<Integer, Order>& f, Polynomial<Integer, Order>& g, Polynomial<Integer, Order>& h,
                     Polynomial<Integer, Order>& s, Polynomial<Integer, Order>& t, const Integer& pk,
                     const Integer& p) {
    using Poly = Polynomial<Integer, Order>;

    // factor correction
    Poly E = (f - g * h) / pk;
    Poly u = s * E;
    Poly v = t * E;
    DivisionResult<Integer, Order> divU = divideMod(u, {h}, p);
    DivisionResult<Integer, Order> divV = divideMod(v, {g}, p);
    Poly q = divU.Q[0];
    Poly r = divU.R;
    Poly a = divV.R;
    Poly b = r + q * g;
    g += pk * a;
    h += pk * b;

    // Bézout correction
    Poly one(Monomial(f.leadingMonomial().size()), Integer(1));
    Poly C = (one - s * g - t * h);
    C /= pk;
    Poly uu = s * C;
    Poly vv = t * C;
    DivisionResult<Integer, Order> divUU = divideMod(uu, {h}, p);
    DivisionResult<Integer, Order> divVV = divideMod(vv, {g}, p);
    Poly qq = divUU.Q[0];
    Poly rr = divUU.R;
    Poly ds = divVV.R;
    Poly dt = rr + qq * g;
    s += pk * ds;
    t += pk * dt;
}

template <class Integer, class Order>
void henselPadic(const Polynomial<Integer, Order>& f, Polynomial<Integer, Order>& g, Polynomial<Integer, Order>& h,
                 Polynomial<Integer, Order>& s, Polynomial<Integer, Order>& t, const Integer& p, const Integer& bound) {
    Integer pk = p;
    while (pk < bound) {
        henselPadicStep(f, g, h, s, t, pk);
        pk *= pk;
    }
}

template <class Integer, class Order>
void henselPadicTree(const Polynomial<Integer, Order>& f, FactorNode<Integer, Order>& node, const Integer& pk,
                     const Integer& p) {
    if (node.isLeaf()) return;

    auto& g = node.left->factor;
    auto& h = node.right->factor;

    henselPadicStep(f, g, h, node.s, node.t, pk, p);

    henselPadicTree(g, *node.left, pk, p);
    henselPadicTree(h, *node.right, pk, p);

    node.factor = g * h;
}
}  // namespace polynomial

#endif
