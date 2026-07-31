#ifndef POLYNOMIAL_HENSEL_LIFT_PADIC_HPP
#define POLYNOMIAL_HENSEL_LIFT_PADIC_HPP

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
    // Algorithm 3.5.6 (Hensel Quadratic Lift)
    // Cohen, A Course in Computational Algebraic Number Theory (1993)
    using Poly = Polynomial<Integer, Order>;

    // factor correction
    Poly E = (f - g * h) / pk;
    auto dg = divideMod(t * E, {g}, p).R;
    auto dh = divideMod(s * E, {h}, p).R;
    g += pk * dg;
    h += pk * dh;

    // Bézout correction
    Poly C = (Poly(Monomial(f.leadingMonomial().size()), Integer(1)) - s * g - t * h) / pk;
    auto ds = divideMod(t * C, {g}, p).R;
    auto dt = divideMod(s * C, {h}, p).R;
    s += pk * ds;
    t += pk * dt;

    // std::cout << "[DEBUG] end of henselPadicStep()" << std::endl;
    // std::cout << "[DEBUG] deg(g): " << g.degree() << std::endl;
    // std::cout << "[DEBUG] deg(h): " << h.degree() << std::endl;
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
