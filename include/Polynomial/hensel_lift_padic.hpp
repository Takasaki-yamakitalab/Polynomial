#ifndef POLYNOMIAL_HENSEL_LIFT_PADIC_HPP
#define POLYNOMIAL_HENSEL_LIFT_PADIC_HPP

#include "polynomial.hpp"
#include "factor_node.hpp"

namespace polynomial {

template <class Integer, class Order>
void henselPadicStep(const Polynomial<Integer, Order>& f, Polynomial<Integer, Order>& g, Polynomial<Integer, Order>& h,
                     Polynomial<Integer, Order>& s, Polynomial<Integer, Order>& t, const int64_t& pk, const int& p);
template <class Integer, class Order>
void henselPadic(const Polynomial<Integer, Order>& f, Polynomial<Integer, Order>& g, Polynomial<Integer, Order>& h,
                 Polynomial<Integer, Order>& s, Polynomial<Integer, Order>& t, const int& p, const int64_t& bound);
template <class Integer, class Order>
void henselPadicTree(const Polynomial<Integer, Order>& f, FactorNode<Integer, Order>& node, const int64_t& pk,
                     const int& p);

template <class Integer, class Order>
void henselPadicStep(const Polynomial<Integer, Order>& f, Polynomial<Integer, Order>& g, Polynomial<Integer, Order>& h,
                     Polynomial<Integer, Order>& s, Polynomial<Integer, Order>& t, const int64_t& pk, const int& p) {
    // Algorithm 3.5.6 (Hensel Quadratic Lift)
    // Cohen, A Course in Computational Algebraic Number Theory (1993)
    using Poly = Polynomial<Integer, Order>;

    // factor correction
    Poly E = (f - g * h) / Integer(pk);
    auto dg = divideMod(t * E, {g}, p).R;
    auto dh = divideMod(s * E, {h}, p).R;
    g += Integer(pk) * dg;
    h += Integer(pk) * dh;

    // Bézout correction
    Poly C = (Poly(Monomial(f.leadingMonomial().size()), Integer(1)) - s * g - t * h) / Integer(pk);
    auto ds = divideMod(t * C, {g}, p).R;
    auto dt = divideMod(s * C, {h}, p).R;
    s += Integer(pk) * ds;
    t += Integer(pk) * dt;

    DEBUG_OUT("end of henselPadicStep()");
    DEBUG_OUT("deg(g): ", g.degree());
    DEBUG_OUT("deg(h): ", h.degree());
}

template <class Integer, class Order>
void henselPadic(const Polynomial<Integer, Order>& f, Polynomial<Integer, Order>& g, Polynomial<Integer, Order>& h,
                 Polynomial<Integer, Order>& s, Polynomial<Integer, Order>& t, const int& p, const int64_t& bound) {
    int64_t pk = p;
    while (pk < bound) {
        henselPadicStep(f, g, h, s, t, pk);
        pk *= pk;
    }
}

template <class Integer, class Order>
void henselPadicTree(const Polynomial<Integer, Order>& f, FactorNode<Integer, Order>& node, const int64_t& pk,
                     const int& p) {
    if (node.isLeaf()) return;

    auto& g = node.left->factor;
    auto& h = node.right->factor;

    henselPadicStep(f, g, h, node.s, node.t, pk, p);
    DEBUG_OUT("bezout (mod ", pk, "): ", modPolynomial(node.s * g + node.t * h, pk));
    DEBUG_OUT("factor error (mod ", pk, "): ", modPolynomial(g * h, pk) - modPolynomial(f, pk));

    henselPadicTree(g, *node.left, pk, p);
    henselPadicTree(h, *node.right, pk, p);

    node.factor = g * h;
}
}  // namespace polynomial

#endif
