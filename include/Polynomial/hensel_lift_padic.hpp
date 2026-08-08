#ifndef POLYNOMIAL_HENSEL_LIFT_PADIC_HPP
#define POLYNOMIAL_HENSEL_LIFT_PADIC_HPP

#include "polynomial.hpp"
#include "factor_node.hpp"

namespace polynomial {

template <class Integer, class Order>
void henselPadicStep(const Polynomial<Integer, Order>& f, Polynomial<Integer, Order>& g, Polynomial<Integer, Order>& h,
                     Polynomial<Integer, Order>& s, Polynomial<Integer, Order>& t, const int64_t& pk, const int& p);
template <class Integer, class Order>
void henselPadicTree(const Polynomial<Integer, Order>& target, FactorNode<Integer, Order>& node, const int64_t& pk,
                     const int& p);

template <class Integer, class Order>
void henselPadicStep(const Polynomial<Integer, Order>& f, Polynomial<Integer, Order>& g, Polynomial<Integer, Order>& h,
                     Polynomial<Integer, Order>& s, Polynomial<Integer, Order>& t, const int64_t& pk, const int& p) {
    // Algorithm 3.5.6 (Hensel Quadratic Lift)
    // Cohen, A Course in Computational Algebraic Number Theory (1993)
    using Poly = Polynomial<Integer, Order>;

    // factor correction
    auto fgh = f - g * h;
    Poly E = (f - g * h) / Integer(pk);
    auto div = divideMod(t * E, {g}, pk);
    Poly q = div.Q[0];
    Poly dg = div.R;
    auto dh = mod(s * E + h * q, pk);
    // DEBUG
    {
        auto factorCorrection = mod(dg * h + g * dh - E, p);
        DEBUG_OUT("factor correction residual (mod ", p, "): ", factorCorrection);
        DEBUG_OUT("f - gh: ", fgh);
        DEBUG_OUT("E: ", E);
        DEBUG_OUT("f - gh - pk*E: ", fgh - Integer(pk) * E);
        DEBUG_OUT("(f - gh - pk*E) mod p: ", mod(fgh - Integer(pk) * E, p));
    }
    g += Integer(pk) * dg;
    h += Integer(pk) * dh;
    g = mod(g, pk*pk);
    h = mod(h, pk*pk);
    // DEBUG
    {
        DEBUG_OUT("pk: ", pk);
        DEBUG_OUT("p: ", p);
        DEBUG_OUT("E: ", E);
        DEBUG_OUT("t*E mod p: ", mod(t * E, p));
        DEBUG_OUT("dg: ", dg);
        DEBUG_OUT("q: ", q);
        DEBUG_OUT("dh: ", dh);
        DEBUG_OUT("linear correction residual: ", mod(dg * h + g * dh - E, pk));
    }

    // Bézout correction
    Poly C = mod((Poly(Monomial(1), Integer(1)) - s * g - t * h) / Integer(pk), pk);
    div = divideMod(t * C, {g}, pk);
    q = div.Q[0];
    auto dt = div.R;
    auto ds = mod(s * C + h * q, pk);
    // DEBUG
    {
        auto residual = mod(ds * g + dt * h - C, pk);
        DEBUG_OUT("bezout correction residual (mod ", pk, "): ", residual);
    }
    s += Integer(pk) * ds;
    t += Integer(pk) * dt;
    s = mod(s, pk*pk);
    t = mod(t, pk*pk);

    DEBUG_OUT("end of henselPadicStep()");
    DEBUG_OUT("deg(g): ", g.degree());
    DEBUG_OUT("deg(h): ", h.degree());
}

template <class Integer, class Order>
void henselPadicTree(const Polynomial<Integer, Order>& target, FactorNode<Integer, Order>& node, const int64_t& pk,
                     const int& p) {
    if (node.isLeaf()) return;

    auto& g = node.left->factor;
    auto& h = node.right->factor;

    DEBUG_OUT("pre-lift bezout (mod ", pk, "): ", modPolynomial(node.s * g + node.t * h, pk));
    DEBUG_OUT("pre-lift factor error (mod ", pk, "): ", modPolynomial(g * h - target, pk));
    henselPadicStep(target, g, h, node.s, node.t, pk, p);
    DEBUG_OUT("post-lift bezout (mod ", pk * pk, "): ", modPolynomial(node.s * g + node.t * h, pk * pk));
    DEBUG_OUT("post-lift factor error (mod ", pk * pk, "): ", modPolynomial(g * h - target, pk * pk));
    node.factor = g * h;

    henselPadicTree(g, *node.left, pk, p);
    henselPadicTree(h, *node.right, pk, p);
}
}  // namespace polynomial

#endif
