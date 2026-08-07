#ifndef POLYNOMIAL_HENSEL_LIFT_MULTIVAR_HPP
#define POLYNOMIAL_HENSEL_LIFT_MULTIVAR_HPP

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

template <class Coef, class Order>
HenselFactor<Coef, Order> henselLiftMultiVarStep(const Polynomial<Coef, Order>& f, Polynomial<Coef, Order>& g,
                                         Polynomial<Coef, Order>& h, const Polynomial<Coef, Order>& s,
                                         const Polynomial<Coef, Order>& t, int var, Coef a, int64_t precision) {
    using Poly = Polynomial<Coef, Order>;

    size_t dim = f.leadingMonomial().size();
    Monomial mono(dim);
    Poly factor;
    for (int64_t n = precision; n >= 0; n--) {
        Monomial m = mono;
        m.add_exp(var, n);
        Coef c = combination(precision, n) * pow(-a, precision - n);
        factor.addTerm(m, c);
    }

    Poly E = f - g * h;
    for (int64_t i = 0; i < precision; i++) {
        auto qr = divideByLinearFactor(E, var, a);
        assert(qr.R.isZero());
        E = qr.Q[0];
    }
    g += (t * E) * factor;
    h += (s * E) * factor;
    return {g, h};
}
template <class Coef, class Order>
void henselLiftMultiVar(const Polynomial<Coef, Order>& f, Polynomial<Coef, Order>& g, Polynomial<Coef, Order>& h,
                        int var, const Coef& a) {
    auto ee = extendedEuclid(substitute(g, var, a), substitute(h, var, a));

    int precision = 1;
    while (precision <= f.degree(var)) {
        henselLiftMultiVarStep(f, g, h, ee.s, ee.t, var, a, precision);
        precision *= 2;
    }
}
template <class Coef, class Order>
void henselLiftMultiVarTree(const Polynomial<Coef, Order>& f, FactorNode<Coef, Order>& node, int var, const Coef& a) {
    using FNode = FactorNode<Coef, Order>;
    using Poly = Polynomial<Coef, Order>;

    if (node.isLeaf()) return;

    Poly g = node.left->factor;
    Poly h = node.right->factor;

    // henselLiftMultiVar(f, g, h, var, a);
    int precision = 1;
    while (precision <= f.degree(var)) {
        henselLiftMultiVarStep(f, g, h, node.s, node.t, var, a, precision);
        precision *= 2;
    }

    node.left->factor = g;
    node.right->factor = h;

    henselLiftMultiVarTree(g, *node->left, var, a);
    henselLiftMultiVarTree(h, *node->right, var, a);

    // 念のため再計算
    node.factor = node.left->factor * node.right->factor;
}
}  // namespace polynomial
#endif
