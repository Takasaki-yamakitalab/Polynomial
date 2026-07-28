#include "Polynomial/hensel_lift_padic.hpp"
#include "Polynomial/polynomial.hpp"
#include "Polynomial/factor_node.hpp"

namespace polynomial {
template <class Coef, class Order>
void henselPadicStep(const Polynomial<Coef, Order>& f, Polynomial<Coef, Order>& g, Polynomial<Coef, Order>& h,
                     Polynomial<Coef, Order>& s, Polynomial<Coef, Order>& t, const Coef& pk) {
    using Poly = Polynomial<Coef, Order>;

    // factor correction
    Poly E = (f - g * h) / pk;
    Poly u = s * E;
    Poly v = t * E;
    DivisionResult<Coef, Order> divU = divide(u, {h});
    DivisionResult<Coef, Order> divV = divide(v, {g});
    Poly q = divU.Q[0];
    Poly r = divU.R;
    Poly a = divV.R;
    Poly b = r + q * g;
    g += pk * a;
    h += pk * b;

    // Bézout correction
    Poly one(Monomial(f.leadingMonomial().size()), Coef(1));
    Poly C = (one - s * g - t * h) / pk;
    Poly uu = s * C;
    Poly vv = t * C;
    DivisionResult<Coef, Order> divUU = divide(uu, {h});
    DivisionResult<Coef, Order> divVV = divide(vv, {g});
    Poly qq = divUU.Q[0];
    Poly rr = divUU.R;
    Poly ds = divVV.R;
    Poly dt = rr + qq * g;
    s += pk * ds;
    t += pk * dt;
}

template <class Coef, class Order>
void henselPadic(const Polynomial<Coef, Order>& f, Polynomial<Coef, Order>& g, Polynomial<Coef, Order>& h,
                 Polynomial<Coef, Order>& s, Polynomial<Coef, Order>& t, const Coef& p, const Coef& bound) {
    Coef pk = p;
    while (pk < bound) {
        henselPadicStep(f, g, h, s, t, pk);
        pk *= pk;
    }
}

template <class Coef, class Order>
void henselPadicTree(const Polynomial<Coef, Order>& f, FactorNode<Coef, Order>& node, const Coef& pk) {
    if (node.isLeaf()) return;

    auto& g = node.left->factor;
    auto& h = node.right->factor;

    henselPadicStep(f, g, h, node.s, node.t, pk);

    henselPadicTree(g, *node.left, pk);
    henselPadicTree(h, *node.right, pk);

    node.factor = g * h;
}

}  // namespace polynomial
