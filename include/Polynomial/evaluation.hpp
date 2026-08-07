#ifndef POLYNOMIAL_EVALUATION_HPP
#define POLYNOMIAL_EVALUATION_HPP

#include "polynomial.hpp"
#include <random>
#include "polynomial_gcd.hpp"

namespace polynomial {

template <class RNG, class Coef, class Order>
Polynomial<Coef, Order> randomEvaluation(const Polynomial<Coef, Order>& f, std::vector<Coef>& values, RNG& rng,
                                         int bound = 10);
template <class Coef, class Order>
bool evaluationIsGood(const Polynomial<Coef, Order>& f, const Polynomial<Coef, Order>& g);

template <class RNG, class Coef, class Order>
Polynomial<Coef, Order> randomEvaluation(const Polynomial<Coef, Order>& f, std::vector<Coef>& values, RNG& rng,
                                         int bound) {
    using Poly = Polynomial<Coef, Order>;

    const int vars = f.leadingMonomial().size();

    values.resize(vars);

    values[0] = Coef(0);

    std::uniform_int_distribution<int> dist(-bound, bound);

    Poly g = f;

    for (int i = 1; i < vars; i++) {
        Coef a;

        do {
            a = Coef(dist(rng));
        } while (a == Coef(0));

        values[i] = a;

        g = substitute(g, i, a);
    }

    return g;
}
template <class Coef, class Order>
bool evaluationIsGood(const Polynomial<Coef, Order>& f, const Polynomial<Coef, Order>& g) {
    if (g.degree() != f.degree()) return false;

    auto d = differential(g, 0);

    if (!uniPolynomialGCD(g, d).isConstant()) return false;

    return true;
}

}  // namespace polynomial

#endif
