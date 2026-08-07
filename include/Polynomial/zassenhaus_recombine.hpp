#ifndef POLYNOMIAL_ZASSENHAUS_RECOMBINE_HPP
#define POLYNOMIAL_ZASSENHAUS_RECOMBINE_HPP

#include "polynomial.hpp"
#include "polynomial_gcd.hpp"
#include "enumerate_subset.hpp"

namespace polynomial {

template <class Integer, class Order>
Polynomial<Integer, Order> centerLift(const Polynomial<Integer, Order>& f, int64_t modulus);

template <class Integer, class Order>
std::vector<Polynomial<Integer, Order>> centerLiftAll(const std::vector<Polynomial<Integer, Order>>& F,
                                                      int64_t modulus);

template <class Coef, class Order>
Polynomial<Coef, Order> productSubset(const std::vector<Polynomial<Coef, Order>>& factors,
                                      const std::vector<int>& subset);

template <class Coef, class Order>
std::vector<Polynomial<Coef, Order>> zassenhausRecombine(Polynomial<Coef, Order> f,
                                                         std::vector<Polynomial<Coef, Order>> lifted);

template <class Integer, class Order>
Polynomial<Integer, Order> centerLift(const Polynomial<Integer, Order>& f, int64_t modulus) {
    Polynomial<Integer, Order> g;
    for (auto it = f.terms.begin(); it != f.terms.end(); it++) {
        Integer c = it->second % Integer(modulus);
        if (c > modulus / 2) c -= Integer(modulus);
        if (c < -modulus / 2) c += Integer(modulus);
        g.addTerm(it->first, c);
    }
    return g;
}
template <class Integer, class Order>
std::vector<Polynomial<Integer, Order>> centerLiftAll(const std::vector<Polynomial<Integer, Order>>& F,
                                                      int64_t modulus) {
    std::vector<Polynomial<Integer, Order>> ret;
    for (auto& f : F) {
        ret.push_back(centerLift(f, modulus));
    }
    return ret;
}

template <class Coef, class Order>
Polynomial<Coef, Order> productSubset(const std::vector<Polynomial<Coef, Order>>& factors,
                                      const std::vector<int>& subset) {
    Polynomial<Coef, Order> p;
    p.addTerm(Monomial(factors[0].leadingMonomial().size()), Coef(1));
    for (size_t i = 0; i < subset.size(); i++) {
        p *= factors[subset[i]];
    }
    return p;
}

template <class Coef, class Order>
std::vector<Polynomial<Coef, Order>> zassenhausRecombine(Polynomial<Coef, Order> f,
                                                         std::vector<Polynomial<Coef, Order>> lifted) {
    using Poly = Polynomial<Coef, Order>;
    std::vector<Poly> result;
    while (!lifted.empty()) {
        bool found = false;
        int m = (int)lifted.size();
        for (int k = 1; k <= m / 2 && !found; k++) {
            std::vector<std::vector<int>> subsets;
            std::vector<bool> exclude(m, false);
            enumerateSubset(subsets, exclude, k);
            for (size_t i = 0; i < subsets.size(); i++) {
                Poly g = productSubset(lifted, subsets[i]);
                Poly d = primitivePart(uniPolynomialGCD(f, g));
                if (d.degree() <= 0 || d == f) continue;
                result.push_back(d);
                f = primitivePart(divide(f, {d}).Q[0]);
                std::vector<Poly> remain;
                for (int j = 0; j < m; j++) {
                    if (std::find(subsets[i].begin(), subsets[i].end(), j) == subsets[i].end())
                        remain.push_back(lifted[j]);
                }
                lifted.swap(remain);
                found = true;
                break;
            }
        }
        if (!found) break;
    }
    result.push_back(f);
    return result;
}
}  // namespace polynomial

#endif
