#ifndef POLYNOMIAL_KRONECKER_HPP
#define POLYNOMIAL_KRONECKER_HPP

#include "monomial_order.hpp"
#include "polynomial.hpp"
#include <algorithm>
#include "zassenhaus_recombine.hpp"
#include "debug_messager.hpp"

namespace polynomial {

template <class Coef, class Order>
int64_t kroneckerBase(const Polynomial<Coef, Order>& f);
template <class Coef, class Order>
Polynomial<Coef, LexOrder> kroneckerSubstittution(const Polynomial<Coef, Order>& f, int64_t B);
template <class Coef>
Polynomial<Coef, LexOrder> inverseKronecker(const Polynomial<Coef, LexOrder>& f, int64_t B, size_t var_num);
template <class Coef>
std::vector<Polynomial<Coef, LexOrder>> inverseKroneckerAll(const std::vector<Polynomial<Coef, LexOrder>>& f, int64_t B,
                                                            size_t var_num);
template <class Integer, class Order>
std::vector<Polynomial<Integer, Order>> kroneckerRecombine(const Polynomial<Integer, Order>& origin,
                                                           const std::vector<Polynomial<Integer, Order>>& uniFactors,
                                                           std::vector<bool>& is_used, int64_t B, size_t var_num);

template <class Coef, class Order>
int64_t kroneckerBase(const Polynomial<Coef, Order>& f) {
    int64_t B = 1;
    size_t n = f.leadingMonomial().size();
    for (int i = 0; i < n; i++) {
        B = std::max(B, f.degree(i) + 1);
    }
    return B;
}

template <class Coef, class Order>
Polynomial<Coef, LexOrder> kroneckerSubstittution(const Polynomial<Coef, Order>& f, int64_t B) {
    using UniPoly = Polynomial<Coef, LexOrder>;
    UniPoly ret;
    for (auto it = f.terms.begin(); it != f.terms.end(); it++) {
        int64_t e = 0;
        int64_t base = 1;
        for (size_t i = 0; i < it->first.size(); i++) {
            e += it->first.exp(i) * base;
            base *= B;
        }
        ret.addTerm(Monomial(1, e), it->second);
    }
    return ret;
}

template <class Coef>
Polynomial<Coef, LexOrder> inverseKronecker(const Polynomial<Coef, LexOrder>& f, int64_t B, size_t var_num) {
    using Poly = Polynomial<Coef, LexOrder>;
    Poly ret;
    for (auto it = f.terms.begin(); it != f.terms.end(); it++) {
        int64_t e = it->first.exp(0);
        std::vector<int64_t> exps(var_num, 0);
        for (int i = 0; i < var_num; i++) {
            exps[i] = e % B;
            e /= B;
        }
        Monomial m(exps);
        ret.addTerm(m, it->second);
    }
    return ret;
}

template <class Coef>
std::vector<Polynomial<Coef, LexOrder>> inverseKroneckerAll(const std::vector<Polynomial<Coef, LexOrder>>& f, int64_t B,
                                                            size_t var_num) {
    std::vector<Polynomial<Coef, LexOrder>> ret;
    for (auto& p : f) {
        ret.push_back(inverseKronecker(p, B, var_num));
    }
    return ret;
}

template <class Integer, class Order>
std::vector<Polynomial<Integer, Order>> kroneckerRecombine(const Polynomial<Integer, Order>& origin,
                                                           const std::vector<Polynomial<Integer, Order>>& factors,
                                                           std::vector<bool>& is_used, int64_t B, size_t var_num) {
    // if (var_num == 1) return inverseKroneckerAll(factors, B, var_num);
    using Poly = Polynomial<Integer, Order>;
    std::vector<Poly> ret;
    Poly f = origin;
    Integer multiplier(1);
    DEBUG_OUT("origin: ", origin);

    for (int k = 1; k < (int)factors.size(); k++) {
        std::vector<std::vector<int>> subsets;
        enumerateSubset(subsets, is_used, k);
        for (auto& s : subsets) {
            if (s.empty()) continue;

            DEBUG_OUT("s: ", s);

            Poly uf(Monomial(1), Integer(1));
            for (auto& i : s)
                if (!is_used[i]) uf = uf * factors[i];  // <-ここ
            if (uf == Poly(Monomial(var_num), Integer(1))) continue;

            auto candidate = inverseKronecker(uf, B, var_num);
            auto div = pseudoDivide(f, candidate);
            if (!div.R.isZero()) continue;

            for (auto& i : s) is_used[i] = true;
            ret.push_back(candidate);
            f = div.Q[0];
            multiplier *= div.multiplier;
        }
    }
    if (!f.isConstant())
        ret.push_back(f);
    else {
        Integer c = f.leadingCoefficient() / multiplier;
        if (c != Integer(1)) ret.push_back(Poly(Monomial(var_num), c));
    }
    return ret;
}
}  // namespace polynomial

#endif
