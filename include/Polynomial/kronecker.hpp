#ifndef KRONECKER_HPP
#define KRONECKER_HPP

#include "polynomial.hpp"
#include <algorithm>

namespace polynomial {

template <class Coef, class Order>
int64_t kroneckerBase(const Polynomial<Coef, Order>& f);
template <class Coef, class Order>
Polynomial<Coef, LexOrder> kroneckerSubstittution(const Polynomial<Coef, Order>& f, int64_t B);
template <class Coef>
Polynomial<Coef, LexOrder> inverseKronecker(const Polynomial<Coef, LexOrder>& f, int64_t B, size_t var_num);
template <class Coef>
std::vector<Polynomial<Coef, LexOrder>> inverseKroneckerAll(const std::vector<Polynomial<Coef, LexOrder>>& f, int64_t B, size_t var_num);


template <class Coef, class Order>
int64_t kroneckerBase(const Polynomial<Coef, Order>& f) {
    int64_t B = 1;
    size_t n = f.leadingMonomial().size();
    for (size_t i = 0; i < n; i++) {
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
        Monomial m(std::vector<int64_t>{e});

        ret.addTerm(m, it->second);
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
std::vector<Polynomial<Coef, LexOrder>> inverseKroneckerAll(const std::vector<Polynomial<Coef, LexOrder>>& f, int64_t B, size_t var_num) {
    std::vector<Polynomial<Coef, LexOrder>> ret;
    for (auto& p : f) {
        ret.push_back(inverseKronecker(p, B, var_num));
    }
    return ret;
}
}

#endif
