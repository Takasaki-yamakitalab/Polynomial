#ifndef FACTOR_HPP
#define FACTOR_HPP

#include <vector>
#include "modint.hpp"
#include "polynomial.hpp"
#include "factor_node.hpp"
#include "hensel_lift_padic.hpp"
#include "kronecker.hpp"
#include "zassenhaus_recombine.hpp"
#include <random>
#include <iostream>

namespace polynomial {

std::vector<std::vector<ModInt>> nullSpace(std::vector<std::vector<ModInt>> A);

template <class Coef, class Order>
std::vector<SquareFreeFactor<Coef, Order>> squareFreeDecomposition(Polynomial<Coef, Order> f);
template <class Order>
std::vector<std::vector<ModInt>> berlekampMatrix(const Polynomial<ModInt, Order>& f);
template <class Order>
std::vector<Polynomial<ModInt, Order>> basisToPolynomial(const std::vector<std::vector<ModInt>>& basis);
template <class Order, class RNG>
Polynomial<ModInt, Order> randomBerlekampPolynomial(const std::vector<Polynomial<ModInt, Order>>& basis, RNG& rng);
template <class Order, class RNG>
std::vector<Polynomial<ModInt, Order>> berlekampFactor(Polynomial<ModInt, Order> f, RNG& rng);
template <class Integer>
std::vector<Polynomial<Integer, LexOrder>> factorUniVar(Polynomial<Integer, LexOrder> f);
template <class Integer, class Order>
std::vector<Polynomial<Integer, Order>> factorSquareFree(Polynomial<Integer, Order> f);
template <class Integer, class Order>
std::vector<std::pair<Polynomial<Integer, Order>, int64_t>> factorize(Polynomial<Integer, Order> f);

template <class Coef, class Order>
std::vector<SquareFreeFactor<Coef, Order>> squareFreeDecomposition(Polynomial<Coef, Order> f) {
    using Poly = Polynomial<Coef, Order>;
    std::vector<SquareFreeFactor<Coef, Order>> ans;

    Poly g = multiDerivativeGCD(f);
    Poly w = divide(f, {g}).Q[0];

    uint64_t i = 1;
    while (!w.isConstant()) {
        Poly y = polynomialGCD(w, g);
        Poly z = divide(w, {y}).Q[0];
        if (!z.isConstant()) ans.push_back({z, i});
        w = y;
        g = divide(g, {y}).Q[0];
        i++;
    }

    if (!g.isConstant()) {
        auto rest = squareFreeDecomposition(g);
        for (auto& s : rest) {
            s.multiplicity++;
            ans.push_back(s);
        }
    }

    return ans;
}
template <class Order>
std::vector<std::vector<ModInt>> berlekampMatrix(const Polynomial<ModInt, Order>& f) {
    int64_t n = f.degree(0);

    std::vector<std::vector<ModInt>> Q(n, std::vector<ModInt>(n));

    Polynomial<ModInt, Order> x;
    x.addTerm(Monomial(std::vector<int64_t>{1}), ModInt(1));

    Polynomial<ModInt, Order> g = x;
    for (int j = 0; j < n; j++) {
        for (auto it = g.terms.begin(); it != g.terms.end(); it++) {
            int i = it->first.exp(0);
            Q[i][j] = it->second;
        }
        g = frobenius(g, f);
    }

    // Q-I
    for (int i = 0; i < n; i++) {
        Q[i][i] -= 1;
    }

    return Q;
}
template <class Order>
std::vector<Polynomial<ModInt, Order>> basisToPolynomial(const std::vector<std::vector<ModInt>>& basis) {
    using Poly = Polynomial<ModInt, Order>;
    std::vector<Poly> ret;
    for (const auto& v : basis) {
        Poly h;
        for (size_t i = 0; i < v.size(); i++) {
            if (v[i].isZero()) continue;
            h.addTerm(Monomial(1, i), v[i]);
        }
        ret.push_back(h);
    }
    return ret;
}
template <class Order, class RNG>
Polynomial<ModInt, Order> randomBerlekampPolynomial(const std::vector<Polynomial<ModInt, Order>>& basis, RNG& rng) {
    using Poly = Polynomial<ModInt, Order>;
    Poly h;
    int p = ModContext::modulus();
    std::uniform_int_distribution<int> dist(0, p - 1);
    for (size_t i = 1; i < basis.size(); i++) {
        int c = dist(rng);
        if (c == 0) continue;
        h += basis[i] * ModInt(c);
    }
    return h;
}
template <class Order, class RNG>
std::vector<Polynomial<ModInt, Order>> berlekampFactor(Polynomial<ModInt, Order> f, RNG& rng) {
    std::vector<Polynomial<ModInt, Order>> ret;
    if (f.degree(0) <= 1) return {f};

    auto basis = nullSpace(berlekampMatrix(f));

    // 定数 1 を追加
    std::vector<ModInt> one(f.degree(0));
    one[0] = 1;
    basis.insert(basis.begin(), one);
    int p = ModContext::modulus();

    if (basis.size() == 1) {
        ret.push_back(f);
        return ret;
    }

    auto basisPoly = basisToPolynomial<Order>(basis);
    for (int trial = 0; trial < 20; trial++) {
        Polynomial<ModInt, Order> h;
        do {
            h = randomBerlekampPolynomial(basisPoly, rng);
        } while (h.isZero());
        for (int a = 0; a < ModContext::modulus(); a++) {
            auto t = h;
            t.addTerm(Monomial(1), -ModInt(a));
            if (t.isConstant()) continue;
            auto g = polynomialGCD(f, t);
            if (!g.isConstant() && g != f) {
                auto div = divide(f, {g});
                assert(div.R.isZero());
                auto left = berlekampFactor(g, rng);
                auto right = berlekampFactor(div.Q[0], rng);
                left.insert(left.end(), right.begin(), right.end());
                return left;
            }
        }
    }
    ret.push_back(f);
    return ret;
}
template <class Integer>
std::vector<Polynomial<Integer, LexOrder>> factorUniVar(Polynomial<Integer, LexOrder> f) {
    using Poly = Polynomial<Integer, LexOrder>;
    if (f.degree() == 1) return {f};

    int p = findGoodPrime(f);
    Polynomial<ModInt, LexOrder> fp = modPolynomial(f, p);

    // 有限体因数分解 (Berlekamp)
    std::mt19937 rng(std::random_device{}());
    std::vector<Polynomial<ModInt, LexOrder>> modFactors = berlekampFactor(fp, rng);

    // Hensel lifting
    auto root = buildFactorTree<Integer, LexOrder>(modFactors, 0, modFactors.size());
    const int bound = henselBound(f);
    int pk = p;
    while (pk < bound) {
        henselPadicTree(f, *root, pk, p);
        pk *= pk;
    }
    std::vector<Poly> lifted;
    collectLeaves(*root, lifted);

    // Recombine (Zassenhaus method)
    lifted = centerLiftAll(lifted, pk);
    return zassenhausRecombine(f, lifted);
}

template <class Integer, class Order>
std::vector<Polynomial<Integer, Order>> factorSquareFree(Polynomial<Integer, Order> f) {
    using Poly = Polynomial<Integer, LexOrder>;

    if (f.isConstant()) return {f};
    if (f.degree() == 1) return {f};

    int64_t B = kroneckerBase(f);
    Polynomial<Integer, LexOrder> univarPoly = kroneckerSubstittution(f, B);

    std::vector<Polynomial<Integer, LexOrder>> univarFactors = factorUniVar(univarPoly);

    univarFactors = inverseKroneckerAll(univarFactors, B, f.leadingMonomial().size());

    return univarFactors;
}
template <class Integer, class Order>
std::vector<std::pair<Polynomial<Integer, Order>, int64_t>> factorize(Polynomial<Integer, Order> f) {
    using Poly = Polynomial<Integer, Order>;

    std::vector<std::pair<Poly, int64_t>> ans;

    // 定数なら終了
    if (f.isConstant()) return ans;

    // content 除去
    Integer c = content(f);
    f /= c;
    if (c != Integer(1)) {
        std::cout << "[DEBUG] coeff: " << c << std::endl;
        Poly p;
        p.addTerm(Monomial(f.leadingMonomial().size(), 0), c);
        ans.push_back({p, 1});
    }

    // square-free decomposition
    auto sqfs = squareFreeDecomposition(f);

    // 各 square-free 部分を因数分解
    for (auto& sqf : sqfs) {
        auto factors = factorSquareFree(sqf.poly);
        for (auto& h : factors) {
            ans.push_back({h, sqf.multiplicity});
            std::cout << "[DEBUG] factor: " << h << std::endl;
        }
    }
    return ans;
}
}  // namespace polynomial

#endif
