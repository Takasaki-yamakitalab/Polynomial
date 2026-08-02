#ifndef POLYNOMIAL_FACTOR_HPP
#define POLYNOMIAL_FACTOR_HPP

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
    // std::cout << "[DEBUG] squareFreeDecomposition(" << f << ")" << std::endl;

    // std::cout << "[DEBUG] call multiDerivativeGCD() from squareFreeDecomposition()" << std::endl;
    Poly g = multiDerivativeGCD(f);
    // std::cout << "[DEBUG] return squareFreeDecomposition()" << std::endl;
    // std::cout << "[DEBUG] call divide() from squareFreeDecomposition()" << std::endl;
    Poly w = divide(f, {g}).Q[0];
    // std::cout << "[DEBUG] return squareFreeDecomposition()" << std::endl;

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

    for (int j = 0; j < n; j++) {
        Polynomial<ModInt, Order> g(Monomial(1, j), ModInt(1));
        g = frobenius(g, f);
        for (auto it = g.terms.begin(); it != g.terms.end(); it++) {
            int i = (int)it->first.exp(0);
            Q[i][j] = it->second;
        }
    }

    // Q-I
    for (int i = 0; i < n; i++) {
        Q[i][i] -= 1;
    }
    // std::cout << "[DEBUG] Q-I:\n";
    // for (int i = 0; i < n; ++i) {
    //     std::cout << "  ";
    //     for (int j = 0; j < n; ++j) std::cout << Q[i][j] << ' ';
    //     std::cout << std::endl;
    // }

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
    // std::cout << "[DEBUG] factoring: " << f << std::endl;
    // std::cout << "[DEBUG] basis size: " << basis.size() << std::endl;
    // for (auto& v : basis) {
    //     std::cout << "[DEBUG] basis: ";
    //     for (auto& x : v) std::cout << x << " ";
    //     std::cout << std::endl;
    // }

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
        // std::cout << "[DEBUG] h: " << h << " (at trial: " << trial << ")" << std::endl;
        for (int a = 0; a < ModContext::modulus(); a++) {
            auto t = h;
            t.addTerm(Monomial(1), -ModInt(a));
            if (t.isConstant()) continue;
            auto g = polynomialGCD(f, t);
            // std::cout << "[DEBUG] a: " << a << " gcd: " << g << std::endl;
            if (!g.isConstant() && g != f) {
                auto div = divide(f, {g});
                // std::cout << "[DEBUG] Q: " << div.Q[0] << std::endl;
                // std::cout << "[DEBUG] R: " << div.R << std::endl;
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

    Poly cont(monomialContent(f), content(f));
    if (cont != Poly(Monomial(f.leadingMonomial().size()), Integer(1))) {
        f = pseudoDivide(f, cont).Q[0];
        std::cout << "[DEBUG] content: " << cont << std::endl;
    }
    if (f.leadingCoefficient() < Integer(0)) {
        cont *= Integer(-1);
        f *= Integer(-1);
    }

    // std::cout << "[DEBUG] factorUniVar(" << f << ")" << std::endl;
    // std::cout << "[DEBUG] before findGoodPrime" << std::endl;
    int p = findGoodPrime(f);
    // std::cout << "[DEBUG] p: " << p << std::endl;
    Polynomial<ModInt, LexOrder> fp = modPolynomial(f, p);
    // std::cout << "[DEBUG] fp: " << fp << std::endl;

    // 有限体因数分解 (Berlekamp)
    std::mt19937 rng(std::random_device{}());
    // std::cout << "[DEBUG] call barlekampFactor() from factorUniVar()" << std::endl;
    std::vector<Polynomial<ModInt, LexOrder>> modFactors = berlekampFactor(fp, rng);
    // std::cout << "[DEBUG] return factorUniVar() from barlekampFactor()" << std::endl;

    // Hensel lifting
    // std::cout << "[DEBUG] call buildFactorTree() from factorUniVar()" << std::endl;
    auto root = buildFactorTree<Integer>(modFactors, 0, modFactors.size());
    // std::cout << "[DEBUG] return factorUniVar() from buildFactorTree()" << std::endl;
    const int64_t bound = henselBound(f);
    int64_t pk = p;
    while (pk < bound) {
        // std::cout << "[DEBUG] call henselPadicTree() from factorUniVar()" << std::endl;
        henselPadicTree(f, *root, pk, p);
        // std::cout << "[DEBUG] return factorUniVar() from henselPadicTree()" << std::endl;
        pk *= pk;
    }
    std::vector<Poly> lifted;
    // std::cout << "[DEBUG] call collectLeaves() from factorUniVar()" << std::endl;
    collectLeaves(*root, lifted);
    // std::cout << "[DEBUG] return factorUniVar() from collectLeaves()" << std::endl;

    // DEBUG
    // {
    //     Polynomial<Integer, LexOrder> product(Monomial(1), Integer(1));
    //     for (auto& f : lifted) product *= f;
    //     std::cout << "lift product = " << product << "\n";
    //     std::cout << "original = " << f << "\n";
    // }

    // Recombine (Zassenhaus method)
    // std::cout << "[DEBUG] call centerLiftAll() from factorUniVar()" << std::endl;
    lifted = centerLiftAll(lifted, pk);
    // std::cout << "[DEBUG] return factorUniVar() from centerLiftAll()" << std::endl;
    // std::cout << "[DEBUG] lifted factors:" << std::endl;
    // for (auto& g : lifted) std::cout << "  " << g << std::endl;
    // std::cout << "[DEBUG] call zassenhausRecombine() from factorUniVar()" << std::endl;
    auto ret = zassenhausRecombine(f, lifted);
    // std::cout << "[DEBUG] return factorUniVar() from zassenhausRecombine()" << std::endl;
    ret.push_back(cont);
    return ret;
}

template <class Integer, class Order>
std::vector<Polynomial<Integer, Order>> factorSquareFree(Polynomial<Integer, Order> f) {
    using Poly = Polynomial<Integer, LexOrder>;
    // std::cout << "[DEBUG] factorSquareFree(" << f << ")" << std::endl;

    if (f.isConstant()) return {f};
    if (f.degree() == 1) return {f};

    int64_t B = kroneckerBase(f);
    Poly univarPoly = kroneckerSubstittution(f, B);
    int64_t min_n = univarPoly.degree();
    for (auto it = univarPoly.terms.begin(); it != univarPoly.terms.end(); it++) {
        if (it->first.exp(0) < min_n) min_n = it->first.exp(0);
    }
    Polynomial<Integer, LexOrder> kroneckerShift(Monomial(1, min_n), Integer(1));
    univarPoly = divide(univarPoly, {kroneckerShift}).Q[0];
    // std::cout << "[DEBUG] B: " << B << std::endl;
    // std::cout << "[DEBUG] g: " << univarPoly << std::endl;
    // std::cout << "[DEBUG] Kronecker shift: " << kroneckerShift << std::endl;

    std::vector<Polynomial<Integer, LexOrder>> univarFactors = factorUniVar(univarPoly);
    for (int i = 0; i < min_n; i++) {
        univarFactors.push_back(Poly(Monomial(1, 1), Integer(1)));
    }
    // for (auto& fact : univarFactors) std::cout << "[DEBUG] univarFactors: " << fact << std::endl;

    // univarFactors = inverseKroneckerAll(univarFactors, B, f.leadingMonomial().size());
    std::vector<bool> is_used(univarFactors.size(), false);
    return kroneckerRecombine(f, univarFactors, is_used, B, f.leadingMonomial().size());
}

template <class Integer, class Order>
std::vector<std::pair<Polynomial<Integer, Order>, int64_t>> factorize(Polynomial<Integer, Order> f) {
    using Poly = Polynomial<Integer, Order>;

    std::vector<std::pair<Poly, int64_t>> ans;

    // 定数なら終了
    if (f.isConstant()) return ans;

    // content 除去
    Poly cont(monomialContent(f), content(f));
    if (f.leadingCoefficient() < Integer(0)) {
        cont *= Integer(-1);
        f *= Integer(-1);
    }
    if (cont != Poly(Monomial(f.leadingMonomial().size()), Integer(1))) {
        f = pseudoDivide(f, cont).Q[0];
        std::cout << "[DEBUG] content: " << cont << std::endl;
        ans.push_back({cont, 1});
    }

    // square-free decomposition
    auto sqfs = squareFreeDecomposition(f);

    // 各 square-free 部分を因数分解
    for (auto& sqf : sqfs) {
        auto factors = factorSquareFree(sqf.poly);
        for (auto& h : factors) {
            if (h.isOne()) continue;
            ans.push_back({h, sqf.multiplicity});
            std::cout << "[DEBUG] factor: " << h << std::endl;
        }
    }
    return ans;
}
}  // namespace polynomial

#endif
