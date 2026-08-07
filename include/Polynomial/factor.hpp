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
#include "debug_messager.hpp"
#include "polynomial_gcd.hpp"

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
    size_t var_num = f.leadingMonomial().size();
    DEBUG_OUT("squareFreeDecomposition(", f, ")");

    DEBUG_OUT("call multiDerivativeGCD() from squareFreeDecomposition()");
    Poly g = multiDerivativeGCD(f);
    DEBUG_OUT("return squareFreeDecomposition()");
    DEBUG_OUT("call divide() from squareFreeDecomposition()");
    Poly w = divide(f, {g}).Q[0];
    DEBUG_OUT("return squareFreeDecomposition()");

    uint64_t i = 1;
    while (!w.isConstant()) {
        Poly y = polynomialGCD(w, g, var_num);
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

    DEBUG_OUT("Q-I:");
    for (int i = 0; i < n; ++i) {
        DEBUG_OUT("  ", Q[i]);
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

    DEBUG_OUT("factoring: ", f);
    DEBUG_OUT("basis size: ", basis.size());
    for (auto& v : basis) {
        DEBUG_OUT("basis: ", v);
    }

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
        DEBUG_OUT("h: ", h, " (at trial: ", trial, ")");
        for (int a = 0; a < ModContext::modulus(); a++) {
            auto t = h;
            t.addTerm(Monomial(1), -ModInt(a));
            if (t.isConstant()) continue;
            auto g = uniPolynomialGCD(f, t);
            DEBUG_OUT("a: ", a, " gcd: ", g);
            if (!g.isConstant() && g != f) {
                auto div = divide(f, {g});
                DEBUG_OUT("Q: ", div.Q[0]);
                DEBUG_OUT("R: ", div.R);
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
        DEBUG_OUT("content: ", cont);
    }
    if (f.leadingCoefficient() < Integer(0)) {
        cont *= Integer(-1);
        f *= Integer(-1);
    }

    DEBUG_OUT("factorUniVar(", f, ")");
    DEBUG_OUT("before findGoodPrime");
    int p = findGoodPrime(f);
    DEBUG_OUT("p: ", p);
    Polynomial<ModInt, LexOrder> fp = modPolynomial(f, p);
    DEBUG_OUT("fp: ", fp);

    // 有限体因数分解 (Berlekamp)
    std::mt19937 rng(std::random_device{}());
    DEBUG_OUT("call barlekampFactor() from factorUniVar()");
    std::vector<Polynomial<ModInt, LexOrder>> modFactors = berlekampFactor(fp, rng);
    DEBUG_OUT("return factorUniVar() from barlekampFactor()");

    // Hensel lifting
    DEBUG_OUT("call buildFactorTree() from factorUniVar()");
    auto root = buildFactorTree<Integer>(modFactors, 0, modFactors.size());
    DEBUG_OUT("return factorUniVar() from buildFactorTree()");
    DEBUG_OUT("s * left + t * right: ", modPolynomial(root->s * root->left->factor + root->t * root->right->factor, p));
    const int64_t bound = henselBound(f);
    int64_t pk = p;
    int64_t last_pk = pk;
    while (pk < bound) {
        DEBUG_OUT("call henselPadicTree() from factorUniVar()");
        last_pk = pk;
        henselPadicTree(f, *root, pk, p);
        DEBUG_OUT("return factorUniVar() from henselPadicTree()");
        pk *= pk;
    }
    std::vector<Poly> lifted;
    DEBUG_OUT("call collectLeaves() from factorUniVar()");
    collectLeaves(*root, lifted);
    DEBUG_OUT("return factorUniVar() from collectLeaves()");

    // DEBUG
    {
        Polynomial<Integer, LexOrder> product(Monomial(1), Integer(1));
        for (auto& lift : lifted) product *= lift;
        DEBUG_OUT("lift product: ", product);
        DEBUG_OUT("original: ", f);
        DEBUG_OUT("lift product - original (mod ", last_pk, "): ", modPolynomial(product - f, last_pk));
    }

    // Recombine (Zassenhaus method)
    DEBUG_OUT("call centerLiftAll() from factorUniVar()");
    lifted = centerLiftAll(lifted, last_pk);
    DEBUG_OUT("return factorUniVar() from centerLiftAll()");
    DEBUG_OUT("lifted factors:");
    for (auto& g : lifted) DEBUG_OUT("  ", g);
    DEBUG_OUT("lifted factors (mod ", p, "): ");
    for (auto& g : lifted) DEBUG_OUT("  ", modPolynomial(g, p));
    DEBUG_OUT("call zassenhausRecombine() from factorUniVar()");
    auto ret = zassenhausRecombine(f, lifted);
    DEBUG_OUT("return factorUniVar() from zassenhausRecombine()");
    ret.push_back(cont);
    return ret;
}

template <class Integer, class Order>
std::vector<Polynomial<Integer, Order>> factorSquareFree(Polynomial<Integer, Order> f) {
    using Poly = Polynomial<Integer, LexOrder>;
    DEBUG_OUT("factorSquareFree(", f, ")");

    if (f.isConstant()) return {f};
    if (f.degree() == 1) return {f};

    int64_t B = kroneckerBase(f);
    Poly univarPoly = kroneckerSubstitution(f, B);
    int64_t min_n = univarPoly.degree();
    for (auto it = univarPoly.terms.begin(); it != univarPoly.terms.end(); it++) {
        if (it->first.exp(0) < min_n) min_n = it->first.exp(0);
    }
    Polynomial<Integer, LexOrder> kroneckerShift(Monomial(1, min_n), Integer(1));
    univarPoly = divide(univarPoly, {kroneckerShift}).Q[0];
    DEBUG_OUT("B: ", B);
    DEBUG_OUT("g: ", univarPoly);
    DEBUG_OUT("Kronecker shift: ", kroneckerShift);

    std::vector<Polynomial<Integer, LexOrder>> univarFactors = factorUniVar(univarPoly);
    for (int i = 0; i < min_n; i++) {
        univarFactors.push_back(Poly(Monomial(1, 1), Integer(1)));
    }
    for (auto& fact : univarFactors) DEBUG_OUT("univarFactors: ", fact);

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
    // square-free decomposition
    auto sqfs = squareFreeDecomposition(f);
    for (auto& sqf : sqfs) DEBUG_OUT("square-free: ", sqf.poly, " (", sqf.multiplicity, ")");

    // 各 square-free 部分を因数分解
    for (auto& sqf : sqfs) {
        auto factors = factorSquareFree(sqf.poly);
        for (auto& h : factors) {
            if (h.leadingCoefficient() < Integer(0)) {
                cont *= Integer(-1);
                h *= Integer(-1);
            }
            if (h.isOne()) continue;
            ans.push_back({h, sqf.multiplicity});
            DEBUG_OUT("factor: ", h);
        }
    }
    if (cont != Poly(Monomial(f.leadingMonomial().size()), Integer(1))) {
        f = pseudoDivide(f, cont).Q[0];
        DEBUG_OUT("content: ", cont);
        ans.insert(ans.begin(), {cont, 1});
    }
    return ans;
}
}  // namespace polynomial

#endif
