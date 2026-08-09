#ifndef POLYNOMIAL_POLYNOMIAL_GCD_HPP
#define POLYNOMIAL_POLYNOMIAL_GCD_HPP

#include <limits>
#include <type_traits>
#include "modint.hpp"
#include "polynomial.hpp"
#include "kronecker.hpp"
#include "coefficient_traits.hpp"
#include "crt.hpp"
#include <set>

namespace polynomial {

template <class Integer, class Order>
Integer findGoodPrime(const Polynomial<Integer, Order>& f, const std::vector<Integer>& exclude = {});

template <class Coef, class Order>
Polynomial<Coef, Order> uniPolynomialGCD(const Polynomial<Coef, Order>& f, const Polynomial<Coef, Order>& g);
template <class Coef, class Order>
Polynomial<Coef, Order> polynomialGCDKronecker(const Polynomial<Coef, Order>& f, const Polynomial<Coef, Order>& g,
                                               const size_t& var_num);
template <class Order>
Polynomial<ModInt, Order> polyModGCDImpl(const Polynomial<ModInt, Order>& f, const Polynomial<ModInt, Order>& g,
                                         const size_t& var_num, const size_t& orig_var_num);
template <class Order>
Polynomial<ModInt, Order> polynomialModGCD(const Polynomial<ModInt, Order>& f, const Polynomial<ModInt, Order>& g,
                                           const size_t& var_num);
template <class Coef, class Order>
Polynomial<Coef, Order> polynomialGCD(const Polynomial<Coef, Order>& f, const Polynomial<Coef, Order>& g,
                                      const size_t& var_num);
template <class Coef, class Order>
Polynomial<Coef, Order> polynomialGCD(const std::vector<Polynomial<Coef, Order>>& F, const size_t& var_num);
template <class Coef, class Order>
Polynomial<Coef, Order> multiDerivativeGCD(const Polynomial<Coef, Order>& f);

template <class Integer, class Order>
Integer findGoodPrime(const Polynomial<Integer, Order>& f, const std::vector<Integer>& exclude) {
    const Integer primes[] = {2,   3,   5,   7,   11,  13,  17,  19,  23,  29,  31,  37,  41,  43,  47,
                              53,  59,  61,  67,  71,  73,  79,  83,  89,  97,  101, 103, 107, 109, 113,
                              127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 199};

    DEBUG_OUT("findGoodPrime(", f, ")");
    for (Integer p : primes) {
        for (auto& ex : exclude)
            if (p == ex) continue;

        auto fp = modPolynomial(f, p);

        // 次数が落ちる素数は除外
        if (fp.degree(0) != f.degree(0)) continue;

        auto df = derivativeMod(f, p);

        DEBUG_OUT("p: ", p);
        DEBUG_OUT("f (mod ", p, "): ", fp);
        DEBUG_OUT("df (mod ", p, "): ", df);
        auto g = uniPolynomialGCD(fp, df);
        DEBUG_OUT("gcd: ", g);
        DEBUG_OUT("TEST: fp/gcd: ", divide(fp, {g}).Q[0]);
        DEBUG_OUT("TEST: df/gcd: ", divide(df, {g}).Q[0]);

        // square-free
        if (g.isConstant()) return p;
    }

    throw std::runtime_error("no suitable prime");
}

template <class Coef, class Order>
static Polynomial<Coef, Order> uniPolyGCDImpl(Polynomial<Coef, Order> f, Polynomial<Coef, Order> g, std::true_type) {
    // 体係数
    using Poly = Polynomial<Coef, Order>;
    while (!g.isZero()) {
        Poly r = polynomialRemainder(f, g);
        f = g;
        g = r;
    }
    if (!f.isZero()) f *= Coef(1) / f.leadingCoefficient();
    return f;
}
template <class Coef, class Order>
static Polynomial<Coef, Order> uniPolyGCDImpl(Polynomial<Coef, Order> f, Polynomial<Coef, Order> g, std::false_type) {
    // 整数係数
    using Poly = Polynomial<Coef, Order>;
    if (g.isZero()) return primitivePart(f);

    Coef c = coefficientGCD(content(f), content(g));
    while (!g.isZero()) {
        Poly r = pseudoRemainder(f, g);
        r = primitivePart(r);
        f = g;
        g = r;
    }
    return c * primitivePart(f);
}

template <class Coef, class Order>
Polynomial<Coef, Order> uniPolynomialGCD(const Polynomial<Coef, Order>& f, const Polynomial<Coef, Order>& g) {
    if (f.isZero()) return g;
    if (g.isZero()) return f;
    return uniPolyGCDImpl(f, g, std::integral_constant<bool, coefficient_traits<Coef>::is_field>());
}

template <class Coef, class Order>
Polynomial<Coef, Order> polynomialGCDKronecker(const Polynomial<Coef, Order>& f, const Polynomial<Coef, Order>& g,
                                               const size_t& var_num) {
    using Poly = Polynomial<Coef, Order>;
    if (f.isZero()) return g;
    if (g.isZero()) return f;
    int64_t B = std::max(kroneckerBase(f), kroneckerBase(g)) * 2;
    Poly f_uni = kroneckerSubstitution(f, B);
    Poly g_uni = kroneckerSubstitution(g, B);
    Poly gcd_uni = uniPolyGCDImpl(f_uni, g_uni, std::integral_constant<bool, coefficient_traits<Coef>::is_field>());
    Poly gcd = inverseKronecker(gcd_uni, B, var_num);
    return gcd;
}

template <class Order>
Polynomial<ModInt, Order> polyModGCDImpl(const Polynomial<ModInt, Order>& f, const Polynomial<ModInt, Order>& g,
                                         const size_t& var_num, const size_t& orig_var_num) {
    using Poly = Polynomial<ModInt, Order>;
    if (var_num == 1) return uniPolynomialGCD(f, g);
    const size_t var = var_num - 1;
    const int64_t P = ModContext::modulus();
    std::vector<ModInt> Alpha;
    std::vector<Poly> H;
    int64_t best_deg = std::numeric_limits<int64_t>::max();
    for (int a = 0; a < P; a++) {
        Poly F = substitute(f, var, ModInt(a));
        Poly G = substitute(g, var, ModInt(a));
        bool bad_eval = false;
        for (size_t i = 0; i < var; i++) {
            if (F.degree(i) < f.degree(i) || G.degree(i) < g.degree(i)) {
                bad_eval = true;
                break;
            }
        }
        if (bad_eval) continue;
        Poly h = polyModGCDImpl(F, G, var_num - 1, orig_var_num);
        if (h.isZero()) return h;
        h = primitivePart(h);
        if (h.degree(0) > best_deg) {
            continue;
        }
        if (h.degree(0) < best_deg) {
            best_deg = h.degree(0);
            Alpha = {a};
            H = {h};
            continue;
        }
        Alpha.push_back(a);
        H.push_back(h);
        if ((int64_t)H.size() < std::min(f.degree(var), g.degree(var)) + 1) continue;
        // H からラグランジュ補間
        Poly gcd;
        Monomial y(orig_var_num);
        y.set_exp(var, 1);
        for (int64_t x_deg = 0; x_deg <= best_deg; x_deg++) {
            std::vector<Poly> coeff;
            for (auto& hh : H) {
                coeff.push_back(hh.coeffPoly(0, x_deg));
            }
            // (x, y) = (Alpha[j], coeff[j]) をラグランジュ補間
            Poly c;
            for (int i = 0; i < Alpha.size(); i++) {
                Poly num(Monomial(orig_var_num), ModInt(1));
                ModInt den(1);
                for (int k = 0; k < Alpha.size(); k++) {
                    if (k == i) continue;
                    num *= (Poly(y, ModInt(1)) + Poly(Monomial(orig_var_num), -Alpha[k]));
                    den *= (Alpha[i] - Alpha[k]);
                }
                c += coeff[i] * num * den.inverse();
            }
            Monomial xn(orig_var_num);
            xn.set_exp(0, x_deg);
            gcd += c * Poly(xn, ModInt(1));
        }
        if (divide(f, {gcd}).R.isZero() && divide(g, {gcd}).R.isZero()) return gcd;
    }
    return Poly();
}

template <class Order>
Polynomial<ModInt, Order> polynomialModGCD(const Polynomial<ModInt, Order>& f, const Polynomial<ModInt, Order>& g,
                                           const size_t& var_num) {
    return polyModGCDImpl(f, g, var_num, var_num);
}

template <class Coef, class Order>
Polynomial<Coef, Order> polynomialGCD(const Polynomial<Coef, Order>& f, const Polynomial<Coef, Order>& g,
                                      const size_t& var_num) {
    if (var_num == 1) return uniPolynomialGCD(f, g);
    if (f.isOne()) return f;
    if (g.isOne()) return g;
    using Poly = Polynomial<Coef, Order>;
    using ModPoly = Polynomial<ModInt, Order>;

    const int64_t default_P = ModContext::modulus();

    // https://mathlandscape.com/prime-table/
    static const int primes[] = {7927, 7933, 7937, 7949, 7951, 7963, 7993, 8009, 8011, 8017,
                                 8039, 8053, 8059, 8069, 8081, 8087, 8089, 8093, 8101, 8111,
                                 8117, 8123, 8147, 8161, 8167, 8171, 8179, 8191, 8209, 8219};

    std::vector<ModPoly> H;
    std::vector<int64_t> P;
    Poly gcd;
    int64_t pp = 1;
    for (int p : primes) {
        ModPoly fp = modPolynomial(f, p);
        ModPoly gp = modPolynomial(g, p);
        bool restart = false;
        for (int i = 0; i < var_num; i++) {
            if (fp.degree(i) < f.degree(i) || gp.degree(i) < g.degree(i)) {
                restart = true;
                break;
            }
        }
        if (restart) continue;
        ModPoly hp = polynomialModGCD(fp, gp, var_num);
        if (hp.isZero()) continue;
        hp = primitivePart(hp);
        H.push_back(hp);
        P.push_back(p);
        pp *= p;
        if (H.size() < 2) continue;
        std::set<Monomial, Order> monos;
        for (auto& h : H) {
            for (auto it = h.terms.begin(); it != h.terms.end(); it++) {
                monos.insert(it->first);
            }
        }
        gcd.clear();
        for (const auto& m : monos) {
            std::vector<int64_t> coeffs;
            for (auto& h : H) {
                coeffs.push_back((int64_t)h.coefficient(m).value());
            }
            auto crt_result = crt(coeffs, P);
            int64_t c = crt_result.first;
            if (c > pp / 2) c -= pp;
            gcd.addTerm(m, Coef(c));
        }
        if (divide(f, {gcd}).R.isZero() && divide(g, {gcd}).R.isZero()) break;
    }
    ModContext::setMod(default_P);
    return gcd;
}

template <class Coef, class Order>
Polynomial<Coef, Order> polynomialGCD(const std::vector<Polynomial<Coef, Order>>& F, const size_t& var_num) {
    using Poly = Polynomial<Coef, Order>;
    assert(!F.empty());
    Poly gcd = F[0];
    for (auto& f : F) {
        // gcd = polynomialGCDKronecker(gcd, f, var_num);
        gcd = polynomialGCD(gcd, f, var_num);
    }
    DEBUG_OUT("Calculate GCD of:");
    for (auto& f : F) DEBUG_OUT("  ", f);
    DEBUG_OUT("result: ", gcd);
    return gcd;
}

template <class Coef, class Order>
Polynomial<Coef, Order> multiDerivativeGCD(const Polynomial<Coef, Order>& f) {
    DEBUG_OUT("multiDerivativeGCD(", f, ")");
    using Poly = Polynomial<Coef, Order>;
    std::vector<Poly> F = {f};
    size_t var_num = (int)f.leadingMonomial().size();
    for (int i = 0; i < var_num; i++) {
        Poly df = differential(f, i);
        if (df.isZero()) continue;
        F.push_back(df);
    }
    Poly gcd = polynomialGCD(F, var_num);
    return gcd;
}

}  // namespace polynomial

#endif
