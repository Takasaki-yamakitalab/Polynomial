#include <cstdint>
#include <cstdio>
#include <random>
#include <type_traits>
#include "Polynomial/coefficient.hpp"
#include "Polynomial/coefficient_traits.hpp"
#include "Polynomial/polynomial.hpp"
#include "Polynomial/modint.hpp"

namespace polynomial {

// Polynomial member functions
template <class Coef, class Order>
bool Polynomial<Coef, Order>::isZero() const {
    return terms.empty();
}
template <class Coef, class Order>
bool Polynomial<Coef, Order>::isConstant() const {
    if (terms.empty()) return true;
    if (terms.size() != 1) return false;
    return leadingMonomial().degree() == 0;
}
template <class Coef, class Order>
bool Polynomial<Coef, Order>::isMonomial() const {
    return terms.size() == 1;
}
template <class Coef, class Order>
int64_t Polynomial<Coef, Order>::degree() const {
    if (isZero()) return -1;
    int64_t deg = 0;
    for (auto it = terms.begin(); it != terms.end(); it++) {
        int64_t n = terms.mono.degree();
        if (n > deg) deg = n;
    }
    return deg;
}
template <class Coef, class Order>
int64_t Polynomial<Coef, Order>::degree(const int& var) const {
    if (isZero()) return -1;
    int64_t deg = 0;
    for (auto it = terms.begin(); it != terms.end(); it++) {
        int64_t n = terms.mono.exp(var);
        if (n > deg) deg = n;
    }
    return deg;
}
template <class Coef, class Order>
void Polynomial<Coef, Order>::normalize() {
    for (auto it = terms.begin(); it != terms.end();) {
        if (it->second == Coef(0)) {
            it = terms.erase(it);
        } else {
            it++;
        }
    }
}
template <class Coef, class Order>
int64_t Polynomial<Coef, Order>::leadingDegree() const {
    return leadingMonomial().degree();
}
template <class Coef, class Order>
Polynomial<Coef, Order> Polynomial<Coef, Order>::leadingPolynomial() const {
    return Polynomial(leadingTerm());
}
template <class Coef, class Order>
Coef Polynomial<Coef, Order>::coefficient(const Monomial& m) const {
    auto it = terms.find(m);
    if (it == terms.end()) return Coef(0);
    return it->second;
}
template <class Coef, class Order>
void Polynomial<Coef, Order>::clear() {
    terms.clear();
}
template <class Coef, class Order>
void Polynomial<Coef, Order>::addTerm(const Monomial& m, const Coef& c) {
    if (c == Coef(0)) return;
    auto it = terms.find(m);
    if (it == terms.end()) {
        terms.emplace(m, c);
        return;
    }
    it->second += c;
    if (it->second == Coef(0)) terms.erase(it);
}
template <class Coef, class Order>
void Polynomial<Coef, Order>::addTerm(const Term<Coef>& t) {
    addTerm(t.mono, t.coef);
}
template <class Coef, class Order>
bool Polynomial<Coef, Order>::reducibleBy(const Polynomial& g) const {
    if (isZero()) return false;
    if (g.isZero()) return false;
    return g.leadingMonomial().divides(leadingMonomial());
}
template <class Coef, class Order>
Polynomial<Coef, Order>& Polynomial<Coef, Order>::reduceOnceBy(const Polynomial<Coef, Order>& g) {
    assert(reducibleBy(g));
    auto t = leadingTerm() / g.leadingTerm();
    *this -= g * t;
    return *this;
}
template <class Coef, class Order>
Polynomial<Coef, Order>& Polynomial<Coef, Order>::reduceBy(const Polynomial<Coef, Order>& g) {
    while (reducibleBy(g)) {
        reduceOnceBy(g);
    }
    return *this;
}
template <class Coef, class Order>
Polynomial<Coef, Order>& Polynomial<Coef, Order>::normalForm(const std::vector<Polynomial<Coef, Order>>& G) {
    auto tmp = divide(*this, G).R;
    terms.swap(tmp.terms);
    return *this;
}
template <class Coef, class Order>
bool Polynomial<Coef, Order>::operator==(const Polynomial<Coef, Order>& rhs) const {
    return terms == rhs.terms;
}
template <class Coef, class Order>
Polynomial<Coef, Order>& Polynomial<Coef, Order>::operator+=(const Polynomial<Coef, Order>& rhs) {
    for (auto& t : rhs.terms) {
        addTerm(t.first, t.second);
    }
    return *this;
}
template <class Coef, class Order>
Polynomial<Coef, Order>& Polynomial<Coef, Order>::operator+=(const Term<Coef>& rhs) {
    addTerm(rhs);
    return *this;
}
template <class Coef, class Order>
Polynomial<Coef, Order>& Polynomial<Coef, Order>::operator-=(const Polynomial<Coef, Order>& rhs) {
    for (auto& t : rhs.terms) {
        addTerm(t.first, -t.second);
    }
    return *this;
}
template <class Coef, class Order>
Polynomial<Coef, Order>& Polynomial<Coef, Order>::operator*=(const Coef& c) {
    if (c == Coef(0)) {
        clear();
        return *this;
    }
    for (auto& term : terms) {
        term.second *= c;
    }
    normalize();
    return *this;
}
template <class Coef, class Order>
Polynomial<Coef, Order>& Polynomial<Coef, Order>::operator*=(const Term<Coef>& t) {
    Polynomial<Coef, Order> tmp;
    for (auto& term : terms) {
        tmp.addTerm(term.first * t.mono, term.second * t.coef);
    }
    terms.swap(tmp.terms);
    return *this;
}
template <class Coef, class Order>
Polynomial<Coef, Order>& Polynomial<Coef, Order>::operator*=(const Polynomial<Coef, Order>& p) {
    auto temp = *this * p;
    terms.swap(temp.terms);
    return *this;
}
template <class Coef, class Order>
Polynomial<Coef, Order>& Polynomial<Coef, Order>::operator/=(const Coef& c) {
    assert(c != Coef(0));
    for (auto& term : terms) {
        term.second /= c;
    }
    normalize();
    return *this;
}
template <class Coef, class Order>
Polynomial<Coef, Order> Polynomial<Coef, Order>::operator-() const {
    Polynomial<Coef, Order> ret = *this;
    for (auto& t : ret.terms) {
        t.second = -t.second;
    }
    return ret;
}

template <class Coef, class Order>
Polynomial<Coef, Order> operator+(const Polynomial<Coef, Order>& a, const Polynomial<Coef, Order>& b) {
    Polynomial<Coef, Order> ret = a;
    for (auto itr = b.terms.begin(); itr != b.terms.end(); itr++) {
        ret.addTerm(itr->first, itr->second);
    }
    return ret;
}
template <class Coef, class Order>
Polynomial<Coef, Order> operator-(const Polynomial<Coef, Order>& a, const Polynomial<Coef, Order>& b) {
    Polynomial<Coef, Order> ret = a;
    for (auto itr = b.terms.begin(); itr != b.terms.end(); itr++) {
        ret.addTerm(itr->first, -itr->second);
    }
    return ret;
}

template <class Coef, class Order>
Polynomial<Coef, Order> operator*(const Polynomial<Coef, Order>& a, const Polynomial<Coef, Order>& b) {
    Polynomial<Coef, Order> ret;
    for (auto itr1 = a.terms.begin(); itr1 != a.terms.end(); itr1++) {
        for (auto itr2 = b.terms.begin(); itr2 != b.terms.end(); itr2++) {
            ret.addTerm(itr1->first * itr2->first, itr1->second * itr2->second);
        }
    }
    return ret;
}
template <class Coef, class Order>
Polynomial<Coef, Order> operator*(Polynomial<Coef, Order> p, const Term<Coef>& t) {
    p *= t;
    return p;
}
template <class Coef, class Order>
Polynomial<Coef, Order> operator*(Polynomial<Coef, Order> p, const Coef& c) {
    Polynomial<Coef, Order> ret = p;
    ret *= c;
    return ret;
}
template <class Coef, class Order>
Polynomial<Coef, Order> operator*(const Coef& c, Polynomial<Coef, Order> p) {
    Polynomial<Coef, Order> ret = p;
    ret *= c;
    return ret;
}

template <class Coef, class Order>
Polynomial<Coef, Order> differential(const Polynomial<Coef, Order>& p, int var) {
    Polynomial<Coef, Order> ret;
    for (auto& term : p.terms) {
        if (term.first.exp[var] == 0) continue;
        Monomial mono = term.first;
        Coef coeff = term.second * Coef(mono.exp(var));
        mono.add_exp(var, -1);
        ret.addTerm(mono, coeff);
    }
    return ret;
}

template <class Coef, class Order>
DivisionResult<Coef, Order> divide(Polynomial<Coef, Order> f, const std::vector<Polynomial<Coef, Order>>& G) {
    std::vector<Polynomial<Coef, Order>> Q(G.size());
    Polynomial<Coef, Order> R;

    while (!f.isZero()) {
        bool reduced = false;
        for (size_t i = 0; i < G.size(); i++) {
            if (G[i].isZero()) continue;
            if (!f.reducibleBy(G[i])) continue;

            Q[i] += f.leadingTerm() / G[i].leadingTerm();
            f.reduceOnceBy(G[i]);
            reduced = true;
            break;
        }
        if (!reduced) {
            R += f.leadingTerm();
            f.eraseLeadingTerm();
        }
    }
    return {Q, R};
}

template <class Coef, class Order>
static Polynomial<Coef, Order> primitivePartImpl(Polynomial<Coef, Order> p, std::true_type) {
    // 体系数
    if (p.isZero()) return p;
    p /= p.leadingCoefficient();
    return p;
}
template <class Coef, class Order>
static Polynomial<Coef, Order> primitivePartImpl(Polynomial<Coef, Order> p, std::false_type) {
    // 整数係数
    if (p.isZero()) return p;
    Coef c = content(p);
    if (c != Coef(1)) p /= c;
    return p;
}
template <class Coef, class Order>
Polynomial<Coef, Order> primitivePart(Polynomial<Coef, Order> p) {
    return primitivePartImpl(p, std::integral_constant<bool, coefficient_traits<Coef>::is_field>());
}

template <class Coef, class Order>
DivisionResult<Coef, Order> pseudoDivide(Polynomial<Coef, Order> f, const Polynomial<Coef, Order>& g) {
    using Poly = Polynomial<Coef, Order>;
    DivisionResult<Coef, Order> ret;
    assert(!g.isZero());

    Poly R = f;
    Poly Q;

    int64_t m = f.degree();
    int64_t n = g.degree();

    if (m < n) {
        ret.R = f;
        ret.Q.push_back(Q);
        return ret;
    }

    Coef c = g.leadingCoefficient();
    while (!R.isZero() && R.degree() >= g.degree()) {
        Monomial u = R.leadingMonomial() / g.leadingMonomial();
        Coef a = R.leadingCoefficient();
        Term<Coef> t{u, a};
        Q *= c;
        Q += t;
        R *= c;
        R -= g * t;
    }
    ret.Q.push_back(Q);
    ret.R = R;
    return ret;
}
template <class Coef, class Order>
Polynomial<Coef, Order> pseudoRemainder(const Polynomial<Coef, Order>& f, const Polynomial<Coef, Order>& g) {
    return pseudoDivide(f, g).R;
}
template <class Coef, class Order>
Polynomial<Coef, Order> polynomialRemainder(Polynomial<Coef, Order> f, const Polynomial<Coef, Order>& g) {
    assert(!g.isZero());
    while (!f.isZero()) {
        if (!f.reducibleBy(g)) break;
        f.reduceOnceBy(g);
    }
    return f;
}

template <class Coef, class Order>
static Polynomial<Coef, Order> polynomialGCDImpl(Polynomial<Coef, Order> f, Polynomial<Coef, Order> g, std::true_type) {
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
static Polynomial<Coef, Order> polynomialGCDImpl(Polynomial<Coef, Order> f, Polynomial<Coef, Order> g, std::false_type) {
    // 整数係数
    using Poly = Polynomial<Coef, Order>;
    while (!g.isZero()) {
        Poly r = primitivePart(pseudoRemainder(f, g));
        f = g;
        g = r;
    }
    return primitivePart(f);
}
template <class Coef, class Order>
Polynomial<Coef, Order> polynomialGCD(Polynomial<Coef, Order> f, Polynomial<Coef, Order> g) {
    if (f.isZero()) return g;
    if (g.isZero()) return f;
    return polynomialGCDImpl(f, g, std::integral_constant<bool, coefficient_traits<Coef>::is_field>());
}

template <class Coef, class Order>
Polynomial<Coef, Order> pow(Polynomial<Coef, Order> p, int64_t e) {
    using Poly = Polynomial<Coef, Order>;
    assert(e >= 0);

    Poly r;
    r.addTerm(Monomial(std::vector<int64_t>(p.leadingMonomial().size(), 0)), Coef(1));
    while (e) {
        if (e & 1) r = r * p;
        p = p * p;
        e >>= 1;
    }
    return r;
}

template <class Coef, class Order>
Polynomial<Coef, Order> substitute(const Polynomial<Coef, Order>& f, int var, const Polynomial<Coef, Order>& g) {
    using Poly = Polynomial<Coef, Order>;

    Poly ans;

    for (auto it = f.terms.begin(); it != f.terms.end(); ++it) {
        const Monomial& m = it->first;
        Coef c = it->second;

        Poly term;
        term.addTerm(Monomial(std::vector<int>(m.size(), 0)), c);

        for (size_t i = 0; i < m.size(); i++) {
            if (m.exp(i) == 0) continue;

            if ((int)i == var) {
                term = term * pow(g, m.exp(i));
            } else {
                std::vector<int64_t> e(m.size(), 0);
                e[i] = m.exp(i);

                Poly t;
                t.addTerm(Monomial(e), Coef(1));

                term = term * t;
            }
        }
        ans += term;
    }
    return ans;
}
template <class Coef, class Order>
Polynomial<Coef, Order> substitute(const Polynomial<Coef, Order>& f, int var, const Coef& value) {
    using Poly = Polynomial<Coef, Order>;

    Poly ans;

    for (auto it = f.terms.begin(); it != f.terms.end(); ++it) {
        Monomial mono = it->first;
        Coef coef = it->second;

        int64_t e = mono.exp(var);

        Coef p = Coef(1);
        for (int i = 0; i < e; i++) p *= value;

        coef *= p;

        mono.add_exp(var, -e);
        ans.addTerm(mono, coef);
    }
    return ans;
}
template <class Coef, class Order>
Polynomial<Coef, Order> substitute(Polynomial<Coef, Order> f, const std::vector<Polynomial<Coef, Order>>& g) {
    assert(g.size() == f.leadingMonomial().size());

    for (int i = (int)g.size() - 1; i >= 0; i--) f = substitute(f, i, g[i]);

    return f;
}

template <class Coef, class Order>
Polynomial<Coef, Order> sPolynomial(const Polynomial<Coef, Order>& f, const Polynomial<Coef, Order>& g) {
    assert(!f.isZero() && !g.isZero());
    Monomial l = lcm(f.leadingMonomial(), g.leadingMonomial());
    Term<Coef> common{l, Coef(1)};
    return f * (common / f.leadingTerm()) - g * (common / g.leadingTerm());
}

template <class Coef, class Order>
int findReducer(const std::vector<Polynomial<Coef, Order>>& G, const Monomial& m) {
    for (int i = 0; i < (int)G.size(); i++) {
        if (G[i].leadingMonomial().divides(m)) return i;
    }
    return -1;
}

template <class Coef, class Order>
Polynomial<Coef, Order> multiDerivativeGCD(Polynomial<Coef, Order> f) {
    using Poly = Polynomial<Coef, Order>;
    Poly g = f;
    int n = f.leadingMonomial().size();
    for (int i = 0; i < n; i++) {
        Poly df = differential(f, i);

        if (df.isZero()) continue;

        g = polynomialGCD(g, df);

        if (g.isConstant()) break;
    }
    return g;
}

template <class Coef, class Order>
DivisionResult<Coef, Order> divideByLinearFactor(const Polynomial<Coef, Order>& p, int var, Coef a) {
    using T = Term<Coef>;
    std::vector<T> termsVec;
    int64_t max_n = 0;
    for (auto it = p.terms.begin(); it != p.terms.end(); it++) {
        termsVec.push_back({it->first, it->second});
        if (it->first.exp(var) > max_n) max_n = it->first.exp(var);
    }

    std::vector<std::vector<T>> termsMat(max_n + 1);
    for (auto& t : termsVec) {
        int64_t n = t.mono.exp(var);
        termsMat[max_n - n].push_back(t);
    }

    for (int i = 0; i < (int)termsMat.size(); i++) {
        if (i != (int)termsMat.size() - 1) {
            for (int j = 0; j < (int)termsMat[i].size(); j++) {
                termsMat[i][j].mono.add_exp(var, -1);
            }
        }
        if (i == 0) continue;
        for (int j = 0; j < (int)termsMat[i - 1].size(); j++) {
            T t = termsMat[i - 1][j];
            if (i != (int)termsMat.size() - 1) t.mono.add_exp(var, -1);
            termsMat[i].push_back(a * t);
        }
    }

    Polynomial<Coef, Order> Q;
    Polynomial<Coef, Order> R;
    for (int i = 0; i < (int)termsMat.size() - 1; i++) {
        for (int j = 0; j < (int)termsMat[i].size(); j++) {
            Q.addTerm(termsMat[i][j]);
        }
    }
    for (int j = 0; j < (int)termsMat[termsMat.size() - 1].size(); j++) {
        R.addTerm(termsMat[termsMat.size() - 1][j]);
    }
    return DivisionResult<Coef, Order>{{Q}, R};
}
template <class Coef, class Order>
ExtendEuclidResult<Coef, Order> extendedEuclid(Polynomial<Coef, Order> f, Polynomial<Coef, Order> g) {
    using Poly = Polynomial<Coef, Order>;

    Poly s0, s1;
    Poly t0, t1;
    // s0 = 1
    s0.addTerm(Monomial(std::vector<int>(f.leadingMonomial().size(), 0)), Coef(1));
    // s1 = 0
    // t0 = 0
    // t1 = 1
    t1.addTerm(Monomial(std::vector<int>(f.leadingMonomial().size(), 0)), Coef(1));

    while (!g.isZero()) {
        DivisionResult<Coef, Order> div = divide(f, {g});
        Poly q = div.Q[0];
        Poly r = div.R;
        Poly s = s0 - q * s1;
        Poly t = t0 - q * t1;
        f = g;
        g = r;
        s0 = s1;
        s1 = s;
        t0 = t1;
        t1 = t;
    }

    Coef c = f.leadingCoefficient();
    if (c != Coef(1)) {
        f *= Coef(1) / c;
        s0 *= Coef(1) / c;
        t0 *= Coef(1) / c;
    }

    return {f, s0, t0};
}
template <class Coef, class Order>
Polynomial<Coef, Order> shiftVariable(const Polynomial<Coef, Order>& f, int var, const Coef& a) {
    using Poly = Polynomial<Coef, Order>;
    Poly ret;
    for (auto it = f.terms.begin(); it != f.terms.end(); it++) {
        Monomial mono = it->first;
        Coef c = it->second;
        int64_t e = mono.exp(var);
        mono.add_exp(var, -e);
        Poly base;
        base.addTerm(mono, c);
        Poly expand;
        for (int k = 0; k <= e; k++) {
            Monomial m = mono;
            m.add_exp(var, k);
            Coef coeff = c * Coef(combination(e, k)) * pow(a, e - k);
            expand.addTerm(m, coeff);
        }
        ret += expand;
    }
    return ret;
}
template <class Coef, class Order>
Polynomial<Coef, Order> truncateByPower(const Polynomial<Coef, Order>& f, int var, int k) {
    Polynomial<Coef, Order> r;
    for (auto it = f.terms.begin(); it != f.terms.end(); it++) {
        if (it->first.exp(var) >= k) continue;
        r.addTerm(it->first, it->second);
    }
    return r;
}

template <class Coef, class Order>
Polynomial<Coef, Order> product(const std::vector<Polynomial<Coef, Order>>& factors) {
    using Poly = Polynomial<Coef, Order>;
    if (factors.empty()) return Poly();
    Poly result = factors[0];
    for (size_t i = 1; i < factors.size(); i++) {
        result *= factors;
    }
    return result;
}
template <class Coef, class Order>
Polynomial<Coef, Order> product(const std::vector<Polynomial<Coef, Order>>& factors, size_t begin, size_t end) {
    using Poly = Polynomial<Coef, Order>;
    if (factors.empty()) return Poly();
    Poly result = factors[begin];
    for (size_t i = begin + 1; i < end; i++) {
        result *= factors;
    }
    return result;
}

template <class Coef, class Order>
int64_t landauMignotteBound(const Polynomial<Coef, Order>& f) {
    int64_t n = f.degree();
    int64_t infNorm = 0;
    for (auto it = f.terms.begin(); it != f.terms.end(); it++) {
        infNorm = std::max(infNorm, std::abs(it->second));
    }
    double B = std::pow(2.0, n) * std::sqrt(double(n + 1)) * double(infNorm);
    return (int64_t)std::ceil(B);
}
template <class Coef, class Order>
int64_t henselBound(const Polynomial<Coef, Order>& f) {
    return 2 * landauMignotteBound(f);
}


template <class Coef, class Order>
bool isFactor(const Polynomial<Coef, Order>& f, const Polynomial<Coef, Order>& g) {
    auto div = divide(f, {g});
    if (!div.R.isZero()) return false;
    return div.Q[0] * g == f;
}

template <class Coef, class Order>
Coef contentImpl(const Polynomial<Coef, Order>& f, std::true_type) {
    return Coef(1);
}
template <class Coef, class Order>
Coef contentImpl(const Polynomial<Coef, Order>& f, std::false_type) {
    if (f.isZero()) return Coef(0);
    auto it = f.terms.begin();
    Coef g = abs(it->second);
    it++;
    for (; it != f.terms.end(); it++) {
        g = coefficientGCD(g, abs(it->second));
    }
    return g;
}
template <class Coef, class Order>
Coef content(const Polynomial<Coef, Order>& f) {
    return contentImpl(f, std::integral_constant<bool, coefficient_traits<Coef>::is_field>());
}

template <class Coef, class Order>
Polynomial<ModInt, Order> modPolynomial(const Polynomial<Coef, Order>& f, int p) {
    ModInt::setMod(p);
    Polynomial<ModInt, Order> ret;
    for (auto it = f.terms.begin(); it != f.terms.end(); it++) {
        ModInt v((int)(it->second + 0.5));
        if (v.value() != 0) ret.addTerm(it->first, v);
    }
    return ret;
}
template <class Order>
Polynomial<ModInt, Order> derivativeMod(const Polynomial<int, Order>& f, int p) {
    auto df = differential(f, 0);
    return modPolynomial(df, p);
}
template <class Order>
int findGoodPrime(const Polynomial<int, Order>& f) {
    static const int primes[] = {2,   3,   5,   7,   11,  13,  17,  19,  23,  29,  31,  37,  41,  43,  47,
                                 53,  59,  61,  67,  71,  73,  79,  83,  89,  97,  101, 103, 107, 109, 113,
                                 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 199};

    for (int p : primes) {
        auto fp = modPolynomial(f, p);

        // 次数が落ちる素数は除外
        if (fp.degree(0) != f.degree(0)) continue;

        auto df = derivativeMod(f, p);

        auto g = polynomialGCD(fp, df);

        // square-free
        if (g.isConstant()) return p;
    }

    throw std::runtime_error("no suitable prime");
}
template <class Order>
Polynomial<ModInt, Order> powModPolynomial(Polynomial<ModInt, Order> a, int64_t n,
                                           const Polynomial<ModInt, Order>& mod) {
    Polynomial<int, Order> ret;
    ret.addTerm(Monomial(1), ModInt(1));

    while (n > 0) {
        if (n & 1) {
            ret = divide(ret * a, {mod}).R;
        }
        a = divide(a * a, {mod}).R;
        n >>= 1;
    }
    // return modPolynomial(result, p);
    return ret;
}
template <class Order>
Polynomial<ModInt, Order> frobenius(const Polynomial<ModInt, Order>& g, const Polynomial<ModInt, Order>& f) {
    // only for univariate polynomials
    Polynomial<ModInt, Order> h;
    int p = ModInt::modulus();
    for (auto it = g.terms.begin(); it != g.terms.end(); it++) {
        Monomial m = it->first;
        m.set_exp(0, m.exp(0) * p);
        h.addTerm(m, it->second);
    }
    return polynomialRemainder(h, f);
}

}  // namespace polynomial
