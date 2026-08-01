#ifndef POLYNOMIAL_POLYNOMIAL_HPP
#define POLYNOMIAL_POLYNOMIAL_HPP

#include <map>
#include "monomial.hpp"
#include "monomial_order.hpp"
#include "term.hpp"
#include "modint.hpp"
#include "coefficient.hpp"
#include <algorithm>
#include <iostream>

namespace polynomial {

template <class Coef, class Order = LexOrder>
struct Polynomial {
    using Map = std::map<Monomial, Coef, Order>;
    Map terms;

    Polynomial() = default;
    Polynomial(const Term<Coef>& t) { addTerm(t); }
    Polynomial(const Monomial& m, const Coef& c) { addTerm(m, c); }

    bool isZero() const;
    bool isConstant() const;
    bool isMonomial() const;
    int64_t degree() const;
    int64_t degree(const int& var) const;
    Term<Coef> leadingTerm() const;
    Term<Coef> leadingTerm(const int& var) const;
    const Monomial& leadingMonomial() const;
    const Monomial& leadingMonomial(const int& var) const;
    const Coef& leadingCoefficient() const;
    const Coef& leadingCoefficient(const int& var) const;
    int64_t leadingDegree() const;
    Coef coefficient(const Monomial& m) const;
    void eraseLeadingTerm();
    void eraseLeadingTerm(const int& var);
    Term<Coef> popLeadingTerm();
    Term<Coef> popLeadingTerm(const int& var);
    Polynomial<Coef, Order> leadingPolynomial() const;
    void normalize();
    void clear();
    void addTerm(const Monomial& m, const Coef& c);
    void addTerm(const Term<Coef>& t);
    bool reducibleBy(const Polynomial& g) const;
    Polynomial& reduceOnceBy(const Polynomial& g);
    Polynomial& reduceBy(const Polynomial& g);
    Polynomial& normalForm(const std::vector<Polynomial>& G);

    bool operator==(const Polynomial& rhs) const;
    bool operator!=(const Polynomial& rhs) const;
    Polynomial& operator+=(const Polynomial& rhs);
    Polynomial& operator+=(const Term<Coef>& rhs);
    Polynomial& operator-=(const Polynomial& rhs);
    Polynomial& operator*=(const Coef& c);
    Polynomial& operator*=(const Term<Coef>& t);
    Polynomial& operator*=(const Polynomial& rhs);
    Polynomial& operator/=(const Coef& c);
    Polynomial operator-() const;

    friend std::ostream& operator<<(std::ostream& os, const Polynomial<Coef, Order>& f) {
        for (auto it = f.terms.begin(); it != f.terms.end();) {
            os << it->second << "*" << it->first;
            if (++it != f.terms.end()) os << "+";
        }
        return os;
    }
};

template <class Coef, class Order = LexOrder>
struct DivisionResult {
    std::vector<Polynomial<Coef, Order>> Q;
    Polynomial<Coef, Order> R;
    Coef multiplier = Coef(1);
};

template <class Coef, class Order>
struct SquareFreeFactor {
    Polynomial<Coef, Order> poly;
    uint64_t multiplicity;
};

template <class Coef, class Order>
struct ExtendEuclidResult {
    Polynomial<Coef, Order> gcd;
    Polynomial<Coef, Order> s;
    Polynomial<Coef, Order> t;
};

template <class Coef, class Order>
Polynomial<Coef, Order> operator+(const Polynomial<Coef, Order>& a, const Polynomial<Coef, Order>& b);
template <class Coef, class Order>
Polynomial<Coef, Order> operator-(const Polynomial<Coef, Order>& a, const Polynomial<Coef, Order>& b);
template <class Coef, class Order>
Polynomial<Coef, Order> operator*(const Polynomial<Coef, Order>& a, const Polynomial<Coef, Order>& b);
template <class Coef, class Order>
Polynomial<Coef, Order> operator*(Polynomial<Coef, Order> p, const Term<Coef>& t);
template <class Coef, class Order>
Polynomial<Coef, Order> operator*(Polynomial<Coef, Order> p, const Coef& c);
template <class Coef, class Order>
Polynomial<Coef, Order> operator*(const Coef& c, Polynomial<Coef, Order> p);
template <class Coef, class Order>
Polynomial<Coef, Order> operator/(Polynomial<Coef, Order> p, const Coef& c);

template <class Coef, class Order>
Polynomial<Coef, Order> differential(const Polynomial<Coef, Order>& p, int var);
template <class Coef, class Order>
DivisionResult<Coef, Order> divide(Polynomial<Coef, Order> f, const std::vector<Polynomial<Coef, Order>>& G);
template <class Integer, class Order>
DivisionResult<Integer, Order> divideMod(Polynomial<Integer, Order> f, const std::vector<Polynomial<Integer, Order>>& G,
                                         const int p);
template <class Coef, class Order>
Polynomial<Coef, Order> primitivePart(Polynomial<Coef, Order> p);
template <class Integer, class Order>
Polynomial<Integer, Order> primitivePartMultiVar(Polynomial<Integer, Order> p, const int& var);
template <class Coef, class Order>
Polynomial<Coef, Order> pow(Polynomial<Coef, Order> p, int64_t e);
template <class Coef, class Order>
Polynomial<Coef, Order> substitute(const Polynomial<Coef, Order>& f, int var, const Polynomial<Coef, Order>& g);
template <class Coef, class Order>
Polynomial<Coef, Order> substitute(const Polynomial<Coef, Order>& f, int var, const Coef& value);
template <class Coef, class Order>
Polynomial<Coef, Order> substitute(Polynomial<Coef, Order> f, const std::vector<Polynomial<Coef, Order>>& g);
template <class Integer, class Order>
DivisionResult<Integer, Order> pseudoDivide(Polynomial<Integer, Order> f, const Polynomial<Integer, Order>& g);

template <class Integer, class Order>
Polynomial<Integer, Order> pseudoRemainder(const Polynomial<Integer, Order>& f, const Polynomial<Integer, Order>& g);
template <class Integer, class Order>
Polynomial<Integer, Order> pseudoRemainder(const Polynomial<Integer, Order>& f, const Polynomial<Integer, Order>& g,
                                           const int& var);

template <class Coef, class Order>
Polynomial<Coef, Order> polynomialRemainder(Polynomial<Coef, Order> f, const Polynomial<Coef, Order>& g);

template <class Coef, class Order>
size_t selectLowDegVar(const std::vector<Polynomial<Coef, Order>>& F);
template <class Coef, class Order>
Polynomial<Coef, Order> polynomialGCD(Polynomial<Coef, Order> f, Polynomial<Coef, Order> g);

template <class Coef, class Order>
Polynomial<Coef, Order> sPolynomial(const Polynomial<Coef, Order>& f, const Polynomial<Coef, Order>& g);

template <class Coef, class Order>
int findReducer(const std::vector<Polynomial<Coef, Order>>& G, const Monomial& m);

template <class Coef, class Order>
Polynomial<Coef, Order> multiDerivativeGCD(Polynomial<Coef, Order> f);

template <class Coef, class Order>
DivisionResult<Coef, Order> divideByLinearFactor(const Polynomial<Coef, Order>& p, int var, Coef a);
template <class Coef, class Order>
ExtendEuclidResult<Coef, Order> extendedEuclid(Polynomial<Coef, Order> f, Polynomial<Coef, Order> g);
template <class Integer, class Order>
ExtendEuclidResult<Integer, Order> extendedEuclidMod(Polynomial<Integer, Order> f, Polynomial<Integer, Order> g, int p);

template <class Coef, class Order>
Polynomial<Coef, Order> product(const std::vector<Polynomial<Coef, Order>>& factors);
template <class Coef, class Order>
Polynomial<Coef, Order> product(const std::vector<Polynomial<Coef, Order>>& factors, size_t begin, size_t end);

template <class Coef, class Order>
int64_t landauMignotteBound(const Polynomial<Coef, Order>& f);
template <class Coef, class Order>
int64_t henselBound(const Polynomial<Coef, Order>& f);

template <class Coef, class Order>
bool isFactor(const Polynomial<Coef, Order>& f, const Polynomial<Coef, Order>& g);

template <class Coef, class Order>
Coef content(const Polynomial<Coef, Order>& f);
template <class Integer, class Order>
Polynomial<Integer, Order> contentMultiVar(const Polynomial<Integer, Order>& f, size_t var);

template <class Coef, class Order>
Polynomial<ModInt, Order> modPolynomial(const Polynomial<Coef, Order>& f, int p);
template <class Order>
Polynomial<ModInt, Order> modPolynomial(const Polynomial<ModInt, Order>& f, int p);
template <class Order>
Polynomial<ModInt, Order> derivativeMod(const Polynomial<int, Order>& f, int p);
template <class Order>
int findGoodPrime(const Polynomial<int, Order>& f);
template <class Order>
Polynomial<ModInt, Order> powModPolynomial(Polynomial<ModInt, Order> a, int64_t n,
                                           const Polynomial<ModInt, Order>& mod);
template <class Order>
Polynomial<ModInt, Order> frobenius(const Polynomial<ModInt, Order>& g, const Polynomial<ModInt, Order>& f);

template <class Coef, class Order>
Monomial monomialContent(const Polynomial<Coef, Order>& f);

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
        int64_t n = it->first.degree();
        if (n > deg) deg = n;
    }
    return deg;
}
template <class Coef, class Order>
int64_t Polynomial<Coef, Order>::degree(const int& var) const {
    if (isZero()) return -1;
    int64_t deg = 0;
    for (auto it = terms.begin(); it != terms.end(); it++) {
        int64_t n = it->first.exp(var);
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
bool Polynomial<Coef, Order>::operator!=(const Polynomial<Coef, Order>& rhs) const {
    return terms != rhs.terms;
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
    p *= c;
    return p;
}
template <class Coef, class Order>
Polynomial<Coef, Order> operator*(const Coef& c, Polynomial<Coef, Order> p) {
    p *= c;
    return p;
}
template <class Coef, class Order>
Polynomial<Coef, Order> operator/(Polynomial<Coef, Order> p, const Coef& c) {
    p /= c;
    return p;
}

template <class Coef, class Order>
Polynomial<Coef, Order> differential(const Polynomial<Coef, Order>& p, int var) {
    Polynomial<Coef, Order> ret;
    for (auto& term : p.terms) {
        if (term.first.exp(var) == 0) continue;
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
template <class Integer, class Order>
DivisionResult<Integer, Order> divideMod(Polynomial<Integer, Order> f, const std::vector<Polynomial<Integer, Order>>& G,
                                         const int p) {
    std::vector<Polynomial<Integer, Order>> Q;
    Polynomial<Integer, Order> R;

    int mod = ModContext::modulus();
    Polynomial<ModInt, Order> fp = modPolynomial(f, p);
    std::vector<Polynomial<ModInt, Order>> Gp;
    for (auto& g : G) {
        Gp.push_back(modPolynomial(g, p));
    }

    auto result = divide(fp, Gp);
    for (auto& q : result.Q) {
        Polynomial<Integer, Order> poly;
        for (auto it = q.terms.begin(); it != q.terms.end(); it++) {
            poly.addTerm(it->first, it->second.value());
        }
        Q.push_back(poly);
    }
    for (auto it = result.R.terms.begin(); it != result.R.terms.end(); it++) {
        R.addTerm(it->first, it->second.value());
    }

    ModContext::setMod(mod);

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
    if (p.isZero()) return p;
    p = primitivePartImpl(p, std::integral_constant<bool, coefficient_traits<Coef>::is_field>());
    if (p.leadingCoefficient() < Coef(0)) p = -p;
    return p;
}
template <class Integer, class Order>
Polynomial<Integer, Order> primitivePartMultiVar(Polynomial<Integer, Order> p, const int& var) {
    if (p.isZero()) return p;
    auto cont = contentMultiVar(p, var);
    // if (cont.divides(p)) return p / cont;
    // return p;
    return divide(p, {cont}).Q[0];
}

template <class Integer, class Order>
DivisionResult<Integer, Order> pseudoDivide(Polynomial<Integer, Order> f, const Polynomial<Integer, Order>& g) {
    using Poly = Polynomial<Integer, Order>;
    DivisionResult<Integer, Order> ret;
    assert(!g.isZero());
    std::cout << "[DEBUG] pseudoDivide(" << f << ", " << g << ")" << std::endl;

    Poly R = f;
    Poly Q;

    if (f.degree() < g.degree()) {
        ret.R = f;
        ret.Q.push_back(Q);
        return ret;
    }

    Integer c = g.leadingCoefficient();
    Integer multiplier = Integer(1);
    while (!R.isZero() && R.degree() >= g.degree()) {
        if (!g.leadingMonomial().divides(R.leadingMonomial())) break;
        Monomial u = R.leadingMonomial() / g.leadingMonomial();
        Integer a = R.leadingCoefficient();
        Term<Integer> t{u, a};
        Q *= c;
        Q += t;
        R *= c;
        R -= g * t;
        multiplier *= c;
    }
    ret.Q.push_back(Q);
    ret.R = R;
    ret.multiplier = multiplier;
    return ret;
}
template <class Integer, class Order>
Polynomial<Integer, Order> pseudoRemainder(const Polynomial<Integer, Order>& f, const Polynomial<Integer, Order>& g) {
    return pseudoDivide(f, g).R;
}
template <class Integer, class Order>
Polynomial<Integer, Order> pseudoRemainder(const Polynomial<Integer, Order>& f, const Polynomial<Integer, Order>& g,
                                           const int& var) {
    using Poly = Polynomial<Integer, Order>;
    assert(!g.isZero());
    std::cout << "[DEBUG] pseudoRemainder(" << f << ", " << g << ", " << var << ")" << std::endl;

    Poly R = f;

    if (f.degree(var) < g.degree(var)) {
        return f;
    }

    Term<Integer> c = g.leadingTerm(var);
    c.mono.set_exp(var, 0);
    while (!R.isZero() && R.degree(var) >= g.degree(var)) {
        // if (!g.leadingMonomial(var).divides(R.leadingMonomial(var))) break;
        // Monomial u = R.leadingMonomial(var) / g.leadingMonomial(var);
        // Integer a = R.leadingCoefficient(var);
        // Term<Integer> t{u, a};
        Term<Integer> t = R.leadingTerm(var);
        t.mono.add_exp(var, -g.degree(var));
        R *= c;
        R -= g * t;
    }
    return R;
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
size_t selectLowDegVar(const std::vector<Polynomial<Coef, Order>>& F) {
    size_t var = 0;
    int64_t min_score = 0;

    if (F.empty()) return 0;

    size_t n = F[0].leadingMonomial().size();

    for (size_t k = 0; k < F.size(); k++) {
        min_score += F[k].degree(0);
    }
    for (size_t i = 1; i < n; i++) {
        size_t score = 0;
        for (auto& f : F) score += f.degree(i);
        if (score < min_score) {
            min_score = score;
            var = i;
        }
    }
    return var;
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
static Polynomial<Coef, Order> polynomialGCDImpl(Polynomial<Coef, Order> f, Polynomial<Coef, Order> g,
                                                 std::false_type) {
    // 整数係数
    using Poly = Polynomial<Coef, Order>;
    size_t var = selectLowDegVar(std::vector{f, g});
    if (g.isZero()) return primitivePartMultiVar(f, var);

    auto c = coefficientGCD(content(f), content(g));
    while (!g.isZero()) {
        Poly r = pseudoRemainder(f, g, var);
        std::cout << "[DEBUG] r: " << r << "\n";
        std::cout << "[DEBUG] isZero: " << r.isZero() << "\n";
        std::cout << "[DEBUG] size: " << r.terms.size() << "\n";
        r = primitivePartMultiVar(r, var);
        f = g;
        g = r;
    }
    return c * primitivePartMultiVar(f, var);
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
    s0.addTerm(Monomial(f.leadingMonomial().size()), Coef(1));
    // s1 = 0
    // t0 = 0
    // t1 = 1
    t1.addTerm(Monomial(f.leadingMonomial().size()), Coef(1));

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
template <class Integer, class Order>
ExtendEuclidResult<Integer, Order> extendedEuclidMod(Polynomial<Integer, Order> f, Polynomial<Integer, Order> g,
                                                     int p) {
    using Poly = Polynomial<Integer, Order>;
    Polynomial<ModInt, Order> fp = modPolynomial(f, p);
    Polynomial<ModInt, Order> gp = modPolynomial(g, p);
    ExtendEuclidResult<ModInt, Order> result = extendedEuclid(fp, gp);
    ExtendEuclidResult<Integer, Order> ret;
    for (auto it = result.gcd.terms.begin(); it != result.gcd.terms.end(); it++) {
        ret.gcd.addTerm(it->first, Integer(it->second.value()));
    }
    for (auto it = result.s.terms.begin(); it != result.s.terms.end(); it++) {
        ret.s.addTerm(it->first, Integer(it->second.value()));
    }
    for (auto it = result.t.terms.begin(); it != result.t.terms.end(); it++) {
        ret.t.addTerm(it->first, Integer(it->second.value()));
    }
    return ret;
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
        infNorm = std::max(infNorm, int64_t(std::abs(it->second)));
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
template <class Integer, class Order>
Polynomial<Integer, Order> contentMultiVar(const Polynomial<Integer, Order>& f, size_t var) {
    if (f.isZero()) return f;
    Integer c = content(f);
    const size_t N = f.leadingMonomial().size();
    Monomial gcd(f.leadingMonomial().size());
    for (size_t i = 0; i < N; i++) {
        if (i == var) continue;
        int64_t e = f.degree(i);
        for (auto it = f.terms.begin(); it != f.terms.end(); it++) {
            if (it->first[i] < e) e = it->first[i];
        }
        gcd.set_exp(i, e);
    }
    Polynomial<Integer, Order> cont(gcd, c);
    return cont;
}

template <class Coef, class Order>
Polynomial<ModInt, Order> modPolynomial(const Polynomial<Coef, Order>& f, int p) {
    ModContext::setMod(p);
    Polynomial<ModInt, Order> ret;
    for (auto it = f.terms.begin(); it != f.terms.end(); it++) {
        ModInt v((int)std::round(it->second));
        if (v.value() != 0) ret.addTerm(it->first, v);
    }
    return ret;
}
template <class Order>
Polynomial<ModInt, Order> modPolynomial(const Polynomial<ModInt, Order>& f, int p) {
    ModContext::setMod(p);
    Polynomial<ModInt, Order> ret;
    for (auto it = f.terms.begin(); it != f.terms.end(); it++) {
        ModInt v((int)std::round(it->second.value() + 0.5));
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

        std::cout << "[DEBUG] p: " << p << std::endl;
        std::cout << "[DEBUG] fp: " << fp << std::endl;
        std::cout << "[DEBUG] df: " << df << std::endl;
        auto g = polynomialGCD(fp, df);
        std::cout << "[DEBUG] gcd: " << g << std::endl;

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
    int p = ModContext::modulus();
    for (auto it = g.terms.begin(); it != g.terms.end(); it++) {
        Monomial m = it->first;
        m.set_exp(0, m.exp(0) * p);
        h.addTerm(m, it->second);
    }
    return polynomialRemainder(h, f);
}

template <class Coef, class Order>
Term<Coef> Polynomial<Coef, Order>::leadingTerm() const {
    assert(!terms.empty());
    auto it = std::prev(terms.end());
    if (std::is_same<Order, GrevLexOrder>::value) it = terms.begin();
    return {it->first, it->second};
}
template <class Coef, class Order>
Term<Coef> Polynomial<Coef, Order>::leadingTerm(const int& var) const {
    assert(!terms.empty());
    int64_t max_n = 0;
    auto ret = terms.begin();
    for (auto it = ret; it != terms.end(); it++) {
        if (it->first[var] > max_n) {
            ret = it;
            max_n = it->first[var];
        }
    }
    return {ret->first, ret->second};
}
template <class Coef, class Order>
const Monomial& Polynomial<Coef, Order>::leadingMonomial() const {
    assert(!terms.empty());
    if (std::is_same<Order, GrevLexOrder>::value) return terms.begin()->first;
    return std::prev(terms.end())->first;
}
template <class Coef, class Order>
const Monomial& Polynomial<Coef, Order>::leadingMonomial(const int& var) const {
    assert(!terms.empty());
    int64_t max_n = 0;
    auto ret = terms.begin();
    for (auto it = ret; it != terms.end(); it++) {
        if (it->first[var] > max_n) {
            ret = it;
            max_n = it->first[var];
        }
    }
    return ret->first;
}
template <class Coef, class Order>
const Coef& Polynomial<Coef, Order>::leadingCoefficient() const {
    assert(!terms.empty());
    if (std::is_same<Order, GrevLexOrder>::value) return terms.begin()->second;
    return std::prev(terms.end())->second;
}
template <class Coef, class Order>
const Coef& Polynomial<Coef, Order>::leadingCoefficient(const int& var) const {
    assert(!terms.empty());
    int64_t max_n = 0;
    auto ret = terms.begin();
    for (auto it = ret; it != terms.end(); it++) {
        if (it->first[var] > max_n) {
            ret = it;
            max_n = it->first[var];
        }
    }
    return ret->second;
}
template <class Coef, class Order>
void Polynomial<Coef, Order>::eraseLeadingTerm() {
    assert(!terms.empty());
    if (std::is_same<Order, GrevLexOrder>::value)
        terms.erase(terms.begin());
    else
        terms.erase(std::prev(terms.end()));
}
template <class Coef, class Order>
void Polynomial<Coef, Order>::eraseLeadingTerm(const int& var) {
    assert(!terms.empty());
    int64_t max_n = 0;
    auto ret = terms.begin();
    for (auto it = ret; it != terms.end(); it++) {
        if (it->first[var] > max_n) {
            ret = it;
            max_n = it->first[var];
        }
    }
    terms.erase(ret);
}
template <class Coef, class Order>
Term<Coef> Polynomial<Coef, Order>::popLeadingTerm() {
    auto it = std::prev(terms.end());
    if (std::is_same<Order, GrevLexOrder>::value) it = terms.begin();
    Term<Coef> t{it->first, it->second};
    terms.erase(it);
    return t;
}
template <class Coef, class Order>
Term<Coef> Polynomial<Coef, Order>::popLeadingTerm(const int& var) {
    int64_t max_n = 0;
    auto ret = terms.begin();
    for (auto it = ret; it != terms.end(); it++) {
        if (it->first[var] > max_n) {
            ret = it;
            max_n = it->first[var];
        }
    }
    Term<Coef> t{ret->first, ret->second};
    terms.erase(ret);
    return t;
}

template <class Coef, class Order>
Monomial monomialContent(const Polynomial<Coef, Order>& f) {
    Monomial ret(f.leadingMonomial().size());
    for (size_t i = 0; i < f.leadingMonomial().size(); i++) {
        int64_t min_n = f.degree();
        for (auto it = f.terms.begin(); it != f.terms.end(); it++) {
            if (it->first.exp(i) < min_n) min_n = it->first.exp(i);
        }
        ret.add_exp(i, min_n);
    }
    return ret;
}
}  // namespace polynomial

#endif
