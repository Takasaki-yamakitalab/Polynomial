#ifndef POLYNOMIAL_HPP
#define POLYNOMIAL_HPP

#include <map>
#include "monomial.hpp"
#include "monomial_order.hpp"
#include "term.hpp"
#include "modint.hpp"

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
    const Monomial& leadingMonomial() const;
    const Coef& leadingCoefficient() const;
    int64_t leadingDegree() const;
    Coef coefficient(const Monomial& m) const;
    void eraseLeadingTerm();
    Term<Coef> popLeadingTerm();
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
    Polynomial& operator+=(const Polynomial& rhs);
    Polynomial& operator+=(const Term<Coef>& rhs);
    Polynomial& operator-=(const Polynomial& rhs);
    Polynomial& operator*=(const Coef& c);
    Polynomial& operator*=(const Term<Coef>& t);
    Polynomial& operator*=(const Polynomial& rhs);
    Polynomial& operator/=(const Coef& c);
    Polynomial operator-() const;
};

template <class Coef, class Order = GrevLexOrder>
struct DivisionResult {
    std::vector<Polynomial<Coef, Order>> Q;
    Polynomial<Coef, Order> R;
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
Polynomial<Coef, Order> differential(const Polynomial<Coef, Order>& p, int var);
template <class Coef, class Order>
DivisionResult<Coef, Order> divide(Polynomial<Coef, Order> f, const std::vector<Polynomial<Coef, Order>>& G);
template <class Coef, class Order>
Polynomial<Coef, Order> primitivePart(Polynomial<Coef, Order> p);
template <class Coef, class Order>
Polynomial<Coef, Order> pow(Polynomial<Coef, Order> p, int64_t e);
template <class Coef, class Order>
Polynomial<Coef, Order> substitute(const Polynomial<Coef, Order>& f, int var, const Polynomial<Coef, Order>& g);
template <class Coef, class Order>
Polynomial<Coef, Order> substitute(const Polynomial<Coef, Order>& f, int var, const Coef& value);
template <class Coef, class Order>
Polynomial<Coef, Order> substitute(Polynomial<Coef, Order> f, const std::vector<Polynomial<Coef, Order>>& g);
template <class Coef, class Order>
DivisionResult<Coef, Order> pseudoDivide(Polynomial<Coef, Order> f, const Polynomial<Coef, Order>& g);
template <class Coef, class Order>
Polynomial<Coef, Order> pseudoRemainder(const Polynomial<Coef, Order>& f, const Polynomial<Coef, Order>& g);
template <class Coef, class Order>
Polynomial<Coef, Order> polynomialRemainder(Polynomial<Coef, Order> f, const Polynomial<Coef, Order>& g);

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

template <class Coef, class Order>
Polynomial<ModInt, Order> modPolynomial(const Polynomial<Coef, Order>& f, int p);
template <class Order>
Polynomial<ModInt, Order> derivativeMod(const Polynomial<int, Order>& f, int p);
template <class Order>
int findGoodPrime(const Polynomial<int, Order>& f);
template <class Order>
Polynomial<ModInt, Order> powModPolynomial(Polynomial<ModInt, Order> a, int64_t n,
                                           const Polynomial<ModInt, Order>& mod);
template <class Order>
Polynomial<ModInt, Order> frobenius(const Polynomial<ModInt, Order>& g, const Polynomial<ModInt, Order>& f);

}  // namespace polynomial

#endif
