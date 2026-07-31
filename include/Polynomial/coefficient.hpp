#ifndef POLYNOMIAL_COEFFICIENT_HPP
#define POLYNOMIAL_COEFFICIENT_HPP

#include <vector>

namespace polynomial {

template <class Coef>
Coef coefficientGCD(Coef a, Coef b);
template <class Coef>
Coef content(const std::vector<Coef>& c);

int64_t permutation(const int64_t& n, const int64_t& r);
int64_t combination(const int64_t& n, const int64_t& r);

template <class Coef>
Coef pow(Coef c, int64_t n);

// ===================================================

template <class Coef>
Coef coefficientGCD(Coef a, Coef b) {
    a = abs(a);
    b = abs(b);
    while (b != Coef(0)) {
        Coef r = fmod(a, b);
        a = b;
        b = r;
    }
    return a;
}
template <class Coef>
Coef content(const std::vector<Coef>& c) {
    assert(!c.empty());
    Coef g = abs(c.front());
    for (size_t i = 1; i < c.size(); i++) {
        g = coefficientGCD(g, abs(c[i]));
    }
    return g;
}
template <class Coef>
Coef pow(Coef c, int64_t n) {
    if (n == 0) return Coef(0);
    if (n == 1) return c;
    for (size_t i = 0; i < n; i++) {
        c *= c;
    }
    return c;
}

}  // namespace polynomial

#endif
