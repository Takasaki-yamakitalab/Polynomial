#ifndef COEFFICIENT_HPP
#define COEFFICIENT_HPP

#include <vector>

namespace polynomial {

template <class Coef>
Coef coefficientGCD(const Coef& a, const Coef& b);
template <class Coef>
Coef content(const std::vector<Coef>& c);

int64_t permutation(const int64_t& n, const int64_t& r);
int64_t combination(const int64_t& n, const int64_t& r);

template <class Coef>
Coef pow(Coef c, int64_t n);

}  // namespace polynomial

#endif
