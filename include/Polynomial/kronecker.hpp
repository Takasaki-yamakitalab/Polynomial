#ifndef KRONECKER_HPP
#define KRONECKER_HPP

#include "polynomial.hpp"

namespace polynomial {

template <class Coef, class Order>
int64_t kroneckerBase(const Polynomial<Coef, Order>& f);

template <class Coef, class Order>
Polynomial<Coef, LexOrder> kroneckerSubstittution(const Polynomial<Coef, Order>& f, int64_t B);

template <class Coef>
Polynomial<Coef, LexOrder> inverseKronecker(const Polynomial<Coef, LexOrder>& f, int64_t B, size_t var_num);

template <class Coef>
std::vector<Polynomial<Coef, LexOrder>> inverseKroneckerAll(const std::vector<Polynomial<Coef, LexOrder>>& f, int64_t B, size_t var_num);

}

#endif
