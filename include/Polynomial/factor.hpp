#ifndef FACTOR_HPP
#define FACTOR_HPP

#include <vector>
#include "polynomial.hpp"

namespace polynomial {

template <class Coef, class Order>
std::vector<SquareFreeFactor<Coef, Order>> squareFreeDecomposition(Polynomial<Coef, Order> f);

template <class Order>
std::vector<std::vector<ModInt>> berlekampMatrix(const Polynomial<ModInt, Order>& f);

std::vector<std::vector<ModInt>> nullSpace(std::vector<std::vector<ModInt>> A);

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

}  // namespace polynomial

#endif
