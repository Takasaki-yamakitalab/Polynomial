#ifndef ZASSENHAUS_RECOMBINE_HPP
#define ZASSENHAUS_RECOMBINE_HPP

#include "polynomial.hpp"

namespace polynomial {

template <class Integer, class Order>
Polynomial<Integer, Order> centerLift(const Polynomial<Integer, Order>& f, int64_t modulus);

template <class Integer, class Order>
std::vector<Polynomial<Integer, Order>> centerLiftAll(const std::vector<Polynomial<Integer, Order>>& F,
                                                      int64_t modulus);

void enumerateSubset(int pos, int need, std::vector<int>& cur, std::vector<std::vector<int>>& ans, int n);

template <class Coef, class Order>
Polynomial<Coef, Order> productSubset(const std::vector<Polynomial<Coef, Order>>& factors,
                                      const std::vector<int>& subset);

template <class Coef, class Order>
std::vector<Polynomial<Coef, Order>> zassenhausRecombine(Polynomial<Coef, Order> f,
                                                         std::vector<Polynomial<Coef, Order>> lifted);

}  // namespace polynomial

#endif
