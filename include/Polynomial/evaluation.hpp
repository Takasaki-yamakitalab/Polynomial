#ifndef EVALUATION_HPP
#define EVALUATION_HPP

#include "polynomial.hpp"

namespace polynomial {

template <class RNG, class Coef, class Order>
Polynomial<Coef, Order> randomEvaluation(const Polynomial<Coef, Order>& f, std::vector<Coef>& values, RNG& rng,
                                         int bound = 10);
template <class Coef, class Order>
bool evaluationIsGood(const Polynomial<Coef, Order>& f, const Polynomial<Coef, Order>& g);

}  // namespace polynomial

#endif
