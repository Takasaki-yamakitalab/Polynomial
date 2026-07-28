#ifndef groebner_basis_hpp
#define groebner_basis_hpp

#include "polynomial.hpp"

namespace polynomial {

template <class Coef, class Order>
std::vector<Polynomial<Coef, Order>> buchberger(std::vector<Polynomial<Coef, Order>> g);

template <class Coef, class Order>
std::vector<Polynomial<Coef, Order>> reduceGroebnerBasis(std::vector<Polynomial<Coef, Order>> g);

template <class Coef, class Order>
std::vector<Polynomial<Coef, Order>> groebnerBasis(std::vector<Polynomial<Coef, Order>> g);

}

#endif
