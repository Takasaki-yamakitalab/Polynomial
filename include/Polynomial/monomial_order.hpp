#ifndef POLYNOMIAL_MONOMIAL_ORDER_HPP
#define POLYNOMIAL_MONOMIAL_ORDER_HPP

#include "monomial.hpp"

namespace polynomial {

struct LexOrder {
    bool operator()(const Monomial& a, const Monomial& b) const;
};
struct GrLexOrder {
    bool operator()(const Monomial& a, const Monomial& b) const;
};
struct GrevLexOrder {
    bool operator()(const Monomial& a, const Monomial& b) const;
};

}  // namespace polynomial

#endif
