#ifndef POLYNOMIAL_COEFFICIENT_TRAITS_HPP
#define POLYNOMIAL_COEFFICIENT_TRAITS_HPP

#include <cstdint>

namespace polynomial {

template <class Coef>
struct coefficient_traits {
    static const bool is_field = false;
};

template <>
struct coefficient_traits<double> {
    static const bool is_field = true;
};

template <>
struct coefficient_traits<int> {
    static const bool is_field = false;
};

template <>
struct coefficient_traits<int64_t> {
    static const bool is_field = false;
};

}  // namespace polynomial
#endif
