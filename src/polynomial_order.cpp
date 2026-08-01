#include "Polynomial/monomial_order.hpp"

namespace polynomial {

bool LexOrder::operator()(const Monomial& a, const Monomial& b) const {
    assert(a.size() == b.size());
    for (size_t i = 0; i < a.size(); i++) {
        if (a.exp(i) != b.exp(i)) return a.exp(i) < b.exp(i);
    }
    return false;
}
bool GrLexOrder::operator()(const Monomial& a, const Monomial& b) const {
    assert(a.size() == b.size());
    auto da = a.degree();
    auto db = b.degree();
    if (da != db) return da < db;
    for (size_t i = 0; i < a.size(); i++) {
        if (a.exp(i) != b.exp(i)) return a.exp(i) < b.exp(i);
    }
    return false;
}
bool GrevLexOrder::operator()(const Monomial& a, const Monomial& b) const {
    assert(a.size() == b.size());
    auto da = a.degree();
    auto db = b.degree();
    if (da != db) return da < db;
    for (int i = (int)a.size() - 1; i >= 0; i--) {
        if (a.exp(i) != b.exp(i)) return a.exp(i) > b.exp(i);
    }
    return false;
}

}  // namespace polynomial
