#ifndef FACTOR_NODE_HPP
#define FACTOR_NODE_HPP

#include "polynomial.hpp"

namespace polynomial {

template <class Coef, class Order>
struct FactorNode {
    Polynomial<Coef, Order> factor;
    std::unique_ptr<FactorNode> left;
    std::unique_ptr<FactorNode> right;
    Polynomial<Coef, Order> s;
    Polynomial<Coef, Order> t;

    bool isLeaf() const { return !left && !right; }
};

template <class Coef, class Order>
std::unique_ptr<FactorNode<Coef, Order>> makeLeaf(const Polynomial<Coef, Order>& f);
template <class Coef, class Order>
std::unique_ptr<FactorNode<Coef, Order>> buildFactorTree(const std::vector<Polynomial<Coef, Order>>& factors,
                                                         size_t begin, size_t end);
template <class Order>
std::unique_ptr<FactorNode<int64_t, Order>> buildFactorTree(const std::vector<Polynomial<ModInt, Order>>& factors,
                                                             size_t begin, size_t end);
template <class Coef, class Order>
void collectLeaves(FactorNode<Coef, Order>& node, std::vector<Polynomial<Coef, Order>>& ret);

}  // namespace polynomial

#endif
