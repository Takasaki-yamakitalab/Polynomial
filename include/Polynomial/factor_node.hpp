#ifndef POLYNOMIAL_FACTOR_NODE_HPP
#define POLYNOMIAL_FACTOR_NODE_HPP

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
template <class Integer, class Order>
std::unique_ptr<FactorNode<Integer, Order>> buildFactorTree(const std::vector<Polynomial<ModInt, Order>>& factors,
                                                            size_t begin, size_t end);
template <class Int1, class Int2, class Order>
void collectLeaves(FactorNode<Int1, Order>& node, std::vector<Polynomial<Int2, Order>>& ret);
template <class Coef, class Order>
void collectLeaves(FactorNode<ModInt, Order>& node, std::vector<Polynomial<Coef, Order>>& ret);


template <class Coef, class Order>
std::unique_ptr<FactorNode<Coef, Order>> makeLeaf(const Polynomial<Coef, Order>& f) {
    auto node = std::make_unique<FactorNode<Coef, Order>>();
    node->factor = f;
    return node;
}
template <class Coef, class Order>
std::unique_ptr<FactorNode<Coef, Order>> buildFactorTree(const std::vector<Polynomial<Coef, Order>>& factors,
                                                         size_t begin, size_t end) {
    // std::cout << "[DEBUG] buildFactorTree(" << begin << ", " << end << ")" << std::endl;
    if (begin >= end) {
        throw std::logic_error("buildFactorTree: empty range");
    }
    assert(begin < end);
    using FNode = FactorNode<Coef, Order>;
    auto node = std::make_unique<FNode>();
    if (end - begin == 1) {
        node->factor = factors[begin];
        return node;
    }
    size_t mid = (begin + end) / 2;
    node->left = buildFactorTree(factors, begin, mid);
    node->right = buildFactorTree(factors, mid, end);
    // std::cout << "[DEBUG] before multiply" << std::endl;
    node->factor = node->left->factor * node->right->factor;
    // std::cout << "[DEBUG] before extendedEuclid" << std::endl;
    auto ee = extendedEuclid(node->left->factor, node->right->factor);
    // std::cout << "[DEBUG] after extendedEuclid" << std::endl;
    node->s = ee.s;
    node->t = ee.t;
    return node;
}
template <class Integer, class Order>
std::unique_ptr<FactorNode<Integer, Order>> buildFactorTree(const std::vector<Polynomial<ModInt, Order>>& factors,
                                                            size_t begin, size_t end) {
    // std::cout << "[DEBUG] buildFactorTree(" << begin << ", " << end << ")" << std::endl;
    if (begin >= end) {
        throw std::logic_error("buildFactorTree: empty range");
    }
    assert(begin < end);
    using FNode = FactorNode<Integer, Order>;
    auto node = std::make_unique<FNode>();
    if (end - begin == 1) {
        node->factor.clear();
        for (auto it = factors[begin].terms.begin(); it != factors[begin].terms.end(); it++) {
            node->factor.addTerm(it->first, it->second.value());
        }
        return node;
    }
    size_t mid = (begin + end + 1) / 2;
    int p = factors[0].leadingCoefficient().modulus();
    node->left = buildFactorTree<Integer, Order>(factors, begin, mid);
    node->right = buildFactorTree<Integer, Order>(factors, mid, end);
    node->factor = node->left->factor * node->right->factor;
    auto ee = extendedEuclidMod(node->left->factor, node->right->factor, p);
    node->s = ee.s;
    node->t = ee.t;
    return node;
}
template <class Int1, class Int2, class Order>
void collectLeaves(FactorNode<Int1, Order>& node, std::vector<Polynomial<Int2, Order>>& ret) {
    if (node.isLeaf()) {
        Polynomial<Int2, Order> poly;
        for (auto it = node.factor.terms.begin(); it != node.factor.terms.end(); it++) {
            poly.addTerm(it->first, Int2(it->second));
        }
        ret.push_back(poly);
        return;
    }
    collectLeaves(*node.left, ret);
    collectLeaves(*node.right, ret);
}
template <class Integer, class Order>
void collectLeaves(FactorNode<ModInt, Order>& node, std::vector<Polynomial<Integer, Order>>& ret) {
    if (node.isLeaf()) {
        Polynomial<Integer, Order> poly;
        for (auto it = node.factor.terms.begin(); it != node.factor.terms.end(); it++) {
            poly.addTerm(it->first, Integer(it->second.value()));
        }
        ret.push_back(poly);
        return;
    }
    collectLeaves(*node.left, ret);
    collectLeaves(*node.right, ret);
}
}  // namespace polynomial

#endif
