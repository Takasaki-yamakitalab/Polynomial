#include "Polynomial/factor_node.hpp"

namespace polynomial {

template <class Coef, class Order>
std::unique_ptr<FactorNode<Coef, Order>> makeLeaf(const Polynomial<Coef, Order>& f) {
    auto node = std::make_unique<FactorNode<Coef, Order>>();
    node->factor = f;
    return node;
}

template <class Coef, class Order>
std::unique_ptr<FactorNode<Coef, Order>> buildFactorTree(const std::vector<Polynomial<Coef, Order>>& factors,
                                                            size_t begin, size_t end) {
    using FNode = FactorNode<Coef, Order>;
    auto node = std::make_unique<FNode>();
    if (end - begin == 1) {
        node->factor = factors[begin];
        return node;
    }
    size_t mid = (begin + end) / 2;
    node->left = buildFactorTree(factors, begin, mid);
    node->right = buildFactorTree(factors, mid, end);
    node->factor = node->left->factor * node->right->factor;
    auto ee = extendedEuclid(node->left->factor, node->right->factor);
    node->s = ee.s;
    node->t = ee.t;
    return node;
}
template <class Order>
std::unique_ptr<FactorNode<int64_t, Order>> buildFactorTree(const std::vector<Polynomial<ModInt, Order>>& factors,
                                                            size_t begin, size_t end) {
    using FNode = FactorNode<int64_t, Order>;
    auto node = std::make_unique<FNode>();
    if (end - begin == 1) {
        node->factor.clear();
        for (auto it = factors[begin].terms.begin(); it != factors[begin].terms.end(); it++) {
            node->factor.addTerm(it->first, it->second.value());
        }
        return node;
    }
    size_t mid = (begin + end) / 2;
    node->left = buildFactorTree(factors, begin, mid);
    node->right = buildFactorTree(factors, mid, end);
    node->factor = node->left->factor * node->right->factor;
    auto ee = extendedEuclid(node->left->factor, node->right->factor);
    node->s = ee.s;
    node->t = ee.t;
    return node;
}

template <class Coef, class Order>
void collectLeaves(FactorNode<Coef, Order>& node, std::vector<Polynomial<Coef, Order>>& ret) {
    if (node.isLeaf()) {
        ret.push_back(node.factor);
        return;
    }
    collectLeaves(*node.left, ret);
    collectLeaves(*node.right, ret);
}

}  // namespace polynomial
