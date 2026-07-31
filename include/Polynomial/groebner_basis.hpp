#ifndef POLYNOMIAL_GROEBNER_BASIS_HPP
#define POLYNOMIAL_GROEBNER_BASIS_HPP

#include "polynomial.hpp"
#include <deque>

namespace polynomial {

template <class Coef, class Order>
std::vector<Polynomial<Coef, Order>> buchberger(std::vector<Polynomial<Coef, Order>> g);

template <class Coef, class Order>
std::vector<Polynomial<Coef, Order>> reduceGroebnerBasis(std::vector<Polynomial<Coef, Order>> g);

template <class Coef, class Order>
std::vector<Polynomial<Coef, Order>> groebnerBasis(std::vector<Polynomial<Coef, Order>> g);

template <class Coef, class Order>
std::vector<Polynomial<Coef, Order>> buchberger(std::vector<Polynomial<Coef, Order>> G) {
    using Poly = Polynomial<Coef, Order>;
    std::deque<std::pair<int, int>> pairs;
    for (int i = 0; i < (int)G.size(); i++) {
        for (int j = i + 1; j < (int)G.size(); j++) {
            pairs.emplace_back(i, j);
        }
    }
    while (!pairs.empty()) {
        int i = pairs.front().first;
        int j = pairs.front().second;
        pairs.pop_front();

        if (gcd(G[i].leadingMonomial(), G[j].leadingMonomial()).isOne()) continue;

        Poly S = sPolynomial(G[i], G[j]);
        S.normalForm(G);
        if (S.isZero()) continue;

        bool exist = false;
        for (auto& g : G) {
            if (g == S) exist = true;
        }
        if (exist) continue;

        int k = (int)G.size();
        G.push_back(S);
        for (int t = 0; t < k; t++) {
            pairs.emplace_back(t, k);
        }
    }
    return G;
}

template <class Coef, class Order>
std::vector<Polynomial<Coef, Order>> reduceGroebnerBasis(std::vector<Polynomial<Coef, Order>> G) {
    using Poly = Polynomial<Coef, Order>;

    for (size_t i = 0; i < G.size();) {
        std::vector<Poly> others;
        for (size_t j = 0; j < G.size(); j++) {
            if (i != j) others.push_back(G[j]);
        }
        G[i].normalForm(others);
        if (G[i].isZero()) {
            G.erase(G.begin() + i);
            continue;
        }
        i++;
    }

    // make monic
    for (auto& g : G) {
        g *= Coef(1) / g.leadingCoefficient();
    }
    return G;
}

template <class Coef, class Order>
std::vector<Polynomial<Coef, Order>> groebnerBasis(std::vector<Polynomial<Coef, Order>> G) {
    G = buchberger(G);
    G = reduceGroebnerBasis(G);
    return G;
}
}

#endif
