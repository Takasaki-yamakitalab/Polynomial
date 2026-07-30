#include <Polynomial/polynomial.hpp>
#include <Polynomial/factor.hpp>

#include <iostream>

int main() {
    polynomial::Polynomial<int> poly;
    poly.addTerm(polynomial::Monomial(1, 3), 1);
    poly.addTerm(polynomial::Monomial(1, 2), 6);
    poly.addTerm(polynomial::Monomial(1, 1), 11);
    poly.addTerm(polynomial::Monomial(1, 0), 6);
    // poly.addTerm(polynomial::Monomial(1, 2), 1);
    // poly.addTerm(polynomial::Monomial(1, 1), 5);
    // poly.addTerm(polynomial::Monomial(1, 0), 6);
    // poly.addTerm(polynomial::Monomial(std::vector<int64_t>{2, 0}), 2);
    // poly.addTerm(polynomial::Monomial(std::vector<int64_t>{1, 1}), 7);
    // poly.addTerm(polynomial::Monomial(std::vector<int64_t>{0, 2}), 3);
    std::cout << poly << std::endl;
    auto ans = polynomial::factorize(poly);
    for (size_t i = 0; i < ans.size(); i++) {
        std::cout << "(" << ans[i].first << ")^(" << ans[i].second << ")";
        if (i + 1 < ans.size()) std::cout << "*";
    }
    std::cout << std::endl;
    return 0;
}
