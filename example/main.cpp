#include <Polynomial/polynomial.hpp>
#include <Polynomial/factor.hpp>

#include <iostream>

int main() {
    polynomial::Polynomial<int> poly1;
    poly1.addTerm(polynomial::Monomial(1, 2), 1);
    poly1.addTerm(polynomial::Monomial(1, 1), 5);
    poly1.addTerm(polynomial::Monomial(1, 0), 6);
    polynomial::Polynomial<int> poly2;
    poly2.addTerm(polynomial::Monomial(1, 3), 1);
    poly2.addTerm(polynomial::Monomial(1, 2), 6);
    poly2.addTerm(polynomial::Monomial(1, 1), 11);
    poly2.addTerm(polynomial::Monomial(1, 0), 6);
    polynomial::Polynomial<int> poly3;
    poly3.addTerm(polynomial::Monomial(std::vector<int64_t>{2, 0}), 2);
    poly3.addTerm(polynomial::Monomial(std::vector<int64_t>{1, 1}), 7);
    poly3.addTerm(polynomial::Monomial(std::vector<int64_t>{0, 2}), 3);

    std::cout << poly1 << std::endl;
    auto ans1 = polynomial::factorize(poly1);
    for (size_t i = 0; i < ans1.size(); i++) {
        std::cout << "(" << ans1[i].first << ")^(" << ans1[i].second << ")";
        if (i + 1 < ans1.size()) std::cout << "*";
    }
    std::cout << std::endl;

    // std::cout << poly2 << std::endl;
    // auto ans2 = polynomial::factorize(poly2);
    // for (size_t i = 0; i < ans2.size(); i++) {
    //     std::cout << "(" << ans2[i].first << ")^(" << ans2[i].second << ")";
    //     if (i + 1 < ans2.size()) std::cout << "*";
    // }
    // std::cout << std::endl;
    //
    // std::cout << poly3 << std::endl;
    // auto ans3 = polynomial::factorize(poly3);
    // for (size_t i = 0; i < ans3.size(); i++) {
    //     std::cout << "(" << ans3[i].first << ")^(" << ans3[i].second << ")";
    //     if (i + 1 < ans3.size()) std::cout << "*";
    // }
    // std::cout << std::endl;

    return 0;
}
