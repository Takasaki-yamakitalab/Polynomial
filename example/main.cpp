#include <Polynomial/polynomial.hpp>
#include <Polynomial/factor.hpp>
#include <Polynomial/zassenhaus_recombine.hpp>

#include <iostream>
#include "Polynomial/monomial.hpp"

int main() {
    // std::vector<std::vector<int>> dst;
    // std::vector<bool> exclude(10, false);
    // exclude[4] = true;
    // exclude[5] = true;
    // exclude[6] = true;
    // std::cout << "exclude: ";
    // for (auto b : exclude) {
    //     std::cout << b << " ";
    // }
    // std::cout << std::endl;
    // polynomial::enumerateSubset(dst, exclude, 3);
    // for (auto& v : dst) {
    //     std::cout << "subset: ";
    //     for (auto& i : v) std::cout << i << " ";
    //     std::cout << std::endl;
    // }
    // return 0;
    polynomial::Polynomial<int> poly1;
    poly1.addTerm(polynomial::Monomial(1, 4), -1);
    poly1.addTerm(polynomial::Monomial(1, 0), 1);
    polynomial::Polynomial<int> poly2;
    poly2.addTerm(polynomial::Monomial(1, 3), 1);
    poly2.addTerm(polynomial::Monomial(1, 2), 6);
    poly2.addTerm(polynomial::Monomial(1, 1), 11);
    poly2.addTerm(polynomial::Monomial(1, 0), 6);
    polynomial::Polynomial<int> poly3;
    poly3.addTerm(polynomial::Monomial(std::vector<int64_t>{2, 0}), 1);
    poly3.addTerm(polynomial::Monomial(std::vector<int64_t>{0, 2}), -1);
    polynomial::Polynomial<int> poly4;
    poly4.addTerm(polynomial::Monomial(std::vector<int64_t>{1, 0, 1, 0}), 1);
    poly4.addTerm(polynomial::Monomial(std::vector<int64_t>{1, 0, 0, 1}), 1);
    poly4.addTerm(polynomial::Monomial(std::vector<int64_t>{0, 1, 1, 0}), 1);
    poly4.addTerm(polynomial::Monomial(std::vector<int64_t>{0, 1, 0, 1}), 1);
    // poly4.addTerm(polynomial::Monomial(1, 5), 1);
    // poly4.addTerm(polynomial::Monomial(1, 4), 1);
    // poly4.addTerm(polynomial::Monomial(1, 1), 1);
    // poly4.addTerm(polynomial::Monomial(1, 0), 1);

    std::cout << "TEST 1" << std::endl;
    auto ans1 = polynomial::factorize(poly1);
    std::cout << "ORIGIN: " << poly1 << std::endl;
    std::cout << "FACTORIZED: ";
    for (size_t i = 0; i < ans1.size(); i++) {
        std::cout << "(" << ans1[i].first << ")^(" << ans1[i].second << ")";
        if (i + 1 < ans1.size()) std::cout << "*";
    }
    std::cout << std::endl;

    std::cout << "TEST 2" << std::endl;
    auto ans2 = polynomial::factorize(poly2);
    std::cout << "ORIGIN: " << poly2 << std::endl;
    std::cout << "FACTORIZED: ";
    for (size_t i = 0; i < ans2.size(); i++) {
        std::cout << "(" << ans2[i].first << ")^(" << ans2[i].second << ")";
        if (i + 1 < ans2.size()) std::cout << "*";
    }
    std::cout << std::endl;

    std::cout << "TEST 3" << std::endl;
    auto ans3 = polynomial::factorize(poly3);
    std::cout << "ORIGIN: " << poly3 << std::endl;
    std::cout << "FACTORIZED: ";
    for (size_t i = 0; i < ans3.size(); i++) {
        std::cout << "(" << ans3[i].first << ")^(" << ans3[i].second << ")";
        if (i + 1 < ans3.size()) std::cout << "*";
    }
    std::cout << std::endl;

    std::cout << "TEST 4" << std::endl;
    auto ans4 = polynomial::factorize(poly4);
    std::cout << "ORIGIN: " << poly4 << std::endl;
    std::cout << "FACTORIZED: ";
    for (size_t i = 0; i < ans4.size(); i++) {
        std::cout << "(" << ans4[i].first << ")^(" << ans4[i].second << ")";
        if (i + 1 < ans4.size()) std::cout << "*";
    }
    std::cout << std::endl;
    return 0;
}
