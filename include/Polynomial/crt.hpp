#ifndef POLYNOMIAL_CRT_HPP
#define POLYNOMIAL_CRT_HPP

#include <utility>
#include <vector>

namespace polynomial {

// 中国剰余定理
// リターン値を (r, m) とすると解は x ≡ r (mod. m)
// 解なしの場合は (0, -1) をリターン
std::pair<int64_t, int64_t> crt(int64_t b1, int64_t m1, int64_t b2, int64_t m2);

// 中国剰余定理
// リターン値を (r, m) とすると解は x ≡ r (mod. m)
// 解なしの場合は (0, -1) をリターン
std::pair<int64_t, int64_t> crt(const std::vector<int64_t>& B, const std::vector<int64_t>& M);
}  // namespace polynomial

#endif
