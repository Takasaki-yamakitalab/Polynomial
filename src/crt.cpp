#include "Polynomial/crt.hpp"

namespace polynomial {

inline int64_t mod(int64_t a, int64_t m) { return (a % m + m) % m; }

// 拡張 Euclid の互除法
// ap + bq = gcd(a, b) となる (p, q) を求め、d = gcd(a, b) を返す
int64_t extendedEuclid(int64_t a, int64_t b, int64_t& p, int64_t& q) {
    if (b == 0) {
        p = 1;
        q = 0;
        return a;
    }
    int64_t d = extendedEuclid(b, a % b, q, p);
    q -= a / b * p;
    return d;
}

// 中国剰余定理
// リターン値を (r, m) とすると解は x ≡ r (mod. m)
// 解なしの場合は (0, -1) をリターン
std::pair<int64_t, int64_t> crt(int64_t b1, int64_t m1, int64_t b2, int64_t m2) {
    int64_t p, q;
    int64_t d = extendedEuclid(m1, m2, p, q);  // p is inv of m1/d (mod. m2/d)
    if ((b2 - b1) % d != 0) return std::make_pair(0, -1);
    int64_t m = m1 * (m2 / d);  // lcm of (m1, m2)
    int64_t tmp = (b2 - b1) / d * p % (m2 / d);
    int64_t r = mod(b1 + m1 * tmp, m);
    return std::make_pair(r, m);
}

// 中国剰余定理
// リターン値を (r, m) とすると解は x ≡ r (mod. m)
// 解なしの場合は (0, -1) をリターン
std::pair<int64_t, int64_t> crt(const std::vector<int64_t>& B, const std::vector<int64_t>& M) {
    int64_t r = 0, m = 1;
    for (size_t i = 0; i < B.size(); i++) {
        int64_t p, q;
        int64_t d = extendedEuclid(m, M[i], p, q);  // p is inv of M/d (mod. m[i]/d)
        if ((B[i] - r) % d != 0) return std::make_pair(0, -1);
        int64_t tmp = (B[i] - r) / d * p % (M[i] / d);
        r += m * tmp;
        m *= M[i] / d;
    }
    return std::make_pair(mod(r, m), m);
}
}  // namespace polynomial
