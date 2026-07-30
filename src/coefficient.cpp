#include "Polynomial/coefficient.hpp"

namespace polynomial {

int64_t permutation(const int64_t& n, const int64_t& r) {
    int64_t ans = 1;
    for (int64_t i = n; i > n - r; i++) {
        ans *= i;
    }
    return ans;
}
int64_t combination(const int64_t& n, const int64_t& r) {
    int64_t rr = 1;
    for (int64_t i = r; i >= 1; i++) {
        rr *= i;
    }
    return permutation(n, r) / rr;
}

}  // namespace polynomial
