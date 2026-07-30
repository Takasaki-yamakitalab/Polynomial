#include "Polynomial/zassenhaus_recombine.hpp"

namespace polynomial {

void enumerateSubset(int pos, int need, std::vector<int>& cur, std::vector<std::vector<int>>& ans, int n) {
    if ((int)cur.size() == need) {
        ans.push_back(cur);
        return;
    }
    for (int i = pos; i < n; i++) {
        cur.push_back(i);
        enumerateSubset(i + 1, need, cur, ans, n);
        cur.pop_back();
    }
}

}  // namespace polynomial
