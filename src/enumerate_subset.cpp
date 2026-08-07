#include "Polynomial/debug_messager.hpp"

namespace polynomial {

void enumerateSubsetImpl(int pos, int need, std::vector<int>& cur, std::vector<std::vector<int>>& dst, int n, const std::vector<bool>& exclude) {
    if ((int)cur.size() == need) {
        dst.push_back(cur);
        return;
    }
    for (int i = pos; i < exclude.size(); i++) {
        if (exclude[i]) continue;
        cur.push_back(i);
        enumerateSubsetImpl(i + 1, need, cur, dst, n, exclude);
        cur.pop_back();
    }
}
void enumerateSubset(std::vector<std::vector<int>>& dst, const std::vector<bool>& exclude, int need) {
    dst.clear();
    int available = 0;
    for (size_t i = 0; i < exclude.size(); i++) {
        if (!exclude[i]) available++;
    }
    if (need >= available) {
        dst.emplace_back();
        for (int i = 0; i < (int)exclude.size(); i++) {
            if (!exclude[i]) dst.back().push_back(i);
        }
        return;
    }
    DEBUG_OUT("need: ", need);
    DEBUG_OUT("available: ", available);
    std::vector<int> cur;
    enumerateSubsetImpl(0, need, cur, dst, available, exclude);
    return;
}

}  // namespace polynomial
