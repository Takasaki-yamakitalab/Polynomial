#ifndef POLYNOMIAL_ENUMERATE_SUBSET_HPP
#define POLYNOMIAL_ENUMERATE_SUBSET_HPP

#include <vector>

namespace polynomial {

void enumerateSubsetImpl(int pos, int need, std::vector<int>& cur, std::vector<std::vector<int>>& ans, int n,
                         const std::vector<bool>& exclude);
void enumerateSubset(std::vector<std::vector<int>>& dst, const std::vector<bool>& exclude, int need);

}
#endif
