#include "Polynomial/factor.hpp"

namespace polynomial {
std::vector<std::vector<ModInt>> nullSpace(std::vector<std::vector<ModInt>> A) {
    // mod p Gauss elimination
    size_t m = A.size();
    size_t n = A[0].size();
    std::vector<int> pivot(n, -1);

    int row = 0;
    for (int col = 0; col < n && row < m; col++) {
        int sel = row;
        while (sel < m && A[sel][col] == 0) sel++;
        if (sel == m) continue;

        std::swap(A[sel], A[row]);

        ModInt inv = A[row][col].inverse();
        for (int j = 0; j < n; j++) A[row][j] = A[row][j] * inv;

        for (int i = 0; i < m; i++) {
            if (i == row) continue;
            ModInt c = A[i][col];
            for (int j = 0; j < n; j++) {
                A[i][j] = (A[i][j] - c * A[row][j]);
            }
        }
        pivot[col] = row;
        row++;
    }

    std::vector<std::vector<ModInt>> basis;
    for (int freeCol = 0; freeCol < n; freeCol++) {
        if (pivot[freeCol] != -1) continue;
        std::vector<ModInt> v(n);
        v[freeCol] = 1;
        for (int col = 0; col < n; col++) {
            if (pivot[col] != -1) {
                v[col] = -A[pivot[col]][freeCol];
            }
        }
        basis.push_back(v);
    }

    // DEBUG
    for (auto& v : basis) {
        std::cout << "[DEBUG] basis:";
        for (auto& x : v) std::cout << " " << x;
        std::cout << '\n';
    }

    return basis;
}
}  // namespace polynomial
