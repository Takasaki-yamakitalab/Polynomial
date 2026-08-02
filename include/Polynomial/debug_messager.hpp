#ifndef POLYNOMIAL_DEBUG_MESSAGER_HPP
#define POLYNOMIAL_DEBUG_MESSAGER_HPP

#include <vector>
#include <iostream>

#ifdef _DEBUG
// #define DEBUG_OUT(...)                     \
//     do {                                   \
//         std::cout << "[DEBUG] ";           \
//         ((std::cout << __VA_ARGS__), ...); \
//         std::cout << std::endl;            \
//     } while (0)
#define DEBUG_OUT(...)                                         \
    do {                                                       \
        debug_messager::debug_print_helper(__VA_ARGS__); \
    } while (0)
#else
#define DEBUG_OUT(...)
#endif

template <typename T>
std::ostream& operator<<(std::ostream& os, const std::vector<T>& vec) {
    os << "[";
    for (size_t i = 0; i < vec.size(); i++) {
        os << vec[i];
        if (i + 1 < vec.size()) {
            os << ", ";
        }
    }
    os << "]";
    return os;
}

namespace debug_messager {
template <typename... Args>
void debug_print_helper(Args&&... args);

template <typename... Args>
void debug_print_helper(Args&&... args) {
    std::cout << "[DEBUG] ";
    ((std::cout << std::forward<Args>(args)), ...);
    std::cout << std::endl;
}
}  // namespace debug_messager

#endif
