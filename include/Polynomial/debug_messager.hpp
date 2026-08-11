#ifndef POLYNOMIAL_DEBUG_MESSAGER_HPP
#define POLYNOMIAL_DEBUG_MESSAGER_HPP

#include <vector>
#include <iostream>

#ifdef ANDROID
#include <android/log.h>
#include <sstream>
#endif

#define DEBUG_OUT(...)                                   \
    do {                                                 \
        debug_messager::debug_print_helper(__VA_ARGS__); \
    } while (0)

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

#ifdef ANDROID

template <typename... Args>
void debug_print_helper(Args&&... args) {
#ifdef _DEBUG
    std::ostringstream oss;
    oss << "[DEBUG] ";
    ((oss << std::forward<Args>(args)), ...);
    oss << std::endl;
    __android_log_print(ANDROID_LOG_INFO, "DEBUG_OUT", "%s", oss.str().c_str());
#else
    (((void)std::forward<Args>(args)), ...);
#endif
}

#else

template <typename... Args>
void debug_print_helper(Args&&... args) {
#ifdef _DEBUG
    std::cout << "[DEBUG] ";
    ((std::cout << std::forward<Args>(args)), ...);
    std::cout << std::endl;
#else
    (((void)std::forward<Args>(args)), ...);
#endif
}

#endif

}  // namespace debug_messager

#endif
