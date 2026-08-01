#ifndef POLYNOMIAL_MODINT_HPP
#define POLYNOMIAL_MODINT_HPP

#include <ostream>
#include <cassert>
#include "coefficient_traits.hpp"

namespace polynomial {

struct ModContext {
    inline static int P;

    static int modulus() {
        return P;
    };
    static void setMod(int p) {
        assert(p > 1);
        P = p;
    }
};

// P must be prime
class ModInt {
    int v;
    int P;

   public:
    void setMod(int p) {
        assert(p > 1);
        P = p;
        v = v % P;
        if (v < 0) v += P;
    }
    ModInt() : v(0), P(ModContext::modulus()) {}
    template <class Integer>
    ModInt(Integer value) : P(ModContext::modulus()) {
        int64_t x = value % P;
        if (x < 0) x += P;
        v = (int)x;
    }

    int value() const { return v; }
    int modulus() const { return P; }
    bool isZero() const { return v == 0; }
    bool isOne() const { return v == 1; }
    static int modPow(int a, int64_t e, int p) {
        int64_t r = 1;
        int64_t x = a;
        while (e > 0) {
            if (e & 1) r = r * x % p;
            x = x * x % p;
            e >>= 1;
        }
        return (int)r;
    }
    ModInt inverse() const {
        assert(v != 0);
        return ModInt(modPow(v, P - 2, P));
    }

    ModInt operator-() const { return ModInt(-v); }

    ModInt operator+(const ModInt& b) const { assert(P == b.P); return ModInt(v + b.v); }
    ModInt operator-(const ModInt& b) const { assert(P == b.P); return ModInt(v - b.v); }
    ModInt operator*(const ModInt& b) const { assert(P == b.P); return ModInt((int64_t)v * b.v); }
    ModInt operator/(const ModInt& b) const { assert(P == b.P); return *this * b.inverse(); }

    friend ModInt operator+(const int& a, const ModInt& b) { return ModInt(a) + b; }
    friend ModInt operator-(const int& a, const ModInt& b) { return ModInt(a) - b; }
    friend ModInt operator*(const int& a, const ModInt& b) { return ModInt(a) * b; }
    friend ModInt operator/(const int& a, const ModInt& b) { return ModInt(a) / b; }

    ModInt& operator+=(const ModInt& b) {
        *this = *this + b;
        return *this;
    }
    ModInt& operator-=(const ModInt& b) {
        *this = *this - b;
        return *this;
    }
    ModInt& operator*=(const ModInt& b) {
        *this = *this * b;
        return *this;
    }
    ModInt& operator/=(const ModInt& b) {
        *this = *this / b;
        return *this;
    }

    bool operator==(const ModInt& b) const { return P == b.P && v == b.v; }
    bool operator!=(const ModInt& b) const { return P != b.P || v != b.v; }
    bool operator<(const ModInt& b) const { return P == b.P && v < b.v; }
    bool operator>(const ModInt& b) const { return P == b.P && v > b.v; }
    bool operator<=(const ModInt& b) const { return P == b.P && v <= b.v; }
    bool operator>=(const ModInt& b) const { return P == b.P && v >= b.v; }

    friend bool operator==(const ModInt& a, const int& b) { return a == ModInt(b); }
    friend bool operator==(const int& a, const ModInt& b) { return ModInt(a) == b; }
    friend bool operator!=(const ModInt& a, const int& b) { return a != ModInt(b); }
    friend bool operator!=(const int& a, const ModInt& b) { return ModInt(a) != b; }

    friend std::ostream& operator<<(std::ostream& os, const ModInt& a) { return os << a.v; }
};

template <>
struct coefficient_traits<ModInt> {
    static const bool is_field = true;
};

}  // namespace polynomial

#endif
