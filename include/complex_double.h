#ifndef COMPLEX_DOUBLE_H
#define COMPLEX_DOUBLE_H

#include <cmath>
#include <ostream>
#include "backend.h"

/*
    Double precision complex number that can be used both on the CPU and inside GPU kernels
    (std::complex cannot be used in CUDA device code).
    It has the same memory layout as std::complex<double>, so it is communicated with MPI_DOUBLE_COMPLEX.
    The usual functions are available: real, imag, conj, norm (= |z|^2, as std::norm), abs, arg and exp.

    Arguments are passed by value on purpose: this allows to use constexpr constants, such as I_number,
    inside GPU kernels.
*/
struct alignas(16) c_double {
    double re;
    double im;

    HOST_DEVICE constexpr c_double(double r = 0.0, double i = 0.0) : re(r), im(i) {}

    HOST_DEVICE constexpr double real() const { return re; }
    HOST_DEVICE constexpr double imag() const { return im; }

    HOST_DEVICE constexpr c_double& operator+=(c_double z) { re += z.re; im += z.im; return *this; }
    HOST_DEVICE constexpr c_double& operator-=(c_double z) { re -= z.re; im -= z.im; return *this; }
    HOST_DEVICE constexpr c_double& operator*=(double a) { re *= a; im *= a; return *this; }
    HOST_DEVICE constexpr c_double& operator*=(c_double z) {
        const double r = re * z.re - im * z.im;
        im = re * z.im + im * z.re;
        re = r;
        return *this;
    }
};

HOST_DEVICE constexpr c_double operator+(c_double a) { return a; }
HOST_DEVICE constexpr c_double operator-(c_double a) { return c_double(-a.re, -a.im); }

HOST_DEVICE constexpr c_double operator+(c_double a, c_double b) { return c_double(a.re + b.re, a.im + b.im); }
HOST_DEVICE constexpr c_double operator+(c_double a, double b)   { return c_double(a.re + b, a.im); }
HOST_DEVICE constexpr c_double operator+(double a, c_double b)   { return c_double(a + b.re, b.im); }

HOST_DEVICE constexpr c_double operator-(c_double a, c_double b) { return c_double(a.re - b.re, a.im - b.im); }
HOST_DEVICE constexpr c_double operator-(c_double a, double b)   { return c_double(a.re - b, a.im); }
HOST_DEVICE constexpr c_double operator-(double a, c_double b)   { return c_double(a - b.re, -b.im); }

HOST_DEVICE constexpr c_double operator*(c_double a, c_double b) {
    return c_double(a.re * b.re - a.im * b.im, a.re * b.im + a.im * b.re);
}
HOST_DEVICE constexpr c_double operator*(c_double a, double b) { return c_double(a.re * b, a.im * b); }
HOST_DEVICE constexpr c_double operator*(double a, c_double b) { return c_double(a * b.re, a * b.im); }

HOST_DEVICE constexpr c_double operator/(c_double a, double b) { return c_double(a.re / b, a.im / b); }
HOST_DEVICE constexpr c_double operator/(c_double a, c_double b) {
    const double den = b.re * b.re + b.im * b.im;
    return c_double((a.re * b.re + a.im * b.im) / den, (a.im * b.re - a.re * b.im) / den);
}
HOST_DEVICE constexpr c_double operator/(double a, c_double b) { return c_double(a) / b; }

HOST_DEVICE constexpr bool operator==(c_double a, c_double b) { return a.re == b.re && a.im == b.im; }
HOST_DEVICE constexpr bool operator!=(c_double a, c_double b) { return !(a == b); }

HOST_DEVICE constexpr double real(c_double z) { return z.re; }
HOST_DEVICE constexpr double imag(c_double z) { return z.im; }
HOST_DEVICE constexpr c_double conj(c_double z) { return c_double(z.re, -z.im); }
HOST_DEVICE constexpr double norm(c_double z) { return z.re * z.re + z.im * z.im; } //|z|^2
HOST_DEVICE inline double abs(c_double z) { return hypot(z.re, z.im); }
HOST_DEVICE inline double arg(c_double z) { return atan2(z.im, z.re); }
HOST_DEVICE inline c_double exp(c_double z) {
    const double modulus = exp(z.re);
    return c_double(modulus * cos(z.im), modulus * sin(z.im));
}

inline std::ostream& operator<<(std::ostream& os, c_double z) {
    return os << '(' << z.re << ',' << z.im << ')';
}

constexpr c_double I_number(0.0, 1.0); //imaginary unit

#endif
