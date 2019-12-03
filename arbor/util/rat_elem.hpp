#pragma once

// An element representing a segment of a rational polynomial function
// of order p, q, as determined by its values on n = p+q+1 nodes at
// [0, 1/n, ..., 1].

#include <array>
#include <functional>
#include <type_traits>
#include <utility>

#include <arbor/math.hpp>

namespace arb {

namespace impl {

template <unsigned n, unsigned sz>
struct array_init_n {
    template <typename A, typename X, typename... Tail>
    static void set(A& array, X value, Tail... tail) {
        array[sz-n] = std::move(value);
        init_n<n-1, sz>::set(array, std::forward<Tail>(tail)...);
    }
};

template <unsigned sz>
struct array_init_n<0, sz> {
    template <typename A>
    static void set(A& array) {}
};

template <unsigned p, unsigned q>
struct rat_eval {
    // (Linear only for now; TODO add generic).
    static double apply(const std::array<double, 1+p+q>& data, double x) {
        if (p!=1 || q!=0) throw std::logic_error("woah, not implemented!");
        return lerp(data[0], data[p+q], x);
    }
};

} // namespace impl

template <unsigned p, unsigned q>
struct rat_element {
    template <typename F, typename _ = std::enable_if_t<!std::is_arithmetic<std::decay_t<F>>::value>>
    rat_element(F&& fn, _* = nullptr) {
        for (unsigned i = 0; i<size(); ++i) data[i] = fn(i/(size()-1.0));
    }

    template <typename... Tail>
    rat_element(double y0, Tail... tail) {
        impl::array_init_n<p+q+1, p+q+1>::set(data, y0, tail...);
    }

    constexpr unsigned size() const { return 1+p+q; }

    double operator()(double x) const {
        return impl::rat_eval<p, q>::apply(data, x);
    }

    double operator[](unsigned i) const {
        return data.at(i);
    }

    double& operator[](unsigned i) {
        return data.at(i);
    }

private:
    std::array<double, 1+p+q> data;
};

} // namespace arb
