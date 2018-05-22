#pragma once

#include <iostream>
#include <random>
#include <vector>

// Represent vector of ordered numeric values as a size and range.
// Where relevant, values are presumed to be equally spaced across range.

template <typename X>
struct mock_vector {
    std::size_t n_;
    X min_, max_;

    mock_vector() { clear(); }
    mock_vector(const mock_vector&) = default;
    mock_vector& operator=(const mock_vector&) = default;

    void clear() { n_ = 0; min_ = X{}; max_ = X{}; }
    bool empty() const { return n_==0; }
    std::size_t size() const { return n_; }
    void resize(std::size_t n) { n_ = n; }

    void push_back(X v) {
        min_ = empty()? v: std::min(v, min_);
        max_ = empty()? v: std::max(v, max_);
        ++n_;
    }

    void insert(std::size_t nelem) {
        n_ += nelem;
    }

    X front() const { return min_; }
    X back()  const { return max_; }

    void append(const mock_vector& m) {
        if (m.empty()) return;
        min_ = empty()? m.min_: std::min(min_, m.min_);
        max_ = empty()? m.max_: std::max(max_, m.max_);
        n_ += m.n_;
    }

    mock_vector take_upto(X v) {
        mock_vector taken;

        if (!empty() && v>=min_) {
            if (v>=max_ || n_==1 || min_==max_) {
                taken = *this;
                n_ = 0;
                min_ = max_;
            }
            else {
                double delta = double(max_-min_)/(n_-1);
                taken.min_ = min_;
                taken.n_ = std::min(n_, std::size_t((v-min_)/delta));
                taken.max_ = std::min(max_, X(min_ + taken.n_*delta));

                n_ -= taken.n_;
                min_ = std::min(max_, X(min_ + (taken.n_+1)*delta));
            }
        }
        return taken;
    }

    friend std::ostream& operator<<(std::ostream& out, const mock_vector& m) {
        if (m.empty()) return out << "[]";
        else return out << '[' << m.n_ << "; " << m.min_ << " -- " << m.max_ << ']';
    }
};

// Make a mock vector with values in range [x0, x1) by a Poisson
// sampling process with given rate.

template <typename X, typename Rng>
mock_vector<X> mock_poisson(X x0, X x1, double rate, Rng& rng) {
    mock_vector<X> v;

    std::exponential_distribution<double> E(rate);
    std::poisson_distribution<std::ptrdiff_t> P;
    using poisson_rate = decltype(P)::param_type;

    double pmin = double(x0) + E(rng);
    double pmax = double(x1) - E(rng);
    std::ptrdiff_t n = 0;

    X xmin(pmin);
    X xmax(pmax);

    if (xmin<x1) {
        v.push_back(xmin);
        if (xmax>x0) {
            v.push_back(xmax);
            v.insert(P(rng, poisson_rate{rate*double(pmax-pmin)}));
        }
    }

    return v;
}

