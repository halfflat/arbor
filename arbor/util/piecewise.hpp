#pragma once

// Create/manipulate 1-d piece-wise defined objects.
//
// Using vectors everywhere here for ease; consider making
// something more container/sequence-generic later.

#include <initializer_list>
#include <vector>

#include "util/partition.hpp"

namespace arb {

using pw_size_type = unsigned;
constexpr pw_size_type pw_npos = -1;

template <typename X = void>
struct pw_elements {
    using size_type = pw_size_type;
    static constexpr size_type npos = pw_npos;

    // Consistency requirements:
    // 1. empty() || element.size()+1 = vertex.size()
    // 2. vertex[i]<=vertex[j] for i<=j.

    std::vector<double> vertex_;
    std::vector<X> element_;

    // ctors and assignment:

    pw_elements() = default;

    template <typename VSeq, typename ESeq>
    pw_elements(const VSeq& vs, const ESeq& es) {
        assign(vs, es);
    }

    pw_elements(std::initializer_list<double> vs, std::initializer_list<X> es) {
        assign(vs, es);
    }

    pw_elements(pw_elements&&) = default;
    pw_elements(const pw_elements&) = default;

    template <typename Y>
    explicit pw_elements(const pw_elements<Y>& from):
        vertex_(from.vertex_), element_(from.element_.begin(), from.element_.end())
    {}

    pw_elements& operator=(pw_elements&&) = default;
    pw_elements& operator=(const pw_elements&) = default;

    // access:

    auto intervals() const { return util::partition_view(vertex_); }
    auto interval(size_type i) const { return intervals()[i]; }

    auto bounds() const { return intervals().bounds(); }

    size_type size() const { return element_.size(); }
    bool empty() const { return size()==0; }

    bool operator==(const pw_elements& x) const {
        return vertex_==x.vertex_ && element_==x.element_;
    }

    bool operator!=(const pw_elements& x) const { return !(*this==x); }

    const auto& elements() const { return element_; }
    const auto& vertices() const { return vertex_; }

    X& element(size_type i) & { return element_[i]; }
    const X& element(size_type i) const & { return element_[i]; }

    // TODO: change this (and add const iterators) to return by
    // value a pair of element and bounds?

    X& operator[](size_type i) { return element(i); }
    const X& operator[](size_type i) const { return element(i); }

    size_type index_of(double x) const {
        if (empty()) return npos;

        auto partn = intervals();
        if (x == partn.bounds().second) return size()-1;
        else return partn.index(x);
    }

    // mutating operations:

    void reserve(size_type n) {
        vertex_.reserve(n+1);
        element_.reserve(n);
    }

    void clear() {
        vertex_.clear();
        element_.clear();
    }

    template <typename U>
    void push_back(double left, double right, U&& elem) {
        if (!empty() && left!=vertex_.back()) {
            throw std::runtime_error("noncontiguous element");
        }

        if (right<left) {
            throw std::runtime_error("inverted element");
        }

        // Extend element_ first in case a conversion/copy/move throws.

        element_.push_back(std::forward<U>(elem));
        if (vertex_.empty()) vertex_.push_back(left);
        vertex_.push_back(right);
    }

    template <typename U>
    void push_back(double right, U&& elem) {
        if (empty()) {
            throw std::runtime_error("require initial left vertex for element");
        }

        push_back(vertex_.back(), right, elem);
    }

    void assign(std::initializer_list<double> vs, std::initializer_list<X> es) {
        using util::make_range;
        assign(make_range(vs.begin(), vs.end()), make_range(es.begin(), es.end()));
    }

    template <typename Seq1, typename Seq2>
    void assign(const Seq1& vertices, const Seq2& elements) {
        using std::begin;
        using std::end;

        auto vi = begin(vertices);
        auto ve = end(vertices);

        auto ei = begin(elements);
        auto ee = end(elements);

        if (ei==ee) { // empty case
            if (vi!=ve) {
                throw std::runtime_error("vertex list too long");
            }
            clear();
            return;
        }

        double left = *vi++;
        if (vi==ve) {
            throw std::runtime_error("vertex list too short");
        }

        clear();

        double right = *vi++;
        push_back(left, right, *ei++);

        while (ei!=ee) {
            if (vi==ve) {
                throw std::runtime_error("vertex list too short");
            }
            double right = *vi++;
            push_back(right, *ei++);
        }

        if (vi!=ve) {
            throw std::runtime_error("vertex list too long");
        }
    }
};

// With X = void, present the element intervals only,
// keeping othewise the same interface.

template <> struct pw_elements<void> {
    using size_type = pw_size_type;
    static constexpr size_type npos = pw_npos;

    std::vector<double> vertex_;

    // ctors and assignment:

    template <typename VSeq>
    explicit pw_elements(const VSeq& vs) { assign(vs); }

    pw_elements(std::initializer_list<double> vs) { assign(vs); }

    pw_elements() = default;
    pw_elements(pw_elements&&) = default;
    pw_elements(const pw_elements&) = default;

    template <typename Y>
    explicit pw_elements(const pw_elements<Y>& from):
        vertex_(from.vertex_) {}

    pw_elements& operator=(pw_elements&&) = default;
    pw_elements& operator=(const pw_elements&) = default;

    // access:

    auto intervals() const { return util::partition_view(vertex_); }
    auto interval(size_type i) const { return intervals()[i]; }

    auto bounds() const { return intervals().bounds(); }

    size_type size() const { return vertex_.empty()? 0: vertex_.size()-1; }
    bool empty() const { return vertex_.empty(); }

    bool operator==(const pw_elements& x) const { return vertex_==x.vertex_; }
    bool operator!=(const pw_elements& x) const { return !(*this==x); }

    const auto& vertices() const { return vertex_; }

    size_type index_of(double x) const {
        if (empty()) return npos;

        auto partn = intervals();
        if (x == partn.bounds().second) return size()-1;
        else return partn.index(x);
    }

    // mutating operations:

    void reserve(size_type n) { vertex_.reserve(n+1); }
    void clear() { vertex_.clear(); }

    void push_back(double left, double right) {
        if (!empty() && left!=vertex_.back()) {
            throw std::runtime_error("noncontiguous element");
        }

        if (right<left) {
            throw std::runtime_error("inverted element");
        }

        if (vertex_.empty()) vertex_.push_back(left);
        vertex_.push_back(right);
    }

    void push_back(double right) {
        if (empty()) {
            throw std::runtime_error("require initial left vertex for element");
        }
        vertex_.push_back(right);
    }

    void assign(std::initializer_list<double> vs) {
        assign(util::make_range(vs.begin(), vs.end()));
    }

    template <typename Seq1>
    void assign(const Seq1& vertices) {
        using std::begin;
        using std::end;

        auto vi = begin(vertices);
        auto ve = end(vertices);

        if (vi==ve) {
            clear();
            return;
        }

        double left = *vi++;
        if (vi==ve) {
            throw std::runtime_error("vertex list too short");
        }

        clear();

        double right = *vi++;
        push_back(left, right);

        while (vi!=ve) {
            double right = *vi++;
            push_back(right);
        }
    }
};


namespace impl {
    template <typename A, typename B> 
    struct pair_type_map { using type = std::pair<A, B>; };

    template <typename A>
    struct pair_type_map<A, void> { using type = A; };

    template <typename B>
    struct pair_type_map<void, B> { using type = B; };

    template <>
    struct pair_type_map<void, void> { using type = void; };

    template <typename A, typename B>
    using pair_type = typename pair_type_map<A, B>::type;

    template <typename A, typename B>
    void general_pw_push_pair(
        pw_elements<std::pair<A, B>>& out,
        double left, double right,
        const pw_elements<A>& a, pw_size_type ai,
        const pw_elements<B>& b, pw_size_type bi)
    {
        out.push_back(left, right, std::pair<A, B>{a[ai], b[bi]});
    }

    template <typename A>
    void general_pw_push_pair(
        pw_elements<A>& out,
        double left, double right,
        const pw_elements<A>& a, pw_size_type ai,
        const pw_elements<void>& b, pw_size_type bi)
    {
        out.push_back(left, right, a[ai]);
    }

    template <typename B>
    void general_pw_push_pair(
        pw_elements<B>& out,
        double left, double right,
        const pw_elements<void>& a, pw_size_type ai,
        const pw_elements<B>& b, pw_size_type bi)
    {
        out.push_back(left, right, b[bi]);
    }

    void general_pw_push_pair(
        pw_elements<void>& out,
        double left, double right,
        const pw_elements<void>& a, pw_size_type ai,
        const pw_elements<void>& b, pw_size_type bi)
    {
        out.push_back(left, right);
    }
}

// TODO: Consider making a lazy `meet_view` version of meet.

template <typename A, typename B>
pw_elements<impl::pair_type<A, B>> meet(const pw_elements<A>& a, const pw_elements<B>& b) {
    pw_elements<impl::pair_type<A, B>> m;

    if (a.empty() || b.empty()) return m;

    double lmax = std::max(a.bounds().first, b.bounds().first);
    double rmin = std::min(a.bounds().second, b.bounds().second);
    if (rmin<lmax) return m;

    double left = lmax;
    pw_size_type ai = a.intervals().index(left);
    pw_size_type bi = b.intervals().index(left);

    if (rmin==left) {
        impl::general_pw_push_pair(m, left, left, a, ai, b, bi);
        return m;
    }

    double a_right = a.interval(ai).second;
    double b_right = b.interval(bi).second;

    while (left<rmin) {
        double right = std::min(a_right, b_right);
        right = std::min(right, rmin);

        impl::general_pw_push_pair(m, left, right, a, ai, b, bi);
        if (a_right<=right) {
            a_right = a.interval(++ai).second;
        }
        if (b_right<=right) {
            b_right = b.interval(++bi).second;
        }

        left = right;
    }
    return m;
}

} // namespace arb
