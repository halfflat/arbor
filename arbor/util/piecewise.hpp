#pragma once

// Create/manipulate 1-d piece-wise defined objects.
//
// Using vectors everywhere here for ease; consider making
// something more container/sequence-generic later.

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

    auto intervals() const { return util::partition_view(vertex_); }
    auto interval(size_type i) const { return intervals()[i]; }

    auto bounds() const { return intervals().bounds(); }

    size_type size() const { return element_.size(); }
    bool empty() const { return size()==0; }

    void reserve(size_type n) {
        vertex_.reserve(n+1);
        element_.reserve(n);
    }

    void clear() {
        vertex_.clear();
        element_.clear();
    }

    const auto& elements() const { return element_; }
    const auto& vertices() const { return vertex_; }

    X& element(size_type i) & { return element_[i]; }
    const X& element(size_type i) const & { return element_[i]; }

    X& operator[](size_type i) { return element(i); }
    const X& operator[](size_type i) const { return element(i); }

    size_type index_of(double x) {
        if (empty()) return npos;

        auto partn = intervals();
        if (x == partn.bounds().second) return size()-1;
        else return partn.index(x);
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

    template <typename Seq1, typename Seq2>
    void assign(const Seq1& vertices, const Seq2& elements) {
        using std::begin;
        using std::end;

        auto vi = begin(vertices);
        auto ve = end(vertices);

        auto ei = begin(elements);
        auto ee = end(elements);

        if (ei==ee) { // empty case
            if (vi==ve) {
                throw std::runtime_error("vertex list too long");
            }
            clear();
            return;
        }

        double left = *vi++;
        if (vi==ve) {
            throw std::runtime_error("vertex list too short");
        }

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

    auto intervals() const { return util::partition_view(vertex_); }
    auto interval(size_type i) const { return intervals()[i]; }

    auto bounds() const { return intervals().bounds(); }

    size_type size() const { return vertex_.empty()? 0: vertex_.size()-1; }
    bool empty() const { return vertex_.empty(); }

    const auto& vertices() const { return vertex_; }

    size_type index_of(double x) {
        if (empty()) return npos;

        auto partn = intervals();
        if (x == partn.bounds().second) return size()-1;
        else return partn.index(x);
    }

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



#if 0

// View of a restricted or subdivided pw_elements object.
// Overrides element()/operator[] to reference elements in the viewed object.
// Provides index() for access to to the indices to the viewed object.

template <typename X>
struct pw_subview: pw_elements<pw_size_type> {
    using base = pw_elements<pw_size_type>;
    using base::size_type;

    const pw_elements<X>& data;

    size_type index(size_type i) const { return element_[i]; }

    X& element(size_type i) { return data.element(index(i)); }
    const X& element(size_type i) const { return data.element(index(i)); }

    X& operator[](size_type i) { return element(index(i)); }
    const X& operator[](size_type i) const { return element(index(i)); }

    explicit pw_subview(const pw_elements<X>& els): data(els) {}

    void push_back(double left, double right) {
        size_type ridx = data.index_of(right);
        if (ridx == npos) {
            throw std::range_error("span outside support");
        }

        size_type lidx = empty()? data.index_of(left): element_.back();
        if (lidx == npos) {
            throw std::range_error("span outside support");
        }
        else if (lidx != ridx) {
            throw std::range_error("span crossed multiple elements");
        }

        base::push_back(left, right, ridx);
    }

    void push_back(double right) {
        if (empty()) {
            throw std::runtime_error("require initial left vertex for element");
        }

        push_back(vertex_.back(), right);
    }

    template <typename Seq1>
    void assign(const Seq1& vertices) {
        impl::generic_vertices_assign(*this, vertices);
    }
};

template <typename X>
pw_subview<X> restrict(const pw_elements<X>& pw, double left, double right) {
    constexpr pw_size_type npos = pw_npos;

    pw_subview<X> v(pw);
    if (pw.empty() || right<left) return v;

    left = std::max(left, pw.bounds().first);
    right = std::min(right, pw.bounds().second);

    pw_size_type lidx = pw.intervals().index(left);
    pw_size_type ridx = pw.intervals().index(right);
    if (ridx==npos) ridx = pw.size()-1;

    if (lidx==ridx) {
        v.push_back(left, right);
        return v;
    }

    v.push_back(left, pw.interval(lidx).second);
    while (lidx<ridx) {
        v.push_back(pw.interval(++lidx).second);
    }
    v.push_back(right);
    return v;
}

// Create a view of `pw` with the same vertices in addition to
// any vertices from `by` that lie within the bounds of `pw`.

template <typename X, typename Seq>
pw_subview<X> subdivide(const pw_elements<X>& pw, const Seq& by) {
    using std::begin;
    using std::end;

    auto bi = begin(by);
    auto be = end(by);

    double left = pw.bounds().first;
    while (bi!=be && *bi<=left) ++bi;

    // TODO ...
}
#endif

} // namespace arb
