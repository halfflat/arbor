#pragma once

template <typename I>
struct count_range {
    struct proxy {
        I counter;
        proxy(I x): counter(x) {}
        I operator*() const { return counter; }
        proxy& operator++() { ++counter; return *this; }
        bool operator==(const proxy& p) const { return counter==p.counter; }
        bool operator!=(const proxy& p) const { return counter!=p.counter; }
    };

    proxy begin() const { return proxy(i0); }
    proxy end() const { return proxy(i1); }

    count_range() = default;
    explicit count_range(I n): i0(I{}), i1(n) {}
    explicit count_range(I i0, I i1): i0(i0), i1(i1) {}

    I i0, i1;
};

template <typename C>
auto size(const C& c) { return c.size(); }

template <typename X, std::size_t n>
auto size(const X (&c)[n]) { return n; }

template <typename I>
auto count(I n) { return count_range<I>(n); }

template <typename C>
auto count_along(const C& c) { return count(size(c)); }

template <typename I, typename J>
auto span(I i0, J i1) { return count_range<decltype(true? i0: i1)>(i0, i1); }
