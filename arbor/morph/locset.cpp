#include <algorithm>
#include <iostream>
#include <numeric>

#include <arbor/morph/morphexcept.hpp>
#include <arbor/morph/morphology.hpp>
#include <arbor/morph/locset.hpp>
#include <arbor/morph/mprovider.hpp>

#include "util/rangeutil.hpp"
#include "util/strprintf.hpp"

namespace arb {
namespace ls {

// Throw on invalid mlocation.
void assert_valid(mlocation x) {
    if (!test_invariants(x)) {
        throw invalid_mlocation(x);
    }
}

// Advance an iterator to the first value that is not equal to its current
// value, or end, whichever comes first.
template <typename T>
T next_unique(T& it, T end) {
    const auto& x = *it;
    ++it;
    while (it!=end && *it==x) ++it;
    return it;
};

// Return the number of times that the value at it is repeated. Advances the
// iterator to the first value not equal to its current value, or end,
// whichever comse first.
template <typename T>
int multiplicity(T& it, T end) {
    const auto b = it;
    return std::distance(b, next_unique(it, end));
};

mlocation_list sum(const mlocation_list& lhs, const mlocation_list& rhs) {
    mlocation_list v;
    v.resize(lhs.size() + rhs.size());
    std::merge(lhs.begin(), lhs.end(), rhs.begin(), rhs.end(), v.begin());
    return v;
}

mlocation_list join(const mlocation_list& lhs, const mlocation_list& rhs) {
    mlocation_list L;
    L.reserve(lhs.size()+rhs.size());

    auto l    = lhs.begin();
    auto lend = lhs.end();
    auto r    = rhs.begin();
    auto rend = rhs.end();

    auto at_end = [&]() { return l==lend || r==rend; };
    while (!at_end()) {
        auto x = (*l<*r) ? *l: *r;
        auto count = (*l<*r)? multiplicity(l, lend):
                     (*r<*l)? multiplicity(r, rend):
                     std::max(multiplicity(l, lend), multiplicity(r, rend));
        L.insert(L.end(), count, x);
    }
    L.insert(L.end(), l, lend);
    L.insert(L.end(), r, rend);

    return L;
}

mlocation_list intersection(const mlocation_list& lhs, const mlocation_list& rhs) {
    mlocation_list L;
    L.reserve(lhs.size()+rhs.size());

    auto l    = lhs.begin();
    auto lend = lhs.end();
    auto r    = rhs.begin();
    auto rend = rhs.end();

    auto at_end = [&]() { return l==lend || r==rend; };
    while (!at_end()) {
        if (*l==*r) {
            auto x = *l;
            auto count = std::min(multiplicity(l, lend), multiplicity(r, rend));
            L.insert(L.end(), count, x);
        }
        else if (*l<*r) {
            next_unique(l, lend);
        }
        else {
            next_unique(r, rend);
        }
    }

    return L;
}


// Empty locset.

struct nil_ {};

locset nil() {
    return locset{nil_{}};
}

mlocation_list thingify_(const nil_& x, const mprovider&) {
    return {};
}

std::ostream& operator<<(std::ostream& o, const nil_& x) {
    return o << "nil";
}

// An explicit location.

struct location_ {
    mlocation loc;
};

locset location(mlocation loc) {
    assert_valid(loc);
    return locset{location_{loc}};
}

mlocation_list thingify_(const location_& x, const mprovider& p) {
    assert_valid(x.loc);
    if (x.loc.branch>=p.morphology().num_branches()) {
        throw no_such_branch(x.loc.branch);
    }
    return {x.loc};
}

std::ostream& operator<<(std::ostream& o, const location_& x) {
    return o << "(location " << x.loc.branch << " " << x.loc.pos << ")";
}


// Location corresponding to a sample id.

struct sample_ {
    msize_t index;
};

locset sample(msize_t index) {
    return locset{sample_{index}};
}

mlocation_list thingify_(const sample_& x, const mprovider& p) {
    return {p.embedding().sample_location(x.index)};
}

std::ostream& operator<<(std::ostream& o, const sample_& x) {
    return o << "(sample " << x.index << ")";
}

// Set of terminal points (most distal points).

struct terminal_ {};

locset terminal() {
    return locset{terminal_{}};
}

mlocation_list thingify_(const terminal_&, const mprovider& p) {
    mlocation_list locs;
    util::assign(locs, util::transform_view(p.morphology().terminal_branches(),
        [](msize_t bid) { return mlocation{bid, 1.}; }));

    return locs;
}

std::ostream& operator<<(std::ostream& o, const terminal_& x) {
    return o << "terminal";
}

// Root location (most proximal point).

struct root_ {};

locset root() {
    return locset{root_{}};
}

mlocation_list thingify_(const root_&, const mprovider& p) {
    return {mlocation{0, 0.}};
}

std::ostream& operator<<(std::ostream& o, const root_& x) {
    return o << "root";
}

// Named locset.

struct named_ {
    std::string name;
};

locset named(std::string name) {
    return locset(named_{std::move(name)});
}

mlocation_list thingify_(const named_& n, const mprovider& p) {
    return p.locset(n.name);
}

std::ostream& operator<<(std::ostream& o, const named_& x) {
    return o << "(named \"" << x.name << "\")";
}


// Intersection of two point sets.

struct land {
    locset lhs;
    locset rhs;
    land(locset lhs, locset rhs): lhs(std::move(lhs)), rhs(std::move(rhs)) {}
};

mlocation_list thingify_(const land& P, const mprovider& p) {
    return intersection(thingify(P.lhs, p), thingify(P.rhs, p));
}

std::ostream& operator<<(std::ostream& o, const land& x) {
    return o << "(intersect " << x.lhs << " " << x.rhs << ")";
}

// Union of two point sets.

struct lor {
    locset lhs;
    locset rhs;
    lor(locset lhs, locset rhs): lhs(std::move(lhs)), rhs(std::move(rhs)) {}
};

mlocation_list thingify_(const lor& P, const mprovider& p) {
    return join(thingify(P.lhs, p), thingify(P.rhs, p));
}

std::ostream& operator<<(std::ostream& o, const lor& x) {
    return o << "(join " << x.lhs << " " << x.rhs << ")";
}

// Sum of two point sets.

struct lsum {
    locset lhs;
    locset rhs;
    lsum(locset lhs, locset rhs): lhs(std::move(lhs)), rhs(std::move(rhs)) {}
};

mlocation_list thingify_(const lsum& P, const mprovider& p) {
    return sum(thingify(P.lhs, p), thingify(P.rhs, p));
}

std::ostream& operator<<(std::ostream& o, const lsum& x) {
    return o << "(sum " << x.lhs << " " << x.rhs << ")";
}

} // namespace ls

// The intersect and join operations in the arb:: namespace with locset so that
// ADL allows for construction of expressions with locsets without having
// to namespace qualify the intersect/join.

locset intersect(locset lhs, locset rhs) {
    return locset(ls::land(std::move(lhs), std::move(rhs)));
}

locset join(locset lhs, locset rhs) {
    return locset(ls::lor(std::move(lhs), std::move(rhs)));
}

locset sum(locset lhs, locset rhs) {
    return locset(ls::lsum(std::move(lhs), std::move(rhs)));
}

locset::locset() {
    *this = ls::nil();
}

locset::locset(mlocation loc) {
    *this = ls::location(loc);
}

locset::locset(std::string name) {
    *this = ls::named(std::move(name));
}

} // namespace arb
