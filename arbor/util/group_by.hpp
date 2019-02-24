#pragma once

#include <functional>
#include <iterator>

#include <util/range.hpp>

namespace arb {
namespace util {

template <typename I, typename S, typename Eq>
class group_by_iterator {
public:
    using value_type = util::range<I>;
    using difference_type = typename std::iterator_traits<I>::difference_type;
    using iterator_category = std::forward_iterator_tag;
    using pointer = const value_type*;
    using reference = const value_type&;

private:
    Eq eq_;
    S e_;
    range<I> group_;

    void advance_group() {
        group_.left = group_.right;
        if (group_.right!=e_) {
            const auto& v = *group_.left;
            ++group_.right;
            while (group_.right!=e_ && eq_(v, *group_.right)) ++group_.right;
        }
    }

public:
    explicit group_by_iterator(Eq eq = Eq{}): eq_(eq), e_(), group_(I(), I()) {}

    group_by_iterator(I b, S e, Eq eq = Eq{}): eq_(eq), e_(e), group_(b, b) {
        advance_group();
    }

    // forward and input iterator requirements

    pointer operator->() const { return &group_; }
    reference operator*() const { return group_; }

    group_by_iterator& operator++() {
        advance_group();
        return *this;
    }

    group_by_iterator operator++(int) {
        group_by_iterator c(*this);
        ++*this;
        return c;
    }

    bool operator==(const group_by_iterator& x) const {
        return (group_.empty() && x.group_.empty()) || group_.left==x.group_.left;
    }

    bool operator!=(const group_by_iterator& x) const {
        return !(*this==x);
    }
};

template <typename Seq, typename Eq = std::equal_to<>>
auto group_by(Seq&& seq, Eq eq = Eq{}) {
    using std::begin;
    using std::end;

    using I = decltype(std::begin(seq));
    using S = decltype(std::end(seq));

    return make_range(
            group_by_iterator<I, S, Eq>(std::begin(seq), std::end(seq), eq),
            group_by_iterator<I, S, Eq>(eq));
}

} // namespace util
} // namespace arb
