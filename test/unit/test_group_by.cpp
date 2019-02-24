#include "../gtest.h"

#include <forward_list>
#include <functional>
#include <string>
#include <vector>

#include <arbor/assert.hpp>

#include "util/group_by.hpp"
#include "common.hpp"

using namespace arb;
using namespace std::literals;

TEST(group_by, forward) {
    std::forward_list<int> fl = {1, 1, 3, 4, 5, 6, 6};
    std::forward_list<int>::iterator expected[] = {
        fl.begin(),
        std::next(fl.begin(), 2),
        std::next(fl.begin(), 3),
        std::next(fl.begin(), 4),
        std::next(fl.begin(), 5),
        fl.end()
    };

    auto g1 = util::group_by(fl, std::equal_to<>{});
    int k = 0;
    for (auto i = g1.begin(); i!=g1.end(); ++i) {
        ASSERT_LE(k, 4);
        EXPECT_EQ(expected[k], i->begin());
        EXPECT_EQ(expected[k+1], i->end());
        ++k;
    }
    ASSERT_EQ(5, k);
}

TEST(group_by, sentinel_range) {
    const char* cstr = "aaacccc";
    auto cstr_range = util::make_range(cstr, testing::null_terminated);

    std::vector<std::string> groups;
    for (auto& g: util::group_by(cstr_range)) {
        groups.emplace_back(g.begin(), g.end());
    }

    ASSERT_EQ(2u, groups.size());
    EXPECT_EQ("aaa"s, groups[0]);
    EXPECT_EQ("cccc"s, groups[1]);
}

TEST(group_by, array) {
    int a[] = {0, 1, 3, 2, 4, 1, 5};
    auto same_parity = [](int x, int y) { return ((x^y)&1)==0; };

    using ivec = std::vector<int>;
    std::vector<ivec> groups;
    for (auto& g: util::group_by(a, same_parity)) {
        groups.emplace_back(g.begin(), g.end());
    }

    ASSERT_EQ(4u, groups.size());
    EXPECT_EQ((ivec{0}), groups[0]);
    EXPECT_EQ((ivec{1, 3}), groups[1]);
    EXPECT_EQ((ivec{2, 4}), groups[2]);
    EXPECT_EQ((ivec{1, 5}), groups[3]);
}

TEST(group_by, empty) {
    std::vector<int> v;
    auto g = util::group_by(v);

    EXPECT_TRUE(g.empty());
}
