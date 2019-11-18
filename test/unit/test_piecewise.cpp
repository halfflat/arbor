#include "../gtest.h"

#include "util/piecewise.hpp"

using namespace arb;

TEST(piecewise, access) {
    pw_elements<int> p;

    double v[5] = {1., 1.5, 2., 2.5, 3.};
    int x[4] = {10, 8, 9, 4};
    p.assign(v, x);

    ASSERT_EQ(4u, p.size());

    EXPECT_EQ(10, p[0]);
    EXPECT_EQ( 8, p[1]);
    EXPECT_EQ( 9, p[2]);
    EXPECT_EQ( 4, p[3]);

    using dp = std::pair<double, double>;
    EXPECT_EQ(dp(1.0, 1.5), p.interval(0));
    EXPECT_EQ(dp(1.5, 2.0), p.interval(1));
    EXPECT_EQ(dp(2.0, 2.5), p.interval(2));
    EXPECT_EQ(dp(2.5, 3.0), p.interval(3));
}
