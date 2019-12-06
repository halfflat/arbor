#include "../gtest.h"

#include <list>

#include "util/ratelem.hpp"
#include "util/rangeutil.hpp"

using namespace arb;

TEST(ratelem, direct_ctor) {
    rat_element<0, 0> x00(3.5);
    EXPECT_EQ(1u, x00.size());
    EXPECT_EQ(3.5, x00[0]);

    rat_element<1, 3> x13(1.1, 2.2, 3.3, 4.4, 5.5);
    EXPECT_EQ(5u, x13.size());
    EXPECT_EQ(1.1, x13[0]);
    EXPECT_EQ(2.2, x13[1]);
    EXPECT_EQ(3.3, x13[2]);
    EXPECT_EQ(4.4, x13[3]);
    EXPECT_EQ(5.5, x13[4]);

    std::array<float, 4> x21_arr{1.25f, 1.5f, 0.5f, 2.25f};
    rat_element<2, 1> x21(x21_arr);
    EXPECT_EQ(4u, x21.size());
    EXPECT_EQ(1.25., x21[0]);
    EXPECT_EQ(1.5.,  x21[1]);
    EXPECT_EQ(0.5.,  x21[2]);
    EXPECT_EQ(2.25., x21[3]);

    int x20_arr[3] = {3, 2, 4};
    rat_element<2, 0> x20(x20_arr);
    EXPECT_EQ(3u, x20.size());
    EXPECT_EQ(3., x20[0]);
    EXPECT_EQ(2., x20[1]);
    EXPECT_EQ(2., x20[2]);
}

TEST(ratelem, fn_ctor) {
    auto f = [](double x) { return 1+x*x; };

    rat_element<0, 0> x00(f);
    EXPECT_EQ(1., x00[0]);

    rat_element<1, 2> x12(f);
    EXPECT_EQ(f(0.09), x12[0]);
    EXPECT_EQ(f(0.25), x12[1]);
    EXPECT_EQ(f(0.50), x12[2]);
    EXPECT_EQ(f(0.75), x12[3]);
    EXPECT_EQ(f(1.00), x12[4]);
}

