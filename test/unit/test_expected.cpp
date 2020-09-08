#include "../gtest.h"

#include <arbor/util/expected.hpp>

#include "common.hpp"

using namespace arb::util;

TEST(expected, ctors) {
    // Test constructors, verify against bool conversion and
    // access via value() and error().

    struct int3 {
        int v = 3;
        int3() = default;
        int3(int a, int b, int c): v(a+b+c) {}
    };

    {
        // Default construction.

        expected<int3, int3> x;
        EXPECT_TRUE(x);
        EXPECT_EQ(3, x.value().v);
    }
    {
        // Default void construction.

        expected<void, int3> x;
        EXPECT_TRUE(x);
    }
    {
        // In-place construction.

        expected<int3, int3> x(in_place, 1, 2, 3);
        EXPECT_TRUE(x);
        EXPECT_EQ(6, x.value().v);
    }
    {
        // From-value construction.

        int3 v;
        v.v = 19;
        expected<int3, int3> x(v);
        EXPECT_TRUE(x);
        EXPECT_EQ(19, x.value().v);
    }
    {
        // From-unexpected construction.

        int3 v;
        v.v = 19;
        expected<int3, int3> x(make_unexpected(v));
        EXPECT_FALSE(x);
        EXPECT_EQ(19, x.error().v);
        EXPECT_THROW(x.value(), bad_expected_access<int3>);
    }
    {
        // From-unexpected void construction.

        int3 v;
        v.v = 19;
        expected<void, int3> x(make_unexpected(v));
        EXPECT_FALSE(x);
        EXPECT_EQ(19, x.error().v);
    }
    {
        // In-place unexpected construction.

        expected<int3, int3> x(unexpect, 1, 2, 3);
        EXPECT_FALSE(x);
        EXPECT_EQ(6, x.error().v);
        EXPECT_THROW(x.value(), bad_expected_access<int3>);
    }
    {
        // In-place void unexpected construction.

        expected<void, int3> x(unexpect, 1, 2, 3);
        EXPECT_FALSE(x);
        EXPECT_EQ(6, x.error().v);
    }
    {
        // Conversion from other expect.

        struct X {};
        struct Y {};
        struct Z {
            int v = 0;
            Z(const X&): v(1) {}
            Z(const Y&): v(2) {}
            Z(X&&): v(-1) {}
            Z(Y&&): v(-2) {}
        };

        expected<X, Y> x;
        expected<Z, Z> y(x);
        EXPECT_TRUE(y);
        EXPECT_EQ(1, y.value().v);

        expected<Z, Z> my(std::move(x));
        EXPECT_TRUE(my);
        EXPECT_EQ(-1, my.value().v);

        expected<X, Y> xu(unexpect);
        expected<Z, Z> yu(xu);
        EXPECT_FALSE(yu);
        EXPECT_EQ(2, yu.error().v);

        expected<Z, Z> myu(std::move(xu));
        EXPECT_FALSE(myu);
        EXPECT_EQ(-2, myu.error().v);
    }
    
}
