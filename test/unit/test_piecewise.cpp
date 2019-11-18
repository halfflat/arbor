#include "../gtest.h"

#include "util/piecewise.hpp"

using namespace arb;

TEST(piecewise, assign) {
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

    pw_elements<int> q1(p);
    pw_elements<int> q2;
    q2 = p;

    EXPECT_EQ((std::vector<double>{1.0, 1.5, 2.0, 2.5, 3.0}), p.vertices());
    EXPECT_EQ((std::vector<int>{10, 8, 9, 4}), p.elements());

    EXPECT_EQ(q1.vertices(), p.vertices());
    EXPECT_EQ(q2.vertices(), p.vertices());

    EXPECT_EQ(q1.elements(), p.elements());
    EXPECT_EQ(q2.elements(), p.elements());
}

TEST(piecewise, bounds) {
    // TODO
}

TEST(piecewise, push) {
    pw_elements<int> q;
    using dp = std::pair<double, double>;

    // Need left hand side!
    EXPECT_THROW(q.push_back(3.1, 4), std::runtime_error);

    q.clear();
    q.push_back(1.1, 3.1, 4);
    q.push_back(3.1, 4.3, 5);
    EXPECT_EQ(dp(1.1, 3.1), q.interval(0));
    EXPECT_EQ(dp(3.1, 4.3), q.interval(1));
    EXPECT_EQ(4, q[0]);
    EXPECT_EQ(5, q[1]);

    q.push_back(7.2, 6);
    EXPECT_EQ(dp(4.3, 7.2), q.interval(2));
    EXPECT_EQ(6, q[2]);

    // Supplied left side doesn't match current right.
    EXPECT_THROW(q.push_back(7.4, 9.1, 7), std::runtime_error);
}

TEST(piecewise, pwvoid) {
    pw_elements<> p;
    using dp = std::pair<double, double>;

    p.push_back(0.1, 0.2);
    p.push_back(0.3);
    p.push_back(0.4);

    EXPECT_EQ(3u, p.size());
    EXPECT_EQ((std::vector<double>{0.1, 0.2, 0.3, 0.4}), p.vertices());
    EXPECT_EQ(dp(0.2,0.3), p.interval(1));

    pw_elements<> q;
    q = p;
    EXPECT_EQ(3u, q.size());
    EXPECT_EQ(p.vertices(), q.vertices());
    EXPECT_EQ(dp(0.2,0.3), q.interval(1));
}

TEST(piecewise, meet) {
    pw_elements<int> p03;
    p03.assign((double [3]){0., 1.5, 3.}, (int [2]){10, 11});

    pw_elements<int> p14;
    p14.assign((double [5]){1, 2.25, 3., 3.5, 4.}, (int [4]){3, 4, 5, 6});

    pw_elements<int> p12;
    p12.assign((double [4]){1, 1.25, 1.75, 2.}, (int [3]){7, 8, 9});

    pw_elements<std::pair<int, int>> p03_14 = meet(p03, p14);
    EXPECT_EQ(1., p03_14.bounds().first);
    EXPECT_EQ(3., p03_14.bounds().second);

    EXPECT_EQ((std::vector<double>{1., 1.5, 2.25, 3.}), p03_14.vertices());
    using ii = std::pair<int, int>;
    EXPECT_EQ((std::vector<ii>{ii(10, 3), ii(11, 3), ii(11, 4)}), p03_14.elements());

    // more TODO.
}

TEST(piecewise, meet_void) {

}
