#include <arbor/morph/morphology.hpp>
#include <arbor/morph/mprovider.hpp>
#include <arbor/morph/primitives.hpp>
#include <arbor/morph/stitch.hpp>

#include "../test/gtest.h"
#include "morph_pred.hpp"

using namespace arb;
using testing::region_eq;

TEST(morph, stitch_none_or_one) {
    stitch_builder B;

    stitched_morphology sm0(B);
    EXPECT_TRUE(sm0.morphology().empty());

    mpoint p1{1, 2, 3, 0.5}, p2{2, 4, 5, 1.};
    B.add({"first", p1, p2, 3});
    stitched_morphology sm1(B);
    morphology m1 = sm1.morphology();

    msegment seg0 = m1.branch_segments(0).front();
    EXPECT_EQ(3, seg0.tag);
    EXPECT_EQ(p1, seg0.prox);
    EXPECT_EQ(p2, seg0.dist);

    mprovider p(m1, sm1.labels("stitch:"));
    EXPECT_TRUE(region_eq(p, "stitch:first", reg::segment(0)));
}

TEST(morph, stitch_two) {
}

TEST(morph, stitch_complex) {
}
