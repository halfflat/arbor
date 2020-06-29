#include <iterator>
#include <numeric>
#include <utility>
#include <vector>

#include <arbor/util/optional.hpp>
#include <arbor/cable_cell.hpp>
#include <arbor/cable_cell_param.hpp>
#include <arbor/morph/morphology.hpp>
#include <arbor/morph/locset.hpp>
#include <arbor/morph/mprovider.hpp>
#include <arbor/morph/region.hpp>

#include "util/filter.hpp"
#include "util/rangeutil.hpp"
#include "util/span.hpp"

#include "common_morphologies.hpp"
#include "morph_pred.hpp"

using namespace arb;
using util::make_span;
using testing::locset_eq;
using testing::mlocationlist_eq;

using namespace common_morphology;

namespace {
    template <typename... A>
    locset as_locset(mlocation head, A... tail) {
        return join(locset(head), locset(tail)...);
    }
}

TEST(cv_policy, explicit_policy) {
    using L = mlocation;
    locset lset = as_locset(L{0, 0},  L{0, 0.5},  L{0, 1.},  L{1,  0.5},  L{4, 0.2});

    cv_policy pol = cv_policy_explicit(lset);
    for (auto& m: {m_reg_b6, m_mlt_b6}) {
        cable_cell cell(m);

        locset result = pol.cv_boundary_points(cell);
        locset expected = join(ls::boundary(reg::all()), lset);
        EXPECT_TRUE(locset_eq(cell.provider(), expected, result));
    }

    // With cables 1 and 2, expect to pick up (1, 0.5) from locset,
    // and cable ends (1, 0), (1, 1), (2, 0), (2, 1), as the two
    // cables constitute two components.

    region b12 = join(reg::branch(1), reg::branch(2));
    pol = cv_policy_explicit(lset, b12);
    for (auto& m: {m_reg_b6, m_mlt_b6}) {
        cable_cell cell(m);

        locset result = pol.cv_boundary_points(cell);
        locset expected = as_locset(L{1, 0}, L{1, 0.5}, L{1, 1}, L{2, 0}, L{2, 1});
        EXPECT_TRUE(locset_eq(cell.provider(), expected, result));
    }

    // Taking the completion of the two cables, the boundary of the region
    // will be (0, 1), (1, 1), (2, 1) for m_mlt_b6.

    pol = cv_policy_explicit(lset, reg::complete(b12));
    for (auto& m: {m_mlt_b6}) {
        cable_cell cell(m);

        locset result = pol.cv_boundary_points(cell);
        locset expected = as_locset(L{0, 1}, L{1, 0.5}, L{1, 1}, L{2, 1});
        EXPECT_TRUE(locset_eq(cell.provider(), expected, result));
    }
}

TEST(cv_policy, empty_morphology) {
    // Any policy applied to an empty morphology should give an empty locset.

    using namespace cv_policy_flag;

    cv_policy policies[] = {
        cv_policy_fixed_per_branch(3),
        cv_policy_fixed_per_branch(3, interior_forks),
        cv_policy_max_extent(0.234),
        cv_policy_max_extent(0.234, interior_forks),
        cv_policy_single(),
        cv_policy_single(reg::all()),
        cv_policy_explicit(ls::location(0, 0))
    };

    cable_cell cell(m_empty);

    for (auto& pol: policies) {
        EXPECT_TRUE(locset_eq(cell.provider(), ls::nil(), pol.cv_boundary_points(cell)));
    }
}

TEST(cv_policy, fixed_per_branch) {
    using namespace cv_policy_flag;
    using L = mlocation;

    // Root branch only:
    {
        cable_cell cell(m_reg_b1);
        {
            // boundary fork points
            cv_policy pol = cv_policy_fixed_per_branch(4);
            locset expected = as_locset(L{0, 0}, L{0, 0.25}, L{0, 0.5}, L{0, 0.75}, L{0, 1});
            EXPECT_TRUE(locset_eq(cell.provider(), expected, pol.cv_boundary_points(cell)));
        }
        {
            // interior fork points
            cv_policy pol = cv_policy_fixed_per_branch(4, interior_forks);
            locset points = pol.cv_boundary_points(cell);
            locset expected = as_locset(L{0, 0}, L{0, 0.125}, L{0, 0.375}, L{0, 0.625}, L{0, 0.875}, L{0, 1});
            EXPECT_TRUE(locset_eq(cell.provider(), expected, pol.cv_boundary_points(cell)));
        }
    }

    // Multiple top level branches:
    // top level branches are 0 and 3, terminal branches are 1, 2, 4 and 5.
    {
        cable_cell cell(m_mlt_b6);

        {
            // With boundary fork points:
            cv_policy pol = cv_policy_fixed_per_branch(2);
            locset expected = as_locset(
                L{0, 0}, L{0, 0.5}, L{0,1}, L{1, 0}, L{1, 0.5}, L{1,1}, L{2, 0}, L{2, 0.5}, L{2,1},
                L{3, 0}, L{3, 0.5}, L{3,1}, L{4, 0}, L{4, 0.5}, L{4,1}, L{5, 0}, L{5, 0.5}, L{5,1}
            );
            EXPECT_TRUE(locset_eq(cell.provider(), expected, pol.cv_boundary_points(cell)));
        }
        {
            // With interior fork points:
            cv_policy pol = cv_policy_fixed_per_branch(2, interior_forks);
            locset expected = as_locset(
                L{0, 0}, L{0, 0.25}, L{0, 0.75},
                L{1, 0.25}, L{1, 0.75}, L{1, 1.0},
                L{2, 0.25}, L{2, 0.75}, L{2, 1.0},
                L{3, 0}, L{3, 0.25}, L{3, 0.75},
                L{4, 0.25}, L{4, 0.75}, L{4, 1.0},
                L{5, 0.25}, L{5, 0.75}, L{5, 1.0}
            );
            EXPECT_TRUE(locset_eq(cell.provider(), expected, pol.cv_boundary_points(cell)));
        }
    }

    // Restrict to an incomplete subtree (distal half of branch 0 and all of branch 2)
    // in m_mlt_b6 morphology.
    //
    // With two per branch and interior forks, expect to see:
    //      (0, 0.5), (0, 0.625), (0.0875) on branch 0;
    //      (2, 0.25), (2, 0.75), (2, 1.) on branch 2;
    //      (1, 0) on branch 1.
    {
        cable_cell cell(m_mlt_b6);

        region reg = mcable_list{{0, 0.5, 1.}, {2, 0., 1.}};
        cv_policy pol = cv_policy_fixed_per_branch(2, reg, interior_forks);
        locset expected = as_locset(
            L{0, 0.5}, L{0, 0.625}, L{0, 0.875},
            L{1, 0},
            L{2, 0.25}, L{2, 0.75}, L{2, 1}
        );
        EXPECT_TRUE(locset_eq(cell.provider(), expected, pol.cv_boundary_points(cell)));
    }
}

TEST(cv_policy, max_extent) {
    using namespace cv_policy_flag;
    using L = mlocation;

    // Root branch only:
    {
        cable_cell cell(m_reg_b1);
        ASSERT_EQ(1.0, cell.embedding().branch_length(0));

        {
            // extent of 0.25 should give exact fp calculation, giving
            // 4 CVs on the root branch.
            cv_policy pol = cv_policy_max_extent(0.25);
            locset expected = as_locset(L{0, 0}, L{0, 0.25}, L{0, 0.5}, L{0, 0.75}, L{0, 1});
            EXPECT_TRUE(locset_eq(cell.provider(), expected, pol.cv_boundary_points(cell)));
        }
        {
            cv_policy pol = cv_policy_max_extent(0.25, interior_forks);
            locset expected = as_locset(L{0, 0}, L{0, 0.125}, L{0, 0.375}, L{0, 0.625}, L{0, 0.875}, L{0, 1});
            EXPECT_TRUE(locset_eq(cell.provider(), expected, pol.cv_boundary_points(cell)));
        }
    }

    // Cell with varying branch lengths; extent not exact fraction:
    {
        cable_cell cell(m_mlt_b6);
        ASSERT_EQ(1.0, cell.embedding().branch_length(0));
        ASSERT_EQ(1.0, cell.embedding().branch_length(1));
        ASSERT_EQ(2.0, cell.embedding().branch_length(2));
        ASSERT_EQ(4.0, cell.embedding().branch_length(3));
        ASSERT_EQ(1.0, cell.embedding().branch_length(4));
        ASSERT_EQ(2.0, cell.embedding().branch_length(5));

        {
            // Max extent of 0.6 should give two CVs on branches of length 1,
            // four CVs on branches of length 2, and seven CVs on the branch of length 4.
            cv_policy pol = cv_policy_max_extent(0.6);
            mlocation_list points = thingify(pol.cv_boundary_points(cell), cell.provider());

            mlocation_list points_b012 = util::assign_from(util::filter(points, [](mlocation l) { return l.branch<3; }));
            mlocation_list expected_b012 = {
                {0, 0},  {0, 0.5},  {0, 1},
                {1, 0},  {1, 0.5},  {1, 1},
                {2, 0},  {2, 0.25}, {2, 0.5}, {2, 0.75}, {2, 1}
            };
            EXPECT_TRUE(mlocationlist_eq(expected_b012, points_b012));

            mlocation_list points_b3 = util::assign_from(util::filter(points, [](mlocation l) { return l.branch==3; }));
            EXPECT_EQ(8u, points_b3.size());
        }
    }
}

TEST(cv_policy, every_sample) {
    using namespace cv_policy_flag;

    // Cell with root branch and two child branches, with multiple samples per branch.
    // Fork is at (0., 0., 4.0).

    std::vector<msample> ms;

    ms.push_back({{  0.,   0., 0., 0.5}, 5});
    for (auto i: make_span(4)) ms.push_back({{  0.,   0., i+1., 0.5}, 5});
    for (auto i: make_span(4)) ms.push_back({{  0., i+1.,  4.0, 0.5}, 5});
    for (auto i: make_span(4)) ms.push_back({{i+1.,    0,  4.0, 0.5}, 5});

    std::vector<msize_t> parents = {mnpos, 0, 1, 2, 3, 4, 5, 6, 7, 4, 9, 10, 11 };
    morphology m{sample_tree(ms, parents), false};

    // Including all samples:
    {
        cable_cell cell(m);
        cv_policy pol = cv_policy_every_sample();

        mlocation_list expected = {
            {0, 0}, {0, 0.25}, {0, 0.5}, {0, 0.75}, {0, 1.},
            {1, 0}, {1, 0.25}, {1, 0.5}, {1, 0.75}, {1, 1.},
            {2, 0}, {2, 0.25}, {2, 0.5}, {2, 0.75}, {2, 1.}
        };

        EXPECT_TRUE(locset_eq(cell.provider(), locset(expected), pol.cv_boundary_points(cell)));
    }
}

TEST(cv_policy, domain) {
    

}
