#include <iterator>
#include <numeric>
#include <utility>
#include <vector>

#include <arbor/cable_cell.hpp>
#include <arbor/morph/morphology.hpp>
#include <arbor/morph/locset.hpp>

#include "morph/em_morphology.hpp"
#include "util/filter.hpp"
#include "util/rangeutil.hpp"
#include "util/span.hpp"

#include "common.hpp"
#include "common_morphologies.hpp"
#include "../common_cells.hpp"

using namespace arb;
using namespace common_morphology;
using util::make_span;

namespace {

template <typename... A>
locset as_locset(mlocation head, A... tail) {
    return join(ls::location(head), ls::location(tail)...);
}

template <typename Seq>
locset as_locset(const Seq& seq) {
    using std::begin;
    using std::end;
    return std::accumulate(begin(seq), end(seq), ls::nil(),
        [](locset ls, const mlocation& p) { return join(std::move(ls), ls::location(p)); });
}

} // anonymous namespace

TEST(cv_policy, explicit_policy) {
    using L = mlocation;
    locset lset = as_locset(L{0, 0},  L{0, 0.5},  L{0, 1.},  L{1,  0.5},  L{4, 0.2});

    cv_policy pol = cv_policy_explicit(lset);
    for (auto& m: {m_sph_b6, m_reg_b6, m_mlt_b6}) {
        cable_cell cell(m);

        locset result = pol.cv_boundary_points(cell);
        EXPECT_EQ(thingify(lset, *cell.morphology()), thingify(result, *cell.morphology()));
    }
}

TEST(cv_policy, empty_morphology) {
    // Any policy applied to an empty morphology should give an empty locset,
    // with the exception of cv_policy_explicit (this is still being debated).

    using namespace cv_policy_flag;

    cv_policy policies[] = {
        cv_policy_fixed_per_branch(3),
        cv_policy_fixed_per_branch(3, single_root_cv|interior_forks),
        cv_policy_max_extent(0.234),
        cv_policy_max_extent(0.234, single_root_cv|interior_forks)
    };

    cable_cell cell(m_empty);
    auto empty_loclist = thingify(ls::nil(), *cell.morphology());

    for (auto& pol: policies) {
        EXPECT_EQ(empty_loclist, thingify(pol.cv_boundary_points(cell), *cell.morphology()));
    }
}

TEST(cv_policy, single_root_cv) {
    // For policies that respect the single_root_cv flag, the boundary points should
    // be the same as if the flag were not provided, except that the points on branch 0
    // should only ever be (0, 0) and (0, 1) if the morphology is not empty.

    using namespace cv_policy_flag;

    std::pair<cv_policy, cv_policy> policy_pairs[] = {
        {cv_policy_fixed_per_branch(3),                 cv_policy_fixed_per_branch(3, single_root_cv)},
        {cv_policy_fixed_per_branch(3, interior_forks), cv_policy_fixed_per_branch(3, single_root_cv|interior_forks)},
        {cv_policy_max_extent(0.234),                   cv_policy_max_extent(0.234, single_root_cv)},
        {cv_policy_max_extent(0.234, interior_forks),   cv_policy_max_extent(0.234, single_root_cv|interior_forks)}
    };

    for (auto& morph: {m_sph_b1, m_reg_b1, m_sph_b6, m_reg_b6, m_mlt_b6}) {
        cable_cell cell(morph);

        for (auto& polpair: policy_pairs) {
            mlocation_list p1 = thingify(polpair.first.cv_boundary_points(cell), *cell.morphology());
            mlocation_list p2 = thingify(polpair.second.cv_boundary_points(cell), *cell.morphology());

            auto p1_no_b0 = util::filter(p1, [](mlocation l) { return l.branch>0; });
            mlocation_list expected = {{0,0}, {0,1}};
            expected.insert(expected.end(), p1_no_b0.begin(), p1_no_b0.end());

            EXPECT_EQ(expected, p2);
        }
    }
}

TEST(cv_policy, fixed_per_branch) {
    using namespace cv_policy_flag;
    using L = mlocation;

    // root branch only
    for (auto& morph: {m_sph_b1, m_reg_b1}) {
        cable_cell cell(morph);
        {
            // boundary fork points
            cv_policy pol = cv_policy_fixed_per_branch(4);
            locset points = pol.cv_boundary_points(cell);
            locset expected = as_locset(L{0, 0}, L{0, 0.25}, L{0, 0.5}, L{0, 0.75}, L{0, 1});
            EXPECT_EQ(thingify(expected, *cell.morphology()), thingify(points, *cell.morphology()));
        }
        {
            // interior fork points
            cv_policy pol = cv_policy_fixed_per_branch(4, interior_forks);
            locset points = pol.cv_boundary_points(cell);
            locset expected = as_locset(L{0, 0.125}, L{0, 0.375}, L{0, 0.625}, L{0, 0.875});
            EXPECT_EQ(thingify(expected, *cell.morphology()), thingify(points, *cell.morphology()));
        }
    }

    // spherical root, six branches and multiple top level branches cases:
    // expected points are the same.
    for (auto& morph: {m_sph_b6, m_mlt_b6}) {
        cable_cell cell(morph);

        {
            // boundary fork points
            cv_policy pol = cv_policy_fixed_per_branch(2);
            locset points = pol.cv_boundary_points(cell);
            locset expected = as_locset(
                L{0, 0}, L{0, 0.5}, L{0,1}, L{1, 0}, L{1, 0.5}, L{1,1}, L{2, 0}, L{2, 0.5}, L{2,1},
                L{3, 0}, L{3, 0.5}, L{3,1}, L{4, 0}, L{4, 0.5}, L{4,1}, L{5, 0}, L{5, 0.5}, L{5,1}
            );
            EXPECT_EQ(thingify(expected, *cell.morphology()), thingify(points, *cell.morphology()));
        }
        {
            // interior fork points
            cv_policy pol = cv_policy_fixed_per_branch(2, interior_forks);
            locset points = pol.cv_boundary_points(cell);
            locset expected = as_locset(
                L{0, 0.25}, L{0, 0.75}, L{1, 0.25}, L{1, 0.75}, L{2, 0.25}, L{2, 0.75},
                L{3, 0.25}, L{3, 0.75}, L{4, 0.25}, L{4, 0.75}, L{5, 0.25}, L{5, 0.75}
            );
            EXPECT_EQ(thingify(expected, *cell.morphology()), thingify(points, *cell.morphology()));
        }
    }
}

TEST(cv_policy, max_extent) {
    using namespace cv_policy_flag;
    using L = mlocation;

    // root branch only
    for (auto& morph: {m_sph_b1, m_reg_b1}) {
        cable_cell cell(morph);
        ASSERT_EQ(1.0, cell.morphology()->branch_length(0));

        {
            // extent of 0.25 should give exact fp calculation, giving
            // 4 CVs on the root branch.
            cv_policy pol = cv_policy_max_extent(0.25);
            locset points = pol.cv_boundary_points(cell);
            locset expected = as_locset(L{0, 0}, L{0, 0.25}, L{0, 0.5}, L{0, 0.75}, L{0, 1});
            EXPECT_EQ(thingify(expected, *cell.morphology()), thingify(points, *cell.morphology()));
        }
        {
            cv_policy pol = cv_policy_max_extent(0.25, interior_forks);
            locset points = pol.cv_boundary_points(cell);
            locset expected = as_locset(L{0, 0.125}, L{0, 0.375}, L{0, 0.625}, L{0, 0.875});
            EXPECT_EQ(thingify(expected, *cell.morphology()), thingify(points, *cell.morphology()));
        }
    }

    // cell with varying branch lengths; extent not exact fraction.
    {
        cable_cell cell(m_mlt_b6);
        ASSERT_EQ(1.0, cell.morphology()->branch_length(0));
        ASSERT_EQ(1.0, cell.morphology()->branch_length(1));
        ASSERT_EQ(2.0, cell.morphology()->branch_length(2));
        ASSERT_EQ(4.0, cell.morphology()->branch_length(3));
        ASSERT_EQ(1.0, cell.morphology()->branch_length(4));
        ASSERT_EQ(2.0, cell.morphology()->branch_length(5));

        {
            // max extent of 0.6 should give two CVs on branches of length 1,
            // four CVs on branches of length 2, and seven CVs on the branch of length 4.
            cv_policy pol = cv_policy_max_extent(0.6);
            mlocation_list points = thingify(pol.cv_boundary_points(cell), *cell.morphology());

            mlocation_list points_b012 = util::assign_from(util::filter(points, [](mlocation l) { return l.branch<3; }));
            mlocation_list expected_b012 = {
                {0, 0},  {0, 0.5},  {0, 1},
                {1, 0},  {1, 0.5},  {1, 1},
                {2, 0},  {2, 0.25}, {2, 0.5}, {2, 0.75}, {2, 1}
            };
            EXPECT_EQ(expected_b012, points_b012);

            mlocation_list points_b3 = util::assign_from(util::filter(points, [](mlocation l) { return l.branch==3; }));
            EXPECT_EQ(8u, points_b3.size());
        }
    }
}

