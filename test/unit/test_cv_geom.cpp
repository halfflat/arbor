#include <algorithm>
#include <utility>

#include <arbor/util/optional.hpp>
#include <arbor/cable_cell.hpp>
#include <arbor/morph/morphology.hpp>
#include <arbor/morph/locset.hpp>

#include "morph/em_morphology.hpp"
#include "fvm_layout.hpp"

#include "common.hpp"
#include "common_morphologies.hpp"
#include "../common_cells.hpp"

using namespace arb;
using util::make_span;

TEST(cv_layout, empty) {
    using namespace common_morphology;

    cable_cell empty_cell{m_empty};
    cv_geometry geom = cv_geometry_from_ends(empty_cell, ls::nil());

    EXPECT_TRUE(geom.cv_parent.empty());
    EXPECT_TRUE(geom.cv_cables.empty());
    EXPECT_TRUE(geom.cv_cables_divs.empty());
    EXPECT_EQ(0u, geom.size());
}

TEST(cv_layout, trivial) {
    using namespace common_morphology;

    for (auto& p: test_morphologies) {
        if (p.second.empty()) continue;

        SCOPED_TRACE(p.first);
        cable_cell cell{p.second};
        auto em = *cell.morphology();


        // Equivalent ways of specifying one CV comprising whole cell:
        cv_geometry geom1 = cv_geometry_from_ends(cell, ls::nil());
        cv_geometry geom2 = cv_geometry_from_ends(cell, ls::terminal());

        EXPECT_EQ(1u, geom1.size());
        EXPECT_EQ(geom1.cv_cables, geom2.cv_cables);

        // These are equivalent too, if there is a single root branch.
        cv_geometry geom3 = cv_geometry_from_ends(cell, ls::root());
        cv_geometry geom4 = cv_geometry_from_ends(cell, join(ls::root(), ls::terminal()));

        EXPECT_EQ(geom3.cv_cables, geom4.cv_cables);
        if (em.branch_children(mnpos).size()==1) {
            EXPECT_EQ(geom1.cv_cables, geom4.cv_cables);
        }

        mcable_list all_cables = thingify(reg::all(), em);
        EXPECT_TRUE(testing::seq_eq(all_cables, geom1.cables(0)));
    }
}

TEST(cv_layout, one_cv_per_branch) {
    using namespace common_morphology;

    for (auto& p: test_morphologies) {
        if (p.second.empty()) continue;

        SCOPED_TRACE(p.first);
        cable_cell cell{p.second};
        auto em = *cell.morphology();

        cv_geometry geom = cv_geometry_from_ends(cell, ls::on_branches(0));

        // Expect trivial CVs at every fork point, and single-cable CVs for each branch.
        std::vector<unsigned> seen_branches(em.num_branches(), 0);
        auto n_branch_child = [&em](msize_t b) { return em.branch_children(b).size(); };
        for (auto i: make_span(geom.size())) {
            auto cables = geom.cables(i);

            ASSERT_EQ(1u, cables.size());
            auto c = cables.front();

            if (c.prox_pos==c.dist_pos) {
                if (c.branch==0) {
                    EXPECT_EQ(0., c.prox_pos);
                    EXPECT_TRUE(n_branch_child(mnpos)>1);
                }
                else {
                    EXPECT_EQ(1., c.prox_pos);
                    EXPECT_TRUE(n_branch_child(c.branch)>1);
                }
            }
            else {
                ++seen_branches[c.branch];
                EXPECT_EQ(1., seen_branches[c.branch]);
                EXPECT_EQ(0., c.prox_pos);
                EXPECT_EQ(1., c.dist_pos);
            }
        }

        EXPECT_TRUE(std::find(seen_branches.begin(), seen_branches.end(), 0)==seen_branches.end());
    }
}

#if 0
TEST(cv_layout, branch_handling) {
    // CVs with differing treatments of a fork point with 3 children (branches 2, 3, 4 and 5 in m_sph_b6):
    // a) Four CVs: one per branch and a zero-volume CV at fork.
    // b) Two CVs: one covering branches 2 and 4; the second covering branches 3 and 5.
    // In addition there will be a single CV for the root branch and branch 1.

    using L = mlocation;
    locset case_a = as_locset(L{1, 1}, L{2, 0}, L{2, 1}, L{3,0}, L{4, 0}, L{4, 0}, L
#endif
