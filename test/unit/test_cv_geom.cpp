#include <utility>

#include <arbor/util/optional.hpp>
#include <arbor/cable_cell.hpp>
#include <arbor/morph/morphology.hpp>
#include <arbor/morph/locset.hpp>

#include "morph/em_morphology.hpp"
#include "fvm_layout.hpp"

#include "common.hpp"
#include "../common_cells.hpp"

using namespace arb;
using util::make_span;

// TODO: factor out test morphologies from this and test_cv_policy.

namespace {
    std::vector<msample> make_samples(unsigned n) {
        std::vector<msample> ms;
        for (auto i: make_span(n)) ms.push_back({{0., 0., (double)i, 0.5}, 5});
        return ms;
    }

    // Test morphologies for CV determination:
    // Samples points have radius 0.5, giving an initial branch length of 1.0
    // for morphologies with spherical roots.

    const morphology m_empty;

    // spherical root, one branch
    const morphology m_sph_b1{sample_tree(make_samples(1), {mnpos}), true};

    // regular root, one branch
    const morphology m_reg_b1{sample_tree(make_samples(2), {mnpos, 0u}), false};

    // spherical root, six branches
    const morphology m_sph_b6{sample_tree(make_samples(8), {mnpos, 0u, 1u, 0u, 3u, 4u, 4u, 4u}), true};

    // regular root, six branches
    const morphology m_reg_b6{sample_tree(make_samples(7), {mnpos, 0u, 1u, 1u, 2u, 2u, 2u}), false};

    // regular root, six branches, mutiple top level branches.
    const morphology m_mlt_b6{sample_tree(make_samples(7), {mnpos, 0u, 1u, 1u, 0u, 4u, 4u}), false};
}

TEST(cv_layout, empty) {
    cable_cell empty_cell{m_empty};
    cv_geometry geom = cv_geometry_from_ends(empty_cell, ls::nil());

    EXPECT_TRUE(geom.cv_parent.empty());
    EXPECT_TRUE(geom.cv_cables.empty());
    EXPECT_TRUE(geom.cv_cables_divs.empty());
    EXPECT_EQ(0u, geom.size());
}

TEST(cv_layout, trivial) {
    std::pair<const char*, morphology> cases[] = {
        {"m_sph_b1", m_sph_b1}, {"m_reg_b1", m_reg_b1}, {"m_sph_b6", m_sph_b6}, {"m_mlt_b6", m_mlt_b6}
    };
    for (auto& p: cases) {
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

#if 0
TEST(cv_layout, branch_handling) {
    // CVs with differing treatments of a fork point with 3 children (branches 2, 3, 4 and 5 in m_sph_b6):
    // a) Four CVs: one per branch and a zero-volume CV at fork.
    // b) Two CVs: one covering branches 2 and 4; the second covering branches 3 and 5.
    // In addition there will be a single CV for the root branch and branch 1.

    using L = mlocation;
    locset case_a = as_locset(L{1, 1}, L{2, 0}, L{2, 1}, L{3,0}, L{4, 0}, L{4, 0}, L
#endif
