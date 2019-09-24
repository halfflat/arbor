#include <vector>

#include <arbor/util/optional.hpp>
#include <arbor/cable_cell.hpp>
#include <arbor/cable_cell_param.hpp>
#include <arbor/morph/morphology.hpp>

#include "util/span.hpp"

#include "common.hpp"
#include "unit_test_catalogue.hpp"
#include "../common_cells.hpp"

using namespace arb;
using util::make_span;

namespace {
    std::vector<msample> make_samples(unsigned n) {
        std::vector<msample> ms;
        for (auto i: make_span(n)) ms.push_back({{0., 0., (double)i, 1.}, 5});
        return ms;
    }

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

    template <typename... A>
    locset as_locset(A... as) {
        return join(ls::location(as)...);
    }
}

TEST(cv_policy, explicit_policy) {
    using L = mlocation;
    locset lset = as_locset(L{0,0}, L{0,0.5}, L{0,1.}, L{1, 0.5}, L{4,0.2});

    cv_policy pol = cv_policy_explicit(lset);
    for (auto& m: {m_sph_b6, m_reg_b6, m_mlt_b6}) {
        cable_cell cell = make_cable_cell(m);

        locset result = pol.cv_boundary_points(cell);
        EXPECT_EQ(thingify(lset, *cell.morphology()), thingify(result, *cell.morphology()));
    }
}

