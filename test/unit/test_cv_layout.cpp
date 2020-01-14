#include <algorithm>
#include <utility>

#include <arbor/util/optional.hpp>
#include <arbor/cable_cell.hpp>
#include <arbor/morph/morphology.hpp>
#include <arbor/morph/locset.hpp>

#include "fvm_layout.hpp"
#include "util/span.hpp"

#include "common.hpp"
#include "common_morphologies.hpp"
#include "../common_cells.hpp"

using namespace arb;
using util::make_span;

TEST(cv_layout, empty) {
    using namespace common_morphology;

    cable_cell empty_cell{m_empty};
    fvm_cv_discretization D = fvm_cv_discretize(empty_cell, neuron_parameter_defaults);

    EXPECT_TRUE(D.empty());
    EXPECT_EQ(0u, D.size());
    EXPECT_EQ(1u, D.n_cell());

    EXPECT_EQ(0u, D.face_conductance.size());
    EXPECT_EQ(0u, D.cv_area.size());
    EXPECT_EQ(0u, D.cv_capacitance.size());
    EXPECT_EQ(0u, D.init_membrane_potential.size());
    EXPECT_EQ(0u, D.temperature_K.size());
    EXPECT_EQ(0u, D.diam_um.size());
}

TEST(cv_layout, trivial) {
    using namespace common_morphology;

    auto params = neuron_parameter_defaults;
    params.discretization = cv_policy_explicit(ls::nil());

    // For each cell, check size, confirm area is morphological area from
    // embedding, and that membrane-properties are equal to defaults.

    std::vector<cable_cell> cells;
    unsigned n_cv = 0;
    for (auto& p: test_morphologies) {
        cells.emplace_back(p.second);
        n_cv += !p.second.empty(); // one cv per non-empty cell
    }

    auto n_cells = cells.size();
    fvm_cv_discretization D = fvm_cv_discretize(cells, params);

    EXPECT_EQ(n_cv, D.size());
    for (unsigned i = 0; i<n_cells; ++i) {
        auto cv_indices = util::make_span(D.geometry.cell_cv_interval(i));
        if (test_morphologies[i].second.empty()) {
            ASSERT_TRUE(cv_indices.empty());
            continue;
        }
        else {
            ASSERT_EQ(1u, cv_indices.size());
        }

        auto cv = cv_indices.front();

        EXPECT_DOUBLE_EQ(params.temperature_K.value(), D.temperature_K[cv]);
        EXPECT_DOUBLE_EQ(params.init_membrane_potential.value(), D.init_membrane_potential[cv]);

        double total_area = 0;
        unsigned n_branch = cells[i].num_branches();
        const auto& embedding = cells[i].embedding();
        for (unsigned b = 0; b<n_branch; ++b) {
            total_area += embedding.integrate_area(mcable{b, 0., 1.});
        }

        EXPECT_DOUBLE_EQ(total_area, D.cv_area[cv]);
        EXPECT_DOUBLE_EQ(total_area*params.membrane_capacitance.value(), D.cv_capacitance[cv]);
    }
}

TEST(cv_layout, cable) {
    auto morph = common_morphology::m_reg_b1; // one branch, cable constant radius.

    auto params = neuron_parameter_defaults;
    params.init_membrane_potential = 0;

    cable_cell c(morph);
    c.allow_partial_paint();
    c.paint(reg::cable({0, 0.0, 0.2}), init_membrane_potential{10});
    c.paint(reg::cable({0, 0.2, 0.7}), init_membrane_potential{20});
    c.paint(reg::cable({0, 0.7, 1.0}), init_membrane_potential{30});

    params.discretization = cv_policy_explicit(ls::nil());
    fvm_cv_discretization D = fvm_cv_discretize(c, params);

    ASSERT_EQ(1u, D.size());
    EXPECT_DOUBLE_EQ(0.2*10+0.5*20+0.3*30, D.init_membrane_potential[0]);

    params.discretization = cv_policy_explicit(ls::location({0, 0.3}));
    D = fvm_cv_discretize(c, params);

    ASSERT_EQ(2u, D.size());
    EXPECT_DOUBLE_EQ((0.2*10+0.1*20)/0.3, D.init_membrane_potential[0]);
    EXPECT_DOUBLE_EQ((0.4*20+0.3*30)/0.7, D.init_membrane_potential[1]);
}

TEST(cv_layout, cable_conductance) {
    auto morph = common_morphology::m_reg_b1; // one branch, cable constant radius.
    const double rho = 5.; // [Ω·cm]

    auto params = neuron_parameter_defaults;
    params.axial_resistivity = rho;

    cable_cell c(morph);
    double radius = c.embedding().radius(mlocation{0, 0.5});
    double length = c.embedding().branch_length(0);

    params.discretization = cv_policy_explicit(ls::location({0, 0.3}));
    fvm_cv_discretization D = fvm_cv_discretize(c, params);

    ASSERT_EQ(2u, D.size());

    // Face conductance should be conductance between (relative) points 0.15 and 0.65.
    double xa = math::pi<double>*radius*radius; // [µm^2]
    double l = (0.65-0.15)*length; // [µm]
    double sigma = 100 * xa/(l*rho); // [µS]

    EXPECT_DOUBLE_EQ(0., D.face_conductance[0]);
    EXPECT_DOUBLE_EQ(sigma, D.face_conductance[1]);
}
