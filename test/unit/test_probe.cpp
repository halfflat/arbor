#include "../gtest.h"

#include <arbor/common_types.hpp>
#include <arbor/mc_cell.hpp>
#include <arbor/util/any_ptr.hpp>

#include "backends/event.hpp"
#include "backends/multicore/fvm.hpp"
#include "fvm_lowered_cell_impl.hpp"
#include "util/rangeutil.hpp"

#include "common.hpp"
#include "../common_cells.hpp"
#include "../simple_recipes.hpp"

using namespace arb;
using fvm_cell = fvm_lowered_cell_impl<multicore::backend>;
using shared_state = multicore::backend::shared_state;

ACCESS_BIND(std::unique_ptr<shared_state> fvm_cell::*, fvm_state_ptr, &fvm_cell::state_);

TEST(probe, fvm) {
    execution_context context;

    mc_cell bs = make_cell_ball_and_stick(false);

    i_clamp stim(0, 100, 0.3);
    bs.add_stimulus({1, 1}, stim);

    cable1d_recipe rec(bs);

    segment_location loc0{0, 0};
    segment_location loc1{1, 1};
    segment_location loc2{1, 0.3};

    rec.add_probe(0, 10, cell_probe_address{mc_cell_probe_kind::voltage, loc0});
    rec.add_probe(0, 20, cell_probe_address{mc_cell_probe_kind::voltage, loc1});
    rec.add_probe(0, 30, cell_probe_address{mc_cell_probe_kind::current_density, loc2});

    std::vector<target_handle> targets;
    probe_association_map<probe_handle> probe_map;

    fvm_cell lcell(context);
    lcell.initialize({0}, rec, targets, probe_map);

    EXPECT_EQ(3u, rec.num_probes(0));
    EXPECT_EQ(3u, probe_map.size());

    EXPECT_EQ(10, probe_map.at({0, 0}).tag);
    EXPECT_EQ(20, probe_map.at({0, 1}).tag);
    EXPECT_EQ(30, probe_map.at({0, 2}).tag);

    probe_handle p0 = probe_map.at({0, 0}).handle;
    probe_handle p1 = probe_map.at({0, 1}).handle;
    probe_handle p2 = probe_map.at({0, 2}).handle;

    // Check that this probes have no associated weights, and counts of one.
    EXPECT_EQ(nullptr, p0.weight);
    EXPECT_EQ(nullptr, p1.weight);
    EXPECT_EQ(nullptr, p2.weight);

    EXPECT_EQ(1u, p0.count);
    EXPECT_EQ(1u, p1.count);
    EXPECT_EQ(1u, p2.count);

    // Expect initial probe values to be the resting potential
    // for the voltage probes (cell membrane potential should
    // be constant), and zero for the current probe.

    auto& state = *(lcell.*fvm_state_ptr).get();
    auto& voltage = state.voltage;

    auto resting = voltage[0];
    EXPECT_NE(0.0, resting);

    EXPECT_EQ(resting, *p0.data);
    EXPECT_EQ(resting, *p1.data);
    EXPECT_EQ(0.0, *p2.data);

    // After an integration step, expect voltage probe values
    // to differ from resting, and between each other, and
    // for there to be a non-zero current.
    //
    // First probe, at (0,0), should match voltage in first
    // compartment.

    lcell.integrate(0.01, 0.0025, {}, {});

    EXPECT_NE(resting, *p0.data);
    EXPECT_NE(resting, *p1.data);
    EXPECT_NE(*p0.data, *p1.data);
    EXPECT_NE(0.0, *p2.data);

    EXPECT_EQ(voltage[0], *p0.data);
}

TEST(probe, fvm_wide) {
    // Test sampling of area-integrated current from each CV.

    execution_context context;

    // Make a ball and 3-stick cell with just passive mechanisms and stimulus.
    // 4 compartments per dendrite, 1 for soma.
    mc_cell c;
    c.add_soma(6.0);
    c.add_cable(0, section_kind::dendrite, 0.5, 0.3, 100);
    c.add_cable(1, section_kind::dendrite, 0.6, 0.4, 100);
    c.add_cable(2, section_kind::dendrite, 0.7, 0.5, 100);

    for (auto& seg: c.segments()) {
        seg->add_mechanism("pas");
        if (seg->is_dendrite()) seg->set_compartments(4);
    }
    // Stimulus on soma: 0.3 nA = 300 pA.
    c.add_stimulus({0, 0.5}, i_clamp(0, 40, 0.3));
    for (unsigned i = 0; i<c.num_segments(); ++i) {
        //c.segment(i)->cm = 1e-12; // Set tiny capacitance: 1pF/m²
    }
    cable1d_recipe rec(c);

    rec.add_probe(0, 0, cell_probe_address{mc_cell_probe_kind::cv_currents, {}});

    std::vector<target_handle> targets;
    probe_association_map<probe_handle> probe_map;

    fvm_cell lcell(context);
    lcell.initialize({0}, rec, targets, probe_map);

    probe_handle handle = probe_map.at({0, 0}).handle;
    auto& state = *(lcell.*fvm_state_ptr).get();
    const fvm_value_type* J = state.current_density.data();
    const fvm_value_type* A = state.cv_area.data();

    EXPECT_EQ(J, handle.data);
    EXPECT_EQ(A, handle.weight);
    EXPECT_EQ(13u, handle.count); // expect one value per CV, = 13 in total for this cell.

    auto meta = util::any_cast<const mc_cell_probe_metadata*>(probe_map.at({0, 0}).metadata);
    ASSERT_NE(nullptr, meta);
    EXPECT_EQ(mc_cell_probe_kind::cv_currents, meta->kind);
    EXPECT_EQ(13u, meta->locations.size());

    std::vector<cell_lid_type> expected_segs = {0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3};
    std::vector<cell_lid_type> location_segs;
    util::assign(location_segs,
        util::transform_view(meta->locations, [](auto sl) { return sl.segment; }));

    EXPECT_EQ(expected_segs, location_segs);

    lcell.integrate(10.0, 0.025, {}, {});

    // Expect non-zero currents per CV, but sum should be close to zero.
    for (unsigned i = 0; i<handle.count; ++i) {
        EXPECT_NE(0.0, handle.data[i]);
        std::cout << "I[" << i << "]: " << handle.data[i]*handle.weight[i] << "\n";
    }
    fvm_value_type I = std::inner_product(handle.data, handle.data+handle.count, handle.weight, 0.);
    EXPECT_DOUBLE_EQ(0.0, I);

}
