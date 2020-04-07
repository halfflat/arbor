#include "../gtest.h"

#include <arbor/cable_cell.hpp>
#include <arbor/common_types.hpp>
#include <arbor/version.hpp>
#include <arborenv/gpu_env.hpp>

#include "backends/event.hpp"
#include "backends/multicore/fvm.hpp"
#ifdef ARB_GPU_ENABLED
#include "backends/gpu/fvm.hpp"
#endif
#include "fvm_lowered_cell_impl.hpp"
#include "memory/cuda_wrappers.hpp"
#include "util/rangeutil.hpp"

#include "common.hpp"
#include "../common_cells.hpp"
#include "../simple_recipes.hpp"

using namespace arb;

using multicore_fvm_cell = fvm_lowered_cell_impl<multicore::backend>;
using multicore_shared_state = multicore::backend::shared_state;
ACCESS_BIND(std::unique_ptr<multicore_shared_state> multicore_fvm_cell::*, multicore_fvm_state_ptr, &multicore_fvm_cell::state_);

template <typename Backend>
struct backend_access {
    using fvm_cell = multicore_fvm_cell;

    static multicore_shared_state& state(fvm_cell& cell) {
        return *(cell.*multicore_fvm_state_ptr).get();
    }

    static fvm_value_type deref(const fvm_value_type* p) { return *p; }
};

#ifdef ARB_GPU_ENABLED

using gpu_fvm_cell = fvm_lowered_cell_impl<gpu::backend>;
using gpu_shared_state = gpu::backend::shared_state;
ACCESS_BIND(std::unique_ptr<gpu_shared_state> gpu_fvm_cell::*, gpu_fvm_state_ptr, &gpu_fvm_cell::state_);

template <>
struct backend_access<gpu::backend> {
    using fvm_cell = gpu_fvm_cell;

    static gpu_shared_state& state(fvm_cell& cell) {
        return *(cell.*gpu_fvm_state_ptr).get();
    }

    static fvm_value_type deref(const fvm_value_type* p) {
        fvm_value_type r;
        memory::cuda_memcpy_d2h(&r, p, sizeof(r));
        return r;
    }
};

#endif

template <typename Backend>
void run_v_i_probe_test(const context& ctx) {
    using fvm_cell = typename backend_access<Backend>::fvm_cell;
    auto deref = [](const fvm_value_type* p) { return backend_access<Backend>::deref(p); };

    cable_cell bs = make_cell_ball_and_stick(false);

    i_clamp stim(0, 100, 0.3);
    bs.place(mlocation{1, 1}, stim);

    cable1d_recipe rec(bs);

    mlocation loc0{0, 0};
    mlocation loc1{1, 1};
    mlocation loc2{1, 0.3};

    rec.add_probe(0, 10, cell_probe_membrane_voltage{loc0});
    rec.add_probe(0, 20, cell_probe_membrane_voltage{loc1});
    rec.add_probe(0, 30, cell_probe_total_ionic_current_density{loc2});

    std::vector<target_handle> targets;
    std::vector<fvm_index_type> cell_to_intdom;
    probe_association_map<probe_handle> probe_map;

    fvm_cell lcell(*ctx);
    lcell.initialize({0}, rec, cell_to_intdom, targets, probe_map);

    EXPECT_EQ(3u, rec.num_probes(0));
    EXPECT_EQ(3u, probe_map.size());

    EXPECT_EQ(10, probe_map.at({0, 0}).tag);
    EXPECT_EQ(20, probe_map.at({0, 1}).tag);
    EXPECT_EQ(30, probe_map.at({0, 2}).tag);

    probe_handle p0 = probe_map.at({0, 0}).handle;
    probe_handle p1 = probe_map.at({0, 1}).handle;
    probe_handle p2 = probe_map.at({0, 2}).handle;

    // Expect initial probe values to be the resting potential
    // for the voltage probes (cell membrane potential should
    // be constant), and zero for the current probe.

    auto& state = backend_access<Backend>::state(lcell);
    auto& voltage = state.voltage;

    fvm_value_type resting = voltage[0];
    EXPECT_NE(0.0, resting);

    // (Probe handles are just pointers in this implementation).
    EXPECT_EQ(resting, deref(p0));
    EXPECT_EQ(resting, deref(p1));
    EXPECT_EQ(0.0, deref(p2));

    // After an integration step, expect voltage probe values
    // to differ from resting, and between each other, and
    // for there to be a non-zero current.
    //
    // First probe, at (0,0), should match voltage in first
    // compartment.

    lcell.integrate(0.01, 0.0025, {}, {});

    EXPECT_NE(resting, deref(p0));
    EXPECT_NE(resting, deref(p1));
    EXPECT_NE(deref(p0), deref(p1));
    EXPECT_NE(0.0, deref(p2));

    fvm_value_type v = voltage[0];
    EXPECT_EQ(v, deref(p0));
}

#if 0
struct cable1d_with_eventgen: cable1d_recipe {
    template <typename Seq>
    cable1d_with_eventgen(const Seq& cells, bool coalesce):
        cable1d_recipe(cells, coalesce) {}

    cable1d_with_eventgen(const cable_cell& cell, bool coalesce):
        cable1d_recipe(cell, coalesce) {}

    std::vector<event_generator> event_generators(cell_gid_type gid) {
        if (generators_on_.count(gid)) {
            return generators_on_.at(gid);
        }
        return {}
    }

    void add_event_generator(cell_gid_type gid, event_generator eg) {
        generators_on_[gid].push_back(std::move(eg));
    }

private:
    std::unordered_map<cell_gid_type, std::vector<event_generator>> generators_on_;
};
#endif

template <typename Backend>
void run_expsyn_g_probe_test(const context& ctx, bool coalesce_synapses = false) {
    using fvm_cell = typename backend_access<Backend>::fvm_cell;
    auto deref = [](const fvm_value_type* p) { return backend_access<Backend>::deref(p); };

    const double tau = 2.0;
    EXPECT_EQ(tau, global_default_catalogue()["expsyn"].parameters.at("tau").default_value);

    // Ball-and-stick cell, two synapses, both in same CV.
    mlocation loc0{1, 0.8};
    mlocation loc1{1, 1.0};

    cable_cell bs = make_cell_ball_and_stick(false);
    bs.place(loc0, "expsyn");
    bs.place(loc1, "expsyn");
    bs.default_parameters.discretization = cv_policy_fixed_per_branch(2);

    cable1d_recipe rec(bs, coalesce_synapses);
    rec.add_probe(0, 10, cell_probe_point_state{0u, "expsyn", "g"});
    rec.add_probe(0, 20, cell_probe_point_state{1u, "expsyn", "g"});

    std::vector<target_handle> targets;
    std::vector<fvm_index_type> cell_to_intdom;
    probe_association_map<probe_handle> probe_map;

    fvm_cell lcell(*ctx);
    lcell.initialize({0}, rec, cell_to_intdom, targets, probe_map);

    EXPECT_EQ(2u, rec.num_probes(0));
    EXPECT_EQ(2u, probe_map.size());

    EXPECT_EQ(10, probe_map.at({0, 0}).tag);
    EXPECT_EQ(20, probe_map.at({0, 1}).tag);

    probe_handle p0 = probe_map.at({0, 0}).handle;
    probe_handle p1 = probe_map.at({0, 1}).handle;

    // Expect initial probe values to be intial synapse g == 0.

    EXPECT_EQ(0.0, deref(p0));
    EXPECT_EQ(0.0, deref(p1));

    if (coalesce_synapses) {
        // Should be the same raw pointer!
        EXPECT_EQ(p0, p1);
    }

    // Integrate to 3 ms, with one event at 1ms to first expsyn weight 0.5,
    // and another at 2ms to second, weight 1.

    std::vector<deliverable_event> evs = {
        {1.0, targets[0], 0.5},
        {2.0, targets[1], 1.0}
    };
    const double tfinal = 3.;
    const double dt = 0.001;
    lcell.integrate(tfinal, dt, evs, {});

    fvm_value_type g0 = deref(p0);
    fvm_value_type g1 = deref(p1);

    // Expected value: weight*exp(-(t_final-t_event)/tau).
    double expected_g0 = 0.5*std::exp(-(tfinal-1.0)/tau);
    double expected_g1 = 1.0*std::exp(-(tfinal-2.0)/tau);

    const double rtol = 1e-6;
    if (coalesce_synapses) {
        EXPECT_TRUE(testing::near_relative(expected_g0+expected_g1, g0, rtol));
        EXPECT_TRUE(testing::near_relative(expected_g0+expected_g1, g1, rtol));
    }
    else {
        EXPECT_TRUE(testing::near_relative(expected_g0, g0, rtol));
        EXPECT_TRUE(testing::near_relative(expected_g1, g1, rtol));
    }
}

TEST(probe, multicore_v_i) {
    context ctx = make_context();
    run_v_i_probe_test<multicore::backend>(ctx);
}

TEST(probe, multicore_expsyn_g) {
    context ctx = make_context();
    SCOPED_TRACE("uncoalesced synapses");
    run_expsyn_g_probe_test<multicore::backend>(ctx, false);
    SCOPED_TRACE("coalesced synapses");
    run_expsyn_g_probe_test<multicore::backend>(ctx, true);
}

#ifdef ARB_GPU_ENABLED
TEST(probe, gpu_v_i) {
    context ctx = make_context(proc_allocation{1, arbenv::default_gpu()});
    if (has_gpu(ctx)) {
        run_v_i_probe_test<gpu::backend>(ctx);
    }
}

TEST(probe, gpu_expsyn_g) {
    context ctx = make_context(proc_allocation{1, arbenv::default_gpu()});
    if (has_gpu(ctx)) {
        SCOPED_TRACE("uncoalesced synapses");
        run_expsyn_g_probe_test<gpu::backend>(ctx, false);
        SCOPED_TRACE("coalesced synapses");
        run_expsyn_g_probe_test<gpu::backend>(ctx, true);
    }
}
#endif

