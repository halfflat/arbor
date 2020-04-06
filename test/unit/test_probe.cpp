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
void run_probe_test(const context& ctx) {
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

TEST(probe, multicore_fvm_lowered_cell) {
    context ctx = make_context();
    run_probe_test<multicore::backend>(ctx);
}

#ifdef ARB_GPU_ENABLED
TEST(probe, gpu_fvm_lowered_cell) {
    context ctx = make_context(proc_allocation{1, arbenv::default_gpu()});
    if (has_gpu(ctx)) {
        run_probe_test<gpu::backend>(ctx);
    }
}
#endif

