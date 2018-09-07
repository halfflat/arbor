#pragma once

#include <arbor/common_types.hpp>
#include <arbor/fvm_types.hpp>

// Structures for the representation of event delivery targets and
// staged events.

namespace arb {

// Post-synaptic spike events

struct target_handle {
    cell_local_size_type mech_id;    // mechanism type identifier (per cell group).
    cell_local_size_type mech_index; // instance of the mechanism
    cell_size_type cell_index;       // which cell (acts as index into e.g. vec_t)

    target_handle() {}
    target_handle(cell_local_size_type mech_id, cell_local_size_type mech_index, cell_size_type cell_index):
        mech_id(mech_id), mech_index(mech_index), cell_index(cell_index) {}
};

struct deliverable_event {
    time_type time;
    target_handle handle;
    float weight;

    deliverable_event() {}
    deliverable_event(time_type time, target_handle handle, float weight):
        time(time), handle(handle), weight(weight)
    {}
};

// Stream index accessor function for multi_event_stream:
inline cell_size_type event_index(const deliverable_event& ev) {
    return ev.handle.cell_index;
}

// Subset of event information required for mechanism delivery.
struct deliverable_event_data {
    cell_local_size_type mech_id;    // same as target_handle::mech_id
    cell_local_size_type mech_index; // same as target_handle::mech_index
    float weight;
};

// Delivery data accessor function for multi_event_stream:
inline deliverable_event_data event_data(const deliverable_event& ev) {
    return {ev.handle.mech_id, ev.handle.mech_index, ev.weight};
}


// Sample events (scalar values)

struct probe_handle {
    const fvm_value_type* data;
    const fvm_value_type* weight; // nullptr => no weights to apply
    unsigned count;
};

struct raw_probe_info {
    probe_handle handle;
    sample_size_type t_offset;  // offset into array to store sample time
    sample_size_type v_offset;  // offset into array to store raw probed values
};

struct sample_event {
    time_type time;
    cell_size_type cell_index;  // which cell probe is on
    raw_probe_info raw;         // event payload: what gets put where on sample
};

inline raw_probe_info event_data(const sample_event& ev) {
    return ev.raw;
}

inline cell_size_type event_index(const sample_event& ev) {
    return ev.cell_index;
}


} // namespace arb
