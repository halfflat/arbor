#pragma once

#include <arbor/fvm_types.hpp>
#include <arbor/cable_cell.hpp>

#include "util/partition.hpp"
#include "util/range.hpp"
#include "util/rangeutil.hpp"

namespace arb {

struct fvm_policy {
    // Discretization policy. Can be overridden by explicit
    // CV end points given on a cable cell. (This will get
    // moved into cable cell global properties.)
    enum {
        cvs_per_branch,
        fixed_dx
    };
    double value;
};

// Discretized cell geometry.

struct mpoint {
    fvm_size_type branch;
    fvm_value_type pos;
};

struct cell_mpoint {
    fvm_index_type cell; // index wrt some vector, not gid.
    mpoint point;
};

struct cv_geometry {
    using size_type = fvm_size_type;

    // In the partiion of CV end points given below, the first
    // point in the interal for a given CV is the most proximal end point,
    // while the remainder constitute the distal end points.

    std::vector<mpoint> cv_ends; // boundary point list for CVs.
    std::vector<size_type> cv_ends_divs; // partions cv_ends by CV index on this cell.

    // Return CV index (on this cell) containing point.
    size_type point_to_cv(mpoint) const;

    // Return end point set for given CV index (on this cell).
    auto end_points(size_type) const {
        auto part = util::partition_view(cv_ends_divs);
        return util::subrange_view(cv_ends, part[i]);
    }
};

struct fvm_discretization2 {
    using value_type = fvm_value_type;
    using size_type = fvm_size_type;
    using index_type = fvm_index_type; // In particular, used for CV indices.

    size_type ncell;
    size_type ncv;

    // Forest of CVs: parent_cv[i] == i implies i is a root.
    std::vector<index_type> parent_cv;

    // CV to cell index (not gid) mapping.
    std::vector<index_type> cv_to_cell;

    // Partitions CV indices by cell index.
    std::vector<index_type> cell_cv_divs;

    // Return half-open interval for CVs by cell index.
    std::pair<index_type, index_type> cell_cvs(size_type i) const {
        auto part = util::partition_view(cell_cv_divs);
        return part[i];
    }

    // Per-cell CV geometry.
    std::vector<cv_geometry> cell_cv_geometry;
};

struct fvm_phys_config {
    // All fields indexed by CV.
    std::vector<value_type> face_conductance; // [µS]
    std::vector<value_type> cv_area;          // [µm²]
    std::vector<value_type> cv_capacitance;   // [pF]
    std::vector<value_type> init_membrane_potential; // [mV]
    std::vector<value_type> temperature_K;    // [K]
};

#if 0
// These are the same as defined in fvm_layout.hpp

struct fvm_mechanism_config {
    using value_type = fvm_value_type;
    using index_type = fvm_index_type;

    mechanismKind kind;

    // Ordered CV indices where mechanism is present; may contain
    // duplicates for point mechanisms.
    std::vector<index_type> cv;

    // Coalesced synapse multiplier (point mechanisms only).
    std::vector<index_type> multiplicity;

    // Normalized area contribution in corresponding CV (density mechanisms only).
    std::vector<value_type> norm_area;

    // Synapse target number (point mechanisms only).
    std::vector<index_type> target;

    // (Non-global) parameters and parameter values across the mechanism instance.
    std::vector<std::pair<std::string, std::vector<value_type>>> param_values;
};

// Post-discretization data for ion channel state.

struct fvm_ion_config {
    using value_type = fvm_value_type;
    using index_type = fvm_index_type;

    // Ordered CV indices where ion must be present.
    std::vector<index_type> cv;

    // Normalized area contribution of default concentration contribution in corresponding CV.
    std::vector<value_type> init_iconc;
    std::vector<value_type> init_econc;

    // Ion-specific (initial) reversal potential per CV.
    std::vector<value_type> init_revpot;
};

struct fvm_mechanism_data {
    // Mechanism config, indexed by mechanism name.
    std::unordered_map<std::string, fvm_mechanism_config> mechanisms;

    // Ion config, indexed by ion name.
    std::unordered_map<std::string, fvm_ion_config> ions;

    // Total number of targets (point-mechanism points)
    std::size_t ntarget = 0;
};
#endif

struct fvm_layout_config2 {
    // Physical properties, all fields indexed by CV.
    struct fvm_phys_config phys;

    // Mechanism config, indexed by mechanism name.
    std::unordered_map<std::string, fvm_mechanism_config> mechanisms;

    // Ion config, indexed by ion name.
    std::unordered_map<std::string, fvm_ion_config> ions;

    // Total number of targets (point-mechanism points)
    std::size_t ntarget = 0;
};

// Discretization procedures.

fvm_discretization2 fvm_discretize(const std::vector<cable_cell>& cells, const cable_cell_parameter_set& params, const fvm_policy&);

fvm_layout_config2 fvm_build_layout(const cable_cell_global_properties& gprop, const std::vector<cable_cell>& cells, const fvm_discretization2& D);


} // namespace arb
