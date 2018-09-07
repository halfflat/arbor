#pragma once

#include <arbor/fvm_types.hpp>
#include <arbor/mc_cell.hpp>
#include <arbor/mechanism.hpp>
#include <arbor/mechinfo.hpp>
#include <arbor/mechcat.hpp>

#include "util/partition.hpp"
#include "util/span.hpp"

namespace arb {

// Discretization data for an unbranched segment.

struct segment_info {
    using value_type = fvm_value_type;
    using index_type = fvm_index_type;

    value_type parent_cv_area = 0;
    value_type distal_cv_area = 0;

    static constexpr index_type npos = -1;

    index_type parent_cv = npos; // npos => no parent.
    index_type proximal_cv = 0;  // First CV in segment, excluding parent.
    index_type distal_cv = 0;    // Last CV in segment (may be shared with other segments).

    bool has_parent() const { return parent_cv!=npos; }

    // Range of CV-indices for segment, excluding parent.
    std::pair<index_type, index_type> cv_range() const {
        return {proximal_cv, 1+distal_cv};
    }

    // Position is proportional distal distance along segment, in [0, 1).
    index_type cv_by_position(double pos) const {
        index_type n = distal_cv+1-proximal_cv;
        index_type i = static_cast<index_type>(n*pos+0.5);
        if (i>0) {
            return proximal_cv+(i-1);
        }
        else {
            return parent_cv==npos? proximal_cv: parent_cv;
        }
    }

    // Convert cv in segment to proportional distal distance.
    // Returns 0 if not in segment.
    double cv_distal_distance(index_type cv) {
        if (cv<proximal_cv || cv>distal_cv) {
            return 0;
        }

        // Special case for soma: report centre rather than
        // point used for flux approximation.
        if (cv==proximal_cv && !has_parent()) {
            return 0.5;
        }

        return double(cv+1-proximal_cv)/(distal_cv+1-proximal_cv);
    }
};

// Discretization of morphologies and electrical properties for
// cells in a cell group.

struct fvm_discretization {
    using value_type = fvm_value_type;
    using size_type = fvm_size_type;
    using index_type = fvm_index_type; // In particular, used for CV indices.

    size_type ncell;
    size_type ncomp;

    // Note: if CV j has no parent, parent_cv[j] = j. TODO: confirm!
    std::vector<index_type> parent_cv;
    std::vector<index_type> cv_to_cell;

    std::vector<value_type> face_conductance; // [µS]
    std::vector<value_type> cv_area;          // [µm²]
    std::vector<value_type> cv_capacitance;   // [pF]

    std::vector<segment_info> segments;
    std::vector<size_type> cell_segment_bounds; // Partitions segment indices by cell.
    std::vector<index_type> cell_cv_bounds;      // Partitions CV indices by cell.

    auto cell_segment_part() const {
        return util::partition_view(cell_segment_bounds);
    }

    auto cell_cv_part() const {
        return util::partition_view(cell_cv_bounds);
    }

    size_type segment_location_to_cv(size_type cell_index, segment_location segloc) const {
        auto cell_segs = cell_segment_part()[cell_index];

        size_type seg = segloc.segment+cell_segs.first;
        arb_assert(seg<cell_segs.second);
        return segments[seg].cv_by_position(segloc.position);
    }
};

fvm_discretization fvm_discretize(const std::vector<mc_cell>& cells);


// Post-discretization data for point and density mechanism instantiation.

struct fvm_mechanism_config {
    using value_type = fvm_value_type;
    using index_type = fvm_index_type;

    mechanismKind kind;

    // Ordered CV indices where mechanism is present; may contain
    // duplicates for point mechanisms.
    std::vector<index_type> cv;

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
    std::vector<value_type> iconc_norm_area;
    std::vector<value_type> econc_norm_area;
};

struct fvm_mechanism_data {
    // Mechanism config, indexed by mechanism name.
    std::unordered_map<std::string, fvm_mechanism_config> mechanisms;

    // Ion config, indexed by ionKind.
    std::unordered_map<ionKind, fvm_ion_config> ions;

    // Total number of targets (point-mechanism points)
    std::size_t ntarget = 0;
};

fvm_mechanism_data fvm_build_mechanism_data(const mechanism_catalogue& catalogue, const std::vector<mc_cell>& cells, const fvm_discretization& D);

} // namespace arb
