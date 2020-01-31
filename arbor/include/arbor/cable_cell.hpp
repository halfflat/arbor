#pragma once

#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include <arbor/arbexcept.hpp>
#include <arbor/cable_cell_param.hpp>
#include <arbor/common_types.hpp>
#include <arbor/constants.hpp>
#include <arbor/mechcat.hpp>
#include <arbor/morph/label_dict.hpp>
#include <arbor/morph/mcable_map.hpp>
#include <arbor/morph/mprovider.hpp>
#include <arbor/morph/morphology.hpp>
#include <arbor/morph/primitives.hpp>
#include <arbor/util/typed_map.hpp>

namespace arb {

// Pair of indexes that describe range of local indices.
// Returned by cable_cell::place() calls, so that the caller can
// refer to targets, detectors, etc on the cell.
struct lid_range {
    cell_lid_type begin;
    cell_lid_type end;
    lid_range(cell_lid_type b, cell_lid_type e):
        begin(b), end(e) {}
};

// Probe type for cell descriptions.
//
// All probes currently return data as `const double`, but this may
// change in the future. Probe kinds are interpreted as follows.
//
// membrane_voltage:
//     Voltage estimate [mV] at <location>,  possibly interpolated.
//
// axial_current:
//     Intracellular current [A] at <location>, in distal direction.
//
// mechanism_state:
//     Value of state variable <key> in mechanism <source> in CV
//     at <location>. (Peeking into point mechanisms is not yet
//     supported.)
//
// total_ionic_current_density:
//     Total current density [A/m²] across membrane excluding capacitive
//     current at <location>.
//
// ionic_current_density:
//     Current density [A/m²] across membrane attributed to the
//     ion <source> at <location>.
//
// ion_int_concentration:
// ion_ext_concentration:
//     Ionic internal (resp. external) concentration [mmol/L] of
//     ion <source> at <location>, possibly interpolated.

enum class cable_cell_probe_kind {
    membrane_voltage,
    axial_current, //  -- not supported yet!
    mechanism_state,
    total_ionic_current_density,
    ionic_current_density,
    ion_int_concentration,
    ion_ext_concentration
};

struct cell_probe_address {
    cable_cell_probe_kind kind;

    mlocation location;
    std::string source; // ion name for ion data, mechanism name for mechanism state.
    std::string key;    // mechanism state variable name.
};

// Forward declare the implementation, for PIMPL.
struct cable_cell_impl;

// Typed maps for access to painted and placed assignments:
//
// Mechanisms and initial ion data are further keyed by
// mechanism name and ion name respectively.

template <typename T>
using region_assignment =
    std::conditional_t<
        std::is_same<T, mechanism_desc>::value || std::is_same<T, initial_ion_data>::value,
        std::unordered_map<std::string, mcable_map<T>>,
        mcable_map<T>>;

template <typename T>
struct placed {
    mlocation loc;
    cell_lid_type lid;
    T item;
};

template <typename T>
using mlocation_map = std::vector<placed<T>>;

template <typename T>
using location_assignment =
    std::conditional_t<
        std::is_same<T, mechanism_desc>::value,
        std::unordered_map<std::string, mlocation_map<T>>,
        mlocation_map<T>>;

using cable_cell_region_map = static_typed_map<region_assignment,
    mechanism_desc, init_membrane_potential, axial_resistivity,
    temperature_K, membrane_capacitance, initial_ion_data>;

using cable_cell_location_map = static_typed_map<location_assignment,
    mechanism_desc, i_clamp, gap_junction_site, threshold_detector>;

// High-level abstract representation of a cell.
class cable_cell {
public:
    using index_type = cell_lid_type;
    using size_type = cell_local_size_type;
    using value_type = double;
    using point_type = point<value_type>;

    using gap_junction_instance = mlocation;

    cable_cell_parameter_set default_parameters;

    /// Default constructor
    cable_cell();

    /// Copy constructor
    cable_cell(const cable_cell& other);

    /// Move constructor
    cable_cell(cable_cell&& other) = default;

    /// construct from morphology
    cable_cell(const class morphology& m, const label_dict& dictionary={});

    /// Access to morphology and embedding
    const concrete_embedding& embedding() const;
    const arb::morphology& morphology() const;
    const mprovider& provider() const;

    // Set cell-wide default physical and ion parameters.

    void set_default(init_membrane_potential prop) {
        default_parameters.init_membrane_potential = prop.value;
    }

    void set_default(axial_resistivity prop) {
        default_parameters.axial_resistivity = prop.value;
    }

    void set_default(temperature_K prop) {
        default_parameters.temperature_K = prop.value;
    }

    void set_default(membrane_capacitance prop) {
        default_parameters.membrane_capacitance = prop.value;
    }

    void set_default(initial_ion_data prop) {
        default_parameters.ion_data[prop.ion] = prop.initial;
    }

    void set_default(ion_reversal_potential_method prop) {
        default_parameters.reversal_potential_method[prop.ion] = prop.method;
    }

    // Painters and placers.
    //
    // Used to describe regions and locations where density channels, stimuli,
    // synapses, gap juncitons and detectors are located.

    // Density channels.
    void paint(const region&, mechanism_desc);

    // Properties.
    void paint(const region&, init_membrane_potential);
    void paint(const region&, axial_resistivity);
    void paint(const region&, temperature_K);
    void paint(const region&, membrane_capacitance);
    void paint(const region&, initial_ion_data);

    // Synapses.
    lid_range place(const locset&, mechanism_desc);

    // Stimuli.
    lid_range place(const locset&, i_clamp);

    // Gap junctions.
    lid_range place(const locset&, gap_junction_site);

    // Spike detectors.
    lid_range place(const locset&, threshold_detector);

    // Convenience access to placed items.

    const std::unordered_map<std::string, mlocation_map<mechanism_desc>>& synapses() const {
        return location_assignments().get<mechanism_desc>();
    }

    const mlocation_map<gap_junction_site>& gap_junction_sites() const {
        return location_assignments().get<gap_junction_site>();
    }

    const mlocation_map<threshold_detector>& detectors() const {
        return location_assignments().get<threshold_detector>();
    }

    const mlocation_map<i_clamp>& stimuli() const {
        return location_assignments().get<i_clamp>();
    }

    // Generic access to painted and placed items.
    const cable_cell_region_map& region_assignments() const;
    const cable_cell_location_map& location_assignments() const;

private:
    std::unique_ptr<cable_cell_impl, void (*)(cable_cell_impl*)> impl_;
};

} // namespace arb
