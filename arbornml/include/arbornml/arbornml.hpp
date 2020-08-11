#pragma once

#include <memory>
#include <string>
#include <vector>

#include <arbor/morph/label_dict.hpp>
#include <arbor/morph/morphology.hpp>
#include <arbor/morph/stitch.hpp>

namespace arbnml {

// Collect and parse morphology elements from XML.
// No validation is performed against the NeuroML v2 schema.

struct morphology_data {
    // Cell id, or empty if morphology was taken from a top-level <morphology> element.
    std::string cell_id;

    // Morphology id.
    std::string id;

    // Morphology constructed from a signle NeuroML <morphology> element.
    arb::morphology morphology;

    // One region expression for each segment id.
    arb::label_dict segments;

    // One region expression for each name applied to one or more segments.
    arb::label_dict named_segments;

    // One region expression for each segmentGroup id.
    arb::label_dict groups;
};

// Represent NeuroML data determined by provided string.

struct neuroml_impl;

struct neuroml {
    neuroml();
    explicit neuroml(std::string nml_document);

    neuroml(neuroml&&);
    neuroml(const neuroml&) = delete;

    neuroml& operator=(neuroml&&);
    neuroml& operator=(const neuroml&) = delete;

    // Query top-level cells and (standalone) morphologies.

    std::vector<std::string> cell_ids() const;
    std::vector<std::string> morphology_ids() const;

    // Parse and retrieve top-level morphology or morphology associated with a cell.
    // Return nullopt if not found.

    arb::util::optional<morphology_data> morphology(const std::string& morph_id) const;
    arb::util::optional<morphology_data> cell_morphology(const std::string& cell_id) const;

    ~neuroml();
private:
    std::unique_ptr<neuroml_impl> impl_;
};

} // namespace arbnml
