#pragma once

#include <memory>
#include <string>

#include <arbor/morph/morphology.hpp>
#include <arbor/morph/primitives.hpp>
#include <arbor/morph/label_dict.hpp>
#include <arbor/morph/region.hpp>
#include <arbor/util/optional.hpp>

namespace arb {

// Stiches represent an alternative building block for morphologies.
//
// A stitch describes a portion of the morphology delimited by two
// `mpoint`s. Stitches can be attached to a parent stich at any
// point along the parent, interpolated linearly from the end points.
// Each stitch is associated with a unique string label, and optionally
// an integer tag value.
//
// The stitch builder collects stitches and produces the corresponding
// morphology and region/location labels.

struct mstitch {
    std::string id;
    util::optional<mpoint> prox;
    mpoint dist;
    int tag;

    mstitch(std::string id, mpoint prox, mpoint dist, int tag = 0):
        id(std::move(id)), prox(std::move(prox)), dist(std::move(dist)), tag(tag)
    {}

    mstitch(std::string id, mpoint dist, int tag = 0):
        id(std::move(id)), dist(std::move(dist)), tag(tag)
    {}
};

struct stitch_builder_impl;
struct stitch_tree;

struct stitch_builder {
    stitch_builder();

    stitch_builder(const stitch_builder&) = delete;
    stitch_builder(stitch_builder&&) = default;

    stitch_builder& operator=(const stitch_builder&) = delete;
    stitch_builder& operator=(stitch_builder&&) = default;

    // Make a new stitch in the morphology, return reference to self.
    //
    // If the stitch does not contained a proximal point, it will be
    // inferred from the point where it attaches to the parent stitch.
    // If the parent is omitted, it will be taken to be the last stitch
    // added.

    stitch_builder& add(mstitch f, const std::string& parent_id, double along = 1.);
    stitch_builder& add(mstitch f, double along = 1.);

private:
    friend stitch_tree;
    std::unique_ptr<stitch_builder_impl> impl_;
};

// From stitch builder construct morphology, region expressions.

struct stitch_tree_impl;

struct stitch_tree {
    stitch_tree() = delete;
    stitch_tree(const stitch_builder&); // implicit
    stitch_tree(stitch_builder&&); // implicit

    arb::morphology morphology() const;

    // Create labeled regions for each stitch with label equal to the stitch id, prepended by `prefix`.
    label_dict labels(const std::string& prefix="") const;

    region stitch(const std::string& id) const;

private:
    std::unique_ptr<stitch_tree_impl> impl_;
};

} // namesapce arb


