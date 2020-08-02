#pragma once

#include <memory>
#include <string>

#include <arbor/morph/morphology.hpp>
#include <arbor/morph/primitives.hpp>
#include <arbor/util/optional.hpp>

namespace arb {

// Fragments represent an alternative building block for morphologies.
//
// A fragment describes a portion of the morphology delimited by two
// `mpoint`s. Fragments can be attached to a parent fragment at any
// point along the parent, interpolated linearly from the end points.
// Fragments are associated with a unique string label, and optionally
// an integer tag value.
//
// A fragment tree comprises a tree of labelled fragments attached as
// above. The fragment tree can generate both the corresponding
// morphology and morphology expressions describing individual
// fragments or points within them.

struct mfragment {
    std::string id;
    util::optional<mpoint> prox;
    mpoint dist;
    int tag;

    mfragment(std::string id, mpoint prox, mpoint dist, int tag = 0):
        id(std::move(id)), prox(std::move(prox)), dist(std::move(dist)), tag(tag)
    {}

    mfragment(std::string id, mpoint dist, int tag = 0):
        id(std::move(id)), dist(std::move(dist)), tag(tag)
    {}
};

struct fragment_builder_impl;
struct fragment_tree_impl;
struct fragment_tree;

// Build a fragment tree incrementally.

struct fragment_tree_impl;

struct fragment_builder {
    fragment_builder();

    fragment_builder(const fragment_builder&) = delete;
    fragment_builder(fragment_builder&&) = default;

    fragment_builder& operator=(const fragment_builder&) = delete;
    fragment_builder& operator=(fragment_builder&&) = default;

    // Attach a fragment to the tree, returning reference to self.
    //
    // The fragment may be given explicitly, or as a distal point and id.
    // In this latter circumstance, the proximal point will be the attachment
    // point on the parent fragment.
    //
    // If no parent id is given, use the last added fragment as parent.
    //
    // Attaching a fragment with the same id as one already added will throw
    // arb::duplicate_fragment_id.

    fragment_builder& add(mfragment f, const std::string& parent_id, double along = 1.);
    fragment_builder& add(mfragment f, double along = 1.);

private:
    std::unique_ptr<fragment_builder_impl> impl_;
    friend fragment_tree;
};

// From fragment builder construct morphology, region expressions.

struct fragment_tree {
    fragment_tree() = delete;
    fragment_tree(const fragment_builder&); // implicit
    fragment_tree(fragment_builder&&); // implicit

    arb::morphology morphology() const;
    //arb::label_dict labels(const std::string& prefix="") const;

    //locset fragment_along(const std::string& label, double along) const;
    //region fragment(const std::string& label) const;

private:
    std::unique_ptr<fragment_tree_impl> impl_;
};


} // namesapce arb


