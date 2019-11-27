#pragma once

#include <arbor/morph/primitives.hpp>

namespace arb {

// A concrete embedding encodes the geometry of a morphology.

class concrete_embedding {
public:
    virtual mlocation sample_location(msize_t sid) const = 0;
    virtual double branch_length(msize_t bid) const = 0;

    // Integrating over geometry:
    // length; membrane area; inverse cross-sectional area.

    virtual double integrate_length(mcable) = 0;
    virtual double integrate_area(mcable) = 0;
    virtual double integrate_ixa(mcable) = 0;

    // Access underlying morphology.
    virtual const arb::morphology& morphology() = 0;
};

// An mprovider wraps a concrete embedding and also provides
// a lookup facility for named regions and locsets.

class mprovider: public concrete_embedding {
public:
    virtual mlocation_list named_locset(const std::string&) const = 0;
    virtual mcable_list named_region(const std::string&) const = 0;
};

// (Sketching stuff out loud ...)

struct pw1d_embedding: public concrete_embedding {
    explicit pw1d_embedding(const arb::morphology& m):
        morphology(m)
    {
        // init...
    }

    mlocation sample_location(msize_t sid) const override {
        // ...
    }

    double branch_length(msize_t bid) const override {
        // ...
    }

    double integrate_length(mcable) const override {
        // ...
    }

    double integrate_area(mcable) const override {
        // ...
    }

    double integrate_ixa(mcable) const override {
        // ...
    }
};

// Cached, lazily evaluated named expressions.

using region_map = std::unordered_map<std::string, mcable_list>;
using locset_map = std::unordered_map<std::string, mlocation_list>;

struct lazy_mprovider: public mprovider {
    lazy_mprovider(const concrete_embedding& embed, const label_dict& dict):
        embed(embed), dict(dict)
    {}

    // Forward geometry queries to embed:

    mlocation sample_location(msize_t sid) const override { return embed.sample_location(id); }
    double branch_length(msize_t bid) const override { return embed.branch_length(bid); }

    double integrate_length(mcable c) override { return embed.integrate_length(c); }
    double integrate_area(mcable c) override { return embed.integrate_area(c); }
    double integrate_ixa(mcable c) override { return embed.integrate_ixa(c); }

    const arb::morphology& morphology() override { return embed.xx; }
};

} // namespace arb
