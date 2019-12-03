#pragma once

#include <string>
#include <unordered_map>

#include <arbor/morph/embed_pwlin1d.hpp>
#include <arbor/morph/primitives.hpp>
#include <arbor/morph/label_dict.hpp>
#include <arbor/util/either.hpp>

namespace arb {

using concrete_embedding = embed_pwlin1d;

struct mprovider: public  {
    morphology morpho;
    concrete_embedding embed;

    // Throw exception on missing or recursive definition.
    const mcable_list& region(const std::string& name);
    const mlocation_list& locset(const std::string& name);

    mprovider(morphology m, const label_dict& dict):
        morpho(m), embed(m), label_dict_ptr(&dict)
    {
        init();
    }

    explicit mprovider(morphology m):
        morpho(m), embed(m), label_dict_ptr(nullptr)
    {
        init();
    }

private:
    struct circular_def {};

    // Maps are mutated only during initialization phase of mprovider.
    mutable std::unordered_map<std::string, util::either<mcable_list, circular_def>> regions_;
    mutable std::unordered_map<std::string, util::either<mlocation_list, circular_def>> locsets_;

    // Non-null only during initialization phase.
    mutable const label_dict* label_dict_ptr;

    // Perform greedy initialization of concrete region, locset maps.
    void init();
};

} // namespace arb
