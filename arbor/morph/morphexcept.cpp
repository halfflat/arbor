#include <string>
#include <sstream>

#include <arbor/morph/primitives.hpp>
#include <arbor/morph/morphexcept.hpp>

#include "util/strprintf.hpp"

namespace arb {

using arb::util::pprintf;

invalid_mlocation::invalid_mlocation(mlocation loc):
    morphology_error(pprintf("invalid mlocation {}", loc)),
    loc(loc)
{}

no_such_branch::no_such_branch(msize_t bid):
    morphology_error(pprintf("no such branch id {}", branch)),
    bid(bid)
{}

invalid_mcable::invalid_mcable(mcable cable):
    morphology_error(pprintf("invalid mcable {}", cable)),
    cable(cable)
{}

invalid_sample_parent::invalid_sample_parent(msize_t parent, msize_t tree_size):
    morphology_error(pprintf("invalid sample parent {} for a sample tree of size {}", parent, tree_size)),
    parent(parent),
    tree_size(tree_size)
{}

label_type_mismatch::label_type_mispatch(const std::string& label):
    morphology_error(pprintf("label \"{}\" is already bound to a different type of object", label)),
    label(label)
{}

incomplete_branch::incomplete_branch(msize_t bid):
    morphology_error(pprintf("insufficent samples to define branch id {}", bid)),
    bid(bid)
{}

} // namespace arb
