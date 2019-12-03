#include <string>

#include <arbor/morph/label_dict.hpp>
#include <arbor/morph/locset.hpp>
#include <arbor/morph/mprovider.hpp>
#include <arbor/morph/primitives.hpp>
#include <arbor/morph/region.hpp>

namespace arb {

void mprovider::init() {
    if (!label_dict_ptr) return;

    for (const auto& pair: label_dict_ptr->region()) {
        (void)region(pair.first);
    }

    for (const auto& pair: label_dict_ptr->region()) {
        (void)locset(pair.first);
    }

    label_dict_ptr = nullptr;
}

template <typename RegOrLocMap, typename Err>
static const auto& try_lookup(const std::string& name, RegOrLocMap& map, const label_dict* dict_ptr, Err errval) {
    auto it = map.find(name);
    if (it==map.end()) {
        if (dict_ptr) {
            map.emplace(name, errval);

            auto it = dict_ptr->find(name);
            if (it==dict_ptr.end()) {
                throw unbound_name(name);
            }

            return *map.emplace(name, thingify(*it))->first;
        }
        else {
            throw unbound_name(name);
        }
    }
    else if (!*it) {
        throw circular_definition(name);
    }
    else {
        return get<0>(*it);
    }
}

const mcable_list& mprovider::region(const std::string& name) {
    return try_lookup(name, regions_, label_dict_ptr, circular_def{});
}

const mcable_list& mprovider::locset(const std::string& name) {
    return try_lookup(name, locsets_, label_dict_ptr, circular_def{});
}


} // namespace arb
