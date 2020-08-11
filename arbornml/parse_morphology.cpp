#include <stack>
#include <string>
#include <vector>

#include <arbor/assert.hpp>
#include <arbor/morph/primitives.hpp>
#include <arbor/morph/stitch.hpp>
#include <arbor/util/optional.hpp>

#include <arbornml/arbornml.hpp>
#include <arbornml/nmlexcept.hpp>

#include "parse_morphology.hpp"
#include "xmlwrap.hpp"

using ::arb::util::optional;
using ::arb::util::nullopt;
using ::arb::util::just;

namespace arbnml {

struct neuroml_segment {
    // Morhpological data:
    non_negative id;
    std::string name;
    optional<arb::mpoint> proximal;
    arb::mpoint distal;
    optional<non_negative> parent_id;
    double along = 1;

    // Data for error reporting:
    unsigned line = 0;

    // Data for topological sort:
    enum { cyclic = -2, unknown = -1 };
    std::ptrdiff_t tdepth = unknown;
};

static arb::stitched_morphology construct_morphology(std::vector<neuroml_segment>& segs) {
    arb::stitch_builder builder;
    if (segs.empty()) return arb::stitched_morphology{builder};

    // Sort segments topologically with respect to parent relationship.

    std::unordered_map<non_negative, non_negative> id_to_index;
    for (std::size_t i = 0; i<segs.size(); ++i) {
        auto iter_success = id_to_index.insert({segs[i].id, (non_negative)i});
        if (!iter_success.second) {
            throw bad_segment(segs[i].id, segs[i].line); // Duplicate id.
        }
    }

    unsigned n_root = 0;
    std::stack<std::size_t> trace;
    for (std::size_t i = 0; i<segs.size(); ++i) {
        std::size_t j = i;
        for (;;) {
            if (segs[j].tdepth>=0) break;

            if (!segs[j].parent_id) {
                if (++n_root>1) {
                    // Duplicate root.
                    throw bad_segment(segs[j].id, segs[i].line);
                }
                segs[j].tdepth = 0;
                break;
            }

            if (segs[j].tdepth==neuroml_segment::cyclic) {
                throw parse_error("cyclic dependency in segments");
            }

            segs[j].tdepth = neuroml_segment::cyclic;
            trace.push(j);

            auto next = id_to_index.find(segs[j].parent_id.value());
            j = next!=id_to_index.end()? next->second: throw bad_segment(segs[j].id, segs[i].line); // Missing parent.
        }

        int d = segs[j].tdepth;
        while (!trace.empty()) {
            segs[trace.top()].tdepth = ++d;
            trace.pop();
        }
        arb_assert(segs[i].tdepth >= 0);
    }

    arb_assert(segs.empty() || n_root==1);
    std::sort(segs.begin(), segs.end(), [](auto& a, auto& b) { return a.tdepth<b.tdepth; });

    // Construct result from sorted segments.

    for (auto& s: segs) {
        arb::mstitch stitch(std::to_string(s.id), s.distal);
        stitch.prox = s.proximal;

        if (s.parent_id) {
            builder.add(stitch, std::to_string(s.parent_id.value()), s.along);
        }
        else {
            builder.add(stitch);
        }
    }

    return arb::stitched_morphology(std::move(builder));
}

morphology_data parse_morphology_element(xml_xpathctx ctx, xml_node morph) {
    morphology_data M;

    M.id = morph.prop<std::string>("id", std::string{});

    // NeuroML schema specifies that segment ids must be non-negative integers
    // of arbitrary magnitude, but we attempt to parse them as unsigned long long
    // values, and will throw an exception if they do not fit in the range.
    std::vector<neuroml_segment> segments;

    // TODO: precompile xpath queries for nml:distal, nml:proximal, nml:parent.
    const char* q_parent = "./nml:parent";
    const char* q_proximal = "./nml:proximal";
    const char* q_distal = "./nml:distal";

    for (auto n: ctx.query(morph, "./nml:segment")) {
        neuroml_segment seg;

        try {
            seg.id = -1;
            seg.id = n.prop<non_negative>("id");
            std::string name = n.prop<std::string>("name", std::string{});

            auto result = ctx.query(n, q_parent);
            if (!result.empty()) {
                seg.parent_id = result[0].prop<non_negative>("segment");
                seg.along = result[0].prop<double>("fractionAlong", 1.0);
            }

            result = ctx.query(n, q_proximal);
            if (!result.empty()) {
                double x = result[0].prop<double>("x");
                double y = result[0].prop<double>("y");
                double z = result[0].prop<double>("z");
                double diameter = result[0].prop<double>("diameter");
                if (diameter<0) throw bad_segment(seg.id, n.line());

                seg.proximal = arb::mpoint{x, y, z, diameter/2};
            }

            if (!seg.parent_id && !seg.proximal) throw bad_segment(seg.id, n.line());

            result = ctx.query(n, q_distal);
            if (!result.empty()) {
                double x = result[0].prop<double>("x");
                double y = result[0].prop<double>("y");
                double z = result[0].prop<double>("z");
                double diameter = result[0].prop<double>("diameter");
                if (diameter<0) throw bad_segment(seg.id, n.line());

                seg.distal = arb::mpoint{x, y, z, diameter/2};
            }
            else {
                throw bad_segment(seg.id, n.line());
            }
        }
        catch (parse_error& e) {
            throw bad_segment(seg.id, n.line());
        }

        seg.line = n.line();
        segments.push_back(std::move(seg));
    }

    if (segments.empty()) return M;

    arb::stitched_morphology stitched = construct_morphology(segments);
    M.morphology = stitched.morphology();
    M.segments = stitched.labels();

    // Collate 'name' attributes for segments: associate each name
    // with the corresponding segments in a region expression.

    std::unordered_multimap<std::string, non_negative> name_to_ids;
    std::unordered_set<std::string> names;

    for (auto& s: segments) {
        if (!s.name.empty()) {
            name_to_ids.insert({s.name, s.id});
            names.insert(s.name);
        }
    }

    for (const auto& n: names) {
        arb::region r;
        auto ids = name_to_ids.equal_range(n);
        for (auto i = ids.first; i!=ids.second; ++i) {
            r = join(std::move(r), M.segments.regions().at(std::to_string(i->second)));
        }
        M.named_segments.set(n, std::move(r));
    }

    // TODO: segmentGroups.
#if 0
    stuct group_info {
        std::vector<non_negative> segments;
        std::vector<std::string> includes;
        // paths
        // subtree
        
    };
#endif
    return M;
}

} // namespace arbnml
