#include <numeric>
#include <stack>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include <arbor/assert.hpp>
#include <arbor/morph/locset.hpp>
#include <arbor/morph/primitives.hpp>
#include <arbor/morph/region.hpp>
#include <arbor/morph/stitch.hpp>
#include <arbor/util/optional.hpp>
#include <arbor/util/variant.hpp>

#include <arbornml/arbornml.hpp>
#include <arbornml/nmlexcept.hpp>

#include "parse_morphology.hpp"
#include "xmlwrap.hpp"

using ::arb::util::optional;
using ::arb::util::nullopt;
using ::arb::util::just;
using ::arb::util::variant;
using ::arb::region;

namespace arbnml {

// Q&D error return class.

template <typename V, typename Err>
struct either: ::arb::util::variant<V, Err> {
    template <typename X, typename = std::enable_if_t<!std::is_same<either, std::decay_t<X>>::value>>
    either(X&& x): ::arb::util::variant<V, Err>(std::forward<X>(x)) {}

    explicit operator bool() const { return ::arb::util::get_if<0>(*this); }

    const V& get() const { return ::arb::util::get<0>(*this); }
    V& get() { return ::arb::util::get<0>(*this); }

    const Err& error() const { return get<1>(*this); }
    Err& error() { return ::arb::util::get<1>(*this); }
};

// Box is a container of size 0 or 1.

template <typename X>
struct box {
    optional<X> x;

    X* begin() { return x? &(*x): nullptr; }
    X* end() { return x? &(*x)+1: nullptr; }

    const X* begin() const { return x? &(*x): nullptr; }
    const X* end() const { return x? &(*x)+1: nullptr; }

    box() = default;
    box(const X& x): x(x) {}

    std::size_t size() { return !!x; }
};

// Return vector of depths; sorting object collection by depth will
// give a topological order.
//
// The functional Inset takes a reference to an object and returns
// a range or collection of indices to that object's precedessors.
//
// If a cycle is encountered, return detected_cycle{i} where i
// is the index of an item in the cycle.

struct cycle_detected { std::size_t index; };

template <typename X, typename Inset>
either<std::vector<std::size_t>, cycle_detected> topological_sort(const std::vector<X>& objects, Inset inset) {
    using std::begin;
    using std::end;

    constexpr std::size_t unknown = -1;
    constexpr std::size_t cycle = -2;

    std::size_t n = objects.size();
    std::vector<std::size_t> depth(n, unknown);
    std::stack<std::size_t> stack;

    for (std::size_t i = 0; i<n; ++i) {
        if (depth[i]!=unknown) continue;

        depth[i] = cycle;
        stack.push(i);

        while (!stack.empty()) {
            std::size_t j = stack.top();
            std::size_t d = 0;
            bool resolve = true;

            auto in = inset(objects[j]);
            for (auto k = begin(in); k!=end(in); ++k) {
                switch (depth[*k]) {
                case cycle:
                    return cycle_detected{*k};
                case unknown:
                    depth[*k] = cycle;
                    stack.push(*k);
                    resolve = false;
                    break;
                default:
                    d = std::max(d, 1+depth[*k]);
                }
            }

            if (resolve) {
                depth[j] = d;
                stack.pop();
            }
        }
    }

    return depth;
}


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

    // Topological depth:
    std::size_t tdepth = 0;
};

struct group_info {
    std::string id;
    std::vector<non_negative> segments;
    std::vector<std::string> includes;
    std::vector<std::pair<non_negative, non_negative>> paths;
    std::vector<non_negative> subtrees;

    // Data for error reporting:
    unsigned line = 0;

    // Topological depth:
    std::size_t tdepth = 0;
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

    for (const auto& s: segs) {
        if (s.parent_id && !id_to_index.count(*s.parent_id)) {
            throw bad_segment(s.id, s.line); // No such parent id.
        }
    }

    auto ts_result = topological_sort(segs,
        [&id_to_index](const neuroml_segment& s) {
            return s.parent_id? box<non_negative>{id_to_index.at(*s.parent_id)}: box<non_negative>{}; });

    if (!ts_result) {
        auto& seg = segs[ts_result.error().index];
        throw cyclic_dependency(std::to_string(seg.id), seg.line);
    }
    const std::vector<std::size_t>& depth = ts_result.get();

    for (std::size_t i = 0; i<segs.size(); ++i) {
        segs[i].tdepth = depth[i];
    }
    std::sort(segs.begin(), segs.end(), [](auto& a, auto& b) { return a.tdepth<b.tdepth; });

    // Check for multiple roots:
    if (segs.size()>1 && segs[1].tdepth==0) throw bad_segment(segs[1].id, segs[1].line);

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

void build_segment_groups(morphology_data& M, std::vector<group_info>& groups) {
    std::vector<region> group_region(groups.size());

    // Sort groups topologically with respect to include relationship.

    std::unordered_map<std::string, std::size_t> id_to_index;
    std::vector<std::vector<std::size_t>> includes_by_index(groups.size());

    for (std::size_t i = 0; i<groups.size(); ++i) {
        auto iter_success = id_to_index.insert({groups[i].id, i});
        if (!iter_success.second) {
            throw bad_segment_group(groups[i].id, groups[i].line); // Duplicate id.
        }
    }
    for (std::size_t i = 0; i<groups.size(); ++i) {
        for (const std::string& include: groups[i].includes) {
            auto iter = id_to_index.find(include);
            if (iter==id_to_index.end()) throw bad_segment_group(include, groups[i].line);
            includes_by_index[i].push_back(iter->second);
        }
    }

    auto ts_result = topological_sort(groups,
        [&id_to_index, includes_by_index](const group_info& g) {
            return includes_by_index.at(id_to_index.at(g.id)); });

    if (!ts_result) {
        auto& group = groups[ts_result.error().index];
        throw cyclic_dependency(group.id, group.line);
    }

    const std::vector<std::size_t>& depth = ts_result.get();

    std::vector<std::size_t> topo_order(groups.size());
    std::iota(topo_order.begin(), topo_order.end(), std::size_t(0));
    std::sort(topo_order.begin(), topo_order.end(), [&depth](auto& a, auto& b) { return depth[a]<depth[b]; });

    for (auto index: topo_order) {
        const auto& group = groups[index];
        region r;

        for (auto seg_id: group.segments) {
            auto opt_reg = M.segments.region(std::to_string(seg_id));
            if (!opt_reg) throw bad_segment_group(group.id, group.line);
            r = join(std::move(r), *opt_reg);
        }

        for (auto path: group.paths) {
            auto opt_from = M.segments.region(std::to_string(path.first));
            if (!opt_from) throw bad_segment_group(group.id, group.line);
            auto opt_to = M.segments.region(std::to_string(path.second));
            if (!opt_to) throw bad_segment_group(group.id, group.line);

            r = join(std::move(r), arb::reg::between(arb::ls::most_proximal(*opt_from), arb::ls::most_distal(*opt_to)));
        }

        for (auto seg: group.subtrees) {
            auto opt_from = M.segments.region(std::to_string(seg));
            if (!opt_from) throw bad_segment_group(group.id, group.line);

            r = join(std::move(r), arb::reg::distal_interval(arb::ls::most_distal(*opt_from), INFINITY));
        }

        for (auto j: includes_by_index.at(index)) {
            r = join(std::move(r), group_region.at(j));
        }

        group_region[index] = r;
        M.groups.set(group.id, std::move(r));
    }
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
        int line = n.line(); // for error context!

        try {
            seg.id = -1;
            seg.id = n.prop<non_negative>("id");
            std::string name = n.prop<std::string>("name", std::string{});

            auto result = ctx.query(n, q_parent);
            if (!result.empty()) {
                line = result[0].line();
                seg.parent_id = result[0].prop<non_negative>("segment");
                seg.along = result[0].prop<double>("fractionAlong", 1.0);
            }

            result = ctx.query(n, q_proximal);
            if (!result.empty()) {
                line = result[0].line();
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
                line = result[0].line();
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
            throw bad_segment(seg.id, line);
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

    // TODO: precompile xpath queries for following:
    const char* q_member = "./nml:member";
    const char* q_include = "./nml:include";
    const char* q_path = "./nml:path";
    const char* q_from = "./nml:from";
    const char* q_to = "./nml:to";
    const char* q_subtree = "./nml:subTree";

    std::vector<group_info> groups;

    for (auto n: ctx.query(morph, "./nml:segmentGroup")) {
        group_info group;
        int line = n.line(); // for error context!

        try {
            group.id = n.prop<std::string>("id");
            for (auto elem: ctx.query(n, q_member)) {
                line = elem.line();
                group.segments.push_back(elem.prop<non_negative>("segment"));
            }
            for (auto elem: ctx.query(n, q_include)) {
                line = elem.line();
                group.includes.push_back(elem.prop<std::string>("segmentGroup"));
            }
            for (auto elem: ctx.query(n, q_path)) {
                line = elem.line();
                auto froms = ctx.query(elem, q_from);
                auto tos = ctx.query(elem, q_to);
                // Schema says its okay to have zero 'from' or 'to' elements in
                // a path, but what does that even mean?
                if (froms.empty() || tos.empty()) {
                    throw bad_segment_group(group.id, line);
                }

                line = froms[0].line();
                non_negative seg_from = froms[0].prop<non_negative>("segment");

                line = tos[0].line();
                non_negative seg_to = tos[0].prop<non_negative>("segment");

                group.paths.push_back({seg_from, seg_to});
            }
            for (auto elem: ctx.query(n, q_subtree)) {
                // Schema says we can have a subTree _to_ a segment instead.
                // It is not documented what this would mean, so we're going
                // to ignore that option for now.
                line = elem.line();
                auto froms = ctx.query(elem, q_from);
                if (froms.empty()) throw bad_segment_group(group.id, line);

                line = froms[0].line();
                non_negative seg_from = froms[0].prop<non_negative>("segment");
                group.subtrees.push_back(seg_from);
            }
        }
        catch (parse_error& e) {
            throw bad_segment_group(group.id, line);
        }

        group.line = n.line();
        groups.push_back(std::move(group));
    }

    build_segment_groups(M, groups);
    return M;
}

} // namespace arbnml
