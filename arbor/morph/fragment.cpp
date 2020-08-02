#include <memory>
#include <unordered_map>

#include <arbor/morph/fragment.hpp>
#include <arbor/morph/morphexcept.hpp>
#include <arbor/util/optional.hpp>

#include "util/ordered_forest.hpp"
#include "util/maputil.hpp"

namespace arb {

struct fragment_builder_impl {
    struct fragment_segment {
        double along_prox;
        double along_dist;

        mpoint prox;
        mpoint dist;
        int tag;
        msize_t seg_id;
    };

    using forest_type = util::ordered_forest<fragment_segment>;

    forest_type forest;
    std::unordered_map<std::string, forest_type::iterator> id_to_node;
    std::string last_id;

    void add(mfragment f, const std::string& parent, double along, bool infer_prox = false) {
        if (id_to_node.count(f.id)) throw duplicate_fragment_id(f.id);

        forest_type::iterator p;

        if (!(parent.empty() && forest.empty())) {
            p = find_fragment_along(parent, along);
            arb_assert(p);

            if (along==p->along_prox) {
                if (!f.prox) f.prox = p->prox;
                p = p.parent();
            }
            else if (along<p->along_dist) {
                // Split parent node p at along.
                auto split = *p;

                mpoint point = lerp(p->prox, p->dist, (along-p->along_prox)/(p->along_dist-p->along_prox));
                if (!f.prox) f.prox = point;

                p->dist = point;
                p->along_dist = along;
                split.prox = point;
                split.along_prox = along;

                auto i = forest.push_child(p, split);
                while (i.next()) forest.graft_child(i, forest.prune_after(i));
            }
            else {
                if (!f.prox) f.prox = p->dist;
            }
        }
        if (!f.prox) throw missing_fragment_start(f.id);

        fragment_segment n{0., 1., f.prox.value(), f.dist, f.tag, msize_t(-1)};
        id_to_node[f.id] = p? forest.push_child(p, n): forest.push_front(n);
        last_id = f.id;
    }

    forest_type::iterator find_fragment_along(const std::string& id, double along) {
        if (along<0 || along>1) throw invalid_fragment_position(id, along);

        auto map_it = id_to_node.find(id);
        if (map_it==id_to_node.end()) throw no_such_fragment(id);

        auto i = map_it->second;
        arb_assert(i->along_prox==0);
        arb_assert(i->along_dist==1 || i.child());

        while (along>i->along_dist) {
            // Continuation is last child.
            i = i.child();
            arb_assert(i);
            while (i.next()) i = i.next();
        }
        return i;
    }
};

fragment_builder::fragment_builder(): impl_(new fragment_builder_impl) {}

fragment_builder& fragment_builder::add(mfragment f, const std::string& parent_id, double along) {
    impl_->add(std::move(f), parent_id, along);
    return *this;
}

fragment_builder& fragment_builder::add(mfragment f, double along) {
    if (impl_->forest.empty()) throw no_such_fragment{""};
    return add(std::move(f), impl_->last_id, along);
}

struct fragment_tree_impl {
    std::unordered_multimap<std::string, msize_t> id_to_segs;
    segment_tree stree;

    fragment_tree_impl(fragment_builder_impl bimpl) {
        auto iter = bimpl.forest.preorder_begin();
        auto end = bimpl.forest.preorder_end();

        while (iter!=end) {
            msize_t seg_parent_id = iter.parent()? iter.parent()->seg_id: mnpos;
            iter->seg_id = stree.append(seg_parent_id, iter->prox, iter->dist, iter->tag);
        }

        for (const auto& id_node: bimpl.id_to_node) {
            const std::string& id = id_node.first;
            auto iter = id_node.second;

            while (iter) {
                id_to_segs.insert({id, iter->seg_id});
                iter = iter.child();
                while (iter.next()) {
                    iter = iter.next();
                }
            }
        }
    }
};

fragment_tree::fragment_tree(fragment_builder&& builder):
    impl_(new fragment_tree_impl(std::move(*builder.impl_)))
{}

fragment_tree::fragment_tree(const fragment_builder& builder):
    impl_(new fragment_tree_impl(*builder.impl_))
{}

arb::morphology fragment_tree::morphology() const {
    return arb::morphology(impl_->stree);
}

} // namespace arb
