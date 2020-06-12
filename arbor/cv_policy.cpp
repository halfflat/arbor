#include <vector>

#include <arbor/cable_cell.hpp>
#include <arbor/cv_policy.hpp>

#include "util/rangeutil.hpp"
#include "util/span.hpp"

// Discretization policy implementations:

namespace arb {

locset cv_policy_max_extent::cv_boundary_points(const cable_cell& cell) const {
    const unsigned nbranch = cell.morphology().num_branches();
    const auto& embed = cell.embedding();
    if (!nbranch || max_extent_<=0) return ls::nil();

    std::vector<mlocation> points;
    double oomax_extent = 1./max_extent_;
    auto comps = components(cell.morphology(), thingify(domain_, cell.provider()));

    for (auto& comp: comps) {
        for (mcable c: comp) {
            double cable_length = embed.integrate_length(c);
            unsigned ncv = std::ceil(cable_length*oomax_extent);
            double scale = (c.dist_pos-c.prox_pos)/ncv;

            if (flags_&cv_policy_flag::interior_forks) {
                for (unsigned i = 0; i<ncv; ++i) {
                    points.push_back({c.branch, c.prox_pos+(1+2*i)*scale/2});
                }
            }
            else {
                for (unsigned i = 0; i<ncv; ++i) {
                    points.push_back({c.branch, i*scale});
                }
                points.push_back({c.branch, c.dist_pos});
            }
        }
    }

    util::sort(points);
    return join(locset(std::move(points)), ls::boundary(domain_));
}

locset cv_policy_fixed_per_branch::cv_boundary_points(const cable_cell& cell) const {
    const unsigned nbranch = cell.morphology().num_branches();
    if (!nbranch) return ls::nil();

    std::vector<mlocation> points;
    double ooncv = 1./cv_per_branch_;
    auto comps = components(cell.morphology(), thingify(domain_, cell.provider()));

    for (auto& comp: comps) {
        for (mcable c: comp) {
            double scale = (c.dist_pos-c.prox_pos)*ooncv;

            if (flags_&cv_policy_flag::interior_forks) {
                for (unsigned i = 0; i<cv_per_branch_; ++i) {
                    points.push_back({c.branch, (1+2*i)*scale/2});
                }
            }
            else {
                for (unsigned i = 0; i<cv_per_branch_; ++i) {
                    points.push_back({c.branch, i*scale});
                }
                points.push_back({c.branch, c.dist_pos});
            }
        }
    }

    util::sort(points);
    return join(locset(std::move(points)), ls::boundary(domain_));
}

locset cv_policy_every_sample::cv_boundary_points(const cable_cell& cell) const {
    const unsigned nbranch = cell.morphology().num_branches();
    if (!nbranch) return ls::nil();

    // Always include branch proximal points, so that forks are trivial.
    // Ignores interior_forks flag.

    auto sample_indices = util::make_span(cell.morphology().num_samples());
    return
        join(
            ls::boundary(domain_),
            ls::restrict(
                join(ls::on_branches(0.),
                    std::accumulate(sample_indices.begin(), sample_indices.end(), ls::nil(),
                        [](auto&& l, auto&& r) { return sum(std::move(l), ls::sample(r)); })),
                domain_));
}

} // namespace arb
