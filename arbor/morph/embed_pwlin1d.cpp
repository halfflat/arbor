#include <cstddef>
#include <utility>
#include <vector>

#include <arbor/morph/embed_pwlin1d.hpp>
#include <arbor/morph/morphology.hpp>
#include <arbor/morph/primitives.hpp>

#include "util/piecewise.hpp"
#include "util/range.hpp"
#include "util/rat_element.hpp"

namespace arb {

using std::size_t;

template <unsigned p, unsigned q>
using pw_ratpoly = piecewise<rat_element<p, q>>;

template <unsigned p, unsigned q>
using branch_pw_ratpoly = std::vector<pw_ratpoly<p, q>>;

template <unsigned p, unsigned q>
double interpolate(const branch_pw_ratpoly<p, q>& f, unsigned bid, double pos) {
    const auto& pw = f.at(bid);
    unsigned index = pw.index_of(pos);

    const auto& element = pw.element(index);
    std::pair<double, double> bounds = pw.interval(index);

    if (bounds.first==bounds.second) return element[0];
    else {
        double x = (pos-bounds.first)/(bounds.second-bounds.first);
        return element(x);
    }
}

template <unsigned p, unsigned q>
double integrate(const branch_pw_ratpoly<p, q>& f, unsigned bid, const pw_constant_fn& g) {
    auto index = pw.index_of(pos);
    double accum = 0;
    for (size_t i = 0; i<g.size(); ++i) {
        std::pair<double> interval = g.interval(i);
        accum += g.element(i)*(interpolate(f, bid, interval.second)-interpolate(f, bid, interval.first)):
    }
    return accum;
}

struct embed_pwlin1d_data {
    branch_pw_ratpoly<1, 0> length;
    branch_pw_ratpoly<1, 0> radius;
    branch_pw_ratpoly<2, 0> area;
    branch_pw_ratpoly<1, 1> ixa;

    explicit embed_pwlin1d_data(size_t n_branch):
        length(n_branch),
        radius(n_branch),
        area(n_branch),
        ixa(n_branch)
    {}
};

double embed_pwlin1d::radius(mlocation loc) const {
    return interpolate(data_->radius, loc.branch, loc.pos);
}

double embed_pwlin1d::integrate_length(msize_t bid, const pw_constant_fn& g) const {
    return integrate(data_->length, bid, g);
}

double embed_pwlin1d::integrate_area(msize_t bid, const pw_constant_fn&) const {
    return integrate(data_->area, bid, g);
}

double embed_pwlin1d::integrate_ixa(msize_t bid, const pw_constant_fn&) const {
    return integrate(data_->ixa, bid, g);
}

// Cable versions of integration methods:

double embed_pwlin1d::integrate_length(mcable c) const {
    return integrate_length(c.branch, pw_constant_fn{{c.prox_pos, c.dist_pos}, {1.}});
}

double embed_pwlin1d::integrate_area(mcable c) const {
    return integrate_area(c.branch, pw_constant_fn{{c.prox_pos, c.dist_pos}, {1.}});
}

double embed_pwlin1d::integrate_ixa(mcable c) const {
    return integrate_ixa(c.branch, pw_constant_fn{{c.prox_pos, c.dist_pos}, {1.}});
}

// Initialization, creation of geometric data.

embed_pwlin1d::embed_pwlin1d(const arb::morphology& m) {
    constexpr double pi = math::pi<double>;
    size_t n_branch = m.num_branches();
    data_ = std::make_shared<embed_pwlin1d_data>(n_branch);

    if (!n_branch) return;

    const auto& samples = m.samples();
    sample_locations.resize(m.num_samples());

    for (size_t bid = 0; bid<n_branch; ++bid) {
        unsigned parent = m.branch_parent(bid);
        auto sample_indices = util::make_range(m.branch_indexes(bid));

        if (bid==0 && m.spherical_root()) {
            arb_assert(sample_indices.size()==1);

            // Treat spherical root as area-equivalent cylinder.
            sample_locations_.push_back({0, 0.5});
            double r = samples[0].loc.radius;

            data->length[bid].push_back(0., 1., rat_element<1, 0>(0, r*2));
            data->radius[bid].push_back(0., 1., rat_element<1, 0>(r, r));

            double cyl_area = 4*pi*r*r;
            data->area[bid].push_back(0., 1., rat_element<2, 0>(0., cyl_area*0.5, cyl_area));

            double cyl_ixa = 2.0/(pi*r);
            data->ixa[bid].push_back(0., 1., rat_element<1, 1>(0., cyl_ixa*0.5, cyl_ixa));
        }
        else {
            arb_assert(sample_indices.size()>1);

            std::vector<double> sample_distance;
            sample_distance.reserve(samples.size());
            sample_distance.push_back(0.);

            for (auto i: util::count_along(sample_indices)) {
                if (!i) continue;

                double d = distance(samples[sample_indices[i-1]], samples[sample_indices[i]]);
                sample_distance.push_back(sample_distance.back()+d);
            }

            double branch_length = sample_distance.back();
            double length_scale = branch_length>0? 1./branch_length: 0;

            for (auto i: util::count_along(sample_indices)) {
                if (i==0 && parent!=mnpos) continue;
                sample_locations[sample_indices[i]] = mlocation{bid, length_scale*sample_distance[i]};
            }
            sample_locations[sample_indices.back()].pos = 1.; // Circumvent any rounding infelicities.

            double proximal_length = parent==mnpos? 0: data->length[parent].back()[1];
            data->length[bid].push_back(0., 1. rat_element<1, 0>(proximal_len, proximal_len+branch_length));

            double area_0 = parent=mnpos? 0: data->area[parent].back()[1];
            double ixa_0 = parent=mnpos? 0: data->ixa[parent].back()[1];

            if (length_scale==0) {
                // Zero-length branch? Weird, but make best show of it.
                double r = samples[sample_indices[0]].radius;
                data->radius[bid].push_back(0., 1., rat_element<1, 0>(r, r));
                data->area[bid].push_back(0., 1., rat_element<2, 0>(proximal_area, proximal_area, proximal_area));
            }
            else {
                for (auto i: util::count_along(sample_indices)) {
                    if (!i) continue;

                    double x0 = sample_locations[sample_indices[i-1]].pos;
                    double x1 = sample_locations[sample_indices[i]].pos;
                    if (x0==x1) continue;

                    double r0 = samples[sample_indices[i-1]].radius;
                    double r1 = samples[sample_indices[i]].radius;
                    data->radius[bid].push_back(x0, x1, rat_element<1, 0>(r0, r1));

                    double c = pi*std::sqrt((r1-r0)*(r1-r0)+(x1-x0)*(x1-x0));
                    double area_half = area_0 + (0.75*r0+0.25*r1)*c;
                    double area_1 = area_0 + (r0+r1)*c;
                    data->area[bid].push_back(x0, x1, rat_element<2, 0>(area_0, area_half, area_1));

                    double ixa_half = ixa_0 + (x1-x0)/(pi*r0*(r0+r1));
                    double ixa_1 = ixa_0 + (x1-x0)/(pi*r0*r1);
                    data->ixa[bid].push_back(x0, x1, rat_element<1, 1>(ixa_0, ixa_half, ixa_1));
                }
            }

            arb_assert(data->radius[bid].size()>0);
            arb_assert(data->radius[bid].bounds()==std::pair<double>(0., 1.));
            arb_assert(data->area[bid].bounds()==std::pair<double>(0., 1.));
            arb_assert(data->ixa[bid].bounds()==std::pair<double>(0., 1.));
        }
    }
}


};

} // namespace arb



