#include <cfloat>
#include <iostream>
#include <vector>

#include <arbor/simulation.hpp>
#include <arbor/load_balance.hpp>
#include <arbor/mc_cell.hpp>
#include <arbor/recipe.hpp>

using namespace arb;

// 1. stability with passive channel.

struct pas1comp: recipe {
    // One compartment parameter set; circuit behaviour as below.
    // Current source from time zero.
    // Initial condition v = erev.
    //
    //        cm
    // +------||------+
    // |              |
    // +-vVVV^---+|---+
    // |   rm    erev |
    // |              |
    // +-----(\)------+
    //       iinj

    struct param {
        double cm    = 0.01;   // [nF]
        double rm    = 100;    // [MΩ]
        double erev  = -65;    // [mV]
        double iinj  = 0.1;    // [nA]
    };

    explicit pas1comp(std::vector<param> param_set):
        param_set_(std::move(param_set))
    {}

    explicit pas1comp(param p):
        param_set_({p})
    {}

    cell_size_type num_cells() const override { return param_set_.size(); }
    cell_size_type num_targets(cell_gid_type) const override { return 0; }
    cell_size_type num_probes(cell_gid_type) const override { return 1; }
    cell_kind get_cell_kind(cell_gid_type) const override { return cell_kind::cable1d_neuron; }

    util::unique_any get_cell_description(cell_gid_type gid) const override {
        auto& p = param_set_.at(gid);
        mc_cell c;

        double r = 9e-6; // [m]
        double area = r*r*4*pi*1e-12; // [m²]
        mechanism_desc pas("pas");
        pas["g"] = 1e-10/(p.rm*area); // [S/cm^2]
        pas["e"] = p.erev;

        auto soma = c.add_soma(r*1e6);
        soma->cm = p.cm*1e-9/area;
        soma->add_mechanism(pas);

        c.add_stimulus({0,0}, {0, FLT_MAX, p.iinj});
        return c;
    }

    probe_info get_probe(cell_member_type probe_id) const override {
        return probe_info{probe_id, 0, cell_probe_address{{0, 0.}, cell_probe_address::membrane_voltage}};
    }

    util::any get_global_properties(cell_kind) const override {
        mc_cell_global_properties props;
        props.init_membrane_potential_mV = 77; //param_set_[0].erev; // pick first!
        return props;
    }

    static constexpr double pi = 3.141592653589793238462643383279502884;
    std::vector<param> param_set_;
};

struct pas1comp_result {
    pas1comp::param p;
    double dt;
    double t;
    double v;
    double v_err;
};

std::vector<pas1comp_result> run_pas1comp(pas1comp::param p, int nstep_min, int nstep_max) {
    std::vector<pas1comp_result> results;

    auto context = make_context();

    double tau = p.rm*p.cm; // [ms]
    double vinf = p.erev+p.iinj*p.rm; // [mV]

    for (int n = nstep_min; n<=nstep_max; ++n) {
        double dt = tau/n;
        double t_end = tau + dt;

        pas1comp recipe(p);
        simulation sim(recipe, partition_load_balance(recipe, context), context);


        double t = -1;
        double v = 0;
        auto sample_once = [&t,&v](cell_member_type, probe_tag, std::size_t n, const sample_record* rec) {
            const double* ptr = nullptr;
            if (!n || !(ptr = util::any_cast<const double*>(rec[0].data))) {
                throw std::runtime_error("sampling error");
            }
            t = rec[0].time;
            v = *ptr;
        };

        sim.add_sampler(all_probes, explicit_schedule({(float)tau}), sample_once);


        auto dump = [](cell_member_type, probe_tag, std::size_t n, const sample_record* rec) {
            for (unsigned i = 0; i<n; ++i) {
                auto ptr = util::any_cast<const double*>(rec[i].data);
                if (!ptr) throw std::runtime_error("sampling error");
                std::cout << "debug: t=" << rec[i].time << "; v=" << *ptr << "\n";
            }
        };
        sim.add_sampler(all_probes, regular_schedule(dt), dump);

        sim.run(t_end, dt);

        if (t<0) {
            throw std::runtime_error("sampler never triggered?!");
        }

        double v0 = p.erev;
        double v_analytic = (v0-vinf)*std::exp(-t/tau) + vinf;
        double v_err = std::abs(v_analytic-v);

        results.push_back({p, dt, t, v, v_err});
    }

    return results;
}

int main() {
    pas1comp::param p; // default
    //auto results = run_pas1comp(p, 1, 100);
    auto results = run_pas1comp(p, 100, 100);

    std::cout << "t, v, dt, verr\n";
    for (auto& r: results) {
   //     std::cout << r.t << ", " << r.v << ", " << r.dt << ", " << r.v_err << "\n";
    }
}
