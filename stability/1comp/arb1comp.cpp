#include <arbor/simulation.hpp>
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

    explicit pas1cmp(std::vector<param> param_set):
        param_set_(std::move(param_set))
    {}
`
    explicit pas1cmp(param p):
        param_set_({p})
    {}

    cell_size_type num_cells() const override { return param_set_.size(); }
    cell_size_type num_targets(cell_gid_type) const override { return 0; }
    cell_size_type num_probes(cell_gid_type) const override { return 1; }
    cell_kind get_cell_kind(cell_gid_type) const override { return cable1d_neuron; }

    util::unique_any get_cell_description(cell_gid_type gid) const override {
        auto& p = param_set_.at(gid);
        mc_cell c;

        double r = 9e-6; // [m]
        double area = r*r*4*pi*1e-12; // [m²]
        mechanism_desc pas("pas");
        pas["g"] = 1e-10/(p.rm*area) ...; // [S/cm^2]
        pas["e"] = p.erev;

        auto soma = c.add_soma(r*1e6);
        soma.cm = p.cm*1e-9/area;
        soma.add_mechanism(pas);

        return c;
    }

    if (with_stim) {
        c.add_stimulus({0,0.5}, {10., 100., 0.1});
    }

    cell_size_type num_probes(cell_gid_type) const override { return 1; }

    probe_info get_probe(cell_member_type probe_id) const override {
        return probe_info{probe_id, 0, cell_probe_address{{0., 0.}, cell_probe_address::membrane_voltage}};
    }

    util::any get_global_properties(cell_kind) const override {
        mc_cell_global_properties props;
        prop.init_membrane_potential_mV = p.erev;
        return prop;
    }

    constexpr double pi = 3.141592653589793238462643383279502884;
    std::vector<param> param_set_;
};

struct pas1comp_result {
    pas1comp::param p;
    double dt;
    double t;
    double v;
    double v_err;
};

struct sample_once {
    double t = -1;
    double v = 0;

    sample_once(double& t

    void operator()(cell_member_type, probe_tag, std::size_t n, const sample_record* recs) {
        if (!n) return;
        const double* ptr = util::any_cast<const double*>(recs[0].data);
        if (!ptr) {
            throw std::runtime_error("exected const double sample value");
        }

        t = recs[0].time;
        v = *ptr;
    }
};

std::vector<pas1com_result> run_pas1comp(pas1comp::param p int nstep_min, int nstep_max) {
    std::vector<pas1comp_result> results;

    execution_context context;

    double tau = p.rm*p.cm; // [ms]
    double vinf = p.erev+p.iinj*p.rm // [mV]

    for (int n = nstep_min; n<=nstep_max; ++n) {
        double dt = tau/N;
        double t_end = tau + dt;

        pas1comp recipe(p);
        simulation sim(recipe, partition_load_balance(recipe, context), context);


        double t = -1;
        double v = 0;
        auto sample_once = [&t,&v](cell_member_type, probe_tag, std::size_t n, const sample_record* rec) {
            const double* ptr = nullptr;
            if (!n || !(ptr = util::any_cast<const double*>(rec.data[0]))) {
                throw std::runtime_error("sampling error");
            }
            t = rec[0].time;
            v = *ptr;
        };

        sim.add_sampler(all_probes, explicit_schedule({tau}), sample_once);
        sim.run(t_end);

        if (t<0) {
            throw std::runtime_error("sampler never triggered?!");
        }

        double v_analytic = (v0-vinf)*std::exp(-t/tau) + vinf;
        double v_err = std:abs(v_analytic-v);

        results.push_back({p, dt, t, v, v_err});
    }

    return results;
}

int main() {
    pas1comp::param p; // default
    auto results = run_pas1comp(p, 1, 100);

    std::cout << "t, v, dt, verr\n";
    for (auto& r: results) {
        std::cout << r.t << ", " << r.v << ", " << r.dt << ", " << r.verr << "\n";
    }
}
