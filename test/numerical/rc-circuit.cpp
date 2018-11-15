/*
 * For a passive one-compartment cell with current injection,
 * arbor will be simulating the behaviour of a simple RC
 * circuit.
 *
 * For fixed RC (and reversal potential) parameters, perform
 * a parameter sweep over injected current and integration
 * time step (as a fraction of the RC time constant).
 */

#include <cfloat>
#include <iostream>
#include <utility>
#include <vector>

#include <arbor/simulation.hpp>
#include <arbor/load_balance.hpp>
#include <arbor/mc_cell.hpp>
#include <arbor/recipe.hpp>
#include <sup/tinyopt.hpp>

#include "backends/multicore/fvm.hpp"
#include "mechanisms/decay_inj.hpp"

const char* usage_long =
    "[OPTION]...\n"
    "Arbor integration error in a simple RC model, parameterized\n"
    "over injection current and integration dt.\n"
    "\n"
    "  -R, --resistance=R    total membrane resistance R [MΩ]\n"
    "  -C, --capacitance=C   total membrane capacitance C [nF]\n"
    "  -E, --reversal=EREV   reversal potential EREV [mV]\n"
    "  -I, --injection=IINJ  peak injection current IINJ [nA]\n"
    "  -c, --decay=T         injection current decay time constant [ms]\n"
    "  -d, --dt=MIN[,MAX]    integration time step [ms]\n"
    "  -T, --time=TEND       integration stop time [ms]\n"
    "  -n, --nsteps=N        number of time steps to test if range is given\n"
//    "  -f, --format=FMT      output format (csv or json)\n"
    "  -s, --show            show parameter values and exit\n"
    "  -h, --help            display usage information and exit\n"
    "\n"
    "Initial voltage is set the the reversal potential.\n"
    "Injection current is either constant (if decay time T is zero,\n"
    "the default), or decays exponentially with time constant T.\n"
    "Integration time step (dt) values are interpolated geometrically\n"
    "between dt MIN and MAX.\n"

const char* usage_short =
    "[OPTION]...\n"
    "Use the --help option for detailed usage information.";

struct rc_param_set {
    double iinj0 = 0.1;    // [nA]
    double injtau = 0;     // [ms]
    double rm = 100;       // [MΩ]
    double cm = 0.01;      // [nF]
    double erev = -65;     // [mV]
};

struct range_spec {
    double min, max;
    friend std::istream& operator>>(std::istream& s, range_spec& r) {
        std::istream::sentry sentry(s);

        s >> r.min;
        r.max = r.min;
        if (s && s.peek()==',') {
            s >> r.max;
        }
        return s;
    }
};

struct run_param_set {
    range_spec dt{0.025, 0.025}; // [ms]
    double t_end = 10;     // [ms]
    unsigned nsteps = 10;
};

struct run_rc_result {
    double dt;      // [ms]
    double t_end;   // [ms]
    double iinj;    // [nA]
    double v;       // [mV]
    double v_exact; // [mV]
};

std::vector<run_rc_result> run_rc(run_param_set, rc_param_set);

template <typename T, typename... Args>
bool set_from_opt(char**& argp,T& value, Args&&... args) {
   if (auto optv = to::parse_opt<T>(argp, std::forward<Args>(args)...)) {
       value = optv.value();
       return true;
   }
   return false;
}

template <typename... Args>
bool flag_from_opt(char**& argp, bool& value, Args&&... args) {
   if (auto optv = to::parse_opt(argp, std::forward<Args>(args)...)) {
       value = true;
       return true;
   }
   value = false;
   return false;
}

int main(int argc, char** argv) {
    rc_param_set rc_params;
    run_param_set run_params;

    try {
        bool opt_help = false;
        bool opt_show = false;

        for (char** arg = argv+1; *arg; ) {
            if (set_from_opt(arg, rc_params.rm, 'R', "resistance") ||
                set_from_opt(arg, rc_params.cm, 'C', "capacitance") ||
                set_from_opt(arg, rc_params.erev, 'E', "reversal") ||
                set_from_opt(arg, rc_params.iinj0, 'I', "injection") ||
                set_from_opt(arg, rc_params.injtau, 'c', "decay") ||
                set_from_opt(arg, run_params.dt, 'd', "dt") ||
                set_from_opt(arg, run_params.t_end, 'T', "time") ||
                set_from_opt(arg, run_params.nsteps, 'n', "nsteps") ||
                flag_from_opt(arg, opt_help, 'h', "help") ||
                flag_from_opt(arg, opt_show, 's', "show")) continue;

            throw to::parse_opt_error(*arg, "unrecognized argument");
        }

        if (opt_help) {
            to::usage(argv[0], usage_long);
            return 0;
        }

        if (rc_params.dt.min>=rc_params.dt.max) {
            rc_params.dt.max = rc_params.dt.min;
            rc_params.nsteps = 1;
        }

        if (opt_show) {
            std::cout <<
                "membrane resistance " << rc_params.rm << " MΩ\n" <<
                "membrane capacitance " << rc_params.cm << " nF\n" <<
                "time constant (τ) "  << rc_params.rm*rc_params.cm << " ms\n" <<
                "reversal potential " << rc_params.erev << " mV\n" <<
                "max injected current " << rc_params.iinj_max << " nA\n" <<
                "min integration time step " << rc_params.dt.min << " ms\n" <<
                "max integration time step " << rc_params.dt.max << " ms\n" <<
                "integration end time " << run_params.t_end_over_tau << " τ\n" <<
                "number of currents " << run_params.ncurrents << "\n" <<
                "number of dts " << run_params.nsteps << "\n";
            return 0;
        }
    }
    catch (to::parse_opt_error& e) {
        to::usage(argv[0], usage_short, e.what());
        return 1;
    }

    auto results = run_rc(run_params, rc_params);
    std::cout << "dt, Iinj, t_end, R, C, Erev, v, v_exact, v_err\n";
    for (auto& r: results) {
        double v_err = std::abs(r.v-r.v_exact);
        std::cout
            << r.dt << ", " << r.iinj << ", " << r.t_end << ", "
            << rc_params.rm << ", " << rc_params.cm << ", "
            << rc_params.erev << ", "
            << r.v << ", " << r.v_exact << ", " << v_err << "\n";
    }
}

using namespace arb;

struct pas1comp: recipe {
    // One compartment parameter set; circuit behaviour as below.
    // Current source from time zero.
    // Initial condition v = erev.
    //
    //        cm
    // +------||------+
    // |              |
    // +-\/\/\---+|---+
    // |   rm    erev |
    // |              |
    // +-----(\)------+
    //       iinj

    explicit pas1comp(rc_param_set rc, std::vector<double> iinj):
        rc_(std::move(rc)),
        iinj_(std::move(iinj))
    {}

    cell_size_type num_cells() const override { return iinj_.size(); }
    cell_size_type num_targets(cell_gid_type) const override { return 0; }
    cell_size_type num_probes(cell_gid_type) const override { return 1; }
    cell_kind get_cell_kind(cell_gid_type) const override { return cell_kind::cable1d_neuron; }

    util::unique_any get_cell_description(cell_gid_type gid) const override {
        mc_cell c;

        double r = 9e-6; // [m]
        double area = r*r*4*pi; // [m²]
        mechanism_desc pas("pas");
        pas["g"] = 1e-10/(rc_.rm*area); // [S/cm^2]
        pas["e"] = rc_.erev;

        auto soma = c.add_soma(r*1e6);
        soma->cm = rc_.cm*1e-9/area;
        soma->add_mechanism(pas);

        if (iinjtau>0) {
            mechanism_desc decay_inj("decay_inj");
            decay_inj["iinj0"] = rc_.iinj0;
            decay_inj["tau"] = rc_.iinjtau;
            soma->add_mechanism(decay_inj);
        }
        else {
            c.add_stimulus({0,0.5}, {0, FLT_MAX, (float)iinj_.at(gid)});
        }
        return c;
    }

    probe_info get_probe(cell_member_type probe_id) const override {
        return probe_info{probe_id, 0, cell_probe_address{{0, 0.}, cell_probe_address::membrane_voltage}};
    }

    util::any get_global_properties(cell_kind) const override {
        mc_cell_global_properties props;
        props.init_membrane_potential_mV = rc_.erev;
        props.catalogue = custom_catalogue();

        return props;
    }

private:
    static mechanism_catalogue* custom_catalogue() {
        static mechanism_catalogue cat = make_custom_catalogue();
        return &cat;
    }

    static mechanism_catalogue make_custom_catalogue() {
        mechanism_catalogue cat;
        cat.add("decay_inj", mechanism_decay_inj_info());
        cat.register_implementation("decay_inj", make_mechanism_decay_ink<multicore::backend>());
        return cat;
    }

    static constexpr double pi = 3.141592653589793238462643383279502884;
    rc_param_set rc_;
    std::vector<double> iinj_;
};

std::vector<run_rc_result> run_rc(run_param_set p, rc_param_set rc) {
    std::vector<run_rc_result> results;
    unsigned n_dt = p.nsteps;
    unsigned n_inj = p.ncurrents;

    std::vector<double> iinj(n_inj);
    std::vector<double> vinf(n_inj);

    std::vector<std::vector<double>> sample_t(n_dt);
    std::vector<std::vector<double>> sample_v(n_dt);

    double tau = rc.rm*rc.cm; // [ms]

    results.reserve(n_dt*n_inj);

    auto sample = [&sample_t, &sample_v](cell_member_type m, probe_tag, std::size_t n, const sample_record* rec) {
        for (std::size_t i = 0; i<n; ++i) {
            auto vptr = util::any_cast<const double*>(rec[i].data);
            if (!vptr) throw std::runtime_error("sampling error");

            sample_t.at(m.gid).push_back(rec[i].time);
            sample_v.at(m.gid).push_back(*ptr);
        }
    };

    auto context = make_context();

    for (unsigned i = 0; i<n_dt; ++i) {
        double dt = p.dt.min;
        if (n_dt>1) dt *= std::pow(p.dt_max/p.dt_min, i/(n_dt-1.));

// work from here...
        pas1comp recipe(rc, iinj);
        simulation sim(recipe, partition_load_balance(recipe, context), context);

        std::fill(sample_t.begin(), sample_t.end(), -1.f);
        std::fill(sample_v.begin(), sample_v.end(), 0.f);

        sim.add_sampler(all_probes, explicit_schedule({(float)tau}), sample_once);
        sim.run(t_end, dt);

        for (unsigned j = 0; j<n_inj; ++j) {
            if (sample_t[j]<0) throw std::runtime_error("sampler never triggered?!");

            double v_analytic = (rc.erev-vinf[j])*std::exp(-sample_t[j]/tau) + vinf[j];
            results.push_back({dt, sample_t[j], iinj[j], sample_v[j], v_analytic});
        }
    }

    return results;
}

