#include <iostream>
#include <type_traits>
#include <vector>

#include "mocksim.h"
#include "tinyopt.h"

const std::string usage_str =
"mocksim [OPTION]\n"
"  -n, --cells=N          simulate N cells per (virtual) rank [default: 1]\n"
"  -t, --time=T           run for T ms of simulated time [default: 1000]\n"
"  -g, --group-size=LIST  take group sizes from values in comma-separated\n"
"                           LIST [default: 1]\n"
"  -N, --ranks=N          simulate N ranks [default: 1]\n"
"  -F, --fanout=N         simulate N events per spike [default: 1000]\n"
"  -M, --min-delay=T      set minimum event delivery delay to T ms [default: 10]\n"
"  -r, --spike-rate=R     set cell spike rate to R kHz [default: 0.3]\n"
"  -s, --seed=K           set RNG seed to integer K [default: 10000]\n"
"  --advance-time=T       busy-wait T ms per cell in advance method [default: 0]\n"
"  --exchange-time=T      busy-wait T ms in exchange task [default: 0]\n"
"\n"
"  -h, --help          display this help and exit\n";

struct global_options {
    int n_cell = 1;
    double sim_time = 1000.;
    std::vector<int> group_sizes = {1};
    bool help = false;
    bool verbose = false;
};

template <typename T>
struct delimited_wrapper {
    const T& items;
    const std::string& delimiter;

    delimited_wrapper(const T& items, const std::string& delimiter):
        items(items), delimiter(delimiter) {}

    friend std::ostream& operator<<(std::ostream& out, const delimited_wrapper& w) {
        bool tail = false;
        for (auto& item: w.items) {
            if (tail) out << w.delimiter;
            out << item;
            tail = true;
        }
        return out;
    }
};

template <typename T>
delimited_wrapper<T> delimited(const T& item, const std::string& delimiter) {
    return delimited_wrapper<T>(item, delimiter);
}

struct sputs: public std::ostream {
    sputs(): sputs(std::cout) {}
    sputs(std::ostream& o): std::ostream(o.rdbuf()) {}
    ~sputs() { *this << '\n'; }
};

int main(int argc, char** argv) {
    try {
        using namespace to;
        mock_parameters mparam;
        global_options opt;
        delimited_parser<int> csv(",");

        for (auto arg = argv+1; *arg; ) {
            [&]() {opt.help = true; }    << parse_opt(arg, 'h', "help") ||
            [&]() {opt.verbose = true; } << parse_opt(arg, 'v', "verbose") ||
            opt.n_cell                   << parse_opt<int>(arg, 'n', "cells") ||
            opt.sim_time                 << parse_opt<double>(arg, 't', "time") ||
            opt.group_sizes              << parse_opt<std::vector<int>>(arg, 'g', "group-size", csv) ||
            mparam.n_rank                << parse_opt<int>(arg, 'N', "ranks") ||
            mparam.fanout                << parse_opt<int>(arg, 'F', "fanout") ||
            mparam.min_delay             << parse_opt<double>(arg, 'M', "min-delay") ||
            mparam.mean_spike_rate       << parse_opt<double>(arg, 'r', "spike-rate") ||
            mparam.rng_seed              << parse_opt<int>(arg, 's', "seed") ||
            mparam.busy_wait_advance     << parse_opt<double>(arg, 0, "advance-time") ||
            mparam.busy_wait_exchange    << parse_opt<double>(arg, 0, "exchange-time") ||
            (throw to::parse_opt_error(*arg, "unrecognized option"), false);
        }

        if (opt.help) {
            to::usage(argv[0], usage_str);
            return 0;
        }

        if (opt.verbose) {
            std::cout
                << "option summary:\n"
                << "number of cells:         " << opt.n_cell << "\n"
                << "simulation time:         " << opt.sim_time << "\n"
                << "cell group sizes:        " << delimited(opt.group_sizes, ", ") << "\n"
                << "virtual ranks:           " << mparam.n_rank << "\n"
                << "spike fanout:            " << mparam.fanout << "\n"
                << "min spike delay:         " << mparam.min_delay << "\n"
                << "mean spike rate:         " << mparam.mean_spike_rate << "\n"
                << "rng seed:                " << mparam.rng_seed << "\n"
                << "advance busy-wait time:  " << mparam.busy_wait_advance << "\n"
                << "exchange busy-wait time: " << mparam.busy_wait_exchange << "\n"
                << "\n";
        }

        std::vector<gid_type> groups;
        if (opt.group_sizes.empty()) {
            groups.assign(opt.n_cell, 1);
        }
        else {
            for (int n = 0, k = 0; n<opt.n_cell; ++k) {
                if (k>=opt.group_sizes.size()) k = 0;

                int s = opt.group_sizes[k];
                if (s<=0) s = 1;
                if (n+s>opt.n_cell) s = opt.n_cell-n;

                groups.push_back(s);
                n += s;
            }
        }

        cell_group_partition partn(groups);
        auto sim = std::unique_ptr<simulation>(new serial_simulation(partn, mparam));
        sim->run(opt.sim_time);

        auto group_time_span = sim->time_minmax();
        if (opt.verbose) {
            std::cout
                << "simulation summary:\n"
                << "events delivered: " << sim->n_ev_delivered() << "\n"
                << "events queued:    " << sim->n_ev_queued() << "\n"
                << "spikes generated: " << sim->n_spike() << "\n"
                << "spikes received:  " << sim->n_recv_spike() << "\n"
                << "cell group times: " << group_time_span.first << "--" << group_time_span.second << "\n";
        }

        // Check event counts:
        auto n_ev_total = sim->n_ev_delivered()+sim->n_ev_queued();
        auto n_spike_out = sim->n_spike();
        auto n_spike_in = sim->n_recv_spike();
        if (n_spike_in!=n_spike_out) {
            std::cout << "error: spike discrepancy (in/out): " << n_spike_in << '/' << n_spike_out << "\n";
        }
        auto n_ev_expected = n_spike_out*mparam.fanout;
        if (n_ev_total!=n_ev_expected) {
            std::cout << "error: spike--event discrepancy (events/spike*fanout): "
                << n_ev_total << '/' << n_ev_expected << "\n";
        }
        if (group_time_span.first!=group_time_span.second) {
            std::cout << "error: cell group time discrepancy (min/max): "
                << group_time_span.first << '/' << group_time_span.second << "\n";
        }
        if (group_time_span.second!=opt.sim_time) {
            std::cout << "error: cell group time--sim time discrepancy (group/sim): "
                << group_time_span.second << '/' << opt.sim_time << "\n";
        }
        return 0;
    }
    catch (to::parse_opt_error& e) {
        to::usage(argv[0], usage_str, e.what());
    }
    catch (std::exception& e) {
        std::cerr << "caught exception: " << e.what() << "\n";
    }
    return 1;
}
