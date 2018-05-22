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
};

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

        auto help = [&opt]() { opt.help = true; };

        for (auto arg = argv+1; *arg; ) {
            help              << parse_opt(arg, 'h', "help") ||
            opt.n_cell        << parse_opt<int>(arg, 'n', "cells") ||
            opt.sim_time      << parse_opt<double>(arg, 't', "time") ||
            opt.group_sizes   << parse_opt<std::vector<int>>(arg, 'g', "group-size", csv) ||
            mparam.n_rank     << parse_opt<int>(arg, 'N', "ranks") ||
            mparam.fanout     << parse_opt<int>(arg, 'F', "fanout") ||
            mparam.min_delay  << parse_opt<double>(arg, 'M', "min-delay") ||
            mparam.mean_spike_rate    << parse_opt<double>(arg, 'r', "spike-rate") ||
            mparam.busy_wait_advance  << parse_opt<double>(arg, 0, "advance-time") ||
            mparam.busy_wait_exchange << parse_opt<double>(arg, 0, "exchange-time") ||
            (throw to::parse_opt_error(*arg, "unrecognized option"), false);
        }

        if (opt.help) {
            to::usage(argv[0], usage_str);
            return 0;
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

        sputs() << "events delivered: " << sim->n_ev_delivered();
        sputs() << "events queued:    " << sim->n_ev_queued();
        sputs() << "spikes generated: " << sim->n_spike();
        sputs() << "spikes received: " << sim->n_recv_spike();
        auto mm = sim->time_minmax();
        sputs() << "cell group times: " << mm.first << ',' << mm.second;
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
