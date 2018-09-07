#include <cmath>
#include <exception>
#include <iostream>
#include <fstream>
#include <memory>
#include <vector>

#include <arbor/context.hpp>
#include <arbor/load_balance.hpp>
#include <arbor/mc_cell.hpp>
#include <arbor/profile/meter_manager.hpp>
#include <arbor/profile/profiler.hpp>
#include <arbor/sampling.hpp>
#include <arbor/schedule.hpp>
#include <arbor/simulation.hpp>
#include <arbor/util/any.hpp>
#include <arbor/version.hpp>


#include <aux/ioutil.hpp>
#include <aux/json_meter.hpp>
#include <aux/path.hpp>
#include <aux/spike_emitter.hpp>
#include <aux/strsub.hpp>
#ifdef ARB_MPI_ENABLED
#include <aux/with_mpi.hpp>
#include <mpi.h>
#endif

#include "io.hpp"
#include "miniapp_recipes.hpp"
#include "trace.hpp"

using namespace arb;

using util::any_cast;

void banner(const context&);
std::unique_ptr<recipe> make_recipe(const io::cl_options&, const probe_distribution&);
sample_trace make_trace(const probe_info& probe);
std::fstream& open_or_throw(std::fstream& file, const aux::path& p, bool exclusive = false);
void report_compartment_stats(const recipe&);

int main(int argc, char** argv) {
    bool root = true;
    int rank = 0;

    try {
#ifdef ARB_MPI_ENABLED
        aux::with_mpi guard(argc, argv, false);
        auto context = arb::make_context(arb::proc_allocation(), MPI_COMM_WORLD);
        {
            MPI_Comm_rank(MPI_COMM_WORLD, &rank);
            root = rank==0;
        }
#else
        auto context = arb::make_context();
#endif
#ifdef ARB_PROFILE_ENABLED
        profile::profiler_initialize(context);
#endif
        profile::meter_manager meters;
        meters.start(context);

        std::cout << aux::mask_stream(root);

        // read parameters
        io::cl_options options = io::read_options(argc, argv, root);

        // TODO: add dry run mode

        // Use a node description that uses the number of threads used by the
        // threading back end, and 1 gpu if available.
        banner(context);

        meters.checkpoint("setup", context);

        // determine what to attach probes to
        probe_distribution pdist;
        pdist.proportion = options.probe_ratio;
        pdist.all_segments = !options.probe_soma_only;

        auto recipe = make_recipe(options, pdist);
        if (options.report_compartments) {
            report_compartment_stats(*recipe);
        }

        auto decomp = partition_load_balance(*recipe, context);
        simulation sim(*recipe, decomp, context);

        // Set up samplers for probes on local cable cells, as requested
        // by command line options.
        std::vector<sample_trace> sample_traces;
        for (const auto& g: decomp.groups) {
            if (g.kind==cell_kind::cable1d_neuron) {
                for (auto gid: g.gids) {
                    if (options.trace_max_gid && gid>*options.trace_max_gid) {
                        continue;
                    }

                    for (cell_lid_type j = 0; j<recipe->num_probes(gid); ++j) {
                        sample_traces.push_back(make_trace(recipe->get_probe({gid, j})));
                    }
                }
            }
        }

        auto ssched = regular_schedule(options.sample_dt);
        for (auto& trace: sample_traces) {
            sim.add_sampler(one_probe(trace.probe_id), ssched, make_simple_sampler(trace.samples));
        }

        // Specify event binning/coalescing.
        auto binning_policy =
            options.bin_dt==0? binning_kind::none:
            options.bin_regular? binning_kind::regular:
            binning_kind::following;

        sim.set_binning_policy(binning_policy, options.bin_dt);

        // Initialize the spike exporting interface
        std::fstream spike_out;
        if (options.spike_file_output) {
            using std::ios_base;

            aux::path p = options.output_path;
            p /= aux::strsub("%_%.%", options.file_name, rank, options.file_extension);

            if (options.single_file_per_rank) {
                spike_out = aux::open_or_throw(p, ios_base::out, !options.over_write);
                sim.set_local_spike_callback(aux::spike_emitter(spike_out));
            }
            else if (rank==0) {
                spike_out = aux::open_or_throw(p, ios_base::out, !options.over_write);
                sim.set_global_spike_callback(aux::spike_emitter(spike_out));
            }
        }

        meters.checkpoint("model-init", context);

        // run model
        sim.run(options.tfinal, options.dt);

        meters.checkpoint("model-simulate", context);

        // output profile and diagnostic feedback
        auto profile = profile::profiler_summary();
        std::cout << profile << "\n";
        std::cout << "\nthere were " << sim.num_spikes() << " spikes\n";

        // save traces
        auto write_trace = options.trace_format=="json"? write_trace_json: write_trace_csv;
        for (const auto& trace: sample_traces) {
            write_trace(trace, options.trace_prefix);
        }

        auto report = profile::make_meter_report(meters, context);
        std::cout << report;
        if (root) {
            std::ofstream fid;
            fid.exceptions(std::ios_base::badbit | std::ios_base::failbit);
            fid.open("meters.json");
            fid << std::setw(1) << aux::to_json(report) << "\n";
        }
    }
    catch (io::usage_error& e) {
        // only print usage/startup errors on master
        std::cerr << aux::mask_stream(root);
        std::cerr << e.what() << "\n";
        return 1;
    }
    catch (std::exception& e) {
        std::cerr << e.what() << "\n";
        return 2;
    }
    return 0;
}

void banner(const context& ctx) {
    std::cout << "==========================================\n";
    std::cout << "  Arbor miniapp\n";
    std::cout << "  - distributed : " << arb::num_ranks(ctx)
              << (arb::has_mpi(ctx)? " (mpi)": " (serial)") << "\n";
    std::cout << "  - threads     : " << arb::num_threads(ctx) << "\n";
    std::cout << "  - gpus        : " << (arb::has_gpu(ctx)? "yes": "no") << "\n";
    std::cout << "==========================================\n";
}

std::unique_ptr<recipe> make_recipe(const io::cl_options& options, const probe_distribution& pdist) {
    basic_recipe_param p;

    if (options.morphologies) {
        std::cout << "loading morphologies...\n";
        p.morphologies.clear();
        load_swc_morphology_glob(p.morphologies, options.morphologies.value());
        std::cout << "loading morphologies: " << p.morphologies.size() << " loaded.\n";
    }
    p.morphology_round_robin = options.morph_rr;

    p.num_compartments = options.compartments_per_segment;

    // TODO: Put all recipe parameters in the recipes file
    p.num_synapses = options.all_to_all? options.cells-1: options.synapses_per_cell;
    p.synapse_type = options.syn_type;

    if (options.all_to_all) {
        return make_basic_kgraph_recipe(options.cells, p, pdist);
    }
    else if (options.ring) {
        return make_basic_ring_recipe(options.cells, p, pdist);
    }
    else {
        return make_basic_rgraph_recipe(options.cells, p, pdist);
    }
}

sample_trace make_trace(const probe_info& probe) {
    std::string name = "";
    std::string units = "";

    auto addr = any_cast<cell_probe_address>(probe.address);
    switch (addr.kind) {
    case mc_cell_probe_kind::voltage:
        name = "v";
        units = "mV";
        break;
    case mc_cell_probe_kind::current_density:
        name = "i";
        units = "mA/cm²";
        break;
    default: ;
    }
    name += addr.location.segment? "dend" : "soma";

    return sample_trace{probe.id, name, units};
}

void report_compartment_stats(const recipe& rec) {
    std::size_t ncell = rec.num_cells();
    std::size_t ncomp_total = 0;
    std::size_t ncomp_min = std::numeric_limits<std::size_t>::max();
    std::size_t ncomp_max = 0;

    for (std::size_t i = 0; i<ncell; ++i) {
        std::size_t ncomp = 0;
        auto c = rec.get_cell_description(i);
        if (auto ptr = any_cast<mc_cell>(&c)) {
            ncomp = ptr->num_compartments();
        }
        ncomp_total += ncomp;
        ncomp_min = std::min(ncomp_min, ncomp);
        ncomp_max = std::max(ncomp_max, ncomp);
    }

    std::cout << "compartments/cell: min=" << ncomp_min <<"; max=" << ncomp_max << "; mean=" << (double)ncomp_total/ncell << "\n";
}

