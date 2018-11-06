/*
 * For a passive one-compartment cell with current injection,
 * arbor will be simulating the behaviour of a simple RC
 * circuit.
 *
 * For fixed RC (and reversal potential) parameters, perform
 * a parameter sweep over injected current and integration
 * time step (as a fraction of the RC time constant).
 */

#include <iostream>
#include <utility>

#include <sup/tinyopt.hpp>

const char* usage_long =
    "[OPTION]...\n"
    "Arbor integration error in a simple RC model, parameterized\n"
    "over injection current and integration dt.\n"
    "\n"
    "  -R, --resistance=R    total membrane resistance R [MΩ]\n"
    "  -C, --capacitance=C   total membrane capacitance C [nF]\n"
    "  -E, --reversal=EREV   reversal potential EREV [mV]\n"
    "  -I, --injection=IINJ  max injection current IINJ [nA]\n"
    "  -d, --dt=DTMIN        minimum time step as proportion of\n"
    "                           time constant τ=RC.\n"
    "  -n, --steps=N         number of currents and time steps to test\n"
    "  -s, --show            show parameter values and exit\n"
    "  -h, --help            display usage information and exit\n"
    "\n"
    "When run with N steps, the test is performed with N currents,\n"
    "selected linearly from I/N to I, and N time steps, selected\n"
    "selected geometrically from DTMIN·τ to τ, the RC time constant.";

const char* usage_short =
    "[OPTION]...\n"
    "Use the --help option for detailed usage information.";

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
    struct param_set {
        double rm = 100;       // [MΩ]
        double cm = 0.01;      // [nF]
        double erev = -65;     // [mV]
        double iinj_max = 0.1; // [nA]
        double dt_min = 0.001; // as proportion of rm*cm
        unsigned n = 5;
        bool help = false;
        bool show = false;
    } params;

    try {
        for (char** arg = argv+1; *arg; ) {
            if (set_from_opt(arg, params.rm, 'R', "resistance") ||
                set_from_opt(arg, params.cm, 'C', "capacitance") ||
                set_from_opt(arg, params.erev, 'E', "reversal") ||
                set_from_opt(arg, params.iinj_max, 'I', "injection") ||
                set_from_opt(arg, params.dt_min, 'd', "dt") ||
                set_from_opt(arg, params.n, 'n', "steps") ||
                flag_from_opt(arg, params.help, 'h', "help") ||
                flag_from_opt(arg, params.show, 's', "show")) continue;

            throw to::parse_opt_error(*arg, "unrecognized argument");
        }

        if (params.help) {
            to::usage(argv[0], usage_long);
            return 0;
        }

        if (params.show) {
            std::cout <<
                "membrane resistance " << params.rm << " MΩ\n" <<
                "membrane capacitance " << params.cm << " nF\n" <<
                "time constant (τ) "  << params.rm*params.cm << " ms\n" <<
                "reversal potential " << params.erev << " mV\n" <<
                "max injected current " << params.iinj_max << " nA\n" <<
                "min integration time step " << params.dt_min << " τ\n" <<
                "number of parameter steps " << params.n << "\n";
            return 0;
        }
    }
    catch (to::parse_opt_error& e) {
        to::usage(argv[0], usage_short, e.what());
        return 1;
    }
}

