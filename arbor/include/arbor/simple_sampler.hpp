#pragma once

/*
 * Simple(st?) implementation of a recorder of scalar
 * trace data from a cell probe, with some metadata.
 */

#include <stdexcept>
#include <type_traits>
#include <vector>

#include <arbor/common_types.hpp>
#include <arbor/sampling.hpp>
#include <arbor/util/any_ptr.hpp>
#include <arbor/util/optional.hpp>

namespace arb {

template <typename V>
struct trace_entry {
    time_type t;
    V v;
};

// `trace_data` wraps a std::vector of trace_entry with
// a copy of the probe-specific metadata associated with a probe.
//
// If `Meta` is void, ignore any metadta.

template <typename V, typename Meta = void>
struct trace_data: public std::vector<trace_entry<V>> {
    Meta meta;
    explicit operator bool() const { return !this->empty(); }
};

template <typename V>
struct trace_data<V, void>: public std::vector<trace_entry<V>> {
    explicit operator bool() const { return !this->empty(); }
};

// `trace_vector` is simply a vector of `trace_data`.
// When there are multiple probes associated with a probe id,
// the ith element will correspond to the sample data obtained
// from the probe with index i.

template <typename V, typename Meta = void>
using trace_vector = std::vector<trace_data<V, Meta>>;

// Add a bit of smarts for collecting variable-length samples which are
// passed back as a pair of pointers describing a range; these can
// be used to populate a trace of vectors.

template <typename V>
struct trace_push_back {
    template <typename Meta>
    static bool push_back(trace_data<V, Meta>& trace, const sample_record& rec) {
        if (auto p = util::any_cast<const V*>(rec.data)) {
            trace.push_back({rec.time, *p});
            return true;
        }
        return false;
    }
};

template <typename V>
struct trace_push_back<std::vector<V>> {
    template <typename Meta>
    static bool push_back(trace_data<std::vector<V>, Meta>& trace, const sample_record& rec) {
        if (auto p = util::any_cast<const std::vector<V>*>(rec.data)) {
            trace.push_back({rec.time, *p});
            return true;
        }
        else if (auto p = util::any_cast<const std::pair<const V*, const V*>*>(rec.data)) {
            trace.push_back({rec.time, std::vector<V>(p->first, p->second)});
            return true;
        }
        return false;
    }
};

template <typename V, typename Meta = void>
class simple_sampler {
public:
    explicit simple_sampler(trace_vector<V, Meta>& trace): trace_(trace) {}

    // TODO: C++17 use if constexpr to test for Meta = void case.
    void operator()(probe_metadata pm, std::size_t n, const sample_record* recs) {
        const Meta* m = util::any_cast<const Meta*>(pm.meta);
        if (!m) {
            throw std::runtime_error("unexpected metadata type in simple_sampler");
        }

        if (trace_.size()<=pm.index) {
            trace_.resize(pm.index+1);
        }

        if (trace_[pm.index].empty()) {
            trace_[pm.index].meta = *m;
        }

        for (std::size_t i = 0; i<n; ++i) {
            if (!trace_push_back<V>::push_back(trace_[pm.index], recs[i])) {
                throw std::runtime_error("unexpected sample type in simple_sampler");
            }
        }
    }

private:
    trace_vector<V, Meta>& trace_;
};

template <typename V>
class simple_sampler<V, void> {
public:
    explicit simple_sampler(trace_vector<V, void>& trace): trace_(trace) {}

    void operator()(probe_metadata pm, std::size_t n, const sample_record* recs) {
        if (trace_.size()<=pm.index) {
            trace_.resize(pm.index+1);
        }

        for (std::size_t i = 0; i<n; ++i) {
            if (!trace_push_back<V>::push_back(trace_[pm.index], recs[i])) {
                throw std::runtime_error("unexpected sample type in simple_sampler");
            }
        }
    }

private:
    trace_vector<V, void>& trace_;
};

template <typename V, typename Meta>
inline simple_sampler<V, Meta> make_simple_sampler(trace_vector<V, Meta>& trace) {
    return simple_sampler<V, Meta>(trace);
}

} // namespace arb
