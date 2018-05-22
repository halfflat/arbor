#pragma once

#include <cfloat>
#include <cstddef>
#include <cassert>
#include <memory>
#include <random>
#include <stdexcept>
#include <utility>
#include <vector>

#include "mock_vector.h"

using time_type = float;
using gid_type = unsigned;

using spike_vector = mock_vector<time_type>;
using pse_vector = mock_vector<time_type>;

struct epoch {
    std::ptrdiff_t id = 0;
    time_type t0 = 0;
    time_type t1 = 0;

    epoch() = default;
    epoch(const epoch&) = default;
    epoch(std::ptrdiff_t id, float t0, float t1): id(id), t0(t0), t1(t1) {}

    epoch& operator=(const epoch&) = default;

    epoch advance(float t_step) const {
        assert(t_step>0);
        return epoch(id+1, t1, t1+t_step);
    }
};

template <typename X>
struct ptr_range {
    using iterator = X*;
    using const_iterator = const X*;
    using difference_type = std::ptrdiff_t;
    using size_type = std::size_t;
    using value_type = X;
    using reference = X&;
    using const_reference = const X&;

    X* begin_ = nullptr;
    X* end_ = nullptr;

    ptr_range() = default;
    ptr_range(X* b, X* e): begin_(b), end_(e) {}
    ptr_range(const ptr_range& pr): begin_(pr.begin_), end_(pr.end_) {}

    ptr_range& operator=(const ptr_range&) = default;

    X* begin() const { return begin_; }
    X* end() const { return end_; }
    X* cbegin() const { return begin_; }
    X* cend() const { return end_; }
    X* data() const { return begin_; }

    void swap(ptr_range& other) {
        std::swap(begin_, other.begin_);
        std::swap(end_, other.end_);
    }

    X& operator[](ptrdiff_t i) const { return begin_[i]; }
    X& at(ptrdiff_t i) const {
        return size_type(i)<size()? (*this)[i]: throw std::out_of_range();
    }

    std::size_t size() { return end()-begin(); }
};

struct mock_parameters {
    int n_rank = 1; // multiplier for spike exchange;
    int fanout = 1000; // events per spike;
    time_type min_delay = 10; // [ms]
    time_type mean_spike_rate = .3; // [kHz]

    // busy wait times are pulled from an exponential
    // distribution (or equivalent) with given mean in ms.
    float busy_wait_advance = 0.; // [ms]
    float busy_wait_exchange = 0.; // [ms]

    int rng_seed = 10000;
};

struct cell_group {
    cell_group(std::pair<gid_type, gid_type> gids, const mock_parameters&);
    spike_vector advance(epoch p, const ptr_range<pse_vector>& events);

    std::pair<gid_type, gid_type> gids_;
    time_type rate_;
    time_type t_;

    std::minstd_rand rng_;
    double wait_per_gid_ = 0;

    std::size_t n_delivered = 0;
    std::size_t n_spike = 0;
};

using cell_group_ptr = std::unique_ptr<cell_group>;

struct cell_group_partition {
    std::vector<gid_type> group_divs = {0};

    explicit cell_group_partition(const std::vector<gid_type>& sizes) {
        group_divs.push_back(0);
        std::partial_sum(sizes.begin(), sizes.end(), std::back_inserter(group_divs));
    }

    std::pair<gid_type, gid_type> operator[](int i) {
        return {group_divs.at(i), group_divs.at(i+1)};
    }

    gid_type n_gid() const { return group_divs.back(); }
    gid_type n_cell_groups() const { return group_divs.size()-1; }
};

struct simulation {
    simulation(const cell_group_partition& p, const mock_parameters&);
    virtual void run(time_type tfinal) = 0;

    virtual std::size_t n_ev_queued() const = 0;
    virtual std::size_t n_recv_spike() const = 0;

    std::size_t n_ev_delivered() const {
        std::size_t n = 0;
        for (auto& g: cell_groups_) n += g->n_delivered;
        return n*n_rank_;
    }

    std::size_t n_spike() const {
        std::size_t n = 0;
        for (auto& g: cell_groups_) n += g->n_spike;
        return n*n_rank_;
    }

    std::pair<time_type, time_type> time_minmax() const {
        std::pair<time_type, time_type> minmax(INFINITY, -INFINITY);
        for (auto& g: cell_groups_) {
            minmax.first = std::min(minmax.first, g->t_);
            minmax.second = std::max(minmax.second, g->t_);
        }
        return minmax;
    }

    virtual ~simulation() {}

    time_type min_delay_;
    cell_group_partition cell_group_gids_;
    std::vector<cell_group_ptr> cell_groups_;

    int n_rank_; // multiplier for spike exchange.
    int fanout_; // one spike -> fanout_ events.
    std::minstd_rand rng_;
    double wait_exchange_ = 0;
};


struct serial_simulation: simulation {
    serial_simulation(const cell_group_partition& p, const mock_parameters&);
    void run(time_type tfinal) override;

    void make_event_queues(const spike_vector&, std::vector<pse_vector>&);
    void setup_events(time_type t0, time_type t1);
    spike_vector spike_exchange(const spike_vector&);

    std::size_t n_ev_queued() const override {
        std::size_t n = 0;
        for (auto& evs: event_queues_) n += evs.size();
        return n;
    }

    std::size_t n_recv_spike() const override { return n_recv_spike_; }

    std::size_t n_recv_spike_ = 0;
    std::vector<pse_vector> pending_events_; // events from spike exchange per gid
    std::vector<pse_vector> event_queues_; // events to deliver per gid
};
