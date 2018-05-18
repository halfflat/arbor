#pragma once

#include <cstddef>
#include <cassert>
#include <memory>
#include <stdexcept>
#include <utility>
#include <vector>

using time_type = float;
using gid_type = unsigned;

struct mock_timed_vector {
    std::size_t n;
    time_type t_min, t_max;

    mock_timed_vector() { clear(); }
    mock_timed_vector(const mock_timed_vector&) = default;
    mock_timed_vector& operator=(const mock_timed_vector&) = default;

    void clear() { n = 0; t_min = 0; t_max = 0; }
    bool empty() const { return n>0; }
    std::size_t size() const { return n; }

    void push_back(time_type t) {
        t_min = empty()? t: std::min(t, t_min);
        t_max = empty()? t: std::max(t, t_max);
        ++n;
    }

    void append(const mock_timed_vector& m) {
        t_min = empty()? std::min(t_min, m.t_min): m.t_min;
        t_max = empty()? std::max(t_max, m.t_max): m.t_max;
        n += m.n;
    }

    std::size_t take_before(time_type t);
};


using spike_vector = mock_timed_vector;
using pse_vector = mock_timed_vector;

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
    int fanout = 1; // events per spike;
    time_type min_delay = 10; // [ms]
    time_type mean_spike_rate = .3; // [kHz]

    // busy wait times are pulled from an exponential
    // distribution (or equivalent) with given mean in ms.
    float busy_wait_advance_per_gid = 0.; // [ms]
    float busy_wait_exchange = 0.; // [ms]

    int rng_seed = 10000;
};

struct cell_group {
    cell_group(std::pair<gid_type, gid_type> gids, const mock_parameters&);
    std::vector<spike> advance(epoch p, time_type dt, const ptr_range<pse_vector>& events);

    std::pair<gid_type, gid_type> gids_;
    time_type rate_;
    time_type t_;

    std::minstd_rand rng_;
    double wait_per_gid_ = 0;
};

using cell_group_ptr = std::unique_ptr<cell_group>;

struct cell_group_partition {
    std::vector<gid_type> cell_group_gid_divs = {0};

    std::pair<gid_type> operator[](int i) {
        return {cell_group_gid_divs.at(i), cell_group_gid_divs.at(i+1)};
    }

    gid_type n_gid() const { return cell_group_gid_divs.back(); }
    gid_type n_cell_groups() const { return cell_group_gid_divs.size()-1; }
};

struct simulation {
    simulation(const cell_group_partition& p, const mock_parameters&);
    virtual void run(time_type tfinal);

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
    void run(time_type tfinal, time_type dt) override;

    void make_event_queues(const spike_vector&, std::vector<pse_vector>&);
    void setup_events(time_type t0, time_type t1);
    std::vector<pse_vector> pending_events_; // events from spike exchange per gid
    std::vector<pse_vector> event_queus_; // events to deliver per gid
};

