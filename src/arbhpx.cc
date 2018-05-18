#include <cassert>
#include <cstddef>
#include <memory>
#include <random>

#include "arbhpx.h"
#include "busywait.h"

mock_timed_vector::take_before(time_type t) {
    if (empty()) return 0;

    std::size_t n_taken = 0;

    if (t>=t_max) {
        n_taken = n;
    }
    else if (t>=t_min) {
        n_taken = std::size_t((n*(t-t_min))/(t_max-t_min));
        if (n_taken>n) n_taken = n;
    }

    t_min = std::min(std::max(t_min, t), t_max);
    n -= n_taken;

    return n_taken;
}

template <typename R>
mock_timed_vector mock_poisson(time_type t0, time_type t1, time_type rate, R& rng) {
    mock_timed_vector v;

    std::exponential_distribution<float> E(lambda);
    std::poisson_distribution<int> P;

    float ev_t0 = t0 + E(rng);
    float ev_t1 = t1 - E(rng);
    int n = 0;

    if (ev_t0<t1) {
        v.push_back(ev_t0);
        if (ev_t1>ev_t0) {
            v.push_back(ev_t1);
            ev.n += P(rng, lambda*(ev_t1-ev_t0));
        }
    }

    return v;
}

template <typename I>
struct count_range {
    struct proxy {
        I counter;
        proxy(I x): counter(x) {}
        const I& operator*() const { return counter; }
        proxy& operator++() { ++counter; return *this; }
        bool operator==(const proxy& p) const { return counter==p.counter; }
    };

    I begin() const { return proxy(i0); }
    I end() const { return proxy(i1); }

    count_range() = default;
    explicit count_range(I n): i0(I{}), i1(n) {}
    explicit count_range(I i0, I i1): i0(i0), i1(i1) {}

    I i0, i1;
};

template <typename C>
auto size(const C& c) { return c.size(); }

template <typename X, std::size_t n>
auto size(const X (&c)[n]) { return n; }

template <typename I>
auto count(I n) { return count_range<I>(n); }

template <typename C>
auto count_along(const C& c) { return count(size(c)); }

template <typename I, typename J>
auto span(I i0, J i1) { return count_range<decltype(true? i0: i1)>(i0, i1); }

// Busy wait

template <typename R>
void busy_wait_exp(double t_mean_ms, R& rng, int k = 1) {
    if (t_mean_ms<=0) return;

    double t_ms =
        k==1?
            std::exponential_distribution<double>(t_mean_ms)(rng):
            std::gamma_distribution<double>(k, t_mean_ms)(rng);

    busy_wait_ns(t_ms*1e6);
}

// Mock cell group implementation:

cell_group::cell_group(std::pair<gid_type, gid_type> gids, const mock_parameters& mp):
    gids_(gids),
    rate_(rate),
    t_(0),
    rng_(mp.rng_seed+gids.first*12345),
    wait_per_gid_(mp.busy_wait_advance_per_gid)
{}

spike_vector cell_group::advance(epoch p, time_type dt, const ptr_range<pse_vector>& events) {
    assert(p.t0==t_);

    for (const pse_vector& ev_vec: events) {
        assert(ev_vec.t_min>=p.t0);
    }

    busy_wait_exp(wait_per_gid_);
    t_ = p.t1;

    float lambda = rate_*(gids_.second-gids_.first);
    return mock_poisson(p.t0, p.t1, lambda, rng_);
}

// Mock simulation implementation:

// Base class common code:
simulation::simulation(cell_group_partition p, const mock_parameters& mp):
    n_rank_(mp.n_rank),
    min_delay_(mp.min_delay),
    cell_group_gids_(p),
    cell_groups_(p.n_cell_groups()),
    rng_(mp.rng_seed+11),
    wait_exchange(mp.busy_wait_exchange)
{
    for (auto i: count_along(cell_groups_)) {
        cell_groups_[i] = std::make_unique(new cell_group(cell_group_gids_[i], mp.mean_spike_rate));
    }
}

// Serial implementation:

void spike_accumulate(spike_vector& a, const spike_vector& l) {
    a.append(l);
}

spike_vector spike_exchange(const spike_vector& local, int n_rank) {
    spike_vector global(local);
    global.n *= n_rank;

    busy_wait_exp(wait_exchange_);
    return global;
}

serial_simulation::serial_simulation(const cell_group_partition& p, const mock_parameters& mp):
    simulation(p, mp),
    pending_events_(p.n_gid()),
    event_queues_(p.n_gid())
{}

void serial_simulation::run(time_type t_final, time_type dt) {
    spike_vector local_spikes, global_spikes;
    time_type t_step = min_delay_/2;
    time_type t_until = std::min(t_final, t_step);

    epoch ep(0, 0.f, t_until);

    while (ep.t0<t_final) {
        epoch next_ep = ep.advance(t_step);
        setup_events(next_ep.t0, next_ep.t1);

        for (auto& i: count_along(cell_groups_)) {
            auto gid_ival = cell_group_gids_[i];
            auto spikes = cell_groups_[i]->run(ep, dt,
                    {&event_queues_[gid_ival.first], &event_queues_[gid_ival.second]});
            spike_accumulate(local_spikes, std::move(spikes));
        }

        global_spikes = spike_exchange(local_spikes, n_rank_);
        local_spikes.clear();
        make_event_queues(global_spikes, pending_events_);

        ep = next_ep;
    }
}

void serial_simulation::make_event_queues(const spike_vector& global, std::vector<pse_vector>& queues) {
    std::exponential_distribution<time_type> E(min_delay_);
    std::uniform_real_distribution<time_type> U(global.t_min, global.t_max);
    std::uniform_int_distribution<gid_type> G(0, cell_group_gids_.n_gid()-1);

    // each spike -> n_fanout events, randomly picked across gids.
    for (auto i: count(global.n)) {
        for (auto j: count(fanout_)) {
            time_type delay = min_delay_ + E(rng_);
            queues.at(G(rng_)).push_back(U(rng_)+delay);
        }
    }
}

void merge_events_serial(time_type t0, time_type t1, pse_vector& queue, pse_vector& pending) {
    queue.take_until(t0);
    queue.append(pending);
}

void serial_simulation::setup_events(time_type t0, time_type t1) {
    for (auto i: count_along(event_queues_)) {
        merge_events_serial(t0, t1, event_queues_[i], pending_events_[i]);
    }
}
