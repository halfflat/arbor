#include <cassert>
#include <cstddef>
#include <memory>
#include <random>

#include "counter.h"
#include "mocksim.h"
#include "mock_vector.h"
#include "busywait.h"

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
    rate_(mp.mean_spike_rate),
    t_(0),
    rng_(mp.rng_seed+gids.first*12345),
    wait_per_gid_(mp.busy_wait_advance)
{}

spike_vector cell_group::advance(epoch p, const ptr_range<pse_vector>& events) {
    assert(p.t0==t_);
    int n_gid = gids_.second-gids_.first;

    for (pse_vector ev_vec: events) {
        if (!ev_vec.empty()) {
            assert(ev_vec.front()>=p.t0);
            n_delivered += ev_vec.take_upto(p.t1).size();
        }
    }

    busy_wait_exp(wait_per_gid_, rng_, n_gid);
    t_ = p.t1;

    auto spikes = mock_poisson<time_type>(p.t0, p.t1, rate_*n_gid, rng_);
    n_spike += spikes.size();
    return spikes;
}

// Mock simulation implementation:

// Base class common code:
simulation::simulation(const cell_group_partition& p, const mock_parameters& mp):
    n_rank_(mp.n_rank),
    fanout_(mp.fanout),
    min_delay_(mp.min_delay),
    cell_group_gids_(p),
    cell_groups_(p.n_cell_groups()),
    rng_(mp.rng_seed+11),
    wait_exchange_(mp.busy_wait_exchange)
{
    for (auto i: count_along(cell_groups_)) {
        cell_groups_[i] = std::make_unique<cell_group>(cell_group_gids_[i], mp);
    }
}

// Serial implementation:

void spike_accumulate(spike_vector& a, const spike_vector& l) {
    a.append(l);
}

spike_vector serial_simulation::spike_exchange(const spike_vector& local) {
    spike_vector global(local);
    global.resize(global.size()*n_rank_);

    busy_wait_exp(wait_exchange_, rng_);
    n_recv_spike_ += global.size();
    return global;
}

serial_simulation::serial_simulation(const cell_group_partition& p, const mock_parameters& mp):
    simulation(p, mp),
    pending_events_(p.n_gid()),
    event_queues_(p.n_gid())
{}

void serial_simulation::run(time_type t_final) {
    spike_vector local_spikes, global_spikes;
    time_type t_step = min_delay_/2;
    time_type t_until = std::min(t_final, t_step);

    epoch ep(0, 0.f, t_until);

    while (ep.t0<t_final) {
        setup_events(ep.t0, ep.t1);

        for (auto i: count_along(cell_groups_)) {
            auto gid_ival = cell_group_gids_[i];
            auto spikes = cell_groups_[i]->advance(ep,
                    {&event_queues_[gid_ival.first], &event_queues_[gid_ival.second]});
            spike_accumulate(local_spikes, std::move(spikes));
        }

        global_spikes = spike_exchange(local_spikes);
        local_spikes.clear();
        make_event_queues(global_spikes, pending_events_);

        ep = ep.advance(t_step, t_final);
    }

    setup_events(ep.t0, ep.t1);
}

void serial_simulation::make_event_queues(const spike_vector& global, std::vector<pse_vector>& queues) {
    // Each spike -> n_fanout events, randomly picked across gids. Only 1/n_rank_ of these gids
    // will be local; presume exactly even distribution (we sneakily know that global.size() has been
    // premultiplied by n_rank_).
 
    std::exponential_distribution<time_type> E(min_delay_);
    std::uniform_real_distribution<time_type> U(global.front(), global.back());
    std::uniform_int_distribution<gid_type> G(0, cell_group_gids_.n_gid()-1);

    for (auto i: count(global.size()/n_rank_)) {
        for (auto j: count(fanout_)) {
            time_type delay = min_delay_ + E(rng_);

            queues.at(G(rng_)).push_back(U(rng_)+delay);
        }
    }
}

void merge_events_serial(time_type t0, time_type t1, pse_vector& queue, pse_vector& pending) {
    queue.take_upto(t0);
    queue.append(pending);
    pending.clear();
}

void serial_simulation::setup_events(time_type t0, time_type t1) {
    for (auto i: count_along(event_queues_)) {
        merge_events_serial(t0, t1, event_queues_[i], pending_events_[i]);
    }
}
