#ifndef FACTORY_TRAJECTORY_SEARCH_HPP
#define FACTORY_TRAJECTORY_SEARCH_HPP

// Search around an executable incumbent rather than rebuilding every plan from
// zero. All checkpoints are actual inventories, including purchased seeds.
namespace factory_trajectory {

struct Event {
    int position, molecule;
    bool purchase;
    double weight;
};

struct Trace {
    vector<vector<ll>> stock, required;
    vector<ll> cost;
    vector<Event> events;
    vector<int> policy;
};

static bool makeTrace(const Plan& plan, int target, ll quantity,
                      Trace& trace, double deadline) {
    const int n = (int)plan.seq.size();
    trace.stock.assign(n + 1, vector<ll>(M));
    trace.cost.assign(n + 1, 0);
    trace.required.assign(n + 1, vector<ll>(M));
    trace.events.clear();
    vector<unordered_map<int, ll>> contribution(M);
    vector<ll> bought(M);
    for (int i = 0; i < n; ++i) {
        if ((i & 31) == 0 && nowSec() >= deadline) return false;
        const auto& reaction = R[plan.seq[i]];
        trace.stock[i + 1] = trace.stock[i];
        trace.cost[i + 1] = trace.cost[i];
        auto& stock = trace.stock[i + 1];
        for (auto [m, q] : reaction.in) {
            ll missing = max(0LL, q - stock[m]);
            if (missing && price[m] < 0) return false;
            if (missing) {
                __int128 paid = (__int128)missing * price[m];
                __int128 next = (__int128)trace.cost[i + 1] + paid;
                if (next >= LINF) return false;
                trace.cost[i + 1] = (ll)next;
                bought[m] += missing;
                trace.events.push_back({i, m, true, (double)paid});
            }
            stock[m] = max(0LL, stock[m] - q);
        }
        for (auto [m, q] : reaction.out) {
            if ((__int128)stock[m] + q >= LINF) return false;
            stock[m] += q;
        }
        for (auto [m, q] : reaction.prod) {
            contribution[m][plan.seq[i]] += q;
            double worth = v0[m] < DINF / 2 ? v0[m] : 0;
            trace.events.push_back({i, m, false, max(0.001, worth * q)});
        }
    }
    if (trace.cost.back() != plan.cost || trace.stock.back()[target] < quantity)
        return false;
    trace.required[n][target] = quantity;
    for (int i = n - 1; i >= 0; --i) {
        trace.required[i] = trace.required[i + 1];
        for (auto [m, q] : R[plan.seq[i]].out)
            trace.required[i][m] = max(0LL, trace.required[i][m] - q);
        for (auto [m, q] : R[plan.seq[i]].in) {
            if ((__int128)trace.required[i][m] + q >= LINF) return false;
            trace.required[i][m] += q;
        }
    }
    trace.policy = plan.choice.size() == (size_t)M ? plan.choice : baseChoice;
    for (int m = 0; m < M; ++m) {
        ll most = bought[m];
        int source = price[m] >= 0 ? CH_BUY : trace.policy[m];
        for (auto [r, q] : contribution[m]) {
            if (q > most || (q == most && r < source)) { most = q; source = r; }
        }
        if (most) trace.policy[m] = source;
    }
    stable_sort(trace.events.begin(), trace.events.end(), [](const Event& a, const Event& b) {
        if (a.purchase != b.purchase) return a.purchase > b.purchase;
        if (a.weight != b.weight) return a.weight > b.weight;
        return a.position < b.position;
    });
    return true;
}

static bool accept(Plan& best, vector<int> candidate, int target, ll quantity) {
    ll cost = 0;
    if (!verifySeq(candidate, target, quantity, cost) || cost >= best.cost) return false;
    // A changed downstream recipe can make an entire old precursor branch dead.
    auto live = factory_sequence_polish::liveSequence(candidate, target, quantity);
    ll liveCost = 0;
    if (live.size() < candidate.size() && verifySeq(live, target, quantity, liveCost)
            && liveCost < cost) {
        candidate = move(live);
        cost = liveCost;
    }
    best.cost = cost;
    best.seq = move(candidate);
    best.proven = false;
    best.dem.clear();
    vector<unsigned char> used(M);
    used[target] = 1;
    for (int r : best.seq) for (auto [m, q] : R[r].in) used[m] = 1;
    for (int m = 0; m < M; ++m) if (used[m]) best.dem.push_back(m);
    return true;
}

static bool replaySuffix(adaptive_detail::Planner& planner, const Plan& plan,
                         const Trace& trace, int begin, int target, ll quantity) {
    if (planner.stock[target] >= quantity) return true;
    for (int i = begin; i < (int)plan.seq.size(); ++i) {
        if (planner.expired()) return false;
        int r = plan.seq[i];
        // Exact sufficient test for redundancy in the unchanged suffix. The
        // reverse requirements include gross catalyst requirements as well.
        bool useful = false;
        for (auto [m, q] : R[r].prod)
            if (planner.stock[m] < trace.required[i + 1][m]) useful = true;
        if (!useful) continue;
        bool ready = false;
        for (int guard = 0; guard < 16; ++guard) {
            ready = true;
            for (auto [m, q] : R[r].in) {
                if (price[m] >= 0 || planner.stock[m] >= q) continue;
                ready = false;
                if (!planner.ensure(m, q, q, false)) return false;
            }
            if (ready) break;
        }
        if (!ready || !planner.fire(r)) return false;
        if (planner.stock[target] >= quantity) return true;
    }
    return planner.stock[target] >= quantity;
}

}  // namespace factory_trajectory

static void trajectoryImprove(Plan& best, int target, ll quantity, double deadline) {
    if (best.cost <= 0 || best.cost >= LINF || best.seq.empty() || nowSec() >= deadline)
        return;
    using namespace factory_trajectory;
    uint64_t seed = 0x6a09e667f3bcc909ULL ^ ((uint64_t)searchSeed << 23)
                  ^ (uint64_t)target * 998244353ULL ^ (uint64_t)quantity;
    adaptive_detail::Planner prepared(best, deadline, seed);
    mt19937 rng((unsigned)(seed ^ (seed >> 32)));
    Trace trace;
    bool dirty = true;
    size_t trial = 0, eventCursor = 0;
    while (nowSec() < deadline) {
        if (dirty) {
            if (!makeTrace(best, target, quantity, trace, deadline)) return;
            trace.events.erase(remove_if(trace.events.begin(), trace.events.end(), [&](const Event& event) {
                if (!event.purchase && price[event.molecule] >= 0) return false;
                for (int r : prepared.source[event.molecule])
                    if (event.purchase || r != best.seq[event.position]) return false;
                return true;
            }), trace.events.end());
            eventCursor = 0;
            dirty = false;
        }
        ++trial;
        const int n = (int)best.seq.size();
        int checkpoint = 0;
        Event event{0, target, false, 0};
        bool targeted = !trace.events.empty() && trial % 3 != 0;
        if (targeted) {
            // Start systematically with the most expensive purchases. A later
            // randomized pass also explores production choices in the tail.
            size_t index = eventCursor++;
            if (index >= trace.events.size()) index = rng() % trace.events.size();
            event = trace.events[index];
            checkpoint = event.position;
        } else if (trial <= 12) {
            checkpoint = (int)((trial % 5) * (size_t)n / 5);
        } else {
            checkpoint = rng() % n;
        }
        if (trace.cost[checkpoint] >= best.cost) continue;
        adaptive_detail::Planner planner = prepared;
        planner.bound = best.cost;
        planner.cost = trace.cost[checkpoint];
        planner.stock = trace.stock[checkpoint];
        planner.sequence.assign(best.seq.begin(), best.seq.begin() + checkpoint);
        planner.cap = min<ll>(500000, max<ll>(10000, best.seq.size() * 4 + 500));
        planner.deadline = min(deadline, nowSec() + min(0.12, 0.012 + 0.000015 * (n - checkpoint)));
        planner.randomState = seed + trial * 11995408973635179863ULL;
        planner.preferred = trial % 4 == 0 && best.choice.size() == (size_t)M
                          ? best.choice : trace.policy;
        planner.preference = trial % 7 == 0 ? 0.6 : 0.005 + 0.04 * planner.uniform();
        planner.noise = trial < 8 ? 0 : 0.05 + 0.25 * planner.uniform();
        planner.credit = (trial % 5) * 0.25;
        planner.ordering = trial < 8 ? 3 : planner.random32() % 4;
        planner.protectAncestors = trial % 5 != 0;
        bool solved = false;
        if (targeted) {
            int m = event.molecule;
            vector<int> sources = planner.source[m];
            int old = best.seq[checkpoint];
            if (!event.purchase)
                sources.erase(remove(sources.begin(), sources.end(), old), sources.end());
            if (!event.purchase && price[m] >= 0) sources.push_back(CH_BUY);
            if (!sources.empty()) {
                int alternative = sources[(trial / 2) % sources.size()];
                if (alternative == CH_BUY) {
                    // Omit this production batch and buy only the eventual
                    // remainder at its consuming reactions.
                    solved = replaySuffix(planner, best, trace, checkpoint + 1, target, quantity);
                } else {
                planner.source[m] = {alternative};
                planner.preferred[m] = alternative;
                ll increase = event.purchase
                    ? max(1LL, qtyOf(R[old].in, m) - planner.stock[m])
                    : max(1LL, min(qtyOf(R[old].prod, m),
                                   max(0LL, trace.required[checkpoint + 1][m] - planner.stock[m])));
                if ((__int128)planner.stock[m] + increase < LINF) {
                    ll need = planner.stock[m] + increase;
                    if (planner.ensure(m, need, need, true) == 1) {
                        int begin = checkpoint + (event.purchase ? 0 : 1);
                        solved = replaySuffix(planner, best, trace, begin, target, quantity);
                    }
                }
                }
            }
        } else {
            if (trial > 8 && planner.uniform() < 0.7) {
                int changes = 1 + rng() % 3;
                for (int j = 0; j < changes && !best.dem.empty(); ++j) {
                    int m = best.dem[rng() % best.dem.size()];
                    if (!planner.source[m].empty())
                        planner.preferred[m] = planner.source[m][rng() % planner.source[m].size()];
                }
            }
            solved = planner.ensure(target, quantity, quantity, true) == 1;
        }
        if (solved && !planner.stopped && planner.cost < best.cost &&
                accept(best, move(planner.sequence), target, quantity)) dirty = true;
        if (best.cost == 0) return;
    }
}

#endif  // FACTORY_TRAJECTORY_SEARCH_HPP
