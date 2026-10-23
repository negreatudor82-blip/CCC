#pragma once

// This is a heuristic incumbent improver. Every accepted result is replayed by
// verifySeq; the estimates below are never used as correctness bounds.
namespace adaptive_detail {

struct Planner {
    vector<ll> stock, request;
    vector<unsigned char> active, touched;
    vector<int> sequence, preferred;
    vector<vector<int>> source;
    vector<double> value;
    ll cost = 0, bound, cap;
    double deadline, credit = 0, noise = 0, preference = 1;
    int ordering = 0, calls = 0, depth = 0;
    uint64_t randomState;
    bool stopped = false, protectAncestors = true;

    uint32_t random32() {
        randomState ^= randomState << 13;
        randomState ^= randomState >> 7;
        randomState ^= randomState << 17;
        return (uint32_t)(randomState >> 16);
    }
    double uniform() { return random32() * (1.0 / 4294967296.0); }

    Planner(const Plan& incumbent, double until, uint64_t seed)
        : stock(M, 0), request(M, 0), active(M, 0), touched(M, 0),
          preferred(incumbent.choice.size() == (size_t)M ? incumbent.choice : baseChoice),
          source(M), value(M), bound(incumbent.cost), deadline(until), randomState(seed) {
        cap = incumbent.cost < LINF
            ? min<ll>(500000, max<ll>(5000, 3 * (ll)incumbent.seq.size() + 500))
            : 200000;
        vector<unsigned char> reachable(M, 0);
        for (int m = 0; m < M; ++m) reachable[m] = price[m] >= 0;
        for (int pass = 0; pass <= M; ++pass) {
            bool changed = false;
            for (const auto& reaction : R) {
                bool feasible = true;
                for (auto [m, q] : reaction.in)
                    if (!reachable[m]) { feasible = false; break; }
                if (!feasible) continue;
                for (auto [m, q] : reaction.out)
                    if (!reachable[m]) { reachable[m] = 1; changed = true; }
            }
            if (!changed) break;
        }
        for (int r = 0; r < (int)R.size(); ++r) {
            bool feasible = true;
            for (auto [m, q] : R[r].in)
                if (!reachable[m]) { feasible = false; break; }
            if (feasible)
                for (auto [m, q] : R[r].prod) source[m].push_back(r);
        }
        for (int m = 0; m < M; ++m) {
            value[m] = v0[m] < DINF / 2 ? max(0.0, v0[m]) : 1e12;
            if (price[m] >= 0) value[m] = min(value[m], (double)price[m]);
        }
    }

    bool expired() {
        if (stopped) return true;
        if (++calls > cap * 24 || (calls % 64 == 0 && nowSec() >= deadline))
            stopped = true;
        return stopped;
    }

    // Seed molecules bought by a self-dependent reaction count as stock for its
    // first run. They are still paid for only by fire(), never speculatively.
    ll runsFor(int r, int m, ll horizon) const {
        ll input = qtyOf(R[r].in, m), output = qtyOf(R[r].out, m);
        ll yield = output - input;
        ll initial = stock[m];
        if (price[m] >= 0) initial = max(initial, input);
        if (yield <= 0) return output >= horizon && initial <= output ? 1 : LINF;
        ll missing = max<ll>(0, horizon - initial);
        return max<ll>(1, missing / yield + (missing % yield != 0));
    }

    double estimate(int r, int m, ll horizon) const {
        if (r == CH_BUY) return (double)max<ll>(0, horizon - stock[m]) * price[m];
        for (auto [a, q] : R[r].cons)
            if (protectAncestors && a != m && active[a] && stock[a] - min(stock[a], q) < request[a]) return DINF;
        ll runs = runsFor(r, m, horizon);
        if (runs >= LINF) return DINF;
        double total = 0;
        for (auto [a, q] : R[r].in) {
            ll consumed = max<ll>(0, q - qtyOf(R[r].out, a));
            long double needed = (long double)q + (long double)(runs - 1) * consumed;
            long double missing = max((long double)0, needed - stock[a]);
            total += (double)missing * value[a];
        }
        double usefulCredit = 0;
        for (auto [a, q] : R[r].prod) {
            if (a == m) continue;
            double made = (double)runs * q;
            double wanted = (double)max<ll>(0, request[a] - stock[a]);
            // Conservative variants only credit live demand. Other variants
            // deliberately explore recipes with useful future coproducts.
            double usable = min(made, wanted) + 0.15 * max(0.0, made - wanted);
            usefulCredit += usable * value[a];
        }
        return max(0.0, total - credit * usefulCredit);
    }

    vector<pair<double, int>> candidates(int m, ll horizon, bool top) {
        vector<pair<double, int>> options;
        if (!top && price[m] >= 0) {
            double score = estimate(CH_BUY, m, horizon);
            if (preferred[m] == CH_BUY) score *= preference;
            if (noise) score = score * (1 + noise * (2 * uniform() - 1)) + 1e-6 * uniform();
            options.push_back({score, CH_BUY});
        }
        for (int r : source[m]) {
            double score = estimate(r, m, horizon);
            if (score >= DINF / 2) continue;
            if (preferred[m] == r) score *= preference;
            if (noise) score = score * (1 + noise * (2 * uniform() - 1)) + 1e-6 * uniform();
            options.push_back({score, r});
        }
        sort(options.begin(), options.end());
        return options;
    }

    bool fire(int r) {
        if (expired() || (ll)sequence.size() >= cap) return false;
        __int128 nextCost = cost;
        for (auto [m, q] : R[r].in) {
            ll missing = max<ll>(0, q - stock[m]);
            if (missing && price[m] < 0) return false;
            if (missing) nextCost += (__int128)missing * price[m];
        }
        if (nextCost >= bound || nextCost >= LINF) return false;
        for (auto [m, q] : R[r].in) stock[m] = max<ll>(0, stock[m] - q);
        for (auto [m, q] : R[r].out) {
            if (stock[m] > LINF - q) return false;
            stock[m] += q;
        }
        cost = (ll)nextCost;
        sequence.push_back(r);
        return true;
    }

    struct RequestScope {
        vector<ll>& all;
        vector<pair<int, ll>> old;
        explicit RequestScope(vector<ll>& r) : all(r) {}
        void raise(int m, ll quantity) {
            for (auto [a, q] : old) if (a == m) { all[m] = max(all[m], quantity); return; }
            old.push_back({m, all[m]});
            all[m] = max(all[m], quantity);
        }
        ~RequestScope() { for (auto [m, q] : old) all[m] = q; }
    };

    // Result: 0 = failed, 1 = physically in stock, 2 = defer a purchase until
    // the consuming reaction. The latter avoids inventing an output operation.
    int ensure(int m, ll needed, ll horizon, bool top = false) {
        if (stock[m] >= needed) return 1;
        if (expired() || depth >= min(192, M + 12)) return 0;
        touched[m] = 1;
        if (active[m]) return !top && price[m] >= 0 ? 2 : 0;
        active[m] = 1;
        ++depth;
        RequestScope scope(request);
        horizon = max(needed, horizon);
        scope.raise(m, horizon);
        int result = 0, iterations = 0;
        while (stock[m] < needed && !expired()) {
            if (++iterations > cap) break;
            auto options = candidates(m, horizon, top);
            bool progressed = false;
            int attempts = 0;
            for (auto [score, r] : options) {
                if (++attempts > (depth <= 2 ? 8 : 4)) break;
                if (r == CH_BUY) { result = 2; break; }
                vector<ll> savedStock = stock;
                ll savedCost = cost;
                size_t savedLength = sequence.size();
                bool ok = gather(r, m, horizon) && fire(r);
                if (ok && stock[m] > savedStock[m]) {
                    progressed = true;
                    break;
                }
                stock = std::move(savedStock);
                cost = savedCost;
                sequence.resize(savedLength);
                if (stopped) break;
            }
            if (result == 2 || !progressed) break;
        }
        if (stock[m] >= needed) result = 1;
        --depth;
        active[m] = 0;
        return result;
    }

    bool gather(int r, int output, ll horizon) {
        const auto& inputs = R[r].in;
        vector<unsigned char> deferred(inputs.size(), 0);
        vector<ll> future(inputs.size());
        ll runs = runsFor(r, output, horizon);
        if (runs >= LINF) return false;
        RequestScope scope(request);
        for (int t = 0; t < (int)inputs.size(); ++t) {
            auto [m, q] = inputs[t];
            ll net = max<ll>(0, q - qtyOf(R[r].out, m));
            __int128 total = (__int128)q + (__int128)(runs - 1) * net;
            future[t] = (ll)min<__int128>(LINF / 2, total);
            scope.raise(m, future[t]);
        }
        for (int guard = 0; guard < 96 + 12 * (int)inputs.size(); ++guard) {
            if (expired()) return false;
            int pick = -1;
            double bestScore = DINF;
            for (int t = 0; t < (int)inputs.size(); ++t) {
                auto [m, q] = inputs[t];
                if (stock[m] >= q || deferred[t]) continue;
                double score = 0;
                if (ordering == 0 || ordering == 1) {
                    auto options = candidates(m, future[t], false);
                    score = options.empty() ? DINF / 4 : options[0].first;
                    score /= max(1.0, (double)future[t] - stock[m]);
                    if (ordering == 1) score = -score;
                    // Postpone buying until coproduct-producing ingredients
                    // have had a chance to provide this molecule for free.
                    if (price[m] >= 0) score += 0.1;
                } else if (ordering == 2) {
                    score = uniform();
                } else {
                    score = t;
                }
                if (pick < 0 || score < bestScore) { pick = t; bestScore = score; }
            }
            if (pick < 0) return true;
            int answer = ensure(inputs[pick].first, inputs[pick].second, future[pick]);
            if (!answer) return false;
            if (answer == 2) deferred[pick] = 1;
        }
        return false;
    }
};

} // namespace adaptive_detail

static Plan adaptiveSolve(ll K, int target, const Plan& incumbent, double deadline) {
    Plan best = incumbent;
    if (best.cost == 0 || nowSec() >= deadline) return best;
    uint64_t seed = 0x9e3779b97f4a7c15ULL ^ (uint64_t)target * 1000003ULL ^ (uint64_t)K ^ ((uint64_t)searchSeed << 17);
    adaptive_detail::Planner prepared(incumbent, deadline, seed);
    int attempt = 0;
    while (nowSec() < deadline) {
        // A runaway recursive policy must leave time for other policies. Large
        // valid incumbents get a longer individual rollout budget.
        double slice = best.cost < LINF ? min(0.35, 0.035 + best.seq.size() * 0.000006) : 0.15;
        adaptive_detail::Planner planner = prepared;
        planner.bound = best.cost;
        planner.preferred = best.choice.size() == (size_t)M ? best.choice : baseChoice;
        planner.deadline = min(deadline, nowSec() + slice);
        planner.randomState = seed + ++attempt * 11995408973635179863ULL;
        planner.credit = (attempt % 5) * 0.25;
        planner.protectAncestors = attempt % 4 != 0;
        planner.noise = attempt <= 5 ? 0 : 0.05 + 0.6 * planner.uniform();
        planner.preference = attempt <= 5 ? 1.0 : (attempt % 2 ? 0.05 + 0.25 * planner.uniform() : 0.55 + 0.9 * planner.uniform());
        planner.ordering = attempt <= 5 ? 0 : planner.random32() % 4;
        if (attempt > 8) {
            for (int m = 0; m < M; ++m) {
                if (planner.uniform() > 0.12 || planner.source[m].empty()) continue;
                planner.preferred[m] = planner.source[m][planner.random32() % planner.source[m].size()];
            }
        }
        // Net-zero target reactions can still leave a purchased catalyst in
        // inventory, a valid operation under the reaction-only output format.
        for (int r = 0; r < (int)R.size(); ++r)
            if (qtyOf(R[r].out, target) >= K &&
                qtyOf(R[r].out, target) <= qtyOf(R[r].in, target) && price[target] >= 0)
                planner.source[target].push_back(r);
        int answer = planner.ensure(target, K, K, true);
        if (answer == 1 && planner.cost < best.cost) {
            ll checked = 0;
            if (verifySeq(planner.sequence, target, K, checked) && checked == planner.cost) {
                best.cost = checked;
                best.seq = std::move(planner.sequence);
                best.proven = false;
                if (best.choice.size() != (size_t)M) best.choice = baseChoice;
                best.dem.clear();
                for (int m = 0; m < M; ++m) if (planner.touched[m]) best.dem.push_back(m);
                if (best.cost == 0) break;
            }
        }
    }
    return best;
}
