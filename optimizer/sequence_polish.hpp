#ifndef FACTORY_SEQUENCE_POLISH_HPP
#define FACTORY_SEQUENCE_POLISH_HPP

// Included after the original globals and verifySeq. This optimizer changes
// only the incumbent; every accepted edit is replayed from an empty inventory.
namespace factory_sequence_polish {

struct Prefix {
    vector<vector<ll>> inventory;
    vector<ll> cost;
};

static bool step(int r, vector<ll>& inventory, ll& cost, ll bound) {
    for (auto e : R[r].in) {
        ll missing = max(0LL, e.second - inventory[e.first]);
        if (missing && price[e.first] < 0) return false;
        if (missing) {
            __int128 next = (__int128)cost + (__int128)missing * price[e.first];
            if (next >= bound || next > numeric_limits<ll>::max()) return false;
            cost = (ll)next;
        }
        inventory[e.first] = max(0LL, inventory[e.first] - e.second);
    }
    for (auto e : R[r].out) {
        __int128 next = (__int128)inventory[e.first] + e.second;
        if (next > numeric_limits<ll>::max()) return false;
        inventory[e.first] = (ll)next;
    }
    return true;
}

static bool makePrefix(const Plan& best, Prefix& prefix, double deadline) {
    prefix.inventory.clear();
    prefix.cost.clear();
    prefix.inventory.reserve(best.seq.size() + 1);
    prefix.cost.reserve(best.seq.size() + 1);
    prefix.inventory.push_back(vector<ll>(M, 0));
    prefix.cost.push_back(0);
    vector<ll> inventory(M, 0);
    ll cost = 0;
    for (size_t i = 0; i < best.seq.size(); ++i) {
        if ((i & 63) == 0 && nowSec() >= deadline) return false;
        if (!step(best.seq[i], inventory, cost, LINF)) return false;
        prefix.inventory.push_back(inventory);
        prefix.cost.push_back(cost);
    }
    return cost == best.cost;
}

static bool accept(Plan& best, vector<int> candidate, size_t unchanged,
                   const Prefix& prefix, int target, ll K, double deadline) {
    if (nowSec() >= deadline || unchanged >= prefix.inventory.size()) return false;
    vector<ll> inventory = prefix.inventory[unchanged];
    ll cost = prefix.cost[unchanged];
    if (cost >= best.cost) return false;
    if (inventory[target] >= K) candidate.resize(unchanged);
    else {
        for (size_t i = unchanged; i < candidate.size(); ++i) {
            if ((i & 31) == 0 && nowSec() >= deadline) return false;
            if (!step(candidate[i], inventory, cost, best.cost)) return false;
            if (inventory[target] >= K) {
                candidate.resize(i + 1);
                break;
            }
        }
    }
    if (inventory[target] < K || cost >= best.cost || nowSec() >= deadline) return false;
    ll verified = 0;
    if (!verifySeq(candidate, target, K, verified) || verified != cost || verified >= best.cost)
        return false;
    best.cost = verified;
    best.seq = move(candidate);
    best.proven = false;
    best.dem.clear();
    vector<unsigned char> demanded(M, 0);
    demanded[target] = 1;
    for (int r : best.seq) for (auto e : R[r].in) demanded[e.first] = 1;
    for (int m = 0; m < M; ++m) if (demanded[m]) best.dem.push_back(m);
    return true;
}

// Removing a consumer can make a complete upstream branch redundant. Reverse
// liveness collects those branches, including gross catalyst requirements.
static vector<int> liveSequence(const vector<int>& sequence, int target, ll K) {
    vector<__int128> need(M, 0);
    need[target] = K;
    vector<int> kept;
    for (int i = (int)sequence.size() - 1; i >= 0; --i) {
        int r = sequence[i];
        bool useful = false;
        for (auto e : R[r].out) if (need[e.first] > 0) useful = true;
        if (!useful) continue;
        kept.push_back(r);
        for (auto e : R[r].out) need[e.first] = max((__int128)0, need[e.first] - e.second);
        for (auto e : R[r].in) need[e.first] += e.second;
    }
    reverse(kept.begin(), kept.end());
    return kept;
}

// Pull a later producer, and any nonbuyable prerequisites it needs, before an
// earlier purchase. This handles indirect byproducts and seed recycling that
// adjacent swaps cannot cross. The remaining actions preserve their order.
static bool reorder(Plan& best, const Prefix& prefix, int target, ll K, double deadline) {
    struct Event { ll value; int position, molecule; };
    vector<Event> events;
    const int n = (int)best.seq.size();
    for (int i = 0; i < n; ++i) {
        if ((i & 63) == 0 && nowSec() >= deadline) return false;
        for (auto e : R[best.seq[i]].in) {
            if (price[e.first] <= 0) continue;
            ll missing = max(0LL, e.second - prefix.inventory[i][e.first]);
            if (missing) events.push_back({(ll)((__int128)missing * price[e.first]), i, e.first});
        }
    }
    sort(events.begin(), events.end(), [](const Event& a, const Event& b) {
        if (a.value != b.value) return a.value > b.value;
        return a.position < b.position;
    });
    for (const Event& event : events) {
        if (nowSec() >= deadline) return false;
        const int begin = event.position;
        int examined = 0;
        for (int end = begin + 1; end < n; ++end) {
            if ((end & 31) == 0 && nowSec() >= deadline) return false;
            bool produces = false;
            for (auto e : R[best.seq[end]].prod)
                if (e.first == event.molecule && e.second > 0) produces = true;
            if (!produces) continue;
            // Several nearby copies provide diversity without quadratic runs
            // of the same producer for a single purchase event.
            if (++examined > 24) break;
            vector<unsigned char> take(end - begin + 1, 0);
            vector<__int128> need(M, 0);
            auto include = [&](int position) {
                take[position - begin] = 1;
                const Reaction& reaction = R[best.seq[position]];
                for (auto e : reaction.out) if (price[e.first] < 0)
                    need[e.first] = max((__int128)0, need[e.first] - e.second);
                for (auto e : reaction.in) if (price[e.first] < 0) need[e.first] += e.second;
            };
            include(end);
            for (int j = end - 1; j >= begin; --j) {
                bool helps = false;
                for (auto e : R[best.seq[j]].out)
                    if (price[e.first] < 0 && need[e.first] > prefix.inventory[begin][e.first])
                        helps = true;
                if (helps) include(j);
            }
            bool feasible = true;
            for (int m = 0; m < M; ++m)
                if (need[m] > prefix.inventory[begin][m]) feasible = false;
            if (!feasible) continue;
            int selected = accumulate(take.begin(), take.end(), 0);
            if (selected == end - begin + 1) continue;
            vector<int> candidate;
            candidate.reserve(n);
            candidate.insert(candidate.end(), best.seq.begin(), best.seq.begin() + begin);
            for (int j = begin; j <= end; ++j) if (take[j - begin]) candidate.push_back(best.seq[j]);
            for (int j = begin; j <= end; ++j) if (!take[j - begin]) candidate.push_back(best.seq[j]);
            candidate.insert(candidate.end(), best.seq.begin() + end + 1, best.seq.end());
            if (accept(best, move(candidate), begin, prefix, target, K, deadline)) return true;
        }
    }
    return false;
}

static bool deleteGroups(Plan& best, const Prefix& prefix, int target, ll K, double deadline) {
    const int n = (int)best.seq.size();
    vector<vector<int>> occurrences(R.size());
    for (int i = 0; i < n; ++i) occurrences[best.seq[i]].push_back(i);
    vector<int> order;
    for (int r = 0; r < (int)R.size(); ++r) if (!occurrences[r].empty()) order.push_back(r);
    sort(order.begin(), order.end(), [&](int a, int b) {
        ll ca = 0, cb = 0;
        for (int i : occurrences[a]) ca += prefix.cost[i + 1] - prefix.cost[i];
        for (int i : occurrences[b]) cb += prefix.cost[i + 1] - prefix.cost[i];
        if (ca != cb) return ca > cb;
        return occurrences[a].size() > occurrences[b].size();
    });
    for (int r : order) {
        if (nowSec() >= deadline) return false;
        const auto& positions = occurrences[r];
        for (int mode = 0; mode < 3; ++mode) {
            if (mode && positions.size() < 2) continue;
            vector<unsigned char> remove(n, 0);
            size_t lo = mode == 2 ? positions.size() / 2 : 0;
            size_t hi = mode == 1 ? (positions.size() + 1) / 2 : positions.size();
            for (size_t j = lo; j < hi; ++j) remove[positions[j]] = 1;
            vector<int> candidate;
            candidate.reserve(n);
            for (int i = 0; i < n; ++i) if (!remove[i]) candidate.push_back(best.seq[i]);
            if (accept(best, candidate, positions[lo], prefix, target, K, deadline)) return true;
            if (nowSec() >= deadline) return false;
            auto live = liveSequence(candidate, target, K);
            if (live.size() < candidate.size() && accept(best, move(live), 0, prefix, target, K, deadline))
                return true;
        }
    }
    // Chunk deletions remove unused multistep branches; the final pass tests
    // every single occurrence rather than fixing one count for each recipe.
    for (int width = max(1, n / 2); ; width = max(1, width / 2)) {
        for (int begin = 0; begin < n; begin += width) {
            if (nowSec() >= deadline) return false;
            int end = min(n, begin + width);
            vector<int> candidate;
            candidate.reserve(n - (end - begin));
            candidate.insert(candidate.end(), best.seq.begin(), best.seq.begin() + begin);
            candidate.insert(candidate.end(), best.seq.begin() + end, best.seq.end());
            if (accept(best, move(candidate), begin, prefix, target, K, deadline)) return true;
        }
        if (width == 1) break;
    }
    return false;
}

}  // namespace factory_sequence_polish

static void polishSequence(Plan& best, int target, ll K, double deadline) {
    if (best.cost >= LINF || best.cost <= 0 || best.seq.empty()) return;
    while (nowSec() < deadline) {
        factory_sequence_polish::Prefix prefix;
        if (!factory_sequence_polish::makePrefix(best, prefix, deadline)) return;
        auto live = factory_sequence_polish::liveSequence(best.seq, target, K);
        if (live.size() < best.seq.size() && factory_sequence_polish::accept(
                best, move(live), 0, prefix, target, K, deadline)) continue;
        if (factory_sequence_polish::reorder(best, prefix, target, K, deadline)) continue;
        if (factory_sequence_polish::deleteGroups(best, prefix, target, K, deadline)) continue;
        break;
    }
}

#endif  // FACTORY_SEQUENCE_POLISH_HPP
