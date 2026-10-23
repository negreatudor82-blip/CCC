#ifndef FACTORY_FAMILY_SEARCH_HPP
#define FACTORY_FAMILY_SEARCH_HPP

// Share valid production prefixes and repeatable batches between independent
// cases. Cases still start empty: inventory is never transferred between them.
namespace factory_family {
struct Block { vector<int> seq; ll yield, cost; };
struct State {
    vector<ll> stock;
    ll cost;
    int parent, block, length;
    bool active = true;
};

static bool step(int r, vector<ll>& stock, ll& cost, ll bound) {
    __int128 next = cost;
    for (auto e : R[r].in) {
        ll missing = max(0LL, e.second - stock[e.first]);
        if (missing && price[e.first] < 0) return false;
        if (missing) next += (__int128)missing * price[e.first];
    }
    if (next >= bound || next >= LINF) return false;
    for (auto e : R[r].in) stock[e.first] = max(0LL, stock[e.first] - e.second);
    for (auto e : R[r].out) {
        if ((__int128)stock[e.first] + e.second >= LINF) return false;
        stock[e.first] += e.second;
    }
    cost = (ll)next;
    return true;
}

static vector<Block> blocksFor(const vector<Plan>& bank, int target, ll K, double deadline) {
    map<ll, vector<Block>> byYield;
    for (const Plan& plan : bank) {
        if (nowSec() >= deadline) break;
        if (plan.cost >= LINF) continue;
        vector<ll> stock(M, 0);
        ll cost = 0, largest = 0;
        for (int i = 0; i < (int)plan.seq.size(); ++i) {
            if ((i & 127) == 0 && nowSec() >= deadline) break;
            if (!step(plan.seq[i], stock, cost, LINF)) break;
            if (stock[target] <= largest) continue;
            largest = stock[target];
            Block block{vector<int>(plan.seq.begin(), plan.seq.begin() + i + 1), largest, cost};
            auto& pool = byYield[largest];
            bool duplicate = false;
            for (const auto& old : pool) if (old.seq == block.seq) duplicate = true;
            if (!duplicate) pool.push_back(move(block));
            stable_sort(pool.begin(), pool.end(), [](const Block& a, const Block& b) {
                if (a.cost != b.cost) return a.cost < b.cost;
                return a.seq.size() < b.seq.size();
            });
            if (pool.size() > 6) pool.resize(6);
            if (largest >= K) break;
        }
    }
    vector<Block> blocks;
    for (auto& entry : byYield) for (auto& block : entry.second) blocks.push_back(move(block));
    return blocks;
}

static void improve(Plan& best, const vector<Plan>& bank, int target, ll K, double deadline) {
    if (best.cost <= 0 || best.cost >= LINF || nowSec() >= deadline) return;
    auto blocks = blocksFor(bank, target, K, deadline);
    if (blocks.empty()) return;
    vector<State> states;
    states.push_back({vector<ll>(M, 0), 0, -1, -1, 0, true});
    using Entry = pair<ll, int>;
    priority_queue<Entry, vector<Entry>, greater<Entry>> queue;
    queue.push({0, 0});
    map<ll, vector<int>> labels;
    labels[0].push_back(0);
    const int maxStates = max(300, min(12000, 1600000 / max(1, M)));
    while (!queue.empty() && (int)states.size() < maxStates && nowSec() < deadline) {
        auto [queuedCost, index] = queue.top(); queue.pop();
        if (!states[index].active || queuedCost != states[index].cost || queuedCost >= best.cost) continue;
        // Own copies survive vector growth while new states are appended.
        const vector<ll> startStock = states[index].stock;
        for (int b = 0; b < (int)blocks.size() && nowSec() < deadline; ++b) {
            vector<ll> stock = startStock;
            ll cost = queuedCost;
            int length = 0;
            bool valid = true;
            for (int r : blocks[b].seq) {
                if ((length & 127) == 0 && nowSec() >= deadline) { valid = false; break; }
                if (!step(r, stock, cost, best.cost)) { valid = false; break; }
                ++length;
                if (stock[target] >= K) break;
            }
            if (!valid || stock[target] <= startStock[target]) continue;
            if (stock[target] >= K) {
                vector<pair<int, int>> path{{b, length}};
                for (int p = index; states[p].parent >= 0; p = states[p].parent)
                    path.push_back({states[p].block, states[p].length});
                reverse(path.begin(), path.end());
                vector<int> sequence;
                for (auto entry : path)
                    sequence.insert(sequence.end(), blocks[entry.first].seq.begin(),
                                    blocks[entry.first].seq.begin() + entry.second);
                ll checked;
                if (verifySeq(sequence, target, K, checked) && checked == cost && cost < best.cost) {
                    best.seq = move(sequence); best.cost = cost; best.proven = false;
                }
                continue;
            }
            auto& peers = labels[stock[target]];
            bool dominated = false;
            for (int p : peers) if (states[p].active && states[p].cost <= cost) {
                bool covers = true;
                for (int m = 0; m < M && covers; ++m) covers = states[p].stock[m] >= stock[m];
                if (covers) { dominated = true; break; }
            }
            if (dominated) continue;
            // Bounded labels are a heuristic only. Stock itself is never clipped.
            if (peers.size() >= 40) {
                auto worst = max_element(peers.begin(), peers.end(), [&](int a, int c) {
                    return states[a].cost < states[c].cost;
                });
                if (states[*worst].cost <= cost) continue;
                states[*worst].active = false;
                peers.erase(worst);
            }
            int next = (int)states.size();
            states.push_back({move(stock), cost, index, b, length, true});
            peers.push_back(next); queue.push({cost, next});
            if ((int)states.size() >= maxStates) break;
        }
    }
}
} // namespace factory_family

static void improveCaseFamilies(vector<Plan>& plans, const vector<pair<ll, int>>& cases,
                                const vector<Plan>& previous, double deadline) {
    vector<Plan> bank = previous;
    bank.insert(bank.end(), plans.begin(), plans.end());
    for (int pass = 0; pass < 2 && nowSec() < deadline; ++pass) {
        bool improved = false;
        for (int i = 0; i < (int)cases.size() && nowSec() < deadline; ++i) {
            ll old = plans[i].cost;
            double slice = max(0.002, (deadline - nowSec()) / max(1, (int)cases.size() - i));
            factory_family::improve(plans[i], bank, cases[i].second, cases[i].first,
                                    min(deadline, nowSec() + slice));
            if (plans[i].cost < old) { bank.push_back(plans[i]); improved = true; }
        }
        if (!improved) break;
    }
}

#endif
