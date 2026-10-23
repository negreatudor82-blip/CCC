#ifndef FACTORY_INTEGER_SEARCH_HPP
#define FACTORY_INTEGER_SEARCH_HPP

// Included after optimization.hpp's scheduling helpers. Linear models guide
// the search only; every accepted sequence is replayed independently.
namespace factory_integer_detail {

static vector<vector<ll>> conservedAtoms() {
    vector<map<string, ll>> formulas(M);
    set<string> elements;
    for (int m = 0; m < M; ++m) {
        const string& name = names[m];
        for (size_t i = 0; i < name.size();) {
            if (!isupper((unsigned char)name[i])) return {};
            string e(1, name[i++]);
            while (i < name.size() && islower((unsigned char)name[i])) e += name[i++];
            ll q = 0;
            while (i < name.size() && isdigit((unsigned char)name[i])) {
                if (q > 100000000) return {};
                q = q * 10 + name[i++] - '0';
            }
            if (q == 0) q = 1;
            formulas[m][e] += q;
            elements.insert(e);
        }
    }
    vector<vector<ll>> result;
    for (const string& element : elements) {
        vector<ll> atom(M);
        for (int m = 0; m < M; ++m) {
            auto it = formulas[m].find(element);
            if (it != formulas[m].end()) atom[m] = it->second;
        }
        bool conserved = true;
        for (const auto& reaction : R) {
            __int128 balance = 0;
            for (auto [m, q] : reaction.out) balance += (__int128)atom[m] * q;
            for (auto [m, q] : reaction.in) balance -= (__int128)atom[m] * q;
            if (balance != 0) { conserved = false; break; }
        }
        if (conserved) result.push_back(move(atom));
    }
    return result;
}

struct CountScheduler {
    vector<ll> inventory, left;
    vector<int> sequence, guide, cursor;
    vector<vector<int>> positions;
    vector<unsigned char> active;
    vector<double> priority;
    ll cap;
    double deadline;
    bool failed = false;
    int calls = 0;

    CountScheduler(const vector<ll>& counts, const vector<ll>& budget,
                   const Plan& incumbent, const vector<ll>& initialStock,
                   const vector<int>& prefix, double until, unsigned seed)
        : inventory(budget), left(counts), guide(M, CH_NONE), cursor(R.size(), 0),
          positions(R.size()), active(M, 0), priority(R.size()), deadline(until) {
        cap = min<ll>(100000, max<ll>(10000, incumbent.seq.size() * 8 + 1000));
        sequence = prefix;
        for (int m = 0; m < M; ++m) {
            __int128 total = (__int128)inventory[m] + initialStock[m];
            if (total >= LINF) { failed = true; break; }
            inventory[m] = (ll)total;
        }
        for (int m = 0; m < M; ++m) if (price[m] >= 0) guide[m] = CH_BUY;
        for (int i = 0; i < (int)incumbent.seq.size(); ++i) {
            int r = incumbent.seq[i];
            positions[r].push_back(i);
            for (auto [m, q] : R[r].prod) if (guide[m] == CH_NONE) guide[m] = r;
        }
        for (int m = 0; m < M; ++m) if (guide[m] == CH_NONE) guide[m] = baseChoice[m];
        mt19937 gen(seed);
        for (int r : prefix) if (cursor[r] < (int)positions[r].size()) ++cursor[r];
        for (double& value : priority) value = (gen() % 1000000) * 1e-6;
    }

    bool expired() {
        if (failed) return true;
        if (++calls > cap * 30 || (calls % 32 == 0 && nowSec() >= deadline)) failed = true;
        return failed;
    }

    bool ready(int r) const {
        for (auto [m, q] : R[r].in) if (inventory[m] < q) return false;
        return true;
    }

    bool execute(int r) {
        if (expired() || (ll)sequence.size() >= cap || !ready(r)) return false;
        for (auto [m, q] : R[r].in) inventory[m] -= q;
        for (auto [m, q] : R[r].out) {
            if ((__int128)inventory[m] + q >= LINF) return false;
            inventory[m] += q;
        }
        sequence.push_back(r);
        if (left[r] > 0) --left[r];
        if (cursor[r] < (int)positions[r].size()) ++cursor[r];
        return true;
    }

    int nextPosition(int r) const {
        return cursor[r] < (int)positions[r].size()
            ? positions[r][cursor[r]] : 1000000000 + r;
    }

    bool ensure(int m, ll quantity, int depth = 0) {
        if (inventory[m] >= quantity) return true;
        if (expired() || depth > M + 2) return false;
        if (price[m] >= 0) { inventory[m] = quantity; return true; }
        if (active[m] || guide[m] < 0) return false;
        active[m] = 1;
        int r = guide[m];
        bool ok = true;
        while (inventory[m] < quantity && ok) {
            ll before = inventory[m];
            // Producing one ingredient can consume a previously gathered one.
            for (int retry = 0; retry < M + 8 && !ready(r); ++retry) {
                for (auto [a, q] : R[r].in) {
                    if (!ensure(a, q, depth + 1)) { ok = false; break; }
                }
                if (!ok || expired()) { ok = false; break; }
            }
            if (!ok || !execute(r) || inventory[m] <= before) { ok = false; break; }
        }
        active[m] = 0;
        return ok && inventory[m] >= quantity;
    }

    bool run(int target, ll quantity, int mode) {
        for (ll iteration = 0; iteration < cap && !expired(); ++iteration) {
            if (inventory[target] >= quantity) return true;
            int pick = -1, blocked = -1;
            double score = -1e300;
            vector<ll> demand(M);
            if (mode == 2) for (int r = 0; r < (int)R.size(); ++r) if (left[r])
                for (auto [m, q] : R[r].in) {
                    __int128 next = (__int128)demand[m] + (__int128)q * left[r];
                    demand[m] = (ll)min<__int128>(LINF / 2, next);
                }
            for (int r = 0; r < (int)R.size(); ++r) if (left[r] > 0) {
                if (blocked < 0 || nextPosition(r) < nextPosition(blocked)) blocked = r;
                if (!ready(r)) continue;
                double value = -nextPosition(r);
                if (mode == 1) value = priority[r] * 10000;
                if (mode == 2) {
                    value = 0;
                    for (auto [m, q] : R[r].out)
                        value += min(q, max(0LL, demand[m] - inventory[m])) * max(0.01, wgt(m));
                }
                if (pick < 0 || value > score) { pick = r; score = value; }
            }
            if (pick >= 0) {
                if (!execute(pick)) return false;
            } else {
                if (blocked < 0) return ensure(target, quantity);
                int ingredient = -1;
                double cheapest = DINF;
                for (auto [m, q] : R[blocked].in) if (inventory[m] < q) {
                    double value = (q - inventory[m]) * max(0.01, wgt(m));
                    if (ingredient < 0 || value < cheapest) { ingredient = m; cheapest = value; }
                }
                if (ingredient < 0 || !ensure(ingredient, qtyOf(R[blocked].in, ingredient))) return false;
            }
        }
        return false;
    }
};

static void tryCounts(const vector<double>& x, const vector<int>& buy, Plan& best,
                      int target, ll quantity, const vector<ll>& initialStock,
                      const vector<int>& prefix, double deadline, unsigned seed) {
    const int nr = R.size();
    vector<ll> counts(nr), budget(M);
    ll total = 0;
    for (int r = 0; r < nr; ++r) {
        if (!isfinite(x[r]) || x[r] > 1000000) return;
        counts[r] = max(0LL, (ll)ceil(x[r] - 1e-7));
        total += counts[r];
    }
    if (total <= 0 || total > max<ll>(10000, best.seq.size() * 10)) return;
    for (int j = 0; j < (int)buy.size(); ++j) {
        double amount = x[nr + j];
        if (!isfinite(amount) || amount > 1e12) return;
        budget[buy[j]] = max(0LL, (ll)ceil(amount - 1e-7));
    }
    for (int mode = 0; mode < 3 && nowSec() < deadline; ++mode) {
        CountScheduler scheduler(counts, budget, best, initialStock, prefix, deadline, seed + mode * 17389u);
        if (!scheduler.run(target, quantity, mode)) continue;
        ll actual;
        if (verifySeq(scheduler.sequence, target, quantity, actual) && actual < best.cost)
            considerSequence(best, scheduler.sequence, target, quantity,
                             best.choice.empty() ? baseChoice : best.choice);
    }
}

} // namespace factory_integer_detail

static void integerImprove(Plan& best, int target, ll K, double deadline) {
    using namespace factory_integer_detail;
    if (best.cost >= LINF || best.cost <= 0 || nowSec() >= deadline) return;
    const int nr = R.size();
    vector<int> buy;
    for (int m = 0; m < M; ++m) if (price[m] >= 0) buy.push_back(m);
    const int nv = nr + buy.size();
    vector<vector<double>> basic(M, vector<double>(nv));
    vector<double> rhs(M), objective(nv);
    rhs[target] = -(double)K;
    for (int r = 0; r < nr; ++r) {
        for (auto [m, q] : R[r].in) basic[m][r] += q;
        for (auto [m, q] : R[r].out) basic[m][r] -= q;
    }
    for (int j = 0; j < (int)buy.size(); ++j) {
        basic[buy[j]][nr + j] = -1;
        objective[nr + j] = -(double)price[buy[j]];
    }
    const auto atoms = conservedAtoms();
    mt19937 rng(searchSeed ^ (unsigned)target * 9781u ^ (unsigned)K * 1723u);
    // Retained prefixes provide real catalyst stock and support correlated
    // recipe substitutions in the remainder without paying startup twice.
    const vector<int> original = best.seq;
    vector<int> cuts = {(int)original.size() / 4, (int)original.size() / 2,
                        0, (int)original.size() * 3 / 4,
                        (int)original.size() / 8, (int)original.size() * 5 / 8};
    for (int neighborhood = 0; neighborhood < (int)cuts.size() && nowSec() < deadline; ++neighborhood) {
        const int cut = cuts[neighborhood];
        const double localDeadline = min(deadline, nowSec() +
            (deadline - nowSec()) / max(1, (int)cuts.size() - neighborhood));
        vector<int> prefix(original.begin(), original.begin() + cut);
        vector<ll> stock(M), incumbentCounts(nr);
        ll spent = 0;
        bool valid = true;
        for (int r : prefix) {
            for (auto [m, q] : R[r].in) {
                ll shortage = max(0LL, q - stock[m]);
                if (shortage && price[m] < 0) { valid = false; break; }
                __int128 nextCost = (__int128)spent + (__int128)shortage * max(0LL, price[m]);
                if (nextCost >= best.cost) { valid = false; break; }
                spent = (ll)nextCost;
                stock[m] = max(0LL, stock[m] - q);
            }
            if (!valid) break;
            for (auto [m, q] : R[r].out) {
                if ((__int128)stock[m] + q >= LINF) { valid = false; break; }
                stock[m] += q;
            }
            if (!valid) break;
        }
        if (!valid) continue;
        if (stock[target] >= K) {
            considerSequence(best, prefix, target, K, best.choice);
            continue;
        }
        for (int i = cut; i < (int)original.size(); ++i) ++incumbentCounts[original[i]];
        auto matrix = basic;
        auto b = rhs;
        for (int m = 0; m < M; ++m) b[m] += stock[m];
        // Restrict to the incumbent support initially. Other neighborhoods
        // introduce alternate reactions while bounding search resources.
        for (int r = 0; r < nr; ++r) {
            bool allowed = incumbentCounts[r] > 0 || neighborhood == 5;
            if (neighborhood == 1 && !allowed) allowed = rng() % 4 == 0;
            if (neighborhood == 3 && !allowed) allowed = rng() % 2 == 0;
            if (neighborhood == 4 && !allowed) allowed = rng() % 3 == 0;
            vector<double> row(nv); row[r] = 1;
            matrix.push_back(move(row));
            b.push_back(allowed ? max<ll>(12, incumbentCounts[r] * 2 + 6) : 0);
        }
        // Conserved-atom stock cuts bound seed needs for the incumbent
        // neighborhood, with already-purchased atoms credited exactly.
        if (neighborhood == 0 || neighborhood == 2) for (const auto& atom : atoms) {
            ll required = 0;
            for (int r = 0; r < nr; ++r) if (incumbentCounts[r]) {
                __int128 need = 0;
                for (auto [m, q] : R[r].in) need += (__int128)q * atom[m];
                if (need >= LINF) { required = 0; break; }
                required = max(required, (ll)need);
            }
            __int128 available = 0;
            for (int m = 0; m < M; ++m) available += (__int128)atom[m] * stock[m];
            if (available >= required || required <= 0) continue;
            vector<double> row(nv);
            for (int j = 0; j < (int)buy.size(); ++j) row[nr + j] = -(double)atom[buy[j]];
            matrix.push_back(move(row)); b.push_back(-(double)(required - available));
        }
        struct Branch { int variable; double bound; bool lower; };
        vector<vector<Branch>> stack(1);
        int nodes = 0;
        while (!stack.empty() && ++nodes <= 18 && nowSec() < localDeadline) {
            auto branches = move(stack.back()); stack.pop_back();
            auto a = matrix; auto bb = b;
            for (const auto& branch : branches) {
                vector<double> row(nv);
                row[branch.variable] = branch.lower ? -1 : 1;
                a.push_back(move(row)); bb.push_back(branch.lower ? -branch.bound : branch.bound);
            }
            auto lp = factory_lp::solveLP(a, bb, objective, clockLimit(localDeadline));
            if (lp.status != factory_lp::LPStatus::Optimal) continue;
            if (-lp.objective + spent >= best.cost - 1e-7) continue;
            tryCounts(lp.x, buy, best, target, K, stock, prefix, localDeadline, rng());
            int fractional = -1;
            double fractionality = 0;
            for (int r = 0; r < nr; ++r) {
                double f = lp.x[r] - floor(lp.x[r]);
                double score = min(f, 1 - f);
                if (score > fractionality + 1e-6) { fractionality = score; fractional = r; }
            }
            if (fractional < 0) continue;
            auto down = branches, up = branches;
            down.push_back({fractional, floor(lp.x[fractional]), false});
            up.push_back({fractional, ceil(lp.x[fractional]), true});
            // Explore the direction nearer the incumbent first.
            if (incumbentCounts[fractional] > lp.x[fractional]) {
                stack.push_back(move(down)); stack.push_back(move(up));
            } else {
                stack.push_back(move(up)); stack.push_back(move(down));
            }
        }
    }
}

#endif
