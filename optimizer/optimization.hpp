#ifndef FACTORY_OPTIMIZATION_HPP
#define FACTORY_OPTIMIZATION_HPP

// Included after the original simulator. Every heuristic candidate is replayed
// from an empty inventory before it can replace the incumbent.
#include "adaptive_planner.hpp"
#include "sequence_polish.hpp"
#include "trajectory_search.hpp"
#include "family_search.hpp"

static chrono::steady_clock::time_point clockLimit(double deadline) {
    return T0 + chrono::duration_cast<chrono::steady_clock::duration>(chrono::duration<double>(deadline));
}

static void considerSequence(Plan& best, const vector<int>& seq, int target, ll K,
                             const vector<int>& choice) {
    ll cost;
    if (!verifySeq(seq, target, K, cost) || cost >= best.cost) return;
    best.cost = cost;
    best.seq = seq;
    best.choice = choice;
    best.proven = false;
    best.dem.clear();
    vector<char> used(M, false);
    used[target] = true;
    for (int r : seq) for (auto e : R[r].in) used[e.first] = true;
    for (int m = 0; m < M; ++m) if (used[m]) best.dem.push_back(m);
}

// Reverse liveness removes operations whose outputs cannot help the target.
// A replay is essential: removing a byproduct producer can increase purchases.
static void trimSequence(Plan& best, int target, ll K) {
    vector<ll> need(M, 0);
    need[target] = K;
    vector<int> keep;
    for (int i = (int)best.seq.size() - 1; i >= 0; --i) {
        int r = best.seq[i];
        bool useful = false;
        for (auto e : R[r].out) if (need[e.first] > 0) useful = true;
        if (!useful) continue;
        keep.push_back(r);
        for (auto e : R[r].out) need[e.first] = max(0LL, need[e.first] - e.second);
        for (auto e : R[r].in) need[e.first] += e.second;
    }
    reverse(keep.begin(), keep.end());
    considerSequence(best, keep, target, K, best.choice);
}

static void scheduleCounts(const vector<ll>& counts, const vector<ll>& budget,
                           Plan& best, int target, ll K, double deadline, int mode) {
    if (nowSec() >= deadline) return;
    ll count = accumulate(counts.begin(), counts.end(), 0LL);
    if (count <= 0 || count > max<ll>(10000, 12LL * best.seq.size())) return;
    vector<ll> left = counts;
    vector<int> ch = best.choice.empty() ? baseChoice : best.choice;
    inv = budget; // Conceptual prepurchases determine the order; output uses reactions only.
    fill(inProg.begin(), inProg.end(), 0);
    fill(dem.begin(), dem.end(), 0);
    seqBuf.clear();
    ChoiceScope choices(ch);
    curCost = 0; costBound = best.cost; stepCnt = 0;
    stepCap = max<ll>(20000, count * 6); aborted = false;
    trialDeadline = deadline;
    mt19937 gen(712367u + (unsigned)mode * 991u);
    vector<double> tie(R.size());
    for (double& t : tie) t = (double)(gen() % 10000) / 10000.0;
    for (ll guard = 0; guard < count * 3 + 100 && nowSec() < deadline; ++guard) {
        if (inv[target] >= K) break;
        int pick = -1, blocked = -1;
        double bestScore = -1e300, seedCost = 1e300;
        for (int r = 0; r < (int)R.size(); ++r) if (left[r]) {
            bool ready = true;
            double missing = 0, value = 0;
            for (auto e : R[r].in) {
                ll deficit = max(0LL, e.second - inv[e.first]);
                if (deficit) ready = false;
                missing += deficit * max(0.01, wgt(e.first));
            }
            for (auto e : R[r].prod)
                value += e.second * max(0.01, wgt(e.first));
            double score = value / (1.0 + R[r].in.size());
            if (mode % 3 == 1) score = -r;
            if (mode % 3 == 2) score = tie[r] * 1000;
            if (ready && score > bestScore) { bestScore = score; pick = r; }
            missing *= 1.0 + 0.05 * tie[r];
            if (missing < seedCost) { seedCost = missing; blocked = r; }
        }
        if (pick < 0) {
            if (blocked < 0) break;
            // Break a startup dependency using the incumbent's known recipes.
            bool ok = true;
            for (auto e : R[blocked].in) {
                if (inv[e.first] >= e.second) continue;
                if (price[e.first] >= 0) {
                    __int128 extra = (__int128)(e.second - inv[e.first]) * price[e.first];
                    if (extra >= LINF || extra + curCost >= costBound) { ok = false; break; }
                    curCost += (ll)extra;
                    inv[e.first] = e.second;
                } else if (!ensureMol(e.first, e.second, false)) { ok = false; break; }
            }
            if (!ok || aborted) break;
            // A seed can consume another ingredient, so check all inputs again.
            bool ready = true;
            for (auto e : R[blocked].in) if (inv[e.first] < e.second) ready = false;
            if (!ready) continue;
            pick = blocked;
        }
        if (!execReaction(pick)) break;
        --left[pick];
    }
    if (!aborted) considerSequence(best, seqBuf, target, K, ch);
}

static void linearImprove(Plan& best, int target, ll K, double deadline) {
    if (best.cost >= LINF || nowSec() >= deadline) return;
    const int nr = (int)R.size();
    vector<int> buy;
    for (int m = 0; m < M; ++m) if (price[m] >= 0) buy.push_back(m);
    const int nv = nr + (int)buy.size();
    vector<vector<double>> A(M, vector<double>(nv));
    vector<double> rhs(M), obj(nv);
    rhs[target] = -(double)K;
    for (int r = 0; r < nr; ++r) {
        for (auto e : R[r].in) A[e.first][r] += e.second;
        for (auto e : R[r].out) A[e.first][r] -= e.second;
    }
    for (int j = 0; j < (int)buy.size(); ++j) {
        A[buy[j]][nr + j] = -1;
        obj[nr + j] = -(double)price[buy[j]];
    }
    vector<ll> incumbentCounts(nr);
    for (int r : best.seq) ++incumbentCounts[r];
    // The full relaxation and incumbent support are complementary: the latter
    // avoids many unseeded cycles in the full stoichiometric relaxation.
    for (int variant = 0; variant < 3 && nowSec() < deadline; ++variant) {
        auto mat = A;
        auto b = rhs;
        auto c = obj;
        if (variant == 1) {
            for (int r = 0; r < nr; ++r) if (!incumbentCounts[r]) {
                vector<double> row(nv); row[r] = 1;
                mat.push_back(move(row)); b.push_back(0);
            }
        }
        if (variant == 2) {
            for (int r = 0; r < nr; ++r) {
                double fee = 0;
                for (auto e : R[r].in) {
                    if (price[e.first] >= 0 && qtyOf(R[r].out, e.first) >= e.second)
                        fee += price[e.first] * e.second;
                }
                c[r] = -0.01 - fee / max(1LL, incumbentCounts[r]);
            }
        }
        auto result = factory_lp::solveLP(mat, b, c, clockLimit(deadline));
        if (result.status != factory_lp::LPStatus::Optimal) continue;
        vector<ll> counts(nr), budget(M);
        for (int r = 0; r < nr; ++r) {
            if (!isfinite(result.x[r]) || result.x[r] > 1000000) return;
            counts[r] = max(0LL, (ll)ceil(result.x[r] - 1e-7));
        }
        for (int j = 0; j < (int)buy.size(); ++j)
            budget[buy[j]] = max(0LL, (ll)ceil(result.x[nr + j] - 1e-7));
        for (int mode = 0; mode < 3 && nowSec() < deadline; ++mode)
            scheduleCounts(counts, budget, best, target, K, deadline, mode);
    }
}

static vector<vector<int>> policyBank;

#include "integer_search.hpp"

static Plan solveBetter(ll K, int target, double deadline, const Plan& previous) {
    const double start = nowSec();
    const double available = max(0.001, deadline - start);
    trialDeadline = start + available * 0.12;
    Plan best = solveCase(K, target, trialDeadline);
    if (previous.cost < best.cost) {
        vector<int> choice = best.choice.empty() ? baseChoice : best.choice;
        considerSequence(best, previous.seq, target, K, choice);
    }
    bool small = M <= 14 && R.size() <= 18;
    if (small) {
        vector<inventory_search::Reaction> reactions;
        for (const auto& r : R) reactions.push_back({r.in, r.out});
        auto result = inventory_search::search(price, reactions, target, K, best.cost,
                        best.seq, clockLimit(deadline), 200000);
        if (result.cost < best.cost) {
            considerSequence(best, result.sequence, target, K, baseChoice);
        }
        best.proven = result.exhausted && best.cost < LINF;
        return best;
    }
    trialDeadline = start + available * 0.20;
    for (const auto& policy : policyBank) {
        if (nowSec() >= trialDeadline) break;
        ll cost = runTrial(policy, target, K, best.cost, CAP0);
        if (cost < best.cost) considerSequence(best, seqBuf, target, K, policy);
    }
    if (best.cost < LINF && nowSec() < trialDeadline) {
        mt19937 rng(392473u + target);
        descend(best, target, K, trialDeadline, rng);
    }
    if (nowSec() < start + available * 0.28)
        best = adaptiveSolve(K, target, best, start + available * 0.28);
    if (best.cost < LINF && nowSec() < start + available * 0.34)
        linearImprove(best, target, K, start + available * 0.34);
    if (best.cost < LINF && nowSec() < start + available * 0.56)
        integerImprove(best, target, K, start + available * 0.56);
    if (best.cost < LINF && nowSec() < start + available * 0.90)
        trajectoryImprove(best, target, K, start + available * 0.90);
    if (best.cost < LINF && nowSec() < deadline)
        polishSequence(best, target, K, deadline);
    if (best.cost < LINF) trimSequence(best, target, K);
    if (best.choice.size() == (size_t)M) {
        if (find(policyBank.begin(), policyBank.end(), best.choice) == policyBank.end())
            policyBank.push_back(best.choice);
        if (policyBank.size() > 16) policyBank.erase(policyBank.begin());
    }
    trialDeadline = 1e100;
    return best;
}

#endif
