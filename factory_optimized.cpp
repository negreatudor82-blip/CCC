// Molecule Factory, Level 3 -- standalone C++17 optimizer.
// Build: g++ -O2 -std=c++17 factory_optimized.cpp -o factory
// Run: ./factory INPUT [OUTPUT] [SECONDS=60] [SEED=987654321]
// A matching existing OUTPUT is validated and used as an incumbent.
// Small networks use exact inventory search; larger networks use heuristics.
// Every accepted solution is independently replayed from empty inventory.
// Larger instances report best-found costs, not unproved global optima.

#include <bits/stdc++.h>
// ---- linear_program.hpp ----
#ifndef FACTORY_LINEAR_PROGRAM_HPP
#define FACTORY_LINEAR_PROGRAM_HPP

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <limits>
#include <numeric>
#include <vector>

// A self-contained two-phase simplex solver.  The caller must still validate
// integer plans; a floating-point LP result is not an exact certificate.
namespace factory_lp {

enum class LPStatus { Optimal, Infeasible, Unbounded, TimedOut, NumericalFailure };

struct LPSolution {
    LPStatus status = LPStatus::NumericalFailure;
    double objective = 0;
    std::vector<double> x;
    std::size_t pivots = 0;
};

namespace detail {
using Real = long double;
using Clock = std::chrono::steady_clock;

class Simplex {
    static constexpr Real eps = 1e-11L;
    int m, n;
    std::vector<int> basic, nonbasic;
    std::vector<std::vector<Real>> tab;
    Clock::time_point deadline;
    std::size_t pivots = 0;

    bool expired() const { return Clock::now() >= deadline; }

    bool pivot(int r, int s) {
        if (expired()) return false;
        const Real p = tab[r][s];
        if (!std::isfinite(p) || std::fabs(p) < 1e-25L) return false;
        const Real inverse = 1 / p;
        for (int i = 0; i < m + 2; ++i) {
            if (i == r) continue;
            const Real factor = tab[i][s] * inverse;
            if (factor != 0) {
                for (int j = 0; j < n + 2; ++j)
                    if (j != s) tab[i][j] -= tab[r][j] * factor;
            }
        }
        for (int j = 0; j < n + 2; ++j)
            if (j != s) tab[r][j] *= inverse;
        for (int i = 0; i < m + 2; ++i)
            if (i != r) tab[i][s] *= -inverse;
        tab[r][s] = inverse;
        std::swap(basic[r], nonbasic[s]);
        ++pivots;
        return true;
    }

    LPStatus simplex(int phase) {
        const int objective = phase == 1 ? m + 1 : m;
        while (true) {
            if (expired()) return LPStatus::TimedOut;
            // Bland's rule avoids cycling on degenerate zero-RHS networks.
            int entering = -1;
            for (int j = 0; j <= n; ++j) {
                if (phase == 2 && nonbasic[j] == -1) continue;
                if (tab[objective][j] < -eps &&
                    (entering < 0 || nonbasic[j] < nonbasic[entering]))
                    entering = j;
            }
            if (entering < 0) return LPStatus::Optimal;
            int leaving = -1;
            for (int i = 0; i < m; ++i) {
                if (tab[i][entering] <= eps) continue;
                if (leaving < 0) {
                    leaving = i;
                    continue;
                }
                Real lhs = tab[i][n + 1] / tab[i][entering];
                Real rhs = tab[leaving][n + 1] / tab[leaving][entering];
                Real tolerance = eps * (1 + std::fabs(lhs) + std::fabs(rhs));
                if (lhs < rhs - tolerance ||
                    (std::fabs(lhs - rhs) <= tolerance && basic[i] < basic[leaving]))
                    leaving = i;
            }
            if (leaving < 0) return LPStatus::Unbounded;
            if (!pivot(leaving, entering))
                return expired() ? LPStatus::TimedOut : LPStatus::NumericalFailure;
        }
    }

public:
    Simplex(const std::vector<std::vector<double>>& a,
            const std::vector<double>& b, const std::vector<double>& c,
            Clock::time_point limit)
        : m(static_cast<int>(b.size())), n(static_cast<int>(c.size())),
          basic(m), nonbasic(n + 1),
          tab(m + 2, std::vector<Real>(n + 2)), deadline(limit) {
        for (int i = 0; i < m; ++i) {
            for (int j = 0; j < n; ++j) tab[i][j] = a[i][j];
            basic[i] = n + i;
            tab[i][n] = -1;
            tab[i][n + 1] = b[i];
        }
        for (int j = 0; j < n; ++j) {
            nonbasic[j] = j;
            tab[m][j] = -c[j];
        }
        nonbasic[n] = -1;
        tab[m + 1][n] = 1;
    }

    LPSolution solve() {
        LPSolution answer;
        if (expired()) {
            answer.status = LPStatus::TimedOut;
            return answer;
        }
        int row = -1;
        for (int i = 0; i < m; ++i)
            if (row < 0 || tab[i][n + 1] < tab[row][n + 1]) row = i;
        if (row >= 0 && tab[row][n + 1] < -eps) {
            if (!pivot(row, n)) {
                answer.status = expired() ? LPStatus::TimedOut : LPStatus::NumericalFailure;
                return answer;
            }
            LPStatus phase1 = simplex(1);
            if (phase1 != LPStatus::Optimal) {
                answer.status = phase1 == LPStatus::Unbounded
                    ? LPStatus::NumericalFailure : phase1;
                answer.pivots = pivots;
                return answer;
            }
            if (tab[m + 1][n + 1] < -eps) {
                answer.status = LPStatus::Infeasible;
                answer.pivots = pivots;
                return answer;
            }
            if (std::fabs(tab[m + 1][n + 1]) > 1e-7L) {
                answer.status = LPStatus::NumericalFailure;
                answer.pivots = pivots;
                return answer;
            }
            // Remove an artificial basic variable if its row is not redundant.
            for (int i = 0; i < m; ++i) {
                if (basic[i] != -1) continue;
                int column = -1;
                for (int j = 0; j <= n; ++j)
                    if (nonbasic[j] != -1 && std::fabs(tab[i][j]) > eps &&
                        (column < 0 || nonbasic[j] < nonbasic[column])) column = j;
                if (column >= 0 && !pivot(i, column)) {
                    answer.status = expired() ? LPStatus::TimedOut : LPStatus::NumericalFailure;
                    answer.pivots = pivots;
                    return answer;
                }
            }
        }
        answer.status = simplex(2);
        answer.pivots = pivots;
        if (answer.status != LPStatus::Optimal) return answer;
        answer.x.assign(n, 0);
        for (int i = 0; i < m; ++i)
            if (basic[i] >= 0 && basic[i] < n)
                answer.x[basic[i]] = static_cast<double>(tab[i][n + 1]);
        answer.objective = static_cast<double>(tab[m][n + 1]);
        return answer;
    }
};
} // namespace detail

// Maximize c*x subject to A*x <= b and x >= 0.
// TimedOut results intentionally have no partial solution or usable bound.
inline LPSolution solveLP(
    const std::vector<std::vector<double>>& a, const std::vector<double>& b,
    const std::vector<double>& c,
    std::chrono::steady_clock::time_point deadline =
        std::chrono::steady_clock::time_point::max()) {
    LPSolution failure;
    if (a.size() != b.size()) return failure;
    for (std::size_t i = 0; i < a.size(); ++i) {
        if (a[i].size() != c.size() || !std::isfinite(b[i])) return failure;
        for (double value : a[i]) if (!std::isfinite(value)) return failure;
    }
    for (double value : c) if (!std::isfinite(value)) return failure;
    detail::Simplex solver(a, b, c, deadline);
    LPSolution answer = solver.solve();
    if (answer.status != LPStatus::Optimal) return answer;
    for (double& value : answer.x) {
        if (!std::isfinite(value) || value < -1e-7) {
            answer.status = LPStatus::NumericalFailure;
            return answer;
        }
        if (value < 0) value = 0;
    }
    for (std::size_t i = 0; i < a.size(); ++i) {
        long double lhs = 0, magnitude = std::fabs(b[i]);
        for (std::size_t j = 0; j < c.size(); ++j) {
            long double term = static_cast<long double>(a[i][j]) * answer.x[j];
            lhs += term;
            magnitude += std::fabs(term);
        }
        if (lhs > b[i] + 1e-7L * (1 + magnitude)) {
            answer.status = LPStatus::NumericalFailure;
            return answer;
        }
    }
    long double objective = 0;
    for (std::size_t j = 0; j < c.size(); ++j)
        objective += static_cast<long double>(c[j]) * answer.x[j];
    if (!std::isfinite(objective)) answer.status = LPStatus::NumericalFailure;
    else answer.objective = static_cast<double>(objective);
    return answer;
}

} // namespace factory_lp
#endif

// ---- inventory_search.hpp ----
#ifndef MOLECULE_FACTORY_INVENTORY_SEARCH_HPP
#define MOLECULE_FACTORY_INVENTORY_SEARCH_HPP

#include <algorithm>
#include <chrono>
#include <cstddef>
#include <functional>
#include <limits>
#include <queue>
#include <unordered_map>
#include <utility>
#include <vector>

// Exact uniform-cost search for small networks. Purchases are made only to
// cover the gross inputs of a reaction. Inventories are never clipped: hitting
// either resource limit returns the incumbent without claiming optimality.
namespace inventory_search {

using ll = long long;

struct Reaction {
    std::vector<std::pair<int, ll>> inputs, outputs;
};

struct Result {
    ll cost;
    std::vector<int> sequence;
    bool exhausted;
};

namespace detail {
struct InventoryHash {
    std::size_t operator()(const std::vector<ll>& inventory) const noexcept {
        std::size_t h = inventory.size();
        for (ll amount : inventory) {
            const std::size_t x = std::hash<ll>{}(amount);
            h ^= x + std::size_t(0x9e3779b97f4a7c15ULL) + (h << 6) + (h >> 2);
        }
        return h;
    }
};

struct Node {
    std::vector<ll> inventory;
    ll cost;
    std::size_t parent;
    int reaction;
};

struct CompactReaction {
    int original;
    std::vector<std::pair<int, ll>> inputs, outputs;
};
}  // namespace detail

inline Result search(
    const std::vector<ll>& prices,
    const std::vector<Reaction>& reactions,
    int target,
    ll needed,
    ll incumbentCost,
    const std::vector<int>& incumbentSequence,
    std::chrono::steady_clock::time_point deadline,
    std::size_t maxStates = 100000) {

    Result answer{incumbentCost, incumbentSequence, false};
    const std::size_t noParent = std::numeric_limits<std::size_t>::max();
    const ll limit = std::numeric_limits<ll>::max();
    const int moleculeCount = static_cast<int>(prices.size());
    if (target < 0 || target >= moleculeCount) return answer;
    if (needed <= 0) return Result{0, {}, true};
    if (incumbentCost == 0 && !incumbentSequence.empty()) {
        answer.exhausted = true;
        return answer;
    }
    if (maxStates == 0 || std::chrono::steady_clock::now() >= deadline)
        return answer;

    // Backward closure is safe even with catalysts and cycles: whenever a
    // reaction may produce a required molecule, retain it and all its gross
    // inputs. Outputs outside this closure cannot enable any retained reaction.
    std::vector<std::vector<int>> producers(moleculeCount);
    for (int r = 0; r < static_cast<int>(reactions.size()); ++r)
        for (const auto& output : reactions[r].outputs)
            if (output.first >= 0 && output.first < moleculeCount && output.second > 0)
                producers[output.first].push_back(r);

    std::vector<unsigned char> relevantMolecule(moleculeCount, 0);
    std::vector<unsigned char> relevantReaction(reactions.size(), 0);
    std::vector<int> pending{target};
    relevantMolecule[target] = 1;
    for (std::size_t p = 0; p < pending.size(); ++p) {
        for (int r : producers[pending[p]]) {
            if (relevantReaction[r]) continue;
            relevantReaction[r] = 1;
            for (const auto& input : reactions[r].inputs) {
                const int m = input.first;
                if (m >= 0 && m < moleculeCount && input.second > 0 &&
                    !relevantMolecule[m]) {
                    relevantMolecule[m] = 1;
                    pending.push_back(m);
                }
            }
        }
    }
    std::vector<int> compactId(moleculeCount, -1);
    std::vector<ll> compactPrices;
    for (int m = 0; m < moleculeCount; ++m)
        if (relevantMolecule[m]) {
            compactId[m] = static_cast<int>(compactPrices.size());
            compactPrices.push_back(prices[m]);
        }
    const int compactTarget = compactId[target];
    const int activeCount = static_cast<int>(compactPrices.size());
    bool arithmeticLimited = false;
    std::vector<detail::CompactReaction> activeReactions;
    for (int r = 0; r < static_cast<int>(reactions.size()); ++r) {
        if (!relevantReaction[r]) continue;
        detail::CompactReaction reaction;
        reaction.original = r;
        bool valid = true;
        auto normalize = [&](const std::vector<std::pair<int, ll>>& source,
                             std::vector<std::pair<int, ll>>& destination,
                             bool isInput) {
            std::vector<__int128> totals(activeCount, 0);
            for (const auto& item : source) {
                if (item.first < 0 || item.first >= moleculeCount || item.second < 0) {
                    valid = false;
                    return;
                }
                int m = compactId[item.first];
                if (m < 0) {
                    if (isInput && item.second > 0) valid = false;
                    continue;
                }
                totals[m] += item.second;
                if (totals[m] > limit) {
                    valid = false;
                    arithmeticLimited = true;
                    return;
                }
            }
            for (int m = 0; m < activeCount; ++m)
                if (totals[m] > 0)
                    destination.emplace_back(m, static_cast<ll>(totals[m]));
        };
        normalize(reactions[r].inputs, reaction.inputs, true);
        normalize(reactions[r].outputs, reaction.outputs, false);
        if (valid) activeReactions.push_back(std::move(reaction));
    }

    std::vector<detail::Node> nodes;
    nodes.reserve(std::min<std::size_t>(maxStates, 100000));
    std::unordered_map<std::vector<ll>, std::size_t, detail::InventoryHash> seen;
    seen.reserve(std::min<std::size_t>(maxStates, 100000));
    std::vector<ll> empty(activeCount, 0);
    nodes.push_back(detail::Node{empty, 0, noParent, -1});
    seen.emplace(empty, 0);
    using QueueEntry = std::pair<ll, std::size_t>;
    std::priority_queue<QueueEntry, std::vector<QueueEntry>, std::greater<QueueEntry>> queue;
    queue.emplace(0, 0);
    auto reconstruct = [&](std::size_t id, int lastReaction) {
        std::vector<int> result;
        if (lastReaction >= 0) result.push_back(lastReaction);
        while (id != noParent && nodes[id].reaction >= 0) {
            result.push_back(nodes[id].reaction);
            id = nodes[id].parent;
        }
        std::reverse(result.begin(), result.end());
        return result;
    };

    while (!queue.empty()) {
        if (std::chrono::steady_clock::now() >= deadline) return answer;
        const QueueEntry entry = queue.top();
        queue.pop();
        const ll currentCost = entry.first;
        const std::size_t currentId = entry.second;
        if (currentCost != nodes[currentId].cost) continue;
        if (currentCost >= answer.cost) continue;

        // Copy before insertion: nodes may reallocate while children are added.
        const std::vector<ll> inventory = nodes[currentId].inventory;
        for (const auto& reaction : activeReactions) {
            if (std::chrono::steady_clock::now() >= deadline) return answer;
            __int128 nextCost = currentCost;
            bool legal = true;
            for (const auto& input : reaction.inputs) {
                const int m = input.first;
                const ll deficit = input.second > inventory[m] ? input.second - inventory[m] : 0;
                if (deficit == 0) continue;
                if (compactPrices[m] < 0) { legal = false; break; }
                nextCost += static_cast<__int128>(deficit) * compactPrices[m];
                if (nextCost >= answer.cost) { legal = false; break; }
            }
            if (!legal) continue;
            if (nextCost > limit) { arithmeticLimited = true; continue; }
            std::vector<ll> nextInventory = inventory;
            for (const auto& input : reaction.inputs)
                nextInventory[input.first] = std::max<ll>(0, nextInventory[input.first] - input.second);
            for (const auto& output : reaction.outputs) {
                const __int128 amount = static_cast<__int128>(nextInventory[output.first]) + output.second;
                if (amount > limit) { legal = false; arithmeticLimited = true; break; }
                nextInventory[output.first] = static_cast<ll>(amount);
            }
            if (!legal) continue;
            const ll candidateCost = static_cast<ll>(nextCost);
            if (nextInventory[compactTarget] >= needed) {
                answer.cost = candidateCost;
                answer.sequence = reconstruct(currentId, reaction.original);
                if (candidateCost == 0) { answer.exhausted = true; return answer; }
                continue;
            }
            auto old = seen.find(nextInventory);
            if (old != seen.end()) {
                const std::size_t id = old->second;
                if (candidateCost >= nodes[id].cost) continue;
                nodes[id].cost = candidateCost;
                nodes[id].parent = currentId;
                nodes[id].reaction = reaction.original;
                queue.emplace(candidateCost, id);
                continue;
            }
            if (nodes.size() >= maxStates) return answer;
            const std::size_t id = nodes.size();
            seen.emplace(nextInventory, id);
            nodes.push_back(detail::Node{std::move(nextInventory), candidateCost,
                                         currentId, reaction.original});
            queue.emplace(candidateCost, id);
        }
    }
    answer.exhausted = !arithmeticLimited;
    return answer;
}

}  // namespace inventory_search

#endif  // MOLECULE_FACTORY_INVENTORY_SEARCH_HPP

using namespace std;
typedef long long ll;

static const int CH_BUY = -1, CH_NONE = -2;
static const double DINF = 1e30;
static const ll LINF = (ll)4e18, CAP0 = 20000000;

struct Reaction {
    string text;
    vector<pair<int, ll>> in, out;     // gross amounts of one run
    vector<pair<int, ll>> cons, prod;  // net consumed / net produced (> 0 only)
};
struct Plan {
    ll cost = LINF;
    vector<int> choice;  // per molecule: CH_BUY, CH_NONE or a reaction index
    vector<int> seq;     // executed reactions, in order
    vector<int> dem;     // molecules this plan actually needed
    bool proven = false;
};

// ---------------- global data ----------------
static vector<string> names;
static unordered_map<string, int> idOf;
static vector<ll> price;                // -1 = cannot be bought
static vector<Reaction> R;
static vector<vector<int>> producers;   // producers[m] = reactions with net output of m
static int M = 0;
static unsigned searchSeed = 987654321u;
static chrono::steady_clock::time_point T0;
static double nowSec() {
    return chrono::duration<double>(chrono::steady_clock::now() - T0).count();
}

static int getId(const string& name) {
    auto it = idOf.find(name);
    if (it != idOf.end()) return it->second;
    int id = (int)names.size();
    idOf[name] = id;
    names.push_back(name);
    price.push_back(-1);
    return id;
}
static void addTo(vector<pair<int, ll>>& v, int id, ll q) {
    for (auto& e : v)
        if (e.first == id) { e.second += q; return; }
    v.push_back(make_pair(id, q));
}
static ll qtyOf(const vector<pair<int, ll>>& v, int id) {
    for (const auto& e : v)
        if (e.first == id) return e.second;
    return 0;
}
static bool nextLine(istream& in, string& line) {
    string s;
    while (getline(in, s)) {
        size_t a = 0, b = s.size();
        while (a < b && isspace((unsigned char)s[a])) ++a;
        while (b > a && isspace((unsigned char)s[b - 1])) --b;
        if (a < b) { line = s.substr(a, b - a); return true; }
    }
    return false;
}
static string deriveOutName(const string& inName) {
    size_t sl = inName.find_last_of("/\\");
    string dir = (sl == string::npos) ? "" : inName.substr(0, sl + 1);
    string base = (sl == string::npos) ? inName : inName.substr(sl + 1);
    if (base.compare(0, 3, "in_") == 0) return dir + "out_" + base.substr(3);
    return inName + ".out";
}

// ---------------- static analysis ----------------
static vector<double> v0;                  // static unit value of each molecule
static vector<int> baseChoice, molOrder, visStamp;
static int stampCtr = 0;

static bool reaches(const vector<int>& ch, int a, int m) {
    if (a == m) return true;
    if (visStamp[a] == stampCtr) return false;
    visStamp[a] = stampCtr;
    int r = ch[a];
    if (r < 0) return false;
    for (const auto& e : R[r].in)
        if (reaches(ch, e.first, m)) return true;
    return false;
}
// would using reaction r as the source of m create a dependency cycle?
static bool wouldCycle(const vector<int>& ch, int m, int r) {
    ++stampCtr;
    for (const auto& e : R[r].in)
        if (reaches(ch, e.first, m)) return true;
    return false;
}
static bool obtainable(int r) {
    for (const auto& e : R[r].in)
        if (v0[e.first] >= DINF / 2) return false;
    return true;
}
static double wgt(int x) { return v0[x] >= DINF / 2 ? 0.0 : v0[x]; }

static void computeStatic() {
    v0.assign(M, DINF);
    baseChoice.assign(M, CH_NONE);
    for (int m = 0; m < M; ++m)
        if (price[m] >= 0) { v0[m] = (double)price[m]; baseChoice[m] = CH_BUY; }
    int nR = (int)R.size();
    for (int iter = 0; iter < 5000; ++iter) {
        bool changed = false;
        for (int r = 0; r < nR; ++r) {
            if (!obtainable(r)) continue;
            double cons = 0;
            for (const auto& e : R[r].cons) cons += (double)e.second * v0[e.first];
            for (const auto& e : R[r].prod) {
                int o = e.first;
                double c = cons / (double)e.second;
                if (c + 1e-9 < v0[o]) {
                    if (wouldCycle(baseChoice, o, r)) continue;
                    v0[o] = c;
                    baseChoice[o] = r;
                    changed = true;
                }
            }
        }
        if (!changed) break;
    }
    molOrder.resize(M);
    iota(molOrder.begin(), molOrder.end(), 0);
    sort(molOrder.begin(), molOrder.end(), [&](int a, int b) {
        if (v0[a] != v0[b]) return v0[a] < v0[b];
        return a < b;
    });
}

// Estimated cost per unit of m when made with reaction r.
// alpha = share of the reaction cost split over all outputs by value (by-product credit).
static double unitCost(int r, int m, double alpha) {
    const Reaction& rc = R[r];
    double cons = 0, W = 0, qm = 1;
    for (const auto& e : rc.cons) cons += (double)e.second * wgt(e.first);
    for (const auto& e : rc.prod) {
        W += (double)e.second * wgt(e.first);
        if (e.first == m) qm = (double)e.second;
    }
    double full = cons / qm;
    if (alpha <= 0 || W <= 1e-12) return full;
    double share = cons * wgt(m) / W;
    return (1.0 - alpha) * full + alpha * share;
}

// Start plan: baseline recipes, then switch to cheaper (credited) options when acyclic.
static vector<int> makeCredited(double alpha, int target) {
    vector<int> ch = baseChoice;
    for (int m : molOrder) {
        if (v0[m] >= DINF / 2) continue;
        int cur = ch[m];
        double curC = DINF;
        if (cur == CH_BUY) { if (m != target) curC = (double)price[m]; }
        else if (cur >= 0) curC = unitCost(cur, m, alpha);
        int bestO = cur;
        double bestC = curC;
        if (m != target && price[m] >= 0 && cur != CH_BUY && (double)price[m] < bestC - 1e-9) {
            bestC = (double)price[m];
            bestO = CH_BUY;
        }
        for (int r : producers[m]) {
            if (r == cur || !obtainable(r)) continue;
            double c = unitCost(r, m, alpha);
            if (c < bestC - 1e-9 && !wouldCycle(ch, m, r)) { bestC = c; bestO = r; }
        }
        ch[m] = bestO;
    }
    return ch;
}

// ---------------- exact simulation of one plan ----------------
static vector<ll> inv;
static vector<char> inProg, dem;
static vector<int> seqBuf;
static const vector<int>* chPtr = nullptr;
static ll curCost = 0, costBound = LINF, stepCap = 0, stepCnt = 0;
static bool aborted = false;
static double trialDeadline = 1e100;

// will the molecule simply be bought by the reaction that needs it?
static inline bool buyOK(int a) {
    return price[a] >= 0 && ((*chPtr)[a] == CH_BUY || inProg[a]);
}
static bool provides(int r, int b) {
    for (const auto& e : R[r].prod)
        if (e.first == b) return true;
    return false;
}

static bool execReaction(int r) {
    if ((stepCnt & 127) == 0 && nowSec() >= trialDeadline) { aborted = true; return false; }
    const Reaction& rc = R[r];
    for (const auto& e : rc.in) {
        ll have = inv[e.first];
        if (have >= e.second) {
            inv[e.first] = have - e.second;
        } else {
            if (price[e.first] < 0) return false;
            __int128 next = (__int128)curCost + (__int128)(e.second - have) * price[e.first];
            if (next >= costBound || next >= LINF) { aborted = true; return false; }
            curCost = (ll)next;
            inv[e.first] = 0;
            if (curCost >= costBound) { aborted = true; return false; }
        }
    }
    for (const auto& e : rc.out) {
        if ((__int128)inv[e.first] + e.second >= LINF) { aborted = true; return false; }
        inv[e.first] += e.second;
    }
    seqBuf.push_back(r);
    if (++stepCnt > stepCap) { aborted = true; return false; }
    return true;
}

static bool ensureMol(int m, ll need, bool top);

// make sure all ingredients of one run of reaction r are available
static bool gather(int r) {
    const vector<pair<int, ll>>& ins = R[r].in;
    const int k = (int)ins.size();
    for (int guard = 0; guard < 1000; ++guard) {
        int cntUns = 0, firstUns = -1;
        for (int t = 0; t < k; ++t) {
            int a = ins[t].first;
            if (inv[a] >= ins[t].second) continue;
            dem[a] = 1;
            if (buyOK(a)) continue;
            if (firstUns < 0) firstUns = t;
            ++cntUns;
        }
        if (cntUns == 0) return true;
        int pick = firstUns;
        if (cntUns > 1) {
            // prefer the ingredient whose reaction also yields another missing ingredient
            int bestScore = -1;
            for (int t = 0; t < k; ++t) {
                int a = ins[t].first;
                if (inv[a] >= ins[t].second || buyOK(a)) continue;
                int score = 0;
                int ra = (*chPtr)[a];
                if (ra >= 0) {
                    for (int t2 = 0; t2 < k; ++t2) {
                        if (t2 == t) continue;
                        int b = ins[t2].first;
                        if (inv[b] < ins[t2].second && !buyOK(b) && provides(ra, b)) ++score;
                    }
                }
                if (score > bestScore) { bestScore = score; pick = t; }
            }
        }
        if (!ensureMol(ins[pick].first, ins[pick].second, false) || aborted) return false;
    }
    return false;
}

// produce molecule m until the inventory holds at least `need`
static bool ensureMol(int m, ll need, bool top) {
    if (aborted || nowSec() >= trialDeadline) { aborted = true; return false; }
    if (inv[m] >= need) return true;
    dem[m] = 1;
    int opt = (*chPtr)[m];
    if (opt == CH_BUY) return !top && price[m] >= 0;
    if (opt < 0) return false;
    if (inProg[m]) return price[m] >= 0;
    inProg[m] = 1;
    bool ok = true;
    while (inv[m] < need) {
        if (!gather(opt) || !execReaction(opt)) { ok = false; break; }
    }
    inProg[m] = 0;
    return ok;
}

// returns total purchase cost, or LINF if invalid / aborted / not cheaper than `bound`
static ll runTrial(const vector<int>& ch, int target, ll K, ll bound, ll cap) {
    chPtr = &ch;
    fill(inv.begin(), inv.end(), 0);
    fill(inProg.begin(), inProg.end(), 0);
    fill(dem.begin(), dem.end(), 0);
    seqBuf.clear();
    curCost = 0; costBound = bound; stepCap = cap; stepCnt = 0; aborted = false;
    bool ok = ensureMol(target, K, true);
    if (!ok || aborted || inv[target] < K) return LINF;
    return curCost;
}

// ---------------- search ----------------
static vector<vector<int>> opts;   // possible sources per molecule

static void buildOpts(int target) {
    opts.assign(M, vector<int>());
    for (int m = 0; m < M; ++m) {
        if (v0[m] >= DINF / 2) continue;
        if (m != target && price[m] >= 0) opts[m].push_back(CH_BUY);
        for (int r : producers[m])
            if (obtainable(r)) opts[m].push_back(r);
    }
}
static void collectDem(vector<int>& out) {
    out.clear();
    for (int m = 0; m < M; ++m)
        if (dem[m]) out.push_back(m);
}

// hill climbing: switch one molecule's source, keep only real improvements
static void descend(Plan& p, int target, ll K, double deadline, mt19937& rng) {
    ll cap = max<ll>(200000, 10LL * (ll)p.seq.size());
    bool improved = true;
    while (improved && nowSec() < deadline) {
        improved = false;
        vector<int> order = p.dem;
        shuffle(order.begin(), order.end(), rng);
        for (int m : order) {
            if (nowSec() >= deadline) break;
            const vector<int>& op = opts[m];
            if (op.size() < 2) continue;
            for (size_t i = 0; i < op.size(); ++i) {
                int o = op[i];
                if (o == p.choice[m]) continue;
                if (o >= 0 && wouldCycle(p.choice, m, o)) continue;
                int old = p.choice[m];
                p.choice[m] = o;
                ll c = runTrial(p.choice, target, K, p.cost, cap);
                if (c < p.cost) {
                    p.cost = c;
                    p.seq = seqBuf;
                    collectDem(p.dem);
                    improved = true;
                } else {
                    p.choice[m] = old;
                }
            }
        }
    }
}

static Plan solveCase(ll K, int target, double deadline) {
    Plan best;
    buildOpts(target);

    vector<Plan> starts;
    vector<vector<int>> tried;
    const double alphas[3] = {0.0, 0.5, 1.0};
    for (int ai = 0; ai < 3; ++ai) {
        vector<int> ch = makeCredited(alphas[ai], target);
        bool dup = false;
        for (const auto& t : tried)
            if (t == ch) { dup = true; break; }
        if (dup) continue;
        tried.push_back(ch);
        Plan p;
        p.choice = ch;
        p.cost = runTrial(p.choice, target, K, LINF, CAP0);
        if (p.cost >= LINF) continue;
        p.seq = seqBuf;
        collectDem(p.dem);
        starts.push_back(std::move(p));
    }
    if (starts.empty()) return best;

    sort(starts.begin(), starts.end(), [](const Plan& a, const Plan& b) { return a.cost < b.cost; });
    mt19937 rng(searchSeed);
    for (Plan& s : starts) {
        if (nowSec() >= deadline) break;
        descend(s, target, K, deadline, rng);
    }
    best = starts[0];
    for (const Plan& s : starts)
        if (s.cost < best.cost) best = s;

    // iterated local search: random kicks + descent
    int fails = 0;
    while (nowSec() < deadline) {
        vector<int> pool;
        for (int m : best.dem)
            if (opts[m].size() >= 2) pool.push_back(m);
        if (pool.empty()) break;
        if (fails > 30 + 2 * (int)pool.size()) break;
        Plan cand;
        cand.choice = best.choice;
        int nChg = 1 + (int)(rng() % 3);
        bool changed = false;
        for (int t = 0; t < nChg; ++t) {
            int m = pool[rng() % pool.size()];
            int o = opts[m][rng() % opts[m].size()];
            if (o == cand.choice[m]) continue;
            if (o >= 0 && wouldCycle(cand.choice, m, o)) continue;
            cand.choice[m] = o;
            changed = true;
        }
        if (!changed) { ++fails; continue; }
        ll cap = max<ll>(200000, 10LL * (ll)best.seq.size());
        cand.cost = runTrial(cand.choice, target, K, LINF, cap);
        if (cand.cost >= LINF) { ++fails; continue; }
        cand.seq = seqBuf;
        collectDem(cand.dem);
        descend(cand, target, K, deadline, rng);
        if (cand.cost < best.cost) { best = cand; fails = 0; }
        else ++fails;
    }
    return best;
}

// independent re-check of a finished sequence
static bool verifySeq(const vector<int>& seq, int target, ll K, ll& cost) {
    vector<ll> iv(M, 0);
    cost = 0;
    for (int r : seq) {
        for (const auto& e : R[r].in) {
            ll have = iv[e.first];
            if (have >= e.second) {
                iv[e.first] = have - e.second;
            } else {
                if (price[e.first] < 0) return false;
                __int128 next = (__int128)cost + (__int128)(e.second - have) * price[e.first];
                if (next >= LINF) return false;
                cost = (ll)next;
                iv[e.first] = 0;
            }
        }
        for (const auto& e : R[r].out) {
            if ((__int128)iv[e.first] + e.second >= LINF) return false;
            iv[e.first] += e.second;
        }
    }
    return iv[target] >= K;
}

// ---- optimization.hpp ----
#ifndef FACTORY_OPTIMIZATION_HPP
#define FACTORY_OPTIMIZATION_HPP

// Included after the original simulator. Every heuristic candidate is replayed
// from an empty inventory before it can replace the incumbent.
// ---- adaptive_planner.hpp ----

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

// ---- sequence_polish.hpp ----
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
    chPtr = &ch;
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

static Plan solveBetter(ll K, int target, double deadline, const Plan& previous) {
    const double start = nowSec();
    const double available = max(0.001, deadline - start);
    trialDeadline = start + available * 0.30;
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
    trialDeadline = start + available * 0.40;
    for (const auto& policy : policyBank) {
        if (nowSec() >= trialDeadline) break;
        ll cost = runTrial(policy, target, K, best.cost, CAP0);
        if (cost < best.cost) considerSequence(best, seqBuf, target, K, policy);
    }
    if (best.cost < LINF && nowSec() < trialDeadline) {
        mt19937 rng(392473u + target);
        descend(best, target, K, trialDeadline, rng);
    }
    if (nowSec() < start + available * 0.52)
        best = adaptiveSolve(K, target, best, start + available * 0.52);
    if (best.cost < LINF && nowSec() < start + available * 0.75)
        linearImprove(best, target, K, start + available * 0.75);
    if (nowSec() < start + available * 0.88)
        best = adaptiveSolve(K, target, best, start + available * 0.88);
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


int main(int argc, char** argv) {
    T0 = chrono::steady_clock::now();
    if (argc < 2) { cerr << "Usage: factory INPUT [OUTPUT] [SECONDS=60] [SEED=987654321]\n"; return 1; }
    string inName = argv[1];
    string outName = (argc > 2) ? string(argv[2]) : deriveOutName(inName);
    double totalTime = (argc > 3) ? atof(argv[3]) : 60.0;
    if (argc > 4) searchSeed = (unsigned)strtoul(argv[4], nullptr, 10);

    if (!isfinite(totalTime) || totalTime <= 0) { cerr << "SECONDS must be positive\n"; return 1; }
    if (inName == outName) { cerr << "Input and output paths must differ\n"; return 1; }
    ifstream fin(inName);
    if (!fin) {
        cerr << "Cannot open input file: " << inName << "\n";
        return 1;
    }

    string line;
    if (!nextLine(fin, line)) return 0;
    int nBuy = stoi(line);
    for (int i = 0; i < nBuy; ++i) {
        if (!nextLine(fin, line)) break;
        stringstream ss(line);
        ll p = 0;
        string arrow, nm;
        ss >> p >> arrow >> nm;
        int id = getId(nm);
        price[id] = p;
    }

    nextLine(fin, line);
    int nReact = stoi(line);
    R.assign(nReact, Reaction());
    for (int r = 0; r < nReact; ++r) {
        if (!nextLine(fin, line)) break;
        R[r].text = line;
        stringstream ss(line);
        string tok;
        bool after = false;
        while (ss >> tok) {
            if (tok == "->") { after = true; continue; }
            ll q = stoll(tok);
            string nm;
            ss >> nm;
            while (!nm.empty() && nm.back() == ',') nm.pop_back();
            int id = getId(nm);
            addTo(after ? R[r].out : R[r].in, id, q);
        }
    }

    nextLine(fin, line);
    int nCases = stoi(line);
    vector<pair<ll, int>> cases;
    for (int i = 0; i < nCases; ++i) {
        if (!nextLine(fin, line)) break;
        stringstream ss(line);
        ll cnt = 0;
        string nm;
        ss >> cnt >> nm;
        int id = getId(nm);
        cases.push_back(make_pair(cnt, id));
    }
    nCases = (int)cases.size();
    fin.close();

    M = (int)names.size();
    producers.assign(M, vector<int>());
    for (int r = 0; r < nReact; ++r) {
        Reaction& rc = R[r];
        for (const auto& e : rc.in) {
            ll q = e.second - qtyOf(rc.out, e.first);
            if (q > 0) rc.cons.push_back(make_pair(e.first, q));
        }
        for (const auto& e : rc.out) {
            ll q = e.second - qtyOf(rc.in, e.first);
            if (q > 0) {
                rc.prod.push_back(make_pair(e.first, q));
                producers[e.first].push_back(r);
            }
        }
    }

    inv.assign(M, 0);
    inProg.assign(M, 0);
    dem.assign(M, 0);
    visStamp.assign(M, 0);
    seqBuf.reserve(1 << 16);
    computeStatic();

    vector<Plan> plans(nCases);
    // Resume a matching output only after independently validating every case.
    vector<Plan> previous(nCases);
    {
        ifstream old(outName);
        if (old) {
            unordered_map<string, int> lookup;
            for (int r = 0; r < nReact; ++r) lookup[R[r].text] = r;
            vector<vector<int>> blocks(1);
            string s;
            bool valid = true;
            while (getline(old, s)) {
                while (!s.empty() && isspace((unsigned char)s.back())) s.pop_back();
                if (s.empty()) { if (!blocks.back().empty()) blocks.emplace_back(); continue; }
                auto it = lookup.find(s);
                if (it == lookup.end()) { valid = false; break; }
                blocks.back().push_back(it->second);
            }
            if (!blocks.empty() && blocks.back().empty()) blocks.pop_back();
            if (valid && blocks.size() == cases.size()) {
                for (int i = 0; i < nCases; ++i) {
                    ll cost;
                    if (verifySeq(blocks[i], cases[i].second, cases[i].first, cost)) {
                        previous[i].cost = cost;
                        previous[i].seq = move(blocks[i]);
                        previous[i].choice = baseChoice;
                    }
                }
            }
        }
    }
    ll total = 0;
    for (int i = 0; i < nCases; ++i) {
        double remaining = totalTime - nowSec();
        double share = max(0.2, remaining / (double)(nCases - i));
        plans[i] = solveBetter(cases[i].first, cases[i].second, nowSec() + share, previous[i]);
        if (plans[i].cost >= LINF) {
            cerr << "Case " << (i + 1) << ": NO VALID PLAN FOUND (" << cases[i].first << " "
                 << names[cases[i].second] << ")\n";
            return 2;
        }
        ll vc = 0;
        bool ok = verifySeq(plans[i].seq, cases[i].second, cases[i].first, vc);
        if (!ok || vc != plans[i].cost) {
            cerr << "ERROR: verification mismatch in case " << (i + 1) << "\n"; return 3;
        }
        total += plans[i].cost;
        cout << "Case " << (i + 1) << "/" << nCases << ": " << cases[i].first << " "
             << names[cases[i].second] << "  cost = " << plans[i].cost
             << "  reactions = " << plans[i].seq.size()
             << (plans[i].proven ? "  [optimal]" : "  [best found]") << "\n";
    }
    cout << "Total cost = " << total << "   (" << nowSec() << " s)\n";

    string outBuf;
    for (int i = 0; i < nCases; ++i) {
        for (int r : plans[i].seq) {
            outBuf += R[r].text;
            outBuf += '\n';
        }
        if (i + 1 < nCases) outBuf += '\n';
    }
    ofstream fout(outName);
    if (!fout) { cerr << "Cannot open output: " << outName << "\n"; return 4; }
    fout << outBuf;
    fout.close();
    if (!fout) { cerr << "Cannot write output\n"; return 4; }
    cout << "Written to " << outName << "\n";
    return 0;
}
