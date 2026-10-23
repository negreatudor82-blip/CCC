// Molecule factory - cost optimizer (Level 3)
// Build: g++ -O2 -std=c++17 -o factory factory.cpp
// Run:   ./factory <input> [output] [max_search_seconds]
//   e.g. ./factory in_level-3_0-example.txt out_level-3_0-example.txt 30
//
// Idea:
//  * A plan gives every molecule one source: buy it, or make it with one reaction.
//  * A plan is executed just-in-time on an exact inventory simulation (inventory is
//    used first, only the shortfall is bought), so by-products are reused.
//  * Start plans come from static unit costs with 0 / 50 / 100 % by-product credit.
//  * Local search + random kicks switch one molecule's source at a time and keep a
//    change only if the simulated total cost drops.
#include <bits/stdc++.h>
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
};

// ---------------- global data ----------------
static vector<string> names;
static unordered_map<string, int> idOf;
static vector<ll> price;                // -1 = cannot be bought
static vector<Reaction> R;
static vector<vector<int>> producers;   // producers[m] = reactions with net output of m
static int M = 0;
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
    const Reaction& rc = R[r];
    for (const auto& e : rc.in) {
        ll have = inv[e.first];
        if (have >= e.second) {
            inv[e.first] = have - e.second;
        } else {
            if (price[e.first] < 0) return false;
            curCost += (e.second - have) * price[e.first];
            inv[e.first] = 0;
            if (curCost >= costBound) { aborted = true; return false; }
        }
    }
    for (const auto& e : rc.out) inv[e.first] += e.second;
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
    if (aborted) return false;
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
    mt19937 rng(987654321u);
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
                cost += (e.second - have) * price[e.first];
                iv[e.first] = 0;
            }
        }
        for (const auto& e : R[r].out) iv[e.first] += e.second;
    }
    return iv[target] >= K;
}

int main(int argc, char** argv) {
    T0 = chrono::steady_clock::now();
    string inName = (argc > 1) ? string(argv[1]) : string("in_level-3_4-enormous.txt");
    string outName = (argc > 2) ? string(argv[2]) : deriveOutName(inName);
    double totalTime = (argc > 3) ? atof(argv[3]) : 30.0;

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
    ll total = 0;
    for (int i = 0; i < nCases; ++i) {
        double remaining = totalTime - nowSec();
        double share = max(0.2, remaining / (double)(nCases - i));
        plans[i] = solveCase(cases[i].first, cases[i].second, nowSec() + share);
        if (plans[i].cost >= LINF) {
            cerr << "Case " << (i + 1) << ": NO VALID PLAN FOUND (" << cases[i].first << " "
                 << names[cases[i].second] << ")\n";
            continue;
        }
        ll vc = 0;
        bool ok = verifySeq(plans[i].seq, cases[i].second, cases[i].first, vc);
        if (!ok || vc != plans[i].cost)
            cerr << "WARNING: verification mismatch in case " << (i + 1) << "\n";
        total += plans[i].cost;
        cout << "Case " << (i + 1) << "/" << nCases << ": " << cases[i].first << " "
             << names[cases[i].second] << "  cost = " << plans[i].cost
             << "  reactions = " << plans[i].seq.size() << "\n";
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
    fout << outBuf;
    fout.close();
    cout << "Written to " << outName << "\n";
    return 0;
}