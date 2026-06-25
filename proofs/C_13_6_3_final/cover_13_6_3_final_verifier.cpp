#include <algorithm>
#include <array>
#include <chrono>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <string>
#include <vector>
using namespace std;

// ---------- small triple bitsets on C(12,3)=220 ----------
using u64 = uint64_t;
struct TMask { u64 w[4]; };
static inline void clearT(TMask& x){ x.w[0]=x.w[1]=x.w[2]=x.w[3]=0; }
static inline void setT(TMask& x,int i){ x.w[i>>6] |= 1ULL << (i&63); }
static inline bool testT(const TMask& x,int i){ return (x.w[i>>6] >> (i&63)) & 1ULL; }

// ---------- local incidence generation ----------
// A degree-d point gives d columns of size 5 on the other 12 points.  Equivalently
// the other 12 points are rows: pairwise-intersecting subsets of a d-set, each
// of size at least 3, with every column sum 5.
struct LocalGenerator {
    int D;
    vector<int> chosen;
    int cap[10]{};
    int targetN[11]{}, usedN[11]{};
    vector<int> candBySize[11];
    vector<array<int,12>> out;

    explicit LocalGenerator(int d): D(d) {
        for (int m=0; m<(1<<D); ++m) {
            int s = __builtin_popcount((unsigned)m);
            if (s >= 3) candBySize[s].push_back(m);
        }
        for (int s=0; s<=D; ++s) sort(candBySize[s].begin(), candBySize[s].end());
    }

    bool residualColumnsSorted() const {
        int k = (int)chosen.size();
        int sig[10]{};
        for (int r=0; r<k; ++r) {
            int bit = 1 << (k-1-r);
            int m = chosen[r];
            for (int c=0; c<D; ++c) if ((m>>c)&1) sig[c] |= bit;
        }
        // After the first row is fixed to columns {0,1,2}, columns 0..2 and
        // 3..D-1 are separately permutable.  Keep the current prefix sorted.
        for (int c=0; c+1<3; ++c) if (sig[c] < sig[c+1]) return false;
        for (int c=3; c+1<D; ++c) if (sig[c] < sig[c+1]) return false;
        return true;
    }

    bool columnsDistinctAtLeaf() const {
        int col[10]{};
        for (int r=0; r<12; ++r) {
            int m = chosen[r];
            for (int c=0; c<D; ++c) if ((m>>c)&1) col[c] |= (1 << r);
        }
        for (int i=0; i<D; ++i) for (int j=i+1; j<D; ++j)
            if (col[i] == col[j]) return false;
        return true;
    }

    void leaf() {
        if (!columnsDistinctAtLeaf()) return;
        array<int,12> a{};
        for (int i=0; i<12; ++i) a[i] = chosen[i];
        out.push_back(a);
    }

    void rec(int sizePhase, int startIdx) {
        int rem = 12 - (int)chosen.size();
        int sumCaps = 0;
        for (int c=0; c<D; ++c) {
            if (cap[c] < 0 || cap[c] > rem) return;
            sumCaps += cap[c];
        }
        if (rem == 0) {
            if (sumCaps == 0) leaf();
            return;
        }
        int expected = 0, remRows = 0;
        for (int s=3; s<=D; ++s) {
            int left = targetN[s] - usedN[s];
            if (left < 0) return;
            expected += s * left;
            remRows += left;
        }
        if (remRows != rem || expected != sumCaps) return;
        while (sizePhase <= D && usedN[sizePhase] >= targetN[sizePhase]) {
            ++sizePhase;
            startIdx = 0;
        }
        if (sizePhase > D) return;

        int forced = 0;
        for (int c=0; c<D; ++c) if (cap[c] == rem) forced |= 1 << c;

        const auto& vec = candBySize[sizePhase];
        for (int idx=startIdx; idx<(int)vec.size(); ++idx) {
            int m = vec[idx];
            if ((m & forced) != forced) continue;
            bool ok = true;
            for (int old : chosen) if ((m & old) == 0) { ok = false; break; }
            if (!ok) continue;
            for (int c=0; c<D; ++c) {
                if ((m>>c)&1) {
                    if (cap[c] <= 0) { ok = false; break; }
                } else {
                    if (cap[c] > rem-1) { ok = false; break; }
                }
            }
            if (!ok) continue;
            for (int c=0; c<D; ++c) if ((m>>c)&1) --cap[c];
            chosen.push_back(m);
            ++usedN[sizePhase];
            if (residualColumnsSorted()) rec(sizePhase, idx);
            --usedN[sizePhase];
            chosen.pop_back();
            for (int c=0; c<D; ++c) if ((m>>c)&1) ++cap[c];
        }
    }

    void runDistribution(const vector<int>& target) {
        memset(targetN,0,sizeof targetN);
        memset(usedN,0,sizeof usedN);
        for (int s=3; s<=D; ++s) targetN[s] = target[s];
        if (targetN[3] == 0) return; // total weight forces at least one size-3 row
        for (int c=0; c<D; ++c) cap[c] = 5;
        chosen.clear();
        int first = 0b111;
        chosen.push_back(first);
        ++usedN[3];
        for (int c=0; c<3; ++c) --cap[c];
        int firstIdx = (int)(lower_bound(candBySize[3].begin(), candBySize[3].end(), first) - candBySize[3].begin());
        rec(3, firstIdx);
    }

    void genTargets(int s, int remainingRows, int remainingSum, vector<int>& target) {
        if (s == D+1) {
            if (remainingRows == 0 && remainingSum == 0) runDistribution(target);
            return;
        }
        for (int n=0; n<=remainingRows; ++n) {
            if (n*s > remainingSum) break;
            target[s] = n;
            genTargets(s+1, remainingRows-n, remainingSum-n*s, target);
            target[s] = 0;
        }
    }

    vector<array<int,12>> generate() {
        vector<int> target(D+1,0);
        genTargets(3, 12, 5*D, target);
        return out;
    }
};

static bool isPermutation(const vector<int>& p, int n) {
    vector<int> seen(n,0);
    if ((int)p.size()!=n) return false;
    for (int x:p) { if (x<0 || x>=n || seen[x]) return false; seen[x]=1; }
    return true;
}

static array<int,12> transformRows(const array<int,12>& rows, const vector<int>& pmap, const vector<int>& cmap) {
    array<int,12> out{};
    for (int i=0; i<12; ++i) {
        int ni = pmap[i];
        int nr = 0;
        for (int c=0; c<9; ++c) if ((rows[i]>>c)&1) nr |= 1 << cmap[c];
        out[ni] = nr;
    }
    return out;
}

// ---------- CDCL SAT solver for the relaxed outside-extension CNF ----------
struct CDCL {
    int nvars=0, branchVars=0;
    vector<vector<int>> cls, watches;
    vector<int> val, lev, reason, trail, trail_lim;
    vector<double> act;
    int qhead=0;
    double var_inc=1.0;

    int newVar(){
        ++nvars; val.push_back(0); lev.push_back(0); reason.push_back(-1); act.push_back(0.0);
        watches.resize(2*nvars+2);
        return nvars;
    }
    int litIndex(int lit) const { int v=abs(lit)-1; return 2*v + (lit<0); }
    int negIndex(int lit) const { return litIndex(-lit); }
    bool litTrue(int lit) const { int x=val[abs(lit)-1]; return lit>0 ? x==1 : x==-1; }
    bool litFalse(int lit) const { int x=val[abs(lit)-1]; return lit>0 ? x==-1 : x==1; }
    int level() const { return (int)trail_lim.size(); }

    bool enqueue(int lit, int rsn) {
        int v=abs(lit), want = lit>0 ? 1 : -1;
        if (val[v-1] != 0) return val[v-1] == want;
        val[v-1] = want; lev[v-1] = level(); reason[v-1] = rsn; trail.push_back(lit);
        return true;
    }

    int addClause(vector<int> c) {
        sort(c.begin(), c.end(), [](int a,int b){ return abs(a)==abs(b) ? a<b : abs(a)<abs(b); });
        vector<int> d; int last = 0;
        for (int l:c) {
            if (l == last) continue;
            if (last == -l) return -2; // tautology
            d.push_back(l); last = l;
        }
        int id = (int)cls.size(); cls.push_back(d);
        if (d.size() >= 2) {
            watches[litIndex(d[0])].push_back(id);
            watches[litIndex(d[1])].push_back(id);
        }
        for (int l:d) if (abs(l) <= branchVars) act[abs(l)-1] += 1.0;
        return id;
    }

    int propagate() {
        while (qhead < (int)trail.size()) {
            int lit = trail[qhead++];
            int idx = negIndex(lit);
            auto& ws = watches[idx];
            int i = 0;
            for (int k=0; k<(int)ws.size(); ++k) {
                int ci = ws[k]; auto& c = cls[ci];
                if (c.size() == 1) {
                    ws[i++] = ci;
                    if (litFalse(c[0])) { for (++k; k<(int)ws.size(); ++k) ws[i++] = ws[k]; ws.resize(i); return ci; }
                    continue;
                }
                if (!litFalse(c[0])) swap(c[0], c[1]);
                if (!litFalse(c[0])) { ws[i++] = ci; continue; }
                if (litTrue(c[1])) { ws[i++] = ci; continue; }
                bool found = false;
                for (int j=2; j<(int)c.size(); ++j) if (!litFalse(c[j])) {
                    swap(c[0], c[j]); watches[litIndex(c[0])].push_back(ci); found = true; break;
                }
                if (found) continue;
                ws[i++] = ci;
                if (litFalse(c[1])) { for (++k; k<(int)ws.size(); ++k) ws[i++] = ws[k]; ws.resize(i); return ci; }
                if (!enqueue(c[1], ci)) { for (++k; k<(int)ws.size(); ++k) ws[i++] = ws[k]; ws.resize(i); return ci; }
            }
            ws.resize(i);
        }
        return -1;
    }

    void newDecisionLevel(){ trail_lim.push_back((int)trail.size()); }
    void cancelUntil(int lvl) {
        if (level() <= lvl) return;
        int lim = (lvl == 0 ? 0 : trail_lim[lvl]);
        for (int i=(int)trail.size()-1; i>=lim; --i) {
            int v = abs(trail[i])-1;
            val[v]=0; reason[v]=-1; lev[v]=0;
        }
        trail.resize(lim); trail_lim.resize(lvl); qhead = min(qhead, lim);
    }
    int pickBranch() {
        int best = 0; double ba = -1.0;
        for (int v=1; v<=branchVars; ++v) if (val[v-1]==0 && act[v-1] > ba) { ba = act[v-1]; best = v; }
        return best;
    }
    bool analyze(int confl, vector<int>& learnt, int& bt) {
        vector<char> seen(nvars+1,0); learnt.clear(); learnt.push_back(0);
        int pathC=0, p=0, idx=(int)trail.size()-1;
        vector<int>* cp = &cls[confl];
        do {
            for (int q:*cp) {
                int v=abs(q); if (v==abs(p)) continue;
                if (!seen[v] && lev[v-1] > 0) {
                    seen[v]=1;
                    if (v <= branchVars) act[v-1] += var_inc;
                    if (lev[v-1] == level()) ++pathC; else learnt.push_back(q);
                }
            }
            while (idx>=0 && (!seen[abs(trail[idx])] || lev[abs(trail[idx])-1] != level())) --idx;
            if (idx < 0) return false;
            p = trail[idx--]; seen[abs(p)] = 0; --pathC;
            int r = reason[abs(p)-1];
            if (pathC > 0) { if (r < 0) return false; cp = &cls[r]; }
        } while (pathC > 0);
        learnt[0] = -p;
        bt = 0;
        for (int i=1; i<(int)learnt.size(); ++i) bt = max(bt, lev[abs(learnt[i])-1]);
        var_inc *= 1.01;
        return true;
    }
    bool solve() {
        for (int i=0; i<(int)cls.size(); ++i) {
            if (cls[i].empty()) return false;
            if (cls[i].size()==1 && !enqueue(cls[i][0], i)) return false;
        }
        int confl = propagate(); if (confl != -1) return false;
        int conflicts = 0;
        while (true) {
            confl = propagate();
            if (confl != -1) {
                ++conflicts;
                if (level() == 0) return false;
                vector<int> learnt; int bt=0;
                if (!analyze(confl, learnt, bt)) return false;
                cancelUntil(bt);
                int ci = addClause(learnt);
                if (ci >= 0 && !enqueue(cls[ci][0], ci)) return false;
                if (conflicts % 100000 == 0) { cancelUntil(0); qhead = (int)trail.size(); }
            } else {
                int v = pickBranch();
                if (v == 0) return true;
                newDecisionLevel();
                // Branch true first: selected blocks are usually needed by cover clauses.
                if (!enqueue(v, -1)) return false;
            }
        }
    }
};

static void addAtMost(CDCL& S, const vector<int>& lits, int K) {
    int n = (int)lits.size();
    if (K >= n) return;
    if (K < 0) { S.addClause({}); return; }
    if (K == 0) { for (int l:lits) S.addClause({-l}); return; }
    vector<vector<int>> seq(n, vector<int>(K+1,0));
    for (int i=0; i<n; ++i) for (int j=1; j<=K; ++j) seq[i][j] = S.newVar();
    for (int i=0; i<n; ++i) {
        S.addClause({-lits[i], seq[i][1]});
        if (i > 0) S.addClause({-seq[i-1][1], seq[i][1]});
        for (int j=2; j<=K; ++j) if (i > 0) {
            S.addClause({-seq[i-1][j], seq[i][j]});
            S.addClause({-lits[i], -seq[i-1][j-1], seq[i][j]});
        }
        if (i > 0) S.addClause({-lits[i], -seq[i-1][K]});
    }
}

// ---------- relaxed outside extension instance ----------
struct OutsideData {
    vector<array<int,3>> triples;
    int tid[12][12][12];
    vector<array<int,6>> blocks;
    vector<int> pointVars[12];
    vector<int> triVars[220];

    OutsideData() {
        memset(tid,-1,sizeof tid);
        for (int a=0; a<12; ++a) for (int b=a+1; b<12; ++b) for (int c=b+1; c<12; ++c) {
            tid[a][b][c] = (int)triples.size(); triples.push_back({a,b,c});
        }
        for (int a=0; a<12; ++a) for (int b=a+1; b<12; ++b) for (int c=b+1; c<12; ++c)
        for (int d=c+1; d<12; ++d) for (int e=d+1; e<12; ++e) for (int f=e+1; f<12; ++f) {
            array<int,6> B = {a,b,c,d,e,f};
            int var = (int)blocks.size() + 1;
            blocks.push_back(B);
            for (int p:B) pointVars[p].push_back(var);
            for (int i=0; i<6; ++i) for (int j=i+1; j<6; ++j) for (int k=j+1; k<6; ++k)
                triVars[tid[B[i]][B[j]][B[k]]].push_back(var);
        }
    }
};

static void localStatsFromRows(const array<int,12>& rows, const OutsideData& OD, int sdeg[12], bool localCovered[220]) {
    memset(sdeg, 0, 12*sizeof(int));
    memset(localCovered, 0, 220*sizeof(bool));
    for (int c=0; c<9; ++c) {
        vector<int> pts;
        for (int p=0; p<12; ++p) if ((rows[p]>>c)&1) { pts.push_back(p); ++sdeg[p]; }
        if (pts.size() != 5) { cerr << "bad generated local column size\n"; exit(2); }
        for (int i=0; i<5; ++i) for (int j=i+1; j<5; ++j) for (int k=j+1; k<5; ++k) {
            int a=pts[i], b=pts[j], c3=pts[k];
            int ar[3] = {a,b,c3}; sort(ar,ar+3);
            localCovered[OD.tid[ar[0]][ar[1]][ar[2]]] = true;
        }
    }
}

static bool relaxedExtensionSAT(const array<int,12>& rows, const OutsideData& OD, vector<array<int,6>>* witness=nullptr) {
    int sdeg[12]; bool covered[220];
    localStatsFromRows(rows, OD, sdeg, covered);
    CDCL S;
    for (int i=0; i<924; ++i) S.newVar();
    S.branchVars = 924;
    vector<int> all(924); iota(all.begin(), all.end(), 1);
    addAtMost(S, all, 11);
    for (int p=0; p<12; ++p) {
        int ub = 12 - sdeg[p]; // original global degree is at most 12 when all degrees are >= 9 and sum is 120
        addAtMost(S, OD.pointVars[p], ub);
    }
    int leaveCount = 0;
    for (int ti=0; ti<220; ++ti) if (!covered[ti]) {
        ++leaveCount;
        S.addClause(OD.triVars[ti]);
    }
    bool sat = S.solve();
    if (sat && witness) {
        witness->clear();
        for (int i=0; i<924; ++i) if (S.val[i] == 1) witness->push_back(OD.blocks[i]);
    }
    return sat;
}

int main(int argc, char** argv) {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    string certPath = (argc >= 2 ? argv[1] : "cover_13_6_3_d9_iso_certificate.txt");
    ifstream test(certPath);
    if (!test && argc < 2) certPath = "/mnt/data/cover_13_6_3_d9_iso_certificate.txt";

    cerr << "Generating degree-8 local cases...\n";
    LocalGenerator gen8(8);
    auto d8 = gen8.generate();
    cout << "degree-8 local cases: " << d8.size() << "\n";
    if (!d8.empty()) { cerr << "degree-8 local contradiction failed\n"; return 1; }

    cerr << "Generating degree-9 local cases...\n";
    LocalGenerator gen9(9);
    auto rows = gen9.generate();
    cout << "degree-9 local cases generated: " << rows.size() << "\n";
    if (rows.size() != 3385) { cerr << "unexpected degree-9 local case count\n"; return 1; }

    ifstream in(certPath);
    if (!in) { cerr << "could not open iso certificate: " << certPath << "\n"; return 1; }
    int nClasses; in >> nClasses;
    vector<int> reps(nClasses);
    for (int i=0; i<nClasses; ++i) in >> reps[i];
    if (nClasses <= 0) { cerr << "bad class count\n"; return 1; }
    vector<int> caseClass(rows.size(), -1);
    cout << "iso certificate classes: " << nClasses << "\n";

    for (int i=0; i<(int)rows.size(); ++i) {
        int cls; in >> cls;
        if (!in || cls < 0 || cls >= nClasses) { cerr << "bad class entry at case " << i << "\n"; return 1; }
        vector<int> pmap(12), cmap(9);
        for (int j=0; j<12; ++j) in >> pmap[j];
        string semi; in >> semi;
        for (int j=0; j<9; ++j) in >> cmap[j];
        if (!isPermutation(pmap,12) || !isPermutation(cmap,9) || semi != ";") {
            cerr << "bad permutation certificate at case " << i << "\n"; return 1;
        }
        if (reps[cls] < 0 || reps[cls] >= (int)rows.size()) { cerr << "bad rep index\n"; return 1; }
        auto tr = transformRows(rows[i], pmap, cmap);
        if (tr != rows[reps[cls]]) {
            cerr << "isomorphism check failed at case " << i << " rep " << reps[cls] << "\n";
            return 1;
        }
        caseClass[i] = cls;
    }

    vector<char> seen(nClasses,0);
    for (int c:caseClass) seen[c]=1;
    for (int c=0; c<nClasses; ++c) if (!seen[c]) { cerr << "unused class " << c << "\n"; return 1; }
    cout << "iso certificate verified for all 3385 cases\n";

    OutsideData OD;
    cerr << "Checking relaxed outside-extension SAT instances for representatives...\n";
    auto st = chrono::steady_clock::now();
    for (int c=0; c<nClasses; ++c) {
        int idx = reps[c];
        vector<array<int,6>> wit;
        bool sat = relaxedExtensionSAT(rows[idx], OD, &wit);
        cerr << "rep " << c+1 << "/" << nClasses << " case " << idx << " sat=" << sat << "\n";
        if (sat) {
            cout << "RELAXED EXTENSION FOUND for representative " << c << " case " << idx << "\n";
            for (auto B:wit) {
                for (int j=0; j<6; ++j) { if (j) cout << ' '; cout << B[j]; }
                cout << "\n";
            }
            return 1;
        }
    }
    double sec = chrono::duration<double>(chrono::steady_clock::now()-st).count();
    cout << "all representative relaxed outside-extension instances are UNSAT\n";
    cout << "VERIFIED: no 20-block 3-cover of 13 points by 6-subsets exists\n";
    cout << "Therefore C(13,6,3) >= 21. rep_check_seconds=" << fixed << setprecision(3) << sec << "\n";
    return 0;
}
