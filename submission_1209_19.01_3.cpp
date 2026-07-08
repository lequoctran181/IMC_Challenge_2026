#include <bits/stdc++.h>
using namespace std;

struct Vec3 {
    double x = 0, y = 0, z = 0;
};
static inline Vec3 operator+(const Vec3 &a, const Vec3 &b) { return {a.x + b.x, a.y + b.y, a.z + b.z}; }
static inline Vec3 operator-(const Vec3 &a, const Vec3 &b) { return {a.x - b.x, a.y - b.y, a.z - b.z}; }
static inline Vec3 operator*(const Vec3 &a, double s) { return {a.x * s, a.y * s, a.z * s}; }
static inline Vec3 operator/(const Vec3 &a, double s) { return {a.x / s, a.y / s, a.z / s}; }
static inline double dot3(const Vec3 &a, const Vec3 &b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
static inline Vec3 cross3(const Vec3 &a, const Vec3 &b) {
    return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
}
static inline double norm2(const Vec3 &a) { return dot3(a, a); }
static inline double norm3(const Vec3 &a) { return sqrt(norm2(a)); }
static inline Vec3 unit3(Vec3 a) { double n = norm3(a); return n > 1e-300 ? a * (1.0 / n) : Vec3{0,0,0}; }
static inline double clampd(double x, double lo, double hi) { return x < lo ? lo : (x > hi ? hi : x); }

struct Face { int v[3]; };

static int N0 = 0, M0 = 0;
static vector<Vec3> origP;
static vector<Face> origF;
static Vec3 bbMin, bbMax, centerP;
static double diagLen = 1.0;
static double orientSign = 1.0;

static int curN = 0;
static vector<Vec3> P;
static vector<Face> F;
static vector<unsigned char> aliveV, aliveF;
static vector<vector<int>> adj;
static vector<double> coverRad;
static int aliveFaceCount = 0;
static chrono::steady_clock::time_point T0;

static inline double elapsed() {
    return chrono::duration<double>(chrono::steady_clock::now() - T0).count();
}

struct FastInput {
    vector<char> buf;
    char *p = nullptr;
    FastInput() {
        buf.reserve(1 << 27);
        char chunk[1 << 16];
        size_t n;
        while ((n = fread(chunk, 1, sizeof(chunk), stdin)) > 0) buf.insert(buf.end(), chunk, chunk + n);
        buf.push_back('\0');
        p = buf.data();
    }
    inline void skip() { while (*p == ' ' || *p == '\n' || *p == '\r' || *p == '\t') ++p; }
    long nextLong() { skip(); return strtol(p, &p, 10); }
    double nextDouble() { skip(); return strtod(p, &p); }
    char nextChar() { skip(); return *p++; }
};

static void readInput() {
    FastInput in;
    N0 = (int)in.nextLong();
    M0 = (int)in.nextLong();
    origP.resize(N0);
    P.resize(N0);
    bbMin = {1e100,1e100,1e100};
    bbMax = {-1e100,-1e100,-1e100};
    centerP = {0,0,0};
    for (int i = 0; i < N0; ++i) {
        (void)in.nextChar();
        Vec3 q{in.nextDouble(), in.nextDouble(), in.nextDouble()};
        origP[i] = P[i] = q;
        centerP = centerP + q;
        bbMin.x = min(bbMin.x, q.x); bbMin.y = min(bbMin.y, q.y); bbMin.z = min(bbMin.z, q.z);
        bbMax.x = max(bbMax.x, q.x); bbMax.y = max(bbMax.y, q.y); bbMax.z = max(bbMax.z, q.z);
    }
    if (N0 > 0) centerP = centerP / (double)N0;
    diagLen = norm3(bbMax - bbMin); if (!(diagLen > 0)) diagLen = 1.0;
    origF.resize(M0); F.resize(M0);
    for (int i = 0; i < M0; ++i) {
        (void)in.nextChar();
        int a = (int)in.nextLong() - 1, b = (int)in.nextLong() - 1, c = (int)in.nextLong() - 1;
        origF[i].v[0] = F[i].v[0] = a;
        origF[i].v[1] = F[i].v[1] = b;
        origF[i].v[2] = F[i].v[2] = c;
    }
    double s = 0;
    int stride = max(1, M0 / 200000);
    for (int i = 0; i < M0; i += stride) {
        const Face &f = origF[i];
        Vec3 cr = cross3(origP[f.v[1]] - origP[f.v[0]], origP[f.v[2]] - origP[f.v[0]]);
        Vec3 ce = (origP[f.v[0]] + origP[f.v[1]] + origP[f.v[2]]) / 3.0;
        s += dot3(cr, ce - centerP);
    }
    orientSign = (s >= 0 ? 1.0 : -1.0);
    curN = N0;
    aliveV.assign(curN, 1);
    aliveF.assign(M0, 1);
    coverRad.assign(curN, 0.0);
    aliveFaceCount = M0;
}

struct U64Hash {
    size_t operator()(uint64_t x) const noexcept {
        x += 0x9e3779b97f4a7c15ULL;
        x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
        x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
        return (size_t)(x ^ (x >> 31));
    }
};

static inline uint64_t edgeKey(int a, int b) {
    if (a > b) swap(a, b);
    return (uint64_t)(uint32_t)a << 32 | (uint32_t)b;
}
static inline bool faceHas(const Face &f, int v) { return f.v[0] == v || f.v[1] == v || f.v[2] == v; }
static inline int thirdOf(const Face &f, int a, int b) {
    for (int k = 0; k < 3; ++k) { int x = f.v[k]; if (x != a && x != b) return x; }
    return -1;
}
static inline Vec3 faceCrossCur(const Face &f) { return cross3(P[f.v[1]] - P[f.v[0]], P[f.v[2]] - P[f.v[0]]); }
static inline Vec3 faceCrossVec(const vector<Vec3> &V, const Face &f) { return cross3(V[f.v[1]] - V[f.v[0]], V[f.v[2]] - V[f.v[0]]); }

static void rebuildAdj() {
    vector<int> deg(curN, 0);
    int cntAlive = 0;
    for (int i = 0; i < (int)F.size(); ++i) {
        if (!aliveF[i]) continue;
        Face f = F[i];
        bool ok = f.v[0] >= 0 && f.v[0] < curN && f.v[1] >= 0 && f.v[1] < curN && f.v[2] >= 0 && f.v[2] < curN &&
                  aliveV[f.v[0]] && aliveV[f.v[1]] && aliveV[f.v[2]] && f.v[0] != f.v[1] && f.v[0] != f.v[2] && f.v[1] != f.v[2];
        if (!ok) { aliveF[i] = 0; continue; }
        deg[f.v[0]]++; deg[f.v[1]]++; deg[f.v[2]]++; cntAlive++;
    }
    aliveFaceCount = cntAlive;
    adj.assign(curN, {});
    for (int i = 0; i < curN; ++i) adj[i].reserve(deg[i]);
    for (int i = 0; i < (int)F.size(); ++i) if (aliveF[i]) {
        Face f = F[i];
        adj[f.v[0]].push_back(i); adj[f.v[1]].push_back(i); adj[f.v[2]].push_back(i);
    }
}

static int activeVertices() { int c = 0; for (unsigned char b : aliveV) c += b != 0; return c; }

static bool sameFaceUnordered(const Face &f, int a, int b, int c) {
    int x[3] = {f.v[0], f.v[1], f.v[2]}; int y[3] = {a,b,c};
    sort(x, x + 3); sort(y, y + 3);
    return x[0] == y[0] && x[1] == y[1] && x[2] == y[2];
}

static bool findEdgeFaces(int u, int v, int ef[2], int opp[2]) {
    if (u < 0 || v < 0 || u >= curN || v >= curN || !aliveV[u] || !aliveV[v]) return false;
    int cnt = 0;
    const vector<int> &A = adj[u].size() < adj[v].size() ? adj[u] : adj[v];
    for (int fid : A) {
        if (!aliveF[fid]) continue;
        const Face &f = F[fid];
        if (!faceHas(f, u) || !faceHas(f, v)) continue;
        if (cnt >= 2) return false;
        ef[cnt] = fid; opp[cnt] = thirdOf(f, u, v); cnt++;
    }
    if (cnt != 2) return false;
    if (opp[0] < 0 || opp[1] < 0 || opp[0] == opp[1]) return false;
    return true;
}

static vector<int> markA, markB;
static int stampA = 1, stampB = 1;

static bool linkOK(int u, int v, const int opp[2]) {
    if (++stampA > 2000000000) { fill(markA.begin(), markA.end(), 0); stampA = 1; }
    if (++stampB > 2000000000) { fill(markB.begin(), markB.end(), 0); stampB = 1; }
    for (int fid : adj[u]) {
        if (!aliveF[fid]) continue;
        const Face &f = F[fid]; if (!faceHas(f, u)) continue;
        for (int k = 0; k < 3; ++k) { int x = f.v[k]; if (x != u && x != v) markA[x] = stampA; }
    }
    int common = 0;
    for (int fid : adj[v]) {
        if (!aliveF[fid]) continue;
        const Face &f = F[fid]; if (!faceHas(f, v)) continue;
        for (int k = 0; k < 3; ++k) {
            int x = f.v[k];
            if (x == u || x == v || markA[x] != stampA) continue;
            if (x != opp[0] && x != opp[1]) return false;
            if (markB[x] != stampB) { markB[x] = stampB; common++; }
        }
    }
    return common == 2 && markB[opp[0]] == stampB && markB[opp[1]] == stampB;
}

struct CollapseParams {
    double epsR = 0.049;
    double planeR = 0.020;
    double minCos = 0.90;
    double areaR = 1e-24;
};
struct EvalCollapse {
    bool ok = false;
    double cost = 1e100;
    double newRad = 0;
    int keep = -1, rem = -1;
};

static bool duplicateFaceAfter(int fid, int rem, int ef0, int ef1, int a, int b, int c) {
    int probe = a;
    if ((int)adj[b].size() < (int)adj[probe].size()) probe = b;
    if ((int)adj[c].size() < (int)adj[probe].size()) probe = c;
    for (int g : adj[probe]) {
        if (!aliveF[g] || g == fid || g == ef0 || g == ef1) continue;
        if (faceHas(F[g], rem)) continue;
        if (sameFaceUnordered(F[g], a, b, c)) return true;
    }
    return false;
}

static EvalCollapse evalEndpoint(int keep, int rem, const int ef[2], const CollapseParams &par) {
    EvalCollapse r; r.keep = keep; r.rem = rem;
    const double eps = par.epsR * diagLen;
    const double planeTol = par.planeR * diagLen;
    double dkr = norm3(P[keep] - P[rem]);
    r.newRad = max(coverRad[keep], coverRad[rem] + dkr);
    if (r.newRad > eps + 1e-12) return r;
    const double minArea2 = max(1e-300, par.areaR * diagLen * diagLen * diagLen * diagLen);
    double maxPlane = 0.0, maxAng = 0.0, areaLoss = 0.0;
    int changed = 0;
    for (int fid : adj[rem]) {
        if (!aliveF[fid] || !faceHas(F[fid], rem)) continue;
        if (fid == ef[0] || fid == ef[1]) continue;
        if (faceHas(F[fid], keep)) return r;
        Face old = F[fid], nf = old;
        for (int k = 0; k < 3; ++k) if (nf.v[k] == rem) nf.v[k] = keep;
        if (nf.v[0] == nf.v[1] || nf.v[0] == nf.v[2] || nf.v[1] == nf.v[2]) return r;
        Vec3 co = faceCrossCur(old), cn = faceCrossCur(nf);
        double ao = norm3(co), an = norm3(cn);
        if (!(ao > 0) || !(an > 0) || norm2(cn) <= minArea2) return r;
        double nd = dot3(co, cn) / (ao * an); nd = clampd(nd, -1.0, 1.0);
        if (nd < par.minCos) return r;
        Vec3 no = co * (1.0 / ao);
        double pd = fabs(dot3(no, P[keep] - P[old.v[0]]));
        if (pd > planeTol) return r;
        if (duplicateFaceAfter(fid, rem, ef[0], ef[1], nf.v[0], nf.v[1], nf.v[2])) return r;
        maxPlane = max(maxPlane, pd);
        maxAng = max(maxAng, 1.0 - nd);
        areaLoss = max(areaLoss, max(0.0, 1.0 - an / ao));
        changed++;
    }
    if (changed == 0) return r;
    r.ok = true;
    r.cost = 1.1 * (r.newRad / (eps + 1e-300)) + 0.9 * (maxPlane / (planeTol + 1e-300)) +
             180.0 * maxAng + 0.04 * areaLoss + 0.0002 * changed;
    return r;
}

static void applyCollapse(const EvalCollapse &e, const int ef[2]) {
    int keep = e.keep, rem = e.rem;
    for (int i = 0; i < 2; ++i) if (ef[i] >= 0 && aliveF[ef[i]]) { aliveF[ef[i]] = 0; aliveFaceCount--; }
    for (int fid : adj[rem]) {
        if (!aliveF[fid] || !faceHas(F[fid], rem)) continue;
        for (int k = 0; k < 3; ++k) if (F[fid].v[k] == rem) F[fid].v[k] = keep;
    }
    aliveV[rem] = 0;
    coverRad[keep] = e.newRad;
    vector<int> merged;
    merged.reserve(adj[keep].size() + adj[rem].size());
    for (int fid : adj[keep]) if (aliveF[fid] && faceHas(F[fid], keep)) merged.push_back(fid);
    for (int fid : adj[rem]) if (aliveF[fid] && faceHas(F[fid], keep)) merged.push_back(fid);
    sort(merged.begin(), merged.end()); merged.erase(unique(merged.begin(), merged.end()), merged.end());
    adj[keep].swap(merged);
    vector<int>().swap(adj[rem]);
}

static bool tryCollapseEdge(int u, int v, const CollapseParams &par) {
    if (u == v || u < 0 || v < 0 || u >= curN || v >= curN || !aliveV[u] || !aliveV[v]) return false;
    int ef[2] = {-1,-1}, opp[2] = {-1,-1};
    if (!findEdgeFaces(u, v, ef, opp)) return false;
    if (!linkOK(u, v, opp)) return false;
    EvalCollapse a = evalEndpoint(u, v, ef, par);
    EvalCollapse b = evalEndpoint(v, u, ef, par);
    if (!a.ok && !b.ok) return false;
    if (b.ok && (!a.ok || b.cost < a.cost)) applyCollapse(b, ef);
    else applyCollapse(a, ef);
    return true;
}

struct SmoothStats { double smooth = 0, coarse = 0, sharp = 0, verySharp = 0, bad = 0; int samples = 0; };
static SmoothStats analyzeSmooth() {
    SmoothStats s;
    if (N0 < 50 || M0 < 50) return s;
    rebuildAdj();
    int stride = max(1, M0 / 60000), limit = 140000;
    const double c10 = cos(10.0 * M_PI / 180.0), c30 = cos(30.0 * M_PI / 180.0);
    const double c22 = cos(22.0 * M_PI / 180.0), c45 = cos(45.0 * M_PI / 180.0);
    int sm = 0, co = 0, sh = 0, vs = 0, bad = 0, tot = 0;
    for (int fid = 0; fid < M0 && tot < limit; fid += stride) {
        const Face &f = F[fid];
        int e[3][2] = {{f.v[0],f.v[1]}, {f.v[1],f.v[2]}, {f.v[2],f.v[0]}};
        for (int k = 0; k < 3 && tot < limit; ++k) {
            int ef[2], op[2];
            if (!findEdgeFaces(e[k][0], e[k][1], ef, op)) { bad++; continue; }
            Vec3 n0 = unit3(faceCrossCur(F[ef[0]]));
            Vec3 n1 = unit3(faceCrossCur(F[ef[1]]));
            double d = clampd(dot3(n0, n1), -1.0, 1.0);
            tot++;
            if (d > c10) sm++;
            if (d > c30) co++;
            if (d < c22) sh++;
            if (d < c45) vs++;
        }
    }
    s.samples = tot;
    if (tot) { s.smooth = (double)sm/tot; s.coarse = (double)co/tot; s.sharp = (double)sh/tot; s.verySharp = (double)vs/tot; }
    if (tot + bad) s.bad = (double)bad / (tot + bad);
    return s;
}

static void collectEdges(vector<uint64_t> &keys, double maxL2) {
    keys.clear(); keys.reserve((size_t)aliveFaceCount * 3 + 16);
    for (int i = 0; i < (int)F.size(); ++i) {
        if (!aliveF[i]) continue;
        const Face &f = F[i];
        if (f.v[0] < 0 || f.v[1] < 0 || f.v[2] < 0 || f.v[0] >= curN || f.v[1] >= curN || f.v[2] >= curN) continue;
        if (!aliveV[f.v[0]] || !aliveV[f.v[1]] || !aliveV[f.v[2]]) continue;
        uint64_t e0 = edgeKey(f.v[0], f.v[1]), e1 = edgeKey(f.v[1], f.v[2]), e2 = edgeKey(f.v[2], f.v[0]);
        keys.push_back(e0); keys.push_back(e1); keys.push_back(e2);
    }
    sort(keys.begin(), keys.end()); keys.erase(unique(keys.begin(), keys.end()), keys.end());
    if (maxL2 < 1e99) {
        size_t w = 0;
        for (uint64_t k : keys) {
            int a = (int)(k >> 32), b = (int)(uint32_t)k;
            if (a >= 0 && b >= 0 && a < curN && b < curN && aliveV[a] && aliveV[b] && norm2(P[a] - P[b]) <= maxL2) keys[w++] = k;
        }
        keys.resize(w);
    }
}
struct EdgeRec { float l2; uint64_t k; };

static bool edgePass(const CollapseParams &par, double timeLimit, int maxRounds, int targetVertices) {
    bool any = false;
    vector<uint64_t> keys;
    vector<EdgeRec> edges;
    double eps = par.epsR * diagLen;
    double maxL2 = (eps * 1.0005) * (eps * 1.0005);
    for (int round = 0; round < maxRounds && elapsed() < timeLimit && activeVertices() > targetVertices; ++round) {
        collectEdges(keys, maxL2);
        if (keys.empty()) break;
        edges.clear(); edges.reserve(keys.size());
        for (uint64_t k : keys) {
            int a = (int)(k >> 32), b = (int)(uint32_t)k;
            if (a >= 0 && b >= 0 && a < curN && b < curN && aliveV[a] && aliveV[b]) edges.push_back({(float)norm2(P[a] - P[b]), k});
        }
        sort(edges.begin(), edges.end(), [](const EdgeRec &a, const EdgeRec &b){ return a.l2 < b.l2; });
        int hit = 0;
        for (size_t i = 0; i < edges.size() && elapsed() < timeLimit && activeVertices() > targetVertices; ++i) {
            int a = (int)(edges[i].k >> 32), b = (int)(uint32_t)edges[i].k;
            if (a >= 0 && b >= 0 && a < curN && b < curN && aliveV[a] && aliveV[b]) if (tryCollapseEdge(a,b,par)) { hit++; any = true; }
        }
        if (hit == 0) break;
        if (round + 1 < maxRounds) rebuildAdj();
    }
    return any;
}

struct PointGrid {
    Vec3 mn, mx;
    double cell = 1.0;
    int nx = 1, ny = 1, nz = 1;
    vector<vector<int>> buck;
    const vector<Vec3> *pts = nullptr;
    int clampi(int x, int n) const { return x < 0 ? 0 : (x >= n ? n - 1 : x); }
    int ix(double x) const { return clampi((int)floor((x - mn.x) / cell), nx); }
    int iy(double y) const { return clampi((int)floor((y - mn.y) / cell), ny); }
    int iz(double z) const { return clampi((int)floor((z - mn.z) / cell), nz); }
    int key(int x, int y, int z) const { return (z * ny + y) * nx + x; }
    void build(const vector<Vec3> &p, double c) {
        pts = &p; cell = max(c, 1e-9);
        mn = {1e100,1e100,1e100}; mx = {-1e100,-1e100,-1e100};
        for (const Vec3 &q : p) { mn.x = min(mn.x,q.x); mn.y = min(mn.y,q.y); mn.z = min(mn.z,q.z); mx.x = max(mx.x,q.x); mx.y = max(mx.y,q.y); mx.z = max(mx.z,q.z); }
        mn.x -= cell; mn.y -= cell; mn.z -= cell; mx.x += cell; mx.y += cell; mx.z += cell;
        nx = max(1, (int)((mx.x - mn.x) / cell) + 1); ny = max(1, (int)((mx.y - mn.y) / cell) + 1); nz = max(1, (int)((mx.z - mn.z) / cell) + 1);
        long long tot = 1LL * nx * ny * nz;
        if (tot > 2500000) { nx = ny = nz = 1; cell = max({mx.x-mn.x, mx.y-mn.y, mx.z-mn.z, 1e-9}); }
        buck.assign((size_t)nx * ny * nz, {});
        for (int i = 0; i < (int)p.size(); ++i) buck[key(ix(p[i].x), iy(p[i].y), iz(p[i].z))].push_back(i);
    }
    bool nearPoint(const Vec3 &q, double r2) const {
        int X = ix(q.x), Y = iy(q.y), Z = iz(q.z);
        for (int dz = -1; dz <= 1; ++dz) { int z = Z + dz; if (z < 0 || z >= nz) continue;
            for (int dy = -1; dy <= 1; ++dy) { int y = Y + dy; if (y < 0 || y >= ny) continue;
                for (int dx = -1; dx <= 1; ++dx) { int x = X + dx; if (x < 0 || x >= nx) continue;
                    for (int id : buck[key(x,y,z)]) if (norm2((*pts)[id] - q) <= r2) return true;
                }
            }
        }
        return false;
    }
    int nearestPoint(const Vec3 &q, double r2) const {
        int X = ix(q.x), Y = iy(q.y), Z = iz(q.z);
        int best = -1; double bd = r2;
        for (int rad = 1; rad <= 3 && best < 0; ++rad) {
            for (int dz = -rad; dz <= rad; ++dz) { int z = Z + dz; if (z < 0 || z >= nz) continue;
                for (int dy = -rad; dy <= rad; ++dy) { int y = Y + dy; if (y < 0 || y >= ny) continue;
                    for (int dx = -rad; dx <= rad; ++dx) { int x = X + dx; if (x < 0 || x >= nx) continue;
                        for (int id : buck[key(x,y,z)]) { double d = norm2((*pts)[id] - q); if (d <= bd) { bd = d; best = id; } }
                    }
                }
            }
        }
        return best;
    }
};

static PointGrid origGridCache;
static bool origGridReady = false;
static double origGridCell = -1;
static void ensureOrigGrid(double eps) {
    if (!origGridReady || fabs(origGridCell - eps) > eps * 1e-6) { origGridCache.build(origP, eps); origGridReady = true; origGridCell = eps; }
}

static bool coverageOK(const vector<Vec3> &cand, double scale = 0.0495, double timeLimit = 18.9) {
    if (cand.empty() || (int)cand.size() > N0) return false;
    double eps = scale * diagLen;
    double r2 = eps * eps * (1.0 + 1e-12);
    PointGrid gCand; gCand.build(cand, eps);
    for (int i = 0; i < N0; ++i) {
        if ((i & 32767) == 0 && elapsed() > timeLimit) return false;
        if (!gCand.nearPoint(origP[i], r2)) return false;
    }
    ensureOrigGrid(eps);
    for (int i = 0; i < (int)cand.size(); ++i) {
        if (!origGridCache.nearPoint(cand[i], r2)) return false;
    }
    return true;
}

static bool candidateFacesBasicOK(const vector<Vec3> &X, const vector<Face> &Y) {
    if (X.empty() || Y.empty() || (int)X.size() > N0) return false;
    double minA2 = max(1e-300, 1e-28 * diagLen * diagLen * diagLen * diagLen);
    unordered_map<uint64_t, int, U64Hash> ec;
    ec.reserve(Y.size() * 3 * 2 + 16);
    for (const Face &f : Y) {
        int a = f.v[0], b = f.v[1], c = f.v[2];
        if (a < 0 || b < 0 || c < 0 || a >= (int)X.size() || b >= (int)X.size() || c >= (int)X.size()) return false;
        if (a == b || a == c || b == c) return false;
        if (norm2(faceCrossVec(X, f)) <= minA2) return false;
        ec[edgeKey(a,b)]++; ec[edgeKey(b,c)]++; ec[edgeKey(c,a)]++;
    }
    for (auto &kv : ec) if (kv.second != 2) return false;
    return true;
}

static void installCandidate(const vector<Vec3> &X, const vector<Face> &Y) {
    curN = (int)X.size(); P = X; F = Y;
    aliveV.assign(curN, 1); aliveF.assign(F.size(), 1); coverRad.assign(curN, 0.0); aliveFaceCount = (int)F.size();
    rebuildAdj(); markA.assign(curN, 0); markB.assign(curN, 0); stampA = stampB = 1;
}

struct RenderMap { vector<float> dep, nx, ny, nz; vector<unsigned char> fg; };
static unordered_map<int, vector<RenderMap>> origRenderCache;

static inline void projectPoint(const Vec3 &p, int view, int res, double &u, double &v, double &z) {
    const double D = 2.5;
    double f = 800.0 * (double)res / 1024.0;
    double c = 0.5 * (double)res;
    double sx, sy;
    if (view == 0) { sx = p.y; sy = p.z; z = D - p.x; }
    else if (view == 1) { sx = -p.y; sy = p.z; z = D + p.x; }
    else if (view == 2) { sx = -p.x; sy = p.z; z = D - p.y; }
    else if (view == 3) { sx = p.x; sy = p.z; z = D + p.y; }
    else if (view == 4) { sx = p.x; sy = p.y; z = D - p.z; }
    else { sx = -p.x; sy = p.y; z = D + p.z; }
    u = f * sx / z + c; v = f * sy / z + c;
}

static void rasterTri(RenderMap &rm, int res, const Vec3 &a, const Vec3 &b, const Vec3 &c, const Vec3 &n, int view) {
    double x0,y0,z0,x1,y1,z1,x2,y2,z2;
    projectPoint(a, view, res, x0,y0,z0); projectPoint(b, view, res, x1,y1,z1); projectPoint(c, view, res, x2,y2,z2);
    if (z0 <= 0 || z1 <= 0 || z2 <= 0) return;
    int xmin = max(0, (int)floor(min(x0, min(x1,x2)) - 0.5));
    int xmax = min(res - 1, (int)ceil(max(x0, max(x1,x2)) + 0.5));
    int ymin = max(0, (int)floor(min(y0, min(y1,y2)) - 0.5));
    int ymax = min(res - 1, (int)ceil(max(y0, max(y1,y2)) + 0.5));
    if (xmin > xmax || ymin > ymax) return;
    double den = (y1-y2)*(x0-x2) + (x2-x1)*(y0-y2);
    if (fabs(den) < 1e-18) return;
    for (int yy = ymin; yy <= ymax; ++yy) {
        double py = yy + 0.5;
        for (int xx = xmin; xx <= xmax; ++xx) {
            double px = xx + 0.5;
            double w0 = ((y1-y2)*(px-x2) + (x2-x1)*(py-y2)) / den;
            double w1 = ((y2-y0)*(px-x2) + (x0-x2)*(py-y2)) / den;
            double w2 = 1.0 - w0 - w1;
            if (w0 < -1e-9 || w1 < -1e-9 || w2 < -1e-9) continue;
            double zp = 1.0 / (w0 / z0 + w1 / z1 + w2 / z2);
            int id = yy * res + xx;
            if (zp < rm.dep[id]) { rm.dep[id] = (float)zp; rm.nx[id] = (float)n.x; rm.ny[id] = (float)n.y; rm.nz[id] = (float)n.z; rm.fg[id] = 1; }
        }
    }
}

static RenderMap renderMesh(const vector<Vec3> &V, const vector<Face> &Y, int res, int view) {
    RenderMap rm; int S = res * res;
    rm.dep.assign(S, 255.0f); rm.nx.assign(S, 0); rm.ny.assign(S, 0); rm.nz.assign(S, 0); rm.fg.assign(S, 0);
    for (const Face &f : Y) {
        if (f.v[0] < 0 || f.v[1] < 0 || f.v[2] < 0 || f.v[0] >= (int)V.size() || f.v[1] >= (int)V.size() || f.v[2] >= (int)V.size()) continue;
        Vec3 cr = cross3(V[f.v[1]] - V[f.v[0]], V[f.v[2]] - V[f.v[0]]);
        double l = norm3(cr); if (!(l > 0)) continue;
        rasterTri(rm, res, V[f.v[0]], V[f.v[1]], V[f.v[2]], cr * (1.0/l), view);
    }
    return rm;
}

static const vector<RenderMap> &getOrigRenders(int res) {
    auto it = origRenderCache.find(res); if (it != origRenderCache.end()) return it->second;
    vector<RenderMap> v; v.reserve(6);
    for (int view = 0; view < 6; ++view) v.push_back(renderMesh(origP, origF, res, view));
    auto ins = origRenderCache.emplace(res, std::move(v)); return ins.first->second;
}

static inline double valAt(const RenderMap &m, int idx, int channel) {
    if (channel == 0) return (m.nx[idx] + 1.0) * 127.5;
    if (channel == 1) return (m.ny[idx] + 1.0) * 127.5;
    if (channel == 2) return (m.nz[idx] + 1.0) * 127.5;
    return m.dep[idx];
}

static double ssimSlow(const RenderMap &a, const RenderMap &b, const vector<unsigned char> &fg, int res, int channel) {
    const int rad = 5;
    const double c1 = (0.01 * 255.0) * (0.01 * 255.0);
    const double c2 = (0.03 * 255.0) * (0.03 * 255.0);
    double total = 0; int cnt = 0;
    for (int y = 0; y < res; ++y) for (int x = 0; x < res; ++x) {
        int center = y * res + x; if (!fg[center]) continue;
        double sx=0, sy=0, sxx=0, syy=0, sxy=0; int n=0;
        for (int dy = -rad; dy <= rad; ++dy) {
            int yy = y + dy; if (yy < 0) yy = 0; if (yy >= res) yy = res - 1;
            for (int dx = -rad; dx <= rad; ++dx) {
                int xx = x + dx; if (xx < 0) xx = 0; if (xx >= res) xx = res - 1;
                int id = yy * res + xx; double vx = valAt(a, id, channel), vy = valAt(b, id, channel);
                sx += vx; sy += vy; sxx += vx*vx; syy += vy*vy; sxy += vx*vy; n++;
            }
        }
        double inv = 1.0 / n, ux = sx * inv, uy = sy * inv;
        double vx = max(0.0, sxx * inv - ux * ux), vy = max(0.0, syy * inv - uy * uy), cov = sxy * inv - ux * uy;
        double den = (ux*ux + uy*uy + c1) * (vx + vy + c2);
        if (den <= 0) continue;
        total += ((2*ux*uy + c1) * (2*cov + c2)) / den; cnt++;
    }
    return cnt ? total / cnt : 1.0;
}

static double visualScoreCandidate(const vector<Vec3> &X, const vector<Face> &Y, int res, double timeLimit = 19.2) {
    if (elapsed() > timeLimit) return 0.0;
    const vector<RenderMap> &orig = getOrigRenders(res);
    double total = 0;
    for (int view = 0; view < 6; ++view) {
        if (elapsed() > timeLimit) return 0.0;
        RenderMap simp = renderMesh(X, Y, res, view);
        vector<unsigned char> fg(res * res, 0);
        for (int i = 0; i < res * res; ++i) fg[i] = (orig[view].fg[i] || simp.fg[i]) ? 1 : 0;
        double ns = 0;
        for (int ch = 0; ch < 3; ++ch) ns += ssimSlow(orig[view], simp, fg, res, ch);
        ns /= 3.0;
        double ds = ssimSlow(orig[view], simp, fg, res, 3);
        total += 0.5 * ns + 0.5 * ds;
    }
    return total / 6.0;
}

static vector<Vec3> originalVertexNormals() {
    vector<Vec3> n(N0, {0,0,0});
    for (const Face &f : origF) {
        Vec3 cr = cross3(origP[f.v[1]] - origP[f.v[0]], origP[f.v[2]] - origP[f.v[0]]);
        n[f.v[0]] = n[f.v[0]] + cr; n[f.v[1]] = n[f.v[1]] + cr; n[f.v[2]] = n[f.v[2]] + cr;
    }
    for (int i = 0; i < N0; ++i) {
        if (norm2(n[i]) < 1e-30) n[i] = unit3(origP[i] - centerP) * orientSign;
        else n[i] = unit3(n[i]);
    }
    return n;
}

static void orientFaceByRef(vector<Face> &Y, const vector<Vec3> &X, Face f, Vec3 ref) {
    Vec3 cr = faceCrossVec(X, f);
    if (dot3(cr, ref) < 0) swap(f.v[1], f.v[2]);
    Y.push_back(f);
}
static void addOriented(vector<Face> &Y, const vector<Vec3> &X, Face f, Vec3 outward) {
    Vec3 cr = faceCrossVec(X, f);
    if (dot3(cr, outward) * orientSign < 0) swap(f.v[1], f.v[2]);
    Y.push_back(f);
}

static bool acceptCandidate(const vector<Vec3> &X, const vector<Face> &Y, double visualNeed, int res, double covScale = 0.0495) {
    if (X.empty() || (int)X.size() >= activeVertices() || (int)X.size() > N0) return false;
    if (!candidateFacesBasicOK(X, Y)) return false;
    if (!coverageOK(X, covScale, 18.6)) return false;
    if (visualNeed > 0 && elapsed() < 17.8) {
        double q = visualScoreCandidate(X, Y, res, 19.1);
        if (q < visualNeed) return false;
    }
    installCandidate(X, Y);
    return true;
}

static bool adj3small(const int a[3], int m, int &base) {
    for (int t = 0; t < 3; ++t) for (int s = 0; s < 2; ++s) {
        int x = (a[t] - s + m) % m; bool ok = true;
        for (int i = 0; i < 3; ++i) { int d = (a[i] - x + m) % m; if (d != 0 && d != 1) { ok = false; break; } }
        if (ok) { base = x; return true; }
    }
    return false;
}

static bool tubeFaceCompatible(const Face &f, int S) {
    if (S < 8 || N0 % S) return false;
    int U = N0 / S; if (U < 8) return false;
    int a[3] = {f.v[0] / S, f.v[1] / S, f.v[2] / S};
    int c[3] = {f.v[0] % S, f.v[1] % S, f.v[2] % S};
    int ra = 0, ca = 0; if (!adj3small(a, U, ra) || !adj3small(c, S, ca)) return false;
    int mask = 0;
    for (int i = 0; i < 3; ++i) {
        int x = (a[i] - ra + U) % U, y = (c[i] - ca + S) % S;
        if (x > 1 || y > 1) return false;
        mask |= 1 << (x * 2 + y);
    }
    return __builtin_popcount((unsigned)mask) == 3;
}

static vector<int> tubeSCandidates() {
    map<int,int> mp;
    int st = max(1, M0 / 120000);
    for (int i = 0; i < M0; i += st) {
        const Face &f = origF[i]; int a[3] = {f.v[0], f.v[1], f.v[2]};
        for (int k = 0; k < 3; ++k) {
            int d = abs(a[k] - a[(k+1)%3]); if (!d) continue; d = min(d, N0 - d);
            if (d >= 6 && d <= N0 / 3) mp[d]++;
        }
    }
    vector<pair<int,int>> q; for (auto &p : mp) q.push_back({p.second,p.first}); sort(q.rbegin(), q.rend());
    vector<int> r;
    auto add = [&](int s){ if (s >= 8 && s <= N0 / 4 && N0 % s == 0 && find(r.begin(), r.end(), s) == r.end()) r.push_back(s); };
    for (int i = 0; i < (int)q.size() && i < 20; ++i) { int d = q[i].second; for (int e = -4; e <= 4; ++e) add(d + e); if (d) add(N0 / d); }
    for (int s : {16,24,32,40,48,64,80,96,100,112,120,128,144,160,192,200,224,240,256,320,384,512}) add(s);
    return r;
}
static bool tubeGoodS(int S) {
    if (N0 % S) return false;
    int U = N0 / S; if (U < 8 || S < 8) return false;
    int st = max(1, M0 / 180000), tot = 0, ok = 0;
    for (int i = 0; i < M0; i += st) { tot++; if (tubeFaceCompatible(origF[i], S)) ok++; }
    return tot > 80 && ok * 1000 >= tot * 995;
}

static bool makeTubeGrid(int S, int U2, int S2, vector<Vec3> &X, vector<Face> &Y, const vector<Vec3> &vn) {
    if (N0 % S) return false; int U = N0 / S;
    if (U2 < 8 || S2 < 8 || U2 > U || S2 > S) return false;
    X.clear(); Y.clear(); X.reserve((size_t)U2 * S2);
    vector<int> src(U2 * S2);
    for (int i = 0; i < U2; ++i) { int oi = (long long)i * U / U2;
        for (int j = 0; j < S2; ++j) { int oj = (long long)j * S / S2; int sid = oi * S + oj; src[i*S2+j] = sid; X.push_back(origP[sid]); }
    }
    auto id = [&](int i, int j){ return ((i % U2 + U2) % U2) * S2 + ((j % S2 + S2) % S2); };
    Y.reserve((size_t)2 * U2 * S2);
    for (int i = 0; i < U2; ++i) for (int j = 0; j < S2; ++j) {
        int a = id(i,j), b = id(i+1,j), c = id(i+1,j+1), d = id(i,j+1);
        orientFaceByRef(Y, X, Face{{a,b,d}}, vn[src[a]] + vn[src[b]] + vn[src[d]]);
        orientFaceByRef(Y, X, Face{{b,c,d}}, vn[src[b]] + vn[src[c]] + vn[src[d]]);
    }
    return true;
}

static bool tryTubeStructured() {
    if (!(M0 == 2 * N0 && N0 >= 8000) || elapsed() > 10.5) return false;
    vector<int> ss = tubeSCandidates(); int S = 0;
    for (int s : ss) if (tubeGoodS(s)) { S = s; break; }
    if (!S) return false;
    int U = N0 / S;
    vector<Vec3> vn = originalVertexNormals();
    vector<pair<int,int>> dims;
    double eps = 0.049 * diagLen;
    double approxU = 0, approxS = 0;
    int sampleU = min(U, 512), sampleS = min(S, 512);
    for (int i = 0; i < sampleU; ++i) {
        int a = (long long)i * U / sampleU * S;
        int b = ((long long)(i+1) * U / sampleU % U) * S;
        approxU += norm3(origP[a] - origP[b]);
    }
    for (int j = 0; j < sampleS; ++j) approxS += norm3(origP[j] - origP[(j+1)%S]);
    int baseU = max(8, (int)ceil(approxU / (eps * 1.25)));
    int baseS = max(8, (int)ceil(approxS / (eps * 1.25)));
    for (double m : {0.65,0.8,1.0,1.25,1.6,2.0,2.6,3.2}) {
        int u2 = min(U, max(8, (int)ceil(baseU * m)));
        int s2 = min(S, max(8, (int)ceil(baseS * m)));
        if (u2 * s2 < N0) dims.push_back({u2,s2});
    }
    sort(dims.begin(), dims.end(), [](auto a, auto b){ return a.first*a.second < b.first*b.second; });
    dims.erase(unique(dims.begin(), dims.end()), dims.end());
    vector<Vec3> X; vector<Face> Y;
    for (auto [u2,s2] : dims) {
        if (elapsed() > 17.2) break;
        if (!makeTubeGrid(S, u2, s2, X, Y, vn)) continue;
        if (acceptCandidate(X, Y, 0.935, N0 < 70000 ? 160 : 112, 0.0492)) return true;
    }
    return false;
}

static inline bool sameOrigFace(int fid, int a, int b, int c) { return sameFaceUnordered(origF[fid], a, b, c); }
static bool detectPolarGrid(int &R, int &V) {
    if (N0 < 300 || M0 != 2 * (N0 - 2)) return false;
    vector<int> cand;
    for (int v = 8; v <= min(4096, N0 - 2); ++v) if ((N0 - 2) % v == 0) cand.push_back(v);
    sort(cand.begin(), cand.end(), [&](int a, int b){ return abs(a - (int)sqrt(N0)) < abs(b - (int)sqrt(N0)); });
    for (int v : cand) {
        int r = (N0 - 2) / v; if (r < 3) continue;
        int step = max(1, v / 96); bool ok = true;
        for (int j = 0; j < v && ok; j += step) {
            int a = 1 + j, b = 1 + (j + 1) % v;
            if (!sameOrigFace(j, 0, b, a)) ok = false;
            int off = v + 2 * (r - 1) * v + j;
            int c = 1 + (r - 1) * v + j, d = 1 + (r - 1) * v + (j + 1) % v;
            if (ok && !sameOrigFace(off, N0 - 1, c, d)) ok = false;
        }
        int span = max(1, (r - 1) * v / 300);
        for (int q = 0; q < (r - 1) * v && ok; q += span) {
            int rr = q / v, j = q - rr * v;
            int a = 1 + rr * v + j, b = 1 + rr * v + (j + 1) % v;
            int c = 1 + (rr + 1) * v + j, d = 1 + (rr + 1) * v + (j + 1) % v;
            int f = v + 2 * (rr * v + j);
            if (!sameOrigFace(f, a,b,c) || !sameOrigFace(f+1, b,d,c)) ok = false;
        }
        if (ok) { R = r; V = v; return true; }
    }
    return false;
}

static bool makePolarGrid(int R, int V, int R2, int V2, vector<Vec3> &X, vector<Face> &Y, const vector<Vec3> &vn) {
    if (R2 < 3 || V2 < 8 || 2 + R2 * V2 >= N0) return false;
    X.clear(); Y.clear(); X.reserve(2 + (size_t)R2 * V2); vector<int> src; src.reserve(2 + (size_t)R2 * V2);
    X.push_back(origP[0]); src.push_back(0);
    for (int i = 0; i < R2; ++i) {
        int oi = 1 + (int)((long long)i * (R - 1) / max(1, R2 - 1));
        for (int j = 0; j < V2; ++j) { int oj = (long long)j * V / V2; int s = 1 + (oi - 1) * V + oj; X.push_back(origP[s]); src.push_back(s); }
    }
    int bot = (int)X.size(); X.push_back(origP[N0 - 1]); src.push_back(N0 - 1);
    auto rid = [&](int r, int j){ return 1 + (r - 1) * V2 + ((j % V2 + V2) % V2); };
    auto add = [&](int a, int b, int c){ orientFaceByRef(Y, X, Face{{a,b,c}}, vn[src[a]] + vn[src[b]] + vn[src[c]]); };
    for (int j = 0; j < V2; ++j) add(0, rid(1,j+1), rid(1,j));
    for (int r = 1; r < R2; ++r) for (int j = 0; j < V2; ++j) {
        int a = rid(r,j), b = rid(r,j+1), c = rid(r+1,j), d = rid(r+1,j+1); add(a,b,c); add(b,d,c);
    }
    for (int j = 0; j < V2; ++j) add(bot, rid(R2,j), rid(R2,j+1));
    return true;
}

static bool tryPolarStructured() {
    if (elapsed() > 10.5) return false;
    int R=0,V=0; if (!detectPolarGrid(R,V)) return false;
    vector<Vec3> vn = originalVertexNormals();
    vector<pair<int,int>> dims;
    double eps = 0.049 * diagLen;
    int rBase = max(3, (int)ceil(M_PI * max({bbMax.x-bbMin.x, bbMax.y-bbMin.y, bbMax.z-bbMin.z}) * 0.5 / (eps * 1.25)));
    int vBase = max(8, (int)ceil(2.0 * M_PI * max({bbMax.x-bbMin.x, bbMax.y-bbMin.y, bbMax.z-bbMin.z}) * 0.5 / (eps * 1.25)));
    for (double m : {0.55,0.7,0.9,1.15,1.45,1.8,2.3,3.0}) {
        int r2 = min(R, max(3, (int)ceil(rBase * m)));
        int v2 = min(V, max(8, (int)ceil(vBase * m)));
        if (2 + r2 * v2 < N0) dims.push_back({r2,v2});
    }
    for (double q : {16.0,12.0,10.0,8.0,6.0,4.0,3.0}) dims.push_back({max(3,(int)ceil(R/q)), max(8,(int)ceil(V/q))});
    sort(dims.begin(), dims.end(), [](auto a, auto b){ return a.first*a.second < b.first*b.second; });
    dims.erase(unique(dims.begin(), dims.end()), dims.end());
    vector<Vec3> X; vector<Face> Y;
    for (auto [r2,v2] : dims) {
        if (elapsed() > 17.4) break;
        if (!makePolarGrid(R,V,r2,v2,X,Y,vn)) continue;
        if (acceptCandidate(X, Y, 0.935, N0 < 70000 ? 160 : 112, 0.0492)) return true;
    }
    return false;
}

static bool tryBoxGrid() {
    if (N0 < 600 || elapsed() > 9.5) return false;
    double Lx = bbMax.x - bbMin.x, Ly = bbMax.y - bbMin.y, Lz = bbMax.z - bbMin.z;
    if (min({Lx,Ly,Lz}) < 1e-7) return false;
    int stride = max(1, N0 / 80000), nearCnt = 0, tot = 0;
    double tol = max(1e-8, 0.0065 * diagLen);
    for (int i = 0; i < N0; i += stride) {
        const Vec3 &p = origP[i];
        double d = min({fabs(p.x-bbMin.x), fabs(p.x-bbMax.x), fabs(p.y-bbMin.y), fabs(p.y-bbMax.y), fabs(p.z-bbMin.z), fabs(p.z-bbMax.z)});
        if (d <= tol) nearCnt++; tot++;
    }
    if (tot < 10 || nearCnt * 1000 < tot * 990) return false;
    double eps = 0.049 * diagLen;
    vector<double> factors = {1.35, 1.18, 1.02, 0.88};
    vector<Vec3> bestX; vector<Face> bestY;
    for (double fac : factors) {
        if (elapsed() > 16.8) break;
        double step = eps * fac;
        int nx = max(1, (int)ceil(Lx / step)), ny = max(1, (int)ceil(Ly / step)), nz = max(1, (int)ceil(Lz / step));
        long long full = 1LL*(nx+1)*(ny+1)*(nz+1);
        if (full > 2000000) continue;
        vector<int> id((size_t)full, -1);
        vector<Vec3> X; vector<Face> Y;
        auto key3 = [&](int i,int j,int k){ return (int)(((long long)i*(ny+1)+j)*(nz+1)+k); };
        auto get = [&](int i,int j,int k)->int{
            int key = key3(i,j,k); int &r = id[key]; if (r >= 0) return r;
            double x = bbMin.x + Lx * (double)i / nx;
            double y = bbMin.y + Ly * (double)j / ny;
            double z = bbMin.z + Lz * (double)k / nz;
            r = (int)X.size(); X.push_back({x,y,z}); return r;
        };
        auto addTri = [&](int a,int b,int c, Vec3 out){ addOriented(Y, X, Face{{a,b,c}}, out); };
        auto addQuad = [&](int a,int b,int c,int d, Vec3 out){ addTri(a,b,c,out); addTri(a,c,d,out); };
        for (int j=0;j<ny;++j) for (int k=0;k<nz;++k) {
            addQuad(get(0,j,k), get(0,j,k+1), get(0,j+1,k+1), get(0,j+1,k), Vec3{-1,0,0});
            addQuad(get(nx,j,k), get(nx,j+1,k), get(nx,j+1,k+1), get(nx,j,k+1), Vec3{1,0,0});
        }
        for (int i=0;i<nx;++i) for (int k=0;k<nz;++k) {
            addQuad(get(i,0,k), get(i+1,0,k), get(i+1,0,k+1), get(i,0,k+1), Vec3{0,-1,0});
            addQuad(get(i,ny,k), get(i,ny,k+1), get(i+1,ny,k+1), get(i+1,ny,k), Vec3{0,1,0});
        }
        for (int i=0;i<nx;++i) for (int j=0;j<ny;++j) {
            addQuad(get(i,j,0), get(i,j+1,0), get(i+1,j+1,0), get(i+1,j,0), Vec3{0,0,-1});
            addQuad(get(i,j,nz), get(i+1,j,nz), get(i+1,j+1,nz), get(i,j+1,nz), Vec3{0,0,1});
        }
        if ((int)X.size() >= N0) continue;
        if (candidateFacesBasicOK(X,Y) && coverageOK(X, 0.0494, 18.2)) {
            double q = 1.0;
            if (elapsed() < 16.8) q = visualScoreCandidate(X,Y,N0 < 80000 ? 160 : 112, 19.0);
            if (q >= 0.94) { installCandidate(X,Y); return true; }
        }
    }
    return false;
}

static void jacobiEigen(double a[3][3], Vec3 axis[3]) {
    double v[3][3] = {{1,0,0},{0,1,0},{0,0,1}};
    for (int it=0; it<45; ++it) {
        int p=0,q=1; double best=fabs(a[0][1]);
        if (fabs(a[0][2]) > best) { p=0; q=2; best=fabs(a[0][2]); }
        if (fabs(a[1][2]) > best) { p=1; q=2; best=fabs(a[1][2]); }
        if (best < 1e-18) break;
        double app=a[p][p], aqq=a[q][q], apq=a[p][q];
        double tau=(aqq-app)/(2*apq);
        double t=(tau>=0?1:-1)/(fabs(tau)+sqrt(1+tau*tau));
        double c=1/sqrt(1+t*t), s=t*c;
        for (int k=0;k<3;++k) if (k!=p && k!=q) { double akp=a[k][p], akq=a[k][q]; a[k][p]=a[p][k]=c*akp-s*akq; a[k][q]=a[q][k]=s*akp+c*akq; }
        a[p][p]=c*c*app-2*s*c*apq+s*s*aqq; a[q][q]=s*s*app+2*s*c*apq+c*c*aqq; a[p][q]=a[q][p]=0;
        for (int k=0;k<3;++k) { double vp=v[k][p], vq=v[k][q]; v[k][p]=c*vp-s*vq; v[k][q]=s*vp+c*vq; }
    }
    int ord[3]={0,1,2}; sort(ord,ord+3,[&](int i,int j){return a[i][i]>a[j][j];});
    for (int j=0;j<3;++j) { int col=ord[j]; axis[j]=unit3(Vec3{v[0][col],v[1][col],v[2][col]}); }
    if (dot3(cross3(axis[0],axis[1]),axis[2]) < 0) axis[2] = axis[2] * -1;
}

struct EllFit { bool ok=false; Vec3 cen, ax[3]; double r[3]; double rms=1e9, maxe=1e9; };
static EllFit fitEllipsoidPCA() {
    EllFit fit; if (N0 < 800) return fit;
    Vec3 mean{}; for (const Vec3 &p: origP) mean = mean + p; mean = mean / (double)N0;
    double cov[3][3]{}; int stride=max(1,N0/250000), sampled=0;
    for (int i=0;i<N0;i+=stride) { Vec3 q=origP[i]-mean; double x[3]={q.x,q.y,q.z}; for(int a=0;a<3;++a) for(int b=0;b<3;++b) cov[a][b]+=x[a]*x[b]; sampled++; }
    for(int a=0;a<3;++a) for(int b=0;b<3;++b) cov[a][b]/=max(1,sampled);
    jacobiEigen(cov, fit.ax);
    double lo[3]={1e100,1e100,1e100}, hi[3]={-1e100,-1e100,-1e100};
    for (const Vec3 &p: origP) for(int k=0;k<3;++k) { double t=dot3(p,fit.ax[k]); lo[k]=min(lo[k],t); hi[k]=max(hi[k],t); }
    fit.cen={0,0,0}; for(int k=0;k<3;++k) { double mid=(lo[k]+hi[k])*0.5; fit.cen=fit.cen+fit.ax[k]*mid; fit.r[k]=(hi[k]-lo[k])*0.5; if(fit.r[k]<1e-12) return fit; }
    double sum=0,mx=0; sampled=0;
    for(int i=0;i<N0;i+=stride) { Vec3 q=origP[i]-fit.cen; double rr=0; for(int k=0;k<3;++k){ double u=dot3(q,fit.ax[k])/fit.r[k]; rr+=u*u; } double e=fabs(sqrt(max(0.0,rr))-1.0); sum+=e*e; mx=max(mx,e); sampled++; }
    fit.rms=sqrt(sum/max(1,sampled)); fit.maxe=mx;
    fit.ok = (fit.rms < (N0 < 5000 ? 0.010 : 0.0075) && fit.maxe < (N0 < 5000 ? 0.050 : 0.035));
    return fit;
}
static bool makeEllipsoid(const EllFit &fit, int lat, int lon, vector<Vec3> &X, vector<Face> &Y) {
    if (!fit.ok || lat < 4 || lon < 8) return false;
    X.clear(); Y.clear(); X.reserve(2 + (lat-1)*lon); Y.reserve(2*lon*(lat-1));
    auto pt = [&](double x,double y,double z){ return fit.cen + fit.ax[0]*(fit.r[0]*x) + fit.ax[1]*(fit.r[1]*y) + fit.ax[2]*(fit.r[2]*z); };
    X.push_back(pt(0,0,1)); X.push_back(pt(0,0,-1));
    auto rid=[&](int r,int j){ return 2+(r-1)*lon+((j%lon+lon)%lon); };
    for(int r=1;r<=lat-1;++r){ double th=M_PI*r/lat, st=sin(th), ct=cos(th); for(int j=0;j<lon;++j){ double ph=2*M_PI*j/lon; X.push_back(pt(st*cos(ph),st*sin(ph),ct)); }}
    auto add=[&](int a,int b,int c){ addOriented(Y,X,Face{{a,b,c}}, (X[a]+X[b]+X[c])/3.0 - fit.cen); };
    for(int j=0;j<lon;++j) add(0,rid(1,j),rid(1,j+1));
    for(int r=1;r<=lat-2;++r) for(int j=0;j<lon;++j){ int a=rid(r,j), b=rid(r+1,j), c=rid(r+1,j+1), d=rid(r,j+1); add(a,b,c); add(a,c,d); }
    for(int j=0;j<lon;++j) add(1,rid(lat-1,j+1),rid(lat-1,j));
    return true;
}
static bool tryEllipsoid() {
    if (elapsed() > 9.5 || N0 < 900) return false;
    EllFit fit = fitEllipsoidPCA(); if (!fit.ok) return false;
    vector<pair<int,int>> dims; double eps=0.049*diagLen, rmax=max({fit.r[0],fit.r[1],fit.r[2]});
    int baseLat=max(5,(int)ceil(M_PI*rmax/(eps*1.20))), baseLon=max(10,(int)ceil(2*M_PI*rmax/(eps*1.20)));
    for(double m:{0.65,0.8,1.0,1.25,1.55,2.0}) dims.push_back({max(4,(int)ceil(baseLat*m)), max(8,(int)ceil(baseLon*m))});
    sort(dims.begin(), dims.end(), [](auto a, auto b){ return 2+(a.first-1)*a.second < 2+(b.first-1)*b.second; });
    dims.erase(unique(dims.begin(), dims.end()), dims.end());
    vector<Vec3>X; vector<Face>Y;
    for(auto [lat,lon]:dims){ if(elapsed()>17.0) break; if(!makeEllipsoid(fit,lat,lon,X,Y)) continue; if((int)X.size()>=N0) continue; if(acceptCandidate(X,Y,0.94,N0<70000?160:112,0.0490)) return true; }
    return false;
}

struct TriKey {
    int a,b,c;
    bool operator==(const TriKey &o) const { return a==o.a && b==o.b && c==o.c; }
};
struct TriKeyHash {
    size_t operator()(const TriKey &t) const noexcept {
        uint64_t x = ((uint64_t)(uint32_t)t.a * 11995408973635179863ULL) ^ ((uint64_t)(uint32_t)t.b * 10150724397891781847ULL) ^ (uint64_t)(uint32_t)t.c;
        return U64Hash{}(x);
    }
};

static bool greedyCover(double r, vector<int> &sel, double timeLimit) {
    PointGrid g; g.build(origP, r);
    vector<unsigned char> cov(N0, 0); sel.clear(); sel.reserve(max(16, N0/32));
    double r2 = r*r;
    int covered=0;
    for(int i=0;i<N0;++i){
        if ((i&32767)==0 && elapsed()>timeLimit) return false;
        if(cov[i]) continue;
        int id=i; sel.push_back(id);
        int X=g.ix(origP[id].x), Y=g.iy(origP[id].y), Z=g.iz(origP[id].z);
        for(int dz=-1;dz<=1;++dz){int z=Z+dz;if(z<0||z>=g.nz)continue;
            for(int dy=-1;dy<=1;++dy){int y=Y+dy;if(y<0||y>=g.ny)continue;
                for(int dx=-1;dx<=1;++dx){int x=X+dx;if(x<0||x>=g.nx)continue;
                    for(int q:g.buck[g.key(x,y,z)]) if(!cov[q] && norm2(origP[q]-origP[id])<=r2){cov[q]=1;covered++;}
                }
            }
        }
    }
    return covered==N0;
}

static bool buildClusterCandidate(const vector<int> &sel, double r, vector<Vec3> &X, vector<Face> &Y, double timeLimit) {
    X.clear(); Y.clear();
    if(sel.empty() || (int)sel.size()>=N0) return false;
    X.reserve(sel.size()); for(int id:sel) X.push_back(origP[id]);
    PointGrid gs; gs.build(X, r);
    vector<int> mp(N0,-1); double r2 = r*r*1.000001;
    for(int i=0;i<N0;++i){ if((i&32767)==0 && elapsed()>timeLimit) return false; int k=gs.nearestPoint(origP[i], r2); if(k<0) return false; mp[i]=k; }
    unordered_set<TriKey, TriKeyHash> seen; seen.reserve((size_t)M0*2/3+16);
    Y.reserve(min(M0, (int)sel.size()*8));
    for(int i=0;i<M0;++i){
        int a=mp[origF[i].v[0]], b=mp[origF[i].v[1]], c=mp[origF[i].v[2]];
        if(a<0||b<0||c<0||a==b||a==c||b==c) continue;
        int s[3]={a,b,c}; sort(s,s+3); TriKey tk{s[0],s[1],s[2]};
        if(seen.insert(tk).second) Y.push_back(Face{{a,b,c}});
    }
    return !Y.empty();
}

static bool tryClusterRemap() {
    if (N0 > 180000 || N0 < 1000 || elapsed() > 11.0) return false;
    vector<double> scales = {0.0490,0.0460,0.0420,0.0380};
    for(double sc:scales){
        if(elapsed()>16.3) break;
        vector<int> sel;
        double r = sc * diagLen;
        if(!greedyCover(r, sel, 15.8)) return false;
        if((int)sel.size() >= activeVertices()*85/100) continue;
        vector<Vec3>X; vector<Face>Y;
        if(!buildClusterCandidate(sel, r*1.02, X, Y, 17.0)) continue;
        if(acceptCandidate(X,Y,0.945,N0<70000?144:96,0.0492)) return true;
    }
    return false;
}

static bool meshManifoldQuick() {
    unordered_map<uint64_t,int,U64Hash> ec;
    ec.reserve((size_t)aliveFaceCount * 3 * 2 + 16);
    double minA2 = 1e-30 * diagLen * diagLen * diagLen * diagLen;
    int facesSeen=0;
    for (int i=0;i<(int)F.size();++i) if(aliveF[i]) {
        Face f=F[i];
        if(f.v[0]<0||f.v[1]<0||f.v[2]<0||f.v[0]>=curN||f.v[1]>=curN||f.v[2]>=curN) return false;
        if(!aliveV[f.v[0]]||!aliveV[f.v[1]]||!aliveV[f.v[2]]) return false;
        if(f.v[0]==f.v[1]||f.v[0]==f.v[2]||f.v[1]==f.v[2]) return false;
        if(norm2(faceCrossCur(f)) <= minA2) return false;
        ec[edgeKey(f.v[0],f.v[1])]++; ec[edgeKey(f.v[1],f.v[2])]++; ec[edgeKey(f.v[2],f.v[0])]++; facesSeen++;
    }
    if(facesSeen==0) return false;
    for(auto &kv:ec) if(kv.second!=2) return false;
    return true;
}

struct MeshState { vector<Vec3>P; vector<Face>F; vector<unsigned char>aliveV,aliveF; vector<vector<int>>adj; vector<double>coverRad; int aliveFaceCount; };
static MeshState saveState(){ return MeshState{P,F,aliveV,aliveF,adj,coverRad,aliveFaceCount}; }
static void loadState(const MeshState&s){ P=s.P; F=s.F; aliveV=s.aliveV; aliveF=s.aliveF; adj=s.adj; coverRad=s.coverRad; aliveFaceCount=s.aliveFaceCount; }

static void collectActiveMesh(vector<Vec3> &X, vector<Face> &Y) {
    vector<int> mp(curN, -1); X.clear(); Y.clear(); X.reserve(activeVertices()); Y.reserve(aliveFaceCount);
    for(int i=0;i<curN;++i) if(aliveV[i]) { mp[i]=(int)X.size(); X.push_back(P[i]); }
    for(int i=0;i<(int)F.size();++i) if(aliveF[i]) {
        Face f=F[i];
        if(f.v[0]>=0&&f.v[1]>=0&&f.v[2]>=0&&f.v[0]<curN&&f.v[1]<curN&&f.v[2]<curN&&mp[f.v[0]]>=0&&mp[f.v[1]]>=0&&mp[f.v[2]]>=0&&f.v[0]!=f.v[1]&&f.v[0]!=f.v[2]&&f.v[1]!=f.v[2]) Y.push_back(Face{{mp[f.v[0]],mp[f.v[1]],mp[f.v[2]]}});
    }
}

static void generalSimplify() {
    rebuildAdj(); markA.assign(curN,0); markB.assign(curN,0); stampA=stampB=1;
    SmoothStats st = analyzeSmooth();
    bool verySmooth = st.samples > 500 && st.smooth >= 0.970 && st.sharp <= 0.035 && st.bad <= 0.025;
    bool smooth = st.samples > 500 && st.coarse >= 0.925 && st.verySharp <= 0.130 && st.bad <= 0.035;
    vector<CollapseParams> passes;
    if (verySmooth) {
        passes.push_back({0.0490,0.042,0.52,1e-25});
        passes.push_back({0.0490,0.034,0.66,1e-25});
        passes.push_back({0.0485,0.027,0.78,1e-25});
        passes.push_back({0.0480,0.020,0.88,1e-25});
    } else if (smooth) {
        passes.push_back({0.0480,0.032,0.70,1e-25});
        passes.push_back({0.0470,0.024,0.84,1e-25});
        passes.push_back({0.0455,0.017,0.92,1e-24});
        passes.push_back({0.0430,0.012,0.965,1e-24});
    } else {
        passes.push_back({0.0440,0.017,0.92,1e-24});
        passes.push_back({0.0415,0.012,0.965,1e-24});
        passes.push_back({0.0385,0.007,0.990,1e-24});
    }
    if (N0 < 20) { passes.clear(); passes.push_back({0.049,0.050,0.40,1e-25}); }
    int target = verySmooth ? max(8, (int)(N0 * 0.050)) : (smooth ? max(8, (int)(N0 * 0.085)) : max(8, (int)(N0 * 0.155)));
    double timeLimit = N0 > 300000 ? 19.15 : 18.35;
    for(size_t i=0;i<passes.size() && elapsed()<timeLimit;++i){ int rounds = (N0 > 300000 ? 1 : 2); edgePass(passes[i], timeLimit, rounds, target); rebuildAdj(); }
}

static void optionalAggressivePost() {
    if (elapsed() > 15.8 || N0 > 180000 || N0 < 4000) return;
    MeshState base = saveState(); int baseV = activeVertices();
    vector<CollapseParams> ps = {{0.0490,0.045,0.62,1e-25},{0.0490,0.038,0.72,1e-25},{0.0485,0.030,0.82,1e-24}};
    MeshState best; int bestV = baseV; bool ok=false;
    for(auto cp:ps){
        if(elapsed()>18.0) break;
        loadState(base);
        edgePass(cp,18.2,2,max(8,N0/28)); rebuildAdj();
        int now=activeVertices(); if(now>=bestV || now>=baseV*98/100) continue;
        if(!meshManifoldQuick()) continue;
        vector<Vec3>X; vector<Face>Y; collectActiveMesh(X,Y);
        double q=visualScoreCandidate(X,Y,N0<70000?128:96,19.15);
        if(q>=0.94){ best=saveState(); bestV=now; ok=true; }
    }
    if(ok) loadState(best); else loadState(base);
}

static void outputMesh() {
    vector<int> mp(curN, -1); int nout=0, fout=0;
    for(int i=0;i<curN;++i) if(aliveV[i]) mp[i]=nout++;
    for(int i=0;i<(int)F.size();++i) if(aliveF[i]) {
        Face f=F[i];
        if(f.v[0]>=0&&f.v[1]>=0&&f.v[2]>=0&&f.v[0]<curN&&f.v[1]<curN&&f.v[2]<curN&&mp[f.v[0]]>=0&&mp[f.v[1]]>=0&&mp[f.v[2]]>=0&&f.v[0]!=f.v[1]&&f.v[0]!=f.v[2]&&f.v[1]!=f.v[2]) fout++;
    }
    string out; out.reserve((size_t)nout*42 + (size_t)fout*26 + 64); char line[128];
    out.append(line, snprintf(line,sizeof(line), "%d %d\n", nout, fout));
    for(int i=0;i<curN;++i) if(aliveV[i]) out.append(line, snprintf(line,sizeof(line), "v %.10g %.10g %.10g\n", P[i].x,P[i].y,P[i].z));
    for(int i=0;i<(int)F.size();++i) if(aliveF[i]) {
        Face f=F[i];
        if(f.v[0]>=0&&f.v[1]>=0&&f.v[2]>=0&&f.v[0]<curN&&f.v[1]<curN&&f.v[2]<curN&&mp[f.v[0]]>=0&&mp[f.v[1]]>=0&&mp[f.v[2]]>=0&&f.v[0]!=f.v[1]&&f.v[0]!=f.v[2]&&f.v[1]!=f.v[2]) out.append(line, snprintf(line,sizeof(line), "f %d %d %d\n", mp[f.v[0]]+1, mp[f.v[1]]+1, mp[f.v[2]]+1));
    }
    fwrite(out.data(),1,out.size(),stdout);
}

int main(){
    T0 = chrono::steady_clock::now();
    readInput();
    bool done=false;
    if(!done) done = tryTubeStructured();
    if(!done) done = tryPolarStructured();
    if(!done) done = tryBoxGrid();
    if(!done) done = tryEllipsoid();
    if(!done) done = tryClusterRemap();
    if(!done) { generalSimplify(); optionalAggressivePost(); }
    if(elapsed() < 19.4) {
        bool ok = meshManifoldQuick();
        if(!ok) { curN=N0; P=origP; F=origF; aliveV.assign(curN,1); aliveF.assign(F.size(),1); aliveFaceCount=(int)F.size(); }
    }
    outputMesh();
    return 0;
}
