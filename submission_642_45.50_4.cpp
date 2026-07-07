#include <bits/stdc++.h>
using namespace std;

struct Vec3 {
    double x, y, z;
    Vec3() : x(0), y(0), z(0) {}
    Vec3(double X, double Y, double Z) : x(X), y(Y), z(Z) {}
    Vec3 operator + (const Vec3& o) const { return Vec3(x+o.x, y+o.y, z+o.z); }
    Vec3 operator - (const Vec3& o) const { return Vec3(x-o.x, y-o.y, z-o.z); }
    Vec3 operator * (double s) const { return Vec3(x*s, y*s, z*s); }
    Vec3 operator / (double s) const { return Vec3(x/s, y/s, z/s); }
};
static inline double dot3(const Vec3& a, const Vec3& b) { return a.x*b.x + a.y*b.y + a.z*b.z; }
static inline Vec3 cross3(const Vec3& a, const Vec3& b) {
    return Vec3(a.y*b.z-a.z*b.y, a.z*b.x-a.x*b.z, a.x*b.y-a.y*b.x);
}
static inline double norm2(const Vec3& a) { return dot3(a,a); }
static inline double norm3(const Vec3& a) { return sqrt(norm2(a)); }
static inline Vec3 normalized(const Vec3& a) {
    double n = norm3(a);
    if (n <= 0) return Vec3(0,0,0);
    return a / n;
}

struct Face {
    int a, b, c;
    Vec3 n;
    double area;
    double w;
    unsigned char alive;
};

struct Vertex {
    Vec3 p;
    double q[10];
    vector<int> inc;
    int head, tail, csz;
    int version;
    unsigned char alive;
};

struct Snapshot {
    vector<Vec3> V;
    vector<array<int,3>> F;
    double ratio;
};

static int N0, F0;
static vector<Vec3> origP;
static vector<Vertex> Vv;
static vector<Face> Ff;
static vector<int> nextMember;
static int activeV = 0, activeF = 0;
static double epsH = 0.0, epsH2 = 0.0;
static chrono::steady_clock::time_point timeStart;

static inline double elapsed_sec() {
    using namespace chrono;
    return duration<double>(steady_clock::now() - timeStart).count();
}

// q storage: [xx,xy,xz,xw, yy,yz,yw, zz,zw, ww]
static inline void q_clear(double q[10]) { for (int i=0;i<10;i++) q[i]=0.0; }
static inline void q_add(double q[10], const double r[10]) { for (int i=0;i<10;i++) q[i]+=r[i]; }
static inline void q_plane_add(double q[10], double a, double b, double c, double d, double w) {
    q[0] += w*a*a; q[1] += w*a*b; q[2] += w*a*c; q[3] += w*a*d;
    q[4] += w*b*b; q[5] += w*b*c; q[6] += w*b*d;
    q[7] += w*c*c; q[8] += w*c*d;
    q[9] += w*d*d;
}
static inline double q_eval_arr(const double q[10], const Vec3& p) {
    double x=p.x,y=p.y,z=p.z;
    return q[0]*x*x + 2*q[1]*x*y + 2*q[2]*x*z + 2*q[3]*x
         + q[4]*y*y + 2*q[5]*y*z + 2*q[6]*y
         + q[7]*z*z + 2*q[8]*z + q[9];
}

static inline uint64_t edge_key(int a, int b) {
    if (a > b) swap(a,b);
    return (uint64_t)(uint32_t)a << 32 | (uint32_t)b;
}
static inline int key_u(uint64_t k) { return (int)(k >> 32); }
static inline int key_v(uint64_t k) { return (int)(k & 0xffffffffu); }

static vector<char> slurp_stdin() {
    vector<char> buf;
    buf.reserve(1 << 27);
    char chunk[1 << 16];
    size_t n;
    while ((n = fread(chunk, 1, sizeof(chunk), stdin)) > 0) {
        buf.insert(buf.end(), chunk, chunk + n);
    }
    buf.push_back('\0');
    return buf;
}

static void load_input() {
    vector<char> buf = slurp_stdin();
    char* p = buf.data();
    N0 = (int)strtol(p, &p, 10);
    F0 = (int)strtol(p, &p, 10);
    origP.resize(N0);
    Vv.resize(N0);
    nextMember.assign(N0, -1);
    for (int i=0;i<N0;i++) {
        while (*p==' ' || *p=='\n' || *p=='\r' || *p=='\t') ++p;
        if (*p=='v') ++p;
        double x = strtod(p, &p);
        double y = strtod(p, &p);
        double z = strtod(p, &p);
        origP[i] = Vec3(x,y,z);
        Vv[i].p = origP[i];
        q_clear(Vv[i].q);
        Vv[i].head = Vv[i].tail = i;
        Vv[i].csz = 1;
        Vv[i].version = 0;
        Vv[i].alive = 1;
    }
    Ff.resize(F0);
    vector<int> deg(N0, 0);
    for (int i=0;i<F0;i++) {
        while (*p==' ' || *p=='\n' || *p=='\r' || *p=='\t') ++p;
        if (*p=='f') ++p;
        int a = (int)strtol(p, &p, 10) - 1;
        int b = (int)strtol(p, &p, 10) - 1;
        int c = (int)strtol(p, &p, 10) - 1;
        Ff[i].a=a; Ff[i].b=b; Ff[i].c=c;
        Ff[i].alive = 1;
        deg[a]++; deg[b]++; deg[c]++;
    }
    for (int i=0;i<N0;i++) Vv[i].inc.reserve(deg[i] + 4);
    for (int i=0;i<F0;i++) {
        Vv[Ff[i].a].inc.push_back(i);
        Vv[Ff[i].b].inc.push_back(i);
        Vv[Ff[i].c].inc.push_back(i);
    }
    Vec3 mn(1e100,1e100,1e100), mx(-1e100,-1e100,-1e100);
    for (const Vec3& v: origP) {
        mn.x=min(mn.x,v.x); mn.y=min(mn.y,v.y); mn.z=min(mn.z,v.z);
        mx.x=max(mx.x,v.x); mx.y=max(mx.y,v.y); mx.z=max(mx.z,v.z);
    }
    double diag = norm3(mx-mn);
    epsH = 0.05 * diag * 0.999999; // small numeric safety margin
    epsH2 = epsH * epsH;
    activeV = N0;
    activeF = F0;
}

static inline bool contains_vertex(const Face& f, int v) {
    return f.a==v || f.b==v || f.c==v;
}
static inline bool contains_edge(const Face& f, int u, int v) {
    return contains_vertex(f,u) && contains_vertex(f,v);
}
static inline int third_of_edge(const Face& f, int u, int v) {
    if (f.a!=u && f.a!=v) return f.a;
    if (f.b!=u && f.b!=v) return f.b;
    return f.c;
}

static bool recompute_face(int fid) {
    Face& f = Ff[fid];
    Vec3 a = Vv[f.a].p, b = Vv[f.b].p, c = Vv[f.c].p;
    Vec3 cr = cross3(b-a, c-a);
    double l = norm3(cr);
    if (!(l > 1e-18)) return false;
    f.n = cr / l;
    f.area = 0.5 * l;
    f.w = f.area * (0.25 + fabs(f.n.x) + fabs(f.n.y) + fabs(f.n.z));
    return true;
}

static void initialize_quadrics() {
    for (int i=0;i<F0;i++) {
        recompute_face(i);
        Face& f = Ff[i];
        Vec3 p = Vv[f.a].p;
        double d = -dot3(f.n, p);
        double w = max(f.w, 1e-18);
        q_plane_add(Vv[f.a].q, f.n.x, f.n.y, f.n.z, d, w);
        q_plane_add(Vv[f.b].q, f.n.x, f.n.y, f.n.z, d, w);
        q_plane_add(Vv[f.c].q, f.n.x, f.n.y, f.n.z, d, w);
    }
}

static void clean_incident(int v) {
    if (!Vv[v].alive) return;
    vector<int>& L = Vv[v].inc;
    if (L.empty()) return;
    // For small lists, the stale overhead is lower than rewriting it repeatedly.
    if (L.size() < 96) return;
    size_t w = 0;
    for (size_t i=0;i<L.size();i++) {
        int fid = L[i];
        if (fid >= 0 && fid < (int)Ff.size() && Ff[fid].alive && contains_vertex(Ff[fid], v)) {
            L[w++] = fid;
        }
    }
    L.resize(w);
}

static void get_edge_faces(int u, int v, vector<int>& efaces, vector<int>& opps) {
    efaces.clear(); opps.clear();
    if (!Vv[u].alive || !Vv[v].alive) return;
    // Scan the shorter incident list; stale entries are filtered.
    int s = (Vv[u].inc.size() <= Vv[v].inc.size()) ? u : v;
    for (int fid : Vv[s].inc) {
        if (!Ff[fid].alive) continue;
        const Face& f = Ff[fid];
        if (contains_edge(f,u,v)) {
            efaces.push_back(fid);
            opps.push_back(third_of_edge(f,u,v));
            if (efaces.size() > 2) return;
        }
    }
}

static void get_neighbors(int u, vector<int>& nb) {
    nb.clear();
    if (!Vv[u].alive) return;
    clean_incident(u);
    for (int fid : Vv[u].inc) {
        if (!Ff[fid].alive) continue;
        const Face& f = Ff[fid];
        if (f.a==u) { if (Vv[f.b].alive) nb.push_back(f.b); if (Vv[f.c].alive) nb.push_back(f.c); }
        else if (f.b==u) { if (Vv[f.a].alive) nb.push_back(f.a); if (Vv[f.c].alive) nb.push_back(f.c); }
        else if (f.c==u) { if (Vv[f.a].alive) nb.push_back(f.a); if (Vv[f.b].alive) nb.push_back(f.b); }
    }
    sort(nb.begin(), nb.end());
    nb.erase(unique(nb.begin(), nb.end()), nb.end());
}

static bool link_condition_ok(int u, int v, const vector<int>& opps) {
    static vector<int> nu, nv, inter;
    if (opps.size() != 2 || opps[0] == opps[1]) return false;
    get_neighbors(u, nu);
    get_neighbors(v, nv);
    inter.clear();
    size_t i=0,j=0;
    while (i<nu.size() && j<nv.size()) {
        if (nu[i] == nv[j]) { inter.push_back(nu[i]); ++i; ++j; }
        else if (nu[i] < nv[j]) ++i;
        else ++j;
    }
    if (inter.size() != 2) return false;
    int a=opps[0], b=opps[1]; if (a>b) swap(a,b);
    int c=inter[0], d=inter[1]; if (c>d) swap(c,d);
    return a==c && b==d;
}

static bool solve_qem_position(const double q[10], Vec3& out) {
    double A[3][4] = {
        {q[0], q[1], q[2], -q[3]},
        {q[1], q[4], q[5], -q[6]},
        {q[2], q[5], q[7], -q[8]}
    };
    // Gaussian elimination with partial pivoting.
    for (int col=0; col<3; col++) {
        int piv = col;
        for (int r=col+1; r<3; r++) if (fabs(A[r][col]) > fabs(A[piv][col])) piv = r;
        if (fabs(A[piv][col]) < 1e-14) return false;
        if (piv != col) for (int k=col;k<4;k++) swap(A[piv][k], A[col][k]);
        double div = A[col][col];
        for (int k=col;k<4;k++) A[col][k] /= div;
        for (int r=0;r<3;r++) if (r != col) {
            double m = A[r][col];
            if (m == 0) continue;
            for (int k=col;k<4;k++) A[r][k] -= m*A[col][k];
        }
    }
    out = Vec3(A[0][3], A[1][3], A[2][3]);
    return isfinite(out.x) && isfinite(out.y) && isfinite(out.z);
}

static bool coverage_ok(int u, int v, const Vec3& cand) {
    for (int h = Vv[u].head; h != -1; h = nextMember[h]) {
        if (norm2(origP[h] - cand) > epsH2) return false;
    }
    for (int h = Vv[v].head; h != -1; h = nextMember[h]) {
        if (norm2(origP[h] - cand) > epsH2) return false;
    }
    return true;
}

struct EdgeEval {
    double cost;
    Vec3 cand;
};

static bool evaluate_local_change(int u, int v, const Vec3& cand, double& normalPenalty) {
    static vector<int> local;
    static vector<array<int,3>> keys;
    local.clear(); keys.clear();
    for (int fid : Vv[u].inc) if (Ff[fid].alive) local.push_back(fid);
    for (int fid : Vv[v].inc) if (Ff[fid].alive) local.push_back(fid);
    sort(local.begin(), local.end());
    local.erase(unique(local.begin(), local.end()), local.end());
    normalPenalty = 0.0;
    keys.reserve(local.size());
    for (int fid : local) {
        const Face& f = Ff[fid];
        bool hu = contains_vertex(f,u), hv = contains_vertex(f,v);
        if (!hu && !hv) continue;
        if (hu && hv) continue; // the two edge-adjacent faces disappear
        int ia=f.a, ib=f.b, ic=f.c;
        Vec3 pa = Vv[ia].p, pb = Vv[ib].p, pc = Vv[ic].p;
        if (ia==u || ia==v) pa = cand;
        if (ib==u || ib==v) pb = cand;
        if (ic==u || ic==v) pc = cand;
        // The topological test should prevent equal indices. Keep this guard for safety.
        int ka=ia, kb=ib, kc=ic;
        if (ka==u || ka==v) ka = -1;
        if (kb==u || kb==v) kb = -1;
        if (kc==u || kc==v) kc = -1;
        // Use a local synthetic id for the merged vertex while detecting duplicate faces.
        const int MID = INT_MAX;
        ka = (ia==u || ia==v) ? MID : ia;
        kb = (ib==u || ib==v) ? MID : ib;
        kc = (ic==u || ic==v) ? MID : ic;
        if (ka==kb || kb==kc || ka==kc) return false;
        array<int,3> kk = {ka,kb,kc};
        sort(kk.begin(), kk.end());
        keys.push_back(kk);
        Vec3 cr = cross3(pb-pa, pc-pa);
        double l2 = norm2(cr);
        if (!(l2 > 1e-28)) return false;
        double l = sqrt(l2);
        Vec3 nn = cr / l;
        double d = dot3(nn, f.n);
        if (d < -0.10) return false; // reject actual local inversions, allow large simplification bends
        double oneMinus = max(0.0, 1.0 - d);
        normalPenalty += f.w * oneMinus * oneMinus;
    }
    sort(keys.begin(), keys.end());
    for (size_t i=1;i<keys.size();i++) if (keys[i] == keys[i-1]) return false;
    return true;
}

static bool compute_edge_eval(int u, int v, bool fullTopo, EdgeEval& ev) {
    if (u == v || !Vv[u].alive || !Vv[v].alive) return false;
    static vector<int> efaces, opps;
    get_edge_faces(u, v, efaces, opps);
    if (efaces.size() != 2) return false;
    if (fullTopo) {
        if (activeV <= 4) return false;
        if (!link_condition_ok(u, v, opps)) return false;
    }
    double qsum[10];
    for (int i=0;i<10;i++) qsum[i] = Vv[u].q[i] + Vv[v].q[i];

    Vec3 pu = Vv[u].p, pv = Vv[v].p;
    Vec3 mid = (pu + pv) * 0.5;
    Vec3 weighted = (pu * (double)Vv[u].csz + pv * (double)Vv[v].csz) / (double)(Vv[u].csz + Vv[v].csz);
    Vec3 opt;
    bool hasOpt = solve_qem_position(qsum, opt);

    Vec3 cands[6];
    int cc = 0;
    if (hasOpt) cands[cc++] = opt;
    cands[cc++] = mid;
    cands[cc++] = weighted;
    cands[cc++] = pu;
    cands[cc++] = pv;
    // A very conservative fallback candidate: average of original cluster endpoints' current reps.
    cands[cc++] = (mid + weighted) * 0.5;

    const Face& f0 = Ff[efaces[0]];
    const Face& f1 = Ff[efaces[1]];
    double crease = max(0.0, 1.0 - dot3(f0.n, f1.n));
    double creaseMul = 1.0 + 1.5 * crease * crease;
    double len2 = norm2(pu - pv);

    bool ok = false;
    double best = 1e300;
    Vec3 bestc = mid;
    for (int i=0;i<cc;i++) {
        Vec3 cand = cands[i];
        if (!isfinite(cand.x) || !isfinite(cand.y) || !isfinite(cand.z)) continue;
        if (!coverage_ok(u, v, cand)) continue;
        double npen = 0.0;
        if (!evaluate_local_change(u, v, cand, npen)) continue;
        double qerr = max(0.0, q_eval_arr(qsum, cand));
        // The normal term is deliberately modest: QEM preserves large-scale shape while the
        // local test prevents true flips. This ordering reaches much lower vertex counts.
        double cost = qerr * creaseMul + 0.035 * npen + 1e-13 * len2;
        // Prefer candidates near the current edge when their QEM is comparable; this reduces shrinkage.
        Vec3 em = mid;
        cost += 1e-7 * norm2(cand - em);
        if (cost < best) {
            best = cost;
            bestc = cand;
            ok = true;
        }
    }
    if (!ok) return false;
    ev.cost = best;
    ev.cand = bestc;
    return true;
}

struct PQEdge {
    double cost;
    int u, v;
    int vu, vv;
    bool operator < (const PQEdge& o) const { return cost > o.cost; }
};

static void push_edge(priority_queue<PQEdge>& pq, int u, int v) {
    if (u==v || !Vv[u].alive || !Vv[v].alive) return;
    EdgeEval ev;
    if (!compute_edge_eval(u, v, false, ev)) return;
    pq.push(PQEdge{ev.cost, u, v, Vv[u].version, Vv[v].version});
}

static Snapshot make_snapshot() {
    Snapshot s;
    s.ratio = (double)activeV / (double)N0;
    vector<int> mp(N0, -1);
    s.V.reserve(activeV);
    for (int i=0;i<N0;i++) if (Vv[i].alive) {
        mp[i] = (int)s.V.size();
        s.V.push_back(Vv[i].p);
    }
    s.F.reserve(activeF);
    for (int i=0;i<F0;i++) if (Ff[i].alive) {
        int a=mp[Ff[i].a], b=mp[Ff[i].b], c=mp[Ff[i].c];
        if (a<0 || b<0 || c<0 || a==b || b==c || a==c) continue;
        s.F.push_back({a,b,c});
    }
    return s;
}

static int apply_collapse(int u, int v, const Vec3& cand) {
    // Keep the endpoint with the larger incident list to limit total update cost.
    int keep = u, lose = v;
    if (Vv[keep].inc.size() < Vv[lose].inc.size()) swap(keep, lose);
    vector<int> loseList;
    loseList.reserve(Vv[lose].inc.size());
    for (int fid : Vv[lose].inc) {
        if (Ff[fid].alive && contains_vertex(Ff[fid], lose)) loseList.push_back(fid);
    }

    Vv[keep].p = cand;
    for (int i=0;i<10;i++) Vv[keep].q[i] += Vv[lose].q[i];
    // O(1) exact cluster merge through a linked list of original vertices.
    if (Vv[keep].tail != -1) nextMember[Vv[keep].tail] = Vv[lose].head;
    else Vv[keep].head = Vv[lose].head;
    Vv[keep].tail = Vv[lose].tail;
    Vv[keep].csz += Vv[lose].csz;

    vector<int> addToKeep;
    addToKeep.reserve(loseList.size());
    for (int fid : loseList) {
        if (!Ff[fid].alive) continue;
        Face& f = Ff[fid];
        if (contains_vertex(f, keep)) {
            f.alive = 0;
            activeF--;
            continue;
        }
        if (f.a == lose) f.a = keep;
        if (f.b == lose) f.b = keep;
        if (f.c == lose) f.c = keep;
        // This was validated before, but recompute current normal for future local checks/rendering.
        if (!recompute_face(fid)) {
            f.alive = 0;
            activeF--;
            continue;
        }
        addToKeep.push_back(fid);
    }
    Vv[keep].inc.insert(Vv[keep].inc.end(), addToKeep.begin(), addToKeep.end());
    vector<int>().swap(Vv[lose].inc);
    Vv[lose].alive = 0;
    Vv[lose].version++;
    Vv[keep].version++;
    activeV--;
    if (Vv[keep].inc.size() > 512) clean_incident(keep);
    return keep;
}

static priority_queue<PQEdge> build_initial_queue() {
    vector<uint64_t> edges;
    edges.reserve((size_t)F0 * 3);
    for (int i=0;i<F0;i++) {
        if (!Ff[i].alive) continue;
        const Face& f = Ff[i];
        edges.push_back(edge_key(f.a, f.b));
        edges.push_back(edge_key(f.b, f.c));
        edges.push_back(edge_key(f.c, f.a));
    }
    sort(edges.begin(), edges.end());
    edges.erase(unique(edges.begin(), edges.end()), edges.end());
    priority_queue<PQEdge> pq;
    for (uint64_t k : edges) {
        int u = key_u(k), v = key_v(k);
        push_edge(pq, u, v);
    }
    vector<uint64_t>().swap(edges);
    return pq;
}


static void simplify_mesh_fast_to(double stopRatio) {
    int targetV = max(4, (int)ceil(stopRatio * N0));
    vector<int> used(N0, 0);
    int pass = 1;
    double stopTime = (N0 > 700000 ? 7.5 : 6.0);
    while (activeV > targetV && elapsed_sec() < stopTime) {
        int collapsedThisPass = 0;
        for (int fid=0; fid<F0 && activeV>targetV; fid++) {
            if ((fid & 16383) == 0 && elapsed_sec() > stopTime) break;
            if (!Ff[fid].alive) continue;
            int a=Ff[fid].a, b=Ff[fid].b, c=Ff[fid].c;
            if (!Vv[a].alive || !Vv[b].alive || !Vv[c].alive) continue;
            pair<double, pair<int,int>> e[3] = {
                {norm2(Vv[a].p - Vv[b].p), {a,b}},
                {norm2(Vv[b].p - Vv[c].p), {b,c}},
                {norm2(Vv[c].p - Vv[a].p), {c,a}}
            };
            sort(e, e+3, [](const auto& x, const auto& y){ return x.first < y.first; });
            for (int k=0;k<3;k++) {
                int u=e[k].second.first, v=e[k].second.second;
                if (u==v || !Vv[u].alive || !Vv[v].alive) continue;
                if (used[u] == pass || used[v] == pass) continue;
                EdgeEval ev;
                if (!compute_edge_eval(u, v, true, ev)) continue;
                int keep = apply_collapse(u, v, ev.cand);
                used[u] = used[v] = used[keep] = pass;
                collapsedThisPass++;
                break;
            }
        }
        if (collapsedThisPass == 0) break;
        pass++;
        if (pass == INT_MAX) { fill(used.begin(), used.end(), 0); pass = 1; }
    }
}

// ------------------------ Low-resolution self-check renderer ------------------------
struct ImagePack {
    int R;
    vector<float> depth, nx, ny, nz;
    vector<unsigned char> fg;
};

static inline void camera_basis(int view, Vec3& eye, Vec3& right, Vec3& up, Vec3& forward) {
    const double D = 2.5;
    if (view == 0) { eye=Vec3( D,0,0); forward=Vec3(-1,0,0); right=Vec3(0,0,-1); up=Vec3(0,1,0); }
    if (view == 1) { eye=Vec3(-D,0,0); forward=Vec3( 1,0,0); right=Vec3(0,0, 1); up=Vec3(0,1,0); }
    if (view == 2) { eye=Vec3(0, D,0); forward=Vec3(0,-1,0); right=Vec3(1,0,0); up=Vec3(0,0,-1); }
    if (view == 3) { eye=Vec3(0,-D,0); forward=Vec3(0, 1,0); right=Vec3(1,0,0); up=Vec3(0,0, 1); }
    if (view == 4) { eye=Vec3(0,0, D); forward=Vec3(0,0,-1); right=Vec3(1,0,0); up=Vec3(0,1,0); }
    if (view == 5) { eye=Vec3(0,0,-D); forward=Vec3(0,0, 1); right=Vec3(-1,0,0); up=Vec3(0,1,0); }
}

static inline double edge2d(double ax,double ay,double bx,double by,double cx,double cy) {
    return (cx-ax)*(by-ay) - (cy-ay)*(bx-ax);
}

static void render_mesh_arrays(const vector<Vec3>& verts, const vector<array<int,3>>& faces, int R, int view, ImagePack& img) {
    img.R = R;
    int P = R*R;
    img.depth.assign(P, 255.0f);
    img.nx.assign(P, 127.5f);
    img.ny.assign(P, 127.5f);
    img.nz.assign(P, 127.5f);
    img.fg.assign(P, 0);
    vector<float> zbuf(P, 1e30f);
    Vec3 eye, right, up, forward;
    camera_basis(view, eye, right, up, forward);
    double fscale = 800.0 * (double)R / 1024.0;
    double c = 0.5 * (double)R;
    for (const auto& tri : faces) {
        const Vec3& A = verts[tri[0]];
        const Vec3& B = verts[tri[1]];
        const Vec3& C = verts[tri[2]];
        Vec3 crw = cross3(B-A, C-A);
        double nl = norm3(crw);
        if (!(nl > 1e-18)) continue;
        Vec3 nw = crw / nl;
        Vec3 pts[3] = {A,B,C};
        double sx[3], sy[3], sz[3];
        bool bad=false;
        for (int i=0;i<3;i++) {
            Vec3 q = pts[i] - eye;
            double xc = dot3(q, right);
            double yc = dot3(q, up);
            double zc = dot3(q, forward);
            if (zc <= 1e-6) { bad=true; break; }
            sz[i] = zc;
            sx[i] = fscale * (xc / zc) + c;
            sy[i] = fscale * (yc / zc) + c;
        }
        if (bad) continue;
        double minx = min(sx[0], min(sx[1], sx[2]));
        double maxx = max(sx[0], max(sx[1], sx[2]));
        double miny = min(sy[0], min(sy[1], sy[2]));
        double maxy = max(sy[0], max(sy[1], sy[2]));
        int x0 = max(0, (int)floor(minx - 0.5));
        int x1 = min(R-1, (int)floor(maxx - 0.5));
        int y0 = max(0, (int)floor(miny - 0.5));
        int y1 = min(R-1, (int)floor(maxy - 0.5));
        if (x0>x1 || y0>y1) continue;
        double area = edge2d(sx[0],sy[0], sx[1],sy[1], sx[2],sy[2]);
        if (fabs(area) < 1e-12) continue;
        double invArea = 1.0 / area;
        float cnx = (float)((nw.x + 1.0) * 127.5);
        float cny = (float)((nw.y + 1.0) * 127.5);
        float cnz = (float)((nw.z + 1.0) * 127.5);
        for (int yy=y0; yy<=y1; yy++) {
            double py = yy + 0.5;
            int row = yy * R;
            for (int xx=x0; xx<=x1; xx++) {
                double px = xx + 0.5;
                double w0 = edge2d(sx[1],sy[1], sx[2],sy[2], px,py) * invArea;
                double w1 = edge2d(sx[2],sy[2], sx[0],sy[0], px,py) * invArea;
                double w2 = 1.0 - w0 - w1;
                if (w0 >= -1e-9 && w1 >= -1e-9 && w2 >= -1e-9) {
                    double iz = w0 / sz[0] + w1 / sz[1] + w2 / sz[2];
                    if (iz <= 0) continue;
                    float zz = (float)(1.0 / iz);
                    int id = row + xx;
                    if (zz < zbuf[id]) {
                        zbuf[id] = zz;
                        img.depth[id] = zz;
                        img.nx[id] = cnx; img.ny[id] = cny; img.nz[id] = cnz;
                        img.fg[id] = 1;
                    }
                }
            }
        }
    }
}

static void render_current_original(int R, int view, ImagePack& img) {
    img.R = R;
    int P = R*R;
    img.depth.assign(P, 255.0f);
    img.nx.assign(P, 127.5f);
    img.ny.assign(P, 127.5f);
    img.nz.assign(P, 127.5f);
    img.fg.assign(P, 0);
    vector<float> zbuf(P, 1e30f);
    Vec3 eye, right, up, forward;
    camera_basis(view, eye, right, up, forward);
    double fscale = 800.0 * (double)R / 1024.0;
    double c = 0.5 * (double)R;
    for (int fi=0; fi<F0; fi++) {
        const Face& ff = Ff[fi];
        const Vec3& A = origP[ff.a];
        const Vec3& B = origP[ff.b];
        const Vec3& C = origP[ff.c];
        Vec3 crw = cross3(B-A, C-A);
        double nl = norm3(crw);
        if (!(nl > 1e-18)) continue;
        Vec3 nw = crw / nl;
        Vec3 pts[3] = {A,B,C};
        double sx[3], sy[3], sz[3];
        bool bad=false;
        for (int i=0;i<3;i++) {
            Vec3 q = pts[i] - eye;
            double xc = dot3(q, right);
            double yc = dot3(q, up);
            double zc = dot3(q, forward);
            if (zc <= 1e-6) { bad=true; break; }
            sz[i] = zc;
            sx[i] = fscale * (xc / zc) + c;
            sy[i] = fscale * (yc / zc) + c;
        }
        if (bad) continue;
        double minx = min(sx[0], min(sx[1], sx[2]));
        double maxx = max(sx[0], max(sx[1], sx[2]));
        double miny = min(sy[0], min(sy[1], sy[2]));
        double maxy = max(sy[0], max(sy[1], sy[2]));
        int x0 = max(0, (int)floor(minx - 0.5));
        int x1 = min(R-1, (int)floor(maxx - 0.5));
        int y0 = max(0, (int)floor(miny - 0.5));
        int y1 = min(R-1, (int)floor(maxy - 0.5));
        if (x0>x1 || y0>y1) continue;
        double area = edge2d(sx[0],sy[0], sx[1],sy[1], sx[2],sy[2]);
        if (fabs(area) < 1e-12) continue;
        double invArea = 1.0 / area;
        float cnx = (float)((nw.x + 1.0) * 127.5);
        float cny = (float)((nw.y + 1.0) * 127.5);
        float cnz = (float)((nw.z + 1.0) * 127.5);
        for (int yy=y0; yy<=y1; yy++) {
            double py = yy + 0.5;
            int row = yy * R;
            for (int xx=x0; xx<=x1; xx++) {
                double px = xx + 0.5;
                double w0 = edge2d(sx[1],sy[1], sx[2],sy[2], px,py) * invArea;
                double w1 = edge2d(sx[2],sy[2], sx[0],sy[0], px,py) * invArea;
                double w2 = 1.0 - w0 - w1;
                if (w0 >= -1e-9 && w1 >= -1e-9 && w2 >= -1e-9) {
                    double iz = w0 / sz[0] + w1 / sz[1] + w2 / sz[2];
                    if (iz <= 0) continue;
                    float zz = (float)(1.0 / iz);
                    int id = row + xx;
                    if (zz < zbuf[id]) {
                        zbuf[id] = zz;
                        img.depth[id] = zz;
                        img.nx[id] = cnx; img.ny[id] = cny; img.nz[id] = cnz;
                        img.fg[id] = 1;
                    }
                }
            }
        }
    }
}

static inline double rect_sum(const vector<double>& I, int W, int x0, int y0, int x1, int y1) {
    // inclusive coordinates, image width W-1 for original R; integral stride W
    int A = y0*W + x0;
    int B = y0*W + (x1+1);
    int C = (y1+1)*W + x0;
    int D = (y1+1)*W + (x1+1);
    return I[D] - I[B] - I[C] + I[A];
}

static void build_integral(const vector<float>& X, const vector<float>& Y, vector<double>& SX, vector<double>& SY, vector<double>& SX2, vector<double>& SY2, vector<double>& SXY, int R) {
    int W = R + 1;
    int SZ = W * W;
    SX.assign(SZ, 0.0); SY.assign(SZ, 0.0); SX2.assign(SZ, 0.0); SY2.assign(SZ, 0.0); SXY.assign(SZ, 0.0);
    for (int y=0; y<R; y++) {
        double rx=0, ry=0, rx2=0, ry2=0, rxy=0;
        for (int x=0; x<R; x++) {
            int id = y*R+x;
            double a=X[id], b=Y[id];
            rx += a; ry += b; rx2 += a*a; ry2 += b*b; rxy += a*b;
            int ii = (y+1)*W + (x+1);
            int upi = y*W + (x+1);
            SX[ii] = SX[upi] + rx;
            SY[ii] = SY[upi] + ry;
            SX2[ii] = SX2[upi] + rx2;
            SY2[ii] = SY2[upi] + ry2;
            SXY[ii] = SXY[upi] + rxy;
        }
    }
}

static double channel_ssim(const vector<float>& X, const vector<float>& Y, const vector<unsigned char>& fgX, const vector<unsigned char>& fgY, int R) {
    static vector<double> SX, SY, SX2, SY2, SXY;
    build_integral(X, Y, SX, SY, SX2, SY2, SXY, R);
    int rad = (R <= 128 ? 1 : 2);
    const double c1 = (0.01 * 255.0) * (0.01 * 255.0);
    const double c2 = (0.03 * 255.0) * (0.03 * 255.0);
    int W = R + 1;
    double acc = 0.0;
    int cnt = 0;
    for (int y=0; y<R; y++) for (int x=0; x<R; x++) {
        int id = y*R+x;
        if (!fgX[id] && !fgY[id]) continue;
        int x0=max(0,x-rad), x1=min(R-1,x+rad);
        int y0=max(0,y-rad), y1=min(R-1,y+rad);
        double n = (double)(x1-x0+1)*(y1-y0+1);
        double sx = rect_sum(SX,W,x0,y0,x1,y1);
        double sy = rect_sum(SY,W,x0,y0,x1,y1);
        double sx2 = rect_sum(SX2,W,x0,y0,x1,y1);
        double sy2 = rect_sum(SY2,W,x0,y0,x1,y1);
        double sxy = rect_sum(SXY,W,x0,y0,x1,y1);
        double mux=sx/n, muy=sy/n;
        double vx=max(0.0, sx2/n - mux*mux);
        double vy=max(0.0, sy2/n - muy*muy);
        double cov=sxy/n - mux*muy;
        double den = (mux*mux + muy*muy + c1) * (vx + vy + c2);
        double ssim = 1.0;
        if (den > 1e-30) ssim = ((2*mux*muy + c1) * (2*cov + c2)) / den;
        if (!isfinite(ssim)) ssim = 0.0;
        acc += max(-1.0, min(1.0, ssim));
        cnt++;
    }
    if (cnt == 0) return 1.0;
    return acc / cnt;
}

static double approx_final_ssim(const vector<ImagePack>& origImgs, const Snapshot& s, int R) {
    double total = 0.0;
    ImagePack sim;
    for (int view=0; view<6; view++) {
        render_mesh_arrays(s.V, s.F, R, view, sim);
        const ImagePack& o = origImgs[view];
        double snx = channel_ssim(o.nx, sim.nx, o.fg, sim.fg, R);
        double sny = channel_ssim(o.ny, sim.ny, o.fg, sim.fg, R);
        double snz = channel_ssim(o.nz, sim.nz, o.fg, sim.fg, R);
        double sd  = channel_ssim(o.depth, sim.depth, o.fg, sim.fg, R);
        total += 0.5 * ((snx + sny + snz) / 3.0) + 0.5 * sd;
    }
    return total / 6.0;
}

static void simplify_mesh(vector<Snapshot>& snapshots) {
    priority_queue<PQEdge> pq = build_initial_queue();

    vector<double> ratios;
    if (N0 < 50) ratios = {0.85, 0.75, 0.65, 0.55};
    else if (N0 <= 6000) ratios = {0.40, 0.32, 0.25, 0.22, 0.18, 0.14, 0.10, 0.070, 0.043};
    else if (N0 < 20000) ratios = {0.20, 0.16, 0.125, 0.10, 0.090, 0.073, 0.060, 0.045};
    else if (N0 < 50000) ratios = {0.16, 0.125, 0.10, 0.082, 0.067, 0.052, 0.043};
    else ratios = {0.16, 0.125, 0.100, 0.082, 0.067, 0.052, 0.043};

    vector<int> thresholds;
    thresholds.reserve(ratios.size());
    for (double r: ratios) thresholds.push_back(max(4, (int)ceil(r * N0)));
    int nextSnap = 0;
    int targetV = thresholds.empty() ? max(4, (int)ceil(0.08*N0)) : thresholds.back();
    targetV = max(4, targetV);

    // Leave enough time for optional snapshot scoring and output.
    double decimLimit = (N0 > 300000 ? 17.6 : 18.8);
    long long pops = 0, collapses = 0;
    while (activeV > targetV && !pq.empty()) {
        if ((pops & 8191LL) == 0 && elapsed_sec() > decimLimit) break;
        PQEdge e = pq.top(); pq.pop();
        pops++;
        if (e.u<0 || e.v<0 || e.u>=N0 || e.v>=N0) continue;
        if (!Vv[e.u].alive || !Vv[e.v].alive) continue;
        if (Vv[e.u].version != e.vu || Vv[e.v].version != e.vv) continue;
        EdgeEval ev;
        if (!compute_edge_eval(e.u, e.v, true, ev)) continue;
        int keep = apply_collapse(e.u, e.v, ev.cand);
        collapses++;
        static vector<int> nb;
        get_neighbors(keep, nb);
        for (int w : nb) push_edge(pq, keep, w);

        while (nextSnap < (int)thresholds.size() && activeV <= thresholds[nextSnap]) {
            snapshots.push_back(make_snapshot());
            nextSnap++;
        }
    }
    // If no threshold was reached, or the final mesh is smaller than the last saved one, keep a final copy.
    if (snapshots.empty() || snapshots.back().V.size() != (size_t)activeV) {
        snapshots.push_back(make_snapshot());
    }
}

static const Snapshot& choose_snapshot(const vector<Snapshot>& snapshots, const vector<ImagePack>& origImgs, int R, bool haveRender) {
    int bestIndex = 0;
    // A balanced fallback around the known competitive target. It is safer than blindly emitting 4-5%.
    double fallbackTarget;
    if (N0 <= 6000) fallbackTarget = 0.25;
    else if (N0 < 20000) fallbackTarget = 0.090;
    else if (N0 < 50000) fallbackTarget = 0.083;
    else fallbackTarget = 0.083;
    int fallback = -1;
    for (int i=0;i<(int)snapshots.size();i++) {
        if (snapshots[i].ratio <= fallbackTarget) { fallback = i; break; }
    }
    if (fallback == -1) fallback = (int)snapshots.size() - 1;
    bestIndex = fallback;

    if (!haveRender || snapshots.empty() || elapsed_sec() > 19.4) return snapshots[bestIndex];

    vector<int> idx(snapshots.size());
    iota(idx.begin(), idx.end(), 0);
    sort(idx.begin(), idx.end(), [&](int a, int b){ return snapshots[a].V.size() < snapshots[b].V.size(); });

    double threshold = (R >= 384 ? 0.905 : (R <= 128 ? 0.947 : 0.925));
    int safest = idx.back();
    int chosen = safest;
    for (int id : idx) {
        if (elapsed_sec() > 20.0) break;
        double sc = approx_final_ssim(origImgs, snapshots[id], R);
        if (sc >= threshold) { chosen = id; break; }
    }
    return snapshots[chosen];
}

static void save_snapshot(const Snapshot& s) {
    string out;
    out.reserve((size_t)s.V.size()*42 + (size_t)s.F.size()*28 + 64);
    char line[128];
    out.append(line, snprintf(line, sizeof(line), "%d %d\n", (int)s.V.size(), (int)s.F.size()));
    for (const Vec3& p : s.V) {
        out.append(line, snprintf(line, sizeof(line), "v %.10g %.10g %.10g\n", p.x, p.y, p.z));
    }
    for (const auto& f : s.F) {
        out.append(line, snprintf(line, sizeof(line), "f %d %d %d\n", f[0]+1, f[1]+1, f[2]+1));
    }
    fwrite(out.data(), 1, out.size(), stdout);
}

int main() {
    timeStart = chrono::steady_clock::now();
    load_input();
    initialize_quadrics();

    int R = 0;
    bool doRender = (N0 >= 2000);
    if (doRender) {
        if (N0 <= 20000) R = 512;
        else if (N0 <= 60000) R = 384;
        else if (N0 > 700000) R = 96;
        else if (N0 > 150000) R = 112;
        else R = 128;
    }
    vector<ImagePack> origImgs;
    bool haveRender = false;
    if (doRender) {
        origImgs.resize(6);
        for (int view=0; view<6; view++) {
            if (elapsed_sec() > 4.5 && N0 > 300000) { haveRender = false; break; }
            render_current_original(R, view, origImgs[view]);
            haveRender = true;
        }
    }

    vector<Snapshot> snapshots;
    if (N0 > 250000) simplify_mesh_fast_to(0.16);
    simplify_mesh(snapshots);
    const Snapshot& chosen = choose_snapshot(snapshots, origImgs, R, haveRender);
    save_snapshot(chosen);
    return 0;
}
