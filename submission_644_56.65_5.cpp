#include <bits/stdc++.h>
using namespace std;

/*
  IMC Challenge 2026 - simplifygeometry
  C++17 single-file solution.

  Strategy:
    * Topology-preserving batched edge collapses on the input triangular 2-manifold.
    * Endpoint/segment QEM placement with conservative vertex-Hausdorff certificates.
    * Dihedral/screen-aware penalties and a low-resolution six-view SSIM guard.

  The program is deliberately self-contained: no Eigen dependency is used.
*/

struct FastInput {
    vector<char> buf;
    char *p;
    FastInput() {
        char tmp[1 << 16];
        size_t n;
        while ((n = fread(tmp, 1, sizeof(tmp), stdin)) > 0) buf.insert(buf.end(), tmp, tmp + n);
        buf.push_back('\0');
        p = buf.data();
    }
    inline void skip_ws() {
        while (*p == ' ' || *p == '\n' || *p == '\r' || *p == '\t') ++p;
    }
    inline int next_int() {
        skip_ws();
        int sign = 1;
        if (*p == '-') { sign = -1; ++p; }
        int x = 0;
        while (*p >= '0' && *p <= '9') { x = x * 10 + (*p - '0'); ++p; }
        return sign * x;
    }
    inline double next_double() {
        skip_ws();
        char *q;
        double x = strtod(p, &q);
        p = q;
        return x;
    }
    inline char next_char() {
        skip_ws();
        return *p++;
    }
};

struct Vec3 {
    double x, y, z;
    Vec3() : x(0), y(0), z(0) {}
    Vec3(double X, double Y, double Z) : x(X), y(Y), z(Z) {}
    inline Vec3 operator + (const Vec3& o) const { return Vec3(x + o.x, y + o.y, z + o.z); }
    inline Vec3 operator - (const Vec3& o) const { return Vec3(x - o.x, y - o.y, z - o.z); }
    inline Vec3 operator * (double s) const { return Vec3(x * s, y * s, z * s); }
    inline Vec3 operator / (double s) const { return Vec3(x / s, y / s, z / s); }
    inline Vec3& operator += (const Vec3& o) { x += o.x; y += o.y; z += o.z; return *this; }
};

static inline double dotv(const Vec3& a, const Vec3& b) { return a.x*b.x + a.y*b.y + a.z*b.z; }
static inline Vec3 crossv(const Vec3& a, const Vec3& b) {
    return Vec3(a.y*b.z - a.z*b.y, a.z*b.x - a.x*b.z, a.x*b.y - a.y*b.x);
}
static inline double norm2v(const Vec3& a) { return dotv(a, a); }
static inline double normv(const Vec3& a) { return sqrt(norm2v(a)); }
static inline double distv(const Vec3& a, const Vec3& b) { return normv(a - b); }
static inline Vec3 lerpv(const Vec3& a, const Vec3& b, double t) { return a + (b - a) * t; }

struct Vertex {
    Vec3 p;
    double rad;   // conservative max distance from original vertices in this cluster to p
    double near;  // conservative distance from p to at least one original vertex
};

struct Face {
    int a, b, c;
};

struct Quadric {
    // Symmetric 4x4 matrix, packed upper triangle:
    // 0 xx, 1 xy, 2 xz, 3 xw, 4 yy, 5 yz, 6 yw, 7 zz, 8 zw, 9 ww
    double q[10];
    Quadric() { memset(q, 0, sizeof(q)); }
    inline void clear() { memset(q, 0, sizeof(q)); }
    inline void add_plane(const Vec3& n, double d, double w) {
        const double a = n.x, b = n.y, c = n.z;
        q[0] += w*a*a; q[1] += w*a*b; q[2] += w*a*c; q[3] += w*a*d;
        q[4] += w*b*b; q[5] += w*b*c; q[6] += w*b*d;
        q[7] += w*c*c; q[8] += w*c*d; q[9] += w*d*d;
    }
    inline void operator += (const Quadric& o) {
        for (int i = 0; i < 10; ++i) q[i] += o.q[i];
    }
    inline double eval(const Vec3& p) const {
        const double x = p.x, y = p.y, z = p.z;
        return q[0]*x*x + 2*q[1]*x*y + 2*q[2]*x*z + 2*q[3]*x
             + q[4]*y*y + 2*q[5]*y*z + 2*q[6]*y
             + q[7]*z*z + 2*q[8]*z + q[9];
    }
    inline double line_A(const Vec3& d) const {
        const double x = d.x, y = d.y, z = d.z;
        return q[0]*x*x + 2*q[1]*x*y + 2*q[2]*x*z
             + q[4]*y*y + 2*q[5]*y*z + q[7]*z*z;
    }
    inline double line_B(const Vec3& p, const Vec3& d) const {
        // derivative's linear coefficient for eval(p + t*d): A*t^2 + B*t + C
        const double x = p.x, y = p.y, z = p.z;
        const double dx = d.x, dy = d.y, dz = d.z;
        double v = 0.0;
        v += dx * (q[0]*x + q[1]*y + q[2]*z + q[3]);
        v += dy * (q[1]*x + q[4]*y + q[5]*z + q[6]);
        v += dz * (q[2]*x + q[5]*y + q[7]*z + q[8]);
        return 2.0 * v;
    }
};

struct EdgeRec {
    int u, v, opp, fid;
};

struct Candidate {
    float cost;
    int u, v;       // endpoints in current compact mesh; v will be removed into u
    float t;        // new position = p[u] + t * (p[v] - p[u])
    float newRad;
    float newNear;
};

static inline bool edge_less(const EdgeRec& a, const EdgeRec& b) {
    if (a.u != b.u) return a.u < b.u;
    if (a.v != b.v) return a.v < b.v;
    return a.opp < b.opp;
}

static inline uint64_t edge_key(int a, int b) {
    if (a > b) swap(a, b);
    return (uint64_t)(uint32_t)a << 32 | (uint32_t)b;
}

class Solver {
public:
    vector<Vertex> V;
    vector<Face> F;
    int originalN = 0;
    double bboxDiag = 0.0;
    double hausTol = 0.0;
    chrono::steady_clock::time_point t0;

    // Low-resolution SSIM guard data.
    bool guardReady = false;
    int guardR = 192;

    struct Rendered {
        int R = 0;
        // Six views; each view has 4 scalar images: normal x/y/z in [0,255], depth in [0,255-ish].
        vector<array<vector<float>, 4>> im;
        vector<vector<unsigned char>> fg;
    } originalRender;

    Solver() { t0 = chrono::steady_clock::now(); }

    inline double elapsed() const {
        return chrono::duration<double>(chrono::steady_clock::now() - t0).count();
    }

    void read() {
        FastInput in;
        int n = in.next_int();
        int m = in.next_int();
        originalN = n;
        V.resize(n);
        Vec3 mn(1e100, 1e100, 1e100), mx(-1e100, -1e100, -1e100);
        for (int i = 0; i < n; ++i) {
            (void)in.next_char();
            double x = in.next_double();
            double y = in.next_double();
            double z = in.next_double();
            V[i].p = Vec3(x, y, z);
            V[i].rad = 0.0;
            V[i].near = 0.0;
            mn.x = min(mn.x, x); mn.y = min(mn.y, y); mn.z = min(mn.z, z);
            mx.x = max(mx.x, x); mx.y = max(mx.y, y); mx.z = max(mx.z, z);
        }
        F.resize(m);
        for (int i = 0; i < m; ++i) {
            (void)in.next_char();
            int a = in.next_int() - 1;
            int b = in.next_int() - 1;
            int c = in.next_int() - 1;
            F[i] = {a, b, c};
        }
        bboxDiag = normv(mx - mn);
        hausTol = 0.05 * bboxDiag;
        if (!(hausTol > 0)) hausTol = 1e-9;
    }

    void write() {
        string out;
        out.reserve(1 << 20);
        char line[128];
        auto flush = [&]() {
            if (!out.empty()) { fwrite(out.data(), 1, out.size(), stdout); out.clear(); }
        };
        out.append(line, snprintf(line, sizeof(line), "%d %d\n", (int)V.size(), (int)F.size()));
        for (const auto& v : V) {
            out.append(line, snprintf(line, sizeof(line), "v %.10g %.10g %.10g\n", v.p.x, v.p.y, v.p.z));
            if (out.size() > (1 << 20)) flush();
        }
        for (const auto& f : F) {
            out.append(line, snprintf(line, sizeof(line), "f %d %d %d\n", f.a + 1, f.b + 1, f.c + 1));
            if (out.size() > (1 << 20)) flush();
        }
        flush();
    }

    static void add_edge_rec(vector<EdgeRec>& er, int a, int b, int opp, int fid) {
        if (a > b) swap(a, b);
        er.push_back({a, b, opp, fid});
    }

    bool validate_mesh(const vector<Vertex>& verts, const vector<Face>& faces, bool fullEdgeCheck=true) {
        const int n = (int)verts.size();
        if (n <= 0 || n > originalN) return false;
        vector<uint64_t> edges;
        if (fullEdgeCheck) edges.reserve((size_t)faces.size() * 3);
        const double eps2 = 1e-24;
        for (const Face& f : faces) {
            if (f.a < 0 || f.a >= n || f.b < 0 || f.b >= n || f.c < 0 || f.c >= n) return false;
            if (f.a == f.b || f.b == f.c || f.c == f.a) return false;
            Vec3 cr = crossv(verts[f.b].p - verts[f.a].p, verts[f.c].p - verts[f.a].p);
            if (!(norm2v(cr) > eps2)) return false;
            if (fullEdgeCheck) {
                edges.push_back(edge_key(f.a, f.b));
                edges.push_back(edge_key(f.b, f.c));
                edges.push_back(edge_key(f.c, f.a));
            }
        }
        if (!fullEdgeCheck) return true;
        sort(edges.begin(), edges.end());
        for (size_t i = 0; i < edges.size();) {
            size_t j = i + 1;
            while (j < edges.size() && edges[j] == edges[i]) ++j;
            if (j - i != 2) return false;
            i = j;
        }
        return true;
    }

    bool build_quadrics_and_topology(
            vector<Quadric>& Q,
            vector<Vec3>& faceN,
            vector<int>& faceStart,
            vector<int>& faceInc,
            vector<int>& adjStart,
            vector<int>& adj,
            vector<EdgeRec>& edgeRecs,
            vector<int>& degreeOut) {
        const int n = (int)V.size();
        const int m = (int)F.size();
        Q.assign(n, Quadric());
        faceN.assign(m, Vec3());
        vector<int> fdeg(n, 0);
        edgeRecs.clear();
        edgeRecs.reserve((size_t)m * 3);

        const double areaEps = 1e-30;
        for (int i = 0; i < m; ++i) {
            const Face& f = F[i];
            Vec3 p0 = V[f.a].p, p1 = V[f.b].p, p2 = V[f.c].p;
            Vec3 cr = crossv(p1 - p0, p2 - p0);
            double len = normv(cr);
            if (!(len > areaEps)) return false;
            Vec3 nrm = cr / len;
            faceN[i] = nrm;
            double d = -dotv(nrm, p0);
            // Area weighting is stable across irregular tessellations. A small floor prevents
            // vanishing quadrics on tiny but visually important triangles.
            double w = max(0.5 * len, 1e-12);
            Quadric q;
            q.add_plane(nrm, d, w);
            Q[f.a] += q; Q[f.b] += q; Q[f.c] += q;
            ++fdeg[f.a]; ++fdeg[f.b]; ++fdeg[f.c];
            add_edge_rec(edgeRecs, f.a, f.b, f.c, i);
            add_edge_rec(edgeRecs, f.b, f.c, f.a, i);
            add_edge_rec(edgeRecs, f.c, f.a, f.b, i);
        }

        faceStart.assign(n + 1, 0);
        for (int i = 0; i < n; ++i) faceStart[i+1] = faceStart[i] + fdeg[i];
        faceInc.assign(faceStart[n], 0);
        vector<int> cur = faceStart;
        for (int i = 0; i < m; ++i) {
            const Face& f = F[i];
            faceInc[cur[f.a]++] = i;
            faceInc[cur[f.b]++] = i;
            faceInc[cur[f.c]++] = i;
        }

        sort(edgeRecs.begin(), edgeRecs.end(), edge_less);
        vector<int> deg(n, 0);
        for (size_t i = 0; i < edgeRecs.size();) {
            size_t j = i + 1;
            while (j < edgeRecs.size() && edgeRecs[j].u == edgeRecs[i].u && edgeRecs[j].v == edgeRecs[i].v) ++j;
            if (j - i >= 1) { ++deg[edgeRecs[i].u]; ++deg[edgeRecs[i].v]; }
            i = j;
        }
        adjStart.assign(n + 1, 0);
        for (int i = 0; i < n; ++i) adjStart[i+1] = adjStart[i] + deg[i];
        adj.assign(adjStart[n], 0);
        cur = adjStart;
        for (size_t i = 0; i < edgeRecs.size();) {
            size_t j = i + 1;
            while (j < edgeRecs.size() && edgeRecs[j].u == edgeRecs[i].u && edgeRecs[j].v == edgeRecs[i].v) ++j;
            int u = edgeRecs[i].u, v = edgeRecs[i].v;
            adj[cur[u]++] = v;
            adj[cur[v]++] = u;
            i = j;
        }
        degreeOut.swap(deg);
        return true;
    }

    bool link_condition(int u, int v, int o1, int o2,
                        const vector<int>& adjStart, const vector<int>& adj,
                        vector<int>& mark, int stamp) {
        for (int i = adjStart[u]; i < adjStart[u+1]; ++i) mark[adj[i]] = stamp;
        int cnt = 0;
        bool has1 = false, has2 = false;
        for (int i = adjStart[v]; i < adjStart[v+1]; ++i) {
            int x = adj[i];
            if (mark[x] == stamp) {
                ++cnt;
                if (x == o1) has1 = true;
                if (x == o2) has2 = true;
                if (cnt > 2) return false;
            }
        }
        return cnt == 2 && has1 && has2 && o1 != o2;
    }

    bool feasible_segment_interval(int u, int v, double& lo, double& hi) {
        const double L = distv(V[u].p, V[v].p);
        if (!(L > 1e-15)) return false;
        const double allow = hausTol * 0.985;
        // rad constraints: rad[u] + t L <= allow, rad[v] + (1-t)L <= allow
        double au = (allow - V[u].rad) / L;
        double av = (allow - V[v].rad) / L;
        if (au < -1e-12 || av < -1e-12) return false;
        lo = max(0.0, 1.0 - av);
        hi = min(1.0, au);
        if (lo > hi + 1e-12) return false;
        // simplified vertex must remain near some original vertex in either cluster.
        // It is sufficient if one side certifies it.
        double nu = (allow - V[u].near) / L;
        double nv = (allow - V[v].near) / L;
        double nearLo1 = 0.0, nearHi1 = min(1.0, nu);
        double nearLo2 = max(0.0, 1.0 - nv), nearHi2 = 1.0;
        bool intersects1 = max(lo, nearLo1) <= min(hi, nearHi1) + 1e-12;
        bool intersects2 = max(lo, nearLo2) <= min(hi, nearHi2) + 1e-12;
        if (!intersects1 && !intersects2) return false;
        return true;
    }

    bool local_geometry_ok(int u, int v, const Vec3& np,
                           const vector<Vec3>& faceN,
                           const vector<int>& faceStart,
                           const vector<int>& faceInc) {
        const double eps2 = 1e-22;
        auto face_has_both = [&](const Face& f) {
            bool hu = (f.a == u || f.b == u || f.c == u);
            bool hv = (f.a == v || f.b == v || f.c == v);
            return hu && hv;
        };
        auto check_face = [&](int fid) -> bool {
            const Face& f = F[fid];
            if (face_has_both(f)) return true; // this face disappears after the collapse
            Vec3 p0 = (f.a == u || f.a == v) ? np : V[f.a].p;
            Vec3 p1 = (f.b == u || f.b == v) ? np : V[f.b].p;
            Vec3 p2 = (f.c == u || f.c == v) ? np : V[f.c].p;
            Vec3 cr = crossv(p1 - p0, p2 - p0);
            double a2 = norm2v(cr);
            if (!(a2 > eps2)) return false;
            Vec3 nn = cr / sqrt(a2);
            // Prevent strong local foldovers. Allow moderate changes on high-curvature zones.
            if (dotv(nn, faceN[fid]) < -0.05) return false;
            return true;
        };
        for (int i = faceStart[u]; i < faceStart[u+1]; ++i) if (!check_face(faceInc[i])) return false;
        for (int i = faceStart[v]; i < faceStart[v+1]; ++i) {
            int fid = faceInc[i];
            const Face& f = F[fid];
            if (f.a == u || f.b == u || f.c == u) continue;
            if (!check_face(fid)) return false;
        }
        return true;
    }

    void build_candidates(const vector<Quadric>& Q,
                          const vector<Vec3>& faceN,
                          const vector<int>& adjStart,
                          const vector<int>& adj,
                          const vector<EdgeRec>& edgeRecs,
                          const vector<int>& deg,
                          vector<Candidate>& cand) {
        const int n = (int)V.size();
        cand.clear();
        cand.reserve(edgeRecs.size() / 3);
        vector<int> mark(n, 0);
        int stamp = 1;
        const double tinyLen = 1e-20;

        for (size_t i = 0; i < edgeRecs.size();) {
            size_t j = i + 1;
            while (j < edgeRecs.size() && edgeRecs[j].u == edgeRecs[i].u && edgeRecs[j].v == edgeRecs[i].v) ++j;
            if (j - i == 2) {
                int u0 = edgeRecs[i].u, v0 = edgeRecs[i].v;
                int o1 = edgeRecs[i].opp, o2 = edgeRecs[i+1].opp;
                if (stamp == INT_MAX) { fill(mark.begin(), mark.end(), 0); stamp = 1; }
                if (link_condition(u0, v0, o1, o2, adjStart, adj, mark, stamp++)) {
                    double len = distv(V[u0].p, V[v0].p);
                    if (len > tinyLen) {
                        double lo, hi;
                        if (feasible_segment_interval(u0, v0, lo, hi)) {
                            Quadric qsum = Q[u0];
                            qsum += Q[v0];
                            Vec3 p0 = V[u0].p;
                            Vec3 d = V[v0].p - p0;
                            double A = qsum.line_A(d);
                            double B = qsum.line_B(p0, d);
                            double t = 0.5 * (lo + hi);
                            if (A > 1e-30) t = -B / (2.0 * A);
                            if (t < lo) t = lo;
                            if (t > hi) t = hi;

                            // Try a few feasible locations; this is cheap and helps on planar/crease regions.
                            double bestT = t;
                            Vec3 bestP = lerpv(V[u0].p, V[v0].p, bestT);
                            double bestVal = qsum.eval(bestP);
                            double probes[5] = {lo, hi, 0.5*(lo+hi), t, (fabs(t-0.5*(lo+hi)) > 1e-9 ? t : lo)};
                            for (double tp : probes) {
                                if (tp < lo - 1e-12 || tp > hi + 1e-12) continue;
                                tp = min(hi, max(lo, tp));
                                Vec3 pp = lerpv(V[u0].p, V[v0].p, tp);
                                double val = qsum.eval(pp);
                                if (val < bestVal) { bestVal = val; bestT = tp; bestP = pp; }
                            }

                            double rad = max(V[u0].rad + len * bestT, V[v0].rad + len * (1.0 - bestT));
                            double nearv = min(V[u0].near + len * bestT, V[v0].near + len * (1.0 - bestT));
                            if (rad <= hausTol * 0.986 && nearv <= hausTol * 0.986) {
                                double dihedral = 0.0;
                                if (edgeRecs[i].fid >= 0 && edgeRecs[i+1].fid >= 0) {
                                    double cd = dotv(faceN[edgeRecs[i].fid], faceN[edgeRecs[i+1].fid]);
                                    cd = max(-1.0, min(1.0, cd));
                                    dihedral = 1.0 - fabs(cd);
                                }
                                // Prefer removing the lower-valence endpoint only slightly; the placement is identical.
                                int u = u0, v = v0;
                                double tt = bestT;
                                if (deg[v0] > deg[u0]) { u = v0; v = u0; tt = 1.0 - bestT; }

                                // The base QEM may be exactly zero on planes/creases. The tie-breakers keep
                                // collapses short and postpone sharp visual features until necessary.
                                double cost = bestVal;
                                cost += 1e-12 * len * len;
                                cost *= (1.0 + 75.0 * dihedral * dihedral);
                                cost += 2e-7 * dihedral * len * len;
                                if (rad > hausTol * 0.80) cost *= 1.0 + 4.0 * (rad / max(hausTol, 1e-12) - 0.80);
                                if (!isfinite(cost)) cost = 1e100;
                                Candidate c;
                                c.cost = (float)min(cost, 3.3e38);
                                c.u = u; c.v = v; c.t = (float)tt;
                                c.newRad = (float)rad;
                                c.newNear = (float)nearv;
                                cand.push_back(c);
                            }
                        }
                    }
                }
            }
            i = j;
        }
        sort(cand.begin(), cand.end(), [](const Candidate& a, const Candidate& b) {
            if (a.cost != b.cost) return a.cost < b.cost;
            if (a.u != b.u) return a.u < b.u;
            return a.v < b.v;
        });
    }

    bool apply_pass(bool conservativeLock, double targetRatioThisPass) {
        const int n = (int)V.size();
        if (n <= 4) return false;
        vector<Quadric> Q;
        vector<Vec3> faceN;
        vector<int> faceStart, faceInc, adjStart, adj, deg;
        vector<EdgeRec> edgeRecs;
        if (!build_quadrics_and_topology(Q, faceN, faceStart, faceInc, adjStart, adj, edgeRecs, deg)) return false;

        vector<Candidate> cand;
        build_candidates(Q, faceN, adjStart, adj, edgeRecs, deg, cand);
        if (cand.empty()) return false;

        int wantN = max(4, (int)ceil(originalN * targetRatioThisPass));
        int maxCollapse = max(1, n - wantN);
        // Avoid over-large batches; smaller batches make QEM/SSIM guard more reliable.
        int batchCap = conservativeLock ? max(1, n / 7) : max(1, n / 3);
        maxCollapse = min(maxCollapse, batchCap);

        vector<unsigned char> locked(n, 0), removed(n, 0), touched(n, 0);
        vector<int> to(n, -1);
        vector<Vec3> newP(n);
        vector<double> newRad(n, 0.0), newNear(n, 0.0);
        int selected = 0;

        auto lock_one_ring = [&](int x) {
            locked[x] = 1;
            for (int i = adjStart[x]; i < adjStart[x+1]; ++i) locked[adj[i]] = 1;
        };

        for (const Candidate& c : cand) {
            if (selected >= maxCollapse) break;
            int u = c.u, v = c.v;
            if ((unsigned)u >= (unsigned)n || (unsigned)v >= (unsigned)n) continue;
            if (removed[u] || removed[v] || locked[u] || locked[v]) continue;
            Vec3 np = lerpv(V[u].p, V[v].p, c.t);
            if (!local_geometry_ok(u, v, np, faceN, faceStart, faceInc)) continue;
            removed[v] = 1;
            to[v] = u;
            touched[u] = 1;
            newP[u] = np;
            newRad[u] = c.newRad;
            newNear[u] = c.newNear;
            if (conservativeLock) {
                lock_one_ring(u);
                lock_one_ring(v);
            } else {
                locked[u] = locked[v] = 1;
            }
            ++selected;
        }
        if (selected == 0) return false;

        vector<int> root(n);
        for (int i = 0; i < n; ++i) root[i] = removed[i] ? to[i] : i;
        vector<int> nid(n, -1);
        vector<Vertex> nV;
        nV.reserve(n - selected);
        for (int i = 0; i < n; ++i) if (!removed[i]) {
            nid[i] = (int)nV.size();
            Vertex vv = V[i];
            if (touched[i]) { vv.p = newP[i]; vv.rad = newRad[i]; vv.near = newNear[i]; }
            nV.push_back(vv);
        }
        vector<Face> nF;
        nF.reserve(F.size() - 2 * selected);
        for (const Face& f : F) {
            int a = root[f.a], b = root[f.b], c = root[f.c];
            if (a == b || b == c || c == a) continue;
            a = nid[a]; b = nid[b]; c = nid[c];
            if (a == b || b == c || c == a) continue;
            nF.push_back({a, b, c});
        }

        if (!validate_mesh(nV, nF, true)) {
            if (!conservativeLock) {
                return apply_pass(true, targetRatioThisPass);
            }
            return false;
        }
        V.swap(nV);
        F.swap(nF);
        return true;
    }

    // --------------------- approximate six-view renderer/SSIM guard ---------------------

    struct ProjP { float sx, sy, z; };

    inline void project_point(const Vec3& p, int view, int R, ProjP& out) const {
        const float D = 2.5f;
        float cx = 0, cy = 0, cz = 1;
        switch (view) {
            case 0: // +X camera
                cx = (float)p.y; cy = (float)p.z; cz = D - (float)p.x; break;
            case 1: // -X
                cx = -(float)p.y; cy = (float)p.z; cz = D + (float)p.x; break;
            case 2: // +Y
                cx = -(float)p.x; cy = (float)p.z; cz = D - (float)p.y; break;
            case 3: // -Y
                cx = (float)p.x; cy = (float)p.z; cz = D + (float)p.y; break;
            case 4: // +Z
                cx = (float)p.x; cy = (float)p.y; cz = D - (float)p.z; break;
            default: // -Z
                cx = -(float)p.x; cy = (float)p.y; cz = D + (float)p.z; break;
        }
        float f = 800.0f * (float)R / 1024.0f;
        out.sx = f * cx / cz + 0.5f * R;
        out.sy = f * cy / cz + 0.5f * R;
        out.z = cz;
    }

    inline Vec3 normal_to_view(const Vec3& n, int view) const {
        switch (view) {
            case 0: return Vec3(n.y, n.z, -n.x);
            case 1: return Vec3(-n.y, n.z, n.x);
            case 2: return Vec3(-n.x, n.z, -n.y);
            case 3: return Vec3(n.x, n.z, n.y);
            case 4: return Vec3(n.x, n.y, -n.z);
            default: return Vec3(-n.x, n.y, n.z);
        }
    }

    Rendered render_mesh(const vector<Vertex>& verts, const vector<Face>& faces, int R) const {
        Rendered rr;
        rr.R = R;
        rr.im.resize(6);
        rr.fg.resize(6);
        const int P = R * R;
        for (int view = 0; view < 6; ++view) {
            for (int k = 0; k < 4; ++k) rr.im[view][k].assign(P, k == 3 ? 255.0f : 127.5f);
            rr.fg[view].assign(P, 0);
            vector<float> zbuf(P, 1e30f);
            for (const Face& f : faces) {
                const Vec3& p0w = verts[f.a].p;
                const Vec3& p1w = verts[f.b].p;
                const Vec3& p2w = verts[f.c].p;
                Vec3 crw = crossv(p1w - p0w, p2w - p0w);
                double clen = normv(crw);
                if (!(clen > 1e-30)) continue;
                Vec3 nw = crw / clen;
                Vec3 nv = normal_to_view(nw, view);
                double nl = normv(nv);
                if (nl > 1e-30) nv = nv / nl;

                ProjP p0, p1, p2;
                project_point(p0w, view, R, p0);
                project_point(p1w, view, R, p1);
                project_point(p2w, view, R, p2);
                if (p0.z <= 1e-6f || p1.z <= 1e-6f || p2.z <= 1e-6f) continue;

                float minx = min(p0.sx, min(p1.sx, p2.sx));
                float maxx = max(p0.sx, max(p1.sx, p2.sx));
                float miny = min(p0.sy, min(p1.sy, p2.sy));
                float maxy = max(p0.sy, max(p1.sy, p2.sy));
                int x0 = max(0, (int)floor(minx));
                int x1 = min(R - 1, (int)floor(maxx));
                int y0 = max(0, (int)floor(miny));
                int y1 = min(R - 1, (int)floor(maxy));
                if (x0 > x1 || y0 > y1) continue;

                float ax = p0.sx, ay = p0.sy;
                float bx = p1.sx, by = p1.sy;
                float cx = p2.sx, cy = p2.sy;
                float den = (by - cy) * (ax - cx) + (cx - bx) * (ay - cy);
                if (fabs(den) < 1e-12f) continue;
                float invDen = 1.0f / den;
                float invz0 = 1.0f / p0.z, invz1 = 1.0f / p1.z, invz2 = 1.0f / p2.z;
                float r0 = (float)((nv.x + 1.0) * 127.5);
                float r1 = (float)((nv.y + 1.0) * 127.5);
                float r2 = (float)((nv.z + 1.0) * 127.5);
                for (int y = y0; y <= y1; ++y) {
                    float py = y + 0.5f;
                    int base = y * R;
                    for (int x = x0; x <= x1; ++x) {
                        float px = x + 0.5f;
                        float w0 = ((by - cy) * (px - cx) + (cx - bx) * (py - cy)) * invDen;
                        float w1 = ((cy - ay) * (px - cx) + (ax - cx) * (py - cy)) * invDen;
                        float w2 = 1.0f - w0 - w1;
                        const float eps = -1e-5f;
                        if (w0 >= eps && w1 >= eps && w2 >= eps) {
                            float iz = w0 * invz0 + w1 * invz1 + w2 * invz2;
                            if (iz <= 0.0f) continue;
                            float zz = 1.0f / iz;
                            int id = base + x;
                            if (zz < zbuf[id]) {
                                zbuf[id] = zz;
                                rr.fg[view][id] = 1;
                                rr.im[view][0][id] = r0;
                                rr.im[view][1][id] = r1;
                                rr.im[view][2][id] = r2;
                                rr.im[view][3][id] = zz;
                            }
                        }
                    }
                }
            }
        }
        return rr;
    }

    static double ssim_channel(const vector<float>& A, const vector<float>& B,
                               const vector<unsigned char>& fgA, const vector<unsigned char>& fgB,
                               int R) {
        const int W = R + 1;
        const int N = R * R;
        vector<double> IA((size_t)W * W, 0), IB((size_t)W * W, 0), IA2((size_t)W * W, 0), IB2((size_t)W * W, 0), IAB((size_t)W * W, 0);
        for (int y = 0; y < R; ++y) {
            double sa = 0, sb = 0, sa2 = 0, sb2 = 0, sab = 0;
            for (int x = 0; x < R; ++x) {
                int id = y * R + x;
                double a = A[id], b = B[id];
                sa += a; sb += b; sa2 += a*a; sb2 += b*b; sab += a*b;
                size_t ii = (size_t)(y + 1) * W + (x + 1);
                size_t up = (size_t)y * W + (x + 1);
                IA[ii] = IA[up] + sa;
                IB[ii] = IB[up] + sb;
                IA2[ii] = IA2[up] + sa2;
                IB2[ii] = IB2[up] + sb2;
                IAB[ii] = IAB[up] + sab;
            }
        }
        auto rect = [&](const vector<double>& I, int x0, int y0, int x1, int y1) -> double {
            return I[(size_t)y1 * W + x1] - I[(size_t)y0 * W + x1] - I[(size_t)y1 * W + x0] + I[(size_t)y0 * W + x0];
        };
        const double c1 = (0.01 * 255.0) * (0.01 * 255.0);
        const double c2 = (0.03 * 255.0) * (0.03 * 255.0);
        double sum = 0.0;
        int cnt = 0;
        for (int y = 0; y < R; ++y) {
            for (int x = 0; x < R; ++x) {
                int id = y * R + x;
                if (!fgA[id] && !fgB[id]) continue;
                int x0 = max(0, x - 5), y0 = max(0, y - 5);
                int x1 = min(R, x + 6), y1 = min(R, y + 6);
                double n = (double)(x1 - x0) * (y1 - y0);
                double ma = rect(IA, x0, y0, x1, y1) / n;
                double mb = rect(IB, x0, y0, x1, y1) / n;
                double va = rect(IA2, x0, y0, x1, y1) / n - ma * ma;
                double vb = rect(IB2, x0, y0, x1, y1) / n - mb * mb;
                double co = rect(IAB, x0, y0, x1, y1) / n - ma * mb;
                if (va < 0 && va > -1e-9) va = 0;
                if (vb < 0 && vb > -1e-9) vb = 0;
                double num = (2 * ma * mb + c1) * (2 * co + c2);
                double den = (ma * ma + mb * mb + c1) * (va + vb + c2);
                double s = (den != 0.0 ? num / den : 1.0);
                if (s < -1.0) s = -1.0;
                if (s > 1.0) s = 1.0;
                sum += s;
                ++cnt;
            }
        }
        return cnt ? sum / cnt : 1.0;
    }

    static double final_ssim_approx(const Rendered& A, const Rendered& B) {
        int R = A.R;
        double total = 0.0;
        for (int view = 0; view < 6; ++view) {
            double ns = 0.0;
            for (int k = 0; k < 3; ++k) ns += ssim_channel(A.im[view][k], B.im[view][k], A.fg[view], B.fg[view], R);
            ns /= 3.0;
            double ds = ssim_channel(A.im[view][3], B.im[view][3], A.fg[view], B.fg[view], R);
            total += 0.5 * ns + 0.5 * ds;
        }
        return total / 6.0;
    }

    void prepare_guard() {
        // The guard is most useful on large/visual cases. On tiny meshes it is cheap anyway.
        if (originalN <= 80000) guardR = 1024;
        else if (originalN <= 250000) guardR = 512;
        else if (originalN <= 700000) guardR = 384;
        else guardR = 288;
        if (elapsed() < 3.0 || originalN <= 200000) {
            originalRender = render_mesh(V, F, guardR);
            guardReady = true;
        } else {
            // For very large inputs, still try the smaller guard if there is enough time left.
            originalRender = render_mesh(V, F, guardR);
            guardReady = true;
        }
    }

    double guard_score_current() {
        if (!guardReady) return 1.0;
        Rendered cur = render_mesh(V, F, guardR);
        return final_ssim_approx(originalRender, cur);
    }

    double target_ratio() const {
        // Hard target: the SSIM guard normally stops earlier when the visual threshold is approached.
        if (originalN <= 20) return 0.80;         // sample/cube-like tiny case
        if (originalN <= 6000) return 0.055;
        if (originalN <= 30000) return 0.050;
        if (originalN <= 70000) return 0.047;
        if (originalN <= 450000) return 0.045;
        return 0.050;                            // keep time/visual safety on the largest case
    }

    double guard_threshold() const {
        // Low-resolution SSIM can be optimistic on very fine details and pessimistic on pixel-level aliases.
        // These margins aim to stay safely above the official 0.9 while still allowing strong compression.
        if (originalN <= 20) return 0.999;        // sample should stay visually exact after planar cleanup
        if (originalN <= 6000) return 0.906;
        if (originalN <= 30000) return 0.908;
        if (originalN <= 70000) return 0.912;
        if (originalN <= 450000) return 0.910;
        return 0.912;
    }

    void simplify() {
        if (originalN <= 0) return;
        prepare_guard();

        vector<Vertex> safeV = V;
        vector<Face> safeF = F;
        double bestScore = 1.0;
        const double hardTarget = target_ratio();
        const double ssimNeed = guard_threshold();
        int stagnant = 0;

        for (int pass = 0; pass < 40; ++pass) {
            double ratio = (double)V.size() / max(1, originalN);
            if (ratio <= hardTarget) break;
            if (elapsed() > 18.2) break;

            // Use large steps while the image score has wide slack, then take finer
            // steps near the visual threshold to avoid leaving easy vertices on the table.
            double slack = bestScore - ssimNeed;
            double shrink;
            if (pass == 0) shrink = 0.60;
            else if (slack < 0.006) shrink = 0.95;
            else if (slack < 0.015) shrink = 0.90;
            else if (slack < 0.035) shrink = 0.84;
            else shrink = 0.72;
            double passTarget = max(hardTarget, ratio * shrink);
            bool conservative = (pass >= 2 && (pass % 3 == 2));
            size_t before = V.size();
            bool ok = apply_pass(conservative, passTarget);
            if (!ok) {
                if (!conservative) ok = apply_pass(true, passTarget);
            }
            if (!ok) break;
            if (V.size() >= before) {
                if (++stagnant >= 2) break;
            } else stagnant = 0;

            bool doGuard = guardReady;
            // On very large meshes, guard every other early pass to save time, then every pass near the threshold.
            if (originalN > 250000 && pass < 3 && (pass & 1)) doGuard = false;
            if (doGuard && elapsed() < 19.2) {
                double sc = guard_score_current();
                if (sc + 1e-9 >= ssimNeed) {
                    bestScore = sc;
                    safeV = V;
                    safeF = F;
                } else {
                    // If the first failed guard is close, retain current only when the compression gain is large
                    // and the approximate score is still above a risky-but-plausible official margin.
                    if (sc >= 0.904 && (double)V.size() / max(1, originalN) < (double)safeV.size() / max(1, originalN) * 0.82) {
                        bestScore = sc;
                        safeV = V;
                        safeF = F;
                    }
                    break;
                }
            } else {
                if (!guardReady) {
                    safeV = V;
                    safeF = F;
                }
            }
        }

        // Final safety: never return a mesh that fails our own topological validator.
        if (validate_mesh(safeV, safeF, true)) {
            V.swap(safeV);
            F.swap(safeF);
        }
        (void)bestScore;
    }
};

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    Solver solver;
    solver.read();
    solver.simplify();
    solver.write();
    return 0;
}
