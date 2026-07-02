#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <numeric>
#include <string>
#include <vector>

using namespace std;

struct Vec3 {
    double x, y, z;
};

static inline Vec3 operator+(const Vec3& a, const Vec3& b) { return {a.x + b.x, a.y + b.y, a.z + b.z}; }
static inline Vec3 operator-(const Vec3& a, const Vec3& b) { return {a.x - b.x, a.y - b.y, a.z - b.z}; }
static inline Vec3 operator*(const Vec3& a, double s) { return {a.x * s, a.y * s, a.z * s}; }
static inline Vec3 operator/(const Vec3& a, double s) { return {a.x / s, a.y / s, a.z / s}; }
static inline double dot3(const Vec3& a, const Vec3& b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
static inline Vec3 cross3(const Vec3& a, const Vec3& b) {
    return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
}
static inline double norm2(const Vec3& a) { return dot3(a, a); }
static inline double norm3(const Vec3& a) { return sqrt(norm2(a)); }
static inline Vec3 normalize3(const Vec3& a) {
    double n = norm3(a);
    if (n <= 0.0) return {0.0, 0.0, 0.0};
    return a / n;
}

struct Normal3f {
    float x, y, z;
};

static inline Vec3 to_vec3(const Normal3f& n) { return {(double)n.x, (double)n.y, (double)n.z}; }

struct Face {
    int a, b, c;
};

struct Quadric {
    // xx xy xz xw yy yz yw zz zw ww
    double q[10];
    Quadric() { memset(q, 0, sizeof(q)); }
    void add_plane(const Vec3& n, double d, double w) {
        const double a = n.x, b = n.y, c = n.z;
        q[0] += w * a * a;
        q[1] += w * a * b;
        q[2] += w * a * c;
        q[3] += w * a * d;
        q[4] += w * b * b;
        q[5] += w * b * c;
        q[6] += w * b * d;
        q[7] += w * c * c;
        q[8] += w * c * d;
        q[9] += w * d * d;
    }
};

static inline double eval_q(const Quadric& q, const Vec3& p) {
    return q.q[0] * p.x * p.x + 2.0 * q.q[1] * p.x * p.y + 2.0 * q.q[2] * p.x * p.z +
           2.0 * q.q[3] * p.x + q.q[4] * p.y * p.y + 2.0 * q.q[5] * p.y * p.z +
           2.0 * q.q[6] * p.y + q.q[7] * p.z * p.z + 2.0 * q.q[8] * p.z + q.q[9];
}

static inline double eval_sum_q(const Quadric& a, const Quadric& b, const Vec3& p) {
    return eval_q(a, p) + eval_q(b, p);
}

static vector<Vec3> V;
static vector<Face> F;
static vector<double> cluster_radius;
static vector<int> cluster_size;
static double haus_eps = 0.0;

struct FastInput {
    vector<char> buf;
    char* p;

    FastInput() : p(nullptr) {
        char chunk[1 << 16];
        size_t n;
        while ((n = fread(chunk, 1, sizeof(chunk), stdin)) > 0) {
            buf.insert(buf.end(), chunk, chunk + n);
        }
        buf.push_back('\0');
        p = buf.data();
    }

    inline void skip_ws() {
        while (*p == ' ' || *p == '\n' || *p == '\r' || *p == '\t') ++p;
    }

    int next_int() {
        skip_ws();
        int sign = 1;
        if (*p == '-') {
            sign = -1;
            ++p;
        }
        int x = 0;
        while (*p >= '0' && *p <= '9') {
            x = x * 10 + (*p - '0');
            ++p;
        }
        return sign * x;
    }

    double next_double() {
        skip_ws();
        char* e = nullptr;
        double x = strtod(p, &e);
        p = e;
        return x;
    }

    char next_char() {
        skip_ws();
        return *p++;
    }
};

static void load_mesh() {
    FastInput in;
    int n = in.next_int();
    int m = in.next_int();
    V.resize(n);
    F.resize(m);
    cluster_radius.assign(n, 0.0);
    cluster_size.assign(n, 1);

    for (int i = 0; i < n; ++i) {
        (void)in.next_char();
        V[i].x = in.next_double();
        V[i].y = in.next_double();
        V[i].z = in.next_double();
    }
    for (int i = 0; i < m; ++i) {
        (void)in.next_char();
        F[i].a = in.next_int() - 1;
        F[i].b = in.next_int() - 1;
        F[i].c = in.next_int() - 1;
    }
}

static void save_mesh() {
    string out;
    out.reserve(1 << 22);
    auto flush = [&]() {
        if (!out.empty()) {
            fwrite(out.data(), 1, out.size(), stdout);
            out.clear();
        }
    };
    char line[128];
    out.append(line, snprintf(line, sizeof(line), "%d %d\n", (int)V.size(), (int)F.size()));
    for (int i = 0; i < (int)V.size(); ++i) {
        out.append(line, snprintf(line, sizeof(line), "v %.9g %.9g %.9g\n", V[i].x, V[i].y, V[i].z));
        if (out.size() > (1 << 21)) flush();
    }
    for (int i = 0; i < (int)F.size(); ++i) {
        out.append(line, snprintf(line, sizeof(line), "f %d %d %d\n", F[i].a + 1, F[i].b + 1, F[i].c + 1));
        if (out.size() > (1 << 21)) flush();
    }
    flush();
}

struct EdgeHalf {
    unsigned long long key;
    int opp;
    int fid;
    bool operator<(const EdgeHalf& other) const {
        if (key != other.key) return key < other.key;
        return fid < other.fid;
    }
};

struct Edge {
    int a, b;
    int o1, o2;
    int f1, f2;
};

struct BuildData {
    bool valid = true;
    vector<Edge> edges;
    vector<int> voff, neigh;
    vector<int> foff, incident;
    vector<Quadric> quad;
    vector<Normal3f> fnormal;
    vector<float> farea;
    vector<float> curvature;
    double feature_fraction = 0.0;
};

static inline unsigned long long edge_key(int u, int v) {
    if (u > v) swap(u, v);
    return (unsigned long long)(unsigned int)u << 32 | (unsigned int)v;
}

static BuildData build_data() {
    const int n = (int)V.size();
    const int m = (int)F.size();
    BuildData d;
    d.quad.assign(n, Quadric());
    d.fnormal.resize(m);
    d.farea.resize(m);
    d.curvature.assign(n, 0.0f);

    vector<Vec3> nsum(n, {0.0, 0.0, 0.0});
    vector<double> area_sum(n, 0.0);
    vector<int> fdeg(n, 0);
    vector<EdgeHalf> halves;
    halves.reserve((size_t)m * 3);

    for (int i = 0; i < m; ++i) {
        const Face& f = F[i];
        if (f.a < 0 || f.a >= n || f.b < 0 || f.b >= n || f.c < 0 || f.c >= n ||
            f.a == f.b || f.b == f.c || f.c == f.a) {
            d.valid = false;
            return d;
        }
        Vec3 p0 = V[f.a], p1 = V[f.b], p2 = V[f.c];
        Vec3 cr = cross3(p1 - p0, p2 - p0);
        double area2 = norm3(cr);
        if (area2 <= 1e-18) {
            d.valid = false;
            return d;
        }
        Vec3 nrm = cr / area2;
        d.fnormal[i] = {(float)nrm.x, (float)nrm.y, (float)nrm.z};
        d.farea[i] = (float)area2;
        double plane_d = -dot3(nrm, p0);
        double w = area2;
        d.quad[f.a].add_plane(nrm, plane_d, w);
        d.quad[f.b].add_plane(nrm, plane_d, w);
        d.quad[f.c].add_plane(nrm, plane_d, w);

        nsum[f.a] = nsum[f.a] + nrm * area2;
        nsum[f.b] = nsum[f.b] + nrm * area2;
        nsum[f.c] = nsum[f.c] + nrm * area2;
        area_sum[f.a] += area2;
        area_sum[f.b] += area2;
        area_sum[f.c] += area2;
        ++fdeg[f.a];
        ++fdeg[f.b];
        ++fdeg[f.c];

        halves.push_back({edge_key(f.a, f.b), f.c, i});
        halves.push_back({edge_key(f.b, f.c), f.a, i});
        halves.push_back({edge_key(f.c, f.a), f.b, i});
    }

    for (int i = 0; i < n; ++i) {
        if (area_sum[i] > 0.0) nsum[i] = normalize3(nsum[i]);
    }
    for (int i = 0; i < m; ++i) {
        const Face& f = F[i];
        Vec3 nrm = to_vec3(d.fnormal[i]);
        double ar = d.farea[i];
        d.curvature[f.a] += (float)(ar * max(0.0, 1.0 - dot3(nsum[f.a], nrm)));
        d.curvature[f.b] += (float)(ar * max(0.0, 1.0 - dot3(nsum[f.b], nrm)));
        d.curvature[f.c] += (float)(ar * max(0.0, 1.0 - dot3(nsum[f.c], nrm)));
    }
    for (int i = 0; i < n; ++i) {
        if (area_sum[i] > 0.0) d.curvature[i] = (float)(d.curvature[i] / area_sum[i]);
    }

    d.foff.assign(n + 1, 0);
    for (int i = 0; i < n; ++i) d.foff[i + 1] = d.foff[i] + fdeg[i];
    d.incident.resize((size_t)m * 3);
    vector<int> fcur = d.foff;
    for (int i = 0; i < m; ++i) {
        const Face& f = F[i];
        d.incident[fcur[f.a]++] = i;
        d.incident[fcur[f.b]++] = i;
        d.incident[fcur[f.c]++] = i;
    }

    sort(halves.begin(), halves.end());
    vector<int> deg(n, 0);
    d.edges.reserve(halves.size() / 2);
    long long feature_edges = 0;
    for (size_t i = 0; i < halves.size();) {
        size_t j = i + 1;
        while (j < halves.size() && halves[j].key == halves[i].key) ++j;
        if (j - i != 2) {
            d.valid = false;
            return d;
        }
        int a = (int)(halves[i].key >> 32);
        int b = (int)(halves[i].key & 0xffffffffu);
        Edge e{a, b, halves[i].opp, halves[i + 1].opp, halves[i].fid, halves[i + 1].fid};
        d.edges.push_back(e);
        ++deg[a];
        ++deg[b];

        Vec3 n1 = to_vec3(d.fnormal[e.f1]);
        Vec3 n2 = to_vec3(d.fnormal[e.f2]);
        double co = max(-1.0, min(1.0, dot3(n1, n2)));
        bool silhouette = (n1.x * n2.x < -1e-7) || (n1.y * n2.y < -1e-7) || (n1.z * n2.z < -1e-7);
        if (1.0 - co > 0.018 || silhouette) ++feature_edges;
        i = j;
    }
    if (!d.edges.empty()) d.feature_fraction = (double)feature_edges / (double)d.edges.size();

    d.voff.assign(n + 1, 0);
    for (int i = 0; i < n; ++i) d.voff[i + 1] = d.voff[i] + deg[i];
    d.neigh.resize((size_t)d.voff[n]);
    vector<int> cur = d.voff;
    for (const Edge& e : d.edges) {
        d.neigh[cur[e.a]++] = e.b;
        d.neigh[cur[e.b]++] = e.a;
    }
    return d;
}

static inline double det3(double a00, double a01, double a02,
                          double a10, double a11, double a12,
                          double a20, double a21, double a22) {
    return a00 * (a11 * a22 - a12 * a21) -
           a01 * (a10 * a22 - a12 * a20) +
           a02 * (a10 * a21 - a11 * a20);
}

static bool qem_optimum(const Quadric& qa, const Quadric& qb, Vec3& out) {
    double q[10];
    for (int i = 0; i < 10; ++i) q[i] = qa.q[i] + qb.q[i];
    double a00 = q[0], a01 = q[1], a02 = q[2];
    double a10 = q[1], a11 = q[4], a12 = q[5];
    double a20 = q[2], a21 = q[5], a22 = q[7];
    double b0 = -q[3], b1 = -q[6], b2 = -q[8];
    double d = det3(a00, a01, a02, a10, a11, a12, a20, a21, a22);
    if (fabs(d) < 1e-16) return false;
    double dx = det3(b0, a01, a02, b1, a11, a12, b2, a21, a22);
    double dy = det3(a00, b0, a02, a10, b1, a12, a20, b2, a22);
    double dz = det3(a00, a01, b0, a10, a11, b1, a20, a21, b2);
    out = {dx / d, dy / d, dz / d};
    return isfinite(out.x) && isfinite(out.y) && isfinite(out.z);
}

struct BestPoint {
    Vec3 p;
    double radius;
    double qcost;
};

static bool choose_best_point(int a, int b, const BuildData& data, BestPoint& best) {
    const Vec3 pa = V[a], pb = V[b];
    const Quadric& qa = data.quad[a];
    const Quadric& qb = data.quad[b];
    double best_cost = numeric_limits<double>::infinity();
    bool ok = false;

    auto try_point = [&](const Vec3& p) {
        if (!isfinite(p.x) || !isfinite(p.y) || !isfinite(p.z)) return;
        double r = max(cluster_radius[a] + norm3(p - pa), cluster_radius[b] + norm3(p - pb));
        if (r > haus_eps) return;
        double c = eval_sum_q(qa, qb, p);
        if (c < best_cost) {
            best_cost = c;
            best = {p, r, c};
            ok = true;
        }
    };

    try_point(pa);
    try_point(pb);
    try_point((pa + pb) * 0.5);
    double wa = (double)cluster_size[a];
    double wb = (double)cluster_size[b];
    try_point((pa * wa + pb * wb) / (wa + wb));
    Vec3 opt;
    if (qem_optimum(qa, qb, opt)) try_point(opt);
    return ok;
}

static bool link_condition(const Edge& e, const BuildData& data, vector<int>& mark, int& stamp) {
    if (++stamp == numeric_limits<int>::max()) {
        fill(mark.begin(), mark.end(), 0);
        stamp = 1;
    }
    for (int it = data.voff[e.a]; it < data.voff[e.a + 1]; ++it) mark[data.neigh[it]] = stamp;
    int common = 0;
    bool seen1 = false, seen2 = false;
    for (int it = data.voff[e.b]; it < data.voff[e.b + 1]; ++it) {
        int x = data.neigh[it];
        if (mark[x] == stamp) {
            ++common;
            if (x == e.o1) seen1 = true;
            if (x == e.o2) seen2 = true;
        }
    }
    return common == 2 && seen1 && seen2 && e.o1 != e.o2;
}

static bool face_has_vertex(const Face& f, int v) {
    return f.a == v || f.b == v || f.c == v;
}

static bool local_geometry_ok(int a, int b, const Vec3& p, const BuildData& data) {
    auto check_face = [&](int fid) -> bool {
        const Face& f = F[fid];
        bool ha = face_has_vertex(f, a);
        bool hb = face_has_vertex(f, b);
        if (ha && hb) return true; // The two incident edge faces disappear.
        if (!ha && !hb) return true;
        Vec3 p0 = (f.a == a || f.a == b) ? p : V[f.a];
        Vec3 p1 = (f.b == a || f.b == b) ? p : V[f.b];
        Vec3 p2 = (f.c == a || f.c == b) ? p : V[f.c];
        Vec3 cr = cross3(p1 - p0, p2 - p0);
        double ar = norm3(cr);
        if (ar <= 1e-16) return false;
        Vec3 oldn = to_vec3(data.fnormal[fid]);
        if (dot3(cr, oldn) <= ar * 1e-5) return false;
        return true;
    };

    for (int it = data.foff[a]; it < data.foff[a + 1]; ++it) {
        if (!check_face(data.incident[it])) return false;
    }
    for (int it = data.foff[b]; it < data.foff[b + 1]; ++it) {
        if (!check_face(data.incident[it])) return false;
    }
    return true;
}

struct Candidate {
    float cost;
    int edge_id;
    bool operator<(const Candidate& other) const { return cost < other.cost; }
};

struct Collapse {
    int keep, remove;
    Vec3 p;
    double radius;
};

static double collapse_cost(const Edge& e, const BuildData& data, const BestPoint& bp) {
    Vec3 pa = V[e.a], pb = V[e.b];
    double len2 = norm2(pa - pb);
    Vec3 n1 = to_vec3(data.fnormal[e.f1]);
    Vec3 n2 = to_vec3(data.fnormal[e.f2]);
    double co = max(-1.0, min(1.0, dot3(n1, n2)));
    double crease = max(0.0, 1.0 - co);
    double silhouette = 0.0;
    if (n1.x * n2.x < -1e-7) silhouette += 1.0;
    if (n1.y * n2.y < -1e-7) silhouette += 1.0;
    if (n1.z * n2.z < -1e-7) silhouette += 1.0;
    double curv = max((double)data.curvature[e.a], (double)data.curvature[e.b]);
    double feature_penalty = (1e-10 + 1.8 * crease * crease + 0.16 * silhouette + 0.07 * curv) * len2;
    double radius_penalty = 2e-8 * (bp.radius / max(haus_eps, 1e-12)) * (bp.radius / max(haus_eps, 1e-12));
    return bp.qcost + feature_penalty + radius_penalty;
}

static int choose_target_vertices(int original_n, double feature_fraction) {
    if (original_n <= 8) return original_n;
    if (original_n <= 20) return original_n - 1;

    double base;
    if (original_n <= 6000) base = 0.105;
    else if (original_n <= 30000) base = 0.085;
    else if (original_n <= 70000) base = 0.075;
    else if (original_n <= 500000) base = 0.064;
    else base = 0.058;

    double ratio = base + min(0.085, 0.34 * feature_fraction);
    ratio = min(0.165, max(0.055, ratio));
    int target = (int)ceil(original_n * ratio);
    return max(4, min(original_n, target));
}

static void apply_collapses(const vector<Collapse>& collapses) {
    const int n = (int)V.size();
    vector<int> redirect(n);
    iota(redirect.begin(), redirect.end(), 0);
    vector<unsigned char> removed(n, 0);

    for (const Collapse& c : collapses) {
        redirect[c.remove] = c.keep;
        removed[c.remove] = 1;
        V[c.keep] = c.p;
        cluster_radius[c.keep] = c.radius;
        cluster_size[c.keep] += cluster_size[c.remove];
    }

    vector<Face> nf;
    nf.reserve(F.size());
    vector<unsigned char> used(n, 0);
    for (const Face& f : F) {
        int a = redirect[f.a], b = redirect[f.b], c = redirect[f.c];
        if (a == b || b == c || c == a) continue;
        Vec3 cr = cross3(V[b] - V[a], V[c] - V[a]);
        if (norm2(cr) <= 1e-28) continue;
        nf.push_back({a, b, c});
        used[a] = used[b] = used[c] = 1;
    }

    vector<int> id(n, -1);
    int nn = 0;
    for (int i = 0; i < n; ++i) {
        if (!removed[i] && used[i]) id[i] = nn++;
    }

    vector<Vec3> nv(nn);
    vector<double> nr(nn);
    vector<int> ns(nn);
    for (int i = 0; i < n; ++i) {
        if (id[i] >= 0) {
            nv[id[i]] = V[i];
            nr[id[i]] = cluster_radius[i];
            ns[id[i]] = cluster_size[i];
        }
    }
    for (Face& f : nf) {
        f.a = id[f.a];
        f.b = id[f.b];
        f.c = id[f.c];
    }
    V.swap(nv);
    F.swap(nf);
    cluster_radius.swap(nr);
    cluster_size.swap(ns);
}

static bool light_validate() {
    const int n = (int)V.size();
    if (n <= 0 || F.empty()) return false;
    vector<unsigned long long> keys;
    keys.reserve((size_t)F.size() * 3);
    for (const Face& f : F) {
        if (f.a < 0 || f.a >= n || f.b < 0 || f.b >= n || f.c < 0 || f.c >= n) return false;
        if (f.a == f.b || f.b == f.c || f.c == f.a) return false;
        if (norm2(cross3(V[f.b] - V[f.a], V[f.c] - V[f.a])) <= 1e-28) return false;
        keys.push_back(edge_key(f.a, f.b));
        keys.push_back(edge_key(f.b, f.c));
        keys.push_back(edge_key(f.c, f.a));
    }
    sort(keys.begin(), keys.end());
    for (size_t i = 0; i < keys.size();) {
        size_t j = i + 1;
        while (j < keys.size() && keys[j] == keys[i]) ++j;
        if (j - i != 2) return false;
        i = j;
    }
    return true;
}

static void simplify() {
    const int original_n = (int)V.size();
    if (original_n <= 1) return;

    double xmin = V[0].x, xmax = V[0].x;
    double ymin = V[0].y, ymax = V[0].y;
    double zmin = V[0].z, zmax = V[0].z;
    for (const Vec3& p : V) {
        xmin = min(xmin, p.x); xmax = max(xmax, p.x);
        ymin = min(ymin, p.y); ymax = max(ymax, p.y);
        zmin = min(zmin, p.z); zmax = max(zmax, p.z);
    }
    double diag = sqrt((xmax - xmin) * (xmax - xmin) +
                       (ymax - ymin) * (ymax - ymin) +
                       (zmax - zmin) * (zmax - zmin));
    haus_eps = max(1e-12, diag * 0.05 * 0.992);

    auto start = chrono::steady_clock::now();
    auto elapsed = [&]() {
        return chrono::duration<double>(chrono::steady_clock::now() - start).count();
    };

    int target = -1;
    vector<Vec3> prevV;
    vector<Face> prevF;
    vector<double> prevR;
    vector<int> prevS;
    bool have_prev = false;

    for (int pass = 0; pass < 30; ++pass) {
        if (elapsed() > 18.7) break;
        BuildData data = build_data();
        if (!data.valid) {
            if (have_prev) {
                V.swap(prevV);
                F.swap(prevF);
                cluster_radius.swap(prevR);
                cluster_size.swap(prevS);
            }
            break;
        }
        if (target < 0) target = choose_target_vertices(original_n, data.feature_fraction);
        if ((int)V.size() <= target) break;

        vector<Candidate> cand;
        cand.reserve(data.edges.size());
        for (int i = 0; i < (int)data.edges.size(); ++i) {
            BestPoint bp;
            const Edge& e = data.edges[i];
            if (!choose_best_point(e.a, e.b, data, bp)) continue;
            double c = collapse_cost(e, data, bp);
            if (!isfinite(c)) continue;
            cand.push_back({(float)c, i});
        }
        sort(cand.begin(), cand.end());

        vector<int> mark(V.size(), 0);
        int stamp = 0;
        vector<unsigned char> blocked(V.size(), 0);
        vector<Collapse> chosen;
        int need_remove = (int)V.size() - target;
        chosen.reserve(min(need_remove, (int)V.size() / 3 + 1));

        for (const Candidate& ca : cand) {
            if ((int)chosen.size() >= need_remove) break;
            const Edge& e = data.edges[ca.edge_id];
            if (blocked[e.a] || blocked[e.b]) continue;
            if (!link_condition(e, data, mark, stamp)) continue;
            BestPoint bp;
            if (!choose_best_point(e.a, e.b, data, bp)) continue;
            if (!local_geometry_ok(e.a, e.b, bp.p, data)) continue;

            int keep = e.a, rem = e.b;
            if (cluster_size[keep] < cluster_size[rem]) swap(keep, rem);
            blocked[e.a] = blocked[e.b] = 1;
            for (int it = data.voff[e.a]; it < data.voff[e.a + 1]; ++it) blocked[data.neigh[it]] = 1;
            for (int it = data.voff[e.b]; it < data.voff[e.b + 1]; ++it) blocked[data.neigh[it]] = 1;
            chosen.push_back({keep, rem, bp.p, bp.radius});
        }

        if (chosen.empty()) break;

        prevV = V;
        prevF = F;
        prevR = cluster_radius;
        prevS = cluster_size;
        have_prev = true;

        apply_collapses(chosen);
    }

    if (!light_validate() && have_prev) {
        V.swap(prevV);
        F.swap(prevF);
        cluster_radius.swap(prevR);
        cluster_size.swap(prevS);
    }
}

int main() {
    load_mesh();
    simplify();
    save_mesh();
    return 0;
}
