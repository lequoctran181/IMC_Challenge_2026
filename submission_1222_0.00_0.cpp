#include <bits/stdc++.h>
using namespace std;

struct Vec3 {
    double x, y, z;
    Vec3() : x(0), y(0), z(0) {}
    Vec3(double X, double Y, double Z) : x(X), y(Y), z(Z) {}
    Vec3 operator + (const Vec3& o) const { return Vec3(x + o.x, y + o.y, z + o.z); }
    Vec3 operator - (const Vec3& o) const { return Vec3(x - o.x, y - o.y, z - o.z); }
    Vec3 operator * (double s) const { return Vec3(x * s, y * s, z * s); }
    Vec3& operator += (const Vec3& o) { x += o.x; y += o.y; z += o.z; return *this; }
};

static inline double dot3(const Vec3& a, const Vec3& b) {
    return a.x*b.x + a.y*b.y + a.z*b.z;
}
static inline Vec3 cross3(const Vec3& a, const Vec3& b) {
    return Vec3(a.y*b.z - a.z*b.y, a.z*b.x - a.x*b.z, a.x*b.y - a.y*b.x);
}
static inline double norm2(const Vec3& a) { return dot3(a, a); }
static inline double norm(const Vec3& a) { return sqrt(max(0.0, norm2(a))); }
static inline double dist3(const Vec3& a, const Vec3& b) { return norm(a - b); }
static inline Vec3 normalized(const Vec3& a) {
    double n = norm(a);
    if (n <= 0) return Vec3(0,0,0);
    return a * (1.0 / n);
}

struct Face {
    int v[3];
    unsigned char alive;
};

static inline uint64_t edge_key_int(int a, int b) {
    if (a > b) swap(a, b);
    return (uint64_t)(uint32_t)a << 32 | (uint32_t)b;
}

struct Candidate {
    float cost;
    int a, b;
    bool operator < (const Candidate& o) const { return cost > o.cost; }
};

class Simplifier {
public:
    int N0 = 0, M0 = 0;
    vector<Vec3> p;
    vector<array<int,3>> original_faces;
    vector<Face> faces;
    vector<unsigned char> alive_v;
    vector<vector<int>> inc;
    vector<double> rad;
    vector<double> importance;
    vector<Vec3> initial_face_n;
    vector<int> mark;
    vector<int> mark2;
    int stamp = 1;
    int alive_vertices = 0;
    int alive_faces = 0;
    double diag = 1.0;
    double haus_eps = 0.05;
    double area_eps2 = 1e-30;
    int stage = 0;
    chrono::steady_clock::time_point t0;
    priority_queue<Candidate> pq;

    explicit Simplifier(vector<Vec3> P, vector<array<int,3>> F)
        : N0((int)P.size()), M0((int)F.size()), p(std::move(P)), original_faces(std::move(F)) {
        faces.resize(M0);
        for (int i = 0; i < M0; ++i) {
            faces[i].v[0] = original_faces[i][0];
            faces[i].v[1] = original_faces[i][1];
            faces[i].v[2] = original_faces[i][2];
            faces[i].alive = 1;
        }
        alive_v.assign(N0, 1);
        inc.assign(N0, {});
        rad.assign(N0, 0.0);
        importance.assign(N0, 1.0);
        mark.assign(N0, 0);
        mark2.assign(N0, 0);
        alive_vertices = N0;
        alive_faces = M0;
        compute_bbox();
        area_eps2 = max(1e-36, pow(max(diag, 1e-12), 4.0) * 1e-28);
        haus_eps = 0.0492 * max(diag, 1e-12);
        build_incidence();
        compute_initial_importance();
    }

    void compute_bbox() {
        if (p.empty()) { diag = 1.0; return; }
        Vec3 mn = p[0], mx = p[0];
        for (const auto& q : p) {
            mn.x = min(mn.x, q.x); mn.y = min(mn.y, q.y); mn.z = min(mn.z, q.z);
            mx.x = max(mx.x, q.x); mx.y = max(mx.y, q.y); mx.z = max(mx.z, q.z);
        }
        diag = norm(mx - mn);
        if (!(diag > 0) || !isfinite(diag)) diag = 1.0;
    }

    void build_incidence() {
        for (int i = 0; i < M0; ++i) {
            if (!faces[i].alive) continue;
            for (int j = 0; j < 3; ++j) {
                int v = faces[i].v[j];
                if (0 <= v && v < N0) inc[v].push_back(i);
            }
        }
    }

    bool face_has(int fid, int v) const {
        const Face& f = faces[fid];
        return f.alive && (f.v[0] == v || f.v[1] == v || f.v[2] == v);
    }

    int other_in_face(int fid, int a, int b) const {
        const Face& f = faces[fid];
        for (int i = 0; i < 3; ++i) if (f.v[i] != a && f.v[i] != b) return f.v[i];
        return -1;
    }

    vector<int>& clean_inc(int v) {
        vector<int>& L = inc[v];
        int w = 0;
        for (int fid : L) {
            if (fid >= 0 && fid < (int)faces.size() && face_has(fid, v)) L[w++] = fid;
        }
        L.resize(w);
        return L;
    }

    Vec3 raw_face_normal_from_vertices(int a, int b, int c) const {
        return cross3(p[b] - p[a], p[c] - p[a]);
    }

    Vec3 raw_face_normal(int fid) const {
        const Face& f = faces[fid];
        return raw_face_normal_from_vertices(f.v[0], f.v[1], f.v[2]);
    }

    void compute_initial_importance() {
        initial_face_n.assign(M0, Vec3(0,0,0));
        vector<Vec3> acc(N0, Vec3(0,0,0));
        vector<double> dev(N0, 0.0), sil(N0, 0.0);
        vector<uint64_t> edge_keys;
        vector<int> edge_faces;
        edge_keys.reserve((size_t)M0 * 3);
        edge_faces.reserve((size_t)M0 * 3);

        for (int i = 0; i < M0; ++i) {
            Vec3 rn = raw_face_normal(i);
            Vec3 n = normalized(rn);
            initial_face_n[i] = n;
            double ar = max(1e-300, norm(rn));
            for (int k = 0; k < 3; ++k) acc[faces[i].v[k]] += n * ar;
            for (int e = 0; e < 3; ++e) {
                int a = faces[i].v[e], b = faces[i].v[(e+1)%3];
                edge_keys.push_back(edge_key_int(a, b));
                edge_faces.push_back(i);
            }
        }

        vector<Vec3> vn(N0);
        for (int i = 0; i < N0; ++i) vn[i] = normalized(acc[i]);
        for (int i = 0; i < M0; ++i) {
            Vec3 n = initial_face_n[i];
            for (int k = 0; k < 3; ++k) {
                int v = faces[i].v[k];
                double d = 1.0 - dot3(n, vn[v]);
                if (isfinite(d)) dev[v] = max(dev[v], max(0.0, d));
            }
        }

        vector<int> order(edge_keys.size());
        iota(order.begin(), order.end(), 0);
        sort(order.begin(), order.end(), [&](int A, int B) {
            if (edge_keys[A] != edge_keys[B]) return edge_keys[A] < edge_keys[B];
            return edge_faces[A] < edge_faces[B];
        });
        for (size_t i = 0; i < order.size();) {
            size_t j = i + 1;
            while (j < order.size() && edge_keys[order[j]] == edge_keys[order[i]]) ++j;
            if (j - i == 2) {
                int f1 = edge_faces[order[i]], f2 = edge_faces[order[i+1]];
                Vec3 n1 = initial_face_n[f1], n2 = initial_face_n[f2];
                double dih = max(0.0, 1.0 - dot3(n1, n2));
                uint64_t key = edge_keys[order[i]];
                int a = (int)(key >> 32), b = (int)(key & 0xffffffffu);
                double s = 0.0;
                if (n1.x * n2.x <= 0.0) s += 0.35;
                if (n1.y * n2.y <= 0.0) s += 0.35;
                if (n1.z * n2.z <= 0.0) s += 0.35;
                double add = 0.55 * sqrt(min(4.0, dih)) + s;
                sil[a] += add;
                sil[b] += add;
            }
            i = j;
        }

        for (int i = 0; i < N0; ++i) {
            double d = min(2.0, max(0.0, dev[i]));
            double s = min(4.0, sil[i]);
            importance[i] = 1.0 + 1.8 * sqrt(d) + 0.45 * s;
            if (!isfinite(importance[i])) importance[i] = 1.0;
        }
    }

    void reset_stamps_if_needed() {
        if (stamp < 2000000000) return;
        fill(mark.begin(), mark.end(), 0);
        fill(mark2.begin(), mark2.end(), 0);
        stamp = 1;
    }

    bool edge_common(int a, int b, int common[2], int opp[2]) {
        auto& La = clean_inc(a);
        int cnt = 0;
        for (int fid : La) {
            if (!faces[fid].alive) continue;
            const Face& f = faces[fid];
            if (f.v[0] == b || f.v[1] == b || f.v[2] == b) {
                if (cnt >= 2) return false;
                common[cnt] = fid;
                opp[cnt] = other_in_face(fid, a, b);
                ++cnt;
            }
        }
        return cnt == 2 && opp[0] >= 0 && opp[1] >= 0 && opp[0] != opp[1];
    }

    void collect_neighbors_marked(int v, int token, int skip) {
        auto& L = clean_inc(v);
        for (int fid : L) {
            const Face& f = faces[fid];
            for (int k = 0; k < 3; ++k) {
                int w = f.v[k];
                if (w != v && w != skip && alive_v[w]) mark[w] = token;
            }
        }
    }

    bool link_condition(int a, int b, int oa, int ob) {
        reset_stamps_if_needed();
        int tok = ++stamp;
        int tok2 = ++stamp;
        collect_neighbors_marked(a, tok, b);
        vector<int> inter;
        auto& Lb = clean_inc(b);
        for (int fid : Lb) {
            const Face& f = faces[fid];
            for (int k = 0; k < 3; ++k) {
                int w = f.v[k];
                if (w == a || w == b || !alive_v[w]) continue;
                if (mark[w] == tok && mark2[w] != tok2) {
                    mark2[w] = tok2;
                    inter.push_back(w);
                    if (inter.size() > 2) return false;
                }
            }
        }
        if (inter.size() != 2) return false;
        return ((inter[0] == oa && inter[1] == ob) ||
                (inter[0] == ob && inter[1] == oa));
    }

    bool new_face_vertices_after_replace(int fid, int drop, int keep, int out[3]) const {
        const Face& f = faces[fid];
        for (int i = 0; i < 3; ++i) out[i] = (f.v[i] == drop ? keep : f.v[i]);
        return out[0] != out[1] && out[1] != out[2] && out[0] != out[2];
    }

    bool would_duplicate_local(int fid, const int nv[3], int keep, int c0, int c1) {
        int s0[3] = {nv[0], nv[1], nv[2]};
        sort(s0, s0 + 3);
        auto& L = clean_inc(keep);
        for (int ofid : L) {
            if (ofid == fid || ofid == c0 || ofid == c1 || !faces[ofid].alive) continue;
            int s1[3] = {faces[ofid].v[0], faces[ofid].v[1], faces[ofid].v[2]};
            sort(s1, s1 + 3);
            if (s0[0] == s1[0] && s0[1] == s1[1] && s0[2] == s1[2]) return true;
        }
        return false;
    }

    double local_radius_limit(int keep, int drop) const {
        if (stage >= 2) return haus_eps;
        double imp = max(importance[keep], importance[drop]);
        double factor;
        if (stage == 0) {
            factor = 1.0 / (1.0 + 0.075 * max(0.0, imp - 1.0));
            factor = max(0.70, min(1.0, factor));
        } else {
            factor = 1.0 / (1.0 + 0.035 * max(0.0, imp - 1.0));
            factor = max(0.84, min(1.0, factor));
        }
        return haus_eps * factor;
    }

    double normal_threshold(int keep, int drop) const {
        double imp = max(importance[keep], importance[drop]);
        double base = (stage == 0 ? 0.015 : (stage == 1 ? -0.045 : -0.12));
        double add = min(0.18, 0.018 * max(0.0, imp - 1.0));
        return base + add;
    }

    bool valid_direction(int keep, int drop, int common[2], int opp[2]) {
        if (!alive_v[keep] || !alive_v[drop] || keep == drop) return false;
        double d = dist3(p[keep], p[drop]);
        if (!isfinite(d)) return false;
        if (max(rad[keep], rad[drop] + d) > local_radius_limit(keep, drop) + 1e-12) return false;
        if (!link_condition(keep, drop, opp[0], opp[1])) return false;

        double nmin = normal_threshold(keep, drop);
        vector<int> L = clean_inc(drop);
        for (int fid : L) {
            if (!faces[fid].alive || fid == common[0] || fid == common[1]) continue;
            int nv[3];
            if (!new_face_vertices_after_replace(fid, drop, keep, nv)) return false;
            Vec3 rn = raw_face_normal_from_vertices(nv[0], nv[1], nv[2]);
            double a2 = norm2(rn);
            if (!(a2 > area_eps2) || !isfinite(a2)) return false;
            Vec3 oldn = normalized(raw_face_normal(fid));
            Vec3 newn = normalized(rn);
            double nd = dot3(oldn, newn);
            if (!isfinite(nd) || nd < nmin) return false;
            if (would_duplicate_local(fid, nv, keep, common[0], common[1])) return false;
        }
        return true;
    }

    bool choose_direction(int a, int b, int& keep, int& drop, int common[2], int opp[2]) {
        if (!alive_v[a] || !alive_v[b]) return false;
        if (!edge_common(a, b, common, opp)) return false;

        bool va = valid_direction(a, b, common, opp);
        bool vb = valid_direction(b, a, common, opp);
        if (!va && !vb) return false;
        if (va && !vb) { keep = a; drop = b; return true; }
        if (!va && vb) { keep = b; drop = a; return true; }

        double da = dist3(p[a], p[b]);
        double scoreA = max(rad[a], rad[b] + da) / max(1e-30, local_radius_limit(a, b));
        double scoreB = max(rad[b], rad[a] + da) / max(1e-30, local_radius_limit(b, a));
        scoreA -= 0.045 * importance[a];
        scoreB -= 0.045 * importance[b];
        if (scoreA < scoreB) { keep = a; drop = b; }
        else { keep = b; drop = a; }
        return true;
    }

    void push_candidate(int a, int b) {
        if (a == b || a < 0 || b < 0 || a >= N0 || b >= N0) return;
        if (!alive_v[a] || !alive_v[b]) return;
        double len = dist3(p[a], p[b]);
        if (!isfinite(len)) return;
        double imp = 0.5 * (importance[a] + importance[b]);
        double rr = (rad[a] + rad[b]) / max(1e-30, haus_eps);
        double c = (len / max(1e-30, diag)) * (1.0 + 0.22 * imp) + 0.025 * rr;
        if (!isfinite(c)) return;
        pq.push(Candidate{(float)c, a, b});
    }

    vector<uint64_t> current_edges() {
        vector<uint64_t> e;
        e.reserve((size_t)alive_faces * 3);
        for (int fid = 0; fid < (int)faces.size(); ++fid) {
            if (!faces[fid].alive) continue;
            const Face& f = faces[fid];
            e.push_back(edge_key_int(f.v[0], f.v[1]));
            e.push_back(edge_key_int(f.v[1], f.v[2]));
            e.push_back(edge_key_int(f.v[2], f.v[0]));
        }
        sort(e.begin(), e.end());
        e.erase(unique(e.begin(), e.end()), e.end());
        return e;
    }

    void rebuild_queue() {
        priority_queue<Candidate> empty;
        pq.swap(empty);
        vector<uint64_t> e = current_edges();
        for (uint64_t key : e) {
            int a = (int)(key >> 32), b = (int)(key & 0xffffffffu);
            push_candidate(a, b);
        }
    }

    void push_neighbors_of(int v) {
        reset_stamps_if_needed();
        int tok = ++stamp;
        auto& L = clean_inc(v);
        for (int fid : L) {
            const Face& f = faces[fid];
            for (int k = 0; k < 3; ++k) {
                int w = f.v[k];
                if (w != v && alive_v[w] && mark[w] != tok) {
                    mark[w] = tok;
                    push_candidate(v, w);
                }
            }
        }
    }

    void collapse(int keep, int drop, int common[2]) {
        double d = dist3(p[keep], p[drop]);
        vector<int> L = clean_inc(drop);
        faces[common[0]].alive = 0;
        faces[common[1]].alive = 0;
        alive_faces -= 2;
        for (int fid : L) {
            if (!faces[fid].alive) continue;
            Face& f = faces[fid];
            for (int k = 0; k < 3; ++k) if (f.v[k] == drop) f.v[k] = keep;
            inc[keep].push_back(fid);
        }
        alive_v[drop] = 0;
        --alive_vertices;
        rad[keep] = max(rad[keep], rad[drop] + d);
        if (importance[drop] > importance[keep]) importance[keep] = max(importance[keep], 0.92 * importance[drop]);
        clean_inc(keep);
        inc[drop].clear();
        push_neighbors_of(keep);
    }

    double target_ratio() const {
        // The public score is 100 * (1 - V_out / V_in).  These targets are
        // intentionally aggressive: for the main-sized cases they aim around
        // 91.5--92.2 before validity/visual guards stop the collapse stream.
        if (N0 >= 700000) return 0.0775;
        if (N0 >= 200000) return 0.0785;
        if (N0 >= 80000) return 0.0795;
        if (N0 >= 35000) return 0.0805;
        if (N0 >= 12000) return 0.0875;
        if (N0 >= 4000) return 0.1100;
        return 0.18;
    }

    bool time_limit_near() const {
        auto now = chrono::steady_clock::now();
        double sec = chrono::duration<double>(now - t0).count();
        return sec > 19.25;
    }

    void run() {
        t0 = chrono::steady_clock::now();
        int target = max(4, (int)ceil(N0 * target_ratio()));
        if (N0 <= 4) return;

        for (stage = 0; stage <= 2 && alive_vertices > target; ++stage) {
            rebuild_queue();
            long long pops = 0, collapses = 0;
            while (!pq.empty() && alive_vertices > target) {
                Candidate c = pq.top(); pq.pop();
                if (!alive_v[c.a] || !alive_v[c.b] || c.a == c.b) continue;
                int common[2], opp[2], keep = -1, drop = -1;
                if (choose_direction(c.a, c.b, keep, drop, common, opp)) {
                    collapse(keep, drop, common);
                    ++collapses;
                }
                if ((++pops & 8191LL) == 0 && time_limit_near()) return;
                if ((pops & 262143LL) == 0 && collapses == 0 && stage < 2) break;
            }
            if (time_limit_near()) break;
        }
    }

    bool validate_active() {
        vector<uint64_t> edge_keys;
        edge_keys.reserve((size_t)alive_faces * 3);
        vector<array<int,3>> tris;
        tris.reserve(alive_faces);
        vector<int> used;
        used.reserve(alive_vertices);
        vector<unsigned char> is_used(N0, 0);

        for (int fid = 0; fid < (int)faces.size(); ++fid) {
            if (!faces[fid].alive) continue;
            const Face& f = faces[fid];
            int a = f.v[0], b = f.v[1], c = f.v[2];
            if (a < 0 || b < 0 || c < 0 || a >= N0 || b >= N0 || c >= N0) return false;
            if (!alive_v[a] || !alive_v[b] || !alive_v[c]) return false;
            if (a == b || b == c || a == c) return false;
            Vec3 rn = raw_face_normal_from_vertices(a, b, c);
            if (!(norm2(rn) > area_eps2) || !isfinite(norm2(rn))) return false;
            array<int,3> t = {a,b,c};
            array<int,3> st = t;
            sort(st.begin(), st.end());
            tris.push_back(st);
            edge_keys.push_back(edge_key_int(a,b));
            edge_keys.push_back(edge_key_int(b,c));
            edge_keys.push_back(edge_key_int(c,a));
            if (!is_used[a]) { is_used[a] = 1; used.push_back(a); }
            if (!is_used[b]) { is_used[b] = 1; used.push_back(b); }
            if (!is_used[c]) { is_used[c] = 1; used.push_back(c); }
        }
        if (tris.empty()) return false;
        sort(tris.begin(), tris.end());
        for (size_t i = 1; i < tris.size(); ++i) if (tris[i] == tris[i-1]) return false;
        sort(edge_keys.begin(), edge_keys.end());
        for (size_t i = 0; i < edge_keys.size();) {
            size_t j = i + 1;
            while (j < edge_keys.size() && edge_keys[j] == edge_keys[i]) ++j;
            if (j - i != 2) return false;
            i = j;
        }

        // Local vertex-link check.  This catches pinched vertices that can pass the edge-count test.
        for (int v : used) {
            auto& L = clean_inc(v);
            if (L.empty()) return false;
            vector<pair<int,int>> link_edges;
            vector<int> neigh;
            link_edges.reserve(L.size());
            neigh.reserve(L.size() * 2);
            for (int fid : L) {
                const Face& f = faces[fid];
                int a = -1, b = -1;
                for (int k = 0; k < 3; ++k) if (f.v[k] != v) {
                    if (a < 0) a = f.v[k]; else b = f.v[k];
                }
                if (a < 0 || b < 0 || a == b) return false;
                link_edges.emplace_back(a, b);
                neigh.push_back(a);
                neigh.push_back(b);
            }
            sort(neigh.begin(), neigh.end());
            neigh.erase(unique(neigh.begin(), neigh.end()), neigh.end());
            int k = (int)neigh.size();
            if (k < 3) return false;
            vector<int> deg(k, 0);
            vector<vector<int>> adj(k);
            auto get_id = [&](int x) {
                return (int)(lower_bound(neigh.begin(), neigh.end(), x) - neigh.begin());
            };
            for (auto [a,b] : link_edges) {
                int ia = get_id(a), ib = get_id(b);
                if (ia < 0 || ia >= k || ib < 0 || ib >= k || ia == ib) return false;
                deg[ia]++; deg[ib]++;
                adj[ia].push_back(ib);
                adj[ib].push_back(ia);
            }
            for (int d : deg) if (d != 2) return false;
            vector<char> seen(k, 0);
            int cnt = 0;
            vector<int> st = {0};
            seen[0] = 1;
            while (!st.empty()) {
                int x = st.back(); st.pop_back(); ++cnt;
                for (int y : adj[x]) if (!seen[y]) { seen[y] = 1; st.push_back(y); }
            }
            if (cnt != k) return false;
        }
        return true;
    }

    void output_original(ostream& out) const {
        out << N0 << ' ' << M0 << '\n';
        out.setf(ios::fmtflags(0), ios::floatfield);
        out << setprecision(15);
        for (const auto& q : p) out << q.x << ' ' << q.y << ' ' << q.z << '\n';
        for (const auto& f : original_faces) out << f[0] + 1 << ' ' << f[1] + 1 << ' ' << f[2] + 1 << '\n';
    }

    void output_active(ostream& out) {
        vector<int> id(N0, -1);
        vector<int> verts;
        verts.reserve(alive_vertices);
        vector<array<int,3>> out_faces;
        out_faces.reserve(alive_faces);
        for (int fid = 0; fid < (int)faces.size(); ++fid) {
            if (!faces[fid].alive) continue;
            int ov[3];
            for (int k = 0; k < 3; ++k) {
                int v = faces[fid].v[k];
                if (id[v] < 0) { id[v] = (int)verts.size(); verts.push_back(v); }
                ov[k] = id[v];
            }
            out_faces.push_back({ov[0], ov[1], ov[2]});
        }
        out << verts.size() << ' ' << out_faces.size() << '\n';
        out.setf(ios::fmtflags(0), ios::floatfield);
        out << setprecision(15);
        for (int v : verts) out << p[v].x << ' ' << p[v].y << ' ' << p[v].z << '\n';
        for (const auto& f : out_faces) out << f[0] + 1 << ' ' << f[1] + 1 << ' ' << f[2] + 1 << '\n';
    }
};

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    int V, F;
    if (!(cin >> V >> F)) return 0;
    if (V <= 0 || F <= 0) return 0;
    vector<Vec3> p(V);
    for (int i = 0; i < V; ++i) cin >> p[i].x >> p[i].y >> p[i].z;
    vector<array<int,3>> faces(F);
    bool ok = true;
    for (int i = 0; i < F; ++i) {
        long long a,b,c;
        cin >> a >> b >> c;
        --a; --b; --c;
        if (a < 0 || b < 0 || c < 0 || a >= V || b >= V || c >= V) ok = false;
        faces[i] = {(int)a, (int)b, (int)c};
    }
    if (!ok) {
        cout << V << ' ' << F << '\n';
        cout << setprecision(15);
        for (const auto& q : p) cout << q.x << ' ' << q.y << ' ' << q.z << '\n';
        for (const auto& f : faces) cout << f[0] + 1 << ' ' << f[1] + 1 << ' ' << f[2] + 1 << '\n';
        return 0;
    }

    Simplifier sim(std::move(p), std::move(faces));
    sim.run();
    if (sim.validate_active()) sim.output_active(cout);
    else sim.output_original(cout);
    return 0;
}
