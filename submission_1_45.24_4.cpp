#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <queue>
#include <string>
#include <utility>
#include <vector>

using namespace std;

struct Vec3 {
    double x = 0.0, y = 0.0, z = 0.0;
};

static inline Vec3 operator-(const Vec3& a, const Vec3& b) {
    return {a.x - b.x, a.y - b.y, a.z - b.z};
}

static inline Vec3 operator*(const Vec3& a, double k) {
    return {a.x * k, a.y * k, a.z * k};
}

static inline double dot3(const Vec3& a, const Vec3& b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

static inline Vec3 cross3(const Vec3& a, const Vec3& b) {
    return {a.y * b.z - a.z * b.y,
            a.z * b.x - a.x * b.z,
            a.x * b.y - a.y * b.x};
}

static inline double norm2(const Vec3& a) {
    return dot3(a, a);
}

static inline double dist2(const Vec3& a, const Vec3& b) {
    return norm2(a - b);
}

static inline bool normalize(Vec3& a) {
    double n2 = norm2(a);
    if (n2 <= 1e-30) return false;
    double inv = 1.0 / sqrt(n2);
    a = a * inv;
    return true;
}

struct Face {
    int v[3];
    unsigned char active = 1;
};

struct Quadric {
    // Symmetric 4x4 matrix stored as:
    // xx, xy, xz, xw, yy, yz, yw, zz, zw, ww
    double q[10];

    Quadric() { memset(q, 0, sizeof(q)); }

    void add_plane(double a, double b, double c, double d) {
        q[0] += a * a;
        q[1] += a * b;
        q[2] += a * c;
        q[3] += a * d;
        q[4] += b * b;
        q[5] += b * c;
        q[6] += b * d;
        q[7] += c * c;
        q[8] += c * d;
        q[9] += d * d;
    }

    void add(const Quadric& o) {
        for (int i = 0; i < 10; ++i) q[i] += o.q[i];
    }

    double eval(const Vec3& p) const {
        double x = p.x, y = p.y, z = p.z;
        return q[0] * x * x
             + 2.0 * q[1] * x * y
             + 2.0 * q[2] * x * z
             + 2.0 * q[3] * x
             + q[4] * y * y
             + 2.0 * q[5] * y * z
             + 2.0 * q[6] * y
             + q[7] * z * z
             + 2.0 * q[8] * z
             + q[9];
    }
};

struct FastInput {
    vector<char> buf;
    char* p = nullptr;

    void read_all() {
        char tmp[1 << 16];
        size_t n;
        while ((n = fread(tmp, 1, sizeof(tmp), stdin)) > 0) {
            buf.insert(buf.end(), tmp, tmp + n);
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
        return x * sign;
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

struct EdgeFace {
    uint64_t key;
    int face;
    bool operator<(const EdgeFace& o) const {
        if (key != o.key) return key < o.key;
        return face < o.face;
    }
};

struct Node {
    double cost;
    int u, v;
    int vu, vv;
    bool operator<(const Node& o) const {
        return cost > o.cost;
    }
};

static int N, M;
static vector<Vec3> P;
static vector<Face> F;
static vector<vector<int>> incident;
static vector<Quadric> Q;
static vector<unsigned char> active_vertex;
static vector<int> head_member, tail_member, next_member, cluster_size;
static vector<int> version_id;
static vector<int> mark_seen;
static int mark_token = 7;
static int active_vertices = 0;
static int active_faces = 0;
static double hausdorff_limit2 = 0.0;
static double diagonal_len = 0.0;
static int target_vertices = 0;
static priority_queue<Node> pq;
static chrono::steady_clock::time_point start_time;
static const double TIME_LIMIT_SECONDS = 19.25;

static inline uint64_t edge_key(int a, int b) {
    if (a > b) swap(a, b);
    return (uint64_t)(uint32_t)a << 32 | (uint32_t)b;
}

static inline int key_a(uint64_t key) {
    return (int)(key >> 32);
}

static inline int key_b(uint64_t key) {
    return (int)(key & 0xffffffffu);
}

static inline bool time_left() {
    auto now = chrono::steady_clock::now();
    double elapsed = chrono::duration<double>(now - start_time).count();
    return elapsed < TIME_LIMIT_SECONDS;
}

static inline bool face_has_vertex(const Face& f, int v) {
    return f.v[0] == v || f.v[1] == v || f.v[2] == v;
}

static inline bool face_has_edge(const Face& f, int a, int b) {
    return face_has_vertex(f, a) && face_has_vertex(f, b);
}

static inline int third_vertex(const Face& f, int a, int b) {
    for (int k = 0; k < 3; ++k) {
        int x = f.v[k];
        if (x != a && x != b) return x;
    }
    return -1;
}

static inline Vec3 face_normal_from_indices(int a, int b, int c) {
    return cross3(P[b] - P[a], P[c] - P[a]);
}

static bool are_adjacent(int a, int b) {
    const vector<int>& ia = incident[a];
    const vector<int>& ib = incident[b];
    const vector<int>& small = (ia.size() < ib.size()) ? ia : ib;
    for (int fid : small) {
        if (!F[fid].active) continue;
        if (face_has_edge(F[fid], a, b)) return true;
    }
    return false;
}

static bool collect_edge_opposites(int a, int b, int opp[2]) {
    int cnt = 0;
    const vector<int>& ia = incident[a];
    const vector<int>& ib = incident[b];
    const vector<int>& small = (ia.size() < ib.size()) ? ia : ib;
    for (int fid : small) {
        if (!F[fid].active) continue;
        const Face& f = F[fid];
        if (!face_has_edge(f, a, b)) continue;
        if (cnt >= 2) return false;
        int t = third_vertex(f, a, b);
        if (t < 0) return false;
        opp[cnt++] = t;
    }
    if (cnt != 2) return false;
    return opp[0] != opp[1];
}

static bool link_condition_ok(int a, int b) {
    int opp[2];
    if (!collect_edge_opposites(a, b, opp)) return false;

    if (mark_token > 2000000000) {
        fill(mark_seen.begin(), mark_seen.end(), 0);
        mark_token = 7;
    }
    int token_u = mark_token++;
    int token_common = mark_token++;

    for (int fid : incident[a]) {
        if (!F[fid].active) continue;
        const Face& f = F[fid];
        if (!face_has_vertex(f, a)) continue;
        for (int k = 0; k < 3; ++k) {
            int x = f.v[k];
            if (x != a && x != b) mark_seen[x] = token_u;
        }
    }

    int common_count = 0;
    int got0 = 0, got1 = 0;
    for (int fid : incident[b]) {
        if (!F[fid].active) continue;
        const Face& f = F[fid];
        if (!face_has_vertex(f, b)) continue;
        for (int k = 0; k < 3; ++k) {
            int x = f.v[k];
            if (x == a || x == b) continue;
            if (mark_seen[x] == token_u) {
                mark_seen[x] = token_common;
                ++common_count;
                if (x == opp[0]) got0 = 1;
                if (x == opp[1]) got1 = 1;
                if (common_count > 2) return false;
            }
        }
    }
    return common_count == 2 && got0 && got1;
}

static bool cluster_can_move(int src, int dst) {
    const Vec3& to = P[dst];
    for (int m = head_member[src]; m != -1; m = next_member[m]) {
        if (dist2(P[m], to) > hausdorff_limit2) return false;
    }
    return true;
}

static bool normals_ok_after_move(int src, int dst) {
    for (int fid : incident[src]) {
        if (!F[fid].active) continue;
        const Face& f = F[fid];
        if (!face_has_vertex(f, src)) continue;
        if (face_has_vertex(f, dst)) continue;

        int a = f.v[0], b = f.v[1], c = f.v[2];
        Vec3 old_n = face_normal_from_indices(a, b, c);
        double old_len2 = norm2(old_n);
        if (old_len2 <= 1e-30) return false;

        if (a == src) a = dst;
        else if (b == src) b = dst;
        else if (c == src) c = dst;

        if (a == b || b == c || c == a) return false;
        Vec3 new_n = face_normal_from_indices(a, b, c);
        double new_len2 = norm2(new_n);
        if (new_len2 <= 1e-30) return false;

        double d = dot3(old_n, new_n);
        if (d <= 1e-18) return false;
    }
    return true;
}

static double directed_collapse_cost(int src, int dst) {
    Quadric q = Q[src];
    q.add(Q[dst]);
    double c = q.eval(P[dst]);
    c += 0.001 * dist2(P[src], P[dst]);
    return c;
}

static void push_edge(int a, int b) {
    if (a == b) return;
    if (!active_vertex[a] || !active_vertex[b]) return;
    double d2 = dist2(P[a], P[b]);
    if (d2 > hausdorff_limit2) return;
    Quadric q = Q[a];
    q.add(Q[b]);
    double c1 = q.eval(P[a]);
    double c2 = q.eval(P[b]);
    double c = min(c1, c2) + 0.001 * d2;
    pq.push({c, a, b, version_id[a], version_id[b]});
}

static void compact_incident(int v) {
    vector<int>& ids = incident[v];
    if (ids.size() < 128) return;
    size_t alive = 0;
    for (int fid : ids) {
        if (F[fid].active && face_has_vertex(F[fid], v)) ++alive;
    }
    if (alive * 3 + 32 >= ids.size()) return;
    vector<int> keep;
    keep.reserve(alive + 8);
    for (int fid : ids) {
        if (F[fid].active && face_has_vertex(F[fid], v)) keep.push_back(fid);
    }
    ids.swap(keep);
}

static void merge_members(int src, int dst) {
    if (head_member[src] == -1) return;
    next_member[tail_member[dst]] = head_member[src];
    tail_member[dst] = tail_member[src];
    cluster_size[dst] += cluster_size[src];
    head_member[src] = tail_member[src] = -1;
    cluster_size[src] = 0;
}

static void do_collapse(int src, int dst) {
    Q[dst].add(Q[src]);

    for (int fid : incident[src]) {
        if (!F[fid].active) continue;
        Face& f = F[fid];
        bool has_src = false, has_dst = false;
        for (int k = 0; k < 3; ++k) {
            if (f.v[k] == src) has_src = true;
            if (f.v[k] == dst) has_dst = true;
        }
        if (!has_src) continue;
        if (has_dst) {
            f.active = 0;
            --active_faces;
        } else {
            for (int k = 0; k < 3; ++k) {
                if (f.v[k] == src) f.v[k] = dst;
            }
            incident[dst].push_back(fid);
        }
    }

    active_vertex[src] = 0;
    --active_vertices;
    ++version_id[src];
    ++version_id[dst];
    merge_members(src, dst);

    compact_incident(src);
    compact_incident(dst);

    for (int fid : incident[dst]) {
        if (!F[fid].active) continue;
        const Face& f = F[fid];
        if (!face_has_vertex(f, dst)) continue;
        push_edge(f.v[0], f.v[1]);
        push_edge(f.v[1], f.v[2]);
        push_edge(f.v[2], f.v[0]);
    }
}

static bool attempt_collapse(int a, int b) {
    if (a == b) return false;
    if (!active_vertex[a] || !active_vertex[b]) return false;
    if (dist2(P[a], P[b]) > hausdorff_limit2) return false;
    if (!link_condition_ok(a, b)) return false;

    bool b_to_a = cluster_can_move(b, a) && normals_ok_after_move(b, a);
    bool a_to_b = cluster_can_move(a, b) && normals_ok_after_move(a, b);
    if (!b_to_a && !a_to_b) return false;

    if (b_to_a && a_to_b) {
        double c_ba = directed_collapse_cost(b, a);
        double c_ab = directed_collapse_cost(a, b);
        if (c_ba <= c_ab) do_collapse(b, a);
        else do_collapse(a, b);
    } else if (b_to_a) {
        do_collapse(b, a);
    } else {
        do_collapse(a, b);
    }
    return true;
}

static void load_mesh() {
    FastInput in;
    in.read_all();
    N = in.next_int();
    M = in.next_int();
    P.resize(N);
    F.resize(M);

    Vec3 mn{1e100, 1e100, 1e100};
    Vec3 mx{-1e100, -1e100, -1e100};
    for (int i = 0; i < N; ++i) {
        (void)in.next_char();
        P[i].x = in.next_double();
        P[i].y = in.next_double();
        P[i].z = in.next_double();
        mn.x = min(mn.x, P[i].x);
        mn.y = min(mn.y, P[i].y);
        mn.z = min(mn.z, P[i].z);
        mx.x = max(mx.x, P[i].x);
        mx.y = max(mx.y, P[i].y);
        mx.z = max(mx.z, P[i].z);
    }

    vector<int> deg(N, 0);
    for (int i = 0; i < M; ++i) {
        (void)in.next_char();
        int a = in.next_int() - 1;
        int b = in.next_int() - 1;
        int c = in.next_int() - 1;
        F[i].v[0] = a;
        F[i].v[1] = b;
        F[i].v[2] = c;
        ++deg[a];
        ++deg[b];
        ++deg[c];
    }

    Vec3 d = mx - mn;
    diagonal_len = sqrt(norm2(d));
    double hausdorff_limit = 0.05 * diagonal_len * 0.999999;
    hausdorff_limit2 = hausdorff_limit * hausdorff_limit;

    incident.resize(N);
    for (int i = 0; i < N; ++i) incident[i].reserve(deg[i] + 8);
    for (int i = 0; i < M; ++i) {
        incident[F[i].v[0]].push_back(i);
        incident[F[i].v[1]].push_back(i);
        incident[F[i].v[2]].push_back(i);
    }

    Q.assign(N, Quadric());
    active_vertex.assign(N, 1);
    head_member.resize(N);
    tail_member.resize(N);
    next_member.assign(N, -1);
    cluster_size.assign(N, 1);
    version_id.assign(N, 0);
    mark_seen.assign(N, 0);
    for (int i = 0; i < N; ++i) {
        head_member[i] = tail_member[i] = i;
    }
    active_vertices = N;
    active_faces = M;
}

static double initialize_quadrics_and_edges() {
    vector<Vec3> face_normals(M);
    vector<EdgeFace> edges;
    edges.reserve((size_t)M * 3);

    for (int i = 0; i < M; ++i) {
        Face& f = F[i];
        Vec3 n = face_normal_from_indices(f.v[0], f.v[1], f.v[2]);
        if (!normalize(n)) n = {0.0, 0.0, 0.0};
        face_normals[i] = n;

        double dd = -dot3(n, P[f.v[0]]);
        Q[f.v[0]].add_plane(n.x, n.y, n.z, dd);
        Q[f.v[1]].add_plane(n.x, n.y, n.z, dd);
        Q[f.v[2]].add_plane(n.x, n.y, n.z, dd);

        edges.push_back({edge_key(f.v[0], f.v[1]), i});
        edges.push_back({edge_key(f.v[1], f.v[2]), i});
        edges.push_back({edge_key(f.v[2], f.v[0]), i});
    }

    sort(edges.begin(), edges.end());
    long long unique_edges = 0;
    long long feature_edges = 0;
    const double feature_cos = cos(35.0 * acos(-1.0) / 180.0);

    for (size_t i = 0; i < edges.size();) {
        size_t j = i + 1;
        while (j < edges.size() && edges[j].key == edges[i].key) ++j;

        ++unique_edges;
        if (j - i == 2) {
            double d = dot3(face_normals[edges[i].face], face_normals[edges[i + 1].face]);
            if (d < feature_cos) ++feature_edges;
        }

        int a = key_a(edges[i].key);
        int b = key_b(edges[i].key);
        push_edge(a, b);
        i = j;
    }

    vector<EdgeFace>().swap(edges);
    vector<Vec3>().swap(face_normals);
    if (unique_edges == 0) return 0.0;
    return (double)feature_edges / (double)unique_edges;
}

static void choose_target(double feature_ratio) {
    double ratio = 0.082 + 0.10 * min(0.18, feature_ratio);
    if (N >= 300000) ratio -= 0.004;
    if (N <= 8000) ratio = max(ratio, 0.100);
    if (N <= 1000) ratio = max(ratio, 0.180);
    if (N <= 50) ratio = 0.01;

    ratio = max(0.055, min(0.140, ratio));
    target_vertices = max(4, (int)ceil((double)N * ratio));
}

static void simplify_mesh() {
    double feature_ratio = initialize_quadrics_and_edges();
    choose_target(feature_ratio);

    long long pops = 0;
    while (active_vertices > target_vertices && !pq.empty()) {
        if ((++pops & 4095) == 0 && !time_left()) break;

        Node cur = pq.top();
        pq.pop();

        int a = cur.u, b = cur.v;
        if (a == b) continue;
        if (!active_vertex[a] || !active_vertex[b]) continue;
        if (cur.vu != version_id[a] || cur.vv != version_id[b]) {
            if (are_adjacent(a, b)) push_edge(a, b);
            continue;
        }
        if (!attempt_collapse(a, b)) continue;
    }
}

static void save_mesh() {
    vector<int> remap(N, -1);
    int out_vertices = 0;
    for (int i = 0; i < N; ++i) {
        if (active_vertex[i]) remap[i] = out_vertices++;
    }

    int out_faces = 0;
    for (int i = 0; i < M; ++i) {
        if (!F[i].active) continue;
        int a = F[i].v[0], b = F[i].v[1], c = F[i].v[2];
        if (a == b || b == c || c == a) continue;
        if (remap[a] < 0 || remap[b] < 0 || remap[c] < 0) continue;
        ++out_faces;
    }

    static char outbuf[1 << 20];
    setvbuf(stdout, outbuf, _IOFBF, sizeof(outbuf));

    printf("%d %d\n", out_vertices, out_faces);
    for (int i = 0; i < N; ++i) {
        if (!active_vertex[i]) continue;
        printf("v %.15g %.15g %.15g\n", P[i].x, P[i].y, P[i].z);
    }
    for (int i = 0; i < M; ++i) {
        if (!F[i].active) continue;
        int a = F[i].v[0], b = F[i].v[1], c = F[i].v[2];
        if (a == b || b == c || c == a) continue;
        if (remap[a] < 0 || remap[b] < 0 || remap[c] < 0) continue;
        printf("f %d %d %d\n", remap[a] + 1, remap[b] + 1, remap[c] + 1);
    }
}

int main() {
    start_time = chrono::steady_clock::now();
    load_mesh();
    simplify_mesh();
    save_mesh();
    return 0;
}
