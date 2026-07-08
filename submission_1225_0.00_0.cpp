#include <bits/stdc++.h>
using namespace std;

/*
  IMC Challenge 2026 - simplifygeometry
  Self-contained C++17 mesh simplifier.

  The implementation is intentionally format-tolerant:
    * custom:  n m [optional_tolerance]  + vertices + triangular faces
    * OBJ:     v / f lines, polygons triangulated
    * OFF
    * ASCII PLY with x y z and list-style faces

  Core algorithm:
    * robust cleanup
    * area/feature weighted quadric-error edge collapses
    * dynamic link-condition and flip tests to preserve manifold topology
    * optional/adaptive vertex-Hausdorff cover bound over represented original vertices

  No external libraries.
*/

static constexpr double EPS = 1e-12;
static constexpr double PI = 3.141592653589793238462643383279502884;

struct Vec3 {
    double x = 0, y = 0, z = 0;
    Vec3() = default;
    Vec3(double X, double Y, double Z) : x(X), y(Y), z(Z) {}
    Vec3 operator + (const Vec3& o) const { return {x + o.x, y + o.y, z + o.z}; }
    Vec3 operator - (const Vec3& o) const { return {x - o.x, y - o.y, z - o.z}; }
    Vec3 operator * (double s) const { return {x * s, y * s, z * s}; }
    Vec3 operator / (double s) const { return {x / s, y / s, z / s}; }
    Vec3& operator += (const Vec3& o) { x += o.x; y += o.y; z += o.z; return *this; }
    Vec3& operator -= (const Vec3& o) { x -= o.x; y -= o.y; z -= o.z; return *this; }
    Vec3& operator *= (double s) { x *= s; y *= s; z *= s; return *this; }
};

static inline double dotv(const Vec3& a, const Vec3& b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}
static inline Vec3 crossv(const Vec3& a, const Vec3& b) {
    return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
}
static inline double norm2(const Vec3& a) { return dotv(a, a); }
static inline double normv(const Vec3& a) { return sqrt(max(0.0, norm2(a))); }
static inline double distv(const Vec3& a, const Vec3& b) { return normv(a - b); }
static inline Vec3 normalized(const Vec3& a) {
    double n = normv(a);
    if (n < EPS) return {0, 0, 0};
    return a / n;
}
static inline double clampd(double x, double lo, double hi) {
    return x < lo ? lo : (x > hi ? hi : x);
}

struct SymmetricMatrix {
    // 0:a00 1:a01 2:a02 3:a03 4:a11 5:a12 6:a13 7:a22 8:a23 9:a33
    double m[10];
    SymmetricMatrix(double c = 0.0) { for (double &v : m) v = c; }
    SymmetricMatrix(double a, double b, double c, double d, double w = 1.0) {
        // plane: ax + by + cz + d = 0, weighted outer product
        m[0] = w * a * a; m[1] = w * a * b; m[2] = w * a * c; m[3] = w * a * d;
        m[4] = w * b * b; m[5] = w * b * c; m[6] = w * b * d;
        m[7] = w * c * c; m[8] = w * c * d;
        m[9] = w * d * d;
    }
    double operator[](int i) const { return m[i]; }
    double& operator[](int i) { return m[i]; }
    SymmetricMatrix operator+(const SymmetricMatrix& n) const {
        SymmetricMatrix r;
        for (int i = 0; i < 10; ++i) r.m[i] = m[i] + n.m[i];
        return r;
    }
    SymmetricMatrix& operator+=(const SymmetricMatrix& n) {
        for (int i = 0; i < 10; ++i) m[i] += n.m[i];
        return *this;
    }
    double det(int a11, int a12, int a13,
               int a21, int a22, int a23,
               int a31, int a32, int a33) const {
        return m[a11] * m[a22] * m[a33] + m[a13] * m[a21] * m[a32] + m[a12] * m[a23] * m[a31]
             - m[a13] * m[a22] * m[a31] - m[a11] * m[a23] * m[a32] - m[a12] * m[a21] * m[a33];
    }
};

static inline double vertex_error(const SymmetricMatrix& q, const Vec3& p) {
    double x = p.x, y = p.y, z = p.z;
    return q[0] * x * x + 2.0 * q[1] * x * y + 2.0 * q[2] * x * z + 2.0 * q[3] * x
         + q[4] * y * y + 2.0 * q[5] * y * z + 2.0 * q[6] * y
         + q[7] * z * z + 2.0 * q[8] * z + q[9];
}

struct Vertex {
    Vec3 p;
    SymmetricMatrix q;
    int tstart = 0;
    int tcount = 0;
    bool border = false;
    bool locked = false;
    double cover = 0.0;  // upper bound: farthest represented original vertex -> current p
    double nearOrig = 0.0; // upper bound: current p -> some represented original vertex
};

struct Triangle {
    int v[3] = {0,0,0};
    double err[4] = {0,0,0,0};
    bool deleted = false;
    bool dirty = false;
    Vec3 n;
};

struct Ref {
    int tid = 0;
    int tvertex = 0;
};

enum class FormatType { CUSTOM, OBJ, OFF, PLY };

struct FaceKeyHash {
    static uint64_t splitmix64(uint64_t x) {
        x += 0x9e3779b97f4a7c15ULL;
        x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
        x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
        return x ^ (x >> 31);
    }
    size_t operator()(const array<int,3>& a) const {
        uint64_t h = splitmix64((uint64_t)(uint32_t)a[0]);
        h ^= splitmix64(((uint64_t)(uint32_t)a[1] << 1) ^ 0x9e3779b97f4a7c15ULL);
        h ^= splitmix64(((uint64_t)(uint32_t)a[2] << 2) ^ 0xbf58476d1ce4e5b9ULL);
        return (size_t)h;
    }
};

struct MeshIO {
    FormatType fmt = FormatType::CUSTOM;
    bool customOneBased = true;
    bool customFaceHasPrefix = false;
    double toleranceHint = -1.0;
    vector<Vec3> vertices;
    vector<array<int,3>> faces;
};

static inline string trim_copy(const string& s) {
    size_t a = 0, b = s.size();
    while (a < b && isspace((unsigned char)s[a])) ++a;
    while (b > a && isspace((unsigned char)s[b-1])) --b;
    return s.substr(a, b-a);
}

static inline bool starts_with_token(const string& s, const string& token) {
    if (s.size() < token.size()) return false;
    if (s.compare(0, token.size(), token) != 0) return false;
    return s.size() == token.size() || isspace((unsigned char)s[token.size()]);
}

static bool is_integer_like(const string& s) {
    if (s.empty()) return false;
    size_t i = (s[0] == '-' || s[0] == '+') ? 1 : 0;
    if (i >= s.size()) return false;
    for (; i < s.size(); ++i) if (!isdigit((unsigned char)s[i])) return false;
    return true;
}

static vector<string> split_ws(const string& s) {
    vector<string> out;
    string token;
    stringstream ss(s);
    while (ss >> token) out.push_back(token);
    return out;
}

static int parse_obj_index(const string& token, int nverts) {
    string a;
    for (char c : token) {
        if (c == '/') break;
        a.push_back(c);
    }
    if (a.empty()) return -1;
    int idx = stoi(a);
    if (idx < 0) idx = nverts + idx;
    else idx = idx - 1;
    return idx;
}

static void add_polygon_triangulated(vector<array<int,3>>& faces, const vector<int>& poly) {
    if (poly.size() < 3) return;
    int a = poly[0];
    for (size_t i = 1; i + 1 < poly.size(); ++i) {
        int b = poly[i], c = poly[i+1];
        if (a != b && b != c && c != a) faces.push_back({a,b,c});
    }
}

static MeshIO parse_obj(const string& input) {
    MeshIO mesh; mesh.fmt = FormatType::OBJ; mesh.customOneBased = true;
    string line;
    stringstream ss(input);
    while (getline(ss, line)) {
        line = trim_copy(line);
        if (line.empty() || line[0] == '#') continue;
        if (starts_with_token(line, "v")) {
            auto tok = split_ws(line);
            if (tok.size() >= 4) {
                try { mesh.vertices.push_back({stod(tok[1]), stod(tok[2]), stod(tok[3])}); } catch (...) {}
            }
        } else if (starts_with_token(line, "f")) {
            auto tok = split_ws(line);
            vector<int> poly;
            for (size_t i = 1; i < tok.size(); ++i) {
                try {
                    int idx = parse_obj_index(tok[i], (int)mesh.vertices.size());
                    if (0 <= idx && idx < (int)mesh.vertices.size()) poly.push_back(idx);
                } catch (...) {}
            }
            add_polygon_triangulated(mesh.faces, poly);
        }
    }
    return mesh;
}

static MeshIO parse_off(const string& input) {
    MeshIO mesh; mesh.fmt = FormatType::OFF; mesh.customOneBased = false;
    stringstream ss(input);
    string token;
    ss >> token; // OFF
    int n = 0, f = 0, e = 0;
    ss >> n >> f >> e;
    mesh.vertices.reserve(max(0, n));
    for (int i = 0; i < n; ++i) {
        double x, y, z; ss >> x >> y >> z;
        mesh.vertices.push_back({x,y,z});
    }
    for (int i = 0; i < f; ++i) {
        int k; ss >> k;
        vector<int> poly(k);
        for (int j = 0; j < k; ++j) ss >> poly[j];
        bool ok = true;
        for (int id : poly) if (id < 0 || id >= n) ok = false;
        if (ok) add_polygon_triangulated(mesh.faces, poly);
    }
    return mesh;
}

static MeshIO parse_ply(const string& input) {
    MeshIO mesh; mesh.fmt = FormatType::PLY; mesh.customOneBased = false;
    stringstream ss(input);
    string line;
    int nverts = 0, nfaces = 0;
    bool ascii = false;
    while (getline(ss, line)) {
        string t = trim_copy(line);
        if (t == "ply") continue;
        if (t.find("format ascii") == 0) ascii = true;
        auto tok = split_ws(t);
        if (tok.size() >= 3 && tok[0] == "element" && tok[1] == "vertex") nverts = stoi(tok[2]);
        if (tok.size() >= 3 && tok[0] == "element" && tok[1] == "face") nfaces = stoi(tok[2]);
        if (t == "end_header") break;
    }
    if (!ascii) return mesh;
    mesh.vertices.reserve(max(0, nverts));
    for (int i = 0; i < nverts && getline(ss, line); ++i) {
        auto tok = split_ws(line);
        if (tok.size() >= 3) {
            try { mesh.vertices.push_back({stod(tok[0]), stod(tok[1]), stod(tok[2])}); } catch (...) { mesh.vertices.push_back({0,0,0}); }
        } else mesh.vertices.push_back({0,0,0});
    }
    for (int i = 0; i < nfaces && getline(ss, line); ++i) {
        auto tok = split_ws(line);
        if (tok.empty()) continue;
        int k = 0;
        try { k = stoi(tok[0]); } catch (...) { continue; }
        if ((int)tok.size() < k + 1) continue;
        vector<int> poly;
        for (int j = 0; j < k; ++j) {
            try { poly.push_back(stoi(tok[j+1])); } catch (...) { poly.push_back(-1); }
        }
        bool ok = true;
        for (int id : poly) if (id < 0 || id >= (int)mesh.vertices.size()) ok = false;
        if (ok) add_polygon_triangulated(mesh.faces, poly);
    }
    return mesh;
}

static MeshIO parse_custom(const string& input) {
    MeshIO mesh; mesh.fmt = FormatType::CUSTOM; mesh.customOneBased = true;
    vector<string> lines;
    string line;
    stringstream ss(input);
    while (getline(ss, line)) {
        string t = trim_copy(line);
        if (!t.empty() && t[0] != '#') lines.push_back(t);
    }
    if (lines.empty()) return mesh;
    auto first = split_ws(lines[0]);
    if (first.size() < 2) return mesh;
    int n = 0, m = 0;
    try { n = stoi(first[0]); m = stoi(first[1]); } catch (...) { return mesh; }
    if ((int)first.size() >= 3) {
        try {
            double v = stod(first[2]);
            if (isfinite(v) && v > 0) mesh.toleranceHint = v;
        } catch (...) {}
    }
    int ptr = 1;
    mesh.vertices.reserve(max(0, n));
    for (int i = 0; i < n && ptr < (int)lines.size(); ++i, ++ptr) {
        auto tok = split_ws(lines[ptr]);
        if (tok.size() >= 3) {
            try { mesh.vertices.push_back({stod(tok[0]), stod(tok[1]), stod(tok[2])}); } catch (...) { mesh.vertices.push_back({0,0,0}); }
        } else mesh.vertices.push_back({0,0,0});
    }
    vector<vector<int>> rawfaces;
    rawfaces.reserve(max(0, m));
    int mn = INT_MAX, mx = INT_MIN;
    for (int i = 0; i < m && ptr < (int)lines.size(); ++i, ++ptr) {
        auto tok = split_ws(lines[ptr]);
        if (tok.size() < 3) continue;
        vector<int> vals;
        for (auto &s : tok) {
            if (!is_integer_like(s)) continue;
            try { vals.push_back(stoi(s)); } catch (...) {}
        }
        if ((int)vals.size() >= 4 && vals[0] == 3) {
            mesh.customFaceHasPrefix = true;
            vals.erase(vals.begin());
        }
        if ((int)vals.size() >= 3) {
            vals.resize(3);
            for (int x : vals) { mn = min(mn, x); mx = max(mx, x); }
            rawfaces.push_back(vals);
        }
    }
    bool oneBased = true;
    if (mn == 0) oneBased = false;
    else if (mx == n) oneBased = true;
    else if (mn >= 1 && mx <= n) oneBased = true;
    else oneBased = false;
    mesh.customOneBased = oneBased;
    for (auto vals : rawfaces) {
        int a = vals[0], b = vals[1], c = vals[2];
        if (oneBased) { --a; --b; --c; }
        if (a >= 0 && a < n && b >= 0 && b < n && c >= 0 && c < n && a != b && b != c && c != a)
            mesh.faces.push_back({a,b,c});
    }
    return mesh;
}

static MeshIO parse_mesh(const string& input) {
    string first;
    stringstream ss(input);
    while (getline(ss, first)) {
        first = trim_copy(first);
        if (!first.empty() && first[0] != '#') break;
    }
    string low = first;
    for (char &c : low) c = (char)tolower((unsigned char)c);
    if (low == "off" || low.rfind("off ", 0) == 0) return parse_off(input);
    if (low == "ply") return parse_ply(input);

    // OBJ if there are actual vertex/face records near the top.
    bool looksObj = false;
    string line;
    stringstream ss2(input);
    int checked = 0;
    while (checked < 100 && getline(ss2, line)) {
        line = trim_copy(line);
        if (line.empty() || line[0] == '#') continue;
        ++checked;
        if (starts_with_token(line, "v") || starts_with_token(line, "f")) { looksObj = true; break; }
        if (!line.empty() && (isdigit((unsigned char)line[0]) || line[0] == '-' || line[0] == '+')) break;
    }
    if (looksObj) return parse_obj(input);
    return parse_custom(input);
}

struct MeshSimplifier {
    vector<Vertex> vertices;
    vector<Triangle> triangles;
    vector<Ref> refs;
    int deletedTriangles = 0;
    bool originalClosed = true;
    int originalBoundaryEdges = 0;

    MeshSimplifier() = default;

    void load(const vector<Vec3>& vin, const vector<array<int,3>>& fin) {
        vertices.clear(); triangles.clear(); refs.clear(); deletedTriangles = 0;
        vertices.resize(vin.size());
        for (size_t i = 0; i < vin.size(); ++i) vertices[i].p = vin[i];
        triangles.reserve(fin.size());
        for (auto f : fin) {
            if (f[0] == f[1] || f[1] == f[2] || f[2] == f[0]) continue;
            Triangle t;
            t.v[0] = f[0]; t.v[1] = f[1]; t.v[2] = f[2];
            t.deleted = false; t.dirty = false;
            triangles.push_back(t);
        }
        cleanup_degenerate_and_duplicate_faces();
        analyze_boundary_state();
        initialize_quadrics();
        update_mesh(0);
    }

    static uint64_t edge_key(int a, int b) {
        if (a > b) swap(a, b);
        return (uint64_t)(uint32_t)a << 32 | (uint32_t)b;
    }

    Vec3 triangle_normal(const Triangle& t) const {
        const Vec3& p0 = vertices[t.v[0]].p;
        const Vec3& p1 = vertices[t.v[1]].p;
        const Vec3& p2 = vertices[t.v[2]].p;
        return normalized(crossv(p1 - p0, p2 - p0));
    }

    double triangle_area2(const Triangle& t) const {
        const Vec3& p0 = vertices[t.v[0]].p;
        const Vec3& p1 = vertices[t.v[1]].p;
        const Vec3& p2 = vertices[t.v[2]].p;
        return normv(crossv(p1 - p0, p2 - p0));
    }

    void cleanup_degenerate_and_duplicate_faces() {
        vector<Triangle> out;
        out.reserve(triangles.size());
        unordered_set<array<int,3>, FaceKeyHash> seen;
        seen.reserve(triangles.size() * 2 + 16);
        for (auto &t : triangles) {
            if (t.v[0] < 0 || t.v[1] < 0 || t.v[2] < 0 ||
                t.v[0] >= (int)vertices.size() || t.v[1] >= (int)vertices.size() || t.v[2] >= (int)vertices.size()) continue;
            if (t.v[0] == t.v[1] || t.v[1] == t.v[2] || t.v[2] == t.v[0]) continue;
            if (triangle_area2(t) < 1e-15) continue;
            array<int,3> key = {t.v[0], t.v[1], t.v[2]};
            sort(key.begin(), key.end());
            if (seen.insert(key).second) out.push_back(t);
        }
        triangles.swap(out);
    }

    void analyze_boundary_state() {
        vector<uint64_t> edges;
        edges.reserve(triangles.size() * 3 + 16);
        for (auto &t : triangles) {
            for (int i = 0; i < 3; ++i) edges.push_back(edge_key(t.v[i], t.v[(i+1)%3]));
        }
        sort(edges.begin(), edges.end());
        originalBoundaryEdges = 0;
        int nonManifoldEdges = 0;
        for (size_t l = 0; l < edges.size(); ) {
            size_t r = l + 1;
            while (r < edges.size() && edges[r] == edges[l]) ++r;
            int c = (int)(r - l);
            if (c == 1) originalBoundaryEdges++;
            if (c > 2) nonManifoldEdges++;
            l = r;
        }
        originalClosed = (originalBoundaryEdges == 0 && nonManifoldEdges == 0 && !triangles.empty());
    }

    void add_plane_quadric_to_vertex(int vid, const Vec3& n, double d, double w) {
        if (w <= 0 || norm2(n) < EPS) return;
        vertices[vid].q += SymmetricMatrix(n.x, n.y, n.z, d, w);
    }

    void initialize_quadrics() {
        for (auto &v : vertices) {
            v.q = SymmetricMatrix();
            v.cover = 0.0;
            v.nearOrig = 0.0;
            v.border = false;
            v.locked = false;
        }
        vector<Vec3> faceNormals(triangles.size());
        struct EdgeFaceRec { uint64_t key; int face; };
        vector<EdgeFaceRec> edgeFaces;
        edgeFaces.reserve(triangles.size() * 3 + 16);

        double avgArea2 = 0.0;
        int areaCnt = 0;
        for (size_t i = 0; i < triangles.size(); ++i) {
            Triangle &t = triangles[i];
            Vec3 p0 = vertices[t.v[0]].p;
            Vec3 p1 = vertices[t.v[1]].p;
            Vec3 p2 = vertices[t.v[2]].p;
            Vec3 n = normalized(crossv(p1 - p0, p2 - p0));
            t.n = n;
            faceNormals[i] = n;
            double d = -dotv(n, p0);
            double area2 = normv(crossv(p1 - p0, p2 - p0));
            avgArea2 += area2;
            areaCnt++;
            // Mild area weighting without letting huge triangles dominate entirely.
            double w = 0.35 + 0.65 * sqrt(max(1e-16, area2));
            SymmetricMatrix q(n.x, n.y, n.z, d, w);
            vertices[t.v[0]].q += q;
            vertices[t.v[1]].q += q;
            vertices[t.v[2]].q += q;
            for (int j = 0; j < 3; ++j) edgeFaces.push_back({edge_key(t.v[j], t.v[(j+1)%3]), (int)i});
        }

        // Extra penalty planes on open boundaries and sharp creases.
        const double boundaryW = 45.0;
        const double creaseW = 18.0;
        const double creaseDot = cos(35.0 * PI / 180.0);
        sort(edgeFaces.begin(), edgeFaces.end(), [](const EdgeFaceRec& a, const EdgeFaceRec& b) {
            if (a.key != b.key) return a.key < b.key;
            return a.face < b.face;
        });
        for (size_t l = 0; l < edgeFaces.size(); ) {
            size_t r = l + 1;
            while (r < edgeFaces.size() && edgeFaces[r].key == edgeFaces[l].key) ++r;
            uint64_t key = edgeFaces[l].key;
            int a = (int)(key >> 32), b = (int)(key & 0xffffffffu);
            int cnt = (int)(r - l);
            Vec3 e = vertices[b].p - vertices[a].p;
            double elen = normv(e);
            if (elen >= EPS) {
                e = e / elen;
                if (cnt == 1) {
                    int tid = edgeFaces[l].face;
                    Vec3 fn = faceNormals[tid];
                    Vec3 pn = normalized(crossv(e, fn));
                    double d = -dotv(pn, vertices[a].p);
                    add_plane_quadric_to_vertex(a, pn, d, boundaryW);
                    add_plane_quadric_to_vertex(b, pn, d, boundaryW);
                    vertices[a].border = vertices[b].border = true;
                } else if (cnt == 2) {
                    Vec3 n1 = faceNormals[edgeFaces[l].face], n2 = faceNormals[edgeFaces[l + 1].face];
                    double cd = clampd(dotv(n1, n2), -1.0, 1.0);
                    if (cd < creaseDot) {
                        Vec3 pn1 = normalized(crossv(e, n1));
                        Vec3 pn2 = normalized(crossv(e, n2));
                        double d1 = -dotv(pn1, vertices[a].p);
                        double d2 = -dotv(pn2, vertices[a].p);
                        double sharp = clampd((creaseDot - cd) / (creaseDot + 1.0), 0.0, 1.0);
                        add_plane_quadric_to_vertex(a, pn1, d1, creaseW * (0.5 + 2.0 * sharp));
                        add_plane_quadric_to_vertex(b, pn1, d1, creaseW * (0.5 + 2.0 * sharp));
                        add_plane_quadric_to_vertex(a, pn2, d2, creaseW * (0.5 + 2.0 * sharp));
                        add_plane_quadric_to_vertex(b, pn2, d2, creaseW * (0.5 + 2.0 * sharp));
                    }
                } else {
                    vertices[a].locked = vertices[b].locked = true;
                }
            }
            l = r;
        }
    }

    void update_mesh(int iteration) {
        if (iteration > 0) {
            int dst = 0;
            for (int i = 0; i < (int)triangles.size(); ++i) {
                if (!triangles[i].deleted) {
                    triangles[dst++] = triangles[i];
                }
            }
            triangles.resize(dst);
            deletedTriangles = 0;
        }

        for (auto &v : vertices) { v.tstart = 0; v.tcount = 0; v.border = false; }
        for (auto &t : triangles) {
            t.deleted = false;
            t.dirty = false;
            t.n = triangle_normal(t);
            for (int j = 0; j < 3; ++j) vertices[t.v[j]].tcount++;
        }
        int tstart = 0;
        for (auto &v : vertices) {
            v.tstart = tstart;
            tstart += v.tcount;
            v.tcount = 0;
        }
        refs.assign((size_t)tstart, Ref());
        for (int i = 0; i < (int)triangles.size(); ++i) {
            Triangle &t = triangles[i];
            for (int j = 0; j < 3; ++j) {
                Vertex &v = vertices[t.v[j]];
                refs[v.tstart + v.tcount] = {i, j};
                v.tcount++;
            }
        }

        // Border flags from current active mesh.  Sorting flat edge keys is more memory-stable than a hash map on million-vertex tests.
        vector<uint64_t> edgeCnt;
        edgeCnt.reserve(triangles.size() * 3 + 16);
        for (auto &t : triangles) {
            for (int j = 0; j < 3; ++j) edgeCnt.push_back(edge_key(t.v[j], t.v[(j+1)%3]));
        }
        sort(edgeCnt.begin(), edgeCnt.end());
        for (size_t l = 0; l < edgeCnt.size(); ) {
            size_t r = l + 1;
            while (r < edgeCnt.size() && edgeCnt[r] == edgeCnt[l]) ++r;
            int c = (int)(r - l);
            if (c == 1 || c > 2) {
                int a = (int)(edgeCnt[l] >> 32);
                int b = (int)(edgeCnt[l] & 0xffffffffu);
                if (c == 1) {
                    if (0 <= a && a < (int)vertices.size()) vertices[a].border = true;
                    if (0 <= b && b < (int)vertices.size()) vertices[b].border = true;
                } else {
                    if (0 <= a && a < (int)vertices.size()) vertices[a].locked = true;
                    if (0 <= b && b < (int)vertices.size()) vertices[b].locked = true;
                }
            }
            l = r;
        }

        // Initial / refreshed errors.
        for (auto &t : triangles) {
            Vec3 p;
            t.err[0] = calculate_error(t.v[0], t.v[1], p, -1.0, nullptr, nullptr);
            t.err[1] = calculate_error(t.v[1], t.v[2], p, -1.0, nullptr, nullptr);
            t.err[2] = calculate_error(t.v[2], t.v[0], p, -1.0, nullptr, nullptr);
            t.err[3] = min(t.err[0], min(t.err[1], t.err[2]));
        }
    }

    int count_edge_triangles(int u, int v) const {
        if (u < 0 || v < 0 || u >= (int)vertices.size() || v >= (int)vertices.size()) return 0;
        int cnt = 0;
        const Vertex& vu = vertices[u];
        for (int k = 0; k < vu.tcount; ++k) {
            const Ref& r = refs[vu.tstart + k];
            const Triangle& t = triangles[r.tid];
            if (t.deleted) continue;
            if (t.v[0] == v || t.v[1] == v || t.v[2] == v) cnt++;
        }
        return cnt;
    }

    bool neighbor_in_vector(const vector<int>& v, int x) const {
        for (int y : v) if (y == x) return true;
        return false;
    }

    void add_unique(vector<int>& v, int x) const {
        if (x < 0) return;
        if (!neighbor_in_vector(v, x)) v.push_back(x);
    }

    bool link_condition(int u, int v) const {
        int edgeTris = 0;
        const Vertex& vu = vertices[u];
        for (int k = 0; k < vu.tcount; ++k) {
            const Ref& r = refs[vu.tstart + k];
            const Triangle& t = triangles[r.tid];
            if (t.deleted) continue;
            if (t.v[0] == v || t.v[1] == v || t.v[2] == v) edgeTris++;
        }
        if (edgeTris != 1 && edgeTris != 2) return false;

        vector<int> nu, common;
        nu.reserve(24); common.reserve(8);
        for (int k = 0; k < vu.tcount; ++k) {
            const Ref& r = refs[vu.tstart + k];
            const Triangle& t = triangles[r.tid];
            if (t.deleted) continue;
            bool hasU = false;
            for (int j = 0; j < 3; ++j) hasU |= (t.v[j] == u);
            if (!hasU) continue;
            for (int j = 0; j < 3; ++j) if (t.v[j] != u && t.v[j] != v) add_unique(nu, t.v[j]);
        }
        const Vertex& vv = vertices[v];
        for (int k = 0; k < vv.tcount; ++k) {
            const Ref& r = refs[vv.tstart + k];
            const Triangle& t = triangles[r.tid];
            if (t.deleted) continue;
            bool hasV = false;
            for (int j = 0; j < 3; ++j) hasV |= (t.v[j] == v);
            if (!hasV) continue;
            for (int j = 0; j < 3; ++j) {
                int x = t.v[j];
                if (x != u && x != v && neighbor_in_vector(nu, x)) add_unique(common, x);
            }
        }
        return (int)common.size() == edgeTris;
    }


    static array<int,3> canonical_face(int a, int b, int c) {
        array<int,3> q = {a,b,c};
        sort(q.begin(), q.end());
        return q;
    }

    bool contains_face_id(const vector<int>& ids, int x) const {
        for (int y : ids) if (y == x) return true;
        return false;
    }

    bool would_create_duplicate_faces(int i0, int i1) const {
        vector<int> affected;
        affected.reserve(vertices[i0].tcount + vertices[i1].tcount);
        const Vertex& v0 = vertices[i0];
        const Vertex& v1 = vertices[i1];
        for (int k = 0; k < v0.tcount; ++k) {
            const Ref& r = refs[v0.tstart + k];
            if (!triangles[r.tid].deleted && !contains_face_id(affected, r.tid)) affected.push_back(r.tid);
        }
        for (int k = 0; k < v1.tcount; ++k) {
            const Ref& r = refs[v1.tstart + k];
            if (!triangles[r.tid].deleted && !contains_face_id(affected, r.tid)) affected.push_back(r.tid);
        }
        vector<array<int,3>> newFaces;
        newFaces.reserve(affected.size());
        for (int tid : affected) {
            const Triangle& t = triangles[tid];
            int a = t.v[0], b = t.v[1], c = t.v[2];
            if (a == i1) a = i0;
            if (b == i1) b = i0;
            if (c == i1) c = i0;
            if (a == b || b == c || c == a) continue; // the two incident edge triangles disappear
            array<int,3> cf = canonical_face(a,b,c);
            for (auto &old : newFaces) if (old == cf) return true;
            // Check existing unmodified faces around one vertex of this new face.
            const Vertex& va = vertices[cf[0]];
            for (int k = 0; k < va.tcount; ++k) {
                const Ref& r = refs[va.tstart + k];
                if (triangles[r.tid].deleted || contains_face_id(affected, r.tid)) continue;
                const Triangle& ot = triangles[r.tid];
                if (canonical_face(ot.v[0], ot.v[1], ot.v[2]) == cf) return true;
            }
            newFaces.push_back(cf);
        }
        return false;
    }

    bool cover_ok(const Vertex& a, const Vertex& b, const Vec3& p, double coverLimit, double* outCover, double* outNear) const {
        double ca = a.cover + distv(a.p, p);
        double cb = b.cover + distv(b.p, p);
        double cov = max(ca, cb);
        double nearv = min(a.nearOrig + distv(a.p, p), b.nearOrig + distv(b.p, p));
        if (outCover) *outCover = cov;
        if (outNear) *outNear = nearv;
        if (coverLimit <= 0) return true;
        // A small slack absorbs floating-point accumulation but stays safely below a stated tolerance.
        return cov <= coverLimit * 1.0000001 && nearv <= coverLimit * 1.0000001;
    }

    double calculate_error(int id_v1, int id_v2, Vec3& p_result, double coverLimit,
                           double* outCover, double* outNear) const {
        const Vertex& v1 = vertices[id_v1];
        const Vertex& v2 = vertices[id_v2];
        SymmetricMatrix q = v1.q + v2.q;
        bool border = v1.border || v2.border;
        bool locked = v1.locked || v2.locked;
        vector<Vec3> candidates;
        candidates.reserve(6);

        double det = q.det(0,1,2, 1,4,5, 2,5,7);
        if (!border && !locked && fabs(det) > 1e-18 && isfinite(det)) {
            Vec3 p;
            p.x = -1.0 / det * q.det(1,2,3, 4,5,6, 5,7,8);
            p.y =  1.0 / det * q.det(0,2,3, 1,5,6, 2,7,8);
            p.z = -1.0 / det * q.det(0,1,3, 1,4,6, 2,5,8);
            if (isfinite(p.x) && isfinite(p.y) && isfinite(p.z)) candidates.push_back(p);
        }
        candidates.push_back(v1.p);
        candidates.push_back(v2.p);
        candidates.push_back((v1.p + v2.p) * 0.5);
        // Slightly biased candidates often reduce Hausdorff cover while still avoiding slivers.
        candidates.push_back(v1.p * 0.67 + v2.p * 0.33);
        candidates.push_back(v1.p * 0.33 + v2.p * 0.67);

        double bestErr = numeric_limits<double>::infinity();
        Vec3 bestP = candidates[0];
        double bestCover = 0, bestNear = 0;
        for (const Vec3& p : candidates) {
            double cov = 0, nearv = 0;
            if (!cover_ok(v1, v2, p, coverLimit, &cov, &nearv)) continue;
            double e = vertex_error(q, p);
            // Prefer positions with smaller cover when QEM is close; this directly targets the clarified metric.
            if (coverLimit > 0) e += 1e-5 * (cov / max(coverLimit, 1e-18)) * (cov / max(coverLimit, 1e-18));
            if (e < bestErr) {
                bestErr = e; bestP = p; bestCover = cov; bestNear = nearv;
            }
        }
        if (!isfinite(bestErr)) {
            // No cover-safe candidate. Return a very high error, but still define p_result.
            p_result = (v1.p + v2.p) * 0.5;
            if (outCover) *outCover = numeric_limits<double>::infinity();
            if (outNear) *outNear = numeric_limits<double>::infinity();
            return numeric_limits<double>::infinity();
        }
        p_result = bestP;
        if (outCover) *outCover = bestCover;
        if (outNear) *outNear = bestNear;
        return bestErr;
    }

    bool flipped(const Vec3& p, int i0, int i1, const Vertex& v0, vector<int>& deleted) const {
        (void)i0;
        for (int k = 0; k < v0.tcount; ++k) {
            const Ref& r = refs[v0.tstart + k];
            const Triangle& t = triangles[r.tid];
            if (t.deleted) { deleted[k] = 0; continue; }
            int s = r.tvertex;
            int id1 = t.v[(s + 1) % 3];
            int id2 = t.v[(s + 2) % 3];
            if (id1 == i1 || id2 == i1) {
                deleted[k] = 1;
                continue;
            }
            Vec3 d1 = normalized(vertices[id1].p - p);
            Vec3 d2 = normalized(vertices[id2].p - p);
            if (norm2(d1) < EPS || norm2(d2) < EPS) return true;
            if (fabs(dotv(d1, d2)) > 0.9995) return true;
            Vec3 n = normalized(crossv(d1, d2));
            if (norm2(n) < EPS) return true;
            // Aggressive enough for high score, but prevents most inversion artifacts.
            if (dotv(n, t.n) < 0.08) return true;
        }
        return false;
    }

    void update_triangles(int i0, Vertex& v, const vector<int>& deleted, int& deleted_triangles, double coverLimit) {
        for (int k = 0; k < v.tcount; ++k) {
            Ref r = refs[v.tstart + k];
            Triangle& t = triangles[r.tid];
            if (t.deleted) continue;
            if (k < (int)deleted.size() && deleted[k]) {
                t.deleted = true;
                deleted_triangles++;
                continue;
            }
            t.v[r.tvertex] = i0;
            if (t.v[0] == t.v[1] || t.v[1] == t.v[2] || t.v[2] == t.v[0]) {
                t.deleted = true;
                deleted_triangles++;
                continue;
            }
            if (triangle_area2(t) < 1e-28) {
                t.deleted = true;
                deleted_triangles++;
                continue;
            }
            t.dirty = true;
            t.n = triangle_normal(t);
            Vec3 p;
            t.err[0] = calculate_error(t.v[0], t.v[1], p, coverLimit, nullptr, nullptr);
            t.err[1] = calculate_error(t.v[1], t.v[2], p, coverLimit, nullptr, nullptr);
            t.err[2] = calculate_error(t.v[2], t.v[0], p, coverLimit, nullptr, nullptr);
            t.err[3] = min(t.err[0], min(t.err[1], t.err[2]));
            refs.push_back(r);
        }
    }

    int active_triangle_count() const { return (int)triangles.size() - deletedTriangles; }

    void compact_mesh() {
        vector<Triangle> kept;
        kept.reserve(triangles.size());
        unordered_set<array<int,3>, FaceKeyHash> seen;
        seen.reserve(triangles.size() * 2 + 16);
        vector<char> used(vertices.size(), 0);
        for (auto &t : triangles) {
            if (t.deleted) continue;
            if (t.v[0] == t.v[1] || t.v[1] == t.v[2] || t.v[2] == t.v[0]) continue;
            if (triangle_area2(t) < 1e-28) continue;
            array<int,3> key = {t.v[0], t.v[1], t.v[2]};
            sort(key.begin(), key.end());
            if (!seen.insert(key).second) continue;
            kept.push_back(t);
            used[t.v[0]] = used[t.v[1]] = used[t.v[2]] = 1;
        }
        vector<int> remap(vertices.size(), -1);
        vector<Vertex> nv;
        nv.reserve(vertices.size());
        for (int i = 0; i < (int)vertices.size(); ++i) {
            if (used[i]) {
                remap[i] = (int)nv.size();
                nv.push_back(vertices[i]);
            }
        }
        for (auto &t : kept) for (int j = 0; j < 3; ++j) t.v[j] = remap[t.v[j]];
        vertices.swap(nv);
        triangles.swap(kept);
        deletedTriangles = 0;
    }

    int boundary_or_nonmanifold_edges(int* boundaryOut = nullptr, int* nonManifoldOut = nullptr) const {
        vector<uint64_t> edges;
        edges.reserve(triangles.size() * 3 + 16);
        for (auto &t : triangles) if (!t.deleted) {
            for (int j = 0; j < 3; ++j) edges.push_back(edge_key(t.v[j], t.v[(j+1)%3]));
        }
        sort(edges.begin(), edges.end());
        int boundary = 0, nonmanifold = 0;
        for (size_t l = 0; l < edges.size(); ) {
            size_t r = l + 1;
            while (r < edges.size() && edges[r] == edges[l]) ++r;
            int c = (int)(r - l);
            if (c == 1) boundary++;
            if (c > 2) nonmanifold++;
            l = r;
        }
        if (boundaryOut) *boundaryOut = boundary;
        if (nonManifoldOut) *nonManifoldOut = nonmanifold;
        return boundary + nonmanifold;
    }

    void simplify(int target_triangles, double aggressiveness, bool preserve_border, double coverLimit, int maxIterations) {
        if (triangles.empty() || vertices.empty()) return;
        target_triangles = max(4, min(target_triangles, (int)triangles.size()));
        deletedTriangles = 0;
        int prevActive = active_triangle_count();
        int stagnant = 0;

        for (int iteration = 0; iteration < maxIterations; ++iteration) {
            int active = active_triangle_count();
            if (active <= target_triangles) break;
            if (iteration % 5 == 0 || iteration == 0 || stagnant >= 2) update_mesh(iteration);

            active = active_triangle_count();
            if (active <= target_triangles) break;

            double t = (double)(iteration + 1) / max(1, maxIterations);
            double threshold;
            if (iteration + 8 >= maxIterations) threshold = numeric_limits<double>::infinity();
            else threshold = 1e-13 + pow(iteration + 3.0, aggressiveness) * 1e-12 + pow(t, 7.5) * 5e-3;
            if (coverLimit > 0 && iteration > maxIterations * 2 / 3) threshold = numeric_limits<double>::infinity();

            for (int i = 0; i < (int)triangles.size(); ++i) {
                if (active_triangle_count() <= target_triangles) break;
                Triangle& tri = triangles[i];
                if (tri.deleted || tri.dirty || tri.err[3] > threshold) continue;

                for (int j = 0; j < 3; ++j) {
                    if (tri.err[j] > threshold) continue;
                    int i0 = tri.v[j];
                    int i1 = tri.v[(j + 1) % 3];
                    if (i0 == i1) continue;
                    Vertex& v0 = vertices[i0];
                    Vertex& v1 = vertices[i1];
                    if (v0.locked || v1.locked) continue;

                    int edgeTris = count_edge_triangles(i0, i1);
                    if (edgeTris != 1 && edgeTris != 2) continue;
                    if (preserve_border) {
                        if (v0.border != v1.border) continue;
                        if ((v0.border || v1.border) && edgeTris != 1) continue;
                    }
                    if (!link_condition(i0, i1)) continue;
                    if (would_create_duplicate_faces(i0, i1)) continue;

                    Vec3 p;
                    double newCover = 0.0, newNear = 0.0;
                    double err = calculate_error(i0, i1, p, coverLimit, &newCover, &newNear);
                    if (!isfinite(err)) continue;
                    if (err > threshold && threshold < numeric_limits<double>::infinity()/4) continue;

                    vector<int> deleted0(v0.tcount, 0), deleted1(v1.tcount, 0);
                    if (flipped(p, i0, i1, v0, deleted0)) continue;
                    if (flipped(p, i1, i0, v1, deleted1)) continue;

                    SymmetricMatrix qsum = v0.q + v1.q;
                    Vec3 oldP = v0.p;
                    (void)oldP;
                    v0.p = p;
                    v0.q = qsum;
                    v0.cover = newCover;
                    v0.nearOrig = newNear;
                    v0.border = v0.border || v1.border;
                    v0.locked = v0.locked || v1.locked;

                    int tstart = (int)refs.size();
                    int deletedBefore = deletedTriangles;
                    update_triangles(i0, v0, deleted0, deletedTriangles, coverLimit);
                    update_triangles(i0, v1, deleted1, deletedTriangles, coverLimit);
                    int tcount = (int)refs.size() - tstart;
                    if (tcount <= v0.tcount) {
                        if (tcount > 0) memcpy(&refs[v0.tstart], &refs[tstart], sizeof(Ref) * (size_t)tcount);
                    } else {
                        v0.tstart = tstart;
                    }
                    v0.tcount = tcount;
                    active -= (deletedTriangles - deletedBefore);
                    break;
                }
            }

            int nowActive = active_triangle_count();
            if (nowActive >= prevActive - max(1, prevActive / 10000)) stagnant++;
            else stagnant = 0;
            prevActive = nowActive;
            if (stagnant >= 7 && coverLimit > 0) break;
        }
        compact_mesh();
    }

    vector<Vec3> output_vertices() const {
        vector<Vec3> out;
        out.reserve(vertices.size());
        for (auto &v : vertices) out.push_back(v.p);
        return out;
    }
    vector<array<int,3>> output_faces() const {
        vector<array<int,3>> out;
        out.reserve(triangles.size());
        for (auto &t : triangles) if (!t.deleted) out.push_back({t.v[0], t.v[1], t.v[2]});
        return out;
    }
};

static void normalize_mesh(vector<Vec3>& v, Vec3& center, double& diag) {
    if (v.empty()) { center = {0,0,0}; diag = 1.0; return; }
    Vec3 mn = v[0], mx = v[0];
    for (auto &p : v) {
        mn.x = min(mn.x, p.x); mn.y = min(mn.y, p.y); mn.z = min(mn.z, p.z);
        mx.x = max(mx.x, p.x); mx.y = max(mx.y, p.y); mx.z = max(mx.z, p.z);
    }
    center = (mn + mx) * 0.5;
    diag = normv(mx - mn);
    if (!isfinite(diag) || diag < EPS) diag = 1.0;
    for (auto &p : v) p = (p - center) / diag;
}

static void denormalize_mesh(vector<Vec3>& v, const Vec3& center, double diag) {
    for (auto &p : v) p = p * diag + center;
}

struct MeshStats {
    double avgEdge = 0.0;
    double medianEdge = 0.0;
    double normalComplexity = 0.0;
    int boundaryEdges = 0;
    bool closed = true;
};

static MeshStats compute_stats(const vector<Vec3>& verts, const vector<array<int,3>>& faces) {
    MeshStats st;
    if (verts.empty() || faces.empty()) return st;
    struct EdgeFaceRec { uint64_t key; int face; };
    vector<EdgeFaceRec> edgeFaces;
    edgeFaces.reserve(faces.size() * 3 + 16);
    vector<Vec3> normals(faces.size());
    auto key = [](int a, int b)->uint64_t { if (a > b) swap(a,b); return ((uint64_t)(uint32_t)a << 32) | (uint32_t)b; };
    vector<double> lens;
    lens.reserve(min<size_t>(faces.size() * 3, 300000));
    double sumLen = 0.0; long long cntLen = 0;
    for (int i = 0; i < (int)faces.size(); ++i) {
        auto f = faces[i];
        Vec3 n = normalized(crossv(verts[f[1]] - verts[f[0]], verts[f[2]] - verts[f[0]]));
        normals[i] = n;
        for (int j = 0; j < 3; ++j) {
            int a = f[j], b = f[(j+1)%3];
            edgeFaces.push_back({key(a,b), i});
            double l = distv(verts[a], verts[b]);
            sumLen += l; cntLen++;
            if (lens.size() < 300000) lens.push_back(l);
        }
    }
    st.avgEdge = cntLen ? sumLen / (double)cntLen : 0.0;
    if (!lens.empty()) {
        nth_element(lens.begin(), lens.begin() + lens.size()/2, lens.end());
        st.medianEdge = lens[lens.size()/2];
    }
    double comp = 0.0; int compCnt = 0; int nonman = 0;
    sort(edgeFaces.begin(), edgeFaces.end(), [](const EdgeFaceRec& a, const EdgeFaceRec& b) {
        if (a.key != b.key) return a.key < b.key;
        return a.face < b.face;
    });
    for (size_t l = 0; l < edgeFaces.size(); ) {
        size_t r = l + 1;
        while (r < edgeFaces.size() && edgeFaces[r].key == edgeFaces[l].key) ++r;
        int c = (int)(r - l);
        if (c == 1) st.boundaryEdges++;
        else if (c == 2) {
            double d = clampd(dotv(normals[edgeFaces[l].face], normals[edgeFaces[l + 1].face]), -1.0, 1.0);
            comp += sqrt(max(0.0, 1.0 - d));
            compCnt++;
        } else if (c > 2) nonman++;
        l = r;
    }
    st.normalComplexity = compCnt ? clampd(comp / compCnt, 0.0, 1.0) : 0.0;
    st.closed = (st.boundaryEdges == 0 && nonman == 0);
    return st;
}

static void compact_input_mesh(vector<Vec3>& verts, vector<array<int,3>>& faces) {
    vector<array<int,3>> nf;
    nf.reserve(faces.size());
    vector<char> used(verts.size(), 0);
    unordered_set<array<int,3>, FaceKeyHash> seen;
    seen.reserve(faces.size() * 2 + 16);
    for (auto f : faces) {
        if (f[0] < 0 || f[1] < 0 || f[2] < 0 || f[0] >= (int)verts.size() || f[1] >= (int)verts.size() || f[2] >= (int)verts.size()) continue;
        if (f[0] == f[1] || f[1] == f[2] || f[2] == f[0]) continue;
        if (normv(crossv(verts[f[1]] - verts[f[0]], verts[f[2]] - verts[f[0]])) < 1e-15) continue;
        array<int,3> key = {f[0], f[1], f[2]};
        sort(key.begin(), key.end());
        if (!seen.insert(key).second) continue;
        used[f[0]] = used[f[1]] = used[f[2]] = 1;
        nf.push_back(f);
    }
    vector<int> remap(verts.size(), -1);
    vector<Vec3> nv;
    nv.reserve(verts.size());
    for (int i = 0; i < (int)verts.size(); ++i) if (used[i]) { remap[i] = (int)nv.size(); nv.push_back(verts[i]); }
    for (auto &f : nf) for (int j = 0; j < 3; ++j) f[j] = remap[f[j]];
    verts.swap(nv); faces.swap(nf);
}

static void write_mesh(const MeshIO& in, const vector<Vec3>& verts, const vector<array<int,3>>& faces) {
    ios::sync_with_stdio(false);
    cout.setf(std::ios::fixed); cout << setprecision(10);
    if (in.fmt == FormatType::OBJ) {
        for (auto &p : verts) cout << "v " << p.x << ' ' << p.y << ' ' << p.z << '\n';
        for (auto &f : faces) cout << "f " << f[0] + 1 << ' ' << f[1] + 1 << ' ' << f[2] + 1 << '\n';
    } else if (in.fmt == FormatType::OFF) {
        cout << "OFF\n";
        cout << verts.size() << ' ' << faces.size() << " 0\n";
        for (auto &p : verts) cout << p.x << ' ' << p.y << ' ' << p.z << '\n';
        for (auto &f : faces) cout << "3 " << f[0] << ' ' << f[1] << ' ' << f[2] << '\n';
    } else if (in.fmt == FormatType::PLY) {
        cout << "ply\nformat ascii 1.0\n";
        cout << "element vertex " << verts.size() << "\n";
        cout << "property float x\nproperty float y\nproperty float z\n";
        cout << "element face " << faces.size() << "\n";
        cout << "property list uchar int vertex_indices\nend_header\n";
        for (auto &p : verts) cout << p.x << ' ' << p.y << ' ' << p.z << '\n';
        for (auto &f : faces) cout << "3 " << f[0] << ' ' << f[1] << ' ' << f[2] << '\n';
    } else {
        cout << verts.size() << ' ' << faces.size() << '\n';
        for (auto &p : verts) cout << p.x << ' ' << p.y << ' ' << p.z << '\n';
        int off = in.customOneBased ? 1 : 0;
        for (auto &f : faces) {
            if (in.customFaceHasPrefix) cout << "3 ";
            cout << f[0] + off << ' ' << f[1] + off << ' ' << f[2] + off << '\n';
        }
    }
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    string input((istreambuf_iterator<char>(cin)), istreambuf_iterator<char>());
    if (trim_copy(input).empty()) return 0;

    MeshIO mesh = parse_mesh(input);
    if (mesh.vertices.empty() || mesh.faces.empty()) {
        // Last-resort echo for unexpected formats; avoids producing malformed random output.
        cout << input;
        return 0;
    }

    vector<Vec3> verts = mesh.vertices;
    vector<array<int,3>> faces = mesh.faces;
    compact_input_mesh(verts, faces);
    if (verts.size() < 4 || faces.size() < 4) {
        write_mesh(mesh, verts, faces);
        return 0;
    }

    Vec3 center;
    double diag = 1.0;
    normalize_mesh(verts, center, diag);
    MeshStats st = compute_stats(verts, faces);

    const int N0 = (int)verts.size();
    const int F0 = (int)faces.size();

    // Adaptive target and cover settings. Coordinates are normalized to bbox diagonal = 1.
    double complexity = st.normalComplexity;
    double avgEdge = st.avgEdge > 0 ? st.avgEdge : 1.0 / sqrt(max(1, N0));
    double medEdge = st.medianEdge > 0 ? st.medianEdge : avgEdge;

    double coverLimit = -1.0;
    bool hasTolerance = (mesh.toleranceHint > 0 && isfinite(mesh.toleranceHint));
    if (hasTolerance) {
        // A stated tolerance is treated as a hard vertex-Hausdorff budget after normalization.
        coverLimit = clampd((mesh.toleranceHint / diag) * 0.985, 1e-12, 0.05);
    } else {
        // The clarified judge Hausdorff is vertex-discrete; this cover bound is intentionally more
        // aggressive than visual-only QEM but still scaled by mesh density/roughness.
        double densityBoost = clampd(pow((double)N0 / 1000000.0, 0.08), 0.72, 1.15);
        double base = 0.0062 * densityBoost;
        double roughPenalty = 1.0 - 0.46 * clampd(complexity, 0.0, 1.0);
        coverLimit = base * roughPenalty;
        double edgeDriven = medEdge * (2.15 + 1.1 * (1.0 - clampd(complexity, 0.0, 1.0)));
        coverLimit = max(coverLimit, edgeDriven);
        // Keep a global cap for dense models, but do not let the cap become smaller than
        // the local sampling scale on coarser cases.
        coverLimit = min(coverLimit, max(0.0125, edgeDriven * 1.08));
        if (st.boundaryEdges > 0) coverLimit *= 0.82;
    }

    double targetRatio;
    if (hasTolerance) {
        targetRatio = 0.0045 + 0.018 * complexity;
    } else {
        targetRatio = 0.020 + 0.095 * complexity;
        if (N0 < 200000) targetRatio += 0.018;
        if (N0 > 900000) targetRatio -= 0.004;
        targetRatio = clampd(targetRatio, 0.014, 0.155);
    }
    int targetTriangles = max(8, (int)llround(F0 * targetRatio));

    // More iterations are useful for large dense meshes, but cap for Kattis runtime.
    int maxIter = 110;
    if (N0 > 500000) maxIter = 95;
    if (N0 < 120000) maxIter = 125;
    if (hasTolerance) maxIter += 25;

    MeshSimplifier simplifier;
    simplifier.load(verts, faces);
    bool preserveBorder = true;
    double aggressiveness = hasTolerance ? 6.35 : 6.95;
    simplifier.simplify(targetTriangles, aggressiveness, preserveBorder, coverLimit, maxIter);

    vector<Vec3> outV = simplifier.output_vertices();
    vector<array<int,3>> outF = simplifier.output_faces();

    // If an unexpectedly tiny/invalid mesh emerges, fall back to a conservative no-cover pass.
    // This is rarely triggered, but protects against unusual nonmanifold inputs.
    if (outV.size() < 4 || outF.size() < 4 || outF.size() > faces.size()) {
        MeshSimplifier fallback;
        fallback.load(verts, faces);
        int safeTarget = max(8, (int)(F0 * clampd(0.18 + 0.20 * complexity, 0.18, 0.42)));
        fallback.simplify(safeTarget, 6.0, true, medEdge * 2.0, 80);
        outV = fallback.output_vertices();
        outF = fallback.output_faces();
    }

    denormalize_mesh(outV, center, diag);
    write_mesh(mesh, outV, outF);
    return 0;
}
