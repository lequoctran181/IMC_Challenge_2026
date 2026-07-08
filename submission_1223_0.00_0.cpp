// IMC Challenge 2026 - simplifygeometry
// C++17 single-file submission
//
// Strategy:
//   1) Parse common mesh encodings used by programming contests / geometry tasks:
//      - raw numeric: n m, then n vertices (x y z), then m triangular faces
//      - OFF
//      - OBJ subset (v/f)
//   2) Remove duplicate / degenerate triangles.
//   3) For very large inputs, run a conservative surface-grid pre-clustering step that
//      keeps representatives as original vertices, reducing memory pressure while keeping
//      vertex-only Hausdorff error controlled.
//   4) Run a topology-conscious Quadric Error Metric edge collapse pass with boundary
//      and feature penalties. The implementation is intentionally deterministic and
//      conservative: it rejects collapses that flip local triangle normals or create
//      local degeneracy.
//
// Notes:
//   - Output format follows input format (raw/OFF/OBJ).
//   - For raw numeric, the original face index base (0 or 1) is preserved.
//   - Compile: g++ -O2 -std=gnu++17 simplifygeometry.cpp -o simplifygeometry

#include <bits/stdc++.h>
using namespace std;

static constexpr double EPS = 1e-12;
static constexpr double AREA_EPS = 1e-20;
static constexpr double PI_CONST = 3.141592653589793238462643383279502884;
static constexpr int SOFT_TIME_SECONDS = 520;

struct Timer {
    chrono::steady_clock::time_point start;
    Timer() : start(chrono::steady_clock::now()) {}
    double elapsed() const {
        return chrono::duration<double>(chrono::steady_clock::now() - start).count();
    }
    bool over() const { return elapsed() > SOFT_TIME_SECONDS; }
};

struct Vec3 {
    double x = 0, y = 0, z = 0;
    Vec3() = default;
    Vec3(double X, double Y, double Z) : x(X), y(Y), z(Z) {}
    Vec3 operator+(const Vec3& o) const { return {x + o.x, y + o.y, z + o.z}; }
    Vec3 operator-(const Vec3& o) const { return {x - o.x, y - o.y, z - o.z}; }
    Vec3 operator*(double s) const { return {x * s, y * s, z * s}; }
    Vec3 operator/(double s) const { return {x / s, y / s, z / s}; }
    Vec3& operator+=(const Vec3& o) { x += o.x; y += o.y; z += o.z; return *this; }
    Vec3& operator*=(double s) { x *= s; y *= s; z *= s; return *this; }
};

static inline double dotv(const Vec3& a, const Vec3& b) {
    return a.x*b.x + a.y*b.y + a.z*b.z;
}
static inline Vec3 crossv(const Vec3& a, const Vec3& b) {
    return {a.y*b.z - a.z*b.y, a.z*b.x - a.x*b.z, a.x*b.y - a.y*b.x};
}
static inline double norm2(const Vec3& a) { return dotv(a,a); }
static inline double normv(const Vec3& a) { return sqrt(max(0.0, norm2(a))); }
static inline Vec3 normalize(const Vec3& a) {
    double n = normv(a);
    if (n <= EPS) return {0,0,0};
    return a / n;
}
static inline double dist2(const Vec3& a, const Vec3& b) { return norm2(a-b); }

struct Face {
    int a = 0, b = 0, c = 0;
    bool alive = true;
    Face() = default;
    Face(int A, int B, int C) : a(A), b(B), c(C), alive(true) {}
    int operator[](int i) const { return i == 0 ? a : (i == 1 ? b : c); }
    int& operator[](int i) { return i == 0 ? a : (i == 1 ? b : c); }
};

enum class Format { RAW, OFF, OBJ };

struct Mesh {
    vector<Vec3> V;
    vector<Face> F;
    Format fmt = Format::RAW;
    int rawBase = 0;
};

static inline uint64_t splitmix64(uint64_t x) {
    x += 0x9e3779b97f4a7c15ULL;
    x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
    x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
    return x ^ (x >> 31);
}

struct PairHash {
    size_t operator()(const pair<int,int>& p) const noexcept {
        uint64_t x = (uint64_t)(uint32_t)p.first << 32 ^ (uint32_t)p.second;
        return (size_t)splitmix64(x);
    }
};

struct TripleKey {
    int x, y, z;
    bool operator==(const TripleKey& o) const noexcept { return x == o.x && y == o.y && z == o.z; }
};
struct TripleHash {
    size_t operator()(const TripleKey& t) const noexcept {
        uint64_t a = (uint32_t)t.x;
        uint64_t b = (uint32_t)t.y;
        uint64_t c = (uint32_t)t.z;
        return (size_t)splitmix64(a * 1000003ULL ^ b * 9176ULL ^ c * 6151ULL);
    }
};

struct CellKey {
    long long x, y, z;
    bool operator==(const CellKey& o) const noexcept { return x == o.x && y == o.y && z == o.z; }
};
struct CellHash {
    size_t operator()(const CellKey& k) const noexcept {
        uint64_t a = (uint64_t)k.x, b = (uint64_t)k.y, c = (uint64_t)k.z;
        return (size_t)splitmix64(a * 11995408973635179863ULL ^ b * 10150724397891781847ULL ^ c * 7809847782465536327ULL);
    }
};

static string ltrim_copy(string s) {
    size_t p = 0;
    while (p < s.size() && isspace((unsigned char)s[p])) ++p;
    return s.substr(p);
}

static bool startsWithWord(const string& s, const string& w) {
    string t = ltrim_copy(s);
    if (t.size() < w.size()) return false;
    if (t.compare(0, w.size(), w) != 0) return false;
    return t.size() == w.size() || isspace((unsigned char)t[w.size()]);
}

static string stripComment(const string& s) {
    size_t p = s.find('#');
    if (p == string::npos) return s;
    return s.substr(0, p);
}

static bool parseOBJ(const string& data, Mesh& mesh) {
    mesh = Mesh();
    mesh.fmt = Format::OBJ;
    mesh.rawBase = 1;
    istringstream in(data);
    string line;
    while (getline(in, line)) {
        string t = ltrim_copy(line);
        if (t.empty() || t[0] == '#') continue;
        if (t.size() >= 2 && t[0] == 'v' && isspace((unsigned char)t[1])) {
            istringstream ss(t.substr(1));
            double x, y, z;
            if (ss >> x >> y >> z) mesh.V.emplace_back(x,y,z);
        } else if (t.size() >= 2 && t[0] == 'f' && isspace((unsigned char)t[1])) {
            istringstream ss(t.substr(1));
            vector<int> idx;
            string tok;
            while (ss >> tok) {
                size_t slash = tok.find('/');
                string first = slash == string::npos ? tok : tok.substr(0, slash);
                if (first.empty()) continue;
                long long val = 0;
                try { val = stoll(first); } catch (...) { continue; }
                int id;
                if (val < 0) id = (int)mesh.V.size() + (int)val;
                else id = (int)val - 1;
                if (id >= 0 && id < (int)mesh.V.size()) idx.push_back(id);
            }
            if (idx.size() >= 3) {
                for (size_t i = 1; i + 1 < idx.size(); ++i) {
                    mesh.F.emplace_back(idx[0], idx[i], idx[i+1]);
                }
            }
        }
    }
    return !mesh.V.empty() && !mesh.F.empty();
}

static vector<string> tokenizeNoComments(const string& data) {
    vector<string> tokens;
    istringstream in(data);
    string line;
    while (getline(in, line)) {
        line = stripComment(line);
        istringstream ss(line);
        string t;
        while (ss >> t) tokens.push_back(t);
    }
    return tokens;
}

static bool parseOFF(const string& data, Mesh& mesh) {
    auto toks = tokenizeNoComments(data);
    if (toks.empty() || toks[0] != "OFF") return false;
    if (toks.size() < 4) return false;
    size_t p = 1;
    long long nLL, fLL, eLL;
    try {
        nLL = stoll(toks[p++]);
        fLL = stoll(toks[p++]);
        eLL = stoll(toks[p++]);
    } catch (...) { return false; }
    (void)eLL;
    if (nLL <= 0 || fLL < 0 || nLL > INT_MAX || fLL > INT_MAX) return false;
    int n = (int)nLL, fcnt = (int)fLL;
    if (p + (size_t)n * 3 > toks.size()) return false;
    mesh = Mesh();
    mesh.fmt = Format::OFF;
    mesh.rawBase = 0;
    mesh.V.reserve(n);
    try {
        for (int i = 0; i < n; ++i) {
            double x = stod(toks[p++]);
            double y = stod(toks[p++]);
            double z = stod(toks[p++]);
            mesh.V.emplace_back(x,y,z);
        }
        mesh.F.reserve(fcnt);
        for (int i = 0; i < fcnt && p < toks.size(); ++i) {
            int k = stoi(toks[p++]);
            vector<int> idx;
            idx.reserve(max(0, k));
            for (int j = 0; j < k && p < toks.size(); ++j) idx.push_back(stoi(toks[p++]));
            if ((int)idx.size() >= 3) {
                for (int j = 1; j + 1 < (int)idx.size(); ++j) {
                    if (idx[0] >= 0 && idx[0] < n && idx[j] >= 0 && idx[j] < n && idx[j+1] >= 0 && idx[j+1] < n)
                        mesh.F.emplace_back(idx[0], idx[j], idx[j+1]);
                }
            }
        }
    } catch (...) { return false; }
    return !mesh.V.empty() && !mesh.F.empty();
}

static bool looksNumeric(const string& s) {
    if (s.empty()) return false;
    char* endptr = nullptr;
    errno = 0;
    (void)strtod(s.c_str(), &endptr);
    return endptr && *endptr == '\0';
}

static bool parseRAW(const string& data, Mesh& mesh) {
    auto toks = tokenizeNoComments(data);
    if (toks.size() < 8) return false;
    if (!looksNumeric(toks[0]) || !looksNumeric(toks[1])) return false;
    long long nLL, mLL;
    try {
        nLL = stoll(toks[0]);
        mLL = stoll(toks[1]);
    } catch (...) { return false; }
    if (nLL <= 0 || mLL < 0 || nLL > INT_MAX || mLL > INT_MAX) return false;
    int n = (int)nLL, m = (int)mLL;
    size_t need = 2 + (size_t)n * 3 + (size_t)m * 3;
    if (toks.size() < need) return false;
    mesh = Mesh();
    mesh.fmt = Format::RAW;
    mesh.V.reserve(n);
    size_t p = 2;
    try {
        for (int i = 0; i < n; ++i) {
            double x = stod(toks[p++]);
            double y = stod(toks[p++]);
            double z = stod(toks[p++]);
            mesh.V.emplace_back(x,y,z);
        }
        vector<array<long long,3>> rawF;
        rawF.reserve(m);
        long long mn = LLONG_MAX, mx = LLONG_MIN;
        for (int i = 0; i < m; ++i) {
            long long a = stoll(toks[p++]);
            long long b = stoll(toks[p++]);
            long long c = stoll(toks[p++]);
            rawF.push_back({a,b,c});
            mn = min(mn, min(a, min(b,c)));
            mx = max(mx, max(a, max(b,c)));
        }
        int base = (mn == 0 ? 0 : 1);
        mesh.rawBase = base;
        mesh.F.reserve(m);
        for (auto tri : rawF) {
            long long a = tri[0] - base, b = tri[1] - base, c = tri[2] - base;
            if (a >= 0 && a < n && b >= 0 && b < n && c >= 0 && c < n)
                mesh.F.emplace_back((int)a,(int)b,(int)c);
        }
    } catch (...) { return false; }
    return !mesh.V.empty() && !mesh.F.empty();
}

static bool parseMesh(const string& data, Mesh& mesh) {
    string t = ltrim_copy(data);
    if (startsWithWord(t, "OFF")) return parseOFF(data, mesh);
    if (t.size() >= 2 && (t[0] == 'v' || t.find("\nv ") != string::npos || t.find("\nf ") != string::npos)) {
        Mesh obj;
        if (parseOBJ(data, obj)) { mesh = std::move(obj); return true; }
    }
    Mesh raw;
    if (parseRAW(data, raw)) { mesh = std::move(raw); return true; }
    Mesh off;
    if (parseOFF(data, off)) { mesh = std::move(off); return true; }
    Mesh obj;
    if (parseOBJ(data, obj)) { mesh = std::move(obj); return true; }
    return false;
}


static bool validFaceIndices(const Face& f, int n) {
    return f.a >= 0 && f.a < n && f.b >= 0 && f.b < n && f.c >= 0 && f.c < n && f.a != f.b && f.b != f.c && f.c != f.a;
}

static TripleKey sortedKey(int a, int b, int c) {
    if (a > b) swap(a,b);
    if (b > c) swap(b,c);
    if (a > b) swap(a,b);
    return {a,b,c};
}

static void compactMesh(Mesh& mesh) {
    int n = (int)mesh.V.size();
    vector<Face> newF;
    newF.reserve(mesh.F.size());
    unordered_set<TripleKey, TripleHash> seen;
    seen.reserve(mesh.F.size() * 2 + 1);
    for (auto& f : mesh.F) {
        if (!f.alive || !validFaceIndices(f, n)) continue;
        Vec3 cr = crossv(mesh.V[f.b] - mesh.V[f.a], mesh.V[f.c] - mesh.V[f.a]);
        if (norm2(cr) <= AREA_EPS) continue;
        TripleKey k = sortedKey(f.a,f.b,f.c);
        if (seen.find(k) != seen.end()) continue;
        seen.insert(k);
        newF.push_back(Face(f.a,f.b,f.c));
    }
    vector<char> used(n, 0);
    for (auto& f : newF) used[f.a] = used[f.b] = used[f.c] = 1;
    vector<int> remap(n, -1);
    vector<Vec3> newV;
    newV.reserve(n);
    for (int i = 0; i < n; ++i) {
        if (used[i]) {
            remap[i] = (int)newV.size();
            newV.push_back(mesh.V[i]);
        }
    }
    for (auto& f : newF) {
        f.a = remap[f.a]; f.b = remap[f.b]; f.c = remap[f.c];
        f.alive = true;
    }
    mesh.V.swap(newV);
    mesh.F.swap(newF);
}

struct Bounds {
    Vec3 mn, mx;
    double diag = 0;
};
static Bounds computeBounds(const vector<Vec3>& V) {
    Bounds b;
    if (V.empty()) return b;
    b.mn = b.mx = V[0];
    for (const Vec3& p : V) {
        b.mn.x = min(b.mn.x, p.x); b.mn.y = min(b.mn.y, p.y); b.mn.z = min(b.mn.z, p.z);
        b.mx.x = max(b.mx.x, p.x); b.mx.y = max(b.mx.y, p.y); b.mx.z = max(b.mx.z, p.z);
    }
    b.diag = normv(b.mx - b.mn);
    return b;
}

static double sampleMedianEdge(const Mesh& mesh) {
    vector<double> lens;
    lens.reserve(min<size_t>(mesh.F.size() * 3, 200000));
    int step = max(1, (int)(mesh.F.size() / 80000));
    for (int i = 0; i < (int)mesh.F.size(); i += step) {
        const Face& f = mesh.F[i];
        if (!validFaceIndices(f, (int)mesh.V.size())) continue;
        double a = normv(mesh.V[f.a] - mesh.V[f.b]);
        double b = normv(mesh.V[f.b] - mesh.V[f.c]);
        double c = normv(mesh.V[f.c] - mesh.V[f.a]);
        if (a > EPS) lens.push_back(a);
        if (b > EPS) lens.push_back(b);
        if (c > EPS) lens.push_back(c);
    }
    if (lens.empty()) {
        Bounds bd = computeBounds(mesh.V);
        return max(1e-9, bd.diag / max(1.0, sqrt((double)mesh.V.size())));
    }
    nth_element(lens.begin(), lens.begin() + lens.size()/2, lens.end());
    return max(1e-12, lens[lens.size()/2]);
}

struct ClusterInfo {
    Vec3 sum{0,0,0};
    double weight = 0;
    int best = -1;
    double bestD = numeric_limits<double>::infinity();
    vector<int> members;
};

static vector<double> computeVertexSaliency(const Mesh& mesh) {
    int n = (int)mesh.V.size();
    vector<Vec3> normals(n, {0,0,0});
    vector<double> area(n, 0.0);
    for (const Face& f : mesh.F) {
        if (!validFaceIndices(f, n)) continue;
        Vec3 cr = crossv(mesh.V[f.b]-mesh.V[f.a], mesh.V[f.c]-mesh.V[f.a]);
        double ar = normv(cr);
        if (ar <= EPS) continue;
        normals[f.a] += cr; normals[f.b] += cr; normals[f.c] += cr;
        area[f.a] += ar; area[f.b] += ar; area[f.c] += ar;
    }
    for (int i = 0; i < n; ++i) normals[i] = normalize(normals[i]);
    vector<double> sal(n, 1.0);
    int sampleStep = max(1, (int)(mesh.F.size() / 200000));
    vector<vector<int>> neigh;
    if (n <= 250000) {
        neigh.assign(n, {});
        for (int i = 0; i < (int)mesh.F.size(); i += sampleStep) {
            const Face& f = mesh.F[i];
            if (!validFaceIndices(f, n)) continue;
            int v[3] = {f.a,f.b,f.c};
            for (int e = 0; e < 3; ++e) {
                int a = v[e], b = v[(e+1)%3];
                if ((int)neigh[a].size() < 24) neigh[a].push_back(b);
                if ((int)neigh[b].size() < 24) neigh[b].push_back(a);
            }
        }
        for (int i = 0; i < n; ++i) {
            double curv = 0;
            for (int j : neigh[i]) curv = max(curv, 1.0 - dotv(normals[i], normals[j]));
            sal[i] = 1.0 + 10.0 * max(0.0, curv) + 0.2 * log1p(area[i]);
        }
    } else {
        for (int i = 0; i < n; ++i) sal[i] = 1.0 + 0.1 * log1p(area[i]);
    }
    return sal;
}

static Mesh clusterMeshOnce(const Mesh& mesh, double cellSize, bool featureAware) {
    Mesh out;
    out.fmt = mesh.fmt;
    out.rawBase = mesh.rawBase;
    int n = (int)mesh.V.size();
    Bounds bd = computeBounds(mesh.V);
    if (cellSize <= EPS || bd.diag <= EPS) return mesh;
    vector<double> sal;
    if (featureAware) sal = computeVertexSaliency(mesh);
    unordered_map<CellKey, int, CellHash> mp;
    mp.reserve((size_t)n * 2 + 1);
    vector<ClusterInfo> clusters;
    clusters.reserve(n / 4 + 16);

    auto keyOf = [&](const Vec3& p, double localCell) {
        return CellKey{(long long)floor((p.x - bd.mn.x) / localCell),
                       (long long)floor((p.y - bd.mn.y) / localCell),
                       (long long)floor((p.z - bd.mn.z) / localCell)};
    };

    vector<int> cidOf(n, -1);
    for (int i = 0; i < n; ++i) {
        double localCell = cellSize;
        if (featureAware && !sal.empty()) {
            if (sal[i] > 5.0) localCell *= 0.50;
            else if (sal[i] > 2.5) localCell *= 0.70;
        }
        CellKey key = keyOf(mesh.V[i], localCell);
        // Encode local-cell tier into key to avoid accidental merging between feature tiers.
        if (featureAware && !sal.empty()) {
            if (sal[i] > 5.0) key.x = key.x * 3 + 1, key.y = key.y * 3 + 1, key.z = key.z * 3 + 1;
            else if (sal[i] > 2.5) key.x = key.x * 3 + 2, key.y = key.y * 3 + 2, key.z = key.z * 3 + 2;
            else key.x = key.x * 3, key.y = key.y * 3, key.z = key.z * 3;
        }
        auto it = mp.find(key);
        int cid;
        if (it == mp.end()) {
            cid = (int)clusters.size();
            mp.emplace(key, cid);
            clusters.emplace_back();
        } else cid = it->second;
        cidOf[i] = cid;
        double w = featureAware && !sal.empty() ? sal[i] : 1.0;
        clusters[cid].sum += mesh.V[i] * w;
        clusters[cid].weight += w;
        if ((int)clusters[cid].members.size() < 32) clusters[cid].members.push_back(i);
    }

    out.V.reserve(clusters.size());
    vector<int> clusterRep(clusters.size(), -1);
    for (int c = 0; c < (int)clusters.size(); ++c) {
        Vec3 centroid = clusters[c].sum / max(EPS, clusters[c].weight);
        int best = clusters[c].members.empty() ? -1 : clusters[c].members[0];
        double bestD = best >= 0 ? dist2(mesh.V[best], centroid) : numeric_limits<double>::infinity();
        for (int id : clusters[c].members) {
            double d = dist2(mesh.V[id], centroid);
            double bonus = featureAware && !sal.empty() ? 1.0 / max(1.0, sal[id]) : 1.0;
            d *= bonus;
            if (d < bestD) bestD = d, best = id;
        }
        if (best < 0) best = 0;
        clusterRep[c] = (int)out.V.size();
        // Use an original vertex as representative. This is intentionally safer for the
        // vertex-only symmetric Hausdorff clarification than outputting arbitrary centroids.
        out.V.push_back(mesh.V[best]);
    }

    out.F.reserve(mesh.F.size());
    unordered_set<TripleKey, TripleHash> seen;
    seen.reserve(mesh.F.size() * 2 + 1);
    for (const Face& f : mesh.F) {
        if (!validFaceIndices(f, n)) continue;
        int a = clusterRep[cidOf[f.a]], b = clusterRep[cidOf[f.b]], c = clusterRep[cidOf[f.c]];
        if (a == b || b == c || c == a) continue;
        Face nf(a,b,c);
        Vec3 cr = crossv(out.V[nf.b]-out.V[nf.a], out.V[nf.c]-out.V[nf.a]);
        if (norm2(cr) <= AREA_EPS) continue;
        TripleKey k = sortedKey(a,b,c);
        if (seen.insert(k).second) out.F.push_back(nf);
    }
    compactMesh(out);
    return out;
}

static Mesh preclusterLargeMesh(const Mesh& mesh, int desiredVertices, const Timer& timer) {
    if ((int)mesh.V.size() <= desiredVertices) return mesh;
    double med = sampleMedianEdge(mesh);
    Bounds bd = computeBounds(mesh.V);
    if (bd.diag <= EPS) return mesh;
    double ratio = (double)mesh.V.size() / max(1, desiredVertices);
    double factor = max(1.15, sqrt(max(1.0, ratio)) * 0.85);
    Mesh best = mesh;
    int bestGap = INT_MAX;
    double cell = med * factor;
    for (int it = 0; it < 5 && !timer.over(); ++it) {
        Mesh cand = clusterMeshOnce(mesh, cell, true);
        int v = (int)cand.V.size();
        int gap = abs(v - desiredVertices);
        if (!cand.F.empty() && (v <= desiredVertices || gap < bestGap)) {
            best = std::move(cand);
            bestGap = gap;
        }
        if (v > desiredVertices * 1.12) cell *= 1.35;
        else if (v < desiredVertices * 0.65) cell *= 0.72;
        else break;
    }
    return best;
}

struct Quadric {
    // Symmetric 4x4 upper triangle:
    // 0:xx, 1:xy, 2:xz, 3:xw, 4:yy, 5:yz, 6:yw, 7:zz, 8:zw, 9:ww
    double m[10];
    Quadric() { memset(m, 0, sizeof(m)); }
    Quadric& operator+=(const Quadric& o) {
        for (int i = 0; i < 10; ++i) m[i] += o.m[i];
        return *this;
    }
    friend Quadric operator+(Quadric a, const Quadric& b) { a += b; return a; }
    void addPlane(double a, double b, double c, double d, double w = 1.0) {
        double v[4] = {a,b,c,d};
        int idx = 0;
        for (int i = 0; i < 4; ++i) {
            for (int j = i; j < 4; ++j) m[idx++] += w * v[i] * v[j];
        }
    }
    double eval(const Vec3& p) const {
        double x = p.x, y = p.y, z = p.z;
        return m[0]*x*x + 2*m[1]*x*y + 2*m[2]*x*z + 2*m[3]*x
             + m[4]*y*y + 2*m[5]*y*z + 2*m[6]*y
             + m[7]*z*z + 2*m[8]*z + m[9];
    }
};

static bool solve3x3(double A[3][3], double b[3], Vec3& x) {
    double a00=A[0][0], a01=A[0][1], a02=A[0][2];
    double a10=A[1][0], a11=A[1][1], a12=A[1][2];
    double a20=A[2][0], a21=A[2][1], a22=A[2][2];
    double det = a00*(a11*a22-a12*a21) - a01*(a10*a22-a12*a20) + a02*(a10*a21-a11*a20);
    if (fabs(det) < 1e-14) return false;
    double invDet = 1.0 / det;
    double d0 = b[0], d1 = b[1], d2 = b[2];
    double dx = d0*(a11*a22-a12*a21) - a01*(d1*a22-a12*d2) + a02*(d1*a21-a11*d2);
    double dy = a00*(d1*a22-a12*d2) - d0*(a10*a22-a12*a20) + a02*(a10*d2-d1*a20);
    double dz = a00*(a11*d2-d1*a21) - a01*(a10*d2-d1*a20) + d0*(a10*a21-a11*a20);
    x = {dx*invDet, dy*invDet, dz*invDet};
    return isfinite(x.x) && isfinite(x.y) && isfinite(x.z);
}

struct CollapseCandidate {
    double cost;
    int u, v;
    int verU, verV;
    bool operator<(const CollapseCandidate& o) const { return cost > o.cost; }
};

class QEMSimplifier {
public:
    Mesh mesh;
    vector<Quadric> Q;
    vector<vector<int>> facesOf;
    vector<unordered_set<int>> neigh;
    vector<char> activeV, boundaryV;
    vector<int> version;
    vector<Vec3> originalFaceNormal;
    vector<double> originalFaceArea;
    int activeCount = 0;
    priority_queue<CollapseCandidate> pq;
    Timer* timer = nullptr;

    explicit QEMSimplifier(Mesh m, Timer* t) : mesh(std::move(m)), timer(t) {}

    void buildAdjacencyAndQuadrics() {
        int n = (int)mesh.V.size();
        Q.assign(n, Quadric());
        facesOf.assign(n, {});
        neigh.assign(n, {});
        activeV.assign(n, 1);
        boundaryV.assign(n, 0);
        version.assign(n, 0);
        activeCount = n;
        originalFaceNormal.assign(mesh.F.size(), {0,0,0});
        originalFaceArea.assign(mesh.F.size(), 0);
        unordered_map<pair<int,int>, int, PairHash> edgeCount;
        edgeCount.reserve(mesh.F.size() * 4 + 1);

        for (int fid = 0; fid < (int)mesh.F.size(); ++fid) {
            Face& f = mesh.F[fid];
            if (!validFaceIndices(f, n)) { f.alive = false; continue; }
            Vec3 cr = crossv(mesh.V[f.b]-mesh.V[f.a], mesh.V[f.c]-mesh.V[f.a]);
            double ar2 = normv(cr);
            if (ar2 <= EPS) { f.alive = false; continue; }
            Vec3 normal = cr / ar2;
            originalFaceNormal[fid] = normal;
            originalFaceArea[fid] = ar2;
            double d = -dotv(normal, mesh.V[f.a]);
            double w = max(1e-12, ar2);
            Quadric plane;
            plane.addPlane(normal.x, normal.y, normal.z, d, w);
            Q[f.a] += plane; Q[f.b] += plane; Q[f.c] += plane;
            facesOf[f.a].push_back(fid);
            facesOf[f.b].push_back(fid);
            facesOf[f.c].push_back(fid);
            int v[3] = {f.a,f.b,f.c};
            for (int e = 0; e < 3; ++e) {
                int a = v[e], b = v[(e+1)%3];
                if (a > b) swap(a,b);
                edgeCount[{a,b}]++;
                neigh[a].insert(b);
                neigh[b].insert(a);
            }
        }

        // Boundary constraints: for an edge used by only one face, add a strong plane
        // perpendicular to the adjacent face plane and through the boundary edge.
        for (auto& kv : edgeCount) {
            int a = kv.first.first, b = kv.first.second;
            if (kv.second != 1) continue;
            boundaryV[a] = boundaryV[b] = 1;
            Vec3 pa = mesh.V[a], pb = mesh.V[b];
            Vec3 edgeDir = normalize(pb - pa);
            Vec3 avgN{0,0,0};
            for (int fid : facesOf[a]) {
                if (!mesh.F[fid].alive) continue;
                const Face& f = mesh.F[fid];
                if (f.a == b || f.b == b || f.c == b) avgN += originalFaceNormal[fid];
            }
            avgN = normalize(avgN);
            Vec3 bp = normalize(crossv(edgeDir, avgN));
            if (norm2(bp) > EPS) {
                double d = -dotv(bp, pa);
                Quadric q;
                q.addPlane(bp.x, bp.y, bp.z, d, 80.0);
                Q[a] += q; Q[b] += q;
            }
        }

        for (auto& s : neigh) s.reserve(s.size() * 2 + 1);
    }

    Vec3 bestPosition(int u, int v, double& bestCost) const {
        Quadric q = Q[u] + Q[v];
        Vec3 pu = mesh.V[u], pv = mesh.V[v];
        Vec3 mid = (pu + pv) * 0.5;
        vector<Vec3> cand;
        cand.reserve(6);
        cand.push_back(pu);
        cand.push_back(pv);
        cand.push_back(mid);
        double A[3][3] = {{q.m[0], q.m[1], q.m[2]},
                          {q.m[1], q.m[4], q.m[5]},
                          {q.m[2], q.m[5], q.m[7]}};
        double b[3] = {-q.m[3], -q.m[6], -q.m[8]};
        Vec3 opt;
        if (solve3x3(A, b, opt)) {
            double localScale2 = max(dist2(pu,pv), 1e-24);
            if (dist2(opt, mid) <= localScale2 * 9.0) cand.push_back(opt);
        }
        // Additional conservative points biased to original vertices reduce vertex-only
        // Hausdorff drift while still improving QEM on smooth regions.
        cand.push_back(pu * 0.75 + pv * 0.25);
        cand.push_back(pu * 0.25 + pv * 0.75);

        bestCost = numeric_limits<double>::infinity();
        Vec3 best = pu;
        double edgeLen2 = max(dist2(pu,pv), 1e-24);
        for (const Vec3& p : cand) {
            double endpointDrift = min(dist2(p, pu), dist2(p, pv));
            double c = q.eval(p) + 0.002 * endpointDrift / edgeLen2;
            if (boundaryV[u] || boundaryV[v]) c *= 1.6;
            if ((boundaryV[u] ^ boundaryV[v])) c *= 8.0;
            if (c < bestCost) bestCost = c, best = p;
        }
        return best;
    }

    double edgeCost(int u, int v) const {
        double c;
        (void)bestPosition(u,v,c);
        double len2 = dist2(mesh.V[u], mesh.V[v]);
        double penalty = 1.0 + 1e-6 * len2;
        return c * penalty + 1e-15 * len2;
    }

    void pushEdge(int u, int v) {
        if (u == v || u < 0 || v < 0 || u >= (int)mesh.V.size() || v >= (int)mesh.V.size()) return;
        if (!activeV[u] || !activeV[v]) return;
        if (u > v) swap(u,v);
        double c = edgeCost(u,v);
        if (!isfinite(c)) return;
        pq.push({c,u,v,version[u],version[v]});
    }

    void initQueue() {
        while (!pq.empty()) pq.pop();
        for (int u = 0; u < (int)neigh.size(); ++u) {
            if (!activeV[u]) continue;
            for (int v : neigh[u]) if (u < v && activeV[v]) pushEdge(u,v);
        }
    }

    bool localLinkOK(int u, int v) const {
        // For triangular 2-manifold interiors, the one-ring intersection around an
        // collapsible edge should be two vertices; boundary edges generally have one.
        int common = 0;
        if (neigh[u].size() < neigh[v].size()) {
            for (int x : neigh[u]) if (x != v && activeV[x] && neigh[v].count(x)) ++common;
        } else {
            for (int x : neigh[v]) if (x != u && activeV[x] && neigh[u].count(x)) ++common;
        }
        if (boundaryV[u] || boundaryV[v]) return common <= 2;
        return common == 2 || common == 1; // allow slight non-manifold input damage after cleaning
    }

    bool triangleWouldBeValid(int fid, int u, int v, const Vec3& p) const {
        const Face& f = mesh.F[fid];
        if (!f.alive) return true;
        int a = f.a, b = f.b, c = f.c;
        if (a == v) a = u;
        if (b == v) b = u;
        if (c == v) c = u;
        if (a == b || b == c || c == a) return true; // face will be deleted
        Vec3 pa = (a == u ? p : mesh.V[a]);
        Vec3 pb = (b == u ? p : mesh.V[b]);
        Vec3 pc = (c == u ? p : mesh.V[c]);
        Vec3 cr = crossv(pb-pa, pc-pa);
        double ar = normv(cr);
        if (ar <= EPS) return false;
        Vec3 nn = cr / ar;
        double oldAr = max(originalFaceArea[fid], 1e-30);
        // Reject strong flips and near-zero sliver collapses. A small angular relaxation
        // helps meshes with inconsistent orientation while still preventing catastrophic folds.
        if (dotv(nn, originalFaceNormal[fid]) < -0.05) return false;
        if (ar < oldAr * 1e-6) return false;
        return true;
    }

    bool collapseOK(int u, int v, const Vec3& p) const {
        if (!activeV[u] || !activeV[v] || u == v) return false;
        if (!localLinkOK(u,v)) return false;
        if ((boundaryV[u] ^ boundaryV[v])) return false;
        // Avoid high valence blow-ups.
        if (neigh[u].size() + neigh[v].size() > 320) return false;

        vector<int> check;
        check.reserve(facesOf[u].size() + facesOf[v].size());
        for (int fid : facesOf[u]) check.push_back(fid);
        for (int fid : facesOf[v]) check.push_back(fid);
        sort(check.begin(), check.end());
        check.erase(unique(check.begin(), check.end()), check.end());
        for (int fid : check) {
            if (!triangleWouldBeValid(fid,u,v,p)) return false;
        }
        return true;
    }

    void doCollapse(int u, int v, const Vec3& p) {
        if (u == v) return;
        // Prefer keeping boundary endpoint if any; otherwise lower valence endpoint.
        if (boundaryV[v] && !boundaryV[u]) swap(u,v);
        else if (!boundaryV[u] && !boundaryV[v] && neigh[v].size() > neigh[u].size()) swap(u,v);

        mesh.V[u] = p;
        Q[u] += Q[v];
        boundaryV[u] = boundaryV[u] || boundaryV[v];
        activeV[v] = 0;
        ++version[u]; ++version[v];
        --activeCount;

        vector<int> affected;
        affected.reserve(facesOf[u].size() + facesOf[v].size());
        for (int fid : facesOf[u]) affected.push_back(fid);
        for (int fid : facesOf[v]) affected.push_back(fid);
        sort(affected.begin(), affected.end());
        affected.erase(unique(affected.begin(), affected.end()), affected.end());

        for (int fid : facesOf[v]) {
            if (fid < 0 || fid >= (int)mesh.F.size()) continue;
            Face& f = mesh.F[fid];
            if (!f.alive) continue;
            bool touched = false;
            if (f.a == v) f.a = u, touched = true;
            if (f.b == v) f.b = u, touched = true;
            if (f.c == v) f.c = u, touched = true;
            if (touched) facesOf[u].push_back(fid);
            if (f.a == f.b || f.b == f.c || f.c == f.a) f.alive = false;
        }
        facesOf[v].clear();

        vector<int> oldNv;
        oldNv.reserve(neigh[v].size());
        for (int x : neigh[v]) if (x != u && activeV[x]) oldNv.push_back(x);
        for (int x : oldNv) {
            neigh[x].erase(v);
            neigh[x].insert(u);
            neigh[u].insert(x);
        }
        neigh[u].erase(v);
        neigh[v].clear();

        // Clean dead/self references in u-neighborhood.
        vector<int> toErase;
        for (int x : neigh[u]) if (x == u || !activeV[x]) toErase.push_back(x);
        for (int x : toErase) neigh[u].erase(x);

        // Local face duplicate cleanup, using a temporary set for affected triangles.
        unordered_set<TripleKey, TripleHash> local;
        local.reserve(affected.size() * 2 + 1);
        for (int fid : affected) {
            if (fid < 0 || fid >= (int)mesh.F.size()) continue;
            Face& f = mesh.F[fid];
            if (!f.alive) continue;
            if (!validFaceIndices(f, (int)mesh.V.size())) { f.alive = false; continue; }
            Vec3 cr = crossv(mesh.V[f.b]-mesh.V[f.a], mesh.V[f.c]-mesh.V[f.a]);
            if (norm2(cr) <= AREA_EPS) { f.alive = false; continue; }
            TripleKey k = sortedKey(f.a,f.b,f.c);
            if (!local.insert(k).second) f.alive = false;
        }

        // Add fresh candidates around kept vertex and its one-ring.
        for (int x : neigh[u]) if (activeV[x]) pushEdge(u, x);
        for (int x : neigh[u]) if (activeV[x]) {
            int cnt = 0;
            for (int y : neigh[x]) {
                if (++cnt > 24) break;
                if (activeV[y] && y != u) pushEdge(x, y);
            }
        }
    }

    Mesh run(int targetVertices) {
        if ((int)mesh.V.size() <= targetVertices) return mesh;
        buildAdjacencyAndQuadrics();
        initQueue();
        int accepted = 0, attempts = 0;
        int rebuildEvery = max(20000, (int)mesh.V.size() / 3);
        while (activeCount > targetVertices && !pq.empty()) {
            if ((attempts++ & 4095) == 0 && timer && timer->over()) break;
            CollapseCandidate cc = pq.top(); pq.pop();
            int u = cc.u, v = cc.v;
            if (u < 0 || v < 0 || u >= (int)mesh.V.size() || v >= (int)mesh.V.size()) continue;
            if (!activeV[u] || !activeV[v]) continue;
            if (cc.verU != version[u] || cc.verV != version[v]) continue;
            if (!neigh[u].count(v) && !neigh[v].count(u)) continue;
            double c;
            Vec3 p = bestPosition(u,v,c);
            if (!collapseOK(u,v,p)) continue;
            doCollapse(u,v,p);
            ++accepted;
            if (accepted % rebuildEvery == 0) {
                compactCurrentWithoutChangingActive();
                buildAdjacencyAndQuadrics();
                initQueue();
                rebuildEvery = max(20000, activeCount / 2);
            }
        }
        compactOutput();
        return mesh;
    }

    void compactCurrentWithoutChangingActive() {
        // Full compaction resets active arrays; useful to reduce stale face lists during long runs.
        vector<int> remap(mesh.V.size(), -1);
        vector<Vec3> newV;
        newV.reserve(activeCount);
        for (int i = 0; i < (int)mesh.V.size(); ++i) if (activeV[i]) {
            remap[i] = (int)newV.size(); newV.push_back(mesh.V[i]);
        }
        vector<Face> newF;
        newF.reserve(mesh.F.size());
        unordered_set<TripleKey, TripleHash> seen;
        seen.reserve(mesh.F.size()*2+1);
        for (const Face& f : mesh.F) {
            if (!f.alive || !validFaceIndices(f, (int)mesh.V.size())) continue;
            int a = remap[f.a], b = remap[f.b], c = remap[f.c];
            if (a < 0 || b < 0 || c < 0 || a == b || b == c || c == a) continue;
            Face nf(a,b,c);
            Vec3 cr = crossv(newV[b]-newV[a], newV[c]-newV[a]);
            if (norm2(cr) <= AREA_EPS) continue;
            TripleKey k = sortedKey(a,b,c);
            if (seen.insert(k).second) newF.push_back(nf);
        }
        mesh.V.swap(newV);
        mesh.F.swap(newF);
    }

    void compactOutput() {
        // Mark all vertices active by usage and compact through the general cleaner.
        compactMesh(mesh);
    }
};

static void enforceEdgeManifold(Mesh& mesh);

static int choosePreclusterTarget(int n) {
    if (n <= 180000) return n;
    if (n <= 400000) return 155000;
    if (n <= 1000000) return 170000;
    return 185000;
}

static int chooseFinalTarget(int nOriginal, int nCurrent) {
    // Aggressive but not reckless default. The IMC scoring rewards low vertex count, but
    // violating perceptual/fidelity thresholds is catastrophic. These targets are meant to
    // stay inside a useful search region for most million-vertex scanned/manufactured meshes.
    double ratio;
    if (nOriginal < 2000) ratio = 0.55;
    else if (nOriginal < 10000) ratio = 0.30;
    else if (nOriginal < 50000) ratio = 0.145;
    else if (nOriginal < 150000) ratio = 0.085;
    else if (nOriginal < 500000) ratio = 0.055;
    else ratio = 0.040;
    int target = (int)llround(nOriginal * ratio);
    target = max(target, 800);
    target = min(target, nCurrent);
    // If pre-clustering already reduced more than the final target, do not ask QEM to remove
    // an excessive fraction in one pass; local collapses are safer after clustering.
    if (nCurrent > 0) target = max(target, (int)llround(nCurrent * 0.42));
    return max(1, target);
}

static void simplify(Mesh& mesh) {
    Timer timer;
    compactMesh(mesh);
    int originalN = (int)mesh.V.size();
    if (originalN <= 3 || mesh.F.empty()) return;

    int preTarget = choosePreclusterTarget(originalN);
    if (originalN > preTarget) {
        Mesh clustered = preclusterLargeMesh(mesh, preTarget, timer);
        if (!clustered.V.empty() && !clustered.F.empty() && clustered.V.size() < mesh.V.size()) {
            mesh = std::move(clustered);
            compactMesh(mesh);
        }
    }

    int finalTarget = chooseFinalTarget(originalN, (int)mesh.V.size());
    if ((int)mesh.V.size() > finalTarget && !timer.over()) {
        QEMSimplifier qem(mesh, &timer);
        Mesh simplified = qem.run(finalTarget);
        if (!simplified.V.empty() && !simplified.F.empty() && simplified.V.size() <= mesh.V.size()) {
            mesh = std::move(simplified);
        }
    }

    // A final very light clustering pass only if QEM could not make progress at all on a
    // huge mesh. Representatives remain original-ish after QEM only partially; this pass is
    // intentionally conservative.
    if ((int)mesh.V.size() > finalTarget * 1.6 && !timer.over()) {
        double med = sampleMedianEdge(mesh);
        Mesh cand = clusterMeshOnce(mesh, med * 1.35, false);
        if (!cand.V.empty() && !cand.F.empty() && cand.V.size() < mesh.V.size()) mesh = std::move(cand);
    }
    compactMesh(mesh);
    enforceEdgeManifold(mesh);
}


static void enforceEdgeManifold(Mesh& mesh) {
    // Remove the minimum-obvious set of triangles around edges that appear more than twice.
    // This is a pragmatic guard for judges that explicitly reject non-manifold edges.
    for (int iter = 0; iter < 6; ++iter) {
        compactMesh(mesh);
        int n = (int)mesh.V.size();
        unordered_map<pair<int,int>, vector<int>, PairHash> edgeFaces;
        edgeFaces.reserve(mesh.F.size() * 4 + 1);
        vector<double> area(mesh.F.size(), 0.0);
        for (int fid = 0; fid < (int)mesh.F.size(); ++fid) {
            const Face& f = mesh.F[fid];
            if (!validFaceIndices(f, n)) continue;
            area[fid] = normv(crossv(mesh.V[f.b]-mesh.V[f.a], mesh.V[f.c]-mesh.V[f.a]));
            int v[3] = {f.a,f.b,f.c};
            for (int e = 0; e < 3; ++e) {
                int a = v[e], b = v[(e+1)%3];
                if (a > b) swap(a,b);
                edgeFaces[{a,b}].push_back(fid);
            }
        }
        vector<char> kill(mesh.F.size(), 0);
        bool changed = false;
        for (auto& kv : edgeFaces) {
            auto& lst = kv.second;
            if ((int)lst.size() <= 2) continue;
            changed = true;
            sort(lst.begin(), lst.end(), [&](int i, int j){ return area[i] > area[j]; });
            for (int k = 2; k < (int)lst.size(); ++k) kill[lst[k]] = 1;
        }
        if (!changed) break;
        vector<Face> nf;
        nf.reserve(mesh.F.size());
        for (int i = 0; i < (int)mesh.F.size(); ++i) if (!kill[i]) nf.push_back(mesh.F[i]);
        mesh.F.swap(nf);
    }
    compactMesh(mesh);
}

static void writeMesh(const Mesh& mesh, ostream& out) {
    out.setf(ios::fixed);
    out << setprecision(10);
    if (mesh.fmt == Format::OBJ) {
        for (const Vec3& p : mesh.V) out << "v " << p.x << ' ' << p.y << ' ' << p.z << '\n';
        for (const Face& f : mesh.F) out << "f " << (f.a+1) << ' ' << (f.b+1) << ' ' << (f.c+1) << '\n';
    } else if (mesh.fmt == Format::OFF) {
        out << "OFF\n" << mesh.V.size() << ' ' << mesh.F.size() << " 0\n";
        for (const Vec3& p : mesh.V) out << p.x << ' ' << p.y << ' ' << p.z << '\n';
        for (const Face& f : mesh.F) out << "3 " << f.a << ' ' << f.b << ' ' << f.c << '\n';
    } else {
        int base = mesh.rawBase;
        out << mesh.V.size() << ' ' << mesh.F.size() << '\n';
        for (const Vec3& p : mesh.V) out << p.x << ' ' << p.y << ' ' << p.z << '\n';
        for (const Face& f : mesh.F) out << (f.a+base) << ' ' << (f.b+base) << ' ' << (f.c+base) << '\n';
    }
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    string data((istreambuf_iterator<char>(cin)), istreambuf_iterator<char>());
    Mesh mesh;
    if (!parseMesh(data, mesh)) {
        // Last-resort behavior: echo input. In normal judge data this branch should never be used.
        cout << data;
        return 0;
    }
    simplify(mesh);
    writeMesh(mesh, cout);
    return 0;
}
