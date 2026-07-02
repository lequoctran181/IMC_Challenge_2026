#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <queue>
#include <string>
#include <vector>

using namespace std;

struct Vec3 {
    double x, y, z;
};

static inline Vec3 operator+(const Vec3& a, const Vec3& b) {
    return {a.x + b.x, a.y + b.y, a.z + b.z};
}
static inline Vec3 operator-(const Vec3& a, const Vec3& b) {
    return {a.x - b.x, a.y - b.y, a.z - b.z};
}
static inline Vec3 operator*(const Vec3& a, double s) {
    return {a.x * s, a.y * s, a.z * s};
}
static inline double dotp(const Vec3& a, const Vec3& b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}
static inline Vec3 crossp(const Vec3& a, const Vec3& b) {
    return {a.y * b.z - a.z * b.y,
            a.z * b.x - a.x * b.z,
            a.x * b.y - a.y * b.x};
}
static inline double norm2(const Vec3& a) {
    return dotp(a, a);
}
static inline double normv(const Vec3& a) {
    return sqrt(norm2(a));
}

struct Face {
    int a, b, c;
    unsigned char alive;
};

struct Quadric {
    double q00 = 0, q01 = 0, q02 = 0, q03 = 0;
    double q11 = 0, q12 = 0, q13 = 0;
    double q22 = 0, q23 = 0, q33 = 0;

    void addPlane(double a, double b, double c, double d, double w) {
        q00 += w * a * a;
        q01 += w * a * b;
        q02 += w * a * c;
        q03 += w * a * d;
        q11 += w * b * b;
        q12 += w * b * c;
        q13 += w * b * d;
        q22 += w * c * c;
        q23 += w * c * d;
        q33 += w * d * d;
    }

    void add(const Quadric& o) {
        q00 += o.q00; q01 += o.q01; q02 += o.q02; q03 += o.q03;
        q11 += o.q11; q12 += o.q12; q13 += o.q13;
        q22 += o.q22; q23 += o.q23; q33 += o.q33;
    }

    double eval(const Vec3& p) const {
        const double x = p.x, y = p.y, z = p.z;
        return q00 * x * x + 2.0 * q01 * x * y + 2.0 * q02 * x * z + 2.0 * q03 * x
             + q11 * y * y + 2.0 * q12 * y * z + 2.0 * q13 * y
             + q22 * z * z + 2.0 * q23 * z
             + q33;
    }
};

struct FastInput {
    vector<char> buf;
    char* p = nullptr;

    void readAll() {
        char chunk[1 << 16];
        size_t n;
        while ((n = fread(chunk, 1, sizeof(chunk), stdin)) > 0) {
            buf.insert(buf.end(), chunk, chunk + n);
        }
        buf.push_back('\0');
        p = buf.data();
    }

    inline void skipWs() {
        while (*p == ' ' || *p == '\n' || *p == '\r' || *p == '\t') ++p;
    }

    long readLong() {
        skipWs();
        return strtol(p, &p, 10);
    }

    double readDouble() {
        skipWs();
        return strtod(p, &p);
    }

    void skipLetter() {
        skipWs();
        ++p;
    }
};

static int N0, M0;
static vector<Vec3> P;
static vector<Vec3> OrigP;
static vector<Face> F;
static vector<vector<int>> incident;
static vector<Quadric> quadric;
static vector<Vec3> boxMinV;
static vector<Vec3> boxMaxV;
static vector<int> memberHead;
static vector<int> memberTail;
static vector<int> memberNext;
static vector<int> memberSize;
static vector<float> importance;
static vector<unsigned char> aliveV;
static vector<int> versionV;
static int aliveVertices = 0;
static int aliveFaces = 0;

static vector<int> markSeen;
static int markStamp = 1;

static double coverLimit = 0.0;
static double coverLimit2 = 0.0;
static double normalCosLimit = 0.65;

static inline unsigned long long edgeKey(int a, int b) {
    if (a > b) swap(a, b);
    return (static_cast<unsigned long long>(static_cast<unsigned int>(a)) << 32)
         | static_cast<unsigned int>(b);
}

static inline int keyA(unsigned long long k) {
    return static_cast<int>(k >> 32);
}

static inline int keyB(unsigned long long k) {
    return static_cast<int>(k & 0xffffffffu);
}

static inline bool faceContains(const Face& f, int v) {
    return f.a == v || f.b == v || f.c == v;
}

static inline bool faceContains2(const Face& f, int u, int v) {
    return faceContains(f, u) && faceContains(f, v);
}

static inline Vec3 minVec(const Vec3& a, const Vec3& b) {
    return {min(a.x, b.x), min(a.y, b.y), min(a.z, b.z)};
}

static inline Vec3 maxVec(const Vec3& a, const Vec3& b) {
    return {max(a.x, b.x), max(a.y, b.y), max(a.z, b.z)};
}

static bool mergedBoxFits(int keep, int other) {
    Vec3 mn = minVec(boxMinV[keep], boxMinV[other]);
    Vec3 mx = maxVec(boxMaxV[keep], boxMaxV[other]);
    const Vec3& p = P[keep];
    double dx = max(fabs(mn.x - p.x), fabs(mx.x - p.x));
    double dy = max(fabs(mn.y - p.y), fabs(mx.y - p.y));
    double dz = max(fabs(mn.z - p.z), fabs(mx.z - p.z));
    return dx * dx + dy * dy + dz * dz <= coverLimit2;
}

static bool clusterFits(int keep, int other) {
    if (mergedBoxFits(keep, other)) return true;
    const Vec3& p = P[keep];
    for (int v = memberHead[other]; v != -1; v = memberNext[v]) {
        if (norm2(OrigP[v] - p) > coverLimit2) return false;
    }
    return true;
}

static bool clusterFitsPoint(int root, const Vec3& p) {
    double dx = max(fabs(boxMinV[root].x - p.x), fabs(boxMaxV[root].x - p.x));
    double dy = max(fabs(boxMinV[root].y - p.y), fabs(boxMaxV[root].y - p.y));
    double dz = max(fabs(boxMinV[root].z - p.z), fabs(boxMaxV[root].z - p.z));
    if (dx * dx + dy * dy + dz * dz <= coverLimit2) return true;

    for (int v = memberHead[root]; v != -1; v = memberNext[v]) {
        if (norm2(OrigP[v] - p) > coverLimit2) return false;
    }
    return true;
}

static bool mergedClustersFitPoint(int a, int b, const Vec3& p) {
    return clusterFitsPoint(a, p) && clusterFitsPoint(b, p);
}

static bool faceNormalCurrent(int fid, Vec3& n, double& len) {
    const Face& f = F[fid];
    Vec3 cr = crossp(P[f.b] - P[f.a], P[f.c] - P[f.a]);
    len = normv(cr);
    if (len <= 1e-30) return false;
    n = cr * (1.0 / len);
    return true;
}

static bool faceNormalAfterReplace(int fid, int from, int to, Vec3& n, double& len) {
    Face f = F[fid];
    if (f.a == from) f.a = to;
    if (f.b == from) f.b = to;
    if (f.c == from) f.c = to;
    if (f.a == f.b || f.b == f.c || f.c == f.a) return false;
    Vec3 cr = crossp(P[f.b] - P[f.a], P[f.c] - P[f.a]);
    len = normv(cr);
    if (len <= 1e-30) return false;
    n = cr * (1.0 / len);
    return true;
}

static bool faceNormalAfterCollapsePoint(int fid, int from, int to, const Vec3& p, Vec3& n, double& len) {
    Face f = F[fid];
    if (f.a == from) f.a = to;
    if (f.b == from) f.b = to;
    if (f.c == from) f.c = to;
    if (f.a == f.b || f.b == f.c || f.c == f.a) return false;

    auto coord = [&](int v) -> Vec3 {
        return v == to ? p : P[v];
    };

    Vec3 cr = crossp(coord(f.b) - coord(f.a), coord(f.c) - coord(f.a));
    len = normv(cr);
    if (len <= 1e-30) return false;
    n = cr * (1.0 / len);
    return true;
}

static void compactIncident(int v) {
    if (!aliveV[v]) return;
    vector<int>& lst = incident[v];
    int w = 0;
    for (int fid : lst) {
        if (F[fid].alive && faceContains(F[fid], v)) {
            lst[w++] = fid;
        }
    }
    lst.resize(w);
}

static void collectNeighbors(int v, vector<int>& out) {
    out.clear();
    if (!aliveV[v]) return;
    if (++markStamp == 0) {
        fill(markSeen.begin(), markSeen.end(), 0);
        markStamp = 1;
    }
    int liveHits = 0;
    const int stamp = markStamp;
    vector<int>& lst = incident[v];
    for (int fid : lst) {
        if (!F[fid].alive) continue;
        const Face& f = F[fid];
        if (!faceContains(f, v)) continue;
        ++liveHits;
        const int vv[3] = {f.a, f.b, f.c};
        for (int i = 0; i < 3; ++i) {
            int w = vv[i];
            if (w == v || !aliveV[w]) continue;
            if (markSeen[w] != stamp) {
                markSeen[w] = stamp;
                out.push_back(w);
            }
        }
    }
    if ((int)lst.size() > liveHits * 4 + 64) compactIncident(v);
}

static int countEdgeFaces(int u, int v) {
    if (incident[u].size() > incident[v].size()) swap(u, v);
    int cnt = 0;
    for (int fid : incident[u]) {
        if (F[fid].alive && faceContains2(F[fid], u, v)) ++cnt;
    }
    return cnt;
}

struct Candidate {
    float cost;
    int from, to;
    int verFrom, verTo;
    double x, y, z;
};

struct CandidateGreater {
    bool operator()(const Candidate& a, const Candidate& b) const {
        return a.cost > b.cost;
    }
};

static priority_queue<Candidate, vector<Candidate>, CandidateGreater> pq;

static inline Vec3 candPos(const Candidate& c) {
    return {c.x, c.y, c.z};
}

static bool solveQuadricPoint(const Quadric& q, Vec3& out) {
    double a[3][4] = {
        {q.q00, q.q01, q.q02, -q.q03},
        {q.q01, q.q11, q.q12, -q.q13},
        {q.q02, q.q12, q.q22, -q.q23},
    };

    for (int col = 0; col < 3; ++col) {
        int pivot = col;
        for (int r = col + 1; r < 3; ++r) {
            if (fabs(a[r][col]) > fabs(a[pivot][col])) pivot = r;
        }
        if (fabs(a[pivot][col]) < 1e-14) return false;
        if (pivot != col) {
            for (int c = col; c < 4; ++c) swap(a[pivot][c], a[col][c]);
        }

        double inv = 1.0 / a[col][col];
        for (int c = col; c < 4; ++c) a[col][c] *= inv;
        for (int r = 0; r < 3; ++r) {
            if (r == col) continue;
            double factor = a[r][col];
            if (factor == 0.0) continue;
            for (int c = col; c < 4; ++c) a[r][c] -= factor * a[col][c];
        }
    }

    out = {a[0][3], a[1][3], a[2][3]};
    return isfinite(out.x) && isfinite(out.y) && isfinite(out.z);
}

static bool makeCandidate(int a, int b, Candidate& cand) {
    if (a == b || !aliveV[a] || !aliveV[b]) return false;
    Vec3 ab = P[b] - P[a];
    const double d2 = norm2(ab);

    Quadric q = quadric[a];
    q.add(quadric[b]);
    const double imp = 1.0 + 6.0 * max(importance[a], importance[b]);
    const double tie = 1e-10 * d2 * imp;
    double bestCost = 1e300;
    Vec3 bestP = P[a];

    auto consider = [&](const Vec3& p) {
        if (!mergedClustersFitPoint(a, b, p)) return;
        double c = q.eval(p) * imp + tie;
        if (c < bestCost) {
            bestCost = c;
            bestP = p;
        }
    };

    consider(P[a]);
    consider(P[b]);
    consider((P[a] + P[b]) * 0.5);

    Vec3 opt;
    if (d2 > 1e-24 && solveQuadricPoint(q, opt)) {
        double t = dotp(opt - P[a], ab) / d2;
        t = max(0.0, min(1.0, t));
        consider(P[a] + ab * t);
        consider(opt);
    }

    if (bestCost >= 1e250) return false;

    int bestFrom = a;
    int bestTo = b;
    if (incident[a].size() > incident[b].size()) {
        bestFrom = b;
        bestTo = a;
    }

    cand.cost = static_cast<float>(bestCost);
    cand.from = bestFrom;
    cand.to = bestTo;
    cand.verFrom = versionV[bestFrom];
    cand.verTo = versionV[bestTo];
    cand.x = bestP.x;
    cand.y = bestP.y;
    cand.z = bestP.z;
    return true;
}

static void pushCandidate(int a, int b) {
    Candidate cand;
    if (makeCandidate(a, b, cand)) pq.push(cand);
}

static bool collapseLegal(const Candidate& cand) {
    int from = cand.from;
    int to = cand.to;
    Vec3 p = candPos(cand);
    if (from == to || !aliveV[from] || !aliveV[to]) return false;

    if (!mergedClustersFitPoint(to, from, p)) return false;

    int edgeFaces = countEdgeFaces(from, to);
    if (edgeFaces != 2) return false;

    vector<int> nf, nt;
    collectNeighbors(from, nf);
    collectNeighbors(to, nt);

    if (++markStamp == 0) {
        fill(markSeen.begin(), markSeen.end(), 0);
        markStamp = 1;
    }
    const int stamp = markStamp;
    for (int v : nf) {
        if (v != to) markSeen[v] = stamp;
    }
    int common = 0;
    for (int v : nt) {
        if (v != from && markSeen[v] == stamp) ++common;
    }
    if (common != 2) return false;

    compactIncident(from);
    compactIncident(to);

    auto faceOk = [&](int fid) {
        Vec3 oldN, newN;
        double oldLen = 0, newLen = 0;
        if (!faceNormalCurrent(fid, oldN, oldLen)) return false;
        if (!faceNormalAfterCollapsePoint(fid, from, to, p, newN, newLen)) return false;
        if (dotp(oldN, newN) < normalCosLimit) return false;
        if (newLen < 1e-20) return false;
        return true;
    };

    for (int fid : incident[from]) {
        if (!F[fid].alive) continue;
        const Face& f = F[fid];
        if (!faceContains(f, from)) continue;
        if (faceContains(f, to)) continue;
        if (!faceOk(fid)) return false;
    }
    for (int fid : incident[to]) {
        if (!F[fid].alive) continue;
        const Face& f = F[fid];
        if (!faceContains(f, to)) continue;
        if (faceContains(f, from)) continue;
        if (!faceOk(fid)) return false;
    }
    return true;
}

static void doCollapse(const Candidate& cand) {
    int from = cand.from;
    int to = cand.to;
    Vec3 p = candPos(cand);

    compactIncident(from);
    compactIncident(to);

    P[to] = p;
    boxMinV[to] = minVec(boxMinV[to], boxMinV[from]);
    boxMaxV[to] = maxVec(boxMaxV[to], boxMaxV[from]);
    memberNext[memberTail[to]] = memberHead[from];
    memberTail[to] = memberTail[from];
    memberSize[to] += memberSize[from];
    memberHead[from] = memberTail[from] = -1;
    memberSize[from] = 0;
    importance[to] = max(importance[to], importance[from]);
    quadric[to].add(quadric[from]);

    for (int fid : incident[from]) {
        if (!F[fid].alive) continue;
        Face& f = F[fid];
        if (!faceContains(f, from)) continue;

        if (faceContains(f, to)) {
            f.alive = 0;
            --aliveFaces;
            continue;
        }

        if (f.a == from) f.a = to;
        if (f.b == from) f.b = to;
        if (f.c == from) f.c = to;
        incident[to].push_back(fid);
    }

    aliveV[from] = 0;
    --aliveVertices;
    ++versionV[from];
    ++versionV[to];
    vector<int>().swap(incident[from]);

    vector<int> neigh;
    collectNeighbors(to, neigh);
    for (int v : neigh) pushCandidate(to, v);
}

static int chooseTargetVertexCount(double avgImportance) {
    if (N0 <= 4) return N0;
    if (N0 <= 20) return max(4, N0 - 1);

    double ratio;
    if (N0 < 10000) ratio = 0.103;
    else if (N0 < 30000) ratio = 0.088;
    else if (N0 < 70000) ratio = 0.084;
    else if (N0 < 500000) ratio = 0.081;
    else ratio = 0.078;

    ratio += min(0.020, avgImportance * 0.10);
    ratio = max(0.078, min(0.120, ratio));

    int target = static_cast<int>(ceil(N0 * ratio));
    target = max(4, min(target, N0));
    return target;
}

static void loadMesh() {
    FastInput in;
    in.readAll();
    N0 = static_cast<int>(in.readLong());
    M0 = static_cast<int>(in.readLong());

    P.resize(N0);
    OrigP.resize(N0);
    F.resize(M0);
    incident.assign(N0, {});
    quadric.assign(N0, {});
    boxMinV.resize(N0);
    boxMaxV.resize(N0);
    memberHead.resize(N0);
    memberTail.resize(N0);
    memberNext.assign(N0, -1);
    memberSize.assign(N0, 1);
    importance.assign(N0, 0.0f);
    aliveV.assign(N0, 1);
    versionV.assign(N0, 0);
    markSeen.assign(N0, 0);

    Vec3 mn = {1e100, 1e100, 1e100};
    Vec3 mx = {-1e100, -1e100, -1e100};
    for (int i = 0; i < N0; ++i) {
        in.skipLetter();
        P[i].x = in.readDouble();
        P[i].y = in.readDouble();
        P[i].z = in.readDouble();
        OrigP[i] = P[i];
        mn.x = min(mn.x, P[i].x); mn.y = min(mn.y, P[i].y); mn.z = min(mn.z, P[i].z);
        mx.x = max(mx.x, P[i].x); mx.y = max(mx.y, P[i].y); mx.z = max(mx.z, P[i].z);
        boxMinV[i] = P[i];
        boxMaxV[i] = P[i];
        memberHead[i] = memberTail[i] = i;
    }

    for (int i = 0; i < M0; ++i) {
        in.skipLetter();
        int a = static_cast<int>(in.readLong()) - 1;
        int b = static_cast<int>(in.readLong()) - 1;
        int c = static_cast<int>(in.readLong()) - 1;
        F[i] = {a, b, c, 1};
        incident[a].push_back(i);
        incident[b].push_back(i);
        incident[c].push_back(i);
    }

    Vec3 diagv = mx - mn;
    double diag = normv(diagv);
    coverLimit = 0.0498 * diag;
    coverLimit2 = coverLimit * coverLimit;
    aliveVertices = N0;
    aliveFaces = M0;
}

static double buildQuadricsAndImportance() {
    vector<Vec3> normalSum(N0, {0.0, 0.0, 0.0});
    vector<int> normalCnt(N0, 0);

    for (int i = 0; i < M0; ++i) {
        const Face& f = F[i];
        Vec3 cr = crossp(P[f.b] - P[f.a], P[f.c] - P[f.a]);
        double len = normv(cr);
        if (len <= 1e-30) continue;
        Vec3 n = cr * (1.0 / len);
        double d = -dotp(n, P[f.a]);
        double w = max(0.5 * len, 1e-12);

        quadric[f.a].addPlane(n.x, n.y, n.z, d, w);
        quadric[f.b].addPlane(n.x, n.y, n.z, d, w);
        quadric[f.c].addPlane(n.x, n.y, n.z, d, w);

        normalSum[f.a] = normalSum[f.a] + n;
        normalSum[f.b] = normalSum[f.b] + n;
        normalSum[f.c] = normalSum[f.c] + n;
        ++normalCnt[f.a];
        ++normalCnt[f.b];
        ++normalCnt[f.c];
    }

    double avgImp = 0.0;
    for (int i = 0; i < N0; ++i) {
        if (normalCnt[i] == 0) continue;
        double meanLen = normv(normalSum[i]) / max(1, normalCnt[i]);
        double curv = 1.0 - max(0.0, min(1.0, meanLen));
        importance[i] = static_cast<float>(max(0.0, min(1.0, curv)));
        avgImp += importance[i];
    }
    return avgImp / max(1, N0);
}

static void initializeQueue() {
    vector<unsigned long long> edges;
    edges.reserve(static_cast<size_t>(M0) * 3);
    for (int i = 0; i < M0; ++i) {
        const Face& f = F[i];
        edges.push_back(edgeKey(f.a, f.b));
        edges.push_back(edgeKey(f.b, f.c));
        edges.push_back(edgeKey(f.c, f.a));
    }
    sort(edges.begin(), edges.end());
    edges.erase(unique(edges.begin(), edges.end()), edges.end());

    for (unsigned long long k : edges) {
        pushCandidate(keyA(k), keyB(k));
    }
}

static void simplify() {
    const auto start = chrono::steady_clock::now();
    double avgImp = buildQuadricsAndImportance();
    int targetVertices = chooseTargetVertexCount(avgImp);
    initializeQueue();

    long long iter = 0;
    while (aliveVertices > targetVertices && !pq.empty()) {
        Candidate cand = pq.top();
        pq.pop();
        ++iter;

        if (!aliveV[cand.from] || !aliveV[cand.to]) continue;
        if (versionV[cand.from] != cand.verFrom || versionV[cand.to] != cand.verTo) continue;

        if (collapseLegal(cand)) {
            doCollapse(cand);
        }

        if ((iter & 4095) == 0) {
            double elapsed = chrono::duration<double>(chrono::steady_clock::now() - start).count();
            if (elapsed > 19.25) break;
        }
    }
}

static void saveMesh() {
    vector<unsigned char> used(N0, 0);
    for (const Face& f : F) {
        if (!f.alive) continue;
        used[f.a] = used[f.b] = used[f.c] = 1;
    }

    vector<int> remap(N0, -1);
    int outV = 0;
    for (int i = 0; i < N0; ++i) {
        if (used[i]) remap[i] = outV++;
    }

    int outF = 0;
    for (const Face& f : F) {
        if (!f.alive) continue;
        if (remap[f.a] < 0 || remap[f.b] < 0 || remap[f.c] < 0) continue;
        if (f.a == f.b || f.b == f.c || f.c == f.a) continue;
        ++outF;
    }

    string out;
    out.reserve(1 << 20);
    char line[128];

    auto flushIfLarge = [&]() {
        if (out.size() >= (1 << 20)) {
            fwrite(out.data(), 1, out.size(), stdout);
            out.clear();
        }
    };

    int n = snprintf(line, sizeof(line), "%d %d\n", outV, outF);
    out.append(line, line + n);

    for (int i = 0; i < N0; ++i) {
        if (remap[i] < 0) continue;
        n = snprintf(line, sizeof(line), "v %.10g %.10g %.10g\n", P[i].x, P[i].y, P[i].z);
        out.append(line, line + n);
        flushIfLarge();
    }

    for (const Face& f : F) {
        if (!f.alive) continue;
        int a = remap[f.a], b = remap[f.b], c = remap[f.c];
        if (a < 0 || b < 0 || c < 0 || a == b || b == c || c == a) continue;
        n = snprintf(line, sizeof(line), "f %d %d %d\n", a + 1, b + 1, c + 1);
        out.append(line, line + n);
        flushIfLarge();
    }

    if (!out.empty()) fwrite(out.data(), 1, out.size(), stdout);
}

int main() {
    loadMesh();
    simplify();
    saveMesh();
    return 0;
}
