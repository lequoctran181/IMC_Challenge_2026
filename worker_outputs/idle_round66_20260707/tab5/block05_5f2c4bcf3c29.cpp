#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <iostream>
#include <limits>
#include <queue>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>
using namespace std;

struct Vec3 {
    double x, y, z;
};

static inline Vec3 operator+(const Vec3& a, const Vec3& b) { return {a.x + b.x, a.y + b.y, a.z + b.z}; }
static inline Vec3 operator-(const Vec3& a, const Vec3& b) { return {a.x - b.x, a.y - b.y, a.z - b.z}; }
static inline Vec3 operator*(const Vec3& a, double s) { return {a.x * s, a.y * s, a.z * s}; }
static inline Vec3 operator/(const Vec3& a, double s) { return {a.x / s, a.y / s, a.z / s}; }
static inline double dot3(const Vec3& a, const Vec3& b) { return a.x*b.x + a.y*b.y + a.z*b.z; }
static inline Vec3 cross3(const Vec3& a, const Vec3& b) {
    return {a.y*b.z - a.z*b.y, a.z*b.x - a.x*b.z, a.x*b.y - a.y*b.x};
}
static inline double norm2(const Vec3& a) { return dot3(a,a); }
static inline double norm3(const Vec3& a) { return sqrt(max(0.0, norm2(a))); }
static inline Vec3 normalized(const Vec3& a) {
    double n = norm3(a);
    if (!(n > 0)) return {0,0,0};
    return a / n;
}
static inline double coord(const Vec3& p, int axis) {
    return axis==0 ? p.x : (axis==1 ? p.y : p.z);
}
static inline pair<double,double> screen2(const Vec3& p, int axis) {
    if (axis == 0) return {p.y, p.z};
    if (axis == 1) return {p.x, p.z};
    return {p.x, p.y};
}

struct Face { int a, b, c; };
struct CandidateEval {
    bool ok=false;
    int v=-1;
    double maxRedistribute=0;
    vector<pair<int,int>> redistribute;
    double cost=1e100;
    vector<int> ring;
    vector<Face> add;
};

static int N, M0;
static vector<Vec3> P;
static vector<Face> F;
static vector<unsigned char> aliveV, aliveF;
static vector<vector<int>> VF;
static vector<unsigned> vstamp;
static vector<vector<int>> ownedOriginals;
static vector<int> faceMark;
static int markToken = 1;
static int liveV = 0, liveF = 0;
static double diagLen = 1.0, epsArea = 1e-20, radiusLimit = 1.0;

static inline uint64_t edgeKey(int a, int b) {
    if (a > b) swap(a,b);
    return (uint64_t)(uint32_t)a << 32 | (uint32_t)b;
}
static inline array<int,3> faceKey(int a, int b, int c) {
    array<int,3> r{a,b,c};
    sort(r.begin(), r.end());
    return r;
}
static inline bool hasVertex(const Face& f, int v) { return f.a==v || f.b==v || f.c==v; }
static inline Vec3 faceCross(const Face& f) { return cross3(P[f.b]-P[f.a], P[f.c]-P[f.a]); }
static inline double faceArea2(const Face& f) { return norm3(faceCross(f)); }

static void rebuildVF() {
    VF.assign(N, {});
    for (int i=0; i<(int)F.size(); ++i) if (aliveF[i]) {
        VF[F[i].a].push_back(i);
        VF[F[i].b].push_back(i);
        VF[F[i].c].push_back(i);
    }
}

static vector<int> incidentFaces(int v) {
    vector<int> out;
    if (v < 0 || v >= N || !aliveV[v]) return out;
    out.reserve(VF[v].size());
    for (int id: VF[v]) {
        if (id >= 0 && id < (int)F.size() && aliveF[id] && hasVertex(F[id], v)) out.push_back(id);
    }
    sort(out.begin(), out.end());
    out.erase(unique(out.begin(), out.end()), out.end());
    return out;
}

static void markFaces(const vector<int>& ids) {
    ++markToken;
    if (markToken > 2000000000) {
        fill(faceMark.begin(), faceMark.end(), 0);
        markToken = 1;
    }
    for (int id: ids) if (id >= 0 && id < (int)faceMark.size()) faceMark[id] = markToken;
}

static int edgeCountOutside(int a, int b) {
    int cnt = 0;
    const vector<int>& L = (VF[a].size() < VF[b].size()) ? VF[a] : VF[b];
    for (int id: L) {
        if (id < 0 || id >= (int)F.size() || !aliveF[id] || faceMark[id] == markToken) continue;
        if (hasVertex(F[id], a) && hasVertex(F[id], b)) ++cnt;
    }
    return cnt;
}

static bool duplicateFaceOutside(int a, int b, int c) {
    auto key = faceKey(a,b,c);
    int best = a;
    if (VF[b].size() < VF[best].size()) best = b;
    if (VF[c].size() < VF[best].size()) best = c;
    for (int id: VF[best]) {
        if (id < 0 || id >= (int)F.size() || !aliveF[id] || faceMark[id] == markToken) continue;
        const Face& g = F[id];
        if (faceKey(g.a,g.b,g.c) == key) return true;
    }
    return false;
}

static int indexInSmall(const vector<int>& a, int x) {
    for (int i=0; i<(int)a.size(); ++i) if (a[i] == x) return i;
    return -1;
}

static bool buildLinkCycle(int v, const vector<int>& inc, vector<int>& ring) {
    vector<int> verts;
    vector<pair<int,int>> edges;
    verts.reserve(inc.size()*2);
    edges.reserve(inc.size());
    for (int id: inc) {
        Face f = F[id];
        int arr[3] = {f.a, f.b, f.c};
        int o[2], k=0;
        for (int i=0; i<3; ++i) {
            if (arr[i] != v) {
                if (k >= 2) return false;
                o[k++] = arr[i];
            }
        }
        if (k != 2 || o[0] == o[1] || !aliveV[o[0]] || !aliveV[o[1]]) return false;
        edges.push_back({o[0], o[1]});
        verts.push_back(o[0]); verts.push_back(o[1]);
    }
    sort(verts.begin(), verts.end());
    verts.erase(unique(verts.begin(), verts.end()), verts.end());
    int k = (int)verts.size();
    if (k != (int)inc.size()) return false;
    vector<array<int,2>> adj(k, {-1,-1});
    vector<int> deg(k,0);
    for (auto &e: edges) {
        int a = indexInSmall(verts, e.first), b = indexInSmall(verts, e.second);
        if (a < 0 || b < 0 || a == b) return false;
        if (deg[a] >= 2 || deg[b] >= 2) return false;
        adj[a][deg[a]++] = b;
        adj[b][deg[b]++] = a;
    }
    for (int i=0; i<k; ++i) if (deg[i] != 2) return false;
    ring.clear(); ring.reserve(k);
    vector<unsigned char> seen(k,0);
    int start = 0, prev = -1, cur = start;
    for (int step=0; step<k; ++step) {
        if (seen[cur]) return false;
        seen[cur] = 1;
        ring.push_back(verts[cur]);
        int nxt = (adj[cur][0] == prev ? adj[cur][1] : adj[cur][0]);
        prev = cur; cur = nxt;
    }
    return cur == start;
}

struct Basis2 {
    Vec3 n, u, w;
};
static Basis2 makeBasis(const Vec3& normal) {
    Basis2 b;
    b.n = normalized(normal);
    Vec3 a = (fabs(b.n.x) < 0.75) ? Vec3{1,0,0} : Vec3{0,1,0};
    b.u = normalized(cross3(b.n, a));
    if (norm3(b.u) <= 0) b.u = {0,0,1};
    b.w = cross3(b.u, b.n);
    return b;
}
static inline pair<double,double> projectBasis(const Vec3& p, const Basis2& b) {
    return {dot3(p, b.u), dot3(p, b.w)};
}
static inline double orient2(pair<double,double> a, pair<double,double> b, pair<double,double> c) {
    return (b.first-a.first)*(c.second-a.second) - (b.second-a.second)*(c.first-a.first);
}
static bool insideTri2(pair<double,double> p, pair<double,double> a, pair<double,double> b, pair<double,double> c, double sgn) {
    double o1 = orient2(a,b,p) * sgn;
    double o2 = orient2(b,c,p) * sgn;
    double o3 = orient2(c,a,p) * sgn;
    double tol = -1e-14 * diagLen * diagLen;
    return o1 >= tol && o2 >= tol && o3 >= tol;
}
static bool originalBoundaryEdge(const vector<int>& ring, int a, int b) {
    int k = (int)ring.size();
    for (int i=0; i<k; ++i) {
        int x=ring[i], y=ring[(i+1)%k];
        if ((x==a && y==b) || (x==b && y==a)) return true;
    }
    return false;
}

static bool basicTriangleOK(const Face& tf, const Vec3& avgN, double cosMin, double& area, double& loss) {
    if (tf.a==tf.b || tf.a==tf.c || tf.b==tf.c) return false;
    Vec3 cr = faceCross(tf);
    area = norm3(cr);
    if (!(area > epsArea)) return false;
    double nd = dot3(cr, avgN) / area;
    if (nd < 0) nd = -nd;
    if (nd < cosMin) return false;
    loss = 1.0 - min(1.0, nd);
    if (duplicateFaceOutside(tf.a, tf.b, tf.c)) return false;
    return true;
}

static bool verifyFill(const vector<int>& ring, const vector<Face>& rawTris, const Vec3& avgN,
                       double cosMin, double oldArea, double& areaRel, double& worstN) {
    if (rawTris.empty()) return false;
    unordered_map<uint64_t,int> newEC;
    newEC.reserve(rawTris.size()*3+4);
    double newArea = 0.0;
    worstN = 0.0;
    for (Face tf: rawTris) {
        Vec3 cr = faceCross(tf);
        double ar = norm3(cr);
        if (!(ar > epsArea)) return false;
        if (dot3(cr, avgN) < 0) swap(tf.b, tf.c);
        cr = faceCross(tf); ar = norm3(cr);
        double nd = dot3(cr, avgN) / ar;
        if (nd < cosMin) return false;
        if (duplicateFaceOutside(tf.a, tf.b, tf.c)) return false;
        newArea += ar;
        worstN = max(worstN, 1.0 - nd);
        int vv[3]={tf.a,tf.b,tf.c};
        for (int i=0;i<3;++i) newEC[edgeKey(vv[i],vv[(i+1)%3])]++;
    }
    for (auto &p: newEC) {
        int a = (int)(p.first >> 32);
        int b = (int)(p.first & 0xffffffffu);
        bool bd = originalBoundaryEdge(ring, a, b);
        int outside = edgeCountOutside(a,b);
        if (bd) {
            if (p.second != 1 || outside != 1) return false;
        } else {
            if (p.second != 2 || outside != 0) return false;
        }
    }
    areaRel = fabs(newArea - oldArea) / (oldArea + 1e-30);
    return areaRel <= 0.75;
}

static bool triangulateEar(const vector<int>& ringInput, const Vec3& avgN, double cosMin,
                           vector<Face>& tris, double& score) {
    int K = (int)ringInput.size();
    if (K < 3) return false;
    vector<int> ring = ringInput;
    Basis2 basis = makeBasis(avgN);
    vector<pair<double,double>> pt(K);
    for (int i=0;i<K;++i) pt[i] = projectBasis(P[ring[i]], basis);
    auto posInRing = [&](int v)->int {
        for (int i=0;i<K;++i) if (ring[i] == v) return i;
        return -1;
    };
    auto P2 = [&](int v)->pair<double,double> { return pt[posInRing(v)]; };
    double polyA = 0.0;
    for (int i=0;i<K;++i) polyA += orient2({0,0}, pt[i], pt[(i+1)%K]);
    if (fabs(polyA) <= 1e-18*diagLen*diagLen) return false;
    double sgn = (polyA > 0) ? 1.0 : -1.0;
    tris.clear();
    vector<int> cur = ring;
    score = 0.0;
    int guard = 0;
    while ((int)cur.size() > 3 && guard++ < K*K) {
        int m = (int)cur.size();
        int bestPos = -1;
        double bestCost = 1e100;
        Face bestTri{0,0,0};
        for (int i=0;i<m;++i) {
            int a = cur[(i+m-1)%m], b = cur[i], c = cur[(i+1)%m];
            auto pa = P2(a), pb = P2(b), pc = P2(c);
            double o = orient2(pa, pb, pc) * sgn;
            if (o <= 1e-18*diagLen*diagLen) continue;
            bool contains = false;
            for (int j=0;j<m;++j) {
                int q = cur[j];
                if (q==a || q==b || q==c) continue;
                if (insideTri2(P2(q), pa, pb, pc, sgn)) { contains = true; break; }
            }
            if (contains) continue;
            if (!originalBoundaryEdge(ringInput, a, c) && edgeCountOutside(a,c) != 0) continue;
            Face tf{a,b,c};
            double ar, loss;
            if (!basicTriangleOK(tf, avgN, cosMin, ar, loss)) continue;
            double diag = norm3(P[a]-P[c]) / (diagLen + 1e-30);
            double local = diag*diag + 25.0*loss + 0.00001*ar/(diagLen*diagLen+1e-30);
            if (local < bestCost) { bestCost = local; bestPos = i; bestTri = tf; }
        }
        if (bestPos < 0) return false;
        tris.push_back(bestTri);
        score += bestCost;
        cur.erase(cur.begin()+bestPos);
    }
    if (cur.size() != 3) return false;
    Face tf{cur[0],cur[1],cur[2]};
    double ar, loss;
    if (!basicTriangleOK(tf, avgN, cosMin, ar, loss)) return false;
    tris.push_back(tf);
    score += 25.0*loss;
    return true;
}

static bool triangulateBest(const vector<int>& ring, const Vec3& avgN, double cosMin,
                            vector<Face>& bestTris, double& bestScore) {
    bestScore = 1e100;
    bestTris.clear();
    for (int pass=0; pass<2; ++pass) {
        vector<int> r = ring;
        if (pass) reverse(r.begin(), r.end());
        vector<Face> tris;
        double sc = 0.0;
        if (!triangulateEar(r, avgN, cosMin, tris, sc)) continue;
        if (sc < bestScore) { bestScore = sc; bestTris = tris; }
    }
    return !bestTris.empty();
}

static bool baryDepthInTri(const Vec3& q, const Face& t, int axis, double& depthOut) {
    auto p = screen2(q, axis);
    auto a = screen2(P[t.a], axis);
    auto b = screen2(P[t.b], axis);
    auto c = screen2(P[t.c], axis);
    double den = orient2(a,b,c);
    if (fabs(den) <= 1e-18*diagLen*diagLen) return false;
    double w0 = orient2(b,c,p) / den;
    double w1 = orient2(c,a,p) / den;
    double w2 = 1.0 - w0 - w1;
    double tol = -1e-9;
    if (w0 < tol || w1 < tol || w2 < tol) return false;
    depthOut = w0*coord(P[t.a], axis) + w1*coord(P[t.b], axis) + w2*coord(P[t.c], axis);
    return true;
}

static bool visualDepthGuard(int v, const vector<Face>& tris, double depthLimit, double& worstDepth) {
    worstDepth = 0.0;
    const Vec3 q = P[v];
    for (int axis=0; axis<3; ++axis) {
        double best = 1e100;
        for (const Face& t: tris) {
            double d;
            if (baryDepthInTri(q, t, axis, d)) {
                best = min(best, fabs(d - coord(q, axis)));
            }
        }
        if (best < 1e90) {
            worstDepth = max(worstDepth, best);
            if (best > depthLimit) return false;
        }
    }
    return true;
}

static CandidateEval evaluateVertex(int v, double cosMin, int maxValence, double depthLimit, double areaLimit) {
    CandidateEval ev;
    ev.v = v;
    if (!aliveV[v]) return ev;
    vector<int> inc = incidentFaces(v);
    int k = (int)inc.size();
    if (k < 3 || k > maxValence) return ev;
    markFaces(inc);

    vector<int> ring;
    if (!buildLinkCycle(v, inc, ring)) return ev;
    if ((int)ring.size() != k) return ev;

    double nearestGeom = 1e100;
    for (int u: ring) {
        if (!aliveV[u]) return ev;
        nearestGeom = min(nearestGeom, norm3(P[v]-P[u]));
    }
    vector<pair<int,int>> redist;
    redist.reserve(ownedOriginals[v].size());
    double bestRedist = 0.0;
    for (int orig: ownedOriginals[v]) {
        int bestU = -1;
        double bestD = 1e100;
        for (int u: ring) {
            double d = norm3(P[orig] - P[u]);
            if (d < bestD) { bestD = d; bestU = u; }
        }
        if (bestU < 0 || bestD > radiusLimit) return ev;
        bestRedist = max(bestRedist, bestD);
        redist.push_back({orig, bestU});
    }

    Vec3 oldSum{0,0,0};
    double oldArea = 0.0;
    double patchNormalSpread = 0.0;
    for (int id: inc) {
        Vec3 cr = faceCross(F[id]);
        double ar = norm3(cr);
        if (!(ar > epsArea)) return ev;
        oldArea += ar;
        oldSum = oldSum + cr;
    }
    double sumLen = norm3(oldSum);
    if (!(sumLen > epsArea)) return ev;
    Vec3 avgN = oldSum / sumLen;
    for (int id: inc) {
        Vec3 cr = faceCross(F[id]);
        double ar = norm3(cr);
        double nd = dot3(cr, avgN) / ar;
        if (nd < cosMin) return ev;
        patchNormalSpread = max(patchNormalSpread, 1.0 - nd);
    }

    vector<Face> tris;
    double triangScore;
    if (!triangulateBest(ring, avgN, cosMin, tris, triangScore)) return ev;

    for (Face& t: tris) {
        if (dot3(faceCross(t), avgN) < 0) swap(t.b, t.c);
    }

    double areaRel = 0.0, worstN = 0.0;
    if (!verifyFill(ring, tris, avgN, cosMin, oldArea, areaRel, worstN)) return ev;
    if (areaRel > areaLimit) return ev;

    double worstDepth = 0.0;
    if (!visualDepthGuard(v, tris, depthLimit, worstDepth)) return ev;

    double sharpPenalty = 0.0;
    vector<int> local = inc;
    for (int u: ring) {
        int cnt = 0;
        Vec3 ns[2];
        for (int id: local) if (hasVertex(F[id],u)) {
            Vec3 cr = faceCross(F[id]); double ar = norm3(cr);
            if (ar > epsArea && cnt < 2) ns[cnt++] = cr / ar;
        }
        if (cnt == 2) sharpPenalty = max(sharpPenalty, 1.0 - dot3(ns[0], ns[1]));
    }

    ev.ok = true;
    ev.maxRedistribute = bestRedist;
    ev.redistribute = redist;
    ev.ring = ring;
    ev.add = tris;
    double density = nearestGeom / (radiusLimit + 1e-30);
    double radUse = bestRedist / (radiusLimit + 1e-30);
    ev.cost = 0.20*k + 0.45*density + 0.25*radUse + 24.0*worstN + 0.65*areaRel
            + 8.0*(worstDepth/(depthLimit + 1e-30)) + 5.0*patchNormalSpread + 1.5*sharpPenalty
            + 0.02*triangScore;
    return ev;
}

struct QueueItem {
    double cost;
    int v;
    unsigned stamp;
    bool operator<(const QueueItem& o) const { return cost > o.cost; }
};
static priority_queue<QueueItem> pq;

static void pushVertex(int v, double cosMin, int maxValence, double depthLimit, double areaLimit) {
    if (!aliveV[v]) return;
    CandidateEval ev = evaluateVertex(v, cosMin, maxValence, depthLimit, areaLimit);
    if (ev.ok) pq.push({ev.cost, v, vstamp[v]});
}

static void applyDelete(const CandidateEval& ev) {
    int v = ev.v;
    vector<int> inc = incidentFaces(v);
    for (int id: inc) if (aliveF[id]) {
        aliveF[id] = 0;
        --liveF;
    }
    aliveV[v] = 0;
    --liveV;
    for (auto &mv: ev.redistribute) {
        int orig = mv.first, to = mv.second;
        if (to >= 0 && to < N && aliveV[to]) ownedOriginals[to].push_back(orig);
    }
    vector<int>().swap(ownedOriginals[v]);

    for (const Face& t: ev.add) {
        int id = (int)F.size();
        F.push_back(t);
        aliveF.push_back(1);
        faceMark.push_back(0);
        ++liveF;
        VF[t.a].push_back(id);
        VF[t.b].push_back(id);
        VF[t.c].push_back(id);
    }
    VF[v].clear();
    ++vstamp[v];
    for (int u: ev.ring) ++vstamp[u];
}

static void runPass(double cosMin, int maxValence, double depthFrac, double areaLimit, int targetV, int opLimit) {
    double depthLimit = depthFrac * diagLen;
    while (!pq.empty()) pq.pop();
    for (int v=0; v<N; ++v) if (aliveV[v]) pushVertex(v, cosMin, maxValence, depthLimit, areaLimit);
    int ops = 0;
    while (!pq.empty() && liveV > targetV && ops < opLimit) {
        QueueItem q = pq.top(); pq.pop();
        if (!aliveV[q.v] || vstamp[q.v] != q.stamp) continue;
        CandidateEval ev = evaluateVertex(q.v, cosMin, maxValence, depthLimit, areaLimit);
        if (!ev.ok) continue;
        if (ev.cost > q.cost * 1.50 + 1e-12) {
            pushVertex(q.v, cosMin, maxValence, depthLimit, areaLimit);
            continue;
        }
        applyDelete(ev);
        ++ops;
        for (int u: ev.ring) {
            pushVertex(u, cosMin, maxValence, depthLimit, areaLimit);
            vector<int> in = incidentFaces(u);
            for (int id: in) {
                int a[3] = {F[id].a, F[id].b, F[id].c};
                for (int j=0; j<3; ++j) pushVertex(a[j], cosMin, maxValence, depthLimit, areaLimit);
            }
        }
        if ((ops & 511) == 0) rebuildVF();
    }
    rebuildVF();
}

static bool closedManifoldCheck() {
    unordered_map<uint64_t,int> ec;
    ec.reserve((size_t)liveF*3 + 16);
    vector<int> used(N,0);
    for (int id=0; id<(int)F.size(); ++id) if (aliveF[id]) {
        Face f = F[id];
        if (f.a<0||f.a>=N||f.b<0||f.b>=N||f.c<0||f.c>=N) return false;
        if (!aliveV[f.a] || !aliveV[f.b] || !aliveV[f.c]) return false;
        if (f.a==f.b || f.a==f.c || f.b==f.c) return false;
        if (faceArea2(f) <= epsArea) return false;
        int a[3]={f.a,f.b,f.c};
        for (int i=0;i<3;++i) ec[edgeKey(a[i],a[(i+1)%3])]++;
        used[f.a]=used[f.b]=used[f.c]=1;
    }
    for (auto &p: ec) if (p.second != 2) return false;
    for (int v=0; v<N; ++v) if (aliveV[v]) {
        if (!used[v]) return false;
        vector<int> in = incidentFaces(v);
        if (in.empty()) return false;
        vector<int> verts;
        vector<pair<int,int>> ed;
        for (int id: in) {
            int other[2], k=0;
            int a[3]={F[id].a,F[id].b,F[id].c};
            for (int i=0;i<3;++i) if (a[i] != v) {
                if (k>=2) return false;
                other[k++] = a[i];
            }
            if (k != 2 || other[0] == other[1]) return false;
            verts.push_back(other[0]); verts.push_back(other[1]);
            ed.push_back({other[0],other[1]});
        }
        sort(verts.begin(), verts.end());
        verts.erase(unique(verts.begin(), verts.end()), verts.end());
        vector<int> deg(verts.size(),0);
        for (auto &e: ed) {
            int ia = lower_bound(verts.begin(), verts.end(), e.first) - verts.begin();
            int ib = lower_bound(verts.begin(), verts.end(), e.second) - verts.begin();
            if (ia<0||ib<0||ia==(int)verts.size()||ib==(int)verts.size()) return false;
            deg[ia]++; deg[ib]++;
        }
        for (int d: deg) if (d != 2) return false;
    }
    return true;
}

static void resetOriginal(const vector<Face>& origF) {
    F = origF;
    aliveV.assign(N,1);
    aliveF.assign(F.size(),1);
    VF.assign(N,{});
    vstamp.assign(N,1);
    ownedOriginals.assign(N,{});
    for (int i=0;i<N;++i) ownedOriginals[i].push_back(i);
    faceMark.assign(F.size(),0);
    markToken = 1;
    liveV = N;
    liveF = (int)F.size();
    rebuildVF();
}

static void emitMesh() {
    vector<int> id(N,-1);
    vector<Vec3> outV;
    outV.reserve(liveV);
    for (int i=0;i<N;++i) if (aliveV[i]) {
        id[i] = (int)outV.size()+1;
        outV.push_back(P[i]);
    }
    vector<Face> outF;
    outF.reserve(liveF);
    for (int i=0;i<(int)F.size();++i) if (aliveF[i]) {
        Face f=F[i];
        if (id[f.a]>0 && id[f.b]>0 && id[f.c]>0 && id[f.a]!=id[f.b] && id[f.a]!=id[f.c] && id[f.b]!=id[f.c])
            outF.push_back({id[f.a], id[f.b], id[f.c]});
    }
    cout.setf(ios::fixed);
    cout.precision(10);
    cout << outV.size() << ' ' << outF.size() << '\n';
    for (const Vec3& p: outV) cout << "v " << p.x << ' ' << p.y << ' ' << p.z << '\n';
    for (const Face& f: outF) cout << "f " << f.a << ' ' << f.b << ' ' << f.c << '\n';
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    if (!(cin >> N >> M0)) return 0;
    P.resize(N);
    Vec3 mn{1e100,1e100,1e100}, mx{-1e100,-1e100,-1e100};
    string tag;
    for (int i=0;i<N;++i) {
        cin >> tag >> P[i].x >> P[i].y >> P[i].z;
        mn.x = min(mn.x, P[i].x); mn.y = min(mn.y, P[i].y); mn.z = min(mn.z, P[i].z);
        mx.x = max(mx.x, P[i].x); mx.y = max(mx.y, P[i].y); mx.z = max(mx.z, P[i].z);
    }
    vector<Face> origF(M0);
    for (int i=0;i<M0;++i) {
        cin >> tag >> origF[i].a >> origF[i].b >> origF[i].c;
        --origF[i].a; --origF[i].b; --origF[i].c;
    }

    diagLen = norm3(mx-mn);
    if (!(diagLen > 0)) diagLen = 1.0;
    epsArea = 1e-20 * diagLen * diagLen;
    radiusLimit = 0.0490 * diagLen;

    resetOriginal(origF);

    if (N <= 8) {
        emitMesh();
        return 0;
    }

    int t1 = max(8, (int)ceil(N * 0.55));
    int t2 = max(8, (int)ceil(N * 0.34));
    int t3 = max(8, (int)ceil(N * 0.21));
    int t4 = max(8, (int)ceil(N * 0.13));
    int t5 = max(8, (int)ceil(N * 0.085));

    runPass(0.995, 6, 0.0045, 0.18, t1, N);
    runPass(0.985, 7, 0.0080, 0.25, t2, N);
    runPass(0.970, 8, 0.0140, 0.35, t3, 2*N);
    runPass(0.945,10, 0.0220, 0.50, t4, 3*N);
    runPass(0.920,12, 0.0300, 0.65, t5, 4*N);

    if (!closedManifoldCheck()) {
        resetOriginal(origF);
    }

    emitMesh();
    return 0;
}