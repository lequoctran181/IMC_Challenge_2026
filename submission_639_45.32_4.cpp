#include <bits/stdc++.h>
using namespace std;

// IMC simplifygeometry - topology-preserving QEM/normal-aware edge collapse.
// C++17, no external dependencies.

struct Vec3 {
    double x, y, z;
    Vec3() : x(0), y(0), z(0) {}
    Vec3(double X, double Y, double Z) : x(X), y(Y), z(Z) {}
    inline Vec3 operator+(const Vec3& o) const { return Vec3(x+o.x,y+o.y,z+o.z); }
    inline Vec3 operator-(const Vec3& o) const { return Vec3(x-o.x,y-o.y,z-o.z); }
    inline Vec3 operator*(double s) const { return Vec3(x*s,y*s,z*s); }
    inline Vec3 operator/(double s) const { return Vec3(x/s,y/s,z/s); }
    inline Vec3& operator+=(const Vec3& o){ x+=o.x; y+=o.y; z+=o.z; return *this; }
};
static inline double dotv(const Vec3& a, const Vec3& b){ return a.x*b.x + a.y*b.y + a.z*b.z; }
static inline Vec3 crossv(const Vec3& a, const Vec3& b){ return Vec3(a.y*b.z-a.z*b.y, a.z*b.x-a.x*b.z, a.x*b.y-a.y*b.x); }
static inline double norm2(const Vec3& a){ return dotv(a,a); }
static inline double normv(const Vec3& a){ return sqrt(norm2(a)); }

struct Quadric {
    // symmetric 4x4 matrix: [x y z 1] Q [x y z 1]^T
    double a00=0,a01=0,a02=0,a03=0,a11=0,a12=0,a13=0,a22=0,a23=0,a33=0;
    inline void add(const Quadric& q){
        a00+=q.a00; a01+=q.a01; a02+=q.a02; a03+=q.a03;
        a11+=q.a11; a12+=q.a12; a13+=q.a13; a22+=q.a22; a23+=q.a23; a33+=q.a33;
    }
    inline void addPlane(const Vec3& n, double d, double w){
        double nx=n.x, ny=n.y, nz=n.z;
        a00 += w*nx*nx; a01 += w*nx*ny; a02 += w*nx*nz; a03 += w*nx*d;
        a11 += w*ny*ny; a12 += w*ny*nz; a13 += w*ny*d;
        a22 += w*nz*nz; a23 += w*nz*d; a33 += w*d*d;
    }
    inline double eval(const Vec3& p) const {
        double x=p.x,y=p.y,z=p.z;
        return a00*x*x + 2*a01*x*y + 2*a02*x*z + 2*a03*x
             + a11*y*y + 2*a12*y*z + 2*a13*y
             + a22*z*z + 2*a23*z + a33;
    }
    bool solve(Vec3& out) const {
        // Solve A p = -b for the upper-left 3x3 block.
        double A00=a00, A01=a01, A02=a02;
        double A10=a01, A11=a11, A12=a12;
        double A20=a02, A21=a12, A22=a22;
        double b0=-a03, b1=-a13, b2=-a23;
        double det = A00*(A11*A22 - A12*A21) - A01*(A10*A22 - A12*A20) + A02*(A10*A21 - A11*A20);
        if (fabs(det) < 1e-15) return false;
        double id = 1.0/det;
        double dx = b0*(A11*A22 - A12*A21) - A01*(b1*A22 - A12*b2) + A02*(b1*A21 - A11*b2);
        double dy = A00*(b1*A22 - A12*b2) - b0*(A10*A22 - A12*A20) + A02*(A10*b2 - b1*A20);
        double dz = A00*(A11*b2 - b1*A21) - A01*(A10*b2 - b1*A20) + b0*(A10*A21 - A11*A20);
        out = Vec3(dx*id, dy*id, dz*id);
        return isfinite(out.x) && isfinite(out.y) && isfinite(out.z);
    }
};

struct Vertex {
    Vec3 p;
    Vec3 mn, mx;       // conservative bbox of original vertices absorbed here
    Quadric q;
    int cnt = 1;
    int ver = 1;
    bool alive = true;
};

struct Face {
    int v[3];
    Vec3 n;
    double area2 = 0.0;
    bool alive = true;
};

struct EdgeItem {
    double key;
    int u, v;
    int vu, vv;
    bool operator<(const EdgeItem& other) const { return key > other.key; } // min-heap
};

static int N0, F0;
static vector<Vertex> V;
static vector<Face> F;
static vector<vector<int>> Inc;
static priority_queue<EdgeItem> PQ;
static int aliveVertices = 0, aliveFaces = 0;
static double bboxDiag = 0.0, hausTol = 0.0, hausTolSafe = 0.0;
static double g_minDot = 0.70;
static chrono::steady_clock::time_point g_start;
static double g_timeLimit = 19.2;
static bool g_fastMode = false;
static bool g_queueActive = true;
static double g_detailP90 = 0.0;

static inline bool timeExceeded(){
    double t = chrono::duration<double>(chrono::steady_clock::now() - g_start).count();
    return t > g_timeLimit;
}

static vector<char> slurp_stdin() {
    vector<char> buf;
    buf.reserve(1 << 26);
    char chunk[1 << 16];
    size_t n;
    while ((n = fread(chunk, 1, sizeof(chunk), stdin)) > 0) buf.insert(buf.end(), chunk, chunk+n);
    buf.push_back('\0');
    return buf;
}

static inline void skip_ws(char*& p){ while (*p==' ' || *p=='\n' || *p=='\r' || *p=='\t') ++p; }

static void recomputeFace(int id){
    Face &f = F[id];
    Vec3 a = V[f.v[0]].p, b = V[f.v[1]].p, c = V[f.v[2]].p;
    Vec3 cr = crossv(b-a, c-a);
    double len = normv(cr);
    f.area2 = len;
    if (len > 0) f.n = cr / len;
    else f.n = Vec3(0,0,0);
}

static void loadMesh(){
    vector<char> buf = slurp_stdin();
    char* p = buf.data();
    N0 = (int)strtol(p, &p, 10);
    F0 = (int)strtol(p, &p, 10);
    V.assign(N0, Vertex());
    F.assign(F0, Face());
    Vec3 mn(1e100,1e100,1e100), mx(-1e100,-1e100,-1e100);
    for (int i=0;i<N0;i++){
        skip_ws(p); if (*p=='v') ++p;
        double x = strtod(p,&p), y=strtod(p,&p), z=strtod(p,&p);
        V[i].p = Vec3(x,y,z);
        V[i].mn = V[i].mx = V[i].p;
        mn.x=min(mn.x,x); mn.y=min(mn.y,y); mn.z=min(mn.z,z);
        mx.x=max(mx.x,x); mx.y=max(mx.y,y); mx.z=max(mx.z,z);
    }
    vector<int> deg(N0,0);
    for (int i=0;i<F0;i++){
        skip_ws(p); if (*p=='f') ++p;
        int a=(int)strtol(p,&p,10)-1, b=(int)strtol(p,&p,10)-1, c=(int)strtol(p,&p,10)-1;
        F[i].v[0]=a; F[i].v[1]=b; F[i].v[2]=c;
        deg[a]++; deg[b]++; deg[c]++;
    }
    bboxDiag = sqrt((mx.x-mn.x)*(mx.x-mn.x)+(mx.y-mn.y)*(mx.y-mn.y)+(mx.z-mn.z)*(mx.z-mn.z));
    hausTol = 0.05 * bboxDiag;
    hausTolSafe = hausTol * 0.985;
    Inc.resize(N0);
    for (int i=0;i<N0;i++) Inc[i].reserve((size_t)deg[i] + 8);
    for (int i=0;i<F0;i++){
        recomputeFace(i);
        for (int k=0;k<3;k++) Inc[F[i].v[k]].push_back(i);
        double area = 0.5 * F[i].area2;
        if (area <= 0) continue;
        Vec3 n = F[i].n;
        double d = -dotv(n, V[F[i].v[0]].p);
        // Area weighting is stable, but put a small floor so tiny feature triangles still count.
        double w = sqrt(area) + 1e-8;
        Quadric q; q.addPlane(n, d, w);
        V[F[i].v[0]].q.add(q); V[F[i].v[1]].q.add(q); V[F[i].v[2]].q.add(q);
    }
    aliveVertices = N0; aliveFaces = F0;
}

static inline uint64_t packEdge(int a, int b){
    if (a>b) swap(a,b);
    return (uint64_t)(uint32_t)a << 32 | (uint32_t)b;
}
static inline int eA(uint64_t e){ return (int)(e >> 32); }
static inline int eB(uint64_t e){ return (int)(e & 0xffffffffu); }

static inline bool faceHasVertex(int fid, int x){
    const Face& f = F[fid];
    return f.alive && (f.v[0]==x || f.v[1]==x || f.v[2]==x);
}

static void filteredInc(int v, vector<int>& out){
    out.clear();
    if (!V[v].alive) return;
    vector<int>& src = Inc[v];
    for (int fid : src) if ((unsigned)fid < (unsigned)F.size() && faceHasVertex(fid, v)) out.push_back(fid);
    if (out.size() > 1) {
        sort(out.begin(), out.end());
        out.erase(unique(out.begin(), out.end()), out.end());
    }
}

static inline double farthestBBoxDist2(const Vec3& p, const Vec3& mn, const Vec3& mx){
    double dx = max(fabs(p.x-mn.x), fabs(p.x-mx.x));
    double dy = max(fabs(p.y-mn.y), fabs(p.y-mx.y));
    double dz = max(fabs(p.z-mn.z), fabs(p.z-mx.z));
    return dx*dx + dy*dy + dz*dz;
}

struct EvalRes { bool ok=false; double key=0; Vec3 p; };

static bool collectTopology(int u, int v, vector<int>& incU, vector<int>& incV, vector<int>& common, int opp[2]){
    filteredInc(u, incU); filteredInc(v, incV);
    if (incU.empty() || incV.empty()) return false;
    common.clear();
    int oi=0;
    // degrees are usually tiny; two-pointer after sorting.
    size_t i=0,j=0;
    while (i<incU.size() && j<incV.size()){
        if (incU[i]==incV[j]){
            int fid = incU[i];
            const Face& f=F[fid];
            if (faceHasVertex(fid,u) && faceHasVertex(fid,v)){
                common.push_back(fid);
                int o=-1;
                for (int k=0;k<3;k++) if (f.v[k]!=u && f.v[k]!=v) o=f.v[k];
                if (oi<2) opp[oi]=o;
                oi++;
            }
            ++i; ++j;
        } else if (incU[i]<incV[j]) ++i; else ++j;
    }
    if ((int)common.size()!=2 || oi!=2 || opp[0]<0 || opp[1]<0 || opp[0]==opp[1]) return false;

    vector<int> nu, nv;
    nu.reserve(incU.size()*2); nv.reserve(incV.size()*2);
    for (int fid: incU){
        const Face& f=F[fid];
        for (int k=0;k<3;k++) if (f.v[k]!=u && f.v[k]!=v) nu.push_back(f.v[k]);
    }
    for (int fid: incV){
        const Face& f=F[fid];
        for (int k=0;k<3;k++) if (f.v[k]!=u && f.v[k]!=v) nv.push_back(f.v[k]);
    }
    sort(nu.begin(), nu.end()); nu.erase(unique(nu.begin(), nu.end()), nu.end());
    sort(nv.begin(), nv.end()); nv.erase(unique(nv.begin(), nv.end()), nv.end());
    int inter=0; bool has0=false, has1=false;
    i=0; j=0;
    while (i<nu.size() && j<nv.size()){
        if (nu[i]==nv[j]){
            int x=nu[i]; inter++;
            if (x==opp[0]) has0=true;
            if (x==opp[1]) has1=true;
            ++i; ++j;
        } else if (nu[i]<nv[j]) ++i; else ++j;
    }
    return inter==2 && has0 && has1;
}

static inline bool sameVertexAfterReplace(const Face& f, int u, int v){
    int a=f.v[0], b=f.v[1], c=f.v[2];
    if (a==v) a=u; if (b==v) b=u; if (c==v) c=u;
    return a==b || b==c || a==c;
}

static EvalRes evaluateEdge(int u, int v){
    EvalRes best;
    if (u==v || u<0 || v<0 || u>=N0 || v>=N0 || !V[u].alive || !V[v].alive) return best;
    static thread_local vector<int> incU, incV, common;
    int opp[2] = {-1,-1};
    if (!collectTopology(u, v, incU, incV, common, opp)) return best;

    Vec3 mn(min(V[u].mn.x,V[v].mn.x), min(V[u].mn.y,V[v].mn.y), min(V[u].mn.z,V[v].mn.z));
    Vec3 mx(max(V[u].mx.x,V[v].mx.x), max(V[u].mx.y,V[v].mx.y), max(V[u].mx.z,V[v].mx.z));
    Quadric q = V[u].q; q.add(V[v].q);
    vector<Vec3> cand;
    cand.reserve(6);
    cand.push_back((V[u].p * (double)V[u].cnt + V[v].p * (double)V[v].cnt) / (double)(V[u].cnt + V[v].cnt));
    cand.push_back((V[u].p + V[v].p) * 0.5);
    cand.push_back(V[u].p);
    cand.push_back(V[v].p);
    Vec3 opt;
    if (q.solve(opt)) cand.push_back(opt);
    // Clamp an over-eager QEM optimum back to the merged cluster bbox enlarged by tolerance.
    if (cand.size()==5){
        Vec3 c = cand.back();
        double pad = hausTolSafe * 0.25;
        c.x = min(max(c.x, mn.x-pad), mx.x+pad);
        c.y = min(max(c.y, mn.y-pad), mx.y+pad);
        c.z = min(max(c.z, mn.z-pad), mx.z+pad);
        cand.back() = c;
    }
    
    // common face ids for O(1) skip
    int c0=common[0], c1=common[1];
    double tol2 = hausTolSafe * hausTolSafe;
    double edgeLen = normv(V[u].p - V[v].p);
    double localScale = max(1e-12, edgeLen);

    for (const Vec3& p : cand){
        if (!isfinite(p.x) || !isfinite(p.y) || !isfinite(p.z)) continue;
        if (farthestBBoxDist2(p, mn, mx) > tol2) continue;
        double qerr = max(0.0, q.eval(p));
        double sumArea = 1e-18, normalCost = 0.0, maxBend = 0.0, skinnyPenalty = 0.0;
        bool ok = true;
        auto scanFace = [&](int fid){
            if (!ok || fid==c0 || fid==c1) return;
            const Face& f = F[fid];
            if (!f.alive) return;
            if (!(f.v[0]==u || f.v[1]==u || f.v[2]==u || f.v[0]==v || f.v[1]==v || f.v[2]==v)) return;
            if (sameVertexAfterReplace(f,u,v)) { ok=false; return; }
            Vec3 pts[3];
            for (int k=0;k<3;k++){
                int id=f.v[k];
                pts[k] = (id==u || id==v) ? p : V[id].p;
            }
            Vec3 cr = crossv(pts[1]-pts[0], pts[2]-pts[0]);
            double len = normv(cr);
            if (!(len > 1e-18)) { ok=false; return; }
            Vec3 nn = cr / len;
            double d = dotv(nn, f.n);
            if (d < g_minDot) { ok=false; return; }
            double bend = max(0.0, 1.0-d);
            double a = max(1e-18, f.area2);
            sumArea += a;
            normalCost += a * bend * bend;
            maxBend = max(maxBend, bend);
            // Discourage extremely skinny triangles; they are valid but usually bad for raster stability.
            double l0=norm2(pts[1]-pts[0]), l1=norm2(pts[2]-pts[1]), l2=norm2(pts[0]-pts[2]);
            double lm=max(l0,max(l1,l2));
            if (lm > 0) {
                double quality = (len*len) / (lm*lm + 1e-30);
                if (quality < 1e-8) skinnyPenalty += (1e-8-quality)*1e5;
            }
        };
        for (int fid: incU) scanFace(fid);
        for (int fid: incV) scanFace(fid);
        if (!ok) continue;
        double radius = sqrt(max(0.0, farthestBBoxDist2(p,mn,mx))) / max(hausTolSafe,1e-12);
        double avgNorm = sqrt(normalCost / sumArea);
        // Key combines quadric distance, normal-map perturbation, cluster radius and local compactness.
        // It is deliberately scale-normalized because all inputs live in a unit sphere.
        double key = sqrt(qerr / (sumArea + 1e-18)) * 16.0
                   + avgNorm * 140.0
                   + maxBend * 35.0
                   + radius * 0.20
                   + edgeLen * 0.02
                   + skinnyPenalty;
        // Endpoint collapses are useful on creases, but prefer the endpoint that moves less visually.
        if (&p == &cand[2] || &p == &cand[3]) key += 1e-7;
        if (!best.ok || key < best.key){ best.ok=true; best.key=key; best.p=p; }
    }
    return best;
}

static void pushEdge(int a, int b){
    if (a==b || a<0 || b<0 || a>=N0 || b>=N0) return;
    if (!V[a].alive || !V[b].alive) return;
    EvalRes e = evaluateEdge(a,b);
    if (!e.ok) return;
    PQ.push(EdgeItem{e.key, a, b, V[a].ver, V[b].ver});
}

static void rebuildQueue(bool compactInc){
    priority_queue<EdgeItem> empty; PQ.swap(empty);
    if (compactInc){
        for (int i=0;i<N0;i++) if (V[i].alive) Inc[i].clear();
        for (int i=0;i<F0;i++) if (F[i].alive){
            Inc[F[i].v[0]].push_back(i); Inc[F[i].v[1]].push_back(i); Inc[F[i].v[2]].push_back(i);
        }
    }
    vector<uint64_t> edges;
    edges.reserve((size_t)aliveFaces * 3);
    for (int i=0;i<F0;i++) if (F[i].alive){
        int a=F[i].v[0], b=F[i].v[1], c=F[i].v[2];
        if (V[a].alive && V[b].alive && V[c].alive && a!=b && b!=c && a!=c){
            edges.push_back(packEdge(a,b)); edges.push_back(packEdge(b,c)); edges.push_back(packEdge(c,a));
        }
    }
    sort(edges.begin(), edges.end());
    edges.erase(unique(edges.begin(), edges.end()), edges.end());
    for (uint64_t e: edges){
        if (timeExceeded()) break;
        pushEdge(eA(e), eB(e));
    }
}

static void collapseEdge(int u, int v, const Vec3& p){
    if (!V[u].alive || !V[v].alive) return;
    static thread_local vector<int> incU, incV, common;
    int opp[2];
    if (!collectTopology(u, v, incU, incV, common, opp)) return;
    int c0=common[0], c1=common[1];

    vector<int> affectedFaces;
    affectedFaces.reserve(incU.size()+incV.size());
    for (int fid: incU) if (fid!=c0 && fid!=c1) affectedFaces.push_back(fid);
    for (int fid: incV) if (fid!=c0 && fid!=c1) affectedFaces.push_back(fid);
    sort(affectedFaces.begin(), affectedFaces.end());
    affectedFaces.erase(unique(affectedFaces.begin(), affectedFaces.end()), affectedFaces.end());

    vector<int> affectedVerts;
    affectedVerts.reserve(affectedFaces.size()*3 + 8);
    affectedVerts.push_back(u); affectedVerts.push_back(v); affectedVerts.push_back(opp[0]); affectedVerts.push_back(opp[1]);
    for (int fid: affectedFaces){
        for (int k=0;k<3;k++) affectedVerts.push_back(F[fid].v[k]);
    }

    // remove the two triangles adjacent to collapsed edge
    if (F[c0].alive){ F[c0].alive=false; aliveFaces--; }
    if (F[c1].alive){ F[c1].alive=false; aliveFaces--; }

    // merge vertex v into u
    V[u].p = p;
    V[u].mn = Vec3(min(V[u].mn.x,V[v].mn.x), min(V[u].mn.y,V[v].mn.y), min(V[u].mn.z,V[v].mn.z));
    V[u].mx = Vec3(max(V[u].mx.x,V[v].mx.x), max(V[u].mx.y,V[v].mx.y), max(V[u].mx.z,V[v].mx.z));
    V[u].cnt += V[v].cnt;
    V[u].q.add(V[v].q);
    V[v].alive = false;
    aliveVertices--;

    for (int fid: affectedFaces){
        Face& f=F[fid];
        if (!f.alive) continue;
        for (int k=0;k<3;k++) if (f.v[k]==v) f.v[k]=u;
        // Degenerate faces should be absent by evaluation/link condition, but guard anyway.
        if (f.v[0]==f.v[1] || f.v[1]==f.v[2] || f.v[0]==f.v[2]){
            f.alive=false; aliveFaces--;
        } else recomputeFace(fid);
    }

    Inc[u].insert(Inc[u].end(), Inc[v].begin(), Inc[v].end());
    Inc[v].clear();
    sort(affectedVerts.begin(), affectedVerts.end());
    affectedVerts.erase(unique(affectedVerts.begin(), affectedVerts.end()), affectedVerts.end());
    for (int x: affectedVerts) if ((unsigned)x < (unsigned)N0 && V[x].alive) V[x].ver++;

    // Push every edge of the one-ring affected by the new vertex position in priority-queue mode.
    if (g_queueActive) {
        static thread_local vector<int> newInc;
        filteredInc(u, newInc);
        for (int fid: newInc){
            if (!F[fid].alive) continue;
            int a=F[fid].v[0], b=F[fid].v[1], c=F[fid].v[2];
            pushEdge(a,b); pushEdge(b,c); pushEdge(c,a);
        }
    }
}


static double estimateDetailP90(){
    // Area-weighted one-ring normal variance.  It is near zero on smooth/planar oversampled
    // meshes, but rises on high-frequency relief where the normal map needs more vertices.
    vector<float> vals;
    vals.reserve(N0);
    for (int i=0;i<N0;i++){
        const vector<int>& fs = Inc[i];
        if (fs.empty()) { vals.push_back(0.0f); continue; }
        Vec3 s(0,0,0); double sw=0.0;
        for (int fid: fs){
            const Face& f = F[fid];
            if (!f.alive) continue;
            double w = max(1e-18, f.area2);
            s += f.n * w;
            sw += w;
        }
        double l = normv(s);
        if (l <= 1e-30 || sw <= 0) { vals.push_back(0.0f); continue; }
        Vec3 m = s / l;
        double var=0.0;
        for (int fid: fs){
            const Face& f = F[fid];
            if (!f.alive) continue;
            double w = max(1e-18, f.area2);
            double d = max(-1.0, min(1.0, dotv(m, f.n)));
            var += w * max(0.0, 1.0 - d);
        }
        vals.push_back((float)(var / sw));
    }
    if (vals.empty()) return 0.0;
    size_t k = (size_t)(0.90 * (double)(vals.size()-1));
    nth_element(vals.begin(), vals.begin()+k, vals.end());
    return vals[k];
}

static int chooseTarget(){
    if (N0 <= 9) return max(4, N0-1);
    double r;
    if (N0 <= 6000) r = 0.085;
    else if (N0 <= 30000) r = 0.082;
    else if (N0 <= 60000) r = 0.080;
    else if (N0 <= 450000) r = 0.080;
    else r = 0.080;
    // Extra budget for meshes whose normal map contains widespread high-frequency detail.
    // A dense smooth sphere/torus has p90 variance below about 0.001; bumpy relief is far higher.
    double extra = 0.0;
    if (g_detailP90 > 0.0012) extra = min(0.085, (g_detailP90 - 0.0012) * 8.0);
    r += extra;
    int t = (int)ceil(N0 * r);
    // Do not chase impossible tiny triangulations on high-genus inputs.
    t = max(t, 4);
    return t;
}


static void rebuildIncOnly(){
    for (int i=0;i<N0;i++) if (V[i].alive) Inc[i].clear();
    for (int i=0;i<F0;i++) if (F[i].alive){
        Inc[F[i].v[0]].push_back(i);
        Inc[F[i].v[1]].push_back(i);
        Inc[F[i].v[2]].push_back(i);
    }
}

static void greedyLargeSimplify(int finalTarget){
    g_queueActive = false;
    vector<int> touched(N0, 0);
    int sweepId = 1;
    struct Stage { double ratio; double minDot; };
    vector<Stage> stages;
    stages.push_back({max((double)finalTarget/N0, 0.300), 0.965});
    stages.push_back({max((double)finalTarget/N0, 0.200), 0.900});
    stages.push_back({max((double)finalTarget/N0, 0.135), 0.800});
    stages.push_back({max((double)finalTarget/N0, 0.100), 0.680});
    stages.push_back({(double)finalTarget/N0, 0.520});

    for (const Stage& st: stages){
        if (aliveVertices <= finalTarget || timeExceeded()) break;
        g_minDot = st.minDot;
        int stageTarget = max(finalTarget, (int)ceil(N0 * st.ratio));
        int stagnant = 0;
        while (aliveVertices > stageTarget && !timeExceeded()){
            rebuildIncOnly();
            if (++sweepId == INT_MAX) { fill(touched.begin(), touched.end(), 0); sweepId=1; }
            int before = aliveVertices;
            int step = 1000003 + 131071 * (sweepId % 97);
            if ((step & 1)==0) ++step;
            while (std::gcd(step, max(1, F0)) != 1) step += 2;
            int start = (int)((uint64_t)2654435761u * (uint32_t)sweepId % (uint32_t)max(1, F0));
            for (int iter=0; iter<F0 && aliveVertices>stageTarget; ++iter){
                if ((iter & 65535)==0 && timeExceeded()) break;
                int fid = (int)(((uint64_t)start + (uint64_t)iter * (uint64_t)step) % (uint64_t)F0);
                if (!F[fid].alive) continue;
                int a=F[fid].v[0], b=F[fid].v[1], c=F[fid].v[2];
                if (!V[a].alive || !V[b].alive || !V[c].alive) continue;
                struct E { int u,v; double l; } es[3] = {
                    {a,b,norm2(V[a].p-V[b].p)}, {b,c,norm2(V[b].p-V[c].p)}, {c,a,norm2(V[c].p-V[a].p)}
                };
                if (es[1].l < es[0].l) swap(es[0],es[1]);
                if (es[2].l < es[1].l) swap(es[1],es[2]);
                if (es[1].l < es[0].l) swap(es[0],es[1]);
                for (int ei=0; ei<3; ++ei){
                    int u=es[ei].u, v=es[ei].v;
                    if (!V[u].alive || !V[v].alive) continue;
                    if (touched[u]==sweepId || touched[v]==sweepId) continue;
                    EvalRes e = evaluateEdge(u,v);
                    if (!e.ok) continue;
                    collapseEdge(u,v,e.p);
                    if (V[u].alive) touched[u]=sweepId;
                    // Mark the immediate triangle vertices too; this keeps each sweep close to a matching.
                    if ((unsigned)a<(unsigned)N0 && V[a].alive) touched[a]=sweepId;
                    if ((unsigned)b<(unsigned)N0 && V[b].alive) touched[b]=sweepId;
                    if ((unsigned)c<(unsigned)N0 && V[c].alive) touched[c]=sweepId;
                    break;
                }
            }
            int removed = before - aliveVertices;
            if (removed < max(10, before/2000)) {
                if (++stagnant >= 2) break;
            } else stagnant = 0;
        }
    }
    rebuildIncOnly();
    g_queueActive = true;
}


static void greedySafeBonus(){
    if (aliveVertices <= 8 || timeExceeded()) return;
    bool oldQueue = g_queueActive;
    g_queueActive = false;
    priority_queue<EdgeItem> empty; PQ.swap(empty);
    int bonusTarget = max(4, (int)ceil(N0 * 0.003));
    // On ordinary smooth/detail meshes this strict dot limit stops quickly.  On exactly planar
    // tessellations it removes the large amount of redundancy that a ratio target would leave.
    g_minDot = 0.9997;
    vector<int> touched(N0, 0);
    int sweepId = 12345;
    int stagnant = 0;
    while (aliveVertices > bonusTarget && !timeExceeded() && stagnant < 3){
        rebuildIncOnly();
        if (++sweepId == INT_MAX) { fill(touched.begin(), touched.end(), 0); sweepId=1; }
        int before = aliveVertices;
        int step = 917513 + 65537 * (sweepId % 89);
        if ((step & 1)==0) ++step;
        while (std::gcd(step, max(1, F0)) != 1) step += 2;
        int start = (int)((uint64_t)2246822519u * (uint32_t)sweepId % (uint32_t)max(1, F0));
        for (int iter=0; iter<F0 && aliveVertices>bonusTarget; ++iter){
            if ((iter & 65535)==0 && timeExceeded()) break;
            int fid = (int)(((uint64_t)start + (uint64_t)iter * (uint64_t)step) % (uint64_t)F0);
            if (!F[fid].alive) continue;
            int a=F[fid].v[0], b=F[fid].v[1], c=F[fid].v[2];
            if (!V[a].alive || !V[b].alive || !V[c].alive) continue;
            struct E { int u,v; double l; } es[3] = {
                {a,b,norm2(V[a].p-V[b].p)}, {b,c,norm2(V[b].p-V[c].p)}, {c,a,norm2(V[c].p-V[a].p)}
            };
            if (es[1].l < es[0].l) swap(es[0],es[1]);
            if (es[2].l < es[1].l) swap(es[1],es[2]);
            if (es[1].l < es[0].l) swap(es[0],es[1]);
            for (int ei=0; ei<3; ++ei){
                int u=es[ei].u, v=es[ei].v;
                if (!V[u].alive || !V[v].alive) continue;
                if (touched[u]==sweepId || touched[v]==sweepId) continue;
                EvalRes e = evaluateEdge(u,v);
                if (!e.ok) continue;
                // Keep the bonus pass truly conservative: no high-cost skinny/large-radius surprises.
                if (e.key > 0.80) continue;
                collapseEdge(u,v,e.p);
                if (V[u].alive) touched[u]=sweepId;
                if ((unsigned)a<(unsigned)N0 && V[a].alive) touched[a]=sweepId;
                if ((unsigned)b<(unsigned)N0 && V[b].alive) touched[b]=sweepId;
                if ((unsigned)c<(unsigned)N0 && V[c].alive) touched[c]=sweepId;
                break;
            }
        }
        int removed = before - aliveVertices;
        if (removed < max(4, before/5000)) stagnant++; else stagnant=0;
    }
    rebuildIncOnly();
    g_queueActive = oldQueue;
}

static void simplify(){
    if (N0 <= 4) return;
    int finalTarget = chooseTarget();
    g_fastMode = (N0 > 250000);
    if (g_fastMode) {
        greedyLargeSimplify(finalTarget);
        greedySafeBonus();
        return;
    }

    struct Stage { double ratio; double minDot; };
    vector<Stage> stages;
    if (N0 <= 12) {
        stages.push_back({(double)finalTarget/N0, 0.05});
    } else {
        stages.push_back({max((double)finalTarget/N0, 0.300), 0.965});
        stages.push_back({max((double)finalTarget/N0, 0.200), 0.910});
        stages.push_back({max((double)finalTarget/N0, 0.135), 0.820});
        stages.push_back({max((double)finalTarget/N0, 0.102), 0.700});
        stages.push_back({(double)finalTarget/N0, 0.560});
    }

    for (size_t si=0; si<stages.size() && aliveVertices>finalTarget; ++si){
        if (timeExceeded()) break;
        g_minDot = stages[si].minDot;
        int stageTarget = max(finalTarget, (int)ceil(N0 * stages[si].ratio));
        rebuildQueue(true);
        long long pops=0, collapsesAtRebuild=0;
        int lastRebuildAlive = aliveVertices;
        while (aliveVertices > stageTarget && !PQ.empty()){
            if ((pops++ & 2047LL)==0 && timeExceeded()) break;
            EdgeItem it = PQ.top(); PQ.pop();
            if (it.u<0 || it.v<0 || it.u>=N0 || it.v>=N0) continue;
            if (!V[it.u].alive || !V[it.v].alive) continue;
            if (it.vu != V[it.u].ver || it.vv != V[it.v].ver){
                pushEdge(it.u, it.v);
                continue;
            }
            EvalRes e = evaluateEdge(it.u, it.v);
            if (!e.ok) continue;
            // If local cost has changed materially, reinsert once with the current key.
            if (e.key > it.key * 1.35 + 1e-9) {
                PQ.push(EdgeItem{e.key, it.u, it.v, V[it.u].ver, V[it.v].ver});
                continue;
            }
            collapseEdge(it.u, it.v, e.p);
            collapsesAtRebuild++;
            // Large stale heaps are a real time killer; rebuild after substantial progress.
            if ((N0 > 100000 && collapsesAtRebuild >= 70000) || PQ.size() > (size_t)max(2000000, aliveFaces*8)){
                if (aliveVertices <= stageTarget || timeExceeded()) break;
                rebuildQueue(true);
                collapsesAtRebuild=0;
                lastRebuildAlive = aliveVertices;
                if (PQ.empty() || lastRebuildAlive <= stageTarget) break;
            }
        }
    }
    greedySafeBonus();
}

static bool finalSanitySmall(){
    // Full check on the final mesh. If this fails, we still print current mesh; this function is mainly
    // a guard against accidental zero-area faces after floating-point rounding in development.
    for (int i=0;i<F0;i++) if (F[i].alive){
        int a=F[i].v[0], b=F[i].v[1], c=F[i].v[2];
        if (a==b || b==c || a==c || !V[a].alive || !V[b].alive || !V[c].alive) return false;
        Vec3 cr = crossv(V[b].p-V[a].p, V[c].p-V[a].p);
        if (!(normv(cr) > 1e-18)) return false;
    }
    return true;
}

static void saveMesh(){
    finalSanitySmall();
    vector<int> id(N0, -1);
    int nv=0, nf=0;
    for (int i=0;i<N0;i++) if (V[i].alive) id[i]=++nv;
    for (int i=0;i<F0;i++) if (F[i].alive){
        int a=F[i].v[0], b=F[i].v[1], c=F[i].v[2];
        if (a!=b && b!=c && a!=c && V[a].alive && V[b].alive && V[c].alive) nf++;
    }
    string out;
    out.reserve(min<size_t>((size_t)120000000, (size_t)nv*42 + (size_t)nf*24 + 64));
    char line[128];
    int len = snprintf(line, sizeof(line), "%d %d\n", nv, nf);
    out.append(line, line+len);
    for (int i=0;i<N0;i++) if (V[i].alive){
        len = snprintf(line, sizeof(line), "v %.10g %.10g %.10g\n", V[i].p.x, V[i].p.y, V[i].p.z);
        out.append(line, line+len);
        if (out.size() > (1<<20)) { fwrite(out.data(),1,out.size(),stdout); out.clear(); }
    }
    for (int i=0;i<F0;i++) if (F[i].alive){
        int a=F[i].v[0], b=F[i].v[1], c=F[i].v[2];
        if (a==b || b==c || a==c || !V[a].alive || !V[b].alive || !V[c].alive) continue;
        len = snprintf(line, sizeof(line), "f %d %d %d\n", id[a], id[b], id[c]);
        out.append(line, line+len);
        if (out.size() > (1<<20)) { fwrite(out.data(),1,out.size(),stdout); out.clear(); }
    }
    if (!out.empty()) fwrite(out.data(),1,out.size(),stdout);
}

int main(){
    g_start = chrono::steady_clock::now();
    loadMesh();
    g_detailP90 = estimateDetailP90();
    simplify();
    saveMesh();
    return 0;
}
