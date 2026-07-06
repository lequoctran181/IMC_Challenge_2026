#include <bits/stdc++.h>
using namespace std;

struct FastInput {
    static const size_t BUFSIZE = 1 << 20;
    int idx = 0, size = 0;
    char buf[BUFSIZE];
    inline char getch() {
        if (idx >= size) {
            size = (int)fread(buf, 1, BUFSIZE, stdin);
            idx = 0;
            if (size == 0) return 0;
        }
        return buf[idx++];
    }
    inline void skipws() {
        char c;
        do { c = getch(); } while (c && (c==' ' || c=='\n' || c=='\r' || c=='\t'));
        if (c) idx--;
    }
    int nextInt() {
        skipws();
        int sgn = 1, x = 0;
        char c = getch();
        if (c == '-') { sgn = -1; c = getch(); }
        while (c >= '0' && c <= '9') {
            x = x * 10 + (c - '0');
            c = getch();
        }
        return x * sgn;
    }
    double nextDouble() {
        skipws();
        char tmp[64]; int n = 0;
        char c = getch();
        while (c && !(c==' ' || c=='\n' || c=='\r' || c=='\t')) {
            if (n < 63) tmp[n++] = c;
            c = getch();
        }
        tmp[n] = 0;
        return strtod(tmp, nullptr);
    }
    char nextChar() { skipws(); return getch(); }
};

struct Vec3 {
    double x, y, z;
    Vec3() : x(0), y(0), z(0) {}
    Vec3(double X, double Y, double Z) : x(X), y(Y), z(Z) {}
};
static inline Vec3 operator+(const Vec3& a, const Vec3& b){ return Vec3(a.x+b.x,a.y+b.y,a.z+b.z); }
static inline Vec3 operator-(const Vec3& a, const Vec3& b){ return Vec3(a.x-b.x,a.y-b.y,a.z-b.z); }
static inline Vec3 operator*(const Vec3& a, double s){ return Vec3(a.x*s,a.y*s,a.z*s); }
static inline double dot3(const Vec3& a, const Vec3& b){ return a.x*b.x+a.y*b.y+a.z*b.z; }
static inline Vec3 cross3(const Vec3& a, const Vec3& b){ return Vec3(a.y*b.z-a.z*b.y, a.z*b.x-a.x*b.z, a.x*b.y-a.y*b.x); }
static inline double norm2(const Vec3& a){ return dot3(a,a); }
static inline double dist2v(const Vec3& a, const Vec3& b){ return norm2(a-b); }

struct Quadric {
    double q[10];
    Quadric() { memset(q, 0, sizeof(q)); }
    inline void addPlane(double a, double b, double c, double d, double w) {
        q[0] += w*a*a; q[1] += w*a*b; q[2] += w*a*c; q[3] += w*a*d;
        q[4] += w*b*b; q[5] += w*b*c; q[6] += w*b*d;
        q[7] += w*c*c; q[8] += w*c*d; q[9] += w*d*d;
    }
    inline void add(const Quadric& o) { for (int i=0;i<10;i++) q[i] += o.q[i]; }
    inline double eval(const Vec3& p) const {
        double x=p.x,y=p.y,z=p.z;
        double v = q[0]*x*x + 2*q[1]*x*y + 2*q[2]*x*z + 2*q[3]*x
                 + q[4]*y*y + 2*q[5]*y*z + 2*q[6]*y
                 + q[7]*z*z + 2*q[8]*z + q[9];
        return v < 0 && v > -1e-18 ? 0 : v;
    }
};

struct Face {
    int v[3];
    unsigned char alive;
};

struct Candidate {
    float cost;
    int u, v;
    int vu, vv;
};
struct CandCmp {
    bool operator()(const Candidate& a, const Candidate& b) const {
        return a.cost > b.cost;
    }
};

static int N, M;
static vector<Vec3> Orig, P;
static vector<Face> F;
static vector<int> headNode, nextNode;
static vector<unsigned char> vAlive;
static vector<int> ver;
static vector<Quadric> Q;
static vector<int> memHead, memTail, memNext, memSize;
static vector<double> rad2;
static vector<int> markA, markB;
static int stampA = 1, stampB = 1000000007;
static int aliveV, aliveF;
static double epsDist, eps2Dist;
static chrono::steady_clock::time_point tStart;

static inline double elapsedSec() {
    return chrono::duration<double>(chrono::steady_clock::now() - tStart).count();
}

static inline uint64_t edgeKey(int a, int b) {
    if (a > b) swap(a,b);
    return (uint64_t)(uint32_t)a << 32 | (uint32_t)b;
}

static inline bool faceHas(const Face& f, int x) {
    return f.v[0] == x || f.v[1] == x || f.v[2] == x;
}
static inline int thirdOf(const Face& f, int a, int b) {
    for (int i=0;i<3;i++) if (f.v[i] != a && f.v[i] != b) return f.v[i];
    return -1;
}

static void cleanList(int u) {
    if (!vAlive[u]) return;
    int nh = -1;
    for (int node = headNode[u]; node != -1; ) {
        int nx = nextNode[node];
        int fid = node / 3, k = node - fid*3;
        if (F[fid].alive && F[fid].v[k] == u) {
            nextNode[node] = nh;
            nh = node;
        }
        node = nx;
    }
    headNode[u] = nh;
}

static inline double faceArea2(int a, int b, int c) {
    return norm2(cross3(P[b]-P[a], P[c]-P[a]));
}

static bool radiusOK(int from, int to, double& nr2) {
    nr2 = rad2[to];
    const Vec3& p = P[to];
    const double lim = eps2Dist * (1.0 + 2e-12);
    for (int m = memHead[from]; m != -1; m = memNext[m]) {
        double d2 = dist2v(Orig[m], p);
        if (d2 > lim) return false;
        if (d2 > nr2) nr2 = d2;
    }
    return true;
}

static inline Vec3 posInFace(int idx, int replFrom, int replTo, const Vec3& newPos) {
    if (idx == replFrom) idx = replTo;
    if (idx == replTo) return newPos;
    return P[idx];
}

static double normalPenaltyNew(int from, int to, const Vec3& newPos) {
    const double tiny = 1e-32;
    double pen = 0.0;
    bool toMoves = dist2v(newPos, P[to]) > 1e-28;
    for (int node = headNode[from]; node != -1; node = nextNode[node]) {
        int fid = node / 3, k = node - fid*3;
        Face &f = F[fid];
        if (!f.alive || f.v[k] != from) continue;
        if (faceHas(f, to)) continue;
        int a=f.v[0], b=f.v[1], c=f.v[2];
        Vec3 oldn = cross3(P[b]-P[a], P[c]-P[a]);
        double oldl2 = norm2(oldn);
        if (oldl2 <= tiny) continue;
        Vec3 pa = posInFace(a, from, to, newPos);
        Vec3 pb = posInFace(b, from, to, newPos);
        Vec3 pc = posInFace(c, from, to, newPos);
        Vec3 newn = cross3(pb-pa, pc-pa);
        double newl2 = norm2(newn);
        if (newl2 <= tiny) return 1e100;
        double d = dot3(oldn, newn);
        if (d <= 1e-14 * sqrt(max(oldl2, tiny) * newl2)) return 1e100;
        double co = d / sqrt(oldl2 * newl2);
        if (co > 1.0) co = 1.0;
        pen += sqrt(oldl2) * (1.0 - co);
    }
    if (toMoves) {
        for (int node = headNode[to]; node != -1; node = nextNode[node]) {
            int fid = node / 3, k = node - fid*3;
            Face &f = F[fid];
            if (!f.alive || f.v[k] != to) continue;
            if (faceHas(f, from)) continue;
            int a=f.v[0], b=f.v[1], c=f.v[2];
            Vec3 oldn = cross3(P[b]-P[a], P[c]-P[a]);
            double oldl2 = norm2(oldn);
            if (oldl2 <= tiny) continue;
            Vec3 pa = (a == to ? newPos : P[a]);
            Vec3 pb = (b == to ? newPos : P[b]);
            Vec3 pc = (c == to ? newPos : P[c]);
            Vec3 newn = cross3(pb-pa, pc-pa);
            double newl2 = norm2(newn);
            if (newl2 <= tiny) return 1e100;
            double d = dot3(oldn, newn);
            if (d <= 1e-14 * sqrt(max(oldl2, tiny) * newl2)) return 1e100;
            double co = d / sqrt(oldl2 * newl2);
            if (co > 1.0) co = 1.0;
            pen += sqrt(oldl2) * (1.0 - co);
        }
    }
    return pen;
}

static bool radiusOKPoint(int a, int b, const Vec3& p, double& nr2) {
    nr2 = 0.0;
    const double lim = eps2Dist * (1.0 + 2e-12);
    for (int m = memHead[a]; m != -1; m = memNext[m]) {
        double d2 = dist2v(Orig[m], p);
        if (d2 > lim) return false;
        if (d2 > nr2) nr2 = d2;
    }
    for (int m = memHead[b]; m != -1; m = memNext[m]) {
        double d2 = dist2v(Orig[m], p);
        if (d2 > lim) return false;
        if (d2 > nr2) nr2 = d2;
    }
    return true;
}

static bool solveOptimal(const Quadric& q, Vec3& p) {
    double a00=q.q[0], a01=q.q[1], a02=q.q[2];
    double a11=q.q[4], a12=q.q[5], a22=q.q[7];
    double b0=-q.q[3], b1=-q.q[6], b2=-q.q[8];
    double det = a00*(a11*a22-a12*a12) - a01*(a01*a22-a12*a02) + a02*(a01*a12-a11*a02);
    if (fabs(det) < 1e-18) return false;
    double c00=(a11*a22-a12*a12), c01=-(a01*a22-a12*a02), c02=(a01*a12-a11*a02);
    double c10=c01, c11=(a00*a22-a02*a02), c12=-(a00*a12-a01*a02);
    double c20=c02, c21=c12, c22=(a00*a11-a01*a01);
    p.x = (c00*b0 + c01*b1 + c02*b2) / det;
    p.y = (c10*b0 + c11*b1 + c12*b2) / det;
    p.z = (c20*b0 + c21*b1 + c22*b2) / det;
    return isfinite(p.x) && isfinite(p.y) && isfinite(p.z);
}

static bool chooseDirection(int u, int v, int& from, int& to, Vec3& newPos, double& newR2) {
    Quadric q = Q[u]; q.add(Q[v]);
    struct Opt { double c; int from, to; Vec3 p; double r2; } opt[7];
    int cnt = 0;
    double r2a;
    if (radiusOK(v, u, r2a)) {
        Vec3 p = P[u];
        double pen = normalPenaltyNew(v, u, p);
        if (isfinite(pen) && pen < 1e90) opt[cnt++] = {q.eval(p) + 0.02 * pen, v, u, p, r2a};
    }
    double r2b;
    if (radiusOK(u, v, r2b)) {
        Vec3 p = P[v];
        double pen = normalPenaltyNew(u, v, p);
        if (isfinite(pen) && pen < 1e90) opt[cnt++] = {q.eval(p) + 0.02 * pen, u, v, p, r2b};
    }
    Vec3 mid = (P[u] + P[v]) * 0.5;
    double r2m;
    if (radiusOKPoint(u, v, mid, r2m)) {
        int keep = (memSize[u] >= memSize[v] ? u : v);
        int rem = (keep == u ? v : u);
        double pen = normalPenaltyNew(rem, keep, mid);
        if (isfinite(pen) && pen < 1e90) opt[cnt++] = {q.eval(mid) + 0.02 * pen + 1e-13, rem, keep, mid, r2m};
        pen = normalPenaltyNew(keep, rem, mid);
        if (isfinite(pen) && pen < 1e90) opt[cnt++] = {q.eval(mid) + 0.02 * pen + 1e-13, keep, rem, mid, r2m};
    }
    Vec3 op;
    if (solveOptimal(q, op)) {
        double r2o;
        if (radiusOKPoint(u, v, op, r2o)) {
            int keep = (memSize[u] >= memSize[v] ? u : v);
            int rem = (keep == u ? v : u);
            double pen = normalPenaltyNew(rem, keep, op);
            if (isfinite(pen) && pen < 1e90) opt[cnt++] = {q.eval(op) + 0.02 * pen, rem, keep, op, r2o};
            pen = normalPenaltyNew(keep, rem, op);
            if (isfinite(pen) && pen < 1e90) opt[cnt++] = {q.eval(op) + 0.02 * pen, keep, rem, op, r2o};
        }
    }
    if (cnt == 0) return false;
    sort(opt, opt + cnt, [](const Opt& a, const Opt& b){ return a.c < b.c; });
    for (int i=0;i<cnt;i++) {
        if (isfinite(opt[i].c)) {
            from = opt[i].from; to = opt[i].to; newPos = opt[i].p; newR2 = opt[i].r2;
            return true;
        }
    }
    return false;
}

static bool linkOK(int u, int v) {
    cleanList(u); cleanList(v);
    if (++stampA == 0) { fill(markA.begin(), markA.end(), 0); stampA = 1; }
    int st = stampA;
    int edgeCnt = 0, opp[3] = {-1,-1,-1};
    for (int node = headNode[u]; node != -1; node = nextNode[node]) {
        int fid = node / 3, k = node - fid*3;
        Face &f = F[fid];
        if (!f.alive || f.v[k] != u) continue;
        bool hasV = false;
        for (int i=0;i<3;i++) if (f.v[i] == v) hasV = true;
        if (hasV) {
            if (edgeCnt < 3) opp[edgeCnt] = thirdOf(f, u, v);
            edgeCnt++;
            if (edgeCnt > 2) return false;
        }
        for (int i=0;i<3;i++) {
            int w = f.v[i];
            if (w != u) markA[w] = st;
        }
    }
    if (edgeCnt != 2 || opp[0] < 0 || opp[1] < 0 || opp[0] == opp[1]) return false;
    if (++stampB == 0) { fill(markB.begin(), markB.end(), 0); stampB = 1000000007; }
    int st2 = stampB;
    int common = 0;
    for (int node = headNode[v]; node != -1; node = nextNode[node]) {
        int fid = node / 3, k = node - fid*3;
        Face &f = F[fid];
        if (!f.alive || f.v[k] != v) continue;
        for (int i=0;i<3;i++) {
            int w = f.v[i];
            if (w == v) continue;
            if (markA[w] == st && markB[w] != st2) {
                markB[w] = st2;
                common++;
                if (common > 2) return false;
            }
        }
    }
    if (common != 2) return false;
    return markB[opp[0]] == st2 && markB[opp[1]] == st2;
}

static double edgeCostRaw(int u, int v) {
    if (!vAlive[u] || !vAlive[v] || u == v) return 1e100;
    double l2 = dist2v(P[u], P[v]);
    if (l2 > 4.0 * eps2Dist * (1.0 + 1e-10)) return 1e100;
    Quadric q = Q[u]; q.add(Q[v]);
    double cu = q.eval(P[u]);
    double cv = q.eval(P[v]);
    double c = min(cu, cv);
    if (!isfinite(c)) return 1e100;
    if (c < 0) c = 0;
    c += 1e-12 * l2;
    return c;
}

static Candidate makeCand(int u, int v) {
    if (u > v) swap(u,v);
    double c = edgeCostRaw(u,v);
    Candidate r;
    r.cost = (float)min(c, 3.3e38);
    r.u = u; r.v = v; r.vu = ver[u]; r.vv = ver[v];
    return r;
}

static void pushIncident(int u, priority_queue<Candidate, vector<Candidate>, CandCmp>& pq) {
    if (!vAlive[u]) return;
    cleanList(u);
    if (++stampA == 0) { fill(markA.begin(), markA.end(), 0); stampA = 1; }
    int st = stampA;
    for (int node = headNode[u]; node != -1; node = nextNode[node]) {
        int fid = node / 3, k = node - fid*3;
        Face &f = F[fid];
        if (!f.alive || f.v[k] != u) continue;
        for (int i=0;i<3;i++) {
            int w = f.v[i];
            if (w == u || !vAlive[w]) continue;
            if (markA[w] == st) continue;
            markA[w] = st;
            Candidate c = makeCand(u,w);
            if (isfinite((double)c.cost) && c.cost < 3e38f) pq.push(c);
        }
    }
}

static bool collapseEdge(int u, int v, priority_queue<Candidate, vector<Candidate>, CandCmp>& pq) {
    if (aliveV <= 4) return false;
    if (!vAlive[u] || !vAlive[v] || u == v) return false;
    if (dist2v(P[u],P[v]) > 4.0 * eps2Dist * (1.0 + 1e-10)) return false;
    if (!linkOK(u,v)) return false;
    int from=-1, to=-1; double nr2=0; Vec3 chosenPos;
    if (!chooseDirection(u,v,from,to,chosenPos,nr2)) return false;

    vAlive[from] = 0;
    aliveV--;
    ver[from]++;
    ver[to]++;
    Q[to].add(Q[from]);
    P[to] = chosenPos;
    rad2[to] = nr2;
    if (memHead[from] != -1) {
        memNext[memTail[to]] = memHead[from];
        memTail[to] = memTail[from];
        memSize[to] += memSize[from];
    }
    memHead[from] = memTail[from] = -1;
    memSize[from] = 0;

    int node = headNode[from];
    while (node != -1) {
        int nx = nextNode[node];
        int fid = node / 3, k = node - fid*3;
        Face &f = F[fid];
        if (f.alive && f.v[k] == from) {
            if (faceHas(f, to)) {
                f.alive = 0;
                aliveF--;
            } else {
                f.v[k] = to;
                nextNode[node] = headNode[to];
                headNode[to] = node;
            }
        }
        node = nx;
    }
    headNode[from] = -1;
    cleanList(to);
    pushIncident(to, pq);
    return true;
}

struct Snapshot {
    vector<Face> faces;
    vector<unsigned char> alive;
    vector<Vec3> pos;
    int av = 0, af = 0;
    bool valid = false;
};
static Snapshot bestSnap;
static bool outputSnap = false;

static void takeSnapshot(Snapshot& s) {
    s.faces = F;
    s.alive = vAlive;
    s.pos = P;
    s.av = aliveV;
    s.af = aliveF;
    s.valid = true;
}

struct RenderMaps {
    int R = 0;
    vector<float> depth;
    vector<float> normal;
    vector<unsigned char> fg;
};

struct ViewDef { int axis; int sign; };

static inline void cameraCoords(const Vec3& p, int view, float& sx, float& sy, float& dep, int R) {
    const double D = 2.5;
    double x=0,y=0,d=0;
    switch(view) {
        case 0: d = D - p.x; x = p.y;  y = p.z; break;
        case 1: d = D + p.x; x = -p.y; y = p.z; break;
        case 2: d = D - p.y; x = -p.x; y = p.z; break;
        case 3: d = D + p.y; x = p.x;  y = p.z; break;
        case 4: d = D - p.z; x = p.x;  y = p.y; break;
        default:d = D + p.z; x = -p.x; y = p.y; break;
    }
    double f = 800.0 * (double)R / 1024.0;
    sx = (float)(f * x / d + 0.5 * R);
    sy = (float)(f * y / d + 0.5 * R);
    dep = (float)d;
}

static inline Vec3 normalToCamera(const Vec3& n, int view) {
    switch(view) {
        case 0: return Vec3(n.y, n.z, n.x);
        case 1: return Vec3(-n.y, n.z, -n.x);
        case 2: return Vec3(-n.x, n.z, n.y);
        case 3: return Vec3(n.x, n.z, -n.y);
        case 4: return Vec3(n.x, n.y, n.z);
        default:return Vec3(-n.x, n.y, -n.z);
    }
}

static inline float edge2d(float ax,float ay,float bx,float by,float cx,float cy) {
    return (cx-ax)*(by-ay) - (cy-ay)*(bx-ax);
}

static RenderMaps renderMesh(const vector<Face>& faces, int R) {
    RenderMaps out;
    out.R = R;
    int pix = R*R, total = 6*pix;
    out.depth.assign(total, 255.0f);
    out.normal.assign((size_t)total*3, 127.5f);
    out.fg.assign(total, 0);
    vector<float> sx(N), sy(N), dep(N);
    for (int view=0; view<6; ++view) {
        for (int i=0;i<N;i++) cameraCoords(P[i], view, sx[i], sy[i], dep[i], R);
        int base = view * pix;
        for (int fid=0; fid<M; ++fid) {
            const Face& f = faces[fid];
            if (!f.alive) continue;
            int ia=f.v[0], ib=f.v[1], ic=f.v[2];
            float x0=sx[ia], y0=sy[ia], z0=dep[ia];
            float x1=sx[ib], y1=sy[ib], z1=dep[ib];
            float x2=sx[ic], y2=sy[ic], z2=dep[ic];
            float area = edge2d(x0,y0,x1,y1,x2,y2);
            if (fabs(area) < 1e-12f) continue;
            float minxf = floorf(min(x0, min(x1,x2)) - 0.5f);
            float maxxf = ceilf (max(x0, max(x1,x2)) - 0.5f);
            float minyf = floorf(min(y0, min(y1,y2)) - 0.5f);
            float maxyf = ceilf (max(y0, max(y1,y2)) - 0.5f);
            int xmin = max(0, (int)minxf), xmax = min(R-1, (int)maxxf);
            int ymin = max(0, (int)minyf), ymax = min(R-1, (int)maxyf);
            if (xmin > xmax || ymin > ymax) continue;
            Vec3 wn = cross3(P[ib]-P[ia], P[ic]-P[ia]);
            double nl = sqrt(norm2(wn));
            if (nl <= 1e-30) continue;
            wn = wn * (1.0 / nl);
            Vec3 cn = normalToCamera(wn, view);
            float nr = (float)((cn.x + 1.0) * 127.5);
            float ng = (float)((cn.y + 1.0) * 127.5);
            float nb = (float)((cn.z + 1.0) * 127.5);
            float invArea = 1.0f / area;
            float iz0 = 1.0f / z0, iz1 = 1.0f / z1, iz2 = 1.0f / z2;
            for (int yy=ymin; yy<=ymax; ++yy) {
                float py = yy + 0.5f;
                int row = base + yy*R;
                for (int xx=xmin; xx<=xmax; ++xx) {
                    float px = xx + 0.5f;
                    float w0 = edge2d(x1,y1,x2,y2,px,py) * invArea;
                    float w1 = edge2d(x2,y2,x0,y0,px,py) * invArea;
                    float w2 = 1.0f - w0 - w1;
                    if (w0 >= -1e-6f && w1 >= -1e-6f && w2 >= -1e-6f) {
                        float zz = 1.0f / (w0*iz0 + w1*iz1 + w2*iz2);
                        int id = row + xx;
                        if (zz < out.depth[id]) {
                            out.depth[id] = zz;
                            out.fg[id] = 1;
                            size_t ni = (size_t)id * 3;
                            out.normal[ni] = nr; out.normal[ni+1] = ng; out.normal[ni+2] = nb;
                        }
                    }
                }
            }
        }
    }
    return out;
}

static double ssimChannel(const RenderMaps& A, const RenderMaps& B, int view, int channel) {
    int R = A.R;
    int pix = R*R;
    int base = view * pix;
    int W = R + 1;
    vector<double> ix((size_t)W*W), iy((size_t)W*W), ix2((size_t)W*W), iy2((size_t)W*W), ixy((size_t)W*W);
    for (int y=0; y<R; ++y) {
        double sxr=0, syr=0, sx2r=0, sy2r=0, sxyr=0;
        for (int x=0; x<R; ++x) {
            int id = base + y*R + x;
            double xv, yv;
            if (channel == 0) { xv = A.depth[id]; yv = B.depth[id]; }
            else {
                size_t ni = (size_t)id*3 + (channel-1);
                xv = A.normal[ni]; yv = B.normal[ni];
            }
            sxr += xv; syr += yv; sx2r += xv*xv; sy2r += yv*yv; sxyr += xv*yv;
            size_t p = (size_t)(y+1)*W + (x+1);
            size_t up = (size_t)y*W + (x+1);
            ix[p] = ix[up] + sxr;
            iy[p] = iy[up] + syr;
            ix2[p] = ix2[up] + sx2r;
            iy2[p] = iy2[up] + sy2r;
            ixy[p] = ixy[up] + sxyr;
        }
    }
    auto rect = [W](const vector<double>& ii, int x0, int y0, int x1, int y1)->double {
        return ii[(size_t)y1*W + x1] - ii[(size_t)y0*W + x1] - ii[(size_t)y1*W + x0] + ii[(size_t)y0*W + x0];
    };
    const double C1 = (0.01*255.0)*(0.01*255.0);
    const double C2 = (0.03*255.0)*(0.03*255.0);
    double sum = 0.0; long long cnt = 0;
    for (int y=0; y<R; ++y) {
        int y0 = max(0, y-5), y1 = min(R, y+6);
        for (int x=0; x<R; ++x) {
            int id = base + y*R + x;
            if (!A.fg[id] && !B.fg[id]) continue;
            int x0 = max(0, x-5), x1 = min(R, x+6);
            double n = (double)(x1-x0)*(y1-y0);
            double sx = rect(ix,x0,y0,x1,y1), sy = rect(iy,x0,y0,x1,y1);
            double sx2 = rect(ix2,x0,y0,x1,y1), sy2 = rect(iy2,x0,y0,x1,y1), sxy = rect(ixy,x0,y0,x1,y1);
            double mux = sx/n, muy = sy/n;
            double vx = max(0.0, sx2/n - mux*mux);
            double vy = max(0.0, sy2/n - muy*muy);
            double cov = sxy/n - mux*muy;
            double num = (2*mux*muy + C1) * (2*cov + C2);
            double den = (mux*mux + muy*muy + C1) * (vx + vy + C2);
            sum += den != 0.0 ? num/den : 1.0;
            cnt++;
        }
    }
    if (cnt == 0) return 1.0;
    return sum / (double)cnt;
}

static double approxSSIM(const RenderMaps& A, const RenderMaps& B) {
    double total = 0.0;
    for (int v=0; v<6; ++v) {
        double dn = ssimChannel(A,B,v,0);
        double n0 = ssimChannel(A,B,v,1);
        double n1 = ssimChannel(A,B,v,2);
        double n2 = ssimChannel(A,B,v,3);
        double ns = (n0+n1+n2)/3.0;
        total += 0.5*dn + 0.5*ns;
    }
    return total / 6.0;
}

static void buildInitial() {
    headNode.assign(N, -1);
    nextNode.assign((size_t)3*M, -1);
    for (int fid=0; fid<M; ++fid) {
        for (int k=0;k<3;k++) {
            int node = 3*fid + k;
            int v = F[fid].v[k];
            nextNode[node] = headNode[v];
            headNode[v] = node;
        }
    }
    Q.assign(N, Quadric());
    for (int i=0;i<M;i++) {
        int a=F[i].v[0], b=F[i].v[1], c=F[i].v[2];
        Vec3 n = cross3(P[b]-P[a], P[c]-P[a]);
        double len = sqrt(norm2(n));
        if (len <= 1e-30) continue;
        double area = 0.5 * len;
        n = n * (1.0 / len);
        double d = -dot3(n, P[a]);
        double w = area;
        Q[a].addPlane(n.x,n.y,n.z,d,w);
        Q[b].addPlane(n.x,n.y,n.z,d,w);
        Q[c].addPlane(n.x,n.y,n.z,d,w);
    }
    vAlive.assign(N, 1);
    ver.assign(N, 0);
    memHead.resize(N); memTail.resize(N); memNext.assign(N, -1); memSize.assign(N, 1); rad2.assign(N, 0.0);
    for (int i=0;i<N;i++) memHead[i] = memTail[i] = i;
    markA.assign(N, 0); markB.assign(N, 0);
    aliveV = N; aliveF = M;
}

static priority_queue<Candidate, vector<Candidate>, CandCmp> buildQueue() {
    vector<uint64_t> keys;
    keys.reserve((size_t)3*M);
    for (int i=0;i<M;i++) {
        int a=F[i].v[0], b=F[i].v[1], c=F[i].v[2];
        keys.push_back(edgeKey(a,b)); keys.push_back(edgeKey(b,c)); keys.push_back(edgeKey(c,a));
    }
    sort(keys.begin(), keys.end());
    vector<Candidate> cs;
    cs.reserve(keys.size()/2 + 1);
    uint64_t last = UINT64_MAX;
    for (uint64_t k : keys) {
        if (k == last) continue;
        last = k;
        int a = (int)(k >> 32), b = (int)(uint32_t)k;
        Candidate c = makeCand(a,b);
        if (isfinite((double)c.cost) && c.cost < 3e38f) cs.push_back(c);
    }
    vector<uint64_t>().swap(keys);
    priority_queue<Candidate, vector<Candidate>, CandCmp> pq(CandCmp(), std::move(cs));
    return pq;
}

static void simplifyMesh() {
    if (N <= 4) {
        vAlive.assign(N, 1);
        aliveV = N;
        aliveF = M;
        return;
    }
    double xmin=Orig[0].x, xmax=Orig[0].x, ymin=Orig[0].y, ymax=Orig[0].y, zmin=Orig[0].z, zmax=Orig[0].z;
    for (int i=1;i<N;i++) {
        xmin=min(xmin,Orig[i].x); xmax=max(xmax,Orig[i].x);
        ymin=min(ymin,Orig[i].y); ymax=max(ymax,Orig[i].y);
        zmin=min(zmin,Orig[i].z); zmax=max(zmax,Orig[i].z);
    }
    double diag = sqrt((xmax-xmin)*(xmax-xmin)+(ymax-ymin)*(ymax-ymin)+(zmax-zmin)*(zmax-zmin));
    epsDist = 0.05 * diag * 0.999;
    eps2Dist = epsDist * epsDist;

    buildInitial();

    bool useEval = (N >= 2000);
    int R = 0;
    RenderMaps originalMap;
    if (useEval) {
        R = (N <= 7000 ? 1024 : (N < 30000 ? 640 : (N < 70000 ? 512 : (N < 300000 ? 384 : 256))));
        if (elapsedSec() < 3.0) originalMap = renderMesh(F, R);
        else useEval = false;
        if (useEval) takeSnapshot(bestSnap);
    }

    auto pq = buildQueue();

    vector<double> checks;
    if (useEval) {
        double arr[] = {0.20,0.15,0.12,0.105,0.095,0.090,0.086,0.083,0.080,0.077,0.074,0.071,0.068,0.064,0.060,0.056,0.052,0.048,0.044,0.040,0.036,0.032,0.028,0.024,0.020};
        checks.assign(arr, arr + sizeof(arr)/sizeof(arr[0]));
    }
    int checkIdx = 0;
    int nextCheck = useEval ? max(4, (int)ceil(N * checks[0])) : -1;
    double fixedRatio = 0.078;
    if (N < 1000) fixedRatio = 0.04;
    int hardTarget = max(4, (int)ceil(N * (useEval ? 0.020 : fixedRatio)));
    const double timeStop = 19.2;
    int popCounter = 0;
    bool stoppedBySSIM = false;

    while (aliveV > hardTarget && !pq.empty()) {
        if ((++popCounter & 4095) == 0 && elapsedSec() > timeStop) break;
        Candidate c = pq.top(); pq.pop();
        int u=c.u, v=c.v;
        if (u<0 || v<0 || u>=N || v>=N) continue;
        if (!vAlive[u] || !vAlive[v]) continue;
        if (ver[u] != c.vu || ver[v] != c.vv) continue;
        collapseEdge(u,v,pq);

        if (useEval && aliveV <= nextCheck) {
            if (elapsedSec() > 18.0) break;
            RenderMaps cur = renderMesh(F, R);
            double sc = approxSSIM(originalMap, cur);
            double ratio = (double)aliveV / (double)N;
            double threshold;
            if (R >= 1024) threshold = 0.9004;
            else if (R >= 640) threshold = (ratio > 0.080 ? 0.9030 : 0.9050);
            else if (R >= 512) threshold = (ratio > 0.080 ? 0.9050 : 0.9075);
            else if (R >= 384) threshold = (ratio > 0.080 ? 0.9070 : 0.9100);
            else threshold = (ratio > 0.080 ? 0.9090 : 0.9120);
            if (sc >= threshold) {
                takeSnapshot(bestSnap);
                checkIdx++;
                if (checkIdx >= (int)checks.size()) {
                    nextCheck = hardTarget;
                } else {
                    nextCheck = max(4, (int)ceil(N * checks[checkIdx]));
                    if (nextCheck >= aliveV) nextCheck = aliveV - 1;
                }
            } else {
                stoppedBySSIM = true;
                break;
            }
        }
    }
    if (useEval && bestSnap.valid) {
        outputSnap = true;
    } else {
        outputSnap = false;
    }
    (void)stoppedBySSIM;
}

static void saveMesh(const vector<Face>& faces, const vector<unsigned char>& aliveFlags, const vector<Vec3>& pos) {
    vector<int> mapIdx(N, -1);
    vector<unsigned char> used(N, 0);
    int nf = 0;
    for (int i=0;i<M;i++) if (faces[i].alive) {
        int a=faces[i].v[0], b=faces[i].v[1], c=faces[i].v[2];
        if (a==b || b==c || a==c) continue;
        used[a]=used[b]=used[c]=1;
        nf++;
    }
    int nv = 0;
    for (int i=0;i<N;i++) if (used[i] && aliveFlags[i]) mapIdx[i] = ++nv;
    if (nv < 4 || nf < 4) {
        printf("%d %d\n", N, M);
        for (int i=0;i<N;i++) printf("v %.10g %.10g %.10g\n", pos[i].x, pos[i].y, pos[i].z);
        for (int i=0;i<M;i++) printf("f %d %d %d\n", F[i].v[0]+1, F[i].v[1]+1, F[i].v[2]+1);
        return;
    }
    string out;
    out.reserve((size_t)nv * 44 + (size_t)nf * 28 + 64);
    char line[128];
    int len = snprintf(line, sizeof(line), "%d %d\n", nv, nf);
    out.append(line, len);
    for (int i=0;i<N;i++) if (mapIdx[i] > 0) {
        len = snprintf(line, sizeof(line), "v %.10g %.10g %.10g\n", pos[i].x, pos[i].y, pos[i].z);
        out.append(line, len);
    }
    for (int i=0;i<M;i++) if (faces[i].alive) {
        int a=faces[i].v[0], b=faces[i].v[1], c=faces[i].v[2];
        if (a==b || b==c || a==c) continue;
        int ia=mapIdx[a], ib=mapIdx[b], ic=mapIdx[c];
        if (ia<=0 || ib<=0 || ic<=0) continue;
        len = snprintf(line, sizeof(line), "f %d %d %d\n", ia, ib, ic);
        out.append(line, len);
    }
    fwrite(out.data(), 1, out.size(), stdout);
}

int main() {
    tStart = chrono::steady_clock::now();
    FastInput in;
    N = in.nextInt(); M = in.nextInt();
    Orig.resize(N); P.resize(N);
    F.resize(M);
    for (int i=0;i<N;i++) {
        char c = in.nextChar(); (void)c;
        double x = in.nextDouble(), y = in.nextDouble(), z = in.nextDouble();
        Orig[i] = P[i] = Vec3(x,y,z);
    }
    for (int i=0;i<M;i++) {
        char c = in.nextChar(); (void)c;
        int a = in.nextInt() - 1, b = in.nextInt() - 1, cidx = in.nextInt() - 1;
        F[i].v[0]=a; F[i].v[1]=b; F[i].v[2]=cidx; F[i].alive=1;
    }
    simplifyMesh();
    if (outputSnap && bestSnap.valid) saveMesh(bestSnap.faces, bestSnap.alive, bestSnap.pos);
    else saveMesh(F, vAlive, P);
    return 0;
}
