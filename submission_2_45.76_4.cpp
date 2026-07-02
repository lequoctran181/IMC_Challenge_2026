#include <bits/stdc++.h>
using namespace std;

struct Vec3 {
    double x, y, z;
    Vec3() : x(0), y(0), z(0) {}
    Vec3(double _x, double _y, double _z) : x(_x), y(_y), z(_z) {}
};
static inline Vec3 operator+(const Vec3& a, const Vec3& b){ return Vec3(a.x+b.x,a.y+b.y,a.z+b.z); }
static inline Vec3 operator-(const Vec3& a, const Vec3& b){ return Vec3(a.x-b.x,a.y-b.y,a.z-b.z); }
static inline Vec3 operator*(const Vec3& a, double s){ return Vec3(a.x*s,a.y*s,a.z*s); }
static inline Vec3 operator/(const Vec3& a, double s){ return Vec3(a.x/s,a.y/s,a.z/s); }
static inline double dotv(const Vec3& a, const Vec3& b){ return a.x*b.x+a.y*b.y+a.z*b.z; }
static inline Vec3 crossv(const Vec3& a, const Vec3& b){ return Vec3(a.y*b.z-a.z*b.y, a.z*b.x-a.x*b.z, a.x*b.y-a.y*b.x); }
static inline double norm2(const Vec3& a){ return dotv(a,a); }
static inline double normv(const Vec3& a){ return sqrt(norm2(a)); }
static inline double distv(const Vec3& a, const Vec3& b){ return normv(a-b); }

struct Quadric {
    double q[10];
    Quadric(){ memset(q, 0, sizeof(q)); }
    inline void addPlane(double a, double b, double c, double d, double w){
        q[0] += w*a*a; q[1] += w*a*b; q[2] += w*a*c; q[3] += w*a*d;
        q[4] += w*b*b; q[5] += w*b*c; q[6] += w*b*d;
        q[7] += w*c*c; q[8] += w*c*d;
        q[9] += w*d*d;
    }
    inline void add(const Quadric& o){ for(int i=0;i<10;i++) q[i]+=o.q[i]; }
};
static inline double qeval(const Quadric& Q, const Vec3& p){
    const double x=p.x,y=p.y,z=p.z;
    return Q.q[0]*x*x + 2.0*Q.q[1]*x*y + 2.0*Q.q[2]*x*z + 2.0*Q.q[3]*x
         + Q.q[4]*y*y + 2.0*Q.q[5]*y*z + 2.0*Q.q[6]*y
         + Q.q[7]*z*z + 2.0*Q.q[8]*z + Q.q[9];
}

struct Face {
    int v[3];
    unsigned char alive;
    Vec3 n;
    double area;
};

struct Vertex {
    Vec3 p;
    Quadric Q;
    double rad;
    int ver;
    unsigned char alive;
    vector<int> faces;
};

static vector<Vertex> Vtx;
static vector<Face> Fac;
static int aliveV = 0;
static int aliveF = 0;
static double hausEps = 0.0;
static double minArea2 = 1e-28;
static const double INF = 1e100;

struct FastScanner {
    vector<char> buf;
    char* p;
    FastScanner(){
        const size_t CH = 1<<20;
        char tmp[CH];
        size_t n;
        while((n=fread(tmp,1,CH,stdin))>0) buf.insert(buf.end(), tmp, tmp+n);
        buf.push_back('\0');
        p = buf.data();
    }
    inline void skip(){ while(*p==' ' || *p=='\n' || *p=='\r' || *p=='\t') ++p; }
    int nextInt(){ skip(); int sign=1; if(*p=='-'){sign=-1;++p;} int x=0; while(*p>='0'&&*p<='9'){ x=x*10+(*p-'0'); ++p; } return x*sign; }
    double nextDouble(){ skip(); char* e; double x=strtod(p,&e); p=e; return x; }
    char nextChar(){ skip(); return *p++; }
};

static inline bool hasVertex(const Face& f, int a){ return f.v[0]==a || f.v[1]==a || f.v[2]==a; }
static inline bool hasBoth(const Face& f, int a, int b){ return hasVertex(f,a) && hasVertex(f,b); }
static inline int otherInFace(const Face& f, int a, int b){
    for(int i=0;i<3;i++) if(f.v[i]!=a && f.v[i]!=b) return f.v[i];
    return -1;
}

static bool recomputeFace(int fid){
    Face& f = Fac[fid];
    const Vec3 &a = Vtx[f.v[0]].p, &b = Vtx[f.v[1]].p, &c = Vtx[f.v[2]].p;
    Vec3 cr = crossv(b-a, c-a);
    double l2 = norm2(cr);
    if(l2 <= minArea2) return false;
    double inv = 1.0 / sqrt(l2);
    f.n = cr * inv;
    f.area = 0.5 / inv;
    return true;
}

static Vec3 faceNormalWithReplacement(const Face& f, int keep, int kill, const Vec3& np, double &l2){
    Vec3 p[3];
    for(int i=0;i<3;i++){
        int id=f.v[i];
        p[i] = (id==keep || id==kill) ? np : Vtx[id].p;
    }
    Vec3 cr = crossv(p[1]-p[0], p[2]-p[0]);
    l2 = norm2(cr);
    if(l2 <= minArea2) return Vec3(0,0,0);
    return cr / sqrt(l2);
}

static void cleanFaces(int v){
    if(!Vtx[v].alive) return;
    vector<int>& a = Vtx[v].faces;
    size_t w=0;
    for(size_t i=0;i<a.size();++i){
        int fid = a[i];
        if(fid>=0 && fid<(int)Fac.size() && Fac[fid].alive && hasVertex(Fac[fid], v)) a[w++] = fid;
    }
    a.resize(w);
}

static vector<int> markArr;
static int markToken = 1;

static void collectNeighbors(int v, vector<int>& out){
    out.clear();
    cleanFaces(v);
    ++markToken;
    if(markToken == INT_MAX){ fill(markArr.begin(), markArr.end(), 0); markToken=1; }
    for(int fid: Vtx[v].faces){
        const Face& f=Fac[fid];
        for(int k=0;k<3;k++){
            int u=f.v[k];
            if(u!=v && Vtx[u].alive && markArr[u]!=markToken){
                markArr[u]=markToken;
                out.push_back(u);
            }
        }
    }
}

static bool solve3x3(const Quadric& Q, Vec3& ans){
    double a00=Q.q[0], a01=Q.q[1], a02=Q.q[2];
    double a11=Q.q[4], a12=Q.q[5];
    double a22=Q.q[7];
    double b0=-Q.q[3], b1=-Q.q[6], b2=-Q.q[8];
    double det = a00*(a11*a22-a12*a12) - a01*(a01*a22-a12*a02) + a02*(a01*a12-a11*a02);
    if(fabs(det) < 1e-13) return false;
    double invDet = 1.0/det;
    double dx = b0*(a11*a22-a12*a12) - a01*(b1*a22-a12*b2) + a02*(b1*a12-a11*b2);
    double dy = a00*(b1*a22-a12*b2) - b0*(a01*a22-a12*a02) + a02*(a01*b2-b1*a02);
    double dz = a00*(a11*b2-b1*a12) - a01*(a01*b2-b1*a02) + b0*(a01*a12-a11*a02);
    ans = Vec3(dx*invDet, dy*invDet, dz*invDet);
    if(!isfinite(ans.x) || !isfinite(ans.y) || !isfinite(ans.z)) return false;
    return true;
}

struct Candidate {
    Vec3 p;
    double cost;
    double newRad;
};

static bool bestGeometricCandidate(int a, int b, Candidate& best){
    if(!Vtx[a].alive || !Vtx[b].alive) return false;
    Quadric Q = Vtx[a].Q; Q.add(Vtx[b].Q);
    const Vec3 pa = Vtx[a].p, pb = Vtx[b].p;
    const double ra = Vtx[a].rad, rb = Vtx[b].rad;
    best.cost = INF;
    best.newRad = INF;
    bool ok=false;
    auto consider = [&](const Vec3& p){
        double nr = max(ra + distv(p, pa), rb + distv(p, pb));
        if(nr > hausEps) return;
        double c = qeval(Q, p);
        if(c < 0 && c > -1e-18) c = 0;
        double key = c + 1e-18 * (nr / max(hausEps, 1e-12));
        if(key < best.cost){ best.cost = key; best.p = p; best.newRad = nr; ok=true; }
    };
    Vec3 opt;
    if(solve3x3(Q, opt)) consider(opt);
    consider(pa);
    consider(pb);
    consider((pa+pb)*0.5);
    double d = distv(pa,pb);
    if(d > 1e-18){
        double t = (rb + d - ra) / (2.0*d);
        if(t < 0) t = 0; if(t > 1) t = 1;
        consider(pa*(1.0-t) + pb*t);
    }
    return ok;
}

struct HeapItem {
    double cost;
    int a,b,va,vb;
    bool operator<(const HeapItem& o) const { return cost > o.cost; } // min heap
};
static priority_queue<HeapItem> heapq;

static void pushEdge(int a, int b){
    if(a==b || !Vtx[a].alive || !Vtx[b].alive) return;
    if(a>b) swap(a,b);
    Candidate c;
    if(!bestGeometricCandidate(a,b,c)) return;
    heapq.push({c.cost, a, b, Vtx[a].ver, Vtx[b].ver});
}

static bool linkCondition(int a, int b){
    if(!Vtx[a].alive || !Vtx[b].alive || a==b) return false;
    cleanFaces(a); cleanFaces(b);

    int edgeFaces = 0;
    if(Vtx[a].faces.size() < Vtx[b].faces.size()){
        for(int fid: Vtx[a].faces) if(Fac[fid].alive && hasBoth(Fac[fid],a,b)) ++edgeFaces;
    } else {
        for(int fid: Vtx[b].faces) if(Fac[fid].alive && hasBoth(Fac[fid],a,b)) ++edgeFaces;
    }
    if(edgeFaces != 2) return false;

    vector<int> nbA, nbB;
    collectNeighbors(a, nbA);
    collectNeighbors(b, nbB);
    ++markToken;
    if(markToken == INT_MAX){ fill(markArr.begin(), markArr.end(), 0); markToken=1; }
    for(int u: nbA) markArr[u] = markToken;
    int common = 0;
    for(int u: nbB) if(markArr[u] == markToken) ++common;
    return common == 2;
}

static bool geometryOK(int keep, int kill, const Vec3& np){
    const double minDot = -1e-6;
    cleanFaces(keep); cleanFaces(kill);

    for(int pass=0; pass<2; ++pass){
        int v = pass ? kill : keep;
        for(int fid: Vtx[v].faces){
            Face& f = Fac[fid];
            if(!f.alive) continue;
            bool hk = hasVertex(f, keep), hl = hasVertex(f, kill);
            if(hk && hl) continue; 
            double l2;
            Vec3 nn = faceNormalWithReplacement(f, keep, kill, np, l2);
            if(l2 <= minArea2) return false;
            if(dotv(nn, f.n) < minDot) return false;
        }
    }
    return true;
}

static bool collapseEdge(int a, int b, const Candidate& cand){
    if(!linkCondition(a,b)) return false;

    cleanFaces(a); cleanFaces(b);
    int keep=a, kill=b;
    if(Vtx[a].faces.size() < Vtx[b].faces.size()) { keep=b; kill=a; }

    if(!geometryOK(keep, kill, cand.p)) return false;

    vector<int> killFaces = Vtx[kill].faces; 

    for(int fid: killFaces){
        Face& f = Fac[fid];
        if(!f.alive) continue;
        bool hk = hasVertex(f, keep), hl = hasVertex(f, kill);
        if(!hl) continue;
        if(hk){
            f.alive = 0;
            --aliveF;
        } else {
            for(int i=0;i<3;i++) if(f.v[i] == kill) f.v[i] = keep;
            Vtx[keep].faces.push_back(fid);
        }
    }

    Vtx[keep].p = cand.p;
    Vtx[keep].rad = cand.newRad;
    Vtx[keep].Q.add(Vtx[kill].Q);
    Vtx[keep].ver++;
    Vtx[kill].alive = 0;
    Vtx[kill].ver++;
    Vtx[kill].faces.clear();
    --aliveV;

    cleanFaces(keep);
    for(int fid: Vtx[keep].faces){
        if(Fac[fid].alive) recomputeFace(fid);
    }
    return true;
}

static void pushEdgesAround(int v){
    static vector<int> nb;
    collectNeighbors(v, nb);
    for(int u: nb) pushEdge(v,u);
}

static uint64_t edgeKey(int a, int b){
    if(a>b) swap(a,b);
    return (uint64_t)(uint32_t)a << 32 | (uint32_t)b;
}

static void loadMesh(){
    FastScanner fs;
    int n = fs.nextInt();
    int m = fs.nextInt();
    Vtx.assign(n, Vertex());
    Fac.assign(m, Face());
    aliveV = n; aliveF = m;

    double xmin=1e100, ymin=1e100, zmin=1e100, xmax=-1e100, ymax=-1e100, zmax=-1e100;
    for(int i=0;i<n;i++){
        char ch = fs.nextChar(); (void)ch;
        double x=fs.nextDouble(), y=fs.nextDouble(), z=fs.nextDouble();
        Vtx[i].p = Vec3(x,y,z);
        Vtx[i].rad = 0.0;
        Vtx[i].ver = 0;
        Vtx[i].alive = 1;
        xmin=min(xmin,x); ymin=min(ymin,y); zmin=min(zmin,z);
        xmax=max(xmax,x); ymax=max(ymax,y); zmax=max(zmax,z);
    }
    double dx=xmax-xmin, dy=ymax-ymin, dz=zmax-zmin;
    double diag = sqrt(dx*dx+dy*dy+dz*dz);
    hausEps = 0.05 * diag * 0.999999; 
    minArea2 = max(1e-30, diag*diag*diag*diag * 1e-28);

    for(int i=0;i<m;i++){
        char ch = fs.nextChar(); (void)ch;
        int a=fs.nextInt()-1, b=fs.nextInt()-1, c=fs.nextInt()-1;
        Fac[i].v[0]=a; Fac[i].v[1]=b; Fac[i].v[2]=c; Fac[i].alive=1;
        Vtx[a].faces.push_back(i); Vtx[b].faces.push_back(i); Vtx[c].faces.push_back(i);
    }

    for(int i=0;i<m;i++){
        recomputeFace(i);
        Face& f=Fac[i];
        const Vec3& p = Vtx[f.v[0]].p;
        double d = -dotv(f.n, p);
        double w = max(f.area, 1e-14);
        for(int k=0;k<3;k++) Vtx[f.v[k]].Q.addPlane(f.n.x, f.n.y, f.n.z, d, w);
    }
}

static double chooseTargetRatio(int n, int m){
    if(n <= 12) return 0.01;
    if(n <= 1000) return 0.18;
    if(n <= 6000) return 0.103;
    if(n <= 30000) return 0.088;
    if(n <= 70000) return 0.082;
    if(n <= 500000) return 0.078;
    return 0.073;
}

static void simplify(){
    const int n0 = (int)Vtx.size();
    const int m0 = (int)Fac.size();
    if(n0 <= 4) return;

    markArr.assign(n0, 0);

    vector<uint64_t> edges;
    edges.reserve((size_t)m0*3);
    for(int i=0;i<m0;i++){
        edges.push_back(edgeKey(Fac[i].v[0], Fac[i].v[1]));
        edges.push_back(edgeKey(Fac[i].v[1], Fac[i].v[2]));
        edges.push_back(edgeKey(Fac[i].v[2], Fac[i].v[0]));
    }
    sort(edges.begin(), edges.end());
    edges.erase(unique(edges.begin(), edges.end()), edges.end());
    for(uint64_t e: edges){
        int a = (int)(e >> 32);
        int b = (int)(e & 0xffffffffu);
        pushEdge(a,b);
    }
    vector<uint64_t>().swap(edges);

    double ratio = chooseTargetRatio(n0, m0);
    int target = max(4, (int)ceil(n0 * ratio));

    long long popCount = 0, failCount = 0;
    const long long popLimit = max<long long>(2000000LL, (long long)n0 * 90LL);

    while(aliveV > target && !heapq.empty() && popCount < popLimit){
        HeapItem it = heapq.top(); heapq.pop();
        ++popCount;
        int a=it.a, b=it.b;
        if(a<0 || b<0 || a>=n0 || b>=n0) continue;
        if(!Vtx[a].alive || !Vtx[b].alive) continue;
        if(Vtx[a].ver != it.va || Vtx[b].ver != it.vb) continue;

        Candidate cand;
        if(!bestGeometricCandidate(a,b,cand)) { ++failCount; continue; }
        if(cand.cost > it.cost * 1.000001 + 1e-18){
            heapq.push({cand.cost, min(a,b), max(a,b), Vtx[min(a,b)].ver, Vtx[max(a,b)].ver});
            continue;
        }
        if(collapseEdge(a,b,cand)){
            int keep = Vtx[a].alive ? a : b;
            if(Vtx[a].alive && Vtx[b].alive) keep = a;
            if(!Vtx[keep].alive) keep = (keep==a?b:a);
            pushEdgesAround(keep);
        } else {
            ++failCount;
        }
    }
}

static void saveMesh(){
    int n = (int)Vtx.size();
    vector<int> id(n, -1);
    vector<unsigned char> used(n, 0);
    int newN=0, newF=0;
    for(const Face& f: Fac) if(f.alive){
        ++newF;
        used[f.v[0]] = used[f.v[1]] = used[f.v[2]] = 1;
    }
    for(int i=0;i<n;i++) if(used[i]) id[i] = ++newN;

    string out;
    out.reserve((size_t)newN*44 + (size_t)newF*28 + 64);
    char line[128];
    out.append(line, snprintf(line, sizeof(line), "%d %d\n", newN, newF));
    for(int i=0;i<n;i++) if(used[i]){
        const Vec3& p = Vtx[i].p;
        out.append(line, snprintf(line, sizeof(line), "v %.10g %.10g %.10g\n", p.x, p.y, p.z));
    }
    for(const Face& f: Fac) if(f.alive){
        out.append(line, snprintf(line, sizeof(line), "f %d %d %d\n", id[f.v[0]], id[f.v[1]], id[f.v[2]]));
    }
    fwrite(out.data(), 1, out.size(), stdout);
}

int main(){
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    loadMesh();
    simplify();
    saveMesh();
    return 0;
}
