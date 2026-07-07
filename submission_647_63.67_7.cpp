#include <bits/stdc++.h>
using namespace std;

// IMC 2026 simplifygeometry - aggressive perception-aware QEM/r-net edge contraction.
// C++17, no external dependencies.

struct Vec3 {
    double x, y, z;
    Vec3() : x(0), y(0), z(0) {}
    Vec3(double X,double Y,double Z):x(X),y(Y),z(Z){}
    Vec3 operator+(const Vec3& o) const { return Vec3(x+o.x,y+o.y,z+o.z); }
    Vec3 operator-(const Vec3& o) const { return Vec3(x-o.x,y-o.y,z-o.z); }
    Vec3 operator*(double s) const { return Vec3(x*s,y*s,z*s); }
    Vec3 operator/(double s) const { return Vec3(x/s,y/s,z/s); }
    Vec3& operator+=(const Vec3& o){ x+=o.x; y+=o.y; z+=o.z; return *this; }
};
static inline double dotv(const Vec3&a,const Vec3&b){ return a.x*b.x+a.y*b.y+a.z*b.z; }
static inline Vec3 crossv(const Vec3&a,const Vec3&b){ return Vec3(a.y*b.z-a.z*b.y, a.z*b.x-a.x*b.z, a.x*b.y-a.y*b.x); }
static inline double norm2(const Vec3&a){ return dotv(a,a); }
static inline double normv(const Vec3&a){ return sqrt(norm2(a)); }

struct Quadric {
    double q00,q01,q02,q03,q11,q12,q13,q22,q23,q33;
    Quadric(){ clear(); }
    void clear(){ q00=q01=q02=q03=q11=q12=q13=q22=q23=q33=0.0; }
    Quadric& operator+=(const Quadric& o){
        q00+=o.q00; q01+=o.q01; q02+=o.q02; q03+=o.q03;
        q11+=o.q11; q12+=o.q12; q13+=o.q13; q22+=o.q22; q23+=o.q23; q33+=o.q33;
        return *this;
    }
};
static inline Quadric operator+(Quadric a, const Quadric& b){ a += b; return a; }
static inline double qcost(const Quadric& q, const Vec3& p){
    const double x=p.x,y=p.y,z=p.z;
    return q.q00*x*x + 2*q.q01*x*y + 2*q.q02*x*z + 2*q.q03*x
         + q.q11*y*y + 2*q.q12*y*z + 2*q.q13*y
         + q.q22*z*z + 2*q.q23*z + q.q33;
}
static inline void addPlane(Quadric& q, const Vec3& n, double d, double w){
    const double a=n.x,b=n.y,c=n.z;
    q.q00 += w*a*a; q.q01 += w*a*b; q.q02 += w*a*c; q.q03 += w*a*d;
    q.q11 += w*b*b; q.q12 += w*b*c; q.q13 += w*b*d;
    q.q22 += w*c*c; q.q23 += w*c*d; q.q33 += w*d*d;
}

struct Face {
    int a,b,c;
    float nx,ny,nz;
    unsigned char alive;
};

struct FastInput {
    vector<char> buf; char* p;
    FastInput(){
        char tmp[1<<16]; size_t n;
        while((n=fread(tmp,1,sizeof(tmp),stdin))>0) buf.insert(buf.end(),tmp,tmp+n);
        buf.push_back('\0'); p=buf.data();
    }
    inline void skip(){ while(*p==' '||*p=='\n'||*p=='\r'||*p=='\t') ++p; }
    long nextLong(){ skip(); long s=1; if(*p=='-'){s=-1; ++p;} long v=0; while(*p>='0'&&*p<='9'){v=v*10+(*p-'0'); ++p;} return v*s; }
    double nextDouble(){ skip(); char* e; double v=strtod(p,&e); p=e; return v; }
    char nextChar(){ skip(); return *p++; }
};

struct OutBuf {
    static const int SZ = 1<<20;
    char buf[SZ]; int n=0;
    ~OutBuf(){ flush(); }
    inline void flush(){ if(n){ fwrite(buf,1,n,stdout); n=0; } }
    inline void put(const char* s, int len){
        if(len >= SZ){ flush(); fwrite(s,1,len,stdout); return; }
        if(n+len > SZ) flush();
        memcpy(buf+n,s,len); n += len;
    }
    template<class... Args> inline void printfmt(const char* fmt, Args... args){
        char tmp[160]; int len = snprintf(tmp,sizeof(tmp),fmt,args...); put(tmp,len);
    }
};

static int NV, NF;
static vector<Vec3> P;
static vector<Face> F;
static vector<vector<int>> Inc;
static vector<Quadric> Q;
static vector<Vec3> BLo, BHi, SumOrig;
static vector<int> Cnt, Ver;
static vector<unsigned char> AliveV;
static int liveV, liveF;
static double diagLen, hausTol, hausTol2;
static double geomComplexity = 0.0, highCurvFraction = 0.0;
static chrono::steady_clock::time_point tStart;
static double timeLimitSec = 19.2;

static inline bool faceContains(const Face& f, int v){ return f.a==v || f.b==v || f.c==v; }
static inline bool faceContains2(const Face& f, int u, int v){ return faceContains(f,u) && faceContains(f,v); }

static bool recomputeFaceNormal(int id){
    Face &f = F[id];
    Vec3 ab = P[f.b]-P[f.a], ac = P[f.c]-P[f.a];
    Vec3 cr = crossv(ab,ac);
    double l = normv(cr);
    if(l <= 1e-30) return false;
    cr = cr / l;
    f.nx = (float)cr.x; f.ny = (float)cr.y; f.nz = (float)cr.z;
    return true;
}

static void cleanInc(int v){
    if(!AliveV[v]) { Inc[v].clear(); return; }
    auto &a = Inc[v];
    int m=0;
    for(int id: a){
        if(id>=0 && id<NF && F[id].alive && faceContains(F[id],v)) a[m++] = id;
    }
    a.resize(m);
    if(a.size() > 1){ sort(a.begin(), a.end()); a.erase(unique(a.begin(), a.end()), a.end()); }
}

static vector<int> getNeighbors(int v){
    cleanInc(v);
    vector<int> nb; nb.reserve(16);
    for(int id: Inc[v]){
        const Face& f = F[id];
        if(f.a!=v && AliveV[f.a]) nb.push_back(f.a);
        if(f.b!=v && AliveV[f.b]) nb.push_back(f.b);
        if(f.c!=v && AliveV[f.c]) nb.push_back(f.c);
    }
    sort(nb.begin(), nb.end()); nb.erase(unique(nb.begin(), nb.end()), nb.end());
    return nb;
}

static void edgeFaces(int u, int v, vector<int>& out){
    out.clear();
    if(!AliveV[u] || !AliveV[v]) return;
    cleanInc(u); cleanInc(v);
    int scan = (Inc[u].size() <= Inc[v].size() ? u : v);
    for(int id: Inc[scan]) if(F[id].alive && faceContains2(F[id],u,v)) out.push_back(id);
}

static inline uint64_t edgeKey(int a, int b){
    if(a>b) swap(a,b);
    return (uint64_t)(uint32_t)a << 32 | (uint32_t)b;
}
static inline uint64_t faceKey3(int a, int b, int c){
    if(a>b) swap(a,b); if(b>c) swap(b,c); if(a>b) swap(a,b);
    return (uint64_t)a | ((uint64_t)b << 21) | ((uint64_t)c << 42);
}

static inline double maxBoxDist2(const Vec3& lo, const Vec3& hi, const Vec3& p){
    double md=0.0;
    for(int ix=0; ix<2; ++ix){ double x = ix?hi.x:lo.x;
        for(int iy=0; iy<2; ++iy){ double y = iy?hi.y:lo.y;
            for(int iz=0; iz<2; ++iz){ double z = iz?hi.z:lo.z;
                double dx=x-p.x, dy=y-p.y, dz=z-p.z;
                md = max(md, dx*dx+dy*dy+dz*dz);
            }
        }
    }
    return md;
}

static inline void projectAxis(const Vec3& p, int view, double& u, double& v){
    constexpr double D=2.5, f=800.0;
    double dep;
    switch(view){
        case 0: dep = D - p.x; u = p.y/dep*f; v = p.z/dep*f; break; // +X
        case 1: dep = D + p.x; u = p.y/dep*f; v = p.z/dep*f; break; // -X
        case 2: dep = D - p.y; u = p.x/dep*f; v = p.z/dep*f; break; // +Y
        case 3: dep = D + p.y; u = p.x/dep*f; v = p.z/dep*f; break; // -Y
        case 4: dep = D - p.z; u = p.x/dep*f; v = p.y/dep*f; break; // +Z
        default: dep = D + p.z; u = p.x/dep*f; v = p.y/dep*f; break; // -Z
    }
}
static inline double maxProjDist2(const Vec3& a, const Vec3& b){
    double best=0.0;
    for(int k=0;k<6;k++){
        double au,av,bu,bv; projectAxis(a,k,au,av); projectAxis(b,k,bu,bv);
        double du=au-bu, dv=av-bv; best=max(best,du*du+dv*dv);
    }
    return best;
}
static inline double normalCosLimit(double pixMove){
    // Tiny screen-space moves may alter a local face normal without visible impact.
    // Larger moves must preserve normals much more tightly.
    if(pixMove < 1.5) return 0.20;
    if(pixMove < 3.0) return 0.45;
    if(pixMove < 6.0) return 0.68;
    if(pixMove < 12.0) return 0.82;
    if(pixMove < 24.0) return 0.90;
    return 0.945;
}

static bool solveOptimal(const Quadric& q, Vec3& out){
    double a00=q.q00, a01=q.q01, a02=q.q02;
    double a11=q.q11, a12=q.q12, a22=q.q22;
    double b0=-q.q03, b1=-q.q13, b2=-q.q23;
    double det = a00*(a11*a22-a12*a12) - a01*(a01*a22-a12*a02) + a02*(a01*a12-a11*a02);
    double scale = fabs(a00)+fabs(a11)+fabs(a22)+fabs(a01)+fabs(a02)+fabs(a12)+1e-30;
    if(fabs(det) < 1e-12 * scale*scale*scale) return false;
    double c00 = (a11*a22-a12*a12);
    double c01 = -(a01*a22-a12*a02);
    double c02 = (a01*a12-a11*a02);
    double c10 = c01;
    double c11 = (a00*a22-a02*a02);
    double c12 = -(a00*a12-a01*a02);
    double c20 = c02;
    double c21 = c12;
    double c22 = (a00*a11-a01*a01);
    out.x = (c00*b0 + c01*b1 + c02*b2)/det;
    out.y = (c10*b0 + c11*b1 + c12*b2)/det;
    out.z = (c20*b0 + c21*b1 + c22*b2)/det;
    return isfinite(out.x) && isfinite(out.y) && isfinite(out.z);
}

struct Cand {
    int u=-1, v=-1, keep=-1, rem=-1;
    Vec3 p;
    double pri=0.0;
};

static bool choosePosition(int u, int v, Cand& cand){
    Quadric qq = Q[u] + Q[v];
    Vec3 lo(min(BLo[u].x,BLo[v].x), min(BLo[u].y,BLo[v].y), min(BLo[u].z,BLo[v].z));
    Vec3 hi(max(BHi[u].x,BHi[v].x), max(BHi[u].y,BHi[v].y), max(BHi[u].z,BHi[v].z));
    vector<Vec3> ps; ps.reserve(8);
    Vec3 opt;
    if(solveOptimal(qq,opt)) ps.push_back(opt);
    ps.push_back((P[u]+P[v])*0.5);
    ps.push_back((SumOrig[u]+SumOrig[v]) / double(Cnt[u]+Cnt[v]));
    ps.push_back(Vec3((lo.x+hi.x)*0.5, (lo.y+hi.y)*0.5, (lo.z+hi.z)*0.5));
    ps.push_back(P[u]); ps.push_back(P[v]);

    double tr = qq.q00 + qq.q11 + qq.q22 + 1e-30;
    double best = 1e300; Vec3 bestp;
    bool ok=false;
    double edgeLen2 = norm2(P[u]-P[v]);
    double edgePix2 = maxProjDist2(P[u],P[v]);
    for(const Vec3& p: ps){
        if(!isfinite(p.x) || !isfinite(p.y) || !isfinite(p.z)) continue;
        double bd2 = maxBoxDist2(lo,hi,p);
        if(bd2 > hausTol2) continue;
        double qc = max(0.0, qcost(qq,p)) / tr;
        // The bbox term breaks planar ties in favor of centered r-net representatives.
        double score = qc + 1e-7*edgeLen2 + 5e-10*edgePix2 + 1e-8*bd2;
        if(score < best){ best=score; bestp=p; ok=true; }
    }
    if(!ok) return false;
    cand.p=bestp; cand.pri=best;
    // Keep the larger dynamic cluster to reduce incident-list churn.
    if(Cnt[u] > Cnt[v] || (Cnt[u]==Cnt[v] && Inc[u].size() >= Inc[v].size())) { cand.keep=u; cand.rem=v; }
    else { cand.keep=v; cand.rem=u; }
    return true;
}

static bool fullCandidate(int u, int v, Cand& cand){
    if(u==v || !AliveV[u] || !AliveV[v]) return false;
    cand.u=u; cand.v=v;
    vector<int> efs; edgeFaces(u,v,efs);
    if(efs.size()!=2) return false;
    if(!choosePosition(u,v,cand)) return false;

    vector<int> nu = getNeighbors(u), nv = getNeighbors(v);
    int common=0;
    size_t i=0,j=0;
    while(i<nu.size() && j<nv.size()){
        if(nu[i]==nv[j]){ if(nu[i]!=u && nu[i]!=v) ++common; ++i; ++j; }
        else if(nu[i]<nv[j]) ++i; else ++j;
    }
    if(common != 2) return false;

    int keep=cand.keep, rem=cand.rem;
    vector<int> incident;
    incident.reserve(Inc[u].size()+Inc[v].size());
    incident.insert(incident.end(), Inc[u].begin(), Inc[u].end());
    incident.insert(incident.end(), Inc[v].begin(), Inc[v].end());
    sort(incident.begin(), incident.end()); incident.erase(unique(incident.begin(), incident.end()), incident.end());

    vector<uint64_t> keys; keys.reserve(incident.size());
    double pixMove = sqrt(max(maxProjDist2(P[u], cand.p), maxProjDist2(P[v], cand.p)));
    double cLim = normalCosLimit(pixMove);
    // Be slightly more permissive when the merge remains far inside the Hausdorff ball.
    Vec3 lo(min(BLo[u].x,BLo[v].x), min(BLo[u].y,BLo[v].y), min(BLo[u].z,BLo[v].z));
    Vec3 hi(max(BHi[u].x,BHi[v].x), max(BHi[u].y,BHi[v].y), max(BHi[u].z,BHi[v].z));
    double fill = sqrt(maxBoxDist2(lo,hi,cand.p)) / max(hausTol, 1e-30);
    if(fill < 0.55) cLim = min(cLim, 0.72);

    for(int id: incident){
        Face &f = F[id];
        if(!f.alive) continue;
        bool hu = faceContains(f,u), hv = faceContains(f,v);
        if(hu && hv) continue; // two faces removed by the contraction
        int a=f.a, b=f.b, c=f.c;
        if(a==u || a==v) a=keep;
        if(b==u || b==v) b=keep;
        if(c==u || c==v) c=keep;
        if(a==b || b==c || a==c) return false;
        keys.push_back(faceKey3(a,b,c));

        Vec3 pa = (f.a==u || f.a==v) ? cand.p : P[f.a];
        Vec3 pb = (f.b==u || f.b==v) ? cand.p : P[f.b];
        Vec3 pc = (f.c==u || f.c==v) ? cand.p : P[f.c];
        Vec3 cr = crossv(pb-pa, pc-pa);
        double l = normv(cr);
        if(l <= 1e-24 * max(1.0, diagLen*diagLen)) return false;
        Vec3 nn = cr / l;
        double d = nn.x*f.nx + nn.y*f.ny + nn.z*f.nz;
        if(d < cLim) return false;
    }
    if(keys.size() > 1){
        sort(keys.begin(), keys.end());
        for(size_t k=1;k<keys.size();++k) if(keys[k]==keys[k-1]) return false;
    }
    return true;
}

struct PQNode {
    float pri;
    int u,v,vu,vv;
    bool operator<(const PQNode& o) const { return pri > o.pri; }
};
static priority_queue<PQNode> pq;

static void pushEdge(int u, int v){
    if(u==v || !AliveV[u] || !AliveV[v]) return;
    Cand c; c.u=u; c.v=v;
    if(!choosePosition(u,v,c)) return;
    double p = c.pri;
    if(!isfinite(p)) return;
    if(p > 1e20) return;
    pq.push(PQNode{(float)p,u,v,Ver[u],Ver[v]});
}

static void applyCollapse(const Cand& cand){
    int keep=cand.keep, rem=cand.rem;
    int u=cand.u, v=cand.v;
    cleanInc(u); cleanInc(v);
    vector<int> incident;
    incident.reserve(Inc[u].size()+Inc[v].size());
    incident.insert(incident.end(), Inc[u].begin(), Inc[u].end());
    incident.insert(incident.end(), Inc[v].begin(), Inc[v].end());
    sort(incident.begin(), incident.end()); incident.erase(unique(incident.begin(), incident.end()), incident.end());

    P[keep] = cand.p;
    Q[keep] += Q[rem];
    BLo[keep].x = min(BLo[keep].x, BLo[rem].x); BLo[keep].y = min(BLo[keep].y, BLo[rem].y); BLo[keep].z = min(BLo[keep].z, BLo[rem].z);
    BHi[keep].x = max(BHi[keep].x, BHi[rem].x); BHi[keep].y = max(BHi[keep].y, BHi[rem].y); BHi[keep].z = max(BHi[keep].z, BHi[rem].z);
    SumOrig[keep] += SumOrig[rem]; Cnt[keep] += Cnt[rem];
    AliveV[rem]=0; ++Ver[keep]; ++Ver[rem]; --liveV;

    vector<int> newIncKeep; newIncKeep.reserve(incident.size());
    for(int id: incident){
        Face &f = F[id];
        if(!f.alive) continue;
        bool hKeep = faceContains(f, keep), hRem = faceContains(f, rem);
        if(hKeep && hRem){ f.alive=0; --liveF; continue; }
        if(hRem){
            if(f.a==rem) f.a=keep;
            if(f.b==rem) f.b=keep;
            if(f.c==rem) f.c=keep;
        }
        // The survivor moved, so every remaining incident face needs a fresh normal.
        if(faceContains(f, keep)){
            recomputeFaceNormal(id);
            newIncKeep.push_back(id);
        }
    }
    sort(newIncKeep.begin(), newIncKeep.end()); newIncKeep.erase(unique(newIncKeep.begin(), newIncKeep.end()), newIncKeep.end());
    Inc[keep].swap(newIncKeep);
    Inc[rem].clear();

    vector<int> nb = getNeighbors(keep);
    for(int x: nb) pushEdge(keep,x);
}

static double targetFraction(){
    if(NV <= 20) return 0.80;      // sample/sanity meshes; constraints decide exact count
    double base;
    if(NV <= 4000) base = 0.135;   // very low-res cases need more visual support
    else if(NV <= 6000) base = 0.120;
    else if(NV <= 30000) base = 0.083;
    else if(NV <= 70000) base = 0.074;
    else if(NV <= 500000) base = 0.060;
    else base = 0.055;

    // If the input has dense high-frequency normal variation, the normal-map SSIM,
    // not Hausdorff, is the limiting constraint. Keep more vertices only then.
    if(geomComplexity > 0.035 || highCurvFraction > 0.20) base = max(base, 0.50);
    else if(geomComplexity > 0.025 || highCurvFraction > 0.12) base = max(base, 0.38);
    else if(geomComplexity > 0.015 || highCurvFraction > 0.07) base = max(base, 0.26);
    else if(geomComplexity > 0.008 || highCurvFraction > 0.03) base = max(base, 0.16);
    return min(base, 0.60);
}

static void buildInitialQuadrics(){
    Q.assign(NV, Quadric());
    vector<Vec3> vN(NV);
    for(int i=0;i<NF;i++){
        recomputeFaceNormal(i);
        Face &f = F[i];
        Vec3 n(f.nx,f.ny,f.nz);
        Vec3 cr = crossv(P[f.b]-P[f.a], P[f.c]-P[f.a]);
        double area2 = max(1e-30, normv(cr));
        double area = 0.5 * area2;
        double d = -dotv(n, P[f.a]);
        // Area weighting follows QEM; the small floor keeps ultra-thin features visible.
        double w = area + 1e-12;
        addPlane(Q[f.a], n, d, w);
        addPlane(Q[f.b], n, d, w);
        addPlane(Q[f.c], n, d, w);
        Vec3 aw = n * area2;
        vN[f.a] += aw; vN[f.b] += aw; vN[f.c] += aw;
    }
    for(int i=0;i<NV;i++){
        double l = normv(vN[i]);
        if(l > 1e-30) vN[i] = vN[i] / l;
    }
    double acc = 0.0; long long cnt = 0, high = 0;
    for(int i=0;i<NF;i++){
        const Face &f = F[i];
        Vec3 n(f.nx,f.ny,f.nz);
        int vs[3] = {f.a,f.b,f.c};
        for(int k=0;k<3;k++){
            double val = max(0.0, 1.0 - dotv(n, vN[vs[k]]));
            acc += val; ++cnt;
            if(val > 0.02) ++high;
        }
    }
    if(cnt){ geomComplexity = acc / double(cnt); highCurvFraction = double(high) / double(cnt); }
}

static void simplify(){
    liveV = NV; liveF = NF;
    if(NV <= 1) return;
    buildInitialQuadrics();

    vector<uint64_t> edges; edges.reserve((size_t)NF*3);
    for(int i=0;i<NF;i++){
        const Face& f=F[i];
        edges.push_back(edgeKey(f.a,f.b)); edges.push_back(edgeKey(f.b,f.c)); edges.push_back(edgeKey(f.c,f.a));
    }
    sort(edges.begin(), edges.end()); edges.erase(unique(edges.begin(), edges.end()), edges.end());
    for(uint64_t k: edges){
        int u = (int)(k>>32), v = (int)(k & 0xffffffffu);
        pushEdge(u,v);
    }
    vector<uint64_t>().swap(edges);

    int target = max(4, (int)ceil(NV * targetFraction()));
    if(NV <= 20) target = max(1, NV-1);

    int iter=0;
    while(liveV > target && !pq.empty()){
        if((++iter & 4095)==0){
            double elapsed = chrono::duration<double>(chrono::steady_clock::now()-tStart).count();
            if(elapsed > timeLimitSec) break;
        }
        PQNode n = pq.top(); pq.pop();
        if(n.u==n.v || !AliveV[n.u] || !AliveV[n.v]) continue;
        if(Ver[n.u]!=n.vu || Ver[n.v]!=n.vv) continue;
        Cand c;
        if(!fullCandidate(n.u,n.v,c)) continue;
        applyCollapse(c);
    }
}

static void readInput(){
    FastInput in;
    NV = (int)in.nextLong(); NF = (int)in.nextLong();
    P.resize(NV); F.resize(NF); Inc.assign(NV, {});
    BLo.resize(NV); BHi.resize(NV); SumOrig.resize(NV); Cnt.assign(NV,1); Ver.assign(NV,0); AliveV.assign(NV,1);
    Vec3 mn(1e100,1e100,1e100), mx(-1e100,-1e100,-1e100);
    for(int i=0;i<NV;i++){
        char ch = in.nextChar(); (void)ch;
        double x=in.nextDouble(), y=in.nextDouble(), z=in.nextDouble();
        P[i]=Vec3(x,y,z); BLo[i]=BHi[i]=SumOrig[i]=P[i];
        mn.x=min(mn.x,x); mn.y=min(mn.y,y); mn.z=min(mn.z,z);
        mx.x=max(mx.x,x); mx.y=max(mx.y,y); mx.z=max(mx.z,z);
    }
    diagLen = normv(mx-mn); hausTol = 0.05 * diagLen * 0.995; hausTol2 = hausTol*hausTol;
    for(int i=0;i<NF;i++){
        char ch = in.nextChar(); (void)ch;
        int a=(int)in.nextLong()-1, b=(int)in.nextLong()-1, c=(int)in.nextLong()-1;
        F[i].a=a; F[i].b=b; F[i].c=c; F[i].alive=1; F[i].nx=F[i].ny=F[i].nz=0;
        Inc[a].push_back(i); Inc[b].push_back(i); Inc[c].push_back(i);
    }
}

static void writeOutput(){
    vector<int> mapv(NV, -1);
    int outV=0, outF=0;
    for(int i=0;i<NF;i++) if(F[i].alive) ++outF;
    for(int i=0;i<NV;i++) if(AliveV[i]) mapv[i]=outV++;
    // In the unlikely event that an isolated live vertex exists, it is harmless for validity only if V'>=1;
    // edge contractions on closed meshes should not create isolated live vertices before all faces disappear.
    OutBuf out;
    out.printfmt("%d %d\n", outV, outF);
    for(int i=0;i<NV;i++) if(mapv[i]>=0){
        out.printfmt("v %.10g %.10g %.10g\n", P[i].x, P[i].y, P[i].z);
    }
    for(int i=0;i<NF;i++) if(F[i].alive){
        out.printfmt("f %d %d %d\n", mapv[F[i].a]+1, mapv[F[i].b]+1, mapv[F[i].c]+1);
    }
    out.flush();
}

int main(){
    tStart = chrono::steady_clock::now();
    readInput();
    simplify();
    writeOutput();
    return 0;
}
