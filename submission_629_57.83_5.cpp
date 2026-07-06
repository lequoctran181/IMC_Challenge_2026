#include <bits/stdc++.h>
using namespace std;

struct Vec3 {
    double x, y, z;
    Vec3() : x(0), y(0), z(0) {}
    Vec3(double X, double Y, double Z) : x(X), y(Y), z(Z) {}
    Vec3 operator+(const Vec3& o) const { return Vec3(x+o.x, y+o.y, z+o.z); }
    Vec3 operator-(const Vec3& o) const { return Vec3(x-o.x, y-o.y, z-o.z); }
    Vec3 operator*(double s) const { return Vec3(x*s, y*s, z*s); }
    Vec3 operator/(double s) const { return Vec3(x/s, y/s, z/s); }
    Vec3& operator+=(const Vec3& o){ x+=o.x; y+=o.y; z+=o.z; return *this; }
};
static inline double dotv(const Vec3& a, const Vec3& b){ return a.x*b.x + a.y*b.y + a.z*b.z; }
static inline Vec3 crossv(const Vec3& a, const Vec3& b){ return Vec3(a.y*b.z-a.z*b.y, a.z*b.x-a.x*b.z, a.x*b.y-a.y*b.x); }
static inline double norm2(const Vec3& a){ return dotv(a,a); }
static inline double normv(const Vec3& a){ return sqrt(norm2(a)); }
static inline bool finitev(const Vec3& a){ return isfinite(a.x) && isfinite(a.y) && isfinite(a.z); }

struct Face { int v[3]; unsigned char active; };

struct Quadric {
    double a00,a01,a02,a03,a11,a12,a13,a22,a23,a33;
    Quadric(){ clear(); }
    void clear(){ a00=a01=a02=a03=a11=a12=a13=a22=a23=a33=0.0; }
    Quadric& operator+=(const Quadric& q){
        a00+=q.a00; a01+=q.a01; a02+=q.a02; a03+=q.a03;
        a11+=q.a11; a12+=q.a12; a13+=q.a13; a22+=q.a22; a23+=q.a23; a33+=q.a33; return *this;
    }
};
static inline Quadric qadd(const Quadric& a, const Quadric& b){ Quadric r=a; r+=b; return r; }
static inline void addPlane(Quadric& q, double a, double b, double c, double d, double w){
    q.a00 += w*a*a; q.a01 += w*a*b; q.a02 += w*a*c; q.a03 += w*a*d;
    q.a11 += w*b*b; q.a12 += w*b*c; q.a13 += w*b*d;
    q.a22 += w*c*c; q.a23 += w*c*d; q.a33 += w*d*d;
}
static inline double qeval(const Quadric& q, const Vec3& p){
    double x=p.x,y=p.y,z=p.z;
    double r = 0.0;
    r += q.a00*x*x + 2.0*q.a01*x*y + 2.0*q.a02*x*z + 2.0*q.a03*x;
    r += q.a11*y*y + 2.0*q.a12*y*z + 2.0*q.a13*y;
    r += q.a22*z*z + 2.0*q.a23*z + q.a33;
    return r < 0.0 && r > -1e-12 ? 0.0 : r;
}

struct BBox {
    double mn[3], mx[3];
};
static inline BBox mergeBox(const BBox& a, const BBox& b){
    BBox r;
    for(int i=0;i<3;i++){ r.mn[i]=min(a.mn[i],b.mn[i]); r.mx[i]=max(a.mx[i],b.mx[i]); }
    return r;
}

static int N, M;
static vector<Vec3> P;
static vector<Face> F;
static vector<Quadric> Q;
static vector<BBox> Box;
static vector<vector<int>> Inc;
static vector<int> Ver;
static vector<unsigned char> VActive;
static vector<int> faceMark, vMark, cMark;
static int faceStamp = 1, vStamp = 1, cStamp = 1;
static double diagLen = 1.0, hausEps = 0.05, hausEps2 = 0.0025;
static int activeV = 0, activeF = 0;
static chrono::steady_clock::time_point startTime;

static inline uint64_t edgeKey(int a, int b){ if(a>b) swap(a,b); return (uint64_t)(uint32_t)a << 32 | (uint32_t)b; }
static inline uint64_t faceKey3(int a, int b, int c){
    if(a>b) swap(a,b); if(b>c) swap(b,c); if(a>b) swap(a,b);
    return ((uint64_t)(uint32_t)a << 42) ^ ((uint64_t)(uint32_t)b << 21) ^ (uint64_t)(uint32_t)c;
}
static inline bool faceContains(const Face& f, int x){ return f.v[0]==x || f.v[1]==x || f.v[2]==x; }
static inline int thirdVertex(const Face& f, int a, int b){ for(int i=0;i<3;i++) if(f.v[i]!=a && f.v[i]!=b) return f.v[i]; return -1; }

struct EdgeRef { uint64_t key; int fid; bool operator<(const EdgeRef& o) const { return key < o.key || (key==o.key && fid < o.fid); } };
struct HeapNode {
    double cost;
    int u, v, vu, vv;
    bool operator<(const HeapNode& o) const { return cost > o.cost; }
};
static priority_queue<HeapNode> heapq;

static vector<char> inputBuf;
static void slurpInput(){
    const size_t CH = 1<<20;
    inputBuf.clear(); inputBuf.reserve(1<<26);
    char tmp[CH]; size_t n;
    while((n=fread(tmp,1,CH,stdin))>0) inputBuf.insert(inputBuf.end(), tmp, tmp+n);
    inputBuf.push_back('\0');
}
static inline void skipws(char*& p){ while(*p==' '||*p=='\n'||*p=='\r'||*p=='\t') ++p; }
static void loadMesh(vector<EdgeRef>& edgeRefs, vector<int>& deg){
    slurpInput(); char* p = inputBuf.data();
    N = (int)strtol(p,&p,10); M = (int)strtol(p,&p,10);
    P.resize(N); F.resize(M); Q.assign(N, Quadric()); Box.resize(N); Ver.assign(N,0); VActive.assign(N,1); deg.assign(N,0);
    double mnx=1e100,mny=1e100,mnz=1e100,mxx=-1e100,mxy=-1e100,mxz=-1e100;
    for(int i=0;i<N;i++){
        skipws(p); if(*p=='v') ++p;
        double x=strtod(p,&p), y=strtod(p,&p), z=strtod(p,&p);
        P[i]=Vec3(x,y,z);
        Box[i].mn[0]=Box[i].mx[0]=x; Box[i].mn[1]=Box[i].mx[1]=y; Box[i].mn[2]=Box[i].mx[2]=z;
        mnx=min(mnx,x); mny=min(mny,y); mnz=min(mnz,z); mxx=max(mxx,x); mxy=max(mxy,y); mxz=max(mxz,z);
    }
    edgeRefs.reserve((size_t)M*3);
    for(int i=0;i<M;i++){
        skipws(p); if(*p=='f') ++p;
        int a=(int)strtol(p,&p,10)-1, b=(int)strtol(p,&p,10)-1, c=(int)strtol(p,&p,10)-1;
        F[i].v[0]=a; F[i].v[1]=b; F[i].v[2]=c; F[i].active=1;
        deg[a]++; deg[b]++; deg[c]++;
        edgeRefs.push_back({edgeKey(a,b), i});
        edgeRefs.push_back({edgeKey(b,c), i});
        edgeRefs.push_back({edgeKey(c,a), i});
    }
    double dx=mxx-mnx, dy=mxy-mny, dz=mxz-mnz;
    diagLen=sqrt(dx*dx+dy*dy+dz*dz);
    if(diagLen <= 0 || !isfinite(diagLen)) diagLen = 1.0;
    hausEps = 0.05 * diagLen * 0.997;
    hausEps2 = hausEps * hausEps;
    activeV=N; activeF=M;
}

static bool solve3(const Quadric& q, Vec3& out){
    double a=q.a00,b=q.a01,c=q.a02,d=q.a11,e=q.a12,f=q.a22;
    double r0=-q.a03, r1=-q.a13, r2=-q.a23;
    double det = a*(d*f-e*e) - b*(b*f-e*c) + c*(b*e-d*c);
    double scale = fabs(a)+fabs(d)+fabs(f)+2.0*(fabs(b)+fabs(c)+fabs(e));
    if(fabs(det) < 1e-13 * max(1.0, scale*scale*scale)) return false;
    double dx = r0*(d*f-e*e) - b*(r1*f-e*r2) + c*(r1*e-d*r2);
    double dy = a*(r1*f-e*r2) - r0*(b*f-e*c) + c*(b*r2-r1*c);
    double dz = a*(d*r2-r1*e) - b*(b*r2-r1*c) + r0*(b*e-d*c);
    out = Vec3(dx/det, dy/det, dz/det);
    return finitev(out);
}

static inline bool bboxOK(const BBox& b, const Vec3& p){
    double xs[2]={b.mn[0], b.mx[0]}, ys[2]={b.mn[1], b.mx[1]}, zs[2]={b.mn[2], b.mx[2]};
    for(int ix=0;ix<2;ix++) for(int iy=0;iy<2;iy++) for(int iz=0;iz<2;iz++){
        double dx=xs[ix]-p.x, dy=ys[iy]-p.y, dz=zs[iz]-p.z;
        if(dx*dx+dy*dy+dz*dz > hausEps2) return false;
    }
    return true;
}

struct PosCost { Vec3 p; double qerr; double cost; };
static vector<PosCost> genPositions(int u, int v){
    Quadric q = qadd(Q[u], Q[v]);
    BBox b = mergeBox(Box[u], Box[v]);
    Vec3 a=P[u], c=P[v], mid=(a+c)*0.5;
    Vec3 opt; vector<Vec3> cand; cand.reserve(7);
    cand.push_back(a); cand.push_back(c); cand.push_back(mid);
    if(solve3(q,opt)) cand.push_back(opt);
    cand.push_back(a*0.75 + c*0.25);
    cand.push_back(a*0.25 + c*0.75);
    vector<PosCost> res; res.reserve(cand.size());
    for(const Vec3& p: cand){
        if(!finitev(p) || !bboxOK(b,p)) continue;
        double qe = qeval(q,p); if(!isfinite(qe)) continue;
        if(qe < 0 && qe > -1e-9) qe = 0;
        double len2 = norm2(a-c);
        double cst = max(0.0, qe) + len2 * 1e-16;
        bool dup=false;
        for(auto &r: res) if(norm2(r.p-p) < 1e-28){ if(cst < r.cost){ r.p=p; r.qerr=qe; r.cost=cst; } dup=true; break; }
        if(!dup) res.push_back({p, qe, cst});
    }
    sort(res.begin(), res.end(), [](const PosCost& A, const PosCost& B){ return A.cost < B.cost; });
    return res;
}
static bool cheapCost(int u, int v, double& cost){
    if(u==v || !VActive[u] || !VActive[v]) return false;
    vector<PosCost> r = genPositions(u,v);
    if(r.empty()) return false;
    cost = r[0].cost;
    return isfinite(cost);
}
static void pushEdge(int u, int v){
    if(u==v || !VActive[u] || !VActive[v]) return;
    double c;
    if(cheapCost(u,v,c)) heapq.push({c,u,v,Ver[u],Ver[v]});
}

static inline Vec3 faceNormalCurrent(const Face& f){
    return crossv(P[f.v[1]]-P[f.v[0]], P[f.v[2]]-P[f.v[0]]);
}
static inline Vec3 newPosFor(int idx, int keep, int kill, const Vec3& p){ return (idx==keep || idx==kill) ? p : P[idx]; }

static bool gatherAndValidate(int keep, int kill, const Vec3& np, vector<int>& star, double& qerrOut){
    if(!VActive[keep] || !VActive[kill] || keep==kill) return false;
    star.clear();
    if(++faceStamp == INT_MAX){ fill(faceMark.begin(), faceMark.end(), 0); faceStamp=1; }
    auto addList = [&](int v){
        for(int fid: Inc[v]){
            if(fid<0 || fid>=M) continue;
            if(!F[fid].active) continue;
            if(faceMark[fid]==faceStamp) continue;
            if(!faceContains(F[fid], v)) continue;
            faceMark[fid]=faceStamp; star.push_back(fid);
        }
    };
    addList(keep); addList(kill);
    int shared=0, opp[2]={-1,-1};
    for(int fid: star){
        Face &f=F[fid]; bool hk=faceContains(f,keep), hl=faceContains(f,kill);
        if(hk && hl){ if(shared<2) opp[shared]=thirdVertex(f,keep,kill); shared++; }
    }
    if(shared != 2 || opp[0] < 0 || opp[1] < 0 || opp[0] == opp[1]) return false;

    if(++vStamp == INT_MAX){ fill(vMark.begin(), vMark.end(), 0); vStamp=1; }
    if(++cStamp == INT_MAX){ fill(cMark.begin(), cMark.end(), 0); cStamp=1; }
    for(int fid: star){
        Face &f=F[fid]; if(!faceContains(f,keep)) continue;
        for(int i=0;i<3;i++){ int w=f.v[i]; if(w!=keep && w!=kill) vMark[w]=vStamp; }
    }
    int commonCnt=0; int commons[4]={-1,-1,-1,-1};
    for(int fid: star){
        Face &f=F[fid]; if(!faceContains(f,kill)) continue;
        for(int i=0;i<3;i++){
            int w=f.v[i]; if(w==keep || w==kill) continue;
            if(vMark[w]==vStamp && cMark[w]!=cStamp){ cMark[w]=cStamp; if(commonCnt<4) commons[commonCnt]=w; commonCnt++; }
        }
    }
    if(commonCnt != 2) return false;
    bool okOpp = ((commons[0]==opp[0] && commons[1]==opp[1]) || (commons[0]==opp[1] && commons[1]==opp[0]));
    if(!okOpp) return false;

    vector<uint64_t> keys; keys.reserve(star.size());
    const double areaTiny = max(1e-30, diagLen*diagLen*1e-28);
    for(int fid: star){
        Face &f=F[fid]; bool hk=false, hl=false;
        int nv[3];
        for(int i=0;i<3;i++){ hk = hk || (f.v[i]==keep); hl = hl || (f.v[i]==kill); nv[i] = (f.v[i]==kill ? keep : f.v[i]); }
        if(hk && hl) continue;
        if(nv[0]==nv[1] || nv[1]==nv[2] || nv[2]==nv[0]) return false;
        keys.push_back(faceKey3(nv[0],nv[1],nv[2]));
        Vec3 oldn = faceNormalCurrent(f);
        Vec3 a = newPosFor(f.v[0],keep,kill,np), b = newPosFor(f.v[1],keep,kill,np), c = newPosFor(f.v[2],keep,kill,np);
        Vec3 newn = crossv(b-a,c-a);
        double oldl = normv(oldn), newl = normv(newn);
        if(!(newl > areaTiny) || !isfinite(newl)) return false;
        if(oldl > areaTiny){
            double d = dotv(oldn,newn) / (oldl*newl);
            if(d < -1e-4) return false;
        }
    }
    sort(keys.begin(), keys.end());
    for(size_t i=1;i<keys.size();i++) if(keys[i]==keys[i-1]) return false;
    Quadric q=qadd(Q[keep],Q[kill]); qerrOut = qeval(q,np); if(qerrOut < 0 && qerrOut > -1e-9) qerrOut=0;
    return true;
}

static bool chooseCollapse(int& keep, int& kill, Vec3& np, vector<int>& star, double& qerr){
    if(Inc[keep].size() < Inc[kill].size()) swap(keep, kill);
    vector<PosCost> cand = genPositions(keep, kill);
    if(cand.empty()) return false;
    for(const auto& pc: cand){
        if(gatherAndValidate(keep, kill, pc.p, star, qerr)){
            np = pc.p;
            return true;
        }
    }
    return false;
}

static void doCollapse(int keep, int kill, const Vec3& np, const vector<int>& star){
    Quadric nq = qadd(Q[keep], Q[kill]);
    BBox nb = mergeBox(Box[keep], Box[kill]);
    for(int fid: star){
        Face &f=F[fid]; if(!f.active) continue;
        bool hk=faceContains(f,keep), hl=faceContains(f,kill);
        if(hk && hl){ f.active=0; activeF--; }
        else if(hl){
            for(int i=0;i<3;i++) if(f.v[i]==kill) f.v[i]=keep;
            Inc[keep].push_back(fid);
        }
    }
    P[keep]=np; Q[keep]=nq; Box[keep]=nb;
    VActive[kill]=0; activeV--; Ver[keep]++; Ver[kill]++;

    if(++faceStamp == INT_MAX){ fill(faceMark.begin(), faceMark.end(), 0); faceStamp=1; }
    for(int fid: Inc[keep]){
        if(!F[fid].active) continue;
        if(faceMark[fid]==faceStamp) continue;
        if(!faceContains(F[fid], keep)) continue;
        faceMark[fid]=faceStamp;
        for(int i=0;i<3;i++){ int w=F[fid].v[i]; if(w!=keep && VActive[w]) pushEdge(keep,w); }
    }
}

static void buildIncidence(const vector<int>& deg){
    Inc.clear(); Inc.resize(N);
    for(int i=0;i<N;i++) Inc[i].reserve((size_t)deg[i] + 8);
    for(int i=0;i<M;i++) for(int j=0;j<3;j++) Inc[F[i].v[j]].push_back(i);
    faceMark.assign(M,0); vMark.assign(N,0); cMark.assign(N,0);
}

static void buildQuadricsAndEdges(vector<EdgeRef>& edgeRefs, double& sharpFrac, double& turnAvg){
    vector<Vec3> fn(M);
    vector<double> fa(M);
    double totalArea=0.0;
    for(int i=0;i<M;i++){
        Vec3 n = crossv(P[F[i].v[1]]-P[F[i].v[0]], P[F[i].v[2]]-P[F[i].v[0]]);
        double l = normv(n);
        if(l > 0){ fn[i]=n/l; fa[i]=0.5*l; totalArea += fa[i]; }
        else { fn[i]=Vec3(0,0,1); fa[i]=0; }
    }
    double avgArea = (M>0 && totalArea>0) ? totalArea / (double)M : 1.0;
    for(int i=0;i<M;i++){
        Face& f=F[i];
        Vec3 n = fn[i];
        double d = -dotv(n, P[f.v[0]]);
        double w = fa[i] / avgArea;
        if(!isfinite(w) || w<=0) w = 1.0;
        if(w < 0.02) w = 0.02; else if(w > 50.0) w = 50.0;
        addPlane(Q[f.v[0]], n.x,n.y,n.z,d,w);
        addPlane(Q[f.v[1]], n.x,n.y,n.z,d,w);
        addPlane(Q[f.v[2]], n.x,n.y,n.z,d,w);
    }
    sort(edgeRefs.begin(), edgeRefs.end());
    long long edgeCnt=0, sharpCnt=0; double turnSum=0.0;
    for(size_t i=0;i<edgeRefs.size();){
        size_t j=i+1; while(j<edgeRefs.size() && edgeRefs[j].key==edgeRefs[i].key) j++;
        uint64_t key=edgeRefs[i].key; int u=(int)(key>>32), v=(int)(key & 0xffffffffu);
        pushEdge(u,v);
        if(j-i>=2){
            int f0=edgeRefs[i].fid, f1=edgeRefs[i+1].fid;
            double d = dotv(fn[f0], fn[f1]); if(d<-1)d=-1; if(d>1)d=1;
            double t = sqrt(max(0.0, 1.0-d));
            turnSum += t;
            if(d < 0.965) sharpCnt++;
        }
        edgeCnt++;
        i=j;
    }
    sharpFrac = edgeCnt ? (double)sharpCnt / (double)edgeCnt : 0.0;
    turnAvg = edgeCnt ? turnSum / (double)edgeCnt : 0.0;
    vector<EdgeRef>().swap(edgeRefs);
}

static int chooseTarget(double sharpFrac, double turnAvg){
    if(N <= 4) return N;
    if(N <= 20) return max(4, N-1);
    double ratio;
    if(N < 8000) ratio = 0.080;
    else if(N < 30000) ratio = 0.050;
    else if(N < 80000) ratio = 0.047;
    else if(N < 500000) ratio = 0.044;
    else ratio = 0.041;
    if(sharpFrac > 0.38 || turnAvg > 0.18) ratio = max(ratio, 0.225);
    else if(sharpFrac > 0.28 || turnAvg > 0.14) ratio = max(ratio, 0.180);
    else if(sharpFrac > 0.16 || turnAvg > 0.095) ratio = max(ratio, 0.125);
    else if(sharpFrac > 0.07 || turnAvg > 0.055) ratio = max(ratio, 0.090);
    if(sharpFrac < 0.012 && turnAvg < 0.030) ratio -= 0.006;
    if(ratio < 0.035) ratio = 0.035;
    if(ratio > 0.240) ratio = 0.240;
    int t = (int)ceil((double)N * ratio);
    if(t < 4) t = 4;
    if(t > N) t = N;
    return t;
}

static void simplify(){
    startTime = chrono::steady_clock::now();
    vector<EdgeRef> edgeRefs; vector<int> deg;
    loadMesh(edgeRefs, deg);
    vector<char>().swap(inputBuf);
    buildIncidence(deg);
    vector<int>().swap(deg);
    double sharpFrac=0.0, turnAvg=0.0;
    buildQuadricsAndEdges(edgeRefs, sharpFrac, turnAvg);
    int target = chooseTarget(sharpFrac, turnAvg);
    const double planarLimit = 1e-11 * max(1.0, (double)M / max(1, N));
    vector<int> star; star.reserve(256);
    long long attempts=0, collapses=0;
    while(!heapq.empty()){
        if((attempts++ & 4095LL) == 0){
            double elapsed = chrono::duration<double>(chrono::steady_clock::now()-startTime).count();
            if(elapsed > 19.3) break;
        }
        HeapNode hn = heapq.top(); heapq.pop();
        if(hn.u<0 || hn.u>=N || hn.v<0 || hn.v>=N) continue;
        if(!VActive[hn.u] || !VActive[hn.v]) continue;
        if(Ver[hn.u] != hn.vu || Ver[hn.v] != hn.vv) continue;
        int keep=hn.u, kill=hn.v; Vec3 np; double qerr=0;
        if(!chooseCollapse(keep, kill, np, star, qerr)) continue;
        if(activeV <= target && qerr > planarLimit) break;
        doCollapse(keep, kill, np, star);
        collapses++;
    }
}

static void saveMesh(){
    vector<int> id(N, -1); int nv=0, nf=0;
    for(int i=0;i<N;i++) if(VActive[i]) id[i]=nv++;
    for(int i=0;i<M;i++) if(F[i].active) nf++;
    static char outbuf[1<<20]; setvbuf(stdout, outbuf, _IOFBF, sizeof(outbuf));
    printf("%d %d\n", nv, nf);
    for(int i=0;i<N;i++) if(VActive[i]) printf("v %.10g %.10g %.10g\n", P[i].x, P[i].y, P[i].z);
    for(int i=0;i<M;i++) if(F[i].active){
        int a=id[F[i].v[0]], b=id[F[i].v[1]], c=id[F[i].v[2]];
        if(a>=0 && b>=0 && c>=0) printf("f %d %d %d\n", a+1,b+1,c+1);
    }
}

int main(){
    simplify();
    saveMesh();
    return 0;
}
