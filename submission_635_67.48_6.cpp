#include <bits/stdc++.h>
using namespace std;

struct Vec3 {
    double x,y,z;
    Vec3():x(0),y(0),z(0){}
    Vec3(double X,double Y,double Z):x(X),y(Y),z(Z){}
    Vec3 operator+(const Vec3& o) const { return {x+o.x,y+o.y,z+o.z}; }
    Vec3 operator-(const Vec3& o) const { return {x-o.x,y-o.y,z-o.z}; }
    Vec3 operator*(double s) const { return {x*s,y*s,z*s}; }
    Vec3 operator/(double s) const { return {x/s,y/s,z/s}; }
    Vec3& operator+=(const Vec3& o){x+=o.x;y+=o.y;z+=o.z;return *this;}
};
static inline double dotv(const Vec3&a,const Vec3&b){return a.x*b.x+a.y*b.y+a.z*b.z;}
static inline Vec3 crossv(const Vec3&a,const Vec3&b){return {a.y*b.z-a.z*b.y,a.z*b.x-a.x*b.z,a.x*b.y-a.y*b.x};}
static inline double norm2(const Vec3&a){return dotv(a,a);} 
static inline double normv(const Vec3&a){return sqrt(norm2(a));}

struct Face { int a,b,c; };
struct Quadric {
    double q[10];
    Quadric(){ memset(q,0,sizeof(q)); }
    inline void addPlane(double a,double b,double c,double d,double w){
        q[0]+=w*a*a; q[1]+=w*a*b; q[2]+=w*a*c; q[3]+=w*a*d;
        q[4]+=w*b*b; q[5]+=w*b*c; q[6]+=w*b*d;
        q[7]+=w*c*c; q[8]+=w*c*d; q[9]+=w*d*d;
    }
    inline void add(const Quadric& o){ for(int i=0;i<10;i++) q[i]+=o.q[i]; }
    inline double eval(const Vec3&p) const {
        double x=p.x,y=p.y,z=p.z;
        return q[0]*x*x + 2*q[1]*x*y + 2*q[2]*x*z + 2*q[3]*x
             + q[4]*y*y + 2*q[5]*y*z + 2*q[6]*y
             + q[7]*z*z + 2*q[8]*z + q[9];
    }
};
static inline Quadric sumQ(const Quadric&a,const Quadric&b){ Quadric r=a; r.add(b); return r; }

struct Ball {
    Vec3 c; double r;
};
static Ball mergeBall(const Ball&A,const Ball&B){
    Vec3 d = B.c - A.c; double dist = normv(d);
    if(A.r >= dist + B.r) return A;
    if(B.r >= dist + A.r) return B;
    if(dist < 1e-30) return {A.c, max(A.r,B.r)};
    double nr = 0.5*(dist + A.r + B.r);
    double t = (nr - A.r) / dist;
    return {A.c + d*t, nr};
}

struct Operation { int dead, keep; Vec3 p; };

static int N, M, ORIGN;
static vector<Vec3> P, OrigP, BaseP;
static vector<Face> F, Forig, OrigF;
static vector<Ball> InitialBS;
static vector<int> InitialCnt;
static vector<unsigned char> faceAlive, vertAlive;
static vector<Quadric> Q;
static vector<Ball> BS;
static vector<Vec3> Cent;
static vector<int> Cnt;
static vector<vector<int>> adjv, incf;
static vector<int> ver;
static vector<Operation> ops;
static int aliveV, aliveF;
static double bboxDiag, hausR, hausR2, RoughEst=-1.0;
static chrono::steady_clock::time_point startTime;

static inline bool faceContains(const Face& f,int v){ return f.a==v || f.b==v || f.c==v; }
static inline bool faceContains2(const Face& f,int u,int v){ return faceContains(f,u) && faceContains(f,v); }
static inline int faceOpposite(const Face& f,int u,int v){ if(f.a!=u && f.a!=v) return f.a; if(f.b!=u && f.b!=v) return f.b; return f.c; }
static inline bool isAliveVertex(int v){ return v>=0 && vertAlive[v]; }

static inline bool hasAdj(int u,int v){
    const auto &a = adjv[u];
    for(int x:a) if(x==v) return true;
    return false;
}
static inline void addAdj(int u,int v){
    if(u==v || !vertAlive[u] || !vertAlive[v]) return;
    auto &a=adjv[u];
    for(int x:a) if(x==v) return;
    a.push_back(v);
}
static inline void remAdj(int u,int v){
    auto &a=adjv[u];
    for(size_t i=0;i<a.size();++i) if(a[i]==v){ a[i]=a.back(); a.pop_back(); return; }
}

static vector<char> slurp_stdin() {
    vector<char> buf;
    buf.reserve(1<<27);
    char chunk[1<<16]; size_t n;
    while((n=fread(chunk,1,sizeof(chunk),stdin))>0) buf.insert(buf.end(), chunk, chunk+n);
    buf.push_back('\0'); return buf;
}
static void load_obj(){
    vector<char> buf=slurp_stdin(); char* p=buf.data();
    N=(int)strtol(p,&p,10); M=(int)strtol(p,&p,10); ORIGN=N;
    P.resize(N); F.resize(M); Forig.resize(M);
    double mnx=1e100,mny=1e100,mnz=1e100,mxx=-1e100,mxy=-1e100,mxz=-1e100;
    for(int i=0;i<N;i++){
        while(*p==' '||*p=='\n'||*p=='\r'||*p=='\t') ++p; ++p;
        double x=strtod(p,&p), y=strtod(p,&p), z=strtod(p,&p);
        P[i]={x,y,z};
        mnx=min(mnx,x); mny=min(mny,y); mnz=min(mnz,z);
        mxx=max(mxx,x); mxy=max(mxy,y); mxz=max(mxz,z);
    }
    for(int i=0;i<M;i++){
        while(*p==' '||*p=='\n'||*p=='\r'||*p=='\t') ++p; ++p;
        int a=(int)strtol(p,&p,10)-1, b=(int)strtol(p,&p,10)-1, c=(int)strtol(p,&p,10)-1;
        F[i]={a,b,c}; Forig[i]=F[i];
    }
    OrigP=P; OrigF=F;
    double dx=mxx-mnx,dy=mxy-mny,dz=mxz-mnz;
    bboxDiag=sqrt(dx*dx+dy*dy+dz*dz);
    hausR=0.05*bboxDiag*0.9990;
    hausR2=hausR*hausR;
}

static void build_structures(){
    faceAlive.assign(M,1); vertAlive.assign(N,1); aliveV=N; aliveF=M;
    Q.assign(N, Quadric()); BS.resize(N); Cent.resize(N); Cnt.assign(N,1);
    vector<int> deg(N,0), ideg(N,0);
    for(const auto& f:F){
        ideg[f.a]++; ideg[f.b]++; ideg[f.c]++;
        deg[f.a]+=2; deg[f.b]+=2; deg[f.c]+=2;
    }
    adjv.resize(N); incf.resize(N);
    for(int i=0;i<N;i++){ adjv[i].reserve(deg[i]); incf[i].reserve(ideg[i]+4); BS[i]=(InitialBS.size()==(size_t)N?InitialBS[i]:Ball{P[i],0.0}); Cent[i]=P[i]; if(InitialCnt.size()==(size_t)N) Cnt[i]=InitialCnt[i]; }
    for(int i=0;i<M;i++){
        Face f=F[i];
        incf[f.a].push_back(i); incf[f.b].push_back(i); incf[f.c].push_back(i);
        adjv[f.a].push_back(f.b); adjv[f.a].push_back(f.c);
        adjv[f.b].push_back(f.a); adjv[f.b].push_back(f.c);
        adjv[f.c].push_back(f.a); adjv[f.c].push_back(f.b);
        Vec3 e1=P[f.b]-P[f.a], e2=P[f.c]-P[f.a];
        Vec3 n=crossv(e1,e2); double len=normv(n);
        if(len>0){
            Vec3 un=n/len; double d=-dotv(un,P[f.a]);
            double w=max(len*0.5, 1e-18);
            Q[f.a].addPlane(un.x,un.y,un.z,d,w);
            Q[f.b].addPlane(un.x,un.y,un.z,d,w);
            Q[f.c].addPlane(un.x,un.y,un.z,d,w);
        }
    }
    for(int i=0;i<N;i++){
        auto &a=adjv[i]; sort(a.begin(), a.end()); a.erase(unique(a.begin(), a.end()), a.end());
    }
    ver.assign(N,0);
}

struct PQItem {
    double cost; int u,v,vu,vv;
    bool operator<(const PQItem& o) const { return cost > o.cost; }
};
static priority_queue<PQItem> pq;

static bool solveOptimal(const Quadric& q, Vec3& out){
    double a00=q.q[0], a01=q.q[1], a02=q.q[2];
    double a11=q.q[4], a12=q.q[5];
    double a22=q.q[7];
    double b0=-q.q[3], b1=-q.q[6], b2=-q.q[8];
    double det = a00*(a11*a22-a12*a12) - a01*(a01*a22-a12*a02) + a02*(a01*a12-a11*a02);
    double scale = fabs(a00)+fabs(a11)+fabs(a22)+fabs(a01)+fabs(a02)+fabs(a12)+1.0;
    if(fabs(det) < 1e-14*scale*scale*scale) return false;
    double c00=(a11*a22-a12*a12), c01=-(a01*a22-a02*a12), c02=(a01*a12-a02*a11);
    double c10=c01, c11=(a00*a22-a02*a02), c12=-(a00*a12-a01*a02);
    double c20=c02, c21=c12, c22=(a00*a11-a01*a01);
    out.x=(c00*b0+c01*b1+c02*b2)/det;
    out.y=(c10*b0+c11*b1+c12*b2)/det;
    out.z=(c20*b0+c21*b1+c22*b2)/det;
    return isfinite(out.x)&&isfinite(out.y)&&isfinite(out.z);
}

static inline bool coveredByBall(const Ball& b,const Vec3&p){
    double d2=norm2(p-b.c);
    double lim=hausR-b.r;
    return lim>=-1e-12 && d2 <= lim*lim + 1e-18;
}

static double quickEdgeCost(int u,int v){
    if(!vertAlive[u]||!vertAlive[v]||u==v) return 1e300;
    Ball mb=mergeBall(BS[u],BS[v]);
    if(mb.r > hausR*1.000001) return 1e300;
    Quadric q=sumQ(Q[u],Q[v]);
    Vec3 best; double bestc=1e300;
    Vec3 cand[5]; int nc=0;
    Vec3 opt; if(solveOptimal(q,opt)) cand[nc++]=opt;
    cand[nc++]=P[u]; cand[nc++]=P[v]; cand[nc++]=(P[u]+P[v])*0.5; cand[nc++]=(Cent[u]*Cnt[u] + Cent[v]*Cnt[v])/(double)(Cnt[u]+Cnt[v]);
    double len2=norm2(P[u]-P[v]);
    for(int i=0;i<nc;i++) if(coveredByBall(mb,cand[i])){
        double c=q.eval(cand[i]); if(c<bestc) bestc=c;
    }
    if(bestc>=1e299) return 1e300;
    return bestc + 1e-12*len2;
}
static inline void pushEdge(int u,int v){
    if(u==v || !vertAlive[u] || !vertAlive[v]) return;
    double c=quickEdgeCost(u,v);
    if(c<1e250) pq.push({c,u,v,ver[u],ver[v]});
}

struct CollapseChoice { int dead, keep; Vec3 p; double cost; vector<int> edgeFaces; };

static void collectEdgeFaces(int u,int v, vector<int>& efaces, vector<int>& opps){
    efaces.clear(); opps.clear();
    const vector<int> &lst = (incf[u].size() < incf[v].size()? incf[u] : incf[v]);
    for(int fid: lst){
        if(fid<0 || fid>=M || !faceAlive[fid]) continue;
        Face f=F[fid];
        if(faceContains(f,u) && faceContains(f,v)){
            efaces.push_back(fid); opps.push_back(faceOpposite(f,u,v));
            if(efaces.size()>2) return;
        }
    }
}

static bool duplicateWouldOccur(int dead,int keep,const vector<int>& edgeFaces){
    // For every face (dead,a,b) that will become (keep,a,b), ensure such face is not already present at keep.
    int ef0=edgeFaces[0], ef1=edgeFaces[1];
    for(int fid: incf[dead]){
        if(!faceAlive[fid] || fid==ef0 || fid==ef1) continue;
        Face f=F[fid]; if(!faceContains(f,dead)) continue;
        int a,b;
        if(f.a==dead){a=f.b;b=f.c;} else if(f.b==dead){a=f.a;b=f.c;} else {a=f.a;b=f.b;}
        if(a==keep || b==keep || a==b) return true;
        for(int gid: incf[keep]){
            if(!faceAlive[gid] || gid==ef0 || gid==ef1 || gid==fid) continue;
            Face g=F[gid];
            if(faceContains(g,a) && faceContains(g,b) && faceContains(g,keep)) return true;
        }
    }
    return false;
}

static bool localGeometryOK(int u,int v,const Vec3&p,const vector<int>& edgeFaces,double &normalPenalty){
    normalPenalty=0.0;
    int ef0=edgeFaces[0], ef1=edgeFaces[1];
    vector<int> touched;
    touched.reserve(incf[u].size()+incf[v].size());
    for(int fid: incf[u]) if(faceAlive[fid] && fid!=ef0 && fid!=ef1) touched.push_back(fid);
    for(int fid: incf[v]) if(faceAlive[fid] && fid!=ef0 && fid!=ef1){
        bool seen=false; for(int x:touched) if(x==fid){seen=true;break;} if(!seen) touched.push_back(fid);
    }
    const double epsArea2 = 1e-28;
    for(int fid: touched){
        Face f=F[fid];
        if(!faceContains(f,u) && !faceContains(f,v)) continue;
        Vec3 A=P[f.a], B=P[f.b], C=P[f.c];
        Vec3 oldN=crossv(B-A,C-A); double oldL2=norm2(oldN);
        if(f.a==u || f.a==v) A=p;
        if(f.b==u || f.b==v) B=p;
        if(f.c==u || f.c==v) C=p;
        Vec3 newN=crossv(B-A,C-A); double newL2=norm2(newN);
        if(newL2 <= epsArea2) return false;
        double dl=dotv(oldN,newN);
        if(dl <= -1e-18) return false;
        double denom=sqrt(max(oldL2*newL2,1e-300));
        double cs=dl/denom;
        if(cs < 0.03) return false;
        if(cs<1.0){
            double area=sqrt(oldL2)*0.5;
            normalPenalty += area * (1.0 - max(-1.0,min(1.0,cs)));
        }
    }
    return true;
}

static bool findBestCollapse(int u,int v, CollapseChoice &res){
    if(!vertAlive[u]||!vertAlive[v]||u==v||!hasAdj(u,v)) return false;
    vector<int> efaces, opps; collectEdgeFaces(u,v,efaces,opps);
    if(efaces.size()!=2) return false;
    if(opps[0]==opps[1] || !vertAlive[opps[0]] || !vertAlive[opps[1]]) return false;
    int common=0;
    for(int n: adjv[u]) if(n!=v && vertAlive[n] && hasAdj(v,n)){
        if(n==opps[0] || n==opps[1]) common++;
        else return false;
    }
    if(common!=2) return false;

    int dead=u, keep=v;
    if(incf[dead].size() > incf[keep].size()) { dead=v; keep=u; }
    if(duplicateWouldOccur(dead,keep,efaces)) return false;

    Ball mb=mergeBall(BS[u],BS[v]);
    if(mb.r > hausR*1.0000001) return false;
    Quadric q=sumQ(Q[u],Q[v]);
    vector<Vec3> cand; cand.reserve(8);
    Vec3 opt; if(solveOptimal(q,opt)) cand.push_back(opt);
    Vec3 pu=P[u], pv=P[v];
    cand.push_back((pu+pv)*0.5);
    // optimal point constrained to the current edge segment
    Vec3 d=pv-pu;
    double dd=norm2(d);
    if(dd>1e-30){
        // q(pu + t d) = at^2 + bt + c
        double eps=1e-6;
        double c0=q.eval(pu), c1=q.eval(pu+d*eps), c_1=q.eval(pu-d*eps);
        double b=(c1-c_1)/(2*eps);
        double a=(c1+c_1-2*c0)/(2*eps*eps);
        if(a>1e-30){ double t=-b/(2*a); if(t>0&&t<1) cand.push_back(pu+d*t); }
    }
    cand.push_back(P[u]); cand.push_back(P[v]);
    cand.push_back((Cent[u]*Cnt[u] + Cent[v]*Cnt[v])/(double)(Cnt[u]+Cnt[v]));
    cand.push_back(mb.c);

    double best=1e300; Vec3 bestp;
    for(const Vec3& p: cand){
        if(!isfinite(p.x)||!isfinite(p.y)||!isfinite(p.z)) continue;
        if(!coveredByBall(mb,p)) continue;
        double np=0;
        if(!localGeometryOK(u,v,p,efaces,np)) continue;
        double len2=norm2(P[u]-P[v]);
        double c=q.eval(p);
        if(c<0 && c>-1e-18) c=0;
        // normal-map preservation is deliberately important; depth/geometric error is in q.
        c += 0.035*np + 1e-12*len2;
        // prefer positions near the cluster center to use Hausdorff budget safely
        double slack=max(0.0, hausR - mb.r);
        if(slack>1e-12) c += 1e-10 * norm2(p-mb.c) / (slack*slack);
        if(c<best){ best=c; bestp=p; }
    }
    if(best>=1e250) return false;
    res.dead=dead; res.keep=keep; res.p=bestp; res.cost=best; res.edgeFaces=efaces;
    return true;
}

static void cleanIncident(int v){
    auto &lst=incf[v];
    size_t w=0;
    for(size_t i=0;i<lst.size();++i){ int fid=lst[i]; if(fid>=0&&fid<M&&faceAlive[fid]&&faceContains(F[fid],v)) lst[w++]=fid; }
    lst.resize(w);
    sort(lst.begin(), lst.end()); lst.erase(unique(lst.begin(), lst.end()), lst.end());
}
static void cleanAdj(int v){
    auto &a=adjv[v]; size_t w=0;
    for(size_t i=0;i<a.size();++i){ int n=a[i]; if(n!=v && n>=0 && n<N && vertAlive[n]) a[w++]=n; }
    a.resize(w); sort(a.begin(),a.end()); a.erase(unique(a.begin(),a.end()),a.end());
}

static void applyCollapse(const CollapseChoice& ch){
    int dead=ch.dead, keep=ch.keep;
    int ef0=ch.edgeFaces[0], ef1=ch.edgeFaces[1];
    if(faceAlive[ef0]){ faceAlive[ef0]=0; aliveF--; }
    if(faceAlive[ef1]){ faceAlive[ef1]=0; aliveF--; }

    // rewrite faces incident to dead (except the two deleted edge faces)
    for(int fid: incf[dead]){
        if(!faceAlive[fid] || fid==ef0 || fid==ef1) continue;
        Face &f=F[fid];
        if(f.a==dead) f.a=keep;
        if(f.b==dead) f.b=keep;
        if(f.c==dead) f.c=keep;
        incf[keep].push_back(fid);
    }

    // update one-ring adjacency
    vector<int> dnei = adjv[dead];
    for(int n: dnei){
        if(n==keep || !vertAlive[n]) continue;
        remAdj(n,dead);
        addAdj(n,keep);
        addAdj(keep,n);
    }
    remAdj(keep,dead);
    adjv[dead].clear(); incf[dead].clear();

    P[keep]=ch.p;
    Q[keep].add(Q[dead]);
    BS[keep]=mergeBall(BS[keep],BS[dead]);
    Cent[keep]=(Cent[keep]*Cnt[keep] + Cent[dead]*Cnt[dead])/(double)(Cnt[keep]+Cnt[dead]);
    Cnt[keep]+=Cnt[dead];
    vertAlive[dead]=0; aliveV--; ver[keep]++; ver[dead]++;
    ops.push_back({dead,keep,ch.p});

    if(incf[keep].size() > (adjv[keep].size()+8)*8) cleanIncident(keep);
    if(adjv[keep].size() > 2000) cleanAdj(keep);

    for(int n: adjv[keep]) if(vertAlive[n]) pushEdge(keep,n);
}

static double targetRatioForN(int n){
    if(n<=10) return 8.0/(double)n;
    if(n<=6000) return 0.060;
    if(n<=30000) return 0.065;
    if(n<=60000) return 0.070;
    if(n<=450000) return 0.072;
    return 0.075;
}
static double minRatioForN(int n){
    if(n<=10) return 8.0/(double)n;
    if(n<=6000) return 0.010;
    if(n<=30000) return 0.010;
    if(n<=60000) return 0.010;
    double r = (RoughEst<0?0.0:RoughEst);
    if(n<=450000){
        if(r>0.0080) return 0.18;
        if(r>0.0025) return 0.060;
        if(r>0.0001) return 0.015;
        return 0.0010;
    }
    if(r>0.0080) return 0.20;
    if(r>0.0025) return 0.070;
    if(r>0.0001) return 0.018;
    return 0.0010;
}


struct ClusterResult {
    vector<Vec3> V;
    vector<Face> F;
    vector<Ball> balls;
    vector<int> cnt;
};
static inline uint64_t edgeKey32(int a,int b){ if(a>b) swap(a,b); return (uint64_t)(uint32_t)a<<32 | (uint32_t)b; }
static double estimateRoughness(){
    size_t FM=OrigF.size();
    if(FM==0) return 0.0;
    vector<Vec3> norms(FM);
    for(size_t i=0;i<FM;i++){
        Face f=OrigF[i];
        Vec3 A=OrigP[f.a], B=OrigP[f.b], C=OrigP[f.c];
        Vec3 n=crossv(B-A,C-A); double l=normv(n);
        norms[i]=(l>0? n/l : Vec3(0,0,1));
    }
    unordered_map<uint64_t,int> em;
    em.reserve(FM*3/2+16);
    double sum=0.0; long long cnt=0;
    for(size_t i=0;i<FM;i++){
        Face f=OrigF[i];
        int vs[3]={f.a,f.b,f.c};
        for(int e=0;e<3;e++){
            uint64_t key=edgeKey32(vs[e],vs[(e+1)%3]);
            auto it=em.find(key);
            if(it==em.end()) em.emplace(key,(int)i);
            else {
                double d=dotv(norms[i], norms[it->second]);
                if(d<-1) d=-1; if(d>1) d=1;
                sum += 1.0 - d;
                cnt++;
            }
        }
    }
    return cnt? sum/(double)cnt : 0.0;
}
static bool buildClusterCandidate(double h,double ox,double oy,double oz,ClusterResult &cr){
    if(h<=0) return false;
    const double base = 1.0000000001;
    vector<int> rid(ORIGN);
    struct Acc { double sx,sy,sz; int cnt; };
    vector<Acc> acc; acc.reserve(min(ORIGN, (int)(ORIGN*0.6)+10));
    unordered_map<uint64_t,int> mp; mp.reserve((size_t)(ORIGN*1.25));
    auto keyOf=[&](const Vec3&p)->uint64_t{
        int ix=(int)floor((p.x + base + ox*h)/h);
        int iy=(int)floor((p.y + base + oy*h)/h);
        int iz=(int)floor((p.z + base + oz*h)/h);
        return ((uint64_t)(uint32_t)ix<<42) ^ ((uint64_t)(uint32_t)iy<<21) ^ (uint64_t)(uint32_t)iz;
    };
    for(int i=0;i<ORIGN;i++){
        uint64_t key=keyOf(OrigP[i]);
        auto it=mp.find(key); int id;
        if(it==mp.end()){
            id=(int)acc.size(); mp.emplace(key,id); acc.push_back({0,0,0,0});
        } else id=it->second;
        rid[i]=id; acc[id].sx+=OrigP[i].x; acc[id].sy+=OrigP[i].y; acc[id].sz+=OrigP[i].z; acc[id].cnt++;
    }
    int K=(int)acc.size();
    vector<Vec3> np(K);
    vector<int> nc(K);
    for(int i=0;i<K;i++){ double inv=1.0/acc[i].cnt; np[i]={acc[i].sx*inv,acc[i].sy*inv,acc[i].sz*inv}; nc[i]=acc[i].cnt; }
    vector<double> rad2(K,0.0);
    for(int i=0;i<ORIGN;i++){
        int id=rid[i]; double d2=norm2(OrigP[i]-np[id]); if(d2>rad2[id]) rad2[id]=d2;
    }
    double coverLimit=hausR*0.999;
    for(int i=0;i<K;i++) if(rad2[i] > coverLimit*coverLimit + 1e-18) return false;

    vector<Face> nf; nf.reserve(OrigF.size());
    vector<int> used(K,0);
    unordered_map<uint64_t, unsigned char> ec; ec.reserve((size_t)OrigF.size()*2);
    const double epsArea2=1e-28;
    for(size_t i=0;i<OrigF.size();i++){
        Face of=OrigF[i]; int a=rid[of.a], b=rid[of.b], c=rid[of.c];
        if(a==b || b==c || a==c) continue;
        Vec3 A=np[a], B=np[b], C=np[c];
        Vec3 nn=crossv(B-A,C-A); double n2=norm2(nn);
        if(n2<=epsArea2) return false;
        Vec3 OA=OrigP[of.a], OB=OrigP[of.b], OC=OrigP[of.c];
        Vec3 on=crossv(OB-OA,OC-OA);
        if(dotv(on,nn) <= 0.0) return false;
        uint64_t e1=edgeKey32(a,b), e2=edgeKey32(b,c), e3=edgeKey32(c,a);
        unsigned char &c1=ec[e1]; if(++c1>2) return false;
        unsigned char &c2=ec[e2]; if(++c2>2) return false;
        unsigned char &c3=ec[e3]; if(++c3>2) return false;
        nf.push_back({a,b,c}); used[a]++; used[b]++; used[c]++;
    }
    if(nf.empty()) return false;
    for(auto &kv: ec) if(kv.second!=2) return false;
    for(int i=0;i<K;i++) if(used[i]==0) return false;
    cr.V.swap(np); cr.F.swap(nf); cr.cnt.swap(nc); cr.balls.resize(K);
    for(int i=0;i<K;i++) cr.balls[i]={cr.V[i], sqrt(rad2[i])};
    return true;
}

static bool preprocessCluster(){
    if(ORIGN < 80000) return false;
    vector<double> factors={0.70,0.60,0.50,0.42,0.35,0.29,0.24,0.19,0.15,0.12,0.095,0.075,0.058,0.048,0.040,0.034,0.029,0.025};
    vector<array<double,3>> offs={{0,0,0},{0.5,0.5,0.5},{0,0.5,0.5},{0.5,0,0.5},{0.5,0.5,0},{0.25,0.25,0.25},{0.75,0.25,0.5},{0.25,0.75,0.5},{0.5,0.25,0.75}};
    ClusterResult chosen; bool have=false;
    RoughEst = estimateRoughness();
    double rough = RoughEst;
    double minStartRatio;
    if(rough < 1e-4) minStartRatio = 0.08;
    else if(rough < 0.0025) minStartRatio = 0.15;
    else if(rough < 0.008) minStartRatio = 0.24;
    else minStartRatio = 0.32;
    double minAcceptRatio = max(0.045, minStartRatio*0.72);
    const double maxStartRatio = 0.58;
    for(double fac: factors){
        double h=hausR*fac;
        if(h<1e-7) break;
        for(auto of: offs){
            ClusterResult cr;
            if(!buildClusterCandidate(h,of[0],of[1],of[2],cr)) continue;
            double ratio=(double)cr.V.size()/(double)ORIGN;
            if(ratio < minStartRatio) { have=true; chosen=std::move(cr); continue; }
            if(ratio <= maxStartRatio){ chosen=std::move(cr); have=true; goto done; }
        }
    }
done:
    if(!have) return false;
    double ratio=(double)chosen.V.size()/(double)ORIGN;
    if(ratio>maxStartRatio || ratio<minAcceptRatio) return false;
    P=std::move(chosen.V); F=std::move(chosen.F); Forig=F; InitialBS=std::move(chosen.balls); InitialCnt=std::move(chosen.cnt);
    N=(int)P.size(); M=(int)F.size();
    return true;
}

static void simplify_mesh(){
    BaseP=P;
    if(N<=3) return;
    build_structures();
    for(int u=0;u<N;u++) for(int v: adjv[u]) if(u<v) pushEdge(u,v);
    int minAlive = max(4, (int)ceil(minRatioForN(ORIGN)*ORIGN));
    if(ORIGN<=10) minAlive=8;
    if(minAlive>=N) minAlive=N;
    const double timeLimit = (N>700000? 16.2 : (N>200000? 17.0 : 18.2));
    long long pops=0;
    while(aliveV>minAlive && !pq.empty()){
        PQItem it=pq.top(); pq.pop(); pops++;
        if(!vertAlive[it.u]||!vertAlive[it.v]||ver[it.u]!=it.vu||ver[it.v]!=it.vv) continue;
        if(!hasAdj(it.u,it.v)) continue;
        CollapseChoice ch;
        if(!findBestCollapse(it.u,it.v,ch)) continue;
        if(ch.cost > it.cost*1.0001 + 1e-16){
            pq.push({ch.cost,it.u,it.v,ver[it.u],ver[it.v]});
            continue;
        }
        applyCollapse(ch);
        if((pops & 4095)==0){
            double elapsed=chrono::duration<double>(chrono::steady_clock::now()-startTime).count();
            if(elapsed>timeLimit) break;
        }
    }
}

// ---------- optional low-resolution perceptual guard ----------
struct Image {
    int S;
    vector<float> r,g,b,d;
    vector<unsigned char> fg;
    Image(int s=0):S(s),r(s*s,127.5f),g(s*s,127.5f),b(s*s,127.5f),d(s*s,255.0f),fg(s*s,0){}
};
struct DSUState {
    vector<int> parent;
    vector<Vec3> pos;
    DSUState(int n=0):parent(n),pos(BaseP.empty()?P:BaseP){ iota(parent.begin(),parent.end(),0); }
    int find(int x){ while(parent[x]!=x){ parent[x]=parent[parent[x]]; x=parent[x]; } return x; }
    void applyOp(const Operation& op){ int d=find(op.dead), k=find(op.keep); if(d!=k){ parent[d]=k; pos[k]=op.p; } }
};

static inline void viewBasis(int view, Vec3& dir, Vec3& right, Vec3& up){
    switch(view){
        case 0: dir={1,0,0}; right={0,1,0}; up={0,0,1}; break;
        case 1: dir={-1,0,0}; right={0,-1,0}; up={0,0,1}; break;
        case 2: dir={0,1,0}; right={-1,0,0}; up={0,0,1}; break;
        case 3: dir={0,-1,0}; right={1,0,0}; up={0,0,1}; break;
        case 4: dir={0,0,1}; right={1,0,0}; up={0,1,0}; break;
        default: dir={0,0,-1}; right={1,0,0}; up={0,-1,0}; break;
    }
}

static void rasterFace(Image& img,const Vec3&A,const Vec3&B,const Vec3&C,const Vec3&Nrm,int view){
    const double D=2.5; int S=img.S; double scale=S/1024.0; double fx=800.0*scale, fy=800.0*scale, cc=0.5*S;
    Vec3 dir,right,up; viewBasis(view,dir,right,up);
    auto proj=[&](const Vec3&p,double&u,double&v,double&dep){
        dep = D - dotv(p,dir);
        double sx=dotv(p,right), sy=dotv(p,up);
        u=fx*sx/dep + cc; v=fy*sy/dep + cc;
    };
    double u0,v0,z0,u1,v1,z1,u2,v2,z2;
    proj(A,u0,v0,z0); proj(B,u1,v1,z1); proj(C,u2,v2,z2);
    if(z0<=0||z1<=0||z2<=0) return;
    double area=(u1-u0)*(v2-v0)-(v1-v0)*(u2-u0);
    if(fabs(area)<1e-12) return;
    int minx=max(0,(int)floor(min({u0,u1,u2})-0.5));
    int maxx=min(S-1,(int)ceil(max({u0,u1,u2})-0.5));
    int miny=max(0,(int)floor(min({v0,v1,v2})-0.5));
    int maxy=min(S-1,(int)ceil(max({v0,v1,v2})-0.5));
    if(minx>maxx||miny>maxy) return;
    double invArea=1.0/area;
    float nr=(float)((Nrm.x+1.0)*127.5), ng=(float)((Nrm.y+1.0)*127.5), nb=(float)((Nrm.z+1.0)*127.5);
    for(int y=miny;y<=maxy;y++){
        double py=y+0.5;
        for(int x=minx;x<=maxx;x++){
            double px=x+0.5;
            double w0=((u1-px)*(v2-py)-(v1-py)*(u2-px))*invArea;
            double w1=((u2-px)*(v0-py)-(v2-py)*(u0-px))*invArea;
            double w2=1.0-w0-w1;
            if(w0>=-1e-9 && w1>=-1e-9 && w2>=-1e-9){
                double dep=1.0/(w0/z0+w1/z1+w2/z2);
                int id=y*S+x;
                if(dep<img.d[id]){ img.d[id]=(float)dep; img.r[id]=nr; img.g[id]=ng; img.b[id]=nb; img.fg[id]=1; }
            }
        }
    }
}

static Image renderOriginal(int S,int view){
    Image img(S);
    for(size_t i=0;i<OrigF.size();i++){
        Face f=OrigF[i]; Vec3 A=OrigP[f.a],B=OrigP[f.b],C=OrigP[f.c];
        Vec3 n=crossv(B-A,C-A); double l=normv(n); if(l<=0) continue; n=n/l;
        rasterFace(img,A,B,C,n,view);
    }
    return img;
}
static Image renderCandidate(int S,int view,DSUState& ds){
    Image img(S);
    for(int i=0;i<M;i++){
        Face f=Forig[i]; int a=ds.find(f.a), b=ds.find(f.b), c=ds.find(f.c);
        if(a==b||b==c||a==c) continue;
        Vec3 A=ds.pos[a],B=ds.pos[b],C=ds.pos[c];
        Vec3 n=crossv(B-A,C-A); double l=normv(n); if(l<=1e-30) continue; n=n/l;
        rasterFace(img,A,B,C,n,view);
    }
    return img;
}

static double ssimOne(const vector<float>& X,const vector<float>& Y,const vector<unsigned char>& fgX,const vector<unsigned char>& fgY,int S,int Rw){
    int Np=S*S; (void)Np;
    vector<double> ix((S+1)*(S+1),0), iy((S+1)*(S+1),0), ix2((S+1)*(S+1),0), iy2((S+1)*(S+1),0), ixy((S+1)*(S+1),0);
    auto at=[&](vector<double>& I,int y,int x)->double&{ return I[y*(S+1)+x]; };
    for(int y=0;y<S;y++){
        double sx=0,sy=0,sx2=0,sy2=0,sxy=0;
        for(int x=0;x<S;x++){
            int id=y*S+x; double vx=X[id], vy=Y[id];
            sx+=vx; sy+=vy; sx2+=vx*vx; sy2+=vy*vy; sxy+=vx*vy;
            at(ix,y+1,x+1)=at(ix,y,x+1)+sx;
            at(iy,y+1,x+1)=at(iy,y,x+1)+sy;
            at(ix2,y+1,x+1)=at(ix2,y,x+1)+sx2;
            at(iy2,y+1,x+1)=at(iy2,y,x+1)+sy2;
            at(ixy,y+1,x+1)=at(ixy,y,x+1)+sxy;
        }
    }
    auto sumRect=[&](vector<double>& I,int y0,int x0,int y1,int x1){
        return I[y1*(S+1)+x1]-I[y0*(S+1)+x1]-I[y1*(S+1)+x0]+I[y0*(S+1)+x0];
    };
    const double c1=6.5025, c2=58.5225;
    double total=0; int cnt=0;
    for(int y=0;y<S;y++) for(int x=0;x<S;x++){
        int id=y*S+x; if(!fgX[id] && !fgY[id]) continue;
        int y0=max(0,y-Rw), y1=min(S,y+Rw+1), x0=max(0,x-Rw), x1=min(S,x+Rw+1);
        double area=(y1-y0)*(x1-x0);
        double sx=sumRect(ix,y0,x0,y1,x1), sy=sumRect(iy,y0,x0,y1,x1);
        double sx2=sumRect(ix2,y0,x0,y1,x1), sy2=sumRect(iy2,y0,x0,y1,x1), sxy=sumRect(ixy,y0,x0,y1,x1);
        double mx=sx/area, my=sy/area;
        double vx=max(0.0,sx2/area-mx*mx), vy=max(0.0,sy2/area-my*my), cov=sxy/area-mx*my;
        double val=((2*mx*my+c1)*(2*cov+c2))/((mx*mx+my*my+c1)*(vx+vy+c2));
        if(isfinite(val)){ total+=val; cnt++; }
    }
    if(cnt==0) return 1.0;
    return max(0.0,min(1.0,total/cnt));
}
static double compareSSIM(const vector<Image>& orig, vector<Image>& cand){
    double final=0;
    for(int v=0;v<6;v++){
        int rw=max(1, (int)floor(5.0*orig[v].S/1024.0));
        double nr=ssimOne(orig[v].r,cand[v].r,orig[v].fg,cand[v].fg,orig[v].S,rw);
        double ng=ssimOne(orig[v].g,cand[v].g,orig[v].fg,cand[v].fg,orig[v].S,rw);
        double nb=ssimOne(orig[v].b,cand[v].b,orig[v].fg,cand[v].fg,orig[v].S,rw);
        double nd=ssimOne(orig[v].d,cand[v].d,orig[v].fg,cand[v].fg,orig[v].S,rw);
        final += 0.5*((nr+ng+nb)/3.0) + 0.5*nd;
    }
    return final/6.0;
}

static int choosePrefix(){
    int maxOps=ops.size();
    if(ORIGN<=10) return maxOps;
    double elapsed=chrono::duration<double>(chrono::steady_clock::now()-startTime).count();
    if(elapsed > 18.0) return maxOps;
    int S = (ORIGN>700000? 448 : (ORIGN>120000? 512 : 512));
    if(ORIGN>120000 && RoughEst>0.0080) S=384;
    vector<Image> orig; orig.reserve(6);
    for(int v=0;v<6;v++) orig.push_back(renderOriginal(S,v));
    vector<double> ratios;
    // from safer to more aggressive; the last passing level is selected
    if(ORIGN<=6000) ratios={0.600,0.450,0.320,0.240,0.180,0.135,0.105,0.083,0.068,0.058,0.050,0.044,0.036,0.030,0.024,0.018,0.012};
    else if(ORIGN<=30000) ratios={0.600,0.450,0.320,0.240,0.180,0.135,0.105,0.085,0.073,0.066,0.064,0.062,0.060,0.058,0.056,0.054,0.052,0.050,0.044,0.038,0.032,0.026,0.020,0.014};
    else if(ORIGN<=60000) ratios={0.600,0.450,0.350,0.300,0.250,0.200,0.160,0.130,0.105,0.090,0.078,0.068,0.060,0.054,0.048,0.042,0.036,0.030,0.024,0.018,0.012};
    else if(ORIGN<=450000) ratios={0.600,0.450,0.350,0.280,0.220,0.170,0.135,0.110,0.095,0.083,0.073,0.064,0.057,0.052,0.047,0.042,0.037,0.032,0.028,0.024,0.021,0.018,0.015,0.012,0.010,0.008,0.006,0.0045,0.0035,0.0025,0.0018,0.0012};
    else ratios={0.600,0.450,0.350,0.280,0.220,0.170,0.140,0.115,0.100,0.088,0.078,0.068,0.060,0.054,0.048,0.042,0.036,0.032,0.028,0.024,0.021,0.018,0.015,0.012,0.010,0.008,0.006,0.0045,0.0035,0.0025,0.0018,0.0012};
    DSUState ds(N);
    int applied=0, bestK=0, lastRendered=-1;
    bool anyPass=false;
    double threshold = 0.905;
    for(double r: ratios){
        int target=max(4,(int)ceil(r*ORIGN));
        int k=min(maxOps, max(0,N-target));
        if(k<applied) continue;
        if(k==lastRendered) continue;
        lastRendered=k;
        while(applied<k){ ds.applyOp(ops[applied]); applied++; }
        vector<Image> cand; cand.reserve(6);
        for(int v=0;v<6;v++) cand.push_back(renderCandidate(S,v,ds));
        double sc=compareSSIM(orig,cand);
        if(sc>=threshold){ bestK=k; anyPass=true; } else break;
        double now=chrono::duration<double>(chrono::steady_clock::now()-startTime).count();
        if(now>19.0) break;
    }
    if(!anyPass){
        int safeTarget=max(4,(int)ceil(0.60*ORIGN));
        bestK=min(maxOps,max(0,N-safeTarget));
    }
    return bestK;
}

static void outputPrefix(int k){
    DSUState ds(N);
    vector<Ball> rb(N);
    bool hasInit = (InitialBS.size()==(size_t)N);
    for(int i=0;i<N;i++) rb[i]=hasInit?InitialBS[i]:Ball{(BaseP.empty()?P[i]:BaseP[i]),0.0};
    for(int i=0;i<k;i++){
        int d=ds.find(ops[i].dead), kk=ds.find(ops[i].keep);
        if(d!=kk){ rb[kk]=mergeBall(rb[kk],rb[d]); ds.parent[d]=kk; ds.pos[kk]=ops[i].p; }
    }
    vector<int> root(N);
    for(int i=0;i<N;i++) root[i]=ds.find(i);
    vector<int> outId(N,-1);
    vector<int> outRoot; outRoot.reserve(N-k);
    vector<Vec3> outV; outV.reserve(N-k);
    // Include all surviving representatives to preserve the vertex-set Hausdorff cover.
    for(int i=0;i<N;i++) if(root[i]==i){ outId[i]=(int)outV.size(); outRoot.push_back(i); outV.push_back(ds.pos[i]); }
    if(ORIGN>10){
        double baseJ = max(1e-15, bboxDiag*1e-10);
        for(size_t i=0;i<outV.size();++i){
            uint64_t h=(uint64_t)(i+1)*11995408973635179863ull;
            double dx=((int)((h>> 8)&1023)-511.5)/511.5;
            double dy=((int)((h>>24)&1023)-511.5)/511.5;
            double dz=((int)((h>>40)&1023)-511.5)/511.5;
            int r=outRoot[i];
            double slack = hausR - rb[r].r - normv(outV[i]-rb[r].c);
            if(slack>0){
                double amp = min(baseJ, slack*0.20/1.7320508075688772);
                outV[i].x += amp*dx; outV[i].y += amp*dy; outV[i].z += amp*dz;
            }
        }
    }
    vector<Face> outF; outF.reserve(max(0,M-2*k));
    for(int i=0;i<M;i++){
        Face f=Forig[i]; int a=root[f.a], b=root[f.b], c=root[f.c];
        if(a==b||b==c||a==c) continue;
        outF.push_back({outId[a],outId[b],outId[c]});
    }
    auto makeOutput=[&](int prec){
        string out; out.reserve((size_t)outV.size()*(prec==15?58:42) + (size_t)outF.size()*28 + 64);
        char line[160];
        out.append(line, snprintf(line,sizeof(line), "%zu %zu\n", outV.size(), outF.size()));
        if(prec==15){
            for(const Vec3& p: outV) out.append(line, snprintf(line,sizeof(line), "v %.15g %.15g %.15g\n", p.x,p.y,p.z));
        }else{
            for(const Vec3& p: outV) out.append(line, snprintf(line,sizeof(line), "v %.10g %.10g %.10g\n", p.x,p.y,p.z));
        }
        for(const Face& f: outF) out.append(line, snprintf(line,sizeof(line), "f %d %d %d\n", f.a+1,f.b+1,f.c+1));
        return out;
    };
    string out = makeOutput(15);
    if(out.size() > 100ull*1024ull*1024ull) out = makeOutput(10);
    fwrite(out.data(),1,out.size(),stdout);
}

int main(){
    startTime=chrono::steady_clock::now();
    load_obj();
    preprocessCluster();
    simplify_mesh();
    int k=choosePrefix();
    outputPrefix(k);
    return 0;
}
