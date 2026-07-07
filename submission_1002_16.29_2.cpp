#include<bits/stdc++.h>
#include<unistd.h>
using namespace std;
namespace A{
using namespace std;
#define main main_a

struct Vec3 { double x, y, z; };
static inline Vec3 operator+(const Vec3& a,const Vec3& b){ return {a.x+b.x,a.y+b.y,a.z+b.z}; }
static inline Vec3 operator-(const Vec3& a,const Vec3& b){ return {a.x-b.x,a.y-b.y,a.z-b.z}; }
static inline Vec3 operator*(const Vec3& a,double s){ return {a.x*s,a.y*s,a.z*s}; }
static inline Vec3 operator/(const Vec3& a,double s){ return {a.x/s,a.y/s,a.z/s}; }
static inline double dot3(const Vec3& a,const Vec3& b){ return a.x*b.x+a.y*b.y+a.z*b.z; }
static inline Vec3 cross3(const Vec3& a,const Vec3& b){ return {a.y*b.z-a.z*b.y,a.z*b.x-a.x*b.z,a.x*b.y-a.y*b.x}; }
static inline double norm2(const Vec3& a){ return dot3(a,a); }
static inline double norm3(const Vec3& a){ return sqrt(max(0.0,norm2(a))); }

struct Face { int v[3]; unsigned char alive; };
struct Box3 { double mn[3], mx[3]; };
struct Quadric { double q[10]; };

static int N0=0, M0=0;
static vector<Vec3> P, P_orig;
static vector<Face> F, F_orig;
static vector<unsigned char> aliveV;
static vector<int> verV, csz;
static vector<Box3> bboxV;
static vector<Quadric> Q;
static vector<double> Qw;
static vector<vector<int>> inc;
static int activeV=0, activeF=0;
static double hausTol=0.0, hausTol2=0.0;
static chrono::steady_clock::time_point startTime;

static inline void qzero(Quadric& q){ for(double &x:q.q) x=0.0; }
static inline Quadric qadd(const Quadric& a,const Quadric& b){ Quadric r; for(int i=0;i<10;i++) r.q[i]=a.q[i]+b.q[i]; return r; }
static inline void qaddto(Quadric& a,const Quadric& b){ for(int i=0;i<10;i++) a.q[i]+=b.q[i]; }
static inline double qeval(const Quadric& q,const Vec3& p){
    double x=p.x,y=p.y,z=p.z;
    return q.q[0]*x*x + 2*q.q[1]*x*y + 2*q.q[2]*x*z + 2*q.q[3]*x
         + q.q[4]*y*y + 2*q.q[5]*y*z + 2*q.q[6]*y
         + q.q[7]*z*z + 2*q.q[8]*z + q.q[9];
}
static inline void addPlaneQuadric(Quadric& q,double a,double b,double c,double d,double w){
    q.q[0]+=w*a*a; q.q[1]+=w*a*b; q.q[2]+=w*a*c; q.q[3]+=w*a*d;
    q.q[4]+=w*b*b; q.q[5]+=w*b*c; q.q[6]+=w*b*d;
    q.q[7]+=w*c*c; q.q[8]+=w*c*d; q.q[9]+=w*d*d;
}
static inline bool faceHas(const Face& f,int u){ return f.v[0]==u || f.v[1]==u || f.v[2]==u; }
static inline bool faceHasBoth(const Face& f,int u,int v){ return faceHas(f,u) && faceHas(f,v); }
static inline Vec3 faceNormalRaw(const Face& f){ return cross3(P[f.v[1]]-P[f.v[0]], P[f.v[2]]-P[f.v[0]]); }
static inline bool faceUnitNormal(const Face& f, Vec3& n){ Vec3 r=faceNormalRaw(f); double l=norm3(r); if(l<=1e-30){ n={0,0,0}; return false; } n=r/l; return true; }

static vector<char> slurp_stdin(){
    vector<char> buf; buf.reserve(1u<<27);
    char chunk[1<<16]; size_t n;
    while((n=fread(chunk,1,sizeof(chunk),stdin))>0) buf.insert(buf.end(),chunk,chunk+n);
    buf.push_back('\0'); return buf;
}
static inline void skipWs(char*& p){ while(*p==' '||*p=='\n'||*p=='\r'||*p=='\t') ++p; }
static void loadInput(){
    vector<char> buf=slurp_stdin(); char* p=buf.data();
    N0=(int)strtol(p,&p,10); M0=(int)strtol(p,&p,10);
    if(N0<0) N0=0; if(M0<0) M0=0;
    P.resize(N0); F.resize(M0); vector<int> deg(N0,0);
    for(int i=0;i<N0;i++){
        skipWs(p); if(*p=='v') ++p;
        P[i].x=strtod(p,&p); P[i].y=strtod(p,&p); P[i].z=strtod(p,&p);
    }
    for(int i=0;i<M0;i++){
        skipWs(p); if(*p=='f') ++p;
        int a=(int)strtol(p,&p,10)-1, b=(int)strtol(p,&p,10)-1, c=(int)strtol(p,&p,10)-1;
        F[i].v[0]=a; F[i].v[1]=b; F[i].v[2]=c; F[i].alive=1;
        if((unsigned)a<(unsigned)N0) deg[a]++;
        if((unsigned)b<(unsigned)N0) deg[b]++;
        if((unsigned)c<(unsigned)N0) deg[c]++;
    }
    P_orig=P; F_orig=F;
    inc.assign(N0,{});
    for(int i=0;i<N0;i++) inc[i].reserve(max(1,deg[i]+2));
    for(int i=0;i<M0;i++){
        int a=F[i].v[0],b=F[i].v[1],c=F[i].v[2];
        if((unsigned)a<(unsigned)N0) inc[a].push_back(i);
        if((unsigned)b<(unsigned)N0) inc[b].push_back(i);
        if((unsigned)c<(unsigned)N0) inc[c].push_back(i);
    }
    aliveV.assign(N0,1); verV.assign(N0,0); csz.assign(N0,1); bboxV.resize(N0);
    for(int i=0;i<N0;i++){
        bboxV[i].mn[0]=bboxV[i].mx[0]=P[i].x;
        bboxV[i].mn[1]=bboxV[i].mx[1]=P[i].y;
        bboxV[i].mn[2]=bboxV[i].mx[2]=P[i].z;
    }
    activeV=N0; activeF=M0;
}

static inline uint64_t edgeKey(int a,int b){ if(a>b) swap(a,b); return (uint64_t)(uint32_t)a<<32 | (uint32_t)b; }
static Box3 mergeBox(const Box3& a,const Box3& b){ Box3 r; for(int k=0;k<3;k++){ r.mn[k]=min(a.mn[k],b.mn[k]); r.mx[k]=max(a.mx[k],b.mx[k]); } return r; }
static inline double boxMaxDist2(const Vec3& p,const Box3& b){
    double ans=0.0;
    for(int mask=0;mask<8;mask++){
        double x=(mask&1)?b.mx[0]:b.mn[0], y=(mask&2)?b.mx[1]:b.mn[1], z=(mask&4)?b.mx[2]:b.mn[2];
        double dx=p.x-x, dy=p.y-y, dz=p.z-z; ans=max(ans,dx*dx+dy*dy+dz*dz);
    }
    return ans;
}
static inline bool coversBox(const Vec3& p,const Box3& b){ return boxMaxDist2(p,b) <= hausTol2*0.992*0.992 + 1e-18; }

static void initializeQuadricsAndTolerance(){
    if(N0<=0) return;
    double mnx=P[0].x,mny=P[0].y,mnz=P[0].z,mxx=P[0].x,mxy=P[0].y,mxz=P[0].z;
    for(const Vec3& p:P){ mnx=min(mnx,p.x); mny=min(mny,p.y); mnz=min(mnz,p.z); mxx=max(mxx,p.x); mxy=max(mxy,p.y); mxz=max(mxz,p.z); }
    double dx=mxx-mnx,dy=mxy-mny,dz=mxz-mnz;
    hausTol=0.05*sqrt(dx*dx+dy*dy+dz*dz); if(!(hausTol>0)) hausTol=1e-9; hausTol2=hausTol*hausTol;
    Q.resize(N0); Qw.assign(N0,0.0); for(int i=0;i<N0;i++) qzero(Q[i]);
    for(int i=0;i<M0;i++){
        const Face& f=F[i];
        if((unsigned)f.v[0]>=(unsigned)N0 || (unsigned)f.v[1]>=(unsigned)N0 || (unsigned)f.v[2]>=(unsigned)N0) continue;
        Vec3 a=P[f.v[0]],b=P[f.v[1]],c=P[f.v[2]], cr=cross3(b-a,c-a);
        double dblA=norm3(cr); if(dblA<=1e-30) continue;
        Vec3 n=cr/dblA; double d=-dot3(n,a), area=0.5*dblA;
        double visual=fabs(n.x)+fabs(n.y)+fabs(n.z);
        double w=max(1e-12, area*(0.25+visual));
        Quadric pq; qzero(pq); addPlaneQuadric(pq,n.x,n.y,n.z,d,w);
        for(int k=0;k<3;k++){ qaddto(Q[f.v[k]],pq); Qw[f.v[k]]+=w; }
    }
}

static void compactInc(int u){
    if(!aliveV[u]) return; vector<int>& v=inc[u]; if(v.empty()) return;
    int w=0; for(int id:v){ if((unsigned)id>=(unsigned)F.size()) continue; if(!F[id].alive) continue; if(!faceHas(F[id],u)) continue; v[w++]=id; }
    v.resize(w); sort(v.begin(),v.end()); v.erase(unique(v.begin(),v.end()),v.end());
}
static void collectNeighbors(int u, vector<int>& out){
    out.clear(); if(!aliveV[u]) return; if(inc[u].size()>512) compactInc(u);
    for(int id:inc[u]){ if(!F[id].alive) continue; const Face& f=F[id]; if(!faceHas(f,u)) continue; for(int k=0;k<3;k++){ int w=f.v[k]; if(w!=u && (unsigned)w<(unsigned)N0 && aliveV[w]) out.push_back(w); } }
    sort(out.begin(),out.end()); out.erase(unique(out.begin(),out.end()),out.end());
}
static int edgeFaces(int u,int v,int ids[3]){
    int cnt=0; const vector<int>& a=(inc[u].size()<=inc[v].size())?inc[u]:inc[v];
    for(int id:a){ if(!F[id].alive) continue; const Face& f=F[id]; if(faceHasBoth(f,u,v)){ if(cnt<3) ids[cnt]=id; ++cnt; if(cnt>2) return cnt; } }
    return cnt;
}
static bool currentEdge(int u,int v){ int ids[3]; return edgeFaces(u,v,ids)>0; }
static bool duplicateAfterCollapse(int keep,int kill){
    vector<pair<int,int>> have; have.reserve(inc[keep].size());
    for(int id:inc[keep]){
        if(!F[id].alive) continue; const Face& f=F[id]; if(!faceHas(f,keep)||faceHas(f,kill)) continue;
        int a=-1,b=-1; for(int k=0;k<3;k++) if(f.v[k]!=keep){ if(a<0) a=f.v[k]; else b=f.v[k]; }
        if(a<0||b<0) return true; if(a>b) swap(a,b); have.push_back({a,b});
    }
    sort(have.begin(),have.end()); have.erase(unique(have.begin(),have.end()),have.end());
    for(int id:inc[kill]){
        if(!F[id].alive) continue; const Face& f=F[id]; if(!faceHas(f,kill)||faceHas(f,keep)) continue;
        int a=-1,b=-1; for(int k=0;k<3;k++) if(f.v[k]!=kill){ if(a<0) a=f.v[k]; else b=f.v[k]; }
        if(a==b||a<0||b<0) return true; if(a>b) swap(a,b); if(binary_search(have.begin(),have.end(),make_pair(a,b))) return true;
    }
    return false;
}
static bool topologyOK(int u,int v){
    if(u==v||!aliveV[u]||!aliveV[v]) return false; int ids[3]; if(edgeFaces(u,v,ids)!=2) return false;
    static thread_local vector<int> nu,nv; collectNeighbors(u,nu); collectNeighbors(v,nv);
    int common=0; size_t i=0,j=0; while(i<nu.size()&&j<nv.size()){
        if(nu[i]==nv[j]){ ++common; ++i; ++j; if(common>2) return false; }
        else if(nu[i]<nv[j]) ++i; else ++j;
    }
    if(common!=2) return false;
    if(duplicateAfterCollapse(u,v)) return false;
    if(duplicateAfterCollapse(v,u)) return false;
    return true;
}
static bool localGeometryOK(int u,int v,const Vec3& p){
    auto scan=[&](int x)->bool{
        for(int id:inc[x]){
            if(!F[id].alive) continue; const Face& f=F[id]; bool hu=faceHas(f,u), hv=faceHas(f,v); if(!hu&&!hv) continue; if(hu&&hv) continue;
            Vec3 a=(f.v[0]==u||f.v[0]==v)?p:P[f.v[0]], b=(f.v[1]==u||f.v[1]==v)?p:P[f.v[1]], c=(f.v[2]==u||f.v[2]==v)?p:P[f.v[2]];
            Vec3 nr=cross3(b-a,c-a); if(!(norm2(nr)>1e-28)) return false; Vec3 old=faceNormalRaw(f); if(dot3(nr,old)<=1e-24) return false;
        }
        return true;
    };
    return scan(u)&&scan(v);
}

struct CandidateInfo { bool ok=false; Vec3 p{}; Box3 b{}; Quadric q{}; double qw=0.0; double cost=0.0; };
static bool solveOptimal(const Quadric& q,Vec3& p){
    double a00=q.q[0],a01=q.q[1],a02=q.q[2],a10=q.q[1],a11=q.q[4],a12=q.q[5],a20=q.q[2],a21=q.q[5],a22=q.q[7];
    double b0=-q.q[3],b1=-q.q[6],b2=-q.q[8];
    double det=a00*(a11*a22-a12*a21)-a01*(a10*a22-a12*a20)+a02*(a10*a21-a11*a20);
    if(fabs(det)<1e-14) return false; double inv=1.0/det;
    double d0=b0*(a11*a22-a12*a21)-a01*(b1*a22-a12*b2)+a02*(b1*a21-a11*b2);
    double d1=a00*(b1*a22-a12*b2)-b0*(a10*a22-a12*a20)+a02*(a10*b2-b1*a20);
    double d2=a00*(a11*b2-b1*a21)-a01*(a10*b2-b1*a20)+b0*(a10*a21-a11*a20);
    p={d0*inv,d1*inv,d2*inv}; return isfinite(p.x)&&isfinite(p.y)&&isfinite(p.z);
}
static double edgeFeaturePenalty(int u,int v){
    int ids[3]; if(edgeFaces(u,v,ids)!=2) return 10.0; Vec3 n1,n2; if(!faceUnitNormal(F[ids[0]],n1)||!faceUnitNormal(F[ids[1]],n2)) return 10.0;
    double d=max(-1.0,min(1.0,dot3(n1,n2))); return max(0.0,1.0-d);
}
static CandidateInfo makeCandidate(int u,int v,bool checkTopo){
    CandidateInfo ci; if(u==v||!aliveV[u]||!aliveV[v]) return ci; if(checkTopo && !topologyOK(u,v)) return ci;
    ci.b=mergeBox(bboxV[u],bboxV[v]); ci.q=qadd(Q[u],Q[v]); ci.qw=Qw[u]+Qw[v];
    Vec3 cand[8]; int cn=0; cand[cn++]=P[u]; cand[cn++]=P[v]; cand[cn++]=(P[u]+P[v])*0.5;
    cand[cn++]=(P[u]*(double)csz[u]+P[v]*(double)csz[v])/(double)(csz[u]+csz[v]);
    Vec3 opt; if(solveOptimal(ci.q,opt)) cand[cn++]=opt;
    cand[cn++]={(ci.b.mn[0]+ci.b.mx[0])*0.5,(ci.b.mn[1]+ci.b.mx[1])*0.5,(ci.b.mn[2]+ci.b.mx[2])*0.5};
    double best=numeric_limits<double>::infinity(); Vec3 bestp=cand[0]; double len2=norm2(P[u]-P[v]); double feature=checkTopo?edgeFeaturePenalty(u,v):0.0;
    for(int i=0;i<cn;i++){
        Vec3 p=cand[i]; if(!isfinite(p.x)||!isfinite(p.y)||!isfinite(p.z)) continue; if(!coversBox(p,ci.b)) continue; if(checkTopo && !localGeometryOK(u,v,p)) continue;
        double qe=max(0.0,qeval(ci.q,p)); double rms=sqrt(qe/max(1e-18,ci.qw));
        double c=rms*rms + 2.5e-4*len2 + feature*(0.12*len2+0.015*hausTol2) + 2.0e-5*boxMaxDist2(p,ci.b);
        if(c<best){ best=c; bestp=p; }
    }
    if(!isfinite(best)) return ci;
    double rms=sqrt(max(0.0,qeval(ci.q,bestp))/max(1e-18,ci.qw)); double limit=hausTol*(N0>=300000?0.62:(N0>=40000?0.55:0.48));
    if(checkTopo && rms>limit) return ci;
    ci.ok=true; ci.p=bestp; ci.cost=best; return ci;
}

struct PQItem { double cost; int u,v,vu,vv; unsigned char precise; bool operator<(const PQItem& o) const { return cost>o.cost; } };
static void pushEdge(priority_queue<PQItem>& pq,int u,int v,bool precise=false){ if(u==v||!aliveV[u]||!aliveV[v]) return; CandidateInfo ci=makeCandidate(u,v,false); if(!ci.ok) return; pq.push({ci.cost,u,v,verV[u],verV[v],(unsigned char)(precise?1:0)}); }
static int chooseKeep(int u,int v){
    size_t iu=inc[u].size(), iv=inc[v].size(); if(iu>iv*11/10) return u; if(iv>iu*11/10) return v; if(csz[u]!=csz[v]) return csz[u]>csz[v]?u:v; return min(u,v);
}
static int performCollapseCore(int u,int v,const CandidateInfo& ci){
    int keep=chooseKeep(u,v), kill=(keep==u?v:u);
    P[keep]=ci.p; bboxV[keep]=ci.b; Q[keep]=ci.q; Qw[keep]=ci.qw; csz[keep]+=csz[kill];
    aliveV[kill]=0; --activeV; ++verV[keep]; ++verV[kill];
    vector<int> killedInc; killedInc.swap(inc[kill]);
    for(int id:killedInc){
        if(!F[id].alive) continue; Face& f=F[id]; bool hasK=faceHas(f,keep), hasD=faceHas(f,kill); if(!hasD) continue;
        if(hasK){ f.alive=0; --activeF; continue; }
        for(int k=0;k<3;k++) if(f.v[k]==kill) f.v[k]=keep;
        if(f.v[0]==f.v[1]||f.v[1]==f.v[2]||f.v[2]==f.v[0]){ f.alive=0; --activeF; }
        else inc[keep].push_back(id);
    }
    if(inc[keep].size()>768) compactInc(keep); return keep;
}
static void performCollapse(int u,int v,const CandidateInfo& ci,priority_queue<PQItem>& pq){ int keep=performCollapseCore(u,v,ci); static thread_local vector<int> nb; collectNeighbors(keep,nb); for(int w:nb) pushEdge(pq,keep,w,false); }

static double desiredRatio(){
    if(N0<=20) return 0.88;
    if(N0<=6000) return 0.115;
    if(N0<=30000) return 0.090;
    if(N0<=70000) return 0.075;
    if(N0<=500000) return 0.060;
    return 0.050;
}
static int targetVertexCount(){ int target=max(4,(int)ceil(N0*desiredRatio())); if(N0>100000) target=max(target,18000); if(N0>800000) target=max(target,42000); return target; }
static int greedyLargePass(int target){
    vector<uint64_t> keys; keys.reserve((size_t)activeF*3);
    for(const Face& f:F) if(f.alive){ int a=f.v[0],b=f.v[1],c=f.v[2]; if(a==b||b==c||c==a) continue; keys.push_back(edgeKey(a,b)); keys.push_back(edgeKey(b,c)); keys.push_back(edgeKey(c,a)); }
    sort(keys.begin(),keys.end()); keys.erase(unique(keys.begin(),keys.end()),keys.end());
    struct EdgeOrd{ float len2; uint64_t key; }; vector<EdgeOrd> ord; ord.reserve(keys.size());
    for(uint64_t k:keys){ int u=(int)(k>>32), v=(int)(k&0xffffffffu); if((unsigned)u<(unsigned)N0&&(unsigned)v<(unsigned)N0&&aliveV[u]&&aliveV[v]) ord.push_back({(float)norm2(P[u]-P[v]),k}); }
    vector<uint64_t>().swap(keys); sort(ord.begin(),ord.end(),[](const EdgeOrd& a,const EdgeOrd& b){ if(a.len2!=b.len2) return a.len2<b.len2; return a.key<b.key; });
    int collapsed=0; for(const EdgeOrd& e:ord){
        if(activeV<=target) break; if((collapsed&4095)==0){ double elapsed=chrono::duration<double>(chrono::steady_clock::now()-startTime).count(); if(elapsed>18.3) break; }
        int u=(int)(e.key>>32), v=(int)(e.key&0xffffffffu); if(!aliveV[u]||!aliveV[v]) continue; if(!currentEdge(u,v)) continue; CandidateInfo ci=makeCandidate(u,v,true); if(!ci.ok) continue; performCollapseCore(u,v,ci); ++collapsed;
    }
    return collapsed;
}
static bool decimateLargeGreedy(int target){
    if(N0<160000) return false; int lastV=activeV;
    for(int pass=0; pass<16 && activeV>target; ++pass){ double elapsed=chrono::duration<double>(chrono::steady_clock::now()-startTime).count(); if(elapsed>18.3) break; int c=greedyLargePass(target); if(c<512) break; if(activeV>=lastV-max(512,lastV/200)) break; lastV=activeV; }
    return true;
}
static void decimate(){
    if(N0<=4) return; initializeQuadricsAndTolerance();
    vector<uint64_t> keys; keys.reserve((size_t)M0*3);
    for(int i=0;i<M0;i++){ int a=F[i].v[0],b=F[i].v[1],c=F[i].v[2]; if((unsigned)a<(unsigned)N0&&(unsigned)b<(unsigned)N0&&(unsigned)c<(unsigned)N0){ keys.push_back(edgeKey(a,b)); keys.push_back(edgeKey(b,c)); keys.push_back(edgeKey(c,a)); } }
    sort(keys.begin(),keys.end()); keys.erase(unique(keys.begin(),keys.end()),keys.end());
    priority_queue<PQItem> pq; for(uint64_t k:keys){ int u=(int)(k>>32), v=(int)(k&0xffffffffu); pushEdge(pq,u,v,false); } vector<uint64_t>().swap(keys);
    int target=targetVertexCount(); if(decimateLargeGreedy(target)) return;
    int refined=0, collapsed=0;
    while(activeV>target && !pq.empty()){
        if((collapsed&4095)==0){ double elapsed=chrono::duration<double>(chrono::steady_clock::now()-startTime).count(); if(elapsed>18.6) break; }
        PQItem it=pq.top(); pq.pop(); int u=it.u,v=it.v;
        if((unsigned)u>=(unsigned)N0||(unsigned)v>=(unsigned)N0) continue; if(!aliveV[u]||!aliveV[v]) continue;
        if(it.vu!=verV[u]||it.vv!=verV[v]){ if(currentEdge(u,v)) pushEdge(pq,u,v,false); continue; }
        if(!currentEdge(u,v)) continue; CandidateInfo ci=makeCandidate(u,v,true); if(!ci.ok) continue;
        if(!it.precise){ pq.push({ci.cost,u,v,verV[u],verV[v],1}); if(++refined>6000000 && activeV>target*12/10) refined=0; continue; }
        performCollapse(u,v,ci,pq); ++collapsed;
    }
}

static bool trySampleLikeSimplify(){ return false; }
static bool buildCompact(vector<Vec3>& outP, vector<array<int,3>>& outF){
    vector<int> mp(N0,-1); outF.clear(); outF.reserve(max(0,activeF));
    for(const Face& f:F) if(f.alive){
        int a=f.v[0],b=f.v[1],c=f.v[2]; if((unsigned)a>=(unsigned)N0||(unsigned)b>=(unsigned)N0||(unsigned)c>=(unsigned)N0) continue; if(a==b||b==c||c==a) continue;
        Vec3 cr=cross3(P[b]-P[a],P[c]-P[a]); if(norm2(cr)<=1e-30) continue;
        if(mp[a]<0) mp[a]=0; if(mp[b]<0) mp[b]=0; if(mp[c]<0) mp[c]=0; outF.push_back({a,b,c});
    }
    int cnt=0; for(int i=0;i<N0;i++) if(mp[i]==0) mp[i]=cnt++;
    outP.resize(cnt); for(int i=0;i<N0;i++) if(mp[i]>=0) outP[mp[i]]=P[i];
    for(auto& t:outF){ t[0]=mp[t[0]]; t[1]=mp[t[1]]; t[2]=mp[t[2]]; }
    return !outP.empty() && !outF.empty();
}
static bool validateMesh(const vector<Vec3>& V,const vector<array<int,3>>& G){
    if(V.empty()||V.size()>(size_t)N0) return false; vector<uint64_t> edges; edges.reserve(G.size()*3);
    for(auto t:G){ int a=t[0],b=t[1],c=t[2]; if((unsigned)a>=V.size()||(unsigned)b>=V.size()||(unsigned)c>=V.size()) return false; if(a==b||b==c||c==a) return false; if(norm2(cross3(V[b]-V[a],V[c]-V[a]))<=1e-28) return false; edges.push_back(edgeKey(a,b)); edges.push_back(edgeKey(b,c)); edges.push_back(edgeKey(c,a)); }
    sort(edges.begin(),edges.end()); for(size_t i=0;i<edges.size();){ size_t j=i+1; while(j<edges.size()&&edges[j]==edges[i]) ++j; if(j-i!=2) return false; i=j; } return true;
}
static void makeOriginalOutput(vector<Vec3>& outP, vector<array<int,3>>& outF){ outP=P_orig; outF.clear(); outF.reserve(F_orig.size()); for(const Face& f:F_orig) outF.push_back({f.v[0],f.v[1],f.v[2]}); }
static void writeMesh(const vector<Vec3>& V,const vector<array<int,3>>& G){
    printf("%zu %zu\n",V.size(),G.size());
    for(const Vec3& p:V) printf("v %.10g %.10g %.10g\n",p.x,p.y,p.z);
    for(auto t:G) printf("f %d %d %d\n",t[0]+1,t[1]+1,t[2]+1);
}
int main(){
    startTime=chrono::steady_clock::now(); loadInput();
    if(!trySampleLikeSimplify()) decimate();
    vector<Vec3> outP; vector<array<int,3>> outF; bool ok=buildCompact(outP,outF); if(ok) ok=validateMesh(outP,outF); if(!ok) makeOriginalOutput(outP,outF); writeMesh(outP,outF); return 0;
}
#undef main
}
namespace B{
using namespace std;
#define main main_b

struct Vec3{double x,y,z;};
struct Face{int a,b,c;};
static inline Vec3 operator+(Vec3 a,Vec3 b){return {a.x+b.x,a.y+b.y,a.z+b.z};}
static inline Vec3 operator-(Vec3 a,Vec3 b){return {a.x-b.x,a.y-b.y,a.z-b.z};}
static inline Vec3 operator*(Vec3 a,double s){return {a.x*s,a.y*s,a.z*s};}
static inline Vec3 operator/(Vec3 a,double s){return {a.x/s,a.y/s,a.z/s};}
static inline double dotv(Vec3 a,Vec3 b){return a.x*b.x+a.y*b.y+a.z*b.z;}
static inline Vec3 crossv(Vec3 a,Vec3 b){return {a.y*b.z-a.z*b.y,a.z*b.x-a.x*b.z,a.x*b.y-a.y*b.x};}
static inline double norm2(Vec3 a){return dotv(a,a);} 
static inline double dist2(Vec3 a,Vec3 b){return norm2(a-b);} 
static inline double normv(Vec3 a){return sqrt(norm2(a));}
static inline Vec3 unitv(Vec3 a){double n=normv(a);return n>1e-300?a/n:Vec3{0,0,1};}
static inline uint64_t ekey(int a,int b){if(a>b)swap(a,b);return (uint64_t)(uint32_t)a<<32 | (uint32_t)b;}
static inline int eka(uint64_t k){return (int)(k>>32);} static inline int ekb(uint64_t k){return (int)(k&0xffffffffu);} 

struct Inp{vector<char>b;char*p;Inp(){char tmp[1<<16];size_t n;while((n=fread(tmp,1,sizeof(tmp),stdin))>0)b.insert(b.end(),tmp,tmp+n);b.push_back(0);p=b.data();}inline void ws(){while(*p==' '||*p=='\n'||*p=='\r'||*p=='\t')++p;}int ni(){ws();int s=1;if(*p=='-')s=-1,++p;int x=0;while(*p>='0'&&*p<='9')x=x*10+*p++-'0';return s*x;}double nd(){ws();char*e;double x=strtod(p,&e);p=e;return x;}char nc(){ws();return *p++;}};

struct Quadric{double q[10];Quadric(){memset(q,0,sizeof(q));}void addp(double a,double b,double c,double d,double w){q[0]+=w*a*a;q[1]+=w*a*b;q[2]+=w*a*c;q[3]+=w*a*d;q[4]+=w*b*b;q[5]+=w*b*c;q[6]+=w*b*d;q[7]+=w*c*c;q[8]+=w*c*d;q[9]+=w*d*d;}void add(const Quadric&o){for(int i=0;i<10;i++)q[i]+=o.q[i];}double eval(Vec3 p)const{double x=p.x,y=p.y,z=p.z;return q[0]*x*x+2*q[1]*x*y+2*q[2]*x*z+2*q[3]*x+q[4]*y*y+2*q[5]*y*z+2*q[6]*y+q[7]*z*z+2*q[8]*z+q[9];}};

struct Mesh{vector<Vec3>V;vector<Face>F;void swap(Mesh&o){V.swap(o.V);F.swap(o.F);}};
static vector<Vec3> OrigV; static vector<Face> OrigF; static Vec3 BB0,BB1; static double DIAG=1, EPS=0.05; static chrono::steady_clock::time_point T0; static double TL=19.2;
static bool tout(){return chrono::duration<double>(chrono::steady_clock::now()-T0).count()>TL;}

struct TriKey{int a,b,c;bool operator==(TriKey const&o)const{return a==o.a&&b==o.b&&c==o.c;}};
struct TriHash{size_t operator()(TriKey const&t)const{return ((uint64_t)t.a*11995408973635179863ull)^((uint64_t)t.b*10150724397891781847ull)^((uint64_t)t.c*1442695040888963407ull);}};
static TriKey tkey(Face f){int a[3]={f.a,f.b,f.c};sort(a,a+3);return {a[0],a[1],a[2]};}

static void bounds(const vector<Vec3>&V){BB0={1e100,1e100,1e100};BB1={-1e100,-1e100,-1e100};for(auto&p:V){BB0.x=min(BB0.x,p.x);BB0.y=min(BB0.y,p.y);BB0.z=min(BB0.z,p.z);BB1.x=max(BB1.x,p.x);BB1.y=max(BB1.y,p.y);BB1.z=max(BB1.z,p.z);}DIAG=normv(BB1-BB0);if(!(DIAG>0))DIAG=1;EPS=0.05*DIAG*0.999;}

static bool valid_closed(const Mesh&m){int n=m.V.size();if(n<4||m.F.size()<4)return false;vector<uint64_t>ed;ed.reserve(m.F.size()*3);unordered_set<TriKey,TriHash>ts;ts.reserve(m.F.size()*2+16);for(auto f:m.F){if(f.a<0||f.a>=n||f.b<0||f.b>=n||f.c<0||f.c>=n)return false;if(f.a==f.b||f.a==f.c||f.b==f.c)return false;Vec3 cr=crossv(m.V[f.b]-m.V[f.a],m.V[f.c]-m.V[f.a]);if(norm2(cr)<=1e-28*DIAG*DIAG*DIAG*DIAG)return false;TriKey tk=tkey(f);if(!ts.insert(tk).second)return false;ed.push_back(ekey(f.a,f.b));ed.push_back(ekey(f.b,f.c));ed.push_back(ekey(f.c,f.a));}sort(ed.begin(),ed.end());for(size_t i=0;i<ed.size();){size_t j=i+1;while(j<ed.size()&&ed[j]==ed[i])++j;if(j-i!=2)return false;i=j;}return true;}

struct PGrid{double cell;Vec3 mn;unordered_map<long long,vector<int>>g;static long long key(long long x,long long y,long long z){return (x*73856093ll)^(y*19349663ll)^(z*83492791ll);}array<long long,3> idx(Vec3 p)const{return {(long long)floor((p.x-mn.x)/cell),(long long)floor((p.y-mn.y)/cell),(long long)floor((p.z-mn.z)/cell)};}void build(const vector<Vec3>&V,double c){cell=max(c,1e-12);mn=BB0;g.clear();g.reserve(V.size()*2+16);for(int i=0;i<(int)V.size();++i){auto a=idx(V[i]);g[key(a[0],a[1],a[2])].push_back(i);}}bool near(Vec3 p,const vector<Vec3>&V,double r,int*best=nullptr)const{double r2=r*r;auto a=idx(p);double bd=r2;int bi=-1;for(long long dx=-1;dx<=1;dx++)for(long long dy=-1;dy<=1;dy++)for(long long dz=-1;dz<=1;dz++){auto it=g.find(key(a[0]+dx,a[1]+dy,a[2]+dz));if(it==g.end())continue;for(int id:it->second){double d=dist2(p,V[id]);if(d<=bd){bd=d;bi=id;}}}if(best)*best=bi;return bi>=0;}};

static bool covers_orig(const vector<Vec3>&sel,double r){PGrid pg;pg.build(sel,r);for(auto&p:OrigV)if(!pg.near(p,sel,r))return false;return true;}

static void orient_add(Mesh&out,int a,int b,int c,Vec3 cen){Vec3 cr=crossv(out.V[b]-out.V[a],out.V[c]-out.V[a]);Vec3 m=(out.V[a]+out.V[b]+out.V[c])/3.0;if(dotv(cr,m-cen)<0)swap(b,c);out.F.push_back({a,b,c});}

static bool try_box(Mesh&out){int n=OrigV.size();if(n<200)return false;double dx=BB1.x-BB0.x,dy=BB1.y-BB0.y,dz=BB1.z-BB0.z;double mind=max(1e-12,DIAG);double surf=0;for(auto&p:OrigV){double d=min({fabs(p.x-BB0.x),fabs(p.x-BB1.x),fabs(p.y-BB0.y),fabs(p.y-BB1.y),fabs(p.z-BB0.z),fabs(p.z-BB1.z)});if(d<EPS*0.18)surf++;}if(surf<n*0.985)return false;double step=EPS*0.82;int nx=max(1,(int)ceil(dx/step)),ny=max(1,(int)ceil(dy/step)),nz=max(1,(int)ceil(dz/step));long long nodes=(long long)(nx+1)*(ny+1)*(nz+1);if(nodes>250000||nodes>n*0.85)return false;PGrid og;og.build(OrigV,EPS);unordered_map<long long,int>id;id.reserve((nx+1)*(ny+1)*2+(nx+1)*(nz+1)*2+(ny+1)*(nz+1)*2);vector<int>origId;Mesh m;auto key=[&](int i,int j,int k){return ((long long)i<<42)^((long long)j<<21)^k;};auto addnode=[&](int i,int j,int k)->int{long long kk=key(i,j,k);auto it=id.find(kk);if(it!=id.end())return it->second;Vec3 q{BB0.x+dx*(double)i/nx,BB0.y+dy*(double)j/ny,BB0.z+dz*(double)k/nz};int bi=-1;if(!og.near(q,OrigV,EPS*0.96,&bi))return -1;int t=m.V.size();id[kk]=t;m.V.push_back(OrigV[bi]);origId.push_back(bi);return t;};auto quad=[&](int a,int b,int c,int d,Vec3 want){if(a<0||b<0||c<0||d<0)return false;if(a==b||a==c||a==d||b==c||b==d||c==d)return false;Vec3 cr=crossv(m.V[b]-m.V[a],m.V[c]-m.V[a]);if(dotv(cr,want)<0){swap(b,d);}m.F.push_back({a,b,c});m.F.push_back({a,c,d});return true;};for(int i=0;i<nx;i++)for(int j=0;j<ny;j++){if(!quad(addnode(i,j,0),addnode(i+1,j,0),addnode(i+1,j+1,0),addnode(i,j+1,0),{0,0,-1}))return false;if(!quad(addnode(i,j,nz),addnode(i,j+1,nz),addnode(i+1,j+1,nz),addnode(i+1,j,nz),{0,0,1}))return false;}for(int i=0;i<nx;i++)for(int k=0;k<nz;k++){if(!quad(addnode(i,0,k),addnode(i,0,k+1),addnode(i+1,0,k+1),addnode(i+1,0,k),{0,-1,0}))return false;if(!quad(addnode(i,ny,k),addnode(i+1,ny,k),addnode(i+1,ny,k+1),addnode(i,ny,k+1),{0,1,0}))return false;}for(int j=0;j<ny;j++)for(int k=0;k<nz;k++){if(!quad(addnode(0,j,k),addnode(0,j+1,k),addnode(0,j+1,k+1),addnode(0,j,k+1),{-1,0,0}))return false;if(!quad(addnode(nx,j,k),addnode(nx,j,k+1),addnode(nx,j+1,k+1),addnode(nx,j+1,k),{1,0,0}))return false;}if((int)m.V.size()>=n||!valid_closed(m))return false;if(!covers_orig(m.V,EPS*0.985))return false;out.swap(m);return true;}

static int genus0_possible(){Mesh om{OrigV,OrigF};vector<uint64_t>ed;ed.reserve(OrigF.size()*3);vector<int>used(OrigV.size());for(auto f:OrigF){used[f.a]=used[f.b]=used[f.c]=1;ed.push_back(ekey(f.a,f.b));ed.push_back(ekey(f.b,f.c));ed.push_back(ekey(f.c,f.a));}sort(ed.begin(),ed.end());int E=0;for(size_t i=0;i<ed.size();){size_t j=i+1;while(j<ed.size()&&ed[j]==ed[i])++j;if(j-i!=2)return 0;E++;i=j;}int V=accumulate(used.begin(),used.end(),0),F=OrigF.size();int chi=V-E+F;return chi==2;}

static bool try_star(Mesh&out){int n=OrigV.size();if(n<3000||!genus0_possible())return false;Vec3 cen{0,0,0};for(auto&p:OrigV)cen=cen+p;cen=cen/(double)n;vector<double>rs;rs.reserve(n);for(auto&p:OrigV){double r=normv(p-cen);if(r>1e-9)rs.push_back(r);}if(rs.size()<n*0.98)return false;nth_element(rs.begin(),rs.begin()+rs.size()*9/10,rs.end());double r90=rs[rs.size()*9/10];if(!(r90>0))return false;double step=min(0.35,max(0.045,EPS/(r90*0.72)));int lon=max(12,min(192,(int)ceil(2*M_PI/step)));int lat=max(6,min(96,(int)ceil(M_PI/step)));if((long long)(lat-1)*lon+2>=n*0.72)return false;auto dir=[&](int i,int j){double ph=M_PI*i/lat,th=2*M_PI*j/lon;return Vec3{sin(ph)*cos(th),sin(ph)*sin(th),cos(ph)};};vector<int>best((lat-1)*lon+2,-1);vector<double>bs(best.size(),-2);auto slot=[&](int i,int j){return 1+(i-1)*lon+(j%lon+lon)%lon;};for(int idx=0;idx<n;idx++){Vec3 q=OrigV[idx]-cen;double r=normv(q);if(r<=1e-12)continue;q=q/r;double ph=acos(max(-1.0,min(1.0,q.z)));double th=atan2(q.y,q.x);if(th<0)th+=2*M_PI;if(ph<M_PI/(2*lat)){double sc=q.z;if(sc>bs[0])bs[0]=sc,best[0]=idx;continue;}if(M_PI-ph<M_PI/(2*lat)){double sc=-q.z;int s=best.size()-1;if(sc>bs[s])bs[s]=sc,best[s]=idx;continue;}int ii=min(lat-1,max(1,(int)floor(ph/M_PI*lat+0.5)));int jj=((int)floor(th/(2*M_PI)*lon+0.5))%lon;for(int di=-1;di<=1;di++)for(int dj=-1;dj<=1;dj++){int ni=ii+di;if(ni<1||ni>=lat)continue;int s=slot(ni,jj+dj);double sc=dotv(q,dir(ni,jj+dj));if(sc>bs[s])bs[s]=sc,best[s]=idx;}}for(int x:best)if(x<0)return false;unordered_set<int>uniq(best.begin(),best.end());if(uniq.size()!=best.size())return false;Mesh m;for(int id:best)m.V.push_back(OrigV[id]);Vec3 zero{0,0,0};for(int j=0;j<lon;j++){orient_add(m,0,slot(1,j),slot(1,j+1),cen);}for(int i=1;i<lat-1;i++)for(int j=0;j<lon;j++){int a=slot(i,j),b=slot(i,j+1),c=slot(i+1,j+1),d=slot(i+1,j);orient_add(m,a,b,c,cen);orient_add(m,a,c,d,cen);}int south=best.size()-1;for(int j=0;j<lon;j++){orient_add(m,south,slot(lat-1,j+1),slot(lat-1,j),cen);}if(!valid_closed(m))return false;if(!covers_orig(m.V,EPS*0.985))return false;if((int)m.V.size()>=n)return false;out.swap(m);return true;}

struct EdgeRec{uint64_t k;int opp;bool operator<(EdgeRec const&o)const{return k<o.k||(k==o.k&&opp<o.opp);}};
struct Cand{float cost;int u,v,o1,o2;bool operator<(Cand const&o)const{return cost<o.cost;}};

static bool common_ok(const vector<int>&a,const vector<int>&b,int o1,int o2){int i=0,j=0,c=0;bool h1=0,h2=0;while(i<(int)a.size()&&j<(int)b.size()){if(a[i]==b[j]){int x=a[i];if(x!=o1&&x!=o2)return false;c++;if(x==o1)h1=1;if(x==o2)h2=1;i++;j++;}else if(a[i]<b[j])i++;else j++;}return c==2&&h1&&h2;}

static bool batch_simplify(Mesh&mesh){int n=mesh.V.size();if(n<5)return false;vector<double>rad(n,0);int maxPass=n>400000?4:(n>120000?5:7);bool anyAll=false;for(int pass=0;pass<maxPass && !tout();++pass){n=mesh.V.size();int mf=mesh.F.size();vector<Vec3>fn(mf);vector<double>fa(mf);vector<vector<int>>inc(n),nei(n);vector<EdgeRec>er;er.reserve((size_t)mf*3);vector<Quadric>Q(n);vector<int>deg(n);for(auto f:mesh.F){deg[f.a]++;deg[f.b]++;deg[f.c]++;}for(int i=0;i<n;i++){inc[i].reserve(deg[i]);nei[i].reserve(min(32,deg[i]*2+4));}for(int i=0;i<mf;i++){auto f=mesh.F[i];Vec3 cr=crossv(mesh.V[f.b]-mesh.V[f.a],mesh.V[f.c]-mesh.V[f.a]);double l=normv(cr);if(!(l>1e-30))continue;fn[i]=cr/l;fa[i]=0.5*l;double d=-dotv(fn[i],mesh.V[f.a]);double w=max(1e-12,fa[i]);Q[f.a].addp(fn[i].x,fn[i].y,fn[i].z,d,w);Q[f.b].addp(fn[i].x,fn[i].y,fn[i].z,d,w);Q[f.c].addp(fn[i].x,fn[i].y,fn[i].z,d,w);inc[f.a].push_back(i);inc[f.b].push_back(i);inc[f.c].push_back(i);nei[f.a].push_back(f.b);nei[f.a].push_back(f.c);nei[f.b].push_back(f.a);nei[f.b].push_back(f.c);nei[f.c].push_back(f.a);nei[f.c].push_back(f.b);er.push_back({ekey(f.a,f.b),f.c});er.push_back({ekey(f.b,f.c),f.a});er.push_back({ekey(f.c,f.a),f.b});}for(auto&v:nei){sort(v.begin(),v.end());v.erase(unique(v.begin(),v.end()),v.end());}sort(er.begin(),er.end());vector<Cand>cs;cs.reserve(er.size()/3);double lim=EPS*(pass==0?0.88:0.96);double lim2=lim*lim;for(size_t i=0;i<er.size();){size_t j=i+1;while(j<er.size()&&er[j].k==er[i].k)++j;if(j==i+2){int u=eka(er[i].k),v=ekb(er[i].k);double l2=dist2(mesh.V[u],mesh.V[v]);if(l2<=lim2*1.08 && min(rad[u]+sqrt(l2),rad[v]+sqrt(l2))<=lim){double nd=dotv(fn[er[i].opp<0?0:0],fn[er[i].opp<0?0:0]);(void)nd;float cost=(float)(l2/(lim2+1e-30)+0.05*(nei[u].size()+nei[v].size()));cs.push_back({cost,u,v,er[i].opp,er[i+1].opp});}}i=j;}if(cs.empty())break;sort(cs.begin(),cs.end());vector<char>used(n,0),dead(n,0);vector<int>par(n);iota(par.begin(),par.end(),0);vector<double>nrad=rad;int accepted=0;auto normal_ok=[&](int rem,int keep)->bool{double ndot=pass==0?0.999999:(pass<3?0.965:0.92);for(int fid:inc[rem]){Face f=mesh.F[fid];bool hk=(f.a==keep||f.b==keep||f.c==keep);if(hk)continue;Vec3 a=mesh.V[f.a],b=mesh.V[f.b],c=mesh.V[f.c];if(f.a==rem)a=mesh.V[keep];if(f.b==rem)b=mesh.V[keep];if(f.c==rem)c=mesh.V[keep];Vec3 cr=crossv(b-a,c-a);double l=normv(cr);if(!(l>1e-28*DIAG*DIAG))return false;double d=dotv(cr/l,fn[fid]);if(d<ndot)return false;}return true;};for(const auto&c:cs){if(tout())break;int u=c.u,v=c.v;if(used[u]||used[v])continue;if(!common_ok(nei[u],nei[v],c.o1,c.o2))continue;double duv=normv(mesh.V[u]-mesh.V[v]);int keep=-1,rem=-1;double ru=max(rad[u],rad[v]+duv),rv=max(rad[v],rad[u]+duv);if(ru<=lim&&(rv>lim||Q[u].eval(mesh.V[u])<=Q[v].eval(mesh.V[v])))keep=u,rem=v;else if(rv<=lim)keep=v,rem=u;else continue;if(!normal_ok(rem,keep))continue;used[u]=used[v]=1;dead[rem]=1;par[rem]=keep;nrad[keep]=max(nrad[keep],max(rad[keep],rad[rem]+duv));accepted++;}if(accepted==0)break;Mesh nm;vector<int>id(n,-1);for(int i=0;i<n;i++)if(!dead[i]){id[i]=nm.V.size();nm.V.push_back(mesh.V[i]);}vector<double>rr(nm.V.size());for(int i=0;i<n;i++)if(id[i]>=0)rr[id[i]]=nrad[i];unordered_set<TriKey,TriHash>seen;seen.reserve(mesh.F.size()*2+16);for(auto f:mesh.F){int a=par[f.a],b=par[f.b],c=par[f.c];if(a==b||a==c||b==c)continue;Face g{id[a],id[b],id[c]};if(g.a<0||g.b<0||g.c<0)continue;Vec3 cr=crossv(nm.V[g.b]-nm.V[g.a],nm.V[g.c]-nm.V[g.a]);if(norm2(cr)<=1e-28*DIAG*DIAG*DIAG*DIAG)continue;TriKey tk=tkey(g);if(seen.insert(tk).second)nm.F.push_back(g);}if((int)nm.V.size()>=n||!valid_closed(nm))break;mesh.swap(nm);rad.swap(rr);anyAll=true;if(accepted<max(4,n/500))break;}return anyAll;}


static bool try_ordered_torus(Mesh&out){
    int N=OrigV.size(), M=OrigF.size();
    if(N<400 || M<700) return false;
    unordered_set<TriKey,TriHash> fs; fs.reserve((size_t)M*2+16);
    for(auto f:OrigF) fs.insert(tkey(f));
    auto hastri=[&](int a,int b,int c){return fs.find(tkey({a,b,c}))!=fs.end();};
    Mesh best; int bestn=N+1;
    vector<int> divs;
    for(int c=8;c*c<=N;c++) if(N%c==0){divs.push_back(c); if(c*c!=N) divs.push_back(N/c);} 
    sort(divs.begin(),divs.end());
    for(int C:divs){
        int R=N/C; if(R<8||C<8||R>20000||C>20000) continue;
        long long cells=1LL*R*C; if(cells*2 > (long long)M*13/10 || cells*2 < (long long)M*55/100) continue;
        int probes=0, ok=0; int step=max(1,(int)(cells/6000));
        for(long long t=0;t<cells;t+=step){int i=t/C,j=t%C;int a=i*C+j,b=((i+1)%R)*C+j,c=((i+1)%R)*C+(j+1)%C,d=i*C+(j+1)%C;probes++; if((hastri(a,b,c)&&hastri(a,c,d))||(hastri(a,c,b)&&hastri(a,d,c))||(hastri(a,b,d)&&hastri(b,c,d))||(hastri(a,d,b)&&hastri(b,d,c))) ok++;}
        if(probes<20 || ok*100<probes*70) continue;
        int maxSr=min(R/4,16), maxSc=min(C/4,16); if(maxSr<1||maxSc<1) continue;
        for(int sr=1;sr<=maxSr;sr++)for(int sc=1;sc<=maxSc;sc++){
            int R2=(R+sr-1)/sr, C2=(C+sc-1)/sc; if(R2<4||C2<4) continue; if(R2*C2>=bestn||R2*C2>=N) continue;
            Mesh m; m.V.reserve(R2*C2); vector<int> id(R2*C2);
            for(int i=0;i<R2;i++){int oi=(int)((long long)i*R/R2); for(int j=0;j<C2;j++){int oj=(int)((long long)j*C/C2); id[i*C2+j]=m.V.size(); m.V.push_back(OrigV[oi*C+oj]);}}
            auto ID=[&](int i,int j){i=(i%R2+R2)%R2;j=(j%C2+C2)%C2;return i*C2+j;};
            Vec3 cen{0,0,0}; for(auto&p:m.V) cen=cen+p; cen=cen/(double)m.V.size();
            for(int i=0;i<R2;i++)for(int j=0;j<C2;j++){int a=ID(i,j),b=ID(i+1,j),c=ID(i+1,j+1),d=ID(i,j+1);orient_add(m,a,b,c,cen);orient_add(m,a,c,d,cen);} 
            if((int)m.V.size()<bestn && valid_closed(m) && covers_orig(m.V,EPS*0.985)){best.swap(m);bestn=best.V.size();}
        }
    }
    if(bestn<N){out.swap(best);return true;}return false;
}

static bool try_ordered_uvsphere(Mesh&out){
    int N=OrigV.size(), M=OrigF.size(); if(N<500) return false;
    unordered_set<TriKey,TriHash> fs; fs.reserve((size_t)M*2+16); for(auto f:OrigF) fs.insert(tkey(f));
    auto hastri=[&](int a,int b,int c){return fs.find(tkey({a,b,c}))!=fs.end();};
    Mesh best; int bestn=N+1;
    for(int V=8;V<=min(1000,N/3);V++) if((N-2)%V==0){
        int R=(N-2)/V+1; if(R<4||R>20000) continue; long long cells=1LL*(R-1)*V; if(cells*2 > (long long)M*15/10 || cells*2 < (long long)M*40/100) continue;
        int probes=0,ok=0; int step=max(1,(int)(cells/5000));
        auto rid=[&](int r,int j){return 1+(r-1)*V+(j%V+V)%V;}; int bot=N-1;
        for(long long t=0;t<cells;t+=step){int r=t/V+1,j=t%V;probes++; if(r==1){if(hastri(0,rid(1,j),rid(1,j+1))||hastri(0,rid(1,j+1),rid(1,j)))ok++;}else if(r==R-1){if(hastri(bot,rid(R-1,j+1),rid(R-1,j))||hastri(bot,rid(R-1,j),rid(R-1,j+1)))ok++;}else{int a=rid(r-1,j),b=rid(r-1,j+1),c=rid(r,j+1),d=rid(r,j); if((hastri(a,b,c)&&hastri(a,c,d))||(hastri(a,c,b)&&hastri(a,d,c))) ok++;}}
        if(probes<20||ok*100<probes*58) continue;
        int maxSr=min((R-1)/3,14), maxSc=min(V/4,16);
        for(int sr=1;sr<=maxSr;sr++)for(int sc=1;sc<=maxSc;sc++){
            int R2=(R+sr-1)/sr, V2=(V+sc-1)/sc; if(R2<3||V2<6) continue; int nv=2+(R2-1)*V2; if(nv>=bestn||nv>=N) continue;
            Mesh m; m.V.reserve(nv); m.V.push_back(OrigV[0]);
            for(int r=1;r<=R2-1;r++){int orr=1+(int)((long long)(r-1)*(R-1)/max(1,R2-1)); if(orr>=R)orr=R-1; for(int j=0;j<V2;j++){int oj=(int)((long long)j*V/V2); m.V.push_back(OrigV[rid(orr,oj)]);}}
            int south=m.V.size(); m.V.push_back(OrigV[bot]); auto ID=[&](int r,int j){return 1+(r-1)*V2+(j%V2+V2)%V2;}; Vec3 cen{0,0,0}; for(auto&p:m.V)cen=cen+p; cen=cen/(double)m.V.size();
            for(int j=0;j<V2;j++)orient_add(m,0,ID(1,j+1),ID(1,j),cen); for(int r=1;r<R2-1;r++)for(int j=0;j<V2;j++){int a=ID(r,j),b=ID(r,j+1),c=ID(r+1,j+1),d=ID(r+1,j);orient_add(m,a,b,c,cen);orient_add(m,a,c,d,cen);} for(int j=0;j<V2;j++)orient_add(m,south,ID(R2-1,j),ID(R2-1,j+1),cen);
            if((int)m.V.size()<bestn&&valid_closed(m)&&covers_orig(m.V,EPS*0.985)){best.swap(m);bestn=best.V.size();}
        }
    }
    if(bestn<N){out.swap(best);return true;} return false;
}

static void write_mesh(const Mesh&m){string out;out.reserve(m.V.size()*42+m.F.size()*26+64);char line[160];out.append(line,snprintf(line,sizeof(line),"%d %d\n",(int)m.V.size(),(int)m.F.size()));bool hp=m.V.size()*2<OrigV.size();for(auto&p:m.V){if(hp)out.append(line,snprintf(line,sizeof(line),"v %.15g %.15g %.15g\n",p.x,p.y,p.z));else out.append(line,snprintf(line,sizeof(line),"v %.10g %.10g %.10g\n",p.x,p.y,p.z));}for(auto f:m.F)out.append(line,snprintf(line,sizeof(line),"f %d %d %d\n",f.a+1,f.b+1,f.c+1));fwrite(out.data(),1,out.size(),stdout);} 

int main(){T0=chrono::steady_clock::now();Inp in;int N=in.ni(),M=in.ni();if(N<=0||M<=0)return 0;OrigV.resize(N);for(int i=0;i<N;i++){in.ws();if(*in.p=='v'||*in.p=='V')++in.p;OrigV[i].x=in.nd();OrigV[i].y=in.nd();OrigV[i].z=in.nd();}OrigF.resize(M);for(int i=0;i<M;i++){in.ws();if(*in.p=='f'||*in.p=='F')++in.p;OrigF[i].a=in.ni()-1;OrigF[i].b=in.ni()-1;OrigF[i].c=in.ni()-1;}bounds(OrigV);Mesh best{OrigV,OrigF},cand;
    auto keep=[&](Mesh&m){if(m.V.size()<best.V.size()&&valid_closed(m)&&covers_orig(m.V,EPS*0.999))best.swap(m);};
    Mesh q=best; if(valid_closed(q)&&batch_simplify(q)) keep(q);
    if(try_box(cand)){keep(cand); Mesh t=cand; if(batch_simplify(t)) keep(t);} cand=Mesh();
    if(try_ordered_torus(cand)){keep(cand); Mesh t=cand; if(batch_simplify(t)&&covers_orig(t.V,EPS*0.999)) keep(t);} cand=Mesh();
    if(try_ordered_uvsphere(cand)){keep(cand); Mesh t=cand; if(batch_simplify(t)&&covers_orig(t.V,EPS*0.999)) keep(t);} cand=Mesh();
    if(try_star(cand)){keep(cand); Mesh t=cand; if(batch_simplify(t)&&covers_orig(t.V,EPS*0.999)) keep(t);} 
    write_mesh(best);return 0;}
#undef main
}
namespace C{
using namespace std;
struct V{double x,y,z;};
struct F{int a,b,c;};
static inline V operator+(V a,V b){return {a.x+b.x,a.y+b.y,a.z+b.z};}
static inline V operator-(V a,V b){return {a.x-b.x,a.y-b.y,a.z-b.z};}
static inline V operator*(V a,double s){return {a.x*s,a.y*s,a.z*s};}
static inline double dotp(V a,V b){return a.x*b.x+a.y*b.y+a.z*b.z;}
static inline V crossp(V a,V b){return {a.y*b.z-a.z*b.y,a.z*b.x-a.x*b.z,a.x*b.y-a.y*b.x};}
static inline double n2(V a){return dotp(a,a);} static inline double d2(V a,V b){return n2(a-b);} static inline double nr(V a){return sqrt(max(0.0,n2(a)));}
static inline bool norm(V &a){double l=nr(a); if(l<=1e-300) return false; a=a*(1.0/l); return true;}
static vector<V>P; static vector<F>Fin; static int N,M; static V mnv,mxv; static double diagv,epsv,eps2v; static chrono::steady_clock::time_point T0;
static inline bool time_left(double margin=0){return chrono::duration<double>(chrono::steady_clock::now()-T0).count()+margin<18.8;}
struct In{vector<char>b;char*p;In(){char q[1<<16];size_t n;while((n=fread(q,1,sizeof(q),stdin)))b.insert(b.end(),q,q+n);b.push_back(0);p=b.data();}void ws(){while(*p==' '||*p=='\n'||*p=='\r'||*p=='\t')++p;}int ni(){ws();int s=1;if(*p=='-')s=-1,++p;int x=0;while(*p>='0'&&*p<='9')x=x*10+*p++-'0';return s*x;}double nd(){ws();char*e;double x=strtod(p,&e);p=e;return x;}};
static void bounds(){mnv={1e100,1e100,1e100};mxv={-1e100,-1e100,-1e100};for(auto&p:P){mnv.x=min(mnv.x,p.x);mnv.y=min(mnv.y,p.y);mnv.z=min(mnv.z,p.z);mxv.x=max(mxv.x,p.x);mxv.y=max(mxv.y,p.y);mxv.z=max(mxv.z,p.z);}diagv=nr(mxv-mnv);if(!(diagv>0))diagv=1;epsv=.05*diagv*.9992;eps2v=epsv*epsv;}
static inline unsigned long long ekey(int a,int b){if(a>b)swap(a,b);return ((unsigned long long)(unsigned)a<<32)|(unsigned)b;}
struct TK{int a,b,c;bool operator==(TK const&o)const{return a==o.a&&b==o.b&&c==o.c;}};
struct TH{size_t operator()(TK const&t)const{return (uint64_t)t.a*11995408973635179863ull ^ (uint64_t)t.b*10150724397891781847ull ^ (uint64_t)t.c*1442695040888963407ull;}};
static TK tk(F f){int a[3]={f.a,f.b,f.c};sort(a,a+3);return {a[0],a[1],a[2]};}
struct Mesh{vector<V>v;vector<F>f;};
static bool valid(const Mesh&m){int n=m.v.size();if(n<4||m.f.size()<4)return false;double mina=max(1e-32,1e-30*diagv*diagv*diagv*diagv);vector<unsigned long long>ed;ed.reserve(m.f.size()*3);unordered_set<TK,TH>ts;ts.reserve(m.f.size()*2+11);for(auto f:m.f){if(f.a<0||f.b<0||f.c<0||f.a>=n||f.b>=n||f.c>=n||f.a==f.b||f.a==f.c||f.b==f.c)return false;if(n2(crossp(m.v[f.b]-m.v[f.a],m.v[f.c]-m.v[f.a]))<=mina)return false;if(!ts.insert(tk(f)).second)return false;ed.push_back(ekey(f.a,f.b));ed.push_back(ekey(f.b,f.c));ed.push_back(ekey(f.c,f.a));}sort(ed.begin(),ed.end());for(size_t i=0;i<ed.size();){size_t j=i+1;while(j<ed.size()&&ed[j]==ed[i])++j;if(j-i!=2)return false;i=j;}return true;}
struct Grid{double c;V o;unordered_map<long long,vector<int>>h;static long long key(long long x,long long y,long long z){return x*73856093ll^y*19349663ll^z*83492791ll;}array<long long,3> id(V p)const{return {(long long)floor((p.x-o.x)/c),(long long)floor((p.y-o.y)/c),(long long)floor((p.z-o.z)/c)};}void init(double cc,size_t cap){c=max(cc,1e-300);o=mnv;h.clear();h.reserve(cap*2+17);}void add(int i,const vector<V>&A){auto a=id(A[i]);h[key(a[0],a[1],a[2])].push_back(i);}bool near(V p,const vector<V>&A,double r2)const{auto a=id(p);for(long long dx=-1;dx<=1;dx++)for(long long dy=-1;dy<=1;dy++)for(long long dz=-1;dz<=1;dz++){auto it=h.find(key(a[0]+dx,a[1]+dy,a[2]+dz));if(it==h.end())continue;for(int j:it->second)if(d2(p,A[j])<=r2)return true;}return false;}bool near_id(V p,const vector<int>&ids,double r2)const{auto a=id(p);for(long long dx=-1;dx<=1;dx++)for(long long dy=-1;dy<=1;dy++)for(long long dz=-1;dz<=1;dz++){auto it=h.find(key(a[0]+dx,a[1]+dy,a[2]+dz));if(it==h.end())continue;for(int k:it->second)if(d2(p,P[ids[k]])<=r2)return true;}return false;}};
static bool cover(const vector<V>&A,const vector<V>&B,double r){Grid g;g.init(r,B.size());for(int i=0;i<(int)B.size();++i)g.add(i,B);double r2=r*r;for(auto&p:A)if(!g.near(p,B,r2))return false;return true;}
static vector<int> greedy(long long step,long long off,int cutoff){vector<int>s;s.reserve(min(N,cutoff));Grid g;g.init(epsv,sqrt((double)N)+16);double r2=eps2v;for(int t=0;t<N;t++){int i=(int)((off+step*(long long)t)%N);if(i<0)i+=N;if(!g.near_id(P[i],s,r2)){if((int)s.size()>=cutoff)return {};s.push_back(i);g.h[Grid::key((long long)floor((P[i].x-mnv.x)/g.c),(long long)floor((P[i].y-mnv.y)/g.c),(long long)floor((P[i].z-mnv.z)/g.c))].push_back((int)s.size()-1);}}return s;}
static vector<int> best_net(){vector<int>best;int cutoff=N;vector<long long>steps={1,N>1?N-1:1,17,97,997,7919,104729,1000003};int tries=0;for(long long st:steps){if(!time_left(.35))break;st%=N;if(st==0)st=1;while(std::gcd((long long)N,st)!=1)++st;vector<int>s=greedy(st,(st==N-1?N-1:0),cutoff);++tries;if(!s.empty()&&(best.empty()||s.size()<best.size())){best.swap(s);cutoff=max(4,(int)best.size());}if(N>300000&&tries>=5)break;if(N>800000&&tries>=4)break;}return best;}
static bool make_bipyramid(const vector<int>&ids,Mesh&out,bool jit=false){if(ids.size()<4)return false;vector<int>cand;auto addc=[&](int x){if(x>=0&&find(cand.begin(),cand.end(),x)==cand.end())cand.push_back(x);};int ix0=ids[0],ix1=ids[0],iy0=ids[0],iy1=ids[0],iz0=ids[0],iz1=ids[0];for(int i:ids){if(P[i].x<P[ix0].x)ix0=i;if(P[i].x>P[ix1].x)ix1=i;if(P[i].y<P[iy0].y)iy0=i;if(P[i].y>P[iy1].y)iy1=i;if(P[i].z<P[iz0].z)iz0=i;if(P[i].z>P[iz1].z)iz1=i;}for(int x:{ix0,ix1,iy0,iy1,iz0,iz1,ids[0],ids[(int)ids.size()/2]})addc(x);int a=ids[0],b=ids.back();double bd=-1;for(int x:cand)for(int y:cand)if(x!=y&&d2(P[x],P[y])>bd){bd=d2(P[x],P[y]);a=x;b=y;}vector<pair<int,int>>pairs={{a,b},{ix0,ix1},{iy0,iy1},{iz0,iz1}};for(auto pr:pairs){a=pr.first;b=pr.second;if(a==b)continue;vector<int>ring;ring.reserve(ids.size()-2);for(int x:ids)if(x!=a&&x!=b)ring.push_back(x);if(ring.size()<3)continue;V axis=P[b]-P[a];if(!norm(axis))continue;V tmp=fabs(axis.x)<.7?V{1,0,0}:V{0,1,0};V u=crossp(axis,tmp);if(!norm(u))continue;V v=crossp(axis,u);V cen{0,0,0};for(int x:ring)cen=cen+P[x];cen=cen*(1.0/ring.size());sort(ring.begin(),ring.end(),[&](int r1,int r2){V p1=P[r1]-cen,p2=P[r2]-cen;double a1=atan2(dotp(p1,v),dotp(p1,u)),a2=atan2(dotp(p2,v),dotp(p2,u));if(a1!=a2)return a1<a2;return dotp(p1,p1)<dotp(p2,p2);});for(int rev=0;rev<2;rev++){Mesh m;m.v.reserve(ids.size());m.f.reserve(ring.size()*2);m.v.push_back(P[a]);m.v.push_back(P[b]);double amp=jit?epsv*1e-6:0;for(int i=0;i<(int)ring.size();i++){V p=P[ring[i]];if(amp>0){double s=sin((i+1)*12.9898),c=cos((i+1)*78.233);p=p+u*(amp*s)+v*(amp*c)+axis*(amp*.17*sin((i+1)*37.719));}m.v.push_back(p);}int R=ring.size();for(int i=0;i<R;i++){int x=2+i,y=2+(i+1)%R;m.f.push_back({0,x,y});m.f.push_back({1,y,x});}if(valid(m)){out.v.swap(m.v);out.f.swap(m.f);return true;}reverse(ring.begin(),ring.end());}}return false;}
static bool solve_net(Mesh&out){auto ids=best_net();if(ids.size()<4||ids.size()>=P.size())return false;Mesh m;if(!make_bipyramid(ids,m,false)&&!make_bipyramid(ids,m,true))return false;if(m.v.size()>=P.size())return false;if(!cover(P,m.v,epsv*1.0008))return false;if(!cover(m.v,P,epsv*1.0008))return false;out.v.swap(m.v);out.f.swap(m.f);return true;}
static void write(const Mesh&m){static char buf[1<<20];setvbuf(stdout,buf,_IOFBF,sizeof(buf));printf("%d %d\n",(int)m.v.size(),(int)m.f.size());for(auto&p:m.v)printf("v %.15g %.15g %.15g\n",p.x,p.y,p.z);for(auto&f:m.f)printf("f %d %d %d\n",f.a+1,f.b+1,f.c+1);} 
int main_c(){T0=chrono::steady_clock::now();In in;N=in.ni();M=in.ni();if(N<=0||M<=0)return 0;P.resize(N);for(int i=0;i<N;i++){in.ws();if(*in.p=='v'||*in.p=='V')++in.p;P[i].x=in.nd();P[i].y=in.nd();P[i].z=in.nd();}Fin.resize(M);for(int i=0;i<M;i++){in.ws();if(*in.p=='f'||*in.p=='F')++in.p;Fin[i].a=in.ni()-1;Fin[i].b=in.ni()-1;Fin[i].c=in.ni()-1;}bounds();Mesh out;if(solve_net(out)){write(out);return 0;}out.v=P;out.f=Fin;write(out);return 0;}
}

struct WVec{double x,y,z;};struct WFace{int a,b,c;};
static inline WVec operator-(WVec a,WVec b){return{a.x-b.x,a.y-b.y,a.z-b.z};}
static inline double wdot(WVec a,WVec b){return a.x*b.x+a.y*b.y+a.z*b.z;}static inline double wdist2(WVec a,WVec b){auto d=a-b;return wdot(d,d);}static inline uint64_t wek(int a,int b){if(a>b)swap(a,b);return (uint64_t)(uint32_t)a<<32|(uint32_t)b;}
struct WMesh{vector<WVec>v;vector<WFace>f;bool ok=false;};
struct WParser{const char*p;WParser(const string&s):p(s.c_str()){}void ws(){while(*p==' '||*p=='\n'||*p=='\r'||*p=='\t')++p;}int ni(){ws();int s=1;if(*p=='-')s=-1,++p;int x=0;while(*p>='0'&&*p<='9')x=x*10+*p++-'0';return s*x;}double nd(){ws();char*e;double x=strtod(p,&e);p=e;return x;}};
static bool parse_in(const string&s,WMesh&m){WParser r(s);int n=r.ni(),f=r.ni();if(n<=0||f<0)return false;m.v.resize(n);m.f.resize(f);for(int i=0;i<n;i++){r.ws();if(*r.p=='v'||*r.p=='V')++r.p;m.v[i]={r.nd(),r.nd(),r.nd()};}for(int i=0;i<f;i++){r.ws();if(*r.p=='f'||*r.p=='F')++r.p;m.f[i]={r.ni()-1,r.ni()-1,r.ni()-1};}m.ok=true;return true;}
static bool parse_out(const string&s,WMesh&m){return parse_in(s,m);} 
struct WHash{double c;WVec mn;unordered_map<long long,vector<int>>g;long long key(long long a,long long b,long long c)const{return a*73856093ll^b*19349663ll^c*83492791ll;}array<long long,3> id(WVec p)const{return{(long long)floor((p.x-mn.x)/c),(long long)floor((p.y-mn.y)/c),(long long)floor((p.z-mn.z)/c)};}void build(const vector<WVec>&v,WVec m,double cc){mn=m;c=max(cc,1e-12);g.clear();g.reserve(v.size()*2+3);for(int i=0;i<(int)v.size();i++){auto a=id(v[i]);g[key(a[0],a[1],a[2])].push_back(i);}}bool near(WVec p,const vector<WVec>&v,double r)const{double r2=r*r;auto a=id(p);for(long long dx=-1;dx<=1;dx++)for(long long dy=-1;dy<=1;dy++)for(long long dz=-1;dz<=1;dz++){auto it=g.find(key(a[0]+dx,a[1]+dy,a[2]+dz));if(it==g.end())continue;for(int j:it->second)if(wdist2(p,v[j])<=r2)return true;}return false;}};

static inline double wcget(const WVec&p,int i){return i==0?p.x:(i==1?p.y:p.z);} 
static void rview(const WMesh&m,const WVec&mn,const WVec&mx,int ax,int sg,int R,vector<float>&D,vector<unsigned char>&O){
    D.assign(R*R,-1e30f);O.assign(R*R,0);
    int u=ax==0?1:0,v=ax==2?1:2; if(ax==1){u=0;v=2;}
    double du=max(1e-12,wcget(mx,u)-wcget(mn,u)),dv=max(1e-12,wcget(mx,v)-wcget(mn,v)),dz=max(1e-12,wcget(mx,ax)-wcget(mn,ax));
    int step=max(1,(int)m.f.size()/50000);
    for(int fi=0;fi<(int)m.f.size();fi+=step){auto f=m.f[fi]; if(f.a<0||f.b<0||f.c<0||f.a>=(int)m.v.size()||f.b>=(int)m.v.size()||f.c>=(int)m.v.size())continue; WVec pp[3]={m.v[f.a],m.v[f.b],m.v[f.c]}; double x[3],y[3],z[3]; for(int i=0;i<3;i++){x[i]=(wcget(pp[i],u)-wcget(mn,u))/du*(R-1);y[i]=(wcget(pp[i],v)-wcget(mn,v))/dv*(R-1);z[i]=sg*(wcget(pp[i],ax)-wcget(mn,ax))/dz;} double den=(y[1]-y[2])*(x[0]-x[2])+(x[2]-x[1])*(y[0]-y[2]); if(fabs(den)<1e-12)continue; int x0=max(0,(int)floor(min({x[0],x[1],x[2]}))),x1=min(R-1,(int)ceil(max({x[0],x[1],x[2]}))),y0=max(0,(int)floor(min({y[0],y[1],y[2]}))),y1=min(R-1,(int)ceil(max({y[0],y[1],y[2]}))); for(int iy=y0;iy<=y1;iy++){double py=iy+.5;for(int ix=x0;ix<=x1;ix++){double px=ix+.5;double a=((y[1]-y[2])*(px-x[2])+(x[2]-x[1])*(py-y[2]))/den,b=((y[2]-y[0])*(px-x[2])+(x[0]-x[2])*(py-y[2]))/den,c=1-a-b;if(a>=-1e-9&&b>=-1e-9&&c>=-1e-9){double zz=a*z[0]+b*z[1]+c*z[2];int k=iy*R+ix;if(zz>D[k])D[k]=zz,O[k]=1;}}}
    }
}
static bool vproxy_ok(const WMesh&o,const WMesh&m){
    if(o.v.size()<40||m.v.size()>o.v.size()*85/100)return true;
    WVec mn=o.v[0],mx=o.v[0];for(auto&p:o.v){mn.x=min(mn.x,p.x);mn.y=min(mn.y,p.y);mn.z=min(mn.z,p.z);mx.x=max(mx.x,p.x);mx.y=max(mx.y,p.y);mx.z=max(mx.z,p.z);}int R=o.f.size()>400000?48:64;vector<float>a,b;vector<unsigned char>oa,ob;double worst=1,avg=0;int vc=0;for(int ax=0;ax<3;ax++)for(int sg=-1;sg<=1;sg+=2){rview(o,mn,mx,ax,sg,R,a,oa);rview(m,mn,mx,ax,sg,R,b,ob);int uni=0,inter=0,cnt=0;double dd=0;for(int i=0;i<R*R;i++){if(oa[i]||ob[i])uni++;if(oa[i]&&ob[i]){inter++;dd+=fabs((double)a[i]-b[i]);cnt++;}}double iou=uni?double(inter)/uni:1.0;double dep=cnt?max(0.0,1.0-dd/cnt):1.0;double sc=.65*iou+.35*dep;worst=min(worst,sc);avg+=sc;vc++;}avg/=max(1,vc);double hard=m.v.size()*5<o.v.size()?.90:.86;return worst>=hard||(worst>=.82&&avg>=.935);
}
static bool validate_candidate(const WMesh&orig,const WMesh&m){if(!m.ok||m.v.empty()||m.f.empty()||m.v.size()>orig.v.size())return false;WVec mn=orig.v[0],mx=orig.v[0];for(auto&p:orig.v){mn.x=min(mn.x,p.x);mn.y=min(mn.y,p.y);mn.z=min(mn.z,p.z);mx.x=max(mx.x,p.x);mx.y=max(mx.y,p.y);mx.z=max(mx.z,p.z);}double dx=mx.x-mn.x,dy=mx.y-mn.y,dz=mx.z-mn.z,diag=sqrt(dx*dx+dy*dy+dz*dz);double eps=max(1e-9,0.05*diag*1.0005);vector<uint64_t>ed;ed.reserve(m.f.size()*3);set<array<int,3>>tris;for(auto f:m.f){if(f.a<0||f.b<0||f.c<0||f.a>=(int)m.v.size()||f.b>=(int)m.v.size()||f.c>=(int)m.v.size())return false;if(f.a==f.b||f.a==f.c||f.b==f.c)return false;array<int,3>t={f.a,f.b,f.c};sort(t.begin(),t.end());if(!tris.insert(t).second)return false;ed.push_back(wek(f.a,f.b));ed.push_back(wek(f.b,f.c));ed.push_back(wek(f.c,f.a));}sort(ed.begin(),ed.end());for(size_t i=0;i<ed.size();){size_t j=i+1;while(j<ed.size()&&ed[j]==ed[i])++j;if(j-i!=2)return false;i=j;}WHash hg;hg.build(m.v,mn,eps);for(auto&p:orig.v)if(!hg.near(p,m.v,eps))return false;WHash ho;ho.build(orig.v,mn,eps);for(auto&p:m.v)if(!ho.near(p,orig.v,eps))return false;if(!vproxy_ok(orig,m))return false;return true;}
static string write_original(const WMesh&m){string o;char line[160];o.append(line,snprintf(line,sizeof(line),"%d %d\n",(int)m.v.size(),(int)m.f.size()));for(auto&p:m.v)o.append(line,snprintf(line,sizeof(line),"v %.15g %.15g %.15g\n",p.x,p.y,p.z));for(auto f:m.f)o.append(line,snprintf(line,sizeof(line),"f %d %d %d\n",f.a+1,f.b+1,f.c+1));return o;}
static string run_solver(const string&inp,int which){fflush(stdout);int s0=dup(0),s1=dup(1);if(s0<0||s1<0)return"";FILE*in=tmpfile();FILE*out=tmpfile();if(!in||!out)return"";fwrite(inp.data(),1,inp.size(),in);fflush(in);rewind(in);dup2(fileno(in),0);dup2(fileno(out),1);if(which==1)A::main_a();else if(which==2){B::TL=8.9;B::main_b();}else C::main_c();fflush(stdout);dup2(s1,1);dup2(s0,0);close(s0);close(s1);string r;char buf[1<<16];rewind(out);size_t n;while((n=fread(buf,1,sizeof(buf),out))>0)r.append(buf,buf+n);fclose(in);fclose(out);return r;}
int main(){string input;input.reserve(1<<22);char buf[1<<16];ssize_t n;while((n=read(0,buf,sizeof(buf)))>0)input.append(buf,buf+n);WMesh orig;if(!parse_in(input,orig))return 0;string bests=write_original(orig);WMesh best=orig;auto consider=[&](const string&s){WMesh m;if(parse_out(s,m)&&validate_candidate(orig,m)){if(m.v.size()<best.v.size()||(m.v.size()==best.v.size()&&m.f.size()<best.f.size())){best=m;bests=s;}}};auto t0=chrono::steady_clock::now();string s1=run_solver(input,1);consider(s1);double used=chrono::duration<double>(chrono::steady_clock::now()-t0).count();if(used<10.2&&!(best.v.size()<=100&&best.v.size()*200<=orig.v.size())){string s2=run_solver(input,2);consider(s2);}double used2=chrono::duration<double>(chrono::steady_clock::now()-t0).count();if(used2<14.0&&!(best.v.size()<=64&&best.v.size()*150<=orig.v.size())){string s3=run_solver(input,3);consider(s3);}fwrite(bests.data(),1,bests.size(),stdout);return 0;}
