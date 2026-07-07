#include <bits/stdc++.h>
using namespace std;

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
    static char buf[1<<20]; setvbuf(stdout,buf,_IOFBF,sizeof(buf));
    printf("%zu %zu\n",V.size(),G.size());
    for(const Vec3& p:V) printf("v %.10g %.10g %.10g\n",p.x,p.y,p.z);
    for(auto t:G) printf("f %d %d %d\n",t[0]+1,t[1]+1,t[2]+1);
}
int main(){
    startTime=chrono::steady_clock::now(); loadInput();
    if(!trySampleLikeSimplify()) decimate();
    vector<Vec3> outP; vector<array<int,3>> outF; bool ok=buildCompact(outP,outF); if(ok) ok=validateMesh(outP,outF); if(!ok) makeOriginalOutput(outP,outF); writeMesh(outP,outF); return 0;
}
