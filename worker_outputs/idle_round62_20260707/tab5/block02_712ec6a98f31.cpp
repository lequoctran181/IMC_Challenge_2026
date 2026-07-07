#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <queue>
#include <unordered_set>
#include <vector>
using namespace std;

struct Vec3{double x,y,z;};
static inline Vec3 operator+(const Vec3&a,const Vec3&b){return {a.x+b.x,a.y+b.y,a.z+b.z};}
static inline Vec3 operator-(const Vec3&a,const Vec3&b){return {a.x-b.x,a.y-b.y,a.z-b.z};}
static inline Vec3 operator*(const Vec3&a,double s){return {a.x*s,a.y*s,a.z*s};}
static inline double dotp(const Vec3&a,const Vec3&b){return a.x*b.x+a.y*b.y+a.z*b.z;}
static inline Vec3 crossp(const Vec3&a,const Vec3&b){return {a.y*b.z-a.z*b.y,a.z*b.x-a.x*b.z,a.x*b.y-a.y*b.x};}
static inline double n2(const Vec3&a){return dotp(a,a);} 
static inline double nv(const Vec3&a){return sqrt(max(0.0,n2(a)));}
struct Face{int v[3];};

static int N,M;
static vector<Vec3> P, Pin;
static vector<Face> F, Fin;
static vector<unsigned char> VA, FA;
static vector<vector<int>> VF;
static vector<double> R;
static vector<unsigned> ST;
static int aliveV, aliveF;
static double diagL, epsA, radMax;

static inline uint64_t keyEdge(int a,int b){ if(a>b) swap(a,b); return (uint64_t)(uint32_t)a<<32 | (uint32_t)b; }
static inline bool hasv(const Face&f,int x){return f.v[0]==x||f.v[1]==x||f.v[2]==x;}
static inline Vec3 fcross(const Face&f){return crossp(P[f.v[1]]-P[f.v[0]],P[f.v[2]]-P[f.v[0]]);} 
static inline array<int,3> sface(int a,int b,int c){array<int,3>s={a,b,c}; sort(s.begin(),s.end()); return s;}
static inline int third(const Face&f,int a,int b){for(int i=0;i<3;i++) if(f.v[i]!=a&&f.v[i]!=b) return f.v[i]; return -1;}

static void rebuildVF(){
    VF.assign(N,{});
    for(int i=0;i<M;i++) if(FA[i]) for(int j=0;j<3;j++) VF[F[i].v[j]].push_back(i);
}
static vector<int> inc(int v){
    vector<int> r; if(v<0||v>=N||!VA[v]) return r;
    for(int id:VF[v]) if(id>=0&&id<M&&FA[id]&&hasv(F[id],v)) r.push_back(id);
    sort(r.begin(),r.end()); r.erase(unique(r.begin(),r.end()),r.end()); return r;
}
static vector<int> nbr(int v){
    vector<int> r;
    for(int id:inc(v)){
        for(int k=0;k<3;k++){int x=F[id].v[k]; if(x!=v&&VA[x]) r.push_back(x);} 
    }
    sort(r.begin(),r.end()); r.erase(unique(r.begin(),r.end()),r.end()); return r;
}
static bool edgeFaces(int a,int b,int ef[2],int op[2]){
    int c=0;
    for(int id:inc(a)) if(hasv(F[id],b)){
        if(c>=2) return false;
        ef[c]=id; op[c]=third(F[id],a,b); c++;
    }
    return c==2 && op[0]>=0 && op[1]>=0 && op[0]!=op[1];
}
static bool linkOK(int a,int b,const int op[2]){
    vector<int>A=nbr(a),B=nbr(b),I;
    size_t i=0,j=0;
    while(i<A.size()&&j<B.size()){
        if(A[i]==B[j]){I.push_back(A[i]); i++; j++;}
        else if(A[i]<B[j]) i++; else j++;
    }
    sort(I.begin(),I.end()); I.erase(unique(I.begin(),I.end()),I.end());
    if(I.size()!=2) return false;
    return (I[0]==op[0]&&I[1]==op[1])||(I[0]==op[1]&&I[1]==op[0]);
}
static bool duplicateFaceAfter(int fid,int keep,int rem,int ef0,int ef1){
    Face nf=F[fid];
    for(int i=0;i<3;i++) if(nf.v[i]==rem) nf.v[i]=keep;
    auto S=sface(nf.v[0],nf.v[1],nf.v[2]);
    vector<int> pool=inc(keep), q=inc(rem); pool.insert(pool.end(),q.begin(),q.end());
    sort(pool.begin(),pool.end()); pool.erase(unique(pool.begin(),pool.end()),pool.end());
    for(int id:pool) if(id!=fid&&id!=ef0&&id!=ef1&&FA[id]){
        if(hasv(F[id],rem)) continue;
        if(sface(F[id].v[0],F[id].v[1],F[id].v[2])==S) return true;
    }
    return false;
}

struct Eval{bool ok=false; double cost=1e100, nr=0; int keep=-1,rem=-1;};
static Eval evalDir(int keep,int rem,const int ef[2],double cosMin,double planeMax){
    Eval e; e.keep=keep; e.rem=rem;
    double d=nv(P[keep]-P[rem]);
    e.nr=max(R[keep],R[rem]+d);
    if(e.nr>radMax) return e;
    vector<int> affected=inc(rem); vector<int> k=inc(keep); affected.insert(affected.end(),k.begin(),k.end());
    sort(affected.begin(),affected.end()); affected.erase(unique(affected.begin(),affected.end()),affected.end());
    double worstN=0, worstP=0, shrink=0; int changed=0;
    for(int id:affected) if(id!=ef[0]&&id!=ef[1]&&FA[id]){
        Face old=F[id], nf=old;
        bool touchesRem=false;
        for(int t=0;t<3;t++) if(nf.v[t]==rem){nf.v[t]=keep; touchesRem=true;}
        if(nf.v[0]==nf.v[1]||nf.v[0]==nf.v[2]||nf.v[1]==nf.v[2]) return e;
        Vec3 co=fcross(old); double ao=nv(co); if(!(ao>epsA)) return e;
        Vec3 cn=crossp(P[nf.v[1]]-P[nf.v[0]],P[nf.v[2]]-P[nf.v[0]]); double an=nv(cn); if(!(an>epsA)) return e;
        double nd=dotp(co,cn)/(ao*an); nd=max(-1.0,min(1.0,nd));
        if(nd<cosMin) return e;
        Vec3 n=co*(1.0/ao);
        double pd=0.0;
        if(touchesRem) pd=fabs(dotp(n,P[keep]-P[rem]));
        if(pd>planeMax) return e;
        if(touchesRem && duplicateFaceAfter(id,keep,rem,ef[0],ef[1])) return e;
        worstN=max(worstN,1.0-nd); worstP=max(worstP,pd/(diagL+1e-30));
        if(an<ao) shrink=max(shrink,1.0-an/ao);
        changed++;
    }
    if(changed==0) return e;
    Vec3 c0=fcross(F[ef[0]]), c1=fcross(F[ef[1]]); double a0=nv(c0),a1=nv(c1);
    double sharp=0.0; if(a0>epsA&&a1>epsA){double q=dotp(c0,c1)/(a0*a1); q=max(-1.0,min(1.0,q)); sharp=1.0-q;}
    e.cost=(d/(diagL+1e-30))*0.45 + (e.nr/(radMax+1e-30))*0.95 + 18.0*worstN + 3.0*worstP + 0.25*shrink + 0.12*sharp + 0.00015*changed;
    e.ok=true; return e;
}
static Eval evalEdge(int a,int b,double cosMin,double planeMax){
    Eval bad; if(a==b||!VA[a]||!VA[b]) return bad;
    int ef[2],op[2]; if(!edgeFaces(a,b,ef,op)) return bad; if(!linkOK(a,b,op)) return bad;
    Eval x=evalDir(a,b,ef,cosMin,planeMax), y=evalDir(b,a,ef,cosMin,planeMax);
    if(x.ok&&(!y.ok||x.cost<=y.cost)) return x;
    if(y.ok) return y;
    return bad;
}
struct Cand{double cost; int a,b; unsigned sa,sb; bool operator<(const Cand&o)const{return cost>o.cost;}};
static priority_queue<Cand> PQ;
static void pushEdge(int a,int b,double cosMin,double planeMax){
    if(a==b||!VA[a]||!VA[b]) return;
    Eval e=evalEdge(a,b,cosMin,planeMax); if(!e.ok) return;
    PQ.push({e.cost,a,b,ST[a],ST[b]});
}
static void collapse(const Eval&e){
    int keep=e.keep, rem=e.rem, ef[2],op[2];
    if(!edgeFaces(keep,rem,ef,op)) return;
    for(int i=0;i<2;i++) if(FA[ef[i]]){FA[ef[i]]=0; aliveF--;}
    vector<int> ir=inc(rem);
    for(int id:ir) if(FA[id]){
        for(int t=0;t<3;t++) if(F[id].v[t]==rem) F[id].v[t]=keep;
        VF[keep].push_back(id);
    }
    VA[rem]=0; aliveV--; R[keep]=e.nr; VF[rem].clear(); ST[keep]++; ST[rem]++;
}
static bool closedManifold(){
    unordered_set<uint64_t> once, twice; once.reserve((size_t)aliveF*3+10); twice.reserve((size_t)aliveF*3+10);
    for(int i=0;i<M;i++) if(FA[i]){
        Face f=F[i]; if(!VA[f.v[0]]||!VA[f.v[1]]||!VA[f.v[2]]) return false;
        if(f.v[0]==f.v[1]||f.v[0]==f.v[2]||f.v[1]==f.v[2]) return false;
        if(nv(fcross(f))<=epsA) return false;
        for(int k=0;k<3;k++){uint64_t e=keyEdge(f.v[k],f.v[(k+1)%3]); if(!once.insert(e).second) twice.insert(e);} 
    }
    if(once.size()!=twice.size() || once.size()*2!=(size_t)aliveF*3) return false;
    for(int v=0;v<N;v++) if(VA[v]){
        vector<int> in=inc(v); if(in.empty()) return false;
        vector<pair<int,int>> linkEdges;
        vector<int> verts;
        for(int id:in){int a=-1,b=-1; for(int k=0;k<3;k++){int x=F[id].v[k]; if(x!=v){ if(a<0)a=x; else b=x; }} if(a<0||b<0) return false; linkEdges.push_back({a,b}); verts.push_back(a); verts.push_back(b);} 
        sort(verts.begin(),verts.end()); verts.erase(unique(verts.begin(),verts.end()),verts.end());
        vector<int> deg(verts.size(),0);
        for(auto &p:linkEdges){int ia=lower_bound(verts.begin(),verts.end(),p.first)-verts.begin(); int ib=lower_bound(verts.begin(),verts.end(),p.second)-verts.begin(); deg[ia]++; deg[ib]++;}
        for(int d:deg) if(d!=2) return false;
    }
    return true;
}
static void runPass(double cosMin,double planeFrac,int target,int opLimit){
    double planeMax=planeFrac*diagL;
    while(!PQ.empty()) PQ.pop(); unordered_set<uint64_t> seen; seen.reserve((size_t)aliveF*3+10);
    for(int i=0;i<M;i++) if(FA[i]) for(int k=0;k<3;k++){int a=F[i].v[k],b=F[i].v[(k+1)%3]; uint64_t kk=keyEdge(a,b); if(seen.insert(kk).second) pushEdge(a,b,cosMin,planeMax);} 
    int ops=0;
    while(!PQ.empty()&&aliveV>target&&ops<opLimit){
        Cand c=PQ.top(); PQ.pop();
        if(!VA[c.a]||!VA[c.b]||ST[c.a]!=c.sa||ST[c.b]!=c.sb) continue;
        Eval e=evalEdge(c.a,c.b,cosMin,planeMax); if(!e.ok) continue;
        if(e.cost>c.cost*1.5+1e-12){pushEdge(c.a,c.b,cosMin,planeMax); continue;}
        collapse(e); ops++;
        vector<int> nb=nbr(e.keep);
        for(int x:nb) pushEdge(e.keep,x,cosMin,planeMax);
        if((ops&511)==0) rebuildVF();
    }
    rebuildVF();
}
static void resetOriginal(){P=Pin; F=Fin; VA.assign(N,1); FA.assign(M,1); VF.assign(N,{}); R.assign(N,0.0); ST.assign(N,1); aliveV=N; aliveF=M; rebuildVF();}
static void emit(){
    vector<int> id(N,-1); vector<Vec3> V; V.reserve(aliveV);
    for(int i=0;i<N;i++) if(VA[i]){id[i]=(int)V.size()+1; V.push_back(P[i]);}
    vector<Face> OF; OF.reserve(aliveF);
    for(int i=0;i<M;i++) if(FA[i]){Face f=F[i]; if(id[f.v[0]]>0&&id[f.v[1]]>0&&id[f.v[2]]>0) OF.push_back({id[f.v[0]],id[f.v[1]],id[f.v[2]]});}
    cout.setf(ios::fixed); cout.precision(10);
    cout<<V.size()<<" "<<OF.size()<<"\n";
    for(auto&p:V) cout<<"v "<<p.x<<" "<<p.y<<" "<<p.z<<"\n";
    for(auto&f:OF) cout<<"f "<<f.v[0]<<" "<<f.v[1]<<" "<<f.v[2]<<"\n";
}
int main(){
    ios::sync_with_stdio(false); cin.tie(nullptr);
    if(!(cin>>N>>M)) return 0;
    P.resize(N); Pin.resize(N); F.resize(M); Fin.resize(M);
    Vec3 mn{1e100,1e100,1e100}, mx{-1e100,-1e100,-1e100}; string t;
    for(int i=0;i<N;i++){cin>>t>>P[i].x>>P[i].y>>P[i].z; Pin[i]=P[i]; mn.x=min(mn.x,P[i].x); mn.y=min(mn.y,P[i].y); mn.z=min(mn.z,P[i].z); mx.x=max(mx.x,P[i].x); mx.y=max(mx.y,P[i].y); mx.z=max(mx.z,P[i].z);} 
    for(int i=0;i<M;i++){cin>>t>>F[i].v[0]>>F[i].v[1]>>F[i].v[2]; --F[i].v[0];--F[i].v[1];--F[i].v[2]; Fin[i]=F[i];}
    diagL=nv(mx-mn); if(!(diagL>0)) diagL=1.0; epsA=1e-20*diagL*diagL; radMax=0.0490*diagL;
    resetOriginal();
    if(N<=8){emit(); return 0;}
    int t1=max(8,(int)ceil(N*0.32));
    int t2=max(8,(int)ceil(N*0.18));
    int t3=max(8,(int)ceil(N*0.10));
    int t4=max(8,(int)ceil(N*0.065));
    runPass(0.985,0.018,t1,N);
    runPass(0.955,0.028,t2,N);
    runPass(0.910,0.040,t3,2*N);
    runPass(0.860,0.048,t4,3*N);
    if(!closedManifold()) resetOriginal();
    emit();
    return 0;
}