#include <bits/stdc++.h>
using namespace std;

struct Vec3{double x=0,y=0,z=0;};
struct Face{int v[3];};
static inline Vec3 operator+(const Vec3&a,const Vec3&b){return {a.x+b.x,a.y+b.y,a.z+b.z};}
static inline Vec3 operator-(const Vec3&a,const Vec3&b){return {a.x-b.x,a.y-b.y,a.z-b.z};}
static inline Vec3 operator*(const Vec3&a,double s){return {a.x*s,a.y*s,a.z*s};}
static inline double dot3(const Vec3&a,const Vec3&b){return a.x*b.x+a.y*b.y+a.z*b.z;}
static inline Vec3 cross3(const Vec3&a,const Vec3&b){return {a.y*b.z-a.z*b.y,a.z*b.x-a.x*b.z,a.x*b.y-a.y*b.x};}
static inline double norm2(const Vec3&a){return dot3(a,a);} 
static inline double norm3(const Vec3&a){return sqrt(max(0.0,norm2(a)));}
static inline uint64_t ekey(int a,int b){if(a>b) swap(a,b); return (uint64_t)(uint32_t)a<<32 | (uint32_t)b;}
static inline array<int,3> tkey(int a,int b,int c){array<int,3> r{a,b,c}; sort(r.begin(),r.end()); return r;}

struct FastInput{
    vector<char> buf; char* p=nullptr;
    FastInput(){buf.reserve(1<<27); char tmp[1<<16]; size_t n; while((n=fread(tmp,1,sizeof(tmp),stdin))>0) buf.insert(buf.end(),tmp,tmp+n); buf.push_back('\0'); p=buf.data();}
    inline void skip(){while(*p==' '||*p=='\n'||*p=='\r'||*p=='\t') ++p;}
    inline bool eof(){skip(); return *p=='\0';}
    long nextLong(){skip(); return strtol(p,&p,10);} 
    double nextDouble(){skip(); return strtod(p,&p);} 
    void skipToken(){skip(); if((*p>='a'&&*p<='z')||(*p>='A'&&*p<='Z')){while(*p&&*p!=' '&&*p!='\n'&&*p!='\r'&&*p!='\t') ++p;}}
};

static int N=0,M=0;
static vector<Vec3> P, origP;
static vector<Face> F, origF;
static vector<unsigned char> aliveV, aliveF, fixedV;
static vector<double> coverR;
static vector<vector<int>> adj;
static vector<int> markA, markB;
static int stampA=1, stampB=1;
static int aliveVertices=0, aliveFaces=0;
static Vec3 bbMin, bbMax, bbExt;
static double diagLen=1.0;
static chrono::steady_clock::time_point T0;

static inline double elapsed(){return chrono::duration<double>(chrono::steady_clock::now()-T0).count();}
static inline double coord(const Vec3&p,int a){return a==0?p.x:(a==1?p.y:p.z);} 
static inline bool faceContains(const Face&f,int v){return f.v[0]==v||f.v[1]==v||f.v[2]==v;}
static inline bool containsF(int fid,int v){return faceContains(F[fid],v);} 
static int thirdVertex(int fid,int a,int b){const Face&f=F[fid]; for(int i=0;i<3;i++){int x=f.v[i]; if(x!=a&&x!=b) return x;} return -1;}
static inline Vec3 triCross(const Vec3&a,const Vec3&b,const Vec3&c){return cross3(b-a,c-a);} 
static inline Vec3 currentCross(const Face&f){return triCross(P[f.v[0]],P[f.v[1]],P[f.v[2]]);} 
static inline Vec3 movedPoint(int id,int mv,const Vec3&pos){return id==mv?pos:P[id];}
static inline Vec3 movedCross(int a,int b,int c,int mv,const Vec3&pos){return triCross(movedPoint(a,mv,pos),movedPoint(b,mv,pos),movedPoint(c,mv,pos));}

static double signedVol6(const vector<Vec3>&V,const vector<Face>&G){long double s=0; for(const auto&f:G) s+=(long double)dot3(V[f.v[0]],cross3(V[f.v[1]],V[f.v[2]])); return (double)s;}
static void orientLikeOriginal(vector<Vec3>&V, vector<Face>&G){double a=signedVol6(origP,origF), b=signedVol6(V,G); if(a*b<0) for(auto &f:G) swap(f.v[1],f.v[2]);}

static bool readInput(){
    FastInput in; if(in.eof()) return false; N=(int)in.nextLong(); M=(int)in.nextLong();
    P.resize(N); origP.resize(N); aliveV.assign(N,1); fixedV.assign(N,0); coverR.assign(N,0.0);
    bbMin={1e100,1e100,1e100}; bbMax={-1e100,-1e100,-1e100};
    for(int i=0;i<N;i++){
        in.skipToken(); P[i].x=in.nextDouble(); P[i].y=in.nextDouble(); P[i].z=in.nextDouble(); origP[i]=P[i];
        bbMin.x=min(bbMin.x,P[i].x); bbMin.y=min(bbMin.y,P[i].y); bbMin.z=min(bbMin.z,P[i].z);
        bbMax.x=max(bbMax.x,P[i].x); bbMax.y=max(bbMax.y,P[i].y); bbMax.z=max(bbMax.z,P[i].z);
    }
    bbExt=bbMax-bbMin; diagLen=norm3(bbExt); if(!(diagLen>0)) diagLen=1.0;
    F.resize(M); origF.resize(M); aliveF.assign(M,1); vector<int> deg(N,0);
    for(int i=0;i<M;i++){
        in.skipToken(); int a=(int)in.nextLong()-1,b=(int)in.nextLong()-1,c=(int)in.nextLong()-1;
        if(a<0||b<0||c<0||a>=N||b>=N||c>=N) return false;
        F[i].v[0]=a; F[i].v[1]=b; F[i].v[2]=c; origF[i]=F[i];
        ++deg[a]; ++deg[b]; ++deg[c];
    }
    adj.assign(N,{}); for(int i=0;i<N;i++) adj[i].reserve(deg[i]+4);
    for(int i=0;i<M;i++){adj[F[i].v[0]].push_back(i); adj[F[i].v[1]].push_back(i); adj[F[i].v[2]].push_back(i);} 
    aliveVertices=N; aliveFaces=M; markA.assign(N,0); markB.assign(N,0);
    for(int ax=0;ax<3;ax++){
        int mn=0,mx=0; for(int i=1;i<N;i++){ if(coord(P[i],ax)<coord(P[mn],ax)) mn=i; if(coord(P[i],ax)>coord(P[mx],ax)) mx=i; }
        fixedV[mn]=1; fixedV[mx]=1;
    }
    return true;
}

struct ShapeStats{
    double smooth10=0, smooth30=0, sharp22=0, sharp45=0, bad=0;
};
static Vec3 faceUnit(int fid){Vec3 cr=currentCross(F[fid]); double l=norm3(cr); return l>1e-300?cr*(1.0/l):Vec3{0,0,0};}
static bool edgeFaces(int u,int v,int ef[2],int op[2]){
    if(adj[u].size()>adj[v].size()) swap(u,v);
    int cnt=0; for(int fid:adj[u]) if(aliveF[fid]&&containsF(fid,u)&&containsF(fid,v)){
        if(cnt>=2) return false; ef[cnt]=fid; op[cnt]=thirdVertex(fid,u,v); ++cnt;
    }
    return cnt==2 && op[0]>=0 && op[1]>=0 && op[0]!=op[1];
}
static ShapeStats inspectShape(){
    ShapeStats s; if(N<50||M<30) return s;
    const int target=60000, stride=max(1,M/target), maxs=120000;
    const double c10=cos(10.0*M_PI/180.0), c30=cos(30.0*M_PI/180.0), c22=cos(22.0*M_PI/180.0), c45=cos(45.0*M_PI/180.0);
    int sampled=0, sm10=0, sm30=0, sh22=0, sh45=0, bad=0;
    for(int fid=0;fid<M && sampled<maxs;fid+=stride){
        const Face&f=F[fid]; int e[3][2]={{f.v[0],f.v[1]},{f.v[1],f.v[2]},{f.v[2],f.v[0]}};
        for(int k=0;k<3 && sampled<maxs;k++){
            int ef[2],op[2]; if(!edgeFaces(e[k][0],e[k][1],ef,op)){bad++; continue;}
            Vec3 n0=faceUnit(ef[0]), n1=faceUnit(ef[1]); double d=dot3(n0,n1); d=max(-1.0,min(1.0,d));
            sampled++; if(d>c10) sm10++; if(d>c30) sm30++; if(d<c22) sh22++; if(d<c45) sh45++;
        }
    }
    if(sampled>0){s.smooth10=(double)sm10/sampled; s.smooth30=(double)sm30/sampled; s.sharp22=(double)sh22/sampled; s.sharp45=(double)sh45/sampled; s.bad=(double)bad/(sampled+bad);} return s;
}

static bool linkOK(int u,int v,const int op[2]){
    if(++stampA>2000000000){fill(markA.begin(),markA.end(),0); stampA=1;} if(++stampB>2000000000){fill(markB.begin(),markB.end(),0); stampB=1;}
    for(int fid:adj[u]) if(aliveF[fid]&&containsF(fid,u)){const Face&f=F[fid]; for(int k=0;k<3;k++){int x=f.v[k]; if(x!=u&&x!=v) markA[x]=stampA;}}
    int inter=0; for(int fid:adj[v]) if(aliveF[fid]&&containsF(fid,v)){const Face&f=F[fid]; for(int k=0;k<3;k++){
        int x=f.v[k]; if(x==u||x==v||markA[x]!=stampA) continue; if(x!=op[0]&&x!=op[1]) return false; if(markB[x]!=stampB){markB[x]=stampB; ++inter;}
    }}
    return inter==2 && markB[op[0]]==stampB && markB[op[1]]==stampB;
}

struct Params{double cover=0, flatCover=0, plane=0, minCos=1, flatPlane=0, flatCos=1, area=0; bool midpoint=false;};
struct Eval{bool ok=false, flat=false; double cost=1e100, rad=0; Vec3 pos;};

static bool duplicateFace(int keep,int rem,int fid,int a,int b,int c,int e0,int e1){
    int best=keep; if(adj[a].size()<adj[best].size()) best=a; if(adj[b].size()<adj[best].size()) best=b; if(adj[c].size()<adj[best].size()) best=c;
    array<int,3> key=tkey(a,b,c);
    for(int g:adj[best]){
        if(!aliveF[g]||g==fid||g==e0||g==e1) continue; if(containsF(g,rem)) continue;
        if(tkey(F[g].v[0],F[g].v[1],F[g].v[2])==key) return true;
    }
    return false;
}

static Eval evalEndpoint(int keep,int rem,const int ef[2],const Params&p){
    Eval r; if(fixedV[rem]) return r; r.pos=P[keep]; double ed=norm3(P[keep]-P[rem]); r.rad=max(coverR[keep],coverR[rem]+ed);
    if(r.rad>p.flatCover) return r;
    double maxPlane=0, bend=0, shrink=0; int changed=0; bool flat=true;
    for(int fid:adj[rem]) if(aliveF[fid]&&containsF(fid,rem)){
        if(fid==ef[0]||fid==ef[1]) continue; if(containsF(fid,keep)) return r;
        Face old=F[fid], nf=old; for(int k=0;k<3;k++) if(nf.v[k]==rem) nf.v[k]=keep;
        if(nf.v[0]==nf.v[1]||nf.v[0]==nf.v[2]||nf.v[1]==nf.v[2]) return r;
        Vec3 oc=currentCross(old), nc=currentCross(nf); double oa=norm3(oc), na=norm3(nc);
        if(!(oa>p.area)||!(na>p.area)||na<oa*1e-8) return r;
        double cs=dot3(oc,nc)/(oa*na); cs=max(-1.0,min(1.0,cs)); if(cs<p.minCos) return r;
        Vec3 n=oc*(1.0/oa); double pd=fabs(dot3(n,P[keep]-P[old.v[0]]));
        if(pd>p.plane && r.rad>p.cover) return r;
        if(duplicateFace(keep,rem,fid,nf.v[0],nf.v[1],nf.v[2],ef[0],ef[1])) return r;
        maxPlane=max(maxPlane,pd); bend=max(bend,1.0-cs); shrink=max(shrink,max(0.0,1.0-na/oa));
        if(pd>p.flatPlane || cs<p.flatCos) flat=false; changed++;
    }
    if(changed==0) return r;
    if(!flat){ if(r.rad>p.cover) return r; if(maxPlane>p.plane) return r; }
    r.ok=true; r.flat=flat; double hc=p.cover>0?r.rad/p.cover:0, pc=p.plane>0?maxPlane/p.plane:0; r.cost=(flat?-1000.0:0.0)+0.75*hc+0.85*pc+260.0*bend+0.02*shrink+0.0007*adj[rem].size(); return r;
}

static Eval evalMove(int keep,int rem,const int ef[2],const Params&p,const Vec3&pos){
    Eval r; if(fixedV[keep]||fixedV[rem]) return r; r.pos=pos; r.rad=max(coverR[keep]+norm3(pos-P[keep]),coverR[rem]+norm3(pos-P[rem])); if(r.rad>p.cover) return r;
    static vector<int> ids; ids.clear(); ids.reserve(adj[keep].size()+adj[rem].size());
    for(int fid:adj[keep]) if(aliveF[fid]&&fid!=ef[0]&&fid!=ef[1]&&containsF(fid,keep)) ids.push_back(fid);
    for(int fid:adj[rem]) if(aliveF[fid]&&fid!=ef[0]&&fid!=ef[1]&&containsF(fid,rem)) ids.push_back(fid);
    sort(ids.begin(),ids.end()); ids.erase(unique(ids.begin(),ids.end()),ids.end()); if(ids.empty()) return r;
    vector<array<int,3>> newKeys; newKeys.reserve(ids.size()); double maxPlane=0,bend=0,shrink=0;
    for(int fid:ids){ bool fromRem=containsF(fid,rem); Face old=F[fid], nf=old; for(int k=0;k<3;k++) if(nf.v[k]==rem) nf.v[k]=keep;
        if(nf.v[0]==nf.v[1]||nf.v[0]==nf.v[2]||nf.v[1]==nf.v[2]) return r;
        Vec3 oc=currentCross(old); Vec3 nc=movedCross(nf.v[0],nf.v[1],nf.v[2],keep,pos); double oa=norm3(oc), na=norm3(nc);
        if(!(oa>p.area)||!(na>p.area)||na<oa*1e-8) return r;
        double cs=dot3(oc,nc)/(oa*na); cs=max(-1.0,min(1.0,cs)); if(cs<p.minCos) return r;
        Vec3 n=oc*(1.0/oa); double pd=fabs(dot3(n,pos-P[old.v[0]])); if(pd>p.plane) return r;
        if(fromRem && duplicateFace(keep,rem,fid,nf.v[0],nf.v[1],nf.v[2],ef[0],ef[1])) return r;
        newKeys.push_back(tkey(nf.v[0],nf.v[1],nf.v[2])); maxPlane=max(maxPlane,pd); bend=max(bend,1.0-cs); shrink=max(shrink,max(0.0,1.0-na/oa));
    }
    sort(newKeys.begin(),newKeys.end()); for(size_t i=1;i<newKeys.size();i++) if(newKeys[i]==newKeys[i-1]) return r;
    r.ok=true; r.flat=false; double hc=p.cover>0?r.rad/p.cover:0, pc=p.plane>0?maxPlane/p.plane:0; r.cost=0.35+0.95*hc+0.95*pc+310.0*bend+0.03*shrink+0.0008*ids.size(); return r;
}

static void compactAdj(int v){
    vector<int>&a=adj[v]; if(a.empty()) return; sort(a.begin(),a.end()); a.erase(unique(a.begin(),a.end()),a.end());
    size_t alive=0; for(int fid:a) if(aliveF[fid]&&containsF(fid,v)) ++alive; if(alive*3+64>=a.size()) return;
    vector<int> b; b.reserve(alive+8); for(int fid:a) if(aliveF[fid]&&containsF(fid,v)) b.push_back(fid); a.swap(b);
}

using PQNode=tuple<double,int,int>;
static void pushIncident(int v, priority_queue<PQNode,vector<PQNode>,greater<PQNode>>&pq){
    for(int fid:adj[v]) if(aliveF[fid]&&containsF(fid,v)){
        const Face&f=F[fid]; for(int k=0;k<3;k++){int a=f.v[k],b=f.v[(k+1)%3]; if(a==v||b==v){int u=a^b^v; (void)u; if(aliveV[a]&&aliveV[b]&&a!=b) pq.emplace(norm3(P[a]-P[b]),a,b);}}
    }
}
static void applyCollapse(int keep,int rem,const int ef[2],const Eval&e, priority_queue<PQNode,vector<PQNode>,greater<PQNode>>&pq){
    for(int i=0;i<2;i++) if(aliveF[ef[i]]){aliveF[ef[i]]=0; --aliveFaces;}
    for(int fid:adj[rem]) if(aliveF[fid]&&containsF(fid,rem)){
        Face&f=F[fid]; for(int k=0;k<3;k++) if(f.v[k]==rem) f.v[k]=keep; adj[keep].push_back(fid);
    }
    aliveV[rem]=0; --aliveVertices; P[keep]=e.pos; coverR[keep]=e.rad; adj[rem].clear(); compactAdj(keep); pushIncident(keep,pq);
}

static vector<uint64_t> aliveEdges(){
    vector<uint64_t> es; es.reserve((size_t)aliveFaces*3);
    for(int i=0;i<M;i++) if(aliveF[i]){const Face&f=F[i]; if(aliveV[f.v[0]]&&aliveV[f.v[1]]&&aliveV[f.v[2]]){es.push_back(ekey(f.v[0],f.v[1])); es.push_back(ekey(f.v[1],f.v[2])); es.push_back(ekey(f.v[2],f.v[0]));}}
    sort(es.begin(),es.end()); es.erase(unique(es.begin(),es.end()),es.end()); return es;
}

static bool tryCollapse(int u,int v,const Params&p, priority_queue<PQNode,vector<PQNode>,greater<PQNode>>&pq){
    if(u==v||!aliveV[u]||!aliveV[v]) return false; int ef[2],op[2]; if(!edgeFaces(u,v,ef,op)) return false; if(!linkOK(u,v,op)) return false;
    Eval uv=evalEndpoint(u,v,ef,p), vu=evalEndpoint(v,u,ef,p);
    if(p.midpoint){
        Vec3 mid=(P[u]+P[v])*0.5; Eval um=evalMove(u,v,ef,p,mid), vm=evalMove(v,u,ef,p,mid);
        if(um.ok&&(!uv.ok||um.cost<uv.cost)) uv=um; if(vm.ok&&(!vu.ok||vm.cost<vu.cost)) vu=vm;
        Vec3 q1=P[u]*0.70+P[v]*0.30, q2=P[u]*0.30+P[v]*0.70;
        Eval uq=evalMove(u,v,ef,p,q1), vq=evalMove(v,u,ef,p,q2);
        if(uq.ok&&(!uv.ok||uq.cost<uv.cost)) uv=uq; if(vq.ok&&(!vu.ok||vq.cost<vu.cost)) vu=vq;
    }
    if(!uv.ok&&!vu.ok) return false;
    if(vu.ok&&(!uv.ok||vu.cost<uv.cost)) applyCollapse(v,u,ef,vu,pq); else applyCollapse(u,v,ef,uv,pq); return true;
}

static void runPass(const Params&p,double timeLimit,int target){
    vector<uint64_t> es=aliveEdges(); priority_queue<PQNode,vector<PQNode>,greater<PQNode>> pq;
    for(uint64_t k:es){int a=(int)(k>>32), b=(int)(uint32_t)k; if(aliveV[a]&&aliveV[b]) pq.emplace(norm3(P[a]-P[b]),a,b);} vector<uint64_t>().swap(es);
    long long pops=0; while(!pq.empty()&&aliveVertices>target){
        if((++pops&4095)==0 && elapsed()>timeLimit) break; auto [l,u,v]=pq.top(); pq.pop();
        if(u==v||!aliveV[u]||!aliveV[v]) continue; int ef[2],op[2]; if(!edgeFaces(u,v,ef,op)) continue; // cheap stale check before full
        if(!linkOK(u,v,op)) continue;
        Eval uv=evalEndpoint(u,v,ef,p), vu=evalEndpoint(v,u,ef,p);
        if(p.midpoint){Vec3 mid=(P[u]+P[v])*0.5; Eval um=evalMove(u,v,ef,p,mid), vm=evalMove(v,u,ef,p,mid); if(um.ok&&(!uv.ok||um.cost<uv.cost)) uv=um; if(vm.ok&&(!vu.ok||vm.cost<vu.cost)) vu=vm;}
        if(!uv.ok&&!vu.ok) continue; if(vu.ok&&(!uv.ok||vu.cost<uv.cost)) applyCollapse(v,u,ef,vu,pq); else applyCollapse(u,v,ef,uv,pq);
    }
}

static int countOutputVertices(){int c=0; for(int i=0;i<N;i++) if(aliveV[i]) ++c; return c;}
static void simplify(){
    ShapeStats st=inspectShape();
    double sharp=st.sharp22, verySharp=st.sharp45, smooth=st.smooth10, bad=st.bad;
    double ratio;
    if(N<64) ratio=0.98;
    else if(N<600) ratio=0.30;
    else if(verySharp>0.18||bad>0.03) ratio=0.205;
    else if(sharp>0.08) ratio=0.190;
    else if(smooth>0.985 && sharp<0.012) ratio=0.166;
    else ratio=0.178;
    if(N<5000) ratio=max(ratio,0.176);
    int target=max(8,(int)ceil(N*ratio));
    double A=max(1e-300,diagLen*diagLen*1e-24);
    Params p1{diagLen*0.045,diagLen*0.045,diagLen*0.008,cos(18*M_PI/180.0),diagLen*0.0035,cos(9*M_PI/180.0),A,false};
    Params p2{diagLen*0.0490,diagLen*0.0490,diagLen*0.018,cos(42*M_PI/180.0),diagLen*0.0085,cos(20*M_PI/180.0),A,false};
    Params p3{diagLen*0.0495,diagLen*0.0495,diagLen*0.030,cos(62*M_PI/180.0),diagLen*0.014,cos(32*M_PI/180.0),A,true};
    Params p4{diagLen*0.0497,diagLen*0.0497,diagLen*0.040,cos(72*M_PI/180.0),diagLen*0.021,cos(45*M_PI/180.0),A,true};
    runPass(p1,5.8,target);
    if(elapsed()<12.5 && aliveVertices>target) runPass(p2,12.5,target);
    if(elapsed()<18.0 && aliveVertices>target) runPass(p3,18.0,target);
    if(elapsed()<19.6 && aliveVertices>target && smooth>0.97 && sharp<0.04 && bad<0.02) runPass(p4,19.6,target);
    // round-2 tiny improvement: one extra strictly-flat endpoint-only chase, capped to 1.5% beyond the safe target.
    int chaseTarget=max(8,(int)floor(target*0.985));
    if(elapsed()<20.0 && aliveVertices>chaseTarget && smooth>0.992 && sharp<0.006 && bad<0.01){
        Params pf{diagLen*0.0492,diagLen*0.0492,diagLen*0.016,cos(25*M_PI/180.0),diagLen*0.0045,cos(8*M_PI/180.0),A,false};
        runPass(pf,20.0,chaseTarget);
    }
}

static void outputMesh(){
    vector<int> id(N,-1); vector<Vec3> V; V.reserve(countOutputVertices());
    for(int i=0;i<N;i++) if(aliveV[i]){id[i]=(int)V.size(); V.push_back(P[i]);}
    vector<Face> G; G.reserve(aliveFaces);
    for(int i=0;i<M;i++) if(aliveF[i]){
        int a=id[F[i].v[0]],b=id[F[i].v[1]],c=id[F[i].v[2]]; if(a>=0&&b>=0&&c>=0&&a!=b&&a!=c&&b!=c) G.push_back({{a,b,c}});
    }
    orientLikeOriginal(V,G);
    static char outbuf[1<<20]; setvbuf(stdout,outbuf,_IOFBF,sizeof(outbuf));
    printf("%zu %zu\n",V.size(),G.size());
    for(const auto&p:V) printf("v %.15g %.15g %.15g\n",p.x,p.y,p.z);
    for(const auto&f:G) printf("f %d %d %d\n",f.v[0]+1,f.v[1]+1,f.v[2]+1);
}

int main(){
    T0=chrono::steady_clock::now();
    if(!readInput()) return 0;
    simplify();
    outputMesh();
    return 0;
}
