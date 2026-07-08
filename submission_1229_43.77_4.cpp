#include <bits/stdc++.h>
using namespace std;

struct Vec3{ double x=0,y=0,z=0; };
static inline Vec3 operator+(const Vec3&a,const Vec3&b){return {a.x+b.x,a.y+b.y,a.z+b.z};}
static inline Vec3 operator-(const Vec3&a,const Vec3&b){return {a.x-b.x,a.y-b.y,a.z-b.z};}
static inline Vec3 operator*(const Vec3&a,double s){return {a.x*s,a.y*s,a.z*s};}
static inline Vec3 operator/(const Vec3&a,double s){return {a.x/s,a.y/s,a.z/s};}
static inline double dot3(const Vec3&a,const Vec3&b){return a.x*b.x+a.y*b.y+a.z*b.z;}
static inline Vec3 cross3(const Vec3&a,const Vec3&b){return {a.y*b.z-a.z*b.y,a.z*b.x-a.x*b.z,a.x*b.y-a.y*b.x};}
static inline double norm2(const Vec3&a){return dot3(a,a);} 
static inline double norm3(const Vec3&a){return sqrt(max(0.0,norm2(a)));}
static inline Vec3 unit3(const Vec3&a){double n=norm3(a); return n>1e-300?a*(1.0/n):Vec3{0,0,0};}
static inline double clampd(double x,double a,double b){return x<a?a:(x>b?b:x);} 
static inline uint64_t ekey(int a,int b){ if(a>b) swap(a,b); return ((uint64_t)(uint32_t)a<<32) | (uint32_t)b; }

struct FastInput{
    vector<char> buf; char* p=nullptr;
    FastInput(){ buf.reserve(1<<27); char tmp[1<<16]; size_t n; while((n=fread(tmp,1,sizeof(tmp),stdin))>0) buf.insert(buf.end(),tmp,tmp+n); buf.push_back('\0'); p=buf.data(); }
    inline void skip(){ while(*p==' '||*p=='\n'||*p=='\r'||*p=='\t') ++p; }
    long nextLong(){ skip(); return strtol(p,&p,10); }
    double nextDouble(){ skip(); return strtod(p,&p); }
    char nextChar(){ skip(); return *p++; }
};

struct Face{ int v[3]; unsigned char alive=1; };
struct Tri{ int a,b,c; };
struct MeshOut{ vector<Vec3> V; vector<Tri> T; };

struct Quadric{
    // xx xy xz xw yy yz yw zz zw ww
    double q[10];
    Quadric(){ memset(q,0,sizeof(q)); }
    inline void addPlane(double a,double b,double c,double d,double w){
        q[0]+=w*a*a; q[1]+=w*a*b; q[2]+=w*a*c; q[3]+=w*a*d;
        q[4]+=w*b*b; q[5]+=w*b*c; q[6]+=w*b*d;
        q[7]+=w*c*c; q[8]+=w*c*d; q[9]+=w*d*d;
    }
    inline void add(const Quadric&o){ for(int i=0;i<10;i++) q[i]+=o.q[i]; }
    inline double eval(const Vec3&p)const{
        double x=p.x,y=p.y,z=p.z;
        return q[0]*x*x+2*q[1]*x*y+2*q[2]*x*z+2*q[3]*x+q[4]*y*y+2*q[5]*y*z+2*q[6]*y+q[7]*z*z+2*q[8]*z+q[9];
    }
};

static int N=0,M=0,aliveVcnt=0,aliveFcnt=0;
static vector<Vec3> origP,P;
static vector<Face> origF,F;
static vector<vector<int>> adj;
static vector<Quadric> Q;
static vector<unsigned char> aliveV, fixedV;
static vector<int> headNode,tailNode,nextNode,clusterSize,ver,markA,markB;
static int stampA=1,stampB=1;
static Vec3 bbMin,bbMax,bbCtr;
static double diagLen=1.0, epsH=0.05, epsH2=0.0025;
static chrono::steady_clock::time_point T0;
static inline double elapsed(){ return chrono::duration<double>(chrono::steady_clock::now()-T0).count(); }
static inline bool hasV(const Face&f,int v){ return f.v[0]==v||f.v[1]==v||f.v[2]==v; }
static inline bool hasEdge(const Face&f,int a,int b){ return hasV(f,a)&&hasV(f,b); }
static inline int thirdV(const Face&f,int a,int b){ for(int k=0;k<3;k++){int x=f.v[k]; if(x!=a&&x!=b) return x;} return -1; }
static inline double dist2(const Vec3&a,const Vec3&b){ return norm2(a-b); }
static inline Vec3 faceCrossCur(const Face&f){ return cross3(P[f.v[1]]-P[f.v[0]],P[f.v[2]]-P[f.v[0]]); }
static inline double coordAxis(const Vec3&p,int ax){ return ax==0?p.x:(ax==1?p.y:p.z); }

static bool readInput(){
    FastInput in; N=(int)in.nextLong(); M=(int)in.nextLong(); if(N<=0||M<=0) return false;
    origP.resize(N); P.resize(N); aliveV.assign(N,1); fixedV.assign(N,0);
    bbMin={1e100,1e100,1e100}; bbMax={-1e100,-1e100,-1e100};
    for(int i=0;i<N;i++){
        (void)in.nextChar(); Vec3 p{in.nextDouble(),in.nextDouble(),in.nextDouble()}; origP[i]=P[i]=p;
        bbMin.x=min(bbMin.x,p.x); bbMin.y=min(bbMin.y,p.y); bbMin.z=min(bbMin.z,p.z);
        bbMax.x=max(bbMax.x,p.x); bbMax.y=max(bbMax.y,p.y); bbMax.z=max(bbMax.z,p.z);
    }
    bbCtr=(bbMin+bbMax)*0.5; diagLen=norm3(bbMax-bbMin); if(!(diagLen>0)) diagLen=1.0; epsH=0.05*diagLen*0.999999; epsH2=epsH*epsH;
    F.resize(M); origF.resize(M); vector<int>deg(N,0);
    for(int i=0;i<M;i++){
        (void)in.nextChar(); int a=(int)in.nextLong()-1,b=(int)in.nextLong()-1,c=(int)in.nextLong()-1;
        F[i]={{a,b,c},1}; origF[i]=F[i]; if((unsigned)a>=N||(unsigned)b>=N||(unsigned)c>=N) return false; deg[a]++; deg[b]++; deg[c]++;
    }
    adj.assign(N,{}); for(int i=0;i<N;i++) adj[i].reserve(deg[i]+8);
    for(int i=0;i<M;i++){ adj[F[i].v[0]].push_back(i); adj[F[i].v[1]].push_back(i); adj[F[i].v[2]].push_back(i); }
    Q.assign(N,Quadric()); headNode.resize(N); tailNode.resize(N); nextNode.assign(N,-1); clusterSize.assign(N,1); ver.assign(N,0); markA.assign(N,0); markB.assign(N,0);
    for(int i=0;i<N;i++) headNode[i]=tailNode[i]=i;
    aliveVcnt=N; aliveFcnt=M;
    for(int ax=0;ax<3;ax++){ int mn=0,mx=0; for(int i=1;i<N;i++){ if(coordAxis(P[i],ax)<coordAxis(P[mn],ax)) mn=i; if(coordAxis(P[i],ax)>coordAxis(P[mx],ax)) mx=i; } fixedV[mn]=fixedV[mx]=1; }
    return true;
}

static void initQuadrics(){
    for(int i=0;i<M;i++){
        const Face&f=F[i]; Vec3 cr=cross3(P[f.v[1]]-P[f.v[0]],P[f.v[2]]-P[f.v[0]]); double l=norm3(cr); if(l<=1e-300) continue;
        Vec3 n=cr*(1.0/l); double d=-dot3(n,P[f.v[0]]); double w=sqrt(max(0.0,l*0.5))+1e-9;
        Q[f.v[0]].addPlane(n.x,n.y,n.z,d,w); Q[f.v[1]].addPlane(n.x,n.y,n.z,d,w); Q[f.v[2]].addPlane(n.x,n.y,n.z,d,w);
    }
}

static bool solveOptimal(const Quadric&q,Vec3&out){
    double a00=q.q[0],a01=q.q[1],a02=q.q[2],a11=q.q[4],a12=q.q[5],a22=q.q[7];
    double b0=-q.q[3],b1=-q.q[6],b2=-q.q[8];
    double det=a00*(a11*a22-a12*a12)-a01*(a01*a22-a12*a02)+a02*(a01*a12-a11*a02);
    if(fabs(det)<1e-15) return false;
    double dx=b0*(a11*a22-a12*a12)-a01*(b1*a22-a12*b2)+a02*(b1*a12-a11*b2);
    double dy=a00*(b1*a22-a12*b2)-b0*(a01*a22-a12*a02)+a02*(a01*b2-b1*a02);
    double dz=a00*(a11*b2-b1*a12)-a01*(a01*b2-b1*a02)+b0*(a01*a12-a11*a02);
    out={dx/det,dy/det,dz/det}; return isfinite(out.x)&&isfinite(out.y)&&isfinite(out.z);
}

static bool collectEdgeOpp(int a,int b,int ef[2],int op[2]){
    const vector<int>& small=(adj[a].size()<adj[b].size()?adj[a]:adj[b]); int cnt=0;
    for(int fid:small){ if(!F[fid].alive) continue; const Face&f=F[fid]; if(!hasEdge(f,a,b)) continue; if(cnt>=2) return false; ef[cnt]=fid; op[cnt]=thirdV(f,a,b); if(op[cnt]<0) return false; cnt++; }
    return cnt==2 && op[0]!=op[1];
}
static bool linkOK(int a,int b,const int op[2]){
    if(++stampA>2000000000){ fill(markA.begin(),markA.end(),0); stampA=1; }
    if(++stampB>2000000000){ fill(markB.begin(),markB.end(),0); stampB=1; }
    for(int fid:adj[a]) if(F[fid].alive&&hasV(F[fid],a)){ const Face&f=F[fid]; for(int k=0;k<3;k++){ int x=f.v[k]; if(x!=a&&x!=b) markA[x]=stampA; } }
    int cnt=0;
    for(int fid:adj[b]) if(F[fid].alive&&hasV(F[fid],b)){ const Face&f=F[fid]; for(int k=0;k<3;k++){ int x=f.v[k]; if(x==a||x==b||markA[x]!=stampA) continue; if(x!=op[0]&&x!=op[1]) return false; if(markB[x]!=stampB){ markB[x]=stampB; cnt++; } } }
    return cnt==2 && markB[op[0]]==stampB && markB[op[1]]==stampB;
}

struct Params{ double minCos, planeTol, maxEdgeFactor; bool allowMove; };
static bool clusterCanMove(int v,const Vec3&to){ for(int m=headNode[v];m!=-1;m=nextNode[m]) if(dist2(origP[m],to)>epsH2+1e-14) return false; return true; }
static bool duplicateFaceAfter(int a,int b,int fid,int x,int y,int z,int e0,int e1){
    array<int,3> key{x,y,z}; sort(key.begin(),key.end()); int best=x; if(adj[y].size()<adj[best].size()) best=y; if(adj[z].size()<adj[best].size()) best=z;
    for(int g:adj[best]){ if(!F[g].alive||g==fid||g==e0||g==e1) continue; if(hasV(F[g],a)||hasV(F[g],b)) continue; array<int,3> h{F[g].v[0],F[g].v[1],F[g].v[2]}; sort(h.begin(),h.end()); if(h==key) return true; }
    return false;
}
static bool localGeometryOK(int a,int b,const int ef[2],const Vec3&to,const Params&par){
    const double minA2=max(1e-32,diagLen*diagLen*diagLen*diagLen*1e-29);
    auto scan=[&](int root)->bool{
        for(int fid:adj[root]){
            if(!F[fid].alive) continue; const Face&f=F[fid]; bool ha=hasV(f,a), hb=hasV(f,b); if(!ha&&!hb) continue; if(ha&&hb) continue;
            int nv[3]={f.v[0],f.v[1],f.v[2]}; for(int k=0;k<3;k++) if(nv[k]==a||nv[k]==b) nv[k]=a;
            if(nv[0]==nv[1]||nv[0]==nv[2]||nv[1]==nv[2]) return false;
            Vec3 oldp[3]={P[f.v[0]],P[f.v[1]],P[f.v[2]]}; Vec3 np[3]={oldp[0],oldp[1],oldp[2]}; for(int k=0;k<3;k++) if(f.v[k]==a||f.v[k]==b) np[k]=to;
            Vec3 co=cross3(oldp[1]-oldp[0],oldp[2]-oldp[0]), cn=cross3(np[1]-np[0],np[2]-np[0]); double ao=norm2(co), an=norm2(cn); if(ao<=1e-300||an<=minA2) return false;
            double cs=dot3(co,cn)/sqrt(ao*an); if(cs<par.minCos) return false;
            Vec3 no=co*(1.0/sqrt(ao)); double pd=0.0; for(int k=0;k<3;k++) if(f.v[k]==a||f.v[k]==b) pd=max(pd,fabs(dot3(no,to-oldp[k]))); if(pd>par.planeTol) return false;
            if(duplicateFaceAfter(a,b,fid,nv[0],nv[1],nv[2],ef[0],ef[1])) return false;
        }
        return true;
    };
    return scan(a)&&scan(b);
}
static bool bestPosition(int a,int b,const int ef[2],const Params&par,Vec3&best,double&bestCost){
    if(dist2(P[a],P[b])>par.maxEdgeFactor*par.maxEdgeFactor*epsH2) return false;
    Quadric q=Q[a]; q.add(Q[b]); Vec3 cand[10]; int cnt=0,optok=0; Vec3 opt;
    if(par.allowMove && solveOptimal(q,opt)){ cand[cnt++]=opt; optok=1; }
    if(par.allowMove){ cand[cnt++]=(P[a]+P[b])*0.5; cand[cnt++]=P[a]*0.75+P[b]*0.25; cand[cnt++]=P[a]*0.25+P[b]*0.75; cand[cnt++]=(P[a]*(double)clusterSize[a]+P[b]*(double)clusterSize[b])/(double)(clusterSize[a]+clusterSize[b]); }
    cand[cnt++]=P[a]; cand[cnt++]=P[b]; bestCost=1e100; bool ok=false;
    for(int i=0;i<cnt;i++){
        Vec3 p=cand[i]; if(!isfinite(p.x)||!isfinite(p.y)||!isfinite(p.z)) continue;
        if(!clusterCanMove(a,p)||!clusterCanMove(b,p)) continue;
        if(!localGeometryOK(a,b,ef,p,par)) continue;
        double c=q.eval(p)+1e-4*(dist2(p,P[a])+dist2(p,P[b]))+1e-7*(clusterSize[a]+clusterSize[b]);
        if(i==0&&optok) c*=0.85; if(c<bestCost){bestCost=c; best=p; ok=true;}
    }
    return ok;
}

struct Node{ double c; int a,b,va,vb; bool operator<(const Node&o)const{return c>o.c;} };
static priority_queue<Node> pq;
static double cheapCost(int a,int b){ Quadric q=Q[a]; q.add(Q[b]); double c=min(q.eval(P[a]),q.eval(P[b])); Vec3 m=(P[a]+P[b])*0.5; c=min(c,q.eval(m)); return c+1e-6*dist2(P[a],P[b])+1e-9*(clusterSize[a]+clusterSize[b]); }
static void pushEdge(int a,int b,const Params&par){ if(a==b||a<0||b<0||a>=N||b>=N||!aliveV[a]||!aliveV[b]) return; if(dist2(P[a],P[b])>par.maxEdgeFactor*par.maxEdgeFactor*epsH2) return; pq.push({cheapCost(a,b),a,b,ver[a],ver[b]}); }
static void compactAdj(int v,bool force=false){ vector<int>&x=adj[v]; if(!force&&x.size()<160) return; size_t good=0; for(int fid:x) if(F[fid].alive&&hasV(F[fid],v)) good++; if(!force&&good*3+64>=x.size()) return; vector<int> y; y.reserve(good+8); for(int fid:x) if(F[fid].alive&&hasV(F[fid],v)) y.push_back(fid); x.swap(y); }
static void mergeCluster(int src,int dst){ if(headNode[src]<0) return; nextNode[tailNode[dst]]=headNode[src]; tailNode[dst]=tailNode[src]; clusterSize[dst]+=clusterSize[src]; headNode[src]=tailNode[src]=-1; clusterSize[src]=0; }
static void doCollapse(int src,int dst,const int ef[2],const Vec3&pos,const Params&par){
    Q[dst].add(Q[src]); P[dst]=pos;
    for(int k=0;k<2;k++) if(F[ef[k]].alive){ F[ef[k]].alive=0; aliveFcnt--; }
    for(int fid:adj[src]){
        if(!F[fid].alive) continue; Face&f=F[fid]; bool hs=false; for(int k=0;k<3;k++) if(f.v[k]==src) hs=true; if(!hs) continue;
        for(int k=0;k<3;k++) if(f.v[k]==src) f.v[k]=dst; adj[dst].push_back(fid);
    }
    aliveV[src]=0; aliveVcnt--; ver[src]++; ver[dst]++; mergeCluster(src,dst); compactAdj(src,true); compactAdj(dst);
    for(int fid:adj[dst]) if(F[fid].alive&&hasV(F[fid],dst)){ const Face&f=F[fid]; pushEdge(f.v[0],f.v[1],par); pushEdge(f.v[1],f.v[2],par); pushEdge(f.v[2],f.v[0],par); }
}
static bool attemptEdge(int a,int b,const Params&par){
    if(a==b||!aliveV[a]||!aliveV[b]) return false; if(fixedV[a]&&fixedV[b]) return false;
    int ef[2],op[2]; if(!collectEdgeOpp(a,b,ef,op)) return false; if(!linkOK(a,b,op)) return false;
    Vec3 pos; double cost; if(!bestPosition(a,b,ef,par,pos,cost)) return false;
    int src=a,dst=b; long long wa=(long long)adj[a].size()+3LL*clusterSize[a]+(fixedV[a]?1000000000LL:0); long long wb=(long long)adj[b].size()+3LL*clusterSize[b]+(fixedV[b]?1000000000LL:0);
    if(wa>wb){ src=b; dst=a; }
    // if a fixed vertex participates, keep it when possible
    if(fixedV[a]&&!fixedV[b]){src=b;dst=a;} else if(fixedV[b]&&!fixedV[a]){src=a;dst=b;}
    doCollapse(src,dst,ef,pos,par); return true;
}
static void buildQueue(const Params&par){
    vector<uint64_t> keys; keys.reserve((size_t)aliveFcnt*3);
    for(int i=0;i<M;i++) if(F[i].alive){ const Face&f=F[i]; keys.push_back(ekey(f.v[0],f.v[1])); keys.push_back(ekey(f.v[1],f.v[2])); keys.push_back(ekey(f.v[2],f.v[0])); }
    sort(keys.begin(),keys.end()); keys.erase(unique(keys.begin(),keys.end()),keys.end()); priority_queue<Node> empty; pq.swap(empty); for(uint64_t k:keys){ pushEdge((int)(k>>32),(int)(uint32_t)k,par); }
}
static void simplifyTo(int target,double timeLimit,const Params&par){
    if(aliveVcnt<=target||elapsed()>timeLimit) return; buildQueue(par); long long pops=0;
    while(!pq.empty()&&aliveVcnt>target){
        if((++pops&8191)==0&&elapsed()>timeLimit) break; Node nd=pq.top(); pq.pop(); int a=nd.a,b=nd.b;
        if(a<0||b<0||a>=N||b>=N||!aliveV[a]||!aliveV[b]||nd.va!=ver[a]||nd.vb!=ver[b]) continue;
        attemptEdge(a,b,par);
    }
    for(int i=0;i<N;i++) if(aliveV[i]&&adj[i].size()>512) compactAdj(i,true);
}

struct SmoothInfo{ double smooth10=0,smooth30=0,sharp45=0,bad=0; };
static SmoothInfo analyzeSmooth(){
    SmoothInfo s; int stride=max(1,M/60000),tot=0,sm10=0,sm30=0,sh45=0,bad=0; double c10=cos(10*M_PI/180.0),c30=cos(30*M_PI/180.0),c45=cos(45*M_PI/180.0);
    vector<uint64_t> keys; keys.reserve((size_t)min(M,60000)*3);
    for(int i=0;i<M;i+=stride){ const Face&f=F[i]; keys.push_back(ekey(f.v[0],f.v[1])); keys.push_back(ekey(f.v[1],f.v[2])); keys.push_back(ekey(f.v[2],f.v[0])); }
    sort(keys.begin(),keys.end()); keys.erase(unique(keys.begin(),keys.end()),keys.end());
    for(uint64_t k:keys){ int a=(int)(k>>32),b=(int)(uint32_t)k,ef[2],op[2]; if(!collectEdgeOpp(a,b,ef,op)){bad++; continue;} Vec3 n0=unit3(faceCrossCur(F[ef[0]])), n1=unit3(faceCrossCur(F[ef[1]])); double d=clampd(dot3(n0,n1),-1,1); tot++; if(d>c10) sm10++; if(d>c30) sm30++; if(d<c45) sh45++; }
    if(tot){s.smooth10=(double)sm10/tot; s.smooth30=(double)sm30/tot; s.sharp45=(double)sh45/tot; s.bad=(double)bad/(tot+bad);} return s;
}

static MeshOut makeCurrent(){ MeshOut o; vector<int> id(N,-1); o.V.reserve(aliveVcnt); for(int i=0;i<N;i++) if(aliveV[i]){id[i]=(int)o.V.size(); o.V.push_back(P[i]);} o.T.reserve(aliveFcnt); for(int i=0;i<M;i++) if(F[i].alive){int a=id[F[i].v[0]],b=id[F[i].v[1]],c=id[F[i].v[2]]; if(a>=0&&b>=0&&c>=0&&a!=b&&a!=c&&b!=c&&norm2(cross3(o.V[b]-o.V[a],o.V[c]-o.V[a]))>1e-30) o.T.push_back({a,b,c});} return o; }
static MeshOut makeOriginal(){ MeshOut o; o.V=origP; o.T.reserve(M); for(const Face&f:origF)o.T.push_back({f.v[0],f.v[1],f.v[2]}); return o; }

static bool validateOut(const MeshOut&o){
    if(o.V.empty()||o.V.size()>(size_t)N||o.T.size()<4) return false; vector<uint64_t> edges; edges.reserve(o.T.size()*3);
    double epsA=max(1e-32,diagLen*diagLen*diagLen*diagLen*1e-30);
    vector<int> used(o.V.size(),0);
    for(const Tri&t:o.T){ if((unsigned)t.a>=o.V.size()||(unsigned)t.b>=o.V.size()||(unsigned)t.c>=o.V.size()) return false; if(t.a==t.b||t.a==t.c||t.b==t.c) return false; if(norm2(cross3(o.V[t.b]-o.V[t.a],o.V[t.c]-o.V[t.a]))<=epsA) return false; edges.push_back(ekey(t.a,t.b)); edges.push_back(ekey(t.b,t.c)); edges.push_back(ekey(t.c,t.a)); used[t.a]=used[t.b]=used[t.c]=1; }
    sort(edges.begin(),edges.end()); for(size_t i=0;i<edges.size();){ size_t j=i+1; while(j<edges.size()&&edges[j]==edges[i]) j++; if(j-i!=2) return false; i=j; }
    return true;
}

struct Image{ int R=0; vector<float>d,nx,ny,nz; vector<unsigned char>m; };
static inline void viewCoord(int view,const Vec3&p,double&x,double&y,double&z){ const double D=2.5; if(view==0){z=D-p.x;x=p.y;y=p.z;} else if(view==1){z=D+p.x;x=-p.y;y=p.z;} else if(view==2){z=D-p.y;x=-p.x;y=p.z;} else if(view==3){z=D+p.y;x=p.x;y=p.z;} else if(view==4){z=D-p.z;x=p.x;y=p.y;} else {z=D+p.z;x=p.x;y=-p.y;} }
static void renderMesh(const vector<Vec3>&V,const vector<Tri>&T,int R,vector<Image>&imgs,double stopT){
    imgs.assign(6,{}); double f=800.0*R/1024.0,c=0.5*R,depthScale=60.0;
    for(int vi=0;vi<6;vi++){
        Image im; im.R=R; int S=R*R; im.d.assign(S,255.0f); im.nx.assign(S,127.5f); im.ny.assign(S,127.5f); im.nz.assign(S,127.5f); im.m.assign(S,0);
        int fc=0; for(const Tri&t:T){ if(((++fc)&4095)==0&&elapsed()>stopT){imgs.clear();return;} const Vec3&p0=V[t.a],&p1=V[t.b],&p2=V[t.c]; Vec3 n=unit3(cross3(p1-p0,p2-p0)); if(norm2(n)<=0) continue; double x0,y0,z0,x1,y1,z1,x2,y2,z2; viewCoord(vi,p0,x0,y0,z0); viewCoord(vi,p1,x1,y1,z1); viewCoord(vi,p2,x2,y2,z2); if(z0<=1e-8||z1<=1e-8||z2<=1e-8) continue; double u0=f*x0/z0+c,v0=f*y0/z0+c,u1=f*x1/z1+c,v1=f*y1/z1+c,u2=f*x2/z2+c,v2=f*y2/z2+c; int ix0=max(0,(int)floor(min({u0,u1,u2})-1)), ix1=min(R-1,(int)ceil(max({u0,u1,u2})+1)); int iy0=max(0,(int)floor(min({v0,v1,v2})-1)), iy1=min(R-1,(int)ceil(max({v0,v1,v2})+1)); if(ix0>ix1||iy0>iy1) continue; double den=(v1-v2)*(u0-u2)+(u2-u1)*(v0-v2); if(fabs(den)<1e-14) continue; float rx=(float)((n.x+1.0)*127.5), ry=(float)((n.y+1.0)*127.5), rz=(float)((n.z+1.0)*127.5); for(int yy=iy0;yy<=iy1;yy++){ double py=yy+0.5; for(int xx=ix0;xx<=ix1;xx++){ double px=xx+0.5; double w0=((v1-v2)*(px-u2)+(u2-u1)*(py-v2))/den; double w1=((v2-v0)*(px-u2)+(u0-u2)*(py-v2))/den; double w2=1.0-w0-w1; if(w0>=-1e-9&&w1>=-1e-9&&w2>=-1e-9){ double iz=w0/z0+w1/z1+w2/z2; if(iz<=0) continue; float dep=(float)clampd((1.0/iz)*depthScale,0,255); int id=yy*R+xx; if(!im.m[id]||dep<im.d[id]){im.d[id]=dep; im.nx[id]=rx; im.ny[id]=ry; im.nz[id]=rz; im.m[id]=1;}} }} }
        imgs[vi]=std::move(im);
    }
}
static double ssimLocal(const vector<float>&A,const vector<float>&B,const vector<unsigned char>&fg,int R,int win){ int rad=win/2; const double C1=6.5025,C2=58.5225; double total=0; int cnt=0; for(int y=0;y<R;y++)for(int x=0;x<R;x++){int id=y*R+x; if(!fg[id]) continue; double sx=0,sy=0,sxx=0,syy=0,sxy=0; int n=0; for(int dy=-rad;dy<=rad;dy++){int yy=min(R-1,max(0,y+dy)); int row=yy*R; for(int dx=-rad;dx<=rad;dx++){int xx=min(R-1,max(0,x+dx)); int k=row+xx; double vx=A[k],vy=B[k]; sx+=vx; sy+=vy; sxx+=vx*vx; syy+=vy*vy; sxy+=vx*vy; n++;}} double inv=1.0/n,mx=sx*inv,my=sy*inv,vx=max(0.0,sxx*inv-mx*mx),vy=max(0.0,syy*inv-my*my),cov=sxy*inv-mx*my; double den=(mx*mx+my*my+C1)*(vx+vy+C2); double r=den?((2*mx*my+C1)*(2*cov+C2))/den:1.0; if(isfinite(r)) total+=clampd(r,-1,1); cnt++; } return cnt?total/cnt:1.0; }
static double visualProxy(const MeshOut&o,int R,double stopT){
    if(elapsed()>stopT) return 0; MeshOut org=makeOriginal(); vector<Image>A,B; renderMesh(org.V,org.T,R,A,stopT); if(A.size()!=6||elapsed()>stopT) return 0; renderMesh(o.V,o.T,R,B,stopT); if(B.size()!=6) return 0; double tot=0; for(int v=0;v<6;v++){ int S=R*R; vector<unsigned char> fg(S); for(int i=0;i<S;i++) fg[i]=A[v].m[i]||B[v].m[i]; double ns=(ssimLocal(A[v].nx,B[v].nx,fg,R,11)+ssimLocal(A[v].ny,B[v].ny,fg,R,11)+ssimLocal(A[v].nz,B[v].nz,fg,R,11))/3.0; double ds=ssimLocal(A[v].d,B[v].d,fg,R,11); tot += 0.5*ns+0.5*ds; if(elapsed()>stopT) return 0; } return tot/6.0;
}

static bool trySampleBox(MeshOut&o){
    if(N>30) return false; vector<Vec3> V={{bbMin.x,bbMin.y,bbMin.z},{bbMax.x,bbMin.y,bbMin.z},{bbMax.x,bbMax.y,bbMin.z},{bbMin.x,bbMax.y,bbMin.z},{bbMin.x,bbMin.y,bbMax.z},{bbMax.x,bbMin.y,bbMax.z},{bbMax.x,bbMax.y,bbMax.z},{bbMin.x,bbMax.y,bbMax.z}}; for(const Vec3&p:origP){ double bd=1e100; for(const Vec3&q:V) bd=min(bd,dist2(p,q)); if(bd>epsH2+1e-12) return false; } int tr[12][3]={{0,2,1},{0,3,2},{4,5,6},{4,6,7},{0,1,5},{0,5,4},{1,2,6},{1,6,5},{2,3,7},{2,7,6},{3,0,4},{3,4,7}}; o.V=V; o.T.clear(); for(auto&t:tr)o.T.push_back({t[0],t[1],t[2]}); return validateOut(o); }

static void saveMesh(const MeshOut&o){ static char obuf[1<<20]; setvbuf(stdout,obuf,_IOFBF,sizeof(obuf)); printf("%zu %zu\n",o.V.size(),o.T.size()); for(const Vec3&p:o.V) printf("v %.15g %.15g %.15g\n",p.x,p.y,p.z); for(const Tri&t:o.T) printf("f %d %d %d\n",t.a+1,t.b+1,t.c+1); }

int main(){
    T0=chrono::steady_clock::now(); if(!readInput()) return 0; MeshOut box; if(trySampleBox(box)){saveMesh(box);return 0;} initQuadrics(); SmoothInfo si=analyzeSmooth();
    double smooth=si.smooth10, smooth30=si.smooth30, sharp=si.sharp45;
    Params p1; double r1;
    if(sharp>0.20){ p1={0.88,0.014*diagLen,1.50,false}; r1=0.245; }
    else if(smooth>0.985&&sharp<0.015){ p1={0.52,0.036*diagLen,2.05,true}; r1=0.118; }
    else if(smooth30>0.86){ p1={0.62,0.030*diagLen,1.92,true}; r1=0.145; }
    else { p1={0.74,0.022*diagLen,1.72,true}; r1=0.185; }
    if(N<2500) r1=max(r1,0.34); else if(N<8000) r1=max(r1,0.19);
    int t1=max(8,(int)ceil(N*r1)); double T1=N>250000?9.7:8.7; simplifyTo(t1,T1,p1);
    MeshOut best=makeCurrent(); if(!validateOut(best)) best=makeOriginal();
    // The first pass is deliberately conservative and is always kept as fallback; the proxy only guards the second pass.
    // Second, more speculative QEM pass. It is only used with an official-style local SSIM proxy margin.
    if(elapsed()<15.0 && N>1800){
        Params p2; double r2; if(sharp>0.20){ p2={0.78,0.022*diagLen,1.74,true}; r2=0.17; }
        else if(smooth>0.985&&sharp<0.015){ p2={0.26,0.058*diagLen,2.35,true}; r2=0.062; }
        else if(smooth30>0.86){ p2={0.42,0.046*diagLen,2.18,true}; r2=0.087; }
        else { p2={0.56,0.036*diagLen,1.98,true}; r2=0.118; }
        if(N<2500) r2=max(r2,0.14); else if(N<8000) r2=max(r2,0.105); int t2=max(8,(int)ceil(N*r2)); simplifyTo(t2,17.5,p2); MeshOut cand=makeCurrent();
        if(cand.V.size()<best.V.size() && validateOut(cand)){ bool take=false; if(N>300000 && sharp<0.10 && smooth30>0.82) take=true; else if(elapsed()<18.4){ int R=N>260000?80:(N>70000?96:128); double vp=visualProxy(cand,R,20.55); double th=(smooth>0.985&&sharp<0.015)?(N<10000?0.900:0.925):(sharp>0.20?0.955:(N<10000?0.922:0.940)); if(vp>=th) take=true; } if(take) best=std::move(cand); }
    }
    if(!validateOut(best)) best=makeOriginal(); saveMesh(best); return 0;
}
