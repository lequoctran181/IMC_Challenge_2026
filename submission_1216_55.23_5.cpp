#include <bits/stdc++.h>
using namespace std;

struct Vec3{ double x,y,z; };
static inline Vec3 operator+(const Vec3&a,const Vec3&b){return {a.x+b.x,a.y+b.y,a.z+b.z};}
static inline Vec3 operator-(const Vec3&a,const Vec3&b){return {a.x-b.x,a.y-b.y,a.z-b.z};}
static inline Vec3 operator*(const Vec3&a,double s){return {a.x*s,a.y*s,a.z*s};}
static inline Vec3 operator/(const Vec3&a,double s){return {a.x/s,a.y/s,a.z/s};}
static inline double dot3(const Vec3&a,const Vec3&b){return a.x*b.x+a.y*b.y+a.z*b.z;}
static inline Vec3 cross3(const Vec3&a,const Vec3&b){return {a.y*b.z-a.z*b.y,a.z*b.x-a.x*b.z,a.x*b.y-a.y*b.x};}
static inline double norm2(const Vec3&a){return dot3(a,a);} 
static inline double norm3(const Vec3&a){return sqrt(norm2(a));}
static inline Vec3 unit3(Vec3 a){ double n=norm3(a); return n>1e-300 ? a*(1.0/n) : Vec3{0,0,0}; }
static inline double clampd(double x,double a,double b){ return x<a?a:(x>b?b:x); }

struct Face{ int v[3]; unsigned char active; };

struct FastInput{
    vector<char> buf; char* p=nullptr;
    FastInput(){
        buf.reserve(1<<27); char tmp[1<<16]; size_t n;
        while((n=fread(tmp,1,sizeof(tmp),stdin))>0) buf.insert(buf.end(),tmp,tmp+n);
        buf.push_back('\0'); p=buf.data();
    }
    inline void skip(){ while(*p==' '||*p=='\n'||*p=='\r'||*p=='\t') ++p; }
    int nextInt(){ skip(); int s=1; if(*p=='-'){s=-1; ++p;} int x=0; while(*p>='0'&&*p<='9'){x=x*10+(*p-'0'); ++p;} return x*s; }
    double nextDouble(){ skip(); char* e=nullptr; double x=strtod(p,&e); p=e; return x; }
    char nextChar(){ skip(); return *p++; }
};

struct Quadric{
    // symmetric 4x4 plane quadric, stored as xx,xy,xz,xw, yy,yz,yw, zz,zw, ww
    double q[10];
    Quadric(){ memset(q,0,sizeof(q)); }
    void addPlane(double a,double b,double c,double d,double w=1.0){
        q[0]+=w*a*a; q[1]+=w*a*b; q[2]+=w*a*c; q[3]+=w*a*d;
        q[4]+=w*b*b; q[5]+=w*b*c; q[6]+=w*b*d;
        q[7]+=w*c*c; q[8]+=w*c*d; q[9]+=w*d*d;
    }
    void add(const Quadric&o){ for(int i=0;i<10;i++) q[i]+=o.q[i]; }
    double eval(const Vec3&p) const{
        double x=p.x,y=p.y,z=p.z;
        return q[0]*x*x + 2*q[1]*x*y + 2*q[2]*x*z + 2*q[3]*x
             + q[4]*y*y + 2*q[5]*y*z + 2*q[6]*y
             + q[7]*z*z + 2*q[8]*z + q[9];
    }
};

static int N=0,M=0;
static vector<Vec3> origP, P;
static vector<Face> origF, F;
static vector<vector<int>> adj;
static vector<Quadric> Q;
static vector<unsigned char> alive;
static vector<int> headNode, tailNode, nextNode, clusterSize, ver, mark;
static int markToken=10;
static int activeVertices=0, activeFaces=0;
static Vec3 bbMin, bbMax;
static double diagLen=1.0, epsH=0.05, epsH2=0.0025;
static chrono::steady_clock::time_point T0;

static inline double elapsed(){ return chrono::duration<double>(chrono::steady_clock::now()-T0).count(); }
static inline uint64_t edgeKey(int a,int b){ if(a>b) swap(a,b); return (uint64_t)(uint32_t)a<<32 | (uint32_t)b; }
static inline bool hasV(const Face&f,int v){ return f.v[0]==v || f.v[1]==v || f.v[2]==v; }
static inline bool hasEdge(const Face&f,int a,int b){ return hasV(f,a) && hasV(f,b); }
static inline int thirdV(const Face&f,int a,int b){ for(int k=0;k<3;k++){int x=f.v[k]; if(x!=a && x!=b) return x;} return -1; }
static inline Vec3 faceCross(const Face&f){ return cross3(P[f.v[1]]-P[f.v[0]], P[f.v[2]]-P[f.v[0]]); }
static inline double dist2(const Vec3&a,const Vec3&b){ return norm2(a-b); }

static void readInput(){
    FastInput in;
    N=in.nextInt(); M=in.nextInt();
    origP.resize(N); P.resize(N);
    bbMin={1e100,1e100,1e100}; bbMax={-1e100,-1e100,-1e100};
    for(int i=0;i<N;i++){
        (void)in.nextChar();
        Vec3 p{in.nextDouble(),in.nextDouble(),in.nextDouble()};
        origP[i]=P[i]=p;
        bbMin.x=min(bbMin.x,p.x); bbMin.y=min(bbMin.y,p.y); bbMin.z=min(bbMin.z,p.z);
        bbMax.x=max(bbMax.x,p.x); bbMax.y=max(bbMax.y,p.y); bbMax.z=max(bbMax.z,p.z);
    }
    diagLen=norm3(bbMax-bbMin); if(!(diagLen>0)) diagLen=1.0;
    epsH=0.05*diagLen*0.999999; epsH2=epsH*epsH;
    origF.resize(M); F.resize(M);
    vector<int> deg(N,0);
    for(int i=0;i<M;i++){
        (void)in.nextChar();
        int a=in.nextInt()-1,b=in.nextInt()-1,c=in.nextInt()-1;
        origF[i]={{a,b,c},1}; F[i]={{a,b,c},1};
        deg[a]++; deg[b]++; deg[c]++;
    }
    adj.resize(N); for(int i=0;i<N;i++) adj[i].reserve(deg[i]+8);
    for(int i=0;i<M;i++){ adj[F[i].v[0]].push_back(i); adj[F[i].v[1]].push_back(i); adj[F[i].v[2]].push_back(i); }
    Q.assign(N,Quadric()); alive.assign(N,1); activeVertices=N; activeFaces=M;
    headNode.resize(N); tailNode.resize(N); nextNode.assign(N,-1); clusterSize.assign(N,1); ver.assign(N,0); mark.assign(N,0);
    for(int i=0;i<N;i++) headNode[i]=tailNode[i]=i;
}

static void initQuadrics(){
    for(int i=0;i<M;i++){
        const Face&f=F[i];
        Vec3 cr=cross3(P[f.v[1]]-P[f.v[0]],P[f.v[2]]-P[f.v[0]]);
        double l=norm3(cr); if(!(l>1e-300)) continue;
        Vec3 n=cr*(1.0/l);
        double d=-dot3(n,P[f.v[0]]);
        // Area weight is deliberately mild: it favours visual planes but still preserves small details through normal/cluster checks.
        double w=sqrt(l*0.5)+1e-9;
        Q[f.v[0]].addPlane(n.x,n.y,n.z,d,w);
        Q[f.v[1]].addPlane(n.x,n.y,n.z,d,w);
        Q[f.v[2]].addPlane(n.x,n.y,n.z,d,w);
    }
}

static bool solveOptimal(const Quadric&q, Vec3&out){
    double a00=q.q[0],a01=q.q[1],a02=q.q[2];
    double a10=q.q[1],a11=q.q[4],a12=q.q[5];
    double a20=q.q[2],a21=q.q[5],a22=q.q[7];
    double b0=-q.q[3],b1=-q.q[6],b2=-q.q[8];
    double det=a00*(a11*a22-a12*a21)-a01*(a10*a22-a12*a20)+a02*(a10*a21-a11*a20);
    if(fabs(det)<1e-14) return false;
    double dx=b0*(a11*a22-a12*a21)-a01*(b1*a22-a12*b2)+a02*(b1*a21-a11*b2);
    double dy=a00*(b1*a22-a12*b2)-b0*(a10*a22-a12*a20)+a02*(a10*b2-b1*a20);
    double dz=a00*(a11*b2-b1*a21)-a01*(a10*b2-b1*a20)+b0*(a10*a21-a11*a20);
    out={dx/det,dy/det,dz/det};
    return isfinite(out.x)&&isfinite(out.y)&&isfinite(out.z);
}

static bool collectEdgeOpp(int a,int b,int opp[2],int ef[2]){
    int cnt=0;
    const vector<int>& A=adj[a]; const vector<int>& B=adj[b];
    const vector<int>& small=(A.size()<B.size()?A:B);
    for(int fid:small){
        if(!F[fid].active) continue;
        const Face&f=F[fid];
        if(!hasEdge(f,a,b)) continue;
        if(cnt>=2) return false;
        int t=thirdV(f,a,b); if(t<0) return false;
        opp[cnt]=t; ef[cnt]=fid; cnt++;
    }
    return cnt==2 && opp[0]!=opp[1];
}

static bool linkOK(int a,int b,const int opp[2]){
    if(markToken>2000000000){ fill(mark.begin(),mark.end(),0); markToken=10; }
    int tokA=markToken++, tokC=markToken++;
    for(int fid:adj[a]){
        if(!F[fid].active) continue; const Face&f=F[fid]; if(!hasV(f,a)) continue;
        for(int k=0;k<3;k++){ int x=f.v[k]; if(x!=a && x!=b) mark[x]=tokA; }
    }
    int cnt=0, g0=0, g1=0;
    for(int fid:adj[b]){
        if(!F[fid].active) continue; const Face&f=F[fid]; if(!hasV(f,b)) continue;
        for(int k=0;k<3;k++){
            int x=f.v[k]; if(x==a || x==b) continue;
            if(mark[x]==tokA){
                mark[x]=tokC; cnt++;
                if(x==opp[0]) g0=1; if(x==opp[1]) g1=1;
                if(cnt>2) return false;
            }
        }
    }
    return cnt==2 && g0 && g1;
}

static bool clusterCanMove(int v,const Vec3&to){
    for(int m=headNode[v]; m!=-1; m=nextNode[m]) if(dist2(origP[m],to)>epsH2+1e-15) return false;
    return true;
}

struct Params{ double minCos, planeTol, maxEdgeFactor; };

static bool localGeometryOK(int a,int b,const Vec3&to,const Params&par){
    const double minArea2=max(1e-32,diagLen*diagLen*diagLen*diagLen*1e-28);
    auto scan=[&](int root)->bool{
        for(int fid:adj[root]){
            if(!F[fid].active) continue;
            const Face&f=F[fid];
            bool ha=hasV(f,a), hb=hasV(f,b);
            if(!ha && !hb) continue;
            if(ha && hb) continue; // the two edge faces disappear
            Vec3 oldp[3]={P[f.v[0]],P[f.v[1]],P[f.v[2]]};
            Vec3 np[3]={oldp[0],oldp[1],oldp[2]};
            for(int k=0;k<3;k++) if(f.v[k]==a || f.v[k]==b) np[k]=to;
            Vec3 co=cross3(oldp[1]-oldp[0],oldp[2]-oldp[0]);
            Vec3 cn=cross3(np[1]-np[0],np[2]-np[0]);
            double lo=norm2(co), ln=norm2(cn);
            if(lo<=1e-300 || ln<=minArea2) return false;
            double d=dot3(co,cn)/sqrt(lo*ln);
            if(d<par.minCos) return false;
            Vec3 no=co*(1.0/sqrt(lo));
            double pd=1e100;
            for(int k=0;k<3;k++) if(f.v[k]==a || f.v[k]==b) pd=min(pd,fabs(dot3(no,to-oldp[k])));
            if(pd>par.planeTol) return false;
        }
        return true;
    };
    return scan(a) && scan(b);
}

static bool bestPosition(int a,int b,const Params&par,Vec3&best,double&bestCost){
    if(dist2(P[a],P[b]) > par.maxEdgeFactor*par.maxEdgeFactor*epsH2) return false;
    Quadric q=Q[a]; q.add(Q[b]);
    Vec3 cand[8]; int cnt=0; Vec3 opt;
    if(solveOptimal(q,opt)) cand[cnt++]=opt;
    cand[cnt++]=(P[a]+P[b])*0.5;
    cand[cnt++]=P[a]; cand[cnt++]=P[b];
    cand[cnt++]=P[a]*0.75+P[b]*0.25; cand[cnt++]=P[a]*0.25+P[b]*0.75;
    // small tangent-free inward bias candidates are often useful for noisy scans while still bounded by cluster tests
    cand[cnt++]=(P[a]*clusterSize[a]+P[b]*clusterSize[b])*(1.0/(clusterSize[a]+clusterSize[b]));
    bestCost=1e100; bool ok=false;
    for(int i=0;i<cnt;i++){
        const Vec3&p=cand[i];
        if(!isfinite(p.x)||!isfinite(p.y)||!isfinite(p.z)) continue;
        if(!clusterCanMove(a,p) || !clusterCanMove(b,p)) continue;
        if(!localGeometryOK(a,b,p,par)) continue;
        double c=q.eval(p) + 2e-4*(dist2(p,P[a])+dist2(p,P[b])) + 1e-7*(clusterSize[a]+clusterSize[b]);
        if(c<bestCost){ bestCost=c; best=p; ok=true; }
    }
    return ok;
}

struct Node{ double cost; int a,b,va,vb; bool operator<(const Node&o)const{ return cost>o.cost; } };
static priority_queue<Node> pq;

static double cheapCost(int a,int b){
    Quadric q=Q[a]; q.add(Q[b]); Vec3 opt; double best=1e100;
    if(solveOptimal(q,opt) && clusterCanMove(a,opt) && clusterCanMove(b,opt)) best=min(best,q.eval(opt));
    Vec3 m=(P[a]+P[b])*0.5;
    if(clusterCanMove(a,m) && clusterCanMove(b,m)) best=min(best,q.eval(m));
    best=min(best,q.eval(P[a])); best=min(best,q.eval(P[b]));
    return best + 1e-5*dist2(P[a],P[b]);
}

static void pushEdge(int a,int b,const Params&par){
    if(a==b || a<0 || b<0 || a>=N || b>=N || !alive[a] || !alive[b]) return;
    if(dist2(P[a],P[b]) > par.maxEdgeFactor*par.maxEdgeFactor*epsH2) return;
    pq.push({cheapCost(a,b),a,b,ver[a],ver[b]});
}

static void compactAdj(int v){
    vector<int>& x=adj[v];
    if(x.size()<96) return;
    size_t good=0; for(int fid:x) if(F[fid].active && hasV(F[fid],v)) good++;
    if(good*3+64>=x.size()) return;
    vector<int> y; y.reserve(good+8);
    for(int fid:x) if(F[fid].active && hasV(F[fid],v)) y.push_back(fid);
    x.swap(y);
}

static void mergeCluster(int src,int dst){
    if(headNode[src]<0) return;
    nextNode[tailNode[dst]]=headNode[src];
    tailNode[dst]=tailNode[src];
    clusterSize[dst]+=clusterSize[src];
    headNode[src]=tailNode[src]=-1; clusterSize[src]=0;
}

static void doCollapse(int src,int dst,const Vec3&pos,const Params&par){
    Q[dst].add(Q[src]); P[dst]=pos;
    for(int fid:adj[src]){
        if(!F[fid].active) continue;
        Face&f=F[fid];
        bool hs=false, hd=false;
        for(int k=0;k<3;k++){ if(f.v[k]==src) hs=true; if(f.v[k]==dst) hd=true; }
        if(!hs) continue;
        if(hd){ f.active=0; activeFaces--; }
        else{
            for(int k=0;k<3;k++) if(f.v[k]==src) f.v[k]=dst;
            adj[dst].push_back(fid);
        }
    }
    alive[src]=0; activeVertices--; ver[src]++; ver[dst]++;
    mergeCluster(src,dst);
    compactAdj(src); compactAdj(dst);
    for(int fid:adj[dst]){
        if(!F[fid].active) continue; const Face&f=F[fid]; if(!hasV(f,dst)) continue;
        pushEdge(f.v[0],f.v[1],par); pushEdge(f.v[1],f.v[2],par); pushEdge(f.v[2],f.v[0],par);
    }
}

static bool attemptEdge(int a,int b,const Params&par){
    if(a==b || !alive[a] || !alive[b]) return false;
    int opp[2], ef[2]; if(!collectEdgeOpp(a,b,opp,ef)) return false;
    if(!linkOK(a,b,opp)) return false;
    Vec3 pos; double c; if(!bestPosition(a,b,par,pos,c)) return false;
    // Keep the side with larger cluster/adjacency to reduce future adjacency copying.
    int src=a,dst=b;
    long long wa=(long long)adj[a].size()+3LL*clusterSize[a];
    long long wb=(long long)adj[b].size()+3LL*clusterSize[b];
    if(wa>wb){ src=b; dst=a; }
    doCollapse(src,dst,pos,par);
    return true;
}

static void buildInitialQueue(const Params&par){
    vector<uint64_t> keys; keys.reserve((size_t)M*3);
    for(int i=0;i<M;i++) if(F[i].active){
        const Face&f=F[i]; keys.push_back(edgeKey(f.v[0],f.v[1])); keys.push_back(edgeKey(f.v[1],f.v[2])); keys.push_back(edgeKey(f.v[2],f.v[0]));
    }
    sort(keys.begin(),keys.end()); keys.erase(unique(keys.begin(),keys.end()),keys.end());
    priority_queue<Node> empty; pq.swap(empty);
    for(uint64_t k:keys){ int a=(int)(k>>32), b=(int)(uint32_t)k; pushEdge(a,b,par); }
}

static void simplifyTo(int target,double until,const Params&par){
    if(activeVertices<=target) return;
    buildInitialQueue(par);
    long long pops=0, hits=0;
    while(!pq.empty() && activeVertices>target && elapsed()<until){
        Node nd=pq.top(); pq.pop(); pops++;
        int a=nd.a,b=nd.b;
        if(a<0||b<0||a>=N||b>=N||!alive[a]||!alive[b]) continue;
        if(nd.va!=ver[a] || nd.vb!=ver[b]) continue;
        if(attemptEdge(a,b,par)) hits++;
        if((pops&262143)==0 && hits==0 && elapsed()>until-0.2) break;
    }
}

struct MeshOut{ vector<Vec3> V; vector<array<int,3>> T; };
static MeshOut makeOutputFromCurrent(){
    MeshOut out; vector<int> id(N,-1); out.V.reserve(activeVertices);
    for(int i=0;i<N;i++) if(alive[i]){ id[i]=(int)out.V.size(); out.V.push_back(P[i]); }
    out.T.reserve(activeFaces);
    for(int i=0;i<M;i++) if(F[i].active){
        int a=F[i].v[0],b=F[i].v[1],c=F[i].v[2];
        if(a<0||b<0||c<0||a>=N||b>=N||c>=N||!alive[a]||!alive[b]||!alive[c]) continue;
        int ia=id[a],ib=id[b],ic=id[c]; if(ia<0||ib<0||ic<0||ia==ib||ia==ic||ib==ic) continue;
        Vec3 cr=cross3(out.V[ib]-out.V[ia],out.V[ic]-out.V[ia]); if(norm2(cr)<=1e-30) continue;
        out.T.push_back({ia,ib,ic});
    }
    return out;
}
static MeshOut makeOriginalOutput(){
    MeshOut out; out.V=origP; out.T.reserve(M);
    for(const Face&f:origF) out.T.push_back({f.v[0],f.v[1],f.v[2]});
    return out;
}

// Low-resolution visual proxy. It is a guard for aggressive passes, not a replacement for the judge.
struct Image{ int R; vector<float>d,nx,ny,nz; vector<unsigned char>m; };
static inline void viewCoord(int view,const Vec3&p,double&x,double&y,double&z){
    const double D=2.5;
    if(view==0){ z=D-p.x; x=p.y; y=p.z; }       // +X
    else if(view==1){ z=D+p.x; x=-p.y; y=p.z; } // -X
    else if(view==2){ z=D-p.y; x=-p.x; y=p.z; } // +Y
    else if(view==3){ z=D+p.y; x=p.x; y=p.z; }  // -Y
    else if(view==4){ z=D-p.z; x=p.x; y=p.y; }  // +Z
    else { z=D+p.z; x=p.x; y=-p.y; }            // -Z
}
static void renderMesh(const vector<Vec3>&V,const vector<array<int,3>>&T,int R,vector<Image>&imgs,double stopTime=1e100){
    imgs.assign(6,{});
    const float INF=1e30f; double f=800.0*R/1024.0, c=0.5*R;
    for(int vi=0;vi<6;vi++){
        Image im; im.R=R; im.d.assign(R*R,INF); im.nx.assign(R*R,0); im.ny.assign(R*R,0); im.nz.assign(R*R,0); im.m.assign(R*R,0);
        int faceCounter=0;
        for(const auto&t:T){
            if(((++faceCounter)&4095)==0 && elapsed()>stopTime){ imgs.clear(); return; }
            Vec3 p0=V[t[0]],p1=V[t[1]],p2=V[t[2]];
            Vec3 n=unit3(cross3(p1-p0,p2-p0)); if(norm2(n)<=0) continue;
            double x0,y0,z0,x1,y1,z1,x2,y2,z2;
            viewCoord(vi,p0,x0,y0,z0); viewCoord(vi,p1,x1,y1,z1); viewCoord(vi,p2,x2,y2,z2);
            if(z0<=1e-6||z1<=1e-6||z2<=1e-6) continue;
            double u0=f*x0/z0+c, v0=f*y0/z0+c, u1=f*x1/z1+c, v1=f*y1/z1+c, u2=f*x2/z2+c, v2=f*y2/z2+c;
            double minx=floor(min(u0,min(u1,u2))-0.5), maxx=ceil(max(u0,max(u1,u2))+0.5);
            double miny=floor(min(v0,min(v1,v2))-0.5), maxy=ceil(max(v0,max(v1,v2))+0.5);
            int ix0=max(0,(int)minx), ix1=min(R-1,(int)maxx), iy0=max(0,(int)miny), iy1=min(R-1,(int)maxy);
            if(ix0>ix1||iy0>iy1) continue;
            double area=(u1-u0)*(v2-v0)-(v1-v0)*(u2-u0); if(fabs(area)<1e-12) continue;
            double invA=1.0/area;
            for(int yy=iy0; yy<=iy1; yy++){
                double py=yy+0.5;
                for(int xx=ix0; xx<=ix1; xx++){
                    double px=xx+0.5;
                    double w0=((u1-px)*(v2-py)-(v1-py)*(u2-px))*invA;
                    double w1=((u2-px)*(v0-py)-(v2-py)*(u0-px))*invA;
                    double w2=1.0-w0-w1;
                    if(w0>=-1e-9 && w1>=-1e-9 && w2>=-1e-9){
                        double iz=w0/z0+w1/z1+w2/z2; if(iz<=0) continue;
                        float dep=(float)(1.0/iz); int idx=yy*R+xx;
                        if(dep<im.d[idx]){ im.d[idx]=dep; im.nx[idx]=(float)n.x; im.ny[idx]=(float)n.y; im.nz[idx]=(float)n.z; im.m[idx]=1; }
                    }
                }
            }
        }
        imgs[vi] = std::move(im);
    }
}
static double ssimGlobalChannel(const vector<double>&A,const vector<double>&B){
    int n=(int)A.size(); if(n<=2) return 1.0;
    long double sx=0,sy=0; for(int i=0;i<n;i++){sx+=A[i]; sy+=B[i];}
    long double mx=sx/n,my=sy/n,vx=0,vy=0,cov=0;
    for(int i=0;i<n;i++){ long double dx=A[i]-mx, dy=B[i]-my; vx+=dx*dx; vy+=dy*dy; cov+=dx*dy; }
    vx/=n; vy/=n; cov/=n;
    const long double C1=6.5025L, C2=58.5225L;
    long double num=(2*mx*my+C1)*(2*cov+C2);
    long double den=(mx*mx+my*my+C1)*(vx+vy+C2);
    if(den==0) return 1.0;
    double r=(double)(num/den); if(!isfinite(r)) r=0; return clampd(r,0,1);
}
static double visualScoreApprox(const MeshOut&out,int R,double timeLimit){
    if(elapsed()>timeLimit) return 0.0;
    MeshOut orig=makeOriginalOutput();
    vector<Image>A,B; renderMesh(orig.V,orig.T,R,A,timeLimit); if(A.size()!=6||elapsed()>timeLimit) return 0.0; renderMesh(out.V,out.T,R,B,timeLimit); if(B.size()!=6) return 0.0;
    double total=0; int comps=0;
    for(int vi=0;vi<6;vi++){
        vector<double> ax,ay,az,bx,by,bz,ad,bd;
        int RR=A[vi].R, NN=RR*RR; ax.reserve(NN/2);
        for(int i=0;i<NN;i++){
            if(!A[vi].m[i] && !B[vi].m[i]) continue;
            double anx=A[vi].m[i]? (A[vi].nx[i]+1.0)*127.5 : 127.5;
            double any=A[vi].m[i]? (A[vi].ny[i]+1.0)*127.5 : 127.5;
            double anz=A[vi].m[i]? (A[vi].nz[i]+1.0)*127.5 : 127.5;
            double bnx=B[vi].m[i]? (B[vi].nx[i]+1.0)*127.5 : 127.5;
            double bny=B[vi].m[i]? (B[vi].ny[i]+1.0)*127.5 : 127.5;
            double bnz=B[vi].m[i]? (B[vi].nz[i]+1.0)*127.5 : 127.5;
            ax.push_back(anx); ay.push_back(any); az.push_back(anz); bx.push_back(bnx); by.push_back(bny); bz.push_back(bnz);
            double da=A[vi].m[i]? A[vi].d[i]*60.0 : 255.0; double db=B[vi].m[i]? B[vi].d[i]*60.0 : 255.0;
            ad.push_back(clampd(da,0,255)); bd.push_back(clampd(db,0,255));
        }
        if(ax.empty()){ total+=1; comps++; continue; }
        double sn=(ssimGlobalChannel(ax,bx)+ssimGlobalChannel(ay,by)+ssimGlobalChannel(az,bz))/3.0;
        double sd=ssimGlobalChannel(ad,bd);
        total += 0.5*sn + 0.5*sd; comps++;
    }
    return comps?total/comps:1.0;
}

struct SmoothInfo{ double smooth=0, sharp=0; };
static SmoothInfo analyzeSmooth(){
    SmoothInfo s; if(M<10) return s;
    int stride=max(1,M/60000); int tot=0, sm=0, sh=0; double c15=cos(15.0*M_PI/180.0), c55=cos(55.0*M_PI/180.0);
    vector<uint64_t> keys; keys.reserve((size_t)min(M,60000)*3);
    for(int i=0;i<M;i+=stride){ const Face&f=F[i]; keys.push_back(edgeKey(f.v[0],f.v[1])); keys.push_back(edgeKey(f.v[1],f.v[2])); keys.push_back(edgeKey(f.v[2],f.v[0])); }
    sort(keys.begin(),keys.end()); keys.erase(unique(keys.begin(),keys.end()),keys.end());
    for(uint64_t k:keys){ int a=(int)(k>>32), b=(int)(uint32_t)k; int opp[2],ef[2]; if(!collectEdgeOpp(a,b,opp,ef)) continue; Vec3 n0=unit3(faceCross(F[ef[0]])), n1=unit3(faceCross(F[ef[1]])); double d=dot3(n0,n1); tot++; if(d>c15) sm++; if(d<c55) sh++; }
    if(tot){ s.smooth=(double)sm/tot; s.sharp=(double)sh/tot; }
    return s;
}

static bool tryTinySampleBox(MeshOut&out){
    if(N>30) return false;
    vector<Vec3> V={
        {bbMin.x,bbMin.y,bbMin.z},{bbMax.x,bbMin.y,bbMin.z},{bbMax.x,bbMax.y,bbMin.z},{bbMin.x,bbMax.y,bbMin.z},
        {bbMin.x,bbMin.y,bbMax.z},{bbMax.x,bbMin.y,bbMax.z},{bbMax.x,bbMax.y,bbMax.z},{bbMin.x,bbMax.y,bbMax.z}
    };
    for(const Vec3&p:origP){
        double best=1e300;
        for(const Vec3&q:V) best=min(best,dist2(p,q));
        if(best>epsH2+1e-12) return false;
    }
    out.V=V;
    int tri[12][3]={{0,1,2},{0,2,3},{4,6,5},{4,7,6},{0,4,5},{0,5,1},{1,5,6},{1,6,2},{2,6,7},{2,7,3},{3,7,4},{3,4,0}};
    out.T.clear(); out.T.reserve(12);
    Vec3 C=(bbMin+bbMax)*0.5;
    for(auto &t:tri){
        array<int,3> a={t[0],t[1],t[2]};
        Vec3 cr=cross3(out.V[a[1]]-out.V[a[0]],out.V[a[2]]-out.V[a[0]]);
        Vec3 cen=(out.V[a[0]]+out.V[a[1]]+out.V[a[2]])/3.0;
        if(dot3(cr,cen-C)<0) swap(a[1],a[2]);
        out.T.push_back(a);
    }
    return true;
}


struct PointHash{
    double cell=0; Vec3 base; const vector<Vec3>*pts=nullptr; unordered_map<long long, vector<int>> mp;
    static long long key(int x,int y,int z){ return ((long long)(x+1048576)<<42) ^ ((long long)(y+1048576)<<21) ^ (long long)(z+1048576); }
    void build(const vector<Vec3>&P0,double c){
        pts=&P0; cell=max(c,1e-9); base={bbMin.x-2*cell,bbMin.y-2*cell,bbMin.z-2*cell};
        mp.clear(); mp.reserve(P0.size()*2+16);
        for(int i=0;i<(int)P0.size();++i){ int ix=(int)floor((P0[i].x-base.x)/cell), iy=(int)floor((P0[i].y-base.y)/cell), iz=(int)floor((P0[i].z-base.z)/cell); mp[key(ix,iy,iz)].push_back(i); }
    }
    bool nearPoint(const Vec3&q,double r2)const{
        if(!pts) return false; int ix=(int)floor((q.x-base.x)/cell), iy=(int)floor((q.y-base.y)/cell), iz=(int)floor((q.z-base.z)/cell);
        for(int dx=-1;dx<=1;dx++)for(int dy=-1;dy<=1;dy++)for(int dz=-1;dz<=1;dz++){
            auto it=mp.find(key(ix+dx,iy+dy,iz+dz)); if(it==mp.end()) continue;
            for(int id:it->second) if(dist2(q,(*pts)[id])<=r2) return true;
        }
        return false;
    }
};

static bool exactVertexHausdorffOK(const MeshOut&out,double stopTime){
    if(out.V.empty()||out.V.size()>(size_t)N) return false;
    PointHash Horig; Horig.build(origP,epsH);
    for(size_t i=0;i<out.V.size();++i){ if((i&4095)==0 && elapsed()>stopTime) return false; if(!Horig.nearPoint(out.V[i],epsH2+1e-12)) return false; }
    PointHash Hout; Hout.build(out.V,epsH);
    for(int i=0;i<N;++i){ if((i&16383)==0 && elapsed()>stopTime) return false; if(!Hout.nearPoint(origP[i],epsH2+1e-12)) return false; }
    return true;
}

static bool tryRadialCubeShell(MeshOut&out,int S,double stopTime){
    if(N<1200 || S<6 || elapsed()>stopTime) return false;
    unordered_map<long long,int> id; id.reserve((size_t)6*(S+1)*(S+1));
    vector<array<int,3>> coord; vector<Vec3> dir;
    auto ckey=[&](int x,int y,int z)->long long{ int B=2*S+1; return ((long long)(x+S)*B + (y+S))*B + (z+S); };
    auto get=[&](int x,int y,int z)->int{
        long long k=ckey(x,y,z); auto it=id.find(k); if(it!=id.end()) return it->second;
        int q=(int)coord.size(); id[k]=q; coord.push_back({x,y,z}); Vec3 d{(double)x,(double)y,(double)z}; dir.push_back(unit3(d)); return q;
    };
    vector<vector<int>> faceGrid(6, vector<int>((S+1)*(S+1)));
    auto at=[&](int f,int i,int j)->int&{ return faceGrid[f][i*(S+1)+j]; };
    for(int sg=0;sg<2;sg++){
        int sgn=sg?1:-1, f=sg;
        for(int i=0;i<=S;i++)for(int j=0;j<=S;j++) at(f,i,j)=get(sgn*S,-S+2*i,-S+2*j);
    }
    for(int sg=0;sg<2;sg++){
        int sgn=sg?1:-1, f=2+sg;
        for(int i=0;i<=S;i++)for(int j=0;j<=S;j++) at(f,i,j)=get(-S+2*i,sgn*S,-S+2*j);
    }
    for(int sg=0;sg<2;sg++){
        int sgn=sg?1:-1, f=4+sg;
        for(int i=0;i<=S;i++)for(int j=0;j<=S;j++) at(f,i,j)=get(-S+2*i,-S+2*j,sgn*S);
    }
    int G=(int)coord.size(); vector<double> best(G,1e300), rad(G,0.0); vector<unsigned char> have(G,0);
    auto upd=[&](int vid,const Vec3&p){ double r=dot3(p,dir[vid]); if(r<=1e-9) return; Vec3 q=dir[vid]*r; double e=dist2(p,q); if(e<best[vid]){ best[vid]=e; rad[vid]=r; have[vid]=1; } };
    for(int idx=0;idx<N;idx++){
        if((idx&32767)==0 && elapsed()>stopTime) return false;
        Vec3 p=origP[idx]; double ax=fabs(p.x),ay=fabs(p.y),az=fabs(p.z); if(max({ax,ay,az})<1e-12) continue;
        int f; double u,v,den;
        if(ax>=ay&&ax>=az){ f=(p.x>=0)?1:0; den=ax; u=p.y/den; v=p.z/den; }
        else if(ay>=az){ f=(p.y>=0)?3:2; den=ay; u=p.x/den; v=p.z/den; }
        else { f=(p.z>=0)?5:4; den=az; u=p.x/den; v=p.y/den; }
        int ci=(int)floor((u+1)*0.5*S+0.5), cj=(int)floor((v+1)*0.5*S+0.5); ci=max(0,min(S,ci)); cj=max(0,min(S,cj));
        for(int di=-2;di<=2;di++)for(int dj=-2;dj<=2;dj++){ int ii=ci+di,jj=cj+dj; if(ii>=0&&ii<=S&&jj>=0&&jj<=S) upd(at(f,ii,jj),p); }
    }
    vector<vector<int>> nb(G);
    auto addEdge=[&](int a,int b){ if(a==b) return; nb[a].push_back(b); nb[b].push_back(a); };
    out.T.clear(); out.T.reserve((size_t)12*S*S);
    for(int f=0;f<6;f++)for(int i=0;i<S;i++)for(int j=0;j<S;j++){
        int a=at(f,i,j),b=at(f,i+1,j),c=at(f,i+1,j+1),d=at(f,i,j+1);
        addEdge(a,b); addEdge(b,c); addEdge(c,d); addEdge(d,a); addEdge(a,c);
        array<int,3> t1={a,b,c},t2={a,c,d};
        Vec3 cen1=dir[a]+dir[b]+dir[c]; if(dot3(cross3(dir[b]-dir[a],dir[c]-dir[a]),cen1)<0) swap(t1[1],t1[2]);
        Vec3 cen2=dir[a]+dir[c]+dir[d]; if(dot3(cross3(dir[c]-dir[a],dir[d]-dir[a]),cen2)<0) swap(t2[1],t2[2]);
        out.T.push_back(t1); out.T.push_back(t2);
    }
    queue<int>q; for(int i=0;i<G;i++) if(have[i]) q.push(i);
    if(q.empty()) return false;
    while(!q.empty()){
        int v=q.front(); q.pop();
        for(int to:nb[v]) if(!have[to]){ have[to]=1; rad[to]=rad[v]; q.push(to); }
    }
    out.V.resize(G); for(int i=0;i<G;i++) out.V[i]=dir[i]*rad[i];
    for(const auto&t:out.T){ Vec3 cr=cross3(out.V[t[1]]-out.V[t[0]],out.V[t[2]]-out.V[t[0]]); if(norm2(cr)<=1e-24) return false; }
    if(!exactVertexHausdorffOK(out,stopTime)) return false;
    int R=(N>250000?72:(N>60000?96:128));
    double sc=visualScoreApprox(out,R,stopTime);
    return sc>0.982;
}

static bool tryBestRadialShell(MeshOut&out,size_t limitVerts,double stopTime){
    if(elapsed()>stopTime||N<1200) return false;
    vector<int> cand={14,18,24,32,44,60,80,104,132};
    for(int S:cand){
        long long approx=6LL*(S+1)*(S+1)-12LL*(S+1)+8;
        if(approx<=0||approx>(long long)limitVerts||approx>=N) continue;
        MeshOut tmp;
        if(tryRadialCubeShell(tmp,S,stopTime)){ out=std::move(tmp); return true; }
        if(elapsed()>stopTime) break;
    }
    return false;
}

static void saveObj(const MeshOut&out){
    string s; s.reserve(out.V.size()*42 + out.T.size()*24 + 64); char buf[128];
    s.append(buf, snprintf(buf,sizeof(buf), "%d %d\n", (int)out.V.size(), (int)out.T.size()));
    for(const Vec3&p:out.V) s.append(buf, snprintf(buf,sizeof(buf), "v %.10g %.10g %.10g\n", p.x,p.y,p.z));
    for(const auto&t:out.T) s.append(buf, snprintf(buf,sizeof(buf), "f %d %d %d\n", t[0]+1,t[1]+1,t[2]+1));
    fwrite(s.data(),1,s.size(),stdout);
}

int main(){
    ios::sync_with_stdio(false);
    T0=chrono::steady_clock::now();
    readInput();
    MeshOut tiny; if(tryTinySampleBox(tiny)){ saveObj(tiny); return 0; }
    initQuadrics();
    SmoothInfo info=analyzeSmooth();

    // Conservative accepted-core pass. It is endpoint/QEM bounded by exact vertex-cluster Hausdorff and local normal/plane tests.
    Params safe;
    if(info.sharp>0.22){ safe={0.86,0.016*diagLen,1.55}; }
    else if(info.smooth>0.70){ safe={0.62,0.032*diagLen,1.95}; }
    else { safe={0.74,0.024*diagLen,1.75}; }
    double safeRatio = (info.sharp>0.22?0.24:(info.smooth>0.70?0.13:0.18));
    if(N<2000) safeRatio=min(safeRatio,0.70); // tiny cases should not over-churn after sample collapse
    int safeTarget=max(8,(int)ceil(N*safeRatio));
    double safeTime = (N>250000?13.6:12.4);
    simplifyTo(safeTarget,safeTime,safe);
    MeshOut best=makeOutputFromCurrent();

    MeshOut radial;
    if(elapsed()<14.8 && tryBestRadialShell(radial,best.V.size()?best.V.size()-1:0,16.75) && radial.V.size()<best.V.size())
        best=std::move(radial);

    // A second visual-aware speculative pass. Kept only if a low-res renderer proxy says it remains comfortably above 0.9.
    if(elapsed()<16.2 && N>2000){
        Params aggr;
        if(info.sharp>0.22) aggr={0.72,0.026*diagLen,1.85};
        else if(info.smooth>0.70) aggr={0.38,0.050*diagLen,2.20};
        else aggr={0.55,0.038*diagLen,2.05};
        double ratio=(info.sharp>0.22?0.17:(info.smooth>0.70?0.075:0.115));
        int target=max(8,(int)ceil(N*ratio));
        simplifyTo(target,18.15,aggr);
        MeshOut cand=makeOutputFromCurrent();
        bool take = cand.V.size() < best.V.size();
        if(take && elapsed()<18.55){
            int R = (N>250000?72:(N>60000?96:128));
            double sc=visualScoreApprox(cand,R,19.45);
            // The proxy is global/low-res, so require a margin above the official 0.9 threshold.
            if(sc < (info.smooth>0.70?0.955:0.968)) take=false;
        }
        if(take) best = std::move(cand);
    }

    // If all else fails, best is still a manifold produced by valid edge collapses. Output compactly.
    saveObj(best);
    return 0;
}
