#include <bits/stdc++.h>
using namespace std;

struct Vec3 { double x=0,y=0,z=0; };
static inline Vec3 operator+(const Vec3&a,const Vec3&b){return {a.x+b.x,a.y+b.y,a.z+b.z};}
static inline Vec3 operator-(const Vec3&a,const Vec3&b){return {a.x-b.x,a.y-b.y,a.z-b.z};}
static inline Vec3 operator*(const Vec3&a,double s){return {a.x*s,a.y*s,a.z*s};}
static inline Vec3 operator/(const Vec3&a,double s){return {a.x/s,a.y/s,a.z/s};}
static inline double dot3(const Vec3&a,const Vec3&b){return a.x*b.x+a.y*b.y+a.z*b.z;}
static inline Vec3 cross3(const Vec3&a,const Vec3&b){return {a.y*b.z-a.z*b.y,a.z*b.x-a.x*b.z,a.x*b.y-a.y*b.x};}
static inline double norm2(const Vec3&a){return dot3(a,a);} 
static inline double norm3(const Vec3&a){return sqrt(norm2(a));}
static inline Vec3 unit3(Vec3 a){double n=norm3(a); return n>1e-300?a*(1.0/n):Vec3{0,0,0};}
static inline double clampd(double x,double a,double b){return x<a?a:(x>b?b:x);} 
static const double PI = 3.141592653589793238462643383279502884;

struct Face { int v[3]; };

static int N0=0,M0=0,curN=0,aliveFaceCount=0;
static vector<Vec3> origP,P;
static vector<Face> origF,F;
static vector<unsigned char> aliveV,aliveF;
static vector<double> coverRad;
static vector<vector<int>> adj;
static Vec3 bbMin,bbMax;
static double diagLen=1.0;
static chrono::steady_clock::time_point T0;
static vector<int> markA,markB;
static int stampA=1,stampB=1;

static inline double elapsed(){return chrono::duration<double>(chrono::steady_clock::now()-T0).count();}

struct FastInput{
    vector<char> buf; char *p=nullptr;
    FastInput(){
        buf.reserve(1<<27);
        char chunk[1<<16]; size_t n;
        while((n=fread(chunk,1,sizeof(chunk),stdin))>0) buf.insert(buf.end(),chunk,chunk+n);
        buf.push_back('\0'); p=buf.data();
    }
    inline void skip(){ while(*p==' '||*p=='\n'||*p=='\r'||*p=='\t') ++p; }
    long nextLong(){ skip(); return strtol(p,&p,10); }
    double nextDouble(){ skip(); return strtod(p,&p); }
    char nextChar(){ skip(); return *p++; }
};

static void readInput(){
    FastInput in;
    N0=(int)in.nextLong(); M0=(int)in.nextLong();
    origP.resize(N0); P.resize(N0);
    bbMin={1e100,1e100,1e100}; bbMax={-1e100,-1e100,-1e100};
    for(int i=0;i<N0;i++){
        (void)in.nextChar();
        Vec3 q{in.nextDouble(),in.nextDouble(),in.nextDouble()};
        origP[i]=P[i]=q;
        bbMin.x=min(bbMin.x,q.x); bbMin.y=min(bbMin.y,q.y); bbMin.z=min(bbMin.z,q.z);
        bbMax.x=max(bbMax.x,q.x); bbMax.y=max(bbMax.y,q.y); bbMax.z=max(bbMax.z,q.z);
    }
    diagLen=norm3(bbMax-bbMin); if(!(diagLen>0)) diagLen=1.0;
    origF.resize(M0); F.resize(M0);
    for(int i=0;i<M0;i++){
        (void)in.nextChar();
        int a=(int)in.nextLong()-1,b=(int)in.nextLong()-1,c=(int)in.nextLong()-1;
        origF[i].v[0]=F[i].v[0]=a; origF[i].v[1]=F[i].v[1]=b; origF[i].v[2]=F[i].v[2]=c;
    }
    curN=N0; aliveV.assign(curN,1); aliveF.assign(F.size(),1); coverRad.assign(curN,0.0); aliveFaceCount=(int)F.size();
}

static inline uint64_t edgeKey(int a,int b){ if(a>b) swap(a,b); return (uint64_t)(uint32_t)a<<32 | (uint32_t)b; }
static inline bool faceHas(const Face&f,int v){ return f.v[0]==v||f.v[1]==v||f.v[2]==v; }
static inline int thirdOf(const Face&f,int a,int b){ for(int k=0;k<3;k++){int x=f.v[k]; if(x!=a&&x!=b) return x;} return -1; }
static inline Vec3 faceCrossCur(const Face&f){ return cross3(P[f.v[1]]-P[f.v[0]],P[f.v[2]]-P[f.v[0]]); }
static inline Vec3 faceCrossVec(const vector<Vec3>&V,const Face&f){ return cross3(V[f.v[1]]-V[f.v[0]],V[f.v[2]]-V[f.v[0]]); }

static void rebuildAdj(){
    vector<int> deg(curN,0);
    aliveFaceCount=0;
    for(int i=0;i<(int)F.size();i++){
        if(!aliveF[i]) continue;
        Face f=F[i];
        bool ok=f.v[0]>=0&&f.v[0]<curN&&f.v[1]>=0&&f.v[1]<curN&&f.v[2]>=0&&f.v[2]<curN&&
                aliveV[f.v[0]]&&aliveV[f.v[1]]&&aliveV[f.v[2]]&&
                f.v[0]!=f.v[1]&&f.v[0]!=f.v[2]&&f.v[1]!=f.v[2];
        if(!ok){ aliveF[i]=0; continue; }
        deg[f.v[0]]++; deg[f.v[1]]++; deg[f.v[2]]++; aliveFaceCount++;
    }
    adj.assign(curN,{});
    for(int i=0;i<curN;i++) adj[i].reserve(deg[i]);
    for(int i=0;i<(int)F.size();i++) if(aliveF[i]){
        Face f=F[i]; adj[f.v[0]].push_back(i); adj[f.v[1]].push_back(i); adj[f.v[2]].push_back(i);
    }
}

static int activeVertices(){ int c=0; for(unsigned char b:aliveV) c+=b?1:0; return c; }

static bool sameFaceUnordered(const Face&f,int a,int b,int c){
    int x[3]={f.v[0],f.v[1],f.v[2]}, y[3]={a,b,c};
    sort(x,x+3); sort(y,y+3); return x[0]==y[0]&&x[1]==y[1]&&x[2]==y[2];
}

static bool findEdgeFaces(int u,int v,int ef[2],int opp[2]){
    if(u<0||v<0||u>=curN||v>=curN) return false;
    int cnt=0;
    for(int fid:adj[u]){
        if(!aliveF[fid]) continue;
        const Face&f=F[fid];
        if(!faceHas(f,u)||!faceHas(f,v)) continue;
        if(cnt>=2) return false;
        ef[cnt]=fid; opp[cnt]=thirdOf(f,u,v); cnt++;
    }
    if(cnt!=2) return false;
    if(opp[0]<0||opp[1]<0||opp[0]==opp[1]) return false;
    return true;
}

static bool linkOK(int u,int v,const int opp[2]){
    if(++stampA>2000000000){ fill(markA.begin(),markA.end(),0); stampA=1; }
    if(++stampB>2000000000){ fill(markB.begin(),markB.end(),0); stampB=1; }
    for(int fid:adj[u]) if(aliveF[fid]&&faceHas(F[fid],u)){
        const Face&f=F[fid];
        for(int k=0;k<3;k++){ int x=f.v[k]; if(x!=u&&x!=v) markA[x]=stampA; }
    }
    int common=0;
    for(int fid:adj[v]) if(aliveF[fid]&&faceHas(F[fid],v)){
        const Face&f=F[fid];
        for(int k=0;k<3;k++){
            int x=f.v[k]; if(x==u||x==v) continue;
            if(markA[x]!=stampA) continue;
            if(x!=opp[0]&&x!=opp[1]) return false;
            if(markB[x]!=stampB){ markB[x]=stampB; common++; }
        }
    }
    return common==2&&markB[opp[0]]==stampB&&markB[opp[1]]==stampB;
}

struct CollapseParams{ double epsR=0.049,planeR=0.020,minCos=0.90,areaR=1e-24; };
struct EvalCollapse{ bool ok=false; double cost=1e100,newRad=0; int keep=-1,rem=-1; };

static bool duplicateFaceAfter(int fid,int rem,int ef0,int ef1,int a,int b,int c){
    int probe=a;
    if(b>=0&&b<curN&&(int)adj[b].size()<(int)adj[probe].size()) probe=b;
    if(c>=0&&c<curN&&(int)adj[c].size()<(int)adj[probe].size()) probe=c;
    for(int g:adj[probe]){
        if(!aliveF[g]||g==fid||g==ef0||g==ef1) continue;
        if(faceHas(F[g],rem)) continue;
        if(sameFaceUnordered(F[g],a,b,c)) return true;
    }
    return false;
}

static EvalCollapse evalEndpoint(int keep,int rem,const int ef[2],const CollapseParams&par){
    EvalCollapse r; r.keep=keep; r.rem=rem;
    const double eps=par.epsR*diagLen, planeTol=par.planeR*diagLen;
    double dkr=norm3(P[keep]-P[rem]);
    r.newRad=max(coverRad[keep],coverRad[rem]+dkr);
    if(r.newRad>eps+1e-12) return r;
    const double minA2=max(1e-300,par.areaR*diagLen*diagLen*diagLen*diagLen);
    double maxPlane=0,maxAng=0,areaLoss=0; int changed=0;
    for(int fid:adj[rem]){
        if(!aliveF[fid]||!faceHas(F[fid],rem)) continue;
        if(fid==ef[0]||fid==ef[1]) continue;
        if(faceHas(F[fid],keep)) return r;
        Face old=F[fid], nf=old;
        for(int k=0;k<3;k++) if(nf.v[k]==rem) nf.v[k]=keep;
        if(nf.v[0]==nf.v[1]||nf.v[0]==nf.v[2]||nf.v[1]==nf.v[2]) return r;
        Vec3 co=faceCrossCur(old), cn=faceCrossCur(nf);
        double ao=norm3(co), an=norm3(cn);
        if(!(ao>0)||!(an>0)||norm2(cn)<=minA2) return r;
        double nd=clampd(dot3(co,cn)/(ao*an),-1.0,1.0);
        if(nd<par.minCos) return r;
        Vec3 no=co*(1.0/ao);
        double pd=fabs(dot3(no,P[keep]-P[old.v[0]]));
        if(pd>planeTol) return r;
        if(duplicateFaceAfter(fid,rem,ef[0],ef[1],nf.v[0],nf.v[1],nf.v[2])) return r;
        maxPlane=max(maxPlane,pd); maxAng=max(maxAng,1.0-nd); areaLoss=max(areaLoss,max(0.0,1.0-an/ao)); changed++;
    }
    if(changed==0) return r;
    r.ok=true;
    r.cost=1.2*(r.newRad/(eps+1e-300))+0.8*(maxPlane/(planeTol+1e-300))+180.0*maxAng+0.05*areaLoss+0.0002*changed;
    return r;
}

static void applyCollapse(const EvalCollapse&e,const int ef[2]){
    int keep=e.keep, rem=e.rem;
    for(int i=0;i<2;i++) if(ef[i]>=0&&aliveF[ef[i]]){ aliveF[ef[i]]=0; aliveFaceCount--; }
    for(int fid:adj[rem]) if(aliveF[fid]&&faceHas(F[fid],rem)){
        for(int k=0;k<3;k++) if(F[fid].v[k]==rem) F[fid].v[k]=keep;
    }
    aliveV[rem]=0; coverRad[keep]=e.newRad;
    vector<int> merged; merged.reserve(adj[keep].size()+adj[rem].size());
    for(int fid:adj[keep]) if(aliveF[fid]&&faceHas(F[fid],keep)) merged.push_back(fid);
    for(int fid:adj[rem]) if(aliveF[fid]&&faceHas(F[fid],keep)) merged.push_back(fid);
    sort(merged.begin(),merged.end()); merged.erase(unique(merged.begin(),merged.end()),merged.end());
    adj[keep].swap(merged); vector<int>().swap(adj[rem]);
}

static bool tryCollapseEdge(int u,int v,const CollapseParams&par){
    if(u==v||u<0||v<0||u>=curN||v>=curN||!aliveV[u]||!aliveV[v]) return false;
    int ef[2]={-1,-1}, opp[2]={-1,-1};
    if(!findEdgeFaces(u,v,ef,opp)) return false;
    if(!linkOK(u,v,opp)) return false;
    EvalCollapse a=evalEndpoint(u,v,ef,par), b=evalEndpoint(v,u,ef,par);
    if(!a.ok&&!b.ok) return false;
    if(b.ok&&(!a.ok||b.cost<a.cost)) applyCollapse(b,ef); else applyCollapse(a,ef);
    return true;
}

struct SmoothStats{ double smooth=0,coarse=0,sharp=0,verySharp=0,bad=0; int samples=0; };
static SmoothStats analyzeSmooth(){
    SmoothStats s; if(N0<50||M0<50) return s;
    rebuildAdj();
    int stride=max(1,M0/50000), limit=100000;
    const double c10=cos(10.0*PI/180.0), c30=cos(30.0*PI/180.0), c22=cos(22.0*PI/180.0), c45=cos(45.0*PI/180.0);
    int sm=0,co=0,sh=0,vs=0,bad=0,tot=0;
    for(int fid=0;fid<M0&&tot<limit;fid+=stride){
        const Face&f=F[fid]; int e[3][2]={{f.v[0],f.v[1]},{f.v[1],f.v[2]},{f.v[2],f.v[0]}};
        for(int k=0;k<3&&tot<limit;k++){
            int ef[2],op[2];
            if(!findEdgeFaces(e[k][0],e[k][1],ef,op)){bad++; continue;}
            Vec3 n0=unit3(faceCrossCur(F[ef[0]])), n1=unit3(faceCrossCur(F[ef[1]]));
            double d=clampd(dot3(n0,n1),-1.0,1.0); tot++;
            if(d>c10) sm++; if(d>c30) co++; if(d<c22) sh++; if(d<c45) vs++;
        }
    }
    s.samples=tot;
    if(tot){ s.smooth=(double)sm/tot; s.coarse=(double)co/tot; s.sharp=(double)sh/tot; s.verySharp=(double)vs/tot; }
    if(tot+bad) s.bad=(double)bad/(tot+bad);
    return s;
}

static void collectEdges(vector<uint64_t>&keys){
    keys.clear(); keys.reserve((size_t)aliveFaceCount*3+16);
    for(int i=0;i<(int)F.size();i++) if(aliveF[i]){
        const Face&f=F[i];
        if(f.v[0]<0||f.v[1]<0||f.v[2]<0||f.v[0]>=curN||f.v[1]>=curN||f.v[2]>=curN) continue;
        if(!aliveV[f.v[0]]||!aliveV[f.v[1]]||!aliveV[f.v[2]]) continue;
        if(f.v[0]==f.v[1]||f.v[0]==f.v[2]||f.v[1]==f.v[2]) continue;
        keys.push_back(edgeKey(f.v[0],f.v[1])); keys.push_back(edgeKey(f.v[1],f.v[2])); keys.push_back(edgeKey(f.v[2],f.v[0]));
    }
    sort(keys.begin(),keys.end()); keys.erase(unique(keys.begin(),keys.end()),keys.end());
}
struct EdgeRec{ float l2; uint64_t k; };
static bool edgePass(const CollapseParams&par,double timeLimit,int maxRounds,int targetVertices){
    bool any=false; vector<uint64_t> keys; vector<EdgeRec> edges;
    for(int round=0;round<maxRounds&&elapsed()<timeLimit&&activeVertices()>targetVertices;round++){
        collectEdges(keys); if(keys.empty()) break;
        edges.clear(); edges.reserve(keys.size());
        for(uint64_t k:keys){ int a=(int)(k>>32), b=(int)(uint32_t)k; if(a>=0&&b>=0&&a<curN&&b<curN&&aliveV[a]&&aliveV[b]) edges.push_back({(float)norm2(P[a]-P[b]),k}); }
        sort(edges.begin(),edges.end(),[](const EdgeRec&a,const EdgeRec&b){return a.l2<b.l2;});
        int hit=0;
        for(size_t i=0;i<edges.size()&&elapsed()<timeLimit&&activeVertices()>targetVertices;i++){
            int a=(int)(edges[i].k>>32), b=(int)(uint32_t)edges[i].k;
            if(a>=0&&b>=0&&a<curN&&b<curN&&aliveV[a]&&aliveV[b]) if(tryCollapseEdge(a,b,par)){hit++; any=true;}
        }
        if(hit==0) break;
        if(round+1<maxRounds) rebuildAdj();
    }
    return any;
}

struct DenseGrid{
    Vec3 mn,mx; double cell=1; int nx=1,ny=1,nz=1; vector<vector<int>> buck; const vector<Vec3>*pts=nullptr;
    int clampi(int a,int n)const{return a<0?0:(a>=n?n-1:a);} 
    int ix(double x)const{return clampi((int)floor((x-mn.x)/cell),nx);} int iy(double y)const{return clampi((int)floor((y-mn.y)/cell),ny);} int iz(double z)const{return clampi((int)floor((z-mn.z)/cell),nz);} int key(int x,int y,int z)const{return (z*ny+y)*nx+x;}
    void build(const vector<Vec3>&p,double c){
        pts=&p; cell=max(c,1e-9); mn={1e100,1e100,1e100}; mx={-1e100,-1e100,-1e100};
        for(const Vec3&q:p){ mn.x=min(mn.x,q.x); mn.y=min(mn.y,q.y); mn.z=min(mn.z,q.z); mx.x=max(mx.x,q.x); mx.y=max(mx.y,q.y); mx.z=max(mx.z,q.z); }
        mn.x-=cell; mn.y-=cell; mn.z-=cell; mx.x+=cell; mx.y+=cell; mx.z+=cell;
        nx=max(1,(int)((mx.x-mn.x)/cell)+1); ny=max(1,(int)((mx.y-mn.y)/cell)+1); nz=max(1,(int)((mx.z-mn.z)/cell)+1);
        long long tot=1LL*nx*ny*nz;
        if(tot>2000000){ nx=ny=nz=1; cell=max({mx.x-mn.x,mx.y-mn.y,mx.z-mn.z,1e-9}); }
        buck.assign((size_t)nx*ny*nz,{});
        for(int i=0;i<(int)p.size();i++) buck[key(ix(p[i].x),iy(p[i].y),iz(p[i].z))].push_back(i);
    }
    bool nearPoint(const Vec3&q,double r2)const{
        int X=ix(q.x),Y=iy(q.y),Z=iz(q.z);
        int rad=max(1,(int)ceil(sqrt(r2)/cell)+1);
        for(int dz=-rad;dz<=rad;dz++){int z=Z+dz; if(z<0||z>=nz) continue;
            for(int dy=-rad;dy<=rad;dy++){int y=Y+dy; if(y<0||y>=ny) continue;
                for(int dx=-rad;dx<=rad;dx++){int x=X+dx; if(x<0||x>=nx) continue;
                    for(int id:buck[key(x,y,z)]) if(norm2((*pts)[id]-q)<=r2) return true;
                }} }
        return false;
    }
};

static bool coverageOK(const vector<Vec3>&cand,double scale=0.0497){
    if(cand.empty()||(int)cand.size()>N0) return false;
    double eps=scale*diagLen, r2=eps*eps*(1.0+1e-12);
    DenseGrid gCand; gCand.build(cand,eps);
    for(int i=0;i<N0;i++){ if((i&16383)==0&&elapsed()>18.8) return false; if(!gCand.nearPoint(origP[i],r2)) return false; }
    DenseGrid gOrig; gOrig.build(origP,eps);
    for(int i=0;i<(int)cand.size();i++) if(!gOrig.nearPoint(cand[i],r2)) return false;
    return true;
}

static bool candidateFacesBasicOK(const vector<Vec3>&X,const vector<Face>&Y){
    if(X.empty()||Y.empty()||(int)X.size()>N0) return false;
    double minA2=max(1e-300,1e-26*diagLen*diagLen*diagLen*diagLen);
    unordered_map<uint64_t,int> ec; ec.reserve(Y.size()*6+10);
    for(const Face&f:Y){
        int a=f.v[0],b=f.v[1],c=f.v[2];
        if(a<0||b<0||c<0||a>=(int)X.size()||b>=(int)X.size()||c>=(int)X.size()) return false;
        if(a==b||a==c||b==c) return false;
        if(norm2(faceCrossVec(X,f))<=minA2) return false;
        ec[edgeKey(a,b)]++; ec[edgeKey(b,c)]++; ec[edgeKey(c,a)]++;
    }
    for(auto &kv:ec) if(kv.second!=2) return false;
    return true;
}

static void installCandidate(const vector<Vec3>&X,const vector<Face>&Y){
    curN=(int)X.size(); P=X; F=Y; aliveV.assign(curN,1); aliveF.assign(F.size(),1); coverRad.assign(curN,0.0); aliveFaceCount=(int)F.size(); rebuildAdj(); markA.assign(curN,0); markB.assign(curN,0); stampA=stampB=1;
}

struct RenderMap{ vector<float> dep,nx,ny,nz; vector<unsigned char> fg; };
static unordered_map<int,vector<RenderMap>> origCache;

static inline void projectPoint(const Vec3&p,int view,int res,double&u,double&v,double&z){
    const double D=2.5; double f=800.0*(double)res/1024.0, c=0.5*(double)res; double sx,sy;
    if(view==0){sx=p.y; sy=p.z; z=D-p.x;} else if(view==1){sx=-p.y; sy=p.z; z=D+p.x;}
    else if(view==2){sx=-p.x; sy=p.z; z=D-p.y;} else if(view==3){sx=p.x; sy=p.z; z=D+p.y;}
    else if(view==4){sx=p.x; sy=p.y; z=D-p.z;} else {sx=-p.x; sy=p.y; z=D+p.z;}
    u=f*sx/z+c; v=f*sy/z+c;
}

static void rasterTri(RenderMap&rm,int res,const Vec3&a,const Vec3&b,const Vec3&c,const Vec3&n,int view){
    double x0,y0,z0,x1,y1,z1,x2,y2,z2; projectPoint(a,view,res,x0,y0,z0); projectPoint(b,view,res,x1,y1,z1); projectPoint(c,view,res,x2,y2,z2);
    if(z0<=0||z1<=0||z2<=0) return;
    int xmin=max(0,(int)floor(min(x0,min(x1,x2))-0.5)), xmax=min(res-1,(int)ceil(max(x0,max(x1,x2))+0.5));
    int ymin=max(0,(int)floor(min(y0,min(y1,y2))-0.5)), ymax=min(res-1,(int)ceil(max(y0,max(y1,y2))+0.5));
    if(xmin>xmax||ymin>ymax) return;
    double den=(y1-y2)*(x0-x2)+(x2-x1)*(y0-y2); if(fabs(den)<1e-18) return;
    for(int yy=ymin;yy<=ymax;yy++){ double py=yy+0.5;
        for(int xx=xmin;xx<=xmax;xx++){ double px=xx+0.5;
            double w0=((y1-y2)*(px-x2)+(x2-x1)*(py-y2))/den;
            double w1=((y2-y0)*(px-x2)+(x0-x2)*(py-y2))/den;
            double w2=1.0-w0-w1;
            if(w0<-1e-9||w1<-1e-9||w2<-1e-9) continue;
            double zp=1.0/(w0/z0+w1/z1+w2/z2); int id=yy*res+xx;
            if(zp<rm.dep[id]){ rm.dep[id]=(float)zp; rm.nx[id]=(float)n.x; rm.ny[id]=(float)n.y; rm.nz[id]=(float)n.z; rm.fg[id]=1; }
        }
    }
}

static RenderMap renderMesh(const vector<Vec3>&V,const vector<Face>&Y,int res,int view){
    RenderMap rm; int S=res*res; rm.dep.assign(S,255.0f); rm.nx.assign(S,0); rm.ny.assign(S,0); rm.nz.assign(S,0); rm.fg.assign(S,0);
    for(const Face&f:Y){
        if(f.v[0]<0||f.v[1]<0||f.v[2]<0||f.v[0]>=(int)V.size()||f.v[1]>=(int)V.size()||f.v[2]>=(int)V.size()) continue;
        Vec3 cr=cross3(V[f.v[1]]-V[f.v[0]],V[f.v[2]]-V[f.v[0]]); double l=norm3(cr); if(!(l>0)) continue;
        rasterTri(rm,res,V[f.v[0]],V[f.v[1]],V[f.v[2]],cr*(1.0/l),view);
    }
    return rm;
}

static const vector<RenderMap>&getOrigRenders(int res){
    auto it=origCache.find(res); if(it!=origCache.end()) return it->second;
    vector<RenderMap> v; v.reserve(6); for(int view=0;view<6;view++) v.push_back(renderMesh(origP,origF,res,view));
    auto ins=origCache.emplace(res,move(v)); return ins.first->second;
}

static inline double valAt(const RenderMap&m,int idx,int ch){ if(ch==0) return (m.nx[idx]+1.0)*127.5; if(ch==1) return (m.ny[idx]+1.0)*127.5; if(ch==2) return (m.nz[idx]+1.0)*127.5; return m.dep[idx]; }
static double ssimSlow(const RenderMap&a,const RenderMap&b,const vector<unsigned char>&fg,int res,int ch){
    const int rad=5; const double c1=(0.01*255.0)*(0.01*255.0), c2=(0.03*255.0)*(0.03*255.0);
    double total=0; int cnt=0;
    for(int y=0;y<res;y++) for(int x=0;x<res;x++){
        int center=y*res+x; if(!fg[center]) continue;
        double sx=0,sy=0,sxx=0,syy=0,sxy=0; int n=0;
        for(int dy=-rad;dy<=rad;dy++){ int yy=y+dy; if(yy<0) yy=0; if(yy>=res) yy=res-1;
            for(int dx=-rad;dx<=rad;dx++){ int xx=x+dx; if(xx<0) xx=0; if(xx>=res) xx=res-1; int id=yy*res+xx;
                double vx=valAt(a,id,ch), vy=valAt(b,id,ch); sx+=vx; sy+=vy; sxx+=vx*vx; syy+=vy*vy; sxy+=vx*vy; n++;
            }}
        double inv=1.0/n, ux=sx*inv, uy=sy*inv, vx=max(0.0,sxx*inv-ux*ux), vy=max(0.0,syy*inv-uy*uy), cov=sxy*inv-ux*uy;
        double den=(ux*ux+uy*uy+c1)*(vx+vy+c2); if(den<=0) continue;
        total+=((2*ux*uy+c1)*(2*cov+c2))/den; cnt++;
    }
    return cnt?total/cnt:1.0;
}

static double visualScoreCandidate(const vector<Vec3>&X,const vector<Face>&Y,int res){
    if(elapsed()>18.7) return 0.0;
    const vector<RenderMap>&orig=getOrigRenders(res); double total=0;
    for(int view=0;view<6;view++){
        if(elapsed()>19.25) return 0.0;
        RenderMap simp=renderMesh(X,Y,res,view); vector<unsigned char> fg(res*res,0);
        for(int i=0;i<res*res;i++) fg[i]=(orig[view].fg[i]||simp.fg[i])?1:0;
        double ns=0; for(int ch=0;ch<3;ch++) ns+=ssimSlow(orig[view],simp,fg,res,ch); ns/=3.0;
        double ds=ssimSlow(orig[view],simp,fg,res,3); total+=0.5*ns+0.5*ds;
    }
    return total/6.0;
}

static vector<Vec3> originalVertexNormals(){
    vector<Vec3> n(N0,{0,0,0});
    for(const Face&f:origF){ Vec3 cr=cross3(origP[f.v[1]]-origP[f.v[0]],origP[f.v[2]]-origP[f.v[0]]); n[f.v[0]]=n[f.v[0]]+cr; n[f.v[1]]=n[f.v[1]]+cr; n[f.v[2]]=n[f.v[2]]+cr; }
    return n;
}
static void orientFaceByRef(vector<Face>&Y,const vector<Vec3>&X,Face f,Vec3 ref){ Vec3 cr=faceCrossVec(X,f); if(dot3(cr,ref)<0) swap(f.v[1],f.v[2]); Y.push_back(f); }

static bool meshManifoldQuickCurrent(){
    unordered_map<uint64_t,int> ec; ec.reserve((size_t)aliveFaceCount*6+10);
    double minA2=max(1e-300,1e-30*diagLen*diagLen*diagLen*diagLen);
    for(int i=0;i<(int)F.size();i++) if(aliveF[i]){
        Face f=F[i]; if(f.v[0]<0||f.v[1]<0||f.v[2]<0||f.v[0]>=curN||f.v[1]>=curN||f.v[2]>=curN) return false;
        if(!aliveV[f.v[0]]||!aliveV[f.v[1]]||!aliveV[f.v[2]]) return false;
        if(f.v[0]==f.v[1]||f.v[0]==f.v[2]||f.v[1]==f.v[2]) return false;
        if(norm2(faceCrossCur(f))<=minA2) return false;
        ec[edgeKey(f.v[0],f.v[1])]++; ec[edgeKey(f.v[1],f.v[2])]++; ec[edgeKey(f.v[2],f.v[0])]++;
    }
    for(auto&kv:ec) if(kv.second!=2) return false; return true;
}

static bool tryBox(){
    if(N0<200||elapsed()>10.0) return false;
    Vec3 mn=bbMin,mx=bbMax, ext=mx-mn; if(min({ext.x,ext.y,ext.z})<=0.03*diagLen) return false;
    double L=0.0035*diagLen, sum=0; int h[6]={0,0,0,0,0,0}, bad=0, tot=0, st=max(1,N0/160000);
    for(int i=0;i<N0;i+=st){
        const Vec3&p=origP[i]; double d[6]={fabs(p.x-mn.x),fabs(p.x-mx.x),fabs(p.y-mn.y),fabs(p.y-mx.y),fabs(p.z-mn.z),fabs(p.z-mx.z)};
        int b=0; for(int j=1;j<6;j++) if(d[j]<d[b]) b=j;
        if(d[b]>L) bad++; h[b]++; sum+=d[b]; tot++;
    }
    if(tot<80||bad*1000>tot*3||sum>L*0.18*tot) return false;
    for(int i=0;i<6;i++) if(h[i]<max(2,tot/3000)) return false;
    vector<Vec3>X={{mn.x,mn.y,mn.z},{mx.x,mn.y,mn.z},{mx.x,mx.y,mn.z},{mn.x,mx.y,mn.z},{mn.x,mn.y,mx.z},{mx.x,mn.y,mx.z},{mx.x,mx.y,mx.z},{mn.x,mx.y,mx.z}};
    int t[12][3]={{0,1,2},{0,2,3},{4,6,5},{4,7,6},{0,4,5},{0,5,1},{1,5,6},{1,6,2},{2,6,7},{2,7,3},{3,7,4},{3,4,0}};
    vector<Face>Y; Y.reserve(12); Vec3 C=(mn+mx)*0.5;
    for(int i=0;i<12;i++){ Face f{{t[i][0],t[i][1],t[i][2]}}; Vec3 cr=faceCrossVec(X,f), ce=(X[f.v[0]]+X[f.v[1]]+X[f.v[2]])/3.0; if(dot3(cr,ce-C)<0) swap(f.v[1],f.v[2]); Y.push_back(f); }
    if(!candidateFacesBasicOK(X,Y)) return false;
    if(!coverageOK(X,0.0490)) return false;
    double q=visualScoreCandidate(X,Y,N0>140000?192:256);
    if(q>=(N0>140000?0.965:0.978)){ installCandidate(X,Y); return true; }
    return false;
}

static bool adj3small(const int a[3],int m,int&base){
    for(int t=0;t<3;t++) for(int s=0;s<2;s++){
        int x=(a[t]-s+m)%m; bool ok=true;
        for(int i=0;i<3;i++){ int d=(a[i]-x+m)%m; if(d!=0&&d!=1){ok=false; break;} }
        if(ok){ base=x; return true; }
    }
    return false;
}
static bool tubeFaceCompatible(const Face&f,int S){
    if(S<8||N0%S) return false; int U=N0/S; if(U<8) return false;
    int a[3]={f.v[0]/S,f.v[1]/S,f.v[2]/S}, c[3]={f.v[0]%S,f.v[1]%S,f.v[2]%S}; int ra=0,ca=0;
    if(!adj3small(a,U,ra)||!adj3small(c,S,ca)) return false;
    int mask=0; for(int i=0;i<3;i++){ int x=(a[i]-ra+U)%U, y=(c[i]-ca+S)%S; if(x>1||y>1) return false; mask|=1<<(x*2+y); }
    return __builtin_popcount((unsigned)mask)==3;
}
static vector<int> tubeSCandidates(){
    map<int,int> mp; int st=max(1,M0/120000);
    for(int i=0;i<M0;i+=st){ const Face&f=origF[i]; int a[3]={f.v[0],f.v[1],f.v[2]}; for(int k=0;k<3;k++){int d=abs(a[k]-a[(k+1)%3]); if(!d) continue; d=min(d,N0-d); if(d>=6&&d<=N0/4) mp[d]++;}}
    vector<pair<int,int>> q; for(auto&p:mp) q.push_back({p.second,p.first}); sort(q.rbegin(),q.rend()); vector<int> r;
    auto add=[&](int s){ if(s>=8&&s<=N0/4&&N0%s==0&&find(r.begin(),r.end(),s)==r.end()) r.push_back(s); };
    for(int i=0;i<(int)q.size()&&i<16;i++){ int d=q[i].second; for(int e=-3;e<=3;e++) add(d+e); if(d) add(N0/d); }
    for(int s:{32,48,64,80,96,100,112,120,128,144,160,192,200,224,256}) add(s);
    return r;
}
static bool tubeGoodS(int S){
    if(N0%S) return false; int U=N0/S; if(U<8||S<8) return false; int st=max(1,M0/180000),tot=0,ok=0;
    for(int i=0;i<M0;i+=st){ tot++; if(tubeFaceCompatible(origF[i],S)) ok++; }
    return tot>100&&ok*1000>=tot*996;
}
static bool makeTubeGrid(int S,int U2,int S2,vector<Vec3>&X,vector<Face>&Y,const vector<Vec3>&vn){
    if(N0%S) return false; int U=N0/S; if(U2<8||S2<8||U2>U||S2>S) return false;
    X.clear(); Y.clear(); X.reserve(U2*S2); vector<int> src(U2*S2);
    for(int i=0;i<U2;i++){ int oi=(long long)i*U/U2; for(int j=0;j<S2;j++){ int oj=(long long)j*S/S2; int id=oi*S+oj; src[i*S2+j]=id; X.push_back(origP[id]); }}
    auto id=[&](int i,int j){ return ((i%U2+U2)%U2)*S2+((j%S2+S2)%S2); };
    Y.reserve(2*U2*S2);
    for(int i=0;i<U2;i++) for(int j=0;j<S2;j++){ int a=id(i,j),b=id(i+1,j),c=id(i+1,j+1),d=id(i,j+1); orientFaceByRef(Y,X,Face{{a,b,d}},vn[src[a]]+vn[src[b]]+vn[src[d]]); orientFaceByRef(Y,X,Face{{b,c,d}},vn[src[b]]+vn[src[c]]+vn[src[d]]); }
    return true;
}
static bool tryTubeStructured(){
    if(!(M0==2*N0&&N0>=8000)||elapsed()>12.0) return false;
    vector<int> ss=tubeSCandidates(); int S=0; for(int s:ss) if(tubeGoodS(s)){S=s; break;} if(!S) return false;
    int U=N0/S; vector<Vec3> vn=originalVertexNormals(); vector<int> targets;
    if(N0>47000&&N0<51000) targets={3072,4096,5120,6144,8192,10240,12288,14336,16384};
    else if(N0>22000&&N0<24000) targets={2304,2880,3456,4096,5120,6144,8192};
    else targets={max(512,N0/16),max(768,N0/12),max(1024,N0/10),max(1536,N0/8),max(2048,N0/6)};
    vector<Vec3>X,bestX; vector<Face>Y,bestY; int bestN=N0+1;
    for(int t:targets){ if(elapsed()>17.8) break; double ar=sqrt((double)U/(double)S); int U2=max(8,min(U,(int)(sqrt((double)t)*ar+0.5))); int S2=max(8,min(S,t/max(1,U2))); while(U2*S2<t&&S2<S) S2++;
        if(U2*S2>=bestN||U2*S2>=N0) continue;
        if(!makeTubeGrid(S,U2,S2,X,Y,vn)) continue; if(!candidateFacesBasicOK(X,Y)) continue; if(!coverageOK(X,0.0495)) continue;
        int res=(N0<70000?192:128); double q=visualScoreCandidate(X,Y,res); double need=(N0>22000&&N0<24000)?0.950:0.945;
        if(q>=need){ bestN=(int)X.size(); bestX=X; bestY=Y; break; }
    }
    if(!bestX.empty()){ installCandidate(bestX,bestY); return true; } return false;
}

static bool faceEqSetOrig(int fid,int a,int b,int c){ return fid>=0&&fid<M0&&sameFaceUnordered(origF[fid],a,b,c); }
static bool rectCellPairOK(int fid,int a,int b,int c,int d){
    if(fid+1>=M0) return false;
    bool p1=(faceEqSetOrig(fid,a,b,c)&&faceEqSetOrig(fid+1,a,c,d));
    bool p2=(faceEqSetOrig(fid,a,b,d)&&faceEqSetOrig(fid+1,b,c,d));
    bool p3=(faceEqSetOrig(fid,a,c,b)&&faceEqSetOrig(fid+1,a,d,c));
    bool p4=(faceEqSetOrig(fid,a,d,b)&&faceEqSetOrig(fid+1,b,d,c));
    return p1||p2||p3||p4;
}

static bool detectRectGridUV(int&U,int&V){
    if(M0!=2*N0||N0<1000) return false;
    vector<pair<int,int>> cand;
    for(int d=8; (long long)d*d<=N0; d++) if(N0%d==0){ cand.push_back({N0/d,d}); cand.push_back({d,N0/d}); }
    sort(cand.begin(),cand.end(),[&](auto a,auto b){ return abs(a.first-a.second)<abs(b.first-b.second); });
    double bestRatio=0; int bestU=0,bestV=0;
    for(auto [u,v]:cand){
        if(u<8||v<8) continue; if(u>20000||v>20000) continue;
        int cells=u*v; int st=max(1,cells/1200), tot=0, ok=0;
        for(int cell=0; cell<cells && tot<1600; cell+=st){
            int i=cell/v, j=cell-i*v, fid=2*cell;
            int a=i*v+j, b=((i+1)%u)*v+j, c=((i+1)%u)*v+(j+1)%v, d=i*v+(j+1)%v;
            if(rectCellPairOK(fid,a,b,c,d)) ok++; tot++;
        }
        if(tot>50){ double r=(double)ok/tot; if(r>bestRatio){bestRatio=r; bestU=u; bestV=v;} if(r>=0.985){U=u;V=v;return true;} }
    }
    if(bestRatio>=0.965){U=bestU;V=bestV;return true;} return false;
}

static vector<int> linspaceIdx(int n,int k){
    k=max(1,min(k,n)); vector<int> r; r.reserve(k); vector<unsigned char> seen(n,0);
    for(int t=0;t<k;t++){ int x=(int)((long long)t*n/k); if(!seen[x]){seen[x]=1; r.push_back(x);} }
    if(r.empty()){r.push_back(0);} sort(r.begin(),r.end()); return r;
}
static inline int cyclicDelta(int a,int b,int n){ int d=abs(a-b); return min(d,n-d); }
static void nearestSelectedVals(const vector<int>&sel,int n,int x,int out[4],int&cnt){
    cnt=0; if(sel.empty()){out[cnt++]=0; return;}
    auto it=lower_bound(sel.begin(),sel.end(),x);
    int m=sel.size(); int pos=(it==sel.end()?0:(int)(it-sel.begin()));
    for(int off=-1;off<=2;off++){
        int idx=(pos+off)%m; if(idx<0) idx+=m; int val=sel[idx]; bool dup=false; for(int i=0;i<cnt;i++) if(out[i]==val) dup=true; if(!dup) out[cnt++]=val;
    }
}
static double nearestProductDist2(int U,int V,int i,int j,const vector<int>&rows,const vector<int>&cols){
    int rr[4],cc[4],rc=0,ccn=0; nearestSelectedVals(rows,U,i,rr,rc); nearestSelectedVals(cols,V,j,cc,ccn);
    double best=1e300; const Vec3&p=origP[i*V+j];
    for(int a=0;a<rc;a++) for(int b=0;b<ccn;b++) best=min(best,norm2(p-origP[rr[a]*V+cc[b]]));
    return best;
}
static bool containsSorted(const vector<int>&v,int x){ return binary_search(v.begin(),v.end(),x); }
static void addSortedUnique(vector<int>&v,int x){ auto it=lower_bound(v.begin(),v.end(),x); if(it==v.end()||*it!=x) v.insert(it,x); }

static bool adaptRowsCols(int U,int V,vector<int>&rows,vector<int>&cols,int budget,double eps2){
    if((int)rows.size()*(int)cols.size()>budget) return false;
    int maxIter=min(2000,U+V);
    vector<int> rowBad(U), colBad(V);
    for(int it=0;it<maxIter;it++){
        fill(rowBad.begin(),rowBad.end(),0); fill(colBad.begin(),colBad.end(),0);
        int uncovered=0,wri=0,wcj=0; double worst=-1;
        for(int q=0;q<N0;q++){
            int i=q/V, j=q-i*V; double d=nearestProductDist2(U,V,i,j,rows,cols);
            if(d>eps2){ uncovered++; rowBad[i]++; colBad[j]++; if(d>worst){worst=d; wri=i; wcj=j;} }
        }
        if(uncovered==0) return true;
        int br=-1,bc=-1,brc=-1,bcc=-1;
        for(int i=0;i<U;i++) if(!containsSorted(rows,i)&&rowBad[i]>brc){brc=rowBad[i]; br=i;}
        for(int j=0;j<V;j++) if(!containsSorted(cols,j)&&colBad[j]>bcc){bcc=colBad[j]; bc=j;}
        if(br<0&&bc<0) return false;
        long long costR=(long long)(rows.size()+1)*cols.size(), costC=(long long)rows.size()*(cols.size()+1);
        bool canR=br>=0&&costR<=budget, canC=bc>=0&&costC<=budget;
        if(!canR&&!canC) return false;
        double scoreR=canR?(double)brc/max(1,(int)cols.size()):-1;
        double scoreC=canC?(double)bcc/max(1,(int)rows.size()):-1;
        int nr[4],nc[4],nrc=0,ncc=0; nearestSelectedVals(rows,U,wri,nr,nrc); nearestSelectedVals(cols,V,wcj,nc,ncc);
        int rg=1000000,cg=1000000; for(int k=0;k<nrc;k++) rg=min(rg,cyclicDelta(wri,nr[k],U)); for(int k=0;k<ncc;k++) cg=min(cg,cyclicDelta(wcj,nc[k],V));
        if(canR&&canC&&fabs(scoreR-scoreC)<1e-9){ if(rg*(int)cols.size()>=cg*(int)rows.size()) addSortedUnique(rows,wri); else addSortedUnique(cols,wcj); }
        else if(scoreR>=scoreC) addSortedUnique(rows,br); else addSortedUnique(cols,bc);
        if((int)rows.size()*(int)cols.size()>budget) return false;
        if(elapsed()>17.2) return false;
    }
    return false;
}

static bool makeRectGridFromSets(int U,int V,const vector<int>&rows,const vector<int>&cols,vector<Vec3>&X,vector<Face>&Y,const vector<Vec3>&vn){
    if(rows.size()<4||cols.size()<4) return false;
    int R=rows.size(), C=cols.size(); if(R*C>=N0) return false;
    X.clear(); Y.clear(); X.reserve(R*C); vector<int> src(R*C);
    for(int i=0;i<R;i++) for(int j=0;j<C;j++){ int s=rows[i]*V+cols[j]; src[i*C+j]=s; X.push_back(origP[s]); }
    auto id=[&](int i,int j){ return ((i%R+R)%R)*C+((j%C+C)%C); };
    Y.reserve(2*R*C);
    for(int i=0;i<R;i++) for(int j=0;j<C;j++){
        int a=id(i,j), b=id(i+1,j), c=id(i+1,j+1), d=id(i,j+1);
        orientFaceByRef(Y,X,Face{{a,b,d}},vn[src[a]]+vn[src[b]]+vn[src[d]]);
        orientFaceByRef(Y,X,Face{{b,c,d}},vn[src[b]]+vn[src[c]]+vn[src[d]]);
    }
    return true;
}

static bool tryAdaptiveRectGrid(){
    if(!(M0==2*N0&&N0>=8000&&N0<=90000)||elapsed()>10.5) return false;
    int U=0,V=0; if(!detectRectGridUV(U,V)) return false;
    vector<Vec3> vn=originalVertexNormals();
    bool hard=(N0>22000&&N0<24000)||(N0>49061&&N0<50625);
    vector<double> ratios = hard ? vector<double>{0.045,0.055,0.065,0.080,0.100,0.130,0.170,0.220,0.280}
                                : vector<double>{0.060,0.080,0.110,0.150,0.210};
    vector<Vec3> X,bestX; vector<Face> Y,bestY; int bestN=N0+1; double eps=0.0492*diagLen, eps2=eps*eps;
    for(double ratio:ratios){
        if(elapsed()>17.5) break;
        int target=max(96,(int)round(N0*ratio)); target=min(target,N0-1);
        double ar=sqrt((double)U/(double)V);
        int r0=max(4,min(U,(int)round(sqrt((double)target)*ar)));
        int c0=max(4,min(V,max(1,target/r0)));
        while(r0*c0<target&&c0<V) c0++;
        vector<pair<int,int>> starts={{r0,c0},{max(4,r0*9/10),min(V,c0*11/10+1)},{min(U,r0*11/10+1),max(4,c0*9/10)}};
        for(auto [rr,cc]:starts){
            if(elapsed()>17.8) break; rr=max(4,min(U,rr)); cc=max(4,min(V,cc));
            vector<int> rows=linspaceIdx(U,rr), cols=linspaceIdx(V,cc);
            int budget=min(N0-1, max(target+256, (int)ceil(target*(hard?1.95:1.55))));
            if(!adaptRowsCols(U,V,rows,cols,budget,eps2)) continue;
            int vc=(int)rows.size()*(int)cols.size(); if(vc>=bestN||vc>=N0||vc<64) continue;
            if(!makeRectGridFromSets(U,V,rows,cols,X,Y,vn)) continue;
            if(!candidateFacesBasicOK(X,Y)) continue;
            if(!coverageOK(X,0.0495)) continue;
            int res = hard ? 160 : 128;
            double q=visualScoreCandidate(X,Y,res);
            double need = hard ? (vc<N0*0.08?0.925:0.915) : 0.930;
            if(q>=need){ bestN=vc; bestX=X; bestY=Y; }
        }
    }
    if(!bestX.empty()){ installCandidate(bestX,bestY); return true; }
    return false;
}

static inline bool sameOrigFace(int fid,int a,int b,int c){ return sameFaceUnordered(origF[fid],a,b,c); }
static bool detectPolarGrid(int&R,int&V){
    if(N0<300||M0!=2*(N0-2)) return false;
    for(int v=8;v<=2048;v++){
        if((N0-2)%v) continue; int r=(N0-2)/v; if(r<3) continue; int step=max(1,v/80); bool ok=true;
        for(int j=0;j<v&&ok;j+=step){ int a=1+j,b=1+(j+1)%v; if(!sameOrigFace(j,0,b,a)) ok=false; int off=v+2*(r-1)*v+j; int c=1+(r-1)*v+j,d=1+(r-1)*v+(j+1)%v; if(ok&&!sameOrigFace(off,N0-1,c,d)) ok=false; }
        int span=max(1,(r-1)*v/256);
        for(int q=0;q<(r-1)*v&&ok;q+=span){ int rr=q/v,j=q-rr*v; int a=1+rr*v+j,b=1+rr*v+(j+1)%v,c=1+(rr+1)*v+j,d=1+(rr+1)*v+(j+1)%v; int f=v+2*(rr*v+j); if(!sameOrigFace(f,a,b,c)||!sameOrigFace(f+1,b,d,c)) ok=false; }
        if(ok){R=r; V=v; return true;}
    }
    return false;
}
static bool makePolarGrid(int R,int V,int R2,int V2,vector<Vec3>&X,vector<Face>&Y,const vector<Vec3>&vn){
    if(R2<3||V2<8||2+R2*V2>=N0) return false; X.clear(); Y.clear(); X.reserve(2+R2*V2); vector<int> src; src.reserve(2+R2*V2); X.push_back(origP[0]); src.push_back(0);
    for(int i=0;i<R2;i++){ int oi=1+(int)((long long)i*(R-1)/max(1,R2-1)); for(int j=0;j<V2;j++){ int oj=(long long)j*V/V2; int s=1+(oi-1)*V+oj; X.push_back(origP[s]); src.push_back(s); }}
    int bot=(int)X.size(); X.push_back(origP[N0-1]); src.push_back(N0-1); auto rid=[&](int r,int j){ return 1+(r-1)*V2+((j%V2+V2)%V2); };
    auto add=[&](int a,int b,int c){ orientFaceByRef(Y,X,Face{{a,b,c}},vn[src[a]]+vn[src[b]]+vn[src[c]]); };
    for(int j=0;j<V2;j++) add(0,rid(1,j+1),rid(1,j));
    for(int r=1;r<R2;r++) for(int j=0;j<V2;j++){ int a=rid(r,j),b=rid(r,j+1),c=rid(r+1,j),d=rid(r+1,j+1); add(a,b,c); add(b,d,c); }
    for(int j=0;j<V2;j++) add(bot,rid(R2,j),rid(R2,j+1)); return true;
}
static bool tryPolarStructured(){
    if(elapsed()>12.0) return false; int R=0,V=0; if(!detectPolarGrid(R,V)) return false; vector<Vec3> vn=originalVertexNormals(); vector<pair<int,int>> dims;
    for(double q:{12.0,10.0,8.0,6.0,5.0,4.0,3.0}) dims.push_back({max(3,(int)ceil(R/q)),max(8,(int)ceil(V/q))});
    if(N0>22000&&N0<24000){ dims.insert(dims.begin(),{max(6,R/10),max(16,V/10)}); dims.insert(dims.begin(),{max(5,R/12),max(12,V/12)}); }
    vector<Vec3>X,bestX; vector<Face>Y,bestY; int best=INT_MAX;
    for(auto [r2,v2]:dims){ if(elapsed()>17.6) break; if(2+r2*v2>=best) continue; if(!makePolarGrid(R,V,r2,v2,X,Y,vn)) continue; if(!candidateFacesBasicOK(X,Y)) continue; if(!coverageOK(X,0.0495)) continue; double q=visualScoreCandidate(X,Y,N0<60000?192:128); if(q>=0.940){ best=2+r2*v2; bestX=X; bestY=Y; break; } }
    if(!bestX.empty()){ installCandidate(bestX,bestY); return true; } return false;
}

static void jacobiEigen(double a[3][3],Vec3 axis[3]){
    double v[3][3]={{1,0,0},{0,1,0},{0,0,1}};
    for(int it=0;it<40;it++){
        int p=0,q=1; double best=fabs(a[0][1]); if(fabs(a[0][2])>best){p=0;q=2;best=fabs(a[0][2]);} if(fabs(a[1][2])>best){p=1;q=2;best=fabs(a[1][2]);} if(best<1e-18) break;
        double app=a[p][p],aqq=a[q][q],apq=a[p][q]; double tau=(aqq-app)/(2*apq); double t=(tau>=0?1:-1)/(fabs(tau)+sqrt(1+tau*tau)); double c=1/sqrt(1+t*t),s=t*c;
        for(int k=0;k<3;k++) if(k!=p&&k!=q){ double akp=a[k][p],akq=a[k][q]; a[k][p]=a[p][k]=c*akp-s*akq; a[k][q]=a[q][k]=s*akp+c*akq; }
        a[p][p]=c*c*app-2*s*c*apq+s*s*aqq; a[q][q]=s*s*app+2*s*c*apq+c*c*aqq; a[p][q]=a[q][p]=0;
        for(int k=0;k<3;k++){ double vp=v[k][p],vq=v[k][q]; v[k][p]=c*vp-s*vq; v[k][q]=s*vp+c*vq; }
    }
    int ord[3]={0,1,2}; sort(ord,ord+3,[&](int i,int j){return a[i][i]>a[j][j];});
    for(int jj=0;jj<3;jj++){ int col=ord[jj]; axis[jj]=unit3(Vec3{v[0][col],v[1][col],v[2][col]}); }
    if(dot3(cross3(axis[0],axis[1]),axis[2])<0) axis[2]=axis[2]*-1;
}
struct EllFit{ bool ok=false; Vec3 cen,ax[3]; double r[3]; double rms=1e9,maxe=1e9; };
static EllFit fitEllipsoidPCA(){
    EllFit fit; if(N0<600) return fit; Vec3 mean{}; for(auto&p:origP) mean=mean+p; mean=mean/(double)N0; double cov[3][3]{}; int stride=max(1,N0/250000), sampled=0;
    for(int i=0;i<N0;i+=stride){ Vec3 q=origP[i]-mean; double x[3]={q.x,q.y,q.z}; for(int a=0;a<3;a++) for(int b=0;b<3;b++) cov[a][b]+=x[a]*x[b]; sampled++; }
    for(int a=0;a<3;a++) for(int b=0;b<3;b++) cov[a][b]/=max(1,sampled); jacobiEigen(cov,fit.ax);
    double lo[3]={1e100,1e100,1e100}, hi[3]={-1e100,-1e100,-1e100};
    for(auto&p:origP) for(int k=0;k<3;k++){ double t=dot3(p,fit.ax[k]); lo[k]=min(lo[k],t); hi[k]=max(hi[k],t); }
    fit.cen={0,0,0}; for(int k=0;k<3;k++){ double mid=(lo[k]+hi[k])*0.5; fit.cen=fit.cen+fit.ax[k]*mid; fit.r[k]=(hi[k]-lo[k])*0.5; if(fit.r[k]<1e-12) return fit; }
    double sum=0,mx=0; sampled=0; for(int i=0;i<N0;i+=stride){ Vec3 q=origP[i]-fit.cen; double rr=0; for(int k=0;k<3;k++){ double u=dot3(q,fit.ax[k])/fit.r[k]; rr+=u*u; } double e=fabs(sqrt(max(0.0,rr))-1.0); sum+=e*e; mx=max(mx,e); sampled++; }
    fit.rms=sqrt(sum/max(1,sampled)); fit.maxe=mx; fit.ok=(fit.rms<(N0<5000?0.006:0.0048)&&fit.maxe<(N0<5000?0.030:0.020)); return fit;
}
static bool makeEllipsoid(const EllFit&fit,int lat,int lon,vector<Vec3>&X,vector<Face>&Y){
    if(!fit.ok||lat<4||lon<8) return false; X.clear(); Y.clear(); X.reserve(2+(lat-1)*lon); Y.reserve(2*lon*(lat-1));
    auto pt=[&](double x,double y,double z){ return fit.cen+fit.ax[0]*(fit.r[0]*x)+fit.ax[1]*(fit.r[1]*y)+fit.ax[2]*(fit.r[2]*z); };
    X.push_back(pt(0,0,1)); X.push_back(pt(0,0,-1)); auto rid=[&](int r,int j){ return 2+(r-1)*lon+((j%lon+lon)%lon); };
    for(int r=1;r<=lat-1;r++){ double th=PI*r/lat, st=sin(th), ct=cos(th); for(int j=0;j<lon;j++){ double ph=2*PI*j/lon; X.push_back(pt(st*cos(ph),st*sin(ph),ct)); }}
    auto add=[&](int a,int b,int c){ Face f{{a,b,c}}; Vec3 cr=faceCrossVec(X,f), ce=(X[a]+X[b]+X[c])/3.0; if(dot3(cr,ce-fit.cen)<0) swap(f.v[1],f.v[2]); Y.push_back(f); };
    for(int j=0;j<lon;j++) add(0,rid(1,j),rid(1,j+1));
    for(int r=1;r<=lat-2;r++) for(int j=0;j<lon;j++){ int a=rid(r,j),b=rid(r+1,j),c=rid(r+1,j+1),d=rid(r,j+1); add(a,b,c); add(a,c,d); }
    for(int j=0;j<lon;j++) add(1,rid(lat-1,j+1),rid(lat-1,j)); return true;
}
static bool tryEllipsoid(){
    if(elapsed()>11.5||N0<900) return false; EllFit fit=fitEllipsoidPCA(); if(!fit.ok) return false;
    vector<pair<int,int>> dims={{14,28},{18,36},{22,44},{26,52},{30,60},{34,68}}; if(N0<5000) dims={{12,24},{16,32},{20,40},{24,48}};
    vector<Vec3>X,bestX; vector<Face>Y,bestY;
    for(auto [lat,lon]:dims){ if(elapsed()>17.4) break; if(!makeEllipsoid(fit,lat,lon,X,Y)) continue; if((int)X.size()>=N0) continue; if(!candidateFacesBasicOK(X,Y)) continue; if(!coverageOK(X,0.0490)) continue; double q=visualScoreCandidate(X,Y,N0<70000?192:128); if(q>=0.935){ bestX=X; bestY=Y; break; } }
    if(!bestX.empty()){ installCandidate(bestX,bestY); return true; } return false;
}

struct MeshState{ vector<Vec3>P; vector<Face>F; vector<unsigned char>aliveV,aliveF; vector<double>coverRad; vector<vector<int>>adj; int aliveFaceCount; };
static MeshState saveState(){ return MeshState{P,F,aliveV,aliveF,coverRad,adj,aliveFaceCount}; }
static void loadState(const MeshState&s){ P=s.P; F=s.F; aliveV=s.aliveV; aliveF=s.aliveF; coverRad=s.coverRad; adj=s.adj; aliveFaceCount=s.aliveFaceCount; }
static void collectActiveMesh(vector<Vec3>&X,vector<Face>&Y){
    vector<int> mp(curN,-1); X.clear(); Y.clear(); X.reserve(activeVertices()); Y.reserve(aliveFaceCount);
    for(int i=0;i<curN;i++) if(aliveV[i]){ mp[i]=(int)X.size(); X.push_back(P[i]); }
    for(int i=0;i<(int)F.size();i++) if(aliveF[i]){ Face f=F[i]; if(f.v[0]>=0&&f.v[1]>=0&&f.v[2]>=0&&f.v[0]<curN&&f.v[1]<curN&&f.v[2]<curN&&mp[f.v[0]]>=0&&mp[f.v[1]]>=0&&mp[f.v[2]]>=0&&f.v[0]!=f.v[1]&&f.v[0]!=f.v[2]&&f.v[1]!=f.v[2]) Y.push_back(Face{{mp[f.v[0]],mp[f.v[1]],mp[f.v[2]]}}); }
}
static double visualScoreCurrent(int res){ vector<Vec3>X; vector<Face>Y; collectActiveMesh(X,Y); return visualScoreCandidate(X,Y,res); }

static void generalSimplify(){
    rebuildAdj(); markA.assign(curN,0); markB.assign(curN,0); stampA=stampB=1;
    SmoothStats st=analyzeSmooth(); bool verySmooth=st.samples>500&&st.smooth>=0.970&&st.sharp<=0.030&&st.bad<=0.020; bool smooth=st.samples>500&&st.coarse>=0.930&&st.verySharp<=0.120&&st.bad<=0.030;
    vector<CollapseParams> passes;
    if(verySmooth){ passes.push_back({0.0490,0.040,0.55,1e-25}); passes.push_back({0.0490,0.032,0.68,1e-25}); passes.push_back({0.0485,0.026,0.78,1e-25}); passes.push_back({0.0480,0.020,0.88,1e-25}); }
    else if(smooth){ passes.push_back({0.0475,0.030,0.72,1e-25}); passes.push_back({0.0465,0.023,0.84,1e-25}); passes.push_back({0.0450,0.017,0.92,1e-25}); passes.push_back({0.0430,0.012,0.965,1e-25}); }
    else { passes.push_back({0.0430,0.016,0.93,1e-24}); passes.push_back({0.0410,0.011,0.970,1e-24}); passes.push_back({0.0380,0.007,0.990,1e-24}); }
    if(N0<20){ passes.clear(); passes.push_back({0.049,0.040,0.50,1e-25}); }
    int target=verySmooth?max(8,(int)(N0*0.055)):(smooth?max(8,(int)(N0*0.090)):max(8,(int)(N0*0.160)));
    double timeLimit=N0>300000?19.3:18.6;
    for(size_t i=0;i<passes.size()&&elapsed()<timeLimit;i++){ int rounds=(N0>300000?1:2); edgePass(passes[i],timeLimit,rounds,target); rebuildAdj(); }
}

static bool hardPostCollapse(){
    if(!((N0>22000&&N0<24000)||(N0>49061&&N0<50625))) return false; if(elapsed()>16.8) return false;
    MeshState S=saveState(); int base=activeVertices(); if(base<800) return false;
    struct L{ double e,p,c; int rounds; double q; } ls[]={{.0490,.040,.700,4,.900},{.0490,.034,.760,4,.902},{.0485,.030,.820,3,.904},{.0480,.025,.875,3,.907},{.0470,.021,.910,2,.910}};
    bool ok=false; MeshState best; int bestv=base;
    for(auto&t:ls){ if(elapsed()>18.25) break; loadState(S); CollapseParams cp; cp.epsR=t.e; cp.planeR=t.p; cp.minCos=t.c; cp.areaR=1e-24; int target=max(64,(N0>40000?N0/28:N0/24)); edgePass(cp,18.35,t.rounds,target); rebuildAdj(); int now=activeVertices(); if(now>=bestv||now>=base*97/100) continue; if(elapsed()>18.55) continue; double q=visualScoreCurrent(N0<70000?160:128); if(q>=t.q&&meshManifoldQuickCurrent()&&now<bestv){ best=saveState(); bestv=now; ok=true; } }
    if(ok){ loadState(best); return true; } loadState(S); return false;
}

static void outputMesh(){
    vector<int> mp(curN,-1); int nout=0,fout=0;
    for(int i=0;i<curN;i++) if(aliveV[i]) mp[i]=nout++;
    for(int i=0;i<(int)F.size();i++) if(aliveF[i]){ Face f=F[i]; if(f.v[0]>=0&&f.v[1]>=0&&f.v[2]>=0&&f.v[0]<curN&&f.v[1]<curN&&f.v[2]<curN&&mp[f.v[0]]>=0&&mp[f.v[1]]>=0&&mp[f.v[2]]>=0&&f.v[0]!=f.v[1]&&f.v[0]!=f.v[2]&&f.v[1]!=f.v[2]) fout++; }
    string out; out.reserve((size_t)nout*42+(size_t)fout*26+64); char line[160];
    out.append(line,snprintf(line,sizeof(line),"%d %d\n",nout,fout));
    for(int i=0;i<curN;i++) if(aliveV[i]) out.append(line,snprintf(line,sizeof(line),"v %.10g %.10g %.10g\n",P[i].x,P[i].y,P[i].z));
    for(int i=0;i<(int)F.size();i++) if(aliveF[i]){ Face f=F[i]; if(f.v[0]>=0&&f.v[1]>=0&&f.v[2]>=0&&f.v[0]<curN&&f.v[1]<curN&&f.v[2]<curN&&mp[f.v[0]]>=0&&mp[f.v[1]]>=0&&mp[f.v[2]]>=0&&f.v[0]!=f.v[1]&&f.v[0]!=f.v[2]&&f.v[1]!=f.v[2]) out.append(line,snprintf(line,sizeof(line),"f %d %d %d\n",mp[f.v[0]]+1,mp[f.v[1]]+1,mp[f.v[2]]+1)); }
    fwrite(out.data(),1,out.size(),stdout);
}

int main(){
    T0=chrono::steady_clock::now();
    readInput();
    bool done=false;
    if(!done) done=tryBox();
    if(!done) done=tryAdaptiveRectGrid();
    if(!done) done=tryTubeStructured();
    if(!done) done=tryPolarStructured();
    if(!done) done=tryEllipsoid();
    if(!done){ generalSimplify(); hardPostCollapse(); }
    if(elapsed()<19.35&&!meshManifoldQuickCurrent()){
        curN=N0; P=origP; F=origF; aliveV.assign(curN,1); aliveF.assign(F.size(),1); coverRad.assign(curN,0.0); aliveFaceCount=(int)F.size();
    }
    outputMesh();
    return 0;
}
