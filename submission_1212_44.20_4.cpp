#include <bits/stdc++.h>
using namespace std;

struct Vec3 {
    double x=0, y=0, z=0;
};
static inline Vec3 operator+(const Vec3&a,const Vec3&b){return {a.x+b.x,a.y+b.y,a.z+b.z};}
static inline Vec3 operator-(const Vec3&a,const Vec3&b){return {a.x-b.x,a.y-b.y,a.z-b.z};}
static inline Vec3 operator*(const Vec3&a,double s){return {a.x*s,a.y*s,a.z*s};}
static inline Vec3 operator/(const Vec3&a,double s){return {a.x/s,a.y/s,a.z/s};}
static inline double dot3(const Vec3&a,const Vec3&b){return a.x*b.x+a.y*b.y+a.z*b.z;}
static inline Vec3 cross3(const Vec3&a,const Vec3&b){return {a.y*b.z-a.z*b.y, a.z*b.x-a.x*b.z, a.x*b.y-a.y*b.x};}
static inline double norm2(const Vec3&a){return dot3(a,a);} 
static inline double norm3(const Vec3&a){return sqrt(norm2(a));}
static inline Vec3 unit3(Vec3 a){double n=norm3(a); return n>1e-300 ? a*(1.0/n) : Vec3{0,0,0};}
static inline double clampd(double x,double lo,double hi){return x<lo?lo:(x>hi?hi:x);} 

struct Face { int v[3]; };

static int N0=0, M0=0;
static vector<Vec3> origP, P;
static vector<Face> origF, F;
static vector<unsigned char> aliveV, aliveF;
static vector<double> coverR;
static vector<vector<int>> adj;
static int aliveN=0, aliveM=0;
static Vec3 bbMin, bbMax;
static double diagLen=1.0;
static chrono::steady_clock::time_point T0;
static const double PI = acos(-1.0);

static inline double elapsed(){ return chrono::duration<double>(chrono::steady_clock::now()-T0).count(); }

struct FastInput {
    vector<char> buf; char* p=nullptr;
    FastInput(){
        buf.reserve(1<<27);
        char chunk[1<<16]; size_t n;
        while((n=fread(chunk,1,sizeof(chunk),stdin))>0) buf.insert(buf.end(), chunk, chunk+n);
        buf.push_back('\0'); p=buf.data();
    }
    inline void skip(){ while(*p==' '||*p=='\n'||*p=='\r'||*p=='\t') ++p; }
    long nextLong(){ skip(); return strtol(p,&p,10); }
    double nextDouble(){ skip(); return strtod(p,&p); }
    char nextChar(){ skip(); return *p++; }
};

static inline uint64_t edgeKey(int a,int b){ if(a>b) swap(a,b); return (uint64_t)(uint32_t)a<<32 | (uint32_t)b; }
static inline bool faceHas(const Face&f,int v){ return f.v[0]==v || f.v[1]==v || f.v[2]==v; }
static inline int thirdOf(const Face&f,int a,int b){ for(int k=0;k<3;k++){int x=f.v[k]; if(x!=a && x!=b) return x;} return -1; }
static inline Vec3 faceCrossCur(const Face&f){ return cross3(P[f.v[1]]-P[f.v[0]], P[f.v[2]]-P[f.v[0]]); }
static inline Vec3 faceCrossOrig(const Face&f){ return cross3(origP[f.v[1]]-origP[f.v[0]], origP[f.v[2]]-origP[f.v[0]]); }
static inline array<int,3> sortedFace3(int a,int b,int c){ array<int,3> t{a,b,c}; sort(t.begin(),t.end()); return t; }
static inline bool sameFaceUnordered(const Face&f,int a,int b,int c){ return sortedFace3(f.v[0],f.v[1],f.v[2]) == sortedFace3(a,b,c); }

static void readInput(){
    FastInput in;
    N0=(int)in.nextLong(); M0=(int)in.nextLong();
    origP.resize(N0); P.resize(N0);
    aliveV.assign(N0,1); coverR.assign(N0,0.0);
    bbMin={1e100,1e100,1e100}; bbMax={-1e100,-1e100,-1e100};
    for(int i=0;i<N0;i++){
        (void)in.nextChar();
        Vec3 q{in.nextDouble(), in.nextDouble(), in.nextDouble()};
        origP[i]=P[i]=q;
        bbMin.x=min(bbMin.x,q.x); bbMin.y=min(bbMin.y,q.y); bbMin.z=min(bbMin.z,q.z);
        bbMax.x=max(bbMax.x,q.x); bbMax.y=max(bbMax.y,q.y); bbMax.z=max(bbMax.z,q.z);
    }
    diagLen=norm3(bbMax-bbMin); if(!(diagLen>0)) diagLen=1.0;
    origF.resize(M0); F.resize(M0); aliveF.assign(M0,1);
    vector<int> deg(N0,0);
    for(int i=0;i<M0;i++){
        (void)in.nextChar();
        int a=(int)in.nextLong()-1, b=(int)in.nextLong()-1, c=(int)in.nextLong()-1;
        origF[i].v[0]=F[i].v[0]=a;
        origF[i].v[1]=F[i].v[1]=b;
        origF[i].v[2]=F[i].v[2]=c;
        if(a>=0&&a<N0) deg[a]++; if(b>=0&&b<N0) deg[b]++; if(c>=0&&c<N0) deg[c]++;
    }
    adj.assign(N0,{});
    for(int i=0;i<N0;i++) adj[i].reserve(deg[i]);
    for(int i=0;i<M0;i++){
        adj[F[i].v[0]].push_back(i); adj[F[i].v[1]].push_back(i); adj[F[i].v[2]].push_back(i);
    }
    aliveN=N0; aliveM=M0;
}

static void rebuildAdj(){
    vector<int> deg(N0,0);
    int nf=0;
    for(int i=0;i<(int)F.size();i++){
        if(!aliveF[i]) continue;
        Face f=F[i];
        bool ok=f.v[0]>=0&&f.v[0]<N0&&f.v[1]>=0&&f.v[1]<N0&&f.v[2]>=0&&f.v[2]<N0&&
                f.v[0]!=f.v[1]&&f.v[0]!=f.v[2]&&f.v[1]!=f.v[2]&&
                aliveV[f.v[0]]&&aliveV[f.v[1]]&&aliveV[f.v[2]];
        if(!ok){ aliveF[i]=0; continue; }
        deg[f.v[0]]++; deg[f.v[1]]++; deg[f.v[2]]++; nf++;
    }
    adj.assign(N0,{});
    for(int i=0;i<N0;i++) if(aliveV[i]) adj[i].reserve(deg[i]);
    for(int i=0;i<(int)F.size();i++) if(aliveF[i]){
        Face f=F[i]; adj[f.v[0]].push_back(i); adj[f.v[1]].push_back(i); adj[f.v[2]].push_back(i);
    }
    aliveM=nf;
    int nv=0; for(unsigned char b:aliveV) if(b) nv++; aliveN=nv;
}

struct State {
    vector<Vec3> P;
    vector<Face> F;
    vector<unsigned char> aliveV, aliveF;
    vector<double> coverR;
    vector<vector<int>> adj;
    int aliveN=0, aliveM=0;
};
static State saveState(){ return {P,F,aliveV,aliveF,coverR,adj,aliveN,aliveM}; }
static void restoreState(const State&s){ P=s.P; F=s.F; aliveV=s.aliveV; aliveF=s.aliveF; coverR=s.coverR; adj=s.adj; aliveN=s.aliveN; aliveM=s.aliveM; }

static bool findEdgeFaces(int u,int v,int ef[2],int opp[2]){
    if(u<0||v<0||u>=N0||v>=N0||!aliveV[u]||!aliveV[v]) return false;
    int cnt=0;
    const vector<int>& au = adj[u].size() <= adj[v].size() ? adj[u] : adj[v];
    for(int fid: au){
        if(!aliveF[fid]) continue;
        const Face& f=F[fid];
        if(!faceHas(f,u)||!faceHas(f,v)) continue;
        if(cnt>=2) return false;
        ef[cnt]=fid; opp[cnt]=thirdOf(f,u,v); cnt++;
    }
    if(cnt!=2) return false;
    if(opp[0]<0||opp[1]<0||opp[0]==opp[1]) return false;
    return true;
}

static vector<int> markA, markB; static int stampA=1, stampB=1;
static bool linkOK(int u,int v,const int opp[2]){
    if(markA.empty()){ markA.assign(N0,0); markB.assign(N0,0); }
    if(++stampA>2000000000){ fill(markA.begin(),markA.end(),0); stampA=1; }
    if(++stampB>2000000000){ fill(markB.begin(),markB.end(),0); stampB=1; }
    for(int fid: adj[u]) if(aliveF[fid] && faceHas(F[fid],u)){
        const Face& f=F[fid];
        for(int k=0;k<3;k++){ int x=f.v[k]; if(x!=u && x!=v) markA[x]=stampA; }
    }
    int common=0;
    for(int fid: adj[v]) if(aliveF[fid] && faceHas(F[fid],v)){
        const Face& f=F[fid];
        for(int k=0;k<3;k++){
            int x=f.v[k]; if(x==u||x==v) continue;
            if(markA[x]!=stampA) continue;
            if(x!=opp[0] && x!=opp[1]) return false;
            if(markB[x]!=stampB){ markB[x]=stampB; common++; }
        }
    }
    return common==2 && markB[opp[0]]==stampB && markB[opp[1]]==stampB;
}

static bool duplicateFaceAfter(int fid,int rem,int ef0,int ef1,int a,int b,int c){
    int probe=a;
    if(adj[b].size()<adj[probe].size()) probe=b;
    if(adj[c].size()<adj[probe].size()) probe=c;
    for(int g: adj[probe]){
        if(!aliveF[g]||g==fid||g==ef0||g==ef1) continue;
        if(faceHas(F[g],rem)) continue;
        if(sameFaceUnordered(F[g],a,b,c)) return true;
    }
    return false;
}

struct CollapseParams {
    double epsR=0.0490;
    double planeR=0.010;
    double minCos=0.96;
    double strictPlaneR=0.004;
    double strictCos=0.995;
    double areaR=1e-24;
    bool allowMid=false;
};
struct EvalCollapse {
    bool ok=false; bool strict=false;
    double cost=1e100, newRad=0; int keep=-1, rem=-1; Vec3 pos;
};

static EvalCollapse evalEndpoint(int keep,int rem,const int ef[2],const CollapseParams&par){
    EvalCollapse r; r.keep=keep; r.rem=rem; r.pos=P[keep];
    const double eps=par.epsR*diagLen;
    const double planeTol=par.planeR*diagLen;
    const double strictPlane=par.strictPlaneR*diagLen;
    const double d=norm3(P[keep]-P[rem]);
    r.newRad=max(coverR[keep], coverR[rem]+d);
    if(r.newRad>eps+1e-12) return r;
    const double minArea2=max(1e-300, par.areaR*diagLen*diagLen*diagLen*diagLen);
    double maxPlane=0, maxAng=0, maxLoss=0; int changed=0; bool strict=true;
    for(int fid: adj[rem]){
        if(!aliveF[fid]||!faceHas(F[fid],rem)) continue;
        if(fid==ef[0]||fid==ef[1]) continue;
        if(faceHas(F[fid],keep)) return r;
        Face old=F[fid], nf=old;
        for(int k=0;k<3;k++) if(nf.v[k]==rem) nf.v[k]=keep;
        if(nf.v[0]==nf.v[1]||nf.v[0]==nf.v[2]||nf.v[1]==nf.v[2]) return r;
        Vec3 co=faceCrossCur(old), cn=faceCrossCur(nf);
        double ao=norm3(co), an=norm3(cn);
        if(!(ao>0)||!(an>0)||norm2(cn)<=minArea2) return r;
        double nd=clampd(dot3(co,cn)/(ao*an),-1,1);
        if(nd<par.minCos) return r;
        Vec3 no=co*(1.0/ao);
        double pd=fabs(dot3(no, P[keep]-P[old.v[0]]));
        if(pd>planeTol) return r;
        if(duplicateFaceAfter(fid,rem,ef[0],ef[1],nf.v[0],nf.v[1],nf.v[2])) return r;
        maxPlane=max(maxPlane,pd); maxAng=max(maxAng,1.0-nd); maxLoss=max(maxLoss,max(0.0,1.0-an/ao)); changed++;
        if(pd>strictPlane || nd<par.strictCos) strict=false;
    }
    if(changed==0) return r;
    r.ok=true; r.strict=strict;
    r.cost=(strict?-1000.0:0.0) + 0.9*(r.newRad/(eps+1e-300)) + 0.9*(maxPlane/(planeTol+1e-300)) + 240.0*maxAng + 0.02*maxLoss + 0.0005*changed;
    return r;
}
static EvalCollapse evalMidpoint(int keep,int rem,const int ef[2],const CollapseParams&par,const Vec3&pos){
    EvalCollapse r; r.keep=keep; r.rem=rem; r.pos=pos;
    const double eps=par.epsR*diagLen, planeTol=par.planeR*diagLen;
    r.newRad=max(coverR[keep]+norm3(P[keep]-pos), coverR[rem]+norm3(P[rem]-pos));
    if(r.newRad>eps+1e-12) return r;
    const double minArea2=max(1e-300, par.areaR*diagLen*diagLen*diagLen*diagLen);
    static vector<int> inc; inc.clear(); inc.reserve(adj[keep].size()+adj[rem].size());
    for(int fid: adj[keep]) if(aliveF[fid]&&fid!=ef[0]&&fid!=ef[1]&&faceHas(F[fid],keep)) inc.push_back(fid);
    for(int fid: adj[rem]) if(aliveF[fid]&&fid!=ef[0]&&fid!=ef[1]&&faceHas(F[fid],rem)) inc.push_back(fid);
    sort(inc.begin(),inc.end()); inc.erase(unique(inc.begin(),inc.end()),inc.end());
    if(inc.empty()) return r;
    double maxPlane=0,maxAng=0,maxLoss=0; int changed=0;
    for(int fid: inc){
        Face old=F[fid], nf=old;
        for(int k=0;k<3;k++) if(nf.v[k]==rem) nf.v[k]=keep;
        if(nf.v[0]==nf.v[1]||nf.v[0]==nf.v[2]||nf.v[1]==nf.v[2]) return r;
        Vec3 a = (nf.v[0]==keep?pos:P[nf.v[0]]);
        Vec3 b = (nf.v[1]==keep?pos:P[nf.v[1]]);
        Vec3 c = (nf.v[2]==keep?pos:P[nf.v[2]]);
        Vec3 co=faceCrossCur(old), cn=cross3(b-a,c-a);
        double ao=norm3(co), an=norm3(cn);
        if(!(ao>0)||!(an>0)||norm2(cn)<=minArea2) return r;
        double nd=clampd(dot3(co,cn)/(ao*an),-1,1);
        if(nd<par.minCos) return r;
        Vec3 no=co*(1.0/ao);
        double pd=fabs(dot3(no, pos-P[old.v[0]]));
        if(pd>planeTol) return r;
        if(faceHas(old,rem) && duplicateFaceAfter(fid,rem,ef[0],ef[1],nf.v[0],nf.v[1],nf.v[2])) return r;
        maxPlane=max(maxPlane,pd); maxAng=max(maxAng,1.0-nd); maxLoss=max(maxLoss,max(0.0,1.0-an/ao)); changed++;
    }
    if(changed==0) return r;
    r.ok=true; r.strict=false;
    r.cost=0.35 + 1.05*(r.newRad/(eps+1e-300)) + 1.05*(maxPlane/(planeTol+1e-300)) + 300.0*maxAng + 0.03*maxLoss + 0.0007*changed;
    return r;
}
static void applyCollapse(const EvalCollapse&e,const int ef[2]){
    int keep=e.keep, rem=e.rem;
    for(int i=0;i<2;i++) if(ef[i]>=0 && ef[i]<(int)F.size() && aliveF[ef[i]]){ aliveF[ef[i]]=0; aliveM--; }
    for(int fid: adj[rem]){
        if(!aliveF[fid]||!faceHas(F[fid],rem)) continue;
        for(int k=0;k<3;k++) if(F[fid].v[k]==rem) F[fid].v[k]=keep;
    }
    aliveV[rem]=0; aliveN--;
    P[keep]=e.pos; coverR[keep]=e.newRad;
    vector<int> merged; merged.reserve(adj[keep].size()+adj[rem].size());
    for(int fid: adj[keep]) if(aliveF[fid]&&faceHas(F[fid],keep)) merged.push_back(fid);
    for(int fid: adj[rem]) if(aliveF[fid]&&faceHas(F[fid],keep)) merged.push_back(fid);
    sort(merged.begin(),merged.end()); merged.erase(unique(merged.begin(),merged.end()),merged.end());
    adj[keep].swap(merged); vector<int>().swap(adj[rem]);
}
static bool tryCollapseEdge(int u,int v,const CollapseParams&par){
    if(u==v||u<0||v<0||u>=N0||v>=N0||!aliveV[u]||!aliveV[v]) return false;
    int ef[2]={-1,-1}, opp[2]={-1,-1};
    if(!findEdgeFaces(u,v,ef,opp)) return false;
    if(!linkOK(u,v,opp)) return false;
    EvalCollapse best=evalEndpoint(u,v,ef,par);
    EvalCollapse b=evalEndpoint(v,u,ef,par);
    if(b.ok && (!best.ok || b.cost<best.cost)) best=b;
    if(par.allowMid){
        Vec3 mid=(P[u]+P[v])*0.5;
        EvalCollapse m1=evalMidpoint(u,v,ef,par,mid), m2=evalMidpoint(v,u,ef,par,mid);
        if(m1.ok && (!best.ok || m1.cost<best.cost)) best=m1;
        if(m2.ok && (!best.ok || m2.cost<best.cost)) best=m2;
    }
    if(!best.ok) return false;
    applyCollapse(best,ef); return true;
}

struct EdgeRec { float l2; uint64_t k; };
static void collectEdges(vector<uint64_t>&keys){
    keys.clear(); keys.reserve((size_t)aliveM*3+16);
    for(int i=0;i<(int)F.size();i++) if(aliveF[i]){
        const Face&f=F[i];
        if(f.v[0]<0||f.v[1]<0||f.v[2]<0||f.v[0]>=N0||f.v[1]>=N0||f.v[2]>=N0) continue;
        if(!aliveV[f.v[0]]||!aliveV[f.v[1]]||!aliveV[f.v[2]]) continue;
        if(f.v[0]==f.v[1]||f.v[0]==f.v[2]||f.v[1]==f.v[2]) continue;
        keys.push_back(edgeKey(f.v[0],f.v[1])); keys.push_back(edgeKey(f.v[1],f.v[2])); keys.push_back(edgeKey(f.v[2],f.v[0]));
    }
    sort(keys.begin(),keys.end()); keys.erase(unique(keys.begin(),keys.end()),keys.end());
}
static int edgePass(const CollapseParams&par,int target,double timeLimit,int rounds){
    vector<uint64_t> keys; vector<EdgeRec> edges; int removed=0;
    for(int rd=0; rd<rounds && elapsed()<timeLimit && aliveN>target; rd++){
        collectEdges(keys); if(keys.empty()) break;
        edges.clear(); edges.reserve(keys.size());
        for(uint64_t k: keys){
            int a=(int)(k>>32), b=(int)(uint32_t)k;
            if(a>=0&&b>=0&&a<N0&&b<N0&&aliveV[a]&&aliveV[b]) edges.push_back({(float)norm2(P[a]-P[b]),k});
        }
        sort(edges.begin(),edges.end(),[](const EdgeRec&a,const EdgeRec&b){ return a.l2<b.l2; });
        int hit=0;
        for(size_t i=0;i<edges.size() && elapsed()<timeLimit && aliveN>target;i++){
            int a=(int)(edges[i].k>>32), b=(int)(uint32_t)edges[i].k;
            if(a>=0&&b>=0&&a<N0&&b<N0&&aliveV[a]&&aliveV[b] && tryCollapseEdge(a,b,par)){ hit++; removed++; }
        }
        if(hit==0) break;
        if(rd+1<rounds) rebuildAdj();
    }
    return removed;
}

// ---------- low-resolution visual proxy ----------
struct PixN { float x=0,y=0,z=0; };
struct RenderMap { vector<float> depth; vector<PixN> normal; vector<unsigned char> mask; int res=0; };
static int cacheRes=0; static vector<RenderMap> origCache;

static inline void projectPoint(const Vec3&p,int view,int res,double&u,double&v,double&z){
    const double D=2.5; const double f=800.0*((double)res/1024.0); const double c=0.5*res;
    double sx=0, sy=0;
    if(view==0){ sx=p.y; sy=p.z; z=D-p.x; }
    else if(view==1){ sx=-p.y; sy=p.z; z=D+p.x; }
    else if(view==2){ sx=-p.x; sy=p.z; z=D-p.y; }
    else if(view==3){ sx=p.x; sy=p.z; z=D+p.y; }
    else if(view==4){ sx=p.x; sy=p.y; z=D-p.z; }
    else { sx=-p.x; sy=p.y; z=D+p.z; }
    u=f*sx/z+c; v=f*sy/z+c;
}
static void rasterTri(RenderMap&rm,const Vec3&a,const Vec3&b,const Vec3&c,const Vec3&n,int view){
    int res=rm.res; double x0,y0,z0,x1,y1,z1,x2,y2,z2;
    projectPoint(a,view,res,x0,y0,z0); projectPoint(b,view,res,x1,y1,z1); projectPoint(c,view,res,x2,y2,z2);
    if(z0<=0||z1<=0||z2<=0) return;
    int xmin=max(0,(int)floor(min(x0,min(x1,x2))-0.5));
    int xmax=min(res-1,(int)ceil(max(x0,max(x1,x2))+0.5));
    int ymin=max(0,(int)floor(min(y0,min(y1,y2))-0.5));
    int ymax=min(res-1,(int)ceil(max(y0,max(y1,y2))+0.5));
    if(xmin>xmax||ymin>ymax) return;
    double den=(y1-y2)*(x0-x2)+(x2-x1)*(y0-y2);
    if(fabs(den)<1e-14) return;
    PixN pn{(float)n.x,(float)n.y,(float)n.z};
    for(int yy=ymin; yy<=ymax; yy++){
        double py=yy+0.5;
        for(int xx=xmin; xx<=xmax; xx++){
            double px=xx+0.5;
            double w0=((y1-y2)*(px-x2)+(x2-x1)*(py-y2))/den;
            double w1=((y2-y0)*(px-x2)+(x0-x2)*(py-y2))/den;
            double w2=1.0-w0-w1;
            if(w0<-1e-9||w1<-1e-9||w2<-1e-9) continue;
            double zp=1.0/(w0/z0+w1/z1+w2/z2);
            int idx=yy*res+xx;
            if(zp<rm.depth[idx]){ rm.depth[idx]=(float)zp; rm.normal[idx]=pn; rm.mask[idx]=1; }
        }
    }
}
static RenderMap renderMesh(bool original,int view,int res){
    RenderMap rm; rm.res=res; int n=res*res;
    rm.depth.assign(n,255.0f); rm.normal.assign(n,PixN{}); rm.mask.assign(n,0);
    if(original){
        for(const Face&f: origF){
            Vec3 a=origP[f.v[0]], b=origP[f.v[1]], c=origP[f.v[2]];
            Vec3 cr=cross3(b-a,c-a); double l=norm3(cr); if(l<=0) continue;
            rasterTri(rm,a,b,c,cr*(1.0/l),view);
        }
    }else{
        for(int i=0;i<(int)F.size();i++) if(aliveF[i]){
            const Face&f=F[i]; if(!aliveV[f.v[0]]||!aliveV[f.v[1]]||!aliveV[f.v[2]]) continue;
            Vec3 cr=faceCrossCur(f); double l=norm3(cr); if(l<=0) continue;
            rasterTri(rm,P[f.v[0]],P[f.v[1]],P[f.v[2]],cr*(1.0/l),view);
        }
    }
    return rm;
}
static void ensureOrigCache(int res){
    if(cacheRes==res && (int)origCache.size()==6) return;
    cacheRes=res; origCache.clear(); origCache.reserve(6);
    for(int v=0; v<6; v++) origCache.push_back(renderMesh(true,v,res));
}
static inline double normChannelValue(const PixN&n,int ch){ return ((ch==0?n.x:(ch==1?n.y:n.z))+1.0)*127.5; }
static double ssimChannelNaive(const RenderMap&a,const RenderMap&b,const vector<unsigned char>&fg,int ch,bool depth){
    int res=a.res, rad=5; const double c1=(0.01*255)*(0.01*255), c2=(0.03*255)*(0.03*255);
    double total=0; int cnt=0;
    for(int y=0;y<res;y++) for(int x=0;x<res;x++){
        int center=y*res+x; if(!fg[center]) continue;
        double sx=0,sy=0,sxx=0,syy=0,sxy=0; int n=0;
        for(int dy=-rad;dy<=rad;dy++){
            int yy=y+dy; if(yy<0) yy=0; else if(yy>=res) yy=res-1;
            for(int dx=-rad;dx<=rad;dx++){
                int xx=x+dx; if(xx<0) xx=0; else if(xx>=res) xx=res-1;
                int idx=yy*res+xx; double vx,vy;
                if(depth){ vx=a.depth[idx]; vy=b.depth[idx]; }
                else { vx=normChannelValue(a.normal[idx],ch); vy=normChannelValue(b.normal[idx],ch); }
                sx+=vx; sy+=vy; sxx+=vx*vx; syy+=vy*vy; sxy+=vx*vy; n++;
            }
        }
        double inv=1.0/n, ux=sx*inv, uy=sy*inv;
        double varx=max(0.0,sxx*inv-ux*ux), vary=max(0.0,syy*inv-uy*uy), cov=sxy*inv-ux*uy;
        double num=(2*ux*uy+c1)*(2*cov+c2), den=(ux*ux+uy*uy+c1)*(varx+vary+c2);
        total += den!=0 ? num/den : 1.0; cnt++;
    }
    return cnt? total/cnt : 1.0;
}
static double visualProxyScore(int res){
    if(elapsed()>19.8) return 0.0;
    ensureOrigCache(res);
    double total=0;
    for(int view=0; view<6; view++){
        RenderMap simp=renderMesh(false,view,res);
        vector<unsigned char> fg(res*res,0);
        for(int i=0;i<res*res;i++) fg[i]=(origCache[view].mask[i]||simp.mask[i])?1:0;
        double ns=0; for(int ch=0;ch<3;ch++) ns+=ssimChannelNaive(origCache[view],simp,fg,ch,false); ns/=3.0;
        double ds=ssimChannelNaive(origCache[view],simp,fg,0,true);
        total += 0.5*ns + 0.5*ds;
        if(elapsed()>20.2) break;
    }
    return total/6.0;
}

// ---------- exact-ish vertex Hausdorff check for generated structural candidates ----------
struct PointGrid {
    const vector<Vec3>* pts=nullptr; double cell=1; Vec3 mn{}; unordered_map<long long, vector<int>> mp;
    static long long key3(int x,int y,int z){ const long long O=1048576LL, M=(1LL<<21)-1; return (((long long)(x+O)&M)<<42) ^ (((long long)(y+O)&M)<<21) ^ ((long long)(z+O)&M); }
    void build(const vector<Vec3>&p,double c){
        pts=&p; cell=max(c,1e-9); mn={1e100,1e100,1e100};
        for(const Vec3&q:p){ mn.x=min(mn.x,q.x); mn.y=min(mn.y,q.y); mn.z=min(mn.z,q.z); }
        mn.x-=cell*2; mn.y-=cell*2; mn.z-=cell*2;
        mp.clear(); mp.reserve(p.size()*2+16);
        for(int i=0;i<(int)p.size();i++){
            int ix=(int)floor((p[i].x-mn.x)/cell), iy=(int)floor((p[i].y-mn.y)/cell), iz=(int)floor((p[i].z-mn.z)/cell);
            mp[key3(ix,iy,iz)].push_back(i);
        }
    }
    bool hasNear(const Vec3&q,double eps2) const{
        int ix=(int)floor((q.x-mn.x)/cell), iy=(int)floor((q.y-mn.y)/cell), iz=(int)floor((q.z-mn.z)/cell);
        for(int dz=-1;dz<=1;dz++) for(int dy=-1;dy<=1;dy++) for(int dx=-1;dx<=1;dx++){
            auto it=mp.find(key3(ix+dx,iy+dy,iz+dz)); if(it==mp.end()) continue;
            for(int id: it->second) if(norm2((*pts)[id]-q)<=eps2) return true;
        }
        return false;
    }
};
static bool vertexHausdorffCheck(const vector<Vec3>&cand,double eps){
    if(cand.empty() || (int)cand.size()>N0) return false;
    double eps2=eps*eps;
    PointGrid gc; gc.build(cand,eps);
    for(int i=0;i<N0;i++){
        if((i&8191)==0 && elapsed()>18.6) return false;
        if(!gc.hasNear(origP[i],eps2)) return false;
    }
    PointGrid go; go.build(origP,eps);
    for(const Vec3&q: cand){ if(!go.hasNear(q,eps2)) return false; }
    return true;
}
static bool facesManifoldOK(const vector<Vec3>&verts,const vector<Face>&fs){
    if(verts.empty()||fs.empty()||(int)verts.size()>N0) return false;
    double minA=max(1e-30,1e-24*diagLen*diagLen*diagLen*diagLen);
    vector<uint64_t> edges; edges.reserve(fs.size()*3);
    for(const Face&f: fs){
        for(int k=0;k<3;k++) if(f.v[k]<0||f.v[k]>=(int)verts.size()) return false;
        if(f.v[0]==f.v[1]||f.v[0]==f.v[2]||f.v[1]==f.v[2]) return false;
        Vec3 cr=cross3(verts[f.v[1]]-verts[f.v[0]],verts[f.v[2]]-verts[f.v[0]]);
        if(!(norm2(cr)>minA)) return false;
        edges.push_back(edgeKey(f.v[0],f.v[1])); edges.push_back(edgeKey(f.v[1],f.v[2])); edges.push_back(edgeKey(f.v[2],f.v[0]));
    }
    sort(edges.begin(),edges.end());
    for(size_t i=0;i<edges.size();){ size_t j=i+1; while(j<edges.size()&&edges[j]==edges[i]) j++; if(j-i!=2) return false; i=j; }
    return true;
}
static void setCandidateMesh(const vector<Vec3>&verts,const vector<Face>&fs){
    P.assign(N0,Vec3{}); aliveV.assign(N0,0); coverR.assign(N0,0.0);
    for(int i=0;i<(int)verts.size();i++){ P[i]=verts[i]; aliveV[i]=1; }
    F=fs; aliveF.assign(F.size(),1); aliveN=(int)verts.size(); aliveM=(int)F.size(); rebuildAdj();
}
static bool applyCandidateWithGate(const vector<Vec3>&verts,const vector<Face>&fs,double proxyNeed,int res){
    if((int)verts.size()>=aliveN || (int)verts.size()>N0) return false;
    if(!facesManifoldOK(verts,fs)) return false;
    const double eps=0.0492*diagLen;
    if(!vertexHausdorffCheck(verts,eps)) return false;
    State st=saveState(); setCandidateMesh(verts,fs);
    bool keep=true;
    if(res>0 && elapsed()<18.8){ double q=visualProxyScore(res); keep=q>=proxyNeed; }
    if(!keep){ restoreState(st); return false; }
    return true;
}

// ---------- orientation helpers and structural mesh generation ----------
static double orientationSignFromOriginal(const Vec3&center){
    int stride=max(1,M0/100000); double s=0; int c=0;
    for(int i=0;i<M0;i+=stride){ const Face&f=origF[i]; Vec3 cr=faceCrossOrig(f); Vec3 ce=(origP[f.v[0]]+origP[f.v[1]]+origP[f.v[2]])/3.0; double v=dot3(cr,ce-center); if(fabs(v)>1e-18){ s+=v; c++; } }
    return s>=0?1.0:-1.0;
}
static void orientFace(vector<Vec3>&v,Face&f,const Vec3&center,double sign){
    Vec3 cr=cross3(v[f.v[1]]-v[f.v[0]], v[f.v[2]]-v[f.v[0]]);
    Vec3 ce=(v[f.v[0]]+v[f.v[1]]+v[f.v[2]])/3.0;
    if(dot3(cr,ce-center)*sign<0) swap(f.v[1],f.v[2]);
}
static bool sameOrigFaceUnordered(int fid,int a,int b,int c){ if(fid<0||fid>=M0) return false; return sortedFace3(origF[fid].v[0],origF[fid].v[1],origF[fid].v[2])==sortedFace3(a,b,c); }
static Vec3 origFaceNormalSafe(int fid){ if(fid<0||fid>=M0) return {0,0,0}; return unit3(faceCrossOrig(origF[fid])); }
static void orientToNormal(const vector<Vec3>&v,Face&f,const Vec3&n){
    Vec3 cr=cross3(v[f.v[1]]-v[f.v[0]], v[f.v[2]]-v[f.v[0]]);
    if(dot3(cr,n)<0) swap(f.v[1],f.v[2]);
}

static bool detectOrderedTorus(int&U,int&V){
    if(N0<400 || M0!=2*N0) return false;
    vector<pair<int,int>> cand;
    for(int d=8;(long long)d*d<=N0;d++) if(N0%d==0){ cand.push_back({N0/d,d}); cand.push_back({d,N0/d}); }
    sort(cand.begin(),cand.end(),[](auto&a,auto&b){ return abs(a.first-a.second)<abs(b.first-b.second); });
    for(auto pr:cand){
        U=pr.first; V=pr.second; if(U<8||V<8||U>5000||V>5000) continue;
        int total=0, ok=0, step=max(1,N0/900);
        for(int cell=0; cell<N0 && total<900; cell+=step){
            int i=cell/V, j=cell%V, fid=2*cell; if(fid+1>=M0) break;
            int a=i*V+j, b=((i+1)%U)*V+j, c=((i+1)%U)*V+(j+1)%V, d=i*V+(j+1)%V;
            if(sameOrigFaceUnordered(fid,a,b,c)&&sameOrigFaceUnordered(fid+1,a,c,d)) ok++;
            total++;
        }
        if(total>100 && ok*100>=98*total) return true;
    }
    return false;
}
static bool buildTorusSubsample(int U,int V,int U2,int V2,vector<Vec3>&verts,vector<Face>&fs){
    if(U2<3||V2<3||U2*V2>=N0) return false;
    verts.clear(); fs.clear(); verts.reserve(U2*V2); fs.reserve(2*U2*V2);
    for(int i=0;i<U2;i++){ int oi=(long long)i*U/U2; for(int j=0;j<V2;j++){ int oj=(long long)j*V/V2; verts.push_back(origP[oi*V+oj]); } }
    auto id=[&](int i,int j){ return ((i%U2+U2)%U2)*V2 + ((j%V2+V2)%V2); };
    for(int i=0;i<U2;i++) for(int j=0;j<V2;j++){
        int a=id(i,j), b=id(i+1,j), c=id(i+1,j+1), d=id(i,j+1);
        Face f1{{a,b,c}}, f2{{a,c,d}};
        int oi=(long long)i*U/U2, oj=(long long)j*V/V2;
        orientToNormal(verts,f1,origFaceNormalSafe(2*(oi*V+oj)));
        orientToNormal(verts,f2,origFaceNormalSafe(2*(oi*V+oj)+1));
        fs.push_back(f1); fs.push_back(f2);
    }
    return true;
}
static bool tryOrderedTorus(){
    if(elapsed()>4.5) return false;
    int U=0,V=0; if(!detectOrderedTorus(U,V)) return false;
    State st=saveState();
    vector<pair<int,int>> tries;
    for(double div: {10.0,8.0,6.0,5.0,4.0,3.0,2.5}){
        int u=max(6,(int)ceil(U/div)); int v=max(6,(int)ceil(V/div));
        u=min(u,U); v=min(v,V); tries.push_back({u,v});
    }
    sort(tries.begin(),tries.end(),[](auto&a,auto&b){return a.first*a.second < b.first*b.second;});
    tries.erase(unique(tries.begin(),tries.end()),tries.end());
    for(auto [u2,v2]:tries){
        if(elapsed()>16.5) break; restoreState(st);
        vector<Vec3> vs; vector<Face> fs; if(!buildTorusSubsample(U,V,u2,v2,vs,fs)) continue;
        int res=N0>120000?64:128; double need=N0>120000?0.930:0.940;
        if(applyCandidateWithGate(vs,fs,need,res)) return true;
    }
    restoreState(st); return false;
}

static bool detectOrderedLatLong(int&R,int&V){
    if(N0<300 || M0!=2*(N0-2)) return false;
    int base=N0-2;
    for(int v=8; v<=min(4096,base); v++){
        if(base%v) continue; int r=base/v; if(r<2) continue;
        int step=max(1,v/80); bool ok=true;
        for(int j=0;j<v && ok;j+=step){ int a=1+j,b=1+(j+1)%v; if(!sameOrigFaceUnordered(j,0,b,a)) ok=false; }
        int span=max(1,(r-1)*v/260);
        for(int q=0;q<(r-1)*v && ok;q+=span){
            int rr=q/v,j=q%v; int a=1+rr*v+j,b=1+rr*v+(j+1)%v,c=1+(rr+1)*v+j,d=1+(rr+1)*v+(j+1)%v; int fid=v+2*(rr*v+j);
            if(!sameOrigFaceUnordered(fid,a,b,c)||!sameOrigFaceUnordered(fid+1,b,d,c)) ok=false;
        }
        if(ok){ R=r; V=v; return true; }
    }
    return false;
}
static bool buildLatLongSubsample(int R,int V,int R2,int V2,vector<Vec3>&verts,vector<Face>&fs){
    if(R2<2||V2<8||2+R2*V2>=N0) return false;
    verts.clear(); fs.clear(); verts.reserve(2+R2*V2); fs.reserve(2*R2*V2);
    verts.push_back(origP[0]);
    for(int r=0;r<R2;r++){ int orr=(long long)r*(R-1)/max(1,R2-1); for(int j=0;j<V2;j++){ int oj=(long long)j*V/V2; verts.push_back(origP[1+orr*V+oj]); } }
    int bot=(int)verts.size(); verts.push_back(origP[N0-1]);
    auto rid=[&](int r,int j){ return 1+r*V2+((j%V2+V2)%V2); };
    for(int j=0;j<V2;j++){ Face f{{0,rid(0,j+1),rid(0,j)}}; int oj=(long long)j*V/V2; orientToNormal(verts,f,origFaceNormalSafe(oj)); fs.push_back(f); }
    for(int r=0;r<R2-1;r++) for(int j=0;j<V2;j++){
        int a=rid(r,j), b=rid(r,j+1), c=rid(r+1,j), d=rid(r+1,j+1);
        Face f1{{a,b,c}}, f2{{b,d,c}}; int orr=(long long)r*(R-1)/max(1,R2-1), oj=(long long)j*V/V2; int fid=V+2*(orr*V+oj);
        orientToNormal(verts,f1,origFaceNormalSafe(fid)); orientToNormal(verts,f2,origFaceNormalSafe(fid+1)); fs.push_back(f1); fs.push_back(f2);
    }
    int bottomStart=V+2*(R-1)*V;
    for(int j=0;j<V2;j++){ Face f{{bot,rid(R2-1,j),rid(R2-1,j+1)}}; int oj=(long long)j*V/V2; orientToNormal(verts,f,origFaceNormalSafe(bottomStart+oj)); fs.push_back(f); }
    return true;
}
static bool tryOrderedLatLong(){
    if(elapsed()>5.5) return false;
    int R=0,V=0; if(!detectOrderedLatLong(R,V)) return false;
    State st=saveState(); vector<pair<int,int>> tries;
    for(double div: {10.0,8.0,6.0,5.0,4.0,3.0,2.4}){
        int r=max(2,(int)ceil(R/div)); int v=max(8,(int)ceil(V/div)); r=min(r,R); v=min(v,V); tries.push_back({r,v});
    }
    sort(tries.begin(),tries.end(),[](auto&a,auto&b){return (2+a.first*a.second)<(2+b.first*b.second);}); tries.erase(unique(tries.begin(),tries.end()),tries.end());
    for(auto [r2,v2]:tries){
        if(elapsed()>16.5) break; restoreState(st);
        vector<Vec3> vs; vector<Face> fs; if(!buildLatLongSubsample(R,V,r2,v2,vs,fs)) continue;
        int res=N0>120000?64:128; double need=N0>120000?0.930:0.940;
        if(applyCandidateWithGate(vs,fs,need,res)) return true;
    }
    restoreState(st); return false;
}

// Jacobi eigenvectors for symmetric 3x3. Output unit axes in descending variance order approximately.
static void jacobiAxes(double a[3][3], Vec3 axes[3]){
    double v[3][3]={{1,0,0},{0,1,0},{0,0,1}};
    for(int it=0;it<50;it++){
        int p=0,q=1; double best=fabs(a[0][1]);
        if(fabs(a[0][2])>best){p=0;q=2;best=fabs(a[0][2]);}
        if(fabs(a[1][2])>best){p=1;q=2;best=fabs(a[1][2]);}
        if(best<1e-18) break;
        double app=a[p][p], aqq=a[q][q], apq=a[p][q];
        double tau=(aqq-app)/(2*apq); double t=(tau>=0?1:-1)/(fabs(tau)+sqrt(1+tau*tau));
        double c=1/sqrt(1+t*t), s=t*c;
        for(int k=0;k<3;k++) if(k!=p&&k!=q){ double akp=a[k][p], akq=a[k][q]; a[k][p]=a[p][k]=c*akp-s*akq; a[k][q]=a[q][k]=s*akp+c*akq; }
        a[p][p]=c*c*app-2*s*c*apq+s*s*aqq; a[q][q]=s*s*app+2*s*c*apq+c*c*aqq; a[p][q]=a[q][p]=0;
        for(int k=0;k<3;k++){ double vip=v[k][p], viq=v[k][q]; v[k][p]=c*vip-s*viq; v[k][q]=s*vip+c*viq; }
    }
    array<int,3> idx{0,1,2}; sort(idx.begin(),idx.end(),[&](int i,int j){return a[i][i]>a[j][j];});
    for(int r=0;r<3;r++){ int i=idx[r]; axes[r]=unit3({v[0][i],v[1][i],v[2][i]}); }
}
struct EllFit { bool ok=false; Vec3 center; Vec3 axis[3]; double rad[3]; double rms=1e9, mx=1e9; };
static EllFit fitEllipsoidPCA(){
    EllFit fit; if(N0<800) return fit;
    Vec3 mean{}; int stride=max(1,N0/250000), cnt=0;
    for(int i=0;i<N0;i+=stride){ mean=mean+origP[i]; cnt++; } mean=mean/(double)max(1,cnt);
    double cov[3][3]{};
    for(int i=0;i<N0;i+=stride){ Vec3 q=origP[i]-mean; double x[3]={q.x,q.y,q.z}; for(int a=0;a<3;a++) for(int b=0;b<3;b++) cov[a][b]+=x[a]*x[b]; }
    for(int a=0;a<3;a++) for(int b=0;b<3;b++) cov[a][b]/=max(1,cnt);
    jacobiAxes(cov,fit.axis);
    double lo[3]={1e100,1e100,1e100}, hi[3]={-1e100,-1e100,-1e100};
    for(const Vec3&p: origP) for(int k=0;k<3;k++){ double t=dot3(p,fit.axis[k]); lo[k]=min(lo[k],t); hi[k]=max(hi[k],t); }
    fit.center={0,0,0};
    for(int k=0;k<3;k++){ double mid=0.5*(lo[k]+hi[k]); fit.center=fit.center+fit.axis[k]*mid; fit.rad[k]=0.5*(hi[k]-lo[k]); if(!(fit.rad[k]>1e-8)) return fit; }
    double ss=0, ma=0; int checked=0;
    for(int i=0;i<N0;i+=stride){ Vec3 q=origP[i]-fit.center; double r2=0; for(int k=0;k<3;k++){ double u=dot3(q,fit.axis[k])/fit.rad[k]; r2+=u*u; } double e=fabs(sqrt(max(0.0,r2))-1.0); ss+=e*e; ma=max(ma,e); checked++; }
    fit.rms=sqrt(ss/max(1,checked)); fit.mx=ma;
    double rmax=max(fit.rad[0],max(fit.rad[1],fit.rad[2])), rmin=min(fit.rad[0],min(fit.rad[1],fit.rad[2]));
    if(rmin<rmax*0.12) return fit;
    double rmsLim=N0<5000?0.010:0.018, maxLim=N0<5000?0.040:0.080;
    if(fit.rms<=rmsLim && fit.mx<=maxLim) fit.ok=true;
    return fit;
}
static bool buildEllipsoid(const EllFit&fit,int lat,int lon,vector<Vec3>&verts,vector<Face>&fs){
    if(!fit.ok||lat<4||lon<8) return false; int vc=2+(lat-1)*lon; if(vc>=N0) return false;
    verts.clear(); fs.clear(); verts.reserve(vc); fs.reserve(2*lon*(lat-1));
    auto makep=[&](double x,double y,double z){ Vec3 p=fit.center; p=p+fit.axis[0]*(fit.rad[0]*x)+fit.axis[1]*(fit.rad[1]*y)+fit.axis[2]*(fit.rad[2]*z); return p; };
    verts.push_back(makep(0,0,1)); verts.push_back(makep(0,0,-1));
    auto rid=[&](int r,int j){ return 2+(r-1)*lon+((j%lon+lon)%lon); };
    for(int r=1;r<=lat-1;r++){ double th=PI*r/lat, st=sin(th), ct=cos(th); for(int j=0;j<lon;j++){ double ph=2*PI*j/lon; verts.push_back(makep(st*cos(ph),st*sin(ph),ct)); } }
    double sign=orientationSignFromOriginal(fit.center);
    auto add=[&](int a,int b,int c){ Face f{{a,b,c}}; orientFace(verts,f,fit.center,sign); fs.push_back(f); };
    for(int j=0;j<lon;j++) add(0,rid(1,j),rid(1,j+1));
    for(int r=1;r<=lat-2;r++) for(int j=0;j<lon;j++){ int a=rid(r,j), b=rid(r+1,j), c=rid(r+1,j+1), d=rid(r,j+1); add(a,b,c); add(a,c,d); }
    for(int j=0;j<lon;j++) add(1,rid(lat-1,j+1),rid(lat-1,j));
    return true;
}
static bool tryEllipsoid(){
    if(elapsed()>6.0 || N0<900) return false;
    EllFit fit=fitEllipsoidPCA(); if(!fit.ok) return false;
    double eps=0.0492*diagLen, rmax=max(fit.rad[0],max(fit.rad[1],fit.rad[2]));
    int baseLon=max(12,(int)ceil(2*PI*rmax/(eps*0.82))); int baseLat=max(6,(int)ceil(PI*rmax/(eps*0.82)));
    baseLon=min(baseLon,128); baseLat=min(baseLat,80);
    vector<pair<int,int>> tries;
    for(double mul: {0.70,0.85,1.0,1.18,1.4}) tries.push_back({max(4,(int)ceil(baseLat*mul)), max(8,(int)ceil(baseLon*mul))});
    sort(tries.begin(),tries.end(),[](auto&a,auto&b){return 2+(a.first-1)*a.second < 2+(b.first-1)*b.second;}); tries.erase(unique(tries.begin(),tries.end()),tries.end());
    State st=saveState();
    for(auto [lat,lon]: tries){ if(elapsed()>17.0) break; restoreState(st); vector<Vec3> vs; vector<Face> fs; if(!buildEllipsoid(fit,lat,lon,vs,fs)) continue; int res=N0>150000?64:128; double need=N0>150000?0.928:0.938; if(applyCandidateWithGate(vs,fs,need,res)) return true; }
    restoreState(st); return false;
}

// ---------- one-ring patch removal ----------
static bool vecContains(const vector<int>&v,int x){ return find(v.begin(),v.end(),x)!=v.end(); }
static int edgeCountOutside(int a,int b,const vector<int>&skip){ int c=0; for(int fid: adj[a]) if(aliveF[fid]&&!vecContains(skip,fid)&&faceHas(F[fid],a)&&faceHas(F[fid],b)) if(++c>2) return c; return c; }
static bool faceExistsOutside(const Face&g,const vector<int>&skip){ int a=g.v[0],b=g.v[1],c=g.v[2], probe=a; if(adj[b].size()<adj[probe].size()) probe=b; if(adj[c].size()<adj[probe].size()) probe=c; for(int fid: adj[probe]) if(aliveF[fid]&&!vecContains(skip,fid)&&sameFaceUnordered(F[fid],a,b,c)) return true; return false; }
struct Patch { int v=-1; vector<int> inc, ring; vector<Face> nf; double score=1e100; };
static bool inspectPatch(int v,Patch&pc){
    pc=Patch(); pc.v=v; if(v<0||v>=N0||!aliveV[v]) return false;
    for(int fid: adj[v]) if(aliveF[fid]&&faceHas(F[fid],v)) pc.inc.push_back(fid);
    int k=(int)pc.inc.size(); if(k<3||k>7) return false;
    vector<int> nb; vector<pair<int,int>> bedges; Vec3 nsum{};
    for(int fid: pc.inc){
        const Face&f=F[fid]; int oth[2], t=0;
        for(int i=0;i<3;i++) if(f.v[i]!=v) oth[t++]=f.v[i];
        if(t!=2||oth[0]==oth[1]||!aliveV[oth[0]]||!aliveV[oth[1]]) return false;
        for(int q=0;q<2;q++) if(!vecContains(nb,oth[q])) nb.push_back(oth[q]);
        bedges.push_back({oth[0],oth[1]}); Vec3 cr=faceCrossCur(f); double l=norm3(cr); if(l<=0) return false; nsum=nsum+cr*(1.0/l);
    }
    if((int)nb.size()!=k) return false;
    Vec3 n=unit3(nsum); if(norm2(n)<0.5) return false;
    double oldCos=cos((N0>50000?13.0:10.0)*PI/180.0);
    for(int fid: pc.inc) if(dot3(unit3(faceCrossCur(F[fid])),n)<oldCos) return false;
    auto idx=[&](int x){ for(int i=0;i<k;i++) if(nb[i]==x) return i; return -1; };
    vector<vector<int>> cyc(k); vector<pair<int,int>> boundary;
    for(auto e: bedges){ int a=idx(e.first), b=idx(e.second); if(a<0||b<0||a==b) return false; int aa=a,bb=b; if(aa>bb) swap(aa,bb); if(find(boundary.begin(),boundary.end(),make_pair(aa,bb))!=boundary.end()) return false; boundary.push_back({aa,bb}); cyc[a].push_back(b); cyc[b].push_back(a); }
    for(int i=0;i<k;i++) if(cyc[i].size()!=2) return false;
    vector<int> ord; int cur=0, prev=-1;
    for(int step=0;step<k;step++){ ord.push_back(cur); int nx=(cyc[cur][0]==prev?cyc[cur][1]:cyc[cur][0]); prev=cur; cur=nx; }
    if(cur!=0) return false;
    for(int id: ord) pc.ring.push_back(nb[id]);
    Vec3 ctr=P[v]; for(int u: pc.ring) ctr=ctr+P[u]; ctr=ctr/(double)(k+1);
    double maxPlane=fabs(dot3(P[v]-ctr,n)), nearest=1e100;
    for(int u: pc.ring){ maxPlane=max(maxPlane,fabs(dot3(P[u]-ctr,n))); nearest=min(nearest,norm3(P[u]-P[v])+coverR[v]); }
    // Vertex-only Hausdorff: the removed original cluster must remain near at least one boundary vertex.
    if(nearest>0.0492*diagLen) return false;
    if(maxPlane>(N0>50000?0.008:0.006)*diagLen) return false;
    auto isBoundary=[&](int a,int b){ int ia=idx(a), ib=idx(b); if(ia<0||ib<0) return false; if(ia>ib) swap(ia,ib); return find(boundary.begin(),boundary.end(),make_pair(ia,ib))!=boundary.end(); };
    pc.nf.clear(); pc.nf.reserve(k-2); double newCos=cos((N0>50000?24.0:18.0)*PI/180.0); double minA=max(1e-30,1e-24*diagLen*diagLen*diagLen*diagLen);
    for(int i=1;i+1<k;i++){
        Face g{{pc.ring[0],pc.ring[i],pc.ring[i+1]}}; Vec3 cr=faceCrossCur(g); if(dot3(cr,n)<0){ swap(g.v[1],g.v[2]); cr=cr*(-1); }
        double a2=norm2(cr), al=sqrt(a2); if(!(a2>minA)||dot3(cr*(1.0/al),n)<newCos) return false;
        if(faceExistsOutside(g,pc.inc)) return false;
        int ee[3][2]={{g.v[0],g.v[1]},{g.v[1],g.v[2]},{g.v[2],g.v[0]}};
        for(auto &e: ee){ int a=e[0], b=e[1]; int cnt=edgeCountOutside(a,b,pc.inc); if(isBoundary(a,b)){ if(cnt!=1) return false; } else { if(cnt!=0) return false; } }
        pc.nf.push_back(g);
    }
    pc.score=nearest/diagLen + 0.3*maxPlane/diagLen + 0.002*(k-3);
    return true;
}
static bool applyPatch(const Patch&pc){
    if(!aliveV[pc.v]) return false; Patch cur; if(!inspectPatch(pc.v,cur)) return false; if(cur.ring!=pc.ring) return false;
    for(int fid: cur.inc) if(aliveF[fid]){ aliveF[fid]=0; aliveM--; }
    aliveV[cur.v]=0; aliveN--; adj[cur.v].clear();
    for(const Face&g: cur.nf){ int id=(int)F.size(); F.push_back(g); aliveF.push_back(1); aliveM++; for(int k=0;k<3;k++) adj[g.v[k]].push_back(id); }
    return true;
}
static int patchPass(int cap,double timeLimit){
    vector<Patch> cand; cand.reserve(min(N0,50000));
    for(int v=0; v<N0 && elapsed()<timeLimit-0.25; v++) if(aliveV[v]){ Patch p; if(inspectPatch(v,p)) cand.push_back(std::move(p)); }
    sort(cand.begin(),cand.end(),[](const Patch&a,const Patch&b){ if(a.score!=b.score) return a.score<b.score; return a.v<b.v; });
    int rm=0;
    for(const Patch&p:cand){ if(rm>=cap||elapsed()>timeLimit) break; if(applyPatch(p)) rm++; }
    return rm;
}

// ---------- smoothness and gates ----------
struct SmoothStats { double smooth=0, coarse=0, sharp=0, verySharp=0, bad=0; int samples=0; };
static SmoothStats analyzeSmooth(){
    SmoothStats s; int stride=max(1,(int)F.size()/50000), limit=120000;
    double c10=cos(10*PI/180.0), c30=cos(30*PI/180.0), c22=cos(22*PI/180.0), c45=cos(45*PI/180.0);
    int sm=0,co=0,sh=0,vs=0,bad=0,tot=0;
    for(int fid=0; fid<(int)F.size() && tot<limit; fid+=stride){
        if(!aliveF[fid]) continue; const Face&f=F[fid]; int e[3][2]={{f.v[0],f.v[1]},{f.v[1],f.v[2]},{f.v[2],f.v[0]}};
        for(int k=0;k<3 && tot<limit;k++){ int ef[2],op[2]; if(!findEdgeFaces(e[k][0],e[k][1],ef,op)){ bad++; continue; } Vec3 n0=unit3(faceCrossCur(F[ef[0]])), n1=unit3(faceCrossCur(F[ef[1]])); double d=clampd(dot3(n0,n1),-1,1); tot++; if(d>c10) sm++; if(d>c30) co++; if(d<c22) sh++; if(d<c45) vs++; }
    }
    s.samples=tot; if(tot){ s.smooth=(double)sm/tot; s.coarse=(double)co/tot; s.sharp=(double)sh/tot; s.verySharp=(double)vs/tot; } if(tot+bad) s.bad=(double)bad/(tot+bad); return s;
}
static bool gatedEdgePhase(const CollapseParams&par,int target,double tl,int rounds,int proxyRes,double need){
    if(aliveN<=target || elapsed()>tl-0.5) return false;
    State st=saveState(); int before=aliveN; int rm=edgePass(par,target,tl,rounds); rebuildAdj();
    if(rm<=0 || aliveN>=before){ restoreState(st); return false; }
    bool keep=true;
    if(proxyRes>0 && elapsed()<18.7){ double q=visualProxyScore(proxyRes); keep=q>=need; }
    if(!keep){ restoreState(st); return false; }
    return true;
}
static bool gatedPatchPhase(int cap,double tl,int proxyRes,double need){
    State st=saveState(); int before=aliveN; int rm=patchPass(cap,tl); rebuildAdj();
    if(rm<=0 || aliveN>=before){ restoreState(st); return false; }
    bool keep=true; if(proxyRes>0 && elapsed()<18.8){ double q=visualProxyScore(proxyRes); keep=q>=need; }
    if(!keep){ restoreState(st); return false; } return true;
}

static void genericSimplify(){
    rebuildAdj();
    if(N0<=12){ // sample and tiny meshes: one conservative pass is enough.
        CollapseParams p; p.epsR=0.049; p.planeR=0.02; p.minCos=cos(25*PI/180.0); p.strictPlaneR=0.01; p.strictCos=cos(10*PI/180.0); p.allowMid=false;
        edgePass(p,1,2.0,3); rebuildAdj(); return;
    }
    SmoothStats ss=analyzeSmooth();
    bool verySmooth = ss.samples>200 && ss.smooth>0.985 && ss.sharp<0.015 && ss.bad<0.02;
    bool smooth = ss.samples>200 && ss.coarse>0.92 && ss.verySharp<0.08 && ss.bad<0.03;
    bool sharp = ss.samples>200 && ss.verySharp>0.18;

    int targetSafe = max(20, (int)(N0 * (sharp?0.28:(smooth?0.18:0.23))));
    int targetMed  = max(20, (int)(N0 * (sharp?0.22:(verySmooth?0.10:(smooth?0.13:0.18)))));
    int targetAgg  = max(20, (int)(N0 * (verySmooth?0.075:(smooth?0.105:0.15))));

    CollapseParams safe; safe.epsR=0.0445; safe.planeR=0.0065; safe.minCos=cos((sharp?7.0:10.0)*PI/180.0); safe.strictPlaneR=0.003; safe.strictCos=cos(4*PI/180.0); safe.allowMid=false;
    CollapseParams med; med.epsR=0.0486; med.planeR=smooth?0.012:0.009; med.minCos=cos((smooth?16.0:12.0)*PI/180.0); med.strictPlaneR=0.0055; med.strictCos=cos(7*PI/180.0); med.allowMid=true;
    CollapseParams ag; ag.epsR=0.0492; ag.planeR=smooth?0.018:0.013; ag.minCos=cos((smooth?24.0:16.0)*PI/180.0); ag.strictPlaneR=0.009; ag.strictCos=cos(10*PI/180.0); ag.allowMid=true;

    // Conservative reduction, kept without proxy because it obeys local angle/coverage guards.
    edgePass(safe,targetSafe, min(11.5, N0>250000?10.0:11.5), N0>250000?1:2); rebuildAdj();

    int resMed = N0>180000?64:(N0>70000?96:128);
    int resAgg = N0>180000?64:(N0>70000?96:128);
    double needMed = N0>180000?0.928:(smooth?0.936:0.944);
    double needAgg = N0>180000?(smooth?0.932:0.945):(smooth?0.944:0.955);
    // Medium and aggressive phases are proxy-gated. If the proxy is too slow, the state rolls back.
    gatedEdgePhase(med,targetMed, min(16.0, N0>250000?14.8:16.0), N0>250000?1:2, resMed, needMed);
    if(smooth || verySmooth || N0<80000) gatedEdgePhase(ag,targetAgg, min(18.2, N0>250000?16.8:18.2), 1, resAgg, needAgg);

    int patchCap=max(20,min(N0>250000?9000:5000, aliveN/18));
    if(elapsed()<18.0) gatedPatchPhase(patchCap,19.0,N0>180000?64:128,N0>180000?0.932:0.946);
}

static bool finalStateBasicOK(bool fullEdges){
    if(aliveN<=0||aliveN>N0||aliveM<=0) return false;
    double minA=max(1e-30,1e-24*diagLen*diagLen*diagLen*diagLen);
    vector<uint64_t> edges; if(fullEdges) edges.reserve((size_t)aliveM*3);
    int fc=0;
    for(int i=0;i<(int)F.size();i++) if(aliveF[i]){
        const Face&f=F[i]; for(int k=0;k<3;k++) if(f.v[k]<0||f.v[k]>=N0||!aliveV[f.v[k]]) return false;
        if(f.v[0]==f.v[1]||f.v[0]==f.v[2]||f.v[1]==f.v[2]) return false;
        if(norm2(faceCrossCur(f))<=minA) return false;
        fc++;
        if(fullEdges){ edges.push_back(edgeKey(f.v[0],f.v[1])); edges.push_back(edgeKey(f.v[1],f.v[2])); edges.push_back(edgeKey(f.v[2],f.v[0])); }
    }
    if(fc!=aliveM) aliveM=fc;
    if(fullEdges){ sort(edges.begin(),edges.end()); for(size_t i=0;i<edges.size();){ size_t j=i+1; while(j<edges.size()&&edges[j]==edges[i]) j++; if(j-i!=2) return false; i=j; } }
    return true;
}
static void outputOriginal(){
    string out; out.reserve(1<<20); char line[128];
    auto emit=[&](const char*s,int len){ if(out.size()+len>(1<<20)){ fwrite(out.data(),1,out.size(),stdout); out.clear(); } out.append(s,s+len); };
    int len=snprintf(line,sizeof(line),"%d %d\n",N0,M0); emit(line,len);
    for(const Vec3&p: origP){ len=snprintf(line,sizeof(line),"v %.10g %.10g %.10g\n",p.x,p.y,p.z); emit(line,len); }
    for(const Face&f: origF){ len=snprintf(line,sizeof(line),"f %d %d %d\n",f.v[0]+1,f.v[1]+1,f.v[2]+1); emit(line,len); }
    if(!out.empty()) fwrite(out.data(),1,out.size(),stdout);
}
static void outputState(){
    bool fullEdges = aliveM<350000 && elapsed()<19.2;
    if(!finalStateBasicOK(fullEdges)){ outputOriginal(); return; }
    vector<int> id(N0,-1), old; old.reserve(aliveN);
    for(int i=0;i<N0;i++) if(aliveV[i]){ id[i]=(int)old.size(); old.push_back(i); }
    vector<Face> fs; fs.reserve(aliveM);
    for(int i=0;i<(int)F.size();i++) if(aliveF[i]){
        Face g{{id[F[i].v[0]],id[F[i].v[1]],id[F[i].v[2]]}};
        if(g.v[0]<0||g.v[1]<0||g.v[2]<0) { outputOriginal(); return; }
        fs.push_back(g);
    }
    string out; out.reserve(1<<20); char line[160];
    auto emit=[&](const char*s,int len){ if(out.size()+len>(1<<20)){ fwrite(out.data(),1,out.size(),stdout); out.clear(); } out.append(s,s+len); };
    int len=snprintf(line,sizeof(line),"%d %d\n",(int)old.size(),(int)fs.size()); emit(line,len);
    bool highPrec = (int)old.size()*2 <= N0;
    for(int oi: old){ const Vec3&p=P[oi]; len=snprintf(line,sizeof(line), highPrec?"v %.15g %.15g %.15g\n":"v %.10g %.10g %.10g\n", p.x,p.y,p.z); emit(line,len); }
    for(const Face&f: fs){ len=snprintf(line,sizeof(line),"f %d %d %d\n",f.v[0]+1,f.v[1]+1,f.v[2]+1); emit(line,len); }
    if(!out.empty()) fwrite(out.data(),1,out.size(),stdout);
}

int main(){
    T0=chrono::steady_clock::now();
    readInput();
    // Structural routes first. They are guarded by exact vertex coverage, manifold checks, and a renderer proxy.
    if(N0>300 && elapsed()<1.0){
        if(tryOrderedTorus()){ outputState(); return 0; }
        if(tryOrderedLatLong()){ outputState(); return 0; }
        if(tryEllipsoid()){ outputState(); return 0; }
    }
    genericSimplify();
    outputState();
    return 0;
}
