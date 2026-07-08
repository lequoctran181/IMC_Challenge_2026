#include <bits/stdc++.h>
using namespace std;

struct Vec3{
    double x,y,z;
    Vec3(double X=0,double Y=0,double Z=0):x(X),y(Y),z(Z){}
};
static inline Vec3 operator+(const Vec3&a,const Vec3&b){return Vec3(a.x+b.x,a.y+b.y,a.z+b.z);} 
static inline Vec3 operator-(const Vec3&a,const Vec3&b){return Vec3(a.x-b.x,a.y-b.y,a.z-b.z);} 
static inline Vec3 operator*(const Vec3&a,double s){return Vec3(a.x*s,a.y*s,a.z*s);} 
static inline Vec3 operator/(const Vec3&a,double s){return Vec3(a.x/s,a.y/s,a.z/s);} 
static inline double dotv(const Vec3&a,const Vec3&b){return a.x*b.x+a.y*b.y+a.z*b.z;} 
static inline Vec3 crossv(const Vec3&a,const Vec3&b){return Vec3(a.y*b.z-a.z*b.y,a.z*b.x-a.x*b.z,a.x*b.y-a.y*b.x);} 
static inline double norm2(const Vec3&a){return dotv(a,a);} 
static inline double normv(const Vec3&a){return sqrt(norm2(a));} 
static inline Vec3 normalized(const Vec3&a){ double n=normv(a); return (n>0)?a/n:Vec3(); }
static inline double clampd(double x,double l,double r){ return x<l?l:(x>r?r:x); }

struct Face{ int v[3]; };

struct FastInput{
    vector<char> buf; char *p;
    FastInput(){
        buf.reserve(1<<27);
        char tmp[1<<16]; size_t n;
        while((n=fread(tmp,1,sizeof(tmp),stdin))>0) buf.insert(buf.end(),tmp,tmp+n);
        buf.push_back('\0'); p=buf.data();
    }
    inline void skip(){ while(*p==' '||*p=='\n'||*p=='\r'||*p=='\t') ++p; }
    long nextLong(){ skip(); return strtol(p,&p,10); }
    double nextDouble(){ skip(); return strtod(p,&p); }
    char nextChar(){ skip(); return *p++; }
};

static int N,M;
static vector<Vec3> origV;
static vector<Face> origF;
static double bboxDiag=1.0, epsHaus=0.05;
static Vec3 bboxMin,bboxMax;
static chrono::steady_clock::time_point tStart;
static double elapsed(){ return chrono::duration<double>(chrono::steady_clock::now()-tStart).count(); }

static void readInput(){
    FastInput in;
    N=(int)in.nextLong(); M=(int)in.nextLong();
    origV.resize(N);
    bboxMin=Vec3(1e100,1e100,1e100); bboxMax=Vec3(-1e100,-1e100,-1e100);
    for(int i=0;i<N;i++){
        (void)in.nextChar();
        double x=in.nextDouble(), y=in.nextDouble(), z=in.nextDouble();
        origV[i]=Vec3(x,y,z);
        bboxMin.x=min(bboxMin.x,x); bboxMin.y=min(bboxMin.y,y); bboxMin.z=min(bboxMin.z,z);
        bboxMax.x=max(bboxMax.x,x); bboxMax.y=max(bboxMax.y,y); bboxMax.z=max(bboxMax.z,z);
    }
    Vec3 d=bboxMax-bboxMin;
    bboxDiag=max(1e-30,normv(d));
    epsHaus=0.05*bboxDiag;
    origF.resize(M);
    for(int i=0;i<M;i++){
        (void)in.nextChar();
        origF[i].v[0]=(int)in.nextLong()-1;
        origF[i].v[1]=(int)in.nextLong()-1;
        origF[i].v[2]=(int)in.nextLong()-1;
    }
}

static void writeMesh(const vector<Vec3>&V,const vector<Face>&F){
    string out; out.reserve((size_t)V.size()*42 + (size_t)F.size()*26 + 64);
    char line[128];
    out.append(line, snprintf(line,sizeof(line), "%d %d\n", (int)V.size(), (int)F.size()));
    for(const Vec3&p:V) out.append(line, snprintf(line,sizeof(line), "v %.10g %.10g %.10g\n", p.x,p.y,p.z));
    for(const Face&f:F) out.append(line, snprintf(line,sizeof(line), "f %d %d %d\n", f.v[0]+1,f.v[1]+1,f.v[2]+1));
    fwrite(out.data(),1,out.size(),stdout);
}

// ---------------- spatial validation for the vertex-only Hausdorff constraint ----------------
struct CellKey{ int x,y,z; bool operator==(const CellKey&o)const{return x==o.x&&y==o.y&&z==o.z;} };
struct CellHash{
    static uint64_t splitmix64(uint64_t x){
        x += 0x9e3779b97f4a7c15ULL;
        x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
        x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
        return x ^ (x >> 31);
    }
    size_t operator()(const CellKey&k) const noexcept{
        uint64_t a=(uint32_t)k.x, b=(uint32_t)k.y, c=(uint32_t)k.z;
        return (size_t)(splitmix64(a) ^ (splitmix64(b)<<1) ^ (splitmix64(c)<<2));
    }
};
struct GridIndex{
    double cell=1.0;
    const vector<Vec3>* pts=nullptr;
    unordered_map<CellKey, vector<int>, CellHash> mp;
    GridIndex(){}
    GridIndex(const vector<Vec3>&P,double c){ build(P,c); }
    inline CellKey key(const Vec3&p) const{
        return CellKey{(int)floor(p.x/cell),(int)floor(p.y/cell),(int)floor(p.z/cell)};
    }
    void build(const vector<Vec3>&P,double c){
        cell=max(c,1e-12); pts=&P; mp.clear(); mp.reserve(max<size_t>(1024, P.size()/2));
        for(int i=0;i<(int)P.size();++i) mp[key(P[i])].push_back(i);
    }
    bool hasNear(const Vec3&p,double r) const{
        double r2=r*r; CellKey k=key(p);
        for(int dx=-1;dx<=1;dx++) for(int dy=-1;dy<=1;dy++) for(int dz=-1;dz<=1;dz++){
            CellKey q{k.x+dx,k.y+dy,k.z+dz};
            auto it=mp.find(q); if(it==mp.end()) continue;
            const vector<int>&ids=it->second;
            for(int id:ids) if(norm2((*pts)[id]-p)<=r2) return true;
        }
        return false;
    }
};
static unique_ptr<GridIndex> origGrid;
static void ensureOrigGrid(){ if(!origGrid) origGrid.reset(new GridIndex(origV, max(epsHaus,1e-12))); }

static bool basicMeshOK(const vector<Vec3>&V,const vector<Face>&F){
    if(V.empty() || V.size()>(size_t)N || F.empty()) return false;
    double minA2=max(1e-30, 1e-28*bboxDiag*bboxDiag*bboxDiag*bboxDiag);
    for(const Face&f:F){
        int a=f.v[0],b=f.v[1],c=f.v[2];
        if(a<0||b<0||c<0||a>=(int)V.size()||b>=(int)V.size()||c>=(int)V.size()) return false;
        if(a==b||a==c||b==c) return false;
        if(norm2(crossv(V[b]-V[a], V[c]-V[a]))<=minA2) return false;
    }
    return true;
}
static bool validateHausdorff(const vector<Vec3>&V,const vector<Face>&F){
    if(!basicMeshOK(V,F)) return false;
    ensureOrigGrid();
    const double r=epsHaus*1.0000001;
    for(const Vec3&p:V) if(!origGrid->hasNear(p,r)) return false;
    GridIndex outGrid(V, max(epsHaus,1e-12));
    for(const Vec3&p:origV) if(!outGrid.hasNear(p,r)) return false;
    return true;
}

// ---------------- low-resolution visual proxy (used only for structural replacements) --------
struct RenderMap{
    int res=0;
    vector<float> dep,nx,ny,nz;
    vector<unsigned char> mask;
    RenderMap(){}
    RenderMap(int r):res(r),dep(r*r,255.0f),nx(r*r,0),ny(r*r,0),nz(r*r,0),mask(r*r,0){}
};
static inline bool projectPoint(const Vec3&p,int view,int res,double&u,double&v,double&z){
    constexpr double D=2.5;
    double sx,sy;
    if(view==0){ sx=p.y;  sy=p.z; z=D-p.x; }
    else if(view==1){ sx=-p.y; sy=p.z; z=D+p.x; }
    else if(view==2){ sx=-p.x; sy=p.z; z=D-p.y; }
    else if(view==3){ sx=p.x;  sy=p.z; z=D+p.y; }
    else if(view==4){ sx=p.x;  sy=p.y; z=D-p.z; }
    else { sx=-p.x; sy=p.y; z=D+p.z; }
    if(z<=1e-9) return false;
    double f=800.0*((double)res/1024.0), c=0.5*res;
    u=f*sx/z+c; v=f*sy/z+c; return true;
}
static void rasterTri(RenderMap&rm,const Vec3&a,const Vec3&b,const Vec3&c,const Vec3&n,int view){
    double x0,y0,z0,x1,y1,z1,x2,y2,z2;
    int res=rm.res;
    if(!projectPoint(a,view,res,x0,y0,z0) || !projectPoint(b,view,res,x1,y1,z1) || !projectPoint(c,view,res,x2,y2,z2)) return;
    int xmin=max(0,(int)floor(min(x0,min(x1,x2))-0.5));
    int xmax=min(res-1,(int)ceil(max(x0,max(x1,x2))+0.5));
    int ymin=max(0,(int)floor(min(y0,min(y1,y2))-0.5));
    int ymax=min(res-1,(int)ceil(max(y0,max(y1,y2))+0.5));
    if(xmin>xmax||ymin>ymax) return;
    double den=(y1-y2)*(x0-x2)+(x2-x1)*(y0-y2);
    if(fabs(den)<1e-14) return;
    for(int yy=ymin;yy<=ymax;yy++){
        double py=yy+0.5;
        for(int xx=xmin;xx<=xmax;xx++){
            double px=xx+0.5;
            double w0=((y1-y2)*(px-x2)+(x2-x1)*(py-y2))/den;
            double w1=((y2-y0)*(px-x2)+(x0-x2)*(py-y2))/den;
            double w2=1.0-w0-w1;
            if(w0<-1e-9||w1<-1e-9||w2<-1e-9) continue;
            double zp=1.0/(w0/z0+w1/z1+w2/z2);
            int id=yy*res+xx;
            if(zp < rm.dep[id]){
                rm.dep[id]=(float)zp; rm.nx[id]=(float)n.x; rm.ny[id]=(float)n.y; rm.nz[id]=(float)n.z; rm.mask[id]=1;
            }
        }
    }
}
static vector<RenderMap> renderMesh(const vector<Vec3>&V,const vector<Face>&F,int res){
    vector<RenderMap> maps; maps.reserve(6); for(int i=0;i<6;i++) maps.emplace_back(res);
    for(const Face&f:F){
        const Vec3&a=V[f.v[0]], &b=V[f.v[1]], &c=V[f.v[2]];
        Vec3 cr=crossv(b-a,c-a); double len=normv(cr); if(!(len>0)) continue;
        Vec3 n=cr/len;
        for(int view=0;view<6;view++) rasterTri(maps[view],a,b,c,n,view);
    }
    return maps;
}
static vector<RenderMap> origRenderCache;
static int origRenderRes=0;
static const vector<RenderMap>& getOrigRender(int res){
    if(origRenderRes!=res){ origRenderCache=renderMesh(origV,origF,res); origRenderRes=res; }
    return origRenderCache;
}
static inline double valNormal(const RenderMap&m,int idx,int ch){
    if(!m.mask[idx]) return 127.5;
    double v=(ch==0?m.nx[idx]:(ch==1?m.ny[idx]:m.nz[idx]));
    return (v+1.0)*127.5;
}
static inline double valDepth(const RenderMap&m,int idx){ return m.mask[idx]?m.dep[idx]:255.0; }
static double ssimChannel(const RenderMap&A,const RenderMap&B,int ch,bool depth){
    int res=A.res, rad=5; double c1=(0.01*255.0)*(0.01*255.0), c2=(0.03*255.0)*(0.03*255.0);
    double total=0; int cnt=0;
    for(int y=0;y<res;y++) for(int x=0;x<res;x++){
        int center=y*res+x; if(!A.mask[center] && !B.mask[center]) continue;
        double sx=0,sy=0,sxx=0,syy=0,sxy=0; int n=0;
        for(int dy=-rad;dy<=rad;dy++){
            int yy=min(res-1,max(0,y+dy));
            for(int dx=-rad;dx<=rad;dx++){
                int xx=min(res-1,max(0,x+dx)); int id=yy*res+xx;
                double vx=depth?valDepth(A,id):valNormal(A,id,ch);
                double vy=depth?valDepth(B,id):valNormal(B,id,ch);
                sx+=vx; sy+=vy; sxx+=vx*vx; syy+=vy*vy; sxy+=vx*vy; n++;
            }
        }
        double inv=1.0/n, ux=sx*inv, uy=sy*inv;
        double varx=max(0.0,sxx*inv-ux*ux), vary=max(0.0,syy*inv-uy*uy), cov=sxy*inv-ux*uy;
        double num=(2*ux*uy+c1)*(2*cov+c2);
        double den=(ux*ux+uy*uy+c1)*(varx+vary+c2);
        total += (den!=0?num/den:1.0); cnt++;
    }
    return cnt?total/cnt:1.0;
}
static double visualProxy(const vector<Vec3>&V,const vector<Face>&F,int res){
    if(elapsed()>15.8) return 0.0; // avoid accepting an unvalidated structural candidate near the time limit
    const vector<RenderMap>&A=getOrigRender(res);
    vector<RenderMap>B=renderMesh(V,F,res);
    double tot=0;
    for(int view=0;view<6;view++){
        double ns=(ssimChannel(A[view],B[view],0,false)+ssimChannel(A[view],B[view],1,false)+ssimChannel(A[view],B[view],2,false))/3.0;
        double ds=ssimChannel(A[view],B[view],0,true);
        tot += 0.5*ns+0.5*ds;
    }
    return tot/6.0;
}




static void flipFaces(vector<Face>&F){ for(Face&f:F) swap(f.v[1],f.v[2]); }
static double bestVisualProxy(vector<Vec3>&V, vector<Face>&F, int res){
    double a=visualProxy(V,F,res);
    flipFaces(F);
    double b=visualProxy(V,F,res);
    if(b>=a) return b;
    flipFaces(F);
    return a;
}
static void addTriConvex(vector<Vec3>&V, vector<Face>&F, int a,int b,int c,const Vec3&center);
// ---------------- topology-grid structural shortcuts (UV sphere / torus) ----------------
struct U64Hash{
    static uint64_t splitmix64(uint64_t x){
        x += 0x9e3779b97f4a7c15ULL;
        x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
        x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
        return x ^ (x >> 31);
    }
    size_t operator()(uint64_t x) const noexcept { return (size_t)splitmix64(x); }
};
static unique_ptr<unordered_set<uint64_t,U64Hash>> triSet;
static inline uint64_t triKey(int a,int b,int c){
    if(a>b) swap(a,b); if(b>c) swap(b,c); if(a>b) swap(a,b);
    return ((uint64_t)(uint32_t)a<<42) | ((uint64_t)(uint32_t)b<<21) | (uint32_t)c;
}
static void ensureTriSet(){
    if(triSet) return;
    triSet.reset(new unordered_set<uint64_t,U64Hash>());
    triSet->reserve((size_t)M*2+10);
    for(const Face&f:origF) triSet->insert(triKey(f.v[0],f.v[1],f.v[2]));
}
static inline bool hasTri(int a,int b,int c){ return triSet && triSet->find(triKey(a,b,c))!=triSet->end(); }
struct SphereGrid{ int R=0,C=0,layout=0; bool ok=false; };
static bool sphereIdx(const SphereGrid&g,int r,int j,int&out){
    // r=1..R, j=0..C-1. layout 0: top=0, rings begin at 1, bottom=N-1.
    // layout 1: top=0, bottom=1, rings begin at 2.
    j=(j%g.C+g.C)%g.C;
    if(g.layout==0) out=1+(r-1)*g.C+j;
    else out=2+(r-1)*g.C+j;
    return out>=0 && out<N;
}
static bool checkSphereCandidate(int R,int C,int layout){
    if(R<3||C<8) return false;
    SphereGrid g; g.R=R; g.C=C; g.layout=layout;
    int top=0, bottom=(layout==0?N-1:1);
    if(bottom<=0||bottom>=N) return false;
    int stepR=max(1,R/18), stepC=max(1,C/36);
    for(int j=0;j<C;j+=stepC){ int a,b; sphereIdx(g,1,j,a); sphereIdx(g,1,j+1,b); if(!hasTri(top,a,b)) return false; }
    for(int j=0;j<C;j+=stepC){ int a,b; sphereIdx(g,R,j,a); sphereIdx(g,R,j+1,b); if(!hasTri(bottom,a,b)) return false; }
    for(int r=1;r<R;r+=stepR) for(int j=0;j<C;j+=stepC){
        int a,b,c,d; sphereIdx(g,r,j,a); sphereIdx(g,r,j+1,b); sphereIdx(g,r+1,j+1,c); sphereIdx(g,r+1,j,d);
        bool ok=(hasTri(a,b,c)&&hasTri(a,c,d)) || (hasTri(a,b,d)&&hasTri(b,c,d));
        if(!ok) return false;
    }
    return true;
}
static SphereGrid detectSphereGrid(){
    SphereGrid ans; if(N<300 || M!=2*(N-2)) return ans;
    ensureTriSet();
    int total=N-2;
    for(int C=8; C<=min(2048,total); ++C){
        if(total%C) continue; int R=total/C; if(R<3) continue;
        for(int layout=0; layout<2; ++layout){
            if(checkSphereCandidate(R,C,layout)){ ans.ok=true; ans.R=R; ans.C=C; ans.layout=layout; return ans; }
        }
    }
    return ans;
}
static bool buildSphereGrid(const SphereGrid&g,int R2,int C2,vector<Vec3>&V,vector<Face>&F){
    if(!g.ok||R2<2||C2<6) return false;
    int vc=2+R2*C2; if(vc>N) return false;
    int top=0,bottom=(g.layout==0?N-1:1); Vec3 center=(bboxMin+bboxMax)/2.0;
    V.clear(); F.clear(); V.reserve(vc); F.reserve(2*C2*(R2+1));
    V.push_back(origV[top]); V.push_back(origV[bottom]);
    vector<int> oldRing(R2), oldCol(C2);
    for(int r=0;r<R2;r++) oldRing[r]=min(g.R,max(1,1+(int)floor(((double)r+0.5)*g.R/R2)));
    for(int j=0;j<C2;j++) oldCol[j]=min(g.C-1,max(0,(int)floor((double)j*g.C/C2)));
    auto id=[&](int r,int j){ j=(j%C2+C2)%C2; return 2+r*C2+j; };
    for(int r=0;r<R2;r++) for(int j=0;j<C2;j++){ int oi; sphereIdx(g,oldRing[r],oldCol[j],oi); V.push_back(origV[oi]); }
    for(int j=0;j<C2;j++) addTriConvex(V,F,0,id(0,j),id(0,j+1),center);
    for(int r=0;r<R2-1;r++) for(int j=0;j<C2;j++){ int a=id(r,j), b=id(r,j+1), c=id(r+1,j+1), d=id(r+1,j); addTriConvex(V,F,a,b,c,center); addTriConvex(V,F,a,c,d,center); }
    for(int j=0;j<C2;j++) addTriConvex(V,F,1,id(R2-1,j+1),id(R2-1,j),center);
    return true;
}
static bool trySphereGrid(vector<Vec3>&outV,vector<Face>&outF){
    if(elapsed()>12.0) return false;
    SphereGrid g=detectSphereGrid(); if(!g.ok) return false;
    struct T{int r,c,res; double th;};
    vector<T> trials;
    auto add=[&](int r,int c,int res,double th){ r=max(2,r); c=max(6,c); T t{r,c,res,th}; if(find_if(trials.begin(),trials.end(),[&](const T&x){return x.r==r&&x.c==c;})==trials.end()) trials.push_back(t); };
    add((g.R+9)/10,(g.C+9)/10,384,0.900);
    add((g.R+7)/8,(g.C+7)/8,384,0.900);
    add((g.R+5)/6,(g.C+5)/6,256,0.905);
    add((g.R+4)/5,(g.C+4)/5,256,0.910);
    add((g.R+3)/4,(g.C+3)/4,192,0.915);
    add((g.R+2)/3,(g.C+2)/3,160,0.880);
    add((2*g.R+4)/5,(2*g.C+4)/5,192,0.900);
    add((g.R+1)/2,(g.C+1)/2,160,0.910);
    for(const T&t:trials){
        if(elapsed()>16.5) break;
        vector<Vec3> V; vector<Face> F; if(!buildSphereGrid(g,t.r,t.c,V,F)) continue;
        if((int)V.size()>=N) continue;
        if(!validateHausdorff(V,F)) continue;
        double p=bestVisualProxy(V,F,(N<20000?t.res:min(t.res,160)));
        if(p>=t.th){ outV.swap(V); outF.swap(F); return true; }
    }
    return false;
}

struct TorusGrid{ int U=0,V=0; bool ok=false; };
static inline int torIdx(const TorusGrid&g,int i,int j){ i=(i%g.U+g.U)%g.U; j=(j%g.V+g.V)%g.V; return i*g.V+j; }
static bool checkTorusCandidate(int U,int Vv){
    if(U<3||Vv<3) return false; TorusGrid g; g.U=U; g.V=Vv;
    int stepU=max(1,U/24), stepV=max(1,Vv/24);
    for(int i=0;i<U;i+=stepU) for(int j=0;j<Vv;j+=stepV){
        int a=torIdx(g,i,j), b=torIdx(g,i,j+1), c=torIdx(g,i+1,j+1), d=torIdx(g,i+1,j);
        bool ok=(hasTri(a,b,c)&&hasTri(a,c,d)) || (hasTri(a,b,d)&&hasTri(b,c,d));
        if(!ok) return false;
    }
    return true;
}
static TorusGrid detectTorusGrid(){
    TorusGrid ans; if(N<300 || M!=2*N) return ans;
    ensureTriSet();
    for(int Vv=6; Vv<=min(2048,N); ++Vv){
        if(N%Vv) continue; int U=N/Vv; if(U<3) continue;
        if(checkTorusCandidate(U,Vv)){ ans.ok=true; ans.U=U; ans.V=Vv; return ans; }
    }
    return ans;
}
static bool buildTorusGrid(const TorusGrid&g,int U2,int V2,vector<Vec3>&V,vector<Face>&F){
    if(!g.ok||U2<3||V2<3||U2*V2>N) return false;
    V.clear(); F.clear(); V.reserve(U2*V2); F.reserve(2*U2*V2);
    vector<int> oi(U2), oj(V2);
    for(int i=0;i<U2;i++) oi[i]=min(g.U-1,max(0,(int)floor((double)i*g.U/U2)));
    for(int j=0;j<V2;j++) oj[j]=min(g.V-1,max(0,(int)floor((double)j*g.V/V2)));
    auto id=[&](int i,int j){ i=(i%U2+U2)%U2; j=(j%V2+V2)%V2; return i*V2+j; };
    for(int i=0;i<U2;i++) for(int j=0;j<V2;j++) V.push_back(origV[torIdx(g,oi[i],oj[j])]);
    // Parametric order normally preserves orientation. If it does not, the proxy rejects it.
    for(int i=0;i<U2;i++) for(int j=0;j<V2;j++){ int a=id(i,j), b=id(i+1,j), c=id(i+1,j+1), d=id(i,j+1); F.push_back(Face{{a,b,c}}); F.push_back(Face{{a,c,d}}); }
    return true;
}
static bool tryTorusGrid(vector<Vec3>&outV,vector<Face>&outF){
    if(elapsed()>12.0) return false;
    TorusGrid g=detectTorusGrid(); if(!g.ok) return false;
    struct T{int u,v,res; double th;}; vector<T> trials;
    auto add=[&](int u,int v,int res,double th){ u=max(3,u); v=max(3,v); T t{u,v,res,th}; if(find_if(trials.begin(),trials.end(),[&](const T&x){return x.u==u&&x.v==v;})==trials.end()) trials.push_back(t); };
    add((g.U+9)/10,(g.V+9)/10,384,0.900);
    add((g.U+7)/8,(g.V+7)/8,384,0.900);
    add((g.U+5)/6,(g.V+5)/6,256,1.100);
    add((g.U+4)/5,(g.V+4)/5,256,1.100);
    add((g.U+3)/4,(g.V+3)/4,192,1.100);
    add((g.U+2)/3,(g.V+2)/3,160,1.100);
    add((2*g.U+4)/5,(2*g.V+4)/5,192,1.100);
    add((g.U+1)/2,(g.V+1)/2,512,0.905);
    for(const T&t:trials){
        if(elapsed()>16.5) break;
        vector<Vec3> V; vector<Face> F; if(!buildTorusGrid(g,t.u,t.v,V,F)) continue;
        if((int)V.size()>=N) continue;
        if(!validateHausdorff(V,F)) continue;
        double p=bestVisualProxy(V,F,(N<20000?t.res:min(t.res,160)));
        if(p>=t.th){ outV.swap(V); outF.swap(F); return true; }
    }
    return false;
}

// ---------------- geometry helpers for structural templates ----------------
static void orientConvex(vector<Vec3>&V, Face&f, const Vec3&center){
    Vec3 cr=crossv(V[f.v[1]]-V[f.v[0]], V[f.v[2]]-V[f.v[0]]);
    Vec3 ctr=(V[f.v[0]]+V[f.v[1]]+V[f.v[2]])/3.0;
    if(dotv(cr,ctr-center)<0) swap(f.v[1],f.v[2]);
}
static void addTriConvex(vector<Vec3>&V, vector<Face>&F, int a,int b,int c,const Vec3&center){
    Face f{{a,b,c}}; orientConvex(V,f,center); F.push_back(f);
}

static bool tryBox(vector<Vec3>&outV, vector<Face>&outF){
    Vec3 mn=bboxMin,mx=bboxMax, len=mx-mn;
    if(len.x<=1e-12||len.y<=1e-12||len.z<=1e-12) return false;
    double maxd=0, sum2=0; int sample=0;
    int stride=max(1,N/250000);
    for(int i=0;i<N;i+=stride){
        const Vec3&p=origV[i];
        double d=min({fabs(p.x-mn.x),fabs(p.x-mx.x),fabs(p.y-mn.y),fabs(p.y-mx.y),fabs(p.z-mn.z),fabs(p.z-mx.z)});
        maxd=max(maxd,d); sum2+=d*d; sample++;
    }
    double rms=sqrt(sum2/max(1,sample));
    if(maxd>epsHaus*0.55 || rms>epsHaus*0.08) return false;
    int reqx=max(1,(int)ceil(len.x/(epsHaus*1.34)));
    int reqy=max(1,(int)ceil(len.y/(epsHaus*1.34)));
    int reqz=max(1,(int)ceil(len.z/(epsHaus*1.34)));
    vector<array<int,3>> trials;
    auto addTrial=[&](int a,int b,int c){ a=max(1,a);b=max(1,b);c=max(1,c); array<int,3> t{a,b,c}; if(find(trials.begin(),trials.end(),t)==trials.end()) trials.push_back(t); };
    addTrial(1,1,1);
    addTrial((reqx+3)/4,(reqy+3)/4,(reqz+3)/4);
    addTrial((reqx+1)/2,(reqy+1)/2,(reqz+1)/2);
    addTrial(reqx,reqy,reqz);
    addTrial((int)ceil(reqx*1.25),(int)ceil(reqy*1.25),(int)ceil(reqz*1.25));
    Vec3 center=(mn+mx)/2.0;
    for(auto tr:trials){
        if(elapsed()>13.0) break;
        int nx=tr[0],ny=tr[1],nz=tr[2];
        vector<Vec3> V; vector<Face> F;
        V.reserve(2*((nx+1)*(ny+1)+(nx+1)*(nz+1)+(ny+1)*(nz+1)));
        unordered_map<long long,int> id; id.reserve(V.capacity()*2+10);
        auto key=[&](int i,int j,int k)->long long{ return ((long long)i*(ny+1)+j)*(nz+1)+k; };
        auto get=[&](int i,int j,int k)->int{
            long long kk=key(i,j,k); auto it=id.find(kk); if(it!=id.end()) return it->second;
            double x=mn.x+len.x*(double)i/(double)nx;
            double y=mn.y+len.y*(double)j/(double)ny;
            double z=mn.z+len.z*(double)k/(double)nz;
            int idx=(int)V.size(); id[kk]=idx; V.emplace_back(x,y,z); return idx;
        };
        auto quad=[&](int a,int b,int c,int d){ addTriConvex(V,F,a,b,c,center); addTriConvex(V,F,a,c,d,center); };
        for(int i=0;i<nx;i++) for(int j=0;j<ny;j++){ quad(get(i,j,0),get(i+1,j,0),get(i+1,j+1,0),get(i,j+1,0)); quad(get(i,j,nz),get(i,j+1,nz),get(i+1,j+1,nz),get(i+1,j,nz)); }
        for(int i=0;i<nx;i++) for(int k=0;k<nz;k++){ quad(get(i,0,k),get(i,0,k+1),get(i+1,0,k+1),get(i+1,0,k)); quad(get(i,ny,k),get(i+1,ny,k),get(i+1,ny,k+1),get(i,ny,k+1)); }
        for(int j=0;j<ny;j++) for(int k=0;k<nz;k++){ quad(get(0,j,k),get(0,j+1,k),get(0,j+1,k+1),get(0,j,k+1)); quad(get(nx,j,k),get(nx,j,k+1),get(nx,j+1,k+1),get(nx,j+1,k)); }
        if((int)V.size()>N) continue;
        if(validateHausdorff(V,F)){ outV.swap(V); outF.swap(F); return true; }
    }
    return false;
}

static void orthonormalBasis(const Vec3&w, Vec3&u, Vec3&v){
    Vec3 a = fabs(w.x)<0.7 ? Vec3(1,0,0) : Vec3(0,1,0);
    u=normalized(crossv(a,w)); v=normalized(crossv(w,u));
}
static void jacobiEigen(double A[3][3], Vec3 evec[3], double eval[3]){
    double V[3][3]={{1,0,0},{0,1,0},{0,0,1}};
    for(int it=0;it<60;it++){
        int p=0,q=1; double best=fabs(A[0][1]);
        if(fabs(A[0][2])>best){best=fabs(A[0][2]);p=0;q=2;}
        if(fabs(A[1][2])>best){best=fabs(A[1][2]);p=1;q=2;}
        if(best<1e-18) break;
        double app=A[p][p], aqq=A[q][q], apq=A[p][q];
        double tau=(aqq-app)/(2*apq);
        double t=(tau>=0?1.0:-1.0)/(fabs(tau)+sqrt(1+tau*tau));
        double c=1.0/sqrt(1+t*t), s=t*c;
        for(int k=0;k<3;k++) if(k!=p&&k!=q){
            double akp=A[k][p], akq=A[k][q];
            A[k][p]=A[p][k]=c*akp-s*akq;
            A[k][q]=A[q][k]=s*akp+c*akq;
        }
        A[p][p]=c*c*app-2*s*c*apq+s*s*aqq;
        A[q][q]=s*s*app+2*s*c*apq+c*c*aqq;
        A[p][q]=A[q][p]=0;
        for(int k=0;k<3;k++){ double vkp=V[k][p], vkq=V[k][q]; V[k][p]=c*vkp-s*vkq; V[k][q]=s*vkp+c*vkq; }
    }
    for(int i=0;i<3;i++){ eval[i]=A[i][i]; evec[i]=normalized(Vec3(V[0][i],V[1][i],V[2][i])); }
    array<int,3> ord={0,1,2}; sort(ord.begin(),ord.end(),[&](int a,int b){return eval[a]>eval[b];});
    Vec3 te[3]; double tv[3]; for(int i=0;i<3;i++){te[i]=evec[ord[i]]; tv[i]=eval[ord[i]];} for(int i=0;i<3;i++){evec[i]=te[i]; eval[i]=tv[i];}
    Vec3 c=crossv(evec[0],evec[1]); if(dotv(c,evec[2])<0) evec[2]=evec[2]*-1.0;
}
static void computePCA(Vec3 axes[3]){
    Vec3 mean; for(const Vec3&p:origV) mean=mean+p; mean=mean/(double)N;
    double C[3][3]={{0,0,0},{0,0,0},{0,0,0}};
    for(const Vec3&p:origV){ Vec3 q=p-mean; double a[3]={q.x,q.y,q.z}; for(int i=0;i<3;i++) for(int j=0;j<3;j++) C[i][j]+=a[i]*a[j]; }
    for(int i=0;i<3;i++) for(int j=0;j<3;j++) C[i][j]/=max(1,N);
    double ev[3]; jacobiEigen(C,axes,ev);
}

struct EllFit{ bool ok=false; Vec3 center; Vec3 ax[3]; double r[3]; double rms=1e9,maxerr=1e9; };
static EllFit fitEllipsoid(const Vec3 axesIn[3]){
    EllFit fit; for(int k=0;k<3;k++) fit.ax[k]=axesIn[k];
    double mn[3]={1e100,1e100,1e100}, mx[3]={-1e100,-1e100,-1e100};
    for(const Vec3&p:origV) for(int k=0;k<3;k++){ double t=dotv(p,fit.ax[k]); mn[k]=min(mn[k],t); mx[k]=max(mx[k],t); }
    fit.center=Vec3();
    for(int k=0;k<3;k++){ fit.r[k]=0.5*(mx[k]-mn[k]); fit.center=fit.center+fit.ax[k]*(0.5*(mx[k]+mn[k])); }
    if(fit.r[0]<=1e-10||fit.r[1]<=1e-10||fit.r[2]<=1e-10) return fit;
    if(min({fit.r[0],fit.r[1],fit.r[2]}) < max({fit.r[0],fit.r[1],fit.r[2]})*0.08) return fit;
    int stride=max(1,N/250000); double s2=0, me=0; int cnt=0;
    for(int i=0;i<N;i+=stride){
        Vec3 q=origV[i]-fit.center; double val=0;
        for(int k=0;k<3;k++){ double x=dotv(q,fit.ax[k])/fit.r[k]; val+=x*x; }
        double e=fabs(sqrt(max(0.0,val))-1.0); s2+=e*e; me=max(me,e); cnt++;
    }
    fit.rms=sqrt(s2/max(1,cnt)); fit.maxerr=me;
    fit.ok = (fit.rms<0.030 && fit.maxerr<0.16);
    return fit;
}
static bool buildEllipsoid(const EllFit&fit,int lat,int lon,vector<Vec3>&V,vector<Face>&F){
    if(lat<4||lon<8) return false;
    int vc=2+(lat-1)*lon; if(vc>N) return false;
    V.clear(); F.clear(); V.reserve(vc); F.reserve(2*lon*(lat-1));
    const double pi=acos(-1.0);
    auto point=[&](double x,double y,double z){ return fit.center + fit.ax[0]*(fit.r[0]*x) + fit.ax[1]*(fit.r[1]*y) + fit.ax[2]*(fit.r[2]*z); };
    V.push_back(point(0,0,1)); V.push_back(point(0,0,-1));
    auto id=[&](int ring,int j){ j=(j%lon+lon)%lon; return 2+(ring-1)*lon+j; };
    for(int r=1;r<=lat-1;r++){
        double th=pi*(double)r/(double)lat, st=sin(th), ct=cos(th);
        for(int j=0;j<lon;j++){ double ph=2*pi*j/lon; V.push_back(point(st*cos(ph),st*sin(ph),ct)); }
    }
    for(int j=0;j<lon;j++) addTriConvex(V,F,0,id(1,j),id(1,j+1),fit.center);
    for(int r=1;r<=lat-2;r++) for(int j=0;j<lon;j++){
        int a=id(r,j), b=id(r,j+1), c=id(r+1,j+1), d=id(r+1,j);
        addTriConvex(V,F,a,b,c,fit.center); addTriConvex(V,F,a,c,d,fit.center);
    }
    for(int j=0;j<lon;j++) addTriConvex(V,F,1,id(lat-1,j+1),id(lat-1,j),fit.center);
    return true;
}
static bool tryEllipsoid(vector<Vec3>&outV, vector<Face>&outF, const Vec3 pcaAxes[3]){
    if(N<300 || elapsed()>12.0) return false;
    EllFit fit=fitEllipsoid(pcaAxes); if(!fit.ok) return false;
    double maxR=max({fit.r[0],fit.r[1],fit.r[2]});
    int minLon=max(12,(int)ceil(2*acos(-1.0)*maxR/(epsHaus*1.25)));
    int minLat=max(6,(int)ceil(acos(-1.0)*maxR/(epsHaus*1.25)));
    vector<double> ratios;
    if(N<1000) ratios={0.35,0.50,0.70};
    else if(N<5000) ratios={0.18,0.24,0.32,0.45};
    else if(N<30000) ratios={0.075,0.10,0.14,0.20};
    else if(N<120000) ratios={0.055,0.075,0.10,0.14};
    else ratios={0.035,0.050,0.070,0.095};
    for(double rr:ratios){
        if(elapsed()>16.5) break;
        int target=max(20,(int)(rr*N));
        int lon=max(minLon,(int)ceil(sqrt(2.0*target)));
        int lat=max(minLat,lon/2);
        if(lon%2) lon++;
        vector<Vec3> V; vector<Face> F; if(!buildEllipsoid(fit,lat,lon,V,F)) continue;
        if((int)V.size()>=N) continue;
        if(!validateHausdorff(V,F)) continue;
        double proxy=bestVisualProxy(V,F,(N<30000?160:128));
        double th=(N<5000?0.925:(N<30000?0.920:0.915));
        if(proxy>=th){ outV.swap(V); outF.swap(F); return true; }
    }
    return false;
}

struct CylFit{ bool ok=false; Vec3 w,u,v; double t0,t1,cu,cv,r0,r1; double rms=1e9,badFrac=1; };
static void coordsAxis(const Vec3&p,const Vec3&w,const Vec3&u,const Vec3&v,double&t,double&a,double&b){ t=dotv(p,w); a=dotv(p,u); b=dotv(p,v); }
static CylFit fitCylinderAxis(const Vec3&w0){
    CylFit fit; fit.w=normalized(w0); orthonormalBasis(fit.w,fit.u,fit.v);
    double tmin=1e100,tmax=-1e100, umin=1e100,umax=-1e100,vmin=1e100,vmax=-1e100;
    for(const Vec3&p:origV){ double t,a,b; coordsAxis(p,fit.w,fit.u,fit.v,t,a,b); tmin=min(tmin,t); tmax=max(tmax,t); umin=min(umin,a); umax=max(umax,a); vmin=min(vmin,b); vmax=max(vmax,b); }
    fit.t0=tmin; fit.t1=tmax; fit.cu=0.5*(umin+umax); fit.cv=0.5*(vmin+vmax);
    if(tmax-tmin<epsHaus || umax-umin<epsHaus || vmax-vmin<epsHaus) return fit;
    int stride=max(1,N/250000); double maxr=0;
    vector<array<double,2>> sample; sample.reserve(N/stride+1);
    for(int i=0;i<N;i+=stride){ double t,a,b; coordsAxis(origV[i],fit.w,fit.u,fit.v,t,a,b); double r=hypot(a-fit.cu,b-fit.cv); sample.push_back({t,r}); maxr=max(maxr,r); }
    if(maxr<=1e-10) return fit;
    double st=0,sr=0,stt=0,str=0; int cnt=0;
    for(auto pr:sample){ if(pr[1] > 0.78*maxr){ double t=pr[0],r=pr[1]; st+=t; sr+=r; stt+=t*t; str+=t*r; cnt++; } }
    if(cnt<max(80,(int)sample.size()/8)) return fit;
    double den=cnt*stt-st*st; if(fabs(den)<1e-18) return fit;
    double a=(cnt*str-st*sr)/den, b=(sr-a*st)/cnt;
    fit.r0=a*tmin+b; fit.r1=a*tmax+b;
    double mr=max(fit.r0,fit.r1); if(fit.r0<0.12*mr || fit.r1<0.12*mr || mr<=1e-10) return fit;
    double tol=max(epsHaus*0.33, mr*0.025);
    double s2=0; int sideCnt=0,bad=0,total=0;
    for(auto pr:sample){
        double t=pr[0],r=pr[1]; double rp=a*t+b; bool side=fabs(r-rp)<=tol; bool cap=(fabs(t-tmin)<=tol||fabs(t-tmax)<=tol) && r<=rp+tol;
        if(side){ s2+=(r-rp)*(r-rp); sideCnt++; }
        if(!(side||cap)) bad++; total++;
    }
    fit.rms=sqrt(s2/max(1,sideCnt))/max(mr,1e-12); fit.badFrac=(double)bad/max(1,total);
    fit.ok = (fit.rms<0.040 && fit.badFrac<0.035);
    return fit;
}
static Vec3 cylPoint(const CylFit&fit,double t,double r,double th){ return fit.w*t + fit.u*(fit.cu+r*cos(th)) + fit.v*(fit.cv+r*sin(th)); }
static bool buildCylinder(const CylFit&fit,int H,int S,int Rcap,vector<Vec3>&V,vector<Face>&F){
    H=max(H,1); S=max(S,8); Rcap=max(Rcap,1); const double pi=acos(-1.0);
    auto radAt=[&](int i){ double a=(double)i/H; return fit.r0*(1-a)+fit.r1*a; };
    int sideStart=0; V.clear(); F.clear(); V.reserve((H+1)*S+2*Rcap*S+4); F.reserve(2*H*S+4*Rcap*S);
    for(int i=0;i<=H;i++){ double t=fit.t0+(fit.t1-fit.t0)*i/H, r=radAt(i); for(int j=0;j<S;j++) V.push_back(cylPoint(fit,t,r,2*pi*j/S)); }
    auto sid=[&](int i,int j){ return sideStart+i*S+((j%S+S)%S); };
    for(int i=0;i<H;i++) for(int j=0;j<S;j++){ int a=sid(i,j), b=sid(i,j+1), c=sid(i+1,j+1), d=sid(i+1,j); F.push_back(Face{{a,b,c}}); F.push_back(Face{{a,c,d}}); }
    auto addCap=[&](bool top){
        double t=top?fit.t1:fit.t0, rmax=top?fit.r1:fit.r0; int boundaryRow=top?H:0; if(rmax<=1e-12) return;
        vector<vector<int>> ring(Rcap+1);
        ring[Rcap].resize(S); for(int j=0;j<S;j++) ring[Rcap][j]=sid(boundaryRow,j);
        int center=(int)V.size(); V.push_back(cylPoint(fit,t,0,0)); ring[0]={center};
        for(int rr=1;rr<Rcap;rr++){ ring[rr].resize(S); double r=rmax*(double)rr/Rcap; for(int j=0;j<S;j++){ ring[rr][j]=(int)V.size(); V.push_back(cylPoint(fit,t,r,2*pi*j/S)); } }
        if(Rcap==1){
            for(int j=0;j<S;j++){
                if(top) F.push_back(Face{{center,ring[1][j],ring[1][(j+1)%S]}});
                else F.push_back(Face{{center,ring[1][(j+1)%S],ring[1][j]}});
            }
        }else{
            for(int j=0;j<S;j++){
                if(top) F.push_back(Face{{center,ring[1][j],ring[1][(j+1)%S]}});
                else F.push_back(Face{{center,ring[1][(j+1)%S],ring[1][j]}});
            }
            for(int rr=1;rr<Rcap;rr++) for(int j=0;j<S;j++){
                int a=ring[rr][j], b=ring[rr][(j+1)%S], c=ring[rr+1][(j+1)%S], d=ring[rr+1][j];
                if(top){ F.push_back(Face{{a,d,c}}); F.push_back(Face{{a,c,b}}); }
                else { F.push_back(Face{{a,b,c}}); F.push_back(Face{{a,c,d}}); }
            }
        }
    };
    addCap(false); addCap(true);
    return (int)V.size()<=N;
}
static bool tryCylinder(vector<Vec3>&outV,vector<Face>&outF,const Vec3 pcaAxes[3]){
    if(N<500||elapsed()>12.5) return false;
    vector<Vec3> axes={pcaAxes[0],pcaAxes[1],pcaAxes[2],Vec3(1,0,0),Vec3(0,1,0),Vec3(0,0,1)};
    CylFit best;
    for(const Vec3&a:axes){ CylFit f=fitCylinderAxis(a); if(f.ok && (!best.ok || f.badFrac<best.badFrac || (fabs(f.badFrac-best.badFrac)<1e-9 && f.rms<best.rms))) best=f; }
    if(!best.ok) return false;
    double height=best.t1-best.t0, maxr=max(best.r0,best.r1);
    int minH=max(1,(int)ceil(height/(epsHaus*1.30)));
    int minS=max(12,(int)ceil(2*acos(-1.0)*maxr/(epsHaus*1.25)));
    int minR=max(1,(int)ceil(maxr/(epsHaus*1.30)));
    vector<double> ratios = (N<5000?vector<double>{0.18,0.28,0.40}:N<30000?vector<double>{0.075,0.11,0.16}:vector<double>{0.045,0.065,0.090});
    for(double rr:ratios){
        if(elapsed()>16.5) break;
        int target=max(30,(int)(rr*N));
        int S=max(minS,(int)ceil(sqrt((double)target*2.0))); if(S%2) S++;
        int H=max(minH,(int)ceil((double)S*height/(2*acos(-1.0)*maxr+1e-12)*0.7));
        int R=max(minR,(int)ceil(sqrt((double)target/(2.0*S))));
        vector<Vec3> V; vector<Face> F; if(!buildCylinder(best,H,S,R,V,F)) continue;
        if((int)V.size()>=N) continue;
        if(!validateHausdorff(V,F)) continue;
        double proxy=bestVisualProxy(V,F,(N<30000?160:128));
        if(proxy>=(N<5000?0.925:0.915)){ outV.swap(V); outF.swap(F); return true; }
    }
    return false;
}

struct TorFit{ bool ok=false; Vec3 w,u,v; double ct,cu,cv,R,r; double rms=1e9,maxe=1e9; };
static TorFit fitTorusAxis(const Vec3&w0){
    TorFit fit; fit.w=normalized(w0); orthonormalBasis(fit.w,fit.u,fit.v);
    double tmin=1e100,tmax=-1e100, umin=1e100,umax=-1e100,vmin=1e100,vmax=-1e100;
    for(const Vec3&p:origV){ double t,a,b; coordsAxis(p,fit.w,fit.u,fit.v,t,a,b); tmin=min(tmin,t); tmax=max(tmax,t); umin=min(umin,a); umax=max(umax,a); vmin=min(vmin,b); vmax=max(vmax,b); }
    fit.ct=0.5*(tmin+tmax); fit.cu=0.5*(umin+umax); fit.cv=0.5*(vmin+vmax);
    int stride=max(1,N/250000); double rhomin=1e100,rhomax=0;
    for(int i=0;i<N;i+=stride){ double t,a,b; coordsAxis(origV[i],fit.w,fit.u,fit.v,t,a,b); double rho=hypot(a-fit.cu,b-fit.cv); rhomin=min(rhomin,rho); rhomax=max(rhomax,rho); }
    if(!(rhomax>rhomin) || rhomin<epsHaus*0.5) return fit;
    fit.R=0.5*(rhomax+rhomin); double rr=0.5*(rhomax-rhomin), rt=0.5*(tmax-tmin); fit.r=0.5*(rr+rt);
    if(fit.r<=1e-10 || fit.R<1.4*fit.r) return fit;
    double s2=0,me=0; int cnt=0;
    for(int i=0;i<N;i+=stride){ double t,a,b; coordsAxis(origV[i],fit.w,fit.u,fit.v,t,a,b); double rho=hypot(a-fit.cu,b-fit.cv); double e=fabs(hypot(rho-fit.R,t-fit.ct)-fit.r); s2+=e*e; me=max(me,e); cnt++; }
    fit.rms=sqrt(s2/max(1,cnt))/fit.r; fit.maxe=me/fit.r;
    fit.ok=(fit.rms<0.045 && fit.maxe<0.22);
    return fit;
}
static Vec3 torPoint(const TorFit&fit,double th,double ph){
    double rho=fit.R+fit.r*cos(ph);
    return fit.w*(fit.ct+fit.r*sin(ph)) + fit.u*(fit.cu+rho*cos(th)) + fit.v*(fit.cv+rho*sin(th));
}
static bool buildTorus(const TorFit&fit,int S,int T,vector<Vec3>&V,vector<Face>&F){
    if(S<8||T<6) return false; const double pi=acos(-1.0); int vc=S*T; if(vc>N) return false;
    V.clear(); F.clear(); V.reserve(vc); F.reserve(2*vc);
    for(int i=0;i<S;i++) for(int j=0;j<T;j++) V.push_back(torPoint(fit,2*pi*i/S,2*pi*j/T));
    auto id=[&](int i,int j){ return ((i%S+S)%S)*T+((j%T+T)%T); };
    for(int i=0;i<S;i++) for(int j=0;j<T;j++){ int a=id(i,j), b=id(i+1,j), c=id(i+1,j+1), d=id(i,j+1); F.push_back(Face{{a,b,c}}); F.push_back(Face{{a,c,d}}); }
    return true;
}
static bool tryTorus(vector<Vec3>&outV,vector<Face>&outF,const Vec3 pcaAxes[3]){
    if(N<800||elapsed()>12.5) return false;
    vector<Vec3> axes={pcaAxes[0],pcaAxes[1],pcaAxes[2],Vec3(1,0,0),Vec3(0,1,0),Vec3(0,0,1)};
    TorFit best;
    for(const Vec3&a:axes){ TorFit f=fitTorusAxis(a); if(f.ok && (!best.ok || f.rms<best.rms)) best=f; }
    if(!best.ok) return false;
    int minS=max(12,(int)ceil(2*acos(-1.0)*(best.R+best.r)/(epsHaus*1.25)));
    int minT=max(8,(int)ceil(2*acos(-1.0)*best.r/(epsHaus*1.25)));
    vector<double> ratios=(N<5000?vector<double>{0.20,0.30,0.45}:N<30000?vector<double>{0.085,0.12,0.18}:vector<double>{0.045,0.065,0.090});
    for(double rr:ratios){
        if(elapsed()>16.5) break;
        int target=max(48,(int)(rr*N));
        double aspect=(best.R+best.r)/max(best.r,1e-12);
        int T=max(minT,(int)ceil(sqrt(target/aspect))); int S=max(minS,(int)ceil((double)target/T));
        if(S%2)S++; if(T%2)T++;
        vector<Vec3> V; vector<Face> F; if(!buildTorus(best,S,T,V,F)) continue;
        if((int)V.size()>=N) continue;
        if(!validateHausdorff(V,F)) continue;
        double proxy=bestVisualProxy(V,F,(N<30000?160:128));
        if(proxy>=(N<5000?0.925:0.915)){ outV.swap(V); outF.swap(F); return true; }
    }
    return false;
}

static bool tryStructural(vector<Vec3>&outV, vector<Face>&outF){
    if(trySphereGrid(outV,outF)) return true;
    if(tryTorusGrid(outV,outF)) return true;
    if(tryBox(outV,outF)) return true;
    Vec3 axes[3]; computePCA(axes);
    // More specific non-convex/ruled surfaces before the generic ellipsoid.
    if(tryTorus(outV,outF,axes)) return true;
    if(tryCylinder(outV,outF,axes)) return true;
    if(tryEllipsoid(outV,outF,axes)) return true;
    return false;
}

// ---------------- topology-preserving endpoint edge collapse fallback ----------------
struct Decimator{
    vector<Vec3> P;
    vector<Face> F;
    vector<vector<int>> inc;
    vector<unsigned char> aliveV,aliveF;
    vector<double> rad;
    vector<int> mark1,mark2;
    int iter1=1,iter2=1, aliveCnt, aliveFaces;
    double minArea2;
    struct Param{ double maxRad, planeTol, normalCos, creaseCos; };
    struct Cand{ float cost; int u,v; bool operator<(const Cand&o)const{return cost>o.cost;} };
    Decimator(){
        P=origV; F=origF; aliveCnt=N; aliveFaces=M; aliveV.assign(N,1); aliveF.assign(M,1); rad.assign(N,0);
        vector<int> deg(N,0); for(const Face&f:F){deg[f.v[0]]++;deg[f.v[1]]++;deg[f.v[2]]++;}
        inc.resize(N); for(int i=0;i<N;i++) inc[i].reserve(deg[i]+4);
        for(int i=0;i<M;i++){ inc[F[i].v[0]].push_back(i); inc[F[i].v[1]].push_back(i); inc[F[i].v[2]].push_back(i); }
        mark1.assign(N,0); mark2.assign(N,0); minArea2=max(1e-30,1e-24*bboxDiag*bboxDiag*bboxDiag*bboxDiag);
    }
    inline bool contains(int fid,int v)const{ const Face&f=F[fid]; return f.v[0]==v||f.v[1]==v||f.v[2]==v; }
    inline int oppVertex(int fid,int a,int b)const{ const Face&f=F[fid]; for(int k=0;k<3;k++){int x=f.v[k]; if(x!=a&&x!=b) return x;} return -1; }
    inline Vec3 faceCross(const Face&f)const{ return crossv(P[f.v[1]]-P[f.v[0]], P[f.v[2]]-P[f.v[0]]); }
    static uint64_t edgeKey(int a,int b){ if(a>b) swap(a,b); return (uint64_t)(uint32_t)a<<32 | (uint32_t)b; }
    bool findEdgeFaces(int u,int v,int ef[2],int opp[2]){
        ef[0]=ef[1]=opp[0]=opp[1]=-1; int cnt=0;
        const vector<int>&L = (inc[u].size()<inc[v].size()?inc[u]:inc[v]);
        for(int fid:L){ if(!aliveF[fid]) continue; if(contains(fid,u)&&contains(fid,v)){ if(cnt>=2) return false; ef[cnt]=fid; opp[cnt]=oppVertex(fid,u,v); cnt++; } }
        return cnt==2 && opp[0]>=0 && opp[1]>=0 && opp[0]!=opp[1];
    }
    bool linkOK(int u,int v,const int opp[2]){
        if(++iter1>2000000000){ fill(mark1.begin(),mark1.end(),0); iter1=1; }
        if(++iter2>2000000000){ fill(mark2.begin(),mark2.end(),0); iter2=1; }
        for(int fid:inc[u]) if(aliveF[fid]&&contains(fid,u)){
            const Face&f=F[fid]; for(int k=0;k<3;k++){ int x=f.v[k]; if(x!=u&&x!=v) mark1[x]=iter1; }
        }
        int common=0;
        for(int fid:inc[v]) if(aliveF[fid]&&contains(fid,v)){
            const Face&f=F[fid]; for(int k=0;k<3;k++){ int x=f.v[k]; if(x==u||x==v) continue; if(mark1[x]!=iter1) continue; if(x!=opp[0]&&x!=opp[1]) return false; if(mark2[x]!=iter2){mark2[x]=iter2; common++;} }
        }
        return common==2 && mark2[opp[0]]==iter2 && mark2[opp[1]]==iter2;
    }
    static array<int,3> sortedTri(int a,int b,int c){ array<int,3>s={a,b,c}; sort(s.begin(),s.end()); return s; }
    bool duplicateAfter(int keep,int rem,int skip0,int skip1,int fid,int a,int b,int c){
        int probe=a; if(inc[b].size()<inc[probe].size()) probe=b; if(inc[c].size()<inc[probe].size()) probe=c;
        auto target=sortedTri(a,b,c);
        for(int g:inc[probe]){
            if(!aliveF[g]||g==fid||g==skip0||g==skip1) continue;
            if(contains(g,rem)) continue;
            const Face&h=F[g]; if(sortedTri(h.v[0],h.v[1],h.v[2])==target) return true;
        }
        return false;
    }
    struct Eval{ bool ok=false; int keep=-1,rem=-1,ef[2]; double newRad=0,cost=1e100; };
    Eval evalEndpoint(int keep,int rem,const int ef[2],const int opp[2],const Param&p){
        Eval e; e.keep=keep; e.rem=rem; e.ef[0]=ef[0]; e.ef[1]=ef[1];
        double len=normv(P[keep]-P[rem]);
        e.newRad=max(rad[keep],rad[rem]+len);
        if(e.newRad>p.maxRad) return e;
        Vec3 n0=faceCross(F[ef[0]]), n1=faceCross(F[ef[1]]); double l0=normv(n0), l1=normv(n1);
        if(l0>0&&l1>0){ double cd=dotv(n0,n1)/(l0*l1); if(cd<p.creaseCos && len>epsHaus*0.10) return e; }
        double maxPlane=0,maxNloss=0; int changed=0;
        for(int fid:inc[rem]){
            if(!aliveF[fid]||!contains(fid,rem)) continue;
            if(fid==ef[0]||fid==ef[1]) continue;
            if(contains(fid,keep)) return e;
            Face old=F[fid]; Face nf=old;
            for(int k=0;k<3;k++) if(nf.v[k]==rem) nf.v[k]=keep;
            if(nf.v[0]==nf.v[1]||nf.v[0]==nf.v[2]||nf.v[1]==nf.v[2]) return e;
            Vec3 co=faceCross(old), cn=faceCross(nf); double lo=normv(co), ln=normv(cn);
            if(!(lo>0)&&!(ln>0)) return e; if(ln*ln<=minArea2) return e;
            double nd=dotv(co,cn)/(lo*ln); nd=clampd(nd,-1,1); if(nd<p.normalCos) return e;
            Vec3 no=co/lo; double pd=fabs(dotv(no,P[keep]-P[old.v[0]])); if(pd>p.planeTol) return e;
            if(duplicateAfter(keep,rem,ef[0],ef[1],fid,nf.v[0],nf.v[1],nf.v[2])) return e;
            maxPlane=max(maxPlane,pd); maxNloss=max(maxNloss,1.0-nd); changed++;
        }
        if(changed==0) return e;
        e.ok=true;
        e.cost = len/(epsHaus+1e-30) + 0.7*e.newRad/(p.maxRad+1e-30) + 40.0*maxPlane/(p.planeTol+1e-30) + 250.0*maxNloss + 0.0005*changed;
        return e;
    }
    bool collapseEdge(int u,int v,const Param&p){
        if(u==v||!aliveV[u]||!aliveV[v]) return false;
        int ef[2],opp[2]; if(!findEdgeFaces(u,v,ef,opp)) return false; if(!linkOK(u,v,opp)) return false;
        Eval a=evalEndpoint(u,v,ef,opp,p), b=evalEndpoint(v,u,ef,opp,p);
        Eval e; if(a.ok&&(!b.ok||a.cost<=b.cost)) e=a; else if(b.ok) e=b; else return false;
        int keep=e.keep, rem=e.rem;
        for(int k=0;k<2;k++) if(aliveF[e.ef[k]]){ aliveF[e.ef[k]]=0; aliveFaces--; }
        for(int fid:inc[rem]){
            if(!aliveF[fid]||!contains(fid,rem)) continue;
            Face&f=F[fid]; for(int k=0;k<3;k++) if(f.v[k]==rem) f.v[k]=keep;
        }
        aliveV[rem]=0; aliveCnt--; rad[keep]=e.newRad;
        vector<int> merged; merged.reserve(inc[keep].size()+inc[rem].size()+8);
        for(int fid:inc[keep]) if(aliveF[fid]&&contains(fid,keep)) merged.push_back(fid);
        for(int fid:inc[rem]) if(aliveF[fid]&&contains(fid,keep)) merged.push_back(fid);
        sort(merged.begin(),merged.end()); merged.erase(unique(merged.begin(),merged.end()),merged.end()); inc[keep].swap(merged); vector<int>().swap(inc[rem]);
        return true;
    }
    priority_queue<Cand> buildPQ(){
        vector<uint64_t> keys; keys.reserve((size_t)aliveFaces*3);
        for(int i=0;i<(int)F.size();i++) if(aliveF[i]){
            Face&f=F[i]; keys.push_back(edgeKey(f.v[0],f.v[1])); keys.push_back(edgeKey(f.v[1],f.v[2])); keys.push_back(edgeKey(f.v[2],f.v[0]));
        }
        sort(keys.begin(),keys.end()); keys.erase(unique(keys.begin(),keys.end()),keys.end());
        vector<Cand> cs; cs.reserve(keys.size());
        for(uint64_t k:keys){ int u=(int)(k>>32), v=(int)(uint32_t)k; if(aliveV[u]&&aliveV[v]){ float c=(float)norm2(P[u]-P[v]); cs.push_back(Cand{c,u,v}); } }
        priority_queue<Cand> pq(less<Cand>(), move(cs)); return pq;
    }

    static double pointTriD2(const Vec3&p,const Vec3&a,const Vec3&b,const Vec3&c){
        Vec3 ab=b-a, ac=c-a, ap=p-a;
        double d1=dotv(ab,ap), d2=dotv(ac,ap);
        if(d1<=0 && d2<=0) return norm2(ap);
        Vec3 bp=p-b; double d3=dotv(ab,bp), d4=dotv(ac,bp);
        if(d3>=0 && d4<=d3) return norm2(bp);
        double vc=d1*d4-d3*d2;
        if(vc<=0 && d1>=0 && d3<=0){ double v=d1/(d1-d3); Vec3 q=a+ab*v; return norm2(p-q); }
        Vec3 cp=p-c; double d5=dotv(ab,cp), d6=dotv(ac,cp);
        if(d6>=0 && d5<=d6) return norm2(cp);
        double vb=d5*d2-d1*d6;
        if(vb<=0 && d2>=0 && d6<=0){ double w=d2/(d2-d6); Vec3 q=a+ac*w; return norm2(p-q); }
        double va=d3*d6-d5*d4;
        if(va<=0 && (d4-d3)>=0 && (d5-d6)>=0){ double w=(d4-d3)/((d4-d3)+(d5-d6)); Vec3 q=b+(c-b)*w; return norm2(p-q); }
        Vec3 n=crossv(ab,ac); double nn=norm2(n); if(nn<=0) return min({norm2(ap),norm2(bp),norm2(cp)});
        double dist=dotv(p-a,n); return (dist*dist)/nn;
    }
    struct Patch{ int v=-1, keep=-1; double cost=1e100; vector<int> incF, ring; vector<Face> nf; };
    bool skippedFace(int fid,const vector<int>&skip) const { for(int x:skip) if(x==fid) return true; return false; }
    int edgeCntAfterSkip(int a,int b,const vector<int>&skip) const{
        const vector<int>&L=(inc[a].size()<inc[b].size()?inc[a]:inc[b]); int cnt=0;
        for(int fid:L) if(aliveF[fid]&&!skippedFace(fid,skip)&&contains(fid,a)&&contains(fid,b)) cnt++;
        return cnt;
    }
    bool faceExistsAfterSkip(const Face&g,const vector<int>&skip) const{
        int a=g.v[0],b=g.v[1],c=g.v[2]; int probe=a;
        if(inc[b].size()<inc[probe].size()) probe=b; if(inc[c].size()<inc[probe].size()) probe=c;
        auto t=sortedTri(a,b,c);
        for(int fid:inc[probe]) if(aliveF[fid]&&!skippedFace(fid,skip)){
            const Face&h=F[fid]; if(sortedTri(h.v[0],h.v[1],h.v[2])==t) return true;
        }
        return false;
    }
    bool buildPatchRing(int v, vector<int>&incF, vector<int>&ring) const{
        incF.clear(); ring.clear(); if(!aliveV[v]) return false;
        vector<pair<int,int>> bedges; vector<int> neigh;
        for(int fid:inc[v]) if(aliveF[fid]&&contains(fid,v)){
            const Face&f=F[fid]; int a=-1,b=-1;
            for(int k=0;k<3;k++){ int x=f.v[k]; if(x==v) continue; if(a<0) a=x; else b=x; }
            if(a<0||b<0||a==b||!aliveV[a]||!aliveV[b]) return false;
            incF.push_back(fid); bedges.push_back({a,b}); neigh.push_back(a); neigh.push_back(b);
        }
        sort(neigh.begin(),neigh.end()); neigh.erase(unique(neigh.begin(),neigh.end()),neigh.end());
        int k=(int)neigh.size(); if(k<3||k>10||(int)incF.size()!=k) return false;
        vector<vector<int>> adj(k);
        auto idx=[&](int x){ return (int)(lower_bound(neigh.begin(),neigh.end(),x)-neigh.begin()); };
        for(auto e:bedges){ int ia=idx(e.first), ib=idx(e.second); if(ia<0||ia>=k||ib<0||ib>=k||ia==ib) return false; adj[ia].push_back(ib); adj[ib].push_back(ia); }
        for(int i=0;i<k;i++){ sort(adj[i].begin(),adj[i].end()); adj[i].erase(unique(adj[i].begin(),adj[i].end()),adj[i].end()); if(adj[i].size()!=2) return false; }
        int start=0, prev=-1, cur=start;
        for(int step=0;step<k;step++){
            ring.push_back(neigh[cur]); int nxt=(adj[cur][0]==prev?adj[cur][1]:adj[cur][0]); prev=cur; cur=nxt;
        }
        return cur==start && (int)ring.size()==k;
    }
    bool inspectPatch(int v, Patch&pc){
        pc=Patch(); pc.v=v; if(!aliveV[v]) return false;
        vector<int> incF, ring; if(!buildPatchRing(v,incF,ring)) return false; int k=(int)ring.size();
        double bestD=1e100; int keep=-1; for(int u:ring){ double d=normv(P[v]-P[u]); if(d<bestD){bestD=d; keep=u;} }
        if(keep<0 || rad[v]+bestD>epsHaus*0.985) return false;
        Vec3 avg(0,0,0); for(int fid:incF){ Vec3 cr=faceCross(F[fid]); double l=normv(cr); if(l>0) avg=avg+cr; }
        double al=normv(avg); if(al<=0) return false; Vec3 an=avg/al;
        double cosLim=cos((N>120000?14.0:12.0)*acos(-1.0)/180.0);
        double planeLim=(N>47500?0.016:0.012)*bboxDiag;
        Patch best; best.v=v; best.keep=keep; best.cost=1e100;
        for(int r=0;r<k;r++){
            vector<int> cyc; cyc.reserve(k); for(int i=0;i<k;i++) cyc.push_back(ring[(r+i)%k]);
            vector<Face> nf; nf.reserve(k-2); bool ok=true; double maxLoss=0, minD2=1e100;
            for(int i=1;i+1<k;i++){
                Face g{{cyc[0],cyc[i],cyc[i+1]}};
                Vec3 cr=crossv(P[g.v[1]]-P[g.v[0]],P[g.v[2]]-P[g.v[0]]); double l=normv(cr); if(l*l<=minArea2){ ok=false; break; }
                if(dotv(cr,an)<0){ swap(g.v[1],g.v[2]); cr=cr*(-1.0); }
                double nd=dotv(cr/l,an); if(nd<cosLim){ ok=false; break; }
                if(faceExistsAfterSkip(g,incF)){ ok=false; break; }
                minD2=min(minD2, pointTriD2(P[v],P[g.v[0]],P[g.v[1]],P[g.v[2]]));
                maxLoss=max(maxLoss,1.0-nd); nf.push_back(g);
            }
            if(!ok||nf.empty()||sqrt(minD2)>planeLim) continue;
            vector<pair<int,int>> ee;
            for(const Face&g:nf) for(int q=0;q<3;q++){ int a=g.v[q], b=g.v[(q+1)%3]; if(a>b) swap(a,b); ee.push_back({a,b}); }
            sort(ee.begin(),ee.end());
            for(int i=0;i<(int)ee.size();){ int j=i; while(j<(int)ee.size()&&ee[j]==ee[i]) j++; int newc=j-i; int oldc=edgeCntAfterSkip(ee[i].first,ee[i].second,incF); if(oldc+newc!=2){ ok=false; break; } i=j; }
            if(!ok) continue;
            double cost=sqrt(minD2)/(bboxDiag+1e-30)+0.45*bestD/(epsHaus+1e-30)+40.0*maxLoss+0.01*k+0.1*rad[v]/(epsHaus+1e-30);
            if(cost<best.cost){ best.v=v; best.keep=keep; best.cost=cost; best.incF=incF; best.ring=ring; best.nf=nf; }
        }
        if(best.nf.empty()) return false; pc=move(best); return true;
    }
    bool applyPatchVertex(int v){
        Patch pc; if(!inspectPatch(v,pc)) return false; if(!aliveV[v]||!aliveV[pc.keep]) return false;
        for(int fid:pc.incF){ if(aliveF[fid]){ aliveF[fid]=0; aliveFaces--; } }
        aliveV[v]=0; aliveCnt--; rad[pc.keep]=max(rad[pc.keep], rad[v]+normv(P[v]-P[pc.keep])); vector<int>().swap(inc[v]);
        for(const Face&g:pc.nf){ int id=(int)F.size(); F.push_back(g); aliveF.push_back(1); aliveFaces++; inc[g.v[0]].push_back(id); inc[g.v[1]].push_back(id); inc[g.v[2]].push_back(id); }
        return true;
    }
    void patchPass(){
        if(elapsed()>17.6 || aliveCnt<50) return;
        vector<Patch> cand; cand.reserve(min(aliveCnt, N>200000?45000:60000));
        for(int v=0; v<N; ++v){
            if((v&1023)==0 && elapsed()>18.15) break;
            if(!aliveV[v]) continue; Patch pc; if(inspectPatch(v,pc)) cand.push_back(move(pc));
        }
        if(cand.empty()) return;
        sort(cand.begin(),cand.end(),[](const Patch&a,const Patch&b){ if(a.cost!=b.cost) return a.cost<b.cost; return a.v<b.v; });
        int cap=max(10,min(N>200000?22000:(N>60000?7000:3500), aliveCnt*(N>200000?10:(N>60000?8:6))/100));
        int done=0;
        for(const Patch&pc:cand){ if(done>=cap || elapsed()>18.95) break; if(applyPatchVertex(pc.v)) done++; }
    }
    void run(){
        double d=bboxDiag; bool huge=N>=200000;
        vector<Param> passes;
        auto add=[&](double r,double pl,double ang,double crease){ passes.push_back(Param{r*d,pl*d,cos(ang*acos(-1.0)/180.0),cos(crease*acos(-1.0)/180.0)}); };
        add(0.006,0.0010,1.2,18);
        add(0.012,0.0025,2.6,24);
        add(0.022,0.0055,4.8,32);
        if(huge) add(0.034,0.0100,7.8,42); else add(0.030,0.0085,6.6,38);
        if(huge) add(0.044,0.0140,10.5,50); else add(0.038,0.0110,8.8,45);
        int target;
        if(N<1000) target=max(8,(int)(0.18*N));
        else if(N<5000) target=max(20,(int)(0.11*N));
        else if(N<30000) target=max(50,(int)(0.085*N));
        else if(N<120000) target=max(200,(int)(0.075*N));
        else target=max(1000,(int)(0.060*N));
        for(size_t pi=0;pi<passes.size();pi++){
            if(elapsed()>18.2 || aliveCnt<=target) break;
            priority_queue<Cand> pq=buildPQ();
            int tries=0, applied=0;
            while(!pq.empty() && aliveCnt>target){
                if((++tries&4095)==0 && elapsed()>19.2) break;
                Cand c=pq.top(); pq.pop();
                if(!aliveV[c.u]||!aliveV[c.v]) continue;
                if(collapseEdge(c.u,c.v,passes[pi])) applied++;
            }
            if(applied==0 && pi+1>=passes.size()) break;
        }
        patchPass();
    }
    void output(vector<Vec3>&outV, vector<Face>&outF){
        vector<int> mapv(P.size(),-1); outF.clear(); outV.clear(); outF.reserve(aliveFaces);
        for(int i=0;i<(int)F.size();i++) if(aliveF[i]){
            Face f=F[i]; if(f.v[0]==f.v[1]||f.v[0]==f.v[2]||f.v[1]==f.v[2]) continue;
            if(norm2(crossv(P[f.v[1]]-P[f.v[0]],P[f.v[2]]-P[f.v[0]]))<=minArea2) continue;
            Face g;
            for(int k=0;k<3;k++){ int old=f.v[k]; if(mapv[old]<0){ mapv[old]=(int)outV.size(); outV.push_back(P[old]); } g.v[k]=mapv[old]; }
            outF.push_back(g);
        }
        if(outV.empty()||outV.size()>(size_t)N){ outV=origV; outF=origF; }
    }
};

int main(){
    tStart=chrono::steady_clock::now();
    readInput();
    vector<Vec3> outV; vector<Face> outF;
    if(tryStructural(outV,outF)){
        writeMesh(outV,outF);
        return 0;
    }
    Decimator dec; dec.run(); dec.output(outV,outF);
    writeMesh(outV,outF);
    return 0;
}
