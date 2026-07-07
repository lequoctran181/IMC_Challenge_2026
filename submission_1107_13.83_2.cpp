#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cstdint>
#include <queue>
#include <vector>
using namespace std;

struct Vec3{ double x,y,z; };
static inline Vec3 operator+(const Vec3&a,const Vec3&b){ return {a.x+b.x,a.y+b.y,a.z+b.z}; }
static inline Vec3 operator-(const Vec3&a,const Vec3&b){ return {a.x-b.x,a.y-b.y,a.z-b.z}; }
static inline Vec3 operator*(const Vec3&a,double s){ return {a.x*s,a.y*s,a.z*s}; }
static inline double dot3(const Vec3&a,const Vec3&b){ return a.x*b.x+a.y*b.y+a.z*b.z; }
static inline Vec3 cross3(const Vec3&a,const Vec3&b){ return {a.y*b.z-a.z*b.y,a.z*b.x-a.x*b.z,a.x*b.y-a.y*b.x}; }
static inline double norm2(const Vec3&a){ return dot3(a,a); }
static inline double norm3(const Vec3&a){ return sqrt(norm2(a)); }

struct Face{ int v[3]; unsigned char on=1; };
struct FastInput{
    vector<char> buf; char *p=nullptr;
    FastInput(){ char tmp[1<<16]; size_t n; while((n=fread(tmp,1,sizeof(tmp),stdin))>0) buf.insert(buf.end(),tmp,tmp+n); buf.push_back('\0'); p=buf.data(); }
    inline void skip(){ while(*p==' '||*p=='\n'||*p=='\r'||*p=='\t') ++p; }
    int nextInt(){ skip(); int s=1; if(*p=='-'){s=-1;++p;} int x=0; while(*p>='0'&&*p<='9'){ x=x*10+(*p-'0'); ++p; } return x*s; }
    double nextDouble(){ skip(); char*e=nullptr; double x=strtod(p,&e); p=e; return x; }
    char nextChar(){ skip(); return *p++; }
};

static int N,M;
static vector<Vec3> P,Orig;
static vector<Face> F;
static vector<vector<int>> inc;
static vector<unsigned char> alive,prot;
static vector<double> rad;
static vector<int> markv,markf;
static int tokenv=11, tokenf=17;
static int aliveV, aliveF;
static Vec3 bb0,bb1;
static double diagLen=1.0, Hlim=1.0, H2=1.0, minArea2=1e-30;
static chrono::steady_clock::time_point T0;
static double TL=20.15;

static inline double elapsed(){ return chrono::duration<double>(chrono::steady_clock::now()-T0).count(); }
static inline bool hasv(const Face&f,int v){ return f.v[0]==v||f.v[1]==v||f.v[2]==v; }
static inline uint64_t keyEdge(int a,int b){ if(a>b) swap(a,b); return (uint64_t)(uint32_t)a<<32 | (uint32_t)b; }
static inline uint64_t keyTri(int a,int b,int c){ if(a>b) swap(a,b); if(b>c) swap(b,c); if(a>b) swap(a,b); return ((uint64_t)(uint32_t)a<<42) ^ ((uint64_t)(uint32_t)b<<21) ^ (uint32_t)c; }
static inline int otherInFace(const Face&f,int a,int b){ for(int k=0;k<3;k++){ int x=f.v[k]; if(x!=a&&x!=b) return x; } return -1; }
static inline Vec3 faceNormal(const Face&f){ return cross3(P[f.v[1]]-P[f.v[0]],P[f.v[2]]-P[f.v[0]]); }

static void readInput(){
    FastInput in; N=in.nextInt(); M=in.nextInt();
    P.resize(N); Orig.resize(N); alive.assign(N,1); prot.assign(N,0); rad.assign(N,0.0);
    bb0={1e100,1e100,1e100}; bb1={-1e100,-1e100,-1e100};
    for(int i=0;i<N;i++){
        (void)in.nextChar(); P[i].x=in.nextDouble(); P[i].y=in.nextDouble(); P[i].z=in.nextDouble(); Orig[i]=P[i];
        bb0.x=min(bb0.x,P[i].x); bb0.y=min(bb0.y,P[i].y); bb0.z=min(bb0.z,P[i].z);
        bb1.x=max(bb1.x,P[i].x); bb1.y=max(bb1.y,P[i].y); bb1.z=max(bb1.z,P[i].z);
    }
    Vec3 d=bb1-bb0; diagLen=max(1e-12,norm3(d)); Hlim=0.0492*diagLen; H2=Hlim*Hlim; minArea2=1e-28*diagLen*diagLen*diagLen*diagLen;
    double eps=diagLen*1e-12;
    for(int i=0;i<N;i++) if(fabs(P[i].x-bb0.x)<=eps||fabs(P[i].x-bb1.x)<=eps||fabs(P[i].y-bb0.y)<=eps||fabs(P[i].y-bb1.y)<=eps||fabs(P[i].z-bb0.z)<=eps||fabs(P[i].z-bb1.z)<=eps) prot[i]=1;
    F.resize(M); vector<int> deg(N,0);
    for(int i=0;i<M;i++){
        (void)in.nextChar(); int a=in.nextInt()-1,b=in.nextInt()-1,c=in.nextInt()-1;
        F[i].v[0]=a; F[i].v[1]=b; F[i].v[2]=c; ++deg[a]; ++deg[b]; ++deg[c];
    }
    inc.resize(N); for(int i=0;i<N;i++) inc[i].reserve(deg[i]+8);
    for(int i=0;i<M;i++){ inc[F[i].v[0]].push_back(i); inc[F[i].v[1]].push_back(i); inc[F[i].v[2]].push_back(i); }
    aliveV=N; aliveF=M; markv.assign(N,0); markf.assign(M,0);
}

static void compactInc(int v){
    if(v<0 || !alive[v]) return;
    vector<int>& a=inc[v];
    if(a.size()<96) return;
    int cnt=0; for(int id:a) if(id>=0&&id<(int)F.size()&&F[id].on&&hasv(F[id],v)) ++cnt;
    if(cnt*3+32>=(int)a.size()) return;
    vector<int> b; b.reserve(cnt+8); for(int id:a) if(id>=0&&id<(int)F.size()&&F[id].on&&hasv(F[id],v)) b.push_back(id); a.swap(b);
}

static bool edgeFaces(int a,int b,int ef[2],int op[2]){
    ef[0]=ef[1]=op[0]=op[1]=-1; int cnt=0;
    const vector<int>& s = inc[a].size()<inc[b].size()?inc[a]:inc[b];
    for(int id:s){
        if(!F[id].on) continue; const Face&f=F[id];
        if(!hasv(f,a)||!hasv(f,b)) continue;
        if(cnt>=2) return false;
        ef[cnt]=id; op[cnt]=otherInFace(f,a,b); if(op[cnt]<0) return false; ++cnt;
    }
    return cnt==2 && op[0]!=op[1];
}

static bool linkOK(int a,int b,const int op[2]){
    if(tokenv>2000000000){ fill(markv.begin(),markv.end(),0); tokenv=11; }
    int ta=tokenv++, tb=tokenv++;
    for(int id:inc[a]) if(F[id].on&&hasv(F[id],a)){
        const Face&f=F[id]; for(int k=0;k<3;k++){ int x=f.v[k]; if(x!=a&&x!=b) markv[x]=ta; }
    }
    int common=0, g0=0,g1=0;
    for(int id:inc[b]) if(F[id].on&&hasv(F[id],b)){
        const Face&f=F[id]; for(int k=0;k<3;k++){ int x=f.v[k]; if(x==a||x==b) continue; if(markv[x]==ta){ markv[x]=tb; ++common; if(x==op[0]) g0=1; if(x==op[1]) g1=1; if(common>2) return false; } }
    }
    return common==2 && g0 && g1;
}

static bool duplicateFaceAfter(int keep,int rem,int fid,int a,int b,int c,int e0,int e1){
    int s=a; if(inc[b].size()<inc[s].size()) s=b; if(inc[c].size()<inc[s].size()) s=c;
    uint64_t kt=keyTri(a,b,c);
    for(int id:inc[s]){
        if(id==fid||id==e0||id==e1||!F[id].on) continue;
        const Face&g=F[id]; if(hasv(g,rem)) continue;
        if(keyTri(g.v[0],g.v[1],g.v[2])==kt) return true;
    }
    return false;
}

struct DirEval{ bool ok=false; double cost=1e100, nr=0; };
static DirEval evalDir(int keep,int rem,const int ef[2],double minCos,double planeLim){
    DirEval r; if(prot[rem]) return r;
    double nd=norm3(P[keep]-P[rem]); if(max(rad[keep],rad[rem]+nd)>Hlim) return r;
    double worstAng=0, worstPlane=0, areaLoss=0; int changed=0;
    for(int id:inc[rem]){
        if(!F[id].on||!hasv(F[id],rem)) continue;
        if(id==ef[0]||id==ef[1]) continue;
        Face f=F[id]; int nv[3]={f.v[0],f.v[1],f.v[2]};
        for(int k=0;k<3;k++) if(nv[k]==rem) nv[k]=keep;
        if(nv[0]==nv[1]||nv[0]==nv[2]||nv[1]==nv[2]) return r;
        Vec3 oldn=faceNormal(f);
        Vec3 newn=cross3(P[nv[1]]-P[nv[0]],P[nv[2]]-P[nv[0]]);
        double o2=norm2(oldn), n2=norm2(newn); if(o2<=minArea2||n2<=minArea2) return r;
        double cs=dot3(oldn,newn)/sqrt(o2*n2); if(cs<minCos) return r;
        Vec3 un=oldn*(1.0/sqrt(o2)); double pl=fabs(dot3(un,P[keep]-P[f.v[0]])); if(pl>planeLim) return r;
        if(duplicateFaceAfter(keep,rem,id,nv[0],nv[1],nv[2],ef[0],ef[1])) return r;
        worstAng=max(worstAng,1.0-cs); worstPlane=max(worstPlane,pl); areaLoss=max(areaLoss,max(0.0,1.0-sqrt(n2/o2))); ++changed;
    }
    if(changed==0) return r;
    r.ok=true; r.nr=max(rad[keep],rad[rem]+nd);
    r.cost = (nd/(Hlim+1e-300))*0.7 + (worstPlane/(planeLim+1e-300))*0.9 + 90.0*worstAng + 0.03*changed + 0.02*areaLoss;
    if(prot[keep]) r.cost-=0.05;
    return r;
}

static bool collapseDir(int keep,int rem,const int ef[2],double newRad){
    for(int z=0;z<2;z++) if(F[ef[z]].on){ F[ef[z]].on=0; --aliveF; }
    for(int id:inc[rem]){
        if(!F[id].on||!hasv(F[id],rem)) continue;
        Face&f=F[id]; for(int k=0;k<3;k++) if(f.v[k]==rem) f.v[k]=keep;
        inc[keep].push_back(id);
    }
    alive[rem]=0; --aliveV; rad[keep]=newRad;
    vector<int> merged; merged.reserve(inc[keep].size()+inc[rem].size());
    if(tokenf>2000000000){ fill(markf.begin(),markf.end(),0); tokenf=17; }
    int tk=tokenf++;
    for(int id:inc[keep]) if(F[id].on&&hasv(F[id],keep)&&markf[id]!=tk){ markf[id]=tk; merged.push_back(id); }
    for(int id:inc[rem]) if(F[id].on&&hasv(F[id],keep)&&markf[id]!=tk){ markf[id]=tk; merged.push_back(id); }
    inc[keep].swap(merged); vector<int>().swap(inc[rem]);
    compactInc(keep); return true;
}

static bool tryEdge(int a,int b,double minCos,double planeLim){
    if(a==b||!alive[a]||!alive[b]) return false;
    int ef[2],op[2]; if(!edgeFaces(a,b,ef,op)) return false;
    if(!linkOK(a,b,op)) return false;
    DirEval ab=evalDir(a,b,ef,minCos,planeLim), ba=evalDir(b,a,ef,minCos,planeLim);
    if(!ab.ok&&!ba.ok) return false;
    if(ba.ok&&(!ab.ok||ba.cost<ab.cost)) return collapseDir(b,a,ef,ba.nr);
    return collapseDir(a,b,ef,ab.nr);
}

static double estimateSmoothness(){
    int stride=max(1,M/90000), samples=0, smooth=0, sharp=0, bad=0;
    double c20=cos(20.0*acos(-1.0)/180.0), c45=cos(45.0*acos(-1.0)/180.0);
    for(int i=0;i<M && samples<90000;i+=stride){
        if(!F[i].on) continue; const Face&f=F[i]; int e[3][2]={{f.v[0],f.v[1]},{f.v[1],f.v[2]},{f.v[2],f.v[0]}};
        for(int k=0;k<3;k++){
            int ef[2],op[2]; if(!edgeFaces(e[k][0],e[k][1],ef,op)){ ++bad; continue; }
            Vec3 n0=faceNormal(F[ef[0]]), n1=faceNormal(F[ef[1]]); double l0=norm3(n0),l1=norm3(n1); if(l0<=0||l1<=0){++bad;continue;}
            double c=dot3(n0,n1)/(l0*l1); if(c>c20) ++smooth; if(c<c45) ++sharp; ++samples;
        }
    }
    if(samples==0) return 0.0;
    double sr=(double)smooth/samples, hr=(double)sharp/samples, br=(double)bad/(samples+bad+1);
    return max(0.0,min(1.0,sr - 0.35*hr - 0.25*br));
}

static vector<uint64_t> collectEdges(){
    vector<uint64_t> e; e.reserve((size_t)aliveF*3);
    for(int i=0;i<M;i++) if(F[i].on){ const Face&f=F[i]; e.push_back(keyEdge(f.v[0],f.v[1])); e.push_back(keyEdge(f.v[1],f.v[2])); e.push_back(keyEdge(f.v[2],f.v[0])); }
    sort(e.begin(),e.end()); e.erase(unique(e.begin(),e.end()),e.end());
    sort(e.begin(),e.end(),[](uint64_t A,uint64_t B){
        int a=A>>32,b=(int)(A&0xffffffffu), c=B>>32,d=(int)(B&0xffffffffu);
        double da=norm2(P[a]-P[b]), db=norm2(P[c]-P[d]); return da<db;
    });
    return e;
}


static double originalOrientationSign(){
    Vec3 c=(bb0+bb1)*0.5; double sum=0.0; int stride=max(1,M/200000);
    for(int i=0;i<M;i+=stride){ const Face&f=F[i]; Vec3 a=Orig[f.v[0]],b=Orig[f.v[1]],d=Orig[f.v[2]]; Vec3 cr=cross3(b-a,d-a); Vec3 ce=(a+b+d)*(1.0/3.0); sum+=dot3(cr,ce-c); }
    return sum>=0.0?1.0:-1.0;
}

static void setCandidateMesh(const vector<Vec3>&V,const vector<Face>&G){
    vector<Vec3> oldP=P; vector<Face> oldF=F; vector<unsigned char> oldAlive=alive;
    P.assign(N,Vec3{0,0,0}); alive.assign(N,0); rad.assign(N,0.0);
    for(int i=0;i<(int)V.size();++i){ P[i]=V[i]; alive[i]=1; }
    F=G; M=(int)F.size(); aliveV=(int)V.size(); aliveF=M;
}

static bool tryEllipsoidPrimitive(){
    if(N<850 || elapsed()>2.0) return false;
    Vec3 c=(bb0+bb1)*0.5; double rx=(bb1.x-bb0.x)*0.5, ry=(bb1.y-bb0.y)*0.5, rz=(bb1.z-bb0.z)*0.5;
    double mx=max(rx,max(ry,rz)), mn=min(rx,min(ry,rz));
    if(!(mx>1e-12) || mn<mx*0.22) return false;
    int stride=max(1,N/220000), cnt=0; double ss=0.0, ma=0.0;
    for(int i=0;i<N;i+=stride){ double x=(Orig[i].x-c.x)/rx, y=(Orig[i].y-c.y)/ry, z=(Orig[i].z-c.z)/rz; double e=fabs(sqrt(max(0.0,x*x+y*y+z*z))-1.0); ss+=e*e; ma=max(ma,e); ++cnt; }
    if(cnt<200) return false;
    double rms=sqrt(ss/cnt);
    double rmsLim=(N<5000?0.0065:0.0048), maxLim=(N<5000?0.030:0.022);
    if(rms>rmsLim || ma>maxLim) return false;
    int lat,lon;
    if(N<1500){lat=16;lon=32;} else if(N<5000){lat=18;lon=36;} else if(N<50000){lat=22;lon=44;} else {lat=26;lon=52;}
    int VC=2+(lat-1)*lon; if(VC>=N || VC>N*92/100) return false;
    vector<Vec3> V; vector<Face> G; V.reserve(VC); G.reserve(2*lat*lon);
    const double pi=acos(-1.0);
    V.push_back({c.x,c.y,c.z+rz});
    for(int i=1;i<lat;i++){ double ph=pi*i/lat, sp=sin(ph), cp=cos(ph); for(int j=0;j<lon;j++){ double th=2*pi*j/lon; V.push_back({c.x+rx*sp*cos(th),c.y+ry*sp*sin(th),c.z+rz*cp}); } }
    V.push_back({c.x,c.y,c.z-rz});
    auto id=[&](int r,int j){ return 1+(r-1)*lon+((j%lon+lon)%lon); };
    auto add=[&](int a,int b,int cc){ Face f; f.v[0]=a; f.v[1]=b; f.v[2]=cc; f.on=1; G.push_back(f); };
    int bottom=VC-1;
    for(int j=0;j<lon;j++) add(0,id(1,j),id(1,j+1));
    for(int r=1;r<lat-1;r++) for(int j=0;j<lon;j++){ int a=id(r,j), b=id(r,j+1), cc=id(r+1,j), d=id(r+1,j+1); add(a,cc,b); add(b,cc,d); }
    for(int j=0;j<lon;j++) add(bottom,id(lat-1,j+1),id(lat-1,j));
    double want=originalOrientationSign(), got=0.0;
    for(const Face&f:G){ Vec3 a=V[f.v[0]],b=V[f.v[1]],d=V[f.v[2]]; got+=dot3(cross3(b-a,d-a),((a+b+d)*(1.0/3.0))-c); }
    if((got>=0?1.0:-1.0)!=want) for(Face&f:G) swap(f.v[1],f.v[2]);
    setCandidateMesh(V,G);
    return true;
}

static void simplify(){
    double smooth=estimateSmoothness();
    double targetRatio;
    double minCos, planeLim;
    if(N<600){ targetRatio=0.55; minCos=0.78; planeLim=0.010*diagLen; }
    else if(smooth>0.965){ targetRatio=0.070; minCos=0.22; planeLim=0.038*diagLen; }
    else if(smooth>0.88){ targetRatio=0.095; minCos=0.32; planeLim=0.032*diagLen; }
    else if(smooth>0.68){ targetRatio=0.125; minCos=0.44; planeLim=0.026*diagLen; }
    else { targetRatio=0.170; minCos=0.60; planeLim=0.018*diagLen; }
    if(N<2500) targetRatio=max(targetRatio,0.20); else if(N<9000) targetRatio=max(targetRatio,0.14);
    int target=max(12,(int)ceil(N*targetRatio));
    for(int pass=0;pass<6 && aliveV>target && elapsed()<TL;pass++){
        vector<uint64_t> edges=collectEdges();
        int changed=0; double passCos=min(0.84,minCos+0.055*pass); double passPlane=planeLim*(1.0-0.04*min(pass,3));
        for(size_t i=0;i<edges.size() && aliveV>target;i++){
            if((i&4095)==0 && elapsed()>TL) break;
            int a=(int)(edges[i]>>32), b=(int)(edges[i]&0xffffffffu);
            if(tryEdge(a,b,passCos,passPlane)) ++changed;
        }
        if(changed==0) break;
    }
}

static bool validOutput(vector<int>& remap, vector<int>& outFaces){
    remap.assign(N,-1); outFaces.clear(); outFaces.reserve(aliveF);
    vector<unsigned char> used(N,0);
    for(int i=0;i<M;i++) if(F[i].on){
        int a=F[i].v[0],b=F[i].v[1],c=F[i].v[2];
        if(a<0||b<0||c<0||a>=N||b>=N||c>=N||!alive[a]||!alive[b]||!alive[c]) return false;
        if(a==b||a==c||b==c) return false;
        if(norm2(faceNormal(F[i]))<=minArea2) return false;
        used[a]=used[b]=used[c]=1; outFaces.push_back(i);
    }
    if(outFaces.empty()) return false;
    int V=0; Vec3 mn{1e100,1e100,1e100}, mx{-1e100,-1e100,-1e100};
    for(int i=0;i<N;i++) if(used[i]){ remap[i]=V++; mn.x=min(mn.x,P[i].x); mn.y=min(mn.y,P[i].y); mn.z=min(mn.z,P[i].z); mx.x=max(mx.x,P[i].x); mx.y=max(mx.y,P[i].y); mx.z=max(mx.z,P[i].z); }
    double eps=diagLen*1e-10;
    if(fabs(mn.x-bb0.x)>eps||fabs(mn.y-bb0.y)>eps||fabs(mn.z-bb0.z)>eps||fabs(mx.x-bb1.x)>eps||fabs(mx.y-bb1.y)>eps||fabs(mx.z-bb1.z)>eps) return false;
    vector<uint64_t> tris; tris.reserve(outFaces.size());
    for(int id:outFaces) tris.push_back(keyTri(remap[F[id].v[0]],remap[F[id].v[1]],remap[F[id].v[2]]));
    sort(tris.begin(),tris.end()); for(size_t i=1;i<tris.size();i++) if(tris[i]==tris[i-1]) return false;
    return V>0 && V<=N && (int)outFaces.size()>0;
}

static void outputOriginal(){
    static char obuf[1<<20]; setvbuf(stdout,obuf,_IOFBF,sizeof(obuf));
    printf("%d %d\n",N,M);
    for(int i=0;i<N;i++) printf("v %.17g %.17g %.17g\n",Orig[i].x,Orig[i].y,Orig[i].z);
    for(int i=0;i<M;i++) printf("f %d %d %d\n",F[i].v[0]+1,F[i].v[1]+1,F[i].v[2]+1);
}

static void outputSimplified(){
    vector<int> remap, ids; if(!validOutput(remap,ids)){ outputOriginal(); return; }
    int V=0; for(int x:remap) if(x>=0) ++V;
    if(V>=N){ outputOriginal(); return; }
    static char obuf[1<<20]; setvbuf(stdout,obuf,_IOFBF,sizeof(obuf));
    printf("%d %d\n",V,(int)ids.size());
    vector<int> old(V,-1); for(int i=0;i<N;i++) if(remap[i]>=0) old[remap[i]]=i;
    for(int j=0;j<V;j++){ int i=old[j]; printf("v %.17g %.17g %.17g\n",P[i].x,P[i].y,P[i].z); }
    for(int id:ids) printf("f %d %d %d\n",remap[F[id].v[0]]+1,remap[F[id].v[1]]+1,remap[F[id].v[2]]+1);
}

int main(){
    T0=chrono::steady_clock::now();
    readInput();
    if(!tryEllipsoidPrimitive() && N>20 && M>20) simplify();
    outputSimplified();
    return 0;
}
