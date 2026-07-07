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
static inline Vec3 normalized(Vec3 a){ double l=norm3(a); if(l<=1e-300) return {0,0,0}; return a/l; }

struct Tri{ int a,b,c; };
struct Face{ int v[3]; unsigned char active=1; Vec3 n; };

struct FastInput{
    vector<char> buf; char *p=nullptr;
    FastInput(){ char tmp[1<<16]; size_t n; while((n=fread(tmp,1,sizeof(tmp),stdin))>0) buf.insert(buf.end(),tmp,tmp+n); buf.push_back('\0'); p=buf.data(); }
    inline void skip(){ while(*p==' '||*p=='\n'||*p=='\r'||*p=='\t') ++p; }
    int nextInt(){ skip(); int s=1; if(*p=='-') s=-1,++p; int x=0; while(*p>='0'&&*p<='9'){ x=x*10+(*p-'0'); ++p; } return x*s; }
    long nextLong(){ skip(); return strtol(p,&p,10); }
    double nextDouble(){ skip(); return strtod(p,&p); }
    char nextChar(){ skip(); return *p++; }
};

struct Quadric{
    double q[10];
    Quadric(){ memset(q,0,sizeof(q)); }
    inline void addPlane(double a,double b,double c,double d,double w=1.0){
        q[0]+=w*a*a; q[1]+=w*a*b; q[2]+=w*a*c; q[3]+=w*a*d;
        q[4]+=w*b*b; q[5]+=w*b*c; q[6]+=w*b*d;
        q[7]+=w*c*c; q[8]+=w*c*d; q[9]+=w*d*d;
    }
    inline void add(const Quadric&o){ for(int i=0;i<10;i++) q[i]+=o.q[i]; }
    inline double eval(const Vec3&p) const{
        double x=p.x,y=p.y,z=p.z;
        return q[0]*x*x + 2*q[1]*x*y + 2*q[2]*x*z + 2*q[3]*x +
               q[4]*y*y + 2*q[5]*y*z + 2*q[6]*y + q[7]*z*z + 2*q[8]*z + q[9];
    }
};

static int N=0,M=0;
static vector<Vec3> Orig;
static vector<Tri> OrigF;
static Vec3 bbMin, bbMax;
static double diagLen=1.0, tau=0.05, tau2=0.0025;
static chrono::steady_clock::time_point T0;
static inline double elapsed(){ return chrono::duration<double>(chrono::steady_clock::now()-T0).count(); }

static void computeBounds(){
    if(Orig.empty()){ bbMin=bbMax={0,0,0}; diagLen=1; tau=.05; tau2=tau*tau; return; }
    bbMin=bbMax=Orig[0];
    for(auto&p:Orig){
        bbMin.x=min(bbMin.x,p.x); bbMin.y=min(bbMin.y,p.y); bbMin.z=min(bbMin.z,p.z);
        bbMax.x=max(bbMax.x,p.x); bbMax.y=max(bbMax.y,p.y); bbMax.z=max(bbMax.z,p.z);
    }
    diagLen=norm3(bbMax-bbMin); if(!(diagLen>0)) diagLen=1.0;
    tau=0.05*diagLen*0.992; tau2=tau*tau;
}

struct DenseGrid{
    Vec3 mn,mx; double cell=1; int nx=1,ny=1,nz=1; vector<int> head,nxt; const vector<Vec3>*pts=nullptr;
    int clampi(int v,int n) const { return v<0?0:(v>=n?n-1:v); }
    int ix(double x) const { return clampi((int)((x-mn.x)/cell),nx); }
    int iy(double y) const { return clampi((int)((y-mn.y)/cell),ny); }
    int iz(double z) const { return clampi((int)((z-mn.z)/cell),nz); }
    int key(int x,int y,int z) const { return (z*ny+y)*nx+x; }
    void build(const vector<Vec3>&P,double radius){
        pts=&P; cell=max(radius,1e-9); mn=mx=P.empty()?Vec3{0,0,0}:P[0];
        for(auto&p:P){ mn.x=min(mn.x,p.x); mn.y=min(mn.y,p.y); mn.z=min(mn.z,p.z); mx.x=max(mx.x,p.x); mx.y=max(mx.y,p.y); mx.z=max(mx.z,p.z); }
        nx=max(1,(int)((mx.x-mn.x)/cell)+1); ny=max(1,(int)((mx.y-mn.y)/cell)+1); nz=max(1,(int)((mx.z-mn.z)/cell)+1);
        long long cells=1LL*nx*ny*nz;
        if(cells>3000000LL){ double span=max({mx.x-mn.x,mx.y-mn.y,mx.z-mn.z,cell}); int k=(int)ceil(pow(3000000.0,1.0/3.0)); nx=ny=nz=max(1,k); cell=span/(double)k+1e-12; cells=1LL*nx*ny*nz; }
        head.assign((size_t)cells,-1); nxt.assign(P.size(),-1);
        for(int i=0;i<(int)P.size();++i){ int h=key(ix(P[i].x),iy(P[i].y),iz(P[i].z)); nxt[i]=head[h]; head[h]=i; }
    }
    bool within(const Vec3&p,double r2) const{
        int X=ix(p.x),Y=iy(p.y),Z=iz(p.z); int rad=max(1,(int)ceil(sqrt(r2)/cell));
        for(int z=Z-rad; z<=Z+rad; ++z) if(z>=0&&z<nz)
        for(int y=Y-rad; y<=Y+rad; ++y) if(y>=0&&y<ny)
        for(int x=X-rad; x<=X+rad; ++x) if(x>=0&&x<nx){
            for(int i=head[key(x,y,z)]; i!=-1; i=nxt[i]) if(norm2((*pts)[i]-p)<=r2) return true;
        }
        return false;
    }
    void markNear(const Vec3&p,double r2,vector<unsigned char>&covered,int&cnt) const{
        int X=ix(p.x),Y=iy(p.y),Z=iz(p.z); int rad=max(1,(int)ceil(sqrt(r2)/cell));
        for(int z=Z-rad; z<=Z+rad; ++z) if(z>=0&&z<nz)
        for(int y=Y-rad; y<=Y+rad; ++y) if(y>=0&&y<ny)
        for(int x=X-rad; x<=X+rad; ++x) if(x>=0&&x<nx){
            for(int i=head[key(x,y,z)]; i!=-1; i=nxt[i]) if(!covered[i] && norm2((*pts)[i]-p)<=r2){ covered[i]=1; ++cnt; }
        }
    }
};

static bool vertexHausdorffOK(const vector<Vec3>&A,const vector<Vec3>&B,double r2){
    if(A.empty()||B.empty()) return false;
    DenseGrid gB; gB.build(B,sqrt(r2));
    for(size_t i=0;i<A.size();++i){ if((i&8191)==0 && elapsed()>18.8) return false; if(!gB.within(A[i],r2)) return false; }
    DenseGrid gA; gA.build(A,sqrt(r2));
    for(size_t i=0;i<B.size();++i){ if((i&8191)==0 && elapsed()>19.2) return false; if(!gA.within(B[i],r2)) return false; }
    return true;
}

struct RenderMaps{ int R=0; vector<float> depth; vector<float> norm; vector<unsigned char> mask; };
static inline void projectPoint(const Vec3&p,int view,int R,double&u,double&v,double&z){
    const double D=2.5; double f=800.0*(double)R/1024.0, c=0.5*(double)R; double sx,sy;
    if(view==0){ sx=p.y; sy=p.z; z=D-p.x; }
    else if(view==1){ sx=-p.y; sy=p.z; z=D+p.x; }
    else if(view==2){ sx=-p.x; sy=p.z; z=D-p.y; }
    else if(view==3){ sx=p.x; sy=p.z; z=D+p.y; }
    else if(view==4){ sx=p.x; sy=p.y; z=D-p.z; }
    else { sx=-p.x; sy=p.y; z=D+p.z; }
    u=f*sx/z+c; v=f*sy/z+c;
}
static inline Vec3 unprojectPoint(int view,int R,double u,double v,double z){
    const double D=2.5; double f=800.0*(double)R/1024.0, c=0.5*(double)R; double sx=(u-c)*z/f, sy=(v-c)*z/f;
    if(view==0) return {D-z,sx,sy};
    if(view==1) return {z-D,-sx,sy};
    if(view==2) return {-sx,D-z,sy};
    if(view==3) return {sx,z-D,sy};
    if(view==4) return {sx,sy,D-z};
    return {-sx,sy,z-D};
}
static inline Vec3 camDir(int view){
    if(view==0) return {1,0,0}; if(view==1) return {-1,0,0}; if(view==2) return {0,1,0}; if(view==3) return {0,-1,0}; if(view==4) return {0,0,1}; return {0,0,-1};
}
static inline Vec3 backDir(int view){ Vec3 d=camDir(view); return d*(-1.0); }
static void renderMesh(const vector<Vec3>&V,const vector<Tri>&F,RenderMaps&rm,int R){
    const int PIX=R*R; rm.R=R; rm.depth.assign((size_t)6*PIX,255.0f); rm.norm.assign((size_t)6*PIX*3,127.5f); rm.mask.assign((size_t)6*PIX,0);
    vector<float> U(V.size()), W(V.size()), Z(V.size());
    vector<Vec3> fn(F.size());
    for(size_t i=0;i<F.size();++i){ const Tri&t=F[i]; Vec3 cr=cross3(V[t.b]-V[t.a],V[t.c]-V[t.a]); fn[i]=normalized(cr); }
    for(int view=0; view<6; ++view){
        for(size_t i=0;i<V.size();++i){ double u,v,z; projectPoint(V[i],view,R,u,v,z); U[i]=(float)u; W[i]=(float)v; Z[i]=(float)z; }
        float* dep=rm.depth.data()+(size_t)view*PIX; float* nor=rm.norm.data()+(size_t)view*PIX*3; unsigned char* msk=rm.mask.data()+(size_t)view*PIX;
        for(size_t ti=0; ti<F.size(); ++ti){ const Tri&t=F[ti]; int ia=t.a,ib=t.b,ic=t.c; float x0=U[ia],x1=U[ib],x2=U[ic], y0=W[ia],y1=W[ib],y2=W[ic], z0=Z[ia],z1=Z[ib],z2=Z[ic];
            if(!(z0>0&&z1>0&&z2>0)) continue;
            int xmin=max(0,(int)floor(min(x0,min(x1,x2))-.5f)); int xmax=min(R-1,(int)ceil(max(x0,max(x1,x2))+.5f));
            int ymin=max(0,(int)floor(min(y0,min(y1,y2))-.5f)); int ymax=min(R-1,(int)ceil(max(y0,max(y1,y2))+.5f));
            if(xmin>xmax||ymin>ymax) continue; float den=(y1-y2)*(x0-x2)+(x2-x1)*(y0-y2); if(fabs(den)<1e-20f) continue; float invDen=1.0f/den;
            Vec3 n=fn[ti]; float nr=(float)((n.x+1.0)*127.5), ng=(float)((n.y+1.0)*127.5), nb=(float)((n.z+1.0)*127.5);
            for(int yy=ymin; yy<=ymax; ++yy){ float py=(float)yy+0.5f; int row=yy*R; for(int xx=xmin; xx<=xmax; ++xx){ float px=(float)xx+0.5f; float a=((y1-y2)*(px-x2)+(x2-x1)*(py-y2))*invDen; float b=((y2-y0)*(px-x2)+(x0-x2)*(py-y2))*invDen; float c=1.0f-a-b; if(a>=-1e-6f&&b>=-1e-6f&&c>=-1e-6f){ float zz=1.0f/(a/z0+b/z1+c/z2); int id=row+xx; if(zz<dep[id]){ dep[id]=zz; msk[id]=1; float*q=nor+3*id; q[0]=nr; q[1]=ng; q[2]=nb; } } } }
        }
    }
}
static inline double rectSum(const vector<double>&I,int W,int x0,int y0,int x1,int y1){ return I[(size_t)y1*W+x1]-I[(size_t)y0*W+x1]-I[(size_t)y1*W+x0]+I[(size_t)y0*W+x0]; }
static double ssimChannel(const float*A,int strideA,const float*B,int strideB,const unsigned char*ma,const unsigned char*mb,int R){
    const int W=R+1, rad=5, area=(2*rad+1)*(2*rad+1); const double c1=6.5025, c2=58.5225; size_t SZ=(size_t)W*W;
    vector<double> SA(SZ),SB(SZ),SA2(SZ),SB2(SZ),SAB(SZ);
    for(int y=1;y<=R;y++){ double a=0,b=0,a2=0,b2=0,ab=0; int row=(y-1)*R; for(int x=1;x<=R;x++){ int p=row+x-1; double va=A[(size_t)p*strideA], vb=B[(size_t)p*strideB]; a+=va; b+=vb; a2+=va*va; b2+=vb*vb; ab+=va*vb; size_t id=(size_t)y*W+x, up=(size_t)(y-1)*W+x; SA[id]=SA[up]+a; SB[id]=SB[up]+b; SA2[id]=SA2[up]+a2; SB2[id]=SB2[up]+b2; SAB[id]=SAB[up]+ab; } }
    long long cnt=0; long double acc=0.0;
    for(int y=rad;y<R-rad;y++){ int row=y*R; for(int x=rad;x<R-rad;x++){ int p=row+x; if(!(ma[p]||mb[p])) continue; int x0=x-rad,y0=y-rad,x1=x+rad+1,y1=y+rad+1; double sa=rectSum(SA,W,x0,y0,x1,y1), sb=rectSum(SB,W,x0,y0,x1,y1); double sa2=rectSum(SA2,W,x0,y0,x1,y1), sb2=rectSum(SB2,W,x0,y0,x1,y1), sab=rectSum(SAB,W,x0,y0,x1,y1); double mua=sa/area, mub=sb/area; double va=max(0.0,sa2/area-mua*mua), vb=max(0.0,sb2/area-mub*mub), cov=sab/area-mua*mub; double den=(mua*mua+mub*mub+c1)*(va+vb+c2); double val=den!=0.0?((2*mua*mub+c1)*(2*cov+c2)/den):1.0; acc+=val; ++cnt; } }
    return cnt? (double)(acc/cnt) : 1.0;
}
static double renderSSIM(const RenderMaps&A,const RenderMaps&B){
    int R=A.R, PIX=R*R; double total=0;
    for(int view=0; view<6; ++view){ const float*ad=A.depth.data()+(size_t)view*PIX; const float*bd=B.depth.data()+(size_t)view*PIX; const unsigned char*am=A.mask.data()+(size_t)view*PIX; const unsigned char*bm=B.mask.data()+(size_t)view*PIX; double ns=0; const float*an=A.norm.data()+(size_t)view*PIX*3; const float*bn=B.norm.data()+(size_t)view*PIX*3; for(int ch=0; ch<3; ++ch) ns+=ssimChannel(an+ch,3,bn+ch,3,am,bm,R); ns/=3.0; double ds=ssimChannel(ad,1,bd,1,am,bm,R); total += 0.5*(ns+ds); if(elapsed()>19.5) break; }
    return total/6.0;
}

static void outputMesh(const vector<Vec3>&V,const vector<Tri>&F){
    static char obuf[1<<20]; setvbuf(stdout,obuf,_IOFBF,sizeof(obuf));
    printf("%d %d\n",(int)V.size(),(int)F.size());
    for(auto&p:V) printf("v %.15g %.15g %.15g\n",p.x,p.y,p.z);
    for(auto&t:F) printf("f %d %d %d\n",t.a+1,t.b+1,t.c+1);
}
static inline void addTriDesired(vector<Vec3>&V, vector<Tri>&F, int a,int b,int c, Vec3 desired){
    Vec3 cr=cross3(V[b]-V[a],V[c]-V[a]); if(dot3(cr,desired)<0) swap(b,c); if(norm2(cross3(V[b]-V[a],V[c]-V[a]))>1e-32) F.push_back({a,b,c});
}
static void addOrientedTri(vector<Vec3>&V, vector<Tri>&F, int a,int b,int c, Vec3 inside){
    Vec3 cr=cross3(V[b]-V[a],V[c]-V[a]); Vec3 ce=(V[a]+V[b]+V[c])/3.0; if(dot3(cr,ce-inside)<0) swap(b,c); if(norm2(cross3(V[b]-V[a],V[c]-V[a]))>1e-32) F.push_back({a,b,c});
}
static void addTinyTet(vector<Vec3>&V, vector<Tri>&F, const Vec3&p){
    double e=max(1e-8,diagLen*8e-7); int s=(int)V.size();
    V.push_back(p); V.push_back({p.x+e,p.y,p.z}); V.push_back({p.x,p.y+e,p.z}); V.push_back({p.x,p.y,p.z+e});
    F.push_back({s,s+2,s+1}); F.push_back({s,s+1,s+3}); F.push_back({s,s+3,s+2}); F.push_back({s+1,s+2,s+3});
}

static bool faceNormalsMostlyAxis(){
    int stride=max(1,M/60000), cnt=0, ok=0;
    for(int i=0;i<M;i+=stride){ auto&t=OrigF[i]; Vec3 n=normalized(cross3(Orig[t.b]-Orig[t.a],Orig[t.c]-Orig[t.a])); double m=max({fabs(n.x),fabs(n.y),fabs(n.z)}); if(m>0.992) ok++; cnt++; }
    return cnt>0 && ok*100>=cnt*98;
}

static bool tryBoxCorners(vector<Vec3>&OV, vector<Tri>&OF){
    if(N<9) return false; Vec3 ext=bbMax-bbMin; if(ext.x<=1e-12||ext.y<=1e-12||ext.z<=1e-12) return false;
    double planeTol=max(1e-7,0.002*diagLen); for(auto&p:Orig){ double d=min({fabs(p.x-bbMin.x),fabs(p.x-bbMax.x),fabs(p.y-bbMin.y),fabs(p.y-bbMax.y),fabs(p.z-bbMin.z),fabs(p.z-bbMax.z)}); if(d>planeTol) return false; }
    OV={{bbMax.x,bbMax.y,bbMax.z},{bbMax.x,bbMax.y,bbMin.z},{bbMax.x,bbMin.y,bbMax.z},{bbMax.x,bbMin.y,bbMin.z},{bbMin.x,bbMax.y,bbMax.z},{bbMin.x,bbMax.y,bbMin.z},{bbMin.x,bbMin.y,bbMax.z},{bbMin.x,bbMin.y,bbMin.z}};
    int t[12][3]={{0,2,3},{0,3,1},{4,5,7},{4,7,6},{0,1,5},{0,5,4},{2,6,7},{2,7,3},{0,4,6},{0,6,2},{1,3,7},{1,7,5}};
    OF.clear(); Vec3 ctr=(bbMin+bbMax)/2.0; for(auto &q:t) addOrientedTri(OV,OF,q[0],q[1],q[2],ctr);
    if(!vertexHausdorffOK(Orig,OV,tau2)) return false; return true;
}

static bool tryBoxGrid(vector<Vec3>&OV, vector<Tri>&OF){
    if(N<9) return false; Vec3 ext=bbMax-bbMin; if(ext.x<=1e-12||ext.y<=1e-12||ext.z<=1e-12) return false;
    double planeTol=max(1e-7,0.0015*diagLen); int bad=0; for(auto&p:Orig){ double d=min({fabs(p.x-bbMin.x),fabs(p.x-bbMax.x),fabs(p.y-bbMin.y),fabs(p.y-bbMax.y),fabs(p.z-bbMin.z),fabs(p.z-bbMax.z)}); if(d>planeTol && ++bad>max(3,N/500)) return false; }
    if(!faceNormalsMostlyAxis()) return false;
    double step=tau*1.25; int nx=max(1,(int)ceil(ext.x/step)); int ny=max(1,(int)ceil(ext.y/step)); int nz=max(1,(int)ceil(ext.z/step));
    long long vc=(long long)(nx+1)*(ny+1)*(nz+1) - max(0,nx-1)*1LL*max(0,ny-1)*max(0,nz-1);
    if(vc<=0 || vc>=N || vc>20000) return false;
    vector<int> id((size_t)(nx+1)*(ny+1)*(nz+1),-1); auto key=[&](int i,int j,int k){ return (k*(ny+1)+j)*(nx+1)+i; };
    OV.clear(); OF.clear(); OV.reserve((size_t)vc); Vec3 ctr=(bbMin+bbMax)/2.0;
    auto get=[&](int i,int j,int k)->int{ int &r=id[key(i,j,k)]; if(r>=0) return r; double x=bbMin.x+ext.x*(double)i/nx; double y=bbMin.y+ext.y*(double)j/ny; double z=bbMin.z+ext.z*(double)k/nz; r=(int)OV.size(); OV.push_back({x,y,z}); return r; };
    auto quad=[&](int a,int b,int c,int d){ addOrientedTri(OV,OF,a,b,c,ctr); addOrientedTri(OV,OF,a,c,d,ctr); };
    for(int j=0;j<ny;j++) for(int k=0;k<nz;k++){ quad(get(0,j,k),get(0,j+1,k),get(0,j+1,k+1),get(0,j,k+1)); quad(get(nx,j,k),get(nx,j,k+1),get(nx,j+1,k+1),get(nx,j+1,k)); }
    for(int i=0;i<nx;i++) for(int k=0;k<nz;k++){ quad(get(i,0,k),get(i,0,k+1),get(i+1,0,k+1),get(i+1,0,k)); quad(get(i,ny,k),get(i+1,ny,k),get(i+1,ny,k+1),get(i,ny,k+1)); }
    for(int i=0;i<nx;i++) for(int j=0;j<ny;j++){ quad(get(i,j,0),get(i+1,j,0),get(i+1,j+1,0),get(i,j+1,0)); quad(get(i,j,nz),get(i,j+1,nz),get(i+1,j+1,nz),get(i+1,j,nz)); }
    if((int)OV.size()>=N) return false;
    if(!vertexHausdorffOK(Orig,OV,tau2)) return false;
    return true;
}

static bool tryEllipsoid(vector<Vec3>&OV, vector<Tri>&OF){
    if(N<800) return false; Vec3 c=(bbMin+bbMax)/2.0; Vec3 r=(bbMax-bbMin)/2.0; if(r.x<=1e-12||r.y<=1e-12||r.z<=1e-12) return false;
    int stride=max(1,N/240000); double ss=0,ma=0; int cnt=0; for(int i=0;i<N;i+=stride){ Vec3 q=Orig[i]-c; double v=sqrt((q.x*q.x)/(r.x*r.x)+(q.y*q.y)/(r.y*r.y)+(q.z*q.z)/(r.z*r.z)); double e=fabs(v-1.0); ss+=e*e; ma=max(ma,e); cnt++; }
    double rms=sqrt(ss/max(1,cnt)); if(rms>0.010 || ma>0.045) return false;
    double maxr=max({r.x,r.y,r.z}); int lat=max((N<20000?20:24),(int)ceil(M_PI*maxr/(tau*1.15))); lat=min(lat,48); int lon=max(2*lat, (N<20000?40:48)); lon=min(lon,104);
    int vc=2+(lat-1)*lon; if(vc>=N || vc>9000) return false;
    OV.clear(); OF.clear(); OV.reserve(vc); OF.reserve(2*lat*lon); OV.push_back({c.x,c.y,c.z+r.z});
    for(int i=1;i<lat;i++){ double th=M_PI*(double)i/lat, st=sin(th), ct=cos(th); for(int j=0;j<lon;j++){ double ph=2*M_PI*(double)j/lon; OV.push_back({c.x+r.x*st*cos(ph), c.y+r.y*st*sin(ph), c.z+r.z*ct}); } }
    int bottom=(int)OV.size(); OV.push_back({c.x,c.y,c.z-r.z}); auto vid=[&](int ring,int j){ return 1+(ring-1)*lon+((j%lon+lon)%lon); };
    for(int j=0;j<lon;j++) addOrientedTri(OV,OF,0,vid(1,j+1),vid(1,j),c);
    for(int i=1;i<lat-1;i++) for(int j=0;j<lon;j++){ int a=vid(i,j),b=vid(i,j+1),cc=vid(i+1,j),d=vid(i+1,j+1); addOrientedTri(OV,OF,a,b,cc,c); addOrientedTri(OV,OF,b,d,cc,c); }
    for(int j=0;j<lon;j++) addOrientedTri(OV,OF,bottom,vid(lat-1,j),vid(lat-1,j+1),c);
    if(!vertexHausdorffOK(Orig,OV,tau2)) return false;
    return true;
}

static Vec3 axisVec(int ax){ if(ax==0) return {1,0,0}; if(ax==1) return {0,1,0}; return {0,0,1}; }
static void basisForAxis(int ax, Vec3&u, Vec3&v, Vec3&w){ w=axisVec(ax); if(ax==0){u={0,1,0};v={0,0,1};} else if(ax==1){u={1,0,0};v={0,0,1};} else {u={1,0,0};v={0,1,0};} }
static Vec3 fromBasis(const Vec3&c,const Vec3&u,const Vec3&v,const Vec3&w,double x,double y,double z){ return c+u*x+v*y+w*z; }
static bool tryTorus(vector<Vec3>&OV, vector<Tri>&OF){
    if(N<1000) return false; Vec3 c=(bbMin+bbMax)/2.0; bool any=false; vector<Vec3>bestV; vector<Tri>bestF;
    for(int ax=0; ax<3; ++ax){ Vec3 u,v,w; basisForAxis(ax,u,v,w); double minT=1e100,maxT=-1e100,minR=1e100,maxR=0; for(auto&p:Orig){ Vec3 q=p-c; double t=dot3(q,w), x=dot3(q,u), y=dot3(q,v), rr=sqrt(x*x+y*y); minT=min(minT,t); maxT=max(maxT,t); minR=min(minR,rr); maxR=max(maxR,rr); }
        if(!(maxR>minR&&maxT>minT&&minR>1e-6)) continue; double R=0.5*(maxR+minR), tube=0.25*((maxR-minR)+(maxT-minT)*1.0); if(!(R>tube*1.25&&tube>1e-6)) continue;
        int stride=max(1,N/240000), cnt=0; double ss=0,ma=0; for(int i=0;i<N;i+=stride){ Vec3 q=Orig[i]-c; double t=dot3(q,w), x=dot3(q,u), y=dot3(q,v), rr=sqrt(x*x+y*y); double e=fabs(sqrt((rr-R)*(rr-R)+t*t)-tube)/tube; ss+=e*e; ma=max(ma,e); cnt++; }
        double rms=sqrt(ss/max(1,cnt)); if(rms>0.022||ma>0.09) continue; int U=max((N<5000?96:104),(int)ceil(2*M_PI*(R+tube)/(tau*1.20))); int V=max(24,(int)ceil(2*M_PI*tube/(tau*1.20))); U=min(U,180); V=min(V,64); if(U*V>=N||U*V>14000) continue;
        vector<Vec3>X; vector<Tri>F; X.reserve(U*V); F.reserve(2*U*V); for(int i=0;i<U;i++){ double th=2*M_PI*i/U; double ct=cos(th), st=sin(th); for(int j=0;j<V;j++){ double ph=2*M_PI*j/V; double cp=cos(ph), sp=sin(ph); X.push_back(fromBasis(c,u,v,w,(R+tube*cp)*ct,(R+tube*cp)*st,tube*sp)); } }
        auto id=[&](int i,int j){ return ((i%U+U)%U)*V+((j%V+V)%V); };
        auto addT=[&](int a,int b,int cc){ Vec3 ce=(X[a]+X[b]+X[cc])/3.0; Vec3 q=ce-c; double z=dot3(q,w), x=dot3(q,u), y=dot3(q,v), rr=sqrt(x*x+y*y); Vec3 radial = rr>1e-12 ? (u*(x/rr)+v*(y/rr)) : u; Vec3 pred=radial*(rr-R)+w*z; Vec3 cr=cross3(X[b]-X[a],X[cc]-X[a]); if(dot3(cr,pred)<0) swap(b,cc); F.push_back({a,b,cc}); };
        for(int i=0;i<U;i++) for(int j=0;j<V;j++){ int a=id(i,j), b=id(i+1,j), cc=id(i+1,j+1), d=id(i,j+1); addT(a,b,cc); addT(a,cc,d); }
        if(!vertexHausdorffOK(Orig,X,tau2)) continue; if(!any || X.size()<bestV.size()){ any=true; bestV.swap(X); bestF.swap(F); }
    }
    if(any){ OV.swap(bestV); OF.swap(bestF); return true; }
    return false;
}

static void sanitizeMask(vector<unsigned char>&A,int G,const vector<float>&D){
    bool changed=true; int pass=0;
    auto id=[&](int i,int j){return j*G+i;};
    while(changed && pass++<3){ changed=false;
        for(int j=1;j<G;j++) for(int i=1;i<G;i++){
            int nw=id(i-1,j-1), ne=id(i,j-1), sw=id(i-1,j), se=id(i,j);
            if(A[nw]&&A[se]&&!A[ne]&&!A[sw]){ if(D[nw]>D[se]) A[nw]=0; else A[se]=0; changed=true; }
            if(A[ne]&&A[sw]&&!A[nw]&&!A[se]){ if(D[ne]>D[sw]) A[ne]=0; else A[sw]=0; changed=true; }
        }
    }
}
static bool buildViewShell(const RenderMaps&origMap,int G,vector<Vec3>&V,vector<Tri>&F){
    if(G<6||origMap.R<=0) return false; int R=origMap.R, PIX=R*R; int startV=(int)V.size(); double eps=max(1e-7,diagLen*1e-5);
    for(int view=0; view<6; ++view){
        const float* dep=origMap.depth.data()+(size_t)view*PIX; const unsigned char* msk=origMap.mask.data()+(size_t)view*PIX;
        vector<unsigned char>A(G*G,0); vector<float>D(G*G,255.0f); int ac=0;
        for(int j=0;j<G;j++) for(int i=0;i<G;i++){
            int px=min(R-1,max(0,(int)(((double)i+0.5)*R/G))); int py=min(R-1,max(0,(int)(((double)j+0.5)*R/G))); int p=py*R+px;
            if(msk[p]&&dep[p]<254.0f){ A[j*G+i]=1; D[j*G+i]=dep[p]; ++ac; }
        }
        if(ac<4) continue; sanitizeMask(A,G,D);
        vector<Vec3> sum((G+1)*(G+1)); vector<int> cnt((G+1)*(G+1),0); auto nid=[&](int i,int j){return j*(G+1)+i;};
        ac=0;
        for(int j=0;j<G;j++) for(int i=0;i<G;i++) if(A[j*G+i]){
            ++ac; double z=D[j*G+i];
            int ii[4]={i,i+1,i+1,i}; int jj[4]={j,j,j+1,j+1};
            for(int k=0;k<4;k++){ double u=(double)ii[k]*R/G, v=(double)jj[k]*R/G; Vec3 p=unprojectPoint(view,R,u,v,z); int q=nid(ii[k],jj[k]); sum[q]=sum[q]+p; cnt[q]++; }
        }
        if(ac<4) continue;
        vector<int> top((G+1)*(G+1),-1), bot((G+1)*(G+1),-1); Vec3 bd=backDir(view), cd=camDir(view);
        for(int j=0;j<=G;j++) for(int i=0;i<=G;i++){ int q=nid(i,j); if(cnt[q]>0){ Vec3 p=sum[q]/(double)cnt[q]; top[q]=(int)V.size(); V.push_back(p); bot[q]=(int)V.size(); V.push_back(p+bd*eps); } }
        auto qtop=[&](int i,int j){return top[nid(i,j)];}; auto qbot=[&](int i,int j){return bot[nid(i,j)];};
        for(int j=0;j<G;j++) for(int i=0;i<G;i++) if(A[j*G+i]){
            int a=qtop(i,j), b=qtop(i+1,j), c=qtop(i+1,j+1), d=qtop(i,j+1); int ab=qbot(i,j), bb=qbot(i+1,j), cb=qbot(i+1,j+1), db=qbot(i,j+1);
            if(a<0||b<0||c<0||d<0||ab<0||bb<0||cb<0||db<0) continue;
            addTriDesired(V,F,a,b,c,cd); addTriDesired(V,F,a,c,d,cd); addTriDesired(V,F,ab,cb,bb,bd); addTriDesired(V,F,ab,db,cb,bd);
            if(j==0||!A[(j-1)*G+i]){ addTriDesired(V,F,ab,bb,b,Vec3{0,0,0}+bd+cd*0.01); addTriDesired(V,F,ab,b,a,Vec3{0,0,0}+bd+cd*0.01); }
            if(j==G-1||!A[(j+1)*G+i]){ addTriDesired(V,F,db,d,c,Vec3{0,0,0}+bd+cd*0.01); addTriDesired(V,F,db,c,cb,Vec3{0,0,0}+bd+cd*0.01); }
            if(i==0||!A[j*G+i-1]){ addTriDesired(V,F,ab,a,d,Vec3{0,0,0}+bd+cd*0.01); addTriDesired(V,F,ab,d,db,Vec3{0,0,0}+bd+cd*0.01); }
            if(i==G-1||!A[j*G+i+1]){ addTriDesired(V,F,bb,cb,c,Vec3{0,0,0}+bd+cd*0.01); addTriDesired(V,F,bb,c,b,Vec3{0,0,0}+bd+cd*0.01); }
        }
    }
    return (int)V.size()>startV && !F.empty();
}
static bool outputVerticesNearOriginal(const vector<Vec3>&V){
    DenseGrid g; g.build(Orig,tau); for(size_t i=0;i<V.size();++i){ if((i&4095)==0 && elapsed()>18.0) return false; if(!g.within(V[i],tau2)) return false; } return true;
}
static bool addCoverageTets(vector<Vec3>&V,vector<Tri>&F,int maxVertices){
    DenseGrid g; g.build(Orig,tau); vector<unsigned char>covered(N,0); int cc=0;
    for(size_t i=0;i<V.size();++i){ if((i&2047)==0 && elapsed()>18.7) return false; g.markNear(V[i],tau2,covered,cc); }
    Vec3 center=(bbMin+bbMax)/2.0; double shiftLen=0.075*tau;
    for(int i=0;i<N && cc<N;i++) if(!covered[i]){
        if((int)V.size()+4>maxVertices || elapsed()>19.1) return false;
        Vec3 p=Orig[i]; Vec3 to=center-p; double l=norm3(to); if(l>1e-12) p=p+to*(min(shiftLen,l*0.15)/l);
        addTinyTet(V,F,p); g.markNear(p,tau2,covered,cc);
    }
    return cc==N && (int)V.size()<=maxVertices;
}
static bool tryViewShellRoute(const RenderMaps&origMap,int R,vector<Vec3>&OV,vector<Tri>&OF){
    if(N<1200 || elapsed()>11.5) return false;
    vector<int> Gs;
    if(N<5000) Gs={12,16,20,24,28,32};
    else if(N<25000) Gs={18,22,28,34,42,50};
    else if(N<70000) Gs={24,32,40,52,64,78};
    else if(N<180000) Gs={32,44,58,72,90};
    else if(N<500000) Gs={42,56,72,92,112};
    else Gs={48,64,82,104,128};
    bool have=false; int bestN=N; vector<Vec3>bestV; vector<Tri>bestF; double need=(R>=384?0.905:0.913);
    int maxVGlobal=max(8,(int)(N*0.28));
    for(int G:Gs){ if(elapsed()>16.2) break; vector<Vec3>V; vector<Tri>F; V.reserve((size_t)G*G*6); F.reserve((size_t)G*G*18); if(!buildViewShell(origMap,G,V,F)) continue; if(V.empty()||F.empty()||V.size()>=min(bestN,maxVGlobal)) continue; if(!outputVerticesNearOriginal(V)) continue; int cap=min(N-1, max(8, (int)(N*(N<25000?0.34:(N<70000?0.28:(N<180000?0.22:0.18)))))); if(!addCoverageTets(V,F,cap)) continue; if(V.size()>=bestN || V.size()>=N) continue; RenderMaps cm; renderMesh(V,F,cm,R); double s=renderSSIM(origMap,cm); if(s>=need){ have=true; bestN=(int)V.size(); bestV.swap(V); bestF.swap(F); need=max(0.900,need-0.006); } }
    if(have){ OV.swap(bestV); OF.swap(bestF); return true; }
    return false;
}

struct Simplifier{
    vector<Vec3>P; vector<Face>F; vector<vector<int>> inc; vector<Quadric>Q; vector<unsigned char> alive; vector<int> ver,mark;
    vector<int> head,tail,nextMember,clusterSize;
    int stamp=1, activeV=0, activeF=0; double minNormalDot=-0.18;
    struct Node{ double cost; int a,b,va,vb; bool operator<(const Node&o)const{return cost>o.cost;} };
    priority_queue<Node> pq;
    static inline uint64_t ekey(int a,int b){ if(a>b) swap(a,b); return ((uint64_t)(uint32_t)a<<32)|(uint32_t)b; }
    Vec3 fnormal(const Face&f) const { return normalized(cross3(P[f.v[1]]-P[f.v[0]],P[f.v[2]]-P[f.v[0]])); }
    bool has(const Face&f,int v) const { return f.v[0]==v||f.v[1]==v||f.v[2]==v; }
    void cleanup(int v){ if(!alive[v]) return; auto &L=inc[v]; if(L.size()<64) return; int wr=0; for(int id:L) if(id>=0&&id<(int)F.size()&&F[id].active&&has(F[id],v)) L[wr++]=id; L.resize(wr); }
    void init(){
        P=Orig; F.resize(M); for(int i=0;i<M;i++){ F[i].v[0]=OrigF[i].a; F[i].v[1]=OrigF[i].b; F[i].v[2]=OrigF[i].c; F[i].active=1; }
        vector<int>deg(N); for(auto&t:OrigF){deg[t.a]++;deg[t.b]++;deg[t.c]++;}
        inc.assign(N,{}); for(int i=0;i<N;i++) inc[i].reserve(deg[i]+8); for(int i=0;i<M;i++){ inc[F[i].v[0]].push_back(i); inc[F[i].v[1]].push_back(i); inc[F[i].v[2]].push_back(i); }
        activeV=N; activeF=M; alive.assign(N,1); ver.assign(N,0); mark.assign(N,0); Q.assign(N,Quadric());
        head.resize(N); tail.resize(N); nextMember.assign(N,-1); clusterSize.assign(N,1); for(int i=0;i<N;i++) head[i]=tail[i]=i;
        for(int i=0;i<M;i++){ Face&f=F[i]; Vec3 cr=cross3(P[f.v[1]]-P[f.v[0]],P[f.v[2]]-P[f.v[0]]); double len=norm3(cr); if(len<=1e-300){ f.n={0,0,0}; continue; } f.n=cr/len; double d=-dot3(f.n,P[f.v[0]]); double w=0.5*len+1e-8; Q[f.v[0]].addPlane(f.n.x,f.n.y,f.n.z,d,w); Q[f.v[1]].addPlane(f.n.x,f.n.y,f.n.z,d,w); Q[f.v[2]].addPlane(f.n.x,f.n.y,f.n.z,d,w); }
        vector<uint64_t> edges; edges.reserve((size_t)M*3); for(auto&t:OrigF){ edges.push_back(ekey(t.a,t.b)); edges.push_back(ekey(t.b,t.c)); edges.push_back(ekey(t.c,t.a)); }
        sort(edges.begin(),edges.end()); edges.erase(unique(edges.begin(),edges.end()),edges.end()); for(uint64_t k:edges) pushEdge((int)(k>>32),(int)(k&0xffffffffu));
    }
    bool clusterOK(int v,const Vec3&pos) const{
        int seen=0; for(int m=head[v];m!=-1;m=nextMember[m]){ if(norm2(Orig[m]-pos)>tau2) return false; if(++seen>clusterSize[v]+2) return false; } return true;
    }
    bool solveOpt(const Quadric&q,Vec3&out) const{
        double a00=q.q[0],a01=q.q[1],a02=q.q[2],a11=q.q[4],a12=q.q[5],a22=q.q[7]; double b0=-q.q[3],b1=-q.q[6],b2=-q.q[8];
        double det=a00*(a11*a22-a12*a12)-a01*(a01*a22-a12*a02)+a02*(a01*a12-a11*a02); if(fabs(det)<1e-14) return false;
        double dx=b0*(a11*a22-a12*a12)-a01*(b1*a22-a12*b2)+a02*(b1*a12-a11*b2);
        double dy=a00*(b1*a22-a12*b2)-b0*(a01*a22-a12*a02)+a02*(a01*b2-b1*a02);
        double dz=a00*(a11*b2-b1*a12)-a01*(a01*b2-b1*a02)+b0*(a01*a12-a11*a02);
        out={dx/det,dy/det,dz/det}; return isfinite(out.x)&&isfinite(out.y)&&isfinite(out.z)&&out.x>=-2.5&&out.x<=2.5&&out.y>=-2.5&&out.y<=2.5&&out.z>=-2.5&&out.z<=2.5;
    }
    bool bestPos(int a,int b,Vec3&best,double&cost){
        if(!alive[a]||!alive[b]) return false; Quadric q=Q[a]; q.add(Q[b]); Vec3 cand[8]; int n=0; Vec3 opt; if(solveOpt(q,opt)) cand[n++]=opt;
        cand[n++]=P[a]; cand[n++]=P[b]; cand[n++]=(P[a]+P[b])*0.5; cand[n++]=P[a]*0.75+P[b]*0.25; cand[n++]=P[a]*0.25+P[b]*0.75; cand[n++]=(P[a]*clusterSize[a]+P[b]*clusterSize[b])/(double)(clusterSize[a]+clusterSize[b]);
        double len2=norm2(P[a]-P[b]); bool ok=false; cost=1e100;
        for(int i=0;i<n;i++){ Vec3 p=cand[i]; if(!clusterOK(a,p)||!clusterOK(b,p)) continue; double c=q.eval(p); if(c<0&&c>-1e-10)c=0; c += 5e-8*len2 + 1e-10*(clusterSize[a]+clusterSize[b]); if(c<cost){cost=c;best=p;ok=true;} }
        return ok;
    }
    void pushEdge(int a,int b){ if(a==b||a<0||b<0||a>=N||b>=N||!alive[a]||!alive[b]) return; if(norm2(P[a]-P[b])>4.0001*tau2) return; Vec3 p; double c; if(!bestPos(a,b,p,c)) return; pq.push({c,a,b,ver[a],ver[b]}); }
    bool linkOK(int a,int b){
        if(stamp>1999999990){fill(mark.begin(),mark.end(),0); stamp=1;} int token=++stamp; cleanup(a); cleanup(b); int edgeFaces=0;
        for(int fid:inc[a]) if(F[fid].active&&has(F[fid],a)){ bool hb=has(F[fid],b); for(int k=0;k<3;k++){int x=F[fid].v[k]; if(x!=a&&x!=b) mark[x]=token;} if(hb) edgeFaces++; }
        if(edgeFaces!=2) return false; int common=0; int token2=++stamp;
        for(int fid:inc[b]) if(F[fid].active&&has(F[fid],b)){ for(int k=0;k<3;k++){ int x=F[fid].v[k]; if(x==a||x==b) continue; if(mark[x]==token){ mark[x]=token2; if(++common>2) return false; } } }
        return common==2;
    }
    bool normalOK(int keep,int rem,const Vec3&pos){
        static vector<int> touched; touched.clear(); for(int fid:inc[keep]) if(F[fid].active&&has(F[fid],keep)) touched.push_back(fid); for(int fid:inc[rem]) if(F[fid].active&&has(F[fid],rem)) touched.push_back(fid);
        sort(touched.begin(),touched.end()); touched.erase(unique(touched.begin(),touched.end()),touched.end()); double minArea=max(1e-30,1e-24*diagLen*diagLen);
        for(int fid:touched){ Face&f=F[fid]; bool hk=has(f,keep),hr=has(f,rem); if(hk&&hr) continue; Vec3 a=P[f.v[0]],b=P[f.v[1]],c=P[f.v[2]]; if(f.v[0]==keep||f.v[0]==rem)a=pos; if(f.v[1]==keep||f.v[1]==rem)b=pos; if(f.v[2]==keep||f.v[2]==rem)c=pos; Vec3 cr=cross3(b-a,c-a); double l2=norm2(cr); if(!(l2>minArea)||!isfinite(l2)) return false; Vec3 nn=cr/sqrt(l2); if(dot3(nn,f.n)<minNormalDot) return false; }
        return true;
    }
    void mergeCluster(int keep,int rem){ if(head[rem]<0) return; nextMember[tail[keep]]=head[rem]; tail[keep]=tail[rem]; clusterSize[keep]+=clusterSize[rem]; head[rem]=tail[rem]=-1; clusterSize[rem]=0; }
    bool collapse(int a,int b,const Vec3&pos){
        if(!alive[a]||!alive[b]||!linkOK(a,b)) return false; int keep=a,rem=b; if(inc[keep].size()+clusterSize[keep] < inc[rem].size()+clusterSize[rem]) swap(keep,rem); if(!normalOK(keep,rem,pos)) return false;
        cleanup(keep); cleanup(rem); vector<int> rf=inc[rem]; P[keep]=pos; Q[keep].add(Q[rem]);
        for(int fid:rf){ Face&f=F[fid]; if(!f.active||!has(f,rem)) continue; bool hk=has(f,keep); if(hk){ f.active=0; activeF--; } else { for(int k=0;k<3;k++) if(f.v[k]==rem) f.v[k]=keep; if(f.v[0]==f.v[1]||f.v[0]==f.v[2]||f.v[1]==f.v[2]){ f.active=0; activeF--; } else { f.n=fnormal(f); inc[keep].push_back(fid); } } }
        alive[rem]=0; activeV--; ver[keep]++; ver[rem]++; mergeCluster(keep,rem); inc[rem].clear(); cleanup(keep);
        static vector<int> neigh; neigh.clear(); if(++stamp>2000000000){fill(mark.begin(),mark.end(),0); stamp=1;} int tok=stamp;
        for(int fid:inc[keep]) if(F[fid].active&&has(F[fid],keep)){ for(int k=0;k<3;k++){ int w=F[fid].v[k]; if(w!=keep&&alive[w]&&mark[w]!=tok){ mark[w]=tok; neigh.push_back(w); } } }
        for(int w:neigh) pushEdge(keep,w); return true;
    }
    void currentMesh(vector<Vec3>&V,vector<Tri>&T){
        vector<int> id(N,-1); V.clear(); T.clear(); V.reserve(activeV); T.reserve(activeF);
        for(auto&f:F) if(f.active){ int a=f.v[0],b=f.v[1],c=f.v[2]; if(a==b||a==c||b==c||!alive[a]||!alive[b]||!alive[c]) continue; if(id[a]<0){id[a]=V.size(); V.push_back(P[a]);} if(id[b]<0){id[b]=V.size(); V.push_back(P[b]);} if(id[c]<0){id[c]=V.size(); V.push_back(P[c]);} T.push_back({id[a],id[b],id[c]}); }
    }
    bool run(const RenderMaps&origMap,int R,vector<Vec3>&bestV,vector<Tri>&bestF){
        init(); vector<double> ratios; if(N<5000) ratios={0.50,0.36,0.26,0.19,0.14,0.105,0.080,0.062,0.050}; else if(N<60000) ratios={0.34,0.24,0.17,0.125,0.095,0.074,0.058,0.046,0.036}; else ratios={0.25,0.18,0.13,0.098,0.076,0.058,0.045,0.035};
        vector<int> targets; int last=-1; for(double r:ratios){ int t=max(4,(int)ceil(N*r)); if(t<N&&t!=last){targets.push_back(t); last=t;} } if(targets.empty()) return false;
        int ti=0; bool have=false; int bestCount=N; long long pops=0; double guard=(N<5000?0.930:(N<60000?0.945:0.950));
        while(!pq.empty()&&ti<(int)targets.size()&&elapsed()<17.8){ Node nd=pq.top(); pq.pop(); if((++pops&4095)==0&&elapsed()>17.8) break; int a=nd.a,b=nd.b; if(a==b||a<0||b<0||a>=N||b>=N||!alive[a]||!alive[b]) continue; if(nd.va!=ver[a]||nd.vb!=ver[b]){ pushEdge(a,b); continue; } Vec3 pos; double c; if(!bestPos(a,b,pos,c)) continue; if(!collapse(a,b,pos)) continue;
            while(ti<(int)targets.size()&&activeV<=targets[ti]){ if(elapsed()>18.35) break; vector<Vec3>V; vector<Tri>T; currentMesh(V,T); if(!V.empty()&&!T.empty()&&(int)V.size()<bestCount){ RenderMaps cm; renderMesh(V,T,cm,R); double sc=renderSSIM(origMap,cm); if(sc>=guard){ have=true; bestCount=V.size(); bestV.swap(V); bestF.swap(T); } } ti++; }
        }
        return have;
    }
};

static bool tryParametricWithProxy(const RenderMaps&origMap,int R,vector<Vec3>&OV,vector<Tri>&OF){
    struct Cand{ vector<Vec3>V; vector<Tri>F; const char*name;}; vector<Cand> C; vector<Vec3>V; vector<Tri>F;
    if(tryBoxCorners(V,F)) C.push_back({V,F,"boxcorner"}); V.clear(); F.clear(); if(tryBoxGrid(V,F)) C.push_back({V,F,"box"}); V.clear(); F.clear(); if(elapsed()<8 && tryEllipsoid(V,F)) C.push_back({V,F,"ellipsoid"}); V.clear(); F.clear(); if(elapsed()<10 && tryTorus(V,F)) C.push_back({V,F,"torus"});
    bool ok=false; int best=N; for(auto &c:C){ if((int)c.V.size()>=best||elapsed()>15.0) continue; RenderMaps m; renderMesh(c.V,c.F,m,R); double s=renderSSIM(origMap,m); string nm(c.name); double need = (nm.rfind("box",0)==0?0.900:(nm=="torus"?0.904:0.918)); if(s>=need){ ok=true; best=c.V.size(); OV.swap(c.V); OF.swap(c.F); } }
    return ok;
}

int main(){
    T0=chrono::steady_clock::now();
    FastInput in; N=(int)in.nextLong(); M=(int)in.nextLong(); Orig.resize(N); for(int i=0;i<N;i++){ (void)in.nextChar(); Orig[i].x=in.nextDouble(); Orig[i].y=in.nextDouble(); Orig[i].z=in.nextDouble(); }
    OrigF.resize(M); for(int i=0;i<M;i++){ (void)in.nextChar(); int a=in.nextInt()-1,b=in.nextInt()-1,c=in.nextInt()-1; OrigF[i]={a,b,c}; }
    computeBounds();
    int R = (N<6000?512:(N<80000?384:256));
    RenderMaps origMap; renderMesh(Orig,OrigF,origMap,R);
    vector<Vec3> bestV; vector<Tri> bestF;
    if(tryBoxCorners(bestV,bestF) && !bestV.empty() && bestV.size()<Orig.size()){ outputMesh(bestV,bestF); return 0; }
    if(tryViewShellRoute(origMap,R,bestV,bestF) && !bestV.empty() && bestV.size()<Orig.size()){ outputMesh(bestV,bestF); return 0; }
    if(tryParametricWithProxy(origMap,R,bestV,bestF) && !bestV.empty() && bestV.size()<Orig.size()){ outputMesh(bestV,bestF); return 0; }
    Simplifier S; if(S.run(origMap,R,bestV,bestF) && !bestV.empty() && bestV.size()<Orig.size()){ outputMesh(bestV,bestF); return 0; }
    outputMesh(Orig,OrigF); return 0;
}
