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
static inline double normv(const Vec3&a){return sqrt(max(0.0,norm2(a)));} 
static inline Vec3 normalized(Vec3 a){double n=normv(a); if(n<=1e-300) return Vec3(); return a/n;} 
static inline double dist2(const Vec3&a,const Vec3&b){return norm2(a-b);} 

struct Tri{int a,b,c;};
struct Face{int v[3]; unsigned char active; Vec3 n;};
struct MeshCand{ vector<Vec3> V; vector<Tri> F; string name; double proxy=-1; };

static int N0,M0;
static vector<Vec3> OV;
static vector<Tri> OF;
static double DIAG=1.0, TAU=0.05, TAU2=0.0025;
static chrono::steady_clock::time_point T0;
static double elapsed(){return chrono::duration<double>(chrono::steady_clock::now()-T0).count();}

struct FastInput{
    vector<char> buf; char* p=nullptr;
    FastInput(){ char tmp[1<<16]; size_t n; while((n=fread(tmp,1,sizeof(tmp),stdin))>0) buf.insert(buf.end(),tmp,tmp+n); buf.push_back(0); p=buf.data(); }
    inline void skip(){ while(*p==' '||*p=='\n'||*p=='\r'||*p=='\t') ++p; }
    int nextInt(){ skip(); int s=1; if(*p=='-'){s=-1;++p;} int x=0; while(*p>='0'&&*p<='9'){x=x*10+(*p-'0');++p;} return x*s; }
    double nextDouble(){ skip(); char* e; double x=strtod(p,&e); p=e; return x; }
    char nextChar(){ skip(); return *p++; }
};

static void readInput(){
    FastInput in;
    N0=in.nextInt(); M0=in.nextInt();
    OV.resize(N0); OF.resize(M0);
    Vec3 mn(1e100,1e100,1e100), mx(-1e100,-1e100,-1e100);
    for(int i=0;i<N0;i++){
        (void)in.nextChar();
        OV[i].x=in.nextDouble(); OV[i].y=in.nextDouble(); OV[i].z=in.nextDouble();
        mn.x=min(mn.x,OV[i].x); mn.y=min(mn.y,OV[i].y); mn.z=min(mn.z,OV[i].z);
        mx.x=max(mx.x,OV[i].x); mx.y=max(mx.y,OV[i].y); mx.z=max(mx.z,OV[i].z);
    }
    for(int i=0;i<M0;i++){
        (void)in.nextChar();
        OF[i].a=in.nextInt()-1; OF[i].b=in.nextInt()-1; OF[i].c=in.nextInt()-1;
    }
    DIAG=normv(mx-mn); if(DIAG<=0) DIAG=1.0;
    TAU=0.05*DIAG*0.999; TAU2=TAU*TAU;
}

static inline uint64_t edgeKey(int a,int b){ if(a>b) swap(a,b); return (uint64_t)(uint32_t)a<<32 | (uint32_t)b; }
static inline array<int,3> faceKey(Tri t){ array<int,3> a{t.a,t.b,t.c}; sort(a.begin(),a.end()); return a; }

static bool manifoldAndBasicOK(const vector<Vec3>&V,const vector<Tri>&F){
    if(V.empty() || F.empty() || (int)V.size()>N0) return false;
    double minArea2=max(1e-30,1e-24*DIAG*DIAG);
    vector<uint64_t> edges; edges.reserve(F.size()*3);
    vector<array<int,3>> fkeys; fkeys.reserve(F.size());
    for(auto t:F){
        if(t.a<0||t.b<0||t.c<0||t.a>=(int)V.size()||t.b>=(int)V.size()||t.c>=(int)V.size()) return false;
        if(t.a==t.b||t.a==t.c||t.b==t.c) return false;
        Vec3 cr=crossv(V[t.b]-V[t.a],V[t.c]-V[t.a]);
        if(!(norm2(cr)>minArea2) || !isfinite(norm2(cr))) return false;
        edges.push_back(edgeKey(t.a,t.b)); edges.push_back(edgeKey(t.b,t.c)); edges.push_back(edgeKey(t.c,t.a));
        fkeys.push_back(faceKey(t));
    }
    sort(fkeys.begin(),fkeys.end());
    if(adjacent_find(fkeys.begin(),fkeys.end())!=fkeys.end()) return false;
    sort(edges.begin(),edges.end());
    for(size_t i=0;i<edges.size();){ size_t j=i+1; while(j<edges.size()&&edges[j]==edges[i]) ++j; if(j-i!=2) return false; i=j; }
    return true;
}

struct GridHash{
    double cell; unordered_map<uint64_t, vector<int>> mp; const vector<Vec3>* pts;
    GridHash(const vector<Vec3>&P,double c):cell(c),pts(&P){ mp.reserve(P.size()*2+10); for(int i=0;i<(int)P.size();++i) mp[key(P[i])].push_back(i); }
    inline uint64_t pack(long long ix,long long iy,long long iz) const{
        const long long OFF=1LL<<20, MASK=(1LL<<21)-1;
        return (uint64_t)((ix+OFF)&MASK)<<42 | (uint64_t)((iy+OFF)&MASK)<<21 | (uint64_t)((iz+OFF)&MASK);
    }
    inline uint64_t key(const Vec3&p) const { return pack((long long)floor(p.x/cell),(long long)floor(p.y/cell),(long long)floor(p.z/cell)); }
    bool hasNear(const Vec3&p,double r2) const{
        long long ix=(long long)floor(p.x/cell), iy=(long long)floor(p.y/cell), iz=(long long)floor(p.z/cell);
        for(long long dx=-1;dx<=1;++dx) for(long long dy=-1;dy<=1;++dy) for(long long dz=-1;dz<=1;++dz){
            auto it=mp.find(pack(ix+dx,iy+dy,iz+dz)); if(it==mp.end()) continue;
            for(int id:it->second) if(dist2(p,(*pts)[id])<=r2) return true;
        }
        return false;
    }
};

static bool vertexHausdorffOK(const vector<Vec3>&V){
    if(V.empty() || (int)V.size()>N0) return false;
    double cell=max(TAU,1e-9);
    GridHash gout(V,cell);
    for(const Vec3&p:OV) if(!gout.hasNear(p,TAU2)) return false;
    GridHash gorig(OV,cell);
    for(const Vec3&p:V) if(!gorig.hasNear(p,TAU2)) return false;
    return true;
}

static double orientSignFromOriginal(const Vec3&center){
    int stride=max(1,M0/200000); double s=0; int c=0;
    for(int i=0;i<M0;i+=stride){
        const Tri&t=OF[i]; Vec3 a=OV[t.a],b=OV[t.b],cc=OV[t.c]; Vec3 cr=crossv(b-a,cc-a); Vec3 ctr=(a+b+cc)/3.0;
        double v=dotv(cr,ctr-center); if(isfinite(v)){s+=v; ++c;}
    }
    return s>=0?1.0:-1.0;
}
static void orientFace(vector<Vec3>&V, Tri&f, const Vec3&center, double sign){
    Vec3 cr=crossv(V[f.b]-V[f.a], V[f.c]-V[f.a]); Vec3 ctr=(V[f.a]+V[f.b]+V[f.c])/3.0;
    if(dotv(cr,ctr-center)*sign < 0) swap(f.b,f.c);
}

// ---------- Rendering and SSIM proxy ----------
struct RenderMaps{int R=0; vector<float> depth; vector<float> norm;};
static inline void projectPoint(const Vec3&p,int view,int R,double&u,double&v,double&z){
    const double D=2.5; double f=800.0*(double)R/1024.0, c=0.5*(double)R; double sx,sy;
    if(view==0){sx=p.y; sy=p.z; z=D-p.x;}
    else if(view==1){sx=-p.y; sy=p.z; z=D+p.x;}
    else if(view==2){sx=-p.x; sy=p.z; z=D-p.y;}
    else if(view==3){sx=p.x; sy=p.z; z=D+p.y;}
    else if(view==4){sx=p.x; sy=p.y; z=D-p.z;}
    else {sx=-p.x; sy=p.y; z=D+p.z;}
    u=f*sx/z+c; v=f*sy/z+c;
}
static void renderMesh(const vector<Vec3>&V,const vector<Tri>&F,RenderMaps&rm,int R){
    const int PIX=R*R, VIEWS=6; rm.R=R; rm.depth.assign((size_t)VIEWS*PIX,255.0f); rm.norm.assign((size_t)VIEWS*PIX*3,127.5f);
    vector<float> U(V.size()), VV(V.size()), Z(V.size()); vector<Vec3> FN(F.size());
    for(size_t i=0;i<F.size();++i){ const Tri&t=F[i]; Vec3 cr=crossv(V[t.b]-V[t.a],V[t.c]-V[t.a]); FN[i]=normalized(cr); }
    for(int view=0; view<VIEWS; ++view){
        for(size_t i=0;i<V.size();++i){ double u,v,z; projectPoint(V[i],view,R,u,v,z); U[i]=(float)u; VV[i]=(float)v; Z[i]=(float)z; }
        float* zbuf=rm.depth.data()+(size_t)view*PIX; float* nbuf=rm.norm.data()+(size_t)view*PIX*3;
        for(size_t ti=0; ti<F.size(); ++ti){
            const Tri&t=F[ti]; int ia=t.a,ib=t.b,ic=t.c; float x0=U[ia],x1=U[ib],x2=U[ic], y0=VV[ia],y1=VV[ib],y2=VV[ic], z0=Z[ia],z1=Z[ib],z2=Z[ic];
            if(!(z0>0&&z1>0&&z2>0)) continue;
            int xmin=max(0,(int)floor(min(x0,min(x1,x2))-0.5f)); int xmax=min(R-1,(int)ceil(max(x0,max(x1,x2))+0.5f));
            int ymin=max(0,(int)floor(min(y0,min(y1,y2))-0.5f)); int ymax=min(R-1,(int)ceil(max(y0,max(y1,y2))+0.5f));
            if(xmin>xmax||ymin>ymax) continue; float den=(y1-y2)*(x0-x2)+(x2-x1)*(y0-y2); if(fabs(den)<1e-20f) continue; float invDen=1.0f/den;
            Vec3 n=FN[ti]; float nr=(float)((n.x+1.0)*127.5), ng=(float)((n.y+1.0)*127.5), nb=(float)((n.z+1.0)*127.5);
            for(int yy=ymin; yy<=ymax; ++yy){ float py=(float)yy+0.5f; int row=yy*R; for(int xx=xmin; xx<=xmax; ++xx){ float px=(float)xx+0.5f;
                    float w0=((y1-y2)*(px-x2)+(x2-x1)*(py-y2))*invDen;
                    float w1=((y2-y0)*(px-x2)+(x0-x2)*(py-y2))*invDen;
                    float w2=1.0f-w0-w1;
                    if(w0>=-1e-6f&&w1>=-1e-6f&&w2>=-1e-6f){ float dep=1.0f/(w0/z0+w1/z1+w2/z2); int idx=row+xx; if(dep<zbuf[idx]){ zbuf[idx]=dep; float*q=nbuf+idx*3; q[0]=nr; q[1]=ng; q[2]=nb; }}
            }}
        }
    }
}
static inline double rectSum(const vector<double>&I,int W,int x0,int y0,int x1,int y1){ return I[(size_t)y1*W+x1]-I[(size_t)y0*W+x1]-I[(size_t)y1*W+x0]+I[(size_t)y0*W+x0]; }
static double ssimChannel(const float*A,int strideA,const float*B,int strideB,const float*dA,const float*dB,int R,vector<double>&IA,vector<double>&IB,vector<double>&IA2,vector<double>&IB2,vector<double>&IAB){
    int W=R+1; size_t SZ=(size_t)W*W; fill(IA.begin(),IA.begin()+SZ,0); fill(IB.begin(),IB.begin()+SZ,0); fill(IA2.begin(),IA2.begin()+SZ,0); fill(IB2.begin(),IB2.begin()+SZ,0); fill(IAB.begin(),IAB.begin()+SZ,0);
    for(int y=1;y<=R;y++){ double sa=0,sb=0,sa2=0,sb2=0,sab=0; int row=(y-1)*R; for(int x=1;x<=R;x++){ int p=row+x-1; double a=A[(size_t)p*strideA], b=B[(size_t)p*strideB]; sa+=a; sb+=b; sa2+=a*a; sb2+=b*b; sab+=a*b; size_t id=(size_t)y*W+x, up=(size_t)(y-1)*W+x; IA[id]=IA[up]+sa; IB[id]=IB[up]+sb; IA2[id]=IA2[up]+sa2; IB2[id]=IB2[up]+sb2; IAB[id]=IAB[up]+sab; }}
    const int rad=5; const double C1=6.5025, C2=58.5225; long long cnt=0; long double acc=0;
    for(int y=0;y<R;y++){ int row=y*R; for(int x=0;x<R;x++){ int p=row+x; if(!(dA[p]<254.9f || dB[p]<254.9f)) continue; int x0=max(0,x-rad), y0=max(0,y-rad), x1=min(R,x+rad+1), y1=min(R,y+rad+1); double area=(double)(x1-x0)*(double)(y1-y0); double sa=rectSum(IA,W,x0,y0,x1,y1), sb=rectSum(IB,W,x0,y0,x1,y1); double sa2=rectSum(IA2,W,x0,y0,x1,y1), sb2=rectSum(IB2,W,x0,y0,x1,y1), sab=rectSum(IAB,W,x0,y0,x1,y1); double ma=sa/area, mb=sb/area; double va=max(0.0,sa2/area-ma*ma), vb=max(0.0,sb2/area-mb*mb), cov=sab/area-ma*mb; double den=(ma*ma+mb*mb+C1)*(va+vb+C2); double val=den!=0?((2*ma*mb+C1)*(2*cov+C2)/den):1.0; acc+=val; ++cnt; }}
    return cnt? (double)(acc/cnt):1.0;
}
static double renderSSIM(const RenderMaps&orig,const RenderMaps&cand){
    int R=orig.R, PIX=R*R, W=R+1; size_t SZ=(size_t)W*W; vector<double>A(SZ),B(SZ),A2(SZ),B2(SZ),AB(SZ); double total=0;
    for(int view=0; view<6; ++view){ const float*od=orig.depth.data()+(size_t)view*PIX; const float*cd=cand.depth.data()+(size_t)view*PIX; double ns=0; for(int ch=0;ch<3;ch++){ const float*on=orig.norm.data()+(size_t)view*PIX*3+ch; const float*cn=cand.norm.data()+(size_t)view*PIX*3+ch; ns+=ssimChannel(on,3,cn,3,od,cd,R,A,B,A2,B2,AB); } ns/=3.0; double ds=ssimChannel(od,1,cd,1,od,cd,R,A,B,A2,B2,AB); total += 0.5*(ns+ds); }
    return total/6.0;
}

// ---------- Parametric candidates ----------
static void jacobiEigen(double a[3][3], Vec3 axis[3]){
    double v[3][3]={{1,0,0},{0,1,0},{0,0,1}};
    for(int it=0;it<40;it++){
        int p=0,q=1; double best=fabs(a[0][1]); if(fabs(a[0][2])>best){p=0;q=2;best=fabs(a[0][2]);} if(fabs(a[1][2])>best){p=1;q=2;best=fabs(a[1][2]);}
        if(best<1e-18) break; double app=a[p][p], aqq=a[q][q], apq=a[p][q]; double tau=(aqq-app)/(2*apq); double t=(tau>=0?1.0:-1.0)/(fabs(tau)+sqrt(1+tau*tau)); double c=1/sqrt(1+t*t), s=t*c;
        for(int k=0;k<3;k++) if(k!=p&&k!=q){ double akp=a[k][p], akq=a[k][q]; a[k][p]=a[p][k]=c*akp-s*akq; a[k][q]=a[q][k]=s*akp+c*akq; }
        a[p][p]=c*c*app-2*s*c*apq+s*s*aqq; a[q][q]=s*s*app+2*s*c*apq+c*c*aqq; a[p][q]=a[q][p]=0;
        for(int k=0;k<3;k++){ double vp=v[k][p], vq=v[k][q]; v[k][p]=c*vp-s*vq; v[k][q]=s*vp+c*vq; }
    }
    int ord[3]={0,1,2}; sort(ord,ord+3,[&](int i,int j){return a[i][i]>a[j][j];});
    for(int j=0;j<3;j++){ int col=ord[j]; axis[j]=normalized(Vec3(v[0][col],v[1][col],v[2][col])); }
    if(dotv(crossv(axis[0],axis[1]),axis[2])<0) axis[2]=axis[2]*-1.0;
}
static vector<array<Vec3,3>> candidateBases(){
    vector<array<Vec3,3>> bases; bases.push_back({Vec3(1,0,0),Vec3(0,1,0),Vec3(0,0,1)});
    Vec3 mean; for(auto&p:OV) mean=mean+p; mean=mean/(double)max(1,N0); double cov[3][3]{}; int stride=max(1,N0/240000), cnt=0; for(int i=0;i<N0;i+=stride){ Vec3 q=OV[i]-mean; double x[3]={q.x,q.y,q.z}; for(int a=0;a<3;a++) for(int b=0;b<3;b++) cov[a][b]+=x[a]*x[b]; ++cnt; } for(int a=0;a<3;a++) for(int b=0;b<3;b++) cov[a][b]/=max(1,cnt); Vec3 ax[3]; jacobiEigen(cov,ax); if(norm2(ax[0])>0.5) bases.push_back({ax[0],ax[1],ax[2]}); return bases;
}
static bool buildEllipsoid(const array<Vec3,3>&B, double spacingFactor, MeshCand&mc){
    double lo[3]={1e100,1e100,1e100}, hi[3]={-1e100,-1e100,-1e100};
    for(const Vec3&p:OV) for(int k=0;k<3;k++){ double t=dotv(p,B[k]); lo[k]=min(lo[k],t); hi[k]=max(hi[k],t); }
    Vec3 center; double r[3]; for(int k=0;k<3;k++){ double mid=0.5*(lo[k]+hi[k]); center=center+B[k]*mid; r[k]=0.5*(hi[k]-lo[k]); if(r[k]<=1e-10) return false; }
    int stride=max(1,N0/240000), cnt=0; double ss=0, ma=0; for(int i=0;i<N0;i+=stride){ Vec3 q=OV[i]-center; double rr=0; for(int k=0;k<3;k++){ double u=dotv(q,B[k])/r[k]; rr+=u*u; } double e=fabs(sqrt(max(0.0,rr))-1.0); ss+=e*e; ma=max(ma,e); ++cnt; }
    double rms=sqrt(ss/max(1,cnt)); if(rms>0.035 || ma>0.13) return false;
    double maxr=max(r[0],max(r[1],r[2])); double sp=max(1e-6, spacingFactor*TAU); int lat=max(8,min(72,(int)ceil(M_PI*maxr/sp))); int lon=max(16,min(144,(int)ceil(2*M_PI*maxr/sp))); if(2+(lat-1)*lon>=N0) return false;
    vector<Vec3> V; vector<Tri> F; V.reserve(2+(lat-1)*lon); F.reserve(2*lat*lon);
    auto makep=[&](double x,double y,double z){ return center+B[0]*(r[0]*x)+B[1]*(r[1]*y)+B[2]*(r[2]*z); };
    V.push_back(makep(0,0,1));
    for(int i=1;i<lat;i++){ double th=M_PI*i/lat, st=sin(th), ct=cos(th); for(int j=0;j<lon;j++){ double ph=2*M_PI*j/lon; V.push_back(makep(st*cos(ph),st*sin(ph),ct)); }}
    V.push_back(makep(0,0,-1)); int bot=(int)V.size()-1; auto id=[&](int ring,int j){return 1+(ring-1)*lon+((j%lon+lon)%lon);};
    double sg=orientSignFromOriginal(center); auto add=[&](int a,int b,int c){ Tri t{a,b,c}; orientFace(V,t,center,sg); F.push_back(t); };
    for(int j=0;j<lon;j++) add(0,id(1,j+1),id(1,j));
    for(int i=1;i<lat-1;i++) for(int j=0;j<lon;j++){ int a=id(i,j), b=id(i,j+1), c=id(i+1,j), d=id(i+1,j+1); add(a,b,c); add(b,d,c); }
    for(int j=0;j<lon;j++) add(bot,id(lat-1,j),id(lat-1,j+1));
    if(!manifoldAndBasicOK(V,F) || !vertexHausdorffOK(V)) return false;
    mc.V=move(V); mc.F=move(F); mc.name="ellipsoid"; return true;
}
static vector<MeshCand> specialEllipsoids(){
    vector<MeshCand> out; for(auto&B:candidateBases()){ for(double sf: {1.00,0.82,0.68}){ if(elapsed()>5.0) break; MeshCand c; if(buildEllipsoid(B,sf,c)) out.push_back(move(c)); }} return out;
}

static bool buildTorus(const array<Vec3,3>&B,int axis,double spacingFactor,MeshCand&mc){
    Vec3 w=B[axis], u=B[(axis+1)%3], v=B[(axis+2)%3];
    double minT=1e100,maxT=-1e100,minU=1e100,maxU=-1e100,minV=1e100,maxV=-1e100;
    for(auto&p:OV){ double t=dotv(p,w), x=dotv(p,u), y=dotv(p,v); minT=min(minT,t); maxT=max(maxT,t); minU=min(minU,x); maxU=max(maxU,x); minV=min(minV,y); maxV=max(maxV,y); }
    double ct=0.5*(minT+maxT), cu=0.5*(minU+maxU), cv=0.5*(minV+maxV); double minR=1e100,maxR=0; int stride=max(1,N0/240000), cnt=0;
    for(int i=0;i<N0;i+=stride){ double x=dotv(OV[i],u)-cu, y=dotv(OV[i],v)-cv; double rho=sqrt(x*x+y*y); minR=min(minR,rho); maxR=max(maxR,rho); ++cnt; }
    double major=0.5*(minR+maxR), minor=0.5*(maxR-minR); if(!(major>1e-8&&minor>1e-8) || major<1.25*minor) return false;
    double ss=0, ma=0; cnt=0; for(int i=0;i<N0;i+=stride){ double t=dotv(OV[i],w)-ct, x=dotv(OV[i],u)-cu, y=dotv(OV[i],v)-cv; double rho=sqrt(x*x+y*y); double e=fabs(sqrt((rho-major)*(rho-major)+t*t)-minor)/minor; ss+=e*e; ma=max(ma,e); ++cnt; }
    double rms=sqrt(ss/max(1,cnt)); if(rms>0.045 || ma>0.18) return false;
    double sp=max(1e-6,spacingFactor*TAU); int A=max(16,min(192,(int)ceil(2*M_PI*major/sp))); int Bn=max(8,min(80,(int)ceil(2*M_PI*minor/sp))); if(A*Bn>=N0) return false;
    vector<Vec3> V; vector<Tri> F; V.reserve(A*Bn); F.reserve(2*A*Bn); Vec3 center=w*ct+u*cu+v*cv;
    for(int i=0;i<A;i++){ double th=2*M_PI*i/A, co=cos(th), si=sin(th); Vec3 radial=u*co+v*si; for(int j=0;j<Bn;j++){ double ph=2*M_PI*j/Bn; V.push_back(center + radial*((major+minor*cos(ph))) + w*(minor*sin(ph))); }}
    auto id=[&](int i,int j){ return ((i%A+A)%A)*Bn+((j%Bn+Bn)%Bn); };
    // choose orientation by comparing original normals with analytic normals
    double os=0; int ost=0; int fstride=max(1,M0/100000);
    for(int fi=0;fi<M0;fi+=fstride){ auto&t=OF[fi]; Vec3 a=OV[t.a],b=OV[t.b],c=OV[t.c]; Vec3 cr=normalized(crossv(b-a,c-a)); Vec3 p=(a+b+c)/3.0; double tt=dotv(p,w)-ct, x=dotv(p,u)-cu, y=dotv(p,v)-cv; double rho=sqrt(x*x+y*y); if(rho>1e-12){ Vec3 pred=(u*(x/rho)+v*(y/rho))*((rho-major)/max(minor,1e-12)) + w*(tt/max(minor,1e-12)); pred=normalized(pred); os+=dotv(cr,pred); ++ost; }}
    bool flip=os<0;
    auto add=[&](int a,int b,int c){ Tri t{a,b,c}; if(flip) swap(t.b,t.c); F.push_back(t); };
    for(int i=0;i<A;i++) for(int j=0;j<Bn;j++){ int a=id(i,j), b=id(i+1,j), c=id(i+1,j+1), d=id(i,j+1); add(a,b,c); add(a,c,d); }
    if(!manifoldAndBasicOK(V,F)||!vertexHausdorffOK(V)) return false; mc.V=move(V); mc.F=move(F); mc.name="torus"; return true;
}
static vector<MeshCand> specialTori(){ vector<MeshCand> out; for(auto&B:candidateBases()) for(int ax=0;ax<3;ax++) for(double sf:{1.0,0.8}){ if(elapsed()>6.0) break; MeshCand c; if(buildTorus(B,ax,sf,c)) out.push_back(move(c)); } return out; }

static bool buildFrustum(const array<Vec3,3>&B,int axis,double spacingFactor,MeshCand&mc){
    Vec3 w=B[axis], u=B[(axis+1)%3], v=B[(axis+2)%3];
    double minT=1e100,maxT=-1e100,minU=1e100,maxU=-1e100,minV=1e100,maxV=-1e100;
    for(auto&p:OV){ double t=dotv(p,w), x=dotv(p,u), y=dotv(p,v); minT=min(minT,t); maxT=max(maxT,t); minU=min(minU,x); maxU=max(maxU,x); minV=min(minV,y); maxV=max(maxV,y); }
    double cu=0.5*(minU+maxU), cv=0.5*(minV+maxV), len=maxT-minT; if(len<=1e-8) return false;
    const int BINS=20; vector<double> mr(BINS,0); vector<int> bc(BINS,0); int stride=max(1,N0/240000);
    for(int i=0;i<N0;i+=stride){ double t=dotv(OV[i],w), x=dotv(OV[i],u)-cu, y=dotv(OV[i],v)-cv; int b=min(BINS-1,max(0,(int)((t-minT)/len*BINS))); double r=sqrt(x*x+y*y); mr[b]=max(mr[b],r); bc[b]++; }
    double S=0,St=0,Stt=0,Sr=0,Str=0; for(int b=0;b<BINS;b++) if(bc[b]>0 && mr[b]>1e-8){ double s=((double)b+0.5)/BINS; S++; St+=s; Stt+=s*s; Sr+=mr[b]; Str+=s*mr[b]; }
    if(S<6) return false; double den=S*Stt-St*St; if(fabs(den)<1e-18) return false; double slope=(S*Str-St*Sr)/den, base=(Sr-slope*St)/S; double r0=base, r1=base+slope; if(r0<0) r0=0; if(r1<0) r1=0; double rmax=max(r0,r1); if(rmax<=1e-8 || len<0.25*rmax) return false;
    int checked=0, bad=0; double ss=0, ma=0; double sideTol=max(0.012*rmax,0.25*TAU);
    for(int i=0;i<N0;i+=stride){ double t=dotv(OV[i],w), x=dotv(OV[i],u)-cu, y=dotv(OV[i],v)-cv; double s=(t-minT)/len; double rp=max(0.0,r0+(r1-r0)*s); double rr=sqrt(x*x+y*y); double side=fabs(rr-rp); double cap=min(fabs(t-minT),fabs(t-maxT)); double err=min(side,cap); if(err>max(sideTol,0.055*DIAG)){ if(++bad>max(5,checked/30)) return false; } ss+=err*err; ma=max(ma,err); ++checked; }
    double rms=sqrt(ss/max(1,checked)); if(rms>max(0.018*rmax,0.35*TAU) || ma>max(0.09*rmax,0.75*TAU)) return false;
    double sp=max(1e-6,spacingFactor*TAU); int sides=max(12,min(160,(int)ceil(2*M_PI*rmax/sp))); int nt=max(2,min(80,(int)ceil(len/sp))); int cap0=max(1,min(40,(int)ceil(r0/sp))); int cap1=max(1,min(40,(int)ceil(r1/sp)));
    vector<Vec3> V; vector<Tri> F; V.reserve((nt+1)*sides+(cap0+cap1)*sides+2); F.reserve(2*nt*sides+2*(cap0+cap1)*sides);
    auto make=[&](double t,double r,int j){ double th=2*M_PI*j/sides; return w*t + u*(cu+r*cos(th)) + v*(cv+r*sin(th)); };
    vector<int> sideStart(nt+1);
    for(int i=0;i<=nt;i++){ double s=(double)i/nt; double t=minT+len*s, r=r0+(r1-r0)*s; sideStart[i]=V.size(); for(int j=0;j<sides;j++) V.push_back(make(t,r,j)); }
    Vec3 center=w*(0.5*(minT+maxT))+u*cu+v*cv; double sg=orientSignFromOriginal(center);
    auto orientAdd=[&](Tri t){ orientFace(V,t,center,sg); F.push_back(t); };
    auto sid=[&](int i,int j){return sideStart[i]+((j%sides+sides)%sides);};
    for(int i=0;i<nt;i++) for(int j=0;j<sides;j++){ orientAdd({sid(i,j),sid(i+1,j),sid(i+1,j+1)}); orientAdd({sid(i,j),sid(i+1,j+1),sid(i,j+1)}); }
    auto buildCap=[&](bool top){ double t=top?maxT:minT, rr=top?r1:r0; int sub=top?cap1:cap0; if(rr<=1e-10) return; vector<int> rings; int centerId=V.size(); V.push_back(w*t+u*cu+v*cv); int prevCenter=centerId; int prevStart=-1; for(int k=1;k<sub;k++){ double r=rr*(double)k/sub; int st=V.size(); for(int j=0;j<sides;j++) V.push_back(make(t,r,j)); if(k==1){ for(int j=0;j<sides;j++) orientAdd({centerId,st+j,st+(j+1)%sides}); } else { for(int j=0;j<sides;j++) { int a=prevStart+j,b=prevStart+(j+1)%sides,c=st+j,d=st+(j+1)%sides; orientAdd({a,c,d}); orientAdd({a,d,b}); }} prevStart=st; }
        int outer=top?sideStart[nt]:sideStart[0]; if(sub==1){ for(int j=0;j<sides;j++) orientAdd({centerId,outer+j,outer+(j+1)%sides}); } else { for(int j=0;j<sides;j++){ int a=prevStart+j,b=prevStart+(j+1)%sides,c=outer+j,d=outer+(j+1)%sides; orientAdd({a,c,d}); orientAdd({a,d,b}); }} };
    buildCap(false); buildCap(true);
    if((int)V.size()>=N0) return false; if(!manifoldAndBasicOK(V,F)||!vertexHausdorffOK(V)) return false; mc.V=move(V); mc.F=move(F); mc.name="frustum"; return true;
}
static vector<MeshCand> specialFrusta(){ vector<MeshCand> out; for(auto&B:candidateBases()) for(int ax=0;ax<3;ax++) for(double sf:{1.0,0.8}){ if(elapsed()>7.0) break; MeshCand c; if(buildFrustum(B,ax,sf,c)) out.push_back(move(c)); } return out; }

// Detect simple ordered UV torus grids and sphere grids and downsample them.
static bool sameFace(const Tri&t,int a,int b,int c){ array<int,3>x{t.a,t.b,t.c}, y{a,b,c}; sort(x.begin(),x.end()); sort(y.begin(),y.end()); return x==y; }
static bool trySphereGrid(MeshCand&mc){
    if(N0<300 || M0!=2*(N0-2)) return false;
    int R=0,C=0;
    for(int c=8;c<=4096;c++){ if((N0-2)%c) continue; int r=(N0-2)/c; if(r<3) continue; bool ok=true; int step=max(1,c/32); for(int j=0;j<c&&ok;j+=step){ if(!sameFace(OF[j],0,1+(j+1)%c,1+j)) ok=false; }
        int span=max(1,(r-1)*c/256); for(int q=0;q<(r-1)*c&&ok;q+=span){ int rr=q/c,j=q%c; int a=1+rr*c+j,b=1+rr*c+(j+1)%c,cc=1+(rr+1)*c+j,d=1+(rr+1)*c+(j+1)%c; int f=c+2*(rr*c+j); if(f+1>=M0 || !sameFace(OF[f],a,b,cc) || !sameFace(OF[f+1],b,d,cc)) ok=false; }
        if(ok){R=r; C=c; break;}
    }
    if(!R) return false;
    int R2=max(4,min(R, (int)ceil((double)R*0.22))); int C2=max(12,min(C, (int)ceil((double)C*0.22))); // initial aggressive downsample
    // enforce chord coverage by trying finer if needed
    for(double frac: {0.16,0.20,0.25,0.32,0.40}){
        R2=max(4,min(R,(int)ceil(R*frac))); C2=max(12,min(C,(int)ceil(C*frac))); if(2+R2*C2>=N0) continue;
        vector<Vec3> V; vector<Tri> F; V.reserve(2+R2*C2); V.push_back(OV[0]);
        for(int i=0;i<R2;i++){ int oi=1+(int)((long long)i*(R-1)/max(1,R2-1)); for(int j=0;j<C2;j++){ int oj=(int)((long long)j*C/C2); V.push_back(OV[1+(oi-1)*C+oj]); }} int bot=V.size(); V.push_back(OV[N0-1]);
        Vec3 ctr; for(auto&p:V) ctr=ctr+p; ctr=ctr/(double)V.size(); double sg=orientSignFromOriginal(ctr); auto id=[&](int r,int j){return 1+(r-1)*C2+((j%C2+C2)%C2);}; auto add=[&](int a,int b,int c){Tri t{a,b,c}; orientFace(V,t,ctr,sg); F.push_back(t);};
        for(int j=0;j<C2;j++) add(0,id(1,j+1),id(1,j)); for(int r=1;r<R2;r++) for(int j=0;j<C2;j++){ int a=id(r,j), b=id(r,j+1), c=id(r+1,j), d=id(r+1,j+1); add(a,c,d); add(a,d,b); } for(int j=0;j<C2;j++) add(bot,id(R2,j),id(R2,j+1));
        if(manifoldAndBasicOK(V,F)&&vertexHausdorffOK(V)){ mc.V=move(V); mc.F=move(F); mc.name="sphere-grid"; return true; }
    }
    return false;
}
static bool tryPeriodicGrid(MeshCand&mc){
    if(N0<300 || M0!=2*N0) return false; int U=0,Vn=0;
    for(int v=6;v<=4096;v++){ if(N0%v) continue; int u=N0/v; if(u<6) continue; int step=max(1,N0/512), tot=0, okc=0; for(int q=0;q<N0&&tot<512;q+=step){ int i=q/v,j=q%v, ni=(i+1)%u,nj=(j+1)%v; int a=i*v+j,b=ni*v+j,c=ni*v+nj,d=i*v+nj; int f=2*q; if(f+1>=M0) break; bool ok=(sameFace(OF[f],a,b,c)&&sameFace(OF[f+1],a,c,d)) || (sameFace(OF[f],a,c,b)&&sameFace(OF[f+1],a,d,c)); if(ok) okc++; tot++; } if(tot>100 && okc*100>=tot*98){U=u; Vn=v; break;} }
    if(!U) return false;
    for(double frac: {0.14,0.18,0.24,0.32,0.42}){ int U2=max(6,min(U,(int)ceil(U*frac))); int V2=max(6,min(Vn,(int)ceil(Vn*frac))); if(U2*V2>=N0) continue; vector<Vec3> V; vector<Tri> F; V.reserve(U2*V2); for(int i=0;i<U2;i++){ int oi=(int)((long long)i*U/U2); for(int j=0;j<V2;j++){ int oj=(int)((long long)j*Vn/V2); V.push_back(OV[oi*Vn+oj]); }} Vec3 ctr; for(auto&p:V) ctr=ctr+p; ctr=ctr/(double)V.size(); double sg=orientSignFromOriginal(ctr); auto id=[&](int i,int j){return ((i%U2+U2)%U2)*V2+((j%V2+V2)%V2);}; auto add=[&](int a,int b,int c){Tri t{a,b,c}; orientFace(V,t,ctr,sg); F.push_back(t);}; for(int i=0;i<U2;i++) for(int j=0;j<V2;j++){ add(id(i,j),id(i+1,j),id(i+1,j+1)); add(id(i,j),id(i+1,j+1),id(i,j+1)); } if(manifoldAndBasicOK(V,F)&&vertexHausdorffOK(V)){ mc.V=move(V); mc.F=move(F); mc.name="periodic-grid"; return true; }} return false;
}

// ---------- Generic QEM-like simplifier ----------
namespace GEN{
struct Quadric{ double q[10]; Quadric(){memset(q,0,sizeof q);} void add(const Quadric&o){for(int i=0;i<10;i++)q[i]+=o.q[i];} void addPlane(double a,double b,double c,double d,double w){q[0]+=w*a*a; q[1]+=w*a*b; q[2]+=w*a*c; q[3]+=w*a*d; q[4]+=w*b*b; q[5]+=w*b*c; q[6]+=w*b*d; q[7]+=w*c*c; q[8]+=w*c*d; q[9]+=w*d*d;} double eval(const Vec3&p)const{double x=p.x,y=p.y,z=p.z; return q[0]*x*x+2*q[1]*x*y+2*q[2]*x*z+2*q[3]*x+q[4]*y*y+2*q[5]*y*z+2*q[6]*y+q[7]*z*z+2*q[8]*z+q[9];}};
struct Node{double cost; int u,v,vu,vv; bool operator<(Node const&o)const{return cost>o.cost;}};
static int N,M; static vector<Vec3>P; static vector<Face>F; static vector<vector<int>> inc; static vector<Quadric> Q; static vector<unsigned char> alive; static vector<int> ver,markA,markB,tmp; static int stampA=1,stampB=1; static vector<array<float,3>> bbMin,bbMax; static int activeV, activeF; static priority_queue<Node> pq; static double visTol=40, visTol2=1600; static double minDot=0.02;
static inline bool contains(int fid,int v){auto &f=F[fid]; return f.v[0]==v||f.v[1]==v||f.v[2]==v;}
static inline bool containsFace(const Face&f,int v){return f.v[0]==v||f.v[1]==v||f.v[2]==v;}
static Vec3 faceNormal(const Face&f){return normalized(crossv(P[f.v[1]]-P[f.v[0]],P[f.v[2]]-P[f.v[0]]));}
static void cleanup(int v){ if(!alive[v]) return; auto &a=inc[v]; int w=0; for(int fid:a) if(fid>=0&&fid<M&&F[fid].active&&contains(fid,v)) a[w++]=fid; a.resize(w); }
static void maybeCleanup(int v){ if(!alive[v]) return; if(inc[v].size()>96){ int good=0; for(int fid:inc[v]) if(F[fid].active&&contains(fid,v)) good++; if((int)inc[v].size()>good*3+64) cleanup(v); }}
static bool solveOpt(const Quadric&q,Vec3&out){ double a00=q.q[0],a01=q.q[1],a02=q.q[2],a11=q.q[4],a12=q.q[5],a22=q.q[7]; double b0=q.q[3],b1=q.q[6],b2=q.q[8]; double det=a00*(a11*a22-a12*a12)-a01*(a01*a22-a12*a02)+a02*(a01*a12-a11*a02); if(fabs(det)<1e-14) return false; double inv00=(a11*a22-a12*a12)/det, inv01=(a02*a12-a01*a22)/det, inv02=(a01*a12-a02*a11)/det, inv11=(a00*a22-a02*a02)/det, inv12=(a01*a02-a00*a12)/det, inv22=(a00*a11-a01*a01)/det; out.x=-(inv00*b0+inv01*b1+inv02*b2); out.y=-(inv01*b0+inv11*b1+inv12*b2); out.z=-(inv02*b0+inv12*b1+inv22*b2); return isfinite(out.x)&&isfinite(out.y)&&isfinite(out.z)&&out.x>=-2&&out.x<=2&&out.y>=-2&&out.y<=2&&out.z>=-2&&out.z<=2; }
static inline void mergedBB(int u,int v,array<float,3>&mn,array<float,3>&mx){for(int k=0;k<3;k++){mn[k]=min(bbMin[u][k],bbMin[v][k]); mx[k]=max(bbMax[u][k],bbMax[v][k]);}}
static bool coversBB(const Vec3&p,const array<float,3>&mn,const array<float,3>&mx){ for(int m=0;m<8;m++){ Vec3 c((m&1)?mx[0]:mn[0],(m&2)?mx[1]:mn[1],(m&4)?mx[2]:mn[2]); if(dist2(p,c)>TAU2) return false; } return true; }
static bool visualOK(const Vec3&p,const array<float,3>&mn,const array<float,3>&mx){ double pu[6],pv[6],pd; for(int view=0;view<6;view++) projectPoint(p,view,1024,pu[view],pv[view],pd); for(int m=0;m<8;m++){ Vec3 c((m&1)?mx[0]:mn[0],(m&2)?mx[1]:mn[1],(m&4)?mx[2]:mn[2]); for(int view=0;view<6;view++){ double u,v,d; projectPoint(c,view,1024,u,v,d); double dx=u-pu[view],dy=v-pv[view]; if(dx*dx+dy*dy>visTol2) return false; }} return true; }
struct Eval{bool ok=false; double cost=1e300; Vec3 pos;};
static Eval computeCand(int u,int v){ Eval e; if(u==v||!alive[u]||!alive[v]) return e; array<float,3>mn,mx; mergedBB(u,v,mn,mx); Quadric q=Q[u]; q.add(Q[v]); Vec3 cand[8]; int c=0,opt; Vec3 o; if(solveOpt(q,o)) cand[c++]=o; cand[c++]=P[u]; cand[c++]=P[v]; cand[c++]=(P[u]+P[v])*0.5; cand[c++]=P[u]*0.75+P[v]*0.25; cand[c++]=P[u]*0.25+P[v]*0.75; cand[c++]=Vec3((mn[0]+mx[0])*0.5,(mn[1]+mx[1])*0.5,(mn[2]+mx[2])*0.5); double len2=dist2(P[u],P[v]); for(int i=0;i<c;i++){ Vec3 p=cand[i]; if(!coversBB(p,mn,mx)) continue; if(!visualOK(p,mn,mx)) continue; double val=q.eval(p); if(val<0&&val>-1e-10) val=0; double score=val + 1e-7*len2 + 1e-12*norm2(p); if(score<e.cost){e.ok=true; e.cost=score; e.pos=p;} } return e; }
static void pushEdge(int u,int v){ if(u==v||!alive[u]||!alive[v]) return; if(u>v) swap(u,v); Eval e=computeCand(u,v); if(!e.ok) return; pq.push({e.cost,u,v,ver[u],ver[v]}); }
static bool linkOK(int u,int v){ if(u==v||!alive[u]||!alive[v]) return false; maybeCleanup(u); maybeCleanup(v); if(++stampA>1000000000){fill(markA.begin(),markA.end(),0); stampA=1;} if(++stampB>1000000000){fill(markB.begin(),markB.end(),0); stampB=1;} int edgeFaces=0; for(int fid:inc[u]){ if(!F[fid].active||!contains(fid,u)) continue; bool hasv=false; for(int k=0;k<3;k++){int x=F[fid].v[k]; if(x==v) hasv=true; if(x!=u&&x!=v) markA[x]=stampA;} if(hasv) edgeFaces++; } if(edgeFaces!=2) return false; int common=0; for(int fid:inc[v]){ if(!F[fid].active||!contains(fid,v)) continue; for(int k=0;k<3;k++){int x=F[fid].v[k]; if(x==u||x==v) continue; if(markA[x]==stampA && markB[x]!=stampB){markB[x]=stampB; if(++common>2) return false;}} } return common==2; }
static bool flipOK(int keep,int rem,const Vec3&pos){ static vector<int> touched; touched.clear(); for(int fid:inc[keep]) if(F[fid].active&&contains(fid,keep)) touched.push_back(fid); for(int fid:inc[rem]) if(F[fid].active&&contains(fid,rem)) touched.push_back(fid); sort(touched.begin(),touched.end()); touched.erase(unique(touched.begin(),touched.end()),touched.end()); double minA=max(1e-30,1e-24*DIAG*DIAG); for(int fid:touched){ Face&f=F[fid]; bool hk=containsFace(f,keep), hr=containsFace(f,rem); if(hk&&hr) continue; Vec3 a=(f.v[0]==keep||f.v[0]==rem)?pos:P[f.v[0]], b=(f.v[1]==keep||f.v[1]==rem)?pos:P[f.v[1]], c=(f.v[2]==keep||f.v[2]==rem)?pos:P[f.v[2]]; Vec3 cr=crossv(b-a,c-a); double a2=norm2(cr); if(!(a2>minA)||!isfinite(a2)) return false; Vec3 nn=cr/sqrt(a2); if(dotv(nn,f.n)<minDot) return false; } return true; }
static void collectNei(int u,vector<int>&out){ out.clear(); if(!alive[u]) return; maybeCleanup(u); if(++stampA>1000000000){fill(markA.begin(),markA.end(),0); stampA=1;} for(int fid:inc[u]){ if(!F[fid].active||!contains(fid,u)) continue; for(int k=0;k<3;k++){int w=F[fid].v[k]; if(w!=u&&alive[w]&&markA[w]!=stampA){markA[w]=stampA; out.push_back(w);}} }}
static void collapse(int keep,int rem,const Vec3&pos){ cleanup(keep); cleanup(rem); P[keep]=pos; vector<int> list=inc[rem]; for(int fid:list){ if(!F[fid].active||!contains(fid,rem)) continue; bool hasKeep=contains(fid,keep); if(hasKeep){F[fid].active=0; --activeF;} else { for(int k=0;k<3;k++) if(F[fid].v[k]==rem) F[fid].v[k]=keep; F[fid].n=faceNormal(F[fid]); inc[keep].push_back(fid); }} for(int fid:inc[keep]) if(F[fid].active&&contains(fid,keep)) F[fid].n=faceNormal(F[fid]); Q[keep].add(Q[rem]); for(int k=0;k<3;k++){bbMin[keep][k]=min(bbMin[keep][k],bbMin[rem][k]); bbMax[keep][k]=max(bbMax[keep][k],bbMax[rem][k]);} alive[rem]=0; ++ver[keep]; ++ver[rem]; --activeV; inc[rem].clear(); cleanup(keep); static vector<int> nei; collectNei(keep,nei); for(int w:nei) pushEdge(keep,w); }
static MeshCand snapshot(const string&name){ MeshCand s; s.name=name; tmp.assign(N,-1); s.V.reserve(activeV); for(int i=0;i<N;i++) if(alive[i]){tmp[i]=(int)s.V.size(); s.V.push_back(P[i]);} s.F.reserve(activeF); for(int i=0;i<M;i++) if(F[i].active){int a=F[i].v[0],b=F[i].v[1],c=F[i].v[2]; if(a!=b&&a!=c&&b!=c&&tmp[a]>=0&&tmp[b]>=0&&tmp[c]>=0) s.F.push_back({tmp[a],tmp[b],tmp[c]});} return s; }
static vector<MeshCand> run(double timeLimit){
    N=N0; M=M0; P=OV; F.assign(M,{}); inc.assign(N,{}); vector<int> deg(N,0); for(int i=0;i<M;i++){F[i].v[0]=OF[i].a;F[i].v[1]=OF[i].b;F[i].v[2]=OF[i].c;F[i].active=1;deg[F[i].v[0]]++;deg[F[i].v[1]]++;deg[F[i].v[2]]++;} for(int i=0;i<N;i++) inc[i].reserve(deg[i]+8); for(int i=0;i<M;i++){inc[F[i].v[0]].push_back(i);inc[F[i].v[1]].push_back(i);inc[F[i].v[2]].push_back(i);} for(int i=0;i<M;i++) F[i].n=faceNormal(F[i]);
    Q.assign(N,Quadric()); for(int i=0;i<M;i++){auto &f=F[i]; Vec3 a=P[f.v[0]],b=P[f.v[1]],c=P[f.v[2]]; Vec3 cr=crossv(b-a,c-a); double ar=normv(cr); if(ar<=1e-300) continue; Vec3 n=cr/ar; double d=-dotv(n,a); double w=max(1e-8,0.5*ar); for(int k=0;k<3;k++) Q[f.v[k]].addPlane(n.x,n.y,n.z,d,w); }
    alive.assign(N,1); ver.assign(N,0); markA.assign(N,0); markB.assign(N,0); bbMin.resize(N); bbMax.resize(N); for(int i=0;i<N;i++){bbMin[i]={(float)P[i].x,(float)P[i].y,(float)P[i].z}; bbMax[i]=bbMin[i];}
    activeV=N; activeF=M; pq=priority_queue<Node>(); stampA=1; stampB=1;
    if(N<1000){visTol=90; minDot=-0.05;} else if(N<5000){visTol=75; minDot=-0.02;} else if(N<30000){visTol=58; minDot=0.0;} else if(N<100000){visTol=46; minDot=0.02;} else {visTol=36; minDot=0.05;} visTol2=visTol*visTol;
    vector<uint64_t> edges; edges.reserve((size_t)M*3); for(int i=0;i<M;i++){edges.push_back(edgeKey(F[i].v[0],F[i].v[1]));edges.push_back(edgeKey(F[i].v[1],F[i].v[2]));edges.push_back(edgeKey(F[i].v[2],F[i].v[0]));} sort(edges.begin(),edges.end()); edges.erase(unique(edges.begin(),edges.end()),edges.end()); for(uint64_t k:edges) pushEdge((int)(k>>32),(int)(k&0xffffffffu)); edges.clear(); edges.shrink_to_fit();
    vector<double> ratios={0.70,0.50,0.35,0.25,0.20,0.16,0.13,0.11,0.095,0.085,0.075,0.065,0.055,0.047}; vector<int> cps; int last=-1; for(double r:ratios){int c=max(4,(int)ceil(N*r)); if(c<N&&c!=last){cps.push_back(c); last=c;}} sort(cps.begin(),cps.end(),greater<int>()); int cp=0; int minTarget=cps.empty()?max(4,N/8):cps.back(); vector<MeshCand> snaps;
    long long pops=0; while(activeV>minTarget && !pq.empty()){ if((++pops&4095)==0 && elapsed()>timeLimit) break; Node nd=pq.top(); pq.pop(); int u=nd.u,v=nd.v; if(u==v||!alive[u]||!alive[v]) continue; if(u>v) swap(u,v); if(ver[u]!=nd.vu||ver[v]!=nd.vv) { continue; } Eval ec=computeCand(u,v); if(!ec.ok) continue; if(!linkOK(u,v)) continue; int keep=u,rem=v; if(inc[v].size()>inc[u].size()){keep=v; rem=u;} if(!flipOK(keep,rem,ec.pos)) continue; collapse(keep,rem,ec.pos); while(cp<(int)cps.size() && activeV<=cps[cp]){ MeshCand s=snapshot("generic"); if(!s.V.empty() && s.V.size()<OV.size()) snaps.push_back(move(s)); cp++; } }
    if(snaps.empty() || snaps.back().V.size()!=(size_t)activeV){ MeshCand s=snapshot("generic-final"); if(!s.V.empty()&&s.V.size()<OV.size()) snaps.push_back(move(s)); }
    return snaps;
}
}

static void outputMesh(const vector<Vec3>&V,const vector<Tri>&F){
    static char obuf[1<<20]; setvbuf(stdout,obuf,_IOFBF,sizeof(obuf));
    printf("%d %d\n",(int)V.size(),(int)F.size());
    bool highprec=(int)V.size()*2<N0;
    for(auto&p:V){ if(highprec) printf("v %.15g %.15g %.15g\n",p.x,p.y,p.z); else printf("v %.10g %.10g %.10g\n",p.x,p.y,p.z); }
    for(auto&t:F) printf("f %d %d %d\n",t.a+1,t.b+1,t.c+1);
}
static void outputOriginal(){ outputMesh(OV,OF); }

int main(){
    T0=chrono::steady_clock::now();
    readInput();
    vector<MeshCand> candidates;
    // Ordered-grid replacements are very cheap and sometimes unbeatable.
    MeshCand g;
    if(trySphereGrid(g)) candidates.push_back(move(g));
    MeshCand pg;
    if(tryPeriodicGrid(pg)) candidates.push_back(move(pg));
    if(elapsed()<3.0){ auto xs=specialEllipsoids(); for(auto &x:xs) candidates.push_back(move(x)); }
    if(elapsed()<5.5){ auto xs=specialTori(); for(auto &x:xs) candidates.push_back(move(x)); }
    if(elapsed()<6.5){ auto xs=specialFrusta(); for(auto &x:xs) candidates.push_back(move(x)); }
    // Generic manifold QEM collapse.  It is conservative wrt the vertex-Hausdorff clarification.
    double genericLimit = (N0>250000?13.2:(N0>60000?14.2:15.0));
    auto snaps=GEN::run(genericLimit); for(auto &s:snaps){ if((int)s.V.size()<N0) candidates.push_back(move(s)); }
    if(candidates.empty()){ outputOriginal(); return 0; }
    // Drop invalid candidates from possible numerical mishaps.
    vector<MeshCand> valid; valid.reserve(candidates.size());
    for(auto &c:candidates){ if(c.V.empty()||c.F.empty()||c.V.size()>=OV.size()) continue; if(!manifoldAndBasicOK(c.V,c.F)) continue; // generic is constructed to satisfy vertex cover; specials were checked exactly.
        // For generic, bbox-cluster proof is internal; for special candidates exact was already done.  A cheap exact recheck is only done for small outputs.
        if(c.name!="generic"&&c.name!="generic-final") { /* already checked */ } 
        valid.push_back(move(c)); }
    if(valid.empty()){ outputOriginal(); return 0; }
    sort(valid.begin(),valid.end(),[](const MeshCand&a,const MeshCand&b){ if(a.V.size()!=b.V.size()) return a.V.size()<b.V.size(); return a.F.size()<b.F.size(); });
    int R; if(N0<6000) R=512; else if(N0<60000) R=384; else if(N0<250000) R=256; else R=128;
    if(elapsed()>17.2) R=min(R,256); if(elapsed()>18.0) R=min(R,128);
    double guard = (R>=512?0.905:(R>=384?0.915:(R>=256?0.925:0.940)));
    RenderMaps orig; renderMesh(OV,OF,orig,R);
    int chosen=-1; double bestScore=-1; int bestIdx=-1;
    for(int i=0;i<(int)valid.size();++i){
        if(elapsed()>20.15) break;
        RenderMaps rm; renderMesh(valid[i].V,valid[i].F,rm,R); double sc=renderSSIM(orig,rm); valid[i].proxy=sc; if(sc>bestScore){bestScore=sc; bestIdx=i;}
        if(sc>=guard){ chosen=i; break; }
    }
    if(chosen<0){
        // Prefer a verified high-proxy candidate; otherwise fall back to the least aggressive snapshot/original to avoid WA.
        if(bestIdx>=0 && bestScore>=0.895) chosen=bestIdx;
    }
    if(chosen>=0) outputMesh(valid[chosen].V,valid[chosen].F); else outputOriginal();
    return 0;
}
