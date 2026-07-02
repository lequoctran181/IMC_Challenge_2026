#include <bits/stdc++.h>
using namespace std;

struct Vec3 {
    double x, y, z;
    Vec3() : x(0), y(0), z(0) {}
    Vec3(double X, double Y, double Z) : x(X), y(Y), z(Z) {}
};
static inline Vec3 operator+(const Vec3& a, const Vec3& b){ return Vec3(a.x+b.x,a.y+b.y,a.z+b.z); }
static inline Vec3 operator-(const Vec3& a, const Vec3& b){ return Vec3(a.x-b.x,a.y-b.y,a.z-b.z); }
static inline Vec3 operator*(const Vec3& a, double s){ return Vec3(a.x*s,a.y*s,a.z*s); }
static inline Vec3 operator/(const Vec3& a, double s){ return Vec3(a.x/s,a.y/s,a.z/s); }
static inline double dotv(const Vec3& a, const Vec3& b){ return a.x*b.x+a.y*b.y+a.z*b.z; }
static inline Vec3 crossv(const Vec3& a, const Vec3& b){ return Vec3(a.y*b.z-a.z*b.y, a.z*b.x-a.x*b.z, a.x*b.y-a.y*b.x); }
static inline double norm2(const Vec3& a){ return dotv(a,a); }
static inline double normv(const Vec3& a){ return sqrt(norm2(a)); }
static inline double distv(const Vec3& a, const Vec3& b){ return normv(a-b); }

struct Quadric {
    double q[10];
    Quadric(){ memset(q, 0, sizeof(q)); }
    inline void addPlane(double a, double b, double c, double d, double w){
        q[0]+=w*a*a; q[1]+=w*a*b; q[2]+=w*a*c; q[3]+=w*a*d;
        q[4]+=w*b*b; q[5]+=w*b*c; q[6]+=w*b*d;
        q[7]+=w*c*c; q[8]+=w*c*d; q[9]+=w*d*d;
    }
    inline void add(const Quadric& o){ for(int i=0;i<10;i++) q[i]+=o.q[i]; }
};
static inline double qeval(const Quadric& Q, const Vec3& p){
    const double x=p.x,y=p.y,z=p.z;
    return Q.q[0]*x*x + 2.0*Q.q[1]*x*y + 2.0*Q.q[2]*x*z + 2.0*Q.q[3]*x
         + Q.q[4]*y*y + 2.0*Q.q[5]*y*z + 2.0*Q.q[6]*y
         + Q.q[7]*z*z + 2.0*Q.q[8]*z + Q.q[9];
}

struct Face { int v[3]; unsigned char alive; Vec3 n; double area; };
struct Vertex { Vec3 p; Quadric Q; double rad; int ver; unsigned char alive; float imp; vector<int> faces; };

static vector<Vertex> Vtx;
static vector<Face> Fac;
static int aliveV = 0, aliveF = 0;
static double hausEps = 0.0, minArea2 = 1e-28;
static const double INF = 1e100;
static chrono::steady_clock::time_point gStart;

struct OutMesh {
    vector<Vec3> V;
    vector<array<int,3>> F;
};
static OutMesh originalMesh;
static vector<OutMesh> snapshots;
static OutMesh answerMesh;

struct FastScanner {
    vector<char> buf; char* p;
    FastScanner(){
        const size_t CH=1<<20; char tmp[CH]; size_t n;
        while((n=fread(tmp,1,CH,stdin))>0) buf.insert(buf.end(), tmp, tmp+n);
        buf.push_back('\0'); p=buf.data();
    }
    inline void skip(){ while(*p==' '||*p=='\n'||*p=='\r'||*p=='\t') ++p; }
    int nextInt(){ skip(); int s=1; if(*p=='-'){s=-1;++p;} int x=0; while(*p>='0'&&*p<='9'){x=x*10+(*p-'0');++p;} return x*s; }
    double nextDouble(){ skip(); char* e; double x=strtod(p,&e); p=e; return x; }
    char nextChar(){ skip(); return *p++; }
};
static inline bool timeExceeded(double s){ return chrono::duration<double>(chrono::steady_clock::now()-gStart).count() > s; }
static inline bool hasVertex(const Face& f,int a){ return f.v[0]==a||f.v[1]==a||f.v[2]==a; }
static inline bool hasBoth(const Face& f,int a,int b){ return hasVertex(f,a)&&hasVertex(f,b); }

static bool recomputeFace(int fid){
    Face& f=Fac[fid];
    Vec3 cr=crossv(Vtx[f.v[1]].p - Vtx[f.v[0]].p, Vtx[f.v[2]].p - Vtx[f.v[0]].p);
    double l2=norm2(cr);
    if(l2<=minArea2) return false;
    double inv=1.0/sqrt(l2);
    f.n=cr*inv; f.area=0.5/inv; return true;
}
static Vec3 faceNormalWithReplacement(const Face& f,int keep,int kill,const Vec3& np,double& l2){
    Vec3 p[3];
    for(int i=0;i<3;i++){ int id=f.v[i]; p[i]=(id==keep||id==kill)?np:Vtx[id].p; }
    Vec3 cr=crossv(p[1]-p[0], p[2]-p[0]); l2=norm2(cr);
    if(l2<=minArea2) return Vec3();
    return cr/sqrt(l2);
}
static void cleanFaces(int v){
    if(!Vtx[v].alive) return;
    vector<int>& a=Vtx[v].faces; size_t w=0;
    for(size_t i=0;i<a.size();++i){ int fid=a[i]; if(fid>=0&&fid<(int)Fac.size()&&Fac[fid].alive&&hasVertex(Fac[fid],v)) a[w++]=fid; }
    a.resize(w);
}
static vector<int> markArr; static int markToken=1;
static void collectNeighbors(int v, vector<int>& out){
    out.clear(); if(!Vtx[v].alive) return; cleanFaces(v);
    ++markToken; if(markToken==INT_MAX){ fill(markArr.begin(),markArr.end(),0); markToken=1; }
    for(int fid: Vtx[v].faces){ const Face& f=Fac[fid]; for(int k=0;k<3;k++){ int u=f.v[k]; if(u!=v&&Vtx[u].alive&&markArr[u]!=markToken){ markArr[u]=markToken; out.push_back(u); } } }
}
static bool solve3x3(const Quadric& Q, Vec3& ans){
    double a00=Q.q[0], a01=Q.q[1], a02=Q.q[2];
    double a11=Q.q[4], a12=Q.q[5], a22=Q.q[7];
    double b0=-Q.q[3], b1=-Q.q[6], b2=-Q.q[8];
    double det=a00*(a11*a22-a12*a12)-a01*(a01*a22-a12*a02)+a02*(a01*a12-a11*a02);
    if(fabs(det)<1e-13) return false;
    double inv=1.0/det;
    double dx=b0*(a11*a22-a12*a12)-a01*(b1*a22-a12*b2)+a02*(b1*a12-a11*b2);
    double dy=a00*(b1*a22-a12*b2)-b0*(a01*a22-a12*a02)+a02*(a01*b2-b1*a02);
    double dz=a00*(a11*b2-b1*a12)-a01*(a01*b2-b1*a02)+b0*(a01*a12-a11*a02);
    ans=Vec3(dx*inv,dy*inv,dz*inv);
    return isfinite(ans.x)&&isfinite(ans.y)&&isfinite(ans.z)&&fabs(ans.x)<=3.0&&fabs(ans.y)<=3.0&&fabs(ans.z)<=3.0;
}
struct Candidate { Vec3 p; double cost; double newRad; };
static bool bestGeometricCandidate(int a,int b,Candidate& best){
    if(!Vtx[a].alive||!Vtx[b].alive) return false;
    Quadric Q=Vtx[a].Q; Q.add(Vtx[b].Q);
    const Vec3 pa=Vtx[a].p,pb=Vtx[b].p; const double ra=Vtx[a].rad, rb=Vtx[b].rad;
    const double imp = 1.0 + 2.0 * max(Vtx[a].imp, Vtx[b].imp);
    best.cost=INF; best.newRad=INF; bool ok=false;
    auto consider=[&](const Vec3& p){
        double nr=max(ra+distv(p,pa), rb+distv(p,pb));
        if(nr>hausEps) return;
        double c=qeval(Q,p); if(c<0&&c>-1e-18) c=0;
        double key=c*imp + 2e-12*norm2(pa-pb) + 2e-14*(nr/max(hausEps,1e-12));
        if(key<best.cost){ best.cost=key; best.p=p; best.newRad=nr; ok=true; }
    };
    Vec3 opt; if(solve3x3(Q,opt)) consider(opt);
    consider(pa); consider(pb); consider((pa+pb)*0.5);
    double d=distv(pa,pb);
    if(d>1e-18){ double t=(rb+d-ra)/(2.0*d); if(t<0)t=0; if(t>1)t=1; consider(pa*(1.0-t)+pb*t); }
    return ok;
}
struct HeapItem { double cost; int a,b,va,vb; bool operator<(const HeapItem& o) const { return cost>o.cost; } };
static priority_queue<HeapItem> heapq;
static void pushEdge(int a,int b){ if(a==b||!Vtx[a].alive||!Vtx[b].alive) return; if(a>b) swap(a,b); Candidate c; if(bestGeometricCandidate(a,b,c)) heapq.push({c.cost,a,b,Vtx[a].ver,Vtx[b].ver}); }
static bool linkCondition(int a,int b){
    if(!Vtx[a].alive||!Vtx[b].alive||a==b) return false; cleanFaces(a); cleanFaces(b);
    int edgeFaces=0;
    if(Vtx[a].faces.size()<Vtx[b].faces.size()){ for(int fid:Vtx[a].faces) if(Fac[fid].alive&&hasBoth(Fac[fid],a,b)) ++edgeFaces; }
    else { for(int fid:Vtx[b].faces) if(Fac[fid].alive&&hasBoth(Fac[fid],a,b)) ++edgeFaces; }
    if(edgeFaces!=2) return false;
    vector<int> na,nb; collectNeighbors(a,na); collectNeighbors(b,nb);
    ++markToken; if(markToken==INT_MAX){ fill(markArr.begin(),markArr.end(),0); markToken=1; }
    for(int u:na) markArr[u]=markToken;
    int common=0; for(int u:nb) if(markArr[u]==markToken) ++common;
    return common==2;
}
static bool geometryOK(int keep,int kill,const Vec3& np){
    const double minDot = -1e-7;
    cleanFaces(keep); cleanFaces(kill);
    for(int pass=0; pass<2; ++pass){
        int v=pass?kill:keep;
        for(int fid: Vtx[v].faces){
            Face& f=Fac[fid]; if(!f.alive) continue;
            bool hk=hasVertex(f,keep), hl=hasVertex(f,kill);
            if(hk&&hl) continue;
            double l2; Vec3 nn=faceNormalWithReplacement(f,keep,kill,np,l2);
            if(l2<=minArea2) return false;
            if(dotv(nn,f.n)<minDot) return false;
        }
    }
    return true;
}
static bool collapseEdge(int a,int b,const Candidate& cand){
    if(!linkCondition(a,b)) return false;
    cleanFaces(a); cleanFaces(b);
    int keep=a,kill=b; if(Vtx[a].faces.size()<Vtx[b].faces.size()) { keep=b; kill=a; }
    if(!geometryOK(keep,kill,cand.p)) return false;
    vector<int> killFaces=Vtx[kill].faces;
    for(int fid:killFaces){
        Face& f=Fac[fid]; if(!f.alive) continue;
        bool hk=hasVertex(f,keep), hl=hasVertex(f,kill); if(!hl) continue;
        if(hk){ f.alive=0; --aliveF; }
        else { for(int i=0;i<3;i++) if(f.v[i]==kill) f.v[i]=keep; Vtx[keep].faces.push_back(fid); }
    }
    Vtx[keep].p=cand.p; Vtx[keep].rad=cand.newRad; Vtx[keep].Q.add(Vtx[kill].Q); Vtx[keep].imp=max(Vtx[keep].imp,Vtx[kill].imp); Vtx[keep].ver++;
    Vtx[kill].alive=0; Vtx[kill].ver++; Vtx[kill].faces.clear(); --aliveV;
    cleanFaces(keep); for(int fid:Vtx[keep].faces) if(Fac[fid].alive) recomputeFace(fid);
    return true;
}
static void pushEdgesAround(int v){ static vector<int> nb; collectNeighbors(v,nb); for(int u:nb) pushEdge(v,u); }
static uint64_t edgeKey(int a,int b){ if(a>b) swap(a,b); return (uint64_t)(uint32_t)a<<32 | (uint32_t)b; }

static void compactCurrent(OutMesh& out){
    int n=(int)Vtx.size(); vector<int> id(n,-1); vector<unsigned char> used(n,0);
    int nf=0; for(const Face& f:Fac) if(f.alive){ ++nf; used[f.v[0]]=used[f.v[1]]=used[f.v[2]]=1; }
    int nv=0; for(int i=0;i<n;i++) if(used[i]) id[i]=nv++;
    out.V.clear(); out.F.clear(); out.V.reserve(nv); out.F.reserve(nf);
    for(int i=0;i<n;i++) if(used[i]) out.V.push_back(Vtx[i].p);
    for(const Face& f:Fac) if(f.alive) out.F.push_back({id[f.v[0]],id[f.v[1]],id[f.v[2]]});
}

// Exact-resolution internal evaluator. It mirrors the statement: 6 axial views,
// flat world-space face normals, perspective-correct depth, 11x11 SSIM windows.
static const int RENDER_R = 1024;
static const double BG_DEPTH = 255.0;
static const double BG_NORM = 127.5;
struct RenderMaps {
    vector<float> data; // [view][channel 0..2 normal RGB, 3 depth][pixel]
    RenderMaps(){ data.assign((size_t)6*4*RENDER_R*RENDER_R, 0.0f); }
    inline float* ptr(int view,int ch){ return data.data()+((size_t)view*4+ch)*(size_t)RENDER_R*RENDER_R; }
    inline const float* ptr(int view,int ch) const { return data.data()+((size_t)view*4+ch)*(size_t)RENDER_R*RENDER_R; }
};
static RenderMaps origMaps;
static bool haveOrigMaps=false;

static inline void projectPoint(const Vec3& p,int view,double& u,double& v,double& z){
    constexpr double D=2.5, f=800.0, c=512.0;
    double sx=0, sy=0;
    if(view==0){ sx=p.y;  sy=p.z; z=D-p.x; }
    else if(view==1){ sx=-p.y; sy=p.z; z=D+p.x; }
    else if(view==2){ sx=-p.x; sy=p.z; z=D-p.y; }
    else if(view==3){ sx=p.x;  sy=p.z; z=D+p.y; }
    else if(view==4){ sx=p.x;  sy=p.y; z=D-p.z; }
    else { sx=-p.x; sy=p.y; z=D+p.z; }
    u=f*sx/z+c; v=f*sy/z+c;
}
static void renderMesh(const OutMesh& mesh, RenderMaps& rm){
    const int R=RENDER_R, P=R*R;
    for(int view=0; view<6; ++view){
        for(int ch=0; ch<3; ++ch) fill(rm.ptr(view,ch), rm.ptr(view,ch)+P, (float)BG_NORM);
        fill(rm.ptr(view,3), rm.ptr(view,3)+P, (float)BG_DEPTH);
    }
    vector<float> zbuf(P);
    const double epsInside=-1e-9;
    for(int view=0; view<6; ++view){
        fill(zbuf.begin(), zbuf.end(), 1e30f);
        float *n0=rm.ptr(view,0), *n1=rm.ptr(view,1), *n2=rm.ptr(view,2), *dd=rm.ptr(view,3);
        for(const auto& af: mesh.F){
            const Vec3 &a=mesh.V[af[0]], &b=mesh.V[af[1]], &c=mesh.V[af[2]];
            Vec3 cr=crossv(b-a,c-a); double len=normv(cr); if(!(len>0)) continue;
            Vec3 n=cr/len;
            double x0,y0,z0,x1,y1,z1,x2,y2,z2;
            projectPoint(a,view,x0,y0,z0); projectPoint(b,view,x1,y1,z1); projectPoint(c,view,x2,y2,z2);
            if(z0<=1e-12||z1<=1e-12||z2<=1e-12) continue;
            int xmin=max(0,(int)floor(min(x0,min(x1,x2))-0.5));
            int xmax=min(R-1,(int)ceil(max(x0,max(x1,x2))+0.5));
            int ymin=max(0,(int)floor(min(y0,min(y1,y2))-0.5));
            int ymax=min(R-1,(int)ceil(max(y0,max(y1,y2))+0.5));
            if(xmin>xmax||ymin>ymax) continue;
            double den=(y1-y2)*(x0-x2)+(x2-x1)*(y0-y2); if(fabs(den)<1e-18) continue;
            double invDen=1.0/den, iz0=1.0/z0, iz1=1.0/z1, iz2=1.0/z2;
            float cn0=(float)((n.x+1.0)*127.5), cn1=(float)((n.y+1.0)*127.5), cn2=(float)((n.z+1.0)*127.5);
            for(int yy=ymin; yy<=ymax; ++yy){
                double py=yy+0.5; int base=yy*R;
                for(int xx=xmin; xx<=xmax; ++xx){
                    double px=xx+0.5;
                    double w0=((y1-y2)*(px-x2)+(x2-x1)*(py-y2))*invDen;
                    double w1=((y2-y0)*(px-x2)+(x0-x2)*(py-y2))*invDen;
                    double w2=1.0-w0-w1;
                    if(w0<epsInside||w1<epsInside||w2<epsInside) continue;
                    double invz=w0*iz0+w1*iz1+w2*iz2; if(invz<=0) continue;
                    float dep=(float)(1.0/invz); int idx=base+xx;
                    if(dep<zbuf[idx]){ zbuf[idx]=dep; dd[idx]=dep; n0[idx]=cn0; n1[idx]=cn1; n2[idx]=cn2; }
                }
            }
        }
    }
}
static double ssimChannel(const float* A,const float* B,const unsigned char* mask){
    const int R=RENDER_R, PAD=5, PR=R+2*PAD, W=PR+1;
    static vector<double> IA,IB,IA2,IB2,IAB;
    size_t sz=(size_t)W*W; IA.assign(sz,0.0); IB.assign(sz,0.0); IA2.assign(sz,0.0); IB2.assign(sz,0.0); IAB.assign(sz,0.0);
    for(int py=1; py<=PR; ++py){
        int sy=py-1-PAD; if(sy<0) sy=0; else if(sy>=R) sy=R-1;
        double sa=0,sb=0,sa2=0,sb2=0,sab=0;
        int row=sy*R;
        for(int px=1; px<=PR; ++px){
            int sx=px-1-PAD; if(sx<0) sx=0; else if(sx>=R) sx=R-1;
            double a=A[row+sx], b=B[row+sx];
            sa+=a; sb+=b; sa2+=a*a; sb2+=b*b; sab+=a*b;
            size_t idx=(size_t)py*W+px, up=(size_t)(py-1)*W+px;
            IA[idx]=IA[up]+sa; IB[idx]=IB[up]+sb; IA2[idx]=IA2[up]+sa2; IB2[idx]=IB2[up]+sb2; IAB[idx]=IAB[up]+sab;
        }
    }
    auto rect=[&](const vector<double>& I,int x0,int y0,int x1,int y1)->double{
        return I[(size_t)y1*W+x1]-I[(size_t)y0*W+x1]-I[(size_t)y1*W+x0]+I[(size_t)y0*W+x0];
    };
    const double C1=(0.01*255.0)*(0.01*255.0), C2=(0.03*255.0)*(0.03*255.0), invN=1.0/121.0;
    double sum=0.0; long long cnt=0;
    for(int y=0; y<R; ++y){ int row=y*R; int y0=y, y1=y+11;
        for(int x=0; x<R; ++x){ if(!mask[row+x]) continue; int x0=x, x1=x+11;
            double sA=rect(IA,x0,y0,x1,y1), sB=rect(IB,x0,y0,x1,y1);
            double sA2=rect(IA2,x0,y0,x1,y1), sB2=rect(IB2,x0,y0,x1,y1), sAB=rect(IAB,x0,y0,x1,y1);
            double muA=sA*invN, muB=sB*invN;
            double varA=max(0.0,sA2*invN-muA*muA), varB=max(0.0,sB2*invN-muB*muB), cov=sAB*invN-muA*muB;
            double den=(muA*muA+muB*muB+C1)*(varA+varB+C2);
            if(den>0.0){ double val=((2.0*muA*muB+C1)*(2.0*cov+C2))/den; sum+=val; }
            ++cnt;
        }
    }
    if(cnt==0) return 1.0;
    double r=sum/(double)cnt; if(r<0) r=0; if(r>1) r=1; return r;
}
static double evaluateVisual(const OutMesh& mesh){
    if(!haveOrigMaps) return 1.0;
    RenderMaps cur; renderMesh(mesh,cur);
    const int R=RENDER_R, P=R*R; vector<unsigned char> mask(P);
    double total=0.0;
    for(int view=0; view<6; ++view){
        const float* od=origMaps.ptr(view,3); const float* cd=cur.ptr(view,3);
        for(int i=0;i<P;i++) mask[i]=(od[i]!=(float)BG_DEPTH || cd[i]!=(float)BG_DEPTH);
        double ns=0.0; for(int ch=0; ch<3; ++ch) ns += ssimChannel(origMaps.ptr(view,ch), cur.ptr(view,ch), mask.data());
        ns/=3.0; double ds=ssimChannel(od,cd,mask.data());
        total += 0.5*ns + 0.5*ds;
    }
    return total/6.0;
}

static void loadMesh(){
    FastScanner fs; int n=fs.nextInt(), m=fs.nextInt();
    Vtx.assign(n,Vertex()); Fac.assign(m,Face()); aliveV=n; aliveF=m;
    double xmin=1e100,ymin=1e100,zmin=1e100,xmax=-1e100,ymax=-1e100,zmax=-1e100;
    for(int i=0;i<n;i++){
        (void)fs.nextChar(); double x=fs.nextDouble(), y=fs.nextDouble(), z=fs.nextDouble();
        Vtx[i].p=Vec3(x,y,z); Vtx[i].rad=0.0; Vtx[i].ver=0; Vtx[i].alive=1; Vtx[i].imp=0.0f;
        xmin=min(xmin,x); ymin=min(ymin,y); zmin=min(zmin,z); xmax=max(xmax,x); ymax=max(ymax,y); zmax=max(zmax,z);
    }
    double dx=xmax-xmin,dy=ymax-ymin,dz=zmax-zmin; double diag=sqrt(dx*dx+dy*dy+dz*dz); if(!(diag>0)) diag=1.0;
    hausEps=0.05*diag*0.999999; minArea2=max(1e-30, diag*diag*diag*diag*1e-28);
    for(int i=0;i<m;i++){
        (void)fs.nextChar(); int a=fs.nextInt()-1,b=fs.nextInt()-1,c=fs.nextInt()-1;
        Fac[i].v[0]=a; Fac[i].v[1]=b; Fac[i].v[2]=c; Fac[i].alive=1;
        Vtx[a].faces.push_back(i); Vtx[b].faces.push_back(i); Vtx[c].faces.push_back(i);
    }
    vector<Vec3> nsum(n,Vec3()); vector<int> ncnt(n,0);
    for(int i=0;i<m;i++){
        recomputeFace(i); Face& f=Fac[i]; const Vec3& p=Vtx[f.v[0]].p; double d=-dotv(f.n,p);
        double vieww=fabs(f.n.x)+fabs(f.n.y)+fabs(f.n.z);
        double w=max(f.area,1e-14)*(0.80+0.38*vieww);
        for(int k=0;k<3;k++) Vtx[f.v[k]].Q.addPlane(f.n.x,f.n.y,f.n.z,d,w);
        for(int k=0;k<3;k++){ int v=f.v[k]; nsum[v]=nsum[v]+f.n; ++ncnt[v]; }
    }
    for(int i=0;i<n;i++){
        double curv=0.0; if(ncnt[i]>0){ double ml=normv(nsum[i])/(double)ncnt[i]; curv=1.0-max(0.0,min(1.0,ml)); }
        double ext=0.0, near=hausEps*1.5; const Vec3& p=Vtx[i].p;
        if(fabs(p.x-xmin)<near||fabs(p.x-xmax)<near||fabs(p.y-ymin)<near||fabs(p.y-ymax)<near||fabs(p.z-zmin)<near||fabs(p.z-zmax)<near) ext=0.25;
        Vtx[i].imp=(float)max(curv,ext);
    }
}
static void buildInitialQueue(){
    vector<uint64_t> edges; edges.reserve((size_t)Fac.size()*3);
    for(const Face& f:Fac){ edges.push_back(edgeKey(f.v[0],f.v[1])); edges.push_back(edgeKey(f.v[1],f.v[2])); edges.push_back(edgeKey(f.v[2],f.v[0])); }
    sort(edges.begin(),edges.end()); edges.erase(unique(edges.begin(),edges.end()),edges.end());
    for(uint64_t e:edges) pushEdge((int)(e>>32),(int)(e&0xffffffffu));
}
static vector<double> makeCheckpoints(int n){
    if(n<=20) return {0.90,0.80,0.70,0.60,0.50,0.45,0.40};
    if(n<=6000) return {0.30,0.22,0.17,0.135,0.115,0.102,0.094,0.088,0.082,0.076,0.070,0.064};
    return {0.25,0.18,0.14,0.115,0.102,0.094,0.088,0.082,0.076,0.070,0.064,0.060};
}
static double hardTargetRatio(int n){
    if(n<=20) return 0.01;
    if(n<=6000) return 0.060;
    if(n<=50000) return 0.058;
    return 0.056;
}
static void simplify(){
    const int n0=(int)Vtx.size();
    compactCurrent(originalMesh);
    renderMesh(originalMesh, origMaps); haveOrigMaps=true;
    if(n0<=4){ answerMesh=originalMesh; return; }
    markArr.assign(n0,0); buildInitialQueue();
    vector<double> cp=makeCheckpoints(n0); size_t cpIdx=0;
    int hardTarget=max(4,(int)ceil(n0*hardTargetRatio(n0)));
    long long popCount=0; const long long popLimit=max<long long>(2500000LL,(long long)n0*130LL);
    while(aliveV>hardTarget && !heapq.empty() && popCount<popLimit){
        if((popCount&8191LL)==0 && timeExceeded(15.7)) break;
        HeapItem it=heapq.top(); heapq.pop(); ++popCount;
        int a=it.a,b=it.b; if(a<0||b<0||a>=n0||b>=n0) continue;
        if(!Vtx[a].alive||!Vtx[b].alive) continue;
        if(Vtx[a].ver!=it.va||Vtx[b].ver!=it.vb) continue;
        Candidate cand; if(!bestGeometricCandidate(a,b,cand)) continue;
        if(cand.cost > it.cost*1.000001 + 1e-18){ heapq.push({cand.cost,min(a,b),max(a,b),Vtx[min(a,b)].ver,Vtx[max(a,b)].ver}); continue; }
        if(collapseEdge(a,b,cand)){
            int keep = Vtx[a].alive ? a : b; if(!Vtx[keep].alive) keep=(keep==a?b:a); pushEdgesAround(keep);
        }
        double ratio=(double)aliveV/(double)n0;
        while(cpIdx<cp.size() && ratio<=cp[cpIdx]){ OutMesh s; compactCurrent(s); snapshots.push_back(std::move(s)); ++cpIdx; }
    }
    OutMesh cur; compactCurrent(cur);
    if(snapshots.empty() || snapshots.back().V.size()!=cur.V.size()) snapshots.push_back(std::move(cur));

    // Start from original as guaranteed valid fallback; exact-test a conservative snapshot,
    // then search from most aggressive to least aggressive.
    answerMesh=originalMesh;
    const double PASS=0.9020;
    vector<int> order;
    int safeIdx=-1;
    for(int i=0;i<(int)snapshots.size();++i){
        double r=(double)snapshots[i].V.size()/max(1,n0);
        if(r>=0.16){ safeIdx=i; break; }
    }
    if(safeIdx<0 && !snapshots.empty()) safeIdx=0;
    if(safeIdx>=0 && !timeExceeded(17.8)){
        double s=evaluateVisual(snapshots[safeIdx]);
        if(s>=PASS) answerMesh=snapshots[safeIdx];
    }
    for(int i=(int)snapshots.size()-1; i>=0; --i){
        if(i==safeIdx) continue;
        if(timeExceeded(20.15)) break;
        double s=evaluateVisual(snapshots[i]);
        if(s>=PASS){ answerMesh=snapshots[i]; break; }
    }
}
static void saveOut(const OutMesh& om){
    string out; out.reserve((size_t)om.V.size()*44+(size_t)om.F.size()*28+64); char line[128];
    out.append(line,snprintf(line,sizeof(line),"%d %d\n",(int)om.V.size(),(int)om.F.size()));
    for(const Vec3& p:om.V) out.append(line,snprintf(line,sizeof(line),"v %.10g %.10g %.10g\n",p.x,p.y,p.z));
    for(const auto& f:om.F) out.append(line,snprintf(line,sizeof(line),"f %d %d %d\n",f[0]+1,f[1]+1,f[2]+1));
    fwrite(out.data(),1,out.size(),stdout);
}
int main(){
    gStart=chrono::steady_clock::now();
    loadMesh();
    simplify();
    saveOut(answerMesh);
    return 0;
}
