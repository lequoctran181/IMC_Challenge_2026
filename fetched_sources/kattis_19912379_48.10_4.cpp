#include <bits/stdc++.h>
using namespace std;

// W5 guarded full-source variant.  It keeps the accepted-family philosophy:
// conservative topology-preserving edge collapses, vertex-only Hausdorff cluster
// guards, and an internal six-view SSIM proxy to choose the most aggressive safe
// snapshot.  No external dependencies; C++17.

struct Vec3{double x=0,y=0,z=0;};
static inline Vec3 operator+(const Vec3&a,const Vec3&b){return {a.x+b.x,a.y+b.y,a.z+b.z};}
static inline Vec3 operator-(const Vec3&a,const Vec3&b){return {a.x-b.x,a.y-b.y,a.z-b.z};}
static inline Vec3 operator*(const Vec3&a,double s){return {a.x*s,a.y*s,a.z*s};}
static inline Vec3 operator/(const Vec3&a,double s){return {a.x/s,a.y/s,a.z/s};}
static inline double dot3(const Vec3&a,const Vec3&b){return a.x*b.x+a.y*b.y+a.z*b.z;}
static inline Vec3 cross3(const Vec3&a,const Vec3&b){return {a.y*b.z-a.z*b.y,a.z*b.x-a.x*b.z,a.x*b.y-a.y*b.x};}
static inline double norm2(const Vec3&a){return dot3(a,a);} 
static inline double norm3(const Vec3&a){return sqrt(max(0.0,norm2(a)));}
static inline Vec3 normalize(Vec3 a){double n=norm3(a);return n>1e-300?a/n:Vec3{};}

struct Tri{int a=0,b=0,c=0;};
struct Face{int v[3]={0,0,0}; unsigned char active=1; Vec3 n{};};

struct FastScanner{
    static constexpr int SZ=1<<20; int idx=0,sz=0; char buf[SZ];
    inline char get(){ if(idx>=sz){ sz=(int)fread(buf,1,SZ,stdin); idx=0; if(!sz) return 0; } return buf[idx++]; }
    inline char prefix(){ char c=get(); while(c==' '||c=='\n'||c=='\r'||c=='\t') c=get(); return c; }
    bool readInt(int&out){ char c=prefix(); if(!c) return false; int s=1; if(c=='-') s=-1,c=get(); int x=0; while(c>='0'&&c<='9'){ x=x*10+(c-'0'); c=get(); } out=x*s; return true; }
    bool readDouble(double&out){ char c=prefix(); if(!c) return false; int s=1; if(c=='-') s=-1,c=get(); double x=0; while(c>='0'&&c<='9'){ x=x*10+(c-'0'); c=get(); } if(c=='.'){ double f=.1; c=get(); while(c>='0'&&c<='9'){ x+=f*(c-'0'); f*=.1; c=get(); }} if(c=='e'||c=='E'){ int es=1,ev=0; c=get(); if(c=='-') es=-1,c=get(); else if(c=='+') c=get(); while(c>='0'&&c<='9'){ ev=ev*10+(c-'0'); c=get(); } x*=pow(10.0,es*ev); } out=s*x; return true; }
};

struct Timer{ chrono::steady_clock::time_point st; Timer(){st=chrono::steady_clock::now();} double elapsed()const{return chrono::duration<double>(chrono::steady_clock::now()-st).count();}};
static Timer GTIMER;
static int ORIGN=0,ORIGM=0;
static vector<Vec3> ORIGV;
static vector<Tri> ORIGF;
static double GLOBAL_DIAG=1.0;

static inline uint64_t ekey(int a,int b){ if(a>b) swap(a,b); return (uint64_t)(uint32_t)a<<32 | (uint32_t)b; }
static inline uint64_t fkey(int a,int b,int c){ if(a>b) swap(a,b); if(b>c) swap(b,c); if(a>b) swap(a,b); return ((uint64_t)(uint32_t)a<<42)^((uint64_t)(uint32_t)b<<21)^(uint32_t)c; }

struct RenderMaps{ int R=0; vector<float> depth; vector<float> normal; vector<unsigned char> mask; };

static inline void project_point(const Vec3&p,int view,int R,double&u,double&v,double&z){
    const double D=2.5; double f=800.0*(double)R/1024.0, c=.5*(double)R; double sx=0,sy=0;
    if(view==0){sx=p.y; sy=p.z; z=D-p.x;}
    else if(view==1){sx=-p.y; sy=p.z; z=D+p.x;}
    else if(view==2){sx=-p.x; sy=p.z; z=D-p.y;}
    else if(view==3){sx=p.x; sy=p.z; z=D+p.y;}
    else if(view==4){sx=p.x; sy=p.y; z=D-p.z;}
    else {sx=-p.x; sy=p.y; z=D+p.z;}
    u=f*sx/z+c; v=f*sy/z+c;
}

static void render_mesh(const vector<Vec3>&V,const vector<Tri>&F,RenderMaps&rm,int R){
    const int PIX=R*R, VIEWS=6; rm.R=R; rm.depth.assign((size_t)VIEWS*PIX,255.f); rm.normal.assign((size_t)VIEWS*PIX*3,127.5f); rm.mask.assign((size_t)VIEWS*PIX,0);
    if(V.empty()||F.empty()) return;
    vector<float> U(V.size()),W(V.size()),Z(V.size());
    vector<Vec3> FN(F.size());
    for(size_t i=0;i<F.size();++i){ const Tri&t=F[i]; if(t.a<0||t.b<0||t.c<0||t.a>=(int)V.size()||t.b>=(int)V.size()||t.c>=(int)V.size()) continue; Vec3 n=cross3(V[t.b]-V[t.a],V[t.c]-V[t.a]); FN[i]=normalize(n); }
    for(int view=0;view<VIEWS;++view){
        for(size_t i=0;i<V.size();++i){ double u,v,z; project_point(V[i],view,R,u,v,z); U[i]=(float)u; W[i]=(float)v; Z[i]=(float)z; }
        float*zbuf=rm.depth.data()+(size_t)view*PIX; float*nbuf=rm.normal.data()+(size_t)view*PIX*3; unsigned char*mbuf=rm.mask.data()+(size_t)view*PIX;
        for(size_t ti=0;ti<F.size();++ti){ const Tri&t=F[ti]; int ia=t.a,ib=t.b,ic=t.c; if(ia<0||ib<0||ic<0||ia>=(int)V.size()||ib>=(int)V.size()||ic>=(int)V.size()) continue;
            float x0=U[ia],x1=U[ib],x2=U[ic], y0=W[ia],y1=W[ib],y2=W[ic], z0=Z[ia],z1=Z[ib],z2=Z[ic];
            if(!(z0>0&&z1>0&&z2>0)) continue;
            int xmin=max(0,(int)floor(min(x0,min(x1,x2))-.5f)), xmax=min(R-1,(int)ceil(max(x0,max(x1,x2))+.5f));
            int ymin=max(0,(int)floor(min(y0,min(y1,y2))-.5f)), ymax=min(R-1,(int)ceil(max(y0,max(y1,y2))+.5f));
            if(xmin>xmax||ymin>ymax) continue; float den=(y1-y2)*(x0-x2)+(x2-x1)*(y0-y2); if(fabs(den)<1e-20f) continue; float inv=1.f/den;
            Vec3 n=FN[ti]; float nr=(float)((n.x+1.)*127.5), ng=(float)((n.y+1.)*127.5), nb=(float)((n.z+1.)*127.5);
            for(int y=ymin;y<=ymax;++y){ float py=(float)y+.5f; int row=y*R; for(int x=xmin;x<=xmax;++x){ float px=(float)x+.5f; float a=((y1-y2)*(px-x2)+(x2-x1)*(py-y2))*inv; float b=((y2-y0)*(px-x2)+(x0-x2)*(py-y2))*inv; float c=1.f-a-b; if(a>=-1e-6f&&b>=-1e-6f&&c>=-1e-6f){ float zz=1.f/(a/z0+b/z1+c/z2); int id=row+x; if(zz<zbuf[id]){ zbuf[id]=zz; mbuf[id]=1; float*q=nbuf+id*3; q[0]=nr; q[1]=ng; q[2]=nb; } } }}
        }
    }
}

static inline double rectsum(const vector<double>&I,int W,int x0,int y0,int x1,int y1){ return I[(size_t)y1*W+x1]-I[(size_t)y0*W+x1]-I[(size_t)y1*W+x0]+I[(size_t)y0*W+x0]; }

static double ssim_channel(const float*A,int strideA,const float*B,int strideB,const unsigned char*fg,int R,vector<double>&IA,vector<double>&IB,vector<double>&IA2,vector<double>&IB2,vector<double>&IAB){
    int W=R+1; size_t SZ=(size_t)W*W; fill(IA.begin(),IA.begin()+SZ,0); fill(IB.begin(),IB.begin()+SZ,0); fill(IA2.begin(),IA2.begin()+SZ,0); fill(IB2.begin(),IB2.begin()+SZ,0); fill(IAB.begin(),IAB.begin()+SZ,0);
    for(int y=1;y<=R;++y){ double sa=0,sb=0,sa2=0,sb2=0,sab=0; int row=(y-1)*R; for(int x=1;x<=R;++x){ int p=row+x-1; double a=A[(size_t)p*strideA], b=B[(size_t)p*strideB]; sa+=a; sb+=b; sa2+=a*a; sb2+=b*b; sab+=a*b; size_t id=(size_t)y*W+x, up=(size_t)(y-1)*W+x; IA[id]=IA[up]+sa; IB[id]=IB[up]+sb; IA2[id]=IA2[up]+sa2; IB2[id]=IB2[up]+sb2; IAB[id]=IAB[up]+sab; }}
    const int rad=5, area=121; const double c1=6.5025,c2=58.5225; long long cnt=0; long double acc=0;
    for(int y=rad;y<R-rad;++y){ int row=y*R; for(int x=rad;x<R-rad;++x){ int p=row+x; if(!fg[p]) continue; int x0=x-rad,y0=y-rad,x1=x+rad+1,y1=y+rad+1; double sa=rectsum(IA,W,x0,y0,x1,y1), sb=rectsum(IB,W,x0,y0,x1,y1); double sa2=rectsum(IA2,W,x0,y0,x1,y1), sb2=rectsum(IB2,W,x0,y0,x1,y1), sab=rectsum(IAB,W,x0,y0,x1,y1); double ma=sa/area, mb=sb/area; double va=max(0.0,sa2/area-ma*ma), vb=max(0.0,sb2/area-mb*mb), cov=sab/area-ma*mb; double den=(ma*ma+mb*mb+c1)*(va+vb+c2); double val=den!=0?((2*ma*mb+c1)*(2*cov+c2)/den):1.0; acc+=val; ++cnt; }}
    return cnt? (double)(acc/cnt):1.0;
}

static double proxy_score(const RenderMaps&orig,const vector<Vec3>&V,const vector<Tri>&F,int R){
    RenderMaps cand; render_mesh(V,F,cand,R); int PIX=R*R,W=R+1; size_t SZ=(size_t)W*W; vector<double>IA(SZ),IB(SZ),IA2(SZ),IB2(SZ),IAB(SZ); double total=0; vector<unsigned char>fg(PIX);
    for(int view=0;view<6;++view){ const float*od=orig.depth.data()+(size_t)view*PIX; const float*cd=cand.depth.data()+(size_t)view*PIX; const float*on=orig.normal.data()+(size_t)view*PIX*3; const float*cn=cand.normal.data()+(size_t)view*PIX*3; const unsigned char*om=orig.mask.data()+(size_t)view*PIX; const unsigned char*cm=cand.mask.data()+(size_t)view*PIX; for(int i=0;i<PIX;++i) fg[i]=(om[i]||cm[i])?1:0; double ns=0; for(int ch=0;ch<3;++ch) ns+=ssim_channel(on+ch,3,cn+ch,3,fg.data(),R,IA,IB,IA2,IB2,IAB); ns/=3.; double ds=ssim_channel(od,1,cd,1,fg.data(),R,IA,IB,IA2,IB2,IAB); total+=0.5*(ns+ds); }
    return total/6.0;
}

static bool validate_mesh(const vector<Vec3>&V,const vector<Tri>&F){
    if(V.empty()||F.empty()||V.size()>(size_t)ORIGN) return false; double eps=max(1e-30,1e-24*pow(GLOBAL_DIAG,4)); vector<uint64_t> edges; edges.reserve(F.size()*3); vector<uint64_t> faces; faces.reserve(F.size());
    vector<vector<int>> adj(V.size());
    for(const auto&t:F){ if(t.a<0||t.b<0||t.c<0||t.a>=(int)V.size()||t.b>=(int)V.size()||t.c>=(int)V.size()) return false; if(t.a==t.b||t.a==t.c||t.b==t.c) return false; if(norm2(cross3(V[t.b]-V[t.a],V[t.c]-V[t.a]))<=eps) return false; edges.push_back(ekey(t.a,t.b)); edges.push_back(ekey(t.b,t.c)); edges.push_back(ekey(t.c,t.a)); faces.push_back(fkey(t.a,t.b,t.c)); adj[t.a].push_back(t.b); adj[t.a].push_back(t.c); adj[t.b].push_back(t.a); adj[t.b].push_back(t.c); adj[t.c].push_back(t.a); adj[t.c].push_back(t.b); }
    sort(faces.begin(),faces.end()); if(adjacent_find(faces.begin(),faces.end())!=faces.end()) return false; sort(edges.begin(),edges.end()); for(size_t i=0;i<edges.size();){ size_t j=i+1; while(j<edges.size()&&edges[j]==edges[i]) ++j; if(j-i!=2) return false; i=j; }
    vector<char> vis(V.size(),0); queue<int>q; int start=-1; for(size_t i=0;i<V.size();++i) if(!adj[i].empty()){start=(int)i;break;} if(start<0) return false; vis[start]=1; q.push(start); int cnt=0; while(!q.empty()){int u=q.front();q.pop();++cnt; for(int v:adj[u]) if(!vis[v]) vis[v]=1,q.push(v);} for(size_t i=0;i<V.size();++i) if(!adj[i].empty()&&!vis[i]) return false; return true;
}

struct Quadric{
    double q[10]; Quadric(){memset(q,0,sizeof(q));}
    void addPlane(double a,double b,double c,double d,double w){q[0]+=w*a*a;q[1]+=w*a*b;q[2]+=w*a*c;q[3]+=w*a*d;q[4]+=w*b*b;q[5]+=w*b*c;q[6]+=w*b*d;q[7]+=w*c*c;q[8]+=w*c*d;q[9]+=w*d*d;}
    void add(const Quadric&o){for(int i=0;i<10;++i) q[i]+=o.q[i];}
    double eval(const Vec3&p)const{double x=p.x,y=p.y,z=p.z; return q[0]*x*x+2*q[1]*x*y+2*q[2]*x*z+2*q[3]*x+q[4]*y*y+2*q[5]*y*z+2*q[6]*y+q[7]*z*z+2*q[8]*z+q[9];}
};

struct Snapshot{ vector<Vec3> V; vector<Tri> F; double ratio=1; double ssim=-2; };

class Simplifier{
public:
    int N=0,M=0,activeV=0,activeF=0; double tau=0,tau2=0,minArea2=0; vector<Vec3>P; vector<Face>F; vector<vector<int>>inc; vector<Quadric>Q; vector<array<double,3>>bbMin,bbMax; vector<unsigned char>alive,locked; vector<int>ver,markA,markB,tmp; int stampA=1,stampB=1; vector<Snapshot>snaps;
    struct Cand{double s; int u,v,vu,vv; bool operator<(const Cand&o)const{return s>o.s;}}; priority_queue<Cand>pq;
    Simplifier(const vector<Vec3>&V,const vector<Tri>&T){build(V,T);}    
    bool timeOK(double lim=17.4)const{return GTIMER.elapsed()<lim;}
    inline bool contains(int fid,int u)const{const Face&f=F[fid];return f.v[0]==u||f.v[1]==u||f.v[2]==u;}
    Vec3 fNormal(const Face&f)const{return normalize(cross3(P[f.v[1]]-P[f.v[0]],P[f.v[2]]-P[f.v[0]]));}
    void build(const vector<Vec3>&V,const vector<Tri>&T){N=(int)V.size();M=(int)T.size();P=V;F.assign(M,Face());vector<int>deg(N); for(int i=0;i<M;++i){F[i].v[0]=T[i].a;F[i].v[1]=T[i].b;F[i].v[2]=T[i].c;F[i].active=1; for(int k=0;k<3;++k) if(F[i].v[k]>=0&&F[i].v[k]<N) ++deg[F[i].v[k]];}inc.assign(N,{});for(int i=0;i<N;++i)inc[i].reserve(deg[i]+8);activeF=M;activeV=N;for(int i=0;i<M;++i){F[i].n=fNormal(F[i]);for(int k=0;k<3;++k)inc[F[i].v[k]].push_back(i);}Q.assign(N,Quadric());for(int i=0;i<M;++i){Face&f=F[i];Vec3 cr=cross3(P[f.v[1]]-P[f.v[0]],P[f.v[2]]-P[f.v[0]]);double a2=norm3(cr);if(a2<=1e-300)continue;Vec3 n=cr/a2;double d=-dot3(n,P[f.v[0]]);double w=max(1e-10,0.5*a2);for(int k=0;k<3;++k)Q[f.v[k]].addPlane(n.x,n.y,n.z,d,w);}bbMin.resize(N);bbMax.resize(N);for(int i=0;i<N;++i){bbMin[i]={P[i].x,P[i].y,P[i].z};bbMax[i]=bbMin[i];}alive.assign(N,1);locked.assign(N,0);{Vec3 gmn=P[0],gmx=P[0];for(auto &pp:P){gmn.x=min(gmn.x,pp.x);gmn.y=min(gmn.y,pp.y);gmn.z=min(gmn.z,pp.z);gmx.x=max(gmx.x,pp.x);gmx.y=max(gmx.y,pp.y);gmx.z=max(gmx.z,pp.z);}double leps=max(1e-12,1e-10*GLOBAL_DIAG);for(int ii=0;ii<N;++ii){if(fabs(P[ii].x-gmn.x)<=leps||fabs(P[ii].x-gmx.x)<=leps||fabs(P[ii].y-gmn.y)<=leps||fabs(P[ii].y-gmx.y)<=leps||fabs(P[ii].z-gmn.z)<=leps||fabs(P[ii].z-gmx.z)<=leps)locked[ii]=1;}}ver.assign(N,0);markA.assign(N,0);markB.assign(N,0);tmp.assign(N,-1);tau=0.05*GLOBAL_DIAG*0.992; tau2=tau*tau; minArea2=max(1e-30,1e-24*pow(GLOBAL_DIAG,4)); initHeap();}
    void initHeap(){vector<uint64_t> E;E.reserve((size_t)M*3);for(const auto&f:F){E.push_back(ekey(f.v[0],f.v[1]));E.push_back(ekey(f.v[1],f.v[2]));E.push_back(ekey(f.v[2],f.v[0]));}sort(E.begin(),E.end());E.erase(unique(E.begin(),E.end()),E.end());for(auto k:E)pushEdge((int)(k>>32),(int)(k&0xffffffffu));}
    void cleanup(int u){if(u<0||u>=N||!alive[u])return;vector<int>&v=inc[u];int w=0;for(int fid:v)if(fid>=0&&fid<M&&F[fid].active&&contains(fid,u))v[w++]=fid;v.resize(w);}    
    void maybeCleanup(int u){if(u<0||u>=N||!alive[u])return;if(inc[u].size()<128)return;int good=0;for(int fid:inc[u])if(F[fid].active&&contains(fid,u))++good;if((int)inc[u].size()>good*3+64)cleanup(u);}    
    void mergeBox(int u,int v,array<double,3>&mn,array<double,3>&mx){for(int k=0;k<3;++k){mn[k]=min(bbMin[u][k],bbMin[v][k]);mx[k]=max(bbMax[u][k],bbMax[v][k]);}}
    bool coverBox(const Vec3&p,const array<double,3>&mn,const array<double,3>&mx){for(int m=0;m<8;++m){double x=(m&1)?mx[0]:mn[0],y=(m&2)?mx[1]:mn[1],z=(m&4)?mx[2]:mn[2];if((p.x-x)*(p.x-x)+(p.y-y)*(p.y-y)+(p.z-z)*(p.z-z)>tau2)return false;}return true;}
    bool solveOpt(const Quadric&q,Vec3&out){double a00=q.q[0],a01=q.q[1],a02=q.q[2],a11=q.q[4],a12=q.q[5],a22=q.q[7];double b0=q.q[3],b1=q.q[6],b2=q.q[8];double det=a00*(a11*a22-a12*a12)-a01*(a01*a22-a12*a02)+a02*(a01*a12-a11*a02);if(fabs(det)<1e-14)return false;double i00=(a11*a22-a12*a12)/det,i01=(a02*a12-a01*a22)/det,i02=(a01*a12-a02*a11)/det,i11=(a00*a22-a02*a02)/det,i12=(a01*a02-a00*a12)/det,i22=(a00*a11-a01*a01)/det;out.x=-(i00*b0+i01*b1+i02*b2);out.y=-(i01*b0+i11*b1+i12*b2);out.z=-(i02*b0+i12*b1+i22*b2);return isfinite(out.x)&&isfinite(out.y)&&isfinite(out.z)&&fabs(out.x)<1e6&&fabs(out.y)<1e6&&fabs(out.z)<1e6;}
    struct Eval{bool ok=false;double score=1e300;Vec3 pos;};
    Eval evalEdge(int u,int v){Eval e;if(u==v||!alive[u]||!alive[v])return e;if(locked[u]&&locked[v])return e;array<double,3>mn,mx;mergeBox(u,v,mn,mx);Quadric q=Q[u];q.add(Q[v]);Vec3 cands[8];int cc=0;if(locked[u]||locked[v]){cands[cc++]=locked[u]?P[u]:P[v];}else{if(solveOpt(q,cands[cc]))++cc;cands[cc++]=(P[u]+P[v])*.5;cands[cc++]=P[u];cands[cc++]=P[v];cands[cc++]={.5*(mn[0]+mx[0]),.5*(mn[1]+mx[1]),.5*(mn[2]+mx[2])};cands[cc++]=P[u]*.75+P[v]*.25;cands[cc++]=P[u]*.25+P[v]*.75;}double l2=norm2(P[u]-P[v]);for(int i=0;i<cc;++i){Vec3 p=cands[i];if(!coverBox(p,mn,mx))continue;double val=q.eval(p);if(val<0&&val>-1e-9)val=0;double sc=val+1e-8*l2+1e-12*norm2(p);if(sc<e.score){e.ok=true;e.score=sc;e.pos=p;}}return e;}
    void pushEdge(int u,int v){if(u==v||u<0||v<0||u>=N||v>=N||!alive[u]||!alive[v])return;if(u>v)swap(u,v);Eval e=evalEdge(u,v);if(e.ok)pq.push({e.score,u,v,ver[u],ver[v]});}
    bool linkOK(int u,int v){if(u==v||!alive[u]||!alive[v])return false;maybeCleanup(u);maybeCleanup(v);if(++stampA>1000000000){fill(markA.begin(),markA.end(),0);stampA=1;}if(++stampB>1000000000){fill(markB.begin(),markB.end(),0);stampB=1;}int edgeFaces=0;for(int fid:inc[u]){Face&f=F[fid];if(!f.active||!contains(fid,u))continue;bool hv=false;for(int k=0;k<3;++k){int x=f.v[k];if(x==v)hv=true;if(x!=u&&x!=v)markA[x]=stampA;}if(hv)++edgeFaces;}if(edgeFaces!=2)return false;int common=0;for(int fid:inc[v]){Face&f=F[fid];if(!f.active||!contains(fid,v))continue;for(int k=0;k<3;++k){int x=f.v[k];if(x==u||x==v)continue;if(markA[x]==stampA&&markB[x]!=stampB){markB[x]=stampB;if(++common>2)return false;}}}return common==2;}
    bool flipOK(int keep,int rem,const Vec3&p,double mindot){static vector<int>touch;touch.clear();for(int fid:inc[keep])if(F[fid].active&&contains(fid,keep))touch.push_back(fid);for(int fid:inc[rem])if(F[fid].active&&contains(fid,rem))touch.push_back(fid);sort(touch.begin(),touch.end());touch.erase(unique(touch.begin(),touch.end()),touch.end());for(int fid:touch){Face&f=F[fid];if(!f.active)continue;bool hk=contains(fid,keep),hr=contains(fid,rem);if(hk&&hr)continue;Vec3 a=P[f.v[0]],b=P[f.v[1]],c=P[f.v[2]];if(f.v[0]==keep||f.v[0]==rem)a=p;if(f.v[1]==keep||f.v[1]==rem)b=p;if(f.v[2]==keep||f.v[2]==rem)c=p;Vec3 cr=cross3(b-a,c-a);double a2=norm2(cr);if(!(a2>minArea2)||!isfinite(a2))return false;Vec3 n=cr/sqrt(a2);if(dot3(n,f.n)<mindot)return false;}return true;}
    void collectNbr(int u,vector<int>&out){out.clear();if(!alive[u])return;maybeCleanup(u);if(++stampA>1000000000){fill(markA.begin(),markA.end(),0);stampA=1;}for(int fid:inc[u]){Face&f=F[fid];if(!f.active||!contains(fid,u))continue;for(int k=0;k<3;++k){int w=f.v[k];if(w!=u&&alive[w]&&markA[w]!=stampA){markA[w]=stampA;out.push_back(w);}}}}
    void collapse(int keep,int rem,const Vec3&p){cleanup(keep);cleanup(rem);vector<int>rf=inc[rem];P[keep]=p;for(int fid:rf){Face&f=F[fid];if(!f.active||!contains(fid,rem))continue;bool hasK=contains(fid,keep);if(hasK){f.active=0;--activeF;}else{for(int k=0;k<3;++k)if(f.v[k]==rem)f.v[k]=keep;if(f.v[0]==f.v[1]||f.v[0]==f.v[2]||f.v[1]==f.v[2]){f.active=0;--activeF;}else{f.n=fNormal(f);inc[keep].push_back(fid);}}}for(int fid:inc[keep])if(F[fid].active&&contains(fid,keep))F[fid].n=fNormal(F[fid]);Q[keep].add(Q[rem]);for(int k=0;k<3;++k){bbMin[keep][k]=min(bbMin[keep][k],bbMin[rem][k]);bbMax[keep][k]=max(bbMax[keep][k],bbMax[rem][k]);}locked[keep]=locked[keep]||locked[rem];alive[rem]=0;++ver[keep];++ver[rem];--activeV;inc[rem].clear();cleanup(keep);static vector<int>nbr;collectNbr(keep,nbr);for(int w:nbr)pushEdge(keep,w);}
    Snapshot snapshot(){Snapshot s;s.ratio=(double)activeV/(double)max(1,N);tmp.assign(N,-1);s.V.reserve(activeV);for(int i=0;i<N;++i)if(alive[i]){tmp[i]=(int)s.V.size();s.V.push_back(P[i]);}s.F.reserve(activeF);for(const Face&f:F)if(f.active){int a=f.v[0],b=f.v[1],c=f.v[2];if(a!=b&&a!=c&&b!=c&&a>=0&&b>=0&&c>=0&&a<N&&b<N&&c<N&&tmp[a]>=0&&tmp[b]>=0&&tmp[c]>=0)s.F.push_back({tmp[a],tmp[b],tmp[c]});}return s;}
    void run(){vector<double>ratios={.55,.40,.30,.24,.20,.17,.15,.13,.115,.105,.095,.085,.078,.070,.064,.058,.052,.047,.042,.038};vector<int>cp;int last=-1;for(double r:ratios){int c=max(4,(int)ceil(N*r));if(c<N&&c!=last){cp.push_back(c);last=c;}}sort(cp.begin(),cp.end(),greater<int>());int idx=0,target=cp.empty()?max(4,N/10):cp.back();long long pops=0;while(activeV>target&&!pq.empty()&&timeOK(16.9)){Cand c=pq.top();pq.pop();if((++pops&4095)==0&&!timeOK(16.9))break;int u=c.u,v=c.v;if(u==v||!alive[u]||!alive[v])continue;if(ver[u]!=c.vu||ver[v]!=c.vv)continue;Eval e=evalEdge(u,v);if(!e.ok)continue;if(!linkOK(u,v))continue;int keep=u,rem=v;if(inc[v].size()>inc[u].size()){keep=v;rem=u;}if(locked[rem]&&!locked[keep])swap(keep,rem);if(locked[keep]&&locked[rem])continue;if(!flipOK(keep,rem,e.pos,-0.04))continue;collapse(keep,rem,e.pos);while(idx<(int)cp.size()&&activeV<=cp[idx]){snaps.push_back(snapshot());++idx;}}
        if(snaps.empty()||snaps.back().V.size()!=(size_t)activeV)snaps.push_back(snapshot());}
};

static void orient_face(vector<Vec3>&V,Tri&f,const Vec3&center){Vec3 cr=cross3(V[f.b]-V[f.a],V[f.c]-V[f.a]);Vec3 ct=(V[f.a]+V[f.b]+V[f.c])/3.0;if(dot3(cr,ct-center)<0) swap(f.b,f.c);} 

static bool make_sphere_mesh(Vec3 c,double r,int lat,int lon,vector<Vec3>&V,vector<Tri>&F){if(lat<4||lon<8)return false;V.clear();F.clear();V.push_back({c.x,c.y,c.z+r});for(int i=1;i<lat;++i){double ph=acos(-1.0)*i/lat,sp=sin(ph),cp=cos(ph);for(int j=0;j<lon;++j){double th=2*acos(-1.0)*j/lon;V.push_back({c.x+r*sp*cos(th),c.y+r*sp*sin(th),c.z+r*cp});}}int bottom=V.size();V.push_back({c.x,c.y,c.z-r});auto id=[&](int ring,int j){return 1+(ring-1)*lon+(j%lon+lon)%lon;};auto add=[&](int a,int b,int cc){Tri t{a,b,cc};orient_face(V,t,c);F.push_back(t);};for(int j=0;j<lon;++j)add(0,id(1,j),id(1,j+1));for(int i=1;i<lat-1;++i)for(int j=0;j<lon;++j){int a=id(i,j),b=id(i,j+1),cc=id(i+1,j),d=id(i+1,j+1);add(a,cc,b);add(b,cc,d);}for(int j=0;j<lon;++j)add(bottom,id(lat-1,j+1),id(lat-1,j));return validate_mesh(V,F);} 

static bool try_sphere(const RenderMaps&orig,int R,vector<Vec3>&bestV,vector<Tri>&bestF){if(ORIGN<600||GTIMER.elapsed()>2.5)return false;Vec3 mn=ORIGV[0],mx=ORIGV[0];for(auto&p:ORIGV){mn.x=min(mn.x,p.x);mn.y=min(mn.y,p.y);mn.z=min(mn.z,p.z);mx.x=max(mx.x,p.x);mx.y=max(mx.y,p.y);mx.z=max(mx.z,p.z);}double ex=mx.x-mn.x,ey=mx.y-mn.y,ez=mx.z-mn.z,hi=max(ex,max(ey,ez)),lo=min(ex,min(ey,ez));if(!(hi>1e-12)||lo<.965*hi)return false;Vec3 c=(mn+mx)*.5;int stride=max(1,ORIGN/200000);double sr=0,sr2=0;int cnt=0;for(int i=0;i<ORIGN;i+=stride){double r=norm3(ORIGV[i]-c);sr+=r;sr2+=r*r;++cnt;}double meanr=sr/max(1,cnt);double r=.5*hi;if(!(r>1e-12))return false;double rms=sqrt(max(0.0,sr2/max(1,cnt)-meanr*meanr))/r,maxdev=0;for(int i=0;i<ORIGN;i+=stride)maxdev=max(maxdev,fabs(norm3(ORIGV[i]-c)-r)/r);if(rms>.012||maxdev>.055)return false;vector<tuple<int,int,double>> trials={{14,28,.925},{18,36,.92},{22,44,.925},{26,52,.935}};bool ok=false;double bestRatio=1;for(auto [lat,lon,thr]:trials){if(GTIMER.elapsed()>5.0)break;vector<Vec3>V;vector<Tri>F;if(!make_sphere_mesh(c,r,lat,lon,V,F))continue;if(V.size()>=bestV.size()||V.size()>=ORIGV.size())continue;double sc=proxy_score(orig,V,F,R);if(sc>=thr&&((double)V.size()/ORIGN<bestRatio)){bestRatio=(double)V.size()/ORIGN;bestV.swap(V);bestF.swap(F);ok=true;}}return ok;}

static bool make_axis_revolve(int axis,double t0,double t1,double cu,double cv,double r0,double r1,int sides,vector<Vec3>&V,vector<Tri>&F){V.clear();F.clear();if(sides<12)return false;auto make=[&](double t,double r,int j){double th=2*acos(-1.0)*j/sides,u=cu+r*cos(th),v=cv+r*sin(th);if(axis==0)return Vec3{t,u,v};if(axis==1)return Vec3{u,t,v};return Vec3{u,v,t};};auto center=[&](double t){if(axis==0)return Vec3{t,cu,cv};if(axis==1)return Vec3{cu,t,cv};return Vec3{cu,cv,t};};Vec3 mid=center(.5*(t0+t1));bool cone0=r0<max(r0,r1)*1e-5, cone1=r1<max(r0,r1)*1e-5;if(cone0&&cone1)return false;int c0=-1,c1=-1,s0=-1,s1=-1;if(!cone0){c0=V.size();V.push_back(center(t0));s0=V.size();for(int j=0;j<sides;++j)V.push_back(make(t0,r0,j));}else{c0=V.size();V.push_back(center(t0));}
    if(!cone1){c1=V.size();V.push_back(center(t1));s1=V.size();for(int j=0;j<sides;++j)V.push_back(make(t1,r1,j));}else{c1=V.size();V.push_back(center(t1));}
    auto ring=[&](int s,int j){return s+(j%sides+sides)%sides;};auto add=[&](int a,int b,int c){Tri t{a,b,c};orient_face(V,t,mid);F.push_back(t);};
    if(!cone0&&!cone1){for(int j=0;j<sides;++j){int k=(j+1)%sides;add(ring(s0,j),ring(s0,k),ring(s1,j));add(ring(s0,k),ring(s1,k),ring(s1,j));add(c0,ring(s0,j),ring(s0,k));add(c1,ring(s1,k),ring(s1,j));}}
    else if(cone0){for(int j=0;j<sides;++j){int k=(j+1)%sides;add(c0,ring(s1,j),ring(s1,k));add(c1,ring(s1,k),ring(s1,j));}}
    else {for(int j=0;j<sides;++j){int k=(j+1)%sides;add(ring(s0,j),ring(s0,k),c1);add(c0,ring(s0,j),ring(s0,k));}}
    return validate_mesh(V,F);
}

static bool try_revolve(const RenderMaps&orig,int R,vector<Vec3>&bestV,vector<Tri>&bestF){if(ORIGN<800||GTIMER.elapsed()>5.5)return false;bool any=false;double br=1;for(int axis=0;axis<3;++axis){double lo=1e100,hi=-1e100,umin=1e100,umax=-1e100,vmin=1e100,vmax=-1e100;auto split=[&](const Vec3&p,double&t,double&u,double&v){if(axis==0){t=p.x;u=p.y;v=p.z;}else if(axis==1){t=p.y;u=p.x;v=p.z;}else{t=p.z;u=p.x;v=p.y;}};for(auto&p:ORIGV){double t,u,v;split(p,t,u,v);lo=min(lo,t);hi=max(hi,t);umin=min(umin,u);umax=max(umax,u);vmin=min(vmin,v);vmax=max(vmax,v);}double cu=.5*(umin+umax),cv=.5*(vmin+vmax),len=hi-lo;if(!(len>1e-12))continue;int stride=max(1,ORIGN/200000);double sr=0,srt=0,st=0,stt=0;int cnt=0;for(int i=0;i<ORIGN;i+=stride){double t,u,v;split(ORIGV[i],t,u,v);double r=hypot(u-cu,v-cv);if(r<1e-6*GLOBAL_DIAG)continue;sr+=r;st+=t;srt+=r*t;stt+=t*t;++cnt;}if(cnt<200)continue;double den=cnt*stt-st*st;if(fabs(den)<1e-18)continue;double a=(cnt*srt-st*sr)/den,b=(sr-a*st)/cnt;double r0=max(0.0,a*lo+b),r1=max(0.0,a*hi+b),mr=max(r0,r1);if(!(mr>1e-10))continue;double ss=0,ma=0;int chk=0;for(int i=0;i<ORIGN;i+=stride){double t,u,v;split(ORIGV[i],t,u,v);double r=hypot(u-cu,v-cv);double pred=max(0.0,a*t+b);double e=fabs(r-pred)/mr;ss+=e*e;ma=max(ma,e);++chk;}double rms=sqrt(ss/max(1,chk));if(rms>.035||ma>.16)continue;for(int sides: {24,32,40,56,72}){if(GTIMER.elapsed()>8.0)break;vector<Vec3>V;vector<Tri>F;if(!make_axis_revolve(axis,lo,hi,cu,cv,r0,r1,sides,V,F))continue;if(V.size()>=bestV.size()||V.size()>=ORIGV.size())continue;double sc=proxy_score(orig,V,F,R);double thr=(sides<=32?.94:.93);if(sc>=thr&&((double)V.size()/ORIGN<br)){br=(double)V.size()/ORIGN;bestV.swap(V);bestF.swap(F);any=true;}}}return any;}

static void output_mesh(const vector<Vec3>&V,const vector<Tri>&F){static char obuf[1<<20];setvbuf(stdout,obuf,_IOFBF,sizeof(obuf));printf("%d %d\n",(int)V.size(),(int)F.size());bool hi=(int)V.size()*2<=ORIGN;for(auto&p:V){if(hi)printf("v %.15g %.15g %.15g\n",p.x,p.y,p.z);else printf("v %.10g %.10g %.10g\n",p.x,p.y,p.z);}for(auto&t:F)printf("f %d %d %d\n",t.a+1,t.b+1,t.c+1);}

int main(){GTIMER=Timer();FastScanner fs;int N,M;if(!fs.readInt(N))return 0;fs.readInt(M);ORIGN=N;ORIGM=M;ORIGV.resize(N);ORIGF.resize(M);Vec3 mn{1e100,1e100,1e100},mx{-1e100,-1e100,-1e100};for(int i=0;i<N;++i){fs.prefix();fs.readDouble(ORIGV[i].x);fs.readDouble(ORIGV[i].y);fs.readDouble(ORIGV[i].z);mn.x=min(mn.x,ORIGV[i].x);mn.y=min(mn.y,ORIGV[i].y);mn.z=min(mn.z,ORIGV[i].z);mx.x=max(mx.x,ORIGV[i].x);mx.y=max(mx.y,ORIGV[i].y);mx.z=max(mx.z,ORIGV[i].z);}for(int i=0;i<M;++i){fs.prefix();int a,b,c;fs.readInt(a);fs.readInt(b);fs.readInt(c);ORIGF[i]={a-1,b-1,c-1};}GLOBAL_DIAG=max(1e-12,norm3(mx-mn));int R=(N<8000?512:(N<80000?384:256));vector<Vec3>bestV=ORIGV;vector<Tri>bestF=ORIGF;double bestRatio=1.0,bestSSIM=1.0;RenderMaps origMaps;render_mesh(ORIGV,ORIGF,origMaps,R);
    if(try_sphere(origMaps,R,bestV,bestF)){bestRatio=(double)bestV.size()/N;bestSSIM=proxy_score(origMaps,bestV,bestF,R);}    
    if(false&&GTIMER.elapsed()<8.0){vector<Vec3>V=bestV;vector<Tri>F=bestF;if(try_revolve(origMaps,R,V,F)&&validate_mesh(V,F)){double sc=proxy_score(origMaps,V,F,R);if(sc>=.925&&(double)V.size()/N<bestRatio){bestV.swap(V);bestF.swap(F);bestRatio=(double)bestV.size()/N;bestSSIM=sc;}}}
    if(GTIMER.elapsed()<16.8){Simplifier sim(ORIGV,ORIGF);sim.run();vector<int>order(sim.snaps.size());iota(order.begin(),order.end(),0);sort(order.begin(),order.end(),[&](int a,int b){return sim.snaps[a].ratio<sim.snaps[b].ratio;});double guard=(R>=512?.908:(R>=384?.914:.922));int chosen=-1;for(int idx:order){if(GTIMER.elapsed()>19.0)break;Snapshot&s=sim.snaps[idx];if(s.V.size()>=bestV.size()||s.V.empty()||s.F.empty())continue;if(!validate_mesh(s.V,s.F))continue;s.ssim=proxy_score(origMaps,s.V,s.F,R);if(s.ssim>=guard){chosen=idx;break;}}if(chosen>=0){bestV.swap(sim.snaps[chosen].V);bestF.swap(sim.snaps[chosen].F);bestRatio=(double)bestV.size()/N;bestSSIM=sim.snaps[chosen].ssim;}else{double score=bestSSIM-.08*bestRatio;for(int idx:order){if(sim.snaps[idx].ssim<-1||sim.snaps[idx].V.empty())continue;double sc=sim.snaps[idx].ssim-.08*sim.snaps[idx].ratio;if(sc>score&&sim.snaps[idx].ssim>.89&&validate_mesh(sim.snaps[idx].V,sim.snaps[idx].F)){score=sc;bestV=sim.snaps[idx].V;bestF=sim.snaps[idx].F;}}}}
    if(bestV.empty()||bestF.empty()||bestV.size()>=ORIGV.size()||!validate_mesh(bestV,bestF)){output_mesh(ORIGV,ORIGF);return 0;}output_mesh(bestV,bestF);return 0;}
