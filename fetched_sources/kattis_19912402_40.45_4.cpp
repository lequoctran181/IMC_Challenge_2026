#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cstdint>
#include <queue>
#include <unordered_set>
#include <vector>
using namespace std;

struct Vec3{double x,y,z;};
static inline Vec3 operator+(const Vec3&a,const Vec3&b){return {a.x+b.x,a.y+b.y,a.z+b.z};}
static inline Vec3 operator-(const Vec3&a,const Vec3&b){return {a.x-b.x,a.y-b.y,a.z-b.z};}
static inline Vec3 operator*(const Vec3&a,double s){return {a.x*s,a.y*s,a.z*s};}
static inline double dot3(const Vec3&a,const Vec3&b){return a.x*b.x+a.y*b.y+a.z*b.z;}
static inline Vec3 cross3(const Vec3&a,const Vec3&b){return {a.y*b.z-a.z*b.y,a.z*b.x-a.x*b.z,a.x*b.y-a.y*b.x};}
static inline double norm2(const Vec3&a){return dot3(a,a);} 
static inline double dist2(const Vec3&a,const Vec3&b){return norm2(a-b);} 
static inline double norm3(const Vec3&a){return sqrt(norm2(a));}

struct Face{int v[3]; unsigned char active;};
struct In{
    vector<char>b; char*p;
    In(){char tmp[1<<16]; size_t n; while((n=fread(tmp,1,sizeof(tmp),stdin))>0)b.insert(b.end(),tmp,tmp+n); b.push_back(0); p=b.data();}
    inline void ws(){while(*p==' '||*p=='\n'||*p=='\r'||*p=='\t')++p;}
    int ni(){ws(); int s=1; if(*p=='-')s=-1,++p; int x=0; while(*p>='0'&&*p<='9')x=x*10+*p++-'0'; return x*s;}
    double nd(){ws(); char*e; double x=strtod(p,&e); p=e; return x;}
    char nc(){ws(); return *p++;}
};
struct Quadric{
    double q[10];
    Quadric(){memset(q,0,sizeof(q));}
    inline void addPlane(double a,double b,double c,double d,double w=1.0){q[0]+=w*a*a;q[1]+=w*a*b;q[2]+=w*a*c;q[3]+=w*a*d;q[4]+=w*b*b;q[5]+=w*b*c;q[6]+=w*b*d;q[7]+=w*c*c;q[8]+=w*c*d;q[9]+=w*d*d;}
    inline void add(const Quadric&o){for(int i=0;i<10;i++)q[i]+=o.q[i];}
    inline double eval(const Vec3&p)const{double x=p.x,y=p.y,z=p.z;return q[0]*x*x+2*q[1]*x*y+2*q[2]*x*z+2*q[3]*x+q[4]*y*y+2*q[5]*y*z+2*q[6]*y+q[7]*z*z+2*q[8]*z+q[9];}
};
struct EdgeFace{uint64_t k; int f; bool operator<(const EdgeFace&o)const{return k<o.k||(k==o.k&&f<o.f);} };
struct Node{double c; int a,b,va,vb; bool operator<(const Node&o)const{return c>o.c;}};

static int N,M,activeV,activeF,targetV;
static vector<Vec3>P,Orig,origP0;
static vector<Face>F,origF0;
static vector<vector<int>>Inc;
static vector<Quadric>Q;
static vector<unsigned char>VAlive,Locked;
static vector<int>Head,Tail,Nxt,Csz,Ver,Mark;
static unordered_set<uint64_t> FaceSet;
static priority_queue<Node> PQ;
static double diagLen=1,haus2=1,areaEps=1e-30,minCos=0.48;
static chrono::steady_clock::time_point T0;

static inline bool aliveFace(int f){return F[f].active;}
static inline bool hasv(const Face&f,int v){return f.v[0]==v||f.v[1]==v||f.v[2]==v;}
static inline uint64_t ekey(int a,int b){if(a>b)swap(a,b);return (uint64_t)(uint32_t)a<<32 | (uint32_t)b;}
static inline uint64_t fkey(int a,int b,int c){if(a>b)swap(a,b); if(b>c)swap(b,c); if(a>b)swap(a,b); return ((uint64_t)a<<42)|((uint64_t)b<<21)|(uint64_t)c;}
static inline uint64_t fkey(const Face&f){return fkey(f.v[0],f.v[1],f.v[2]);}
static inline double elapsed(){return chrono::duration<double>(chrono::steady_clock::now()-T0).count();}
static inline bool timeok(){return elapsed()<19.45;}
static inline Vec3 fnormIdx(int a,int b,int c){return cross3(P[b]-P[a],P[c]-P[a]);}
static inline Vec3 fnorm(const Face&f){return fnormIdx(f.v[0],f.v[1],f.v[2]);}

static void readInput(){
    In in; N=in.ni(); M=in.ni(); P.resize(N); Orig.resize(N); origP0.resize(N);
    Vec3 mn{1e100,1e100,1e100}, mx{-1e100,-1e100,-1e100};
    for(int i=0;i<N;i++){in.nc(); P[i].x=in.nd(); P[i].y=in.nd(); P[i].z=in.nd(); Orig[i]=origP0[i]=P[i]; mn.x=min(mn.x,P[i].x);mn.y=min(mn.y,P[i].y);mn.z=min(mn.z,P[i].z); mx.x=max(mx.x,P[i].x);mx.y=max(mx.y,P[i].y);mx.z=max(mx.z,P[i].z);} 
    F.resize(M); origF0.resize(M); vector<int>deg(N,0);
    for(int i=0;i<M;i++){in.nc(); int a=in.ni()-1,b=in.ni()-1,c=in.ni()-1; F[i]={{a,b,c},1}; origF0[i]=F[i]; deg[a]++;deg[b]++;deg[c]++;}
    diagLen=norm3(mx-mn); if(!(diagLen>0))diagLen=1; haus2=(0.0492*diagLen)*(0.0492*diagLen); areaEps=1e-28*diagLen*diagLen*diagLen*diagLen;
    Inc.resize(N); for(int i=0;i<N;i++)Inc[i].reserve(deg[i]+8); for(int i=0;i<M;i++)for(int k=0;k<3;k++)Inc[F[i].v[k]].push_back(i);
    Q.assign(N,Quadric()); VAlive.assign(N,1); Locked.assign(N,0); Head.resize(N); Tail.resize(N); Nxt.assign(N,-1); Csz.assign(N,1); Ver.assign(N,0); Mark.assign(N,0);
    for(int i=0;i<N;i++)Head[i]=Tail[i]=i;
    activeV=N; activeF=M;
    int ixmn=0,iymn=0,izmn=0,ixmx=0,iymx=0,izmx=0;
    for(int i=1;i<N;i++){if(P[i].x<P[ixmn].x)ixmn=i;if(P[i].y<P[iymn].y)iymn=i;if(P[i].z<P[izmn].z)izmn=i;if(P[i].x>P[ixmx].x)ixmx=i;if(P[i].y>P[iymx].y)iymx=i;if(P[i].z>P[izmx].z)izmx=i;}
    Locked[ixmn]=Locked[iymn]=Locked[izmn]=Locked[ixmx]=Locked[iymx]=Locked[izmx]=1;
    FaceSet.reserve((size_t)M*2+1024); for(int i=0;i<M;i++)FaceSet.insert(fkey(F[i]));
}

static bool normalize(Vec3&v){double n=norm3(v); if(!(n>1e-300))return false; v=v*(1.0/n); return true;}
static bool solveOpt(const Quadric&q,Vec3&out){
    double a00=q.q[0],a01=q.q[1],a02=q.q[2],a11=q.q[4],a12=q.q[5],a22=q.q[7];
    double b0=-q.q[3],b1=-q.q[6],b2=-q.q[8];
    double det=a00*(a11*a22-a12*a12)-a01*(a01*a22-a12*a02)+a02*(a01*a12-a11*a02);
    if(fabs(det)<1e-14)return false;
    double dx=b0*(a11*a22-a12*a12)-a01*(b1*a22-a12*b2)+a02*(b1*a12-a11*b2);
    double dy=a00*(b1*a22-a12*b2)-b0*(a01*a22-a12*a02)+a02*(a01*b2-b1*a02);
    double dz=a00*(a11*b2-b1*a12)-a01*(a01*b2-b1*a02)+b0*(a01*a12-a11*a02);
    out={dx/det,dy/det,dz/det}; return isfinite(out.x)&&isfinite(out.y)&&isfinite(out.z);
}

static bool collectEdge(int a,int b,int ef[2],int op[2]){
    int c=0; const vector<int>&S=Inc[a].size()<Inc[b].size()?Inc[a]:Inc[b];
    for(int id:S){if(!F[id].active)continue; const Face&f=F[id]; if(!hasv(f,a)||!hasv(f,b))continue; if(c>=2)return false; ef[c]=id; for(int k=0;k<3;k++){int x=f.v[k]; if(x!=a&&x!=b){op[c]=x; break;}} c++;}
    return c==2&&op[0]!=op[1]&&op[0]>=0&&op[1]>=0;
}
static bool linkOK(int a,int b){
    int ef[2]={-1,-1},op[2]={-1,-1}; if(!collectEdge(a,b,ef,op))return false;
    static int tok=1; if(tok>2000000000){fill(Mark.begin(),Mark.end(),0); tok=1;} int t=tok++;
    for(int id:Inc[a])if(F[id].active&&hasv(F[id],a)){const Face&f=F[id]; for(int k=0;k<3;k++){int x=f.v[k]; if(x!=a&&x!=b)Mark[x]=t;}}
    int cnt=0,g0=0,g1=0; int t2=tok++;
    for(int id:Inc[b])if(F[id].active&&hasv(F[id],b)){const Face&f=F[id]; for(int k=0;k<3;k++){int x=f.v[k]; if(x==a||x==b)continue; if(Mark[x]==t){ if(x!=op[0]&&x!=op[1])return false; if(Mark[x]!=t2){Mark[x]=t2;cnt++; if(x==op[0])g0=1; if(x==op[1])g1=1;} }}}
    return cnt==2&&g0&&g1;
}
static bool clusterOK(int src,int dst){for(int m=Head[src];m!=-1;m=Nxt[m])if(dist2(Orig[m],P[dst])>haus2)return false; return true;}
static bool canRemove(int src,int dst,double &cost){
    if(src==dst||!VAlive[src]||!VAlive[dst]||Locked[src])return false;
    if(dist2(P[src],P[dst])>haus2*1.00001)return false;
    if(!clusterOK(src,dst))return false;
    vector<uint64_t> nk; nk.reserve(Inc[src].size());
    for(int id:Inc[src]){
        if(!F[id].active||!hasv(F[id],src))continue;
        Face f=F[id]; bool hd=hasv(f,dst);
        if(hd)continue;
        for(int k=0;k<3;k++)if(f.v[k]==src)f.v[k]=dst;
        if(f.v[0]==f.v[1]||f.v[0]==f.v[2]||f.v[1]==f.v[2])return false;
        Vec3 oldn=fnorm(F[id]), newn=fnorm(f); double no=norm2(oldn), nn=norm2(newn); if(!(no>areaEps)||!(nn>areaEps))return false;
        if(dot3(oldn,newn)<=minCos*sqrt(no*nn))return false;
        uint64_t k=fkey(f); if(FaceSet.find(k)!=FaceSet.end())return false; nk.push_back(k);
    }
    sort(nk.begin(),nk.end()); for(size_t i=1;i<nk.size();i++)if(nk[i]==nk[i-1])return false;
    Quadric qq=Q[src]; qq.add(Q[dst]); double e=dist2(P[src],P[dst])/(diagLen*diagLen+1e-300);
    cost=qq.eval(P[dst])+0.02*Csz[src]+0.0001*Csz[dst]+1e-4*e;
    return true;
}
static void compact(int v){auto &a=Inc[v]; if(a.size()<96)return; size_t live=0; for(int id:a)if(F[id].active&&hasv(F[id],v))live++; if(live*3+32>=a.size())return; vector<int>b; b.reserve(live+8); for(int id:a)if(F[id].active&&hasv(F[id],v))b.push_back(id); a.swap(b);} 
static void pushEdge(int a,int b){
    if(a==b||a<0||b<0||a>=N||b>=N||!VAlive[a]||!VAlive[b]||(Locked[a]&&Locked[b]))return;
    if(dist2(P[a],P[b])>haus2*1.00001)return;
    Quadric qq=Q[a]; qq.add(Q[b]); double c=min(qq.eval(P[a]),qq.eval(P[b]))+1e-6*dist2(P[a],P[b])/(diagLen*diagLen+1e-300);
    PQ.push({c,a,b,Ver[a],Ver[b]});
}
static void mergeList(int src,int dst){if(Head[src]<0)return; Nxt[Tail[dst]]=Head[src]; Tail[dst]=Tail[src]; Csz[dst]+=Csz[src]; Head[src]=Tail[src]=-1; Csz[src]=0;}
static void collapseDo(int src,int dst){
    Q[dst].add(Q[src]);
    for(int id:Inc[src])if(F[id].active&&hasv(F[id],src))FaceSet.erase(fkey(F[id]));
    for(int id:Inc[src]){
        if(!F[id].active||!hasv(F[id],src))continue;
        if(hasv(F[id],dst)){F[id].active=0; activeF--; continue;}
        for(int k=0;k<3;k++)if(F[id].v[k]==src)F[id].v[k]=dst;
        Inc[dst].push_back(id); FaceSet.insert(fkey(F[id]));
    }
    VAlive[src]=0; activeV--; Ver[src]++; Ver[dst]++; mergeList(src,dst); Inc[src].clear(); compact(dst);
    for(int id:Inc[dst])if(F[id].active&&hasv(F[id],dst)){Face &f=F[id]; pushEdge(f.v[0],f.v[1]); pushEdge(f.v[1],f.v[2]); pushEdge(f.v[2],f.v[0]);}
}
static bool attempt(int a,int b){
    if(a==b||!VAlive[a]||!VAlive[b]||!linkOK(a,b))return false;
    double ca=1e100,cb=1e100; bool ra=canRemove(a,b,ca), rb=canRemove(b,a,cb);
    if(!ra&&!rb)return false;
    if(rb&&(!ra||cb<ca))collapseDo(b,a); else collapseDo(a,b);
    return true;
}

static double buildQuadricsAndQueue(){
    vector<Vec3> FN(M); vector<EdgeFace> E; E.reserve((size_t)M*3);
    for(int i=0;i<M;i++){
        Vec3 n=fnorm(F[i]); if(!normalize(n))n={0,0,0}; FN[i]=n; double d=-dot3(n,P[F[i].v[0]]);
        for(int k=0;k<3;k++)Q[F[i].v[k]].addPlane(n.x,n.y,n.z,d);
        E.push_back({ekey(F[i].v[0],F[i].v[1]),i}); E.push_back({ekey(F[i].v[1],F[i].v[2]),i}); E.push_back({ekey(F[i].v[2],F[i].v[0]),i});
    }
    sort(E.begin(),E.end()); long long ue=0,feat=0; double fc=cos(38.0*acos(-1.0)/180.0);
    for(size_t i=0;i<E.size();){size_t j=i+1; while(j<E.size()&&E[j].k==E[i].k)j++; ue++; if(j-i==2&&dot3(FN[E[i].f],FN[E[i+1].f])<fc)feat++; int a=(int)(E[i].k>>32),b=(int)(E[i].k&0xffffffffu); pushEdge(a,b); i=j;}
    return ue?double(feat)/double(ue):0.0;
}
static void chooseTarget(double fr){
    double r=0.142+min(0.105,fr*0.38);
    if(N<800)r=0.62; else if(N<2500)r=max(r,0.32); else if(N<8000)r=max(r,0.23); else if(N<30000)r=max(r,0.18);
    if(N>200000&&fr<0.04)r=min(r,0.125);
    r=max(0.115,min(0.30,r)); targetV=max(8,(int)ceil(N*r));
}
static void simplify(){
    double fr=buildQuadricsAndQueue(); chooseTarget(fr); long long pops=0;
    while(activeV>targetV&&!PQ.empty()){
        if((++pops&4095)==0&&!timeok())break;
        Node n=PQ.top(); PQ.pop(); int a=n.a,b=n.b;
        if(a==b||!VAlive[a]||!VAlive[b])continue;
        if(n.va!=Ver[a]||n.vb!=Ver[b]){pushEdge(a,b); continue;}
        attempt(a,b);
    }
}

static bool bboxSame(const vector<Vec3>&V){
    if(V.empty())return false; Vec3 a=V[0],b=V[0],c=origP0[0],d=origP0[0];
    for(auto&p:V){a.x=min(a.x,p.x);a.y=min(a.y,p.y);a.z=min(a.z,p.z);b.x=max(b.x,p.x);b.y=max(b.y,p.y);b.z=max(b.z,p.z);} for(auto&p:origP0){c.x=min(c.x,p.x);c.y=min(c.y,p.y);c.z=min(c.z,p.z);d.x=max(d.x,p.x);d.y=max(d.y,p.y);d.z=max(d.z,p.z);} double e=diagLen*1e-11+1e-12; return fabs(a.x-c.x)<e&&fabs(a.y-c.y)<e&&fabs(a.z-c.z)<e&&fabs(b.x-d.x)<e&&fabs(b.y-d.y)<e&&fabs(b.z-d.z)<e;
}

static double orientSignOrig(const Vec3&c){
    double s=0; int st=max(1,M/120000); for(int i=0;i<M;i+=st){auto &f=origF0[i]; Vec3 a=origP0[f.v[0]],b=origP0[f.v[1]],d=origP0[f.v[2]]; Vec3 cr=cross3(b-a,d-a); Vec3 ce=(a+b+d)*(1.0/3.0); s+=dot3(cr,ce-c);} return s>=0?1.0:-1.0;
}
static bool genOK(const vector<Vec3>&V,const vector<array<int,3>>&G){
    if(V.empty()||G.empty()||V.size()>origP0.size()||!bboxSame(V))return false; double eps=1e-28*diagLen*diagLen*diagLen*diagLen; unordered_set<uint64_t>S; S.reserve(G.size()*2+16);
    for(auto&t:G){int a=t[0],b=t[1],c=t[2]; if(a<0||b<0||c<0||a>=(int)V.size()||b>=(int)V.size()||c>=(int)V.size()||a==b||a==c||b==c)return false; if(norm2(cross3(V[b]-V[a],V[c]-V[a]))<=eps)return false; if(!S.insert(fkey(a,b,c)).second)return false;} return true;
}
static bool tryBox(vector<Vec3>&V,vector<array<int,3>>&G){
    if(N<500)return false; Vec3 mn=origP0[0],mx=origP0[0]; for(auto&p:origP0){mn.x=min(mn.x,p.x);mn.y=min(mn.y,p.y);mn.z=min(mn.z,p.z);mx.x=max(mx.x,p.x);mx.y=max(mx.y,p.y);mx.z=max(mx.z,p.z);} double ex=mx.x-mn.x,ey=mx.y-mn.y,ez=mx.z-mn.z; if(ex<=0||ey<=0||ez<=0)return false; int st=max(1,N/180000),bad=0,samp=0; double tol=diagLen*(N<5000?0.0015:0.001);
    for(int i=0;i<N;i+=st){auto&p=origP0[i]; double d=min({fabs(p.x-mn.x),fabs(p.x-mx.x),fabs(p.y-mn.y),fabs(p.y-mx.y),fabs(p.z-mn.z),fabs(p.z-mx.z)}); if(d>tol&&++bad>2+N/st/400)return false; samp++;}
    if(samp<200)return false; V={{mn.x,mn.y,mn.z},{mx.x,mn.y,mn.z},{mx.x,mx.y,mn.z},{mn.x,mx.y,mn.z},{mn.x,mn.y,mx.z},{mx.x,mn.y,mx.z},{mx.x,mx.y,mx.z},{mn.x,mx.y,mx.z}}; int t[12][3]={{0,2,1},{0,3,2},{4,5,6},{4,6,7},{0,1,5},{0,5,4},{1,2,6},{1,6,5},{2,3,7},{2,7,6},{3,0,4},{3,4,7}}; Vec3 c=(mn+mx)*0.5; double sg=orientSignOrig(c); G.clear(); for(auto &q:t){array<int,3>a={q[0],q[1],q[2]}; Vec3 cr=cross3(V[a[1]]-V[a[0]],V[a[2]]-V[a[0]]); Vec3 ce=(V[a[0]]+V[a[1]]+V[a[2]])*(1.0/3.0); if(dot3(cr,ce-c)*sg<0)swap(a[1],a[2]); G.push_back(a);} return genOK(V,G);
}
static bool tryEllipsoid(vector<Vec3>&V,vector<array<int,3>>&G){
    if(N<900)return false; Vec3 mn=origP0[0],mx=origP0[0]; for(auto&p:origP0){mn.x=min(mn.x,p.x);mn.y=min(mn.y,p.y);mn.z=min(mn.z,p.z);mx.x=max(mx.x,p.x);mx.y=max(mx.y,p.y);mx.z=max(mx.z,p.z);} Vec3 c=(mn+mx)*0.5; double rx=(mx.x-mn.x)*0.5,ry=(mx.y-mn.y)*0.5,rz=(mx.z-mn.z)*0.5; double hi=max(rx,max(ry,rz)),lo=min(rx,min(ry,rz)); if(!(lo>1e-12)||lo<hi*0.72)return false;
    int st=max(1,N/220000),samp=0; double ss=0,sm=0,mxerr=0; for(int i=0;i<N;i+=st){auto&p=origP0[i]; double x=(p.x-c.x)/rx,y=(p.y-c.y)/ry,z=(p.z-c.z)/rz; double e=fabs(sqrt(max(0.0,x*x+y*y+z*z))-1.0); ss+=e*e; sm+=e; mxerr=max(mxerr,e); samp++;}
    if(samp<200)return false; double rms=sqrt(ss/samp), mean=sm/samp; double rl=N<5000?0.0065:0.0050, ml=N<5000?0.028:0.022; if(rms>rl||mean>rl*0.72||mxerr>ml)return false;
    int lat,lon; if(N<3000)lat=14,lon=28; else if(N<9000)lat=16,lon=32; else if(N<40000)lat=18,lon=36; else if(N<120000)lat=20,lon=40; else lat=22,lon=44; if(2+(lat-1)*lon>=N*92/100)return false;
    V.clear();G.clear(); V.reserve(2+(lat-1)*lon); G.reserve(2*lat*lon); const double pi=acos(-1.0); V.push_back({c.x,c.y,c.z+rz}); for(int i=1;i<lat;i++){double th=pi*i/lat,stt=sin(th),ct=cos(th); for(int j=0;j<lon;j++){double ph=2*pi*j/lon; V.push_back({c.x+rx*stt*cos(ph),c.y+ry*stt*sin(ph),c.z+rz*ct});}} int bot=V.size(); V.push_back({c.x,c.y,c.z-rz}); auto id=[&](int r,int j){return 1+(r-1)*lon+((j%lon+lon)%lon);}; double sg=orientSignOrig(c); auto add=[&](int a,int b,int d){array<int,3>f={a,b,d}; Vec3 cr=cross3(V[f[1]]-V[f[0]],V[f[2]]-V[f[0]]); Vec3 ce=(V[f[0]]+V[f[1]]+V[f[2]])*(1.0/3.0); if(dot3(cr,ce-c)*sg<0)swap(f[1],f[2]); G.push_back(f);}; for(int j=0;j<lon;j++)add(0,id(1,j+1),id(1,j)); for(int r=1;r<lat-1;r++)for(int j=0;j<lon;j++){int a=id(r,j),b=id(r,j+1),cc=id(r+1,j),d=id(r+1,j+1); add(a,b,cc); add(b,d,cc);} for(int j=0;j<lon;j++)add(bot,id(lat-1,j),id(lat-1,j+1)); return genOK(V,G);
}

static bool detectWidth(int &W){
    if(N<200000||M<300000)return false; int cand[24]={0},cnt[24]={0}; auto add=[&](int d){for(int q=-1;q<=1;q++){int x=d+q;if(x<8||x>N/8||N%x)continue;int i=0;for(;i<24&&cand[i]&&cand[i]!=x;i++);if(i<24){cand[i]=x;cnt[i]++;}}}; int st=max(1,M/60000); for(int i=0;i<M;i+=st){add(abs(F[i].v[0]-F[i].v[1]));add(abs(F[i].v[1]-F[i].v[2]));add(abs(F[i].v[2]-F[i].v[0]));}
    for(int it=0;it<24;it++){int best=-1;for(int i=0;i<24;i++)if(cnt[i]>0&&(best<0||cnt[i]>cnt[best]))best=i; if(best<0)break; int w=cand[best]; cnt[best]=-1; int H=N/w; if(H<16)continue; long long expect=2LL*(H-1)*(w-1); if(llabs((long long)M-expect)*1000>expect*3)continue; int bad=0,step=max(1,M/1000); for(int i=0;i<M;i+=step){int s[3]={F[i].v[0],F[i].v[1],F[i].v[2]}; sort(s,s+3); int x=s[0],y=s[1],z=s[2],r=x/w,cc=x%w; bool ok=false; if(r<H-1&&cc<w-1&&y==x+1&&(z==x+w||z==x+w+1))ok=true; if(r<H-1&&cc>0&&y==x+w-1&&z==x+w)ok=true; if(!ok&&++bad>3)break; } if(bad<=3){W=w;return true;}}
    return false;
}
static bool tryOrderedGrid(vector<Vec3>&V,vector<array<int,3>>&G){
    int W=0; if(!detectWidth(W))return false; int H=N/W; double keep=0.58, s=sqrt(keep); int nw=max(32,(int)(W*s+0.5)), nh=max(32,(int)(H*s+0.5)); if(nw>=W||nh>=H||1LL*nw*nh>=N*92LL/100)return false;
    vector<int>rs(nh),cs(nw); for(int i=0;i<nh;i++)rs[i]=(long long)i*(H-1)/(nh-1); for(int j=0;j<nw;j++)cs[j]=(long long)j*(W-1)/(nw-1);
    V.reserve((size_t)nh*nw); for(int i=0;i<nh;i++)for(int j=0;j<nw;j++)V.push_back(origP0[rs[i]*W+cs[j]]); if(!bboxSame(V)){V.clear();return false;}
    G.reserve((size_t)(nh-1)*(nw-1)*2); auto add=[&](int a,int b,int c){G.push_back({a,b,c});};
    bool diag=false,orient=true; // infer from first cell if possible
    for(int i=0;i<min(M,256);i++){int s3[3]={F[i].v[0],F[i].v[1],F[i].v[2]};sort(s3,s3+3);int x=s3[0],y=s3[1],z=s3[2],r=x/W,cc=x%W;if(r<H-1&&cc<W-1){if(y==x+1&&z==x+W+1)diag=true;break;}}
    for(int i=0;i<nh-1;i++)for(int j=0;j<nw-1;j++){int a=i*nw+j,b=a+1,c=a+nw,d=c+1; if(!diag){if(orient)add(a,b,c),add(b,d,c);else add(a,c,b),add(b,c,d);}else{if(orient)add(a,b,d),add(a,d,c);else add(a,d,b),add(a,c,d);}}
    return true;
}

static bool validateCurrent(){
    unordered_set<uint64_t>S; S.reserve((size_t)activeF*2+16); vector<unsigned char>ref(N,0); int fc=0;
    for(int i=0;i<M;i++)if(F[i].active){Face &f=F[i]; for(int k=0;k<3;k++)if(f.v[k]<0||f.v[k]>=N||!VAlive[f.v[k]])return false; if(f.v[0]==f.v[1]||f.v[0]==f.v[2]||f.v[1]==f.v[2])return false; if(norm2(fnorm(f))<=areaEps)return false; uint64_t k=fkey(f); if(!S.insert(k).second)return false; ref[f.v[0]]=ref[f.v[1]]=ref[f.v[2]]=1; fc++;}
    for(int i=0;i<N;i++)if(VAlive[i]&&!ref[i])return false; return fc>0;
}
static void outputOriginal(){
    static char ob[1<<20]; setvbuf(stdout,ob,_IOFBF,sizeof(ob)); printf("%d %d\n",N,M); for(auto&p:origP0)printf("v %.17g %.17g %.17g\n",p.x,p.y,p.z); for(auto&f:origF0)printf("f %d %d %d\n",f.v[0]+1,f.v[1]+1,f.v[2]+1);
}
static void outputCurrent(){
    if(!validateCurrent()){outputOriginal();return;} vector<int>mp(N,-1); int n=0,m=0; for(int i=0;i<N;i++)if(VAlive[i])mp[i]=n++; for(int i=0;i<M;i++)if(F[i].active)m++;
    static char ob[1<<20]; setvbuf(stdout,ob,_IOFBF,sizeof(ob)); printf("%d %d\n",n,m); for(int i=0;i<N;i++)if(VAlive[i])printf("v %.17g %.17g %.17g\n",P[i].x,P[i].y,P[i].z); for(int i=0;i<M;i++)if(F[i].active)printf("f %d %d %d\n",mp[F[i].v[0]]+1,mp[F[i].v[1]]+1,mp[F[i].v[2]]+1);
}
static void outputGrid(const vector<Vec3>&V,const vector<array<int,3>>&G){
    static char ob[1<<20]; setvbuf(stdout,ob,_IOFBF,sizeof(ob)); printf("%d %d\n",(int)V.size(),(int)G.size()); for(auto&p:V)printf("v %.17g %.17g %.17g\n",p.x,p.y,p.z); for(auto&f:G)printf("f %d %d %d\n",f[0]+1,f[1]+1,f[2]+1);
}
int main(){
    T0=chrono::steady_clock::now(); readInput();
    vector<Vec3>GV; vector<array<int,3>>GF; if(tryBox(GV,GF)||tryEllipsoid(GV,GF)||tryOrderedGrid(GV,GF)){outputGrid(GV,GF);return 0;}
    simplify(); outputCurrent(); return 0;
}
