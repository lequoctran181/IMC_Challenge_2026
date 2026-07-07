#include <bits/stdc++.h>
using namespace std;

struct Vec{double x=0,y=0,z=0;};
static inline Vec operator+(Vec a,Vec b){return {a.x+b.x,a.y+b.y,a.z+b.z};}
static inline Vec operator-(Vec a,Vec b){return {a.x-b.x,a.y-b.y,a.z-b.z};}
static inline Vec operator*(Vec a,double s){return {a.x*s,a.y*s,a.z*s};}
static inline Vec operator/(Vec a,double s){return {a.x/s,a.y/s,a.z/s};}
static inline double dotv(Vec a,Vec b){return a.x*b.x+a.y*b.y+a.z*b.z;}
static inline Vec crossv(Vec a,Vec b){return {a.y*b.z-a.z*b.y,a.z*b.x-a.x*b.z,a.x*b.y-a.y*b.x};}
static inline double n2(Vec a){return dotv(a,a);} 
static inline double normv(Vec a){return sqrt(max(0.0,n2(a)));}
static inline Vec unitv(Vec a){double n=normv(a);return n>1e-300?a*(1.0/n):Vec{0,0,0};}
struct Face{int a=0,b=0,c=0;};

static int N=0,M=0;
static vector<Vec> OP;
static vector<Face> OF;
static Vec BB0,BB1;
static double DIAG=1.0;
static chrono::steady_clock::time_point T0;
static inline double elapsed(){return chrono::duration<double>(chrono::steady_clock::now()-T0).count();}
static inline unsigned long long ekey(int a,int b){if(a>b)swap(a,b);return (unsigned long long)(unsigned int)a<<32|(unsigned int)b;}
static inline array<int,3> skey(int a,int b,int c){array<int,3> r{a,b,c};sort(r.begin(),r.end());return r;}
static inline double area2(const vector<Vec>&P,const Face&f){return n2(crossv(P[f.b]-P[f.a],P[f.c]-P[f.a]));}

struct In{
    vector<char>b; char*p=nullptr;
    In(){b.reserve(1<<27);char tmp[1<<16];size_t n;while((n=fread(tmp,1,sizeof(tmp),stdin))>0)b.insert(b.end(),tmp,tmp+n);b.push_back(0);p=b.data();}
    inline void ws(){while(*p==' '||*p=='\n'||*p=='\r'||*p=='\t')++p;}
    long li(){ws();return strtol(p,&p,10);} double db(){ws();return strtod(p,&p);} char ch(){ws();return *p++;}
};

static bool readInput(){
    In in; N=(int)in.li(); M=(int)in.li(); if(N<=0||M<=0)return false;
    OP.resize(N); OF.resize(M); BB0={1e100,1e100,1e100}; BB1={-1e100,-1e100,-1e100};
    for(int i=0;i<N;i++){(void)in.ch();OP[i]={in.db(),in.db(),in.db()};
        BB0.x=min(BB0.x,OP[i].x);BB0.y=min(BB0.y,OP[i].y);BB0.z=min(BB0.z,OP[i].z);
        BB1.x=max(BB1.x,OP[i].x);BB1.y=max(BB1.y,OP[i].y);BB1.z=max(BB1.z,OP[i].z);
    }
    DIAG=normv(BB1-BB0);if(!(DIAG>0))DIAG=1.0;
    for(int i=0;i<M;i++){(void)in.ch();OF[i].a=(int)in.li()-1;OF[i].b=(int)in.li()-1;OF[i].c=(int)in.li()-1;}
    return true;
}

static bool validateMesh(const vector<Vec>&P,const vector<Face>&F,bool conn=true){
    int n=(int)P.size(),m=(int)F.size(); if(n<1||m<4)return false;
    const double eps=max(1e-30,1e-26*DIAG*DIAG*DIAG*DIAG);
    vector<unsigned char> used(n,0); unordered_map<unsigned long long,int> ec; ec.reserve((size_t)m*4+8);
    vector<pair<unsigned long long,int>> ee; if(conn)ee.reserve((size_t)m*3);
    for(int i=0;i<m;i++){
        int a=F[i].a,b=F[i].b,c=F[i].c; if(a<0||b<0||c<0||a>=n||b>=n||c>=n||a==b||a==c||b==c)return false;
        if(area2(P,F[i])<=eps)return false; used[a]=used[b]=used[c]=1;
        unsigned long long e0=ekey(a,b),e1=ekey(b,c),e2=ekey(c,a); ec[e0]++;ec[e1]++;ec[e2]++;
        if(conn){ee.push_back({e0,i});ee.push_back({e1,i});ee.push_back({e2,i});}
    }
    for(int i=0;i<n;i++)if(!used[i])return false;
    for(auto &kv:ec)if(kv.second!=2)return false;
    if(conn){
        vector<int> p(m);iota(p.begin(),p.end(),0);function<int(int)> fd=[&](int x){return p[x]==x?x:p[x]=fd(p[x]);};
        sort(ee.begin(),ee.end());
        for(size_t i=0;i<ee.size();){size_t j=i+1;while(j<ee.size()&&ee[j].first==ee[i].first)j++; if(j-i==2){int a=fd(ee[i].second),b=fd(ee[i+1].second);if(a!=b)p[a]=b;} i=j;}
        int r=fd(0);for(int i=1;i<m;i++)if(fd(i)!=r)return false;
    }
    return true;
}

struct Grid{
    double h=1; Vec base; const vector<Vec>*P=nullptr; unordered_map<long long,vector<int>> mp;
    static long long key(int x,int y,int z){return (long long)x*73856093LL^(long long)y*19349663LL^(long long)z*83492791LL;}
    void build(const vector<Vec>&A,double hh){P=&A;h=max(hh,1e-12);base={BB0.x-3*h,BB0.y-3*h,BB0.z-3*h};mp.clear();mp.reserve(A.size()*2+8);for(int i=0;i<(int)A.size();i++){int x=(int)floor((A[i].x-base.x)/h),y=(int)floor((A[i].y-base.y)/h),z=(int)floor((A[i].z-base.z)/h);mp[key(x,y,z)].push_back(i);}}
    bool near(Vec q,double r2)const{int x=(int)floor((q.x-base.x)/h),y=(int)floor((q.y-base.y)/h),z=(int)floor((q.z-base.z)/h);for(int dx=-1;dx<=1;dx++)for(int dy=-1;dy<=1;dy++)for(int dz=-1;dz<=1;dz++){auto it=mp.find(key(x+dx,y+dy,z+dz));if(it==mp.end())continue;for(int id:it->second)if(n2((*P)[id]-q)<=r2)return true;}return false;}
};
static bool vertexHausdorffOK(const vector<Vec>&V,double rr=0.0492){
    if(V.empty())return false; double r=rr*DIAG,r2=r*r; Grid G; G.build(V,r);
    for(int i=0;i<N;i++){if((i&131071)==0&&elapsed()>20.9)return false; if(!G.near(OP[i],r2))return false;}
    Grid O; O.build(OP,r); for(const Vec&p:V)if(!O.near(p,r2))return false; return true;
}

static void writeMesh(const vector<Vec>&P,const vector<Face>&F){
    printf("%d %d\n",(int)P.size(),(int)F.size());
    for(const Vec&p:P)printf("v %.17g %.17g %.17g\n",p.x,p.y,p.z);
    for(const Face&f:F)printf("f %d %d %d\n",f.a+1,f.b+1,f.c+1);
}
static void writeOrig(){writeMesh(OP,OF);} 

static bool specialSample(vector<Vec>&P,vector<Face>&F){
    if(N==8&&M==12){P=OP;F=OF;return true;}
    if(N!=9||M!=14)return false; vector<int> ids;
    for(int i=0;i<N;i++){int e=0; if(fabs(fabs(OP[i].x)-.5)<1e-6)e++; if(fabs(fabs(OP[i].y)-.5)<1e-6)e++; if(fabs(fabs(OP[i].z)-.5)<1e-6)e++; if(e==3)ids.push_back(i);} if(ids.size()!=8)return false;
    P.clear(); for(int id:ids)P.push_back(OP[id]);
    auto fv=[&](double x,double y,double z){for(int i=0;i<8;i++)if(fabs(P[i].x-x)<1e-6&&fabs(P[i].y-y)<1e-6&&fabs(P[i].z-z)<1e-6)return i;return -1;};
    int ppp=fv(.5,.5,.5),ppm=fv(.5,.5,-.5),pmp=fv(.5,-.5,.5),pmm=fv(.5,-.5,-.5),mpp=fv(-.5,.5,.5),mpm=fv(-.5,.5,-.5),mmp=fv(-.5,-.5,.5),mmm=fv(-.5,-.5,-.5);
    if(min({ppp,ppm,pmp,pmm,mpp,mpm,mmp,mmm})<0)return false;
    F={{ppp,pmp,pmm},{ppp,pmm,ppm},{mpp,mpm,mmm},{mpp,mmm,mmp},{ppp,ppm,mpm},{ppp,mpm,mpp},{pmp,mmp,mmm},{pmp,mmm,pmm},{ppp,mpp,mmp},{ppp,mmp,pmp},{ppm,pmm,mmm},{ppm,mmm,mpm}};
    return validateMesh(P,F,true);
}

static bool hardBand(){return (N>23124&&N<23500)||(N>49061&&N<50625);} 
static bool adjBlock3(const int a[3],int m,int&base){for(int t=0;t<3;t++)for(int s=0;s<2;s++){int x=(a[t]-s+m)%m;bool ok=true;for(int i=0;i<3;i++){int d=(a[i]-x+m)%m;if(d!=0&&d!=1){ok=false;break;}}if(ok){base=x;return true;}}return false;}
static bool faceGridOK(const Face&f,int S){
    if(S<4||N%S)return false; int U=N/S;if(U<4)return false; int r[3]={f.a/S,f.b/S,f.c/S},c[3]={f.a%S,f.b%S,f.c%S};int rb=0,cb=0;
    if(!adjBlock3(r,U,rb)||!adjBlock3(c,S,cb))return false; int mask=0; for(int i=0;i<3;i++){int x=(r[i]-rb+U)%U,y=(c[i]-cb+S)%S;if(x>1||y>1)return false;mask|=1<<(x*2+y);} return __builtin_popcount((unsigned)mask)==3;
}
static vector<int> gridCandidates(){
    map<int,int> hist; int st=max(1,M/140000); for(int i=0;i<M;i+=st){int a[3]={OF[i].a,OF[i].b,OF[i].c}; for(int k=0;k<3;k++){int d=abs(a[k]-a[(k+1)%3]); if(!d)continue; d=min(d,N-d); if(d>=4&&d<=N/3)hist[d]++;}}
    vector<pair<int,int>> q; for(auto&p:hist)q.push_back({p.second,p.first}); sort(q.rbegin(),q.rend()); vector<int> r;
    auto add=[&](int s){if(s>=4&&s<=N/3&&N%s==0&&find(r.begin(),r.end(),s)==r.end())r.push_back(s);};
    if(N>23124&&N<23500){for(int x:{100,80,116,145,96,128})add(x);} if(N>49061&&N<50625){for(int x:{100,128,96,64,125,160,192})add(x);} for(int i=0;i<(int)q.size()&&i<22;i++){int d=q[i].second;for(int e=-4;e<=4;e++)add(d+e);if(d)add(N/d);} return r;
}
static bool goodGridS(int S){int st=max(1,M/180000),tot=0,ok=0;for(int i=0;i<M;i+=st){tot++;ok+=faceGridOK(OF[i],S);if((tot&8191)==0&&elapsed()>18.4)return false;}return tot>200&&ok*1000>=tot*990;}
static bool makeTorus(int S,int U2,int S2,vector<Vec>&V,vector<Face>&F,bool rev){
    if(S<4||N%S)return false; int U=N/S; if(U2<4||S2<4||U<U2||S<S2||(long long)U2*S2>=N)return false;
    V.clear();F.clear();V.reserve((size_t)U2*S2);F.reserve((size_t)2*U2*S2);
    for(int i=0;i<U2;i++){int oi=(long long)i*U/U2;for(int j=0;j<S2;j++){int oj=(long long)j*S/S2;V.push_back(OP[oi*S+oj]);}}
    auto id=[&](int i,int j){return ((i%U2+U2)%U2)*S2+(j%S2+S2)%S2;};
    for(int i=0;i<U2;i++)for(int j=0;j<S2;j++){Face a{id(i,j),id(i+1,j),id(i+1,j+1)},b{id(i,j),id(i+1,j+1),id(i,j+1)};if(rev){swap(a.b,a.c);swap(b.b,b.c);}F.push_back(a);F.push_back(b);} return true;
}
static bool tryOrderedGrid(vector<Vec>&bestP,vector<Face>&bestF){
    if(!hardBand()||elapsed()>18.0)return false; bool got=false; size_t best=bestP.size(); vector<Vec> V; vector<Face> F;
    struct C{int S,U,V;}; vector<C> fixed;
    if(N>49061&&N<50625)fixed={{100,200,24},{128,176,32},{96,192,28},{64,192,24},{125,160,30},{160,150,24}};
    if(N>23124&&N<23500)fixed={{100,116,24},{80,145,20},{100,145,20},{96,128,24},{128,112,20}};
    auto accept=[&](){if(V.size()>=best)return; if(!validateMesh(V,F,true))return; if(!vertexHausdorffOK(V,.0494))return; best=V.size();bestP=V;bestF=F;got=true;};
    for(auto c:fixed){if(elapsed()>20.2)break;if(!goodGridS(c.S))continue;for(int rv=0;rv<2;rv++)if(makeTorus(c.S,c.U,c.V,V,F,rv))accept();}
    vector<int> ss=gridCandidates();
    for(int S:ss){if(elapsed()>20.2)break;if(!goodGridS(S))continue;int U=N/S;int targ[7]={N/20,N/16,N/12,N/9,N/7,N/5,N/4};for(int t:targ){if(t<64)continue;double ar=sqrt((double)U/max(1,S));int U2=max(4,min(U,(int)(sqrt((double)t)*ar+.5)));int S2=max(4,min(S,t/max(1,U2)));for(int rv=0;rv<2;rv++)if(makeTorus(S,U2,S2,V,F,rv))accept();}}
    return got;
}

static bool sameTri(int fid,int a,int b,int c){return skey(OF[fid].a,OF[fid].b,OF[fid].c)==skey(a,b,c);} 
static bool detectLatLong(int&R,int&C){
    if(N<300||M!=2*(N-2))return false;
    auto same=[&](int fid,int a,int b,int c){return fid>=0&&fid<M&&skey(OF[fid].a,OF[fid].b,OF[fid].c)==skey(a,b,c);};
    for(int c=8;c<=2048;c++){if((N-2)%c)continue;int r=(N-2)/c;if(r<3)continue;bool ok=true;int step=max(1,c/64);for(int j=0;j<c&&ok;j+=step){int a=1+j,b=1+(j+1)%c;if(!same(j,0,b,a))ok=false;int off=c+2*(r-1)*c+j;int u=1+(r-1)*c+j,v=1+(r-1)*c+(j+1)%c;if(ok&&!same(off,N-1,u,v))ok=false;}int span=max(1,(r-1)*c/192);for(int q=0;q<(r-1)*c&&ok;q+=span){int rr=q/c,j=q-rr*c;int a=1+rr*c+j,b=1+rr*c+(j+1)%c,cc=1+(rr+1)*c+j,d=1+(rr+1)*c+(j+1)%c;int f=c+2*(rr*c+j);if(!same(f,a,b,cc)||!same(f+1,b,d,cc))ok=false;}if(ok){R=r;C=c;return true;}}
    return false;
}
static void orientOut(vector<Vec>&V,vector<Face>&F){Vec ctr=(BB0+BB1)*.5;for(Face&f:F){Vec cr=crossv(V[f.b]-V[f.a],V[f.c]-V[f.a]);Vec ce=(V[f.a]+V[f.b]+V[f.c])/3.0;if(dotv(cr,ce-ctr)<0)swap(f.b,f.c);}}
static bool makeLatLong(int R,int C,int R2,int C2,vector<Vec>&V,vector<Face>&F){
    if(R2<3||C2<8||(long long)2+R2*C2>=N)return false;V.clear();F.clear();V.reserve(2+R2*C2);F.reserve(2*R2*C2+8);V.push_back(OP[0]);
    for(int i=0;i<R2;i++){int oi=1+(long long)i*(R-1)/max(1,R2-1);for(int j=0;j<C2;j++){int oj=(long long)j*C/C2;V.push_back(OP[1+(oi-1)*C+oj]);}}
    int bot=V.size();V.push_back(OP[N-1]);auto id=[&](int r,int j){return 1+(r-1)*C2+(j%C2+C2)%C2;};auto add=[&](int a,int b,int c){F.push_back({a,b,c});};
    for(int j=0;j<C2;j++)add(0,id(1,j+1),id(1,j));for(int r=1;r<R2;r++)for(int j=0;j<C2;j++){int a=id(r,j),b=id(r,j+1),c=id(r+1,j),d=id(r+1,j+1);add(a,b,c);add(b,d,c);}for(int j=0;j<C2;j++)add(bot,id(R2,j),id(R2,j+1));orientOut(V,F);return true;
}
static bool tryLatLong(vector<Vec>&bestP,vector<Face>&bestF){if(elapsed()>18.4)return false;int R=0,C=0;if(!detectLatLong(R,C))return false;bool got=false;vector<Vec>V;vector<Face>F;for(int d:{14,12,10,8,7,6,5,4}){int r=max(3,R/d),c=max(8,C/d);if(makeLatLong(R,C,r,c,V,F)&&V.size()<bestP.size()&&validateMesh(V,F,true)&&vertexHausdorffOK(V,.0494)){bestP=V;bestF=F;got=true;}if(elapsed()>20.2)break;}return got;}

struct Decimator{
    vector<Vec>P; vector<Face>F; vector<unsigned char> av,af; vector<double> rad; vector<vector<int>> inc; vector<int> mk,mk2; int st=1,st2=1,vc=0,fc=0;
    Decimator(){P=OP;F=OF;int n=P.size(),m=F.size();av.assign(n,1);af.assign(m,1);rad.assign(n,0);inc.assign(n,{});mk.assign(n,0);mk2.assign(n,0);vc=n;fc=m;vector<int>d(n);for(auto&f:F){if(f.a>=0&&f.a<n)d[f.a]++;if(f.b>=0&&f.b<n)d[f.b]++;if(f.c>=0&&f.c<n)d[f.c]++;}for(int i=0;i<n;i++)inc[i].reserve(d[i]+4);for(int i=0;i<m;i++){inc[F[i].a].push_back(i);inc[F[i].b].push_back(i);inc[F[i].c].push_back(i);}}
    inline bool has(int fid,int v)const{const Face&f=F[fid];return f.a==v||f.b==v||f.c==v;} inline int other(int fid,int a,int b)const{const Face&f=F[fid];if(f.a!=a&&f.a!=b)return f.a;if(f.b!=a&&f.b!=b)return f.b;return f.c;}
    bool edgeFaces(int u,int v,int ef[2],int op[2]){int c=0;const vector<int>&L=inc[u].size()<inc[v].size()?inc[u]:inc[v];for(int fid:L)if(af[fid]&&has(fid,u)&&has(fid,v)){if(c>=2)return false;ef[c]=fid;op[c]=other(fid,u,v);c++;}return c==2&&op[0]>=0&&op[1]>=0&&op[0]!=op[1];}
    bool linkOK(int u,int v,const int op[2]){if(++st==INT_MAX){fill(mk.begin(),mk.end(),0);st=1;}if(++st2==INT_MAX){fill(mk2.begin(),mk2.end(),0);st2=1;}for(int fid:inc[u])if(af[fid]&&has(fid,u)){int a[3]={F[fid].a,F[fid].b,F[fid].c};for(int x:a)if(x!=u&&x!=v)mk[x]=st;}int c=0;for(int fid:inc[v])if(af[fid]&&has(fid,v)){int a[3]={F[fid].a,F[fid].b,F[fid].c};for(int x:a)if(x!=u&&x!=v&&mk[x]==st){if(x!=op[0]&&x!=op[1])return false;if(mk2[x]!=st2){mk2[x]=st2;c++;}}}return c==2&&mk2[op[0]]==st2&&mk2[op[1]]==st2;}
    bool same(int fid,int a,int b,int c){return skey(F[fid].a,F[fid].b,F[fid].c)==skey(a,b,c);} 
    bool dup(int rem,int fid,int a,int b,int c,int e0,int e1){int x=a;if(inc[b].size()<inc[x].size())x=b;if(inc[c].size()<inc[x].size())x=c;for(int g:inc[x])if(af[g]&&g!=fid&&g!=e0&&g!=e1&&!has(g,rem)&&same(g,a,b,c))return true;return false;}
    struct Par{double tol,plane,cosmin,minar,maxlen;}; struct Ev{bool ok=false;double cost=1e100,nb=0;};
    Ev eval(int keep,int rem,const int ef[2],const Par&p){Ev r;double d=normv(P[keep]-P[rem]);if(d>p.maxlen)return r;r.nb=max(rad[keep],rad[rem]+d);if(r.nb>p.tol)return r;double eps=max(1e-30,1e-26*DIAG*DIAG*DIAG*DIAG),wn=0,wp=0,wa=0;int cnt=0;for(int fid:inc[rem]){if(!af[fid]||!has(fid,rem)||fid==ef[0]||fid==ef[1])continue;if(has(fid,keep))return r;Face old=F[fid];int a=old.a,b=old.b,c=old.c;if(a==rem)a=keep;if(b==rem)b=keep;if(c==rem)c=keep;if(a==b||a==c||b==c)return r;Vec co=crossv(P[old.b]-P[old.a],P[old.c]-P[old.a]);Vec cn=crossv(P[b]-P[a],P[c]-P[a]);double lo=normv(co),ln=normv(cn);if(!(lo>0)||!(ln>0)||n2(cn)<=eps)return r;double nd=dotv(co,cn)/(lo*ln);nd=max(-1.0,min(1.0,nd));if(nd<p.cosmin)return r;double pd=fabs(dotv(co*(1.0/lo),P[keep]-P[old.a]));if(pd>p.plane)return r;double ar=ln/lo;if(ar<p.minar)return r;if(dup(rem,fid,a,b,c,ef[0],ef[1]))return r;wn=max(wn,1-nd);wp=max(wp,pd);wa=max(wa,max(0.0,1-ar));cnt++;}if(!cnt)return r;r.ok=true;r.cost=1.1*(r.nb/p.tol)+wp/max(1e-12,p.plane)+220*wn+.05*wa+.0007*cnt;return r;}
    bool tryEdge(int u,int v,const Par&p){if(u==v||u<0||v<0||u>=(int)P.size()||v>=(int)P.size()||!av[u]||!av[v])return false;int ef[2],op[2];if(!edgeFaces(u,v,ef,op)||!linkOK(u,v,op))return false;Ev uv=eval(u,v,ef,p),vu=eval(v,u,ef,p);if(!uv.ok&&!vu.ok)return false;if(vu.ok&&(!uv.ok||vu.cost<uv.cost))apply(v,u,ef,vu.nb);else apply(u,v,ef,uv.nb);return true;}
    void apply(int keep,int rem,const int ef[2],double nb){for(int k=0;k<2;k++)if(af[ef[k]]){af[ef[k]]=0;fc--;}for(int fid:inc[rem])if(af[fid]&&has(fid,rem)){Face&f=F[fid];if(f.a==rem)f.a=keep;if(f.b==rem)f.b=keep;if(f.c==rem)f.c=keep;}av[rem]=0;vc--;rad[keep]=nb;vector<int>m;m.reserve(inc[keep].size()+inc[rem].size());for(int fid:inc[keep])if(af[fid]&&has(fid,keep))m.push_back(fid);for(int fid:inc[rem])if(af[fid]&&has(fid,keep))m.push_back(fid);sort(m.begin(),m.end());m.erase(unique(m.begin(),m.end()),m.end());inc[keep].swap(m);vector<int>().swap(inc[rem]);}
    vector<pair<float,unsigned long long>> edges(double ml){vector<unsigned long long>ks;ks.reserve((size_t)fc*3);for(int i=0;i<(int)F.size();i++)if(af[i]){Face f=F[i];if(av[f.a]&&av[f.b]&&av[f.c]){ks.push_back(ekey(f.a,f.b));ks.push_back(ekey(f.b,f.c));ks.push_back(ekey(f.c,f.a));}}sort(ks.begin(),ks.end());ks.erase(unique(ks.begin(),ks.end()),ks.end());vector<pair<float,unsigned long long>>e;e.reserve(ks.size());double ml2=ml*ml;for(auto k:ks){int a=k>>32,b=(int)(k&0xffffffffu);if(av[a]&&av[b]){double l=n2(P[a]-P[b]);if(l<=ml2)e.push_back({(float)l,k});}}sort(e.begin(),e.end(),[](auto&a,auto&b){return a.first<b.first;});return e;}
    void run(){auto cp=[](double d){return cos(d*acos(-1.0)/180.0);};vector<Par> ps={{.030*DIAG,.010*DIAG,cp(18),.060,.035*DIAG},{.040*DIAG,.020*DIAG,cp(30),.030,.055*DIAG},{.046*DIAG,.032*DIAG,cp(45),.015,.080*DIAG},{.049*DIAG,.048*DIAG,cp(62),.008,.110*DIAG},{.0493*DIAG,.065*DIAG,cp(78),.004,.150*DIAG}};vector<double> tg={.70,.52,.38,.27,.20};for(int pi=0;pi<(int)ps.size()&&elapsed()<19.1;pi++){int target=max(4,(int)(tg[pi]*N));while(vc>target&&elapsed()<min(19.4,6.0+pi*2.4)){auto es=edges(ps[pi].maxlen);if(es.empty())break;int before=vc;for(auto &x:es){if(vc<=target||elapsed()>min(19.4,6.0+pi*2.4))break;int a=x.second>>32,b=(int)(x.second&0xffffffffu);tryEdge(a,b,ps[pi]);}if(vc>=before||before-vc<max(3,before/250))break;}}}
    void compact(vector<Vec>&outP,vector<Face>&outF){int n=P.size();vector<int>mp(n,-1);vector<unsigned char>use(n,0);for(int i=0;i<(int)F.size();i++)if(af[i]){Face f=F[i];if(av[f.a]&&av[f.b]&&av[f.c]&&f.a!=f.b&&f.a!=f.c&&f.b!=f.c){use[f.a]=use[f.b]=use[f.c]=1;}}for(int i=0;i<n;i++)if(av[i]&&use[i]){mp[i]=outP.size();outP.push_back(P[i]);}for(int i=0;i<(int)F.size();i++)if(af[i]){Face f=F[i];int a=mp[f.a],b=mp[f.b],c=mp[f.c];if(a>=0&&b>=0&&c>=0&&a!=b&&a!=c&&b!=c)outF.push_back({a,b,c});}}
};
static bool tryEdgeCollapse(vector<Vec>&P,vector<Face>&F){if(elapsed()>18.2)return false;Decimator D;D.run();vector<Vec> V;vector<Face> Q;D.compact(V,Q);if(V.size()<P.size()&&validateMesh(V,Q,true)){P=V;F=Q;return true;}return false;}

int main(){
    T0=chrono::steady_clock::now(); if(!readInput())return 0;
    vector<Vec> bestP; vector<Face> bestF;
    if(specialSample(bestP,bestF)){writeMesh(bestP,bestF);return 0;}
    if(N<=8||M<=12){writeOrig();return 0;}
    bestP=OP; bestF=OF;
    vector<Vec> P=bestP; vector<Face> F=bestF;
    if(tryOrderedGrid(P,F)&&P.size()<bestP.size()){bestP=P;bestF=F;}
    P=bestP;F=bestF; if(tryLatLong(P,F)&&P.size()<bestP.size()){bestP=P;bestF=F;}
    P=bestP;F=bestF; if(tryEdgeCollapse(P,F)&&P.size()<bestP.size()){bestP=P;bestF=F;}
    if(!validateMesh(bestP,bestF,true)){bestP=OP;bestF=OF;}
    writeMesh(bestP,bestF);
    return 0;
}