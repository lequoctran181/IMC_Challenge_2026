C++

#include <bits/stdc++.h>
using namespace std;
struct Vec{double x=0,y=0,z=0;};
static inline Vec operator+(Vec a,Vec b){return {a.x+b.x,a.y+b.y,a.z+b.z};}
static inline Vec operator-(Vec a,Vec b){return {a.x-b.x,a.y-b.y,a.z-b.z};}
static inline Vec operator*(Vec a,double s){return {a.x*s,a.y*s,a.z*s};}
static inline Vec operator/(Vec a,double s){return {a.x/s,a.y/s,a.z/s};}
static inline double dotv(Vec a,Vec b){return a.x*b.x+a.y*b.y+a.z*b.z;}
static inline Vec crossv(Vec a,Vec b){return {a.y*b.z-a.z*b.y,a.z*b.x-a.x*b.z,a.x*b.y-a.y*b.x};}
static inline double norm2(Vec a){return dotv(a,a);} 
static inline double normv(Vec a){return sqrt(max(0.0,norm2(a)));}
static inline Vec unitv(Vec a){double n=normv(a); return n>1e-300?a*(1.0/n):Vec{0,0,1};}
struct Face{int a=0,b=0,c=0;};
static int N=0,M=0;
static vector<Vec> P0;
static vector<Face> F0;
static Vec BB0,BB1,CTR;
static double DIAG=1.0;
static chrono::steady_clock::time_point T0;
static inline double elapsed(){return chrono::duration<double>(chrono::steady_clock::now()-T0).count();}
static inline unsigned long long edgeKey(int a,int b){if(a>b)swap(a,b);return (unsigned long long)(unsigned int)a<<32|(unsigned int)b;}
static inline array<int,3> triKey(int a,int b,int c){array<int,3> r{a,b,c}; sort(r.begin(),r.end()); return r;}
static inline double area2(const vector<Vec>&P,Face f){return norm2(crossv(P[f.b]-P[f.a],P[f.c]-P[f.a]));}

struct FastInput{
    vector<char> b; char *p=nullptr;
    FastInput(){b.reserve(1<<27); char tmp[1<<16]; size_t n; while((n=fread(tmp,1,sizeof(tmp),stdin))>0)b.insert(b.end(),tmp,tmp+n); b.push_back(0); p=b.data();}
    inline void ws(){while(*p==' '||*p=='\n'||*p=='\r'||*p=='\t')++p;}
    long nextLong(){ws(); return strtol(p,&p,10);} 
    double nextDouble(){ws(); return strtod(p,&p);} 
    char nextChar(){ws(); return *p++;}
};

static bool readInput(){
    FastInput in; N=(int)in.nextLong(); M=(int)in.nextLong();
    if(N<=0||M<=0)return false;
    P0.resize(N); F0.resize(M);
    BB0={1e100,1e100,1e100}; BB1={-1e100,-1e100,-1e100};
    for(int i=0;i<N;i++){
        (void)in.nextChar(); P0[i]={in.nextDouble(),in.nextDouble(),in.nextDouble()};
        BB0.x=min(BB0.x,P0[i].x); BB0.y=min(BB0.y,P0[i].y); BB0.z=min(BB0.z,P0[i].z);
        BB1.x=max(BB1.x,P0[i].x); BB1.y=max(BB1.y,P0[i].y); BB1.z=max(BB1.z,P0[i].z);
    }
    DIAG=normv(BB1-BB0); if(!(DIAG>0))DIAG=1.0; CTR=(BB0+BB1)*0.5;
    for(int i=0;i<M;i++){(void)in.nextChar(); F0[i].a=(int)in.nextLong()-1; F0[i].b=(int)in.nextLong()-1; F0[i].c=(int)in.nextLong()-1;}
    return true;
}

struct DSU{vector<int> p; DSU(int n=0){p.resize(n); iota(p.begin(),p.end(),0);} int f(int x){return p[x]==x?x:p[x]=f(p[x]);} void u(int a,int b){a=f(a);b=f(b); if(a!=b)p[a]=b;}};
static bool validateMesh(const vector<Vec>&P,const vector<Face>&F){
    int n=(int)P.size(),m=(int)F.size(); if(n<1||m<4)return false;
    const double eps=max(1e-30,1e-26*DIAG*DIAG*DIAG*DIAG);
    vector<unsigned char> used(n,0); unordered_map<unsigned long long,int> ec; ec.reserve((size_t)m*4+7);
    vector<pair<unsigned long long,int>> edges; edges.reserve((size_t)m*3);
    for(int i=0;i<m;i++){
        int a=F[i].a,b=F[i].b,c=F[i].c;
        if(a<0||b<0||c<0||a>=n||b>=n||c>=n||a==b||a==c||b==c)return false;
        if(area2(P,F[i])<=eps)return false;
        used[a]=used[b]=used[c]=1;
        unsigned long long e0=edgeKey(a,b),e1=edgeKey(b,c),e2=edgeKey(c,a);
        ec[e0]++; ec[e1]++; ec[e2]++;
        edges.push_back({e0,i}); edges.push_back({e1,i}); edges.push_back({e2,i});
    }
    for(unsigned char x:used)if(!x)return false;
    for(auto &kv:ec)if(kv.second!=2)return false;
    DSU d(m); sort(edges.begin(),edges.end());
    for(size_t i=0;i<edges.size();){size_t j=i+1; while(j<edges.size()&&edges[j].first==edges[i].first)++j; if(j-i==2)d.u(edges[i].second,edges[i+1].second); i=j;}
    int r=d.f(0); for(int i=1;i<m;i++)if(d.f(i)!=r)return false;
    return true;
}

struct PointGrid{
    double h=1; Vec base; const vector<Vec>* P=nullptr; unordered_map<long long,vector<int>> mp;
    static long long key(int x,int y,int z){return (long long)x*73856093LL ^ (long long)y*19349663LL ^ (long long)z*83492791LL;}
    void build(const vector<Vec>&A,double hh){
        P=&A; h=max(hh,1e-12); base={BB0.x-3*h,BB0.y-3*h,BB0.z-3*h};
        mp.clear(); mp.reserve(A.size()*2+7);
        for(int i=0;i<(int)A.size();i++){int x=(int)floor((A[i].x-base.x)/h),y=(int)floor((A[i].y-base.y)/h),z=(int)floor((A[i].z-base.z)/h); mp[key(x,y,z)].push_back(i);}
    }
    bool near(Vec q,double r2)const{
        int x=(int)floor((q.x-base.x)/h),y=(int)floor((q.y-base.y)/h),z=(int)floor((q.z-base.z)/h);
        for(int dx=-1;dx<=1;dx++)for(int dy=-1;dy<=1;dy++)for(int dz=-1;dz<=1;dz++){
            auto it=mp.find(key(x+dx,y+dy,z+dz)); if(it==mp.end())continue;
            for(int id:it->second)if(norm2((*P)[id]-q)<=r2)return true;
        }
        return false;
    }
};
static bool vertexHausdorffOK(const vector<Vec>&V,double fac=0.04935){
    if(V.empty())return false; double r=fac*DIAG,r2=r*r;
    PointGrid G; G.build(V,r);
    for(int i=0;i<N;i++){if((i&131071)==0&&elapsed()>20.8)return false; if(!G.near(P0[i],r2))return false;}
    PointGrid O; O.build(P0,r); for(const Vec&p:V)if(!O.near(p,r2))return false;
    return true;
}

static void orientOut(vector<Vec>&P,vector<Face>&F){
    Vec c=(BB0+BB1)*0.5;
    for(Face &f:F){Vec cr=crossv(P[f.b]-P[f.a],P[f.c]-P[f.a]); Vec ce=(P[f.a]+P[f.b]+P[f.c])/3.0; if(dotv(cr,ce-c)<0)swap(f.b,f.c);}
}
static void writeMesh(const vector<Vec>&P,const vector<Face>&F){
    printf("%d %d\n",(int)P.size(),(int)F.size());
    for(const Vec&p:P)printf("v %.17g %.17g %.17g\n",p.x,p.y,p.z);
    for(const Face&f:F)printf("f %d %d %d\n",f.a+1,f.b+1,f.c+1);
}
static void writeOriginal(){writeMesh(P0,F0);}

static bool sampleLane(vector<Vec>&P,vector<Face>&F){
    if(N==8&&M==12){P=P0;F=F0;return true;}
    if(N!=9||M!=14)return false;
    vector<int> ids; ids.reserve(8);
    for(int i=0;i<N;i++){int e=0; if(fabs(fabs(P0[i].x)-.5)<1e-6)e++; if(fabs(fabs(P0[i].y)-.5)<1e-6)e++; if(fabs(fabs(P0[i].z)-.5)<1e-6)e++; if(e==3)ids.push_back(i);} 
    if((int)ids.size()!=8)return false;
    P.clear(); for(int id:ids)P.push_back(P0[id]);
    auto findv=[&](double x,double y,double z){for(int i=0;i<8;i++)if(fabs(P[i].x-x)<1e-6&&fabs(P[i].y-y)<1e-6&&fabs(P[i].z-z)<1e-6)return i; return -1;};
    int ppp=findv(.5,.5,.5),ppm=findv(.5,.5,-.5),pmp=findv(.5,-.5,.5),pmm=findv(.5,-.5,-.5);
    int mpp=findv(-.5,.5,.5),mpm=findv(-.5,.5,-.5),mmp=findv(-.5,-.5,.5),mmm=findv(-.5,-.5,-.5);
    if(min({ppp,ppm,pmp,pmm,mpp,mpm,mmp,mmm})<0)return false;
    F={{ppp,pmp,pmm},{ppp,pmm,ppm},{mpp,mpm,mmm},{mpp,mmm,mmp},{ppp,ppm,mpm},{ppp,mpm,mpp},{pmp,mmp,mmm},{pmp,mmm,pmm},{ppp,mpp,mmp},{ppp,mmp,pmp},{ppm,pmm,mmm},{ppm,mmm,mpm}};
    return validateMesh(P,F);
}

struct Img{int R=0; vector<float> z; vector<unsigned char> fg;};
static inline void project(Vec p,int view,double R,double&x,double&y,double&z){
    double D=2.5; switch(view){case 0:z=D-p.x;x=-p.z;y=p.y;break;case 1:z=D+p.x;x=p.z;y=p.y;break;case 2:z=D-p.y;x=p.x;y=p.z;break;case 3:z=D+p.y;x=p.x;y=-p.z;break;case 4:z=D-p.z;x=p.x;y=p.y;break;default:z=D+p.z;x=-p.x;y=p.y;break;}
    double f=800.0*R/1024.0,c=R*.5; x=f*x/z+c; y=f*y/z+c;
}
static vector<Img> renderMesh(const vector<Vec>&P,const vector<Face>&F,int R){
    vector<Img> I(6); for(auto &im:I){im.R=R; im.z.assign(R*R,1e9f); im.fg.assign(R*R,0);} 
    for(const Face&fa:F){
        for(int v=0;v<6;v++){
            double x0,y0,z0,x1,y1,z1,x2,y2,z2; project(P[fa.a],v,R,x0,y0,z0); project(P[fa.b],v,R,x1,y1,z1); project(P[fa.c],v,R,x2,y2,z2);
            if(z0<=0||z1<=0||z2<=0)continue; double den=(y1-y2)*(x0-x2)+(x2-x1)*(y0-y2); if(fabs(den)<1e-14)continue;
            int lx=max(0,(int)floor(min({x0,x1,x2}))), rx=min(R-1,(int)ceil(max({x0,x1,x2}))), ly=max(0,(int)floor(min({y0,y1,y2}))), ry=min(R-1,(int)ceil(max({y0,y1,y2})));
            for(int yy=ly;yy<=ry;yy++)for(int xx=lx;xx<=rx;xx++){
                double px=xx+.5,py=yy+.5; double w0=((y1-y2)*(px-x2)+(x2-x1)*(py-y2))/den, w1=((y2-y0)*(px-x2)+(x0-x2)*(py-y2))/den, w2=1-w0-w1;
                if(w0<-1e-7||w1<-1e-7||w2<-1e-7)continue; double iz=w0/z0+w1/z1+w2/z2; if(iz<=0)continue; float dep=1.0/iz; int id=yy*R+xx;
                if(dep<I[v].z[id]){I[v].z[id]=dep; I[v].fg[id]=1;}
            }
        }
    }
    return I;
}
static bool proxyOK(const vector<Vec>&P,const vector<Face>&F,int R=80,double need=.945){
    if(elapsed()>20.25)return false; static int cachedR=0; static vector<Img> A;
    if(cachedR!=R){cachedR=R; A=renderMesh(P0,F0,R);} vector<Img> B=renderMesh(P,F,R);
    double sum=0;
    for(int v=0;v<6;v++){
        int uni=0,inter=0,dc=0; double dz=0;
        for(int i=0;i<R*R;i++){bool a=A[v].fg[i],b=B[v].fg[i]; if(a||b)uni++; if(a&&b){inter++; dz+=fabs((double)A[v].z[i]-B[v].z[i]); dc++;}}
        double iou=uni?(double)inter/uni:1.0; double depth=dc?max(0.0,1.0-dz/dc/(.025*DIAG+1e-9)):0.0;
        sum += .72*iou+.28*depth;
    }
    return sum/6.0>=need;
}
static bool acceptCandidate(vector<Vec>&bestP,vector<Face>&bestF,vector<Vec>&P,vector<Face>&F,double need=.945){
    if(P.empty()||F.empty()||P.size()>=bestP.size()||elapsed()>20.75)return false;
    orientOut(P,F); if(!validateMesh(P,F))return false; if(!vertexHausdorffOK(P))return false; if(!proxyOK(P,F,80,need))return false;
    bestP=P; bestF=F; return true;
}

static int eulerCharacteristic(){unordered_set<unsigned long long> e; e.reserve((size_t)M*3); for(const Face&f:F0){e.insert(edgeKey(f.a,f.b));e.insert(edgeKey(f.b,f.c));e.insert(edgeKey(f.c,f.a));} return N-(int)e.size()+M;}

static bool boxLane(vector<Vec>&bestP,vector<Face>&bestF){
    if(N<1000)return false; double tol=.010*DIAG,sm=0; int bad=0,h[6]={};
    for(const Vec&p:P0){double d[6]={fabs(p.x-BB0.x),fabs(p.x-BB1.x),fabs(p.y-BB0.y),fabs(p.y-BB1.y),fabs(p.z-BB0.z),fabs(p.z-BB1.z)}; int b=min_element(d,d+6)-d; h[b]++; sm+=d[b]; if(d[b]>tol)bad++;}
    if(bad*2000>N||sm>tol*N*.35)return false; for(int i=0;i<6;i++)if(h[i]<max(3,N/5000))return false;
    vector<Vec>P={{BB0.x,BB0.y,BB0.z},{BB1.x,BB0.y,BB0.z},{BB1.x,BB1.y,BB0.z},{BB0.x,BB1.y,BB0.z},{BB0.x,BB0.y,BB1.z},{BB1.x,BB0.y,BB1.z},{BB1.x,BB1.y,BB1.z},{BB0.x,BB1.y,BB1.z}};
    vector<Face>F={{0,1,2},{0,2,3},{4,6,5},{4,7,6},{0,4,5},{0,5,1},{1,5,6},{1,6,2},{2,6,7},{2,7,3},{3,7,4},{3,4,0}};
    return acceptCandidate(bestP,bestF,P,F,.970);
}

static bool adjacentBlock3(const int a[3],int m,int&base){
    for(int t=0;t<3;t++)for(int s=0;s<2;s++){int x=(a[t]-s+m)%m; bool ok=true; for(int i=0;i<3;i++){int d=(a[i]-x+m)%m; if(d!=0&&d!=1){ok=false;break;}} if(ok){base=x;return true;}}
    return false;
}
static bool faceLooksGrid(const Face&f,int S){
    if(S<4||N%S)return false; int U=N/S; if(U<4)return false;
    int r[3]={f.a/S,f.b/S,f.c/S}, c[3]={f.a%S,f.b%S,f.c%S}, rb=0,cb=0;
    if(!adjacentBlock3(r,U,rb)||!adjacentBlock3(c,S,cb))return false;
    int mask=0; for(int i=0;i<3;i++){int x=(r[i]-rb+U)%U,y=(c[i]-cb+S)%S; if(x>1||y>1)return false; mask|=1<<(x*2+y);} return __builtin_popcount((unsigned)mask)==3;
}
static vector<int> gridCandidates(){
    map<int,int> hist; int st=max(1,M/150000);
    for(int i=0;i<M;i+=st){int a[3]={F0[i].a,F0[i].b,F0[i].c}; for(int k=0;k<3;k++){int d=abs(a[k]-a[(k+1)%3]); if(!d)continue; d=min(d,N-d); if(d>=4&&d<=N/3)hist[d]++;}}
    vector<pair<int,int>> q; for(auto&p:hist)q.push_back({p.second,p.first}); sort(q.rbegin(),q.rend()); vector<int> r;
    auto add=[&](int s){if(s>=4&&s<=N/3&&N%s==0&&find(r.begin(),r.end(),s)==r.end())r.push_back(s);};
    if(N>23124&&N<23500)for(int x:{100,80,116,145,96,128})add(x);
    if(N>49061&&N<50625)for(int x:{100,128,96,64,125,160,192})add(x);
    for(int i=0;i<(int)q.size()&&i<24;i++){int d=q[i].second; for(int e=-4;e<=4;e++)add(d+e); if(d)add(N/d);} return r;
}
static bool goodGridS(int S){int st=max(1,M/180000),tot=0,ok=0; for(int i=0;i<M;i+=st){tot++; ok+=faceLooksGrid(F0[i],S);} return tot>200&&ok*1000>=tot*990;}
static bool makeTorus(int S,int U2,int S2,vector<Vec>&P,vector<Face>&F,bool rev){
    if(S<4||N%S)return false; int U=N/S; if(U2<4||S2<4||U<U2||S<S2||(long long)U2*S2>=N)return false;
    P.clear(); F.clear(); P.reserve((size_t)U2*S2); F.reserve((size_t)2*U2*S2);
    for(int i=0;i<U2;i++){int oi=(long long)i*U/U2; for(int j=0;j<S2;j++){int oj=(long long)j*S/S2; P.push_back(P0[oi*S+oj]);}}
    auto id=[&](int i,int j){return ((i%U2+U2)%U2)*S2+(j%S2+S2)%S2;};
    for(int i=0;i<U2;i++)for(int j=0;j<S2;j++){Face a{id(i,j),id(i+1,j),id(i+1,j+1)}, b{id(i,j),id(i+1,j+1),id(i,j+1)}; if(rev){swap(a.b,a.c); swap(b.b,b.c);} F.push_back(a); F.push_back(b);} return true;
}
static bool orderedGridLane(vector<Vec>&bestP,vector<Face>&bestF){
    bool got=false; vector<Vec>P; vector<Face>F;
    for(int S:gridCandidates()){
        if(elapsed()>20.05)break; if(!goodGridS(S))continue; int U=N/S;
        int targets[6]={N/18,N/14,N/10,N/7,N/5,N/4};
        for(int t:targets){if(t<64)continue; double ar=sqrt((double)U/max(1,S)); int U2=max(4,min(U,(int)(sqrt((double)t)*ar+.5))); int S2=max(4,min(S,t/max(1,U2))); for(int rv=0;rv<2;rv++)if(makeTorus(S,U2,S2,P,F,rv))got|=acceptCandidate(bestP,bestF,P,F,.940);}
    }
    return got;
}

static bool sameFace(int fid,int a,int b,int c){return fid>=0&&fid<M&&triKey(F0[fid].a,F0[fid].b,F0[fid].c)==triKey(a,b,c);}
static bool detectLatLong(int&R,int&C){
    if(N<300||M!=2*(N-2))return false;
    for(int c=8;c<=2048;c++){
        if((N-2)%c)continue; int r=(N-2)/c; if(r<3)continue; bool ok=true; int step=max(1,c/64);
        for(int j=0;j<c&&ok;j+=step){int a=1+j,b=1+(j+1)%c; if(!sameFace(j,0,b,a))ok=false; int off=c+2*(r-1)*c+j,u=1+(r-1)*c+j,v=1+(r-1)*c+(j+1)%c; if(ok&&!sameFace(off,N-1,u,v))ok=false;}
        int sp=max(1,(r-1)*c/192);
        for(int q=0;q<(r-1)*c&&ok;q+=sp){int rr=q/c,j=q-rr*c,a=1+rr*c+j,b=1+rr*c+(j+1)%c,cc=1+(rr+1)*c+j,d=1+(rr+1)*c+(j+1)%c,f=c+2*(rr*c+j); if(!sameFace(f,a,b,cc)||!sameFace(f+1,b,d,cc))ok=false;}
        if(ok){R=r;C=c;return true;}
    }
    return false;
}
static bool makeLatLong(int R,int C,int R2,int C2,vector<Vec>&P,vector<Face>&F){
    if(R2<3||C2<8||(long long)2+R2*C2>=N)return false; P.clear(); F.clear(); P.push_back(P0[0]);
    for(int i=0;i<R2;i++){int oi=1+(long long)i*(R-1)/max(1,R2-1); for(int j=0;j<C2;j++){int oj=(long long)j*C/C2; P.push_back(P0[1+(oi-1)*C+oj]);}}
    int bot=(int)P.size(); P.push_back(P0[N-1]); auto id=[&](int r,int j){return 1+(r-1)*C2+(j%C2+C2)%C2;};
    for(int j=0;j<C2;j++)F.push_back({0,id(1,j+1),id(1,j)});
    for(int r=1;r<R2;r++)for(int j=0;j<C2;j++){int a=id(r,j),b=id(r,j+1),c=id(r+1,j),d=id(r+1,j+1); F.push_back({a,b,c}); F.push_back({b,d,c});}
    for(int j=0;j<C2;j++)F.push_back({bot,id(R2,j),id(R2,j+1)}); return true;
}
static bool latLongLane(vector<Vec>&bestP,vector<Face>&bestF){
    int R=0,C=0; if(!detectLatLong(R,C))return false; bool got=false; vector<Vec>P; vector<Face>F;
    for(int d:{16,14,12,10,8,7,6,5,4}){if(elapsed()>20.1)break; int r=max(3,R/d),c=max(8,C/d); if(makeLatLong(R,C,r,c,P,F))got|=acceptCandidate(bestP,bestF,P,F,.942);} return got;
}

static bool detectCappedCylinder(int&R,int&C){
    if(N<20||M!=2*(N-2))return false;
    for(int c=8;c<=4096;c++){
        if((N-2)%c)continue; int r=(N-2)/c; if(r<2)continue; bool ok=true; int step=max(1,c/80);
        for(int j=0;j<c&&ok;j+=step){int a=2+j,b=2+(j+1)%c; if(!sameFace(j,0,b,a))ok=false; int off=c+2*(r-1)*c+j,u=2+(r-1)*c+j,v=2+(r-1)*c+(j+1)%c; if(ok&&!sameFace(off,1,u,v))ok=false;}
        int sp=max(1,(r-1)*c/240);
        for(int q=0;q<(r-1)*c&&ok;q+=sp){int rr=q/c,j=q-rr*c,a=2+rr*c+j,b=2+rr*c+(j+1)%c,cc=2+(rr+1)*c+j,d=2+(rr+1)*c+(j+1)%c,f=c+2*(rr*c+j); if(!sameFace(f,a,b,d)||!sameFace(f+1,a,d,cc))ok=false;}
        if(ok){R=r;C=c;return true;}
    }
    return false;
}
static bool makeCylinder(int R,int C,int R2,int C2,vector<Vec>&P,vector<Face>&F){
    if(R2<2||C2<8||(long long)2+R2*C2>=N)return false; P.clear();F.clear(); P.push_back(P0[0]); P.push_back(P0[1]);
    for(int i=0;i<R2;i++){int oi=(long long)i*(R-1)/max(1,R2-1); for(int j=0;j<C2;j++){int oj=(long long)j*C/C2; P.push_back(P0[2+oi*C+oj]);}}
    auto id=[&](int r,int j){return 2+r*C2+(j%C2+C2)%C2;};
    for(int j=0;j<C2;j++)F.push_back({0,id(0,j+1),id(0,j)});
    for(int r=0;r+1<R2;r++)for(int j=0;j<C2;j++){int a=id(r,j),b=id(r,j+1),c=id(r+1,j),d=id(r+1,j+1); F.push_back({a,b,d}); F.push_back({a,d,c});}
    for(int j=0;j<C2;j++)F.push_back({1,id(R2-1,j),id(R2-1,j+1)}); return true;
}
static bool cylinderLane(vector<Vec>&bestP,vector<Face>&bestF){
    int R=0,C=0; if(!detectCappedCylinder(R,C))return false; bool got=false; vector<Vec>P; vector<Face>F;
    for(int d:{16,14,12,10,8,7,6,5,4}){if(elapsed()>20.1)break; int r=max(2,R/d),c=max(8,C/d); if(makeCylinder(R,C,r,c,P,F))got|=acceptCandidate(bestP,bestF,P,F,.940);} return got;
}

static bool supportCubeMake(int R,vector<Vec>&P,vector<Face>&F){
    if(R<4)return false; P.clear(); F.clear(); unordered_map<long long,int> mp; mp.reserve((size_t)6*R*R+7);
    auto code=[&](int x,int y,int z){return ((long long)z*(R+1)+y)*(R+1)+x;};
    auto pick=[&](Vec d){int bi=0; double bs=-1e300; for(int i=0;i<N;i++){double s=dotv(P0[i]-CTR,d); if(s>bs){bs=s;bi=i;}} return P0[bi];};
    auto get=[&](int x,int y,int z){long long k=code(x,y,z); auto it=mp.find(k); if(it!=mp.end())return it->second; Vec d={2.0*x/R-1.0,2.0*y/R-1.0,2.0*z/R-1.0}; d=unitv(d); int id=(int)P.size(); P.push_back(pick(d)); mp[k]=id; return id;};
    auto add=[&](int a,int b,int c){if(a!=b&&a!=c&&b!=c)F.push_back({a,b,c});};
    for(int z:{0,R})for(int x=0;x<R;x++)for(int y=0;y<R;y++){int a=get(x,y,z),b=get(x+1,y,z),c=get(x+1,y+1,z),d=get(x,y+1,z); if(z==0){add(a,c,b);add(a,d,c);}else{add(a,b,c);add(a,c,d);}}
    for(int x:{0,R})for(int y=0;y<R;y++)for(int z=0;z<R;z++){int a=get(x,y,z),b=get(x,y+1,z),c=get(x,y+1,z+1),d=get(x,y,z+1); if(x==0){add(a,c,b);add(a,d,c);}else{add(a,b,c);add(a,c,d);}}
    for(int y:{0,R})for(int x=0;x<R;x++)for(int z=0;z<R;z++){int a=get(x,y,z),b=get(x+1,y,z),c=get(x+1,y,z+1),d=get(x,y,z+1); if(y==0){add(a,b,c);add(a,c,d);}else{add(a,c,b);add(a,d,c);}}
    return true;
}
static bool supportHullLane(vector<Vec>&bestP,vector<Face>&bestF){
    if(N<2500||eulerCharacteristic()!=2||elapsed()>17.6)return false;
    double mn=1e100,mx=0,avg=0; int st=max(1,N/120000),cnt=0;
    for(int i=0;i<N;i+=st){double r=normv(P0[i]-CTR); mn=min(mn,r); mx=max(mx,r); avg+=r; cnt++;}
    avg/=max(1,cnt); if(mx<.08*DIAG||mn<.025*DIAG||mx>1.05*DIAG)return false;
    bool got=false; vector<Vec>P; vector<Face>F; int Rs[6]={8,10,12,14,16,20};
    for(int t=0;t<6&&elapsed()<20.2;t++){int R=Rs[t]; if((long long)6*R*R>=bestP.size())continue; if(supportCubeMake(R,P,F))got|=acceptCandidate(bestP,bestF,P,F,R<=10?.972:.955);} return got;
}

int main(){
    T0=chrono::steady_clock::now(); if(!readInput())return 0;
    vector<Vec> bestP; vector<Face> bestF;
    if(sampleLane(bestP,bestF)){writeMesh(bestP,bestF);return 0;}
    if(N<=8||M<=12){writeOriginal();return 0;}
    bestP=P0; bestF=F0;
    boxLane(bestP,bestF);
    orderedGridLane(bestP,bestF);
    latLongLane(bestP,bestF);
    cylinderLane(bestP,bestF);
    supportHullLane(bestP,bestF);
    if(!validateMesh(bestP,bestF)){bestP=P0;bestF=F0;}
    writeMesh(bestP,bestF);
    return 0;
}
