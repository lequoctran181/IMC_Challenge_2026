#include <bits/stdc++.h>
using namespace std;

struct Vec3{double x=0,y=0,z=0;};
static inline Vec3 operator+(const Vec3&a,const Vec3&b){return {a.x+b.x,a.y+b.y,a.z+b.z};}
static inline Vec3 operator-(const Vec3&a,const Vec3&b){return {a.x-b.x,a.y-b.y,a.z-b.z};}
static inline Vec3 operator*(const Vec3&a,double s){return {a.x*s,a.y*s,a.z*s};}
static inline Vec3 operator/(const Vec3&a,double s){return {a.x/s,a.y/s,a.z/s};}
static inline double dot3(const Vec3&a,const Vec3&b){return a.x*b.x+a.y*b.y+a.z*b.z;}
static inline Vec3 cross3(const Vec3&a,const Vec3&b){return {a.y*b.z-a.z*b.y,a.z*b.x-a.x*b.z,a.x*b.y-a.y*b.x};}
static inline double norm2(const Vec3&a){return dot3(a,a);} 
static inline double norm3(const Vec3&a){return sqrt(max(0.0,norm2(a)));}
static inline double dist2(const Vec3&a,const Vec3&b){return norm2(a-b);} 
static inline double dist3(const Vec3&a,const Vec3&b){return sqrt(dist2(a,b));}
static inline Vec3 normalize3(const Vec3&a){double n=norm3(a); return n>1e-300?a/n:Vec3{0,0,0};}

struct Face{int v[3]; unsigned char alive=1;};
struct FastInput{
    vector<char> b; char *p=nullptr;
    FastInput(){ b.reserve(1<<26); char tmp[1<<16]; size_t n; while((n=fread(tmp,1,sizeof(tmp),stdin))>0)b.insert(b.end(),tmp,tmp+n); b.push_back(0); p=b.data(); }
    inline void ws(){while(*p==' '||*p=='\n'||*p=='\r'||*p=='\t')++p;}
    int ni(){ws(); int s=1;if(*p=='-'){s=-1;++p;} int x=0; while(*p>='0'&&*p<='9'){x=x*10+*p-'0';++p;} return x*s;}
    double nd(){ws(); char*e; double x=strtod(p,&e); p=e; return x;}
    char nc(){ws(); return *p++;}
};

static int N=0,M=0,activeV=0,activeF=0;
static vector<Vec3>P,Orig;
static vector<array<int,3>> F0;
static vector<Face> F;
static vector<vector<int>> adj;
static vector<unsigned char> aliveV, emask;
static vector<int> ver, markv;
static vector<double> radErr, nearErr;
static int stampv=1, extCnt[6];
static Vec3 Bmin,Bmax,Center;
static double diagLen=1, H=1, areaEps=1e-30, featRatio=0;
static chrono::steady_clock::time_point T0;
static inline double elapsed(){return chrono::duration<double>(chrono::steady_clock::now()-T0).count();}

struct OutMesh{vector<Vec3> X; vector<array<int,3>> T; string tag;};
static vector<OutMesh> candidates;

static inline bool hasV(const Face&f,int v){return f.v[0]==v||f.v[1]==v||f.v[2]==v;}
static inline bool hasV0(const array<int,3>&f,int v){return f[0]==v||f[1]==v||f[2]==v;}
static inline uint64_t ekey(int a,int b){if(a>b)swap(a,b);return (uint64_t)(uint32_t)a<<32|(uint32_t)b;}
static inline array<int,3> tkey(int a,int b,int c){array<int,3> t{a,b,c};sort(t.begin(),t.end());return t;}
static inline Vec3 faceNormalCur(const Face&f){return cross3(P[f.v[1]]-P[f.v[0]],P[f.v[2]]-P[f.v[0]]);} 
static inline Vec3 faceNormalOrig(const array<int,3>&f){return cross3(Orig[f[1]]-Orig[f[0]],Orig[f[2]]-Orig[f[0]]);} 
static inline int third(const Face&f,int a,int b){for(int k=0;k<3;k++){int x=f.v[k]; if(x!=a&&x!=b)return x;} return -1;}
static inline bool faceHasEdge(const Face&f,int a,int b){return hasV(f,a)&&hasV(f,b);} 

static void readInput(){
    FastInput in; N=in.ni(); M=in.ni();
    P.resize(N); Orig.resize(N); F0.resize(M); F.resize(M);
    Bmin={1e100,1e100,1e100}; Bmax={-1e100,-1e100,-1e100};
    for(int i=0;i<N;i++){
        (void)in.nc(); P[i].x=in.nd(); P[i].y=in.nd(); P[i].z=in.nd(); Orig[i]=P[i];
        Bmin.x=min(Bmin.x,P[i].x); Bmin.y=min(Bmin.y,P[i].y); Bmin.z=min(Bmin.z,P[i].z);
        Bmax.x=max(Bmax.x,P[i].x); Bmax.y=max(Bmax.y,P[i].y); Bmax.z=max(Bmax.z,P[i].z);
    }
    Center=(Bmin+Bmax)/2.0; diagLen=max(1e-12,norm3(Bmax-Bmin)); H=0.0490*diagLen; areaEps=max(1e-30,diagLen*diagLen*diagLen*diagLen*1e-26);
    vector<int> deg(N,0);
    for(int i=0;i<M;i++){
        (void)in.nc(); int a=in.ni()-1,b=in.ni()-1,c=in.ni()-1;
        F0[i]={a,b,c}; F[i].v[0]=a; F[i].v[1]=b; F[i].v[2]=c; F[i].alive=1;
        if(0<=a&&a<N)deg[a]++; if(0<=b&&b<N)deg[b]++; if(0<=c&&c<N)deg[c]++;
    }
    adj.assign(N,{}); for(int i=0;i<N;i++)adj[i].reserve(deg[i]+4);
    for(int i=0;i<M;i++){adj[F[i].v[0]].push_back(i);adj[F[i].v[1]].push_back(i);adj[F[i].v[2]].push_back(i);} 
    aliveV.assign(N,1); emask.assign(N,0); memset(extCnt,0,sizeof(extCnt));
    double eps=diagLen*1e-11+1e-14;
    for(int i=0;i<N;i++){
        unsigned char m=0;
        if(fabs(P[i].x-Bmin.x)<=eps)m|=1; if(fabs(P[i].x-Bmax.x)<=eps)m|=2;
        if(fabs(P[i].y-Bmin.y)<=eps)m|=4; if(fabs(P[i].y-Bmax.y)<=eps)m|=8;
        if(fabs(P[i].z-Bmin.z)<=eps)m|=16; if(fabs(P[i].z-Bmax.z)<=eps)m|=32;
        emask[i]=m; for(int b=0;b<6;b++)if(m&(1u<<b))extCnt[b]++;
    }
    ver.assign(N,0); markv.assign(N,0); radErr.assign(N,0); nearErr.assign(N,0); activeV=N; activeF=M;
}

static void appendLine(string&out,const char*s,int n){ if(out.size()+n>(1u<<20)){ fwrite(out.data(),1,out.size(),stdout); out.clear(); } out.append(s,s+n); }
static void outputRawOriginal(){
    string out; out.reserve(1<<20); char line[160];
    int n=snprintf(line,sizeof(line),"%d %d\n",N,M); appendLine(out,line,n);
    for(auto&p:Orig){n=snprintf(line,sizeof(line),"v %.10g %.10g %.10g\n",p.x,p.y,p.z); appendLine(out,line,n);} 
    for(auto&f:F0){n=snprintf(line,sizeof(line),"f %d %d %d\n",f[0]+1,f[1]+1,f[2]+1); appendLine(out,line,n);} 
    if(!out.empty())fwrite(out.data(),1,out.size(),stdout);
}
static void outputMesh(const OutMesh&mo){
    string out; out.reserve(1<<20); char line[180];
    int n=snprintf(line,sizeof(line),"%d %d\n",(int)mo.X.size(),(int)mo.T.size()); appendLine(out,line,n);
    int prec=(int)mo.X.size()*2<=N?15:10;
    for(auto&p:mo.X){ if(prec==15)n=snprintf(line,sizeof(line),"v %.15g %.15g %.15g\n",p.x,p.y,p.z); else n=snprintf(line,sizeof(line),"v %.10g %.10g %.10g\n",p.x,p.y,p.z); appendLine(out,line,n);} 
    for(auto&t:mo.T){n=snprintf(line,sizeof(line),"f %d %d %d\n",t[0]+1,t[1]+1,t[2]+1); appendLine(out,line,n);} 
    if(!out.empty())fwrite(out.data(),1,out.size(),stdout);
}

struct GridNear{
    double h=1; Vec3 mn; unordered_map<long long, vector<int>> mp; const vector<Vec3>*pts=nullptr;
    static long long key(int x,int y,int z){return ((long long)x*73856093LL)^((long long)y*19349663LL)^((long long)z*83492791LL);} 
    void build(const vector<Vec3>&p,double cell){pts=&p; h=max(cell,1e-12); mn={1e100,1e100,1e100}; for(auto&q:p){mn.x=min(mn.x,q.x);mn.y=min(mn.y,q.y);mn.z=min(mn.z,q.z);} mp.clear(); mp.reserve(p.size()*2+10); for(int i=0;i<(int)p.size();i++){int x=(int)floor((p[i].x-mn.x)/h),y=(int)floor((p[i].y-mn.y)/h),z=(int)floor((p[i].z-mn.z)/h); mp[key(x,y,z)].push_back(i);} }
    bool near(const Vec3&q,double r2)const{int X=(int)floor((q.x-mn.x)/h),Y=(int)floor((q.y-mn.y)/h),Z=(int)floor((q.z-mn.z)/h); for(int x=X-1;x<=X+1;x++)for(int y=Y-1;y<=Y+1;y++)for(int z=Z-1;z<=Z+1;z++){auto it=mp.find(key(x,y,z)); if(it==mp.end())continue; for(int id:it->second) if(dist2((*pts)[id],q)<=r2)return true;} return false;}
};
static bool coverOK(const vector<Vec3>&X,double slack=0.990){
    if(X.empty()||X.size()>(size_t)max(N,4))return false; double r=H*slack,r2=r*r; GridNear gx; gx.build(X,r);
    for(int i=0;i<N;i++){ if((i&8191)==0&&elapsed()>18.8)return false; if(!gx.near(Orig[i],r2))return false; }
    GridNear go; go.build(Orig,r); for(auto&p:X) if(!go.near(p,r2))return false; return true;
}

static bool validateMesh(const OutMesh&mo){
    int n=mo.X.size(), m=mo.T.size(); if(n<1||m<4)return false; vector<unsigned char> used(n,0); vector<uint64_t>E; E.reserve((size_t)m*3);
    double ae=max(1e-30,diagLen*diagLen*diagLen*diagLen*1e-28);
    for(auto&t:mo.T){int a=t[0],b=t[1],c=t[2]; if(a<0||a>=n||b<0||b>=n||c<0||c>=n)return false; if(a==b||a==c||b==c)return false; Vec3 cr=cross3(mo.X[b]-mo.X[a],mo.X[c]-mo.X[a]); if(norm2(cr)<=ae)return false; used[a]=used[b]=used[c]=1; E.push_back(ekey(a,b));E.push_back(ekey(b,c));E.push_back(ekey(c,a));}
    for(int i=0;i<n;i++)if(!used[i])return false; sort(E.begin(),E.end());
    for(size_t i=0;i<E.size();){size_t j=i+1; while(j<E.size()&&E[j]==E[i])j++; if(j-i!=2)return false; i=j;} return true;
}

struct Img{int R=0; vector<float> z; vector<Vec3> n; vector<unsigned char> fg;};
static inline bool camCoord(const Vec3&p,int view,double D,double&sx,double&sy,double&zz,Vec3&worldNbasis){
    (void)worldNbasis; Vec3 q=p-Center; double x=0,y=0,a=0;
    switch(view){
        case 0: a= q.x; x=-q.y; y= q.z; break; // +x camera
        case 1: a=-q.x; x= q.y; y= q.z; break;
        case 2: a= q.y; x= q.x; y= q.z; break;
        case 3: a=-q.y; x=-q.x; y= q.z; break;
        case 4: a= q.z; x= q.x; y= q.y; break;
        default:a=-q.z; x= q.x; y=-q.y; break;
    }
    zz=D-a; if(zz<=0.05)return false; sx=x; sy=y; return true;
}
static Vec3 orientForView(Vec3 n,int view){
    Vec3 d; switch(view){case 0:d={-1,0,0};break;case 1:d={1,0,0};break;case 2:d={0,-1,0};break;case 3:d={0,1,0};break;case 4:d={0,0,-1};break;default:d={0,0,1};break;} if(dot3(n,d)<0)n=n*-1.0; return normalize3(n);
}
static void rasterTri(Img&im,const Vec3&a,const Vec3&b,const Vec3&c,Vec3 nw,int view){
    double D=2.5, f=0.78125*im.R, ax,ay,az,bx,by,bz,cx,cy,cz; Vec3 dummy;
    if(!camCoord(a,view,D,ax,ay,az,dummy)||!camCoord(b,view,D,bx,by,bz,dummy)||!camCoord(c,view,D,cx,cy,cz,dummy))return;
    ax=f*ax/az+im.R*0.5; ay=f*ay/az+im.R*0.5; bx=f*bx/bz+im.R*0.5; by=f*by/bz+im.R*0.5; cx=f*cx/cz+im.R*0.5; cy=f*cy/cz+im.R*0.5;
    double minx=floor(min(ax,min(bx,cx))), maxx=ceil(max(ax,max(bx,cx))), miny=floor(min(ay,min(by,cy))), maxy=ceil(max(ay,max(by,cy)));
    int x0=max(0,(int)minx), x1=min(im.R-1,(int)maxx), y0=max(0,(int)miny), y1=min(im.R-1,(int)maxy); if(x0>x1||y0>y1)return;
    double den=(by-cy)*(ax-cx)+(cx-bx)*(ay-cy); if(fabs(den)<1e-12)return; nw=orientForView(nw,view);
    for(int y=y0;y<=y1;y++)for(int x=x0;x<=x1;x++){
        double px=x+0.5, py=y+0.5; double w0=((by-cy)*(px-cx)+(cx-bx)*(py-cy))/den; double w1=((cy-ay)*(px-cx)+(ax-cx)*(py-cy))/den; double w2=1.0-w0-w1;
        if(w0<-1e-9||w1<-1e-9||w2<-1e-9)continue; double z=w0*az+w1*bz+w2*cz; int id=y*im.R+x; if(z<im.z[id]){im.z[id]=(float)z; im.n[id]=nw; im.fg[id]=1;}
    }
}
static Img renderMesh(const vector<Vec3>&V,const vector<array<int,3>>&T,int R,int view){
    Img im; im.R=R; im.z.assign(R*R,1e30f); im.n.assign(R*R,Vec3{0,0,0}); im.fg.assign(R*R,0);
    for(auto&t:T){Vec3 a=V[t[0]],b=V[t[1]],c=V[t[2]]; Vec3 nw=cross3(b-a,c-a); if(norm2(nw)<=areaEps)continue; rasterTri(im,a,b,c,nw,view);} return im;
}
static Img renderOriginal(int R,int view){
    Img im; im.R=R; im.z.assign(R*R,1e30f); im.n.assign(R*R,Vec3{0,0,0}); im.fg.assign(R*R,0);
    int step=1; (void)step;
    for(auto&t:F0){Vec3 a=Orig[t[0]],b=Orig[t[1]],c=Orig[t[2]]; Vec3 nw=cross3(b-a,c-a); if(norm2(nw)<=areaEps)continue; rasterTri(im,a,b,c,nw,view);} return im;
}
static double proxyScore(const OutMesh&mo,int R){
    if(elapsed()>19.2)return 0.0; double total=0; int views=0;
    for(int v=0;v<6;v++){
        Img A=renderOriginal(R,v); Img B=renderMesh(mo.X,mo.T,R,v);
        long long ca=0,cb=0,ci=0,cu=0; double de=0,na=0; long long both=0;
        for(int i=0;i<R*R;i++){bool fa=A.fg[i],fb=B.fg[i]; ca+=fa; cb+=fb; ci+=(fa&&fb); cu+=(fa||fb); if(fa&&fb){double dz=fabs((double)A.z[i]-(double)B.z[i])/(0.025*diagLen+1e-9); if(dz>1)dz=1; de+=dz; double nd=dot3(A.n[i],B.n[i]); if(nd<0)nd=0; if(nd>1)nd=1; na+=nd; both++;}}
        if(ca==0&&cb==0){total+=1;views++;continue;} double iou=cu?double(ci)/double(cu):0; double dep=both?1.0-de/both:0; double nor=both?na/both:0; double s=0.46*iou+0.34*dep+0.20*nor; total+=s; views++;
    }
    return views?total/views:0;
}

static void maybeAddCandidate(const string&tag){
    OutMesh mo; mo.tag=tag; vector<int> id(N,-1); mo.T.reserve(activeF); mo.X.reserve(activeV);
    for(int i=0;i<M;i++)if(F[i].alive){int a=F[i].v[0],b=F[i].v[1],c=F[i].v[2]; if(a<0||a>=N||b<0||b>=N||c<0||c>=N||!aliveV[a]||!aliveV[b]||!aliveV[c])return; if(a==b||a==c||b==c)return; int vs[3]={a,b,c}; array<int,3> nt; for(int k=0;k<3;k++){int old=vs[k]; if(id[old]<0){id[old]=(int)mo.X.size(); mo.X.push_back(P[old]);} nt[k]=id[old];} mo.T.push_back(nt);} 
    if(mo.X.empty()||mo.T.empty())return; if(!validateMesh(mo))return; candidates.push_back(std::move(mo));
}

static bool tryBoxGrid(){
    if(N<=20)return false; double tol=diagLen*0.0015; int nearPlane=0; 
    for(auto&p:Orig){double d=min({fabs(p.x-Bmin.x),fabs(p.x-Bmax.x),fabs(p.y-Bmin.y),fabs(p.y-Bmax.y),fabs(p.z-Bmin.z),fabs(p.z-Bmax.z)}); if(d<=tol)nearPlane++;}
    if(nearPlane<N*0.995)return false;
    Vec3 L=Bmax-Bmin; double step=max(H*0.78,diagLen*1e-6); int nx=max(1,(int)ceil(L.x/step)),ny=max(1,(int)ceil(L.y/step)),nz=max(1,(int)ceil(L.z/step));
    if(2*(nx+1)*(ny+1)+2*(nx+1)*(nz+1)+2*(ny+1)*(nz+1)>max(40,N*2))return false;
    map<tuple<int,int,int>,int> mp; OutMesh mo; mo.tag="boxgrid";
    auto add=[&](int ix,int iy,int iz){auto key=make_tuple(ix,iy,iz); auto it=mp.find(key); if(it!=mp.end())return it->second; double x=Bmin.x+L.x*(double)ix/nx, y=Bmin.y+L.y*(double)iy/ny, z=Bmin.z+L.z*(double)iz/nz; int id=mo.X.size(); mp[key]=id; mo.X.push_back({x,y,z}); return id;};
    auto quad=[&](int a,int b,int c,int d,bool flip){array<int,3> t1{a,b,c},t2{a,c,d}; if(flip){swap(t1[1],t1[2]);swap(t2[1],t2[2]);} mo.T.push_back(t1);mo.T.push_back(t2);};
    for(int ix=0;ix<nx;ix++)for(int iy=0;iy<ny;iy++){quad(add(ix,iy,0),add(ix+1,iy,0),add(ix+1,iy+1,0),add(ix,iy+1,0),true); quad(add(ix,iy,nz),add(ix,iy+1,nz),add(ix+1,iy+1,nz),add(ix+1,iy,nz),true);} 
    for(int ix=0;ix<nx;ix++)for(int iz=0;iz<nz;iz++){quad(add(ix,0,iz),add(ix,0,iz+1),add(ix+1,0,iz+1),add(ix+1,0,iz),true); quad(add(ix,ny,iz),add(ix+1,ny,iz),add(ix+1,ny,iz+1),add(ix,ny,iz+1),true);} 
    for(int iy=0;iy<ny;iy++)for(int iz=0;iz<nz;iz++){quad(add(0,iy,iz),add(0,iy+1,iz),add(0,iy+1,iz+1),add(0,iy,iz+1),true); quad(add(nx,iy,iz),add(nx,iy,iz+1),add(nx,iy+1,iz+1),add(nx,iy+1,iz),true);} 
    if(!validateMesh(mo))return false; if(!coverOK(mo.X,0.98))return false; int R=N>200000?64:96; double ps=proxyScore(mo,R); if(ps<0.925)return false; outputMesh(mo); return true;
}

static double ellErrSample(Vec3 c,Vec3 r,double &mx){
    int st=max(1,N/250000); double ss=0; int cnt=0; mx=0; for(int i=0;i<N;i+=st){Vec3 p=Orig[i]; double x=(p.x-c.x)/max(r.x,1e-12),y=(p.y-c.y)/max(r.y,1e-12),z=(p.z-c.z)/max(r.z,1e-12); double e=fabs(sqrt(x*x+y*y+z*z)-1.0); ss+=e*e; mx=max(mx,e); cnt++;} return sqrt(ss/max(1,cnt));
}
static bool makeEllipsoid(int lat,int lon,OutMesh&mo){
    Vec3 c=Center,r=(Bmax-Bmin)/2.0; if(r.x<=1e-12||r.y<=1e-12||r.z<=1e-12)return false; mo.X.clear(); mo.T.clear(); mo.tag="ellipsoid"; double pi=acos(-1.0);
    mo.X.push_back({c.x,c.y,c.z+r.z});
    for(int i=1;i<lat;i++){double th=pi*i/lat,st=sin(th),ct=cos(th); for(int j=0;j<lon;j++){double ph=2*pi*j/lon; mo.X.push_back({c.x+r.x*st*cos(ph),c.y+r.y*st*sin(ph),c.z+r.z*ct});}}
    int bot=mo.X.size(); mo.X.push_back({c.x,c.y,c.z-r.z}); auto id=[&](int ring,int j){return 1+(ring-1)*lon+((j%lon+lon)%lon);};
    for(int j=0;j<lon;j++)mo.T.push_back({0,id(1,j+1),id(1,j)});
    for(int i=1;i<lat-1;i++)for(int j=0;j<lon;j++){int a=id(i,j),b=id(i+1,j),c2=id(i+1,j+1),d=id(i,j+1); mo.T.push_back({a,c2,b}); mo.T.push_back({a,d,c2});}
    for(int j=0;j<lon;j++)mo.T.push_back({bot,id(lat-1,j),id(lat-1,j+1)}); return validateMesh(mo);
}
static bool tryEllipsoid(){
    if(N<400||elapsed()>1.8)return false; Vec3 r=(Bmax-Bmin)/2.0; double anis=max({r.x,r.y,r.z})/max(1e-12,min({r.x,r.y,r.z})); if(anis>3.5)return false; double mx=0,rms=ellErrSample(Center,r,mx); if(!(rms<0.020&&mx<0.080))return false;
    int base=max(24,(int)ceil(2*acos(-1.0)*max({r.x,r.y,r.z})/(H*0.72))); base=min(192,base+(base&1));
    vector<int> lons={base, min(224,(int)(base*1.35)), min(256,(int)(base*1.7))};
    for(int lon:lons){int lat=max(12,lon/2); OutMesh mo; if(!makeEllipsoid(lat,lon,mo))continue; if(mo.X.size()>=N)continue; if(!coverOK(mo.X,0.985))continue; int R=N>200000?64:96; if(proxyScore(mo,R)<0.918)continue; outputMesh(mo); return true; if(elapsed()>7.0)break;}
    return false;
}

static bool collectEdge(int a,int b,int ef[2],int op[2]){
    if(a<0||b<0||a>=N||b>=N)return false; int cnt=0; const vector<int>&A=adj[a],&B=adj[b]; const vector<int>&S=A.size()<B.size()?A:B;
    for(int fid:S){if(fid<0||fid>=M||!F[fid].alive)continue; const Face&f=F[fid]; if(!faceHasEdge(f,a,b))continue; if(cnt>=2)return false; ef[cnt]=fid; op[cnt]=third(f,a,b); cnt++;}
    return cnt==2&&op[0]>=0&&op[1]>=0&&op[0]!=op[1];
}
static bool linkOK(int a,int b,const int op[2]){
    if(++stampv>2000000000){fill(markv.begin(),markv.end(),0);stampv=1;} int s=stampv;
    for(int fid:adj[a])if(F[fid].alive&&hasV(F[fid],a)){auto&f=F[fid]; for(int k=0;k<3;k++){int x=f.v[k]; if(x!=a&&x!=b)markv[x]=s;}}
    int cnt=0,g0=0,g1=0;
    for(int fid:adj[b])if(F[fid].alive&&hasV(F[fid],b)){auto&f=F[fid]; for(int k=0;k<3;k++){int x=f.v[k]; if(x==a||x==b)continue; if(markv[x]==s){if(x==op[0])g0=1; else if(x==op[1])g1=1; else return false; cnt++; markv[x]=s+1; if(cnt>2)return false;}}}
    return cnt==2&&g0&&g1;
}
static bool duplicateAfter(int keep,int rem,int skip0,int skip1,int a,int b,int c){
    if(a==rem)a=keep; if(b==rem)b=keep; if(c==rem)c=keep; if(a==b||a==c||b==c)return true; int pick=a; if(adj[b].size()<adj[pick].size())pick=b; if(adj[c].size()<adj[pick].size())pick=c; auto key=tkey(a,b,c);
    for(int fid:adj[pick]){if(fid==skip0||fid==skip1||!F[fid].alive)continue; auto&f=F[fid]; if(hasV(f,rem))continue; if(tkey(f.v[0],f.v[1],f.v[2])==key)return true;} return false;
}
static void compactAdj(int v){
    if(v<0||v>=N||adj[v].size()<192)return; size_t live=0; for(int fid:adj[v])if(F[fid].alive&&hasV(F[fid],v))live++; if(live*3+64>=adj[v].size())return; vector<int>w; w.reserve(live+8); for(int fid:adj[v])if(F[fid].alive&&hasV(F[fid],v))w.push_back(fid); adj[v].swap(w);
}
struct Cand{bool ok=false; double cost=1e100; Vec3 q; double nr=0,nn=0; int keep=-1,rem=-1;};
static bool qPreservesExtreme(unsigned char m,const Vec3&q){
    double lim=H*0.82; if((m&1)&&fabs(q.x-Bmin.x)>lim)return false; if((m&2)&&fabs(q.x-Bmax.x)>lim)return false; if((m&4)&&fabs(q.y-Bmin.y)>lim)return false; if((m&8)&&fabs(q.y-Bmax.y)>lim)return false; if((m&16)&&fabs(q.z-Bmin.z)>lim)return false; if((m&32)&&fabs(q.z-Bmax.z)>lim)return false; return true;
}
static Cand evalCand(int keep,int rem,const int ef[2],const Vec3&q,double minCos,double planeLim){
    Cand r; r.q=q; r.keep=keep; r.rem=rem; if(!aliveV[keep]||!aliveV[rem])return r; unsigned char mm=emask[keep]|emask[rem]; if(!qPreservesExtreme(mm,q))return r;
    double dk=dist3(P[keep],q), dr=dist3(P[rem],q); double nr=max(radErr[keep]+dk,radErr[rem]+dr); double nn=min(nearErr[keep]+dk,nearErr[rem]+dr); if(nr>H*0.995||nn>H*0.995)return r;
    vector<int> touched; touched.reserve(adj[keep].size()+adj[rem].size());
    for(int fid:adj[keep])if(F[fid].alive&&hasV(F[fid],keep)&&fid!=ef[0]&&fid!=ef[1])touched.push_back(fid);
    for(int fid:adj[rem])if(F[fid].alive&&hasV(F[fid],rem)&&fid!=ef[0]&&fid!=ef[1])touched.push_back(fid);
    sort(touched.begin(),touched.end()); touched.erase(unique(touched.begin(),touched.end()),touched.end()); if(touched.empty())return r;
    double worstN=0,worstP=0,worstA=0; int cnt=0; vector<array<int,3>> made; made.reserve(adj[rem].size());
    for(int fid:touched){Face old=F[fid]; if(hasV(old,keep)&&hasV(old,rem))return r; Face nf=old; for(int k=0;k<3;k++)if(nf.v[k]==rem)nf.v[k]=keep; if(nf.v[0]==nf.v[1]||nf.v[0]==nf.v[2]||nf.v[1]==nf.v[2])return r;
        if(hasV(old,rem)){auto key=tkey(nf.v[0],nf.v[1],nf.v[2]); for(auto&t:made)if(t==key)return r; made.push_back(key); if(duplicateAfter(keep,rem,ef[0],ef[1],old.v[0],old.v[1],old.v[2]))return r;}
        Vec3 oa=P[old.v[0]],ob=P[old.v[1]],oc=P[old.v[2]]; Vec3 na=(old.v[0]==keep||old.v[0]==rem)?q:P[old.v[0]], nb=(old.v[1]==keep||old.v[1]==rem)?q:P[old.v[1]], nc=(old.v[2]==keep||old.v[2]==rem)?q:P[old.v[2]];
        Vec3 no=cross3(ob-oa,oc-oa), nnv=cross3(nb-na,nc-na); double lo=norm3(no), ln=norm3(nnv); if(lo*lo<=areaEps||ln*ln<=areaEps)return r; double co=dot3(no,nnv)/(lo*ln); if(co<minCos)return r;
        Vec3 un=no/lo; double pd=0; if(old.v[0]==keep||old.v[0]==rem)pd=max(pd,fabs(dot3(un,q-oa))); if(old.v[1]==keep||old.v[1]==rem)pd=max(pd,fabs(dot3(un,q-ob))); if(old.v[2]==keep||old.v[2]==rem)pd=max(pd,fabs(dot3(un,q-oc))); if(pd>planeLim)return r;
        worstN=max(worstN,1.0-co); worstP=max(worstP,pd/planeLim); worstA=max(worstA,max(0.0,1.0-ln/lo)); cnt++;
    }
    r.ok=true; r.nr=nr; r.nn=nn; double L=dist3(P[keep],P[rem]); r.cost=L*L*(1+80*worstN+10*worstP+3*worstA)+1e-7*(adj[keep].size()+adj[rem].size())+1e-4*nr/H; return r;
}
static void applyCollapse(const Cand&c,const int ef[2]){
    int keep=c.keep,rem=c.rem; unsigned char oldm=emask[keep], rm=emask[rem]; for(int i=0;i<2;i++)if(F[ef[i]].alive){F[ef[i]].alive=0;activeF--;}
    for(int fid:adj[rem]){if(!F[fid].alive)continue; Face&f=F[fid]; if(!hasV(f,rem))continue; for(int k=0;k<3;k++)if(f.v[k]==rem)f.v[k]=keep; adj[keep].push_back(fid);} 
    aliveV[rem]=0; activeV--; for(int b=0;b<6;b++)if(rm&(1u<<b)){ if(oldm&(1u<<b))extCnt[b]--; }
    emask[keep]=oldm|rm; emask[rem]=0; P[keep]=c.q; radErr[keep]=c.nr; nearErr[keep]=c.nn; ver[keep]++; ver[rem]++; compactAdj(rem); compactAdj(keep);
}
struct Node{double c; int a,b,va,vb; bool operator<(const Node&o)const{return c>o.c;}};
static priority_queue<Node> pq;
static void pushEdge(int a,int b){if(a==b||a<0||b<0||a>=N||b>=N||!aliveV[a]||!aliveV[b])return; double d2=dist2(P[a],P[b]); if(d2>H*H*4.1)return; pq.push({d2,a,b,ver[a],ver[b]});}
static bool attemptCollapse(int a,int b,double minCos,double planeLim){
    if(a==b||!aliveV[a]||!aliveV[b])return false; int ef[2],op[2]; if(!collectEdge(a,b,ef,op)||!linkOK(a,b,op))return false;
    Vec3 qa=P[a], qb=P[b], qm=(P[a]+P[b])*0.5; Cand best;
    Cand ab=evalCand(a,b,ef,qa,minCos,planeLim), ba=evalCand(b,a,ef,qb,minCos,planeLim), am=evalCand(a,b,ef,qm,minCos+0.08,planeLim*0.85), bm=evalCand(b,a,ef,qm,minCos+0.08,planeLim*0.85);
    auto take=[&](const Cand&x){if(x.ok&&(!best.ok||x.cost<best.cost))best=x;}; take(ab);take(ba);take(am);take(bm); if(!best.ok)return false; int keep=best.keep; applyCollapse(best,ef);
    for(int fid:adj[keep])if(F[fid].alive&&hasV(F[fid],keep)){Face&f=F[fid]; pushEdge(f.v[0],f.v[1]);pushEdge(f.v[1],f.v[2]);pushEdge(f.v[2],f.v[0]);}
    return true;
}
static double buildFeatureAndPQ(){
    vector<Vec3> fn(M); for(int i=0;i<M;i++){Vec3 n=faceNormalCur(F[i]); fn[i]=normalize3(n);} vector<pair<uint64_t,int>> E; E.reserve((size_t)M*3);
    for(int i=0;i<M;i++){auto&f=F[i]; E.push_back({ekey(f.v[0],f.v[1]),i});E.push_back({ekey(f.v[1],f.v[2]),i});E.push_back({ekey(f.v[2],f.v[0]),i});}
    sort(E.begin(),E.end(),[](const auto&a,const auto&b){return a.first<b.first;}); long long ue=0,sh=0; double c35=cos(35.0*acos(-1.0)/180.0);
    for(size_t i=0;i<E.size();){size_t j=i+1; while(j<E.size()&&E[j].first==E[i].first)j++; int a=(int)(E[i].first>>32),b=(int)(E[i].first&0xffffffffu); if(j-i==2&&dot3(fn[E[i].second],fn[E[i+1].second])<c35)sh++; ue++; pushEdge(a,b); i=j;} return ue?double(sh)/double(ue):0.0;
}
static void runCollapse(){
    featRatio=buildFeatureAndPQ(); double minCos=featRatio>0.10?0.78:(featRatio>0.035?0.66:0.50); double planeLim=(featRatio>0.10?0.010:(featRatio>0.035?0.0135:0.0175))*diagLen; if(N>160000){minCos-=0.05; planeLim*=1.15;} minCos=max(0.40,minCos);
    vector<int> targets; if(N>1000){targets.push_back(max(40,N/7));targets.push_back(max(48,N/13));targets.push_back(max(64,N/24));targets.push_back(max(80,N/42));}
    size_t nextT=0; long long pops=0, succ=0; double tlim=N>300000?14.2:15.8;
    while(!pq.empty()&&elapsed()<tlim){Node nd=pq.top();pq.pop();pops++; if(nd.a<0||nd.b<0||nd.a>=N||nd.b>=N||!aliveV[nd.a]||!aliveV[nd.b]||ver[nd.a]!=nd.va||ver[nd.b]!=nd.vb)continue; if(attemptCollapse(nd.a,nd.b,minCos,planeLim))succ++;
        while(nextT<targets.size()&&activeV<=targets[nextT]){maybeAddCandidate(string("collapse_")+to_string(targets[nextT])); nextT++; if(elapsed()>tlim)break;}
        if(activeV<=90)break; if((pops&65535)==0&&pq.size()>6000000){priority_queue<Node> npq; swap(pq,npq); vector<uint64_t>E; E.reserve((size_t)activeF*3); for(int i=0;i<M;i++)if(F[i].alive){auto&f=F[i]; E.push_back(ekey(f.v[0],f.v[1]));E.push_back(ekey(f.v[1],f.v[2]));E.push_back(ekey(f.v[2],f.v[0]));} sort(E.begin(),E.end()); E.erase(unique(E.begin(),E.end()),E.end()); for(uint64_t e:E)pushEdge((int)(e>>32),(int)(e&0xffffffffu));}
    }
    maybeAddCandidate("collapse_final");
}

int main(){
    T0=chrono::steady_clock::now(); readInput();
    if(N<=20){outputRawOriginal(); return 0;} // guarantees official cube sample first line 8 12
    if(tryBoxGrid())return 0;
    if(tryEllipsoid())return 0;
    runCollapse();
    sort(candidates.begin(),candidates.end(),[](const OutMesh&a,const OutMesh&b){return a.X.size()<b.X.size();});
    vector<OutMesh> uniq; size_t last=0; for(auto &m:candidates){if(!validateMesh(m))continue; if(!uniq.empty()&&m.X.size()==last)continue; last=m.X.size(); uniq.push_back(std::move(m)); if(uniq.size()>=5)break;}
    int R=N>300000?64:96; double bestScore=-1; int best=-1;
    for(int i=0;i<(int)uniq.size();i++){
        if(elapsed()>19.0)break; double s=proxyScore(uniq[i],R); double thresh=0.905; if(uniq[i].X.size()>(size_t)N/6)thresh=0.895; if(s>bestScore){bestScore=s;best=i;} if(s>=thresh){outputMesh(uniq[i]); return 0;}
    }
    if(best>=0&&N>200000){outputMesh(uniq[best]); return 0;}
    if(best>=0&&bestScore>0.86){outputMesh(uniq[best]); return 0;}
    outputRawOriginal();
    return 0;
}
