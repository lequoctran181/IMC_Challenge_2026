#include <bits/stdc++.h>
#ifndef M_PI
#define M_PI 3.141592653589793238462643383279502884
#endif
using namespace std;

struct Vec3{double x,y,z;};
struct Face{int v[3];};
static inline Vec3 operator+(const Vec3&a,const Vec3&b){return {a.x+b.x,a.y+b.y,a.z+b.z};}
static inline Vec3 operator-(const Vec3&a,const Vec3&b){return {a.x-b.x,a.y-b.y,a.z-b.z};}
static inline Vec3 operator*(const Vec3&a,double s){return {a.x*s,a.y*s,a.z*s};}
static inline double dot3(const Vec3&a,const Vec3&b){return a.x*b.x+a.y*b.y+a.z*b.z;}
static inline Vec3 cross3(const Vec3&a,const Vec3&b){return {a.y*b.z-a.z*b.y,a.z*b.x-a.x*b.z,a.x*b.y-a.y*b.x};}
static inline double norm2(const Vec3&a){return dot3(a,a);} 
static inline double norm3(const Vec3&a){return sqrt(max(0.0,norm2(a)));}
static inline uint64_t ekey(int a,int b){if(a>b)swap(a,b);return (uint64_t)(uint32_t)a<<32 | (uint32_t)b;}
static inline array<int,3> tkey(int a,int b,int c){array<int,3> t{a,b,c};sort(t.begin(),t.end());return t;}

struct FastInput{
    vector<char> b; char *p;
    FastInput(){char buf[1<<16];size_t n;while((n=fread(buf,1,sizeof(buf),stdin))) b.insert(b.end(),buf,buf+n);b.push_back(0);p=b.data();}
    inline void skip(){while(*p==' '||*p=='\n'||*p=='\r'||*p=='\t')++p;}
    inline bool eof(){skip();return *p==0;}
    long long nextInt(){skip(); long long s=1,x=0; if(*p=='-')s=-1,++p; while(*p>='0'&&*p<='9'){x=x*10+(*p-'0');++p;} return s*x;}
    double nextDouble(){skip(); char* q; double x=strtod(p,&q); p=q; return x;}
    void skipLabel(){skip(); if((*p>='a'&&*p<='z')||(*p>='A'&&*p<='Z')){while(*p && *p!=' '&&*p!='\n'&&*p!='\r'&&*p!='\t')++p;}}
};

static int N,M;
static vector<Vec3> P, origP;
static vector<Face> F, origF;
static vector<unsigned char> aliveV, aliveF, fixedV;
static vector<double> coverR;
static vector<vector<int>> adj;
static int aliveCnt, aliveFaces;
static Vec3 bbMin,bbMax,bbExt,bbCtr;
static double diagLen=1.0;
static chrono::steady_clock::time_point T0;
static vector<int> mark1,mark2; static int stamp1=1,stamp2=1;

static inline double elapsed(){return chrono::duration<double>(chrono::steady_clock::now()-T0).count();}
static inline double coord(const Vec3&p,int a){return a==0?p.x:(a==1?p.y:p.z);} 
static inline void setcoord(Vec3&p,int a,double v){if(a==0)p.x=v;else if(a==1)p.y=v;else p.z=v;}
static inline bool contains(const Face&f,int x){return f.v[0]==x||f.v[1]==x||f.v[2]==x;}
static inline bool containsF(int fid,int x){return contains(F[fid],x);} 
static int thirdOf(int fid,int a,int b){const Face&f=F[fid];for(int i=0;i<3;i++){int x=f.v[i];if(x!=a&&x!=b)return x;}return -1;}
static Vec3 faceCross(const vector<Vec3>&V,const Face&f){return cross3(V[f.v[1]]-V[f.v[0]],V[f.v[2]]-V[f.v[0]]);} 
static Vec3 curCross(const Face&f){return cross3(P[f.v[1]]-P[f.v[0]],P[f.v[2]]-P[f.v[0]]);} 
static double signedVol6(const vector<Vec3>&V,const vector<Face>&FF){long double s=0;for(auto &f:FF)s+=(long double)dot3(V[f.v[0]],cross3(V[f.v[1]],V[f.v[2]]));return (double)s;}

static bool readInput(){
    FastInput in; if(in.eof()) return false; N=(int)in.nextInt(); M=(int)in.nextInt();
    P.resize(N); origP.resize(N); aliveV.assign(N,1); coverR.assign(N,0); fixedV.assign(N,0);
    bbMin={1e100,1e100,1e100}; bbMax={-1e100,-1e100,-1e100};
    for(int i=0;i<N;i++){
        in.skipLabel(); P[i].x=in.nextDouble(); P[i].y=in.nextDouble(); P[i].z=in.nextDouble(); origP[i]=P[i];
        bbMin.x=min(bbMin.x,P[i].x); bbMin.y=min(bbMin.y,P[i].y); bbMin.z=min(bbMin.z,P[i].z);
        bbMax.x=max(bbMax.x,P[i].x); bbMax.y=max(bbMax.y,P[i].y); bbMax.z=max(bbMax.z,P[i].z);
    }
    bbExt=bbMax-bbMin; bbCtr=(bbMin+bbMax)*0.5; diagLen=norm3(bbExt); if(!(diagLen>0)) diagLen=1.0;
    F.resize(M); origF.resize(M); aliveF.assign(M,1); vector<int> deg(N,0);
    for(int i=0;i<M;i++){
        in.skipLabel(); int a=(int)in.nextInt()-1,b=(int)in.nextInt()-1,c=(int)in.nextInt()-1;
        if(a<0||b<0||c<0||a>=N||b>=N||c>=N){a=max(0,min(N-1,a));b=max(0,min(N-1,b));c=max(0,min(N-1,c));}
        F[i]={{a,b,c}}; origF[i]=F[i]; deg[a]++; deg[b]++; deg[c]++;
    }
    adj.assign(N,{}); for(int i=0;i<N;i++) adj[i].reserve(deg[i]+4);
    for(int i=0;i<M;i++){adj[F[i].v[0]].push_back(i);adj[F[i].v[1]].push_back(i);adj[F[i].v[2]].push_back(i);}    
    aliveCnt=N; aliveFaces=M; mark1.assign(N,0); mark2.assign(N,0);
    int mnid[3]={0,0,0}, mxid[3]={0,0,0};
    for(int a=0;a<3;a++) for(int i=1;i<N;i++){if(coord(P[i],a)<coord(P[mnid[a]],a))mnid[a]=i;if(coord(P[i],a)>coord(P[mxid[a]],a))mxid[a]=i;}
    for(int a=0;a<3;a++){fixedV[mnid[a]]=1; fixedV[mxid[a]]=1;}
    return true;
}

static void orientLikeOriginal(vector<Vec3>&V, vector<Face>&FF){double s0=signedVol6(origP,origF), s1=signedVol6(V,FF); if(s0*s1<0) for(auto &f:FF) swap(f.v[1],f.v[2]);}
static void printMesh(const vector<Vec3>&V,const vector<Face>&FF){
    printf("%zu %zu\n",V.size(),FF.size());
    for(auto&p:V) printf("v %.17g %.17g %.17g\n",p.x,p.y,p.z);
    for(auto&f:FF) printf("f %d %d %d\n",f.v[0]+1,f.v[1]+1,f.v[2]+1);
}

static long long originalChi(){
    vector<uint64_t> es; es.reserve((size_t)M*3);
    for(auto &f:origF){es.push_back(ekey(f.v[0],f.v[1]));es.push_back(ekey(f.v[1],f.v[2]));es.push_back(ekey(f.v[2],f.v[0]));}
    sort(es.begin(),es.end()); es.erase(unique(es.begin(),es.end()),es.end());
    return (long long)N-(long long)es.size()+(long long)M;
}

static int addGridVertex(map<tuple<int,int,int>,int>&mp, vector<Vec3>&V,int i,int j,int k,int nx,int ny,int nz){
    auto key=make_tuple(i,j,k); auto it=mp.find(key); if(it!=mp.end()) return it->second;
    Vec3 p{bbMin.x+bbExt.x*(double)i/nx,bbMin.y+bbExt.y*(double)j/ny,bbMin.z+bbExt.z*(double)k/nz};
    int id=(int)V.size(); mp[key]=id; V.push_back(p); return id;
}
static bool tryBoxGrid(long long chi){
    if(chi!=2 || N<2000) return false;
    double tol=max(1e-12,diagLen*0.0025); long long bad=0; int cnt[6]={0,0,0,0,0,0};
    for(int i=0;i<N;i++){
        bool ok=false; double d[6]={fabs(P[i].x-bbMin.x),fabs(P[i].x-bbMax.x),fabs(P[i].y-bbMin.y),fabs(P[i].y-bbMax.y),fabs(P[i].z-bbMin.z),fabs(P[i].z-bbMax.z)};
        for(int k=0;k<6;k++) if(d[k]<=tol){ok=true;cnt[k]++;}
        if(!ok) bad++;
    }
    if(bad*1000>N) return false; for(int k=0;k<6;k++) if(cnt[k]<max(4,N/10000)) return false;
    double h=diagLen*0.064; int nx=max(1,(int)ceil(bbExt.x/h)),ny=max(1,(int)ceil(bbExt.y/h)),nz=max(1,(int)ceil(bbExt.z/h));
    if(nx>90||ny>90||nz>90) return false;
    vector<Vec3> V; vector<Face> G; map<tuple<int,int,int>,int> mp;
    auto id=[&](int i,int j,int k){return addGridVertex(mp,V,i,j,k,nx,ny,nz);};
    auto tri=[&](int a,int b,int c){if(a!=b&&a!=c&&b!=c)G.push_back({{a,b,c}});};
    for(int j=0;j<ny;j++)for(int k=0;k<nz;k++){int a=id(0,j,k),b=id(0,j+1,k),c=id(0,j+1,k+1),d=id(0,j,k+1);tri(a,c,b);tri(a,d,c);} 
    for(int j=0;j<ny;j++)for(int k=0;k<nz;k++){int a=id(nx,j,k),b=id(nx,j+1,k),c=id(nx,j+1,k+1),d=id(nx,j,k+1);tri(a,b,c);tri(a,c,d);} 
    for(int i=0;i<nx;i++)for(int k=0;k<nz;k++){int a=id(i,0,k),b=id(i+1,0,k),c=id(i+1,0,k+1),d=id(i,0,k+1);tri(a,b,c);tri(a,c,d);} 
    for(int i=0;i<nx;i++)for(int k=0;k<nz;k++){int a=id(i,ny,k),b=id(i+1,ny,k),c=id(i+1,ny,k+1),d=id(i,ny,k+1);tri(a,c,b);tri(a,d,c);} 
    for(int i=0;i<nx;i++)for(int j=0;j<ny;j++){int a=id(i,j,0),b=id(i+1,j,0),c=id(i+1,j+1,0),d=id(i,j+1,0);tri(a,c,b);tri(a,d,c);} 
    for(int i=0;i<nx;i++)for(int j=0;j<ny;j++){int a=id(i,j,nz),b=id(i+1,j,nz),c=id(i+1,j+1,nz),d=id(i,j+1,nz);tri(a,b,c);tri(a,c,d);} 
    if(V.size()>= (size_t)N || V.empty() || G.empty()) return false;
    orientLikeOriginal(V,G); printMesh(V,G); return true;
}


static Vec3 cylPoint(int axis,int a,int b,double lx,double ly,double lz){Vec3 p=bbCtr; setcoord(p,a,lx); setcoord(p,b,ly); setcoord(p,axis,lz); return p;}
static bool tryCylinder(long long chi){
    if(chi!=2 || N<5000) return false;
    double ext[3]={bbExt.x,bbExt.y,bbExt.z}; int bestAxis=-1; double bestBad=1e100;
    for(int ax=0;ax<3;ax++){
        int a=(ax+1)%3,b=(ax+2)%3; double ra=ext[a]*0.5, rb=ext[b]*0.5, r=(ra+rb)*0.5, L=ext[ax]; if(r<=0||L<=0) continue;
        if(fabs(ra-rb)/max(ra,rb)>0.045) continue;
        double tol=max(1e-12,diagLen*0.010); long long bad=0,side=0,cap0=0,cap1=0; int step=max(1,N/350000);
        for(int i=0;i<N;i+=step){double x=coord(P[i],a)-coord(bbCtr,a),y=coord(P[i],b)-coord(bbCtr,b),z=coord(P[i],ax); double rho=sqrt(x*x+y*y); bool cap=(fabs(z-coord(bbMin,ax))<tol||fabs(z-coord(bbMax,ax))<tol)&&rho<=r+1.5*tol; bool sd=fabs(rho-r)<tol&&z>=coord(bbMin,ax)-tol&&z<=coord(bbMax,ax)+tol; if(!cap&&!sd) bad++; if(sd) side++; if(fabs(z-coord(bbMin,ax))<tol) cap0++; if(fabs(z-coord(bbMax,ax))<tol) cap1++;}
        int cnt=(N+step-1)/step; if(bad*1000>cnt || side<cnt/5 || cap0<4 || cap1<4) continue; if((double)bad<bestBad){bestBad=bad; bestAxis=ax;}
    }
    if(bestAxis<0) return false;
    int ax=bestAxis,a=(ax+1)%3,b=(ax+2)%3; double ra=ext[a]*0.5, rb=ext[b]*0.5, r=(ra+rb)*0.5; double z0=coord(bbMin,ax), z1=coord(bbMax,ax); double ca=coord(bbCtr,a), cb=coord(bbCtr,b);
    double h=diagLen*0.062; int U=max(16,(int)ceil(2*M_PI*r/h)); U=(U+3)/4*4; int H=max(1,(int)ceil((z1-z0)/h)); int Rn=max(1,(int)ceil(r/h));
    if(U>256||H>256||Rn>128) return false; int outN=U*(H+1)+2*(1+max(0,Rn-1)*U); if(outN>=N*0.82) return false;
    vector<Vec3> V; vector<Face> G; vector<vector<int>> side(H+1,vector<int>(U));
    auto sideP=[&](int k,int j){double ph=2*M_PI*j/U;double z=z0+(z1-z0)*(double)k/H;return cylPoint(ax,a,b,ca+r*cos(ph),cb+r*sin(ph),z);};
    for(int k=0;k<=H;k++)for(int j=0;j<U;j++){side[k][j]=V.size();V.push_back(sideP(k,j));}
    for(int k=0;k<H;k++)for(int j=0;j<U;j++){int p=side[k][j],q=side[k][(j+1)%U],r2=side[k+1][(j+1)%U],s=side[k+1][j];G.push_back({{p,q,r2}});G.push_back({{p,r2,s}});} 
    auto buildCap=[&](bool top){int center=V.size();V.push_back(cylPoint(ax,a,b,ca,cb,top?z1:z0)); vector<vector<int>> ring(Rn+1); ring[0].push_back(center); for(int rr=1;rr<Rn;rr++){ring[rr].resize(U); double rad=r*(double)rr/Rn; for(int j=0;j<U;j++){double ph=2*M_PI*j/U; ring[rr][j]=V.size(); V.push_back(cylPoint(ax,a,b,ca+rad*cos(ph),cb+rad*sin(ph),top?z1:z0));}} ring[Rn].resize(U); for(int j=0;j<U;j++) ring[Rn][j]=side[top?H:0][j];
        for(int rr=0;rr<Rn;rr++){
            if(rr==0){for(int j=0;j<U;j++){int o=ring[1][j],on=ring[1][(j+1)%U]; if(top) G.push_back({{center,o,on}}); else G.push_back({{center,on,o}});}}
            else for(int j=0;j<U;j++){int i0=ring[rr][j],i1=ring[rr][(j+1)%U],o0=ring[rr+1][j],o1=ring[rr+1][(j+1)%U]; if(top){G.push_back({{i0,o0,o1}});G.push_back({{i0,o1,i1}});} else {G.push_back({{i0,o1,o0}});G.push_back({{i0,i1,o1}});}}
        }};
    buildCap(false); buildCap(true); orientLikeOriginal(V,G); printMesh(V,G); return true;
}

static bool tryEllipsoid(long long chi){
    if(chi!=2 || N<4000) return false;
    double rx=bbExt.x*0.5, ry=bbExt.y*0.5, rz=bbExt.z*0.5; if(rx<=0||ry<=0||rz<=0) return false;
    long double ss=0; double md=0; int step=max(1,N/250000),cnt=0;
    for(int i=0;i<N;i+=step){double x=(P[i].x-bbCtr.x)/rx,y=(P[i].y-bbCtr.y)/ry,z=(P[i].z-bbCtr.z)/rz;double q=sqrt(x*x+y*y+z*z);double d=fabs(q-1.0);ss+=d*d;md=max(md,d);cnt++;}
    double rms=sqrt((double)(ss/max(1,cnt))); if(rms>0.025||md>0.11) return false;
    double h=diagLen*0.062; int U=max(16,(int)ceil(2*M_PI*max(rx,ry)/h)); U=(U+3)/4*4; int R=max(8,(int)ceil(M_PI*max({rx,ry,rz})/h)); if(R&1)R++;
    if(U>256||R>128) return false; int outN=2+(R-1)*U; if(outN>=N*0.82) return false;
    vector<Vec3> V; vector<Face> G; V.reserve(outN); G.reserve(2*U*R);
    V.push_back({bbCtr.x,bbCtr.y,bbMax.z});
    for(int i=1;i<R;i++){double th=M_PI*i/R, st=sin(th), ct=cos(th); for(int j=0;j<U;j++){double ph=2*M_PI*j/U; V.push_back({bbCtr.x+rx*st*cos(ph),bbCtr.y+ry*st*sin(ph),bbCtr.z+rz*ct});}}
    int bot=(int)V.size(); V.push_back({bbCtr.x,bbCtr.y,bbMin.z});
    auto ring=[&](int i,int j){return 1+(i-1)*U+(j%U);};
    for(int j=0;j<U;j++) G.push_back({{0,ring(1,j),ring(1,j+1)}});
    for(int i=1;i<R-1;i++) for(int j=0;j<U;j++){int a=ring(i,j),b=ring(i+1,j),c=ring(i+1,j+1),d=ring(i,j+1);G.push_back({{a,b,c}});G.push_back({{a,c,d}});} 
    for(int j=0;j<U;j++) G.push_back({{ring(R-1,j),bot,ring(R-1,j+1)}});
    if(V.size()>= (size_t)N) return false; orientLikeOriginal(V,G); printMesh(V,G); return true;
}

static bool tryTorus(long long chi){
    if(chi!=0 || N<5000) return false;
    double ext[3]={bbExt.x,bbExt.y,bbExt.z}; int ax=min_element(ext,ext+3)-ext; int a=(ax+1)%3,b=(ax+2)%3;
    double ra=ext[a]*0.5, rb=ext[b]*0.5, r=ext[ax]*0.5; if(r<=0) return false; if(fabs(ra-rb)/max(ra,rb)>0.05) return false;
    double R=(ra+rb)*0.5-r; if(R<1.45*r) return false;
    long double ss=0; double md=0; int step=max(1,N/250000),cnt=0;
    for(int i=0;i<N;i+=step){double x=coord(P[i],a)-coord(bbCtr,a), y=coord(P[i],b)-coord(bbCtr,b), z=coord(P[i],ax)-coord(bbCtr,ax);double rho=sqrt(x*x+y*y);double q=sqrt((rho-R)*(rho-R)+z*z)/r;double d=fabs(q-1.0);ss+=d*d;md=max(md,d);cnt++;}
    double rms=sqrt((double)(ss/max(1,cnt))); if(rms>0.035||md>0.13) return false;
    double h=diagLen*0.062; int U=max(20,(int)ceil(2*M_PI*(R+r)/h)); U=(U+3)/4*4; int Vn=max(8,(int)ceil(2*M_PI*r/h)); if(Vn&1)Vn++;
    if(U>384||Vn>160) return false; if(U*Vn>=N*0.82) return false;
    vector<Vec3> V; vector<Face> G; V.reserve(U*Vn); G.reserve(2*U*Vn);
    for(int i=0;i<U;i++){double u=2*M_PI*i/U,cu=cos(u),su=sin(u);for(int j=0;j<Vn;j++){double v=2*M_PI*j/Vn,cv=cos(v),sv=sin(v);Vec3 p=bbCtr;setcoord(p,a,coord(bbCtr,a)+(R+r*cv)*cu);setcoord(p,b,coord(bbCtr,b)+(R+r*cv)*su);setcoord(p,ax,coord(bbCtr,ax)+r*sv);V.push_back(p);}}
    auto id=[&](int i,int j){i%=U;j%=Vn;return i*Vn+j;};
    for(int i=0;i<U;i++)for(int j=0;j<Vn;j++){int p=id(i,j),q=id(i+1,j),r2=id(i+1,j+1),s=id(i,j+1);G.push_back({{p,q,r2}});G.push_back({{p,r2,s}});} 
    orientLikeOriginal(V,G); printMesh(V,G); return true;
}

struct Params{double cov,plane,ncos,area;};
struct Cand{bool ok=false; double cost=1e100,rad=0;};

static bool edgeFaces(int u,int v,int ef[2],int op[2]){
    if(adj[u].size()>adj[v].size()) swap(u,v);
    int c=0; for(int fid:adj[u]) if(aliveF[fid]&&containsF(fid,u)&&containsF(fid,v)){ if(c>=2) return false; ef[c]=fid; op[c]=thirdOf(fid,u,v); c++; }
    return c==2 && op[0]>=0 && op[1]>=0 && op[0]!=op[1];
}
static bool linkOK(int u,int v,const int op[2]){
    if(++stamp1>2000000000){fill(mark1.begin(),mark1.end(),0);stamp1=1;} if(++stamp2>2000000000){fill(mark2.begin(),mark2.end(),0);stamp2=1;}
    for(int fid:adj[u]) if(aliveF[fid]&&containsF(fid,u)){Face &f=F[fid];for(int k=0;k<3;k++){int x=f.v[k];if(x!=u&&x!=v) mark1[x]=stamp1;}}
    int inter=0; for(int fid:adj[v]) if(aliveF[fid]&&containsF(fid,v)){Face &f=F[fid];for(int k=0;k<3;k++){int x=f.v[k];if(x==u||x==v||mark1[x]!=stamp1) continue; if(x!=op[0]&&x!=op[1]) return false; if(mark2[x]!=stamp2){mark2[x]=stamp2;inter++;}}}
    return inter==2 && mark2[op[0]]==stamp2 && mark2[op[1]]==stamp2;
}
static bool localDuplicate(int keep,int rem,int e0,int e1,const vector<array<int,3>>&nk){
    vector<array<int,3>> tmp=nk; sort(tmp.begin(),tmp.end()); for(size_t i=1;i<tmp.size();i++) if(tmp[i]==tmp[i-1]) return true;
    for(auto key:nk){int verts[3]={key[0],key[1],key[2]}; int best=verts[0]; if(adj[verts[1]].size()<adj[best].size())best=verts[1]; if(adj[verts[2]].size()<adj[best].size())best=verts[2];
        for(int fid:adj[best]){if(!aliveF[fid]||fid==e0||fid==e1||containsF(fid,rem))continue; if(tkey(F[fid].v[0],F[fid].v[1],F[fid].v[2])==key) return true;}
    }
    return false;
}
static Cand evalCollapse(int keep,int rem,const int ef[2],const Params&p){
    Cand c; if(fixedV[rem]) return c; double d=norm3(P[keep]-P[rem]); c.rad=max(coverR[keep],coverR[rem]+d); if(c.rad>p.cov) return c;
    vector<array<int,3>> newKeys; newKeys.reserve(adj[rem].size()); int changed=0; double worst=0,pen=0;
    for(int fid:adj[rem]) if(aliveF[fid]&&containsF(fid,rem)){
        if(fid==ef[0]||fid==ef[1]) continue; if(containsF(fid,keep)) return c; Face old=F[fid], nf=old; for(int k=0;k<3;k++) if(nf.v[k]==rem) nf.v[k]=keep;
        if(nf.v[0]==nf.v[1]||nf.v[0]==nf.v[2]||nf.v[1]==nf.v[2]) return c;
        Vec3 oc=curCross(old), nc=curCross(nf); double oa=norm3(oc), na=norm3(nc); if(oa<=p.area||na<=p.area) return c; double cs=dot3(oc,nc)/(oa*na); if(cs<p.ncos) return c;
        Vec3 n=oc*(1.0/oa); double pd=fabs(dot3(n,P[keep]-P[old.v[0]])); if(pd>p.plane) return c;
        worst=max(worst,pd); pen+=1.0-cs; newKeys.push_back(tkey(nf.v[0],nf.v[1],nf.v[2])); changed++;
    }
    if(changed==0) return c; if(localDuplicate(keep,rem,ef[0],ef[1],newKeys)) return c;
    c.ok=true; c.cost=d/diagLen + 0.7*c.rad/p.cov + 0.4*worst/max(1e-30,p.plane) + 0.03*pen + 1e-5*adj[rem].size(); return c;
}
static void compactAdj(int v){
    vector<int> w; w.reserve(adj[v].size()); sort(adj[v].begin(),adj[v].end()); adj[v].erase(unique(adj[v].begin(),adj[v].end()),adj[v].end());
    for(int fid:adj[v]) if(aliveF[fid]&&containsF(fid,v)) w.push_back(fid); adj[v].swap(w);
}
static void applyCollapse(int keep,int rem,const int ef[2],double newRad, priority_queue<tuple<double,int,int>,vector<tuple<double,int,int>>,greater<tuple<double,int,int>>>&pq){
    for(int k=0;k<2;k++) if(aliveF[ef[k]]){aliveF[ef[k]]=0; aliveFaces--;}
    for(int fid:adj[rem]) if(aliveF[fid]&&containsF(fid,rem)){for(int k=0;k<3;k++) if(F[fid].v[k]==rem) F[fid].v[k]=keep; adj[keep].push_back(fid);}    
    aliveV[rem]=0; aliveCnt--; coverR[keep]=newRad; adj[rem].clear(); compactAdj(keep);
    for(int fid:adj[keep]) if(aliveF[fid]){Face &f=F[fid]; for(int k=0;k<3;k++){int x=f.v[k]; if(x!=keep&&aliveV[x]) pq.emplace(norm3(P[keep]-P[x]),keep,x);}}
}
static vector<uint64_t> buildAliveEdges(){
    vector<uint64_t> es; es.reserve((size_t)aliveFaces*3);
    for(int i=0;i<M;i++) if(aliveF[i]){Face &f=F[i]; if(aliveV[f.v[0]]&&aliveV[f.v[1]]&&aliveV[f.v[2]]){es.push_back(ekey(f.v[0],f.v[1]));es.push_back(ekey(f.v[1],f.v[2]));es.push_back(ekey(f.v[2],f.v[0]));}}
    sort(es.begin(),es.end()); es.erase(unique(es.begin(),es.end()),es.end()); return es;
}
static void runPhase(const Params&p,double until,int target){
    vector<uint64_t> es=buildAliveEdges(); priority_queue<tuple<double,int,int>,vector<tuple<double,int,int>>,greater<tuple<double,int,int>>> pq;
    for(uint64_t k:es){int u=(int)(k>>32),v=(int)(uint32_t)k; if(aliveV[u]&&aliveV[v]) pq.emplace(norm3(P[u]-P[v]),u,v);} vector<uint64_t>().swap(es);
    while(!pq.empty()&&aliveCnt>target){ if((pq.size()&4095)==0 && elapsed()>until) break; auto [l,u,v]=pq.top(); pq.pop(); if(u==v||!aliveV[u]||!aliveV[v]) continue;
        int ef[2],op[2]; if(!edgeFaces(u,v,ef,op)) continue; if(!linkOK(u,v,op)) continue;
        Cand a=evalCollapse(u,v,ef,p), b=evalCollapse(v,u,ef,p); if(!a.ok&&!b.ok) continue; if(b.ok&&(!a.ok||b.cost<a.cost)) applyCollapse(v,u,ef,b.rad,pq); else applyCollapse(u,v,ef,a.rad,pq);
    }
}
static void fallbackSimplify(){
    double A=max(1e-300,diagLen*diagLen*1e-18); int target=max(12,(int)(N*0.055));
    Params p1{diagLen*0.024,diagLen*0.010,cos(28*M_PI/180.0),A};
    Params p2{diagLen*0.038,diagLen*0.022,cos(55*M_PI/180.0),A};
    Params p3{diagLen*0.0490,diagLen*0.036,0.03,A};
    runPhase(p1,6.2,target); if(elapsed()<18.8&&aliveCnt>target) runPhase(p2,13.2,target); if(elapsed()<19.6&&aliveCnt>target) runPhase(p3,19.7,target);
}
static void outputCurrent(){
    vector<int> id(N,-1); vector<Vec3> V; V.reserve(aliveCnt); for(int i=0;i<N;i++) if(aliveV[i]){id[i]=(int)V.size(); V.push_back(P[i]);}
    vector<Face> G; G.reserve(aliveFaces); for(int i=0;i<M;i++) if(aliveF[i]){int a=id[F[i].v[0]],b=id[F[i].v[1]],c=id[F[i].v[2]]; if(a>=0&&b>=0&&c>=0&&a!=b&&a!=c&&b!=c) G.push_back({{a,b,c}});} 
    orientLikeOriginal(V,G); printMesh(V,G);
}

int main(){
    T0=chrono::steady_clock::now(); if(!readInput()) return 0; long long chi=originalChi();
    if(tryBoxGrid(chi)) return 0;
    if(tryCylinder(chi)) return 0;
    if(tryEllipsoid(chi)) return 0;
    if(tryTorus(chi)) return 0;
    fallbackSimplify(); outputCurrent(); return 0;
}
