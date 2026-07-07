#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <queue>
#include <unordered_map>
#include <utility>
#include <vector>
using namespace std;

struct Vec3{double x,y,z;};
static inline Vec3 operator+(const Vec3&a,const Vec3&b){return{a.x+b.x,a.y+b.y,a.z+b.z};}
static inline Vec3 operator-(const Vec3&a,const Vec3&b){return{a.x-b.x,a.y-b.y,a.z-b.z};}
static inline Vec3 operator*(const Vec3&a,double s){return{a.x*s,a.y*s,a.z*s};}
static inline Vec3 operator/(const Vec3&a,double s){return{a.x/s,a.y/s,a.z/s};}
static inline double dot3(const Vec3&a,const Vec3&b){return a.x*b.x+a.y*b.y+a.z*b.z;}
static inline Vec3 cross3(const Vec3&a,const Vec3&b){return{a.y*b.z-a.z*b.y,a.z*b.x-a.x*b.z,a.x*b.y-a.y*b.x};}
static inline double norm2(const Vec3&a){return dot3(a,a);} 
static inline double dist2(const Vec3&a,const Vec3&b){return norm2(a-b);} 
static inline double norm3(const Vec3&a){return sqrt(norm2(a));}

struct Face{int v[3]; unsigned char alive;};
struct FastInput{
    vector<char> buf; char* p=nullptr;
    FastInput(){char tmp[1<<16]; size_t n; while((n=fread(tmp,1,sizeof(tmp),stdin))>0) buf.insert(buf.end(),tmp,tmp+n); buf.push_back('\0'); p=buf.data();}
    inline void ws(){while(*p==' '||*p=='\n'||*p=='\r'||*p=='\t')++p;}
    int nextInt(){ws(); int s=1; if(*p=='-'){s=-1;++p;} int x=0; while(*p>='0'&&*p<='9'){x=x*10+(*p-'0');++p;} return x*s;}
    double nextDouble(){ws(); char* e; double x=strtod(p,&e); p=e; return x;}
    char nextChar(){ws(); return *p++;}
};

static int N,M,activeV,activeF;
static vector<Vec3> P,Orig;
static vector<Face> F;
static vector<vector<int>> adj;
static vector<unsigned char> valive, emask;
static int extCnt[6];
static Vec3 Bmin,Bmax; static double diagLen=1.0, HT=1.0, HT2=1.0, areaEps=1e-30;
static chrono::steady_clock::time_point T0;
static inline double elapsed(){return chrono::duration<double>(chrono::steady_clock::now()-T0).count();}

static vector<int> headM,tailM,nextM,csz,ver,markv; static int stampv=1;
struct Node{double c; int a,b,va,vb; bool operator<(const Node&o)const{return c>o.c;}};
static priority_queue<Node> pq;

static inline bool hasV(const Face&f,int v){return f.v[0]==v||f.v[1]==v||f.v[2]==v;}
static inline uint64_t ekey(int a,int b){if(a>b)swap(a,b); return (uint64_t)(uint32_t)a<<32 | (uint32_t)b;}
static inline array<int,3> triKey(int a,int b,int c){array<int,3> t{a,b,c}; sort(t.begin(),t.end()); return t;}
static inline Vec3 fnormal(const Face&f){return cross3(P[f.v[1]]-P[f.v[0]],P[f.v[2]]-P[f.v[0]]);} 
static inline double faceArea2(const Face&f){return norm2(fnormal(f));}
static inline int third(const Face&f,int a,int b){for(int k=0;k<3;k++){int x=f.v[k]; if(x!=a&&x!=b)return x;} return -1;}
static bool faceHasEdge(const Face&f,int a,int b){return hasV(f,a)&&hasV(f,b);} 

static void readInput(){
    FastInput in; N=in.nextInt(); M=in.nextInt();
    P.resize(N); Orig.resize(N); F.resize(M);
    Bmin={1e100,1e100,1e100}; Bmax={-1e100,-1e100,-1e100};
    for(int i=0;i<N;i++){
        (void)in.nextChar(); P[i].x=in.nextDouble(); P[i].y=in.nextDouble(); P[i].z=in.nextDouble(); Orig[i]=P[i];
        Bmin.x=min(Bmin.x,P[i].x); Bmin.y=min(Bmin.y,P[i].y); Bmin.z=min(Bmin.z,P[i].z);
        Bmax.x=max(Bmax.x,P[i].x); Bmax.y=max(Bmax.y,P[i].y); Bmax.z=max(Bmax.z,P[i].z);
    }
    Vec3 d=Bmax-Bmin; diagLen=max(1e-12,norm3(d)); HT=0.0490*diagLen; HT2=HT*HT; areaEps=max(1e-30,diagLen*diagLen*diagLen*diagLen*1e-24);
    vector<int> deg(N,0);
    for(int i=0;i<M;i++){
        (void)in.nextChar(); int a=in.nextInt()-1,b=in.nextInt()-1,c=in.nextInt()-1;
        F[i].v[0]=a; F[i].v[1]=b; F[i].v[2]=c; F[i].alive=1;
        if(a>=0&&a<N)deg[a]++; if(b>=0&&b<N)deg[b]++; if(c>=0&&c<N)deg[c]++;
    }
    adj.assign(N,{}); for(int i=0;i<N;i++) adj[i].reserve(deg[i]+6);
    for(int i=0;i<M;i++){adj[F[i].v[0]].push_back(i); adj[F[i].v[1]].push_back(i); adj[F[i].v[2]].push_back(i);} 
    valive.assign(N,1); emask.assign(N,0); memset(extCnt,0,sizeof(extCnt));
    const double eps=diagLen*1e-12+1e-14;
    for(int i=0;i<N;i++){
        unsigned char m=0;
        if(fabs(P[i].x-Bmin.x)<=eps)m|=1; if(fabs(P[i].x-Bmax.x)<=eps)m|=2;
        if(fabs(P[i].y-Bmin.y)<=eps)m|=4; if(fabs(P[i].y-Bmax.y)<=eps)m|=8;
        if(fabs(P[i].z-Bmin.z)<=eps)m|=16; if(fabs(P[i].z-Bmax.z)<=eps)m|=32;
        emask[i]=m; for(int b=0;b<6;b++) if(m&(1u<<b)) extCnt[b]++;
    }
    headM.resize(N); tailM.resize(N); nextM.assign(N,-1); csz.assign(N,1); ver.assign(N,0); markv.assign(N,0);
    for(int i=0;i<N;i++) headM[i]=tailM[i]=i;
    activeV=N; activeF=M;
}

static bool canDeleteExt(int v){unsigned char m=emask[v]; for(int b=0;b<6;b++) if((m&(1u<<b)) && extCnt[b]<=1) return false; return true;}
static bool clusterCanMove(int v,const Vec3&to){for(int m=headM[v];m!=-1;m=nextM[m]) if(dist2(Orig[m],to)>HT2) return false; return true;}

static bool collectEdge(int a,int b,int ef[2],int op[2]){
    int cnt=0; const vector<int>&A=adj[a]; const vector<int>&B=adj[b]; const vector<int>&S=(A.size()<B.size()?A:B);
    for(int fid:S){if(fid<0||fid>=M||!F[fid].alive)continue; const Face&f=F[fid]; if(!faceHasEdge(f,a,b))continue; if(cnt>=2)return false; ef[cnt]=fid; op[cnt]=third(f,a,b); cnt++;}
    if(cnt!=2||op[0]<0||op[1]<0||op[0]==op[1])return false; return true;
}
static bool adjacentNow(int a,int b){int ef[2],op[2]; return collectEdge(a,b,ef,op);}

static bool linkOK(int a,int b,const int op[2]){
    if(++stampv>2000000000){fill(markv.begin(),markv.end(),0); stampv=1;} int s=stampv;
    for(int fid:adj[a]){if(!F[fid].alive||!hasV(F[fid],a))continue; const Face&f=F[fid]; for(int k=0;k<3;k++){int x=f.v[k]; if(x!=a&&x!=b) markv[x]=s;}}
    int cnt=0,g0=0,g1=0;
    for(int fid:adj[b]){if(!F[fid].alive||!hasV(F[fid],b))continue; const Face&f=F[fid]; for(int k=0;k<3;k++){int x=f.v[k]; if(x==a||x==b)continue; if(markv[x]==s){ if(x==op[0])g0=1; else if(x==op[1])g1=1; else return false; cnt++; markv[x]=s+1; if(cnt>2)return false; }}}
    return cnt==2&&g0&&g1;
}

static bool duplicateFaceAfter(int keep,int rem,int skip0,int skip1,int a,int b,int c){
    int x=a,y=b,z=c; if(x==rem)x=keep; if(y==rem)y=keep; if(z==rem)z=keep; if(x==y||x==z||y==z)return true;
    int pick=x; if(adj[y].size()<adj[pick].size())pick=y; if(adj[z].size()<adj[pick].size())pick=z;
    auto key=triKey(x,y,z);
    for(int fid:adj[pick]){
        if(fid==skip0||fid==skip1||!F[fid].alive)continue; const Face&f=F[fid];
        if(hasV(f,rem))continue; if(triKey(f.v[0],f.v[1],f.v[2])==key)return true;
    }
    return false;
}

struct Cand{bool ok=false; double cost=1e100;};
static Cand evalCollapse(int keep,int rem,const int ef[2]){
    Cand r; if(!valive[keep]||!valive[rem]||!canDeleteExt(rem))return r; const Vec3&to=P[keep];
    if(!clusterCanMove(rem,to))return r;
    const double minCos=0.62; double worst=0.0; int touched=0; vector<array<int,3>> made; made.reserve(adj[rem].size());
    for(int fid:adj[rem]){
        if(!F[fid].alive||!hasV(F[fid],rem))continue; if(fid==ef[0]||fid==ef[1])continue; const Face&old=F[fid];
        if(hasV(old,keep))return r;
        Face nf=old; for(int k=0;k<3;k++) if(nf.v[k]==rem) nf.v[k]=keep;
        auto nk=triKey(nf.v[0],nf.v[1],nf.v[2]); for(auto&q:made) if(q==nk)return r; made.push_back(nk);
        if(duplicateFaceAfter(keep,rem,ef[0],ef[1],old.v[0],old.v[1],old.v[2]))return r;
        Vec3 no=fnormal(old); Vec3 nn=cross3(P[nf.v[1]]-P[nf.v[0]],P[nf.v[2]]-P[nf.v[0]]);
        double ao=norm2(no), an=norm2(nn); if(ao<=areaEps||an<=areaEps)return r;
        double co=dot3(no,nn)/sqrt(ao*an); if(co<minCos)return r; worst=max(worst,1.0-co); touched++;
    }
    if(touched==0)return r;
    double L2=dist2(P[keep],P[rem]);
    r.ok=true; r.cost=L2*(1.0+45.0*worst)+0.000001*(double)(adj[keep].size()+adj[rem].size())+0.00001*csz[rem]; return r;
}

static void compactAdj(int v){
    if(v<0||v>=N||adj[v].size()<128)return; size_t live=0; for(int fid:adj[v]) if(F[fid].alive&&hasV(F[fid],v)) live++;
    if(live*3+32>=adj[v].size())return; vector<int> w; w.reserve(live+8); for(int fid:adj[v]) if(F[fid].alive&&hasV(F[fid],v)) w.push_back(fid); adj[v].swap(w);
}
static void mergeMembers(int rem,int keep){
    if(headM[rem]<0)return; nextM[tailM[keep]]=headM[rem]; tailM[keep]=tailM[rem]; csz[keep]+=csz[rem]; headM[rem]=tailM[rem]=-1; csz[rem]=0;
}
static void pushEdge(int a,int b){
    if(a==b||a<0||b<0||a>=N||b>=N||!valive[a]||!valive[b])return; double d=dist2(P[a],P[b]); if(d>HT2*4.0001)return; pq.push({d,a,b,ver[a],ver[b]});
}
static void doCollapse(int keep,int rem,const int ef[2]){
    for(int i=0;i<2;i++) if(F[ef[i]].alive){F[ef[i]].alive=0; activeF--;}
    for(int fid:adj[rem]){
        if(!F[fid].alive)continue; Face&f=F[fid]; if(!hasV(f,rem))continue;
        for(int k=0;k<3;k++) if(f.v[k]==rem) f.v[k]=keep; adj[keep].push_back(fid);
    }
    valive[rem]=0; activeV--; for(int b=0;b<6;b++) if(emask[rem]&(1u<<b)) extCnt[b]--;
    mergeMembers(rem,keep); ver[rem]++; ver[keep]++;
    compactAdj(rem); compactAdj(keep);
    for(int fid:adj[keep]) if(F[fid].alive&&hasV(F[fid],keep)){Face&f=F[fid]; pushEdge(f.v[0],f.v[1]); pushEdge(f.v[1],f.v[2]); pushEdge(f.v[2],f.v[0]);}
}
static bool attempt(int a,int b){
    if(a==b||!valive[a]||!valive[b])return false; int ef[2],op[2]; if(!collectEdge(a,b,ef,op))return false; if(!linkOK(a,b,op))return false;
    Cand ab=evalCollapse(a,b,ef), ba=evalCollapse(b,a,ef); if(!ab.ok&&!ba.ok)return false;
    if(ba.ok&&(!ab.ok||ba.cost<ab.cost))doCollapse(b,a,ef); else doCollapse(a,b,ef); return true;
}

static double buildPQAndFeature(){
    vector<Vec3> fn(M); for(int i=0;i<M;i++){Vec3 n=fnormal(F[i]); double l=norm3(n); fn[i]=(l>0?n/l:Vec3{0,0,0});}
    vector<pair<uint64_t,int>> E; E.reserve((size_t)M*3);
    for(int i=0;i<M;i++){Face&f=F[i]; E.push_back({ekey(f.v[0],f.v[1]),i}); E.push_back({ekey(f.v[1],f.v[2]),i}); E.push_back({ekey(f.v[2],f.v[0]),i});}
    sort(E.begin(),E.end(),[](const auto&a,const auto&b){return a.first<b.first;}); long long ue=0,fe=0; const double fcos=cos(35.0*acos(-1.0)/180.0);
    for(size_t i=0;i<E.size();){size_t j=i+1; while(j<E.size()&&E[j].first==E[i].first)j++; ue++; int a=(int)(E[i].first>>32), b=(int)(E[i].first&0xffffffffu); if(j-i==2){double d=dot3(fn[E[i].second],fn[E[i+1].second]); if(d<fcos)fe++;} pushEdge(a,b); i=j;}
    return ue? (double)fe/(double)ue : 0.0;
}

struct Grid{
    double h; Vec3 mn; unordered_map<long long, vector<int>> mp; const vector<Vec3>* pts;
    static long long key(int x,int y,int z){return ((long long)x*73856093LL)^((long long)y*19349663LL)^((long long)z*83492791LL);} 
    void build(const vector<Vec3>&p,double cell){pts=&p; h=max(cell,1e-12); mn={1e100,1e100,1e100}; for(auto&q:p){mn.x=min(mn.x,q.x); mn.y=min(mn.y,q.y); mn.z=min(mn.z,q.z);} mp.reserve(p.size()*2+10); for(int i=0;i<(int)p.size();i++){int x=(int)floor((p[i].x-mn.x)/h),y=(int)floor((p[i].y-mn.y)/h),z=(int)floor((p[i].z-mn.z)/h); mp[key(x,y,z)].push_back(i);} }
    bool near(const Vec3&q,double r2)const{int X=(int)floor((q.x-mn.x)/h),Y=(int)floor((q.y-mn.y)/h),Z=(int)floor((q.z-mn.z)/h); for(int x=X-1;x<=X+1;x++)for(int y=Y-1;y<=Y+1;y++)for(int z=Z-1;z<=Z+1;z++){auto it=mp.find(key(x,y,z)); if(it==mp.end())continue; for(int id:it->second) if(dist2((*pts)[id],q)<=r2) return true;} return false;}
};
static bool coverOK(const vector<Vec3>&X){
    if(X.empty()||X.size()>(size_t)N)return false; double r=HT*0.995, r2=r*r; Grid gx; gx.build(X,r); for(int i=0;i<N;i++){if((i&8191)==0&&elapsed()>4.0)return false; if(!gx.near(Orig[i],r2))return false;} Grid go; go.build(Orig,r); for(const Vec3&p:X) if(!go.near(p,r2))return false; return true;
}
static double orientSign(const vector<Vec3>&X,const vector<array<int,3>>&T,const Vec3&c){double s=0; int step=max(1,(int)T.size()/100000); for(int i=0;i<(int)T.size();i+=step){auto t=T[i]; Vec3 a=X[t[0]],b=X[t[1]],d=X[t[2]]; Vec3 cr=cross3(b-a,d-a),ct=(a+b+d)/3.0; s+=dot3(cr,ct-c);} return s>=0?1.0:-1.0;}
static void outputMesh(const vector<Vec3>&X, vector<array<int,3>> T){
    printf("%d %d\n",(int)X.size(),(int)T.size());
    for(const Vec3&p:X) printf("v %.17g %.17g %.17g\n",p.x,p.y,p.z);
    for(auto&t:T) printf("f %d %d %d\n",t[0]+1,t[1]+1,t[2]+1);
}
static bool makeEllipsoid(vector<Vec3>&X,vector<array<int,3>>&T,int lat,int lon){
    Vec3 c=(Bmin+Bmax)/2.0, r=(Bmax-Bmin)/2.0; if(r.x<=1e-12||r.y<=1e-12||r.z<=1e-12)return false; X.clear(); T.clear(); const double pi=acos(-1.0);
    X.push_back({c.x,c.y,c.z+r.z});
    for(int i=1;i<lat;i++){double th=pi*i/lat, st=sin(th), ct=cos(th); for(int j=0;j<lon;j++){double ph=2*pi*j/lon; X.push_back({c.x+r.x*st*cos(ph),c.y+r.y*st*sin(ph),c.z+r.z*ct});}}
    int bot=X.size(); X.push_back({c.x,c.y,c.z-r.z}); auto id=[&](int ring,int j){return 1+(ring-1)*lon+((j%lon+lon)%lon);};
    for(int j=0;j<lon;j++)T.push_back({0,id(1,j),id(1,j+1)});
    for(int i=1;i<lat-1;i++)for(int j=0;j<lon;j++){int a=id(i,j),b=id(i+1,j),cc=id(i+1,j+1),d=id(i,j+1); T.push_back({a,b,cc}); T.push_back({a,cc,d});}
    for(int j=0;j<lon;j++)T.push_back({bot,id(lat-1,j+1),id(lat-1,j)});
    if(orientSign(X,T,c)<0) for(auto&t:T) swap(t[1],t[2]); return true;
}

static inline void mapAxis(const Vec3&p,int ax,double&t,double&u,double&v){if(ax==0){t=p.x;u=p.y;v=p.z;}else if(ax==1){t=p.y;u=p.x;v=p.z;}else{t=p.z;u=p.x;v=p.y;}}
static inline Vec3 unmapAxis(int ax,double t,double u,double v){if(ax==0)return{t,u,v}; if(ax==1)return{u,t,v}; return{u,v,t};}
struct TorFit{bool ok=false;int ax=2;double ct,cu,cv,R,r,rms,mx;};
static TorFit fitTorusAxis(int ax){TorFit z;z.ax=ax;if(N<700)return z;double t0=1e100,t1=-1e100,u0=1e100,u1=-1e100,v0=1e100,v1=-1e100;for(const Vec3&p:Orig){double t,u,v;mapAxis(p,ax,t,u,v);t0=min(t0,t);t1=max(t1,t);u0=min(u0,u);u1=max(u1,u);v0=min(v0,v);v1=max(v1,v);}z.ct=(t0+t1)*.5;z.cu=(u0+u1)*.5;z.cv=(v0+v1)*.5;double minrho=1e100,maxrho=0;int st=max(1,N/260000);for(int i=0;i<N;i+=st){double t,u,v;mapAxis(Orig[i],ax,t,u,v);double du=u-z.cu,dv=v-z.cv,rho=sqrt(du*du+dv*dv);minrho=min(minrho,rho);maxrho=max(maxrho,rho);}if(!(maxrho>minrho)&&minrho>0)return z;double R=(maxrho+minrho)*.5, rrho=(maxrho-minrho)*.5, rt=(t1-t0)*.5, r=(rrho+rt)*.5;if(!(R>1e-12)||!(r>1e-12)||R<1.25*r||fabs(rrho-rt)>0.28*r)return z;double ss=0,ma=0;int c=0;for(int i=0;i<N;i+=st){double t,u,v;mapAxis(Orig[i],ax,t,u,v);double du=u-z.cu,dv=v-z.cv,rho=sqrt(du*du+dv*dv);double e=fabs(sqrt((rho-R)*(rho-R)+(t-z.ct)*(t-z.ct))-r)/r;ss+=e*e;ma=max(ma,e);c++;}z.rms=sqrt(ss/max(1,c));z.mx=ma;z.R=R;z.r=r;z.ok=z.rms<0.018&&z.mx<0.075;return z;}
static bool makeTorus(const TorFit&f,int A,int B,vector<Vec3>&X,vector<array<int,3>>&T){if(!f.ok||A<12||B<6)return false;X.clear();T.clear();double pi=acos(-1.0);for(int i=0;i<A;i++){double th=2*pi*i/A,ct=cos(th),st=sin(th);for(int j=0;j<B;j++){double ph=2*pi*j/B,cp=cos(ph),sp=sin(ph),rho=f.R+f.r*cp;X.push_back(unmapAxis(f.ax,f.ct+f.r*sp,f.cu+rho*ct,f.cv+rho*st));}}auto id=[&](int i,int j){return ((i%A+A)%A)*B+((j%B+B)%B);};auto pred=[&](const Vec3&p){double t,u,v;mapAxis(p,f.ax,t,u,v);double du=u-f.cu,dv=v-f.cv,rho=sqrt(du*du+dv*dv);if(rho<1e-12)return Vec3{0,0,0};return unmapAxis(f.ax,(t-f.ct)/f.r,(rho-f.R)/f.r*du/rho,(rho-f.R)/f.r*dv/rho);};double os=0;int oc=0,step=max(1,M/100000);for(int q=0;q<M;q+=step){Face &ff=F[q];Vec3 a=Orig[ff.v[0]],b=Orig[ff.v[1]],c=Orig[ff.v[2]],cr=cross3(b-a,c-a),ce=(a+b+c)/3.0,pr=pred(ce);double l=norm3(cr)*norm3(pr);if(l>1e-20){os+=dot3(cr,pr)/l;oc++;}}bool flip=oc&&os<0;for(int i=0;i<A;i++)for(int j=0;j<B;j++){array<int,3>a{id(i,j),id(i+1,j),id(i+1,j+1)},b{id(i,j),id(i+1,j+1),id(i,j+1)};if(flip)swap(a[1],a[2]),swap(b[1],b[2]);T.push_back(a);T.push_back(b);}return true;}
static bool tryTorus(){if(N<700||elapsed()>2.0)return false;TorFit best;for(int ax=0;ax<3;ax++){TorFit f=fitTorusAxis(ax);if(f.ok&&(!best.ok||f.rms<best.rms))best=f;}if(!best.ok)return false;int tr[4][2]={{48,12},{64,16},{80,20},{96,24}};vector<Vec3>X;vector<array<int,3>>T;for(auto&a:tr){if(a[0]*a[1]>=N)continue;if(!makeTorus(best,a[0],a[1],X,T))continue;if(coverOK(X)){outputMesh(X,T);return true;}if(elapsed()>5.3)break;}return false;}
struct RevFit{bool ok=false;int ax=2;double cu,cv,t0,t1,r0,r1,rms,mx;};
static RevFit fitRevAxis(int ax){RevFit f;f.ax=ax;if(N<700)return f;double t0=1e100,t1=-1e100,u0=1e100,u1=-1e100,v0=1e100,v1=-1e100;for(const Vec3&p:Orig){double t,u,v;mapAxis(p,ax,t,u,v);t0=min(t0,t);t1=max(t1,t);u0=min(u0,u);u1=max(u1,u);v0=min(v0,v);v1=max(v1,v);}if(!(t1>t0))return f;f.t0=t0;f.t1=t1;f.cu=(u0+u1)*.5;f.cv=(v0+v1)*.5;int st=max(1,N/240000),cnt=0,axisPts=0;double S=0,St=0,Stt=0,Sr=0,Str=0,maxr=0;for(int i=0;i<N;i+=st){double t,u,v;mapAxis(Orig[i],ax,t,u,v);double r=sqrt((u-f.cu)*(u-f.cu)+(v-f.cv)*(v-f.cv));maxr=max(maxr,r);if(r<1e-9)continue;S++;St+=t;Stt+=t*t;Sr+=r;Str+=t*r;cnt++;}if(cnt<200||maxr<1e-12)return f;double den=S*Stt-St*St;if(fabs(den)<1e-18)return f;double a=(S*Str-St*Sr)/den,b=(Sr-a*St)/S;f.r0=max(0.0,a*t0+b);f.r1=max(0.0,a*t1+b);if(max(f.r0,f.r1)<maxr*.25)return f;double ss=0,ma=0;int c=0;double len=t1-t0;for(int i=0;i<N;i+=st){double t,u,v;mapAxis(Orig[i],ax,t,u,v);double r=sqrt((u-f.cu)*(u-f.cu)+(v-f.cv)*(v-f.cv));if(r<maxr*.04){if(min(fabs(t-t0),fabs(t-t1))>len*.035)axisPts+=5;continue;}double pr=max(0.0,a*t+b);double e=fabs(r-pr)/maxr;ss+=e*e;ma=max(ma,e);c++;}if(axisPts>cnt/3||c<200)return f;f.rms=sqrt(ss/max(1,c));f.mx=ma;f.ok=f.rms<0.010&&f.mx<0.045;return f;}
static bool makeRev(const RevFit&f,int L,int S,vector<Vec3>&X,vector<array<int,3>>&T){
    if(!f.ok||L<2||S<12)return false;X.clear();T.clear();double pi=acos(-1.0),rmax=max(f.r0,f.r1);auto rad=[&](int i){double ss=(double)i/(L-1);return f.r0+(f.r1-f.r0)*ss;};vector<int>ring(L,-1),apex(L,-1);for(int i=0;i<L;i++){double t=f.t0+(f.t1-f.t0)*(double)i/(L-1),r=rad(i);if(r<rmax*1e-6){apex[i]=X.size();X.push_back(unmapAxis(f.ax,t,f.cu,f.cv));}else{ring[i]=X.size();for(int j=0;j<S;j++){double th=2*pi*j/S;X.push_back(unmapAxis(f.ax,t,f.cu+r*cos(th),f.cv+r*sin(th)));}}}
    auto id=[&](int i,int j){return ring[i]+((j%S+S)%S);};
    for(int i=0;i+1<L;i++)for(int j=0;j<S;j++){if(ring[i]>=0&&ring[i+1]>=0){T.push_back({id(i,j),id(i,j+1),id(i+1,j)});T.push_back({id(i,j+1),id(i+1,j+1),id(i+1,j)});}else if(apex[i]>=0&&ring[i+1]>=0){T.push_back({apex[i],id(i+1,j+1),id(i+1,j)});}else if(ring[i]>=0&&apex[i+1]>=0){T.push_back({apex[i+1],id(i,j),id(i,j+1)});}}
    if(ring[0]>=0){int c=X.size();X.push_back(unmapAxis(f.ax,f.t0,f.cu,f.cv));for(int j=0;j<S;j++)T.push_back({c,id(0,j),id(0,j+1)});}if(ring[L-1]>=0){int c=X.size();X.push_back(unmapAxis(f.ax,f.t1,f.cu,f.cv));for(int j=0;j<S;j++)T.push_back({c,id(L-1,j+1),id(L-1,j)});}Vec3 cen=unmapAxis(f.ax,(f.t0+f.t1)*.5,f.cu,f.cv);if(orientSign(X,T,cen)<0)for(auto&t:T)swap(t[1],t[2]);return !X.empty()&&!T.empty();}
static bool tryRevolve(){if(N<700||elapsed()>2.0)return false;RevFit best;for(int ax=0;ax<3;ax++){RevFit f=fitRevAxis(ax);if(f.ok&&(!best.ok||f.rms<best.rms))best=f;}if(!best.ok)return false;int tr[4][2]={{12,32},{18,40},{24,48},{32,64}};vector<Vec3>X;vector<array<int,3>>T;for(auto&a:tr){if(a[0]*a[1]+2>=N)continue;if(!makeRev(best,a[0],a[1],X,T))continue;if(coverOK(X)){outputMesh(X,T);return true;}if(elapsed()>5.3)break;}return false;}

static bool tryEllipsoid(){
    if(N<900||elapsed()>2.0)return false; Vec3 c=(Bmin+Bmax)/2.0, r=(Bmax-Bmin)/2.0; double hi=max(r.x,max(r.y,r.z)), lo=min(r.x,min(r.y,r.z)); if(lo<=1e-12||hi/lo>2.8)return false;
    int stride=max(1,N/260000); double ss=0,ma=0; int cnt=0; for(int i=0;i<N;i+=stride){Vec3 q=Orig[i]-c; double rr=sqrt((q.x*q.x)/(r.x*r.x)+(q.y*q.y)/(r.y*r.y)+(q.z*q.z)/(r.z*r.z)); double e=fabs(rr-1.0); ss+=e*e; ma=max(ma,e); cnt++;}
    double rms=sqrt(ss/max(1,cnt)); if(rms>0.010||ma>0.050)return false;
    vector<Vec3>X; vector<array<int,3>>T; int trials[7][2]={{12,24},{14,28},{16,32},{18,36},{20,40},{24,48},{28,56}};
    for(auto &tr:trials){ if(2+(tr[0]-1)*tr[1]>=N)continue; if(!makeEllipsoid(X,T,tr[0],tr[1]))continue; if(coverOK(X)){outputMesh(X,T); return true;} if(elapsed()>5.2)break; }
    return false;
}

static void simplifyGeneric(){
    double fr=buildPQAndFeature(); double ratio=0.145+0.11*min(0.5,fr); if(N<2000)ratio=max(ratio,0.24); else if(N<10000)ratio=max(ratio,0.19); else if(N>200000&&fr<0.03)ratio=0.125; ratio=max(0.115,min(0.34,ratio)); int target=max(8,(int)ceil(N*ratio));
    long long pops=0,success=0; while(activeV>target&&!pq.empty()){
        if((++pops&4095)==0 && elapsed()>19.35)break; Node cur=pq.top(); pq.pop(); int a=cur.a,b=cur.b; if(a<0||b<0||a>=N||b>=N||!valive[a]||!valive[b])continue; if(cur.va!=ver[a]||cur.vb!=ver[b]){if(adjacentNow(a,b))pushEdge(a,b); continue;} if(attempt(a,b))success++; }
}
static void outputCurrent(){
    vector<int> remap(N,-1); int nv=0; for(int i=0;i<N;i++) if(valive[i]) remap[i]=nv++;
    vector<array<int,3>> T; T.reserve(activeF); for(int i=0;i<M;i++) if(F[i].alive){int a=F[i].v[0],b=F[i].v[1],c=F[i].v[2]; if(a==b||a==c||b==c)continue; if(a<0||b<0||c<0||a>=N||b>=N||c>=N)continue; if(remap[a]<0||remap[b]<0||remap[c]<0)continue; T.push_back({remap[a],remap[b],remap[c]});}
    printf("%d %d\n",nv,(int)T.size());
    for(int i=0;i<N;i++) if(remap[i]>=0) printf("v %.17g %.17g %.17g\n",P[i].x,P[i].y,P[i].z);
    for(auto&t:T) printf("f %d %d %d\n",t[0]+1,t[1]+1,t[2]+1);
}
int main(){T0=chrono::steady_clock::now(); readInput(); if(tryTorus())return 0; if(tryRevolve())return 0; if(tryEllipsoid())return 0; simplifyGeneric(); outputCurrent(); return 0;}
