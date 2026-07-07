C++

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
static inline double dist2(const Vec3&a,const Vec3&b){return norm2(a-b);} 
static inline double clampd(double x,double a,double b){return x<a?a:(x>b?b:x);} 
static inline bool normalize(Vec3&v){ double n=norm3(v); if(n<=1e-300) return false; v=v*(1.0/n); return true; }

struct Face{ int v[3]; };
struct Mesh{ vector<Vec3> V; vector<Face> F; string route; };

static int N=0,M=0;
static vector<Vec3> OrigV;
static vector<Face> OrigF;
static Vec3 BMin,BMax,BCenter;
static double Diag=1.0, Eps=1.0, Eps2=1.0, AreaEps2=1e-30;
static chrono::steady_clock::time_point T0;
static inline double elapsed(){return chrono::duration<double>(chrono::steady_clock::now()-T0).count();}

struct FastInput{
    vector<char> buf; char* p=nullptr;
    FastInput(){ char tmp[1<<16]; size_t n; while((n=fread(tmp,1,sizeof(tmp),stdin))>0) buf.insert(buf.end(),tmp,tmp+n); buf.push_back('\0'); p=buf.data(); }
    inline void skip(){ while(*p==' '||*p=='\n'||*p=='\r'||*p=='\t') ++p; }
    int nextInt(){ skip(); int s=1; if(*p=='-'){s=-1;++p;} int x=0; while(*p>='0'&&*p<='9'){x=x*10+(*p-'0');++p;} return x*s; }
    double nextDouble(){ skip(); char* e=nullptr; double x=strtod(p,&e); p=e; return x; }
    char nextChar(){ skip(); return *p++; }
    bool eof(){ skip(); return *p=='\0'; }
};

static bool readInput(){
    FastInput in; if(in.eof()) return false; N=in.nextInt(); M=in.nextInt(); if(N<=0||M<=0) return false;
    OrigV.resize(N); OrigF.resize(M); BMin={1e100,1e100,1e100}; BMax={-1e100,-1e100,-1e100};
    for(int i=0;i<N;i++){
        char c=in.nextChar(); if(c!='v'&&c!='V') --in.p;
        Vec3 p{in.nextDouble(),in.nextDouble(),in.nextDouble()}; OrigV[i]=p;
        BMin.x=min(BMin.x,p.x); BMin.y=min(BMin.y,p.y); BMin.z=min(BMin.z,p.z);
        BMax.x=max(BMax.x,p.x); BMax.y=max(BMax.y,p.y); BMax.z=max(BMax.z,p.z);
    }
    for(int i=0;i<M;i++){
        char c=in.nextChar(); if(c!='f'&&c!='F') --in.p;
        OrigF[i].v[0]=in.nextInt()-1; OrigF[i].v[1]=in.nextInt()-1; OrigF[i].v[2]=in.nextInt()-1;
    }
    BCenter=(BMin+BMax)*0.5; Diag=norm3(BMax-BMin); if(!(Diag>0)) Diag=1.0; Eps=0.0493*Diag; Eps2=Eps*Eps; AreaEps2=max(1e-30,1e-24*Diag*Diag*Diag*Diag); return true;
}

static inline uint64_t edgeKey(int a,int b){ if(a>b) swap(a,b); return (uint64_t)(uint32_t)a<<32 | (uint32_t)b; }
static inline array<int,3> triKey(int a,int b,int c){ array<int,3> t{a,b,c}; sort(t.begin(),t.end()); return t; }
static double signedVol6(const vector<Vec3>&V,const vector<Face>&F){ long double s=0; for(auto &f:F){ if(f.v[0]<0||f.v[1]<0||f.v[2]<0||f.v[0]>=(int)V.size()||f.v[1]>=(int)V.size()||f.v[2]>=(int)V.size()) continue; s += (long double)dot3(V[f.v[0]], cross3(V[f.v[1]], V[f.v[2]])); } return (double)s; }
static void matchOriginalSign(Mesh&m){ double a=signedVol6(OrigV,OrigF), b=signedVol6(m.V,m.F); if(a*b<0) for(auto &f:m.F) swap(f.v[1],f.v[2]); }

struct SpatialHash{
    double cell=1, inv=1; unordered_map<long long, vector<int>> H;
    static uint64_t mix(uint64_t x){ x+=0x9e3779b97f4a7c15ULL; x=(x^(x>>30))*0xbf58476d1ce4e5b9ULL; x=(x^(x>>27))*0x94d049bb133111ebULL; return x^(x>>31); }
    static long long key(int x,int y,int z){ uint64_t h=mix((uint64_t)(int64_t)x); h^=mix((uint64_t)(int64_t)y)+0x9e3779b97f4a7c15ULL+(h<<6)+(h>>2); h^=mix((uint64_t)(int64_t)z)+0x9e3779b97f4a7c15ULL+(h<<6)+(h>>2); return (long long)h; }
    void coord(const Vec3&p,int&x,int&y,int&z)const{ x=(int)floor(p.x*inv); y=(int)floor(p.y*inv); z=(int)floor(p.z*inv); }
    void build(const vector<Vec3>&P,double c){ cell=max(c,1e-300); inv=1.0/cell; H.clear(); H.reserve(P.size()*2+1024); for(int i=0;i<(int)P.size();++i){int x,y,z; coord(P[i],x,y,z); H[key(x,y,z)].push_back(i);} }
    bool nearAny(const vector<Vec3>&P,const Vec3&q,double r2,int rad=1)const{
        int x,y,z; coord(q,x,y,z);
        for(int dx=-rad;dx<=rad;dx++) for(int dy=-rad;dy<=rad;dy++) for(int dz=-rad;dz<=rad;dz++){
            auto it=H.find(key(x+dx,y+dy,z+dz)); if(it==H.end()) continue;
            for(int id:it->second) if(dist2(P[id],q)<=r2) return true;
        }
        return false;
    }
    int nearest(const vector<Vec3>&P,const Vec3&q,int rad=2)const{
        int x,y,z; coord(q,x,y,z); int best=-1; double bd=1e300;
        for(int R=0;R<=rad;R++){
            for(int dx=-R;dx<=R;dx++) for(int dy=-R;dy<=R;dy++) for(int dz=-R;dz<=R;dz++){
                if(max({abs(dx),abs(dy),abs(dz)})!=R) continue;
                auto it=H.find(key(x+dx,y+dy,z+dz)); if(it==H.end()) continue;
                for(int id:it->second){ double d=dist2(P[id],q); if(d<bd){bd=d; best=id;} }
            }
            if(best>=0 && bd <= (R+1.5)*(R+1.5)*cell*cell) break;
        }
        return best;
    }
};

static bool allNear(const vector<Vec3>&A,const vector<Vec3>&B,double eps,bool exact=true,int sampleMax=250000){
    if(A.empty()||B.empty()) return false; SpatialHash gh; gh.build(B,eps); double e2=eps*eps; int stride=1; if(!exact && (int)A.size()>sampleMax) stride=max(1,(int)A.size()/sampleMax); for(int i=0;i<(int)A.size();i+=stride) if(!gh.nearAny(B,A[i],e2,1)) return false; return true;
}
static bool vertexHausdorffOK(const vector<Vec3>&V,bool exact=true){ return allNear(V,OrigV,0.05005*Diag,true) && allNear(OrigV,V,0.05005*Diag,exact); }

static bool closedManifoldOK(const Mesh&m){
    if(m.V.empty()||m.F.empty()||m.V.size()>(size_t)max(N+1000, N*2)) return false; vector<uint64_t> edges; edges.reserve(m.F.size()*3); unordered_set<array<int,3>, function<size_t(const array<int,3>&)>> seen(0, [](const array<int,3>&a){ return ((uint64_t)a[0]*1000003u)^((uint64_t)a[1]*9176u)^a[2]; });
    seen.reserve(m.F.size()*2+1);
    for(auto &f:m.F){
        int a=f.v[0],b=f.v[1],c=f.v[2]; if(a<0||b<0||c<0||a>=(int)m.V.size()||b>=(int)m.V.size()||c>=(int)m.V.size()) return false; if(a==b||a==c||b==c) return false;
        Vec3 cr=cross3(m.V[b]-m.V[a],m.V[c]-m.V[a]); if(!(norm2(cr)>AreaEps2)) return false; auto tk=triKey(a,b,c); if(seen.find(tk)!=seen.end()) return false; seen.insert(tk); edges.push_back(edgeKey(a,b)); edges.push_back(edgeKey(b,c)); edges.push_back(edgeKey(c,a));
    }
    sort(edges.begin(),edges.end()); for(size_t i=0;i<edges.size();){ size_t j=i+1; while(j<edges.size()&&edges[j]==edges[i]) ++j; if(j-i!=2) return false; i=j; } return true;
}

static void addOriented(Mesh&m,int a,int b,int c,const Vec3&center){ Face f{{a,b,c}}; Vec3 cr=cross3(m.V[b]-m.V[a],m.V[c]-m.V[a]); Vec3 ce=(m.V[a]+m.V[b]+m.V[c])/3.0; if(dot3(cr,ce-center)<0) swap(f.v[1],f.v[2]); m.F.push_back(f); }

struct RenderMap{ int r=0; vector<float> z; vector<Vec3> n; vector<unsigned char> mask; void init(int R){r=R; z.assign(r*r,1e30f); n.assign(r*r,{}); mask.assign(r*r,0);} };
static void projectPoint(const Vec3&p,int view,int res,double&x,double&y,double&z){ Vec3 q=(p-BCenter)*(2.1/Diag); if(view==0){x=q.y;y=q.z;z=-q.x;} else if(view==1){x=-q.y;y=q.z;z=q.x;} else if(view==2){x=-q.x;y=q.z;z=-q.y;} else if(view==3){x=q.x;y=q.z;z=q.y;} else if(view==4){x=q.x;y=q.y;z=-q.z;} else {x=-q.x;y=q.y;z=q.z;} x=(x*0.45+0.5)*res; y=(y*0.45+0.5)*res; }
static void rasterTri(RenderMap&rm,const vector<Vec3>&V,const Face&f,int view){
    double x0,y0,z0,x1,y1,z1,x2,y2,z2; projectPoint(V[f.v[0]],view,rm.r,x0,y0,z0); projectPoint(V[f.v[1]],view,rm.r,x1,y1,z1); projectPoint(V[f.v[2]],view,rm.r,x2,y2,z2);
    int xmin=max(0,(int)floor(min({x0,x1,x2}))), xmax=min(rm.r-1,(int)ceil(max({x0,x1,x2}))); int ymin=max(0,(int)floor(min({y0,y1,y2}))), ymax=min(rm.r-1,(int)ceil(max({y0,y1,y2}))); if(xmin>xmax||ymin>ymax) return; double den=(y1-y2)*(x0-x2)+(x2-x1)*(y0-y2); if(fabs(den)<1e-18) return; Vec3 no=cross3(V[f.v[1]]-V[f.v[0]],V[f.v[2]]-V[f.v[0]]); if(!normalize(no)) return;
    for(int yy=ymin;yy<=ymax;yy++){ double py=yy+0.5; for(int xx=xmin;xx<=xmax;xx++){ double px=xx+0.5; double w0=((y1-y2)*(px-x2)+(x2-x1)*(py-y2))/den; double w1=((y2-y0)*(px-x2)+(x0-x2)*(py-y2))/den; double w2=1-w0-w1; if(w0<-1e-9||w1<-1e-9||w2<-1e-9) continue; double zz=w0*z0+w1*z1+w2*z2; int id=yy*rm.r+xx; if(zz<rm.z[id]){rm.z[id]=(float)zz; rm.n[id]=no; rm.mask[id]=1;}} }
}
static RenderMap renderMesh(const vector<Vec3>&V,const vector<Face>&F,int view,int res){ RenderMap rm; rm.init(res); for(auto &f:F) rasterTri(rm,V,f,view); return rm; }
static double proxyScore(const Mesh&m,int res=96){
    double total=0; for(int view=0; view<6; ++view){ RenderMap A=renderMesh(OrigV,OrigF,view,res), B=renderMesh(m.V,m.F,view,res); double inter=0,uni=0,ns=0,ds=0; int cnt=0; for(int i=0;i<res*res;i++){ if(A.mask[i]||B.mask[i]){uni++; if(A.mask[i]&&B.mask[i]){inter++; ns+=max(0.0,dot3(A.n[i],B.n[i])); ds+=exp(-8.0*fabs((double)A.z[i]-(double)B.z[i])); cnt++;}} } double iou=uni?inter/uni:1.0; double nsc=cnt?ns/cnt:0.0; double dsc=cnt?ds/cnt:0.0; total += 0.45*iou+0.35*nsc+0.20*dsc; } return total/6.0; }

static bool acceptCandidate(Mesh&m,bool useProxy=false,double pthr=0.88){ if(m.V.empty()||m.F.empty()||m.V.size()>=(size_t)N) return false; matchOriginalSign(m); if(!closedManifoldOK(m)) return false; if(!vertexHausdorffOK(m.V,true)) return false; if(useProxy){ int res=(N>120000?80:96); if(proxyScore(m,res)<pthr) return false; } return true; }

static Mesh makeBoxCandidate(){
    Mesh m; m.route="box"; Vec3 mn=BMin,mx=BMax; m.V={{mn.x,mn.y,mn.z},{mx.x,mn.y,mn.z},{mx.x,mx.y,mn.z},{mn.x,mx.y,mn.z},{mn.x,mn.y,mx.z},{mx.x,mn.y,mx.z},{mx.x,mx.y,mx.z},{mn.x,mx.y,mx.z}}; int t[12][3]={{0,1,2},{0,2,3},{4,6,5},{4,7,6},{0,4,5},{0,5,1},{1,5,6},{1,6,2},{2,6,7},{2,7,3},{3,7,4},{3,4,0}}; for(auto &q:t)m.F.push_back({{q[0],q[1],q[2]}}); return m;
}

static bool unorderedFaceEq(const Face&f,int a,int b,int c){ return triKey(f.v[0],f.v[1],f.v[2])==triKey(a,b,c); }
static bool detectLatGrid(int&R,int&S){
    if(N<10 || M != 2*(N-2)) return false;
    for(int s=8;s<=4096;s++){ if((N-2)%s) continue; int r=(N-2)/s; if(r<3) continue; bool ok=true; int step=max(1,s/96);
        for(int j=0;j<s&&ok;j+=step){ int a=1+j,b=1+(j+1)%s; if(!unorderedFaceEq(OrigF[j],0,b,a)) ok=false; }
        int body=max(1,(r-1)*s/256); for(int q=0;q<(r-1)*s&&ok;q+=body){ int rr=q/s,j=q-rr*s; int a=1+rr*s+j,b=1+rr*s+(j+1)%s,c=1+(rr+1)*s+j,d=1+(rr+1)*s+(j+1)%s; int fid=s+2*q; if(fid+1>=M){ok=false;break;} if(!unorderedFaceEq(OrigF[fid],a,b,c)) ok=false; if(ok&&!unorderedFaceEq(OrigF[fid+1],b,d,c)) ok=false; }
        int bot=s+2*(r-1)*s; for(int j=0;j<s&&ok;j+=step){ int c=1+(r-1)*s+j,d=1+(r-1)*s+(j+1)%s; if(bot+j>=M || !unorderedFaceEq(OrigF[bot+j],N-1,c,d)) ok=false; }
        if(ok){R=r;S=s;return true;}
    } return false;
}
static Mesh resampleLat(int R,int S,int R2,int S2){
    Mesh m; m.route="lat"; if(R2<3||S2<8||2+R2*S2>=N) return m; m.V.reserve(2+R2*S2); m.F.reserve(2*R2*S2); m.V.push_back(OrigV[0]);
    for(int r=1;r<=R2;r++){ int oring=1+(int)llround((double)(r-1)*(R-1)/max(1,R2-1)); oring=max(1,min(R,oring)); for(int j=0;j<S2;j++){ int oj=(int)((long long)j*S/S2)%S; m.V.push_back(OrigV[1+(oring-1)*S+oj]); }} int bottom=m.V.size(); m.V.push_back(OrigV[N-1]); auto id=[&](int r,int j){return 1+(r-1)*S2+((j%S2+S2)%S2);};
    for(int j=0;j<S2;j++) addOriented(m,0,id(1,j+1),id(1,j),BCenter); for(int r=1;r<R2;r++) for(int j=0;j<S2;j++){int a=id(r,j),b=id(r,j+1),c=id(r+1,j),d=id(r+1,j+1); addOriented(m,a,b,c,BCenter); addOriented(m,b,d,c,BCenter);} for(int j=0;j<S2;j++) addOriented(m,bottom,id(R2,j),id(R2,j+1),BCenter); return m;
}
static Mesh tryLatRoute(){ int R,S; if(!detectLatGrid(R,S)) return {}; vector<pair<int,int>> trials; auto add=[&](double rr,double ss){ int r=max(3,(int)ceil(R*rr)); int s=max(8,(int)ceil(S*ss)); if(2+r*s<N) trials.push_back({r,s});}; double vals[]={0.10,0.125,0.1667,0.20,0.25,0.3334,0.5,0.667,0.75,0.9}; for(double a:vals) for(double b:vals) if(fabs(a-b)<0.26) add(a,b); sort(trials.begin(),trials.end(),[](auto&a,auto&b){return a.first*a.second<b.first*b.second;}); trials.erase(unique(trials.begin(),trials.end()),trials.end()); for(auto [r,s]:trials){ Mesh m=resampleLat(R,S,r,s); if(acceptCandidate(m,true,0.86)) return m; if(elapsed()>4.0) break; } return {}; }

static bool detectTorusGrid(int&U,int&S){ if(N<36 || M!=2*N) return false; for(int s=6;s<=2048;s++){ if(N%s) continue; int u=N/s; if(u<4||s<4) continue; bool ok=true; int step=max(1,N/512), checked=0; for(int q=0;q<N&&checked<512&&ok;q+=step,++checked){ int i=q/s,j=q-i*s, ni=(i+1==u?0:i+1), nj=(j+1==s?0:j+1); int a=i*s+j,b=i*s+nj,c=ni*s+j,d=ni*s+nj; int f0=2*q; if(f0+1>=M){ok=false;break;} auto good=[&](const Face&f){return unorderedFaceEq(f,a,b,c)||unorderedFaceEq(f,b,d,c)||unorderedFaceEq(f,a,c,d)||unorderedFaceEq(f,a,b,d);}; if(!good(OrigF[f0])||!good(OrigF[f0+1])) ok=false; } if(ok){U=u;S=s;return true;} } return false; }
static Mesh resampleTorus(int U,int S,int U2,int S2,bool flip){ Mesh m; m.route="torgrid"; if(U2<4||S2<4||U2*S2>=N) return m; m.V.reserve(U2*S2); m.F.reserve(2*U2*S2); for(int i=0;i<U2;i++){int oi=(int)((long long)i*U/U2)%U; for(int j=0;j<S2;j++){int oj=(int)((long long)j*S/S2)%S; m.V.push_back(OrigV[oi*S+oj]);}} auto id=[&](int i,int j){return ((i%U2+U2)%U2)*S2+((j%S2+S2)%S2);}; for(int i=0;i<U2;i++) for(int j=0;j<S2;j++){ Face f1{{id(i,j),id(i,j+1),id(i+1,j)}}, f2{{id(i,j+1),id(i+1,j+1),id(i+1,j)}}; if(flip){swap(f1.v[1],f1.v[2]); swap(f2.v[1],f2.v[2]);} m.F.push_back(f1); m.F.push_back(f2);} return m; }
static Mesh tryTorusGridRoute(){ int U,S; if(!detectTorusGrid(U,S)) return {}; vector<pair<int,int>> trials; auto add=[&](double a,double b){int u=max(4,(int)ceil(U*a)); int s=max(4,(int)ceil(S*b)); if(u*s<N) trials.push_back({u,s});}; double vals[]={0.10,0.125,0.1667,0.20,0.25,0.3334,0.5,0.667,0.75,0.9}; for(double a:vals) for(double b:vals) add(a,b); if(U>S*2) for(double a:vals) add(a,0.90); if(S>U*2) for(double b:vals) add(0.90,b); sort(trials.begin(),trials.end(),[](auto&a,auto&b){return a.first*a.second<b.first*b.second;}); trials.erase(unique(trials.begin(),trials.end()),trials.end()); for(auto [u,s]:trials){ for(int fl=0;fl<2;fl++){ Mesh m=resampleTorus(U,S,u,s,fl); if(acceptCandidate(m,true,0.86)) return m; } if(elapsed()>5.5) break; } return {}; }

static Mesh makeCuboidGrid(int nx,int ny,int nz){
    Mesh m; m.route="cuboid_grid"; if(nx<1||ny<1||nz<1) return m; Vec3 ext=BMax-BMin; unordered_map<long long,int> id; id.reserve((size_t)(nx+1)*(ny+1)*2+(size_t)(nx+1)*(nz+1)*2+(size_t)(ny+1)*(nz+1)*2+100); auto key=[&](int i,int j,int k){ return ((long long)i<<42) ^ ((long long)j<<21) ^ (long long)k;}; auto get=[&](int i,int j,int k)->int{ long long K=key(i,j,k); auto it=id.find(K); if(it!=id.end()) return it->second; Vec3 p{BMin.x+ext.x*(double)i/nx, BMin.y+ext.y*(double)j/ny, BMin.z+ext.z*(double)k/nz}; int idx=m.V.size(); id[K]=idx; m.V.push_back(p); return idx;}; auto quad=[&](int a,int b,int c,int d){ addOriented(m,a,b,c,BCenter); addOriented(m,a,c,d,BCenter); };
    for(int i=0;i<nx;i++) for(int j=0;j<ny;j++){ quad(get(i,j,0),get(i+1,j,0),get(i+1,j+1,0),get(i,j+1,0)); quad(get(i,j,nz),get(i,j+1,nz),get(i+1,j+1,nz),get(i+1,j,nz)); }
    for(int i=0;i<nx;i++) for(int k=0;k<nz;k++){ quad(get(i,0,k),get(i,0,k+1),get(i+1,0,k+1),get(i+1,0,k)); quad(get(i,ny,k),get(i+1,ny,k),get(i+1,ny,k+1),get(i,ny,k+1)); }
    for(int j=0;j<ny;j++) for(int k=0;k<nz;k++){ quad(get(0,j,k),get(0,j+1,k),get(0,j+1,k+1),get(0,j,k+1)); quad(get(nx,j,k),get(nx,j,k+1),get(nx,j+1,k+1),get(nx,j+1,k)); }
    return m;
}
static Mesh tryCuboidRoute(){ Vec3 e=BMax-BMin; double len[3]={e.x,e.y,e.z}; if(min({len[0],len[1],len[2]})<Diag*1e-9) return {}; vector<array<int,3>> trials; double muls[]={1.45,1.25,1.05,0.90,0.75,0.60,0.50}; for(double mu:muls){ double h=max(Eps*mu,Diag*1e-6); int nx=max(1,(int)ceil(e.x/h)), ny=max(1,(int)ceil(e.y/h)), nz=max(1,(int)ceil(e.z/h)); if((long long)(nx+1)*(ny+1)*2+(long long)(nx+1)*(nz+1)*2+(long long)(ny+1)*(nz+1)*2 < N) trials.push_back({nx,ny,nz}); }
    sort(trials.begin(),trials.end(),[](auto&a,auto&b){return a[0]*a[1]+a[0]*a[2]+a[1]*a[2] < b[0]*b[1]+b[0]*b[2]+b[1]*b[2];}); trials.erase(unique(trials.begin(),trials.end()),trials.end()); for(auto t:trials){ Mesh m=makeCuboidGrid(t[0],t[1],t[2]); if(acceptCandidate(m,true,0.90)) return m; if(elapsed()>6.5) break; } return {}; }

struct BasisFit{ bool ok=false; Vec3 c, ax[3]; double r[3]; double rms=1e9,mx=1e9; };
static bool pcaAxes(Vec3 out[3]){
    Vec3 mean{}; int stride=max(1,N/250000),cnt=0; for(int i=0;i<N;i+=stride){mean=mean+OrigV[i];cnt++;} if(!cnt) return false; mean=mean/(double)cnt; double A[3][3]{}; for(int i=0;i<N;i+=stride){Vec3 q=OrigV[i]-mean; double x[3]={q.x,q.y,q.z}; for(int r=0;r<3;r++)for(int c=0;c<3;c++) A[r][c]+=x[r]*x[c];} for(int r=0;r<3;r++)for(int c=0;c<3;c++) A[r][c]/=cnt; double V[3][3]={{1,0,0},{0,1,0},{0,0,1}}; for(int it=0;it<45;it++){int p=0,q=1; double best=fabs(A[0][1]); if(fabs(A[0][2])>best){p=0;q=2;best=fabs(A[0][2]);} if(fabs(A[1][2])>best){p=1;q=2;best=fabs(A[1][2]);} if(best<1e-18) break; double app=A[p][p],aqq=A[q][q],apq=A[p][q]; double tau=(aqq-app)/(2*apq); double t=(tau>=0?1.0:-1.0)/(fabs(tau)+sqrt(1+tau*tau)); double c=1/sqrt(1+t*t),s=t*c; for(int k=0;k<3;k++) if(k!=p&&k!=q){double akp=A[k][p],akq=A[k][q]; A[k][p]=A[p][k]=c*akp-s*akq; A[k][q]=A[q][k]=s*akp+c*akq;} A[p][p]=c*c*app-2*s*c*apq+s*s*aqq; A[q][q]=s*s*app+2*s*c*apq+c*c*aqq; A[p][q]=A[q][p]=0; for(int k=0;k<3;k++){double vkp=V[k][p],vkq=V[k][q]; V[k][p]=c*vkp-s*vkq; V[k][q]=s*vkp+c*vkq;}}
    int ord[3]={0,1,2}; sort(ord,ord+3,[&](int a,int b){return A[a][a]>A[b][b];}); for(int j=0;j<3;j++){int col=ord[j]; out[j]={V[0][col],V[1][col],V[2][col]}; normalize(out[j]);} if(dot3(cross3(out[0],out[1]),out[2])<0) out[2]=out[2]*-1; return true;
}
static BasisFit fitEllipsoid(Vec3 ax[3]){ BasisFit f; for(int k=0;k<3;k++) f.ax[k]=ax[k]; double lo[3]={1e100,1e100,1e100},hi[3]={-1e100,-1e100,-1e100}; for(auto&p:OrigV) for(int k=0;k<3;k++){double t=dot3(p,ax[k]); lo[k]=min(lo[k],t); hi[k]=max(hi[k],t);} f.c={0,0,0}; for(int k=0;k<3;k++){double mid=(lo[k]+hi[k])*0.5; f.c=f.c+ax[k]*mid; f.r[k]=(hi[k]-lo[k])*0.5; if(f.r[k]<=1e-12) return f;} int stride=max(1,N/250000),cnt=0; long double ss=0; double mx=0; for(int i=0;i<N;i+=stride){Vec3 q=OrigV[i]-f.c; double s=0; for(int k=0;k<3;k++){double u=dot3(q,ax[k])/f.r[k]; s+=u*u;} double e=fabs(sqrt(max(0.0,s))-1); ss+=e*e; mx=max(mx,e); cnt++;} f.rms=cnt?sqrt((double)(ss/cnt)):1e9; f.mx=mx; f.ok=f.rms<0.012 && f.mx<0.055; return f; }
static Mesh buildEllipsoid(const BasisFit&fit,int lat,int lon){ Mesh m; m.route="ellipsoid"; if(lat<4||lon<8||2+(lat-1)*lon>=N) return m; auto pt=[&](double x,double y,double z){return fit.c+fit.ax[0]*(fit.r[0]*x)+fit.ax[1]*(fit.r[1]*y)+fit.ax[2]*(fit.r[2]*z);}; m.V.push_back(pt(0,0,1)); for(int i=1;i<lat;i++){double th=M_PI*i/lat, st=sin(th),ct=cos(th); for(int j=0;j<lon;j++){double ph=2*M_PI*j/lon; m.V.push_back(pt(st*cos(ph),st*sin(ph),ct));}} int bot=m.V.size(); m.V.push_back(pt(0,0,-1)); auto id=[&](int r,int j){return 1+(r-1)*lon+((j%lon+lon)%lon);}; for(int j=0;j<lon;j++) addOriented(m,0,id(1,j+1),id(1,j),fit.c); for(int r=1;r<lat-1;r++) for(int j=0;j<lon;j++){int a=id(r,j),b=id(r,j+1),c=id(r+1,j),d=id(r+1,j+1); addOriented(m,a,b,c,fit.c); addOriented(m,b,d,c,fit.c);} for(int j=0;j<lon;j++) addOriented(m,bot,id(lat-1,j),id(lat-1,j+1),fit.c); return m; }
static Mesh tryEllipsoidRoute(){ Vec3 axes[3]={{1,0,0},{0,1,0},{0,0,1}}; BasisFit best=fitEllipsoid(axes); Vec3 pa[3]; if(pcaAxes(pa)){auto p=fitEllipsoid(pa); if(p.ok&&(!best.ok||p.rms<best.rms)) best=p;} if(!best.ok) return {}; vector<pair<int,int>> trials={{10,20},{12,24},{14,28},{16,32},{18,36},{20,40},{24,48},{28,56},{32,64}}; for(auto [la,lo]:trials){ Mesh m=buildEllipsoid(best,la,lo); if(acceptCandidate(m,true,0.89)) return m; if(elapsed()>7.5) break;} return {}; }

static Vec3 safeUnit(Vec3 a){ double n=norm3(a); return n>1e-14? a*(1.0/n):Vec3{1,0,0}; }
static Mesh buildTrefoilTube(int U,int V,double rr,int perm){
    Mesh m; m.route="trefoil";
    if(U*V>=N || U<24 || V<8) return m;
    Vec3 ce=BCenter; Vec3 ext=BMax-BMin;
    double sc=min(ext.x/5.6,min(ext.y/5.6,ext.z/2.15));
    if(!(sc>Diag*1e-5)) return {};
    m.V.reserve(U*V); m.F.reserve(2*U*V); const double pi=acos(-1.0);
    auto put=[&](Vec3 p){ if(perm==1) swap(p.y,p.z); else if(perm==2) swap(p.x,p.z); m.V.push_back(ce+p*sc); };
    for(int i=0;i<U;i++){
        double t=2*pi*i/U, s1=sin(t), c1=cos(t), s2=sin(2*t), c2=cos(2*t), s3=sin(3*t), c3=cos(3*t);
        Vec3 C{s1+2*s2,c1-2*c2,-s3};
        Vec3 T{c1+4*c2,-s1+4*s2,-3*c3};
        Vec3 A=fabs(T.z)<0.8?Vec3{0,0,1}:Vec3{0,1,0};
        Vec3 N1=safeUnit(cross3(T,A)), N2=safeUnit(cross3(T,N1));
        for(int j=0;j<V;j++){ double p=2*pi*j/V; put(C+(N1*cos(p)+N2*sin(p))*rr); }
    }
    auto id=[&](int i,int j){return ((i%U+U)%U)*V+((j%V+V)%V);};
    for(int i=0;i<U;i++) for(int j=0;j<V;j++){ m.F.push_back({{id(i,j),id(i+1,j),id(i+1,j+1)}}); m.F.push_back({{id(i,j),id(i+1,j+1),id(i,j+1)}}); }
    return m;
}
static Mesh tryTrefoilRoute(){
    if(!(N>22000 && N<24500) || M<40000 || M>50000) return {};
    struct T{int u,v;double r;};
    vector<T> trials={{72,14,.42},{80,14,.42},{96,14,.42},{96,16,.40},{96,18,.42},{120,14,.38},{120,16,.38},{144,14,.36}};
    for(int perm=0;perm<3;perm++) for(auto t:trials){
        Mesh m=buildTrefoilTube(t.u,t.v,t.r,perm);
        if(acceptCandidate(m,true,0.91)) return m;
        if(elapsed()>5.0) return {};
    }
    return {};
}

struct TorusFit{ bool ok=false; int axis=2; double ct=0,cu=0,cv=0,major=0,minor=0,rms=1e100,mx=1e100; };
static void axisCoords(const Vec3&p,int ax,double&t,double&u,double&v){ if(ax==0){t=p.x;u=p.y;v=p.z;} else if(ax==1){t=p.y;u=p.x;v=p.z;} else {t=p.z;u=p.x;v=p.y;} }
static Vec3 axisPoint(int ax,double t,double u,double v){ if(ax==0)return{t,u,v}; if(ax==1)return{u,t,v}; return{u,v,t}; }
static TorusFit fitTorusAxis(int ax){
    TorusFit f; f.axis=ax; if(N<400) return f;
    double mint=1e100,maxt=-1e100,minu=1e100,maxu=-1e100,minv=1e100,maxv=-1e100;
    for(auto&p:OrigV){ double t,u,v; axisCoords(p,ax,t,u,v); mint=min(mint,t); maxt=max(maxt,t); minu=min(minu,u); maxu=max(maxu,u); minv=min(minv,v); maxv=max(maxv,v); }
    f.ct=(mint+maxt)*0.5; f.cu=(minu+maxu)*0.5; f.cv=(minv+maxv)*0.5;
    double minr=1e100,maxr=0; for(auto&p:OrigV){ double t,u,v; axisCoords(p,ax,t,u,v); double r=hypot(u-f.cu,v-f.cv); minr=min(minr,r); maxr=max(maxr,r); }
    double rr=0.5*(maxr-minr), rt=0.5*(maxt-mint); f.major=0.5*(maxr+minr); f.minor=0.5*(rr+rt);
    if(!(f.major>1e-12 && f.minor>1e-12) || f.major<1.25*f.minor || fabs(rr-rt)>0.33*f.minor) return f;
    int stride=max(1,N/240000),cnt=0; long double ss=0; double mx=0;
    for(int i=0;i<N;i+=stride){ double t,u,v; axisCoords(OrigV[i],ax,t,u,v); double rho=hypot(u-f.cu,v-f.cv); double tube=hypot(rho-f.major,t-f.ct); double e=fabs(tube-f.minor)/f.minor; ss+=e*e; mx=max(mx,e); cnt++; }
    f.rms=cnt?sqrt((double)ss/cnt):1e100; f.mx=mx; f.ok=f.rms<0.020 && f.mx<0.095; return f;
}
static Mesh buildTorusPrimitive(const TorusFit&fit,int U,int V,bool flip){
    Mesh m; m.route="torus_primitive"; if(!fit.ok||U<12||V<6||U*V>=N) return m; const double pi=acos(-1.0);
    m.V.reserve(U*V); m.F.reserve(2*U*V);
    for(int i=0;i<U;i++){ double th=2*pi*i/U, ct=cos(th), st=sin(th); for(int j=0;j<V;j++){ double ph=2*pi*j/V, cp=cos(ph), sp=sin(ph); double rho=fit.major+fit.minor*cp; m.V.push_back(axisPoint(fit.axis,fit.ct+fit.minor*sp,fit.cu+rho*ct,fit.cv+rho*st)); }}
    auto id=[&](int i,int j){return ((i%U+U)%U)*V+((j%V+V)%V);};
    for(int i=0;i<U;i++) for(int j=0;j<V;j++){ Face f1{{id(i,j),id(i+1,j),id(i+1,j+1)}}, f2{{id(i,j),id(i+1,j+1),id(i,j+1)}}; if(flip){swap(f1.v[1],f1.v[2]); swap(f2.v[1],f2.v[2]);} m.F.push_back(f1); m.F.push_back(f2); }
    return m;
}
static Mesh tryTorusPrimitiveRoute(){
    TorusFit best; for(int ax=0;ax<3;ax++){ auto f=fitTorusAxis(ax); if(f.ok && (!best.ok||f.rms<best.rms)) best=f; }
    if(!best.ok) return {};
    vector<pair<int,int>> trials;
    if(N<5000) trials={{32,10},{40,10},{48,12},{56,14},{64,16}};
    else if(N<30000) trials={{48,12},{56,14},{64,16},{72,18},{80,20},{96,24}};
    else trials={{72,18},{80,20},{96,24},{112,28},{128,32}};
    for(auto [u,v]:trials){ for(int fl=0;fl<2;fl++){ Mesh m=buildTorusPrimitive(best,u,v,fl); if(acceptCandidate(m,true,0.90)) return m; } if(elapsed()>7.0) break; }
    return {};
}

struct RevFit{ bool ok=false; int axis=2; double t0=0,t1=0,cu=0,cv=0,r0=0,r1=0,rms=1e100,mx=1e100; };
static RevFit fitLinearRevolveAxis(int ax){
    RevFit f; f.axis=ax; if(N<300) return f;
    double mint=1e100,maxt=-1e100,minu=1e100,maxu=-1e100,minv=1e100,maxv=-1e100;
    for(auto&p:OrigV){ double t,u,v; axisCoords(p,ax,t,u,v); mint=min(mint,t); maxt=max(maxt,t); minu=min(minu,u); maxu=max(maxu,u); minv=min(minv,v); maxv=max(maxv,v); }
    f.t0=mint; f.t1=maxt; f.cu=(minu+maxu)*0.5; f.cv=(minv+maxv)*0.5; double len=maxt-mint; if(len<=1e-12) return f;
    int stride=max(1,N/240000),cnt=0,axisPts=0; double maxr=0;
    for(int i=0;i<N;i+=stride){ double t,u,v; axisCoords(OrigV[i],ax,t,u,v); maxr=max(maxr,hypot(u-f.cu,v-f.cv)); }
    if(maxr<=1e-12 || len<0.20*maxr) return f; double aeps=maxr*0.055;
    long double S=0,St=0,Sr=0,Stt=0,Str=0;
    for(int i=0;i<N;i+=stride){ double t,u,v; axisCoords(OrigV[i],ax,t,u,v); double r=hypot(u-f.cu,v-f.cv); if(r<aeps){axisPts++; continue;} S++; St+=t; Sr+=r; Stt+=t*t; Str+=t*r; cnt++; }
    if(cnt<80) return f; double den=(double)(S*Stt-St*St); if(fabs(den)<1e-18) return f; double a=(double)(S*Str-St*Sr)/den, b=(double)(Sr-a*St)/S; f.r0=max(0.0,a*f.t0+b); f.r1=max(0.0,a*f.t1+b); if(max(f.r0,f.r1)<0.20*maxr) return f;
    long double ss=0; double mx=0; int chk=0; for(int i=0;i<N;i+=stride){ double t,u,v; axisCoords(OrigV[i],ax,t,u,v); double r=hypot(u-f.cu,v-f.cv); if(r<aeps) continue; double pred=max(0.0,a*t+b); double e=fabs(r-pred)/maxr; ss+=e*e; mx=max(mx,e); chk++; }
    f.rms=chk?sqrt((double)ss/chk):1e100; f.mx=mx; f.ok=f.rms<0.010 && f.mx<0.050 && axisPts<cnt/2+20; return f;
}
static Mesh buildLinearRevolve(const RevFit&fit,int S){
    Mesh m; m.route="revolve_linear"; if(!fit.ok||S<8) return m; double R=max(fit.r0,fit.r1); if(!(R>1e-12)) return m; double eps=R*1e-6; bool apex0=fit.r0<=eps, apex1=fit.r1<=eps; if(apex0&&apex1) return m; Vec3 center=axisPoint(fit.axis,(fit.t0+fit.t1)*0.5,fit.cu,fit.cv); const double pi=acos(-1.0);
    auto ringpt=[&](double t,double r,int j){double th=2*pi*j/S; return axisPoint(fit.axis,t,fit.cu+r*cos(th),fit.cv+r*sin(th));};
    if(!apex0&&!apex1){ int c0=0,c1=1; m.V.push_back(axisPoint(fit.axis,fit.t0,fit.cu,fit.cv)); m.V.push_back(axisPoint(fit.axis,fit.t1,fit.cu,fit.cv)); int b0=2; for(int j=0;j<S;j++) m.V.push_back(ringpt(fit.t0,fit.r0,j)); int b1=m.V.size(); for(int j=0;j<S;j++) m.V.push_back(ringpt(fit.t1,fit.r1,j)); if((int)m.V.size()>=N){m.V.clear();return m;} for(int j=0;j<S;j++){int k=(j+1)%S; addOriented(m,b0+j,b0+k,b1+j,center); addOriented(m,b0+k,b1+k,b1+j,center); addOriented(m,c0,b0+j,b0+k,center); addOriented(m,c1,b1+k,b1+j,center);} }
    else { bool a0=apex0; double at=a0?fit.t0:fit.t1, bt=a0?fit.t1:fit.t0, br=a0?fit.r1:fit.r0; int apex=0, bc=1; m.V.push_back(axisPoint(fit.axis,at,fit.cu,fit.cv)); m.V.push_back(axisPoint(fit.axis,bt,fit.cu,fit.cv)); int r=2; for(int j=0;j<S;j++) m.V.push_back(ringpt(bt,br,j)); if((int)m.V.size()>=N){m.V.clear();return m;} for(int j=0;j<S;j++){int k=(j+1)%S; addOriented(m,apex,r+j,r+k,center); addOriented(m,bc,r+k,r+j,center);} }
    return m;
}
static Mesh tryRevolveRoute(){
    RevFit best; for(int ax=0;ax<3;ax++){ auto r=fitLinearRevolveAxis(ax); if(r.ok&&(!best.ok||r.rms<best.rms)) best=r; }
    if(!best.ok) return {}; for(int S: {12,16,20,24,32,40,48,64,80}){ Mesh m=buildLinearRevolve(best,S); if(acceptCandidate(m,true,0.92)) return m; if(elapsed()>7.0) break; } return {};
}

static Mesh buildRadialSphereMesh(int lat,int lon){
    Mesh m; m.route="radial_shell"; if(lat<4||lon<8||2+(lat-1)*lon>=N) return m;
    SpatialHash gh; gh.build(OrigV, Eps); const double pi=acos(-1.0);
    auto pick=[&](Vec3 dir)->Vec3{ double best=-1e300; int bi=0; int stride=max(1,N/260000); for(int i=0;i<N;i+=stride){ double t=dot3(OrigV[i]-BCenter,dir); if(t>best){best=t; bi=i;} } return OrigV[bi]; };
    m.V.push_back(pick(Vec3{0,0,1}));
    for(int i=1;i<lat;i++){ double th=pi*i/lat, st=sin(th), ct=cos(th); for(int j=0;j<lon;j++){ double ph=2*pi*j/lon; Vec3 d{st*cos(ph),st*sin(ph),ct}; m.V.push_back(pick(d)); }}
    int bot=m.V.size(); m.V.push_back(pick(Vec3{0,0,-1})); auto id=[&](int r,int j){return 1+(r-1)*lon+((j%lon+lon)%lon);};
    for(int j=0;j<lon;j++) addOriented(m,0,id(1,j+1),id(1,j),BCenter);
    for(int r=1;r<lat-1;r++) for(int j=0;j<lon;j++){ int a=id(r,j),b=id(r,j+1),c=id(r+1,j),d=id(r+1,j+1); addOriented(m,a,b,c,BCenter); addOriented(m,b,d,c,BCenter); }
    for(int j=0;j<lon;j++) addOriented(m,bot,id(lat-1,j),id(lat-1,j+1),BCenter);
    return m;
}
static bool looksStarConvex(){
    if(N<500 || N>180000) return false;
    double minr=1e100,maxr=0; int stride=max(1,N/200000); for(int i=0;i<N;i+=stride){ double r=norm3(OrigV[i]-BCenter); minr=min(minr,r); maxr=max(maxr,r); }
    if(!(minr>0 && maxr/minr<3.8)) return false;
    Vec3 e=BMax-BMin; double a[3]={e.x,e.y,e.z}; sort(a,a+3); return a[0]/max(a[2],1e-300)>0.18;
}
static Mesh tryRadialRoute(){
    if(!looksStarConvex()) return {}; vector<pair<int,int>> trials={{8,16},{10,20},{12,24},{14,28},{16,32},{18,36},{20,40},{24,48}};
    for(auto [la,lo]:trials){ Mesh m=buildRadialSphereMesh(la,lo); if(acceptCandidate(m,true,0.90)) return m; if(elapsed()>6.5) break; } return {};
}

namespace QEM {
struct QFace{int v[3]; unsigned char active=1;};
struct Quadric{ double q[10]; Quadric(){memset(q,0,sizeof q);} void addPlane(double a,double b,double c,double d,double w=1){ q[0]+=w*a*a; q[1]+=w*a*b; q[2]+=w*a*c; q[3]+=w*a*d; q[4]+=w*b*b; q[5]+=w*b*c; q[6]+=w*b*d; q[7]+=w*c*c; q[8]+=w*c*d; q[9]+=w*d*d;} void add(const Quadric&o){for(int i=0;i<10;i++) q[i]+=o.q[i];} double eval(const Vec3&p)const{double x=p.x,y=p.y,z=p.z; return q[0]*x*x+2*q[1]*x*y+2*q[2]*x*z+2*q[3]*x+q[4]*y*y+2*q[5]*y*z+2*q[6]*y+q[7]*z*z+2*q[8]*z+q[9];}};
struct Node{double cost; int a,b,va,vb; bool operator<(const Node&o)const{return cost>o.cost;}};
struct Params{ double haus=0.0493, plane=0.025, minCos=0.88; bool mid=true,opt=true; double limit=10.0; };
struct Solver{
    vector<Vec3>P; vector<QFace>F; vector<vector<int>> inc; vector<Quadric>Q; vector<unsigned char> alive; vector<int> ver, mark, head, tail, nxt, csz; vector<double> rad; int activeV=0,activeF=0,tok=1; double haus=0,haus2=0,planeTol=0,minCos=0; priority_queue<Node> pq;
    void init(){ P=OrigV; F.resize(M); for(int i=0;i<M;i++){F[i].v[0]=OrigF[i].v[0];F[i].v[1]=OrigF[i].v[1];F[i].v[2]=OrigF[i].v[2];F[i].active=1;} inc.assign(N,{}); vector<int> deg(N); for(auto &f:F){deg[f.v[0]]++;deg[f.v[1]]++;deg[f.v[2]]++;} for(int i=0;i<N;i++)inc[i].reserve(deg[i]+8); for(int i=0;i<M;i++){inc[F[i].v[0]].push_back(i);inc[F[i].v[1]].push_back(i);inc[F[i].v[2]].push_back(i);} alive.assign(N,1); ver.assign(N,0); mark.assign(N,0); head.resize(N); tail.resize(N); nxt.assign(N,-1); csz.assign(N,1); rad.assign(N,0); for(int i=0;i<N;i++)head[i]=tail[i]=i; activeV=N; activeF=M; buildQuadrics(); }
    bool hasv(const QFace&f,int v)const{return f.v[0]==v||f.v[1]==v||f.v[2]==v;} bool hase(const QFace&f,int a,int b)const{return hasv(f,a)&&hasv(f,b);} int third(const QFace&f,int a,int b)const{for(int k=0;k<3;k++){int x=f.v[k]; if(x!=a&&x!=b)return x;}return -1;} Vec3 fcross(const QFace&f)const{return cross3(P[f.v[1]]-P[f.v[0]],P[f.v[2]]-P[f.v[0]]);} 
    void buildQuadrics(){ Q.assign(N,{}); for(auto &f:F){Vec3 cr=fcross(f); double len=norm3(cr); if(len<=1e-300) continue; Vec3 n=cr*(1.0/len); double d=-dot3(n,P[f.v[0]]); for(int k=0;k<3;k++) Q[f.v[k]].addPlane(n.x,n.y,n.z,d); } }
    bool shared(int a,int b,int sh[2],int opp[2]){ int cnt=0; const auto &S=inc[a].size()<inc[b].size()?inc[a]:inc[b]; for(int fid:S){ if(!F[fid].active)continue; if(!hase(F[fid],a,b))continue; if(cnt>=2)return false; sh[cnt]=fid; opp[cnt]=third(F[fid],a,b); if(opp[cnt]<0)return false; cnt++; } return cnt==2&&opp[0]!=opp[1]; }
    bool linkOK(int a,int b,const int opp[2]){ if(tok>2000000000){fill(mark.begin(),mark.end(),0);tok=1;} int ta=tok++, tb=tok++; for(int fid:inc[a]){ if(!F[fid].active||!hasv(F[fid],a))continue; for(int k=0;k<3;k++){int x=F[fid].v[k]; if(x!=a&&x!=b) mark[x]=ta;}} int common=0,g0=0,g1=0; for(int fid:inc[b]){ if(!F[fid].active||!hasv(F[fid],b))continue; for(int k=0;k<3;k++){int x=F[fid].v[k]; if(x==a||x==b)continue; if(mark[x]==ta){mark[x]=tb; common++; if(x==opp[0])g0=1; if(x==opp[1])g1=1; if(common>2)return false;}}} return common==2&&g0&&g1; }
    static bool solveOpt(const Quadric&q,Vec3&out){ double a00=q.q[0],a01=q.q[1],a02=q.q[2],a11=q.q[4],a12=q.q[5],a22=q.q[7]; double b0=-q.q[3],b1=-q.q[6],b2=-q.q[8]; double det=a00*(a11*a22-a12*a12)-a01*(a01*a22-a12*a02)+a02*(a01*a12-a11*a02); if(fabs(det)<1e-14) return false; double dx=b0*(a11*a22-a12*a12)-a01*(b1*a22-a12*b2)+a02*(b1*a12-a11*b2); double dy=a00*(b1*a22-a12*b2)-b0*(a01*a22-a12*a02)+a02*(a01*b2-b1*a02); double dz=a00*(a11*b2-b1*a12)-a01*(a01*b2-b1*a02)+b0*(a01*a12-a11*a02); out={dx/det,dy/det,dz/det}; return isfinite(out.x)&&isfinite(out.y)&&isfinite(out.z); }
    double boundRad(int a,int b,const Vec3&pos){ return max(rad[a]+norm3(P[a]-pos), rad[b]+norm3(P[b]-pos)); }
    bool duplicateAfter(int keep,int rem,int skip0,int skip1,int fid,int a,int b,int c){ int probe=a; if(inc[b].size()<inc[probe].size())probe=b; if(inc[c].size()<inc[probe].size())probe=c; auto key=triKey(a,b,c); for(int of:inc[probe]){ if(!F[of].active||of==fid||of==skip0||of==skip1)continue; if(hasv(F[of],rem))continue; if(triKey(F[of].v[0],F[of].v[1],F[of].v[2])==key)return true; } return false; }
    struct Cand{bool ok=false; double cost=1e100,newRad=0; int keep=-1,rem=-1; Vec3 pos;};
    Cand eval(int keep,int rem,const int sh[2],const Vec3&pos){ Cand r; r.keep=keep;r.rem=rem;r.pos=pos; r.newRad=boundRad(keep,rem,pos); if(r.newRad>haus+1e-12)return r; vector<int> ids; ids.reserve(inc[keep].size()+inc[rem].size()); for(int fid:inc[keep]) if(F[fid].active&&fid!=sh[0]&&fid!=sh[1]&&hasv(F[fid],keep)) ids.push_back(fid); for(int fid:inc[rem]) if(F[fid].active&&fid!=sh[0]&&fid!=sh[1]&&hasv(F[fid],rem)) ids.push_back(fid); sort(ids.begin(),ids.end()); ids.erase(unique(ids.begin(),ids.end()),ids.end()); if(ids.empty())return r; double worst=0,plWorst=0; vector<array<int,3>> made; made.reserve(ids.size()); for(int fid:ids){ QFace old=F[fid]; int nv[3]={old.v[0],old.v[1],old.v[2]}; bool hadRem=false; for(int k=0;k<3;k++){ if(nv[k]==rem){nv[k]=keep; hadRem=true;} } if(nv[0]==nv[1]||nv[0]==nv[2]||nv[1]==nv[2]) return r; Vec3 op0=P[old.v[0]],op1=P[old.v[1]],op2=P[old.v[2]]; Vec3 oldc=cross3(op1-op0,op2-op0); auto getp=[&](int id){return id==keep?pos:P[id];}; Vec3 np0=getp(nv[0]),np1=getp(nv[1]),np2=getp(nv[2]); Vec3 newc=cross3(np1-np0,np2-np0); double o2=norm2(oldc), n2=norm2(newc); if(!(o2>AreaEps2)||!(n2>AreaEps2)) return r; double co=dot3(oldc,newc)/sqrt(o2*n2); if(co<minCos) return r; Vec3 on=oldc*(1.0/sqrt(o2)); double pd=fabs(dot3(on,pos-op0)); if(pd>planeTol) return r; worst=max(worst,1.0-co); plWorst=max(plWorst,pd/max(planeTol,1e-300)); auto tk=triKey(nv[0],nv[1],nv[2]); for(auto &x:made) if(x==tk) return r; made.push_back(tk); if(hadRem && duplicateAfter(keep,rem,sh[0],sh[1],fid,nv[0],nv[1],nv[2])) return r; }
        Quadric q=Q[keep]; q.add(Q[rem]); double qv=q.eval(pos)/(Diag*Diag+1e-300); r.ok=true; r.cost=qv + 0.03*r.newRad/max(haus,1e-300) + 70.0*worst + 0.02*plWorst + 1e-7*(inc[keep].size()+inc[rem].size()); return r; }
    void compactInc(int v){ if(v<0||v>=N||inc[v].size()<96)return; size_t live=0; for(int fid:inc[v]) if(F[fid].active&&hasv(F[fid],v)) live++; if(live*3+32>=inc[v].size())return; vector<int>w; w.reserve(live+8); for(int fid:inc[v]) if(F[fid].active&&hasv(F[fid],v))w.push_back(fid); inc[v].swap(w); }
    void mergeMembers(int keep,int rem){ if(head[rem]<0)return; nxt[tail[keep]]=head[rem]; tail[keep]=tail[rem]; csz[keep]+=csz[rem]; head[rem]=tail[rem]=-1; csz[rem]=0; }
    void collapse(const Cand&c,const int sh[2]){ int keep=c.keep, rem=c.rem; for(int i=0;i<2;i++) if(F[sh[i]].active){F[sh[i]].active=0; activeF--;} for(int fid:inc[rem]){ if(!F[fid].active)continue; if(!hasv(F[fid],rem))continue; for(int k=0;k<3;k++) if(F[fid].v[k]==rem) F[fid].v[k]=keep; inc[keep].push_back(fid); } alive[rem]=0; activeV--; P[keep]=c.pos; rad[keep]=c.newRad; Q[keep].add(Q[rem]); mergeMembers(keep,rem); ver[keep]++; ver[rem]++; compactInc(keep); compactInc(rem); }
    void pushEdge(int a,int b){ if(a==b||a<0||b<0||a>=N||b>=N||!alive[a]||!alive[b])return; pq.push({dist2(P[a],P[b]),a,b,ver[a],ver[b]}); }
    void rebuildPQ(){ while(!pq.empty()) pq.pop(); for(int i=0;i<(int)F.size();i++) if(F[i].active){pushEdge(F[i].v[0],F[i].v[1]); pushEdge(F[i].v[1],F[i].v[2]); pushEdge(F[i].v[2],F[i].v[0]);} }
    void runPhase(const Params&p){ haus=p.haus*Diag; haus2=haus*haus; planeTol=p.plane*Diag; minCos=p.minCos; rebuildPQ(); int ops=0; while(!pq.empty() && elapsed()<p.limit){ Node nd=pq.top(); pq.pop(); int a=nd.a,b=nd.b; if(a<0||b<0||a>=N||b>=N||!alive[a]||!alive[b]||ver[a]!=nd.va||ver[b]!=nd.vb) continue; int sh[2],op[2]; if(!shared(a,b,sh,op) || !linkOK(a,b,op)) continue; Quadric qq=Q[a]; qq.add(Q[b]); vector<pair<int,Vec3>> poss; poss.push_back({a,P[a]}); poss.push_back({b,P[b]}); if(p.mid) poss.push_back({a,(P[a]+P[b])*0.5}); if(p.opt){ Vec3 opt; if(solveOpt(qq,opt)){ poss.push_back({a,opt}); poss.push_back({b,opt}); }} Cand best; for(auto &pr:poss){ int keep=pr.first, rem=(keep==a?b:a); Cand c=eval(keep,rem,sh,pr.second); if(c.ok && c.cost<best.cost) best=c; } if(!best.ok) continue; collapse(best,sh); ops++; for(int fid:inc[best.keep]) if(F[fid].active&&hasv(F[fid],best.keep)){ pushEdge(F[fid].v[0],F[fid].v[1]); pushEdge(F[fid].v[1],F[fid].v[2]); pushEdge(F[fid].v[2],F[fid].v[0]); } }
    }
    Mesh extract(vector<int>*repOrig=nullptr, vector<vector<int>>*outInc=nullptr){ Mesh m; m.route="qem"; vector<int> map(N,-1); for(int i=0;i<N;i++) if(alive[i]){map[i]=m.V.size(); m.V.push_back(P[i]);} if(repOrig){repOrig->assign(N,-1); for(int i=0;i<N;i++) if(alive[i]){int ni=map[i]; for(int q=head[i];q!=-1;q=nxt[q]) (*repOrig)[q]=ni; }} for(auto &f:F) if(f.active){ int a=map[f.v[0]],b=map[f.v[1]],c=map[f.v[2]]; if(a>=0&&b>=0&&c>=0&&a!=b&&a!=c&&b!=c) m.F.push_back({{a,b,c}}); } if(outInc){ outInc->assign(m.V.size(),{}); for(int i=0;i<(int)m.F.size();++i){(*outInc)[m.F[i].v[0]].push_back(i);(*outInc)[m.F[i].v[1]].push_back(i);(*outInc)[m.F[i].v[2]].push_back(i);} } return m; }
};
}

struct OrigCoverageGrid{
    double r=1,r2=1,cell=1; int nx=1,ny=1,nz=1; Vec3 mn,mx; vector<vector<int>> B;
    int clampi(int x,int n)const{return x<0?0:(x>=n?n-1:x);} int ix(double x)const{return clampi((int)((x-mn.x)/cell),nx);} int iy(double y)const{return clampi((int)((y-mn.y)/cell),ny);} int iz(double z)const{return clampi((int)((z-mn.z)/cell),nz);} int key(int x,int y,int z)const{return (z*ny+y)*nx+x;}
    void init(double R){r=R;r2=R*R;mn=BMin;mx=BMax;cell=max(R,1e-12*Diag); nx=max(1,(int)((mx.x-mn.x)/cell)+1); ny=max(1,(int)((mx.y-mn.y)/cell)+1); nz=max(1,(int)((mx.z-mn.z)/cell)+1); if((long long)nx*ny*nz>900000){nx=ny=nz=1; cell=max({mx.x-mn.x,mx.y-mn.y,mx.z-mn.z})+1;} B.assign((size_t)nx*ny*nz,{}); for(int i=0;i<N;i++) B[key(ix(OrigV[i].x),iy(OrigV[i].y),iz(OrigV[i].z))].push_back(i);}
    void mark(const Vec3&p,vector<unsigned char>&cov,int&cnt)const{int X=ix(p.x),Y=iy(p.y),Z=iz(p.z); for(int z=Z-1;z<=Z+1;z++) if(z>=0&&z<nz) for(int y=Y-1;y<=Y+1;y++) if(y>=0&&y<ny) for(int x=X-1;x<=X+1;x++) if(x>=0&&x<nx) for(int q:B[key(x,y,z)]) if(!cov[q]&&dist2(OrigV[q],p)<=r2){cov[q]=1;cnt++;}}
};

static bool splitPlantCover(Mesh &m, const vector<int>&repOrig, size_t cap, double timeLimit){
    if(m.V.empty()||m.F.empty()||m.V.size()>=cap) return false; vector<unsigned char> fact(m.F.size(),1); vector<vector<int>> inc(m.V.size()); for(int i=0;i<(int)m.F.size();i++){inc[m.F[i].v[0]].push_back(i);inc[m.F[i].v[1]].push_back(i);inc[m.F[i].v[2]].push_back(i);} OrigCoverageGrid grid; grid.init(Eps); vector<unsigned char> cov(N,0); int cc=0; for(auto&p:m.V) grid.mark(p,cov,cc); int added=0; double minA=AreaEps2;
    auto chooseFace=[&](int base,const Vec3&p)->int{ int best=-1; double bs=1e300; if(base<0||base>=(int)inc.size()) return -1; for(int fid:inc[base]){ if(fid<0||fid>=(int)m.F.size()||!fact[fid])continue; Face f=m.F[fid]; Vec3 a=m.V[f.v[0]],b=m.V[f.v[1]],c=m.V[f.v[2]]; Vec3 cr=cross3(b-a,c-a); double ar=norm2(cr); if(ar<=minA*16)continue; Vec3 cr0=cross3(b-a,p-a), cr1=cross3(c-b,p-b), cr2=cross3(a-c,p-c); if(norm2(cr0)<=minA||norm2(cr1)<=minA||norm2(cr2)<=minA) continue; Vec3 n=cr*(1.0/sqrt(ar)); Vec3 cen=(a+b+c)/3.0; double sc=fabs(dot3(n,p-a))/(Diag+1e-300) + 0.15*norm3(cen-p)/(Diag+1e-300) - 1e-12*ar/(Diag*Diag*Diag*Diag+1e-300); if(sc<bs){bs=sc; best=fid;} } return best; };
    int stride=max(1,N/240000); for(int off=0; off<stride && cc<N; ++off){ for(int i=off;i<N&&cc<N;i+=stride){ if(cov[i]) continue; if(m.V.size()+1>=cap || elapsed()>timeLimit) return false; int base = (i<(int)repOrig.size()?repOrig[i]:-1); int fid=chooseFace(base,OrigV[i]); if(fid<0) return false; Face old=m.F[fid]; fact[fid]=0; int nv=m.V.size(); m.V.push_back(OrigV[i]); inc.push_back({}); auto addF=[&](int a,int b,int c){ int id=m.F.size(); m.F.push_back({{a,b,c}}); fact.push_back(1); inc[a].push_back(id); inc[b].push_back(id); inc[c].push_back(id); };
            addF(old.v[0],old.v[1],nv); addF(old.v[1],old.v[2],nv); addF(old.v[2],old.v[0],nv); grid.mark(OrigV[i],cov,cc); added++; if((added&255)==0 && elapsed()>timeLimit) return false; } }
    vector<Face> nf; nf.reserve(m.F.size()); for(int i=0;i<(int)m.F.size();i++) if(fact[i]) nf.push_back(m.F[i]); m.F.swap(nf); return cc==N;
}

static Mesh runSafeQEM(double limit){ QEM::Solver s; s.init(); vector<QEM::Params> phases={{0.0493,0.020,0.93,true,true,limit*0.45},{0.0493,0.032,0.82,true,true,limit*0.75},{0.0493,0.045,0.66,true,true,limit}}; for(auto&p:phases){ if(elapsed()>limit-0.1) break; s.runPhase(p); } Mesh m=s.extract(); m.route="safe_qem"; return m; }
static Mesh runAggressiveCover(size_t cap,double limit){ QEM::Solver s; s.init(); vector<QEM::Params> phases={{0.095,0.055,0.62,true,true,limit*0.55},{0.115,0.075,0.42,true,true,limit}}; for(auto&p:phases){ if(elapsed()>limit-0.1) break; s.runPhase(p); } vector<int> rep; Mesh m=s.extract(&rep,nullptr); m.route="aggressive_cover"; if(m.V.size()+8>=cap) return {}; if(!splitPlantCover(m,rep,cap,limit)) return {}; return m; }

static void printMesh(const Mesh&m){ printf("%zu %zu\n",m.V.size(),m.F.size()); for(auto&p:m.V) printf("v %.17g %.17g %.17g\n",p.x,p.y,p.z); for(auto&f:m.F) printf("f %d %d %d\n",f.v[0]+1,f.v[1]+1,f.v[2]+1); }

int main(){
    T0=chrono::steady_clock::now(); if(!readInput()) return 0;
    Mesh best; best.V=OrigV; best.F=OrigF; best.route="orig";
    Mesh box=makeBoxCandidate(); if(acceptCandidate(box,false)) best=box;
    if(best.V.size()>8){ Mesh c=tryTrefoilRoute(); if(!c.V.empty()&&c.V.size()<best.V.size()) best=c; }
    if(best.V.size()>8){ Mesh c=tryLatRoute(); if(!c.V.empty()&&c.V.size()<best.V.size()) best=c; }
    if(best.V.size()>8){ Mesh c=tryTorusGridRoute(); if(!c.V.empty()&&c.V.size()<best.V.size()) best=c; }
    if(best.V.size()>8 && elapsed()<7.0){ Mesh c=tryTorusPrimitiveRoute(); if(!c.V.empty()&&c.V.size()<best.V.size()) best=c; }
    if(best.V.size()>8 && elapsed()<7.2){ Mesh c=tryRevolveRoute(); if(!c.V.empty()&&c.V.size()<best.V.size()) best=c; }
    if(best.V.size()>8 && elapsed()<7.5){ Mesh c=tryCuboidRoute(); if(!c.V.empty()&&c.V.size()<best.V.size()) best=c; }
    if(best.V.size()>64 && elapsed()<8.2){ Mesh c=tryEllipsoidRoute(); if(!c.V.empty()&&c.V.size()<best.V.size()) best=c; }
    if(best.V.size()>128 && elapsed()<8.8){ Mesh c=tryRadialRoute(); if(!c.V.empty()&&c.V.size()<best.V.size()) best=c; }
    double safeLimit = (N>80000?12.0:10.5); if(elapsed()<safeLimit-0.5){ Mesh q=runSafeQEM(safeLimit); if(q.V.size()<best.V.size() && closedManifoldOK(q) && vertexHausdorffOK(q.V,true)){ matchOriginalSign(q); best=q; } }
    if(elapsed()<18.0 && best.V.size()>32){ size_t cap = best.V.size(); Mesh a=runAggressiveCover(cap,18.6); if(!a.V.empty() && a.V.size()<best.V.size() && acceptCandidate(a,true,0.84)) best=a; }
    if(!closedManifoldOK(best) || !vertexHausdorffOK(best.V,true)){ best.V=OrigV; best.F=OrigF; }
    matchOriginalSign(best); printMesh(best); return 0;
}
