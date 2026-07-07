#include <bits/stdc++.h>
using namespace std;

struct Vec3{ double x=0,y=0,z=0; };
struct Face{ int a=0,b=0,c=0; };
static inline Vec3 operator+(const Vec3&a,const Vec3&b){ return {a.x+b.x,a.y+b.y,a.z+b.z}; }
static inline Vec3 operator-(const Vec3&a,const Vec3&b){ return {a.x-b.x,a.y-b.y,a.z-b.z}; }
static inline Vec3 operator*(const Vec3&a,double s){ return {a.x*s,a.y*s,a.z*s}; }
static inline Vec3 operator/(const Vec3&a,double s){ return {a.x/s,a.y/s,a.z/s}; }
static inline double dot3(const Vec3&a,const Vec3&b){ return a.x*b.x+a.y*b.y+a.z*b.z; }
static inline Vec3 cross3(const Vec3&a,const Vec3&b){ return {a.y*b.z-a.z*b.y,a.z*b.x-a.x*b.z,a.x*b.y-a.y*b.x}; }
static inline double norm2(const Vec3&a){ return dot3(a,a); }
static inline double norm3(const Vec3&a){ return sqrt(norm2(a)); }
static inline Vec3 unit(Vec3 v){ double l=norm3(v); return l>0 ? v/l : Vec3{0,0,0}; }
static inline long long edgeKey(int a,int b){ if(a>b) swap(a,b); return ((long long)(unsigned int)a<<32)|(unsigned int)b; }
static inline array<int,3> triKey(Face f){ array<int,3> t{f.a,f.b,f.c}; sort(t.begin(),t.end()); return t; }
struct TriHash{ size_t operator()(const array<int,3>&t) const { uint64_t a=t[0],b=t[1],c=t[2]; return (size_t)(a*11995408973635179863ull ^ b*10150724397891781847ull ^ c*7809847782465536323ull); } };

int N=0,M=0;
vector<Vec3> Orig;
vector<Face> OrigF;
Vec3 bbMin,bbMax,bbCtr;
double diagLen=1.0, HTO=1.0;
chrono::steady_clock::time_point t0;
static inline double elapsed(){ return chrono::duration<double>(chrono::steady_clock::now()-t0).count(); }

static void readInput(){
    ios::sync_with_stdio(false); cin.tie(nullptr);
    cin>>N>>M;
    Orig.resize(N); OrigF.resize(M);
    bbMin={1e100,1e100,1e100}; bbMax={-1e100,-1e100,-1e100};
    char ch;
    for(int i=0;i<N;i++){
        cin>>ch>>Orig[i].x>>Orig[i].y>>Orig[i].z;
        bbMin.x=min(bbMin.x,Orig[i].x); bbMin.y=min(bbMin.y,Orig[i].y); bbMin.z=min(bbMin.z,Orig[i].z);
        bbMax.x=max(bbMax.x,Orig[i].x); bbMax.y=max(bbMax.y,Orig[i].y); bbMax.z=max(bbMax.z,Orig[i].z);
    }
    for(int i=0;i<M;i++){ cin>>ch>>OrigF[i].a>>OrigF[i].b>>OrigF[i].c; --OrigF[i].a; --OrigF[i].b; --OrigF[i].c; }
    bbCtr=(bbMin+bbMax)*0.5; diagLen=norm3(bbMax-bbMin); if(!(diagLen>0)) diagLen=1.0; HTO=0.0490*diagLen;
}

static void printMesh(const vector<Vec3>&V,const vector<Face>&F){
    cout.setf(ios::fixed); cout<<setprecision(10);
    cout<<V.size()<<" "<<F.size()<<"\n";
    for(const Vec3&p:V) cout<<"v "<<p.x<<" "<<p.y<<" "<<p.z<<"\n";
    for(const Face&f:F) cout<<"f "<<f.a+1<<" "<<f.b+1<<" "<<f.c+1<<"\n";
}
static void printOriginal(){ printMesh(Orig,OrigF); }

struct NearGrid{
    Vec3 mn; double cell=1.0; unordered_map<long long, vector<int>> b;
    static long long pack(int x,int y,int z){ const long long B=1048576LL,O=524288LL; return ((long long)(x+O)*B+(y+O))*B+(z+O); }
    NearGrid(){}
    NearGrid(const vector<Vec3>&P,double c){ init(P,c); }
    void init(const vector<Vec3>&P,double c){
        cell=max(c,1e-12); mn={1e100,1e100,1e100};
        for(const Vec3&p:P){ mn.x=min(mn.x,p.x); mn.y=min(mn.y,p.y); mn.z=min(mn.z,p.z); }
        b.reserve(P.size()*2+9);
        for(int i=0;i<(int)P.size();++i){ int x,y,z; idx(P[i],x,y,z); b[pack(x,y,z)].push_back(i); }
    }
    void idx(const Vec3&p,int&x,int&y,int&z) const { x=(int)floor((p.x-mn.x)/cell); y=(int)floor((p.y-mn.y)/cell); z=(int)floor((p.z-mn.z)/cell); }
    bool near(const vector<Vec3>&P,const Vec3&q,double r) const{
        int X,Y,Z; idx(q,X,Y,Z); int R=(int)ceil(r/cell)+1; double r2=r*r;
        for(int dx=-R;dx<=R;dx++) for(int dy=-R;dy<=R;dy++) for(int dz=-R;dz<=R;dz++){
            auto it=b.find(pack(X+dx,Y+dy,Z+dz)); if(it==b.end()) continue;
            for(int id:it->second) if(norm2(P[id]-q)<=r2) return true;
        }
        return false;
    }
};

static bool symmetricOK(const vector<Vec3>&V,double tol){
    if(V.empty() || V.size()>(size_t)N) return false;
    NearGrid gv(V,tol);
    for(const Vec3&p:Orig) if(!gv.near(V,p,tol)) return false;
    NearGrid go(Orig,tol);
    for(const Vec3&p:V) if(!go.near(Orig,p,tol)) return false;
    return true;
}

static bool closedOK(const vector<Vec3>&V,const vector<Face>&F){
    if(V.empty()||F.empty()||V.size()>(size_t)N) return false;
    double eps=max(1e-32,1e-27*diagLen*diagLen);
    unordered_map<long long,int> ec; ec.reserve(F.size()*4+9);
    unordered_set<array<int,3>,TriHash> ts; ts.reserve(F.size()*2+9);
    vector<unsigned char> used(V.size(),0);
    for(const Face&f:F){
        if(f.a<0||f.b<0||f.c<0||f.a>=(int)V.size()||f.b>=(int)V.size()||f.c>=(int)V.size()) return false;
        if(f.a==f.b||f.a==f.c||f.b==f.c) return false;
        if(norm2(cross3(V[f.b]-V[f.a],V[f.c]-V[f.a]))<=eps) return false;
        if(!ts.insert(triKey(f)).second) return false;
        ec[edgeKey(f.a,f.b)]++; ec[edgeKey(f.b,f.c)]++; ec[edgeKey(f.c,f.a)]++;
        used[f.a]=used[f.b]=used[f.c]=1;
    }
    for(auto &kv:ec) if(kv.second!=2) return false;
    for(unsigned char u:used) if(!u) return false;
    return true;
}

static void orientOut(vector<Vec3>&V,Face&f,const Vec3&ref){
    Vec3 n=cross3(V[f.b]-V[f.a],V[f.c]-V[f.a]); Vec3 c=(V[f.a]+V[f.b]+V[f.c])/3.0;
    if(dot3(n,c-ref)<0) swap(f.b,f.c);
}
static void addFace(vector<Vec3>&V,vector<Face>&F,int a,int b,int c,const Vec3&ref){ if(a==b||a==c||b==c) return; Face f{a,b,c}; orientOut(V,f,ref); F.push_back(f); }
static bool finalOK(const vector<Vec3>&V,const vector<Face>&F){ return V.size()<Orig.size() && closedOK(V,F) && symmetricOK(V,HTO*0.998); }

static bool tryDenseBox(vector<Vec3>&V,vector<Face>&F){
    if(N<1200) return false;
    double dx=bbMax.x-bbMin.x,dy=bbMax.y-bbMin.y,dz=bbMax.z-bbMin.z;
    if(min({dx,dy,dz})<0.025*diagLen) return false;
    int st=max(1,N/260000),cnt=0,bad=0,h[6]={}; double sum=0,mx=0;
    for(int i=0;i<N;i+=st){
        const Vec3&p=Orig[i]; double d[6]={fabs(p.x-bbMin.x),fabs(p.x-bbMax.x),fabs(p.y-bbMin.y),fabs(p.y-bbMax.y),fabs(p.z-bbMin.z),fabs(p.z-bbMax.z)};
        int b=0; for(int k=1;k<6;k++) if(d[k]<d[b]) b=k;
        h[b]++; cnt++; sum+=d[b]; mx=max(mx,d[b]); if(d[b]>0.0065*diagLen) bad++;
    }
    if(cnt<60 || bad*1000>cnt*2 || sum>cnt*0.0018*diagLen || mx>0.012*diagLen) return false;
    for(int i=0;i<6;i++) if(h[i]<max(2,cnt/7000)) return false;
    double step=max(1e-12,0.66*HTO);
    int nx=max(1,(int)ceil(dx/step)),ny=max(1,(int)ceil(dy/step)),nz=max(1,(int)ceil(dz/step));
    if(2LL*(nx+1)*(ny+1)+2LL*(nx+1)*(nz+1)+2LL*(ny+1)*(nz+1)>=N) return false;
    V.clear(); F.clear(); unordered_map<long long,int> id; id.reserve(N);
    auto key=[&](int i,int j,int k){ return ((long long)i*(ny+1)+j)*(nz+1)+k; };
    auto get=[&](int i,int j,int k)->int{
        long long kk=key(i,j,k); auto it=id.find(kk); if(it!=id.end()) return it->second;
        Vec3 p{bbMin.x+dx*(double)i/nx,bbMin.y+dy*(double)j/ny,bbMin.z+dz*(double)k/nz}; int r=(int)V.size(); V.push_back(p); id[kk]=r; return r;
    };
    for(int i=0;i<nx;i++) for(int j=0;j<ny;j++){
        int a=get(i,j,0),b=get(i+1,j,0),c=get(i+1,j+1,0),d=get(i,j+1,0); addFace(V,F,a,b,c,bbCtr); addFace(V,F,a,c,d,bbCtr);
        a=get(i,j,nz);b=get(i+1,j,nz);c=get(i+1,j+1,nz);d=get(i,j+1,nz); addFace(V,F,a,c,b,bbCtr); addFace(V,F,a,d,c,bbCtr);
    }
    for(int i=0;i<nx;i++) for(int k=0;k<nz;k++){
        int a=get(i,0,k),b=get(i+1,0,k),c=get(i+1,0,k+1),d=get(i,0,k+1); addFace(V,F,a,b,c,bbCtr); addFace(V,F,a,c,d,bbCtr);
        a=get(i,ny,k);b=get(i+1,ny,k);c=get(i+1,ny,k+1);d=get(i,ny,k+1); addFace(V,F,a,c,b,bbCtr); addFace(V,F,a,d,c,bbCtr);
    }
    for(int j=0;j<ny;j++) for(int k=0;k<nz;k++){
        int a=get(0,j,k),b=get(0,j+1,k),c=get(0,j+1,k+1),d=get(0,j,k+1); addFace(V,F,a,b,c,bbCtr); addFace(V,F,a,c,d,bbCtr);
        a=get(nx,j,k);b=get(nx,j+1,k);c=get(nx,j+1,k+1);d=get(nx,j,k+1); addFace(V,F,a,c,b,bbCtr); addFace(V,F,a,d,c,bbCtr);
    }
    return finalOK(V,F);
}

static int cubeFace(Vec3 q,double&u,double&v){
    double ax=fabs(q.x),ay=fabs(q.y),az=fabs(q.z);
    if(ax>=ay && ax>=az){ u=q.y/ax; v=q.z/ax; return q.x>=0?0:1; }
    if(ay>=az){ u=q.x/ay; v=q.z/ay; return q.y>=0?2:3; }
    u=q.x/az; v=q.y/az; return q.z>=0?4:5;
}
static Vec3 dirFromLattice(int ix,int iy,int iz,int R){
    double x=2.0*ix/R-1.0, y=2.0*iy/R-1.0, z=2.0*iz/R-1.0;
    return unit({x,y,z});
}
static bool tryCubemapShell(int R,vector<Vec3>&V,vector<Face>&F){
    if(N<5000 || elapsed()>14.0) return false;
    vector<vector<int>> bucket(6*R*R);
    for(int i=0;i<N;i++){
        Vec3 q=Orig[i]-bbCtr; if(norm2(q)==0) continue; double u,v; int f=cubeFace(q,u,v);
        int x=min(R-1,max(0,(int)((u+1.0)*0.5*R))), y=min(R-1,max(0,(int)((v+1.0)*0.5*R)));
        bucket[(f*R+y)*R+x].push_back(i);
    }
    unordered_map<long long,int> vid; vid.reserve(6*R*R+9); V.clear(); F.clear();
    auto latticeKey=[&](int x,int y,int z){ return ((long long)x*(R+1)+y)*(R+1)+z; };
    auto selectVertex=[&](int ix,int iy,int iz)->Vec3{
        Vec3 d=dirFromLattice(ix,iy,iz,R); double u,v; int f=cubeFace(d,u,v);
        int cx=min(R-1,max(0,(int)((u+1.0)*0.5*R))), cy=min(R-1,max(0,(int)((v+1.0)*0.5*R)));
        int best=-1; double score=-1e300;
        for(int rad=0;rad<=5 && best<0;rad++){
            for(int yy=max(0,cy-rad);yy<=min(R-1,cy+rad);yy++) for(int xx=max(0,cx-rad);xx<=min(R-1,cx+rad);xx++){
                for(int id: bucket[(f*R+yy)*R+xx]){
                    double sc=dot3(Orig[id]-bbCtr,d) - 0.015*norm2(unit(Orig[id]-bbCtr)-d)*diagLen;
                    if(sc>score){ score=sc; best=id; }
                }
            }
        }
        if(best<0) return bbCtr + d*(0.5*diagLen);
        double rr=max(0.0,dot3(Orig[best]-bbCtr,d));
        return bbCtr + d*rr;
    };
    auto get=[&](int ix,int iy,int iz)->int{
        long long k=latticeKey(ix,iy,iz); auto it=vid.find(k); if(it!=vid.end()) return it->second;
        Vec3 p=selectVertex(ix,iy,iz); int id=(int)V.size(); V.push_back(p); vid[k]=id; return id;
    };
    auto quad=[&](int a,int b,int c,int d){ addFace(V,F,a,b,c,bbCtr); addFace(V,F,a,c,d,bbCtr); };
    for(int i=0;i<R;i++) for(int j=0;j<R;j++){
        quad(get(i,j,0),get(i+1,j,0),get(i+1,j+1,0),get(i,j+1,0));
        quad(get(i,j,R),get(i,j+1,R),get(i+1,j+1,R),get(i+1,j,R));
        quad(get(i,0,j),get(i+1,0,j),get(i+1,0,j+1),get(i,0,j+1));
        quad(get(i,R,j),get(i,R,j+1),get(i+1,R,j+1),get(i+1,R,j));
        quad(get(0,i,j),get(0,i,j+1),get(0,i+1,j+1),get(0,i+1,j));
        quad(get(R,i,j),get(R,i+1,j),get(R,i+1,j+1),get(R,i,j+1));
    }
    return finalOK(V,F);
}

struct Decimator{
    vector<Vec3> V; vector<Face> F; vector<unsigned char> av,af,prot; vector<vector<int>> inc; vector<double> err; int nv,nf; double eps; vector<int> ma,mb; int sta=1,stb=1;
    struct Cand{ double c; int a,b; bool operator<(const Cand&o) const { return c>o.c; } }; priority_queue<Cand> pq;
    Decimator(){
        V=Orig; F=OrigF; nv=N; nf=M; av.assign(N,1); af.assign(M,1); prot.assign(N,0); err.assign(N,0); inc.assign(N,{}); eps=max(1e-32,1e-27*diagLen*diagLen);
        for(int i=0;i<M;i++){ inc[F[i].a].push_back(i); inc[F[i].b].push_back(i); inc[F[i].c].push_back(i); }
        ma.assign(N,0); mb.assign(N,0); protect(); seed();
    }
    bool has(const Face&f,int v) const { return f.a==v||f.b==v||f.c==v; }
    Vec3 fn(int id) const { return cross3(V[F[id].b]-V[F[id].a],V[F[id].c]-V[F[id].a]); }
    void clean(int v){ vector<int> r; for(int id:inc[v]) if(id>=0&&id<(int)F.size()&&af[id]&&has(F[id],v)) r.push_back(id); sort(r.begin(),r.end()); r.erase(unique(r.begin(),r.end()),r.end()); inc[v].swap(r); }
    void protect(){
        vector<Vec3> n(M); for(int i=0;i<M;i++) n[i]=unit(cross3(Orig[OrigF[i].b]-Orig[OrigF[i].a],Orig[OrigF[i].c]-Orig[OrigF[i].a]));
        unordered_map<long long,int> first; first.reserve(3*M+9); double sc=cos(36.0*M_PI/180.0);
        auto edge=[&](int a,int b,int id){ long long k=edgeKey(a,b); auto it=first.find(k); if(it==first.end()) first[k]=id; else if(dot3(n[id],n[it->second])<sc) prot[a]=prot[b]=1; };
        for(int i=0;i<M;i++){ edge(OrigF[i].a,OrigF[i].b,i); edge(OrigF[i].b,OrigF[i].c,i); edge(OrigF[i].c,OrigF[i].a,i); }
        int R=N>150000?72:(N>40000?56:42); double sx=max(1e-30,bbMax.x-bbMin.x),sy=max(1e-30,bbMax.y-bbMin.y),sz=max(1e-30,bbMax.z-bbMin.z);
        auto mark=[&](int ax){ vector<int> lo(R*R,-1),hi(R*R,-1); vector<double> lv(R*R,1e300),hv(R*R,-1e300); for(int i=0;i<N;i++){ double u,v,d; if(ax==0){u=(Orig[i].y-bbMin.y)/sy;v=(Orig[i].z-bbMin.z)/sz;d=Orig[i].x;} else if(ax==1){u=(Orig[i].x-bbMin.x)/sx;v=(Orig[i].z-bbMin.z)/sz;d=Orig[i].y;} else {u=(Orig[i].x-bbMin.x)/sx;v=(Orig[i].y-bbMin.y)/sy;d=Orig[i].z;} int x=min(R-1,max(0,(int)(u*R))),y=min(R-1,max(0,(int)(v*R))),id=x*R+y; if(d<lv[id])lv[id]=d,lo[id]=i; if(d>hv[id])hv[id]=d,hi[id]=i; } for(int id:lo)if(id>=0)prot[id]=1; for(int id:hi)if(id>=0)prot[id]=1; };
        mark(0); mark(1); mark(2);
    }
    void push(int a,int b){ if(a==b||!av[a]||!av[b]) return; double pc=(prot[a]||prot[b])?1.6:0.0; pq.push({norm3(V[a]-V[b])+pc*HTO+0.2*(err[a]+err[b]),a,b}); }
    void seed(){ unordered_set<long long> seen; seen.reserve(3*M+9); for(const Face&f:F){ int a[3]={f.a,f.b,f.c},b[3]={f.b,f.c,f.a}; for(int k=0;k<3;k++) if(seen.insert(edgeKey(a[k],b[k])).second) push(a[k],b[k]); } }
    int opp(const Face&f,int a,int b) const { if(f.a!=a&&f.a!=b) return f.a; if(f.b!=a&&f.b!=b) return f.b; return f.c; }
    bool edgeFaces(int a,int b,int ef[2],int op[2]){ if(inc[a].size()>1500) clean(a); if(inc[b].size()>1500) clean(b); int c=0; const vector<int>&L=inc[a].size()<inc[b].size()?inc[a]:inc[b]; for(int id:L) if(af[id]&&has(F[id],a)&&has(F[id],b)){ if(c>=2) return false; ef[c]=id; op[c]=opp(F[id],a,b); c++; } return c==2&&op[0]!=op[1]; }
    bool linkOK(int a,int b,const int op[2]){ if(++sta==INT_MAX){fill(ma.begin(),ma.end(),0);sta=1;} if(++stb==INT_MAX){fill(mb.begin(),mb.end(),0);stb=1;} for(int id:inc[a]) if(af[id]&&has(F[id],a)){ int xs[3]={F[id].a,F[id].b,F[id].c}; for(int x:xs) if(x!=a&&x!=b) ma[x]=sta; } int cc=0; for(int id:inc[b]) if(af[id]&&has(F[id],b)){ int xs[3]={F[id].a,F[id].b,F[id].c}; for(int x:xs) if(x!=a&&x!=b&&ma[x]==sta){ if(x!=op[0]&&x!=op[1]) return false; if(mb[x]!=stb){mb[x]=stb;cc++;} } } return cc==2&&mb[op[0]]==stb&&mb[op[1]]==stb; }
    bool dup(int fid,int rem,const Face&nf,const int ef[2]){ auto k=triKey(nf); int vs[3]={nf.a,nf.b,nf.c}; int best=vs[0]; if(inc[vs[1]].size()<inc[best].size()) best=vs[1]; if(inc[vs[2]].size()<inc[best].size()) best=vs[2]; for(int id:inc[best]) if(af[id]&&id!=fid&&id!=ef[0]&&id!=ef[1]){ if(has(F[id],rem)) continue; if(triKey(F[id])==k) return true; } return false; }
    bool tryDir(int keep,int rem,const int ef[2],double&ne,double&score){
        if(prot[rem]) return false; ne=max(err[keep],err[rem]+norm3(V[keep]-V[rem])); if(ne>0.982*HTO) return false; double worst=1,pl=0; int cnt=0;
        for(int id:inc[rem]){ if(!af[id]||!has(F[id],rem)||id==ef[0]||id==ef[1]) continue; Face old=F[id],nf=old; if(nf.a==rem)nf.a=keep; if(nf.b==rem)nf.b=keep; if(nf.c==rem)nf.c=keep; if(nf.a==nf.b||nf.a==nf.c||nf.b==nf.c) return false; Vec3 no=fn(id),nn=cross3(V[nf.b]-V[nf.a],V[nf.c]-V[nf.a]); double ao=norm3(no),an=norm3(nn); if(ao*ao<=eps||an*an<=eps) return false; double cs=dot3(no,nn)/(ao*an); worst=min(worst,cs); bool hard=prot[old.a]||prot[old.b]||prot[old.c]; if(cs<(hard?0.990:0.940)) return false; double pd=fabs(dot3(no/ao,V[keep]-V[old.a])); pl=max(pl,pd); if(pd>(hard?0.20:0.48)*HTO) return false; if(dup(id,rem,nf,ef)) return false; cnt++; }
        if(!cnt) return false; score=ne/HTO+0.04*cnt+0.8*(1-worst)+0.35*pl/HTO; return true;
    }
    void collapse(int keep,int rem,const int ef[2],double ne){ for(int i=0;i<2;i++) if(af[ef[i]]){af[ef[i]]=0;nf--;} for(int id:inc[rem]) if(af[id]&&has(F[id],rem)){ if(F[id].a==rem)F[id].a=keep; if(F[id].b==rem)F[id].b=keep; if(F[id].c==rem)F[id].c=keep; inc[keep].push_back(id); } av[rem]=0;nv--;err[keep]=ne;inc[rem].clear();clean(keep); for(int id:inc[keep]) if(af[id]){Face f=F[id];push(f.a,f.b);push(f.b,f.c);push(f.c,f.a);} }
    bool step(){ if(pq.empty()) return false; Cand q=pq.top(); pq.pop(); int a=q.a,b=q.b; if(a==b||!av[a]||!av[b]) return false; int ef[2],op[2]; if(!edgeFaces(a,b,ef,op)||!linkOK(a,b,op)) return false; double ea,eb,sa,sb; bool oa=tryDir(a,b,ef,ea,sa),ob=tryDir(b,a,ef,eb,sb); if(!oa&&!ob) return false; if(ob&&(!oa||sb<sa)) collapse(b,a,ef,eb); else collapse(a,b,ef,ea); return true; }
    void run(){ int bad=0,target=max(64,N/3); double lim=N>100000?18.4:17.0; while(!pq.empty()&&elapsed()<lim&&nv>target){ bool ok=step(); if(ok) bad=0; else if(++bad>200000&&bad>(int)pq.size()/3) break; } }
    void exportMesh(vector<Vec3>&OV,vector<Face>&OF){ vector<int> mp(N,-1); OV.clear(); OF.clear(); for(int i=0;i<N;i++) if(av[i]){mp[i]=(int)OV.size(); OV.push_back(V[i]);} for(int i=0;i<(int)F.size();i++) if(af[i]){ Face f=F[i]; if(mp[f.a]<0||mp[f.b]<0||mp[f.c]<0) continue; Face g{mp[f.a],mp[f.b],mp[f.c]}; if(g.a!=g.b&&g.a!=g.c&&g.b!=g.c) OF.push_back(g); } }
};
static bool tryDecimator(vector<Vec3>&V,vector<Face>&F){ if(N<1500) return false; Decimator D; D.run(); D.exportMesh(V,F); return finalOK(V,F); }

int main(){
    t0=chrono::steady_clock::now(); readInput();
    // Hard sample gate: the official sample cube has 8 vertices, 12 faces. Keep all small meshes unchanged.
    if(N<=1000){ printOriginal(); return 0; }
    vector<Vec3> bestV,V; vector<Face> bestF,F; bool have=false;
    auto consider=[&](bool ok){ if(ok && (!have || V.size()<bestV.size())){ bestV=V; bestF=F; have=true; } };
    consider(tryDenseBox(V,F));
    int R0 = N>250000?72:(N>100000?60:(N>40000?48:36));
    consider(tryCubemapShell(R0,V,F));
    if(elapsed()<16.5) consider(tryCubemapShell(max(20,R0*5/4),V,F));
    if(have){ printMesh(bestV,bestF); return 0; }
    if(tryDecimator(V,F)){ printMesh(V,F); return 0; }
    printOriginal(); return 0;
}