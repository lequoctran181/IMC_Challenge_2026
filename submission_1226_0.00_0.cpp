#include <bits/stdc++.h>
using namespace std;

struct P{double x,y,z;};
static inline P operator+(const P&a,const P&b){return {a.x+b.x,a.y+b.y,a.z+b.z};}
static inline P operator-(const P&a,const P&b){return {a.x-b.x,a.y-b.y,a.z-b.z};}
static inline P operator*(const P&a,double s){return {a.x*s,a.y*s,a.z*s};}
static inline double dotp(const P&a,const P&b){return a.x*b.x+a.y*b.y+a.z*b.z;}
static inline P crossp(const P&a,const P&b){return {a.y*b.z-a.z*b.y,a.z*b.x-a.x*b.z,a.x*b.y-a.y*b.x};}
static inline double norm2(const P&a){return dotp(a,a);} 
static inline double normv(const P&a){return sqrt(norm2(a));}
static inline double coord(const P&p,int ax){return ax==0?p.x:(ax==1?p.y:p.z);} 
static inline double dist2(const P&a,const P&b){return norm2(a-b);} 

struct Tri{int a,b,c;};
struct Mesh{vector<P> v; vector<Tri> f;};

int N,M;
vector<P> V;
vector<Tri> F;
P mnv,mxv,cen;
double DIAG,EPS,EPS2,AREA;

struct H64{static uint64_t smix(uint64_t x){x+=0x9e3779b97f4a7c15ULL;x=(x^(x>>30))*0xbf58476d1ce4e5b9ULL;x=(x^(x>>27))*0x94d049bb133111ebULL;return x^(x>>31);} size_t operator()(uint64_t x)const{return smix(x);} };

struct KD{
    const vector<P>* p; vector<int> id;
    KD(){}
    KD(const vector<P>&pts){init(pts);} 
    void init(const vector<P>&pts){p=&pts; id.resize(pts.size()); iota(id.begin(),id.end(),0); build(0,(int)id.size(),0);} 
    void build(int l,int r,int d){ if(r-l<=1) return; int m=(l+r)>>1, ax=d%3; nth_element(id.begin()+l,id.begin()+m,id.begin()+r,[&](int i,int j){return coord((*p)[i],ax)<coord((*p)[j],ax);}); build(l,m,d+1); build(m+1,r,d+1);} 
    void near_rec(const P&q,int l,int r,int d,double&best)const{ if(l>=r) return; int m=(l+r)>>1, ax=d%3; const P&pt=(*p)[id[m]]; double dd=dist2(q,pt); if(dd<best) best=dd; double df=coord(q,ax)-coord(pt,ax); int l1=l,r1=m,l2=m+1,r2=r; if(df>0){swap(l1,l2);swap(r1,r2);} near_rec(q,l1,r1,d+1,best); if(df*df<best) near_rec(q,l2,r2,d+1,best);} 
    double nearest2(const P&q)const{double b=1e300; near_rec(q,0,(int)id.size(),0,b); return b;}
};

static bool validate_mesh(const Mesh&me){
    int n=me.v.size(), m=me.f.size(); if(n<4||m<4||n>N) return false;
    vector<int> inc(n,0);
    struct EI{int c=0,s=0;};
    unordered_map<uint64_t,EI,H64> ed; ed.reserve((size_t)m*4);
    double tiny=max(1e-30,DIAG*DIAG*1e-24);
    for(auto &t:me.f){
        int a=t.a,b=t.b,c=t.c; if(a<0||b<0||c<0||a>=n||b>=n||c>=n||a==b||b==c||c==a) return false;
        P cr=crossp(me.v[b]-me.v[a],me.v[c]-me.v[a]); if(norm2(cr)<=tiny) return false;
        inc[a]++; inc[b]++; inc[c]++;
        int x[3]={a,b,c}, y[3]={b,c,a};
        for(int i=0;i<3;i++){int u=x[i],v=y[i]; int p=min(u,v),q=max(u,v); uint64_t key=((uint64_t)p<<32)|(uint32_t)q; auto &e=ed[key]; e.c++; e.s += (u==p?1:-1); if(e.c>2) return false;}
    }
    for(int x:inc) if(!x) return false;
    for(auto &kv:ed) if(kv.second.c!=2||kv.second.s!=0) return false;
    return true;
}

static void print_mesh(const Mesh&me){
    cout.setf(ios::fmtflags(0), ios::floatfield); cout<<setprecision(10);
    cout<<me.v.size()<<' '<<me.f.size()<<'\n';
    for(auto &p:me.v) cout<<p.x<<' '<<p.y<<' '<<p.z<<'\n';
    for(auto &t:me.f) cout<<t.a+1<<' '<<t.b+1<<' '<<t.c+1<<'\n';
}
static void print_original(){ Mesh me; me.v=V; me.f=F; print_mesh(me); }

static Mesh gen_box(int nx,int ny,int nz){
    nx=max(1,nx);ny=max(1,ny);nz=max(1,nz);
    Mesh me; vector<int> id((nx+1)*(ny+1)*(nz+1),-1);
    auto lin=[&](int i,int j,int k){return (i*(ny+1)+j)*(nz+1)+k;};
    auto get=[&](int i,int j,int k)->int{
        int &r=id[lin(i,j,k)]; if(r>=0) return r;
        double x=mnv.x+(mxv.x-mnv.x)*i/nx, y=mnv.y+(mxv.y-mnv.y)*j/ny, z=mnv.z+(mxv.z-mnv.z)*k/nz;
        r=(int)me.v.size(); me.v.push_back({x,y,z}); return r;
    };
    auto add=[&](int a,int b,int c,int d){me.f.push_back({a,b,c}); me.f.push_back({a,c,d});};
    for(int i=0;i<nx;i++)for(int j=0;j<ny;j++){int a=get(i,j,nz),b=get(i+1,j,nz),c=get(i+1,j+1,nz),d=get(i,j+1,nz);add(a,b,c,d);} // +z
    for(int i=0;i<nx;i++)for(int j=0;j<ny;j++){int a=get(i,j,0),b=get(i+1,j,0),c=get(i+1,j+1,0),d=get(i,j+1,0);add(a,d,c,b);} // -z
    for(int j=0;j<ny;j++)for(int k=0;k<nz;k++){int a=get(nx,j,k),b=get(nx,j+1,k),c=get(nx,j+1,k+1),d=get(nx,j,k+1);add(a,b,c,d);} // +x
    for(int j=0;j<ny;j++)for(int k=0;k<nz;k++){int a=get(0,j,k),b=get(0,j+1,k),c=get(0,j+1,k+1),d=get(0,j,k+1);add(a,d,c,b);} // -x
    for(int i=0;i<nx;i++)for(int k=0;k<nz;k++){int a=get(i,ny,k),b=get(i,ny,k+1),c=get(i+1,ny,k+1),d=get(i+1,ny,k);add(a,b,c,d);} // +y
    for(int i=0;i<nx;i++)for(int k=0;k<nz;k++){int a=get(i,0,k),b=get(i,0,k+1),c=get(i+1,0,k+1),d=get(i+1,0,k);add(a,d,c,b);} // -y
    return me;
}

static Mesh gen_ellipsoid(int L,int U){
    L=max(6,L); U=max(12,U); Mesh me;
    double rx=(mxv.x-mnv.x)*0.5, ry=(mxv.y-mnv.y)*0.5, rz=(mxv.z-mnv.z)*0.5;
    int north=0,south=1; me.v.push_back({cen.x,cen.y,cen.z+rz}); me.v.push_back({cen.x,cen.y,cen.z-rz});
    vector<vector<int>> ring(L-1, vector<int>(U)); const double PI=acos(-1.0);
    for(int i=1;i<L;i++){ double th=PI*i/L, st=sin(th), ct=cos(th); for(int j=0;j<U;j++){ double ph=2*PI*j/U; ring[i-1][j]=me.v.size(); me.v.push_back({cen.x+rx*st*cos(ph),cen.y+ry*st*sin(ph),cen.z+rz*ct}); }}
    for(int j=0;j<U;j++) me.f.push_back({north,ring[0][j],ring[0][(j+1)%U]});
    for(int i=0;i<L-2;i++)for(int j=0;j<U;j++){int a=ring[i][j],b=ring[i][(j+1)%U],c=ring[i+1][j],d=ring[i+1][(j+1)%U]; me.f.push_back({a,c,d}); me.f.push_back({a,d,b});}
    for(int j=0;j<U;j++) me.f.push_back({ring[L-2][j],south,ring[L-2][(j+1)%U]});
    return me;
}

static Mesh gen_torus(int ax,double R,double r,int U,int W){
    U=max(12,U); W=max(8,W); Mesh me; const double PI=acos(-1.0);
    auto point=[&](double u,double v){ double a=(R+r*cos(v))*cos(u), b=(R+r*cos(v))*sin(u), h=r*sin(v); if(ax==0) return P{cen.x+h,cen.y+a,cen.z+b}; if(ax==1) return P{cen.x+b,cen.y+h,cen.z+a}; return P{cen.x+a,cen.y+b,cen.z+h}; };
    vector<vector<int>> id(U,vector<int>(W));
    for(int i=0;i<U;i++)for(int j=0;j<W;j++){id[i][j]=me.v.size(); me.v.push_back(point(2*PI*i/U,2*PI*j/W));}
    for(int i=0;i<U;i++)for(int j=0;j<W;j++){int a=id[i][j],b=id[i][(j+1)%W],c=id[(i+1)%U][j],d=id[(i+1)%U][(j+1)%W]; me.f.push_back({a,c,d}); me.f.push_back({a,d,b});}
    return me;
}


static Mesh gen_cylinder(int ax,double R,int H,int U,int S){
    H=max(1,H); U=max(12,U); S=max(1,S); Mesh me; const double PI=acos(-1.0);
    double h0=coord(mnv,ax), h1=coord(mxv,ax);
    auto pt=[&](double h,double r,double ph){ double ca=cos(ph),sa=sin(ph); if(ax==0) return P{h,cen.y+r*ca,cen.z+r*sa}; if(ax==1) return P{cen.x+r*ca,h,cen.z+r*sa}; return P{cen.x+r*ca,cen.y+r*sa,h}; };
    vector<vector<int>> side(H+1,vector<int>(U));
    for(int i=0;i<=H;i++){ double h=h0+(h1-h0)*i/H; for(int j=0;j<U;j++){side[i][j]=me.v.size(); me.v.push_back(pt(h,R,2*PI*j/U));}}
    for(int i=0;i<H;i++)for(int j=0;j<U;j++){int a=side[i][j],b=side[i][(j+1)%U],c=side[i+1][j],d=side[i+1][(j+1)%U]; me.f.push_back({a,b,d}); me.f.push_back({a,d,c});}
    auto cap=[&](bool top){
        double h=top?h1:h0; int center=me.v.size(); me.v.push_back(pt(h,0,0)); vector<vector<int>> ring(S+1,vector<int>(U));
        for(int j=0;j<U;j++) ring[S][j]=side[top?H:0][j];
        for(int s=1;s<S;s++) for(int j=0;j<U;j++){ring[s][j]=me.v.size(); me.v.push_back(pt(h,R*s/S,2*PI*j/U));}
        for(int j=0;j<U;j++){ int a=center,b=ring[1][j],c=ring[1][(j+1)%U]; if(top) me.f.push_back({a,b,c}); else me.f.push_back({a,c,b}); }
        for(int s=1;s<S;s++) for(int j=0;j<U;j++){int a=ring[s][j],b=ring[s][(j+1)%U],c=ring[s+1][j],d=ring[s+1][(j+1)%U]; if(top){me.f.push_back({a,c,d}); me.f.push_back({a,d,b});} else {me.f.push_back({a,d,c}); me.f.push_back({a,b,d});}}
    };
    cap(false); cap(true); return me;
}

static bool hausdorff_vertices_ok(const Mesh&me, KD*origKD){
    if(!validate_mesh(me)) return false;
    KD kdC(me.v);
    for(const P&p:V) if(kdC.nearest2(p)>EPS2*1.0001) return false;
    for(const P&p:me.v) if(origKD->nearest2(p)>EPS2*1.0001) return false;
    return true;
}

static bool detect_box(){
    double tol=0.006*DIAG, sum=0; int cnt=0;
    for(auto&p:V){ double d=min({fabs(p.x-mnv.x),fabs(p.x-mxv.x),fabs(p.y-mnv.y),fabs(p.y-mxv.y),fabs(p.z-mnv.z),fabs(p.z-mxv.z)}); sum+=d; if(d<tol) cnt++; }
    return cnt>0.985*N && sum/N < 0.0035*DIAG;
}
static bool detect_ellipsoid(){
    double rx=(mxv.x-mnv.x)*0.5,ry=(mxv.y-mnv.y)*0.5,rz=(mxv.z-mnv.z)*0.5; if(min({rx,ry,rz})<DIAG*1e-8) return false;
    double sum=0; int c1=0,c2=0;
    for(auto&p:V){ double s=sqrt((p.x-cen.x)*(p.x-cen.x)/(rx*rx)+(p.y-cen.y)*(p.y-cen.y)/(ry*ry)+(p.z-cen.z)*(p.z-cen.z)/(rz*rz)); double e=fabs(s-1); sum+=e; if(e<0.035)c1++; if(e<0.060)c2++; }
    return c1>0.965*N && c2>0.990*N && sum/N<0.018;
}
struct TorD{bool ok=false; int ax=2; double R=0,r=0,err=1e9;};
static TorD detect_torus(){
    TorD best; 
    for(int ax=0;ax<3;ax++){
        double sr=0; for(auto&p:V){ double h=coord(p,ax)-coord(cen,ax); int a=(ax+1)%3,b=(ax+2)%3; double u=coord(p,a)-coord(cen,a), v=coord(p,b)-coord(cen,b); (void)h; sr+=sqrt(u*u+v*v); }
        double R=sr/N, rr=0;
        for(auto&p:V){ double h=coord(p,ax)-coord(cen,ax); int a=(ax+1)%3,b=(ax+2)%3; double u=coord(p,a)-coord(cen,a), v=coord(p,b)-coord(cen,b); double d=sqrt((sqrt(u*u+v*v)-R)*(sqrt(u*u+v*v)-R)+h*h); rr+=d; }
        rr/=N; if(rr<DIAG*1e-8||R/rr<1.15) continue;
        double sum=0; int c1=0,c2=0;
        for(auto&p:V){ double h=coord(p,ax)-coord(cen,ax); int a=(ax+1)%3,b=(ax+2)%3; double u=coord(p,a)-coord(cen,a), v=coord(p,b)-coord(cen,b); double d=sqrt((sqrt(u*u+v*v)-R)*(sqrt(u*u+v*v)-R)+h*h); double e=fabs(d/rr-1); sum+=e; if(e<0.055)c1++; if(e<0.10)c2++; }
        double err=sum/N; if(c1>0.94*N && c2>0.985*N && err<0.032 && err<best.err){best={true,ax,R,rr,err};}
    }
    return best;
}


struct CylD{bool ok=false; int ax=2; double R=0,err=1e9;};
static CylD detect_cylinder(){
    CylD best; double tol=max(EPS*0.40,DIAG*0.012);
    for(int ax=0;ax<3;ax++){
        int a=(ax+1)%3,b=(ax+2)%3; double wa=coord(mxv,a)-coord(mnv,a), wb=coord(mxv,b)-coord(mnv,b), h=coord(mxv,ax)-coord(mnv,ax); if(h<DIAG*1e-8) continue;
        if(min(wa,wb)<=DIAG*1e-8 || fabs(wa-wb)>0.10*max(wa,wb)) continue; double R=(wa+wb)*0.25; if(R<=DIAG*1e-8) continue;
        int ok=0,side=0,cap=0; double sum=0;
        for(auto&p:V){ double x=coord(p,a)-coord(cen,a), y=coord(p,b)-coord(cen,b), rr=sqrt(x*x+y*y); double ds=fabs(rr-R); double dc=min(fabs(coord(p,ax)-coord(mnv,ax)),fabs(coord(p,ax)-coord(mxv,ax))); double d=min(ds,dc); sum+=d; if(d<tol) ok++; if(ds<tol) side++; if(dc<tol) cap++; }
        double err=sum/(N*max(DIAG,1e-300)); if(ok>0.985*N && side>0.25*N && cap>0.05*N && err<0.006 && err<best.err) best={true,ax,R,err};
    }
    return best;
}


struct EC{
    vector<P> p; vector<Tri> f; vector<vector<int>> inc; vector<double> rad; vector<int> ver,mark; vector<char> va,fa; int aliveV,aliveF,stamp=1; double LIM,LIM2;
    struct Q{double w; int u,v,vu,vv; bool operator<(Q const&o)const{return w>o.w;}};
    priority_queue<Q> pq;
    EC(double lim):p(V),f(F),rad(N,0),ver(N,0),mark(N,0),va(N,1),fa(M,1),aliveV(N),aliveF(M),LIM(lim),LIM2(lim*lim){
        inc.assign(N,{}); for(int i=0;i<M;i++){inc[f[i].a].push_back(i);inc[f[i].b].push_back(i);inc[f[i].c].push_back(i);} }
    inline bool hasv(const Tri&t,int u)const{return t.a==u||t.b==u||t.c==u;}
    inline int opp(const Tri&t,int u,int v)const{if(t.a!=u&&t.a!=v)return t.a;if(t.b!=u&&t.b!=v)return t.b;return t.c;}
    void push_edge(int a,int b){ if(a==b||a<0||b<0||a>=N||b>=N||!va[a]||!va[b]) return; double d=normv(p[a]-p[b]); if(min(d+rad[a],d+rad[b])>LIM*1.000001) return; pq.push({d,a,b,ver[a],ver[b]}); }
    void push_face_edges(int id){ if(id<0||id>=M||!fa[id])return; Tri&t=f[id]; push_edge(t.a,t.b); push_edge(t.b,t.c); push_edge(t.c,t.a); }
    int edge_faces(int u,int v,int op[2]){
        int c=0; const vector<int>&lst = inc[u].size()<inc[v].size()?inc[u]:inc[v];
        for(int id:lst) if(id>=0&&id<M&&fa[id]){Tri&t=f[id]; if(hasv(t,u)&&hasv(t,v)){ if(c<2) op[c]=opp(t,u,v); c++; if(c>2) return c; }}
        return c;
    }
    void neigh(int u, vector<int>&out){
        out.clear(); if(++stamp==INT_MAX){fill(mark.begin(),mark.end(),0);stamp=1;}
        for(int id:inc[u]) if(id>=0&&id<M&&fa[id]){Tri&t=f[id]; if(!hasv(t,u)) continue; int a[3]={t.a,t.b,t.c}; for(int x:a) if(x!=u&&x>=0&&x<N&&va[x]&&mark[x]!=stamp){mark[x]=stamp; out.push_back(x);}}
    }
    bool link_ok(int u,int v,int op0,int op1){
        static vector<int> a,b; neigh(u,a); neigh(v,b); sort(a.begin(),a.end()); sort(b.begin(),b.end()); int i=0,j=0,cm=0;
        while(i<(int)a.size()&&j<(int)b.size()){ if(a[i]<b[j]) i++; else if(b[j]<a[i]) j++; else { int x=a[i]; if(!((x==op0)||(x==op1))) return false; cm++; i++; j++; } }
        return cm==2;
    }
    bool normals_ok(int u,int v){
        double tiny=max(1e-30,DIAG*DIAG*1e-24);
        for(int id:inc[u]) if(id>=0&&id<M&&fa[id]){ Tri&t=f[id]; if(!hasv(t,u)||hasv(t,v)) continue; int a=t.a,b=t.b,c=t.c; P oldn=crossp(p[b]-p[a],p[c]-p[a]); if(a==u)a=v; if(b==u)b=v; if(c==u)c=v; if(a==b||b==c||c==a) return false; P newn=crossp(p[b]-p[a],p[c]-p[a]); double no=normv(oldn), nn=normv(newn); if(nn*nn<=tiny) return false; if(no>0 && dotp(oldn,newn)<-0.05*no*nn) return false; }
        return true;
    }
    bool collapse(int u,int v){
        if(u==v||u<0||v<0||u>=N||v>=N||!va[u]||!va[v])return false; double d=normv(p[u]-p[v]); if(d+rad[u]>LIM*1.000001) return false;
        int op[2]={-1,-1}; if(edge_faces(u,v,op)!=2||op[0]<0||op[1]<0||op[0]==op[1]) return false; if(!link_ok(u,v,op[0],op[1])) return false; if(!normals_ok(u,v)) return false;
        vector<int> touch; touch.reserve(inc[u].size()+8);
        for(int id:inc[u]) if(id>=0&&id<M&&fa[id]&&hasv(f[id],u)) touch.push_back(id);
        va[u]=0; aliveV--; ver[u]++; ver[v]++; rad[v]=max(rad[v],d+rad[u]);
        for(int id:touch){ if(!fa[id]) continue; Tri&t=f[id]; if(hasv(t,v)){ fa[id]=0; aliveF--; continue; } if(t.a==u)t.a=v; if(t.b==u)t.b=v; if(t.c==u)t.c=v; inc[v].push_back(id); }
        vector<int>().swap(inc[u]);
        for(int id:inc[v]) if(id>=0&&id<M&&fa[id]&&hasv(f[id],v)) push_face_edges(id);
        return true;
    }
    Mesh run(double targetRatio){
        for(int i=0;i<M;i++) push_face_edges(i);
        int target=max(4,(int)ceil(N*targetRatio)); long long pops=0, cap=max(2000000LL,(long long)N*80);
        while(aliveV>target && !pq.empty() && pops<cap){ Q q=pq.top(); pq.pop(); pops++; int u=q.u,v=q.v; if(u<0||v<0||u>=N||v>=N||!va[u]||!va[v]||ver[u]!=q.vu||ver[v]!=q.vv) continue; double d=normv(p[u]-p[v]); if(d>q.w*1.5+1e-15 && min(d+rad[u],d+rad[v])>LIM) continue; bool aok=d+rad[u]<=LIM*1.000001, bok=d+rad[v]<=LIM*1.000001; if(!aok&&!bok) continue; bool done=false; if(aok&&bok){ double ru=max(rad[v],d+rad[u]), rv=max(rad[u],d+rad[v]); if(ru<rv) done=collapse(u,v); else done=collapse(v,u); if(!done){ if(ru<rv) done=collapse(v,u); else done=collapse(u,v); } } else if(aok) done=collapse(u,v); else done=collapse(v,u); }
        Mesh me; vector<int> mp(N,-1); vector<char> used(N,0);
        for(int i=0;i<M;i++) if(fa[i]){Tri&t=f[i]; if(t.a!=t.b&&t.b!=t.c&&t.c!=t.a){used[t.a]=used[t.b]=used[t.c]=1;}}
        for(int i=0;i<N;i++) if(va[i]&&used[i]){mp[i]=me.v.size(); me.v.push_back(p[i]);}
        for(int i=0;i<M;i++) if(fa[i]){Tri t=f[i]; if(t.a==t.b||t.b==t.c||t.c==t.a) continue; if(mp[t.a]<0||mp[t.b]<0||mp[t.c]<0) continue; me.f.push_back({mp[t.a],mp[t.b],mp[t.c]});}
        return me;
    }
};

static bool edge_collapse_candidate(double ratio, Mesh&out){
    if(N<1000) return false;
    EC ec(EPS*0.94); Mesh me=ec.run(ratio);
    if(me.v.size()>=V.size()||me.v.size()<4||me.f.size()<4) return false;
    if(!validate_mesh(me)) return false;
    out=std::move(me); return true;
}


struct Cell{long long x,y,z; bool operator==(Cell const&o)const{return x==o.x&&y==o.y&&z==o.z;}};
struct CellHash{size_t operator()(Cell const&c)const{uint64_t h=H64::smix((uint64_t)c.x); h^=H64::smix((uint64_t)c.y+0x9e3779b97f4a7c15ULL); h^=H64::smix((uint64_t)c.z+0xbf58476d1ce4e5b9ULL); return h;}};
struct Cl{double sx=0,sy=0,sz=0,best=1e300,maxd=0; int cnt=0,rep=-1;};

static bool cluster_candidate(double s,double minRatio,Mesh&out){
    if(s<=0) return false;
    unordered_map<Cell,int,CellHash> mp; mp.reserve((size_t)(N*min(0.8,max(0.02,minRatio*2)))+1000);
    vector<Cl> cl; vector<int> bel(N);
    for(int i=0;i<N;i++){
        Cell c{(long long)floor((V[i].x-mnv.x)/s),(long long)floor((V[i].y-mnv.y)/s),(long long)floor((V[i].z-mnv.z)/s)};
        auto it=mp.find(c); int id;
        if(it==mp.end()){id=cl.size(); mp.emplace(c,id); cl.push_back(Cl());}
        else id=it->second;
        bel[i]=id; auto &q=cl[id]; q.sx+=V[i].x; q.sy+=V[i].y; q.sz+=V[i].z; q.cnt++;
    }
    int K=cl.size(); if(K<4||K>N||K<max(4,(int)(N*minRatio*0.78))) return false;
    for(auto &q:cl){q.sx/=q.cnt; q.sy/=q.cnt; q.sz/=q.cnt;}
    for(int i=0;i<N;i++){auto &q=cl[bel[i]]; double dx=V[i].x-q.sx,dy=V[i].y-q.sy,dz=V[i].z-q.sz,d=dx*dx+dy*dy+dz*dz; if(d<q.best){q.best=d; q.rep=i;}}
    for(int i=0;i<N;i++){auto &q=cl[bel[i]]; double d=dist2(V[i],V[q.rep]); if(d>q.maxd) q.maxd=d;}
    for(auto&q:cl) if(q.maxd>EPS2*1.0001) return false;
    Mesh me; me.v.resize(K); for(int i=0;i<K;i++) me.v[i]=V[cl[i].rep];
    unordered_set<uint64_t,H64> seen; seen.reserve(F.size()); double tiny=max(1e-30,DIAG*DIAG*1e-24);
    for(auto&t:F){int a=bel[t.a],b=bel[t.b],c=bel[t.c]; if(a==b||b==c||c==a) continue; int x=a,y=b,z=c; if(x>y)swap(x,y); if(y>z)swap(y,z); if(x>y)swap(x,y); uint64_t key=((uint64_t)x<<42)^((uint64_t)y<<21)^(uint64_t)z; if(!seen.insert(key).second) continue; P cr=crossp(me.v[b]-me.v[a],me.v[c]-me.v[a]); if(norm2(cr)<=tiny) continue; me.f.push_back({a,b,c});}
    if((int)me.v.size()>N || me.f.empty()) return false;
    if(!validate_mesh(me)) return false;
    out = std::move(me); return true;
}

int main(){
    ios::sync_with_stdio(false); cin.tie(nullptr);
    if(!(cin>>N>>M)) return 0; V.resize(N); F.resize(M);
    for(int i=0;i<N;i++) cin>>V[i].x>>V[i].y>>V[i].z;
    for(int i=0;i<M;i++){cin>>F[i].a>>F[i].b>>F[i].c; --F[i].a;--F[i].b;--F[i].c;}
    mnv=mxv=V[0];
    for(auto&p:V){mnv.x=min(mnv.x,p.x);mnv.y=min(mnv.y,p.y);mnv.z=min(mnv.z,p.z);mxv.x=max(mxv.x,p.x);mxv.y=max(mxv.y,p.y);mxv.z=max(mxv.z,p.z);} 
    cen={(mnv.x+mxv.x)*0.5,(mnv.y+mxv.y)*0.5,(mnv.z+mxv.z)*0.5}; DIAG=normv(mxv-mnv); if(DIAG<=0){print_original();return 0;} EPS=0.049*DIAG; EPS2=EPS*EPS;
    AREA=0; for(auto&t:F) if(t.a>=0&&t.b>=0&&t.c>=0&&t.a<N&&t.b<N&&t.c<N) AREA+=0.5*normv(crossp(V[t.b]-V[t.a],V[t.c]-V[t.a]));

    bool tryPrim = N>=2000;
    if(tryPrim){
        vector<Mesh> cand;
        bool isBox = detect_box();
        bool isEll = detect_ellipsoid();
        TorD td = detect_torus();
        if(isBox){
            for(double ff: {1.70,1.55,1.42,1.30,1.18}){ int nx=(int)ceil((mxv.x-mnv.x)/(EPS*ff)), ny=(int)ceil((mxv.y-mnv.y)/(EPS*ff)), nz=(int)ceil((mxv.z-mnv.z)/(EPS*ff)); nx=min(80,max(1,nx)); ny=min(80,max(1,ny)); nz=min(80,max(1,nz)); cand.push_back(gen_box(nx,ny,nz)); }
        }
        if(isEll){
            double r=max({mxv.x-mnv.x,mxv.y-mnv.y,mxv.z-mnv.z})*0.5; for(double ff: {1.70,1.55,1.40,1.25,1.12}){ int U=(int)ceil(2*acos(-1.0)*r/(EPS*ff)); int L=(int)ceil(acos(-1.0)*r/(EPS*ff)); U=min(120,max(14,U)); L=min(64,max(7,L)); cand.push_back(gen_ellipsoid(L,U)); }
        }
        if(td.ok){ for(double ff: {1.70,1.55,1.40,1.25,1.12}){ int U=(int)ceil(2*acos(-1.0)*(td.R+td.r)/(EPS*ff)); int W=(int)ceil(2*acos(-1.0)*td.r/(EPS*ff)); U=min(180,max(18,U)); W=min(80,max(8,W)); cand.push_back(gen_torus(td.ax,td.R,td.r,U,W)); } }
        CylD cd=detect_cylinder();
        if(cd.ok){ double h=coord(mxv,cd.ax)-coord(mnv,cd.ax); for(double ff: {1.70,1.55,1.40,1.25,1.12}){ int U=(int)ceil(2*acos(-1.0)*cd.R/(EPS*ff)); int H=(int)ceil(h/(EPS*ff)); int S=(int)ceil(cd.R/(EPS*ff)); U=min(160,max(14,U)); H=min(100,max(1,H)); S=min(50,max(1,S)); cand.push_back(gen_cylinder(cd.ax,cd.R,H,U,S)); } }
        if(!cand.empty()){
            KD origKD(V); int best=-1; size_t bv=(size_t)N+1;
            for(int i=0;i<(int)cand.size();i++) if(cand[i].v.size()<bv && hausdorff_vertices_ok(cand[i],&origKD)){best=i; bv=cand[i].v.size();}
            if(best>=0){print_mesh(cand[best]); return 0;}
        }
    }

    double minRatio = (N>250000?0.070:(N>70000?0.080:(N>20000?0.085:(N>8000?0.095:0.11))));
    Mesh best; bool have=false;
    { Mesh me; if(edge_collapse_candidate(minRatio,me)){ best=std::move(me); have=true; } }
    double Ktar=max(100.0,N*minRatio); double base=sqrt(max(AREA,DIAG*DIAG*1e-12)/Ktar);
    vector<double> mul={1.35,1.22,1.10,1.0,0.90,0.80,0.70,0.60,0.50,0.42,0.35};
    for(double m:mul){ Mesh me; if(cluster_candidate(base*m,minRatio,me)){ if(!have||me.v.size()<best.v.size()){best = std::move(me); have=true;} if(have && best.v.size()<=N*minRatio*1.03) break; }}
    if(have){print_mesh(best); return 0;}
    print_original();
    return 0;
}
