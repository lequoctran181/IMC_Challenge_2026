#include <bits/stdc++.h>
using namespace std;

struct P{double x,y,z;};
static inline P operator+(P a,P b){return {a.x+b.x,a.y+b.y,a.z+b.z};}
static inline P operator-(P a,P b){return {a.x-b.x,a.y-b.y,a.z-b.z};}
static inline P operator*(P a,double s){return {a.x*s,a.y*s,a.z*s};}
static inline P operator/(P a,double s){return {a.x/s,a.y/s,a.z/s};}
static inline double dotp(P a,P b){return a.x*b.x+a.y*b.y+a.z*b.z;}
static inline P crossp(P a,P b){return {a.y*b.z-a.z*b.y,a.z*b.x-a.x*b.z,a.x*b.y-a.y*b.x};}
static inline double n2(P a){return dotp(a,a);} 
static inline double norm(P a){return sqrt(n2(a));}
static inline P unit(P a){double l=norm(a); return l>0?a/l:P{0,0,0};}

struct F{int a,b,c;};
struct Mesh{vector<P> p; vector<F> f; vector<double> r;};
struct EdgeRec{int a,b,o,fi;};
struct EdgeInfo{int a,b,o1,o2,f1,f2;};
struct Cand{float cost; int drop,keep; double nr;};

static uint64_t mix64(uint64_t x){
    x += 0x9e3779b97f4a7c15ULL;
    x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
    x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
    return x ^ (x >> 31);
}
static uint64_t key3i(int a,int b,int c){
    uint64_t x=(uint32_t)a, y=(uint32_t)b, z=(uint32_t)c;
    return mix64(x) ^ (mix64(y)<<1) ^ (mix64(z)<<2);
}
static uint64_t key3ll(long long a,long long b,long long c){
    return mix64((uint64_t)(a+10000019LL)) ^ (mix64((uint64_t)(b+20000033LL))<<1) ^ (mix64((uint64_t)(c+30000091LL))<<2);
}

static bool canonicalFace(F t, array<int,3>& q){
    if(t.a==t.b||t.b==t.c||t.c==t.a) return false;
    q={t.a,t.b,t.c}; sort(q.begin(),q.end()); return true;
}

static void compactMesh(Mesh& m){
    int n=m.p.size(); vector<int> used(n,0);
    vector<F> nf; nf.reserve(m.f.size());
    for(auto t:m.f){
        if(t.a==t.b||t.b==t.c||t.c==t.a) continue;
        nf.push_back(t); used[t.a]=used[t.b]=used[t.c]=1;
    }
    vector<int> mp(n,-1); vector<P> np; vector<double> nr;
    np.reserve(n); nr.reserve(n);
    for(int i=0;i<n;i++) if(used[i]){ mp[i]=(int)np.size(); np.push_back(m.p[i]); nr.push_back(m.r.empty()?0:m.r[i]); }
    for(auto &t:nf){ t.a=mp[t.a]; t.b=mp[t.b]; t.c=mp[t.c]; }
    m.p.swap(np); m.r.swap(nr); m.f.swap(nf);
}

static bool buildEdges(const Mesh& m, vector<EdgeInfo>* out=nullptr){
    vector<EdgeRec> e; e.reserve(m.f.size()*3);
    for(int i=0;i<(int)m.f.size();++i){
        F t=m.f[i];
        int a=t.a,b=t.b,c=t.c;
        if(a==b||b==c||c==a) return false;
        auto add=[&](int u,int v,int o){ if(u>v) swap(u,v); e.push_back({u,v,o,i}); };
        add(a,b,c); add(b,c,a); add(c,a,b);
    }
    sort(e.begin(),e.end(),[](const EdgeRec&x,const EdgeRec&y){return x.a==y.a?x.b<y.b:x.a<y.a;});
    if(out) out->clear(), out->reserve(e.size()/2);
    for(size_t i=0;i<e.size();){
        size_t j=i+1; while(j<e.size()&&e[j].a==e[i].a&&e[j].b==e[i].b) ++j;
        if(j-i!=2) return false;
        if(out) out->push_back({e[i].a,e[i].b,e[i].o,e[i+1].o,e[i].fi,e[i+1].fi});
        i=j;
    }
    return true;
}

static bool validateMesh(const Mesh& m){
    int n=m.p.size(); if(n<4||m.f.size()<4) return false;
    unordered_set<uint64_t> seen; seen.reserve(m.f.size()*2+10);
    double box=0; P mn=m.p[0],mx=m.p[0];
    for(auto v:m.p){mn.x=min(mn.x,v.x);mn.y=min(mn.y,v.y);mn.z=min(mn.z,v.z);mx.x=max(mx.x,v.x);mx.y=max(mx.y,v.y);mx.z=max(mx.z,v.z);} box=norm(mx-mn); double amin=max(1e-28, box*box*1e-24);
    for(auto t:m.f){
        if(t.a<0||t.a>=n||t.b<0||t.b>=n||t.c<0||t.c>=n) return false;
        array<int,3> q; if(!canonicalFace(t,q)) return false;
        uint64_t k=key3i(q[0],q[1],q[2]); if(seen.count(k)) return false; seen.insert(k);
        if(n2(crossp(m.p[t.b]-m.p[t.a],m.p[t.c]-m.p[t.a]))<=amin) return false;
    }
    return buildEdges(m,nullptr);
}

struct Grid{
    double h, inv; P base; unordered_map<uint64_t, vector<int>> cell; const vector<P>* pts;
    Grid(const vector<P>& p,double hh):h(hh),inv(1.0/hh),pts(&p){
        if(p.empty()) return;
        base=p[0]; P mx=p[0];
        for(auto v:p){base.x=min(base.x,v.x);base.y=min(base.y,v.y);base.z=min(base.z,v.z);mx.x=max(mx.x,v.x);mx.y=max(mx.y,v.y);mx.z=max(mx.z,v.z);} 
        cell.reserve(p.size()*2+10);
        for(int i=0;i<(int)p.size();++i){ auto [a,b,c]=idx(p[i]); cell[key3ll(a,b,c)].push_back(i); }
    }
    tuple<long long,long long,long long> idx(P p) const{
        return {(long long)floor((p.x-base.x)*inv),(long long)floor((p.y-base.y)*inv),(long long)floor((p.z-base.z)*inv)};
    }
    bool hasNear(P q,double eps2) const{
        auto [a,b,c]=idx(q);
        for(long long dx=-1;dx<=1;dx++) for(long long dy=-1;dy<=1;dy++) for(long long dz=-1;dz<=1;dz++){
            auto it=cell.find(key3ll(a+dx,b+dy,c+dz)); if(it==cell.end()) continue;
            for(int id:it->second) if(n2((*pts)[id]-q)<=eps2) return true;
        }
        return false;
    }
};
static bool hausdorffOK(const vector<P>& A,const vector<P>& B,double eps){
    if(A.empty()||B.empty()) return false;
    double e2=eps*eps;
    Grid gb(B,eps); for(auto q:A) if(!gb.hasNear(q,e2)) return false;
    Grid ga(A,eps); for(auto q:B) if(!ga.hasNear(q,e2)) return false;
    return true;
}

static int eulerChi(const Mesh& m){ vector<EdgeInfo> ed; if(!buildEdges(m,&ed)) return -999999; return (int)m.p.size()-(int)ed.size()+(int)m.f.size(); }

static Mesh makeIcoSphere(P c,double R,int sub){
    Mesh m; double t=(1.0+sqrt(5.0))/2.0;
    vector<P> v={{-1,t,0},{1,t,0},{-1,-t,0},{1,-t,0},{0,-1,t},{0,1,t},{0,-1,-t},{0,1,-t},{t,0,-1},{t,0,1},{-t,0,-1},{-t,0,1}};
    for(auto &x:v){x=unit(x); m.p.push_back(c+x*R);} 
    int ff[20][3]={{0,11,5},{0,5,1},{0,1,7},{0,7,10},{0,10,11},{1,5,9},{5,11,4},{11,10,2},{10,7,6},{7,1,8},{3,9,4},{3,4,2},{3,2,6},{3,6,8},{3,8,9},{4,9,5},{2,4,11},{6,2,10},{8,6,7},{9,8,1}};
    for(auto &q:ff) m.f.push_back({q[0],q[1],q[2]});
    for(int s=0;s<sub;s++){
        unordered_map<unsigned long long,int> mid; mid.reserve(m.f.size()*2);
        auto getmid=[&](int a,int b)->int{
            int x=min(a,b),y=max(a,b); unsigned long long k=((unsigned long long)(unsigned int)x<<32)|(unsigned int)y;
            auto it=mid.find(k); if(it!=mid.end()) return it->second;
            P dir=unit(((m.p[a]-c)+(m.p[b]-c))*0.5); int id=m.p.size(); m.p.push_back(c+dir*R); mid[k]=id; return id;
        };
        vector<F> nf; nf.reserve(m.f.size()*4);
        for(auto f:m.f){int a=f.a,b=f.b,cx=f.c; int ab=getmid(a,b),bc=getmid(b,cx),ca=getmid(cx,a); nf.push_back({a,ab,ca}); nf.push_back({b,bc,ab}); nf.push_back({cx,ca,bc}); nf.push_back({ab,bc,ca});}
        m.f.swap(nf);
    }
    m.r.assign(m.p.size(),0); return m;
}


static bool tryCluster(const Mesh& in, Mesh& out, double eps){
    int n=in.p.size(); if(n<100) return false;
    P mn=in.p[0],mx=in.p[0];
    for(auto p:in.p){mn.x=min(mn.x,p.x);mn.y=min(mn.y,p.y);mn.z=min(mn.z,p.z);mx.x=max(mx.x,p.x);mx.y=max(mx.y,p.y);mx.z=max(mx.z,p.z);} 
    vector<double> facs={0.88,0.82,0.74,0.66,0.58,0.50,0.42};
    vector<P> origN(in.f.size());
    for(int i=0;i<(int)in.f.size();++i){F t=in.f[i]; origN[i]=unit(crossp(in.p[t.b]-in.p[t.a],in.p[t.c]-in.p[t.a]));}
    for(double fac:facs){
        double h=eps*fac; if(h<=0) continue; double inv=1.0/h;
        unordered_map<uint64_t,int> id; id.reserve(n*2+10);
        vector<array<long long,3>> keys; keys.reserve(n);
        vector<int> cid(n);
        struct C{P sum; int cnt,rep; double best;}; vector<C> cs;
        for(int i=0;i<n;i++){
            long long a=(long long)floor((in.p[i].x-mn.x)*inv), b=(long long)floor((in.p[i].y-mn.y)*inv), c=(long long)floor((in.p[i].z-mn.z)*inv);
            uint64_t k=key3ll(a,b,c); auto it=id.find(k); int z;
            if(it==id.end()){z=cs.size(); id[k]=z; cs.push_back({{0,0,0},0,-1,1e100}); keys.push_back({a,b,c});}
            else z=it->second;
            cid[i]=z; cs[z].sum=cs[z].sum+in.p[i]; cs[z].cnt++;
        }
        if(cs.size()>=in.p.size()) continue;
        for(auto &c:cs) c.sum=c.sum/(double)c.cnt;
        for(int i=0;i<n;i++){C &c=cs[cid[i]]; double d=n2(in.p[i]-c.sum); if(d<c.best){c.best=d; c.rep=i;}}
        bool cover=true; double e2=eps*eps*0.9604;
        for(int i=0;i<n&&cover;i++) if(n2(in.p[i]-in.p[cs[cid[i]].rep])>e2) cover=false;
        if(!cover) continue;
        Mesh m; m.p.reserve(cs.size()); m.r.assign(cs.size(),0);
        vector<int> newId(cs.size());
        for(int i=0;i<(int)cs.size();i++){newId[i]=i; m.p.push_back(in.p[cs[i].rep]);}
        unordered_set<uint64_t> seen; seen.reserve(in.f.size()*2+10); vector<unsigned char> used(cs.size(),0);
        bool bad=false; m.f.reserve(in.f.size());
        for(int i=0;i<(int)in.f.size();i++){
            F t=in.f[i]; F q{cid[t.a],cid[t.b],cid[t.c]}; if(q.a==q.b||q.b==q.c||q.c==q.a) continue;
            P cr=crossp(m.p[q.b]-m.p[q.a],m.p[q.c]-m.p[q.a]); double ar=norm(cr); if(ar<=1e-20) {bad=true; break;} P nn=cr/ar;
            if(dotp(nn,origN[i])<0.02){bad=true; break;}
            array<int,3> can={q.a,q.b,q.c}; sort(can.begin(),can.end()); uint64_t kk=key3i(can[0],can[1],can[2]);
            if(seen.count(kk)){bad=true; break;} seen.insert(kk); m.f.push_back(q); used[q.a]=used[q.b]=used[q.c]=1;
        }
        if(bad) continue;
        for(int i=0;i<(int)used.size();i++) if(!used[i]){bad=true; break;}
        if(bad) continue;
        if(validateMesh(m)){out=std::move(m); return true;}
    }
    return false;
}

static bool trySphere(const Mesh& in, Mesh& out, double eps){
    if(eulerChi(in)!=2 || in.p.size()<80) return false;
    P c{0,0,0}; for(auto p:in.p)c=c+p; c=c/(double)in.p.size();
    double R=0; for(auto p:in.p)R+=norm(p-c); R/=in.p.size(); if(R<=eps) return false;
    double mx=0,ss=0; for(auto p:in.p){double d=fabs(norm(p-c)-R); mx=max(mx,d); ss+=d*d;}
    double rms=sqrt(ss/in.p.size()); if(mx>eps*0.70 || rms>eps*0.22) return false;
    for(int sub=1;sub<=6;sub++){
        Mesh cand=makeIcoSphere(c,R,sub);
        if(cand.p.size()>=in.p.size()) return false;
        if(hausdorffOK(in.p,cand.p,eps*0.96) && validateMesh(cand)){ out=std::move(cand); return true; }
    }
    return false;
}

static void jacobi3(double A[3][3], double V[3][3]){
    for(int i=0;i<3;i++)for(int j=0;j<3;j++)V[i][j]=(i==j);
    for(int it=0;it<32;it++){
        int p=0,q=1; double best=fabs(A[0][1]);
        if(fabs(A[0][2])>best) best=fabs(A[0][2]),p=0,q=2;
        if(fabs(A[1][2])>best) best=fabs(A[1][2]),p=1,q=2;
        if(best<1e-14) break;
        double phi=0.5*atan2(2*A[p][q],A[q][q]-A[p][p]); double c=cos(phi),s=sin(phi);
        double app=c*c*A[p][p]-2*s*c*A[p][q]+s*s*A[q][q];
        double aqq=s*s*A[p][p]+2*s*c*A[p][q]+c*c*A[q][q];
        A[p][q]=A[q][p]=0; A[p][p]=app; A[q][q]=aqq;
        for(int k=0;k<3;k++) if(k!=p&&k!=q){double akp=A[k][p],akq=A[k][q]; A[k][p]=A[p][k]=c*akp-s*akq; A[k][q]=A[q][k]=s*akp+c*akq;}
        for(int k=0;k<3;k++){double vkp=V[k][p],vkq=V[k][q]; V[k][p]=c*vkp-s*vkq; V[k][q]=s*vkp+c*vkq;}
    }
}

static Mesh makeIcoEllipsoid(P cen,P e1,P e2,P e3,double a,double b,double c,int sub){
    Mesh u=makeIcoSphere({0,0,0},1.0,sub);
    for(auto &p:u.p){P d=p; p=cen+e1*(a*d.x)+e2*(b*d.y)+e3*(c*d.z);} return u;
}
static bool tryEllipsoid(const Mesh& in, Mesh& out, double eps){
    if(eulerChi(in)!=2 || in.p.size()<120) return false;
    P mean{0,0,0}; for(auto p:in.p) mean=mean+p; mean=mean/(double)in.p.size();
    double A[3][3]={{0,0,0},{0,0,0},{0,0,0}};
    for(auto p:in.p){P q=p-mean; double x[3]={q.x,q.y,q.z}; for(int i=0;i<3;i++)for(int j=0;j<3;j++)A[i][j]+=x[i]*x[j];}
    double V[3][3]; jacobi3(A,V); int id[3]={0,1,2}; sort(id,id+3,[&](int a,int b){return A[a][a]>A[b][b];});
    auto col=[&](int k){return unit(P{V[0][k],V[1][k],V[2][k]});};
    P e1=col(id[0]), e2=col(id[1]); e2=unit(e2-e1*dotp(e1,e2)); P e3=unit(crossp(e1,e2));
    double mnv[3]={1e100,1e100,1e100}, mxv[3]={-1e100,-1e100,-1e100};
    for(auto p:in.p){P q=p-mean; double x[3]={dotp(q,e1),dotp(q,e2),dotp(q,e3)}; for(int k=0;k<3;k++){mnv[k]=min(mnv[k],x[k]); mxv[k]=max(mxv[k],x[k]);}}
    double mid[3],rad[3]; for(int k=0;k<3;k++){mid[k]=(mnv[k]+mxv[k])*0.5; rad[k]=(mxv[k]-mnv[k])*0.5; if(rad[k]<eps*1.5) return false;}
    P cen=mean+e1*mid[0]+e2*mid[1]+e3*mid[2];
    double ss=0,mx=0,scale=(rad[0]+rad[1]+rad[2])/3.0;
    for(auto p:in.p){P q=p-cen; double x=dotp(q,e1)/rad[0],y=dotp(q,e2)/rad[1],z=dotp(q,e3)/rad[2]; double dev=fabs(sqrt(x*x+y*y+z*z)-1.0)*scale; ss+=dev*dev; mx=max(mx,dev);} 
    if(mx>eps*0.80 || sqrt(ss/in.p.size())>eps*0.25) return false;
    for(int sub=1;sub<=6;sub++){Mesh cand=makeIcoEllipsoid(cen,e1,e2,e3,rad[0],rad[1],rad[2],sub); if(cand.p.size()>=in.p.size()) return false; if(hausdorffOK(in.p,cand.p,eps*0.96)&&validateMesh(cand)){out=std::move(cand); return true;}}
    return false;
}

static Mesh makeTorus(P c,P e1,P e2,P ez,double R,double r,int M,int N){
    Mesh m; m.p.reserve(M*N); m.f.reserve(M*N*2);
    const double PI=acos(-1.0);
    for(int i=0;i<M;i++){double u=2*PI*i/M,cu=cos(u),su=sin(u); for(int j=0;j<N;j++){double v=2*PI*j/N,cv=cos(v),sv=sin(v); m.p.push_back(c+e1*((R+r*cv)*cu)+e2*((R+r*cv)*su)+ez*(r*sv));}}
    for(int i=0;i<M;i++){int ni=(i+1)%M; for(int j=0;j<N;j++){int nj=(j+1)%N; int a=i*N+j,b=ni*N+j,cx=ni*N+nj,d=i*N+nj; m.f.push_back({a,b,d}); m.f.push_back({b,cx,d});}}
    m.r.assign(m.p.size(),0); return m;
}
static bool tryTorus(const Mesh& in, Mesh& out, double eps){
    if(eulerChi(in)!=0 || in.p.size()<200) return false;
    P c{0,0,0}; for(auto p:in.p)c=c+p; c=c/(double)in.p.size();
    double A[3][3]={{0,0,0},{0,0,0},{0,0,0}};
    for(auto p:in.p){P q=p-c; double x[3]={q.x,q.y,q.z}; for(int i=0;i<3;i++)for(int j=0;j<3;j++)A[i][j]+=x[i]*x[j];}
    double V[3][3]; jacobi3(A,V); int id[3]={0,1,2}; sort(id,id+3,[&](int a,int b){return A[a][a]<A[b][b];});
    auto col=[&](int k){return unit(P{V[0][k],V[1][k],V[2][k]});};
    P ez=col(id[0]), e1=col(id[2]), e2=unit(crossp(ez,e1)); e1=unit(crossp(e2,ez));
    double R=0; vector<double> rho; rho.reserve(in.p.size());
    for(auto p:in.p){P q=p-c; double x=dotp(q,e1), y=dotp(q,e2); double rr=sqrt(x*x+y*y); rho.push_back(rr); R+=rr;} R/=in.p.size();
    double r=0; for(int i=0;i<(int)in.p.size();i++){P q=in.p[i]-c; double z=dotp(q,ez); r+=sqrt((rho[i]-R)*(rho[i]-R)+z*z);} r/=in.p.size();
    if(R<eps*2||r<eps*0.5||R<r*1.15) return false;
    double mx=0,ss=0; for(int i=0;i<(int)in.p.size();i++){P q=in.p[i]-c; double z=dotp(q,ez); double d=fabs(sqrt((rho[i]-R)*(rho[i]-R)+z*z)-r); mx=max(mx,d); ss+=d*d;}
    double rms=sqrt(ss/in.p.size()); if(mx>eps*0.75||rms>eps*0.24) return false;
    const double PI=acos(-1.0); int M=max(12,min(256,(int)ceil(2*PI*R/(eps*0.72)))); int N=max(8,min(160,(int)ceil(2*PI*r/(eps*0.72))));
    for(int add=0;add<3;add++){
        Mesh cand=makeTorus(c,e1,e2,ez,R,r,M+(add?M/4*add:0),N+(add?N/4*add:0));
        if(cand.p.size()>=in.p.size()) return false;
        if(hausdorffOK(in.p,cand.p,eps*0.96) && validateMesh(cand)){out=std::move(cand); return true;}
    }
    return false;
}

struct PassData{
    vector<EdgeInfo> edges; vector<int> np,adj,fp,vf; vector<P> fn; vector<double> area,sal;
};
static bool prepare(const Mesh& m, PassData& d){
    if(!buildEdges(m,&d.edges)) return false;
    int n=m.p.size();
    d.np.assign(n+1,0); for(auto &e:d.edges){d.np[e.a+1]++; d.np[e.b+1]++;}
    for(int i=1;i<=n;i++) d.np[i]+=d.np[i-1];
    d.adj.assign(d.np[n],0); vector<int> cur=d.np;
    for(auto &e:d.edges){d.adj[cur[e.a]++]=e.b; d.adj[cur[e.b]++]=e.a;}
    d.fp.assign(n+1,0); for(auto f:m.f){d.fp[f.a+1]++; d.fp[f.b+1]++; d.fp[f.c+1]++;}
    for(int i=1;i<=n;i++) d.fp[i]+=d.fp[i-1];
    d.vf.assign(d.fp[n],0); cur=d.fp;
    d.fn.resize(m.f.size()); d.area.resize(m.f.size());
    for(int i=0;i<(int)m.f.size();i++){
        auto t=m.f[i]; P cr=crossp(m.p[t.b]-m.p[t.a],m.p[t.c]-m.p[t.a]); double ar=norm(cr); if(ar==0) return false; d.fn[i]=cr/ar; d.area[i]=ar;
        d.vf[cur[t.a]++]=i; d.vf[cur[t.b]++]=i; d.vf[cur[t.c]++]=i;
    }
    d.sal.assign(n,0); vector<P> av(n,{0,0,0});
    for(int i=0;i<(int)m.f.size();i++){auto t=m.f[i]; av[t.a]=av[t.a]+d.fn[i]; av[t.b]=av[t.b]+d.fn[i]; av[t.c]=av[t.c]+d.fn[i];}
    for(int v=0;v<n;v++){
        P u=unit(av[v]); double mx=0; for(int k=d.fp[v];k<d.fp[v+1];k++) mx=max(mx,1.0-dotp(u,d.fn[d.vf[k]]));
        d.sal[v]=min(1.0,mx/0.18);
    }
    return true;
}

static F repl(F t,int a,int b){ if(t.a==a)t.a=b; if(t.b==a)t.b=b; if(t.c==a)t.c=b; return t; }

static bool linkOK(const PassData& d, int u,int v,int o1,int o2, vector<int>& mark, int token){
    for(int k=d.np[u];k<d.np[u+1];k++) mark[d.adj[k]]=token;
    int cnt=0; bool h1=false,h2=false;
    for(int k=d.np[v];k<d.np[v+1];k++) if(mark[d.adj[k]]==token){cnt++; if(d.adj[k]==o1)h1=true; if(d.adj[k]==o2)h2=true;}
    return cnt==2 && h1 && h2;
}

static bool evalCollapse(const Mesh& m,const PassData& d,int drop,int keep,int f1,int f2,double eps,double aggr,Cand& out){
    double len=norm(m.p[drop]-m.p[keep]); double s=max(d.sal[drop],d.sal[keep]);
    double lim=eps*aggr*(1.0-0.66*pow(s,0.65)); lim=max(lim,eps*0.22);
    if(s>0.92 && len>eps*0.22) return false;
    double nr=max(m.r[keep],m.r[drop]+len); if(nr>lim || nr>eps) return false;
    double pcost=0; int cnt=0; double mindot = 0.18 + 0.55*s;
    for(int k=d.fp[drop];k<d.fp[drop+1];k++){
        int fi=d.vf[k]; if(fi==f1||fi==f2) continue; F nt=repl(m.f[fi],drop,keep); if(nt.a==nt.b||nt.b==nt.c||nt.c==nt.a) return false;
        P cr=crossp(m.p[nt.b]-m.p[nt.a],m.p[nt.c]-m.p[nt.a]); double ar=norm(cr); if(ar<=1e-20) return false; P nn=cr/ar; if(dotp(nn,d.fn[fi])<mindot) return false;
        F ot=m.f[fi]; P base=m.p[ot.a]; double dis=dotp(m.p[keep]-base,d.fn[fi]); pcost += dis*dis*(0.25+d.area[fi]); cnt++;
    }
    for(int k=d.fp[keep];k<d.fp[keep+1];k++){
        int fi=d.vf[k]; if(fi==f1||fi==f2) continue; F nt=repl(m.f[fi],drop,keep); if(nt.a==nt.b||nt.b==nt.c||nt.c==nt.a) continue;
        P cr=crossp(m.p[nt.b]-m.p[nt.a],m.p[nt.c]-m.p[nt.a]); double ar=norm(cr); if(ar<=1e-20) return false; P nn=cr/ar; if(dotp(nn,d.fn[fi])<mindot) return false;
    }
    double th=eps*eps*(0.015+0.10*(1.0-s))*max(1,cnt);
    if(pcost>th) return false;
    out={(float)(pcost + len*len*(0.012+2.5*s) + nr*nr*0.04),drop,keep,nr}; return true;
}

static bool passOnce(Mesh& m,double eps,int passNo,double deadline){
    auto nowsec=[](){return chrono::duration<double>(chrono::steady_clock::now().time_since_epoch()).count();};
    if(nowsec()>deadline) return false;
    PassData d; if(!prepare(m,d)) return false; int n=m.p.size();
    P mn=m.p[0],mx=m.p[0]; for(auto p:m.p){mn.x=min(mn.x,p.x);mn.y=min(mn.y,p.y);mn.z=min(mn.z,p.z);mx.x=max(mx.x,p.x);mx.y=max(mx.y,p.y);mx.z=max(mx.z,p.z);} 
    for(int i=0;i<n;i++){
        int ec=0; if(m.p[i].x<mn.x+eps*0.25||m.p[i].x>mx.x-eps*0.25) ec++; if(m.p[i].y<mn.y+eps*0.25||m.p[i].y>mx.y-eps*0.25) ec++; if(m.p[i].z<mn.z+eps*0.25||m.p[i].z>mx.z-eps*0.25) ec++; double ex=(ec>=3?1.0:(ec==2?0.82:(ec==1?0.12:0.0))); d.sal[i]=max(d.sal[i],ex);
    }
    vector<Cand> cand; cand.reserve(d.edges.size()/3+10); vector<int> mark(n,0); int token=1;
    double aggr=min(1.0,0.58+0.055*passNo);
    for(auto &e:d.edges){
        if(token>INT_MAX-10){fill(mark.begin(),mark.end(),0); token=1;}
        if(!linkOK(d,e.a,e.b,e.o1,e.o2,mark,token++)) continue;
        Cand c1,c2; bool b1=evalCollapse(m,d,e.a,e.b,e.f1,e.f2,eps,aggr,c1); bool b2=evalCollapse(m,d,e.b,e.a,e.f1,e.f2,eps,aggr,c2);
        if(b1&&b2) cand.push_back(c1.cost<c2.cost?c1:c2); else if(b1) cand.push_back(c1); else if(b2) cand.push_back(c2);
    }
    if(cand.empty()) return false;
    sort(cand.begin(),cand.end(),[](const Cand&a,const Cand&b){return a.cost<b.cost;});
    int minSel=max(8,n/200000);
    int capBase=INT_MAX;
    for(int attempt=0;attempt<6;attempt++){
        vector<unsigned char> busy(n,0),dead(n,0); vector<int> rep(n); iota(rep.begin(),rep.end(),0); vector<double> nr=m.r; int sel=0;
        for(auto &c:cand){
            if(sel>=capBase) break;
            if(dead[c.drop]||dead[c.keep]||busy[c.drop]||busy[c.keep]) continue;
            rep[c.drop]=c.keep; dead[c.drop]=1; nr[c.keep]=max(nr[c.keep],c.nr); sel++;
            auto block=[&](int v){ busy[v]=1; for(int k=d.np[v];k<d.np[v+1];k++) busy[d.adj[k]]=1; };
            block(c.drop); block(c.keep);
        }
        if(sel<minSel) return false;
        Mesh nm; nm.p=m.p; nm.r=nr; nm.f.reserve(m.f.size()>2*(size_t)sel?m.f.size()-2*(size_t)sel:m.f.size());
        for(auto t:m.f){t.a=rep[t.a];t.b=rep[t.b];t.c=rep[t.c]; if(t.a==t.b||t.b==t.c||t.c==t.a) continue; nm.f.push_back(t);} 
        compactMesh(nm); if(validateMesh(nm)){m=std::move(nm); return true;}
        capBase = sel/2;
    }
    return false;
}

static void simplify(Mesh& m,double eps){
    auto start=chrono::steady_clock::now(); double st=chrono::duration<double>(start.time_since_epoch()).count(); double deadline=st+19.2;
    if(m.r.size()!=m.p.size()) m.r.assign(m.p.size(),0);
    for(int pass=0;pass<80;pass++){
        int before=m.p.size(); Mesh old=m;
        if(!passOnce(m,eps,pass,deadline)){m=std::move(old); break;}
        if(m.p.size()>=before*0.997) break;
        double now=chrono::duration<double>(chrono::steady_clock::now().time_since_epoch()).count(); if(now>deadline) break;
    }
}

int main(){
    ios::sync_with_stdio(false); cin.tie(nullptr);
    int V,FN; if(!(cin>>V>>FN)) return 0; Mesh in; in.p.resize(V); in.r.assign(V,0); in.f.resize(FN);
    string s; for(int i=0;i<V;i++){cin>>s>>in.p[i].x>>in.p[i].y>>in.p[i].z;}
    for(int i=0;i<FN;i++){cin>>s>>in.f[i].a>>in.f[i].b>>in.f[i].c; --in.f[i].a;--in.f[i].b;--in.f[i].c;}
    P mn=in.p[0],mx=in.p[0]; for(auto p:in.p){mn.x=min(mn.x,p.x);mn.y=min(mn.y,p.y);mn.z=min(mn.z,p.z);mx.x=max(mx.x,p.x);mx.y=max(mx.y,p.y);mx.z=max(mx.z,p.z);} double diag=norm(mx-mn); double eps=max(1e-12,0.0490*diag);
    Mesh best=in, cand;
    if(V>200 && trySphere(in,cand,eps) && cand.p.size()<best.p.size()) best=cand;
    if(V>250 && tryEllipsoid(in,cand,eps) && cand.p.size()<best.p.size()) best=cand;
    if(V>400 && tryTorus(in,cand,eps) && cand.p.size()<best.p.size()) best=cand;
    if(V>300 && tryCluster(in,cand,eps) && cand.p.size()<best.p.size()) best=cand;
    Mesh gen=in; if(validateMesh(gen)){ simplify(gen,eps); if(validateMesh(gen) && gen.p.size()<best.p.size()) best=std::move(gen); }
    cout.setf(ios::fixed); cout<<setprecision(9);
    cout<<best.p.size()<<' '<<best.f.size()<<'\n';
    for(auto p:best.p) cout<<"v "<<p.x<<' '<<p.y<<' '<<p.z<<'\n';
    for(auto f:best.f) cout<<"f "<<f.a+1<<' '<<f.b+1<<' '<<f.c+1<<'\n';
    return 0;
}
