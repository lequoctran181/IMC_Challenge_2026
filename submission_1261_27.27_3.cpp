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

static uint64_t mix64(uint64_t x){
    x += 0x9e3779b97f4a7c15ULL;
    x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
    x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
    return x ^ (x >> 31);
}
static uint64_t key3ll(long long a,long long b,long long c){
    return mix64((uint64_t)(a+10000019LL)) ^ (mix64((uint64_t)(b+20000033LL))<<1) ^ (mix64((uint64_t)(c+30000091LL))<<2);
}
static uint64_t key3i(int a,int b,int c){
    return mix64((uint32_t)a) ^ (mix64((uint32_t)b)<<1) ^ (mix64((uint32_t)c)<<2);
}
static double nowsec(){return chrono::duration<double>(chrono::steady_clock::now().time_since_epoch()).count();}

static bool buildEdges(const Mesh& m, vector<EdgeInfo>* out=nullptr){
    vector<EdgeRec> e; e.reserve(m.f.size()*3);
    int n=(int)m.p.size();
    for(int i=0;i<(int)m.f.size();++i){
        F t=m.f[i];
        if(t.a<0||t.a>=n||t.b<0||t.b>=n||t.c<0||t.c>=n) return false;
        int a=t.a,b=t.b,c=t.c;
        if(a==b||b==c||c==a) return false;
        auto add=[&](int u,int v,int o){ if(u>v) swap(u,v); e.push_back({u,v,o,i}); };
        add(a,b,c); add(b,c,a); add(c,a,b);
    }
    sort(e.begin(),e.end(),[](const EdgeRec&x,const EdgeRec&y){return x.a==y.a?x.b<y.b:x.a<y.a;});
    if(out){out->clear(); out->reserve(e.size()/2+1);}    
    for(size_t i=0;i<e.size();){
        size_t j=i+1; while(j<e.size()&&e[j].a==e[i].a&&e[j].b==e[i].b) ++j;
        if(j-i!=2) return false;
        if(out) out->push_back({e[i].a,e[i].b,e[i].o,e[i+1].o,e[i].fi,e[i+1].fi});
        i=j;
    }
    return true;
}
static void compactMesh(Mesh& m){
    int n=(int)m.p.size(); vector<unsigned char> used(n,0); vector<F> nf; nf.reserve(m.f.size());
    unordered_set<uint64_t> seen; seen.reserve(m.f.size()*2+10);
    for(auto t:m.f){
        if(t.a==t.b||t.b==t.c||t.c==t.a) continue;
        array<int,3> q={t.a,t.b,t.c}; sort(q.begin(),q.end()); uint64_t k=key3i(q[0],q[1],q[2]);
        if(seen.insert(k).second){nf.push_back(t); used[t.a]=used[t.b]=used[t.c]=1;}
    }
    vector<int> mp(n,-1); vector<P> np; vector<double> nr; np.reserve(n); nr.reserve(n);
    for(int i=0;i<n;i++) if(used[i]){mp[i]=(int)np.size(); np.push_back(m.p[i]); nr.push_back(m.r.empty()?0:m.r[i]);}
    for(auto &t:nf){t.a=mp[t.a];t.b=mp[t.b];t.c=mp[t.c];}
    m.p.swap(np); m.f.swap(nf); m.r.swap(nr);
}
static bool validateMesh(const Mesh& m){
    int n=(int)m.p.size(); if(n<4||m.f.size()<4||m.r.size()!=m.p.size()) return false;
    P mn=m.p[0],mx=m.p[0]; for(auto v:m.p){mn.x=min(mn.x,v.x);mn.y=min(mn.y,v.y);mn.z=min(mn.z,v.z);mx.x=max(mx.x,v.x);mx.y=max(mx.y,v.y);mx.z=max(mx.z,v.z);} double diag=max(1e-30,norm(mx-mn));
    unordered_set<uint64_t> seen; seen.reserve(m.f.size()*2+10);
    for(auto t:m.f){
        if(t.a<0||t.a>=n||t.b<0||t.b>=n||t.c<0||t.c>=n) return false;
        if(t.a==t.b||t.b==t.c||t.c==t.a) return false;
        array<int,3> q={t.a,t.b,t.c}; sort(q.begin(),q.end()); uint64_t k=key3i(q[0],q[1],q[2]); if(!seen.insert(k).second) return false;
        if(n2(crossp(m.p[t.b]-m.p[t.a],m.p[t.c]-m.p[t.a]))<=diag*diag*1e-26) return false;
    }
    return buildEdges(m,nullptr);
}
static int eulerChi(const Mesh& m){ vector<EdgeInfo> ed; if(!buildEdges(m,&ed)) return -1000000000; return (int)m.p.size()-(int)ed.size()+(int)m.f.size(); }

static inline unsigned long long edgeKey2(int a,int b){ if(a>b) swap(a,b); return ((unsigned long long)(unsigned int)a<<32)|(unsigned int)b; }
static bool orientMesh(Mesh& m){
    int n=(int)m.p.size(), Fm=(int)m.f.size();
    if(n<4||Fm<4) return false;
    struct Ref{int f; int dir;};
    unordered_map<unsigned long long, vector<Ref>> mp; mp.reserve((size_t)Fm*4+10);
    auto add=[&](int f,int a,int b){ int u=min(a,b),v=max(a,b); int dir=(a==u&&b==v)?1:-1; mp[((unsigned long long)(unsigned int)u<<32)|(unsigned int)v].push_back({f,dir}); };
    for(int i=0;i<Fm;i++){F t=m.f[i]; if(t.a==t.b||t.b==t.c||t.c==t.a) return false; add(i,t.a,t.b); add(i,t.b,t.c); add(i,t.c,t.a);}    
    vector<vector<pair<int,int>>> adj(Fm);
    for(auto &kv:mp){
        auto &v=kv.second; if(v.size()!=2) return false;
        // sign_g = rel * sign_f, where rel = -dir_f*dir_g
        int rel = -v[0].dir*v[1].dir;
        adj[v[0].f].push_back({v[1].f,rel});
        adj[v[1].f].push_back({v[0].f,rel});
    }
    vector<int> sg(Fm,0);
    for(int st=0;st<Fm;st++) if(!sg[st]){
        sg[st]=1; vector<int> q={st};
        for(size_t qi=0;qi<q.size();qi++){
            int u=q[qi];
            for(auto [v,rel]:adj[u]){
                int ns=sg[u]*rel;
                if(!sg[v]){sg[v]=ns; q.push_back(v);} else if(sg[v]!=ns) return false;
            }
        }
    }
    for(int i=0;i<Fm;i++) if(sg[i]<0) swap(m.f[i].b,m.f[i].c);
    return true;
}
static bool finalizeCandidate(Mesh& m){
    compactMesh(m);
    if(!validateMesh(m)) return false;
    if(!orientMesh(m)) return false;
    return validateMesh(m);
}

struct Grid{
    double h,inv; P base; unordered_map<uint64_t, vector<int>> cell; const vector<P>* pts;
    Grid(){}
    Grid(const vector<P>& p,double hh){init(p,hh);}    
    void init(const vector<P>& p,double hh){h=max(hh,1e-30);inv=1.0/h;pts=&p;cell.clear(); if(p.empty())return; base=p[0]; for(auto v:p){base.x=min(base.x,v.x);base.y=min(base.y,v.y);base.z=min(base.z,v.z);} cell.reserve(p.size()*2+10); for(int i=0;i<(int)p.size();i++){auto [a,b,c]=idx(p[i]); cell[key3ll(a,b,c)].push_back(i);} }
    tuple<long long,long long,long long> idx(P p) const{return {(long long)floor((p.x-base.x)*inv),(long long)floor((p.y-base.y)*inv),(long long)floor((p.z-base.z)*inv)};}
    template<class Fn> void eachNear(P q,Fn fn) const{
        auto [a,b,c]=idx(q);
        for(long long dx=-1;dx<=1;dx++)for(long long dy=-1;dy<=1;dy++)for(long long dz=-1;dz<=1;dz++){
            auto it=cell.find(key3ll(a+dx,b+dy,c+dz)); if(it==cell.end()) continue;
            for(int id:it->second) fn(id);
        }
    }
    bool hasNear(P q,double e2) const{bool ok=false; eachNear(q,[&](int id){if(!ok&&n2((*pts)[id]-q)<=e2) ok=true;}); return ok;}
};
static bool hausdorffOK(const vector<P>& A,const vector<P>& B,double eps){
    if(A.empty()||B.empty()) return false; double e2=eps*eps; Grid gb(B,eps); for(auto q:A) if(!gb.hasNear(q,e2)) return false; Grid ga(A,eps); for(auto q:B) if(!ga.hasNear(q,e2)) return false; return true;
}

struct Masks{int R; array<vector<unsigned char>,3> m; P mn,mx;};
static inline pair<double,double> proj2(P p,int v){
    if(v==0) return {p.x,p.y};
    if(v==1) return {p.x,p.z};
    return {p.y,p.z};
}
static inline double ed2(double ax,double ay,double bx,double by,double px,double py){return (px-ax)*(by-ay)-(py-ay)*(bx-ax);} 
static Masks rasterMesh(const Mesh& me,const P& mn,const P& mx,int R){
    Masks M; M.R=R; M.mn=mn; M.mx=mx; for(int i=0;i<3;i++) M.m[i].assign(R*R,0);
    double lo[3]={mn.x,mn.y,mn.z}, hi[3]={mx.x,mx.y,mx.z};
    double range[3]={max(hi[0]-lo[0],1e-12),max(hi[1]-lo[1],1e-12),max(hi[2]-lo[2],1e-12)};
    int a1[3]={0,0,1}, a2[3]={1,2,2};
    for(auto f:me.f){
        P p[3]={me.p[f.a],me.p[f.b],me.p[f.c]};
        for(int v=0;v<3;v++){
            int A=a1[v],B=a2[v]; double x[3],y[3];
            for(int k=0;k<3;k++){double cc[3]={p[k].x,p[k].y,p[k].z}; x[k]=(cc[A]-lo[A])/range[A]*(R-1); y[k]=(cc[B]-lo[B])/range[B]*(R-1);}            
            double area=ed2(x[0],y[0],x[1],y[1],x[2],y[2]); if(fabs(area)<1e-12) continue;
            int xmin=max(0,(int)floor(min({x[0],x[1],x[2]}))-1), xmax=min(R-1,(int)ceil(max({x[0],x[1],x[2]}))+1);
            int ymin=max(0,(int)floor(min({y[0],y[1],y[2]}))-1), ymax=min(R-1,(int)ceil(max({y[0],y[1],y[2]}))+1);
            unsigned char* mask=M.m[v].data(); bool pos=area>0;
            for(int yy=ymin;yy<=ymax;yy++){
                double py=yy+0.5;
                for(int xx=xmin;xx<=xmax;xx++){
                    double px=xx+0.5;
                    double e0=ed2(x[0],y[0],x[1],y[1],px,py), e1=ed2(x[1],y[1],x[2],y[2],px,py), e2=ed2(x[2],y[2],x[0],y[0],px,py);
                    if(pos?(e0>=-1e-9&&e1>=-1e-9&&e2>=-1e-9):(e0<=1e-9&&e1<=1e-9&&e2<=1e-9)) mask[yy*R+xx]=1;
                }
            }
        }
    }
    return M;
}
static double maskScore(const vector<unsigned char>& A,const vector<unsigned char>& B){
    int n=A.size(); double sx=0,sy=0,sxx=0,syy=0,sxy=0; int inter=0,uni=0;
    for(int i=0;i<n;i++){double x=A[i],y=B[i]; sx+=x;sy+=y;sxx+=x*x;syy+=y*y;sxy+=x*y; inter+=(A[i]&B[i]); uni+=(A[i]|B[i]);}
    if(!uni) return 1.0; double mx=sx/n,my=sy/n; double vx=sxx/n-mx*mx,vy=syy/n-my*my,cov=sxy/n-mx*my;
    const double C1=0.0001,C2=0.0009; double ssim=((2*mx*my+C1)*(2*cov+C2))/((mx*mx+my*my+C1)*(vx+vy+C2));
    double iou=(double)inter/(double)uni; return min(ssim, 0.65*ssim+0.35*iou);
}
static bool maskOK(const Masks& orig,const Mesh& cand,double th){
    Masks M=rasterMesh(cand,orig.mn,orig.mx,orig.R); for(int i=0;i<3;i++) if(maskScore(orig.m[i],M.m[i])<th) return false; return true;
}

static Mesh makeIcoSphere(P c,double R,int sub){
    Mesh m; double t=(1.0+sqrt(5.0))/2.0; vector<P> v={{-1,t,0},{1,t,0},{-1,-t,0},{1,-t,0},{0,-1,t},{0,1,t},{0,-1,-t},{0,1,-t},{t,0,-1},{t,0,1},{-t,0,-1},{-t,0,1}};
    for(auto x:v)m.p.push_back(c+unit(x)*R);
    int ff[20][3]={{0,11,5},{0,5,1},{0,1,7},{0,7,10},{0,10,11},{1,5,9},{5,11,4},{11,10,2},{10,7,6},{7,1,8},{3,9,4},{3,4,2},{3,2,6},{3,6,8},{3,8,9},{4,9,5},{2,4,11},{6,2,10},{8,6,7},{9,8,1}};
    for(auto &q:ff)m.f.push_back({q[0],q[1],q[2]});
    for(int s=0;s<sub;s++){
        unordered_map<unsigned long long,int> mid; mid.reserve(m.f.size()*2+10);
        auto gm=[&](int a,int b){int x=min(a,b),y=max(a,b); unsigned long long k=((unsigned long long)(unsigned)x<<32)|(unsigned)y; auto it=mid.find(k); if(it!=mid.end()) return it->second; P d=unit((m.p[a]-c)+(m.p[b]-c)); int id=m.p.size(); m.p.push_back(c+d*R); mid[k]=id; return id;};
        vector<F> nf; nf.reserve(m.f.size()*4); for(auto f:m.f){int a=f.a,b=f.b,cx=f.c; int ab=gm(a,b),bc=gm(b,cx),ca=gm(cx,a); nf.push_back({a,ab,ca}); nf.push_back({b,bc,ab}); nf.push_back({cx,ca,bc}); nf.push_back({ab,bc,ca});} m.f.swap(nf);
    }
    m.r.assign(m.p.size(),0); return m;
}
static bool trySphere(const Mesh& in,Mesh& out,double eps,const Masks& om){
    if(in.p.size()<80) return false; P c{0,0,0}; for(auto p:in.p)c=c+p; c=c/(double)in.p.size(); double R=0; for(auto p:in.p)R+=norm(p-c); R/=in.p.size(); if(R<=eps) return false;
    double mx=0,ss=0; for(auto p:in.p){double d=fabs(norm(p-c)-R);mx=max(mx,d);ss+=d*d;} if(mx>eps*.85||sqrt(ss/in.p.size())>eps*.28) return false;
    for(int sub=1;sub<=6;sub++){Mesh cand=makeIcoSphere(c,R,sub); if(cand.p.size()>=in.p.size()) return false; if(finalizeCandidate(cand)&&hausdorffOK(in.p,cand.p,eps*.995)&&maskOK(om,cand,.92)){out=move(cand);return true;}}
    return false;
}

struct CoverResult{vector<int> centers; vector<int> assign;};
static CoverResult greedyCover(const vector<P>& pts,double eps,double deadline){
    int n=pts.size(); CoverResult R; if(!n) return R; double e2=eps*eps; Grid g(pts,eps); vector<unsigned char> cov(n,0); int left=n; R.centers.reserve(max(4,n/100));
    auto countUncov=[&](int id){int cnt=0; P q=pts[id]; g.eachNear(q,[&](int j){if(!cov[j]&&n2(pts[j]-q)<=e2) cnt++;}); return cnt;};
    auto mark=[&](int id){int cnt=0; P q=pts[id]; g.eachNear(q,[&](int j){if(!cov[j]&&n2(pts[j]-q)<=e2){cov[j]=1;cnt++;}}); return cnt;};
    int ptr=0,loops=0;
    while(left>0){
        if(nowsec()>deadline) {while(ptr<n&&cov[ptr])ptr++; if(ptr>=n)break; R.centers.push_back(ptr); left-=mark(ptr); continue;}
        while(ptr<n&&cov[ptr]) ptr++; if(ptr>=n) break; int seed=ptr,bestSeed=seed,bestSeedCnt=-1;
        int seen=0; for(int j=ptr;j<n&&seen<16;j++) if(!cov[j]){int c=countUncov(j); if(c>bestSeedCnt){bestSeedCnt=c;bestSeed=j;} seen++;}
        seed=bestSeed;
        vector<int> neigh; neigh.reserve(256); P sp=pts[seed]; g.eachNear(sp,[&](int j){if(n2(pts[j]-sp)<=e2) neigh.push_back(j);});
        int best=seed,bestc=-1; int cap=160; int step=max(1,(int)neigh.size()/cap);
        for(int ii=0;ii<(int)neigh.size();ii+=step){int c=countUncov(neigh[ii]); if(c>bestc){bestc=c;best=neigh[ii];}}
        int cs=countUncov(seed); if(cs>bestc){best=seed;bestc=cs;}
        R.centers.push_back(best); int got=mark(best); left-=got; loops++;
        if(got==0){cov[seed]=1;left--;}
    }
    R.assign.assign(n,-1);
    if(!R.centers.empty()){
        vector<P> cp; cp.reserve(R.centers.size()); for(int id:R.centers)cp.push_back(pts[id]); Grid cg(cp,eps);
        for(int i=0;i<n;i++){double bd=1e300; int bi=0; cg.eachNear(pts[i],[&](int j){double d=n2(cp[j]-pts[i]); if(d<bd){bd=d;bi=j;}}); R.assign[i]=bi;}
    }
    return R;
}
static bool mappedCoverMesh(const Mesh& in,const CoverResult& cr,Mesh& out,const Masks& om,double eps){
    int k=cr.centers.size(); if(k<4||cr.assign.size()!=in.p.size()) return false;
    Mesh m; m.p.reserve(k); for(int id:cr.centers)m.p.push_back(in.p[id]); m.r.assign(k,0); m.f.reserve(in.f.size()); unordered_set<uint64_t> seen; seen.reserve(in.f.size()*2+10);
    for(auto t:in.f){int a=cr.assign[t.a],b=cr.assign[t.b],c=cr.assign[t.c]; if(a<0||b<0||c<0||a==b||b==c||c==a) continue; array<int,3> q={a,b,c}; sort(q.begin(),q.end()); uint64_t kk=key3i(q[0],q[1],q[2]); if(seen.insert(kk).second) m.f.push_back({a,b,c});}
    compactMesh(m); if(m.p.size()>=in.p.size()) return false; if(finalizeCandidate(m)&&hausdorffOK(in.p,m.p,eps*.999)&&maskOK(om,m,.92)){out=move(m);return true;} return false;
}
static bool makeBipyramid(const vector<P>& p,Mesh& m){
    int n=p.size(); if(n<4) return false; m=Mesh(); m.p=p; m.r.assign(n,0); if(n==4){m.f={{0,1,2},{0,3,1},{1,3,2},{2,3,0}}; return validateMesh(m);}    
    P mn=p[0],mx=p[0]; for(auto q:p){mn.x=min(mn.x,q.x);mn.y=min(mn.y,q.y);mn.z=min(mn.z,q.z);mx.x=max(mx.x,q.x);mx.y=max(mx.y,q.y);mx.z=max(mx.z,q.z);}    
    vector<array<int,2>> tries; for(int ax=0;ax<3;ax++){int lo=0,hi=0; for(int i=1;i<n;i++){double vi=ax==0?p[i].x:(ax==1?p[i].y:p[i].z), vlo=ax==0?p[lo].x:(ax==1?p[lo].y:p[lo].z), vhi=ax==0?p[hi].x:(ax==1?p[hi].y:p[hi].z); if(vi<vlo)lo=i; if(vi>vhi)hi=i;} if(lo!=hi) tries.push_back({lo,hi});}
    sort(tries.begin(),tries.end(),[&](auto A,auto B){return n2(p[A[0]]-p[A[1]])>n2(p[B[0]]-p[B[1]]);});
    for(auto tr:tries){int top=tr[0],bot=tr[1]; P ez=unit(p[top]-p[bot]); if(norm(ez)==0) continue; P tmp=fabs(ez.x)<.8?P{1,0,0}:P{0,1,0}; P e1=unit(crossp(ez,tmp)), e2=crossp(ez,e1); vector<pair<double,int>> ring; ring.reserve(n-2); P c=(p[top]+p[bot])*.5; for(int i=0;i<n;i++) if(i!=top&&i!=bot){P q=p[i]-c; ring.push_back({atan2(dotp(q,e2),dotp(q,e1)),i});} sort(ring.begin(),ring.end()); int r=ring.size(); if(r<3)continue; m.f.clear(); m.f.reserve(2*r); for(int i=0;i<r;i++){int a=ring[i].second,b=ring[(i+1)%r].second; m.f.push_back({top,a,b}); m.f.push_back({bot,b,a});} if(validateMesh(m)) return true; }
    return false;
}
static bool coverBipyramid(const Mesh& in,const CoverResult& cr,Mesh& out,const Masks& om,double eps){
    vector<P> p; p.reserve(cr.centers.size()); for(int id:cr.centers)p.push_back(in.p[id]); if(p.size()<4) return false; Mesh m; if(!makeBipyramid(p,m)) return false; if(m.p.size()<in.p.size()&&finalizeCandidate(m)&&hausdorffOK(in.p,m.p,eps*.999)&&maskOK(om,m,.985)){out=move(m);return true;} return false;
}


static inline uint64_t faceKeySorted(int a,int b,int c){
    array<int,3> q={a,b,c}; sort(q.begin(),q.end()); return key3i(q[0],q[1],q[2]);
}
static unordered_set<uint64_t> makeFaceSet(const Mesh& m){
    unordered_set<uint64_t> S; S.reserve(m.f.size()*2+10);
    for(auto f:m.f) S.insert(faceKeySorted(f.a,f.b,f.c));
    return S;
}
static bool hasTri(const unordered_set<uint64_t>& S,int a,int b,int c){return S.find(faceKeySorted(a,b,c))!=S.end();}

static bool detectLatLonSphereGrid(const Mesh& in,int& R,int& Vv){
    int N=in.p.size(), M=in.f.size(); if(N<300||M!=2*(N-2)) return false;
    auto S=makeFaceSet(in);
    for(int v=8; v<=1024; ++v){
        if((N-2)%v) continue; int r=(N-2)/v; if(r<3) continue;
        int checks=0, ok=0; int step=max(1,v/80);
        auto id=[&](int rr,int j){j=(j%v+v)%v; return 2+(rr-1)*v+j;};
        for(int j=0;j<v;j+=step){
            checks++; ok += hasTri(S,0,id(1,j),id(1,j+1));
            checks++; ok += hasTri(S,1,id(r,j+1),id(r,j));
        }
        int rstep=max(1,r/50);
        for(int rr=1; rr<r; rr+=rstep){
            for(int j=0;j<v;j+=step){
                int a=id(rr,j), b=id(rr+1,j), c=id(rr+1,j+1), d=id(rr,j+1);
                bool ac=hasTri(S,a,b,c)&&hasTri(S,a,c,d);
                bool bd=hasTri(S,a,b,d)&&hasTri(S,b,d,c);
                checks+=2; ok += ac||bd; ok += ac||bd;
            }
        }
        if(checks>0 && ok*100 >= checks*97){R=r; Vv=v; return true;}
    }
    return false;
}
static bool makeLatLonFromOriginal(const Mesh& in,int R,int Vv,int R2,int V2,Mesh& out){
    if(R2<3||V2<8||2+R2*V2>= (int)in.p.size()) return false;
    Mesh m; m.p.reserve(2+R2*V2); m.f.reserve(2*R2*V2); m.p.push_back(in.p[0]); m.p.push_back(in.p[1]);
    vector<int> rr(R2);
    for(int i=0;i<R2;i++){ int r=(int)llround((double)(i+1)*(R+1)/(R2+1)); r=max(1,min(R,r)); rr[i]=r; }
    auto oid=[&](int r,int j){j=(j%Vv+Vv)%Vv; return 2+(r-1)*Vv+j;};
    auto nid=[&](int i,int j){j=(j%V2+V2)%V2; return 2+i*V2+j;};
    for(int i=0;i<R2;i++) for(int j=0;j<V2;j++){ int oj=(int)((long long)j*Vv/V2); m.p.push_back(in.p[oid(rr[i],oj)]); }
    for(int j=0;j<V2;j++) m.f.push_back({0,nid(0,j),nid(0,j+1)});
    for(int i=0;i+1<R2;i++) for(int j=0;j<V2;j++){int a=nid(i,j),b=nid(i+1,j),c=nid(i+1,j+1),d=nid(i,j+1); m.f.push_back({a,b,c}); m.f.push_back({a,c,d});}
    for(int j=0;j<V2;j++) m.f.push_back({1,nid(R2-1,j+1),nid(R2-1,j)});
    m.r.assign(m.p.size(),0); out=move(m); return true;
}
static bool tryLatLonGrid(const Mesh& in,Mesh& out,double eps,const Masks& om,double deadline){
    int R=0,Vv=0; if(!detectLatLonSphereGrid(in,R,Vv)) return false;
    struct T{int r,v; double th;}; vector<T> trials;
    auto add=[&](int r,int v,double th){r=max(3,min(R,r)); v=max(8,min(Vv,v)); if(2+r*v<(int)in.p.size()) trials.push_back({r,v,th});};
    add(R/10,Vv/10,.92); add(R/8,Vv/8,.925); add(R/6,Vv/6,.93); add(R/5,Vv/5,.935); add(R/4,Vv/4,.94); add(R/3,Vv/3,.945); add(R/2,Vv/2,.95);
    sort(trials.begin(),trials.end(),[](const T&a,const T&b){return a.r*a.v<b.r*b.v;});
    trials.erase(unique(trials.begin(),trials.end(),[](const T&a,const T&b){return a.r==b.r&&a.v==b.v;}),trials.end());
    Mesh best; bool got=false;
    for(auto t:trials){ if(nowsec()>deadline) break; Mesh cand; if(!makeLatLonFromOriginal(in,R,Vv,t.r,t.v,cand)) continue; if(!finalizeCandidate(cand)) continue; if(!hausdorffOK(in.p,cand.p,eps*.999)) continue; if(!maskOK(om,cand,t.th)) continue; if(!got||cand.p.size()<best.p.size()){best=move(cand); got=true;} }
    if(got){out=move(best); return true;} return false;
}

static bool detectPeriodicGrid(const Mesh& in,int& U,int& Vv){
    int N=in.p.size(), M=in.f.size(); if(N<300||M!=2*N) return false;
    auto S=makeFaceSet(in);
    vector<int> divisors;
    for(int v=6; v<=1024; ++v) if(N%v==0) divisors.push_back(v);
    sort(divisors.begin(),divisors.end(),[&](int a,int b){return abs(a-128)<abs(b-128);});
    for(int v:divisors){
        int u=N/v; if(u<3) continue;
        int stepU=max(1,u/70), stepV=max(1,v/70); int checks=0, ok=0;
        auto id=[&](int i,int j){i=(i%u+u)%u; j=(j%v+v)%v; return i*v+j;};
        for(int i=0;i<u;i+=stepU) for(int j=0;j<v;j+=stepV){
            int a=id(i,j), b=id(i+1,j), c=id(i+1,j+1), d=id(i,j+1);
            bool ac=hasTri(S,a,b,c)&&hasTri(S,a,c,d);
            bool bd=hasTri(S,a,b,d)&&hasTri(S,b,c,d);
            checks++; ok += (ac||bd);
        }
        if(checks>0 && ok*100 >= checks*98){U=u; Vv=v; return true;}
    }
    return false;
}
static bool makePeriodicFromOriginal(const Mesh& in,int U,int Vv,int U2,int V2,bool diagBD,Mesh& out){
    if(U2<3||V2<3||U2*V2>= (int)in.p.size()) return false;
    Mesh m; m.p.reserve(U2*V2); m.f.reserve(2*U2*V2);
    auto oid=[&](int i,int j){i=(i%U+U)%U; j=(j%Vv+Vv)%Vv; return i*Vv+j;};
    auto nid=[&](int i,int j){i=(i%U2+U2)%U2; j=(j%V2+V2)%V2; return i*V2+j;};
    for(int i=0;i<U2;i++) for(int j=0;j<V2;j++){ int oi=(int)((long long)i*U/U2), oj=(int)((long long)j*Vv/V2); m.p.push_back(in.p[oid(oi,oj)]); }
    for(int i=0;i<U2;i++) for(int j=0;j<V2;j++){
        int a=nid(i,j), b=nid(i+1,j), c=nid(i+1,j+1), d=nid(i,j+1);
        if(!diagBD){m.f.push_back({a,b,c}); m.f.push_back({a,c,d});}
        else {m.f.push_back({a,b,d}); m.f.push_back({b,c,d});}
    }
    m.r.assign(m.p.size(),0); out=move(m); return true;
}
static bool tryPeriodicGrid(const Mesh& in,Mesh& out,double eps,const Masks& om,double deadline){
    int U=0,Vv=0; if(!detectPeriodicGrid(in,U,Vv)) return false;
    struct T{int u,v; double th;}; vector<T> tr;
    auto add=[&](int u,int v,double th){u=max(3,min(U,u)); v=max(3,min(Vv,v)); if(u*v<(int)in.p.size()) tr.push_back({u,v,th});};
    add(U/6,Vv/6,.91); add(U/5,Vv/5,.92); add(U/4,Vv/4,.93); add(U/4,(Vv*7+19)/20,.935); add(U/3,Vv/3,.935); add(U/3,Vv/2,.94); add(U/2,Vv/2,.945); add(U/2,Vv,.95);
    sort(tr.begin(),tr.end(),[](const T&a,const T&b){return a.u*a.v<b.u*b.v;});
    tr.erase(unique(tr.begin(),tr.end(),[](const T&a,const T&b){return a.u==b.u&&a.v==b.v;}),tr.end());
    Mesh best; bool got=false;
    for(auto t:tr){ if(nowsec()>deadline) break; for(int bd=0;bd<2;bd++){
        Mesh cand; if(!makePeriodicFromOriginal(in,U,Vv,t.u,t.v,bd,cand)) continue; if(!finalizeCandidate(cand)) continue; if(!hausdorffOK(in.p,cand.p,eps*.999)) continue; if(!maskOK(om,cand,t.th)) continue; if(!got||cand.p.size()<best.p.size()){best=move(cand); got=true;}
    }}
    if(got){out=move(best); return true;} return false;
}

static P fromAxis(int ax,double t,double u,double v){ if(ax==0) return {t,u,v}; if(ax==1) return {u,t,v}; return {u,v,t}; }
static void toAxis(P p,int ax,double& t,double& u,double& v){ if(ax==0){t=p.x;u=p.y;v=p.z;} else if(ax==1){t=p.y;u=p.x;v=p.z;} else {t=p.z;u=p.x;v=p.y;} }
struct TorusFit{bool ok=false; int ax=2; double ct=0,cu=0,cv=0,R=0,r=0;};
static TorusFit fitAxisTorus(const Mesh& in,int ax,double eps){
    TorusFit fit; fit.ax=ax; double mnt=1e100,mxt=-1e100,mnu=1e100,mxu=-1e100,mnv=1e100,mxv=-1e100;
    for(auto p:in.p){double t,u,v; toAxis(p,ax,t,u,v); mnt=min(mnt,t);mxt=max(mxt,t);mnu=min(mnu,u);mxu=max(mxu,u);mnv=min(mnv,v);mxv=max(mxv,v);} 
    fit.ct=.5*(mnt+mxt); fit.cu=.5*(mnu+mxu); fit.cv=.5*(mnv+mxv);
    double minrho=1e100,maxrho=0; for(auto p:in.p){double t,u,v; toAxis(p,ax,t,u,v); double rho=hypot(u-fit.cu,v-fit.cv); minrho=min(minrho,rho); maxrho=max(maxrho,rho);} 
    if(!(maxrho>minrho&&minrho>1e-12)) return fit; fit.R=.5*(maxrho+minrho); double rrho=.5*(maxrho-minrho), rt=.5*(mxt-mnt); fit.r=.5*(rrho+rt); if(!(fit.R>fit.r*1.35&&fit.r>eps*.2)) return fit;
    double mx=0,ss=0; int cnt=0; for(auto p:in.p){double t,u,v; toAxis(p,ax,t,u,v); double rho=hypot(u-fit.cu,v-fit.cv); double d=fabs(hypot(rho-fit.R,t-fit.ct)-fit.r); mx=max(mx,d); ss+=d*d; cnt++;}
    double rms=sqrt(ss/max(1,cnt)); if(mx<eps*.86&&rms<eps*.32) fit.ok=true; return fit;
}
static bool makeTorus(const TorusFit& fit,int A,int B,Mesh& out){
    if(!fit.ok||A<12||B<6) return false; Mesh m; m.p.reserve(A*B); m.f.reserve(2*A*B); const double pi=acos(-1.0);
    for(int i=0;i<A;i++){double th=2*pi*i/A, co=cos(th), si=sin(th); for(int j=0;j<B;j++){double ph=2*pi*j/B, cp=cos(ph), sp=sin(ph); double rho=fit.R+fit.r*cp; m.p.push_back(fromAxis(fit.ax,fit.ct+fit.r*sp,fit.cu+rho*co,fit.cv+rho*si));}}
    auto id=[&](int i,int j){i=(i%A+A)%A; j=(j%B+B)%B; return i*B+j;};
    for(int i=0;i<A;i++)for(int j=0;j<B;j++){int a=id(i,j),b=id(i+1,j),c=id(i+1,j+1),d=id(i,j+1); m.f.push_back({a,b,c}); m.f.push_back({a,c,d});}
    m.r.assign(m.p.size(),0); out=move(m); return true;
}
static bool tryAnalyticTorus(const Mesh& in,Mesh& out,double eps,const Masks& om,double deadline){
    if(in.p.size()<600) return false; TorusFit best; for(int ax=0;ax<3;ax++){auto f=fitAxisTorus(in,ax,eps); if(f.ok&&(!best.ok||f.r<best.r)) best=f;} if(!best.ok) return false;
    struct T{int a,b; double th;}; vector<T> tr={{24,10,.92},{32,12,.925},{40,14,.93},{48,16,.935},{64,20,.94},{80,24,.945}};
    Mesh bm; bool got=false; for(auto t:tr){if(nowsec()>deadline)break; Mesh cand; if(!makeTorus(best,t.a,t.b,cand))continue; if(cand.p.size()>=in.p.size())continue; if(!finalizeCandidate(cand))continue; if(!hausdorffOK(in.p,cand.p,eps*.995))continue; if(!maskOK(om,cand,t.th))continue; if(!got||cand.p.size()<bm.p.size()){bm=move(cand);got=true;}}
    if(got){out=move(bm);return true;} return false;
}


static bool makeBoxGrid(const Mesh& in,double eps,Mesh& out){
    if(in.p.size()<20) return false; P mn=in.p[0],mx=in.p[0];
    for(auto p:in.p){mn.x=min(mn.x,p.x);mn.y=min(mn.y,p.y);mn.z=min(mn.z,p.z);mx.x=max(mx.x,p.x);mx.y=max(mx.y,p.y);mx.z=max(mx.z,p.z);} 
    double ex=mx.x-mn.x,ey=mx.y-mn.y,ez=mx.z-mn.z; if(ex<eps*.2||ey<eps*.2||ez<eps*.2) return false;
    int nearCnt=0; double tol=eps*.75;
    for(auto p:in.p){double d=min({fabs(p.x-mn.x),fabs(p.x-mx.x),fabs(p.y-mn.y),fabs(p.y-mx.y),fabs(p.z-mn.z),fabs(p.z-mx.z)}); if(d<=tol) nearCnt++;}
    if(nearCnt*100 < (int)in.p.size()*98) return false;
    double h=eps*.82; int nx=max(1,(int)ceil(ex/h)), ny=max(1,(int)ceil(ey/h)), nz=max(1,(int)ceil(ez/h));
    long long vc=(long long)(nx+1)*(ny+1)*(nz+1) - (long long)max(0,nx-1)*max(0,ny-1)*max(0,nz-1);
    if(vc>= (long long)in.p.size() || vc>250000) return false;
    vector<double> xs(nx+1),ys(ny+1),zs(nz+1); for(int i=0;i<=nx;i++)xs[i]=mn.x+ex*i/nx; for(int j=0;j<=ny;j++)ys[j]=mn.y+ey*j/ny; for(int k=0;k<=nz;k++)zs[k]=mn.z+ez*k/nz;
    vector<int> id((nx+1)*(ny+1)*(nz+1),-1); auto idx=[&](int i,int j,int k){return (i*(ny+1)+j)*(nz+1)+k;};
    Mesh m; m.p.reserve((size_t)vc); m.f.reserve((size_t)2*(nx*ny+nx*nz+ny*nz));
    auto get=[&](int i,int j,int k){int &r=id[idx(i,j,k)]; if(r<0){r=m.p.size(); m.p.push_back({xs[i],ys[j],zs[k]});} return r;};
    auto addq=[&](int a,int b,int c,int d){m.f.push_back({a,b,c}); m.f.push_back({a,c,d});};
    for(int i=0;i<nx;i++)for(int j=0;j<ny;j++){addq(get(i,j,nz),get(i+1,j,nz),get(i+1,j+1,nz),get(i,j+1,nz)); addq(get(i,j,0),get(i,j+1,0),get(i+1,j+1,0),get(i+1,j,0));}
    for(int i=0;i<nx;i++)for(int k=0;k<nz;k++){addq(get(i,ny,k),get(i+1,ny,k),get(i+1,ny,k+1),get(i,ny,k+1)); addq(get(i,0,k),get(i,0,k+1),get(i+1,0,k+1),get(i+1,0,k));}
    for(int j=0;j<ny;j++)for(int k=0;k<nz;k++){addq(get(nx,j,k),get(nx,j+1,k),get(nx,j+1,k+1),get(nx,j,k+1)); addq(get(0,j,k),get(0,j,k+1),get(0,j+1,k+1),get(0,j+1,k));}
    m.r.assign(m.p.size(),0); out=move(m); return true;
}
static bool tryBoxGrid(const Mesh& in,Mesh& out,double eps,const Masks& om){
    Mesh cand; if(!makeBoxGrid(in,eps,cand)) return false; if(!finalizeCandidate(cand)) return false; if(!hausdorffOK(in.p,cand.p,eps*.995)) return false; if(!maskOK(om,cand,.955)) return false; out=move(cand); return true;
}


static void jacobiEigen3(double a[3][3],P e[3]){
    double v[3][3]={{1,0,0},{0,1,0},{0,0,1}};
    for(int it=0;it<60;it++){
        int p=0,q=1; double best=fabs(a[0][1]); if(fabs(a[0][2])>best){p=0;q=2;best=fabs(a[0][2]);} if(fabs(a[1][2])>best){p=1;q=2;best=fabs(a[1][2]);}
        if(best<1e-18) break; double app=a[p][p], aqq=a[q][q], apq=a[p][q]; double tau=(aqq-app)/(2*apq); double tt=(tau>=0?1:-1)/(fabs(tau)+sqrt(1+tau*tau)); double cs=1/sqrt(1+tt*tt), sn=tt*cs;
        for(int k=0;k<3;k++) if(k!=p&&k!=q){double akp=a[k][p],akq=a[k][q]; a[k][p]=a[p][k]=cs*akp-sn*akq; a[k][q]=a[q][k]=sn*akp+cs*akq;}
        a[p][p]=cs*cs*app-2*sn*cs*apq+sn*sn*aqq; a[q][q]=sn*sn*app+2*sn*cs*apq+cs*cs*aqq; a[p][q]=a[q][p]=0;
        for(int k=0;k<3;k++){double vkp=v[k][p],vkq=v[k][q]; v[k][p]=cs*vkp-sn*vkq; v[k][q]=sn*vkp+cs*vkq;}
    }
    int ord[3]={0,1,2}; sort(ord,ord+3,[&](int i,int j){return a[i][i]>a[j][j];});
    for(int j=0;j<3;j++){int c=ord[j]; e[j]=unit({v[0][c],v[1][c],v[2][c]});}
    if(dotp(crossp(e[0],e[1]),e[2])<0) e[2]=e[2]*-1;
}
struct EllFit{bool ok=false; P c; P ax[3]; double r[3]; double rms=1e100,mx=1e100;};
static EllFit fitEllipsoidPCA(const Mesh& in,double eps){
    EllFit fit; if(in.p.size()<100) return fit; P mean{0,0,0}; for(auto p:in.p) mean=mean+p; mean=mean/(double)in.p.size(); double C[3][3]={{0}};
    for(auto p:in.p){P q=p-mean; double x[3]={q.x,q.y,q.z}; for(int i=0;i<3;i++)for(int j=0;j<3;j++) C[i][j]+=x[i]*x[j];}
    for(int i=0;i<3;i++)for(int j=0;j<3;j++) C[i][j]/=max(1,(int)in.p.size()); jacobiEigen3(C,fit.ax);
    double lo[3]={1e100,1e100,1e100}, hi[3]={-1e100,-1e100,-1e100};
    for(auto p:in.p){P q=p-mean; for(int k=0;k<3;k++){double t=dotp(q,fit.ax[k]); lo[k]=min(lo[k],t); hi[k]=max(hi[k],t);}}
    fit.c=mean; for(int k=0;k<3;k++){double mid=.5*(lo[k]+hi[k]); fit.c=fit.c+fit.ax[k]*mid; fit.r[k]=.5*(hi[k]-lo[k]); if(!(fit.r[k]>eps*.3)) return fit;}
    double ss=0,mx=0; int cnt=0; double scale=(fit.r[0]+fit.r[1]+fit.r[2])/3.0;
    for(auto p:in.p){P q=p-fit.c; double s2=0; for(int k=0;k<3;k++){double t=dotp(q,fit.ax[k])/fit.r[k]; s2+=t*t;} double s=sqrt(max(1e-30,s2)); double d=fabs(1.0-1.0/s)*norm(q); mx=max(mx,d); ss+=d*d; cnt++;}
    fit.rms=sqrt(ss/max(1,cnt)); fit.mx=mx; if(mx<eps*.82 && fit.rms<eps*.30) fit.ok=true; return fit;
}
static bool makeEllipsoid(const EllFit& fit,int lat,int lon,Mesh& out){
    if(!fit.ok||lat<4||lon<8) return false; Mesh m; m.p.reserve(2+(lat-1)*lon); m.f.reserve(2*lat*lon); const double pi=acos(-1.0);
    auto pt=[&](double x,double y,double z){return fit.c + fit.ax[0]*(fit.r[0]*x) + fit.ax[1]*(fit.r[1]*y) + fit.ax[2]*(fit.r[2]*z);};
    m.p.push_back(pt(0,0,1)); m.p.push_back(pt(0,0,-1)); auto id=[&](int r,int j){j=(j%lon+lon)%lon; return 2+(r-1)*lon+j;};
    for(int r=1;r<=lat-1;r++){double th=pi*r/lat, st=sin(th), ct=cos(th); for(int j=0;j<lon;j++){double ph=2*pi*j/lon; m.p.push_back(pt(st*cos(ph),st*sin(ph),ct));}}
    for(int j=0;j<lon;j++) m.f.push_back({0,id(1,j),id(1,j+1)});
    for(int r=1;r<lat-1;r++) for(int j=0;j<lon;j++){int a=id(r,j),b=id(r+1,j),c=id(r+1,j+1),d=id(r,j+1); m.f.push_back({a,b,c}); m.f.push_back({a,c,d});}
    for(int j=0;j<lon;j++) m.f.push_back({1,id(lat-1,j+1),id(lat-1,j)}); m.r.assign(m.p.size(),0); out=move(m); return true;
}
static bool tryEllipsoid(const Mesh& in,Mesh& out,double eps,const Masks& om,double deadline){
    EllFit fit=fitEllipsoidPCA(in,eps); if(!fit.ok) return false; struct T{int lat,lon; double th;}; vector<T> tr={{12,24,.91},{16,32,.92},{20,40,.93},{24,48,.94},{28,56,.945},{32,64,.95}};
    Mesh best; bool got=false; for(auto t:tr){if(nowsec()>deadline)break; Mesh cand; if(!makeEllipsoid(fit,t.lat,t.lon,cand))continue; if(cand.p.size()>=in.p.size())continue; if(!finalizeCandidate(cand))continue; if(!hausdorffOK(in.p,cand.p,eps*.995))continue; if(!maskOK(om,cand,t.th))continue; if(!got||cand.p.size()<best.p.size()){best=move(cand);got=true;}}
    if(got){out=move(best);return true;} return false;
}


struct FrustumFit{bool ok=false; P c,w,u,v; double t0=0,t1=0,cu=0,cv=0,rx0=0,rx1=0,ry0=0,ry1=0,mx=1e100,rms=1e100;};
static void basisFromAxis(P w,P hint,P& u,P& v){
    w=unit(w); u=hint-w*dotp(hint,w); if(n2(u)<1e-20){P h=fabs(w.x)<.75?P{1,0,0}:(fabs(w.y)<.75?P{0,1,0}:P{0,0,1}); u=h-w*dotp(h,w);} u=unit(u); v=unit(crossp(w,u));
}
static FrustumFit fitFrustumBasis(const Mesh& in,P w,P hint,double eps){
    FrustumFit f; basisFromAxis(w,hint,f.u,f.v); f.w=unit(w); if(n2(f.w)<.5) return f;
    double lot=1e100,hit=-1e100,lou=1e100,hiu=-1e100,lov=1e100,hiv=-1e100; vector<array<double,3>> pr; pr.reserve(in.p.size());
    for(auto p:in.p){double t=dotp(p,f.w), x=dotp(p,f.u), y=dotp(p,f.v); pr.push_back({t,x,y}); lot=min(lot,t);hit=max(hit,t);lou=min(lou,x);hiu=max(hiu,x);lov=min(lov,y);hiv=max(hiv,y);} 
    if(!(hit>lot+eps)) return f; f.t0=lot; f.t1=hit; f.cu=.5*(lou+hiu); f.cv=.5*(lov+hiv); double len=hit-lot; const int B=36; vector<double> sx(B,0),sy(B,0); vector<int> bc(B,0);
    for(auto q:pr){double s=(q[0]-lot)/len; int b=min(B-1,max(0,(int)floor(s*B))); sx[b]=max(sx[b],fabs(q[1]-f.cu)); sy[b]=max(sy[b],fabs(q[2]-f.cv)); bc[b]++;}
    double S=0,St=0,Stt=0,Sx=0,Stx=0,Sy=0,Sty=0; int bins=0; for(int b=0;b<B;b++) if(bc[b]>0&&sx[b]>eps*.2&&sy[b]>eps*.2){double s=(b+.5)/B; S++;St+=s;Stt+=s*s;Sx+=sx[b];Stx+=s*sx[b];Sy+=sy[b];Sty+=s*sy[b];bins++;}
    if(bins<5) return f; double den=S*Stt-St*St; if(fabs(den)<1e-18) return f; double bx=(S*Stx-St*Sx)/den, ax=(Sx-bx*St)/S; double by=(S*Sty-St*Sy)/den, ay=(Sy-by*St)/S;
    f.rx0=ax; f.rx1=ax+bx; f.ry0=ay; f.ry1=ay+by; double mr=max({f.rx0,f.rx1,f.ry0,f.ry1}), nr=min({f.rx0,f.rx1,f.ry0,f.ry1}); if(!(mr>eps*.8&&nr>eps*.15)) return f; if(mr/nr>20) return f;
    double ss=0,mx=0; int cnt=0,bad=0; for(auto q:pr){double s=(q[0]-lot)/len; s=min(1.0,max(0.0,s)); double rx=f.rx0+(f.rx1-f.rx0)*s, ry=f.ry0+(f.ry1-f.ry0)*s; double X=fabs(q[1]-f.cu), Y=fabs(q[2]-f.cv); double e=max(X/max(rx,1e-30),Y/max(ry,1e-30)); double side=fabs(e-1.0)*min(rx,ry); double cap=1e100; double dt=min(fabs(q[0]-lot),fabs(q[0]-hit)); if(e<=1.04) cap=dt; else cap=hypot(dt,(e-1.0)*min(rx,ry)); double d=min(side,cap); mx=max(mx,d); ss+=d*d; cnt++; if(d>eps*.95) bad++; }
    f.mx=mx; f.rms=sqrt(ss/max(1,cnt)); if(bad*100<=cnt*2 && mx<eps*1.10 && f.rms<eps*.36) f.ok=true; return f;
}
static bool makeFrustum(const FrustumFit& fit,int seg,int sides,Mesh& out){
    if(!fit.ok||seg<1||sides<6) return false; Mesh m; const double pi=acos(-1.0); m.p.reserve((seg+1)*sides+2); m.f.reserve(2*seg*sides+2*sides);
    auto pt=[&](double t,double x,double y){return fit.w*t + fit.u*x + fit.v*y;};
    for(int i=0;i<=seg;i++){double s=(double)i/seg, t=fit.t0+(fit.t1-fit.t0)*s, rx=fit.rx0+(fit.rx1-fit.rx0)*s, ry=fit.ry0+(fit.ry1-fit.ry0)*s; for(int j=0;j<sides;j++){double a=2*pi*j/sides; m.p.push_back(pt(t,fit.cu+rx*cos(a),fit.cv+ry*sin(a)));}}
    int c0=m.p.size(); m.p.push_back(pt(fit.t0,fit.cu,fit.cv)); int c1=m.p.size(); m.p.push_back(pt(fit.t1,fit.cu,fit.cv)); auto id=[&](int i,int j){j=(j%sides+sides)%sides; return i*sides+j;};
    for(int i=0;i<seg;i++) for(int j=0;j<sides;j++){int a=id(i,j),b=id(i+1,j),c=id(i+1,j+1),d=id(i,j+1); m.f.push_back({a,b,c}); m.f.push_back({a,c,d});}
    for(int j=0;j<sides;j++){m.f.push_back({c0,id(0,j+1),id(0,j)}); m.f.push_back({c1,id(seg,j),id(seg,j+1)});} m.r.assign(m.p.size(),0); out=move(m); return true;
}
static bool tryFrustum(const Mesh& in,Mesh& out,double eps,const Masks& om,double deadline){
    if(in.p.size()<300) return false; P mean{0,0,0}; for(auto p:in.p) mean=mean+p; mean=mean/(double)in.p.size(); double C[3][3]={{0}}; for(auto p:in.p){P q=p-mean; double x[3]={q.x,q.y,q.z}; for(int i=0;i<3;i++)for(int j=0;j<3;j++)C[i][j]+=x[i]*x[j];} for(int i=0;i<3;i++)for(int j=0;j<3;j++)C[i][j]/=max(1,(int)in.p.size()); P ax[3]; jacobiEigen3(C,ax);
    vector<pair<P,P>> bases={{ax[0],ax[1]},{ax[1],ax[0]},{ax[2],ax[0]},{P{1,0,0},P{0,1,0}},{P{0,1,0},P{1,0,0}},{P{0,0,1},P{1,0,0}}};
    FrustumFit best; for(auto &bh:bases){auto f=fitFrustumBasis(in,bh.first,bh.second,eps); if(f.ok&&(!best.ok||f.rms<best.rms)) best=f;} if(!best.ok) return false;
    struct T{int seg,sides; double th;}; vector<T> tr={{2,12,.90},{3,16,.91},{4,20,.92},{6,24,.93},{8,32,.94},{10,40,.945},{14,48,.95}};
    Mesh bm; bool got=false; for(auto t:tr){if(nowsec()>deadline)break; Mesh cand; if(!makeFrustum(best,t.seg,t.sides,cand))continue; if(cand.p.size()>=in.p.size())continue; if(!finalizeCandidate(cand))continue; if(!hausdorffOK(in.p,cand.p,eps*.995))continue; if(!maskOK(om,cand,t.th))continue; if(!got||cand.p.size()<bm.p.size()){bm=move(cand);got=true;}}
    if(got){out=move(bm);return true;} return false;
}

struct Q{double a00=0,a01=0,a02=0,a03=0,a11=0,a12=0,a13=0,a22=0,a23=0,a33=0;};
static inline Q operator+(const Q&x,const Q&y){return {x.a00+y.a00,x.a01+y.a01,x.a02+y.a02,x.a03+y.a03,x.a11+y.a11,x.a12+y.a12,x.a13+y.a13,x.a22+y.a22,x.a23+y.a23,x.a33+y.a33};}
static inline void addPlane(Q& q,P n,double d,double w){double a=n.x,b=n.y,c=n.z; q.a00+=w*a*a; q.a01+=w*a*b; q.a02+=w*a*c; q.a03+=w*a*d; q.a11+=w*b*b; q.a12+=w*b*c; q.a13+=w*b*d; q.a22+=w*c*c; q.a23+=w*c*d; q.a33+=w*d*d;}
static inline double qcost(const Q&q,P p){double x=p.x,y=p.y,z=p.z; return q.a00*x*x+2*q.a01*x*y+2*q.a02*x*z+2*q.a03*x+q.a11*y*y+2*q.a12*y*z+2*q.a13*y+q.a22*z*z+2*q.a23*z+q.a33;}
static bool solveQ(const Q&q,P& p){
    double A00=q.a00,A01=q.a01,A02=q.a02,A11=q.a11,A12=q.a12,A22=q.a22; double b0=-q.a03,b1=-q.a13,b2=-q.a23;
    double det=A00*(A11*A22-A12*A12)-A01*(A01*A22-A12*A02)+A02*(A01*A12-A11*A02); if(fabs(det)<1e-16) return false;
    double dx=b0*(A11*A22-A12*A12)-A01*(b1*A22-A12*b2)+A02*(b1*A12-A11*b2);
    double dy=A00*(b1*A22-A12*b2)-b0*(A01*A22-A12*A02)+A02*(A01*b2-b1*A02);
    double dz=A00*(A11*b2-b1*A12)-A01*(A01*b2-b1*A02)+b0*(A01*A12-A11*A02);
    p={dx/det,dy/det,dz/det}; return isfinite(p.x)&&isfinite(p.y)&&isfinite(p.z);
}
struct PassData{vector<EdgeInfo> edges; vector<int> np,adj,fp,vf; vector<P> fn; vector<double> area,sal; vector<Q> quad;};
static bool prepare(const Mesh& m,PassData& d){
    if(!buildEdges(m,&d.edges)) return false; int n=m.p.size(); d.np.assign(n+1,0); for(auto&e:d.edges){d.np[e.a+1]++;d.np[e.b+1]++;} for(int i=1;i<=n;i++)d.np[i]+=d.np[i-1]; d.adj.assign(d.np[n],0); vector<int> cur=d.np; for(auto&e:d.edges){d.adj[cur[e.a]++]=e.b;d.adj[cur[e.b]++]=e.a;}
    d.fp.assign(n+1,0); for(auto f:m.f){d.fp[f.a+1]++;d.fp[f.b+1]++;d.fp[f.c+1]++;} for(int i=1;i<=n;i++)d.fp[i]+=d.fp[i-1]; d.vf.assign(d.fp[n],0); cur=d.fp;
    d.fn.resize(m.f.size()); d.area.resize(m.f.size()); d.quad.assign(n,Q());
    for(int i=0;i<(int)m.f.size();i++){auto t=m.f[i]; P cr=crossp(m.p[t.b]-m.p[t.a],m.p[t.c]-m.p[t.a]); double ar=norm(cr); if(ar<=0)return false; P nn=cr/ar; d.fn[i]=nn; d.area[i]=ar*.5; double dd=-dotp(nn,m.p[t.a]); double w=max(1e-18,d.area[i]); addPlane(d.quad[t.a],nn,dd,w); addPlane(d.quad[t.b],nn,dd,w); addPlane(d.quad[t.c],nn,dd,w); d.vf[cur[t.a]++]=i;d.vf[cur[t.b]++]=i;d.vf[cur[t.c]++]=i;}
    d.sal.assign(n,0); vector<P> av(n,{0,0,0}); vector<double> amin(n,1),amax(n,-1); P ax[3]={{1,0,0},{0,1,0},{0,0,1}};
    for(int i=0;i<(int)m.f.size();i++){auto t=m.f[i]; P nn=d.fn[i]; av[t.a]=av[t.a]+nn; av[t.b]=av[t.b]+nn; av[t.c]=av[t.c]+nn;}
    for(int v=0;v<n;v++){P u=unit(av[v]); double mx=0; for(int k=d.fp[v];k<d.fp[v+1];k++){P nn=d.fn[d.vf[k]]; mx=max(mx,1.0-dotp(u,nn));} d.sal[v]=min(1.0,mx/0.25);}    
    return true;
}
static bool linkOK(const PassData& d,int u,int v,int o1,int o2,vector<int>& mark,int token){
    for(int k=d.np[u];k<d.np[u+1];k++) mark[d.adj[k]]=token; int cnt=0; bool h1=false,h2=false; for(int k=d.np[v];k<d.np[v+1];k++) if(mark[d.adj[k]]==token){cnt++; h1|=d.adj[k]==o1; h2|=d.adj[k]==o2;} return cnt==2&&h1&&h2;
}
struct Cand{float cost; int drop,keep; P pos; double nr;};
static bool checkCandPos(const Mesh& m,const PassData& d,int drop,int keep,int f1,int f2,P x,double eps,double aggr,Cand& out){
    double s=max(d.sal[drop],d.sal[keep]); double lim=eps*(0.985-0.50*pow(s,0.7)); lim=max(lim,eps*(0.20+0.18*aggr)); lim=min(lim,eps*.995);
    double nr=max(m.r[drop]+norm(x-m.p[drop]),m.r[keep]+norm(x-m.p[keep])); if(nr>lim) return false;
    double mindot=-0.08+0.48*s; double pc=0; int cnt=0; Q qq=d.quad[drop]+d.quad[keep];
    auto posOf=[&](int id)->P{return (id==drop||id==keep)?x:m.p[id];};
    auto chkFace=[&](int fi)->bool{F ot=m.f[fi]; F nt=ot; if(nt.a==drop)nt.a=keep; if(nt.b==drop)nt.b=keep; if(nt.c==drop)nt.c=keep; if(nt.a==nt.b||nt.b==nt.c||nt.c==nt.a)return false; P A=posOf(nt.a),B=posOf(nt.b),C=posOf(nt.c); P cr=crossp(B-A,C-A); double ar=norm(cr); if(ar<=1e-20)return false; P nn=cr/ar; if(dotp(nn,d.fn[fi])<mindot)return false; pc+=qcost(qq,x)*(0.05+d.area[fi]); cnt++; return true;};
    for(int k=d.fp[drop];k<d.fp[drop+1];k++){int fi=d.vf[k]; if(fi==f1||fi==f2)continue; if(!chkFace(fi))return false;}
    for(int k=d.fp[keep];k<d.fp[keep+1];k++){int fi=d.vf[k]; if(fi==f1||fi==f2)continue; if(!chkFace(fi))return false;}
    double len2=n2(m.p[drop]-m.p[keep]); double cost=qcost(qq,x)+len2*(0.002+1.4*s)+nr*nr*(0.01+0.4*s)+pc*1e-4; out={(float)cost,drop,keep,x,nr}; return true;
}
static bool evalEdge(const Mesh& m,const PassData& d,const EdgeInfo& e,double eps,double aggr,Cand& best){
    int drop=e.a,keep=e.b; Q q=d.quad[drop]+d.quad[keep]; vector<P> cand; P opt; if(solveQ(q,opt)) cand.push_back(opt); cand.push_back((m.p[drop]+m.p[keep])*.5); cand.push_back(m.p[drop]); cand.push_back(m.p[keep]); cand.push_back((m.p[drop]*.35+m.p[keep]*.65)); cand.push_back((m.p[drop]*.65+m.p[keep]*.35));
    bool ok=false; Cand bc; bc.cost=1e30; for(P x:cand){Cand c; if(checkCandPos(m,d,drop,keep,e.f1,e.f2,x,eps,aggr,c)&&c.cost<bc.cost){bc=c;ok=true;}} if(!ok){swap(drop,keep); cand.clear(); if(solveQ(q,opt)) cand.push_back(opt); cand.push_back((m.p[drop]+m.p[keep])*.5); cand.push_back(m.p[drop]); cand.push_back(m.p[keep]); for(P x:cand){Cand c; if(checkCandPos(m,d,drop,keep,e.f1,e.f2,x,eps,aggr,c)&&c.cost<bc.cost){bc=c;ok=true;}}}
    if(ok)best=bc; return ok;
}
static bool qemPass(Mesh& m,double eps,int passNo,double deadline){
    if(nowsec()>deadline) return false; PassData d; if(!prepare(m,d)) return false; int n=m.p.size();
    P mn=m.p[0],mx=m.p[0]; for(auto p:m.p){mn.x=min(mn.x,p.x);mn.y=min(mn.y,p.y);mn.z=min(mn.z,p.z);mx.x=max(mx.x,p.x);mx.y=max(mx.y,p.y);mx.z=max(mx.z,p.z);} 
    for(int i=0;i<n;i++){int ec=0; if(m.p[i].x<mn.x+eps*.22||m.p[i].x>mx.x-eps*.22)ec++; if(m.p[i].y<mn.y+eps*.22||m.p[i].y>mx.y-eps*.22)ec++; if(m.p[i].z<mn.z+eps*.22||m.p[i].z>mx.z-eps*.22)ec++; if(ec>=2)d.sal[i]=max(d.sal[i],.75); else if(ec)d.sal[i]=max(d.sal[i],.18);}    
    vector<Cand> cands; cands.reserve(d.edges.size()/3+10); vector<int> mark(n,0); int tok=1; double aggr=min(1.0,0.42+0.035*passNo);
    for(auto &e:d.edges){if(nowsec()>deadline)break; if(tok>INT_MAX-10){fill(mark.begin(),mark.end(),0);tok=1;} if(!linkOK(d,e.a,e.b,e.o1,e.o2,mark,tok++))continue; Cand c; if(evalEdge(m,d,e,eps,aggr,c)) cands.push_back(c);} if(cands.empty())return false;
    sort(cands.begin(),cands.end(),[](const Cand&a,const Cand&b){return a.cost<b.cost;}); int minSel=max(4,n/260000); int cap=INT_MAX;
    for(int attempt=0;attempt<6;attempt++){
        vector<unsigned char> busy(n,0),dead(n,0); vector<int> rep(n); iota(rep.begin(),rep.end(),0); vector<P> np=m.p; vector<double> nr=m.r; int sel=0;
        for(auto &c:cands){if(sel>=cap)break; if(dead[c.drop]||dead[c.keep]||busy[c.drop]||busy[c.keep])continue; rep[c.drop]=c.keep; dead[c.drop]=1; np[c.keep]=c.pos; nr[c.keep]=c.nr; sel++; auto block=[&](int v){busy[v]=1; for(int k=d.np[v];k<d.np[v+1];k++)busy[d.adj[k]]=1;}; block(c.drop); block(c.keep);} if(sel<minSel)return false;
        Mesh nm; nm.p=move(np); nm.r=move(nr); nm.f.reserve(m.f.size()>2ull*sel?m.f.size()-2ull*sel:m.f.size());
        for(auto t:m.f){t.a=rep[t.a];t.b=rep[t.b];t.c=rep[t.c]; if(t.a==t.b||t.b==t.c||t.c==t.a)continue; nm.f.push_back(t);} compactMesh(nm); bool rr=true; for(double x:nm.r) if(x>eps*1.0001){rr=false;break;} if(rr&&finalizeCandidate(nm)){m=move(nm);return true;} cap=sel/2;
    }
    return false;
}
static Mesh simplifyQEM(Mesh in,double eps,const Masks& om,double deadline){
    Mesh best=in,cur=in; if(cur.r.size()!=cur.p.size())cur.r.assign(cur.p.size(),0); double th=.895; int stale=0;
    for(int pass=0;pass<90 && nowsec()<deadline;pass++){
        int before=cur.p.size(); Mesh old=cur; if(!qemPass(cur,eps,pass,deadline)){cur=move(old);break;} if(cur.p.size()>=before*.999) stale++; else stale=0;
        if((pass&1)==1||cur.p.size()<best.p.size()*0.94){ if(maskOK(om,cur,th)){ if(cur.p.size()<best.p.size())best=cur; } else {cur=move(old); break;} }
        if(stale>=3)break;
    }
    if(cur.p.size()<best.p.size()&&maskOK(om,cur,th))best=move(cur); return best;
}

int main(){
    ios::sync_with_stdio(false); cin.tie(nullptr);
    int V,FN; if(!(cin>>V>>FN)) return 0; Mesh in; in.p.resize(V); in.f.resize(FN); in.r.assign(V,0); string s;
    for(int i=0;i<V;i++){cin>>s>>in.p[i].x>>in.p[i].y>>in.p[i].z;}
    for(int i=0;i<FN;i++){cin>>s>>in.f[i].a>>in.f[i].b>>in.f[i].c; --in.f[i].a;--in.f[i].b;--in.f[i].c;}
    P mn=in.p[0],mx=in.p[0]; for(auto p:in.p){mn.x=min(mn.x,p.x);mn.y=min(mn.y,p.y);mn.z=min(mn.z,p.z);mx.x=max(mx.x,p.x);mx.y=max(mx.y,p.y);mx.z=max(mx.z,p.z);} double diag=max(1e-12,norm(mx-mn)); double eps=diag*0.05*0.999999;
    int R = (V>400000?80:96); double start=nowsec(); double deadline=start+18.65; Masks om=rasterMesh(in,mn,mx,R);
    Mesh best=in,cand;
    if(validateMesh(in)){
        if(nowsec()<deadline&&trySphere(in,cand,eps,om)&&cand.p.size()<best.p.size()) best=move(cand);
        if(nowsec()<deadline&&tryEllipsoid(in,cand,eps,om,deadline-1.0)&&cand.p.size()<best.p.size()) best=move(cand);
        if(nowsec()<deadline&&tryFrustum(in,cand,eps,om,deadline-1.0)&&cand.p.size()<best.p.size()) best=move(cand);
        if(nowsec()<deadline&&tryBoxGrid(in,cand,eps,om)&&cand.p.size()<best.p.size()) best=move(cand);
        if(nowsec()<deadline-1.0&&tryLatLonGrid(in,cand,eps,om,deadline-1.0)&&cand.p.size()<best.p.size()) best=move(cand);
        if(nowsec()<deadline-1.0&&tryPeriodicGrid(in,cand,eps,om,deadline-1.0)&&cand.p.size()<best.p.size()) best=move(cand);
        if(nowsec()<deadline-1.0&&tryAnalyticTorus(in,cand,eps,om,deadline-1.0)&&cand.p.size()<best.p.size()) best=move(cand);
        if(nowsec()<deadline-2.0 && V>80){
            CoverResult cr=greedyCover(in.p,eps,deadline-5.5);
            if(nowsec()<deadline-4.8 && mappedCoverMesh(in,cr,cand,om,eps)&&cand.p.size()<best.p.size()) best=move(cand);
            if(nowsec()<deadline-4.2 && coverBipyramid(in,cr,cand,om,eps)&&cand.p.size()<best.p.size()) best=move(cand);
        }
        if(nowsec()<deadline){ Mesh q=simplifyQEM(in,eps,om,deadline); if(q.p.size()<best.p.size()) best=move(q); }
    }
    cout.setf(ios::fixed); cout<<setprecision(10);
    cout<<best.p.size()<<' '<<best.f.size()<<'\n';
    for(auto p:best.p) cout<<"v "<<p.x<<' '<<p.y<<' '<<p.z<<'\n';
    for(auto f:best.f) cout<<"f "<<f.a+1<<' '<<f.b+1<<' '<<f.c+1<<'\n';
    return 0;
}
