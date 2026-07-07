#include <bits/stdc++.h>
using namespace std;

struct P{double x,y,z;};
struct F{int a,b,c;};
static inline P operator+(const P&a,const P&b){return {a.x+b.x,a.y+b.y,a.z+b.z};}
static inline P operator-(const P&a,const P&b){return {a.x-b.x,a.y-b.y,a.z-b.z};}
static inline P operator*(const P&a,double s){return {a.x*s,a.y*s,a.z*s};}
static inline double dotp(const P&a,const P&b){return a.x*b.x+a.y*b.y+a.z*b.z;}
static inline P crossp(const P&a,const P&b){return {a.y*b.z-a.z*b.y,a.z*b.x-a.x*b.z,a.x*b.y-a.y*b.x};}
static inline double d2(const P&a,const P&b){double x=a.x-b.x,y=a.y-b.y,z=a.z-b.z; return x*x+y*y+z*z;}
static inline double norm2(const P&a){return dotp(a,a);} 

struct Mesh{ vector<P> v; vector<F> f; };

struct DSU{ vector<int> p,sz; DSU(int n=0){init(n);} void init(int n){p.resize(n); sz.assign(n,1); iota(p.begin(),p.end(),0);} int find(int x){while(p[x]!=x){p[x]=p[p[x]]; x=p[x];}return x;} void unite(int a,int b){a=find(a);b=find(b); if(a==b)return; if(sz[a]<sz[b])swap(a,b); p[b]=a; sz[a]+=sz[b];}};

static bool nondeg(const vector<P>&v,const F&t){
    if(t.a==t.b||t.b==t.c||t.c==t.a) return false;
    P cr=crossp(v[t.b]-v[t.a], v[t.c]-v[t.a]);
    double bb=0; for(auto&p:v) bb=max(bb,norm2(p));
    return norm2(cr) > max(1e-28, bb*1e-26);
}

static unsigned long long ekey(int a,int b){ if(a>b) swap(a,b); return (unsigned long long)(unsigned int)a<<32 | (unsigned int)b; }

static bool validate_manifold(const Mesh&m){
    int n=m.v.size(), nf=m.f.size();
    if(n<4 || nf<4) return false;
    vector<vector<int>> incident(n);
    unordered_map<unsigned long long, vector<int>> edgeFaces;
    edgeFaces.reserve(nf*3+7);
    set<array<int,3>> seen;
    for(int i=0;i<nf;i++){
        F t=m.f[i];
        if(t.a<0||t.a>=n||t.b<0||t.b>=n||t.c<0||t.c>=n) return false;
        if(!nondeg(m.v,t)) return false;
        array<int,3> s={t.a,t.b,t.c}; sort(s.begin(),s.end());
        if(!seen.insert(s).second) return false;
        incident[t.a].push_back(i); incident[t.b].push_back(i); incident[t.c].push_back(i);
        edgeFaces[ekey(t.a,t.b)].push_back(i);
        edgeFaces[ekey(t.b,t.c)].push_back(i);
        edgeFaces[ekey(t.c,t.a)].push_back(i);
    }
    for(auto &kv:edgeFaces) if(kv.second.size()!=2) return false;
    for(int v=0; v<n; ++v){
        auto &inc=incident[v];
        if(inc.empty()) return false;
        unordered_map<int,int> loc; loc.reserve(inc.size()*2+1);
        for(int i=0;i<(int)inc.size();++i) loc[inc[i]]=i;
        DSU d(inc.size());
        for(auto &kv:edgeFaces){
            unsigned long long k=kv.first; int a=(int)(k>>32), b=(int)(k&0xffffffffu);
            if(a!=v && b!=v) continue;
            int x=kv.second[0], y=kv.second[1];
            auto ix=loc.find(x), iy=loc.find(y);
            if(ix!=loc.end() && iy!=loc.end()) d.unite(ix->second, iy->second);
        }
        int r=d.find(0); for(int i=1;i<(int)inc.size();++i) if(d.find(i)!=r) return false;
    }
    return true;
}

struct GridNN{
    double h, inv; P mn; unordered_map<long long, vector<int>> mp; const vector<P>* pts;
    static long long mix(int x,int y,int z){
        long long X=(long long)(x+1048576)&2097151LL, Y=(long long)(y+1048576)&2097151LL, Z=(long long)(z+1048576)&2097151LL;
        return (X<<42) ^ (Y<<21) ^ Z;
    }
    int ix(double x) const { return (int)floor((x-mn.x)*inv); }
    int iy(double y) const { return (int)floor((y-mn.y)*inv); }
    int iz(double z) const { return (int)floor((z-mn.z)*inv); }
    GridNN(const vector<P>&p,double radius){ pts=&p; h=max(radius,1e-12); inv=1.0/h; mn={1e100,1e100,1e100}; for(auto&q:p){mn.x=min(mn.x,q.x);mn.y=min(mn.y,q.y);mn.z=min(mn.z,q.z);} mp.reserve(p.size()*2+1); for(int i=0;i<(int)p.size();++i) mp[mix(ix(p[i].x),iy(p[i].y),iz(p[i].z))].push_back(i); }
    double nearest2(const P&q,double best) const{
        int x=ix(q.x), y=iy(q.y), z=iz(q.z); int R=2;
        for(int dx=-R;dx<=R;dx++)for(int dy=-R;dy<=R;dy++)for(int dz=-R;dz<=R;dz++){
            auto it=mp.find(mix(x+dx,y+dy,z+dz)); if(it==mp.end()) continue;
            for(int id:it->second) best=min(best,d2(q,(*pts)[id]));
        }
        return best;
    }
};

static bool hausdorff_vertices_ok(const Mesh&orig,const Mesh&cand,double eps){
    double e2=eps*eps*(1.0000001);
    if(cand.v.empty()) return false;
    GridNN gC(cand.v, max(eps,1e-9));
    for(const auto&p:orig.v) if(gC.nearest2(p,1e100)>e2) return false;
    GridNN gO(orig.v, max(eps,1e-9));
    for(const auto&p:cand.v) if(gO.nearest2(p,1e100)>e2) return false;
    return true;
}

static Mesh compact_mesh(const Mesh&m){
    int n=m.v.size(); vector<int> used(n,0), mapv(n,-1); for(auto&t:m.f){used[t.a]=used[t.b]=used[t.c]=1;}
    Mesh r; for(int i=0;i<n;i++) if(used[i]){mapv[i]=r.v.size(); r.v.push_back(m.v[i]);}
    r.f.reserve(m.f.size()); for(auto&t:m.f) r.f.push_back({mapv[t.a],mapv[t.b],mapv[t.c]}); return r;
}

static bool bbox_candidate(const Mesh&orig,double eps,Mesh&out){
    P mn={1e100,1e100,1e100}, mx={-1e100,-1e100,-1e100};
    for(auto&p:orig.v){mn.x=min(mn.x,p.x);mn.y=min(mn.y,p.y);mn.z=min(mn.z,p.z);mx.x=max(mx.x,p.x);mx.y=max(mx.y,p.y);mx.z=max(mx.z,p.z);}    
    vector<P> vv={{mn.x,mn.y,mn.z},{mx.x,mn.y,mn.z},{mx.x,mx.y,mn.z},{mn.x,mx.y,mn.z},{mn.x,mn.y,mx.z},{mx.x,mn.y,mx.z},{mx.x,mx.y,mx.z},{mn.x,mx.y,mx.z}};
    Mesh box; box.v=vv;
    box.f={{0,2,1},{0,3,2},{4,5,6},{4,6,7},{0,1,5},{0,5,4},{1,2,6},{1,6,5},{2,3,7},{2,7,6},{3,0,4},{3,4,7}};
    if(validate_manifold(box) && hausdorff_vertices_ok(orig,box,eps)){ out=box; return true; }
    return false;
}

struct Key3{int x,y,z; bool operator==(Key3 const&o)const{return x==o.x&&y==o.y&&z==o.z;}};
struct KeyHash{ size_t operator()(Key3 const&k)const{ uint64_t h=1469598103934665603ull; auto add=[&](uint64_t v){h^=v+0x9e3779b97f4a7c15ull+(h<<6)+(h>>2);}; add((uint32_t)k.x); add((uint32_t)k.y); add((uint32_t)k.z); return (size_t)h; }};

static bool clustered_candidate(const Mesh&orig,double cell,double eps,Mesh&out){
    int n=orig.v.size(); if(n<16) return false;
    P mn={1e100,1e100,1e100}; for(auto&p:orig.v){mn.x=min(mn.x,p.x);mn.y=min(mn.y,p.y);mn.z=min(mn.z,p.z);} double inv=1.0/max(cell,1e-12);
    unordered_map<Key3,int,KeyHash> cid; cid.reserve(n*2+1);
    vector<vector<int>> groups;
    vector<int> cls(n);
    for(int i=0;i<n;i++){
        Key3 k{(int)floor((orig.v[i].x-mn.x)*inv),(int)floor((orig.v[i].y-mn.y)*inv),(int)floor((orig.v[i].z-mn.z)*inv)};
        auto it=cid.find(k); if(it==cid.end()){int id=groups.size(); cid[k]=id; groups.push_back({}); it=cid.find(k);} cls[i]=it->second; groups[cls[i]].push_back(i);
    }
    if((int)groups.size()>=n) return false;
    vector<int> rep(groups.size()); vector<P> nv(groups.size());
    for(int g=0;g<(int)groups.size();++g){
        P avg{0,0,0}; for(int id:groups[g]) avg=avg+orig.v[id]; avg=avg*(1.0/groups[g].size());
        int best=groups[g][0]; double bd=d2(orig.v[best],avg); for(int id:groups[g]){double dd=d2(orig.v[id],avg); if(dd<bd){bd=dd; best=id;}}
        rep[g]=best; nv[g]=orig.v[best];
    }
    Mesh m; m.v=nv; m.f.reserve(orig.f.size());
    set<array<int,3>> used;
    for(auto&t:orig.f){ int a=cls[t.a],b=cls[t.b],c=cls[t.c]; if(a==b||b==c||c==a) continue; array<int,3>s={a,b,c}; sort(s.begin(),s.end()); if(used.insert(s).second) m.f.push_back({a,b,c}); }
    m=compact_mesh(m);
    if(m.v.size()>=orig.v.size() || m.f.empty()) return false;
    if(!validate_manifold(m)) return false;
    if(!hausdorff_vertices_ok(orig,m,eps)) return false;
    out=m; return true;
}

static bool edge_length_filter_candidate(const Mesh&orig,double thresh,double eps,Mesh&out){
    int n=orig.v.size(); vector<int> parent(n); iota(parent.begin(),parent.end(),0);
    function<int(int)> find=[&](int x){return parent[x]==x?x:parent[x]=find(parent[x]);};
    auto uni=[&](int a,int b){a=find(a);b=find(b); if(a!=b) parent[b]=a;};
    for(auto&t:orig.f){
        if(d2(orig.v[t.a],orig.v[t.b])<thresh*thresh) uni(t.a,t.b);
        if(d2(orig.v[t.b],orig.v[t.c])<thresh*thresh) uni(t.b,t.c);
        if(d2(orig.v[t.c],orig.v[t.a])<thresh*thresh) uni(t.c,t.a);
    }
    unordered_map<int,int> rid; vector<vector<int>> groups; vector<int> cls(n);
    for(int i=0;i<n;i++){int r=find(i); auto it=rid.find(r); if(it==rid.end()){int id=groups.size(); rid[r]=id; groups.push_back({}); it=rid.find(r);} cls[i]=it->second; groups[cls[i]].push_back(i);}    
    if((int)groups.size()>=n) return false;
    Mesh m; m.v.resize(groups.size());
    for(int g=0;g<(int)groups.size();g++){ P avg{0,0,0}; for(int id:groups[g]) avg=avg+orig.v[id]; avg=avg*(1.0/groups[g].size()); int best=groups[g][0]; double bd=d2(orig.v[best],avg); for(int id:groups[g]){double dd=d2(orig.v[id],avg); if(dd<bd){bd=dd; best=id;}} m.v[g]=orig.v[best]; }
    set<array<int,3>> used;
    for(auto&t:orig.f){int a=cls[t.a],b=cls[t.b],c=cls[t.c]; if(a==b||b==c||c==a) continue; array<int,3>s={a,b,c}; sort(s.begin(),s.end()); if(used.insert(s).second) m.f.push_back({a,b,c});}
    m=compact_mesh(m);
    if(m.v.size()>=orig.v.size()) return false;
    if(!validate_manifold(m)) return false;
    if(!hausdorff_vertices_ok(orig,m,eps)) return false;
    out=m; return true;
}

int main(){
    ios::sync_with_stdio(false); cin.tie(nullptr);
    int V,FN; if(!(cin>>V>>FN)) return 0;
    Mesh orig; orig.v.reserve(V); orig.f.reserve(FN);
    string tag;
    for(int i=0;i<V;i++){ cin>>tag; P p; cin>>p.x>>p.y>>p.z; orig.v.push_back(p); }
    for(int i=0;i<FN;i++){ cin>>tag; int a,b,c; cin>>a>>b>>c; --a;--b;--c; orig.f.push_back({a,b,c}); }
    P mn={1e100,1e100,1e100}, mx={-1e100,-1e100,-1e100};
    for(auto&p:orig.v){mn.x=min(mn.x,p.x);mn.y=min(mn.y,p.y);mn.z=min(mn.z,p.z);mx.x=max(mx.x,p.x);mx.y=max(mx.y,p.y);mx.z=max(mx.z,p.z);}    
    double diag=sqrt(d2(mn,mx)); double eps=max(1e-12,0.05*diag);
    Mesh best=orig, cand;
    if(bbox_candidate(orig,eps,cand) && cand.v.size()<best.v.size()) best=cand;
    vector<double> factors={1.00,0.90,0.80,0.70,0.62,0.55,0.49,0.43,0.38,0.34,0.30,0.27,0.24,0.21,0.18,0.16,0.14,0.12};
    for(double f:factors){ if(clustered_candidate(orig,eps*f,eps,cand) && cand.v.size()<best.v.size()) best=cand; }
    vector<double> ef={0.80,0.65,0.52,0.42,0.34,0.28,0.22,0.18,0.14,0.11,0.09};
    for(double f:ef){ if(edge_length_filter_candidate(orig,eps*f,eps,cand) && cand.v.size()<best.v.size()) best=cand; }
    cout<<best.v.size()<<" "<<best.f.size()<<"\n";
    cout.setf(ios::fixed); cout<<setprecision(10);
    for(auto&p:best.v) cout<<"v "<<p.x<<" "<<p.y<<" "<<p.z<<"\n";
    for(auto&t:best.f) cout<<"f "<<t.a+1<<" "<<t.b+1<<" "<<t.c+1<<"\n";
    return 0;
}