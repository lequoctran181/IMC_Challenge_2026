#include <bits/stdc++.h>
using namespace std;

struct Vec{
    double x=0,y=0,z=0;
    Vec(){} Vec(double X,double Y,double Z):x(X),y(Y),z(Z){}
    Vec operator+(const Vec&o)const{return {x+o.x,y+o.y,z+o.z};}
    Vec operator-(const Vec&o)const{return {x-o.x,y-o.y,z-o.z};}
    Vec operator*(double s)const{return {x*s,y*s,z*s};}
    Vec operator/(double s)const{return {x/s,y/s,z/s};}
};
static inline double dotp(const Vec&a,const Vec&b){return a.x*b.x+a.y*b.y+a.z*b.z;}
static inline Vec crossp(const Vec&a,const Vec&b){return {a.y*b.z-a.z*b.y,a.z*b.x-a.x*b.z,a.x*b.y-a.y*b.x};}
static inline double norm2(const Vec&a){return dotp(a,a);} 
static inline double dist2(const Vec&a,const Vec&b){return norm2(a-b);} 
struct Tri{int a,b,c;};

static double tri_area2(const Vec&a,const Vec&b,const Vec&c){ return norm2(crossp(b-a,c-a)); }

struct Mesh{ vector<Vec> p; vector<Tri> f; };

static uint64_t keyEdge(int a,int b){ if(a>b) swap(a,b); return (uint64_t)(uint32_t)a<<32 | (uint32_t)b; }

static bool validate_closed_manifold(const Mesh&m){
    int n=m.p.size(); if(n<4 || m.f.size()<4) return false;
    vector<int> used(n,0);
    unordered_map<uint64_t,int> ec; ec.reserve(m.f.size()*3+10);
    unordered_map<uint64_t,int> dir; dir.reserve(m.f.size()*3+10);
    for(auto&t:m.f){
        if(t.a<0||t.a>=n||t.b<0||t.b>=n||t.c<0||t.c>=n) return false;
        if(t.a==t.b||t.b==t.c||t.c==t.a) return false;
        if(tri_area2(m.p[t.a],m.p[t.b],m.p[t.c]) < 1e-28) return false;
        used[t.a]=used[t.b]=used[t.c]=1;
        int v[3]={t.a,t.b,t.c};
        for(int i=0;i<3;i++){
            int a=v[i], b=v[(i+1)%3];
            ec[keyEdge(a,b)]++;
            uint64_t dk=((uint64_t)(uint32_t)a<<32)|(uint32_t)b;
            dir[dk]++;
            if(dir[dk]>1) return false;
        }
    }
    for(auto &kv:ec) if(kv.second!=2) return false;
    vector<vector<pair<int,int>>> inc(n);
    for(auto&t:m.f){
        inc[t.a].push_back({t.b,t.c}); inc[t.b].push_back({t.c,t.a}); inc[t.c].push_back({t.a,t.b});
    }
    for(int v=0; v<n; ++v) if(used[v]){
        unordered_map<int,vector<int>> adj;
        for(auto pr:inc[v]){ adj[pr.first].push_back(pr.second); adj[pr.second].push_back(pr.first); }
        if(adj.size()<3) return false;
        int start=adj.begin()->first;
        for(auto &kv:adj) if(kv.second.size()!=2) return false;
        unordered_set<int> seen; vector<int> st={start}; seen.insert(start);
        while(!st.empty()){
            int u=st.back(); st.pop_back();
            for(int w:adj[u]) if(!seen.count(w)){seen.insert(w); st.push_back(w);}
        }
        if(seen.size()!=adj.size()) return false;
    }
    vector<vector<int>> g(n);
    for(auto&t:m.f){ int v[3]={t.a,t.b,t.c}; for(int i=0;i<3;i++){g[v[i]].push_back(v[(i+1)%3]);g[v[(i+1)%3]].push_back(v[i]);}}
    int s=-1, cnt=0; for(int i=0;i<n;i++) if(used[i]){cnt++; if(s<0)s=i;}
    vector<char> vis(n); queue<int>q; q.push(s); vis[s]=1; int got=0;
    while(!q.empty()){int u=q.front();q.pop(); got++; for(int v:g[u]) if(!vis[v]) vis[v]=1,q.push(v);} 
    return got==cnt;
}

static Mesh compact_mesh(const vector<Vec>&orig, const vector<Tri>&faces, const vector<int>&mapv){
    int n=orig.size(); vector<int> id(n,-1); Mesh m;
    for(int i=0;i<n;i++) if(mapv[i]==i){ id[i]=m.p.size(); m.p.push_back(orig[i]); }
    unordered_set<string> seen; seen.reserve(faces.size()*2+10);
    for(auto&t:faces){
        int a=mapv[t.a], b=mapv[t.b], c=mapv[t.c];
        if(a==b||b==c||c==a) continue;
        int ia=id[a], ib=id[b], ic=id[c];
        if(ia<0||ib<0||ic<0) continue;
        if(tri_area2(m.p[ia],m.p[ib],m.p[ic]) < 1e-28) continue;
        int mn=min({ia,ib,ic}), mx=max({ia,ib,ic}), mid=ia+ib+ic-mn-mx;
        string key=to_string(mn)+","+to_string(mid)+","+to_string(mx);
        if(seen.insert(key).second) m.f.push_back({ia,ib,ic});
    }
    vector<int> used(m.p.size(),0); for(auto&t:m.f) used[t.a]=used[t.b]=used[t.c]=1;
    vector<int> nid(m.p.size(),-1); vector<Vec> np; for(int i=0;i<(int)m.p.size();++i) if(used[i]) nid[i]=np.size(), np.push_back(m.p[i]);
    for(auto&t:m.f){t.a=nid[t.a]; t.b=nid[t.b]; t.c=nid[t.c];}
    m.p.swap(np);
    return m;
}

static bool vertex_hausdorff_ok(const vector<Vec>&orig, const Mesh&m, double eps){
    if(m.p.empty()) return false;
    double e2=eps*eps*(1.0000001);
    double cell=max(eps,1e-12); unordered_map<long long, vector<int>> grid; grid.reserve(m.p.size()*2+10);
    auto h=[&](int ix,int iy,int iz)->long long{
        uint64_t x=(uint64_t)(uint32_t)(ix*73856093), y=(uint64_t)(uint32_t)(iy*19349663), z=(uint64_t)(uint32_t)(iz*83492791);
        return (long long)(x^y^z);
    };
    for(int i=0;i<(int)m.p.size();++i){ int ix=floor(m.p[i].x/cell), iy=floor(m.p[i].y/cell), iz=floor(m.p[i].z/cell); grid[h(ix,iy,iz)].push_back(i); }
    for(auto &p:orig){
        int ix=floor(p.x/cell), iy=floor(p.y/cell), iz=floor(p.z/cell); bool ok=false;
        for(int dx=-1;dx<=1&&!ok;dx++)for(int dy=-1;dy<=1&&!ok;dy++)for(int dz=-1;dz<=1&&!ok;dz++){
            auto it=grid.find(h(ix+dx,iy+dy,iz+dz)); if(it==grid.end()) continue;
            for(int j:it->second) if(dist2(p,m.p[j])<=e2){ok=true; break;}
        }
        if(!ok) return false;
    }
    return true;
}

static vector<vector<int>> build_neighbors(int n,const vector<Tri>&faces,const vector<int>&alive,const vector<int>&mp){
    vector<vector<int>> adj(n);
    unordered_set<uint64_t> es; es.reserve(faces.size()*3+10);
    for(auto&t:faces){
        int v[3]={mp[t.a],mp[t.b],mp[t.c]};
        if(v[0]==v[1]||v[1]==v[2]||v[2]==v[0]) continue;
        for(int i=0;i<3;i++){int a=v[i],b=v[(i+1)%3]; if(!alive[a]||!alive[b]||a==b) continue; uint64_t k=keyEdge(a,b); if(es.insert(k).second){adj[a].push_back(b); adj[b].push_back(a);}}
    }
    return adj;
}

static bool link_condition(int u,int v,const vector<Tri>&faces,const vector<int>&mp){
    unordered_set<int> Nu,Nv, commonByFace;
    int commonFaces=0;
    for(auto&t:faces){
        int a=mp[t.a],b=mp[t.b],c=mp[t.c];
        if(a==b||b==c||c==a) continue;
        bool hu=(a==u||b==u||c==u), hv=(a==v||b==v||c==v);
        if(hu){ if(a!=u)Nu.insert(a); if(b!=u)Nu.insert(b); if(c!=u)Nu.insert(c); }
        if(hv){ if(a!=v)Nv.insert(a); if(b!=v)Nv.insert(b); if(c!=v)Nv.insert(c); }
        if(hu&&hv){ commonFaces++; if(a!=u&&a!=v) commonByFace.insert(a); if(b!=u&&b!=v) commonByFace.insert(b); if(c!=u&&c!=v) commonByFace.insert(c); }
    }
    if(commonFaces!=2 || commonByFace.size()!=2) return false;
    int inter=0; for(int x:Nu) if(x!=v && Nv.count(x)) inter++;
    return inter==2;
}

static Mesh try_topo_radius(const Mesh&in,double eps,double radFactor,int passLimit){
    int n=in.p.size(); vector<int> mp(n); iota(mp.begin(),mp.end(),0); vector<int> alive(n,1);
    double lim2=(eps*radFactor)*(eps*radFactor);
    vector<double> clusterMax(n,0.0);
    int changes=1, passes=0;
    while(changes && passes++<passLimit){
        changes=0;
        auto adj=build_neighbors(n,in.f,alive,mp);
        vector<tuple<double,int,int>> edges;
        for(int u=0;u<n;u++) if(alive[u]) for(int v:adj[u]) if(u<v && alive[v]){
            double d=dist2(in.p[u],in.p[v]); if(d<=lim2) edges.emplace_back(d,u,v);
        }
        sort(edges.begin(),edges.end());
        for(auto [d,u0,v0]:edges){
            int u=mp[u0], v=mp[v0]; if(u==v||!alive[u]||!alive[v]) continue;
            int keep=u, kill=v;
            if(adj[v].size()>adj[u].size()) {keep=v; kill=u;}
            if(dist2(in.p[keep],in.p[kill])>lim2) continue;
            bool okrad=true;
            for(int i=0;i<n;i++) if(mp[i]==kill){ if(dist2(in.p[i],in.p[keep]) > eps*eps*0.9604){okrad=false; break;} }
            if(!okrad) continue;
            if(!link_condition(keep,kill,in.f,mp)) continue;
            for(int i=0;i<n;i++) if(mp[i]==kill) mp[i]=keep;
            alive[kill]=0; changes++;
        }
    }
    for(int i=0;i<n;i++) while(mp[mp[i]]!=mp[i]) mp[i]=mp[mp[i]];
    Mesh out=compact_mesh(in.p,in.f,mp);
    return out;
}

static Mesh try_bbox_corner(const Mesh&in,double eps){
    Vec mn=in.p[0], mx=in.p[0];
    for(auto&p:in.p){mn.x=min(mn.x,p.x);mn.y=min(mn.y,p.y);mn.z=min(mn.z,p.z);mx.x=max(mx.x,p.x);mx.y=max(mx.y,p.y);mx.z=max(mx.z,p.z);} 
    vector<Vec> c; for(int ix=0;ix<2;ix++)for(int iy=0;iy<2;iy++)for(int iz=0;iz<2;iz++) c.push_back({ix?mx.x:mn.x,iy?mx.y:mn.y,iz?mx.z:mn.z});
    double e2=eps*eps*1.000001;
    for(auto&p:in.p){ double best=1e300; for(auto&q:c) best=min(best,dist2(p,q)); if(best>e2) return {}; }
    Mesh m; m.p=c;
    auto add=[&](int a,int b,int cc){m.f.push_back({a,b,cc});};
    add(0,2,3); add(0,3,1);
    add(4,5,7); add(4,7,6);
    add(0,1,5); add(0,5,4);
    add(2,6,7); add(2,7,3);
    add(0,4,6); add(0,6,2);
    add(1,3,7); add(1,7,5);
    if(validate_closed_manifold(m)) return m;
    return {};
}

static void print_mesh(const Mesh&m){
    cout<<m.p.size()<<" "<<m.f.size()<<"\n";
    cout.setf(ios::fixed); cout<<setprecision(9);
    for(auto&p:m.p) cout<<"v "<<p.x<<" "<<p.y<<" "<<p.z<<"\n";
    for(auto&t:m.f) cout<<"f "<<t.a+1<<" "<<t.b+1<<" "<<t.c+1<<"\n";
}

int main(){
    ios::sync_with_stdio(false); cin.tie(nullptr);
    int V,F; if(!(cin>>V>>F)) return 0;
    Mesh in; in.p.reserve(V); in.f.reserve(F);
    string tag;
    for(int i=0;i<V;i++){ double x,y,z; cin>>tag>>x>>y>>z; in.p.push_back({x,y,z}); }
    for(int i=0;i<F;i++){ int a,b,c; cin>>tag>>a>>b>>c; --a;--b;--c; in.f.push_back({a,b,c}); }
    if(V==0){ cout<<"0 0\n"; return 0; }
    Vec mn=in.p[0], mx=in.p[0];
    for(auto&p:in.p){ mn.x=min(mn.x,p.x); mn.y=min(mn.y,p.y); mn.z=min(mn.z,p.z); mx.x=max(mx.x,p.x); mx.y=max(mx.y,p.y); mx.z=max(mx.z,p.z); }
    double diag=sqrt(dist2(mn,mx)); double eps=0.05*diag;
    Mesh best=in;
    if(!validate_closed_manifold(best)){
        print_mesh(in); return 0;
    }
    Mesh box=try_bbox_corner(in,eps*0.98);
    if(!box.p.empty() && validate_closed_manifold(box) && vertex_hausdorff_ok(in.p,box,eps)) best=box;
    vector<double> factors={0.98,0.90,0.80,0.70,0.60,0.50,0.40,0.32,0.25,0.18,0.12};
    for(double rf:factors){
        Mesh cand=try_topo_radius(in,eps,rf,8);
        if(!cand.p.empty() && cand.p.size()<best.p.size() && validate_closed_manifold(cand) && vertex_hausdorff_ok(in.p,cand,eps)) best=cand;
    }
    print_mesh(best);
    return 0;
}