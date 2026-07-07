#include <bits/stdc++.h>
using namespace std;

struct P{double x,y,z;};
struct F{int a,b,c;};
static inline P operator+(const P&a,const P&b){return {a.x+b.x,a.y+b.y,a.z+b.z};}
static inline P operator-(const P&a,const P&b){return {a.x-b.x,a.y-b.y,a.z-b.z};}
static inline P operator*(const P&a,double s){return {a.x*s,a.y*s,a.z*s};}
static inline double dotp(const P&a,const P&b){return a.x*b.x+a.y*b.y+a.z*b.z;}
static inline P crossp(const P&a,const P&b){return {a.y*b.z-a.z*b.y,a.z*b.x-a.x*b.z,a.x*b.y-a.y*b.x};}
static inline double n2(const P&a){return dotp(a,a);} 
static inline double dist2(const P&a,const P&b){return n2(a-b);} 
static inline double norm(const P&a){return sqrt(max(0.0,n2(a)));}
static inline P unit(P a){double l=norm(a); if(l<=1e-300) return {0,0,0}; return a*(1.0/l);} 

struct Mesh{
    vector<P> v;
    vector<F> f;
    vector<vector<int>> mem;
};

static unsigned long long ekey(int a,int b){ if(a>b) swap(a,b); return (unsigned long long)(unsigned int)a<<32 | (unsigned int)b; }
static array<int,3> skey(int a,int b,int c){ array<int,3>s={a,b,c}; sort(s.begin(),s.end()); return s; }

struct GridNN{
    double h, inv; P mn; const vector<P>* pts; unordered_map<long long, vector<int>> mp;
    static long long key(int x,int y,int z){
        uint64_t X=(uint32_t)(x*73856093), Y=(uint32_t)(y*19349663), Z=(uint32_t)(z*83492791);
        return (long long)(X ^ (Y<<1) ^ (Z<<2));
    }
    GridNN(){}
    GridNN(const vector<P>&p,double radius){init(p,radius);} 
    void init(const vector<P>&p,double radius){
        pts=&p; h=max(radius,1e-9); inv=1.0/h; mn={1e100,1e100,1e100};
        for(auto&q:p){mn.x=min(mn.x,q.x); mn.y=min(mn.y,q.y); mn.z=min(mn.z,q.z);} 
        mp.clear(); mp.reserve(p.size()*2+3);
        for(int i=0;i<(int)p.size();++i){ int x=(int)floor((p[i].x-mn.x)*inv), y=(int)floor((p[i].y-mn.y)*inv), z=(int)floor((p[i].z-mn.z)*inv); mp[key(x,y,z)].push_back(i); }
    }
    double nearest2(const P&q,double best=1e300) const{
        int x=(int)floor((q.x-mn.x)*inv), y=(int)floor((q.y-mn.y)*inv), z=(int)floor((q.z-mn.z)*inv);
        int R=2;
        for(int dx=-R;dx<=R;dx++) for(int dy=-R;dy<=R;dy++) for(int dz=-R;dz<=R;dz++){
            auto it=mp.find(key(x+dx,y+dy,z+dz)); if(it==mp.end()) continue;
            for(int id:it->second) best=min(best,dist2(q,(*pts)[id]));
        }
        return best;
    }
};

static bool face_ok(const vector<P>&v,const F&t){
    if(t.a<0||t.b<0||t.c<0||t.a>=(int)v.size()||t.b>=(int)v.size()||t.c>=(int)v.size()) return false;
    if(t.a==t.b||t.b==t.c||t.c==t.a) return false;
    P cr=crossp(v[t.b]-v[t.a],v[t.c]-v[t.a]);
    double scale=1.0; for(auto&p:v) scale=max(scale,n2(p));
    return n2(cr)>scale*1e-26;
}

static bool validate_manifold(const Mesh&m){
    int n=m.v.size(), nf=m.f.size(); if(n<4||nf<4) return false;
    vector<vector<int>> inc(n);
    unordered_map<unsigned long long, array<int,3>> ec; ec.reserve(nf*4+7);
    set<array<int,3>> sf;
    for(int i=0;i<nf;i++){
        F t=m.f[i]; if(!face_ok(m.v,t)) return false;
        auto s=skey(t.a,t.b,t.c); if(!sf.insert(s).second) return false;
        inc[t.a].push_back(i); inc[t.b].push_back(i); inc[t.c].push_back(i);
        unsigned long long ks[3]={ekey(t.a,t.b),ekey(t.b,t.c),ekey(t.c,t.a)};
        for(auto k:ks){ auto &a=ec[k]; if(a[0]==0){a[1]=i;} else if(a[0]==1){a[2]=i;} a[0]++; if(a[0]>2) return false; }
    }
    for(auto &kv:ec) if(kv.second[0]!=2) return false;
    for(int v=0;v<n;v++){
        if(inc[v].empty()) return false;
        unordered_map<int,int> loc; loc.reserve(inc[v].size()*2+1);
        for(int i=0;i<(int)inc[v].size();i++) loc[inc[v][i]]=i;
        vector<int> par(inc[v].size()); iota(par.begin(),par.end(),0);
        function<int(int)> fd=[&](int x){return par[x]==x?x:par[x]=fd(par[x]);};
        auto un=[&](int a,int b){a=fd(a);b=fd(b);if(a!=b)par[b]=a;};
        for(auto &kv:ec){
            int a=(int)(kv.first>>32), b=(int)(kv.first&0xffffffffu); if(a!=v&&b!=v) continue;
            auto x=loc.find(kv.second[1]), y=loc.find(kv.second[2]); if(x==loc.end()||y==loc.end()) return false; un(x->second,y->second);
        }
        int r=fd(0); for(int i=1;i<(int)par.size();i++) if(fd(i)!=r) return false;
    }
    return true;
}

static Mesh compact(const Mesh&m){
    int n=m.v.size(); vector<int> used(n,0), mp(n,-1); for(auto&t:m.f){used[t.a]=used[t.b]=used[t.c]=1;}
    Mesh r; r.v.reserve(n); r.mem.reserve(n);
    for(int i=0;i<n;i++) if(used[i]){mp[i]=r.v.size(); r.v.push_back(m.v[i]); if(!m.mem.empty()) r.mem.push_back(m.mem[i]);}
    r.f.reserve(m.f.size()); for(auto&t:m.f) r.f.push_back({mp[t.a],mp[t.b],mp[t.c]}); return r;
}

static bool hausdorff_vertices_ok(const vector<P>&orig,const Mesh&m,double eps){
    if(m.v.empty()) return false; double e2=eps*eps*1.000001;
    GridNN gm(m.v,max(eps,1e-9)); for(auto&p:orig) if(gm.nearest2(p)>e2) return false;
    GridNN go(orig,max(eps,1e-9)); for(auto&p:m.v) if(go.nearest2(p)>e2) return false;
    return true;
}

static bool bbox_if_exact_or_safe(const vector<P>&orig,const Mesh&in,double eps,Mesh&out){
    P mn={1e100,1e100,1e100}, mx={-1e100,-1e100,-1e100};
    for(auto&p:orig){mn.x=min(mn.x,p.x);mn.y=min(mn.y,p.y);mn.z=min(mn.z,p.z);mx.x=max(mx.x,p.x);mx.y=max(mx.y,p.y);mx.z=max(mx.z,p.z);} 
    Mesh b; b.v={{mn.x,mn.y,mn.z},{mx.x,mn.y,mn.z},{mx.x,mx.y,mn.z},{mn.x,mx.y,mn.z},{mn.x,mn.y,mx.z},{mx.x,mn.y,mx.z},{mx.x,mx.y,mx.z},{mn.x,mx.y,mx.z}};
    b.f={{0,2,1},{0,3,2},{4,5,6},{4,6,7},{0,1,5},{0,5,4},{1,2,6},{1,6,5},{2,3,7},{2,7,6},{3,0,4},{3,4,7}};
    b.mem.resize(8);
    for(int i=0;i<(int)orig.size();i++){
        int best=0; double bd=dist2(orig[i],b.v[0]); for(int j=1;j<8;j++){double dd=dist2(orig[i],b.v[j]); if(dd<bd){bd=dd; best=j;}}
        b.mem[best].push_back(i);
    }
    if(validate_manifold(b) && hausdorff_vertices_ok(orig,b,eps)){ out=b; return true; }
    return false;
}

struct Topo{
    vector<vector<int>> nb;
    unordered_map<unsigned long long, pair<int,int>> opp;
    vector<P> fn;
};
static Topo build_topo(const Mesh&m){
    int n=m.v.size(); Topo T; T.nb.assign(n,{}); T.fn.resize(m.f.size()); unordered_map<unsigned long long, vector<int>> eopp; eopp.reserve(m.f.size()*4+7);
    vector<unordered_set<int>> ns(n);
    for(int i=0;i<(int)m.f.size();i++){
        auto t=m.f[i]; T.fn[i]=unit(crossp(m.v[t.b]-m.v[t.a],m.v[t.c]-m.v[t.a]));
        ns[t.a].insert(t.b); ns[t.a].insert(t.c); ns[t.b].insert(t.a); ns[t.b].insert(t.c); ns[t.c].insert(t.a); ns[t.c].insert(t.b);
        eopp[ekey(t.a,t.b)].push_back(t.c); eopp[ekey(t.b,t.c)].push_back(t.a); eopp[ekey(t.c,t.a)].push_back(t.b);
    }
    for(int i=0;i<n;i++){ T.nb[i].assign(ns[i].begin(),ns[i].end()); sort(T.nb[i].begin(),T.nb[i].end()); }
    T.opp.reserve(eopp.size()*2+1); for(auto &kv:eopp) if(kv.second.size()==2) T.opp[kv.first]={kv.second[0],kv.second[1]};
    return T;
}

static bool link_condition(const Topo&T,int u,int v){
    auto it=T.opp.find(ekey(u,v)); if(it==T.opp.end()) return false;
    int o1=it->second.first,o2=it->second.second; if(o1==o2) return false;
    vector<int> inter; const auto&A=T.nb[u]; const auto&B=T.nb[v];
    set_intersection(A.begin(),A.end(),B.begin(),B.end(),back_inserter(inter));
    sort(inter.begin(),inter.end()); vector<int> need={o1,o2}; sort(need.begin(),need.end());
    return inter==need;
}

static double vertex_flatness(const Mesh&m,const Topo&T,int v){
    P avg{0,0,0}; int cnt=0;
    for(int i=0;i<(int)m.f.size();i++){auto f=m.f[i]; if(f.a==v||f.b==v||f.c==v){avg=avg+T.fn[i]; cnt++;}}
    if(cnt==0) return 1.0; avg=unit(avg); double worst=0;
    for(int i=0;i<(int)m.f.size();i++){auto f=m.f[i]; if(f.a==v||f.b==v||f.c==v) worst=max(worst,1.0-dotp(avg,T.fn[i]));}
    return worst;
}

static bool best_representative(const vector<P>&orig,const Mesh&m,int u,int v,double eps,int &repOrig,P &repP, vector<int>&merged){
    merged.clear(); merged.reserve(m.mem[u].size()+m.mem[v].size());
    merged.insert(merged.end(),m.mem[u].begin(),m.mem[u].end()); merged.insert(merged.end(),m.mem[v].begin(),m.mem[v].end());
    vector<int> trials; trials.push_back(m.mem[u][0]); trials.push_back(m.mem[v][0]);
    if(merged.size()<=48) for(int id:merged) trials.push_back(id);
    else { for(int k=0;k<24;k++) trials.push_back(merged[(long long)k*merged.size()/24]); }
    sort(trials.begin(),trials.end()); trials.erase(unique(trials.begin(),trials.end()),trials.end());
    double best=1e300; int bestid=-1; double lim2=eps*eps*0.999;
    for(int cand:trials){ double mx=0; const P&p=orig[cand]; for(int id:merged){ mx=max(mx,dist2(p,orig[id])); if(mx>lim2 && mx>=best) break; } if(mx<best){best=mx; bestid=cand;} }
    if(bestid<0 || best>lim2) return false; repOrig=bestid; repP=orig[bestid]; return true;
}

static Mesh apply_pairs(const vector<P>&orig,const Mesh&m,const vector<pair<int,int>>&pairs,const vector<int>&repOrigs,const vector<vector<int>>&mergeds){
    int n=m.v.size(); vector<int> root(n); iota(root.begin(),root.end(),0); vector<int> pairIndex(n,-1);
    for(int i=0;i<(int)pairs.size();i++){ int u=pairs[i].first,v=pairs[i].second; root[u]=v; pairIndex[v]=i; pairIndex[u]=i; }
    Mesh r; vector<int> mp(n,-1); r.v.reserve(n-pairs.size()); r.mem.reserve(n-pairs.size());
    vector<char> killed(n,0); for(auto&p:pairs) killed[p.first]=1;
    for(int i=0;i<n;i++) if(!killed[i]){
        mp[i]=r.v.size(); int pi=pairIndex[i];
        if(pi>=0){ r.v.push_back(orig[repOrigs[pi]]); r.mem.push_back(mergeds[pi]); }
        else { r.v.push_back(m.v[i]); r.mem.push_back(m.mem[i]); }
    }
    set<array<int,3>> seen;
    for(auto&t:m.f){
        int a=root[t.a], b=root[t.b], c=root[t.c]; if(a==b||b==c||c==a) continue;
        a=mp[a]; b=mp[b]; c=mp[c]; if(a<0||b<0||c<0||a==b||b==c||c==a) continue;
        F nf{a,b,c}; if(!face_ok(r.v,nf)) continue;
        auto s=skey(a,b,c); if(seen.insert(s).second) r.f.push_back(nf);
    }
    return compact(r);
}

static bool collapse_round(const vector<P>&orig,const Mesh&m,double eps,double lengthFactor,double curveLimit,Mesh&out){
    Topo T=build_topo(m); int n=m.v.size();
    struct E{double score,len; int u,v;}; vector<E> edges; edges.reserve(m.f.size()*3);
    set<unsigned long long> usedE;
    for(auto&t:m.f){ int aa[3]={t.a,t.b,t.c}; for(int k=0;k<3;k++){int u=aa[k],v=aa[(k+1)%3]; auto key=ekey(u,v); if(!usedE.insert(key).second) continue; double l=sqrt(dist2(m.v[u],m.v[v])); if(l<=eps*lengthFactor){ double cu=vertex_flatness(m,T,u), cv=vertex_flatness(m,T,v); edges.push_back({l*(1.0+2.5*max(cu,cv)),l,u,v}); } } }
    sort(edges.begin(),edges.end(),[](const E&a,const E&b){return a.score<b.score;});
    vector<char> taken(n,0); vector<pair<int,int>> pairs; vector<int> reps; vector<vector<int>> mergeds;
    int maxPairs=max(1,n/5);
    for(auto&e:edges){
        int u=e.u,v=e.v; if(taken[u]||taken[v]) continue; if(!link_condition(T,u,v)) continue;
        double curv=max(vertex_flatness(m,T,u),vertex_flatness(m,T,v)); if(curv>curveLimit && e.len>eps*0.20) continue;
        int rep; P rp; vector<int> merged; if(!best_representative(orig,m,u,v,eps,rep,rp,merged)) continue;
        if(T.nb[u].size()>T.nb[v].size()) swap(u,v);
        pairs.push_back({u,v}); reps.push_back(rep); mergeds.push_back(std::move(merged)); taken[u]=taken[v]=1;
        if((int)pairs.size()>=maxPairs) break;
    }
    if(pairs.empty()) return false;
    while(!pairs.empty()){
        Mesh c=apply_pairs(orig,m,pairs,reps,mergeds);
        if(c.v.size()<m.v.size() && validate_manifold(c) && hausdorff_vertices_ok(orig,c,eps)){ out=c; return true; }
        int keep=(int)pairs.size()/2; pairs.resize(keep); reps.resize(keep); mergeds.resize(keep);
    }
    return false;
}

static Mesh edge_collapse_search(const vector<P>&orig,const Mesh&start,double eps){
    Mesh best=start;
    vector<double> lf={0.98,0.82,0.68,0.56,0.46,0.38,0.31,0.25,0.20,0.16,0.13,0.10};
    vector<double> cl={0.08,0.12,0.18,0.25,0.35,0.50,0.75,1.10,2.00};
    for(double curve:cl){
        Mesh cur=start; bool improved=true; int passes=0;
        while(improved && passes<18){
            improved=false; passes++;
            for(double f:lf){
                Mesh nxt; if(collapse_round(orig,cur,eps,f,curve,nxt)){
                    cur=nxt; improved=true;
                    if(cur.v.size()<best.v.size()) best=cur;
                    break;
                }
            }
        }
    }
    return best;
}

static bool voxel_cluster_candidate(const vector<P>&orig,const Mesh&in,double eps,double cell,Mesh&out){
    struct K{int x,y,z; bool operator==(K const&o)const{return x==o.x&&y==o.y&&z==o.z;}};
    struct H{size_t operator()(K const&k)const{return (uint64_t)(uint32_t)(k.x*73856093) ^ ((uint64_t)(uint32_t)(k.y*19349663)<<1) ^ ((uint64_t)(uint32_t)(k.z*83492791)<<2);} };
    P mn={1e100,1e100,1e100}; for(auto&p:in.v){mn.x=min(mn.x,p.x);mn.y=min(mn.y,p.y);mn.z=min(mn.z,p.z);} double inv=1.0/max(cell,1e-12);
    unordered_map<K,int,H> id; id.reserve(in.v.size()*2+1); vector<vector<int>> verts;
    vector<int> cls(in.v.size());
    for(int i=0;i<(int)in.v.size();i++){
        K k{(int)floor((in.v[i].x-mn.x)*inv),(int)floor((in.v[i].y-mn.y)*inv),(int)floor((in.v[i].z-mn.z)*inv)};
        auto it=id.find(k); if(it==id.end()){int nid=verts.size(); id[k]=nid; verts.push_back({}); it=id.find(k);} cls[i]=it->second; verts[cls[i]].push_back(i);
    }
    if(verts.size()>=in.v.size()) return false;
    Mesh m; m.v.resize(verts.size()); m.mem.resize(verts.size());
    for(int g=0;g<(int)verts.size();g++){
        vector<int> all; for(int v:verts[g]) all.insert(all.end(),in.mem[v].begin(),in.mem[v].end());
        vector<int> trials=all; if(trials.size()>40){ vector<int> t; for(int k=0;k<40;k++) t.push_back(trials[(long long)k*trials.size()/40]); trials.swap(t); }
        double best=1e300; int bid=all[0]; for(int cand:trials){double mx=0; for(int oid:all){mx=max(mx,dist2(orig[cand],orig[oid])); if(mx>eps*eps) break;} if(mx<best){best=mx; bid=cand;}}
        if(best>eps*eps*0.999) return false; m.v[g]=orig[bid]; m.mem[g]=std::move(all);
    }
    set<array<int,3>> seen;
    for(auto&t:in.f){int a=cls[t.a],b=cls[t.b],c=cls[t.c]; if(a==b||b==c||c==a) continue; F nf{a,b,c}; if(!face_ok(m.v,nf)) continue; auto s=skey(a,b,c); if(seen.insert(s).second) m.f.push_back(nf);} 
    m=compact(m);
    if(m.v.size()<in.v.size() && validate_manifold(m) && hausdorff_vertices_ok(orig,m,eps)){out=m; return true;} return false;
}

int main(){
    ios::sync_with_stdio(false); cin.tie(nullptr);
    int V,FN; if(!(cin>>V>>FN)) return 0;
    Mesh orig; orig.v.reserve(V); orig.f.reserve(FN); orig.mem.resize(V);
    string tag; for(int i=0;i<V;i++){P p; cin>>tag>>p.x>>p.y>>p.z; orig.v.push_back(p); orig.mem[i]={i};}
    for(int i=0;i<FN;i++){int a,b,c; cin>>tag>>a>>b>>c; --a;--b;--c; orig.f.push_back({a,b,c});}
    P mn={1e100,1e100,1e100}, mx={-1e100,-1e100,-1e100}; for(auto&p:orig.v){mn.x=min(mn.x,p.x);mn.y=min(mn.y,p.y);mn.z=min(mn.z,p.z);mx.x=max(mx.x,p.x);mx.y=max(mx.y,p.y);mx.z=max(mx.z,p.z);} double eps=0.05*sqrt(dist2(mn,mx)); eps=max(eps,1e-12);
    Mesh base=compact(orig); if(!validate_manifold(base)){ base=orig; }
    Mesh best=base, cand;
    if(bbox_if_exact_or_safe(orig.v,base,eps,cand) && cand.v.size()<best.v.size()) best=cand;
    for(double cf: {0.95,0.80,0.66,0.54,0.44,0.36,0.29,0.23,0.18,0.14}){
        Mesh vc; if(voxel_cluster_candidate(orig.v,base,eps,eps*cf,vc)){
            Mesh ec=edge_collapse_search(orig.v,vc,eps); if(ec.v.size()<best.v.size()) best=ec; else if(vc.v.size()<best.v.size()) best=vc;
        }
    }
    Mesh ec=edge_collapse_search(orig.v,base,eps); if(ec.v.size()<best.v.size()) best=ec;
    if(!validate_manifold(best) || !hausdorff_vertices_ok(orig.v,best,eps)) best=base;
    cout<<best.v.size()<<" "<<best.f.size()<<"\n";
    cout.setf(ios::fixed); cout<<setprecision(10);
    for(auto&p:best.v) cout<<"v "<<p.x<<" "<<p.y<<" "<<p.z<<"\n";
    for(auto&t:best.f) cout<<"f "<<t.a+1<<" "<<t.b+1<<" "<<t.c+1<<"\n";
    return 0;
}