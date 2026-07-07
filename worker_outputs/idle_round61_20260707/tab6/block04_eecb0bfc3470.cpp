#include <bits/stdc++.h>
using namespace std;

struct Vec{double x,y,z;};
struct P2{double x,y;};
struct Tri{int a,b,c;};
struct Mesh{vector<Vec> p; vector<Tri> f;};

static inline Vec operator+(const Vec&a,const Vec&b){return {a.x+b.x,a.y+b.y,a.z+b.z};}
static inline Vec operator-(const Vec&a,const Vec&b){return {a.x-b.x,a.y-b.y,a.z-b.z};}
static inline Vec operator*(const Vec&a,double s){return {a.x*s,a.y*s,a.z*s};}
static inline double dotv(const Vec&a,const Vec&b){return a.x*b.x+a.y*b.y+a.z*b.z;}
static inline Vec crossv(const Vec&a,const Vec&b){return {a.y*b.z-a.z*b.y,a.z*b.x-a.x*b.z,a.x*b.y-a.y*b.x};}
static inline double n2(const Vec&a){return dotv(a,a);} 
static inline double dist2v(const Vec&a,const Vec&b){return n2(a-b);} 
static inline double len(const Vec&a){return sqrt(max(0.0,n2(a)));}
static Vec normed(Vec a){double l=len(a); if(l==0) return {0,0,0}; return {a.x/l,a.y/l,a.z/l};}
static uint64_t ekey(int a,int b){ if(a>b) swap(a,b); return (uint64_t)(uint32_t)a<<32 | (uint32_t)b; }
static double area2_3(const Vec&a,const Vec&b,const Vec&c){return n2(crossv(b-a,c-a));}
static double orient2(const P2&a,const P2&b,const P2&c){return (b.x-a.x)*(c.y-a.y)-(b.y-a.y)*(c.x-a.x);} 
static double areaPoly(const vector<P2>&q){double s=0; for(int i=0,n=q.size();i<n;i++){int j=(i+1)%n; s+=q[i].x*q[j].y-q[i].y*q[j].x;} return s*0.5;}

static bool pointInTri(P2 p,P2 a,P2 b,P2 c){
    double o1=orient2(a,b,p), o2=orient2(b,c,p), o3=orient2(c,a,p);
    double eps=1e-12;
    return (o1>=-eps&&o2>=-eps&&o3>=-eps)||(o1<=eps&&o2<=eps&&o3<=eps);
}
static bool pointInPoly(const vector<P2>&poly,P2 p){
    bool in=false; int n=poly.size();
    for(int i=0,j=n-1;i<n;j=i++){
        const P2&a=poly[i], &b=poly[j];
        if(((a.y>p.y)!=(b.y>p.y)) && p.x < (b.x-a.x)*(p.y-a.y)/(b.y-a.y+1e-300)+a.x) in=!in;
    }
    return in;
}

static bool validateClosed(const Mesh&m){
    int n=m.p.size(); if(n<4 || m.f.size()<4) return false;
    unordered_map<uint64_t,int> cnt; cnt.reserve(m.f.size()*4+7);
    unordered_map<uint64_t,int> dir; dir.reserve(m.f.size()*4+7);
    vector<int> used(n,0);
    for(auto&t:m.f){
        if(t.a<0||t.a>=n||t.b<0||t.b>=n||t.c<0||t.c>=n) return false;
        if(t.a==t.b||t.b==t.c||t.c==t.a) return false;
        if(area2_3(m.p[t.a],m.p[t.b],m.p[t.c])<1e-26) return false;
        used[t.a]=used[t.b]=used[t.c]=1;
        int v[3]={t.a,t.b,t.c};
        for(int i=0;i<3;i++){
            int a=v[i],b=v[(i+1)%3]; cnt[ekey(a,b)]++;
            uint64_t dk=(uint64_t)(uint32_t)a<<32 | (uint32_t)b; if(++dir[dk]>1) return false;
        }
    }
    for(auto &kv:cnt) if(kv.second!=2) return false;
    vector<vector<pair<int,int>>> inc(n);
    for(auto&t:m.f){inc[t.a].push_back({t.b,t.c}); inc[t.b].push_back({t.c,t.a}); inc[t.c].push_back({t.a,t.b});}
    for(int v=0; v<n; ++v) if(used[v]){
        unordered_map<int, vector<int>> lk; lk.reserve(inc[v].size()*2+1);
        for(auto pr:inc[v]){lk[pr.first].push_back(pr.second); lk[pr.second].push_back(pr.first);} 
        if(lk.size()<3) return false;
        for(auto &kv:lk) if(kv.second.size()!=2) return false;
        unordered_set<int> seen; vector<int> st; st.push_back(lk.begin()->first); seen.insert(st.back());
        while(!st.empty()){int u=st.back(); st.pop_back(); for(int w:lk[u]) if(!seen.count(w)){seen.insert(w); st.push_back(w);}}
        if(seen.size()!=lk.size()) return false;
    }
    vector<vector<int>> g(n);
    for(auto&t:m.f){int v[3]={t.a,t.b,t.c}; for(int i=0;i<3;i++){g[v[i]].push_back(v[(i+1)%3]); g[v[(i+1)%3]].push_back(v[i]);}}
    int s=-1, uc=0; for(int i=0;i<n;i++) if(used[i]){uc++; if(s<0)s=i;}
    vector<char> vis(n); queue<int>q; q.push(s); vis[s]=1; int got=0;
    while(!q.empty()){int u=q.front();q.pop();got++; for(int v:g[u]) if(!vis[v]) vis[v]=1,q.push(v);} 
    return got==uc;
}

struct GridNN{
    double cell; unordered_map<long long, vector<int>> mp; const vector<Vec>* pts;
    static long long h(int x,int y,int z){
        uint64_t a=(uint32_t)x*11995408973635179863ull;
        uint64_t b=(uint32_t)y*10150724397891781847ull;
        uint64_t c=(uint32_t)z*14932028262755578307ull;
        return (long long)(a^b^c);
    }
    void build(const vector<Vec>&p,double c){pts=&p; cell=max(c,1e-12); mp.clear(); mp.reserve(p.size()*2+7); for(int i=0;i<(int)p.size();++i){int ix=floor(p[i].x/cell),iy=floor(p[i].y/cell),iz=floor(p[i].z/cell);mp[h(ix,iy,iz)].push_back(i);} }
    bool anyWithin(const Vec&p,double r2)const{
        int ix=floor(p.x/cell),iy=floor(p.y/cell),iz=floor(p.z/cell); int R=2;
        for(int dx=-R;dx<=R;dx++)for(int dy=-R;dy<=R;dy++)for(int dz=-R;dz<=R;dz++){
            auto it=mp.find(h(ix+dx,iy+dy,iz+dz)); if(it==mp.end()) continue;
            for(int j:it->second) if(dist2v(p,(*pts)[j])<=r2) return true;
        }
        return false;
    }
};
static bool vertexHausdorff(const vector<Vec>&orig,const Mesh&m,double eps){
    if(m.p.empty()) return false; double r2=eps*eps*1.0000001;
    GridNN a,b; a.build(m.p,eps); for(auto&p:orig) if(!a.anyWithin(p,r2)) return false;
    b.build(orig,eps); for(auto&p:m.p) if(!b.anyWithin(p,r2)) return false;
    return true;
}

static vector<int> earClip(const vector<P2>&poly){
    int n=poly.size(); vector<int> res; if(n<3) return res;
    vector<int> V(n); iota(V.begin(),V.end(),0);
    if(areaPoly(poly)<0) reverse(V.begin(),V.end());
    int guard=0;
    while(V.size()>3 && guard++<n*n){
        bool cut=false; int m=V.size();
        for(int ii=0; ii<m; ++ii){
            int ip=V[(ii+m-1)%m], ic=V[ii], in=V[(ii+1)%m];
            if(orient2(poly[ip],poly[ic],poly[in])<=1e-14) continue;
            bool has=false;
            for(int jj=0;jj<m;jj++){int v=V[jj]; if(v==ip||v==ic||v==in) continue; if(pointInTri(poly[v],poly[ip],poly[ic],poly[in])){has=true;break;}}
            if(has) continue;
            res.push_back(ip);res.push_back(ic);res.push_back(in); V.erase(V.begin()+ii); cut=true; break;
        }
        if(!cut) break;
    }
    if(V.size()==3){res.push_back(V[0]);res.push_back(V[1]);res.push_back(V[2]);}
    if((int)res.size()/3 != n-2) res.clear();
    return res;
}

struct DSU{vector<int> p,sz; DSU(int n=0){init(n);} void init(int n){p.resize(n);sz.assign(n,1);iota(p.begin(),p.end(),0);} int find(int x){return p[x]==x?x:p[x]=find(p[x]);} void unite(int a,int b){a=find(a);b=find(b); if(a==b)return; if(sz[a]<sz[b])swap(a,b); p[b]=a; sz[a]+=sz[b];}};

struct Basis{Vec o,u,v,n; int drop;};
static Basis makeBasis(const Mesh&in,const vector<int>&faces){
    Vec n={0,0,0}, cen={0,0,0}; double asum=0;
    for(int fi:faces){auto t=in.f[fi]; Vec cr=crossv(in.p[t.b]-in.p[t.a],in.p[t.c]-in.p[t.a]); double l=len(cr); if(l>0){n=n+cr; asum+=l;} cen=cen+in.p[t.a]+in.p[t.b]+in.p[t.c];}
    if(len(n)<1e-20){auto t=in.f[faces[0]]; n=crossv(in.p[t.b]-in.p[t.a],in.p[t.c]-in.p[t.a]);}
    n=normed(n); cen=cen*(1.0/(3.0*faces.size()));
    Vec ax = (fabs(n.x)<0.8?Vec{1,0,0}:Vec{0,1,0});
    Vec u=normed(crossv(ax,n)); Vec v=crossv(n,u);
    int drop=0; double axn=fabs(n.x), ayn=fabs(n.y), azn=fabs(n.z); if(ayn>axn&&ayn>=azn) drop=1; else if(azn>axn&&azn>ayn) drop=2;
    return {cen,u,v,n,drop};
}
static P2 proj(const Basis&b,const Vec&p){Vec d=p-b.o; return {dotv(d,b.u), dotv(d,b.v)};}

static bool extractOneLoop(const Mesh&in,const vector<int>&regionFaces, vector<int>&loop, vector<int>&allVerts){
    unordered_map<uint64_t,int> edgeCnt; edgeCnt.reserve(regionFaces.size()*4+7);
    unordered_map<uint64_t,pair<int,int>> ep;
    unordered_set<int> verts;
    for(int fi:regionFaces){auto t=in.f[fi]; int a[3]={t.a,t.b,t.c}; for(int i=0;i<3;i++){verts.insert(a[i]); int u=a[i],v=a[(i+1)%3]; uint64_t k=ekey(u,v); edgeCnt[k]++; ep[k]={u,v};}}
    allVerts.assign(verts.begin(),verts.end());
    vector<pair<int,int>> be;
    for(auto &kv:edgeCnt) if(kv.second==1) be.push_back(ep[kv.first]);
    if(be.size()<3) return false;
    unordered_map<int, vector<int>> adj; adj.reserve(be.size()*2+1);
    for(auto e:be){adj[e.first].push_back(e.second); adj[e.second].push_back(e.first);} 
    for(auto &kv:adj) if(kv.second.size()!=2) return false;
    loop.clear(); int start=adj.begin()->first, prev=-1, cur=start;
    for(int step=0; step<=(int)adj.size(); ++step){
        loop.push_back(cur); int a=adj[cur][0], b=adj[cur][1]; int nx=(a==prev?b:a); prev=cur; cur=nx; if(cur==start) break;
    }
    if(cur!=start || loop.size()!=adj.size()) return false;
    return true;
}

static Mesh planarAtlasRemesh(const Mesh&in,double eps,double normalCos,double planeTolFactor){
    int F=in.f.size(), V=in.p.size(); if(F<4) return {};
    vector<Vec> fn(F); vector<double> fd(F,0);
    for(int i=0;i<F;i++){auto t=in.f[i]; fn[i]=normed(crossv(in.p[t.b]-in.p[t.a],in.p[t.c]-in.p[t.a])); fd[i]=dotv(fn[i],in.p[t.a]);}
    unordered_map<uint64_t, vector<int>> efaces; efaces.reserve(F*3+7);
    for(int i=0;i<F;i++){auto t=in.f[i]; int a[3]={t.a,t.b,t.c}; for(int j=0;j<3;j++) efaces[ekey(a[j],a[(j+1)%3])].push_back(i);}
    DSU dsu(F);
    double ptol=eps*planeTolFactor;
    for(auto &kv:efaces) if(kv.second.size()==2){int a=kv.second[0], b=kv.second[1]; double co=fabs(dotv(fn[a],fn[b])); if(co>=normalCos){
        double da=fabs(dotv(fn[a], in.p[in.f[b].a])-fd[a]);
        double db=fabs(dotv(fn[b], in.p[in.f[a].a])-fd[b]);
        if(max(da,db)<=ptol) dsu.unite(a,b);
    }}
    unordered_map<int, vector<int>> regs; regs.reserve(F);
    for(int i=0;i<F;i++) regs[dsu.find(i)].push_back(i);
    Mesh out; out.p=in.p;
    vector<Tri> nf;
    vector<char> touchedFace(F,0);
    int useful=0;
    for(auto &rv:regs){
        auto &faces=rv.second;
        if((int)faces.size()<3) { for(int fi:faces) nf.push_back(in.f[fi]); continue; }
        vector<int> loop, allv; if(!extractOneLoop(in,faces,loop,allv)){ for(int fi:faces) nf.push_back(in.f[fi]); continue; }
        if(loop.size()<3 || loop.size()>20000){ for(int fi:faces) nf.push_back(in.f[fi]); continue; }
        Basis bas=makeBasis(in,faces);
        vector<P2> poly; poly.reserve(loop.size()); for(int id:loop) poly.push_back(proj(bas,in.p[id]));
        if(fabs(areaPoly(poly))<1e-16){ for(int fi:faces) nf.push_back(in.f[fi]); continue; }
        vector<int> ears=earClip(poly); if(ears.empty()){ for(int fi:faces) nf.push_back(in.f[fi]); continue; }
        vector<int> sel=loop; unordered_set<int> selected(loop.begin(),loop.end());
        vector<int> interior;
        for(int v:allv) if(!selected.count(v)) interior.push_back(v);
        double cover2=eps*eps*0.86*0.86;
        vector<char> covered(interior.size(),0);
        for(int i=0;i<(int)interior.size();++i){ for(int s:sel) if(dist2v(in.p[interior[i]],in.p[s])<=cover2){covered[i]=1;break;} }
        vector<int> order(interior.size()); iota(order.begin(),order.end(),0);
        sort(order.begin(),order.end(),[&](int a,int b){return dist2v(in.p[interior[a]],bas.o)>dist2v(in.p[interior[b]],bas.o);});
        for(int oi:order) if(!covered[oi]){
            int v=interior[oi]; sel.push_back(v); selected.insert(v);
            for(int j=0;j<(int)interior.size();++j) if(!covered[j] && dist2v(in.p[interior[j]],in.p[v])<=cover2) covered[j]=1;
        }
        vector<P2> q; q.reserve(sel.size()); for(int id:sel) q.push_back(proj(bas,in.p[id]));
        vector<Tri> tris;
        for(int k=0;k<(int)ears.size();k+=3) tris.push_back({ears[k],ears[k+1],ears[k+2]});
        bool ok=true;
        for(int ii=(int)loop.size(); ii<(int)sel.size() && ok; ++ii){
            P2 p=q[ii]; if(!pointInPoly(poly,p)) {ok=false; break;}
            int hit=-1;
            for(int ti=0; ti<(int)tris.size(); ++ti){ if(pointInTri(p,q[tris[ti].a],q[tris[ti].b],q[tris[ti].c])){hit=ti;break;} }
            if(hit<0){ok=false; break;}
            Tri old=tris[hit]; tris.erase(tris.begin()+hit);
            tris.push_back({old.a,old.b,ii}); tris.push_back({old.b,old.c,ii}); tris.push_back({old.c,old.a,ii});
        }
        if(!ok){ for(int fi:faces) nf.push_back(in.f[fi]); continue; }
        int before=nf.size();
        for(auto tr:tris){
            int a=sel[tr.a],b=sel[tr.b],c=sel[tr.c];
            if(a==b||b==c||c==a) continue;
            if(area2_3(in.p[a],in.p[b],in.p[c])<1e-26) continue;
            if(dotv(crossv(in.p[b]-in.p[a],in.p[c]-in.p[a]),bas.n)<0) swap(b,c);
            nf.push_back({a,b,c});
        }
        if((int)nf.size()==before){ for(int fi:faces) nf.push_back(in.f[fi]); }
        else useful += (int)faces.size() - ((int)nf.size()-before);
    }
    if(useful<=0) return {};
    vector<int> used(V,0); for(auto&t:nf) used[t.a]=used[t.b]=used[t.c]=1;
    vector<int> id(V,-1); Mesh m; for(int i=0;i<V;i++) if(used[i]){id[i]=m.p.size(); m.p.push_back(in.p[i]);}
    unordered_set<string> seen; seen.reserve(nf.size()*2+7);
    for(auto t:nf){int a=id[t.a],b=id[t.b],c=id[t.c]; if(a<0||b<0||c<0||a==b||b==c||c==a) continue; int mn=min({a,b,c}), mx=max({a,b,c}), md=a+b+c-mn-mx; string key=to_string(mn)+","+to_string(md)+","+to_string(mx); if(seen.insert(key).second) m.f.push_back({a,b,c});}
    return m;
}

static bool linkOK(int u,int v,const Mesh&base,const vector<int>&mp){
    unordered_set<int> Nu,Nv,cf; int shared=0;
    for(auto&t:base.f){int a=mp[t.a],b=mp[t.b],c=mp[t.c]; if(a==b||b==c||c==a) continue; bool hu=(a==u||b==u||c==u), hv=(a==v||b==v||c==v); if(hu){if(a!=u)Nu.insert(a); if(b!=u)Nu.insert(b); if(c!=u)Nu.insert(c);} if(hv){if(a!=v)Nv.insert(a); if(b!=v)Nv.insert(b); if(c!=v)Nv.insert(c);} if(hu&&hv){shared++; if(a!=u&&a!=v)cf.insert(a); if(b!=u&&b!=v)cf.insert(b); if(c!=u&&c!=v)cf.insert(c);}}
    if(shared!=2 || cf.size()!=2) return false; int inter=0; for(int x:Nu) if(x!=v && Nv.count(x)) inter++; return inter==2;
}
static Mesh compactMap(const Mesh&in, vector<int> mp){
    int n=in.p.size(); for(int i=0;i<n;i++){int x=i; while(mp[x]!=x) x=mp[x]; while(mp[i]!=i){int y=mp[i]; mp[i]=x; i=y;} }
    vector<int> rep(n,0); for(int i=0;i<n;i++) rep[mp[i]]=1;
    vector<int> id(n,-1); Mesh m; for(int i=0;i<n;i++) if(rep[i]){id[i]=m.p.size();m.p.push_back(in.p[i]);}
    unordered_set<string> seen; seen.reserve(in.f.size()*2+7);
    for(auto&t:in.f){int a=mp[t.a],b=mp[t.b],c=mp[t.c]; if(a==b||b==c||c==a) continue; a=id[a];b=id[b];c=id[c]; if(area2_3(m.p[a],m.p[b],m.p[c])<1e-26) continue; int mn=min({a,b,c}),mx=max({a,b,c}),md=a+b+c-mn-mx; string key=to_string(mn)+","+to_string(md)+","+to_string(mx); if(seen.insert(key).second) m.f.push_back({a,b,c});}
    vector<int> used(m.p.size(),0); for(auto&t:m.f) used[t.a]=used[t.b]=used[t.c]=1; vector<int> nid(m.p.size(),-1); vector<Vec> np; for(int i=0;i<(int)m.p.size();i++) if(used[i]){nid[i]=np.size();np.push_back(m.p[i]);} for(auto&t:m.f){t.a=nid[t.a];t.b=nid[t.b];t.c=nid[t.c];} m.p.swap(np); return m;
}
static Mesh sourceCollapse(const Mesh&in,double eps,double factor,int rounds,int mode){
    int n=in.p.size(); vector<int> mp(n); iota(mp.begin(),mp.end(),0); vector<int> alive(n,1); double lim2=eps*eps*factor*factor;
    vector<vector<int>> cluster(n); for(int i=0;i<n;i++) cluster[i].push_back(i);
    for(int r=0;r<rounds;r++){
        unordered_set<uint64_t> es; es.reserve(in.f.size()*3+7); vector<pair<double,pair<int,int>>> cand;
        for(auto&t:in.f){int vv[3]={mp[t.a],mp[t.b],mp[t.c]}; if(vv[0]==vv[1]||vv[1]==vv[2]||vv[2]==vv[0]) continue; for(int i=0;i<3;i++){int a=vv[i],b=vv[(i+1)%3]; if(a==b||!alive[a]||!alive[b]) continue; uint64_t k=ekey(a,b); if(es.insert(k).second){double d=dist2v(in.p[a],in.p[b]); if(d<=lim2)cand.push_back({d,{a,b}});}}}
        if(cand.empty()) break;
        if(mode==0) sort(cand.begin(),cand.end());
        else if(mode==1) sort(cand.begin(),cand.end(),greater<pair<double,pair<int,int>>>());
        else sort(cand.begin(),cand.end(),[&](auto&A,auto&B){return (cluster[A.second.first].size()+cluster[A.second.second].size())>(cluster[B.second.first].size()+cluster[B.second.second].size());});
        int changed=0; vector<char> locked(n,0);
        for(auto &cd:cand){int a=mp[cd.second.first],b=mp[cd.second.second]; if(a==b||!alive[a]||!alive[b]||locked[a]||locked[b]) continue; int keep=a,kill=b; if(cluster[b].size()>cluster[a].size()) swap(keep,kill); bool cover=true; for(int x:cluster[kill]) if(dist2v(in.p[x],in.p[keep])>eps*eps*0.9409){cover=false;break;} if(!cover) continue; if(!linkOK(keep,kill,in,mp)) continue; for(int x:cluster[kill]){mp[x]=keep; cluster[keep].push_back(x);} cluster[kill].clear(); alive[kill]=0; locked[keep]=1; locked[kill]=1; changed++;}
        if(!changed) break;
    }
    return compactMap(in,mp);
}

static Mesh bboxCube(const Mesh&in,double eps){
    if(in.p.empty()) return {}; Vec mn=in.p[0],mx=in.p[0]; for(auto&p:in.p){mn.x=min(mn.x,p.x);mn.y=min(mn.y,p.y);mn.z=min(mn.z,p.z);mx.x=max(mx.x,p.x);mx.y=max(mx.y,p.y);mx.z=max(mx.z,p.z);} vector<Vec> c; for(int ix=0;ix<2;ix++)for(int iy=0;iy<2;iy++)for(int iz=0;iz<2;iz++) c.push_back({ix?mx.x:mn.x,iy?mx.y:mn.y,iz?mx.z:mn.z}); double r2=eps*eps*0.98*0.98; for(auto&p:in.p){double best=1e300; for(auto&q:c) best=min(best,dist2v(p,q)); if(best>r2) return {};} Mesh m; m.p=c; auto add=[&](int a,int b,int c){m.f.push_back({a,b,c});}; add(0,2,3);add(0,3,1);add(4,5,7);add(4,7,6);add(0,1,5);add(0,5,4);add(2,6,7);add(2,7,3);add(0,4,6);add(0,6,2);add(1,3,7);add(1,7,5); return m;
}

static void printMesh(const Mesh&m){
    cout<<m.p.size()<<" "<<m.f.size()<<"\n"; cout.setf(ios::fixed); cout<<setprecision(10);
    for(auto&p:m.p) cout<<"v "<<p.x<<" "<<p.y<<" "<<p.z<<"\n";
    for(auto&t:m.f) cout<<"f "<<t.a+1<<" "<<t.b+1<<" "<<t.c+1<<"\n";
}

int main(){
    ios::sync_with_stdio(false); cin.tie(nullptr);
    int V,F; if(!(cin>>V>>F)) return 0; Mesh in; in.p.reserve(V); in.f.reserve(F); string s;
    for(int i=0;i<V;i++){Vec p; cin>>s>>p.x>>p.y>>p.z; in.p.push_back(p);} 
    for(int i=0;i<F;i++){Tri t; cin>>s>>t.a>>t.b>>t.c; --t.a;--t.b;--t.c; in.f.push_back(t);} 
    if(V==0){cout<<"0 0\n";return 0;}
    Vec mn=in.p[0],mx=in.p[0]; for(auto&p:in.p){mn.x=min(mn.x,p.x);mn.y=min(mn.y,p.y);mn.z=min(mn.z,p.z);mx.x=max(mx.x,p.x);mx.y=max(mx.y,p.y);mx.z=max(mx.z,p.z);} double diag=sqrt(dist2v(mn,mx)); double eps=diag*0.05;
    Mesh best=in;
    if(!validateClosed(in)){ printMesh(in); return 0; }
    auto consider=[&](const Mesh&cand){ if(cand.p.empty()) return; if(cand.p.size()>=best.p.size()) return; if(!validateClosed(cand)) return; if(!vertexHausdorff(in.p,cand,eps)) return; best=cand; };
    consider(bboxCube(in,eps));
    vector<pair<double,double>> pars={{0.99995,0.04},{0.9998,0.06},{0.9993,0.08},{0.998,0.10},{0.995,0.12},{0.990,0.15}};
    for(auto pr:pars) consider(planarAtlasRemesh(in,eps,pr.first,pr.second));
    vector<double> fac={0.92,0.80,0.68,0.56,0.44,0.32,0.22};
    for(int mode=0;mode<3;mode++) for(double f:fac) consider(sourceCollapse(in,eps,f,10,mode));
    printMesh(best);
    return 0;
}