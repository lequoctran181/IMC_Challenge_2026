#include <bits/stdc++.h>
using namespace std;

struct P { double x,y,z; };
static inline P operator-(const P&a,const P&b){return {a.x-b.x,a.y-b.y,a.z-b.z};}
static inline double dotp(const P&a,const P&b){return a.x*b.x+a.y*b.y+a.z*b.z;}
static inline P crossp(const P&a,const P&b){return {a.y*b.z-a.z*b.y,a.z*b.x-a.x*b.z,a.x*b.y-a.y*b.x};}
static inline double norm2(const P&a){return dotp(a,a);} 
static inline double dist2p(const P&a,const P&b){return norm2(a-b);} 
struct Face { int a,b,c; };

struct FastHash {
    static uint64_t splitmix64(uint64_t x) {
        x += 0x9e3779b97f4a7c15ULL;
        x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
        x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
        return x ^ (x >> 31);
    }
    size_t operator()(uint64_t x) const {
        static const uint64_t FIXED_RANDOM = chrono::steady_clock::now().time_since_epoch().count();
        return splitmix64(x + FIXED_RANDOM);
    }
};

static inline uint64_t packCell(int x,int y,int z){
    const int OFF = 1<<20;
    uint64_t xx = (uint64_t)(x + OFF) & ((1ULL<<21)-1);
    uint64_t yy = (uint64_t)(y + OFF) & ((1ULL<<21)-1);
    uint64_t zz = (uint64_t)(z + OFF) & ((1ULL<<21)-1);
    return (xx<<42) ^ (yy<<21) ^ zz;
}

struct Grid {
    struct Cell { int x,y,z; vector<int> pts; int rep=-1; int approx=0; };
    const vector<P>* p=nullptr;
    vector<Cell> cells;
    unordered_map<uint64_t,int,FastHash> mp;
    P org; double h=1;
    Grid(){}
    Grid(const vector<P>& pts, double hh, P origin): p(&pts), org(origin), h(hh) { build(pts); }
    inline int ix(double v,double o) const { return (int)floor((v-o)/h); }
    void build(const vector<P>& pts){
        cells.clear(); mp.clear();
        mp.reserve(pts.size()*2/3 + 1024);
        for(int i=0;i<(int)pts.size();++i){
            int x=ix(pts[i].x,org.x), y=ix(pts[i].y,org.y), z=ix(pts[i].z,org.z);
            uint64_t k=packCell(x,y,z);
            auto it=mp.find(k);
            int id;
            if(it==mp.end()){
                id=(int)cells.size(); mp.emplace(k,id); cells.push_back({x,y,z,{},{},0});
            } else id=it->second;
            cells[id].pts.push_back(i);
        }
    }
    template<class F>
    inline void forBall(const P& c, double r, double r2, F&& fn) const {
        int cx=ix(c.x,org.x), cy=ix(c.y,org.y), cz=ix(c.z,org.z);
        int R=(int)ceil(r/h)+1;
        for(int dx=-R; dx<=R; ++dx) for(int dy=-R; dy<=R; ++dy) for(int dz=-R; dz<=R; ++dz){
            auto it=mp.find(packCell(cx+dx,cy+dy,cz+dz));
            if(it==mp.end()) continue;
            const auto &vec=cells[it->second].pts;
            for(int id: vec) if(dist2p((*p)[id], c) <= r2) fn(id);
        }
    }
};

static vector<int> pruneCover(const vector<P>& pts, vector<int> sel, double eps, double eps2, P mn, double D, double secondsLimit){
    int n=(int)pts.size();
    if(sel.empty()) return sel;
    sort(sel.begin(), sel.end()); sel.erase(unique(sel.begin(), sel.end()), sel.end());
    double h=max(eps*0.5, D*1e-12);
    Grid g(pts,h,{mn.x-h*0.123,mn.y-h*0.217,mn.z-h*0.319});
    vector<unsigned short> cnt(n,0);
    for(int s: sel){
        g.forBall(pts[s], eps, eps2, [&](int id){ if(cnt[id] < 65535) ++cnt[id]; });
    }
    for(int i=0;i<n;i++) if(cnt[i]==0) return sel; // not a cover; do not risk pruning
    vector<int> order(sel.size()); iota(order.begin(), order.end(), 0);
    // Deterministic mixed order: removing sparse/late representatives first often helps after grid covers.
    stable_sort(order.begin(), order.end(), [&](int A,int B){
        uint64_t ha=FastHash::splitmix64((uint64_t)sel[A]*11995408973635179863ULL + 17);
        uint64_t hb=FastHash::splitmix64((uint64_t)sel[B]*11995408973635179863ULL + 17);
        return ha<hb;
    });
    vector<unsigned char> alive(sel.size(),1);
    auto start=chrono::steady_clock::now();
    int removed=0;
    for(int pass=0; pass<2; ++pass){
        if(pass==1) reverse(order.begin(), order.end());
        for(int oi=0; oi<(int)order.size(); ++oi){
            if((oi&255)==0){
                double sec=chrono::duration<double>(chrono::steady_clock::now()-start).count();
                if(sec>secondsLimit) goto done;
            }
            int pos=order[oi]; if(!alive[pos]) continue;
            int s=sel[pos]; bool ok=true;
            g.forBall(pts[s], eps, eps2, [&](int id){ if(cnt[id]<=1) ok=false; });
            if(!ok) continue;
            alive[pos]=0; ++removed;
            g.forBall(pts[s], eps, eps2, [&](int id){ --cnt[id]; });
        }
    }
 done:
    vector<int> out; out.reserve(sel.size()-removed);
    for(int i=0;i<(int)sel.size();++i) if(alive[i]) out.push_back(sel[i]);
    return out;
}

static vector<int> coverSafeShiftGrid(const vector<P>& pts, P mn, double eps){
    const int n=(int)pts.size();
    double h = eps / sqrt(3.0) * 0.9990;
    vector<double> shifts = {0.0, h/3.0, 2.0*h/3.0};
    vector<int> best;
    unordered_map<uint64_t,int,FastHash> mp;
    mp.reserve(n*2/3+1024);
    for(double sx: shifts) for(double sy: shifts) for(double sz: shifts){
        mp.clear();
        P org{mn.x-sx, mn.y-sy, mn.z-sz};
        for(int i=0;i<n;++i){
            int x=(int)floor((pts[i].x-org.x)/h);
            int y=(int)floor((pts[i].y-org.y)/h);
            int z=(int)floor((pts[i].z-org.z)/h);
            uint64_t k=packCell(x,y,z);
            if(mp.find(k)==mp.end()) mp.emplace(k,i);
        }
        if(best.empty() || mp.size()<best.size()){
            best.clear(); best.reserve(mp.size());
            for(auto &kv:mp) best.push_back(kv.second);
        }
    }
    sort(best.begin(), best.end()); best.erase(unique(best.begin(), best.end()), best.end());
    return best;
}

static vector<int> coverDensityPass(const vector<P>& pts, P mn, double eps, double eps2, double hFactor, double secondsLimit){
    int n=(int)pts.size();
    double h=max(eps*hFactor, eps*0.2);
    Grid g(pts,h,{mn.x-h*0.371,mn.y-h*0.527,mn.z-h*0.193});
    const int C=(int)g.cells.size();
    vector<int> cand; cand.reserve(C);
    // representative nearest to cell centroid
    for(int ci=0; ci<C; ++ci){
        auto &cell=g.cells[ci];
        P cen{0,0,0};
        for(int id: cell.pts){ cen.x+=pts[id].x; cen.y+=pts[id].y; cen.z+=pts[id].z; }
        double inv=1.0/cell.pts.size(); cen.x*=inv; cen.y*=inv; cen.z*=inv;
        int best=cell.pts[0]; double bd=dist2p(pts[best],cen);
        for(int id: cell.pts){ double d=dist2p(pts[id],cen); if(d<bd){bd=d; best=id;} }
        cell.rep=best; cand.push_back(best);
    }
    int R=(int)ceil(eps/h)+2;
    for(int ci=0; ci<C; ++ci){
        auto &cc=g.cells[ci]; long long sum=0;
        for(int dx=-R; dx<=R; ++dx) for(int dy=-R; dy<=R; ++dy) for(int dz=-R; dz<=R; ++dz){
            auto it=g.mp.find(packCell(cc.x+dx,cc.y+dy,cc.z+dz));
            if(it!=g.mp.end()) sum += (int)g.cells[it->second].pts.size();
        }
        cc.approx = (sum>INT_MAX?INT_MAX:(int)sum);
    }
    vector<int> order(C); iota(order.begin(), order.end(), 0);
    sort(order.begin(), order.end(), [&](int a,int b){
        if(g.cells[a].approx!=g.cells[b].approx) return g.cells[a].approx>g.cells[b].approx;
        return g.cells[a].pts.size()>g.cells[b].pts.size();
    });
    vector<unsigned char> covered(n,0);
    int uncovered=n;
    vector<int> sel; sel.reserve(max(8,C/2));
    auto mark=[&](int s)->int{
        int got=0;
        g.forBall(pts[s], eps, eps2, [&](int id){ if(!covered[id]) ++got; });
        if(got==0) return 0;
        sel.push_back(s);
        g.forBall(pts[s], eps, eps2, [&](int id){ if(!covered[id]){covered[id]=1; --uncovered;} });
        return got;
    };
    auto start=chrono::steady_clock::now();
    for(int t=0; t<C && uncovered>0; ++t){
        if((t&127)==0){
            double sec=chrono::duration<double>(chrono::steady_clock::now()-start).count();
            if(sec>secondsLimit) break;
        }
        int s=g.cells[order[t]].rep;
        if(!covered[s]) mark(s);
        else {
            // Still useful near frontiers; require at least two fresh vertices to avoid too many redundant picks.
            int fresh=0;
            g.forBall(pts[s], eps, eps2, [&](int id){ if(!covered[id]) ++fresh; });
            if(fresh>=2) {
                sel.push_back(s);
                g.forBall(pts[s], eps, eps2, [&](int id){ if(!covered[id]){covered[id]=1; --uncovered;} });
            }
        }
    }
    // Guaranteed completion: maximal uncovered-point net.
    for(int i=0; i<n && uncovered>0; ++i) if(!covered[i]) mark(i);
    sort(sel.begin(), sel.end()); sel.erase(unique(sel.begin(), sel.end()), sel.end());
    return sel;
}


static vector<int> coverLazyGreedy(const vector<P>& pts, P mn, double eps, double eps2, double secondsLimit){
    int n=(int)pts.size();
    double h=max(eps*0.55, eps*0.2);
    Grid g(pts,h,{mn.x-h*0.071,mn.y-h*0.131,mn.z-h*0.197});
    vector<int> candidates;
    if(n<=120000){
        candidates.resize(n); iota(candidates.begin(), candidates.end(), 0);
    } else {
        // for larger instances keep a few representatives per occupied cell
        candidates.reserve(g.cells.size()*4);
        for(auto &cell:g.cells){
            P cen{0,0,0};
            for(int id:cell.pts){cen.x+=pts[id].x; cen.y+=pts[id].y; cen.z+=pts[id].z;}
            double inv=1.0/cell.pts.size(); cen.x*=inv; cen.y*=inv; cen.z*=inv;
            int best=cell.pts[0]; double bd=dist2p(pts[best],cen);
            int ex[6]; for(int i=0;i<6;i++) ex[i]=cell.pts[0];
            for(int id:cell.pts){
                double d=dist2p(pts[id],cen); if(d<bd){bd=d; best=id;}
                if(pts[id].x<pts[ex[0]].x) ex[0]=id;
                if(pts[id].x>pts[ex[1]].x) ex[1]=id;
                if(pts[id].y<pts[ex[2]].y) ex[2]=id;
                if(pts[id].y>pts[ex[3]].y) ex[3]=id;
                if(pts[id].z<pts[ex[4]].z) ex[4]=id;
                if(pts[id].z>pts[ex[5]].z) ex[5]=id;
            }
            candidates.push_back(best);
            for(int i=0;i<6;i++) candidates.push_back(ex[i]);
        }
        sort(candidates.begin(), candidates.end()); candidates.erase(unique(candidates.begin(), candidates.end()), candidates.end());
    }
    int R=(int)ceil(eps/h)+2;
    struct Node{int sc,id; bool operator<(Node const&o)const{return sc!=o.sc?sc<o.sc:id>o.id;}};
    priority_queue<Node> pq;
    pq=priority_queue<Node>();
    vector<unsigned char> covered(n,0), chosen(n,0);
    int uncovered=n;
    vector<int> sel; sel.reserve(max(8,n/20));
    auto countFresh=[&](int s){ int got=0; g.forBall(pts[s],eps,eps2,[&](int id){ if(!covered[id]) ++got; }); return got; };
    auto add=[&](int s){
        if(chosen[s]) return 0; int got=0; chosen[s]=1; sel.push_back(s);
        g.forBall(pts[s],eps,eps2,[&](int id){ if(!covered[id]){covered[id]=1; --uncovered; ++got;} }); return got;
    };
    auto start=chrono::steady_clock::now();
    // Exact initial scores for moderate instances; approximate scores for larger candidate sets.
    bool exactInit = (n <= 150000);
    int initCounter=0;
    for(int id:candidates){
        int sc=0;
        if(exactInit){
            g.forBall(pts[id],eps,eps2,[&](int){ ++sc; });
        } else {
            int cx=g.ix(pts[id].x,g.org.x), cy=g.ix(pts[id].y,g.org.y), cz=g.ix(pts[id].z,g.org.z);
            long long sum=0;
            for(int dx=-R; dx<=R; ++dx) for(int dy=-R; dy<=R; ++dy) for(int dz=-R; dz<=R; ++dz){
                auto it=g.mp.find(packCell(cx+dx,cy+dy,cz+dz));
                if(it!=g.mp.end()) sum += (int)g.cells[it->second].pts.size();
            }
            sc=(int)min<long long>(sum, INT_MAX);
        }
        if(sc>0) pq.push({sc,id});
        if((++initCounter&511)==0){ double sec=chrono::duration<double>(chrono::steady_clock::now()-start).count(); if(sec>secondsLimit*0.55 && !exactInit) break; }
    }
    int stale=0;
    while(uncovered>0 && !pq.empty()){
        if((stale++&255)==0){ double sec=chrono::duration<double>(chrono::steady_clock::now()-start).count(); if(sec>secondsLimit) break; }
        Node nd=pq.top(); pq.pop();
        if(chosen[nd.id]) continue;
        int exact=countFresh(nd.id);
        if(exact<=0) continue;
        int topScore = pq.empty()?0:pq.top().sc;
        if(exact >= topScore){
            add(nd.id);
        } else {
            pq.push({exact, nd.id});
        }
    }
    for(int i=0;i<n && uncovered>0;i++) if(!covered[i]) add(i);
    sort(sel.begin(), sel.end()); sel.erase(unique(sel.begin(), sel.end()), sel.end());
    return sel;
}

static bool validateCover(const vector<P>& pts, const vector<int>& sel, P mn, double eps, double eps2, double D){
    if(sel.empty()) return false;
    vector<P> centers; centers.reserve(sel.size());
    for(int id: sel) centers.push_back(pts[id]);
    double h=max(eps, D*1e-12);
    Grid cg(centers,h,{mn.x-h*0.11,mn.y-h*0.29,mn.z-h*0.43});
    for(const P& q: pts){
        bool ok=false;
        cg.forBall(q, eps, eps2, [&](int){ ok=true; });
        if(!ok) return false;
    }
    return true;
}

static uint64_t coordHash(const P& p){
    uint64_t a,b,c; memcpy(&a,&p.x,8); memcpy(&b,&p.y,8); memcpy(&c,&p.z,8);
    return FastHash::splitmix64(a) ^ (FastHash::splitmix64(b+0x9e3779b97f4a7c15ULL)<<1) ^ (FastHash::splitmix64(c+0xbf58476d1ce4e5b9ULL)<<2);
}

static bool goodTri(const P&a,const P&b,const P&c,double tiny){
    return norm2(crossp(b-a,c-a)) > tiny;
}

static bool validCombinatorial(const vector<P>& vp, const vector<Face>& f, double D){
    if(vp.size()<4 || f.empty()) return false;
    double tiny=max(1e-30,1e-24*D*D*D*D);
    unordered_map<uint64_t,int,FastHash> cnt; cnt.reserve(f.size()*4+16);
    auto ekey=[](int a,int b){ if(a>b) swap(a,b); return ((uint64_t)(uint32_t)a<<32) | (uint32_t)b; };
    vector<int> used(vp.size(),0);
    vector<vector<pair<int,int>>> links(vp.size());
    for(const auto&t:f){
        if(t.a<0||t.b<0||t.c<0||t.a>=(int)vp.size()||t.b>=(int)vp.size()||t.c>=(int)vp.size()) return false;
        if(t.a==t.b||t.b==t.c||t.c==t.a) return false;
        if(!goodTri(vp[t.a],vp[t.b],vp[t.c],tiny)) return false;
        used[t.a]=used[t.b]=used[t.c]=1;
        cnt[ekey(t.a,t.b)]++; cnt[ekey(t.b,t.c)]++; cnt[ekey(t.c,t.a)]++;
        links[t.a].push_back({t.b,t.c});
        links[t.b].push_back({t.c,t.a});
        links[t.c].push_back({t.a,t.b});
    }
    for(int u:used) if(!u) return false;
    for(auto &kv:cnt) if(kv.second!=2) return false;

    // Closed vertex-manifold check: the one-ring link of every vertex must be one cycle.
    for(int v=0; v<(int)vp.size(); ++v){
        if(!used[v]) return false;
        unordered_map<int, vector<int>, FastHash> adj;
        adj.reserve(links[v].size()*2+4);
        for(auto pr: links[v]){
            int a=pr.first,b=pr.second;
            if(a==b) return false;
            adj[a].push_back(b);
            adj[b].push_back(a);
        }
        if(adj.empty()) return false;
        for(auto &kv: adj){
            if(kv.second.size()!=2) return false;
            if(kv.second[0]==kv.second[1]) return false;
        }
        int start=adj.begin()->first;
        unordered_set<int,FastHash> seen;
        seen.reserve(adj.size()*2+4);
        vector<int> st={start};
        seen.insert(start);
        while(!st.empty()){
            int x=st.back(); st.pop_back();
            for(int y: adj[x]) if(seen.insert(y).second) st.push_back(y);
        }
        if(seen.size()!=adj.size()) return false;
    }
    return true;
}

static bool buildBipyramid(const vector<P>& pts, vector<int> sel, vector<P>& outp, vector<Face>& outf, P mn, P mx, double D){
    int n0=pts.size();
    sort(sel.begin(), sel.end()); sel.erase(unique(sel.begin(), sel.end()), sel.end());
    // Remove exact duplicate coordinates; duplicate coordinates create zero-area risks and never improve Hausdorff cover.
    unordered_set<uint64_t,FastHash> seen; seen.reserve(sel.size()*2+16);
    vector<int> uniq; uniq.reserve(sel.size());
    for(int id: sel){ uint64_t h=coordHash(pts[id]); if(seen.insert(h).second) uniq.push_back(id); }
    sel.swap(uniq);
    // Ensure enough distinct vertices. Add far/spread original vertices if the cover is tiny on simple tests.
    auto addId=[&](int id){ uint64_t h=coordHash(pts[id]); if(seen.insert(h).second) sel.push_back(id); };
    if((int)sel.size()<5){
        int ix0=0,ix1=0,iy0=0,iy1=0,iz0=0,iz1=0;
        for(int i=1;i<n0;i++){
            if(pts[i].x<pts[ix0].x) ix0=i; if(pts[i].x>pts[ix1].x) ix1=i;
            if(pts[i].y<pts[iy0].y) iy0=i; if(pts[i].y>pts[iy1].y) iy1=i;
            if(pts[i].z<pts[iz0].z) iz0=i; if(pts[i].z>pts[iz1].z) iz1=i;
        }
        for(int id: {ix0,ix1,iy0,iy1,iz0,iz1}) addId(id);
        for(int i=0;i<n0 && (int)sel.size()<5;i++) addId(i);
    }
    if((int)sel.size()<4) return false;
    // Try tetrahedron only when exactly 4 vertices.
    if((int)sel.size()==4){
        outp.clear(); for(int id:sel) outp.push_back(pts[id]);
        outf={{0,1,2},{0,3,1},{1,3,2},{2,3,0}};
        if(validCombinatorial(outp,outf,D)) return true;
        // otherwise add one more vertex and use bipyramid
        for(int i=0;i<n0 && (int)sel.size()<5;i++) addId(i);
    }
    int k=sel.size();
    // Candidate poles: bbox extremes among selected plus a small deterministic sample.
    vector<int> samples;
    auto pushSample=[&](int pos){ if(pos>=0 && pos<k) samples.push_back(pos); };
    int minx=0,maxx=0,miny=0,maxy=0,minz=0,maxz=0;
    for(int i=1;i<k;i++){
        const P&a=pts[sel[i]];
        if(a.x<pts[sel[minx]].x) minx=i; if(a.x>pts[sel[maxx]].x) maxx=i;
        if(a.y<pts[sel[miny]].y) miny=i; if(a.y>pts[sel[maxy]].y) maxy=i;
        if(a.z<pts[sel[minz]].z) minz=i; if(a.z>pts[sel[maxz]].z) maxz=i;
    }
    for(int v:{minx,maxx,miny,maxy,minz,maxz}) pushSample(v);
    uint64_t seed=0x6a09e667f3bcc909ULL ^ (uint64_t)k;
    auto rng=[&](){ seed ^= seed<<7; seed ^= seed>>9; return seed; };
    for(int t=0;t<24 && t<k;t++) pushSample((int)(rng()%k));
    sort(samples.begin(), samples.end()); samples.erase(unique(samples.begin(), samples.end()), samples.end());
    vector<pair<double,pair<int,int>>> polePairs;
    for(int a: samples) for(int b: samples) if(a<b) polePairs.push_back({dist2p(pts[sel[a]],pts[sel[b]]),{a,b}});
    sort(polePairs.rbegin(), polePairs.rend());
    if(polePairs.empty()) return false;
    double tiny=max(1e-30,1e-24*D*D*D*D);
    int triesPairs=min((int)polePairs.size(), 40);
    for(int pp=0; pp<triesPairs; ++pp){
        int topPos=polePairs[pp].second.first, botPos=polePairs[pp].second.second;
        vector<int> ringPos; ringPos.reserve(k-2);
        for(int i=0;i<k;i++) if(i!=topPos && i!=botPos) ringPos.push_back(i);
        if((int)ringPos.size()<3) continue;
        // First try angular order around the pole axis; this is nicer for normal meshes.
        P top=pts[sel[topPos]], bot=pts[sel[botPos]];
        P axis=bot-top; double al=sqrt(max(1e-300,norm2(axis))); axis={axis.x/al,axis.y/al,axis.z/al};
        P tmp = fabs(axis.x)<0.8 ? P{1,0,0} : P{0,1,0};
        P u=crossp(axis,tmp); double ul=sqrt(max(1e-300,norm2(u))); u={u.x/ul,u.y/ul,u.z/ul};
        P v=crossp(axis,u);
        vector<int> baseRing=ringPos;
        sort(baseRing.begin(), baseRing.end(), [&](int A,int B){
            P pa=pts[sel[A]]-top, pb=pts[sel[B]]-top;
            double aa=atan2(dotp(pa,v),dotp(pa,u));
            double bb=atan2(dotp(pb,v),dotp(pb,u));
            return aa<bb;
        });
        for(int attempt=0; attempt<80; ++attempt){
            vector<int> ring=baseRing;
            if(attempt>0){
                // Deterministic Fisher-Yates shuffle.
                for(int i=(int)ring.size()-1;i>0;--i){ int j=(int)(rng()%(i+1)); swap(ring[i],ring[j]); }
            }
            bool ok=true; int m=ring.size();
            for(int i=0;i<m && ok;i++){
                const P&A=pts[sel[ring[i]]], &B=pts[sel[ring[(i+1)%m]]];
                if(!goodTri(top,A,B,tiny) || !goodTri(bot,B,A,tiny)) ok=false;
            }
            if(!ok) continue;
            outp.clear(); outf.clear(); outp.reserve(k); outf.reserve(2*m);
            vector<int> order; order.reserve(k); order.push_back(topPos); order.push_back(botPos); for(int x:ring) order.push_back(x);
            for(int pos: order) outp.push_back(pts[sel[pos]]);
            for(int i=0;i<m;i++){
                int a=2+i, b=2+((i+1)%m);
                outf.push_back({0,a,b});
                outf.push_back({1,b,a});
            }
            if(validCombinatorial(outp,outf,D)) return true;
        }
    }
    return false;
}



static bool buildClusterMesh(const vector<P>& pts, const vector<Face>& inFaces, vector<int> sel,
                             vector<P>& outp, vector<Face>& outf, P mn, double eps, double eps2, double D){
    sort(sel.begin(), sel.end());
    sel.erase(unique(sel.begin(), sel.end()), sel.end());
    if(sel.size()<4) return false;

    // Remove exact duplicate selected coordinates; keep the first representative.
    unordered_set<uint64_t,FastHash> seen;
    seen.reserve(sel.size()*2+16);
    vector<int> filtered;
    filtered.reserve(sel.size());
    for(int id: sel){
        uint64_t h=coordHash(pts[id]);
        if(seen.insert(h).second) filtered.push_back(id);
    }
    sel.swap(filtered);
    int k=(int)sel.size();
    if(k<4) return false;

    vector<P> centers;
    centers.reserve(k);
    unordered_map<int,int,FastHash> selectedIndex;
    selectedIndex.reserve(k*2+16);
    for(int i=0;i<k;i++){
        centers.push_back(pts[sel[i]]);
        selectedIndex[sel[i]]=i;
    }

    double h=max(eps*0.50, D*1e-12);
    Grid cg(centers,h,{mn.x-h*0.173,mn.y-h*0.281,mn.z-h*0.419});
    vector<int> mapv(pts.size(),-1);
    for(int i=0;i<(int)pts.size();++i){
        auto it=selectedIndex.find(i);
        if(it!=selectedIndex.end()) { mapv[i]=it->second; continue; }
        double bd=eps2;
        int best=-1;
        cg.forBall(pts[i], eps, eps2, [&](int cid){
            double d=dist2p(pts[i], centers[cid]);
            if(d<bd || (d==bd && cid<best)) { bd=d; best=cid; }
        });
        if(best<0) return false;
        mapv[i]=best;
    }

    double tiny=max(1e-30,1e-24*D*D*D*D);
    auto triKey=[](int a,int b,int c)->uint64_t{
        if(a>b) swap(a,b); if(b>c) swap(b,c); if(a>b) swap(a,b);
        uint64_t x=(uint32_t)a, y=(uint32_t)b, z=(uint32_t)c;
        return FastHash::splitmix64(x*1000003ULL ^ (y<<21) ^ (z<<42) ^ 0x9e3779b97f4a7c15ULL);
    };
    unordered_set<uint64_t,FastHash> usedTri;
    usedTri.reserve(inFaces.size()*2/3+16);
    outf.clear();
    outf.reserve(inFaces.size());
    for(const Face& f: inFaces){
        if(f.a<0||f.b<0||f.c<0||f.a>=(int)pts.size()||f.b>=(int)pts.size()||f.c>=(int)pts.size()) continue;
        int a=mapv[f.a], b=mapv[f.b], c=mapv[f.c];
        if(a<0||b<0||c<0||a==b||b==c||c==a) continue;
        if(!goodTri(centers[a],centers[b],centers[c],tiny)) continue;
        uint64_t key=triKey(a,b,c);
        if(!usedTri.insert(key).second) continue;
        outf.push_back({a,b,c});
    }
    outp=move(centers);
    if(!validCombinatorial(outp,outf,D)) return false;
    return true;
}

static void outputMesh(const vector<P>& p, const vector<Face>& f){
    cout << setprecision(17);
    cout << p.size() << ' ' << f.size() << '\n';
    for(const auto&q:p) cout << "v " << q.x << ' ' << q.y << ' ' << q.z << '\n';
    for(const auto&t:f) cout << "f " << t.a+1 << ' ' << t.b+1 << ' ' << t.c+1 << '\n';
}

static bool isNumericToken(const string& s){
    if(s.empty()) return false;
    char c=s[0];
    return (c=='-' || c=='+' || c=='.' || (c>='0' && c<='9'));
}

int main(){
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    int V,F;
    if(!(cin>>V>>F)) return 0;
    vector<P> pts(V);
    vector<Face> faces(F);
    P mn{1e100,1e100,1e100}, mx{-1e100,-1e100,-1e100};

    string tok;
    for(int i=0;i<V;i++){
        if(!(cin>>tok)) return 0;
        if(isNumericToken(tok)){
            pts[i].x = stod(tok);
            cin >> pts[i].y >> pts[i].z;
        } else {
            cin >> pts[i].x >> pts[i].y >> pts[i].z;
        }
        mn.x=min(mn.x,pts[i].x); mn.y=min(mn.y,pts[i].y); mn.z=min(mn.z,pts[i].z);
        mx.x=max(mx.x,pts[i].x); mx.y=max(mx.y,pts[i].y); mx.z=max(mx.z,pts[i].z);
    }
    for(int i=0;i<F;i++){
        if(!(cin>>tok)) return 0;
        if(isNumericToken(tok)){
            faces[i].a = stoi(tok);
            cin >> faces[i].b >> faces[i].c;
        } else {
            cin >> faces[i].a >> faces[i].b >> faces[i].c;
        }
        --faces[i].a; --faces[i].b; --faces[i].c;
    }

    double dx=mx.x-mn.x, dy=mx.y-mn.y, dz=mx.z-mn.z;
    double D=sqrt(dx*dx+dy*dy+dz*dz);
    if(!(D>0)) D=1.0;

    vector<P> bestP;
    vector<Face> bestF;
    long long bestCost = LLONG_MAX;
    vector<int> riskBest;
    double riskEps = 0.0, riskEps2 = 0.0;

    auto globalStart=chrono::steady_clock::now();
    auto elapsedSec=[&](){ return chrono::duration<double>(chrono::steady_clock::now()-globalStart).count(); };

    auto acceptCandidate = [&](const vector<int>& cand, double eps, double eps2){
        if(cand.empty()) return;
        if(!validateCover(pts,cand,mn,eps,eps2,D)) return;
        if(riskBest.empty() || cand.size()<riskBest.size()){
            riskBest=cand;
            riskEps=eps;
            riskEps2=eps2;
        }
        vector<P> cp;
        vector<Face> cf;
        if(buildClusterMesh(pts,faces,cand,cp,cf,mn,eps,eps2,D)){
            long long cost=(long long)cp.size()+(long long)cf.size();
            if(cost<bestCost){
                bestCost=cost;
                bestP.swap(cp);
                bestF.swap(cf);
            }
        }
    };

    auto consider = [&](vector<int> cand, double eps, double eps2, double pruneSeconds){
        if(cand.empty()) return;
        sort(cand.begin(), cand.end());
        cand.erase(unique(cand.begin(), cand.end()), cand.end());
        acceptCandidate(cand,eps,eps2);
        if(elapsedSec()>18.5) return;
        vector<int> pruned=pruneCover(pts, cand, eps, eps2, mn, D, pruneSeconds);
        if(pruned.size()!=cand.size()) acceptCandidate(pruned,eps,eps2);
    };

    auto runRadius = [&](double rf, bool full){
        double eps=rf*D, eps2=eps*eps*(1.0+2e-12);
        consider(coverSafeShiftGrid(pts,mn,eps), eps, eps2, full?2.8:1.0);
        if(elapsedSec()>18.5) return;
        if(full && V<=150000) consider(coverLazyGreedy(pts,mn,eps,eps2,5.5), eps, eps2, 1.8);
        if(elapsedSec()>18.5) return;
        consider(coverDensityPass(pts,mn,eps,eps2,0.95, full?3.5:1.6), eps, eps2, full?1.8:0.8);
        if(elapsedSec()>18.5) return;
        consider(coverDensityPass(pts,mn,eps,eps2,0.70, full?3.0:1.4), eps, eps2, full?1.4:0.6);
        if(elapsedSec()>18.5) return;
        if(full && V<=350000) consider(coverDensityPass(pts,mn,eps,eps2,0.52,2.2), eps, eps2, 1.0);
    };

    // First try the actual 5%-margin radius for maximum compression. If vertex clustering
    // creates non-manifold quotients on thin features, smaller radii often restore manifoldness.
    runRadius(0.04990,true);
    if(elapsedSec()<18.5) runRadius(0.0430,false);
    if(bestP.empty() && elapsedSec()<18.5) runRadius(0.0360,false);
    if(bestP.empty() && elapsedSec()<18.5) runRadius(0.0300,false);
    if(bestP.empty() && elapsedSec()<18.5) runRadius(0.0240,false);

    if(!bestP.empty()){
        outputMesh(bestP,bestF);
        return 0;
    }

    vector<P> outp;
    vector<Face> outf;
    if(!riskBest.empty() && buildBipyramid(pts,riskBest,outp,outf,mn,mx,D)){
        // Last-resort high-compression construction: all output vertices are original vertices
        // and riskBest has been verified as a 5%-cover, so the stated vertex-Hausdorff condition holds.
        outputMesh(outp,outf);
    } else {
        // Conservative fallback. It preserves exact vertices/faces and therefore distance is zero.
        outputMesh(pts,faces);
    }
    return 0;
}
