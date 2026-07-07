#include <bits/stdc++.h>
using namespace std;

struct P{double x,y,z;};
static inline P operator+(const P&a,const P&b){return {a.x+b.x,a.y+b.y,a.z+b.z};}
static inline P operator-(const P&a,const P&b){return {a.x-b.x,a.y-b.y,a.z-b.z};}
static inline P operator*(const P&a,double s){return {a.x*s,a.y*s,a.z*s};}
static inline double dotp(const P&a,const P&b){return a.x*b.x+a.y*b.y+a.z*b.z;}
static inline P crossp(const P&a,const P&b){return {a.y*b.z-a.z*b.y,a.z*b.x-a.x*b.z,a.x*b.y-a.y*b.x};}
static inline double norm2(const P&a){return dotp(a,a);} 
static inline double norm(const P&a){return sqrt(max(0.0,norm2(a)));}

struct Tri{int a,b,c;};
struct EdgeInfo{int cnt=0; int opp1=-1, opp2=-1;};
static inline long long key2(int a,int b){ if(a>b) swap(a,b); return ( (long long)a<<32 ) ^ (unsigned int)b; }

static double triArea2(const vector<P>&v,const Tri&t){ return norm(crossp(v[t.b]-v[t.a], v[t.c]-v[t.a])); }
static P triN(const vector<P>&v,const Tri&t){ P n=crossp(v[t.b]-v[t.a], v[t.c]-v[t.a]); double l=norm(n); if(l==0) return {0,0,0}; return n*(1.0/l); }

static bool validClosed(const vector<P>&v,const vector<Tri>&f){
    unordered_map<long long,int> ec; ec.reserve(f.size()*3+10);
    double bb=0; P mn={1e100,1e100,1e100}, mx={-1e100,-1e100,-1e100};
    for(auto&p:v){mn.x=min(mn.x,p.x);mn.y=min(mn.y,p.y);mn.z=min(mn.z,p.z);mx.x=max(mx.x,p.x);mx.y=max(mx.y,p.y);mx.z=max(mx.z,p.z);} bb=norm(mx-mn);
    double eps=max(1e-18,bb*1e-12);
    for(auto&t:f){
        if(t.a<0||t.b<0||t.c<0||t.a>=(int)v.size()||t.b>=(int)v.size()||t.c>=(int)v.size()) return false;
        if(t.a==t.b||t.b==t.c||t.c==t.a) return false;
        if(triArea2(v,t)<=eps) return false;
        ec[key2(t.a,t.b)]++; ec[key2(t.b,t.c)]++; ec[key2(t.c,t.a)]++;
    }
    for(auto &kv:ec) if(kv.second!=2) return false;
    return true;
}

int main(){
    ios::sync_with_stdio(false); cin.tie(nullptr);
    int V,F; if(!(cin>>V>>F)) return 0;
    vector<P> origV; origV.reserve(V); vector<Tri> origF; origF.reserve(F);
    string s;
    for(int i=0;i<V;i++){ cin>>s; P p; cin>>p.x>>p.y>>p.z; origV.push_back(p); }
    for(int i=0;i<F;i++){ cin>>s; Tri t; cin>>t.a>>t.b>>t.c; --t.a;--t.b;--t.c; origF.push_back(t); }
    if(V==0 || F==0){ cout<<0<<" "<<0<<"\n"; return 0; }

    P mn={1e100,1e100,1e100}, mx={-1e100,-1e100,-1e100};
    for(auto&p:origV){ mn.x=min(mn.x,p.x); mn.y=min(mn.y,p.y); mn.z=min(mn.z,p.z); mx.x=max(mx.x,p.x); mx.y=max(mx.y,p.y); mx.z=max(mx.z,p.z); }
    double diag=norm(mx-mn), H=0.05*diag;
    if(H<=0 || !isfinite(H)){
        cout<<V<<" "<<F<<"\n";
        cout.setf(ios::fixed); cout<<setprecision(10);
        for(auto&p:origV) cout<<"v "<<p.x<<" "<<p.y<<" "<<p.z<<"\n";
        for(auto&t:origF) cout<<"f "<<t.a+1<<" "<<t.b+1<<" "<<t.c+1<<"\n";
        return 0;
    }

    vector<P> v=origV; vector<Tri> f=origF;
    int n=V;
    vector<char> alive(n,1);
    vector<int> sz(n,1);
    vector<double> rad(n,0.0);

    double minArea=max(1e-20, diag*diag*1e-18);
    int rounds = (V>120000?7:(V>50000?9:13));
    for(int round=0; round<rounds; ++round){
        unordered_map<long long,EdgeInfo> em; em.reserve(f.size()*4+10);
        vector<vector<int>> neigh(n), inc(n);
        for(int fi=0; fi<(int)f.size(); ++fi){
            Tri t=f[fi]; int a=t.a,b=t.b,c=t.c;
            inc[a].push_back(fi); inc[b].push_back(fi); inc[c].push_back(fi);
            neigh[a].push_back(b); neigh[a].push_back(c);
            neigh[b].push_back(a); neigh[b].push_back(c);
            neigh[c].push_back(a); neigh[c].push_back(b);
            int e0[3][3]={{a,b,c},{b,c,a},{c,a,b}};
            for(auto &e:e0){ long long k=key2(e[0],e[1]); auto &x=em[k]; x.cnt++; if(x.opp1<0)x.opp1=e[2]; else x.opp2=e[2]; }
        }
        for(int i=0;i<n;i++) if(alive[i]){ auto &w=neigh[i]; sort(w.begin(),w.end()); w.erase(unique(w.begin(),w.end()),w.end()); }
        struct Cand{double d; int a,b; int o1,o2;}; vector<Cand> cand; cand.reserve(em.size());
        double lenLim = H*(round<3?0.55:(round<8?0.72:0.90));
        for(auto &kv:em){
            if(kv.second.cnt!=2) continue;
            int a=(int)(kv.first>>32), b=(int)(kv.first & 0xffffffffu);
            if(!alive[a]||!alive[b]) continue;
            double d=norm(v[a]-v[b]); if(d>lenLim) continue;
            cand.push_back({d,a,b,kv.second.opp1,kv.second.opp2});
        }
        sort(cand.begin(),cand.end(),[](const Cand&x,const Cand&y){return x.d<y.d;});
        vector<int> mark(n,0); int changed=0;
        for(const auto &c:cand){
            int a=c.a,b=c.b; if(!alive[a]||!alive[b]||mark[a]||mark[b]) continue;
            if(sz[a]<sz[b] || (sz[a]==sz[b] && b<a)) swap(a,b);
            vector<int> inter; const auto&A=neigh[a]; const auto&B=neigh[b];
            size_t ia=0,ib=0; while(ia<A.size()&&ib<B.size()){
                if(A[ia]==B[ib]){inter.push_back(A[ia]); ia++; ib++;}
                else if(A[ia]<B[ib]) ia++; else ib++;
            }
            sort(inter.begin(),inter.end()); inter.erase(unique(inter.begin(),inter.end()),inter.end());
            if(inter.size()!=2) continue;
            if(!((inter[0]==min(c.o1,c.o2) && inter[1]==max(c.o1,c.o2)))) continue;
            P np=(v[a]*(double)sz[a]+v[b]*(double)sz[b])*(1.0/(sz[a]+sz[b]));
            double nr=max(norm(np-v[a])+rad[a], norm(np-v[b])+rad[b]);
            if(nr>H*0.985) continue;
            bool ok=true;
            for(int fi:inc[a]){
                Tri t=f[fi]; if(t.a==b||t.b==b||t.c==b) continue;
                Tri nt=t; if(nt.a==a){} if(nt.b==a){} if(nt.c==a){}
                auto getp=[&](int q)->P{return q==a?np:v[q];};
                P oldn=triN(v,t);
                P ncr=crossp(getp(nt.b)-getp(nt.a), getp(nt.c)-getp(nt.a));
                double nl=norm(ncr); if(nl==0){ok=false;break;}
                P nnew=ncr*(1.0/nl);
                if(dotp(oldn,nnew)<0.05 || nl<minArea){ok=false;break;}
            }
            if(!ok) continue;
            for(int fi:inc[b]){
                Tri t=f[fi]; if(t.a==a||t.b==a||t.c==a) continue;
                Tri nt=t; if(nt.a==b)nt.a=a; if(nt.b==b)nt.b=a; if(nt.c==b)nt.c=a;
                auto getp=[&](int q)->P{return q==a?np:v[q];};
                P oldn=triN(v,t);
                P ncr=crossp(getp(nt.b)-getp(nt.a), getp(nt.c)-getp(nt.a));
                double nl=norm(ncr); if(nl==0){ok=false;break;}
                P nnew=ncr*(1.0/nl);
                if(dotp(oldn,nnew)<0.05 || nl<minArea){ok=false;break;}
            }
            if(!ok) continue;
            v[a]=np; rad[a]=nr; sz[a]+=sz[b]; alive[b]=0; mark[a]=mark[b]=1; changed++;
            for(auto &t:f){ if(t.a==b)t.a=a; if(t.b==b)t.b=a; if(t.c==b)t.c=a; }
        }
        vector<Tri> nf; nf.reserve(f.size()); unordered_set<string> seen; seen.reserve(f.size()*2+10);
        for(auto&t:f){
            if(t.a==t.b||t.b==t.c||t.c==t.a) continue;
            if(!alive[t.a]||!alive[t.b]||!alive[t.c]) continue;
            if(triArea2(v,t)<minArea) continue;
            int x[3]={t.a,t.b,t.c}; int y[3]={x[0],x[1],x[2]}; sort(y,y+3);
            string key=to_string(y[0])+","+to_string(y[1])+","+to_string(y[2]);
            if(seen.insert(key).second) nf.push_back(t);
        }
        f.swap(nf);
        if(changed==0) break;
        if((int)f.size()<4) break;
    }

    vector<int> id(n,-1); vector<P> outV; outV.reserve(n);
    for(int i=0;i<n;i++) if(alive[i]){ id[i]=(int)outV.size(); outV.push_back(v[i]); }
    vector<Tri> outF; outF.reserve(f.size());
    for(auto&t:f){ if(id[t.a]<0||id[t.b]<0||id[t.c]<0) continue; Tri u={id[t.a],id[t.b],id[t.c]}; if(u.a!=u.b&&u.b!=u.c&&u.c!=u.a) outF.push_back(u); }

    if(outV.empty() || outF.empty() || !validClosed(outV,outF) || outV.size()>origV.size()){
        outV=origV; outF=origF;
    }

    cout.setf(ios::fixed); cout<<setprecision(10);
    cout<<outV.size()<<" "<<outF.size()<<"\n";
    for(auto&p:outV) cout<<"v "<<p.x<<" "<<p.y<<" "<<p.z<<"\n";
    for(auto&t:outF) cout<<"f "<<t.a+1<<" "<<t.b+1<<" "<<t.c+1<<"\n";
    return 0;
}