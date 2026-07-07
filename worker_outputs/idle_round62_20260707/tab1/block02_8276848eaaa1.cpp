#include <bits/stdc++.h>
using namespace std;

struct P{double x,y,z;};
struct F{int a,b,c;};
static inline P operator+(const P&a,const P&b){return {a.x+b.x,a.y+b.y,a.z+b.z};}
static inline P operator-(const P&a,const P&b){return {a.x-b.x,a.y-b.y,a.z-b.z};}
static inline P operator*(const P&a,double s){return {a.x*s,a.y*s,a.z*s};}
static inline double dotp(const P&a,const P&b){return a.x*b.x+a.y*b.y+a.z*b.z;}
static inline P crossp(const P&a,const P&b){return {a.y*b.z-a.z*b.y,a.z*b.x-a.x*b.z,a.x*b.y-a.y*b.x};}
static inline double n2(const P&a){return dotp(a,a);} static inline double normp(const P&a){return sqrt(max(0.0,n2(a)));}
static inline long long ekey(int a,int b){ if(a>b) swap(a,b); return ((long long)(unsigned)a<<32) | (unsigned)b; }
static inline P face_cross(const vector<P>&V,const F&f){ return crossp(V[f.b]-V[f.a],V[f.c]-V[f.a]); }
static inline double area2(const vector<P>&V,const F&f){ return normp(face_cross(V,f)); }
static inline P unit_or_zero(P q){ double l=normp(q); if(l<=0) return {0,0,0}; return q*(1.0/l); }

struct DSU{ vector<int> p,sz; DSU(int n=0){init(n);} void init(int n){p.resize(n); sz.assign(n,1); iota(p.begin(),p.end(),0);} int find(int x){while(p[x]!=x){p[x]=p[p[x]];x=p[x];}return x;} void unite(int a,int b){a=find(a);b=find(b); if(a==b)return; if(sz[a]<sz[b]) swap(a,b); p[b]=a; sz[a]+=sz[b];}};

static bool validate_closed_manifold(const vector<P>&V,const vector<F>&T){
    int N=(int)V.size(), M=(int)T.size(); if(N<4||M<4) return false;
    P mn{1e100,1e100,1e100}, mx{-1e100,-1e100,-1e100};
    for(auto&p:V){mn.x=min(mn.x,p.x);mn.y=min(mn.y,p.y);mn.z=min(mn.z,p.z);mx.x=max(mx.x,p.x);mx.y=max(mx.y,p.y);mx.z=max(mx.z,p.z);} double eps=max(1e-18,n2(mx-mn)*1e-24);
    unordered_map<long long,int> ec; ec.reserve((size_t)M*4+9);
    vector<int> vdeg(N,0); vector<vector<pair<int,int>>> link(N);
    for(const auto&f:T){
        if(f.a<0||f.b<0||f.c<0||f.a>=N||f.b>=N||f.c>=N) return false;
        if(f.a==f.b||f.b==f.c||f.c==f.a) return false;
        if(area2(V,f)<=eps) return false;
        ec[ekey(f.a,f.b)]++; ec[ekey(f.b,f.c)]++; ec[ekey(f.c,f.a)]++;
        vdeg[f.a]++; vdeg[f.b]++; vdeg[f.c]++;
        link[f.a].push_back({f.b,f.c}); link[f.b].push_back({f.c,f.a}); link[f.c].push_back({f.a,f.b});
    }
    for(auto &kv:ec) if(kv.second!=2) return false;
    for(int v=0; v<N; ++v){
        if(vdeg[v]==0) return false;
        unordered_map<int, array<int,3>> mp; mp.reserve(link[v].size()*3+5);
        for(auto pr:link[v]){
            int a=pr.first,b=pr.second; auto &A=mp[a]; if(A[0]<2) A[++A[0]]=b; else A[0]=3; auto &B=mp[b]; if(B[0]<2) B[++B[0]]=a; else B[0]=3;
        }
        if(mp.size()<3) return false;
        for(auto &kv:mp) if(kv.second[0]!=2 || kv.second[1]==kv.second[2]) return false;
        int start=mp.begin()->first, prev=-1, cur=start; size_t cnt=0;
        do{
            auto it=mp.find(cur); if(it==mp.end()) return false; int n1=it->second[1], n2v=it->second[2]; int nxt=(n1==prev?n2v:n1); prev=cur; cur=nxt; cnt++; if(cnt>mp.size()+1) return false;
        }while(cur!=start);
        if(cnt!=mp.size()) return false;
    }
    return true;
}

struct EdgeRec{long long key; int a,b,f;};
struct Cand{int gain=0; vector<pair<int,int>> maps;};

static bool build_output(const vector<P>&V,const vector<F>&F0,const vector<int>&rep,double H,vector<P>&outV,vector<F>&outF){
    int N=V.size(); vector<int> usedRep(N,0), incidentRep(N,0);
    for(int i=0;i<N;i++) usedRep[rep[i]]=1;
    outF.clear(); outF.reserve(F0.size()); unordered_set<string> seen; seen.reserve(F0.size()*2+7);
    P mn{1e100,1e100,1e100}, mx{-1e100,-1e100,-1e100}; for(auto&p:V){mn.x=min(mn.x,p.x);mn.y=min(mn.y,p.y);mn.z=min(mn.z,p.z);mx.x=max(mx.x,p.x);mx.y=max(mx.y,p.y);mx.z=max(mx.z,p.z);} double eps=max(1e-18,n2(mx-mn)*1e-24);
    vector<F> tmp; tmp.reserve(F0.size());
    for(auto f:F0){
        int a=rep[f.a], b=rep[f.b], c=rep[f.c]; if(a==b||b==c||c==a) continue; F t{a,b,c}; if(area2(V,t)<=eps) continue;
        int x[3]={a,b,c}; sort(x,x+3); string k=to_string(x[0])+"#"+to_string(x[1])+"#"+to_string(x[2]); if(!seen.insert(k).second) continue;
        tmp.push_back(t); incidentRep[a]=incidentRep[b]=incidentRep[c]=1;
    }
    for(int i=0;i<N;i++) if(usedRep[i]&&!incidentRep[i]) return false;
    vector<int> id(N,-1); outV.clear(); outV.reserve(N);
    for(int i=0;i<N;i++) if(usedRep[i]){ id[i]=(int)outV.size(); outV.push_back(V[i]); }
    outF.reserve(tmp.size()); for(auto&t:tmp) outF.push_back({id[t.a],id[t.b],id[t.c]});
    return validate_closed_manifold(outV,outF);
}

int main(){
    ios::sync_with_stdio(false); cin.tie(nullptr);
    int N,M; if(!(cin>>N>>M)) return 0; vector<P> V(N); vector<F> F0(M); string s;
    for(int i=0;i<N;i++){cin>>s>>V[i].x>>V[i].y>>V[i].z;}
    for(int i=0;i<M;i++){cin>>s>>F0[i].a>>F0[i].b>>F0[i].c; --F0[i].a; --F0[i].b; --F0[i].c;}
    auto print=[&](const vector<P>&VV,const vector<F>&FF){
        cout<<VV.size()<<" "<<FF.size()<<"\n"<<fixed<<setprecision(10);
        for(auto&p:VV) cout<<"v "<<p.x<<" "<<p.y<<" "<<p.z<<"\n";
        for(auto&f:FF) cout<<"f "<<f.a+1<<" "<<f.b+1<<" "<<f.c+1<<"\n";
    };
    if(N<=20 || M<20){ print(V,F0); return 0; }
    P mn{1e100,1e100,1e100}, mx{-1e100,-1e100,-1e100};
    for(auto&p:V){mn.x=min(mn.x,p.x);mn.y=min(mn.y,p.y);mn.z=min(mn.z,p.z);mx.x=max(mx.x,p.x);mx.y=max(mx.y,p.y);mx.z=max(mx.z,p.z);} double D=normp(mx-mn); if(!(D>0)) D=1; double H=0.05*D;
    vector<P> fn(M); vector<double> fa(M), fd(M);
    for(int i=0;i<M;i++){P cr=face_cross(V,F0[i]); fa[i]=normp(cr); fn[i]=fa[i]>0?cr*(1.0/fa[i]):P{0,0,0}; fd[i]=dotp(fn[i],V[F0[i].a]);}
    vector<EdgeRec> er; er.reserve((size_t)M*3);
    for(int i=0;i<M;i++){int a=F0[i].a,b=F0[i].b,c=F0[i].c; er.push_back({ekey(a,b),a,b,i}); er.push_back({ekey(b,c),b,c,i}); er.push_back({ekey(c,a),c,a,i});}
    sort(er.begin(),er.end(),[](const EdgeRec&x,const EdgeRec&y){return x.key<y.key;});
    DSU dsu(M); const double ndThr=0.99935, planeThr=H*0.075;
    for(size_t i=0,j;i<er.size();i=j){j=i+1;while(j<er.size()&&er[j].key==er[i].key)j++; if(j-i==2){int f1=er[i].f,f2=er[i+1].f; double nd=dotp(fn[f1],fn[f2]); if(nd>ndThr){double d1=fabs(dotp(fn[f1],V[F0[f2].a])-fd[f1]); double d2=fabs(dotp(fn[f2],V[F0[f1].a])-fd[f2]); if(max(d1,d2)<planeThr) dsu.unite(f1,f2);}}}
    vector<int> root(M); for(int i=0;i<M;i++) root[i]=dsu.find(i);
    unordered_map<int, vector<int>> compFaces; compFaces.reserve(M/4+7); for(int i=0;i<M;i++) compFaces[root[i]].push_back(i);
    unordered_map<int, vector<int>> compBoundary; compBoundary.reserve(compFaces.size()*2+7);
    for(size_t i=0,j;i<er.size();i=j){j=i+1;while(j<er.size()&&er[j].key==er[i].key)j++; if(j-i==2){int f1=er[i].f,f2=er[i+1].f; int r1=root[f1],r2=root[f2]; if(r1!=r2){compBoundary[r1].push_back(er[i].a); compBoundary[r1].push_back(er[i].b); compBoundary[r2].push_back(er[i].a); compBoundary[r2].push_back(er[i].b);}}}
    vector<Cand> cands;
    cands.reserve(compFaces.size());
    vector<char> isBoundary(N,0), inComp(N,0);
    for(auto &kv:compFaces){
        auto &fl=kv.second; if(fl.size()<18) continue;
        vector<int> verts; verts.reserve(fl.size()); P sumN{0,0,0}; double areaSum=0;
        for(int fid:fl){F f=F0[fid]; verts.push_back(f.a); verts.push_back(f.b); verts.push_back(f.c); sumN=sumN+fn[fid]*fa[fid]; areaSum+=fa[fid];}
        sort(verts.begin(),verts.end()); verts.erase(unique(verts.begin(),verts.end()),verts.end()); if(verts.size()<12) continue;
        auto bit=compBoundary.find(kv.first); if(bit==compBoundary.end()) continue; auto bnd=bit->second; sort(bnd.begin(),bnd.end()); bnd.erase(unique(bnd.begin(),bnd.end()),bnd.end()); if(bnd.size()<3 || bnd.size()>=verts.size()) continue;
        P Nrm=unit_or_zero(sumN); if(normp(Nrm)<=0) continue; double d0=0; for(int fid:fl) d0 += fd[fid]*fa[fid]; d0/=max(areaSum,1e-300);
        double maxDev=0, minNd=1; for(int fid:fl){minNd=min(minNd,dotp(Nrm,fn[fid])); F f=F0[fid]; maxDev=max(maxDev,fabs(dotp(Nrm,V[f.a])-d0)); maxDev=max(maxDev,fabs(dotp(Nrm,V[f.b])-d0)); maxDev=max(maxDev,fabs(dotp(Nrm,V[f.c])-d0));}
        if(minNd<0.9965 || maxDev>H*0.22) continue;
        P ax = fabs(Nrm.x)<0.7?P{1,0,0}:(fabs(Nrm.y)<0.7?P{0,1,0}:P{0,0,1}); P U=unit_or_zero(crossp(Nrm,ax)); P W=crossp(Nrm,U);
        for(int x:verts) inComp[x]=1; for(int x:bnd) isBoundary[x]=1;
        double gs=H*0.43; if(!(gs>0)) gs=1;
        unordered_map<long long,int> cellRep; cellRep.reserve(verts.size()*2+7);
        vector<pair<int,int>> maps; maps.reserve(verts.size());
        for(int x:verts){ if(isBoundary[x]) continue; double u=dotp(V[x],U)/gs, w=dotp(V[x],W)/gs; int iu=(int)floor(u), iw=(int)floor(w); long long key=((long long)(unsigned)(iu+1000000000)<<32) ^ (unsigned)(iw+1000000000); auto it=cellRep.find(key); if(it==cellRep.end()) cellRep[key]=x; }
        bool ok=true; for(int x:verts){ if(isBoundary[x]) continue; double u=dotp(V[x],U)/gs,w=dotp(V[x],W)/gs; int iu=(int)floor(u),iw=(int)floor(w); long long key=((long long)(unsigned)(iu+1000000000)<<32) ^ (unsigned)(iw+1000000000); int r=cellRep[key]; if(normp(V[x]-V[r])>H*0.88){ok=false;break;} if(r!=x) maps.push_back({x,r}); }
        int oldInterior=0; for(int x:verts) if(!isBoundary[x]) oldInterior++; int newInterior=(int)cellRep.size(); int gain=oldInterior-newInterior;
        for(int x:verts) inComp[x]=0; for(int x:bnd) isBoundary[x]=0;
        if(!ok || gain<3) continue; Cand c; c.gain=gain; c.maps.swap(maps); cands.push_back(move(c));
    }
    sort(cands.begin(),cands.end(),[](const Cand&a,const Cand&b){return a.gain>b.gain;});
    vector<int> rep(N); iota(rep.begin(),rep.end(),0); vector<P> bestV,trialV; vector<F> bestF,trialF;
    bool have=false; int baseN=N;
    vector<int> allRep=rep;
    for(auto &c:cands) for(auto pr:c.maps) allRep[pr.first]=pr.second;
    if(!cands.empty() && build_output(V,F0,allRep,H,trialV,trialF) && (int)trialV.size()<baseN){ bestV=trialV; bestF=trialF; rep=allRep; have=true; }
    else{
        auto start=chrono::steady_clock::now(); int tests=0;
        for(auto &c:cands){
            bool conflict=false; for(auto pr:c.maps){ if(rep[pr.first]!=pr.first || rep[pr.second]!=pr.second){conflict=true;break;} } if(conflict) continue;
            vector<int> old; old.reserve(c.maps.size()); for(auto pr:c.maps){old.push_back(pr.first); rep[pr.first]=pr.second;}
            tests++; bool ok=build_output(V,F0,rep,H,trialV,trialF) && (int)trialV.size() < (have?(int)bestV.size():baseN);
            if(ok){bestV=trialV; bestF=trialF; have=true;} else for(int x:old) rep[x]=x;
            double el=chrono::duration<double>(chrono::steady_clock::now()-start).count(); if(tests>=45 || el>7.5) break;
        }
    }
    if(!have){ print(V,F0); return 0; }
    print(bestV,bestF); return 0;
}