#include <bits/stdc++.h>
using namespace std;
struct P{double x,y,z;};
struct F{int a,b,c;};
static inline P operator+(P a,P b){return {a.x+b.x,a.y+b.y,a.z+b.z};}
static inline P operator-(P a,P b){return {a.x-b.x,a.y-b.y,a.z-b.z};}
static inline P operator*(P a,double s){return {a.x*s,a.y*s,a.z*s};}
static inline double dotp(P a,P b){return a.x*b.x+a.y*b.y+a.z*b.z;}
static inline P crossp(P a,P b){return {a.y*b.z-a.z*b.y,a.z*b.x-a.x*b.z,a.x*b.y-a.y*b.x};}
static inline double n2(P a){return dotp(a,a);} static inline double nr(P a){return sqrt(max(0.0,n2(a)));}
static inline long long ekey(int a,int b){ if(a>b) swap(a,b); return ((long long)a<<32) ^ (unsigned)b; }
static inline string tkey(int a,int b,int c){ int x[3]={a,b,c}; sort(x,x+3); return to_string(x[0])+"#"+to_string(x[1])+"#"+to_string(x[2]); }
static inline double area2(const vector<P>&v,const F&t){return nr(crossp(v[t.b]-v[t.a],v[t.c]-v[t.a]));}
static inline P normal(const vector<P>&v,const F&t){P q=crossp(v[t.b]-v[t.a],v[t.c]-v[t.a]); double l=nr(q); if(l<=0) return {0,0,0}; return q*(1.0/l);} 
static bool validate(const vector<P>&v,const vector<F>&fs){
    if(v.empty()||fs.empty()) return false; unordered_map<long long,int> cnt; cnt.reserve(fs.size()*4+7);
    P mn={1e100,1e100,1e100},mx={-1e100,-1e100,-1e100}; for(auto&p:v){mn.x=min(mn.x,p.x);mn.y=min(mn.y,p.y);mn.z=min(mn.z,p.z);mx.x=max(mx.x,p.x);mx.y=max(mx.y,p.y);mx.z=max(mx.z,p.z);} double eps=max(1e-18,n2(mx-mn)*1e-24);
    for(auto&t:fs){ if(t.a<0||t.b<0||t.c<0||t.a>=(int)v.size()||t.b>=(int)v.size()||t.c>=(int)v.size()) return false; if(t.a==t.b||t.b==t.c||t.c==t.a) return false; if(area2(v,t)<=eps) return false; cnt[ekey(t.a,t.b)]++; cnt[ekey(t.b,t.c)]++; cnt[ekey(t.c,t.a)]++; }
    for(auto&kv:cnt) if(kv.second!=2) return false; return true;
}

int main(){
    ios::sync_with_stdio(false); cin.tie(nullptr);
    int N,M; if(!(cin>>N>>M)) return 0; vector<P> V(N); vector<F> T(M); string s;
    for(int i=0;i<N;i++){cin>>s>>V[i].x>>V[i].y>>V[i].z;}
    for(int i=0;i<M;i++){cin>>s>>T[i].a>>T[i].b>>T[i].c; --T[i].a;--T[i].b;--T[i].c;}
    vector<P> origV=V; vector<F> origT=T;
    if(N<=20){ cout<<N<<" "<<M<<"\n"<<fixed<<setprecision(10); for(auto&p:V) cout<<"v "<<p.x<<" "<<p.y<<" "<<p.z<<"\n"; for(auto&t:T) cout<<"f "<<t.a+1<<" "<<t.b+1<<" "<<t.c+1<<"\n"; return 0; }
    P mn={1e100,1e100,1e100},mx={-1e100,-1e100,-1e100}; for(auto&p:V){mn.x=min(mn.x,p.x);mn.y=min(mn.y,p.y);mn.z=min(mn.z,p.z);mx.x=max(mx.x,p.x);mx.y=max(mx.y,p.y);mx.z=max(mx.z,p.z);} double diag=nr(mx-mn); if(!(diag>0)) diag=1; double H=.05*diag; double epsA=max(1e-18,diag*diag*1e-18);
    vector<char> aliveV(N,1), aliveF(M,1);
    int maxPass = (N>100000?4:(N>30000?6:8));
    for(int pass=0; pass<maxPass; ++pass){
        vector<vector<int>> inc(N); vector<int> liveDeg(N,0); int liveF=0;
        for(int i=0;i<M;i++) if(aliveF[i]){auto&t=T[i]; if(!aliveV[t.a]||!aliveV[t.b]||!aliveV[t.c]){aliveF[i]=0; continue;} inc[t.a].push_back(i); inc[t.b].push_back(i); inc[t.c].push_back(i); liveDeg[t.a]++; liveDeg[t.b]++; liveDeg[t.c]++; liveF++;}
        unordered_set<string> outsideTri; outsideTri.reserve(liveF*2+7); for(int i=0;i<M;i++) if(aliveF[i]) outsideTri.insert(tkey(T[i].a,T[i].b,T[i].c));
        struct C{double sc; int v; vector<int> cyc; vector<F> add; vector<int> rem;}; vector<C> cand; cand.reserve(N/8+10);
        for(int r=0;r<N;r++) if(aliveV[r]){
            int k=(int)inc[r].size(); if(k<3||k>14) continue;
            unordered_map<int,vector<int>> ladj; ladj.reserve(k*3); vector<int> neigh; neigh.reserve(k); vector<int> rem=inc[r]; bool bad=false; P avn={0,0,0}; double minEdge=1e100, maxBend=0;
            for(int fid:inc[r]){ F t=T[fid]; int a=-1,b=-1; if(t.a==r){a=t.b;b=t.c;} else if(t.b==r){a=t.c;b=t.a;} else if(t.c==r){a=t.a;b=t.b;} else {bad=true;break;} if(a==b||!aliveV[a]||!aliveV[b]){bad=true;break;} ladj[a].push_back(b); ladj[b].push_back(a); neigh.push_back(a); neigh.push_back(b); P nn=normal(V,t); avn=avn+nn; minEdge=min(minEdge,min(nr(V[r]-V[a]),nr(V[r]-V[b]))); }
            if(bad||minEdge>H*.985) continue; sort(neigh.begin(),neigh.end()); neigh.erase(unique(neigh.begin(),neigh.end()),neigh.end()); if((int)neigh.size()!=k) continue; for(int x:neigh){ auto &w=ladj[x]; sort(w.begin(),w.end()); w.erase(unique(w.begin(),w.end()),w.end()); if(w.size()!=2){bad=true;break;} } if(bad) continue;
            double avl=nr(avn); if(avl<=0) continue; avn=avn*(1.0/avl); for(int fid:inc[r]){ P nn=normal(V,T[fid]); maxBend=max(maxBend,1.0-dotp(avn,nn)); } if(maxBend>0.10) continue;
            vector<int> cyc; cyc.reserve(k); int start=*min_element(neigh.begin(),neigh.end()), prev=-1, cur=start; for(int step=0; step<k; ++step){cyc.push_back(cur); auto w=ladj[cur]; int nxt=(w[0]==prev?w[1]:w[0]); prev=cur; cur=nxt;} if(cur!=start){ reverse(cyc.begin()+1,cyc.end()); prev=-1; cur=start; cyc.clear(); for(int step=0; step<k; ++step){cyc.push_back(cur); auto w=ladj[cur]; int nxt=(w[1]==prev?w[0]:w[1]); prev=cur; cur=nxt;} if(cur!=start) continue; }
            for(int fid:rem) outsideTri.erase(tkey(T[fid].a,T[fid].b,T[fid].c));
            vector<F> best; double bestScore=1e100; int kk=cyc.size();
            for(int root=0; root<kk; ++root){
                vector<int> ord; for(int q=0;q<kk;q++) ord.push_back(cyc[(root+q)%kk]);
                for(int flip=0; flip<2; ++flip){ vector<F> add; bool ok=true; double worst=0, flat=0, minar=1e100; for(int i=1;i+1<kk;i++){ F nt = flip?F{ord[0],ord[i+1],ord[i]}:F{ord[0],ord[i],ord[i+1]}; if(nt.a==nt.b||nt.b==nt.c||nt.c==nt.a){ok=false;break;} string tk=tkey(nt.a,nt.b,nt.c); if(outsideTri.count(tk)){ok=false;break;} double ar=area2(V,nt); if(ar<=epsA){ok=false;break;} P nn=normal(V,nt); double d=dotp(nn,avn); if(d<0.35){ok=false;break;} worst=max(worst,1.0-d); minar=min(minar,ar); P cr=crossp(V[nt.b]-V[nt.a],V[nt.c]-V[nt.a]); double cl=nr(cr); double dist = cl>0?fabs(dotp(cr,V[r]-V[nt.a]))/cl:1e100; flat=max(flat,dist); add.push_back(nt); }
                    if(!ok) continue; if(flat>H*0.38 && maxBend>0.025) continue; double sc=500*maxBend+40*worst+flat/H+0.0001*kk-minar/(diag*diag+1e-30); if(sc<bestScore){bestScore=sc; best=add;} }
            }
            for(int fid:rem) outsideTri.insert(tkey(T[fid].a,T[fid].b,T[fid].c));
            if(best.empty()) continue; cand.push_back({bestScore,r,cyc,best,rem});
        }
        sort(cand.begin(),cand.end(),[](const C&a,const C&b){return a.sc<b.sc;});
        vector<char> used(N,0), usedF(M,0), delV(N,0), delF(M,0); vector<F> addAll; int chosen=0;
        for(auto &c:cand){ if(!aliveV[c.v]||used[c.v]) continue; bool ok=true; for(int x:c.cyc) if(used[x]||!aliveV[x]){ok=false;break;} for(int fid:c.rem) if(usedF[fid]||!aliveF[fid]){ok=false;break;} if(!ok) continue; used[c.v]=1; delV[c.v]=1; for(int x:c.cyc) used[x]=1; for(int fid:c.rem){usedF[fid]=1; delF[fid]=1;} for(auto&t:c.add) addAll.push_back(t); chosen++; }
        if(!chosen) break;
        for(int i=0;i<N;i++) if(delV[i]) aliveV[i]=0; for(int i=0;i<M;i++) if(delF[i]) aliveF[i]=0; for(auto&t:addAll){ T.push_back(t); aliveF.push_back(1); M++; }
    }
    vector<int> id(N,-1); vector<P> outV; outV.reserve(N); for(int i=0;i<N;i++) if(aliveV[i]){id[i]=outV.size(); outV.push_back(V[i]);}
    vector<F> outF; outF.reserve(T.size()); unordered_set<string> seen; seen.reserve(T.size()*2+7); for(int i=0;i<(int)T.size();i++) if(aliveF[i]){F t=T[i]; if(id[t.a]<0||id[t.b]<0||id[t.c]<0) continue; F u{id[t.a],id[t.b],id[t.c]}; if(u.a==u.b||u.b==u.c||u.c==u.a) continue; string k=tkey(u.a,u.b,u.c); if(seen.insert(k).second) outF.push_back(u);}
    if(!validate(outV,outF)){ outV=origV; outF=origT; }
    cout<<outV.size()<<" "<<outF.size()<<"\n"<<fixed<<setprecision(10); for(auto&p:outV) cout<<"v "<<p.x<<" "<<p.y<<" "<<p.z<<"\n"; for(auto&t:outF) cout<<"f "<<t.a+1<<" "<<t.b+1<<" "<<t.c+1<<"\n";
    return 0;
}