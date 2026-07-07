#include <bits/stdc++.h>
using namespace std;
struct P{double x,y,z;};
struct T{int a,b,c;};
static inline P operator+(P a,P b){return {a.x+b.x,a.y+b.y,a.z+b.z};}
static inline P operator-(P a,P b){return {a.x-b.x,a.y-b.y,a.z-b.z};}
static inline P operator*(P a,double s){return {a.x*s,a.y*s,a.z*s};}
static inline double dotp(P a,P b){return a.x*b.x+a.y*b.y+a.z*b.z;}
static inline P crossp(P a,P b){return {a.y*b.z-a.z*b.y,a.z*b.x-a.x*b.z,a.x*b.y-a.y*b.x};}
static inline double n2(P a){return dotp(a,a);} static inline double nr(P a){return sqrt(max(0.0,n2(a)));}
static inline long long E(int a,int b){ if(a>b) swap(a,b); return ((long long)a<<32)^(unsigned)b; }
static inline double A2(const vector<P>&v,const T&t){return nr(crossp(v[t.b]-v[t.a],v[t.c]-v[t.a]));}
static bool closedOK(const vector<P>&v,const vector<T>&f){
    if(v.size()<4||f.size()<4) return false;
    P mn{1e100,1e100,1e100},mx{-1e100,-1e100,-1e100}; for(auto&p:v){mn.x=min(mn.x,p.x);mn.y=min(mn.y,p.y);mn.z=min(mn.z,p.z);mx.x=max(mx.x,p.x);mx.y=max(mx.y,p.y);mx.z=max(mx.z,p.z);} double eps=max(1e-18,n2(mx-mn)*1e-24);
    unordered_map<long long,int> ec; ec.reserve(f.size()*4+9);
    for(auto&t:f){ if(t.a<0||t.b<0||t.c<0||t.a>=(int)v.size()||t.b>=(int)v.size()||t.c>=(int)v.size()) return false; if(t.a==t.b||t.b==t.c||t.c==t.a) return false; if(A2(v,t)<=eps) return false; ec[E(t.a,t.b)]++; ec[E(t.b,t.c)]++; ec[E(t.c,t.a)]++; }
    for(auto &kv:ec) if(kv.second!=2) return false; return true;
}
static P Nrm(const vector<P>&v,const T&t){P q=crossp(v[t.b]-v[t.a],v[t.c]-v[t.a]); double l=nr(q); if(l<=0) return {0,0,0}; return q*(1.0/l);} 

int main(){
    ios::sync_with_stdio(false); cin.tie(nullptr);
    int N,M; if(!(cin>>N>>M)) return 0; vector<P> V(N), OV; vector<T> F(M), OF; string tag;
    for(int i=0;i<N;i++){cin>>tag>>V[i].x>>V[i].y>>V[i].z;} for(int i=0;i<M;i++){cin>>tag>>F[i].a>>F[i].b>>F[i].c;--F[i].a;--F[i].b;--F[i].c;} OV=V; OF=F;
    auto print=[&](const vector<P>&v,const vector<T>&f){ cout<<v.size()<<" "<<f.size()<<"\n"<<fixed<<setprecision(10); for(auto&p:v) cout<<"v "<<p.x<<" "<<p.y<<" "<<p.z<<"\n"; for(auto&t:f) cout<<"f "<<t.a+1<<" "<<t.b+1<<" "<<t.c+1<<"\n"; };
    if(N<=20){print(V,F);return 0;}
    P mn{1e100,1e100,1e100},mx{-1e100,-1e100,-1e100},cen{0,0,0}; for(auto&p:V){mn.x=min(mn.x,p.x);mn.y=min(mn.y,p.y);mn.z=min(mn.z,p.z);mx.x=max(mx.x,p.x);mx.y=max(mx.y,p.y);mx.z=max(mx.z,p.z);cen=cen+p;} cen=cen*(1.0/max(1,N)); double D=nr(mx-mn); if(!(D>0)) D=1; double H=.05*D;
    int target = max(64, min(N, (int)(sqrt((double)N)*22.0)));
    vector<int> ord(N); iota(ord.begin(),ord.end(),0);
    sort(ord.begin(),ord.end(),[&](int a,int b){return n2(V[a]-cen)>n2(V[b]-cen);});
    vector<int> land; land.reserve(target+8); vector<double> best(N,1e300);
    for(int step=0; step<target && step<N; ++step){
        int pick=-1; double bd=-1;
        int scan = (step<32?N:min(N, 25000 + step*700));
        for(int ii=0; ii<scan; ++ii){int i=ord[ii]; if(best[i]>bd){bd=best[i]; pick=i;}}
        if(pick<0) break; land.push_back(pick); for(int i=0;i<N;i++){double d=n2(V[i]-V[pick]); if(d<best[i]) best[i]=d;}
        if(step>20 && sqrt(bd)<H*.72) break;
    }
    double cover=0; for(double d:best) cover=max(cover,sqrt(d));
    vector<P> outV; vector<T> outF;
    bool tried=false;
    if(cover < H*0.86 && (int)land.size()>=4 && (int)land.size()<=180){
        tried=true; outV.reserve(land.size()); unordered_map<int,int> mp; for(int id:land){mp[id]=outV.size(); outV.push_back(V[id]);}
        int L=land.size(); vector<array<int,3>> hull; double eps=max(1e-12,D*1e-10);
        for(int i=0;i<L;i++) for(int j=i+1;j<L;j++) for(int k=j+1;k<L;k++){
            P a=outV[i],b=outV[j],c=outV[k]; P q=crossp(b-a,c-a); double ql=nr(q); if(ql<eps) continue;
            int pos=0,neg=0; for(int m=0;m<L;m++) if(m!=i&&m!=j&&m!=k){double s=dotp(q,outV[m]-a); if(s>eps) pos++; else if(s<-eps) neg++; if(pos&&neg) break;}
            if(pos&&neg) continue; if(pos==0&&neg==0) continue; if(pos==0) hull.push_back({i,k,j}); else hull.push_back({i,j,k});
        }
        sort(hull.begin(),hull.end()); hull.erase(unique(hull.begin(),hull.end()),hull.end()); outF.reserve(hull.size()); for(auto h:hull) outF.push_back({h[0],h[1],h[2]});
    }
    if(!tried || outF.empty() || !closedOK(outV,outF)){
        vector<char> av(N,1), af(F.size(),1); int FN=F.size();
        int passes = (N>80000?3:5);
        for(int pass=0; pass<passes; ++pass){
            vector<vector<int>> inc(N); inc.assign(N,{}); for(int i=0;i<FN;i++) if(af[i]){inc[F[i].a].push_back(i);inc[F[i].b].push_back(i);inc[F[i].c].push_back(i);} 
            vector<int> killV,killF; vector<T> add; vector<char> touched(N,0), tf(FN,0);
            struct C{double s;int v;int f0,f1,f2;T nt;}; vector<C> cs;
            for(int v=0; v<N; ++v) if(av[v] && inc[v].size()==3){
                int f0=inc[v][0],f1=inc[v][1],f2=inc[v][2]; vector<int> nb; for(int fid:inc[v]){T t=F[fid]; if(t.a!=v) nb.push_back(t.a); if(t.b!=v) nb.push_back(t.b); if(t.c!=v) nb.push_back(t.c);} sort(nb.begin(),nb.end()); nb.erase(unique(nb.begin(),nb.end()),nb.end()); if(nb.size()!=3) continue;
                if(min({nr(V[v]-V[nb[0]]),nr(V[v]-V[nb[1]]),nr(V[v]-V[nb[2]])})>H*.96) continue;
                P n0=Nrm(V,F[f0]),n1=Nrm(V,F[f1]),n2v=Nrm(V,F[f2]); P nav=n0+n1+n2v; double nl=nr(nav); if(nl<=0) continue; nav=nav*(1.0/nl); if(dotp(nav,n0)<.965||dotp(nav,n1)<.965||dotp(nav,n2v)<.965) continue;
                T nt{nb[0],nb[1],nb[2]}; P ntN=Nrm(V,nt); if(dotp(ntN,nav)<0) swap(nt.b,nt.c),ntN=Nrm(V,nt); if(dotp(ntN,nav)<.90) continue; if(A2(V,nt)<D*D*1e-18) continue;
                double dist=fabs(dotp(crossp(V[nt.b]-V[nt.a],V[nt.c]-V[nt.a]),V[v]-V[nt.a]))/max(1e-300,A2(V,nt)); if(dist>H*.25) continue;
                cs.push_back({dist/H + (1.0-min({dotp(nav,n0),dotp(nav,n1),dotp(nav,n2v)}))*10.0,v,f0,f1,f2,nt});
            }
            sort(cs.begin(),cs.end(),[](const C&a,const C&b){return a.s<b.s;});
            for(auto &c:cs){ if(!av[c.v]||touched[c.v]||!af[c.f0]||!af[c.f1]||!af[c.f2]) continue; T nt=c.nt; if(touched[nt.a]||touched[nt.b]||touched[nt.c]) continue; touched[c.v]=touched[nt.a]=touched[nt.b]=touched[nt.c]=1; av[c.v]=0; af[c.f0]=af[c.f1]=af[c.f2]=0; F.push_back(nt); af.push_back(1); FN++; }
            if(cs.empty()) break;
        }
        vector<int> id(N,-1); outV.clear(); for(int i=0;i<N;i++) if(av[i]){id[i]=outV.size(); outV.push_back(V[i]);}
        outF.clear(); unordered_set<string> seen; seen.reserve(FN*2+7); for(int i=0;i<FN;i++) if(af[i]){T t=F[i]; if(id[t.a]<0||id[t.b]<0||id[t.c]<0) continue; T u{id[t.a],id[t.b],id[t.c]}; int ar[3]={u.a,u.b,u.c}; sort(ar,ar+3); string key=to_string(ar[0])+","+to_string(ar[1])+","+to_string(ar[2]); if(u.a!=u.b&&u.b!=u.c&&u.c!=u.a&&seen.insert(key).second) outF.push_back(u);}        
    }
    if(!closedOK(outV,outF)) { outV=OV; outF=OF; }
    print(outV,outF); return 0;
}
