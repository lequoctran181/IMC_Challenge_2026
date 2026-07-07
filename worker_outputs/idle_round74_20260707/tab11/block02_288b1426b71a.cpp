#include <bits/stdc++.h>
using namespace std;
struct Vec{double x,y,z;};
static inline Vec operator+(Vec a,Vec b){return {a.x+b.x,a.y+b.y,a.z+b.z};}
static inline Vec operator-(Vec a,Vec b){return {a.x-b.x,a.y-b.y,a.z-b.z};}
static inline Vec operator*(Vec a,double s){return {a.x*s,a.y*s,a.z*s};}
static inline Vec operator/(Vec a,double s){return {a.x/s,a.y/s,a.z/s};}
static inline double dot(Vec a,Vec b){return a.x*b.x+a.y*b.y+a.z*b.z;}
static inline Vec cross(Vec a,Vec b){return {a.y*b.z-a.z*b.y,a.z*b.x-a.x*b.z,a.x*b.y-a.y*b.x};}
static inline double n2(Vec a){return dot(a,a);} static inline double norm(Vec a){return sqrt(max(0.0,n2(a)));}
struct Face{int a,b,c;bool alive=true;};
struct Vertex{Vec p;bool alive=true;vector<int> mem;};
struct EK{int a,b;EK(int u=0,int v=0){if(u<v){a=u;b=v;}else{a=v;b=u;}}bool operator==(EK const&o)const{return a==o.a&&b==o.b;}};
struct HK{size_t operator()(EK const&e)const{return (uint64_t(uint32_t(e.a))*11995408973635179863ull)^(uint32_t(e.b)+0x9e3779b9);}};
static inline Vec normal(Vec a,Vec b,Vec c){return cross(b-a,c-a);} 

int main(){
    ios::sync_with_stdio(false);cin.tie(nullptr);
    int V,F;if(!(cin>>V>>F)) return 0;
    vector<Vertex> v(V); vector<Face> f; f.reserve(F*2+100); vector<Vec> origP(V); vector<Face> origF(F);
    string t; Vec mn{1e100,1e100,1e100},mx{-1e100,-1e100,-1e100};
    for(int i=0;i<V;i++){cin>>t>>v[i].p.x>>v[i].p.y>>v[i].p.z;origP[i]=v[i].p;v[i].mem.push_back(i);mn.x=min(mn.x,v[i].p.x);mn.y=min(mn.y,v[i].p.y);mn.z=min(mn.z,v[i].p.z);mx.x=max(mx.x,v[i].p.x);mx.y=max(mx.y,v[i].p.y);mx.z=max(mx.z,v[i].p.z);}    
    for(int i=0;i<F;i++){Face q;cin>>t>>q.a>>q.b>>q.c;--q.a;--q.b;--q.c;f.push_back(q);origF[i]=q;}
    cout.setf(ios::fixed); cout<<setprecision(10);
    auto output_orig=[&](){cout<<V<<" "<<origF.size()<<"\n";for(auto&p:origP)cout<<"v "<<p.x<<" "<<p.y<<" "<<p.z<<"\n";for(auto&q:origF)cout<<"f "<<q.a+1<<" "<<q.b+1<<" "<<q.c+1<<"\n";};
    if(V<=8){output_orig();return 0;}
    double diag=norm(mx-mn), EPS=max(1e-12,0.0491*diag), EPS2=EPS*EPS; double AREA2=max(1e-32,diag*diag*diag*diag*1e-30);

    auto clean=[&](){
        unordered_set<string> seen; seen.reserve(f.size()*2+1);
        for(auto &q:f) if(q.alive){
            if(q.a<0||q.b<0||q.c<0||q.a>=V||q.b>=V||q.c>=V||!v[q.a].alive||!v[q.b].alive||!v[q.c].alive||q.a==q.b||q.b==q.c||q.c==q.a){q.alive=false;continue;}
            if(n2(normal(v[q.a].p,v[q.b].p,v[q.c].p))<AREA2){q.alive=false;continue;}
            int s[3]={q.a,q.b,q.c};sort(s,s+3);string k=to_string(s[0])+","+to_string(s[1])+","+to_string(s[2]);
            if(!seen.insert(k).second) q.alive=false;
        }
    };
    auto build=[&](vector<vector<int>>&vf, vector<unordered_set<int>>&nb, unordered_map<EK,vector<int>,HK>&ef){
        vf.assign(V,{}); nb.assign(V,{}); ef.clear(); ef.reserve(f.size()*4+1);
        for(int i=0;i<(int)f.size();++i) if(f[i].alive){int a=f[i].a,b=f[i].b,c=f[i].c;vf[a].push_back(i);vf[b].push_back(i);vf[c].push_back(i);nb[a].insert(b);nb[a].insert(c);nb[b].insert(a);nb[b].insert(c);nb[c].insert(a);nb[c].insert(b);ef[EK(a,b)].push_back(i);ef[EK(b,c)].push_back(i);ef[EK(c,a)].push_back(i);}    
    };
    auto within_any_live=[&](int orig, const vector<int>&cands){for(int x:cands) if(x>=0&&x<V&&v[x].alive&&n2(origP[orig]-v[x].p)<=EPS2) return true; return false;};
    auto pointOK=[&](int a,int b,Vec p){for(int id:v[a].mem) if(n2(origP[id]-p)>EPS2) return false;for(int id:v[b].mem) if(n2(origP[id]-p)>EPS2) return false;return true;};
    auto chooseCollapsePoint=[&](int a,int b,Vec&out){
        Vec cen{0,0,0};int cnt=0;for(int id:v[a].mem){cen=cen+origP[id];cnt++;}for(int id:v[b].mem){cen=cen+origP[id];cnt++;}
        vector<Vec> cand={v[a].p,v[b].p,(v[a].p+v[b].p)*0.5}; if(cnt)cand.push_back(cen/(double)cnt); cand.push_back(v[a].p*0.67+v[b].p*0.33); cand.push_back(v[a].p*0.33+v[b].p*0.67);
        double best=1e300; Vec bp=v[a].p; for(Vec p:cand) if(pointOK(a,b,p)){double d=0;for(int id:v[a].mem)d=max(d,n2(origP[id]-p));for(int id:v[b].mem)d=max(d,n2(origP[id]-p));if(d<best){best=d;bp=p;}}
        if(best>EPS2) return false; out=bp; return true;
    };
    auto localCollapseOK=[&](int a,int b,Vec np){
        for(auto&q:f) if(q.alive){bool h=(q.a==a||q.b==a||q.c==a||q.a==b||q.b==b||q.c==b); if(!h) continue; bool both=(q.a==a||q.b==a||q.c==a)&&(q.a==b||q.b==b||q.c==b); if(both) continue; Vec A=v[q.a].p,B=v[q.b].p,C=v[q.c].p; Vec n0=normal(A,B,C); if(q.a==a||q.a==b)A=np;if(q.b==a||q.b==b)B=np;if(q.c==a||q.c==b)C=np;Vec n1=normal(A,B,C);double l0=norm(n0),l1=norm(n1); if(n2(n1)<AREA2) return false; if(l0>1e-30&&l1>1e-30&&dot(n0,n1)<-0.10*l0*l1) return false;}
        return true;
    };
    auto simulateManifoldAfterCollapse=[&](int a,int b){
        unordered_map<EK,int,HK> cnt; unordered_set<string> tri; cnt.reserve(f.size()*4+1); tri.reserve(f.size()*2+1);
        for(auto&q:f) if(q.alive){int x=q.a,y=q.b,z=q.c;if(x==b)x=a;if(y==b)y=a;if(z==b)z=a;if(x==y||y==z||z==x) continue;int ss[3]={x,y,z};sort(ss,ss+3);string k=to_string(ss[0])+","+to_string(ss[1])+","+to_string(ss[2]);if(!tri.insert(k).second)return false;cnt[EK(x,y)]++;cnt[EK(y,z)]++;cnt[EK(z,x)]++;}
        for(auto&kv:cnt) if(kv.second!=2) return false; return true;
    };
    auto doCollapse=[&](int a,int b,Vec np){
        v[a].p=np; v[a].mem.insert(v[a].mem.end(),v[b].mem.begin(),v[b].mem.end()); v[b].alive=false;
        for(auto&q:f) if(q.alive){if(q.a==b)q.a=a;if(q.b==b)q.b=a;if(q.c==b)q.c=a;if(q.a==q.b||q.b==q.c||q.c==q.a) q.alive=false;}
    };

    // Phase A: star deletion/retriangulation. This is a hard pivot from pure collapse: remove a vertex if its original cluster is covered by neighbors and fill its one-ring polygon.
    for(int pass=0;pass<80;pass++){
        clean(); vector<vector<int>> vf; vector<unordered_set<int>> nb; unordered_map<EK,vector<int>,HK> ef; build(vf,nb,ef);
        int best=-1; double bestScore=1e300; vector<int> bestRing;
        for(int x=0;x<V;x++) if(v[x].alive){
            int deg=nb[x].size(); if(deg<3||deg>10) continue; if((int)vf[x].size()!=deg) continue;
            vector<int> ring; ring.reserve(deg);
            int start=-1; for(int y:nb[x]){start=y;break;} if(start<0) continue;
            int prev=-1,cur=start; bool ok=true;
            for(int step=0;step<deg;step++){
                ring.push_back(cur); int nxt=-1;
                for(int fi:vf[x]){int a=f[fi].a,b=f[fi].b,c=f[fi].c; if(a!=x&&b!=x&&c!=x) continue; int u=-1,w=-1; if(a==x){u=b;w=c;} else if(b==x){u=c;w=a;} else {u=a;w=b;} if(u==cur&&w!=prev)nxt=w; if(w==cur&&u!=prev)nxt=u;}
                if(nxt<0){ok=false;break;} prev=cur;cur=nxt; if(cur==start&&step+1<deg){ok=false;break;}
            }
            if(!ok||cur!=start) continue; set<int> uniq(ring.begin(),ring.end()); if((int)uniq.size()!=deg) continue;
            bool cover=true; for(int id:v[x].mem) if(!within_any_live(id,ring)){cover=false;break;} if(!cover) continue;
            Vec c{0,0,0}; for(int r:ring)c=c+v[r].p; c=c/(double)deg; Vec nsum{0,0,0}; for(int fi:vf[x]) nsum=nsum+normal(v[f[fi].a].p,v[f[fi].b].p,v[f[fi].c].p); int axis=0; Vec an{fabs(nsum.x),fabs(nsum.y),fabs(nsum.z)}; if(an.y>=an.x&&an.y>=an.z)axis=1; else if(an.z>=an.x&&an.z>=an.y)axis=2;
            auto orient2=[&](Vec p)->pair<double,double>{if(axis==0)return {p.y,p.z}; if(axis==1)return {p.x,p.z}; return {p.x,p.y};};
            double poly=0; for(int i=0;i<deg;i++){auto A=orient2(v[ring[i]].p),B=orient2(v[ring[(i+1)%deg]].p); poly+=A.first*B.second-A.second*B.first;} if(fabs(poly)<1e-30) continue;
            // fan from best ring vertex, no new vertices. Try all anchors and reject duplicates/inversions.
            double score=0; for(int r:ring) score+=n2(v[x].p-v[r].p); score/=deg; score*=deg-2;
            if(score<bestScore){bestScore=score;best=x;bestRing=ring;}
        }
        if(best<0) break;
        int deg=bestRing.size(); vector<int> old; for(int i=0;i<(int)f.size();++i) if(f[i].alive&&(f[i].a==best||f[i].b==best||f[i].c==best)) old.push_back(i);
        if((int)old.size()!=deg) continue;
        Vec nsum{0,0,0}; for(int fi:old) nsum=nsum+normal(v[f[fi].a].p,v[f[fi].b].p,v[f[fi].c].p);
        bool applied=false;
        for(int anchor=0;anchor<deg&&!applied;anchor++){
            vector<array<int,3>> add; bool ok=true; for(int k=1;k<deg-1;k++){array<int,3> tri={bestRing[anchor],bestRing[(anchor+k)%deg],bestRing[(anchor+k+1)%deg]}; Vec nn=normal(v[tri[0]].p,v[tri[1]].p,v[tri[2]].p); if(n2(nn)<AREA2){ok=false;break;} if(dot(nn,nsum)<0) swap(tri[1],tri[2]); add.push_back(tri);} if(!ok) continue;
            unordered_map<EK,int,HK> cnt; unordered_set<string> seen; bool bad=false;
            for(auto&q:f) if(q.alive){bool rem=(q.a==best||q.b==best||q.c==best); if(rem) continue; int ss[3]={q.a,q.b,q.c};sort(ss,ss+3);string key=to_string(ss[0])+","+to_string(ss[1])+","+to_string(ss[2]);seen.insert(key);cnt[EK(q.a,q.b)]++;cnt[EK(q.b,q.c)]++;cnt[EK(q.c,q.a)]++;}
            for(auto&t3:add){int ss[3]={t3[0],t3[1],t3[2]};sort(ss,ss+3);string key=to_string(ss[0])+","+to_string(ss[1])+","+to_string(ss[2]);if(!seen.insert(key).second){bad=true;break;}cnt[EK(t3[0],t3[1])]++;cnt[EK(t3[1],t3[2])]++;cnt[EK(t3[2],t3[0])]++;}
            if(bad) continue; for(auto&kv:cnt) if(kv.second!=2){bad=true;break;} if(bad) continue;
            for(int fi:old) f[fi].alive=false; v[best].alive=false; for(auto&t3:add){Face nf{t3[0],t3[1],t3[2],true}; f.push_back(nf);} applied=true;
        }
        if(!applied) break;
    }

    // Phase B: guarded local collapses for remaining short, smooth edges.
    for(int pass=0;pass<90;pass++){
        clean(); vector<vector<int>> vf; vector<unordered_set<int>> nb; unordered_map<EK,vector<int>,HK> ef; build(vf,nb,ef);
        struct C{double c;int a,b;}; vector<C> cand; cand.reserve(ef.size());
        for(auto&kv:ef){if(kv.second.size()!=2) continue;int a=kv.first.a,b=kv.first.b;if(!v[a].alive||!v[b].alive)continue;int common=0;for(int x:nb[a]) if(x!=b&&nb[b].count(x))common++;if(common!=2)continue;double l2=n2(v[a].p-v[b].p);if(l2>3.80*EPS2)continue;Vec n0=normal(v[f[kv.second[0]].a].p,v[f[kv.second[0]].b].p,v[f[kv.second[0]].c].p);Vec n1=normal(v[f[kv.second[1]].a].p,v[f[kv.second[1]].b].p,v[f[kv.second[1]].c].p);double cr=0;if(norm(n0)>0&&norm(n1)>0)cr=max(0.0,1.0-dot(n0,n1)/(norm(n0)*norm(n1)));cand.push_back({l2*(1+2.5*cr)*(1+0.008*(v[a].mem.size()+v[b].mem.size())),a,b});}
        sort(cand.begin(),cand.end(),[](auto&a,auto&b){return a.c<b.c;}); int did=0,limit=max(6,V/70);
        for(auto&c:cand){int a=c.a,b=c.b;if(!v[a].alive||!v[b].alive)continue;vector<vector<int>> vf2;vector<unordered_set<int>> nb2;unordered_map<EK,vector<int>,HK> ef2;build(vf2,nb2,ef2);auto it=ef2.find(EK(a,b));if(it==ef2.end()||it->second.size()!=2)continue;int common=0;for(int x:nb2[a])if(x!=b&&nb2[b].count(x))common++;if(common!=2)continue;if(v[a].mem.size()<v[b].mem.size())swap(a,b);Vec np;if(!chooseCollapsePoint(a,b,np))continue;if(!localCollapseOK(a,b,np))continue;if(!simulateManifoldAfterCollapse(a,b))continue;doCollapse(a,b,np);did++;if(did>=limit)break;}
        if(!did) break;
    }

    clean(); vector<int> id(V,-1); vector<Vec> ov; for(int i=0;i<V;i++) if(v[i].alive){id[i]=ov.size();ov.push_back(v[i].p);} vector<array<int,3>> of; set<array<int,3>> seen;
    for(auto&q:f) if(q.alive){if(id[q.a]<0||id[q.b]<0||id[q.c]<0)continue;array<int,3> tr={id[q.a],id[q.b],id[q.c]}, st=tr;sort(st.begin(),st.end());if(st[0]==st[1]||st[1]==st[2])continue;if(seen.insert(st).second)of.push_back(tr);}    
    auto finalOK=[&](){if(ov.size()<4||of.size()<4)return false;unordered_map<EK,int,HK> cnt;for(auto&tr:of){if(n2(normal(ov[tr[0]],ov[tr[1]],ov[tr[2]]))<AREA2)return false;cnt[EK(tr[0],tr[1])]++;cnt[EK(tr[1],tr[2])]++;cnt[EK(tr[2],tr[0])]++;}for(auto&kv:cnt)if(kv.second!=2)return false;return true;};
    if(!finalOK()){output_orig();return 0;}
    cout<<ov.size()<<" "<<of.size()<<"\n"; for(auto&p:ov) cout<<"v "<<p.x<<" "<<p.y<<" "<<p.z<<"\n"; for(auto&tr:of) cout<<"f "<<tr[0]+1<<" "<<tr[1]+1<<" "<<tr[2]+1<<"\n";
    return 0;
}
