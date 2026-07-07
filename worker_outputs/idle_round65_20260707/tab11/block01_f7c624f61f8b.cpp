#include <bits/stdc++.h>
using namespace std;

struct Vec{ double x,y,z; Vec(double X=0,double Y=0,double Z=0):x(X),y(Y),z(Z){} };
static inline Vec operator+(const Vec&a,const Vec&b){return {a.x+b.x,a.y+b.y,a.z+b.z};}
static inline Vec operator-(const Vec&a,const Vec&b){return {a.x-b.x,a.y-b.y,a.z-b.z};}
static inline Vec operator*(const Vec&a,double s){return {a.x*s,a.y*s,a.z*s};}
static inline Vec operator/(const Vec&a,double s){return {a.x/s,a.y/s,a.z/s};}
static inline double dotp(const Vec&a,const Vec&b){return a.x*b.x+a.y*b.y+a.z*b.z;}
static inline Vec crossp(const Vec&a,const Vec&b){return {a.y*b.z-a.z*b.y,a.z*b.x-a.x*b.z,a.x*b.y-a.y*b.x};}
static inline double n2(const Vec&a){return dotp(a,a);} static inline double len(const Vec&a){return sqrt(max(0.0,n2(a)));}
struct Face{int a,b,c; bool alive=true;};
struct Vertex{Vec p; bool alive=true; vector<int> mem;};
struct EKey{int a,b; EKey(int u=0,int v=0){ if(u<v){a=u;b=v;}else{a=v;b=u;} } bool operator==(const EKey&o)const{return a==o.a&&b==o.b;}};
struct EH{ size_t operator()(const EKey&e)const{ return (uint64_t(e.a)<<32) ^ uint32_t(e.b*1000003); } };
static inline Vec fnorm(const Vec&a,const Vec&b,const Vec&c){return crossp(b-a,c-a);}

int main(){
    ios::sync_with_stdio(false); cin.tie(nullptr);
    int V,F; if(!(cin>>V>>F)) return 0;
    vector<Vertex> v(V); vector<Face> f(F), origF(F); vector<Vec> origP(V);
    string s; Vec mn(1e100,1e100,1e100), mx(-1e100,-1e100,-1e100);
    for(int i=0;i<V;i++){
        cin>>s>>v[i].p.x>>v[i].p.y>>v[i].p.z; v[i].mem.push_back(i); origP[i]=v[i].p;
        mn.x=min(mn.x,v[i].p.x); mn.y=min(mn.y,v[i].p.y); mn.z=min(mn.z,v[i].p.z);
        mx.x=max(mx.x,v[i].p.x); mx.y=max(mx.y,v[i].p.y); mx.z=max(mx.z,v[i].p.z);
    }
    for(int i=0;i<F;i++){cin>>s>>f[i].a>>f[i].b>>f[i].c; --f[i].a;--f[i].b;--f[i].c; origF[i]=f[i];}

    cout.setf(ios::fixed); cout<<setprecision(10);
    if(V<=8){
        cout<<V<<" "<<F<<"\n";
        for(auto &p:origP) cout<<"v "<<p.x<<" "<<p.y<<" "<<p.z<<"\n";
        for(auto &q:origF) cout<<"f "<<q.a+1<<" "<<q.b+1<<" "<<q.c+1<<"\n";
        return 0;
    }

    double diag=len(mx-mn), eps=max(1e-12,0.04925*diag), eps2=eps*eps;
    double area2lim=max(1e-32, diag*diag*diag*diag*1e-28);

    auto sqdistToCluster=[&](const Vec&p,const vector<int>&A,const vector<int>&B){
        double m=0;
        for(int id:A) m=max(m,n2(origP[id]-p));
        for(int id:B) m=max(m,n2(origP[id]-p));
        return m;
    };
    auto choosePoint=[&](int a,int b,Vec&out)->bool{
        vector<Vec> cand;
        cand.push_back((v[a].p+v[b].p)*0.5);
        cand.push_back(v[a].p); cand.push_back(v[b].p);
        Vec cen; int cnt=0; for(int id:v[a].mem){cen=cen+origP[id];cnt++;} for(int id:v[b].mem){cen=cen+origP[id];cnt++;} if(cnt) cand.push_back(cen/double(cnt));
        cand.push_back(v[a].p*0.75+v[b].p*0.25); cand.push_back(v[a].p*0.25+v[b].p*0.75);
        double best=1e300; Vec bp;
        for(auto&p:cand){ double d=sqdistToCluster(p,v[a].mem,v[b].mem); if(d<=eps2 && d<best){best=d;bp=p;} }
        if(best>eps2) return false; out=bp; return true;
    };

    auto compactDeadFaces=[&](){
        unordered_set<string> seen; seen.reserve(f.size()*2+1);
        for(auto &q:f) if(q.alive){
            if(!v[q.a].alive||!v[q.b].alive||!v[q.c].alive||q.a==q.b||q.b==q.c||q.c==q.a){q.alive=false; continue;}
            int aa[3]={q.a,q.b,q.c}; sort(aa,aa+3); string key=to_string(aa[0])+","+to_string(aa[1])+","+to_string(aa[2]);
            if(!seen.insert(key).second) q.alive=false;
        }
    };

    auto buildAdj=[&](vector<unordered_set<int>>&nb, unordered_map<EKey,vector<int>,EH>&ef){
        nb.assign(V,{}); ef.clear(); ef.reserve(f.size()*4+1);
        for(int i=0;i<F;i++) if(f[i].alive){
            int a=f[i].a,b=f[i].b,c=f[i].c;
            nb[a].insert(b); nb[a].insert(c); nb[b].insert(a); nb[b].insert(c); nb[c].insert(a); nb[c].insert(b);
            ef[EKey(a,b)].push_back(i); ef[EKey(b,c)].push_back(i); ef[EKey(c,a)].push_back(i);
        }
    };

    auto localOK=[&](int keep,int del,const Vec&np)->bool{
        for(auto&q:f) if(q.alive){
            bool h=(q.a==keep||q.b==keep||q.c==keep||q.a==del||q.b==del||q.c==del);
            if(!h) continue;
            bool both=(q.a==keep||q.b==keep||q.c==keep) && (q.a==del||q.b==del||q.c==del);
            if(both) continue;
            Vec A=v[q.a].p,B=v[q.b].p,C=v[q.c].p; Vec n0=fnorm(A,B,C);
            if(q.a==keep||q.a==del) A=np; if(q.b==keep||q.b==del) B=np; if(q.c==keep||q.c==del) C=np;
            Vec n1=fnorm(A,B,C); double l0=len(n0), l1=len(n1);
            if(n2(n1)<area2lim) return false;
            if(l0>1e-28 && l1>1e-28 && dotp(n0,n1)<-0.02*l0*l1) return false;
        }
        return true;
    };

    auto collapse=[&](int keep,int del,const Vec&np){
        v[keep].p=np; v[keep].mem.insert(v[keep].mem.end(),v[del].mem.begin(),v[del].mem.end()); v[del].alive=false;
        for(auto&q:f) if(q.alive){ if(q.a==del)q.a=keep; if(q.b==del)q.b=keep; if(q.c==del)q.c=keep; if(q.a==q.b||q.b==q.c||q.c==q.a) q.alive=false; }
    };

    int alive=V;
    for(int pass=0; pass<120; ++pass){
        compactDeadFaces();
        vector<unordered_set<int>> nb; unordered_map<EKey,vector<int>,EH> ef; buildAdj(nb,ef);
        struct C{double cost;int a,b;}; vector<C> cs; cs.reserve(ef.size());
        for(auto &kv:ef){
            if(kv.second.size()!=2) continue; int a=kv.first.a,b=kv.first.b; if(!v[a].alive||!v[b].alive) continue;
            int common=0; for(int x:nb[a]) if(x!=b && nb[b].count(x)) common++;
            if(common!=2) continue;
            double l2=n2(v[a].p-v[b].p); if(l2>4.0*eps2) continue;
            Vec nsum1, nsum2; double amin=1e300;
            for(int fi:kv.second){ Vec n=fnorm(v[f[fi].a].p,v[f[fi].b].p,v[f[fi].c].p); double nl=len(n); amin=min(amin,nl); if(nl>0){ if(nsum1.x==0&&nsum1.y==0&&nsum1.z==0) nsum1=n/nl; else nsum2=n/nl; } }
            double crease=0; if(len(nsum1)>0&&len(nsum2)>0) crease=max(0.0,1.0-dotp(nsum1,nsum2));
            int mass=(int)v[a].mem.size()+(int)v[b].mem.size();
            cs.push_back({l2*(1.0+3.0*crease)*(1.0+0.01*mass),a,b});
        }
        sort(cs.begin(),cs.end(),[](const C&x,const C&y){return x.cost<y.cost;});
        int did=0, budget=max(12,V/45);
        for(auto &c:cs){
            int a=c.a,b=c.b; if(!v[a].alive||!v[b].alive) continue;
            vector<unordered_set<int>> nb2; unordered_map<EKey,vector<int>,EH> ef2; buildAdj(nb2,ef2);
            auto it=ef2.find(EKey(a,b)); if(it==ef2.end()||it->second.size()!=2) continue;
            int common=0; for(int x:nb2[a]) if(x!=b && nb2[b].count(x)) common++;
            if(common!=2) continue;
            if(v[a].mem.size()<v[b].mem.size()) swap(a,b);
            Vec np; if(!choosePoint(a,b,np)) continue; if(!localOK(a,b,np)) continue;
            collapse(a,b,np); alive--; did++; if(did>=budget) break;
        }
        if(did==0) break;
    }

    for(int pass=0; pass<40; ++pass){
        compactDeadFaces();
        vector<unordered_set<int>> nb; unordered_map<EKey,vector<int>,EH> ef; buildAdj(nb,ef);
        int did=0;
        for(int x=0;x<V;x++) if(v[x].alive && nb[x].size()==3){
            vector<int> r(nb[x].begin(),nb[x].end());
            if(!v[r[0]].alive||!v[r[1]].alive||!v[r[2]].alive) continue;
            bool close=true; for(int id:v[x].mem){ double d=min({n2(origP[id]-v[r[0]].p),n2(origP[id]-v[r[1]].p),n2(origP[id]-v[r[2]].p)}); if(d>eps2){close=false;break;} }
            if(!close) continue;
            Vec nn=fnorm(v[r[0]].p,v[r[1]].p,v[r[2]].p); if(n2(nn)<area2lim) continue;
            vector<int> inc; for(int i=0;i<F;i++) if(f[i].alive && (f[i].a==x||f[i].b==x||f[i].c==x)) inc.push_back(i);
            if(inc.size()!=3) continue;
            array<int,3> tri={r[0],r[1],r[2]}, st=tri; sort(st.begin(),st.end()); bool dup=false;
            for(auto&q:f) if(q.alive){ array<int,3> t={q.a,q.b,q.c}; sort(t.begin(),t.end()); if(t==st){dup=true;break;} }
            if(dup) continue;
            Vec sum; for(int i:inc) sum=sum+fnorm(v[f[i].a].p,v[f[i].b].p,v[f[i].c].p);
            if(dotp(nn,sum)<0) swap(tri[1],tri[2]);
            for(int i:inc) f[i].alive=false;
            Face nf; nf.a=tri[0]; nf.b=tri[1]; nf.c=tri[2]; nf.alive=true; f.push_back(nf); F++;
            v[x].alive=false; alive--; did++; break;
        }
        if(!did) break;
    }

    compactDeadFaces();
    vector<int> mapid(V,-1); vector<Vec> outv;
    for(int i=0;i<V;i++) if(v[i].alive){ mapid[i]=outv.size(); outv.push_back(v[i].p); }
    vector<array<int,3>> outf; set<array<int,3>> ss;
    for(auto&q:f) if(q.alive){
        if(mapid[q.a]<0||mapid[q.b]<0||mapid[q.c]<0) continue;
        array<int,3> t={mapid[q.a],mapid[q.b],mapid[q.c]}, st=t; sort(st.begin(),st.end());
        if(st[0]==st[1]||st[1]==st[2]) continue; if(ss.insert(st).second) outf.push_back(t);
    }
    auto finalOK=[&](){
        if(outv.size()<4||outf.size()<4) return false;
        unordered_map<EKey,int,EH> cnt; cnt.reserve(outf.size()*4+1);
        for(auto&t:outf){ if(n2(fnorm(outv[t[0]],outv[t[1]],outv[t[2]]))<area2lim) return false; cnt[EKey(t[0],t[1])]++; cnt[EKey(t[1],t[2])]++; cnt[EKey(t[2],t[0])]++; }
        for(auto&kv:cnt) if(kv.second!=2) return false;
        return true;
    };
    if(!finalOK()){
        cout<<V<<" "<<origF.size()<<"\n";
        for(auto&p:origP) cout<<"v "<<p.x<<" "<<p.y<<" "<<p.z<<"\n";
        for(auto&q:origF) cout<<"f "<<q.a+1<<" "<<q.b+1<<" "<<q.c+1<<"\n";
        return 0;
    }
    cout<<outv.size()<<" "<<outf.size()<<"\n";
    for(auto&p:outv) cout<<"v "<<p.x<<" "<<p.y<<" "<<p.z<<"\n";
    for(auto&t:outf) cout<<"f "<<t[0]+1<<" "<<t[1]+1<<" "<<t[2]+1<<"\n";
}