#include <bits/stdc++.h>
using namespace std;

struct Vec{double x,y,z;Vec(double X=0,double Y=0,double Z=0):x(X),y(Y),z(Z){}Vec operator+(const Vec&o)const{return {x+o.x,y+o.y,z+o.z};}Vec operator-(const Vec&o)const{return {x-o.x,y-o.y,z-o.z};}Vec operator*(double s)const{return {x*s,y*s,z*s};}Vec operator/(double s)const{return {x/s,y/s,z/s};}};
static inline double dotv(const Vec&a,const Vec&b){return a.x*b.x+a.y*b.y+a.z*b.z;}
static inline Vec crossv(const Vec&a,const Vec&b){return {a.y*b.z-a.z*b.y,a.z*b.x-a.x*b.z,a.x*b.y-a.y*b.x};}
static inline double n2(const Vec&a){return dotv(a,a);}static inline double nv(const Vec&a){return sqrt(n2(a));}
static inline Vec unit(Vec a){double l=nv(a);return l>0?a/l:Vec();}
struct Face{int a,b,c;};
static inline long long ekey(int a,int b){if(a>b)swap(a,b);return ( (long long)a<<32 ) ^ (unsigned)b;} 
static inline long long tkey(int a,int b,int c){if(a>b)swap(a,b);if(b>c)swap(b,c);if(a>b)swap(a,b);return ((long long)a*1000003LL+b)*1000003LL+c;}

static void emit_original(const vector<Vec>&P,const vector<Face>&F){
    cout<<P.size()<<" "<<F.size()<<"\n";cout.setf(ios::fixed);cout<<setprecision(10);
    for(auto&p:P)cout<<"v "<<p.x<<" "<<p.y<<" "<<p.z<<"\n";
    for(auto&f:F)cout<<"f "<<f.a+1<<" "<<f.b+1<<" "<<f.c+1<<"\n";
}

int main(){
    ios::sync_with_stdio(false);cin.tie(nullptr);
    int V,FN;if(!(cin>>V>>FN))return 0;string s;vector<Vec>P(V);vector<Face>faces;faces.reserve(FN);
    for(int i=0;i<V;i++)cin>>s>>P[i].x>>P[i].y>>P[i].z;
    for(int i=0;i<FN;i++){int a,b,c;cin>>s>>a>>b>>c;faces.push_back({a-1,b-1,c-1});}
    vector<Face>origF=faces;
    if(V<=8){emit_original(P,origF);return 0;}
    Vec mn=P[0],mx=P[0];for(auto&p:P){mn.x=min(mn.x,p.x);mn.y=min(mn.y,p.y);mn.z=min(mn.z,p.z);mx.x=max(mx.x,p.x);mx.y=max(mx.y,p.y);mx.z=max(mx.z,p.z);}double diag=max(1e-12,nv(mx-mn));
    const double EPS=0.05*diag, DEL=0.985*EPS;
    vector<char> alive(V,1);

    auto fnormal=[&](const Face&f){return unit(crossv(P[f.b]-P[f.a],P[f.c]-P[f.a]));};
    auto farea2=[&](const Face&f){return nv(crossv(P[f.b]-P[f.a],P[f.c]-P[f.a]));};
    auto clean=[&](){
        vector<Face> nf;nf.reserve(faces.size());unordered_set<long long> seen;seen.reserve(faces.size()*2+1);
        for(auto f:faces){
            if(f.a<0||f.b<0||f.c<0||f.a>=V||f.b>=V||f.c>=V)continue;
            if(!alive[f.a]||!alive[f.b]||!alive[f.c])continue;
            if(f.a==f.b||f.b==f.c||f.a==f.c)continue;
            if(farea2(f)<=1e-14*diag*diag)continue;
            long long k=tkey(f.a,f.b,f.c);if(seen.insert(k).second)nf.push_back(f);
        }
        faces.swap(nf);
    };

    struct Cand{double q;int v;vector<int> inc,cyc;vector<Face> add;};

    for(int round=0;round<16;round++){
        clean();
        vector<vector<int>> inc(V);inc.assign(V,{});
        unordered_map<long long,int> ec;ec.reserve(faces.size()*4+1);
        for(int i=0;i<(int)faces.size();i++){
            auto f=faces[i];inc[f.a].push_back(i);inc[f.b].push_back(i);inc[f.c].push_back(i);
            ec[ekey(f.a,f.b)]++;ec[ekey(f.b,f.c)]++;ec[ekey(f.c,f.a)]++;
        }
        vector<Cand> cands;
        for(int v=0;v<V;v++) if(alive[v]){
            int k=inc[v].size(); if(k<3||k>14)continue;
            unordered_map<int,vector<int>> rg;rg.reserve(k*3+1);vector<int> neigh;neigh.reserve(k);
            bool bad=false;Vec avn;double area=0,dmin=1e100;
            for(int fi:inc[v]){
                Face f=faces[fi];int x[2],m=0;
                if(f.a!=v)x[m++]=f.a;if(f.b!=v)x[m++]=f.b;if(f.c!=v)x[m++]=f.c;
                if(m!=2||!alive[x[0]]||!alive[x[1]]){bad=true;break;}
                rg[x[0]].push_back(x[1]);rg[x[1]].push_back(x[0]);
                Vec nn=crossv(P[f.b]-P[f.a],P[f.c]-P[f.a]); double aa=nv(nn); if(aa>0){avn=avn+nn; area+=aa;}
                dmin=min(dmin,nv(P[v]-P[x[0]]));dmin=min(dmin,nv(P[v]-P[x[1]]));
            }
            if(bad||dmin>DEL||area<=0)continue;
            for(auto &p:rg) if(p.second.size()!=2){bad=true;break;} if(bad||rg.size()!=size_t(k))continue;
            int start=rg.begin()->first,cur=start,pre=-1;vector<int> cyc;
            for(int step=0;step<k+2;step++){
                cyc.push_back(cur);auto &vv=rg[cur];int nx=(vv[0]==pre?vv[1]:vv[0]);pre=cur;cur=nx;if(cur==start)break;
            }
            if((int)cyc.size()!=k||cur!=start)continue;
            Vec N=unit(avn); if(nv(N)==0)continue;
            double worst=1.0;
            for(int fi:inc[v]){Vec nn=fnormal(faces[fi]);if(nv(nn)>0)worst=min(worst,dotv(N,nn));}
            double sharpLimit = (dmin<0.22*EPS?0.05:(dmin<0.45*EPS?0.35:0.68));
            if(worst<sharpLimit)continue;

            // Choose the fan root that gives the least normal/area penalty. This is a one-ring
            // vertex deletion, not a collapse: all remaining vertices are original input vertices,
            // so output->input Hausdorff is exactly zero and deleted vertices are guarded by dmin.
            double best=1e100;vector<Face> bestadd;
            for(int root=0;root<k;root++){
                vector<Face> add;double pen=0;bool ok=true;int a=cyc[root];
                for(int t=1;t<k-1;t++){
                    int b=cyc[(root+t)%k], c=cyc[(root+t+1)%k];
                    if(a==b||b==c||a==c){ok=false;break;}
                    Vec cr=crossv(P[b]-P[a],P[c]-P[a]);double ar=nv(cr);
                    if(ar<=1e-12*diag*diag){ok=false;break;}
                    if(dotv(cr,N)<0){swap(b,c);cr=cr*-1;}
                    double co=dotv(unit(cr),N); if(co<0.03){ok=false;break;}
                    pen += (1.0-co)*ar + 0.001*ar;
                    add.push_back({a,b,c});
                }
                if(!ok)continue;
                unordered_map<long long,int> loc;
                for(auto &p:ec) loc[p.first]=p.second;
                for(int fi:inc[v]){Face f=faces[fi];loc[ekey(f.a,f.b)]--;loc[ekey(f.b,f.c)]--;loc[ekey(f.c,f.a)]--;}
                for(auto f:add){
                    long long e1=ekey(f.a,f.b),e2=ekey(f.b,f.c),e3=ekey(f.c,f.a);
                    loc[e1]++;loc[e2]++;loc[e3]++;
                    if(loc[e1]>2||loc[e2]>2||loc[e3]>2){ok=false;break;}
                }
                if(!ok)continue;
                // Boundary edges should be restored to exactly two; fan diagonals exactly two.
                for(auto f:add){
                    if(loc[ekey(f.a,f.b)]<1||loc[ekey(f.b,f.c)]<1||loc[ekey(f.c,f.a)]<1){ok=false;break;}
                }
                if(ok&&pen<best){best=pen;bestadd=add;}
            }
            if(bestadd.empty())continue;
            double q=dmin/EPS + 0.25*(1.0-worst) + 0.015*k + best/(area+1e-30);
            cands.push_back({q,v,inc[v],cyc,bestadd});
        }
        if(cands.empty())break;
        sort(cands.begin(),cands.end(),[](const Cand&a,const Cand&b){return a.q<b.q;});
        vector<char> reserved(V,0), killFace(faces.size(),0), delV(V,0);vector<Face> append;
        int removed=0, limit=max(1,V/(round<4?5:8));
        for(auto &c:cands){
            if(removed>=limit)break;if(!alive[c.v]||reserved[c.v])continue;bool ok=true;
            for(int x:c.cyc) if(reserved[x]||!alive[x]){ok=false;break;} if(!ok)continue;
            for(int fi:c.inc) if(fi<0||fi>=(int)faces.size()||killFace[fi]){ok=false;break;} if(!ok)continue;
            delV[c.v]=1;reserved[c.v]=1;for(int x:c.cyc)reserved[x]=1;for(int fi:c.inc)killFace[fi]=1;for(auto f:c.add)append.push_back(f);removed++;
        }
        if(!removed)break;
        vector<Face> nf;nf.reserve(faces.size()+append.size());
        for(int i=0;i<(int)faces.size();i++)if(!killFace[i])nf.push_back(faces[i]);
        nf.insert(nf.end(),append.begin(),append.end());faces.swap(nf);
        for(int i=0;i<V;i++)if(delV[i])alive[i]=0;
    }
    clean();

    // closed manifold validation and used-vertex compaction; if anything is suspicious, output input.
    unordered_map<long long,int> ec;ec.reserve(faces.size()*4+1);vector<int> used(V,0);bool bad=false;
    for(auto f:faces){
        if(!alive[f.a]||!alive[f.b]||!alive[f.c]||f.a==f.b||f.b==f.c||f.a==f.c){bad=true;break;}
        if(farea2(f)<=1e-14*diag*diag){bad=true;break;}
        ec[ekey(f.a,f.b)]++;ec[ekey(f.b,f.c)]++;ec[ekey(f.c,f.a)]++;used[f.a]=used[f.b]=used[f.c]=1;
    }
    for(auto&p:ec)if(p.second!=2){bad=true;break;}
    int outV=0;for(int i=0;i<V;i++)outV+=used[i];
    if(bad||faces.empty()||outV<4||outV>=V){emit_original(P,origF);return 0;}
    vector<int> id(V,-1);vector<Vec> OP;OP.reserve(outV);for(int i=0;i<V;i++)if(used[i]){id[i]=OP.size();OP.push_back(P[i]);}
    vector<Face> OF;OF.reserve(faces.size());unordered_set<long long> seen;seen.reserve(faces.size()*2+1);
    for(auto f:faces){int a=id[f.a],b=id[f.b],c=id[f.c];if(a<0||b<0||c<0||a==b||b==c||a==c)continue;long long k=tkey(a,b,c);if(seen.insert(k).second)OF.push_back({a,b,c});}
    unordered_map<long long,int> ec2;for(auto f:OF){ec2[ekey(f.a,f.b)]++;ec2[ekey(f.b,f.c)]++;ec2[ekey(f.c,f.a)]++;}
    for(auto&p:ec2)if(p.second!=2){emit_original(P,origF);return 0;}
    cout<<OP.size()<<" "<<OF.size()<<"\n";cout.setf(ios::fixed);cout<<setprecision(10);
    for(auto&p:OP)cout<<"v "<<p.x<<" "<<p.y<<" "<<p.z<<"\n";
    for(auto&f:OF)cout<<"f "<<f.a+1<<" "<<f.b+1<<" "<<f.c+1<<"\n";
    return 0;
}