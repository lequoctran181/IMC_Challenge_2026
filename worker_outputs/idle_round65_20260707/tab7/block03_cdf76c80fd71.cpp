#include <bits/stdc++.h>
using namespace std;

struct Vec{
    double x,y,z;
    Vec(double X=0,double Y=0,double Z=0):x(X),y(Y),z(Z){}
    Vec operator+(const Vec&o)const{return Vec(x+o.x,y+o.y,z+o.z);} 
    Vec operator-(const Vec&o)const{return Vec(x-o.x,y-o.y,z-o.z);} 
    Vec operator*(double s)const{return Vec(x*s,y*s,z*s);} 
    Vec operator/(double s)const{return Vec(x/s,y/s,z/s);} 
};
static inline double dotv(const Vec&a,const Vec&b){return a.x*b.x+a.y*b.y+a.z*b.z;}
static inline Vec crossv(const Vec&a,const Vec&b){return Vec(a.y*b.z-a.z*b.y,a.z*b.x-a.x*b.z,a.x*b.y-a.y*b.x);} 
static inline double n2(const Vec&a){return dotv(a,a);} 
static inline double nv(const Vec&a){return sqrt(max(0.0,n2(a)));}
static inline Vec unit(const Vec&a){double l=nv(a);return l>0?a/l:Vec();}
struct Face{int a,b,c;};
static inline long long ekey(int a,int b){ if(a>b) swap(a,b); return ( (long long)a<<32 ) ^ (unsigned int)b; }
static inline long long tkey(int a,int b,int c){ if(a>b) swap(a,b); if(b>c) swap(b,c); if(a>b) swap(a,b); return ((long long)a*1000003LL+b)*1000003LL+c; }

static void emit(const vector<Vec>&P,const vector<Face>&F){
    cout<<P.size()<<" "<<F.size()<<"\n";
    cout.setf(ios::fixed); cout<<setprecision(10);
    for(const auto&p:P) cout<<"v "<<p.x<<" "<<p.y<<" "<<p.z<<"\n";
    for(const auto&f:F) cout<<"f "<<f.a+1<<" "<<f.b+1<<" "<<f.c+1<<"\n";
}

int main(){
    ios::sync_with_stdio(false); cin.tie(nullptr);
    int V,FN; if(!(cin>>V>>FN)) return 0;
    string tag; vector<Vec>P(V); vector<Face> orig; orig.reserve(FN);
    for(int i=0;i<V;i++) cin>>tag>>P[i].x>>P[i].y>>P[i].z;
    for(int i=0;i<FN;i++){int a,b,c;cin>>tag>>a>>b>>c;orig.push_back({a-1,b-1,c-1});}
    if(V<=8){ emit(P,orig); return 0; }

    Vec mn=P[0], mx=P[0];
    for(auto&p:P){mn.x=min(mn.x,p.x);mn.y=min(mn.y,p.y);mn.z=min(mn.z,p.z);mx.x=max(mx.x,p.x);mx.y=max(mx.y,p.y);mx.z=max(mx.z,p.z);} 
    double diag=max(1e-12,nv(mx-mn));
    double EPS=0.05*diag, SAFE=0.982*EPS;
    double area_eps=1e-14*diag*diag;

    vector<Face> faces=orig;
    vector<char> alive(V,1);
    vector<double> rad(V,0.0); // max distance from covered original vertices to this surviving output vertex

    auto area2 = [&](const Face&f)->double{return nv(crossv(P[f.b]-P[f.a],P[f.c]-P[f.a]));};
    auto normal = [&](const Face&f)->Vec{return unit(crossv(P[f.b]-P[f.a],P[f.c]-P[f.a]));};

    auto compact_clean = [&](){
        vector<Face> nf; nf.reserve(faces.size()); unordered_set<long long> seen; seen.reserve(faces.size()*2+1);
        for(auto f:faces){
            if(f.a<0||f.b<0||f.c<0||f.a>=V||f.b>=V||f.c>=V) continue;
            if(!alive[f.a]||!alive[f.b]||!alive[f.c]) continue;
            if(f.a==f.b||f.b==f.c||f.a==f.c) continue;
            if(area2(f)<=area_eps) continue;
            long long k=tkey(f.a,f.b,f.c);
            if(seen.insert(k).second) nf.push_back(f);
        }
        faces.swap(nf);
    };

    struct Cand{
        double score; int v, cover;
        vector<int> inc, ring;
        vector<Face> add;
        double newRad;
    };

    auto build_edge_counts = [&](){
        unordered_map<long long,int> ec; ec.reserve(faces.size()*4+1);
        for(auto f:faces){ ec[ekey(f.a,f.b)]++; ec[ekey(f.b,f.c)]++; ec[ekey(f.c,f.a)]++; }
        return ec;
    };

    auto try_candidate = [&](int v, const vector<vector<int>>&inc, const unordered_map<long long,int>&ec)->Cand{
        Cand fail; fail.score=1e100; fail.v=-1;
        if(!alive[v]) return fail;
        int deg=(int)inc[v].size();
        if(deg<3 || deg>22) return fail;

        unordered_map<int, vector<int>> link; link.reserve(deg*3+1);
        Vec sumN; double sumA=0, worstIncident=1.0;
        double nearest=1e100, farthest=0.0;
        vector<int> incident=inc[v];
        for(int fi:incident){
            const Face&f=faces[fi];
            int other[2],m=0;
            if(f.a!=v) other[m++]=f.a; if(f.b!=v) other[m++]=f.b; if(f.c!=v) other[m++]=f.c;
            if(m!=2 || !alive[other[0]] || !alive[other[1]]) return fail;
            link[other[0]].push_back(other[1]); link[other[1]].push_back(other[0]);
            Vec cr=crossv(P[f.b]-P[f.a],P[f.c]-P[f.a]); double a=nv(cr);
            if(a>0){ sumN=sumN+cr; sumA+=a; }
            nearest=min(nearest,nv(P[v]-P[other[0]])); nearest=min(nearest,nv(P[v]-P[other[1]]));
            farthest=max(farthest,nv(P[v]-P[other[0]])); farthest=max(farthest,nv(P[v]-P[other[1]]));
        }
        if((int)link.size()!=deg || sumA<=0) return fail;
        for(auto &kv:link) if(kv.second.size()!=2) return fail;
        Vec N=unit(sumN); if(nv(N)==0) return fail;
        for(int fi:incident){ Vec n=normal(faces[fi]); if(nv(n)>0) worstIncident=min(worstIncident,dotv(N,n)); }

        // Preserve creases/silhouette: allow flatter/lower-risk vertices only. Tiny one-rings can pass with more curvature.
        double curvatureGate = (nearest < 0.24*EPS ? -0.10 : (nearest < 0.48*EPS ? 0.18 : 0.50));
        if(worstIncident < curvatureGate) return fail;
        if(nearest + rad[v] > SAFE) return fail;

        // Reconstruct ordered cycle of neighbors around v.
        int start=link.begin()->first, cur=start, pre=-1;
        vector<int> ring; ring.reserve(deg);
        for(int step=0; step<deg+3; ++step){
            ring.push_back(cur);
            auto &w=link[cur]; int nx=(w[0]==pre?w[1]:w[0]);
            pre=cur; cur=nx; if(cur==start) break;
        }
        if(cur!=start || (int)ring.size()!=deg) return fail;
        vector<int> sortedRing=ring; sort(sortedRing.begin(),sortedRing.end());
        if(unique(sortedRing.begin(),sortedRing.end())!=sortedRing.end()) return fail;

        // Choose a cover vertex for Hausdorff. It must remain alive, and its existing covered set plus v's covered set
        // must all stay within 5% diagonal after v is removed.
        int cover=-1; double bestR=1e100;
        for(int u:ring){
            double nr=max(rad[u], rad[v] + nv(P[v]-P[u]));
            if(nr < bestR){ bestR=nr; cover=u; }
        }
        if(cover<0 || bestR>SAFE) return fail;

        // Try all fan roots; select the one with best normal preservation and shorter diagonals.
        vector<Face> bestAdd; double bestPenalty=1e100;
        for(int root=0; root<deg; ++root){
            int a=ring[root]; vector<Face> add; add.reserve(deg-2);
            double penalty=0; bool ok=true;
            unordered_map<long long,int> loc=ec;
            for(int fi:incident){
                Face f=faces[fi]; loc[ekey(f.a,f.b)]--; loc[ekey(f.b,f.c)]--; loc[ekey(f.c,f.a)]--;
            }
            for(int t=1;t<deg-1;t++){
                int b=ring[(root+t)%deg], c=ring[(root+t+1)%deg];
                if(a==b||b==c||a==c){ok=false;break;}
                Vec cr=crossv(P[b]-P[a],P[c]-P[a]); double ar=nv(cr);
                if(ar<=area_eps){ok=false;break;}
                if(dotv(cr,N)<0){ swap(b,c); cr=cr*(-1); }
                double co=dotv(unit(cr),N);
                if(co < -0.02){ok=false;break;}
                double longest=max(nv(P[a]-P[b]), max(nv(P[b]-P[c]), nv(P[c]-P[a])));
                // Do not bridge across a large visual gap unless the patch is almost planar.
                if(longest > 2.35*EPS && co < 0.82){ok=false;break;}
                loc[ekey(a,b)]++; loc[ekey(b,c)]++; loc[ekey(c,a)]++;
                if(loc[ekey(a,b)]>2 || loc[ekey(b,c)]>2 || loc[ekey(c,a)]>2){ok=false;break;}
                penalty += (1.0-co)*ar + 0.035*longest*longest;
                add.push_back({a,b,c});
            }
            if(!ok) continue;
            // Full local edge accounting: no zero-count or over-two edge in the affected patch.
            for(int x:ring){
                for(int y:link[x]){
                    long long k=ekey(x,y); if(loc[k]!=2){ok=false;break;}
                }
                if(!ok) break;
            }
            if(!ok) continue;
            for(auto f:add){
                if(loc[ekey(f.a,f.b)]!=2 || loc[ekey(f.b,f.c)]!=2 || loc[ekey(f.c,f.a)]!=2){ok=false;break;}
            }
            if(ok && penalty<bestPenalty){bestPenalty=penalty; bestAdd=add;}
        }
        if(bestAdd.empty()) return fail;

        // Keep bbox extremal vertices and strong axial-depth extrema: the six fixed cameras are sensitive to these.
        double margin=0.006*diag;
        if(P[v].x<mn.x+margin || P[v].x>mx.x-margin || P[v].y<mn.y+margin || P[v].y>mx.y-margin || P[v].z<mn.z+margin || P[v].z>mx.z-margin){
            if(nearest>0.18*EPS || worstIncident<0.92) return fail;
        }

        double score = nearest/EPS + 0.18*(1.0-worstIncident) + 0.006*deg + bestPenalty/(sumA+1e-30) + 0.08*(bestR/EPS);
        return {score,v,cover,incident,ring,bestAdd,bestR};
    };

    int aliveCount=V;
    for(int round=0; round<28; ++round){
        compact_clean();
        vector<vector<int>> inc(V);
        for(int i=0;i<(int)faces.size();++i){
            Face f=faces[i]; if(alive[f.a]) inc[f.a].push_back(i); if(alive[f.b]) inc[f.b].push_back(i); if(alive[f.c]) inc[f.c].push_back(i);
        }
        auto ec=build_edge_counts();
        vector<Cand> cand; cand.reserve(V/4+1);
        for(int v=0; v<V; ++v){
            Cand c=try_candidate(v,inc,ec); if(c.v>=0) cand.push_back(std::move(c));
        }
        if(cand.empty()) break;
        sort(cand.begin(),cand.end(),[](const Cand&a,const Cand&b){return a.score<b.score;});
        vector<char> reserved(V,0), killF(faces.size(),0), delV(V,0);
        vector<pair<int,double>> radUpdate;
        vector<Face> addFaces; int did=0;
        int budget=max(1, aliveCount/(round<5?4:(round<14?6:9)));
        for(auto &c:cand){
            if(did>=budget) break;
            if(!alive[c.v] || reserved[c.v] || reserved[c.cover]) continue;
            bool ok=true;
            for(int r:c.ring) if(!alive[r] || reserved[r]) {ok=false;break;}
            if(!ok) continue;
            for(int fi:c.inc) if(fi<0 || fi>=(int)faces.size() || killF[fi]) {ok=false;break;}
            if(!ok) continue;
            // Re-check coverage because cover may have been updated by an earlier accepted deletion.
            double nr=max(rad[c.cover], rad[c.v]+nv(P[c.v]-P[c.cover]));
            if(nr>SAFE) continue;
            reserved[c.v]=1; for(int r:c.ring) reserved[r]=1;
            delV[c.v]=1; for(int fi:c.inc) killF[fi]=1;
            for(auto f:c.add) addFaces.push_back(f);
            radUpdate.push_back({c.cover,nr});
            did++;
        }
        if(!did) break;
        vector<Face> nf; nf.reserve(faces.size()+addFaces.size());
        for(int i=0;i<(int)faces.size();++i) if(!killF[i]) nf.push_back(faces[i]);
        nf.insert(nf.end(),addFaces.begin(),addFaces.end()); faces.swap(nf);
        for(int i=0;i<V;i++) if(delV[i]){ alive[i]=0; aliveCount--; }
        for(auto &u:radUpdate) rad[u.first]=max(rad[u.first],u.second);
        if(aliveCount<=12) break;
    }
    compact_clean();

    // Validate closed two-manifold, nondegenerate faces, and all surviving/covered vertices within EPS.
    bool bad=false; unordered_map<long long,int> finalEC; finalEC.reserve(faces.size()*4+1); vector<int> used(V,0);
    for(auto f:faces){
        if(f.a<0||f.b<0||f.c<0||f.a>=V||f.b>=V||f.c>=V||!alive[f.a]||!alive[f.b]||!alive[f.c]||f.a==f.b||f.b==f.c||f.a==f.c){bad=true;break;}
        if(area2(f)<=area_eps){bad=true;break;}
        finalEC[ekey(f.a,f.b)]++; finalEC[ekey(f.b,f.c)]++; finalEC[ekey(f.c,f.a)]++;
        used[f.a]=used[f.b]=used[f.c]=1;
    }
    for(auto &kv:finalEC) if(kv.second!=2){bad=true;break;}
    for(int i=0;i<V;i++) if(alive[i] && rad[i]>EPS+1e-9){bad=true;break;}
    int outN=0; for(int i=0;i<V;i++) outN+=used[i];
    if(bad || faces.empty() || outN<4 || outN>=V){ emit(P,orig); return 0; }

    vector<int> id(V,-1); vector<Vec> OP; OP.reserve(outN);
    for(int i=0;i<V;i++) if(used[i]){ id[i]=(int)OP.size(); OP.push_back(P[i]); }
    vector<Face> OF; OF.reserve(faces.size()); unordered_set<long long> seen; seen.reserve(faces.size()*2+1);
    for(auto f:faces){
        int a=id[f.a], b=id[f.b], c=id[f.c]; if(a<0||b<0||c<0||a==b||b==c||a==c) continue;
        long long k=tkey(a,b,c); if(seen.insert(k).second) OF.push_back({a,b,c});
    }
    unordered_map<long long,int> ec2; ec2.reserve(OF.size()*4+1); bad=false;
    for(auto f:OF){ if(nv(crossv(OP[f.b]-OP[f.a],OP[f.c]-OP[f.a]))<=area_eps){bad=true;break;} ec2[ekey(f.a,f.b)]++; ec2[ekey(f.b,f.c)]++; ec2[ekey(f.c,f.a)]++; }
    for(auto &kv:ec2) if(kv.second!=2){bad=true;break;}
    if(bad || OF.empty()){ emit(P,orig); return 0; }
    emit(OP,OF);
    return 0;
}