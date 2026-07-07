#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <iostream>
#include <limits>
#include <queue>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>
using namespace std;

struct Vec3{ double x=0,y=0,z=0; };
static inline Vec3 operator+(const Vec3&a,const Vec3&b){ return {a.x+b.x,a.y+b.y,a.z+b.z}; }
static inline Vec3 operator-(const Vec3&a,const Vec3&b){ return {a.x-b.x,a.y-b.y,a.z-b.z}; }
static inline Vec3 operator*(const Vec3&a,double s){ return {a.x*s,a.y*s,a.z*s}; }
static inline double dot3(const Vec3&a,const Vec3&b){ return a.x*b.x+a.y*b.y+a.z*b.z; }
static inline Vec3 cross3(const Vec3&a,const Vec3&b){ return {a.y*b.z-a.z*b.y,a.z*b.x-a.x*b.z,a.x*b.y-a.y*b.x}; }
static inline double norm2(const Vec3&a){ return dot3(a,a); }
static inline double norm3(const Vec3&a){ return sqrt(max(0.0,norm2(a))); }

struct Face{ int a,b,c; };

static int N, M0;
static vector<Vec3> P;
static vector<Face> F;
static vector<unsigned char> aliveV, aliveF;
static vector<vector<int>> VF;
static vector<unsigned> stampV;
static int liveV=0, liveF=0;
static double diagLen=1.0, epsArea=1e-18, deleteRadius=1.0;

static inline uint64_t ekey(int a,int b){ if(a>b) swap(a,b); return (uint64_t)(uint32_t)a<<32 | (uint32_t)b; }
static inline array<int,3> fkey(int a,int b,int c){ array<int,3>s={a,b,c}; sort(s.begin(),s.end()); return s; }
static inline bool hasv(const Face&f,int v){ return f.a==v||f.b==v||f.c==v; }
static inline Vec3 faceCross(const Face&f){ return cross3(P[f.b]-P[f.a], P[f.c]-P[f.a]); }

static vector<int> incident(int v){
    vector<int> r;
    if(v<0||v>=N||!aliveV[v]) return r;
    r.reserve(VF[v].size());
    for(int id:VF[v]) if(id>=0 && id<(int)F.size() && aliveF[id] && hasv(F[id],v)) r.push_back(id);
    sort(r.begin(),r.end());
    r.erase(unique(r.begin(),r.end()),r.end());
    return r;
}

static void rebuildVF(){
    VF.assign(N,{});
    for(int i=0;i<(int)F.size();++i) if(aliveF[i]){
        VF[F[i].a].push_back(i);
        VF[F[i].b].push_back(i);
        VF[F[i].c].push_back(i);
    }
}

static int edgeCountAlive(int a,int b,const vector<unsigned char>*ignore=nullptr){
    int cnt=0;
    const vector<int>& small = (VF[a].size()<VF[b].size()?VF[a]:VF[b]);
    for(int id:small){
        if(id<0||id>=(int)F.size()||!aliveF[id]) continue;
        if(ignore && (*ignore)[id]) continue;
        if(hasv(F[id],a)&&hasv(F[id],b)) ++cnt;
    }
    return cnt;
}

static bool duplicateAliveFace(int a,int b,int c,const vector<unsigned char>&ignore){
    auto K=fkey(a,b,c);
    int vv[3]={a,b,c};
    int best=vv[0];
    if(VF[vv[1]].size()<VF[best].size()) best=vv[1];
    if(VF[vv[2]].size()<VF[best].size()) best=vv[2];
    for(int id:VF[best]){
        if(id<0||id>=(int)F.size()||!aliveF[id]||ignore[id]) continue;
        const Face&g=F[id];
        if(fkey(g.a,g.b,g.c)==K) return true;
    }
    return false;
}

struct Eval{
    bool ok=false;
    double cost=1e100;
    int v=-1;
    vector<int> poly;
    vector<Face> add;
};

static bool buildLinkCycle(int v, const vector<int>&inc, vector<int>&poly){
    unordered_map<int, vector<int>> adj;
    adj.reserve(inc.size()*2+4);
    for(int id:inc){
        Face f=F[id];
        int other[2], k=0;
        int vv[3]={f.a,f.b,f.c};
        for(int i=0;i<3;i++) if(vv[i]!=v) other[k++]=vv[i];
        if(k!=2 || other[0]==other[1]) return false;
        adj[other[0]].push_back(other[1]);
        adj[other[1]].push_back(other[0]);
    }
    if(adj.size()!=inc.size()) return false;
    for(auto &p:adj) if(p.second.size()!=2) return false;
    int start=adj.begin()->first, prev=-1, cur=start;
    poly.clear(); poly.reserve(adj.size());
    unordered_set<int> seen; seen.reserve(adj.size()+2);
    for(size_t step=0; step<adj.size(); ++step){
        if(seen.count(cur)) return false;
        seen.insert(cur);
        poly.push_back(cur);
        const auto &nb=adj[cur];
        int nxt = (nb[0]==prev ? nb[1] : nb[0]);
        prev=cur; cur=nxt;
    }
    return cur==start;
}

static Eval evalVertex(int v, double cosMin, int maxValence){
    Eval e; e.v=v;
    if(!aliveV[v]) return e;
    vector<int> inc=incident(v);
    int k=(int)inc.size();
    if(k<3 || k>maxValence) return e;
    vector<int> poly;
    if(!buildLinkCycle(v,inc,poly)) return e;
    if((int)poly.size()!=k) return e;

    double nearest=1e100;
    for(int u:poly) nearest=min(nearest,norm3(P[v]-P[u]));
    if(nearest>deleteRadius) return e;

    vector<unsigned char> ignore(F.size(),0);
    Vec3 oldSum{0,0,0};
    double oldArea=0.0;
    for(int id:inc){
        ignore[id]=1;
        Vec3 cr=faceCross(F[id]);
        double a=norm3(cr);
        if(a<=epsArea) return e;
        oldSum=oldSum+cr;
        oldArea+=a;
    }
    double oldNorm=norm3(oldSum);
    if(oldNorm<=epsArea) return e;
    Vec3 avg=oldSum*(1.0/oldNorm);

    for(int id:inc){
        Vec3 cr=faceCross(F[id]);
        double a=norm3(cr);
        double d=dot3(cr,avg)/a;
        if(d<cosMin) return e;
    }

    auto makeFan = [&](const vector<int>&order)->vector<Face>{
        vector<Face> tris;
        for(int i=1;i+1<(int)order.size();++i) tris.push_back({order[0],order[i],order[i+1]});
        return tris;
    };
    vector<Face> tris=makeFan(poly);
    Vec3 fanSum{0,0,0};
    for(const Face&t:tris) fanSum=fanSum+faceCross(t);
    if(dot3(fanSum,avg)<0){
        reverse(poly.begin(),poly.end());
        tris=makeFan(poly);
    }

    unordered_map<uint64_t,int> newEdgeCnt;
    newEdgeCnt.reserve(tris.size()*3+4);
    for(const Face&t:tris){
        int vv[3]={t.a,t.b,t.c};
        if(vv[0]==vv[1]||vv[0]==vv[2]||vv[1]==vv[2]) return e;
        Vec3 cr=faceCross(t);
        double ar=norm3(cr);
        if(ar<=epsArea) return e;
        double nd=dot3(cr,avg)/ar;
        if(nd<cosMin) return e;
        if(duplicateAliveFace(t.a,t.b,t.c,ignore)) return e;
        for(int j=0;j<3;j++) newEdgeCnt[ekey(vv[j],vv[(j+1)%3])]++;
    }

    for(int i=0;i<k;i++){
        int a=poly[i], b=poly[(i+1)%k];
        uint64_t ke=ekey(a,b);
        if(newEdgeCnt[ke]!=1) return e;
        if(edgeCountAlive(a,b,&ignore)!=1) return e;
    }
    for(auto &p:newEdgeCnt){
        int cnt=p.second;
        if(cnt==1) continue;
        if(cnt!=2) return e;
        int a=(int)(p.first>>32), b=(int)(p.first&0xffffffffu);
        if(edgeCountAlive(a,b,&ignore)!=0) return e;
    }

    double newArea=0.0, worstN=0.0;
    for(const Face&t:tris){
        Vec3 cr=faceCross(t);
        double ar=norm3(cr);
        newArea+=ar;
        worstN=max(worstN,1.0-dot3(cr,avg)/ar);
    }
    double areaRel=fabs(newArea-oldArea)/(oldArea+1e-30);

    e.cost = 0.55*(nearest/(deleteRadius+1e-30)) + 18.0*worstN + 0.40*areaRel + 0.015*k;
    e.poly=poly;
    e.add=tris;
    e.ok=true;
    return e;
}

struct Cand{
    double cost;
    int v;
    unsigned st;
    bool operator<(const Cand&o)const{return cost>o.cost;}
};

static priority_queue<Cand> pq;

static void pushVertex(int v,double cosMin,int maxValence){
    if(!aliveV[v]) return;
    Eval e=evalVertex(v,cosMin,maxValence);
    if(e.ok) pq.push({e.cost,v,stampV[v]});
}

static void applyEval(const Eval&e){
    int v=e.v;
    if(!aliveV[v]) return;
    vector<int> inc=incident(v);
    for(int id:inc) if(aliveF[id]){
        aliveF[id]=0;
        --liveF;
    }
    aliveV[v]=0;
    --liveV;

    for(const Face&t:e.add){
        int id=(int)F.size();
        F.push_back(t);
        aliveF.push_back(1);
        ++liveF;
        VF[t.a].push_back(id);
        VF[t.b].push_back(id);
        VF[t.c].push_back(id);
    }
    VF[v].clear();
    ++stampV[v];
    for(int u:e.poly) ++stampV[u];
}

static bool closedManifoldCheck(){
    unordered_map<uint64_t,int> ec;
    ec.reserve((size_t)liveF*3+4);
    for(int id=0; id<(int)F.size(); ++id) if(aliveF[id]){
        const Face&f=F[id];
        if(!aliveV[f.a]||!aliveV[f.b]||!aliveV[f.c]) return false;
        if(f.a==f.b||f.a==f.c||f.b==f.c) return false;
        if(norm3(faceCross(f))<=epsArea) return false;
        int vv[3]={f.a,f.b,f.c};
        for(int i=0;i<3;i++) ec[ekey(vv[i],vv[(i+1)%3])]++;
    }
    for(auto &p:ec) if(p.second!=2) return false;

    for(int v=0; v<N; ++v) if(aliveV[v]){
        vector<int> in=incident(v);
        if(in.empty()) return false;
        unordered_map<int, vector<int>> adj;
        for(int id:in){
            int a=-1,b=-1;
            int vv[3]={F[id].a,F[id].b,F[id].c};
            for(int i=0;i<3;i++) if(vv[i]!=v){ if(a<0) a=vv[i]; else b=vv[i]; }
            if(a<0||b<0||a==b) return false;
            adj[a].push_back(b); adj[b].push_back(a);
        }
        for(auto &q:adj) if(q.second.size()!=2) return false;
        int start=adj.begin()->first, prev=-1, cur=start, seen=0;
        unordered_set<int> vis; vis.reserve(adj.size()+2);
        while(true){
            if(vis.count(cur)) break;
            vis.insert(cur); ++seen;
            auto &nb=adj[cur];
            int nxt=(nb[0]==prev?nb[1]:nb[0]);
            prev=cur; cur=nxt;
        }
        if(cur!=start || seen!=(int)adj.size()) return false;
    }
    return true;
}

static void runPass(double cosMin,int maxValence,int targetV,int opLimit){
    while(!pq.empty()) pq.pop();
    for(int v=0; v<N; ++v) if(aliveV[v]) pushVertex(v,cosMin,maxValence);
    int ops=0;
    while(!pq.empty() && liveV>targetV && ops<opLimit){
        Cand c=pq.top(); pq.pop();
        if(!aliveV[c.v]||stampV[c.v]!=c.st) continue;
        Eval e=evalVertex(c.v,cosMin,maxValence);
        if(!e.ok) continue;
        if(e.cost>c.cost*1.35+1e-12){
            pushVertex(c.v,cosMin,maxValence);
            continue;
        }
        applyEval(e);
        ++ops;
        for(int u:e.poly){
            pushVertex(u,cosMin,maxValence);
            vector<int> iu=incident(u);
            for(int id:iu){
                int vv[3]={F[id].a,F[id].b,F[id].c};
                for(int j=0;j<3;j++) pushVertex(vv[j],cosMin,maxValence);
            }
        }
        if((ops&511)==0) rebuildVF();
    }
    rebuildVF();
}

static void emit(){
    vector<int> id(N,-1);
    vector<Vec3> outV; outV.reserve(liveV);
    for(int i=0;i<N;i++) if(aliveV[i]){
        id[i]=(int)outV.size()+1;
        outV.push_back(P[i]);
    }
    vector<Face> outF; outF.reserve(liveF);
    for(int i=0;i<(int)F.size();++i) if(aliveF[i]){
        Face f=F[i];
        if(id[f.a]>0&&id[f.b]>0&&id[f.c]>0&&id[f.a]!=id[f.b]&&id[f.a]!=id[f.c]&&id[f.b]!=id[f.c])
            outF.push_back({id[f.a],id[f.b],id[f.c]});
    }
    cout.setf(ios::fixed);
    cout.precision(10);
    cout << outV.size() << ' ' << outF.size() << '\n';
    for(const Vec3&p:outV) cout << "v " << p.x << ' ' << p.y << ' ' << p.z << '\n';
    for(const Face&f:outF) cout << "f " << f.a << ' ' << f.b << ' ' << f.c << '\n';
}

int main(){
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    if(!(cin>>N>>M0)) return 0;
    P.resize(N);
    Vec3 mn{1e100,1e100,1e100}, mx{-1e100,-1e100,-1e100};
    string tag;
    for(int i=0;i<N;i++){
        cin >> tag >> P[i].x >> P[i].y >> P[i].z;
        mn.x=min(mn.x,P[i].x); mn.y=min(mn.y,P[i].y); mn.z=min(mn.z,P[i].z);
        mx.x=max(mx.x,P[i].x); mx.y=max(mx.y,P[i].y); mx.z=max(mx.z,P[i].z);
    }
    F.resize(M0);
    for(int i=0;i<M0;i++){
        cin >> tag >> F[i].a >> F[i].b >> F[i].c;
        --F[i].a; --F[i].b; --F[i].c;
    }

    Vec3 box=mx-mn;
    diagLen=norm3(box);
    if(!(diagLen>0)) diagLen=1.0;
    epsArea=1e-20*diagLen*diagLen;
    deleteRadius=0.0490*diagLen;

    aliveV.assign(N,1);
    aliveF.assign(F.size(),1);
    stampV.assign(N,1);
    liveV=N;
    liveF=(int)F.size();
    rebuildVF();

    if(N<=8){
        emit();
        return 0;
    }

    int t1=max(8,(int)ceil(N*0.55));
    int t2=max(8,(int)ceil(N*0.34));
    int t3=max(8,(int)ceil(N*0.22));
    int t4=max(8,(int)ceil(N*0.15));

    runPass(0.985,4,t1,N);
    runPass(0.965,5,t2,N);
    runPass(0.940,6,t3,2*N);
    runPass(0.900,7,t4,3*N);

    if(!closedManifoldCheck()){
        aliveV.assign(N,1);
        aliveF.assign(F.size(),0);
        for(int i=0;i<M0;i++) aliveF[i]=1;
        liveV=N; liveF=M0;
        F.resize(M0);
        rebuildVF();
    }

    emit();
    return 0;
}