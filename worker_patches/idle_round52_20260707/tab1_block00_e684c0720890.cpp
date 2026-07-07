#include <bits/stdc++.h>
using namespace std;

struct Vec3{
    double x,y,z;
};
static inline Vec3 operator+(const Vec3&a,const Vec3&b){return {a.x+b.x,a.y+b.y,a.z+b.z};}
static inline Vec3 operator-(const Vec3&a,const Vec3&b){return {a.x-b.x,a.y-b.y,a.z-b.z};}
static inline Vec3 operator*(const Vec3&a,double s){return {a.x*s,a.y*s,a.z*s};}
static inline double dotv(const Vec3&a,const Vec3&b){return a.x*b.x+a.y*b.y+a.z*b.z;}
static inline Vec3 crossv(const Vec3&a,const Vec3&b){return {a.y*b.z-a.z*b.y,a.z*b.x-a.x*b.z,a.x*b.y-a.y*b.x};}
static inline double n2(const Vec3&a){return dotv(a,a);} 
static inline double norm(const Vec3&a){return sqrt(n2(a));}

struct Face{int a,b,c;};
static int N,M;
static vector<Vec3> P, Orig;
static vector<Face> F;
static vector<unsigned char> aliveV, aliveF;
static vector<double> errV;
static vector<vector<int>> inc;
static int aliveFaces;
static double diagLen, EPS;
static chrono::steady_clock::time_point T0;

static inline bool hasv(const Face&f,int v){return f.a==v||f.b==v||f.c==v;}
static inline int thirdv(const Face&f,int u,int v){
    if(f.a!=u&&f.a!=v) return f.a;
    if(f.b!=u&&f.b!=v) return f.b;
    if(f.c!=u&&f.c!=v) return f.c;
    return -1;
}
static inline array<int,3> keytri(int a,int b,int c){array<int,3> r={a,b,c}; sort(r.begin(),r.end()); return r;}
static inline long long edgekey(int a,int b){ if(a>b) swap(a,b); return (long long)(unsigned)a<<32 | (unsigned)b; }
static inline double elapsed(){return chrono::duration<double>(chrono::steady_clock::now()-T0).count();}
static inline Vec3 fnormal(int fid){const Face&f=F[fid]; return crossv(P[f.b]-P[f.a],P[f.c]-P[f.a]);}

static void read_input(){
    ios::sync_with_stdio(false); cin.tie(nullptr);
    cin>>N>>M;
    P.resize(N); Orig.resize(N);
    Vec3 mn{1e100,1e100,1e100}, mx{-1e100,-1e100,-1e100};
    for(int i=0;i<N;i++){
        char ch; cin>>ch>>P[i].x>>P[i].y>>P[i].z; Orig[i]=P[i];
        mn.x=min(mn.x,P[i].x); mn.y=min(mn.y,P[i].y); mn.z=min(mn.z,P[i].z);
        mx.x=max(mx.x,P[i].x); mx.y=max(mx.y,P[i].y); mx.z=max(mx.z,P[i].z);
    }
    F.resize(M); inc.assign(N,{});
    for(int i=0;i<M;i++){
        char ch; int a,b,c; cin>>ch>>a>>b>>c; --a;--b;--c;
        F[i]={a,b,c};
        inc[a].push_back(i); inc[b].push_back(i); inc[c].push_back(i);
    }
    aliveV.assign(N,1); aliveF.assign(M,1); errV.assign(N,0.0); aliveFaces=M;
    diagLen=norm(mx-mn); if(!(diagLen>0)) diagLen=1.0;
    EPS=0.05*diagLen;
}

static void clean_incident(int v){
    auto &q=inc[v]; int w=0;
    for(int id:q) if(id>=0 && id<M && aliveF[id] && hasv(F[id],v)) q[w++]=id;
    q.resize(w);
}

static bool edge_faces(int u,int v,int ef[2],int op[2]){
    clean_incident(u);
    int cnt=0;
    for(int fid:inc[u]) if(aliveF[fid] && hasv(F[fid],v)){
        if(cnt>=2) return false;
        ef[cnt]=fid; op[cnt]=thirdv(F[fid],u,v); cnt++;
    }
    if(cnt!=2) return false;
    if(op[0]<0||op[1]<0||op[0]==op[1]) return false;
    return true;
}

static bool link_ok(int u,int v,int op0,int op1){
    static vector<int> mark; static int stamp=1;
    if((int)mark.size()<N) mark.assign(N,0);
    if(++stamp==INT_MAX){fill(mark.begin(),mark.end(),0); stamp=1;}
    clean_incident(u); clean_incident(v);
    for(int fid:inc[u]) if(aliveF[fid]){
        Face f=F[fid]; int a[3]={f.a,f.b,f.c};
        for(int x:a) if(x!=u && x!=v) mark[x]=stamp;
    }
    int common=0; bool s0=false,s1=false;
    for(int fid:inc[v]) if(aliveF[fid]){
        Face f=F[fid]; int a[3]={f.a,f.b,f.c};
        for(int x:a) if(x!=u && x!=v && mark[x]==stamp){
            if(x==op0) s0=true; else if(x==op1) s1=true; else return false;
        }
    }
    common = (s0?1:0)+(s1?1:0);
    return common==2;
}

static bool duplicate_after(int keep,int rem,int fid,int a,int b,int c,int e0,int e1){
    if(a==b||a==c||b==c) return true;
    auto want=keytri(a,b,c);
    int base=a;
    if(inc[b].size()<inc[base].size()) base=b;
    if(inc[c].size()<inc[base].size()) base=c;
    clean_incident(base);
    for(int g:inc[base]) if(aliveF[g] && g!=fid && g!=e0 && g!=e1){
        if(hasv(F[g],rem)) continue;
        if(keytri(F[g].a,F[g].b,F[g].c)==want) return true;
    }
    return false;
}

struct Eval{bool ok=false; double cost=1e100, newErr=0;};

static Eval eval_collapse(int keep,int rem,const int ef[2],double emax,double ndotMin,double planeMax,double areaRelMin){
    Eval r;
    double d=norm(P[keep]-P[rem]);
    r.newErr=max(errV[keep],errV[rem]+d);
    if(r.newErr>emax) return r;
    clean_incident(rem);
    double worstPlane=0, worstFlip=0, areaLoss=0;
    int touched=0;
    for(int fid:inc[rem]) if(aliveF[fid] && fid!=ef[0] && fid!=ef[1]){
        Face old=F[fid];
        int a=old.a,b=old.b,c=old.c;
        if(a==rem) a=keep; if(b==rem) b=keep; if(c==rem) c=keep;
        if(duplicate_after(keep,rem,fid,a,b,c,ef[0],ef[1])) return r;
        Vec3 nOld=crossv(P[old.b]-P[old.a],P[old.c]-P[old.a]);
        Vec3 nNew=crossv(P[b]-P[a],P[c]-P[a]);
        double lo=norm(nOld), ln=norm(nNew);
        if(!(lo>1e-30)||!(ln>1e-30)) return r;
        if(ln < lo*areaRelMin) return r;
        double co=dotv(nOld,nNew)/(lo*ln);
        if(co<ndotMin) return r;
        Vec3 un=nOld*(1.0/lo);
        double pd=fabs(dotv(un,P[keep]-P[old.a]));
        if(pd>planeMax) return r;
        worstPlane=max(worstPlane,pd);
        worstFlip=max(worstFlip,1.0-co);
        areaLoss=max(areaLoss,max(0.0,1.0-ln/lo));
        touched++;
    }
    if(touched==0) return r;
    r.ok=true;
    double l=d/(EPS+1e-30), er=r.newErr/(EPS+1e-30), pl=worstPlane/(EPS+1e-30);
    r.cost=0.70*l+1.60*er+2.5*pl+180.0*worstFlip+0.07*touched+0.04*areaLoss;
    return r;
}

static bool do_collapse(int keep,int rem,const int ef[2],double newErr){
    if(!aliveV[keep]||!aliveV[rem]) return false;
    for(int k=0;k<2;k++) if(aliveF[ef[k]]){ aliveF[ef[k]]=0; aliveFaces--; }
    clean_incident(rem);
    vector<int> moved=inc[rem];
    for(int fid:moved) if(aliveF[fid] && hasv(F[fid],rem)){
        Face &f=F[fid];
        if(f.a==rem) f.a=keep; if(f.b==rem) f.b=keep; if(f.c==rem) f.c=keep;
        inc[keep].push_back(fid);
    }
    aliveV[rem]=0; errV[keep]=newErr; inc[rem].clear();
    if(inc[keep].size()>512) clean_incident(keep);
    return true;
}

struct Cand{double cost; int u,v;};
struct Cmp{bool operator()(const Cand&a,const Cand&b)const{return a.cost>b.cost;}};

static void push_edges_of_vertex(int v, priority_queue<Cand,vector<Cand>,Cmp>&pq, unordered_set<long long>&seen, double scale){
    if(!aliveV[v]) return;
    clean_incident(v);
    for(int fid:inc[v]) if(aliveF[fid]){
        Face f=F[fid]; int a[3]={f.a,f.b,f.c};
        for(int i=0;i<3;i++){
            int x=a[i], y=a[(i+1)%3];
            if(!aliveV[x]||!aliveV[y]||x==y) continue;
            long long k=edgekey(x,y);
            if(seen.insert(k).second){
                double len=norm(P[x]-P[y]);
                pq.push({len*scale + 1e-9*(inc[x].size()+inc[y].size()), x, y});
            }
        }
    }
}

static void simplify_pass(double efrac,double ndotMin,double planeFrac,double areaRelMin,double timelimit){
    double emax=EPS*efrac, planeMax=EPS*planeFrac;
    priority_queue<Cand,vector<Cand>,Cmp> pq;
    unordered_set<long long> seen; seen.reserve((size_t)M*3/2+10);
    for(int i=0;i<M;i++) if(aliveF[i]){
        int a=F[i].a,b=F[i].b,c=F[i].c;
        long long keys[3]={edgekey(a,b),edgekey(b,c),edgekey(c,a)};
        int xs[3]={a,b,c}, ys[3]={b,c,a};
        for(int e=0;e<3;e++) if(seen.insert(keys[e]).second){
            double len=norm(P[xs[e]]-P[ys[e]]);
            pq.push({len,xs[e],ys[e]});
        }
    }
    int ops=0;
    while(!pq.empty() && elapsed()<timelimit){
        Cand c=pq.top(); pq.pop();
        int u=c.u,v=c.v;
        if(u<0||v<0||u>=N||v>=N||!aliveV[u]||!aliveV[v]) continue;
        int ef[2]={-1,-1}, op[2]={-1,-1};
        if(!edge_faces(u,v,ef,op)) continue;
        if(!link_ok(u,v,op[0],op[1])) continue;
        Eval eu=eval_collapse(u,v,ef,emax,ndotMin,planeMax,areaRelMin);
        Eval ev=eval_collapse(v,u,ef,emax,ndotMin,planeMax,areaRelMin);
        if(!eu.ok && !ev.ok) continue;
        int keep=u,rem=v; Eval best=eu;
        if(ev.ok && (!eu.ok || ev.cost<eu.cost)){keep=v; rem=u; best=ev;}
        if(!do_collapse(keep,rem,ef,best.newErr)) continue;
        unordered_set<long long> loc; loc.reserve(256);
        push_edges_of_vertex(keep,pq,loc,0.85);
        if(++ops%512==0) clean_incident(keep);
    }
}

static bool verify_closed(){
    unordered_map<long long,int> cnt; cnt.reserve((size_t)aliveFaces*3+10);
    for(int i=0;i<M;i++) if(aliveF[i]){
        Face f=F[i];
        if(f.a==f.b||f.a==f.c||f.b==f.c) return false;
        if(!aliveV[f.a]||!aliveV[f.b]||!aliveV[f.c]) return false;
        if(norm(crossv(P[f.b]-P[f.a],P[f.c]-P[f.a]))<=1e-30) return false;
        cnt[edgekey(f.a,f.b)]++; cnt[edgekey(f.b,f.c)]++; cnt[edgekey(f.c,f.a)]++;
    }
    for(auto &kv:cnt) if(kv.second!=2) return false;
    return true;
}

static void output_mesh_original(){
    cout.setf(ios::fixed); cout<<setprecision(10);
    cout<<N<<" "<<M<<"\n";
    for(auto &p:Orig) cout<<"v "<<p.x<<" "<<p.y<<" "<<p.z<<"\n";
    for(auto &f:F) cout<<"f "<<f.a+1<<" "<<f.b+1<<" "<<f.c+1<<"\n";
}

static void output_mesh(){
    if(!verify_closed()){ output_mesh_original(); return; }
    vector<int> id(N,-1); int nv=0,nf=0;
    for(int i=0;i<N;i++) if(aliveV[i]) id[i]=++nv;
    for(int i=0;i<M;i++) if(aliveF[i]) nf++;
    cout.setf(ios::fixed); cout<<setprecision(10);
    cout<<nv<<" "<<nf<<"\n";
    for(int i=0;i<N;i++) if(aliveV[i]) cout<<"v "<<P[i].x<<" "<<P[i].y<<" "<<P[i].z<<"\n";
    for(int i=0;i<M;i++) if(aliveF[i]) cout<<"f "<<id[F[i].a]<<" "<<id[F[i].b]<<" "<<id[F[i].c]<<"\n";
}

int main(){
    T0=chrono::steady_clock::now();
    read_input();
