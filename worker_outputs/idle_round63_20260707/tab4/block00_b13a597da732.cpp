#include <bits/stdc++.h>
using namespace std;

struct Vec3{double x,y,z;};
static inline Vec3 operator+(const Vec3&a,const Vec3&b){return {a.x+b.x,a.y+b.y,a.z+b.z};}
static inline Vec3 operator-(const Vec3&a,const Vec3&b){return {a.x-b.x,a.y-b.y,a.z-b.z};}
static inline Vec3 operator*(const Vec3&a,double s){return {a.x*s,a.y*s,a.z*s};}
static inline double dotp(const Vec3&a,const Vec3&b){return a.x*b.x+a.y*b.y+a.z*b.z;}
static inline Vec3 crossp(const Vec3&a,const Vec3&b){return {a.y*b.z-a.z*b.y,a.z*b.x-a.x*b.z,a.x*b.y-a.y*b.x};}
static inline double n2(const Vec3&a){return dotp(a,a);} 
static inline double norm(const Vec3&a){return sqrt(n2(a));}

struct Face{int a,b,c;};
struct EdgeCand{int u,v; double cost;};

int N,M;
vector<Vec3> P, Pin;
vector<Face> F;
vector<unsigned char> aliveV, aliveF;
vector<double> radv;
vector<vector<int>> vfaces;
double diagLen, EPS;
chrono::steady_clock::time_point T0;

static inline double elapsed(){return chrono::duration<double>(chrono::steady_clock::now()-T0).count();}
static inline unsigned long long ekey(int a,int b){ if(a>b) swap(a,b); return (unsigned long long)(unsigned int)a<<32 | (unsigned int)b; }
static inline array<int,3> skey(int a,int b,int c){ array<int,3>s={a,b,c}; sort(s.begin(),s.end()); return s; }
static inline bool hasv(const Face&f,int x){return f.a==x||f.b==x||f.c==x;}
static inline int thirdv(const Face&f,int u,int v){ if(f.a!=u&&f.a!=v) return f.a; if(f.b!=u&&f.b!=v) return f.b; return f.c; }
static inline Vec3 faceCross(const Face&f, const vector<Vec3>&Q){return crossp(Q[f.b]-Q[f.a], Q[f.c]-Q[f.a]);}

void readInput(){
    ios::sync_with_stdio(false); cin.tie(nullptr);
    cin>>N>>M; P.resize(N); Pin.resize(N); aliveV.assign(N,1); radv.assign(N,0.0);
    Vec3 mn{1e100,1e100,1e100}, mx{-1e100,-1e100,-1e100};
    for(int i=0;i<N;i++){ char ch; cin>>ch>>P[i].x>>P[i].y>>P[i].z; Pin[i]=P[i];
        mn.x=min(mn.x,P[i].x); mn.y=min(mn.y,P[i].y); mn.z=min(mn.z,P[i].z);
        mx.x=max(mx.x,P[i].x); mx.y=max(mx.y,P[i].y); mx.z=max(mx.z,P[i].z);
    }
    diagLen=norm(mx-mn); if(!(diagLen>0)) diagLen=1.0; EPS=0.05*diagLen;
    F.resize(M); aliveF.assign(M,1);
    for(int i=0;i<M;i++){ char ch; cin>>ch; cin>>F[i].a>>F[i].b>>F[i].c; --F[i].a;--F[i].b;--F[i].c; }
}

void rebuildVFaces(){
    vfaces.assign(N,{});
    vector<int> deg(N,0);
    for(int i=0;i<(int)F.size();i++) if(aliveF[i]){deg[F[i].a]++;deg[F[i].b]++;deg[F[i].c]++;}
    for(int i=0;i<N;i++) if(aliveV[i]) vfaces[i].reserve(deg[i]);
    for(int i=0;i<(int)F.size();i++) if(aliveF[i]){vfaces[F[i].a].push_back(i);vfaces[F[i].b].push_back(i);vfaces[F[i].c].push_back(i);}    
}

bool edgeFaces(int u,int v,int ef[2],int op[2]){
    int cnt=0;
    for(int id: vfaces[u]) if(aliveF[id] && hasv(F[id],v)){
        if(cnt>=2) return false;
        ef[cnt]=id; op[cnt]=thirdv(F[id],u,v); cnt++;
    }
    return cnt==2 && op[0]!=op[1] && op[0]>=0 && op[1]>=0;
}

bool linkOK(int u,int v,int op0,int op1){
    static vector<int> mark; static int stamp=1; static vector<int> seen;
    if((int)mark.size()!=N) mark.assign(N,0);
    if(++stamp==INT_MAX){fill(mark.begin(),mark.end(),0); stamp=1;}
    for(int id:vfaces[u]) if(aliveF[id]){
        const Face&f=F[id]; int xs[3]={f.a,f.b,f.c};
        for(int x:xs) if(x!=u && x!=v) mark[x]=stamp;
    }
    int got0=0,got1=0;
    for(int id:vfaces[v]) if(aliveF[id]){
        const Face&f=F[id]; int xs[3]={f.a,f.b,f.c};
        for(int x:xs) if(x!=u && x!=v && mark[x]==stamp){
            if(x==op0) got0=1; else if(x==op1) got1=1; else return false;
        }
    }
    return got0&&got1;
}

bool wouldDuplicate(int keep,int rem,int fidSkip0,int fidSkip1,int a,int b,int c){
    auto target=skey(a,b,c);
    int probe=keep;
    for(int id:vfaces[probe]) if(aliveF[id] && id!=fidSkip0 && id!=fidSkip1){
        Face g=F[id];
        if(g.a==rem) g.a=keep; if(g.b==rem) g.b=keep; if(g.c==rem) g.c=keep;
        if(g.a==g.b||g.a==g.c||g.b==g.c) continue;
        if(skey(g.a,g.b,g.c)==target) return true;
    }
    return false;
}

struct Trial{bool ok=false; int keep=-1,rem=-1; Vec3 pos; double nr=0,cost=1e100;};

Trial testCollapse(int keep,int rem,const Vec3&pos,int ef0,int ef1,double minDot,double areaFrac){
    Trial tr; tr.keep=keep; tr.rem=rem; tr.pos=pos;
    double nr=max(radv[keep]+norm(pos-P[keep]), radv[rem]+norm(pos-P[rem]));
    if(nr>EPS*0.985) return tr;
    vector<int> ids; ids.reserve(vfaces[keep].size()+vfaces[rem].size());
    for(int id:vfaces[keep]) if(aliveF[id] && id!=ef0 && id!=ef1) ids.push_back(id);
    for(int id:vfaces[rem]) if(aliveF[id] && id!=ef0 && id!=ef1) ids.push_back(id);
    sort(ids.begin(),ids.end()); ids.erase(unique(ids.begin(),ids.end()),ids.end());
    for(int id:ids){
        Face old=F[id], nw=old;
        if(nw.a==rem) nw.a=keep; if(nw.b==rem) nw.b=keep; if(nw.c==rem) nw.c=keep;
        if(nw.a==nw.b||nw.a==nw.c||nw.b==nw.c) return tr;
        Vec3 co=faceCross(old,P); double ao=norm(co); if(!(ao>1e-20)) return tr;
        Vec3 pa=(nw.a==keep?pos:P[nw.a]); Vec3 pb=(nw.b==keep?pos:P[nw.b]); Vec3 pc=(nw.c==keep?pos:P[nw.c]);
        Vec3 cn=crossp(pb-pa,pc-pa); double an=norm(cn); if(!(an>1e-20)) return tr;
        double d=dotp(co,cn)/(ao*an); if(d<minDot) return tr;
        if(an<ao*areaFrac) return tr;
        if(hasv(old,rem) && wouldDuplicate(keep,rem,ef0,ef1,nw.a,nw.b,nw.c)) return tr;
    }
    tr.ok=true; tr.nr=nr; tr.cost=nr/EPS + 0.015*ids.size(); return tr;
}

bool collapseEdge(int u,int v,double minDot,double areaFrac){
    if(u==v||!aliveV[u]||!aliveV[v]) return false;
    int ef[2]={-1,-1},op[2]={-1,-1};
    if(!edgeFaces(u,v,ef,op)) return false;
    if(!linkOK(u,v,op[0],op[1])) return false;
    vector<Trial> cand;
    Vec3 mid=(P[u]+P[v])*0.5;
    cand.push_back(testCollapse(u,v,mid,ef[0],ef[1],minDot,areaFrac));
    cand.push_back(testCollapse(v,u,mid,ef[0],ef[1],minDot,areaFrac));
    cand.push_back(testCollapse(u,v,P[u],ef[0],ef[1],minDot,areaFrac));
    cand.push_back(testCollapse(v,u,P[v],ef[0],ef[1],minDot,areaFrac));
    int bi=-1; double bc=1e100;
    for(int i=0;i<(int)cand.size();i++) if(cand[i].ok && cand[i].cost<bc){bc=cand[i].cost;bi=i;}
    if(bi<0) return false;
    Trial t=cand[bi]; int keep=t.keep, rem=t.rem;
    aliveF[ef[0]]=aliveF[ef[1]]=0;
    for(int id:vfaces[rem]) if(aliveF[id]){
        if(F[id].a==rem) F[id].a=keep; if(F[id].b==rem) F[id].b=keep; if(F[id].c==rem) F[id].c=keep;
        if(F[id].a==F[id].b||F[id].a==F[id].c||F[id].b==F[id].c) aliveF[id]=0;
        else vfaces[keep].push_back(id);
    }
    aliveV[rem]=0; P[keep]=t.pos; radv[keep]=t.nr; vector<int>().swap(vfaces[rem]);
    return true;
}

vector<EdgeCand> collectEdges(){
    unordered_set<unsigned long long> seen; seen.reserve(F.size()*2+1);
    vector<EdgeCand> e; e.reserve(F.size()*3/2+1);
    for(int i=0;i<(int)F.size();i++) if(aliveF[i]){
        int a[3]={F[i].a,F[i].b,F[i].c};
        for(int k=0;k<3;k++){
            int u=a[k],v=a[(k+1)%3]; if(!aliveV[u]||!aliveV[v]) continue;
            auto key=ekey(u,v); if(seen.insert(key).second){
                double L=norm(P[u]-P[v]);
                if(L<=2.0*EPS) e.push_back({u,v,L + 0.35*(radv[u]+radv[v])});
            }
        }
    }
    sort(e.begin(),e.end(),[](const EdgeCand&a,const EdgeCand&b){return a.cost<b.cost;});
    return e;
}

void simplify(){
    if(N<=8) return;
    rebuildVFaces();
    const double mdots[]={0.995,0.985,0.965,0.94,0.90,0.86,0.80};
    const double afracs[]={0.25,0.18,0.12,0.08,0.04,0.02,0.005};
    int phaseCount=7;
    int aliveCount=N;
    for(int ph=0; ph<phaseCount && elapsed()<8.2; ph++){
        bool progress=true; int loops=0;
        while(progress && loops<4 && elapsed()<8.2){
            loops++; progress=false;
            auto edges=collectEdges();
            int done=0, tried=0;
            for(const auto &ec:edges){
                if(elapsed()>8.4) break;
                tried++;
                if(collapseEdge(ec.u,ec.v,mdots[ph],afracs[ph])){ progress=true; done++; aliveCount--; if(done%256==0) rebuildVFaces(); }
                if(tried>250000 && done<16) break;
            }
            rebuildVFaces();
            if(done<8) break;
        }
    }
}

void outputMesh(){
    unordered_set<string> used; used.reserve(F.size()*2+1);
    vector<Face> outF; outF.reserve(F.size());
    vector<int> usedV(N,0);
    for(int i=0;i<(int)F.size();i++) if(aliveF[i]){
        Face f=F[i]; if(!aliveV[f.a]||!aliveV[f.b]||!aliveV[f.c]) continue;
        if(f.a==f.b||f.a==f.c||f.b==f.c) continue;
        if(norm(faceCross(f,P))<=1e-18*diagLen*diagLen) continue;
        auto s=skey(f.a,f.b,f.c);
        string key=to_string(s[0])+","+to_string(s[1])+","+to_string(s[2]);
        if(used.insert(key).second){ outF.push_back(f); usedV[f.a]=usedV[f.b]=usedV[f.c]=1; }
    }
    vector<int> id(N,-1); vector<Vec3> outP; outP.reserve(N);
    for(int i=0;i<N;i++) if(aliveV[i]&&usedV[i]){ id[i]=(int)outP.size()+1; outP.push_back(P[i]); }
    cout.setf(ios::fixed); cout<<setprecision(10);
    cout<<outP.size()<<" "<<outF.size()<<"\n";
    for(auto &p:outP) cout<<"v "<<p.x<<" "<<p.y<<" "<<p.z<<"\n";
    for(auto &f:outF) cout<<"f "<<id[f.a]<<" "<<id[f.b]<<" "<<id[f.c]<<"\n";
}

int main(){
    T0=chrono::steady_clock::now();
    readInput();
    simplify();
    outputMesh();
    return 0;
}