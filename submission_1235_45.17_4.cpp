
#include <bits/stdc++.h>
using namespace std;

struct Vec3{double x,y,z;};
static inline Vec3 operator+(const Vec3&a,const Vec3&b){return {a.x+b.x,a.y+b.y,a.z+b.z};}
static inline Vec3 operator-(const Vec3&a,const Vec3&b){return {a.x-b.x,a.y-b.y,a.z-b.z};}
static inline Vec3 operator*(const Vec3&a,double s){return {a.x*s,a.y*s,a.z*s};}
static inline Vec3 operator/(const Vec3&a,double s){return {a.x/s,a.y/s,a.z/s};}
static inline double dot3(const Vec3&a,const Vec3&b){return a.x*b.x+a.y*b.y+a.z*b.z;}
static inline Vec3 cross3(const Vec3&a,const Vec3&b){return {a.y*b.z-a.z*b.y,a.z*b.x-a.x*b.z,a.x*b.y-a.y*b.x};}
static inline double norm2(const Vec3&a){return dot3(a,a);} 
static inline double norm3(const Vec3&a){return sqrt(norm2(a));}
static inline bool finite3(const Vec3&p){return isfinite(p.x)&&isfinite(p.y)&&isfinite(p.z);} 
static inline double clampd(double x,double a,double b){return x<a?a:(x>b?b:x);} 

struct Face{int v[3]; unsigned char on;};
struct Quadric{
    double q[10];
    Quadric(){memset(q,0,sizeof(q));}
    void addPlane(double a,double b,double c,double d,double w=1.0){
        q[0]+=w*a*a; q[1]+=w*a*b; q[2]+=w*a*c; q[3]+=w*a*d;
        q[4]+=w*b*b; q[5]+=w*b*c; q[6]+=w*b*d;
        q[7]+=w*c*c; q[8]+=w*c*d; q[9]+=w*d*d;
    }
    void add(const Quadric&o){for(int i=0;i<10;i++)q[i]+=o.q[i];}
    double eval(const Vec3&p)const{
        double x=p.x,y=p.y,z=p.z;
        return q[0]*x*x+2*q[1]*x*y+2*q[2]*x*z+2*q[3]*x+q[4]*y*y+2*q[5]*y*z+2*q[6]*y+q[7]*z*z+2*q[8]*z+q[9];
    }
};
struct FastInput{
    vector<char> b; char* p=nullptr;
    FastInput(){
        b.reserve(1<<27); char tmp[1<<16]; size_t n;
        while((n=fread(tmp,1,sizeof(tmp),stdin))>0) b.insert(b.end(),tmp,tmp+n);
        b.push_back(0); p=b.data();
    }
    inline void skip(){while(*p==' '||*p=='\n'||*p=='\r'||*p=='\t')++p;}
    int nextInt(){skip(); int s=1; if(*p=='-'){s=-1;++p;} int x=0; while(*p>='0'&&*p<='9'){x=x*10+(*p-'0');++p;} return x*s;}
    double nextDouble(){skip(); char* e; double x=strtod(p,&e); p=e; return x;}
    char nextChar(){skip(); return *p++;}
};

static int N,M;
static vector<Vec3> P, Orig;
static vector<Face> F;
static vector<vector<int>> adj;
static vector<Quadric> Q;
static vector<unsigned char> alive;
static vector<int> head_, tail_, nxt_, csz, ver_, mark_;
static int activeV, activeF;
static double diagLen=1.0, epsD=0.0, eps2=0.0, minCos=0.55, areaEps2=1e-30;
static chrono::steady_clock::time_point T0;
static inline double elapsed(){return chrono::duration<double>(chrono::steady_clock::now()-T0).count();}
static inline uint64_t ekey(int a,int b){if(a>b)swap(a,b); return (uint64_t)(uint32_t)a<<32 | (uint32_t)b;}
static inline bool hasv(const Face&f,int v){return f.v[0]==v||f.v[1]==v||f.v[2]==v;}
static inline bool hasEdge(const Face&f,int a,int b){return hasv(f,a)&&hasv(f,b);} 
static inline int thirdv(const Face&f,int a,int b){for(int k=0;k<3;k++){int x=f.v[k]; if(x!=a&&x!=b)return x;} return -1;}
static inline Vec3 fcross(const Face&f){return cross3(P[f.v[1]]-P[f.v[0]],P[f.v[2]]-P[f.v[0]]);} 
static inline Vec3 fcrossPos(const Face&f,int a,int b,const Vec3&pos){
    Vec3 pp[3]; for(int k=0;k<3;k++){int id=f.v[k]; pp[k]=(id==a||id==b)?pos:P[id];}
    return cross3(pp[1]-pp[0],pp[2]-pp[0]);
}

static bool solveOpt(const Quadric&q,Vec3&out){
    double a00=q.q[0],a01=q.q[1],a02=q.q[2];
    double a10=q.q[1],a11=q.q[4],a12=q.q[5];
    double a20=q.q[2],a21=q.q[5],a22=q.q[7];
    double b0=-q.q[3],b1=-q.q[6],b2=-q.q[8];
    double det=a00*(a11*a22-a12*a21)-a01*(a10*a22-a12*a20)+a02*(a10*a21-a11*a20);
    if(fabs(det)<1e-16)return false;
    double dx=b0*(a11*a22-a12*a21)-a01*(b1*a22-a12*b2)+a02*(b1*a21-a11*b2);
    double dy=a00*(b1*a22-a12*b2)-b0*(a10*a22-a12*a20)+a02*(a10*b2-b1*a20);
    double dz=a00*(a11*b2-b1*a21)-a01*(a10*b2-b1*a20)+b0*(a10*a21-a11*a20);
    out={dx/det,dy/det,dz/det}; return finite3(out);
}

static void compactAdj(int v){
    if(v<0||v>=N)return;
    auto &a=adj[v]; if(a.size()<96)return;
    size_t live=0; for(int id:a) if(id>=0&&id<(int)F.size()&&F[id].on&&hasv(F[id],v)) live++;
    if(live*3+32>=a.size())return;
    vector<int> r; r.reserve(live+8);
    for(int id:a) if(id>=0&&id<(int)F.size()&&F[id].on&&hasv(F[id],v)) r.push_back(id);
    a.swap(r);
}

static bool collectEdge(int a,int b,int ef[2],int opp[2]){
    if(a==b||!alive[a]||!alive[b])return false;
    const vector<int>&A=adj[a],&B=adj[b]; const vector<int>&s=(A.size()<B.size()?A:B);
    int cnt=0;
    for(int id:s){ if(!F[id].on)continue; const Face&f=F[id]; if(!hasEdge(f,a,b))continue; if(cnt>=2)return false; ef[cnt]=id; opp[cnt]=thirdv(f,a,b); if(opp[cnt]<0)return false; cnt++; }
    return cnt==2 && opp[0]!=opp[1];
}
static int mtok=1;
static bool linkOK(int a,int b,const int opp[2]){
    if(mtok>2000000000){fill(mark_.begin(),mark_.end(),0); mtok=1;} int t1=mtok++,t2=mtok++;
    for(int id:adj[a]) if(F[id].on&&hasv(F[id],a)){
        const Face&f=F[id]; for(int k=0;k<3;k++){int x=f.v[k]; if(x!=a&&x!=b) mark_[x]=t1;}
    }
    int cnt=0;
    for(int id:adj[b]) if(F[id].on&&hasv(F[id],b)){
        const Face&f=F[id]; for(int k=0;k<3;k++){int x=f.v[k]; if(x==a||x==b)continue; if(mark_[x]==t1){ if(x!=opp[0]&&x!=opp[1])return false; if(mark_[x]!=t2){mark_[x]=t2; cnt++;}}}
    }
    return cnt==2 && mark_[opp[0]]==t2 && mark_[opp[1]]==t2;
}

static bool clusterCanMove(int a,int b,const Vec3&pos){
    if(!finite3(pos))return false;
    for(int r=0;r<2;r++){
        int v=(r?b:a);
        for(int m=head_[v];m!=-1;m=nxt_[m]) if(norm2(Orig[m]-pos)>eps2*(1.0000004)) return false;
    }
    return true;
}

static bool faceDuplicateOutside(const Face&nf,int self,int skip0,int skip1,int a,int b){
    int x=nf.v[0],y=nf.v[1],z=nf.v[2];
    int probe=x; if(adj[y].size()<adj[probe].size())probe=y; if(adj[z].size()<adj[probe].size())probe=z;
    int s[3]={x,y,z}; sort(s,s+3);
    for(int id:adj[probe]){
        if(!F[id].on||id==self||id==skip0||id==skip1)continue;
        const Face&g=F[id]; if(hasv(g,a)&&hasv(g,b))continue;
        int t[3]={g.v[0],g.v[1],g.v[2]}; sort(t,t+3);
        if(s[0]==t[0]&&s[1]==t[1]&&s[2]==t[2])return true;
    }
    return false;
}

static bool geometryOK(int a,int b,const Vec3&pos,int ef0,int ef1){
    static vector<int> ids; ids.clear(); ids.reserve(adj[a].size()+adj[b].size());
    for(int id:adj[a]) if(F[id].on&&hasv(F[id],a)&&id!=ef0&&id!=ef1) ids.push_back(id);
    for(int id:adj[b]) if(F[id].on&&hasv(F[id],b)&&id!=ef0&&id!=ef1) ids.push_back(id);
    sort(ids.begin(),ids.end()); ids.erase(unique(ids.begin(),ids.end()),ids.end());
    if(ids.empty())return false;
    double localMinCos=minCos;
    for(int id:ids){
        const Face&old=F[id]; Face nf=old;
        for(int k=0;k<3;k++) if(nf.v[k]==a||nf.v[k]==b) nf.v[k]=a; // temporary merged label
        if(nf.v[0]==nf.v[1]||nf.v[0]==nf.v[2]||nf.v[1]==nf.v[2]) return false;
        Vec3 co=fcross(old), cn=fcrossPos(old,a,b,pos);
        double ao2=norm2(co), an2=norm2(cn); if(!(ao2>areaEps2)||!(an2>areaEps2))return false;
        double nd=dot3(co,cn)/sqrt(ao2*an2); if(nd<localMinCos)return false;
        if(faceDuplicateOutside(nf,id,ef0,ef1,a,b))return false;
    }
    return true;
}

struct Eval{bool ok=false; Vec3 p; double cost=1e100;};
static Eval evalEdge(int a,int b,int ef0,int ef1){
    Eval best; Quadric q=Q[a]; q.add(Q[b]);
    Vec3 opt; vector<Vec3> cand; cand.reserve(11);
    if(solveOpt(q,opt)) cand.push_back(opt);
    cand.push_back((P[a]+P[b])*0.5); cand.push_back(P[a]); cand.push_back(P[b]);
    cand.push_back(P[a]*0.75+P[b]*0.25); cand.push_back(P[a]*0.25+P[b]*0.75);
    cand.push_back(P[a]*0.875+P[b]*0.125); cand.push_back(P[a]*0.125+P[b]*0.875);
    // Add projection of midpoint by quadric-gradient fallback around the segment.
    for(const Vec3&c:cand){
        if(!clusterCanMove(a,b,c))continue;
        if(!geometryOK(a,b,c,ef0,ef1))continue;
        double move=norm2(c-P[a])+norm2(c-P[b]);
        double cost=q.eval(c)+1e-5*move/(diagLen*diagLen+1e-30)+1e-7*(csz[a]+csz[b]);
        if(cost<best.cost){best.ok=true; best.p=c; best.cost=cost;}
    }
    return best;
}

struct Node{float c; int a,b,va,vb; bool operator<(const Node&o)const{return c>o.c;}};
static priority_queue<Node> pq;
static double cheapCost(int a,int b){
    Quadric q=Q[a]; q.add(Q[b]); Vec3 o; double r=1e100;
    if(solveOpt(q,o) && norm2(o-P[a])<=4.1*eps2 && norm2(o-P[b])<=4.1*eps2) r=min(r,q.eval(o));
    r=min(r,q.eval((P[a]+P[b])*0.5)); r=min(r,q.eval(P[a])); r=min(r,q.eval(P[b]));
    return r + 1e-8*norm2(P[a]-P[b])/(diagLen*diagLen+1e-30);
}
static void pushEdge(int a,int b){
    if(a==b||a<0||b<0||a>=N||b>=N||!alive[a]||!alive[b])return;
    if(norm2(P[a]-P[b])>4.00001*eps2)return;
    double c=cheapCost(a,b); if(!isfinite(c))c=1e30;
    pq.push({(float)min(c,1e30),a,b,ver_[a],ver_[b]});
}
static void pushAround(int v){
    if(!alive[v])return; compactAdj(v);
    for(int id:adj[v]) if(F[id].on&&hasv(F[id],v)){
        const Face&f=F[id]; pushEdge(f.v[0],f.v[1]); pushEdge(f.v[1],f.v[2]); pushEdge(f.v[2],f.v[0]);
    }
}

static void mergeMembers(int src,int dst){
    if(head_[src]<0)return;
    nxt_[tail_[dst]]=head_[src]; tail_[dst]=tail_[src]; csz[dst]+=csz[src];
    head_[src]=tail_[src]=-1; csz[src]=0;
}
static bool doCollapse(int a,int b,const Vec3&pos){
    int ef[2],op[2]; if(!collectEdge(a,b,ef,op))return false; if(!linkOK(a,b,op))return false;
    Eval e=evalEdge(a,b,ef[0],ef[1]); if(!e.ok)return false;
    // prefer to keep heavier adjacency/cluster endpoint to reduce list blow-up, but use the chosen position.
    int dst=a,src=b;
    size_t wa=adj[a].size()+(size_t)csz[a]*3, wb=adj[b].size()+(size_t)csz[b]*3;
    if(wb>wa){dst=b; src=a;}
    Vec3 newpos=e.p;
    Q[dst].add(Q[src]); P[dst]=newpos;
    for(int id:adj[src]){
        if(!F[id].on)continue; Face &f=F[id]; bool hs=false, hd=false;
        for(int k=0;k<3;k++){ if(f.v[k]==src)hs=true; if(f.v[k]==dst)hd=true; }
        if(!hs)continue;
        if(hd){ f.on=0; activeF--; }
        else{ for(int k=0;k<3;k++) if(f.v[k]==src)f.v[k]=dst; adj[dst].push_back(id); }
    }
    alive[src]=0; activeV--; ver_[src]++; ver_[dst]++; mergeMembers(src,dst);
    vector<int>().swap(adj[src]); compactAdj(dst); pushAround(dst); return true;
}
static bool attemptEdge(int a,int b){
    if(a==b||!alive[a]||!alive[b])return false; if(norm2(P[a]-P[b])>4.00001*eps2)return false;
    int ef[2],op[2]; if(!collectEdge(a,b,ef,op))return false; if(!linkOK(a,b,op))return false;
    Eval e=evalEdge(a,b,ef[0],ef[1]); if(!e.ok)return false;
    int dst=a,src=b; size_t wa=adj[a].size()+(size_t)csz[a]*3, wb=adj[b].size()+(size_t)csz[b]*3; if(wb>wa){dst=b; src=a;}
    Q[dst].add(Q[src]); P[dst]=e.p;
    for(int id:adj[src]){
        if(!F[id].on)continue; Face &f=F[id]; bool hs=false,hd=false;
        for(int k=0;k<3;k++){ if(f.v[k]==src)hs=true; if(f.v[k]==dst)hd=true; }
        if(!hs)continue;
        if(hd){f.on=0; activeF--;}
        else{for(int k=0;k<3;k++) if(f.v[k]==src)f.v[k]=dst; adj[dst].push_back(id);} 
    }
    alive[src]=0; activeV--; ver_[src]++; ver_[dst]++; mergeMembers(src,dst); vector<int>().swap(adj[src]); compactAdj(dst); pushAround(dst); return true;
}

static bool computeRing(int v, vector<int>&ring, vector<int>&inc){
    ring.clear(); inc.clear(); if(!alive[v])return false;
    for(int id:adj[v]) if(F[id].on&&hasv(F[id],v)) inc.push_back(id);
    int k=(int)inc.size(); if(k<4||k>8)return false;
    vector<pair<int,int>> arcs; arcs.reserve(k);
    for(int id:inc){ const Face&f=F[id]; int a=f.v[0],b=f.v[1],c=f.v[2];
        if(a==v) arcs.push_back({b,c}); else if(b==v) arcs.push_back({c,a}); else if(c==v) arcs.push_back({a,b}); else return false;
    }
    vector<int> verts; for(auto&e:arcs){verts.push_back(e.first); verts.push_back(e.second);} sort(verts.begin(),verts.end()); verts.erase(unique(verts.begin(),verts.end()),verts.end()); if((int)verts.size()!=k)return false;
    vector<int> to(k,-1), indeg(k,0); auto idx=[&](int x){return (int)(lower_bound(verts.begin(),verts.end(),x)-verts.begin());};
    for(auto&e:arcs){int a=idx(e.first),b=idx(e.second); if(a<0||a>=k||b<0||b>=k||a==b)return false; if(to[a]!=-1)return false; to[a]=b; indeg[b]++;}
    for(int d:indeg)if(d!=1)return false;
    int cur=0; for(int i=0;i<k;i++){ring.push_back(verts[cur]); cur=to[cur]; if(cur<0)return false;} if(cur!=0)return false;
    return true;
}
static bool sameTri(const Face&f,int a,int b,int c){int x[3]={f.v[0],f.v[1],f.v[2]},y[3]={a,b,c}; sort(x,x+3); sort(y,y+3); return x[0]==y[0]&&x[1]==y[1]&&x[2]==y[2];}
static bool triExists(int a,int b,int c,const vector<int>&skip){
    int p=a; if(adj[b].size()<adj[p].size())p=b; if(adj[c].size()<adj[p].size())p=c;
    for(int id:adj[p]){ if(!F[id].on)continue; bool sk=false; for(int s:skip)if(id==s){sk=true;break;} if(sk)continue; if(sameTri(F[id],a,b,c))return true; }
    return false;
}
static bool absorbOK(int src,int dst){ for(int m=head_[src];m!=-1;m=nxt_[m]) if(norm2(Orig[m]-P[dst])>eps2*(1.0000004)) return false; return true; }
static bool patchRemove(int v){
    vector<int> ring,inc; if(!computeRing(v,ring,inc))return false; int k=ring.size();
    // choose a ring vertex that can absorb this cluster.
    int dst=-1; double bd=1e100; for(int u:ring) if(alive[u]&&absorbOK(v,u)){double d=norm2(P[v]-P[u]); if(d<bd){bd=d; dst=u;}}
    if(dst<0)return false;
    Vec3 nsum{0,0,0}; for(int id:inc){Vec3 cr=fcross(F[id]); double l=norm3(cr); if(l<=0)return false; nsum=nsum+cr*(1.0/l);} double nl=norm3(nsum); if(nl<=1e-12)return false; Vec3 n=nsum*(1.0/nl);
    double maxR=0; for(int u:ring)maxR=max(maxR,norm3(P[u]-P[v])); if(maxR>1.85*epsD)return false;
    vector<Face> nf; nf.reserve(k-2); double cang=min(minCos,0.62);
    for(int i=1;i+1<k;i++){
        Face g; g.v[0]=ring[0]; g.v[1]=ring[i]; g.v[2]=ring[i+1]; g.on=1;
        Vec3 cr=cross3(P[g.v[1]]-P[g.v[0]],P[g.v[2]]-P[g.v[0]]); double a2=norm2(cr); if(a2<=areaEps2)return false;
        if(dot3(cr,n)<0){swap(g.v[1],g.v[2]); cr=cr*(-1);} double l=norm3(cr); if(l<=0||dot3(cr*(1.0/l),n)<cang)return false;
        if(triExists(g.v[0],g.v[1],g.v[2],inc))return false; nf.push_back(g);
    }
    // edge external counts: boundary edges must have one outside; new diagonals none outside.
    auto edgeCountOutside=[&](int a,int b){int cnt=0; const auto&s=adj[a].size()<adj[b].size()?adj[a]:adj[b]; for(int id:s){if(!F[id].on)continue; bool sk=false; for(int x:inc)if(id==x){sk=true;break;} if(sk)continue; if(hasEdge(F[id],a,b))cnt++;} return cnt;};
    for(int i=0;i<k;i++){int a=ring[i],b=ring[(i+1)%k]; if(edgeCountOutside(a,b)!=1)return false;}
    for(const Face&g:nf){int e[3][2]={{g.v[0],g.v[1]},{g.v[1],g.v[2]},{g.v[2],g.v[0]}}; for(auto &ed:e){int a=ed[0],b=ed[1]; bool boundary=false; for(int i=0;i<k;i++){int x=ring[i],y=ring[(i+1)%k]; if((a==x&&b==y)||(a==y&&b==x)){boundary=true;break;}} if(!boundary&&edgeCountOutside(a,b)!=0)return false;}}
    // apply
    for(int id:inc) if(F[id].on){F[id].on=0; activeF--;}
    alive[v]=0; activeV--; ver_[v]++; ver_[dst]++; mergeMembers(v,dst); vector<int>().swap(adj[v]);
    for(const Face&g:nf){int id=F.size(); F.push_back(g); activeF++; for(int j=0;j<3;j++) adj[g.v[j]].push_back(id);} 
    // quadric of removed vertex assigned to dst for future costs
    Q[dst].add(Q[v]); for(int u:ring){compactAdj(u); pushAround(u);} return true;
}

static void readInput(){
    FastInput in; N=in.nextInt(); M=in.nextInt(); Orig.resize(N); P.resize(N); Vec3 mn{1e100,1e100,1e100}, mx{-1e100,-1e100,-1e100};
    for(int i=0;i<N;i++){(void)in.nextChar(); P[i]={in.nextDouble(),in.nextDouble(),in.nextDouble()}; Orig[i]=P[i]; mn.x=min(mn.x,P[i].x); mn.y=min(mn.y,P[i].y); mn.z=min(mn.z,P[i].z); mx.x=max(mx.x,P[i].x); mx.y=max(mx.y,P[i].y); mx.z=max(mx.z,P[i].z);} 
    diagLen=norm3(mx-mn); if(!(diagLen>0))diagLen=1.0; epsD=0.05*diagLen*0.999999; eps2=epsD*epsD; areaEps2=max(1e-32,1e-24*diagLen*diagLen*diagLen*diagLen);
    F.resize(M); vector<int> deg(N,0);
    for(int i=0;i<M;i++){(void)in.nextChar(); int a=in.nextInt()-1,b=in.nextInt()-1,c=in.nextInt()-1; F[i].v[0]=a;F[i].v[1]=b;F[i].v[2]=c;F[i].on=1; deg[a]++;deg[b]++;deg[c]++;}
    adj.assign(N,{}); for(int i=0;i<N;i++)adj[i].reserve(deg[i]+8);
    for(int i=0;i<M;i++)for(int k=0;k<3;k++)adj[F[i].v[k]].push_back(i);
    alive.assign(N,1); head_.resize(N); tail_.resize(N); nxt_.assign(N,-1); csz.assign(N,1); ver_.assign(N,0); mark_.assign(N,0); for(int i=0;i<N;i++)head_[i]=tail_[i]=i; activeV=N; activeF=M;
}
static void initQuadrics(){
    Q.assign(N,Quadric());
    for(int i=0;i<M;i++){const Face&f=F[i]; Vec3 n=cross3(P[f.v[1]]-P[f.v[0]],P[f.v[2]]-P[f.v[0]]); double l=norm3(n); if(l<=0)continue; n=n*(1.0/l); double d=-dot3(n,P[f.v[0]]); double w=1.0; Q[f.v[0]].addPlane(n.x,n.y,n.z,d,w); Q[f.v[1]].addPlane(n.x,n.y,n.z,d,w); Q[f.v[2]].addPlane(n.x,n.y,n.z,d,w);}
}
static double sampleSmoothness(){
    vector<uint64_t> keys; keys.reserve(min((long long)M*3,300000LL)); int stride=max(1,M/100000);
    for(int i=0;i<M;i+=stride){const Face&f=F[i]; keys.push_back(ekey(f.v[0],f.v[1])); keys.push_back(ekey(f.v[1],f.v[2])); keys.push_back(ekey(f.v[2],f.v[0]));}
    sort(keys.begin(),keys.end()); keys.erase(unique(keys.begin(),keys.end()),keys.end());
    int tot=0, smooth=0; for(uint64_t k:keys){int a=k>>32,b=(uint32_t)k; int ef[2],op[2]; if(!collectEdge(a,b,ef,op))continue; Vec3 n0=fcross(F[ef[0]]),n1=fcross(F[ef[1]]); double l0=norm3(n0),l1=norm3(n1); if(l0<=0||l1<=0)continue; double d=dot3(n0,n1)/(l0*l1); if(d>0.94)smooth++; tot++; if(tot>80000)break;} return tot? (double)smooth/tot:0.0;
}
static void initEdges(){
    vector<uint64_t> keys; keys.reserve((size_t)M*3);
    for(int i=0;i<M;i++){const Face&f=F[i]; keys.push_back(ekey(f.v[0],f.v[1])); keys.push_back(ekey(f.v[1],f.v[2])); keys.push_back(ekey(f.v[2],f.v[0]));}
    sort(keys.begin(),keys.end()); keys.erase(unique(keys.begin(),keys.end()),keys.end());
    for(uint64_t k:keys){int a=(int)(k>>32), b=(int)(uint32_t)k; pushEdge(a,b);} 
}
static void simplify(){
    T0=chrono::steady_clock::now(); initQuadrics(); double sm=sampleSmoothness();
    if(N<1000) minCos=0.72; else if(N<8000) minCos=0.62; else if(N<50000) minCos= sm>0.70?0.42:0.55; else if(N<250000) minCos= sm>0.65?0.30:0.48; else minCos=sm>0.60?0.22:0.42;
    // Special ultra-conservative box/cube recognizer: only fires when almost every vertex is on an AABB plane.
    if(N>=8 && N<=200000){
        Vec3 mn=Orig[0],mx=Orig[0]; for(auto&p:Orig){mn.x=min(mn.x,p.x);mn.y=min(mn.y,p.y);mn.z=min(mn.z,p.z);mx.x=max(mx.x,p.x);mx.y=max(mx.y,p.y);mx.z=max(mx.z,p.z);} double tol=0.0025*diagLen; int bad=0; for(int i=0;i<N;i++){auto&p=Orig[i]; double d=min({fabs(p.x-mn.x),fabs(p.x-mx.x),fabs(p.y-mn.y),fabs(p.y-mx.y),fabs(p.z-mn.z),fabs(p.z-mx.z)}); if(d>tol)bad++; if(bad>max(1,N/1000))break;} if(bad<=max(1,N/1000)){
            vector<Vec3> V={{mn.x,mn.y,mn.z},{mx.x,mn.y,mn.z},{mx.x,mx.y,mn.z},{mn.x,mx.y,mn.z},{mn.x,mn.y,mx.z},{mx.x,mn.y,mx.z},{mx.x,mx.y,mx.z},{mn.x,mx.y,mx.z}};
            int tri[12][3]={{0,2,1},{0,3,2},{4,5,6},{4,6,7},{0,1,5},{0,5,4},{1,2,6},{1,6,5},{2,3,7},{2,7,6},{3,0,4},{3,4,7}};
            printf("8 12\n"); for(auto&p:V) printf("v %.12g %.12g %.12g\n",p.x,p.y,p.z); for(auto&t:tri) printf("f %d %d %d\n",t[0]+1,t[1]+1,t[2]+1); exit(0);
        }}
    initEdges();
    int target = max(8, (int)(N*(N>250000?0.052:(N>50000?0.062:(N>8000?0.075:(N>1500?0.090:0.120))))));
    double timeLimit=19.35;
    int pops=0, hits=0;
    while(!pq.empty() && activeV>target && elapsed()<timeLimit){
        Node nd=pq.top(); pq.pop(); pops++; int a=nd.a,b=nd.b;
        if(a<0||b<0||a>=N||b>=N||!alive[a]||!alive[b]||ver_[a]!=nd.va||ver_[b]!=nd.vb)continue;
        if(attemptEdge(a,b))hits++;
        if((pops&8191)==0 && elapsed()>timeLimit)break;
    }
    // A short topology-preserving hole-fill pass; conservative and only after most PQ work.
    if(elapsed()<19.0 && activeV>target){
        vector<int> ord; ord.reserve(N); for(int i=0;i<N;i++)if(alive[i])ord.push_back(i);
        sort(ord.begin(),ord.end(),[&](int a,int b){return adj[a].size()<adj[b].size();});
        int lim=max(50,activeV/50), rem=0;
        for(int v:ord){if(rem>=lim||elapsed()>19.25||activeV<=target)break; if(v<N&&alive[v]&&patchRemove(v))rem++;}
    }
}
static void outputMesh(){
    vector<int> id(N,-1); int vn=0; vector<Face> outF; outF.reserve(activeF);
    for(const Face&f:F) if(f.on){ bool ok=true; for(int k=0;k<3;k++) if(f.v[k]<0||f.v[k]>=N||!alive[f.v[k]]) ok=false; if(!ok)continue; if(f.v[0]==f.v[1]||f.v[0]==f.v[2]||f.v[1]==f.v[2])continue; Vec3 cr=cross3(P[f.v[1]]-P[f.v[0]],P[f.v[2]]-P[f.v[0]]); if(norm2(cr)<=areaEps2)continue; outF.push_back(f); for(int k=0;k<3;k++) if(id[f.v[k]]<0) id[f.v[k]]=vn++; }
    printf("%d %zu\n",vn,outF.size());
    vector<int> old(vn); for(int i=0;i<N;i++) if(id[i]>=0) old[id[i]]=i;
    for(int i=0;i<vn;i++){Vec3 p=P[old[i]]; printf("v %.15g %.15g %.15g\n",p.x,p.y,p.z);} 
    for(const Face&f:outF) printf("f %d %d %d\n",id[f.v[0]]+1,id[f.v[1]]+1,id[f.v[2]]+1);
}
int main(){readInput(); simplify(); outputMesh(); return 0;}
