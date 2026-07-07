#include <bits/stdc++.h>
using namespace std;

struct P{double x,y,z;};
static inline P operator+(P a,P b){return {a.x+b.x,a.y+b.y,a.z+b.z};}
static inline P operator-(P a,P b){return {a.x-b.x,a.y-b.y,a.z-b.z};}
static inline P operator*(P a,double s){return {a.x*s,a.y*s,a.z*s};}
static inline double dotp(P a,P b){return a.x*b.x+a.y*b.y+a.z*b.z;}
static inline P crossp(P a,P b){return {a.y*b.z-a.z*b.y,a.z*b.x-a.x*b.z,a.x*b.y-a.y*b.x};}
static inline double n2(P a){return dotp(a,a);} 
static inline double norm(P a){return sqrt(n2(a));}
static inline P unit(P a){double l=norm(a); return l>0?a*(1.0/l):P{0,0,0};}

struct F{int a,b,c; bool alive; P n;};
struct Cand{double cost; int v,ver; bool operator<(Cand const& o)const{return cost>o.cost;}};

int N0,F0,aliveV,aliveF; double diagLen,tol,areaFloor;
vector<P> V0; vector<F> faces; vector<vector<int>> inc, cluster;
vector<char> alive; vector<int> ver;
unordered_map<unsigned long long,int> ecnt;
unordered_set<unsigned long long> fset;

static inline unsigned long long ekey(int a,int b){ if(a>b) swap(a,b); return (unsigned long long)(unsigned int)a<<32 | (unsigned int)b; }
static inline unsigned long long fkey3(int a,int b,int c){
    if(a>b) swap(a,b); if(b>c) swap(b,c); if(a>b) swap(a,b);
    // 21 bits each, enough for Kattis-sized meshes; fallback hash mix for larger is still unique enough for duplicate detection.
    return ((unsigned long long)(unsigned int)a*11995408973635179863ull) ^ ((unsigned long long)(unsigned int)b*10150724397891781847ull) ^ ((unsigned long long)(unsigned int)c*5714713534750834931ull);
}
static inline double area2tri(P a,P b,P c){return n2(crossp(b-a,c-a));}
static inline void addEdge(int a,int b,int d){ auto k=ekey(a,b); int &x=ecnt[k]; x+=d; if(x==0) ecnt.erase(k); }
static inline void addFaceEdges(const F& f,int d){ addEdge(f.a,f.b,d); addEdge(f.b,f.c,d); addEdge(f.c,f.a,d); }
static inline bool hasv(const F& f,int v){return f.a==v||f.b==v||f.c==v;}
static inline P fnormal(int a,int b,int c){ return unit(crossp(V0[b]-V0[a],V0[c]-V0[a])); }

static vector<char> slurp(){ vector<char>b; char buf[1<<16]; size_t n; while((n=fread(buf,1,sizeof(buf),stdin))) b.insert(b.end(),buf,buf+n); b.push_back(0); return b; }
static inline void sws(char*&p){while(*p==' '||*p=='\n'||*p=='\r'||*p=='\t')++p;}

static void readInput(){
    auto buf=slurp(); char* p=buf.data();
    N0=strtol(p,&p,10); F0=strtol(p,&p,10);
    V0.resize(N0); faces.reserve(F0+N0); inc.assign(N0,{}); cluster.assign(N0,{}); alive.assign(N0,1); ver.assign(N0,0);
    double mnx=1e100,mny=1e100,mnz=1e100,mxx=-1e100,mxy=-1e100,mxz=-1e100;
    for(int i=0;i<N0;i++){
        sws(p); if(*p=='v') ++p; double x=strtod(p,&p), y=strtod(p,&p), z=strtod(p,&p);
        V0[i]={x,y,z}; cluster[i].push_back(i);
        mnx=min(mnx,x); mny=min(mny,y); mnz=min(mnz,z); mxx=max(mxx,x); mxy=max(mxy,y); mxz=max(mxz,z);
    }
    for(int i=0;i<F0;i++){
        sws(p); if(*p=='f') ++p; int a=strtol(p,&p,10)-1,b=strtol(p,&p,10)-1,c=strtol(p,&p,10)-1;
        F f{a,b,c,true,unit(crossp(V0[b]-V0[a],V0[c]-V0[a]))}; faces.push_back(f);
        inc[a].push_back(i); inc[b].push_back(i); inc[c].push_back(i); addFaceEdges(f,1); fset.insert(fkey3(a,b,c));
    }
    aliveV=N0; aliveF=F0; diagLen=sqrt((mxx-mnx)*(mxx-mnx)+(mxy-mny)*(mxy-mny)+(mxz-mnz)*(mxz-mnz));
    tol=0.05*diagLen; areaFloor=max(1e-30,diagLen*diagLen*1e-28);
}

static void cleanInc(int v){
    vector<int> r; r.reserve(inc[v].size());
    for(int id:inc[v]) if(id>=0 && id<(int)faces.size() && faces[id].alive && hasv(faces[id],v)) r.push_back(id);
    inc[v].swap(r);
}

static double pointSegDist2(P p,P a,P b){ P ab=b-a; double d=n2(ab); double t=d>0?dotp(p-a,ab)/d:0; t=max(0.0,min(1.0,t)); P q=a+ab*t; return n2(p-q); }
static double pointTriDist2(P p,P a,P b,P c){
    // Christer Ericson closest-point test.
    P ab=b-a, ac=c-a, ap=p-a; double d1=dotp(ab,ap), d2=dotp(ac,ap);
    if(d1<=0&&d2<=0) return n2(ap);
    P bp=p-b; double d3=dotp(ab,bp), d4=dotp(ac,bp);
    if(d3>=0&&d4<=d3) return n2(bp);
    double vc=d1*d4-d3*d2;
    if(vc<=0&&d1>=0&&d3<=0){ double v=d1/(d1-d3); P q=a+ab*v; return n2(p-q); }
    P cp=p-c; double d5=dotp(ab,cp), d6=dotp(ac,cp);
    if(d6>=0&&d5<=d6) return n2(cp);
    double vb=d5*d2-d1*d6;
    if(vb<=0&&d2>=0&&d6<=0){ double w=d2/(d2-d6); P q=a+ac*w; return n2(p-q); }
    double va=d3*d6-d5*d4;
    if(va<=0&&(d4-d3)>=0&&(d5-d6)>=0){ double w=(d4-d3)/((d4-d3)+(d5-d6)); P q=b+(c-b)*w; return n2(p-q); }
    P n=crossp(ab,ac); double nn=n2(n); if(nn<=0) return min({n2(ap),n2(bp),n2(cp)});
    double dist=dotp(p-a,n); return dist*dist/nn;
}

struct Plan{ vector<int> ring; vector<array<int,3>> tris; double cost; };

static bool makeCycle(int v, vector<int>& ring, vector<int>& ifaces, P& avgN, double& oldArea){
    if(!alive[v]) return false; cleanInc(v); ifaces=inc[v]; int d=ifaces.size(); if(d<3||d>14) return false;
    vector<int> nb; nb.reserve(d); vector<array<int,2>> ledges; ledges.reserve(d); avgN={0,0,0}; oldArea=0;
    for(int id:ifaces){ F &f=faces[id]; int a=f.a,b=f.b,c=f.c; int x,y;
        if(a==v){x=b;y=c;} else if(b==v){x=c;y=a;} else if(c==v){x=a;y=b;} else return false;
        if(!alive[x]||!alive[y]||x==y) return false;
        ledges.push_back({x,y}); nb.push_back(x); nb.push_back(y);
        P cr=crossp(V0[f.b]-V0[f.a],V0[f.c]-V0[f.a]); double ar=norm(cr); oldArea+=0.5*ar; if(ar>0) avgN=avgN+cr;
    }
    sort(nb.begin(),nb.end()); nb.erase(unique(nb.begin(),nb.end()),nb.end()); if((int)nb.size()!=d) return false;
    vector<vector<int>> adj(nb.size());
    auto idx=[&](int x){return int(lower_bound(nb.begin(),nb.end(),x)-nb.begin());};
    for(auto e:ledges){ int i=idx(e[0]),j=idx(e[1]); if(i<0||i>=d||j<0||j>=d||nb[i]!=e[0]||nb[j]!=e[1]) return false; adj[i].push_back(j); adj[j].push_back(i); }
    for(auto &a:adj){ sort(a.begin(),a.end()); a.erase(unique(a.begin(),a.end()),a.end()); if(a.size()!=2) return false; }
    ring.clear(); ring.reserve(d); int start=0, prev=-1, cur=start;
    for(int k=0;k<d;k++){ ring.push_back(nb[cur]); int nx=(adj[cur][0]==prev?adj[cur][1]:adj[cur][0]); prev=cur; cur=nx; }
    if(cur!=start) return false;
    P poly={0,0,0}; P cen={0,0,0}; for(int x:ring) cen=cen+V0[x]; cen=cen*(1.0/d);
    for(int i=0;i<d;i++) poly=poly+crossp(V0[ring[i]]-cen,V0[ring[(i+1)%d]]-cen);
    if(dotp(poly,avgN)<0) reverse(ring.begin(),ring.end());
    avgN=unit(avgN); return norm(avgN)>0;
}

static bool evaluatePlan(int v, const vector<int>& ring, const vector<int>& ifaces, P avgN, double oldArea, Plan& best){
    int d=ring.size(); unordered_map<unsigned long long,int> rem,newc; rem.reserve(ifaces.size()*3+8);
    for(int id:ifaces){ F &f=faces[id]; rem[ekey(f.a,f.b)]++; rem[ekey(f.b,f.c)]++; rem[ekey(f.c,f.a)]++; }
    double coverLim=tol*0.985, surfLim=tol*(N0>50000?0.42:0.34); double cover2=coverLim*coverLim, surf2=surfLim*surfLim;
    // Every represented original point currently owned by v must be close to some surviving link vertex.
    double maxCover=0;
    for(int orig:cluster[v]){
        double md=1e300; for(int r:ring) md=min(md,n2(V0[orig]-V0[r]));
        if(md>cover2) return false; maxCover=max(maxCover,sqrt(md));
    }
    best.cost=1e300; bool ok=false; double minDot=(N0>40000?0.55:(N0>10000?0.66:0.78));
    for(int s=0;s<d;s++){
        vector<int> seq; seq.reserve(d); for(int k=0;k<d;k++) seq.push_back(ring[(s+k)%d]);
        vector<array<int,3>> tr; tr.reserve(d-2); bool bad=false; double newArea=0, ndErr=0;
        newc.clear();
        for(int i=1;i<d-1;i++){
            int a=seq[0],b=seq[i],c=seq[i+1]; if(a==b||b==c||a==c){bad=true;break;}
            if(fset.find(fkey3(a,b,c))!=fset.end()){ bad=true; break; }
            P cr=crossp(V0[b]-V0[a],V0[c]-V0[a]); double ar=norm(cr); if(ar*ar<=areaFloor){bad=true;break;}
            P nn=cr*(1.0/ar); double dp=dotp(nn,avgN);
            if(dp<0){ swap(b,c); cr=crossp(V0[b]-V0[a],V0[c]-V0[a]); ar=norm(cr); if(ar<=0){bad=true;break;} nn=cr*(1.0/ar); dp=dotp(nn,avgN); }
            if(dp<minDot){bad=true;break;}
            ndErr += (1.0-dp)*ar; newArea += 0.5*ar; tr.push_back({a,b,c});
            newc[ekey(a,b)]++; newc[ekey(b,c)]++; newc[ekey(c,a)]++;
        }
        if(bad) continue;
        // local edge counts must finish as closed 2-manifold counts.
        for(auto &kv:newc){ int base=0; auto it=ecnt.find(kv.first); if(it!=ecnt.end()) base=it->second; int r=0; auto ir=rem.find(kv.first); if(ir!=rem.end()) r=ir->second; if(base-r+kv.second!=2){bad=true;break;} }
        if(bad) continue;
        // removed radial edges must vanish, boundary edges must be restored by one new face.
        for(auto &kv:rem){ int add=0; auto ia=newc.find(kv.first); if(ia!=newc.end()) add=ia->second; int base=ecnt.count(kv.first)?ecnt[kv.first]:0; int fin=base-kv.second+add; if(fin!=0 && fin!=2){bad=true;break;} }
        if(bad) continue;
        double maxSurf=0;
        for(int orig:cluster[v]){
            double md=1e300; P p=V0[orig]; for(auto t:tr) md=min(md,pointTriDist2(p,V0[t[0]],V0[t[1]],V0[t[2]]));
            if(md>surf2){bad=true;break;} maxSurf=max(maxSurf,sqrt(md));
        }
        if(bad) continue;
        double areaRatio = oldArea>1e-300 ? fabs(newArea-oldArea)/oldArea : 0;
        if(areaRatio>(N0>40000?0.85:0.55)) continue;
        double cost = maxCover/(tol+1e-300)*0.9 + maxSurf/(surfLim+1e-300)*1.4 + ndErr/(oldArea+1e-300)*0.8 + 0.018*d + areaRatio*0.25;
        if(cost<best.cost){best.ring=ring; best.tris=tr; best.cost=cost; ok=true;}
    }
    return ok;
}

static bool getPlan(int v, Plan& p){ vector<int> ring,ifs; P n; double ar; return makeCycle(v,ring,ifs,n,ar)&&evaluatePlan(v,ring,ifs,n,ar,p); }

static void pushCand(int v, priority_queue<Cand>&pq){ if(!alive[v]) return; Plan p; if(getPlan(v,p)) pq.push({p.cost,v,ver[v]}); }

static bool applyRemove(int v, priority_queue<Cand>&pq){
    vector<int> ring,ifs; P avn; double oldArea; Plan pl; if(!makeCycle(v,ring,ifs,avn,oldArea)||!evaluatePlan(v,ring,ifs,avn,oldArea,pl)) return false;
    for(int id:ifs){ if(!faces[id].alive) continue; F f=faces[id]; faces[id].alive=false; aliveF--; addFaceEdges(f,-1); fset.erase(fkey3(f.a,f.b,f.c)); }
    for(auto t:pl.tris){ F nf{t[0],t[1],t[2],true,fnormal(t[0],t[1],t[2])}; int id=faces.size(); faces.push_back(nf); aliveF++; addFaceEdges(nf,1); fset.insert(fkey3(nf.a,nf.b,nf.c)); inc[nf.a].push_back(id); inc[nf.b].push_back(id); inc[nf.c].push_back(id); }
    // Reassign represented original vertices to their nearest surviving ring vertex.
    for(int orig:cluster[v]){ int br=ring[0]; double bd=n2(V0[orig]-V0[br]); for(int r:ring){ double dd=n2(V0[orig]-V0[r]); if(dd<bd){bd=dd;br=r;} } cluster[br].push_back(orig); }
    vector<int>().swap(cluster[v]); alive[v]=0; aliveV--; ver[v]++;
    vector<int>().swap(inc[v]);
    for(int r:ring){ ver[r]++; cleanInc(r); }
    // refresh candidates in two rings.
    unordered_set<int> todo; todo.reserve(64); for(int r:ring){ todo.insert(r); for(int id:inc[r]){ F &f=faces[id]; if(f.alive){ todo.insert(f.a); todo.insert(f.b); todo.insert(f.c); } } }
    for(int x:todo) if(alive[x]) pushCand(x,pq);
    return true;
}

static bool finalValid(){
    vector<int> id(N0,-1); int nv=0,nf=0; for(int i=0;i<N0;i++) if(alive[i]) id[i]=nv++;
    if(nv<=0) return false; unordered_map<unsigned long long,int> ee; ee.reserve(aliveF*4+10);
    for(auto &f:faces) if(f.alive){ int a=id[f.a],b=id[f.b],c=id[f.c]; if(a<0||b<0||c<0||a==b||b==c||a==c) return false; if(area2tri(V0[f.a],V0[f.b],V0[f.c])<=areaFloor) return false; ee[ekey(a,b)]++; ee[ekey(b,c)]++; ee[ekey(c,a)]++; nf++; }
    if(nf<=0) return false; for(auto &kv:ee) if(kv.second!=2) return false;
    double lim2=tol*tol*1.0000001;
    for(int i=0;i<N0;i++) if(!alive[i]){ double md=1e300; for(int j=0;j<N0;j++) if(alive[j]){ double d=n2(V0[i]-V0[j]); if(d<md) md=d; } if(md>lim2) return false; }
    return true;
}

static void outputIdentity(){
    printf("%d %d\n",N0,F0); for(auto p:V0) printf("v %.17g %.17g %.17g\n",p.x,p.y,p.z);
    for(int i=0;i<F0;i++){F &f=faces[i]; printf("f %d %d %d\n",f.a+1,f.b+1,f.c+1);} }

static void outputMesh(){
    vector<int> id(N0,-1); vector<int> rev; rev.reserve(aliveV); for(int i=0;i<N0;i++) if(alive[i]){id[i]=rev.size(); rev.push_back(i);} 
    vector<array<int,3>> out; out.reserve(aliveF); for(auto &f:faces) if(f.alive) out.push_back({id[f.a],id[f.b],id[f.c]});
    printf("%d %d\n",(int)rev.size(),(int)out.size()); for(int oi:rev){P p=V0[oi]; printf("v %.17g %.17g %.17g\n",p.x,p.y,p.z);} for(auto t:out) printf("f %d %d %d\n",t[0]+1,t[1]+1,t[2]+1);
}

int main(){
    ios::sync_with_stdio(false); cin.tie(nullptr);
    readInput();
    if(N0<=8 || F0<=12 || diagLen<=0){ outputIdentity(); return 0; }
    priority_queue<Cand> pq; for(int i=0;i<N0;i++) pushCand(i,pq);
    clock_t st=clock(); int attempts=0, done=0; double tlim=8.5;
    while(!pq.empty()){
        if((++attempts&255)==0 && double(clock()-st)/CLOCKS_PER_SEC>tlim) break;
        Cand c=pq.top(); pq.pop(); if(!alive[c.v]||ver[c.v]!=c.ver) continue;
        if(applyRemove(c.v,pq)) done++;
        if(aliveV<12) break;
    }
    if(!finalValid()) outputIdentity(); else outputMesh();
    return 0;
}