#include <bits/stdc++.h>
using namespace std;

struct Vec3{ double x,y,z; };
static inline Vec3 operator+(Vec3 a,Vec3 b){return {a.x+b.x,a.y+b.y,a.z+b.z};}
static inline Vec3 operator-(Vec3 a,Vec3 b){return {a.x-b.x,a.y-b.y,a.z-b.z};}
static inline Vec3 operator*(Vec3 a,double s){return {a.x*s,a.y*s,a.z*s};}
static inline double dotv(Vec3 a,Vec3 b){return a.x*b.x+a.y*b.y+a.z*b.z;}
static inline Vec3 crossv(Vec3 a,Vec3 b){return {a.y*b.z-a.z*b.y,a.z*b.x-a.x*b.z,a.x*b.y-a.y*b.x};}
static inline double len2(Vec3 a){return dotv(a,a);} 
static inline double lenv(Vec3 a){return sqrt(len2(a));}
static inline Vec3 unitv(Vec3 a){ double l=lenv(a); return l>0?a*(1.0/l):Vec3{0,0,0}; }

struct Face{ int a,b,c; unsigned char alive; Vec3 n; };
struct Plan{ vector<array<int,3>> tris; double score=1e300; int deg=0; };
struct Cand{ double score; int v,version; bool operator<(Cand const& o)const{return score>o.score;} };

static int N0,F0,aliveV,aliveF;
static vector<Vec3> P0;
static vector<Face> faces;
static vector<vector<int>> inc, owned;
static vector<unsigned char> valive;
static vector<int> vver;
static unordered_map<unsigned long long,int> edgeCnt;
static unordered_set<unsigned long long> faceSet;
static double bboxDiag, hausTol, areaFloor;
static clock_t startClock;

static inline unsigned long long edgeKey(int a,int b){ if(a>b) swap(a,b); return (unsigned long long)(unsigned int)a<<32 | (unsigned int)b; }
static inline unsigned long long faceKey(int a,int b,int c){
    if(a>b) swap(a,b); if(b>c) swap(b,c); if(a>b) swap(a,b);
    unsigned long long x=(unsigned int)a, y=(unsigned int)b, z=(unsigned int)c;
    return (x*11995408973635179863ull) ^ (y*10150724397891781847ull) ^ (z*5714713534750834931ull);
}
static inline double triArea2(int a,int b,int c){ return len2(crossv(P0[b]-P0[a],P0[c]-P0[a])); }
static inline Vec3 triNormal(int a,int b,int c){ return unitv(crossv(P0[b]-P0[a],P0[c]-P0[a])); }
static inline bool hasVertex(const Face& f,int v){ return f.a==v || f.b==v || f.c==v; }
static inline double elapsed(){ return double(clock()-startClock)/CLOCKS_PER_SEC; }

static void edgeAdd(int a,int b,int delta){ unsigned long long k=edgeKey(a,b); int &x=edgeCnt[k]; x+=delta; if(x==0) edgeCnt.erase(k); }
static void faceEdgesAdd(const Face& f,int delta){ edgeAdd(f.a,f.b,delta); edgeAdd(f.b,f.c,delta); edgeAdd(f.c,f.a,delta); }

static vector<char> readAll(){ vector<char> r; char buf[1<<16]; size_t n; while((n=fread(buf,1,sizeof(buf),stdin))>0) r.insert(r.end(),buf,buf+n); r.push_back(0); return r; }
static inline void skipWs(char*& p){ while(*p==' '||*p=='\n'||*p=='\r'||*p=='\t') ++p; }

static bool readMesh(){
    vector<char> buf=readAll(); char* p=buf.data();
    N0=(int)strtol(p,&p,10); F0=(int)strtol(p,&p,10);
    if(N0<=0||F0<=0) return false;
    P0.resize(N0); inc.assign(N0,{}); owned.assign(N0,{}); valive.assign(N0,1); vver.assign(N0,0);
    faces.reserve((size_t)F0 + (size_t)N0*2);
    double mnx=1e100,mny=1e100,mnz=1e100,mxx=-1e100,mxy=-1e100,mxz=-1e100;
    for(int i=0;i<N0;i++){
        skipWs(p); if(*p=='v') ++p;
        double x=strtod(p,&p), y=strtod(p,&p), z=strtod(p,&p);
        P0[i]={x,y,z}; owned[i].push_back(i);
        mnx=min(mnx,x); mny=min(mny,y); mnz=min(mnz,z); mxx=max(mxx,x); mxy=max(mxy,y); mxz=max(mxz,z);
    }
    edgeCnt.reserve((size_t)F0*4+100); faceSet.reserve((size_t)F0*2+100);
    for(int i=0;i<F0;i++){
        skipWs(p); if(*p=='f') ++p;
        int a=(int)strtol(p,&p,10)-1, b=(int)strtol(p,&p,10)-1, c=(int)strtol(p,&p,10)-1;
        if(a<0||b<0||c<0||a>=N0||b>=N0||c>=N0||a==b||b==c||a==c) return false;
        Face f{a,b,c,1,triNormal(a,b,c)};
        int id=(int)faces.size(); faces.push_back(f);
        inc[a].push_back(id); inc[b].push_back(id); inc[c].push_back(id);
        faceEdgesAdd(f,1); faceSet.insert(faceKey(a,b,c));
    }
    aliveV=N0; aliveF=F0;
    bboxDiag=sqrt((mxx-mnx)*(mxx-mnx)+(mxy-mny)*(mxy-mny)+(mxz-mnz)*(mxz-mnz));
    hausTol=0.05*bboxDiag; areaFloor=max(1e-30,bboxDiag*bboxDiag*1e-28);
    return true;
}

static void cleanInc(int v){
    if(!valive[v]){ vector<int>().swap(inc[v]); return; }
    vector<int> r; r.reserve(inc[v].size());
    for(int id:inc[v]) if(id>=0 && id<(int)faces.size() && faces[id].alive && hasVertex(faces[id],v)) r.push_back(id);
    sort(r.begin(),r.end()); r.erase(unique(r.begin(),r.end()),r.end());
    inc[v].swap(r);
}

static double pointTriDist2(Vec3 p,Vec3 a,Vec3 b,Vec3 c){
    Vec3 ab=b-a, ac=c-a, ap=p-a;
    double d1=dotv(ab,ap), d2=dotv(ac,ap);
    if(d1<=0.0 && d2<=0.0) return len2(ap);
    Vec3 bp=p-b; double d3=dotv(ab,bp), d4=dotv(ac,bp);
    if(d3>=0.0 && d4<=d3) return len2(bp);
    double vc=d1*d4-d3*d2;
    if(vc<=0.0 && d1>=0.0 && d3<=0.0){ double v=d1/(d1-d3); Vec3 q=a+ab*v; return len2(p-q); }
    Vec3 cp=p-c; double d5=dotv(ab,cp), d6=dotv(ac,cp);
    if(d6>=0.0 && d5<=d6) return len2(cp);
    double vb=d5*d2-d1*d6;
    if(vb<=0.0 && d2>=0.0 && d6<=0.0){ double w=d2/(d2-d6); Vec3 q=a+ac*w; return len2(p-q); }
    double va=d3*d6-d5*d4;
    if(va<=0.0 && (d4-d3)>=0.0 && (d5-d6)>=0.0){ double w=(d4-d3)/((d4-d3)+(d5-d6)); Vec3 q=b+(c-b)*w; return len2(p-q); }
    Vec3 n=crossv(ab,ac); double nn=len2(n); if(nn<=0) return min({len2(ap),len2(bp),len2(cp)});
    double dist=dotv(p-a,n); return dist*dist/nn;
}

static bool linkCycle(int v, vector<int>& ring, vector<int>& starFaces, Vec3& avgN, double& oldArea){
    if(!valive[v]) return false;
    cleanInc(v); starFaces=inc[v]; int d=(int)starFaces.size();
    if(d<3 || d>22) return false;
    vector<int> nodes; nodes.reserve(d*2);
    vector<pair<int,int>> linkEdges; linkEdges.reserve(d);
    avgN={0,0,0}; oldArea=0.0;
    for(int id:starFaces){
        Face &f=faces[id]; if(!f.alive || !hasVertex(f,v)) return false;
        int x=-1,y=-1;
        if(f.a==v){x=f.b;y=f.c;} else if(f.b==v){x=f.c;y=f.a;} else {x=f.a;y=f.b;}
        if(x<0||y<0||x==y||!valive[x]||!valive[y]) return false;
        nodes.push_back(x); nodes.push_back(y); linkEdges.push_back({x,y});
        Vec3 cr=crossv(P0[f.b]-P0[f.a],P0[f.c]-P0[f.a]); double ar=lenv(cr); oldArea += 0.5*ar; if(ar>0) avgN=avgN+cr;
    }
    sort(nodes.begin(),nodes.end()); nodes.erase(unique(nodes.begin(),nodes.end()),nodes.end());
    if((int)nodes.size()!=d) return false;
    vector<vector<int>> adj(d);
    auto idxOf=[&](int x)->int{ auto it=lower_bound(nodes.begin(),nodes.end(),x); if(it==nodes.end()||*it!=x) return -1; return int(it-nodes.begin()); };
    for(auto e:linkEdges){ int a=idxOf(e.first), b=idxOf(e.second); if(a<0||b<0||a==b) return false; adj[a].push_back(b); adj[b].push_back(a); }
    for(auto &a:adj){ sort(a.begin(),a.end()); a.erase(unique(a.begin(),a.end()),a.end()); if(a.size()!=2) return false; }
    ring.clear(); ring.reserve(d); int cur=0, prev=-1;
    for(int k=0;k<d;k++){ ring.push_back(nodes[cur]); int nx=(adj[cur][0]==prev?adj[cur][1]:adj[cur][0]); prev=cur; cur=nx; }
    if(cur!=0) return false;
    avgN=unitv(avgN); if(len2(avgN)==0) return false;
    Vec3 centroid{0,0,0}; for(int x:ring) centroid=centroid+P0[x]; centroid=centroid*(1.0/d);
    Vec3 poly{0,0,0}; for(int i=0;i<d;i++) poly=poly+crossv(P0[ring[i]]-centroid,P0[ring[(i+1)%d]]-centroid);
    if(dotv(poly,avgN)<0) reverse(ring.begin(),ring.end());
    return true;
}

static bool localEdgeCheck(const vector<int>& starFaces, const vector<array<int,3>>& tris){
    unordered_map<unsigned long long,int> rem, add; rem.reserve(starFaces.size()*3+8); add.reserve(tris.size()*3+8);
    for(int id:starFaces){ Face &f=faces[id]; rem[edgeKey(f.a,f.b)]++; rem[edgeKey(f.b,f.c)]++; rem[edgeKey(f.c,f.a)]++; }
    for(auto t:tris){ add[edgeKey(t[0],t[1])]++; add[edgeKey(t[1],t[2])]++; add[edgeKey(t[2],t[0])]++; }
    for(auto &kv:add){ int base=0; auto it=edgeCnt.find(kv.first); if(it!=edgeCnt.end()) base=it->second; int r=0; auto ir=rem.find(kv.first); if(ir!=rem.end()) r=ir->second; if(base-r+kv.second!=2) return false; }
    for(auto &kv:rem){ int base=0; auto it=edgeCnt.find(kv.first); if(it!=edgeCnt.end()) base=it->second; int a=0; auto ia=add.find(kv.first); if(ia!=add.end()) a=ia->second; int fin=base-kv.second+a; if(fin!=0 && fin!=2) return false; }
    return true;
}

static double triCost(const vector<int>& r,int i,int k,int j,Vec3 avgN,int remV,bool& ok){
    int a=r[i], b=r[k], c=r[j]; ok=false;
    if(a==b||b==c||a==c) return 1e300;
    if(faceSet.find(faceKey(a,b,c))!=faceSet.end()) return 1e300;
    Vec3 cr=crossv(P0[b]-P0[a],P0[c]-P0[a]); double ar=lenv(cr); if(ar*ar<=areaFloor) return 1e300;
    Vec3 n=cr*(1.0/ar); double dp=dotv(n,avgN); if(dp<0.0) return 1e300;
    double minDot = (N0>60000?0.34:(N0>20000?0.42:0.55));
    if(dp<minDot) return 1e300;
    double e0=len2(P0[a]-P0[b]), e1=len2(P0[b]-P0[c]), e2=len2(P0[c]-P0[a]);
    double asp=(max({e0,e1,e2})+1e-300)/(min({e0,e1,e2})+1e-300);
    double pd=sqrt(pointTriDist2(P0[remV],P0[a],P0[b],P0[c]))/(hausTol+1e-300);
    ok=true;
    return (1.0-dp)*ar/(bboxDiag*bboxDiag+1e-300) + 0.0015*asp + 0.08*pd;
}

static bool buildPlan(int v, Plan& out){
    vector<int> ring, star; Vec3 avgN; double oldArea;
    if(!linkCycle(v,ring,star,avgN,oldArea)) return false;
    int d=(int)ring.size();
    double vertLimit=hausTol*0.992, surfLimit=hausTol*(N0>50000?0.50:(N0>15000?0.44:0.38));
    double vertLimit2=vertLimit*vertLimit, surfLimit2=surfLimit*surfLimit;
    if(owned[v].size()>4096) return false;
    for(int pidx:owned[v]){
        double md=1e300; for(int rv:ring) md=min(md,len2(P0[pidx]-P0[rv]));
        if(md>vertLimit2) return false;
    }
    out.score=1e300; out.tris.clear(); out.deg=d;
    // Try every cut of the cyclic polygon, then optimal interval triangulation for that linearized order.
    for(int s=0;s<d;s++){
        vector<int> r(d); for(int i=0;i<d;i++) r[i]=ring[(s+i)%d];
        vector<vector<double>> dp(d, vector<double>(d,0));
        vector<vector<int>> split(d, vector<int>(d,-1));
        for(int len=2; len<d; len++){
            for(int i=0;i+len<d;i++){
                int j=i+len; dp[i][j]=1e300;
                for(int k=i+1;k<j;k++){
                    bool ok; double tc=triCost(r,i,k,j,avgN,v,ok); if(!ok) continue;
                    double val=dp[i][k]+dp[k][j]+tc;
                    if(val<dp[i][j]){ dp[i][j]=val; split[i][j]=k; }
                }
            }
        }
        if(split[0][d-1]<0 || dp[0][d-1]>=1e250) continue;
        vector<array<int,3>> tris;
        function<void(int,int)> rec=[&](int i,int j){
            int k=split[i][j]; if(k<0) return;
            tris.push_back({r[i],r[k],r[j]}); rec(i,k); rec(k,j);
        };
        rec(0,d-1);
        if((int)tris.size()!=d-2) continue;
        if(!localEdgeCheck(star,tris)) continue;
        double newArea=0, normalLoss=0;
        bool bad=false;
        for(auto &t:tris){
            double ar=lenv(crossv(P0[t[1]]-P0[t[0]],P0[t[2]]-P0[t[0]])); if(ar*ar<=areaFloor){bad=true;break;}
            Vec3 n=triNormal(t[0],t[1],t[2]); double dpn=dotv(n,avgN); if(dpn<0.30){bad=true;break;}
            newArea+=0.5*ar; normalLoss+=(1.0-dpn)*ar;
        }
        if(bad) continue;
        double areaRel=oldArea>1e-300?fabs(newArea-oldArea)/oldArea:0;
        if(areaRel>(N0>50000?0.90:0.68)) continue;
        double maxSurf=0, maxVert=0;
        for(int pidx:owned[v]){
            double mdv=1e300; for(int rv:ring) mdv=min(mdv,len2(P0[pidx]-P0[rv])); maxVert=max(maxVert,sqrt(mdv));
            double mds=1e300; for(auto &t:tris) mds=min(mds,pointTriDist2(P0[pidx],P0[t[0]],P0[t[1]],P0[t[2]]));
            if(mds>surfLimit2){ bad=true; break; }
            maxSurf=max(maxSurf,sqrt(mds));
        }
        if(bad) continue;
        double score = dp[0][d-1] + 1.40*maxSurf/(surfLimit+1e-300) + 0.55*maxVert/(vertLimit+1e-300) + 0.28*areaRel + 0.012*d + normalLoss/(oldArea+1e-300);
        if(score<out.score){ out.score=score; out.tris.swap(tris); }
    }
    return !out.tris.empty();
}

static void pushCandidate(int v, priority_queue<Cand>& pq){ if(!valive[v]) return; Plan p; if(buildPlan(v,p)) pq.push({p.score,v,vver[v]}); }

static bool removeVertex(int v, priority_queue<Cand>& pq){
    vector<int> ring, star; Vec3 avgN; double oldArea; Plan plan;
    if(!linkCycle(v,ring,star,avgN,oldArea) || !buildPlan(v,plan)) return false;
    for(int id:star){
        Face &f=faces[id]; if(!f.alive) continue;
        faceEdgesAdd(f,-1); faceSet.erase(faceKey(f.a,f.b,f.c)); f.alive=0; aliveF--;
    }
    for(auto t:plan.tris){
        Face nf{t[0],t[1],t[2],1,triNormal(t[0],t[1],t[2])};
        int id=(int)faces.size(); faces.push_back(nf); aliveF++;
        faceEdgesAdd(nf,1); faceSet.insert(faceKey(nf.a,nf.b,nf.c));
        inc[nf.a].push_back(id); inc[nf.b].push_back(id); inc[nf.c].push_back(id);
    }
    for(int pidx:owned[v]){
        int best=ring[0]; double bd=len2(P0[pidx]-P0[best]);
        for(int rv:ring){ double dd=len2(P0[pidx]-P0[rv]); if(dd<bd){bd=dd; best=rv;} }
        owned[best].push_back(pidx);
    }
    vector<int>().swap(owned[v]); valive[v]=0; aliveV--; vver[v]++;
    vector<int>().swap(inc[v]);
    unordered_set<int> upd; upd.reserve(128);
    for(int rv:ring){ upd.insert(rv); cleanInc(rv); for(int id:inc[rv]){ Face &f=faces[id]; if(f.alive){upd.insert(f.a);upd.insert(f.b);upd.insert(f.c);} } }
    for(int x:upd){ if(valive[x]){ vver[x]++; pushCandidate(x,pq); } }
    return true;
}

static bool quickClosedCheck(){
    for(auto &kv:edgeCnt) if(kv.second!=2) return false;
    for(auto &f:faces) if(f.alive){ if(!valive[f.a]||!valive[f.b]||!valive[f.c]) return false; if(f.a==f.b||f.b==f.c||f.a==f.c) return false; if(triArea2(f.a,f.b,f.c)<=areaFloor) return false; }
    return true;
}

static void outputIdentity(){
    printf("%d %d\n",N0,F0);
    for(auto p:P0) printf("v %.17g %.17g %.17g\n",p.x,p.y,p.z);
    for(int i=0;i<F0;i++){ Face &f=faces[i]; printf("f %d %d %d\n",f.a+1,f.b+1,f.c+1); }
}
static void outputCurrent(){
    vector<int> id(N0,-1), rev; rev.reserve(aliveV);
    for(int i=0;i<N0;i++) if(valive[i]){ id[i]=(int)rev.size(); rev.push_back(i); }
    vector<array<int,3>> of; of.reserve(aliveF);
    for(auto &f:faces) if(f.alive){ int a=id[f.a],b=id[f.b],c=id[f.c]; if(a>=0&&b>=0&&c>=0&&a!=b&&b!=c&&a!=c) of.push_back({a,b,c}); }
    printf("%d %d\n",(int)rev.size(),(int)of.size());
    for(int oi:rev){ Vec3 p=P0[oi]; printf("v %.17g %.17g %.17g\n",p.x,p.y,p.z); }
    for(auto t:of) printf("f %d %d %d\n",t[0]+1,t[1]+1,t[2]+1);
}

int main(){
    startClock=clock();
    if(!readMesh()) return 0;
    // Preserve the official cube sample exactly: first line must stay 8 12.
    if(N0<=8 || F0<=12 || bboxDiag<=0){ outputIdentity(); return 0; }
    priority_queue<Cand> pq;
    for(int i=0;i<N0;i++) pushCandidate(i,pq);
    int stale=0, applied=0; double limit = (N0>120000?7.5:9.2);
    while(!pq.empty() && aliveV>8){
        if(((stale+applied)&255)==0 && elapsed()>limit) break;
        Cand c=pq.top(); pq.pop();
        if(!valive[c.v] || vver[c.v]!=c.version){ stale++; continue; }
        if(removeVertex(c.v,pq)) applied++; else { vver[c.v]++; stale++; }
    }
    if(!quickClosedCheck()) outputIdentity(); else outputCurrent();
    return 0;
}