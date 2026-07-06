#include <bits/stdc++.h>
using namespace std;

struct Vec3{
    double x,y,z;
    Vec3():x(0),y(0),z(0){}
    Vec3(double X,double Y,double Z):x(X),y(Y),z(Z){}
    Vec3 operator+(const Vec3& o) const { return {x+o.x,y+o.y,z+o.z}; }
    Vec3 operator-(const Vec3& o) const { return {x-o.x,y-o.y,z-o.z}; }
    Vec3 operator*(double s) const { return {x*s,y*s,z*s}; }
    Vec3 operator/(double s) const { return {x/s,y/s,z/s}; }
};
static inline double dotv(const Vec3&a,const Vec3&b){return a.x*b.x+a.y*b.y+a.z*b.z;}
static inline Vec3 crossv(const Vec3&a,const Vec3&b){return {a.y*b.z-a.z*b.y,a.z*b.x-a.x*b.z,a.x*b.y-a.y*b.x};}
static inline double norm2(const Vec3&a){return dotv(a,a);} 
static inline double normv(const Vec3&a){return sqrt(norm2(a));}

struct Quadric{
    double q[10];
    Quadric(){ memset(q,0,sizeof(q)); }
    inline void addPlane(double a,double b,double c,double d,double w){
        q[0]+=w*a*a; q[1]+=w*a*b; q[2]+=w*a*c; q[3]+=w*a*d;
        q[4]+=w*b*b; q[5]+=w*b*c; q[6]+=w*b*d;
        q[7]+=w*c*c; q[8]+=w*c*d; q[9]+=w*d*d;
    }
    inline void add(const Quadric& o){ for(int i=0;i<10;i++) q[i]+=o.q[i]; }
    inline double eval(const Vec3&p) const{
        double x=p.x,y=p.y,z=p.z;
        return q[0]*x*x + 2*q[1]*x*y + 2*q[2]*x*z + 2*q[3]*x
             + q[4]*y*y + 2*q[5]*y*z + 2*q[6]*y
             + q[7]*z*z + 2*q[8]*z + q[9];
    }
};

struct Face{ int v[3]; unsigned char alive; };
struct BBox{ double mn[3], mx[3]; };

static int N,FN;
static vector<Vec3> P;
static vector<Face> Faces;
static vector<vector<int>> Inc;
static vector<Quadric> Q;
static vector<BBox> BB;
static vector<unsigned char> VAlive;
static vector<int> Ver;
static int aliveV, aliveF;
static double haus2, diagLen;
static vector<int> markV;
static int curMark=1;

static vector<char> slurp_stdin(){
    vector<char> buf; buf.reserve(1<<26);
    char tmp[1<<16]; size_t n;
    while((n=fread(tmp,1,sizeof(tmp),stdin))>0) buf.insert(buf.end(), tmp, tmp+n);
    buf.push_back('\0'); return buf;
}

static void load_obj(){
    vector<char> buf=slurp_stdin(); char* p=buf.data();
    N=(int)strtol(p,&p,10); FN=(int)strtol(p,&p,10);
    P.resize(N); Faces.resize(FN); Inc.assign(N, {}); Q.assign(N, Quadric()); BB.resize(N); VAlive.assign(N,1); Ver.assign(N,0); markV.assign(N,0);
    double mnx=1e100,mny=1e100,mnz=1e100,mxx=-1e100,mxy=-1e100,mxz=-1e100;
    for(int i=0;i<N;i++){
        while(*p==' '||*p=='\n'||*p=='\r'||*p=='\t') ++p;
        if(*p=='v') ++p;
        double x=strtod(p,&p), y=strtod(p,&p), z=strtod(p,&p);
        P[i]={x,y,z};
        BB[i].mn[0]=BB[i].mx[0]=x; BB[i].mn[1]=BB[i].mx[1]=y; BB[i].mn[2]=BB[i].mx[2]=z;
        mnx=min(mnx,x); mny=min(mny,y); mnz=min(mnz,z); mxx=max(mxx,x); mxy=max(mxy,y); mxz=max(mxz,z);
    }
    for(int i=0;i<FN;i++){
        while(*p==' '||*p=='\n'||*p=='\r'||*p=='\t') ++p;
        if(*p=='f') ++p;
        int a=(int)strtol(p,&p,10)-1, b=(int)strtol(p,&p,10)-1, c=(int)strtol(p,&p,10)-1;
        Faces[i].v[0]=a; Faces[i].v[1]=b; Faces[i].v[2]=c; Faces[i].alive=1;
        Inc[a].push_back(i); Inc[b].push_back(i); Inc[c].push_back(i);
    }
    diagLen=sqrt((mxx-mnx)*(mxx-mnx)+(mxy-mny)*(mxy-mny)+(mxz-mnz)*(mxz-mnz));
    double h=0.0495*diagLen;
    haus2=h*h;
    aliveV=N; aliveF=FN;
}

static void init_quadrics(){
    for(int i=0;i<FN;i++){
        int a=Faces[i].v[0], b=Faces[i].v[1], c=Faces[i].v[2];
        Vec3 e1=P[b]-P[a], e2=P[c]-P[a];
        Vec3 n=crossv(e1,e2);
        double len=normv(n);
        if(len<=0) continue;
        Vec3 un=n/len;
        double d=-dotv(un,P[a]);
        double area=0.5*len;
        double w=max(area, 1e-18);
        Q[a].addPlane(un.x,un.y,un.z,d,w);
        Q[b].addPlane(un.x,un.y,un.z,d,w);
        Q[c].addPlane(un.x,un.y,un.z,d,w);
    }
}

static inline BBox mergeBB(const BBox&a,const BBox&b){
    BBox r;
    for(int k=0;k<3;k++){ r.mn[k]=min(a.mn[k],b.mn[k]); r.mx[k]=max(a.mx[k],b.mx[k]); }
    return r;
}

static inline bool insideHaus(const Vec3&p,const BBox&b){
    double md=0;
    for(int ix=0;ix<2;ix++) for(int iy=0;iy<2;iy++) for(int iz=0;iz<2;iz++){
        double x=ix?b.mx[0]:b.mn[0], y=iy?b.mx[1]:b.mn[1], z=iz?b.mx[2]:b.mn[2];
        double dx=p.x-x, dy=p.y-y, dz=p.z-z;
        double d=dx*dx+dy*dy+dz*dz;
        if(d>md) md=d;
    }
    return md <= haus2*(1.0000001);
}

static bool solveOptimal(const Quadric&q, Vec3& out){
    double a00=q.q[0], a01=q.q[1], a02=q.q[2];
    double a10=q.q[1], a11=q.q[4], a12=q.q[5];
    double a20=q.q[2], a21=q.q[5], a22=q.q[7];
    double b0=-q.q[3], b1=-q.q[6], b2=-q.q[8];
    double det = a00*(a11*a22-a12*a21) - a01*(a10*a22-a12*a20) + a02*(a10*a21-a11*a20);
    if(fabs(det) < 1e-18) return false;
    double invDet=1.0/det;
    double x = (b0*(a11*a22-a12*a21) - a01*(b1*a22-a12*b2) + a02*(b1*a21-a11*b2))*invDet;
    double y = (a00*(b1*a22-a12*b2) - b0*(a10*a22-a12*a20) + a02*(a10*b2-b1*a20))*invDet;
    double z = (a00*(a11*b2-b1*a21) - a01*(a10*b2-b1*a20) + b0*(a10*a21-a11*a20))*invDet;
    if(!isfinite(x)||!isfinite(y)||!isfinite(z)) return false;
    out={x,y,z}; return true;
}

static inline bool faceHas(const Face&f,int v){ return f.v[0]==v || f.v[1]==v || f.v[2]==v; }

static bool edgeFaces(int u,int v,int opp[2],int adjF[2]){
    int cnt=0;
    for(int fid: Inc[u]){
        const Face& f=Faces[fid]; if(!f.alive) continue;
        bool hu=false,hv=false; int op=-1;
        for(int k=0;k<3;k++){ int w=f.v[k]; if(w==u) hu=true; else if(w==v) hv=true; else op=w; }
        if(hu && hv){
            if(cnt<2){ opp[cnt]=op; adjF[cnt]=fid; }
            cnt++;
            if(cnt>2) return false;
        }
    }
    return cnt==2 && opp[0]!=opp[1];
}

static bool linkOK(int u,int v,const int opp[2]){
    if(++curMark > INT_MAX-4){ fill(markV.begin(), markV.end(), 0); curMark=1; }
    int tag=curMark++;
    for(int fid: Inc[u]){
        const Face& f=Faces[fid]; if(!f.alive) continue;
        if(!faceHas(f,u)) continue;
        for(int k=0;k<3;k++){ int w=f.v[k]; if(w!=u && w!=v) markV[w]=tag; }
    }
    int tag2=curMark++;
    int common[4]; int cnt=0;
    for(int fid: Inc[v]){
        const Face& f=Faces[fid]; if(!f.alive) continue;
        if(!faceHas(f,v)) continue;
        for(int k=0;k<3;k++){
            int w=f.v[k];
            if(w==u || w==v) continue;
            if(markV[w]==tag){
                markV[w]=tag2;
                if(cnt<4) common[cnt]=w;
                cnt++;
                if(cnt>2) return false;
            }
        }
    }
    if(cnt!=2) return false;
    return ( (common[0]==opp[0] && common[1]==opp[1]) || (common[0]==opp[1] && common[1]==opp[0]) );
}

static Vec3 faceNormalCurrent(const Face&f){
    return crossv(P[f.v[1]]-P[f.v[0]], P[f.v[2]]-P[f.v[0]]);
}

static bool geometryOK(int u,int v,const Vec3&np,const int adjF[2]){
    const double minArea2 = 1e-28;
    int af0=adjF[0], af1=adjF[1];
    if(++curMark > INT_MAX-4){ fill(markV.begin(), markV.end(), 0); curMark=1; }
    auto checkFace=[&](int fid)->bool{
        if(fid==af0 || fid==af1) return true;
        const Face& f=Faces[fid]; if(!f.alive) return true;
        bool affected=false;
        Vec3 a,b,c;
        int ids[3]={f.v[0],f.v[1],f.v[2]};
        for(int k=0;k<3;k++) if(ids[k]==u || ids[k]==v) affected=true;
        if(!affected) return true;
        a = (ids[0]==u || ids[0]==v) ? np : P[ids[0]];
        b = (ids[1]==u || ids[1]==v) ? np : P[ids[1]];
        c = (ids[2]==u || ids[2]==v) ? np : P[ids[2]];
        Vec3 nn=crossv(b-a,c-a);
        double n2=norm2(nn);
        if(n2 <= minArea2) return false;
        Vec3 on=faceNormalCurrent(f);
        double od=dotv(on,nn);
        if(od <= -1e-18) return false;
        return true;
    };
    for(int fid: Inc[u]) if(!checkFace(fid)) return false;
    for(int fid: Inc[v]) if(!checkFace(fid)) return false;
    return true;
}

struct Cand{ bool ok; Vec3 p; double cost; };

static double edgeSharpness(const int adjF[2]){
    Vec3 n1=faceNormalCurrent(Faces[adjF[0]]), n2=faceNormalCurrent(Faces[adjF[1]]);
    double l1=normv(n1), l2=normv(n2); if(l1<=0||l2<=0) return 10.0;
    n1=n1/l1; n2=n2/l2;
    double d=max(-1.0,min(1.0,dotv(n1,n2)));
    double s=1.0-d;
    double sig=0.0;
    double a1[3]={n1.x,n1.y,n1.z}, a2[3]={n2.x,n2.y,n2.z};
    for(int k=0;k<3;k++){
        if(a1[k]*a2[k] < 0) sig=max(sig, min(1.0, fabs(a1[k]-a2[k])));
    }
    return s + 0.35*sig;
}

static Cand computeCandidate(int u,int v,bool needGeom=false){
    Cand res; res.ok=false; res.cost=1e100; res.p=Vec3();
    if(!VAlive[u]||!VAlive[v]||u==v) return res;
    int opp[2], adjF[2];
    if(!edgeFaces(u,v,opp,adjF)) return res;
    if(needGeom && !linkOK(u,v,opp)) return res;
    BBox mb=mergeBB(BB[u],BB[v]);
    Quadric qq=Q[u]; qq.add(Q[v]);
    Vec3 cand[5]; int cn=0;
    cand[cn++]=P[u]; cand[cn++]=P[v]; cand[cn++]=(P[u]+P[v])*0.5;
    Vec3 opt;
    if(solveOptimal(qq,opt)) cand[cn++]=opt;
    Vec3 e=P[v]-P[u]; double el2=norm2(e);
    if(el2>1e-30){
        double f0=qq.eval(P[u]), f1=qq.eval(P[v]), fm=qq.eval((P[u]+P[v])*0.5);
        double A=2*(f1+f0-2*fm);
        double B=f1-f0-A;
        if(fabs(A)>1e-30){ double t=-B/(2*A); if(t>0 && t<1) cand[cn++]=P[u]+e*t; }
    }
    double sharp=0.0;
    if(adjF[0]>=0) sharp=edgeSharpness(adjF);
    double el=sqrt(max(0.0,norm2(P[u]-P[v])));
    for(int i=0;i<cn;i++){
        Vec3 p=cand[i];
        if(!insideHaus(p,mb)) continue;
        if(needGeom && !geometryOK(u,v,p,adjF)) continue;
        double c=qq.eval(p);
        double penalty = 1.0 + 80.0*sharp*sharp;
        c = c*penalty + 1e-12*el*el*penalty;
        if(c<res.cost){ res.cost=c; res.p=p; res.ok=true; }
    }
    return res;
}

struct Node{
    double cost; int u,v; int vu,vv;
    bool operator<(const Node&o) const { return cost > o.cost; }
};
static priority_queue<Node> PQ;

static void pushEdge(int u,int v){
    if(u==v || !VAlive[u] || !VAlive[v]) return;
    Cand c=computeCandidate(u,v,false);
    if(!c.ok) return;
    PQ.push(Node{c.cost,u,v,Ver[u],Ver[v]});
}

static double initialCurvP90 = 0.0;

static Vec3 unitFaceNormalById(int fid){
    const Face& f=Faces[fid];
    Vec3 n=crossv(P[f.v[1]]-P[f.v[0]], P[f.v[2]]-P[f.v[0]]);
    double l=normv(n);
    if(l<=0) return Vec3(0,0,1);
    return n/l;
}

static inline uint64_t packEdgeFace(int a,int b,int fid){
    if(a>b) swap(a,b);
    return ((uint64_t)(uint32_t)a << 43) | ((uint64_t)(uint32_t)b << 22) | (uint64_t)(uint32_t)fid;
}
static inline int unpackA(uint64_t r){ return (int)(r >> 43); }
static inline int unpackB(uint64_t r){ return (int)((r >> 22) & ((1u<<21)-1)); }
static inline int unpackFid(uint64_t r){ return (int)(r & ((1u<<22)-1)); }
static inline uint64_t packedEdgeOnly(uint64_t r){ return r >> 22; }

static void pushInitialEdges(){
    vector<uint64_t> edges; edges.reserve((size_t)FN*3);
    for(int i=0;i<FN;i++){
        int a=Faces[i].v[0], b=Faces[i].v[1], c=Faces[i].v[2];
        edges.push_back(packEdgeFace(a,b,i));
        edges.push_back(packEdgeFace(b,c,i));
        edges.push_back(packEdgeFace(c,a,i));
    }
    sort(edges.begin(), edges.end());
    const int BINS=2048;
    vector<int> hist(BINS,0);
    long long manifoldEdges=0;
    for(size_t i=0;i<edges.size();){
        size_t j=i+1; uint64_t ek=packedEdgeOnly(edges[i]);
        while(j<edges.size() && packedEdgeOnly(edges[j])==ek) ++j;
        int a=unpackA(edges[i]), b=unpackB(edges[i]);
        pushEdge(a,b);
        if(j==i+2){
            Vec3 n1=unitFaceNormalById(unpackFid(edges[i]));
            Vec3 n2=unitFaceNormalById(unpackFid(edges[i+1]));
            double d=max(-1.0,min(1.0,dotv(n1,n2)));
            double s=1.0-d;
            int bin=(int)floor((min(2.0,max(0.0,s))/2.0)*(BINS-1));
            hist[bin]++; manifoldEdges++;
        }
        i=j;
    }
    if(manifoldEdges>0){
        long long need=(long long)ceil(0.90*manifoldEdges), acc=0;
        int bin=0;
        for(;bin<BINS;bin++){ acc+=hist[bin]; if(acc>=need) break; }
        initialCurvP90 = 2.0*(bin+0.5)/BINS;
    }
}

static void cleanInc(int u){
    vector<int> nv; nv.reserve(Inc[u].size());
    for(int fid: Inc[u]){
        if(!Faces[fid].alive) continue;
        const Face& f=Faces[fid];
        if(f.v[0]==u || f.v[1]==u || f.v[2]==u) nv.push_back(fid);
    }
    Inc[u].swap(nv);
}

static void collectAndPushNeighbors(int u){
    if(++curMark > INT_MAX-4){ fill(markV.begin(), markV.end(), 0); curMark=1; }
    int tag=curMark++;
    vector<int> neigh;
    for(int fid: Inc[u]){
        const Face& f=Faces[fid]; if(!f.alive) continue;
        for(int k=0;k<3;k++){
            int w=f.v[k]; if(w==u || !VAlive[w]) continue;
            if(markV[w]!=tag){ markV[w]=tag; neigh.push_back(w); }
        }
    }
    for(int w: neigh) pushEdge(u,w);
}

static bool doCollapse(int u,int v,const Vec3&np){
    int opp[2], adjF[2];
    if(!edgeFaces(u,v,opp,adjF)) return false;
    if(!linkOK(u,v,opp)) return false;
    if(!geometryOK(u,v,np,adjF)) return false;
    int keep=u, rem=v;
    if(Inc[v].size() > Inc[u].size()) { keep=v; rem=u; }
    BBox mb=mergeBB(BB[keep],BB[rem]);
    if(!insideHaus(np,mb)) return false;
    P[keep]=np;
    BB[keep]=mb;
    Q[keep].add(Q[rem]);
    VAlive[rem]=0; aliveV--; Ver[keep]++; Ver[rem]++;
    for(int fid: Inc[rem]){
        Face& f=Faces[fid]; if(!f.alive) continue;
        bool hasKeep=false, hasRem=false;
        for(int k=0;k<3;k++){ if(f.v[k]==keep) hasKeep=true; if(f.v[k]==rem) hasRem=true; }
        if(!hasRem) continue;
        if(hasKeep){ f.alive=0; aliveF--; }
        else{
            for(int k=0;k<3;k++) if(f.v[k]==rem) f.v[k]=keep;
            Inc[keep].push_back(fid);
        }
    }
    Inc[rem].clear();
    cleanInc(keep);
    collectAndPushNeighbors(keep);
    return true;
}

static double targetRatio(){
    if(N < 20) return 0.05;
    double base;
    if(N < 1000) base = 0.13;
    else if(N <= 6000) base = 0.080;
    else if(N <= 30000) base = 0.050;
    else if(N <= 60000) base = 0.047;
    else if(N <= 150000) base = 0.042;
    else if(N <= 600000) base = 0.035;
    else base = 0.032;
    double add = 0.24 * sqrt(max(0.0, initialCurvP90));
    if(add > 0.080) add = 0.080;
    double r = base + add;
    if(r > 0.18) r = 0.18;
    if(r < 0.025) r = 0.025;
    return r;
}

static void simplify(){
    if(N<=4) return;
    init_quadrics();
    pushInitialEdges();
    int target=max(4, (int)ceil(N*targetRatio()));
    while(aliveV>target && !PQ.empty()){
        Node nd=PQ.top(); PQ.pop();
        if(!VAlive[nd.u]||!VAlive[nd.v]) continue;
        if(Ver[nd.u]!=nd.vu || Ver[nd.v]!=nd.vv) continue;
        Cand c=computeCandidate(nd.u,nd.v,true);
        if(!c.ok) continue;
        doCollapse(nd.u, nd.v, c.p);
    }
}

static void save_obj(){
    vector<int> remap(N,-1); int nv=0;
    for(int i=0;i<N;i++) if(VAlive[i]) remap[i]=nv++;
    vector<array<int,3>> of;
    of.reserve(aliveF>0?aliveF:0);
    for(int i=0;i<FN;i++) if(Faces[i].alive){
        int a=remap[Faces[i].v[0]], b=remap[Faces[i].v[1]], c=remap[Faces[i].v[2]];
        if(a>=0 && b>=0 && c>=0 && a!=b && b!=c && c!=a) of.push_back({a,b,c});
    }
    string out; out.reserve((size_t)nv*48 + (size_t)of.size()*28 + 64);
    char line[128];
    out.append(line, snprintf(line,sizeof(line), "%d %d\n", nv, (int)of.size()));
    for(int i=0;i<N;i++) if(VAlive[i]){
        out.append(line, snprintf(line,sizeof(line), "v %.10g %.10g %.10g\n", P[i].x,P[i].y,P[i].z));
    }
    for(auto &f: of){
        out.append(line, snprintf(line,sizeof(line), "f %d %d %d\n", f[0]+1,f[1]+1,f[2]+1));
    }
    fwrite(out.data(),1,out.size(),stdout);
}

int main(){
    load_obj();
    simplify();
    save_obj();
    return 0;
}
