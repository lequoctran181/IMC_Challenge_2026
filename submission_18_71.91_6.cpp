#include <bits/stdc++.h>
using namespace std;

// Perception-aware mesh simplifier for IMC simplifygeometry.
// QEM manifold edge collapses + internal multi-view render/SSIM checkpoint selection.
// Self-contained C++17, no Eigen dependency.

struct FastScanner {
    static const size_t BUFSIZE = 1 << 20;
    int idx = 0, size = 0;
    char buf[BUFSIZE];
    inline char getch() {
        if (idx >= size) {
            size = (int)fread(buf, 1, BUFSIZE, stdin);
            idx = 0;
            if (!size) return 0;
        }
        return buf[idx++];
    }
    template<class T> bool readInt(T &out) {
        char c = getch();
        while (c && (c==' ' || c=='\n' || c=='\r' || c=='\t')) c = getch();
        if (!c) return false;
        T sign = 1, x = 0;
        if (c=='-') sign = -1, c = getch();
        while (c>='0' && c<='9') { x = x*10 + (c-'0'); c = getch(); }
        out = x * sign;
        return true;
    }
    inline bool readPrefix(char want) {
        char c = getch();
        while (c && (c==' ' || c=='\n' || c=='\r' || c=='\t')) c = getch();
        return c == want;
    }
    bool readDouble(double &out) {
        char c = getch();
        while (c && (c==' ' || c=='\n' || c=='\r' || c=='\t')) c = getch();
        if (!c) return false;
        int sign = 1;
        if (c=='-') sign = -1, c = getch();
        double x = 0.0;
        while (c>='0' && c<='9') { x = x*10.0 + (double)(c-'0'); c = getch(); }
        if (c=='.') {
            double base = 0.1; c = getch();
            while (c>='0' && c<='9') { x += base * (double)(c-'0'); base *= 0.1; c = getch(); }
        }
        if (c=='e' || c=='E') {
            int esign=1, expv=0; c = getch();
            if (c=='-') esign=-1, c=getch(); else if (c=='+') c=getch();
            while (c>='0' && c<='9') { expv = expv*10 + (c-'0'); c=getch(); }
            x *= pow(10.0, esign * expv);
        }
        out = sign * x;
        return true;
    }
};

struct Vec3 {
    double x, y, z;
    Vec3():x(0),y(0),z(0){}
    Vec3(double X,double Y,double Z):x(X),y(Y),z(Z){}
    inline Vec3 operator+(const Vec3& o) const { return Vec3(x+o.x,y+o.y,z+o.z); }
    inline Vec3 operator-(const Vec3& o) const { return Vec3(x-o.x,y-o.y,z-o.z); }
    inline Vec3 operator*(double s) const { return Vec3(x*s,y*s,z*s); }
    inline Vec3 operator/(double s) const { return Vec3(x/s,y/s,z/s); }
};
struct Tri { int a,b,c; };
static inline double dotv(const Vec3& a,const Vec3& b){ return a.x*b.x+a.y*b.y+a.z*b.z; }
static inline Vec3 crossv(const Vec3& a,const Vec3& b){ return Vec3(a.y*b.z-a.z*b.y, a.z*b.x-a.x*b.z, a.x*b.y-a.y*b.x); }
static inline double norm2(const Vec3& a){ return dotv(a,a); }
static inline double normv(const Vec3& a){ return sqrt(norm2(a)); }

struct Quadric {
    // xx xy xz xw yy yz yw zz zw ww
    double q[10];
    Quadric(){ memset(q,0,sizeof(q)); }
    inline void add(const Quadric& o){ for(int i=0;i<10;i++) q[i]+=o.q[i]; }
    static Quadric plane(double a,double b,double c,double d,double w) {
        Quadric Q;
        Q.q[0]=w*a*a; Q.q[1]=w*a*b; Q.q[2]=w*a*c; Q.q[3]=w*a*d;
        Q.q[4]=w*b*b; Q.q[5]=w*b*c; Q.q[6]=w*b*d;
        Q.q[7]=w*c*c; Q.q[8]=w*c*d; Q.q[9]=w*d*d;
        return Q;
    }
    inline double eval(const Vec3& p) const {
        double x=p.x,y=p.y,z=p.z;
        return q[0]*x*x + 2*q[1]*x*y + 2*q[2]*x*z + 2*q[3]*x
             + q[4]*y*y + 2*q[5]*y*z + 2*q[6]*y
             + q[7]*z*z + 2*q[8]*z + q[9];
    }
};
static inline Quadric qsum(const Quadric& a,const Quadric& b){ Quadric r=a; r.add(b); return r; }

struct Face {
    int v[3];
    unsigned char active;
    Vec3 n;
};
struct Candidate {
    double score;
    int u, v, vu, vv;
    bool operator<(Candidate const& o) const { return score > o.score; }
};
struct EvalCand { bool ok; double score; double rawError; Vec3 pos; };
struct Snapshot { vector<Vec3> V; vector<Tri> F; int nv=0; double ratio=1.0; double ssim=-1.0; };
struct RenderMaps { int R=0; vector<float> depth; vector<float> norm; };

static int N, M;
static vector<Vec3> P;
static vector<Face> F;
static vector<vector<int>> inc;
static vector<Quadric> Qv;
static vector<array<float,3>> bbMin, bbMax;
static vector<unsigned char> alive;
static vector<int> ver, markA, markB, tmpId;
static int stampA=1, stampB=1;
static double hausTau=0.0, hausTau2=0.0;
static int activeVertices=0, activeFaces=0;
static double visualPixelTol=18.0, visualPixelTol2=324.0;
static vector<Snapshot> snapshots;
static Snapshot finalSnap;
static bool haveFinalSnap=false;

static inline bool face_contains(int fid,int u){ const Face& f=F[fid]; return f.v[0]==u || f.v[1]==u || f.v[2]==u; }
static inline Vec3 computeFaceNormal(const Face& f){
    Vec3 cr = crossv(P[f.v[1]]-P[f.v[0]], P[f.v[2]]-P[f.v[0]]);
    double l=normv(cr); if(l<=1e-300) return Vec3(0,0,0); return cr/l;
}

static void cleanupIncident(int u){
    if(!alive[u]) return;
    vector<int>& lst=inc[u]; int wr=0;
    for(int fid: lst) if(fid>=0 && fid<M && F[fid].active && face_contains(fid,u)) lst[wr++]=fid;
    lst.resize(wr);
}
static void maybeCleanupIncident(int u){
    if(!alive[u]) return;
    if(inc[u].size()>96){
        int good=0; for(int fid: inc[u]) if(F[fid].active && face_contains(fid,u)) ++good;
        if((int)inc[u].size()>good*3+64) cleanupIncident(u);
    }
}
static inline void mergedBBox(int u,int v,array<float,3>& mn,array<float,3>& mx){
    mn[0]=min(bbMin[u][0],bbMin[v][0]); mn[1]=min(bbMin[u][1],bbMin[v][1]); mn[2]=min(bbMin[u][2],bbMin[v][2]);
    mx[0]=max(bbMax[u][0],bbMax[v][0]); mx[1]=max(bbMax[u][1],bbMax[v][1]); mx[2]=max(bbMax[u][2],bbMax[v][2]);
}
static inline bool coversBBox(const Vec3& p,const array<float,3>& mn,const array<float,3>& mx){
    double maxd2=0.0;
    for(int mask=0;mask<8;mask++){
        double x=(mask&1)?mx[0]:mn[0], y=(mask&2)?mx[1]:mn[1], z=(mask&4)?mx[2]:mn[2];
        double dx=p.x-x,dy=p.y-y,dz=p.z-z; double d2=dx*dx+dy*dy+dz*dz;
        if(d2>maxd2) maxd2=d2;
    }
    return maxd2 <= hausTau2;
}
static inline void project1024(const Vec3& p,int view,double& u,double& v,double& dep){
    int ax=view/2; int sg = (view&1)?-1:1; // +, - for x/y/z
    double sx, sy;
    if(ax==0){ sx=p.y; sy=p.z; dep=2.5 - sg*p.x; if(sg<0) sx=-sx; }
    else if(ax==1){ sx=p.x; sy=p.z; dep=2.5 - sg*p.y; if(sg>0) sx=-sx; }
    else { sx=p.x; sy=p.y; dep=2.5 - sg*p.z; if(sg<0) sx=-sx; }
    u = 800.0*sx/dep + 512.0; v = 800.0*sy/dep + 512.0;
}
static inline bool visualBBoxOK(const Vec3& p,const array<float,3>& mn,const array<float,3>& mx){
    // Conservative screen-space cluster radius in the judge's 1024px camera.  Planar zero-QEM collapses bypass this elsewhere.
    double pu[6], pv[6], pd;
    for(int view=0;view<6;view++) project1024(p,view,pu[view],pv[view],pd);
    double maxd2=0.0;
    for(int mask=0;mask<8;mask++){
        Vec3 c((mask&1)?mx[0]:mn[0], (mask&2)?mx[1]:mn[1], (mask&4)?mx[2]:mn[2]);
        for(int view=0;view<6;view++){
            double u,v,d; project1024(c,view,u,v,d);
            double du=u-pu[view], dv=v-pv[view]; double d2=du*du+dv*dv;
            if(d2>maxd2) maxd2=d2;
            if(maxd2>visualPixelTol2) return false;
        }
    }
    return true;
}

static bool solveOptimal(const Quadric& q,Vec3& out){
    double a00=q.q[0], a01=q.q[1], a02=q.q[2];
    double a11=q.q[4], a12=q.q[5], a22=q.q[7];
    double b0=q.q[3], b1=q.q[6], b2=q.q[8];
    double det=a00*(a11*a22-a12*a12)-a01*(a01*a22-a12*a02)+a02*(a01*a12-a11*a02);
    if(fabs(det)<1e-14) return false;
    double inv00=(a11*a22-a12*a12)/det;
    double inv01=(a02*a12-a01*a22)/det;
    double inv02=(a01*a12-a02*a11)/det;
    double inv11=(a00*a22-a02*a02)/det;
    double inv12=(a01*a02-a00*a12)/det;
    double inv22=(a00*a11-a01*a01)/det;
    out.x=-(inv00*b0+inv01*b1+inv02*b2);
    out.y=-(inv01*b0+inv11*b1+inv12*b2);
    out.z=-(inv02*b0+inv12*b1+inv22*b2);
    if(!isfinite(out.x)||!isfinite(out.y)||!isfinite(out.z)) return false;
    return out.x>=-2.0&&out.x<=2.0&&out.y>=-2.0&&out.y<=2.0&&out.z>=-2.0&&out.z<=2.0;
}
static EvalCand computeCandidate(int u,int v){
    EvalCand ec{false,1e300,1e300,Vec3()};
    if(u==v || !alive[u] || !alive[v]) return ec;
    if(u>v) swap(u,v);
    array<float,3> mn,mx; mergedBBox(u,v,mn,mx);
    Quadric q=qsum(Qv[u],Qv[v]);
    Vec3 cand[7]; int cc=0; Vec3 opt;
    if(solveOptimal(q,opt)) cand[cc++]=opt;
    cand[cc++]=P[u]; cand[cc++]=P[v];
    cand[cc++]=(P[u]+P[v])*0.5;
    cand[cc++]=Vec3((mn[0]+mx[0])*0.5,(mn[1]+mx[1])*0.5,(mn[2]+mx[2])*0.5);
    // two biased midpoints often keep silhouettes better than a pure QEM minimizer on long skinny triangles.
    cand[cc++]=P[u]*0.67 + P[v]*0.33;
    cand[cc++]=P[u]*0.33 + P[v]*0.67;
    double len2=norm2(P[u]-P[v]);
    for(int i=0;i<cc;i++){
        Vec3 p=cand[i];
        if(!coversBBox(p,mn,mx)) continue;
        double e=q.eval(p); if(e<0 && e>-1e-12) e=0;
        // Internal render/SSIM selection is a stronger guard than a hard projected-radius cut.
        // Keep all Hausdorff-valid QEM candidates here; ordering still favors short/local collapses.
        // if(e > 1e-11 && !visualBBoxOK(p,mn,mx)) continue;
        // QEM + projected/locality bias.  The bias is deliberately tiny except on near-planar ties.
        double sc=e + 2e-10*len2 + 5e-15*norm2(p);
        if(sc<ec.score){ ec.ok=true; ec.score=sc; ec.rawError=e; ec.pos=p; }
    }
    return ec;
}
static bool linkConditionOK(int u,int v){
    if(u==v || !alive[u] || !alive[v]) return false;
    maybeCleanupIncident(u); maybeCleanupIncident(v);
    if(++stampA==INT_MAX){ fill(markA.begin(),markA.end(),0); stampA=1; }
    if(++stampB==INT_MAX){ fill(markB.begin(),markB.end(),0); stampB=1; }
    int edgeFaces=0;
    for(int fid: inc[u]){
        Face& f=F[fid]; if(!f.active || !face_contains(fid,u)) continue;
        bool hasv=false;
        for(int k=0;k<3;k++){ int w=f.v[k]; if(w==v) hasv=true; if(w!=u) markA[w]=stampA; }
        if(hasv) ++edgeFaces;
    }
    if(edgeFaces!=2) return false;
    int common=0;
    for(int fid: inc[v]){
        Face& f=F[fid]; if(!f.active || !face_contains(fid,v)) continue;
        for(int k=0;k<3;k++){
            int w=f.v[k]; if(w==u || w==v) continue;
            if(markA[w]==stampA && markB[w]!=stampB){ markB[w]=stampB; if(++common>2) return false; }
        }
    }
    return common==2;
}
static bool flipOK(int keep,int rem,const Vec3& newPos,double minDot){
    static vector<int> touched; touched.clear();
    for(int fid: inc[keep]) if(F[fid].active && face_contains(fid,keep)) touched.push_back(fid);
    for(int fid: inc[rem]) if(F[fid].active && face_contains(fid,rem)) touched.push_back(fid);
    sort(touched.begin(),touched.end()); touched.erase(unique(touched.begin(),touched.end()),touched.end());
    const double minArea2=1e-28;
    for(int fid: touched){
        Face& f=F[fid]; if(!f.active) continue;
        bool hasK=false,hasR=false; for(int k=0;k<3;k++){ if(f.v[k]==keep) hasK=true; if(f.v[k]==rem) hasR=true; }
        if(hasK && hasR) continue;
        int a=f.v[0],b=f.v[1],c=f.v[2];
        Vec3 pa=(a==keep||a==rem)?newPos:P[a];
        Vec3 pb=(b==keep||b==rem)?newPos:P[b];
        Vec3 pc=(c==keep||c==rem)?newPos:P[c];
        Vec3 cr=crossv(pb-pa,pc-pa); double a2=norm2(cr);
        if(!(a2>minArea2) || !isfinite(a2)) return false;
        Vec3 nn=cr/sqrt(a2);
        if(dotv(nn,f.n)<minDot) return false;
    }
    return true;
}
static void pushEdge(priority_queue<Candidate>& pq,int u,int v){
    if(u==v || !alive[u] || !alive[v]) return; if(u>v) swap(u,v);
    EvalCand ec=computeCandidate(u,v); if(!ec.ok) return;
    pq.push({ec.score,u,v,ver[u],ver[v]});
}
static void collectNeighbors(int u,vector<int>& out){
    out.clear(); if(!alive[u]) return; maybeCleanupIncident(u);
    if(++stampA==INT_MAX){ fill(markA.begin(),markA.end(),0); stampA=1; }
    for(int fid: inc[u]){
        Face& f=F[fid]; if(!f.active || !face_contains(fid,u)) continue;
        for(int k=0;k<3;k++){ int w=f.v[k]; if(w!=u && alive[w] && markA[w]!=stampA){ markA[w]=stampA; out.push_back(w); } }
    }
}
static void doCollapse(int keep,int rem,const Vec3& newPos,priority_queue<Candidate>& pq){
    cleanupIncident(keep); cleanupIncident(rem);
    vector<int> remFaces=inc[rem];
    for(int fid: remFaces){
        Face& f=F[fid]; if(!f.active || !face_contains(fid,rem)) continue;
        bool hasKeep=false; for(int k=0;k<3;k++) if(f.v[k]==keep) hasKeep=true;
        if(hasKeep){ f.active=0; --activeFaces; }
        else{
            for(int k=0;k<3;k++) if(f.v[k]==rem) f.v[k]=keep;
            if(f.v[0]==f.v[1] || f.v[1]==f.v[2] || f.v[2]==f.v[0]){ f.active=0; --activeFaces; }
            else{ P[keep]=newPos; f.n=computeFaceNormal(f); inc[keep].push_back(fid); }
        }
    }
    P[keep]=newPos;
    for(int fid: inc[keep]){ Face& f=F[fid]; if(f.active && face_contains(fid,keep)) f.n=computeFaceNormal(f); }
    Qv[keep].add(Qv[rem]);
    for(int d=0;d<3;d++){ bbMin[keep][d]=min(bbMin[keep][d],bbMin[rem][d]); bbMax[keep][d]=max(bbMax[keep][d],bbMax[rem][d]); }
    alive[rem]=0; ++ver[keep]; ++ver[rem]; --activeVertices; inc[rem].clear(); cleanupIncident(keep);
    static vector<int> neigh; collectNeighbors(keep,neigh); for(int w: neigh) pushEdge(pq,keep,w);
}

static Snapshot makeSnapshot(double ratio){
    Snapshot s; s.ratio=ratio; s.nv=activeVertices;
    tmpId.assign(N,-1); s.V.reserve(activeVertices);
    for(int i=0;i<N;i++) if(alive[i]){ tmpId[i]=(int)s.V.size(); s.V.push_back(P[i]); }
    s.F.reserve(activeFaces);
    for(int i=0;i<M;i++) if(F[i].active){
        int a=F[i].v[0],b=F[i].v[1],c=F[i].v[2];
        if(a!=b && b!=c && c!=a && tmpId[a]>=0 && tmpId[b]>=0 && tmpId[c]>=0) s.F.push_back({tmpId[a],tmpId[b],tmpId[c]});
    }
    s.nv=(int)s.V.size();
    return s;
}

// ---------------- internal renderer and SSIM ----------------
static int evalR=512, evalPix=512*512;
static constexpr int VIEWS=6;

static void renderMesh(const vector<Vec3>& verts,const vector<Tri>& tris,RenderMaps& maps,int R){
    int PIX=R*R; double focal=800.0*(double)R/1024.0, center=0.5*(double)R;
    maps.R=R; maps.depth.assign((size_t)VIEWS*PIX,255.0f); maps.norm.assign((size_t)VIEWS*PIX*3,127.5f);
    int nv=(int)verts.size(); int nt=(int)tris.size();
    vector<float> U(nv), VV(nv), Z(nv);
    vector<Vec3> fn(nt);
    for(int i=0;i<nt;i++){
        const Tri& t=tris[i]; Vec3 cr=crossv(verts[t.b]-verts[t.a],verts[t.c]-verts[t.a]); double l=normv(cr);
        fn[i]=(l>1e-300)?cr/l:Vec3(0,0,0);
    }
    for(int view=0;view<VIEWS;view++){
        int ax=view/2; int sg=(view&1)?-1:1;
        for(int i=0;i<nv;i++){
            const Vec3& p=verts[i]; double sx,sy,dep;
            if(ax==0){ sx=p.y; sy=p.z; dep=2.5 - sg*p.x; if(sg<0) sx=-sx; }
            else if(ax==1){ sx=p.x; sy=p.z; dep=2.5 - sg*p.y; if(sg>0) sx=-sx; }
            else { sx=p.x; sy=p.y; dep=2.5 - sg*p.z; if(sg<0) sx=-sx; }
            Z[i]=(float)dep; U[i]=(float)(focal*sx/dep + center); VV[i]=(float)(focal*sy/dep + center);
        }
        float* zbuf = maps.depth.data() + (size_t)view*PIX;
        float* nbuf = maps.norm.data() + (size_t)view*PIX*3;
        for(int ti=0;ti<nt;ti++){
            const Tri& t=tris[ti]; int ia=t.a,ib=t.b,ic=t.c;
            float u0=U[ia],u1=U[ib],u2=U[ic]; float v0=VV[ia],v1=VV[ib],v2=VV[ic];
            float z0=Z[ia],z1=Z[ib],z2=Z[ic];
            if(!(z0>0 && z1>0 && z2>0)) continue;
            float minx=min(u0,min(u1,u2)), maxx=max(u0,max(u1,u2));
            float miny=min(v0,min(v1,v2)), maxy=max(v0,max(v1,v2));
            int x0i=max(0,(int)floor(minx)); int x1i=min(R-1,(int)ceil(maxx));
            int y0i=max(0,(int)floor(miny)); int y1i=min(R-1,(int)ceil(maxy));
            if(x0i>x1i || y0i>y1i) continue;
            float den=(v1-v2)*(u0-u2)+(u2-u1)*(v0-v2);
            if(fabs(den)<1e-20f) continue;
            float invDen=1.0f/den;
            Vec3 n=fn[ti];
            float nr=(float)((n.x+1.0)*127.5), ng=(float)((n.y+1.0)*127.5), nb=(float)((n.z+1.0)*127.5);
            for(int y=y0i;y<=y1i;y++){
                float py=(float)y+0.5f;
                int row=y*R;
                for(int x=x0i;x<=x1i;x++){
                    float px=(float)x+0.5f;
                    float w0=((v1-v2)*(px-u2)+(u2-u1)*(py-v2))*invDen;
                    float w1=((v2-v0)*(px-u2)+(u0-u2)*(py-v2))*invDen;
                    float w2=1.0f-w0-w1;
                    if(w0>=-1e-6f && w1>=-1e-6f && w2>=-1e-6f){
                        float dep=1.0f/(w0/z0+w1/z1+w2/z2);
                        int idx=row+x;
                        if(dep<zbuf[idx]){
                            zbuf[idx]=dep;
                            float* q=nbuf+idx*3; q[0]=nr; q[1]=ng; q[2]=nb;
                        }
                    }
                }
            }
        }
    }
}

static inline double rectSum(const vector<double>& I,int stride,int x0,int y0,int x1,int y1){
    return I[(size_t)y1*stride+x1]-I[(size_t)y0*stride+x1]-I[(size_t)y1*stride+x0]+I[(size_t)y0*stride+x0];
}
static double ssimChannel(const float* X,int strideX,const float* Y,int strideY,const float* dX,const float* dY,int R,
                          vector<double>& IX,vector<double>& IY,vector<double>& IX2,vector<double>& IY2,vector<double>& IXY){
    int W=R+1; size_t SZ=(size_t)W*W;
    fill(IX.begin(),IX.begin()+SZ,0.0); fill(IY.begin(),IY.begin()+SZ,0.0); fill(IX2.begin(),IX2.begin()+SZ,0.0); fill(IY2.begin(),IY2.begin()+SZ,0.0); fill(IXY.begin(),IXY.begin()+SZ,0.0);
    for(int y=1;y<=R;y++){
        double sx=0,sy=0,sx2=0,sy2=0,sxy=0; int row=(y-1)*R;
        for(int x=1;x<=R;x++){
            int p=row+(x-1);
            double xv=X[(size_t)p*strideX], yv=Y[(size_t)p*strideY];
            sx+=xv; sy+=yv; sx2+=xv*xv; sy2+=yv*yv; sxy+=xv*yv;
            size_t id=(size_t)y*W+x, up=(size_t)(y-1)*W+x;
            IX[id]=IX[up]+sx; IY[id]=IY[up]+sy; IX2[id]=IX2[up]+sx2; IY2[id]=IY2[up]+sy2; IXY[id]=IXY[up]+sxy;
        }
    }
    const int rad=5, area=121; const double c1=6.5025, c2=58.5225;
    long long cnt=0; long double acc=0.0L;
    for(int y=rad;y<R-rad;y++){
        int row=y*R;
        for(int x=rad;x<R-rad;x++){
            int p=row+x;
            if(!(dX[p]<254.0f || dY[p]<254.0f)) continue;
            int x0=x-rad, y0=y-rad, x1=x+rad+1, y1=y+rad+1;
            double sx=rectSum(IX,W,x0,y0,x1,y1), sy=rectSum(IY,W,x0,y0,x1,y1);
            double sx2=rectSum(IX2,W,x0,y0,x1,y1), sy2=rectSum(IY2,W,x0,y0,x1,y1), sxy=rectSum(IXY,W,x0,y0,x1,y1);
            double mux=sx/area, muy=sy/area;
            double vx=sx2/area-mux*mux, vy=sy2/area-muy*muy, cov=sxy/area-mux*muy;
            if(vx<0 && vx>-1e-6) vx=0; if(vy<0 && vy>-1e-6) vy=0;
            double denom=(mux*mux+muy*muy+c1)*(vx+vy+c2);
            double val=(denom!=0.0)?((2*mux*muy+c1)*(2*cov+c2)/denom):1.0;
            acc += val; ++cnt;
        }
    }
    return cnt? (double)(acc/cnt) : 1.0;
}
static double renderSSIM(const RenderMaps& orig,const RenderMaps& cand){
    int R=orig.R, PIX=R*R, W=R+1; size_t SZ=(size_t)W*W;
    vector<double> IX(SZ),IY(SZ),IX2(SZ),IY2(SZ),IXY(SZ);
    double total=0.0;
    for(int view=0;view<VIEWS;view++){
        const float* od=orig.depth.data()+(size_t)view*PIX;
        const float* cd=cand.depth.data()+(size_t)view*PIX;
        double ns=0.0;
        for(int ch=0;ch<3;ch++){
            const float* on=orig.norm.data()+(size_t)view*PIX*3+ch;
            const float* cn=cand.norm.data()+(size_t)view*PIX*3+ch;
            ns += ssimChannel(on,3,cn,3,od,cd,R,IX,IY,IX2,IY2,IXY);
        }
        ns/=3.0;
        double ds = ssimChannel(od,1,cd,1,od,cd,R,IX,IY,IX2,IY2,IXY);
        total += 0.5*(ns+ds);
    }
    return total/6.0;
}

static void simplifyCore(){
    if(N<=4 || M==0){ finalSnap.V=P; finalSnap.F.clear(); for(auto& f:F) finalSnap.F.push_back({f.v[0],f.v[1],f.v[2]}); haveFinalSnap=true; return; }
    double xmin=P[0].x,xmax=P[0].x,ymin=P[0].y,ymax=P[0].y,zmin=P[0].z,zmax=P[0].z;
    for(int i=1;i<N;i++){ xmin=min(xmin,P[i].x); xmax=max(xmax,P[i].x); ymin=min(ymin,P[i].y); ymax=max(ymax,P[i].y); zmin=min(zmin,P[i].z); zmax=max(zmax,P[i].z); }
    double diag=sqrt((xmax-xmin)*(xmax-xmin)+(ymax-ymin)*(ymax-ymin)+(zmax-zmin)*(zmax-zmin));
    hausTau=0.05*diag*0.992; hausTau2=hausTau*hausTau;
    if(N<8000) visualPixelTol=28.0; else if(N<60000) visualPixelTol=23.0; else if(N<200000) visualPixelTol=19.0; else visualPixelTol=16.0;
    visualPixelTol2=visualPixelTol*visualPixelTol;

    inc.assign(N,{}); for(int i=0;i<N;i++) inc[i].reserve(8);
    activeFaces=M; activeVertices=N;
    for(int i=0;i<M;i++){
        F[i].active=1; F[i].n=computeFaceNormal(F[i]);
        inc[F[i].v[0]].push_back(i); inc[F[i].v[1]].push_back(i); inc[F[i].v[2]].push_back(i);
    }
    Qv.assign(N,Quadric());
    for(int i=0;i<M;i++){
        Face& f=F[i]; Vec3 p0=P[f.v[0]],p1=P[f.v[1]],p2=P[f.v[2]];
        Vec3 cr=crossv(p1-p0,p2-p0); double area2=normv(cr); if(area2<=1e-300) continue;
        Vec3 n=cr/area2; double d=-dotv(n,p0);
        // Mostly screen/area-weighted, with a tiny uniform component to avoid erasing small visible creases too early.
        double w=max(1e-8,0.5*area2) + 2e-7;
        Quadric q=Quadric::plane(n.x,n.y,n.z,d,w);
        Qv[f.v[0]].add(q); Qv[f.v[1]].add(q); Qv[f.v[2]].add(q);
    }
    bbMin.resize(N); bbMax.resize(N);
    for(int i=0;i<N;i++){ bbMin[i]={(float)P[i].x,(float)P[i].y,(float)P[i].z}; bbMax[i]=bbMin[i]; }
    alive.assign(N,1); ver.assign(N,0); markA.assign(N,0); markB.assign(N,0);

    vector<uint64_t> edges; edges.reserve((size_t)M*3);
    auto addEdge=[&](int a,int b){ if(a>b) swap(a,b); edges.push_back(((uint64_t)(uint32_t)a<<32)|(uint32_t)b); };
    for(int i=0;i<M;i++){ addEdge(F[i].v[0],F[i].v[1]); addEdge(F[i].v[1],F[i].v[2]); addEdge(F[i].v[2],F[i].v[0]); }
    sort(edges.begin(),edges.end()); edges.erase(unique(edges.begin(),edges.end()),edges.end());
    vector<Candidate> heapVec; heapVec.reserve(edges.size());
    for(uint64_t key: edges){ int u=(int)(key>>32), v=(int)(key&0xffffffffu); EvalCand ec=computeCandidate(u,v); if(ec.ok) heapVec.push_back({ec.score,u,v,0,0}); }
    edges.clear(); edges.shrink_to_fit();
    priority_queue<Candidate> pq(less<Candidate>(),move(heapVec));

    vector<double> ratios={0.50,0.35,0.25,0.22,0.20,0.19,0.18,0.17,0.16,0.15,0.14,0.13,0.12,0.115,0.110,0.105,0.100,0.095,0.090,0.085,0.080,0.075,0.070,0.065,0.060,0.055};
    vector<int> cps;
    int last=-1;
    for(double r: ratios){ int c=max(4,(int)ceil(N*r)); if(c<N && c!=last){ cps.push_back(c); last=c; } }
    sort(cps.begin(),cps.end(),greater<int>());
    int cpIdx=0; int minTarget=cps.empty()?max(4,N/10):cps.back();

    while(!pq.empty() && activeVertices>minTarget){
        Candidate c=pq.top(); pq.pop(); int u=c.u,v=c.v;
        if(u==v || !alive[u] || !alive[v]) continue; if(u>v) swap(u,v);
        if(ver[u]!=c.vu || ver[v]!=c.vv) continue;
        EvalCand ec=computeCandidate(u,v); if(!ec.ok) continue;
        if(!linkConditionOK(u,v)) continue;
        int keep=u, rem=v; if(inc[v].size()>inc[u].size()){ keep=v; rem=u; }
        // Local normal-change limit. 0.12 allows strong simplification while rejecting obvious shading disasters.
        if(!flipOK(keep,rem,ec.pos,-0.02)) continue;
        doCollapse(keep,rem,ec.pos,pq);
        while(cpIdx<(int)cps.size() && activeVertices<=cps[cpIdx]){
            snapshots.push_back(makeSnapshot((double)activeVertices/(double)N));
            ++cpIdx;
        }
    }
    // Always have at least one compact candidate.
    if(snapshots.empty() || snapshots.back().nv!=activeVertices) snapshots.push_back(makeSnapshot((double)activeVertices/(double)N));
}

static void chooseByInternalSSIM(const RenderMaps& origMaps){
    if(snapshots.empty()){ finalSnap=makeSnapshot((double)activeVertices/(double)N); haveFinalSnap=true; return; }
    int R=origMaps.R;
    double guard = (R>=1000 ? 0.904 : (R>=700 ? 0.910 : 0.916));
    RenderMaps cand;
    int bestIdx=-1; double bestScore=-1.0;
    auto evalOne = [&](int idx)->double{
        if(idx<0 || idx>=(int)snapshots.size()) return -1.0;
        if(snapshots[idx].ssim >= -0.5) return snapshots[idx].ssim;
        renderMesh(snapshots[idx].V,snapshots[idx].F,cand,R);
        double sc=renderSSIM(origMaps,cand);
        snapshots[idx].ssim=sc;
        if(sc>bestScore){ bestScore=sc; bestIdx=idx; }
        return sc;
    };

    // Coarse-to-fine selection. Sparse probes keep the large tests fast; after a
    // passing bracket is found, all saved ratios inside that bracket are tested.
    vector<double> desired={0.055,0.065,0.075,0.085,0.095,0.105,0.115,0.130,0.160,0.200,0.250,0.350,0.500};
    vector<int> coarse;
    vector<unsigned char> used(snapshots.size(),0);
    for(double d: desired){
        int bi=-1; double bd=1e100;
        for(int i=0;i<(int)snapshots.size();i++){
            double diff=fabs(snapshots[i].ratio-d);
            if(diff<bd){ bd=diff; bi=i; }
        }
        if(bi>=0 && !used[bi]){ used[bi]=1; coarse.push_back(bi); }
    }
    sort(coarse.begin(),coarse.end(),[&](int a,int b){ return snapshots[a].ratio < snapshots[b].ratio; });

    int prevFail=-1;
    for(int idx: coarse){
        double sc=evalOne(idx);
        if(sc>=guard){
            int chosen=idx;
            if(prevFail!=-1){
                // snapshots are stored conservative -> aggressive; walk from just above the last fail
                // toward this pass, selecting the first ratio that satisfies the guard.
                for(int j=prevFail-1;j>=idx;j--){
                    double sj=evalOne(j);
                    if(sj>=guard){ chosen=j; break; }
                }
            }
            finalSnap=std::move(snapshots[chosen]); haveFinalSnap=true; return;
        }
        prevFail=idx;
    }
    if(bestIdx<0){
        for(int i=0;i<(int)snapshots.size();i++) evalOne(i);
    }
    if(bestIdx>=0){ finalSnap=std::move(snapshots[bestIdx]); haveFinalSnap=true; return; }
    finalSnap=makeSnapshot((double)activeVertices/(double)N); haveFinalSnap=true;
}

static void load(){
    FastScanner fs; fs.readInt(N); fs.readInt(M); P.resize(N);
    for(int i=0;i<N;i++){ fs.readPrefix('v'); fs.readDouble(P[i].x); fs.readDouble(P[i].y); fs.readDouble(P[i].z); }
    F.resize(M);
    for(int i=0;i<M;i++){ fs.readPrefix('f'); int a,b,c; fs.readInt(a); fs.readInt(b); fs.readInt(c); F[i].v[0]=a-1; F[i].v[1]=b-1; F[i].v[2]=c-1; F[i].active=1; }
}
static void save(){
    const vector<Vec3>& Vout = haveFinalSnap ? finalSnap.V : P;
    const vector<Tri>& Fout = finalSnap.F;
    string out; out.reserve((size_t)Vout.size()*40 + (size_t)Fout.size()*24 + 64); char line[128];
    int len=snprintf(line,sizeof(line),"%d %d\n",(int)Vout.size(),(int)Fout.size()); out.append(line,len);
    for(const Vec3& p: Vout){ len=snprintf(line,sizeof(line),"v %.10g %.10g %.10g\n",p.x,p.y,p.z); out.append(line,len); }
    for(const Tri& t: Fout){ len=snprintf(line,sizeof(line),"f %d %d %d\n",t.a+1,t.b+1,t.c+1); out.append(line,len); }
    fwrite(out.data(),1,out.size(),stdout);
}
int main(){
    load();
    // Render the original once. Exact 1024 for smaller tests; 512 for large cases to protect runtime.
    evalR = 1024; evalPix=evalR*evalR;
    vector<Tri> origTris; origTris.reserve(M);
    for(int i=0;i<M;i++) origTris.push_back({F[i].v[0],F[i].v[1],F[i].v[2]});
    RenderMaps origMaps; renderMesh(P,origTris,origMaps,evalR);
    origTris.clear(); origTris.shrink_to_fit();

    simplifyCore();

    // Free heavy QEM state before rendering candidate checkpoints.
    inc.clear(); inc.shrink_to_fit(); Qv.clear(); Qv.shrink_to_fit(); bbMin.clear(); bbMin.shrink_to_fit(); bbMax.clear(); bbMax.shrink_to_fit();
    alive.clear(); alive.shrink_to_fit(); ver.clear(); ver.shrink_to_fit(); markA.clear(); markA.shrink_to_fit(); markB.clear(); markB.shrink_to_fit(); tmpId.clear(); tmpId.shrink_to_fit();

    chooseByInternalSSIM(origMaps);
    save();
    return 0;
}
