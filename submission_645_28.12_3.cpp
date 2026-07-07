#include <bits/stdc++.h>
using namespace std;

struct Vec3 {
    double x, y, z;
    Vec3() : x(0), y(0), z(0) {}
    Vec3(double X, double Y, double Z) : x(X), y(Y), z(Z) {}
    Vec3 operator+(const Vec3& o) const { return Vec3(x+o.x, y+o.y, z+o.z); }
    Vec3 operator-(const Vec3& o) const { return Vec3(x-o.x, y-o.y, z-o.z); }
    Vec3 operator*(double s) const { return Vec3(x*s, y*s, z*s); }
    Vec3 operator/(double s) const { return Vec3(x/s, y/s, z/s); }
};
static inline double dotv(const Vec3& a, const Vec3& b){ return a.x*b.x + a.y*b.y + a.z*b.z; }
static inline Vec3 crossv(const Vec3& a, const Vec3& b){ return Vec3(a.y*b.z-a.z*b.y, a.z*b.x-a.x*b.z, a.x*b.y-a.y*b.x); }
static inline double norm2(const Vec3& a){ return dotv(a,a); }
static inline double normv(const Vec3& a){ return sqrt(norm2(a)); }
static inline double distv(const Vec3& a, const Vec3& b){ return normv(a-b); }

struct Quadric {
    // symmetric 4x4: xx,xy,xz,xw, yy,yz,yw, zz,zw, ww
    double q[10];
    Quadric(){ memset(q,0,sizeof(q)); }
    void addPlane(double a,double b,double c,double d,double w){
        q[0]+=w*a*a; q[1]+=w*a*b; q[2]+=w*a*c; q[3]+=w*a*d;
        q[4]+=w*b*b; q[5]+=w*b*c; q[6]+=w*b*d;
        q[7]+=w*c*c; q[8]+=w*c*d; q[9]+=w*d*d;
    }
    void add(const Quadric& o){ for(int i=0;i<10;i++) q[i]+=o.q[i]; }
    double eval(const Vec3& p) const {
        double x=p.x,y=p.y,z=p.z;
        return q[0]*x*x + 2*q[1]*x*y + 2*q[2]*x*z + 2*q[3]*x
             + q[4]*y*y + 2*q[5]*y*z + 2*q[6]*y
             + q[7]*z*z + 2*q[8]*z + q[9];
    }
};

struct Vertex {
    Vec3 p;
    Quadric Q;
    vector<int> faces;
    double rad = 0.0;      // exact sufficient cluster radius to current representative p
    double feature = 0.0;  // original normal variation around this vertex
    int version = 0;
    bool alive = true;
};
struct Face {
    int v[3];
    Vec3 n;
    double area2 = 0.0;
    bool alive = true;
};

struct Candidate {
    double cost;
    int u, v;
    int vu, vv;
    bool operator<(const Candidate& o) const { return cost > o.cost; } // min-heap
};

static vector<Vertex> Vtx;
static vector<Face> Faces;
static int aliveV = 0;
static Vec3 bboxMin, bboxMax;
static double bboxDiag = 0, hausTol = 0, hausTol2 = 0;
static double collapseTol = 0;
static int targetV = 4;
static chrono::steady_clock::time_point startTime;
static double timeLimitSec = 19.2;

static vector<char> slurp_stdin() {
    vector<char> buf;
    buf.reserve(1 << 27);
    char chunk[1 << 16];
    size_t n;
    while ((n=fread(chunk,1,sizeof(chunk),stdin))>0) buf.insert(buf.end(), chunk, chunk+n);
    buf.push_back('\0');
    return buf;
}
static inline void skip_ws(char*& p){ while(*p==' '||*p=='\n'||*p=='\r'||*p=='\t') ++p; }

static void read_input(){
    vector<char> buf=slurp_stdin();
    char* p=buf.data();
    long nv=strtol(p,&p,10), nf=strtol(p,&p,10);
    Vtx.resize(nv); Faces.resize(nf); aliveV=(int)nv;
    Vec3 mn(1e100,1e100,1e100), mx(-1e100,-1e100,-1e100);
    for(long i=0;i<nv;i++){
        skip_ws(p); if(*p=='v') ++p;
        double x=strtod(p,&p), y=strtod(p,&p), z=strtod(p,&p);
        Vtx[i].p=Vec3(x,y,z);
        mn.x=min(mn.x,x); mn.y=min(mn.y,y); mn.z=min(mn.z,z);
        mx.x=max(mx.x,x); mx.y=max(mx.y,y); mx.z=max(mx.z,z);
    }
    for(long i=0;i<nf;i++){
        skip_ws(p); if(*p=='f') ++p;
        int a=(int)strtol(p,&p,10)-1, b=(int)strtol(p,&p,10)-1, c=(int)strtol(p,&p,10)-1;
        Faces[i].v[0]=a; Faces[i].v[1]=b; Faces[i].v[2]=c;
        Vtx[a].faces.push_back((int)i);
        Vtx[b].faces.push_back((int)i);
        Vtx[c].faces.push_back((int)i);
    }
    bboxMin = mn; bboxMax = mx;
    bboxDiag = normv(mx-mn);
    hausTol = 0.05 * bboxDiag;
    hausTol2 = hausTol * hausTol;
    collapseTol = hausTol * 0.985;
}

static bool recompute_face(int fid){
    Face &f=Faces[fid];
    Vec3 a=Vtx[f.v[0]].p, b=Vtx[f.v[1]].p, c=Vtx[f.v[2]].p;
    Vec3 cr=crossv(b-a, c-a);
    double ar=normv(cr);
    f.area2=ar;
    if(ar<=1e-18) { f.n=Vec3(0,0,0); return false; }
    f.n=cr/ar;
    return true;
}

static void init_quadrics_features(){
    int n=(int)Vtx.size();
    vector<Vec3> nsum(n, Vec3(0,0,0));
    for(int i=0;i<(int)Faces.size();i++){
        recompute_face(i);
        Face &f=Faces[i];
        Vec3 p0=Vtx[f.v[0]].p;
        double d=-dotv(f.n,p0);
        // sqrt-area weighting is intentionally mild; huge triangles should not dominate all detail.
        double w=max(1e-12, sqrt(max(0.0, f.area2)));
        for(int k=0;k<3;k++){
            Vtx[f.v[k]].Q.addPlane(f.n.x,f.n.y,f.n.z,d,w);
            nsum[f.v[k]] = nsum[f.v[k]] + f.n * f.area2;
        }
    }
    for(int i=0;i<n;i++){
        double l=normv(nsum[i]);
        Vec3 avg = (l>1e-30 ? nsum[i]/l : Vec3(0,0,0));
        double mind=1.0;
        for(int fid: Vtx[i].faces) if(Faces[fid].alive){
            mind=min(mind, dotv(avg, Faces[fid].n));
        }
        Vtx[i].feature = max(0.0, 1.0 - mind); // 0 smooth/planar, high on creases/corners/noise
    }
}

static inline bool face_contains(const Face& f, int a){ return f.v[0]==a || f.v[1]==a || f.v[2]==a; }
static inline bool face_contains2(const Face& f, int a, int b){ return face_contains(f,a) && face_contains(f,b); }
static inline int other_vertex_in_face(const Face& f, int a, int b){
    for(int k=0;k<3;k++) if(f.v[k]!=a && f.v[k]!=b) return f.v[k];
    return -1;
}

static void cleanup_faces_list(int v){
    if(v<0 || v>=(int)Vtx.size() || !Vtx[v].alive) return;
    vector<int>& L=Vtx[v].faces;
    int wr=0;
    for(int id: L){
        if(id>=0 && id<(int)Faces.size() && Faces[id].alive && face_contains(Faces[id], v)) L[wr++]=id;
    }
    L.resize(wr);
    sort(L.begin(), L.end());
    L.erase(unique(L.begin(), L.end()), L.end());
}

static vector<int> alive_faces_of(int v){
    vector<int> res;
    if(v<0 || v>=(int)Vtx.size() || !Vtx[v].alive) return res;
    auto &L=Vtx[v].faces;
    res.reserve(L.size());
    for(int id: L) if(id>=0 && id<(int)Faces.size() && Faces[id].alive && face_contains(Faces[id], v)) res.push_back(id);
    if((int)L.size() > (int)res.size()*3 + 64) cleanup_faces_list(v);
    return res;
}

static void unique_small(vector<int>& a){ sort(a.begin(), a.end()); a.erase(unique(a.begin(), a.end()), a.end()); }

struct CollapseChoice { int keep=-1, kill=-1; double cost=0; double newRad=0; Vec3 pos; };


static bool solve_qem_position(const Quadric& Q, Vec3& out){
    double a00=Q.q[0], a01=Q.q[1], a02=Q.q[2];
    double a11=Q.q[4], a12=Q.q[5], a22=Q.q[7];
    double b0=-Q.q[3], b1=-Q.q[6], b2=-Q.q[8];
    double det = a00*(a11*a22-a12*a12) - a01*(a01*a22-a12*a02) + a02*(a01*a12-a11*a02);
    if(fabs(det) < 1e-18) return false;
    double id=1.0/det;
    double c00=(a11*a22-a12*a12)*id;
    double c01=(a02*a12-a01*a22)*id;
    double c02=(a01*a12-a02*a11)*id;
    double c11=(a00*a22-a02*a02)*id;
    double c12=(a01*a02-a00*a12)*id;
    double c22=(a00*a11-a01*a01)*id;
    out=Vec3(c00*b0+c01*b1+c02*b2,
             c01*b0+c11*b1+c12*b2,
             c02*b0+c12*b1+c22*b2);
    return isfinite(out.x) && isfinite(out.y) && isfinite(out.z);
}

static double vertex_limit(int a, int b){
    // Preserve very sharp/noisy vertices a little more tightly, but still use most of the official 5% budget.
    double f=max(Vtx[a].feature, Vtx[b].feature);
    double mul=0.985;
    if(f>0.35) mul=0.90;
    else if(f>0.18) mul=0.94;
    else if(f>0.08) mul=0.965;
    return hausTol * mul;
}

static bool choose_collapse(int a, int b, CollapseChoice &ch){
    ch=CollapseChoice();
    if(a==b || a<0 || b<0 || a>=(int)Vtx.size() || b>=(int)Vtx.size()) return false;
    if(!Vtx[a].alive || !Vtx[b].alive) return false;
    double edgeLen=distv(Vtx[a].p, Vtx[b].p);
    double lim=vertex_limit(a,b);
    Quadric Q=Vtx[a].Q; Q.add(Vtx[b].Q);

    vector<Vec3> cand;
    cand.reserve(5);
    cand.push_back(Vtx[a].p);
    cand.push_back(Vtx[b].p);
    cand.push_back((Vtx[a].p + Vtx[b].p) * 0.5);
    double wa=max(1e-9, Vtx[b].rad + edgeLen*0.5), wb=max(1e-9, Vtx[a].rad + edgeLen*0.5);
    cand.push_back((Vtx[a].p*wa + Vtx[b].p*wb) / (wa+wb));
    Vec3 opt;
    if(solve_qem_position(Q,opt)) cand.push_back(opt);

    double best=1e300, bestRad=0.0;
    Vec3 bestPos;
    bool any=false;
    for(const Vec3& p: cand){
        if(!isfinite(p.x)||!isfinite(p.y)||!isfinite(p.z)) continue;
        double r=max(Vtx[a].rad + distv(p,Vtx[a].p), Vtx[b].rad + distv(p,Vtx[b].p));
        if(r > lim) continue;
        double c=Q.eval(p) + 1e-9*edgeLen*edgeLen + 1e-7*r*r;
        // Slightly prefer not moving very sharp representatives too much.
        double feat=max(Vtx[a].feature,Vtx[b].feature);
        if(feat>0.12) c += 1e-5*(distv(p,Vtx[a].p)+distv(p,Vtx[b].p));
        if(c<best){ best=c; bestRad=r; bestPos=p; any=true; }
    }
    if(!any) return false;
    // Keep the higher-valence endpoint to reduce index churn; if the best point is exactly an endpoint,
    // keep that endpoint to avoid unnecessary movement of the other one-ring.
    double da=distv(bestPos,Vtx[a].p), db=distv(bestPos,Vtx[b].p);
    if(da < 1e-14) { ch.keep=a; ch.kill=b; }
    else if(db < 1e-14) { ch.keep=b; ch.kill=a; }
    else if(Vtx[a].faces.size() >= Vtx[b].faces.size()) { ch.keep=a; ch.kill=b; }
    else { ch.keep=b; ch.kill=a; }
    ch.cost=best;
    ch.newRad=bestRad;
    ch.pos=bestPos;
    return true;
}

static bool has_duplicate_face_after_replace(int fid, int keep, int kill){
    int tri[3];
    for(int k=0;k<3;k++) tri[k]=(Faces[fid].v[k]==kill?keep:Faces[fid].v[k]);
    if(tri[0]==tri[1] || tri[0]==tri[2] || tri[1]==tri[2]) return true;
    int keyv=tri[0];
    if(Vtx[tri[1]].faces.size() < Vtx[keyv].faces.size()) keyv=tri[1];
    if(Vtx[tri[2]].faces.size() < Vtx[keyv].faces.size()) keyv=tri[2];
    int s0=tri[0], s1=tri[1], s2=tri[2];
    if(s0>s1) swap(s0,s1); if(s1>s2) swap(s1,s2); if(s0>s1) swap(s0,s1);
    for(int ofid: Vtx[keyv].faces){
        if(ofid==fid || ofid<0 || ofid>=(int)Faces.size() || !Faces[ofid].alive) continue;
        // Faces incident to the collapsing edge are removed, so they cannot duplicate the new face.
        if(face_contains2(Faces[ofid], keep, kill)) continue;
        int t[3]={Faces[ofid].v[0],Faces[ofid].v[1],Faces[ofid].v[2]};
        if(t[0]>t[1]) swap(t[0],t[1]); if(t[1]>t[2]) swap(t[1],t[2]); if(t[0]>t[1]) swap(t[0],t[1]);
        if(t[0]==s0 && t[1]==s1 && t[2]==s2) return true;
    }
    return false;
}

static bool check_link_and_geometry(int keep, int kill, const Vec3& newPos,
                                     vector<int>& removeFaces, vector<int>& changeFaces,
                                     vector<int>& affectedFaces){
    removeFaces.clear(); changeFaces.clear(); affectedFaces.clear();
    if(!Vtx[keep].alive || !Vtx[kill].alive) return false;
    vector<int> fkeep=alive_faces_of(keep);
    vector<int> fkill=alive_faces_of(kill);
    vector<int> nKeep, nKill, opposites;
    nKeep.reserve(fkeep.size()*2); nKill.reserve(fkill.size()*2);
    for(int fid: fkeep){
        Face &f=Faces[fid];
        for(int k=0;k<3;k++) if(f.v[k]!=keep) nKeep.push_back(f.v[k]);
    }
    for(int fid: fkill){
        Face &f=Faces[fid];
        bool both=face_contains(f, keep);
        if(both){ removeFaces.push_back(fid); opposites.push_back(other_vertex_in_face(f,keep,kill)); }
        else changeFaces.push_back(fid);
        for(int k=0;k<3;k++) if(f.v[k]!=kill) nKill.push_back(f.v[k]);
    }
    unique_small(removeFaces); unique_small(changeFaces); unique_small(nKeep); unique_small(nKill); unique_small(opposites);
    if(removeFaces.size()!=2 || opposites.size()!=2) return false;
    vector<int> common;
    size_t i=0,j=0;
    while(i<nKeep.size() && j<nKill.size()){
        if(nKeep[i]==nKill[j]){ if(nKeep[i]!=keep && nKeep[i]!=kill) common.push_back(nKeep[i]); ++i; ++j; }
        else if(nKeep[i]<nKill[j]) ++i; else ++j;
    }
    unique_small(common);
    if(common.size()!=2) return false;
    if(!(common[0]==opposites[0] && common[1]==opposites[1])) return false;

    for(int fid: changeFaces) if(has_duplicate_face_after_replace(fid, keep, kill)) return false;

    // All faces incident to either endpoint are geometrically affected because the survivor may move.
    affectedFaces = fkeep;
    affectedFaces.insert(affectedFaces.end(), changeFaces.begin(), changeFaces.end());
    unique_small(affectedFaces);
    vector<int> remSorted=removeFaces;
    unique_small(remSorted);
    vector<int> tmp; tmp.reserve(affectedFaces.size());
    for(int fid: affectedFaces){
        if(!binary_search(remSorted.begin(), remSorted.end(), fid)) tmp.push_back(fid);
    }
    affectedFaces.swap(tmp);

    const double minArea2 = max(1e-20, bboxDiag*bboxDiag*1e-20);
    double cosHard = -0.05;
    double feat=max(Vtx[keep].feature, Vtx[kill].feature);
    if(feat>0.25) cosHard=0.18;
    else if(feat>0.10) cosHard=0.02;
    for(int fid: affectedFaces){
        Face &old=Faces[fid];
        int ids[3]={old.v[0],old.v[1],old.v[2]};
        for(int k=0;k<3;k++) if(ids[k]==kill) ids[k]=keep;
        if(ids[0]==ids[1] || ids[0]==ids[2] || ids[1]==ids[2]) return false;
        Vec3 p[3];
        for(int k=0;k<3;k++) p[k]=(ids[k]==keep?newPos:Vtx[ids[k]].p);
        Vec3 cr=crossv(p[1]-p[0], p[2]-p[0]);
        double ar=normv(cr);
        if(ar <= minArea2) return false;
        Vec3 nn=cr/ar;
        if(dotv(nn, old.n) < cosHard) return false;
    }
    return true;
}

static bool collapse_edge(int keep, int kill, const Vec3& newPos, double newRad, priority_queue<Candidate>& pq){
    vector<int> rem, chg, affected;
    if(!check_link_and_geometry(keep, kill, newPos, rem, chg, affected)) return false;

    Quadric newQ=Vtx[keep].Q; newQ.add(Vtx[kill].Q);
    for(int fid: rem) Faces[fid].alive=false;

    // Move survivor before recomputing all affected normals.
    Vtx[keep].p = newPos;
    for(int fid: chg){
        Face &f=Faces[fid];
        for(int k=0;k<3;k++) if(f.v[k]==kill) f.v[k]=keep;
        Vtx[keep].faces.push_back(fid);
    }
    for(int fid: affected) if(Faces[fid].alive) recompute_face(fid);

    Vtx[keep].Q = newQ;
    Vtx[keep].rad = newRad;
    Vtx[keep].feature = max(Vtx[keep].feature, Vtx[kill].feature);
    Vtx[keep].version++;
    Vtx[kill].alive=false;
    Vtx[kill].version++;
    aliveV--;

    if(Vtx[keep].faces.size() > 256) cleanup_faces_list(keep);
    vector<int> lf=alive_faces_of(keep);
    for(int fid: lf){
        Face &f=Faces[fid];
        for(int k=0;k<3;k++){
            int nb=f.v[k]; if(nb==keep || !Vtx[nb].alive) continue;
            CollapseChoice cc;
            if(choose_collapse(keep, nb, cc)){
                pq.push({cc.cost, keep, nb, Vtx[keep].version, Vtx[nb].version});
            }
        }
    }
    return true;
}

static void add_initial_candidates(priority_queue<Candidate>& pq){
    vector<unsigned long long> edges;
    edges.reserve((size_t)Faces.size()*3);
    for(int i=0;i<(int)Faces.size();i++){
        int a=Faces[i].v[0], b=Faces[i].v[1], c=Faces[i].v[2];
        int e[3][2]={{a,b},{b,c},{c,a}};
        for(auto &x: e){
            int u=x[0],v=x[1]; if(u>v) swap(u,v);
            edges.push_back((unsigned long long)(unsigned int)u<<32 | (unsigned int)v);
        }
    }
    sort(edges.begin(), edges.end());
    edges.erase(unique(edges.begin(), edges.end()), edges.end());
    for(auto key: edges){
        int u=(int)(key>>32), v=(int)(key & 0xffffffffu);
        CollapseChoice cc;
        if(choose_collapse(u,v,cc)) pq.push({cc.cost,u,v,Vtx[u].version,Vtx[v].version});
    }
}

static void compute_target(){
    int n=(int)Vtx.size();
    double ratio;
    if(n <= 20) ratio = 0.80;            // sample-like: only remove truly redundant vertices
    else if(n <= 1000) ratio = 0.30;
    else if(n <= 6000) ratio = 0.250;
    else if(n <= 30000) ratio = 0.080;
    else if(n <= 70000) ratio = 0.050;
    else if(n <= 500000) ratio = 0.025;
    else ratio = 0.015;

    if(n > 20){
        // Estimate how much sharp detail/noise the mesh has. If many vertices are creases/corners,
        // keep more vertices to avoid SSIM failures on flat-normal images.
        int sampleCnt=0, sharpCnt=0;
        for(int i=0;i<n;i+=(n>20000? n/20000 : 1)){
            sampleCnt++;
            if(Vtx[i].feature > 0.12) sharpCnt++;
        }
        double sharpRatio = sampleCnt? (double)sharpCnt/sampleCnt : 0.0;
        if(sharpRatio > 0.45) ratio *= 1.15;
        else if(sharpRatio > 0.25) ratio *= 1.08;
        else if(sharpRatio > 0.12) ratio *= 1.04;
        ratio = min(ratio, 0.92);
    }

    // If the mesh is already quite small compared with the image resolution, avoid gambling too hard.
    targetV = max(4, (int)ceil(n * ratio));
    targetV = min(targetV, n);
}


static long long cell_key3(long long ix, long long iy, long long iz){
    // Mix three signed 21-bit-ish coordinates into one hashable key.
    const long long B = 1048576LL;
    return (ix+B)*4398046511104LL ^ (iy+B)*2097152LL ^ (iz+B);
}

static bool try_box_remesh(){
    int n=(int)Vtx.size();
    if(n < 100) return false;
    Vec3 L=bboxMax-bboxMin;
    if(L.x<=1e-12 || L.y<=1e-12 || L.z<=1e-12) return false;
    double eps=max(1e-9, bboxDiag*1e-7);
    for(const auto& vt: Vtx){
        const Vec3& p=vt.p;
        bool on = fabs(p.x-bboxMin.x)<=eps || fabs(p.x-bboxMax.x)<=eps ||
                  fabs(p.y-bboxMin.y)<=eps || fabs(p.y-bboxMax.y)<=eps ||
                  fabs(p.z-bboxMin.z)<=eps || fabs(p.z-bboxMax.z)<=eps;
        if(!on) return false;
    }
    double step=max(hausTol*1.25, bboxDiag*1e-6);
    int nx=max(1,(int)ceil(L.x/step));
    int ny=max(1,(int)ceil(L.y/step));
    int nz=max(1,(int)ceil(L.z/step));
    auto surfCount = 1LL*(nx+1)*(ny+1)*(nz+1) - 1LL*max(0,nx-1)*max(0,ny-1)*max(0,nz-1);
    if(surfCount <= 0 || surfCount > n) return false;

    double orient[3][2]; memset(orient,0,sizeof(orient));
    for(const Face& ff: Faces){
        int a=ff.v[0], b=ff.v[1], c=ff.v[2];
        Vec3 cr=crossv(Vtx[b].p - Vtx[a].p, Vtx[c].p - Vtx[a].p);
        double ar=normv(cr); if(ar<=1e-30) continue;
        Vec3 nn=cr/ar;
        Vec3 cen=(Vtx[a].p+Vtx[b].p+Vtx[c].p)/3.0;
        if(fabs(cen.x-bboxMin.x)<=eps*10) orient[0][0]+=nn.x;
        if(fabs(cen.x-bboxMax.x)<=eps*10) orient[0][1]+=nn.x;
        if(fabs(cen.y-bboxMin.y)<=eps*10) orient[1][0]+=nn.y;
        if(fabs(cen.y-bboxMax.y)<=eps*10) orient[1][1]+=nn.y;
        if(fabs(cen.z-bboxMin.z)<=eps*10) orient[2][0]+=nn.z;
        if(fabs(cen.z-bboxMax.z)<=eps*10) orient[2][1]+=nn.z;
    }
    int desiredSign[3][2];
    for(int ax=0; ax<3; ++ax){ desiredSign[ax][0] = (orient[ax][0] >= 0 ? +1 : -1); desiredSign[ax][1] = (orient[ax][1] >= 0 ? +1 : -1); }
    // If a side had no reliable vote, use outward orientation as fallback.
    for(int ax=0; ax<3; ++ax){ if(fabs(orient[ax][0])<1e-12) desiredSign[ax][0]=-1; if(fabs(orient[ax][1])<1e-12) desiredSign[ax][1]=+1; }

    vector<Vec3> newV;
    vector<array<int,3>> newF;
    newV.reserve((size_t)surfCount);
    unordered_map<long long,int> id;
    id.reserve((size_t)surfCount*2+16);
    auto coord=[&](int ix,int iy,int iz){
        return Vec3(bboxMin.x + L.x * ((double)ix/nx),
                    bboxMin.y + L.y * ((double)iy/ny),
                    bboxMin.z + L.z * ((double)iz/nz));
    };
    auto get_id=[&](int ix,int iy,int iz)->int{
        long long key=cell_key3(ix,iy,iz);
        auto it=id.find(key);
        if(it!=id.end()) return it->second;
        int nid=(int)newV.size(); id[key]=nid; newV.push_back(coord(ix,iy,iz)); return nid;
    };
    int dims[3]={nx,ny,nz};
    auto add_face=[&](int axis, int fixedSide, int nsign){
        int a,b;
        if(axis==0){ if(nsign>0){a=1;b=2;} else {a=2;b=1;} }
        else if(axis==1){ if(nsign>0){a=2;b=0;} else {a=0;b=2;} }
        else { if(nsign>0){a=0;b=1;} else {a=1;b=0;} }
        int fixed = (fixedSide>0 ? dims[axis] : 0);
        for(int ia=0; ia<dims[a]; ++ia){
            for(int ib=0; ib<dims[b]; ++ib){
                int idx00[3], idx10[3], idx01[3], idx11[3];
                idx00[axis]=idx10[axis]=idx01[axis]=idx11[axis]=fixed;
                idx00[a]=ia;   idx00[b]=ib;
                idx10[a]=ia+1; idx10[b]=ib;
                idx01[a]=ia;   idx01[b]=ib+1;
                idx11[a]=ia+1; idx11[b]=ib+1;
                int p00=get_id(idx00[0],idx00[1],idx00[2]);
                int p10=get_id(idx10[0],idx10[1],idx10[2]);
                int p01=get_id(idx01[0],idx01[1],idx01[2]);
                int p11=get_id(idx11[0],idx11[1],idx11[2]);
                newF.push_back({p00,p10,p11});
                newF.push_back({p00,p11,p01});
            }
        }
    };
    for(int ax=0; ax<3; ++ax){ add_face(ax,-1,desiredSign[ax][0]); add_face(ax,+1,desiredSign[ax][1]); }
    if((int)newV.size() > n) return false;

    // Original vertices must be covered by lattice vertices.
    double tol2=hausTol*hausTol*0.998*0.998;
    for(const auto& vt: Vtx){
        const Vec3& p=vt.p;
        int ix=(int)llround((p.x-bboxMin.x)/L.x*nx); ix=max(0,min(nx,ix));
        int iy=(int)llround((p.y-bboxMin.y)/L.y*ny); iy=max(0,min(ny,iy));
        int iz=(int)llround((p.z-bboxMin.z)/L.z*nz); iz=max(0,min(nz,iz));
        Vec3 q=coord(ix,iy,iz);
        if(norm2(p-q) > tol2) return false;
    }

    // Lattice vertices must be close to some original vertex. This avoids failing sparse boxes.
    double cell=hausTol;
    if(cell<=0) return false;
    unordered_map<long long, vector<int>> grid;
    grid.reserve((size_t)n*2+16);
    auto cell_index=[&](const Vec3& p){
        long long ix=(long long)floor((p.x-bboxMin.x)/cell);
        long long iy=(long long)floor((p.y-bboxMin.y)/cell);
        long long iz=(long long)floor((p.z-bboxMin.z)/cell);
        return array<long long,3>{ix,iy,iz};
    };
    for(int i=0;i<n;i++){
        auto ci=cell_index(Vtx[i].p);
        grid[cell_key3(ci[0],ci[1],ci[2])].push_back(i);
    }
    double tol2b=hausTol*hausTol*1.002*1.002;
    for(const Vec3& q: newV){
        auto ci=cell_index(q);
        bool ok=false;
        for(long long dx=-1; dx<=1 && !ok; ++dx) for(long long dy=-1; dy<=1 && !ok; ++dy) for(long long dz=-1; dz<=1 && !ok; ++dz){
            auto it=grid.find(cell_key3(ci[0]+dx,ci[1]+dy,ci[2]+dz));
            if(it==grid.end()) continue;
            for(int oi: it->second){ if(norm2(q-Vtx[oi].p) <= tol2b){ ok=true; break; } }
        }
        if(!ok) return false;
    }

    Vtx.clear(); Faces.clear();
    Vtx.resize(newV.size()); Faces.resize(newF.size());
    for(int i=0;i<(int)newV.size();++i){ Vtx[i].p=newV[i]; Vtx[i].alive=true; Vtx[i].rad=0; }
    for(int i=0;i<(int)newF.size();++i){
        Faces[i].v[0]=newF[i][0]; Faces[i].v[1]=newF[i][1]; Faces[i].v[2]=newF[i][2]; Faces[i].alive=true;
        Vtx[Faces[i].v[0]].faces.push_back(i);
        Vtx[Faces[i].v[1]].faces.push_back(i);
        Vtx[Faces[i].v[2]].faces.push_back(i);
    }
    aliveV=(int)Vtx.size();
    for(int i=0;i<(int)Faces.size();++i) recompute_face(i);
    return true;
}

static bool time_exceeded(){
    using namespace chrono;
    double t=duration<double>(steady_clock::now()-startTime).count();
    return t > timeLimitSec;
}

static void simplify_mesh(){
    if(Vtx.empty()) return;
    if(try_box_remesh()) return;
    init_quadrics_features();
    compute_target();

    priority_queue<Candidate> pq;
    add_initial_candidates(pq);

    const double freeCost = 2.0e-9 * max(1.0, bboxDiag*bboxDiag);
    long long pops=0, successful=0;
    while(!pq.empty()){
        if((pops & 4095LL)==0 && time_exceeded()) break;
        Candidate cand=pq.top(); pq.pop(); pops++;
        if(aliveV <= targetV && cand.cost > freeCost) break;
        int u=cand.u, v=cand.v;
        if(u<0||v<0||u>=(int)Vtx.size()||v>=(int)Vtx.size()) continue;
        if(!Vtx[u].alive || !Vtx[v].alive) continue;
        if(Vtx[u].version!=cand.vu || Vtx[v].version!=cand.vv) continue;
        CollapseChoice cc;
        if(!choose_collapse(u,v,cc)) continue;
        if(aliveV <= targetV && cc.cost > freeCost) break;
        // Lazy priority: if the recomputed cost is much worse, reinsert once with fresh priority.
        if(cc.cost > cand.cost * 1.25 + 1e-14){
            pq.push({cc.cost, u, v, Vtx[u].version, Vtx[v].version});
            continue;
        }
        if(collapse_edge(cc.keep, cc.kill, cc.pos, cc.newRad, pq)) successful++;
    }
}

static void write_output(){
    int n=(int)Vtx.size(), m=(int)Faces.size();
    vector<int> used(n,0);
    int aliveF=0;
    for(int i=0;i<m;i++) if(Faces[i].alive){
        int a=Faces[i].v[0], b=Faces[i].v[1], c=Faces[i].v[2];
        if(a==b||a==c||b==c) { Faces[i].alive=false; continue; }
        Vec3 cr=crossv(Vtx[b].p - Vtx[a].p, Vtx[c].p - Vtx[a].p);
        if(norm2(cr) <= 1e-30) { Faces[i].alive=false; continue; }
        used[a]=used[b]=used[c]=1; aliveF++;
    }
    vector<int> id(n,-1);
    int outV=0;
    for(int i=0;i<n;i++) if(used[i] && Vtx[i].alive) id[i]=outV++;
    // Fallback should never trigger, but preserving a valid output is better than risking empty output.
    if(outV<4 || aliveF<4){
        // Output original mesh compactly.
        outV=n; aliveF=m;
        string out; out.reserve((size_t)n*40 + (size_t)m*24 + 32);
        char line[128];
        out.append(line, snprintf(line,sizeof(line), "%d %d\n", n, m));
        for(int i=0;i<n;i++) out.append(line, snprintf(line,sizeof(line), "v %.10g %.10g %.10g\n", Vtx[i].p.x,Vtx[i].p.y,Vtx[i].p.z));
        for(int i=0;i<m;i++) out.append(line, snprintf(line,sizeof(line), "f %d %d %d\n", Faces[i].v[0]+1,Faces[i].v[1]+1,Faces[i].v[2]+1));
        fwrite(out.data(),1,out.size(),stdout);
        return;
    }
    string out;
    out.reserve((size_t)outV*42 + (size_t)aliveF*26 + 64);
    char line[128];
    out.append(line, snprintf(line,sizeof(line), "%d %d\n", outV, aliveF));
    for(int i=0;i<n;i++) if(id[i]>=0){
        out.append(line, snprintf(line,sizeof(line), "v %.10g %.10g %.10g\n", Vtx[i].p.x,Vtx[i].p.y,Vtx[i].p.z));
    }
    for(int i=0;i<m;i++) if(Faces[i].alive){
        int a=id[Faces[i].v[0]], b=id[Faces[i].v[1]], c=id[Faces[i].v[2]];
        if(a>=0 && b>=0 && c>=0) out.append(line, snprintf(line,sizeof(line), "f %d %d %d\n", a+1,b+1,c+1));
    }
    fwrite(out.data(),1,out.size(),stdout);
}

int main(){
    startTime=chrono::steady_clock::now();
    read_input();
    simplify_mesh();
    write_output();
    return 0;
}
