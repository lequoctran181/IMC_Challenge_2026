#include <bits/stdc++.h>
using namespace std;

// Fixed-target topology-preserving QEM simplifier.
// Target retained vertex ratio is tuned for a 91+ score if all official cases pass FinalSSIM.

struct Vec3 {
    double x, y, z;
    Vec3() : x(0), y(0), z(0) {}
    Vec3(double X, double Y, double Z) : x(X), y(Y), z(Z) {}
};
static inline Vec3 operator+(const Vec3& a, const Vec3& b){ return Vec3(a.x+b.x,a.y+b.y,a.z+b.z); }
static inline Vec3 operator-(const Vec3& a, const Vec3& b){ return Vec3(a.x-b.x,a.y-b.y,a.z-b.z); }
static inline Vec3 operator*(const Vec3& a, double s){ return Vec3(a.x*s,a.y*s,a.z*s); }
static inline Vec3 operator/(const Vec3& a, double s){ return Vec3(a.x/s,a.y/s,a.z/s); }
static inline double dot3(const Vec3& a,const Vec3& b){ return a.x*b.x+a.y*b.y+a.z*b.z; }
static inline Vec3 cross3(const Vec3& a,const Vec3& b){ return Vec3(a.y*b.z-a.z*b.y, a.z*b.x-a.x*b.z, a.x*b.y-a.y*b.x); }
static inline double norm2(const Vec3& a){ return dot3(a,a); }
static inline double norm3(const Vec3& a){ return sqrt(norm2(a)); }

struct Face {
    int v[3];
    float nx, ny, nz;
    unsigned char alive;
};

struct Vertex {
    Vec3 p;
    double q[10];        // xx,xy,xz,xw, yy,yz,yw, zz,zw, ww
    double mn[3], mx[3]; // AABB of all original vertices represented by this vertex
    int head, tail;      // incident face linked-list nodes
    int ver;
    unsigned char alive;
};

struct MeshOut {
    vector<Vec3> V;
    vector<array<int,3>> F;
};

static vector<char> slurp_stdin() {
    vector<char> buf;
    buf.reserve(1 << 27);
    char chunk[1 << 16];
    size_t n;
    while ((n = fread(chunk, 1, sizeof(chunk), stdin)) > 0) buf.insert(buf.end(), chunk, chunk+n);
    buf.push_back('\0');
    return buf;
}
static inline void skip_ws(char*& p) {
    while (*p==' ' || *p=='\n' || *p=='\r' || *p=='\t') ++p;
}
static void load_obj(vector<Vec3>& V, vector<array<int,3>>& F) {
    vector<char> buf = slurp_stdin();
    char* p = buf.data();
    long nv = strtol(p, &p, 10);
    long nf = strtol(p, &p, 10);
    V.resize(nv);
    F.resize(nf);
    for (long i=0;i<nv;++i) {
        skip_ws(p); if (*p=='v') ++p;
        V[i].x = strtod(p,&p); V[i].y = strtod(p,&p); V[i].z = strtod(p,&p);
    }
    for (long i=0;i<nf;++i) {
        skip_ws(p); if (*p=='f') ++p;
        F[i][0] = (int)strtol(p,&p,10)-1;
        F[i][1] = (int)strtol(p,&p,10)-1;
        F[i][2] = (int)strtol(p,&p,10)-1;
    }
}
static void save_obj(const vector<Vec3>& V, const vector<array<int,3>>& F) {
    string out;
    out.reserve((size_t)V.size()*42 + (size_t)F.size()*28 + 64);
    char line[128];
    out.append(line, snprintf(line, sizeof(line), "%zu %zu\n", V.size(), F.size()));
    for (const Vec3& p: V) out.append(line, snprintf(line, sizeof(line), "v %.10g %.10g %.10g\n", p.x,p.y,p.z));
    for (const auto& f: F) out.append(line, snprintf(line, sizeof(line), "f %d %d %d\n", f[0]+1,f[1]+1,f[2]+1));
    fwrite(out.data(), 1, out.size(), stdout);
}

class Simplifier {
public:
    int n0, m0, aliveV, aliveF;
    vector<Vertex> vert;
    vector<Face> face;
    vector<int> nodeFace, nextNode;
    vector<int> markV, markC, markF, mapIndex;
    int stampV=1, stampC=1, stampF=1;
    double diag=0, eps=0, eps2=0;
    chrono::steady_clock::time_point t0;

    struct Item { double cost; int u,v,vu,vv; };
    struct Cmp { bool operator()(const Item& a, const Item& b) const { return a.cost > b.cost; } };
    priority_queue<Item, vector<Item>, Cmp> pq;

    Simplifier(const vector<Vec3>& V, const vector<array<int,3>>& F) {
        t0 = chrono::steady_clock::now();
        n0 = (int)V.size(); m0 = (int)F.size();
        aliveV = n0; aliveF = m0;
        vert.resize(n0); face.resize(m0);
        markV.assign(n0,0); markC.assign(n0,0); markF.assign(m0,0); mapIndex.assign(n0,-1);
        Vec3 mn(1e100,1e100,1e100), mx(-1e100,-1e100,-1e100);
        for (int i=0;i<n0;++i) {
            mn.x=min(mn.x,V[i].x); mn.y=min(mn.y,V[i].y); mn.z=min(mn.z,V[i].z);
            mx.x=max(mx.x,V[i].x); mx.y=max(mx.y,V[i].y); mx.z=max(mx.z,V[i].z);
            vert[i].p=V[i];
            for(int k=0;k<10;++k) vert[i].q[k]=0.0;
            vert[i].mn[0]=vert[i].mx[0]=V[i].x;
            vert[i].mn[1]=vert[i].mx[1]=V[i].y;
            vert[i].mn[2]=vert[i].mx[2]=V[i].z;
            vert[i].head=vert[i].tail=-1;
            vert[i].ver=1; vert[i].alive=1;
        }
        diag = norm3(mx-mn);
        eps = 0.05 * diag;
        eps2 = eps*eps*(1.0-1e-12);
        nodeFace.reserve((size_t)3*m0); nextNode.reserve((size_t)3*m0);
        for (int i=0;i<m0;++i) {
            face[i].v[0]=F[i][0]; face[i].v[1]=F[i][1]; face[i].v[2]=F[i][2]; face[i].alive=1;
            computeFaceNormal(i);
            addInitialQuadric(i);
            addIncident(F[i][0], i); addIncident(F[i][1], i); addIncident(F[i][2], i);
        }
    }

    double elapsed() const {
        return chrono::duration<double>(chrono::steady_clock::now()-t0).count();
    }

    void addIncident(int v, int fid) {
        int id=(int)nodeFace.size(); nodeFace.push_back(fid); nextNode.push_back(-1);
        if (vert[v].head==-1) vert[v].head=vert[v].tail=id;
        else { nextNode[vert[v].tail]=id; vert[v].tail=id; }
    }

    static inline bool contains(const Face& f, int v) { return f.v[0]==v || f.v[1]==v || f.v[2]==v; }

    bool computeFaceNormal(int fid) {
        Face& f=face[fid];
        Vec3 a=vert[f.v[0]].p, b=vert[f.v[1]].p, c=vert[f.v[2]].p;
        Vec3 cr=cross3(b-a,c-a); double l=norm3(cr);
        if (l<=1e-300) { f.nx=f.ny=f.nz=0; return false; }
        cr=cr/l; f.nx=(float)cr.x; f.ny=(float)cr.y; f.nz=(float)cr.z; return true;
    }

    void addPlane(int vi, const Vec3& n, double d, double w) {
        double a=n.x,b=n.y,c=n.z; double* q=vert[vi].q;
        q[0]+=w*a*a; q[1]+=w*a*b; q[2]+=w*a*c; q[3]+=w*a*d;
        q[4]+=w*b*b; q[5]+=w*b*c; q[6]+=w*b*d;
        q[7]+=w*c*c; q[8]+=w*c*d; q[9]+=w*d*d;
    }
    void addInitialQuadric(int fid) {
        Face& f=face[fid];
        Vec3 p0=vert[f.v[0]].p, p1=vert[f.v[1]].p, p2=vert[f.v[2]].p;
        Vec3 cr=cross3(p1-p0,p2-p0); double l=norm3(cr);
        if (l<=1e-300) return;
        Vec3 n=cr/l; double d=-dot3(n,p0);
        // Projected/normal-map error is closer to surface-integral error than to pure triangle count.
        double w=max(1e-12,0.5*l);
        addPlane(f.v[0],n,d,w); addPlane(f.v[1],n,d,w); addPlane(f.v[2],n,d,w);
    }

    static inline double qerr(const double q[10], const Vec3& p) {
        double x=p.x,y=p.y,z=p.z;
        return q[0]*x*x + 2*q[1]*x*y + 2*q[2]*x*z + 2*q[3]*x
             + q[4]*y*y + 2*q[5]*y*z + 2*q[6]*y
             + q[7]*z*z + 2*q[8]*z + q[9];
    }
    void mergedQ(int u,int v,double q[10]) const { for(int i=0;i<10;++i) q[i]=vert[u].q[i]+vert[v].q[i]; }
    void mergedBox(int u,int v,double mn[3],double mx[3]) const {
        for(int k=0;k<3;++k){ mn[k]=min(vert[u].mn[k],vert[v].mn[k]); mx[k]=max(vert[u].mx[k],vert[v].mx[k]); }
    }
    bool bboxOK(const double mn[3], const double mx[3], const Vec3& p) const {
        double dx=max(fabs(p.x-mn[0]),fabs(p.x-mx[0]));
        double dy=max(fabs(p.y-mn[1]),fabs(p.y-mx[1]));
        double dz=max(fabs(p.z-mn[2]),fabs(p.z-mx[2]));
        return dx*dx+dy*dy+dz*dz <= eps2;
    }
    static bool solveOptimal(const double q[10], Vec3& out) {
        double a00=q[0],a01=q[1],a02=q[2];
        double a10=q[1],a11=q[4],a12=q[5];
        double a20=q[2],a21=q[5],a22=q[7];
        double b0=-q[3],b1=-q[6],b2=-q[8];
        double det=a00*(a11*a22-a12*a21)-a01*(a10*a22-a12*a20)+a02*(a10*a21-a11*a20);
        double scale=fabs(a00)+fabs(a11)+fabs(a22)+fabs(a01)+fabs(a02)+fabs(a12);
        if (fabs(det)<=1e-14*max(1.0,scale*scale*scale)) return false;
        double id=1.0/det;
        double dx=b0*(a11*a22-a12*a21)-a01*(b1*a22-a12*b2)+a02*(b1*a21-a11*b2);
        double dy=a00*(b1*a22-a12*b2)-b0*(a10*a22-a12*a20)+a02*(a10*b2-b1*a20);
        double dz=a00*(a11*b2-b1*a21)-a01*(a10*b2-b1*a20)+b0*(a10*a21-a11*a20);
        out=Vec3(dx*id,dy*id,dz*id);
        return isfinite(out.x)&&isfinite(out.y)&&isfinite(out.z);
    }

    struct Eval { bool ok; Vec3 p; double cost; };

    bool affectedOK(int u,int v,const Vec3& cand) {
        const double MIN_AREA2=1e-28;
        const double MIN_DOT=0.030; // protects flat normal maps; still permissive enough for 8.8% target
        int fs=++stampF; if(stampF>INT_MAX-10){ fill(markF.begin(),markF.end(),0); stampF=1; fs=++stampF; }
        auto proc=[&](int head)->bool{
            for(int nd=head; nd!=-1; nd=nextNode[nd]){
                int fid=nodeFace[nd]; Face& f=face[fid];
                if(!f.alive || markF[fid]==fs) continue;
                bool hu=false,hv=false;
                for(int i=0;i<3;++i){ hu|=(f.v[i]==u); hv|=(f.v[i]==v); }
                if(!hu && !hv) continue;
                markF[fid]=fs;
                if(hu && hv) continue;
                Vec3 p[3];
                for(int i=0;i<3;++i){ int vi=f.v[i]; p[i]=(vi==u||vi==v)?cand:vert[vi].p; }
                Vec3 cr=cross3(p[1]-p[0],p[2]-p[0]); double a2=norm2(cr);
                if(!(a2>MIN_AREA2)) return false;
                Vec3 nn=cr*(1.0/sqrt(a2));
                double d=nn.x*f.nx + nn.y*f.ny + nn.z*f.nz;
                if(d<MIN_DOT) return false;
            }
            return true;
        };
        return proc(vert[u].head) && proc(vert[v].head);
    }

    bool candidateEval(int u,int v,bool full,Eval& ev) {
        double q[10], mn[3], mx[3]; mergedQ(u,v,q); mergedBox(u,v,mn,mx);
        Vec3 pu=vert[u].p, pv=vert[v].p;
        Vec3 opt;
        vector<Vec3> cand; cand.reserve(8);
        if(full && solveOptimal(q,opt)) cand.push_back(opt);
        cand.push_back(pu); cand.push_back(pv); cand.push_back((pu+pv)*0.5);
        cand.push_back(pu*0.75+pv*0.25); cand.push_back(pu*0.25+pv*0.75);
        cand.push_back(Vec3((mn[0]+mx[0])*0.5,(mn[1]+mx[1])*0.5,(mn[2]+mx[2])*0.5));
        double best=1e300; Vec3 bestP; bool any=false; double len2=norm2(pu-pv);
        for(size_t i=0;i<cand.size();++i){
            Vec3 p=cand[i]; if(!isfinite(p.x)||!isfinite(p.y)||!isfinite(p.z)) continue;
            bool dup=false; for(size_t j=0;j<i;++j) if(norm2(p-cand[j])<1e-26){ dup=true; break; }
            if(dup || !bboxOK(mn,mx,p)) continue;
            if(full && !affectedOK(u,v,p)) continue;
            double e=qerr(q,p); if(e<0 && e>-1e-14) e=0;
            double cost=e + 1e-11*len2;
            if(cost<best){ best=cost; bestP=p; any=true; }
        }
        ev.ok=any; ev.p=bestP; ev.cost=best; return any;
    }

    bool edgeLinkOK(int u,int v) {
        int sv=++stampV, sc=++stampC;
        if(stampV>INT_MAX-10 || stampC>INT_MAX-10){ fill(markV.begin(),markV.end(),0); fill(markC.begin(),markC.end(),0); stampV=stampC=1; sv=++stampV; sc=++stampC; }
        int edgeFaces=0;
        for(int nd=vert[u].head; nd!=-1; nd=nextNode[nd]){
            int fid=nodeFace[nd]; Face& f=face[fid]; if(!f.alive) continue;
            bool hu=false,hv=false; for(int i=0;i<3;++i){ hu|=(f.v[i]==u); hv|=(f.v[i]==v); }
            if(!hu) continue;
            for(int i=0;i<3;++i){ int w=f.v[i]; if(w!=u && vert[w].alive) markV[w]=sv; }
            if(hv) ++edgeFaces;
        }
        if(edgeFaces!=2) return false;
        int common=0;
        for(int nd=vert[v].head; nd!=-1; nd=nextNode[nd]){
            int fid=nodeFace[nd]; Face& f=face[fid]; if(!f.alive) continue;
            bool hv=false; for(int i=0;i<3;++i) hv|=(f.v[i]==v);
            if(!hv) continue;
            for(int i=0;i<3;++i){
                int w=f.v[i]; if(w==v || w==u || !vert[w].alive) continue;
                if(markV[w]==sv && markC[w]!=sc){ markC[w]=sc; ++common; if(common>2) return false; }
            }
        }
        return common==2;
    }

    static inline uint64_t keyEdge(int a,int b){ if(a>b) swap(a,b); return (uint64_t)(uint32_t)a<<32 | (uint32_t)b; }
    void pushEdge(int u,int v){
        if(u==v || u<0 || v<0) return; if(u>v) swap(u,v);
        if(!vert[u].alive || !vert[v].alive) return;
        Eval ev; if(!candidateEval(u,v,false,ev)) return;
        pq.push(Item{ev.cost,u,v,vert[u].ver,vert[v].ver});
    }
    void buildPQ(){
        vector<uint64_t> edges; edges.reserve((size_t)3*m0);
        for(int i=0;i<m0;++i){ int a=face[i].v[0],b=face[i].v[1],c=face[i].v[2]; edges.push_back(keyEdge(a,b)); edges.push_back(keyEdge(b,c)); edges.push_back(keyEdge(c,a)); }
        sort(edges.begin(),edges.end()); edges.erase(unique(edges.begin(),edges.end()),edges.end());
        for(uint64_t k:edges) pushEdge((int)(k>>32),(int)(k&0xffffffffu));
    }

    void appendList(int keep,int kill){
        if(vert[kill].head==-1) return;
        if(vert[keep].head==-1){ vert[keep].head=vert[kill].head; vert[keep].tail=vert[kill].tail; }
        else { nextNode[vert[keep].tail]=vert[kill].head; vert[keep].tail=vert[kill].tail; }
        vert[kill].head=vert[kill].tail=-1;
    }
    void recomputeAround(int keep){
        int fs=++stampF; if(stampF>INT_MAX-10){ fill(markF.begin(),markF.end(),0); stampF=1; fs=++stampF; }
        for(int nd=vert[keep].head; nd!=-1; nd=nextNode[nd]){
            int fid=nodeFace[nd]; Face& f=face[fid];
            if(!f.alive || markF[fid]==fs) continue; markF[fid]=fs;
            if(contains(f,keep)) computeFaceNormal(fid);
        }
    }
    void pushAdjacent(int keep){
        int sv=++stampV; if(stampV>INT_MAX-10){ fill(markV.begin(),markV.end(),0); stampV=1; sv=++stampV; }
        markV[keep]=sv;
        for(int nd=vert[keep].head; nd!=-1; nd=nextNode[nd]){
            int fid=nodeFace[nd]; Face& f=face[fid]; if(!f.alive || !contains(f,keep)) continue;
            for(int i=0;i<3;++i){ int w=f.v[i]; if(w!=keep && vert[w].alive && markV[w]!=sv){ markV[w]=sv; pushEdge(keep,w); } }
        }
    }
    void collapse(int keep,int kill,const Vec3& p){
        vert[keep].p=p;
        for(int i=0;i<10;++i) vert[keep].q[i]+=vert[kill].q[i];
        for(int k=0;k<3;++k){ vert[keep].mn[k]=min(vert[keep].mn[k],vert[kill].mn[k]); vert[keep].mx[k]=max(vert[keep].mx[k],vert[kill].mx[k]); }
        for(int nd=vert[kill].head; nd!=-1; nd=nextNode[nd]){
            int fid=nodeFace[nd]; Face& f=face[fid]; if(!f.alive) continue;
            bool hk=false,hl=false; for(int i=0;i<3;++i){ hk|=(f.v[i]==keep); hl|=(f.v[i]==kill); }
            if(!hl) continue;
            if(hk){ f.alive=0; --aliveF; }
            else for(int i=0;i<3;++i) if(f.v[i]==kill) f.v[i]=keep;
        }
        appendList(keep,kill);
        vert[kill].alive=0; ++vert[keep].ver; ++vert[kill].ver; --aliveV;
        recomputeAround(keep);
        pushAdjacent(keep);
    }

    MeshOut compact(){
        MeshOut out; out.V.reserve(aliveV); out.F.reserve(aliveF); fill(mapIndex.begin(),mapIndex.end(),-1);
        auto idOf=[&](int old)->int{ int& r=mapIndex[old]; if(r==-1){ r=(int)out.V.size(); out.V.push_back(vert[old].p); } return r; };
        for(int i=0;i<m0;++i){
            Face& f=face[i]; if(!f.alive) continue; int a=f.v[0],b=f.v[1],c=f.v[2];
            if(a==b||b==c||c==a) continue; if(!vert[a].alive||!vert[b].alive||!vert[c].alive) continue;
            Vec3 pa=vert[a].p,pb=vert[b].p,pc=vert[c].p; if(norm2(cross3(pb-pa,pc-pa))<=1e-28) continue;
            out.F.push_back({idOf(a),idOf(b),idOf(c)});
        }
        return out;
    }

    MeshOut runSmallPlanar(){
        buildPQ(); const double ZERO=1e-8; int guard=0;
        while(!pq.empty() && guard<n0*4){
            Item it=pq.top(); pq.pop();
            if(!vert[it.u].alive||!vert[it.v].alive) continue;
            if(vert[it.u].ver!=it.vu||vert[it.v].ver!=it.vv) continue;
            if(!edgeLinkOK(it.u,it.v)) continue;
            Eval ev; if(!candidateEval(it.u,it.v,true,ev)) continue;
            if(ev.cost>ZERO) break;
            collapse(it.u,it.v,ev.p); ++guard;
        }
        return compact();
    }

    MeshOut run(){
        if(n0<=50) return runSmallPlanar();
        buildPQ();
        const double TARGET_RATIO = 0.0880; // retained vertices -> 91.20 compression when reached
        int target=max(4,(int)ceil(n0*TARGET_RATIO));
        const double SOFT_LIMIT=18.6;
        long long pops=0;
        while(aliveV>target && !pq.empty()){
            if((pops++ & 8191LL)==0 && elapsed()>SOFT_LIMIT) break;
            Item it=pq.top(); pq.pop();
            if(!vert[it.u].alive||!vert[it.v].alive) continue;
            if(vert[it.u].ver!=it.vu||vert[it.v].ver!=it.vv) continue;
            if(!edgeLinkOK(it.u,it.v)) continue;
            Eval ev; if(!candidateEval(it.u,it.v,true,ev)) continue;
            collapse(it.u,it.v,ev.p);
        }
        return compact();
    }
};

int main(){
    vector<Vec3> V; vector<array<int,3>> F; load_obj(V,F);
    if(V.empty()||F.empty()){ save_obj(V,F); return 0; }
    Simplifier s(V,F);
    MeshOut out=s.run();
    if(out.V.empty()||out.F.empty()) save_obj(V,F);
    else save_obj(out.V,out.F);
    return 0;
}
