#include <bits/stdc++.h>
using namespace std;

struct FastInput {
    vector<char> buf; char* p;
    FastInput() {
        buf.reserve(1 << 27);
        char chunk[1 << 16]; size_t n;
        while ((n = fread(chunk, 1, sizeof(chunk), stdin)) > 0) buf.insert(buf.end(), chunk, chunk + n);
        buf.push_back('\0'); p = buf.data();
    }
    inline void skipws(){ while(*p==' '||*p=='\n'||*p=='\r'||*p=='\t') ++p; }
    inline long readLong(){ skipws(); return strtol(p, &p, 10); }
    inline double readDouble(){ skipws(); return strtod(p, &p); }
    inline char readChar(){ skipws(); return *p++; }
};

struct Vec3 {
    double x,y,z;
    Vec3():x(0),y(0),z(0){}
    Vec3(double X,double Y,double Z):x(X),y(Y),z(Z){}
    inline Vec3 operator+(const Vec3& o) const { return {x+o.x,y+o.y,z+o.z}; }
    inline Vec3 operator-(const Vec3& o) const { return {x-o.x,y-o.y,z-o.z}; }
    inline Vec3 operator*(double s) const { return {x*s,y*s,z*s}; }
    inline Vec3 operator/(double s) const { return {x/s,y/s,z/s}; }
};
static inline double dotp(const Vec3&a,const Vec3&b){ return a.x*b.x+a.y*b.y+a.z*b.z; }
static inline Vec3 crossp(const Vec3&a,const Vec3&b){ return {a.y*b.z-a.z*b.y, a.z*b.x-a.x*b.z, a.x*b.y-a.y*b.x}; }
static inline double norm2(const Vec3&a){ return dotp(a,a); }
static inline double normv(const Vec3&a){ return sqrt(norm2(a)); }
static inline double distv(const Vec3&a,const Vec3&b){ return normv(a-b); }

struct Quadric {
    // xx,xy,xz,xw,yy,yz,yw,zz,zw,ww
    double q[10];
    Quadric(){ memset(q,0,sizeof(q)); }
    inline void addPlane(double a,double b,double c,double d,double w){
        q[0]+=w*a*a; q[1]+=w*a*b; q[2]+=w*a*c; q[3]+=w*a*d;
        q[4]+=w*b*b; q[5]+=w*b*c; q[6]+=w*b*d;
        q[7]+=w*c*c; q[8]+=w*c*d; q[9]+=w*d*d;
    }
    inline void operator+=(const Quadric& o){ for(int i=0;i<10;i++) q[i]+=o.q[i]; }
    inline double eval(const Vec3& p) const {
        double x=p.x,y=p.y,z=p.z;
        return q[0]*x*x + 2*q[1]*x*y + 2*q[2]*x*z + 2*q[3]*x
             + q[4]*y*y + 2*q[5]*y*z + 2*q[6]*y
             + q[7]*z*z + 2*q[8]*z + q[9];
    }
};
static inline Quadric qsum(const Quadric&a,const Quadric&b){ Quadric r=a; r+=b; return r; }

struct Face { int v[3]; unsigned char alive; };
struct Vertex {
    Vec3 p;
    Quadric q;
    double qw = 0.0;
    double rad = 0.0;
    float imp = 0.0f;
    int ver = 0;
    unsigned char alive = 1;
    vector<int> inc;
};

static vector<Vertex> V;
static vector<Face> F;
static int N0, M0, liveV, liveF;
static double bboxDiag, epsHaus, epsBase, minArea2Global;
static double meanSharp = 0.0, meanSqrtSharp = 0.0;
static vector<int> vMark, fMark;
static int vStamp=1, fStamp=1;

static inline bool face_has(const Face& f,int a){ return f.v[0]==a || f.v[1]==a || f.v[2]==a; }
static inline int face_opp(const Face& f,int a,int b){ for(int i=0;i<3;i++) if(f.v[i]!=a && f.v[i]!=b) return f.v[i]; return -1; }
static inline bool face_degenerate_indices(const Face& f){ return f.v[0]==f.v[1] || f.v[0]==f.v[2] || f.v[1]==f.v[2]; }
static inline Vec3 face_normal_raw(const Face& f){ return crossp(V[f.v[1]].p - V[f.v[0]].p, V[f.v[2]].p - V[f.v[0]].p); }
static inline double tri_area2_pos(const Vec3&a,const Vec3&b,const Vec3&c){ return norm2(crossp(b-a,c-a)); }
static inline uint64_t edge_key(int a,int b){ if(a>b) swap(a,b); return ((uint64_t)(uint32_t)a<<32) | (uint32_t)b; }
static inline int key_a(uint64_t k){ return (int)(k>>32); }
static inline int key_b(uint64_t k){ return (int)(k & 0xffffffffu); }

static bool solve_optimal(const Quadric& Q, Vec3& out){
    double a00=Q.q[0], a01=Q.q[1], a02=Q.q[2];
    double a11=Q.q[4], a12=Q.q[5], a22=Q.q[7];
    double b0=-Q.q[3], b1=-Q.q[6], b2=-Q.q[8];
    double det = a00*(a11*a22-a12*a12) - a01*(a01*a22-a12*a02) + a02*(a01*a12-a11*a02);
    double sc = fabs(a00)+fabs(a11)+fabs(a22)+fabs(a01)+fabs(a02)+fabs(a12)+1e-30;
    if(fabs(det) < 1e-12*sc*sc*sc) return false;
    double id=1.0/det;
    double dx = b0*(a11*a22-a12*a12) - a01*(b1*a22-a12*b2) + a02*(b1*a12-a11*b2);
    double dy = a00*(b1*a22-a12*b2) - b0*(a01*a22-a12*a02) + a02*(a01*b2-b1*a02);
    double dz = a00*(a11*b2-b1*a12) - a01*(a01*b2-b1*a02) + b0*(a01*a12-a11*a02);
    out = {dx*id, dy*id, dz*id};
    if(!isfinite(out.x)||!isfinite(out.y)||!isfinite(out.z)) return false;
    if(norm2(out) > 4.25) return false;
    return true;
}

static inline void compact_inc(int u){
    if(++fStamp == INT_MAX){ fill(fMark.begin(), fMark.end(), 0); fStamp=1; }
    vector<int> nw; nw.reserve(min<size_t>(V[u].inc.size(), 96));
    for(int fid: V[u].inc){
        if(fid<0||fid>=M0) continue;
        if(fMark[fid]==fStamp) continue;
        fMark[fid]=fStamp;
        if(F[fid].alive && face_has(F[fid],u)) nw.push_back(fid);
    }
    V[u].inc.swap(nw);
}

static int edge_faces(int u,int v,int ef[4],int opp[4]){
    const vector<int>& A = (V[u].inc.size() <= V[v].inc.size()) ? V[u].inc : V[v].inc;
    if(++fStamp == INT_MAX){ fill(fMark.begin(), fMark.end(), 0); fStamp=1; }
    int cnt=0;
    for(int fid: A){
        if(fid<0||fid>=M0) continue;
        if(fMark[fid]==fStamp) continue;
        fMark[fid]=fStamp;
        Face &f=F[fid];
        if(!f.alive) continue;
        if(face_has(f,u) && face_has(f,v)){
            if(cnt<4){ ef[cnt]=fid; opp[cnt]=face_opp(f,u,v); }
            ++cnt; if(cnt>3) break;
        }
    }
    return cnt;
}

static bool link_condition(int u,int v,int o0,int o1){
    if(o0<0 || o1<0 || o0==o1) return false;
    if(vStamp > INT_MAX-8){ fill(vMark.begin(), vMark.end(), 0); vStamp=1; }
    int st=vStamp++;
    for(int fid: V[u].inc){
        if(!F[fid].alive || !face_has(F[fid],u)) continue;
        Face &f=F[fid];
        for(int k=0;k<3;k++){ int w=f.v[k]; if(w!=u && w!=v && V[w].alive) vMark[w]=st; }
    }
    int common[8], cnt=0;
    for(int fid: V[v].inc){
        if(!F[fid].alive || !face_has(F[fid],v)) continue;
        Face &f=F[fid];
        for(int k=0;k<3;k++){
            int w=f.v[k];
            if(w==u||w==v||!V[w].alive) continue;
            if(vMark[w]==st){ if(cnt<8) common[cnt]=w; ++cnt; vMark[w] = -st; }
        }
    }
    if(cnt!=2) return false;
    return (common[0]==o0 && common[1]==o1) || (common[0]==o1 && common[1]==o0);
}

struct Candidate { bool ok=false; Vec3 p; double rad=0,qerr=0,ndev=0,cost=0; };

static bool normal_check(int u,int v,const Vec3& np,double& maxDev){
    maxDev=0.0;
    if(++fStamp == INT_MAX){ fill(fMark.begin(), fMark.end(), 0); fStamp=1; }
    auto scan = [&](int src)->bool{
        for(int fid: V[src].inc){
            if(fid<0||fid>=M0) continue;
            if(fMark[fid]==fStamp) continue;
            fMark[fid]=fStamp;
            Face &f=F[fid]; if(!f.alive) continue;
            bool hu=face_has(f,u), hv=face_has(f,v);
            if(!hu && !hv) continue;
            if(hu && hv) continue;
            Vec3 p0=V[f.v[0]].p, p1=V[f.v[1]].p, p2=V[f.v[2]].p;
            Vec3 oldn=crossp(p1-p0,p2-p0);
            double oldl=normv(oldn); if(!(oldl>0)) return false;
            if(f.v[0]==u || f.v[0]==v) p0=np;
            if(f.v[1]==u || f.v[1]==v) p1=np;
            if(f.v[2]==u || f.v[2]==v) p2=np;
            Vec3 newn=crossp(p1-p0,p2-p0);
            double newl2=norm2(newn); if(newl2 <= minArea2Global) return false;
            double newl=sqrt(newl2);
            double d=dotp(oldn,newn)/(oldl*newl);
            if(d < 0.08) return false;
            if(d < 1.0) maxDev=max(maxDev,1.0-d);
        }
        return true;
    };
    if(!scan(u)) return false;
    if(!scan(v)) return false;
    return true;
}

static Candidate compute_candidate(int u,int v,bool fullCheck){
    Candidate best;
    if(u==v || !V[u].alive || !V[v].alive) return best;
    if(fullCheck){
        int ef[4],opp[4];
        int cnt=edge_faces(u,v,ef,opp);
        if(cnt!=2) return best;
        if(!link_condition(u,v,opp[0],opp[1])) return best;
    }
    Quadric Q=qsum(V[u].q,V[v].q);
    double qw=V[u].qw+V[v].qw+1e-30;
    double imp=max((double)V[u].imp,(double)V[v].imp);
    double localLimit=epsBase*(0.990 - 0.46*imp);
    localLimit=max(localLimit, epsBase*0.40);
    localLimit=min(localLimit, epsHaus*0.999);

    Vec3 pts[7]; int nc=0; Vec3 opt;
    if(solve_optimal(Q,opt)) pts[nc++]=opt;
    pts[nc++]=V[u].p; pts[nc++]=V[v].p;
    pts[nc++]=(V[u].p+V[v].p)*0.5;
    pts[nc++]=V[u].p*0.75 + V[v].p*0.25;
    pts[nc++]=V[u].p*0.25 + V[v].p*0.75;
    for(int i=0;i<nc;i++){
        Vec3 p=pts[i];
        double rad=max(distv(p,V[u].p)+V[u].rad, distv(p,V[v].p)+V[v].rad);
        if(rad > localLimit + 1e-12) continue;
        double ev=Q.eval(p);
        if(ev < 0 && ev > -1e-17) ev=0;
        if(ev < -1e-12 || !isfinite(ev)) continue;
        double qerr=sqrt(max(0.0,ev/qw));
        double ndev=0.0;
        if(fullCheck && !normal_check(u,v,p,ndev)) continue;
        double cost = qerr*(1.0+2.6*imp) + 0.012*epsBase*ndev*(1.0+imp) + 1e-7*rad/(epsBase+1e-30);
        if(!best.ok || cost < best.cost){ best.ok=true; best.p=p; best.rad=rad; best.qerr=qerr; best.ndev=ndev; best.cost=cost; }
    }
    return best;
}

struct Item { double cost; int u,v,vu,vv; bool operator<(const Item& o) const { return cost > o.cost; } };
static priority_queue<Item> pq;

static inline void push_edge(int u,int v){
    if(u==v || !V[u].alive || !V[v].alive) return;
    Candidate c=compute_candidate(u,v,false);
    if(c.ok) pq.push({c.cost,u,v,V[u].ver,V[v].ver});
}

static bool collapse_edge(int u,int v,const Candidate& c){
    int ef[4],opp[4];
    int cnt=edge_faces(u,v,ef,opp);
    if(cnt!=2) return false;
    if(!link_condition(u,v,opp[0],opp[1])) return false;
    for(int i=0;i<2;i++) if(F[ef[i]].alive){ F[ef[i]].alive=0; --liveF; }
    for(int fid: V[v].inc){
        if(fid<0||fid>=M0) continue;
        Face &f=F[fid]; if(!f.alive) continue;
        bool changed=false;
        for(int k=0;k<3;k++) if(f.v[k]==v){ f.v[k]=u; changed=true; }
        if(changed){ if(face_degenerate_indices(f)){ f.alive=0; --liveF; } else V[u].inc.push_back(fid); }
    }
    V[u].p=c.p; V[u].rad=c.rad; V[u].q+=V[v].q; V[u].qw+=V[v].qw; V[u].imp=max(V[u].imp,V[v].imp);
    V[v].alive=0; ++V[u].ver; ++V[v].ver; --liveV;
    compact_inc(u);
    if(vStamp > INT_MAX-8){ fill(vMark.begin(), vMark.end(), 0); vStamp=1; }
    int st=vStamp++;
    for(int fid: V[u].inc){
        if(!F[fid].alive) continue;
        Face &f=F[fid];
        for(int k=0;k<3;k++){ int w=f.v[k]; if(w!=u && V[w].alive && vMark[w]!=st){ vMark[w]=st; push_edge(u,w); } }
    }
    return true;
}

static void build_quadrics_and_edges(){
    vector<int> deg(N0,0);
    for(int i=0;i<M0;i++){ deg[F[i].v[0]]++; deg[F[i].v[1]]++; deg[F[i].v[2]]++; }
    for(int i=0;i<N0;i++) V[i].inc.reserve(deg[i]+4);
    for(int i=0;i<M0;i++) for(int k=0;k<3;k++) V[F[i].v[k]].inc.push_back(i);

    vector<Vec3> fNorm(M0); vector<double> fArea(M0);
    double totalArea=0.0;
    for(int i=0;i<M0;i++){
        Vec3 cr=face_normal_raw(F[i]); double len=normv(cr), area=0.5*len; fArea[i]=area; totalArea+=area;
        fNorm[i]=(len>0? cr/len : Vec3());
    }
    double avgArea=totalArea/max(1,M0), wFloor=max(avgArea*0.03,1e-16);
    for(int i=0;i<M0;i++){
        Vec3 n=fNorm[i]; if(norm2(n)==0) continue;
        double d=-dotp(n,V[F[i].v[0]].p), w=max(fArea[i],wFloor);
        Quadric q; q.addPlane(n.x,n.y,n.z,d,w);
        for(int k=0;k<3;k++){ int a=F[i].v[k]; V[a].q+=q; V[a].qw+=w; }
    }
    vector<pair<uint64_t,int>> refs; refs.reserve((size_t)M0*3);
    for(int i=0;i<M0;i++){
        int a=F[i].v[0], b=F[i].v[1], c=F[i].v[2];
        refs.emplace_back(edge_key(a,b),i); refs.emplace_back(edge_key(b,c),i); refs.emplace_back(edge_key(c,a),i);
    }
    sort(refs.begin(), refs.end(), [](const auto&x,const auto&y){return x.first<y.first;});
    const Vec3 axes[3]={{1,0,0},{0,1,0},{0,0,1}};
    double sumSharp=0.0, sumSqrt=0.0; long long eCnt=0;
    for(size_t i=0;i<refs.size();){
        size_t j=i+1; while(j<refs.size() && refs[j].first==refs[i].first) ++j;
        uint64_t key=refs[i].first; int a=key_a(key), b=key_b(key);
        if(j-i==2){
            int f0=refs[i].second, f1=refs[i+1].second;
            Vec3 n0=fNorm[f0], n1=fNorm[f1];
            double d=dotp(n0,n1); if(d>1)d=1; if(d<-1)d=-1;
            double sharp=1.0 - fabs(d);
            if(sharp<0) sharp=0;
            sumSharp += sharp; sumSqrt += sqrt(sharp); ++eCnt;
            double sil=0.0;
            for(const Vec3& ax: axes){ double s0=dotp(n0,ax), s1=dotp(n1,ax); if(s0*s1 < -1e-7) sil += 0.35; }
            double val=min(1.0, 0.70*sharp + sil);
            if(val>0){ V[a].imp=max(V[a].imp,(float)val); V[b].imp=max(V[b].imp,(float)val); }
        } else {
            V[a].imp=V[b].imp=1.0f;
        }
        push_edge(a,b);
        i=j;
    }
    if(eCnt>0){ meanSharp=sumSharp/(double)eCnt; meanSqrtSharp=sumSqrt/(double)eCnt; }
}

static int choose_target_vertices(){
    if(N0<=4) return N0;
    double minAgg;
    if(N0 < 200) minAgg = 0.50;
    else if(N0 < 3000) minAgg = 0.070;
    else if(N0 < 30000) minAgg = 0.050;
    else if(N0 < 100000) minAgg = 0.048;
    else if(N0 < 450000) minAgg = 0.044;
    else minAgg = 0.040;

    // Normal-map complexity proxy: for the same smooth object this value falls as tessellation density grows.
    double visual = 5.75 * meanSqrtSharp - 0.020;
    if(meanSqrtSharp < 0.010) visual *= 0.55;
    double ratio = max(minAgg, visual);
    if(N0 < 3000) ratio = min(ratio, 0.34);
    else if(N0 < 10000) ratio = min(ratio, 0.24);
    else if(N0 < 30000) ratio = min(ratio, 0.22);
    else if(N0 < 100000) ratio = min(ratio, 0.16);
    else if(N0 < 450000) ratio = min(ratio, 0.12);
    else ratio = min(ratio, 0.10);
    int target=(int)ceil(N0*ratio);
    if(target<4) target=4;
    return min(target,N0);
}

static void simplify(){
    liveV=N0; liveF=M0;
    vMark.assign(N0,0); fMark.assign(M0,0);
    build_quadrics_and_edges();
    int target=choose_target_vertices();
    const double hardQErr = epsBase * 0.58;
    const double exactErr = max(1e-12, epsBase*1e-7);
    while(!pq.empty()){
        Item it=pq.top(); pq.pop();
        if(it.u<0||it.v<0||it.u>=N0||it.v>=N0) continue;
        if(!V[it.u].alive || !V[it.v].alive || V[it.u].ver!=it.vu || V[it.v].ver!=it.vv) continue;
        Candidate c=compute_candidate(it.u,it.v,true);
        if(!c.ok) continue;
        if(c.cost > it.cost*1.35 + 1e-12){ pq.push({c.cost,it.u,it.v,V[it.u].ver,V[it.v].ver}); continue; }
        bool exactish = (c.qerr <= exactErr && c.ndev <= 1e-9);
        if(liveV <= target && !exactish) break;
        if(c.qerr > hardQErr && !exactish) {
            if(liveV <= target) break;
            continue;
        }
        collapse_edge(it.u,it.v,c);
    }
}

static void load_obj(){
    FastInput in;
    N0=(int)in.readLong(); M0=(int)in.readLong();
    V.assign(N0, Vertex()); F.assign(M0, Face());
    double xmin=1e100,ymin=1e100,zmin=1e100,xmax=-1e100,ymax=-1e100,zmax=-1e100;
    for(int i=0;i<N0;i++){
        (void)in.readChar(); double x=in.readDouble(), y=in.readDouble(), z=in.readDouble();
        V[i].p={x,y,z}; xmin=min(xmin,x); ymin=min(ymin,y); zmin=min(zmin,z); xmax=max(xmax,x); ymax=max(ymax,y); zmax=max(zmax,z);
    }
    for(int i=0;i<M0;i++){
        (void)in.readChar();
        F[i].v[0]=(int)in.readLong()-1; F[i].v[1]=(int)in.readLong()-1; F[i].v[2]=(int)in.readLong()-1; F[i].alive=1;
    }
    double lx=xmax-xmin, ly=ymax-ymin, lz=zmax-zmin;
    bboxDiag=sqrt(lx*lx+ly*ly+lz*lz); if(!(bboxDiag>0)) bboxDiag=1.0;
    epsHaus=0.05*bboxDiag; epsBase=max(epsHaus,1e-12);
    minArea2Global=max(1e-30, bboxDiag*bboxDiag*bboxDiag*bboxDiag*1e-24);
}

static void save_obj(){
    vector<int> id(N0,-1); int nv=0;
    for(int i=0;i<N0;i++) if(V[i].alive) id[i]=nv++;
    vector<array<int,3>> outF; outF.reserve(liveF);
    for(int i=0;i<M0;i++) if(F[i].alive){
        int a=F[i].v[0], b=F[i].v[1], c=F[i].v[2];
        if(a<0||b<0||c<0||a>=N0||b>=N0||c>=N0) continue;
        if(!V[a].alive||!V[b].alive||!V[c].alive) continue;
        if(a==b||b==c||c==a) continue;
        if(tri_area2_pos(V[a].p,V[b].p,V[c].p) <= minArea2Global) continue;
        outF.push_back({id[a]+1,id[b]+1,id[c]+1});
    }
    string out; out.reserve((size_t)nv*42 + outF.size()*30 + 64);
    char line[128];
    out.append(line, snprintf(line,sizeof(line), "%d %zu\n", nv, outF.size()));
    for(int i=0;i<N0;i++) if(V[i].alive) out.append(line, snprintf(line,sizeof(line), "v %.10g %.10g %.10g\n", V[i].p.x,V[i].p.y,V[i].p.z));
    for(auto &f: outF) out.append(line, snprintf(line,sizeof(line), "f %d %d %d\n", f[0],f[1],f[2]));
    fwrite(out.data(),1,out.size(),stdout);
}

int main(){
    load_obj();
    simplify();
    save_obj();
    return 0;
}
