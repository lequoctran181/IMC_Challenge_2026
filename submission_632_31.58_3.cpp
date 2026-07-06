#include <bits/stdc++.h>
using namespace std;

struct Vec3 {
    float x, y, z;
    Vec3() : x(0), y(0), z(0) {}
    Vec3(float X, float Y, float Z) : x(X), y(Y), z(Z) {}
};
static inline Vec3 operator+(const Vec3& a, const Vec3& b){ return Vec3(a.x+b.x,a.y+b.y,a.z+b.z); }
static inline Vec3 operator-(const Vec3& a, const Vec3& b){ return Vec3(a.x-b.x,a.y-b.y,a.z-b.z); }
static inline Vec3 operator*(const Vec3& a, float s){ return Vec3(a.x*s,a.y*s,a.z*s); }
static inline double dotd(const Vec3& a, const Vec3& b){ return (double)a.x*b.x + (double)a.y*b.y + (double)a.z*b.z; }
static inline Vec3 crossf(const Vec3& a, const Vec3& b){ return Vec3(a.y*b.z-a.z*b.y, a.z*b.x-a.x*b.z, a.x*b.y-a.y*b.x); }
static inline double norm2d(const Vec3& a){ return dotd(a,a); }
static inline double distd(const Vec3& a, const Vec3& b){ return sqrt(norm2d(a-b)); }

struct Quadric {
    double q[10];
    Quadric(){ memset(q,0,sizeof(q)); }
    inline void add(const Quadric& o){ for(int i=0;i<10;i++) q[i]+=o.q[i]; }
    inline double eval(const Vec3& p) const {
        double x=p.x,y=p.y,z=p.z;
        return q[0]*x*x + 2*q[1]*x*y + 2*q[2]*x*z + 2*q[3]*x
             + q[4]*y*y + 2*q[5]*y*z + 2*q[6]*y
             + q[7]*z*z + 2*q[8]*z + q[9];
    }
    inline void addPlane(double a,double b,double c,double d,double w){
        q[0]+=w*a*a; q[1]+=w*a*b; q[2]+=w*a*c; q[3]+=w*a*d;
        q[4]+=w*b*b; q[5]+=w*b*c; q[6]+=w*b*d;
        q[7]+=w*c*c; q[8]+=w*c*d; q[9]+=w*d*d;
    }
    inline void addPoint(const Vec3& p, double w){
        double x=p.x,y=p.y,z=p.z;
        q[0]+=w; q[4]+=w; q[7]+=w;
        q[3]+=-w*x; q[6]+=-w*y; q[8]+=-w*z;
        q[9]+=w*(x*x+y*y+z*z);
    }
};

struct Face {
    int v[3];
    float nx, ny, nz;
    float area;
    float w;
    unsigned char alive;
};

struct Tri { int a,b,c; };

static int N, M;
static vector<Vec3> P;
static vector<Face> faces;
static vector<Tri> origFaces;
static vector<vector<int>> inc;
static vector<Quadric> Q;
static vector<float> radiusRep;
static vector<int> clusterSz;
static vector<unsigned char> vAlive;
static vector<float> pointWeight;
static int aliveVertices, aliveFaces;
static double hausTol;
static double avgAreaGlobal;

static vector<char> slurp_stdin(){
    vector<char> buf;
    buf.reserve(1<<26);
    char tmp[1<<16];
    size_t n;
    while((n=fread(tmp,1,sizeof(tmp),stdin))>0) buf.insert(buf.end(),tmp,tmp+n);
    buf.push_back('\0');
    return buf;
}
static inline void skipws(char*& p){ while(*p==' '||*p=='\n'||*p=='\r'||*p=='\t') ++p; }

static void load_input(){
    vector<char> buf=slurp_stdin();
    char* p=buf.data();
    N=(int)strtol(p,&p,10);
    M=(int)strtol(p,&p,10);
    P.resize(N);
    float xmin=1e9f,ymin=1e9f,zmin=1e9f,xmax=-1e9f,ymax=-1e9f,zmax=-1e9f;
    for(int i=0;i<N;i++){
        skipws(p); if(*p=='v') ++p;
        double x=strtod(p,&p), y=strtod(p,&p), z=strtod(p,&p);
        P[i]=Vec3((float)x,(float)y,(float)z);
        xmin=min(xmin,(float)x); ymin=min(ymin,(float)y); zmin=min(zmin,(float)z);
        xmax=max(xmax,(float)x); ymax=max(ymax,(float)y); zmax=max(zmax,(float)z);
    }
    faces.resize(M);
    origFaces.resize(M);
    vector<int> deg(N,0);
    for(int i=0;i<M;i++){
        skipws(p); if(*p=='f') ++p;
        int a=(int)strtol(p,&p,10)-1;
        int b=(int)strtol(p,&p,10)-1;
        int c=(int)strtol(p,&p,10)-1;
        faces[i].v[0]=a; faces[i].v[1]=b; faces[i].v[2]=c;
        faces[i].nx=faces[i].ny=faces[i].nz=0; faces[i].area=0; faces[i].w=1; faces[i].alive=1;
        origFaces[i]={a,b,c};
        deg[a]++; deg[b]++; deg[c]++;
    }
    inc.resize(N);
    for(int i=0;i<N;i++) inc[i].reserve(deg[i]+8);
    for(int i=0;i<M;i++){
        inc[faces[i].v[0]].push_back(i);
        inc[faces[i].v[1]].push_back(i);
        inc[faces[i].v[2]].push_back(i);
    }
    double dx=(double)xmax-xmin, dy=(double)ymax-ymin, dz=(double)zmax-zmin;
    hausTol=0.05*sqrt(dx*dx+dy*dy+dz*dz)*0.985;
    if(hausTol<=0) hausTol=1e-9;
}

static inline bool face_has_vertex(const Face& f, int v){ return f.v[0]==v || f.v[1]==v || f.v[2]==v; }
static inline bool face_has_edge(const Face& f, int a, int b){ return face_has_vertex(f,a) && face_has_vertex(f,b); }
static inline int face_other_of_edge(const Face& f, int a, int b){
    for(int k=0;k<3;k++) if(f.v[k]!=a && f.v[k]!=b) return f.v[k];
    return -1;
}
static inline bool face_has_two(const Face& f, int a, int b){ return face_has_vertex(f,a) && face_has_vertex(f,b); }

static bool compute_face_geom(int fid){
    Face &f=faces[fid];
    Vec3 a=P[f.v[0]], b=P[f.v[1]], c=P[f.v[2]];
    Vec3 cr=crossf(b-a,c-a);
    double len=sqrt(norm2d(cr));
    if(len<=0){ f.nx=f.ny=f.nz=0; f.area=0; return false; }
    f.nx=(float)(cr.x/len); f.ny=(float)(cr.y/len); f.nz=(float)(cr.z/len);
    f.area=(float)(0.5*len);
    return true;
}

static void compute_visual_marks(vector<unsigned short>& marks){
    marks.assign(N,0);
    if(N==0) return;
    double target = max(4.0, 0.08 * (double)N);
    int R = (int)(sqrt(max(1.0, target*0.55/6.0)) + 0.5);
    R = max(12, min(96, R));
    const double D=2.5, FOC=800.0;
    vector<float> best(R*R);
    vector<int> bestId(R*R);
    auto project = [&](int view, const Vec3& p, double& u, double& v, double& dep)->bool{
        double sx, sy;
        if(view==0){ dep=D-p.x; sx=p.y; sy=p.z; }
        else if(view==1){ dep=D+p.x; sx=-p.y; sy=p.z; }
        else if(view==2){ dep=D-p.y; sx=-p.x; sy=p.z; }
        else if(view==3){ dep=D+p.y; sx=p.x; sy=p.z; }
        else if(view==4){ dep=D-p.z; sx=p.x; sy=p.y; }
        else { dep=D+p.z; sx=-p.x; sy=p.y; }
        if(dep<=1e-9) return false;
        u=FOC*sx/dep + 512.0;
        v=FOC*sy/dep + 512.0;
        return u>=0.0 && u<1024.0 && v>=0.0 && v<1024.0;
    };
    for(int view=0; view<6; ++view){
        fill(best.begin(), best.end(), 1e30f);
        fill(bestId.begin(), bestId.end(), -1);
        for(int i=0;i<N;i++){
            double u,v,d;
            if(!project(view,P[i],u,v,d)) continue;
            int ix=(int)(u * R / 1024.0);
            int iy=(int)(v * R / 1024.0);
            if(ix<0||ix>=R||iy<0||iy>=R) continue;
            int id=iy*R+ix;
            if(d < best[id]){ best[id]=(float)d; bestId[id]=i; }
        }
        for(int id: bestId) if(id>=0 && marks[id] < 60000) marks[id]++;
    }
    int minx=0,maxx=0,miny=0,maxy=0,minz=0,maxz=0;
    for(int i=1;i<N;i++){
        if(P[i].x<P[minx].x) minx=i; if(P[i].x>P[maxx].x) maxx=i;
        if(P[i].y<P[miny].y) miny=i; if(P[i].y>P[maxy].y) maxy=i;
        if(P[i].z<P[minz].z) minz=i; if(P[i].z>P[maxz].z) maxz=i;
    }
    int ex[6]={minx,maxx,miny,maxy,minz,maxz};
    for(int k=0;k<6;k++) marks[ex[k]] = (unsigned short)min<int>(65535, marks[ex[k]] + 40);
}

static void initialize_quadrics(){
    Q.assign(N, Quadric());
    radiusRep.assign(N, 0.0f);
    clusterSz.assign(N, 1);
    vAlive.assign(N, 1);
    pointWeight.assign(N, 0.0f);
    aliveVertices=N; aliveFaces=M;

    vector<Vec3> nsum(N, Vec3(0,0,0));
    vector<double> asum(N, 0.0), curv(N, 0.0);
    double totalArea=0.0;
    for(int i=0;i<M;i++){
        compute_face_geom(i);
        totalArea += faces[i].area;
        Vec3 n(faces[i].nx, faces[i].ny, faces[i].nz);
        for(int k=0;k<3;k++){
            int v=faces[i].v[k];
            nsum[v]=nsum[v] + n*(faces[i].area);
            asum[v]+=faces[i].area;
        }
    }
    avgAreaGlobal = (M>0 && totalArea>0.0) ? totalArea / (double)M : 1.0;
    vector<Vec3> meanN(N);
    for(int i=0;i<N;i++){
        double l=sqrt(norm2d(nsum[i]));
        if(l>1e-30) meanN[i]=nsum[i]*(float)(1.0/l);
        else meanN[i]=Vec3(0,0,1);
    }
    for(int i=0;i<M;i++){
        Vec3 n(faces[i].nx,faces[i].ny,faces[i].nz);
        for(int k=0;k<3;k++){
            int v=faces[i].v[k];
            double d=fabs(dotd(n,meanN[v]));
            if(d>1) d=1;
            curv[v] += (double)faces[i].area * (1.0-d);
        }
    }
    vector<unsigned short> marks;
    compute_visual_marks(marks);

    for(int i=0;i<N;i++){
        double c = (asum[i]>0) ? curv[i]/asum[i] : 0.0;
        double w = 0.00025 + 0.0040 * min(1.0, c*4.0) + 0.0025 * (double)marks[i];
        if(marks[i] >= 40) w += 0.08;
        pointWeight[i]=(float)w;
        Q[i].addPoint(P[i], w);
    }
    for(int i=0;i<M;i++){
        Face& f=faces[i];
        double areaRatio = f.area / max(1e-30, avgAreaGlobal);
        double viewFactor = 0.70 + fabs((double)f.nx) + fabs((double)f.ny) + fabs((double)f.nz);
        double w = max(0.025, areaRatio) * viewFactor;
        f.w = (float)w;
        Vec3 a=P[f.v[0]];
        double d = -((double)f.nx*a.x + (double)f.ny*a.y + (double)f.nz*a.z);
        Quadric qf;
        qf.addPlane(f.nx,f.ny,f.nz,d,w);
        Q[f.v[0]].add(qf); Q[f.v[1]].add(qf); Q[f.v[2]].add(qf);
    }
}

struct Candidate {
    bool ok;
    int keep, del;
    double cost;
};

static double direction_normal_penalty(int keep, int del, bool& ok){
    ok=true;
    double pen=0.0;
    const Vec3 pk=P[keep];
    for(int fid: inc[del]){
        Face& f=faces[fid];
        if(!f.alive) continue;
        if(!face_has_vertex(f,del)) continue;
        if(face_has_vertex(f,keep)) continue;
        int a=f.v[0], b=f.v[1], c=f.v[2];
        if(a==del) a=keep; else if(b==del) b=keep; else c=keep;
        if(a==b || b==c || c==a){ ok=false; return 1e100; }
        Vec3 cr=crossf(P[b]-P[a], P[c]-P[a]);
        double len2=norm2d(cr);
        if(len2 <= 1e-28){ ok=false; return 1e100; }
        double inv=1.0/sqrt(len2);
        double nx=cr.x*inv, ny=cr.y*inv, nz=cr.z*inv;
        double d=nx*f.nx + ny*f.ny + nz*f.nz;
        if(d < -0.02){ ok=false; return 1e100; }
        if(d < 1.0){ double t=1.0-d; pen += (double)f.w * t*t; }
    }
    (void)pk;
    return pen;
}

static Candidate compute_candidate(int u, int v){
    Candidate ret{false,-1,-1,1e100};
    if(u<0||v<0||u==v||u>=N||v>=N||!vAlive[u]||!vAlive[v]) return ret;
    Quadric qs=Q[u]; qs.add(Q[v]);
    double d=distd(P[u],P[v]);
    auto try_dir = [&](int keep, int del){
        if(max((double)radiusRep[keep], d + (double)radiusRep[del]) > hausTol) return;
        bool ok=true;
        double pen=direction_normal_penalty(keep,del,ok);
        if(!ok) return;
        double c=qs.eval(P[keep]) + 0.00020 * pen;
        c += 0.00002 * d*d * (double)(clusterSz[keep] + clusterSz[del]);
        if(c < ret.cost){ ret={true,keep,del,c}; }
    };
    try_dir(u,v);
    try_dir(v,u);
    return ret;
}

struct HeapItem {
    float cost;
    int u, v;
    bool operator<(const HeapItem& o) const { return cost > o.cost; }
};
static priority_queue<HeapItem> heapq;

static inline uint64_t edge_key_int(int a, int b){ if(a>b) swap(a,b); return (uint64_t)(uint32_t)a<<32 | (uint32_t)b; }

static void push_edge(int a, int b){
    if(a==b || !vAlive[a] || !vAlive[b]) return;
    Candidate c=compute_candidate(a,b);
    if(c.ok){
        int u=min(a,b), v=max(a,b);
        double cc=c.cost;
        if(!(cc>=0) || !isfinite(cc)) cc=0;
        if(cc > 3.3e38) cc = 3.3e38;
        heapq.push(HeapItem{(float)cc,u,v});
    }
}

static void build_initial_heap(){
    vector<uint64_t> edges;
    edges.reserve((size_t)M*3);
    for(int i=0;i<M;i++){
        int a=faces[i].v[0], b=faces[i].v[1], c=faces[i].v[2];
        edges.push_back(edge_key_int(a,b));
        edges.push_back(edge_key_int(b,c));
        edges.push_back(edge_key_int(c,a));
    }
    sort(edges.begin(), edges.end());
    edges.erase(unique(edges.begin(), edges.end()), edges.end());
    for(uint64_t k: edges){
        int a=(int)(k>>32), b=(int)(uint32_t)k;
        push_edge(a,b);
    }
}

static vector<int> markArr;
static int stampBase=1;
static inline void reset_stamps_if_needed(){
    if(stampBase > 1000000000){ fill(markArr.begin(), markArr.end(), 0); stampBase=1; }
}

static bool has_face_with_pair(int v, int a, int b, int exclude){
    for(int fid: inc[v]){
        if(fid==exclude) continue;
        Face& f=faces[fid];
        if(!f.alive) continue;
        if(face_has_vertex(f,v) && face_has_vertex(f,a) && face_has_vertex(f,b)) return true;
    }
    return false;
}

static bool link_condition_ok(int u, int v){
    if(aliveVertices <= 4) return false;
    reset_stamps_if_needed();
    int st=stampBase; stampBase+=3;
    int edgeFacesCnt=0;
    int opp[4]={-1,-1,-1,-1};
    for(int fid: inc[u]){
        Face& f=faces[fid];
        if(!f.alive || !face_has_vertex(f,u)) continue;
        for(int k=0;k<3;k++) if(f.v[k]!=u) markArr[f.v[k]]=st;
        if(face_has_vertex(f,v)){
            if(edgeFacesCnt<4) opp[edgeFacesCnt]=face_other_of_edge(f,u,v);
            edgeFacesCnt++;
        }
    }
    if(edgeFacesCnt!=2 || opp[0]<0 || opp[1]<0 || opp[0]==opp[1]) return false;
    markArr[opp[0]]=st+1;
    markArr[opp[1]]=st+1;
    int common=0;
    bool bad=false;
    for(int fid: inc[v]){
        Face& f=faces[fid];
        if(!f.alive || !face_has_vertex(f,v)) continue;
        for(int k=0;k<3;k++){
            int w=f.v[k];
            if(w==v) continue;
            if(markArr[w]==st || markArr[w]==st+1){
                if(markArr[w]!=st+1) bad=true;
                common++;
                markArr[w]=st+2;
            }
        }
    }
    if(bad || common!=2) return false;
    int fuv0=-1, fuv1=-1;
    for(int fid: inc[u]){
        Face& f=faces[fid]; if(!f.alive) continue;
        if(face_has_edge(f,u,v)){ if(fuv0<0) fuv0=fid; else fuv1=fid; }
    }
    bool hu=has_face_with_pair(u,opp[0],opp[1],-1);
    bool hv=has_face_with_pair(v,opp[0],opp[1],-1);
    (void)fuv0; (void)fuv1;
    if(hu && hv) return false;
    return true;
}

static void compact_incident(int v){
    vector<int>& lst=inc[v];
    size_t w=0;
    for(size_t i=0;i<lst.size();++i){
        int fid=lst[i];
        if(fid>=0 && fid<M && faces[fid].alive && face_has_vertex(faces[fid],v)) lst[w++]=fid;
    }
    lst.resize(w);
}

static bool g_push_updates=true;
static bool collapse_edge(int keep, int del){
    if(!vAlive[keep] || !vAlive[del] || keep==del) return false;
    if(max((double)radiusRep[keep], distd(P[keep],P[del]) + (double)radiusRep[del]) > hausTol) return false;
    bool ok=true;
    (void)direction_normal_penalty(keep,del,ok);
    if(!ok) return false;

    vector<int> affected;
    affected.reserve(32);
    affected.push_back(keep);
    for(int fid: inc[del]){
        Face& f=faces[fid];
        if(!f.alive || !face_has_vertex(f,del)) continue;
        int a=f.v[0], b=f.v[1], c=f.v[2];
        if(face_has_vertex(f,keep)){
            f.alive=0;
            aliveFaces--;
            for(int k=0;k<3;k++) if(f.v[k]!=keep && f.v[k]!=del) affected.push_back(f.v[k]);
        } else {
            if(f.v[0]==del) f.v[0]=keep;
            else if(f.v[1]==del) f.v[1]=keep;
            else f.v[2]=keep;
            compute_face_geom(fid);
            inc[keep].push_back(fid);
            affected.push_back(a); affected.push_back(b); affected.push_back(c);
        }
    }
    vAlive[del]=0;
    aliveVertices--;
    radiusRep[keep]=(float)max((double)radiusRep[keep], distd(P[keep],P[del]) + (double)radiusRep[del]);
    clusterSz[keep]+=clusterSz[del];
    Q[keep].add(Q[del]);
    compact_incident(keep);

    sort(affected.begin(), affected.end());
    affected.erase(unique(affected.begin(), affected.end()), affected.end());
    if(g_push_updates){
        int pushed=0;
        for(int a: affected){
            if(a<0 || a>=N || !vAlive[a]) continue;
            if((int)inc[a].size() > 800 && a!=keep) continue;
            int scans=0;
            for(int fid: inc[a]){
                Face& f=faces[fid];
                if(!f.alive || !face_has_vertex(f,a)) continue;
                push_edge(f.v[0],f.v[1]);
                push_edge(f.v[1],f.v[2]);
                push_edge(f.v[2],f.v[0]);
                if(++scans > 1200) break;
            }
            if((pushed+=scans) > 4000) break;
        }
    }
    return true;
}

struct MeshOut {
    vector<Vec3> verts;
    vector<Tri> tris;
    double approxScore=0;
};

static MeshOut make_snapshot(){
    MeshOut mo;
    mo.verts.reserve(aliveVertices);
    vector<int> mp(N,-1);
    for(int i=0;i<N;i++) if(vAlive[i]){ mp[i]=(int)mo.verts.size(); mo.verts.push_back(P[i]); }
    mo.tris.reserve(aliveFaces);
    for(int i=0;i<M;i++) if(faces[i].alive){
        int a=mp[faces[i].v[0]], b=mp[faces[i].v[1]], c=mp[faces[i].v[2]];
        if(a>=0&&b>=0&&c>=0&&a!=b&&b!=c&&c!=a) mo.tris.push_back({a,b,c});
    }
    return mo;
}

static vector<int> target_counts(){
    vector<double> ratios;
    if(N < 20){
        return vector<int>{max(4,N-1)};
    } else if(N < 200){
        ratios={0.82,0.70,0.60,0.50};
    } else if(N <= 6500){
        ratios={0.320,0.260,0.220,0.175,0.140,0.115,0.095,0.080,0.068,0.058};
    } else if(N <= 70000){
        ratios={0.280,0.230,0.180,0.148,0.122,0.102,0.088,0.077,0.067,0.058,0.051};
    } else if(N <= 500000){
        ratios={0.240,0.195,0.160,0.132,0.110,0.094,0.082,0.072,0.063,0.055,0.049};
    } else {
        ratios={0.220,0.180,0.145,0.122,0.104,0.090,0.079,0.069,0.060,0.053,0.047};
    }
    vector<int> res;
    int last=N+1;
    for(double r: ratios){
        int t=(int)ceil(N*r);
        t=max(4,min(N-1,t));
        if(t<last){ res.push_back(t); last=t; }
    }
    if(res.empty()) res.push_back(max(4,N-1));
    return res;
}

struct RenderSet {
    int R=0;
    vector<float> img;
    vector<unsigned char> mask;
};

static inline void view_project_low(int view, const Vec3& p, int R, float& u, float& v, float& dep){
    const float D=2.5f;
    float sx, sy;
    if(view==0){ dep=D-p.x; sx=p.y; sy=p.z; }
    else if(view==1){ dep=D+p.x; sx=-p.y; sy=p.z; }
    else if(view==2){ dep=D-p.y; sx=-p.x; sy=p.z; }
    else if(view==3){ dep=D+p.y; sx=p.x; sy=p.z; }
    else if(view==4){ dep=D-p.z; sx=p.x; sy=p.y; }
    else { dep=D+p.z; sx=-p.x; sy=p.y; }
    float f = 800.0f * (float)R / 1024.0f;
    float c = 0.5f * (float)R;
    u = f * sx / dep + c;
    v = f * sy / dep + c;
}

static RenderSet render_mesh_low(const vector<Vec3>& verts, const vector<Tri>& tris, int R){
    RenderSet rs; rs.R=R;
    int S=R*R;
    rs.img.assign((size_t)6*4*S, 127.5f);
    rs.mask.assign((size_t)6*S, 0);
    vector<float> zbuf(S);
    for(int view=0; view<6; ++view){
        fill(zbuf.begin(), zbuf.end(), 1e30f);
        float* base = rs.img.data() + (size_t)view*4*S;
        fill(base + 3*S, base + 4*S, 255.0f);
        for(const Tri& t: tris){
            const Vec3 &A=verts[t.a], &B=verts[t.b], &C=verts[t.c];
            Vec3 cr=crossf(B-A,C-A);
            double len2=norm2d(cr);
            if(len2<=1e-30) continue;
            double invn=1.0/sqrt(len2);
            float nr=(float)((cr.x*invn + 1.0)*127.5);
            float ng=(float)((cr.y*invn + 1.0)*127.5);
            float nb=(float)((cr.z*invn + 1.0)*127.5);
            float u0,v0,d0,u1,v1,d1,u2,v2,d2;
            view_project_low(view,A,R,u0,v0,d0);
            view_project_low(view,B,R,u1,v1,d1);
            view_project_low(view,C,R,u2,v2,d2);
            if(d0<=0||d1<=0||d2<=0) continue;
            float minx=floorf(min(u0,min(u1,u2))-0.5f), maxx=ceilf(max(u0,max(u1,u2))+0.5f);
            float miny=floorf(min(v0,min(v1,v2))-0.5f), maxy=ceilf(max(v0,max(v1,v2))+0.5f);
            int ix0=max(0,(int)minx), ix1=min(R-1,(int)maxx);
            int iy0=max(0,(int)miny), iy1=min(R-1,(int)maxy);
            if(ix0>ix1||iy0>iy1) continue;
            float denom=(u1-u0)*(v2-v0)-(v1-v0)*(u2-u0);
            if(fabs(denom)<1e-12f) continue;
            float invDen=1.0f/denom;
            for(int y=iy0;y<=iy1;y++){
                float py=y+0.5f;
                for(int x=ix0;x<=ix1;x++){
                    float px=x+0.5f;
                    float w1=((px-u0)*(v2-v0)-(py-v0)*(u2-u0))*invDen;
                    float w2=((u1-u0)*(py-v0)-(v1-v0)*(px-u0))*invDen;
                    float w0=1.0f-w1-w2;
                    if(w0>=-1e-5f && w1>=-1e-5f && w2>=-1e-5f){
                        float invz=w0/d0 + w1/d1 + w2/d2;
                        if(invz<=0) continue;
                        float dep=1.0f/invz;
                        int id=y*R+x;
                        if(dep < zbuf[id]){
                            zbuf[id]=dep;
                            base[id]=nr; base[S+id]=ng; base[2*S+id]=nb; base[3*S+id]=dep;
                        }
                    }
                }
            }
        }
        for(int i=0;i<S;i++) if(zbuf[i]<1e29f) rs.mask[(size_t)view*S+i]=1;
    }
    return rs;
}

static RenderSet render_original_low(int R){
    vector<Tri> tris=origFaces;
    return render_mesh_low(P,tris,R);
}

static double ssim_channel(const float* X, const float* Y, const unsigned char* mx, const unsigned char* my, int R){
    int S=R*R;
    int W=R+1;
    vector<double> ix((size_t)W*W,0), iy((size_t)W*W,0), ix2((size_t)W*W,0), iy2((size_t)W*W,0), ixy((size_t)W*W,0);
    for(int r=0;r<R;r++){
        double sx=0,sy=0,sx2=0,sy2=0,sxy=0;
        for(int c=0;c<R;c++){
            int id=r*R+c;
            sx += X[id]; sy += Y[id]; sx2 += (double)X[id]*X[id]; sy2 += (double)Y[id]*Y[id]; sxy += (double)X[id]*Y[id];
            size_t o=(size_t)(r+1)*W+(c+1), up=(size_t)r*W+(c+1);
            ix[o]=ix[up]+sx; iy[o]=iy[up]+sy; ix2[o]=ix2[up]+sx2; iy2[o]=iy2[up]+sy2; ixy[o]=ixy[up]+sxy;
        }
    }
    auto rect=[&](const vector<double>& I,int y0,int x0,int y1,int x1){
        return I[(size_t)y1*W+x1]-I[(size_t)y0*W+x1]-I[(size_t)y1*W+x0]+I[(size_t)y0*W+x0];
    };
    const double C1=6.5025, C2=58.5225;
    double sum=0.0; int cnt=0;
    int rad = (R>=192?5:4);
    for(int y=0;y<R;y++) for(int x=0;x<R;x++){
        int id=y*R+x;
        if(!mx[id] && !my[id]) continue;
        int y0=max(0,y-rad), y1=min(R,y+rad+1), x0=max(0,x-rad), x1=min(R,x+rad+1);
        double n=(double)(y1-y0)*(x1-x0);
        double sx=rect(ix,y0,x0,y1,x1), sy=rect(iy,y0,x0,y1,x1);
        double sx2=rect(ix2,y0,x0,y1,x1), sy2=rect(iy2,y0,x0,y1,x1), sxy=rect(ixy,y0,x0,y1,x1);
        double mux=sx/n, muy=sy/n;
        double vx=max(0.0, sx2/n - mux*mux);
        double vy=max(0.0, sy2/n - muy*muy);
        double cov=sxy/n - mux*muy;
        double num=(2*mux*muy + C1)*(2*cov + C2);
        double den=(mux*mux + muy*muy + C1)*(vx+vy+C2);
        double s=(den!=0)?num/den:1.0;
        if(s<0) s=0; if(s>1) s=1;
        sum += s; cnt++;
    }
    return cnt?sum/cnt:1.0;
}

static double final_ssim_low(const RenderSet& A, const RenderSet& B){
    int R=A.R, S=R*R;
    double total=0.0;
    for(int view=0; view<6; ++view){
        const float* ax=A.img.data()+(size_t)view*4*S;
        const float* bx=B.img.data()+(size_t)view*4*S;
        const unsigned char* am=A.mask.data()+(size_t)view*S;
        const unsigned char* bm=B.mask.data()+(size_t)view*S;
        double ns=0;
        for(int ch=0; ch<3; ++ch) ns += ssim_channel(ax+ch*S,bx+ch*S,am,bm,R);
        ns/=3.0;
        double ds=ssim_channel(ax+3*S,bx+3*S,am,bm,R);
        total += 0.5*(ns+ds);
    }
    return total/6.0;
}

static bool validate_mesh(const MeshOut& mo){
    int V=(int)mo.verts.size();
    if(V<1 || V>N) return false;
    if(mo.tris.empty()) return false;
    vector<unsigned char> used(V,0);
    vector<uint64_t> edges;
    edges.reserve((size_t)mo.tris.size()*3);
    for(const Tri& t: mo.tris){
        if(t.a<0||t.a>=V||t.b<0||t.b>=V||t.c<0||t.c>=V) return false;
        if(t.a==t.b||t.b==t.c||t.c==t.a) return false;
        Vec3 cr=crossf(mo.verts[t.b]-mo.verts[t.a], mo.verts[t.c]-mo.verts[t.a]);
        if(norm2d(cr) <= 1e-28) return false;
        used[t.a]=used[t.b]=used[t.c]=1;
        edges.push_back(edge_key_int(t.a,t.b)); edges.push_back(edge_key_int(t.b,t.c)); edges.push_back(edge_key_int(t.c,t.a));
    }
    for(int i=0;i<V;i++) if(!used[i]) return false;
    sort(edges.begin(), edges.end());
    for(size_t i=0;i<edges.size();){
        size_t j=i+1; while(j<edges.size() && edges[j]==edges[i]) j++;
        if(j-i != 2) return false;
        i=j;
    }
    return true;
}

static void output_mesh(const MeshOut& mo){
    string out;
    out.reserve((size_t)mo.verts.size()*42 + (size_t)mo.tris.size()*26 + 64);
    char line[128];
    out.append(line, snprintf(line,sizeof(line), "%d %d\n", (int)mo.verts.size(), (int)mo.tris.size()));
    for(const Vec3& p: mo.verts){
        out.append(line, snprintf(line,sizeof(line), "v %.10g %.10g %.10g\n", p.x,p.y,p.z));
    }
    for(const Tri& t: mo.tris){
        out.append(line, snprintf(line,sizeof(line), "f %d %d %d\n", t.a+1,t.b+1,t.c+1));
    }
    fwrite(out.data(),1,out.size(),stdout);
}

static void output_original(){
    string out;
    out.reserve((size_t)N*42 + (size_t)M*28 + 64);
    char line[128];
    out.append(line, snprintf(line,sizeof(line), "%d %d\n", N, M));
    for(const Vec3& p: P) out.append(line, snprintf(line,sizeof(line), "v %.10g %.10g %.10g\n", p.x,p.y,p.z));
    for(const Tri& t: origFaces) out.append(line, snprintf(line,sizeof(line), "f %d %d %d\n", t.a+1,t.b+1,t.c+1));
    fwrite(out.data(),1,out.size(),stdout);
}

struct EdgeCost { float cost; int u,v; bool operator<(const EdgeCost& o) const { return cost<o.cost; } };
static vector<EdgeCost> build_edge_costs_current(){
    vector<uint64_t> keys; keys.reserve((size_t)aliveFaces*3);
    for(int i=0;i<M;i++) if(faces[i].alive){
        int a=faces[i].v[0], b=faces[i].v[1], c=faces[i].v[2];
        keys.push_back(edge_key_int(a,b)); keys.push_back(edge_key_int(b,c)); keys.push_back(edge_key_int(c,a));
    }
    sort(keys.begin(), keys.end()); keys.erase(unique(keys.begin(), keys.end()), keys.end());
    vector<EdgeCost> ec; ec.reserve(keys.size());
    for(uint64_t k: keys){
        int a=(int)(k>>32), b=(int)(uint32_t)k;
        Candidate c=compute_candidate(a,b);
        if(c.ok){
            double cc=c.cost;
            if(!isfinite(cc)||cc<0) cc=0;
            if(cc>3.3e38) cc=3.3e38;
            ec.push_back({(float)cc,a,b});
        }
    }
    sort(ec.begin(), ec.end());
    return ec;
}

static vector<MeshOut> simplify_batch(const vector<int>& targets){
    vector<MeshOut> candidates;
    int nextIdx=0;
    int targetMin=targets.back();
    g_push_updates=false;
    int pass=0;
    long long totalColl=0;
    vector<int> touched(N,0);
    while(aliveVertices>targetMin && pass<30){
        auto edges=build_edge_costs_current();
        if(edges.empty()) break;
        int before=aliveVertices;
        long long coll=0;
        int st=pass+1;
        for(const auto& e: edges){
            if(aliveVertices<=targetMin) break;
            int u=e.u,v=e.v;
            if(!vAlive[u]||!vAlive[v]||u==v) continue;
            if(touched[u]==st || touched[v]==st) continue;
            Candidate c=compute_candidate(u,v);
            if(!c.ok) continue;
            if(!link_condition_ok(u,v)) continue;
            if(collapse_edge(c.keep,c.del)){
                coll++;
                totalColl++;
                touched[c.keep]=st;
                touched[c.del]=st;
                while(nextIdx<(int)targets.size() && aliveVertices <= targets[nextIdx]){
                    MeshOut snap=make_snapshot();
                    if(!snap.tris.empty()) candidates.push_back(std::move(snap));
                    nextIdx++;
                }
            }
        }
        pass++;
        if(coll==0) break;
        if(aliveVertices==before) break;
    }
    (void)totalColl;
    g_push_updates=true;
    return candidates;
}

static void simplify_and_output(){
    if(N<=4){ output_original(); return; }
    initialize_quadrics();
    markArr.assign(N,0);
    vector<int> targets=target_counts();
    vector<MeshOut> candidates = simplify_batch(targets);
    if(candidates.empty()) candidates.push_back(make_snapshot());

    int chosen=-1;
    if(N>=200 && !candidates.empty()){
        sort(candidates.begin(), candidates.end(), [](const MeshOut& a, const MeshOut& b){ return a.verts.size() < b.verts.size(); });
        int R = (N<=7000?256:(N<=80000?256:(N<=500000?176:144)));
        RenderSet ref = render_original_low(R);
        double threshold = (R>=256 ? 0.905 : (R>=176 ? 0.906 : 0.908));
        for(int i=0;i<(int)candidates.size();++i){
            RenderSet rr = render_mesh_low(candidates[i].verts, candidates[i].tris, R);
            candidates[i].approxScore = final_ssim_low(ref, rr);
            if(candidates[i].approxScore >= threshold){ chosen=i; break; }
        }
        if(chosen<0) chosen=(int)candidates.size()-1;
    } else {
        chosen=(int)candidates.size()-1;
    }

    sort(candidates.begin(), candidates.end(), [](const MeshOut& a, const MeshOut& b){ return a.verts.size() < b.verts.size(); });
    if(chosen<0 || chosen>=(int)candidates.size()) chosen=(int)candidates.size()-1;
    for(int i=chosen;i<(int)candidates.size();++i){
        if(validate_mesh(candidates[i])){ output_mesh(candidates[i]); return; }
    }
    for(int i=chosen-1;i>=0;--i){
        if(validate_mesh(candidates[i])){ output_mesh(candidates[i]); return; }
    }
    output_original();
}

int main(){
    load_input();
    simplify_and_output();
    return 0;
}
