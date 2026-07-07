#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <numeric>
#include <unordered_map>
#include <utility>
#include <vector>
using namespace std;

/*
  workerI_case7_fast_safe.cpp
  CASE7 / HIGH-N FAST-SAFE COMPRESSION

  Activation policy:
    * N < HIGH_N_TRIGGER: exact identity output. This file is intended as a
      high-N drop-in/worker; it deliberately avoids changing small and mid cases.
    * N >= HIGH_N_TRIGGER: first try strict row-major UV/grid detection
      (open grid, seam grid, tube, torus). If validated, remesh by sampled bands.
    * If grid fails, try very strict axis-aligned ellipsoid/sphere detector.
    * Otherwise run a bounded topology-preserving endpoint edge-collapse fallback.

  Safety guards:
    * All generic collapses require a 2-face manifold edge, the link condition,
      local normal preservation, duplicate-face prevention, and a cluster radius cap.
    * Every heavy loop checks a wall-clock guard; on timeout it outputs the best
      currently valid mesh instead of continuing.
    * Output is buffered and numeric precision is capped to keep files small.
*/

static constexpr int HIGH_N_TRIGGER = 200000;
static constexpr double TIME_LIMIT_SEC = 18.20;
static constexpr double GRID_STRICT_MATCH = 0.940;
static constexpr int MAX_GRID_VALIDATE_FACES = 220000;
static constexpr int MAX_DELTA_SAMPLE_FACES = 260000;

struct Vec3 {
    double x=0, y=0, z=0;
    Vec3() = default;
    Vec3(double X,double Y,double Z):x(X),y(Y),z(Z){}
};
static inline Vec3 operator+(const Vec3&a,const Vec3&b){return {a.x+b.x,a.y+b.y,a.z+b.z};}
static inline Vec3 operator-(const Vec3&a,const Vec3&b){return {a.x-b.x,a.y-b.y,a.z-b.z};}
static inline Vec3 operator*(const Vec3&a,double s){return {a.x*s,a.y*s,a.z*s};}
static inline Vec3 operator/(const Vec3&a,double s){return {a.x/s,a.y/s,a.z/s};}
static inline double dotv(const Vec3&a,const Vec3&b){return a.x*b.x+a.y*b.y+a.z*b.z;}
static inline Vec3 crossv(const Vec3&a,const Vec3&b){return {a.y*b.z-a.z*b.y,a.z*b.x-a.x*b.z,a.x*b.y-a.y*b.x};}
static inline double norm2v(const Vec3&a){return dotv(a,a);} 
static inline double normv(const Vec3&a){return sqrt(norm2v(a));}
static inline double dist2v(const Vec3&a,const Vec3&b){return norm2v(a-b);} 

struct Face { int a=0,b=0,c=0; unsigned char active=1; };

struct FastInput {
    vector<char> buf;
    char* p=nullptr;
    FastInput(){
        buf.reserve(1u<<26);
        char tmp[1<<16];
        size_t n=0;
        while((n=fread(tmp,1,sizeof(tmp),stdin))>0) buf.insert(buf.end(),tmp,tmp+n);
        buf.push_back('\0');
        p=buf.data();
    }
    inline void skip_ws(){ while(*p==' '||*p=='\n'||*p=='\r'||*p=='\t') ++p; }
    int next_int(){ skip_ws(); int s=1; if(*p=='-'){s=-1;++p;} int x=0; while(*p>='0'&&*p<='9'){x=x*10+(*p-'0');++p;} return x*s; }
    long next_long(){ skip_ws(); return strtol(p,&p,10); }
    double next_double(){ skip_ws(); char* e=nullptr; double x=strtod(p,&e); p=e; return x; }
    char next_char(){ skip_ws(); return *p++; }
};

static int N=0, M=0;
static vector<Vec3> P;
static vector<Face> F;
static Vec3 bboxMin, bboxMax, bboxCtr;
static double bboxDiag = 1.0;
static chrono::steady_clock::time_point startTime;

static inline bool time_ok(double marginSec=TIME_LIMIT_SEC){
    double t=chrono::duration<double>(chrono::steady_clock::now()-startTime).count();
    return t < marginSec;
}
static inline uint64_t edge_key(int a,int b){ if(a>b) swap(a,b); return (uint64_t)(uint32_t)a<<32 | (uint32_t)b; }
static inline int key_a(uint64_t k){ return (int)(k>>32); }
static inline int key_b(uint64_t k){ return (int)(k & 0xffffffffu); }
static inline array<int,3> sorted_tri(int a,int b,int c){ array<int,3> s{a,b,c}; sort(s.begin(),s.end()); return s; }
static inline bool face_contains(const Face&f,int v){ return f.a==v||f.b==v||f.c==v; }
static inline Vec3 tri_normal_idx(int a,int b,int c){ return crossv(P[b]-P[a],P[c]-P[a]); }

static void load_mesh(){
    FastInput in;
    N=(int)in.next_long();
    M=(int)in.next_long();
    P.resize(N);
    F.resize(M);
    bboxMin={1e100,1e100,1e100}; bboxMax={-1e100,-1e100,-1e100};
    for(int i=0;i<N;i++){
        (void)in.next_char();
        P[i].x=in.next_double(); P[i].y=in.next_double(); P[i].z=in.next_double();
        bboxMin.x=min(bboxMin.x,P[i].x); bboxMin.y=min(bboxMin.y,P[i].y); bboxMin.z=min(bboxMin.z,P[i].z);
        bboxMax.x=max(bboxMax.x,P[i].x); bboxMax.y=max(bboxMax.y,P[i].y); bboxMax.z=max(bboxMax.z,P[i].z);
    }
    bboxCtr=(bboxMin+bboxMax)*0.5;
    bboxDiag=normv(bboxMax-bboxMin); if(!(bboxDiag>0)) bboxDiag=1.0;
    for(int i=0;i<M;i++){
        (void)in.next_char();
        int a=(int)in.next_long()-1, b=(int)in.next_long()-1, c=(int)in.next_long()-1;
        F[i].a=a; F[i].b=b; F[i].c=c; F[i].active=1;
    }
}

static void write_mesh_vectors(const vector<Vec3>& V, const vector<array<int,3>>& T){
    static char outbuf[1<<20];
    setvbuf(stdout,outbuf,_IOFBF,sizeof(outbuf));
    printf("%d %d\n", (int)V.size(), (int)T.size());
    for(const Vec3& p: V) printf("v %.12g %.12g %.12g\n", p.x,p.y,p.z);
    for(const auto& t: T) printf("f %d %d %d\n", t[0]+1,t[1]+1,t[2]+1);
}

static void write_identity(){
    static char outbuf[1<<20];
    setvbuf(stdout,outbuf,_IOFBF,sizeof(outbuf));
    printf("%d %d\n", N, M);
    for(int i=0;i<N;i++) printf("v %.12g %.12g %.12g\n", P[i].x,P[i].y,P[i].z);
    for(int i=0;i<M;i++) printf("f %d %d %d\n", F[i].a+1,F[i].b+1,F[i].c+1);
}

// ------------------------- strict index-grid / UV remeshing -------------------------

struct GridInfo {
    bool ok=false;
    int W=0,H=0;
    bool wrapR=false, wrapC=false;
    double match=0.0;
    double faceErr=1.0;
    array<array<int,3>,2> patterns{{array<int,3>{0,1,2}, array<int,3>{1,3,2}}};
};

static inline int grid_idx(int r,int c,int H,int W,bool wrapR,bool wrapC){
    if(wrapR){ r%=H; if(r<0) r+=H; }
    else if(r<0||r>=H) return -1;
    if(wrapC){ c%=W; if(c<0) c+=W; }
    else if(c<0||c>=W) return -1;
    return r*W+c;
}

static bool decode_grid_face(int a,int b,int c,int H,int W,bool wrapR,bool wrapC,array<int,3>* orderedPattern=nullptr){
    int vs[3]={a,b,c};
    int rr[3]={a/W,b/W,c/W};
    int cc[3]={a%W,b%W,c%W};
    for(int t=0;t<3;t++){
        for(int dr=0;dr<=1;dr++) for(int dc=0;dc<=1;dc++){
            int br=rr[t]-dr, bc=cc[t]-dc;
            if(!wrapR && (br<0 || br+1>=H)) continue;
            if(!wrapC && (bc<0 || bc+1>=W)) continue;
            int corner[4];
            corner[0]=grid_idx(br,bc,H,W,wrapR,wrapC);
            corner[1]=grid_idx(br,bc+1,H,W,wrapR,wrapC);
            corner[2]=grid_idx(br+1,bc,H,W,wrapR,wrapC);
            corner[3]=grid_idx(br+1,bc+1,H,W,wrapR,wrapC);
            int code[3]; bool used[4]={false,false,false,false}; bool good=true;
            for(int i=0;i<3;i++){
                code[i]=-1;
                for(int k=0;k<4;k++) if(vs[i]==corner[k]) { code[i]=k; break; }
                if(code[i]<0 || used[code[i]]) { good=false; break; }
                used[code[i]]=true;
            }
            if(!good) continue;
            int r0=(code[0]>=2), r1=(code[1]>=2), r2=(code[2]>=2);
            int c0=(code[0]&1), c1=(code[1]&1), c2=(code[2]&1);
            if((r0==r1 && r1==r2) || (c0==c1 && c1==c2)) continue; // collinear in param cell
            if(orderedPattern) *orderedPattern={code[0],code[1],code[2]};
            return true;
        }
    }
    return false;
}

static bool validate_grid_candidate(int W,bool wrapR,bool wrapC,GridInfo& out){
    if(W<3 || N%W!=0) return false;
    int H=N/W;
    if(H<3) return false;
    long long cellsR = wrapR ? H : H-1;
    long long cellsC = wrapC ? W : W-1;
    if(cellsR<=0 || cellsC<=0) return false;
    long long expectedFaces = 2LL*cellsR*cellsC;
    double faceErr = fabs((double)M - (double)expectedFaces) / max(1.0, (double)M);
    // Strict enough not to fire on arbitrary indexed meshes, but tolerant of a few cap/seam leftovers.
    if(faceErr > 0.075) return false;

    int step=max(1, M / MAX_GRID_VALIDATE_FACES);
    int total=0, matched=0;
    unordered_map<int,int> patCount;
    patCount.reserve(16);
    for(int i=0;i<M;i+=step){
        const Face& f=F[i];
        if(f.a<0||f.a>=N||f.b<0||f.b>=N||f.c<0||f.c>=N) continue;
        total++;
        array<int,3> ptn;
        if(decode_grid_face(f.a,f.b,f.c,H,W,wrapR,wrapC,&ptn)){
            matched++;
            int enc=ptn[0] | (ptn[1]<<2) | (ptn[2]<<4);
            patCount[enc]++;
        }
        if((total&65535)==0 && !time_ok(5.5)) break;
    }
    if(total<1000 && M>1000) return false;
    double ratio = total ? (double)matched/(double)total : 0.0;
    if(ratio < GRID_STRICT_MATCH) return false;

    vector<pair<int,int>> pats;
    pats.reserve(patCount.size());
    for(auto& kv: patCount) pats.push_back({kv.second, kv.first});
    sort(pats.rbegin(), pats.rend());
    array<array<int,3>,2> best{{array<int,3>{0,1,2}, array<int,3>{1,3,2}}};
    if(pats.size()>=2){
        for(int j=0;j<2;j++){
            int enc=pats[j].second;
            best[j]={enc&3,(enc>>2)&3,(enc>>4)&3};
        }
        int cover=pats[0].first+pats[1].first;
        if(cover < (int)(0.70*max(1,matched))){
            // Too many orientation/pattern variants: keep default orientation rather than trusting noise.
            best={{array<int,3>{0,1,2}, array<int,3>{1,3,2}}};
        }
    }
    out.ok=true; out.W=W; out.H=H; out.wrapR=wrapR; out.wrapC=wrapC; out.match=ratio; out.faceErr=faceErr; out.patterns=best;
    return true;
}

static void add_candidate_w(vector<int>& cand, int w){
    if(w>=3 && w<=N/3 && N%w==0) cand.push_back(w);
}

static vector<int> collect_grid_width_candidates(){
    vector<int> cand;
    cand.reserve(128);
    int faceStep=max(1, M / MAX_DELTA_SAMPLE_FACES);
    vector<int> deltas; deltas.reserve((size_t)min(M/faceStep+1, MAX_DELTA_SAMPLE_FACES)*3);
    for(int i=0;i<M;i+=faceStep){
        int a=F[i].a,b=F[i].b,c=F[i].c;
        if(a<0||b<0||c<0||a>=N||b>=N||c>=N) continue;
        int d1=abs(a-b), d2=abs(b-c), d3=abs(c-a);
        if(d1>1) deltas.push_back(d1);
        if(d2>1) deltas.push_back(d2);
        if(d3>1) deltas.push_back(d3);
    }
    sort(deltas.begin(), deltas.end());
    vector<pair<int,int>> top;
    for(size_t i=0;i<deltas.size();){
        size_t j=i+1; while(j<deltas.size() && deltas[j]==deltas[i]) j++;
        int cnt=(int)(j-i); int d=deltas[i];
        top.push_back({cnt,d});
        i=j;
    }
    sort(top.rbegin(), top.rend());
    for(int i=0;i<(int)top.size() && i<24;i++){
        int d=top[i].second;
        for(int e=-2;e<=2;e++) add_candidate_w(cand,d+e);
    }
    // Face-count algebra for common open/tube grids.
    if((M&1)==0){
        long long half=M/2;
        long long tubeW=(long long)N-half; // M/2=(H-1)*W=N-W
        if(tubeW>0 && tubeW<=N) add_candidate_w(cand,(int)tubeW);
        long long sumHW=(long long)N + 1 - half; // open grid: M/2=N-H-W+1
        if(sumHW>0){
            // If x*(sum-x)=N, x is a dimension. Try around the quadratic roots via divisors below.
            (void)sumHW;
        }
    }
    // Divisors near sqrt help torus M=2N where face count alone gives no W.
    int r=(int)sqrt((double)N);
    for(int d=max(3,r-4096); d<=r+4096 && d<=N/d; ++d){
        if(N%d==0){ add_candidate_w(cand,d); add_candidate_w(cand,N/d); }
    }
    sort(cand.begin(), cand.end());
    cand.erase(unique(cand.begin(), cand.end()), cand.end());
    return cand;
}

static GridInfo detect_index_grid(){
    GridInfo best;
    if(N<HIGH_N_TRIGGER || M<2) return best;
    vector<int> cand=collect_grid_width_candidates();
    double bestScore=-1e100;
    for(int W: cand){
        if(!time_ok(6.5)) break;
        for(int wr=0; wr<2; ++wr) for(int wc=0; wc<2; ++wc){
            GridInfo g;
            if(!validate_grid_candidate(W,wr,wc,g)) continue;
            // Score favors match first, then exact expected face count, then balanced dimensions.
            double aspectPenalty=fabs(log((double)g.W/(double)g.H))*0.002;
            double score=g.match - 0.45*g.faceErr - aspectPenalty;
            if(score>bestScore){ bestScore=score; best=g; }
        }
    }
    return best;
}

static vector<int> make_samples_1d(int n, int want, bool periodic){
    vector<int> s;
    if(n<=0) return s;
    want=max(periodic?3:2, min(want,n));
    s.reserve(want);
    if(periodic){
        int last=-1;
        for(int i=0;i<want;i++){
            int id=(int)((long long)i*n / want);
            if(id!=last) { s.push_back(id); last=id; }
        }
    }else{
        if(want>=n){ s.resize(n); iota(s.begin(),s.end(),0); return s; }
        int last=-1;
        for(int i=0;i<want;i++){
            int id=(int)llround((long double)i*(n-1)/(long double)(want-1));
            id=max(0,min(n-1,id));
            if(id!=last){ s.push_back(id); last=id; }
        }
        if(s.front()!=0) s.insert(s.begin(),0);
        if(s.back()!=n-1) s.push_back(n-1);
    }
    s.erase(unique(s.begin(), s.end()), s.end());
    return s;
}

static bool tri_area_ok(const vector<Vec3>& V, int a,int b,int c){
    Vec3 cr=crossv(V[b]-V[a],V[c]-V[a]);
    return norm2v(cr) > bboxDiag*bboxDiag*bboxDiag*bboxDiag*1e-28;
}

static void remesh_grid_and_write(const GridInfo& g){
    double ratio;
    if(g.match>0.985 && g.faceErr<0.015) ratio = (N>=800000 ? 0.026 : 0.034);
    else ratio = 0.052;
    int target=(int)llround(N*ratio);
    target=max(9000, min(target, 70000));
    if(N<260000) target=max(target, 10500);

    double aspect=(double)g.W/(double)g.H;
    int wantH=(int)llround(sqrt(max(1.0,(double)target/aspect)));
    int wantW=(int)llround((double)target/max(1,wantH));
    wantH=max(g.wrapR?8:3, min(wantH,g.H));
    wantW=max(g.wrapC?8:3, min(wantW,g.W));
    // Rebalance after clamps.
    if((long long)wantH*wantW > target*13LL/10 && wantH>3 && wantW>3){
        double shrink=sqrt((double)target/((double)wantH*wantW));
        wantH=max(g.wrapR?8:3, min(g.H,(int)floor(wantH*shrink)));
        wantW=max(g.wrapC?8:3, min(g.W,(int)floor(wantW*shrink)));
    }

    vector<int> rows=make_samples_1d(g.H,wantH,g.wrapR);
    vector<int> cols=make_samples_1d(g.W,wantW,g.wrapC);
    int RH=(int)rows.size(), CW=(int)cols.size();
    vector<Vec3> V; V.reserve((size_t)RH*CW);
    for(int r: rows) for(int c: cols) V.push_back(P[r*g.W+c]);
    auto oid=[&](int r,int c)->int{
        if(g.wrapR){ r%=RH; if(r<0) r+=RH; }
        if(g.wrapC){ c%=CW; if(c<0) c+=CW; }
        return r*CW+c;
    };
    auto code_to_vid=[&](int rr,int cc,int code)->int{
        if(code==0) return oid(rr,cc);
        if(code==1) return oid(rr,cc+1);
        if(code==2) return oid(rr+1,cc);
        return oid(rr+1,cc+1);
    };
    vector<array<int,3>> T;
    int rCells=g.wrapR?RH:RH-1;
    int cCells=g.wrapC?CW:CW-1;
    T.reserve((size_t)max(0,rCells)*max(0,cCells)*2);
    for(int r=0;r<rCells;r++){
        for(int c=0;c<cCells;c++){
            for(int p=0;p<2;p++){
                int a=code_to_vid(r,c,g.patterns[p][0]);
                int b=code_to_vid(r,c,g.patterns[p][1]);
                int d=code_to_vid(r,c,g.patterns[p][2]);
                if(a==b||b==d||a==d) continue;
                if(!tri_area_ok(V,a,b,d)) continue;
                T.push_back({a,b,d});
            }
        }
    }
    write_mesh_vectors(V,T);
}

// ------------------------- strict ellipsoid/sphere remesh -------------------------

static bool detect_ellipsoid(){
    Vec3 rad=(bboxMax-bboxMin)*0.5;
    double rmin=min(rad.x,min(rad.y,rad.z));
    double rmax=max(rad.x,max(rad.y,rad.z));
    if(!(rmin>bboxDiag*1e-9) || rmax/rmin>6.0) return false;
    int step=max(1,N/220000);
    long double sum=0, sum2=0; double maxDev=0; int cnt=0;
    for(int i=0;i<N;i+=step){
        double x=(P[i].x-bboxCtr.x)/rad.x;
        double y=(P[i].y-bboxCtr.y)/rad.y;
        double z=(P[i].z-bboxCtr.z)/rad.z;
        double rr=sqrt(x*x+y*y+z*z);
        double dev=fabs(rr-1.0);
        sum += rr; sum2 += (long double)rr*rr; maxDev=max(maxDev,dev); cnt++;
    }
    if(cnt<1000) return false;
    double mean=(double)(sum/cnt);
    double var=max(0.0,(double)(sum2/cnt)-mean*mean);
    double sd=sqrt(var);
    // Very strict: only near-perfect UV spheres/ellipsoids should fire.
    return fabs(mean-1.0)<0.006 && sd<0.0045 && maxDev<0.024;
}

static void remesh_ellipsoid_and_write(){
    int target=(int)llround(N*(N>=800000?0.020:0.028));
    target=max(8000,min(target,52000));
    int lat=max(36,min(180,(int)llround(sqrt(target/2.0))));
    int lon=max(72,min(360,2*lat));
    Vec3 rad=(bboxMax-bboxMin)*0.5;
    vector<Vec3> V;
    V.reserve(2+(lat-1)*lon);
    V.push_back({bboxCtr.x,bboxCtr.y,bboxCtr.z+rad.z});
    const double PI=acos(-1.0);
    for(int i=1;i<lat;i++){
        double phi=PI*(double)i/(double)lat;
        double sp=sin(phi), cp=cos(phi);
        for(int j=0;j<lon;j++){
            double th=2.0*PI*(double)j/(double)lon;
            V.push_back({bboxCtr.x+rad.x*sp*cos(th), bboxCtr.y+rad.y*sp*sin(th), bboxCtr.z+rad.z*cp});
        }
    }
    int south=(int)V.size();
    V.push_back({bboxCtr.x,bboxCtr.y,bboxCtr.z-rad.z});
    auto ring=[&](int i,int j)->int{ // i in [1,lat-1]
        j%=lon; if(j<0) j+=lon;
        return 1+(i-1)*lon+j;
    };
    vector<array<int,3>> T;
    T.reserve((size_t)2*lon*lat);
    for(int j=0;j<lon;j++) T.push_back({0, ring(1,j), ring(1,j+1)});
    for(int i=1;i<lat-1;i++){
        for(int j=0;j<lon;j++){
            int a=ring(i,j), b=ring(i,j+1), c=ring(i+1,j), d=ring(i+1,j+1);
            T.push_back({a,c,b});
            T.push_back({b,c,d});
        }
    }
    for(int j=0;j<lon;j++) T.push_back({south, ring(lat-1,j+1), ring(lat-1,j)});
    write_mesh_vectors(V,T);
}

// ------------------------- generic safe bounded edge collapse fallback -------------------------

static vector<unsigned char> aliveV;
static vector<vector<int>> inc;
static vector<double> clusR;
static int activeV=0, activeF=0;
static vector<int> markSeen;
static int markTok=1;

static void cleanup_incident(int v){
    if(v<0 || v>=N) return;
    vector<int>& ids=inc[v];
    if(ids.size()<96) return;
    size_t good=0;
    for(int fid: ids) if(F[fid].active && face_contains(F[fid],v)) good++;
    if(good*3+32 >= ids.size()) return;
    vector<int> keep; keep.reserve(good+8);
    for(int fid: ids) if(F[fid].active && face_contains(F[fid],v)) keep.push_back(fid);
    ids.swap(keep);
}

static void init_collapse_state(){
    aliveV.assign(N,1);
    clusR.assign(N,0.0);
    vector<int> deg(N,0);
    activeF=0;
    for(int i=0;i<M;i++){
        Face& f=F[i];
        if(f.a<0||f.b<0||f.c<0||f.a>=N||f.b>=N||f.c>=N||f.a==f.b||f.a==f.c||f.b==f.c){ f.active=0; continue; }
        f.active=1; activeF++;
        deg[f.a]++; deg[f.b]++; deg[f.c]++;
    }
    inc.assign(N,{});
    for(int i=0;i<N;i++) inc[i].reserve(deg[i]);
    for(int i=0;i<M;i++) if(F[i].active){ inc[F[i].a].push_back(i); inc[F[i].b].push_back(i); inc[F[i].c].push_back(i); }
    activeV=N;
    markSeen.assign(N,0);
}

static int third_vertex(const Face& f,int a,int b){
    if(f.a!=a && f.a!=b) return f.a;
    if(f.b!=a && f.b!=b) return f.b;
    if(f.c!=a && f.c!=b) return f.c;
    return -1;
}

static bool collect_edge_faces(int a,int b,int ef[2],int opp[2]){
    if(a==b || !aliveV[a] || !aliveV[b]) return false;
    cleanup_incident(a); cleanup_incident(b);
    const vector<int>& A = inc[a].size() < inc[b].size() ? inc[a] : inc[b];
    int cnt=0;
    for(int fid: A){
        if(!F[fid].active) continue;
        const Face& f=F[fid];
        if(!face_contains(f,a) || !face_contains(f,b)) continue;
        if(cnt>=2) return false;
        int o=third_vertex(f,a,b);
        if(o<0) return false;
        ef[cnt]=fid; opp[cnt]=o; cnt++;
    }
    return cnt==2 && opp[0]!=opp[1];
}

static bool link_ok(int a,int b,const int opp[2]){
    if(markTok > 2000000000){ fill(markSeen.begin(),markSeen.end(),0); markTok=1; }
    int tokA=markTok++, tokC=markTok++;
    for(int fid: inc[a]){
        if(!F[fid].active) continue;
        const Face& f=F[fid]; if(!face_contains(f,a)) continue;
        int vs[3]={f.a,f.b,f.c};
        for(int x: vs) if(x!=a && x!=b) markSeen[x]=tokA;
    }
    int common=0, got0=0, got1=0;
    for(int fid: inc[b]){
        if(!F[fid].active) continue;
        const Face& f=F[fid]; if(!face_contains(f,b)) continue;
        int vs[3]={f.a,f.b,f.c};
        for(int x: vs){
            if(x==a||x==b) continue;
            if(markSeen[x]==tokA){
                markSeen[x]=tokC; common++;
                if(x==opp[0]) got0=1;
                if(x==opp[1]) got1=1;
                if(common>2) return false;
            }
        }
    }
    return common==2 && got0 && got1;
}

static bool existing_duplicate_after_replace(int /*keep*/,int rem,int fidSkip0,int fidSkip1,int changedFid,int a,int b,int c){
    array<int,3> s=sorted_tri(a,b,c);
    int probe=a;
    if(inc[b].size()<inc[probe].size()) probe=b;
    if(inc[c].size()<inc[probe].size()) probe=c;
    for(int fid: inc[probe]){
        if(fid==fidSkip0 || fid==fidSkip1 || fid==changedFid) continue;
        if(!F[fid].active) continue;
        const Face& f=F[fid];
        if(face_contains(f,rem)) continue; // will be rewritten in this operation; checked separately.
        if(sorted_tri(f.a,f.b,f.c)==s) return true;
    }
    return false;
}

struct CollapsePlan { bool ok=false; int keep=-1, rem=-1; double newR=0, cost=1e100; };

static CollapsePlan try_directed_collapse(int keep,int rem,const int ef[2],double errLimit,double cosMin){
    CollapsePlan plan; plan.keep=keep; plan.rem=rem;
    double d=sqrt(dist2v(P[keep],P[rem]));
    double nr=max(clusR[keep], clusR[rem]+d);
    if(nr>errLimit) return plan;
    vector<array<int,3>> changed;
    changed.reserve(inc[rem].size());
    const double tinyArea = bboxDiag*bboxDiag*bboxDiag*bboxDiag*1e-30;
    for(int fid: inc[rem]){
        if(!F[fid].active) continue;
        const Face& f=F[fid];
        if(!face_contains(f,rem)) continue;
        if(fid==ef[0] || fid==ef[1]) continue;
        if(face_contains(f,keep)) return plan;
        int a=f.a,b=f.b,c=f.c;
        if(a==rem) a=keep; else if(b==rem) b=keep; else if(c==rem) c=keep;
        if(a==b||a==c||b==c) return plan;
        Vec3 oldN=tri_normal_idx(f.a,f.b,f.c);
        Vec3 newN=tri_normal_idx(a,b,c);
        double lo=norm2v(oldN), ln=norm2v(newN);
        if(lo<=tinyArea || ln<=tinyArea) return plan;
        double co=dotv(oldN,newN)/sqrt(lo*ln);
        if(co < cosMin) return plan;
        if(ln < lo*1e-5) return plan;
        array<int,3> s=sorted_tri(a,b,c);
        for(const auto& q: changed) if(q==s) return plan;
        if(existing_duplicate_after_replace(keep,rem,ef[0],ef[1],fid,a,b,c)) return plan;
        changed.push_back(s);
    }
    plan.ok=true;
    plan.newR=nr;
    plan.cost=nr + 0.001*d + 0.000002*(double)(inc[keep].size()+inc[rem].size());
    return plan;
}

static void apply_collapse_plan(const CollapsePlan& p,const int ef[2]){
    int keep=p.keep, rem=p.rem;
    for(int i=0;i<2;i++) if(ef[i]>=0 && F[ef[i]].active){ F[ef[i]].active=0; activeF--; }
    vector<int> touched;
    touched.reserve(inc[rem].size()+8);
    for(int fid: inc[rem]){
        if(!F[fid].active) continue;
        Face& f=F[fid];
        if(!face_contains(f,rem)) continue;
        if(f.a==rem) f.a=keep;
        if(f.b==rem) f.b=keep;
        if(f.c==rem) f.c=keep;
        if(f.a==f.b || f.a==f.c || f.b==f.c){ f.active=0; activeF--; continue; }
        inc[keep].push_back(fid);
        touched.push_back(fid);
    }
    aliveV[rem]=0; activeV--; clusR[keep]=p.newR;
    inc[rem].clear();
    cleanup_incident(keep);
}

struct EdgeCand { float d2; uint64_t key; };
static bool by_d2(const EdgeCand&a,const EdgeCand&b){ return a.d2<b.d2; }

static double estimate_mean_edge_len(){
    vector<uint64_t> keys; keys.reserve((size_t)activeF*3);
    for(int i=0;i<M;i++) if(F[i].active){
        const Face& f=F[i];
        keys.push_back(edge_key(f.a,f.b)); keys.push_back(edge_key(f.b,f.c)); keys.push_back(edge_key(f.c,f.a));
    }
    sort(keys.begin(),keys.end()); keys.erase(unique(keys.begin(),keys.end()),keys.end());
    long double sum=0; long long cnt=0;
    int step=max(1,(int)(keys.size()/300000));
    for(size_t i=0;i<keys.size();i+=step){
        int a=key_a(keys[i]), b=key_b(keys[i]);
        if(a>=0&&b>=0&&a<N&&b<N) { sum += sqrt(dist2v(P[a],P[b])); cnt++; }
    }
    if(cnt==0) return bboxDiag / sqrt((double)max(1,N));
    return max(1e-18,(double)(sum/cnt));
}

static void collect_candidates_under(double lenLimit, vector<EdgeCand>& cand){
    vector<uint64_t> keys; keys.reserve((size_t)activeF*3);
    for(int i=0;i<M;i++) if(F[i].active){
        const Face& f=F[i];
        if(!aliveV[f.a]||!aliveV[f.b]||!aliveV[f.c]) continue;
        keys.push_back(edge_key(f.a,f.b)); keys.push_back(edge_key(f.b,f.c)); keys.push_back(edge_key(f.c,f.a));
    }
    sort(keys.begin(),keys.end());
    double lim2=lenLimit*lenLimit;
    cand.clear(); cand.reserve(keys.size()/4+16);
    uint64_t last=~0ull;
    for(uint64_t k: keys){
        if(k==last) continue;
        last=k;
        int a=key_a(k), b=key_b(k);
        if(a<0||b<0||a>=N||b>=N||!aliveV[a]||!aliveV[b]) continue;
        double d2=dist2v(P[a],P[b]);
        if(d2<=lim2) cand.push_back({(float)d2,k});
    }
    vector<uint64_t>().swap(keys);
    sort(cand.begin(),cand.end(),by_d2);
}

static void repair_duplicates_and_edges(){
    struct FT { array<int,3> s; int fid; bool operator<(const FT&o)const{return s<o.s;} };
    vector<FT> ts; ts.reserve(activeF);
    for(int i=0;i<M;i++) if(F[i].active){
        Face& f=F[i];
        if(f.a==f.b||f.a==f.c||f.b==f.c){ f.active=0; activeF--; continue; }
        ts.push_back({sorted_tri(f.a,f.b,f.c),i});
    }
    sort(ts.begin(),ts.end());
    for(size_t i=1;i<ts.size();i++){
        if(ts[i].s==ts[i-1].s && F[ts[i].fid].active){ F[ts[i].fid].active=0; activeF--; }
    }
    struct EF { uint64_t e; int fid; bool operator<(const EF&o)const{ if(e!=o.e) return e<o.e; return fid<o.fid; } };
    vector<EF> es; es.reserve((size_t)activeF*3);
    for(int i=0;i<M;i++) if(F[i].active){
        const Face& f=F[i];
        es.push_back({edge_key(f.a,f.b),i}); es.push_back({edge_key(f.b,f.c),i}); es.push_back({edge_key(f.c,f.a),i});
    }
    sort(es.begin(),es.end());
    for(size_t i=0;i<es.size();){
        size_t j=i+1; while(j<es.size() && es[j].e==es[i].e) j++;
        if(j-i>2){
            // Keep the first two active incident faces; drop extras to avoid non-manifold edges.
            int kept=0;
            for(size_t k=i;k<j;k++) if(F[es[k].fid].active){
                if(kept<2) kept++; else { F[es[k].fid].active=0; activeF--; }
            }
        }
        i=j;
    }
}

static void write_active_collapse_mesh(){
    repair_duplicates_and_edges();
    vector<int> remap(N,-1);
    int nv=0;
    for(int i=0;i<N;i++) if(aliveV[i]) remap[i]=nv++;
    vector<Vec3> V; V.reserve(nv);
    for(int i=0;i<N;i++) if(aliveV[i]) V.push_back(P[i]);
    vector<array<int,3>> T; T.reserve(activeF);
    for(int i=0;i<M;i++) if(F[i].active){
        int a=F[i].a,b=F[i].b,c=F[i].c;
        if(a==b||a==c||b==c) continue;
        if(a<0||b<0||c<0||a>=N||b>=N||c>=N) continue;
        if(remap[a]<0||remap[b]<0||remap[c]<0) continue;
        T.push_back({remap[a],remap[b],remap[c]});
    }
    write_mesh_vectors(V,T);
}

static void safe_collapse_fallback_and_write(){
    init_collapse_state();
    int targetV = min((int)llround(N*0.24), 260000);
    if(N>=650000) targetV=min(targetV,220000);
    if(N>=1200000) targetV=min(targetV,180000);
    targetV=max(8000,targetV);

    double meanEdge=estimate_mean_edge_len();
    double errRel = (N>=800000 ? 0.0085 : 0.0065);
    double errLimit = bboxDiag * errRel;
    double mults[] = {1.55, 2.20, 3.10, 4.50, 6.50, 9.00, 12.50};
    vector<EdgeCand> cand;

    for(int pass=0; pass<7 && activeV>targetV && time_ok(16.9); ++pass){
        double lenLimit=min(errLimit, meanEdge*mults[pass]);
        collect_candidates_under(lenLimit,cand);
        double cosMin = (pass<=2 ? 0.965 : (pass<=4 ? 0.940 : 0.905));
        long long tries=0, done=0;
        for(const EdgeCand& ec: cand){
            if(activeV<=targetV) break;
            if((++tries & 8191)==0 && !time_ok(17.35)) break;
            int a=key_a(ec.key), b=key_b(ec.key);
            if(a<0||b<0||a>=N||b>=N||!aliveV[a]||!aliveV[b]) continue;
            int ef[2]={-1,-1}, opp[2]={-1,-1};
            if(!collect_edge_faces(a,b,ef,opp)) continue;
            if(!link_ok(a,b,opp)) continue;
            CollapsePlan ab=try_directed_collapse(a,b,ef,errLimit,cosMin);
            CollapsePlan ba=try_directed_collapse(b,a,ef,errLimit,cosMin);
            CollapsePlan p;
            if(ab.ok && (!ba.ok || ab.cost<=ba.cost)) p=ab; else if(ba.ok) p=ba; else continue;
            apply_collapse_plan(p,ef);
            done++;
        }
        cand.clear();
        if(done==0 && pass>=3) break;
    }
    write_active_collapse_mesh();
}

int main(){
    startTime=chrono::steady_clock::now();
    load_mesh();

    // Do not perturb small/mid cases. The known high-N bottleneck starts at 200k.
    if(N < HIGH_N_TRIGGER){
        write_identity();
        return 0;
    }

    GridInfo grid=detect_index_grid();
    if(grid.ok && grid.match>=GRID_STRICT_MATCH && time_ok(14.0)){
        remesh_grid_and_write(grid);
        return 0;
    }

    if(time_ok(10.5) && detect_ellipsoid()){
        remesh_ellipsoid_and_write();
        return 0;
    }

    // Generic safe fallback: no identity at high N, but all collapses are local-manifold checked.
    safe_collapse_fallback_and_write();
    return 0;
}
