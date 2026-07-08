#include <bits/stdc++.h>
using namespace std;

// IMC Challenge 2026 - simplifygeometry
// C++17 single-file heuristic: prefix-tolerant parser, topology-safe QEM fallback,
// row-major periodic grid remesher, and strict snapped ellipsoid remesher.

struct Vec3{
    double x=0,y=0,z=0;
};
static inline Vec3 operator+(const Vec3&a,const Vec3&b){return {a.x+b.x,a.y+b.y,a.z+b.z};}
static inline Vec3 operator-(const Vec3&a,const Vec3&b){return {a.x-b.x,a.y-b.y,a.z-b.z};}
static inline Vec3 operator*(const Vec3&a,double s){return {a.x*s,a.y*s,a.z*s};}
static inline Vec3 operator/(const Vec3&a,double s){return {a.x/s,a.y/s,a.z/s};}
static inline double dot3(const Vec3&a,const Vec3&b){return a.x*b.x+a.y*b.y+a.z*b.z;}
static inline Vec3 cross3(const Vec3&a,const Vec3&b){return {a.y*b.z-a.z*b.y,a.z*b.x-a.x*b.z,a.x*b.y-a.y*b.x};}
static inline double norm2(const Vec3&a){return dot3(a,a);} 
static inline double norm3(const Vec3&a){return sqrt(max(0.0,norm2(a)));}
static inline double dist2(const Vec3&a,const Vec3&b){return norm2(a-b);} 
static inline bool finite3(const Vec3&p){return isfinite(p.x)&&isfinite(p.y)&&isfinite(p.z);} 

struct Face{int v[3]; unsigned char active=1;};

struct FastScanner{
    vector<char> buf; char* p=nullptr;
    FastScanner(){
        buf.reserve(1<<26);
        char tmp[1<<16]; size_t n;
        while((n=fread(tmp,1,sizeof(tmp),stdin))>0) buf.insert(buf.end(),tmp,tmp+n);
        buf.push_back('\0'); p=buf.data();
    }
    inline void skip(){ while(*p==' '||*p=='\n'||*p=='\r'||*p=='\t') ++p; }
    string next_token(){ skip(); string s; while(*p && *p!=' '&&*p!='\n'&&*p!='\r'&&*p!='\t') s.push_back(*p++); return s; }
    long long next_ll(){ string s=next_token(); return strtoll(s.c_str(),nullptr,10); }
    double next_double(){ string s=next_token(); return strtod(s.c_str(),nullptr); }
};
static inline bool is_prefix_token(const string&s){ return !s.empty() && isalpha((unsigned char)s[0]); }
static inline int parse_index_token(const string&s){
    const char* c=s.c_str(); char* e=nullptr; long v=strtol(c,&e,10); return (int)v;
}

struct Quadric{
    double q[10];
    Quadric(){ memset(q,0,sizeof(q)); }
    inline void add_plane(double a,double b,double c,double d,double w=1.0){
        q[0]+=w*a*a; q[1]+=w*a*b; q[2]+=w*a*c; q[3]+=w*a*d; q[4]+=w*b*b;
        q[5]+=w*b*c; q[6]+=w*b*d; q[7]+=w*c*c; q[8]+=w*c*d; q[9]+=w*d*d;
    }
    inline void add(const Quadric&o){ for(int i=0;i<10;i++) q[i]+=o.q[i]; }
    inline double eval(const Vec3&p)const{
        double x=p.x,y=p.y,z=p.z;
        return q[0]*x*x+2*q[1]*x*y+2*q[2]*x*z+2*q[3]*x+q[4]*y*y+2*q[5]*y*z+2*q[6]*y+q[7]*z*z+2*q[8]*z+q[9];
    }
};
struct Node{ double cost; int a,b,va,vb; bool operator<(const Node&o)const{return cost>o.cost;} };
struct EdgeRef{ uint64_t key; int f; bool operator<(const EdgeRef&o)const{return key<o.key || (key==o.key&&f<o.f);} };

static int N=0,M=0,activeV=0,activeF=0;
static vector<Vec3> Orig,P;
static vector<Face> F;
static vector<vector<int>> Inc;
static vector<Quadric> Q;
static vector<unsigned char> Alive,Locked;
static vector<int> Head,Tail,NextMem,ClusterSize,Version,MarkV,MarkF;
static vector<double> Radius;
static int tokenV=3, tokenF=5;
static double diagLen=1.0, epsLen=1.0, eps2=1.0, areaEps2=1e-30;
static Vec3 bboxMin,bboxMax;
static bool outPrefixed=false;
static chrono::steady_clock::time_point T0;
static priority_queue<Node> PQ;

static inline double elapsed(){ return chrono::duration<double>(chrono::steady_clock::now()-T0).count(); }
static inline bool hasv(const Face&f,int v){return f.v[0]==v||f.v[1]==v||f.v[2]==v;}
static inline uint64_t ekey(int a,int b){ if(a>b) swap(a,b); return (uint64_t)(uint32_t)a<<32 | (uint32_t)b; }
static inline int ka(uint64_t k){return (int)(k>>32);} static inline int kb(uint64_t k){return (int)(k&0xffffffffu);} 
static inline array<int,3> tri_sorted(int a,int b,int c){ array<int,3> t{a,b,c}; sort(t.begin(),t.end()); return t; }
static inline unsigned long long pack_tri(array<int,3> t){ return ((unsigned long long)(unsigned)t[0]<<42)^((unsigned long long)(unsigned)t[1]<<21)^(unsigned)t[2]; }

static bool solve3(const Quadric&q,Vec3&out){
    double a00=q.q[0],a01=q.q[1],a02=q.q[2],a10=q.q[1],a11=q.q[4],a12=q.q[5],a20=q.q[2],a21=q.q[5],a22=q.q[7];
    double b0=-q.q[3],b1=-q.q[6],b2=-q.q[8];
    double det=a00*(a11*a22-a12*a21)-a01*(a10*a22-a12*a20)+a02*(a10*a21-a11*a20);
    double sc=max(1.0,fabs(a00)+fabs(a11)+fabs(a22)+fabs(a01)+fabs(a02)+fabs(a12));
    if(fabs(det)<1e-18*sc*sc*sc) return false;
    double dx=b0*(a11*a22-a12*a21)-a01*(b1*a22-a12*b2)+a02*(b1*a21-a11*b2);
    double dy=a00*(b1*a22-a12*b2)-b0*(a10*a22-a12*a20)+a02*(a10*b2-b1*a20);
    double dz=a00*(a11*b2-b1*a21)-a01*(a10*b2-b1*a20)+b0*(a10*a21-a11*a20);
    out={dx/det,dy/det,dz/det}; return finite3(out);
}

static void read_input(){
    FastScanner in;
    string sN=in.next_token(); if(sN.empty()) return;
    N=(int)strtol(sN.c_str(),nullptr,10);
    M=(int)in.next_ll();
    Orig.resize(N); P.resize(N);
    bboxMin={1e100,1e100,1e100}; bboxMax={-1e100,-1e100,-1e100};
    for(int i=0;i<N;i++){
        string t=in.next_token(); double x,y,z;
        if(is_prefix_token(t)){ outPrefixed=true; x=in.next_double(); y=in.next_double(); z=in.next_double(); }
        else{ x=strtod(t.c_str(),nullptr); y=in.next_double(); z=in.next_double(); }
        Orig[i]=P[i]={x,y,z};
        bboxMin.x=min(bboxMin.x,x); bboxMin.y=min(bboxMin.y,y); bboxMin.z=min(bboxMin.z,z);
        bboxMax.x=max(bboxMax.x,x); bboxMax.y=max(bboxMax.y,y); bboxMax.z=max(bboxMax.z,z);
    }
    F.assign(M,Face()); vector<int> deg(max(0,N),0);
    for(int i=0;i<M;i++){
        string t=in.next_token(); int a,b,c;
        if(is_prefix_token(t)){ outPrefixed=true; a=parse_index_token(in.next_token()); b=parse_index_token(in.next_token()); c=parse_index_token(in.next_token()); }
        else{ a=parse_index_token(t); b=parse_index_token(in.next_token()); c=parse_index_token(in.next_token()); }
        --a;--b;--c;
        F[i].v[0]=a; F[i].v[1]=b; F[i].v[2]=c; F[i].active=1;
        if(0<=a&&a<N)deg[a]++; if(0<=b&&b<N)deg[b]++; if(0<=c&&c<N)deg[c]++;
    }
    diagLen=max(1e-30,norm3(bboxMax-bboxMin));
    epsLen=diagLen*0.04975; eps2=epsLen*epsLen; areaEps2=max(1e-30,diagLen*diagLen*1e-28);
    Inc.assign(N,{}); for(int i=0;i<N;i++) Inc[i].reserve(deg[i]+4);
    for(int i=0;i<M;i++) for(int k=0;k<3;k++) if(0<=F[i].v[k]&&F[i].v[k]<N) Inc[F[i].v[k]].push_back(i);
    Q.assign(N,Quadric()); Alive.assign(N,1); Locked.assign(N,0); Radius.assign(N,0.0);
    Head.resize(N); Tail.resize(N); NextMem.assign(N,-1); ClusterSize.assign(N,1); Version.assign(N,0); MarkV.assign(N,0); MarkF.assign(max(1,M),0);
    for(int i=0;i<N;i++) Head[i]=Tail[i]=i;
    activeV=N; activeF=M;
    if(N>0){
        int ix0=0,ix1=0,iy0=0,iy1=0,iz0=0,iz1=0;
        for(int i=0;i<N;i++){
            if(Orig[i].x<Orig[ix0].x)ix0=i; if(Orig[i].x>Orig[ix1].x)ix1=i;
            if(Orig[i].y<Orig[iy0].y)iy0=i; if(Orig[i].y>Orig[iy1].y)iy1=i;
            if(Orig[i].z<Orig[iz0].z)iz0=i; if(Orig[i].z>Orig[iz1].z)iz1=i;
        }
        Locked[ix0]=Locked[ix1]=Locked[iy0]=Locked[iy1]=Locked[iz0]=Locked[iz1]=1;
    }
}

static bool validate_mesh_arrays(const vector<Vec3>&X,const vector<Face>&A,bool requireClosed=true){
    if(X.empty()||A.empty()||X.size()>(size_t)max(1,N)) return false;
    unordered_map<uint64_t,int> ec; ec.reserve(A.size()*4+16);
    double ae=max(1e-30,diagLen*diagLen*1e-28);
    for(const Face&f:A){
        int a=f.v[0],b=f.v[1],c=f.v[2];
        if(a<0||b<0||c<0||a>=(int)X.size()||b>=(int)X.size()||c>=(int)X.size()) return false;
        if(a==b||a==c||b==c) return false;
        if(norm2(cross3(X[b]-X[a],X[c]-X[a]))<=ae) return false;
        ec[ekey(a,b)]++; ec[ekey(b,c)]++; ec[ekey(c,a)]++;
    }
    if(requireClosed){ for(auto &kv:ec) if(kv.second!=2) return false; }
    return true;
}

static void print_mesh(const vector<Vec3>&X,const vector<Face>&A){
    static char obuf[1<<20]; setvbuf(stdout,obuf,_IOFBF,sizeof(obuf));
    printf("%d %d\n",(int)X.size(),(int)A.size());
    if(outPrefixed){
        for(const auto&p:X) printf("v %.15g %.15g %.15g\n",p.x,p.y,p.z);
        for(const Face&f:A) printf("f %d %d %d\n",f.v[0]+1,f.v[1]+1,f.v[2]+1);
    }else{
        for(const auto&p:X) printf("%.15g %.15g %.15g\n",p.x,p.y,p.z);
        for(const Face&f:A) printf("%d %d %d\n",f.v[0]+1,f.v[1]+1,f.v[2]+1);
    }
}
static void print_original(){
    vector<Face>A=F; for(auto&f:A) f.active=1; print_mesh(Orig,A);
}

struct KDTree{
    const vector<Vec3>* pts=nullptr; vector<int> id;
    void init(const vector<Vec3>&p){ pts=&p; id.resize(p.size()); iota(id.begin(),id.end(),0); build(0,(int)id.size(),0); }
    static double coord(const Vec3&p,int ax){ return ax==0?p.x:(ax==1?p.y:p.z); }
    void build(int l,int r,int d){ if(r-l<=1) return; int m=(l+r)/2, ax=d%3; nth_element(id.begin()+l,id.begin()+m,id.begin()+r,[&](int a,int b){return coord((*pts)[a],ax)<coord((*pts)[b],ax);}); build(l,m,d+1); build(m+1,r,d+1); }
    void nearest_rec(const Vec3&q,int l,int r,int d,int&best,double&bd2)const{
        if(l>=r) return; int m=(l+r)/2, ax=d%3; int pi=id[m]; double dd=dist2(q,(*pts)[pi]); if(dd<bd2){bd2=dd; best=pi;}
        double df=coord(q,ax)-coord((*pts)[pi],ax); int l1=l,r1=m,l2=m+1,r2=r; if(df>0){swap(l1,l2); swap(r1,r2);} nearest_rec(q,l1,r1,d+1,best,bd2); if(df*df<bd2) nearest_rec(q,l2,r2,d+1,best,bd2);
    }
    int nearest(const Vec3&q,double*bd=nullptr)const{ int b=-1; double d=1e300; nearest_rec(q,0,(int)id.size(),0,b,d); if(bd)*bd=d; return b; }
};

static bool coverage_points(const vector<Vec3>&S){
    if(S.empty()) return false;
    KDTree kd; kd.init(S); double lim2=eps2*0.999;
    for(int i=0;i<N;i++){ double b; kd.nearest(Orig[i],&b); if(b>lim2) return false; if((i&16383)==0 && elapsed()>7.0) return false; }
    return true;
}

static void orient_face_to_center(vector<Vec3>&X,Face&f,const Vec3&center){
    Vec3 n=cross3(X[f.v[1]]-X[f.v[0]],X[f.v[2]]-X[f.v[0]]);
    Vec3 c=(X[f.v[0]]+X[f.v[1]]+X[f.v[2]])/3.0;
    if(dot3(n,c-center)<0) swap(f.v[1],f.v[2]);
}

static bool try_snapped_ellipsoid(){
    if(N<700 || elapsed()>1.5) return false;
    Vec3 cen=(bboxMin+bboxMax)/2.0; Vec3 rad=(bboxMax-bboxMin)/2.0;
    if(rad.x<diagLen*1e-8||rad.y<diagLen*1e-8||rad.z<diagLen*1e-8) return false;
    double maxAbs=0, avgAbs=0; int sample=0;
    int step=max(1,N/60000);
    for(int i=0;i<N;i+=step){
        double sx=(Orig[i].x-cen.x)/rad.x, sy=(Orig[i].y-cen.y)/rad.y, sz=(Orig[i].z-cen.z)/rad.z;
        double r=sqrt(sx*sx+sy*sy+sz*sz); double dev=fabs(r-1.0);
        maxAbs=max(maxAbs,dev); avgAbs+=dev; sample++;
    }
    avgAbs/=max(1,sample);
    // Strict: only true spheres/ellipsoids should enter this branch.
    if(maxAbs>0.040 || avgAbs>0.010) return false;
    KDTree origKD; origKD.init(Orig);
    double rmean=(rad.x+rad.y+rad.z)/3.0;
    double baseStep=max(epsLen*1.35, diagLen*0.012);
    int baseLat=max(8,(int)ceil(M_PI*rmean/baseStep));
    vector<pair<int,int>> tries;
    for(double mul: {0.85,1.0,1.18,1.38,1.62}){
        int lat=max(6,(int)ceil(baseLat*mul)); int lon=max(12,2*lat);
        if(2+(lat-1)*lon >= N) continue;
        tries.push_back({lat,lon});
    }
    for(auto [lat,lon]: tries){
        if(elapsed()>6.5) break;
        vector<int> oldId; oldId.reserve(2+(lat-1)*lon);
        auto makep=[&](double theta,double phi){
            double st=sin(theta), ct=cos(theta), cp=cos(phi), sp=sin(phi);
            return Vec3{cen.x+rad.x*st*cp, cen.y+rad.y*st*sp, cen.z+rad.z*ct};
        };
        double bd; int north=origKD.nearest(Vec3{cen.x,cen.y,cen.z+rad.z},&bd); oldId.push_back(north);
        int south=origKD.nearest(Vec3{cen.x,cen.y,cen.z-rad.z},&bd); if(south==north) continue; oldId.push_back(south);
        bool dup=false; unordered_set<int> used; used.reserve((size_t)(lat-1)*lon*2); used.insert(north); used.insert(south);
        for(int r=1;r<=lat-1 && !dup;r++){
            double theta=M_PI*(double)r/(double)lat;
            for(int j=0;j<lon;j++){
                double phi=2.0*M_PI*(double)j/(double)lon; int id=origKD.nearest(makep(theta,phi),&bd);
                if(!used.insert(id).second){ dup=true; break; }
                oldId.push_back(id);
            }
        }
        if(dup) continue;
        vector<Vec3>X; X.reserve(oldId.size()); for(int id:oldId) X.push_back(Orig[id]);
        if(!coverage_points(X)) continue;
        auto rid=[&](int r,int j){ j=(j%lon+lon)%lon; return 2+(r-1)*lon+j; };
        vector<Face>A; A.reserve((size_t)2*lon*(lat-1));
        auto add=[&](int a,int b,int c){ Face f{{a,b,c},1}; orient_face_to_center(X,f,cen); A.push_back(f); };
        for(int j=0;j<lon;j++) add(0,rid(1,j),rid(1,j+1));
        for(int r=1;r<lat-1;r++) for(int j=0;j<lon;j++){ add(rid(r,j),rid(r+1,j),rid(r,j+1)); add(rid(r,j+1),rid(r+1,j),rid(r+1,j+1)); }
        for(int j=0;j<lon;j++) add(1,rid(lat-1,j+1),rid(lat-1,j));
        if(validate_mesh_arrays(X,A,true)){ print_mesh(X,A); return true; }
    }
    return false;
}

static bool face_set_eq(const Face&f, array<int,3> t){
    array<int,3> a{f.v[0],f.v[1],f.v[2]}; sort(a.begin(),a.end()); sort(t.begin(),t.end()); return a==t;
}
static vector<int> grid_candidates(){
    vector<int> c; if(N<=0) return c;
    int sq=(int)sqrt((double)N);
    for(int d=3; d<=sq+3; ++d) if(d>0 && N%d==0){ c.push_back(d); if(d!=N/d) c.push_back(N/d); }
    sort(c.begin(),c.end(),[&](int a,int b){ return abs(a-sq)<abs(b-sq); });
    c.erase(unique(c.begin(),c.end()),c.end()); return c;
}
static bool good_grid(int V){
    if(V<3||N%V) return false; int U=N/V; if(U<3) return false; if(M<2*N-8||M>2*N+8) return false;
    long long tot=0,ok=0;
    int limCells=min(N,200000);
    int stride=max(1,N/limCells);
    for(int cell=0; cell<N; cell+=stride){
        int r=cell/V, j=cell%V; int fid=2*cell; if(fid+1>=M) break;
        int a=r*V+j, b=((r+1)%U)*V+j, c=((r+1)%U)*V+(j+1)%V, d=r*V+(j+1)%V;
        bool one = face_set_eq(F[fid],{a,b,d}) && face_set_eq(F[fid+1],{b,c,d});
        bool two = face_set_eq(F[fid],{a,b,c}) && face_set_eq(F[fid+1],{a,c,d});
        bool three = face_set_eq(F[fid],{a,c,d}) && face_set_eq(F[fid+1],{a,b,c});
        if(one||two||three) ok++; tot++;
    }
    return tot>100 && ok*1000>=tot*995;
}
static vector<int> make_sel(int full,int want,const vector<int>&must){
    vector<int>s; want=max(4,min(full,want)); s.reserve(want+must.size()+4);
    for(int i=0;i<want;i++) s.push_back((int)((long long)i*full/want));
    for(int x:must) if(0<=x&&x<full) s.push_back(x);
    sort(s.begin(),s.end()); s.erase(unique(s.begin(),s.end()),s.end()); return s;
}
static int cyclic_dist(int a,int b,int n){ int d=abs(a-b); return min(d,n-d); }
static bool coverage_grid(int U,int V,const vector<int>&rs,const vector<int>&cs){
    if(rs.empty()||cs.empty()) return false; double lim2=eps2*0.999;
    for(int i=0;i<U;i++){
        int r0=rs[0],r1=rs[0],bd=INT_MAX;
        for(int r:rs){ int d=cyclic_dist(i,r,U); if(d<bd){bd=d;r0=r;} }
        // second candidate is enough for smooth grids; include neighbor in selected list if any.
        r1=r0;
        for(int j=0;j<V;j++){
            int c0=cs[0],cd=INT_MAX; for(int c:cs){ int d=cyclic_dist(j,c,V); if(d<cd){cd=d;c0=c;} }
            int id=i*V+j; double best=dist2(Orig[id],Orig[r0*V+c0]);
            // Try adjacent selected rows/cols around nearest for robustness.
            for(int rr:rs) if(cyclic_dist(i,rr,U)<=bd+1) for(int cc:cs) if(cyclic_dist(j,cc,V)<=cd+1) best=min(best,dist2(Orig[id],Orig[rr*V+cc]));
            if(best>lim2) return false;
        }
        if((i&255)==0 && elapsed()>6.0) return false;
    }
    return true;
}
static void add_oriented(vector<Vec3>&X,vector<Face>&A,int a,int b,int c,const Vec3&ref){ Face f{{a,b,c},1}; Vec3 n=cross3(X[b]-X[a],X[c]-X[a]); if(dot3(n,ref)<0) swap(f.v[1],f.v[2]); A.push_back(f); }
static bool try_regular_grid(){
    if(N<1500 || elapsed()>2.5) return false;
    int V=0; for(int x:grid_candidates()) if(good_grid(x)){ V=x; break; }
    if(!V) return false; int U=N/V;
    vector<int> mustR,mustC;
    for(int i=0;i<N;i++){
        if(fabs(Orig[i].x-bboxMin.x)<1e-12||fabs(Orig[i].x-bboxMax.x)<1e-12||fabs(Orig[i].y-bboxMin.y)<1e-12||fabs(Orig[i].y-bboxMax.y)<1e-12||fabs(Orig[i].z-bboxMin.z)<1e-12||fabs(Orig[i].z-bboxMax.z)<1e-12){ mustR.push_back(i/V); mustC.push_back(i%V); }
    }
    double ratios[]={0.075,0.090,0.110,0.135,0.165,0.210,0.260};
    for(double rr:ratios){
        if(elapsed()>7.0) break;
        double s=sqrt(rr); int ru=max(6,min(U,(int)ceil(U*s))), cv=max(6,min(V,(int)ceil(V*s)));
        vector<int>rs=make_sel(U,ru,mustR), cs=make_sel(V,cv,mustC);
        if((long long)rs.size()*cs.size()>=N) continue;
        if(!coverage_grid(U,V,rs,cs)) continue;
        vector<Vec3>X; vector<Face>A; X.reserve(rs.size()*cs.size());
        for(int r:rs) for(int c:cs) X.push_back(Orig[r*V+c]);
        auto id=[&](int i,int j){ return ((i+(int)rs.size())%(int)rs.size())*(int)cs.size()+((j+(int)cs.size())%(int)cs.size()); };
        A.reserve(2*X.size());
        for(int i=0;i<(int)rs.size();i++) for(int j=0;j<(int)cs.size();j++){
            int oi=rs[i], oj=cs[j], fid=2*(oi*V+oj); Vec3 ref{0,0,0};
            if(fid<M){ Face fo=F[fid]; ref=ref+cross3(Orig[fo.v[1]]-Orig[fo.v[0]],Orig[fo.v[2]]-Orig[fo.v[0]]); }
            if(fid+1<M){ Face fo=F[fid+1]; ref=ref+cross3(Orig[fo.v[1]]-Orig[fo.v[0]],Orig[fo.v[2]]-Orig[fo.v[0]]); }
            int a=id(i,j),b=id(i+1,j),c=id(i+1,j+1),d=id(i,j+1);
            add_oriented(X,A,a,b,d,ref); add_oriented(X,A,b,c,d,ref);
        }
        if(validate_mesh_arrays(X,A,true)){ print_mesh(X,A); return true; }
    }
    return false;
}

static bool collect_edge(int a,int b,int ef[2],int op[2]){
    int cnt=0;
    const vector<int>&L = Inc[a].size()<Inc[b].size()?Inc[a]:Inc[b];
    for(int fid:L){ if(fid<0||fid>=M||!F[fid].active) continue; Face f=F[fid]; if(hasv(f,a)&&hasv(f,b)){
            if(cnt>=2) return false; ef[cnt]=fid; for(int k=0;k<3;k++){int x=f.v[k]; if(x!=a&&x!=b){op[cnt]=x; break;}} cnt++;
        }}
    return cnt==2 && op[0]!=op[1];
}
static bool link_ok(int a,int b){
    int ef[2],op[2]; if(!collect_edge(a,b,ef,op)) return false;
    if(++tokenV>2000000000){ fill(MarkV.begin(),MarkV.end(),0); tokenV=3; }
    int ta=tokenV++;
    for(int fid:Inc[a]) if(F[fid].active&&hasv(F[fid],a)) for(int k=0;k<3;k++){ int x=F[fid].v[k]; if(x!=a&&x!=b) MarkV[x]=ta; }
    int cnt=0,g0=0,g1=0; int tc=tokenV++;
    for(int fid:Inc[b]) if(F[fid].active&&hasv(F[fid],b)) for(int k=0;k<3;k++){
        int x=F[fid].v[k]; if(x==a||x==b) continue; if(MarkV[x]!=ta) continue;
        if(x!=op[0]&&x!=op[1]) return false;
        if(MarkV[x]!=tc){ MarkV[x]=tc; cnt++; if(x==op[0])g0=1; if(x==op[1])g1=1; if(cnt>2) return false; }
    }
    return cnt==2&&g0&&g1;
}
static void compact_inc(int u){
    if(u<0||u>=N||!Alive[u]) return; vector<int>&L=Inc[u]; int w=0;
    for(int fid:L) if(0<=fid&&fid<M&&F[fid].active&&hasv(F[fid],u)) L[w++]=fid;
    L.resize(w);
}
static bool cluster_fits(int u,const Vec3&pos){ return Radius[u]+sqrt(dist2(P[u],pos)) <= epsLen; }
static bool geometry_ok(int src,int dst,const Vec3&pos,double normalCos){
    if(!finite3(pos)||Locked[src]) return false;
    if(Locked[dst] && dist2(P[dst],pos)>max(1e-30,diagLen*diagLen*1e-30)) return false;
    if(!cluster_fits(src,pos)||!cluster_fits(dst,pos)) return false;
    if(++tokenF>2000000000){ fill(MarkF.begin(),MarkF.end(),0); tokenF=5; }
    int tf=tokenF++; vector<int>aff; aff.reserve(Inc[src].size()+Inc[dst].size());
    for(int fid:Inc[src]) if(F[fid].active&&hasv(F[fid],src)){ MarkF[fid]=tf; aff.push_back(fid); }
    for(int fid:Inc[dst]) if(F[fid].active&&hasv(F[fid],dst)&&MarkF[fid]!=tf){ MarkF[fid]=tf; aff.push_back(fid); }
    unordered_set<unsigned long long> local; local.reserve(aff.size()*2+4); vector<array<int,3>> newTris; newTris.reserve(aff.size());
    for(int fid:aff){
        Face f=F[fid]; bool hs=hasv(f,src), hd=hasv(f,dst); if(hs&&hd) continue;
        Vec3 oldp[3]={P[f.v[0]],P[f.v[1]],P[f.v[2]]}; Vec3 newp[3]={oldp[0],oldp[1],oldp[2]}; int idx[3]={f.v[0],f.v[1],f.v[2]};
        for(int k=0;k<3;k++) if(idx[k]==src||idx[k]==dst){ idx[k]=dst; newp[k]=pos; }
        if(idx[0]==idx[1]||idx[0]==idx[2]||idx[1]==idx[2]) return false;
        Vec3 on=cross3(oldp[1]-oldp[0],oldp[2]-oldp[0]); Vec3 nn=cross3(newp[1]-newp[0],newp[2]-newp[0]); double ol2=norm2(on),nl2=norm2(nn);
        if(ol2<=areaEps2||nl2<=areaEps2) return false;
        if(dot3(on,nn) < normalCos*sqrt(ol2*nl2)) return false;
        auto ts=tri_sorted(idx[0],idx[1],idx[2]); unsigned long long h=pack_tri(ts); if(local.find(h)!=local.end()) return false; local.insert(h); newTris.push_back(ts);
    }
    for(auto ts:newTris){ int probe=ts[0]; for(int fid:Inc[probe]) if(F[fid].active&&MarkF[fid]!=tf){ Face f=F[fid]; if(tri_sorted(f.v[0],f.v[1],f.v[2])==ts) return false; } }
    return true;
}
static double cheap_cost(int a,int b){
    Quadric q=Q[a]; q.add(Q[b]); Vec3 p; double best=1e100; if(solve3(q,p)) best=min(best,q.eval(p)); best=min(best,q.eval((P[a]+P[b])*0.5)); best=min(best,q.eval(P[a])); best=min(best,q.eval(P[b]));
    return best + 1e-7*dist2(P[a],P[b])/(diagLen*diagLen+1e-30);
}
static void push_edge(int a,int b){
    if(a==b||a<0||b<0||a>=N||b>=N||!Alive[a]||!Alive[b]) return;
    if(Locked[a]&&Locked[b]) return;
    if(dist2(P[a],P[b])>4.00001*eps2) return;
    PQ.push({cheap_cost(a,b),a,b,Version[a],Version[b]});
}
static bool best_candidate(int a,int b,int&src,int&dst,Vec3&bestp,double normalCos){
    if(!link_ok(a,b)) return false;
    if(Locked[a]&&!Locked[b]){ dst=a; src=b; }
    else if(Locked[b]&&!Locked[a]){ dst=b; src=a; }
    else if(Inc[a].size()+2ull*ClusterSize[a] <= Inc[b].size()+2ull*ClusterSize[b]){ src=a; dst=b; }
    else{ src=b; dst=a; }
    Quadric q=Q[a]; q.add(Q[b]); Vec3 cand[8]; int n=0; Vec3 opt;
    if(Locked[dst]){ cand[n++]=P[dst]; }
    else{
        if(solve3(q,opt)) cand[n++]=opt;
        cand[n++]=(P[a]+P[b])*0.5; cand[n++]=P[dst]; cand[n++]=P[src]; cand[n++]=P[dst]*0.85+P[src]*0.15; cand[n++]=P[dst]*0.70+P[src]*0.30;
    }
    double best=1e100; bool ok=false;
    for(int i=0;i<n;i++){
        Vec3 p=cand[i]; if(!geometry_ok(src,dst,p,normalCos)) continue;
        double c=q.eval(p)+1e-6*(dist2(p,P[a])+dist2(p,P[b]))/(diagLen*diagLen+1e-30);
        if(i==2) c*=0.985; // prefer an existing endpoint when nearly tied
        if(c<best){ best=c; bestp=p; ok=true; }
    }
    return ok;
}
static void merge_cluster(int src,int dst,const Vec3&oldSrc,const Vec3&oldDst,const Vec3&pos){
    double nr=max(Radius[src]+sqrt(dist2(oldSrc,pos)), Radius[dst]+sqrt(dist2(oldDst,pos)));
    Radius[dst]=nr;
    if(Head[src]!=-1){ NextMem[Tail[dst]]=Head[src]; Tail[dst]=Tail[src]; ClusterSize[dst]+=ClusterSize[src]; Locked[dst]|=Locked[src]; Head[src]=Tail[src]=-1; ClusterSize[src]=0; }
}
static void do_collapse(int src,int dst,const Vec3&pos){
    Vec3 oldSrc=P[src], oldDst=P[dst]; Q[dst].add(Q[src]); P[dst]=pos;
    for(int fid:Inc[src]) if(F[fid].active){ Face&f=F[fid]; bool hs=hasv(f,src), hd=hasv(f,dst); if(!hs) continue; if(hd){ f.active=0; --activeF; } else { for(int k=0;k<3;k++) if(f.v[k]==src) f.v[k]=dst; Inc[dst].push_back(fid); } }
    Alive[src]=0; --activeV; ++Version[src]; ++Version[dst]; merge_cluster(src,dst,oldSrc,oldDst,pos); compact_inc(dst); Inc[src].clear();
    for(int fid:Inc[dst]) if(F[fid].active){ Face f=F[fid]; push_edge(f.v[0],f.v[1]); push_edge(f.v[1],f.v[2]); push_edge(f.v[2],f.v[0]); }
}
static bool attempt(int a,int b,double normalCos){ int src,dst; Vec3 pos; if(!best_candidate(a,b,src,dst,pos,normalCos)) return false; do_collapse(src,dst,pos); return true; }

static double initialize_quadrics(){
    vector<Vec3> fn(M); vector<EdgeRef> er; er.reserve((size_t)M*3);
    for(int i=0;i<M;i++){
        Face&f=F[i]; if(f.v[0]<0||f.v[1]<0||f.v[2]<0||f.v[0]>=N||f.v[1]>=N||f.v[2]>=N){f.active=0; continue;}
        Vec3 cr=cross3(P[f.v[1]]-P[f.v[0]],P[f.v[2]]-P[f.v[0]]); double len=norm3(cr); if(len<=0){f.active=0; continue;} Vec3 n=cr/len; fn[i]=n; double d=-dot3(n,P[f.v[0]]); double w=max(1e-12,len/diagLen);
        Q[f.v[0]].add_plane(n.x,n.y,n.z,d,w); Q[f.v[1]].add_plane(n.x,n.y,n.z,d,w); Q[f.v[2]].add_plane(n.x,n.y,n.z,d,w);
        er.push_back({ekey(f.v[0],f.v[1]),i}); er.push_back({ekey(f.v[1],f.v[2]),i}); er.push_back({ekey(f.v[2],f.v[0]),i});
    }
    sort(er.begin(),er.end()); long long uni=0,feat=0; double fcos=cos(35.0*M_PI/180.0);
    for(size_t i=0;i<er.size();){ size_t j=i+1; while(j<er.size()&&er[j].key==er[i].key) j++; ++uni; int a=ka(er[i].key),b=kb(er[i].key); push_edge(a,b); if(j-i==2 && dot3(fn[er[i].f],fn[er[i+1].f])<fcos) ++feat; i=j; }
    return uni?double(feat)/double(uni):0.0;
}
static void run_qem(){
    double feature=initialize_quadrics();
    double ratio=0.090+min(0.060,feature*0.24);
    if(N<1200) ratio=max(ratio,0.22); else if(N<6000) ratio=max(ratio,0.145); else if(N<25000) ratio=max(ratio,0.112);
    if(feature>0.18) ratio=max(ratio,0.135); ratio=min(ratio,0.20);
    int target=max(8,(int)ceil(N*ratio)); double normalCos=feature>0.18?0.50:0.30; if(N>200000) normalCos=min(normalCos,0.24);
    long long pops=0, stale=0;
    while(activeV>target && !PQ.empty()){
        if((++pops&4095)==0 && elapsed()>19.2) break;
        Node nd=PQ.top(); PQ.pop(); int a=nd.a,b=nd.b;
        if(a==b||a<0||b<0||a>=N||b>=N||!Alive[a]||!Alive[b]) continue;
        if(nd.va!=Version[a]||nd.vb!=Version[b]){ if(++stale<4000000) push_edge(a,b); continue; }
        attempt(a,b,normalCos);
    }
}
static bool build_current(vector<Vec3>&X,vector<Face>&A){
    vector<int> remap(N,-1); X.clear(); X.reserve(activeV);
    for(int i=0;i<N;i++) if(Alive[i]){ remap[i]=(int)X.size(); X.push_back(P[i]); }
    A.clear(); A.reserve(activeF); unordered_set<unsigned long long> seen; seen.reserve((size_t)activeF*2+16);
    for(int i=0;i<M;i++) if(F[i].active){ int a=F[i].v[0],b=F[i].v[1],c=F[i].v[2]; if(a<0||b<0||c<0||a>=N||b>=N||c>=N) continue; int ra=remap[a],rb=remap[b],rc=remap[c]; if(ra<0||rb<0||rc<0||ra==rb||ra==rc||rb==rc) continue; auto ts=tri_sorted(ra,rb,rc); unsigned long long h=pack_tri(ts); if(seen.insert(h).second){ Face f{{ra,rb,rc},1}; A.push_back(f); } }
    return validate_mesh_arrays(X,A,true);
}

int main(){
    T0=chrono::steady_clock::now(); read_input();
    if(N<=0||M<=0){ return 0; }
    // High-upside safe special cases first. They print only after exact manifold + vertex-cover validation.
    if(try_snapped_ellipsoid()) return 0;
    if(try_regular_grid()) return 0;
    run_qem();
    vector<Vec3>X; vector<Face>A;
    if(build_current(X,A)) print_mesh(X,A); else print_original();
    return 0;
}
