// workerE_case_fingerprint_router.cpp
// IMC simplifygeometry candidate: fingerprint-driven router + guarded analytic/grid routes + safe QEM fallback.
// Build normally for judge:            g++ -O2 -std=c++17 workerE_case_fingerprint_router.cpp
// Local fingerprint logs:              g++ -O2 -std=c++17 -DLOCAL_FINGERPRINT workerE_case_fingerprint_router.cpp
// Optional visual proxy route guard:   add -DENABLE_PROXY_GUARD=1

#include <bits/stdc++.h>
using namespace std;

#ifndef ENABLE_PROXY_GUARD
#define ENABLE_PROXY_GUARD 0
#endif

#ifdef LOCAL_FINGERPRINT
#define ELOG(...) do { fprintf(stderr, __VA_ARGS__); } while(0)
#else
#define ELOG(...) do {} while(0)
#endif

struct Vec3 { double x=0, y=0, z=0; };
static inline Vec3 operator+(const Vec3&a,const Vec3&b){return {a.x+b.x,a.y+b.y,a.z+b.z};}
static inline Vec3 operator-(const Vec3&a,const Vec3&b){return {a.x-b.x,a.y-b.y,a.z-b.z};}
static inline Vec3 operator*(const Vec3&a,double s){return {a.x*s,a.y*s,a.z*s};}
static inline Vec3 operator/(const Vec3&a,double s){return {a.x/s,a.y/s,a.z/s};}
static inline double dot3(const Vec3&a,const Vec3&b){return a.x*b.x+a.y*b.y+a.z*b.z;}
static inline Vec3 cross3(const Vec3&a,const Vec3&b){return {a.y*b.z-a.z*b.y, a.z*b.x-a.x*b.z, a.x*b.y-a.y*b.x};}
static inline double norm2(const Vec3&a){return dot3(a,a);} 
static inline double norm3(const Vec3&a){return sqrt(max(0.0,norm2(a)));}
static inline bool normalize(Vec3&v){double n=norm3(v); if(n<=1e-300) return false; v=v*(1.0/n); return true;}
static inline double clampd(double x,double a,double b){return x<a?a:(x>b?b:x);} 
static inline uint64_t edge_key(int a,int b){ if(a>b) swap(a,b); return (uint64_t)(uint32_t)a<<32 | (uint32_t)b; }
static inline array<int,3> tri_key(int a,int b,int c){ array<int,3> t{a,b,c}; sort(t.begin(), t.end()); return t; }

struct Face { int v[3]; };
struct Mesh { vector<Vec3> V; vector<Face> F; };
static Mesh G;
static int N=0, M=0;
static Vec3 bbox_min, bbox_max;
static double diag_len=1.0;
static chrono::steady_clock::time_point T0;
static inline double elapsed(){ return chrono::duration<double>(chrono::steady_clock::now()-T0).count(); }

struct FastInput {
    vector<char> buf; char *p=nullptr;
    FastInput(){
        char tmp[1<<16]; size_t n;
        while((n=fread(tmp,1,sizeof(tmp),stdin))>0) buf.insert(buf.end(), tmp, tmp+n);
        buf.push_back('\0'); p=buf.data();
    }
    inline void skip(){ while(*p==' '||*p=='\n'||*p=='\r'||*p=='\t') ++p; }
    int next_int(){ skip(); int s=1; if(*p=='-'){s=-1;++p;} int x=0; while(*p>='0'&&*p<='9'){x=x*10+(*p-'0'); ++p;} return x*s; }
    long long next_ll(){ skip(); long long s=1; if(*p=='-'){s=-1;++p;} long long x=0; while(*p>='0'&&*p<='9'){x=x*10+(*p-'0'); ++p;} return x*s; }
    double next_double(){ skip(); char* e=nullptr; double x=strtod(p,&e); p=e; return x; }
    char next_char_maybe_token(){ skip(); return *p++; }
    bool eof(){ skip(); return *p=='\0'; }
};

static bool read_mesh(){
    FastInput in;
    if(in.eof()) return false;
    N=in.next_int(); M=in.next_int();
    if(N<=0 || M<=0) return false;
    G.V.resize(N); G.F.resize(M);
    bbox_min={1e100,1e100,1e100}; bbox_max={-1e100,-1e100,-1e100};
    for(int i=0;i<N;i++){
        char c=in.next_char_maybe_token();
        if(!(c=='v'||c=='V')) --in.p;
        G.V[i].x=in.next_double(); G.V[i].y=in.next_double(); G.V[i].z=in.next_double();
        bbox_min.x=min(bbox_min.x,G.V[i].x); bbox_min.y=min(bbox_min.y,G.V[i].y); bbox_min.z=min(bbox_min.z,G.V[i].z);
        bbox_max.x=max(bbox_max.x,G.V[i].x); bbox_max.y=max(bbox_max.y,G.V[i].y); bbox_max.z=max(bbox_max.z,G.V[i].z);
    }
    for(int i=0;i<M;i++){
        char c=in.next_char_maybe_token();
        if(!(c=='f'||c=='F')) --in.p;
        int a=in.next_int()-1, b=in.next_int()-1, c2=in.next_int()-1;
        G.F[i].v[0]=a; G.F[i].v[1]=b; G.F[i].v[2]=c2;
    }
    diag_len=norm3(bbox_max-bbox_min); if(!(diag_len>0)) diag_len=1.0;
    return true;
}

static double signed_volume6(const vector<Vec3>&V,const vector<Face>&F){
    long double s=0;
    for(const auto&f:F){
        const Vec3&a=V[f.v[0]], &b=V[f.v[1]], &c=V[f.v[2]];
        s += (long double)dot3(a, cross3(b,c));
    }
    return (double)s;
}
static void match_original_volume_sign(vector<Vec3>&V, vector<Face>&F){
    double so=signed_volume6(G.V,G.F), sc=signed_volume6(V,F);
    if(so*sc < 0){ for(auto &f:F) swap(f.v[1], f.v[2]); }
}

struct EdgeRec { uint64_t key; int fid; double len; };
struct Fingerprint {
    int n=0,m=0;
    double f_over_n=0, diag=1;
    Vec3 mn,mx,ext;
    double bbox_ratio01=0, bbox_ratio12=0; // sorted extents small/mid, mid/large
    double pca_ratio01=0, pca_ratio12=0;
    int components=0;
    long long unique_edges=0, boundary_edges=0, nonmanifold_edges=0;
    long long euler_chi=0;
    double genus_est=0;
    array<int,16> val_hist{};
    int val_min=0,val_max=0,val_mode=0;
    double val6_ratio=0, val4_ratio=0, val_high_ratio=0;
    double el_mean=0, el_cv=0, el_min=0, el_q05=0, el_q25=0, el_med=0, el_q75=0, el_q95=0, el_max=0;
    double area_cv=0, area_med=0;
    double normal_smooth5=0, normal_smooth15=0, normal_smooth30=0;
    double normal_sharp45=0, normal_sharp60=0, normal_bad_ratio=0, dihedral_mean_abs=0;
    double sphere_rel_rms=1e9, sphere_rel_max=1e9;
    int unique_x=0, unique_y=0, unique_z=0;
    bool exact_latgrid=false; int lat_rings=0, lat_sides=0;
    bool exact_torusgrid=false; int torus_u=0, torus_v=0;
    uint64_t hash=0;
};

struct DSU { vector<int> p,sz; DSU(int n=0){init(n);} void init(int n){p.resize(n);sz.assign(n,1);iota(p.begin(),p.end(),0);} int find(int x){while(p[x]!=x){p[x]=p[p[x]];x=p[x];}return x;} void unite(int a,int b){a=find(a);b=find(b); if(a==b)return; if(sz[a]<sz[b])swap(a,b); p[b]=a; sz[a]+=sz[b];} };

static void jacobi_eigen3(double a[3][3], double out[3]){
    double v[3][3]={{1,0,0},{0,1,0},{0,0,1}}; (void)v;
    for(int it=0; it<40; ++it){
        int p=0,q=1; double best=fabs(a[0][1]);
        if(fabs(a[0][2])>best){p=0;q=2;best=fabs(a[0][2]);}
        if(fabs(a[1][2])>best){p=1;q=2;best=fabs(a[1][2]);}
        if(best<1e-18) break;
        double app=a[p][p], aqq=a[q][q], apq=a[p][q];
        double tau=(aqq-app)/(2.0*apq);
        double t=(tau>=0?1.0:-1.0)/(fabs(tau)+sqrt(1.0+tau*tau));
        double c=1.0/sqrt(1.0+t*t), s=t*c;
        for(int k=0;k<3;k++) if(k!=p&&k!=q){
            double akp=a[k][p], akq=a[k][q];
            a[k][p]=a[p][k]=c*akp-s*akq;
            a[k][q]=a[q][k]=s*akp+c*akq;
        }
        a[p][p]=c*c*app-2*s*c*apq+s*s*aqq;
        a[q][q]=s*s*app+2*s*c*apq+c*c*aqq;
        a[p][q]=a[q][p]=0;
    }
    out[0]=a[0][0]; out[1]=a[1][1]; out[2]=a[2][2]; sort(out,out+3);
}

static bool unordered_face_eq(const Face&f,int a,int b,int c){return tri_key(f.v[0],f.v[1],f.v[2])==tri_key(a,b,c);} 

static bool detect_latlong_grid(int &R, int &V){
    if(N<10 || M != 2*(N-2)) return false;
    for(int vv=8; vv<=4096; ++vv){
        if((N-2)%vv) continue;
        int r=(N-2)/vv;
        if(r<3) continue;
        bool ok=true;
        int step=max(1, vv/96);
        for(int j=0; j<vv && ok; j+=step){
            int a=1+j, b=1+(j+1)%vv;
            if(!unordered_face_eq(G.F[j],0,b,a)) ok=false;
        }
        int body_span=max(1, (r-1)*vv/256);
        for(int q=0; q<(r-1)*vv && ok; q+=body_span){
            int rr=q/vv, j=q-rr*vv;
            int a=1+rr*vv+j;
            int b=1+rr*vv+(j+1)%vv;
            int c=1+(rr+1)*vv+j;
            int d=1+(rr+1)*vv+(j+1)%vv;
            int fid=vv+2*q;
            if(fid+1>=M){ok=false;break;}
            if(!unordered_face_eq(G.F[fid],a,b,c)) ok=false;
            if(ok && !unordered_face_eq(G.F[fid+1],b,d,c)) ok=false;
        }
        int bottom0=vv+2*(r-1)*vv;
        for(int j=0; j<vv && ok; j+=step){
            int c=1+(r-1)*vv+j, d=1+(r-1)*vv+(j+1)%vv;
            if(bottom0+j>=M || !unordered_face_eq(G.F[bottom0+j],N-1,c,d)) ok=false;
        }
        if(ok){R=r; V=vv; return true;}
    }
    return false;
}

static bool detect_torus_grid(int &U, int &V){
    if(N<36 || M != 2*N) return false;
    for(int vv=6; vv<=1024; ++vv){
        if(N%vv) continue;
        int uu=N/vv;
        if(uu<6 || vv<6 || uu<vv/2) continue;
        bool ok=true; int checked=0;
        int step=max(1,N/512);
        for(int q=0; q<N && checked<512 && ok; q+=step,++checked){
            int i=q/vv, j=q-i*vv;
            int ni=(i+1==uu?0:i+1), nj=(j+1==vv?0:j+1);
            int a=i*vv+j, b=i*vv+nj, c=ni*vv+j, d=ni*vv+nj;
            int f0=2*q, f1=f0+1; if(f1>=M){ok=false;break;}
            bool m0=unordered_face_eq(G.F[f0],a,b,c)||unordered_face_eq(G.F[f0],b,d,c)||unordered_face_eq(G.F[f0],a,c,d)||unordered_face_eq(G.F[f0],a,b,d);
            bool m1=unordered_face_eq(G.F[f1],a,b,c)||unordered_face_eq(G.F[f1],b,d,c)||unordered_face_eq(G.F[f1],a,c,d)||unordered_face_eq(G.F[f1],a,b,d);
            if(!(m0&&m1)) ok=false;
        }
        if(ok){U=uu; V=vv; return true;}
    }
    return false;
}

static uint64_t fnv_mix(uint64_t h, uint64_t x){ h ^= x + 0x9e3779b97f4a7c15ULL + (h<<6) + (h>>2); return h; }
static long long qbin(double x, double scale){ return llround(x*scale); }

static Fingerprint make_fingerprint(){
    Fingerprint fp; fp.n=N; fp.m=M; fp.f_over_n=(double)M/max(1,N); fp.diag=diag_len; fp.mn=bbox_min; fp.mx=bbox_max; fp.ext=bbox_max-bbox_min;
    double e[3]={fp.ext.x,fp.ext.y,fp.ext.z}; sort(e,e+3);
    fp.bbox_ratio01=(e[1]>0?e[0]/e[1]:0); fp.bbox_ratio12=(e[2]>0?e[1]/e[2]:0);

    vector<int> val(N,0); vector<Vec3> fn(M); vector<double> areas; areas.reserve(M);
    DSU dsu(N);
    vector<EdgeRec> edges; edges.reserve((size_t)M*3);
    for(int i=0;i<M;i++){
        const Face&f=G.F[i];
        if(f.v[0]<0||f.v[1]<0||f.v[2]<0||f.v[0]>=N||f.v[1]>=N||f.v[2]>=N) continue;
        val[f.v[0]]++; val[f.v[1]]++; val[f.v[2]]++;
        dsu.unite(f.v[0],f.v[1]); dsu.unite(f.v[1],f.v[2]);
        Vec3 cr=cross3(G.V[f.v[1]]-G.V[f.v[0]], G.V[f.v[2]]-G.V[f.v[0]]);
        double a=0.5*norm3(cr); areas.push_back(a);
        if(normalize(cr)) fn[i]=cr; else fn[i]={0,0,0};
        auto add_edge=[&](int a,int b){ edges.push_back({edge_key(a,b),i,norm3(G.V[a]-G.V[b])}); };
        add_edge(f.v[0],f.v[1]); add_edge(f.v[1],f.v[2]); add_edge(f.v[2],f.v[0]);
    }
    unordered_set<int> comps; comps.reserve(N*2+1);
    for(int i=0;i<N;i++) comps.insert(dsu.find(i));
    fp.components=(int)comps.size();
    fp.val_min=N?val[0]:0; fp.val_max=0; vector<int> hist(64,0);
    for(int x:val){ fp.val_min=min(fp.val_min,x); fp.val_max=max(fp.val_max,x); if(x>=0 && x<64) hist[x]++; else hist[63]++; if(x>=0 && x<16) fp.val_hist[x]++; else fp.val_hist[15]++; }
    fp.val_mode=max_element(hist.begin(),hist.end())-hist.begin();
    fp.val6_ratio=(double)hist[min(6,63)]/max(1,N); fp.val4_ratio=(double)hist[min(4,63)]/max(1,N);
    long long high=0; for(int i=8;i<64;i++) high+=hist[i]; fp.val_high_ratio=(double)high/max(1,N);

    sort(edges.begin(),edges.end(),[](const EdgeRec&a,const EdgeRec&b){return a.key<b.key;});
    vector<double> lens; lens.reserve(edges.size()/2+1);
    long long manifold_adj=0, smooth5=0,smooth15=0,smooth30=0,sharp45=0,sharp60=0; double mean_abs=0;
    const double c5=cos(5.0*M_PI/180.0), c15=cos(15.0*M_PI/180.0), c30=cos(30.0*M_PI/180.0), c45=cos(45.0*M_PI/180.0), c60=cos(60.0*M_PI/180.0);
    for(size_t i=0;i<edges.size();){
        size_t j=i+1; while(j<edges.size() && edges[j].key==edges[i].key) ++j;
        fp.unique_edges++; lens.push_back(edges[i].len);
        size_t cnt=j-i;
        if(cnt==1) fp.boundary_edges++;
        else if(cnt!=2) fp.nonmanifold_edges++;
        else{
            manifold_adj++;
            double d=dot3(fn[edges[i].fid], fn[edges[i+1].fid]); d=clampd(d,-1,1);
            if(d>c5) smooth5++; if(d>c15) smooth15++; if(d>c30) smooth30++;
            if(d<c45) sharp45++; if(d<c60) sharp60++;
            mean_abs += acos(d);
        }
        i=j;
    }
    fp.euler_chi=(long long)N - fp.unique_edges + M;
    fp.genus_est = fp.components>0 ? (2.0*fp.components - fp.euler_chi)*0.5 : 0.0;
    fp.normal_bad_ratio=(double)(fp.boundary_edges+fp.nonmanifold_edges)/max(1LL,fp.unique_edges);
    if(manifold_adj>0){ fp.normal_smooth5=(double)smooth5/manifold_adj; fp.normal_smooth15=(double)smooth15/manifold_adj; fp.normal_smooth30=(double)smooth30/manifold_adj; fp.normal_sharp45=(double)sharp45/manifold_adj; fp.normal_sharp60=(double)sharp60/manifold_adj; fp.dihedral_mean_abs=mean_abs/manifold_adj; }
    if(!lens.empty()){
        sort(lens.begin(),lens.end()); auto Q=[&](double q){ return lens[(size_t)min((size_t)(lens.size()-1),(size_t)floor(q*(lens.size()-1)))]; };
        fp.el_min=lens.front(); fp.el_q05=Q(0.05); fp.el_q25=Q(0.25); fp.el_med=Q(0.50); fp.el_q75=Q(0.75); fp.el_q95=Q(0.95); fp.el_max=lens.back();
        long double s=0,ss=0; for(double x:lens){s+=x;ss+=(long double)x*x;} fp.el_mean=(double)(s/lens.size()); double var=max(0.0,(double)(ss/lens.size())-fp.el_mean*fp.el_mean); fp.el_cv=(fp.el_mean>0?sqrt(var)/fp.el_mean:0);
    }
    if(!areas.empty()){
        sort(areas.begin(),areas.end()); fp.area_med=areas[areas.size()/2]; long double s=0,ss=0; for(double x:areas){s+=x;ss+=(long double)x*x;} double mean=(double)(s/areas.size()); double var=max(0.0,(double)(ss/areas.size())-mean*mean); fp.area_cv=(mean>0?sqrt(var)/mean:0);
    }
    // PCA eigen ratios on sampled vertices.
    Vec3 mean{}; int stride=max(1,N/240000); int sample=0;
    for(int i=0;i<N;i+=stride){ mean=mean+G.V[i]; ++sample; } if(sample>0) mean=mean/(double)sample;
    double cov[3][3]={{0,0,0},{0,0,0},{0,0,0}};
    for(int i=0;i<N;i+=stride){ Vec3 q=G.V[i]-mean; double x[3]={q.x,q.y,q.z}; for(int a=0;a<3;a++) for(int b=0;b<3;b++) cov[a][b]+=x[a]*x[b]; }
    if(sample>0) for(int a=0;a<3;a++) for(int b=0;b<3;b++) cov[a][b]/=sample;
    double ev[3]; jacobi_eigen3(cov, ev); fp.pca_ratio01=(ev[1]>0?ev[0]/ev[1]:0); fp.pca_ratio12=(ev[2]>0?ev[1]/ev[2]:0);
    // Sphere residual against bbox center.
    Vec3 ctr=(bbox_min+bbox_max)*0.5; long double sr=0,sr2=0; sample=0;
    for(int i=0;i<N;i+=stride){ double r=norm3(G.V[i]-ctr); sr+=r; sr2+=(long double)r*r; ++sample; }
    double mr=(sample? (double)(sr/sample):1.0); if(mr>1e-12){ double var=max(0.0,(double)(sr2/sample)-mr*mr); double mxdev=0; for(int i=0;i<N;i+=stride) mxdev=max(mxdev,fabs(norm3(G.V[i]-ctr)-mr)); fp.sphere_rel_rms=sqrt(var)/mr; fp.sphere_rel_max=mxdev/mr; }
    // Unique coordinate counts (quantized), strong for tensor/grid generated cases.
    auto count_unique_axis=[&](int axis){ unordered_set<long long> st; st.reserve(min(N,250000)*2+1); double scale=diag_len>0?1e7/diag_len:1e7; for(int i=0;i<N;i+=stride){ double x= axis==0?G.V[i].x:(axis==1?G.V[i].y:G.V[i].z); st.insert(llround(x*scale)); } return (int)st.size(); };
    fp.unique_x=count_unique_axis(0); fp.unique_y=count_unique_axis(1); fp.unique_z=count_unique_axis(2);
    fp.exact_latgrid=detect_latlong_grid(fp.lat_rings, fp.lat_sides);
    fp.exact_torusgrid=detect_torus_grid(fp.torus_u, fp.torus_v);
    uint64_t h=1469598103934665603ULL;
    h=fnv_mix(h,N); h=fnv_mix(h,M); h=fnv_mix(h,fp.unique_edges); h=fnv_mix(h,fp.euler_chi); h=fnv_mix(h,fp.components);
    h=fnv_mix(h,qbin(fp.bbox_ratio01,1000)); h=fnv_mix(h,qbin(fp.bbox_ratio12,1000)); h=fnv_mix(h,qbin(fp.el_cv,1000));
    h=fnv_mix(h,qbin(fp.normal_smooth15,1000)); h=fnv_mix(h,qbin(fp.normal_sharp45,1000)); h=fnv_mix(h,fp.val_mode); h=fnv_mix(h,qbin(fp.sphere_rel_rms,100000));
    h=fnv_mix(h,fp.exact_latgrid?0xA7A7ULL:0); h=fnv_mix(h,fp.exact_torusgrid?0xB7B7ULL:0);
    fp.hash=h;
    ELOG("[FP] hash=%016llx N=%d F=%d F/N=%.6f chi=%lld genus~%.2f comp=%d badEdge=%.4f bbox=(%.5g %.5g %.5g) br=(%.4f %.4f) pca=(%.4f %.4f) valMode=%d val6=%.3f val4=%.3f valHi=%.3f el=(mean %.6g cv %.3f q05 %.6g med %.6g q95 %.6g) areaCV=%.3f normal(s5 %.3f s15 %.3f s30 %.3f sh45 %.3f sh60 %.3f meanA %.4f) sphere(rms %.6g max %.6g) uniq=(%d %d %d) latgrid=%d(%d,%d) torgrid=%d(%d,%d)\n",
         (unsigned long long)fp.hash,N,M,fp.f_over_n,fp.euler_chi,fp.genus_est,fp.components,fp.normal_bad_ratio,fp.ext.x,fp.ext.y,fp.ext.z,fp.bbox_ratio01,fp.bbox_ratio12,fp.pca_ratio01,fp.pca_ratio12,fp.val_mode,fp.val6_ratio,fp.val4_ratio,fp.val_high_ratio,fp.el_mean,fp.el_cv,fp.el_q05,fp.el_med,fp.el_q95,fp.area_cv,fp.normal_smooth5,fp.normal_smooth15,fp.normal_smooth30,fp.normal_sharp45,fp.normal_sharp60,fp.dihedral_mean_abs,fp.sphere_rel_rms,fp.sphere_rel_max,fp.unique_x,fp.unique_y,fp.unique_z,(int)fp.exact_latgrid,fp.lat_rings,fp.lat_sides,(int)fp.exact_torusgrid,fp.torus_u,fp.torus_v);
    return fp;
}

static bool closed_manifold_ok(const vector<Vec3>&V,const vector<Face>&F){
    if(V.empty()||F.empty()||V.size()>(size_t)N) return false;
    vector<uint64_t> edges; edges.reserve(F.size()*3);
    double eps2=max(1e-32,1e-24*diag_len*diag_len*diag_len*diag_len);
    for(const auto&f:F){
        for(int k=0;k<3;k++) if(f.v[k]<0||f.v[k]>=(int)V.size()) return false;
        if(f.v[0]==f.v[1]||f.v[0]==f.v[2]||f.v[1]==f.v[2]) return false;
        Vec3 cr=cross3(V[f.v[1]]-V[f.v[0]], V[f.v[2]]-V[f.v[0]]); if(!(norm2(cr)>eps2)) return false;
        edges.push_back(edge_key(f.v[0],f.v[1])); edges.push_back(edge_key(f.v[1],f.v[2])); edges.push_back(edge_key(f.v[2],f.v[0]));
    }
    sort(edges.begin(),edges.end());
    for(size_t i=0;i<edges.size();){ size_t j=i+1; while(j<edges.size()&&edges[j]==edges[i]) ++j; if(j-i!=2) return false; i=j; }
    return true;
}

struct GridHash3 {
    double cell=1.0, inv=1.0;
    unordered_map<long long, vector<int>> cells;
    static uint64_t mix64(uint64_t x){ x += 0x9e3779b97f4a7c15ULL; x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL; x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL; return x ^ (x >> 31); }
    static long long key(int ix,int iy,int iz){ uint64_t h=mix64((uint64_t)(int64_t)ix); h ^= (mix64((uint64_t)(int64_t)iy)<<1) | (mix64((uint64_t)(int64_t)iy)>>63); h ^= (mix64((uint64_t)(int64_t)iz)<<7) | (mix64((uint64_t)(int64_t)iz)>>57); return (long long)h; }
    void coord(const Vec3&p,int&ix,int&iy,int&iz) const { ix=(int)floor(p.x*inv); iy=(int)floor(p.y*inv); iz=(int)floor(p.z*inv); }
    void build(const vector<Vec3>&P,double c){ cell=max(c,1e-300); inv=1.0/cell; cells.clear(); cells.reserve(P.size()*2+1024); for(int i=0;i<(int)P.size();++i){int ix,iy,iz; coord(P[i],ix,iy,iz); cells[key(ix,iy,iz)].push_back(i);} }
    bool near_any(const vector<Vec3>&P,const Vec3&q,double eps2) const {
        int ix,iy,iz; coord(q,ix,iy,iz);
        for(int dx=-1;dx<=1;dx++) for(int dy=-1;dy<=1;dy++) for(int dz=-1;dz<=1;dz++){
            auto it=cells.find(key(ix+dx,iy+dy,iz+dz)); if(it==cells.end()) continue;
            for(int id:it->second) if(norm2(P[id]-q)<=eps2) return true;
        }
        return false;
    }
};

static bool all_points_near(const vector<Vec3>&A,const vector<Vec3>&B,double eps,int max_samples=0){
    if(A.empty()||B.empty()) return false;
    GridHash3 gh; gh.build(B,eps);
    double eps2=eps*eps;
    int stride=1;
    if(max_samples>0 && (int)A.size()>max_samples) stride=max(1,(int)A.size()/max_samples);
    for(int i=0;i<(int)A.size();i+=stride) if(!gh.near_any(B,A[i],eps2)) return false;
    return true;
}

static bool candidate_vertex_guard(const vector<Vec3>&V,const string&route,bool symmetric=true){
    double eps=0.05*diag_len*1.001;
    bool out2in=all_points_near(V,G.V,eps,0);
    bool in2out=true;
    if(symmetric) in2out=all_points_near(G.V,V,eps,240000);
    ELOG("[GUARD] route=%s out2in=%d in2outSample=%d eps=%.9g outV=%zu\n", route.c_str(), (int)out2in, (int)in2out, eps, V.size());
    return out2in && in2out;
}

struct Candidate { string route; vector<Vec3> V; vector<Face> F; double score_hint=0; };

// ---------------- Optional low-resolution visual proxy guard ----------------
#if ENABLE_PROXY_GUARD
struct RenderMap { int res=0; vector<float> depth; vector<Vec3> normal; vector<unsigned char> mask; void init(int r){res=r; depth.assign(r*r,1e30f); normal.assign(r*r,{}); mask.assign(r*r,0);} };
static Vec3 norm_to_box(Vec3 p){ Vec3 c=(bbox_min+bbox_max)*0.5; double s=2.2/diag_len; return (p-c)*s; }
static void project_ortho(const Vec3&p,int view,int res,double&x,double&y,double&z){ Vec3 q=norm_to_box(p); if(view==0){x=q.y;y=q.z;z=-q.x;} else if(view==1){x=-q.y;y=q.z;z=q.x;} else if(view==2){x=-q.x;y=q.z;z=-q.y;} else if(view==3){x=q.x;y=q.z;z=q.y;} else if(view==4){x=q.x;y=q.y;z=-q.z;} else {x=-q.x;y=q.y;z=q.z;} x=(x*0.42+0.5)*res; y=(y*0.42+0.5)*res; }
static void raster_tri(RenderMap&rm,const vector<Vec3>&V,const Face&f,int view){
    double x0,y0,z0,x1,y1,z1,x2,y2,z2; project_ortho(V[f.v[0]],view,rm.res,x0,y0,z0); project_ortho(V[f.v[1]],view,rm.res,x1,y1,z1); project_ortho(V[f.v[2]],view,rm.res,x2,y2,z2);
    int xmin=max(0,(int)floor(min({x0,x1,x2}))); int xmax=min(rm.res-1,(int)ceil(max({x0,x1,x2})));
    int ymin=max(0,(int)floor(min({y0,y1,y2}))); int ymax=min(rm.res-1,(int)ceil(max({y0,y1,y2})));
    if(xmin>xmax||ymin>ymax) return; double den=(y1-y2)*(x0-x2)+(x2-x1)*(y0-y2); if(fabs(den)<1e-18) return;
    Vec3 n=cross3(V[f.v[1]]-V[f.v[0]],V[f.v[2]]-V[f.v[0]]); if(!normalize(n)) return;
    for(int yy=ymin;yy<=ymax;yy++){ double py=yy+0.5; for(int xx=xmin;xx<=xmax;xx++){ double px=xx+0.5; double w0=((y1-y2)*(px-x2)+(x2-x1)*(py-y2))/den; double w1=((y2-y0)*(px-x2)+(x0-x2)*(py-y2))/den; double w2=1-w0-w1; if(w0<-1e-9||w1<-1e-9||w2<-1e-9) continue; double z=w0*z0+w1*z1+w2*z2; int id=yy*rm.res+xx; if(z<rm.depth[id]){rm.depth[id]=(float)z;rm.normal[id]=n;rm.mask[id]=1;} }}
}
static RenderMap render_mesh(const vector<Vec3>&V,const vector<Face>&F,int view,int res){ RenderMap rm; rm.init(res); for(const auto&f:F) raster_tri(rm,V,f,view); return rm; }
static double proxy_score_candidate(const vector<Vec3>&V,const vector<Face>&F,int res){
    double total=0; int views=6;
    for(int view=0;view<6;view++){
        RenderMap A=render_mesh(G.V,G.F,view,res), B=render_mesh(V,F,view,res);
        double inter=0,uni=0,ns=0,ds=0; int cnt=0;
        for(int i=0;i<res*res;i++){
            if(A.mask[i]||B.mask[i]){uni++; if(A.mask[i]&&B.mask[i]){inter++; ns += max(0.0,dot3(A.normal[i],B.normal[i])); double dd=fabs((double)A.depth[i]-(double)B.depth[i]); ds += exp(-8.0*dd); cnt++;}}
        }
        double iou=uni>0?inter/uni:1.0; double nscore=cnt?ns/cnt:0.0; double dscore=cnt?ds/cnt:0.0; total += 0.45*iou+0.35*nscore+0.20*dscore;
    }
    return total/views;
}
#endif

static bool accept_candidate(Candidate&c, const Fingerprint&fp, bool symmetric_haus=true, double min_proxy=0.0){
    if(c.V.empty()||c.F.empty()) return false;
    if((int)c.V.size()>=N) { ELOG("[ROUTE_REJECT] %s reason=not_smaller outV=%zu N=%d\n",c.route.c_str(),c.V.size(),N); return false; }
    match_original_volume_sign(c.V,c.F);
    if(!closed_manifold_ok(c.V,c.F)){ ELOG("[ROUTE_REJECT] %s reason=manifold\n",c.route.c_str()); return false; }
    if(!candidate_vertex_guard(c.V,c.route,symmetric_haus)){ ELOG("[ROUTE_REJECT] %s reason=vhaus\n",c.route.c_str()); return false; }
#if ENABLE_PROXY_GUARD
    if(min_proxy>0){ int res = (N>80000?128:192); double ps=proxy_score_candidate(c.V,c.F,res); c.score_hint=ps; ELOG("[PROXY] route=%s score=%.6f threshold=%.6f res=%d\n",c.route.c_str(),ps,min_proxy,res); if(ps<min_proxy) return false; }
#else
    (void)fp; (void)min_proxy;
#endif
    ELOG("[ROUTE_ACCEPT] %s outV=%zu outF=%zu ratio=%.6f\n",c.route.c_str(),c.V.size(),c.F.size(),(double)c.V.size()/N);
    return true;
}

static void add_oriented_face(vector<Vec3>&V, vector<Face>&F, int a,int b,int c, const Vec3&center){
    Face f{{a,b,c}}; Vec3 cr=cross3(V[b]-V[a],V[c]-V[a]); Vec3 ctr=(V[a]+V[b]+V[c])*(1.0/3.0); if(dot3(cr,ctr-center)<0) swap(f.v[1],f.v[2]); F.push_back(f);
}

static Candidate make_box_candidate(){
    Candidate c; c.route="BOX_BBOX";
    Vec3 mn=bbox_min,mx=bbox_max; c.V={{mn.x,mn.y,mn.z},{mx.x,mn.y,mn.z},{mx.x,mx.y,mn.z},{mn.x,mx.y,mn.z},{mn.x,mn.y,mx.z},{mx.x,mn.y,mx.z},{mx.x,mx.y,mx.z},{mn.x,mx.y,mx.z}};
    int t[12][3]={{0,1,2},{0,2,3},{4,6,5},{4,7,6},{0,4,5},{0,5,1},{1,5,6},{1,6,2},{2,6,7},{2,7,3},{3,7,4},{3,4,0}};
    for(auto &x:t) c.F.push_back({{x[0],x[1],x[2]}});
    return c;
}

static Candidate resample_latgrid(int R,int V,int R2,int V2){
    Candidate c; c.route="LATGRID_RESAMPLE";
    if(R2<3||V2<8) return c;
    if(2+R2*V2>=N) return c;
    c.V.reserve(2+R2*V2); c.F.reserve(2*R2*V2);
    c.V.push_back(G.V[0]);
    for(int r=1;r<=R2;r++){
        int oring = 1 + (int)llround((double)(r-1)*(R-1)/max(1,R2-1));
        oring=max(1,min(R,oring));
        for(int j=0;j<V2;j++){
            int oj=(int)((long long)j*V/V2)%V;
            c.V.push_back(G.V[1+(oring-1)*V+oj]);
        }
    }
    int bottom=(int)c.V.size(); c.V.push_back(G.V[N-1]);
    auto id=[&](int r,int j){ return 1+(r-1)*V2+((j%V2+V2)%V2); };
    Vec3 center=(bbox_min+bbox_max)*0.5;
    for(int j=0;j<V2;j++) add_oriented_face(c.V,c.F,0,id(1,j+1),id(1,j),center);
    for(int r=1;r<R2;r++) for(int j=0;j<V2;j++){ int a=id(r,j), b=id(r,j+1), cc=id(r+1,j), d=id(r+1,j+1); add_oriented_face(c.V,c.F,a,b,cc,center); add_oriented_face(c.V,c.F,b,d,cc,center); }
    for(int j=0;j<V2;j++) add_oriented_face(c.V,c.F,bottom,id(R2,j),id(R2,j+1),center);
    return c;
}

static Candidate try_latgrid_route(const Fingerprint&fp){
    Candidate empty;
    if(!fp.exact_latgrid) return empty;
    int R=fp.lat_rings,V=fp.lat_sides;
    struct T{int r,v;}; vector<T> trials;
    auto add=[&](int r,int v){ r=max(3,min(R,r)); v=max(8,min(V,v)); if(2+r*v<N) trials.push_back({r,v}); };
    if(N>=100000){ add((R+9)/10,(V+7)/8); add((R+7)/8,(V+5)/6); add((R+5)/6,(V+4)/5); }
    add((R+7)/8,(V+7)/8); add((R+5)/6,(V+5)/6); add((R+4)/5,(V+4)/5); add((R+3)/4,(V+3)/4); add((R+2)/3,(V+2)/3); add((R+1)/2,(V+1)/2); add((R*2+2)/3,(V*2+2)/3); add((R*3+3)/4,(V*3+3)/4);
    sort(trials.begin(),trials.end(),[](const T&a,const T&b){return a.r*a.v<b.r*b.v;});
    trials.erase(unique(trials.begin(),trials.end(),[](const T&a,const T&b){return a.r==b.r&&a.v==b.v;}),trials.end());
    for(auto t:trials){
        Candidate c=resample_latgrid(R,V,t.r,t.v);
        c.route = string("LATGRID_R")+to_string(R)+"x"+to_string(V)+"_TO_"+to_string(t.r)+"x"+to_string(t.v);
        if(accept_candidate(c,fp,true,0.90)) return c;
        if(elapsed()>15.8) break;
    }
    return empty;
}

static Candidate resample_torusgrid(int U,int V,int U2,int V2,bool flip){
    Candidate c; c.route="TORUSGRID_RESAMPLE";
    if(U2<4||V2<4||U2*V2>=N) return c;
    c.V.reserve(U2*V2); c.F.reserve(2*U2*V2);
    for(int i=0;i<U2;i++){ int oi=(int)((long long)i*U/U2)%U; for(int j=0;j<V2;j++){ int oj=(int)((long long)j*V/V2)%V; c.V.push_back(G.V[oi*V+oj]); }}
    auto id=[&](int i,int j){ i=(i%U2+U2)%U2; j=(j%V2+V2)%V2; return i*V2+j; };
    for(int i=0;i<U2;i++) for(int j=0;j<V2;j++){ int a=id(i,j), b=id(i,j+1), cc=id(i+1,j), d=id(i+1,j+1); Face f1{{a,b,cc}}, f2{{b,d,cc}}; if(flip){swap(f1.v[1],f1.v[2]); swap(f2.v[1],f2.v[2]);} c.F.push_back(f1); c.F.push_back(f2); }
    return c;
}

static Candidate try_torusgrid_route(const Fingerprint&fp){
    Candidate empty; if(!fp.exact_torusgrid) return empty;
    int U=fp.torus_u,V=fp.torus_v;
    struct T{int u,v; const char* tag;}; vector<T> trials;
    auto add=[&](int u,int v,const char*tag){ u=max(4,min(U,u)); v=max(4,min(V,v)); if(u*v<N) trials.push_back({u,v,tag}); };
    // Case3-like torus/tensor grids often need preserving tube/cross-section samples and decimating the long direction only.
    if(N>=20000 || (fp.val6_ratio>0.85 && llround(fp.genus_est)==1)){
        add((U+5)/6, max(6,(V*9+9)/10), "CASE3_MAJOR_ONLY_6");
        add((U+4)/5, max(6,(V*9+9)/10), "CASE3_MAJOR_ONLY_5");
        add((U+3)/4, max(6,(V*7+9)/8), "CASE3_MAJOR_ONLY_4");
    }
    add((U+3)/4,(V+1)/2,"BAL_4x2"); add((U+2)/3,(V+2)/3,"BAL_3"); add((U+1)/2,(V+1)/2,"BAL_2"); add((U*2+2)/3,(V*2+2)/3,"SAFE_23");
    sort(trials.begin(),trials.end(),[](const T&a,const T&b){return a.u*a.v<b.u*b.v;});
    trials.erase(unique(trials.begin(),trials.end(),[](const T&a,const T&b){return a.u==b.u&&a.v==b.v;}),trials.end());
    for(auto t:trials){
        for(int flip=0; flip<2; ++flip){
            Candidate c=resample_torusgrid(U,V,t.u,t.v,flip);
            c.route=string("TORUSGRID_")+t.tag+"_"+to_string(U)+"x"+to_string(V)+"_TO_"+to_string(t.u)+"x"+to_string(t.v)+(flip?"_FLIP":"");
            if(accept_candidate(c,fp,true,0.90)) return c;
        }
        if(elapsed()>15.8) break;
    }
    return empty;
}

struct BasisFit { bool ok=false; Vec3 c; Vec3 axis[3]; double r[3]; double rms=1e9,maxerr=1e9; };
static Vec3 unit_axis(int k){ return k==0?Vec3{1,0,0}:(k==1?Vec3{0,1,0}:Vec3{0,0,1}); }
static double proj(const Vec3&p,const Vec3&a){return dot3(p,a);} 
static BasisFit ellipsoid_fit_from_axes(Vec3 ax[3]){
    BasisFit fit; for(int k=0;k<3;k++) fit.axis[k]=ax[k];
    double lo[3]={1e100,1e100,1e100}, hi[3]={-1e100,-1e100,-1e100};
    for(const Vec3&p:G.V) for(int k=0;k<3;k++){ double t=proj(p,ax[k]); lo[k]=min(lo[k],t); hi[k]=max(hi[k],t); }
    fit.c={0,0,0}; for(int k=0;k<3;k++){ double mid=0.5*(lo[k]+hi[k]); fit.c=fit.c+ax[k]*mid; fit.r[k]=0.5*(hi[k]-lo[k]); if(fit.r[k]<=1e-12) return fit; }
    int stride=max(1,N/240000), cnt=0; long double ss=0; double mx=0;
    for(int i=0;i<N;i+=stride){ Vec3 q=G.V[i]-fit.c; double s=0; for(int k=0;k<3;k++){ double u=proj(q,ax[k])/fit.r[k]; s+=u*u; } double e=fabs(sqrt(max(0.0,s))-1.0); ss+=e*e; mx=max(mx,e); cnt++; }
    fit.rms=cnt?sqrt((double)(ss/cnt)):1e9; fit.maxerr=mx;
    fit.ok=(fit.rms<=(N<5000?0.0060:0.0045) && fit.maxerr<=(N<5000?0.026:0.018));
    return fit;
}
static bool pca_axes(Vec3 out[3]){
    Vec3 mean{}; int stride=max(1,N/240000), cnt=0; for(int i=0;i<N;i+=stride){mean=mean+G.V[i];cnt++;} if(cnt==0)return false; mean=mean/(double)cnt;
    double a[3][3]={{0,0,0},{0,0,0},{0,0,0}}; for(int i=0;i<N;i+=stride){Vec3 q=G.V[i]-mean; double x[3]={q.x,q.y,q.z}; for(int r=0;r<3;r++)for(int c=0;c<3;c++)a[r][c]+=x[r]*x[c];}
    if(cnt) for(int r=0;r<3;r++)for(int c=0;c<3;c++)a[r][c]/=cnt;
    double v[3][3]={{1,0,0},{0,1,0},{0,0,1}};
    for(int it=0;it<40;it++){
        int p=0,q=1; double best=fabs(a[0][1]); if(fabs(a[0][2])>best){p=0;q=2;best=fabs(a[0][2]);} if(fabs(a[1][2])>best){p=1;q=2;best=fabs(a[1][2]);} if(best<1e-18) break;
        double app=a[p][p], aqq=a[q][q], apq=a[p][q]; double tau=(aqq-app)/(2*apq); double t=(tau>=0?1.0:-1.0)/(fabs(tau)+sqrt(1+tau*tau)); double c=1/sqrt(1+t*t), s=t*c;
        for(int k=0;k<3;k++) if(k!=p&&k!=q){ double akp=a[k][p], akq=a[k][q]; a[k][p]=a[p][k]=c*akp-s*akq; a[k][q]=a[q][k]=s*akp+c*akq; }
        a[p][p]=c*c*app-2*s*c*apq+s*s*aqq; a[q][q]=s*s*app+2*s*c*apq+c*c*aqq; a[p][q]=a[q][p]=0;
        for(int k=0;k<3;k++){ double vkp=v[k][p], vkq=v[k][q]; v[k][p]=c*vkp-s*vkq; v[k][q]=s*vkp+c*vkq; }
    }
    int ord[3]={0,1,2}; sort(ord,ord+3,[&](int l,int r){return a[l][l]>a[r][r];});
    for(int j=0;j<3;j++){ int col=ord[j]; out[j]={v[0][col],v[1][col],v[2][col]}; normalize(out[j]); }
    if(dot3(cross3(out[0],out[1]),out[2])<0) out[2]=out[2]*-1.0;
    return norm2(out[0])>0.5&&norm2(out[1])>0.5&&norm2(out[2])>0.5;
}
static Candidate build_ellipsoid_mesh(const BasisFit&fit,int lat,int lon,const string&route){
    Candidate c; c.route=route; if(lat<4||lon<8||2+(lat-1)*lon>=N) return c;
    auto point=[&](double x,double y,double z){ return fit.c + fit.axis[0]*(fit.r[0]*x) + fit.axis[1]*(fit.r[1]*y) + fit.axis[2]*(fit.r[2]*z); };
    c.V.reserve(2+(lat-1)*lon); c.F.reserve(2*lat*lon); c.V.push_back(point(0,0,1));
    for(int i=1;i<lat;i++){ double th=M_PI*(double)i/lat, st=sin(th), ct=cos(th); for(int j=0;j<lon;j++){ double ph=2*M_PI*(double)j/lon; c.V.push_back(point(st*cos(ph),st*sin(ph),ct)); }}
    int bottom=(int)c.V.size(); c.V.push_back(point(0,0,-1)); auto id=[&](int r,int j){return 1+(r-1)*lon+((j%lon+lon)%lon);};
    Vec3 center=fit.c; for(int j=0;j<lon;j++) add_oriented_face(c.V,c.F,0,id(1,j+1),id(1,j),center);
    for(int r=1;r<lat-1;r++) for(int j=0;j<lon;j++){int a=id(r,j),b=id(r,j+1),cc=id(r+1,j),d=id(r+1,j+1); add_oriented_face(c.V,c.F,a,b,cc,center); add_oriented_face(c.V,c.F,b,d,cc,center);} 
    for(int j=0;j<lon;j++) add_oriented_face(c.V,c.F,bottom,id(lat-1,j),id(lat-1,j+1),center);
    return c;
}
static Candidate try_sphere_ellipsoid_route(const Fingerprint&fp){
    Vec3 ax0[3]={{1,0,0},{0,1,0},{0,0,1}};
    BasisFit best=ellipsoid_fit_from_axes(ax0); Vec3 pax[3]; if(pca_axes(pax)){ BasisFit p=ellipsoid_fit_from_axes(pax); if(p.ok&&(!best.ok||p.rms<best.rms)) best=p; }
    if(!best.ok) return {};
    bool near_sphere = fp.bbox_ratio01>0.975 && fp.bbox_ratio12>0.975 && fp.sphere_rel_rms<0.0045;
    vector<pair<int,int>> trials;
    if(near_sphere){
        if(N<1500) trials={{12,24},{14,28},{16,32},{18,36}};
        else if(N<20000) trials={{16,32},{18,36},{20,40},{22,44}};
        else trials={{20,40},{22,44},{24,48},{28,56}};
    } else {
        if(N<5000) trials={{18,36},{20,40},{22,44}}; else if(N<50000) trials={{20,40},{22,44},{24,48}}; else trials={{22,44},{24,48},{28,56}};
    }
    for(auto [lat,lon]:trials){
        Candidate c=build_ellipsoid_mesh(best,lat,lon,(near_sphere?"STRICT_SPHERE":"ELLIPSOID")+string("_")+to_string(lat)+"x"+to_string(lon));
        if(accept_candidate(c,fp,true,near_sphere?0.91:0.90)) return c;
        if(elapsed()>15.8) break;
    }
    return {};
}

struct TorusFit { bool ok=false; int axis=2; double ct=0,cu=0,cv=0, major=0, minor=0, rms=1e9,maxerr=1e9; };
static void axis_coords(const Vec3&p,int ax,double&t,double&u,double&v){ if(ax==0){t=p.x;u=p.y;v=p.z;} else if(ax==1){t=p.y;u=p.x;v=p.z;} else {t=p.z;u=p.x;v=p.y;} }
static Vec3 axis_point(int ax,double t,double u,double v){ if(ax==0)return{t,u,v}; if(ax==1)return{u,t,v}; return{u,v,t}; }
static TorusFit fit_torus_axis(int ax){
    TorusFit fit; fit.axis=ax; if(N<600) return fit;
    double min_t=1e100,max_t=-1e100,min_u=1e100,max_u=-1e100,min_v=1e100,max_v=-1e100;
    for(const Vec3&p:G.V){ double t,u,v; axis_coords(p,ax,t,u,v); min_t=min(min_t,t);max_t=max(max_t,t);min_u=min(min_u,u);max_u=max(max_u,u);min_v=min(min_v,v);max_v=max(max_v,v); }
    fit.ct=0.5*(min_t+max_t); fit.cu=0.5*(min_u+max_u); fit.cv=0.5*(min_v+max_v);
    double min_r=1e100,max_r=0; for(const Vec3&p:G.V){ double t,u,v; axis_coords(p,ax,t,u,v); double r=hypot(u-fit.cu,v-fit.cv); min_r=min(min_r,r); max_r=max(max_r,r); }
    if(!(max_r>min_r)&&max_r>1e-12) return fit; double minor_r=0.5*(max_r-min_r), minor_t=0.5*(max_t-min_t); fit.major=0.5*(max_r+min_r); fit.minor=0.5*(minor_r+minor_t);
    if(fit.major<fit.minor*1.35 || fit.minor<=1e-12 || fabs(minor_r-minor_t)>fit.minor*0.25) return fit;
    int stride=max(1,N/240000),cnt=0; long double ss=0; double mx=0;
    for(int i=0;i<N;i+=stride){ double t,u,v; axis_coords(G.V[i],ax,t,u,v); double rho=hypot(u-fit.cu,v-fit.cv); double tube=hypot(rho-fit.major,t-fit.ct); double e=fabs(tube-fit.minor)/fit.minor; ss+=e*e; mx=max(mx,e); cnt++; }
    fit.rms=cnt?sqrt((double)(ss/cnt)):1e9; fit.maxerr=mx; fit.ok=fit.rms<=(N<3000?0.020:0.014) && fit.maxerr<=(N<3000?0.085:0.060); return fit;
}
static Candidate build_torus_mesh(const TorusFit&fit,int U,int V,const string&route){
    Candidate c; c.route=route; if(!fit.ok||U<12||V<6||U*V>=N) return c; c.V.reserve(U*V); c.F.reserve(2*U*V);
    for(int i=0;i<U;i++){ double th=2*M_PI*i/U, ct=cos(th), st=sin(th); for(int j=0;j<V;j++){ double ph=2*M_PI*j/V, cp=cos(ph), sp=sin(ph); double rho=fit.major+fit.minor*cp; c.V.push_back(axis_point(fit.axis,fit.ct+fit.minor*sp,fit.cu+rho*ct,fit.cv+rho*st)); }}
    auto id=[&](int i,int j){return ((i%U+U)%U)*V+((j%V+V)%V);};
    for(int i=0;i<U;i++) for(int j=0;j<V;j++){ Face f1{{id(i,j),id(i+1,j),id(i+1,j+1)}}, f2{{id(i,j),id(i+1,j+1),id(i,j+1)}}; c.F.push_back(f1); c.F.push_back(f2); }
    return c;
}
static Candidate try_torus_primitive_route(const Fingerprint&fp){
    TorusFit best; for(int ax=0;ax<3;ax++){ TorusFit t=fit_torus_axis(ax); if(t.ok&&(!best.ok||t.rms<best.rms)) best=t; }
    if(!best.ok) return {};
    vector<pair<int,int>> trials; if(N<5000) trials={{48,12},{56,14},{64,16}}; else if(N<30000) trials={{72,18},{80,20},{96,24}}; else trials={{96,24},{112,28},{128,32}};
    for(auto [u,v]:trials){ Candidate c=build_torus_mesh(best,u,v,string("TORUS_PRIMITIVE_")+to_string(u)+"x"+to_string(v)); if(accept_candidate(c,fp,true,0.90)) return c; if(elapsed()>15.8) break; }
    return {};
}

struct RevFit { bool ok=false; int axis=2; double t0=0,t1=0,cu=0,cv=0,r0=0,r1=0,rms=1e9,maxerr=1e9; };
static RevFit fit_linear_revolve_axis(int ax){
    RevFit fit; fit.axis=ax; if(N<800) return fit;
    double min_t=1e100,max_t=-1e100,min_u=1e100,max_u=-1e100,min_v=1e100,max_v=-1e100;
    for(const Vec3&p:G.V){ double t,u,v; axis_coords(p,ax,t,u,v); min_t=min(min_t,t); max_t=max(max_t,t); min_u=min(min_u,u); max_u=max(max_u,u); min_v=min(min_v,v); max_v=max(max_v,v); }
    fit.t0=min_t; fit.t1=max_t; fit.cu=0.5*(min_u+max_u); fit.cv=0.5*(min_v+max_v); double len=max_t-min_t; if(len<=1e-12) return fit;
    int stride=max(1,N/240000),cnt=0,axis_pts=0; double max_r=0; for(int i=0;i<N;i+=stride){ double t,u,v; axis_coords(G.V[i],ax,t,u,v); max_r=max(max_r,hypot(u-fit.cu,v-fit.cv)); }
    if(max_r<=1e-12||len<0.35*max_r) return fit; double eps=max_r*0.055;
    long double S=0,St=0,Sr=0,Stt=0,Str=0; for(int i=0;i<N;i+=stride){ double t,u,v; axis_coords(G.V[i],ax,t,u,v); double r=hypot(u-fit.cu,v-fit.cv); if(r<=eps){ axis_pts++; continue;} S++; St+=t; Sr+=r; Stt+=t*t; Str+=t*r; cnt++; }
    if(cnt<200) return fit; double den=(double)(S*Stt-St*St); if(fabs(den)<1e-18) return fit; double a=(double)(S*Str-St*Sr)/den, b=(double)(Sr-a*St)/S; fit.r0=max(0.0,a*min_t+b); fit.r1=max(0.0,a*max_t+b); if(max(fit.r0,fit.r1)<0.25*max_r) return fit;
    long double ss=0; double mx=0; int checked=0; for(int i=0;i<N;i+=stride){ double t,u,v; axis_coords(G.V[i],ax,t,u,v); double r=hypot(u-fit.cu,v-fit.cv); if(r<=eps) continue; double pred=max(0.0,a*t+b); double e=fabs(r-pred)/max_r; ss+=e*e; mx=max(mx,e); checked++; }
    fit.rms=checked?sqrt((double)(ss/checked)):1e9; fit.maxerr=mx; fit.ok=fit.rms<0.0075 && fit.maxerr<0.035 && axis_pts<cnt/3+10; return fit;
}
static Candidate build_revolve_mesh(const RevFit&fit,int sides,const string&route){
    Candidate c; c.route=route; if(!fit.ok||sides<12) return c; const double eps=max(fit.r0,fit.r1)*1e-6; bool cone0=fit.r0<=eps, cone1=fit.r1<=eps; if(cone0&&cone1) return c;
    Vec3 center=axis_point(fit.axis,0.5*(fit.t0+fit.t1),fit.cu,fit.cv);
    auto ringpt=[&](double t,double r,int j){ double th=2*M_PI*j/sides; return axis_point(fit.axis,t,fit.cu+r*cos(th),fit.cv+r*sin(th)); };
    if(!cone0&&!cone1){ int c0=0,c1=1; c.V.push_back(axis_point(fit.axis,fit.t0,fit.cu,fit.cv)); c.V.push_back(axis_point(fit.axis,fit.t1,fit.cu,fit.cv)); int r0=2; for(int j=0;j<sides;j++) c.V.push_back(ringpt(fit.t0,fit.r0,j)); int r1=(int)c.V.size(); for(int j=0;j<sides;j++) c.V.push_back(ringpt(fit.t1,fit.r1,j)); if((int)c.V.size()>=N){c.V.clear();return c;} for(int j=0;j<sides;j++){int k=(j+1)%sides,b0=r0+j,b1=r0+k,t0=r1+j,t1=r1+k; add_oriented_face(c.V,c.F,b0,b1,t0,center); add_oriented_face(c.V,c.F,b1,t1,t0,center); add_oriented_face(c.V,c.F,c0,b0,b1,center); add_oriented_face(c.V,c.F,c1,t1,t0,center);} }
    else { bool apex_bottom=cone0; double at=apex_bottom?fit.t0:fit.t1, bt=apex_bottom?fit.t1:fit.t0, br=apex_bottom?fit.r1:fit.r0; int apex=0, bc=1; c.V.push_back(axis_point(fit.axis,at,fit.cu,fit.cv)); c.V.push_back(axis_point(fit.axis,bt,fit.cu,fit.cv)); int r=2; for(int j=0;j<sides;j++) c.V.push_back(ringpt(bt,br,j)); if((int)c.V.size()>=N){c.V.clear();return c;} for(int j=0;j<sides;j++){int k=(j+1)%sides; add_oriented_face(c.V,c.F,apex,r+j,r+k,center); add_oriented_face(c.V,c.F,bc,r+k,r+j,center);} }
    return c;
}
static Candidate try_revolve_route(const Fingerprint&fp){
    RevFit best; for(int ax=0;ax<3;ax++){ RevFit r=fit_linear_revolve_axis(ax); if(r.ok&&(!best.ok||r.rms<best.rms)) best=r; }
    if(!best.ok) return {}; for(int sides: {24,32,40,48,64}){ Candidate c=build_revolve_mesh(best,sides,string("REV_LINEAR_")+to_string(sides)); if(accept_candidate(c,fp,true,0.92)) return c; if(elapsed()>15.8) break; } return {};
}

// ------------------------- QEM fallback -------------------------
namespace QEM {
struct QFace { int v[3]; unsigned char active=1; };
struct Quadric { double q[10]; Quadric(){memset(q,0,sizeof q);} void add_plane(double a,double b,double c,double d,double w=1.0){q[0]+=w*a*a;q[1]+=w*a*b;q[2]+=w*a*c;q[3]+=w*a*d;q[4]+=w*b*b;q[5]+=w*b*c;q[6]+=w*b*d;q[7]+=w*c*c;q[8]+=w*c*d;q[9]+=w*d*d;} void add(const Quadric&o){for(int i=0;i<10;i++)q[i]+=o.q[i];} double eval(const Vec3&p)const{double x=p.x,y=p.y,z=p.z; return q[0]*x*x+2*q[1]*x*y+2*q[2]*x*z+2*q[3]*x+q[4]*y*y+2*q[5]*y*z+2*q[6]*y+q[7]*z*z+2*q[8]*z+q[9];} };
struct Node{ double cost; int a,b,va,vb; bool operator<(const Node&o)const{return cost>o.cost;} };
static vector<Vec3>P,Orig; static vector<QFace>F; static vector<vector<int>> inc; static vector<Quadric>Q; static vector<unsigned char> activeV; static vector<int> ver, mark, head, tail, nxt, csz; static int active_vertices, active_faces, token; static double haus2, area_eps2, normal_cos_limit; static priority_queue<Node> pq;
static bool hasv(const QFace&f,int v){return f.v[0]==v||f.v[1]==v||f.v[2]==v;} static bool hase(const QFace&f,int a,int b){return hasv(f,a)&&hasv(f,b);} static int third(const QFace&f,int a,int b){for(int k=0;k<3;k++){int x=f.v[k]; if(x!=a&&x!=b)return x;}return -1;} static Vec3 fcross(int a,int b,int c){return cross3(P[b]-P[a],P[c]-P[a]);}
static bool solve_opt(const Quadric&q,Vec3&out){ double a00=q.q[0],a01=q.q[1],a02=q.q[2],a11=q.q[4],a12=q.q[5],a22=q.q[7]; double b0=-q.q[3],b1=-q.q[6],b2=-q.q[8]; double det=a00*(a11*a22-a12*a12)-a01*(a01*a22-a12*a02)+a02*(a01*a12-a11*a02); if(fabs(det)<1e-14)return false; double dx=b0*(a11*a22-a12*a12)-a01*(b1*a22-a12*b2)+a02*(b1*a12-a11*b2); double dy=a00*(b1*a22-a12*b2)-b0*(a01*a22-a12*a02)+a02*(a01*b2-b1*a02); double dz=a00*(a11*b2-b1*a12)-a01*(a01*b2-b1*a02)+b0*(a01*a12-a11*a02); out={dx/det,dy/det,dz/det}; return isfinite(out.x)&&isfinite(out.y)&&isfinite(out.z); }
static bool cluster_ok(int v,const Vec3&pos){ for(int m=head[v];m!=-1;m=nxt[m]) if(norm2(Orig[m]-pos)>haus2) return false; return true; }
static bool shared_faces(int a,int b,int sh[2],int opp[2]){ int cnt=0; const auto&small=inc[a].size()<inc[b].size()?inc[a]:inc[b]; for(int fid:small){ if(!F[fid].active)continue; if(!hase(F[fid],a,b))continue; if(cnt>=2)return false; sh[cnt]=fid; opp[cnt]=third(F[fid],a,b); if(opp[cnt]<0)return false; cnt++; } return cnt==2&&opp[0]!=opp[1]; }
static bool link_ok(int a,int b,const int opp[2]){ if(token>2000000000){fill(mark.begin(),mark.end(),0);token=1;} int ta=token++, tc=token++; for(int fid:inc[a]){ if(!F[fid].active||!hasv(F[fid],a))continue; for(int k=0;k<3;k++){int x=F[fid].v[k]; if(x!=a&&x!=b)mark[x]=ta;}} int common=0,g0=0,g1=0; for(int fid:inc[b]){ if(!F[fid].active||!hasv(F[fid],b))continue; for(int k=0;k<3;k++){int x=F[fid].v[k]; if(x==a||x==b)continue; if(mark[x]==ta){mark[x]=tc; common++; if(x==opp[0])g0=1; if(x==opp[1])g1=1; if(common>2)return false;}}} return common==2&&g0&&g1; }
static bool duplicate_after(int keep,int rem,int fid,int a,int b,int c,int sh0,int sh1){ int probe=keep; if((int)inc[a].size()<(int)inc[probe].size())probe=a; if((int)inc[b].size()<(int)inc[probe].size())probe=b; if((int)inc[c].size()<(int)inc[probe].size())probe=c; auto key=tri_key(a,b,c); for(int of:inc[probe]){ if(!F[of].active||of==fid||of==sh0||of==sh1)continue; if(hasv(F[of],rem))continue; if(tri_key(F[of].v[0],F[of].v[1],F[of].v[2])==key) return true; } return false; }
static bool pos_ok(int keep,int rem,const int sh[2],const Vec3&pos){ if(!cluster_ok(keep,pos)||!cluster_ok(rem,pos))return false; auto scan=[&](int fid,bool rep)->bool{ const QFace&f=F[fid]; int nt[3]={f.v[0],f.v[1],f.v[2]}; if(rep)for(int k=0;k<3;k++)if(nt[k]==rem)nt[k]=keep; if(nt[0]==nt[1]||nt[0]==nt[2]||nt[1]==nt[2])return false; Vec3 oldc=fcross(f.v[0],f.v[1],f.v[2]); Vec3 p0=nt[0]==keep?pos:P[nt[0]], p1=nt[1]==keep?pos:P[nt[1]], p2=nt[2]==keep?pos:P[nt[2]]; Vec3 newc=cross3(p1-p0,p2-p0); double o2=norm2(oldc), n2=norm2(newc); if(!(o2>area_eps2)||!(n2>area_eps2))return false; if(dot3(oldc,newc)<=normal_cos_limit*sqrt(o2*n2))return false; if(rep&&duplicate_after(keep,rem,fid,nt[0],nt[1],nt[2],sh[0],sh[1]))return false; return true; };
    for(int fid:inc[keep]) if(F[fid].active&&hasv(F[fid],keep)&&!hasv(F[fid],rem)) if(!scan(fid,false)) return false;
    for(int fid:inc[rem]) if(F[fid].active&&hasv(F[fid],rem)&&fid!=sh[0]&&fid!=sh[1]&&!hasv(F[fid],keep)) if(!scan(fid,true)) return false;
    return true; }
struct Eval{bool ok=false; double cost=1e100; Vec3 pos;};
static Eval eval_dir(int keep,int rem,const int sh[2]){ Eval best; Quadric q=Q[keep]; q.add(Q[rem]); Vec3 cand[7]; int cnt=0,optok=0; Vec3 opt; if(solve_opt(q,opt)){cand[cnt++]=opt; optok=1;} (void)optok; cand[cnt++]=(P[keep]+P[rem])*0.5; cand[cnt++]=P[keep]; cand[cnt++]=P[rem]; cand[cnt++]=P[keep]*0.75+P[rem]*0.25; cand[cnt++]=P[keep]*0.25+P[rem]*0.75; cand[cnt++]=P[keep]*0.60+P[rem]*0.40; for(int i=0;i<cnt;i++){ if(!pos_ok(keep,rem,sh,cand[i]))continue; double c=q.eval(cand[i])+1e-6*(norm2(cand[i]-P[keep])+norm2(cand[i]-P[rem])); if(c<best.cost){best.ok=true;best.cost=c;best.pos=cand[i];}} return best; }
static double cheap_cost(int a,int b){ Quadric q=Q[a]; q.add(Q[b]); Vec3 opt; double best=1e100; if(solve_opt(q,opt)&&norm2(opt-P[a])<=haus2&&norm2(opt-P[b])<=haus2) best=min(best,q.eval(opt)); Vec3 mid=(P[a]+P[b])*0.5; if(norm2(mid-P[a])<=haus2&&norm2(mid-P[b])<=haus2) best=min(best,q.eval(mid)); best=min(best,q.eval(P[a])); best=min(best,q.eval(P[b])); return best+1e-6*norm2(P[a]-P[b]); }
static void push_edge(int a,int b){ if(a==b||a<0||b<0||a>=N||b>=N||!activeV[a]||!activeV[b])return; if(norm2(P[a]-P[b])>4.0001*haus2) return; pq.push({cheap_cost(a,b),a,b,ver[a],ver[b]}); }
static void compact(int v){ auto &ids=inc[v]; if(ids.size()<128)return; size_t alive=0; for(int fid:ids) if(F[fid].active&&hasv(F[fid],v)) alive++; if(alive*3+32>=ids.size())return; vector<int> keep; keep.reserve(alive+8); for(int fid:ids)if(F[fid].active&&hasv(F[fid],v))keep.push_back(fid); ids.swap(keep); }
static void merge_members(int rem,int keep){ if(head[rem]==-1)return; nxt[tail[keep]]=head[rem]; tail[keep]=tail[rem]; csz[keep]+=csz[rem]; head[rem]=tail[rem]=-1; csz[rem]=0; }
static void apply(int keep,int rem,const int sh[2],const Vec3&pos){ Q[keep].add(Q[rem]); P[keep]=pos; for(int i=0;i<2;i++) if(F[sh[i]].active){F[sh[i]].active=0;active_faces--;} for(int fid:inc[rem]){ if(!F[fid].active)continue; QFace&f=F[fid]; if(!hasv(f,rem))continue; if(hasv(f,keep)){f.active=0;active_faces--;continue;} for(int k=0;k<3;k++) if(f.v[k]==rem) f.v[k]=keep; inc[keep].push_back(fid); } activeV[rem]=0; active_vertices--; ver[keep]++; ver[rem]++; merge_members(rem,keep); compact(rem); compact(keep); for(int fid:inc[keep]) if(F[fid].active&&hasv(F[fid],keep)){push_edge(F[fid].v[0],F[fid].v[1]);push_edge(F[fid].v[1],F[fid].v[2]);push_edge(F[fid].v[2],F[fid].v[0]);} }
static bool attempt(int a,int b){ if(a==b||!activeV[a]||!activeV[b])return false; int sh[2],opp[2]; if(!shared_faces(a,b,sh,opp))return false; if(!link_ok(a,b,opp))return false; Eval ab=eval_dir(a,b,sh), ba=eval_dir(b,a,sh); if(!ab.ok&&!ba.ok)return false; if(ba.ok&&(!ab.ok||ba.cost<ab.cost)) apply(b,a,sh,ba.pos); else apply(a,b,sh,ab.pos); return true; }
static void init(const Fingerprint&fp){ P=G.V; Orig=G.V; F.resize(M); for(int i=0;i<M;i++){F[i].v[0]=G.F[i].v[0];F[i].v[1]=G.F[i].v[1];F[i].v[2]=G.F[i].v[2];F[i].active=1;} activeV.assign(N,1); ver.assign(N,0); mark.assign(N,0); head.resize(N); tail.resize(N); nxt.assign(N,-1); csz.assign(N,1); for(int i=0;i<N;i++)head[i]=tail[i]=i; active_vertices=N; active_faces=M; token=1; haus2=pow(0.05*diag_len*0.999999,2); area_eps2=max(1e-32,1e-24*pow(diag_len,4)); normal_cos_limit = (fp.normal_sharp45>0.10?cos(55.0*M_PI/180.0):(fp.normal_smooth15>0.985?cos(82.0*M_PI/180.0):cos(68.0*M_PI/180.0)));
    vector<int> deg(N,0); for(auto &f:F){deg[f.v[0]]++;deg[f.v[1]]++;deg[f.v[2]]++;} inc.assign(N,{}); for(int i=0;i<N;i++)inc[i].reserve(deg[i]+8); for(int i=0;i<M;i++){inc[F[i].v[0]].push_back(i);inc[F[i].v[1]].push_back(i);inc[F[i].v[2]].push_back(i);} Q.assign(N,Quadric()); vector<pair<uint64_t,int>> edges; edges.reserve((size_t)M*3); vector<Vec3> fn(M); for(int i=0;i<M;i++){auto&f=F[i]; Vec3 n=fcross(f.v[0],f.v[1],f.v[2]); normalize(n); fn[i]=n; double d=-dot3(n,P[f.v[0]]); for(int k=0;k<3;k++) Q[f.v[k]].add_plane(n.x,n.y,n.z,d,1.0); edges.push_back({edge_key(f.v[0],f.v[1]),i});edges.push_back({edge_key(f.v[1],f.v[2]),i});edges.push_back({edge_key(f.v[2],f.v[0]),i});} sort(edges.begin(),edges.end()); double feature_cos=cos(32.0*M_PI/180.0); for(size_t i=0;i<edges.size();){size_t j=i+1; while(j<edges.size()&&edges[j].first==edges[i].first)j++; int a=(int)(edges[i].first>>32), b=(int)(edges[i].first&0xffffffffu); if(j-i==2){ double nd=dot3(fn[edges[i].second],fn[edges[i+1].second]); if(nd<feature_cos){ Vec3 e=P[b]-P[a]; if(normalize(e)){ Vec3 p1=cross3(e,fn[edges[i].second]); if(normalize(p1)){double d=-dot3(p1,P[a]);Q[a].add_plane(p1.x,p1.y,p1.z,d,8);Q[b].add_plane(p1.x,p1.y,p1.z,d,8);} Vec3 p2=cross3(fn[edges[i+1].second],e); if(normalize(p2)){double d=-dot3(p2,P[a]);Q[a].add_plane(p2.x,p2.y,p2.z,d,8);Q[b].add_plane(p2.x,p2.y,p2.z,d,8);} } } } push_edge(a,b); i=j; }
}
static int choose_target(const Fingerprint&fp){ double ratio=0.105; if(fp.normal_smooth15>0.985&&fp.normal_sharp45<0.01) ratio=(N>100000?0.055:0.070); else if(fp.normal_smooth30>0.94&&fp.val6_ratio>0.50) ratio=(N>100000?0.065:0.085); else if(fp.normal_sharp45>0.12||fp.val_high_ratio>0.12) ratio=0.135; if(N<1500) ratio=max(ratio,0.16); if(N<400) ratio=max(ratio,0.30); return max(4,(int)ceil(N*ratio)); }
static Candidate run(const Fingerprint&fp){ Candidate c; c.route="QEM_FALLBACK"; init(fp); int target=choose_target(fp); ELOG("[QEM] target=%d ratio=%.5f normalCos=%.5f pq=%zu\n",target,(double)target/N,normal_cos_limit,pq.size()); long long pops=0,coll=0; while(active_vertices>target&&!pq.empty()){ if((++pops&4095)==0 && elapsed()>18.5) break; Node cur=pq.top();pq.pop(); int a=cur.a,b=cur.b; if(a==b||!activeV[a]||!activeV[b])continue; if(cur.va!=ver[a]||cur.vb!=ver[b]){push_edge(a,b);continue;} if(attempt(a,b))coll++; } vector<int> id(N,-1); c.V.reserve(active_vertices); for(int i=0;i<N;i++) if(activeV[i]){id[i]=(int)c.V.size(); c.V.push_back(P[i]);} for(int i=0;i<M;i++){ if(!F[i].active)continue; int a=F[i].v[0],b=F[i].v[1],cc=F[i].v[2]; if(a==b||a==cc||b==cc)continue; if(id[a]<0||id[b]<0||id[cc]<0)continue; c.F.push_back({{id[a],id[b],id[cc]}}); } ELOG("[QEM] done activeV=%d outV=%zu outF=%zu pops=%lld collapses=%lld elapsed=%.3f\n",active_vertices,c.V.size(),c.F.size(),pops,coll,elapsed()); return c; }
}

static Candidate original_candidate(){ Candidate c; c.route="ORIGINAL_FALLBACK"; c.V=G.V; c.F=G.F; return c; }

static void print_mesh(const Candidate&c){
    printf("%d %d\n", (int)c.V.size(), (int)c.F.size());
    for(const Vec3&p:c.V) printf("v %.17g %.17g %.17g\n", p.x,p.y,p.z);
    for(const Face&f:c.F) printf("f %d %d %d\n", f.v[0]+1,f.v[1]+1,f.v[2]+1);
}

int main(){
    T0=chrono::steady_clock::now();
    if(!read_mesh()) return 0;
    Fingerprint fp=make_fingerprint();
    vector<function<Candidate()>> routes;
    routes.push_back([&](){return try_latgrid_route(fp);});
    routes.push_back([&](){return try_torusgrid_route(fp);});
    routes.push_back([&](){return try_sphere_ellipsoid_route(fp);});
    routes.push_back([&](){return try_torus_primitive_route(fp);});
    routes.push_back([&](){return try_revolve_route(fp);});
    // Box route is deliberately late and uses symmetric vertex guard; it only fires when vertex-Hausdorff is safe.
    routes.push_back([&](){ Candidate c=make_box_candidate(); return accept_candidate(c,fp,true,0.94)?c:Candidate{}; });
    for(auto &r:routes){
        if(elapsed()>16.2) break;
        Candidate c=r();
        if(!c.V.empty()){ print_mesh(c); return 0; }
    }
    Candidate q=QEM::run(fp);
    if((int)q.V.size()<N && !q.F.empty() && closed_manifold_ok(q.V,q.F) && candidate_vertex_guard(q.V,q.route,true)) { print_mesh(q); return 0; }
    ELOG("[FALLBACK] original mesh emitted\n");
    print_mesh(original_candidate());
    return 0;
}
