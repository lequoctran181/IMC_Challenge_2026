#include <bits/stdc++.h>
using namespace std;

struct Vec3{ double x,y,z; };
static inline Vec3 operator+(const Vec3&a,const Vec3&b){return {a.x+b.x,a.y+b.y,a.z+b.z};}
static inline Vec3 operator-(const Vec3&a,const Vec3&b){return {a.x-b.x,a.y-b.y,a.z-b.z};}
static inline Vec3 operator*(const Vec3&a,double s){return {a.x*s,a.y*s,a.z*s};}
static inline Vec3 operator/(const Vec3&a,double s){return {a.x/s,a.y/s,a.z/s};}
static inline double dot3(const Vec3&a,const Vec3&b){return a.x*b.x+a.y*b.y+a.z*b.z;}
static inline Vec3 cross3(const Vec3&a,const Vec3&b){return {a.y*b.z-a.z*b.y, a.z*b.x-a.x*b.z, a.x*b.y-a.y*b.x};}
static inline double norm2(const Vec3&a){return dot3(a,a);} 
static inline double norm3(const Vec3&a){return sqrt(max(0.0,norm2(a)));}
static inline Vec3 unit(Vec3 a){ double n=norm3(a); if(n<=1e-300) return {0,0,0}; return a/n; }
static inline double clampd(double x,double a,double b){return x<a?a:(x>b?b:x);} 

struct Face{ int v[3]; };

struct FastInput{
    vector<char> buf; char* p=nullptr;
    FastInput(){
        buf.reserve(1<<27);
        char chunk[1<<16]; size_t n;
        while((n=fread(chunk,1,sizeof(chunk),stdin))>0) buf.insert(buf.end(),chunk,chunk+n);
        buf.push_back('\0'); p=buf.data();
    }
    inline void skip(){ while(*p==' '||*p=='\n'||*p=='\r'||*p=='\t') ++p; }
    long long next_ll(){ skip(); long long s=1; if(*p=='-'){s=-1;++p;} long long x=0; while(*p>='0'&&*p<='9'){x=x*10+(*p-'0');++p;} return x*s; }
    int next_int(){ return (int)next_ll(); }
    double next_double(){ skip(); char* e=nullptr; double x=strtod(p,&e); p=e; return x; }
    char next_char(){ skip(); return *p++; }
};

static int N=0,M=0;
static vector<Vec3> V0;
static vector<Face> F0;
static Vec3 bbMin, bbMax;
static double diagLen=1.0, haus=0.05;
static chrono::steady_clock::time_point T0;
static inline double elapsed(){ return chrono::duration<double>(chrono::steady_clock::now()-T0).count(); }
static inline uint64_t edge_key(int a,int b){ if(a>b) swap(a,b); return ((uint64_t)(uint32_t)a<<32) | (uint32_t)b; }
static inline array<int,3> tri_key(int a,int b,int c){ array<int,3> t{a,b,c}; sort(t.begin(),t.end()); return t; }

static bool same_unordered(const Face&f,int a,int b,int c){ return tri_key(f.v[0],f.v[1],f.v[2]) == tri_key(a,b,c); }

static void read_input(){
    FastInput in;
    N=in.next_int(); M=in.next_int();
    V0.resize(N); F0.resize(M);
    bbMin={1e100,1e100,1e100}; bbMax={-1e100,-1e100,-1e100};
    for(int i=0;i<N;i++){
        char c=in.next_char(); (void)c;
        V0[i].x=in.next_double(); V0[i].y=in.next_double(); V0[i].z=in.next_double();
        bbMin.x=min(bbMin.x,V0[i].x); bbMin.y=min(bbMin.y,V0[i].y); bbMin.z=min(bbMin.z,V0[i].z);
        bbMax.x=max(bbMax.x,V0[i].x); bbMax.y=max(bbMax.y,V0[i].y); bbMax.z=max(bbMax.z,V0[i].z);
    }
    for(int i=0;i<M;i++){
        char c=in.next_char(); (void)c;
        F0[i].v[0]=in.next_int()-1; F0[i].v[1]=in.next_int()-1; F0[i].v[2]=in.next_int()-1;
    }
    diagLen=norm3(bbMax-bbMin); if(!(diagLen>0)) diagLen=1.0; haus=0.05*diagLen;
}

struct SpatialHash{
    double h=0.1; const vector<Vec3>* pts=nullptr; unordered_map<long long, vector<int>> mp;
    static long long mix(int x,int y,int z){
        uint64_t a=(uint32_t)x, b=(uint32_t)y, c=(uint32_t)z;
        uint64_t v=a*11995408973635179863ull ^ b*10150724397891781847ull ^ c*11056076308931462609ull;
        return (long long)v;
    }
    inline int qi(double x) const { return (int)floor(x/h); }
    void build(const vector<Vec3>&P,double hh){
        pts=&P; h=max(hh,1e-9); mp.clear(); mp.reserve(P.size()*2+1);
        for(int i=0;i<(int)P.size();++i){
            int ix=qi(P[i].x), iy=qi(P[i].y), iz=qi(P[i].z);
            mp[mix(ix,iy,iz)].push_back(i);
        }
    }
    bool any_within(const Vec3&p,double r2) const{
        int ix=qi(p.x), iy=qi(p.y), iz=qi(p.z);
        for(int dx=-1;dx<=1;dx++) for(int dy=-1;dy<=1;dy++) for(int dz=-1;dz<=1;dz++){
            auto it=mp.find(mix(ix+dx,iy+dy,iz+dz)); if(it==mp.end()) continue;
            const vector<int>& ids=it->second;
            for(int id:ids){ if(norm2((*pts)[id]-p) <= r2) return true; }
        }
        return false;
    }
};

static bool coverage_ok(const vector<Vec3>&A,const vector<Vec3>&B,double r){
    if(A.empty()||B.empty()) return false;
    const double r2=(r+1e-10)*(r+1e-10);
    SpatialHash hb; hb.build(B, max(r,1e-8));
    for(const Vec3&p:A) if(!hb.any_within(p,r2)) return false;
    SpatialHash ha; ha.build(A, max(r,1e-8));
    for(const Vec3&p:B) if(!ha.any_within(p,r2)) return false;
    return true;
}

static bool validate_basic(const vector<Vec3>&V,const vector<Face>&F){
    if(V.empty() || F.empty() || (int)V.size()>N) return false;
    double minA2=max(1e-32,1e-28*diagLen*diagLen*diagLen*diagLen);
    vector<uint64_t> e; e.reserve(F.size()*3);
    for(const Face&f:F){
        int a=f.v[0],b=f.v[1],c=f.v[2];
        if(a<0||b<0||c<0||a>=(int)V.size()||b>=(int)V.size()||c>=(int)V.size()) return false;
        if(a==b||a==c||b==c) return false;
        if(norm2(cross3(V[b]-V[a],V[c]-V[a]))<=minA2) return false;
        e.push_back(edge_key(a,b)); e.push_back(edge_key(b,c)); e.push_back(edge_key(c,a));
    }
    sort(e.begin(),e.end());
    for(size_t i=0;i<e.size();){ size_t j=i+1; while(j<e.size()&&e[j]==e[i]) ++j; if(j-i!=2) return false; i=j; }
    return true;
}

static double volume6(const vector<Vec3>&V,const vector<Face>&F){
    long double s=0;
    for(const Face&f:F){ s += (long double)dot3(V[f.v[0]], cross3(V[f.v[1]], V[f.v[2]])); }
    return (double)s;
}
static void orient_like_input(vector<Vec3>&V, vector<Face>&F){
    double a=volume6(V0,F0), b=volume6(V,F);
    if(a*b<0) for(Face&f:F) swap(f.v[1],f.v[2]);
}

static void write_output(const vector<Vec3>&V,const vector<Face>&F){
    string out; out.reserve(V.size()*42 + F.size()*24 + 64);
    char buf[128];
    out.append(buf, snprintf(buf,sizeof(buf), "%d %d\n", (int)V.size(), (int)F.size()));
    for(const Vec3&p:V) out.append(buf, snprintf(buf,sizeof(buf), "v %.10g %.10g %.10g\n", p.x,p.y,p.z));
    for(const Face&f:F) out.append(buf, snprintf(buf,sizeof(buf), "f %d %d %d\n", f.v[0]+1,f.v[1]+1,f.v[2]+1));
    fwrite(out.data(),1,out.size(),stdout);
}

// -----------------------------------------------------------------------------
// Exact/near-exact cuboid surface remesher.  It is deliberately strict: it only
// runs when every original vertex lies on one of the six AABB planes.  For those
// cases the generated mesh has identical silhouette/depth/flat normals, while a
// vertex grid spacing under the vertex-Hausdorff limit covers all input vertices.
// -----------------------------------------------------------------------------
static bool try_box(vector<Vec3>&outV, vector<Face>&outF){
    if(N<1000) return false;
    Vec3 ext=bbMax-bbMin;
    if(ext.x<=1e-12 || ext.y<=1e-12 || ext.z<=1e-12) return false;
    double tol=max(1e-8, diagLen*2e-6);
    int on[6]={0,0,0,0,0,0};
    for(const Vec3&p:V0){
        double d[6]={fabs(p.x-bbMin.x),fabs(p.x-bbMax.x),fabs(p.y-bbMin.y),fabs(p.y-bbMax.y),fabs(p.z-bbMin.z),fabs(p.z-bbMax.z)};
        int k=min_element(d,d+6)-d;
        if(d[k]>tol) return false;
        on[k]++;
    }
    for(int k=0;k<6;k++) if(on[k]==0) return false;
    double step=max(haus*1.30, 1e-9);
    int nx=max(1,(int)ceil(ext.x/step));
    int ny=max(1,(int)ceil(ext.y/step));
    int nz=max(1,(int)ceil(ext.z/step));
    // Slightly densify if a face diagonal cell would exceed Hausdorff.
    auto fix2=[&](int &a,int &b,double A,double B){
        for(int it=0; it<6; ++it){
            double da=A/max(1,a), db=B/max(1,b);
            if(0.5*sqrt(da*da+db*db) <= haus*0.985) break;
            if(da>db) ++a; else ++b;
        }
    };
    fix2(ny,nz,ext.y,ext.z); fix2(nx,nz,ext.x,ext.z); fix2(nx,ny,ext.x,ext.y);
    auto vid = [&](int ix,int iy,int iz, unordered_map<long long,int>&id, vector<Vec3>&V)->int{
        long long key=((long long)ix<<42) ^ ((long long)iy<<21) ^ (long long)iz;
        auto it=id.find(key); if(it!=id.end()) return it->second;
        Vec3 p{bbMin.x + ext.x*ix/nx, bbMin.y + ext.y*iy/ny, bbMin.z + ext.z*iz/nz};
        int v=(int)V.size(); id[key]=v; V.push_back(p); return v;
    };
    vector<Vec3> V; vector<Face> F; unordered_map<long long,int> id; id.reserve((nx+1)*(ny+1)*2+(nx+1)*(nz+1)*2+(ny+1)*(nz+1)*2);
    auto addq=[&](int a,int b,int c,int d, Vec3 outward){
        Face f1{{a,b,c}}, f2{{a,c,d}};
        Vec3 cr1=cross3(V[f1.v[1]]-V[f1.v[0]], V[f1.v[2]]-V[f1.v[0]]);
        if(dot3(cr1,outward)<0) swap(f1.v[1],f1.v[2]);
        Vec3 cr2=cross3(V[f2.v[1]]-V[f2.v[0]], V[f2.v[2]]-V[f2.v[0]]);
        if(dot3(cr2,outward)<0) swap(f2.v[1],f2.v[2]);
        F.push_back(f1); F.push_back(f2);
    };
    for(int side=0; side<6; ++side){
        if(side<2){ int ix = side==0?0:nx; Vec3 n{side==0?-1.0:1.0,0,0};
            for(int iy=0;iy<ny;iy++) for(int iz=0;iz<nz;iz++){
                int a=vid(ix,iy,iz,id,V), b=vid(ix,iy+1,iz,id,V), c=vid(ix,iy+1,iz+1,id,V), d=vid(ix,iy,iz+1,id,V); addq(a,b,c,d,n);
            }
        }else if(side<4){ int iy = side==2?0:ny; Vec3 n{0,side==2?-1.0:1.0,0};
            for(int ix=0;ix<nx;ix++) for(int iz=0;iz<nz;iz++){
                int a=vid(ix,iy,iz,id,V), b=vid(ix,iy,iz+1,id,V), c=vid(ix+1,iy,iz+1,id,V), d=vid(ix+1,iy,iz,id,V); addq(a,b,c,d,n);
            }
        }else{ int iz = side==4?0:nz; Vec3 n{0,0,side==4?-1.0:1.0};
            for(int ix=0;ix<nx;ix++) for(int iy=0;iy<ny;iy++){
                int a=vid(ix,iy,iz,id,V), b=vid(ix+1,iy,iz,id,V), c=vid(ix+1,iy+1,iz,id,V), d=vid(ix,iy+1,iz,id,V); addq(a,b,c,d,n);
            }
        }
    }
    if((int)V.size()>=N) return false;
    if(!validate_basic(V,F)) return false;
    if(!coverage_ok(V0,V,haus*0.999)) return false;
    outV.swap(V); outF.swap(F); return true;
}


static bool ellipsoid_quality(double &ax,double &ay,double &az){
    ax=max(fabs(bbMin.x),fabs(bbMax.x)); ay=max(fabs(bbMin.y),fabs(bbMax.y)); az=max(fabs(bbMin.z),fabs(bbMax.z));
    if(ax<diagLen*0.08||ay<diagLen*0.08||az<diagLen*0.08) return false;
    long double se=0; double mx=0;
    int stride=max(1,N/200000), cnt=0;
    for(int i=0;i<N;i+=stride){
        const Vec3&p=V0[i]; double q=(p.x*p.x)/(ax*ax)+(p.y*p.y)/(ay*ay)+(p.z*p.z)/(az*az); double d=fabs(q-1.0); se += d*d; mx=max(mx,d); cnt++;
    }
    double rms=sqrt((double)(se/max(1,cnt)));
    return rms<0.012 && mx<0.055;
}
static bool try_ellipsoid_like(vector<Vec3>&outV, vector<Face>&outF){
    if(N<3000) return false;
    double ax,ay,az; if(!ellipsoid_quality(ax,ay,az)) return false;
    double ratio=(N>=200000?0.055:(N>=40000?0.070:0.100));
    int target=max(96,(int)(N*ratio));
    for(int attempt=0; attempt<8; ++attempt){
        int body=max(1,target-2);
        int R=max(4,(int)sqrt(body/2.0));
        int S=max(12,body/max(1,R));
        if(S>4*R) { S=4*R; R=max(4,body/max(1,S)); }
        vector<Vec3> V; vector<Face> F; V.reserve(2+R*S);
        V.push_back({0,0,az});
        for(int i=0;i<R;i++){
            double th=M_PI*(i+1)/(R+1); double st=sin(th), ct=cos(th);
            for(int j=0;j<S;j++){ double ph=2*M_PI*j/S; V.push_back({ax*st*cos(ph), ay*st*sin(ph), az*ct}); }
        }
        int bot=(int)V.size(); V.push_back({0,0,-az});
        auto id=[&](int i,int j){return 1+i*S+((j%S+S)%S);};
        auto add=[&](int a,int b,int c){ Face f{{a,b,c}}; Vec3 cr=cross3(V[b]-V[a],V[c]-V[a]); Vec3 cc=(V[a]+V[b]+V[c])/3.0; if(dot3(cr,cc)<0) swap(f.v[1],f.v[2]); F.push_back(f); };
        for(int j=0;j<S;j++) add(0,id(0,j+1),id(0,j));
        for(int i=0;i<R-1;i++) for(int j=0;j<S;j++){ int a=id(i,j),b=id(i,j+1),c=id(i+1,j),d=id(i+1,j+1); add(a,b,c); add(b,d,c); }
        for(int j=0;j<S;j++) add(bot,id(R-1,j),id(R-1,j+1));
        if((int)V.size()<N && validate_basic(V,F) && coverage_ok(V0,V,haus*0.995)){ outV.swap(V); outF.swap(F); return true; }
        target=min(N-1,(int)(target*1.5)+64); if(target>=N) break;
    }
    return false;
}

// -----------------------------------------------------------------------------
// Ordered latitude-longitude sphere recognizer/resampler.
// -----------------------------------------------------------------------------
static bool detect_latlong(int &R,int &S){
    if(N<50 || M != 2*(N-2)) return false;
    int X=N-2;
    vector<int> cand;
    for(int d=3;(long long)d*d<=X;d++) if(X%d==0){ cand.push_back(d); if(d*d!=X) cand.push_back(X/d); }
    sort(cand.begin(),cand.end());
    for(int s:cand){
        int r=X/s; if(s<6||r<3) continue;
        bool ok=true; int step=max(1,s/64);
        for(int j=0;j<s && ok;j+=step){ int a=1+j,b=1+(j+1)%s; if(!same_unordered(F0[j],0,b,a)) ok=false; }
        int span=max(1,(r-1)*s/256);
        for(int q=0;q<(r-1)*s && ok;q+=span){
            int rr=q/s, j=q-rr*s; int a=1+rr*s+j, b=1+rr*s+(j+1)%s, c=1+(rr+1)*s+j, d=1+(rr+1)*s+(j+1)%s; int fid=s+2*q;
            if(fid+1>=M) {ok=false; break;}
            if(!same_unordered(F0[fid],a,b,c)) ok=false;
            if(ok && !same_unordered(F0[fid+1],b,d,c)) ok=false;
        }
        int bottom=s+2*(r-1)*s;
        for(int j=0;j<s && ok;j+=step){ int c=1+(r-1)*s+j,d=1+(r-1)*s+(j+1)%s; if(bottom+j>=M||!same_unordered(F0[bottom+j],N-1,c,d)) ok=false; }
        if(ok){R=r;S=s;return true;}
    }
    return false;
}
static bool sphere_quality_ok(){
    double mean=0, mxdev=0; for(const Vec3&p:V0) mean += norm3(p); mean/=max(1,N);
    if(mean<0.15) return false;
    for(const Vec3&p:V0) mxdev=max(mxdev,fabs(norm3(p)-mean));
    return mxdev <= max(haus*0.35, mean*0.035);
}
static bool try_latlong_sphere(vector<Vec3>&outV, vector<Face>&outF){
    int R,S; if(!detect_latlong(R,S)) return false; if(!sphere_quality_ok()) return false;
    double ratio = (N>=200000?0.070:(N>=40000?0.085:0.110));
    int target=max(32,(int)(N*ratio));
    for(int attempt=0; attempt<8; ++attempt){
        int body=max(1,target-2);
        double ar=(double)R/max(1,S);
        int R2=max(2,min(R,(int)(sqrt(body*ar)+0.5)));
        int S2=max(6,min(S,body/max(1,R2)));
        // keep enough sides for normal-map SSIM; avoid too-coarse facets.
        if(S2<16 && S>=16) S2=min(S,16);
        if(R2<8 && R>=8) R2=min(R,8);
        vector<Vec3> V; vector<Face> F; V.reserve(2+R2*S2);
        vector<int> src; src.reserve(2+R2*S2);
        V.push_back(V0[0]); src.push_back(0);
        for(int i=0;i<R2;i++){
            int oi = (R2<=1?0:(int)llround((double)i*(R-1)/(R2-1))); if(oi<0) oi=0; if(oi>=R) oi=R-1;
            for(int j=0;j<S2;j++){
                int oj=(long long)j*S/S2; if(oj>=S) oj=S-1;
                int id=1+oi*S+oj; V.push_back(V0[id]); src.push_back(id);
            }
        }
        int bot=(int)V.size(); V.push_back(V0[N-1]); src.push_back(N-1);
        auto id=[&](int i,int j){return 1+i*S2+((j%S2+S2)%S2);};
        auto add=[&](int a,int b,int c){ Face f{{a,b,c}}; Vec3 cr=cross3(V[b]-V[a],V[c]-V[a]); Vec3 cc=(V[a]+V[b]+V[c])/3.0; if(dot3(cr,cc)<0) swap(f.v[1],f.v[2]); F.push_back(f); };
        for(int j=0;j<S2;j++) add(0,id(0,j+1),id(0,j));
        for(int i=0;i<R2-1;i++) for(int j=0;j<S2;j++){ int a=id(i,j),b=id(i,j+1),c=id(i+1,j),d=id(i+1,j+1); add(a,b,c); add(b,d,c); }
        for(int j=0;j<S2;j++) add(bot,id(R2-1,j),id(R2-1,j+1));
        if((int)V.size()<N && validate_basic(V,F) && coverage_ok(V0,V,haus*0.999)){ outV.swap(V); outF.swap(F); return true; }
        target = min(N-1, (int)(target*1.45)+16);
        if(target>=N) break;
    }
    return false;
}


static vector<Vec3> vertex_normals();
static bool latlong_revolve_ok(int R,int S){
    if(R<2||S<8) return false;
    vector<Vec3> cen(R,{0,0,0});
    for(int i=0;i<R;i++){
        for(int j=0;j<S;j++) cen[i]=cen[i]+V0[1+i*S+j];
        cen[i]=cen[i]/(double)S;
    }
    Vec3 axis=cen.back()-cen.front(); double al=norm3(axis);
    if(al < diagLen*0.08) return false;
    axis=axis/al;
    int good=0, used=0;
    for(int i=0;i<R;i++){
        double mean=0, ss=0; vector<double> rr; rr.reserve(S);
        for(int j=0;j<S;j++){
            Vec3 d=V0[1+i*S+j]-cen[i];
            Vec3 perp=d-axis*dot3(d,axis);
            double r=norm3(perp); rr.push_back(r); mean+=r;
        }
        mean/=S; if(mean<diagLen*1e-5) continue; used++;
        for(double r:rr) ss+=(r-mean)*(r-mean);
        double cv=sqrt(ss/S)/max(mean,1e-12);
        if(cv<0.035) good++;
    }
    if(used<max(2,R/3)) return false;
    return good*10 >= used*9;
}
static bool try_latlong_revolve(vector<Vec3>&outV, vector<Face>&outF){
    int R,S; if(!detect_latlong(R,S)) return false;
    if(!latlong_revolve_ok(R,S)) return false;
    vector<Vec3> vn=vertex_normals();
    double ratio = (N>=200000?0.080:(N>=40000?0.095:0.125));
    int target=max(48,(int)(N*ratio));
    for(int attempt=0; attempt<8; ++attempt){
        int body=max(1,target-2);
        double ar=(double)R/max(1,S);
        int R2=max(2,min(R,(int)(sqrt(body*ar)+0.5)));
        int S2=max(8,min(S,body/max(1,R2)));
        if(S2<20 && S>=20) S2=min(S,20);
        vector<Vec3> V; vector<Face> F; vector<int> src; V.reserve(2+R2*S2); src.reserve(2+R2*S2);
        V.push_back(V0[0]); src.push_back(0);
        for(int i=0;i<R2;i++){
            int oi=(R2<=1?0:(int)llround((double)i*(R-1)/(R2-1))); if(oi<0) oi=0; if(oi>=R) oi=R-1;
            for(int j=0;j<S2;j++){
                int oj=(long long)j*S/S2; if(oj>=S) oj=S-1;
                int sid=1+oi*S+oj; V.push_back(V0[sid]); src.push_back(sid);
            }
        }
        int bot=(int)V.size(); V.push_back(V0[N-1]); src.push_back(N-1);
        auto id=[&](int i,int j){return 1+i*S2+((j%S2+S2)%S2);};
        auto add=[&](int a,int b,int c){ Face f{{a,b,c}}; Vec3 cr=cross3(V[b]-V[a],V[c]-V[a]); Vec3 ref=vn[src[a]]+vn[src[b]]+vn[src[c]]; if(dot3(cr,ref)<0) swap(f.v[1],f.v[2]); F.push_back(f); };
        for(int j=0;j<S2;j++) add(0,id(0,j+1),id(0,j));
        for(int i=0;i<R2-1;i++) for(int j=0;j<S2;j++){ int a=id(i,j),b=id(i,j+1),c=id(i+1,j),d=id(i+1,j+1); add(a,b,c); add(b,d,c); }
        for(int j=0;j<S2;j++) add(bot,id(R2-1,j),id(R2-1,j+1));
        if((int)V.size()<N && validate_basic(V,F) && coverage_ok(V0,V,haus*0.999)){ outV.swap(V); outF.swap(F); return true; }
        target=min(N-1,(int)(target*1.45)+32); if(target>=N) break;
    }
    return false;
}

// -----------------------------------------------------------------------------
// Ordered periodic UxS grid / torus-tube resampler.  Very strict face-order
// detector; selected vertices are original vertices, so only original->output
// coverage is the limiting Hausdorff direction.
// -----------------------------------------------------------------------------
static bool adj2_all(const int a[3],int mod,int &base){
    for(int t=0;t<3;t++) for(int s=0;s<2;s++){
        int x=(a[t]-s+mod)%mod; bool ok=true;
        for(int i=0;i<3;i++){ int d=(a[i]-x+mod)%mod; if(d!=0&&d!=1){ok=false;break;} }
        if(ok){base=x; return true;}
    }
    return false;
}
static bool torus_face_ok(const Face&f,int S){
    if(S<3||N%S) return false; int U=N/S; if(U<3) return false;
    int r[3]={f.v[0]/S,f.v[1]/S,f.v[2]/S};
    int c[3]={f.v[0]%S,f.v[1]%S,f.v[2]%S};
    int rb=0,cb=0; if(!adj2_all(r,U,rb)||!adj2_all(c,S,cb)) return false;
    int mask=0;
    for(int i=0;i<3;i++){ int x=(r[i]-rb+U)%U, y=(c[i]-cb+S)%S; if(x>1||y>1) return false; mask |= 1<<(x*2+y); }
    return __builtin_popcount((unsigned)mask)==3;
}
static bool detect_torus(int &U,int &S){
    if(N<100 || M!=2*N) return false;
    map<int,int> cnt; int stride=max(1,M/120000);
    for(int i=0;i<M;i+=stride){
        int a[3]={F0[i].v[0],F0[i].v[1],F0[i].v[2]};
        for(int k=0;k<3;k++){ int d=abs(a[k]-a[(k+1)%3]); d=min(d,N-d); if(d>=3&&d<=N/3) cnt[d]++; }
    }
    vector<pair<int,int>> q; for(auto&p:cnt) q.push_back({p.second,p.first}); sort(q.rbegin(),q.rend());
    vector<int> cand;
    auto add=[&](int s){ if(s>=6&&s<=N/3&&N%s==0&&find(cand.begin(),cand.end(),s)==cand.end()) cand.push_back(s); };
    for(int i=0;i<(int)q.size()&&i<16;i++){ int d=q[i].second; for(int e=-2;e<=2;e++) add(d+e); if(d) add(N/d); }
    for(int s:cand){
        int u=N/s; if(u<6) continue; int st=max(1,M/200000), total=0, ok=0;
        for(int i=0;i<M;i+=st){ total++; if(torus_face_ok(F0[i],s)) ok++; }
        if(total>0 && ok*1000>=total*995){ U=u; S=s; return true; }
    }
    return false;
}
static vector<Vec3> vertex_normals(){
    vector<Vec3> n(N,{0,0,0});
    for(const Face&f:F0){ Vec3 cr=cross3(V0[f.v[1]]-V0[f.v[0]], V0[f.v[2]]-V0[f.v[0]]); for(int k=0;k<3;k++) n[f.v[k]]=n[f.v[k]]+cr; }
    return n;
}
static bool try_torus_grid(vector<Vec3>&outV, vector<Face>&outF){
    int U,S; if(!detect_torus(U,S)) return false;
    vector<Vec3> vn=vertex_normals();
    double ratio = (N>=200000?0.075:(N>=40000?0.090:0.120));
    int target=max(64,(int)(N*ratio));
    for(int attempt=0; attempt<8; ++attempt){
        double ar=sqrt((double)U/max(1,S));
        int U2=max(6,min(U,(int)(sqrt((double)target)*ar+0.5)));
        int S2=max(6,min(S,target/max(1,U2)));
        if(U2*S2>=N) break;
        vector<Vec3> V; vector<Face> F; vector<int> src; V.reserve(U2*S2); src.reserve(U2*S2);
        for(int i=0;i<U2;i++){ int oi=(long long)i*U/U2; if(oi>=U) oi=U-1; for(int j=0;j<S2;j++){ int oj=(long long)j*S/S2; if(oj>=S) oj=S-1; int id=oi*S+oj; src.push_back(id); V.push_back(V0[id]); }}
        auto id=[&](int i,int j){return ((i%U2+U2)%U2)*S2 + ((j%S2+S2)%S2);};
        auto add=[&](int a,int b,int c){ Face f{{a,b,c}}; Vec3 cr=cross3(V[b]-V[a],V[c]-V[a]); Vec3 ref=vn[src[a]]+vn[src[b]]+vn[src[c]]; if(dot3(cr,ref)<0) swap(f.v[1],f.v[2]); F.push_back(f); };
        for(int i=0;i<U2;i++) for(int j=0;j<S2;j++){ int a=id(i,j), b=id(i,j+1), c=id(i+1,j), d=id(i+1,j+1); add(a,b,c); add(b,d,c); }
        if(validate_basic(V,F) && coverage_ok(V0,V,haus*0.999)){ outV.swap(V); outF.swap(F); return true; }
        target = min(N-1,(int)(target*1.45)+32);
        if(target>=N) break;
    }
    return false;
}

// -----------------------------------------------------------------------------
// Topology-preserving endpoint edge-collapse fallback.
// -----------------------------------------------------------------------------
struct Simplifier{
    vector<Vec3> V;
    vector<Face> F;
    vector<unsigned char> aliveV, aliveF;
    vector<vector<int>> inc;
    vector<double> cover;
    int vcnt=0, fcnt=0;
    vector<int> markA, markB; int tokenA=1, tokenB=1;
    double minA2=1e-32;
    struct Params{ double maxCover, planeTol, minDot, minDotHard, areaRatio; };
    struct Eval{ bool ok=false; double cost=1e100; double newCover=0; };
    struct PQE{ double key; int u,v; bool operator<(const PQE&o)const{return key>o.key;} };

    Simplifier(){
        V=V0; F=F0; aliveV.assign(N,1); aliveF.assign(M,1); cover.assign(N,0); vcnt=N; fcnt=M; minA2=max(1e-32,1e-30*diagLen*diagLen*diagLen*diagLen);
        vector<int> deg(N,0); for(const Face&f:F){ deg[f.v[0]]++; deg[f.v[1]]++; deg[f.v[2]]++; }
        inc.resize(N); for(int i=0;i<N;i++) inc[i].reserve(deg[i]);
        for(int i=0;i<M;i++){ inc[F[i].v[0]].push_back(i); inc[F[i].v[1]].push_back(i); inc[F[i].v[2]].push_back(i); }
        markA.assign(N,0); markB.assign(N,0);
    }
    inline bool contains(int fid,int v) const { const Face&f=F[fid]; return f.v[0]==v||f.v[1]==v||f.v[2]==v; }
    inline int other_in_face(int fid,int a,int b) const { const Face&f=F[fid]; for(int k=0;k<3;k++){int x=f.v[k]; if(x!=a&&x!=b) return x;} return -1; }
    bool edge_faces(int u,int v,int ef[2],int opp[2]){
        int cnt=0;
        const vector<int>& L = inc[u].size()<inc[v].size()?inc[u]:inc[v];
        for(int fid:L){ if(!aliveF[fid]) continue; if(!contains(fid,u)||!contains(fid,v)) continue; if(cnt>=2) return false; ef[cnt]=fid; opp[cnt]=other_in_face(fid,u,v); cnt++; }
        if(cnt!=2) return false; if(opp[0]<0||opp[1]<0||opp[0]==opp[1]) return false; return true;
    }
    bool link_ok(int u,int v,const int opp[2]){
        if(++tokenA>2000000000){ fill(markA.begin(),markA.end(),0); tokenA=1; }
        if(++tokenB>2000000000){ fill(markB.begin(),markB.end(),0); tokenB=1; }
        for(int fid:inc[u]) if(aliveF[fid]&&contains(fid,u)){
            const Face&f=F[fid]; for(int k=0;k<3;k++){ int x=f.v[k]; if(x!=u&&x!=v) markA[x]=tokenA; }
        }
        int common=0;
        for(int fid:inc[v]) if(aliveF[fid]&&contains(fid,v)){
            const Face&f=F[fid]; for(int k=0;k<3;k++){ int x=f.v[k]; if(x==u||x==v) continue; if(markA[x]!=tokenA) continue; if(x!=opp[0]&&x!=opp[1]) return false; if(markB[x]!=tokenB){ markB[x]=tokenB; common++; } }
        }
        return common==2 && markB[opp[0]]==tokenB && markB[opp[1]]==tokenB;
    }
    bool duplicate_existing(int a,int b,int c,int skip1,int skip2,int skip3) const{
        int pick=a; if(inc[b].size()<inc[pick].size()) pick=b; if(inc[c].size()<inc[pick].size()) pick=c;
        auto key=tri_key(a,b,c);
        for(int fid:inc[pick]){
            if(!aliveF[fid]||fid==skip1||fid==skip2||fid==skip3) continue;
            const Face&f=F[fid]; if(tri_key(f.v[0],f.v[1],f.v[2])==key) return true;
        }
        return false;
    }
    Eval evaluate(int keep,int rem,const int ef[2],const Params&p){
        Eval e; double len=norm3(V[keep]-V[rem]); e.newCover=max(cover[keep], cover[rem]+len); if(e.newCover>p.maxCover) return e;
        vector<array<int,3>> newkeys; newkeys.reserve(inc[rem].size());
        double maxPlane=0, maxFlip=0, maxAreaLoss=0; int changed=0;
        for(int fid:inc[rem]){
            if(!aliveF[fid]||!contains(fid,rem)) continue; if(fid==ef[0]||fid==ef[1]) continue;
            if(contains(fid,keep)) return e;
            Face old=F[fid]; int a=old.v[0],b=old.v[1],c=old.v[2];
            int na=a,nb=b,nc=c; if(na==rem) na=keep; if(nb==rem) nb=keep; if(nc==rem) nc=keep;
            if(na==nb||na==nc||nb==nc) return e;
            Vec3 cr0=cross3(V[b]-V[a],V[c]-V[a]); double ar0=norm3(cr0); if(!(ar0>0)) return e;
            Vec3 cr1=cross3(V[nb]-V[na],V[nc]-V[na]); double ar1=norm3(cr1); if(ar1*ar1<=minA2) return e;
            double nd=dot3(cr0,cr1)/(ar0*ar1); nd=clampd(nd,-1,1); if(nd<p.minDotHard) return e;
            Vec3 n=cr0/ar0; double plane=fabs(dot3(n,V[keep]-V[rem]));
            if(plane>p.planeTol && e.newCover>haus*0.25) return e;
            if(nd<p.minDot && plane>p.planeTol*0.35) return e;
            if(duplicate_existing(na,nb,nc,fid,ef[0],ef[1])) return e;
            newkeys.push_back(tri_key(na,nb,nc));
            maxPlane=max(maxPlane,plane); maxFlip=max(maxFlip,1.0-nd); maxAreaLoss=max(maxAreaLoss,max(0.0,1.0-ar1/ar0)); changed++;
        }
        if(changed==0) return e;
        sort(newkeys.begin(),newkeys.end()); for(int i=1;i<(int)newkeys.size();++i) if(newkeys[i]==newkeys[i-1]) return e;
        e.ok=true; e.cost=len/(haus+1e-12) + 0.9*e.newCover/(haus+1e-12) + 80.0*maxFlip + 5.0*maxPlane/(haus+1e-12) + 0.02*changed + 0.3*maxAreaLoss;
        return e;
    }
    void merge_inc(int keep,int rem){
        vector<int> m; m.reserve(inc[keep].size()+inc[rem].size());
        for(int fid:inc[keep]) if(aliveF[fid]&&contains(fid,keep)) m.push_back(fid);
        for(int fid:inc[rem]) if(aliveF[fid]&&contains(fid,keep)) m.push_back(fid);
        sort(m.begin(),m.end()); m.erase(unique(m.begin(),m.end()),m.end()); inc[keep].swap(m); vector<int>().swap(inc[rem]);
    }
    void apply(int keep,int rem,const int ef[2],double newCover){
        for(int k=0;k<2;k++) if(aliveF[ef[k]]){ aliveF[ef[k]]=0; fcnt--; }
        for(int fid:inc[rem]) if(aliveF[fid]&&contains(fid,rem)){
            Face&f=F[fid]; for(int i=0;i<3;i++) if(f.v[i]==rem) f.v[i]=keep;
        }
        aliveV[rem]=0; vcnt--; cover[keep]=newCover; merge_inc(keep,rem);
    }
    bool collapse_edge(int u,int v,const Params&p){
        if(u==v||u<0||v<0||u>=N||v>=N||!aliveV[u]||!aliveV[v]) return false;
        int ef[2],opp[2]; if(!edge_faces(u,v,ef,opp)) return false; if(!link_ok(u,v,opp)) return false;
        Eval a=evaluate(u,v,ef,p), b=evaluate(v,u,ef,p);
        if(!a.ok&&!b.ok) return false;
        if(b.ok && (!a.ok||b.cost<a.cost)) apply(v,u,ef,b.newCover); else apply(u,v,ef,a.newCover);
        return true;
    }
    void build_edges(priority_queue<PQE>&pq){
        vector<uint64_t> es; es.reserve((size_t)fcnt*3);
        for(int i=0;i<M;i++) if(aliveF[i]){ const Face&f=F[i]; es.push_back(edge_key(f.v[0],f.v[1])); es.push_back(edge_key(f.v[1],f.v[2])); es.push_back(edge_key(f.v[2],f.v[0])); }
        sort(es.begin(),es.end()); es.erase(unique(es.begin(),es.end()),es.end());
        for(uint64_t k:es){ int a=(int)(k>>32), b=(int)(uint32_t)k; if(aliveV[a]&&aliveV[b]){ double l=norm2(V[a]-V[b]); pq.push({l,a,b}); } }
    }
    void push_incident(priority_queue<PQE>&pq,int v){
        if(v<0||v>=N||!aliveV[v]) return;
        for(int fid:inc[v]) if(aliveF[fid]&&contains(fid,v)){
            const Face&f=F[fid]; for(int k=0;k<3;k++){ int x=f.v[k]; if(x!=v&&aliveV[x]) pq.push({norm2(V[v]-V[x]),min(v,x),max(v,x)}); }
        }
    }
    void run_pass(const Params&p,int target,double timeLimit){
        priority_queue<PQE> pq; build_edges(pq);
        int attempts=0, acc=0;
        while(!pq.empty() && vcnt>target && elapsed()<timeLimit){
            PQE e=pq.top(); pq.pop(); attempts++;
            int u=e.u,v=e.v; if(!aliveV[u]||!aliveV[v]) continue;
            int before=vcnt;
            if(collapse_edge(u,v,p)){
                acc++; int keep=-1; if(aliveV[u]&&vcnt==before-1) keep=u; if(aliveV[v]&&vcnt==before-1) keep=v; if(keep>=0) push_incident(pq,keep);
            }
            if((attempts&131071)==0 && elapsed()>timeLimit) break;
        }
    }
    void output(vector<Vec3>&OV, vector<Face>&OF){
        vector<int> mp(N,-1); int n=0;
        for(int i=0;i<M;i++) if(aliveF[i]){ Face&f=F[i]; for(int k=0;k<3;k++){ int x=f.v[k]; if(mp[x]<0) mp[x]=n++; } }
        OV.assign(n,{0,0,0}); for(int i=0;i<N;i++) if(mp[i]>=0) OV[mp[i]]=V[i];
        OF.clear(); OF.reserve(fcnt);
        for(int i=0;i<M;i++) if(aliveF[i]){ Face f=F[i]; f.v[0]=mp[f.v[0]]; f.v[1]=mp[f.v[1]]; f.v[2]=mp[f.v[2]]; if(f.v[0]!=f.v[1]&&f.v[0]!=f.v[2]&&f.v[1]!=f.v[2]) OF.push_back(f); }
    }
};

static bool try_general(vector<Vec3>&outV, vector<Face>&outF){
    Simplifier S;
    // The first pass is nearly planar and is what solves the sample exactly.
    vector<Simplifier::Params> P;
    P.push_back({haus*0.999, diagLen*1e-8, 0.9999995, 0.99999, 0.02});
    P.push_back({haus*0.995, diagLen*0.0015, 0.9990, 0.9960, 0.05});
    P.push_back({haus*0.992, diagLen*0.0050, 0.9950, 0.9850, 0.08});
    P.push_back({haus*0.988, diagLen*0.0120, 0.9850, 0.9550, 0.12});
    if(N>=60000) P.push_back({haus*0.985, diagLen*0.0250, 0.9650, 0.9000, 0.20});
    vector<double> ratios;
    ratios.push_back(0.55);
    ratios.push_back(N<20000?0.28:0.24);
    ratios.push_back(N<20000?0.18:0.15);
    ratios.push_back(N<60000?0.13:0.105);
    if(N>=60000) ratios.push_back(0.082);
    for(int i=0;i<(int)P.size();++i){
        int target=max(4,(int)(N*ratios[min(i,(int)ratios.size()-1)]));
        double tl = (i+1==(int)P.size()?19.2: (9.0 + 2.5*i));
        if(N<1000) tl=19.5;
        S.run_pass(P[i], target, tl);
        if(elapsed()>19.25) break;
    }
    S.output(outV,outF);
    if(!validate_basic(outV,outF)){
        // fail closed: original mesh is always valid and accepted, but sample still
        // should have simplified before this point.
        outV=V0; outF=F0; return false;
    }
    orient_like_input(outV,outF);
    return true;
}

int main(){
    T0=chrono::steady_clock::now();
    read_input();
    vector<Vec3> outV; vector<Face> outF;
    bool ok=false;
    // The structural routes are narrow and validate vertex-Hausdorff exactly.
    if(!ok) ok=try_box(outV,outF);
    if(!ok) ok=try_ellipsoid_like(outV,outF);
    if(!ok) ok=try_latlong_sphere(outV,outF);
    if(!ok) ok=try_latlong_revolve(outV,outF);
    if(!ok) ok=try_torus_grid(outV,outF);
    if(!ok) ok=try_general(outV,outF);
    if(!validate_basic(outV,outF) || (int)outV.size()>N){ outV=V0; outF=F0; }
    write_output(outV,outF);
    return 0;
}
