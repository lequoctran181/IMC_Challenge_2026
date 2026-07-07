#include <bits/stdc++.h>
using namespace std;

struct Vec3 {
    double x=0,y=0,z=0;
};
static inline Vec3 operator+(const Vec3&a,const Vec3&b){return {a.x+b.x,a.y+b.y,a.z+b.z};}
static inline Vec3 operator-(const Vec3&a,const Vec3&b){return {a.x-b.x,a.y-b.y,a.z-b.z};}
static inline Vec3 operator*(const Vec3&a,double s){return {a.x*s,a.y*s,a.z*s};}
static inline Vec3 operator/(const Vec3&a,double s){return {a.x/s,a.y/s,a.z/s};}
static inline double dot3(const Vec3&a,const Vec3&b){return a.x*b.x+a.y*b.y+a.z*b.z;}
static inline Vec3 cross3(const Vec3&a,const Vec3&b){return {a.y*b.z-a.z*b.y, a.z*b.x-a.x*b.z, a.x*b.y-a.y*b.x};}
static inline double norm2(const Vec3&a){return dot3(a,a);} 
static inline double norm3(const Vec3&a){return sqrt(norm2(a));}
static inline Vec3 normalize3(const Vec3&a){ double n=norm3(a); return n>0? a*(1.0/n):Vec3{0,0,0}; }

struct Face { int v[3]; };

struct FastInput {
    vector<char> buf; char *p=nullptr;
    FastInput(){
        buf.reserve(1<<27); char chunk[1<<16]; size_t n;
        while((n=fread(chunk,1,sizeof(chunk),stdin))>0) buf.insert(buf.end(),chunk,chunk+n);
        buf.push_back('\0'); p=buf.data();
    }
    inline void skip(){ while(*p==' '||*p=='\n'||*p=='\r'||*p=='\t') ++p; }
    long next_long(){ skip(); return strtol(p,&p,10); }
    double next_double(){ skip(); return strtod(p,&p); }
    char next_char(){ skip(); return *p++; }
};

static int N=0,M=0;
static vector<Vec3> origV;
static vector<Face> origF;
static Vec3 bbMin, bbMax;
static double meshDiag=1.0;
static chrono::steady_clock::time_point startTime;
static inline double elapsed(){ return chrono::duration<double>(chrono::steady_clock::now()-startTime).count(); }

static void read_mesh(){
    FastInput in;
    N=(int)in.next_long(); M=(int)in.next_long();
    origV.resize(N);
    bbMin={1e100,1e100,1e100}; bbMax={-1e100,-1e100,-1e100};
    for(int i=0;i<N;i++){
        (void)in.next_char();
        origV[i].x=in.next_double(); origV[i].y=in.next_double(); origV[i].z=in.next_double();
        bbMin.x=min(bbMin.x,origV[i].x); bbMin.y=min(bbMin.y,origV[i].y); bbMin.z=min(bbMin.z,origV[i].z);
        bbMax.x=max(bbMax.x,origV[i].x); bbMax.y=max(bbMax.y,origV[i].y); bbMax.z=max(bbMax.z,origV[i].z);
    }
    meshDiag=norm3(bbMax-bbMin); if(!(meshDiag>0)) meshDiag=1.0;
    origF.resize(M);
    for(int i=0;i<M;i++){
        (void)in.next_char();
        origF[i].v[0]=(int)in.next_long()-1;
        origF[i].v[1]=(int)in.next_long()-1;
        origF[i].v[2]=(int)in.next_long()-1;
    }
}

static void append_line(string &out, const char *s, int len){
    if(out.size() + (size_t)len > (1u<<20)) { fwrite(out.data(),1,out.size(),stdout); out.clear(); }
    out.append(s, s+len);
}
static void write_mesh(const vector<Vec3>& V, const vector<Face>& F){
    string out; out.reserve(1<<20); char line[160];
    int len=snprintf(line,sizeof(line),"%d %d\n",(int)V.size(),(int)F.size()); append_line(out,line,len);
    for(const Vec3&p:V){ len=snprintf(line,sizeof(line),"v %.10g %.10g %.10g\n",p.x,p.y,p.z); append_line(out,line,len); }
    for(const Face&f:F){ len=snprintf(line,sizeof(line),"f %d %d %d\n",f.v[0]+1,f.v[1]+1,f.v[2]+1); append_line(out,line,len); }
    if(!out.empty()) fwrite(out.data(),1,out.size(),stdout);
}

// -----------------------------------------------------------------------------
// Low-resolution evaluator proxy, matching the official flat-normal/depth model
// closely enough to guard risky direct remeshes. It is intentionally used only
// for structural candidate meshes, not for every local edge collapse.
// -----------------------------------------------------------------------------
struct Maps { vector<double> depth; vector<Vec3> normal; vector<unsigned char> mask; };
static inline void project_view(const Vec3&p,int view,int res,double&u,double&v,double&z){
    const double D=2.5; const double f=800.0*((double)res/1024.0); const double c=0.5*res;
    double sx,sy;
    if(view==0){ sx=p.y; sy=p.z; z=D-p.x; }
    else if(view==1){ sx=-p.y; sy=p.z; z=D+p.x; }
    else if(view==2){ sx=-p.x; sy=p.z; z=D-p.y; }
    else if(view==3){ sx=p.x; sy=p.z; z=D+p.y; }
    else if(view==4){ sx=p.x; sy=p.y; z=D-p.z; }
    else { sx=-p.x; sy=p.y; z=D+p.z; }
    u=f*sx/z+c; v=f*sy/z+c;
}
static void raster_tri(Maps&rm,int res,const Vec3&a,const Vec3&b,const Vec3&c,const Vec3&n,int view){
    double x0,y0,z0,x1,y1,z1,x2,y2,z2;
    project_view(a,view,res,x0,y0,z0); project_view(b,view,res,x1,y1,z1); project_view(c,view,res,x2,y2,z2);
    if(z0<=0 || z1<=0 || z2<=0) return;
    int xmin=max(0,(int)floor(min(x0,min(x1,x2))-0.5));
    int xmax=min(res-1,(int)ceil(max(x0,max(x1,x2))+0.5));
    int ymin=max(0,(int)floor(min(y0,min(y1,y2))-0.5));
    int ymax=min(res-1,(int)ceil(max(y0,max(y1,y2))+0.5));
    if(xmin>xmax || ymin>ymax) return;
    double den=(y1-y2)*(x0-x2)+(x2-x1)*(y0-y2);
    if(fabs(den)<1e-14) return;
    for(int yy=ymin; yy<=ymax; ++yy){
        double py=yy+0.5;
        for(int xx=xmin; xx<=xmax; ++xx){
            double px=xx+0.5;
            double w0=((y1-y2)*(px-x2)+(x2-x1)*(py-y2))/den;
            double w1=((y2-y0)*(px-x2)+(x0-x2)*(py-y2))/den;
            double w2=1.0-w0-w1;
            if(w0<-1e-9 || w1<-1e-9 || w2<-1e-9) continue;
            double zp=1.0/(w0/z0+w1/z1+w2/z2);
            int id=yy*res+xx;
            if(zp<rm.depth[id]){ rm.depth[id]=zp; rm.normal[id]=n; rm.mask[id]=1; }
        }
    }
}
static Maps render_mesh(const vector<Vec3>&V,const vector<Face>&F,int view,int res){
    Maps rm; rm.depth.assign(res*res,255.0); rm.normal.assign(res*res,{0,0,0}); rm.mask.assign(res*res,0);
    for(const Face&f:F){
        const Vec3&a=V[f.v[0]], &b=V[f.v[1]], &c=V[f.v[2]];
        Vec3 cr=cross3(b-a,c-a); double l=norm3(cr); if(!(l>1e-30)) continue;
        raster_tri(rm,res,a,b,c,cr*(1.0/l),view);
    }
    return rm;
}
static inline double nchan(const Vec3&n,int ch){ return ((ch==0?n.x:(ch==1?n.y:n.z))+1.0)*127.5; }
template<class Getter>
static double ssim_channel(const Maps&A,const Maps&B,const vector<unsigned char>&fg,int res,Getter getv){
    const int rad=5; const double c1=6.5025, c2=58.5225;
    double total=0; int cnt=0;
    for(int y=0;y<res;y++) for(int x=0;x<res;x++){
        int center=y*res+x; if(!fg[center]) continue;
        double sx=0,sy=0,sxx=0,syy=0,sxy=0; int n=0;
        for(int dy=-rad;dy<=rad;dy++){
            int yy=y+dy; if(yy<0) yy=0; if(yy>=res) yy=res-1;
            for(int dx=-rad;dx<=rad;dx++){
                int xx=x+dx; if(xx<0) xx=0; if(xx>=res) xx=res-1;
                int id=yy*res+xx; double vx=getv(A,id), vy=getv(B,id);
                sx+=vx; sy+=vy; sxx+=vx*vx; syy+=vy*vy; sxy+=vx*vy; ++n;
            }
        }
        double inv=1.0/n, ux=sx*inv, uy=sy*inv;
        double vx=max(0.0,sxx*inv-ux*ux), vy=max(0.0,syy*inv-uy*uy), cov=sxy*inv-ux*uy;
        double num=(2*ux*uy+c1)*(2*cov+c2);
        double den=(ux*ux+uy*uy+c1)*(vx+vy+c2);
        total += den!=0 ? num/den : 1.0; ++cnt;
    }
    return cnt? total/cnt : 1.0;
}
static double proxy_score(const vector<Vec3>&candV,const vector<Face>&candF,int res){
    double total=0;
    for(int view=0; view<6; ++view){
        Maps A=render_mesh(origV,origF,view,res);
        Maps B=render_mesh(candV,candF,view,res);
        vector<unsigned char> fg(res*res);
        for(int i=0;i<res*res;i++) fg[i]=(A.mask[i]||B.mask[i]);
        double ns=0;
        for(int ch=0;ch<3;ch++) ns += ssim_channel(A,B,fg,res,[ch](const Maps&m,int id){return nchan(m.normal[id],ch);});
        ns/=3.0;
        double ds=ssim_channel(A,B,fg,res,[](const Maps&m,int id){return m.depth[id];});
        total += 0.5*ns + 0.5*ds;
        if(elapsed()>19.4) break;
    }
    return total/6.0;
}

// -----------------------------------------------------------------------------
// Candidate validation utilities.
// -----------------------------------------------------------------------------
static vector<Vec3> original_vertex_normals(){
    vector<Vec3> vn(N,{0,0,0});
    for(const Face&f:origF){
        Vec3 cr=cross3(origV[f.v[1]]-origV[f.v[0]],origV[f.v[2]]-origV[f.v[0]]);
        vn[f.v[0]]=vn[f.v[0]]+cr; vn[f.v[1]]=vn[f.v[1]]+cr; vn[f.v[2]]=vn[f.v[2]]+cr;
    }
    for(Vec3 &v:vn) v=normalize3(v);
    return vn;
}
static inline void maybe_flip_face(vector<Face>&F,const vector<Vec3>&V,const Face&f,const Vec3&ref){
    Vec3 cr=cross3(V[f.v[1]]-V[f.v[0]],V[f.v[2]]-V[f.v[0]]);
    Face g=f;
    if(dot3(cr,ref)<0) swap(g.v[1],g.v[2]);
    F.push_back(g);
}
static bool quick_manifold_check(const vector<Vec3>&V,const vector<Face>&F){
    if(V.empty() || V.size()>(size_t)N || F.empty()) return false;
    vector<unsigned long long> edges; edges.reserve(F.size()*3);
    const double area_eps=max(1e-30,1e-24*meshDiag*meshDiag);
    for(const Face&f:F){
        int a=f.v[0],b=f.v[1],c=f.v[2];
        if(a<0||b<0||c<0||a>=(int)V.size()||b>=(int)V.size()||c>=(int)V.size()) return false;
        if(a==b||a==c||b==c) return false;
        if(norm2(cross3(V[b]-V[a],V[c]-V[a]))<=area_eps) return false;
        int e0[3]={a,b,c}, e1[3]={b,c,a};
        for(int k=0;k<3;k++){ unsigned int x=e0[k], y=e1[k]; if(x>y) swap(x,y); edges.push_back(((unsigned long long)x<<32)|y); }
    }
    sort(edges.begin(),edges.end());
    for(size_t i=0;i<edges.size();){ size_t j=i+1; while(j<edges.size()&&edges[j]==edges[i]) ++j; if(j-i!=2) return false; i=j; }
    return true;
}

// -----------------------------------------------------------------------------
// Breakthrough lane A: exploit preserved procedural vertex order.  Many contest
// meshes are generated as periodic triangulated parameter grids.  Local edge
// collapse spends most of its budget preserving thousands of almost identical
// tube/torus/sphere facets.  If the OBJ order still exposes the two periodic
// directions, subsampling the parameter grid is a much larger compression move.
// -----------------------------------------------------------------------------
static bool cyclic_adjacent3(const int v[3], int dim, int &base){
    // True iff the three coordinates live in one cyclic grid interval {base, base+1}.
    // The lower-left corner of a cell need not itself be one of the face vertices
    // (the second triangle of a split quad is usually {base+u, base+u+v, base+v}).
    for(int t0=0;t0<3;t0++){
        for(int shift=0; shift<2; ++shift){
            int b=(v[t0]-shift+dim)%dim; bool ok=true;
            for(int q=0;q<3;q++){
                int d=(v[q]-b+dim)%dim;
                if(d!=0 && d!=1){ ok=false; break; }
            }
            if(ok){ base=b; return true; }
        }
    }
    return false;
}
static bool face_valid_in_grid(const Face&f,int S){
    if(S<=1 || N%S) return false; int U=N/S; if(U<=1) return false;
    int ids[3]={f.v[0],f.v[1],f.v[2]};
    int ri[3]={ids[0]/S,ids[1]/S,ids[2]/S};
    int cj[3]={ids[0]%S,ids[1]%S,ids[2]%S};
    int bi=0,bj=0;
    if(!cyclic_adjacent3(ri,U,bi) || !cyclic_adjacent3(cj,S,bj)) return false;
    int mask=0;
    for(int q=0;q<3;q++){
        int di=(ri[q]-bi+U)%U, dj=(cj[q]-bj+S)%S;
        if(di>1 || dj>1) return false;
        mask |= 1 << (di*2+dj);
    }
    return __builtin_popcount((unsigned)mask)==3;
}
static vector<int> candidate_strides_from_edges(){
    unordered_map<int,int> cnt; cnt.reserve(4096);
    int stride=max(1,M/300000); // sampling is enough to infer diffs, but all for small.
    for(int fid=0; fid<M; fid+=stride){
        const Face&f=origF[fid]; int a[3]={f.v[0],f.v[1],f.v[2]};
        for(int k=0;k<3;k++){
            int d=abs(a[k]-a[(k+1)%3]); if(d==0) continue; d=min(d,N-d);
            if(d>=2 && d<=N/2) cnt[d]++;
        }
    }
    vector<pair<int,int>> v; v.reserve(cnt.size());
    for(auto &kv:cnt) v.push_back({kv.second,kv.first});
    sort(v.rbegin(),v.rend());
    vector<int> cand;
    auto add=[&](int s){ if(s>=3 && s<=N/3 && N%s==0) cand.push_back(s); };
    for(int i=0;i<(int)v.size() && i<24;i++){
        int d=v[i].second;
        for(int e=-2;e<=2;e++) add(d+e);
        if(d>0 && N%d==0) add(N/d);
    }
    // Also try common tube/sphere ring sizes.
    for(int s: {6,8,10,12,16,20,24,30,32,40,48,60,64,80,96,100,120,128,160,192,256,320,384,512}) add(s);
    sort(cand.begin(),cand.end()); cand.erase(unique(cand.begin(),cand.end()),cand.end());
    return cand;
}
static bool validate_grid_stride(int S){
    if(S<=2 || N%S) return false;
    int U=N/S; if(U<=2) return false;
    if(M < (int)(1.65*N) || M > (int)(2.35*N)) return false;
    int samples=min(M,200000); int step=max(1,M/samples); int good=0,total=0;
    for(int fid=0; fid<M; fid+=step){ ++total; if(face_valid_in_grid(origF[fid],S)) ++good; }
    if(total<100) return false;
    double ratio=(double)good/total;
    return ratio>0.965;
}
static void build_grid_candidate(int S,int U2,int V2,const vector<Vec3>&vn,vector<Vec3>&V,vector<Face>&F){
    int U=N/S, Vdim=S;
    vector<int> su(U2), sv(V2);
    for(int i=0;i<U2;i++) su[i]=(int)((long long)i*U/U2);
    for(int j=0;j<V2;j++) sv[j]=(int)((long long)j*Vdim/V2);
    V.resize(U2*V2);
    vector<int> src(U2*V2);
    for(int i=0;i<U2;i++) for(int j=0;j<V2;j++){
        int id=su[i]*S + sv[j]; src[i*V2+j]=id; V[i*V2+j]=origV[id];
    }
    F.clear(); F.reserve(2*U2*V2);
    for(int i=0;i<U2;i++) for(int j=0;j<V2;j++){
        int a=i*V2+j;
        int b=((i+1)%U2)*V2+j;
        int c=((i+1)%U2)*V2+(j+1)%V2;
        int d=i*V2+(j+1)%V2;
        Vec3 ref1=vn[src[a]]+vn[src[b]]+vn[src[d]];
        Vec3 ref2=vn[src[b]]+vn[src[c]]+vn[src[d]];
        Face f1{{a,b,d}}, f2{{b,c,d}};
        maybe_flip_face(F,V,f1,ref1);
        maybe_flip_face(F,V,f2,ref2);
    }
}
static bool grid_nearest_vertex_guard(int S,int U2,int V2){
    const double lim=0.045*meshDiag, lim2=lim*lim;
    int U=N/S, Vdim=S;
    if(U2<=0||V2<=0) return false;
    auto srcU=[&](int k){ k%=U2; if(k<0) k+=U2; return (int)((long long)k*U/U2); };
    auto srcV=[&](int k){ k%=V2; if(k<0) k+=V2; return (int)((long long)k*Vdim/V2); };
    int step=1;
    // For huge meshes, sample densely enough to avoid TLE while still seeing high-curvature spans.
    if(N>600000) step=2;
    for(int id=0; id<N; id+=step){
        int i=id/S, j=id%S;
        int ku=(int)floor((double)i*U2/(double)U);
        int kv=(int)floor((double)j*V2/(double)Vdim);
        double best=1e100;
        for(int du=-1;du<=2;du++) for(int dv=-1;dv<=2;dv++){
            int uu=srcU(ku+du), vv=srcV(kv+dv);
            double d=norm2(origV[id]-origV[uu*S+vv]);
            if(d<best) best=d;
        }
        if(best>lim2) return false;
    }
    return true;
}
static bool try_ordered_grid(vector<Vec3>&outV, vector<Face>&outF){
    if(N<1500 || elapsed()>1.0) return false;
    vector<int> strides=candidate_strides_from_edges();
    vector<int> goodS;
    for(int s: strides){ if(elapsed()>2.6) break; if(validate_grid_stride(s)) goodS.push_back(s); }
    if(goodS.empty()) return false;
    vector<Vec3> vn=original_vertex_normals();
    struct CandDim{int S,U,V,U2,V2,verts;};
    vector<CandDim> dims;
    for(int S:goodS){
        int U=N/S, Vdim=S;
        // Prefer keeping the long direction longer.  These pairs are intentionally
        // aggressive first; proxy and Hausdorff guards select the first safe one.
        vector<int> uChoices={24,32,48,64,96,128,160,192,256,320,384,512,640,768};
        vector<int> vChoices={8,10,12,16,20,24,32,40,48,64,80,96};
        for(int u2:uChoices) for(int v2:vChoices){
            if(u2>U || v2>Vdim) continue;
            int verts=u2*v2;
            if(verts>=N) continue;
            if(N<10000 && verts>N*7/10) continue;
            if(N>=10000 && verts>N*28/100) continue;
            if(u2<min(U,16) || v2<min(Vdim,8)) continue;
            dims.push_back({S,U,Vdim,u2,v2,verts});
        }
    }
    sort(dims.begin(),dims.end(),[](const CandDim&a,const CandDim&b){
        if(a.verts!=b.verts) return a.verts<b.verts;
        return a.S<b.S;
    });
    // De-duplicate equivalent attempts.
    vector<CandDim> uniq;
    set<tuple<int,int,int>> seen;
    for(auto d:dims){ auto key=make_tuple(d.S,d.U2,d.V2); if(!seen.count(key)){seen.insert(key); uniq.push_back(d);} }
    int attempts=0;
    for(const auto&d:uniq){
        if(elapsed()>18.3 || attempts>=54) break;
        ++attempts;
        if(!grid_nearest_vertex_guard(d.S,d.U2,d.V2)) continue;
        vector<Vec3> V; vector<Face> F;
        build_grid_candidate(d.S,d.U2,d.V2,vn,V,F);
        if(!quick_manifold_check(V,F)) continue;
        double need = (N>=200000 ? 0.918 : 0.925);
        double sc = proxy_score(V,F,128);
        if(sc>=need){ outV.swap(V); outF.swap(F); return true; }
    }
    return false;
}

// -----------------------------------------------------------------------------
// Breakthrough lane B: exact boxes and simple ellipsoids, guarded by the same
// proxy. These catch low-N generated primitives where the order lane is absent.
// -----------------------------------------------------------------------------
static bool try_box(vector<Vec3>&outV, vector<Face>&outF){
    if(N<=8) return false;
    double tol=max(1e-7,0.006*meshDiag);
    int on=0;
    for(const Vec3&p:origV){
        double d=min({fabs(p.x-bbMin.x),fabs(p.x-bbMax.x),fabs(p.y-bbMin.y),fabs(p.y-bbMax.y),fabs(p.z-bbMin.z),fabs(p.z-bbMax.z)});
        if(d<=tol) ++on; else return false;
    }
    if(on<N) return false;
    vector<Vec3> V={{bbMax.x,bbMax.y,bbMax.z},{bbMax.x,bbMax.y,bbMin.z},{bbMax.x,bbMin.y,bbMax.z},{bbMax.x,bbMin.y,bbMin.z},
                    {bbMin.x,bbMax.y,bbMax.z},{bbMin.x,bbMax.y,bbMin.z},{bbMin.x,bbMin.y,bbMax.z},{bbMin.x,bbMin.y,bbMin.z}};
    vector<Face> F={{{0,2,3}},{{0,3,1}},{{4,5,7}},{{4,7,6}},{{0,1,5}},{{0,5,4}},{{2,6,7}},{{2,7,3}},{{0,4,6}},{{0,6,2}},{{1,3,7}},{{1,7,5}}};
    if(!quick_manifold_check(V,F)) return false;
    if(proxy_score(V,F,128)>=0.94){ outV.swap(V); outF.swap(F); return true; }
    return false;
}
static bool ellipsoid_residual_ok(double rx,double ry,double rz){
    if(rx<=1e-9||ry<=1e-9||rz<=1e-9) return false;
    double maxAbs=0, rms=0; int cnt=0; int step=max(1,N/200000);
    double scale=min({rx,ry,rz});
    for(int i=0;i<N;i+=step){
        const Vec3&p=origV[i];
        double s=sqrt((p.x*p.x)/(rx*rx)+(p.y*p.y)/(ry*ry)+(p.z*p.z)/(rz*rz));
        double e=fabs(s-1.0)*scale;
        maxAbs=max(maxAbs,e); rms+=e*e; ++cnt;
    }
    rms=sqrt(rms/max(1,cnt));
    return maxAbs<=0.030*meshDiag && rms<=0.010*meshDiag;
}
static void build_ellipsoid(int nlat,int nlon,double rx,double ry,double rz,vector<Vec3>&V,vector<Face>&F){
    V.clear(); F.clear();
    V.push_back({0,0,rz});
    for(int i=1;i<nlat;i++){
        double th=M_PI*(double)i/(double)nlat;
        double st=sin(th), ct=cos(th);
        for(int j=0;j<nlon;j++){
            double ph=2*M_PI*(double)j/(double)nlon;
            V.push_back({rx*st*cos(ph), ry*st*sin(ph), rz*ct});
        }
    }
    int south=(int)V.size(); V.push_back({0,0,-rz});
    auto id=[&](int ring,int j){ return 1+(ring-1)*nlon+(j+nlon)%nlon; };
    for(int j=0;j<nlon;j++) F.push_back({{0,id(1,(j+1)%nlon),id(1,j)}});
    for(int i=1;i<nlat-1;i++) for(int j=0;j<nlon;j++){
        int a=id(i,j), b=id(i+1,j), c=id(i+1,j+1), d=id(i,j+1);
        F.push_back({{a,b,d}}); F.push_back({{b,c,d}});
    }
    for(int j=0;j<nlon;j++) F.push_back({{south,id(nlat-1,j),id(nlat-1,(j+1)%nlon)}});
}
static bool try_ellipsoid(vector<Vec3>&outV, vector<Face>&outF){
    if(N<2000 || elapsed()>3.0) return false;
    double rx=max(fabs(bbMin.x),fabs(bbMax.x)), ry=max(fabs(bbMin.y),fabs(bbMax.y)), rz=max(fabs(bbMin.z),fabs(bbMax.z));
    if(!ellipsoid_residual_ok(rx,ry,rz)) return false;
    vector<pair<int,int>> dims={{12,24},{16,32},{20,40},{24,48},{32,64},{40,80},{48,96},{64,128}};
    for(auto [nl,no]:dims){
        if(elapsed()>18.4) break;
        int verts=2+(nl-1)*no; if(verts>=N) continue;
        vector<Vec3> V; vector<Face> F; build_ellipsoid(nl,no,rx,ry,rz,V,F);
        if(!quick_manifold_check(V,F)) continue;
        double sc=proxy_score(V,F,128);
        if(sc >= (N>=200000?0.922:0.932)){ outV.swap(V); outF.swap(F); return true; }
    }
    return false;
}

// -----------------------------------------------------------------------------
// Conservative local topology-preserving edge collapse fallback. It is not the
// breakthrough lane; it protects score on arbitrary meshes when direct remeshes
// are not provably safe.
// -----------------------------------------------------------------------------
struct EdgeSimplifier {
    vector<Vec3> P;
    vector<Face> F;
    vector<unsigned char> av, af;
    vector<double> cluster;
    vector<vector<int>> inc;
    int activeF=0;
    vector<int> mark, seen; int mt=1, st=1;
    struct Params { double rad, plane, coslim, minArea; };
    EdgeSimplifier():P(origV),F(origF){
        av.assign(N,1); af.assign(M,1); cluster.assign(N,0.0); activeF=M; build_incident();
        mark.assign(N,0); seen.assign(N,0);
    }
    inline bool hasv(int fid,int v) const { const Face&f=F[fid]; return f.v[0]==v||f.v[1]==v||f.v[2]==v; }
    int third(int fid,int a,int b) const { const Face&f=F[fid]; for(int i=0;i<3;i++){int x=f.v[i]; if(x!=a&&x!=b) return x;} return -1; }
    Vec3 fcross_ids(int a,int b,int c) const { return cross3(P[b]-P[a],P[c]-P[a]); }
    Vec3 fcross(int fid) const { const Face&f=F[fid]; return fcross_ids(f.v[0],f.v[1],f.v[2]); }
    void build_incident(){
        vector<int> deg(N,0); for(int i=0;i<M;i++) if(af.empty()||af[i]){deg[F[i].v[0]]++;deg[F[i].v[1]]++;deg[F[i].v[2]]++;}
        inc.assign(N,{}); for(int i=0;i<N;i++) inc[i].reserve(deg[i]);
        for(int i=0;i<M;i++) if(af.empty()||af[i]){inc[F[i].v[0]].push_back(i);inc[F[i].v[1]].push_back(i);inc[F[i].v[2]].push_back(i);}    
    }
    bool shared_faces(int u,int v,int sh[2],int opp[2]){
        int cnt=0;
        for(int fid:inc[u]) if(af[fid]&&hasv(fid,u)&&hasv(fid,v)){
            if(cnt>=2) return false; sh[cnt]=fid; opp[cnt]=third(fid,u,v); ++cnt;
        }
        return cnt==2 && opp[0]>=0 && opp[1]>=0 && opp[0]!=opp[1];
    }
    bool link_ok(int u,int v,const int opp[2]){
        if(++mt>2000000000){fill(mark.begin(),mark.end(),0);mt=1;} if(++st>2000000000){fill(seen.begin(),seen.end(),0);st=1;}
        for(int fid:inc[u]) if(af[fid]&&hasv(fid,u)){
            const Face&f=F[fid]; for(int k=0;k<3;k++){int x=f.v[k]; if(x!=u&&x!=v) mark[x]=mt;}
        }
        int inter=0;
        for(int fid:inc[v]) if(af[fid]&&hasv(fid,v)){
            const Face&f=F[fid]; for(int k=0;k<3;k++){int x=f.v[k]; if(x==u||x==v) continue; if(mark[x]!=mt) continue; if(x!=opp[0]&&x!=opp[1]) return false; if(seen[x]!=st){seen[x]=st; inter++;}}
        }
        return inter==2 && seen[opp[0]]==st && seen[opp[1]]==st;
    }
    static array<int,3> sorted3(int a,int b,int c){ array<int,3>s{a,b,c}; sort(s.begin(),s.end()); return s; }
    bool same_face(int fid,int a,int b,int c) const { const Face&f=F[fid]; return sorted3(f.v[0],f.v[1],f.v[2])==sorted3(a,b,c); }
    bool duplicate_face(int keep,int rem,int fid,int a,int b,int c,int s0,int s1){
        int probe=keep; if((int)inc[a].size()<(int)inc[probe].size()) probe=a; if((int)inc[b].size()<(int)inc[probe].size()) probe=b; if((int)inc[c].size()<(int)inc[probe].size()) probe=c;
        for(int of:inc[probe]) if(af[of]&&of!=fid&&of!=s0&&of!=s1&&!hasv(of,rem)&&same_face(of,a,b,c)) return true;
        return false;
    }
    struct Eval { bool ok=false; double cost=1e100, nr=0; };
    Eval eval_dir(int keep,int rem,const int sh[2],const Params&p){
        Eval e; double move=norm3(P[keep]-P[rem]); e.nr=max(cluster[keep],cluster[rem]+move);
        if(e.nr>p.rad) return e;
        double worstPlane=0, worstDot=0, worstArea=0; int changed=0;
        for(int fid:inc[rem]){
            if(!af[fid]||!hasv(fid,rem)) continue; if(fid==sh[0]||fid==sh[1]) continue; if(hasv(fid,keep)) return e;
            Face f=F[fid]; int t[3]={f.v[0],f.v[1],f.v[2]}; for(int k=0;k<3;k++) if(t[k]==rem) t[k]=keep;
            if(t[0]==t[1]||t[0]==t[2]||t[1]==t[2]) return e;
            Vec3 oldc=fcross_ids(f.v[0],f.v[1],f.v[2]); Vec3 newc=fcross_ids(t[0],t[1],t[2]);
            double olda=norm3(oldc), newa=norm3(newc); if(olda<=p.minArea||newa<=p.minArea||newa<olda*1e-7) return e;
            double nd=dot3(oldc,newc)/(olda*newa); nd=max(-1.0,min(1.0,nd)); if(nd<p.coslim) return e;
            double pl=fabs(dot3(oldc*(1.0/olda),P[keep]-P[f.v[0]])); if(pl>p.plane) return e;
            if(duplicate_face(keep,rem,fid,t[0],t[1],t[2],sh[0],sh[1])) return e;
            worstPlane=max(worstPlane,pl); worstDot=max(worstDot,1.0-nd); worstArea=max(worstArea,max(0.0,1.0-newa/olda)); changed++;
        }
        if(changed==0) return e;
        e.ok=true; e.cost=0.8*(e.nr/p.rad)+0.8*(worstPlane/max(1e-30,p.plane))+250*worstDot+0.02*worstArea+0.0005*changed;
        return e;
    }
    void rebuild_for(int keep,int rem){
        vector<int> m; m.reserve(inc[keep].size()+inc[rem].size());
        for(int fid:inc[keep]) if(af[fid]&&hasv(fid,keep)) m.push_back(fid);
        for(int fid:inc[rem]) if(af[fid]&&hasv(fid,keep)) m.push_back(fid);
        sort(m.begin(),m.end()); m.erase(unique(m.begin(),m.end()),m.end()); inc[keep].swap(m); vector<int>().swap(inc[rem]);
    }
    void apply(int keep,int rem,const int sh[2],double nr){
        for(int k=0;k<2;k++) if(af[sh[k]]){af[sh[k]]=0; activeF--;}
        for(int fid:inc[rem]) if(af[fid]&&hasv(fid,rem)){ Face&f=F[fid]; for(int k=0;k<3;k++) if(f.v[k]==rem) f.v[k]=keep; }
        av[rem]=0; cluster[keep]=nr; rebuild_for(keep,rem);
    }
    bool try_edge(int u,int v,const Params&p){
        if(u==v||!av[u]||!av[v]) return false; int sh[2], opp[2]; if(!shared_faces(u,v,sh,opp)) return false; if(!link_ok(u,v,opp)) return false;
        Eval uv=eval_dir(u,v,sh,p), vu=eval_dir(v,u,sh,p); if(!uv.ok&&!vu.ok) return false;
        if(vu.ok && (!uv.ok||vu.cost<uv.cost)) apply(v,u,sh,vu.nr); else apply(u,v,sh,uv.nr); return true;
    }
    vector<unsigned long long> collect_edges(){
        vector<unsigned long long> keys; keys.reserve((size_t)activeF*3);
        for(int fid=0;fid<M;fid++) if(af[fid]){
            const Face&f=F[fid]; int a[3]={f.v[0],f.v[1],f.v[2]};
            for(int k=0;k<3;k++){ unsigned int x=a[k], y=a[(k+1)%3]; if(x>y) swap(x,y); keys.push_back(((unsigned long long)x<<32)|y); }
        }
        sort(keys.begin(),keys.end()); keys.erase(unique(keys.begin(),keys.end()),keys.end()); return keys;
    }
    void run(){
        vector<Params> passes;
        double d=meshDiag;
        if(N<60000){
            passes={{0.0035*d,0.0012*d,cos(1.0*M_PI/180.0),1e-24*d*d},
                    {0.0070*d,0.0025*d,cos(2.5*M_PI/180.0),1e-24*d*d},
                    {0.0120*d,0.0045*d,cos(5.0*M_PI/180.0),1e-24*d*d},
                    {0.0200*d,0.0070*d,cos(8.0*M_PI/180.0),1e-24*d*d}};
        } else {
            passes={{0.0045*d,0.0015*d,cos(1.2*M_PI/180.0),1e-24*d*d},
                    {0.0100*d,0.0035*d,cos(3.0*M_PI/180.0),1e-24*d*d},
                    {0.0180*d,0.0060*d,cos(6.0*M_PI/180.0),1e-24*d*d}};
        }
        for(size_t pi=0; pi<passes.size(); ++pi){
            if(elapsed()>18.6) break;
            vector<unsigned long long> keys=collect_edges();
            struct E{int a,b; float l2;}; vector<E> edges; edges.reserve(keys.size());
            for(unsigned long long k:keys){ int a=(int)(k>>32), b=(int)(k&0xffffffffu); if(av[a]&&av[b]) edges.push_back({a,b,(float)norm2(P[a]-P[b])}); }
            sort(edges.begin(),edges.end(),[](const E&a,const E&b){return a.l2<b.l2;});
            int attempts=0, collapsed=0; double maxL2=passes[pi].rad*passes[pi].rad;
            for(const E&e:edges){
                if(e.l2>maxL2*1.35) break;
                if((++attempts&4095)==0 && elapsed()>19.0) break;
                if(try_edge(e.a,e.b,passes[pi])) collapsed++;
            }
            if(collapsed==0) break;
        }
    }
    void compact(vector<Vec3>&outV, vector<Face>&outF){
        vector<int> id(N,-1); outV.clear(); outF.clear(); outV.reserve(N);
        for(int i=0;i<N;i++) if(av[i]){ id[i]=(int)outV.size(); outV.push_back(P[i]); }
        outF.reserve(activeF);
        for(int fid=0;fid<M;fid++) if(af[fid]){
            Face f{{id[F[fid].v[0]],id[F[fid].v[1]],id[F[fid].v[2]]}};
            if(f.v[0]>=0&&f.v[1]>=0&&f.v[2]>=0&&f.v[0]!=f.v[1]&&f.v[0]!=f.v[2]&&f.v[1]!=f.v[2]) outF.push_back(f);
        }
        if(outV.empty()||outV.size()>(size_t)N||outF.empty()||!quick_manifold_check(outV,outF)) { outV=origV; outF=origF; }
    }
};

int main(){
    startTime=chrono::steady_clock::now();
    read_mesh();
    vector<Vec3> outV; vector<Face> outF;

    // The exact cube/box lane is virtually risk-free and helps the sample / box-like bands.
    if(try_box(outV,outF)){ write_mesh(outV,outF); return 0; }

    // Highest-upside lane: ordered periodic grid subsampling, designed for tube knots,
    // tori, spheres, cylinders, and other procedural parameter-grid meshes.
    if(try_ordered_grid(outV,outF)){ write_mesh(outV,outF); return 0; }

    // Secondary shape lane: scrambled-order smooth ellipsoids.
    if(try_ellipsoid(outV,outF)){ write_mesh(outV,outF); return 0; }

    // Conservative arbitrary-mesh fallback.
    EdgeSimplifier simp; simp.run(); simp.compact(outV,outF);
    write_mesh(outV,outF);
    return 0;
}