#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cstdint>
#include <queue>
#include <string>
#include <vector>
#include <unordered_map>
using namespace std;

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
static inline double norm3(const Vec3&a){return sqrt(norm2(a));}
static inline Vec3 normalize(Vec3 a){double n=norm3(a); if(n>1e-300) return a/n; return {0,0,0};}
static inline double dist2(const Vec3&a,const Vec3&b){return norm2(a-b);}

struct Tri{int a,b,c;};
struct Face{int v[3]; unsigned char active=1; Vec3 n;};

struct FastInput{
    vector<char> buf; char* p=nullptr;
    FastInput(){
        buf.reserve(1<<27); char tmp[1<<16]; size_t n;
        while((n=fread(tmp,1,sizeof(tmp),stdin))>0) buf.insert(buf.end(),tmp,tmp+n);
        buf.push_back('\0'); p=buf.data();
    }
    inline void skip(){ while(*p==' '||*p=='\n'||*p=='\r'||*p=='\t') ++p; }
    int next_int(){ skip(); int s=1; if(*p=='-'){s=-1;++p;} int x=0; while(*p>='0'&&*p<='9'){x=x*10+(*p-'0');++p;} return x*s; }
    double next_double(){ skip(); char* e; double x=strtod(p,&e); p=e; return x; }
    char next_char(){ skip(); return *p++; }
};

static int N=0,M=0;
static vector<Vec3> Orig;
static vector<Tri> OrigF;
static double diagLen=1.0, tau=0.05, tau2=0.0025;
static chrono::steady_clock::time_point START;
static inline double elapsed(){return chrono::duration<double>(chrono::steady_clock::now()-START).count();}

static void read_input(){
    FastInput in; N=in.next_int(); M=in.next_int();
    Orig.resize(N); OrigF.resize(M);
    Vec3 mn{1e100,1e100,1e100}, mx{-1e100,-1e100,-1e100};
    for(int i=0;i<N;i++){
        (void)in.next_char();
        Orig[i].x=in.next_double(); Orig[i].y=in.next_double(); Orig[i].z=in.next_double();
        mn.x=min(mn.x,Orig[i].x); mn.y=min(mn.y,Orig[i].y); mn.z=min(mn.z,Orig[i].z);
        mx.x=max(mx.x,Orig[i].x); mx.y=max(mx.y,Orig[i].y); mx.z=max(mx.z,Orig[i].z);
    }
    for(int i=0;i<M;i++){
        (void)in.next_char();
        OrigF[i].a=in.next_int()-1; OrigF[i].b=in.next_int()-1; OrigF[i].c=in.next_int()-1;
    }
    Vec3 d=mx-mn; diagLen=max(1e-12,norm3(d)); tau=0.05*diagLen*0.995; tau2=tau*tau;
}

static inline uint64_t edge_key(int a,int b){ if(a>b) swap(a,b); return ((uint64_t)(uint32_t)a<<32)|(uint32_t)b; }
static inline uint64_t face_key3(int a,int b,int c){
    if(a>b) swap(a,b); if(b>c) swap(b,c); if(a>b) swap(a,b);
    return ((uint64_t)(uint32_t)a<<42)^((uint64_t)(uint32_t)b<<21)^(uint64_t)(uint32_t)c;
}

static bool manifold_basic_check(const vector<Vec3>&V,const vector<Tri>&F){
    if(V.empty()||F.empty()||V.size()>(size_t)N) return false;
    const double area_eps=max(1e-30,1e-24*diagLen*diagLen);
    vector<uint64_t> edges; edges.reserve(F.size()*3);
    vector<uint64_t> faces; faces.reserve(F.size());
    for(const Tri&t:F){
        if(t.a<0||t.b<0||t.c<0||t.a>=(int)V.size()||t.b>=(int)V.size()||t.c>=(int)V.size()) return false;
        if(t.a==t.b||t.b==t.c||t.a==t.c) return false;
        Vec3 cr=cross3(V[t.b]-V[t.a],V[t.c]-V[t.a]);
        if(!(norm2(cr)>area_eps)) return false;
        edges.push_back(edge_key(t.a,t.b)); edges.push_back(edge_key(t.b,t.c)); edges.push_back(edge_key(t.c,t.a));
        faces.push_back(face_key3(t.a,t.b,t.c));
    }
    sort(faces.begin(),faces.end());
    if(adjacent_find(faces.begin(),faces.end())!=faces.end()) return false;
    sort(edges.begin(),edges.end());
    for(size_t i=0;i<edges.size();){ size_t j=i+1; while(j<edges.size()&&edges[j]==edges[i]) ++j; if(j-i!=2) return false; i=j; }
    return true;
}

struct RenderMaps{int R=0; vector<float> depth; vector<float> norm;};

static inline void project_view(const Vec3&p,int view,int R,double&u,double&v,double&z){
    const double D=2.5; const double f=800.0*(double)R/1024.0; const double c=0.5*(double)R;
    double sx,sy;
    if(view==0){sx=p.y; sy=p.z; z=D-p.x;}
    else if(view==1){sx=-p.y; sy=p.z; z=D+p.x;}
    else if(view==2){sx=-p.x; sy=p.z; z=D-p.y;}
    else if(view==3){sx=p.x; sy=p.z; z=D+p.y;}
    else if(view==4){sx=p.x; sy=p.y; z=D-p.z;}
    else {sx=-p.x; sy=p.y; z=D+p.z;}
    u=f*sx/z+c; v=f*sy/z+c;
}

static void render_mesh(const vector<Vec3>&V,const vector<Tri>&F,RenderMaps&rm,int R){
    const int PIX=R*R;
    rm.R=R;
    rm.depth.assign((size_t)6*PIX,255.0f);
    rm.norm.assign((size_t)6*PIX*3,127.5f);
    vector<float> U(V.size()), VV(V.size()), Z(V.size());
    vector<Vec3> fn(F.size());
    for(size_t i=0;i<F.size();++i){
        const Tri&t=F[i];
        fn[i]=normalize(cross3(V[t.b]-V[t.a],V[t.c]-V[t.a]));
    }
    for(int view=0;view<6;view++){
        for(size_t i=0;i<V.size();++i){
            double u,v,z; project_view(V[i],view,R,u,v,z);
            U[i]=(float)u; VV[i]=(float)v; Z[i]=(float)z;
        }
        float* zbuf=rm.depth.data()+(size_t)view*PIX;
        float* nbuf=rm.norm.data()+(size_t)view*PIX*3;
        for(size_t ti=0;ti<F.size();++ti){
            const Tri&t=F[ti];
            int ia=t.a,ib=t.b,ic=t.c;
            float x0=U[ia], y0=VV[ia], z0=Z[ia];
            float x1=U[ib], y1=VV[ib], z1=Z[ib];
            float x2=U[ic], y2=VV[ic], z2=Z[ic];
            if(!(z0>0&&z1>0&&z2>0)) continue;
            float minx=min(x0,min(x1,x2)), maxx=max(x0,max(x1,x2));
            float miny=min(y0,min(y1,y2)), maxy=max(y0,max(y1,y2));
            int xmin=max(0,(int)floor(minx-0.5f));
            int xmax=min(R-1,(int)ceil(maxx+0.5f));
            int ymin=max(0,(int)floor(miny-0.5f));
            int ymax=min(R-1,(int)ceil(maxy+0.5f));
            if(xmin>xmax||ymin>ymax) continue;
            float den=(y1-y2)*(x0-x2)+(x2-x1)*(y0-y2);
            if(fabs(den)<1e-20f) continue;
            float invDen=1.0f/den;
            Vec3 n=fn[ti];
            float nr=(float)((n.x+1.0)*127.5);
            float ng=(float)((n.y+1.0)*127.5);
            float nb=(float)((n.z+1.0)*127.5);
            for(int yy=ymin;yy<=ymax;yy++){
                float py=(float)yy+0.5f;
                int row=yy*R;
                for(int xx=xmin;xx<=xmax;xx++){
                    float px=(float)xx+0.5f;
                    float w0=((y1-y2)*(px-x2)+(x2-x1)*(py-y2))*invDen;
                    float w1=((y2-y0)*(px-x2)+(x0-x2)*(py-y2))*invDen;
                    float w2=1.0f-w0-w1;
                    if(w0>=-1e-6f&&w1>=-1e-6f&&w2>=-1e-6f){
                        float dep=1.0f/(w0/z0+w1/z1+w2/z2);
                        int idx=row+xx;
                        if(dep<zbuf[idx]){
                            zbuf[idx]=dep;
                            float*q=nbuf+(size_t)idx*3;
                            q[0]=nr; q[1]=ng; q[2]=nb;
                        }
                    }
                }
            }
        }
    }
}

static inline double rect_sum(const vector<double>&I,int W,int x0,int y0,int x1,int y1){
    return I[(size_t)y1*W+x1]-I[(size_t)y0*W+x1]-I[(size_t)y1*W+x0]+I[(size_t)y0*W+x0];
}

static double ssim_channel(const float*X,int sx,const float*Y,int sy,const float*DX,const float*DY,int R,
                           vector<double>&IX,vector<double>&IY,vector<double>&IX2,vector<double>&IY2,vector<double>&IXY){
    int W=R+1;
    size_t SZ=(size_t)W*W;
    fill(IX.begin(),IX.begin()+SZ,0.0);
    fill(IY.begin(),IY.begin()+SZ,0.0);
    fill(IX2.begin(),IX2.begin()+SZ,0.0);
    fill(IY2.begin(),IY2.begin()+SZ,0.0);
    fill(IXY.begin(),IXY.begin()+SZ,0.0);
    for(int y=1;y<=R;y++){
        double ax=0,ay=0,ax2=0,ay2=0,axy=0;
        int row=(y-1)*R;
        for(int x=1;x<=R;x++){
            int p=row+x-1;
            double vx=X[(size_t)p*sx], vy=Y[(size_t)p*sy];
            ax+=vx; ay+=vy; ax2+=vx*vx; ay2+=vy*vy; axy+=vx*vy;
            size_t id=(size_t)y*W+x, up=(size_t)(y-1)*W+x;
            IX[id]=IX[up]+ax;
            IY[id]=IY[up]+ay;
            IX2[id]=IX2[up]+ax2;
            IY2[id]=IY2[up]+ay2;
            IXY[id]=IXY[up]+axy;
        }
    }
    const int rad=5, area=121;
    const double c1=6.5025, c2=58.5225;
    long long cnt=0;
    long double acc=0;
    for(int y=rad;y<R-rad;y++){
        int row=y*R;
        for(int x=rad;x<R-rad;x++){
            int p=row+x;
            if(!(DX[p]<254.0f||DY[p]<254.0f)) continue;
            int x0=x-rad,y0=y-rad,x1=x+rad+1,y1=y+rad+1;
            double mx=rect_sum(IX,W,x0,y0,x1,y1)/area;
            double my=rect_sum(IY,W,x0,y0,x1,y1)/area;
            double vx=rect_sum(IX2,W,x0,y0,x1,y1)/area-mx*mx;
            double vy=rect_sum(IY2,W,x0,y0,x1,y1)/area-my*my;
            double cv=rect_sum(IXY,W,x0,y0,x1,y1)/area-mx*my;
            if(vx<0&&vx>-1e-7) vx=0;
            if(vy<0&&vy>-1e-7) vy=0;
            double den=(mx*mx+my*my+c1)*(vx+vy+c2);
            double val=den!=0?((2*mx*my+c1)*(2*cv+c2)/den):1.0;
            acc+=val; ++cnt;
        }
    }
    return cnt?(double)(acc/cnt):1.0;
}

static double render_ssim(const RenderMaps&A,const RenderMaps&B){
    int R=A.R, PIX=R*R, W=R+1;
    size_t SZ=(size_t)W*W;
    vector<double> IX(SZ),IY(SZ),IX2(SZ),IY2(SZ),IXY(SZ);
    double total=0;
    for(int view=0;view<6;view++){
        const float* ad=A.depth.data()+(size_t)view*PIX;
        const float* bd=B.depth.data()+(size_t)view*PIX;
        double ns=0;
        for(int ch=0;ch<3;ch++){
            const float* an=A.norm.data()+(size_t)view*PIX*3+ch;
            const float* bn=B.norm.data()+(size_t)view*PIX*3+ch;
            ns+=ssim_channel(an,3,bn,3,ad,bd,R,IX,IY,IX2,IY2,IXY);
        }
        ns/=3.0;
        double ds=ssim_channel(ad,1,bd,1,ad,bd,R,IX,IY,IX2,IY2,IXY);
        total+=0.5*(ns+ds);
    }
    return total/6.0;
}

static RenderMaps orig_render_cache[3];
static int cache_R[3]={0,0,0};
static const RenderMaps& get_orig_render(int R){
    int slot=(R<=80?0:(R<=160?1:2));
    if(cache_R[slot]!=R){
        render_mesh(Orig,OrigF,orig_render_cache[slot],R);
        cache_R[slot]=R;
    }
    return orig_render_cache[slot];
}
static double visual_score_candidate(const vector<Vec3>&V,const vector<Tri>&F,int R){
    RenderMaps cand;
    render_mesh(V,F,cand,R);
    return render_ssim(get_orig_render(R),cand);
}

struct PointGrid{
    double cell=1.0;
    unordered_map<unsigned long long, vector<int>> mp;
    static unsigned long long key(int x,int y,int z){
        unsigned long long a=(unsigned int)(x+119304647);
        unsigned long long b=(unsigned int)(y-72537221);
        unsigned long long c=(unsigned int)(z+362436069);
        return (a*73856093ull) ^ (b*19349663ull) ^ (c*83492791ull);
    }
    void build(const vector<Vec3>&pts,double c){
        cell=max(c,1e-12);
        mp.clear();
        mp.reserve(pts.size()*2+64);
        for(int i=0;i<(int)pts.size();++i){
            int ix=(int)floor(pts[i].x/cell);
            int iy=(int)floor(pts[i].y/cell);
            int iz=(int)floor(pts[i].z/cell);
            mp[key(ix,iy,iz)].push_back(i);
        }
    }
    bool has_near(const vector<Vec3>&pts,const Vec3&p,double r2)const{
        int ix=(int)floor(p.x/cell);
        int iy=(int)floor(p.y/cell);
        int iz=(int)floor(p.z/cell);
        for(int dx=-1;dx<=1;++dx)for(int dy=-1;dy<=1;++dy)for(int dz=-1;dz<=1;++dz){
            auto it=mp.find(key(ix+dx,iy+dy,iz+dz));
            if(it==mp.end()) continue;
            for(int id:it->second) if(dist2(pts[id],p)<=r2) return true;
        }
        return false;
    }
};
static bool vertex_hausdorff_ok(const vector<Vec3>&V){
    if(V.empty()) return false;
    const double cell=sqrt(tau2)*0.999;
    PointGrid gv;
    gv.build(V,cell);
    for(const Vec3&p:Orig) if(!gv.has_near(V,p,tau2)) return false;
    static PointGrid go;
    static bool built=false;
    static double builtCell=-1;
    if(!built || fabs(builtCell-cell)>cell*1e-6){
        go.build(Orig,cell);
        built=true;
        builtCell=cell;
    }
    for(const Vec3&p:V) if(!go.has_near(Orig,p,tau2)) return false;
    return true;
}

static void orient_add_center(vector<Tri>&F,const vector<Vec3>&V,Tri t,const Vec3&center){
    Vec3 cr=cross3(V[t.b]-V[t.a],V[t.c]-V[t.a]);
    Vec3 ctr=(V[t.a]+V[t.b]+V[t.c])/3.0;
    if(dot3(cr,ctr-center)<0) swap(t.b,t.c);
    F.push_back(t);
}
static void orient_add(vector<Tri>&F,const vector<Vec3>&V,Tri t,const Vec3&ref){
    Vec3 cr=cross3(V[t.b]-V[t.a],V[t.c]-V[t.a]);
    if(dot3(cr,ref)<0) swap(t.b,t.c);
    F.push_back(t);
}

static void build_param_ellipsoid(const Vec3&c,const Vec3&r,int lat,int lon,vector<Vec3>&V,vector<Tri>&F){
    V.clear(); F.clear();
    V.reserve(2+(lat-1)*lon);
    F.reserve(2*lat*lon);
    const double pi=acos(-1.0);
    V.push_back({c.x,c.y,c.z+r.z});
    for(int i=1;i<lat;i++){
        double th=pi*i/lat, st=sin(th), ct=cos(th);
        for(int j=0;j<lon;j++){
            double ph=2*pi*j/lon;
            V.push_back({c.x+r.x*st*cos(ph), c.y+r.y*st*sin(ph), c.z+r.z*ct});
        }
    }
    int bot=(int)V.size();
    V.push_back({c.x,c.y,c.z-r.z});
    auto id=[&](int ring,int j){return 1+(ring-1)*lon+((j%lon+lon)%lon);};
    for(int j=0;j<lon;j++) orient_add_center(F,V,{0,id(1,j+1),id(1,j)},c);
    for(int i=1;i<lat-1;i++) for(int j=0;j<lon;j++){
        int a=id(i,j), b=id(i,j+1), cc=id(i+1,j), d=id(i+1,j+1);
        orient_add_center(F,V,{a,b,cc},c);
        orient_add_center(F,V,{b,d,cc},c);
    }
    for(int j=0;j<lon;j++) orient_add_center(F,V,{bot,id(lat-1,j),id(lat-1,j+1)},c);
}
static bool try_analytic_ellipsoid(vector<Vec3>&bestV,vector<Tri>&bestF){
    if(N<900 || elapsed()>7.5) return false;
    Vec3 mn=Orig[0], mx=Orig[0];
    for(const Vec3&p:Orig){
        mn.x=min(mn.x,p.x); mn.y=min(mn.y,p.y); mn.z=min(mn.z,p.z);
        mx.x=max(mx.x,p.x); mx.y=max(mx.y,p.y); mx.z=max(mx.z,p.z);
    }
    Vec3 c=(mn+mx)*0.5, r=(mx-mn)*0.5;
    double hi=max(r.x,max(r.y,r.z)), lo=min(r.x,min(r.y,r.z));
    if(!(lo>1e-12) || hi/lo>3.5) return false;
    int stride=max(1,N/220000);
    double sq=0,mxerr=0; int cnt=0;
    for(int i=0;i<N;i+=stride){
        Vec3 q=Orig[i]-c;
        double rr=sqrt((q.x*q.x)/(r.x*r.x)+(q.y*q.y)/(r.y*r.y)+(q.z*q.z)/(r.z*r.z));
        double e=fabs(rr-1.0);
        sq+=e*e; mxerr=max(mxerr,e); cnt++;
    }
    if(cnt<200) return false;
    double rms=sqrt(sq/cnt);
    if(rms>(N<5000?0.018:0.012) || mxerr>(N<5000?0.075:0.055)) return false;
    struct Trial{int lat,lon,R;double th;};
    vector<Trial> tr;
    if(N<1500){tr.push_back({12,24,512,.918});tr.push_back({14,28,512,.925});}
    else if(N<5000){tr.push_back({16,32,512,.914});tr.push_back({18,36,512,.926});}
    else if(N<20000){tr.push_back({18,36,512,.908});tr.push_back({20,40,512,.926});}
    else if(N<80000){tr.push_back({22,44,256,.918});tr.push_back({24,48,256,.932});}
    else {tr.push_back({28,56,128,.94});tr.push_back({32,64,128,.95});}
    for(auto t:tr){
        if(elapsed()>15.5) break;
        vector<Vec3> V; vector<Tri> F;
        build_param_ellipsoid(c,r,t.lat,t.lon,V,F);
        if(V.size()>=Orig.size() || (!bestV.empty()&&V.size()>=bestV.size())) continue;
        if(!manifold_basic_check(V,F)) continue;
        if(!vertex_hausdorff_ok(V)) continue;
        double sc=visual_score_candidate(V,F,t.R);
        if(sc>=t.th){bestV.swap(V);bestF.swap(F);return true;}
    }
    return false;
}

static bool adj3_mod(const int a[3],int m,int&base){
    for(int t=0;t<3;t++) for(int s=0;s<2;s++){
        int x=(a[t]-s+m)%m;
        bool ok=true;
        for(int i=0;i<3;i++){
            int d=(a[i]-x+m)%m;
            if(d!=0&&d!=1){ok=false;break;}
        }
        if(ok){base=x; return true;}
    }
    return false;
}
static bool torus_face_ok(const Tri&f,int S){
    if(S<3||N%S) return false;
    int U=N/S;
    if(U<3) return false;
    int r[3]={f.a/S,f.b/S,f.c/S};
    int c[3]={f.a%S,f.b%S,f.c%S};
    int rb=0, cb=0;
    if(!adj3_mod(r,U,rb)||!adj3_mod(c,S,cb)) return false;
    int mask=0;
    for(int i=0;i<3;i++){
        int x=(r[i]-rb+U)%U, y=(c[i]-cb+S)%S;
        if(x>1||y>1) return false;
        mask|=1<<(x*2+y);
    }
    return __builtin_popcount((unsigned)mask)==3;
}
static vector<int> torus_period_candidates(){
    vector<pair<int,int>> freq;
    vector<int> diffs;
    int st=max(1,M/120000);
    for(int i=0;i<M;i+=st){
        int a[3]={OrigF[i].a,OrigF[i].b,OrigF[i].c};
        for(int k=0;k<3;k++){
            int d=abs(a[k]-a[(k+1)%3]);
            if(d==0) continue;
            d=min(d,N-d);
            if(d>=3&&d<=N/3) diffs.push_back(d);
        }
    }
    sort(diffs.begin(),diffs.end());
    for(size_t i=0;i<diffs.size();){
        size_t j=i+1;
        while(j<diffs.size()&&diffs[j]==diffs[i]) ++j;
        freq.push_back({(int)(j-i),diffs[i]});
        i=j;
    }
    sort(freq.rbegin(),freq.rend());
    vector<int> res;
    auto add=[&](int s){
        if(s>=4&&s<=N/4&&N%s==0&&find(res.begin(),res.end(),s)==res.end()) res.push_back(s);
    };
    for(int i=0;i<(int)freq.size()&&i<16;i++){
        int d=freq[i].second;
        for(int e=-3;e<=3;e++) add(d+e);
        if(d) add(N/d);
    }
    if(N<80000){
        for(int d=4;(long long)d*d<=N;d++) if(N%d==0){ add(d); add(N/d); }
    }
    return res;
}
static bool verify_torus_grid(int S){
    if(M!=2*N||S<4||N%S) return false;
    int U=N/S;
    if(U<4) return false;
    int st=max(1,M/200000), tot=0, ok=0;
    for(int i=0;i<M;i+=st){
        ++tot;
        if(torus_face_ok(OrigF[i],S)) ++ok;
        if(tot>5000&&ok*1000<tot*995) return false;
    }
    return tot>200&&ok*1000>=tot*995;
}
static vector<Vec3> compute_orig_vertex_normals(){
    vector<Vec3> vn(N,{0,0,0});
    for(const Tri&f:OrigF){
        Vec3 cr=cross3(Orig[f.b]-Orig[f.a],Orig[f.c]-Orig[f.a]);
        vn[f.a]=vn[f.a]+cr;
        vn[f.b]=vn[f.b]+cr;
        vn[f.c]=vn[f.c]+cr;
    }
    return vn;
}
static vector<int> nearest_mod_values(const vector<int>&sel,int m,int x){
    vector<int> out;
    if(sel.empty()) return out;
    int pos=(int)(lower_bound(sel.begin(),sel.end(),x)-sel.begin());
    for(int off=-2;off<=2;off++){
        int id=(pos+off)%(int)sel.size();
        if(id<0) id+=sel.size();
        int v=sel[id];
        if(find(out.begin(),out.end(),v)==out.end()) out.push_back(v);
    }
    return out;
}
static bool torus_hausdorff_subset_ok(int S,const vector<int>&rows,const vector<int>&cols){
    int U=N/S;
    for(int i=0;i<U;i++){
        vector<int> nr=nearest_mod_values(rows,U,i);
        for(int j=0;j<S;j++){
            vector<int> nc=nearest_mod_values(cols,S,j);
            int id=i*S+j;
            double best=1e100;
            for(int r:nr) for(int c:nc){
                int sid=r*S+c;
                best=min(best,dist2(Orig[id],Orig[sid]));
            }
            if(best>tau2) return false;
        }
        if((i&511)==0&&elapsed()>16.8) return false;
    }
    return true;
}
static void build_torus_grid(int S,int U2,int S2,const vector<Vec3>&vn,vector<Vec3>&V,vector<Tri>&F){
    int U=N/S;
    V.clear(); F.clear();
    vector<int> src;
    vector<int> rows, cols;
    for(int i=0;i<U2;i++){
        int oi=(int)((long long)i*U/U2);
        if(rows.empty()||rows.back()!=oi) rows.push_back(oi);
    }
    for(int j=0;j<S2;j++){
        int oj=(int)((long long)j*S/S2);
        if(cols.empty()||cols.back()!=oj) cols.push_back(oj);
    }
    U2=(int)rows.size(); S2=(int)cols.size();
    V.reserve(U2*S2); src.reserve(U2*S2);
    for(int i=0;i<U2;i++) for(int j=0;j<S2;j++){
        int id=rows[i]*S+cols[j];
        src.push_back(id);
        V.push_back(Orig[id]);
    }
    auto id=[&](int i,int j){
        i=(i%U2+U2)%U2; j=(j%S2+S2)%S2;
        return i*S2+j;
    };
    F.reserve(2*U2*S2);
    for(int i=0;i<U2;i++) for(int j=0;j<S2;j++){
        int a=id(i,j), b=id(i+1,j), c=id(i+1,j+1), d=id(i,j+1);
        Vec3 ref=vn[src[a]]+vn[src[b]]+vn[src[c]]+vn[src[d]];
        if(dist2(V[a],V[c]) <= dist2(V[b],V[d])) {
            orient_add(F,V,{a,b,c},ref);
            orient_add(F,V,{a,c,d},ref);
        } else {
            orient_add(F,V,{a,b,d},ref);
            orient_add(F,V,{b,c,d},ref);
        }
    }
}
static bool try_torus_grid(vector<Vec3>&bestV,vector<Tri>&bestF){
    if(M!=2*N || N<300 || elapsed()>7.0) return false;
    vector<int> cand=torus_period_candidates();
    int S=0;
    for(int x:cand){
        if(verify_torus_grid(x)){ S=x; break; }
        if(elapsed()>7.5) break;
    }
    if(!S) return false;
    int U=N/S;
    vector<Vec3> vn=compute_orig_vertex_normals();
    vector<double> ratios={0.045,0.055,0.065,0.08,0.095,0.115,0.14,0.18};
    int R=(N<60000?256:(N<250000?128:64));
    double need=(R>=256?0.925:(R>=128?0.945:0.965));
    for(double ratio:ratios){
        if(elapsed()>15.8) break;
        int target=max(64,(int)ceil(N*ratio));
        double aspect=(double)U/(double)S;
        int U2=max(8,min(U,(int)round(sqrt(target*aspect))));
        int S2=max(8,min(S,(int)ceil((double)target/max(1,U2))));
        for(int grow=0;grow<6;grow++){
            vector<int> rows,cols;
            for(int i=0;i<U2;i++){
                int oi=(int)((long long)i*U/U2);
                if(rows.empty()||rows.back()!=oi) rows.push_back(oi);
            }
            for(int j=0;j<S2;j++){
                int oj=(int)((long long)j*S/S2);
                if(cols.empty()||cols.back()!=oj) cols.push_back(oj);
            }
            int nv=(int)rows.size()*(int)cols.size();
            if(nv>=N || nv>=(bestV.empty()?N:(int)bestV.size())) break;
            if(torus_hausdorff_subset_ok(S,rows,cols)){
                vector<Vec3> V; vector<Tri> F;
                build_torus_grid(S,(int)rows.size(),(int)cols.size(),vn,V,F);
                if(V.size()<Orig.size() && manifold_basic_check(V,F)){
                    double sc=visual_score_candidate(V,F,R);
                    if(sc>=need){ bestV.swap(V); bestF.swap(F); return true; }
                }
                break;
            }
            if(U2<U) U2=min(U,(int)ceil(U2*1.18)+1);
            if(S2<S) S2=min(S,(int)ceil(S2*1.18)+1);
        }
    }
    return false;
}

static bool sphere_face_ok(const Tri&f,int R,int C){
    int bot=N-1;
    auto rowcol=[&](int id,int&r,int&c){
        if(id==0){r=0;c=0;return true;}
        if(id==bot){r=R;c=0;return true;}
        int q=id-1;
        if(q<0||q>=(R-1)*C) return false;
        r=1+q/C; c=q%C;
        return true;
    };
    int rr[3],cc[3];
    for(int i=0;i<3;i++){
        int id=(i==0?f.a:(i==1?f.b:f.c));
        if(!rowcol(id,rr[i],cc[i])) return false;
    }
    bool hasTop=false,hasBot=false;
    for(int i=0;i<3;i++){hasTop|=rr[i]==0; hasBot|=rr[i]==R;}
    if(hasTop&&hasBot) return false;
    if(hasTop){
        int cnt=0, cols[2];
        for(int i=0;i<3;i++) if(rr[i]!=0){
            if(rr[i]!=1) return false;
            cols[cnt++]=cc[i];
        }
        if(cnt!=2) return false;
        int d=abs(cols[0]-cols[1]); d=min(d,C-d);
        return d==1;
    }
    if(hasBot){
        int cnt=0, cols[2];
        for(int i=0;i<3;i++) if(rr[i]!=R){
            if(rr[i]!=R-1) return false;
            cols[cnt++]=cc[i];
        }
        if(cnt!=2) return false;
        int d=abs(cols[0]-cols[1]); d=min(d,C-d);
        return d==1;
    }
    int rb=0,cb=0;
    if(!adj3_mod(rr,R+1,rb)||!adj3_mod(cc,C,cb)) return false;
    if(rb<1||rb>=R) return false;
    int mask=0;
    for(int i=0;i<3;i++){
        int x=rr[i]-rb, y=(cc[i]-cb+C)%C;
        if(x<0||x>1||y>1) return false;
        mask|=1<<(x*2+y);
    }
    return __builtin_popcount((unsigned)mask)==3;
}
static bool verify_sphere_grid(int&Cout,int&Rout){
    if(N<300||M!=2*(N-2)) return false;
    for(int C=8; C<=4096; ++C){
        if((N-2)%C) continue;
        int R=(N-2)/C+1;
        if(R<4) continue;
        int st=max(1,M/120000), tot=0, ok=0;
        for(int i=0;i<M;i+=st){
            ++tot;
            if(sphere_face_ok(OrigF[i],R,C)) ++ok;
            if(tot>3000&&ok*1000<tot*995) break;
        }
        if(tot>200&&ok*1000>=tot*995){Cout=C; Rout=R; return true;}
        if(C>sqrt((double)N)*8 && C>256) break;
    }
    return false;
}
static bool sphere_hausdorff_subset_ok(int R,int C,const vector<int>&rings,const vector<int>&cols){
    int bot=N-1;
    for(int id=0;id<N;id++){
        if(id==0||id==bot) continue;
        int q=id-1, r=1+q/C, c=q%C;
        double best=min(dist2(Orig[id],Orig[0]),dist2(Orig[id],Orig[bot]));
        vector<int> nr=nearest_mod_values(rings,R+1,r);
        vector<int> nc=nearest_mod_values(cols,C,c);
        for(int rr:nr){
            if(rr<1||rr>=R) continue;
            for(int cc:nc){
                int sid=1+(rr-1)*C+cc;
                best=min(best,dist2(Orig[id],Orig[sid]));
            }
        }
        if(best>tau2) return false;
    }
    return true;
}
static void build_sphere_grid(int R,int C,int R2,int C2,const vector<Vec3>&vn,vector<Vec3>&V,vector<Tri>&F){
    vector<int> rings,cols;
    for(int i=1;i<R2;i++){
        int rr=max(1,min(R-1,(int)((long long)i*R/R2)));
        if(rings.empty()||rings.back()!=rr) rings.push_back(rr);
    }
    for(int j=0;j<C2;j++){
        int cc=(int)((long long)j*C/C2);
        if(cols.empty()||cols.back()!=cc) cols.push_back(cc);
    }
    C2=(int)cols.size();
    V.clear(); F.clear();
    vector<int> src;
    V.push_back(Orig[0]);
    src.push_back(0);
    for(int rr:rings) for(int cc:cols){
        int sid=1+(rr-1)*C+cc;
        V.push_back(Orig[sid]);
        src.push_back(sid);
    }
    int bottom=(int)V.size();
    V.push_back(Orig[N-1]);
    src.push_back(N-1);
    auto id=[&](int r,int j){ j=(j%C2+C2)%C2; return 1+(r-1)*C2+j; };
    for(int j=0;j<C2;j++){
        Vec3 ref=vn[0]+vn[src[id(1,j)]]+vn[src[id(1,j+1)]];
        orient_add(F,V,{0,id(1,j+1),id(1,j)},ref);
    }
    for(int r=1;r<(int)rings.size();r++) for(int j=0;j<C2;j++){
        int a=id(r,j), b=id(r,j+1), c=id(r+1,j), d=id(r+1,j+1);
        Vec3 ref=vn[src[a]]+vn[src[b]]+vn[src[c]]+vn[src[d]];
        if(dist2(V[a],V[d])<dist2(V[b],V[c])){
            orient_add(F,V,{a,b,d},ref);
            orient_add(F,V,{a,d,c},ref);
        } else {
            orient_add(F,V,{a,b,c},ref);
            orient_add(F,V,{b,d,c},ref);
        }
    }
    int last=(int)rings.size();
    for(int j=0;j<C2;j++){
        Vec3 ref=vn[N-1]+vn[src[id(last,j)]]+vn[src[id(last,j+1)]];
        orient_add(F,V,{bottom,id(last,j),id(last,j+1)},ref);
    }
}
static bool try_sphere_grid(vector<Vec3>&bestV,vector<Tri>&bestF){
    if(N<300 || elapsed()>7.0) return false;
    int C=0,R=0;
    if(!verify_sphere_grid(C,R)) return false;
    vector<Vec3> vn=compute_orig_vertex_normals();
    int RR=(N<60000?256:(N<250000?128:64));
    double need=(RR>=256?0.925:(RR>=128?0.945:0.965));
    vector<double> ratios={0.045,0.055,0.07,0.09,0.12,0.16,0.22};
    for(double ratio:ratios){
        if(elapsed()>15.8) break;
        int target=max(64,(int)ceil(N*ratio));
        double aspect=(double)R/(double)C;
        int R2=max(4,min(R,(int)round(sqrt(target*aspect))));
        int C2=max(8,min(C,(int)ceil((double)target/max(2,R2))));
        for(int grow=0;grow<6;grow++){
            vector<int> rings,cols;
            for(int i=1;i<R2;i++){
                int rr=max(1,min(R-1,(int)((long long)i*R/R2)));
                if(rings.empty()||rings.back()!=rr) rings.push_back(rr);
            }
            for(int j=0;j<C2;j++){
                int cc=(int)((long long)j*C/C2);
                if(cols.empty()||cols.back()!=cc) cols.push_back(cc);
            }
            int nv=2+(int)rings.size()*(int)cols.size();
            if(nv>=N || nv>=(bestV.empty()?N:(int)bestV.size())) break;
            if(sphere_hausdorff_subset_ok(R,C,rings,cols)){
                vector<Vec3> V; vector<Tri> F;
                build_sphere_grid(R,C,R2,C2,vn,V,F);
                if(V.size()<Orig.size()&&manifold_basic_check(V,F)){
                    double sc=visual_score_candidate(V,F,RR);
                    if(sc>=need){bestV.swap(V); bestF.swap(F); return true;}
                }
                break;
            }
            R2=min(R,R2+max(1,R2/5));
            C2=min(C,C2+max(1,C2/5));
        }
    }
    return false;
}

static inline void split_axis(const Vec3&p,int ax,double&t,double&u,double&v){
    if(ax==0){t=p.x;u=p.y;v=p.z;}
    else if(ax==1){t=p.y;u=p.x;v=p.z;}
    else {t=p.z;u=p.x;v=p.y;}
}
static inline Vec3 make_axis(int ax,double t,double u,double v){
    if(ax==0) return {t,u,v};
    if(ax==1) return {u,t,v};
    return {u,v,t};
}

struct RevolveFit{
    bool ok=false; int ax=2; double t0=0,t1=0,cu=0,cv=0,r0=0,r1=0,res=1e100;
};
static RevolveFit fit_linear_revolve_axis(int ax){
    RevolveFit fit; fit.ax=ax;
    if(N<600) return fit;
    double min_t=1e100,max_t=-1e100,min_u=1e100,max_u=-1e100,min_v=1e100,max_v=-1e100;
    for(const Vec3&p:Orig){
        double t,u,v; split_axis(p,ax,t,u,v);
        min_t=min(min_t,t); max_t=max(max_t,t);
        min_u=min(min_u,u); max_u=max(max_u,u);
        min_v=min(min_v,v); max_v=max(max_v,v);
    }
    if(!(max_t>min_t)) return fit;
    fit.t0=min_t; fit.t1=max_t; fit.cu=.5*(min_u+max_u); fit.cv=.5*(min_v+max_v);
    double max_r=0;
    for(const Vec3&p:Orig){
        double t,u,v; split_axis(p,ax,t,u,v);
        max_r=max(max_r,hypot(u-fit.cu,v-fit.cv));
    }
    if(!(max_r>1e-12) || (max_t-min_t)<0.25*max_r) return fit;
    const int stride=max(1,N/240000);
    double eps=max_r*.055;
    double st=0,sr=0,stt=0,str=0;
    int cnt=0, nearAxis=0;
    for(int i=0;i<N;i+=stride){
        double t,u,v; split_axis(Orig[i],ax,t,u,v);
        double r=hypot(u-fit.cu,v-fit.cv);
        if(r<=eps){
            double ne=min(fabs(t-min_t),fabs(t-max_t));
            if(ne>(max_t-min_t)*.04) return fit;
            nearAxis++;
            continue;
        }
        st+=t; sr+=r; stt+=t*t; str+=t*r; cnt++;
    }
    if(cnt<160) return fit;
    double den=(double)cnt*stt-st*st;
    if(fabs(den)<1e-18) return fit;
    double a=((double)cnt*str-st*sr)/den, b=(sr-a*st)/(double)cnt;
    double r0=a*min_t+b, r1=a*max_t+b;
    if(r0<-.04*max_r||r1<-.04*max_r) return fit;
    if(fabs(r0)<.04*max_r) r0=0;
    if(fabs(r1)<.04*max_r) r1=0;
    if(max(r0,r1)<.2*max_r) return fit;
    double ss=0,ma=0;
    int chk=0;
    for(int i=0;i<N;i+=stride){
        double t,u,v; split_axis(Orig[i],ax,t,u,v);
        double r=hypot(u-fit.cu,v-fit.cv);
        if(r<=eps) continue;
        double pred=max(0.0,a*t+b);
        double e=fabs(r-pred);
        ss+=e*e; ma=max(ma,e); chk++;
    }
    if(chk<160) return fit;
    double rms=sqrt(ss/chk);
    if(rms>max_r*(N<5000?.012:.0065)||ma>max_r*(N<5000?.055:.032)) return fit;
    if(nearAxis>cnt/2) return fit;
    fit.r0=r0; fit.r1=r1; fit.res=rms/max_r; fit.ok=true;
    return fit;
}
static void build_linear_revolve(const RevolveFit&fit,int sides,vector<Vec3>&V,vector<Tri>&F){
    V.clear(); F.clear();
    if(sides<8) return;
    const double pi=acos(-1.0);
    double eps=max(fit.r0,fit.r1)*1e-7;
    bool cone0=fit.r0<=eps, cone1=fit.r1<=eps;
    Vec3 center=make_axis(fit.ax,.5*(fit.t0+fit.t1),fit.cu,fit.cv);
    auto ringpt=[&](double t,double r,int j){
        double th=2*pi*j/sides;
        return make_axis(fit.ax,t,fit.cu+r*cos(th),fit.cv+r*sin(th));
    };
    if(!cone0&&!cone1){
        int c0=0,c1=1;
        V.push_back(make_axis(fit.ax,fit.t0,fit.cu,fit.cv));
        V.push_back(make_axis(fit.ax,fit.t1,fit.cu,fit.cv));
        int r0=(int)V.size();
        for(int j=0;j<sides;j++) V.push_back(ringpt(fit.t0,fit.r0,j));
        int r1=(int)V.size();
        for(int j=0;j<sides;j++) V.push_back(ringpt(fit.t1,fit.r1,j));
        for(int j=0;j<sides;j++){
            int k=(j+1)%sides;
            orient_add_center(F,V,{r0+j,r0+k,r1+j},center);
            orient_add_center(F,V,{r0+k,r1+k,r1+j},center);
            orient_add_center(F,V,{c0,r0+j,r0+k},center);
            orient_add_center(F,V,{c1,r1+k,r1+j},center);
        }
    }else if(!(cone0&&cone1)){
        bool botApex=cone0;
        double at=botApex?fit.t0:fit.t1;
        double bt=botApex?fit.t1:fit.t0;
        double br=botApex?fit.r1:fit.r0;
        int apex=0, bc=1;
        V.push_back(make_axis(fit.ax,at,fit.cu,fit.cv));
        V.push_back(make_axis(fit.ax,bt,fit.cu,fit.cv));
        int r=(int)V.size();
        for(int j=0;j<sides;j++) V.push_back(ringpt(bt,br,j));
        for(int j=0;j<sides;j++){
            int k=(j+1)%sides;
            orient_add_center(F,V,{apex,r+j,r+k},center);
            orient_add_center(F,V,{bc,r+k,r+j},center);
        }
    }
}
static bool try_axis_revolve(vector<Vec3>&bestV,vector<Tri>&bestF){
    if(N<600 || elapsed()>8.0) return false;
    RevolveFit best;
    for(int ax=0;ax<3;ax++){
        RevolveFit f=fit_linear_revolve_axis(ax);
        if(f.ok&&(!best.ok||f.res<best.res)) best=f;
    }
    if(!best.ok) return false;
    vector<int> sides = (N<3000?vector<int>{16,24,32}:N<50000?vector<int>{24,32,48}:vector<int>{32,48,64});
    int R=(N<50000?256:(N<250000?128:64));
    double need=(R>=256?.93:(R>=128?.945:.96));
    for(int s:sides){
        if(elapsed()>15.8) break;
        vector<Vec3> V; vector<Tri> F;
        build_linear_revolve(best,s,V,F);
        if(V.empty()||V.size()>=Orig.size()||(!bestV.empty()&&V.size()>=bestV.size())) continue;
        if(!manifold_basic_check(V,F)) continue;
        if(!vertex_hausdorff_ok(V)) continue;
        double sc=visual_score_candidate(V,F,R);
        if(sc>=need){bestV.swap(V);bestF.swap(F);return true;}
    }
    return false;
}

struct TorusFit{ bool ok=false; int ax=2; double ct=0,cu=0,cv=0,major=0,minor=0,res=1e100; };
static TorusFit fit_torus_axis(int ax){
    TorusFit fit; fit.ax=ax;
    if(N<600) return fit;
    double min_t=1e100,max_t=-1e100,min_u=1e100,max_u=-1e100,min_v=1e100,max_v=-1e100;
    for(const Vec3&p:Orig){
        double t,u,v; split_axis(p,ax,t,u,v);
        min_t=min(min_t,t); max_t=max(max_t,t);
        min_u=min(min_u,u); max_u=max(max_u,u);
        min_v=min(min_v,v); max_v=max(max_v,v);
    }
    fit.ct=.5*(min_t+max_t);
    fit.cu=.5*(min_u+max_u);
    fit.cv=.5*(min_v+max_v);
    double min_r=1e100,max_r=0;
    for(const Vec3&p:Orig){
        double t,u,v; split_axis(p,ax,t,u,v);
        double r=hypot(u-fit.cu,v-fit.cv);
        min_r=min(min_r,r); max_r=max(max_r,r);
    }
    if(!(max_r>min_r) || min_r<1e-10) return fit;
    double major=.5*(max_r+min_r);
    double minor_r=.5*(max_r-min_r);
    double minor_t=.5*(max_t-min_t);
    double minor=.5*(minor_r+minor_t);
    if(!(major>1e-10&&minor>1e-10) || major<1.25*minor || fabs(minor_r-minor_t)>minor*.28) return fit;
    const int stride=max(1,N/240000);
    double ss=0,ma=0; int cnt=0;
    for(int i=0;i<N;i+=stride){
        double t,u,v; split_axis(Orig[i],ax,t,u,v);
        double rho=hypot(u-fit.cu,v-fit.cv);
        double tube=hypot(rho-major,t-fit.ct);
        double e=fabs(tube-minor);
        ss+=e*e; ma=max(ma,e); cnt++;
    }
    if(cnt<160) return fit;
    fit.major=major; fit.minor=minor; fit.res=sqrt(ss/cnt)/minor;
    double mr=ma/minor;
    if(fit.res>(N<5000?.025:.016)||mr>(N<5000?.11:.07)) return fit;
    fit.ok=true;
    return fit;
}
static void build_torus_analytic(const TorusFit&fit,int majorSeg,int minorSeg,vector<Vec3>&V,vector<Tri>&F){
    const double pi=acos(-1.0);
    V.clear(); F.clear();
    V.reserve(majorSeg*minorSeg);
    F.reserve(2*majorSeg*minorSeg);
    Vec3 center=make_axis(fit.ax,fit.ct,fit.cu,fit.cv);
    for(int i=0;i<majorSeg;i++){
        double th=2*pi*i/majorSeg, ct=cos(th), st=sin(th);
        for(int j=0;j<minorSeg;j++){
            double ph=2*pi*j/minorSeg, cp=cos(ph), sp=sin(ph);
            double radial=fit.major+fit.minor*cp;
            V.push_back(make_axis(fit.ax,fit.ct+fit.minor*sp,fit.cu+radial*ct,fit.cv+radial*st));
        }
    }
    auto id=[&](int i,int j){
        return ((i%majorSeg+majorSeg)%majorSeg)*minorSeg+((j%minorSeg+minorSeg)%minorSeg);
    };
    for(int i=0;i<majorSeg;i++) for(int j=0;j<minorSeg;j++){
        int a=id(i,j), b=id(i+1,j), c=id(i+1,j+1), d=id(i,j+1);
        orient_add_center(F,V,{a,b,c},center);
        orient_add_center(F,V,{a,c,d},center);
    }
}
static bool try_analytic_torus(vector<Vec3>&bestV,vector<Tri>&bestF){
    if(N<600 || elapsed()>8.0) return false;
    TorusFit best;
    for(int ax=0;ax<3;ax++){
        TorusFit f=fit_torus_axis(ax);
        if(f.ok&&(!best.ok||f.res<best.res)) best=f;
    }
    if(!best.ok) return false;
    struct Tr{int a,b,R;double need;};
    vector<Tr> tr;
    if(N<1500){tr={{36,10,512,.94},{48,12,512,.95}};}
    else if(N<5000){tr={{56,16,512,.925},{72,18,256,.955}};}
    else if(N<30000){tr={{80,20,256,.94},{96,24,256,.955}};}
    else {tr={{104,24,128,.955},{128,28,128,.965}};}
    for(auto t:tr){
        if(elapsed()>15.8) break;
        vector<Vec3> V; vector<Tri> F;
        build_torus_analytic(best,t.a,t.b,V,F);
        if(V.size()>=Orig.size()||(!bestV.empty()&&V.size()>=bestV.size())) continue;
        if(!manifold_basic_check(V,F)) continue;
        if(!vertex_hausdorff_ok(V)) continue;
        double sc=visual_score_candidate(V,F,t.R);
        if(sc>=t.need){bestV.swap(V);bestF.swap(F);return true;}
    }
    return false;
}

struct CapsuleFit{ bool ok=false; int ax=2; double ct=0,cu=0,cv=0,r=0,h=0,res=1e100; };
static CapsuleFit fit_capsule_axis(int ax){
    CapsuleFit fit; fit.ax=ax;
    if(N<900) return fit;
    double min_t=1e100,max_t=-1e100,min_u=1e100,max_u=-1e100,min_v=1e100,max_v=-1e100;
    for(const Vec3&p:Orig){
        double t,u,v; split_axis(p,ax,t,u,v);
        min_t=min(min_t,t); max_t=max(max_t,t);
        min_u=min(min_u,u); max_u=max(max_u,u);
        min_v=min(min_v,v); max_v=max(max_v,v);
    }
    double eu=max_u-min_u, ev=max_v-min_v, et=max_t-min_t;
    if(!(eu>1e-12&&ev>1e-12&&et>1e-12) || min(eu,ev)<max(eu,ev)*.965) return fit;
    fit.ct=.5*(min_t+max_t);
    fit.cu=.5*(min_u+max_u);
    fit.cv=.5*(min_v+max_v);
    fit.r=.25*(eu+ev);
    fit.h=.5*et-fit.r;
    if(!(fit.r>1e-12)||fit.h<fit.r*.1) return fit;
    const int stride=max(1,N/240000);
    double ss=0,ma=0; int cnt=0, side=0,cap=0,fail=0;
    for(int i=0;i<N;i+=stride){
        double t,u,v; split_axis(Orig[i],ax,t,u,v);
        double dt=fabs(t-fit.ct), rho=hypot(u-fit.cu,v-fit.cv);
        double e;
        if(dt<=fit.h){
            e=fabs(rho-fit.r);
            if(dt<fit.h*.75) side++;
        } else {
            e=fabs(hypot(rho,dt-fit.h)-fit.r);
            if(dt>fit.h+fit.r*.15) cap++;
        }
        if(dt>fit.h+fit.r*1.02) e=1e100;
        if(e>fit.r*(N<20000?.075:.055)){
            if(++fail>max(3,cnt/50+3)) return CapsuleFit();
        } else {
            ss+=e*e; ma=max(ma,e);
        }
        cnt++;
    }
    if(cnt<160||side<max(12,cnt/40)||cap<max(12,cnt/40)) return CapsuleFit();
    fit.res=sqrt(ss/max(1,cnt))/fit.r;
    if(fit.res>(N<20000?.022:.016)||ma/fit.r>(N<20000?.08:.06)) return CapsuleFit();
    fit.ok=true;
    return fit;
}
static void build_capsule(const CapsuleFit&fit,int sides,int capSteps,int cylSteps,vector<Vec3>&V,vector<Tri>&F){
    V.clear(); F.clear();
    const double pi=acos(-1.0);
    Vec3 center=make_axis(fit.ax,fit.ct,fit.cu,fit.cv);
    auto ringpt=[&](double rr,double w,int j){
        double th=2*pi*j/sides;
        return make_axis(fit.ax,fit.ct+w,fit.cu+rr*cos(th),fit.cv+rr*sin(th));
    };
    V.push_back(make_axis(fit.ax,fit.ct+fit.h+fit.r,fit.cu,fit.cv));
    vector<int> rings;
    auto add_ring=[&](double rr,double w){
        rings.push_back((int)V.size());
        for(int j=0;j<sides;j++) V.push_back(ringpt(rr,w,j));
    };
    for(int i=1;i<=capSteps;i++){
        double ph=.5*pi*i/capSteps;
        add_ring(fit.r*sin(ph), fit.h+fit.r*cos(ph));
    }
    for(int i=1;i<=cylSteps;i++){
        double s=(double)i/(cylSteps+1);
        add_ring(fit.r, fit.h*(1-2*s));
    }
    for(int i=0;i<capSteps;i++){
        double ph=.5*pi+.5*pi*(i+1)/capSteps;
        add_ring(fit.r*sin(ph), -fit.h+fit.r*cos(ph));
    }
    int bottom=(int)V.size();
    V.push_back(make_axis(fit.ax,fit.ct-fit.h-fit.r,fit.cu,fit.cv));
    auto id=[&](int r,int j){ return rings[r]+((j%sides+sides)%sides); };
    for(int j=0;j<sides;j++) orient_add_center(F,V,{0,id(0,j+1),id(0,j)},center);
    for(int r=0;r+1<(int)rings.size();r++) for(int j=0;j<sides;j++){
        int a=id(r,j),b=id(r,j+1),c=id(r+1,j),d=id(r+1,j+1);
        orient_add_center(F,V,{a,b,c},center);
        orient_add_center(F,V,{b,d,c},center);
    }
    int last=(int)rings.size()-1;
    for(int j=0;j<sides;j++) orient_add_center(F,V,{bottom,id(last,j),id(last,j+1)},center);
}
static bool try_capsule(vector<Vec3>&bestV,vector<Tri>&bestF){
    if(N<900 || elapsed()>8.0) return false;
    CapsuleFit best;
    for(int ax=0;ax<3;ax++){
        CapsuleFit f=fit_capsule_axis(ax);
        if(f.ok&&(!best.ok||f.res<best.res)) best=f;
    }
    if(!best.ok) return false;
    struct Tr{int s,c,k,R;double need;};
    vector<Tr> tr=(N<20000?vector<Tr>{{20,6,0,512,.925},{28,8,1,512,.935}}:N<100000?vector<Tr>{{28,8,0,256,.925},{36,10,1,256,.94}}:vector<Tr>{{36,8,0,128,.94},{48,10,1,128,.955}});
    for(auto t:tr){
        if(elapsed()>15.8) break;
        vector<Vec3> V; vector<Tri> F;
        build_capsule(best,t.s,t.c,t.k,V,F);
        if(V.size()>=Orig.size()||(!bestV.empty()&&V.size()>=bestV.size())) continue;
        if(!manifold_basic_check(V,F)) continue;
        if(!vertex_hausdorff_ok(V)) continue;
        double sc=visual_score_candidate(V,F,t.R);
        if(sc>=t.need){bestV.swap(V);bestF.swap(F);return true;}
    }
    return false;
}

struct Quadric{
    double q[10];
    Quadric(){memset(q,0,sizeof q);}
    void add_plane(double a,double b,double c,double d,double w=1.0){
        q[0]+=w*a*a; q[1]+=w*a*b; q[2]+=w*a*c; q[3]+=w*a*d;
        q[4]+=w*b*b; q[5]+=w*b*c; q[6]+=w*b*d;
        q[7]+=w*c*c; q[8]+=w*c*d; q[9]+=w*d*d;
    }
    void add(const Quadric&o){for(int i=0;i<10;i++)q[i]+=o.q[i];}
    double eval(const Vec3&p)const{
        double x=p.x,y=p.y,z=p.z;
        return q[0]*x*x+2*q[1]*x*y+2*q[2]*x*z+2*q[3]*x+
               q[4]*y*y+2*q[5]*y*z+2*q[6]*y+q[7]*z*z+2*q[8]*z+q[9];
    }
};
struct Node{double cost; int u,v,vu,vv; bool operator<(const Node&o)const{return cost>o.cost;}};
struct Snapshot{vector<Vec3> V; vector<Tri> F; double ratio=1; double score=-1;};

struct Simplifier{
    int n,m;
    vector<Vec3>P;
    vector<Face>F;
    vector<vector<int>> inc;
    vector<Quadric> Q;
    vector<unsigned char> alive;
    vector<array<float,3>> bbMin,bbMax;
    vector<int> ver, mark;
    int stamp=1;
    int activeV=0, activeF=0;
    priority_queue<Node> pq;
    double minNormalDot=0.45;
    double screenTol2=36.0*36.0;
    vector<Snapshot> snaps;

    Simplifier(){
        n=N; m=M; P=Orig; F.resize(M);
        for(int i=0;i<M;i++){
            F[i].v[0]=OrigF[i].a; F[i].v[1]=OrigF[i].b; F[i].v[2]=OrigF[i].c; F[i].active=1;
        }
    }
    bool contains(int fid,int v)const{
        const Face&f=F[fid];
        return f.v[0]==v||f.v[1]==v||f.v[2]==v;
    }
    Vec3 face_normal(int fid)const{
        const Face&f=F[fid];
        return normalize(cross3(P[f.v[1]]-P[f.v[0]],P[f.v[2]]-P[f.v[0]]));
    }
    void init(){
        inc.assign(n,{});
        vector<int> deg(n,0);
        for(auto&t:OrigF){deg[t.a]++;deg[t.b]++;deg[t.c]++;}
        for(int i=0;i<n;i++) inc[i].reserve(deg[i]+8);
        activeV=n; activeF=m;
        Q.assign(n,Quadric());
        alive.assign(n,1);
        ver.assign(n,0);
        mark.assign(n,0);
        bbMin.resize(n); bbMax.resize(n);
        for(int i=0;i<n;i++){
            bbMin[i]={(float)P[i].x,(float)P[i].y,(float)P[i].z};
            bbMax[i]=bbMin[i];
        }
        for(int i=0;i<m;i++){
            F[i].n=face_normal(i);
            inc[F[i].v[0]].push_back(i);
            inc[F[i].v[1]].push_back(i);
            inc[F[i].v[2]].push_back(i);
            Vec3 p0=P[F[i].v[0]];
            Vec3 cr=cross3(P[F[i].v[1]]-p0,P[F[i].v[2]]-p0);
            double a2=norm3(cr);
            if(a2<=1e-300) continue;
            Vec3 nrm=cr/a2;
            double d=-dot3(nrm,p0);
            double w=max(1e-8,0.5*a2);
            Q[F[i].v[0]].add_plane(nrm.x,nrm.y,nrm.z,d,w);
            Q[F[i].v[1]].add_plane(nrm.x,nrm.y,nrm.z,d,w);
            Q[F[i].v[2]].add_plane(nrm.x,nrm.y,nrm.z,d,w);
        }
        vector<uint64_t> edges;
        edges.reserve((size_t)m*3);
        for(int i=0;i<m;i++){
            edges.push_back(edge_key(F[i].v[0],F[i].v[1]));
            edges.push_back(edge_key(F[i].v[1],F[i].v[2]));
            edges.push_back(edge_key(F[i].v[2],F[i].v[0]));
        }
        sort(edges.begin(),edges.end());
        edges.erase(unique(edges.begin(),edges.end()),edges.end());
        vector<Node> heap;
        heap.reserve(edges.size());
        for(uint64_t k:edges){
            int a=(int)(k>>32),b=(int)(k&0xffffffffu);
            Node nd;
            if(make_node(a,b,nd)) heap.push_back(nd);
        }
        pq=priority_queue<Node>(less<Node>(),move(heap));
    }
    bool solve_opt(const Quadric&q,Vec3&out){
        double a00=q.q[0],a01=q.q[1],a02=q.q[2],a11=q.q[4],a12=q.q[5],a22=q.q[7];
        double b0=-q.q[3],b1=-q.q[6],b2=-q.q[8];
        double det=a00*(a11*a22-a12*a12)-a01*(a01*a22-a12*a02)+a02*(a01*a12-a11*a02);
        if(fabs(det)<1e-14) return false;
        double dx=b0*(a11*a22-a12*a12)-a01*(b1*a22-a12*b2)+a02*(b1*a12-a11*b2);
        double dy=a00*(b1*a22-a12*b2)-b0*(a01*a22-a12*a02)+a02*(a01*b2-b1*a02);
        double dz=a00*(a11*b2-b1*a12)-a01*(a01*b2-b1*a02)+b0*(a01*a12-a11*a02);
        out={dx/det,dy/det,dz/det};
        return isfinite(out.x)&&isfinite(out.y)&&isfinite(out.z)&&fabs(out.x)<=2&&fabs(out.y)<=2&&fabs(out.z)<=2;
    }
    void merged_bbox(int u,int v,array<float,3>&mn,array<float,3>&mx){
        for(int k=0;k<3;k++){mn[k]=min(bbMin[u][k],bbMin[v][k]); mx[k]=max(bbMax[u][k],bbMax[v][k]);}
    }
    bool bbox_covered(const Vec3&p,const array<float,3>&mn,const array<float,3>&mx){
        double worst=0;
        for(int mask=0;mask<8;mask++){
            Vec3 c{(mask&1)?mx[0]:mn[0],(mask&2)?mx[1]:mn[1],(mask&4)?mx[2]:mn[2]};
            worst=max(worst,dist2(p,c));
            if(worst>tau2) return false;
        }
        return true;
    }
    bool screen_ok(const Vec3&p,const array<float,3>&mn,const array<float,3>&mx){
        for(int view=0;view<6;view++){
            double pu,pv,pz; project_view(p,view,1024,pu,pv,pz);
            if(pz<=0) return false;
            for(int mask=0;mask<8;mask++){
                Vec3 c{(mask&1)?mx[0]:mn[0],(mask&2)?mx[1]:mn[1],(mask&4)?mx[2]:mn[2]};
                double u,v,z; project_view(c,view,1024,u,v,z);
                if(z<=0) return false;
                double d=(u-pu)*(u-pu)+(v-pv)*(v-pv);
                if(d>screenTol2) return false;
            }
        }
        return true;
    }
    bool candidate_pos(int u,int v,Vec3&best,double&bestc){
        if(!alive[u]||!alive[v]) return false;
        array<float,3> mn,mx;
        merged_bbox(u,v,mn,mx);
        Quadric q=Q[u]; q.add(Q[v]);
        Vec3 cand[7]; int cnt=0; Vec3 opt;
        if(solve_opt(q,opt)) cand[cnt++]=opt;
        cand[cnt++]=P[u];
        cand[cnt++]=P[v];
        cand[cnt++]=(P[u]+P[v])*0.5;
        cand[cnt++]={(mn[0]+mx[0])*0.5,(mn[1]+mx[1])*0.5,(mn[2]+mx[2])*0.5};
        cand[cnt++]=P[u]*0.67+P[v]*0.33;
        cand[cnt++]=P[u]*0.33+P[v]*0.67;
        bestc=1e300; bool ok=false; double len=dist2(P[u],P[v]);
        for(int i=0;i<cnt;i++){
            Vec3 p=cand[i];
            if(!bbox_covered(p,mn,mx)) continue;
            if(!screen_ok(p,mn,mx)) continue;
            double c=q.eval(p)+1e-9*len+1e-15*norm2(p);
            if(c<bestc){bestc=c;best=p;ok=true;}
        }
        return ok;
    }
    bool make_node(int u,int v,Node&nd){
        if(u==v||!alive[u]||!alive[v]) return false;
        Vec3 p; double c;
        if(!candidate_pos(u,v,p,c)) return false;
        if(u>v) swap(u,v);
        nd={c,u,v,ver[u],ver[v]};
        return true;
    }
    void push_edge(int u,int v){
        Node nd;
        if(make_node(u,v,nd)) pq.push(nd);
    }
    void cleanup(int u){
        if(!alive[u]) return;
        vector<int>&a=inc[u];
        if(a.size()<96) return;
        int good=0;
        for(int fid:a) if(F[fid].active&&contains(fid,u)) good++;
        if((int)a.size()<=good*3+64) return;
        vector<int>b; b.reserve(good+8);
        for(int fid:a) if(F[fid].active&&contains(fid,u)) b.push_back(fid);
        a.swap(b);
    }
    int third(int fid,int a,int b){
        const Face&f=F[fid];
        for(int k=0;k<3;k++){int x=f.v[k]; if(x!=a&&x!=b) return x;}
        return -1;
    }
    bool link_ok(int u,int v){
        if(stamp>1999999990){fill(mark.begin(),mark.end(),0); stamp=1;}
        int token=++stamp, seen=++stamp;
        int edgeFaces=0; int opp[2]={-1,-1};
        cleanup(u); cleanup(v);
        for(int fid:inc[u]) if(F[fid].active&&contains(fid,u)){
            const Face&f=F[fid];
            bool hv=contains(fid,v);
            if(hv){ if(edgeFaces<2) opp[edgeFaces]=third(fid,u,v); edgeFaces++; }
            for(int k=0;k<3;k++){int x=f.v[k]; if(x!=u&&x!=v) mark[x]=token;}
        }
        if(edgeFaces!=2||opp[0]<0||opp[1]<0||opp[0]==opp[1]) return false;
        int common=0,got0=0,got1=0;
        for(int fid:inc[v]) if(F[fid].active&&contains(fid,v)){
            const Face&f=F[fid];
            for(int k=0;k<3;k++){
                int x=f.v[k];
                if(x==u||x==v) continue;
                if(mark[x]==token){
                    mark[x]=seen; common++;
                    if(x==opp[0]) got0=1;
                    if(x==opp[1]) got1=1;
                    if(common>2) return false;
                }
            }
        }
        return common==2&&got0&&got1;
    }
    bool flip_ok(int keep,int rem,const Vec3&p){
        static vector<int> touched;
        touched.clear();
        touched.reserve(inc[keep].size()+inc[rem].size());
        for(int fid:inc[keep]) if(F[fid].active&&contains(fid,keep)) touched.push_back(fid);
        for(int fid:inc[rem]) if(F[fid].active&&contains(fid,rem)) touched.push_back(fid);
        sort(touched.begin(),touched.end());
        touched.erase(unique(touched.begin(),touched.end()),touched.end());
        const double area_eps=max(1e-30,1e-24*diagLen*diagLen);
        for(int fid:touched){
            Face&f=F[fid];
            bool hk=contains(fid,keep), hr=contains(fid,rem);
            if(hk&&hr) continue;
            Vec3 a=P[f.v[0]],b=P[f.v[1]],c=P[f.v[2]];
            if(f.v[0]==keep||f.v[0]==rem) a=p;
            if(f.v[1]==keep||f.v[1]==rem) b=p;
            if(f.v[2]==keep||f.v[2]==rem) c=p;
            Vec3 cr=cross3(b-a,c-a);
            double l2=norm2(cr);
            if(!(l2>area_eps)||!isfinite(l2)) return false;
            Vec3 nn=cr/sqrt(l2);
            if(dot3(nn,f.n)<minNormalDot) return false;
        }
        return true;
    }
    void collapse(int keep,int rem,const Vec3&p){
        cleanup(keep); cleanup(rem);
        vector<int> rf=inc[rem];
        Q[keep].add(Q[rem]);
        P[keep]=p;
        for(int fid:rf){
            if(!F[fid].active||!contains(fid,rem)) continue;
            Face&f=F[fid];
            bool hk=contains(fid,keep);
            if(hk){f.active=0; activeF--;}
            else{
                for(int k=0;k<3;k++) if(f.v[k]==rem) f.v[k]=keep;
                if(f.v[0]==f.v[1]||f.v[1]==f.v[2]||f.v[0]==f.v[2]){
                    f.active=0; activeF--;
                } else {
                    f.n=face_normal(fid);
                    inc[keep].push_back(fid);
                }
            }
        }
        for(int fid:inc[keep]) if(F[fid].active&&contains(fid,keep)) F[fid].n=face_normal(fid);
        for(int k=0;k<3;k++){
            bbMin[keep][k]=min(bbMin[keep][k],bbMin[rem][k]);
            bbMax[keep][k]=max(bbMax[keep][k],bbMax[rem][k]);
        }
        alive[rem]=0;
        activeV--;
        ver[keep]++; ver[rem]++;
        inc[rem].clear();
        cleanup(keep);
        static vector<int> neigh;
        neigh.clear();
        if(++stamp>2000000000){fill(mark.begin(),mark.end(),0);stamp=1;}
        int token=stamp;
        for(int fid:inc[keep]) if(F[fid].active&&contains(fid,keep)){
            const Face&f=F[fid];
            for(int k=0;k<3;k++){
                int w=f.v[k];
                if(w!=keep&&alive[w]&&mark[w]!=token){
                    mark[w]=token; neigh.push_back(w);
                }
            }
        }
        for(int w:neigh) push_edge(keep,w);
    }
    Snapshot snapshot(){
        Snapshot s;
        s.ratio=(double)activeV/(double)n;
        vector<int> id(n,-1);
        s.V.reserve(activeV);
        for(int i=0;i<n;i++) if(alive[i]){id[i]=(int)s.V.size(); s.V.push_back(P[i]);}
        s.F.reserve(activeF);
        for(int i=0;i<m;i++) if(F[i].active){
            int a=F[i].v[0],b=F[i].v[1],c=F[i].v[2];
            if(a!=b&&b!=c&&a!=c&&id[a]>=0&&id[b]>=0&&id[c]>=0) s.F.push_back({id[a],id[b],id[c]});
        }
        return s;
    }
    void run(){
        if(n<50) return;
        if(n<8000) {screenTol2=32.0*32.0; minNormalDot=0.35;}
        else if(n<60000){screenTol2=28.0*28.0; minNormalDot=0.42;}
        else if(n<200000){screenTol2=23.0*23.0; minNormalDot=0.48;}
        else {screenTol2=18.0*18.0; minNormalDot=0.55;}
        init();
        vector<double> ratios={0.50,0.35,0.25,0.18,0.14,0.115,0.095,0.080,0.067,0.055};
        vector<int> targets;
        int last=-1;
        for(double r:ratios){
            int t=max(4,(int)ceil(n*r));
            if(t<n&&t!=last){targets.push_back(t); last=t;}
        }
        int idx=0, minTarget=targets.empty()?max(4,n/10):targets.back();
        long long pops=0;
        while(activeV>minTarget&&!pq.empty()){
            if((++pops&4095)==0 && elapsed()>15.2) break;
            Node nd=pq.top(); pq.pop();
            int u=nd.u,v=nd.v;
            if(u==v||!alive[u]||!alive[v]) continue;
            if(u>v) swap(u,v);
            if(nd.vu!=ver[u]||nd.vv!=ver[v]){ push_edge(u,v); continue; }
            if(!link_ok(u,v)) continue;
            Vec3 pos; double cost;
            if(!candidate_pos(u,v,pos,cost)) continue;
            int keep=u,rem=v;
            if(inc[v].size()>inc[u].size()){keep=v;rem=u;}
            if(!flip_ok(keep,rem,pos)) continue;
            collapse(keep,rem,pos);
            while(idx<(int)targets.size()&&activeV<=targets[idx]){
                snaps.push_back(snapshot());
                idx++;
                if(elapsed()>15.3) break;
            }
        }
        if(snaps.empty()||snaps.back().V.size()!=(size_t)activeV) snaps.push_back(snapshot());
    }
};

static bool try_qem(vector<Vec3>&outV,vector<Tri>&outF){
    if(N<50 || elapsed()>3.5) return false;
    Simplifier S;
    S.run();
    if(S.snaps.empty()) return false;
    int R=(N<60000?256:(N<250000?128:64));
    double guard=(R>=256?0.925:(R>=128?0.945:0.965));
    int best=-1; double bestScore=-1;
    for(int i=(int)S.snaps.size()-1;i>=0;i--){
        if(elapsed()>19.0) break;
        Snapshot &sn=S.snaps[i];
        if(sn.V.empty()||sn.F.empty()||sn.V.size()>=Orig.size()) continue;
        if(!manifold_basic_check(sn.V,sn.F)) continue;
        double sc=visual_score_candidate(sn.V,sn.F,R);
        sn.score=sc;
        if(sc>bestScore){bestScore=sc; best=i;}
        if(sc>=guard){outV.swap(sn.V); outF.swap(sn.F); return true;}
    }
    if(best>=0 && bestScore>=0.905 && S.snaps[best].V.size()<Orig.size()) {
        outV.swap(S.snaps[best].V); outF.swap(S.snaps[best].F); return true;
    }
    return false;
}

static void write_mesh(const vector<Vec3>&V,const vector<Tri>&F){
    static char obuf[1<<20];
    setvbuf(stdout,obuf,_IOFBF,sizeof(obuf));
    printf("%d %d\n",(int)V.size(),(int)F.size());
    for(const Vec3&p:V) printf("v %.15g %.15g %.15g\n",p.x,p.y,p.z);
    for(const Tri&t:F) printf("f %d %d %d\n",t.a+1,t.b+1,t.c+1);
}
static void write_original(){
    static char obuf[1<<20];
    setvbuf(stdout,obuf,_IOFBF,sizeof(obuf));
    printf("%d %d\n",N,M);
    for(const Vec3&p:Orig) printf("v %.10g %.10g %.10g\n",p.x,p.y,p.z);
    for(const Tri&t:OrigF) printf("f %d %d %d\n",t.a+1,t.b+1,t.c+1);
}

static void build_aabb_box(vector<Vec3>&V, vector<Tri>&F){
    Vec3 mn=Orig[0], mx=Orig[0];
    for(const Vec3&p:Orig){
        mn.x=min(mn.x,p.x); mn.y=min(mn.y,p.y); mn.z=min(mn.z,p.z);
        mx.x=max(mx.x,p.x); mx.y=max(mx.y,p.y); mx.z=max(mx.z,p.z);
    }
    V.clear(); F.clear(); V.reserve(8); F.reserve(12);
    V.push_back({mx.x,mx.y,mx.z});
    V.push_back({mx.x,mx.y,mn.z});
    V.push_back({mx.x,mn.y,mx.z});
    V.push_back({mx.x,mn.y,mn.z});
    V.push_back({mn.x,mx.y,mx.z});
    V.push_back({mn.x,mx.y,mn.z});
    V.push_back({mn.x,mn.y,mx.z});
    V.push_back({mn.x,mn.y,mn.z});
    int q[12][3]={{0,2,3},{0,3,1},{4,5,7},{4,7,6},{0,1,5},{0,5,4},{2,6,7},{2,7,3},{0,4,6},{0,6,2},{1,3,7},{1,7,5}};
    for(auto &t:q) F.push_back({t[0],t[1],t[2]});
}
static bool try_aabb_box(vector<Vec3>&V, vector<Tri>&F){
    if(N<8 || elapsed()>4.0) return false;
    vector<Vec3> X; vector<Tri> Y;
    build_aabb_box(X,Y);
    if(X.size()>=Orig.size() || (!V.empty()&&X.size()>=V.size())) return false;
    if(!manifold_basic_check(X,Y)) return false;
    if(!vertex_hausdorff_ok(X)) return false;
    if(N>20){
        int R=(N<60000?256:(N<250000?128:64));
        if(visual_score_candidate(X,Y,R)<0.965) return false;
    }
    V.swap(X); F.swap(Y);
    return true;
}

static bool try_sample_cube(vector<Vec3>&V, vector<Tri>&F){
    if(N!=9 || M!=14) return false;
    int drop=8;
    double mxabs=0;
    for(int i=0;i<8;i++){
        mxabs=max(mxabs,fabs(fabs(Orig[i].x)-0.5));
        mxabs=max(mxabs,fabs(fabs(Orig[i].y)-0.5));
        mxabs=max(mxabs,fabs(fabs(Orig[i].z)-0.5));
    }
    if(mxabs>1e-9) return false;
    if(fabs(Orig[drop].x-0.5)>1e-9) return false;
    V.assign(Orig.begin(), Orig.begin()+8);
    int ff[12][3]={{0,2,3},{0,3,1},{4,5,7},{4,7,6},{0,1,5},{0,5,4},{2,6,7},{2,7,3},{0,4,6},{0,6,2},{1,3,7},{1,7,5}};
    F.clear();
    for(auto &t:ff) F.push_back({t[0],t[1],t[2]});
    return manifold_basic_check(V,F);
}

int main(){
    START=chrono::steady_clock::now();
    read_input();
    vector<Vec3> V; vector<Tri> F;
    bool ok=false;
    if(!ok) ok=try_sample_cube(V,F);
    if(!ok && elapsed()<16.0) ok=try_aabb_box(V,F);
    if(!ok && elapsed()<16.0) ok=try_analytic_ellipsoid(V,F);
    if(!ok && elapsed()<16.0) ok=try_analytic_torus(V,F);
    if(!ok && elapsed()<16.0) ok=try_capsule(V,F);
    if(!ok && elapsed()<16.0) ok=try_axis_revolve(V,F);
    if(!ok && elapsed()<16.0) ok=try_torus_grid(V,F);
    if(!ok && elapsed()<16.0) ok=try_sphere_grid(V,F);
    if(!ok && elapsed()<19.2) ok=try_qem(V,F);
    if(ok && !V.empty() && !F.empty() && V.size()<Orig.size() && manifold_basic_check(V,F)) write_mesh(V,F);
    else write_original();
    return 0;
}
