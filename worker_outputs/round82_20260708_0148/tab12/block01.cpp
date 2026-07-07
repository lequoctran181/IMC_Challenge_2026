C++

#include <bits/stdc++.h>
using namespace std;

struct Vec3{ double x=0,y=0,z=0; };
static inline Vec3 operator+(const Vec3&a,const Vec3&b){return {a.x+b.x,a.y+b.y,a.z+b.z};}
static inline Vec3 operator-(const Vec3&a,const Vec3&b){return {a.x-b.x,a.y-b.y,a.z-b.z};}
static inline Vec3 operator*(const Vec3&a,double s){return {a.x*s,a.y*s,a.z*s};}
static inline Vec3 operator/(const Vec3&a,double s){return {a.x/s,a.y/s,a.z/s};}
static inline double dot3(const Vec3&a,const Vec3&b){return a.x*b.x+a.y*b.y+a.z*b.z;}
static inline Vec3 cross3(const Vec3&a,const Vec3&b){return {a.y*b.z-a.z*b.y,a.z*b.x-a.x*b.z,a.x*b.y-a.y*b.x};}
static inline double norm2(const Vec3&a){return dot3(a,a);} 
static inline double norm3(const Vec3&a){return sqrt(max(0.0,norm2(a)));}
static inline bool norm_inplace(Vec3&v){double n=norm3(v); if(n<=1e-300) return false; v=v/n; return true;}
static inline double clampd(double x,double a,double b){return x<a?a:(x>b?b:x);} 
static inline uint64_t edge_key(int a,int b){ if(a>b) swap(a,b); return (uint64_t)(uint32_t)a<<32 | (uint32_t)b; }
static inline array<int,3> tri_key(int a,int b,int c){ array<int,3> t{a,b,c}; sort(t.begin(),t.end()); return t; }

struct Face{ int v[3]; };
struct Mesh{ vector<Vec3> V; vector<Face> F; string route; };
static int N=0,M=0; static vector<Vec3> P; static vector<Face> Fin; static Vec3 bbmin,bbmax; static double DIAG=1.0, EPSH=0.05; static chrono::steady_clock::time_point T0;
static inline double elapsed(){ return chrono::duration<double>(chrono::steady_clock::now()-T0).count(); }

struct FastInput{
    vector<char> b; char* p=nullptr;
    FastInput(){ char tmp[1<<16]; size_t n; while((n=fread(tmp,1,sizeof(tmp),stdin))>0) b.insert(b.end(),tmp,tmp+n); b.push_back('\0'); p=b.data(); }
    inline void skip(){ while(*p==' '||*p=='\n'||*p=='\r'||*p=='\t') ++p; }
    int next_int(){ skip(); int s=1; if(*p=='-'){s=-1;++p;} int x=0; while(*p>='0'&&*p<='9'){x=x*10+(*p-'0');++p;} return x*s; }
    double next_double(){ skip(); char* e=nullptr; double x=strtod(p,&e); p=e; return x; }
    char next_char(){ skip(); return *p++; }
    bool eof(){ skip(); return *p=='\0'; }
};
static bool read_mesh(){
    FastInput in; if(in.eof()) return false; N=in.next_int(); M=in.next_int(); if(N<=0||M<=0) return false;
    P.resize(N); Fin.resize(M); bbmin={1e100,1e100,1e100}; bbmax={-1e100,-1e100,-1e100};
    for(int i=0;i<N;i++){ char c=in.next_char(); if(c!='v'&&c!='V') --in.p; P[i].x=in.next_double(); P[i].y=in.next_double(); P[i].z=in.next_double(); bbmin.x=min(bbmin.x,P[i].x); bbmin.y=min(bbmin.y,P[i].y); bbmin.z=min(bbmin.z,P[i].z); bbmax.x=max(bbmax.x,P[i].x); bbmax.y=max(bbmax.y,P[i].y); bbmax.z=max(bbmax.z,P[i].z); }
    for(int i=0;i<M;i++){ char c=in.next_char(); if(c!='f'&&c!='F') --in.p; int a=in.next_int()-1,b=in.next_int()-1,c2=in.next_int()-1; Fin[i].v[0]=a; Fin[i].v[1]=b; Fin[i].v[2]=c2; }
    DIAG=norm3(bbmax-bbmin); if(!(DIAG>0)) DIAG=1.0; EPSH=0.05*DIAG*0.999;
    return true;
}

static double signed_vol6(const vector<Vec3>&V,const vector<Face>&F){ long double s=0; for(auto &f:F) s += (long double)dot3(V[f.v[0]], cross3(V[f.v[1]], V[f.v[2]])); return (double)s; }
static void orient_like_input(Mesh &m){ double a=signed_vol6(P,Fin), b=signed_vol6(m.V,m.F); if(a*b<0) for(auto &f:m.F) swap(f.v[1],f.v[2]); }
static void add_face_oriented(vector<Vec3>&V, vector<Face>&F, int a,int b,int c, const Vec3&ctr){ if(a==b||a==c||b==c) return; Vec3 cr=cross3(V[b]-V[a],V[c]-V[a]); if(norm2(cr)<=1e-30*DIAG*DIAG*DIAG*DIAG) return; Face f{{a,b,c}}; Vec3 fc=(V[a]+V[b]+V[c])/3.0; if(dot3(cr,fc-ctr)<0) swap(f.v[1],f.v[2]); F.push_back(f); }

struct GridHash3{
    double cell=1,inv=1; unordered_map<long long, vector<int>> mp;
    static uint64_t mix(uint64_t x){ x+=0x9e3779b97f4a7c15ULL; x=(x^(x>>30))*0xbf58476d1ce4e5b9ULL; x=(x^(x>>27))*0x94d049bb133111ebULL; return x^(x>>31); }
    static long long key(int x,int y,int z){ uint64_t h=mix((uint64_t)(int64_t)x); h^=(mix((uint64_t)(int64_t)y)<<1)|(mix((uint64_t)(int64_t)y)>>63); h^=(mix((uint64_t)(int64_t)z)<<7)|(mix((uint64_t)(int64_t)z)>>57); return (long long)h; }
    void coord(const Vec3&p,int&x,int&y,int&z)const{ x=(int)floor(p.x*inv); y=(int)floor(p.y*inv); z=(int)floor(p.z*inv); }
    void build(const vector<Vec3>&A,double c){ cell=max(c,1e-300); inv=1.0/cell; mp.clear(); mp.reserve(A.size()*2+16); for(int i=0;i<(int)A.size();++i){int x,y,z; coord(A[i],x,y,z); mp[key(x,y,z)].push_back(i);} }
    bool near_any(const vector<Vec3>&A,const Vec3&q,double eps2)const{ int x,y,z; coord(q,x,y,z); for(int dx=-1;dx<=1;dx++) for(int dy=-1;dy<=1;dy++) for(int dz=-1;dz<=1;dz++){ auto it=mp.find(key(x+dx,y+dy,z+dz)); if(it==mp.end()) continue; for(int id:it->second) if(norm2(A[id]-q)<=eps2) return true; } return false; }
};
static bool all_near(const vector<Vec3>&A,const vector<Vec3>&B,double eps){ if(A.empty()||B.empty()) return false; GridHash3 h; h.build(B,eps); double e2=eps*eps; for(const Vec3&q:A) if(!h.near_any(B,q,e2)) return false; return true; }
static bool closed_manifold_ok(const vector<Vec3>&V,const vector<Face>&F){
    if(V.empty()||F.empty()) return false; vector<uint64_t> e; e.reserve(F.size()*3); double eps2=max(1e-32,1e-24*DIAG*DIAG*DIAG*DIAG);
    for(auto &f:F){ for(int k=0;k<3;k++) if(f.v[k]<0||f.v[k]>=(int)V.size()) return false; if(f.v[0]==f.v[1]||f.v[0]==f.v[2]||f.v[1]==f.v[2]) return false; if(norm2(cross3(V[f.v[1]]-V[f.v[0]],V[f.v[2]]-V[f.v[0]]))<=eps2) return false; e.push_back(edge_key(f.v[0],f.v[1])); e.push_back(edge_key(f.v[1],f.v[2])); e.push_back(edge_key(f.v[2],f.v[0])); }
    sort(e.begin(),e.end()); for(size_t i=0;i<e.size();){ size_t j=i+1; while(j<e.size()&&e[j]==e[i]) ++j; if(j-i!=2) return false; i=j; } return true;
}
static bool accept_candidate(Mesh &m, bool do_guard=true){
    if(m.V.empty()||m.F.empty()||(int)m.V.size()>=N) return false; orient_like_input(m); if(!closed_manifold_ok(m.V,m.F)) return false; if(do_guard){ double eps=0.05*DIAG*1.001; if(!all_near(m.V,P,eps)) return false; if(!all_near(P,m.V,eps)) return false; } return true;
}

static Mesh original_mesh(){ Mesh m; m.route="ORIG"; m.V=P; m.F=Fin; return m; }
static Mesh make_bbox_box(){
    Mesh m; m.route="BBOX8"; Vec3 mn=bbmin,mx=bbmax; m.V={{mn.x,mn.y,mn.z},{mx.x,mn.y,mn.z},{mx.x,mx.y,mn.z},{mn.x,mx.y,mn.z},{mn.x,mn.y,mx.z},{mx.x,mn.y,mx.z},{mx.x,mx.y,mx.z},{mn.x,mx.y,mx.z}};
    int t[12][3]={{0,1,2},{0,2,3},{4,6,5},{4,7,6},{0,4,5},{0,5,1},{1,5,6},{1,6,2},{2,6,7},{2,7,3},{3,7,4},{3,4,0}}; for(auto &x:t) m.F.push_back({{x[0],x[1],x[2]}}); return m;
}

struct TopoInfo{ long long E=0,boundary=0,nonman=0,chi=0; int comps=1; double smooth15=0,sharp45=0; };
static TopoInfo topo_info(){
    TopoInfo ti; vector<uint64_t> edges; edges.reserve((size_t)M*3); vector<Vec3> fn(M); vector<int> p(N),sz(N,1); iota(p.begin(),p.end(),0); function<int(int)> fd=[&](int x){while(p[x]!=x){p[x]=p[p[x]];x=p[x];}return x;}; auto un=[&](int a,int b){a=fd(a);b=fd(b); if(a==b)return; if(sz[a]<sz[b])swap(a,b);p[b]=a;sz[a]+=sz[b];};
    for(int i=0;i<M;i++){ auto &f=Fin[i]; un(f.v[0],f.v[1]); un(f.v[1],f.v[2]); Vec3 n=cross3(P[f.v[1]]-P[f.v[0]],P[f.v[2]]-P[f.v[0]]); norm_inplace(n); fn[i]=n; edges.push_back(edge_key(f.v[0],f.v[1])); edges.push_back(edge_key(f.v[1],f.v[2])); edges.push_back(edge_key(f.v[2],f.v[0])); }
    unordered_set<int> cs; cs.reserve(N*2+1); for(int i=0;i<N;i++) cs.insert(fd(i)); ti.comps=cs.size(); sort(edges.begin(),edges.end()); long long adj=0,s15=0,sh45=0; double c15=cos(15*M_PI/180.0), c45=cos(45*M_PI/180.0);
    struct ER{uint64_t k;int fid;}; vector<ER> er; er.reserve((size_t)M*3); for(int i=0;i<M;i++){auto&f=Fin[i]; er.push_back({edge_key(f.v[0],f.v[1]),i});er.push_back({edge_key(f.v[1],f.v[2]),i});er.push_back({edge_key(f.v[2],f.v[0]),i});} sort(er.begin(),er.end(),[](const ER&a,const ER&b){return a.k<b.k;});
    for(size_t i=0;i<er.size();){ size_t j=i+1; while(j<er.size()&&er[j].k==er[i].k) ++j; ti.E++; if(j-i==1) ti.boundary++; else if(j-i!=2) ti.nonman++; else {adj++; double d=clampd(dot3(fn[er[i].fid],fn[er[i+1].fid]),-1,1); if(d>c15)s15++; if(d<c45)sh45++;} i=j; }
    ti.chi=(long long)N-ti.E+M; if(adj){ti.smooth15=(double)s15/adj; ti.sharp45=(double)sh45/adj;} return ti;
}

static inline void coords_axis(const Vec3&p,int ax,double&t,double&u,double&v){ if(ax==0){t=p.x;u=p.y;v=p.z;} else if(ax==1){t=p.y;u=p.x;v=p.z;} else {t=p.z;u=p.x;v=p.y;} }
static inline Vec3 point_axis(int ax,double t,double u,double v){ if(ax==0)return {t,u,v}; if(ax==1)return {u,t,v}; return {u,v,t}; }
struct Ring{ vector<int> ids; double t=0,cu=0,cv=0,r=0,rcv=0; bool cap=false; vector<pair<double,int>> ang; };
static bool circular_gap_ok(vector<pair<double,int>>&a){ if(a.size()<6) return false; sort(a.begin(),a.end()); double mg=0; for(int i=0;i<(int)a.size();++i){ double x=a[(i+1)%a.size()].first-a[i].first; if(i+1==(int)a.size()) x+=2*M_PI; mg=max(mg,x); } return mg < M_PI*0.85; }
static Mesh build_ring_candidate(const vector<Ring>&rings,int ax,const vector<int>&sel,int sides,const Vec3&center){
    Mesh m; m.route="RINGREV"; if(sel.size()<2||sides<6) return m; vector<int> base(sel.size(),-1), cnt(sel.size(),0); unordered_set<int> used; used.reserve(sel.size()*sides*2+8);
    for(int li=0;li<(int)sel.size();++li){ const Ring&r=rings[sel[li]]; if(r.cap){ int best=r.ids[0]; double bestd=1e300; for(int id:r.ids){double t,u,v; coords_axis(P[id],ax,t,u,v); double d=(u-r.cu)*(u-r.cu)+(v-r.cv)*(v-r.cv); if(d<bestd){bestd=d;best=id;}} if(used.count(best)) return {}; used.insert(best); base[li]=m.V.size(); cnt[li]=1; m.V.push_back(P[best]); }
        else{ base[li]=m.V.size(); cnt[li]=sides; vector<pair<double,int>> a=r.ang; sort(a.begin(),a.end()); vector<int> chosen; chosen.reserve(sides); for(int j=0;j<sides;j++){ double target=2*M_PI*j/sides; auto it=lower_bound(a.begin(),a.end(),make_pair(target,-1)); int cand=-1; double bd=1e100; for(int z=-1;z<=1;z++){ int k; if(it==a.end()) k=0; else k=it-a.begin(); k=(k+z+(int)a.size())%(int)a.size(); double da=fabs(a[k].first-target); da=min(da,2*M_PI-da); if(da<bd){bd=da;cand=a[k].second;} } if(cand<0||used.count(cand)) return {}; used.insert(cand); chosen.push_back(cand); }
            for(int id:chosen) m.V.push_back(P[id]); }
    }
    auto vid=[&](int li,int j){ if(cnt[li]==1) return base[li]; j=(j%sides+sides)%sides; return base[li]+j; };
    for(int li=0;li+1<(int)sel.size();++li){ if(cnt[li]==1&&cnt[li+1]==1) return {}; if(cnt[li]==1){ for(int j=0;j<sides;j++) add_face_oriented(m.V,m.F,vid(li,0),vid(li+1,j),vid(li+1,j+1),center); }
        else if(cnt[li+1]==1){ for(int j=0;j<sides;j++) add_face_oriented(m.V,m.F,vid(li+1,0),vid(li,j+1),vid(li,j),center); }
        else{ for(int j=0;j<sides;j++){ int a=vid(li,j), b=vid(li,j+1), c=vid(li+1,j), d=vid(li+1,j+1); add_face_oriented(m.V,m.F,a,b,c,center); add_face_oriented(m.V,m.F,b,d,c,center); } }
    }
    if(cnt[0]>1){ int li=0; for(int j=1;j+1<sides;j++) add_face_oriented(m.V,m.F,vid(li,0),vid(li,j),vid(li,j+1),center); }
    int last=sel.size()-1; if(cnt[last]>1){ for(int j=1;j+1<sides;j++) add_face_oriented(m.V,m.F,vid(last,0),vid(last,j+1),vid(last,j),center); }
    return m;
}
static Mesh try_ring_revolution(){
    if(N<200 || elapsed()>5.0) return {};
    for(int ax=0; ax<3; ++ax){
        double tol=max(DIAG*1e-8,1e-12);
        if(N>5000){
            int stride=max(1,N/120000); unordered_set<long long> uq; uq.reserve(N/stride*2+16);
            double scale=1.0/tol; for(int i=0;i<N;i+=stride){double t,u,v; coords_axis(P[i],ax,t,u,v); uq.insert((long long)llround(t*scale));}
            int samples=(N+stride-1)/stride; if((int)uq.size()>samples*7/10) continue;
        }
        vector<int> ord(N); iota(ord.begin(),ord.end(),0); sort(ord.begin(),ord.end(),[&](int a,int b){double ta,ua,va,tb,ub,vb; coords_axis(P[a],ax,ta,ua,va); coords_axis(P[b],ax,tb,ub,vb); return ta<tb;});
        vector<Ring> rings;
        for(int s=0;s<N;){ double t0,u,v; coords_axis(P[ord[s]],ax,t0,u,v); int e=s+1; while(e<N){double tt,uu,vv; coords_axis(P[ord[e]],ax,tt,uu,vv); if(fabs(tt-t0)>tol) break; ++e;} Ring r; r.ids.assign(ord.begin()+s,ord.begin()+e); rings.push_back(move(r)); s=e; }
        if((int)rings.size()<5 || (int)rings.size()>N/3) continue;
        long long largePts=0; int largeCnt=0; for(auto&r:rings) if((int)r.ids.size()>=8){largePts+=r.ids.size();largeCnt++;}
        if(largeCnt<3 || largePts < (long long)(0.65*N)) continue;
        double gcu=0,gcv=0,gw=0,tmin=1e100,tmax=-1e100;
        for(auto &r:rings){ double st=0,su=0,sv=0; for(int id:r.ids){double t,u,v; coords_axis(P[id],ax,t,u,v); st+=t;su+=u;sv+=v;} int c=r.ids.size(); r.t=st/c; r.cu=su/c; r.cv=sv/c; if(c>=8){gcu+=r.cu*c; gcv+=r.cv*c; gw+=c;} tmin=min(tmin,r.t);tmax=max(tmax,r.t); }
        if(gw<=0) continue; gcu/=gw; gcv/=gw;
        double centerDrift=0, driftW=0, rmax=0, cvSum=0, cvW=0; int goodR=0;
        for(auto &r:rings){ long double sr=0,sr2=0; vector<pair<double,int>> ang; ang.reserve(r.ids.size()); for(int id:r.ids){double t,u,v; coords_axis(P[id],ax,t,u,v); double du=u-gcu,dv=v-gcv; double rr=hypot(du,dv); sr+=rr; sr2+=(long double)rr*rr; double a=atan2(dv,du); if(a<0)a+=2*M_PI; ang.push_back({a,id});}
            int c=r.ids.size(); r.r=(double)(sr/c); double var=max(0.0,(double)(sr2/c)-r.r*r.r); r.rcv = r.r>1e-12?sqrt(var)/r.r:0; r.cap=(c<6 || r.r<0.018*DIAG); r.ang.swap(ang); rmax=max(rmax,r.r); if(c>=8&&!r.cap){ centerDrift += hypot(r.cu-gcu,r.cv-gcv)*c; driftW+=c; cvSum+=r.rcv*c; cvW+=c; if(r.rcv<0.08 && circular_gap_ok(r.ang)) goodR++; } }
        if(driftW>0 && centerDrift/driftW > 0.015*DIAG) continue; if(cvW>0 && cvSum/cvW > 0.045) continue; if(goodR < max(3,largeCnt*2/3)) continue;
        Vec3 ctr=point_axis(ax,0.5*(tmin+tmax),gcu,gcv);
        int sideReq=(int)ceil(2*M_PI*rmax/(0.78*0.05*DIAG)); sideReq=max(8,min(sideReq,160));
        int ringReq=(int)ceil((tmax-tmin)/(0.78*0.05*DIAG))+1; ringReq=max(3,min(ringReq,(int)rings.size()));
        vector<pair<int,int>> trials;
        auto addTrial=[&](int rr,int ss){ rr=max(3,min(rr,(int)rings.size())); ss=max(6,min(ss,sideReq)); if(rr*ss+4<N) trials.push_back({rr,ss}); };
        addTrial((ringReq*55+99)/100, (sideReq*65+99)/100); addTrial((ringReq*70+99)/100, (sideReq*75+99)/100); addTrial(ringReq, sideReq); addTrial((ringReq*120+99)/100, sideReq); addTrial(ringReq, min(160,(sideReq*120+99)/100));
        sort(trials.begin(),trials.end(),[](auto&a,auto&b){return a.first*a.second<b.first*b.second;}); trials.erase(unique(trials.begin(),trials.end()),trials.end());
        for(auto [R,S]:trials){ vector<int> sel; sel.reserve(R+4); sel.push_back(0); for(int k=1;k<R-1;k++){ int id=(int)llround((double)k*(rings.size()-1)/(R-1)); sel.push_back(id); } sel.push_back((int)rings.size()-1); sort(sel.begin(),sel.end()); sel.erase(unique(sel.begin(),sel.end()),sel.end()); if(sel.size()<2) continue; Mesh m=build_ring_candidate(rings,ax,sel,S,ctr); m.route="RINGREV"; if(accept_candidate(m,true)) return m; if(elapsed()>8.0) break; }
    }
    return {};
}

struct TorFit{bool ok=false; int ax=2; double ct=0,cu=0,cv=0,R=0,r=0,err=1e9;};
static TorFit fit_torus(int ax){ TorFit f; f.ax=ax; if(N<200) return f; double lt=1e100,ht=-1e100,lu=1e100,hu=-1e100,lv=1e100,hv=-1e100; for(auto&p:P){double t,u,v;coords_axis(p,ax,t,u,v);lt=min(lt,t);ht=max(ht,t);lu=min(lu,u);hu=max(hu,u);lv=min(lv,v);hv=max(hv,v);} f.ct=(lt+ht)/2; f.cu=(lu+hu)/2; f.cv=(lv+hv)/2; double minrho=1e100,maxrho=0; for(auto&p:P){double t,u,v;coords_axis(p,ax,t,u,v); double rho=hypot(u-f.cu,v-f.cv); minrho=min(minrho,rho); maxrho=max(maxrho,rho);} double rt=(ht-lt)/2, rr=(maxrho-minrho)/2; f.R=(maxrho+minrho)/2; f.r=(rt+rr)/2; if(!(f.r>1e-12)||f.R<1.25*f.r||fabs(rt-rr)>0.25*f.r) return f; int stride=max(1,N/250000),cnt=0; long double ss=0; double mx=0; for(int i=0;i<N;i+=stride){double t,u,v;coords_axis(P[i],ax,t,u,v); double rho=hypot(u-f.cu,v-f.cv); double tube=hypot(rho-f.R,t-f.ct); double e=fabs(tube-f.r)/f.r; ss+=e*e; mx=max(mx,e);cnt++;} f.err=cnt?sqrt((double)(ss/cnt)):1e9; f.ok=f.err<0.045 && mx<0.18; return f; }
static Mesh build_torus_param(const TorFit&tf,int U,int Vv){ Mesh m; m.route="TORPARAM"; if(U<8||Vv<6||U*Vv>=N) return m; struct Cell{int id=-1; double score=1e100;}; vector<Cell> cell(U*Vv); for(int id=0; id<N; ++id){ double t,u,v; coords_axis(P[id],tf.ax,t,u,v); double A=atan2(v-tf.cv,u-tf.cu); if(A<0)A+=2*M_PI; double rho=hypot(u-tf.cu,v-tf.cv); double B=atan2(t-tf.ct,rho-tf.R); if(B<0)B+=2*M_PI; int ia=(int)llround(A/(2*M_PI)*U)%U; int ib=(int)llround(B/(2*M_PI)*Vv)%Vv; double ta=2*M_PI*ia/U, tb=2*M_PI*ib/Vv; double da=fabs(A-ta); da=min(da,2*M_PI-da); double db=fabs(B-tb); db=min(db,2*M_PI-db); double sc=da*da+db*db; int k=ia*Vv+ib; if(sc<cell[k].score){cell[k].score=sc; cell[k].id=id;} }
    unordered_set<int> used; used.reserve(U*Vv*2+1); for(auto &c:cell){ if(c.id<0||used.count(c.id)) return {}; used.insert(c.id); m.V.push_back(P[c.id]); }
    auto idx=[&](int i,int j){return ((i%U+U)%U)*Vv+((j%Vv+Vv)%Vv);}; Vec3 ctr=point_axis(tf.ax,tf.ct,tf.cu,tf.cv); for(int i=0;i<U;i++) for(int j=0;j<Vv;j++){ int a=idx(i,j),b=idx(i+1,j),c=idx(i+1,j+1),d=idx(i,j+1); add_face_oriented(m.V,m.F,a,b,c,ctr); add_face_oriented(m.V,m.F,a,c,d,ctr);} return m; }
static Mesh try_torus_param(){ if(N<500||elapsed()>8.5) return {}; for(int ax=0; ax<3; ++ax){ TorFit tf=fit_torus(ax); if(!tf.ok) continue; int Ureq=(int)ceil(2*M_PI*(tf.R+tf.r)/(0.78*0.05*DIAG)); int Vreq=(int)ceil(2*M_PI*tf.r/(0.78*0.05*DIAG)); Ureq=max(12,min(192,Ureq)); Vreq=max(8,min(96,Vreq)); vector<pair<int,int>> trials={{(Ureq*65+99)/100,(Vreq*70+99)/100},{(Ureq*80+99)/100,(Vreq*85+99)/100},{Ureq,Vreq},{min(220,(Ureq*115+99)/100),Vreq}}; sort(trials.begin(),trials.end(),[](auto&a,auto&b){return a.first*a.second<b.first*b.second;}); trials.erase(unique(trials.begin(),trials.end()),trials.end()); for(auto [U,Vv]:trials){ Mesh m=build_torus_param(tf,U,Vv); if(accept_candidate(m,true)) return m; if(elapsed()>11.0) break; } } return {}; }

static Mesh build_star_latlon(const Vec3&ctr,int L,int S){ Mesh m; m.route="STARNET"; if(L<4||S<8||2+(L-1)*S>=N) return m; int nodes=2+(L-1)*S; vector<int> best(nodes,-1); vector<double> bscore(nodes,-1e100); auto node_id=[&](int r,int j){ if(r==0) return 0; if(r==L) return nodes-1; return 1+(r-1)*S+((j%S+S)%S);}; vector<Vec3> dirs(nodes); dirs[0]={0,0,1}; dirs[nodes-1]={0,0,-1}; for(int r=1;r<L;r++){ double th=M_PI*r/L, st=sin(th), ct=cos(th); for(int j=0;j<S;j++){ double ph=2*M_PI*j/S; dirs[node_id(r,j)]={st*cos(ph),st*sin(ph),ct}; }}
    for(int i=0;i<N;i++){ Vec3 q=P[i]-ctr; double rr=norm3(q); if(rr<=1e-12) continue; q=q/rr; double th=acos(clampd(q.z,-1,1)); int r=(int)llround(th/M_PI*L); r=max(0,min(L,r)); int j=0; if(r>0&&r<L){ double ph=atan2(q.y,q.x); if(ph<0)ph+=2*M_PI; j=(int)llround(ph/(2*M_PI)*S)%S; } int id=node_id(r,j); double sc=dot3(q,dirs[id]) + 1e-9*rr/DIAG; if(sc>bscore[id]){bscore[id]=sc;best[id]=i;} }
    for(int id=0;id<nodes;id++) if(best[id]<0) return {}; unordered_set<int> used; used.reserve(nodes*2+1); for(int id=0;id<nodes;id++){ if(used.count(best[id])) return {}; used.insert(best[id]); m.V.push_back(P[best[id]]); }
    Vec3 center=ctr; for(int j=0;j<S;j++) add_face_oriented(m.V,m.F,node_id(0,0),node_id(1,j+1),node_id(1,j),center); for(int r=1;r<L-1;r++) for(int j=0;j<S;j++){ int a=node_id(r,j),b=node_id(r,j+1),c=node_id(r+1,j),d=node_id(r+1,j+1); add_face_oriented(m.V,m.F,a,b,c,center); add_face_oriented(m.V,m.F,b,d,c,center);} for(int j=0;j<S;j++) add_face_oriented(m.V,m.F,node_id(L,0),node_id(L-1,j),node_id(L-1,j+1),center); return m; }
static bool radial_shell_plausible(const Vec3&ctr){ int B=18,S=36; vector<double> mn(B*S,1e100), mx(B*S,0); vector<int> cnt(B*S,0); int stride=max(1,N/300000); for(int i=0;i<N;i+=stride){ Vec3 q=P[i]-ctr; double r=norm3(q); if(r<=1e-12) continue; q=q/r; int bi=min(B-1,max(0,(int)(acos(clampd(q.z,-1,1))/M_PI*B))); double ph=atan2(q.y,q.x); if(ph<0)ph+=2*M_PI; int sj=min(S-1,max(0,(int)(ph/(2*M_PI)*S))); int id=bi*S+sj; mn[id]=min(mn[id],r); mx[id]=max(mx[id],r); cnt[id]++; }
    int good=0,bad=0; for(int i=0;i<B*S;i++) if(cnt[i]>=3){ if(mx[i]-mn[i] <= 0.095*DIAG) good++; else bad++; } return good>=40 && bad*4<good; }
static Mesh try_star_net(const TopoInfo&ti){ if(N<1000||elapsed()>11.5) return {}; if(ti.boundary||ti.nonman||ti.comps!=1) return {};
    vector<Vec3> centers; centers.push_back((bbmin+bbmax)/2.0); Vec3 mean{}; int stride=max(1,N/250000),cnt=0; for(int i=0;i<N;i+=stride){mean=mean+P[i];cnt++;} if(cnt) centers.push_back(mean/(double)cnt);
    double rmax=0; for(int i=0;i<N;i+=stride) rmax=max(rmax,norm3(P[i]-centers[0])); int Lreq=(int)ceil(M_PI*max(rmax,DIAG*0.25)/(0.72*0.05*DIAG)); Lreq=max(8,min(54,Lreq)); vector<pair<int,int>> trials={{max(8,(Lreq*65+99)/100),max(16,2*((Lreq*65+99)/100))},{max(10,(Lreq*80+99)/100),max(20,2*((Lreq*80+99)/100))},{Lreq,2*Lreq},{min(64,(Lreq*115+99)/100),2*min(64,(Lreq*115+99)/100)}};
    for(Vec3 ctr:centers){ if(!radial_shell_plausible(ctr)) continue; for(auto [L,S]:trials){ if(2+(L-1)*S>=N) continue; Mesh m=build_star_latlon(ctr,L,S); if(accept_candidate(m,true)) return m; if(elapsed()>15.0) return {}; } }
    return {};
}

namespace QEM{
struct QFace{int v[3]; unsigned char active=1;};
struct Quadric{ double q[10]; Quadric(){memset(q,0,sizeof q);} void plane(double a,double b,double c,double d,double w=1){q[0]+=w*a*a;q[1]+=w*a*b;q[2]+=w*a*c;q[3]+=w*a*d;q[4]+=w*b*b;q[5]+=w*b*c;q[6]+=w*b*d;q[7]+=w*c*c;q[8]+=w*c*d;q[9]+=w*d*d;} void add(const Quadric&o){for(int i=0;i<10;i++)q[i]+=o.q[i];} double eval(const Vec3&p)const{double x=p.x,y=p.y,z=p.z;return q[0]*x*x+2*q[1]*x*y+2*q[2]*x*z+2*q[3]*x+q[4]*y*y+2*q[5]*y*z+2*q[6]*y+q[7]*z*z+2*q[8]*z+q[9];}};
struct Node{double c;int a,b,va,vb; bool operator<(const Node&o)const{return c>o.c;}};
static vector<Vec3> V,Orig; static vector<QFace> F; static vector<vector<int>> inc; static vector<Quadric> Q; static vector<unsigned char> av; static vector<int> ver,mark,head,tail,nxt,csz; static int activeV,activeF,tok; static double eps2,area_eps2,ncos; static priority_queue<Node> pq;
static bool hasv(const QFace&f,int v){return f.v[0]==v||f.v[1]==v||f.v[2]==v;} static bool hase(const QFace&f,int a,int b){return hasv(f,a)&&hasv(f,b);} static int third(const QFace&f,int a,int b){for(int k=0;k<3;k++){int x=f.v[k]; if(x!=a&&x!=b)return x;}return -1;} static Vec3 fcross(int a,int b,int c){return cross3(V[b]-V[a],V[c]-V[a]);}
static bool solve_opt(const Quadric&q,Vec3&out){ double a00=q.q[0],a01=q.q[1],a02=q.q[2],a11=q.q[4],a12=q.q[5],a22=q.q[7]; double b0=-q.q[3],b1=-q.q[6],b2=-q.q[8]; double det=a00*(a11*a22-a12*a12)-a01*(a01*a22-a12*a02)+a02*(a01*a12-a11*a02); if(fabs(det)<1e-14)return false; double dx=b0*(a11*a22-a12*a12)-a01*(b1*a22-a12*b2)+a02*(b1*a12-a11*b2); double dy=a00*(b1*a22-a12*b2)-b0*(a01*a22-a12*a02)+a02*(a01*b2-b1*a02); double dz=a00*(a11*b2-b1*a12)-a01*(a01*b2-b1*a02)+b0*(a01*a12-a11*a02); out={dx/det,dy/det,dz/det}; return isfinite(out.x)&&isfinite(out.y)&&isfinite(out.z); }
static bool cluster_ok(int v,const Vec3&p){ for(int m=head[v];m!=-1;m=nxt[m]) if(norm2(Orig[m]-p)>eps2) return false; return true; }
static bool shared(int a,int b,int sh[2],int op[2]){ int cnt=0; const auto&small=inc[a].size()<inc[b].size()?inc[a]:inc[b]; for(int id:small){ if(!F[id].active)continue; if(!hase(F[id],a,b))continue; if(cnt>=2)return false; sh[cnt]=id; op[cnt]=third(F[id],a,b); if(op[cnt]<0)return false; cnt++; } return cnt==2&&op[0]!=op[1]; }
static bool link_ok(int a,int b,const int op[2]){ if(tok>2000000000){fill(mark.begin(),mark.end(),0);tok=1;} int ta=tok++, tb=tok++; for(int id:inc[a]) if(F[id].active&&hasv(F[id],a)){ for(int k=0;k<3;k++){int x=F[id].v[k]; if(x!=a&&x!=b)mark[x]=ta;}} int com=0,g0=0,g1=0; for(int id:inc[b]) if(F[id].active&&hasv(F[id],b)){ for(int k=0;k<3;k++){int x=F[id].v[k]; if(x==a||x==b)continue; if(mark[x]==ta){mark[x]=tb; com++; if(x==op[0])g0=1; if(x==op[1])g1=1; if(com>2)return false;}}} return com==2&&g0&&g1; }
static bool dup_after(int keep,int rem,int fid,int a,int b,int c,int s0,int s1){ auto key=tri_key(a,b,c); int probe=keep; if(inc[a].size()<inc[probe].size())probe=a; if(inc[b].size()<inc[probe].size())probe=b; if(inc[c].size()<inc[probe].size())probe=c; for(int of:inc[probe]){ if(!F[of].active||of==fid||of==s0||of==s1)continue; if(hasv(F[of],rem))continue; if(tri_key(F[of].v[0],F[of].v[1],F[of].v[2])==key)return true;} return false; }
static bool pos_ok(int keep,int rem,const int sh[2],const Vec3&pos){ if(!cluster_ok(keep,pos)||!cluster_ok(rem,pos))return false; auto scan=[&](int id,bool rep){ QFace f=F[id]; int nt[3]={f.v[0],f.v[1],f.v[2]}; if(rep)for(int k=0;k<3;k++)if(nt[k]==rem)nt[k]=keep; if(nt[0]==nt[1]||nt[0]==nt[2]||nt[1]==nt[2])return false; Vec3 oldc=fcross(f.v[0],f.v[1],f.v[2]); Vec3 p0=nt[0]==keep?pos:V[nt[0]],p1=nt[1]==keep?pos:V[nt[1]],p2=nt[2]==keep?pos:V[nt[2]]; Vec3 newc=cross3(p1-p0,p2-p0); double o2=norm2(oldc),n2=norm2(newc); if(!(o2>area_eps2)||!(n2>area_eps2))return false; if(dot3(oldc,newc)<=ncos*sqrt(o2*n2))return false; if(rep&&dup_after(keep,rem,id,nt[0],nt[1],nt[2],sh[0],sh[1]))return false; return true; };
    for(int id:inc[keep]) if(F[id].active&&hasv(F[id],keep)&&!hasv(F[id],rem)) if(!scan(id,false)) return false; for(int id:inc[rem]) if(F[id].active&&hasv(F[id],rem)&&id!=sh[0]&&id!=sh[1]&&!hasv(F[id],keep)) if(!scan(id,true)) return false; return true; }
struct Ev{bool ok=false; double c=1e100; Vec3 p;};
static Ev eval_dir(int keep,int rem,const int sh[2]){ Ev best; Quadric q=Q[keep]; q.add(Q[rem]); Vec3 cand[8]; int n=0; Vec3 opt; if(solve_opt(q,opt)) cand[n++]=opt; cand[n++]=(V[keep]+V[rem])*0.5; cand[n++]=V[keep]; cand[n++]=V[rem]; cand[n++]=V[keep]*0.75+V[rem]*0.25; cand[n++]=V[keep]*0.25+V[rem]*0.75; cand[n++]=V[keep]*0.60+V[rem]*0.40; cand[n++]=V[keep]*0.40+V[rem]*0.60; for(int i=0;i<n;i++){ if(!pos_ok(keep,rem,sh,cand[i]))continue; double c=q.eval(cand[i])+1e-8*(norm2(cand[i]-V[keep])+norm2(cand[i]-V[rem])); if(c<best.c){best.ok=true;best.c=c;best.p=cand[i];}} return best; }
static double cheap(int a,int b){ Quadric q=Q[a]; q.add(Q[b]); double best=1e100; Vec3 opt; if(solve_opt(q,opt)&&norm2(opt-V[a])<=eps2&&norm2(opt-V[b])<=eps2) best=min(best,q.eval(opt)); Vec3 mid=(V[a]+V[b])*0.5; if(norm2(mid-V[a])<=eps2&&norm2(mid-V[b])<=eps2) best=min(best,q.eval(mid)); best=min(best,q.eval(V[a])+1e-9*norm2(V[a]-V[b])); best=min(best,q.eval(V[b])+1e-9*norm2(V[a]-V[b])); return best; }
static void push_edge(int a,int b){ if(a==b||a<0||b<0||a>=N||b>=N||!av[a]||!av[b])return; if(norm2(V[a]-V[b])>4.0001*eps2)return; pq.push({cheap(a,b),a,b,ver[a],ver[b]}); }
static void compact(int v){ auto &x=inc[v]; if(x.size()<128)return; size_t alive=0; for(int id:x) if(F[id].active&&hasv(F[id],v)) alive++; if(alive*3+32>=x.size())return; vector<int> y; y.reserve(alive+4); for(int id:x) if(F[id].active&&hasv(F[id],v)) y.push_back(id); x.swap(y); }
static void merge_members(int rem,int keep){ if(head[rem]==-1)return; nxt[tail[keep]]=head[rem]; tail[keep]=tail[rem]; csz[keep]+=csz[rem]; head[rem]=tail[rem]=-1; csz[rem]=0; }
static void apply(int keep,int rem,const int sh[2],const Vec3&pos){ Q[keep].add(Q[rem]); V[keep]=pos; for(int i=0;i<2;i++) if(F[sh[i]].active){F[sh[i]].active=0;activeF--;} for(int id:inc[rem]){ if(!F[id].active)continue; if(!hasv(F[id],rem))continue; if(hasv(F[id],keep)){F[id].active=0;activeF--;continue;} for(int k=0;k<3;k++) if(F[id].v[k]==rem) F[id].v[k]=keep; inc[keep].push_back(id); } av[rem]=0; activeV--; ver[keep]++; ver[rem]++; merge_members(rem,keep); compact(rem); compact(keep); for(int id:inc[keep]) if(F[id].active&&hasv(F[id],keep)){push_edge(F[id].v[0],F[id].v[1]);push_edge(F[id].v[1],F[id].v[2]);push_edge(F[id].v[2],F[id].v[0]);} }
static bool attempt(int a,int b){ if(a==b||!av[a]||!av[b])return false; int sh[2],op[2]; if(!shared(a,b,sh,op))return false; if(!link_ok(a,b,op))return false; Ev ab=eval_dir(a,b,sh), ba=eval_dir(b,a,sh); if(!ab.ok&&!ba.ok)return false; if(ba.ok&&(!ab.ok||ba.c<ab.c)) apply(b,a,sh,ba.p); else apply(a,b,sh,ab.p); return true; }
static void init(const TopoInfo&ti){ V=P; Orig=P; F.resize(M); for(int i=0;i<M;i++){F[i].v[0]=Fin[i].v[0];F[i].v[1]=Fin[i].v[1];F[i].v[2]=Fin[i].v[2];F[i].active=1;} av.assign(N,1); ver.assign(N,0); mark.assign(N,0); head.resize(N);tail.resize(N);nxt.assign(N,-1);csz.assign(N,1); for(int i=0;i<N;i++)head[i]=tail[i]=i; activeV=N; activeF=M; tok=1; eps2=pow(0.05*DIAG*0.999,2); area_eps2=max(1e-32,1e-24*pow(DIAG,4)); ncos=(ti.sharp45>0.10?cos(55*M_PI/180.0):(ti.smooth15>0.985?cos(82*M_PI/180.0):cos(68*M_PI/180.0)));
    vector<int>d(N,0); for(auto&f:F){d[f.v[0]]++;d[f.v[1]]++;d[f.v[2]]++;} inc.assign(N,{}); for(int i=0;i<N;i++)inc[i].reserve(d[i]+8); for(int i=0;i<M;i++){inc[F[i].v[0]].push_back(i);inc[F[i].v[1]].push_back(i);inc[F[i].v[2]].push_back(i);} Q.assign(N,Quadric()); struct ER{uint64_t k;int fid;}; vector<ER> er; er.reserve((size_t)M*3); vector<Vec3> fn(M); for(int i=0;i<M;i++){auto&f=F[i]; Vec3 n=fcross(f.v[0],f.v[1],f.v[2]); norm_inplace(n); fn[i]=n; double dd=-dot3(n,V[f.v[0]]); for(int k=0;k<3;k++) Q[f.v[k]].plane(n.x,n.y,n.z,dd,1.0); er.push_back({edge_key(f.v[0],f.v[1]),i});er.push_back({edge_key(f.v[1],f.v[2]),i});er.push_back({edge_key(f.v[2],f.v[0]),i});}
    sort(er.begin(),er.end(),[](const ER&a,const ER&b){return a.k<b.k;}); double feature=cos(32*M_PI/180.0); while(!pq.empty())pq.pop(); for(size_t i=0;i<er.size();){size_t j=i+1; while(j<er.size()&&er[j].k==er[i].k)j++; int a=(int)(er[i].k>>32),b=(int)(er[i].k&0xffffffffu); if(j-i==2){ double nd=dot3(fn[er[i].fid],fn[er[i+1].fid]); if(nd<feature){ Vec3 e=V[b]-V[a]; if(norm_inplace(e)){ Vec3 p1=cross3(e,fn[er[i].fid]); if(norm_inplace(p1)){double d0=-dot3(p1,V[a]); Q[a].plane(p1.x,p1.y,p1.z,d0,8); Q[b].plane(p1.x,p1.y,p1.z,d0,8);} Vec3 p2=cross3(fn[er[i+1].fid],e); if(norm_inplace(p2)){double d0=-dot3(p2,V[a]); Q[a].plane(p2.x,p2.y,p2.z,d0,8); Q[b].plane(p2.x,p2.y,p2.z,d0,8);} } } } push_edge(a,b); i=j; }
}
static int target_vertices(const TopoInfo&ti){ double r=0.105; if(ti.smooth15>0.985&&ti.sharp45<0.01) r=(N>100000?0.052:0.066); else if(ti.smooth15>0.94) r=(N>100000?0.064:0.083); else if(ti.sharp45>0.12) r=0.135; if(N<1500) r=max(r,0.16); if(N<400) r=max(r,0.30); return max(4,(int)ceil(N*r)); }
static Mesh run(const TopoInfo&ti){ Mesh m; m.route="QEM"; init(ti); int target=target_vertices(ti); long long pops=0; while(activeV>target&&!pq.empty()){ if((++pops&4095)==0 && elapsed()>19.0) break; Node cur=pq.top();pq.pop(); int a=cur.a,b=cur.b; if(a==b||!av[a]||!av[b])continue; if(cur.va!=ver[a]||cur.vb!=ver[b]){push_edge(a,b);continue;} attempt(a,b); }
    vector<int> id(N,-1); m.V.reserve(activeV); for(int i=0;i<N;i++) if(av[i]){id[i]=m.V.size(); m.V.push_back(V[i]);} for(int i=0;i<M;i++) if(F[i].active){int a=F[i].v[0],b=F[i].v[1],c=F[i].v[2]; if(a==b||a==c||b==c)continue; if(id[a]<0||id[b]<0||id[c]<0)continue; m.F.push_back({{id[a],id[b],id[c]}});} return m; }
}

static void print_mesh(const Mesh&m){ printf("%d %d\n",(int)m.V.size(),(int)m.F.size()); for(auto&p:m.V) printf("v %.17g %.17g %.17g\n",p.x,p.y,p.z); for(auto&f:m.F) printf("f %d %d %d\n",f.v[0]+1,f.v[1]+1,f.v[2]+1); }
int main(){ T0=chrono::steady_clock::now(); if(!read_mesh()) return 0; TopoInfo ti=topo_info();
    { Mesh box=make_bbox_box(); if(accept_candidate(box,true)){ print_mesh(box); return 0; } }
    Mesh m;
    m=try_ring_revolution(); if(!m.V.empty()){ print_mesh(m); return 0; }
    m=try_torus_param(); if(!m.V.empty()){ print_mesh(m); return 0; }
    m=try_star_net(ti); if(!m.V.empty()){ print_mesh(m); return 0; }
    m=QEM::run(ti); if((int)m.V.size()<N && !m.F.empty() && closed_manifold_ok(m.V,m.F) && all_near(m.V,P,0.05*DIAG*1.001) && all_near(P,m.V,0.05*DIAG*1.001)){ print_mesh(m); return 0; }
    print_mesh(original_mesh()); return 0; }
