#include <bits/stdc++.h>
using namespace std;

struct Vec3{double x,y,z;};
static inline Vec3 operator+(const Vec3&a,const Vec3&b){return {a.x+b.x,a.y+b.y,a.z+b.z};}
static inline Vec3 operator-(const Vec3&a,const Vec3&b){return {a.x-b.x,a.y-b.y,a.z-b.z};}
static inline Vec3 operator*(const Vec3&a,double s){return {a.x*s,a.y*s,a.z*s};}
static inline Vec3 operator/(const Vec3&a,double s){return {a.x/s,a.y/s,a.z/s};}
static inline double dotv(const Vec3&a,const Vec3&b){return a.x*b.x+a.y*b.y+a.z*b.z;}
static inline Vec3 crossv(const Vec3&a,const Vec3&b){return {a.y*b.z-a.z*b.y,a.z*b.x-a.x*b.z,a.x*b.y-a.y*b.x};}
static inline double norm2(const Vec3&a){return dotv(a,a);} 
static inline double normv(const Vec3&a){return sqrt(max(0.0,norm2(a)));}
static inline Vec3 unit(Vec3 a){double n=normv(a);return n>1e-300?a/n:Vec3{0,0,0};}
struct Tri{int a,b,c;};
struct Mesh{vector<Vec3> v; vector<Tri> f;};

struct FastInput{
    vector<char> buf; char *p;
    FastInput(){ char tmp[1<<16]; size_t n; while((n=fread(tmp,1,sizeof(tmp),stdin))>0) buf.insert(buf.end(),tmp,tmp+n); buf.push_back(0); p=buf.data(); }
    inline void skip(){while(*p==' '||*p=='\n'||*p=='\r'||*p=='\t') ++p;}
    int nextInt(){ skip(); int s=1; if(*p=='-') s=-1,++p; int x=0; while(*p>='0'&&*p<='9') x=x*10+(*p++-'0'); return x*s;}
    double nextDouble(){ skip(); char *e; double x=strtod(p,&e); p=e; return x;}
    char nextChar(){ skip(); return *p++;}
};

static Mesh originalMesh;
static int ORIGN=0, ORIGM=0;
static double DIAG=1.0, TAU=0.05, TAU2=0.0025;
static chrono::steady_clock::time_point T0;
static inline double elapsed(){return chrono::duration<double>(chrono::steady_clock::now()-T0).count();}

static uint64_t edgeKey(int a,int b){ if(a>b) swap(a,b); return (uint64_t)(uint32_t)a<<32 | (uint32_t)b; }
static uint64_t triKey(int a,int b,int c){ if(a>b) swap(a,b); if(b>c) swap(b,c); if(a>b) swap(a,b); return ((uint64_t)(uint32_t)a<<42) ^ ((uint64_t)(uint32_t)b<<21) ^ (uint32_t)c; }

static bool basicValid(const Mesh&m){
    if(m.v.empty()||m.f.empty()||m.v.size()>(size_t)ORIGN) return false;
    const double eps=max(1e-30,1e-24*DIAG*DIAG);
    vector<uint64_t> faces; faces.reserve(m.f.size());
    vector<uint64_t> edges; edges.reserve(m.f.size()*3);
    vector<int> used(m.v.size(),0);
    for(const auto&t:m.f){
        if(t.a<0||t.b<0||t.c<0||t.a>=(int)m.v.size()||t.b>=(int)m.v.size()||t.c>=(int)m.v.size()) return false;
        if(t.a==t.b||t.b==t.c||t.a==t.c) return false;
        Vec3 cr=crossv(m.v[t.b]-m.v[t.a],m.v[t.c]-m.v[t.a]);
        if(!(norm2(cr)>eps)) return false;
        faces.push_back(triKey(t.a,t.b,t.c));
        edges.push_back(edgeKey(t.a,t.b)); edges.push_back(edgeKey(t.b,t.c)); edges.push_back(edgeKey(t.c,t.a));
        used[t.a]=used[t.b]=used[t.c]=1;
    }
    for(int u:used) if(!u) return false;
    sort(faces.begin(),faces.end());
    if(adjacent_find(faces.begin(),faces.end())!=faces.end()) return false;
    sort(edges.begin(),edges.end());
    for(size_t i=0;i<edges.size();){ size_t j=i+1; while(j<edges.size()&&edges[j]==edges[i]) ++j; if(j-i!=2) return false; i=j; }
    return true;
}

struct PointGrid{
    double h,rad2; const vector<Vec3>* pts; unordered_map<long long, vector<int>> mp;
    static long long key(int x,int y,int z){ const long long B=1048576, M=(1LL<<21)-1; return (((long long)(x+B)&M)<<42) | (((long long)(y+B)&M)<<21) | ((long long)(z+B)&M); }
    void build(const vector<Vec3>&p,double r){ pts=&p; h=max(r,1e-12); rad2=r*r*(1.0000001); mp.clear(); mp.reserve(p.size()*2+16); for(int i=0;i<(int)p.size();++i){ int x=floor(p[i].x/h), y=floor(p[i].y/h), z=floor(p[i].z/h); mp[key(x,y,z)].push_back(i);} }
    bool near(const Vec3&q)const{ int x=floor(q.x/h), y=floor(q.y/h), z=floor(q.z/h); int rr=2; for(int dx=-rr;dx<=rr;++dx)for(int dy=-rr;dy<=rr;++dy)for(int dz=-rr;dz<=rr;++dz){ auto it=mp.find(key(x+dx,y+dy,z+dz)); if(it==mp.end()) continue; for(int id:it->second) if(norm2((*pts)[id]-q)<=rad2) return true; } return false; }
};
static bool vertexHausdorffOK(const Mesh&m,double timeLimit=17.5){
    PointGrid g; g.build(m.v,TAU);
    for(int i=0;i<ORIGN;++i){ if((i&8191)==0 && elapsed()>timeLimit) return false; if(!g.near(originalMesh.v[i])) return false; }
    PointGrid go; go.build(originalMesh.v,TAU);
    for(int i=0;i<(int)m.v.size();++i){ if((i&8191)==0 && elapsed()>timeLimit) return false; if(!go.near(m.v[i])) return false; }
    return true;
}

struct Maps{int R=0; vector<float>d,n;};
static inline void projectPoint(const Vec3&p,int view,int R,double&u,double&v,double&z){
    const double D=2.5, f=800.0*R/1024.0, c=0.5*R;
    double sx,sy;
    if(view==0){sx=p.y; sy=p.z; z=D-p.x;}
    else if(view==1){sx=-p.y; sy=p.z; z=D+p.x;}
    else if(view==2){sx=-p.x; sy=p.z; z=D-p.y;}
    else if(view==3){sx=p.x; sy=p.z; z=D+p.y;}
    else if(view==4){sx=p.x; sy=p.y; z=D-p.z;}
    else {sx=-p.x; sy=p.y; z=D+p.z;}
    u=f*sx/z+c; v=f*sy/z+c;
}
static Maps renderMesh(const Mesh&m,int R){
    Maps mp; mp.R=R; int P=R*R; mp.d.assign((size_t)6*P,255.0f); mp.n.assign((size_t)6*P*3,127.5f);
    vector<float> U(m.v.size()), V(m.v.size()), Z(m.v.size());
    vector<Vec3> fn(m.f.size());
    for(int i=0;i<(int)m.f.size();++i){ const Tri&t=m.f[i]; fn[i]=unit(crossv(m.v[t.b]-m.v[t.a],m.v[t.c]-m.v[t.a])); }
    for(int view=0;view<6;++view){
        for(int i=0;i<(int)m.v.size();++i){ double u,v,z; projectPoint(m.v[i],view,R,u,v,z); U[i]=u; V[i]=v; Z[i]=z; }
        float* depth=mp.d.data()+(size_t)view*P; float* norm=mp.n.data()+(size_t)view*P*3;
        for(int ti=0;ti<(int)m.f.size();++ti){ const Tri&t=m.f[ti]; int ia=t.a,ib=t.b,ic=t.c; float z0=Z[ia],z1=Z[ib],z2=Z[ic]; if(!(z0>0&&z1>0&&z2>0)) continue;
            float u0=U[ia],u1=U[ib],u2=U[ic], v0=V[ia],v1=V[ib],v2=V[ic];
            int xmin=max(0,(int)floor(min(u0,min(u1,u2)))); int xmax=min(R-1,(int)ceil(max(u0,max(u1,u2))));
            int ymin=max(0,(int)floor(min(v0,min(v1,v2)))); int ymax=min(R-1,(int)ceil(max(v0,max(v1,v2)))); if(xmin>xmax||ymin>ymax) continue;
            float den=(v1-v2)*(u0-u2)+(u2-u1)*(v0-v2); if(fabs(den)<1e-20f) continue; float inv=1.0f/den;
            Vec3 nn=fn[ti]; float nr=(nn.x+1.0)*127.5, ng=(nn.y+1.0)*127.5, nb=(nn.z+1.0)*127.5;
            for(int y=ymin;y<=ymax;++y){ float py=y+0.5f; int row=y*R; for(int x=xmin;x<=xmax;++x){ float px=x+0.5f; float w0=((v1-v2)*(px-u2)+(u2-u1)*(py-v2))*inv; float w1=((v2-v0)*(px-u2)+(u0-u2)*(py-v2))*inv; float w2=1.0f-w0-w1; if(w0>=-1e-6f&&w1>=-1e-6f&&w2>=-1e-6f){ float dep=1.0f/(w0/z0+w1/z1+w2/z2); int idx=row+x; if(dep<depth[idx]){ depth[idx]=dep; float*q=norm+3*idx; q[0]=nr; q[1]=ng; q[2]=nb; } } } }
        }
    }
    return mp;
}
static inline double rect(const vector<double>&I,int W,int x0,int y0,int x1,int y1){ return I[(size_t)y1*W+x1]-I[(size_t)y0*W+x1]-I[(size_t)y1*W+x0]+I[(size_t)y0*W+x0]; }
static double ssimChannel(const float*A,int strideA,const float*B,int strideB,const float*DA,const float*DB,int R){
    int W=R+1; vector<double> IA((size_t)W*W),IB((size_t)W*W),IA2((size_t)W*W),IB2((size_t)W*W),IAB((size_t)W*W);
    for(int y=1;y<=R;++y){ double sa=0,sb=0,sa2=0,sb2=0,sab=0; int row=(y-1)*R; for(int x=1;x<=R;++x){ int p=row+x-1; double a=A[(size_t)p*strideA],b=B[(size_t)p*strideB]; sa+=a; sb+=b; sa2+=a*a; sb2+=b*b; sab+=a*b; size_t id=(size_t)y*W+x, up=(size_t)(y-1)*W+x; IA[id]=IA[up]+sa; IB[id]=IB[up]+sb; IA2[id]=IA2[up]+sa2; IB2[id]=IB2[up]+sb2; IAB[id]=IAB[up]+sab; } }
    const int rad=5, area=121; const double c1=6.5025,c2=58.5225; long long cnt=0; long double acc=0;
    for(int y=rad;y<R-rad;++y){ int row=y*R; for(int x=rad;x<R-rad;++x){ int p=row+x; if(!(DA[p]<254.0f||DB[p]<254.0f)) continue; int x0=x-rad,y0=y-rad,x1=x+rad+1,y1=y+rad+1; double sa=rect(IA,W,x0,y0,x1,y1), sb=rect(IB,W,x0,y0,x1,y1); double sa2=rect(IA2,W,x0,y0,x1,y1), sb2=rect(IB2,W,x0,y0,x1,y1), sab=rect(IAB,W,x0,y0,x1,y1); double ma=sa/area,mb=sb/area; double va=sa2/area-ma*ma, vb=sb2/area-mb*mb, cov=sab/area-ma*mb; if(va<0&&va>-1e-7)va=0; if(vb<0&&vb>-1e-7)vb=0; double den=(ma*ma+mb*mb+c1)*(va+vb+c2); double val=den?((2*ma*mb+c1)*(2*cov+c2)/den):1.0; acc+=val; ++cnt; } }
    return cnt?(double)(acc/cnt):1.0;
}
static double ssimMaps(const Maps&A,const Maps&B){ int R=A.R, P=R*R; double total=0; for(int view=0;view<6;++view){ const float*ad=A.d.data()+(size_t)view*P; const float*bd=B.d.data()+(size_t)view*P; double ns=0; for(int ch=0;ch<3;++ch){ ns+=ssimChannel(A.n.data()+(size_t)view*P*3+ch,3,B.n.data()+(size_t)view*P*3+ch,3,ad,bd,R); } ns/=3.0; double ds=ssimChannel(ad,1,bd,1,ad,bd,R); total+=0.5*(ns+ds); } return total/6.0; }

static void orientFace(Mesh&m,Tri&t,const Vec3&center){ Vec3 cr=crossv(m.v[t.b]-m.v[t.a],m.v[t.c]-m.v[t.a]); Vec3 ctr=(m.v[t.a]+m.v[t.b]+m.v[t.c])/3.0; if(dotv(cr,ctr-center)<0) swap(t.b,t.c); }
static Mesh makeBox(){
    Vec3 mn=originalMesh.v[0], mx=mn; for(auto&p:originalMesh.v){mn.x=min(mn.x,p.x);mn.y=min(mn.y,p.y);mn.z=min(mn.z,p.z);mx.x=max(mx.x,p.x);mx.y=max(mx.y,p.y);mx.z=max(mx.z,p.z);} Mesh m; m.v={{mn.x,mn.y,mn.z},{mx.x,mn.y,mn.z},{mx.x,mx.y,mn.z},{mn.x,mx.y,mn.z},{mn.x,mn.y,mx.z},{mx.x,mn.y,mx.z},{mx.x,mx.y,mx.z},{mn.x,mx.y,mx.z}}; int F[12][3]={{0,2,1},{0,3,2},{4,5,6},{4,6,7},{0,1,5},{0,5,4},{1,2,6},{1,6,5},{2,3,7},{2,7,6},{3,0,4},{3,4,7}}; for(auto &q:F) m.f.push_back({q[0],q[1],q[2]}); return m; }
static Mesh makeEllipsoid(Vec3 c, Vec3 r, int lat,int lon){ Mesh m; m.v.reserve(2+(lat-1)*lon); m.f.reserve(2*lat*lon); const double PI=acos(-1); m.v.push_back({c.x,c.y,c.z+r.z}); for(int i=1;i<lat;++i){ double th=PI*i/lat, st=sin(th), ct=cos(th); for(int j=0;j<lon;++j){ double ph=2*PI*j/lon; m.v.push_back({c.x+r.x*st*cos(ph),c.y+r.y*st*sin(ph),c.z+r.z*ct}); } } int bot=m.v.size(); m.v.push_back({c.x,c.y,c.z-r.z}); auto id=[&](int ring,int j){return 1+(ring-1)*lon+((j%lon+lon)%lon);}; auto add=[&](int a,int b,int c0){Tri t{a,b,c0}; orientFace(m,t,c); m.f.push_back(t);}; for(int j=0;j<lon;++j)add(0,id(1,j+1),id(1,j)); for(int i=1;i<lat-1;++i)for(int j=0;j<lon;++j){int a=id(i,j),b=id(i,j+1),cc=id(i+1,j),d=id(i+1,j+1); add(a,b,cc); add(b,d,cc);} for(int j=0;j<lon;++j)add(bot,id(lat-1,j),id(lat-1,j+1)); return m; }
static bool ellipsoidFit(Vec3&c,Vec3&r){ Vec3 mn=originalMesh.v[0], mx=mn; for(auto&p:originalMesh.v){mn.x=min(mn.x,p.x);mn.y=min(mn.y,p.y);mn.z=min(mn.z,p.z);mx.x=max(mx.x,p.x);mx.y=max(mx.y,p.y);mx.z=max(mx.z,p.z);} c=(mn+mx)*0.5; r=(mx-mn)*0.5; if(r.x<1e-9||r.y<1e-9||r.z<1e-9) return false; int step=max(1,ORIGN/240000); double ss=0,mxerr=0; int n=0; for(int i=0;i<ORIGN;i+=step){ Vec3 q=originalMesh.v[i]-c; double rr=sqrt((q.x*q.x)/(r.x*r.x)+(q.y*q.y)/(r.y*r.y)+(q.z*q.z)/(r.z*r.z)); double e=fabs(rr-1); ss+=e*e; mxerr=max(mxerr,e); ++n; } double rms=sqrt(ss/max(1,n)); return rms<0.018 && mxerr<0.075; }

static bool acceptSpecial(const Mesh&cand,const Maps&origMaps,double guard,Mesh&best){
    if(cand.v.size()>=best.v.size()||cand.v.size()>=originalMesh.v.size()) return false;
    if(elapsed()>16.6) return false;
    if(!basicValid(cand)) return false;
    if(!vertexHausdorffOK(cand,17.0)) return false;
    if(elapsed()>17.2) return false;
    Maps cm=renderMesh(cand,origMaps.R); double s=ssimMaps(origMaps,cm);
    if(s>=guard){ best=cand; return true; }
    return false;
}

struct Simplifier{
    struct Face{int v[3]; unsigned char active; Vec3 n;};
    struct Quad{double q[10]; Quad(){memset(q,0,sizeof(q));} void add(const Quad&o){for(int i=0;i<10;++i)q[i]+=o.q[i];} void plane(double a,double b,double c,double d,double w){q[0]+=w*a*a;q[1]+=w*a*b;q[2]+=w*a*c;q[3]+=w*a*d;q[4]+=w*b*b;q[5]+=w*b*c;q[6]+=w*b*d;q[7]+=w*c*c;q[8]+=w*c*d;q[9]+=w*d*d;} double eval(Vec3 p)const{double x=p.x,y=p.y,z=p.z;return q[0]*x*x+2*q[1]*x*y+2*q[2]*x*z+2*q[3]*x+q[4]*y*y+2*q[5]*y*z+2*q[6]*y+q[7]*z*z+2*q[8]*z+q[9];}};
    struct Node{double c; int u,v,vu,vv; bool operator<(Node const&o)const{return c>o.c;}};
    int n,m,activeV,activeF; vector<Vec3>P; vector<Face>F; vector<vector<int>> inc; vector<Quad> Q; vector<array<float,3>> mn,mx; vector<unsigned char> alive; vector<int> ver,markA,markB,tmp; int stampA=1,stampB=1; priority_queue<Node> pq; vector<Mesh> snaps;
    Simplifier(const Mesh&src){ n=src.v.size(); m=src.f.size(); P=src.v; F.resize(m); for(int i=0;i<m;++i){F[i].v[0]=src.f[i].a;F[i].v[1]=src.f[i].b;F[i].v[2]=src.f[i].c;F[i].active=1;} }
    bool has(int fid,int v)const{const Face&f=F[fid]; return f.v[0]==v||f.v[1]==v||f.v[2]==v;}
    Vec3 normal(const Face&f)const{return unit(crossv(P[f.v[1]]-P[f.v[0]],P[f.v[2]]-P[f.v[0]]));}
    void clean(int v){ if(!alive[v])return; auto&L=inc[v]; int w=0; for(int fid:L) if(F[fid].active&&has(fid,v)) L[w++]=fid; L.resize(w); }
    void maybe(int v){ if(alive[v]&&inc[v].size()>128){int g=0; for(int fid:inc[v]) if(F[fid].active&&has(fid,v))++g; if((int)inc[v].size()>g*3+64) clean(v);} }
    bool solveOpt(const Quad&q,Vec3&out){ double a00=q.q[0],a01=q.q[1],a02=q.q[2],a11=q.q[4],a12=q.q[5],a22=q.q[7]; double b0=-q.q[3],b1=-q.q[6],b2=-q.q[8]; double det=a00*(a11*a22-a12*a12)-a01*(a01*a22-a12*a02)+a02*(a01*a12-a11*a02); if(fabs(det)<1e-14) return false; double dx=b0*(a11*a22-a12*a12)-a01*(b1*a22-a12*b2)+a02*(b1*a12-a11*b2); double dy=a00*(b1*a22-a12*b2)-b0*(a01*a22-a12*a02)+a02*(a01*b2-b1*a02); double dz=a00*(a11*b2-b1*a12)-a01*(a01*b2-b1*a02)+b0*(a01*a12-a11*a02); out={dx/det,dy/det,dz/det}; return isfinite(out.x)&&isfinite(out.y)&&isfinite(out.z)&&fabs(out.x)<5&&fabs(out.y)<5&&fabs(out.z)<5; }
    void mergedBox(int u,int v,array<float,3>&a,array<float,3>&b){ for(int k=0;k<3;++k){a[k]=min(mn[u][k],mn[v][k]); b[k]=max(mx[u][k],mx[v][k]);} }
    bool coverBox(const Vec3&p,const array<float,3>&a,const array<float,3>&b){ for(int mask=0;mask<8;++mask){ double x=(mask&1)?b[0]:a[0], y=(mask&2)?b[1]:a[1], z=(mask&4)?b[2]:a[2]; if((p.x-x)*(p.x-x)+(p.y-y)*(p.y-y)+(p.z-z)*(p.z-z)>TAU2*0.999) return false; } return true; }
    bool candidate(int u,int v,Vec3&pos,double&cost){ if(u==v||!alive[u]||!alive[v]) return false; array<float,3>a,b; mergedBox(u,v,a,b); Quad q=Q[u]; q.add(Q[v]); Vec3 cand[7]; int cc=0,optok=0; Vec3 opt; if(solveOpt(q,opt)) cand[cc++]=opt; cand[cc++]=(P[u]+P[v])*0.5; cand[cc++]=P[u]; cand[cc++]=P[v]; cand[cc++]={0.5*(a[0]+b[0]),0.5*(a[1]+b[1]),0.5*(a[2]+b[2])}; cand[cc++]=P[u]*0.75+P[v]*0.25; cand[cc++]=P[u]*0.25+P[v]*0.75; bool ok=false; cost=1e300; for(int i=0;i<cc;++i){ if(!coverBox(cand[i],a,b)) continue; double c=max(0.0,q.eval(cand[i]))+1e-8*norm2(P[u]-P[v]); if(c<cost){cost=c; pos=cand[i]; ok=true;} } return ok; }
    void push(int u,int v){ if(u==v||!alive[u]||!alive[v]) return; Vec3 p; double c; if(candidate(u,v,p,c)) pq.push({c,min(u,v),max(u,v),ver[min(u,v)],ver[max(u,v)]}); }
    bool linkOK(int u,int v){ if(!alive[u]||!alive[v]||u==v) return false; maybe(u); maybe(v); int opp[2],ec=0; const vector<int>&small=inc[u].size()<inc[v].size()?inc[u]:inc[v]; for(int fid:small){ if(!F[fid].active) continue; bool hu=has(fid,u),hv=has(fid,v); if(hu&&hv){ if(ec>=2) return false; int o=-1; for(int k=0;k<3;++k){int x=F[fid].v[k]; if(x!=u&&x!=v)o=x;} if(o<0)return false; opp[ec++]=o; } } if(ec!=2||opp[0]==opp[1]) return false; if(++stampA>1000000000){fill(markA.begin(),markA.end(),0);stampA=1;} if(++stampB>1000000000){fill(markB.begin(),markB.end(),0);stampB=1;} for(int fid:inc[u]) if(F[fid].active&&has(fid,u)){ for(int k=0;k<3;++k){int x=F[fid].v[k]; if(x!=u&&x!=v) markA[x]=stampA; }} int cnt=0,got0=0,got1=0; for(int fid:inc[v]) if(F[fid].active&&has(fid,v)){ for(int k=0;k<3;++k){int x=F[fid].v[k]; if(x==u||x==v)continue; if(markA[x]==stampA&&markB[x]!=stampB){ markB[x]=stampB; ++cnt; if(x==opp[0])got0=1; if(x==opp[1])got1=1; if(cnt>2)return false; } }} return cnt==2&&got0&&got1; }
    bool sameTri(const Face&f,int a,int b,int c){ int x[3]={f.v[0],f.v[1],f.v[2]}, y[3]={a,b,c}; sort(x,x+3); sort(y,y+3); return x[0]==y[0]&&x[1]==y[1]&&x[2]==y[2]; }
    bool existsFace(int a,int b,int c,const vector<int>&skip){ int ix=a; if(inc[b].size()<inc[ix].size()) ix=b; if(inc[c].size()<inc[ix].size()) ix=c; for(int fid:inc[ix]){ if(!F[fid].active) continue; if(find(skip.begin(),skip.end(),fid)!=skip.end()) continue; if(sameTri(F[fid],a,b,c)) return true; } return false; }
    bool flipOK(int keep,int rem,const Vec3&np){ vector<int> touched; for(int fid:inc[keep]) if(F[fid].active&&has(fid,keep)) touched.push_back(fid); for(int fid:inc[rem]) if(F[fid].active&&has(fid,rem)) touched.push_back(fid); sort(touched.begin(),touched.end()); touched.erase(unique(touched.begin(),touched.end()),touched.end()); double minDot = ORIGN<3000?-0.05:0.05; double eps=max(1e-30,1e-24*DIAG*DIAG); for(int fid:touched){ Face&f=F[fid]; bool hk=has(fid,keep),hr=has(fid,rem); if(hk&&hr) continue; int a=f.v[0],b=f.v[1],c=f.v[2]; if(a==rem)a=keep; if(b==rem)b=keep; if(c==rem)c=keep; if(a==b||a==c||b==c) return false; Vec3 pa=(a==keep?np:P[a]), pb=(b==keep?np:P[b]), pc=(c==keep?np:P[c]); Vec3 cr=crossv(pb-pa,pc-pa); if(!(norm2(cr)>eps)) return false; if(dotv(unit(cr),f.n)<minDot) return false; if(hr&&!hk&&existsFace(a,b,c,touched)) return false; } return true; }
    void collectN(int u,vector<int>&out){ out.clear(); if(!alive[u]) return; maybe(u); if(++stampA>1000000000){fill(markA.begin(),markA.end(),0);stampA=1;} for(int fid:inc[u]) if(F[fid].active&&has(fid,u)){ for(int k=0;k<3;++k){int w=F[fid].v[k]; if(w!=u&&alive[w]&&markA[w]!=stampA){markA[w]=stampA;out.push_back(w);}} } }
    void collapse(int keep,int rem,const Vec3&np){ clean(keep); clean(rem); vector<int> R=inc[rem]; for(int fid:R){ Face&f=F[fid]; if(!f.active||!has(fid,rem)) continue; bool hk=has(fid,keep); if(hk){f.active=0;--activeF;} else {for(int k=0;k<3;++k) if(f.v[k]==rem) f.v[k]=keep; if(f.v[0]==f.v[1]||f.v[0]==f.v[2]||f.v[1]==f.v[2]){f.active=0;--activeF;} else inc[keep].push_back(fid);} }
        P[keep]=np; Q[keep].add(Q[rem]); for(int k=0;k<3;++k){mn[keep][k]=min(mn[keep][k],mn[rem][k]); mx[keep][k]=max(mx[keep][k],mx[rem][k]);}
        alive[rem]=0; --activeV; ++ver[keep]; ++ver[rem]; inc[rem].clear(); clean(keep); for(int fid:inc[keep]) if(F[fid].active&&has(fid,keep)) F[fid].n=normal(F[fid]); vector<int> nb; collectN(keep,nb); for(int w:nb) push(keep,w); }
    Mesh snapshot(){ Mesh s; tmp.assign(n,-1); vector<int> use(n,0); for(int i=0;i<m;++i) if(F[i].active){for(int k=0;k<3;++k) if(alive[F[i].v[k]]) use[F[i].v[k]]=1;} for(int i=0;i<n;++i) if(use[i]){tmp[i]=s.v.size(); s.v.push_back(P[i]);} for(int i=0;i<m;++i) if(F[i].active){ int a=F[i].v[0],b=F[i].v[1],c=F[i].v[2]; if(a!=b&&a!=c&&b!=c&&tmp[a]>=0&&tmp[b]>=0&&tmp[c]>=0) s.f.push_back({tmp[a],tmp[b],tmp[c]}); } return s; }
    void init(){ inc.assign(n,{}); vector<int> deg(n); for(auto&f:F){++deg[f.v[0]];++deg[f.v[1]];++deg[f.v[2]];} for(int i=0;i<n;++i) inc[i].reserve(deg[i]+8); for(int i=0;i<m;++i){F[i].active=1; F[i].n=normal(F[i]); inc[F[i].v[0]].push_back(i); inc[F[i].v[1]].push_back(i); inc[F[i].v[2]].push_back(i);} activeV=n; activeF=m; Q.assign(n,Quad()); for(int i=0;i<m;++i){ Face&f=F[i]; Vec3 cr=crossv(P[f.v[1]]-P[f.v[0]],P[f.v[2]]-P[f.v[0]]); double ar=normv(cr); if(ar<=1e-300) continue; Vec3 nn=cr/ar; double d=-dotv(nn,P[f.v[0]]); double w=max(1e-8,0.5*ar); for(int k=0;k<3;++k) Q[f.v[k]].plane(nn.x,nn.y,nn.z,d,w); } mn.resize(n); mx.resize(n); for(int i=0;i<n;++i) mn[i]=mx[i]={(float)P[i].x,(float)P[i].y,(float)P[i].z}; alive.assign(n,1); ver.assign(n,0); markA.assign(n,0); markB.assign(n,0); vector<uint64_t> ed; ed.reserve((size_t)m*3); for(auto&f:F){ed.push_back(edgeKey(f.v[0],f.v[1])); ed.push_back(edgeKey(f.v[1],f.v[2])); ed.push_back(edgeKey(f.v[2],f.v[0]));} sort(ed.begin(),ed.end()); ed.erase(unique(ed.begin(),ed.end()),ed.end()); for(uint64_t k:ed) push(k>>32,k&0xffffffffu); }
    vector<Mesh> run(){ init(); vector<double> ratios={0.50,0.35,0.25,0.22,0.20,0.18,0.16,0.14,0.12,0.10,0.09,0.08}; vector<int> cps; for(double r:ratios){int c=max(4,(int)ceil(n*r)); if(c<n&&(cps.empty()||cps.back()!=c)) cps.push_back(c);} int cp=0; int target=cps.empty()?max(4,n/5):cps.back(); long long pops=0; while(activeV>target&&!pq.empty()){ if((++pops&4095)==0&&elapsed()>14.2) break; Node nd=pq.top(); pq.pop(); int u=nd.u,v=nd.v; if(!alive[u]||!alive[v]||ver[u]!=nd.vu||ver[v]!=nd.vv) continue; Vec3 pos; double c; if(!candidate(u,v,pos,c)) continue; if(!linkOK(u,v)) continue; int keep=u,rem=v; if(inc[v].size()>inc[u].size()) keep=v,rem=u; if(!flipOK(keep,rem,pos)) continue; collapse(keep,rem,pos); while(cp<(int)cps.size()&&activeV<=cps[cp]){ snaps.push_back(snapshot()); ++cp; } }
        if(snaps.empty()||snaps.back().v.size()!=(size_t)activeV) snaps.push_back(snapshot()); return snaps; }
};

static bool detectLatLong(Mesh&out){
    if(ORIGN<50 || ORIGM!=2*(ORIGN-2)) return false; int lon=0, rings=0; auto same=[&](const Tri&t,int a,int b,int c){int x[3]={t.a,t.b,t.c},y[3]={a,b,c}; sort(x,x+3); sort(y,y+3); return x[0]==y[0]&&x[1]==y[1]&&x[2]==y[2];};
    for(int L=8;L<=1024;++L){ if((ORIGN-2)%L) continue; int R=(ORIGN-2)/L; if(R<3) continue; bool ok=true; int st=max(1,L/64); for(int j=0;j<L&&ok;j+=st){ int a=1+j,b=1+(j+1)%L; if(!same(originalMesh.f[j],0,a,b)) ok=false; } if(ok){lon=L;rings=R;break;} }
    if(!lon) return false; int R2=max(4,min(rings,rings/6)), L2=max(12,min(lon,lon/6)); if(2+R2*L2>=ORIGN) return false; out.v.clear(); out.f.clear(); out.v.push_back(originalMesh.v[0]); for(int i=0;i<R2;++i){ int oi=1+(long long)i*(rings-1)/max(1,R2-1); for(int j=0;j<L2;++j){ int oj=(long long)j*lon/L2; out.v.push_back(originalMesh.v[1+(oi-1)*lon+oj]); }} int bot=out.v.size(); out.v.push_back(originalMesh.v.back()); auto id=[&](int r,int j){return 1+(r-1)*L2+((j%L2+L2)%L2);}; Vec3 ctr{}; for(auto&p:out.v) ctr=ctr+p; ctr=ctr/(double)out.v.size(); auto add=[&](int a,int b,int c){Tri t{a,b,c}; orientFace(out,t,ctr); out.f.push_back(t);}; for(int j=0;j<L2;++j)add(0,id(1,j),id(1,j+1)); for(int r=1;r<R2;++r)for(int j=0;j<L2;++j){int a=id(r,j),b=id(r,j+1),c=id(r+1,j),d=id(r+1,j+1); add(a,c,b); add(b,c,d);} for(int j=0;j<L2;++j)add(bot,id(R2,j+1),id(R2,j)); return true;
}

static bool detectPeriodicGrid(Mesh&out){
    if(ORIGN<100||ORIGM!=2*ORIGN) return false; auto same=[&](const Tri&t,int a,int b,int c){int x[3]={t.a,t.b,t.c},y[3]={a,b,c}; sort(x,x+3); sort(y,y+3); return x[0]==y[0]&&x[1]==y[1]&&x[2]==y[2];}; int U=0,V=0; for(int v=6;v<=512;++v){ if(ORIGN%v) continue; int u=ORIGN/v; if(u<6) continue; int st=max(1,ORIGN/400),tot=0,ok=0; for(int cell=0;cell<ORIGN&&tot<400;cell+=st){int i=cell/v,j=cell%v; int a=i*v+j,b=i*v+(j+1)%v,c=((i+1)%u)*v+j,d=((i+1)%u)*v+(j+1)%v; int f=2*cell; if(f+1>=ORIGM) break; ok += (same(originalMesh.f[f],a,b,c)||same(originalMesh.f[f],a,c,b)||same(originalMesh.f[f],a,b,d)||same(originalMesh.f[f],b,c,d)); ++tot;} if(tot>100&&ok*100>=tot*92){U=u;V=v;break;} }
    if(!U) return false; int U2=max(8,min(U,U/4)), V2=max(8,min(V,V/4)); if(U2*V2>=ORIGN) return false; out.v.clear(); out.f.clear(); for(int i=0;i<U2;++i){int oi=(long long)i*U/U2; for(int j=0;j<V2;++j){int oj=(long long)j*V/V2; out.v.push_back(originalMesh.v[oi*V+oj]);}} auto id=[&](int i,int j){return ((i%U2+U2)%U2)*V2+((j%V2+V2)%V2);}; Vec3 ctr{}; for(auto&p:out.v)ctr=ctr+p; ctr=ctr/(double)out.v.size(); auto add=[&](int a,int b,int c){Tri t{a,b,c}; orientFace(out,t,ctr); out.f.push_back(t);}; for(int i=0;i<U2;++i)for(int j=0;j<V2;++j){int a=id(i,j),b=id(i,j+1),c=id(i+1,j),d=id(i+1,j+1); add(a,c,b); add(b,c,d);} return true;
}


static inline void coordAxis(const Vec3&p,int ax,double&t,double&u,double&v){ if(ax==0){t=p.x;u=p.y;v=p.z;} else if(ax==1){t=p.y;u=p.x;v=p.z;} else {t=p.z;u=p.x;v=p.y;} }
static inline Vec3 makeAxis(int ax,double t,double u,double v){ if(ax==0)return {t,u,v}; if(ax==1)return {u,t,v}; return {u,v,t}; }

struct TorusFit{bool ok=false;int ax=2;double ct=0,cu=0,cv=0,R=1,r=0.2,err=1e9;};
static TorusFit fitTorusAxis(int ax){
    TorusFit f; f.ax=ax; if(ORIGN<200) return f; double t0=1e100,t1=-1e100,u0=1e100,u1=-1e100,v0=1e100,v1=-1e100;
    for(auto&p:originalMesh.v){double t,u,v; coordAxis(p,ax,t,u,v); t0=min(t0,t);t1=max(t1,t);u0=min(u0,u);u1=max(u1,u);v0=min(v0,v);v1=max(v1,v);} f.ct=(t0+t1)*.5; f.cu=(u0+u1)*.5; f.cv=(v0+v1)*.5; double minrho=1e100,maxrho=0; for(auto&p:originalMesh.v){double t,u,v; coordAxis(p,ax,t,u,v); double rho=hypot(u-f.cu,v-f.cv); minrho=min(minrho,rho); maxrho=max(maxrho,rho);} if(!(maxrho>minrho&&minrho>1e-8)) return f; double R=.5*(maxrho+minrho), rr=.25*((maxrho-minrho)+(t1-t0)); if(!(R>1e-8&&rr>1e-8&&R>1.25*rr)) return f; if(fabs(.5*(maxrho-minrho)-.5*(t1-t0))>0.35*rr) return f; int step=max(1,ORIGN/240000); double ss=0,mxerr=0; int n=0; for(int i=0;i<ORIGN;i+=step){double t,u,v; coordAxis(originalMesh.v[i],ax,t,u,v); double rho=hypot(u-f.cu,v-f.cv); double e=fabs(hypot(rho-R,t-f.ct)-rr)/rr; ss+=e*e; mxerr=max(mxerr,e); ++n;} double rms=sqrt(ss/max(1,n)); if(rms>.022||mxerr>.09) return f; f.ok=true; f.R=R; f.r=rr; f.err=rms; return f;
}
static Mesh makeTorus(const TorusFit&f,int U,int V){ Mesh m; const double PI=acos(-1); m.v.reserve(U*V); m.f.reserve(2*U*V); for(int i=0;i<U;++i){double th=2*PI*i/U, ct=cos(th),st=sin(th); for(int j=0;j<V;++j){double ph=2*PI*j/V, cp=cos(ph),sp=sin(ph); double rho=f.R+f.r*cp; m.v.push_back(makeAxis(f.ax,f.ct+f.r*sp,f.cu+rho*ct,f.cv+rho*st));}} auto id=[&](int i,int j){return ((i%U+U)%U)*V+((j%V+V)%V);}; Vec3 cen=makeAxis(f.ax,f.ct,f.cu,f.cv); auto add=[&](int a,int b,int c){Tri t{a,b,c}; orientFace(m,t,cen); m.f.push_back(t);}; for(int i=0;i<U;++i)for(int j=0;j<V;++j){int a=id(i,j),b=id(i+1,j),c=id(i+1,j+1),d=id(i,j+1); add(a,b,c); add(a,c,d);} return m; }
static bool tryTorus(const Maps&origMaps,Mesh&best){ if(ORIGN<400||elapsed()>13.5) return false; TorusFit bestf; for(int ax=0;ax<3;++ax){auto f=fitTorusAxis(ax); if(f.ok&&(!bestf.ok||f.err<bestf.err)) bestf=f;} if(!bestf.ok) return false; int U=max(24,min(160,(int)ceil(2*acos(-1)*(bestf.R+bestf.r)/max(TAU*.92,1e-9)))); int V=max(8,min(64,(int)ceil(2*acos(-1)*bestf.r/max(TAU*.92,1e-9)))); bool ok=false; for(int scale=0;scale<3&&elapsed()<16.2;++scale){int u=min(192,(int)(U*(1.0+0.18*scale))), v=min(80,(int)(V*(1.0+0.18*scale))); Mesh m=makeTorus(bestf,u,v); if(m.v.size()<best.v.size()) ok|=acceptSpecial(m,origMaps,0.918,best);} return ok; }

struct RevFit{bool ok=false;int ax=2;double t0,t1,cu,cv,r0,r1,err=1e9;};
static RevFit fitRevAxis(int ax){ RevFit f; f.ax=ax; if(ORIGN<200) return f; double t0=1e100,t1=-1e100,u0=1e100,u1=-1e100,v0=1e100,v1=-1e100; for(auto&p:originalMesh.v){double t,u,v; coordAxis(p,ax,t,u,v); t0=min(t0,t);t1=max(t1,t);u0=min(u0,u);u1=max(u1,u);v0=min(v0,v);v1=max(v1,v);} if(!(t1>t0))return f; f.t0=t0;f.t1=t1;f.cu=(u0+u1)*.5;f.cv=(v0+v1)*.5; double maxr=0; for(auto&p:originalMesh.v){double t,u,v;coordAxis(p,ax,t,u,v);maxr=max(maxr,hypot(u-f.cu,v-f.cv));} if(maxr<1e-8||(t1-t0)<0.25*maxr)return f; int step=max(1,ORIGN/240000); double eps=maxr*.06,S=0,St=0,Sr=0,Stt=0,Str=0; int cnt=0,axisCnt=0; for(int i=0;i<ORIGN;i+=step){double t,u,v;coordAxis(originalMesh.v[i],ax,t,u,v);double r=hypot(u-f.cu,v-f.cv); if(r<eps){ if(min(fabs(t-t0),fabs(t-t1))>(t1-t0)*.06) return f; ++axisCnt; continue;} S++; St+=t; Sr+=r; Stt+=t*t; Str+=t*r; ++cnt;} if(cnt<100) return f; double den=S*Stt-St*St; if(fabs(den)<1e-18)return f; double a=(S*Str-St*Sr)/den,b=(Sr-a*St)/S; double r0=max(0.0,a*t0+b),r1=max(0.0,a*t1+b); if(max(r0,r1)<0.2*maxr) return f; double ss=0,mxerr=0; int n=0; for(int i=0;i<ORIGN;i+=step){double t,u,v;coordAxis(originalMesh.v[i],ax,t,u,v);double r=hypot(u-f.cu,v-f.cv); double pred=max(0.0,a*t+b); if(r<eps) continue; double e=fabs(r-pred)/maxr; ss+=e*e; mxerr=max(mxerr,e); ++n;} double rms=sqrt(ss/max(1,n)); if(rms>.018||mxerr>.07) return f; f.ok=true;f.r0=r0;f.r1=r1;f.err=rms;return f; }
static Mesh makeRevolution(const RevFit&f,int sides,int rings){ Mesh m; const double PI=acos(-1); bool cone0=f.r0<TAU*.2, cone1=f.r1<TAU*.2; vector<int> start; Vec3 cen=makeAxis(f.ax,(f.t0+f.t1)*.5,f.cu,f.cv); auto addRing=[&](double t,double r){int st=m.v.size(); start.push_back(st); for(int j=0;j<sides;++j){double th=2*PI*j/sides; m.v.push_back(makeAxis(f.ax,t,f.cu+r*cos(th),f.cv+r*sin(th)));}}; if(cone0){m.v.push_back(makeAxis(f.ax,f.t0,f.cu,f.cv));} else addRing(f.t0,f.r0); for(int i=1;i<rings-1;++i){double s=(double)i/(rings-1), t=f.t0+(f.t1-f.t0)*s, r=f.r0+(f.r1-f.r0)*s; if(r>TAU*.15) addRing(t,r);} int topApex=-1; if(cone1){topApex=m.v.size(); m.v.push_back(makeAxis(f.ax,f.t1,f.cu,f.cv));} else addRing(f.t1,f.r1); auto add=[&](int a,int b,int c){Tri t{a,b,c}; orientFace(m,t,cen); m.f.push_back(t);}; auto rid=[&](int r,int j){return start[r]+((j%sides+sides)%sides);}; int rs=start.size(); if(cone0&&rs>0){for(int j=0;j<sides;++j)add(0,rid(0,j+1),rid(0,j));} for(int r=0;r+1<rs;++r)for(int j=0;j<sides;++j){int a=rid(r,j),b=rid(r,j+1),c=rid(r+1,j),d=rid(r+1,j+1);add(a,c,b);add(b,c,d);} if(cone1&&rs>0){for(int j=0;j<sides;++j)add(topApex,rid(rs-1,j),rid(rs-1,j+1));} if(!cone0&&rs>0){int c0=m.v.size(); m.v.push_back(makeAxis(f.ax,f.t0,f.cu,f.cv)); for(int j=0;j<sides;++j)add(c0,rid(0,j),rid(0,j+1));} if(!cone1&&rs>0){int c1=m.v.size(); m.v.push_back(makeAxis(f.ax,f.t1,f.cu,f.cv)); for(int j=0;j<sides;++j)add(c1,rid(rs-1,j+1),rid(rs-1,j));} return m; }

struct CapFit{bool ok=false;int ax=2;double ct,cu,cv,r,h,err=1e9;};
static CapFit fitCapsuleAxis(int ax){ CapFit f; f.ax=ax; if(ORIGN<300)return f; double t0=1e100,t1=-1e100,u0=1e100,u1=-1e100,v0=1e100,v1=-1e100; for(auto&p:originalMesh.v){double t,u,v;coordAxis(p,ax,t,u,v);t0=min(t0,t);t1=max(t1,t);u0=min(u0,u);u1=max(u1,u);v0=min(v0,v);v1=max(v1,v);} double eu=u1-u0,ev=v1-v0,et=t1-t0; if(!(eu>1e-9&&ev>1e-9&&et>1e-9))return f; if(min(eu,ev)<max(eu,ev)*.94)return f; double r=.25*(eu+ev), h=.5*et-r; if(!(r>1e-9&&h>0.15*r))return f; f.ct=.5*(t0+t1);f.cu=.5*(u0+u1);f.cv=.5*(v0+v1);f.r=r;f.h=h; int step=max(1,ORIGN/240000); double ss=0,mxe=0; int n=0,side=0,cap=0; for(int i=0;i<ORIGN;i+=step){double t,u,v;coordAxis(originalMesh.v[i],ax,t,u,v);double at=fabs(t-f.ct),rho=hypot(u-f.cu,v-f.cv),e; if(at<=h){e=fabs(rho-r)/r; ++side;} else {e=fabs(hypot(rho,at-h)-r)/r; ++cap;} ss+=e*e;mxe=max(mxe,e);++n;} double rms=sqrt(ss/max(1,n)); if(side<max(20,n/20)||cap<max(20,n/20)||rms>.025||mxe>.09)return f; f.ok=true;f.err=rms;return f; }
static Mesh makeCapsule(const CapFit&f,int sides,int capN,int cylN){ Mesh m; const double PI=acos(-1); vector<int> st; Vec3 cen=makeAxis(f.ax,f.ct,f.cu,f.cv); m.v.push_back(makeAxis(f.ax,f.ct+f.h+f.r,f.cu,f.cv)); auto ring=[&](double t,double rr){int s=m.v.size();st.push_back(s);for(int j=0;j<sides;++j){double th=2*PI*j/sides;m.v.push_back(makeAxis(f.ax,t,f.cu+rr*cos(th),f.cv+rr*sin(th)));}}; for(int i=1;i<=capN;++i){double phi=(PI/2)*i/capN; ring(f.ct+f.h+f.r*cos(phi),f.r*sin(phi));} for(int i=1;i<cylN;++i){double a=(double)i/cylN; ring(f.ct+f.h*(1-2*a),f.r);} for(int i=1;i<capN;++i){double phi=PI/2+(PI/2)*i/capN; ring(f.ct-f.h+f.r*cos(phi),f.r*sin(phi));} int bot=m.v.size(); m.v.push_back(makeAxis(f.ax,f.ct-f.h-f.r,f.cu,f.cv)); auto id=[&](int r,int j){return st[r]+((j%sides+sides)%sides);}; auto add=[&](int a,int b,int c){Tri t{a,b,c};orientFace(m,t,cen);m.f.push_back(t);}; int R=st.size(); for(int j=0;j<sides;++j)add(0,id(0,j+1),id(0,j)); for(int r=0;r+1<R;++r)for(int j=0;j<sides;++j){int a=id(r,j),b=id(r,j+1),c=id(r+1,j),d=id(r+1,j+1);add(a,c,b);add(b,c,d);} for(int j=0;j<sides;++j)add(bot,id(R-1,j),id(R-1,j+1)); return m; }
static bool tryCapsule(const Maps&origMaps,Mesh&best){ if(ORIGN<500||elapsed()>14.8)return false; CapFit bf; for(int ax=0;ax<3;++ax){auto f=fitCapsuleAxis(ax); if(f.ok&&(!bf.ok||f.err<bf.err))bf=f;} if(!bf.ok)return false; int sides=max(18,min(128,(int)ceil(2*acos(-1)*bf.r/max(TAU*.9,1e-9)))); int capN=max(4,min(32,(int)ceil((acos(-1)/2*bf.r)/max(TAU*.9,1e-9)))); int cylN=max(1,min(80,(int)ceil((2*bf.h)/max(TAU*.9,1e-9)))); Mesh m=makeCapsule(bf,sides,capN,cylN); if(m.v.size()>=best.v.size())return false; return acceptSpecial(m,origMaps,0.925,best); }

static bool tryRevolution(const Maps&origMaps,Mesh&best){ if(ORIGN<400||elapsed()>14.0)return false; RevFit bf; for(int ax=0;ax<3;++ax){auto f=fitRevAxis(ax); if(f.ok&&(!bf.ok||f.err<bf.err))bf=f;} if(!bf.ok)return false; double maxr=max(bf.r0,bf.r1), len=bf.t1-bf.t0; int sides=max(18,min(128,(int)ceil(2*acos(-1)*maxr/max(TAU*.9,1e-9)))); int rings=max(2,min(80,(int)ceil(len/max(TAU*.9,1e-9))+1)); Mesh m=makeRevolution(bf,sides,rings); if(m.v.size()>=best.v.size())return false; return acceptSpecial(m,origMaps,0.925,best); }

static void writeMesh(const Mesh&m){
    static char obuf[1<<20]; setvbuf(stdout,obuf,_IOFBF,sizeof(obuf));
    printf("%d %d\n",(int)m.v.size(),(int)m.f.size());
    for(auto&p:m.v) printf("v %.12g %.12g %.12g\n",p.x,p.y,p.z);
    for(auto&t:m.f) printf("f %d %d %d\n",t.a+1,t.b+1,t.c+1);
}

int main(){
    T0=chrono::steady_clock::now();
    FastInput in; ORIGN=in.nextInt(); ORIGM=in.nextInt(); originalMesh.v.resize(ORIGN); originalMesh.f.resize(ORIGM); Vec3 mn{1e100,1e100,1e100},mx{-1e100,-1e100,-1e100};
    for(int i=0;i<ORIGN;++i){ in.nextChar(); originalMesh.v[i].x=in.nextDouble(); originalMesh.v[i].y=in.nextDouble(); originalMesh.v[i].z=in.nextDouble(); mn.x=min(mn.x,originalMesh.v[i].x); mn.y=min(mn.y,originalMesh.v[i].y); mn.z=min(mn.z,originalMesh.v[i].z); mx.x=max(mx.x,originalMesh.v[i].x); mx.y=max(mx.y,originalMesh.v[i].y); mx.z=max(mx.z,originalMesh.v[i].z); }
    for(int i=0;i<ORIGM;++i){ in.nextChar(); int a=in.nextInt()-1,b=in.nextInt()-1,c=in.nextInt()-1; originalMesh.f[i]={a,b,c}; }
    DIAG=normv(mx-mn); if(!(DIAG>0)) DIAG=1; TAU=0.05*DIAG*0.999999; TAU2=TAU*TAU;
    Mesh best=originalMesh;
    int R=(ORIGN<20000?256:(ORIGN<120000?192:128));
    Maps origMaps=renderMesh(originalMesh,R);
    if(elapsed()<12.0){ Mesh box=makeBox(); acceptSpecial(box,origMaps,ORIGN<2000?0.955:0.965,best); }
    if(elapsed()<13.0){ Vec3 c,r; if(ellipsoidFit(c,r)){ vector<pair<int,int>> trials={{12,24},{16,32},{20,40},{24,48},{28,56},{32,64}}; for(auto [la,lo]:trials){ if(elapsed()>15.8) break; Mesh e=makeEllipsoid(c,r,la,lo); if(e.v.size()<best.v.size()) acceptSpecial(e,origMaps,0.925,best); } } }
    if(elapsed()<13.2){ Mesh g; if(detectLatLong(g)) acceptSpecial(g,origMaps,0.92,best); }
    if(elapsed()<13.5){ Mesh g; if(detectPeriodicGrid(g)) acceptSpecial(g,origMaps,0.92,best); }
    if(elapsed()<14.0) tryTorus(origMaps,best);
    if(elapsed()<14.4) tryRevolution(origMaps,best);
    if(elapsed()<14.8) tryCapsule(origMaps,best);
    if(best.v.size()==originalMesh.v.size() && elapsed()<15.0){ Simplifier sim(originalMesh); vector<Mesh> snaps=sim.run(); double guard=(R>=256?0.918:0.925); Mesh chosen=best; double bestRatio=1.0; for(auto &s:snaps){ if(elapsed()>18.4) break; if(s.v.empty()||s.f.empty()||s.v.size()>=chosen.v.size()) continue; if(!basicValid(s)) continue; Maps cm=renderMesh(s,R); double sc=ssimMaps(origMaps,cm); if(sc>=guard && (double)s.v.size()/ORIGN<bestRatio){ chosen=move(s); bestRatio=(double)chosen.v.size()/ORIGN; } }
        best=chosen;
    }
    if(!basicValid(best)) best=originalMesh;
    writeMesh(best);
    return 0;
}
