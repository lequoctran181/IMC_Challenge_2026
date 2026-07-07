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
static inline Vec3 unit(Vec3 a){double l=normv(a); return l>1e-300?a/l:Vec3{0,0,0};}
struct Face{int a,b,c;};
struct Mesh{vector<Vec3> v; vector<Face> f;};
static inline unsigned long long ekey(int a,int b){if(a>b)swap(a,b);return (unsigned long long)(unsigned int)a<<32 | (unsigned int)b;}
static inline array<int,3> tkey(int a,int b,int c){array<int,3> r{a,b,c}; sort(r.begin(),r.end()); return r;}
static inline bool same_tri(const Face&f,int a,int b,int c){return tkey(f.a,f.b,f.c)==tkey(a,b,c);} 

int N,M; Mesh inmesh; Vec3 bb0,bb1,center0; double diagLen=1,epsH=1,eps2=1;

struct FastIn{
    vector<char> b; char *p;
    FastIn(){char tmp[1<<16]; size_t n; while((n=fread(tmp,1,sizeof(tmp),stdin))>0) b.insert(b.end(),tmp,tmp+n); b.push_back(0); p=b.data();}
    void skip(){while(*p==' '||*p=='\n'||*p=='\r'||*p=='\t')++p;}
    int nextInt(){skip(); int s=1; if(*p=='-')s=-1,++p; int x=0; while(*p>='0'&&*p<='9') x=x*10+(*p++-'0'); return x*s;}
    double nextDouble(){skip(); char*e; double x=strtod(p,&e); p=e; return x;}
    char nextChar(){skip(); return *p++;}
};

bool readInput(){
    FastIn in; N=in.nextInt(); M=in.nextInt(); if(N<=0||M<=0) return false;
    inmesh.v.resize(N); inmesh.f.resize(M);
    bb0={1e100,1e100,1e100}; bb1={-1e100,-1e100,-1e100}; center0={0,0,0};
    for(int i=0;i<N;i++){
        char c=in.nextChar(); if(c!='v'&&c!='V') --in.p;
        Vec3 p{in.nextDouble(),in.nextDouble(),in.nextDouble()}; inmesh.v[i]=p; center0=center0+p;
        bb0.x=min(bb0.x,p.x); bb0.y=min(bb0.y,p.y); bb0.z=min(bb0.z,p.z);
        bb1.x=max(bb1.x,p.x); bb1.y=max(bb1.y,p.y); bb1.z=max(bb1.z,p.z);
    }
    center0=center0/(double)N;
    for(int i=0;i<M;i++){
        char c=in.nextChar(); if(c!='f'&&c!='F') --in.p;
        int a=in.nextInt()-1,b=in.nextInt()-1,c2=in.nextInt()-1; inmesh.f[i]={a,b,c2};
    }
    diagLen=normv(bb1-bb0); if(!(diagLen>0)) diagLen=1; epsH=0.05*diagLen*0.999; eps2=epsH*epsH;
    return true;
}

void printMesh(const Mesh&m){
    printf("%d %d\n",(int)m.v.size(),(int)m.f.size());
    for(auto &p:m.v) printf("v %.12g %.12g %.12g\n",p.x,p.y,p.z);
    for(auto &f:m.f) printf("f %d %d %d\n",f.a+1,f.b+1,f.c+1);
}

void orientOutward(Mesh&m, Vec3 cen){
    for(auto &f:m.f){
        Vec3 cr=crossv(m.v[f.b]-m.v[f.a],m.v[f.c]-m.v[f.a]);
        Vec3 cc=(m.v[f.a]+m.v[f.b]+m.v[f.c])/3.0;
        if(dotv(cr,cc-cen)<0) swap(f.b,f.c);
    }
}

bool manifoldOK(const Mesh&m){
    if(m.v.size()<4 || m.f.size()<4 || m.v.size()>(size_t)N) return false;
    vector<unsigned long long> edges; edges.reserve(m.f.size()*3); set<array<int,3>> tris; vector<char> used(m.v.size(),0);
    double minArea=max(1e-30,1e-24*diagLen*diagLen);
    for(auto &f:m.f){
        if(f.a<0||f.b<0||f.c<0||f.a>=(int)m.v.size()||f.b>=(int)m.v.size()||f.c>=(int)m.v.size())return false;
        if(f.a==f.b||f.a==f.c||f.b==f.c)return false;
        if(norm2(crossv(m.v[f.b]-m.v[f.a],m.v[f.c]-m.v[f.a]))<=minArea)return false;
        if(!tris.insert(tkey(f.a,f.b,f.c)).second)return false;
        used[f.a]=used[f.b]=used[f.c]=1; edges.push_back(ekey(f.a,f.b)); edges.push_back(ekey(f.b,f.c)); edges.push_back(ekey(f.c,f.a));
    }
    for(char c:used) if(!c) return false;
    sort(edges.begin(),edges.end());
    for(size_t i=0;i<edges.size();){size_t j=i+1; while(j<edges.size()&&edges[j]==edges[i])++j; if(j-i!=2)return false; i=j;}
    return true;
}

static inline long long cellKey(int x,int y,int z){
    long long X=(long long)(x+1048576),Y=(long long)(y+1048576),Z=(long long)(z+1048576);
    return (X<<42) ^ (Y<<21) ^ Z;
}
struct RadiusGrid{
    double h,lim2; const vector<Vec3>*pts; unordered_map<long long, vector<int>> mp;
    void build(const vector<Vec3>&p,double r){pts=&p; h=max(r,1e-12); lim2=r*r*1.000004; mp.reserve(p.size()*2+7); for(int i=0;i<(int)p.size();++i){int x=floor(p[i].x/h),y=floor(p[i].y/h),z=floor(p[i].z/h);mp[cellKey(x,y,z)].push_back(i);}}
    bool near(Vec3 q)const{int x=floor(q.x/h),y=floor(q.y/h),z=floor(q.z/h); for(int dx=-1;dx<=1;dx++)for(int dy=-1;dy<=1;dy++)for(int dz=-1;dz<=1;dz++){auto it=mp.find(cellKey(x+dx,y+dy,z+dz)); if(it==mp.end())continue; for(int id:it->second) if(norm2((*pts)[id]-q)<=lim2) return true;} return false;}
};
bool hausOK(const Mesh&m){
    if(m.v.empty()||m.v.size()>(size_t)N)return false;
    RadiusGrid gs; gs.build(m.v,epsH); for(auto &p:inmesh.v) if(!gs.near(p)) return false;
    RadiusGrid go; go.build(inmesh.v,epsH); for(auto &p:m.v) if(!go.near(p)) return false;
    return true;
}

struct Maps{int R; vector<float> d,n;};
static void project(const Vec3&p,int view,int R,double&u,double&v,double&z){
    const double D=2.5, F=800.0*R/1024.0, C=0.5*R; double sx,sy;
    if(view==0){sx=p.y;sy=p.z;z=D-p.x;} else if(view==1){sx=-p.y;sy=p.z;z=D+p.x;}
    else if(view==2){sx=-p.x;sy=p.z;z=D-p.y;} else if(view==3){sx=p.x;sy=p.z;z=D+p.y;}
    else if(view==4){sx=p.x;sy=p.y;z=D-p.z;} else {sx=-p.x;sy=p.y;z=D+p.z;}
    u=F*sx/z+C; v=F*sy/z+C;
}
Maps renderMesh(const Mesh&m,int R){
    Maps mp; mp.R=R; int P=R*R; mp.d.assign((size_t)6*P,255.f); mp.n.assign((size_t)6*P*3,127.5f);
    vector<float> U(m.v.size()),Vv(m.v.size()),Z(m.v.size()); vector<Vec3> FN(m.f.size());
    for(int i=0;i<(int)m.f.size();++i){auto f=m.f[i];FN[i]=unit(crossv(m.v[f.b]-m.v[f.a],m.v[f.c]-m.v[f.a]));}
    for(int view=0;view<6;++view){
        for(int i=0;i<(int)m.v.size();++i){double u,v,z;project(m.v[i],view,R,u,v,z);U[i]=u;Vv[i]=v;Z[i]=z;}
        float*D=mp.d.data()+(size_t)view*P; float*Nn=mp.n.data()+(size_t)view*P*3;
        for(int ti=0;ti<(int)m.f.size();++ti){auto f=m.f[ti];int a=f.a,b=f.b,c=f.c;float z0=Z[a],z1=Z[b],z2=Z[c]; if(!(z0>0&&z1>0&&z2>0))continue;float u0=U[a],u1=U[b],u2=U[c],v0=Vv[a],v1=Vv[b],v2=Vv[c];
            int xmin=max(0,(int)floor(min({u0,u1,u2})-.5)),xmax=min(R-1,(int)ceil(max({u0,u1,u2})+.5));
            int ymin=max(0,(int)floor(min({v0,v1,v2})-.5)),ymax=min(R-1,(int)ceil(max({v0,v1,v2})+.5));
            if(xmin>xmax||ymin>ymax)continue; float den=(v1-v2)*(u0-u2)+(u2-u1)*(v0-v2); if(fabs(den)<1e-20)continue; float inv=1.f/den;
            Vec3 nn=FN[ti]; float nr=(nn.x+1)*127.5f,ng=(nn.y+1)*127.5f,nb=(nn.z+1)*127.5f;
            for(int y=ymin;y<=ymax;++y){float py=y+.5f;int row=y*R;for(int x=xmin;x<=xmax;++x){float px=x+.5f;float w0=((v1-v2)*(px-u2)+(u2-u1)*(py-v2))*inv;float w1=((v2-v0)*(px-u2)+(u0-u2)*(py-v2))*inv;float w2=1-w0-w1; if(w0>=-1e-6&&w1>=-1e-6&&w2>=-1e-6){float dep=1.f/(w0/z0+w1/z1+w2/z2);int id=row+x;if(dep<D[id]){D[id]=dep;float*q=Nn+3*id;q[0]=nr;q[1]=ng;q[2]=nb;}}}}
        }
    }
    return mp;
}
static inline double rect(const vector<double>&I,int W,int x0,int y0,int x1,int y1){return I[(size_t)y1*W+x1]-I[(size_t)y0*W+x1]-I[(size_t)y1*W+x0]+I[(size_t)y0*W+x0];}
double ssimChan(const float*A,int sa,const float*B,int sb,const float*DA,const float*DB,int R){
    int W=R+1; vector<double> IA((size_t)W*W),IB(IA),IA2(IA),IB2(IA),IAB(IA);
    for(int y=1;y<=R;++y){double a=0,b=0,a2=0,b2=0,ab=0;int row=(y-1)*R;for(int x=1;x<=R;++x){int p=row+x-1;double va=A[(size_t)p*sa],vb=B[(size_t)p*sb];a+=va;b+=vb;a2+=va*va;b2+=vb*vb;ab+=va*vb;size_t id=(size_t)y*W+x,up=(size_t)(y-1)*W+x;IA[id]=IA[up]+a;IB[id]=IB[up]+b;IA2[id]=IA2[up]+a2;IB2[id]=IB2[up]+b2;IAB[id]=IAB[up]+ab;}}
    int rad=5,area=121; const double c1=6.5025,c2=58.5225; long long cnt=0; long double acc=0;
    for(int y=rad;y<R-rad;++y){int row=y*R;for(int x=rad;x<R-rad;++x){int p=row+x;if(!(DA[p]<254.f||DB[p]<254.f))continue;int x0=x-rad,y0=y-rad,x1=x+rad+1,y1=y+rad+1;double sa0=rect(IA,W,x0,y0,x1,y1),sb0=rect(IB,W,x0,y0,x1,y1),sa2=rect(IA2,W,x0,y0,x1,y1),sb2=rect(IB2,W,x0,y0,x1,y1),sab=rect(IAB,W,x0,y0,x1,y1);double ma=sa0/area,mb=sb0/area,va=max(0.0,sa2/area-ma*ma),vb=max(0.0,sb2/area-mb*mb),cov=sab/area-ma*mb,den=(ma*ma+mb*mb+c1)*(va+vb+c2);acc+=den?((2*ma*mb+c1)*(2*cov+c2)/den):1;cnt++;}}
    return cnt?(double)(acc/cnt):1.0;
}
double approxSSIM(const Maps&a,const Maps&b){int R=a.R,P=R*R;double tot=0;for(int view=0;view<6;++view){const float*ad=a.d.data()+(size_t)view*P;const float*bd=b.d.data()+(size_t)view*P;double ns=0;for(int c=0;c<3;++c)ns+=ssimChan(a.n.data()+(size_t)view*P*3+c,3,b.n.data()+(size_t)view*P*3+c,3,ad,bd,R);ns/=3;double ds=ssimChan(ad,1,bd,1,ad,bd,R);tot+=0.5*(ns+ds);}return tot/6;}
bool acceptCandidate(Mesh&m,const Maps&origMaps,int R,double th){if(m.v.size()>=inmesh.v.size())return false;if(!manifoldOK(m))return false;if(!hausOK(m))return false;Maps rm=renderMesh(m,R);double s=approxSSIM(origMaps,rm);return s>=th;}

bool detectLatLong(int&rings,int&sides){
    if(N<100||M!=2*(N-2))return false;
    for(int s=8;s<=4096;++s){if((N-2)%s)continue;int r=(N-2)/s;if(r<3)continue;bool ok=true;int step=max(1,s/64);for(int j=0;j<s&&ok;j+=step){int a=1+j,b=1+(j+1)%s;if(!same_tri(inmesh.f[j],0,b,a))ok=false;}int span=max(1,(r-1)*s/256);for(int q=0;q<(r-1)*s&&ok;q+=span){int rr=q/s,j=q-rr*s;int a=1+rr*s+j,b=1+rr*s+(j+1)%s,c=1+(rr+1)*s+j,d=1+(rr+1)*s+(j+1)%s;int fid=s+2*q;if(fid+1>=M||!same_tri(inmesh.f[fid],a,b,c)||!same_tri(inmesh.f[fid+1],b,d,c))ok=false;}if(ok){rings=r;sides=s;return true;}}
    return false;
}
bool latLongRemesh(Mesh&out,const Maps&origMaps,int R){
    int rings,sides;if(!detectLatLong(rings,sides))return false;
    vector<pair<int,int>> trials; for(int d:{16,12,10,8,6,5,4,3}) trials.push_back({max(3,rings/d),max(8,sides/d)});
    for(auto [rr,ss]:trials){if(2+rr*ss>=N)continue;Mesh m;m.v.push_back(inmesh.v[0]);for(int i=0;i<rr;i++){int oi=1+(long long)i*(rings-1)/max(1,rr-1);for(int j=0;j<ss;j++){int oj=(long long)j*sides/ss;m.v.push_back(inmesh.v[1+(oi-1)*sides+oj]);}}int bot=m.v.size();m.v.push_back(inmesh.v[N-1]);auto id=[&](int r,int j){return 1+r*ss+(j%ss+ss)%ss;};for(int j=0;j<ss;j++)m.f.push_back({0,id(0,j),id(0,j+1)});for(int i=0;i<rr-1;i++)for(int j=0;j<ss;j++){int a=id(i,j),b=id(i,j+1),c=id(i+1,j),d=id(i+1,j+1);m.f.push_back({a,b,c});m.f.push_back({b,d,c});}for(int j=0;j<ss;j++)m.f.push_back({bot,id(rr-1,j+1),id(rr-1,j)});orientOutward(m,center0);if(acceptCandidate(m,origMaps,R,0.915)){out=m;return true;}}
    return false;
}

bool detectTorusGrid(int&U,int&Vn){
    if(N<100||M!=2*N)return false;
    for(int v=6;v<=4096;++v){if(N%v)continue;int u=N/v;if(u<6)continue;int step=max(1,N/650),tot=0,ok=0;for(int q=0;q<N&&tot<650;q+=step){int i=q/v,j=q%v;int a=i*v+j,b=i*v+(j+1)%v,c=((i+1)%u)*v+j,d=((i+1)%u)*v+(j+1)%v;int fid=2*q;if(fid+1>=M)break;bool good=(same_tri(inmesh.f[fid],a,b,c)&&same_tri(inmesh.f[fid+1],b,d,c))||(same_tri(inmesh.f[fid],a,c,b)&&same_tri(inmesh.f[fid+1],b,c,d))||(same_tri(inmesh.f[fid],a,b,d)&&same_tri(inmesh.f[fid+1],a,d,c));ok+=good;tot++;}if(tot>100&&ok*100>=tot*96){U=u;Vn=v;return true;}}
    return false;
}
bool torusGridRemesh(Mesh&out,const Maps&origMaps,int R){
    int U,Vn;if(!detectTorusGrid(U,Vn))return false;
    for(int d:{18,14,11,9,7,6,5,4,3}){int u=max(6,U/d),v=max(6,Vn/d);if(u*v>=N)continue;Mesh m;for(int i=0;i<u;i++){int oi=(long long)i*U/u;for(int j=0;j<v;j++){int oj=(long long)j*Vn/v;m.v.push_back(inmesh.v[oi*Vn+oj]);}}auto id=[&](int i,int j){return ((i%u+u)%u)*v+(j%v+v)%v;};for(int i=0;i<u;i++)for(int j=0;j<v;j++){int a=id(i,j),b=id(i,j+1),c=id(i+1,j),dd=id(i+1,j+1);m.f.push_back({a,b,c});m.f.push_back({b,dd,c});}orientOutward(m,center0);if(acceptCandidate(m,origMaps,R,0.915)){out=m;return true;}}
    return false;
}

struct DirKey{int x,y,z;bool operator<(const DirKey&o)const{return x!=o.x?x<o.x:(y!=o.y?y<o.y:z<o.z);} };
static DirKey keyFace(int face,int i,int j,int n){int a=-n+2*i,b=-n+2*j;if(face==0)return { n,a,b};if(face==1)return {-n,a,b};if(face==2)return {a, n,b};if(face==3)return {a,-n,b};if(face==4)return {a,b, n};return {a,b,-n};}
static Vec3 dirFromKey(const DirKey&k){return unit(Vec3{(double)k.x,(double)k.y,(double)k.z});}
static DirKey mapPointToKey(Vec3 q,int n){double ax=fabs(q.x),ay=fabs(q.y),az=fabs(q.z);auto quant=[&](double t){int q=(int)llround((t+1)*0.5*n);q=max(0,min(n,q));return -n+2*q;};if(ax>=ay&&ax>=az){double s=max(ax,1e-300);return {q.x>=0?n:-n,quant(q.y/s),quant(q.z/s)};}if(ay>=ax&&ay>=az){double s=max(ay,1e-300);return {quant(q.x/s),q.y>=0?n:-n,quant(q.z/s)};}double s=max(az,1e-300);return {quant(q.x/s),quant(q.y/s),q.z>=0?n:-n};}
bool cubeStarRemesh(Mesh&out,const Maps&origMaps,int R){
    if(N<1200)return false; Vec3 cen=center0; double bdiag=diagLen;
    vector<int> ns={6,8,10,12,14,16,18,20,24,28,32,36,44,52};
    for(int n:ns){
        map<DirKey,int> id; vector<DirKey> keys; vector<Vec3> dirs;
        auto getId=[&](DirKey k)->int{auto it=id.find(k);if(it!=id.end())return it->second;int r=id.size();id[k]=r;keys.push_back(k);dirs.push_back(dirFromKey(k));return r;};
        Mesh m;
        for(int face=0;face<6;face++)for(int i=0;i<n;i++)for(int j=0;j<n;j++){int a=getId(keyFace(face,i,j,n)),b=getId(keyFace(face,i+1,j,n)),c=getId(keyFace(face,i+1,j+1,n)),d=getId(keyFace(face,i,j+1,n));m.f.push_back({a,b,c});m.f.push_back({a,c,d});}
        int K=keys.size(); if(K>=N||K>30000)continue; vector<int> pick(K,-1); vector<double> best(K,-1e300);
        for(int i=0;i<N;i++){Vec3 q=inmesh.v[i]-cen; double l=normv(q); if(l<=1e-12)continue; DirKey kk=mapPointToKey(q,n); auto it=id.find(kk); if(it==id.end())continue; int h=it->second; double sc=dotv(q,dirs[h])+0.02*l; if(sc>best[h]){best[h]=sc;pick[h]=i;}}
        int empty=0; for(int i=0;i<K;i++)empty+=pick[i]<0; if(empty>K/3)continue;
        for(int h=0;h<K;h++) if(pick[h]<0){double bs=-1e300;int bi=-1;for(int i=0;i<N;i+=max(1,N/250000)){Vec3 q=inmesh.v[i]-cen;double l=normv(q);if(l<=1e-12)continue;double sc=dotv(q/l,dirs[h])+0.003*l/bdiag;if(sc>bs){bs=sc;bi=i;}}pick[h]=bi;}
        vector<int> tmp=pick; sort(tmp.begin(),tmp.end()); if(unique(tmp.begin(),tmp.end())!=tmp.end())continue;
        m.v.resize(K); for(int i=0;i<K;i++)m.v[i]=inmesh.v[pick[i]]; orientOutward(m,cen);
        if(acceptCandidate(m,origMaps,R,0.93)){out=m;return true;}
    }
    return false;
}

bool ellipsoidRemesh(Mesh&out,const Maps&origMaps,int R){
    if(N<900)return false; Vec3 c=(bb0+bb1)*0.5, rr=(bb1-bb0)*0.5; double rmax=max({rr.x,rr.y,rr.z}); if(rmax<=1e-12||min({rr.x,rr.y,rr.z})<0.15*rmax)return false;
    int stride=max(1,N/240000),cnt=0; double rms=0,ma=0; for(int i=0;i<N;i+=stride){Vec3 q=inmesh.v[i]-c; double s=sqrt(q.x*q.x/(rr.x*rr.x)+q.y*q.y/(rr.y*rr.y)+q.z*q.z/(rr.z*rr.z)); double e=fabs(s-1); rms+=e*e; ma=max(ma,e); cnt++;} rms=sqrt(rms/max(1,cnt)); if(rms>0.010||ma>0.050)return false;
    for(auto [lat,lon]:vector<pair<int,int>>{{12,24},{16,32},{20,40},{24,48},{28,56},{32,64},{40,80}}){if(2+(lat-1)*lon>=N)continue; Mesh m; m.v.push_back({c.x,c.y,c.z+rr.z}); for(int i=1;i<lat;i++){double th=M_PI*i/lat;for(int j=0;j<lon;j++){double ph=2*M_PI*j/lon;m.v.push_back({c.x+rr.x*sin(th)*cos(ph),c.y+rr.y*sin(th)*sin(ph),c.z+rr.z*cos(th)});}}int bot=m.v.size();m.v.push_back({c.x,c.y,c.z-rr.z});auto id=[&](int ring,int j){return 1+(ring-1)*lon+(j%lon+lon)%lon;};for(int j=0;j<lon;j++)m.f.push_back({0,id(1,j+1),id(1,j)});for(int i=1;i<lat-1;i++)for(int j=0;j<lon;j++){int a=id(i,j),b=id(i,j+1),cc=id(i+1,j),d=id(i+1,j+1);m.f.push_back({a,b,cc});m.f.push_back({b,d,cc});}for(int j=0;j<lon;j++)m.f.push_back({bot,id(lat-1,j),id(lat-1,j+1)});orientOutward(m,c); if(acceptCandidate(m,origMaps,R,0.925)){out=m;return true;}}
    return false;
}

struct Decimator{
    struct DF{int a,b,c;bool on;Vec3 n;};
    struct Edge{double c;int u,v,vu,vv;bool operator<(const Edge&o)const{return c>o.c;}};
    vector<Vec3>P; vector<DF>F; vector<vector<int>>inc; vector<char>alive; vector<int>ver,mark,head,tail,nxt,sz; vector<double>rad; int stamp=1,aliveN;
    bool has(int fid,int v){auto &f=F[fid];return f.a==v||f.b==v||f.c==v;}
    Vec3 fn(const DF&f){return unit(crossv(P[f.b]-P[f.a],P[f.c]-P[f.a]));}
    void clean(int v){if(!alive[v])return;auto &L=inc[v];int w=0;for(int id:L)if(F[id].on&&has(id,v))L[w++]=id;L.resize(w);} 
    Decimator(){P=inmesh.v;F.resize(M);for(int i=0;i<M;i++){F[i]={inmesh.f[i].a,inmesh.f[i].b,inmesh.f[i].c,true,{}};}aliveN=N;inc.assign(N,{});vector<int>d(N);for(auto &f:F){d[f.a]++;d[f.b]++;d[f.c]++;}for(int i=0;i<N;i++)inc[i].reserve(d[i]+8);for(int i=0;i<M;i++){inc[F[i].a].push_back(i);inc[F[i].b].push_back(i);inc[F[i].c].push_back(i);F[i].n=fn(F[i]);}alive.assign(N,1);ver.assign(N,0);mark.assign(N,0);rad.assign(N,0);head.resize(N);tail.resize(N);nxt.assign(N,-1);sz.assign(N,1);iota(head.begin(),head.end(),0);iota(tail.begin(),tail.end(),0);} 
    bool edgeFaces(int a,int b,int ids[2],int op[2]){int c=0;if(inc[a].size()>inc[b].size())swap(a,b);for(int id:inc[a])if(F[id].on&&has(id,a)&&has(id,b)){if(c>=2)return false;ids[c]=id;auto &f=F[id];op[c]=-1;if(f.a!=a&&f.a!=b)op[c]=f.a;if(f.b!=a&&f.b!=b)op[c]=f.b;if(f.c!=a&&f.c!=b)op[c]=f.c;c++;}return c==2&&op[0]>=0&&op[1]>=0&&op[0]!=op[1];}
    bool link(int a,int b){int ids[2],op[2];if(!edgeFaces(a,b,ids,op))return false;if(++stamp>1000000000){fill(mark.begin(),mark.end(),0);stamp=1;}int s1=stamp++,s2=stamp++;for(int id:inc[a])if(F[id].on&&has(id,a)){int vs[3]={F[id].a,F[id].b,F[id].c};for(int x:vs)if(x!=a&&x!=b)mark[x]=s1;}int cnt=0,g0=0,g1=0;for(int id:inc[b])if(F[id].on&&has(id,b)){int vs[3]={F[id].a,F[id].b,F[id].c};for(int x:vs)if(x!=a&&x!=b&&mark[x]==s1){mark[x]=s2;cnt++;if(x==op[0])g0=1;if(x==op[1])g1=1;if(cnt>2)return false;}}return cnt==2&&g0&&g1;}
    bool canMoveCluster(int kill,int keep){for(int x=head[kill];x!=-1;x=nxt[x])if(norm2(inmesh.v[x]-P[keep])>eps2)return false;return true;}
    bool normalOK(int keep,int kill){vector<int>ids=inc[keep];ids.insert(ids.end(),inc[kill].begin(),inc[kill].end());sort(ids.begin(),ids.end());ids.erase(unique(ids.begin(),ids.end()),ids.end());double coslim=cos(42.0*M_PI/180.0);double mina=max(1e-30,1e-24*diagLen*diagLen);for(int id:ids)if(F[id].on){auto f=F[id];bool hk=has(id,keep),hr=has(id,kill);if(hk&&hr)continue;int a=f.a,b=f.b,c=f.c;if(a==kill)a=keep;if(b==kill)b=keep;if(c==kill)c=keep;Vec3 cr=crossv(P[b]-P[a],P[c]-P[a]);if(norm2(cr)<=mina)return false;if(dotv(unit(cr),f.n)<coslim)return false;}return true;}
    bool duplicateAfter(int keep,int kill,int s0,int s1){set<array<int,3>> seen;for(int i=0;i<(int)F.size();i++)if(F[i].on&&i!=s0&&i!=s1){int a=F[i].a,b=F[i].b,c=F[i].c;if(a==kill)a=keep;if(b==kill)b=keep;if(c==kill)c=keep;if(a==b||a==c||b==c)return true;auto k=tkey(a,b,c);if(!seen.insert(k).second)return true;}return false;}
    bool collapse(int a,int b,priority_queue<Edge>&pq){int ids[2],op[2];if(!edgeFaces(a,b,ids,op)||!link(a,b))return false;int keep=a,kill=b;if(sz[a]<sz[b]){keep=b;kill=a;} if(!canMoveCluster(kill,keep))return false;if(!normalOK(keep,kill))return false;if(duplicateAfter(keep,kill,ids[0],ids[1]))return false;for(int id:ids)F[id].on=false;vector<int> L=inc[kill];for(int id:L)if(F[id].on&&has(id,kill)){if(F[id].a==kill)F[id].a=keep;if(F[id].b==kill)F[id].b=keep;if(F[id].c==kill)F[id].c=keep;inc[keep].push_back(id);}alive[kill]=0;aliveN--;ver[keep]++;ver[kill]++;nxt[tail[keep]]=head[kill];tail[keep]=tail[kill];sz[keep]+=sz[kill];inc[kill].clear();clean(keep);for(int id:inc[keep])if(F[id].on)F[id].n=fn(F[id]);if(++stamp>1000000000){fill(mark.begin(),mark.end(),0);stamp=1;}int st=stamp++;for(int id:inc[keep])if(F[id].on){int vs[3]={F[id].a,F[id].b,F[id].c};for(int x:vs)if(x!=keep&&alive[x]&&mark[x]!=st){mark[x]=st;push(keep,x,pq);}}return true;}
    void push(int a,int b,priority_queue<Edge>&pq){if(a==b||!alive[a]||!alive[b])return;double c=norm2(P[a]-P[b])+1e-9*(inc[a].size()+inc[b].size());pq.push({c,a,b,ver[a],ver[b]});}
    Mesh snap(){Mesh m;vector<int>id(N,-1);for(int i=0;i<N;i++)if(alive[i]){id[i]=m.v.size();m.v.push_back(P[i]);}for(auto &f:F)if(f.on){int a=id[f.a],b=id[f.b],c=id[f.c];if(a>=0&&b>=0&&c>=0&&a!=b&&a!=c&&b!=c)m.f.push_back({a,b,c});}return m;}
    vector<Mesh> run(){priority_queue<Edge>pq;vector<unsigned long long>ed;for(auto &f:F){ed.push_back(ekey(f.a,f.b));ed.push_back(ekey(f.b,f.c));ed.push_back(ekey(f.c,f.a));}sort(ed.begin(),ed.end());ed.erase(unique(ed.begin(),ed.end()),ed.end());for(auto k:ed)push(k>>32,k&0xffffffffu,pq);vector<int>targets;for(double r:{.55,.42,.32,.25,.20,.16,.13,.11,.095,.083}){int t=max(4,(int)ceil(N*r));if(t<N)targets.push_back(t);}int ti=0;vector<Mesh>res;int finalT=targets.empty()?max(4,N/8):targets.back();while(aliveN>finalT&&!pq.empty()&&clock()<CLOCKS_PER_SEC*18){auto e=pq.top();pq.pop();if(e.u<0||e.v<0||e.u>=N||e.v>=N||!alive[e.u]||!alive[e.v]||ver[e.u]!=e.vu||ver[e.v]!=e.vv)continue;if(collapse(e.u,e.v,pq)){while(ti<(int)targets.size()&&aliveN<=targets[ti]){res.push_back(snap());ti++;}}}if(res.empty())res.push_back(snap());return res;}
};

bool endpointDecimator(Mesh&out,const Maps&origMaps,int R){Decimator d;auto snaps=d.run();sort(snaps.begin(),snaps.end(),[](const Mesh&a,const Mesh&b){return a.v.size()<b.v.size();});for(auto &m:snaps){if(m.v.size()>=inmesh.v.size())continue;if(!manifoldOK(m))continue;Maps rm=renderMesh(m,R);if(approxSSIM(origMaps,rm)>=0.912){out=m;return true;}}return false;}

int main(){
    if(!readInput())return 0;
    if(N==8&&M==12){
        printf("8 12\n");
        for(auto &p:inmesh.v)printf("v %.9f %.9f %.9f\n",p.x,p.y,p.z);
        for(auto &f:inmesh.f)printf("f %d %d %d\n",f.a+1,f.b+1,f.c+1);
        return 0;
    }
    int R = N<25000?192:(N<150000?128:96);
    Maps origMaps=renderMesh(inmesh,R);
    Mesh best=inmesh,cand; bool have=false;
    auto take=[&](Mesh&m){if(m.v.size()<best.v.size()){best=m;have=true;}};
    if(latLongRemesh(cand,origMaps,R))take(cand);
    if(torusGridRemesh(cand,origMaps,R))take(cand);
    if(cubeStarRemesh(cand,origMaps,R))take(cand);
    if(ellipsoidRemesh(cand,origMaps,R))take(cand);
    if(endpointDecimator(cand,origMaps,R))take(cand);
    if(!have||!manifoldOK(best))best=inmesh;
    printMesh(best);
    return 0;
}