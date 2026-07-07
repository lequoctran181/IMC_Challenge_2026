#include <bits/stdc++.h>
using namespace std;

struct Vec3{double x=0,y=0,z=0;};
static inline Vec3 operator+(const Vec3&a,const Vec3&b){return {a.x+b.x,a.y+b.y,a.z+b.z};}
static inline Vec3 operator-(const Vec3&a,const Vec3&b){return {a.x-b.x,a.y-b.y,a.z-b.z};}
static inline Vec3 operator*(const Vec3&a,double s){return {a.x*s,a.y*s,a.z*s};}
static inline Vec3 operator/(const Vec3&a,double s){return {a.x/s,a.y/s,a.z/s};}
static inline double dot3(const Vec3&a,const Vec3&b){return a.x*b.x+a.y*b.y+a.z*b.z;}
static inline Vec3 cross3(const Vec3&a,const Vec3&b){return {a.y*b.z-a.z*b.y,a.z*b.x-a.x*b.z,a.x*b.y-a.y*b.x};}
static inline double norm2(const Vec3&a){return dot3(a,a);} 
static inline double norm3(const Vec3&a){return sqrt(max(0.0,norm2(a)));}
static inline Vec3 unit3(Vec3 a){double n=norm3(a);return n>1e-300?a*(1.0/n):Vec3{0,0,0};}
static inline double clampd(double x,double a,double b){return x<a?a:(x>b?b:x);} 

struct Face{int v[3];};
struct Mesh{vector<Vec3> V;vector<Face> F;};

struct FastInput{
    vector<char> b; char *p=nullptr;
    FastInput(){b.reserve(1<<27);char t[1<<16];size_t n;while((n=fread(t,1,sizeof(t),stdin))>0)b.insert(b.end(),t,t+n);b.push_back(0);p=b.data();}
    inline void skip(){while(*p==' '||*p=='\n'||*p=='\r'||*p=='\t')++p;}
    int nextInt(){skip();int s=1;if(*p=='-'){s=-1;++p;}int x=0;while(*p>='0'&&*p<='9'){x=x*10+(*p-'0');++p;}return x*s;}
    double nextDouble(){skip();char*e;double x=strtod(p,&e);p=e;return x;}
    char nextChar(){skip();return *p?*p++:0;}
};

static Mesh ORG;
static int N0=0,M0=0;
static Vec3 BBmin,BBmax;
static double DIAG=1.0, EPS=0.0495, EPS2=0.0495*0.0495;
static chrono::steady_clock::time_point T0;
static inline double elapsed(){return chrono::duration<double>(chrono::steady_clock::now()-T0).count();}

static inline uint64_t ekey(int a,int b){if(a>b)swap(a,b);return (uint64_t)(uint32_t)a<<32 | (uint32_t)b;}
static inline array<int,3> triKey(int a,int b,int c){array<int,3> r{a,b,c};sort(r.begin(),r.end());return r;}
static inline bool sameTri(const Face&f,int a,int b,int c){return triKey(f.v[0],f.v[1],f.v[2])==triKey(a,b,c);} 
static inline Vec3 faceCross(const vector<Vec3>&V,const Face&f){return cross3(V[f.v[1]]-V[f.v[0]],V[f.v[2]]-V[f.v[0]]);} 

static bool readInput(){
    FastInput in;N0=in.nextInt();M0=in.nextInt();
    if(N0<=0||M0<=0)return false;
    ORG.V.resize(N0);ORG.F.resize(M0);
    BBmin={1e100,1e100,1e100};BBmax={-1e100,-1e100,-1e100};
    for(int i=0;i<N0;i++){
        char c=in.nextChar(); if(c!='v'&&c!='V') --in.p;
        Vec3 p{in.nextDouble(),in.nextDouble(),in.nextDouble()};
        ORG.V[i]=p;
        BBmin.x=min(BBmin.x,p.x);BBmin.y=min(BBmin.y,p.y);BBmin.z=min(BBmin.z,p.z);
        BBmax.x=max(BBmax.x,p.x);BBmax.y=max(BBmax.y,p.y);BBmax.z=max(BBmax.z,p.z);
    }
    for(int i=0;i<M0;i++){
        char c=in.nextChar(); if(c!='f'&&c!='F') --in.p;
        int a=in.nextInt()-1,b=in.nextInt()-1,c3=in.nextInt()-1;
        ORG.F[i]={{a,b,c3}};
    }
    DIAG=norm3(BBmax-BBmin); if(!(DIAG>0))DIAG=1.0;
    EPS=0.0495*DIAG; EPS2=EPS*EPS;
    return true;
}

static double vol6(const Mesh&m){long double s=0;for(const Face&f:m.F){const Vec3&a=m.V[f.v[0]],&b=m.V[f.v[1]],&c=m.V[f.v[2]];s+=dot3(a,cross3(b,c));}return (double)s;}
static void orientLikeOriginal(Mesh&m){double a=vol6(ORG),b=vol6(m);if(a*b<0)for(Face&f:m.F)swap(f.v[1],f.v[2]);}
static void outputMesh(const Mesh&m){
    static char outbuf[1<<20];setvbuf(stdout,outbuf,_IOFBF,sizeof(outbuf));
    printf("%d %d\n",(int)m.V.size(),(int)m.F.size());
    for(const Vec3&p:m.V)printf("v %.17g %.17g %.17g\n",p.x,p.y,p.z);
    for(const Face&f:m.F)printf("f %d %d %d\n",f.v[0]+1,f.v[1]+1,f.v[2]+1);
}
static void outputOriginal(){outputMesh(ORG);} 

static bool validateMesh(const Mesh&m){
    int n=(int)m.V.size(),mf=(int)m.F.size();
    if(n<=0||mf<=0||n>N0)return false;
    vector<unsigned char> used(n,0);
    vector<uint64_t> ed;ed.reserve((size_t)mf*3);
    vector<array<int,3>> tk;tk.reserve(mf);
    double minA2=max(1e-300,1e-26*DIAG*DIAG*DIAG*DIAG);
    for(const Face&f:m.F){
        int a=f.v[0],b=f.v[1],c=f.v[2];
        if(a<0||b<0||c<0||a>=n||b>=n||c>=n||a==b||a==c||b==c)return false;
        double a2=norm2(faceCross(m.V,f));
        if(!(a2>minA2)||!isfinite(a2))return false;
        used[a]=used[b]=used[c]=1;
        ed.push_back(ekey(a,b));ed.push_back(ekey(b,c));ed.push_back(ekey(c,a));
        tk.push_back(triKey(a,b,c));
    }
    for(int i=0;i<n;i++)if(!used[i])return false;
    sort(tk.begin(),tk.end()); if(adjacent_find(tk.begin(),tk.end())!=tk.end())return false;
    sort(ed.begin(),ed.end());
    for(size_t i=0;i<ed.size();){size_t j=i+1;while(j<ed.size()&&ed[j]==ed[i])++j;if(j-i!=2)return false;i=j;}
    return true;
}

struct PointGrid{
    Vec3 mn,mx;double cell=1;int nx=1,ny=1,nz=1;vector<vector<int>> bucket;const vector<Vec3>*P=nullptr;
    int ci(int a,int n)const{return a<0?0:(a>=n?n-1:a);} 
    int ix(double x)const{return ci((int)floor((x-mn.x)/cell),nx);}int iy(double y)const{return ci((int)floor((y-mn.y)/cell),ny);}int iz(double z)const{return ci((int)floor((z-mn.z)/cell),nz);}int key(int x,int y,int z)const{return (z*ny+y)*nx+x;}
    void build(const vector<Vec3>&pts,double c){
        P=&pts;cell=max(c,1e-12);mn={1e100,1e100,1e100};mx={-1e100,-1e100,-1e100};
        for(const Vec3&p:pts){mn.x=min(mn.x,p.x);mn.y=min(mn.y,p.y);mn.z=min(mn.z,p.z);mx.x=max(mx.x,p.x);mx.y=max(mx.y,p.y);mx.z=max(mx.z,p.z);} 
        mn=mn-Vec3{cell,cell,cell};mx=mx+Vec3{cell,cell,cell};
        nx=max(1,(int)floor((mx.x-mn.x)/cell)+1);ny=max(1,(int)floor((mx.y-mn.y)/cell)+1);nz=max(1,(int)floor((mx.z-mn.z)/cell)+1);
        long long tot=1LL*nx*ny*nz;
        if(tot>3000000LL){double side=max({mx.x-mn.x,mx.y-mn.y,mx.z-mn.z,cell});cell=side;nx=ny=nz=1;}
        bucket.assign((size_t)nx*ny*nz,{});
        for(int i=0;i<(int)pts.size();i++)bucket[key(ix(pts[i].x),iy(pts[i].y),iz(pts[i].z))].push_back(i);
    }
    bool near(const Vec3&q,double r2)const{
        int X=ix(q.x),Y=iy(q.y),Z=iz(q.z);
        int rad=max(1,(int)ceil(sqrt(r2)/cell));
        for(int dz=-rad;dz<=rad;dz++){int z=Z+dz;if(z<0||z>=nz)continue;
            for(int dy=-rad;dy<=rad;dy++){int y=Y+dy;if(y<0||y>=ny)continue;
                for(int dx=-rad;dx<=rad;dx++){int x=X+dx;if(x<0||x>=nx)continue;
                    for(int id:bucket[key(x,y,z)])if(norm2((*P)[id]-q)<=r2)return true;
                }} }
        return false;
    }
};
static bool coverageOK(const vector<Vec3>&cand,double scale=0.0495){
    if(cand.empty()||(int)cand.size()>N0)return false;
    double e=scale*DIAG,r2=e*e*(1.0+1e-11);
    PointGrid g;g.build(cand,e);
    for(int i=0;i<N0;i++){if((i&8191)==0&&elapsed()>19.25)return false;if(!g.near(ORG.V[i],r2))return false;}
    PointGrid go;go.build(ORG.V,e);
    for(int i=0;i<(int)cand.size();i++)if(!go.near(cand[i],r2))return false;
    return true;
}

struct RenderMap{int R=0;vector<float>d,n;};
static unordered_map<int,RenderMap> origRenderCache;
static inline void projectPoint(const Vec3&p,int view,int R,double&u,double&v,double&z){
    const double D=2.5;double f=800.0*(double)R/1024.0,c=0.5*(double)R,sx,sy;
    if(view==0){sx=p.y;sy=p.z;z=D-p.x;}else if(view==1){sx=-p.y;sy=p.z;z=D+p.x;}
    else if(view==2){sx=-p.x;sy=p.z;z=D-p.y;}else if(view==3){sx=p.x;sy=p.z;z=D+p.y;}
    else if(view==4){sx=p.x;sy=p.y;z=D-p.z;}else{sx=-p.x;sy=p.y;z=D+p.z;}
    u=f*sx/z+c;v=f*sy/z+c;
}
static RenderMap renderMesh(const Mesh&m,int R){
    RenderMap rm;rm.R=R;int P=R*R;rm.d.assign((size_t)6*P,255.0f);rm.n.assign((size_t)6*P*3,127.5f);
    vector<Vec3> fn(m.F.size());
    for(size_t i=0;i<m.F.size();i++)fn[i]=unit3(faceCross(m.V,m.F[i]));
    vector<float> U(m.V.size()),V(m.V.size()),Z(m.V.size());
    for(int view=0;view<6;view++){
        for(size_t i=0;i<m.V.size();i++){double u,v,z;projectPoint(m.V[i],view,R,u,v,z);U[i]=(float)u;V[i]=(float)v;Z[i]=(float)z;}
        float*db=rm.d.data()+(size_t)view*P;float*nb=rm.n.data()+(size_t)view*P*3;
        for(size_t ti=0;ti<m.F.size();ti++){
            const Face&f=m.F[ti];int ia=f.v[0],ib=f.v[1],ic=f.v[2];
            float x0=U[ia],x1=U[ib],x2=U[ic],y0=V[ia],y1=V[ib],y2=V[ic],z0=Z[ia],z1=Z[ib],z2=Z[ic];
            if(!(z0>0&&z1>0&&z2>0))continue;
            int xmin=max(0,(int)floor(min(x0,min(x1,x2))-0.5f)),xmax=min(R-1,(int)ceil(max(x0,max(x1,x2))+0.5f));
            int ymin=max(0,(int)floor(min(y0,min(y1,y2))-0.5f)),ymax=min(R-1,(int)ceil(max(y0,max(y1,y2))+0.5f));
            if(xmin>xmax||ymin>ymax)continue;
            float den=(y1-y2)*(x0-x2)+(x2-x1)*(y0-y2);if(fabs(den)<1e-20f)continue;float inv=1.0f/den;
            Vec3 nn=fn[ti];float nr=(float)((nn.x+1.0)*127.5),ng=(float)((nn.y+1.0)*127.5),nzv=(float)((nn.z+1.0)*127.5);
            for(int yy=ymin;yy<=ymax;yy++){float py=yy+0.5f;int row=yy*R;for(int xx=xmin;xx<=xmax;xx++){
                float px=xx+0.5f;float a=((y1-y2)*(px-x2)+(x2-x1)*(py-y2))*inv;float b=((y2-y0)*(px-x2)+(x0-x2)*(py-y2))*inv;float c=1.0f-a-b;
                if(a>=-1e-7f&&b>=-1e-7f&&c>=-1e-7f){float dep=1.0f/(a/z0+b/z1+c/z2);int id=row+xx;if(dep<db[id]){db[id]=dep;float*q=nb+(size_t)id*3;q[0]=nr;q[1]=ng;q[2]=nzv;}}
            }}
        }
    }
    return rm;
}
static const RenderMap& origRender(int R){auto it=origRenderCache.find(R);if(it!=origRenderCache.end())return it->second;auto p=origRenderCache.emplace(R,renderMesh(ORG,R));return p.first->second;}
static inline double rectSum(const vector<double>&I,int W,int x0,int y0,int x1,int y1){return I[(size_t)y1*W+x1]-I[(size_t)y0*W+x1]-I[(size_t)y1*W+x0]+I[(size_t)y0*W+x0];}
static double ssimChannel(const float*A,int sa,const float*B,int sb,const float*DA,const float*DB,int R,vector<double>&IA,vector<double>&IB,vector<double>&IA2,vector<double>&IB2,vector<double>&IAB){
    int W=R+1;size_t SZ=(size_t)W*W;fill(IA.begin(),IA.begin()+SZ,0);fill(IB.begin(),IB.begin()+SZ,0);fill(IA2.begin(),IA2.begin()+SZ,0);fill(IB2.begin(),IB2.begin()+SZ,0);fill(IAB.begin(),IAB.begin()+SZ,0);
    for(int y=1;y<=R;y++){double a=0,b=0,a2=0,b2=0,ab=0;int row=(y-1)*R;for(int x=1;x<=R;x++){int p=row+x-1;double av=A[(size_t)p*sa],bv=B[(size_t)p*sb];a+=av;b+=bv;a2+=av*av;b2+=bv*bv;ab+=av*bv;size_t id=(size_t)y*W+x,up=(size_t)(y-1)*W+x;IA[id]=IA[up]+a;IB[id]=IB[up]+b;IA2[id]=IA2[up]+a2;IB2[id]=IB2[up]+b2;IAB[id]=IAB[up]+ab;}}
    const int rad=5;const double area=(2*rad+1)*(2*rad+1),c1=6.5025,c2=58.5225;long long cnt=0;long double acc=0;
    for(int y=0;y<R;y++){int y0=max(0,y-rad),y1=min(R,y+rad+1),row=y*R;for(int x=0;x<R;x++){int p=row+x;if(!(DA[p]<254.0f||DB[p]<254.0f))continue;int x0=max(0,x-rad),x1=min(R,x+rad+1);double ar=(y1-y0)*(x1-x0);double sx=rectSum(IA,W,x0,y0,x1,y1),sy=rectSum(IB,W,x0,y0,x1,y1),sx2=rectSum(IA2,W,x0,y0,x1,y1),sy2=rectSum(IB2,W,x0,y0,x1,y1),sxy=rectSum(IAB,W,x0,y0,x1,y1);double mx=sx/ar,my=sy/ar,vx=max(0.0,sx2/ar-mx*mx),vy=max(0.0,sy2/ar-my*my),cov=sxy/ar-mx*my;double den=(mx*mx+my*my+c1)*(vx+vy+c2);acc+=den?((2*mx*my+c1)*(2*cov+c2)/den):1.0;++cnt;}}
    (void)area;return cnt?(double)(acc/cnt):1.0;
}
static double visualProxy(const Mesh&cand,int R){
    if(elapsed()>18.7)return 0.0;const RenderMap&a=origRender(R);RenderMap b=renderMesh(cand,R);int P=R*R,W=R+1;vector<double>IA((size_t)W*W),IB((size_t)W*W),IA2((size_t)W*W),IB2((size_t)W*W),IAB((size_t)W*W);double tot=0;
    for(int v=0;v<6;v++){const float*ad=a.d.data()+(size_t)v*P;const float*bd=b.d.data()+(size_t)v*P;double ns=0;for(int c=0;c<3;c++)ns+=ssimChannel(a.n.data()+(size_t)v*P*3+c,3,b.n.data()+(size_t)v*P*3+c,3,ad,bd,R,IA,IB,IA2,IB2,IAB);ns/=3.0;double ds=ssimChannel(ad,1,bd,1,ad,bd,R,IA,IB,IA2,IB2,IAB);tot+=0.5*(ns+ds);}return tot/6.0;
}
static int proxyRes(){if(N0<4000)return 384;if(N0<30000)return 224;if(N0<90000)return 160;return 112;}
static double proxyNeed(){if(N0<4000)return .935;if(N0<30000)return .925;if(N0<90000)return .918;return .912;}

static vector<Vec3> vertexNormals(){vector<Vec3> n(N0,{0,0,0});for(const Face&f:ORG.F){Vec3 cr=faceCross(ORG.V,f);n[f.v[0]]=n[f.v[0]]+cr;n[f.v[1]]=n[f.v[1]]+cr;n[f.v[2]]=n[f.v[2]]+cr;}return n;}
static void addOriented(vector<Face>&F,const vector<Vec3>&V,Face f,const Vec3&ref){Vec3 cr=faceCross(V,f);if(dot3(cr,ref)<0)swap(f.v[1],f.v[2]);F.push_back(f);} 
static bool adjSmall(const int a[3],int m,int&base){for(int t=0;t<3;t++)for(int s=0;s<2;s++){int x=(a[t]-s+m)%m;bool ok=true;for(int i=0;i<3;i++){int d=(a[i]-x+m)%m;if(d!=0&&d!=1){ok=false;break;}}if(ok){base=x;return true;}}return false;}
static bool tubeFaceCompatible(const Face&f,int S){if(S<8||N0%S)return false;int U=N0/S;if(U<8)return false;int r[3]={f.v[0]/S,f.v[1]/S,f.v[2]/S};int c[3]={f.v[0]%S,f.v[1]%S,f.v[2]%S};int rb=0,cb=0;if(!adjSmall(r,U,rb)||!adjSmall(c,S,cb))return false;int mask=0;for(int i=0;i<3;i++){int x=(r[i]-rb+U)%U,y=(c[i]-cb+S)%S;if(x>1||y>1)return false;mask|=1<<(x*2+y);}return __builtin_popcount((unsigned)mask)==3;}
static vector<int> tubeCandidates(){
    map<int,int> cnt;int st=max(1,M0/100000);
    for(int i=0;i<M0;i+=st){const Face&f=ORG.F[i];int a[3]={f.v[0],f.v[1],f.v[2]};for(int k=0;k<3;k++){int d=abs(a[k]-a[(k+1)%3]);if(!d)continue;d=min(d,N0-d);if(d>=8&&d<=N0/3)cnt[d]++;}}
    vector<pair<int,int>> q;for(auto&p:cnt)q.push_back({p.second,p.first});sort(q.rbegin(),q.rend());vector<int> r;
    auto add=[&](int s){if(s>=8&&s<=N0/4&&N0%s==0&&find(r.begin(),r.end(),s)==r.end())r.push_back(s);};
    for(int i=0;i<(int)q.size()&&i<20;i++){int d=q[i].second;for(int e=-4;e<=4;e++)add(d+e);if(d)add(N0/d);} 
    for(int s:{16,20,24,28,32,40,48,56,64,72,80,96,100,112,120,128,144,160,192,200,224,256,320,384})add(s);
    return r;
}
static bool goodTubeS(int S){if(N0%S)return false;int st=max(1,M0/120000),tot=0,ok=0;for(int i=0;i<M0;i+=st){tot++;if(tubeFaceCompatible(ORG.F[i],S))ok++;}return tot>50&&ok*1000>=tot*995;}
static bool makeTube(int S,int U2,int S2,Mesh&out,const vector<Vec3>&vn){
    if(N0%S)return false;int U=N0/S;if(U2<4||S2<4||U2>U||S2>S)return false;out.V.clear();out.F.clear();out.V.reserve(U2*S2);vector<int> src(U2*S2);
    for(int i=0;i<U2;i++){int oi=(long long)i*U/U2;for(int j=0;j<S2;j++){int oj=(long long)j*S/S2;int s=oi*S+oj;src[i*S2+j]=s;out.V.push_back(ORG.V[s]);}}
    auto id=[&](int i,int j){return ((i%U2+U2)%U2)*S2+((j%S2+S2)%S2);};
    out.F.reserve(2*U2*S2);for(int i=0;i<U2;i++)for(int j=0;j<S2;j++){int a=id(i,j),b=id(i+1,j),c=id(i+1,j+1),d=id(i,j+1);addOriented(out.F,out.V,Face{{a,b,d}},vn[src[a]]+vn[src[b]]+vn[src[d]]);addOriented(out.F,out.V,Face{{b,c,d}},vn[src[b]]+vn[src[c]]+vn[src[d]]);} 
    orientLikeOriginal(out);return true;
}
static bool tryTube(Mesh&best){
    if(M0!=2*N0||N0<2000||elapsed()>10.0)return false;vector<int> ss=tubeCandidates();int S=0;for(int s:ss)if(goodTubeS(s)){S=s;break;}if(!S)return false;int U=N0/S;vector<Vec3> vn=vertexNormals();
    vector<pair<int,int>> dims;auto add=[&](int u,int v){u=max(4,min(U,u));v=max(4,min(S,v));if(u*v<N0&&find(dims.begin(),dims.end(),make_pair(u,v))==dims.end())dims.push_back({u,v});};
    for(double q:{18.0,16.0,14.0,12.0,10.0,9.0,8.0,7.0,6.0,5.0,4.0,3.0,2.5,2.0})add((int)ceil(U/q),(int)ceil(S/q));
    for(int t:{768,1024,1536,2048,2560,3072,4096,5120,6144,8192,10240,12288,16384}){double ar=sqrt((double)U/max(1,S));int u=max(4,(int)(sqrt((double)t)*ar+0.5));int v=max(4,t/max(1,u));add(u,v);}
    sort(dims.begin(),dims.end(),[](auto&a,auto&b){return a.first*a.second<b.first*b.second;});
    Mesh cand;int R=proxyRes();double need=proxyNeed()+0.015;
    for(auto [u,v]:dims){if(elapsed()>18.0)break;if(u*v>=(int)best.V.size())continue;if(!makeTube(S,u,v,cand,vn))continue;if(!validateMesh(cand))continue;if(!coverageOK(cand.V,0.0492))continue;double q=visualProxy(cand,R);if(q>=need){best=cand;return true;}}
    return false;
}
static bool sameOrigFace(int fid,int a,int b,int c){return fid>=0&&fid<M0&&sameTri(ORG.F[fid],a,b,c);} 
static bool detectPolar(int&R,int&Vn){
    if(N0<300||M0!=2*(N0-2))return false;
    for(int v=8;v<=4096;v++){if((N0-2)%v)continue;int r=(N0-2)/v;if(r<3)continue;bool ok=true;int step=max(1,v/96);
        for(int j=0;j<v&&ok;j+=step){int a=1+j,b=1+(j+1)%v;if(!sameOrigFace(j,0,b,a))ok=false;int off=v+2*(r-1)*v+j,c=1+(r-1)*v+j,d=1+(r-1)*v+(j+1)%v;if(ok&&!sameOrigFace(off,N0-1,c,d))ok=false;}
        int span=max(1,(r-1)*v/400);for(int q=0;q<(r-1)*v&&ok;q+=span){int rr=q/v,j=q%v;int a=1+rr*v+j,b=1+rr*v+(j+1)%v,c=1+(rr+1)*v+j,d=1+(rr+1)*v+(j+1)%v;int f=v+2*(rr*v+j);if(!sameOrigFace(f,a,b,c)||!sameOrigFace(f+1,b,d,c))ok=false;}
        if(ok){R=r;Vn=v;return true;}
    }return false;
}
static bool makePolar(int R,int Vn,int R2,int V2,Mesh&out,const vector<Vec3>&vn){
    if(R2<2||V2<6||2+R2*V2>=N0)return false;out.V.clear();out.F.clear();vector<int> src;out.V.reserve(2+R2*V2);src.reserve(2+R2*V2);out.V.push_back(ORG.V[0]);src.push_back(0);
    for(int i=0;i<R2;i++){int oi=1+(int)((long long)i*(R-1)/max(1,R2-1));for(int j=0;j<V2;j++){int oj=(long long)j*Vn/V2;int s=1+(oi-1)*Vn+oj;src.push_back(s);out.V.push_back(ORG.V[s]);}}
    int bot=(int)out.V.size();out.V.push_back(ORG.V[N0-1]);src.push_back(N0-1);auto rid=[&](int r,int j){return 1+(r-1)*V2+((j%V2+V2)%V2);};auto add=[&](int a,int b,int c){addOriented(out.F,out.V,Face{{a,b,c}},vn[src[a]]+vn[src[b]]+vn[src[c]]);};
    for(int j=0;j<V2;j++)add(0,rid(1,j+1),rid(1,j));for(int r=1;r<R2;r++)for(int j=0;j<V2;j++){int a=rid(r,j),b=rid(r,j+1),c=rid(r+1,j),d=rid(r+1,j+1);add(a,b,c);add(b,d,c);}for(int j=0;j<V2;j++)add(bot,rid(R2,j),rid(R2,j+1));orientLikeOriginal(out);return true;
}
static bool tryPolar(Mesh&best){
    if(elapsed()>10.0)return false;int R=0,Vn=0;if(!detectPolar(R,Vn))return false;vector<Vec3> vn=vertexNormals();vector<pair<int,int>> dims;auto add=[&](int r,int v){r=max(2,min(R,r));v=max(6,min(Vn,v));if(2+r*v<N0&&find(dims.begin(),dims.end(),make_pair(r,v))==dims.end())dims.push_back({r,v});};
    for(double q:{20.0,18.0,16.0,14.0,12.0,10.0,9.0,8.0,7.0,6.0,5.0,4.0,3.0,2.5,2.0})add((int)ceil(R/q),(int)ceil(Vn/q));
    for(int t:{258,386,514,770,1026,1538,2050,3074,4098,6146,8194}){double ar=sqrt((double)R/max(1,Vn));int r=max(2,(int)(sqrt((double)t)*ar+0.5));int v=max(6,(t-2)/max(1,r));add(r,v);}sort(dims.begin(),dims.end(),[](auto&a,auto&b){return 2+a.first*a.second<2+b.first*b.second;});
    Mesh cand;int Rr=proxyRes();double need=proxyNeed()+0.018;for(auto [r,v]:dims){if(elapsed()>18.0)break;if(2+r*v>=(int)best.V.size())continue;if(!makePolar(R,Vn,r,v,cand,vn))continue;if(!validateMesh(cand))continue;if(!coverageOK(cand.V,0.0492))continue;double q=visualProxy(cand,Rr);if(q>=need){best=cand;return true;}}
    return false;
}

static bool tryCluster(Mesh&best){
    if(N0<1000||elapsed()>8.0)return false;
    vector<double> cells={0.044,0.040,0.036,0.032,0.028,0.024};
    int R=proxyRes();double need=proxyNeed()+0.025;
    for(double cs:cells){if(elapsed()>16.5)break;double cell=cs*DIAG;unordered_map<long long,int> mp;mp.reserve(N0*2);vector<Vec3> sum;vector<int> cnt;vector<int> mapv(N0,-1);
        auto hkey=[&](const Vec3&p){long long ix=(long long)floor((p.x-BBmin.x)/cell),iy=(long long)floor((p.y-BBmin.y)/cell),iz=(long long)floor((p.z-BBmin.z)/cell);return (ix*73856093LL)^(iy*19349663LL)^(iz*83492791LL);};
        for(int i=0;i<N0;i++){long long k=hkey(ORG.V[i]);auto it=mp.find(k);int id;if(it==mp.end()){id=(int)sum.size();mp.emplace(k,id);sum.push_back(ORG.V[i]);cnt.push_back(1);}else{id=it->second;sum[id]=sum[id]+ORG.V[i];cnt[id]++;}mapv[i]=id;}
        if((int)sum.size()>=(int)best.V.size()||sum.empty())continue;Mesh cand;cand.V.resize(sum.size());for(int i=0;i<(int)sum.size();i++)cand.V[i]=sum[i]/(double)cnt[i];
        vector<array<int,3>> seen;seen.reserve(M0);cand.F.reserve(M0);for(const Face&f:ORG.F){int a=mapv[f.v[0]],b=mapv[f.v[1]],c=mapv[f.v[2]];if(a==b||a==c||b==c)continue;array<int,3> tk=triKey(a,b,c);seen.push_back(tk);cand.F.push_back(Face{{a,b,c}});} 
        if(cand.F.empty())continue;sort(seen.begin(),seen.end());vector<Face> nf;nf.reserve(cand.F.size());unordered_set<unsigned long long> used;used.reserve(cand.F.size()*2+1);auto pack=[&](array<int,3> t)->unsigned long long{return ((unsigned long long)(uint32_t)t[0]<<42)^((unsigned long long)(uint32_t)t[1]<<21)^(uint32_t)t[2];};
        for(const Face&f:cand.F){array<int,3> t=triKey(f.v[0],f.v[1],f.v[2]);unsigned long long k=pack(t);if(used.insert(k).second)nf.push_back(f);}cand.F.swap(nf);
        vector<int> usedv(cand.V.size(),0);for(const Face&f:cand.F)usedv[f.v[0]]=usedv[f.v[1]]=usedv[f.v[2]]=1;vector<int> rem(cand.V.size(),-1);Mesh comp;for(int i=0;i<(int)cand.V.size();i++)if(usedv[i]){rem[i]=comp.V.size();comp.V.push_back(cand.V[i]);}for(const Face&f:cand.F)comp.F.push_back(Face{{rem[f.v[0]],rem[f.v[1]],rem[f.v[2]]}});orientLikeOriginal(comp);
        if((int)comp.V.size()>=(int)best.V.size())continue;if(!validateMesh(comp))continue;if(!coverageOK(comp.V,0.0490))continue;double q=visualProxy(comp,R);if(q>=need){best=comp;return true;}
    }
    return false;
}

static vector<Vec3> P;static vector<Face> F;static vector<unsigned char> aliveV,aliveF;static vector<double> coverR;static vector<vector<int>> adj;static int aliveFaces=0,aliveVerts=0;static vector<int> markA,markB;static int stampA=1,stampB=1;
static inline bool faceHas(const Face&f,int v){return f.v[0]==v||f.v[1]==v||f.v[2]==v;}
static inline int thirdOf(const Face&f,int a,int b){for(int k=0;k<3;k++){int x=f.v[k];if(x!=a&&x!=b)return x;}return -1;}
static inline Vec3 curCross(const Face&f){return cross3(P[f.v[1]]-P[f.v[0]],P[f.v[2]]-P[f.v[0]]);} 
static void initWork(){P=ORG.V;F=ORG.F;aliveV.assign(N0,1);aliveF.assign(M0,1);coverR.assign(N0,0);aliveFaces=M0;aliveVerts=N0;markA.assign(N0,0);markB.assign(N0,0);stampA=stampB=1;}
static void rebuildAdj(){vector<int> deg(P.size(),0);for(int i=0;i<(int)F.size();i++)if(aliveF[i]){Face f=F[i];if(f.v[0]>=0&&f.v[1]>=0&&f.v[2]>=0&&f.v[0]<(int)P.size()&&f.v[1]<(int)P.size()&&f.v[2]<(int)P.size()){deg[f.v[0]]++;deg[f.v[1]]++;deg[f.v[2]]++;}}
    adj.assign(P.size(),{});for(int i=0;i<(int)P.size();i++)adj[i].reserve(deg[i]);for(int i=0;i<(int)F.size();i++)if(aliveF[i]){Face f=F[i];if(f.v[0]>=0&&f.v[1]>=0&&f.v[2]>=0&&f.v[0]<(int)P.size()&&f.v[1]<(int)P.size()&&f.v[2]<(int)P.size()){adj[f.v[0]].push_back(i);adj[f.v[1]].push_back(i);adj[f.v[2]].push_back(i);}}}
static bool findEdgeFaces(int u,int v,int ef[2],int opp[2]){int cnt=0;if(u<0||v<0||u>=(int)P.size()||v>=(int)P.size())return false;for(int fid:adj[u]){if(!aliveF[fid])continue;const Face&f=F[fid];if(!faceHas(f,u)||!faceHas(f,v))continue;if(cnt>=2)return false;ef[cnt]=fid;opp[cnt]=thirdOf(f,u,v);cnt++;}return cnt==2&&opp[0]>=0&&opp[1]>=0&&opp[0]!=opp[1];}
static bool linkOK(int u,int v,const int opp[2]){if(++stampA>2000000000){fill(markA.begin(),markA.end(),0);stampA=1;}if(++stampB>2000000000){fill(markB.begin(),markB.end(),0);stampB=1;}for(int fid:adj[u])if(aliveF[fid]&&faceHas(F[fid],u)){const Face&f=F[fid];for(int k=0;k<3;k++){int x=f.v[k];if(x!=u&&x!=v)markA[x]=stampA;}}
    int common=0;for(int fid:adj[v])if(aliveF[fid]&&faceHas(F[fid],v)){const Face&f=F[fid];for(int k=0;k<3;k++){int x=f.v[k];if(x==u||x==v)continue;if(markA[x]!=stampA)continue;if(x!=opp[0]&&x!=opp[1])return false;if(markB[x]!=stampB){markB[x]=stampB;common++;}}}return common==2&&markB[opp[0]]==stampB&&markB[opp[1]]==stampB;}
static bool duplicateAfter(int fid,int rem,int ef0,int ef1,int a,int b,int c){int probe=a;if(adj[b].size()<adj[probe].size())probe=b;if(adj[c].size()<adj[probe].size())probe=c;for(int g:adj[probe]){if(!aliveF[g]||g==fid||g==ef0||g==ef1)continue;if(faceHas(F[g],rem))continue;if(sameTri(F[g],a,b,c))return true;}return false;}
struct Par{double epsR,planeR,minCos,areaR;};
struct Eval{bool ok=false;double cost=1e100,newR=0;int keep=-1,rem=-1;};
static Eval evalEndpoint(int keep,int rem,const int ef[2],const Par&par){Eval r;r.keep=keep;r.rem=rem;double eps=par.epsR*DIAG,plane=par.planeR*DIAG;double d=norm3(P[keep]-P[rem]);r.newR=max(coverR[keep],coverR[rem]+d);if(r.newR>eps+1e-12)return r;double minA2=max(1e-300,par.areaR*DIAG*DIAG*DIAG*DIAG);double maxPlane=0,maxAng=0,areaLoss=0;int changed=0;for(int fid:adj[rem]){if(!aliveF[fid]||!faceHas(F[fid],rem))continue;if(fid==ef[0]||fid==ef[1])continue;if(faceHas(F[fid],keep))return r;Face old=F[fid],nf=old;for(int k=0;k<3;k++)if(nf.v[k]==rem)nf.v[k]=keep;if(nf.v[0]==nf.v[1]||nf.v[0]==nf.v[2]||nf.v[1]==nf.v[2])return r;Vec3 co=curCross(old),cn=curCross(nf);double ao=norm3(co),an=norm3(cn);if(!(ao>0&&an>0)||norm2(cn)<=minA2)return r;double nd=clampd(dot3(co,cn)/(ao*an),-1,1);if(nd<par.minCos)return r;Vec3 no=co*(1.0/ao);double pd=fabs(dot3(no,P[keep]-P[old.v[0]]));if(pd>plane)return r;if(duplicateAfter(fid,rem,ef[0],ef[1],nf.v[0],nf.v[1],nf.v[2]))return r;maxPlane=max(maxPlane,pd);maxAng=max(maxAng,1.0-nd);areaLoss=max(areaLoss,max(0.0,1.0-an/ao));changed++;}if(!changed)return r;r.ok=true;r.cost=1.0*(r.newR/(eps+1e-300))+0.7*(maxPlane/(plane+1e-300))+180.0*maxAng+0.04*areaLoss+0.0002*changed;return r;}
static void applyCollapse(const Eval&e,const int ef[2]){int keep=e.keep,rem=e.rem;for(int i=0;i<2;i++)if(ef[i]>=0&&aliveF[ef[i]]){aliveF[ef[i]]=0;aliveFaces--;}for(int fid:adj[rem]){if(!aliveF[fid]||!faceHas(F[fid],rem))continue;for(int k=0;k<3;k++)if(F[fid].v[k]==rem)F[fid].v[k]=keep;}aliveV[rem]=0;aliveVerts--;coverR[keep]=e.newR;vector<int> merged;merged.reserve(adj[keep].size()+adj[rem].size());for(int fid:adj[keep])if(aliveF[fid]&&faceHas(F[fid],keep))merged.push_back(fid);for(int fid:adj[rem])if(aliveF[fid]&&faceHas(F[fid],keep))merged.push_back(fid);sort(merged.begin(),merged.end());merged.erase(unique(merged.begin(),merged.end()),merged.end());adj[keep].swap(merged);vector<int>().swap(adj[rem]);}
static bool tryCollapseEdge(int a,int b,const Par&par){if(a==b||a<0||b<0||a>=(int)P.size()||b>=(int)P.size()||!aliveV[a]||!aliveV[b])return false;int ef[2]={-1,-1},opp[2]={-1,-1};if(!findEdgeFaces(a,b,ef,opp)||!linkOK(a,b,opp))return false;Eval x=evalEndpoint(a,b,ef,par),y=evalEndpoint(b,a,ef,par);if(!x.ok&&!y.ok)return false;if(y.ok&&(!x.ok||y.cost<x.cost))applyCollapse(y,ef);else applyCollapse(x,ef);return true;}
static void collectEdges(vector<uint64_t>&keys){keys.clear();keys.reserve((size_t)aliveFaces*3+16);for(int i=0;i<(int)F.size();i++)if(aliveF[i]){const Face&f=F[i];if(aliveV[f.v[0]]&&aliveV[f.v[1]]&&aliveV[f.v[2]]){keys.push_back(ekey(f.v[0],f.v[1]));keys.push_back(ekey(f.v[1],f.v[2]));keys.push_back(ekey(f.v[2],f.v[0]));}}sort(keys.begin(),keys.end());keys.erase(unique(keys.begin(),keys.end()),keys.end());}
struct ERec{float l2;uint64_t k;};
static bool edgePass(const Par&par,double tim,int rounds,int target){bool any=false;vector<uint64_t>keys;vector<ERec>edges;for(int r=0;r<rounds&&elapsed()<tim&&aliveVerts>target;r++){collectEdges(keys);if(keys.empty())break;edges.clear();edges.reserve(keys.size());for(uint64_t k:keys){int a=(int)(k>>32),b=(int)(uint32_t)k;if(aliveV[a]&&aliveV[b])edges.push_back({(float)norm2(P[a]-P[b]),k});}sort(edges.begin(),edges.end(),[](const ERec&a,const ERec&b){return a.l2<b.l2;});int hit=0;for(const ERec&e:edges){if(elapsed()>=tim||aliveVerts<=target)break;int a=(int)(e.k>>32),b=(int)(uint32_t)e.k;if(a>=0&&b>=0&&a<(int)P.size()&&b<(int)P.size()&&aliveV[a]&&aliveV[b]&&tryCollapseEdge(a,b,par)){hit++;any=true;}}if(!hit)break;if(r+1<rounds)rebuildAdj();}return any;}
struct Smooth{double sm=0,co=0,sh=0;int tot=0;};
static Smooth smoothStats(){Smooth s;rebuildAdj();int stride=max(1,M0/60000);double c10=cos(10.0*M_PI/180.0),c30=cos(30.0*M_PI/180.0),c45=cos(45.0*M_PI/180.0);int sm=0,co=0,sh=0,tot=0;for(int fid=0;fid<M0;fid+=stride){const Face&f=F[fid];int e[3][2]={{f.v[0],f.v[1]},{f.v[1],f.v[2]},{f.v[2],f.v[0]}};for(int k=0;k<3;k++){int ef[2],op[2];if(!findEdgeFaces(e[k][0],e[k][1],ef,op))continue;Vec3 n0=unit3(curCross(F[ef[0]])),n1=unit3(curCross(F[ef[1]]));double d=clampd(dot3(n0,n1),-1,1);if(d>c10)sm++;if(d>c30)co++;if(d<c45)sh++;tot++;}}s.tot=tot;if(tot){s.sm=(double)sm/tot;s.co=(double)co/tot;s.sh=(double)sh/tot;}return s;}
static Mesh currentMesh(){vector<int> rem(P.size(),-1);Mesh m;m.V.reserve(aliveVerts);for(int i=0;i<(int)P.size();i++)if(aliveV[i]){rem[i]=(int)m.V.size();m.V.push_back(P[i]);}m.F.reserve(aliveFaces);for(int i=0;i<(int)F.size();i++)if(aliveF[i]){Face f=F[i];if(rem[f.v[0]]>=0&&rem[f.v[1]]>=0&&rem[f.v[2]]>=0)m.F.push_back(Face{{rem[f.v[0]],rem[f.v[1]],rem[f.v[2]]}});}return m;}
static Mesh simplifyGeneric(){
    initWork();rebuildAdj();Smooth st=smoothStats();int target=max(8,N0/200);vector<Par> ps;
    ps.push_back({0.030,0.006,0.9990,1e-24});ps.push_back({0.040,0.012,0.9950,1e-24});ps.push_back({0.046,0.020,0.9850,1e-24});
    if(st.co>0.55||N0>20000)ps.push_back({0.0490,0.030,0.9650,1e-24});
    if(st.sm>0.70&&st.sh<0.15)ps.push_back({0.0493,0.045,0.9250,1e-24});
    double limit=19.2;for(size_t i=0;i<ps.size()&&elapsed()<limit;i++){edgePass(ps[i],limit,2,target);rebuildAdj();}
    Mesh m=currentMesh();if(!validateMesh(m)||!coverageOK(m.V,0.0495))return ORG;orientLikeOriginal(m);return m;
}

static bool tryStructured(Mesh&out){
    Mesh best=ORG;bool ok=false;
    Mesh cand;
    if(tryTube(cand)&&cand.V.size()<best.V.size()){best=cand;ok=true;}
    if(elapsed()<16.0&&tryPolar(cand)&&cand.V.size()<best.V.size()){best=cand;ok=true;}
    if(elapsed()<16.0&&tryCluster(cand)&&cand.V.size()<best.V.size()){best=cand;ok=true;}
    if(ok){out=best;return true;}return false;
}

int main(){
    T0=chrono::steady_clock::now();
    if(!readInput())return 0;
    if(N0==8&&M0==12){outputOriginal();return 0;}
    Mesh ans;
    if(tryStructured(ans)&&validateMesh(ans)&&coverageOK(ans.V,0.0495)){outputMesh(ans);return 0;}
    ans=simplifyGeneric();
    if(!validateMesh(ans))ans=ORG;
    outputMesh(ans);
    return 0;
}
