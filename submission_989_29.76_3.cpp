#include <bits/stdc++.h>
using namespace std;

struct Vec3{ double x,y,z; };
static inline Vec3 operator+(const Vec3&a,const Vec3&b){return {a.x+b.x,a.y+b.y,a.z+b.z};}
static inline Vec3 operator-(const Vec3&a,const Vec3&b){return {a.x-b.x,a.y-b.y,a.z-b.z};}
static inline Vec3 operator*(const Vec3&a,double s){return {a.x*s,a.y*s,a.z*s};}
static inline Vec3 operator/(const Vec3&a,double s){return {a.x/s,a.y/s,a.z/s};}
static inline double dot3(const Vec3&a,const Vec3&b){return a.x*b.x+a.y*b.y+a.z*b.z;}
static inline Vec3 cross3(const Vec3&a,const Vec3&b){return {a.y*b.z-a.z*b.y,a.z*b.x-a.x*b.z,a.x*b.y-a.y*b.x};}
static inline double norm2(const Vec3&a){return dot3(a,a);}
static inline double norm3(const Vec3&a){return sqrt(max(0.0,norm2(a)));}
static inline Vec3 safe_unit(Vec3 v){double n=norm3(v); if(n<=1e-300) return {0,0,0}; return v*(1.0/n);}
static inline double clampd(double x,double a,double b){return x<a?a:(x>b?b:x);}

struct Face{ int v[3]; };
struct FastInput{
    vector<char> buf; char *p=nullptr;
    FastInput(){
        buf.reserve(1<<27); char tmp[1<<16]; size_t n;
        while((n=fread(tmp,1,sizeof(tmp),stdin))>0) buf.insert(buf.end(),tmp,tmp+n);
        buf.push_back('\0'); p=buf.data();
    }
    inline void skip(){ while(*p==' '||*p=='\n'||*p=='\r'||*p=='\t') ++p; }
    int nextInt(){ skip(); int s=1; if(*p=='-'){s=-1;++p;} int x=0; while(*p>='0'&&*p<='9'){x=x*10+(*p-'0');++p;} return x*s; }
    double nextDouble(){ skip(); char *e=nullptr; double x=strtod(p,&e); p=e; return x; }
    char nextChar(){ skip(); return *p++; }
};

static int N,M;
static vector<Vec3> V0;
static vector<Face> F0;
static Vec3 BMin,BMax;
static double DIAG=1.0;
static chrono::steady_clock::time_point T0;
static inline double elapsed(){ return chrono::duration<double>(chrono::steady_clock::now()-T0).count(); }

static inline array<int,3> key3(int a,int b,int c){ array<int,3> k{a,b,c}; sort(k.begin(),k.end()); return k; }
static inline bool sameTri(const Face&f,int a,int b,int c){ return key3(f.v[0],f.v[1],f.v[2])==key3(a,b,c); }

static void readInput(){
    FastInput in;
    N=in.nextInt(); M=in.nextInt();
    V0.resize(N); F0.resize(M);
    BMin={1e100,1e100,1e100}; BMax={-1e100,-1e100,-1e100};
    for(int i=0;i<N;i++){
        char c=in.nextChar(); (void)c;
        V0[i].x=in.nextDouble(); V0[i].y=in.nextDouble(); V0[i].z=in.nextDouble();
        BMin.x=min(BMin.x,V0[i].x); BMin.y=min(BMin.y,V0[i].y); BMin.z=min(BMin.z,V0[i].z);
        BMax.x=max(BMax.x,V0[i].x); BMax.y=max(BMax.y,V0[i].y); BMax.z=max(BMax.z,V0[i].z);
    }
    for(int i=0;i<M;i++){
        char c=in.nextChar(); (void)c;
        F0[i].v[0]=in.nextInt()-1; F0[i].v[1]=in.nextInt()-1; F0[i].v[2]=in.nextInt()-1;
    }
    DIAG=norm3(BMax-BMin); if(!(DIAG>0)) DIAG=1.0;
}

static vector<Vec3> vertexNormals(){
    vector<Vec3> n(N,{0,0,0});
    for(const Face&f:F0){
        Vec3 cr=cross3(V0[f.v[1]]-V0[f.v[0]], V0[f.v[2]]-V0[f.v[0]]);
        n[f.v[0]]=n[f.v[0]]+cr; n[f.v[1]]=n[f.v[1]]+cr; n[f.v[2]]=n[f.v[2]]+cr;
    }
    return n;
}

static void orientFace(vector<Face>&F,const vector<Vec3>&X,Face f,const Vec3&ref){
    Vec3 cr=cross3(X[f.v[1]]-X[f.v[0]], X[f.v[2]]-X[f.v[0]]);
    if(dot3(cr,ref)<0) swap(f.v[1],f.v[2]);
    F.push_back(f);
}

static void writeMesh(const vector<Vec3>&X,const vector<Face>&F){
    printf("%d %d\n",(int)X.size(),(int)F.size());
    for(const Vec3&p:X) printf("v %.10g %.10g %.10g\n",p.x,p.y,p.z);
    for(const Face&f:F) printf("f %d %d %d\n",f.v[0]+1,f.v[1]+1,f.v[2]+1);
}

static bool adj3cyc(const int a[3],int m,int&base){
    for(int t=0;t<3;t++) for(int s=0;s<2;s++){
        int b=(a[t]-s+m)%m; bool ok=true;
        for(int i=0;i<3;i++){ int d=(a[i]-b+m)%m; if(d!=0&&d!=1){ ok=false; break; } }
        if(ok){ base=b; return true; }
    }
    return false;
}
static bool torusFaceOK(const Face&f,int S){
    if(S<3 || N%S) return false; int U=N/S; if(U<3) return false;
    int r[3]={f.v[0]/S,f.v[1]/S,f.v[2]/S};
    int c[3]={f.v[0]%S,f.v[1]%S,f.v[2]%S};
    int rb=0,cb=0; if(!adj3cyc(r,U,rb)||!adj3cyc(c,S,cb)) return false;
    int mask=0;
    for(int i=0;i<3;i++){
        int x=(r[i]-rb+U)%U, y=(c[i]-cb+S)%S;
        if(x>1||y>1) return false;
        mask |= 1<<(x*2+y);
    }
    return __builtin_popcount((unsigned)mask)==3;
}
static bool detectTorus(int&S,int&U){
    if(N<36 || M!=2*N) return false;
    unordered_map<int,int> freq; freq.reserve(256);
    int st=max(1,M/3000);
    for(int i=0;i<M;i+=st){
        int a[3]={F0[i].v[0],F0[i].v[1],F0[i].v[2]};
        for(int k=0;k<3;k++){
            int d=abs(a[k]-a[(k+1)%3]); if(!d) continue; d=min(d,N-d);
            if(d>=3 && d<=N/3) freq[d]++;
        }
    }
    vector<pair<int,int>> q; q.reserve(freq.size());
    for(auto &p:freq) q.push_back({p.second,p.first});
    sort(q.rbegin(),q.rend());
    vector<int> cand;
    auto add=[&](int s){ if(s>=3 && s<=N/3 && N%s==0 && find(cand.begin(),cand.end(),s)==cand.end()) cand.push_back(s); };
    for(int i=0;i<(int)q.size() && i<20;i++){
        int d=q[i].second; for(int e=-3;e<=3;e++) add(d+e); if(d) add(N/d);
    }
    for(int s=3;s*s<=N;s++) if(N%s==0){ if((int)cand.size()<80){add(s); add(N/s);} }
    for(int s:cand){
        int u=N/s; if(u<3) continue;
        int sample=max(1,M/8000), tot=0, ok=0;
        for(int i=0;i<M;i+=sample){ tot++; if(torusFaceOK(F0[i],s)) ok++; }
        if(tot==0 || ok*1000 < tot*995) continue;
        bool all=true; int checkStep=max(1,M/250000);
        for(int i=0;i<M;i+=checkStep){ if(!torusFaceOK(F0[i],s)){ all=false; break; } }
        if(!all) continue;
        S=s; U=u; return true;
    }
    return false;
}

static bool latFaceOK(const Face&f,int S,int R){
    const int south=N-1;
    int a=f.v[0],b=f.v[1],c=f.v[2];
    int arr[3]={a,b,c};
    auto ring=[&](int id)->int{ if(id==0) return -1; if(id==south) return R; return (id-1)/S; };
    auto col=[&](int id)->int{ if(id==0||id==south) return -1; return (id-1)%S; };
    bool hasN=(a==0||b==0||c==0), hasS=(a==south||b==south||c==south);
    if(hasN&&hasS) return false;
    if(hasN){
        int cs[2], t=0; for(int id:arr) if(id!=0){ if(id<=0||id>=south||ring(id)!=0) return false; cs[t++]=col(id); }
        if(t!=2) return false; int d=(cs[0]-cs[1]+S)%S; return d==1||d==S-1;
    }
    if(hasS){
        int cs[2], t=0; for(int id:arr) if(id!=south){ if(id<=0||id>=south||ring(id)!=R-1) return false; cs[t++]=col(id); }
        if(t!=2) return false; int d=(cs[0]-cs[1]+S)%S; return d==1||d==S-1;
    }
    int rr[3]={ring(a),ring(b),ring(c)}, cc[3]={col(a),col(b),col(c)};
    for(int i=0;i<3;i++) if(rr[i]<0||rr[i]>=R) return false;
    int rb=*min_element(rr,rr+3); if(rb<0||rb+1>=R) return false;
    bool rOK=true; for(int i=0;i<3;i++) if(rr[i]!=rb&&rr[i]!=rb+1) rOK=false; if(!rOK) return false;
    int cb=0; if(!adj3cyc(cc,S,cb)) return false;
    int mask=0; for(int i=0;i<3;i++){ int x=rr[i]-rb,y=(cc[i]-cb+S)%S; if(x>1||y>1)return false; mask|=1<<(x*2+y); }
    return __builtin_popcount((unsigned)mask)==3;
}
static bool detectLatLong(int&R,int&S){
    if(N<10 || M!=2*(N-2)) return false;
    vector<int> cand;
    for(int s=3; (long long)s*s<=N-2; ++s) if((N-2)%s==0){ cand.push_back(s); if(s!=(N-2)/s) cand.push_back((N-2)/s); }
    sort(cand.begin(),cand.end(),[](int a,int b){return a>b;});
    for(int s:cand){
        int r=(N-2)/s; if(r<2||s<3) continue;
        int sample=max(1,M/8000), tot=0, ok=0;
        for(int i=0;i<M;i+=sample){ tot++; if(latFaceOK(F0[i],s,r)) ok++; }
        if(tot==0 || ok*1000 < tot*995) continue;
        bool all=true; int checkStep=max(1,M/250000);
        for(int i=0;i<M;i+=checkStep){ if(!latFaceOK(F0[i],s,r)){ all=false; break; } }
        if(!all) continue;
        R=r; S=s; return true;
    }
    return false;
}

static vector<int> chooseCyclic(int n,int k){
    k=max(3,min(k,n)); vector<int> r; r.reserve(k);
    for(int i=0;i<k;i++){ int v=(int)llround((double)i*n/k); if(v>=n) v=n-1; if(!r.empty()&&v==r.back()) v=min(n-1,v+1); r.push_back(v); }
    sort(r.begin(),r.end()); r.erase(unique(r.begin(),r.end()),r.end());
    while((int)r.size()<k){ for(int i=0;i<n && (int)r.size()<k;i++) if(!binary_search(r.begin(),r.end(),i)){ r.insert(lower_bound(r.begin(),r.end(),i),i); } }
    return r;
}
static vector<int> chooseLinear(int n,int k){
    k=max(1,min(k,n)); vector<int> r; r.reserve(k);
    if(k==1){ r.push_back(n/2); return r; }
    for(int i=0;i<k;i++){ int v=(int)llround((double)i*(n-1)/(k-1)); r.push_back(min(n-1,max(0,v))); }
    sort(r.begin(),r.end()); r.erase(unique(r.begin(),r.end()),r.end());
    while((int)r.size()<k){ for(int i=0;i<n && (int)r.size()<k;i++) if(!binary_search(r.begin(),r.end(),i)){ r.insert(lower_bound(r.begin(),r.end(),i),i); } }
    return r;
}
static pair<int,int> nearest2linear(const vector<int>&sel,int x){
    int p=lower_bound(sel.begin(),sel.end(),x)-sel.begin();
    int a=sel[min(p,(int)sel.size()-1)], b=sel[max(0,p-1)];
    return {a,b};
}
static pair<int,int> nearest2cyclic(const vector<int>&sel,int x,int n){
    int p=lower_bound(sel.begin(),sel.end(),x)-sel.begin();
    int cand[4]={sel[p%(int)sel.size()], sel[(p-1+(int)sel.size())%(int)sel.size()], sel[(p+1)%(int)sel.size()], sel[(p-2+2*(int)sel.size())%(int)sel.size()]};
    int best=cand[0], second=cand[1]; int bd=INT_MAX, sd=INT_MAX;
    for(int t=0;t<4;t++){ int d=abs(cand[t]-x); d=min(d,n-d); if(d<bd){sd=bd; second=best; bd=d; best=cand[t];} else if(cand[t]!=best && d<sd){sd=d; second=cand[t];} }
    return {best,second};
}

static double torusCoverage(int S,int U,const vector<int>&rows,const vector<int>&cols){
    double mx=0;
    for(int i=0;i<U;i++){
        auto rr=nearest2cyclic(rows,i,U);
        for(int j=0;j<S;j++){
            auto cc=nearest2cyclic(cols,j,S); int id=i*S+j;
            int ids[4]={rr.first*S+cc.first, rr.first*S+cc.second, rr.second*S+cc.first, rr.second*S+cc.second};
            double bd=1e100; for(int q=0;q<4;q++) bd=min(bd,norm2(V0[id]-V0[ids[q]]));
            mx=max(mx,bd);
        }
    }
    return sqrt(mx);
}
static double latCoverage(int R,int S,const vector<int>&rows,const vector<int>&cols){
    double mx=0; const int south=N-1;
    for(int r=0;r<R;r++){
        auto rr=nearest2linear(rows,r);
        for(int j=0;j<S;j++){
            auto cc=nearest2cyclic(cols,j,S); int id=1+r*S+j;
            int ids[6]={1+rr.first*S+cc.first,1+rr.first*S+cc.second,1+rr.second*S+cc.first,1+rr.second*S+cc.second,0,south};
            double bd=1e100; for(int q=0;q<6;q++) bd=min(bd,norm2(V0[id]-V0[ids[q]]));
            mx=max(mx,bd);
        }
    }
    return sqrt(mx);
}

static bool generatedVerticesNearOriginal(const vector<Vec3>&X,double lim){
    if(X.empty()) return false;
    const double cell=max(lim,1e-12); const double lim2=lim*lim;
    struct H{ size_t operator()(const long long&a)const{return (size_t)(a^(a>>32));} };
    unordered_map<long long, vector<int>, H> buckets; buckets.reserve((size_t)N*2+1);
    auto pack=[&](int ix,int iy,int iz){ const long long O=1LL<<20; return ((ix+O)<<42)^((iy+O)<<21)^(iz+O); };
    auto get=[&](const Vec3&p,int&ix,int&iy,int&iz){ ix=(int)floor((p.x-BMin.x)/cell); iy=(int)floor((p.y-BMin.y)/cell); iz=(int)floor((p.z-BMin.z)/cell); };
    for(int i=0;i<N;i++){int ix,iy,iz; get(V0[i],ix,iy,iz); buckets[pack(ix,iy,iz)].push_back(i);}
    for(const Vec3&p:X){
        int ix,iy,iz; get(p,ix,iy,iz); double best=1e100;
        for(int dx=-1;dx<=1;dx++)for(int dy=-1;dy<=1;dy++)for(int dz=-1;dz<=1;dz++){
            auto it=buckets.find(pack(ix+dx,iy+dy,iz+dz)); if(it==buckets.end()) continue;
            for(int id:it->second) best=min(best,norm2(p-V0[id]));
        }
        if(best>lim2) return false;
    }
    return true;
}

static void makeEllipsoidMesh(const Vec3&c,const Vec3&ax,int lat,int lon,vector<Vec3>&X,vector<Face>&OF){
    X.clear(); OF.clear(); const double pi=acos(-1.0);
    X.reserve(2+(lat-1)*lon); OF.reserve(2*lat*lon);
    auto pt=[&](double th,double ph){return Vec3{c.x+ax.x*sin(th)*cos(ph),c.y+ax.y*sin(th)*sin(ph),c.z+ax.z*cos(th)};};
    X.push_back({c.x,c.y,c.z+ax.z});
    for(int r=1;r<lat;r++){ double th=pi*(double)r/(double)lat; for(int j=0;j<lon;j++){ double ph=2*pi*(double)j/(double)lon; X.push_back(pt(th,ph)); }}
    int bot=(int)X.size(); X.push_back({c.x,c.y,c.z-ax.z});
    auto id=[&](int r,int j){return 1+(r-1)*lon+((j%lon+lon)%lon);};
    auto add=[&](int a,int b,int d){ Face f{{a,b,d}}; Vec3 cen=(X[a]+X[b]+X[d])/3.0; orientFace(OF,X,f,cen-c); };
    for(int j=0;j<lon;j++) add(0,id(1,j),id(1,j+1));
    for(int r=1;r<lat-1;r++) for(int j=0;j<lon;j++){ int a=id(r,j), b=id(r+1,j), cc=id(r+1,j+1), d=id(r,j+1); add(a,b,d); add(b,cc,d); }
    for(int j=0;j<lon;j++) add(bot,id(lat-1,j+1),id(lat-1,j));
}

static bool ellipsoidCoversOriginal(const Vec3&c,const Vec3&ax,int lat,int lon,double lim){
    const double pi=acos(-1.0), lim2=lim*lim;
    auto ep=[&](int r,int j)->Vec3{
        if(r<=0) return {c.x,c.y,c.z+ax.z};
        if(r>=lat) return {c.x,c.y,c.z-ax.z};
        double th=pi*(double)r/(double)lat, ph=2*pi*(double)((j%lon+lon)%lon)/(double)lon;
        return {c.x+ax.x*sin(th)*cos(ph),c.y+ax.y*sin(th)*sin(ph),c.z+ax.z*cos(th)};
    };
    for(const Vec3&p:V0){
        double x=(p.x-c.x)/ax.x, y=(p.y-c.y)/ax.y, z=(p.z-c.z)/ax.z;
        double rr=sqrt(max(0.0,x*x+y*y+z*z)); if(rr<1e-12) return false;
        x/=rr; y/=rr; z/=rr;
        double th=acos(clampd(z,-1.0,1.0)); double ph=atan2(y,x); if(ph<0) ph+=2*pi;
        int ri=(int)llround(th/pi*lat); int ji=(int)llround(ph/(2*pi)*lon);
        double bd=1e100;
        for(int dr=-2;dr<=2;dr++) for(int dj=-2;dj<=2;dj++){ int r=ri+dr; if(r<0)r=0; if(r>lat)r=lat; Vec3 q=ep(r,ji+dj); bd=min(bd,norm2(p-q)); }
        if(bd>lim2) return false;
    }
    return true;
}

static bool tryEllipsoidRoute(vector<Vec3>&X, vector<Face>&OF){
    if(N<900) return false;
    Vec3 c=(BMin+BMax)*0.5; Vec3 ax=(BMax-BMin)*0.5;
    if(ax.x<=1e-12||ax.y<=1e-12||ax.z<=1e-12) return false;
    double mx=max(ax.x,max(ax.y,ax.z)), mn=min(ax.x,min(ax.y,ax.z));
    if(mn<mx*0.12) return false;
    int stride=max(1,N/240000); int cnt=0; double ss=0, ma=0;
    for(int i=0;i<N;i+=stride){
        double x=(V0[i].x-c.x)/ax.x, y=(V0[i].y-c.y)/ax.y, z=(V0[i].z-c.z)/ax.z;
        double e=fabs(sqrt(max(0.0,x*x+y*y+z*z))-1.0); ss+=e*e; ma=max(ma,e); cnt++;
    }
    if(cnt<200) return false;
    double rms=sqrt(ss/cnt);
    double maLim=(N<5000?0.010:0.018), rmsLim=(N<5000?0.0038:0.0065);
    if(ma>maLim || rms>rmsLim) return false;
    struct T{int lat,lon;};
    vector<T> trials;
    if(N<3000) trials={{12,24},{14,28},{16,32}};
    else if(N<20000) trials={{14,28},{16,32},{18,36},{20,40}};
    else if(N<100000) trials={{18,36},{20,40},{24,48},{28,56}};
    else trials={{24,48},{28,56},{32,64},{36,72}};
    double lim=0.0475*DIAG;
    for(const T&t:trials){
        int nv=2+(t.lat-1)*t.lon; if(nv>=N) continue;
        vector<Vec3>C; vector<Face>F; makeEllipsoidMesh(c,ax,t.lat,t.lon,C,F);
        if(ellipsoidCoversOriginal(c,ax,t.lat,t.lon,lim) && generatedVerticesNearOriginal(C,lim)){ X.swap(C); OF.swap(F); return true; }
    }
    return false;
}

static bool tryCapsuleRoute(vector<Vec3>&X, vector<Face>&OF){
    if(N<1500) return false;
    Vec3 ext=BMax-BMin; int ax=0; if(ext.y>ext.x&&ext.y>=ext.z) ax=1; else if(ext.z>ext.x&&ext.z>=ext.y) ax=2;
    int u=(ax==0?1:0), v=(ax==2?1:2); if(ax==1){u=0;v=2;}
    double c0=(ax==0?(BMin.x+BMax.x):(ax==1?(BMin.y+BMax.y):(BMin.z+BMax.z)))*0.5;
    double cu=(u==0?(BMin.x+BMax.x):(u==1?(BMin.y+BMax.y):(BMin.z+BMax.z)))*0.5;
    double cv=(v==0?(BMin.x+BMax.x):(v==1?(BMin.y+BMax.y):(BMin.z+BMax.z)))*0.5;
    double eu=(u==0?ext.x:(u==1?ext.y:ext.z)), ev=(v==0?ext.x:(v==1?ext.y:ext.z)), et=(ax==0?ext.x:(ax==1?ext.y:ext.z));
    if(min(eu,ev)<max(eu,ev)*0.86 || et<max(eu,ev)*1.35) return false;
    double R=0.25*(eu+ev); double H=0.5*et-R; if(!(R>1e-9)||H<=R*0.05) return false;
    auto comp=[&](const Vec3&p,int a){return a==0?p.x:(a==1?p.y:p.z);};
    int stride=max(1,N/240000), cnt=0; double ss=0, ma=0; int side=0,cap=0;
    for(int i=0;i<N;i+=stride){
        double t=comp(V0[i],ax)-c0, x=comp(V0[i],u)-cu, y=comp(V0[i],v)-cv; double at=fabs(t);
        double dist;
        if(at<=H){ dist=fabs(sqrt(x*x+y*y)-R); side++; }
        else{ double dt=at-H; dist=fabs(sqrt(x*x+y*y+dt*dt)-R); cap++; }
        double rel=dist/R; ss+=rel*rel; ma=max(ma,rel); cnt++;
    }
    if(cnt<200||side<cnt/5||cap<cnt/12||ma>0.045||sqrt(ss/cnt)>0.012) return false;
    int sides=(N<5000?32:(N<100000?48:64)), rings=(N<5000?12:(N<100000?16:20));
    int zrings=2*rings+3; int nv=2+zrings*sides; if(nv>=N) return false;
    X.clear(); OF.clear(); X.reserve(nv); OF.reserve(2*zrings*sides);
    auto make=[&](double t,double rr,int j){ double ph=2*acos(-1.0)*(double)j/sides; double a[3]={0,0,0}; a[ax]=c0+t; a[u]=cu+rr*cos(ph); a[v]=cv+rr*sin(ph); return Vec3{a[0],a[1],a[2]}; };
    double top[3]={0,0,0}, bot[3]={0,0,0}; top[ax]=c0+H+R; bot[ax]=c0-H-R; top[u]=bot[u]=cu; top[v]=bot[v]=cv;
    X.push_back({top[0],top[1],top[2]});
    vector<pair<double,double>> rr;
    for(int i=1;i<=rings;i++){ double phi=0.5*acos(-1.0)*(double)i/(rings+1); rr.push_back({H+R*cos(phi),R*sin(phi)}); }
    rr.push_back({H,R}); rr.push_back({0,R}); rr.push_back({-H,R});
    for(int i=rings;i>=1;i--){ double phi=0.5*acos(-1.0)*(double)i/(rings+1); rr.push_back({-H-R*cos(phi),R*sin(phi)}); }
    for(auto pr:rr) for(int j=0;j<sides;j++) X.push_back(make(pr.first,pr.second,j));
    int bottom=X.size(); X.push_back({bot[0],bot[1],bot[2]});
    auto id=[&](int r,int j){ return 1+r*sides+((j%sides+sides)%sides); };
    auto add=[&](int a,int b,int cidx){ Face f{{a,b,cidx}}; Vec3 cen=(X[a]+X[b]+X[cidx])/3.0; Vec3 ref=cen; double t=comp(cen,ax)-c0; double q[3]={0,0,0}; q[ax]=c0+clampd(t,-H,H); q[u]=cu; q[v]=cv; ref=cen-Vec3{q[0],q[1],q[2]}; orientFace(OF,X,f,ref); };
    int Rn=rr.size();
    for(int j=0;j<sides;j++) add(0,id(0,j+1),id(0,j));
    for(int r=0;r+1<Rn;r++)for(int j=0;j<sides;j++){int a=id(r,j),b=id(r+1,j),c=id(r+1,j+1),d=id(r,j+1);add(a,b,d);add(b,c,d);}
    for(int j=0;j<sides;j++) add(bottom,id(Rn-1,j),id(Rn-1,j+1));
    double lim=0.047*DIAG;
    if(!generatedVerticesNearOriginal(X,lim)) return false;
    return true;
}

static bool tryBoxGridRoute(vector<Vec3>&X, vector<Face>&OF){
    if(N<1200) return false;
    Vec3 L=BMax-BMin; double tol=max(1e-10,1e-7*DIAG);
    int oncnt=0;
    for(const Vec3&p:V0){
        double d=min({fabs(p.x-BMin.x),fabs(p.x-BMax.x),fabs(p.y-BMin.y),fabs(p.y-BMax.y),fabs(p.z-BMin.z),fabs(p.z-BMax.z)});
        if(d>tol) return false; oncnt++;
    }
    double step=max(1e-9,0.064*DIAG);
    int nx=max(1,(int)ceil(L.x/step)), ny=max(1,(int)ceil(L.y/step)), nz=max(1,(int)ceil(L.z/step));
    long long est=2LL*((long long)(nx+1)*(ny+1)+(long long)(nx+1)*(nz+1)+(long long)(ny+1)*(nz+1));
    if(est>=N || est>70000) return false;
    auto coord=[&](int axis,int i)->double{
        if(axis==0) return BMin.x + (BMax.x-BMin.x)*(double)i/max(1,nx);
        if(axis==1) return BMin.y + (BMax.y-BMin.y)*(double)i/max(1,ny);
        return BMin.z + (BMax.z-BMin.z)*(double)i/max(1,nz);
    };
    vector<Vec3>C; C.reserve((size_t)est);
    unordered_map<unsigned long long,int> mp; mp.reserve((size_t)est*2+10);
    auto key=[&](int i,int j,int k)->unsigned long long{ return ((unsigned long long)(unsigned)i<<42)^((unsigned long long)(unsigned)j<<21)^(unsigned)k; };
    auto addv=[&](int i,int j,int k)->int{
        unsigned long long kk=key(i,j,k); auto it=mp.find(kk); if(it!=mp.end()) return it->second;
        int id=(int)C.size(); mp[kk]=id; C.push_back({coord(0,i),coord(1,j),coord(2,k)}); return id;
    };
    vector<Face> CF;
    auto addFace=[&](int a,int b,int c,Vec3 ref){ Face f{{a,b,c}}; Vec3 cr=cross3(C[b]-C[a],C[c]-C[a]); if(dot3(cr,ref)<0) swap(f.v[1],f.v[2]); CF.push_back(f); };
    for(int side=0;side<6;side++){
        int ax=side/2, hi=side&1; Vec3 ref{0,0,0}; if(ax==0) ref.x=hi?1:-1; if(ax==1) ref.y=hi?1:-1; if(ax==2) ref.z=hi?1:-1;
        int A=(ax==0?1:0), B=(ax==2?1:2); if(ax==1){A=0;B=2;}
        int na=(A==0?nx:(A==1?ny:nz)), nb=(B==0?nx:(B==1?ny:nz));
        for(int i=0;i<na;i++) for(int j=0;j<nb;j++){
            int q[4][3];
            for(int t=0;t<4;t++){ q[t][0]=0; q[t][1]=0; q[t][2]=0; q[t][ax]=hi?(ax==0?nx:(ax==1?ny:nz)):0; }
            q[0][A]=i; q[0][B]=j; q[1][A]=i+1; q[1][B]=j; q[2][A]=i+1; q[2][B]=j+1; q[3][A]=i; q[3][B]=j+1;
            int a=addv(q[0][0],q[0][1],q[0][2]), b=addv(q[1][0],q[1][1],q[1][2]), c=addv(q[2][0],q[2][1],q[2][2]), d=addv(q[3][0],q[3][1],q[3][2]);
            addFace(a,b,d,ref); addFace(b,c,d,ref);
        }
    }
    double lim=0.047*DIAG, lim2=lim*lim;
    for(const Vec3&p:V0){
        double best=1e100;
        for(int ax=0;ax<3;ax++) for(int hi=0;hi<2;hi++){
            double plane=(ax==0?(hi?BMax.x:BMin.x):(ax==1?(hi?BMax.y:BMin.y):(hi?BMax.z:BMin.z)));
            double pd=fabs((ax==0?p.x:(ax==1?p.y:p.z))-plane); if(pd>tol*10) continue;
            double qx=p.x,qy=p.y,qz=p.z;
            if(ax!=0){ int ix=(int)llround((p.x-BMin.x)/max(1e-300,L.x)*nx); ix=max(0,min(nx,ix)); qx=coord(0,ix); } else qx=plane;
            if(ax!=1){ int iy=(int)llround((p.y-BMin.y)/max(1e-300,L.y)*ny); iy=max(0,min(ny,iy)); qy=coord(1,iy); } else qy=plane;
            if(ax!=2){ int iz=(int)llround((p.z-BMin.z)/max(1e-300,L.z)*nz); iz=max(0,min(nz,iz)); qz=coord(2,iz); } else qz=plane;
            best=min(best,(p.x-qx)*(p.x-qx)+(p.y-qy)*(p.y-qy)+(p.z-qz)*(p.z-qz));
        }
        if(best>lim2) return false;
    }
    double cell=lim;
    struct H{ size_t operator()(const long long&a)const{return (size_t)(a^(a>>32));} };
    unordered_map<long long, vector<int>, H> buckets; buckets.reserve(N*2u);
    auto packCell=[&](int ix,int iy,int iz){ const long long O=1LL<<20; return ((ix+O)<<42) ^ ((iy+O)<<21) ^ (iz+O); };
    auto cid=[&](const Vec3&p){ int ix=(int)floor((p.x-BMin.x)/cell), iy=(int)floor((p.y-BMin.y)/cell), iz=(int)floor((p.z-BMin.z)/cell); return packCell(ix,iy,iz); };
    auto cijk=[&](const Vec3&p,int&ix,int&iy,int&iz){ ix=(int)floor((p.x-BMin.x)/cell); iy=(int)floor((p.y-BMin.y)/cell); iz=(int)floor((p.z-BMin.z)/cell); };
    for(int i=0;i<N;i++) buckets[cid(V0[i])].push_back(i);
    for(const Vec3&p:C){
        int ix,iy,iz; cijk(p,ix,iy,iz); double bd=1e100;
        for(int dx=-1;dx<=1;dx++) for(int dy=-1;dy<=1;dy++) for(int dz=-1;dz<=1;dz++){
            long long kk=packCell(ix+dx,iy+dy,iz+dz);
            auto it=buckets.find(kk); if(it==buckets.end()) continue;
            for(int id:it->second) bd=min(bd,norm2(p-V0[id]));
        }
        if(bd>lim2) return false;
    }
    X.swap(C); OF.swap(CF); return X.size()>=4 && OF.size()>=4 && X.size()<(size_t)N;
}

static bool tryTorusRoute(vector<Vec3>&X,vector<Face>&OF){
    int S=0,U=0; if(!detectTorus(S,U)) return false;
    if(N<2000) return false;
    vector<Vec3> vn=vertexNormals();
    vector<int> targets;
    if(N<15000) targets={max(64,N/8),max(96,N/6),max(128,N/5),max(192,N/4)};
    else if(N<70000) targets={3072,4096,6144,8192,10240,12288,14336,16384};
    else targets={8192,12288,16384,24576,32768,49152};
    double covLimit = (N<70000?0.030:0.026)*DIAG;
    vector<int> bestR,bestC; int bestNV=N+1;
    for(int tg:targets){
        if(tg>=N) continue;
        double ar=sqrt((double)U/max(1,S));
        int U2=(int)llround(sqrt((double)tg)*ar); U2=max(6,min(U,U2));
        int S2=max(6,min(S,tg/max(1,U2)));
        for(int du=-2;du<=2;du++) for(int ds=-2;ds<=2;ds++){
            int u2=max(6,min(U,U2+du)), s2=max(6,min(S,S2+ds)); int nv=u2*s2; if(nv>=N||nv>=bestNV) continue;
            vector<int> rr=chooseCyclic(U,u2), cc=chooseCyclic(S,s2);
            double cov=torusCoverage(S,U,rr,cc);
            if(cov<=covLimit){ bestNV=nv; bestR.swap(rr); bestC.swap(cc); }
        }
    }
    if(bestNV>=N) return false;
    int U2=bestR.size(), S2=bestC.size();
    X.clear(); OF.clear(); X.reserve(U2*S2); vector<int> src(U2*S2);
    for(int i=0;i<U2;i++) for(int j=0;j<S2;j++){ int oid=bestR[i]*S+bestC[j]; src[i*S2+j]=oid; X.push_back(V0[oid]); }
    auto id=[&](int i,int j){ i=(i%U2+U2)%U2; j=(j%S2+S2)%S2; return i*S2+j; };
    OF.reserve(2*U2*S2);
    for(int i=0;i<U2;i++) for(int j=0;j<S2;j++){
        int a=id(i,j), b=id(i+1,j), c=id(i+1,j+1), d=id(i,j+1);
        Vec3 ref=vn[src[a]]+vn[src[b]]+vn[src[c]]+vn[src[d]];
        orientFace(OF,X,Face{{a,b,d}},ref); orientFace(OF,X,Face{{b,c,d}},ref);
    }
    return X.size()>=4 && OF.size()>=4;
}

static bool tryLatRoute(vector<Vec3>&X,vector<Face>&OF){
    int R=0,S=0; if(!detectLatLong(R,S)) return false;
    if(N<1000) return false;
    vector<Vec3> vn=vertexNormals();
    vector<int> targets;
    if(N<12000) targets={256,384,512,768,1024,1536};
    else targets={1024,1536,2048,3072,4096,6144,8192,10240};
    double covLimit=(N<40000?0.028:0.025)*DIAG;
    vector<int> bestR,bestC; int bestNV=N+1;
    for(int tg:targets){
        if(tg>=N) continue;
        double ar=sqrt((double)R/max(1,S));
        int R2=max(1,min(R,(int)llround(sqrt((double)tg)*ar)));
        int S2=max(3,min(S,tg/max(1,R2)));
        for(int dr=-2;dr<=2;dr++) for(int ds=-3;ds<=3;ds++){
            int r2=max(1,min(R,R2+dr)), s2=max(3,min(S,S2+ds)); int nv=2+r2*s2; if(nv>=N||nv>=bestNV) continue;
            vector<int> rr=chooseLinear(R,r2), cc=chooseCyclic(S,s2);
            double cov=latCoverage(R,S,rr,cc);
            if(cov<=covLimit){ bestNV=nv; bestR.swap(rr); bestC.swap(cc); }
        }
    }
    if(bestNV>=N) return false;
    int R2=bestR.size(), S2=bestC.size();
    X.clear(); OF.clear(); X.reserve(2+R2*S2); vector<int> src; src.reserve(2+R2*S2);
    X.push_back(V0[0]); src.push_back(0);
    for(int r=0;r<R2;r++) for(int j=0;j<S2;j++){ int oid=1+bestR[r]*S+bestC[j]; src.push_back(oid); X.push_back(V0[oid]); }
    int south=(int)X.size(); X.push_back(V0[N-1]); src.push_back(N-1);
    auto rid=[&](int r,int j){ return 1+r*S2+((j%S2+S2)%S2); };
    OF.reserve(2*S2*max(1,R2));
    for(int j=0;j<S2;j++){
        int a=rid(0,j), b=rid(0,j+1); Vec3 ref=vn[0]+vn[src[a]]+vn[src[b]];
        orientFace(OF,X,Face{{0,a,b}},ref);
    }
    for(int r=0;r+1<R2;r++) for(int j=0;j<S2;j++){
        int a=rid(r,j), b=rid(r+1,j), c=rid(r+1,j+1), d=rid(r,j+1);
        Vec3 ref=vn[src[a]]+vn[src[b]]+vn[src[c]]+vn[src[d]];
        orientFace(OF,X,Face{{a,b,d}},ref); orientFace(OF,X,Face{{b,c,d}},ref);
    }
    for(int j=0;j<S2;j++){
        int a=rid(R2-1,j), b=rid(R2-1,j+1); Vec3 ref=vn[N-1]+vn[src[a]]+vn[src[b]];
        orientFace(OF,X,Face{{a,south,b}},ref);
    }
    return X.size()>=4 && OF.size()>=4;
}

struct Simplifier{
    vector<Vec3> P;
    vector<Face> F;
    vector<unsigned char> av,af;
    vector<vector<int>> inc;
    vector<double> err;
    vector<int> mark;
    int aliveV=0, aliveF=0, stamp=1;
    double maxErr, maxPlane, minDot, minAreaRatio;
    struct Node{ double cost; int u,v; bool operator<(const Node&o)const{return cost>o.cost;} };
    priority_queue<Node> pq;

    Simplifier(){
        P=V0; F=F0; av.assign(N,1); af.assign(M,1); err.assign(N,0); inc.assign(N,{}); mark.assign(N,0); aliveV=N; aliveF=M;
        vector<int> deg(N,0); for(auto &f:F){deg[f.v[0]]++;deg[f.v[1]]++;deg[f.v[2]]++;}
        for(int i=0;i<N;i++) inc[i].reserve(deg[i]+8);
        for(int i=0;i<M;i++){inc[F[i].v[0]].push_back(i);inc[F[i].v[1]].push_back(i);inc[F[i].v[2]].push_back(i);}
        maxErr=0.049*DIAG;
        if(N<1000){ maxPlane=0.0015*DIAG; minDot=0.9990; }
        else if(N<20000){ maxPlane=0.0055*DIAG; minDot=0.985; }
        else if(N<80000){ maxPlane=0.0075*DIAG; minDot=0.972; }
        else { maxPlane=0.0090*DIAG; minDot=0.955; }
        minAreaRatio=1e-8;
    }
    inline bool hasv(int fid,int v)const{ const Face&f=F[fid]; return f.v[0]==v||f.v[1]==v||f.v[2]==v; }
    inline int third(int fid,int a,int b)const{ const Face&f=F[fid]; for(int i=0;i<3;i++){int x=f.v[i]; if(x!=a&&x!=b) return x;} return -1; }
    bool edgeFaces(int u,int v,int ef[2],int opp[2]){
        if(!av[u]||!av[v]) return false; int cnt=0;
        const vector<int>&L=(inc[u].size()<inc[v].size()?inc[u]:inc[v]);
        for(int fid:L){ if(!af[fid]) continue; if(hasv(fid,u)&&hasv(fid,v)){ if(cnt>=2) return false; ef[cnt]=fid; opp[cnt]=third(fid,u,v); cnt++; } }
        return cnt==2 && opp[0]>=0 && opp[1]>=0 && opp[0]!=opp[1];
    }
    bool linkOK(int u,int v,const int opp[2]){
        if(++stamp>2000000000){ fill(mark.begin(),mark.end(),0); stamp=1; }
        int st=stamp;
        for(int fid:inc[u]) if(af[fid]&&hasv(fid,u)){
            const Face&f=F[fid]; for(int k=0;k<3;k++){ int x=f.v[k]; if(x!=u&&x!=v) mark[x]=st; }
        }
        int common=0; bool seen0=false,seen1=false;
        if(++stamp>2000000000){ fill(mark.begin(),mark.end(),0); stamp=1; }
        int st2=stamp;
        for(int fid:inc[v]) if(af[fid]&&hasv(fid,v)){
            const Face&f=F[fid]; for(int k=0;k<3;k++){
                int x=f.v[k]; if(x==u||x==v) continue;
                if(mark[x]==st){ if(x!=opp[0]&&x!=opp[1]) return false; if(mark[x]!=st2){ mark[x]=st2; common++; if(x==opp[0])seen0=true; if(x==opp[1])seen1=true; } }
            }
        }
        return common==2 && seen0 && seen1;
    }
    Vec3 fNormalCurrent(int fid,int moved,const Vec3&np)const{
        const Face&f=F[fid]; Vec3 a=(f.v[0]==moved?np:P[f.v[0]]), b=(f.v[1]==moved?np:P[f.v[1]]), c=(f.v[2]==moved?np:P[f.v[2]]);
        return cross3(b-a,c-a);
    }
    bool duplicateAfter(int fid,int a,int b,int c,int skip0,int skip1,int rem)const{
        int verts[3]={a,b,c}; int base=verts[0];
        for(int t=1;t<3;t++) if(inc[verts[t]].size()<inc[base].size()) base=verts[t];
        auto key=key3(a,b,c);
        for(int g:inc[base]){
            if(!af[g]||g==fid||g==skip0||g==skip1) continue;
            int x[3]={F[g].v[0],F[g].v[1],F[g].v[2]};
            if(x[0]==rem||x[1]==rem||x[2]==rem) continue;
            if(key3(x[0],x[1],x[2])==key) return true;
        }
        return false;
    }
    struct Cand{ bool ok=false; double cost=1e100,newErr=0; Vec3 np{}; int keep=-1,rem=-1; };
    Cand eval(int keep,int rem,const int ef[2],bool midpoint){
        Cand c; c.keep=keep; c.rem=rem; c.np = midpoint? (P[keep]+P[rem])*0.5 : P[keep];
        c.newErr=max(err[keep]+norm3(P[keep]-c.np), err[rem]+norm3(P[rem]-c.np));
        if(c.newErr>maxErr) return c;
        vector<int> affected; affected.reserve(inc[keep].size()+inc[rem].size());
        for(int fid:inc[keep]) if(af[fid]&&hasv(fid,keep)&&fid!=ef[0]&&fid!=ef[1]) affected.push_back(fid);
        for(int fid:inc[rem]) if(af[fid]&&hasv(fid,rem)&&fid!=ef[0]&&fid!=ef[1]) affected.push_back(fid);
        sort(affected.begin(),affected.end()); affected.erase(unique(affected.begin(),affected.end()),affected.end());
        double worstN=0,worstP=0,worstA=0; int cnt=0;
        for(int fid:affected){
            Face old=F[fid]; int nv[3]={old.v[0],old.v[1],old.v[2]};
            for(int k=0;k<3;k++) if(nv[k]==rem) nv[k]=keep;
            int a=nv[0],b=nv[1],d=nv[2];
            if(a==b||a==d||b==d) return c;
            Vec3 oldCr=cross3(P[old.v[1]]-P[old.v[0]],P[old.v[2]]-P[old.v[0]]);
            Vec3 pa=(a==keep?c.np:P[a]), pb=(b==keep?c.np:P[b]), pc=(d==keep?c.np:P[d]);
            Vec3 newCr=cross3(pb-pa,pc-pa);
            double oa=norm3(oldCr), na=norm3(newCr); if(!(oa>1e-18&&na>1e-18)) return c;
            double nd=dot3(oldCr,newCr)/(oa*na); if(nd<minDot) return c;
            if(na < oa*minAreaRatio) return c;
            Vec3 un=oldCr*(1.0/oa);
            double pd=fabs(dot3(un,c.np-P[old.v[0]]));
            if(pd>maxPlane) return c;
            if(duplicateAfter(fid,a,b,d,ef[0],ef[1],rem)) return c;
            worstN=max(worstN,1.0-nd); worstP=max(worstP,pd/max(1e-12,maxPlane)); worstA=max(worstA,max(0.0,1.0-na/oa)); cnt++;
        }
        if(cnt==0) return c;
        c.ok=true;
        c.cost = c.newErr/maxErr + 2.0*worstP + 400.0*worstN + 0.05*worstA + (midpoint?0.05:0.0) + 0.0001*cnt;
        return c;
    }
    bool collapse(int u,int v){
        if(u==v||!av[u]||!av[v]) return false;
        int ef[2]={-1,-1}, opp[2]={-1,-1}; if(!edgeFaces(u,v,ef,opp)) return false; if(!linkOK(u,v,opp)) return false;
        Cand best=eval(u,v,ef,false), b2=eval(v,u,ef,false);
        if(b2.ok && (!best.ok || b2.cost<best.cost)) best=b2;
        if(N>300){ Cand b3=eval(u,v,ef,true), b4=eval(v,u,ef,true); if(b3.ok&&(!best.ok||b3.cost<best.cost))best=b3; if(b4.ok&&(!best.ok||b4.cost<best.cost))best=b4; }
        if(!best.ok) return false;
        int keep=best.keep, rem=best.rem;
        for(int t=0;t<2;t++) if(af[ef[t]]){ af[ef[t]]=0; aliveF--; }
        vector<int> touched;
        for(int fid:inc[rem]) if(af[fid]&&hasv(fid,rem)){
            for(int k=0;k<3;k++) if(F[fid].v[k]==rem) F[fid].v[k]=keep;
            touched.push_back(fid);
        }
        av[rem]=0; aliveV--; P[keep]=best.np; err[keep]=best.newErr;
        inc[keep].insert(inc[keep].end(),inc[rem].begin(),inc[rem].end()); vector<int>().swap(inc[rem]);
        sort(inc[keep].begin(),inc[keep].end()); inc[keep].erase(unique(inc[keep].begin(),inc[keep].end()),inc[keep].end());
        for(int fid:inc[keep]) if(af[fid]){
            const Face&f=F[fid]; pushEdge(f.v[0],f.v[1]); pushEdge(f.v[1],f.v[2]); pushEdge(f.v[2],f.v[0]);
        }
        return true;
    }
    double edgeCost(int u,int v){
        if(!av[u]||!av[v]) return 1e100;
        int ef[2],op[2]; double l2=norm2(P[u]-P[v]);
        if(edgeFaces(u,v,ef,op)){
            Vec3 n0=safe_unit(cross3(P[F[ef[0]].v[1]]-P[F[ef[0]].v[0]],P[F[ef[0]].v[2]]-P[F[ef[0]].v[0]]));
            Vec3 n1=safe_unit(cross3(P[F[ef[1]].v[1]]-P[F[ef[1]].v[0]],P[F[ef[1]].v[2]]-P[F[ef[1]].v[0]]));
            double cur=max(0.0,1.0-dot3(n0,n1)); return l2*(1.0+80.0*cur)+1e-18*(u+v);
        }
        return l2+1.0;
    }
    void pushEdge(int u,int v){ if(u==v||!av[u]||!av[v]) return; pq.push({edgeCost(u,v),u,v}); }
    void initializePQ(){
        for(int i=0;i<M;i++){ const Face&f=F[i]; pushEdge(f.v[0],f.v[1]); pushEdge(f.v[1],f.v[2]); pushEdge(f.v[2],f.v[0]); }
    }
    void run(){
        initializePQ();
        int target=(N<1000?max(4,N/2):(N<10000?max(24,N/35):(N<80000?max(160,N/22):max(900,N/16))));
        long long attempts=0, success=0;
        while(!pq.empty() && elapsed()<18.2 && aliveV>target){
            Node nd=pq.top(); pq.pop();
            if(!av[nd.u]||!av[nd.v]) continue;
            if(++attempts%200000==0 && elapsed()>17.2) break;
            if(collapse(nd.u,nd.v)) success++;
            if(success>0 && success%50000==0 && elapsed()>16.8) break;
        }
        if(elapsed()<18.4){
            double oldPlane=maxPlane, oldDot=minDot; maxPlane=min(maxPlane,0.0020*DIAG); minDot=max(minDot,0.9992);
            int lim=200000;
            while(!pq.empty() && elapsed()<18.7 && lim--){ Node nd=pq.top(); pq.pop(); if(av[nd.u]&&av[nd.v]) collapse(nd.u,nd.v); }
            maxPlane=oldPlane; minDot=oldDot;
        }
    }
    void build(vector<Vec3>&X,vector<Face>&OF){
        vector<int> mp(N,-1); X.clear(); OF.clear(); X.reserve(aliveV); OF.reserve(aliveF);
        for(int i=0;i<N;i++) if(av[i]){ mp[i]=X.size(); X.push_back(P[i]); }
        for(int i=0;i<M;i++) if(af[i]){
            int a=mp[F[i].v[0]],b=mp[F[i].v[1]],c=mp[F[i].v[2]];
            if(a>=0&&b>=0&&c>=0&&a!=b&&a!=c&&b!=c) OF.push_back(Face{{a,b,c}});
        }
    }
};

int main(){
    T0=chrono::steady_clock::now();
    readInput();
    vector<Vec3>X; vector<Face>OF;
    if(tryEllipsoidRoute(X,OF)){ writeMesh(X,OF); return 0; }
    if(tryCapsuleRoute(X,OF)){ writeMesh(X,OF); return 0; }
    if(tryBoxGridRoute(X,OF)){ writeMesh(X,OF); return 0; }
    if(tryLatRoute(X,OF)){ writeMesh(X,OF); return 0; }
    if(tryTorusRoute(X,OF)){ writeMesh(X,OF); return 0; }
    Simplifier s; s.run(); s.build(X,OF);
    if(X.empty()||OF.empty()){ writeMesh(V0,F0); return 0; }
    writeMesh(X,OF);
    return 0;
}
