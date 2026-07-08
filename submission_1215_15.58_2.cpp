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
static inline double norm3(const Vec3&a){return sqrt(norm2(a));}
static inline Vec3 unit3(Vec3 a){ double n=norm3(a); return n>1e-300? a*(1.0/n):Vec3{0,0,0}; }
static inline double clampd(double x,double l,double r){return x<l?l:(x>r?r:x);} 
struct Face{ int v[3]; };

static int N0=0,M0=0;
static vector<Vec3> origV;
static vector<Face> origF;
static Vec3 bbMin,bbMax;
static double diagLen=1.0, epsH=0.05;
static chrono::steady_clock::time_point startTime;
static inline double elapsed(){return chrono::duration<double>(chrono::steady_clock::now()-startTime).count();}

struct FastInput{
    vector<char> buf; char* p=nullptr;
    FastInput(){
        buf.reserve(1<<27);
        char tmp[1<<16]; size_t n;
        while((n=fread(tmp,1,sizeof(tmp),stdin))>0) buf.insert(buf.end(),tmp,tmp+n);
        buf.push_back('\0'); p=buf.data();
    }
    inline void skip(){while(*p==' '||*p=='\n'||*p=='\r'||*p=='\t') ++p;}
    long nextLong(){skip(); return strtol(p,&p,10);} 
    double nextDouble(){skip(); return strtod(p,&p);} 
    char nextChar(){skip(); return *p++;}
};

static void readInput(){
    FastInput in;
    N0=(int)in.nextLong(); M0=(int)in.nextLong();
    origV.resize(N0); origF.resize(M0);
    bbMin={1e100,1e100,1e100}; bbMax={-1e100,-1e100,-1e100};
    for(int i=0;i<N0;i++){
        (void)in.nextChar();
        Vec3 p{in.nextDouble(),in.nextDouble(),in.nextDouble()};
        origV[i]=p;
        bbMin.x=min(bbMin.x,p.x); bbMin.y=min(bbMin.y,p.y); bbMin.z=min(bbMin.z,p.z);
        bbMax.x=max(bbMax.x,p.x); bbMax.y=max(bbMax.y,p.y); bbMax.z=max(bbMax.z,p.z);
    }
    diagLen=norm3(bbMax-bbMin); if(!(diagLen>0)) diagLen=1.0; epsH=0.05*diagLen;
    for(int i=0;i<M0;i++){
        (void)in.nextChar();
        origF[i].v[0]=(int)in.nextLong()-1;
        origF[i].v[1]=(int)in.nextLong()-1;
        origF[i].v[2]=(int)in.nextLong()-1;
    }
}

static inline uint64_t edgeKey(int a,int b){ if(a>b) swap(a,b); return (uint64_t)(uint32_t)a<<32 | (uint32_t)b; }
static inline Vec3 faceCross(const vector<Vec3>&V,const Face&f){ return cross3(V[f.v[1]]-V[f.v[0]], V[f.v[2]]-V[f.v[0]]); }

struct HashGrid{
    double cell=1.0; Vec3 mn{}; unordered_map<long long, vector<int>> mp; const vector<Vec3>* pts=nullptr;
    static long long pack(int x,int y,int z){
        const long long B=1048576LL;
        return ((long long)(x+B)<<42) ^ ((long long)(y+B)<<21) ^ (long long)(z+B);
    }
    void build(const vector<Vec3>&p,double c){
        pts=&p; cell=max(c,1e-9); mn=bbMin; mn.x-=cell*2; mn.y-=cell*2; mn.z-=cell*2;
        mp.clear(); mp.reserve(p.size()*2+10);
        for(int i=0;i<(int)p.size();++i){
            int ix=(int)floor((p[i].x-mn.x)/cell), iy=(int)floor((p[i].y-mn.y)/cell), iz=(int)floor((p[i].z-mn.z)/cell);
            mp[pack(ix,iy,iz)].push_back(i);
        }
    }
    bool nearPoint(const Vec3&q,double r2) const{
        int ix=(int)floor((q.x-mn.x)/cell), iy=(int)floor((q.y-mn.y)/cell), iz=(int)floor((q.z-mn.z)/cell);
        for(int dz=-1;dz<=1;dz++) for(int dy=-1;dy<=1;dy++) for(int dx=-1;dx<=1;dx++){
            auto it=mp.find(pack(ix+dx,iy+dy,iz+dz)); if(it==mp.end()) continue;
            for(int id: it->second) if(norm2((*pts)[id]-q)<=r2) return true;
        }
        return false;
    }
};

static bool coverageOrigByCand(const vector<Vec3>&cand,double scale=0.997){
    if(cand.empty() || (int)cand.size()>N0) return false;
    double r=epsH*scale, r2=r*r*(1.0+1e-12);
    HashGrid g; g.build(cand,r);
    for(int i=0;i<N0;i++){
        if((i&32767)==0 && elapsed()>18.5) return false;
        if(!g.nearPoint(origV[i],r2)) return false;
    }
    return true;
}

static bool basicMeshOK(const vector<Vec3>&V,const vector<Face>&F){
    if(V.empty()||F.empty()||(int)V.size()>N0) return false;
    double minA2=max(1e-300,1e-28*diagLen*diagLen*diagLen*diagLen);
    unordered_map<uint64_t,int> cnt; cnt.reserve((size_t)F.size()*4+10);
    for(const Face&f:F){
        int a=f.v[0],b=f.v[1],c=f.v[2];
        if(a<0||b<0||c<0||a>=(int)V.size()||b>=(int)V.size()||c>=(int)V.size()) return false;
        if(a==b||a==c||b==c) return false;
        if(norm2(faceCross(V,f))<=minA2) return false;
        cnt[edgeKey(a,b)]++; cnt[edgeKey(b,c)]++; cnt[edgeKey(c,a)]++;
    }
    for(auto &kv:cnt) if(kv.second!=2) return false;
    return true;
}

struct RenderMap{ vector<float>d,nx,ny,nz; vector<unsigned char>fg; };
static unordered_map<int, array<RenderMap,6>> origRenderCache;

static inline void projectPoint(const Vec3&p,int view,int res,double&u,double&v,double&z){
    const double D=2.5; double f=800.0*(double)res/1024.0, c=0.5*(double)res; double sx,sy;
    if(view==0){ sx=p.y; sy=p.z; z=D-p.x; }
    else if(view==1){ sx=-p.y; sy=p.z; z=D+p.x; }
    else if(view==2){ sx=-p.x; sy=p.z; z=D-p.y; }
    else if(view==3){ sx=p.x; sy=p.z; z=D+p.y; }
    else if(view==4){ sx=p.x; sy=p.y; z=D-p.z; }
    else { sx=-p.x; sy=p.y; z=D+p.z; }
    u=f*sx/z+c; v=f*sy/z+c;
}
static void rasterTri(RenderMap&rm,int res,const Vec3&a,const Vec3&b,const Vec3&c,const Vec3&n,int view){
    double x0,y0,z0,x1,y1,z1,x2,y2,z2;
    projectPoint(a,view,res,x0,y0,z0); projectPoint(b,view,res,x1,y1,z1); projectPoint(c,view,res,x2,y2,z2);
    if(z0<=0||z1<=0||z2<=0) return;
    int xmin=max(0,(int)floor(min({x0,x1,x2})-0.5));
    int xmax=min(res-1,(int)ceil(max({x0,x1,x2})+0.5));
    int ymin=max(0,(int)floor(min({y0,y1,y2})-0.5));
    int ymax=min(res-1,(int)ceil(max({y0,y1,y2})+0.5));
    if(xmin>xmax||ymin>ymax) return;
    double den=(y1-y2)*(x0-x2)+(x2-x1)*(y0-y2); if(fabs(den)<1e-20) return;
    for(int yy=ymin; yy<=ymax; ++yy){ double py=yy+0.5;
        for(int xx=xmin; xx<=xmax; ++xx){ double px=xx+0.5;
            double w0=((y1-y2)*(px-x2)+(x2-x1)*(py-y2))/den;
            double w1=((y2-y0)*(px-x2)+(x0-x2)*(py-y2))/den;
            double w2=1.0-w0-w1;
            if(w0<-1e-9||w1<-1e-9||w2<-1e-9) continue;
            double zp=1.0/(w0/z0+w1/z1+w2/z2);
            int id=yy*res+xx;
            if(zp<rm.d[id]){rm.d[id]=(float)zp; rm.nx[id]=(float)n.x; rm.ny[id]=(float)n.y; rm.nz[id]=(float)n.z; rm.fg[id]=1;}
        }
    }
}
static RenderMap renderMesh(const vector<Vec3>&V,const vector<Face>&F,int res,int view){
    int S=res*res; RenderMap rm; rm.d.assign(S,255.0f); rm.nx.assign(S,0); rm.ny.assign(S,0); rm.nz.assign(S,0); rm.fg.assign(S,0);
    for(const Face&f:F){
        if(f.v[0]<0||f.v[1]<0||f.v[2]<0||f.v[0]>=(int)V.size()||f.v[1]>=(int)V.size()||f.v[2]>=(int)V.size()) continue;
        Vec3 cr=faceCross(V,f); double l=norm3(cr); if(!(l>0)) continue;
        rasterTri(rm,res,V[f.v[0]],V[f.v[1]],V[f.v[2]],cr*(1.0/l),view);
    }
    return rm;
}
static const array<RenderMap,6>& getOrigRenders(int res){
    auto it=origRenderCache.find(res); if(it!=origRenderCache.end()) return it->second;
    array<RenderMap,6> arr;
    for(int v=0;v<6;v++) arr[v]=renderMesh(origV,origF,res,v);
    auto ins=origRenderCache.emplace(res, std::move(arr)); return ins.first->second;
}
static double visualProxyScore(const vector<Vec3>&V,const vector<Face>&F,int res=96){
    const auto& O=getOrigRenders(res);
    double total=0; int views=0;
    for(int view=0; view<6; ++view){
        RenderMap R=renderMesh(V,F,res,view);
        const RenderMap&A=O[view]; double sc=0; int cnt=0; int S=res*res;
        for(int i=0;i<S;i++){
            if(!A.fg[i] && !R.fg[i]) continue;
            cnt++;
            if(!A.fg[i] || !R.fg[i]){ sc += 0.05; continue; }
            double nd=A.nx[i]*R.nx[i]+A.ny[i]*R.ny[i]+A.nz[i]*R.nz[i];
            double ns=clampd((nd+1.0)*0.5,0.0,1.0);
            double dz=fabs((double)A.d[i]-(double)R.d[i]);
            double ds=exp(-(dz*dz)/(0.018*0.018));
            sc += 0.58*ns + 0.42*ds;
        }
        if(cnt==0) return 0; total += sc/cnt; views++;
    }
    return views? total/views:0;
}

static void saveMesh(const vector<Vec3>&V,const vector<Face>&F){
    string out; out.reserve((size_t)V.size()*45+(size_t)F.size()*28+64);
    char buf[128];
    out.append(buf, snprintf(buf,sizeof(buf), "%d %d\n", (int)V.size(), (int)F.size()));
    for(const Vec3&p:V) out.append(buf, snprintf(buf,sizeof(buf), "v %.12g %.12g %.12g\n", p.x,p.y,p.z));
    for(const Face&f:F) out.append(buf, snprintf(buf,sizeof(buf), "f %d %d %d\n", f.v[0]+1,f.v[1]+1,f.v[2]+1));
    fwrite(out.data(),1,out.size(),stdout);
}

// ---------------- regular ordered-grid and radial candidates ----------------
static bool triLooksLocal(int a,int b,int c,int U,int Vc){
    int ids[3]={a,b,c}; int u[3],v[3]; for(int k=0;k<3;k++){u[k]=ids[k]/Vc; v[k]=ids[k]%Vc;}
    for(int base=0;base<3;base++){
        int mnU=0,mxU=0,mnV=0,mxV=0; bool ok=true;
        for(int k=0;k<3;k++){
            int du=(u[k]-u[base]+U)%U; if(du>U/2) du-=U;
            int dv=(v[k]-v[base]+Vc)%Vc; if(dv>Vc/2) dv-=Vc;
            if(du<-1||du>1||dv<-1||dv>1){ok=false; break;}
            mnU=min(mnU,du); mxU=max(mxU,du); mnV=min(mnV,dv); mxV=max(mxV,dv);
        }
        if(ok && mxU-mnU<=1 && mxV-mnV<=1 && (mxU-mnU+mxV-mnV)>=1) return true;
    }
    return false;
}
static double gridOrderScore(int U,int Vc){
    if(U<4||Vc<4||1LL*U*Vc!=N0) return 0;
    int stride=max(1,M0/25000), cnt=0, good=0;
    for(int i=0;i<M0;i+=stride){
        const Face&f=origF[i]; cnt++;
        if(triLooksLocal(f.v[0],f.v[1],f.v[2],U,Vc)) good++;
    }
    return cnt? (double)good/cnt:0;
}
static vector<int> circularIndices(int n,int k){
    vector<int>a; a.reserve(k);
    int last=-1;
    for(int i=0;i<k;i++){ int x=(int)floor((double)i*n/k+1e-9); if(x!=last){a.push_back(x); last=x;} }
    return a;
}
static bool torusParamCoverage(const vector<int>&Us,const vector<int>&Vs,int U,int Vc){
    double r2=(epsH*0.985)*(epsH*0.985);
    for(int i=0;i<N0;i++){
        int u=i/Vc, v=i%Vc;
        auto iu=lower_bound(Us.begin(),Us.end(),u); int cu[2];
        cu[0]=(iu==Us.end()?Us[0]:*iu); cu[1]=(iu==Us.begin()?Us.back():*(iu-1));
        auto iv=lower_bound(Vs.begin(),Vs.end(),v); int cv[2];
        cv[0]=(iv==Vs.end()?Vs[0]:*iv); cv[1]=(iv==Vs.begin()?Vs.back():*(iv-1));
        bool ok=false;
        for(int a=0;a<2&&!ok;a++) for(int b=0;b<2&&!ok;b++){
            int id=((cu[a]%U+U)%U)*Vc + ((cv[b]%Vc+Vc)%Vc);
            if(norm2(origV[id]-origV[i])<=r2) ok=true;
        }
        if(!ok) return false;
    }
    return true;
}
static void buildTorusCandidate(const vector<int>&Us,const vector<int>&Vs,int U,int Vc,int pat,vector<Vec3>&Vout,vector<Face>&Fout){
    int Ku=Us.size(), Kv=Vs.size(); Vout.clear(); Fout.clear(); Vout.reserve(Ku*Kv); Fout.reserve(Ku*Kv*2);
    auto id=[&](int i,int j){return ((i%Ku+Ku)%Ku)*Kv + ((j%Kv+Kv)%Kv);};
    for(int i=0;i<Ku;i++) for(int j=0;j<Kv;j++) Vout.push_back(origV[Us[i]*Vc+Vs[j]]);
    for(int i=0;i<Ku;i++) for(int j=0;j<Kv;j++){
        int a=id(i,j), b=id(i+1,j), c=id(i,j+1), d=id(i+1,j+1);
        if(pat==0){ Fout.push_back({{a,b,c}}); Fout.push_back({{b,d,c}}); }
        else if(pat==1){ Fout.push_back({{a,c,b}}); Fout.push_back({{b,c,d}}); }
        else if(pat==2){ Fout.push_back({{a,b,d}}); Fout.push_back({{a,d,c}}); }
        else { Fout.push_back({{a,d,b}}); Fout.push_back({{a,c,d}}); }
    }
}
static bool tryOrderedTorus(vector<Vec3>&bestV, vector<Face>&bestF){
    if(N0<3000 || elapsed()>2.0) return false;
    if(llabs((long long)M0-2LL*N0)>8) return false;
    vector<pair<int,int>> cand;
    for(int d=4; (long long)d*d<=N0; ++d) if(N0%d==0){ cand.push_back({d,N0/d}); cand.push_back({N0/d,d}); }
    double bestS=0; int bestU=0,bestW=0;
    for(auto [U,Vc]:cand){
        if(U<8||Vc<8) continue;
        double asp=(double)max(U,Vc)/min(U,Vc); if(asp>100) continue;
        double s=gridOrderScore(U,Vc);
        if(s>bestS){bestS=s; bestU=U; bestW=Vc;}
    }
    if(bestS<0.985) return false;
    int U=bestU,Vc=bestW;
    vector<double> ratios={0.035,0.05,0.07,0.09,0.12,0.16,0.22,0.30};
    double bestProxy=-1; vector<Vec3> tv; vector<Face> tf;
    for(double ratio: ratios){
        if(elapsed()>7.5) break;
        int target=max(16,(int)ceil(N0*ratio));
        int Ku=max(4,min(U,(int)round(sqrt(target*(double)U/Vc))));
        int Kv=max(4,min(Vc,(int)ceil((double)target/Ku)));
        Ku=min(Ku,U); Kv=min(Kv,Vc);
        vector<int> Us=circularIndices(U,Ku), Vs=circularIndices(Vc,Kv);
        if((int)Us.size()<4||(int)Vs.size()<4) continue;
        if(!torusParamCoverage(Us,Vs,U,Vc)) continue;
        for(int pat=0;pat<4;pat++){
            buildTorusCandidate(Us,Vs,U,Vc,pat,tv,tf);
            if(!basicMeshOK(tv,tf)) continue;
            double proxy=visualProxyScore(tv,tf,88);
            if(proxy>bestProxy){bestProxy=proxy; bestV=tv; bestF=tf;}
        }
        if(bestProxy>0.91) break;
    }
    return bestProxy>0.885 && !bestV.empty();
}

static bool buildRadialCandidate(int lat,int lon,vector<Vec3>&Vout,vector<Face>&Fout){
    const double PI=acos(-1.0);
    int rings=lat-1; if(rings<2||lon<6) return false;
    int top=-1,bot=-1; double zmax=-1e100,zmin=1e100;
    vector<int> bin(rings*lon,-1); vector<double> bscore(rings*lon,-1e100);
    vector<Vec3> dirs(rings*lon);
    for(int i=1;i<lat;i++){
        double phi=PI*i/lat; double sp=sin(phi), cp=cos(phi);
        for(int j=0;j<lon;j++){ double th=2*PI*j/lon; dirs[(i-1)*lon+j]={cos(th)*sp,sin(th)*sp,cp}; }
    }
    for(int id=0;id<N0;id++){
        Vec3 p=origV[id]; double r=norm3(p); if(r<1e-12) continue; Vec3 u=p/r;
        if(p.z>zmax){zmax=p.z; top=id;} if(p.z<zmin){zmin=p.z; bot=id;}
        double phi=acos(clampd(u.z,-1.0,1.0)); int ii=(int)floor(phi/PI*lat+0.5); ii=max(1,min(lat-1,ii));
        double th=atan2(u.y,u.x); if(th<0) th+=2*PI; int jj=(int)floor(th/(2*PI)*lon+0.5); if(jj>=lon) jj-=lon;
        int b=(ii-1)*lon+jj; double sc=dot3(u,dirs[b]);
        if(sc>bscore[b]){bscore[b]=sc; bin[b]=id;}
    }
    if(top<0||bot<0||top==bot) return false;
    for(int i=0;i<rings;i++) for(int j=0;j<lon;j++) if(bin[i*lon+j]<0){
        int best=-1; double bs=-1e100;
        Vec3 d=dirs[i*lon+j];
        int maxRad=max(rings,lon/2);
        for(int rr=1;rr<=maxRad && best<0;rr++){
            for(int di=-rr;di<=rr;di++) for(int dj=-rr;dj<=rr;dj++){
                if(abs(di)!=rr && abs(dj)!=rr) continue;
                int ni=i+di; if(ni<0||ni>=rings) continue;
                int nj=(j+dj)%lon; if(nj<0) nj+=lon;
                int id=bin[ni*lon+nj]; if(id<0) continue;
                double sc=dot3(unit3(origV[id]),d);
                if(sc>bs){bs=sc; best=id;}
            }
        }
        if(best<0) return false;
        bin[i*lon+j]=best;
    }
    Vout.clear(); Fout.clear(); Vout.reserve(2+rings*lon); Fout.reserve(2*lat*lon);
    Vout.push_back(origV[top]);
    for(int i=0;i<rings;i++) for(int j=0;j<lon;j++) Vout.push_back(origV[bin[i*lon+j]]);
    int south=(int)Vout.size(); Vout.push_back(origV[bot]);
    auto rid=[&](int i,int j){return 1+i*lon+((j%lon+lon)%lon);};
    for(int j=0;j<lon;j++){ Fout.push_back({{0,rid(0,j+1),rid(0,j)}}); }
    for(int i=0;i<rings-1;i++) for(int j=0;j<lon;j++){
        int a=rid(i,j), b=rid(i+1,j), c=rid(i,j+1), d=rid(i+1,j+1);
        Fout.push_back({{a,b,c}}); Fout.push_back({{b,d,c}});
    }
    for(int j=0;j<lon;j++){ Fout.push_back({{south,rid(rings-1,j),rid(rings-1,j+1)}}); }
    return true;
}
static bool tryRadialGenus0(vector<Vec3>&bestV, vector<Face>&bestF){
    if(N0<6000 || elapsed()>8.0) return false;
    if(llabs((long long)M0 - (2LL*N0-4))>16) return false;
    vector<double> ratios={0.045,0.07,0.10,0.14,0.20,0.28};
    double bestProxy=-1; vector<Vec3> rv; vector<Face> rf;
    for(double ratio:ratios){
        if(elapsed()>11.5) break;
        int target=max(100,(int)ceil(N0*ratio));
        int lat=max(8,(int)round(sqrt(target/2.0))); int lon=max(16,2*lat);
        if(lat*lon+2>N0) continue;
        if(!buildRadialCandidate(lat,lon,rv,rf)) continue;
        if(!basicMeshOK(rv,rf)) continue;
        if(!coverageOrigByCand(rv,0.985)) continue;
        vector<Face> origFaces=rf;
        double p0=visualProxyScore(rv,rf,88);
        for(auto &f:rf) swap(f.v[1],f.v[2]);
        double p1=visualProxyScore(rv,rf,88);
        if(p0>=p1) rf=origFaces; else p0=p1;
        if(p0>bestProxy){bestProxy=p0; bestV=rv; bestF=rf;}
        if(bestProxy>0.92) break;
    }
    return bestProxy>0.895 && !bestV.empty();
}

// ---------------- QEM/edge-collapse fallback ----------------
struct Quadric{ double q[10]; Quadric(){memset(q,0,sizeof(q));} };
static inline Quadric qadd(const Quadric&a,const Quadric&b){Quadric r; for(int i=0;i<10;i++) r.q[i]=a.q[i]+b.q[i]; return r;}
static inline void qacc(Quadric&Q,double a,double b,double c,double d,double w){
    Q.q[0]+=w*a*a; Q.q[1]+=w*a*b; Q.q[2]+=w*a*c; Q.q[3]+=w*a*d;
    Q.q[4]+=w*b*b; Q.q[5]+=w*b*c; Q.q[6]+=w*b*d; Q.q[7]+=w*c*c; Q.q[8]+=w*c*d; Q.q[9]+=w*d*d;
}
static inline double qeval(const Quadric&Q,const Vec3&p){
    double x=p.x,y=p.y,z=p.z;
    return Q.q[0]*x*x+2*Q.q[1]*x*y+2*Q.q[2]*x*z+2*Q.q[3]*x+Q.q[4]*y*y+2*Q.q[5]*y*z+2*Q.q[6]*y+Q.q[7]*z*z+2*Q.q[8]*z+Q.q[9];
}
static bool qsolve(const Quadric&Q,Vec3&out){
    double a00=Q.q[0], a01=Q.q[1], a02=Q.q[2], b0=-Q.q[3];
    double a10=Q.q[1], a11=Q.q[4], a12=Q.q[5], b1=-Q.q[6];
    double a20=Q.q[2], a21=Q.q[5], a22=Q.q[7], b2=-Q.q[8];
    double det=a00*(a11*a22-a12*a21)-a01*(a10*a22-a12*a20)+a02*(a10*a21-a11*a20);
    if(fabs(det)<1e-14) return false;
    double dx=b0*(a11*a22-a12*a21)-a01*(b1*a22-a12*b2)+a02*(b1*a21-a11*b2);
    double dy=a00*(b1*a22-a12*b2)-b0*(a10*a22-a12*a20)+a02*(a10*b2-b1*a20);
    double dz=a00*(a11*b2-b1*a21)-a01*(a10*b2-b1*a20)+b0*(a10*a21-a11*a20);
    out={dx/det,dy/det,dz/det};
    if(!isfinite(out.x)||!isfinite(out.y)||!isfinite(out.z)) return false;
    return true;
}

static vector<Vec3> P;
static vector<Face> WF;
static vector<unsigned char> aliveV, aliveF;
static vector<double> coverR;
static vector<Quadric> Qs;
static vector<vector<int>> adj;
static vector<int> markA,markB; static int stampA=1, stampB=1;
static int aliveN=0, aliveM=0;

static inline bool faceHas(const Face&f,int v){return f.v[0]==v||f.v[1]==v||f.v[2]==v;}
static inline int thirdOf(const Face&f,int a,int b){for(int k=0;k<3;k++){int x=f.v[k]; if(x!=a&&x!=b) return x;} return -1;}
static inline Vec3 curCross(const Face&f){return cross3(P[f.v[1]]-P[f.v[0]],P[f.v[2]]-P[f.v[0]]);} 

static void initWork(){
    P=origV; WF=origF; aliveV.assign(N0,1); aliveF.assign(M0,1); coverR.assign(N0,0.0); Qs.assign(N0,Quadric()); aliveN=N0; aliveM=M0;
    for(const Face&f:origF){
        Vec3 cr=faceCross(origV,f); double l=norm3(cr); if(!(l>0)) continue;
        Vec3 n=cr*(1.0/l); double d=-dot3(n,origV[f.v[0]]); double w=max(1e-12,l);
        qacc(Qs[f.v[0]],n.x,n.y,n.z,d,w); qacc(Qs[f.v[1]],n.x,n.y,n.z,d,w); qacc(Qs[f.v[2]],n.x,n.y,n.z,d,w);
    }
}
static void rebuildAdj(){
    vector<int> deg(N0,0);
    for(int i=0;i<M0;i++) if(aliveF[i]){
        Face f=WF[i];
        if(f.v[0]<0||f.v[1]<0||f.v[2]<0||f.v[0]>=N0||f.v[1]>=N0||f.v[2]>=N0||!aliveV[f.v[0]]||!aliveV[f.v[1]]||!aliveV[f.v[2]]||f.v[0]==f.v[1]||f.v[0]==f.v[2]||f.v[1]==f.v[2]){aliveF[i]=0; aliveM--; continue;}
        deg[f.v[0]]++; deg[f.v[1]]++; deg[f.v[2]]++;
    }
    adj.assign(N0,{}); for(int i=0;i<N0;i++) if(aliveV[i]) adj[i].reserve(deg[i]);
    for(int i=0;i<M0;i++) if(aliveF[i]){Face f=WF[i]; adj[f.v[0]].push_back(i); adj[f.v[1]].push_back(i); adj[f.v[2]].push_back(i);}    
    markA.assign(N0,0); markB.assign(N0,0); stampA=stampB=1;
}
static bool findEdgeFaces(int u,int v,int ef[2],int opp[2]){
    int cnt=0; if(u<0||v<0||u>=N0||v>=N0) return false;
    for(int fid:adj[u]) if(aliveF[fid]&&faceHas(WF[fid],u)&&faceHas(WF[fid],v)){
        if(cnt>=2) return false; ef[cnt]=fid; opp[cnt]=thirdOf(WF[fid],u,v); cnt++;
    }
    return cnt==2 && opp[0]>=0 && opp[1]>=0 && opp[0]!=opp[1];
}
static bool linkOK(int u,int v,const int opp[2]){
    if(++stampA>2000000000){fill(markA.begin(),markA.end(),0);stampA=1;}
    if(++stampB>2000000000){fill(markB.begin(),markB.end(),0);stampB=1;}
    for(int fid:adj[u]) if(aliveF[fid]&&faceHas(WF[fid],u)){
        const Face&f=WF[fid]; for(int k=0;k<3;k++){int x=f.v[k]; if(x!=u&&x!=v) markA[x]=stampA;}
    }
    int common=0;
    for(int fid:adj[v]) if(aliveF[fid]&&faceHas(WF[fid],v)){
        const Face&f=WF[fid]; for(int k=0;k<3;k++){int x=f.v[k]; if(x==u||x==v) continue; if(markA[x]!=stampA) continue; if(x!=opp[0]&&x!=opp[1]) return false; if(markB[x]!=stampB){markB[x]=stampB; common++;}}
    }
    return common==2 && markB[opp[0]]==stampB && markB[opp[1]]==stampB;
}
static bool sameUnordered(const Face&f,int a,int b,int c){int x[3]={f.v[0],f.v[1],f.v[2]}, y[3]={a,b,c}; sort(x,x+3); sort(y,y+3); return x[0]==y[0]&&x[1]==y[1]&&x[2]==y[2];}
static bool duplicateAfter(int rem,int fid,int ef0,int ef1,int a,int b,int c){
    int probe=a; if(adj[b].size()<adj[probe].size()) probe=b; if(adj[c].size()<adj[probe].size()) probe=c;
    for(int g:adj[probe]) if(aliveF[g]&&g!=fid&&g!=ef0&&g!=ef1&&!faceHas(WF[g],rem)&&sameUnordered(WF[g],a,b,c)) return true;
    return false;
}
struct Params{ double coverScale, planeScale, minCos, areaScale; bool allowMove; int rounds; double timeLimit; };
struct Eval{ bool ok=false; int keep=-1,rem=-1; Vec3 pos{}; double newR=0,cost=1e100; };
static Vec3 posOf(int id,int keep,const Vec3&pos){return id==keep?pos:P[id];}
static Eval evalCollapse(int keep,int rem,const int ef[2],const Params&par,const Vec3&pos){
    Eval e; e.keep=keep; e.rem=rem; e.pos=pos;
    double nr=max(coverR[keep]+norm3(P[keep]-pos), coverR[rem]+norm3(P[rem]-pos));
    if(nr>par.coverScale*diagLen+1e-12) return e; e.newR=nr;
    static vector<int> affected; affected.clear(); affected.reserve(adj[keep].size()+adj[rem].size());
    for(int fid:adj[keep]) if(aliveF[fid]&&fid!=ef[0]&&fid!=ef[1]&&faceHas(WF[fid],keep)) affected.push_back(fid);
    for(int fid:adj[rem]) if(aliveF[fid]&&fid!=ef[0]&&fid!=ef[1]&&faceHas(WF[fid],rem)) affected.push_back(fid);
    sort(affected.begin(),affected.end()); affected.erase(unique(affected.begin(),affected.end()),affected.end());
    if(affected.empty()) return e;
    double minA2=max(1e-300,par.areaScale*diagLen*diagLen*diagLen*diagLen);
    double maxAng=0,maxPlane=0,maxAreaLoss=0; int changed=0;
    for(int fid:affected){
        Face old=WF[fid], nf=old; bool hadRem=false;
        for(int k=0;k<3;k++) if(nf.v[k]==rem){nf.v[k]=keep; hadRem=true;}
        if(nf.v[0]==nf.v[1]||nf.v[0]==nf.v[2]||nf.v[1]==nf.v[2]) return e;
        Vec3 co=curCross(old);
        Vec3 a=posOf(nf.v[0],keep,pos), b=posOf(nf.v[1],keep,pos), c=posOf(nf.v[2],keep,pos);
        Vec3 cn=cross3(b-a,c-a); double ao=norm3(co), an=norm3(cn);
        if(!(ao>0)||!(an>0)||norm2(cn)<=minA2) return e;
        double nd=clampd(dot3(co,cn)/(ao*an),-1.0,1.0); if(nd<par.minCos) return e;
        Vec3 no=co*(1.0/ao); double pd=0;
        Vec3 op0=P[old.v[0]];
        pd=max(pd,fabs(dot3(no,a-op0))); pd=max(pd,fabs(dot3(no,b-op0))); pd=max(pd,fabs(dot3(no,c-op0)));
        if(pd>par.planeScale*diagLen) return e;
        if(hadRem && duplicateAfter(rem,fid,ef[0],ef[1],nf.v[0],nf.v[1],nf.v[2])) return e;
        maxAng=max(maxAng,1.0-nd); maxPlane=max(maxPlane,pd/(par.planeScale*diagLen+1e-30)); maxAreaLoss=max(maxAreaLoss,max(0.0,1.0-an/ao)); changed++;
    }
    Quadric q=qadd(Qs[keep],Qs[rem]); double qe=max(0.0,qeval(q,pos))/(diagLen*diagLen+1e-30);
    e.cost=qe + 120.0*maxAng + 0.5*maxPlane + 0.05*maxAreaLoss + 0.02*(nr/(par.coverScale*diagLen+1e-30)) + 1e-5*changed;
    e.ok=true; return e;
}
static void applyCollapse(const Eval&e,const int ef[2]){
    int keep=e.keep, rem=e.rem;
    for(int i=0;i<2;i++) if(ef[i]>=0&&aliveF[ef[i]]){aliveF[ef[i]]=0; aliveM--;}
    for(int fid:adj[rem]) if(aliveF[fid]&&faceHas(WF[fid],rem)) for(int k=0;k<3;k++) if(WF[fid].v[k]==rem) WF[fid].v[k]=keep;
    aliveV[rem]=0; aliveN--; P[keep]=e.pos; coverR[keep]=e.newR; Qs[keep]=qadd(Qs[keep],Qs[rem]);
    vector<int> merged; merged.reserve(adj[keep].size()+adj[rem].size());
    for(int fid:adj[keep]) if(aliveF[fid]&&faceHas(WF[fid],keep)) merged.push_back(fid);
    for(int fid:adj[rem]) if(aliveF[fid]&&faceHas(WF[fid],keep)) merged.push_back(fid);
    sort(merged.begin(),merged.end()); merged.erase(unique(merged.begin(),merged.end()),merged.end()); adj[keep].swap(merged); vector<int>().swap(adj[rem]);
}
static bool tryEdge(int u,int v,const Params&par){
    if(u==v||u<0||v<0||u>=N0||v>=N0||!aliveV[u]||!aliveV[v]) return false;
    int ef[2]={-1,-1}, opp[2]={-1,-1}; if(!findEdgeFaces(u,v,ef,opp)) return false; if(!linkOK(u,v,opp)) return false;
    Quadric q=qadd(Qs[u],Qs[v]); vector<pair<int,Vec3>> cand;
    cand.push_back({0,P[u]}); cand.push_back({1,P[v]});
    if(par.allowMove){ cand.push_back({0,(P[u]+P[v])*0.5}); Vec3 opt; if(qsolve(q,opt)) cand.push_back({0,opt}); }
    Eval best;
    for(auto &cp:cand){
        if(cp.first==0){ Eval e=evalCollapse(u,v,ef,par,cp.second); if(e.ok&&e.cost<best.cost) best=e; }
        else { Eval e=evalCollapse(v,u,ef,par,cp.second); if(e.ok&&e.cost<best.cost) best=e; }
    }
    if(!best.ok) return false; applyCollapse(best,ef); return true;
}
struct EdgeRec{ float c; uint64_t k; };
static bool collapseRound(const Params&par,int targetN){
    vector<uint64_t> keys; keys.reserve((size_t)aliveM*3+16);
    for(int i=0;i<M0;i++) if(aliveF[i]){
        const Face&f=WF[i];
        if(aliveV[f.v[0]]&&aliveV[f.v[1]]&&aliveV[f.v[2]]){keys.push_back(edgeKey(f.v[0],f.v[1])); keys.push_back(edgeKey(f.v[1],f.v[2])); keys.push_back(edgeKey(f.v[2],f.v[0]));}
    }
    sort(keys.begin(),keys.end()); keys.erase(unique(keys.begin(),keys.end()),keys.end());
    vector<EdgeRec> edges; edges.reserve(keys.size()); double e2=(par.coverScale*diagLen)*(par.coverScale*diagLen)*4.05;
    for(uint64_t k:keys){ int a=(int)(k>>32), b=(int)(uint32_t)k; if(aliveV[a]&&aliveV[b]){ double l2=norm2(P[a]-P[b]); if(l2<=e2 || par.allowMove) edges.push_back({(float)l2,k}); } }
    sort(edges.begin(),edges.end(),[](const EdgeRec&a,const EdgeRec&b){return a.c<b.c;});
    int hit=0;
    for(size_t i=0;i<edges.size() && elapsed()<par.timeLimit && aliveN>targetN;i++){
        int a=(int)(edges[i].k>>32), b=(int)(uint32_t)edges[i].k;
        if(a>=0&&b>=0&&a<N0&&b<N0&&aliveV[a]&&aliveV[b]) if(tryEdge(a,b,par)) hit++;
    }
    return hit>0;
}
static double analyzeSmoothness(){
    int stride=max(1,M0/40000), cnt=0, smooth=0;
    for(int i=0;i<M0 && cnt<80000;i+=stride){
        const Face&f=origF[i]; Vec3 n0=unit3(faceCross(origV,f));
        for(int e=0;e<3;e++){
            int a=f.v[e], b=f.v[(e+1)%3];
            // approximate by scanning small adjacency after it is built
            int ef[2],op[2]; if(findEdgeFaces(a,b,ef,op)){
                Vec3 n1=unit3(curCross(WF[ef[0]])), n2=unit3(curCross(WF[ef[1]])); double d=dot3(n1,n2); cnt++; if(d>0.94) smooth++;
            }
        }
    }
    return cnt? (double)smooth/cnt:0.5;
}
static void runEdgeCollapse(){
    initWork(); rebuildAdj();
    double sm=analyzeSmoothness();
    int target=max(8,(int)(N0*(sm>0.75?0.045:(sm>0.45?0.065:0.09))));
    vector<Params> phases;
    phases.push_back({0.04945,0.014,0.965,1e-28,true,2,16.8});
    phases.push_back({0.04965,0.026,0.925,1e-29,true,2,18.4});
    if(sm>0.55) phases.push_back({0.04975,0.040,0.865,1e-30,true,2,19.5});
    else phases.push_back({0.04970,0.030,0.905,1e-30,true,1,19.2});
    for(const Params&ph:phases){
        for(int r=0;r<ph.rounds && elapsed()<ph.timeLimit && aliveN>target;r++){
            bool any=collapseRound(ph,target);
            rebuildAdj();
            if(!any) break;
        }
    }
}
static void compactWork(vector<Vec3>&Vout,vector<Face>&Fout){
    vector<int> mp(N0,-1); Vout.clear(); Fout.clear(); Vout.reserve(aliveN); Fout.reserve(aliveM);
    for(int i=0;i<N0;i++) if(aliveV[i]){mp[i]=(int)Vout.size(); Vout.push_back(P[i]);}
    for(int i=0;i<M0;i++) if(aliveF[i]){
        Face f=WF[i]; int a=mp[f.v[0]],b=mp[f.v[1]],c=mp[f.v[2]];
        if(a>=0&&b>=0&&c>=0&&a!=b&&a!=c&&b!=c) Fout.push_back({{a,b,c}});
    }
}

int main(){
    startTime=chrono::steady_clock::now();
    readInput();
    vector<Vec3> candV; vector<Face> candF;
    bool installed=false;
    if(tryOrderedTorus(candV,candF)) installed=true;
    if(!installed && tryRadialGenus0(candV,candF)) installed=true;
    if(!installed){ runEdgeCollapse(); compactWork(candV,candF); }
    if(!basicMeshOK(candV,candF) || (int)candV.size()>N0){ candV=origV; candF=origF; }
    saveMesh(candV,candF);
    return 0;
}
