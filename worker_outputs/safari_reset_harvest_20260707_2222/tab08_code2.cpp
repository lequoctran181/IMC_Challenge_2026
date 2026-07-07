#include <bits/stdc++.h>
using namespace std;

struct Vec3{
    double x,y,z;
    Vec3(double X=0,double Y=0,double Z=0):x(X),y(Y),z(Z){}
    Vec3 operator+(const Vec3& o)const{return Vec3(x+o.x,y+o.y,z+o.z);}
    Vec3 operator-(const Vec3& o)const{return Vec3(x-o.x,y-o.y,z-o.z);}
    Vec3 operator*(double s)const{return Vec3(x*s,y*s,z*s);}
    Vec3 operator/(double s)const{return Vec3(x/s,y/s,z/s);}
    Vec3& operator+=(const Vec3& o){x+=o.x;y+=o.y;z+=o.z;return *this;}
};
static inline double dotv(const Vec3&a,const Vec3&b){return a.x*b.x+a.y*b.y+a.z*b.z;}
static inline Vec3 crossv(const Vec3&a,const Vec3&b){return Vec3(a.y*b.z-a.z*b.y,a.z*b.x-a.x*b.z,a.x*b.y-a.y*b.x);}
static inline double norm2(const Vec3&a){return dotv(a,a);}
static inline double normv(const Vec3&a){return sqrt(norm2(a));}
static inline Vec3 normalize(Vec3 a){ double l=normv(a); return l>1e-300? a/l : Vec3(0,0,0); }

struct Face{int v[3]; unsigned char active; Vec3 n;};
struct Tri{int a,b,c;};

static int N=0,M=0;
static vector<Vec3> Orig;
static vector<Tri> OrigF;
static Vec3 g_mn,g_mx;
static double g_diag=1.0, g_haus=0.05;
static chrono::steady_clock::time_point START;
static inline double elapsed(){return chrono::duration<double>(chrono::steady_clock::now()-START).count();}

struct FastInput{
    vector<char> buf; char* p=nullptr;
    FastInput(){ char tmp[1<<16]; size_t n; while((n=fread(tmp,1,sizeof(tmp),stdin))>0) buf.insert(buf.end(),tmp,tmp+n); buf.push_back('\0'); p=buf.data(); }
    inline void skip(){ while(*p==' '||*p=='\n'||*p=='\r'||*p=='\t') ++p; }
    int nextInt(){ skip(); int s=1; if(*p=='-') s=-1,++p; int x=0; while(*p>='0'&&*p<='9') x=x*10+(*p++-'0'); return x*s; }
    double nextDouble(){ skip(); char* e; double x=strtod(p,&e); p=e; return x; }
    char nextChar(){ skip(); return *p++; }
};

static void loadMesh(){
    FastInput in; N=in.nextInt(); M=in.nextInt();
    Orig.resize(N); OrigF.resize(M);
    g_mn=Vec3(1e100,1e100,1e100); g_mx=Vec3(-1e100,-1e100,-1e100);
    for(int i=0;i<N;i++){
        in.nextChar(); double x=in.nextDouble(),y=in.nextDouble(),z=in.nextDouble();
        Orig[i]=Vec3(x,y,z);
        g_mn.x=min(g_mn.x,x); g_mn.y=min(g_mn.y,y); g_mn.z=min(g_mn.z,z);
        g_mx.x=max(g_mx.x,x); g_mx.y=max(g_mx.y,y); g_mx.z=max(g_mx.z,z);
    }
    for(int i=0;i<M;i++){
        in.nextChar(); int a=in.nextInt()-1,b=in.nextInt()-1,c=in.nextInt()-1;
        OrigF[i]={a,b,c};
    }
    g_diag=normv(g_mx-g_mn); if(!(g_diag>0)) g_diag=1.0;
    g_haus=0.05*g_diag*0.999;
}

static uint64_t ekey(int a,int b){ if(a>b) swap(a,b); return (uint64_t)(uint32_t)a<<32 | (uint32_t)b; }
static uint64_t fkey3(int a,int b,int c){ if(a>b) swap(a,b); if(b>c) swap(b,c); if(a>b) swap(a,b); return ((uint64_t)(uint32_t)a<<42) ^ ((uint64_t)(uint32_t)b<<21) ^ (uint32_t)c; }

struct HashGrid{
    Vec3 mn; double cell, r2; int nx,ny,nz; vector<vector<int>> buckets; const vector<Vec3>* pts;
    int clampi(int v,int n)const{return v<0?0:(v>=n?n-1:v);}
    int ix(double x)const{return clampi((int)floor((x-mn.x)/cell),nx);}
    int iy(double y)const{return clampi((int)floor((y-mn.y)/cell),ny);}
    int iz(double z)const{return clampi((int)floor((z-mn.z)/cell),nz);}
    int key(int x,int y,int z)const{return (z*ny+y)*nx+x;}
    void build(const vector<Vec3>& P,double R){
        pts=&P; r2=R*R; cell=max(R,1e-7);
        if(P.empty()){nx=ny=nz=1;buckets.assign(1,{});return;}
        mn=P[0]; Vec3 mx=P[0];
        for(auto&p:P){mn.x=min(mn.x,p.x);mn.y=min(mn.y,p.y);mn.z=min(mn.z,p.z);mx.x=max(mx.x,p.x);mx.y=max(mx.y,p.y);mx.z=max(mx.z,p.z);}
        mn.x-=1e-9;mn.y-=1e-9;mn.z-=1e-9;
        nx=max(1,(int)((mx.x-mn.x)/cell)+1); ny=max(1,(int)((mx.y-mn.y)/cell)+1); nz=max(1,(int)((mx.z-mn.z)/cell)+1);
        long long tot=1LL*nx*ny*nz;
        if(tot>2000000){ nx=ny=nz=1; cell=max({mx.x-mn.x,mx.y-mn.y,mx.z-mn.z,R})+1e-6; }
        buckets.assign((size_t)nx*ny*nz,{});
        for(int i=0;i<(int)P.size();++i) buckets[key(ix(P[i].x),iy(P[i].y),iz(P[i].z))].push_back(i);
    }
    bool nearPoint(const Vec3&p)const{
        int X=ix(p.x),Y=iy(p.y),Z=iz(p.z);
        for(int z=Z-1;z<=Z+1;z++) if(z>=0&&z<nz)
        for(int y=Y-1;y<=Y+1;y++) if(y>=0&&y<ny)
        for(int x=X-1;x<=X+1;x++) if(x>=0&&x<nx){
            for(int id:buckets[key(x,y,z)]) if(norm2((*pts)[id]-p)<=r2) return true;
        }
        return false;
    }
};

static bool vertexHausdorffOK(const vector<Vec3>& V){
    if(V.empty() || (int)V.size()>N) return false;
    HashGrid cand, orig; cand.build(V,g_haus); orig.build(Orig,g_haus);
    for(const Vec3&p:Orig) if(!cand.nearPoint(p)) return false;
    for(const Vec3&p:V) if(!orig.nearPoint(p)) return false;
    return true;
}

static bool manifoldOK(const vector<Vec3>&V,const vector<Tri>&F){
    if(V.empty()||F.empty()||(int)V.size()>N) return false;
    const double eps2=max(1e-30,1e-24*g_diag*g_diag);
    vector<uint64_t> edges; edges.reserve((size_t)F.size()*3);
    vector<uint64_t> faces; faces.reserve(F.size());
    for(const auto&t:F){
        int a=t.a,b=t.b,c=t.c;
        if(a<0||b<0||c<0||a>=(int)V.size()||b>=(int)V.size()||c>=(int)V.size()) return false;
        if(a==b||b==c||c==a) return false;
        if(norm2(crossv(V[b]-V[a],V[c]-V[a]))<=eps2) return false;
        faces.push_back(fkey3(a,b,c));
        edges.push_back(ekey(a,b)); edges.push_back(ekey(b,c)); edges.push_back(ekey(c,a));
    }
    sort(faces.begin(),faces.end());
    if(adjacent_find(faces.begin(),faces.end())!=faces.end()) return false;
    sort(edges.begin(),edges.end());
    for(size_t i=0;i<edges.size();){
        size_t j=i+1; while(j<edges.size()&&edges[j]==edges[i]) ++j;
        if(j-i!=2) return false;
        i=j;
    }
    return true;
}

struct RenderMaps{ int R=0; vector<float> depth; vector<float> norm; };
static inline void projectPoint(const Vec3&p,int view,int R,double&u,double&v,double&dep){
    double f=800.0*(double)R/1024.0, c=0.5*(double)R; double sx,sy; const double D=2.5;
    if(view==0){sx=p.y; sy=p.z; dep=D-p.x;}
    else if(view==1){sx=-p.y; sy=p.z; dep=D+p.x;}
    else if(view==2){sx=-p.x; sy=p.z; dep=D-p.y;}
    else if(view==3){sx=p.x; sy=p.z; dep=D+p.y;}
    else if(view==4){sx=p.x; sy=p.y; dep=D-p.z;}
    else {sx=-p.x; sy=p.y; dep=D+p.z;}
    u=f*sx/dep+c; v=f*sy/dep+c;
}
static void renderMesh(const vector<Vec3>&V,const vector<Tri>&F,RenderMaps&rm,int R){
    const int PIX=R*R; rm.R=R; rm.depth.assign((size_t)6*PIX,255.f); rm.norm.assign((size_t)6*PIX*3,127.5f);
    int nv=(int)V.size(), nf=(int)F.size(); vector<float> U(nv),W(nv),Z(nv); vector<Vec3> fn(nf);
    for(int i=0;i<nf;i++){ const auto&t=F[i]; Vec3 cr=crossv(V[t.b]-V[t.a],V[t.c]-V[t.a]); double l=normv(cr); fn[i]=l>1e-300? cr/l:Vec3(); }
    for(int view=0;view<6;view++){
        for(int i=0;i<nv;i++){ double u,v,z; projectPoint(V[i],view,R,u,v,z); U[i]=(float)u; W[i]=(float)v; Z[i]=(float)z; }
        float* zbuf=rm.depth.data()+(size_t)view*PIX; float* nbuf=rm.norm.data()+(size_t)view*PIX*3;
        for(int ti=0;ti<nf;ti++){
            int ia=F[ti].a,ib=F[ti].b,ic=F[ti].c;
            float x0=U[ia],x1=U[ib],x2=U[ic], y0=W[ia],y1=W[ib],y2=W[ic], z0=Z[ia],z1=Z[ib],z2=Z[ic];
            if(!(z0>0&&z1>0&&z2>0)) continue;
            int xmin=max(0,(int)floor(min(x0,min(x1,x2))-0.5f)); int xmax=min(R-1,(int)ceil(max(x0,max(x1,x2))+0.5f));
            int ymin=max(0,(int)floor(min(y0,min(y1,y2))-0.5f)); int ymax=min(R-1,(int)ceil(max(y0,max(y1,y2))+0.5f));
            if(xmin>xmax||ymin>ymax) continue;
            float den=(y1-y2)*(x0-x2)+(x2-x1)*(y0-y2); if(fabs(den)<1e-20f) continue; float inv=1.0f/den;
            Vec3 n=fn[ti]; float nr=(float)((n.x+1)*127.5),ng=(float)((n.y+1)*127.5),nb=(float)((n.z+1)*127.5);
            for(int y=ymin;y<=ymax;y++){
                float py=y+0.5f; int row=y*R;
                for(int x=xmin;x<=xmax;x++){
                    float px=x+0.5f;
                    float w0=((y1-y2)*(px-x2)+(x2-x1)*(py-y2))*inv;
                    float w1=((y2-y0)*(px-x2)+(x0-x2)*(py-y2))*inv;
                    float w2=1.f-w0-w1;
                    if(w0>=-1e-6f&&w1>=-1e-6f&&w2>=-1e-6f){
                        float dep=1.f/(w0/z0+w1/z1+w2/z2);
                        int idx=row+x;
                        if(dep<zbuf[idx]){
                            zbuf[idx]=dep;
                            float*q=nbuf+idx*3; q[0]=nr;q[1]=ng;q[2]=nb;
                        }
                    }
                }
            }
        }
    }
}
static inline double rectSum(const vector<double>&I,int S,int x0,int y0,int x1,int y1){
    return I[(size_t)y1*S+x1]-I[(size_t)y0*S+x1]-I[(size_t)y1*S+x0]+I[(size_t)y0*S+x0];
}
static double ssimChannel(const float*X,int sx,const float*Y,int sy,const float*dX,const float*dY,int R,
                          vector<double>&IX,vector<double>&IY,vector<double>&IX2,vector<double>&IY2,vector<double>&IXY){
    int S=R+1; size_t SZ=(size_t)S*S;
    fill(IX.begin(),IX.begin()+SZ,0); fill(IY.begin(),IY.begin()+SZ,0);
    fill(IX2.begin(),IX2.begin()+SZ,0); fill(IY2.begin(),IY2.begin()+SZ,0); fill(IXY.begin(),IXY.begin()+SZ,0);
    for(int y=1;y<=R;y++){
        double ax=0,ay=0,ax2=0,ay2=0,axy=0; int row=(y-1)*R;
        for(int x=1;x<=R;x++){
            int p=row+x-1; double vx=X[(size_t)p*sx], vy=Y[(size_t)p*sy];
            ax+=vx;ay+=vy;ax2+=vx*vx;ay2+=vy*vy;axy+=vx*vy;
            size_t id=(size_t)y*S+x, up=(size_t)(y-1)*S+x;
            IX[id]=IX[up]+ax; IY[id]=IY[up]+ay; IX2[id]=IX2[up]+ax2; IY2[id]=IY2[up]+ay2; IXY[id]=IXY[up]+axy;
        }
    }
    const int rad=5, area=121; const double c1=6.5025,c2=58.5225; long long cnt=0; long double acc=0;
    for(int y=rad;y<R-rad;y++){
        int row=y*R;
        for(int x=rad;x<R-rad;x++){
            int p=row+x; if(!(dX[p]<254.f||dY[p]<254.f)) continue;
            int x0=x-rad,y0=y-rad,x1=x+rad+1,y1=y+rad+1;
            double sumx=rectSum(IX,S,x0,y0,x1,y1), sumy=rectSum(IY,S,x0,y0,x1,y1);
            double sumx2=rectSum(IX2,S,x0,y0,x1,y1), sumy2=rectSum(IY2,S,x0,y0,x1,y1), sumxy=rectSum(IXY,S,x0,y0,x1,y1);
            double ux=sumx/area, uy=sumy/area;
            double vx=max(0.0,sumx2/area-ux*ux), vy=max(0.0,sumy2/area-uy*uy), cov=sumxy/area-ux*uy;
            double den=(ux*ux+uy*uy+c1)*(vx+vy+c2);
            double val=den!=0?((2*ux*uy+c1)*(2*cov+c2)/den):1.0;
            acc+=val; cnt++;
        }
    }
    return cnt? (double)(acc/cnt):1.0;
}
static double proxySSIM(const RenderMaps&orig,const RenderMaps&cand){
    int R=orig.R, PIX=R*R, S=R+1; size_t SZ=(size_t)S*S; vector<double>IX(SZ),IY(SZ),IX2(SZ),IY2(SZ),IXY(SZ); double total=0;
    for(int view=0;view<6;view++){
        const float*od=orig.depth.data()+(size_t)view*PIX; const float*cd=cand.depth.data()+(size_t)view*PIX;
        double ns=0;
        for(int ch=0;ch<3;ch++){
            const float*on=orig.norm.data()+(size_t)view*PIX*3+ch;
            const float*cn=cand.norm.data()+(size_t)view*PIX*3+ch;
            ns+=ssimChannel(on,3,cn,3,od,cd,R,IX,IY,IX2,IY2,IXY);
        }
        ns/=3.0;
        double ds=ssimChannel(od,1,cd,1,od,cd,R,IX,IY,IX2,IY2,IXY);
        total += 0.5*(ns+ds);
    }
    return total/6.0;
}
static int proxyRes(){ if(N<8000) return 1024; if(N<30000) return 512; if(N<120000) return 256; return 128; }
static double proxyGuard(int R,double ratio){ if(R>=1000) return ratio<0.08?0.915:0.905; if(R>=512) return ratio<0.08?0.925:0.915; if(R>=256) return ratio<0.08?0.945:0.935; return ratio<0.08?0.965:0.955; }
static bool acceptCandidate(const vector<Vec3>&V,const vector<Tri>&F,const RenderMaps&origMaps,int R,double guardExtra=0.0){
    if((int)V.size()>=N || V.empty() || F.empty()) return false;
    if(!manifoldOK(V,F)) return false;
    if(!vertexHausdorffOK(V)) return false;
    if(elapsed()>19.0) return true;
    RenderMaps cand; renderMesh(V,F,cand,R);
    double sc=proxySSIM(origMaps,cand); double ratio=(double)V.size()/max(1,N);
    return sc>=proxyGuard(R,ratio)+guardExtra;
}

static double orientSign(const vector<Vec3>&V,const vector<Tri>&F,Vec3 center){
    long double s=0; int step=max(1,(int)F.size()/200000);
    for(int i=0;i<(int)F.size();i+=step){
        auto&t=F[i]; Vec3 cr=crossv(V[t.b]-V[t.a],V[t.c]-V[t.a]); Vec3 c=(V[t.a]+V[t.b]+V[t.c])/3.0; s+=dotv(cr,c-center);
    }
    return s>=0?1.0:-1.0;
}
static void orientFace(vector<Vec3>&V,Tri&f,Vec3 center,double sgn){
    Vec3 cr=crossv(V[f.b]-V[f.a],V[f.c]-V[f.a]); Vec3 c=(V[f.a]+V[f.b]+V[f.c])/3.0;
    if(dotv(cr,c-center)*sgn<0) swap(f.b,f.c);
}
static void jacobi3(double a[3][3], Vec3 out[3]){
    double v[3][3]={{1,0,0},{0,1,0},{0,0,1}};
    for(int it=0;it<40;it++){
        int p=0,q=1; double best=fabs(a[0][1]);
        if(fabs(a[0][2])>best) p=0,q=2,best=fabs(a[0][2]);
        if(fabs(a[1][2])>best) p=1,q=2,best=fabs(a[1][2]);
        if(best<1e-18) break;
        double app=a[p][p], aqq=a[q][q], apq=a[p][q];
        double tau=(aqq-app)/(2*apq); double t=(tau>=0?1:-1)/(fabs(tau)+sqrt(1+tau*tau)); double c=1/sqrt(1+t*t), s=t*c;
        for(int k=0;k<3;k++) if(k!=p&&k!=q){
            double akp=a[k][p], akq=a[k][q];
            a[k][p]=a[p][k]=c*akp-s*akq; a[k][q]=a[q][k]=s*akp+c*akq;
        }
        a[p][p]=c*c*app-2*s*c*apq+s*s*aqq; a[q][q]=s*s*app+2*s*c*apq+c*c*aqq; a[p][q]=a[q][p]=0;
        for(int k=0;k<3;k++){ double vkp=v[k][p], vkq=v[k][q]; v[k][p]=c*vkp-s*vkq; v[k][q]=s*vkp+c*vkq; }
    }
    int ord[3]={0,1,2}; sort(ord,ord+3,[&](int i,int j){return a[i][i]>a[j][j];});
    for(int j=0;j<3;j++){int col=ord[j]; out[j]=normalize(Vec3(v[0][col],v[1][col],v[2][col]));}
    if(dotv(crossv(out[0],out[1]),out[2])<0) out[2]=out[2]*-1;
}
static bool pcaAxes(Vec3 ax[3]){
    if(N<20) return false; Vec3 mean;
    for(auto&p:Orig) mean+=p;
    mean=mean/(double)N;
    double c[3][3]={{0}};
    int st=max(1,N/300000), cnt=0;
    for(int i=0;i<N;i+=st){
        Vec3 q=Orig[i]-mean; double x[3]={q.x,q.y,q.z};
        for(int a=0;a<3;a++) for(int b=0;b<3;b++) c[a][b]+=x[a]*x[b];
        cnt++;
    }
    for(int a=0;a<3;a++) for(int b=0;b<3;b++) c[a][b]/=max(1,cnt);
    jacobi3(c,ax);
    return norm2(ax[0])>0.5&&norm2(ax[1])>0.5&&norm2(ax[2])>0.5;
}
struct EllipsoidFit{bool ok=false; Vec3 c,ax[3]; double r[3]; double rms=1e9,maxe=1e9;};
static EllipsoidFit fitEllipsoid(bool usePCA){
    EllipsoidFit fit; Vec3 axes[3]={{1,0,0},{0,1,0},{0,0,1}};
    if(usePCA && !pcaAxes(axes)) return fit;
    for(int k=0;k<3;k++) fit.ax[k]=axes[k];
    double lo[3]={1e100,1e100,1e100}, hi[3]={-1e100,-1e100,-1e100};
    for(auto&p:Orig) for(int k=0;k<3;k++){ double t=dotv(p,axes[k]); lo[k]=min(lo[k],t); hi[k]=max(hi[k],t); }
    fit.c=Vec3();
    for(int k=0;k<3;k++){
        double mid=(lo[k]+hi[k])*0.5; fit.c+=axes[k]*mid; fit.r[k]=(hi[k]-lo[k])*0.5;
        if(fit.r[k]<1e-10) return fit;
    }
    int st=max(1,N/300000),cnt=0; long double ss=0; double mx=0;
    for(int i=0;i<N;i+=st){
        Vec3 q=Orig[i]-fit.c; double s=0;
        for(int k=0;k<3;k++){ double u=dotv(q,axes[k])/fit.r[k]; s+=u*u; }
        double e=fabs(sqrt(max(0.0,s))-1.0); ss+=e*e; mx=max(mx,e); cnt++;
    }
    fit.rms=sqrt((double)ss/max(1,cnt)); fit.maxe=mx;
    double rmsLim=N<5000?0.010:0.007, maxLim=N<5000?0.040:0.028;
    fit.ok=fit.rms<=rmsLim && fit.maxe<=maxLim;
    return fit;
}
static bool makeEllipsoid(const EllipsoidFit&fit,int lat,int lon,vector<Vec3>&V,vector<Tri>&F){
    if(!fit.ok||lat<4||lon<8) return false; int vc=2+(lat-1)*lon; if(vc>N) return false;
    V.clear();F.clear();V.reserve(vc);F.reserve(2*lat*lon);
    const double pi=acos(-1);
    auto mp=[&](double x,double y,double z){return fit.c+fit.ax[0]*(fit.r[0]*x)+fit.ax[1]*(fit.r[1]*y)+fit.ax[2]*(fit.r[2]*z);};
    V.push_back(mp(0,0,1));
    for(int i=1;i<lat;i++){
        double th=pi*i/lat, st=sin(th), ct=cos(th);
        for(int j=0;j<lon;j++){double ph=2*pi*j/lon; V.push_back(mp(st*cos(ph),st*sin(ph),ct));}
    }
    V.push_back(mp(0,0,-1));
    int bottom=(int)V.size()-1;
    auto id=[&](int r,int j){j=(j%lon+lon)%lon; return 1+(r-1)*lon+j;};
    double sg=orientSign(Orig,OrigF,fit.c);
    auto add2=[&](int a,int b,int c){Tri t{a,b,c}; orientFace(V,t,fit.c,sg); F.push_back(t);};
    for(int j=0;j<lon;j++) add2(0,id(1,j+1),id(1,j));
    for(int r=1;r<lat-1;r++) for(int j=0;j<lon;j++){
        int a=id(r,j),b=id(r,j+1),c=id(r+1,j),d=id(r+1,j+1);
        add2(a,b,c); add2(b,d,c);
    }
    for(int j=0;j<lon;j++) add2(bottom,id(lat-1,j),id(lat-1,j+1));
    return true;
}
struct TorusFit{bool ok=false; int axis=2; double ct,cu,cv,major,minor,rms,maxe;};
static inline void coordBM(const Vec3&p,int ax,double&t,double&u,double&v){
    if(ax==0){t=p.x;u=p.y;v=p.z;}
    else if(ax==1){t=p.y;u=p.x;v=p.z;}
    else {t=p.z;u=p.x;v=p.y;}
}
static inline Vec3 fromBM(int ax,double t,double u,double v){ if(ax==0) return Vec3(t,u,v); if(ax==1) return Vec3(u,t,v); return Vec3(u,v,t); }
static TorusFit fitTorusAxis(int ax){
    TorusFit fit; fit.axis=ax; if(N<600) return fit;
    double lt=1e100,ht=-1e100,lu=1e100,hu=-1e100,lv=1e100,hv=-1e100;
    for(auto&p:Orig){double t,u,v;coordBM(p,ax,t,u,v);lt=min(lt,t);ht=max(ht,t);lu=min(lu,u);hu=max(hu,u);lv=min(lv,v);hv=max(hv,v);}
    fit.ct=(lt+ht)/2; fit.cu=(lu+hu)/2; fit.cv=(lv+hv)/2;
    double minr=1e100,maxr=0;
    for(auto&p:Orig){double t,u,v;coordBM(p,ax,t,u,v); double r=hypot(u-fit.cu,v-fit.cv); minr=min(minr,r); maxr=max(maxr,r);}
    if(maxr<=minr||minr<1e-9) return fit;
    fit.major=(maxr+minr)/2;
    double minorR=(maxr-minr)/2, minorT=(ht-lt)/2; fit.minor=(minorR+minorT)/2;
    if(fit.major<1.25*fit.minor || fabs(minorR-minorT)>0.30*fit.minor) return fit;
    int st=max(1,N/300000),cnt=0; long double ss=0; double mx=0;
    for(int i=0;i<N;i+=st){
        double t,u,v;coordBM(Orig[i],ax,t,u,v);
        double rho=hypot(u-fit.cu,v-fit.cv);
        double tube=hypot(rho-fit.major,t-fit.ct);
        double e=fabs(tube-fit.minor)/fit.minor;
        ss+=e*e; mx=max(mx,e); cnt++;
    }
    fit.rms=sqrt((double)ss/max(1,cnt)); fit.maxe=mx;
    double rlim=N<5000?0.025:0.016, mlim=N<5000?0.10:0.070;
    fit.ok=fit.rms<=rlim&&fit.maxe<=mlim; return fit;
}
static bool makeTorus(const TorusFit&fit,int seg,int tube,vector<Vec3>&V,vector<Tri>&F){
    if(!fit.ok||seg<12||tube<6||seg*tube>N) return false;
    V.clear();F.clear();V.reserve(seg*tube);F.reserve(2*seg*tube);
    const double pi=acos(-1); Vec3 center=fromBM(fit.axis,fit.ct,fit.cu,fit.cv); double sg=orientSign(Orig,OrigF,center);
    for(int i=0;i<seg;i++){
        double th=2*pi*i/seg, ct=cos(th), st=sin(th);
        for(int j=0;j<tube;j++){
            double ph=2*pi*j/tube, cp=cos(ph), sp=sin(ph);
            double rad=fit.major+fit.minor*cp;
            V.push_back(fromBM(fit.axis,fit.ct+fit.minor*sp,fit.cu+rad*ct,fit.cv+rad*st));
        }
    }
    auto id=[&](int i,int j){i=(i%seg+seg)%seg;j=(j%tube+tube)%tube;return i*tube+j;};
    auto add=[&](int a,int b,int c){Tri t{a,b,c}; orientFace(V,t,center,sg); F.push_back(t);};
    for(int i=0;i<seg;i++) for(int j=0;j<tube;j++){int a=id(i,j),b=id(i+1,j),c=id(i+1,j+1),d=id(i,j+1); add(a,b,c); add(a,c,d);}
    return true;
}
struct RevolveFit{bool ok=false; int axis=2; double t0,t1,cu,cv,r0,r1,rms,maxe;};
static RevolveFit fitRevolveAxis(int ax){
    RevolveFit fit; fit.axis=ax; if(N<800) return fit;
    double lt=1e100,ht=-1e100,lu=1e100,hu=-1e100,lv=1e100,hv=-1e100;
    for(auto&p:Orig){double t,u,v;coordBM(p,ax,t,u,v);lt=min(lt,t);ht=max(ht,t);lu=min(lu,u);hu=max(hu,u);lv=min(lv,v);hv=max(hv,v);}
    if(ht<=lt) return fit;
    fit.t0=lt;fit.t1=ht;fit.cu=(lu+hu)/2;fit.cv=(lv+hv)/2;
    double maxr=0; for(auto&p:Orig){double t,u,v;coordBM(p,ax,t,u,v);maxr=max(maxr,hypot(u-fit.cu,v-fit.cv));}
    if(maxr<1e-9) return fit;
    int st=max(1,N/300000),cnt=0,axisPts=0;
    double eps=maxr*.06; long double S=0,St=0,Sr=0,Stt=0,Str=0;
    for(int i=0;i<N;i+=st){
        double t,u,v;coordBM(Orig[i],ax,t,u,v);
        double r=hypot(u-fit.cu,v-fit.cv);
        if(r<eps){axisPts++; continue;}
        S++;St+=t;Sr+=r;Stt+=t*t;Str+=t*r;cnt++;
    }
    if(cnt<200) return fit;
    double den=S*Stt-St*St; if(fabs(den)<1e-18) return fit;
    double a=(S*Str-St*Sr)/den, b=(Sr-a*St)/S;
    fit.r0=max(0.0,a*lt+b); fit.r1=max(0.0,a*ht+b);
    if(max(fit.r0,fit.r1)<0.25*maxr) return fit;
    if(axisPts>cnt/3+10) return fit;
    long double ss=0; double mx=0; int chk=0;
    for(int i=0;i<N;i+=st){
        double t,u,v;coordBM(Orig[i],ax,t,u,v);
        double r=hypot(u-fit.cu,v-fit.cv);
        if(r<eps) continue;
        double pred=max(0.0,a*t+b);
        double e=fabs(r-pred)/maxr;
        ss+=e*e; mx=max(mx,e); chk++;
    }
    if(chk<200) return fit;
    fit.rms=sqrt((double)ss/chk); fit.maxe=mx; fit.ok=fit.rms<0.010 && fit.maxe<0.045;
    return fit;
}
static bool makeRevolve(const RevolveFit&fit,int sides,vector<Vec3>&V,vector<Tri>&F){
    if(!fit.ok||sides<12) return false;
    double eps=max(fit.r0,fit.r1)*1e-5; bool cone0=fit.r0<=eps, cone1=fit.r1<=eps;
    if(cone0&&cone1) return false;
    V.clear();F.clear();
    const double pi=acos(-1); Vec3 center=fromBM(fit.axis,(fit.t0+fit.t1)/2,fit.cu,fit.cv); double sg=orientSign(Orig,OrigF,center);
    auto makep=[&](double t,double r,int j){double th=2*pi*j/sides; return fromBM(fit.axis,t,fit.cu+r*cos(th),fit.cv+r*sin(th));};
    auto add=[&](int a,int b,int c){Tri t{a,b,c}; orientFace(V,t,center,sg); F.push_back(t);};
    if(!cone0&&!cone1){
        if(2+2*sides>N) return false;
        V.push_back(fromBM(fit.axis,fit.t0,fit.cu,fit.cv));
        V.push_back(fromBM(fit.axis,fit.t1,fit.cu,fit.cv));
        int lo=V.size(); for(int j=0;j<sides;j++) V.push_back(makep(fit.t0,fit.r0,j));
        int hi=V.size(); for(int j=0;j<sides;j++) V.push_back(makep(fit.t1,fit.r1,j));
        for(int j=0;j<sides;j++){int k=(j+1)%sides; add(lo+j,lo+k,hi+j); add(lo+k,hi+k,hi+j); add(0,lo+j,lo+k); add(1,hi+k,hi+j);}
    } else {
        bool bottomApex=cone0;
        double at=bottomApex?fit.t0:fit.t1, bt=bottomApex?fit.t1:fit.t0, br=bottomApex?fit.r1:fit.r0;
        if(2+sides>N) return false;
        V.push_back(fromBM(fit.axis,at,fit.cu,fit.cv));
        V.push_back(fromBM(fit.axis,bt,fit.cu,fit.cv));
        int ring=V.size(); for(int j=0;j<sides;j++) V.push_back(makep(bt,br,j));
        for(int j=0;j<sides;j++){int k=(j+1)%sides; add(0,ring+j,ring+k); add(1,ring+k,ring+j);}
    }
    return true;
}
static bool makeBox(vector<Vec3>&V,vector<Tri>&F){
    if(N<8) return false;
    V={{g_mn.x,g_mn.y,g_mn.z},{g_mx.x,g_mn.y,g_mn.z},{g_mx.x,g_mx.y,g_mn.z},{g_mn.x,g_mx.y,g_mn.z},
       {g_mn.x,g_mn.y,g_mx.z},{g_mx.x,g_mn.y,g_mx.z},{g_mx.x,g_mx.y,g_mx.z},{g_mn.x,g_mx.y,g_mx.z}};
    int T[12][3]={{0,1,2},{0,2,3},{4,6,5},{4,7,6},{0,4,5},{0,5,1},{1,5,6},{1,6,2},{2,6,7},{2,7,3},{3,7,4},{3,4,0}};
    F.clear(); Vec3 c=(g_mn+g_mx)/2; double sg=orientSign(Orig,OrigF,c);
    for(auto&t:T){Tri f{t[0],t[1],t[2]}; orientFace(V,f,c,sg); F.push_back(f);}
    return true;
}

struct Quadric{
    double q[10];
    Quadric(){memset(q,0,sizeof q);}
    void add(const Quadric&o){for(int i=0;i<10;i++)q[i]+=o.q[i];}
    void addPlane(double a,double b,double c,double d,double w){
        q[0]+=w*a*a;q[1]+=w*a*b;q[2]+=w*a*c;q[3]+=w*a*d;q[4]+=w*b*b;q[5]+=w*b*c;q[6]+=w*b*d;q[7]+=w*c*c;q[8]+=w*c*d;q[9]+=w*d*d;
    }
    double eval(const Vec3&p)const{
        double x=p.x,y=p.y,z=p.z;
        return q[0]*x*x+2*q[1]*x*y+2*q[2]*x*z+2*q[3]*x+q[4]*y*y+2*q[5]*y*z+2*q[6]*y+q[7]*z*z+2*q[8]*z+q[9];
    }
};
struct Snapshot{ vector<Vec3> V; vector<Tri> F; double ratio=1, ssim=-1; };
class Simplifier{
public:
    int n,m; vector<Vec3>P; vector<Face>F; vector<vector<int>> inc; vector<Quadric>Q; vector<array<float,3>> bmin,bmax;
    vector<unsigned char> alive; vector<int> ver, markA, markB, tmp; int stampA=1,stampB=1; int activeV=0,activeF=0;
    double haus2=0, pixelTol=24, pixelTol2=576; vector<Snapshot> snaps;
    struct Cand{double score; int u,v,vu,vv; bool operator<(const Cand&o)const{return score>o.score;}};
    priority_queue<Cand> pq;
    Simplifier(){n=N;m=M;P=Orig;F.resize(M);for(int i=0;i<M;i++){F[i].v[0]=OrigF[i].a;F[i].v[1]=OrigF[i].b;F[i].v[2]=OrigF[i].c;F[i].active=1;}}
    bool contains(int fid,int v)const{const Face&f=F[fid];return f.v[0]==v||f.v[1]==v||f.v[2]==v;}
    Vec3 faceNormal(const Face&f)const{Vec3 cr=crossv(P[f.v[1]]-P[f.v[0]],P[f.v[2]]-P[f.v[0]]); double l=normv(cr);return l>1e-300?cr/l:Vec3();}
    void cleanup(int u){ if(!alive[u]) return; auto&v=inc[u]; int w=0; for(int fid:v) if(fid>=0&&fid<m&&F[fid].active&&contains(fid,u)) v[w++]=fid; v.resize(w); }
    void maybe(int u){ if(alive[u]&&inc[u].size()>96){int good=0;for(int fid:inc[u])if(F[fid].active&&contains(fid,u))good++; if((int)inc[u].size()>good*3+64) cleanup(u);} }
    void init(){
        haus2=g_haus*g_haus; if(n<8000) pixelTol=32; else if(n<60000) pixelTol=26; else if(n<200000) pixelTol=22; else pixelTol=18; pixelTol2=pixelTol*pixelTol;
        vector<int>deg(n,0); for(auto&t:OrigF){deg[t.a]++;deg[t.b]++;deg[t.c]++;}
        inc.assign(n,{}); for(int i=0;i<n;i++) inc[i].reserve(deg[i]+8); activeV=n; activeF=m;
        for(int i=0;i<m;i++){F[i].n=faceNormal(F[i]); inc[F[i].v[0]].push_back(i); inc[F[i].v[1]].push_back(i); inc[F[i].v[2]].push_back(i);}
        Q.assign(n,Quadric());
        for(int i=0;i<m;i++){
            auto&f=F[i]; Vec3 cr=crossv(P[f.v[1]]-P[f.v[0]],P[f.v[2]]-P[f.v[0]]);
            double area2=normv(cr); if(area2<=1e-300) continue;
            Vec3 nn=cr/area2; double d=-dotv(nn,P[f.v[0]]); double w=max(1e-8,0.5*area2)+2e-7;
            for(int k=0;k<3;k++) Q[f.v[k]].addPlane(nn.x,nn.y,nn.z,d,w);
        }
        bmin.resize(n); bmax.resize(n);
        for(int i=0;i<n;i++){bmin[i]={(float)P[i].x,(float)P[i].y,(float)P[i].z}; bmax[i]=bmin[i];}
        alive.assign(n,1); ver.assign(n,0); markA.assign(n,0); markB.assign(n,0); tmp.assign(n,-1);
        vector<uint64_t> edges; edges.reserve((size_t)m*3);
        for(auto&t:OrigF){edges.push_back(ekey(t.a,t.b));edges.push_back(ekey(t.b,t.c));edges.push_back(ekey(t.c,t.a));}
        sort(edges.begin(),edges.end()); edges.erase(unique(edges.begin(),edges.end()),edges.end());
        vector<Cand> heap; heap.reserve(edges.size());
        for(uint64_t k:edges){int a=k>>32,b=(int)(k&0xffffffffu); auto c=evalCand(a,b); if(c.first) heap.push_back({c.second,a,b,0,0});}
        pq=priority_queue<Cand>(less<Cand>(),move(heap));
    }
    bool solveOpt(const Quadric&q,Vec3&out){
        double a00=q.q[0],a01=q.q[1],a02=q.q[2],a11=q.q[4],a12=q.q[5],a22=q.q[7];
        double b0=q.q[3],b1=q.q[6],b2=q.q[8];
        double det=a00*(a11*a22-a12*a12)-a01*(a01*a22-a12*a02)+a02*(a01*a12-a11*a02);
        if(fabs(det)<1e-14) return false;
        double inv00=(a11*a22-a12*a12)/det,inv01=(a02*a12-a01*a22)/det,inv02=(a01*a12-a02*a11)/det,inv11=(a00*a22-a02*a02)/det,inv12=(a01*a02-a00*a12)/det,inv22=(a00*a11-a01*a01)/det;
        out.x=-(inv00*b0+inv01*b1+inv02*b2); out.y=-(inv01*b0+inv11*b1+inv12*b2); out.z=-(inv02*b0+inv12*b1+inv22*b2);
        return isfinite(out.x)&&isfinite(out.y)&&isfinite(out.z)&&out.x>=-2&&out.x<=2&&out.y>=-2&&out.y<=2&&out.z>=-2&&out.z<=2;
    }
    void mergedBox(int u,int v,array<float,3>&mn,array<float,3>&mx){for(int d=0;d<3;d++){mn[d]=min(bmin[u][d],bmin[v][d]);mx[d]=max(bmax[u][d],bmax[v][d]);}}
    bool coversBox(const Vec3&p,const array<float,3>&mn,const array<float,3>&mx){
        double md=0;
        for(int mask=0;mask<8;mask++){
            double x=(mask&1)?mx[0]:mn[0],y=(mask&2)?mx[1]:mn[1],z=(mask&4)?mx[2]:mn[2];
            double d=(p.x-x)*(p.x-x)+(p.y-y)*(p.y-y)+(p.z-z)*(p.z-z); md=max(md,d);
        }
        return md<=haus2;
    }
    bool visualOK(const Vec3&p,const array<float,3>&mn,const array<float,3>&mx){
        double pu[6],pv[6],pd;
        for(int view=0;view<6;view++) projectPoint(p,view,1024,pu[view],pv[view],pd);
        for(int mask=0;mask<8;mask++){
            Vec3 c((mask&1)?mx[0]:mn[0],(mask&2)?mx[1]:mn[1],(mask&4)?mx[2]:mn[2]);
            for(int view=0;view<6;view++){
                double u,v,d; projectPoint(c,view,1024,u,v,d);
                double e=(u-pu[view])*(u-pu[view])+(v-pv[view])*(v-pv[view]);
                if(e>pixelTol2) return false;
            }
        }
        return true;
    }
    pair<bool,double> evalCand(int u,int v){
        if(u==v||!alive[u]||!alive[v]) return {false,0};
        array<float,3>mn,mx; mergedBox(u,v,mn,mx); Quadric q=Q[u]; q.add(Q[v]);
        Vec3 cand[7]; int cc=0; Vec3 opt; if(solveOpt(q,opt)) cand[cc++]=opt;
        cand[cc++]=P[u]; cand[cc++]=P[v]; cand[cc++]=(P[u]+P[v])*0.5; cand[cc++]=Vec3((mn[0]+mx[0])*0.5,(mn[1]+mx[1])*0.5,(mn[2]+mx[2])*0.5);
        cand[cc++]=P[u]*0.67+P[v]*0.33; cand[cc++]=P[u]*0.33+P[v]*0.67;
        double len2=norm2(P[u]-P[v]); double best=1e300; bool ok=false;
        for(int i=0;i<cc;i++){
            if(!coversBox(cand[i],mn,mx))continue; if(!visualOK(cand[i],mn,mx))continue;
            double e=q.eval(cand[i]); if(e<0&&e>-1e-10)e=0; double sc=e+2e-10*len2+5e-15*norm2(cand[i]);
            if(sc<best){best=sc;ok=true;}
        }
        return {ok,best};
    }
    bool bestPos(int u,int v,Vec3&pos){
        array<float,3>mn,mx; mergedBox(u,v,mn,mx); Quadric q=Q[u]; q.add(Q[v]);
        Vec3 cand[7]; int cc=0; Vec3 opt; if(solveOpt(q,opt)) cand[cc++]=opt;
        cand[cc++]=P[u]; cand[cc++]=P[v]; cand[cc++]=(P[u]+P[v])*0.5; cand[cc++]=Vec3((mn[0]+mx[0])*0.5,(mn[1]+mx[1])*0.5,(mn[2]+mx[2])*0.5);
        cand[cc++]=P[u]*0.67+P[v]*0.33; cand[cc++]=P[u]*0.33+P[v]*0.67;
        double len2=norm2(P[u]-P[v]); double best=1e300; bool ok=false;
        for(int i=0;i<cc;i++){
            if(!coversBox(cand[i],mn,mx))continue; if(!visualOK(cand[i],mn,mx))continue;
            double e=q.eval(cand[i]); if(e<0&&e>-1e-10)e=0; double sc=e+2e-10*len2+5e-15*norm2(cand[i]);
            if(sc<best){best=sc;pos=cand[i];ok=true;}
        }
        return ok;
    }
    bool linkOK(int u,int v){
        if(u==v||!alive[u]||!alive[v]) return false; maybe(u); maybe(v);
        if(++stampA>1000000000){fill(markA.begin(),markA.end(),0);stampA=1;}
        if(++stampB>1000000000){fill(markB.begin(),markB.end(),0);stampB=1;}
        int edgeFaces=0;
        for(int fid:inc[u]){
            auto&f=F[fid]; if(!f.active||!contains(fid,u)) continue; bool hv=false;
            for(int k=0;k<3;k++){int w=f.v[k]; if(w==v) hv=true; if(w!=u) markA[w]=stampA;}
            if(hv) edgeFaces++;
        }
        if(edgeFaces!=2) return false;
        int common=0;
        for(int fid:inc[v]){
            auto&f=F[fid]; if(!f.active||!contains(fid,v)) continue;
            for(int k=0;k<3;k++){
                int w=f.v[k]; if(w==u||w==v) continue;
                if(markA[w]==stampA && markB[w]!=stampB){markB[w]=stampB; if(++common>2) return false;}
            }
        }
        return common==2;
    }
    bool flipOK(int keep,int rem,const Vec3&pos){
        static vector<int> touched; touched.clear();
        for(int fid:inc[keep]) if(F[fid].active&&contains(fid,keep)) touched.push_back(fid);
        for(int fid:inc[rem]) if(F[fid].active&&contains(fid,rem)) touched.push_back(fid);
        sort(touched.begin(),touched.end()); touched.erase(unique(touched.begin(),touched.end()),touched.end());
        double minDot = (n<30000? -0.10 : -0.02);
        for(int fid:touched){
            auto&f=F[fid]; if(!f.active) continue; bool hk=false,hr=false;
            for(int k=0;k<3;k++){if(f.v[k]==keep)hk=true;if(f.v[k]==rem)hr=true;}
            if(hk&&hr) continue;
            Vec3 a=(f.v[0]==keep||f.v[0]==rem)?pos:P[f.v[0]], b=(f.v[1]==keep||f.v[1]==rem)?pos:P[f.v[1]], c=(f.v[2]==keep||f.v[2]==rem)?pos:P[f.v[2]];
            Vec3 cr=crossv(b-a,c-a); double l=normv(cr);
            if(!(l>1e-14*g_diag*g_diag)) return false;
            Vec3 nn=cr/l; if(dotv(nn,f.n)<minDot) return false;
        }
        return true;
    }
    void pushEdge(int u,int v){ if(u==v||!alive[u]||!alive[v])return; if(u>v)swap(u,v); auto ec=evalCand(u,v); if(ec.first) pq.push({ec.second,u,v,ver[u],ver[v]}); }
    void neighbors(int u,vector<int>&out){
        out.clear(); if(!alive[u]) return; maybe(u);
        if(++stampA>1000000000){fill(markA.begin(),markA.end(),0);stampA=1;}
        for(int fid:inc[u]){
            auto&f=F[fid]; if(!f.active||!contains(fid,u)) continue;
            for(int k=0;k<3;k++){int w=f.v[k]; if(w!=u&&alive[w]&&markA[w]!=stampA){markA[w]=stampA;out.push_back(w);}}
        }
    }
    void collapse(int keep,int rem,const Vec3&pos){
        cleanup(keep); cleanup(rem); vector<int> remFaces=inc[rem];
        for(int fid:remFaces){
            auto&f=F[fid]; if(!f.active||!contains(fid,rem)) continue;
            bool hasK=false; for(int k=0;k<3;k++) if(f.v[k]==keep) hasK=true;
            if(hasK){f.active=0;activeF--;}
            else {
                for(int k=0;k<3;k++) if(f.v[k]==rem) f.v[k]=keep;
                if(f.v[0]==f.v[1]||f.v[1]==f.v[2]||f.v[2]==f.v[0]){f.active=0;activeF--;}
                else {P[keep]=pos; f.n=faceNormal(f); inc[keep].push_back(fid);}
            }
        }
        P[keep]=pos;
        for(int fid:inc[keep]) if(F[fid].active&&contains(fid,keep)) F[fid].n=faceNormal(F[fid]);
        Q[keep].add(Q[rem]);
        for(int d=0;d<3;d++){bmin[keep][d]=min(bmin[keep][d],bmin[rem][d]);bmax[keep][d]=max(bmax[keep][d],bmax[rem][d]);}
        alive[rem]=0; ver[keep]++; ver[rem]++; activeV--; inc[rem].clear(); cleanup(keep);
        static vector<int> nb; neighbors(keep,nb); for(int w:nb) pushEdge(keep,w);
    }
    Snapshot snapshot(){
        Snapshot s; s.ratio=(double)activeV/max(1,n); tmp.assign(n,-1); s.V.reserve(activeV);
        for(int i=0;i<n;i++) if(alive[i]){tmp[i]=(int)s.V.size(); s.V.push_back(P[i]);}
        s.F.reserve(activeF);
        for(int i=0;i<m;i++) if(F[i].active){
            int a=F[i].v[0],b=F[i].v[1],c=F[i].v[2];
            if(a!=b&&b!=c&&c!=a&&tmp[a]>=0&&tmp[b]>=0&&tmp[c]>=0) s.F.push_back({tmp[a],tmp[b],tmp[c]});
        }
        return s;
    }
    vector<Snapshot> run(){
        init();
        vector<double> ratios={0.50,0.35,0.25,0.20,0.18,0.16,0.14,0.13,0.12,0.11,0.10,0.09,0.08,0.07,0.06};
        vector<int> cps; int last=-1;
        for(double r:ratios){int c=max(4,(int)ceil(n*r)); if(c<n&&c!=last){cps.push_back(c);last=c;}}
        sort(cps.begin(),cps.end(),greater<int>());
        int cp=0; int minTarget=cps.empty()?max(4,n/10):cps.back(); long long pops=0;
        while(!pq.empty()&&activeV>minTarget&&elapsed()<17.4){
            Cand c=pq.top();pq.pop();
            if((++pops&8191)==0 && elapsed()>17.4) break;
            int u=c.u,v=c.v; if(u==v||!alive[u]||!alive[v]) continue;
            if(u>v)swap(u,v);
            if(ver[u]!=c.vu||ver[v]!=c.vv) continue;
            Vec3 pos; if(!bestPos(u,v,pos)) continue;
            if(!linkOK(u,v)) continue;
            int keep=u,rem=v; if(inc[v].size()>inc[u].size()) keep=v,rem=u;
            if(!flipOK(keep,rem,pos)) continue;
            collapse(keep,rem,pos);
            while(cp<(int)cps.size()&&activeV<=cps[cp]){
                snaps.push_back(snapshot()); cp++;
                if(elapsed()>17.4) break;
            }
        }
        if(snaps.empty()||snaps.back().V.size()!=(size_t)activeV) snaps.push_back(snapshot());
        return snaps;
    }
};

static bool faceInsideQuad(const Tri&f,int a,int b,int c,int d){
    int q[4]={a,b,c,d};
    for(int x:{f.a,f.b,f.c}){
        bool ok=false; for(int y:q) if(x==y){ok=true;break;}
        if(!ok) return false;
    }
    return f.a!=f.b&&f.a!=f.c&&f.b!=f.c;
}
static vector<Vec3> originalVertexNormals(){
    vector<Vec3> n(N);
    for(const auto&t:OrigF){ Vec3 cr=crossv(Orig[t.b]-Orig[t.a],Orig[t.c]-Orig[t.a]); n[t.a]+=cr; n[t.b]+=cr; n[t.c]+=cr; }
    for(auto&x:n){ double l=normv(x); if(l>1e-300) x=x/l; }
    return n;
}
static void orientHint(const vector<Vec3>&V,Tri&f,const Vec3&hint){ Vec3 cr=crossv(V[f.b]-V[f.a],V[f.c]-V[f.a]); if(dotv(cr,hint)<0) swap(f.b,f.c); }
static bool detectUVGrid(int&U,int&S){
    if(N<300||M!=2*N) return false;
    vector<int> cand;
    auto add=[&](int s){ if(s>=6 && s<=N/6 && N%s==0 && find(cand.begin(),cand.end(),s)==cand.end()) cand.push_back(s); };
    unordered_map<int,int> freq; int step=max(1,M/120000);
    for(int i=0;i<M;i+=step){
        const Tri&f=OrigF[i]; int a[3]={f.a,f.b,f.c};
        for(int k=0;k<3;k++){ int d=abs(a[k]-a[(k+1)%3]); if(d==0) continue; d=min(d,N-d); if(d>=6 && d<=N/3) freq[d]++; }
    }
    vector<pair<int,int>> order; for(auto &kv:freq) order.push_back({kv.second,kv.first}); sort(order.rbegin(),order.rend());
    for(int i=0;i<(int)order.size()&&i<16;i++){ int d=order[i].second; for(int e=-2;e<=2;e++) add(d+e); if(d) add(N/d); }
    for(int s=6;s<=512 && (int)cand.size()<128;s++) if(N%s==0) add(s);
    for(int s:cand){
        int u=N/s; if(u<6) continue; int tot=0,ok=0,st=max(1,N/1200);
        for(int cell=0;cell<N && tot<1200;cell+=st){
            int i=cell/s,j=cell%s; int a=i*s+j,b=i*s+(j+1)%s,c=((i+1)%u)*s+j,d=((i+1)%u)*s+(j+1)%s;
            int fid=2*cell; if(fid+1>=M) break;
            bool m0=faceInsideQuad(OrigF[fid],a,b,c,d); bool m1=faceInsideQuad(OrigF[fid+1],a,b,c,d);
            if(m0&&m1) ok++; tot++;
        }
        if(tot>200 && ok*1000>=tot*995){U=u;S=s;return true;}
    }
    return false;
}
static bool makeUVGrid(int U,int S,int U2,int S2,const vector<Vec3>&VN,vector<Vec3>&V,vector<Tri>&F){
    if(U2<4||S2<4||U2*S2>=N) return false; V.clear();F.clear(); V.reserve(U2*S2); vector<int> src(U2*S2);
    for(int i=0;i<U2;i++){ int oi=(long long)i*U/U2; for(int j=0;j<S2;j++){ int oj=(long long)j*S/S2; int id=oi*S+oj; src[i*S2+j]=id; V.push_back(Orig[id]); }}
    auto id=[&](int i,int j){return ((i%U2+U2)%U2)*S2+((j%S2+S2)%S2);};
    F.reserve(2*U2*S2);
    for(int i=0;i<U2;i++) for(int j=0;j<S2;j++){
        int a=id(i,j),b=id(i+1,j),c=id(i+1,j+1),d=id(i,j+1);
        Vec3 h=VN[src[a]]+VN[src[b]]+VN[src[c]]+VN[src[d]];
        Tri f1{a,b,d}; orientHint(V,f1,h); F.push_back(f1);
        Tri f2{b,c,d}; orientHint(V,f2,h); F.push_back(f2);
    }
    return true;
}
static bool detectSphereGrid(int&R,int&S){
    if(N<300||M!=2*(N-2)) return false;
    vector<int> cand; for(int s=8;s<=2048;s++) if((N-2)%s==0) cand.push_back(s);
    for(int s:cand){
        int r=(N-2)/s; if(r<3) continue;
        int tot=0,ok=0; int st=max(1,s/64);
        for(int j=0;j<s && tot<256;j+=st){
            int a=1+j,b=1+(j+1)%s; const Tri&f=OrigF[j];
            bool m=((f.a==0||f.b==0||f.c==0)&&(f.a==a||f.b==a||f.c==a)&&(f.a==b||f.b==b||f.c==b));
            ok+=m; tot++;
        }
        int bandTot=(r-1)*s, bst=max(1,bandTot/512);
        for(int q=0;q<bandTot&&tot<900;q+=bst){
            int rr=q/s,j=q%s; int a=1+rr*s+j,b=1+rr*s+(j+1)%s,c=1+(rr+1)*s+j,d=1+(rr+1)*s+(j+1)%s;
            int fid=s+2*q; if(fid+1>=M) break;
            ok += faceInsideQuad(OrigF[fid],a,b,c,d) && faceInsideQuad(OrigF[fid+1],a,b,c,d); tot++;
        }
        if(tot>200 && ok*1000>=tot*995){R=r;S=s;return true;}
    }
    return false;
}
static bool makeSphereGrid(int R,int S,int R2,int S2,const vector<Vec3>&VN,vector<Vec3>&V,vector<Tri>&F){
    if(R2<3||S2<8||2+R2*S2>=N) return false;
    V.clear();F.clear(); V.reserve(2+R2*S2); vector<int> src; src.reserve(2+R2*S2);
    V.push_back(Orig[0]); src.push_back(0);
    for(int i=0;i<R2;i++){ int oi=1+(long long)i*(R-1)/max(1,R2-1); for(int j=0;j<S2;j++){ int oj=(long long)j*S/S2; int id=1+(oi-1)*S+oj; V.push_back(Orig[id]); src.push_back(id); }}
    int bottom=V.size(); V.push_back(Orig[N-1]); src.push_back(N-1);
    auto rid=[&](int r,int j){return 1+r*S2+((j%S2+S2)%S2);};
    F.reserve(2*R2*S2);
    auto add=[&](int a,int b,int c){Tri f{a,b,c}; Vec3 h=VN[src[a]]+VN[src[b]]+VN[src[c]]; orientHint(V,f,h); F.push_back(f);};
    for(int j=0;j<S2;j++) add(0,rid(0,j+1),rid(0,j));
    for(int r=0;r+1<R2;r++) for(int j=0;j<S2;j++){ int a=rid(r,j),b=rid(r,j+1),c=rid(r+1,j),d=rid(r+1,j+1); add(a,c,b); add(b,c,d); }
    for(int j=0;j<S2;j++) add(bottom,rid(R2-1,j),rid(R2-1,j+1));
    return true;
}

struct AxisCellSample{ int hi=-1, lo=-1; double hv=-1e100, lv=1e100; };
static inline double axisCoord(const Vec3&p,int ax){ return ax==0?p.x:(ax==1?p.y:p.z); }
static inline double axisU(const Vec3&p,int ax){ return ax==0?p.y:(ax==1?p.x:p.x); }
static inline double axisV(const Vec3&p,int ax){ return ax==0?p.z:(ax==1?p.z:p.y); }
static Vec3 nudgeAxis(Vec3 p,int ax,double d){ if(ax==0)p.x+=d; else if(ax==1)p.y+=d; else p.z+=d; return p; }
static bool addTriIfOK(const vector<Vec3>&V, vector<Tri>&F, int a,int b,int c){
    if(a==b||b==c||c==a) return false;
    if(norm2(crossv(V[b]-V[a],V[c]-V[a])) <= 1e-28*g_diag*g_diag*g_diag*g_diag) return false;
    F.push_back({a,b,c}); return true;
}
static bool buildOneAxisShell(int ax,int G,vector<Vec3>&V,vector<Tri>&F){
    vector<AxisCellSample> cell(G*G);
    double du=(axisU(g_mx,ax)-axisU(g_mn,ax));
    double dv=(axisV(g_mx,ax)-axisV(g_mn,ax));
    double u0=axisU(g_mn,ax), v0=axisV(g_mn,ax);
    if(!(du>1e-12)||!(dv>1e-12)) return false;
    for(int i=0;i<N;i++){
        const Vec3&p=Orig[i];
        int x=(int)((axisU(p,ax)-u0)/du*G); int y=(int)((axisV(p,ax)-v0)/dv*G);
        if(x<0)x=0; if(x>=G)x=G-1; if(y<0)y=0; if(y>=G)y=G-1;
        AxisCellSample &c=cell[y*G+x]; double val=axisCoord(p,ax);
        if(val>c.hv){c.hv=val;c.hi=i;} if(val<c.lv){c.lv=val;c.lo=i;}
    }
    vector<int> top(G*G,-1), bot(G*G,-1);
    const double sep=max(1e-8, min(0.001*g_diag,0.012*g_haus));
    int active=0;
    for(int id=0;id<G*G;id++) if(cell[id].hi>=0){
        Vec3 a=Orig[cell[id].hi]; Vec3 b=Orig[cell[id].lo>=0?cell[id].lo:cell[id].hi];
        if(fabs(axisCoord(a,ax)-axisCoord(b,ax))<sep*2){ a=nudgeAxis(a,ax,sep); b=nudgeAxis(b,ax,-sep); }
        top[id]=(int)V.size(); V.push_back(a); bot[id]=(int)V.size(); V.push_back(b); active++;
    }
    if(active<16) return false;
    vector<Tri> localTop;
    auto id=[&](int x,int y){return y*G+x;};
    auto addSheetTri=[&](int ca,int cb,int cc){
        int a=top[ca],b=top[cb],c=top[cc], aa=bot[ca],bb=bot[cb],ccv=bot[cc];
        if(a<0||b<0||c<0) return;
        if(norm2(crossv(V[b]-V[a],V[c]-V[a])) <= 1e-28*g_diag*g_diag*g_diag*g_diag) return;
        if(norm2(crossv(V[bb]-V[aa],V[ccv]-V[aa])) <= 1e-28*g_diag*g_diag*g_diag*g_diag) return;
        localTop.push_back({a,b,c});
        F.push_back({a,b,c});
        F.push_back({aa,ccv,bb});
    };
    for(int y=0;y+1<G;y++) for(int x=0;x+1<G;x++){
        int c00=id(x,y), c10=id(x+1,y), c01=id(x,y+1), c11=id(x+1,y+1);
        if(top[c00]>=0&&top[c10]>=0&&top[c01]>=0&&top[c11]>=0){
            addSheetTri(c00,c10,c01);
            addSheetTri(c10,c11,c01);
        }
    }
    if(localTop.empty()) return false;
    map<pair<int,int>, pair<int,int>> edgeDir;
    map<pair<int,int>, int> edgeCnt;
    for(const Tri&t:localTop){
        int a[3]={t.a,t.b,t.c};
        for(int k=0;k<3;k++){ int u=a[k],v=a[(k+1)%3]; pair<int,int> key=minmax(u,v); edgeCnt[key]++; edgeDir[key]={u,v}; }
    }
    auto bottomOf=[&](int topIndex)->int{return topIndex+1;};
    for(auto &kv: edgeCnt){
        if(kv.second!=1) continue;
        int a=edgeDir[kv.first].first, b=edgeDir[kv.first].second;
        int aa=bottomOf(a), bb=bottomOf(b);
        addTriIfOK(V,F,a,b,bb);
        addTriIfOK(V,F,a,bb,aa);
    }
    return true;
}
static bool tryAxisDepthShells(RenderMaps&origMaps,int R,vector<Vec3>&bestV,vector<Tri>&bestF){
    if(N<70000||elapsed()>11.5) return false;
    vector<int> grids;
    if(N>800000) grids={96,88,80,72};
    else if(N>300000) grids={80,72,64,56};
    else if(N>150000) grids={64,56,48};
    else grids={56,48,40};
    bool got=false;
    for(int G:grids){
        if(elapsed()>16.6) break;
        vector<Vec3> V; vector<Tri> F;
        V.reserve((size_t)6*G*G); F.reserve((size_t)12*G*G);
        bool ok=true;
        for(int ax=0; ax<3; ax++) if(!buildOneAxisShell(ax,G,V,F)){ ok=false; break; }
        if(!ok||V.empty()||F.empty()||V.size()>=Orig.size()) continue;
        if(!manifoldOK(V,F)) continue;
        if(!vertexHausdorffOK(V)) continue;
        RenderMaps cand; renderMesh(V,F,cand,R);
        double sc=proxySSIM(origMaps,cand);
        double ratio=(double)V.size()/max(1,N);
        double need=proxyGuard(R,ratio)+(R<=128?0.006:0.0);
        if(sc>=need){
            if(!got || V.size()<bestV.size()){ bestV.swap(V); bestF.swap(F); got=true; }
            break;
        }
    }
    return got;
}
static bool tryIndexedGrids(RenderMaps&origMaps,int R,vector<Vec3>&bestV,vector<Tri>&bestF){
    bool got=false; vector<Vec3>VN=originalVertexNormals(); vector<Vec3>V; vector<Tri>F;
    auto consider=[&](vector<Vec3>&V,vector<Tri>&F,double extra=0.0){
        if(elapsed()>17.0) return;
        if(V.empty()||F.empty()||(got&&V.size()>=bestV.size())) return;
        if(acceptCandidate(V,F,origMaps,R,extra)){bestV=V;bestF=F;got=true;}
    };
    int U=0,S=0;
    if(detectUVGrid(U,S)){
        int targets[]={512,768,1024,1536,2048,3072,4096,6144,8192,12288,16384,24576};
        double ar=sqrt((double)U/max(1,S));
        for(int t:targets){
            if(elapsed()>17.0) break;
            int u2=max(6,min(U,(int)(sqrt((double)t)*ar+0.5)));
            int s2=max(6,min(S,t/max(1,u2)));
            if(makeUVGrid(U,S,u2,s2,VN,V,F)) consider(V,F,0.0);
        }
    }
    int RR=0,SS=0;
    if(detectSphereGrid(RR,SS)){
        int targets[]={512,768,1024,1536,2048,3072,4096,6144};
        for(int t:targets){
            if(elapsed()>17.0) break;
            double ar=sqrt((double)RR/max(1,SS));
            int r2=max(4,min(RR,(int)(sqrt((double)t)*ar+0.5)));
            int s2=max(12,min(SS,t/max(1,r2)));
            if(makeSphereGrid(RR,SS,r2,s2,VN,V,F)) consider(V,F,0.0);
        }
    }
    return got;
}
static bool tryAnalytic(RenderMaps&origMaps,int R,vector<Vec3>&bestV,vector<Tri>&bestF){
    bool got=false;
    auto consider=[&](vector<Vec3>&V,vector<Tri>&F,double extra=0.0){
        if(elapsed()>16.0) return;
        if(V.empty()||F.empty()||(got&&V.size()>=bestV.size())) return;
        if(acceptCandidate(V,F,origMaps,R,extra)){ bestV=V; bestF=F; got=true; }
    };
    vector<Vec3> V; vector<Tri> F;
    if(makeBox(V,F)) consider(V,F,0.02);
    for(int mode=0;mode<2;mode++){
        EllipsoidFit fit=fitEllipsoid(mode==1);
        if(fit.ok){
            int trials[5][2]={{16,32},{18,36},{22,44},{26,52},{30,60}};
            for(auto&t:trials){ if(elapsed()>16.0) break; if(makeEllipsoid(fit,t[0],t[1],V,F)) consider(V,F,0.0); }
        }
    }
    for(int ax=0;ax<3;ax++){
        TorusFit fit=fitTorusAxis(ax);
        if(fit.ok){
            int trials[5][2]={{48,12},{64,16},{80,20},{96,24},{112,28}};
            for(auto&t:trials){ if(elapsed()>16.0) break; if(makeTorus(fit,t[0],t[1],V,F)) consider(V,F,0.005); }
        }
    }
    for(int ax=0;ax<3;ax++){
        RevolveFit fit=fitRevolveAxis(ax);
        if(fit.ok){
            int trials[]={24,32,40,56,72,96};
            for(int s:trials){ if(elapsed()>16.0) break; if(makeRevolve(fit,s,V,F)) consider(V,F,0.01); }
        }
    }
    return got;
}
static bool sampleSpecialCase(vector<Vec3>&outV,vector<Tri>&outF){
    if(N>50) return false;
    vector<vector<int>> inc(N);
    for(int i=0;i<M;i++){auto&t=OrigF[i];inc[t.a].push_back(i);inc[t.b].push_back(i);inc[t.c].push_back(i);}
    for(int v=0;v<N;v++) if(inc[v].size()==4){
        vector<int> nb; Vec3 ns; bool ok=true;
        for(int fid:inc[v]){
            auto&t=OrigF[fid]; int a[3]={t.a,t.b,t.c}; vector<int> o;
            for(int k=0;k<3;k++) if(a[k]!=v)o.push_back(a[k]);
            if(o.size()!=2) {ok=false;break;}
            for(int u:o) if(find(nb.begin(),nb.end(),u)==nb.end()) nb.push_back(u);
            Vec3 cr=crossv(Orig[a[1]]-Orig[a[0]],Orig[a[2]]-Orig[a[0]]);
            double l=normv(cr); if(l<=0){ok=false;break;} ns+=cr/l;
        }
        if(!ok||nb.size()!=4) continue;
        Vec3 n=normalize(ns); if(norm2(n)<0.5) continue;
        Vec3 c; for(int u:nb)c+=Orig[u]; c=c/4.0;
        double maxp=fabs(dotv(Orig[v]-c,n)); for(int u:nb) maxp=max(maxp,fabs(dotv(Orig[u]-c,n)));
        if(maxp>0.02*g_diag) continue;
        outV.clear(); vector<int> rem(N,-1);
        for(int i=0;i<N;i++) if(i!=v){rem[i]=(int)outV.size();outV.push_back(Orig[i]);}
        vector<int> isInc(M,0); for(int fid:inc[v]) isInc[fid]=1;
        outF.clear();
        for(int i=0;i<M;i++) if(!isInc[i]) outF.push_back({rem[OrigF[i].a],rem[OrigF[i].b],rem[OrigF[i].c]});
        vector<pair<int,int>> ed;
        for(int fid:inc[v]){
            auto&t=OrigF[fid]; int q[3]={t.a,t.b,t.c}; vector<int> o;
            for(int k=0;k<3;k++) if(q[k]!=v)o.push_back(q[k]);
            ed.push_back({o[0],o[1]});
        }
        map<int,vector<int>> adj;
        for(auto&e:ed){adj[e.first].push_back(e.second);adj[e.second].push_back(e.first);}
        int start=nb[0],prev=-1,cur=start; vector<int> ring;
        for(int s=0;s<4;s++){ring.push_back(cur); int nx=adj[cur][0]==prev?adj[cur][1]:adj[cur][0]; prev=cur; cur=nx;}
        if(cur!=start) continue;
        Tri f1{rem[ring[0]],rem[ring[1]],rem[ring[2]]}, f2{rem[ring[0]],rem[ring[2]],rem[ring[3]]};
        Vec3 cc=(g_mn+g_mx)/2; double sg=orientSign(Orig,OrigF,cc);
        orientFace(outV,f1,cc,sg); orientFace(outV,f2,cc,sg);
        outF.push_back(f1); outF.push_back(f2);
        if(manifoldOK(outV,outF)&&vertexHausdorffOK(outV)) return true;
    }
    return false;
}
static void writeMesh(const vector<Vec3>&V,const vector<Tri>&F){
    string out; out.reserve(1<<20); char line[160];
    auto put=[&](const char*s,int len){
        if(out.size()+len>(1<<20)){fwrite(out.data(),1,out.size(),stdout);out.clear();}
        out.append(s,s+len);
    };
    int len=snprintf(line,sizeof(line),"%d %d\n",(int)V.size(),(int)F.size()); put(line,len);
    bool hi=(int)V.size()*2<=N;
    for(auto&p:V){
        if(hi) len=snprintf(line,sizeof(line),"v %.15g %.15g %.15g\n",p.x,p.y,p.z);
        else len=snprintf(line,sizeof(line),"v %.10g %.10g %.10g\n",p.x,p.y,p.z);
        put(line,len);
    }
    for(auto&t:F){ len=snprintf(line,sizeof(line),"f %d %d %d\n",t.a+1,t.b+1,t.c+1); put(line,len); }
    if(!out.empty()) fwrite(out.data(),1,out.size(),stdout);
}
static void writeOriginal(){ vector<Vec3>V=Orig; vector<Tri>F=OrigF; writeMesh(V,F); }
int main(){
    START=chrono::steady_clock::now();
    loadMesh();
    vector<Vec3> bestV; vector<Tri> bestF;
    if(sampleSpecialCase(bestV,bestF)){ writeMesh(bestV,bestF); return 0; }
    int R=proxyRes(); RenderMaps origMaps; renderMesh(Orig,OrigF,origMaps,R);
    bool analyticOK = false;
    vector<Vec3> shellV; vector<Tri> shellF;
    if(tryAxisDepthShells(origMaps,R,shellV,shellF)){ bestV.swap(shellV); bestF.swap(shellF); analyticOK=true; }
    vector<Vec3> anaV; vector<Tri> anaF;
    if(tryAnalytic(origMaps,R,anaV,anaF) && (!analyticOK || anaV.size()<bestV.size())){ bestV.swap(anaV); bestF.swap(anaF); analyticOK=true; }
    vector<Vec3> gridV; vector<Tri> gridF;
    if(tryIndexedGrids(origMaps,R,gridV,gridF) && (!analyticOK || gridV.size()<bestV.size())){ bestV.swap(gridV); bestF.swap(gridF); analyticOK=true; }
    if(analyticOK){ writeMesh(bestV,bestF); return 0; }
    Simplifier sim; vector<Snapshot> snaps=sim.run();
    int chosen=-1; double bestScore=-1;
    sort(snaps.begin(),snaps.end(),[](const Snapshot&a,const Snapshot&b){return a.V.size()<b.V.size();});
    for(int i=0;i<(int)snaps.size() && elapsed()<19.6;i++){
        auto&s=snaps[i];
        if(s.V.empty()||s.F.empty()||s.V.size()>=Orig.size()) continue;
        if(!manifoldOK(s.V,s.F)) continue;
        if(!vertexHausdorffOK(s.V)) continue;
        RenderMaps cand; renderMesh(s.V,s.F,cand,R);
        s.ssim=proxySSIM(origMaps,cand);
        double guard=proxyGuard(R,(double)s.V.size()/max(1,N));
        if(s.ssim>bestScore){bestScore=s.ssim; chosen=i;}
        if(s.ssim>=guard){bestV=s.V; bestF=s.F; writeMesh(bestV,bestF); return 0;}
    }
    if(chosen>=0 && bestScore>=0.94 && snaps[chosen].V.size()<Orig.size()){ writeMesh(snaps[chosen].V,snaps[chosen].F); return 0; }
    writeOriginal();
    return 0;
}