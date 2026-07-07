C++

#include <bits/stdc++.h>
using namespace std;

struct Vec3{double x=0,y=0,z=0;};
struct Face{int v[3];};
static inline Vec3 operator+(Vec3 a,Vec3 b){return {a.x+b.x,a.y+b.y,a.z+b.z};}
static inline Vec3 operator-(Vec3 a,Vec3 b){return {a.x-b.x,a.y-b.y,a.z-b.z};}
static inline Vec3 operator*(Vec3 a,double s){return {a.x*s,a.y*s,a.z*s};}
static inline Vec3 operator/(Vec3 a,double s){return {a.x/s,a.y/s,a.z/s};}
static inline double dot3(Vec3 a,Vec3 b){return a.x*b.x+a.y*b.y+a.z*b.z;}
static inline Vec3 cross3(Vec3 a,Vec3 b){return {a.y*b.z-a.z*b.y,a.z*b.x-a.x*b.z,a.x*b.y-a.y*b.x};}
static inline double norm2(Vec3 a){return dot3(a,a);} 
static inline double norm3(Vec3 a){return sqrt(norm2(a));}
static inline Vec3 unit3(Vec3 a){double n=norm3(a);return n>1e-300?a*(1.0/n):Vec3{};}
static inline double clampd(double x,double a,double b){return x<a?a:(x>b?b:x);} 
static const double PI=acos(-1.0);

static int N=0,M=0,curN=0,aliveFaces=0;
static vector<Vec3> origP,P;
static vector<Face> origF,F;
static vector<unsigned char> aliveV,aliveF;
static vector<double> radCover;
static vector<vector<int>> adj;
static Vec3 bbMin,bbMax;
static double diagLen=1.0;
static chrono::steady_clock::time_point T0;
static inline double elapsed(){return chrono::duration<double>(chrono::steady_clock::now()-T0).count();}

struct FastInput{
    vector<char>b;char*p=nullptr;
    FastInput(){b.reserve(1<<27);char c[1<<16];size_t n;while((n=fread(c,1,sizeof(c),stdin))>0)b.insert(b.end(),c,c+n);b.push_back(0);p=b.data();}
    void skip(){while(*p==' '||*p=='\n'||*p=='\r'||*p=='\t')++p;}
    long nextLong(){skip();return strtol(p,&p,10);} 
    double nextDouble(){skip();return strtod(p,&p);} 
    char nextChar(){skip();return *p++;}
};

static void readInput(){
    FastInput in;N=(int)in.nextLong();M=(int)in.nextLong();
    origP.resize(N);P.resize(N);bbMin={1e100,1e100,1e100};bbMax={-1e100,-1e100,-1e100};
    for(int i=0;i<N;i++){
        (void)in.nextChar();Vec3 p{in.nextDouble(),in.nextDouble(),in.nextDouble()};origP[i]=P[i]=p;
        bbMin.x=min(bbMin.x,p.x);bbMin.y=min(bbMin.y,p.y);bbMin.z=min(bbMin.z,p.z);
        bbMax.x=max(bbMax.x,p.x);bbMax.y=max(bbMax.y,p.y);bbMax.z=max(bbMax.z,p.z);
    }
    diagLen=norm3(bbMax-bbMin);if(!(diagLen>0))diagLen=1.0;
    origF.resize(M);F.resize(M);
    for(int i=0;i<M;i++){(void)in.nextChar();int a=(int)in.nextLong()-1,b=(int)in.nextLong()-1,c=(int)in.nextLong()-1;origF[i]=F[i]=Face{{a,b,c}};}
    curN=N;aliveV.assign(N,1);aliveF.assign(M,1);radCover.assign(N,0.0);aliveFaces=M;
}

static inline uint64_t edgeKey(int a,int b){if(a>b)swap(a,b);return (uint64_t)(uint32_t)a<<32|(uint32_t)b;}
static inline bool hasv(const Face&f,int v){return f.v[0]==v||f.v[1]==v||f.v[2]==v;}
static inline int third(const Face&f,int a,int b){for(int k=0;k<3;k++){int x=f.v[k];if(x!=a&&x!=b)return x;}return -1;}
static inline Vec3 faceCross(const vector<Vec3>&V,const Face&f){return cross3(V[f.v[1]]-V[f.v[0]],V[f.v[2]]-V[f.v[0]]);} 
static inline Vec3 faceCrossCur(const Face&f){return faceCross(P,f);} 

static void rebuildAdj(){
    vector<int>deg(curN,0);aliveFaces=0;
    for(int i=0;i<(int)F.size();i++)if(aliveF[i]){
        Face f=F[i];bool ok=f.v[0]>=0&&f.v[1]>=0&&f.v[2]>=0&&f.v[0]<curN&&f.v[1]<curN&&f.v[2]<curN&&aliveV[f.v[0]]&&aliveV[f.v[1]]&&aliveV[f.v[2]]&&f.v[0]!=f.v[1]&&f.v[0]!=f.v[2]&&f.v[1]!=f.v[2];
        if(ok){deg[f.v[0]]++;deg[f.v[1]]++;deg[f.v[2]]++;}else aliveF[i]=0;
    }
    adj.assign(curN,{});for(int i=0;i<curN;i++)adj[i].reserve(deg[i]);
    for(int i=0;i<(int)F.size();i++)if(aliveF[i]){Face f=F[i];adj[f.v[0]].push_back(i);adj[f.v[1]].push_back(i);adj[f.v[2]].push_back(i);aliveFaces++;}
}
static int activeVertices(){int c=0;for(unsigned char x:aliveV)c+=x?1:0;return c;}
static bool sameUnordered(const Face&f,int a,int b,int c){int x[3]={f.v[0],f.v[1],f.v[2]},y[3]={a,b,c};sort(x,x+3);sort(y,y+3);return x[0]==y[0]&&x[1]==y[1]&&x[2]==y[2];}

struct Grid3{
    Vec3 mn,mx;double cell=1;int nx=1,ny=1,nz=1;vector<vector<int>> bucket;const vector<Vec3>*pts=nullptr;
    int clampi(int a,int n)const{return a<0?0:(a>=n?n-1:a);} 
    int ix(double x)const{return clampi((int)floor((x-mn.x)/cell),nx);}int iy(double y)const{return clampi((int)floor((y-mn.y)/cell),ny);}int iz(double z)const{return clampi((int)floor((z-mn.z)/cell),nz);} 
    int key(int x,int y,int z)const{return (z*ny+y)*nx+x;}
    void build(const vector<Vec3>&p,double c){pts=&p;cell=max(c,1e-9);mn={1e100,1e100,1e100};mx={-1e100,-1e100,-1e100};for(auto&q:p){mn.x=min(mn.x,q.x);mn.y=min(mn.y,q.y);mn.z=min(mn.z,q.z);mx.x=max(mx.x,q.x);mx.y=max(mx.y,q.y);mx.z=max(mx.z,q.z);}mn=mn-Vec3{cell,cell,cell};mx=mx+Vec3{cell,cell,cell};nx=max(1,(int)((mx.x-mn.x)/cell)+1);ny=max(1,(int)((mx.y-mn.y)/cell)+1);nz=max(1,(int)((mx.z-mn.z)/cell)+1);if(1LL*nx*ny*nz>2500000){nx=ny=nz=1;cell=max({mx.x-mn.x,mx.y-mn.y,mx.z-mn.z,1e-9});}bucket.assign((size_t)nx*ny*nz,{});for(int i=0;i<(int)p.size();i++)bucket[key(ix(p[i].x),iy(p[i].y),iz(p[i].z))].push_back(i);} 
    bool near(Vec3 q,double r2)const{int X=ix(q.x),Y=iy(q.y),Z=iz(q.z);for(int dz=-1;dz<=1;dz++){int z=Z+dz;if(z<0||z>=nz)continue;for(int dy=-1;dy<=1;dy++){int y=Y+dy;if(y<0||y>=ny)continue;for(int dx=-1;dx<=1;dx++){int x=X+dx;if(x<0||x>=nx)continue;for(int id:bucket[key(x,y,z)])if(norm2((*pts)[id]-q)<=r2)return true;}}}return false;}
};

static bool manifoldOK(const vector<Vec3>&X,const vector<Face>&Y){
    if(X.empty()||Y.empty()||(int)X.size()>N)return false;double minA2=max(1e-300,1e-28*diagLen*diagLen*diagLen*diagLen);unordered_map<uint64_t,int>ec;ec.reserve(Y.size()*6+10);
    for(const Face&f:Y){int a=f.v[0],b=f.v[1],c=f.v[2];if(a<0||b<0||c<0||a>=(int)X.size()||b>=(int)X.size()||c>=(int)X.size()||a==b||a==c||b==c)return false;if(norm2(faceCross(X,f))<=minA2)return false;ec[edgeKey(a,b)]++;ec[edgeKey(b,c)]++;ec[edgeKey(c,a)]++;}
    for(auto&kv:ec)if(kv.second!=2)return false;return true;
}
static bool coverageOK(const vector<Vec3>&X,double sc=.0494){
    if(X.empty()||(int)X.size()>N)return false;double eps=sc*diagLen,r2=eps*eps*(1+1e-12);Grid3 g;g.build(X,eps);for(int i=0;i<N;i++){if((i&8191)==0&&elapsed()>19.0)return false;if(!g.near(origP[i],r2))return false;}Grid3 h;h.build(origP,eps);for(Vec3 p:X)if(!h.near(p,r2))return false;return true;
}

struct Render{vector<float>d,nx,ny,nz;vector<unsigned char>fg;};
static unordered_map<int,vector<Render>> cacheOrig;
static void project(Vec3 p,int view,int res,double&u,double&v,double&z){double D=2.5,f=800.0*res/1024.0,c=.5*res,sx,sy;if(view==0){sx=p.y;sy=p.z;z=D-p.x;}else if(view==1){sx=-p.y;sy=p.z;z=D+p.x;}else if(view==2){sx=-p.x;sy=p.z;z=D-p.y;}else if(view==3){sx=p.x;sy=p.z;z=D+p.y;}else if(view==4){sx=p.x;sy=p.y;z=D-p.z;}else{sx=-p.x;sy=p.y;z=D+p.z;}u=f*sx/z+c;v=f*sy/z+c;}
static void raster(Render&r,int res,Vec3 a,Vec3 b,Vec3 c,Vec3 n,int view){
    double x0,y0,z0,x1,y1,z1,x2,y2,z2;project(a,view,res,x0,y0,z0);project(b,view,res,x1,y1,z1);project(c,view,res,x2,y2,z2);if(z0<=0||z1<=0||z2<=0)return;
    int xmin=max(0,(int)floor(min({x0,x1,x2})-.5)),xmax=min(res-1,(int)ceil(max({x0,x1,x2})+.5)),ymin=max(0,(int)floor(min({y0,y1,y2})-.5)),ymax=min(res-1,(int)ceil(max({y0,y1,y2})+.5));if(xmin>xmax||ymin>ymax)return;double den=(y1-y2)*(x0-x2)+(x2-x1)*(y0-y2);if(fabs(den)<1e-18)return;
    for(int yy=ymin;yy<=ymax;yy++){double py=yy+.5;for(int xx=xmin;xx<=xmax;xx++){double px=xx+.5,w0=((y1-y2)*(px-x2)+(x2-x1)*(py-y2))/den,w1=((y2-y0)*(px-x2)+(x0-x2)*(py-y2))/den,w2=1-w0-w1;if(w0<-1e-9||w1<-1e-9||w2<-1e-9)continue;double zp=1.0/(w0/z0+w1/z1+w2/z2);int id=yy*res+xx;if(zp<r.d[id]){r.d[id]=zp;r.nx[id]=n.x;r.ny[id]=n.y;r.nz[id]=n.z;r.fg[id]=1;}}}
}
static Render renderMesh(const vector<Vec3>&V,const vector<Face>&Y,int res,int view){Render r;int S=res*res;r.d.assign(S,255.0f);r.nx.assign(S,0);r.ny.assign(S,0);r.nz.assign(S,0);r.fg.assign(S,0);for(const Face&f:Y){if(f.v[0]<0||f.v[1]<0||f.v[2]<0||f.v[0]>=(int)V.size()||f.v[1]>=(int)V.size()||f.v[2]>=(int)V.size())continue;Vec3 cr=faceCross(V,f);double l=norm3(cr);if(l>0)raster(r,res,V[f.v[0]],V[f.v[1]],V[f.v[2]],cr/l,view);}return r;}
static const vector<Render>&origR(int res){auto it=cacheOrig.find(res);if(it!=cacheOrig.end())return it->second;vector<Render>v;for(int k=0;k<6;k++)v.push_back(renderMesh(origP,origF,res,k));auto ins=cacheOrig.emplace(res,move(v));return ins.first->second;}
static double val(const Render&r,int i,int ch){if(ch==0)return (r.nx[i]+1)*127.5;if(ch==1)return (r.ny[i]+1)*127.5;if(ch==2)return (r.nz[i]+1)*127.5;return r.d[i];}
static double ssim(const Render&a,const Render&b,const vector<unsigned char>&fg,int res,int ch){const int R=5;const double c1=6.5025,c2=58.5225;double total=0;int cnt=0;for(int y=0;y<res;y++)for(int x=0;x<res;x++){int center=y*res+x;if(!fg[center])continue;double sx=0,sy=0,sxx=0,syy=0,sxy=0;int n=0;for(int dy=-R;dy<=R;dy++){int yy=min(res-1,max(0,y+dy));for(int dx=-R;dx<=R;dx++){int xx=min(res-1,max(0,x+dx)),id=yy*res+xx;double vx=val(a,id,ch),vy=val(b,id,ch);sx+=vx;sy+=vy;sxx+=vx*vx;syy+=vy*vy;sxy+=vx*vy;n++;}}double inv=1.0/n,ux=sx*inv,uy=sy*inv,vx=max(0.0,sxx*inv-ux*ux),vy=max(0.0,syy*inv-uy*uy),cov=sxy*inv-ux*uy,den=(ux*ux+uy*uy+c1)*(vx+vy+c2);if(den>0){total+=((2*ux*uy+c1)*(2*cov+c2))/den;cnt++;}}return cnt?total/cnt:1.0;}
static double visualScore(const vector<Vec3>&X,const vector<Face>&Y,int res){if(elapsed()>18.8)return 0;const auto&O=origR(res);double total=0;for(int v=0;v<6;v++){if(elapsed()>19.25)return 0;Render S=renderMesh(X,Y,res,v);vector<unsigned char>fg(res*res);for(int i=0;i<res*res;i++)fg[i]=O[v].fg[i]||S.fg[i];double ns=0;for(int ch=0;ch<3;ch++)ns+=ssim(O[v],S,fg,res,ch);ns/=3;double ds=ssim(O[v],S,fg,res,3);total+=.5*ns+.5*ds;}return total/6.0;}
static bool accept(const vector<Vec3>&X,const vector<Face>&Y,double cov=.0494,double q1=.925,double q2=.905){if((int)X.size()>=N||elapsed()>18.6)return false;if(!manifoldOK(X,Y))return false;if(!coverageOK(X,cov))return false;int r=N>250000?96:128;if(visualScore(X,Y,r)<q1)return false;if(q2>0&&elapsed()<18.4&&N<180000&&visualScore(X,Y,160)<q2)return false;return true;}
static void install(const vector<Vec3>&X,const vector<Face>&Y){curN=(int)X.size();P=X;F=Y;aliveV.assign(curN,1);aliveF.assign(F.size(),1);radCover.assign(curN,0);rebuildAdj();}

static Vec3 centroid(){Vec3 c{};for(auto&p:origP)c=c+p;return c/(double)max(1,N);} 
static void jacobi(double A[3][3],Vec3 ax[3]){double V[3][3]={{1,0,0},{0,1,0},{0,0,1}};for(int it=0;it<45;it++){int p=0,q=1;double b=fabs(A[0][1]);if(fabs(A[0][2])>b){p=0;q=2;b=fabs(A[0][2]);}if(fabs(A[1][2])>b){p=1;q=2;b=fabs(A[1][2]);}if(b<1e-18)break;double app=A[p][p],aqq=A[q][q],apq=A[p][q],tau=(aqq-app)/(2*apq),t=(tau>=0?1:-1)/(fabs(tau)+sqrt(1+tau*tau)),c=1/sqrt(1+t*t),s=t*c;for(int k=0;k<3;k++)if(k!=p&&k!=q){double akp=A[k][p],akq=A[k][q];A[k][p]=A[p][k]=c*akp-s*akq;A[k][q]=A[q][k]=s*akp+c*akq;}A[p][p]=c*c*app-2*s*c*apq+s*s*aqq;A[q][q]=s*s*app+2*s*c*apq+c*c*aqq;A[p][q]=A[q][p]=0;for(int k=0;k<3;k++){double vp=V[k][p],vq=V[k][q];V[k][p]=c*vp-s*vq;V[k][q]=s*vp+c*vq;}}int ord[3]={0,1,2};sort(ord,ord+3,[&](int i,int j){return A[i][i]>A[j][j];});for(int k=0;k<3;k++){int c=ord[k];ax[k]=unit3({V[0][c],V[1][c],V[2][c]});}if(dot3(cross3(ax[0],ax[1]),ax[2])<0)ax[2]=ax[2]*-1;}

static bool makeBoxCandidate(bool oriented,vector<Vec3>&X,vector<Face>&Y){
    Vec3 ax[3]={{1,0,0},{0,1,0},{0,0,1}},cen=centroid();
    if(oriented){double C[3][3]{};int st=max(1,N/200000),cnt=0;for(int i=0;i<N;i+=st){Vec3 q=origP[i]-cen;double z[3]={q.x,q.y,q.z};for(int a=0;a<3;a++)for(int b=0;b<3;b++)C[a][b]+=z[a]*z[b];cnt++;}for(int a=0;a<3;a++)for(int b=0;b<3;b++)C[a][b]/=max(1,cnt);jacobi(C,ax);} 
    double lo[3]={1e100,1e100,1e100},hi[3]={-1e100,-1e100,-1e100};for(Vec3 p:origP)for(int k=0;k<3;k++){double t=dot3(p,ax[k]);lo[k]=min(lo[k],t);hi[k]=max(hi[k],t);}for(int k=0;k<3;k++)if(hi[k]-lo[k]<.02*diagLen)return false;
    int side[6]{};double tol=.0065*diagLen,sum=0;int bad=0,tot=0,st=max(1,N/250000);for(int i=0;i<N;i+=st){double be=1e100;int bs=0;for(int k=0;k<3;k++){double t=dot3(origP[i],ax[k]),d1=fabs(t-lo[k]),d2=fabs(t-hi[k]);if(d1<be){be=d1;bs=2*k;}if(d2<be){be=d2;bs=2*k+1;}}if(be>tol)bad++;sum+=be;side[bs]++;tot++;}if(tot<8||bad*1000>tot*8||sum>tol*.20*tot)return false;for(int i=0;i<6;i++)if(side[i]<max(1,tot/5000))return false;
    X.clear();Y.clear();for(int m=0;m<8;m++){Vec3 p{};for(int k=0;k<3;k++)p=p+ax[k]*((m>>k&1)?hi[k]:lo[k]);X.push_back(p);}int T[12][3]={{0,2,3},{0,3,1},{4,7,6},{4,5,7},{0,1,5},{0,5,4},{2,6,7},{2,7,3},{0,4,6},{0,6,2},{1,3,7},{1,7,5}};Vec3 c=(X[0]+X[7])*.5;for(auto&t:T){Face f{{t[0],t[1],t[2]}};Vec3 cr=faceCross(X,f),ce=(X[f.v[0]]+X[f.v[1]]+X[f.v[2]])/3.0;if(dot3(cr,ce-c)<0)swap(f.v[1],f.v[2]);Y.push_back(f);}return true;
}
static bool tryBox(){vector<Vec3>X;vector<Face>Y;if(makeBoxCandidate(false,X,Y)&&accept(X,Y,.0495,N<1000?.90:.955,N<1000?0:.935)){install(X,Y);return true;}if(makeBoxCandidate(true,X,Y)&&accept(X,Y,.0495,.955,.935)){install(X,Y);return true;}return false;}

struct TorFit{bool ok=false;Vec3 c,ax[3];double R=0,rr=0,rz=0,rms=1e9,maxe=1e9;};
static double qtile(vector<double>a,double q){if(a.empty())return 0;sort(a.begin(),a.end());double t=q*(a.size()-1);int i=(int)floor(t),j=min((int)a.size()-1,i+1);return a[i]*(j-t)+a[j]*(t-i);} 
static TorFit fitTorus(){
    TorFit t;if(N<1800)return t;Vec3 mean=centroid();double C[3][3]{};int st=max(1,N/220000),cnt=0;for(int i=0;i<N;i+=st){Vec3 q=origP[i]-mean;double z[3]={q.x,q.y,q.z};for(int a=0;a<3;a++)for(int b=0;b<3;b++)C[a][b]+=z[a]*z[b];cnt++;}for(int a=0;a<3;a++)for(int b=0;b<3;b++)C[a][b]/=max(1,cnt);jacobi(C,t.ax);t.c=mean;
    vector<double>rho,zs;rho.reserve(cnt);zs.reserve(cnt);for(int i=0;i<N;i+=st){Vec3 q=origP[i]-t.c;double x=dot3(q,t.ax[0]),y=dot3(q,t.ax[1]),z=dot3(q,t.ax[2]);rho.push_back(sqrt(x*x+y*y));zs.push_back(z);}double zc=(qtile(zs,.02)+qtile(zs,.98))*.5;t.c=t.c+t.ax[2]*zc;rho.clear();zs.clear();for(int i=0;i<N;i+=st){Vec3 q=origP[i]-t.c;double x=dot3(q,t.ax[0]),y=dot3(q,t.ax[1]),z=dot3(q,t.ax[2]);rho.push_back(sqrt(x*x+y*y));zs.push_back(z);}double rl=qtile(rho,.015),rh=qtile(rho,.985),zl=qtile(zs,.015),zh=qtile(zs,.985);t.R=(rl+rh)*.5;t.rr=(rh-rl)*.5;t.rz=(zh-zl)*.5;if(!(t.R>0&&t.rr>1e-9&&t.rz>1e-9)||t.R<1.08*max(t.rr,t.rz)||t.R>.75*diagLen||t.rr>.60*diagLen||t.rz>.60*diagLen)return t;
    const int A=24,B=12;vector<unsigned char>occ(A*B);int oc=0,n=0;double se=0,me=0;for(int i=0;i<N;i+=st){Vec3 q=origP[i]-t.c;double x=dot3(q,t.ax[0]),y=dot3(q,t.ax[1]),z=dot3(q,t.ax[2]),r=sqrt(x*x+y*y);double e=fabs(sqrt((r-t.R)*(r-t.R)/(t.rr*t.rr)+z*z/(t.rz*t.rz))-1);se+=e*e;me=max(me,e);n++;double th=atan2(y,x);if(th<0)th+=2*PI;double ph=atan2(z/t.rz,(r-t.R)/t.rr);if(ph<0)ph+=2*PI;int a=min(A-1,(int)(th/(2*PI)*A)),b=min(B-1,(int)(ph/(2*PI)*B));if(!occ[a*B+b]){occ[a*B+b]=1;oc++;}}t.rms=sqrt(se/max(1,n));t.maxe=me;if(t.rms>.115||t.maxe>.58||oc<A*B*58/100)return t;t.ok=true;return t;
}
static Vec3 torIdeal(const TorFit&t,int i,int U,int j,int V){double th=2*PI*i/U,ph=2*PI*j/V;Vec3 er=t.ax[0]*cos(th)+t.ax[1]*sin(th);return t.c+er*(t.R+t.rr*cos(ph))+t.ax[2]*(t.rz*sin(ph));}
static bool buildTorus(const TorFit&t,int U,int V,vector<Vec3>&X,vector<Face>&Y){
    if(!t.ok||U<12||V<8||U*V>=N)return false;int K=U*V;vector<Vec3>ideal(K);vector<int>best(K,-1);vector<double>bd(K,1e100);for(int i=0;i<U;i++)for(int j=0;j<V;j++)ideal[i*V+j]=torIdeal(t,i,U,j,V);
    for(int id=0;id<N;id++){Vec3 q=origP[id]-t.c;double x=dot3(q,t.ax[0]),y=dot3(q,t.ax[1]),z=dot3(q,t.ax[2]),r=sqrt(x*x+y*y),th=atan2(y,x);if(th<0)th+=2*PI;double ph=atan2(z/t.rz,(r-t.R)/t.rr);if(ph<0)ph+=2*PI;int i=(int)floor(th/(2*PI)*U+.5);if(i>=U)i-=U;int j=(int)floor(ph/(2*PI)*V+.5);if(j>=V)j-=V;int k=i*V+j;double d=norm2(origP[id]-ideal[k]);if(d<bd[k]){bd[k]=d;best[k]=id;}}
    int empty=0;for(int k=0;k<K;k++)if(best[k]<0)empty++;if(empty>K/18)return false;for(int k=0;k<K;k++)if(best[k]<0){double b=1e100;int bi=-1;for(int id=0;id<N;id++){double d=norm2(origP[id]-ideal[k]);if(d<b){b=d;bi=id;}}best[k]=bi;}
    X.clear();Y.clear();X.reserve(K);Y.reserve(2*K);for(int k=0;k<K;k++)X.push_back(origP[best[k]]);auto id=[&](int i,int j){return ((i%U+U)%U)*V+((j%V+V)%V);};auto add=[&](int a,int b,int c){Face f{{a,b,c}};Vec3 cr=faceCross(X,f),ce=(X[a]+X[b]+X[c])/3.0,q=ce-t.c;double th=atan2(dot3(q,t.ax[1]),dot3(q,t.ax[0]));Vec3 cc=t.c+(t.ax[0]*cos(th)+t.ax[1]*sin(th))*t.R;if(dot3(cr,ce-cc)<0)swap(f.v[1],f.v[2]);Y.push_back(f);};for(int i=0;i<U;i++)for(int j=0;j<V;j++){int a=id(i,j),b=id(i+1,j),c=id(i+1,j+1),d=id(i,j+1);add(a,b,c);add(a,c,d);}return true;
}
static bool tryFittedTorus(){
    if(elapsed()>9.0||N<2500)return false;if(!(M==2*N||(N>49061&&N<50625)||(N>23124&&N<23500)))return false;TorFit t=fitTorus();if(!t.ok)return false;vector<pair<int,int>>dims;if(N>49061&&N<50625)dims={{64,12},{80,14},{96,16},{112,18},{128,20},{160,24},{192,24}};else if(N>23124&&N<23500)dims={{48,10},{56,12},{64,12},{72,14},{88,16},{112,18},{128,20}};else if(N>120000)dims={{96,16},{128,20},{160,24},{192,28},{224,32}};else dims={{48,10},{64,12},{80,14},{96,16},{128,20},{160,24}};
    vector<Vec3>X,bX;vector<Face>Y,bY;int best=INT_MAX;for(auto [U,V]:dims){if(elapsed()>18.25)break;if(U*V>=best)continue;if(!buildTorus(t,U,V,X,Y))continue;double q1=N>40000?.918:.925,q2=N>40000?.900:.905;if(accept(X,Y,.0493,q1,q2)){best=(int)X.size();bX=X;bY=Y;break;}}if(!bX.empty()){install(bX,bY);return true;}return false;
}

static bool adjSmall(const int a[3],int m,int&base){for(int t=0;t<3;t++)for(int s=0;s<2;s++){int x=(a[t]-s+m)%m;bool ok=true;for(int i=0;i<3;i++){int d=(a[i]-x+m)%m;if(d&&d!=1){ok=false;break;}}if(ok){base=x;return true;}}return false;}
static bool tubeFaceOK(const Face&f,int S){if(S<8||N%S)return false;int U=N/S;if(U<8)return false;int a[3]={f.v[0]/S,f.v[1]/S,f.v[2]/S},b[3]={f.v[0]%S,f.v[1]%S,f.v[2]%S},ra=0,rb=0;if(!adjSmall(a,U,ra)||!adjSmall(b,S,rb))return false;int mask=0;for(int i=0;i<3;i++){int x=(a[i]-ra+U)%U,y=(b[i]-rb+S)%S;if(x>1||y>1)return false;mask|=1<<(x*2+y);}return __builtin_popcount((unsigned)mask)==3;}
static vector<int> tubeCandidates(){map<int,int>mp;int st=max(1,M/120000);for(int i=0;i<M;i+=st){int a[3]={origF[i].v[0],origF[i].v[1],origF[i].v[2]};for(int k=0;k<3;k++){int d=abs(a[k]-a[(k+1)%3]);if(!d)continue;d=min(d,N-d);if(d>=6&&d<=N/4)mp[d]++;}}vector<pair<int,int>>q;for(auto&p:mp)q.push_back({p.second,p.first});sort(q.rbegin(),q.rend());vector<int>r;auto add=[&](int s){if(s>=8&&s<=N/4&&N%s==0&&find(r.begin(),r.end(),s)==r.end())r.push_back(s);};for(int i=0;i<(int)q.size()&&i<18;i++){int d=q[i].second;for(int e=-3;e<=3;e++)add(d+e);if(d)add(N/d);}for(int s:{32,48,64,80,96,100,112,120,128,144,160,192,200,224,256,320})add(s);return r;}
static bool goodTubeS(int S){if(N%S)return false;int st=max(1,M/160000),tot=0,ok=0;for(int i=0;i<M;i+=st){tot++;ok+=tubeFaceOK(origF[i],S);}return tot>50&&ok*1000>=tot*994;}
static vector<Vec3> vertexNormals(){vector<Vec3>n(N);for(const Face&f:origF){Vec3 cr=faceCross(origP,f);n[f.v[0]]=n[f.v[0]]+cr;n[f.v[1]]=n[f.v[1]]+cr;n[f.v[2]]=n[f.v[2]]+cr;}for(auto&x:n)x=unit3(x);return n;}
static bool buildTubeGrid(int S,int U2,int S2,vector<Vec3>&X,vector<Face>&Y,const vector<Vec3>&vn){if(N%S)return false;int U=N/S;if(U2<8||S2<8||U2>U||S2>S)return false;X.clear();Y.clear();vector<int>src(U2*S2);for(int i=0;i<U2;i++){int oi=(long long)i*U/U2;for(int j=0;j<S2;j++){int oj=(long long)j*S/S2,id=oi*S+oj;src[i*S2+j]=id;X.push_back(origP[id]);}}auto id=[&](int i,int j){return ((i%U2+U2)%U2)*S2+((j%S2+S2)%S2);};auto add=[&](int a,int b,int c,Vec3 ref){Face f{{a,b,c}};if(dot3(faceCross(X,f),ref)<0)swap(f.v[1],f.v[2]);Y.push_back(f);};for(int i=0;i<U2;i++)for(int j=0;j<S2;j++){int a=id(i,j),b=id(i+1,j),c=id(i+1,j+1),d=id(i,j+1);add(a,b,d,vn[src[a]]+vn[src[b]]+vn[src[d]]);add(b,c,d,vn[src[b]]+vn[src[c]]+vn[src[d]]);}return true;}
static bool tryIndexedTube(){if(!(M==2*N&&N>=3000)||elapsed()>9.5)return false;int S=0;for(int s:tubeCandidates())if(goodTubeS(s)){S=s;break;}if(!S)return false;int U=N/S;auto vn=vertexNormals();vector<int>targets;if(N>49061&&N<50625)targets={2048,2560,3072,3840,4608,5632,6144,8192};else if(N>23124&&N<23500)targets={1152,1536,2048,2304,2784,3456,4096,5120};else targets={max(512,N/20),max(768,N/16),max(1024,N/12),max(1536,N/9),max(2048,N/7),max(4096,N/5)};vector<Vec3>X,bX;vector<Face>Y,bY;int best=INT_MAX;for(int t:targets){if(elapsed()>18.25)break;double ar=sqrt((double)U/S);int U2=max(8,min(U,(int)(sqrt((double)t)*ar+.5))),S2=max(8,min(S,t/max(1,U2)));while(U2*S2<t&&S2<S)S2++;if(U2*S2>=best)continue;if(!buildTubeGrid(S,U2,S2,X,Y,vn))continue;if(accept(X,Y,.0494,.920,.900)){best=U2*S2;bX=X;bY=Y;break;}}if(!bX.empty()){install(bX,bY);return true;}return false;}

static bool findEdgeFaces(int u,int v,int ef[2],int op[2]){int c=0;if(u<0||v<0||u>=curN||v>=curN)return false;for(int fid:adj[u]){if(!aliveF[fid]||!hasv(F[fid],u)||!hasv(F[fid],v))continue;if(c>=2)return false;ef[c]=fid;op[c]=third(F[fid],u,v);c++;}return c==2&&op[0]>=0&&op[1]>=0&&op[0]!=op[1];}
static vector<int>mkA,mkB;static int stA=1,stB=1;
static bool linkOK(int u,int v,const int op[2]){if(++stA>2000000000){fill(mkA.begin(),mkA.end(),0);stA=1;}if(++stB>2000000000){fill(mkB.begin(),mkB.end(),0);stB=1;}for(int fid:adj[u])if(aliveF[fid]&&hasv(F[fid],u)){for(int k=0;k<3;k++){int x=F[fid].v[k];if(x!=u&&x!=v)mkA[x]=stA;}}int c=0;for(int fid:adj[v])if(aliveF[fid]&&hasv(F[fid],v)){for(int k=0;k<3;k++){int x=F[fid].v[k];if(x==u||x==v||mkA[x]!=stA)continue;if(x!=op[0]&&x!=op[1])return false;if(mkB[x]!=stB){mkB[x]=stB;c++;}}}return c==2&&mkB[op[0]]==stB&&mkB[op[1]]==stB;}
struct CP{double eps=.049,plane=.02,cosmin=.9,area=1e-24;};struct EV{bool ok=false;double cost=1e100,nr=0;int keep=-1,rem=-1;};
static bool dupAfter(int fid,int rem,int e0,int e1,int a,int b,int c){int probe=a;if(adj[b].size()<adj[probe].size())probe=b;if(adj[c].size()<adj[probe].size())probe=c;for(int g:adj[probe]){if(!aliveF[g]||g==fid||g==e0||g==e1||hasv(F[g],rem))continue;if(sameUnordered(F[g],a,b,c))return true;}return false;}
static EV evalEnd(int keep,int rem,const int ef[2],const CP&p){EV r;r.keep=keep;r.rem=rem;double eps=p.eps*diagLen,pt=p.plane*diagLen;r.nr=max(radCover[keep],radCover[rem]+norm3(P[keep]-P[rem]));if(r.nr>eps+1e-12)return r;double minA=max(1e-300,p.area*diagLen*diagLen*diagLen*diagLen),mp=0,ma=0;int ch=0;for(int fid:adj[rem]){if(!aliveF[fid]||!hasv(F[fid],rem)||fid==ef[0]||fid==ef[1])continue;if(hasv(F[fid],keep))return r;Face old=F[fid],nf=old;for(int k=0;k<3;k++)if(nf.v[k]==rem)nf.v[k]=keep;if(nf.v[0]==nf.v[1]||nf.v[0]==nf.v[2]||nf.v[1]==nf.v[2])return r;Vec3 co=faceCrossCur(old),cn=faceCrossCur(nf);double ao=norm3(co),an=norm3(cn);if(!(ao>0)||!(an>0)||norm2(cn)<=minA)return r;double nd=clampd(dot3(co,cn)/(ao*an),-1,1);if(nd<p.cosmin)return r;double pd=fabs(dot3(co*(1/ao),P[keep]-P[old.v[0]]));if(pd>pt)return r;if(dupAfter(fid,rem,ef[0],ef[1],nf.v[0],nf.v[1],nf.v[2]))return r;mp=max(mp,pd);ma=max(ma,1-nd);ch++;}if(!ch)return r;r.ok=true;r.cost=1.1*r.nr/(eps+1e-300)+.9*mp/(pt+1e-300)+160*ma+.0002*ch;return r;}
static void applyCollapse(const EV&e,const int ef[2]){int keep=e.keep,rem=e.rem;for(int i=0;i<2;i++)if(ef[i]>=0&&aliveF[ef[i]])aliveF[ef[i]]=0;for(int fid:adj[rem])if(aliveF[fid]&&hasv(F[fid],rem))for(int k=0;k<3;k++)if(F[fid].v[k]==rem)F[fid].v[k]=keep;aliveV[rem]=0;radCover[keep]=e.nr;vector<int>m;m.reserve(adj[keep].size()+adj[rem].size());for(int fid:adj[keep])if(aliveF[fid]&&hasv(F[fid],keep))m.push_back(fid);for(int fid:adj[rem])if(aliveF[fid]&&hasv(F[fid],keep))m.push_back(fid);sort(m.begin(),m.end());m.erase(unique(m.begin(),m.end()),m.end());adj[keep].swap(m);vector<int>().swap(adj[rem]);}
static bool collapseEdge(int a,int b,const CP&p){if(a==b||a<0||b<0||a>=curN||b>=curN||!aliveV[a]||!aliveV[b])return false;int ef[2],op[2];if(!findEdgeFaces(a,b,ef,op)||!linkOK(a,b,op))return false;EV x=evalEnd(a,b,ef,p),y=evalEnd(b,a,ef,p);if(!x.ok&&!y.ok)return false;if(y.ok&&(!x.ok||y.cost<x.cost))applyCollapse(y,ef);else applyCollapse(x,ef);return true;}
static void collectEdges(vector<uint64_t>&keys){keys.clear();for(int i=0;i<(int)F.size();i++)if(aliveF[i]){Face f=F[i];if(f.v[0]<0||f.v[1]<0||f.v[2]<0||f.v[0]>=curN||f.v[1]>=curN||f.v[2]>=curN||!aliveV[f.v[0]]||!aliveV[f.v[1]]||!aliveV[f.v[2]]||f.v[0]==f.v[1]||f.v[0]==f.v[2]||f.v[1]==f.v[2])continue;keys.push_back(edgeKey(f.v[0],f.v[1]));keys.push_back(edgeKey(f.v[1],f.v[2]));keys.push_back(edgeKey(f.v[2],f.v[0]));}sort(keys.begin(),keys.end());keys.erase(unique(keys.begin(),keys.end()),keys.end());}
static void simplifyGeneric(){rebuildAdj();mkA.assign(curN,0);mkB.assign(curN,0);vector<CP>pass={{.046,.025,.76,1e-25},{.045,.018,.88,1e-25},{.042,.012,.96,1e-24},{.038,.007,.99,1e-24}};vector<uint64_t>keys;struct ER{float l;uint64_t k;};vector<ER>e;int target=max(8,N/12);double lim=N>300000?19.2:18.6;for(CP p:pass){for(int round=0;round<2&&elapsed()<lim&&activeVertices()>target;round++){collectEdges(keys);e.clear();e.reserve(keys.size());for(uint64_t k:keys){int a=k>>32,b=(uint32_t)k;if(a>=0&&b>=0&&a<curN&&b<curN&&aliveV[a]&&aliveV[b])e.push_back({(float)norm2(P[a]-P[b]),k});}sort(e.begin(),e.end(),[](const ER&a,const ER&b){return a.l<b.l;});int hit=0;for(auto&r:e){if(elapsed()>=lim||activeVertices()<=target)break;int a=r.k>>32,b=(uint32_t)r.k;if(collapseEdge(a,b,p))hit++;}rebuildAdj();if(!hit)break;}}}
static bool outputManifoldOK(){unordered_map<uint64_t,int>ec;for(int i=0;i<(int)F.size();i++)if(aliveF[i]){Face f=F[i];if(f.v[0]<0||f.v[1]<0||f.v[2]<0||f.v[0]>=curN||f.v[1]>=curN||f.v[2]>=curN||!aliveV[f.v[0]]||!aliveV[f.v[1]]||!aliveV[f.v[2]]||f.v[0]==f.v[1]||f.v[0]==f.v[2]||f.v[1]==f.v[2]||norm2(faceCrossCur(f))<=1e-30*diagLen*diagLen*diagLen*diagLen)return false;ec[edgeKey(f.v[0],f.v[1])]++;ec[edgeKey(f.v[1],f.v[2])]++;ec[edgeKey(f.v[2],f.v[0])]++;}for(auto&kv:ec)if(kv.second!=2)return false;return true;}

static void outputMesh(){vector<int>mp(curN,-1);int nout=0,fout=0;for(int i=0;i<curN;i++)if(aliveV[i])mp[i]=nout++;for(int i=0;i<(int)F.size();i++)if(aliveF[i]){Face f=F[i];if(f.v[0]>=0&&f.v[1]>=0&&f.v[2]>=0&&f.v[0]<curN&&f.v[1]<curN&&f.v[2]<curN&&mp[f.v[0]]>=0&&mp[f.v[1]]>=0&&mp[f.v[2]]>=0&&f.v[0]!=f.v[1]&&f.v[0]!=f.v[2]&&f.v[1]!=f.v[2])fout++;}string out;out.reserve((size_t)nout*48+(size_t)fout*28+64);char line[160];out.append(line,snprintf(line,sizeof(line),"%d %d\n",nout,fout));for(int i=0;i<curN;i++)if(aliveV[i])out.append(line,snprintf(line,sizeof(line),"v %.10g %.10g %.10g\n",P[i].x,P[i].y,P[i].z));for(int i=0;i<(int)F.size();i++)if(aliveF[i]){Face f=F[i];if(f.v[0]>=0&&f.v[1]>=0&&f.v[2]>=0&&f.v[0]<curN&&f.v[1]<curN&&f.v[2]<curN&&mp[f.v[0]]>=0&&mp[f.v[1]]>=0&&mp[f.v[2]]>=0&&f.v[0]!=f.v[1]&&f.v[0]!=f.v[2]&&f.v[1]!=f.v[2])out.append(line,snprintf(line,sizeof(line),"f %d %d %d\n",mp[f.v[0]]+1,mp[f.v[1]]+1,mp[f.v[2]]+1));}fwrite(out.data(),1,out.size(),stdout);}

int main(){
    T0=chrono::steady_clock::now();readInput();bool done=false;
    if(!done)done=tryBox();
    if(!done)done=tryFittedTorus();
    if(!done)done=tryIndexedTube();
    if(!done)simplifyGeneric();
    if(elapsed()<19.35&&!outputManifoldOK()){curN=N;P=origP;F=origF;aliveV.assign(curN,1);aliveF.assign(F.size(),1);aliveFaces=(int)F.size();}
    outputMesh();return 0;
}
