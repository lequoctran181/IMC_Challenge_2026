#include <bits/stdc++.h>
using namespace std;

struct Vec{double x,y,z;};
static inline Vec operator+(const Vec&a,const Vec&b){return {a.x+b.x,a.y+b.y,a.z+b.z};}
static inline Vec operator-(const Vec&a,const Vec&b){return {a.x-b.x,a.y-b.y,a.z-b.z};}
static inline Vec operator*(const Vec&a,double s){return {a.x*s,a.y*s,a.z*s};}
static inline Vec operator/(const Vec&a,double s){return {a.x/s,a.y/s,a.z/s};}
static inline double dotv(const Vec&a,const Vec&b){return a.x*b.x+a.y*b.y+a.z*b.z;}
static inline Vec crossv(const Vec&a,const Vec&b){return {a.y*b.z-a.z*b.y,a.z*b.x-a.x*b.z,a.x*b.y-a.y*b.x};}
static inline double norm2(const Vec&a){return dotv(a,a);} 
static inline double normv(const Vec&a){return sqrt(max(0.0,norm2(a)));}
static inline Vec unit(Vec a){double n=normv(a); if(n<=0) return {0,0,0}; return a/n;}
struct Face{int a,b,c;};

struct Scanner{
    string s; size_t p=0;
    Scanner(){ ios::sync_with_stdio(false); cin.tie(nullptr); string part; std::ostringstream oss; oss<<cin.rdbuf(); s=oss.str(); }
    string next(){ while(p<s.size() && isspace((unsigned char)s[p])) p++; size_t q=p; while(q<s.size() && !isspace((unsigned char)s[q])) q++; string t=s.substr(p,q-p); p=q; return t; }
    bool empty(){ while(p<s.size() && isspace((unsigned char)s[p])) p++; return p>=s.size(); }
};
static int parseIndexToken(const string&t){ const char* c=t.c_str(); char* e=nullptr; long v=strtol(c,&e,10); return (int)v; }
static bool isNumTok(const string&t){ if(t.empty()) return false; char c=t[0]; return (c=='-'||c=='+'||c=='.'||(c>='0'&&c<='9')); }

static int N,M; static vector<Vec> origV; static vector<Face> origF; static Vec bbMin,bbMax; static double diagLen, epsH; static chrono::steady_clock::time_point t0;
static double elapsed(){return chrono::duration<double>(chrono::steady_clock::now()-t0).count();}

static void readInput(){
    Scanner sc; string t=sc.next(); if(t.empty()) exit(0); N=stoi(t); M=stoi(sc.next());
    origV.resize(N); origF.resize(M); bbMin={1e100,1e100,1e100}; bbMax={-1e100,-1e100,-1e100};
    for(int i=0;i<N;i++){
        string a=sc.next(); string sx,sy,sz; if(isNumTok(a)){sx=a; sy=sc.next(); sz=sc.next();} else {sx=sc.next(); sy=sc.next(); sz=sc.next();}
        Vec p{stod(sx),stod(sy),stod(sz)}; origV[i]=p;
        bbMin.x=min(bbMin.x,p.x); bbMin.y=min(bbMin.y,p.y); bbMin.z=min(bbMin.z,p.z);
        bbMax.x=max(bbMax.x,p.x); bbMax.y=max(bbMax.y,p.y); bbMax.z=max(bbMax.z,p.z);
    }
    for(int i=0;i<M;i++){
        string a=sc.next(); string sb,scs,sd; if(isNumTok(a)){sb=a; scs=sc.next(); sd=sc.next();} else {sb=sc.next(); scs=sc.next(); sd=sc.next();}
        origF[i]={parseIndexToken(sb)-1,parseIndexToken(scs)-1,parseIndexToken(sd)-1};
    }
    diagLen=normv(bbMax-bbMin); if(!(diagLen>0)) diagLen=1; epsH=0.05*diagLen;
}

static uint64_t ekey(int a,int b){ if(a>b) swap(a,b); return (uint64_t)(uint32_t)a<<32 | (uint32_t)b; }
static bool validateMesh(const vector<Vec>&V,const vector<Face>&F){
    if(V.size()<4 || F.size()<4 || V.size()>(size_t)max(1,N)) return false;
    double minA2=max(1e-32,1e-26*diagLen*diagLen*diagLen*diagLen);
    vector<uint64_t> edges; edges.reserve(F.size()*3);
    for(const Face&f:F){
        if(f.a<0||f.b<0||f.c<0||f.a>=(int)V.size()||f.b>=(int)V.size()||f.c>=(int)V.size()) return false;
        if(f.a==f.b||f.a==f.c||f.b==f.c) return false;
        if(norm2(crossv(V[f.b]-V[f.a],V[f.c]-V[f.a]))<=minA2) return false;
        edges.push_back(ekey(f.a,f.b)); edges.push_back(ekey(f.b,f.c)); edges.push_back(ekey(f.c,f.a));
    }
    sort(edges.begin(),edges.end());
    for(size_t i=0;i<edges.size();){ size_t j=i+1; while(j<edges.size()&&edges[j]==edges[i]) j++; if(j-i!=2) return false; i=j; }
    return true;
}
static void orientOut(vector<Vec>&V, vector<Face>&F){
    Vec c{0,0,0}; for(auto&p:V)c=c+p; c=c/(double)max<size_t>(1,V.size());
    for(auto &f:F){ Vec cr=crossv(V[f.b]-V[f.a],V[f.c]-V[f.a]); Vec fc=(V[f.a]+V[f.b]+V[f.c])/3.0; if(dotv(cr,fc-c)<0) swap(f.b,f.c); }
}
static void compactUsed(vector<Vec>&V, vector<Face>&F){
    vector<int> mp(V.size(),-1); int n=0; for(auto&f:F){ if(mp[f.a]<0) mp[f.a]=n++; if(mp[f.b]<0) mp[f.b]=n++; if(mp[f.c]<0) mp[f.c]=n++; }
    vector<Vec> NV(n); for(size_t i=0;i<V.size();++i) if(mp[i]>=0) NV[mp[i]]=V[i];
    for(auto&f:F){f.a=mp[f.a]; f.b=mp[f.b]; f.c=mp[f.c];}
    V.swap(NV);
}

struct KDTree{
    struct Node{int id,l=-1,r=-1,ax; double split;};
    const vector<Vec>*P=nullptr; vector<int> idx; vector<Node> tr;
    int build(int l,int r){ if(l>=r) return -1; Vec mn{1e100,1e100,1e100}, mx{-1e100,-1e100,-1e100}; for(int i=l;i<r;i++){const Vec&p=(*P)[idx[i]]; mn.x=min(mn.x,p.x); mn.y=min(mn.y,p.y); mn.z=min(mn.z,p.z); mx.x=max(mx.x,p.x); mx.y=max(mx.y,p.y); mx.z=max(mx.z,p.z);} int ax=0; Vec d=mx-mn; if(d.y>d.x&&d.y>=d.z) ax=1; else if(d.z>d.x&&d.z>=d.y) ax=2; int m=(l+r)/2; nth_element(idx.begin()+l,idx.begin()+m,idx.begin()+r,[&](int a,int b){const Vec&A=(*P)[a],&B=(*P)[b]; return ax==0?A.x<B.x:ax==1?A.y<B.y:A.z<B.z;}); int node=tr.size(); tr.push_back({idx[m],-1,-1,ax,0}); const Vec&p=(*P)[idx[m]]; tr[node].split=(ax==0?p.x:ax==1?p.y:p.z); tr[node].l=build(l,m); tr[node].r=build(m+1,r); return node; }
    void init(const vector<Vec>&pts){P=&pts; idx.resize(pts.size()); iota(idx.begin(),idx.end(),0); tr.clear(); tr.reserve(pts.size()); build(0,idx.size());}
    void queryRec(int n,const Vec&q,double&best) const{ if(n<0) return; const Node&nd=tr[n]; const Vec&p=(*P)[nd.id]; double d=norm2(p-q); if(d<best) best=d; double val=nd.ax==0?q.x:nd.ax==1?q.y:q.z; double diff=val-nd.split; int first=diff<0?nd.l:nd.r, second=diff<0?nd.r:nd.l; queryRec(first,q,best); if(diff*diff<best) queryRec(second,q,best); }
    double nearest2(const Vec&q) const{ if(tr.empty()) return 1e300; double b=1e300; queryRec(0,q,b); return b; }
};
static bool hausdorffKD(const vector<Vec>&V,double eps){
    if(V.empty()) return false; double e2=eps*eps;
    KDTree kdO,kdC; kdC.init(V); for(const Vec&p:origV){ if(kdC.nearest2(p)>e2*1.000001) return false; }
    kdO.init(origV); for(const Vec&p:V){ if(kdO.nearest2(p)>e2*1.000001) return false; }
    return true;
}

struct Render{int R; vector<double>d; vector<Vec>n; vector<unsigned char>m;};
static inline double coordAxis(const Vec&p,int ax){return ax==0?p.x:ax==1?p.y:p.z;}
static void projectPoint(const Vec&p,int view,int R,double &x,double&y,double&dep){
    int ax=view/2; bool pos=(view%2==0); int uax=(ax==0?1:0), vax=(ax==2?1:2); if(ax==1){uax=0;vax=2;} if(ax==2){uax=0;vax=1;}
    double u0=coordAxis(bbMin,uax), u1=coordAxis(bbMax,uax), v0=coordAxis(bbMin,vax), v1=coordAxis(bbMax,vax); double ur=max(1e-12,u1-u0), vr=max(1e-12,v1-v0);
    x=(coordAxis(p,uax)-u0)/ur*(R-1); y=(coordAxis(p,vax)-v0)/vr*(R-1);
    double ar=max(1e-12,coordAxis(bbMax,ax)-coordAxis(bbMin,ax)); dep=(pos?(coordAxis(bbMax,ax)-coordAxis(p,ax)):(coordAxis(p,ax)-coordAxis(bbMin,ax)))/ar*255.0; if(dep<0)dep=0; if(dep>255)dep=255;
}
static void rastTri(Render&rm,const Vec&A,const Vec&B,const Vec&C,const Vec&normal,int view){
    double x0,y0,z0,x1,y1,z1,x2,y2,z2; int R=rm.R; projectPoint(A,view,R,x0,y0,z0); projectPoint(B,view,R,x1,y1,z1); projectPoint(C,view,R,x2,y2,z2);
    int minx=max(0,(int)floor(min(x0,min(x1,x2)))); int maxx=min(R-1,(int)ceil(max(x0,max(x1,x2)))); int miny=max(0,(int)floor(min(y0,min(y1,y2)))); int maxy=min(R-1,(int)ceil(max(y0,max(y1,y2)))); if(minx>maxx||miny>maxy) return;
    double den=(y1-y2)*(x0-x2)+(x2-x1)*(y0-y2); if(fabs(den)<1e-12) return;
    for(int yy=miny;yy<=maxy;yy++){ double py=yy+0.5; for(int xx=minx;xx<=maxx;xx++){ double px=xx+0.5; double w0=((y1-y2)*(px-x2)+(x2-x1)*(py-y2))/den; double w1=((y2-y0)*(px-x2)+(x0-x2)*(py-y2))/den; double w2=1-w0-w1; if(w0<-1e-9||w1<-1e-9||w2<-1e-9) continue; double z=w0*z0+w1*z1+w2*z2; int id=yy*R+xx; if(z<rm.d[id]){rm.d[id]=z; rm.n[id]=normal; rm.m[id]=1;} }}
}
static Render renderMesh(const vector<Vec>&V,const vector<Face>&F,int view,int R){
    Render rm; rm.R=R; rm.d.assign(R*R,255.0); rm.n.assign(R*R,{0,0,0}); rm.m.assign(R*R,0);
    for(const Face&f:F){ Vec cr=crossv(V[f.b]-V[f.a],V[f.c]-V[f.a]); double len=normv(cr); if(len<=0) continue; rastTri(rm,V[f.a],V[f.b],V[f.c],cr/len,view); }
    return rm;
}
static double ssimChannel(const Render&A,const Render&B,const vector<unsigned char>&mask,int ch){
    int R=A.R, rad=3; double C1=6.5025, C2=58.5225, sum=0; int cnt=0;
    auto val=[&](const Render&M,int id)->double{ if(ch==0) return M.d[id]; const Vec&n=M.n[id]; double v= ch==1?n.x:ch==2?n.y:n.z; return (v+1.0)*127.5; };
    for(int y=0;y<R;y++) for(int x=0;x<R;x++){ int id=y*R+x; if(!mask[id]) continue; double sx=0,sy=0,sxx=0,syy=0,sxy=0; int n=0; for(int dy=-rad;dy<=rad;dy++){ int yy=min(R-1,max(0,y+dy)); for(int dx=-rad;dx<=rad;dx++){ int xx=min(R-1,max(0,x+dx)); int j=yy*R+xx; double a=val(A,j), b=val(B,j); sx+=a; sy+=b; sxx+=a*a; syy+=b*b; sxy+=a*b; n++; }} double inv=1.0/n, mx=sx*inv,my=sy*inv, vx=max(0.0,sxx*inv-mx*mx), vy=max(0.0,syy*inv-my*my), cv=sxy*inv-mx*my; double num=(2*mx*my+C1)*(2*cv+C2), den=(mx*mx+my*my+C1)*(vx+vy+C2); sum+=num/den; cnt++; }
    return cnt?sum/cnt:1.0;
}
static double proxyScore(const vector<Vec>&V,const vector<Face>&F,int R=64){
    if(elapsed()>15.5) return 0; double sc=0; for(int view=0;view<6;view++){ if(elapsed()>17.2) return 0; Render A=renderMesh(origV,origF,view,R), B=renderMesh(V,F,view,R); vector<unsigned char>mask(R*R); for(int i=0;i<R*R;i++) mask[i]=A.m[i]||B.m[i]; double ns=(ssimChannel(A,B,mask,1)+ssimChannel(A,B,mask,2)+ssimChannel(A,B,mask,3))/3.0; double ds=ssimChannel(A,B,mask,0); sc+=0.5*ns+0.5*ds; } return sc/6.0; }

static bool makeExactAABB(vector<Vec>&V, vector<Face>&F){
    Vec mn=bbMin,mx=bbMax; vector<Vec> C={{mn.x,mn.y,mn.z},{mx.x,mn.y,mn.z},{mx.x,mx.y,mn.z},{mn.x,mx.y,mn.z},{mn.x,mn.y,mx.z},{mx.x,mn.y,mx.z},{mx.x,mx.y,mx.z},{mn.x,mx.y,mx.z}};
    double e2=epsH*epsH*0.98; for(const Vec&p:origV){ double best=1e300; for(auto&q:C) best=min(best,norm2(p-q)); if(best>e2) return false; }
    int t[12][3]={{0,2,1},{0,3,2},{4,5,6},{4,6,7},{0,1,5},{0,5,4},{1,2,6},{1,6,5},{2,3,7},{2,7,6},{3,0,4},{3,4,7}};
    V=C; F.clear(); for(auto&r:t) F.push_back({r[0],r[1],r[2]}); orientOut(V,F); return validateMesh(V,F);
}

static void jacobiEigen(double a[3][3], Vec axes[3]){
    double v[3][3]={{1,0,0},{0,1,0},{0,0,1}};
    for(int it=0;it<60;it++){ int p=0,q=1; double m=fabs(a[0][1]); if(fabs(a[0][2])>m){p=0;q=2;m=fabs(a[0][2]);} if(fabs(a[1][2])>m){p=1;q=2;m=fabs(a[1][2]);} if(m<1e-18) break; double app=a[p][p], aqq=a[q][q], apq=a[p][q]; double phi=0.5*atan2(2*apq,aqq-app); double c=cos(phi),s=sin(phi); for(int k=0;k<3;k++){ double akp=a[k][p], akq=a[k][q]; a[k][p]=akp*c-akq*s; a[k][q]=akp*s+akq*c; } for(int k=0;k<3;k++){ double apk=a[p][k], aqk=a[q][k]; a[p][k]=apk*c-aqk*s; a[q][k]=apk*s+aqk*c; } a[p][q]=a[q][p]=0; for(int k=0;k<3;k++){ double vkp=v[k][p], vkq=v[k][q]; v[k][p]=vkp*c-vkq*s; v[k][q]=vkp*s+vkq*c; }}
    int id[3]={0,1,2}; sort(id,id+3,[&](int i,int j){return a[i][i]>a[j][j];}); for(int r=0;r<3;r++) axes[r]=unit({v[0][id[r]],v[1][id[r]],v[2][id[r]]});
}
static void pcaAxes(Vec&center, Vec axes[3]){
    center={0,0,0}; for(auto&p:origV) center=center+p; center=center/(double)N; double c[3][3]={{0}}; for(auto&p:origV){ Vec q=p-center; double x[3]={q.x,q.y,q.z}; for(int i=0;i<3;i++) for(int j=0;j<3;j++) c[i][j]+=x[i]*x[j]; } for(int i=0;i<3;i++) for(int j=0;j<3;j++) c[i][j]/=max(1,N); jacobiEigen(c,axes);
}
static Vec linComb(const Vec&c,const Vec&a,const Vec&b,const Vec&d,double x,double y,double z){ return c+a*x+b*y+d*z; }

static bool makeOBBGridBox(vector<Vec>&V, vector<Face>&F){
    if(N<8 || M<12 || elapsed()>4.5) return false;
    vector<pair<double,Vec>> ns; ns.reserve(min(M,200000)); int fst=max(1,M/200000);
    for(int i=0;i<M;i+=fst){ const Face&f=origF[i]; Vec cr=crossv(origV[f.b]-origV[f.a],origV[f.c]-origV[f.a]); double ar=normv(cr); if(ar>1e-20) ns.push_back({ar,cr/ar}); }
    if(ns.size()<12) return false; sort(ns.begin(),ns.end(),[](auto&a,auto&b){return a.first>b.first;});
    vector<Vec> sums; vector<double> wt; double cth=cos(8.0*M_PI/180.0);
    for(auto &pr:ns){ Vec n=pr.second; double w=pr.first; int bi=-1; double bd=0; for(int j=0;j<(int)sums.size();j++){ Vec c=unit(sums[j]); double d=fabs(dotv(c,n)); if(d>bd){bd=d;bi=j;} } if(bi>=0&&bd>cth){ Vec c=unit(sums[bi]); if(dotv(c,n)<0) n=n*(-1); sums[bi]=sums[bi]+n*w; wt[bi]+=w; } else { sums.push_back(n*w); wt.push_back(w); if(sums.size()>6) return false; } }
    if(sums.size()!=3) return false; Vec ax[3]={unit(sums[0]),unit(sums[1]),unit(sums[2])};
    for(int i=0;i<3;i++) for(int j=i+1;j<3;j++) if(fabs(dotv(ax[i],ax[j]))>0.10) return false;
    ax[0]=unit(ax[0]); ax[1]=unit(ax[1]-ax[0]*dotv(ax[1],ax[0])); ax[2]=unit(crossv(ax[0],ax[1])); if(fabs(dotv(ax[2],unit(sums[2])))<0.90) return false; if(dotv(ax[2],sums[2])<0) ax[2]=ax[2]*(-1);
    double lo[3]={1e100,1e100,1e100}, hi[3]={-1e100,-1e100,-1e100};
    for(const Vec&p:origV) for(int k=0;k<3;k++){ double t=dotv(p,ax[k]); lo[k]=min(lo[k],t); hi[k]=max(hi[k],t); }
    double len[3]={hi[0]-lo[0],hi[1]-lo[1],hi[2]-lo[2]}; if(len[0]<1e-9*diagLen||len[1]<1e-9*diagLen||len[2]<1e-9*diagLen) return false;
    double tol=0.005*diagLen; int bad=0,tot=0,pc[6]={0}; int st=max(1,N/200000);
    for(int i=0;i<N;i+=st){ const Vec&p=origV[i]; double best=1e100; int bj=0; for(int k=0;k<3;k++){ double t=dotv(p,ax[k]); double d0=fabs(t-lo[k]), d1=fabs(t-hi[k]); if(d0<best){best=d0;bj=2*k;} if(d1<best){best=d1;bj=2*k+1;} } if(best>tol) bad++; else pc[bj]++; tot++; }
    if(tot<8 || bad*200>tot) return false; for(int i=0;i<6;i++) if(pc[i]<max(1,tot/3000)) return false;
    auto pcoord=[&](double x,double y,double z){ return ax[0]*x+ax[1]*y+ax[2]*z; };
    vector<Vec> corners; for(int a=0;a<2;a++) for(int b=0;b<2;b++) for(int c=0;c<2;c++) corners.push_back(pcoord(a?hi[0]:lo[0],b?hi[1]:lo[1],c?hi[2]:lo[2]));
    bool eight=true; double e2=epsH*epsH*0.98; for(const Vec&p:origV){ double best=1e300; for(auto&q:corners) best=min(best,norm2(p-q)); if(best>e2){eight=false;break;} }
    if(eight){ int t[12][3]={{0,1,3},{0,3,2},{4,6,7},{4,7,5},{0,4,5},{0,5,1},{2,3,7},{2,7,6},{0,2,6},{0,6,4},{1,5,7},{1,7,3}}; V=corners; F.clear(); for(auto&r:t)F.push_back({r[0],r[1],r[2]}); orientOut(V,F); return validateMesh(V,F); }
    double step=epsH*1.38; int nx=max(1,(int)ceil(len[0]/step)), ny=max(1,(int)ceil(len[1]/step)), nz=max(1,(int)ceil(len[2]/step));
    unordered_map<uint64_t,int> id; id.reserve((nx+1)*(ny+1)*2+(nx+1)*(nz+1)*2+(ny+1)*(nz+1)*2); vector<Vec> vv; vector<Face> ff;
    auto get=[&](int i,int j,int k)->int{ uint64_t key=((uint64_t)i<<42)|((uint64_t)j<<21)|(uint64_t)k; auto it=id.find(key); if(it!=id.end()) return it->second; int r=vv.size(); id[key]=r; double x=lo[0]+len[0]*i/nx,y=lo[1]+len[1]*j/ny,z=lo[2]+len[2]*k/nz; vv.push_back(pcoord(x,y,z)); return r; };
    auto quad=[&](int a,int b,int c,int d){ ff.push_back({a,b,c}); ff.push_back({a,c,d}); };
    for(int i=0;i<nx;i++) for(int j=0;j<ny;j++){ quad(get(i,j,0),get(i+1,j,0),get(i+1,j+1,0),get(i,j+1,0)); quad(get(i,j,nz),get(i,j+1,nz),get(i+1,j+1,nz),get(i+1,j,nz)); }
    for(int i=0;i<nx;i++) for(int k=0;k<nz;k++){ quad(get(i,0,k),get(i,0,k+1),get(i+1,0,k+1),get(i+1,0,k)); quad(get(i,ny,k),get(i+1,ny,k),get(i+1,ny,k+1),get(i,ny,k+1)); }
    for(int j=0;j<ny;j++) for(int k=0;k<nz;k++){ quad(get(0,j,k),get(0,j+1,k),get(0,j+1,k+1),get(0,j,k+1)); quad(get(nx,j,k),get(nx,j,k+1),get(nx,j+1,k+1),get(nx,j+1,k)); }
    if(vv.size()>=origV.size()*0.92) return false; orientOut(vv,ff); if(!validateMesh(vv,ff)) return false; if(!hausdorffKD(vv,epsH*0.995)) return false; double ps=proxyScore(vv,ff,64); if(ps<0.925) return false; V.swap(vv); F.swap(ff); return true;
}

static bool makeAABBGridBox(vector<Vec>&V, vector<Face>&F){
    if(N<30 || elapsed()>4.0) return false;
    Vec mn=bbMin,mx=bbMax; double lx=mx.x-mn.x, ly=mx.y-mn.y, lz=mx.z-mn.z;
    if(lx<=1e-12*diagLen||ly<=1e-12*diagLen||lz<=1e-12*diagLen) return false;
    double tol=0.0045*diagLen; int cnt[6]={0,0,0,0,0,0}, bad=0; double av=0;
    int st=max(1,N/200000);
    for(int i=0;i<N;i+=st){ const Vec&p=origV[i]; double d[6]={fabs(p.x-mn.x),fabs(p.x-mx.x),fabs(p.y-mn.y),fabs(p.y-mx.y),fabs(p.z-mn.z),fabs(p.z-mx.z)}; int b=0; for(int k=1;k<6;k++) if(d[k]<d[b]) b=k; if(d[b]>tol) bad++; else cnt[b]++; av+=d[b]; }
    int sm=0; for(int k=0;k<6;k++) sm+=cnt[k]; if(sm<20 || bad*100>sm+bad || av/max(1,sm+bad)>tol*0.35) return false; for(int k=0;k<6;k++) if(cnt[k]<max(1,sm/2000)) return false;
    double step=epsH*1.38; int nx=max(1,(int)ceil(lx/step)), ny=max(1,(int)ceil(ly/step)), nz=max(1,(int)ceil(lz/step));
    auto pos=[&](int i,int j,int k){ return Vec{mn.x+lx*i/nx,mn.y+ly*j/ny,mn.z+lz*k/nz}; };
    unordered_map<uint64_t,int> id; id.reserve((nx+1)*(ny+1)*2+(nx+1)*(nz+1)*2+(ny+1)*(nz+1)*2);
    vector<Vec> vv; vector<Face> ff;
    auto get=[&](int i,int j,int k)->int{ uint64_t key=((uint64_t)i<<42)|((uint64_t)j<<21)|(uint64_t)k; auto it=id.find(key); if(it!=id.end()) return it->second; int r=vv.size(); id[key]=r; vv.push_back(pos(i,j,k)); return r; };
    auto quad=[&](int a,int b,int c,int d){ ff.push_back({a,b,c}); ff.push_back({a,c,d}); };
    for(int i=0;i<nx;i++) for(int j=0;j<ny;j++){ quad(get(i,j,0),get(i+1,j,0),get(i+1,j+1,0),get(i,j+1,0)); quad(get(i,j,nz),get(i,j+1,nz),get(i+1,j+1,nz),get(i+1,j,nz)); }
    for(int i=0;i<nx;i++) for(int k=0;k<nz;k++){ quad(get(i,0,k),get(i,0,k+1),get(i+1,0,k+1),get(i+1,0,k)); quad(get(i,ny,k),get(i+1,ny,k),get(i+1,ny,k+1),get(i,ny,k+1)); }
    for(int j=0;j<ny;j++) for(int k=0;k<nz;k++){ quad(get(0,j,k),get(0,j+1,k),get(0,j+1,k+1),get(0,j,k+1)); quad(get(nx,j,k),get(nx,j,k+1),get(nx,j+1,k+1),get(nx,j+1,k)); }
    if(vv.size()>=origV.size()*0.90) return false; orientOut(vv,ff); if(!validateMesh(vv,ff)) return false; if(!hausdorffKD(vv,epsH*0.995)) return false; double ps=proxyScore(vv,ff,64); if(ps<0.925) return false; V.swap(vv); F.swap(ff); return true;
}

static bool makeEllipsoid(vector<Vec>&V, vector<Face>&F){
    if(N<200 || elapsed()>3.0) return false; Vec cen; Vec ax[3]; pcaAxes(cen,ax); double rad[3]={0,0,0}; for(auto&p:origV){ Vec q=p-cen; for(int k=0;k<3;k++) rad[k]=max(rad[k],fabs(dotv(q,ax[k]))); }
    double rmax=max(rad[0],max(rad[1],rad[2])), rmin=min(rad[0],min(rad[1],rad[2])); if(rmin<0.04*diagLen || rmax<=0) return false;
    double sum=0,sq=0,mxdev=0; int cnt=0, st=max(1,N/200000); for(int i=0;i<N;i+=st){ Vec q=origV[i]-cen; double rr=0; for(int k=0;k<3;k++){ double u=dotv(q,ax[k])/rad[k]; rr+=u*u;} double dev=fabs(sqrt(max(0.0,rr))-1.0); sum+=dev; sq+=dev*dev; mxdev=max(mxdev,dev); cnt++; }
    double mean=sum/cnt, rms=sqrt(sq/cnt); if(!(mean<0.018 && rms<0.032 && mxdev<0.105)) return false;
    double step=max(epsH*0.82, rmax*0.035); int nlon=(int)ceil(2*M_PI*rmax/step); int nlat=(int)ceil(M_PI*rmax/step); nlon=max(12,min(96,nlon)); nlat=max(8,min(48,nlat));
    vector<Vec> vv; vector<Face> ff; vv.reserve(2+(nlat-1)*nlon); vv.push_back(linComb(cen,ax[0],ax[1],ax[2],0,0,rad[2])); vv.push_back(linComb(cen,ax[0],ax[1],ax[2],0,0,-rad[2]));
    auto idx=[&](int ilat,int ilon){ ilon=(ilon%nlon+nlon)%nlon; return 2+(ilat-1)*nlon+ilon; };
    for(int il=1;il<nlat;il++){ double th=M_PI*il/nlat, stt=sin(th), ctt=cos(th); for(int j=0;j<nlon;j++){ double ph=2*M_PI*j/nlon; vv.push_back(linComb(cen,ax[0],ax[1],ax[2],rad[0]*stt*cos(ph),rad[1]*stt*sin(ph),rad[2]*ctt)); }}
    for(int j=0;j<nlon;j++) ff.push_back({0,idx(1,j),idx(1,j+1)});
    for(int il=1;il<nlat-1;il++) for(int j=0;j<nlon;j++){ int a=idx(il,j),b=idx(il+1,j),c=idx(il+1,j+1),d=idx(il,j+1); ff.push_back({a,b,c}); ff.push_back({a,c,d}); }
    for(int j=0;j<nlon;j++) ff.push_back({1,idx(nlat-1,j+1),idx(nlat-1,j)});
    if(vv.size()>=origV.size()*0.85) return false; orientOut(vv,ff); if(!validateMesh(vv,ff)) return false; if(!hausdorffKD(vv,epsH*0.985)) return false; double ps=proxyScore(vv,ff,64); if(ps<0.925) return false; V.swap(vv); F.swap(ff); return true;
}

static bool makeTorus(vector<Vec>&V, vector<Face>&F){
    if(N<500 || elapsed()>8.0) return false;
    Vec cen; Vec ax[3]; pcaAxes(cen,ax); Vec u=ax[0], v=ax[1], w=unit(crossv(u,v));
    vector<double> rho; rho.reserve(N); double rmin=1e100,rmax=-1e100,zabs=0;
    int st=max(1,N/200000);
    for(int i=0;i<N;i+=st){ Vec q=origV[i]-cen; double x=dotv(q,u), y=dotv(q,v), z=dotv(q,w); double rr=sqrt(x*x+y*y); rho.push_back(rr); rmin=min(rmin,rr); rmax=max(rmax,rr); zabs=max(zabs,fabs(z)); }
    if(rho.size()<100 || rmin<=0) return false; double R=0.5*(rmin+rmax); double minor0=max(0.5*(rmax-rmin), zabs); if(R<minor0*1.15 || minor0<0.025*diagLen) return false;
    double sumr=0; for(int i=0;i<N;i+=st){ Vec q=origV[i]-cen; double x=dotv(q,u), y=dotv(q,v), z=dotv(q,w); double rr=sqrt(x*x+y*y); sumr+=sqrt((rr-R)*(rr-R)+z*z); } double rrmean=sumr/rho.size(); double sum=0,sq=0,mxd=0; int cnt=0;
    for(int i=0;i<N;i+=st){ Vec q=origV[i]-cen; double x=dotv(q,u), y=dotv(q,v), z=dotv(q,w); double rr=sqrt(x*x+y*y); double d=fabs(sqrt((rr-R)*(rr-R)+z*z)-rrmean)/max(1e-12,rrmean); sum+=d; sq+=d*d; mxd=max(mxd,d); cnt++; }
    if(sum/cnt>0.022 || sqrt(sq/cnt)>0.035 || mxd>0.14) return false;
    double step=epsH*0.82; int nt=max(16,min(160,(int)ceil(2*M_PI*R/step))); int np=max(8,min(80,(int)ceil(2*M_PI*rrmean/step)));
    vector<Vec> vv; vector<Face> ff; vv.reserve(nt*np); auto id=[&](int i,int j){i=(i%nt+nt)%nt;j=(j%np+np)%np;return i*np+j;};
    for(int i=0;i<nt;i++){ double ph=2*M_PI*i/nt; Vec dir=u*cos(ph)+v*sin(ph); for(int j=0;j<np;j++){ double th=2*M_PI*j/np; vv.push_back(cen+dir*(R+rrmean*cos(th))+w*(rrmean*sin(th))); }}
    for(int i=0;i<nt;i++) for(int j=0;j<np;j++){ int a=id(i,j), b=id(i+1,j), c=id(i+1,j+1), d=id(i,j+1); ff.push_back({a,b,c}); ff.push_back({a,c,d}); }
    if(vv.size()>=origV.size()*0.85) return false; orientOut(vv,ff); if(!validateMesh(vv,ff)) return false; if(!hausdorffKD(vv,epsH*0.985)) return false; double ps=proxyScore(vv,ff,64); if(ps<0.92) return false; V.swap(vv); F.swap(ff); return true;
}

static bool makeCylinder(vector<Vec>&V, vector<Face>&F){
    if(N<300 || elapsed()>5.0) return false; Vec cen; Vec ax[3]; pcaAxes(cen,ax); Vec w=ax[0]; double zmin=1e100,zmax=-1e100; vector<double> zs,rs; zs.reserve(N); rs.reserve(N); for(auto&p:origV){ Vec q=p-cen; double z=dotv(q,w); Vec pr=q-w*z; double r=normv(pr); zs.push_back(z); rs.push_back(r); zmin=min(zmin,z); zmax=max(zmax,z); }
    vector<double> sr=rs; nth_element(sr.begin(),sr.begin()+sr.size()*9/10,sr.end()); double R=sr[sr.size()*9/10]; double H=zmax-zmin; if(R<0.04*diagLen || H<0.12*diagLen) return false; double tol=0.018*diagLen; int ok=0,side=0,cap=0; for(int i=0;i<N;i++){ bool s=fabs(rs[i]-R)<tol && zs[i]>zmin-tol && zs[i]<zmax+tol; bool c=(fabs(zs[i]-zmin)<tol||fabs(zs[i]-zmax)<tol) && rs[i]<=R+tol; if(s||c) ok++; if(s) side++; if(c) cap++; }
    if(ok<N*97/100 || side<N/5 || cap<N/10) return false;
    Vec u=ax[1]; Vec v=unit(crossv(w,u)); u=unit(crossv(v,w)); double step=epsH*0.78; int nt=max(12,min(96,(int)ceil(2*M_PI*R/step))); int nz=max(1,min(80,(int)ceil(H/step))); int nr=max(1,min(40,(int)ceil(R/step)));
    vector<Vec> vv; vector<Face> ff; auto sideIdx=[&](int iz,int it){it=(it%nt+nt)%nt; return iz*nt+it;}; vv.reserve((nz+1)*nt+2*(1+max(0,nr-1)*nt));
    for(int iz=0;iz<=nz;iz++){ double z=zmin+H*iz/nz; Vec base=cen+w*z; for(int it=0;it<nt;it++){ double ph=2*M_PI*it/nt; vv.push_back(base+u*(R*cos(ph))+v*(R*sin(ph))); }}
    for(int iz=0;iz<nz;iz++) for(int it=0;it<nt;it++){ int a=sideIdx(iz,it), b=sideIdx(iz+1,it), c=sideIdx(iz+1,it+1), d=sideIdx(iz,it+1); ff.push_back({a,b,c}); ff.push_back({a,c,d}); }
    auto addCap=[&](bool top){ int center=vv.size(); double z=top?zmax:zmin; Vec base=cen+w*z; vv.push_back(base); vector<vector<int>> ring(nr+1); ring[0]={center}; for(int ir=1;ir<nr;ir++){ ring[ir].resize(nt); double rr=R*ir/nr; for(int it=0;it<nt;it++){ double ph=2*M_PI*it/nt; ring[ir][it]=vv.size(); vv.push_back(base+u*(rr*cos(ph))+v*(rr*sin(ph))); }} ring[nr].resize(nt); for(int it=0;it<nt;it++) ring[nr][it]=sideIdx(top?nz:0,it); if(nr==1){ for(int it=0;it<nt;it++) ff.push_back({center,ring[1][it],ring[1][(it+1)%nt]}); } else { for(int it=0;it<nt;it++) ff.push_back({center,ring[1][it],ring[1][(it+1)%nt]}); for(int ir=1;ir<nr;ir++) for(int it=0;it<nt;it++){ int a=ring[ir][it], b=ring[ir+1][it], c=ring[ir+1][(it+1)%nt], d=ring[ir][(it+1)%nt]; ff.push_back({a,b,c}); ff.push_back({a,c,d}); }} };
    addCap(false); addCap(true); if(vv.size()>=origV.size()*0.85) return false; orientOut(vv,ff); if(!validateMesh(vv,ff)) return false; if(!hausdorffKD(vv,epsH*0.985)) return false; double ps=proxyScore(vv,ff,64); if(ps<0.92) return false; V.swap(vv); F.swap(ff); return true;
}

struct CellKey{int x,y,z; bool operator==(const CellKey&o)const{return x==o.x&&y==o.y&&z==o.z;}};
struct CellHash{size_t operator()(const CellKey&k)const{ uint64_t a=(uint32_t)k.x*11995408973635179863ull ^ (uint32_t)k.y*10150724397891781847ull ^ (uint32_t)k.z*780291637; return (size_t)(a^(a>>32)); }};
static bool makeQuant(double q, vector<Vec>&V, vector<Face>&F, bool doProxy){
    unordered_map<CellKey,int,CellHash> mp; mp.reserve(N*2); vector<int> remap(N); vector<Vec> vv; vv.reserve(N);
    for(int i=0;i<N;i++){ Vec p=origV[i]; CellKey k{(int)floor((p.x-bbMin.x)/q),(int)floor((p.y-bbMin.y)/q),(int)floor((p.z-bbMin.z)/q)}; auto it=mp.find(k); if(it==mp.end()){ int id=vv.size(); mp.emplace(k,id); vv.push_back(p); remap[i]=id; } else remap[i]=it->second; }
    if(vv.size()>=origV.size()*0.97 || vv.size()<4) return false;
    vector<Face> ff; ff.reserve(origF.size()); unordered_set<uint64_t> seen; seen.reserve(origF.size()*2);
    const uint64_t LIM=(1ull<<21)-1; if(vv.size()>LIM) return false;
    for(const Face&of:origF){ int a=remap[of.a],b=remap[of.b],c=remap[of.c]; if(a==b||a==c||b==c) continue; int s0=a,s1=b,s2=c; if(s0>s1) swap(s0,s1); if(s1>s2) swap(s1,s2); if(s0>s1) swap(s0,s1); uint64_t key=((uint64_t)s0<<42)|((uint64_t)s1<<21)|(uint64_t)s2; if(seen.insert(key).second) ff.push_back({a,b,c}); }
    if(ff.size()<4) return false; compactUsed(vv,ff); if(!validateMesh(vv,ff)) return false; if(doProxy){ double ps=proxyScore(vv,ff,64); if(ps<0.925) return false; }
    V.swap(vv); F.swap(ff); return true;
}

struct Collapser{
    vector<Vec>P; vector<Face>F; vector<unsigned char> va,fa; vector<double> er; vector<vector<int>> inc; int activeF;
    double maxErr,minArea2,normalDot,maxLen;
    void init(){P=origV; F=origF; va.assign(N,1); fa.assign(M,1); er.assign(N,0); activeF=M; buildInc(); minArea2=max(1e-32,1e-26*diagLen*diagLen*diagLen*diagLen);} 
    void buildInc(){inc.assign(P.size(),{}); vector<int>d(P.size()); for(int i=0;i<(int)F.size();i++) if(fa.empty()||fa[i]){d[F[i].a]++;d[F[i].b]++;d[F[i].c]++;} for(size_t i=0;i<P.size();i++) inc[i].reserve(d[i]); for(int i=0;i<(int)F.size();i++) if(fa.empty()||fa[i]){inc[F[i].a].push_back(i); inc[F[i].b].push_back(i); inc[F[i].c].push_back(i);} }
    bool has(const Face&f,int v){return f.a==v||f.b==v||f.c==v;}
    int otherInFace(const Face&f,int a,int b){ if(f.a!=a&&f.a!=b)return f.a; if(f.b!=a&&f.b!=b)return f.b; return f.c; }
    array<int,3> triKey(int a,int b,int c){ array<int,3> s{a,b,c}; sort(s.begin(),s.end()); return s; }
    bool findEdge(int u,int v,int by[2],int op[2]){ int cnt=0; if(u<0||v<0||u>=(int)P.size()||v>=(int)P.size()||!va[u]||!va[v]) return false; for(int fid:inc[u]) if(fa[fid]&&has(F[fid],v)){ if(cnt>=2) return false; by[cnt]=fid; op[cnt]=otherInFace(F[fid],u,v); cnt++; } return cnt==2 && op[0]!=op[1] && op[0]>=0 && op[1]>=0; }
    bool linkOK(int u,int v,int op[2]){ static vector<int> mark; static int stamp=1; if((int)mark.size()<(int)P.size()) mark.assign(P.size(),0); if(++stamp==INT_MAX){fill(mark.begin(),mark.end(),0);stamp=1;} for(int fid:inc[u]) if(fa[fid]&&has(F[fid],u)){ Face f=F[fid]; int arr[3]={f.a,f.b,f.c}; for(int x:arr) if(x!=u&&x!=v) mark[x]=stamp; } int common=0; bool o0=false,o1=false; for(int fid:inc[v]) if(fa[fid]&&has(F[fid],v)){ Face f=F[fid]; int arr[3]={f.a,f.b,f.c}; for(int x:arr) if(x!=u&&x!=v&&mark[x]==stamp){ if(x!=op[0]&&x!=op[1]) return false; if(x==op[0]) o0=true; if(x==op[1]) o1=true; common++; mark[x]=stamp+1; }} return common==2&&o0&&o1; }
    bool duplicateFace(int self,int rem,int a,int b,int c,int by0,int by1){ auto key=triKey(a,b,c); int cand[3]={a,b,c}; int best=cand[0]; for(int x:cand) if(inc[x].size()<inc[best].size()) best=x; for(int fid:inc[best]) if(fa[fid]&&fid!=self&&fid!=by0&&fid!=by1){ if(has(F[fid],rem)) continue; if(triKey(F[fid].a,F[fid].b,F[fid].c)==key) return true; } return false; }
    bool evalDir(int keep,int rem,int by[2],const Vec&np,double &ne){ ne=max(er[keep]+normv(P[keep]-np), er[rem]+normv(P[rem]-np)); if(ne>maxErr) return false; vector<int> affected; affected.reserve(inc[keep].size()+inc[rem].size()); for(int id:inc[keep]) if(fa[id]&&id!=by[0]&&id!=by[1]) affected.push_back(id); for(int id:inc[rem]) if(fa[id]&&id!=by[0]&&id!=by[1]) affected.push_back(id); sort(affected.begin(),affected.end()); affected.erase(unique(affected.begin(),affected.end()),affected.end()); for(int fid:affected){ Face old=F[fid]; Face nf=old; if(nf.a==rem) nf.a=keep; if(nf.b==rem) nf.b=keep; if(nf.c==rem) nf.c=keep; if(nf.a==nf.b||nf.a==nf.c||nf.b==nf.c) return false; auto getp=[&](int id){return id==keep?np:P[id];}; Vec oldcr=crossv(P[old.b]-P[old.a],P[old.c]-P[old.a]); Vec newcr=crossv(getp(nf.b)-getp(nf.a),getp(nf.c)-getp(nf.a)); double ao=normv(oldcr), an=normv(newcr); if(ao<=0||an<=0||an*an<minArea2) return false; if(dotv(oldcr,newcr)/(ao*an)<normalDot) return false; if(duplicateFace(fid,rem,nf.a,nf.b,nf.c,by[0],by[1])) return false; } return true; }
    bool collapseEdge(int u,int v){ int by[2],op[2]; if(!findEdge(u,v,by,op)) return false; if(!linkOK(u,v,op)) return false; double best=1e100,ne=0; int keep=-1,rem=-1; Vec np{}; Vec cand[3]={P[u],P[v],(P[u]+P[v])*0.5}; int keeps[4]={u,v,u,v}, rems[4]={v,u,v,u}; Vec poss[4]={P[u],P[v],cand[2],cand[2]}; for(int i=0;i<4;i++){ double e; if(evalDir(keeps[i],rems[i],by,poss[i],e)){ double cost=e + (i>=2?0.0001*diagLen:0); if(cost<best){best=cost;keep=keeps[i];rem=rems[i];np=poss[i];ne=e;}} } if(keep<0) return false; for(int i=0;i<2;i++) if(fa[by[i]]){fa[by[i]]=0; activeF--;} for(int fid:inc[rem]) if(fa[fid]&&has(F[fid],rem)){ if(F[fid].a==rem)F[fid].a=keep; if(F[fid].b==rem)F[fid].b=keep; if(F[fid].c==rem)F[fid].c=keep; } va[rem]=0; P[keep]=np; er[keep]=ne; vector<int> merged; merged.reserve(inc[keep].size()+inc[rem].size()); for(int fid:inc[keep]) if(fa[fid]&&has(F[fid],keep)) merged.push_back(fid); for(int fid:inc[rem]) if(fa[fid]&&has(F[fid],keep)) merged.push_back(fid); sort(merged.begin(),merged.end()); merged.erase(unique(merged.begin(),merged.end()),merged.end()); inc[keep].swap(merged); vector<int>().swap(inc[rem]); return true; }
    int activeVerts(){int c=0; for(auto x:va)c+=x; return c;}
    bool run(double seconds){ init(); maxErr=epsH*0.965; vector<pair<double,uint64_t>> edges; double nds[4]={0.995,0.985,0.965,0.92}; double lfs[4]={1.05,1.35,1.70,2.00}; for(int ph=0;ph<4 && elapsed()<seconds;ph++){ normalDot=nds[ph]; maxLen=epsH*lfs[ph]; bool any=true; int iter=0; while(any&&elapsed()<seconds&&iter<3){ any=false; edges.clear(); edges.reserve(activeF*3); for(int i=0;i<(int)F.size();i++) if(fa[i]){ int a[3]={F[i].a,F[i].b,F[i].c}; for(int e=0;e<3;e++){ int u=a[e],v=a[(e+1)%3]; if(u>v) swap(u,v); double len=normv(P[u]-P[v]); if(len<=maxLen) edges.push_back({len,((uint64_t)(uint32_t)u<<32)|(uint32_t)v}); }} sort(edges.begin(),edges.end(),[](auto&a,auto&b){return a.first<b.first;}); edges.erase(unique(edges.begin(),edges.end(),[](auto&a,auto&b){return a.second==b.second;}),edges.end()); int succ=0; for(auto &pr:edges){ if(elapsed()>=seconds) break; int u=pr.second>>32, v=(uint32_t)pr.second; if(collapseEdge(u,v)){succ++; any=true;} } if(succ==0) break; iter++; if(activeVerts()<4) break; }} return true; }
    void output(vector<Vec>&V, vector<Face>&OF){ vector<int> mp(P.size(),-1); V.clear(); OF.clear(); for(size_t i=0;i<P.size();i++) if(va[i]){mp[i]=V.size(); V.push_back(P[i]);} for(size_t i=0;i<F.size();i++) if(fa[i]){ Face f{mp[F[i].a],mp[F[i].b],mp[F[i].c]}; if(f.a>=0&&f.b>=0&&f.c>=0&&f.a!=f.b&&f.a!=f.c&&f.b!=f.c) OF.push_back(f);} }
};

static void printMesh(const vector<Vec>&V,const vector<Face>&F){
    cout.setf(ios::fmtflags(0), ios::floatfield); cout<<V.size()<<" "<<F.size()<<"\n"; cout.setf(ios::fixed); cout<<setprecision(10); for(const Vec&p:V) cout<<"v "<<p.x<<" "<<p.y<<" "<<p.z<<"\n"; for(const Face&f:F) cout<<"f "<<f.a+1<<" "<<f.b+1<<" "<<f.c+1<<"\n";
}

int main(){
    t0=chrono::steady_clock::now(); readInput(); vector<Vec> bestV=origV; vector<Face> bestF=origF; bool have=validateMesh(bestV,bestF);
    vector<Vec> V; vector<Face> F;
    if(makeExactAABB(V,F)){ bestV=V; bestF=F; printMesh(bestV,bestF); return 0; }
    if(makeOBBGridBox(V,F)){ bestV=V; bestF=F; have=true; }
    if(makeAABBGridBox(V,F)){ bestV=V; bestF=F; have=true; }
    double factors[]={0.55,0.45,0.35,0.27};
    for(double f:factors){ if(elapsed()>10.0) break; V.clear(); F.clear(); bool needProxy = (f>=0.35); if(makeQuant(epsH*f,V,F,needProxy)){ if(V.size()<bestV.size()){bestV=V; bestF=F; have=true;} break; }}
    if(elapsed()<11.0){ V.clear(); F.clear(); if(makeEllipsoid(V,F) && V.size()<bestV.size()){bestV=V; bestF=F; have=true;} }
    if(elapsed()<12.0){ V.clear(); F.clear(); if(makeTorus(V,F) && V.size()<bestV.size()){bestV=V; bestF=F; have=true;} }
    if(elapsed()<13.0){ V.clear(); F.clear(); if(makeCylinder(V,F) && V.size()<bestV.size()){bestV=V; bestF=F; have=true;} }
    if(elapsed()<18.0){ Collapser C; C.run(18.2); V.clear(); F.clear(); C.output(V,F); if(validateMesh(V,F) && V.size()<bestV.size()){bestV=V; bestF=F; have=true;} }
    if(!have){bestV=origV; bestF=origF;} printMesh(bestV,bestF); return 0;
}

Sources










Pro Extended


