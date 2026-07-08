#include <bits/stdc++.h>
using namespace std;

struct P{double x,y,z;};
static inline P operator+(P a,P b){return {a.x+b.x,a.y+b.y,a.z+b.z};}
static inline P operator-(P a,P b){return {a.x-b.x,a.y-b.y,a.z-b.z};}
static inline P operator*(P a,double s){return {a.x*s,a.y*s,a.z*s};}
static inline double dotp(P a,P b){return a.x*b.x+a.y*b.y+a.z*b.z;}
static inline P crossp(P a,P b){return {a.y*b.z-a.z*b.y,a.z*b.x-a.x*b.z,a.x*b.y-a.y*b.x};}
static inline double n2(P a){return dotp(a,a);} 
static inline double norm(P a){return sqrt(max(0.0,n2(a)));}
static inline double dist2p(P a,P b){return n2(a-b);} 

struct Face{int a,b,c; unsigned char alive;};
struct Q{double a00,a01,a02,a03,a11,a12,a13,a22,a23,a33;};
static inline Q qzero(){return {0,0,0,0,0,0,0,0,0,0};}
static inline void qadd(Q& A,const Q& B){A.a00+=B.a00;A.a01+=B.a01;A.a02+=B.a02;A.a03+=B.a03;A.a11+=B.a11;A.a12+=B.a12;A.a13+=B.a13;A.a22+=B.a22;A.a23+=B.a23;A.a33+=B.a33;}
static inline Q qsum(Q A,const Q& B){qadd(A,B);return A;}
static inline double qeval(const Q& q,const P& p){double x=p.x,y=p.y,z=p.z;return q.a00*x*x+2*q.a01*x*y+2*q.a02*x*z+2*q.a03*x+q.a11*y*y+2*q.a12*y*z+2*q.a13*y+q.a22*z*z+2*q.a23*z+q.a33;}
static inline void addPlane(Q& q,P a,P b,P c){
    P n=crossp(b-a,c-a); double l=norm(n); if(l==0) return; n=n*(1.0/l); double d=-dotp(n,a);
    double v[4]={n.x,n.y,n.z,d};
    q.a00+=v[0]*v[0]; q.a01+=v[0]*v[1]; q.a02+=v[0]*v[2]; q.a03+=v[0]*v[3];
    q.a11+=v[1]*v[1]; q.a12+=v[1]*v[2]; q.a13+=v[1]*v[3];
    q.a22+=v[2]*v[2]; q.a23+=v[2]*v[3]; q.a33+=v[3]*v[3];
}

static vector<P> V;
static vector<Face> F;
static vector<vector<int>> inc;
static vector<unsigned char> aliveV;
static vector<Q> quad;
static vector<double> radv, qw;
static vector<int> ver, markv;
static int ST=1,N0,M0;
static double EPS,EPS2,DIAG;
static int aliveCnt;

static inline bool fhas(const Face& f,int v){return f.a==v||f.b==v||f.c==v;}
static inline int otherInFace(const Face& f,int u,int v){ if(f.a!=u&&f.a!=v)return f.a; if(f.b!=u&&f.b!=v)return f.b; return f.c; }
static inline unsigned long long ekey(int a,int b){ if(a>b) swap(a,b); return (unsigned long long)(unsigned)a<<32 | (unsigned)b; }
static inline unsigned long long tkey(int a,int b,int c){ if(a>b)swap(a,b); if(b>c)swap(b,c); if(a>b)swap(a,b); return ((unsigned long long)a<<42)^((unsigned long long)b<<21)^(unsigned long long)c; }

static void cleanInc(int v){
    if(v<0||v>=N0) return;
    auto &L=inc[v];
    int w=0;
    if(aliveV[v]){
        for(int id:L) if(id>=0 && id<(int)F.size() && F[id].alive && fhas(F[id],v)) L[w++]=id;
    }
    L.resize(w);
    if(L.size()>1){ sort(L.begin(),L.end()); L.erase(unique(L.begin(),L.end()),L.end()); }
}
static void collectNei(int v, vector<int>& out){
    out.clear(); if(!aliveV[v]) return; cleanInc(v); ++ST; if(ST==INT_MAX){fill(markv.begin(),markv.end(),0); ST=1;}
    for(int id:inc[v]){ const Face &f=F[id]; int a[3]={f.a,f.b,f.c}; for(int k=0;k<3;k++){int x=a[k]; if(x!=v && aliveV[x] && markv[x]!=ST){markv[x]=ST; out.push_back(x);}} }
}
static int edgeOpps(int u,int v,int &o1,int &o2,int &f1,int &f2){
    o1=o2=f1=f2=-1; if(!aliveV[u]||!aliveV[v]) return 0; cleanInc(u); cleanInc(v);
    vector<int> tmp; tmp.reserve(4);
    const vector<int> &L = inc[u].size()<=inc[v].size()?inc[u]:inc[v];
    for(int id:L){ const Face &f=F[id]; if(f.alive && fhas(f,u) && fhas(f,v)) tmp.push_back(id); }
    sort(tmp.begin(),tmp.end()); tmp.erase(unique(tmp.begin(),tmp.end()),tmp.end());
    int cnt=0;
    for(int id:tmp){ int o=otherInFace(F[id],u,v); if(cnt==0){o1=o;f1=id;} else if(cnt==1){o2=o;f2=id;} cnt++; }
    return cnt;
}
static P fnormFace(const Face& f){return crossp(V[f.b]-V[f.a],V[f.c]-V[f.a]);}
static P fnormReplace(const Face& f,int rem,int keep){
    int a=f.a,b=f.b,c=f.c; if(a==rem)a=keep; if(b==rem)b=keep; if(c==rem)c=keep; return crossp(V[b]-V[a],V[c]-V[a]);
}
struct ER{double c; int keep,rem,f1,f2; bool ok;};

static ER evalCollapse(int keep,int rem){
    ER r{1e300,keep,rem,-1,-1,false}; if(!aliveV[keep]||!aliveV[rem]||keep==rem) return r;
    int o1,o2,f1,f2; int ec=edgeOpps(keep,rem,o1,o2,f1,f2); if(ec!=2 || o1<0||o2<0||o1==o2) return r;
    vector<int> Nu,Nv; collectNei(keep,Nu); ++ST; if(ST==INT_MAX){fill(markv.begin(),markv.end(),0); ST=1;}
    for(int x:Nu) markv[x]=ST;
    collectNei(rem,Nv); int com=0, has1=0, has2=0; ++ST; if(ST==INT_MAX){fill(markv.begin(),markv.end(),0); ST=1;}
    // mark neighbors of keep with previous stamp impossible after collectNei changed ST, redo cheaper:
    for(int x:Nu) markv[x]=ST;
    ++ST; if(ST==INT_MAX){fill(markv.begin(),markv.end(),0); ST=1;}
    int seenStamp=ST;
    for(int x:Nv){
        bool in=false; for(int y:Nu) if(y==x){in=true;break;} // degree is small; avoids another array-stamp collision
        if(in){ if(markv[x]==seenStamp) continue; markv[x]=seenStamp; com++; if(x==o1)has1=1; if(x==o2)has2=1; if(com>2) return r; }
    }
    if(com!=2||!has1||!has2) return r;
    double nd=sqrt(dist2p(V[keep],V[rem]));
    double nr=max(radv[keep],radv[rem]+nd); if(nr>EPS*0.992) return r;

    cleanInc(keep); cleanInc(rem);
    vector<unsigned long long> oldK; oldK.reserve(inc[keep].size());
    for(int id:inc[keep]) if(F[id].alive && !fhas(F[id],rem)) oldK.push_back(tkey(F[id].a,F[id].b,F[id].c));
    sort(oldK.begin(),oldK.end()); oldK.erase(unique(oldK.begin(),oldK.end()),oldK.end());
    vector<unsigned long long> newK; newK.reserve(inc[rem].size());
    double minDot=1.0, localArea=0.0;
    for(int id:inc[rem]){
        const Face &ff=F[id]; if(!ff.alive) continue;
        if(fhas(ff,keep)) continue;
        int a=ff.a,b=ff.b,c=ff.c; if(a==rem)a=keep; if(b==rem)b=keep; if(c==rem)c=keep;
        if(a==b||b==c||a==c) return r;
        unsigned long long tk=tkey(a,b,c);
        if(binary_search(oldK.begin(),oldK.end(),tk)) return r;
        newK.push_back(tk);
        P no=fnormFace(ff), nn=fnormReplace(ff,rem,keep); double lo=norm(no), ln=norm(nn);
        if(ln<1e-24*max(1.0,DIAG*DIAG)) return r;
        if(lo>0){ double d=dotp(no,nn)/(lo*ln); minDot=min(minDot,d); if(d<-0.02) return r; }
        localArea += ln;
    }
    sort(newK.begin(),newK.end()); for(size_t i=1;i<newK.size();++i) if(newK[i]==newK[i-1]) return r;
    Q qs=qsum(quad[keep],quad[rem]); double qe=max(0.0,qeval(qs,V[keep])); double w=max(1.0,qw[keep]+qw[rem]);
    if(sqrt(qe/w)>EPS*0.58) return r;
    double len2=dist2p(V[keep],V[rem]);
    r.c=qe + len2*(0.000001 + 0.02*max(0.0,1.0-minDot)) + localArea*1e-18; r.ok=true; r.f1=f1; r.f2=f2; return r;
}

static void doCollapse(const ER& e){
    int keep=e.keep, rem=e.rem; if(!e.ok) return; cleanInc(keep); cleanInc(rem);
    for(int id:inc[rem]){
        Face &ff=F[id]; if(!ff.alive) continue;
        if(fhas(ff,keep)){ ff.alive=0; continue; }
        if(ff.a==rem) ff.a=keep; if(ff.b==rem) ff.b=keep; if(ff.c==rem) ff.c=keep;
        if(ff.a==ff.b||ff.b==ff.c||ff.a==ff.c) ff.alive=0; else inc[keep].push_back(id);
    }
    aliveV[rem]=0; aliveCnt--; radv[keep]=max(radv[keep],radv[rem]+sqrt(dist2p(V[keep],V[rem])));
    qadd(quad[keep],quad[rem]); qw[keep]+=qw[rem]; ver[keep]++; ver[rem]++; inc[rem].clear(); cleanInc(keep);
}

struct Cand{double c; int u,v,vu,vv; bool operator<(Cand const& o)const{return c>o.c;}};
static priority_queue<Cand> pq;
static void pushEdge(int a,int b){ if(a==b||!aliveV[a]||!aliveV[b]) return; if(a>b) swap(a,b); double c=dist2p(V[a],V[b]); pq.push({c,a,b,ver[a],ver[b]}); }

static bool buildOutput(vector<P>& OV, vector<array<int,3>>& OF, bool validate){
    vector<int> mp(N0,-1); OV.clear(); OF.clear(); OV.reserve(aliveCnt);
    for(const Face& ff:F) if(ff.alive){ int a[3]={ff.a,ff.b,ff.c}; for(int k=0;k<3;k++) if(mp[a[k]]<0){mp[a[k]]=OV.size(); OV.push_back(V[a[k]]);} }
    OF.reserve(F.size());
    vector<unsigned long long> ek; if(validate) ek.reserve((size_t)F.size()*3);
    for(const Face& ff:F) if(ff.alive){ int a=mp[ff.a],b=mp[ff.b],c=mp[ff.c]; if(a<0||b<0||c<0||a==b||b==c||a==c) return false; OF.push_back({a,b,c}); if(validate){ek.push_back(ekey(a,b));ek.push_back(ekey(b,c));ek.push_back(ekey(c,a));} }
    if(OV.empty()||OF.empty()||OV.size()>(size_t)N0) return false;
    if(validate){ sort(ek.begin(),ek.end()); for(size_t i=0;i<ek.size();){size_t j=i+1; while(j<ek.size()&&ek[j]==ek[i])j++; if(j-i!=2) return false; i=j;} }
    return true;
}
static void printMesh(const vector<P>& OV,const vector<array<int,3>>& OF){
    cout.setf(ios::fixed); cout<<setprecision(10);
    cout<<OV.size()<<' '<<OF.size()<<'\n';
    for(auto&p:OV) cout<<p.x<<' '<<p.y<<' '<<p.z<<'\n';
    for(auto&t:OF) cout<<t[0]+1<<' '<<t[1]+1<<' '<<t[2]+1<<'\n';
}
static void printOriginal(){
    cout.setf(ios::fixed); cout<<setprecision(10);
    cout<<N0<<' '<<M0<<'\n'; for(auto&p:V) cout<<p.x<<' '<<p.y<<' '<<p.z<<'\n';
    for(auto&f:F) cout<<f.a+1<<' '<<f.b+1<<' '<<f.c+1<<'\n';
}

struct MeshOut{vector<P> v; vector<array<int,3>> f;};
static void addTriO(MeshOut& m,int a,int b,int c,P want){
    P n=crossp(m.v[b]-m.v[a],m.v[c]-m.v[a]); if(dotp(n,want)<0) swap(b,c); if(a!=b&&b!=c&&a!=c)m.f.push_back({a,b,c});
}
static bool validSmall(const MeshOut& m){
    if(m.v.empty()||m.f.empty()||m.v.size()>(size_t)N0) return false; vector<unsigned long long> e; e.reserve(m.f.size()*3);
    for(auto&t:m.f){ if(t[0]==t[1]||t[1]==t[2]||t[0]==t[2]) return false; e.push_back(ekey(t[0],t[1])); e.push_back(ekey(t[1],t[2])); e.push_back(ekey(t[2],t[0])); }
    sort(e.begin(),e.end()); for(size_t i=0;i<e.size();){size_t j=i+1; while(j<e.size()&&e[j]==e[i])j++; if(j-i!=2) return false; i=j;} return true;
}

static bool tryBox(P mn,P mx,MeshOut& out){
    double dx=mx.x-mn.x,dy=mx.y-mn.y,dz=mx.z-mn.z; if(dx<=0||dy<=0||dz<=0) return false;
    double md=0,ad=0; for(auto&p:V){double d=min({p.x-mn.x,mx.x-p.x,p.y-mn.y,mx.y-p.y,p.z-mn.z,mx.z-p.z}); md=max(md,d); ad+=d;} ad/=max(1,N0);
    if(!(md<=EPS*0.18 && ad<=EPS*0.035)) return false;
    double step=max(EPS*1.25,DIAG*1e-9); int nx=max(1,(int)ceil(dx/step)),ny=max(1,(int)ceil(dy/step)),nz=max(1,(int)ceil(dz/step));
    long long rough=2LL*((nx+1LL)*(ny+1)+(nx+1LL)*(nz+1)+(ny+1LL)*(nz+1)); if(rough>N0) return false;
    out=MeshOut(); unordered_map<unsigned long long,int> id; id.reserve((size_t)rough*2+10);
    auto get=[&](int i,int j,int k){ unsigned long long key=((unsigned long long)i<<42)^((unsigned long long)j<<21)^(unsigned long long)k; auto it=id.find(key); if(it!=id.end()) return it->second; P p{mn.x+dx*i/nx,mn.y+dy*j/ny,mn.z+dz*k/nz}; int r=out.v.size(); id[key]=r; out.v.push_back(p); return r; };
    auto quad=[&](int a,int b,int c,int d,P w){addTriO(out,a,b,c,w);addTriO(out,a,c,d,w);};
    for(int i=0;i<nx;i++)for(int j=0;j<ny;j++){quad(get(i,j,nz),get(i+1,j,nz),get(i+1,j+1,nz),get(i,j+1,nz),{0,0,1}); quad(get(i,j,0),get(i,j+1,0),get(i+1,j+1,0),get(i+1,j,0),{0,0,-1});}
    for(int i=0;i<nx;i++)for(int k=0;k<nz;k++){quad(get(i,ny,k),get(i+1,ny,k),get(i+1,ny,k+1),get(i,ny,k+1),{0,1,0}); quad(get(i,0,k),get(i,0,k+1),get(i+1,0,k+1),get(i+1,0,k),{0,-1,0});}
    for(int j=0;j<ny;j++)for(int k=0;k<nz;k++){quad(get(nx,j,k),get(nx,j+1,k),get(nx,j+1,k+1),get(nx,j,k+1),{1,0,0}); quad(get(0,j,k),get(0,j,k+1),get(0,j+1,k+1),get(0,j+1,k),{-1,0,0});}
    return validSmall(out);
}
static bool tryEllipsoid(P mn,P mx,MeshOut& out){
    P c{(mn.x+mx.x)/2,(mn.y+mx.y)/2,(mn.z+mx.z)/2}; double rx=(mx.x-mn.x)/2,ry=(mx.y-mn.y)/2,rz=(mx.z-mn.z)/2; if(rx<=0||ry<=0||rz<=0) return false;
    double rmax=max({rx,ry,rz}), rav=(rx+ry+rz)/3, md=0,ad=0;
    for(auto&p:V){double s=sqrt((p.x-c.x)*(p.x-c.x)/(rx*rx)+(p.y-c.y)*(p.y-c.y)/(ry*ry)+(p.z-c.z)*(p.z-c.z)/(rz*rz)); double d=fabs(s-1)*rav; md=max(md,d); ad+=d;} ad/=max(1,N0);
    if(!(md<=EPS*0.12 && ad<=EPS*0.028)) return false;
    double step=max(EPS*1.20,DIAG*1e-9); int nl=max(8,(int)ceil(M_PI*rmax/step)), no=max(12,(int)ceil(2*M_PI*rmax/step)); nl=min(nl,180); no=min(no,360);
    long long vv=2LL+(long long)(nl-1)*no; if(vv>N0) return false; out=MeshOut(); out.v.reserve(vv);
    int top=out.v.size(); out.v.push_back({c.x,c.y,c.z+rz}); vector<vector<int>> ring(nl-1,vector<int>(no));
    for(int i=1;i<nl;i++){double th=M_PI*i/nl,st=sin(th),ct=cos(th); for(int j=0;j<no;j++){double ph=2*M_PI*j/no; ring[i-1][j]=out.v.size(); out.v.push_back({c.x+rx*st*cos(ph),c.y+ry*st*sin(ph),c.z+rz*ct});}}
    int bot=out.v.size(); out.v.push_back({c.x,c.y,c.z-rz});
    auto want=[&](int a,int b,int cc){P g=(out.v[a]+out.v[b]+out.v[cc])*(1.0/3.0); return g-c;};
    for(int j=0;j<no;j++){int j2=(j+1)%no; addTriO(out,top,ring[0][j],ring[0][j2],want(top,ring[0][j],ring[0][j2]));}
    for(int i=0;i<nl-2;i++)for(int j=0;j<no;j++){int j2=(j+1)%no; int a=ring[i][j],b=ring[i+1][j],cc=ring[i+1][j2],d=ring[i][j2]; addTriO(out,a,b,cc,want(a,b,cc)); addTriO(out,a,cc,d,want(a,cc,d));}
    for(int j=0;j<no;j++){int j2=(j+1)%no; addTriO(out,bot,ring[nl-2][j2],ring[nl-2][j],want(bot,ring[nl-2][j2],ring[nl-2][j]));}
    return validSmall(out);
}
static P setAxis(int ax,double t,double u,double v){ if(ax==0) return {t,u,v}; if(ax==1) return {u,t,v}; return {u,v,t}; }
static bool tryCylinderAxis(int ax,P mn,P mx,MeshOut& out){
    double lo[3]={mn.x,mn.y,mn.z}, hi[3]={mx.x,mx.y,mx.z}; int ay=(ax+1)%3, az=(ax+2)%3; double H=hi[ax]-lo[ax], r1=(hi[ay]-lo[ay])/2, r2=(hi[az]-lo[az])/2; if(H<=0||r1<=0||r2<=0||fabs(r1-r2)>EPS*0.18) return false; double R=(r1+r2)/2, cy=(lo[ay]+hi[ay])/2, cz=(lo[az]+hi[az])/2;
    double md=0,ad=0; for(auto&p:V){double a[3]={p.x,p.y,p.z}; double rr=hypot(a[ay]-cy,a[az]-cz); double se=fabs(rr-R), ce=min(a[ax]-lo[ax],hi[ax]-a[ax]); if(rr>R+EPS*.2) return false; double d=min(se,ce); md=max(md,d); ad+=d;} ad/=max(1,N0); if(!(md<=EPS*0.18&&ad<=EPS*0.04)) return false;
    double step=max(EPS*1.15,DIAG*1e-9); int nt=max(12,(int)ceil(2*M_PI*R/step)), nh=max(1,(int)ceil(H/step)), nr=max(1,(int)ceil(R/step)); nt=min(nt,360); nh=min(nh,260); nr=min(nr,160); long long vv=(long long)(nh+1)*nt+2LL*(1LL+(nr-1LL)*nt); if(vv>N0) return false;
    out=MeshOut(); vector<vector<int>> side(nh+1,vector<int>(nt));
    auto pos=[&](double t,double rr,int j){double ph=2*M_PI*j/nt; return setAxis(ax,t,cy+rr*cos(ph),cz+rr*sin(ph));};
    for(int h=0;h<=nh;h++){double t=lo[ax]+H*h/nh; for(int j=0;j<nt;j++){side[h][j]=out.v.size(); out.v.push_back(pos(t,R,j));}}
    P cen=setAxis(ax,(lo[ax]+hi[ax])/2,cy,cz);
    auto want=[&](int a,int b,int c){P g=(out.v[a]+out.v[b]+out.v[c])*(1.0/3.0); return g-cen;};
    for(int h=0;h<nh;h++)for(int j=0;j<nt;j++){int j2=(j+1)%nt; int a=side[h][j],b=side[h+1][j],c=side[h+1][j2],d=side[h][j2]; addTriO(out,a,b,c,want(a,b,c)); addTriO(out,a,c,d,want(a,c,d));}
    for(int cap=0;cap<2;cap++){
        double t=cap?hi[ax]:lo[ax]; vector<vector<int>> rings(nr+1); rings[0].push_back(out.v.size()); out.v.push_back(setAxis(ax,t,cy,cz));
        for(int k=1;k<=nr;k++){rings[k].resize(nt); if(k==nr){for(int j=0;j<nt;j++) rings[k][j]=side[cap?nh:0][j];} else {double rr=R*k/nr; for(int j=0;j<nt;j++){rings[k][j]=out.v.size(); out.v.push_back(pos(t,rr,j));}}}
        for(int j=0;j<nt;j++){int j2=(j+1)%nt; addTriO(out,rings[0][0],rings[1][cap?j:j2],rings[1][cap?j2:j],setAxis(ax,cap?1:-1,0,0));}
        for(int k=2;k<=nr;k++)for(int j=0;j<nt;j++){int j2=(j+1)%nt; int a=rings[k-1][j],b=rings[k][j],c=rings[k][j2],d=rings[k-1][j2]; addTriO(out,a,b,c,setAxis(ax,cap?1:-1,0,0)); addTriO(out,a,c,d,setAxis(ax,cap?1:-1,0,0));}
    }
    return validSmall(out);
}
static bool tryCylinder(P mn,P mx,MeshOut& out){ for(int ax=0;ax<3;ax++) if(tryCylinderAxis(ax,mn,mx,out)) return true; return false; }

static bool tryFrustumAxis(int ax,P mn,P mx,MeshOut& out){
    double lo[3]={mn.x,mn.y,mn.z}, hi[3]={mx.x,mx.y,mx.z}; int ay=(ax+1)%3, az=(ax+2)%3; double H=hi[ax]-lo[ax]; if(H<=EPS) return false; double cy=(lo[ay]+hi[ay])/2, cz=(lo[az]+hi[az])/2;
    double band=max(EPS*0.35,H*0.015), rlo=0,rhi=0; int clo=0,chi=0;
    for(auto&p:V){double a[3]={p.x,p.y,p.z}; double rr=hypot(a[ay]-cy,a[az]-cz); if(a[ax]-lo[ax]<=band){rlo=max(rlo,rr);clo++;} if(hi[ax]-a[ax]<=band){rhi=max(rhi,rr);chi++;}}
    double rmax=max(rlo,rhi); if(rmax<=EPS*0.5||clo==0||chi==0) return false;
    double md=0,ad=0; for(auto&p:V){double a[3]={p.x,p.y,p.z}; double t=(a[ax]-lo[ax])/H; if(t<-0.02||t>1.02) return false; t=min(1.0,max(0.0,t)); double rr=hypot(a[ay]-cy,a[az]-cz); double re=rlo+(rhi-rlo)*t; double se=fabs(rr-re); double ce=1e100; if(a[ax]-lo[ax]<=band && rr<=rlo+EPS*.2) ce=min(ce,a[ax]-lo[ax]); if(hi[ax]-a[ax]<=band && rr<=rhi+EPS*.2) ce=min(ce,hi[ax]-a[ax]); double d=min(se,ce); md=max(md,d); ad+=d; if(rr>rmax+EPS*.25) return false; }
    ad/=max(1,N0); if(!(md<=EPS*0.18&&ad<=EPS*0.045)) return false;
    double step=max(EPS*1.12,DIAG*1e-9); int nt=max(12,(int)ceil(2*M_PI*rmax/step)), nh=max(1,(int)ceil(sqrt(H*H+(rhi-rlo)*(rhi-rlo))/step)), nr=max(1,(int)ceil(rmax/step)); nt=min(nt,360); nh=min(nh,260); nr=min(nr,160);
    long long vv=0; for(int h=0;h<=nh;h++){double rr=rlo+(rhi-rlo)*(double)h/nh; vv += (rr<EPS*0.08?1:nt);} if(rlo>EPS*0.08) vv += 1LL+(nr-1LL)*nt; if(rhi>EPS*0.08) vv += 1LL+(nr-1LL)*nt; if(vv>N0) return false;
    out=MeshOut(); vector<vector<int>> side(nh+1);
    auto pos=[&](double t,double rr,int j){double ph=2*M_PI*j/nt; return setAxis(ax,t,cy+rr*cos(ph),cz+rr*sin(ph));};
    for(int h=0;h<=nh;h++){double t=lo[ax]+H*(double)h/nh, rr=rlo+(rhi-rlo)*(double)h/nh; if(rr<EPS*0.08){side[h].push_back(out.v.size()); out.v.push_back(setAxis(ax,t,cy,cz));} else {side[h].resize(nt); for(int j=0;j<nt;j++){side[h][j]=out.v.size(); out.v.push_back(pos(t,rr,j));}}}
    P cen=setAxis(ax,(lo[ax]+hi[ax])/2,cy,cz); auto want=[&](int a,int b,int c){P g=(out.v[a]+out.v[b]+out.v[c])*(1.0/3.0); return g-cen;};
    for(int h=0;h<nh;h++){
        bool p0=side[h].size()==1,p1=side[h+1].size()==1;
        if(p0&&p1) continue;
        for(int j=0;j<nt;j++){int j2=(j+1)%nt; if(p0){int a=side[h][0],b=side[h+1][j],c=side[h+1][j2]; addTriO(out,a,b,c,want(a,b,c));} else if(p1){int a=side[h+1][0],b=side[h][j2],c=side[h][j]; addTriO(out,a,b,c,want(a,b,c));} else {int a=side[h][j],b=side[h+1][j],c=side[h+1][j2],d=side[h][j2]; addTriO(out,a,b,c,want(a,b,c)); addTriO(out,a,c,d,want(a,c,d));}}
    }
    for(int cap=0;cap<2;cap++){
        double R=cap?rhi:rlo; if(R<EPS*0.08) continue; double t=cap?hi[ax]:lo[ax]; vector<vector<int>> rings(nr+1); rings[0].push_back(out.v.size()); out.v.push_back(setAxis(ax,t,cy,cz));
        for(int k=1;k<=nr;k++){rings[k].resize(nt); if(k==nr){for(int j=0;j<nt;j++) rings[k][j]=side[cap?nh:0][j];} else {double rr=R*k/nr; for(int j=0;j<nt;j++){rings[k][j]=out.v.size(); out.v.push_back(pos(t,rr,j));}}}
        P w=setAxis(ax,cap?1:-1,0,0); for(int j=0;j<nt;j++){int j2=(j+1)%nt; addTriO(out,rings[0][0],rings[1][cap?j:j2],rings[1][cap?j2:j],w);} for(int k=2;k<=nr;k++)for(int j=0;j<nt;j++){int j2=(j+1)%nt; int a=rings[k-1][j],b=rings[k][j],c=rings[k][j2],d=rings[k-1][j2]; addTriO(out,a,b,c,w); addTriO(out,a,c,d,w);}
    }
    return validSmall(out);
}
static bool tryFrustum(P mn,P mx,MeshOut& out){ for(int ax=0;ax<3;ax++) if(tryFrustumAxis(ax,mn,mx,out)) return true; return false; }

static P torusWantP(P g,int ax,double ca,double cy,double cz,double R){
    double a[3]={g.x,g.y,g.z}; int ay=(ax+1)%3, az=(ax+2)%3; double dy=a[ay]-cy,dz=a[az]-cz; double rr=hypot(dy,dz); if(rr<1e-30) return {0,0,0}; double m[3]; m[ax]=ca; m[ay]=cy+R*dy/rr; m[az]=cz+R*dz/rr; P mp{m[0],m[1],m[2]}; return g-mp;
}
static bool tryTorusAxis(int ax,P mn,P mx,MeshOut& out){
    double lo[3]={mn.x,mn.y,mn.z}, hi[3]={mx.x,mx.y,mx.z}; int ay=(ax+1)%3, az=(ax+2)%3; double ca=(lo[ax]+hi[ax])/2, cy=(lo[ay]+hi[ay])/2, cz=(lo[az]+hi[az])/2; double rax=(hi[ax]-lo[ax])/2; if(rax<=EPS*.2) return false;
    double rmin=1e100,rmax=0; for(auto&p:V){double a[3]={p.x,p.y,p.z}; double rr=hypot(a[ay]-cy,a[az]-cz); rmin=min(rmin,rr); rmax=max(rmax,rr);} if(rmin<=EPS*.25||rmax<=rmin) return false; double R=(rmax+rmin)/2, rrad=(rmax-rmin)/2, rr=(rax+rrad)/2; if(R<=rr+EPS*.1||fabs(rax-rrad)>EPS*.18) return false;
    double md=0,ad=0; for(auto&p:V){double a[3]={p.x,p.y,p.z}; double rho=hypot(a[ay]-cy,a[az]-cz); double d=fabs(hypot(rho-R,a[ax]-ca)-rr); md=max(md,d); ad+=d;} ad/=max(1,N0); if(!(md<=EPS*.13&&ad<=EPS*.035)) return false;
    double step=max(EPS*1.15,DIAG*1e-9); int nu=max(12,(int)ceil(2*M_PI*(R+rr)/step)), nv=max(8,(int)ceil(2*M_PI*rr/step)); nu=min(nu,420); nv=min(nv,220); long long vv=1LL*nu*nv; if(vv>N0) return false;
    out=MeshOut(); out.v.reserve(vv); vector<vector<int>> id(nu,vector<int>(nv));
    for(int i=0;i<nu;i++){double u=2*M_PI*i/nu; for(int j=0;j<nv;j++){double v=2*M_PI*j/nv; double rho=R+rr*cos(v); double arr[3]; arr[ax]=ca+rr*sin(v); arr[ay]=cy+rho*cos(u); arr[az]=cz+rho*sin(u); id[i][j]=out.v.size(); out.v.push_back({arr[0],arr[1],arr[2]});}}
    auto want=[&](int a,int b,int c){P g=(out.v[a]+out.v[b]+out.v[c])*(1.0/3.0); return torusWantP(g,ax,ca,cy,cz,R);};
    for(int i=0;i<nu;i++)for(int j=0;j<nv;j++){int i2=(i+1)%nu,j2=(j+1)%nv; int a=id[i][j],b=id[i2][j],c=id[i2][j2],d=id[i][j2]; addTriO(out,a,b,c,want(a,b,c)); addTriO(out,a,c,d,want(a,c,d));}
    return validSmall(out);
}
static bool tryTorus(P mn,P mx,MeshOut& out){ for(int ax=0;ax<3;ax++) if(tryTorusAxis(ax,mn,mx,out)) return true; return false; }

int main(){
    ios::sync_with_stdio(false); cin.tie(nullptr);
    if(!(cin>>N0>>M0)) return 0; V.resize(N0); for(auto &p:V) cin>>p.x>>p.y>>p.z; F.resize(M0); for(int i=0;i<M0;i++){cin>>F[i].a>>F[i].b>>F[i].c; --F[i].a;--F[i].b;--F[i].c; F[i].alive=1;}
    P mn=V[0], mx=V[0]; for(auto&p:V){mn.x=min(mn.x,p.x);mn.y=min(mn.y,p.y);mn.z=min(mn.z,p.z);mx.x=max(mx.x,p.x);mx.y=max(mx.y,p.y);mx.z=max(mx.z,p.z);} DIAG=norm(mx-mn); EPS=max(1e-12,0.05*DIAG); EPS2=EPS*EPS;
    MeshOut prim; if(N0>1000 && (tryBox(mn,mx,prim)||tryEllipsoid(mn,mx,prim)||tryCylinder(mn,mx,prim)||tryFrustum(mn,mx,prim)||tryTorus(mn,mx,prim))){ printMesh(prim.v,prim.f); return 0; }

    inc.assign(N0,{}); aliveV.assign(N0,1); quad.assign(N0,qzero()); qw.assign(N0,0); radv.assign(N0,0); ver.assign(N0,0); markv.assign(N0,0);
    for(int i=0;i<M0;i++){Face &ff=F[i]; if(ff.a<0||ff.a>=N0||ff.b<0||ff.b>=N0||ff.c<0||ff.c>=N0||ff.a==ff.b||ff.b==ff.c||ff.a==ff.c){printOriginal();return 0;} inc[ff.a].push_back(i); inc[ff.b].push_back(i); inc[ff.c].push_back(i); addPlane(quad[ff.a],V[ff.a],V[ff.b],V[ff.c]); addPlane(quad[ff.b],V[ff.a],V[ff.b],V[ff.c]); addPlane(quad[ff.c],V[ff.a],V[ff.b],V[ff.c]); qw[ff.a]++;qw[ff.b]++;qw[ff.c]++; }
    aliveCnt=N0;
    vector<unsigned long long> ed; ed.reserve((size_t)M0*3); for(auto&ff:F){ed.push_back(ekey(ff.a,ff.b));ed.push_back(ekey(ff.b,ff.c));ed.push_back(ekey(ff.c,ff.a));}
    sort(ed.begin(),ed.end()); ed.erase(unique(ed.begin(),ed.end()),ed.end());
    for(unsigned long long k:ed) pushEdge((int)(k>>32),(int)(k&0xffffffffu)); vector<unsigned long long>().swap(ed);
    int ops=0; const int maxOps=max(0,N0-4); while(!pq.empty() && ops<maxOps){
        Cand c=pq.top(); pq.pop(); if(!aliveV[c.u]||!aliveV[c.v]||ver[c.u]!=c.vu||ver[c.v]!=c.vv) continue;
        ER a=evalCollapse(c.u,c.v), b=evalCollapse(c.v,c.u); ER e; if(a.ok&&(!b.ok||a.c<=b.c)) e=a; else if(b.ok) e=b; else continue;
        doCollapse(e); ops++;
        vector<int> nb; collectNei(e.keep,nb); for(int x:nb) pushEdge(e.keep,x);
    }
    vector<P> OV; vector<array<int,3>> OF; if(buildOutput(OV,OF,true)) printMesh(OV,OF); else printOriginal();
    return 0;
}
