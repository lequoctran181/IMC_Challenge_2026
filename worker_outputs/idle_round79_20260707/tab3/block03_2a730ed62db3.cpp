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
static inline double normv(const Vec&a){return sqrt(norm2(a));}
static inline Vec unit(Vec a){double l=normv(a); if(!(l>0)) return {0,0,0}; return a/l;}
struct Face{int v[3];};
struct DFace{int v[3]; unsigned char alive=1;};

static chrono::steady_clock::time_point T0;
static double elapsed(){return chrono::duration<double>(chrono::steady_clock::now()-T0).count();}
static int N,M;
static vector<Vec> Orig;
static vector<Face> OrigF;
static double DIAG=1.0, EPS=0.05;
static Vec BBmn,BBmx;
static double ORIENT=1.0;
static inline unsigned long long ekey(int a,int b){ if(a>b) swap(a,b); return (unsigned long long)(unsigned int)a<<32 | (unsigned int)b; }
static inline array<int,3> tkey(int a,int b,int c){ array<int,3> r={a,b,c}; sort(r.begin(),r.end()); return r; }

struct MeshOut{ vector<Vec> p; vector<Face> f; string tag; double proxy=-1; };

struct FastIn{
    static const int SZ=1<<20; int idx=0,n=0; char buf[SZ];
    inline char gc(){ if(idx>=n){ n=(int)fread(buf,1,SZ,stdin); idx=0; if(!n) return 0;} return buf[idx++]; }
    inline void skip(){ char c; while((c=gc()) && (c==' '||c=='\n'||c=='\r'||c=='\t')); if(c) --idx; }
    long long nextLong(){ skip(); char c=gc(); int sg=1; if(c=='-'){sg=-1;c=gc();} long long x=0; while(c>='0'&&c<='9'){x=x*10+(c-'0'); c=gc();} return x*sg; }
    double nextDouble(){ skip(); char c=gc(); int sg=1; if(c=='-'){sg=-1;c=gc();} double x=0; while(c>='0'&&c<='9'){x=x*10+(c-'0'); c=gc();} if(c=='.'){ double m=.1; c=gc(); while(c>='0'&&c<='9'){x+=m*(c-'0');m*=.1;c=gc();}} if(c=='e'||c=='E'){ int es=1,ee=0; c=gc(); if(c=='-'){es=-1;c=gc();} else if(c=='+') c=gc(); while(c>='0'&&c<='9'){ee=ee*10+(c-'0');c=gc();} x*=pow(10.0,es*ee);} return sg*x; }
    char nextChar(){ skip(); return gc(); }
};

static void readInput(){
    FastIn in; N=(int)in.nextLong(); M=(int)in.nextLong();
    Orig.resize(N); OrigF.resize(M); BBmn={1e100,1e100,1e100}; BBmx={-1e100,-1e100,-1e100};
    for(int i=0;i<N;i++){ (void)in.nextChar(); Orig[i].x=in.nextDouble(); Orig[i].y=in.nextDouble(); Orig[i].z=in.nextDouble();
        BBmn.x=min(BBmn.x,Orig[i].x); BBmn.y=min(BBmn.y,Orig[i].y); BBmn.z=min(BBmn.z,Orig[i].z);
        BBmx.x=max(BBmx.x,Orig[i].x); BBmx.y=max(BBmx.y,Orig[i].y); BBmx.z=max(BBmx.z,Orig[i].z);
    }
    DIAG=normv(BBmx-BBmn); if(!(DIAG>0)) DIAG=1; EPS=0.05*DIAG*0.999999;
    for(int i=0;i<M;i++){ (void)in.nextChar(); OrigF[i].v[0]=(int)in.nextLong()-1; OrigF[i].v[1]=(int)in.nextLong()-1; OrigF[i].v[2]=(int)in.nextLong()-1; }
    double vol=0; for(auto &f:OrigF){ if(f.v[0]>=0&&f.v[1]>=0&&f.v[2]>=0&&f.v[0]<N&&f.v[1]<N&&f.v[2]<N) vol += dotv(Orig[f.v[0]], crossv(Orig[f.v[1]], Orig[f.v[2]])); }
    ORIENT = (vol>=0?1.0:-1.0);
}

struct Quad{
    double a00=0,a01=0,a02=0,a03=0,a11=0,a12=0,a13=0,a22=0,a23=0,a33=0;
    void addPlane(double a,double b,double c,double d,double w){
        a00+=w*a*a; a01+=w*a*b; a02+=w*a*c; a03+=w*a*d;
        a11+=w*b*b; a12+=w*b*c; a13+=w*b*d; a22+=w*c*c; a23+=w*c*d; a33+=w*d*d;
    }
    void add(const Quad&o){a00+=o.a00;a01+=o.a01;a02+=o.a02;a03+=o.a03;a11+=o.a11;a12+=o.a12;a13+=o.a13;a22+=o.a22;a23+=o.a23;a33+=o.a33;}
    double eval(const Vec&p)const{double x=p.x,y=p.y,z=p.z;return a00*x*x+2*a01*x*y+2*a02*x*z+2*a03*x+a11*y*y+2*a12*y*z+2*a13*y+a22*z*z+2*a23*z+a33;}
    bool optimum(Vec&out)const{
        double A[3][4]={{a00,a01,a02,-a03},{a01,a11,a12,-a13},{a02,a12,a22,-a23}};
        for(int i=0;i<3;i++){
            int piv=i; for(int r=i+1;r<3;r++) if(fabs(A[r][i])>fabs(A[piv][i])) piv=r;
            if(fabs(A[piv][i])<1e-14) return false;
            if(piv!=i) for(int c=i;c<4;c++) swap(A[piv][c],A[i][c]);
            double div=A[i][i]; for(int c=i;c<4;c++) A[i][c]/=div;
            for(int r=0;r<3;r++) if(r!=i){ double f=A[r][i]; for(int c=i;c<4;c++) A[r][c]-=f*A[i][c]; }
        }
        out={A[0][3],A[1][3],A[2][3]};
        return isfinite(out.x)&&isfinite(out.y)&&isfinite(out.z);
    }
};

struct Decimator{
    vector<Vec> P; vector<DFace> F; vector<unsigned char> alive; vector<vector<int>> inc; vector<double> rad; vector<Quad> Q; vector<int> mark,mark2; int stamp=1,stamp2=1; int aliveV=0,aliveF=0; double eps=0;
    Decimator(){ }
    void init(){
        P=Orig; F.resize(M); alive.assign(N,1); rad.assign(N,0); Q.assign(N,Quad()); aliveV=N; aliveF=M; vector<int> deg(N,0);
        for(int i=0;i<M;i++){ for(int k=0;k<3;k++){F[i].v[k]=OrigF[i].v[k]; deg[F[i].v[k]]++;} F[i].alive=1; }
        inc.assign(N,{}); for(int i=0;i<N;i++) inc[i].reserve(deg[i]+8);
        for(int i=0;i<M;i++){ inc[F[i].v[0]].push_back(i); inc[F[i].v[1]].push_back(i); inc[F[i].v[2]].push_back(i); }
        eps=EPS; mark.assign(N,0); mark2.assign(N,0); buildQuadrics();
    }
    inline bool has(int fid,int v)const{ const auto &f=F[fid]; return f.v[0]==v||f.v[1]==v||f.v[2]==v; }
    int third(int fid,int u,int v)const{ const auto&f=F[fid]; for(int k=0;k<3;k++){int x=f.v[k]; if(x!=u&&x!=v) return x;} return -1; }
    Vec faceCross(int fid)const{ const auto&f=F[fid]; return crossv(P[f.v[1]]-P[f.v[0]],P[f.v[2]]-P[f.v[0]]); }
    void buildQuadrics(){
        for(auto &q:Q) q=Quad();
        for(int i=0;i<(int)F.size();i++) if(F[i].alive){ Vec cr=faceCross(i); double ar=normv(cr); if(!(ar>1e-30*DIAG*DIAG)) continue; Vec n=cr/ar; double d=-dotv(n,P[F[i].v[0]]); double w=max(1e-12,ar); for(int k=0;k<3;k++) Q[F[i].v[k]].addPlane(n.x,n.y,n.z,d,w); }
    }
    void rebuildInc(){
        vector<int> deg(N,0); for(int i=0;i<(int)F.size();i++) if(F[i].alive){ for(int k=0;k<3;k++) deg[F[i].v[k]]++; }
        inc.assign(N,{}); for(int i=0;i<N;i++) if(alive[i]) inc[i].reserve(deg[i]+8);
        for(int i=0;i<(int)F.size();i++) if(F[i].alive){ for(int k=0;k<3;k++) inc[F[i].v[k]].push_back(i); }
    }
    bool edgeFaces(int u,int v,int ef[2],int op[2])const{
        if(u==v||u<0||v<0||u>=N||v>=N||!alive[u]||!alive[v]) return false; int c=0;
        for(int id:inc[u]) if(F[id].alive && has(id,v)){ if(c>=2) return false; ef[c]=id; op[c]=third(id,u,v); c++; }
        return c==2 && op[0]>=0 && op[1]>=0 && op[0]!=op[1];
    }
    bool linkOK(int u,int v,const int op[2]){
        if(++stamp==INT_MAX){ fill(mark.begin(),mark.end(),0); stamp=1; }
        if(++stamp2==INT_MAX){ fill(mark2.begin(),mark2.end(),0); stamp2=1; }
        for(int id:inc[u]) if(F[id].alive && has(id,u)){ for(int k=0;k<3;k++){int x=F[id].v[k]; if(x!=u&&x!=v) mark[x]=stamp; } }
        int cnt=0;
        for(int id:inc[v]) if(F[id].alive && has(id,v)){ for(int k=0;k<3;k++){int x=F[id].v[k]; if(x==u||x==v) continue; if(mark[x]==stamp){ if(x!=op[0]&&x!=op[1]) return false; if(mark2[x]!=stamp2){mark2[x]=stamp2; cnt++;} } } }
        return cnt==2 && mark2[op[0]]==stamp2 && mark2[op[1]]==stamp2;
    }
    struct Trial{bool ok=false; int keep=-1,rem=-1; Vec pos{}; double nr=0,cost=1e100;};
    bool dupExisting(const vector<int>&ds,int e0,int e1,int fid,const int nv[3]){
        int a=nv[0],b=nv[1],c=nv[2]; if(a==b||a==c||b==c) return true; array<int,3> key=tkey(a,b,c);
        int best=a; if(inc[b].size()<inc[best].size()) best=b; if(inc[c].size()<inc[best].size()) best=c;
        for(int gid:inc[best]){
            if(!F[gid].alive||gid==fid||gid==e0||gid==e1) continue;
            if(binary_search(ds.begin(),ds.end(),gid)) continue;
            if(tkey(F[gid].v[0],F[gid].v[1],F[gid].v[2])==key) return true;
        }
        return false;
    }
    Trial evalDir(int keep,int rem,const int ef[2],const Vec&pos,double minDot,double minArea,double planeTol){
        Trial tr; tr.keep=keep; tr.rem=rem; tr.pos=pos;
        double nr=max(rad[keep]+normv(P[keep]-pos),rad[rem]+normv(P[rem]-pos)); if(nr>eps) return tr;
        vector<int> ds; ds.reserve(inc[keep].size()+inc[rem].size());
        for(int id:inc[keep]) if(F[id].alive && id!=ef[0] && id!=ef[1] && has(id,keep)) ds.push_back(id);
        for(int id:inc[rem]) if(F[id].alive && id!=ef[0] && id!=ef[1] && has(id,rem)) ds.push_back(id);
        sort(ds.begin(),ds.end()); ds.erase(unique(ds.begin(),ds.end()),ds.end()); if(ds.empty()) return tr;
        set<array<int,3>> newKeys; double worstFlip=0,worstPlane=0,worstArea=0; int touched=0;
        for(int id:ds){
            int nv[3]={F[id].v[0],F[id].v[1],F[id].v[2]}; for(int k=0;k<3;k++) if(nv[k]==rem) nv[k]=keep;
            if(nv[0]==nv[1]||nv[0]==nv[2]||nv[1]==nv[2]) return tr;
            Vec oldp[3]={P[F[id].v[0]],P[F[id].v[1]],P[F[id].v[2]]};
            Vec newp[3]; for(int k=0;k<3;k++){ int vv=nv[k]; newp[k]=(vv==keep?pos:P[vv]); }
            Vec co=crossv(oldp[1]-oldp[0],oldp[2]-oldp[0]); Vec cn=crossv(newp[1]-newp[0],newp[2]-newp[0]);
            double ao=normv(co), an=normv(cn); if(!(ao>1e-24*DIAG*DIAG)||!(an>1e-24*DIAG*DIAG)) return tr;
            double nd=dotv(co,cn)/(ao*an); if(nd>1) nd=1; if(nd<-1) nd=-1; if(nd<minDot) return tr;
            if(an<ao*minArea) return tr;
            Vec no=co/ao; double pd=fabs(dotv(no,pos-oldp[0])); if(pd>planeTol) return tr;
            if(!newKeys.insert(tkey(nv[0],nv[1],nv[2])).second) return tr;
            if(dupExisting(ds,ef[0],ef[1],id,nv)) return tr;
            worstFlip=max(worstFlip,1.0-nd); worstPlane=max(worstPlane,pd); worstArea=max(worstArea,max(0.0,1.0-an/ao)); touched++;
        }
        tr.ok=true; tr.nr=nr; Quad qq=Q[keep]; qq.add(Q[rem]); double qcost=qq.eval(pos)/(DIAG*DIAG+1e-30);
        tr.cost=0.45*(nr/(eps+1e-30)) + 0.000001*qcost + 160.0*worstFlip + 0.7*worstPlane/(eps+1e-30) + 0.02*worstArea + 0.0004*touched + 0.000001*(inc[keep].size()+inc[rem].size());
        return tr;
    }
    bool collapse(int u,int v,double minDot,double minArea,double planeTol){
        int ef[2],op[2]; if(!edgeFaces(u,v,ef,op)) return false; if(!linkOK(u,v,op)) return false;
        Quad qq=Q[u]; qq.add(Q[v]); vector<Vec> cand; Vec opt; if(qq.optimum(opt)) cand.push_back(opt); cand.push_back((P[u]+P[v])*0.5); cand.push_back(P[u]); cand.push_back(P[v]); cand.push_back(P[u]*0.75+P[v]*0.25); cand.push_back(P[u]*0.25+P[v]*0.75);
        Trial best;
        for(Vec p:cand){
            if(!(isfinite(p.x)&&isfinite(p.y)&&isfinite(p.z))) continue;
            Trial a=evalDir(u,v,ef,p,minDot,minArea,planeTol); if(a.ok&&a.cost<best.cost) best=a;
            Trial b=evalDir(v,u,ef,p,minDot,minArea,planeTol); if(b.ok&&b.cost<best.cost) best=b;
        }
        if(!best.ok) return false; int keep=best.keep,rem=best.rem;
        for(int k=0;k<2;k++) if(F[ef[k]].alive){ F[ef[k]].alive=0; aliveF--; }
        for(int id:inc[rem]) if(F[id].alive && has(id,rem)){ for(int k=0;k<3;k++) if(F[id].v[k]==rem) F[id].v[k]=keep; if(F[id].v[0]==F[id].v[1]||F[id].v[0]==F[id].v[2]||F[id].v[1]==F[id].v[2]){F[id].alive=0; aliveF--;} else inc[keep].push_back(id); }
        P[keep]=best.pos; rad[keep]=best.nr; Q[keep].add(Q[rem]); alive[rem]=0; aliveV--; inc[rem].clear();
        if(inc[keep].size()>256){ vector<int> tmp; tmp.reserve(inc[keep].size()); for(int id:inc[keep]) if(F[id].alive&&has(id,keep)) tmp.push_back(id); sort(tmp.begin(),tmp.end()); tmp.erase(unique(tmp.begin(),tmp.end()),tmp.end()); inc[keep].swap(tmp); }
        return true;
    }
    struct Edge{int u,v; double c;};
    vector<Edge> collect(double maxLen){
        unordered_set<unsigned long long> seen; seen.reserve((size_t)aliveF*3+10); vector<Edge> e; e.reserve((size_t)aliveF*3/2+10); double ml2=maxLen*maxLen;
        for(int i=0;i<(int)F.size();i++) if(F[i].alive){ int a[3]={F[i].v[0],F[i].v[1],F[i].v[2]}; for(int k=0;k<3;k++){ int u=a[k],v=a[(k+1)%3]; if(!alive[u]||!alive[v]) continue; auto key=ekey(u,v); if(!seen.insert(key).second) continue; double d2=norm2(P[u]-P[v]); if(d2>ml2) continue; Vec m=(P[u]+P[v])*0.5; Quad qq=Q[u]; qq.add(Q[v]); double q=qq.eval(m)/(DIAG*DIAG+1e-30); e.push_back({u,v,d2+1e-7*q+1e-10*(inc[u].size()+inc[v].size())}); }}
        sort(e.begin(),e.end(),[](const Edge&a,const Edge&b){return a.c<b.c;}); return e;
    }
    MeshOut exportMesh(const string&tag){
        MeshOut out; out.tag=tag; vector<int> used(N,0),id(N,-1); vector<Face> ff; ff.reserve(aliveF);
        unordered_set<string> tri; tri.reserve(aliveF*2+3);
        for(int i=0;i<(int)F.size();i++) if(F[i].alive){ int a=F[i].v[0],b=F[i].v[1],c=F[i].v[2]; if(a==b||a==c||b==c) continue; Vec cr=crossv(P[b]-P[a],P[c]-P[a]); if(norm2(cr)<=1e-28*DIAG*DIAG*DIAG*DIAG) continue; auto tk=tkey(a,b,c); string s=to_string(tk[0])+"/"+to_string(tk[1])+"/"+to_string(tk[2]); if(!tri.insert(s).second) continue; Face f; f.v[0]=a;f.v[1]=b;f.v[2]=c; ff.push_back(f); used[a]=used[b]=used[c]=1; }
        out.p.reserve(aliveV); for(int i=0;i<N;i++) if(alive[i]&&used[i]){ id[i]=(int)out.p.size(); out.p.push_back(P[i]); }
        out.f.reserve(ff.size()); for(auto f:ff){ if(id[f.v[0]]<0||id[f.v[1]]<0||id[f.v[2]]<0) continue; Face g; g.v[0]=id[f.v[0]];g.v[1]=id[f.v[1]];g.v[2]=id[f.v[2]]; out.f.push_back(g); }
        return out;
    }
    vector<MeshOut> run(){
        vector<MeshOut> snaps; init();
        struct Phase{double md,ma,pt,ml;int sw;};
        vector<Phase> ph={{.9995,.35,.06,1.25,1},{.997,.25,.11,1.55,1},{.992,.16,.18,1.85,1},{.982,.09,.30,2.05,1},{.955,.04,.48,2.25,1},{.905,.015,.72,2.45,1},{.84,.004,.98,2.70,1}};
        int phaseNo=0;
        for(auto p:ph){
            if(elapsed()>13.2) break; int totalDone=0;
            for(int sw=0;sw<p.sw && elapsed()<13.6;sw++){
                auto edges=collect(eps*p.ml); int done=0,tried=0;
                for(auto &ed:edges){ if(elapsed()>14.0) break; tried++; if(collapse(ed.u,ed.v,p.md,p.ma,eps*p.pt)){ done++; totalDone++; if((done&1023)==0) rebuildInc(); } if(tried>450000 && done<8) break; }
                rebuildInc(); if(done<4) break;
            }
            phaseNo++; if(totalDone>0){ MeshOut mo=exportMesh("qem_phase_"+to_string(phaseNo)); if(!mo.p.empty()) snaps.push_back(std::move(mo)); }
        }
        return snaps;
    }
};

static bool validateMesh(const MeshOut&mo){
    if(mo.p.empty()||mo.f.empty()||mo.p.size()>(size_t)N) return false;
    unordered_map<unsigned long long,int> ec; ec.reserve(mo.f.size()*3+10); unordered_set<string> tri; tri.reserve(mo.f.size()*2+10); vector<int> use(mo.p.size(),0);
    double minA2=1e-28*DIAG*DIAG*DIAG*DIAG;
    for(auto &f:mo.f){ int a=f.v[0],b=f.v[1],c=f.v[2]; if(a<0||b<0||c<0||a>=(int)mo.p.size()||b>=(int)mo.p.size()||c>=(int)mo.p.size()) return false; if(a==b||a==c||b==c) return false; if(norm2(crossv(mo.p[b]-mo.p[a],mo.p[c]-mo.p[a]))<=minA2) return false; auto tk=tkey(a,b,c); string s=to_string(tk[0])+"/"+to_string(tk[1])+"/"+to_string(tk[2]); if(!tri.insert(s).second) return false; ec[ekey(a,b)]++; ec[ekey(b,c)]++; ec[ekey(c,a)]++; use[a]=use[b]=use[c]=1; }
    for(int u:use) if(!u) return false;
    for(auto &kv:ec) if(kv.second!=2) return false;
    return true;
}

struct GridNear{
    double h; const vector<Vec>*pts; unordered_map<long long, vector<int>> mp;
    static long long key(int x,int y,int z){ return (long long)x*73856093LL ^ (long long)y*19349663LL ^ (long long)z*83492791LL; }
    void build(const vector<Vec>&p,double hh){ h=hh; pts=&p; mp.clear(); mp.reserve(p.size()*2+7); for(int i=0;i<(int)p.size();i++){ int ix=(int)floor(p[i].x/h),iy=(int)floor(p[i].y/h),iz=(int)floor(p[i].z/h); mp[key(ix,iy,iz)].push_back(i); } }
    bool near(const Vec&p,double e2)const{ int ix=(int)floor(p.x/h),iy=(int)floor(p.y/h),iz=(int)floor(p.z/h); for(int dx=-1;dx<=1;dx++)for(int dy=-1;dy<=1;dy++)for(int dz=-1;dz<=1;dz++){ auto it=mp.find(key(ix+dx,iy+dy,iz+dz)); if(it==mp.end()) continue; for(int id:it->second) if(norm2(p-(*pts)[id])<=e2) return true; } return false; }
};
static bool coverAB(const vector<Vec>&A,const vector<Vec>&B,double eps){ if(A.empty()||B.empty()) return false; GridNear g; g.build(B,max(eps,1e-12)); double e2=eps*eps*1.000004; for(const Vec&p:A) if(!g.near(p,e2)) return false; return true; }
static bool vertexHausdorffOK(const MeshOut&mo){ return coverAB(Orig,mo.p,EPS) && coverAB(mo.p,Orig,EPS); }

struct Render{ vector<double>d; vector<Vec>n; vector<unsigned char>m; int R; };
static inline void projectPoint(const Vec&p,int view,int R,double&x,double&y,double&z){
    const double D=2.5; double f=800.0*((double)R/1024.0), c=0.5*R, sx=0,sy=0;
    if(view==0){sx=p.y;sy=p.z;z=D-p.x;} else if(view==1){sx=-p.y;sy=p.z;z=D+p.x;} else if(view==2){sx=-p.x;sy=p.z;z=D-p.y;} else if(view==3){sx=p.x;sy=p.z;z=D+p.y;} else if(view==4){sx=p.x;sy=p.y;z=D-p.z;} else {sx=-p.x;sy=p.y;z=D+p.z;}
    x=f*sx/z+c; y=f*sy/z+c;
}
static void rastTri(Render&r,const Vec&a,const Vec&b,const Vec&c,const Vec&normal,int view){
    double x0,y0,z0,x1,y1,z1,x2,y2,z2; int R=r.R; projectPoint(a,view,R,x0,y0,z0); projectPoint(b,view,R,x1,y1,z1); projectPoint(c,view,R,x2,y2,z2); if(z0<=0||z1<=0||z2<=0) return;
    int minx=max(0,(int)floor(min(x0,min(x1,x2))-.5)), maxx=min(R-1,(int)ceil(max(x0,max(x1,x2))+.5)); int miny=max(0,(int)floor(min(y0,min(y1,y2))-.5)), maxy=min(R-1,(int)ceil(max(y0,max(y1,y2))+.5)); if(minx>maxx||miny>maxy) return;
    double den=(y1-y2)*(x0-x2)+(x2-x1)*(y0-y2); if(fabs(den)<1e-12) return;
    for(int yy=miny;yy<=maxy;yy++){ double py=yy+.5; for(int xx=minx;xx<=maxx;xx++){ double px=xx+.5; double w0=((y1-y2)*(px-x2)+(x2-x1)*(py-y2))/den; double w1=((y2-y0)*(px-x2)+(x0-x2)*(py-y2))/den; double w2=1-w0-w1; if(w0<-1e-9||w1<-1e-9||w2<-1e-9) continue; double z=1.0/(w0/z0+w1/z1+w2/z2); int id=yy*R+xx; if(z<r.d[id]){r.d[id]=z; r.n[id]=normal; r.m[id]=1;} } }
}
static Render renderMesh(const vector<Vec>&p,const vector<Face>&f,int view,int R){
    Render r; r.R=R; r.d.assign(R*R,255.0); r.n.assign(R*R,Vec{0,0,0}); r.m.assign(R*R,0);
    for(const auto &fc:f){ Vec cr=crossv(p[fc.v[1]]-p[fc.v[0]],p[fc.v[2]]-p[fc.v[0]]); double l=normv(cr); if(!(l>0)) continue; rastTri(r,p[fc.v[0]],p[fc.v[1]],p[fc.v[2]],cr/l,view); }
    return r;
}
struct Integral{ int R; vector<double>s; void build(const vector<double>&a,int r){R=r; s.assign((R+1)*(R+1),0); for(int y=0;y<R;y++){ double row=0; for(int x=0;x<R;x++){ row+=a[y*R+x]; s[(y+1)*(R+1)+x+1]=s[y*(R+1)+x+1]+row; } } } double sum(int x0,int y0,int x1,int y1)const{ if(x0<0)x0=0;if(y0<0)y0=0;if(x1>=R)x1=R-1;if(y1>=R)y1=R-1; if(x0>x1||y0>y1) return 0; int W=R+1; return s[(y1+1)*W+x1+1]-s[y0*W+x1+1]-s[(y1+1)*W+x0]+s[y0*W+x0]; }};
static double ssimChannel(const vector<double>&A,const vector<double>&B,const vector<unsigned char>&fg,int R){
    vector<double>A2(R*R),B2(R*R),AB(R*R); for(int i=0;i<R*R;i++){A2[i]=A[i]*A[i];B2[i]=B[i]*B[i];AB[i]=A[i]*B[i];}
    Integral ia,ib,ia2,ib2,iab; ia.build(A,R); ib.build(B,R); ia2.build(A2,R); ib2.build(B2,R); iab.build(AB,R); const int rad=5; const double C1=6.5025,C2=58.5225; double res=0; int cnt=0;
    for(int y=0;y<R;y++) for(int x=0;x<R;x++){ int id=y*R+x; if(!fg[id]) continue; int x0=max(0,x-rad),x1=min(R-1,x+rad),y0=max(0,y-rad),y1=min(R-1,y+rad); double n=(x1-x0+1)*(y1-y0+1); double sx=ia.sum(x0,y0,x1,y1), sy=ib.sum(x0,y0,x1,y1); double ux=sx/n, uy=sy/n; double vx=max(0.0,ia2.sum(x0,y0,x1,y1)/n-ux*ux); double vy=max(0.0,ib2.sum(x0,y0,x1,y1)/n-uy*uy); double cov=iab.sum(x0,y0,x1,y1)/n-ux*uy; double num=(2*ux*uy+C1)*(2*cov+C2); double den=(ux*ux+uy*uy+C1)*(vx+vy+C2); res+=num/den; cnt++; }
    return cnt?res/cnt:1.0;
}
static double proxyScore(const MeshOut&mo,int R){
    double total=0; vector<Vec> op=Orig; vector<Face> of=OrigF;
    for(int view=0;view<6;view++){
        Render A=renderMesh(op,of,view,R), B=renderMesh(mo.p,mo.f,view,R); vector<unsigned char> fg(R*R); int any=0; for(int i=0;i<R*R;i++){fg[i]=(A.m[i]||B.m[i]); any+=fg[i];} if(!any){total+=1; continue;}
        vector<double> va(R*R),vb(R*R); double ns=0;
        for(int ch=0;ch<3;ch++){ for(int i=0;i<R*R;i++){ double ca=(ch==0?A.n[i].x:(ch==1?A.n[i].y:A.n[i].z)); double cb=(ch==0?B.n[i].x:(ch==1?B.n[i].y:B.n[i].z)); va[i]=(ca+1.0)*127.5; vb[i]=(cb+1.0)*127.5; } ns+=ssimChannel(va,vb,fg,R); }
        ns/=3.0; double ds=ssimChannel(A.d,B.d,fg,R); total += 0.5*ns+0.5*ds;
        if(elapsed()>18.0) break;
    }
    return total/6.0;
}

static void orientFace(vector<Vec>&p,Face&f,const Vec&center){ Vec cr=crossv(p[f.v[1]]-p[f.v[0]],p[f.v[2]]-p[f.v[0]]); Vec fc=(p[f.v[0]]+p[f.v[1]]+p[f.v[2]])/3.0; if(dotv(cr,fc-center)*ORIENT<0) swap(f.v[1],f.v[2]); }
static MeshOut makeEllipsoid(const array<Vec,3>&ax,int lat,int lon,string tag){
    double lo[3]={1e100,1e100,1e100},hi[3]={-1e100,-1e100,-1e100}; for(const Vec&p:Orig) for(int k=0;k<3;k++){double t=dotv(p,ax[k]); lo[k]=min(lo[k],t); hi[k]=max(hi[k],t);} Vec cen{0,0,0}; double r[3]; for(int k=0;k<3;k++){ double m=.5*(lo[k]+hi[k]); cen=cen+ax[k]*m; r[k]=.5*(hi[k]-lo[k]); }
    MeshOut mo; mo.tag=tag; if(r[0]<=0||r[1]<=0||r[2]<=0) return mo; mo.p.push_back(cen+ax[2]*r[2]);
    auto idx=[&](int i,int j){return 1+(i-1)*lon+(j%lon+lon)%lon;};
    for(int i=1;i<lat;i++){ double phi=M_PI*i/lat; double sp=sin(phi),cp=cos(phi); for(int j=0;j<lon;j++){ double th=2*M_PI*j/lon; mo.p.push_back(cen+ax[0]*(r[0]*sp*cos(th))+ax[1]*(r[1]*sp*sin(th))+ax[2]*(r[2]*cp)); } }
    int bot=mo.p.size(); mo.p.push_back(cen-ax[2]*r[2]); auto add=[&](int a,int b,int c){ Face f{{a,b,c}}; orientFace(mo.p,f,cen); mo.f.push_back(f); };
    for(int j=0;j<lon;j++) add(0,idx(1,j),idx(1,j+1));
    for(int i=1;i<lat-1;i++) for(int j=0;j<lon;j++){ int a=idx(i,j),b=idx(i,j+1),c=idx(i+1,j),d=idx(i+1,j+1); add(a,c,b); add(b,c,d); }
    for(int j=0;j<lon;j++) add(bot,idx(lat-1,j+1),idx(lat-1,j));
    return mo;
}
static void jacobiAxes(array<Vec,3>&axes){
    Vec c{0,0,0}; for(auto&p:Orig)c=c+p; c=c/(double)max(1,N); double A[3][3]={{0,0,0},{0,0,0},{0,0,0}}; for(auto&p:Orig){ Vec q=p-c; double v[3]={q.x,q.y,q.z}; for(int i=0;i<3;i++)for(int j=0;j<3;j++)A[i][j]+=v[i]*v[j]; }
    double V[3][3]={{1,0,0},{0,1,0},{0,0,1}}; for(int it=0;it<40;it++){ int p=0,q=1; double best=fabs(A[0][1]); if(fabs(A[0][2])>best){p=0;q=2;best=fabs(A[0][2]);} if(fabs(A[1][2])>best){p=1;q=2;best=fabs(A[1][2]);} if(best<1e-18) break; double app=A[p][p],aqq=A[q][q],apq=A[p][q]; double phi=.5*atan2(2*apq,aqq-app); double c0=cos(phi),s0=sin(phi); for(int k=0;k<3;k++){ double akp=A[k][p],akq=A[k][q]; A[k][p]=c0*akp-s0*akq; A[k][q]=s0*akp+c0*akq; } for(int k=0;k<3;k++){ double apk=A[p][k],aqk=A[q][k]; A[p][k]=c0*apk-s0*aqk; A[q][k]=s0*apk+c0*aqk; } for(int k=0;k<3;k++){ double vkp=V[k][p],vkq=V[k][q]; V[k][p]=c0*vkp-s0*vkq; V[k][q]=s0*vkp+c0*vkq; } }
    vector<pair<double,int>> ev={{A[0][0],0},{A[1][1],1},{A[2][2],2}}; sort(ev.rbegin(),ev.rend()); for(int i=0;i<3;i++){ int k=ev[i].second; axes[i]=unit({V[0][k],V[1][k],V[2][k]}); } axes[2]=unit(crossv(axes[0],axes[1])); axes[1]=unit(crossv(axes[2],axes[0]));
}
static MeshOut makeCapsule(const array<Vec,3>&ax,int lon,int capR,string tag){
    Vec e0=ax[0],e1=ax[1],e2=ax[2]; double lo=1e100,hi=-1e100; Vec cen{0,0,0}; for(auto&p:Orig){cen=cen+p; double t=dotv(p,e0); lo=min(lo,t); hi=max(hi,t);} cen=cen/(double)N; double ct=.5*(lo+hi), T=.5*(hi-lo); vector<double> rr; rr.reserve(N); for(auto&p:Orig){Vec q=p-e0*dotv(p,e0); double u=dotv(p,e1),v=dotv(p,e2); double t=fabs(dotv(p,e0)-ct); if(t<T*.65) rr.push_back(sqrt(u*u+v*v));} if(rr.empty()) return MeshOut(); nth_element(rr.begin(),rr.begin()+rr.size()/2,rr.end()); double R=rr[rr.size()/2]; if(!(R>DIAG*1e-6)) return MeshOut(); double H=max(0.0,T-R); MeshOut mo; mo.tag=tag; Vec C=e0*ct; auto loc=[&](double z,double r,double th){return C+e0*z+e1*(r*cos(th))+e2*(r*sin(th));};
    mo.p.push_back(loc(H+R,0,0)); vector<int> ringStart; auto addRing=[&](double z,double r){ ringStart.push_back((int)mo.p.size()); for(int j=0;j<lon;j++) mo.p.push_back(loc(z,r,2*M_PI*j/lon)); };
    for(int i=1;i<=capR;i++){ double ph=(M_PI/2)*i/capR; addRing(H+R*cos(ph),R*sin(ph)); }
    int cyl=max(1,(int)ceil(H/(R+1e-30)*capR)); for(int j=1;j<cyl;j++){ double t=(double)j/cyl; addRing(H*(1-2*t),R); }
    for(int i=1;i<capR;i++){ double ph=M_PI/2+(M_PI/2)*i/capR; addRing(-H+R*cos(ph),R*sin(ph)); }
    int bot=mo.p.size(); mo.p.push_back(loc(-H-R,0,0)); Vec center=C;
    auto ring=[&](int r,int j){return ringStart[r]+(j%lon+lon)%lon;}; auto add=[&](int a,int b,int c){Face f{{a,b,c}}; orientFace(mo.p,f,center); mo.f.push_back(f);}; int Rn=ringStart.size(); if(Rn==0) return MeshOut(); for(int j=0;j<lon;j++) add(0,ring(0,j+1),ring(0,j)); for(int r=0;r+1<Rn;r++) for(int j=0;j<lon;j++){int a=ring(r,j),b=ring(r,j+1),c=ring(r+1,j),d=ring(r+1,j+1); add(a,b,c); add(b,d,c);} for(int j=0;j<lon;j++) add(bot,ring(Rn-1,j),ring(Rn-1,j+1)); return mo;
}
static MeshOut makeTorus(const array<Vec,3>&ax,int maj,int minr,string tag){
    Vec c{0,0,0}; for(auto&p:Orig)c=c+p; c=c/(double)N; Vec e0=ax[0],e1=ax[1],e2=ax[2]; vector<double> qs,ds; qs.reserve(N); for(auto&p:Orig){Vec q=p-c; double u=dotv(q,e0),v=dotv(q,e1); qs.push_back(sqrt(u*u+v*v));} if(qs.empty()) return MeshOut(); sort(qs.begin(),qs.end()); double qlo=qs[qs.size()/20], qhi=qs[qs.size()*19/20]; double R=.5*(qlo+qhi); vector<double> rr; rr.reserve(N); for(auto&p:Orig){Vec q=p-c; double u=dotv(q,e0),v=dotv(q,e1),z=dotv(q,e2); rr.push_back(sqrt((sqrt(u*u+v*v)-R)*(sqrt(u*u+v*v)-R)+z*z));} nth_element(rr.begin(),rr.begin()+rr.size()/2,rr.end()); double r=rr[rr.size()/2]; if(!(R>r*1.3&&r>DIAG*1e-5)) return MeshOut(); MeshOut mo; mo.tag=tag; mo.p.reserve(maj*minr); for(int i=0;i<maj;i++){ double a=2*M_PI*i/maj; Vec radial=e0*cos(a)+e1*sin(a); for(int j=0;j<minr;j++){ double b=2*M_PI*j/minr; mo.p.push_back(c+radial*(R+r*cos(b))+e2*(r*sin(b))); } } auto id=[&](int i,int j){return ((i%maj+maj)%maj)*minr+((j%minr+minr)%minr);}; for(int i=0;i<maj;i++) for(int j=0;j<minr;j++){ Face f1{{id(i,j),id(i+1,j),id(i,j+1)}}; Face f2{{id(i,j+1),id(i+1,j),id(i+1,j+1)}}; orientFace(mo.p,f1,c); orientFace(mo.p,f2,c); mo.f.push_back(f1); mo.f.push_back(f2);} return mo;
}

static vector<MeshOut> primitiveCandidates(){
    vector<MeshOut> out; array<Vec,3> world={Vec{1,0,0},Vec{0,1,0},Vec{0,0,1}}, pca; jacobiAxes(pca);
    if(N>200 && elapsed()<14.5){ out.push_back(makeEllipsoid(world,12,28,"ellipsoid_world_12x28")); out.push_back(makeEllipsoid(pca,12,28,"ellipsoid_pca_12x28")); }
    if(N>600 && elapsed()<14.8){ out.push_back(makeEllipsoid(pca,16,36,"ellipsoid_pca_16x36")); out.push_back(makeCapsule(pca,28,7,"capsule_pca")); array<Vec,3> tor={pca[0],pca[1],pca[2]}; out.push_back(makeTorus(tor,32,12,"torus_pca")); }
    vector<MeshOut> r; for(auto &m:out) if(!m.p.empty()&&!m.f.empty()&&m.p.size()<=(size_t)N) r.push_back(std::move(m)); return r;
}

static MeshOut originalMesh(){ MeshOut mo; mo.p=Orig; mo.f=OrigF; mo.tag="original"; return mo; }
static void printMesh(const MeshOut&mo){
    static char obuf[1<<20]; setvbuf(stdout,obuf,_IOFBF,sizeof(obuf));
    printf("%d %d\n",(int)mo.p.size(),(int)mo.f.size());
    for(auto&p:mo.p) printf("v %.15g %.15g %.15g\n",p.x,p.y,p.z);
    for(auto&f:mo.f) printf("f %d %d %d\n",f.v[0]+1,f.v[1]+1,f.v[2]+1);
}

int main(){
    T0=chrono::steady_clock::now(); readInput();
    if(N<=8 || M<=12){ printMesh(originalMesh()); return 0; }
    vector<MeshOut> candidates;
    Decimator dec; vector<MeshOut> q=dec.run(); for(auto &m:q) candidates.push_back(std::move(m));
    if(elapsed()<14.8){ auto pc=primitiveCandidates(); for(auto &m:pc) candidates.push_back(std::move(m)); }
    MeshOut best=originalMesh(); int R = (N>120000?128:(N>30000?160:192)); double threshold = (R<160?0.895:0.902);
    sort(candidates.begin(),candidates.end(),[](const MeshOut&a,const MeshOut&b){return a.p.size()<b.p.size();});
    for(auto &m:candidates){
        if(elapsed()>18.6) break; if(m.p.empty()||m.f.empty()||m.p.size()>=best.p.size()) continue; if(!validateMesh(m)) continue; if(!vertexHausdorffOK(m)) continue; double pr=proxyScore(m,R); m.proxy=pr; double need=threshold; if(m.tag.find("ellipsoid")!=string::npos||m.tag.find("capsule")!=string::npos||m.tag.find("torus")!=string::npos) need-=0.020; if(m.p.size()*8<best.p.size()) need+=0.010; if(pr>=need){ best=std::move(m); break; }
    }
    if(!validateMesh(best)) best=originalMesh();
    printMesh(best); return 0;
}