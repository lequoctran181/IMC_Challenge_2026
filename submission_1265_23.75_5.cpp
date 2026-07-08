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
static inline Vec3 unit3(const Vec3&a){double l=norm3(a); return l>1e-300 ? a*(1.0/l) : Vec3{};}
static inline uint64_t ekey(int a,int b){ if(a>b) swap(a,b); return (uint64_t)(uint32_t)a<<32 | (uint32_t)b; }

struct Face{int v[3]; unsigned char live=1;};
struct Tri{int a,b,c;};

struct FastInput{
    vector<char> b; char* p=nullptr;
    FastInput(){ b.reserve(1<<27); char buf[1<<16]; size_t n; while((n=fread(buf,1,sizeof(buf),stdin))>0) b.insert(b.end(),buf,buf+n); b.push_back(0); p=b.data(); }
    inline void ws(){ while(*p==' '||*p=='\n'||*p=='\r'||*p=='\t') ++p; }
    int nextInt(){ ws(); int s=1; if(*p=='-') s=-1,++p; int x=0; while(*p>='0'&&*p<='9') x=x*10+(*p++-'0'); return x*s; }
    long nextLong(){ ws(); return strtol(p,&p,10); }
    double nextDouble(){ ws(); char* e; double x=strtod(p,&e); p=e; return x; }
    char nextChar(){ ws(); return *p++; }
};

static int N=0,M=0;
static vector<Vec3> Orig,P;
static vector<Face> F,OrigF;
static vector<vector<int>> Adj;
static vector<unsigned char> Alive;
static vector<int> Head,Tail,Next,CSz,Ver,MarkA,MarkB;
static int aliveV=0, aliveF=0, stampA=1, stampB=1;
static Vec3 BBmin,BBmax,Center;
static double Diag=1.0, Tau=0.05, Tau2=0.0025;
static chrono::steady_clock::time_point T0;
static inline double elapsed(){ return chrono::duration<double>(chrono::steady_clock::now()-T0).count(); }

static bool readInput(){
    FastInput in; N=(int)in.nextLong(); M=(int)in.nextLong();
    if(N<=0||M<=0) return false;
    Orig.resize(N); P.resize(N); Alive.assign(N,1);
    BBmin={1e100,1e100,1e100}; BBmax={-1e100,-1e100,-1e100};
    for(int i=0;i<N;i++){
        (void)in.nextChar(); Vec3 p{in.nextDouble(),in.nextDouble(),in.nextDouble()}; Orig[i]=P[i]=p;
        BBmin.x=min(BBmin.x,p.x); BBmin.y=min(BBmin.y,p.y); BBmin.z=min(BBmin.z,p.z);
        BBmax.x=max(BBmax.x,p.x); BBmax.y=max(BBmax.y,p.y); BBmax.z=max(BBmax.z,p.z);
    }
    Center=(BBmin+BBmax)*0.5; Diag=norm3(BBmax-BBmin); if(!(Diag>0)) Diag=1.0; Tau=0.05*Diag*0.999999; Tau2=Tau*Tau;
    F.resize(M); OrigF.resize(M); vector<int> deg(N,0);
    for(int i=0;i<M;i++){
        (void)in.nextChar(); int a=in.nextInt()-1,b=in.nextInt()-1,c=in.nextInt()-1;
        if((unsigned)a>=(unsigned)N||(unsigned)b>=(unsigned)N||(unsigned)c>=(unsigned)N) return false;
        F[i].v[0]=a; F[i].v[1]=b; F[i].v[2]=c; F[i].live=1; OrigF[i]=F[i];
        deg[a]++; deg[b]++; deg[c]++;
    }
    Adj.assign(N,{}); for(int i=0;i<N;i++) Adj[i].reserve(deg[i]+8);
    for(int i=0;i<M;i++){ Adj[F[i].v[0]].push_back(i); Adj[F[i].v[1]].push_back(i); Adj[F[i].v[2]].push_back(i); }
    Head.resize(N); Tail.resize(N); Next.assign(N,-1); CSz.assign(N,1); Ver.assign(N,0); MarkA.assign(N,0); MarkB.assign(N,0);
    for(int i=0;i<N;i++) Head[i]=Tail[i]=i;
    aliveV=N; aliveF=M;
    return true;
}
static inline bool hasV(const Face&f,int x){return f.v[0]==x||f.v[1]==x||f.v[2]==x;}
static inline bool hasEdge(const Face&f,int a,int b){return hasV(f,a)&&hasV(f,b);} 
static inline int thirdV(const Face&f,int a,int b){for(int k=0;k<3;k++){int x=f.v[k]; if(x!=a&&x!=b) return x;} return -1;}
static inline Vec3 faceCross(const Face&f){return cross3(P[f.v[1]]-P[f.v[0]],P[f.v[2]]-P[f.v[0]]);} 

// exact cluster-Hausdorff check: all original vertices in src cluster must be within Tau from dst representative.
static inline bool clusterCanMoveTo(int src,const Vec3&to){
    for(int m=Head[src]; m!=-1; m=Next[m]) if(norm2(Orig[m]-to)>Tau2+1e-14) return false;
    return true;
}
static bool collectEdge(int a,int b,int ef[2],int op[2]){
    const vector<int>& small = (Adj[a].size()<Adj[b].size()?Adj[a]:Adj[b]); int c=0;
    for(int fid:small){ if(!F[fid].live) continue; const Face&f=F[fid]; if(!hasEdge(f,a,b)) continue; if(c>=2) return false; ef[c]=fid; op[c]=thirdV(f,a,b); if(op[c]<0) return false; ++c; }
    return c==2 && op[0]!=op[1];
}
static bool linkOK(int a,int b,const int op[2]){
    if(++stampA>2000000000){ fill(MarkA.begin(),MarkA.end(),0); stampA=1; }
    if(++stampB>2000000000){ fill(MarkB.begin(),MarkB.end(),0); stampB=1; }
    for(int fid:Adj[a]) if(F[fid].live && hasV(F[fid],a)){ const Face&f=F[fid]; for(int k=0;k<3;k++){ int x=f.v[k]; if(x!=a&&x!=b) MarkA[x]=stampA; } }
    int cnt=0;
    for(int fid:Adj[b]) if(F[fid].live && hasV(F[fid],b)){ const Face&f=F[fid]; for(int k=0;k<3;k++){ int x=f.v[k]; if(x==a||x==b||MarkA[x]!=stampA) continue; if(x!=op[0]&&x!=op[1]) return false; if(MarkB[x]!=stampB){ MarkB[x]=stampB; ++cnt; } } }
    return cnt==2 && MarkB[op[0]]==stampB && MarkB[op[1]]==stampB;
}
static bool sameFace(int fid,int a,int b,int c){ array<int,3>x{F[fid].v[0],F[fid].v[1],F[fid].v[2]}, y{a,b,c}; sort(x.begin(),x.end()); sort(y.begin(),y.end()); return x==y; }
static bool dupAfter(int self,int rem,int e0,int e1,int a,int b,int c){
    int probe=a; if(Adj[b].size()<Adj[probe].size()) probe=b; if(Adj[c].size()<Adj[probe].size()) probe=c;
    for(int fid:Adj[probe]){ if(!F[fid].live||fid==self||fid==e0||fid==e1) continue; if(hasV(F[fid],rem)) continue; if(sameFace(fid,a,b,c)) return true; }
    return false;
}
static bool localFacesOK(int keep,int rem,const int ef[2]){
    const double minA2=max(1e-30,1e-24*Diag*Diag);
    // No normal-angle veto here: global six-view SSIM will choose the safe snapshot.  We only prevent degeneracy and duplicate faces.
    for(int fid:Adj[rem]){
        if(!F[fid].live||!hasV(F[fid],rem)) continue; if(fid==ef[0]||fid==ef[1]) continue; if(hasV(F[fid],keep)) return false;
        Face f=F[fid]; int a=f.v[0],b=f.v[1],c=f.v[2]; if(a==rem)a=keep; if(b==rem)b=keep; if(c==rem)c=keep;
        if(a==b||a==c||b==c) return false;
        Vec3 cr=cross3(P[b]-P[a],P[c]-P[a]); if(norm2(cr)<=minA2) return false;
        if(dupAfter(fid,rem,ef[0],ef[1],a,b,c)) return false;
    }
    return true;
}
struct Node{double cost; int a,b,va,vb; bool operator<(const Node&o)const{return cost>o.cost;}};
static priority_queue<Node> PQ;
static inline double edgeCost(int a,int b){
    double l=norm2(P[a]-P[b]);
    // Prefer removing high-valence/noisy vertices later by adding a tiny degree term, but primarily shortest edges.
    return l*(1.0 + 0.00001*(Adj[a].size()+Adj[b].size()));
}
static void pushEdge(int a,int b){ if(a==b||!Alive[a]||!Alive[b]) return; if(norm2(P[a]-P[b])>4.00001*Tau2) return; PQ.push({edgeCost(a,b),a,b,Ver[a],Ver[b]}); }
static void compactAdj(int v){ if(!Alive[v]||Adj[v].size()<160) return; int good=0; for(int fid:Adj[v]) if(F[fid].live&&hasV(F[fid],v)) ++good; if((int)Adj[v].size()<good*3+64) return; vector<int> nw; nw.reserve(good+8); for(int fid:Adj[v]) if(F[fid].live&&hasV(F[fid],v)) nw.push_back(fid); Adj[v].swap(nw); }
static void mergeCluster(int src,int dst){ if(Head[src]<0) return; Next[Tail[dst]]=Head[src]; Tail[dst]=Tail[src]; CSz[dst]+=CSz[src]; Head[src]=Tail[src]=-1; CSz[src]=0; }
static void doCollapse(int keep,int rem,const int ef[2]){
    for(int k=0;k<2;k++) if(F[ef[k]].live){ F[ef[k]].live=0; --aliveF; }
    vector<int> lst=Adj[rem];
    for(int fid:lst){ if(!F[fid].live||!hasV(F[fid],rem)) continue; Face&f=F[fid]; if(hasV(f,keep)){ f.live=0; --aliveF; continue; } for(int k=0;k<3;k++) if(f.v[k]==rem) f.v[k]=keep; Adj[keep].push_back(fid); }
    Alive[rem]=0; --aliveV; ++Ver[keep]; ++Ver[rem]; mergeCluster(rem,keep); Adj[rem].clear(); compactAdj(keep);
    static vector<int> neigh; neigh.clear();
    for(int fid:Adj[keep]) if(F[fid].live&&hasV(F[fid],keep)){ const Face&f=F[fid]; for(int k=0;k<3;k++){ int x=f.v[k]; if(x!=keep&&Alive[x]) neigh.push_back(x); } }
    sort(neigh.begin(),neigh.end()); neigh.erase(unique(neigh.begin(),neigh.end()),neigh.end()); for(int x:neigh) pushEdge(keep,x);
}
static bool tryCollapse(int a,int b){
    if(a==b||!Alive[a]||!Alive[b]) return false; if(norm2(P[a]-P[b])>4.00001*Tau2) return false;
    int ef[2],op[2]; if(!collectEdge(a,b,ef,op)) return false; if(!linkOK(a,b,op)) return false;
    bool ab=clusterCanMoveTo(b,P[a]) && localFacesOK(a,b,ef);
    bool ba=clusterCanMoveTo(a,P[b]) && localFacesOK(b,a,ef);
    if(!ab&&!ba) return false;
    // Prefer keeping the representative with smaller resulting cluster list scan pressure / degree.
    int keep,rem;
    if(ab&&ba){ double ca=edgeCost(a,b)+1e-12*(Adj[a].size()+CSz[a]); double cb=edgeCost(a,b)+1e-12*(Adj[b].size()+CSz[b]); if(ca<=cb){keep=a;rem=b;} else {keep=b;rem=a;} }
    else if(ab){keep=a;rem=b;} else {keep=b;rem=a;}
    doCollapse(keep,rem,ef); return true;
}

struct MeshOut{vector<Vec3> V; vector<Tri> T;};
static bool buildSnapshot(MeshOut& out){
    vector<int> id(N,-1); out.V.clear(); out.T.clear(); out.V.reserve(aliveV); out.T.reserve(aliveF);
    for(int i=0;i<N;i++) if(Alive[i]){ id[i]=(int)out.V.size(); out.V.push_back(P[i]); }
    const double minA2=max(1e-30,1e-24*Diag*Diag);
    for(int i=0;i<M;i++) if(F[i].live){ int a=F[i].v[0],b=F[i].v[1],c=F[i].v[2]; if(a==b||a==c||b==c) continue; if(id[a]<0||id[b]<0||id[c]<0) continue; Vec3 cr=cross3(P[b]-P[a],P[c]-P[a]); if(norm2(cr)<=minA2) continue; out.T.push_back({id[a],id[b],id[c]}); }
    if(out.V.empty()||out.T.empty()||out.V.size()>(size_t)N) return false;
    vector<uint64_t> ed; ed.reserve(out.T.size()*3); vector<array<int,3>> keys; keys.reserve(out.T.size());
    for(const Tri&t:out.T){ if(t.a==t.b||t.a==t.c||t.b==t.c) return false; ed.push_back(ekey(t.a,t.b)); ed.push_back(ekey(t.b,t.c)); ed.push_back(ekey(t.c,t.a)); array<int,3>k{t.a,t.b,t.c}; sort(k.begin(),k.end()); keys.push_back(k); }
    sort(keys.begin(),keys.end()); if(adjacent_find(keys.begin(),keys.end())!=keys.end()) return false;
    sort(ed.begin(),ed.end()); for(size_t i=0;i<ed.size();){ size_t j=i+1; while(j<ed.size()&&ed[j]==ed[i]) ++j; if(j-i!=2) return false; i=j; }
    return true;
}

// Exact evaluator clone (resolution selectable).  Used only as a guard; accepted fallback is the original mesh.
struct Maps{int R=0,PIX=0; vector<float> D,NM;};
static inline void project(const Vec3&p,int view,int R,double&u,double&v,double&z){
    const double D=2.5, f=800.0*(double)R/1024.0, c=0.5*(double)R; double sx,sy;
    if(view==0){sx=p.y; sy=p.z; z=D-p.x;} else if(view==1){sx=-p.y; sy=p.z; z=D+p.x;} else if(view==2){sx=-p.x; sy=p.z; z=D-p.y;} else if(view==3){sx=p.x; sy=p.z; z=D+p.y;} else if(view==4){sx=p.x; sy=p.y; z=D-p.z;} else {sx=-p.x; sy=p.y; z=D+p.z;}
    u=f*sx/z+c; v=f*sy/z+c;
}
static void renderMesh(const vector<Vec3>&V,const vector<Tri>&T,Maps& m,int R){
    int PIX=R*R; m.R=R; m.PIX=PIX; m.D.assign((size_t)6*PIX,255.f); m.NM.assign((size_t)6*PIX*3,127.5f);
    int nv=(int)V.size(); vector<float> U(nv),VV(nv),Z(nv); vector<Vec3> FN(T.size());
    for(size_t i=0;i<T.size();i++){ const Tri&t=T[i]; FN[i]=unit3(cross3(V[t.b]-V[t.a],V[t.c]-V[t.a])); }
    for(int view=0;view<6;view++){
        for(int i=0;i<nv;i++){ double u,v,z; project(V[i],view,R,u,v,z); U[i]=(float)u; VV[i]=(float)v; Z[i]=(float)z; }
        float*zbuf=m.D.data()+(size_t)view*PIX; float*nbuf=m.NM.data()+(size_t)view*PIX*3;
        for(size_t ti=0;ti<T.size();ti++){
            const Tri&t=T[ti]; int ia=t.a,ib=t.b,ic=t.c; float x0=U[ia],y0=VV[ia],z0=Z[ia],x1=U[ib],y1=VV[ib],z1=Z[ib],x2=U[ic],y2=VV[ic],z2=Z[ic]; if(!(z0>0&&z1>0&&z2>0)) continue;
            int xmin=max(0,(int)floor(min(x0,min(x1,x2))-.5f)); int xmax=min(R-1,(int)ceil(max(x0,max(x1,x2))+.5f)); int ymin=max(0,(int)floor(min(y0,min(y1,y2))-.5f)); int ymax=min(R-1,(int)ceil(max(y0,max(y1,y2))+.5f)); if(xmin>xmax||ymin>ymax) continue;
            float den=(y1-y2)*(x0-x2)+(x2-x1)*(y0-y2); if(fabs(den)<1e-20f) continue; float inv=1.f/den; Vec3 n=FN[ti]; float nr=(float)((n.x+1.)*127.5),ng=(float)((n.y+1.)*127.5),nb=(float)((n.z+1.)*127.5);
            for(int yy=ymin;yy<=ymax;yy++){ float py=(float)yy+.5f; int row=yy*R; for(int xx=xmin;xx<=xmax;xx++){ float px=(float)xx+.5f; float w0=((y1-y2)*(px-x2)+(x2-x1)*(py-y2))*inv; float w1=((y2-y0)*(px-x2)+(x0-x2)*(py-y2))*inv; float w2=1.f-w0-w1; if(w0<-1e-6f||w1<-1e-6f||w2<-1e-6f) continue; float dep=1.f/(w0/z0+w1/z1+w2/z2); int id=row+xx; if(dep<zbuf[id]){ zbuf[id]=dep; float*q=nbuf+id*3; q[0]=nr; q[1]=ng; q[2]=nb; } } }
        }
    }
}
static inline double rectSum(const vector<double>&I,int W,int x0,int y0,int x1,int y1){ return I[(size_t)y1*W+x1]-I[(size_t)y0*W+x1]-I[(size_t)y1*W+x0]+I[(size_t)y0*W+x0]; }
static double ssimChannel(const float*X,int sx,const float*Y,int sy,const float*dX,const float*dY,int R,vector<double>&IX,vector<double>&IY,vector<double>&IX2,vector<double>&IY2,vector<double>&IXY){
    int W=R+1; size_t SZ=(size_t)W*W; fill(IX.begin(),IX.begin()+SZ,0.0); fill(IY.begin(),IY.begin()+SZ,0.0); fill(IX2.begin(),IX2.begin()+SZ,0.0); fill(IY2.begin(),IY2.begin()+SZ,0.0); fill(IXY.begin(),IXY.begin()+SZ,0.0);
    for(int y=1;y<=R;y++){ double sxv=0,syv=0,sx2=0,sy2=0,sxy=0; int row=(y-1)*R; for(int x=1;x<=R;x++){ int p=row+x-1; double xv=X[(size_t)p*sx], yv=Y[(size_t)p*sy]; sxv+=xv; syv+=yv; sx2+=xv*xv; sy2+=yv*yv; sxy+=xv*yv; size_t id=(size_t)y*W+x, up=(size_t)(y-1)*W+x; IX[id]=IX[up]+sxv; IY[id]=IY[up]+syv; IX2[id]=IX2[up]+sx2; IY2[id]=IY2[up]+sy2; IXY[id]=IXY[up]+sxy; } }
    const int rad=5, area=121; const double c1=6.5025, c2=58.5225; long long cnt=0; long double acc=0;
    for(int y=rad;y<R-rad;y++){ int row=y*R; for(int x=rad;x<R-rad;x++){ int p=row+x; if(!(dX[p]<254.f||dY[p]<254.f)) continue; int x0=x-rad,y0=y-rad,x1=x+rad+1,y1=y+rad+1; double ax=rectSum(IX,W,x0,y0,x1,y1), ay=rectSum(IY,W,x0,y0,x1,y1); double ax2=rectSum(IX2,W,x0,y0,x1,y1), ay2=rectSum(IY2,W,x0,y0,x1,y1), axy=rectSum(IXY,W,x0,y0,x1,y1); double mx=ax/area,my=ay/area; double vx=max(0.0,ax2/area-mx*mx), vy=max(0.0,ay2/area-my*my), cov=axy/area-mx*my; double den=(mx*mx+my*my+c1)*(vx+vy+c2); acc += den?((2*mx*my+c1)*(2*cov+c2)/den):1.0; ++cnt; } }
    return cnt?(double)(acc/cnt):1.0;
}
static double visualSSIM(const Maps&A,const Maps&B){
    int R=A.R, PIX=R*R, W=R+1; size_t SZ=(size_t)W*W; vector<double>IX(SZ),IY(SZ),IX2(SZ),IY2(SZ),IXY(SZ); double total=0;
    for(int view=0;view<6;view++){ const float*ad=A.D.data()+(size_t)view*PIX; const float*bd=B.D.data()+(size_t)view*PIX; double ns=0; for(int c=0;c<3;c++){ const float*an=A.NM.data()+(size_t)view*PIX*3+c; const float*bn=B.NM.data()+(size_t)view*PIX*3+c; ns+=ssimChannel(an,3,bn,3,ad,bd,R,IX,IY,IX2,IY2,IXY); } ns/=3.0; double ds=ssimChannel(ad,1,bd,1,ad,bd,R,IX,IY,IX2,IY2,IXY); total+=0.5*(ns+ds); }
    return total/6.0;
}
static MeshOut makeOriginalOut(){ MeshOut o; o.V=Orig; o.T.reserve(M); for(const auto&f:OrigF) o.T.push_back({f.v[0],f.v[1],f.v[2]}); return o; }

static bool validMeshOut(const MeshOut&o){
    if(o.V.empty()||o.T.empty()||o.V.size()>(size_t)N) return false;
    const double minA2=max(1e-30,1e-24*Diag*Diag);
    vector<uint64_t> ed; ed.reserve(o.T.size()*3); vector<array<int,3>> keys; keys.reserve(o.T.size());
    for(const Tri&t:o.T){
        if(t.a<0||t.b<0||t.c<0||t.a>=(int)o.V.size()||t.b>=(int)o.V.size()||t.c>=(int)o.V.size()) return false;
        if(t.a==t.b||t.a==t.c||t.b==t.c) return false;
        if(norm2(cross3(o.V[t.b]-o.V[t.a],o.V[t.c]-o.V[t.a]))<=minA2) return false;
        ed.push_back(ekey(t.a,t.b)); ed.push_back(ekey(t.b,t.c)); ed.push_back(ekey(t.c,t.a));
        array<int,3> k{t.a,t.b,t.c}; sort(k.begin(),k.end()); keys.push_back(k);
    }
    sort(keys.begin(),keys.end()); if(adjacent_find(keys.begin(),keys.end())!=keys.end()) return false;
    sort(ed.begin(),ed.end()); for(size_t i=0;i<ed.size();){ size_t j=i+1; while(j<ed.size()&&ed[j]==ed[i]) ++j; if(j-i!=2) return false; i=j; }
    return true;
}
static double orientSignOriginal(){
    long double s=0; int step=max(1,M/200000), cnt=0;
    for(int i=0;i<M;i+=step){ const Face&f=OrigF[i]; Vec3 cr=cross3(Orig[f.v[1]]-Orig[f.v[0]],Orig[f.v[2]]-Orig[f.v[0]]); Vec3 c=(Orig[f.v[0]]+Orig[f.v[1]]+Orig[f.v[2]])/3.0; s+=dot3(cr,c-Center); if(++cnt>200000) break; }
    return s>=0?1.0:-1.0;
}
static void addFaceOrient(MeshOut&o,int a,int b,int c,double sg){
    Tri t{a,b,c}; Vec3 cr=cross3(o.V[b]-o.V[a],o.V[c]-o.V[a]); Vec3 cc=(o.V[a]+o.V[b]+o.V[c])/3.0; if(dot3(cr,cc-Center)*sg<0) swap(t.b,t.c); o.T.push_back(t);
}
static bool sameTriOrig(int fid,int a,int b,int c){
    if(fid<0||fid>=M) return false; array<int,3>x{OrigF[fid].v[0],OrigF[fid].v[1],OrigF[fid].v[2]},y{a,b,c}; sort(x.begin(),x.end()); sort(y.begin(),y.end()); return x==y;
}
static int nearestCyclic(const vector<int>&s,int n,int x,int delta){
    int m=(int)s.size(); if(!m) return -1; int pos=(int)((long long)x*m/max(1,n)); int best=s[(pos%m+m)%m], bd=INT_MAX;
    for(int dd=-delta;dd<=delta;dd++){ int pp=(pos+dd)%m; if(pp<0) pp+=m; int y=s[pp]; int d=abs(y-x); d=min(d,n-d); if(d<bd){bd=d;best=y;} }
    return best;
}
static bool detectLatLongGrid(int&R,int&Vv){
    if(N<10||M!=2*(N-2)) return false;
    for(int v=8; v<=4096; ++v){ if((N-2)%v) continue; int r=(N-2)/v; if(r<3) continue; bool ok=true; int st=max(1,v/64);
        for(int j=0;j<v&&ok;j+=st){ int a=1+j,b=1+(j+1)%v; if(!sameTriOrig(j,0,b,a)) ok=false; int off=v+2*(r-1)*v+j; int c=1+(r-1)*v+j,d=1+(r-1)*v+(j+1)%v; if(ok&&!sameTriOrig(off,N-1,c,d)) ok=false; }
        int span=max(1,(r-1)*v/256);
        for(int q=0;q<(r-1)*v&&ok;q+=span){ int rr=q/v,j=q%v; int a=1+rr*v+j,b=1+rr*v+(j+1)%v,c=1+(rr+1)*v+j,d=1+(rr+1)*v+(j+1)%v; int fid=v+2*q; if(!sameTriOrig(fid,a,b,c)||!sameTriOrig(fid+1,b,d,c)) ok=false; }
        if(ok){R=r; Vv=v; return true;}
    }
    return false;
}
static bool coverageLatLong(int R,int Vv,const vector<int>&rs,const vector<int>&ls,const MeshOut&o){
    for(int ri=0;ri<R;ri++){
        auto it=lower_bound(rs.begin(),rs.end(),ri); int candR[2], rc=0; if(it!=rs.end()) candR[rc++]=*it; if(it!=rs.begin()) candR[rc++]=*prev(it);
        for(int j=0;j<Vv;j++){ int lj=nearestCyclic(ls,Vv,j,2); int lid=(int)(find(ls.begin(),ls.end(),lj)-ls.begin()); if(lid<0||lid>=(int)ls.size()) return false; double best=1e100;
            for(int z=0;z<rc;z++){ int rr=candR[z]; int rid=(int)(find(rs.begin(),rs.end(),rr)-rs.begin()); for(int dd=-1;dd<=1;dd++){ int ll=(lid+dd+(int)ls.size())%(int)ls.size(); int id=1+rid*(int)ls.size()+ll; best=min(best,norm2(Orig[1+ri*Vv+j]-o.V[id])); } }
            if(best>Tau2*1.0001) return false;
        }
    }
    return true;
}
static bool makeLatLongGrid(int R,int Vv,int R2,int V2,MeshOut&o){
    if(R2<2||V2<8||R2>R||V2>Vv) return false; vector<int>rs,ls;
    for(int i=0;i<R2;i++){ int r=(int)llround((double)i*(R-1)/max(1,R2-1)); if(!rs.empty()&&r<=rs.back()) r=rs.back()+1; if(r>=R) r=R-1; rs.push_back(r); }
    sort(rs.begin(),rs.end()); rs.erase(unique(rs.begin(),rs.end()),rs.end()); if(rs.size()<2) return false; R2=(int)rs.size();
    for(int j=0;j<V2;j++){ int l=(int)((long long)j*Vv/V2); if(ls.empty()||l!=ls.back()) ls.push_back(l); } if(ls.size()<8) return false; V2=(int)ls.size();
    o.V.clear(); o.T.clear(); o.V.reserve(2+(size_t)R2*V2); o.V.push_back(Orig[0]); for(int r:rs) for(int l:ls) o.V.push_back(Orig[1+r*Vv+l]); int bot=(int)o.V.size(); o.V.push_back(Orig[N-1]);
    auto id=[&](int r,int j){return 1+r*V2+((j%V2+V2)%V2);}; double sg=orientSignOriginal();
    for(int j=0;j<V2;j++) addFaceOrient(o,0,id(0,j+1),id(0,j),sg);
    for(int r=0;r+1<R2;r++) for(int j=0;j<V2;j++){ int a=id(r,j),b=id(r,j+1),c=id(r+1,j),d=id(r+1,j+1); addFaceOrient(o,a,b,c,sg); addFaceOrient(o,b,d,c,sg); }
    for(int j=0;j<V2;j++) addFaceOrient(o,bot,id(R2-1,j),id(R2-1,j+1),sg);
    return coverageLatLong(R,Vv,rs,ls,o)&&validMeshOut(o);
}
static bool faceQuadMember(int fid,int a,int b,int c,int d){ return sameTriOrig(fid,a,b,c)||sameTriOrig(fid,a,b,d)||sameTriOrig(fid,a,c,d)||sameTriOrig(fid,b,c,d); }
static bool detectTorusGrid(int&U,int&Vv){
    if(N<36||M!=2*N) return false; vector<pair<int,int>> cand; for(int v=6;(long long)v*v<=N;v++) if(N%v==0){cand.push_back({N/v,v}); cand.push_back({v,N/v});}
    sort(cand.begin(),cand.end(),[](auto&a,auto&b){return abs(a.first-a.second)<abs(b.first-b.second);});
    for(auto [u,v]:cand){ if(u<6||v<6||u>20000||v>20000) continue; int st=max(1,N/1000),tot=0,okc=0; for(int cell=0; cell<N&&tot<1000; cell+=st){ int i=cell/v,j=cell%v,ni=(i+1==u?0:i+1),nj=(j+1==v?0:j+1); int a=i*v+j,b=i*v+nj,c=ni*v+j,d=ni*v+nj; int f0=2*cell,f1=f0+1; if(f1>=M) break; if(faceQuadMember(f0,a,b,c,d)&&faceQuadMember(f1,a,b,c,d)) ++okc; ++tot; } if(tot>100&&okc*100>=tot*98){U=u; Vv=v; return true;} }
    return false;
}
static bool coverageTorus(int U,int Vv,const vector<int>&us,const vector<int>&vs,const MeshOut&o){
    for(int i=0;i<U;i++) for(int j=0;j<Vv;j++){ int bu=nearestCyclic(us,U,i,2), bv=nearestCyclic(vs,Vv,j,2); int uid=(int)(find(us.begin(),us.end(),bu)-us.begin()), vid=(int)(find(vs.begin(),vs.end(),bv)-vs.begin()); if(uid<0||uid>=(int)us.size()||vid<0||vid>=(int)vs.size()) return false; double best=1e100; for(int du=-1;du<=1;du++) for(int dv=-1;dv<=1;dv++){ int uu=(uid+du+(int)us.size())%(int)us.size(); int vv=(vid+dv+(int)vs.size())%(int)vs.size(); best=min(best,norm2(Orig[i*Vv+j]-o.V[uu*(int)vs.size()+vv])); } if(best>Tau2*1.0001) return false; }
    return true;
}
static bool makeTorusGrid(int U,int Vv,int U2,int V2,MeshOut&o,bool flip=false){
    if(U2<4||V2<4||U2>U||V2>Vv||U2*V2>=N) return false; vector<int>us,vs; for(int i=0;i<U2;i++){int x=(int)((long long)i*U/U2); if(us.empty()||x!=us.back()) us.push_back(x);} for(int j=0;j<V2;j++){int x=(int)((long long)j*Vv/V2); if(vs.empty()||x!=vs.back()) vs.push_back(x);} U2=(int)us.size(); V2=(int)vs.size(); if(U2<4||V2<4) return false;
    o.V.clear(); o.T.clear(); o.V.reserve((size_t)U2*V2); o.T.reserve((size_t)U2*V2*2); for(int u:us) for(int v:vs) o.V.push_back(Orig[u*Vv+v]); auto id=[&](int i,int j){return ((i%U2+U2)%U2)*V2+((j%V2+V2)%V2);}; for(int i=0;i<U2;i++) for(int j=0;j<V2;j++){ Tri t{id(i,j),id(i,j+1),id(i+1,j)}; if(flip) swap(t.b,t.c); o.T.push_back(t); t={id(i,j+1),id(i+1,j+1),id(i+1,j)}; if(flip) swap(t.b,t.c); o.T.push_back(t); }
    return coverageTorus(U,Vv,us,vs,o)&&validMeshOut(o);
}
static bool visualPass(const MeshOut&s,int R,double need){
    static int cachedR=0; static Maps O; static MeshOut orig;
    if(cachedR!=R){ orig=makeOriginalOut(); renderMesh(orig.V,orig.T,O,R); cachedR=R; }
    Maps C; renderMesh(s.V,s.T,C,R); return visualSSIM(O,C)>=need;
}
static bool tryStructuredGrid(MeshOut&ans){
    if(elapsed()>2.0) return false;
    int Rg,Vg;
    if(detectLatLongGrid(Rg,Vg)){
        vector<pair<int,int>> cand; for(int r:{3,4,5,6,8,10,12,14,16,18,20,24,28,32,40,48,56,64}) for(int l:{12,16,20,24,28,32,36,40,48,56,64,80,96,128,160}) cand.push_back({2+r*l,r*1000+l}); sort(cand.begin(),cand.end()); MeshOut cur; int R=(N<=10000?1024:512); double need=(R==1024?0.9007:0.912); for(auto pr:cand){ if(elapsed()>15.5) break; int r=pr.second/1000,l=pr.second%1000; if(!makeLatLongGrid(Rg,Vg,r,l,cur)) continue; if(visualPass(cur,R,need)){ ans=std::move(cur); return true; } }
    }
    int U,Vv;
    if(detectTorusGrid(U,Vv)){
        vector<pair<int,int>> cand; for(int v2:{6,8,10,12,14,16,18,20,24,28,32,36,40,48,56,64,80,96,128}){ double base=(double)U/max(1,Vv)*v2; for(double m:{0.45,0.6,0.75,0.9,1.05,1.25,1.5,1.8}){ int u2=max(4,(int)llround(base*m)); u2=min(u2,U); if(v2<=Vv) cand.push_back({u2*v2,u2*1000+v2}); }} sort(cand.begin(),cand.end()); cand.erase(unique(cand.begin(),cand.end()),cand.end()); MeshOut cur,cur2; int R=(N<=10000?1024:512); double need=(R==1024?0.9007:0.912); for(auto pr:cand){ if(elapsed()>15.5) break; int u=pr.second/1000,v=pr.second%1000; if(!makeTorusGrid(U,Vv,u,v,cur,false)) continue; if(visualPass(cur,R,need)){ ans=std::move(cur); return true; } if(makeTorusGrid(U,Vv,u,v,cur2,true) && visualPass(cur2,R,need)){ ans=std::move(cur2); return true; } }
    }
    return false;
}

static void writeOut(const MeshOut&o){
    string out; out.reserve(o.V.size()*44+o.T.size()*26+64); char line[128]; auto add=[&](const char*s,int n){ if(out.size()+n>(1<<20)){ fwrite(out.data(),1,out.size(),stdout); out.clear(); } out.append(s,s+n); };
    int l=snprintf(line,sizeof(line),"%d %d\n",(int)o.V.size(),(int)o.T.size()); add(line,l);
    for(const Vec3&p:o.V){ l=snprintf(line,sizeof(line),"v %.15g %.15g %.15g\n",p.x,p.y,p.z); add(line,l); }
    for(const Tri&t:o.T){ l=snprintf(line,sizeof(line),"f %d %d %d\n",t.a+1,t.b+1,t.c+1); add(line,l); }
    if(!out.empty()) fwrite(out.data(),1,out.size(),stdout);
}

static bool trySampleBox(MeshOut&out){
    if(N>32) return false;
    vector<Vec3> C; C.reserve(8);
    for(int sx=0;sx<2;sx++) for(int sy=0;sy<2;sy++) for(int sz=0;sz<2;sz++) C.push_back({sx?BBmax.x:BBmin.x, sy?BBmax.y:BBmin.y, sz?BBmax.z:BBmin.z});
    vector<int> rep(8,-1); for(int k=0;k<8;k++){ double best=1e100; for(int i=0;i<N;i++){ double d=norm2(Orig[i]-C[k]); if(d<best){best=d;rep[k]=i;} } if(best>Tau2*1.0001) return false; }
    for(const Vec3&p:Orig){ double best=1e100; for(int k=0;k<8;k++) best=min(best,norm2(p-Orig[rep[k]])); if(best>Tau2*1.0001) return false; }
    vector<int> id(N,-1); out.V.clear(); for(int k=0;k<8;k++){ if(id[rep[k]]<0){ id[rep[k]]=(int)out.V.size(); out.V.push_back(Orig[rep[k]]); } }
    if(out.V.size()!=8) return false;
    int t[12][3]={{0,4,6},{0,6,2},{1,3,7},{1,7,5},{0,1,5},{0,5,4},{2,6,7},{2,7,3},{0,2,3},{0,3,1},{4,5,7},{4,7,6}};
    out.T.clear(); for(auto &q:t) out.T.push_back({q[0],q[1],q[2]}); return true;
}

static MeshOut chooseByRendering(const vector<MeshOut>&snaps){
    MeshOut orig=makeOriginalOut(); if(snaps.empty()) return orig;
    // Evaluate from most compressed to least.  High resolution for small/medium; high-N uses stricter low-res guard to stay under 21s.
    int R = (N<=70000?1024:(N<=300000?512:256));
    double need = (R==1024?0.9008:(R==512?0.925:0.948));
    if(elapsed()>18.6) return snaps.back().V.empty()?orig:snaps.back();
    Maps O; renderMesh(orig.V,orig.T,O,R);
    MeshOut best=orig; double bestRatio=1.0;
    for(const MeshOut&s:snaps){
        if(s.V.empty()||s.T.empty()||s.V.size()>=bestRatio*N) continue;
        if(elapsed()>20.1) break;
        Maps C; renderMesh(s.V,s.T,C,R); double sc=visualSSIM(O,C);
        if(sc>=need){ best=s; bestRatio=(double)s.V.size()/max(1,N); break; }
    }
    return best;
}

static vector<MeshOut> aggressiveCollapseSnapshots(){
    vector<double> ratios;
    if(N<1000) ratios={0.20,0.25,0.35,0.50};
    else if(N<8000) ratios={0.055,0.070,0.085,0.10,0.12,0.15,0.20,0.28,0.40,0.55};
    else if(N<80000) ratios={0.045,0.055,0.065,0.075,0.085,0.10,0.12,0.15,0.19,0.25,0.35,0.50};
    else ratios={0.035,0.045,0.055,0.065,0.075,0.09,0.11,0.14,0.18,0.25,0.36,0.52};
    vector<int> targets; for(double r:ratios){ int t=max(4,(int)ceil(N*r)); if(t<N) targets.push_back(t); }
    sort(targets.begin(),targets.end(),greater<int>()); targets.erase(unique(targets.begin(),targets.end()),targets.end());
    vector<uint64_t> edges; edges.reserve((size_t)M*3);
    for(const Face&f:F){ edges.push_back(ekey(f.v[0],f.v[1])); edges.push_back(ekey(f.v[1],f.v[2])); edges.push_back(ekey(f.v[2],f.v[0])); }
    sort(edges.begin(),edges.end()); edges.erase(unique(edges.begin(),edges.end()),edges.end());
    for(uint64_t k:edges) pushEdge((int)(k>>32),(int)(k&0xffffffffu)); vector<uint64_t>().swap(edges);
    vector<MeshOut> snaps; int ti=0; long long pops=0; double hardStop = (N>300000?13.8:14.8);
    int minTarget=targets.empty()?max(4,N/8):targets.back();
    while(aliveV>minTarget && !PQ.empty() && elapsed()<hardStop){
        Node cur=PQ.top(); PQ.pop(); int a=cur.a,b=cur.b; if(a==b||!Alive[a]||!Alive[b]) continue; if(Ver[a]!=cur.va||Ver[b]!=cur.vb){ pushEdge(a,b); continue; }
        tryCollapse(a,b); ++pops;
        while(ti<(int)targets.size() && aliveV<=targets[ti]){ MeshOut s; if(buildSnapshot(s)) snaps.push_back(std::move(s)); ++ti; }
        if((pops&8191)==0 && elapsed()>hardStop) break;
    }
    if(snaps.empty() || (snaps.back().V.size()!=(size_t)aliveV)){ MeshOut s; if(buildSnapshot(s)) snaps.push_back(std::move(s)); }
    sort(snaps.begin(),snaps.end(),[](const MeshOut&a,const MeshOut&b){return a.V.size()<b.V.size();});
    return snaps;
}

int main(){
    T0=chrono::steady_clock::now();
    if(!readInput()) return 0;
    MeshOut sample; if(trySampleBox(sample)){ writeOut(sample); return 0; }
    MeshOut structured; if(tryStructuredGrid(structured)){ writeOut(structured); return 0; }
    vector<MeshOut> snaps=aggressiveCollapseSnapshots();
    MeshOut ans=chooseByRendering(snaps);
    writeOut(ans);
    return 0;
}
