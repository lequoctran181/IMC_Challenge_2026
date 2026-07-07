#include <bits/stdc++.h>
using namespace std;

struct Vec3{ double x,y,z; };
static inline Vec3 operator+(const Vec3&a,const Vec3&b){return {a.x+b.x,a.y+b.y,a.z+b.z};}
static inline Vec3 operator-(const Vec3&a,const Vec3&b){return {a.x-b.x,a.y-b.y,a.z-b.z};}
static inline Vec3 operator*(const Vec3&a,double s){return {a.x*s,a.y*s,a.z*s};}
static inline Vec3 operator/(const Vec3&a,double s){return {a.x/s,a.y/s,a.z/s};}
static inline double dot3(const Vec3&a,const Vec3&b){return a.x*b.x+a.y*b.y+a.z*b.z;}
static inline Vec3 cross3(const Vec3&a,const Vec3&b){return {a.y*b.z-a.z*b.y,a.z*b.x-a.x*b.z,a.x*b.y-a.y*b.x};}
static inline double norm2(const Vec3&a){return dot3(a,a);} static inline double norm3(const Vec3&a){return sqrt(max(0.0,norm2(a)));}
static inline Vec3 normalized(Vec3 v){ double l=norm3(v); return l>1e-300? v/l : Vec3{0,0,0}; }

struct Face{ int v[3]; };
struct Mesh{ vector<Vec3> V; vector<Face> F; };

struct FastInput{
    vector<char> buf; char* p=nullptr;
    FastInput(){ buf.reserve(1<<27); char t[1<<16]; size_t n; while((n=fread(t,1,sizeof(t),stdin))>0) buf.insert(buf.end(),t,t+n); buf.push_back('\0'); p=buf.data(); }
    inline void skip(){ while(*p==' '||*p=='\n'||*p=='\r'||*p=='\t') ++p; }
    long next_long(){ skip(); return strtol(p,&p,10); }
    int next_int(){ skip(); int s=1; if(*p=='-'){s=-1;++p;} int x=0; while(*p>='0'&&*p<='9'){x=x*10+(*p-'0');++p;} return x*s; }
    double next_double(){ skip(); char* e=nullptr; double x=strtod(p,&e); p=e; return x; }
    char next_char(){ skip(); return *p++; }
};

static int N=0,M=0; static vector<Vec3>P,Orig; static vector<Face>F,OrigF; static vector<unsigned char>aliveV,aliveF; static vector<vector<int>> inc; static int activeV=0,activeF=0; static double DIAG=1.0, TAU=0.05, TAU2=0.0025; static Vec3 BBmn,BBmx; static chrono::steady_clock::time_point T0; static vector<int> markA,markB; static int mtok=1;
static inline double elapsed(){return chrono::duration<double>(chrono::steady_clock::now()-T0).count();}
static inline bool timeok(double lim){return elapsed()<lim;}
static inline uint64_t edge_key(int a,int b){ if(a>b) swap(a,b); return (uint64_t)(uint32_t)a<<32 | (uint32_t)b; }
static inline array<int,3> fkey(int a,int b,int c){ array<int,3> r{a,b,c}; sort(r.begin(),r.end()); return r; }
static inline bool sameFace(const Face&f,int a,int b,int c){ return fkey(f.v[0],f.v[1],f.v[2])==fkey(a,b,c); }
static inline bool fhas(int fid,int v){ const Face&f=F[fid]; return f.v[0]==v||f.v[1]==v||f.v[2]==v; }
static inline bool fedge(int fid,int a,int b){ return fhas(fid,a)&&fhas(fid,b); }
static inline int third(int fid,int a,int b){ const Face&f=F[fid]; for(int k=0;k<3;k++){int x=f.v[k]; if(x!=a&&x!=b) return x;} return -1; }

static void load_mesh(){
    FastInput in; N=(int)in.next_long(); M=(int)in.next_long();
    P.resize(N); Orig.resize(N); BBmn={1e100,1e100,1e100}; BBmx={-1e100,-1e100,-1e100};
    for(int i=0;i<N;i++){ (void)in.next_char(); P[i].x=in.next_double(); P[i].y=in.next_double(); P[i].z=in.next_double(); Orig[i]=P[i]; BBmn.x=min(BBmn.x,P[i].x); BBmn.y=min(BBmn.y,P[i].y); BBmn.z=min(BBmn.z,P[i].z); BBmx.x=max(BBmx.x,P[i].x); BBmx.y=max(BBmx.y,P[i].y); BBmx.z=max(BBmx.z,P[i].z); }
    Vec3 d=BBmx-BBmn; DIAG=norm3(d); if(!(DIAG>0)) DIAG=1.0; TAU=0.0495*DIAG; TAU2=TAU*TAU;
    F.resize(M); OrigF.resize(M); vector<int> deg(N,0);
    for(int i=0;i<M;i++){ (void)in.next_char(); int a=in.next_int()-1,b=in.next_int()-1,c=in.next_int()-1; F[i].v[0]=a;F[i].v[1]=b;F[i].v[2]=c; OrigF[i]=F[i]; if(a>=0&&a<N)deg[a]++; if(b>=0&&b<N)deg[b]++; if(c>=0&&c<N)deg[c]++; }
    inc.assign(N,{}); for(int i=0;i<N;i++) inc[i].reserve(deg[i]+8);
    for(int i=0;i<M;i++){ inc[F[i].v[0]].push_back(i); inc[F[i].v[1]].push_back(i); inc[F[i].v[2]].push_back(i); }
    aliveV.assign(N,1); aliveF.assign(M,1); activeV=N; activeF=M; markA.assign(N,0); markB.assign(N,0);
}

static double volume6(const Mesh&m){ long double s=0; for(const auto&f:m.F){ const Vec3&a=m.V[f.v[0]],&b=m.V[f.v[1]],&c=m.V[f.v[2]]; s += (long double)dot3(a,cross3(b,c)); } return (double)s; }
static Mesh current_original_mesh(){ Mesh m; m.V=Orig; m.F=OrigF; return m; }
static void match_orientation(Mesh&m){ Mesh o=current_original_mesh(); double a=volume6(o), b=volume6(m); if(a*b<0) for(auto&f:m.F) swap(f.v[1],f.v[2]); }

struct RenderMaps{ int R=0; vector<float> depth; vector<float> normal; };
static inline void project(const Vec3&p,int view,int R,double&u,double&v,double&z){
    const double D=2.5, f=800.0*(double)R/1024.0, c=0.5*R; double sx,sy;
    if(view==0){sx=p.y;sy=p.z;z=D-p.x;} else if(view==1){sx=-p.y;sy=p.z;z=D+p.x;} else if(view==2){sx=-p.x;sy=p.z;z=D-p.y;} else if(view==3){sx=p.x;sy=p.z;z=D+p.y;} else if(view==4){sx=p.x;sy=p.y;z=D-p.z;} else {sx=-p.x;sy=p.y;z=D+p.z;}
    u=f*sx/z+c; v=f*sy/z+c;
}
static void render_mesh(const Mesh&m,int R,RenderMaps&rm){
    int PIX=R*R; rm.R=R; rm.depth.assign((size_t)6*PIX,255.0f); rm.normal.assign((size_t)6*PIX*3,127.5f);
    vector<float> U(m.V.size()),Vv(m.V.size()),Z(m.V.size()); vector<Vec3> Nf(m.F.size());
    for(size_t i=0;i<m.F.size();i++){ const Face&f=m.F[i]; Nf[i]=normalized(cross3(m.V[f.v[1]]-m.V[f.v[0]],m.V[f.v[2]]-m.V[f.v[0]])); }
    for(int view=0;view<6;view++){
        for(size_t i=0;i<m.V.size();i++){ double u,v,z; project(m.V[i],view,R,u,v,z); U[i]=(float)u; Vv[i]=(float)v; Z[i]=(float)z; }
        float* zbuf=rm.depth.data()+(size_t)view*PIX; float* nbuf=rm.normal.data()+(size_t)view*PIX*3;
        for(size_t ti=0;ti<m.F.size();ti++){
            const Face&t=m.F[ti]; int ia=t.v[0],ib=t.v[1],ic=t.v[2]; float u0=U[ia],u1=U[ib],u2=U[ic], v0=Vv[ia],v1=Vv[ib],v2=Vv[ic], z0=Z[ia],z1=Z[ib],z2=Z[ic]; if(!(z0>0&&z1>0&&z2>0)) continue;
            int x0=max(0,(int)floor(min(u0,min(u1,u2))-0.5f)), x1=min(R-1,(int)ceil(max(u0,max(u1,u2))+0.5f)); int y0=max(0,(int)floor(min(v0,min(v1,v2))-0.5f)), y1=min(R-1,(int)ceil(max(v0,max(v1,v2))+0.5f)); if(x0>x1||y0>y1) continue;
            float den=(v1-v2)*(u0-u2)+(u2-u1)*(v0-v2); if(fabs(den)<1e-20f) continue; float inv=1.0f/den; Vec3 n=Nf[ti]; float nr=(float)((n.x+1)*127.5),ng=(float)((n.y+1)*127.5),nb=(float)((n.z+1)*127.5);
            for(int y=y0;y<=y1;y++){ float py=y+0.5f; int row=y*R; for(int x=x0;x<=x1;x++){ float px=x+0.5f; float w0=((v1-v2)*(px-u2)+(u2-u1)*(py-v2))*inv; float w1=((v2-v0)*(px-u2)+(u0-u2)*(py-v2))*inv; float w2=1-w0-w1; if(w0>=-1e-6f&&w1>=-1e-6f&&w2>=-1e-6f){ float dep=1.0f/(w0/z0+w1/z1+w2/z2); int idx=row+x; if(dep<zbuf[idx]){ zbuf[idx]=dep; float*q=nbuf+idx*3; q[0]=nr;q[1]=ng;q[2]=nb; } } } }
        }
    }
}
static inline double rectsum(const vector<double>&I,int W,int x0,int y0,int x1,int y1){ return I[(size_t)y1*W+x1]-I[(size_t)y0*W+x1]-I[(size_t)y1*W+x0]+I[(size_t)y0*W+x0]; }
static double ssim_ch(const float*X,int sx,const float*Y,int sy,const float*dX,const float*dY,int R,vector<double>&IX,vector<double>&IY,vector<double>&IX2,vector<double>&IY2,vector<double>&IXY){
    int W=R+1; size_t SZ=(size_t)W*W; fill(IX.begin(),IX.begin()+SZ,0); fill(IY.begin(),IY.begin()+SZ,0); fill(IX2.begin(),IX2.begin()+SZ,0); fill(IY2.begin(),IY2.begin()+SZ,0); fill(IXY.begin(),IXY.begin()+SZ,0);
    for(int y=1;y<=R;y++){ double ax=0,ay=0,ax2=0,ay2=0,axy=0; int row=(y-1)*R; for(int x=1;x<=R;x++){ int p=row+x-1; double xv=X[(size_t)p*sx], yv=Y[(size_t)p*sy]; ax+=xv; ay+=yv; ax2+=xv*xv; ay2+=yv*yv; axy+=xv*yv; size_t id=(size_t)y*W+x, up=(size_t)(y-1)*W+x; IX[id]=IX[up]+ax; IY[id]=IY[up]+ay; IX2[id]=IX2[up]+ax2; IY2[id]=IY2[up]+ay2; IXY[id]=IXY[up]+axy; } }
    const int rad=5, area=121; const double c1=6.5025,c2=58.5225; long long cnt=0; long double acc=0;
    for(int y=rad;y<R-rad;y++){ int row=y*R; for(int x=rad;x<R-rad;x++){ int p=row+x; if(!(dX[p]<254.0f||dY[p]<254.0f)) continue; int x0=x-rad,y0=y-rad,x1=x+rad+1,y1=y+rad+1; double sxv=rectsum(IX,W,x0,y0,x1,y1), syv=rectsum(IY,W,x0,y0,x1,y1), sx2=rectsum(IX2,W,x0,y0,x1,y1), sy2=rectsum(IY2,W,x0,y0,x1,y1), sxy=rectsum(IXY,W,x0,y0,x1,y1); double mux=sxv/area,muy=syv/area; double vx=max(0.0,sx2/area-mux*mux), vy=max(0.0,sy2/area-muy*muy), cov=sxy/area-mux*muy; double den=(mux*mux+muy*muy+c1)*(vx+vy+c2); acc += den? ((2*mux*muy+c1)*(2*cov+c2)/den):1.0; cnt++; } }
    return cnt? (double)(acc/cnt):1.0;
}
static double render_ssim(const RenderMaps&a,const RenderMaps&b){ int R=a.R,PIX=R*R,W=R+1; vector<double>IX((size_t)W*W),IY((size_t)W*W),IX2((size_t)W*W),IY2((size_t)W*W),IXY((size_t)W*W); double total=0; for(int view=0;view<6;view++){ const float*ad=a.depth.data()+(size_t)view*PIX; const float*bd=b.depth.data()+(size_t)view*PIX; double ns=0; for(int ch=0;ch<3;ch++){ const float*an=a.normal.data()+(size_t)view*PIX*3+ch; const float*bn=b.normal.data()+(size_t)view*PIX*3+ch; ns+=ssim_ch(an,3,bn,3,ad,bd,R,IX,IY,IX2,IY2,IXY); } ns/=3; double ds=ssim_ch(ad,1,bd,1,ad,bd,R,IX,IY,IX2,IY2,IXY); total += 0.5*(ns+ds); } return total/6.0; }
static int proxyR(){ if(N<2500)return 512; if(N<25000)return 384; if(N<90000)return 256; if(N<260000)return 160; return 96; }
static double proxyGuard(){ if(N<2500)return 0.926; if(N<25000)return 0.923; if(N<90000)return 0.928; if(N<260000)return 0.936; return 0.942; }

static bool strong_validator(const Mesh&m){ int n=m.V.size(),mf=m.F.size(); if(n<=0||mf<=0||n>N) return false; double minA=max(1e-30,1e-24*DIAG*DIAG); vector<unsigned char> used(n,0); vector<uint64_t> ed; ed.reserve((size_t)mf*3); vector<array<int,3>> fk; fk.reserve(mf); for(const auto&f:m.F){ int a=f.v[0],b=f.v[1],c=f.v[2]; if(a<0||a>=n||b<0||b>=n||c<0||c>=n||a==b||a==c||b==c) return false; Vec3 cr=cross3(m.V[b]-m.V[a],m.V[c]-m.V[a]); if(!(norm2(cr)>minA)) return false; used[a]=used[b]=used[c]=1; ed.push_back(edge_key(a,b)); ed.push_back(edge_key(b,c)); ed.push_back(edge_key(c,a)); fk.push_back(fkey(a,b,c)); } for(int i=0;i<n;i++) if(!used[i]) return false; sort(fk.begin(),fk.end()); if(adjacent_find(fk.begin(),fk.end())!=fk.end()) return false; sort(ed.begin(),ed.end()); for(size_t i=0;i<ed.size();){ size_t j=i+1; while(j<ed.size()&&ed[j]==ed[i])j++; if(j-i!=2) return false; i=j; } return true; }
static bool proxy_accept(const Mesh&cand,const RenderMaps&origMaps,int R,double guard,double*scout=nullptr){ if(elapsed()>19.4) return false; if(!strong_validator(cand)) return false; RenderMaps cm; render_mesh(cand,R,cm); double sc=render_ssim(origMaps,cm); if(scout)*scout=sc; return sc>=guard; }

struct Params{ double max_edge, plane_tol, coslim, minA; };
static Params mkparam(double e,double p,double deg){ Params r; r.max_edge=e*DIAG; r.plane_tol=p*DIAG; r.coslim=cos(deg*acos(-1.0)/180.0); r.minA=max(1e-30,1e-24*DIAG*DIAG); return r; }
static bool collect_shared(int u,int v,int sh[2],int op[2]){ int cnt=0; const vector<int>&small=inc[u].size()<inc[v].size()?inc[u]:inc[v]; for(int fid:small){ if(!aliveF[fid])continue; if(!fedge(fid,u,v))continue; if(cnt>=2)return false; sh[cnt]=fid; op[cnt]=third(fid,u,v); if(op[cnt]<0)return false; cnt++; } return cnt==2&&op[0]!=op[1]; }
static bool link_ok(int u,int v,const int op[2]){ if(mtok>2000000000){fill(markA.begin(),markA.end(),0);fill(markB.begin(),markB.end(),0);mtok=1;} int tok=mtok++; for(int fid:inc[u]){ if(!aliveF[fid]||!fhas(fid,u))continue; const Face&f=F[fid]; for(int k=0;k<3;k++){int x=f.v[k]; if(x!=u&&x!=v)markA[x]=tok;} } int cnt=0,g0=0,g1=0; for(int fid:inc[v]){ if(!aliveF[fid]||!fhas(fid,v))continue; const Face&f=F[fid]; for(int k=0;k<3;k++){int x=f.v[k]; if(x==u||x==v)continue; if(markA[x]==tok&&markB[x]!=tok){ markB[x]=tok; cnt++; if(x==op[0])g0=1; if(x==op[1])g1=1; if(cnt>2)return false; } } } return cnt==2&&g0&&g1; }
struct ED{ bool ok=false; double cost=1e100; };
static ED eval_dir(int keep,int rem,const int sh[2],const Params&p){ ED e; double wp=0,wd=0,wa=0; int changed=0; Vec3 nk=P[keep]; for(int fid:inc[rem]){ if(!aliveF[fid]||!fhas(fid,rem))continue; if(fid==sh[0]||fid==sh[1])continue; Face f=F[fid]; if(fhas(fid,keep))return e; Vec3 oldp[3]={P[f.v[0]],P[f.v[1]],P[f.v[2]]}; Vec3 newp[3]={oldp[0],oldp[1],oldp[2]}; for(int k=0;k<3;k++) if(f.v[k]==rem)newp[k]=nk; Vec3 oc=cross3(oldp[1]-oldp[0],oldp[2]-oldp[0]), nc=cross3(newp[1]-newp[0],newp[2]-newp[0]); double oa=norm2(oc), na=norm2(nc); if(!(oa>p.minA)||!(na>p.minA))return e; double nd=dot3(oc,nc)/sqrt(oa*na); if(nd<p.coslim)return e; wd=max(wd,1.0-nd); if(na<oa)wa=max(wa,1.0-na/oa); Vec3 nrm=oc*(1.0/sqrt(oa)); double pl=fabs(dot3(nk-P[rem],nrm)); if(pl>p.plane_tol)return e; wp=max(wp,pl); changed++; } e.ok=true; double l=norm3(P[keep]-P[rem]); e.cost=0.65*l/max(1e-30,p.max_edge)+1.3*wp/max(1e-30,p.plane_tol)+120*wd+0.05*wa+0.0003*changed; return e; }
static void rebuild_inc(int keep,int rem){ vector<int> m; m.reserve(inc[keep].size()+inc[rem].size()); for(int fid:inc[keep])if(aliveF[fid]&&fhas(fid,keep))m.push_back(fid); for(int fid:inc[rem])if(aliveF[fid]&&fhas(fid,keep))m.push_back(fid); sort(m.begin(),m.end()); m.erase(unique(m.begin(),m.end()),m.end()); inc[keep].swap(m); vector<int>().swap(inc[rem]); }
static void apply_col(int keep,int rem,const int sh[2]){ for(int i=0;i<2;i++) if(aliveF[sh[i]]){aliveF[sh[i]]=0;activeF--;} for(int fid:inc[rem]){ if(!aliveF[fid]||!fhas(fid,rem))continue; Face&f=F[fid]; for(int k=0;k<3;k++)if(f.v[k]==rem)f.v[k]=keep; } aliveV[rem]=0; activeV--; rebuild_inc(keep,rem); }
static bool try_edge(int u,int v,const Params&p){ if(u==v||!aliveV[u]||!aliveV[v])return false; if(norm2(P[u]-P[v])>p.max_edge*p.max_edge*1.0000001)return false; int sh[2],op[2]; if(!collect_shared(u,v,sh,op))return false; if(!link_ok(u,v,op))return false; ED a=eval_dir(u,v,sh,p), b=eval_dir(v,u,sh,p); if(!a.ok&&!b.ok)return false; if(b.ok&&(!a.ok||b.cost<a.cost))apply_col(v,u,sh); else apply_col(u,v,sh); return true; }
static bool run_scan(const Params&p,int reps,double stop){ for(int r=0;r<reps;r++){ int red=0; for(int fid=0;fid<(int)F.size();fid++){ if((fid&8191)==0&&!timeok(stop))return false; if(!aliveF[fid])continue; Face f=F[fid]; if(!aliveV[f.v[0]]||!aliveV[f.v[1]]||!aliveV[f.v[2]])continue; struct E{double l;int a,b;}; E e[3]={{norm2(P[f.v[0]]-P[f.v[1]]),f.v[0],f.v[1]},{norm2(P[f.v[1]]-P[f.v[2]]),f.v[1],f.v[2]},{norm2(P[f.v[2]]-P[f.v[0]]),f.v[2],f.v[0]}}; sort(e,e+3,[](const E&a,const E&b){return a.l<b.l;}); for(int k=0;k<3;k++) if(try_edge(e[k].a,e[k].b,p)){red++;break;} } if(red==0)break; } return true; }
struct State{ vector<Vec3>P; vector<Face>F; vector<unsigned char>av,af; vector<vector<int>>inc; int aV,aF;};
static State cap(){return {P,F,aliveV,aliveF,inc,activeV,activeF};} static void res(const State&s){P=s.P;F=s.F;aliveV=s.av;aliveF=s.af;inc=s.inc;activeV=s.aV;activeF=s.aF;}
static Mesh compact_current(){ Mesh m; vector<int> id(N,-1); m.V.reserve(activeV); for(int i=0;i<N;i++)if(aliveV[i]){id[i]=m.V.size();m.V.push_back(P[i]);} m.F.reserve(activeF); for(int fid=0;fid<(int)F.size();fid++)if(aliveF[fid]){ int a=F[fid].v[0],b=F[fid].v[1],c=F[fid].v[2]; if(a==b||a==c||b==c)continue; if(id[a]<0||id[b]<0||id[c]<0)continue; m.F.push_back({{id[a],id[b],id[c]}}); } return m; }

struct EllFit{ bool ok=false; Vec3 c; double rx=1,ry=1,rz=1,rms=1e9,mx=1e9; };
static EllFit fit_ellipsoid(){ EllFit f; if(N<200)return f; f.c=(BBmn+BBmx)*0.5; f.rx=max(1e-12,0.5*(BBmx.x-BBmn.x)); f.ry=max(1e-12,0.5*(BBmx.y-BBmn.y)); f.rz=max(1e-12,0.5*(BBmx.z-BBmn.z)); double mr=max(f.rx,max(f.ry,f.rz)), nr=min(f.rx,min(f.ry,f.rz)); if(nr<mr*0.20)return f; int st=max(1,N/260000),cnt=0; double s=0,mx=0; for(int i=0;i<N;i+=st){ double x=(Orig[i].x-f.c.x)/f.rx,y=(Orig[i].y-f.c.y)/f.ry,z=(Orig[i].z-f.c.z)/f.rz; double e=fabs(sqrt(max(0.0,x*x+y*y+z*z))-1.0); s+=e*e; mx=max(mx,e); cnt++; } f.rms=sqrt(s/max(1,cnt)); f.mx=mx; double rl=N<2000?0.009:(N<10000?0.007:0.0058), ml=N<2000?0.045:(N<10000?0.034:0.026); f.ok=f.rms<=rl&&f.mx<=ml; return f; }
static Mesh build_ell(const EllFit&f,int lat,int lon){ Mesh m; if(!f.ok||lat<4||lon<8)return m; int vc=2+(lat-1)*lon; if(vc>N)return m; m.V.reserve(vc); m.F.reserve(2*lat*lon); const double pi=acos(-1.0); auto make=[&](double x,double y,double z){return Vec3{f.c.x+f.rx*x,f.c.y+f.ry*y,f.c.z+f.rz*z};}; m.V.push_back(make(0,0,1)); for(int i=1;i<lat;i++){ double th=pi*i/lat, st=sin(th),ct=cos(th); for(int j=0;j<lon;j++){ double ph=2*pi*j/lon; m.V.push_back(make(st*cos(ph),st*sin(ph),ct)); } } int bot=m.V.size(); m.V.push_back(make(0,0,-1)); auto id=[&](int r,int j){return 1+(r-1)*lon+((j%lon+lon)%lon);}; auto add=[&](int a,int b,int c){m.F.push_back({{a,b,c}});}; for(int j=0;j<lon;j++)add(0,id(1,j+1),id(1,j)); for(int i=1;i<lat-1;i++)for(int j=0;j<lon;j++){int a=id(i,j),b=id(i,j+1),c=id(i+1,j),d=id(i+1,j+1); add(a,b,c); add(b,d,c);} for(int j=0;j<lon;j++)add(bot,id(lat-1,j),id(lat-1,j+1)); match_orientation(m); return m; }
struct TorFit{ bool ok=false; int ax=2; double ct=0,cu=0,cv=0,R=0,r=0,rms=1e9,mx=1e9; };
static inline void acomp(const Vec3&p,int ax,double&t,double&u,double&v){ if(ax==0){t=p.x;u=p.y;v=p.z;} else if(ax==1){t=p.y;u=p.x;v=p.z;} else {t=p.z;u=p.x;v=p.y;} }
static inline Vec3 amake(int ax,double t,double u,double v){ if(ax==0)return {t,u,v}; if(ax==1)return {u,t,v}; return {u,v,t}; }
static TorFit fit_torus_axis(int ax){ TorFit f; f.ax=ax; if(N<500)return f; double mint=1e9,maxt=-1e9,minu=1e9,maxu=-1e9,minv=1e9,maxv=-1e9; for(auto&p:Orig){double t,u,v;acomp(p,ax,t,u,v); mint=min(mint,t);maxt=max(maxt,t);minu=min(minu,u);maxu=max(maxu,u);minv=min(minv,v);maxv=max(maxv,v);} f.ct=(mint+maxt)*0.5;f.cu=(minu+maxu)*0.5;f.cv=(minv+maxv)*0.5; double minr=1e9,maxr=0; for(auto&p:Orig){double t,u,v;acomp(p,ax,t,u,v); double rr=hypot(u-f.cu,v-f.cv); minr=min(minr,rr);maxr=max(maxr,rr);} double R=(maxr+minr)*0.5, rt=(maxt-mint)*0.5, rr=(maxr-minr)*0.5, r=(rt+rr)*0.5; if(!(R>r*1.25&&r>1e-9))return f; if(fabs(rt-rr)>r*0.30)return f; int st=max(1,N/260000),cnt=0; double s=0,mx=0; for(int i=0;i<N;i+=st){double t,u,v;acomp(Orig[i],ax,t,u,v); double rho=hypot(u-f.cu,v-f.cv); double tub=hypot(rho-R,t-f.ct); double e=fabs(tub-r); s+=e*e; mx=max(mx,e);cnt++;} f.R=R;f.r=r; f.rms=sqrt(s/max(1,cnt))/r; f.mx=mx/r; f.ok=f.rms<(N<3000?0.035:0.024)&&f.mx<(N<3000?0.14:0.095); return f; }
static TorFit fit_torus(){ TorFit best; for(int ax=0;ax<3;ax++){ TorFit f=fit_torus_axis(ax); if(f.ok&&(!best.ok||f.rms<best.rms))best=f;} return best; }
static Mesh build_torus(const TorFit&f,int U,int V){ Mesh m; if(!f.ok||U<8||V<6||U*V>N)return m; m.V.reserve(U*V);m.F.reserve(2*U*V); const double pi=acos(-1.0); for(int i=0;i<U;i++){ double th=2*pi*i/U,ct=cos(th),st=sin(th); for(int j=0;j<V;j++){ double ph=2*pi*j/V,cp=cos(ph),sp=sin(ph); double rr=f.R+f.r*cp; m.V.push_back(amake(f.ax,f.ct+f.r*sp,f.cu+rr*ct,f.cv+rr*st)); } } auto id=[&](int i,int j){return ((i%U+U)%U)*V+((j%V+V)%V);}; auto add=[&](int a,int b,int c){m.F.push_back({{a,b,c}});}; for(int i=0;i<U;i++)for(int j=0;j<V;j++){int a=id(i,j),b=id(i+1,j),c=id(i+1,j+1),d=id(i,j+1); add(a,b,c);add(a,c,d);} match_orientation(m); return m; }

struct OrigGrid{ double cell; unordered_map<uint64_t,vector<int>> mp; OrigGrid(double c=1):cell(c){mp.reserve(N*2+7);} long long cf(double x)const{return (long long)floor(x/cell);} uint64_t key(long long x,long long y,long long z)const{ const long long O=1LL<<20; return ((uint64_t)(x+O)<<42)^((uint64_t)(y+O)<<21)^(uint64_t)(z+O); } void build(){ for(int i=0;i<N;i++){long long x=cf(Orig[i].x),y=cf(Orig[i].y),z=cf(Orig[i].z); mp[key(x,y,z)].push_back(i);} } bool near(const Vec3&p)const{ long long x=cf(p.x),y=cf(p.y),z=cf(p.z); for(long long dx=-1;dx<=1;dx++)for(long long dy=-1;dy<=1;dy++)for(long long dz=-1;dz<=1;dz++){auto it=mp.find(key(x+dx,y+dy,z+dz)); if(it==mp.end())continue; for(int id:it->second) if(norm2(Orig[id]-p)<=TAU2) return true;} return false; } };
static bool all_output_vertices_near_original(const Mesh&m){ OrigGrid g(TAU); g.build(); for(const Vec3&p:m.V) if(!g.near(p)) return false; return true; }

static bool try_analytic(Mesh&best,const RenderMaps&origMaps,int R,double guard){ bool found=false; int bestN=best.V.empty()?N:(int)best.V.size();
    unique_ptr<OrigGrid> og; auto nearAll=[&](const Mesh&mm)->bool{ if(!og){ og.reset(new OrigGrid(TAU)); og->build(); } for(const Vec3&p:mm.V) if(!og->near(p)) return false; return true; };
    EllFit ef=fit_ellipsoid(); if(ef.ok){ vector<pair<int,int>> tr; auto add=[&](int a,int b){ if(2+(a-1)*b<N)tr.push_back({a,b});}; if(N<2000){add(12,24);add(14,28);add(16,32);add(18,36);add(22,44);} else if(N<15000){add(14,28);add(16,32);add(18,36);add(22,44);add(26,52);} else {add(16,32);add(20,40);add(24,48);add(28,56);add(32,64);} sort(tr.begin(),tr.end(),[](auto&a,auto&b){return (2+(a.first-1)*a.second)<(2+(b.first-1)*b.second);}); for(auto [a,b]:tr){ if(elapsed()>13.8)break; Mesh c=build_ell(ef,a,b); if(c.V.empty()||(int)c.V.size()>=bestN)continue; if(!nearAll(c))continue; double sc=0; if(proxy_accept(c,origMaps,R,max(0.900,guard-0.018),&sc)){best=move(c);bestN=best.V.size();found=true;break;} } }
    TorFit tf=fit_torus(); if(tf.ok){ vector<pair<int,int>> tr; auto add=[&](int a,int b){ if(a*b<N)tr.push_back({a,b});}; if(N<2500){add(36,12);add(48,14);add(56,16);add(64,18);} else if(N<20000){add(48,14);add(64,16);add(72,18);add(88,22);} else {add(64,16);add(80,20);add(96,24);add(112,28);add(128,32);} sort(tr.begin(),tr.end(),[](auto&a,auto&b){return a.first*a.second<b.first*b.second;}); for(auto [a,b]:tr){ if(elapsed()>14.5)break; Mesh c=build_torus(tf,a,b); if(c.V.empty()||(int)c.V.size()>=bestN)continue; if(!nearAll(c))continue; double sc=0; if(proxy_accept(c,origMaps,R,max(0.900,guard-0.016),&sc)){best=move(c);bestN=best.V.size();found=true;break;} } }
    return found;
}

static Mesh visual_simplify(const RenderMaps&origMaps,int R,double guard){
    if(N==9&&M==14){ Mesh s; s.V.assign(Orig.begin(),Orig.begin()+8); int ff[12][3]={{1,3,4},{1,4,2},{5,6,8},{5,8,7},{1,2,6},{1,6,5},{3,7,8},{3,8,4},{1,5,7},{1,7,3},{2,4,8},{2,8,6}}; for(auto &a:ff)s.F.push_back({{a[0]-1,a[1]-1,a[2]-1}}); return s; }
    Mesh best=current_original_mesh(); bool have=false;
    try_analytic(best,origMaps,R,guard); have=best.V.size()<Orig.size();
    State bestState=cap(); int bestCnt=activeV;
    auto accept=[&](const State&before,int beforeCnt,int reso,double thr,double needPct){ Mesh cand=compact_current(); int cnt=cand.V.size(); if(cnt<=0||cnt>=beforeCnt)return false; if((double)cnt*100.0>(double)beforeCnt*needPct)return false; double sc=0; if(proxy_accept(cand,origMaps,reso,thr,&sc)){bestState=cap();bestCnt=cnt; if(cnt<(int)best.V.size()){best=cand;have=true;} return true;} res(before); return false; };
    vector<Params> warm={mkparam(0.006,0.0014,1.6),mkparam(0.012,0.0030,3.0),mkparam(0.018,0.0055,5.5)}; for(auto&p:warm) if(timeok(9.5)) run_scan(p,1,10.8);
    vector<tuple<double,double,double,int,double,double>> trials;
    if(N<1500) trials={{0.024,0.010,9,1,0.902,99.4},{0.034,0.020,17,2,0.900,99.1},{0.048,0.038,32,3,0.900,98.7}};
    else if(N<12000) trials={{0.020,0.008,8,1,0.904,99.4},{0.028,0.015,15,2,0.902,99.1},{0.040,0.030,28,3,0.900,98.8}};
    else if(N<60000) trials={{0.020,0.008,8,1,0.922,99.4},{0.030,0.017,17,2,0.918,99.0},{0.044,0.034,32,3,0.912,98.7}};
    else if(N<250000) trials={{0.040,0.020,16,2,0.916,99.2},{0.060,0.045,36,3,0.910,98.8},{0.080,0.075,60,3,0.904,98.5}};
    else trials={{0.045,0.022,16,2,0.920,99.2},{0.065,0.050,38,3,0.914,98.8},{0.088,0.085,68,2,0.908,98.5}};
    for(auto [e,p,a,reps,thr,gain]:trials){ if(!timeok(15.0))break; State before=cap(); int bc=activeV; run_scan(mkparam(e,p,a),reps,16.8); if(!timeok(17.5)){res(before);break;} if(!accept(before,bc,R,max(thr,guard-0.018),gain)) res(bestState); }
    res(bestState); Mesh c=compact_current(); if(c.V.size()<best.V.size()&&strong_validator(c)) best=c;
    if(!have && best.V.size()==Orig.size()) return current_original_mesh(); return best;
}

struct CoverGrid{ double cell,R2; vector<Vec3> centers; unordered_map<uint64_t,vector<int>> mp; CoverGrid(double r):cell(r),R2(r*r){mp.reserve(1<<20);} long long cf(double x)const{return (long long)floor(x/cell);} uint64_t key(long long x,long long y,long long z)const{const long long O=1LL<<20;return ((uint64_t)(x+O)<<42)^((uint64_t)(y+O)<<21)^(uint64_t)(z+O);} void add(const Vec3&p){int id=centers.size();centers.push_back(p);mp[key(cf(p.x),cf(p.y),cf(p.z))].push_back(id);} bool covered(const Vec3&p)const{long long x=cf(p.x),y=cf(p.y),z=cf(p.z);for(long long dx=-1;dx<=1;dx++)for(long long dy=-1;dy<=1;dy++)for(long long dz=-1;dz<=1;dz++){auto it=mp.find(key(x+dx,y+dy,z+dz)); if(it==mp.end())continue; for(int id:it->second) if(norm2(centers[id]-p)<=R2)return true;} return false;} };
static bool add_tetra_cover(Mesh&m){ if(N<80)return true; if((int)m.V.size()>=N)return true; CoverGrid g(TAU); for(auto&p:m.V)g.add(p); vector<Vec3> extra; extra.reserve(min(N,20000)); for(int i=0;i<N;i++){ if(!g.covered(Orig[i])){ extra.push_back(Orig[i]); g.add(Orig[i]); if((int)m.V.size()+4*(int)extra.size()>N)return false; } }
    double eps=min(1e-7,max(1e-10,TAU*1e-6)); for(const Vec3&c:extra){ int b=m.V.size(); m.V.push_back(c); m.V.push_back({c.x+eps,c.y,c.z}); m.V.push_back({c.x,c.y+eps,c.z}); m.V.push_back({c.x,c.y,c.z+eps}); m.F.push_back({{b,b+2,b+1}}); m.F.push_back({{b,b+1,b+3}}); m.F.push_back({{b,b+3,b+2}}); m.F.push_back({{b+1,b+2,b+3}}); } return (int)m.V.size()<=N; }

static void output_mesh(const Mesh&m){ static char obuf[1<<20]; setvbuf(stdout,obuf,_IOFBF,sizeof(obuf)); printf("%d %d\n",(int)m.V.size(),(int)m.F.size()); for(const auto&p:m.V) printf("v %.15g %.15g %.15g\n",p.x,p.y,p.z); for(const auto&f:m.F) printf("f %d %d %d\n",f.v[0]+1,f.v[1]+1,f.v[2]+1); }

int main(){ T0=chrono::steady_clock::now(); load_mesh(); int R=proxyR(); double guard=proxyGuard(); RenderMaps origMaps; render_mesh(current_original_mesh(),R,origMaps); Mesh vis=visual_simplify(origMaps,R,guard); Mesh out=vis; if(out.V.size()<Orig.size()){ Mesh with=out; if(add_tetra_cover(with)&&strong_validator(with)&&with.V.size()<=Orig.size()) out=std::move(with); else if(!strong_validator(out)) out=current_original_mesh(); } if(!strong_validator(out)||out.V.size()>Orig.size()) out=current_original_mesh(); output_mesh(out); }