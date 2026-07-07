#include <bits/stdc++.h>
using namespace std;

struct Vec3{ double x=0,y=0,z=0; };
static inline Vec3 operator+(const Vec3&a,const Vec3&b){ return {a.x+b.x,a.y+b.y,a.z+b.z}; }
static inline Vec3 operator-(const Vec3&a,const Vec3&b){ return {a.x-b.x,a.y-b.y,a.z-b.z}; }
static inline Vec3 operator*(const Vec3&a,double s){ return {a.x*s,a.y*s,a.z*s}; }
static inline Vec3 operator/(const Vec3&a,double s){ return {a.x/s,a.y/s,a.z/s}; }
static inline double dot3(const Vec3&a,const Vec3&b){ return a.x*b.x+a.y*b.y+a.z*b.z; }
static inline Vec3 cross3(const Vec3&a,const Vec3&b){ return {a.y*b.z-a.z*b.y,a.z*b.x-a.x*b.z,a.x*b.y-a.y*b.x}; }
static inline double norm2(const Vec3&a){ return dot3(a,a); }
static inline double norm3(const Vec3&a){ return sqrt(norm2(a)); }
static inline Vec3 unit(Vec3 a){ double n=norm3(a); return n>1e-300?a*(1.0/n):Vec3{0,0,1}; }
static inline double clampd(double x,double lo,double hi){ return x<lo?lo:(x>hi?hi:x); }

struct Face{ int v[3]; unsigned char active=1; };
struct Mesh{ vector<Vec3> V; vector<Face> F; };

struct FastInput{
    vector<char> buf; char* p=nullptr;
    FastInput(){ char tmp[1<<16]; size_t n; buf.reserve(1<<27); while((n=fread(tmp,1,sizeof(tmp),stdin))>0) buf.insert(buf.end(),tmp,tmp+n); buf.push_back('\0'); p=buf.data(); }
    inline void ws(){ while(*p==' '||*p=='\n'||*p=='\r'||*p=='\t') ++p; }
    int next_int(){ ws(); return (int)strtol(p,&p,10); }
    double next_double(){ ws(); return strtod(p,&p); }
    char next_char(){ ws(); return *p++; }
};

static int N=0,M=0;
static vector<Vec3> P, Orig;
static vector<Face> F, OrigF;
static vector<unsigned char> alive;
static vector<unsigned char> factive;
static vector<vector<int>> inc;
static vector<double> errv;
static int activeFaces=0;
static double diagLen=1.0, tau=0.05, tau2=0.0025;
static Vec3 bbMin, bbMax;
static chrono::steady_clock::time_point T0;
static vector<int> markv;
static int markToken=1;

static inline double elapsed(){ return chrono::duration<double>(chrono::steady_clock::now()-T0).count(); }
static inline uint64_t ekey(int a,int b){ if(a>b) swap(a,b); return (uint64_t)(uint32_t)a<<32 | (uint32_t)b; }
static inline bool hasv(const Face&f,int v){ return f.v[0]==v||f.v[1]==v||f.v[2]==v; }
static inline int thirdv(const Face&f,int a,int b){ for(int k=0;k<3;k++){ int x=f.v[k]; if(x!=a&&x!=b) return x; } return -1; }
static inline array<int,3> tri_key(int a,int b,int c){ array<int,3> t{a,b,c}; sort(t.begin(),t.end()); return t; }

static void compute_bbox(){
    bbMin={1e100,1e100,1e100}; bbMax={-1e100,-1e100,-1e100};
    for(auto &p:Orig){ bbMin.x=min(bbMin.x,p.x); bbMin.y=min(bbMin.y,p.y); bbMin.z=min(bbMin.z,p.z); bbMax.x=max(bbMax.x,p.x); bbMax.y=max(bbMax.y,p.y); bbMax.z=max(bbMax.z,p.z); }
    diagLen=norm3(bbMax-bbMin); if(!(diagLen>0)) diagLen=1.0; tau=0.05*diagLen*0.999; tau2=tau*tau;
}

static void read_input(){
    FastInput in; N=in.next_int(); M=in.next_int();
    P.resize(N); Orig.resize(N);
    for(int i=0;i<N;i++){ (void)in.next_char(); P[i].x=in.next_double(); P[i].y=in.next_double(); P[i].z=in.next_double(); Orig[i]=P[i]; }
    F.resize(M); OrigF.resize(M); vector<int> deg(N,0);
    for(int i=0;i<M;i++){ (void)in.next_char(); int a=in.next_int()-1,b=in.next_int()-1,c=in.next_int()-1; F[i].v[0]=a;F[i].v[1]=b;F[i].v[2]=c;F[i].active=1; OrigF[i]=F[i]; if(a>=0&&a<N)deg[a]++; if(b>=0&&b<N)deg[b]++; if(c>=0&&c<N)deg[c]++; }
    compute_bbox(); alive.assign(N,1); factive.assign(M,1); errv.assign(N,0.0); inc.assign(N,{}); for(int i=0;i<N;i++) inc[i].reserve(deg[i]+4);
    for(int i=0;i<M;i++){ inc[F[i].v[0]].push_back(i); inc[F[i].v[1]].push_back(i); inc[F[i].v[2]].push_back(i); }
    activeFaces=M; markv.assign(N,0); T0=chrono::steady_clock::now();
}

static int count_active_vertices(){
    vector<unsigned char> used(N,0); int c=0;
    for(int i=0;i<(int)F.size();i++) if(factive[i]){
        Face &f=F[i]; for(int k=0;k<3;k++){ int v=f.v[k]; if(v>=0&&v<N&&alive[v]&&!used[v]) used[v]=1,c++; }
    }
    return c;
}

static bool face_alive(int fid){ return fid>=0 && fid<(int)F.size() && factive[fid]; }
static Vec3 face_cross(const Face&f){ return cross3(P[f.v[1]]-P[f.v[0]], P[f.v[2]]-P[f.v[0]]); }
static Vec3 face_unit(const Face&f){ return unit(face_cross(f)); }

static void compact_inc(int v){
    if(v<0||v>=N) return; vector<int>& a=inc[v]; if(a.size()<96) return; int good=0; for(int fid:a) if(face_alive(fid)&&hasv(F[fid],v)) good++; if((int)a.size()<=good*3+32) return; vector<int> b; b.reserve(good+4); for(int fid:a) if(face_alive(fid)&&hasv(F[fid],v)) b.push_back(fid); a.swap(b);
}

static bool collect_edge(int u,int v,int ef[2],int opp[2]){
    if(u<0||v<0||u>=N||v>=N||!alive[u]||!alive[v]||u==v) return false;
    compact_inc(u); compact_inc(v);
    const vector<int> &lst = inc[u].size()<=inc[v].size()?inc[u]:inc[v];
    int cnt=0;
    for(int fid:lst){ if(!face_alive(fid)) continue; const Face&f=F[fid]; if(hasv(f,u)&&hasv(f,v)){ if(cnt>=2) return false; ef[cnt]=fid; opp[cnt]=thirdv(f,u,v); cnt++; } }
    return cnt==2 && opp[0]>=0 && opp[1]>=0 && opp[0]!=opp[1];
}

static bool link_ok(int u,int v,const int opp[2]){
    if(markToken>2000000000){ fill(markv.begin(),markv.end(),0); markToken=1; }
    int tok=markToken++;
    for(int fid:inc[u]) if(face_alive(fid)&&hasv(F[fid],u)){
        const Face&f=F[fid]; for(int k=0;k<3;k++){ int x=f.v[k]; if(x!=u&&x!=v) markv[x]=tok; }
    }
    int common=0; bool got0=false,got1=false;
    for(int fid:inc[v]) if(face_alive(fid)&&hasv(F[fid],v)){
        const Face&f=F[fid]; for(int k=0;k<3;k++){ int x=f.v[k]; if(x==u||x==v) continue; if(markv[x]==tok){ if(x!=opp[0]&&x!=opp[1]) return false; if(x==opp[0]) got0=true; if(x==opp[1]) got1=true; markv[x]=0; common++; } }
    }
    return common==2 && got0 && got1;
}

static bool same_tri(const Face&f,int a,int b,int c){ return tri_key(f.v[0],f.v[1],f.v[2])==tri_key(a,b,c); }
static bool duplicate_face_outside(int a,int b,int c,int rem,int ef0,int ef1,int skipfid){
    int ix=a; if(inc[b].size()<inc[ix].size()) ix=b; if(inc[c].size()<inc[ix].size()) ix=c;
    for(int fid:inc[ix]){
        if(!face_alive(fid)||fid==skipfid||fid==ef0||fid==ef1) continue;
        if(hasv(F[fid],rem)) continue;
        if(same_tri(F[fid],a,b,c)) return true;
    }
    return false;
}

struct Param{ double maxErr, plane, cosAng; bool mid=false; };
struct Eval{ bool ok=false; double cost=1e100, newErr=0; Vec3 pos; };

static Eval eval_dir(int keep,int rem,const int ef[2],const Param&p,const Vec3&pos){
    Eval r; r.pos=pos; r.newErr=max(errv[keep]+norm3(pos-P[keep]), errv[rem]+norm3(pos-P[rem]));
    if(r.newErr>p.maxErr) return r;
    if(r.newErr>tau*0.999) return r;
    vector<int> touched; touched.reserve(inc[keep].size()+inc[rem].size());
    for(int fid:inc[rem]) if(face_alive(fid)&&hasv(F[fid],rem)&&fid!=ef[0]&&fid!=ef[1]) touched.push_back(fid);
    if(p.mid){ for(int fid:inc[keep]) if(face_alive(fid)&&hasv(F[fid],keep)&&fid!=ef[0]&&fid!=ef[1]) touched.push_back(fid); sort(touched.begin(),touched.end()); touched.erase(unique(touched.begin(),touched.end()),touched.end()); }
    if(touched.empty()) return r;
    double maxPlane=0, maxFlip=0; int cnt=0;
    const double area_eps=max(1e-30,1e-24*diagLen*diagLen);
    for(int fid:touched){
        Face f=F[fid]; bool hasRem=hasv(f,rem), hasKeep=hasv(f,keep);
        if(hasRem&&hasKeep) continue;
        int nv[3]={f.v[0],f.v[1],f.v[2]}; for(int k=0;k<3;k++) if(nv[k]==rem) nv[k]=keep;
        if(nv[0]==nv[1]||nv[0]==nv[2]||nv[1]==nv[2]) return r;
        Vec3 oldc=face_cross(f); double oldl=norm3(oldc); if(!(oldl>area_eps)) return r;
        Vec3 a=(nv[0]==keep?pos:P[nv[0]]), b=(nv[1]==keep?pos:P[nv[1]]), c=(nv[2]==keep?pos:P[nv[2]]);
        Vec3 newc=cross3(b-a,c-a); double newl=norm3(newc); if(!(newl>area_eps)) return r;
        double nd=dot3(oldc,newc)/(oldl*newl); if(nd>1)nd=1; if(nd<-1)nd=-1; if(nd<p.cosAng) return r;
        Vec3 n=oldc*(1.0/oldl); double pl=fabs(dot3(n,pos-P[f.v[0]])); if(pl>p.plane) return r;
        if(hasRem && duplicate_face_outside(nv[0],nv[1],nv[2],rem,ef[0],ef[1],fid)) return r;
        maxPlane=max(maxPlane,pl); maxFlip=max(maxFlip,1.0-nd); cnt++;
    }
    if(cnt==0) return r;
    r.ok=true; r.cost=0.9*(r.newErr/(p.maxErr+1e-300))+0.7*(maxPlane/(p.plane+1e-300))+200.0*maxFlip+1e-5*cnt;
    return r;
}

static void apply_collapse(int keep,int rem,const int ef[2],const Vec3&pos,double newErr){
    for(int i=0;i<2;i++) if(face_alive(ef[i])){ factive[ef[i]]=0; F[ef[i]].active=0; activeFaces--; }
    vector<int> remList=inc[rem];
    for(int fid:remList){ if(!face_alive(fid)||!hasv(F[fid],rem)) continue; Face &f=F[fid]; bool hasKeep=hasv(f,keep); if(hasKeep){ factive[fid]=0; f.active=0; activeFaces--; continue; } for(int k=0;k<3;k++) if(f.v[k]==rem) f.v[k]=keep; inc[keep].push_back(fid); }
    alive[rem]=0; P[keep]=pos; errv[keep]=newErr; inc[rem].clear(); compact_inc(keep);
}

static bool try_collapse_edge(int u,int v,const Param&p){
    int ef[2],opp[2]; if(!collect_edge(u,v,ef,opp)) return false; if(!link_ok(u,v,opp)) return false;
    Eval a=eval_dir(u,v,ef,p,P[u]); Eval b=eval_dir(v,u,ef,p,P[v]);
    if(p.mid){ Vec3 m=(P[u]+P[v])*0.5; Eval am=eval_dir(u,v,ef,p,m); Eval bm=eval_dir(v,u,ef,p,m); if(am.ok&&(!a.ok||am.cost<a.cost)) a=am; if(bm.ok&&(!b.ok||bm.cost<b.cost)) b=bm; }
    if(!a.ok&&!b.ok) return false;
    if(b.ok&&(!a.ok||b.cost<a.cost)) apply_collapse(v,u,ef,b.pos,b.newErr); else apply_collapse(u,v,ef,a.pos,a.newErr);
    return true;
}

static int face_pass(const Param&p,double stopTime){
    int done=0;
    for(int fid=0;fid<(int)F.size();fid++){
        if((fid&4095)==0 && elapsed()>stopTime) break;
        if(!face_alive(fid)) continue; Face f=F[fid]; if(!alive[f.v[0]]||!alive[f.v[1]]||!alive[f.v[2]]) continue;
        struct E{ double l; int a,b; } e[3]={{norm2(P[f.v[0]]-P[f.v[1]]),f.v[0],f.v[1]},{norm2(P[f.v[1]]-P[f.v[2]]),f.v[1],f.v[2]},{norm2(P[f.v[2]]-P[f.v[0]]),f.v[2],f.v[0]}};
        sort(e,e+3,[](const E&a,const E&b){return a.l<b.l;});
        for(auto &x:e){ if(try_collapse_edge(x.a,x.b,p)){ done++; break; } }
    }
    return done;
}

static int sorted_pass(const Param&p,double stopTime,int maxTry=-1){
    vector<uint64_t> keys; keys.reserve((size_t)activeFaces*3);
    for(int fid=0;fid<(int)F.size();fid++){
        if((fid&8191)==0 && elapsed()>stopTime-0.2) return 0;
        if(!face_alive(fid)) continue; Face &f=F[fid]; if(!alive[f.v[0]]||!alive[f.v[1]]||!alive[f.v[2]]) continue;
        keys.push_back(ekey(f.v[0],f.v[1])); keys.push_back(ekey(f.v[1],f.v[2])); keys.push_back(ekey(f.v[2],f.v[0]));
    }
    sort(keys.begin(),keys.end()); keys.erase(unique(keys.begin(),keys.end()),keys.end());
    struct E{ double l; uint64_t k; }; vector<E> ed; ed.reserve(keys.size());
    for(uint64_t k:keys){ int a=k>>32,b=k&0xffffffffu; if(a>=0&&b>=0&&a<N&&b<N&&alive[a]&&alive[b]) ed.push_back({norm2(P[a]-P[b]),k}); }
    sort(ed.begin(),ed.end(),[](const E&a,const E&b){return a.l<b.l;});
    int done=0,tr=0;
    for(auto &e:ed){ if((tr++&1023)==0 && elapsed()>stopTime) break; if(maxTry>=0&&tr>maxTry) break; int a=e.k>>32,b=e.k&0xffffffffu; if(a>=0&&b>=0&&a<N&&b<N&&alive[a]&&alive[b]) if(try_collapse_edge(a,b,p)) done++; }
    return done;
}

struct Stats{ bool ok=false; double smooth=0, coarse=0, sharp=0, vsharp=0, bad=0; };
static Stats sample_stats(){
    Stats s; if(N<200||M<200) return s; int stride=max(1,M/40000), sampled=0,bad=0,sm=0,co=0,sh=0,vs=0; double c10=cos(10*M_PI/180.0), c30=cos(30*M_PI/180.0), c22=cos(22*M_PI/180.0), c45=cos(45*M_PI/180.0);
    for(int fid=0;fid<M&&sampled<120000;fid+=stride){ Face &f=F[fid]; int ed[3][2]={{f.v[0],f.v[1]},{f.v[1],f.v[2]},{f.v[2],f.v[0]}}; for(int k=0;k<3&&sampled<120000;k++){ int ef[2],op[2]; if(!collect_edge(ed[k][0],ed[k][1],ef,op)){bad++; continue;} Vec3 n0=face_unit(F[ef[0]]),n1=face_unit(F[ef[1]]); double d=clampd(dot3(n0,n1),-1,1); sampled++; if(d>c10)sm++; if(d>c30)co++; if(d<c22)sh++; if(d<c45)vs++; } }
    if(sampled<200) return s; s.ok=true; s.smooth=(double)sm/sampled; s.coarse=(double)co/sampled; s.sharp=(double)sh/sampled; s.vsharp=(double)vs/sampled; s.bad=(double)bad/(sampled+bad+1e-9); return s;
}

static bool strong_validator_current(){
    vector<int> id(N,-1); int vn=0; vector<Face> fs; fs.reserve(activeFaces); double eps=max(1e-30,1e-24*diagLen*diagLen);
    for(int fid=0;fid<(int)F.size();fid++) if(face_alive(fid)){
        Face f=F[fid]; for(int k=0;k<3;k++){ int v=f.v[k]; if(v<0||v>=N||!alive[v]) return false; if(id[v]<0) id[v]=vn++; f.v[k]=id[v]; }
        if(f.v[0]==f.v[1]||f.v[0]==f.v[2]||f.v[1]==f.v[2]) return false; Vec3 cr=cross3(P[F[fid].v[1]]-P[F[fid].v[0]],P[F[fid].v[2]]-P[F[fid].v[0]]); if(norm2(cr)<=eps) return false; fs.push_back(f);
    }
    if(vn<=0||fs.empty()||vn>N) return false;
    vector<uint64_t> edges; edges.reserve(fs.size()*3); vector<array<int,3>> tris; tris.reserve(fs.size());
    for(auto &f:fs){ tris.push_back(tri_key(f.v[0],f.v[1],f.v[2])); edges.push_back(ekey(f.v[0],f.v[1])); edges.push_back(ekey(f.v[1],f.v[2])); edges.push_back(ekey(f.v[2],f.v[0])); }
    sort(tris.begin(),tris.end()); if(adjacent_find(tris.begin(),tris.end())!=tris.end()) return false;
    sort(edges.begin(),edges.end()); for(size_t i=0;i<edges.size();){ size_t j=i+1; while(j<edges.size()&&edges[j]==edges[i]) j++; if(j-i!=2) return false; i=j; }
    return true;
}

static void restore_original(){
    P=Orig; F=OrigF; N=(int)P.size(); M=(int)F.size(); alive.assign(N,1); factive.assign(M,1); errv.assign(N,0.0); inc.assign(N,{}); vector<int>deg(N,0); for(auto &f:F){deg[f.v[0]]++;deg[f.v[1]]++;deg[f.v[2]]++;}
    for(int i=0;i<N;i++) inc[i].reserve(deg[i]+4); for(int i=0;i<M;i++){F[i].active=1;inc[F[i].v[0]].push_back(i);inc[F[i].v[1]].push_back(i);inc[F[i].v[2]].push_back(i);} activeFaces=M; markv.assign(N,0); markToken=1;
}

// -------------------------- rendering / proxy SSIM --------------------------
struct Maps{ int R=0; vector<float> depth; vector<float> norm; };
static inline void project_point(const Vec3&p,int view,int R,double&u,double&v,double&z){
    double f=800.0*(double)R/1024.0, c=0.5*R; double sx,sy; // views: +X,-X,+Y,-Y,+Z,-Z
    if(view==0){ sx=p.y; sy=p.z; z=2.5-p.x; }
    else if(view==1){ sx=-p.y; sy=p.z; z=2.5+p.x; }
    else if(view==2){ sx=-p.x; sy=p.z; z=2.5-p.y; }
    else if(view==3){ sx=p.x; sy=p.z; z=2.5+p.y; }
    else if(view==4){ sx=p.x; sy=p.y; z=2.5-p.z; }
    else { sx=-p.x; sy=p.y; z=2.5+p.z; }
    u=f*sx/z+c; v=f*sy/z+c;
}
static void render_mesh(const vector<Vec3>&V,const vector<Face>&FF,Maps&mp,int R){
    int pix=R*R; mp.R=R; mp.depth.assign((size_t)6*pix,255.f); mp.norm.assign((size_t)6*pix*3,127.5f);
    int nv=V.size(); vector<float> U(nv),VV(nv),Z(nv);
    vector<Vec3> fn(FF.size());
    for(size_t i=0;i<FF.size();i++){ const Face&f=FF[i]; Vec3 cr=cross3(V[f.v[1]]-V[f.v[0]],V[f.v[2]]-V[f.v[0]]); fn[i]=unit(cr); }
    for(int view=0;view<6;view++){
        for(int i=0;i<nv;i++){ double u,v,z; project_point(V[i],view,R,u,v,z); U[i]=u; VV[i]=v; Z[i]=z; }
        float *dptr=mp.depth.data()+(size_t)view*pix; float *nptr=mp.norm.data()+(size_t)view*pix*3;
        for(size_t ti=0;ti<FF.size();ti++){
            const Face&f=FF[ti]; int a=f.v[0],b=f.v[1],c=f.v[2]; if(a<0||b<0||c<0||a>=nv||b>=nv||c>=nv) continue;
            float x0=U[a],y0=VV[a],z0=Z[a], x1=U[b],y1=VV[b],z1=Z[b], x2=U[c],y2=VV[c],z2=Z[c]; if(!(z0>0&&z1>0&&z2>0)) continue;
            int xmin=max(0,(int)floor(min(x0,min(x1,x2))-.5f)), xmax=min(R-1,(int)ceil(max(x0,max(x1,x2))+.5f));
            int ymin=max(0,(int)floor(min(y0,min(y1,y2))-.5f)), ymax=min(R-1,(int)ceil(max(y0,max(y1,y2))+.5f)); if(xmin>xmax||ymin>ymax) continue;
            float den=(y1-y2)*(x0-x2)+(x2-x1)*(y0-y2); if(fabs(den)<1e-20f) continue; float inv=1.0f/den;
            Vec3 n=fn[ti]; float nr=(n.x+1)*127.5f, ng=(n.y+1)*127.5f, nb=(n.z+1)*127.5f;
            for(int yy=ymin;yy<=ymax;yy++){ float py=yy+0.5f; int row=yy*R; for(int xx=xmin;xx<=xmax;xx++){ float px=xx+0.5f; float w0=((y1-y2)*(px-x2)+(x2-x1)*(py-y2))*inv; float w1=((y2-y0)*(px-x2)+(x0-x2)*(py-y2))*inv; float w2=1-w0-w1; if(w0>=-1e-6f&&w1>=-1e-6f&&w2>=-1e-6f){ float dep=1.f/(w0/z0+w1/z1+w2/z2); int idx=row+xx; if(dep<dptr[idx]){ dptr[idx]=dep; float*q=nptr+idx*3; q[0]=nr;q[1]=ng;q[2]=nb; } } } }
        }
    }
}
static inline double rsum(const vector<double>&I,int W,int x0,int y0,int x1,int y1){ return I[(size_t)y1*W+x1]-I[(size_t)y0*W+x1]-I[(size_t)y1*W+x0]+I[(size_t)y0*W+x0]; }
static double ssim_channel(const float*X,int sx,const float*Y,int sy,const float*dX,const float*dY,int R,vector<double>&IX,vector<double>&IY,vector<double>&IX2,vector<double>&IY2,vector<double>&IXY){
    int W=R+1; size_t SZ=(size_t)W*W; fill(IX.begin(),IX.begin()+SZ,0); fill(IY.begin(),IY.begin()+SZ,0); fill(IX2.begin(),IX2.begin()+SZ,0); fill(IY2.begin(),IY2.begin()+SZ,0); fill(IXY.begin(),IXY.begin()+SZ,0);
    for(int y=1;y<=R;y++){ double ax=0,ay=0,ax2=0,ay2=0,axy=0; int row=(y-1)*R; for(int x=1;x<=R;x++){ int p=row+x-1; double xv=X[(size_t)p*sx], yv=Y[(size_t)p*sy]; ax+=xv; ay+=yv; ax2+=xv*xv; ay2+=yv*yv; axy+=xv*yv; size_t id=(size_t)y*W+x, up=(size_t)(y-1)*W+x; IX[id]=IX[up]+ax; IY[id]=IY[up]+ay; IX2[id]=IX2[up]+ax2; IY2[id]=IY2[up]+ay2; IXY[id]=IXY[up]+axy; } }
    const int rad=5, area=121; const double c1=6.5025,c2=58.5225; long long cnt=0; long double acc=0;
    for(int y=rad;y<R-rad;y++){ int row=y*R; for(int x=rad;x<R-rad;x++){ int p=row+x; if(!(dX[p]<254.f||dY[p]<254.f)) continue; int x0=x-rad,y0=y-rad,x1=x+rad+1,y1=y+rad+1; double Sx=rsum(IX,W,x0,y0,x1,y1), Sy=rsum(IY,W,x0,y0,x1,y1); double Sx2=rsum(IX2,W,x0,y0,x1,y1), Sy2=rsum(IY2,W,x0,y0,x1,y1), Sxy=rsum(IXY,W,x0,y0,x1,y1); double mux=Sx/area,muy=Sy/area; double vx=max(0.0,Sx2/area-mux*mux), vy=max(0.0,Sy2/area-muy*muy), cov=Sxy/area-mux*muy; double den=(mux*mux+muy*muy+c1)*(vx+vy+c2); double val=den!=0?((2*mux*muy+c1)*(2*cov+c2)/den):1.0; acc+=val; cnt++; } }
    return cnt?double(acc/cnt):1.0;
}
static double proxy_ssim(const vector<Vec3>&V,const vector<Face>&FF,int R){
    if(elapsed()>18.8) return -1;
    Maps A,B; render_mesh(Orig,OrigF,A,R); if(elapsed()>19.4) return -1; render_mesh(V,FF,B,R); int pix=R*R,W=R+1; vector<double>IX((size_t)W*W),IY((size_t)W*W),IX2((size_t)W*W),IY2((size_t)W*W),IXY((size_t)W*W); double total=0;
    for(int view=0;view<6;view++){ const float*od=A.depth.data()+(size_t)view*pix; const float*cd=B.depth.data()+(size_t)view*pix; double ns=0; for(int ch=0;ch<3;ch++){ const float*on=A.norm.data()+(size_t)view*pix*3+ch; const float*cn=B.norm.data()+(size_t)view*pix*3+ch; ns+=ssim_channel(on,3,cn,3,od,cd,R,IX,IY,IX2,IY2,IXY); } ns/=3.0; double ds=ssim_channel(od,1,cd,1,od,cd,R,IX,IY,IX2,IY2,IXY); total+=0.5*(ns+ds); }
    return total/6.0;
}

// -------------------------- candidate validation / hausdorff --------------------------
struct CellHash{ size_t operator()(const long long&x)const noexcept{ uint64_t z=(uint64_t)x; z^=z>>33; z*=0xff51afd7ed558ccdULL; z^=z>>33; z*=0xc4ceb9fe1a85ec53ULL; z^=z>>33; return (size_t)z; } };
static inline long long cell_key(int x,int y,int z){ return ((long long)(x+1048576)<<42) ^ ((long long)(y+1048576)<<21) ^ (long long)(z+1048576); }
static bool one_way_cover(const vector<Vec3>&A,const vector<Vec3>&B,double r){
    if(A.empty()||B.empty()) return false; double cell=max(r,1e-9), r2=r*r; unordered_map<long long, vector<int>, CellHash> mp; mp.reserve(B.size()*2+10);
    Vec3 mn=bbMin; for(size_t i=0;i<B.size();i++){ int ix=floor((B[i].x-mn.x)/cell), iy=floor((B[i].y-mn.y)/cell), iz=floor((B[i].z-mn.z)/cell); mp[cell_key(ix,iy,iz)].push_back((int)i); }
    for(size_t qi=0;qi<A.size();qi++){ const Vec3&q=A[qi]; int ix=floor((q.x-mn.x)/cell), iy=floor((q.y-mn.y)/cell), iz=floor((q.z-mn.z)/cell); bool ok=false; for(int dx=-1;dx<=1&&!ok;dx++) for(int dy=-1;dy<=1&&!ok;dy++) for(int dz=-1;dz<=1&&!ok;dz++){ auto it=mp.find(cell_key(ix+dx,iy+dy,iz+dz)); if(it==mp.end()) continue; for(int j:it->second) if(norm2(q-B[j])<=r2){ ok=true; break; } } if(!ok) return false; }
    return true;
}
static bool vertex_hausdorff_ok(const vector<Vec3>&V,double r){ return one_way_cover(Orig,V,r) && one_way_cover(V,Orig,r); }
static bool valid_mesh_basic(const vector<Vec3>&V,const vector<Face>&FF){
    if(V.empty()||FF.empty()||V.size()>(size_t)N) return false; double eps=max(1e-30,1e-24*diagLen*diagLen); vector<uint64_t> edges; edges.reserve(FF.size()*3); vector<array<int,3>> tris; tris.reserve(FF.size());
    for(const Face&f:FF){ int a=f.v[0],b=f.v[1],c=f.v[2]; if(a<0||b<0||c<0||a>=(int)V.size()||b>=(int)V.size()||c>=(int)V.size()) return false; if(a==b||a==c||b==c) return false; if(norm2(cross3(V[b]-V[a],V[c]-V[a]))<=eps) return false; edges.push_back(ekey(a,b)); edges.push_back(ekey(b,c)); edges.push_back(ekey(c,a)); tris.push_back(tri_key(a,b,c)); }
    sort(tris.begin(),tris.end()); if(adjacent_find(tris.begin(),tris.end())!=tris.end()) return false; sort(edges.begin(),edges.end()); for(size_t i=0;i<edges.size();){ size_t j=i+1; while(j<edges.size()&&edges[j]==edges[i]) j++; if(j-i!=2) return false; i=j; } return true;
}
static void orient_face(vector<Vec3>&V, Face&f, const Vec3&center){ Vec3 cr=cross3(V[f.v[1]]-V[f.v[0]],V[f.v[2]]-V[f.v[0]]); Vec3 ce=(V[f.v[0]]+V[f.v[1]]+V[f.v[2]])/3.0; if(dot3(cr,ce-center)<0) swap(f.v[1],f.v[2]); }
static bool accept_candidate(const Mesh&mesh,double ssimNeed=0.915,int R=256){
    if(mesh.V.empty()||mesh.F.empty()) return false; int cur=count_active_vertices(); if((int)mesh.V.size()>=cur) return false; if(!valid_mesh_basic(mesh.V,mesh.F)) return false; if(!vertex_hausdorff_ok(mesh.V,tau*0.998)) return false;
    if(elapsed()<15.0){ double sc=proxy_ssim(mesh.V,mesh.F,R); if(sc>=0 && sc<ssimNeed) return false; }
    // Load candidate into current state and finish.
    int n=(int)mesh.V.size(), m=(int)mesh.F.size(); N=n; M=m; P=mesh.V; Orig=P; F=mesh.F; OrigF=F; alive.assign(N,1); factive.assign(M,1); errv.assign(N,0); inc.assign(N,{}); vector<int>deg(N,0); for(auto &f:F){deg[f.v[0]]++;deg[f.v[1]]++;deg[f.v[2]]++;}
    for(int i=0;i<N;i++) inc[i].reserve(deg[i]); for(int i=0;i<M;i++){F[i].active=1; inc[F[i].v[0]].push_back(i); inc[F[i].v[1]].push_back(i); inc[F[i].v[2]].push_back(i);} activeFaces=M; return true;
}

// -------------------------- analytic / structured replacements --------------------------
static Mesh make_uv_sphere(Vec3 c,double rx,double ry,double rz,int lat,int lon){
    Mesh m; lat=max(lat,4); lon=max(lon,8); m.V.reserve(2+(lat-1)*lon); m.F.reserve(2*lat*lon); m.V.push_back({c.x,c.y,c.z+rz}); for(int i=1;i<lat;i++){ double th=M_PI*i/lat; double st=sin(th), ct=cos(th); for(int j=0;j<lon;j++){ double ph=2*M_PI*j/lon; m.V.push_back({c.x+rx*st*cos(ph), c.y+ry*st*sin(ph), c.z+rz*ct}); } } int bot=m.V.size(); m.V.push_back({c.x,c.y,c.z-rz}); auto id=[&](int r,int j){ return 1+(r-1)*lon+((j%lon+lon)%lon); }; auto add=[&](int a,int b,int cc){ Face f{{a,b,cc},1}; orient_face(m.V,f,c); m.F.push_back(f);}; for(int j=0;j<lon;j++) add(0,id(1,j+1),id(1,j)); for(int r=1;r<lat-1;r++) for(int j=0;j<lon;j++){ int a=id(r,j),b=id(r,j+1),cc=id(r+1,j),d=id(r+1,j+1); add(a,b,cc); add(b,d,cc); } for(int j=0;j<lon;j++) add(bot,id(lat-1,j),id(lat-1,j+1)); return m;
}
static bool try_ellipsoid(){
    if(N<700||elapsed()>3.0) return false; Vec3 c=(bbMin+bbMax)*0.5; double rx=(bbMax.x-bbMin.x)*0.5, ry=(bbMax.y-bbMin.y)*0.5, rz=(bbMax.z-bbMin.z)*0.5; double rmax=max(rx,max(ry,rz)), rmin=min(rx,min(ry,rz)); if(!(rmin>1e-12)||rmin<0.35*rmax) return false;
    int stride=max(1,N/240000); double sum=0,sq=0,mx=0; int cnt=0; for(int i=0;i<N;i+=stride){ Vec3 q=Orig[i]-c; double rr=sqrt((q.x*q.x)/(rx*rx)+(q.y*q.y)/(ry*ry)+(q.z*q.z)/(rz*rz)); double e=fabs(rr-1); sum+=e; sq+=e*e; mx=max(mx,e); cnt++; }
    double rms=sqrt(sq/max(1,cnt)); if(!(rms<0.010 && mx<0.055)) return false;
    double step=max(0.045*diagLen, tau*0.82); int lat=(int)ceil(M_PI*rmax/step); lat=max(12,min(56,lat)); int lon=max(24,min(112,2*lat)); Mesh m=make_uv_sphere(c,rx,ry,rz,lat,lon); int R=N<30000?512:(N<120000?256:128); double need=N<5000?0.925:0.915; return accept_candidate(m,need,R);
}

static Mesh make_torus(int axis,double ct,double cu,double cv,double Rmaj,double rminr,int U,int Vv){
    Mesh m; U=max(U,12); Vv=max(Vv,6); m.V.reserve(U*Vv); m.F.reserve(2*U*Vv); auto make=[&](double t,double u,double v){ if(axis==0) return Vec3{t,u,v}; if(axis==1) return Vec3{u,t,v}; return Vec3{u,v,t}; };
    for(int i=0;i<U;i++){ double th=2*M_PI*i/U, co=cos(th), si=sin(th); for(int j=0;j<Vv;j++){ double ph=2*M_PI*j/Vv, cp=cos(ph), sp=sin(ph); double rr=Rmaj+rminr*cp; m.V.push_back(make(ct+rminr*sp, cu+rr*co, cv+rr*si)); } }
    Vec3 center=make(ct,cu,cv); auto id=[&](int i,int j){ return ((i%U+U)%U)*Vv+((j%Vv+Vv)%Vv);}; auto add=[&](int a,int b,int c){ Face f{{a,b,c},1}; orient_face(m.V,f,center); m.F.push_back(f);}; for(int i=0;i<U;i++) for(int j=0;j<Vv;j++){ int a=id(i,j),b=id(i+1,j),c=id(i+1,j+1),d=id(i,j+1); add(a,b,c); add(a,c,d);} return m;
}
static bool try_torus(){
    if(N<600||N>80000||elapsed()>4.5) return false; for(int ax=0;ax<3;ax++){
        double min_t=1e100,max_t=-1e100,min_u=1e100,max_u=-1e100,min_v=1e100,max_v=-1e100; auto coord=[&](const Vec3&p,double&t,double&u,double&v){ if(ax==0){t=p.x;u=p.y;v=p.z;} else if(ax==1){t=p.y;u=p.x;v=p.z;} else {t=p.z;u=p.x;v=p.y;} };
        for(auto&p:Orig){ double t,u,v; coord(p,t,u,v); min_t=min(min_t,t); max_t=max(max_t,t); min_u=min(min_u,u); max_u=max(max_u,u); min_v=min(min_v,v); max_v=max(max_v,v); }
        double ct=(min_t+max_t)*.5, cu=(min_u+max_u)*.5, cv=(min_v+max_v)*.5; double minr=1e100,maxr=0; for(auto&p:Orig){ double t,u,v; coord(p,t,u,v); double r=hypot(u-cu,v-cv); minr=min(minr,r); maxr=max(maxr,r);} if(!(maxr>minr&&minr>1e-8)) continue; double Rmaj=(maxr+minr)*.5, rt=(max_t-min_t)*.5, rr=(maxr-minr)*.5, r=(rt+rr)*.5; if(!(Rmaj>1.25*r&&r>1e-8)) continue; if(fabs(rt-rr)>0.28*r) continue; int stride=max(1,N/220000); double sq=0,mx=0; int cnt=0; for(int i=0;i<N;i+=stride){ double t,u,v; coord(Orig[i],t,u,v); double rho=hypot(u-cu,v-cv); double e=fabs(hypot(rho-Rmaj,t-ct)-r)/r; sq+=e*e; mx=max(mx,e); cnt++; } double rms=sqrt(sq/max(1,cnt)); if(rms>0.035||mx>0.12) continue; double step=tau*0.78; int U=max(24,min(144,(int)ceil(2*M_PI*Rmaj/step))); int Vv=max(10,min(64,(int)ceil(2*M_PI*r/step))); Mesh m=make_torus(ax,ct,cu,cv,Rmaj,r,U,Vv); int Rres=N<12000?512:256; if(accept_candidate(m,0.915,Rres)) return true; }
    return false;
}

static Mesh make_frustum(int axis,double t0,double t1,double cu,double cv,double r0,double r1,int rings,int sides){
    Mesh m; rings=max(rings,2); sides=max(sides,8); auto make=[&](double t,double u,double v){ if(axis==0) return Vec3{t,u,v}; if(axis==1) return Vec3{u,t,v}; return Vec3{u,v,t}; }; Vec3 ctr=make((t0+t1)*.5,cu,cv); vector<int> ringStart; for(int i=0;i<rings;i++){ double s=(double)i/(rings-1), t=t0+(t1-t0)*s, r=r0+(r1-r0)*s; if(r<tau*0.15){ ringStart.push_back((int)m.V.size()); m.V.push_back(make(t,cu,cv)); } else { ringStart.push_back((int)m.V.size()); for(int j=0;j<sides;j++){ double ph=2*M_PI*j/sides; m.V.push_back(make(t,cu+r*cos(ph),cv+r*sin(ph))); } } }
    auto vid=[&](int i,int j){ double s=(double)i/(rings-1); double r=r0+(r1-r0)*s; int st=ringStart[i]; return r<tau*.15?st:st+((j%sides+sides)%sides); };
    auto add=[&](int a,int b,int c){ if(a==b||a==c||b==c) return; Face f{{a,b,c},1}; orient_face(m.V,f,ctr); m.F.push_back(f); };
    for(int i=0;i+1<rings;i++){ bool aA=(vid(i,0)==vid(i,1)), bA=(vid(i+1,0)==vid(i+1,1)); for(int j=0;j<sides;j++){ int a=vid(i,j),b=vid(i,j+1),c=vid(i+1,j),d=vid(i+1,j+1); if(aA&&!bA) add(a,d,c); else if(!aA&&bA) add(c,a,b); else { add(a,b,c); add(b,d,c); } } }
    // caps when ring is not apex
    for(int end=0;end<2;end++){ int i=end?rings-1:0; if(vid(i,0)==vid(i,1)) continue; int cen=m.V.size(); double t=end?t1:t0; m.V.push_back(make(t,cu,cv)); for(int j=0;j<sides;j++){ if(end) add(cen,vid(i,j+1),vid(i,j)); else add(cen,vid(i,j),vid(i,j+1)); } }
    return m;
}
static bool try_frustum(){
    if(N<800||N>160000||elapsed()>5.5) return false; for(int ax=0;ax<3;ax++){
        double mt=1e100,Mt=-1e100,mu=1e100,Mu=-1e100,mv=1e100,Mv=-1e100; auto coord=[&](const Vec3&p,double&t,double&u,double&v){ if(ax==0){t=p.x;u=p.y;v=p.z;} else if(ax==1){t=p.y;u=p.x;v=p.z;} else {t=p.z;u=p.x;v=p.y;} };
        for(auto&p:Orig){ double t,u,v; coord(p,t,u,v); mt=min(mt,t);Mt=max(Mt,t);mu=min(mu,u);Mu=max(Mu,u);mv=min(mv,v);Mv=max(Mv,v);} double len=Mt-mt, cu=(mu+Mu)/2, cv=(mv+Mv)/2; if(!(len>1e-6)) continue; int stride=max(1,N/220000); double S=0,St=0,Stt=0,Sr=0,Str=0,maxr=0; int cnt=0,axisPts=0; for(int i=0;i<N;i+=stride){ double t,u,v; coord(Orig[i],t,u,v); double r=hypot(u-cu,v-cv); maxr=max(maxr,r); if(r<0.04*max(Mu-mu,Mv-mv)){ axisPts++; continue; } S++; St+=t; Stt+=t*t; Sr+=r; Str+=t*r; cnt++; } if(cnt<200||!(maxr>1e-8)) continue; double den=S*Stt-St*St; if(fabs(den)<1e-18) continue; double a=(S*Str-St*Sr)/den, b=(Sr-a*St)/S; double r0=max(0.0,a*mt+b), r1=max(0.0,a*Mt+b); if(max(r0,r1)<0.25*maxr) continue; double sq=0,mx=0; int ck=0; for(int i=0;i<N;i+=stride){ double t,u,v; coord(Orig[i],t,u,v); double r=hypot(u-cu,v-cv), pred=max(0.0,a*t+b); if(r<0.05*maxr) continue; double e=fabs(r-pred)/maxr; sq+=e*e; mx=max(mx,e); ck++; } if(ck<200) continue; double rms=sqrt(sq/ck); if(rms>0.025||mx>0.10) continue; int sides=max(16,min(128,(int)ceil(2*M_PI*max(r0,r1)/(tau*.82)))); int rings=max(2,min(80,(int)ceil(len/(tau*.82))+1)); Mesh m=make_frustum(ax,mt,Mt,cu,cv,r0,r1,rings,sides); int R=N<20000?512:256; if(accept_candidate(m,0.925,R)) return true; }
    return false;
}

// Box boundary grid. Useful for subdivided boxes; Hausdorff verifier rejects non-boxes.

static Mesh make_capsule(int axis,double ct,double cu,double cv,double halfCyl,double r,int sides,int capSeg){
    Mesh m; sides=max(sides,8); capSeg=max(capSeg,3); auto make=[&](double t,double u,double v){ if(axis==0) return Vec3{t,u,v}; if(axis==1) return Vec3{u,t,v}; return Vec3{u,v,t}; };
    Vec3 ctr=make(ct,cu,cv); vector<int> rings; rings.reserve(2*capSeg+2);
    int topPole=m.V.size(); m.V.push_back(make(ct+halfCyl+r,cu,cv));
    auto add_ring=[&](double t,double rr){ rings.push_back((int)m.V.size()); for(int j=0;j<sides;j++){ double ph=2*M_PI*j/sides; m.V.push_back(make(t,cu+rr*cos(ph),cv+rr*sin(ph))); } };
    for(int i=1;i<=capSeg;i++){ double phi=(M_PI/2.0)*(double)i/capSeg; add_ring(ct+halfCyl+r*cos(phi), r*sin(phi)); }
    if(halfCyl>tau*.25) add_ring(ct-halfCyl, r);
    for(int i=capSeg-1;i>=1;i--){ double phi=(M_PI/2.0)*(double)i/capSeg; add_ring(ct-halfCyl-r*cos(phi), r*sin(phi)); }
    int botPole=m.V.size(); m.V.push_back(make(ct-halfCyl-r,cu,cv));
    auto rid=[&](int ri,int j){ return rings[ri]+((j%sides+sides)%sides); };
    auto add=[&](int a,int b,int c){ Face f{{a,b,c},1}; orient_face(m.V,f,ctr); m.F.push_back(f); };
    for(int j=0;j<sides;j++) add(topPole,rid(0,j+1),rid(0,j));
    for(int ri=0;ri+1<(int)rings.size();ri++) for(int j=0;j<sides;j++){ add(rid(ri,j),rid(ri+1,j),rid(ri+1,j+1)); add(rid(ri,j),rid(ri+1,j+1),rid(ri,j+1)); }
    int last=(int)rings.size()-1; for(int j=0;j<sides;j++) add(botPole,rid(last,j),rid(last,j+1));
    return m;
}
static bool try_capsule(){
    if(N<1000||N>300000||elapsed()>5.0) return false;
    for(int ax=0;ax<3;ax++){
        double mt=1e100,Mt=-1e100,mu=1e100,Mu=-1e100,mv=1e100,Mv=-1e100; auto coord=[&](const Vec3&p,double&t,double&u,double&v){ if(ax==0){t=p.x;u=p.y;v=p.z;} else if(ax==1){t=p.y;u=p.x;v=p.z;} else {t=p.z;u=p.x;v=p.y;} };
        for(const Vec3&p:Orig){ double t,u,v; coord(p,t,u,v); mt=min(mt,t); Mt=max(Mt,t); mu=min(mu,u); Mu=max(Mu,u); mv=min(mv,v); Mv=max(Mv,v); }
        double extT=Mt-mt, eu=Mu-mu, ev=Mv-mv; if(!(extT>1e-9&&eu>1e-9&&ev>1e-9)) continue; if(min(eu,ev)<max(eu,ev)*0.96) continue;
        double r=0.25*(eu+ev); double halfCyl=0.5*extT-r; if(!(r>1e-9&&halfCyl>0.12*r)) continue; double ct=(mt+Mt)*.5, cu=(mu+Mu)*.5, cv=(mv+Mv)*.5;
        int stride=max(1,N/240000), cnt=0, side=0, cap=0, bad=0; double sq=0,mx=0;
        for(int i=0;i<N;i+=stride){ double t,u,v; coord(Orig[i],t,u,v); double dt=t-ct, du=u-cu, dv=v-cv; double rho=hypot(du,dv); double e; if(fabs(dt)<=halfCyl){ e=fabs(rho-r)/r; side++; } else { double q=fabs(dt)-halfCyl; e=fabs(hypot(rho,q)-r)/r; cap++; } sq+=e*e; mx=max(mx,e); if(e>0.075) bad++; cnt++; }
        if(cnt<300||side<cnt/8||cap<cnt/12) continue; double rms=sqrt(sq/cnt); if(rms>0.018||mx>0.095||bad>max(4,cnt/50)) continue;
        int sides=max(20,min(96,(int)ceil(2*M_PI*r/(tau*.82)))); int capSeg=max(5,min(20,(int)ceil((M_PI/2*r)/(tau*.82)))); Mesh m=make_capsule(ax,ct,cu,cv,halfCyl,r,sides,capSeg); int R=N<30000?512:256; if(accept_candidate(m,0.925,R)) return true;
    }
    return false;
}

static bool try_box_grid(){
    if(N<200||elapsed()>6.0) return false; double ex=bbMax.x-bbMin.x,ey=bbMax.y-bbMin.y,ez=bbMax.z-bbMin.z; double mx=max(ex,max(ey,ez)); if(!(mx>1e-9)) return false; int stride=max(1,N/200000); int near=0,cnt=0; double eps=max(tau*.18,mx*1e-5); for(int i=0;i<N;i+=stride){ Vec3 p=Orig[i]; double d=min({fabs(p.x-bbMin.x),fabs(p.x-bbMax.x),fabs(p.y-bbMin.y),fabs(p.y-bbMax.y),fabs(p.z-bbMin.z),fabs(p.z-bbMax.z)}); if(d<=eps) near++; cnt++; } if(near*100<cnt*98) return false;
    int nx=max(1,min(80,(int)ceil(ex/(tau*.82)))), ny=max(1,min(80,(int)ceil(ey/(tau*.82)))), nz=max(1,min(80,(int)ceil(ez/(tau*.82))));
    Mesh m; unordered_map<long long,int,CellHash> id; id.reserve((nx+1)*(ny+1)*2+(nx+1)*(nz+1)*2+(ny+1)*(nz+1)*2);
    auto key=[&](int i,int j,int k){ return ((long long)i<<42)^((long long)j<<21)^k; };
    auto get=[&](int i,int j,int k)->int{ long long kk=key(i,j,k); auto it=id.find(kk); if(it!=id.end()) return it->second; double x=bbMin.x+ex*(double)i/nx, y=bbMin.y+ey*(double)j/ny, z=bbMin.z+ez*(double)k/nz; int idx=m.V.size(); id[kk]=idx; m.V.push_back({x,y,z}); return idx; };
    Vec3 ctr=(bbMin+bbMax)*.5; auto add=[&](int a,int b,int c){ Face f{{a,b,c},1}; orient_face(m.V,f,ctr); m.F.push_back(f); };
    auto quad=[&](int a,int b,int c,int d){ add(a,b,c); add(a,c,d); };
    for(int i=0;i<nx;i++) for(int j=0;j<ny;j++){ quad(get(i,j,0),get(i+1,j,0),get(i+1,j+1,0),get(i,j+1,0)); quad(get(i,j,nz),get(i,j+1,nz),get(i+1,j+1,nz),get(i+1,j,nz)); }
    for(int i=0;i<nx;i++) for(int k=0;k<nz;k++){ quad(get(i,0,k),get(i,0,k+1),get(i+1,0,k+1),get(i+1,0,k)); quad(get(i,ny,k),get(i+1,ny,k),get(i+1,ny,k+1),get(i,ny,k+1)); }
    for(int j=0;j<ny;j++) for(int k=0;k<nz;k++){ quad(get(0,j,k),get(0,j+1,k),get(0,j+1,k+1),get(0,j,k+1)); quad(get(nx,j,k),get(nx,j,k+1),get(nx,j+1,k+1),get(nx,j+1,k)); }
    int R=N<20000?512:256; return accept_candidate(m,0.93,R);
}

static bool detect_torus_grid(int&U,int&Vv){
    if(M!=2*N||N<300) return false; vector<pair<int,int>> cand; for(int d=6;(long long)d*d<=N;d++) if(N%d==0){ cand.push_back({N/d,d}); cand.push_back({d,N/d}); }
    auto same=[&](const Face&f,int a,int b,int c){ return tri_key(f.v[0],f.v[1],f.v[2])==tri_key(a,b,c); };
    for(auto pr:cand){ U=pr.first; Vv=pr.second; if(U<6||Vv<6||U>20000||Vv>20000) continue; int ok=0,tot=0,st=max(1,N/700); for(int q=0;q<N&&tot<700;q+=st){ int i=q/Vv,j=q-i*Vv; int ni=(i+1)%U,nj=(j+1)%Vv; int a=i*Vv+j,b=ni*Vv+j,c=ni*Vv+nj,d=i*Vv+nj; int fid=2*q; if(fid+1>=M) break; if((same(OrigF[fid],a,b,c)&&same(OrigF[fid+1],a,c,d))||(same(OrigF[fid],a,b,d)&&same(OrigF[fid+1],b,c,d))) ok++; tot++; } if(tot>100&&ok*100>=tot*98) return true; }
    return false;
}
static bool try_grid_downsample(){
    if(elapsed()>7.0) return false; int U=0,Vv=0; if(!detect_torus_grid(U,Vv)) return false; int cur=N; vector<pair<int,int>> trials; for(double fct:{10.0,8.0,6.0,5.0,4.0,3.0,2.0}){ int u=max(6,(int)round(U/fct)), v=max(6,(int)round(Vv/fct)); if(u*v<cur) trials.push_back({u,v}); } sort(trials.begin(),trials.end(),[](auto a,auto b){return a.first*a.second<b.first*b.second;}); trials.erase(unique(trials.begin(),trials.end()),trials.end());
    for(auto [U2,V2]:trials){ Mesh m; m.V.reserve(U2*V2); m.F.reserve(2*U2*V2); for(int i=0;i<U2;i++){ int oi=(long long)i*U/U2; for(int j=0;j<V2;j++){ int oj=(long long)j*Vv/V2; m.V.push_back(Orig[oi*Vv+oj]); } } auto id=[&](int i,int j){return ((i%U2+U2)%U2)*V2+((j%V2+V2)%V2);}; for(int i=0;i<U2;i++)for(int j=0;j<V2;j++){ Face a{{id(i,j),id(i+1,j),id(i+1,j+1)},1}, b{{id(i,j),id(i+1,j+1),id(i,j+1)},1}; m.F.push_back(a); m.F.push_back(b);} if(accept_candidate(m,0.92,N<30000?512:256)) return true; }
    return false;
}

static bool detect_latlong(int&R,int&Vv){
    if(N<300||M!=2*(N-2)) return false; auto same=[&](int fid,int a,int b,int c){ return tri_key(OrigF[fid].v[0],OrigF[fid].v[1],OrigF[fid].v[2])==tri_key(a,b,c); };
    for(int vv=8;vv<=4096;vv++){ if((N-2)%vv) continue; int r=(N-2)/vv; if(r<3) continue; bool ok=true; int st=max(1,vv/64); for(int j=0;j<vv&&ok;j+=st){ int a=1+j,b=1+(j+1)%vv; if(!same(j,0,b,a)&&!same(j,0,a,b)) ok=false; int off=vv+2*(r-1)*vv+j; int c=1+(r-1)*vv+j,d=1+(r-1)*vv+(j+1)%vv; if(ok&&!same(off,N-1,c,d)&&!same(off,N-1,d,c)) ok=false; } if(ok){R=r;Vv=vv;return true;} }
    return false;
}
static bool try_latlong_downsample(){
    if(elapsed()>7.5) return false; int R0=0,V0=0; if(!detect_latlong(R0,V0)) return false; vector<pair<int,int>> trials; for(double f:{10,8,6,5,4,3,2}){ int r=max(3,(int)round(R0/f)), v=max(8,(int)round(V0/f)); if(2+r*v<N) trials.push_back({r,v}); } sort(trials.begin(),trials.end(),[](auto a,auto b){return 2+a.first*a.second<2+b.first*b.second;});
    for(auto [R2,V2]:trials){ Mesh m; m.V.reserve(2+R2*V2); m.F.reserve(2*R2*V2); m.V.push_back(Orig[0]); for(int i=0;i<R2;i++){ int oi=1+(long long)i*(R0-1)/max(1,R2-1); for(int j=0;j<V2;j++){ int oj=(long long)j*V0/V2; m.V.push_back(Orig[1+(oi-1)*V0+oj]); } } int bot=m.V.size(); m.V.push_back(Orig[N-1]); auto id=[&](int r,int j){return 1+r*V2+((j%V2+V2)%V2);}; Vec3 ctr=(bbMin+bbMax)*.5; auto add=[&](int a,int b,int c){ Face f{{a,b,c},1}; orient_face(m.V,f,ctr); m.F.push_back(f);}; for(int j=0;j<V2;j++) add(0,id(0,j+1),id(0,j)); for(int r=0;r+1<R2;r++)for(int j=0;j<V2;j++){ add(id(r,j),id(r+1,j),id(r+1,j+1)); add(id(r,j),id(r+1,j+1),id(r,j+1)); } for(int j=0;j<V2;j++) add(bot,id(R2-1,j),id(R2-1,j+1)); if(accept_candidate(m,0.92,N<30000?512:256)) return true; }
    return false;
}

static bool try_special_replacements(){
    // Cheap structured topology first, then analytic primitives. Each candidate is guarded by exact vertex-Hausdorff and proxy SSIM.
    if(try_latlong_downsample()) return true;
    if(try_grid_downsample()) return true;
    if(try_ellipsoid()) return true;
    if(try_torus()) return true;
    if(try_frustum()) return true;
    if(try_capsule()) return true;
    if(try_box_grid()) return true;
    return false;
}

static void simplify_edges(){
    Stats st=sample_stats(); double d=diagLen; vector<Param> passes;
    auto add=[&](double a,double pl,double deg,bool mid=false){ passes.push_back({min(a*d,tau*.985),pl*d,cos(deg*M_PI/180.0),mid}); };
    add(.004,.001,1.0); add(.008,.002,1.8); add(.012,.0032,2.8); if(N<1000) add(.018,.005,4.0); else if(N<30000) add(.022,.0065,4.8); else add(.026,.008,5.8);
    for(auto&p:passes){ if(elapsed()>15.5) break; int got=face_pass(p,16.0); if(got==0 && p.maxErr>.015*d) break; }
    bool smooth=st.ok && st.coarse>.965 && st.vsharp<.025 && st.bad<.01;
    bool medium=st.ok && st.coarse>.86 && st.sharp>.02 && st.sharp<.18 && st.vsharp<.08 && st.bad<.012;
    if(elapsed()<16.4 && smooth){ Param p{min(.042*d,tau*.992),.015*d,cos(10.5*M_PI/180.0),true}; sorted_pass(p,18.0); }
    if(elapsed()<16.8 && medium){ Param p{min(.034*d,tau*.99),.010*d,cos(7.8*M_PI/180.0),false}; sorted_pass(p,18.1, max(40000,N/2)); }
    if(elapsed()<17.4 && N>=30000){ Param p{min(.030*d,tau*.985),.009*d,cos(6.5*M_PI/180.0),false}; sorted_pass(p,18.3, max(50000,N/3)); }
    // final tiny planar vertex removal: a safe extra pass on very flat local regions
    if(elapsed()<18.2){ Param p{min(.049*d,tau*.995),.020*d,cos(14.0*M_PI/180.0),true}; sorted_pass(p,18.8, N<50000?80000:160000); }
}

static void write_current(){
    if(!strong_validator_current()) restore_original();
    vector<int> remap(N,-1); vector<int> ids; ids.reserve(N); for(int fid=0;fid<(int)F.size();fid++) if(face_alive(fid)){ Face &f=F[fid]; for(int k=0;k<3;k++){ int v=f.v[k]; if(v>=0&&v<N&&alive[v]&&remap[v]<0){ remap[v]=ids.size(); ids.push_back(v); } } }
    vector<Face> outF; outF.reserve(activeFaces); for(int fid=0;fid<(int)F.size();fid++) if(face_alive(fid)){ Face f=F[fid]; int a=remap[f.v[0]],b=remap[f.v[1]],c=remap[f.v[2]]; if(a>=0&&b>=0&&c>=0&&a!=b&&a!=c&&b!=c){ Face g{{a,b,c},1}; outF.push_back(g); } }
    static char obuf[1<<20]; setvbuf(stdout,obuf,_IOFBF,sizeof(obuf)); printf("%d %d\n",(int)ids.size(),(int)outF.size()); bool highprec=(int)ids.size()*2<= (int)Orig.size();
    for(int old:ids){ const Vec3&p=P[old]; if(highprec) printf("v %.15g %.15g %.15g\n",p.x,p.y,p.z); else printf("v %.10g %.10g %.10g\n",p.x,p.y,p.z); }
    for(auto &f:outF) printf("f %d %d %d\n",f.v[0]+1,f.v[1]+1,f.v[2]+1);
}

int main(){
    read_input();
    // Candidate replacements are intentionally time-capped and fail closed.
    if(!try_special_replacements()){
        simplify_edges();
    }
    write_current();
    return 0;
}
