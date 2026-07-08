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
static inline Vec3 normalize(Vec3 a){ double l=norm3(a); return l>1e-300? a/l : Vec3{0,0,0}; }

struct FastInput{
    vector<char> buf; char *p=nullptr;
    FastInput(){
        buf.reserve(1<<27); char tmp[1<<16]; size_t n;
        while((n=fread(tmp,1,sizeof(tmp),stdin))>0) buf.insert(buf.end(),tmp,tmp+n);
        buf.push_back(0); p=buf.data();
    }
    inline void ws(){ while(*p==' '||*p=='\n'||*p=='\r'||*p=='\t') ++p; }
    inline char ch(){ ws(); return *p++; }
    inline int nextInt(){ ws(); int s=1; if(*p=='-') s=-1,++p; int x=0; while(*p>='0'&&*p<='9') x=x*10+(*p++-'0'); return x*s; }
    inline double nextDouble(){ ws(); char* e; double x=strtod(p,&e); p=e; return x; }
};

struct Face{ int v[3]; unsigned char active=1; Vec3 n{0,0,0}; };
struct Tri{ int a,b,c; };
struct Quadric{
    double q[10];
    Quadric(){ memset(q,0,sizeof(q)); }
    inline void addPlane(double a,double b,double c,double d,double w){
        q[0]+=w*a*a; q[1]+=w*a*b; q[2]+=w*a*c; q[3]+=w*a*d;
        q[4]+=w*b*b; q[5]+=w*b*c; q[6]+=w*b*d;
        q[7]+=w*c*c; q[8]+=w*c*d; q[9]+=w*d*d;
    }
    inline void add(const Quadric&o){ for(int i=0;i<10;i++) q[i]+=o.q[i]; }
    inline double eval(const Vec3&p)const{
        double x=p.x,y=p.y,z=p.z;
        return q[0]*x*x + 2*q[1]*x*y + 2*q[2]*x*z + 2*q[3]*x + q[4]*y*y +
               2*q[5]*y*z + 2*q[6]*y + q[7]*z*z + 2*q[8]*z + q[9];
    }
};

static int N,M;
static vector<Vec3> P0;          // original vertices
static vector<Tri> T0;           // original faces, zero-based
static Vec3 BMin,BMax;
static double diagLen=1.0, tau=0.05, tau2=0.0025;
static chrono::steady_clock::time_point TSTART;
static inline double elapsed(){ return chrono::duration<double>(chrono::steady_clock::now()-TSTART).count(); }

static inline uint64_t edgeKey(int a,int b){ if(a>b) swap(a,b); return (uint64_t)(uint32_t)a<<32 | (uint32_t)b; }
static inline int ea(uint64_t k){ return (int)(k>>32); }
static inline int eb(uint64_t k){ return (int)(k&0xffffffffu); }

struct RenderMaps{
    int R=0, PIX=0;
    vector<float> depth; // 6*PIX
    vector<float> norm;  // 6*PIX*3, normal color [0,255]
};

static inline void projectView(const Vec3&p,int view,int R,double&u,double&v,double&z){
    const double D=2.5; double f=800.0*(double)R/1024.0, c=0.5*(double)R;
    double sx,sy;
    if(view==0){ sx=p.y;  sy=p.z; z=D-p.x; }
    else if(view==1){ sx=-p.y; sy=p.z; z=D+p.x; }
    else if(view==2){ sx=-p.x; sy=p.z; z=D-p.y; }
    else if(view==3){ sx=p.x;  sy=p.z; z=D+p.y; }
    else if(view==4){ sx=p.x;  sy=p.y; z=D-p.z; }
    else { sx=-p.x; sy=p.y; z=D+p.z; }
    u=f*sx/z+c; v=f*sy/z+c;
}

static void renderMesh(const vector<Vec3>&V,const vector<Tri>&F,RenderMaps&rm,int R){
    const int PIX=R*R; rm.R=R; rm.PIX=PIX;
    rm.depth.assign((size_t)6*PIX,255.0f);
    rm.norm.assign((size_t)6*PIX*3,127.5f);
    const int nv=(int)V.size();
    vector<float> U(nv), VV(nv), Z(nv);
    vector<Vec3> FN(F.size());
    for(size_t i=0;i<F.size();++i){
        const Tri&t=F[i]; Vec3 cr=cross3(V[t.b]-V[t.a],V[t.c]-V[t.a]); FN[i]=normalize(cr);
    }
    for(int view=0;view<6;++view){
        for(int i=0;i<nv;++i){ double u,v,z; projectView(V[i],view,R,u,v,z); U[i]=(float)u; VV[i]=(float)v; Z[i]=(float)z; }
        float* zbuf=rm.depth.data()+(size_t)view*PIX;
        float* nbuf=rm.norm.data()+(size_t)view*PIX*3;
        for(size_t ti=0;ti<F.size();++ti){
            const Tri&t=F[ti]; int a=t.a,b=t.b,c=t.c;
            float x0=U[a],y0=VV[a],z0=Z[a], x1=U[b],y1=VV[b],z1=Z[b], x2=U[c],y2=VV[c],z2=Z[c];
            if(!(z0>0&&z1>0&&z2>0)) continue;
            float minx=min(x0,min(x1,x2)), maxx=max(x0,max(x1,x2));
            float miny=min(y0,min(y1,y2)), maxy=max(y0,max(y1,y2));
            int X0=max(0,(int)floor(minx-0.5f)), X1=min(R-1,(int)ceil(maxx+0.5f));
            int Y0=max(0,(int)floor(miny-0.5f)), Y1=min(R-1,(int)ceil(maxy+0.5f));
            if(X0>X1||Y0>Y1) continue;
            float den=(y1-y2)*(x0-x2)+(x2-x1)*(y0-y2);
            if(fabs(den)<1e-20f) continue;
            float inv=1.0f/den;
            Vec3 n=FN[ti]; float nr=(float)((n.x+1.0)*127.5), ng=(float)((n.y+1.0)*127.5), nb=(float)((n.z+1.0)*127.5);
            for(int yy=Y0; yy<=Y1; ++yy){
                float py=(float)yy+0.5f; int row=yy*R;
                for(int xx=X0; xx<=X1; ++xx){
                    float px=(float)xx+0.5f;
                    float w0=((y1-y2)*(px-x2)+(x2-x1)*(py-y2))*inv;
                    float w1=((y2-y0)*(px-x2)+(x0-x2)*(py-y2))*inv;
                    float w2=1.0f-w0-w1;
                    if(w0<-1e-6f||w1<-1e-6f||w2<-1e-6f) continue;
                    float dep=1.0f/(w0/z0+w1/z1+w2/z2);
                    int id=row+xx;
                    if(dep<zbuf[id]){
                        zbuf[id]=dep; float*q=nbuf+id*3; q[0]=nr; q[1]=ng; q[2]=nb;
                    }
                }
            }
        }
    }
}

static inline double rectSum(const vector<double>&I,int W,int x0,int y0,int x1,int y1){
    return I[(size_t)y1*W+x1]-I[(size_t)y0*W+x1]-I[(size_t)y1*W+x0]+I[(size_t)y0*W+x0];
}
static double ssimChannel(const float*X,int sx,const float*Y,int sy,const float*dX,const float*dY,int R,
                          vector<double>&IX,vector<double>&IY,vector<double>&IX2,vector<double>&IY2,vector<double>&IXY){
    const int W=R+1; size_t SZ=(size_t)W*W;
    fill(IX.begin(),IX.begin()+SZ,0); fill(IY.begin(),IY.begin()+SZ,0); fill(IX2.begin(),IX2.begin()+SZ,0); fill(IY2.begin(),IY2.begin()+SZ,0); fill(IXY.begin(),IXY.begin()+SZ,0);
    for(int y=1;y<=R;y++){
        double ax=0,ay=0,ax2=0,ay2=0,axy=0; int row=(y-1)*R;
        for(int x=1;x<=R;x++){
            int p=row+x-1; double vx=X[(size_t)p*sx], vy=Y[(size_t)p*sy];
            ax+=vx; ay+=vy; ax2+=vx*vx; ay2+=vy*vy; axy+=vx*vy;
            size_t id=(size_t)y*W+x, up=(size_t)(y-1)*W+x;
            IX[id]=IX[up]+ax; IY[id]=IY[up]+ay; IX2[id]=IX2[up]+ax2; IY2[id]=IY2[up]+ay2; IXY[id]=IXY[up]+axy;
        }
    }
    const int rad=5, area=121; const double c1=6.5025, c2=58.5225;
    long long cnt=0; long double acc=0;
    for(int y=rad;y<R-rad;y++){
        int row=y*R;
        for(int x=rad;x<R-rad;x++){
            int p=row+x; if(!(dX[p]<254.0f || dY[p]<254.0f)) continue;
            int x0=x-rad,y0=y-rad,x1=x+rad+1,y1=y+rad+1;
            double vx=rectSum(IX,W,x0,y0,x1,y1), vy=rectSum(IY,W,x0,y0,x1,y1);
            double vx2=rectSum(IX2,W,x0,y0,x1,y1), vy2=rectSum(IY2,W,x0,y0,x1,y1), vxy=rectSum(IXY,W,x0,y0,x1,y1);
            double mx=vx/area,my=vy/area;
            double sx2=max(0.0,vx2/area-mx*mx), sy2=max(0.0,vy2/area-my*my), cov=vxy/area-mx*my;
            double den=(mx*mx+my*my+c1)*(sx2+sy2+c2);
            double val=den!=0.0?((2*mx*my+c1)*(2*cov+c2)/den):1.0;
            acc+=val; ++cnt;
        }
    }
    return cnt? (double)(acc/cnt) : 1.0;
}
static double renderSSIM(const RenderMaps&A,const RenderMaps&B){
    int R=A.R, PIX=R*R, W=R+1; size_t SZ=(size_t)W*W;
    vector<double>IX(SZ),IY(SZ),IX2(SZ),IY2(SZ),IXY(SZ);
    double total=0;
    for(int view=0;view<6;view++){
        const float*ad=A.depth.data()+(size_t)view*PIX; const float*bd=B.depth.data()+(size_t)view*PIX;
        double ns=0;
        for(int ch=0;ch<3;ch++){
            const float*an=A.norm.data()+(size_t)view*PIX*3+ch; const float*bn=B.norm.data()+(size_t)view*PIX*3+ch;
            ns += ssimChannel(an,3,bn,3,ad,bd,R,IX,IY,IX2,IY2,IXY);
        }
        ns/=3.0;
        double ds=ssimChannel(ad,1,bd,1,ad,bd,R,IX,IY,IX2,IY2,IXY);
        total += 0.5*(ns+ds);
    }
    return total/6.0;
}

struct Simplifier{
    vector<Vec3> P;
    vector<Face> F;
    vector<vector<int>> inc;
    vector<Quadric> Q;
    vector<array<float,3>> bbMin,bbMax;
    vector<unsigned char> alive;
    vector<int> ver, markA, markB, remap;
    int stampA=1,stampB=1000000;
    int activeV=0, activeF=0;
    double visualTol=24.0, visualTol2=576.0, minDot=0.12;
    struct Cand{ double cost; int u,v,vu,vv; bool operator<(const Cand&o)const{return cost>o.cost;} };
    priority_queue<Cand> pq;
    struct Snapshot{ vector<Vec3> V; vector<Tri> F; double ratio=1.0; double score=-1.0; };
    vector<Snapshot> shots;

    Vec3 faceNormal(int fid) const{
        const Face&f=F[fid]; return normalize(cross3(P[f.v[1]]-P[f.v[0]],P[f.v[2]]-P[f.v[0]]));
    }
    bool hasVert(int fid,int v) const{ const Face&f=F[fid]; return f.v[0]==v||f.v[1]==v||f.v[2]==v; }
    bool solveOptimal(const Quadric&q,Vec3&out){
        double a00=q.q[0],a01=q.q[1],a02=q.q[2],a11=q.q[4],a12=q.q[5],a22=q.q[7];
        double b0=-q.q[3],b1=-q.q[6],b2=-q.q[8];
        double det=a00*(a11*a22-a12*a12)-a01*(a01*a22-a12*a02)+a02*(a01*a12-a11*a02);
        if(fabs(det)<1e-14) return false;
        double dx=b0*(a11*a22-a12*a12)-a01*(b1*a22-a12*b2)+a02*(b1*a12-a11*b2);
        double dy=a00*(b1*a22-a12*b2)-b0*(a01*a22-a12*a02)+a02*(a01*b2-b1*a02);
        double dz=a00*(a11*b2-b1*a12)-a01*(a01*b2-b1*a02)+b0*(a01*a12-a11*a02);
        out={dx/det,dy/det,dz/det};
        return isfinite(out.x)&&isfinite(out.y)&&isfinite(out.z)&&out.x>=-2&&out.x<=2&&out.y>=-2&&out.y<=2&&out.z>=-2&&out.z<=2;
    }
    inline void mergedBox(int u,int v,array<float,3>&mn,array<float,3>&mx){
        for(int k=0;k<3;k++){ mn[k]=min(bbMin[u][k],bbMin[v][k]); mx[k]=max(bbMax[u][k],bbMax[v][k]); }
    }
    bool coversBox(const Vec3&p,const array<float,3>&mn,const array<float,3>&mx){
        for(int mask=0;mask<8;mask++){
            Vec3 c{(mask&1)?mx[0]:mn[0],(mask&2)?mx[1]:mn[1],(mask&4)?mx[2]:mn[2]};
            if(norm2(p-c)>tau2) return false;
        }
        return true;
    }
    bool visualBoxOK(const Vec3&p,const array<float,3>&mn,const array<float,3>&mx){
        double pu[6],pv[6],pd;
        for(int v=0;v<6;v++) projectView(p,v,1024,pu[v],pv[v],pd);
        for(int mask=0;mask<8;mask++){
            Vec3 c{(mask&1)?mx[0]:mn[0],(mask&2)?mx[1]:mn[1],(mask&4)?mx[2]:mn[2]};
            for(int v=0;v<6;v++){
                double u,w,d; projectView(c,v,1024,u,w,d);
                double dx=u-pu[v], dy=w-pv[v]; if(dx*dx+dy*dy>visualTol2) return false;
            }
        }
        return true;
    }
    bool computeCandidate(int u,int v,Vec3&best,double&bestCost){
        if(u==v||!alive[u]||!alive[v]) return false;
        array<float,3> mn,mx; mergedBox(u,v,mn,mx);
        Quadric q=Q[u]; q.add(Q[v]);
        Vec3 cand[8]; int cc=0; Vec3 opt;
        if(solveOptimal(q,opt)) cand[cc++]=opt;
        cand[cc++]=P[u]; cand[cc++]=P[v]; cand[cc++]=(P[u]+P[v])*0.5;
        cand[cc++]={0.5*(mn[0]+mx[0]),0.5*(mn[1]+mx[1]),0.5*(mn[2]+mx[2])};
        cand[cc++]=P[u]*0.75+P[v]*0.25; cand[cc++]=P[u]*0.25+P[v]*0.75;
        bestCost=1e300; bool ok=false; double len2=norm2(P[u]-P[v]);
        for(int i=0;i<cc;i++){
            Vec3 p=cand[i]; if(!coversBox(p,mn,mx)) continue; if(!visualBoxOK(p,mn,mx)) continue;
            double c=q.eval(p); if(c<0&&c>-1e-10)c=0;
            c += 1e-7*len2 + 1e-12*norm2(p);
            if(c<bestCost){ bestCost=c; best=p; ok=true; }
        }
        return ok;
    }
    void cleanup(int u){
        if(!alive[u]) return;
        vector<int>&a=inc[u]; int w=0;
        for(int fid:a) if(fid>=0&&fid<(int)F.size()&&F[fid].active&&hasVert(fid,u)) a[w++]=fid;
        a.resize(w);
    }
    void maybeCleanup(int u){
        if(!alive[u]) return; if(inc[u].size()<96) return;
        int good=0; for(int fid:inc[u]) if(F[fid].active&&hasVert(fid,u)) ++good;
        if((int)inc[u].size()>good*3+64) cleanup(u);
    }
    bool linkOK(int u,int v){
        if(u==v||!alive[u]||!alive[v]) return false;
        maybeCleanup(u); maybeCleanup(v);
        if(++stampA>900000000){ fill(markA.begin(),markA.end(),0); stampA=1; }
        if(++stampB>1900000000){ fill(markB.begin(),markB.end(),0); stampB=1000000; }
        int edgeFaces=0;
        for(int fid:inc[u]){
            if(!F[fid].active||!hasVert(fid,u)) continue;
            bool hasv=false; const Face&f=F[fid];
            for(int k=0;k<3;k++){ int x=f.v[k]; if(x==v) hasv=true; if(x!=u&&x!=v) markA[x]=stampA; }
            if(hasv) ++edgeFaces;
        }
        if(edgeFaces!=2) return false;
        int common=0;
        for(int fid:inc[v]){
            if(!F[fid].active||!hasVert(fid,v)) continue;
            const Face&f=F[fid];
            for(int k=0;k<3;k++){
                int x=f.v[k]; if(x==u||x==v) continue;
                if(markA[x]==stampA && markB[x]!=stampB){ markB[x]=stampB; if(++common>2) return false; }
            }
        }
        return common==2;
    }
    bool normalOK(int u,int v,const Vec3&p){
        static vector<int> touched; touched.clear();
        for(int fid:inc[u]) if(F[fid].active&&hasVert(fid,u)) touched.push_back(fid);
        for(int fid:inc[v]) if(F[fid].active&&hasVert(fid,v)) touched.push_back(fid);
        sort(touched.begin(),touched.end()); touched.erase(unique(touched.begin(),touched.end()),touched.end());
        double eps=max(1e-30,1e-24*diagLen*diagLen);
        for(int fid:touched){
            Face&f=F[fid]; if(!f.active) continue;
            bool hu=hasVert(fid,u), hv=hasVert(fid,v); if(hu&&hv) continue;
            Vec3 a=(f.v[0]==u||f.v[0]==v)?p:P[f.v[0]];
            Vec3 b=(f.v[1]==u||f.v[1]==v)?p:P[f.v[1]];
            Vec3 c=(f.v[2]==u||f.v[2]==v)?p:P[f.v[2]];
            Vec3 cr=cross3(b-a,c-a); double l2=norm2(cr); if(!(l2>eps)||!isfinite(l2)) return false;
            Vec3 nn=cr/sqrt(l2); if(dot3(nn,f.n)<minDot) return false;
        }
        return true;
    }
    void pushEdge(int u,int v){
        if(u==v||!alive[u]||!alive[v]) return; if(u>v) swap(u,v);
        if(norm2(P[u]-P[v])>4.05*tau2) return;
        Vec3 p; double c; if(!computeCandidate(u,v,p,c)) return;
        pq.push({c,u,v,ver[u],ver[v]});
    }
    void collectNeighbors(int u,vector<int>&out){
        out.clear(); if(!alive[u]) return; maybeCleanup(u);
        if(++stampA>900000000){ fill(markA.begin(),markA.end(),0); stampA=1; }
        for(int fid:inc[u]) if(F[fid].active&&hasVert(fid,u)){
            const Face&f=F[fid];
            for(int k=0;k<3;k++){ int w=f.v[k]; if(w!=u&&alive[w]&&markA[w]!=stampA){ markA[w]=stampA; out.push_back(w); } }
        }
    }
    void collapse(int keep,int rem,const Vec3&p){
        cleanup(keep); cleanup(rem);
        vector<int> rfaces=inc[rem];
        for(int fid:rfaces){
            if(!F[fid].active||!hasVert(fid,rem)) continue;
            bool hasK=hasVert(fid,keep); Face&f=F[fid];
            if(hasK){ f.active=0; --activeF; }
            else{
                for(int k=0;k<3;k++) if(f.v[k]==rem) f.v[k]=keep;
                if(f.v[0]==f.v[1]||f.v[0]==f.v[2]||f.v[1]==f.v[2]){ f.active=0; --activeF; }
                else inc[keep].push_back(fid);
            }
        }
        P[keep]=p; Q[keep].add(Q[rem]);
        for(int k=0;k<3;k++){ bbMin[keep][k]=min(bbMin[keep][k],bbMin[rem][k]); bbMax[keep][k]=max(bbMax[keep][k],bbMax[rem][k]); }
        alive[rem]=0; --activeV; ++ver[keep]; ++ver[rem]; inc[rem].clear();
        cleanup(keep);
        for(int fid:inc[keep]) if(F[fid].active&&hasVert(fid,keep)) F[fid].n=faceNormal(fid);
        static vector<int> nei; collectNeighbors(keep,nei); for(int w:nei) pushEdge(keep,w);
    }
    Snapshot makeSnapshot(){
        Snapshot s; s.ratio=(double)activeV/(double)N; remap.assign(N,-1); s.V.reserve(activeV);
        for(int i=0;i<N;i++) if(alive[i]){ remap[i]=(int)s.V.size(); s.V.push_back(P[i]); }
        s.F.reserve(activeF);
        for(int i=0;i<M;i++) if(F[i].active){
            int a=F[i].v[0],b=F[i].v[1],c=F[i].v[2];
            if(a==b||a==c||b==c) continue; if(remap[a]<0||remap[b]<0||remap[c]<0) continue;
            s.F.push_back({remap[a],remap[b],remap[c]});
        }
        return s;
    }
    void init(){
        P=P0; F.assign(M,{}); inc.assign(N,{}); vector<int>deg(N,0);
        for(int i=0;i<M;i++){ F[i].v[0]=T0[i].a; F[i].v[1]=T0[i].b; F[i].v[2]=T0[i].c; F[i].active=1; ++deg[T0[i].a]; ++deg[T0[i].b]; ++deg[T0[i].c]; }
        for(int i=0;i<N;i++) inc[i].reserve(deg[i]+8);
        for(int i=0;i<M;i++){ inc[F[i].v[0]].push_back(i); inc[F[i].v[1]].push_back(i); inc[F[i].v[2]].push_back(i); }
        for(int i=0;i<M;i++) F[i].n=faceNormal(i);
        Q.assign(N,Quadric());
        for(int i=0;i<M;i++){
            Face&f=F[i]; Vec3 cr=cross3(P[f.v[1]]-P[f.v[0]],P[f.v[2]]-P[f.v[0]]); double l=norm3(cr); if(l<=1e-300) continue;
            Vec3 n=cr/l; double d=-dot3(n,P[f.v[0]]); double w=max(1e-10,0.5*l);
            for(int k=0;k<3;k++) Q[f.v[k]].addPlane(n.x,n.y,n.z,d,w);
        }
        bbMin.resize(N); bbMax.resize(N); alive.assign(N,1); ver.assign(N,0); markA.assign(N,0); markB.assign(N,0); remap.resize(N);
        for(int i=0;i<N;i++) bbMin[i]=bbMax[i]={(float)P[i].x,(float)P[i].y,(float)P[i].z};
        activeV=N; activeF=M;
        if(N<5000){ visualTol=34; minDot=0.02; }
        else if(N<60000){ visualTol=27; minDot=0.05; }
        else if(N<200000){ visualTol=22; minDot=0.08; }
        else { visualTol=18; minDot=0.10; }
        visualTol2=visualTol*visualTol;
        vector<uint64_t> edges; edges.reserve((size_t)M*3);
        for(const auto&t:T0){ edges.push_back(edgeKey(t.a,t.b)); edges.push_back(edgeKey(t.b,t.c)); edges.push_back(edgeKey(t.c,t.a)); }
        sort(edges.begin(),edges.end()); edges.erase(unique(edges.begin(),edges.end()),edges.end());
        vector<Cand> base; base.reserve(edges.size());
        for(uint64_t e:edges){ int u=ea(e),v=eb(e); Vec3 p; double c; if(computeCandidate(u,v,p,c)) base.push_back({c,u,v,0,0}); }
        pq=priority_queue<Cand>(less<Cand>(),move(base));
    }
    vector<double> checkpointRatios(){
        vector<double> r;
        if(N<1000) r={0.55,0.45,0.35,0.28,0.22,0.18,0.15,0.13,0.11};
        else if(N<8000) r={0.35,0.28,0.22,0.18,0.15,0.13,0.11,0.095,0.085};
        else if(N<60000) r={0.26,0.22,0.18,0.15,0.13,0.115,0.10,0.09,0.08,0.07};
        else if(N<200000) r={0.20,0.17,0.145,0.125,0.11,0.095,0.085,0.075,0.065};
        else r={0.16,0.135,0.115,0.10,0.088,0.078,0.068,0.058,0.05};
        return r;
    }
    void runCollapse(){
        vector<double> cps=checkpointRatios(); vector<int> targets;
        for(double r:cps){ int t=max(4,(int)ceil(N*r)); if(t<N && (targets.empty()||targets.back()!=t)) targets.push_back(t); }
        int cp=0; int finalTarget=targets.empty()?max(4,N/10):targets.back();
        const double workLimit = (N>200000?14.0:(N>60000?13.5:13.0));
        long long pops=0;
        while(activeV>finalTarget && !pq.empty()){
            if((++pops&2047)==0 && elapsed()>workLimit) break;
            Cand c=pq.top(); pq.pop(); int u=c.u,v=c.v; if(u>v) swap(u,v);
            if(u==v||!alive[u]||!alive[v]) continue;
            if(c.vu!=ver[u]||c.vv!=ver[v]){ pushEdge(u,v); continue; }
            if(!linkOK(u,v)) continue;
            Vec3 p; double cost; if(!computeCandidate(u,v,p,cost)) continue;
            int keep=u,rem=v; if(inc[v].size()>inc[u].size()) keep=v,rem=u;
            if(!normalOK(u,v,p)) continue;
            collapse(keep,rem,p);
            while(cp<(int)targets.size() && activeV<=targets[cp]){ shots.push_back(makeSnapshot()); ++cp; if(shots.size()>12) shots.erase(shots.begin()); }
        }
        if(shots.empty() || (int)shots.back().V.size()!=activeV) shots.push_back(makeSnapshot());
    }
};

static bool validBasic(const vector<Vec3>&V,const vector<Tri>&F){
    if(V.empty()||F.empty()||(int)V.size()>N) return false;
    double eps=max(1e-30,1e-24*diagLen*diagLen);
    vector<uint64_t> edges; edges.reserve(F.size()*3);
    for(auto&t:F){
        if(t.a<0||t.b<0||t.c<0||t.a>=(int)V.size()||t.b>=(int)V.size()||t.c>=(int)V.size()) return false;
        if(t.a==t.b||t.a==t.c||t.b==t.c) return false;
        if(norm2(cross3(V[t.b]-V[t.a],V[t.c]-V[t.a]))<=eps) return false;
        edges.push_back(edgeKey(t.a,t.b)); edges.push_back(edgeKey(t.b,t.c)); edges.push_back(edgeKey(t.c,t.a));
    }
    sort(edges.begin(),edges.end());
    for(size_t i=0;i<edges.size();){ size_t j=i+1; while(j<edges.size()&&edges[j]==edges[i]) ++j; if(j-i!=2) return false; i=j; }
    return true;
}



static int orientSignFromCenter(const vector<Vec3>&V,const vector<Tri>&F,const Vec3&c){
    long double sum=0; int cnt=0;
    int stride=max(1,(int)F.size()/50000);
    for(int i=0;i<(int)F.size();i+=stride){
        const Tri&t=F[i]; Vec3 a=V[t.a],b=V[t.b],cc=V[t.c];
        Vec3 cr=cross3(b-a,cc-a); Vec3 ctr=(a+b+cc)/3.0;
        double v=dot3(cr,ctr-c); if(fabs(v)>1e-20){ sum+=v; ++cnt; }
    }
    return sum>=0?1:-1;
}
static inline void addOrientedTri(vector<Vec3>&V,vector<Tri>&F,int a,int b,int c,const Vec3&center,int sign){
    Vec3 cr=cross3(V[b]-V[a],V[c]-V[a]); Vec3 ctr=(V[a]+V[b]+V[c])/3.0;
    if(dot3(cr,ctr-center)*sign<0) swap(b,c);
    F.push_back({a,b,c});
}
static bool buildEllipsoidMesh(const Vec3&cen,const Vec3 basis[3],const double rad[3],int lat,int lon,vector<Vec3>&V,vector<Tri>&F){
    if(lat<4||lon<8) return false;
    int nv=2+(lat-1)*lon; if(nv>N) return false;
    const double pi=acos(-1.0); V.clear(); F.clear(); V.reserve(nv); F.reserve(2*lat*lon);
    auto make=[&](double x,double y,double z){ return cen + basis[0]*(rad[0]*x) + basis[1]*(rad[1]*y) + basis[2]*(rad[2]*z); };
    V.push_back(make(0,0,1));
    for(int i=1;i<lat;i++){
        double th=pi*(double)i/(double)lat, st=sin(th), ct=cos(th);
        for(int j=0;j<lon;j++){ double ph=2*pi*(double)j/(double)lon; V.push_back(make(st*cos(ph),st*sin(ph),ct)); }
    }
    V.push_back(make(0,0,-1)); int bottom=(int)V.size()-1;
    auto id=[&](int ring,int j){ return 1+(ring-1)*lon+((j%lon+lon)%lon); };
    int sg=orientSignFromCenter(V, vector<Tri>{{0,id(1,1),id(1,0)},{bottom,id(lat-1,0),id(lat-1,1)}}, cen);
    for(int j=0;j<lon;j++) addOrientedTri(V,F,0,id(1,j+1),id(1,j),cen,sg);
    for(int r=1;r<lat-1;r++) for(int j=0;j<lon;j++){
        int a=id(r,j), b=id(r,j+1), c=id(r+1,j), d=id(r+1,j+1);
        addOrientedTri(V,F,a,b,c,cen,sg); addOrientedTri(V,F,b,d,c,cen,sg);
    }
    for(int j=0;j<lon;j++) addOrientedTri(V,F,bottom,id(lat-1,j),id(lat-1,j+1),cen,sg);
    return validBasic(V,F);
}
static void jacobiAxes(double A[3][3],Vec3 ax[3]){
    double V[3][3]={{1,0,0},{0,1,0},{0,0,1}};
    for(int it=0;it<36;it++){
        int p=0,q=1; double best=fabs(A[0][1]);
        if(fabs(A[0][2])>best) p=0,q=2,best=fabs(A[0][2]);
        if(fabs(A[1][2])>best) p=1,q=2,best=fabs(A[1][2]);
        if(best<1e-18) break;
        double app=A[p][p], aqq=A[q][q], apq=A[p][q];
        double tauj=(aqq-app)/(2*apq); double t=(tauj>=0?1.0:-1.0)/(fabs(tauj)+sqrt(1+tauj*tauj));
        double c=1.0/sqrt(1+t*t), s=t*c;
        for(int k=0;k<3;k++) if(k!=p&&k!=q){ double akp=A[k][p], akq=A[k][q]; A[k][p]=A[p][k]=c*akp-s*akq; A[k][q]=A[q][k]=s*akp+c*akq; }
        A[p][p]=c*c*app-2*s*c*apq+s*s*aqq; A[q][q]=s*s*app+2*s*c*apq+c*c*aqq; A[p][q]=A[q][p]=0;
        for(int k=0;k<3;k++){ double vkp=V[k][p], vkq=V[k][q]; V[k][p]=c*vkp-s*vkq; V[k][q]=s*vkp+c*vkq; }
    }
    int ord[3]={0,1,2}; sort(ord,ord+3,[&](int a,int b){return A[a][a]>A[b][b];});
    for(int i=0;i<3;i++){ int c=ord[i]; ax[i]=normalize(Vec3{V[0][c],V[1][c],V[2][c]}); }
    if(norm2(cross3(ax[0],ax[1])-ax[2])>norm2(cross3(ax[0],ax[1])+ax[2])) ax[2]=ax[2]*-1.0;
}
static bool fitEllipsoid(const Vec3 basis[3],Vec3&cen,double rad[3],double&rms,double&mxerr){
    double lo[3]={1e100,1e100,1e100}, hi[3]={-1e100,-1e100,-1e100};
    for(const Vec3&p:P0) for(int k=0;k<3;k++){ double t=dot3(p,basis[k]); lo[k]=min(lo[k],t); hi[k]=max(hi[k],t); }
    cen={0,0,0}; for(int k=0;k<3;k++){ rad[k]=0.5*(hi[k]-lo[k]); if(!(rad[k]>1e-12)) return false; cen=cen+basis[k]*(0.5*(hi[k]+lo[k])); }
    int stride=max(1,N/240000), cnt=0; long double ss=0; mxerr=0;
    for(int i=0;i<N;i+=stride){ Vec3 q=P0[i]-cen; double r2=0; for(int k=0;k<3;k++){ double u=dot3(q,basis[k])/rad[k]; r2+=u*u; } double e=fabs(sqrt(max(0.0,r2))-1.0); ss+=e*e; mxerr=max(mxerr,e); ++cnt; }
    if(cnt<100) return false; rms=sqrt((double)(ss/cnt));
    double minr=min(rad[0],min(rad[1],rad[2])), maxr=max(rad[0],max(rad[1],rad[2]));
    return maxr>1e-12 && minr/maxr>0.08;
}

static inline void projWUV(const Vec3&p,const Vec3&w,const Vec3&u,const Vec3&v,double&t,double&x,double&y);

struct CapsuleFit{bool ok=false;Vec3 w,u,v;double ct=0,cu=0,cv=0,r=0,h=0,err=1e100;};
static CapsuleFit fitCapsuleAxis(Vec3 w,Vec3 hint){
    CapsuleFit fit; fit.w=normalize(w); hint=hint-fit.w*dot3(hint,fit.w); if(norm2(hint)<1e-20) hint=fabs(fit.w.x)<0.8?Vec3{1,0,0}:Vec3{0,1,0}; fit.u=normalize(hint); fit.v=normalize(cross3(fit.w,fit.u));
    if(N<1000||norm2(fit.w)<.5) return fit;
    double t0=1e100,t1=-1e100,x0=1e100,x1=-1e100,y0=1e100,y1=-1e100;
    for(const Vec3&p:P0){double t,x,y;projWUV(p,fit.w,fit.u,fit.v,t,x,y);t0=min(t0,t);t1=max(t1,t);x0=min(x0,x);x1=max(x1,x);y0=min(y0,y);y1=max(y1,y);} double ex=x1-x0, ey=y1-y0, et=t1-t0; if(!(ex>1e-12&&ey>1e-12&&et>1e-12)) return fit; if(min(ex,ey)<max(ex,ey)*.965) return fit;
    fit.ct=.5*(t0+t1); fit.cu=.5*(x0+x1); fit.cv=.5*(y0+y1); fit.r=.25*(ex+ey); fit.h=.5*et-fit.r; if(!(fit.r>1e-12)||fit.h<fit.r*.10) return fit;
    int stride=max(1,N/240000), cnt=0, side=0, cap=0, bad=0; long double se=0; double me=0;
    for(int i=0;i<N;i+=stride){double t,x,y;projWUV(P0[i],fit.w,fit.u,fit.v,t,x,y);double dt=t-fit.ct, xx=x-fit.cu, yy=y-fit.cv, rr=sqrt(xx*xx+yy*yy), ad=fabs(dt); double rel; if(ad<=fit.h){rel=fabs(rr-fit.r)/fit.r; if(ad<fit.h*.7)++side;} else {double q=sqrt(rr*rr+(ad-fit.h)*(ad-fit.h)); rel=fabs(q-fit.r)/fit.r; if(ad>fit.h+fit.r*.2)++cap;} if(ad>fit.h+fit.r*1.03) rel=1e9; if(rel>.09){ if(++bad>max(3,cnt/40+2)) return fit; } else {se+=rel*rel; me=max(me,rel);} ++cnt;}
    if(cnt<120||side<max(15,cnt/40)||cap<max(15,cnt/60)) return fit; fit.err=sqrt((double)(se/cnt)); if(fit.err>.025||me>.09) return fit; fit.ok=true; return fit;
}
static bool buildCapsuleMesh(const CapsuleFit&f,int sides,int cn,int cj,vector<Vec3>&V,vector<Tri>&F){
    if(!f.ok||sides<12||cn<3||cj<1) return false; const double pi=acos(-1.0); int rings=2*cn+cj-1; int nv=2+rings*sides; if(nv>N) return false; V.clear();F.clear();V.reserve(nv);F.reserve(2*rings*sides); Vec3 cen=f.w*f.ct+f.u*f.cu+f.v*f.cv; auto make=[&](double rr,double tt,int j){double th=2*pi*j/sides; return f.w*(f.ct+tt)+f.u*(f.cu+rr*cos(th))+f.v*(f.cv+rr*sin(th));}; V.push_back(make(0,f.h+f.r,0)); vector<int>rs; auto ring=[&](double rr,double tt){int st=V.size();rs.push_back(st);for(int j=0;j<sides;j++)V.push_back(make(rr,tt,j));}; for(int i=1;i<=cn;i++){double ph=.5*pi*i/cn; ring(f.r*sin(ph), f.h+f.r*cos(ph));} for(int i=1;i<=cj;i++){double a=(double)i/cj; ring(f.r, f.h*(1-2*a));} for(int i=1;i<cn;i++){double ph=.5*pi+.5*pi*i/cn; ring(f.r*sin(ph), -f.h+f.r*cos(ph));} int bot=V.size(); V.push_back(make(0,-f.h-f.r,0)); auto rid=[&](int r,int j){return rs[r]+((j%sides+sides)%sides);}; int sg=1; for(int j=0;j<sides;j++) addOrientedTri(V,F,0,rid(0,j+1),rid(0,j),cen,sg); for(int r=0;r+1<rings;r++) for(int j=0;j<sides;j++){int a=rid(r,j),b=rid(r,j+1),c=rid(r+1,j),d=rid(r+1,j+1);addOrientedTri(V,F,a,b,c,cen,sg);addOrientedTri(V,F,b,d,c,cen,sg);} for(int j=0;j<sides;j++) addOrientedTri(V,F,bot,rid(rings-1,j),rid(rings-1,j+1),cen,sg); return validBasic(V,F);
}

struct FrustumFit{bool ok=false;Vec3 w,u,v;double cu=0,cv=0,t0=0,t1=0,r0=0,r1=0,err=1e100;};
static inline void projWUV(const Vec3&p,const Vec3&w,const Vec3&u,const Vec3&v,double&t,double&x,double&y){t=dot3(p,w);x=dot3(p,u);y=dot3(p,v);} 
static inline Vec3 fromWUV(const FrustumFit&f,double t,double x,double y){return f.w*t+f.u*x+f.v*y;}
static FrustumFit fitFrustumAxis(Vec3 w,Vec3 hint){
    FrustumFit fit; fit.w=normalize(w); hint=hint-fit.w*dot3(hint,fit.w); if(norm2(hint)<1e-20) hint=fabs(fit.w.x)<0.8?Vec3{1,0,0}:Vec3{0,1,0}; fit.u=normalize(hint); fit.v=normalize(cross3(fit.w,fit.u));
    if(norm2(fit.w)<.5||norm2(fit.u)<.5||norm2(fit.v)<.5||N<800) return fit;
    double t0=1e100,t1=-1e100,x0=1e100,x1=-1e100,y0=1e100,y1=-1e100;
    for(const Vec3&p:P0){double t,x,y;projWUV(p,fit.w,fit.u,fit.v,t,x,y);t0=min(t0,t);t1=max(t1,t);x0=min(x0,x);x1=max(x1,x);y0=min(y0,y);y1=max(y1,y);} 
    double len=t1-t0; if(!(len>1e-12)) return fit; fit.t0=t0; fit.t1=t1; fit.cu=.5*(x0+x1); fit.cv=.5*(y0+y1);
    int stride=max(1,N/240000), cnt=0, axial=0; double maxr=0;
    for(int i=0;i<N;i+=stride){double t,x,y;projWUV(P0[i],fit.w,fit.u,fit.v,t,x,y);x-=fit.cu;y-=fit.cv;maxr=max(maxr,sqrt(x*x+y*y));++cnt;} if(!(maxr>1e-12)) return fit;
    double eps=maxr*.04; long double S=0,St=0,Stt=0,Sr=0,Str=0; int used=0, centerPts=0;
    for(int i=0;i<N;i+=stride){double t,x,y;projWUV(P0[i],fit.w,fit.u,fit.v,t,x,y);x-=fit.cu;y-=fit.cv;double r=sqrt(x*x+y*y); if(r<=eps){ if(min(fabs(t-t0),fabs(t-t1))>len*.06) return fit; ++centerPts; continue;} double ss=(t-t0)/len; S+=1; St+=ss; Stt+=ss*ss; Sr+=r; Str+=ss*r; ++used; }
    if(used<120) return fit; double den=(double)(S*Stt-St*St); if(fabs(den)<1e-18) return fit;
    double b=(double)((S*Str-St*Sr)/den), a=(double)((Sr-b*St)/S); double r0=a, r1=a+b; if(r0<-.04*maxr||r1<-.04*maxr) return fit; if(fabs(r0)<.04*maxr) r0=0; if(fabs(r1)<.04*maxr) r1=0; if(max(r0,r1)<.2*maxr) return fit;
    long double se=0; double me=0; int checked=0; for(int i=0;i<N;i+=stride){double t,x,y;projWUV(P0[i],fit.w,fit.u,fit.v,t,x,y);x-=fit.cu;y-=fit.cv;double r=sqrt(x*x+y*y); if(r<=eps) continue; double ss=(t-t0)/len; double pred=max(0.0,r0+(r1-r0)*ss); double e=fabs(r-pred)/maxr; se+=e*e; me=max(me,e); ++checked;}
    if(checked<120) return fit; fit.err=sqrt((double)(se/checked)); if(fit.err>.012||me>.055) return fit; fit.r0=r0; fit.r1=r1; fit.ok=true; return fit;
}
static bool buildFrustumMesh(const FrustumFit&f,int sides,vector<Vec3>&V,vector<Tri>&F){
    if(!f.ok||sides<12) return false; const double pi=acos(-1.0); double eps=max(f.r0,f.r1)*1e-7; bool cone0=f.r0<=eps, cone1=f.r1<=eps; if(cone0&&cone1) return false; V.clear(); F.clear(); Vec3 cen=fromWUV(f,.5*(f.t0+f.t1),f.cu,f.cv);
    auto ringPt=[&](double t,double r,int j){double th=2*pi*j/sides; return fromWUV(f,t,f.cu+r*cos(th),f.cv+r*sin(th));};
    int sg=1;
    if(!cone0&&!cone1){int c0=0,c1=1;V.push_back(fromWUV(f,f.t0,f.cu,f.cv));V.push_back(fromWUV(f,f.t1,f.cu,f.cv));int r0s=V.size();for(int j=0;j<sides;j++)V.push_back(ringPt(f.t0,f.r0,j));int r1s=V.size();for(int j=0;j<sides;j++)V.push_back(ringPt(f.t1,f.r1,j)); if((int)V.size()>N)return false; for(int j=0;j<sides;j++){int k=(j+1)%sides;addOrientedTri(V,F,r0s+j,r0s+k,r1s+j,cen,sg);addOrientedTri(V,F,r0s+k,r1s+k,r1s+j,cen,sg);addOrientedTri(V,F,c0,r0s+j,r0s+k,cen,sg);addOrientedTri(V,F,c1,r1s+k,r1s+j,cen,sg);} }
    else{bool bottom=cone0; double at=bottom?f.t0:f.t1, bt=bottom?f.t1:f.t0, br=bottom?f.r1:f.r0; int apex=0,bc=1;V.push_back(fromWUV(f,at,f.cu,f.cv));V.push_back(fromWUV(f,bt,f.cu,f.cv));int rs=V.size();for(int j=0;j<sides;j++)V.push_back(ringPt(bt,br,j)); if((int)V.size()>N)return false; for(int j=0;j<sides;j++){int k=(j+1)%sides;addOrientedTri(V,F,apex,rs+j,rs+k,cen,sg);addOrientedTri(V,F,bc,rs+k,rs+j,cen,sg);} }
    return validBasic(V,F);
}

static bool primitiveCandidates(vector<Vec3>&bestV,vector<Tri>&bestF){
    if(N<600||elapsed()>3.0) return false;
    int R=N<30000?512:(N<150000?256:192); RenderMaps orig,cand; renderMesh(P0,T0,orig,R);
    bool have=false; int bestN=N; double bestSc=-1;
    auto consider=[&](const vector<Vec3>&V,const vector<Tri>&F,double th){
        if(V.empty()||F.empty()||(int)V.size()>=bestN||!validBasic(V,F)) return;
        renderMesh(V,F,cand,R); double sc=renderSSIM(orig,cand);
        if(sc>=th && ((int)V.size()<bestN || sc>bestSc+.015)){ bestV=V; bestF=F; bestN=(int)V.size(); bestSc=sc; have=true; }
    };
    Vec3 mean{0,0,0}; for(const Vec3&p:P0) mean=mean+p; mean=mean/(double)max(1,N);
    double C[3][3]={{0,0,0},{0,0,0},{0,0,0}}; int stride=max(1,N/240000), cnt=0;
    for(int i=0;i<N;i+=stride){ Vec3 q=P0[i]-mean; double x[3]={q.x,q.y,q.z}; for(int a=0;a<3;a++) for(int b=0;b<3;b++) C[a][b]+=x[a]*x[b]; ++cnt; }
    for(int a=0;a<3;a++) for(int b=0;b<3;b++) C[a][b]/=max(1,cnt);
    Vec3 axes[3]; jacobiAxes(C,axes);
    Vec3 identities[3]={{1,0,0},{0,1,0},{0,0,1}};
    vector<array<Vec3,3>> bases; bases.push_back({identities[0],identities[1],identities[2]}); if(norm2(axes[0])>.5) bases.push_back({axes[0],axes[1],axes[2]});
    for(auto &ba:bases){
        Vec3 basis[3]={ba[0],ba[1],ba[2]}, cen; double rad[3],rms,mxerr;
        if(!fitEllipsoid(basis,cen,rad,rms,mxerr)) continue;
        double thFit = N<5000?0.018:0.014, mxFit=N<5000?0.075:0.055;
        if(rms>thFit||mxerr>mxFit) continue;
        vector<pair<int,int>> trials;
        if(N<1500) trials={{12,24},{14,28},{16,32}};
        else if(N<8000) trials={{14,28},{16,32},{18,36}};
        else if(N<40000) trials={{16,32},{18,36},{20,40},{22,44}};
        else trials={{20,40},{24,48},{28,56}};
        for(auto [lat,lon]:trials){ if(elapsed()>6.0) break; vector<Vec3>V; vector<Tri>F; if(buildEllipsoidMesh(cen,basis,rad,lat,lon,V,F)) consider(V,F,(R>=512?0.895:0.918)); }
    }
    vector<pair<Vec3,Vec3>> axlist;
    axlist.push_back({Vec3{1,0,0},Vec3{0,1,0}}); axlist.push_back({Vec3{0,1,0},Vec3{1,0,0}}); axlist.push_back({Vec3{0,0,1},Vec3{1,0,0}});
    if(norm2(axes[0])>.5){ axlist.push_back({axes[0],axes[1]}); axlist.push_back({axes[1],axes[2]}); axlist.push_back({axes[2],axes[0]}); }
    if(elapsed()<6.3){
        for(auto pr:axlist){ if(elapsed()>7.5) break; FrustumFit ff=fitFrustumAxis(pr.first,pr.second); if(!ff.ok) continue; int sidesList[5]={24,32,40,48,64}; for(int sides:sidesList){ vector<Vec3>V; vector<Tri>F; if(buildFrustumMesh(ff,sides,V,F)) consider(V,F,(R>=512?0.895:0.918)); } }
    }
    if(elapsed()<7.8){
        for(auto pr:axlist){ if(elapsed()>8.8) break; CapsuleFit cf=fitCapsuleAxis(pr.first,pr.second); if(!cf.ok) continue; int sidesList[4]={24,32,40,48}; int cnList[2]={6,8}; for(int sides:sidesList) for(int cn:cnList){ vector<Vec3>V; vector<Tri>F; if(buildCapsuleMesh(cf,sides,cn,1,V,F)) consider(V,F,(R>=512?0.895:0.918)); } }
    }
    return have;
}

static bool tryRegularTorusGrid(vector<Vec3>&OV, vector<Tri>&OF){
    // Detect a common UV torus / periodic rectangular grid: N = U*V, M = 2*N.
    if(N<1000 || M!=2*N || elapsed()>1.5) return false;
    int bestU=0,bestV=0;
    vector<pair<int,int>> cand;
    for(int d=6;(long long)d*d<=N;d++) if(N%d==0){ cand.push_back({N/d,d}); cand.push_back({d,N/d}); }
    for(auto [U,V]:cand){
        if(U<8||V<6||U>5000||V>5000) continue;
        int step=max(1,N/600), ok=0, tot=0;
        for(int cell=0; cell<N && tot<600; cell+=step,++tot){
            int i=cell/V,j=cell%V; int a=i*V+j,b=i*V+(j+1)%V,c=((i+1)%U)*V+j,d=((i+1)%U)*V+(j+1)%V;
            int f=2*cell; if(f+1>=M) break;
            auto same=[&](const Tri&t,int x,int y,int z){ array<int,3>A{t.a,t.b,t.c},B{x,y,z}; sort(A.begin(),A.end()); sort(B.begin(),B.end()); return A==B; };
            bool good=(same(T0[f],a,b,c)&&same(T0[f+1],b,d,c))||(same(T0[f],a,b,d)&&same(T0[f+1],a,d,c))||(same(T0[f],a,c,d)&&same(T0[f+1],a,d,b));
            if(good) ++ok;
        }
        if(tot>100 && ok*100>=tot*98){ bestU=U; bestV=V; break; }
    }
    if(!bestU) return false;
    int U=bestU,V=bestV;
    int U2=max(8,min(U, (int)round(sqrt((double)N)*0.42)));
    int V2=max(8,min(V, (int)round(sqrt((double)N)*0.42)));
    if(U2*V2>=N) return false;
    OV.clear(); OF.clear(); OV.reserve(U2*V2); OF.reserve(2*U2*V2);
    for(int i=0;i<U2;i++){ int oi=(long long)i*U/U2; for(int j=0;j<V2;j++){ int oj=(long long)j*V/V2; OV.push_back(P0[oi*V+oj]); } }
    auto id=[&](int i,int j){ return ((i%U2+U2)%U2)*V2 + ((j%V2+V2)%V2); };
    for(int i=0;i<U2;i++) for(int j=0;j<V2;j++){ OF.push_back({id(i,j),id(i+1,j),id(i+1,j+1)}); OF.push_back({id(i,j),id(i+1,j+1),id(i,j+1)}); }
    if(!validBasic(OV,OF)) return false;
    int R=N<40000?512:256; RenderMaps A,B; renderMesh(P0,T0,A,R); renderMesh(OV,OF,B,R); double s=renderSSIM(A,B);
    return s>=(N<40000?0.925:0.935);
}

static void chooseAndWrite(const vector<Vec3>&V,const vector<Tri>&F){
    static char obuf[1<<20]; setvbuf(stdout,obuf,_IOFBF,sizeof(obuf));
    printf("%d %d\n",(int)V.size(),(int)F.size());
    const char* fmt = ((int)V.size()*2<=N) ? "v %.15g %.15g %.15g\n" : "v %.10g %.10g %.10g\n";
    for(const Vec3&p:V) printf(fmt,p.x,p.y,p.z);
    for(const Tri&t:F) printf("f %d %d %d\n",t.a+1,t.b+1,t.c+1);
}

int main(){
    TSTART=chrono::steady_clock::now();
    FastInput in; N=in.nextInt(); M=in.nextInt(); P0.resize(N); T0.resize(M);
    BMin={1e100,1e100,1e100}; BMax={-1e100,-1e100,-1e100};
    for(int i=0;i<N;i++){
        (void)in.ch(); P0[i].x=in.nextDouble(); P0[i].y=in.nextDouble(); P0[i].z=in.nextDouble();
        BMin.x=min(BMin.x,P0[i].x); BMin.y=min(BMin.y,P0[i].y); BMin.z=min(BMin.z,P0[i].z);
        BMax.x=max(BMax.x,P0[i].x); BMax.y=max(BMax.y,P0[i].y); BMax.z=max(BMax.z,P0[i].z);
    }
    for(int i=0;i<M;i++){ (void)in.ch(); T0[i].a=in.nextInt()-1; T0[i].b=in.nextInt()-1; T0[i].c=in.nextInt()-1; }
    diagLen=norm3(BMax-BMin); if(!(diagLen>0)) diagLen=1; tau=0.05*diagLen*0.995; tau2=tau*tau;

    // Official sample and tiny meshes: preserve the known redundant-vertex win using the general path if possible.
    vector<Vec3> specialV; vector<Tri> specialF;
    if(tryRegularTorusGrid(specialV,specialF)){ chooseAndWrite(specialV,specialF); return 0; }
    if(primitiveCandidates(specialV,specialF)){ chooseAndWrite(specialV,specialF); return 0; }

    Simplifier S; S.init(); S.runCollapse();
    vector<Vec3> bestV; vector<Tri> bestF;
    if(S.shots.empty()){
        bestV=P0; bestF=T0;
    }else{
        int R = (N<8000?512:(N<80000?384:(N<250000?256:192)));
        double guard = (R>=512?0.918:(R>=384?0.925:0.935));
        RenderMaps orig; renderMesh(P0,T0,orig,R);
        RenderMaps cand; int chosen=-1; double bestScore=-1; int bestIdx=0;
        // Evaluate from most compressed to safest; keep first passing high guard.
        for(int i=(int)S.shots.size()-1;i>=0;i--){
            if(elapsed()>19.2 && chosen<0) break;
            if(!validBasic(S.shots[i].V,S.shots[i].F)) continue;
            renderMesh(S.shots[i].V,S.shots[i].F,cand,R);
            double sc=renderSSIM(orig,cand); S.shots[i].score=sc;
            if(sc>bestScore){ bestScore=sc; bestIdx=i; }
            if(sc>=guard){ chosen=i; break; }
        }
        if(chosen<0) chosen=bestIdx;
        bestV=move(S.shots[chosen].V); bestF=move(S.shots[chosen].F);
        if(!validBasic(bestV,bestF) || bestV.empty() || bestF.empty() || (int)bestV.size()>N){ bestV=P0; bestF=T0; }
    }
    chooseAndWrite(bestV,bestF);
    return 0;
}
