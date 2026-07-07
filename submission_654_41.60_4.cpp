#include <bits/stdc++.h>
using namespace std;

struct Vec3{
    double x=0,y=0,z=0;
    Vec3(){}
    Vec3(double X,double Y,double Z):x(X),y(Y),z(Z){}
    Vec3 operator+(const Vec3&o)const{return {x+o.x,y+o.y,z+o.z};}
    Vec3 operator-(const Vec3&o)const{return {x-o.x,y-o.y,z-o.z};}
    Vec3 operator*(double s)const{return {x*s,y*s,z*s};}
    Vec3 operator/(double s)const{return {x/s,y/s,z/s};}
};
static inline double dotv(const Vec3&a,const Vec3&b){return a.x*b.x+a.y*b.y+a.z*b.z;}
static inline Vec3 crossv(const Vec3&a,const Vec3&b){return {a.y*b.z-a.z*b.y,a.z*b.x-a.x*b.z,a.x*b.y-a.y*b.x};}
static inline double norm2(const Vec3&a){return dotv(a,a);} 
static inline double normv(const Vec3&a){return sqrt(norm2(a));}
static inline double dist2(const Vec3&a,const Vec3&b){return norm2(a-b);} 
static inline Vec3 normalized(Vec3 a){double l=normv(a); return l>1e-300?a/l:Vec3{0,0,0};}

struct Face{int v[3]; unsigned char active=1; Vec3 n;};
struct Tri{int a,b,c;};

struct FastInput{
    vector<char> buf; char* p=nullptr;
    FastInput(){ char tmp[1<<16]; size_t n; buf.reserve(1<<27); while((n=fread(tmp,1,sizeof(tmp),stdin))>0) buf.insert(buf.end(),tmp,tmp+n); buf.push_back(0); p=buf.data(); }
    inline void skip(){ while(*p==' '||*p=='\n'||*p=='\r'||*p=='\t') ++p; }
    int nextInt(){ skip(); int s=1; if(*p=='-') s=-1,++p; int x=0; while(*p>='0'&&*p<='9'){x=x*10+(*p-'0');++p;} return x*s; }
    double nextDouble(){ skip(); char* e=nullptr; double x=strtod(p,&e); p=e; return x; }
    char nextChar(){ skip(); return *p++; }
};

struct Quadric{
    double q[10];
    Quadric(){ memset(q,0,sizeof(q)); }
    void addPlane(double a,double b,double c,double d,double w=1.0){
        q[0]+=w*a*a; q[1]+=w*a*b; q[2]+=w*a*c; q[3]+=w*a*d;
        q[4]+=w*b*b; q[5]+=w*b*c; q[6]+=w*b*d;
        q[7]+=w*c*c; q[8]+=w*c*d; q[9]+=w*d*d;
    }
    void add(const Quadric&o){ for(int i=0;i<10;i++) q[i]+=o.q[i]; }
    double eval(const Vec3&p)const{
        double x=p.x,y=p.y,z=p.z;
        return q[0]*x*x + 2*q[1]*x*y + 2*q[2]*x*z + 2*q[3]*x +
               q[4]*y*y + 2*q[5]*y*z + 2*q[6]*y + q[7]*z*z + 2*q[8]*z + q[9];
    }
};

static int N=0,M=0;
static vector<Vec3> P, Orig;
static vector<Face> F;
static vector<Tri> OrigTri;
static vector<vector<int>> inc;
static vector<Quadric> Qv;
static vector<array<float,3>> bbMin, bbMax;
static vector<unsigned char> alive;
static vector<int> ver, markA, markB, tempId;
static int stampA=1, stampB=1000007;
static int activeV=0, activeF=0;
static double diagLen=1.0, tau=0.05, tau2=0.0025;
static chrono::steady_clock::time_point T0;
static double timeLimit=19.35;
static double visualTol=22.0, visualTol2=484.0;

static inline double elapsed(){return chrono::duration<double>(chrono::steady_clock::now()-T0).count();}
static inline bool timeOK(double margin=0){return elapsed()+margin<timeLimit;}
static inline uint64_t ekey(int a,int b){ if(a>b) swap(a,b); return (uint64_t)(uint32_t)a<<32 | (uint32_t)b; }
static inline bool faceContains(const Face&f,int u){return f.v[0]==u||f.v[1]==u||f.v[2]==u;}
static Vec3 faceNormal(const Face&f){
    Vec3 cr=crossv(P[f.v[1]]-P[f.v[0]],P[f.v[2]]-P[f.v[0]]);
    double l=normv(cr); return l>1e-300?cr/l:Vec3{0,0,0};
}

static void readInput(){
    FastInput in; N=in.nextInt(); M=in.nextInt();
    P.resize(N); Orig.resize(N); Vec3 mn{1e100,1e100,1e100}, mx{-1e100,-1e100,-1e100};
    for(int i=0;i<N;i++){
        (void)in.nextChar(); P[i].x=in.nextDouble(); P[i].y=in.nextDouble(); P[i].z=in.nextDouble(); Orig[i]=P[i];
        mn.x=min(mn.x,P[i].x); mn.y=min(mn.y,P[i].y); mn.z=min(mn.z,P[i].z);
        mx.x=max(mx.x,P[i].x); mx.y=max(mx.y,P[i].y); mx.z=max(mx.z,P[i].z);
    }
    F.resize(M); OrigTri.reserve(M); vector<int> deg(N,0);
    for(int i=0;i<M;i++){
        (void)in.nextChar(); int a=in.nextInt()-1,b=in.nextInt()-1,c=in.nextInt()-1;
        F[i].v[0]=a; F[i].v[1]=b; F[i].v[2]=c; F[i].active=1;
        OrigTri.push_back({a,b,c}); ++deg[a]; ++deg[b]; ++deg[c];
    }
    Vec3 d=mx-mn; diagLen=normv(d); if(!(diagLen>0)) diagLen=1;
    tau=0.05*diagLen*0.9995; tau2=tau*tau;
    inc.assign(N,{}); for(int i=0;i<N;i++) inc[i].reserve(deg[i]+8);
    for(int i=0;i<M;i++){inc[F[i].v[0]].push_back(i); inc[F[i].v[1]].push_back(i); inc[F[i].v[2]].push_back(i);}    
    activeV=N; activeF=M;
}

static bool validateMesh(const vector<Vec3>&V,const vector<Tri>&T){
    if(V.empty()||T.empty()||V.size()>(size_t)N) return false;
    const double areaEps=max(1e-30,1e-24*diagLen*diagLen);
    vector<uint64_t> edges; edges.reserve(T.size()*3);
    vector<array<int,3>> faces; faces.reserve(T.size());
    for(const auto&t:T){
        if(t.a<0||t.b<0||t.c<0||t.a>=(int)V.size()||t.b>=(int)V.size()||t.c>=(int)V.size()) return false;
        if(t.a==t.b||t.a==t.c||t.b==t.c) return false;
        Vec3 cr=crossv(V[t.b]-V[t.a],V[t.c]-V[t.a]); if(!(norm2(cr)>areaEps)) return false;
        edges.push_back(ekey(t.a,t.b)); edges.push_back(ekey(t.b,t.c)); edges.push_back(ekey(t.c,t.a));
        array<int,3> f{t.a,t.b,t.c}; sort(f.begin(),f.end()); faces.push_back(f);
    }
    sort(faces.begin(),faces.end()); if(adjacent_find(faces.begin(),faces.end())!=faces.end()) return false;
    sort(edges.begin(),edges.end());
    for(size_t i=0;i<edges.size();){ size_t j=i+1; while(j<edges.size()&&edges[j]==edges[i]) ++j; if(j-i!=2) return false; i=j; }
    return true;
}

static void outputMesh(const vector<Vec3>&V,const vector<Tri>&T){
    static char obuf[1<<20]; setvbuf(stdout,obuf,_IOFBF,sizeof(obuf));
    printf("%d %d\n",(int)V.size(),(int)T.size());
    bool highPrec = (int)V.size()*2 <= N;
    for(const auto&p:V){
        if(highPrec) printf("v %.15g %.15g %.15g\n",p.x,p.y,p.z);
        else printf("v %.10g %.10g %.10g\n",p.x,p.y,p.z);
    }
    for(const auto&t:T) printf("f %d %d %d\n",t.a+1,t.b+1,t.c+1);
}
static void outputOriginal(){ outputMesh(Orig, OrigTri); }

struct RenderMaps{int R=0; vector<float> depth; vector<float> norm;};
static inline void projectPoint(const Vec3&p,int view,int R,double&u,double&v,double&dep){
    double f=800.0*(double)R/1024.0, c=0.5*(double)R; double sx,sy;
    if(view==0){sx=p.y; sy=p.z; dep=2.5-p.x;}
    else if(view==1){sx=-p.y; sy=p.z; dep=2.5+p.x;}
    else if(view==2){sx=-p.x; sy=p.z; dep=2.5-p.y;}
    else if(view==3){sx=p.x; sy=p.z; dep=2.5+p.y;}
    else if(view==4){sx=p.x; sy=p.y; dep=2.5-p.z;}
    else {sx=-p.x; sy=p.y; dep=2.5+p.z;}
    u=f*sx/dep+c; v=f*sy/dep+c;
}

static void renderMesh(const vector<Vec3>&V,const vector<Tri>&T,RenderMaps&rm,int R){
    const int PIX=R*R, VIEWS=6; rm.R=R; rm.depth.assign((size_t)VIEWS*PIX,255.f); rm.norm.assign((size_t)VIEWS*PIX*3,127.5f);
    vector<float> U(V.size()), W(V.size()), Z(V.size());
    vector<Vec3> FN(T.size());
    for(size_t i=0;i<T.size();++i){ const auto&t=T[i]; FN[i]=normalized(crossv(V[t.b]-V[t.a],V[t.c]-V[t.a])); }
    double focal=800.0*(double)R/1024.0, center=0.5*(double)R;
    for(int view=0;view<VIEWS;++view){
        for(size_t i=0;i<V.size();++i){ double u,v,z; projectPoint(V[i],view,R,u,v,z); U[i]=(float)u; W[i]=(float)v; Z[i]=(float)z; }
        float*zbuf=rm.depth.data()+(size_t)view*PIX; float*nbuf=rm.norm.data()+(size_t)view*PIX*3;
        for(size_t ti=0;ti<T.size();++ti){
            const auto&t=T[ti]; int ia=t.a, ib=t.b, ic=t.c;
            float x0=U[ia],x1=U[ib],x2=U[ic], y0=W[ia],y1=W[ib],y2=W[ic], z0=Z[ia],z1=Z[ib],z2=Z[ic];
            if(!(z0>0&&z1>0&&z2>0)) continue;
            float minx=min(x0,min(x1,x2)), maxx=max(x0,max(x1,x2));
            float miny=min(y0,min(y1,y2)), maxy=max(y0,max(y1,y2));
            int xi0=max(0,(int)floor(minx-0.5f)), xi1=min(R-1,(int)ceil(maxx+0.5f));
            int yi0=max(0,(int)floor(miny-0.5f)), yi1=min(R-1,(int)ceil(maxy+0.5f));
            if(xi0>xi1||yi0>yi1) continue;
            float den=(y1-y2)*(x0-x2)+(x2-x1)*(y0-y2); if(fabs(den)<1e-20f) continue; float inv=1.0f/den;
            Vec3 n=FN[ti]; float nr=(float)((n.x+1)*127.5), ng=(float)((n.y+1)*127.5), nb=(float)((n.z+1)*127.5);
            for(int y=yi0;y<=yi1;++y){ float py=(float)y+0.5f; int row=y*R; for(int x=xi0;x<=xi1;++x){ float px=(float)x+0.5f;
                float w0=((y1-y2)*(px-x2)+(x2-x1)*(py-y2))*inv;
                float w1=((y2-y0)*(px-x2)+(x0-x2)*(py-y2))*inv;
                float w2=1.f-w0-w1;
                if(w0>=-1e-6f&&w1>=-1e-6f&&w2>=-1e-6f){ float dep=1.f/(w0/z0+w1/z1+w2/z2); int id=row+x; if(dep<zbuf[id]){ zbuf[id]=dep; float*q=nbuf+id*3; q[0]=nr; q[1]=ng; q[2]=nb; } }
            }}
        }
    }
}
static inline double rectSum(const vector<double>&I,int W,int x0,int y0,int x1,int y1){return I[(size_t)y1*W+x1]-I[(size_t)y0*W+x1]-I[(size_t)y1*W+x0]+I[(size_t)y0*W+x0];}
static double ssimChannel(const float*X,int sx,const float*Y,int sy,const float*dX,const float*dY,int R,vector<double>&IX,vector<double>&IY,vector<double>&IX2,vector<double>&IY2,vector<double>&IXY){
    int W=R+1; size_t SZ=(size_t)W*W; fill(IX.begin(),IX.begin()+SZ,0); fill(IY.begin(),IY.begin()+SZ,0); fill(IX2.begin(),IX2.begin()+SZ,0); fill(IY2.begin(),IY2.begin()+SZ,0); fill(IXY.begin(),IXY.begin()+SZ,0);
    for(int y=1;y<=R;y++){ double ax=0,ay=0,ax2=0,ay2=0,axy=0; int row=(y-1)*R; for(int x=1;x<=R;x++){ int p=row+x-1; double xv=X[(size_t)p*sx], yv=Y[(size_t)p*sy]; ax+=xv; ay+=yv; ax2+=xv*xv; ay2+=yv*yv; axy+=xv*yv; size_t id=(size_t)y*W+x, up=(size_t)(y-1)*W+x; IX[id]=IX[up]+ax; IY[id]=IY[up]+ay; IX2[id]=IX2[up]+ax2; IY2[id]=IY2[up]+ay2; IXY[id]=IXY[up]+axy; }}
    const int rad=5, area=121; const double c1=6.5025,c2=58.5225; long long cnt=0; long double acc=0;
    for(int y=rad;y<R-rad;y++){ int row=y*R; for(int x=rad;x<R-rad;x++){ int p=row+x; if(!(dX[p]<254.f||dY[p]<254.f)) continue; int x0=x-rad,y0=y-rad,x1=x+rad+1,y1=y+rad+1; double sumx=rectSum(IX,W,x0,y0,x1,y1), sumy=rectSum(IY,W,x0,y0,x1,y1); double sumx2=rectSum(IX2,W,x0,y0,x1,y1), sumy2=rectSum(IY2,W,x0,y0,x1,y1), sumxy=rectSum(IXY,W,x0,y0,x1,y1); double mx=sumx/area,my=sumy/area; double vx=sumx2/area-mx*mx, vy=sumy2/area-my*my, cov=sumxy/area-mx*my; if(vx<0&&vx>-1e-7)vx=0; if(vy<0&&vy>-1e-7)vy=0; double den=(mx*mx+my*my+c1)*(vx+vy+c2); double val=den?((2*mx*my+c1)*(2*cov+c2)/den):1.0; acc+=val; ++cnt; }}
    return cnt?(double)(acc/cnt):1.0;
}
static double renderSSIM(const RenderMaps&A,const RenderMaps&B){
    int R=A.R, PIX=R*R, W=R+1; size_t SZ=(size_t)W*W; vector<double>IX(SZ),IY(SZ),IX2(SZ),IY2(SZ),IXY(SZ); double tot=0; const int VIEWS=6;
    for(int v=0;v<VIEWS;v++){ const float*ad=A.depth.data()+(size_t)v*PIX; const float*bd=B.depth.data()+(size_t)v*PIX; double ns=0; for(int ch=0;ch<3;ch++){ const float*an=A.norm.data()+(size_t)v*PIX*3+ch; const float*bn=B.norm.data()+(size_t)v*PIX*3+ch; ns+=ssimChannel(an,3,bn,3,ad,bd,R,IX,IY,IX2,IY2,IXY);} ns/=3; double ds=ssimChannel(ad,1,bd,1,ad,bd,R,IX,IY,IX2,IY2,IXY); tot += 0.5*(ns+ds); }
    return tot/6.0;
}

static RenderMaps *cachedOrig=nullptr; static int cachedR=0;
static double proxySSIM(const vector<Vec3>&V,const vector<Tri>&T,int R){
    static RenderMaps orig; if(cachedOrig==nullptr||cachedR!=R){ vector<Vec3> ov=Orig; renderMesh(ov,OrigTri,orig,R); cachedOrig=&orig; cachedR=R; }
    RenderMaps cand; renderMesh(V,T,cand,R); return renderSSIM(*cachedOrig,cand);
}

static void compactIncident(int u){ if(!alive[u]) return; vector<int>&lst=inc[u]; if(lst.size()<96) return; int good=0; for(int id:lst) if(id>=0&&id<M&&F[id].active&&faceContains(F[id],u)) ++good; if((int)lst.size()<=good*3+64) return; vector<int> n; n.reserve(good+8); for(int id:lst) if(F[id].active&&faceContains(F[id],u)) n.push_back(id); lst.swap(n); }
static bool edgeOpposites(int u,int v,int opp[2]){
    compactIncident(u); compactIncident(v); int cnt=0;
    const vector<int>&a = inc[u].size()<inc[v].size()?inc[u]:inc[v];
    for(int fid:a){ if(!F[fid].active) continue; const Face&f=F[fid]; if(!faceContains(f,u)||!faceContains(f,v)) continue; if(cnt>=2) return false; int t=-1; for(int k=0;k<3;k++){int x=f.v[k]; if(x!=u&&x!=v) t=x;} if(t<0) return false; opp[cnt++]=t; }
    return cnt==2 && opp[0]!=opp[1];
}
static bool linkOK(int u,int v){
    int opp[2]; if(!edgeOpposites(u,v,opp)) return false;
    if(++stampA>1000000000){fill(markA.begin(),markA.end(),0); stampA=1;}
    if(++stampB>1000000000){fill(markB.begin(),markB.end(),0); stampB=1000007;}
    for(int fid:inc[u]) if(F[fid].active&&faceContains(F[fid],u)){ const Face&f=F[fid]; for(int k=0;k<3;k++){int x=f.v[k]; if(x!=u&&x!=v) markA[x]=stampA; }}
    int common=0, got0=0, got1=0;
    for(int fid:inc[v]) if(F[fid].active&&faceContains(F[fid],v)){ const Face&f=F[fid]; for(int k=0;k<3;k++){ int x=f.v[k]; if(x==u||x==v) continue; if(markA[x]==stampA && markB[x]!=stampB){ markB[x]=stampB; ++common; if(x==opp[0]) got0=1; if(x==opp[1]) got1=1; if(common>2) return false; if(x!=opp[0]&&x!=opp[1]) return false; } }}
    return common==2&&got0&&got1;
}
static inline void mergedBB(int u,int v,array<float,3>&mn,array<float,3>&mx){ for(int k=0;k<3;k++){mn[k]=min(bbMin[u][k],bbMin[v][k]); mx[k]=max(bbMax[u][k],bbMax[v][k]);} }
static bool coversBBox(const Vec3&p,const array<float,3>&mn,const array<float,3>&mx){
    for(int mask=0;mask<8;mask++){ double x=(mask&1)?mx[0]:mn[0], y=(mask&2)?mx[1]:mn[1], z=(mask&4)?mx[2]:mn[2]; double dx=p.x-x,dy=p.y-y,dz=p.z-z; if(dx*dx+dy*dy+dz*dz>tau2) return false; }
    return true;
}
static bool visualBBoxOK(const Vec3&p,const array<float,3>&mn,const array<float,3>&mx){
    if(N<1000) return true;
    double pu[6],pv[6],pd; for(int v=0;v<6;v++) projectPoint(p,v,1024,pu[v],pv[v],pd);
    for(int mask=0;mask<8;mask++){ Vec3 c{(double)((mask&1)?mx[0]:mn[0]),(double)((mask&2)?mx[1]:mn[1]),(double)((mask&4)?mx[2]:mn[2])}; for(int v=0;v<6;v++){ double u,w,d; projectPoint(c,v,1024,u,w,d); double du=u-pu[v], dv=w-pv[v]; if(du*du+dv*dv>visualTol2) return false; }}
    return true;
}
static bool solveOpt(const Quadric&q,Vec3&out){
    double a00=q.q[0],a01=q.q[1],a02=q.q[2],a11=q.q[4],a12=q.q[5],a22=q.q[7]; double b0=-q.q[3],b1=-q.q[6],b2=-q.q[8];
    double det=a00*(a11*a22-a12*a12)-a01*(a01*a22-a12*a02)+a02*(a01*a12-a11*a02); if(fabs(det)<1e-14) return false;
    double dx=b0*(a11*a22-a12*a12)-a01*(b1*a22-a12*b2)+a02*(b1*a12-a11*b2);
    double dy=a00*(b1*a22-a12*b2)-b0*(a01*a22-a12*a02)+a02*(a01*b2-b1*a02);
    double dz=a00*(a11*b2-b1*a12)-a01*(a01*b2-b1*a02)+b0*(a01*a12-a11*a02);
    out={dx/det,dy/det,dz/det}; return isfinite(out.x)&&isfinite(out.y)&&isfinite(out.z)&&out.x>=-2&&out.x<=2&&out.y>=-2&&out.y<=2&&out.z>=-2&&out.z<=2;
}
struct Eval{bool ok=false; double score=1e300; Vec3 pos;};
static Eval computeCandidate(int u,int v){
    Eval best; if(u==v||!alive[u]||!alive[v]) return best; array<float,3>mn,mx; mergedBB(u,v,mn,mx);
    double bdx=mx[0]-mn[0],bdy=mx[1]-mn[1],bdz=mx[2]-mn[2]; if(bdx*bdx+bdy*bdy+bdz*bdz>4.0001*tau2) return best;
    Quadric q=Qv[u]; q.add(Qv[v]); Vec3 cands[8]; int c=0; Vec3 opt; if(solveOpt(q,opt)) cands[c++]=opt; cands[c++]=P[u]; cands[c++]=P[v]; cands[c++]=(P[u]+P[v])*0.5; cands[c++]={(mn[0]+mx[0])*0.5,(mn[1]+mx[1])*0.5,(mn[2]+mx[2])*0.5}; cands[c++]=P[u]*0.75+P[v]*0.25; cands[c++]=P[u]*0.25+P[v]*0.75;
    double len2=dist2(P[u],P[v]);
    for(int i=0;i<c;i++){ Vec3 p=cands[i]; if(!coversBBox(p,mn,mx)) continue; if(!visualBBoxOK(p,mn,mx)) continue; double e=q.eval(p); if(e<0&&e>-1e-10) e=0; double sc=e+2e-10*len2+1e-14*norm2(p); if(sc<best.score){best.ok=true; best.score=sc; best.pos=p;} }
    return best;
}
static bool flipOK(int keep,int rem,const Vec3&pos,double minDot){
    static vector<int> touched; touched.clear(); for(int fid:inc[keep]) if(F[fid].active&&faceContains(F[fid],keep)) touched.push_back(fid); for(int fid:inc[rem]) if(F[fid].active&&faceContains(F[fid],rem)) touched.push_back(fid); sort(touched.begin(),touched.end()); touched.erase(unique(touched.begin(),touched.end()),touched.end());
    double areaEps=max(1e-30,1e-24*diagLen*diagLen);
    for(int fid:touched){ Face&f=F[fid]; if(!f.active) continue; bool hk=faceContains(f,keep), hr=faceContains(f,rem); if(hk&&hr) continue; Vec3 a=(f.v[0]==keep||f.v[0]==rem)?pos:P[f.v[0]], b=(f.v[1]==keep||f.v[1]==rem)?pos:P[f.v[1]], c=(f.v[2]==keep||f.v[2]==rem)?pos:P[f.v[2]]; Vec3 cr=crossv(b-a,c-a); double l2=norm2(cr); if(!(l2>areaEps)||!isfinite(l2)) return false; Vec3 nn=cr/sqrt(l2); if(dotv(nn,f.n)<minDot) return false; }
    return true;
}
struct Node{double score; int u,v,vu,vv; bool operator<(const Node&o)const{return score>o.score;}};
static void pushEdge(priority_queue<Node>&pq,int u,int v){ if(u==v||!alive[u]||!alive[v]) return; if(u>v) swap(u,v); Eval e=computeCandidate(u,v); if(e.ok) pq.push({e.score,u,v,ver[u],ver[v]}); }
static void collectNeighbors(int u,vector<int>&out){ out.clear(); if(!alive[u]) return; compactIncident(u); if(++stampA>1000000000){fill(markA.begin(),markA.end(),0); stampA=1;} for(int fid:inc[u]) if(F[fid].active&&faceContains(F[fid],u)){ const Face&f=F[fid]; for(int k=0;k<3;k++){int w=f.v[k]; if(w!=u&&alive[w]&&markA[w]!=stampA){markA[w]=stampA; out.push_back(w);}}} }
static void doCollapse(int keep,int rem,const Vec3&pos,priority_queue<Node>&pq){
    compactIncident(keep); compactIncident(rem); vector<int> remFaces=inc[rem];
    for(int fid:remFaces){ if(!F[fid].active||!faceContains(F[fid],rem)) continue; Face&f=F[fid]; bool hasKeep=faceContains(f,keep); if(hasKeep){f.active=0; --activeF;} else { for(int k=0;k<3;k++) if(f.v[k]==rem) f.v[k]=keep; if(f.v[0]==f.v[1]||f.v[0]==f.v[2]||f.v[1]==f.v[2]){f.active=0;--activeF;} else {inc[keep].push_back(fid);} } }
    P[keep]=pos; for(int fid:inc[keep]) if(F[fid].active&&faceContains(F[fid],keep)) F[fid].n=faceNormal(F[fid]); Qv[keep].add(Qv[rem]); for(int k=0;k<3;k++){bbMin[keep][k]=min(bbMin[keep][k],bbMin[rem][k]); bbMax[keep][k]=max(bbMax[keep][k],bbMax[rem][k]);}
    alive[rem]=0; ++ver[keep]; ++ver[rem]; --activeV; inc[rem].clear(); compactIncident(keep); static vector<int> neigh; collectNeighbors(keep,neigh); for(int w:neigh) pushEdge(pq,keep,w);
}
struct Snapshot{vector<Vec3> V; vector<Tri> T; int nv=0; double ratio=1.0, ssim=-1;};
static Snapshot makeSnapshot(){ Snapshot s; s.nv=activeV; s.ratio=(double)activeV/max(1,N); tempId.assign(N,-1); s.V.reserve(activeV); for(int i=0;i<N;i++) if(alive[i]){tempId[i]=(int)s.V.size(); s.V.push_back(P[i]);} s.T.reserve(activeF); for(int i=0;i<M;i++) if(F[i].active){int a=F[i].v[0],b=F[i].v[1],c=F[i].v[2]; if(a!=b&&a!=c&&b!=c&&tempId[a]>=0&&tempId[b]>=0&&tempId[c]>=0) s.T.push_back({tempId[a],tempId[b],tempId[c]});} s.nv=(int)s.V.size(); return s; }
static void initDecimator(){
    activeV=N; activeF=M; alive.assign(N,1); ver.assign(N,0); markA.assign(N,0); markB.assign(N,0); tempId.assign(N,-1); Qv.assign(N,Quadric()); bbMin.resize(N); bbMax.resize(N);
    for(int i=0;i<N;i++){bbMin[i]={(float)P[i].x,(float)P[i].y,(float)P[i].z}; bbMax[i]=bbMin[i];}
    for(int i=0;i<M;i++){ F[i].active=1; F[i].n=faceNormal(F[i]); Vec3 cr=crossv(P[F[i].v[1]]-P[F[i].v[0]],P[F[i].v[2]]-P[F[i].v[0]]); double ar=normv(cr); if(ar>1e-300){ Vec3 n=cr/ar; double d=-dotv(n,P[F[i].v[0]]); double w=max(1e-8,0.5*ar)+2e-7; for(int k=0;k<3;k++) Qv[F[i].v[k]].addPlane(n.x,n.y,n.z,d,w); }}
}
static vector<Snapshot> runDecimator(double minRatio,bool makeSnaps){
    initDecimator(); priority_queue<Node> pq; vector<uint64_t> edges; edges.reserve((size_t)M*3);
    for(int i=0;i<M;i++){edges.push_back(ekey(F[i].v[0],F[i].v[1])); edges.push_back(ekey(F[i].v[1],F[i].v[2])); edges.push_back(ekey(F[i].v[2],F[i].v[0]));}
    sort(edges.begin(),edges.end()); edges.erase(unique(edges.begin(),edges.end()),edges.end()); for(uint64_t k:edges) pushEdge(pq,(int)(k>>32),(int)(k&0xffffffffu)); vector<uint64_t>().swap(edges);
    vector<double> ratios={.50,.35,.25,.20,.17,.15,.13,.115,.105,.095,.085,.075,.065,.055,.048}; vector<int> cps; int last=-1; if(makeSnaps){ for(double r:ratios){int c=max(4,(int)ceil(N*r)); if(c<N&&c!=last){cps.push_back(c); last=c;}} sort(cps.begin(),cps.end(),greater<int>()); }
    int target=max(4,(int)ceil(N*minRatio)); int cp=0; vector<Snapshot> snaps; long long pops=0; double minDot = (N<=70000? -0.04 : cos(38.0*acos(-1)/180.0));
    while(activeV>target&&!pq.empty()&&timeOK(0.35)){
        if((++pops&4095)==0 && !timeOK(0.35)) break;
        Node n=pq.top(); pq.pop(); int u=n.u,v=n.v; if(u==v||!alive[u]||!alive[v]) continue; if(ver[u]!=n.vu||ver[v]!=n.vv){ if(alive[u]&&alive[v]) pushEdge(pq,u,v); continue; }
        Eval e=computeCandidate(u,v); if(!e.ok) continue; if(!linkOK(u,v)) continue; int keep=u,rem=v; if(inc[v].size()>inc[u].size()) keep=v,rem=u; if(!flipOK(keep,rem,e.pos,minDot)) continue; doCollapse(keep,rem,e.pos,pq);
        while(makeSnaps && cp<(int)cps.size() && activeV<=cps[cp]){ snaps.push_back(makeSnapshot()); ++cp; if(!timeOK(0.35)) break; }
    }
    if(makeSnaps && (snaps.empty() || snaps.back().nv!=activeV)) snaps.push_back(makeSnapshot());
    if(!makeSnaps) snaps.push_back(makeSnapshot()); return snaps;
}

static bool tryToroidalGrid(vector<Vec3>&OV,vector<Tri>&OT){
    if(N<300||M!=2*N||elapsed()>2.0) return false;
    auto same=[&](const Tri&t,int a,int b,int c){array<int,3>x{t.a,t.b,t.c},y{a,b,c}; sort(x.begin(),x.end()); sort(y.begin(),y.end()); return x==y;};
    int U=0,V=0; vector<pair<int,int>> cand; for(int d=6;(long long)d*d<=N;d++) if(N%d==0){cand.push_back({N/d,d}); cand.push_back({d,N/d});}
    for(auto [u,v]:cand){ if(u<6||v<6||u>20000||v>20000) continue; int ok=0,tot=0,st=max(1,N/800); for(int cell=0;cell<N&&tot<800;cell+=st){int i=cell/v,j=cell-i*v,fid=2*cell; if(fid+1>=M) break; int a=i*v+j,b=((i+1)%u)*v+j,c=((i+1)%u)*v+(j+1)%v,d=i*v+(j+1)%v; if(same(OrigTri[fid],a,b,c)&&same(OrigTri[fid+1],a,c,d)) ++ok; ++tot;} if(tot>100&&ok*100>=tot*98){U=u;V=v;break;}}
    if(!U) return false;
    auto build=[&](int U2,int V2,vector<Vec3>&X,vector<Tri>&T)->bool{ if(U2<4||V2<4||U2*V2>=N) return false; X.clear(); T.clear(); X.reserve(U2*V2); T.reserve(2*U2*V2); for(int i=0;i<U2;i++){int oi=(long long)i*U/U2; for(int j=0;j<V2;j++){int oj=(long long)j*V/V2; X.push_back(Orig[oi*V+oj]);}} auto id=[&](int i,int j){return ((i%U2+U2)%U2)*V2+((j%V2+V2)%V2);}; for(int i=0;i<U2;i++)for(int j=0;j<V2;j++){T.push_back({id(i,j),id(i+1,j),id(i+1,j+1)}); T.push_back({id(i,j),id(i+1,j+1),id(i,j+1)});} double lim2=(0.0496*diagLen)*(0.0496*diagLen); for(int q=0;q<N;q++){int i=q/V,j=q-i*V, ii=(int)(((long long)i*U2+U/2)/U), jj=(int)(((long long)j*V2+V/2)/V); double best=1e100; for(int di=-1;di<=1;di++)for(int dj=-1;dj<=1;dj++)best=min(best,dist2(Orig[q],X[id(ii+di,jj+dj)])); if(best>lim2) return false;} return validateMesh(X,T);};
    vector<pair<int,int>> trials; for(double r: {10.0,8.0,6.0,5.0,4.0,3.0}) trials.push_back({max(6,(int)round(U/r)),max(6,(int)round(V/r))});
    int R = N<35000?512:384; double guard=N<35000?0.925:0.935; vector<Vec3>X; vector<Tri>T; bool got=false; int bestN=N; double bestScore=-1;
    for(auto [u2,v2]:trials){ if(elapsed()>8.5) break; if(!build(u2,v2,X,T)) continue; if((int)X.size()>=bestN) continue; double s=proxySSIM(X,T,R); if(s>=guard){OV=X;OT=T;bestN=X.size();got=true;bestScore=s;} }
    return got;
}

static bool tryLatLong(vector<Vec3>&OV,vector<Tri>&OT){
    if(N<300||elapsed()>2.0) return false; int lon=0, rings=0; bool ok=false;
    for(int v=8; v<=4096; ++v){ if((N-2)%v) continue; int r=(N-2)/v; if(r<3) continue; if(M!=2*v*r) continue; int good=0,tot=0,st=max(1,v/64); for(int j=0;j<v&&tot<128;j+=st){ int a=1+j,b=1+(j+1)%v; if(OrigTri[j].a==0 || OrigTri[j].b==0 || OrigTri[j].c==0) good++; tot++; } if(tot>0&&good*100>=tot*95){lon=v;rings=r;ok=true;break;} }
    if(!ok) return false;
    auto origId=[&](int r,int j){return 1+(r-1)*lon+((j%lon+lon)%lon);};
    auto build=[&](int R2,int V2,vector<Vec3>&X,vector<Tri>&T)->bool{ if(R2<3||V2<8||2+R2*V2>=N) return false; X.clear();T.clear(); X.push_back(Orig[0]); for(int r=1;r<=R2;r++){int orr=1+(long long)(r-1)*(rings-1)/max(1,R2-1); for(int j=0;j<V2;j++){int oj=(long long)j*lon/V2; X.push_back(Orig[origId(orr,oj)]);}} int bot=X.size(); X.push_back(Orig[N-1]); auto id=[&](int r,int j){return 1+(r-1)*V2+((j%V2+V2)%V2);}; for(int j=0;j<V2;j++) T.push_back({0,id(1,j+1),id(1,j)}); for(int r=1;r<R2;r++)for(int j=0;j<V2;j++){T.push_back({id(r,j),id(r,j+1),id(r+1,j)}); T.push_back({id(r,j+1),id(r+1,j+1),id(r+1,j)});} for(int j=0;j<V2;j++) T.push_back({bot,id(R2,j),id(R2,j+1)}); return validateMesh(X,T);};
    vector<pair<int,int>> trials; for(double d:{8.0,6.0,5.0,4.0,3.0}) trials.push_back({max(4,(int)round(rings/d)),max(12,(int)round(lon/d))});
    int R=N<35000?512:384; double guard=.93; vector<Vec3>X;vector<Tri>T; bool got=false; int best=N; for(auto [r2,v2]:trials){ if(elapsed()>8.5) break; if(!build(r2,v2,X,T)) continue; if((int)X.size()>=best) continue; double s=proxySSIM(X,T,R); if(s>=guard){OV=X;OT=T;best=X.size();got=true;}} return got;
}

static bool chooseSnapshot(vector<Snapshot>&snaps,vector<Vec3>&OV,vector<Tri>&OT){
    if(snaps.empty()) return false; int R=N<12000?512:(N<70000?384:256); double guard=N<12000?0.918:(N<70000?0.925:0.94); int chosen=-1; double best=-1; for(int i=(int)snaps.size()-1;i>=0;i--){ if(elapsed()>18.2) break; if(!validateMesh(snaps[i].V,snaps[i].T)) continue; double s=proxySSIM(snaps[i].V,snaps[i].T,R); snaps[i].ssim=s; if(s>best) best=s, chosen=i; if(s>=guard){OV=move(snaps[i].V); OT=move(snaps[i].T); return true;} }
    if(chosen>=0 && best>=0.89){OV=move(snaps[chosen].V); OT=move(snaps[chosen].T); return true;} return false;
}

int main(){
    T0=chrono::steady_clock::now(); readInput();
    vector<Vec3> OV; vector<Tri> OT;
    if(N<=12){ auto snaps=runDecimator(0.01,false); if(!snaps.empty()&&validateMesh(snaps[0].V,snaps[0].T)){outputMesh(snaps[0].V,snaps[0].T); return 0;} outputOriginal(); return 0; }
    // First try index-structured remeshing. These routes are guarded by vertex-Hausdorff, manifold validation and a renderer proxy.
    if(N<=80000 && timeOK(8.0)){
        if(tryToroidalGrid(OV,OT) || tryLatLong(OV,OT)){ if(validateMesh(OV,OT)){ outputMesh(OV,OT); return 0; } }
    }
    double minRatio;
    if(N<1500) minRatio=.12; else if(N<8000) minRatio=.065; else if(N<70000) minRatio=.050; else if(N<250000) minRatio=.080; else minRatio=.090;
    if(N<8000) visualTol=30; else if(N<70000) visualTol=24; else if(N<250000) visualTol=18; else visualTol=15; visualTol2=visualTol*visualTol;
    bool useSnaps = (N<=70000 && timeOK(5.0));
    auto snaps=runDecimator(minRatio,useSnaps);
    bool ok=false;
    if(useSnaps && timeOK(1.0)) ok=chooseSnapshot(snaps,OV,OT);
    if(!ok && !snaps.empty()){ Snapshot &s=snaps.back(); if(validateMesh(s.V,s.T)){OV=move(s.V); OT=move(s.T); ok=true;} }
    if(ok && OV.size()<Orig.size() && validateMesh(OV,OT)){ outputMesh(OV,OT); return 0; }
    outputOriginal();
    return 0;
}
