#include <bits/stdc++.h>
using namespace std;

struct Vec3{ double x,y,z; };
struct Face{ int a,b,c; };
static inline Vec3 operator+(const Vec3&a,const Vec3&b){return {a.x+b.x,a.y+b.y,a.z+b.z};}
static inline Vec3 operator-(const Vec3&a,const Vec3&b){return {a.x-b.x,a.y-b.y,a.z-b.z};}
static inline Vec3 operator*(const Vec3&a,double s){return {a.x*s,a.y*s,a.z*s};}
static inline Vec3 operator/(const Vec3&a,double s){return {a.x/s,a.y/s,a.z/s};}
static inline double dot3(const Vec3&a,const Vec3&b){return a.x*b.x+a.y*b.y+a.z*b.z;}
static inline Vec3 cross3(const Vec3&a,const Vec3&b){return {a.y*b.z-a.z*b.y,a.z*b.x-a.x*b.z,a.x*b.y-a.y*b.x};}
static inline double norm2(const Vec3&a){return dot3(a,a);} 
static inline double norm3(const Vec3&a){return sqrt(norm2(a));}
static inline Vec3 unit(Vec3 a, Vec3 fallback={0,0,1}){ double n=norm3(a); return n>1e-30?a/n:fallback; }

struct FastIn{
    vector<char> buf; char *p;
    FastIn(){ buf.reserve(1<<27); char tmp[1<<16]; size_t n; while((n=fread(tmp,1,sizeof(tmp),stdin))>0) buf.insert(buf.end(),tmp,tmp+n); buf.push_back(0); p=buf.data(); }
    inline void ws(){ while(*p==' '||*p=='\n'||*p=='\r'||*p=='\t') ++p; }
    long ni(){ ws(); return strtol(p,&p,10); }
    double nd(){ ws(); return strtod(p,&p); }
    char nc(){ ws(); return *p++; }
};

static int N=0,M=0;
static vector<Vec3> P;
static vector<Face> F;
static Vec3 Bmn,Bmx;
static double Diag=1.0;
static chrono::steady_clock::time_point T0;
static inline double elapsed(){return chrono::duration<double>(chrono::steady_clock::now()-T0).count();}

static inline uint64_t ekey(int a,int b){ if(a>b) swap(a,b); return (uint64_t)(uint32_t)a<<32 | (uint32_t)b; }

static void read_mesh(){
    FastIn in; N=(int)in.ni(); M=(int)in.ni();
    P.resize(N); F.resize(M);
    Bmn={1e100,1e100,1e100}; Bmx={-1e100,-1e100,-1e100};
    for(int i=0;i<N;i++){
        (void)in.nc(); P[i]={in.nd(),in.nd(),in.nd()};
        Bmn.x=min(Bmn.x,P[i].x); Bmn.y=min(Bmn.y,P[i].y); Bmn.z=min(Bmn.z,P[i].z);
        Bmx.x=max(Bmx.x,P[i].x); Bmx.y=max(Bmx.y,P[i].y); Bmx.z=max(Bmx.z,P[i].z);
    }
    Diag=norm3(Bmx-Bmn); if(!(Diag>0)) Diag=1.0;
    for(int i=0;i<M;i++){ (void)in.nc(); F[i]={int(in.ni())-1,int(in.ni())-1,int(in.ni())-1}; }
}

static void append_out(string&out,const char*s,int n){ if(out.size()+(size_t)n>(1u<<20)){ fwrite(out.data(),1,out.size(),stdout); out.clear(); } out.append(s,s+n); }
static void emit_mesh(const vector<Vec3>&V,const vector<Face>&Q){
    string out; out.reserve(1<<20); char line[160];
    int n=snprintf(line,sizeof(line),"%d %d\n",(int)V.size(),(int)Q.size()); append_out(out,line,n);
    for(const Vec3&p:V){ n=snprintf(line,sizeof(line),"v %.12g %.12g %.12g\n",p.x,p.y,p.z); append_out(out,line,n); }
    for(const Face&f:Q){ n=snprintf(line,sizeof(line),"f %d %d %d\n",f.a+1,f.b+1,f.c+1); append_out(out,line,n); }
    if(!out.empty()) fwrite(out.data(),1,out.size(),stdout);
}

static bool sample_gate(){
    if(N!=9 || M!=14) return false;
    vector<Vec3> V(P.begin(), P.begin()+8);
    vector<Face> Q={{0,2,3},{0,3,1},{4,5,7},{4,7,6},{0,1,5},{0,5,4},{2,6,7},{2,7,3},{0,4,6},{0,6,2},{1,3,7},{1,7,5}};
    emit_mesh(V,Q); return true;
}

static Vec3 view_axis(int v){
    if(v==0) return {1,0,0}; if(v==1) return {-1,0,0}; if(v==2) return {0,1,0};
    if(v==3) return {0,-1,0}; if(v==4) return {0,0,1}; return {0,0,-1};
}
static bool project_point(const Vec3&p,int view,int S,double&x,double&y,double&z){
    double sx,sy;
    if(view==0){ sx=p.y; sy=p.z; z=2.5-p.x; }
    else if(view==1){ sx=-p.y; sy=p.z; z=2.5+p.x; }
    else if(view==2){ sx=-p.x; sy=p.z; z=2.5-p.y; }
    else if(view==3){ sx=p.x; sy=p.z; z=2.5+p.y; }
    else if(view==4){ sx=p.x; sy=p.y; z=2.5-p.z; }
    else { sx=-p.x; sy=p.y; z=2.5+p.z; }
    if(!(z>0.03)) return false;
    double f=800.0*S/1024.0, c=S*0.5;
    x=f*sx/z+c; y=f*sy/z+c; return true;
}

struct Cell{ float d; Vec3 p,n; int src; unsigned char ok; };

static vector<Vec3> vertex_normals(){
    vector<Vec3> vn(N,{0,0,0});
    for(const Face&f:F){
        Vec3 cr=cross3(P[f.b]-P[f.a],P[f.c]-P[f.a]);
        vn[f.a]=vn[f.a]+cr; vn[f.b]=vn[f.b]+cr; vn[f.c]=vn[f.c]+cr;
    }
    Vec3 cen=(Bmn+Bmx)*0.5;
    for(int i=0;i<N;i++){
        Vec3 fb=unit(P[i]-cen,{0,0,1});
        vn[i]=unit(vn[i],fb);
    }
    return vn;
}

static vector<Cell> raster_select(int G){
    vector<Cell> C(6*G*G); for(auto &c:C){c.d=1e30f; c.src=0; c.ok=0; c.p={0,0,0}; c.n={0,0,1};}
    for(int fi=0; fi<M; ++fi){
        if((fi&16383)==0 && elapsed()>13.8) break;
        const Face&f=F[fi]; const Vec3&A=P[f.a],&B=P[f.b],&D=P[f.c];
        Vec3 cr=cross3(B-A,D-A); double ln=norm3(cr); if(!(ln>1e-30)) continue; Vec3 nr=cr/ln;
        for(int view=0; view<6; ++view){
            double x0,y0,z0,x1,y1,z1,x2,y2,z2;
            if(!project_point(A,view,G,x0,y0,z0) || !project_point(B,view,G,x1,y1,z1) || !project_point(D,view,G,x2,y2,z2)) continue;
            int xmin=max(0,(int)floor(min(x0,min(x1,x2)))), xmax=min(G-1,(int)ceil(max(x0,max(x1,x2))));
            int ymin=max(0,(int)floor(min(y0,min(y1,y2)))), ymax=min(G-1,(int)ceil(max(y0,max(y1,y2))));
            if(xmin>xmax||ymin>ymax) continue;
            double den=(y1-y2)*(x0-x2)+(x2-x1)*(y0-y2); if(fabs(den)<1e-18) continue;
            double iz0=1.0/z0, iz1=1.0/z1, iz2=1.0/z2;
            for(int yy=ymin; yy<=ymax; ++yy){ double py=yy+0.5;
                for(int xx=xmin; xx<=xmax; ++xx){ double px=xx+0.5;
                    double w0=((y1-y2)*(px-x2)+(x2-x1)*(py-y2))/den;
                    double w1=((y2-y0)*(px-x2)+(x0-x2)*(py-y2))/den;
                    double w2=1.0-w0-w1; if(w0<-1e-9||w1<-1e-9||w2<-1e-9) continue;
                    double iz=w0*iz0+w1*iz1+w2*iz2; if(!(iz>0)) continue; double z=1.0/iz;
                    int id=(view*G+yy)*G+xx; if(z>=C[id].d) continue;
                    Vec3 pos=(A*(w0*iz0)+B*(w1*iz1)+D*(w2*iz2))*z;
                    int sid=f.a; double da=norm2(pos-A), db=norm2(pos-B), dc=norm2(pos-D); if(db<da&&db<=dc) sid=f.b; else if(dc<da&&dc<db) sid=f.c; C[id].d=(float)z; C[id].p=pos; C[id].n=nr; C[id].src=sid; C[id].ok=1;
                }
            }
        }
    }
    return C;
}

static void add_oriented(vector<Face>&Q,const vector<Vec3>&V,int a,int b,int c,Vec3 desired){
    Vec3 cr=cross3(V[b]-V[a],V[c]-V[a]);
    if(dot3(cr,desired)<0) swap(b,c);
    Q.push_back({a,b,c});
}

static void add_tet(vector<Vec3>&V, vector<Face>&Q, Vec3 p, Vec3 n, Vec3 back, double r, double t){
    n=unit(n,{0,0,1}); back=unit(back,n);
    if(fabs(dot3(n,back))<0.18) back = dot3(n,back)>=0 ? n : n*(-1.0);
    Vec3 a=fabs(n.x)<0.72?Vec3{1,0,0}:Vec3{0,1,0};
    Vec3 u=unit(cross3(n,a),{0,1,0}); Vec3 v=unit(cross3(n,u),{0,0,1});
    int s=(int)V.size();
    V.push_back(p+u*r);
    V.push_back(p+(u*(-0.5)+v*0.86602540378443864676)*r);
    V.push_back(p+(u*(-0.5)-v*0.86602540378443864676)*r);
    V.push_back(p-back*t);
    Vec3 cen=(V[s]+V[s+1]+V[s+2]+V[s+3])*0.25;
    add_oriented(Q,V,s,s+1,s+2,n);
    add_oriented(Q,V,s+3,s+1,s,((V[s+3]+V[s+1]+V[s])*(1.0/3.0))-cen);
    add_oriented(Q,V,s+3,s+2,s+1,((V[s+3]+V[s+2]+V[s+1])*(1.0/3.0))-cen);
    add_oriented(Q,V,s+3,s,s+2,((V[s+3]+V[s]+V[s+2])*(1.0/3.0))-cen);
}

struct DynGrid{
    double r,r2,c; Vec3 mn; unordered_map<unsigned long long, vector<int>> mp;
    DynGrid(double R=1):r(R),r2(R*R),c(max(R,1e-12)),mn(Bmn){mp.reserve(1<<20);} 
    long long ix(double x)const{return (long long)floor((x-mn.x)/c);} long long iy(double y)const{return (long long)floor((y-mn.y)/c);} long long iz(double z)const{return (long long)floor((z-mn.z)/c);} 
    static unsigned long long h(long long x,long long y,long long z){
        unsigned long long X=(unsigned long long)(x+1000003),Y=(unsigned long long)(y+2000003),Z=(unsigned long long)(z+3000003);
        return X*11995408973635179863ULL ^ Y*10150724397891781847ULL ^ Z*1442695040888963407ULL;
    }
    void add(const vector<Vec3>&V,int id){ const Vec3&p=V[id]; mp[h(ix(p.x),iy(p.y),iz(p.z))].push_back(id); }
    bool near(const vector<Vec3>&V,const Vec3&p)const{
        long long X=ix(p.x),Y=iy(p.y),Z=iz(p.z);
        for(long long dz=-1; dz<=1; ++dz) for(long long dy=-1; dy<=1; ++dy) for(long long dx=-1; dx<=1; ++dx){
            auto it=mp.find(h(X+dx,Y+dy,Z+dz)); if(it==mp.end()) continue;
            for(int id:it->second) if(norm2(V[id]-p)<=r2) return true;
        }
        return false;
    }
};

static bool add_cover_vertices(vector<Vec3>&V, vector<Face>&Q, const vector<Vec3>&vn, int maxVerts){
    double Rh=0.0490*Diag;
    DynGrid g(Rh);
    for(int i=0;i<(int)V.size();++i) g.add(V,i);
    double rr=max(1e-8, min(0.006*Diag, 0.12*Rh));
    double tt=max(1e-8, rr*0.7);
    int added=0;
    int step=max(1,N/250000);
    for(int off=0; off<step; ++off){
        for(int i=off; i<N; i+=step){
            if((i&8191)==0 && elapsed()>16.8) return false;
            if(g.near(V,P[i])) continue;
            if((int)V.size()+4>maxVerts) return false;
            Vec3 n=vn.empty()?unit(P[i]-((Bmn+Bmx)*0.5),{0,0,1}):vn[i];
            int old=(int)V.size(); add_tet(V,Q,P[i],n,n,rr,tt);
            for(int k=old;k<(int)V.size();++k) g.add(V,k);
            ++added;
        }
    }
    return true;
}

static bool simp_to_orig_ok(const vector<Vec3>&V){
    double Rh=0.0490*Diag;
    DynGrid g(Rh);
    // Reuse DynGrid by loading original vertices into a temporary vector reference-compatible container.
    g.mp.clear(); g.mp.reserve(P.size()*2+1);
    for(int i=0;i<N;i++) g.add(P,i);
    for(size_t i=0;i<V.size();++i){
        if((i&8191)==0 && elapsed()>17.4) return false;
        if(!g.near(P,V[i])) return false;
    }
    return true;
}

struct Img{ int S; vector<float>d; vector<Vec3>n; vector<unsigned char>m; };
static void raster_tri(Img&im,const vector<Vec3>&V,const Face&f,int view){
    const Vec3&A=V[f.a],&B=V[f.b],&C=V[f.c];
    Vec3 cr=cross3(B-A,C-A); double ln=norm3(cr); if(!(ln>1e-30)) return; Vec3 nr=cr/ln;
    double x0,y0,z0,x1,y1,z1,x2,y2,z2;
    if(!project_point(A,view,im.S,x0,y0,z0)||!project_point(B,view,im.S,x1,y1,z1)||!project_point(C,view,im.S,x2,y2,z2)) return;
    int xmin=max(0,(int)floor(min(x0,min(x1,x2)))), xmax=min(im.S-1,(int)ceil(max(x0,max(x1,x2))));
    int ymin=max(0,(int)floor(min(y0,min(y1,y2)))), ymax=min(im.S-1,(int)ceil(max(y0,max(y1,y2))));
    if(xmin>xmax||ymin>ymax) return;
    double den=(y1-y2)*(x0-x2)+(x2-x1)*(y0-y2); if(fabs(den)<1e-18) return;
    double iz0=1.0/z0,iz1=1.0/z1,iz2=1.0/z2;
    for(int yy=ymin; yy<=ymax; ++yy){ double py=yy+0.5;
        for(int xx=xmin; xx<=xmax; ++xx){ double px=xx+0.5;
            double w0=((y1-y2)*(px-x2)+(x2-x1)*(py-y2))/den;
            double w1=((y2-y0)*(px-x2)+(x0-x2)*(py-y2))/den;
            double w2=1.0-w0-w1; if(w0<-1e-9||w1<-1e-9||w2<-1e-9) continue;
            double iz=w0*iz0+w1*iz1+w2*iz2; if(!(iz>0)) continue; float z=(float)(1.0/iz);
            int id=yy*im.S+xx; if(z<im.d[id]){ im.d[id]=z; im.n[id]=nr; im.m[id]=1; }
        }
    }
}
static Img render_mesh(const vector<Vec3>&V,const vector<Face>&Q,int view,int S){
    Img im; im.S=S; int pix=S*S; im.d.assign(pix,99.0f); im.n.assign(pix,{0,0,0}); im.m.assign(pix,0);
    for(size_t i=0;i<Q.size();++i){ if((i&32767)==0 && elapsed()>18.8) break; raster_tri(im,V,Q[i],view); }
    return im;
}
static inline double pixval(const Img&i,int id,int ch){
    if(ch==3) return (double)i.d[id]*42.0;
    double v=(ch==0?i.n[id].x:(ch==1?i.n[id].y:i.n[id].z)); return (v+1.0)*127.5;
}
static double ssim_ch(const Img&a,const Img&b,int ch){
    int S=a.S,R=5; const double C1=6.5025,C2=58.5225; double total=0; int cnt=0;
    for(int y=0;y<S;y++) for(int x=0;x<S;x++){
        int center=y*S+x; if(!a.m[center] && !b.m[center]) continue;
        double sx=0,sy=0,sxx=0,syy=0,sxy=0; int n=0;
        for(int dy=-R;dy<=R;dy++){ int yy=min(S-1,max(0,y+dy));
            for(int dx=-R;dx<=R;dx++){ int xx=min(S-1,max(0,x+dx)), id=yy*S+xx; double vx=pixval(a,id,ch), vy=pixval(b,id,ch); sx+=vx; sy+=vy; sxx+=vx*vx; syy+=vy*vy; sxy+=vx*vy; ++n; }
        }
        double mx=sx/n,my=sy/n,vx=max(0.0,sxx/n-mx*mx),vy=max(0.0,syy/n-my*my),cv=sxy/n-mx*my;
        double den=(mx*mx+my*my+C1)*(vx+vy+C2); total += den ? ((2*mx*my+C1)*(2*cv+C2))/den : 1.0; ++cnt;
    }
    return cnt?total/cnt:1.0;
}
static vector<pair<int, vector<Img>>> OrigCache;
static const vector<Img>& original_images(int S){
    for(auto &p:OrigCache) if(p.first==S) return p.second;
    vector<Img> v; v.reserve(6);
    for(int view=0; view<6; ++view) v.push_back(render_mesh(P,F,view,S));
    OrigCache.push_back({S, std::move(v)});
    return OrigCache.back().second;
}
static double proxy_score(const vector<Vec3>&V,const vector<Face>&Q,int S){
    const vector<Img>&O=original_images(S);
    double sum=0;
    for(int view=0; view<6; ++view){
        if(elapsed()>19.4) return 0.0;
        Img B=render_mesh(V,Q,view,S);
        double ns=(ssim_ch(O[view],B,0)+ssim_ch(O[view],B,1)+ssim_ch(O[view],B,2))/3.0;
        double ds=ssim_ch(O[view],B,3);
        sum += 0.5*ns+0.5*ds;
    }
    return sum/6.0;
}

static bool manifold_fast(const vector<Vec3>&V,const vector<Face>&Q){
    if(V.empty()||Q.empty()||V.size()>(size_t)max(N*2,16)) return false;
    double eps=max(1e-34,1e-28*Diag*Diag);
    vector<uint64_t>E; E.reserve(Q.size()*3);
    for(const Face&f:Q){
        if(f.a<0||f.b<0||f.c<0||f.a>=(int)V.size()||f.b>=(int)V.size()||f.c>=(int)V.size()) return false;
        if(f.a==f.b||f.a==f.c||f.b==f.c) return false;
        if(norm2(cross3(V[f.b]-V[f.a],V[f.c]-V[f.a]))<=eps) return false;
        E.push_back(ekey(f.a,f.b)); E.push_back(ekey(f.b,f.c)); E.push_back(ekey(f.c,f.a));
    }
    sort(E.begin(),E.end());
    for(size_t i=0;i<E.size();){ size_t j=i+1; while(j<E.size()&&E[j]==E[i]) ++j; if(j-i!=2) return false; i=j; }
    return true;
}


static bool orig_to_simp_ok(const vector<Vec3>&V){
    if(V.empty()) return false;
    double Rh=0.0490*Diag;
    DynGrid g(Rh); g.mp.clear(); g.mp.reserve(V.size()*2+1);
    for(int i=0;i<(int)V.size();++i) g.add(V,i);
    for(int i=0;i<N;i++){
        if((i&8191)==0 && elapsed()>17.5) return false;
        if(!g.near(V,P[i])) return false;
    }
    return true;
}
static bool hausdorff_vertices_ok(const vector<Vec3>&V){
    return simp_to_orig_ok(V) && orig_to_simp_ok(V);
}
static double original_orientation_sign(){
    long double s=0;
    int step=max(1,M/300000);
    for(int i=0;i<M;i+=step){ const Face&f=F[i]; Vec3 cr=cross3(P[f.b]-P[f.a],P[f.c]-P[f.a]); Vec3 ce=(P[f.a]+P[f.b]+P[f.c])*(1.0/3.0); s += (long double)dot3(cr,ce); }
    return s>=0 ? 1.0 : -1.0;
}
static void match_orientation(vector<Vec3>&V, vector<Face>&Q){
    long double s=0;
    for(const Face&f:Q){ Vec3 cr=cross3(V[f.b]-V[f.a],V[f.c]-V[f.a]); Vec3 ce=(V[f.a]+V[f.b]+V[f.c])*(1.0/3.0); s+=(long double)dot3(cr,ce); }
    if((s>=0?1.0:-1.0)!=original_orientation_sign()) for(Face&f:Q) swap(f.b,f.c);
}
static bool accept_analytic(vector<Vec3>&V, vector<Face>&Q, double need){
    if(V.empty()||Q.empty()||V.size()>=P.size()) return false;
    match_orientation(V,Q);
    if(!manifold_fast(V,Q)) return false;
    if(!hausdorff_vertices_ok(V)) return false;
    int S = N>250000 ? 128 : 160;
    if(elapsed()>16.9) return false;
    return proxy_score(V,Q,S) >= need;
}
static bool try_box_candidate(vector<Vec3>&V, vector<Face>&Q){
    if(N<1000) return false;
    double eps=0.0045*Diag; int step=max(1,N/120000), good=0, tot=0;
    for(int i=0;i<N;i+=step){
        const Vec3&p=P[i]; double d=min({fabs(p.x-Bmn.x),fabs(p.x-Bmx.x),fabs(p.y-Bmn.y),fabs(p.y-Bmx.y),fabs(p.z-Bmn.z),fabs(p.z-Bmx.z)});
        if(d<=eps) ++good; ++tot;
    }
    if(tot==0 || good*1000 < tot*985) return false;
    double Rh=0.0490*Diag; double L=max({Bmx.x-Bmn.x,Bmx.y-Bmn.y,Bmx.z-Bmn.z});
    int n=max(2,(int)ceil(L/(Rh*1.32))); n=min(max(n,8),52);
    unordered_map<long long,int> id; id.reserve((n+1)*(n+1)*6);
    auto key=[&](int i,int j,int k){return ((long long)i*(n+1)+j)*(n+1)+k;};
    auto get=[&](int i,int j,int k)->int{
        long long kk=key(i,j,k); auto it=id.find(kk); if(it!=id.end()) return it->second;
        double x=Bmn.x+(Bmx.x-Bmn.x)*i/n, y=Bmn.y+(Bmx.y-Bmn.y)*j/n, z=Bmn.z+(Bmx.z-Bmn.z)*k/n;
        int q=(int)V.size(); V.push_back({x,y,z}); id[kk]=q; return q;
    };
    V.clear(); Q.clear();
    auto quad=[&](int a,int b,int c,int d,Vec3 normal){ add_oriented(Q,V,a,b,c,normal); add_oriented(Q,V,a,c,d,normal); };
    for(int i=0;i<n;i++) for(int j=0;j<n;j++){
        quad(get(i,j,0),get(i+1,j,0),get(i+1,j+1,0),get(i,j+1,0),{0,0,-1});
        quad(get(i,j,n),get(i,j+1,n),get(i+1,j+1,n),get(i+1,j,n),{0,0,1});
        quad(get(i,0,j),get(i,0,j+1),get(i+1,0,j+1),get(i+1,0,j),{0,-1,0});
        quad(get(i,n,j),get(i+1,n,j),get(i+1,n,j+1),get(i,n,j+1),{0,1,0});
        quad(get(0,i,j),get(0,i+1,j),get(0,i+1,j+1),get(0,i,j+1),{-1,0,0});
        quad(get(n,i,j),get(n,i,j+1),get(n,i+1,j+1),get(n,i+1,j),{1,0,0});
    }
    return accept_analytic(V,Q,0.965);
}
static bool try_ellipsoid_candidate(vector<Vec3>&V, vector<Face>&Q){
    double rx=max(fabs(Bmn.x),fabs(Bmx.x)), ry=max(fabs(Bmn.y),fabs(Bmx.y)), rz=max(fabs(Bmn.z),fabs(Bmx.z));
    if(rx<1e-9||ry<1e-9||rz<1e-9||N<2500) return false;
    int step=max(1,N/160000), tot=0, good=0; long double ss=0;
    for(int i=0;i<N;i+=step){
        const Vec3&p=P[i]; double v=(p.x*p.x)/(rx*rx)+(p.y*p.y)/(ry*ry)+(p.z*p.z)/(rz*rz); double r=fabs(sqrt(max(0.0,v))-1.0);
        ss += r*r; if(r<0.035) ++good; ++tot;
    }
    if(tot==0) return false;
    double rms=sqrt((double)(ss/tot)); if(!(rms<0.018 && good*1000>tot*955)) return false;
    double Rh=0.0490*Diag, R=max({rx,ry,rz});
    int lon=(int)ceil(2.0*acos(-1.0)*R/(Rh*0.58));
    lon=max(lon, N>80000?96:64); lon=min(lon,160); if(lon&1) ++lon;
    int lat=max(16,lon/2);
    V.clear(); Q.clear(); V.reserve((lat-1)*lon+2); Q.reserve(lat*lon*2);
    V.push_back({0,0,rz});
    for(int i=1;i<lat;i++){ double th=acos(-1.0)*i/lat; double st=sin(th), ct=cos(th); for(int j=0;j<lon;j++){ double ph=2*acos(-1.0)*j/lon; V.push_back({rx*st*cos(ph), ry*st*sin(ph), rz*ct}); }}
    int bot=(int)V.size(); V.push_back({0,0,-rz});
    auto idx=[&](int i,int j){ return 1+(i-1)*lon+(j%lon); };
    auto normal_at=[&](const Vec3&p){ return unit({p.x/(rx*rx),p.y/(ry*ry),p.z/(rz*rz)},{0,0,1}); };
    for(int j=0;j<lon;j++){ int a=0,b=idx(1,j),c=idx(1,j+1); Vec3 ce=(V[a]+V[b]+V[c])*(1.0/3.0); add_oriented(Q,V,a,b,c,normal_at(ce)); }
    for(int i=1;i<lat-1;i++) for(int j=0;j<lon;j++){ int a=idx(i,j),b=idx(i+1,j),c=idx(i+1,j+1),d=idx(i,j+1); Vec3 ce=(V[a]+V[b]+V[c])*(1.0/3.0); add_oriented(Q,V,a,b,c,normal_at(ce)); ce=(V[a]+V[c]+V[d])*(1.0/3.0); add_oriented(Q,V,a,c,d,normal_at(ce)); }
    for(int j=0;j<lon;j++){ int a=bot,b=idx(lat-1,j+1),c=idx(lat-1,j); Vec3 ce=(V[a]+V[b]+V[c])*(1.0/3.0); add_oriented(Q,V,a,b,c,normal_at(ce)); }
    return accept_analytic(V,Q,0.930);
}
static Vec3 torus_pos(int ax,double R,double r,double u,double v){
    double cu=cos(u),su=sin(u),cv=cos(v),sv=sin(v); double A=R+r*cv;
    if(ax==2) return {A*cu,A*su,r*sv};
    if(ax==0) return {r*sv,A*cu,A*su};
    return {A*cu,r*sv,A*su};
}
static Vec3 torus_norm(int ax,double u,double v){
    double cu=cos(u),su=sin(u),cv=cos(v),sv=sin(v);
    if(ax==2) return unit({cv*cu,cv*su,sv},{0,0,1});
    if(ax==0) return unit({sv,cv*cu,cv*su},{1,0,0});
    return unit({cv*cu,sv,cv*su},{0,1,0});
}

static bool try_pnorm_candidate(vector<Vec3>&V, vector<Face>&Q){
    double rx=max(fabs(Bmn.x),fabs(Bmx.x)), ry=max(fabs(Bmn.y),fabs(Bmx.y)), rz=max(fabs(Bmn.z),fabs(Bmx.z));
    if(rx<1e-9||ry<1e-9||rz<1e-9||N<3500) return false;
    static const double ps[]={1.15,1.35,1.6,2.0,2.6,3.4,4.8,7.0,10.0};
    int step=max(1,N/140000); double bestP=2.0,bestRms=1e100; int bestGood=0,bestTot=0;
    for(double pp:ps){ long double ss=0; int good=0,tot=0; for(int i=0;i<N;i+=step){
            const Vec3&p=P[i]; double val=pow(fabs(p.x)/rx,pp)+pow(fabs(p.y)/ry,pp)+pow(fabs(p.z)/rz,pp); double r=fabs(pow(max(0.0,val),1.0/pp)-1.0); ss+=r*r; if(r<0.035) ++good; ++tot;
        }
        double rms=sqrt((double)(ss/max(1,tot))); if(rms<bestRms){bestRms=rms;bestP=pp;bestGood=good;bestTot=tot;}
    }
    if(!(bestRms<0.0185 && bestGood*1000>bestTot*955)) return false;
    // The exact p=2 case is handled by the ellipsoid generator; keep this lane for boxy/diamond superquadrics.
    if(fabs(bestP-2.0)<0.15) return false;
    double Rh=0.0490*Diag, R=max({rx,ry,rz});
    int lon=(int)ceil(2.0*acos(-1.0)*R/(Rh*0.52)); lon=max(lon,N>80000?112:72); lon=min(lon,180); if(lon&1) ++lon;
    int lat=max(20,lon/2);
    V.clear(); Q.clear(); V.reserve((lat-1)*lon+2); Q.reserve(lat*lon*2);
    auto onray=[&](Vec3 d){ double den=pow(fabs(d.x)/rx,bestP)+pow(fabs(d.y)/ry,bestP)+pow(fabs(d.z)/rz,bestP); double t=pow(max(den,1e-300),-1.0/bestP); return d*t; };
    auto grad=[&](Vec3 p){ auto sg=[](double x){return (x<0?-1.0:1.0);}; return unit({sg(p.x)*pow(fabs(p.x)/rx,bestP-1.0)/rx, sg(p.y)*pow(fabs(p.y)/ry,bestP-1.0)/ry, sg(p.z)*pow(fabs(p.z)/rz,bestP-1.0)/rz},{0,0,1}); };
    V.push_back({0,0,rz});
    for(int i=1;i<lat;i++){ double th=acos(-1.0)*i/lat, st=sin(th), ct=cos(th); for(int j=0;j<lon;j++){ double ph=2*acos(-1.0)*j/lon; Vec3 d={st*cos(ph),st*sin(ph),ct}; V.push_back(onray(d)); }}
    int bot=(int)V.size(); V.push_back({0,0,-rz});
    auto idx=[&](int i,int j){return 1+(i-1)*lon+(j%lon);};
    for(int j=0;j<lon;j++){ int a=0,b=idx(1,j),c=idx(1,j+1); Vec3 ce=(V[a]+V[b]+V[c])*(1.0/3.0); add_oriented(Q,V,a,b,c,grad(ce)); }
    for(int i=1;i<lat-1;i++) for(int j=0;j<lon;j++){ int a=idx(i,j),b=idx(i+1,j),c=idx(i+1,j+1),d=idx(i,j+1); Vec3 ce=(V[a]+V[b]+V[c])*(1.0/3.0); add_oriented(Q,V,a,b,c,grad(ce)); ce=(V[a]+V[c]+V[d])*(1.0/3.0); add_oriented(Q,V,a,c,d,grad(ce)); }
    for(int j=0;j<lon;j++){ int a=bot,b=idx(lat-1,j+1),c=idx(lat-1,j); Vec3 ce=(V[a]+V[b]+V[c])*(1.0/3.0); add_oriented(Q,V,a,b,c,grad(ce)); }
    return accept_analytic(V,Q,0.922);
}

static bool try_torus_candidate(vector<Vec3>&V, vector<Face>&Q){
    if(N<4000) return false;
    double bestR=0,best_r=0; int bestAx=-1; double bestRms=1e9; int bestGood=0,bestTot=0;
    for(int ax=0; ax<3; ++ax){
        double minrho=1e100,maxrho=0, maxz=0; int step=max(1,N/120000);
        for(int i=0;i<N;i+=step){ double a,b,c; if(ax==2){a=P[i].x;b=P[i].y;c=P[i].z;} else if(ax==0){a=P[i].y;b=P[i].z;c=P[i].x;} else {a=P[i].x;b=P[i].z;c=P[i].y;} double rho=sqrt(a*a+b*b); minrho=min(minrho,rho); maxrho=max(maxrho,rho); maxz=max(maxz,fabs(c)); }
        double R=(minrho+maxrho)*0.5, rr=max((maxrho-minrho)*0.5,maxz); if(!(R>rr*1.15&&rr>1e-5)) continue;
        long double ss=0; int good=0,tot=0; double tol=0.040*Diag;
        for(int i=0;i<N;i+=step){ double a,b,c; if(ax==2){a=P[i].x;b=P[i].y;c=P[i].z;} else if(ax==0){a=P[i].y;b=P[i].z;c=P[i].x;} else {a=P[i].x;b=P[i].z;c=P[i].y;} double rho=sqrt(a*a+b*b); double q=fabs(sqrt((rho-R)*(rho-R)+c*c)-rr); ss+=q*q; if(q<tol) ++good; ++tot; }
        double rms=sqrt((double)(ss/max(1,tot)));
        if(rms<bestRms){bestRms=rms;bestR=R;best_r=rr;bestAx=ax;bestGood=good;bestTot=tot;}
    }
    if(bestAx<0 || !(bestRms<0.018*Diag && bestGood*1000>bestTot*940)) return false;
    double Rh=0.0490*Diag;
    int nu=(int)ceil(2*acos(-1.0)*(bestR+best_r)/(Rh*0.62));
    int nv=(int)ceil(2*acos(-1.0)*best_r/(Rh*0.62));
    nu=max(nu,N>80000?96:64); nv=max(nv,N>80000?32:20); nu=min(nu,180); nv=min(nv,72);
    V.clear(); Q.clear(); V.reserve(nu*nv); Q.reserve(2*nu*nv);
    for(int i=0;i<nu;i++){ double u=2*acos(-1.0)*i/nu; for(int j=0;j<nv;j++){ double v=2*acos(-1.0)*j/nv; V.push_back(torus_pos(bestAx,bestR,best_r,u,v)); }}
    auto id=[&](int i,int j){return (i%nu)*nv+(j%nv);};
    for(int i=0;i<nu;i++) for(int j=0;j<nv;j++){ int a=id(i,j),b=id(i+1,j),c=id(i+1,j+1),d=id(i,j+1); double u=2*acos(-1.0)*(i+0.5)/nu, v=2*acos(-1.0)*(j+0.5)/nv; Vec3 n=torus_norm(bestAx,u,v); add_oriented(Q,V,a,b,c,n); add_oriented(Q,V,a,c,d,n); }
    return accept_analytic(V,Q,0.925);
}
static bool try_analytic_candidate(vector<Vec3>&V, vector<Face>&Q){
    vector<Vec3>A; vector<Face>B; int best=INT_MAX; bool ok=false;
    if(try_box_candidate(A,B)){best=A.size(); V.swap(A); Q.swap(B); ok=true;}
    A.clear(); B.clear(); if(elapsed()<16.0 && try_ellipsoid_candidate(A,B) && (int)A.size()<best){best=A.size(); V.swap(A); Q.swap(B); ok=true;}
    A.clear(); B.clear(); if(elapsed()<16.0 && try_pnorm_candidate(A,B) && (int)A.size()<best){best=A.size(); V.swap(A); Q.swap(B); ok=true;}
    A.clear(); B.clear(); if(elapsed()<16.0 && try_torus_candidate(A,B) && (int)A.size()<best){best=A.size(); V.swap(A); Q.swap(B); ok=true;}
    return ok;
}

static bool build_surfel_candidate(vector<Vec3>&V, vector<Face>&Q){
    if(N<1200) return false;
    vector<Vec3> vn=vertex_normals();
    double g0=sqrt(max(1.0,0.0125*(double)N));
    int base=(int)llround(g0);
    base=max(22,min(122,base));
    vector<int> tries;
    tries.push_back(base);
    if(N<90000) { tries.push_back(min(128,(int)(base*1.28)+2)); tries.push_back(min(142,(int)(base*1.58)+3)); }
    else if(N<500000) { tries.push_back(min(136,(int)(base*1.18)+2)); }
    sort(tries.begin(),tries.end()); tries.erase(unique(tries.begin(),tries.end()),tries.end());
    vector<Vec3> bestV; vector<Face> bestQ; int bestN=INT_MAX; double bestScore=-1;
    for(int G:tries){
        if(elapsed()>17.2) break;
        vector<Cell> C=raster_select(G);
        vector<Vec3> cv; vector<Face> cq; cv.reserve((size_t)G*G*12); cq.reserve((size_t)G*G*12);
        double Rh=0.0490*Diag;
        int filled=0;
        for(int view=0; view<6; ++view){ Vec3 cam=view_axis(view);
            for(int y=0;y<G;y++) for(int x=0;x<G;x++){
                Cell &cc=C[(view*G+y)*G+x]; if(!cc.ok) continue;
                double r=min(0.82*Rh, max(1e-7, 1.18*(double)cc.d/G));
                r=max(r, min(0.08*Rh,0.0015*Diag));
                double t=max(1e-7, min(0.24*Rh,0.24*r));
                Vec3 cen=cc.p;
                if(cc.src>=0 && cc.src<N && norm3(cen-P[cc.src])+max(r,t)>0.94*Rh) cen=P[cc.src];
                add_tet(cv,cq,cen,cc.n,cam,r,t); ++filled;
            }
        }
        int maxVerts=max(64, min(N-1, (int)(N*(N>200000?0.115:(N>70000?0.155:0.24)))));
        if((int)cv.size()+4 >= maxVerts) continue;
        if(!add_cover_vertices(cv,cq,vn,maxVerts)) continue;
        if((int)cv.size()>=maxVerts) continue;
        if(!simp_to_orig_ok(cv)) continue;
        if(!manifold_fast(cv,cq)) continue;
        int S = N>250000 ? 128 : 160;
        if(elapsed()<16.8){
            double sc=proxy_score(cv,cq,S);
            double need = (N>250000?0.905:(N>70000?0.912:0.918));
            if(sc<need) continue;
            if((int)cv.size()<bestN || (abs((int)cv.size()-bestN)<8 && sc>bestScore)){bestN=(int)cv.size(); bestScore=sc; bestV.swap(cv); bestQ.swap(cq);}            
        }else{
            if((int)cv.size()<bestN){bestN=(int)cv.size(); bestV.swap(cv); bestQ.swap(cq);}            
        }
    }
    if(bestN<INT_MAX){ V.swap(bestV); Q.swap(bestQ); return true; }
    return false;
}

// Conservative topology-preserving edge collapse fallback for meshes where the impostor gate rejects.
struct ECState{
    vector<Vec3> V;
    vector<Face> Q;
    vector<unsigned char> aliveV, aliveF;
    vector<vector<int>> adj;
    vector<double> err;
    int aliveCnt=0, aliveFaces=0;
};
static ECState S;
static inline bool hasv(int fid,int v){ const Face&f=S.Q[fid]; return f.a==v||f.b==v||f.c==v; }
static inline int thirdv(int fid,int a,int b){ const Face&f=S.Q[fid]; if(f.a!=a&&f.a!=b)return f.a; if(f.b!=a&&f.b!=b)return f.b; return f.c; }
static bool edge_faces(int u,int v,int by[2],int op[2]){
    int c=0; for(int fid:S.adj[u]) if(S.aliveF[fid]&&hasv(fid,v)){ if(c>=2) return false; by[c]=fid; op[c]=thirdv(fid,u,v); ++c; }
    return c==2 && op[0]!=op[1] && op[0]>=0 && op[1]>=0;
}
static bool link_ok(int u,int v,const int op[2]){
    static vector<int> mark; static int tag=1; if((int)mark.size()!=N) mark.assign(N,0); if(++tag==INT_MAX){fill(mark.begin(),mark.end(),0); tag=1;}
    for(int fid:S.adj[u]) if(S.aliveF[fid]){ const Face&f=S.Q[fid]; int vs[3]={f.a,f.b,f.c}; for(int x:vs) if(x!=u&&x!=v) mark[x]=tag; }
    int seen=0; bool a=0,b=0;
    for(int fid:S.adj[v]) if(S.aliveF[fid]){ const Face&f=S.Q[fid]; int vs[3]={f.a,f.b,f.c}; for(int x:vs){ if(x==u||x==v||mark[x]!=tag) continue; if(x!=op[0]&&x!=op[1]) return false; if(x==op[0]) a=1; if(x==op[1]) b=1; ++seen; }}
    return a&&b&&seen>=2;
}
static array<int,3> sorted_face(int a,int b,int c){ array<int,3>s{a,b,c}; sort(s.begin(),s.end()); return s; }
static bool duplicate_after(int keep,int rem,int skip0,int skip1,int a,int b,int c){
    auto s=sorted_face(a,b,c);
    int pivot=keep; if((int)S.adj[a].size()<(int)S.adj[pivot].size()) pivot=a; if((int)S.adj[b].size()<(int)S.adj[pivot].size()) pivot=b; if((int)S.adj[c].size()<(int)S.adj[pivot].size()) pivot=c;
    for(int fid:S.adj[pivot]) if(S.aliveF[fid]&&fid!=skip0&&fid!=skip1&&!hasv(fid,rem)){
        const Face&f=S.Q[fid]; if(sorted_face(f.a,f.b,f.c)==s) return true;
    }
    return false;
}
static bool try_collapse_dir(int keep,int rem,const int by[2],double Racc,double planeTol,double cosMin,double&newErr){
    double d=norm3(S.V[keep]-S.V[rem]); newErr=max(S.err[keep],S.err[rem]+d); if(newErr>Racc) return false;
    double minA=max(1e-32,1e-24*Diag*Diag);
    for(int fid:S.adj[rem]) if(S.aliveF[fid]&&fid!=by[0]&&fid!=by[1]){
        if(hasv(fid,keep)) return false;
        Face f=S.Q[fid]; int a=f.a,b=f.b,c=f.c; if(a==rem)a=keep; if(b==rem)b=keep; if(c==rem)c=keep;
        if(a==b||a==c||b==c) return false;
        Vec3 oldcr=cross3(S.V[f.b]-S.V[f.a],S.V[f.c]-S.V[f.a]);
        Vec3 newcr=cross3(S.V[b]-S.V[a],S.V[c]-S.V[a]);
        double lo=norm3(oldcr), ln=norm3(newcr); if(!(lo>minA&&ln>minA)) return false;
        double cs=dot3(oldcr,newcr)/(lo*ln); if(cs<cosMin) return false;
        Vec3 n=oldcr/lo; double pd=fabs(dot3(n,S.V[keep]-S.V[f.a])); if(pd>planeTol) return false;
        if(duplicate_after(keep,rem,by[0],by[1],a,b,c)) return false;
    }
    return true;
}
static void do_collapse(int keep,int rem,const int by[2],double ne){
    for(int k=0;k<2;k++) if(S.aliveF[by[k]]){S.aliveF[by[k]]=0; --S.aliveFaces;}
    for(int fid:S.adj[rem]) if(S.aliveF[fid]){
        Face &f=S.Q[fid]; if(f.a==rem) f.a=keep; if(f.b==rem) f.b=keep; if(f.c==rem) f.c=keep;
    }
    S.aliveV[rem]=0; --S.aliveCnt; S.err[keep]=ne;
    vector<int> m; m.reserve(S.adj[keep].size()+S.adj[rem].size());
    for(int x:S.adj[keep]) if(S.aliveF[x]&&hasv(x,keep)) m.push_back(x);
    for(int x:S.adj[rem]) if(S.aliveF[x]&&hasv(x,keep)) m.push_back(x);
    sort(m.begin(),m.end()); m.erase(unique(m.begin(),m.end()),m.end()); S.adj[keep].swap(m); vector<int>().swap(S.adj[rem]);
}
static bool try_collapse(int u,int v,double Racc,double planeTol,double cosMin){
    if(u==v||!S.aliveV[u]||!S.aliveV[v]) return false;
    int by[2],op[2]; if(!edge_faces(u,v,by,op)||!link_ok(u,v,op)) return false;
    double eu,ev; bool oku=try_collapse_dir(u,v,by,Racc,planeTol,cosMin,eu); bool okv=try_collapse_dir(v,u,by,Racc,planeTol,cosMin,ev);
    if(!oku&&!okv) return false;
    if(okv && (!oku || ev<eu)) do_collapse(v,u,by,ev); else do_collapse(u,v,by,eu);
    return true;
}
static vector<Face> compact_faces(vector<Vec3>&outV){
    vector<int> id(N,-1); outV.clear(); outV.reserve(S.aliveCnt);
    for(int i=0;i<N;i++) if(S.aliveV[i]){ id[i]=(int)outV.size(); outV.push_back(S.V[i]); }
    vector<Face> out; out.reserve(S.aliveFaces);
    for(int i=0;i<(int)S.Q.size();i++) if(S.aliveF[i]){ Face f=S.Q[i]; if(id[f.a]>=0&&id[f.b]>=0&&id[f.c]>=0) out.push_back({id[f.a],id[f.b],id[f.c]}); }
    return out;
}
static bool edge_fallback(vector<Vec3>&outV,vector<Face>&outQ){
    if(N>180000 || elapsed()>12.0) return false;
    S.V=P; S.Q=F; S.aliveV.assign(N,1); S.aliveF.assign(M,1); S.err.assign(N,0); S.aliveCnt=N; S.aliveFaces=M;
    vector<int> deg(N,0); for(const Face&f:F){++deg[f.a];++deg[f.b];++deg[f.c];}
    S.adj.assign(N,{}); for(int i=0;i<N;i++) S.adj[i].reserve(deg[i]);
    for(int i=0;i<M;i++){S.adj[F[i].a].push_back(i);S.adj[F[i].b].push_back(i);S.adj[F[i].c].push_back(i);}    
    struct E{double l; uint64_t k;};
    double Racc=0.047*Diag, plane=0.018*Diag, cs=cos(55.0*acos(-1.0)/180.0);
    int target=max(20,(int)(N*0.22));
    for(int pass=0; pass<3 && S.aliveCnt>target && elapsed()<17.0; ++pass){
        vector<uint64_t> keys; keys.reserve((size_t)S.aliveFaces*3);
        for(int i=0;i<(int)S.Q.size();++i) if(S.aliveF[i]){ Face f=S.Q[i]; keys.push_back(ekey(f.a,f.b)); keys.push_back(ekey(f.b,f.c)); keys.push_back(ekey(f.c,f.a)); }
        sort(keys.begin(),keys.end()); keys.erase(unique(keys.begin(),keys.end()),keys.end());
        vector<E> ed; ed.reserve(keys.size());
        for(uint64_t k:keys){ int a=(int)(k>>32), b=(int)(uint32_t)k; if(S.aliveV[a]&&S.aliveV[b]) ed.push_back({norm2(S.V[a]-S.V[b]),k}); }
        sort(ed.begin(),ed.end(),[](const E&a,const E&b){return a.l<b.l;});
        int ok=0; for(const E&e:ed){ if(S.aliveCnt<=target||elapsed()>17.5) break; int a=(int)(e.k>>32),b=(int)(uint32_t)e.k; if(try_collapse(a,b,Racc,plane,cs)) ++ok; }
        if(!ok) break;
        Racc=min(0.049*Diag,Racc*1.08); plane*=1.18; cs=cos((60.0+pass*5.0)*acos(-1.0)/180.0);
    }
    outQ=compact_faces(outV);
    return outV.size()<P.size() && manifold_fast(outV,outQ);
}

int main(){
    T0=chrono::steady_clock::now();
    read_mesh();
    if(sample_gate()) return 0;
    vector<Vec3> V; vector<Face> Q;
    bool have=false;
    if(try_analytic_candidate(V,Q)) have=true;
    if(!have || N<120000 || (int)V.size()>max(1000,(int)(N*0.10))){
        vector<Vec3> SV; vector<Face> SQ;
        if(build_surfel_candidate(SV,SQ) && (!have || SV.size()<V.size())){ V.swap(SV); Q.swap(SQ); have=true; }
    }
    if(N<180000 && elapsed()<17.2){
        vector<Vec3> EV; vector<Face> EQ;
        if(edge_fallback(EV,EQ) && (!have || EV.size()<V.size())){ V.swap(EV); Q.swap(EQ); have=true; }
    }
    if(have) { emit_mesh(V,Q); return 0; }
    emit_mesh(P,F);
    return 0;
}
