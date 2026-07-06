#include <bits/stdc++.h>
using namespace std;

struct Vec3{
    double x,y,z;
    Vec3():x(0),y(0),z(0){}
    Vec3(double x_,double y_,double z_):x(x_),y(y_),z(z_){}
    Vec3 operator+(const Vec3& o) const {return Vec3(x+o.x,y+o.y,z+o.z);}
    Vec3 operator-(const Vec3& o) const {return Vec3(x-o.x,y-o.y,z-o.z);}
    Vec3 operator*(double s) const {return Vec3(x*s,y*s,z*s);}
    Vec3 operator/(double s) const {return Vec3(x/s,y/s,z/s);}
    Vec3& operator+=(const Vec3& o){x+=o.x;y+=o.y;z+=o.z;return *this;}
};
static inline double dotv(const Vec3&a,const Vec3&b){return a.x*b.x+a.y*b.y+a.z*b.z;}
static inline Vec3 crossv(const Vec3&a,const Vec3&b){return Vec3(a.y*b.z-a.z*b.y,a.z*b.x-a.x*b.z,a.x*b.y-a.y*b.x);}
static inline double norm2(const Vec3&a){return dotv(a,a);}
static inline double normv(const Vec3&a){return sqrt(norm2(a));}

struct Quadric{
    double q[10];
    Quadric(){ memset(q,0,sizeof(q)); }
    void clear(){ memset(q,0,sizeof(q)); }
    void addPlane(double a,double b,double c,double d,double w){
        q[0]+=w*a*a; q[1]+=w*a*b; q[2]+=w*a*c; q[3]+=w*a*d;
        q[4]+=w*b*b; q[5]+=w*b*c; q[6]+=w*b*d;
        q[7]+=w*c*c; q[8]+=w*c*d; q[9]+=w*d*d;
    }
    Quadric operator+(const Quadric&o) const { Quadric r; for(int i=0;i<10;i++) r.q[i]=q[i]+o.q[i]; return r; }
    Quadric& operator+=(const Quadric&o){ for(int i=0;i<10;i++) q[i]+=o.q[i]; return *this; }
    double eval(const Vec3&p) const{
        double x=p.x,y=p.y,z=p.z;
        return q[0]*x*x + 2*q[1]*x*y + 2*q[2]*x*z + 2*q[3]*x
             + q[4]*y*y + 2*q[5]*y*z + 2*q[6]*y
             + q[7]*z*z + 2*q[8]*z + q[9];
    }
};

struct Face{
    int a,b,c;
    Vec3 n;
    bool active;
};
struct Vertex{
    Vec3 p;
    Vec3 sum;
    Vec3 sc;
    double sr;
    int cnt;
    int ver;
    bool active;
    Quadric q;
    vector<int> inc;
};

static vector<Vertex> V;
static vector<Face> F;
static int NV, NF;
static int activeVerts;
static double hausTol, hausTolSafe;
static vector<int> vmark, fmark;
static int vstamp=1, fstamp=1;
static chrono::steady_clock::time_point startTime;

static vector<char> slurp_stdin(){
    vector<char> buf;
    buf.reserve(1<<27);
    char tmp[1<<16]; size_t n;
    while((n=fread(tmp,1,sizeof(tmp),stdin))>0) buf.insert(buf.end(),tmp,tmp+n);
    buf.push_back(0); return buf;
}
static void load_obj(){
    vector<char> buf=slurp_stdin(); char *p=buf.data();
    long nv=strtol(p,&p,10), nf=strtol(p,&p,10);
    NV=(int)nv; NF=(int)nf; V.resize(NV); F.resize(NF);
    for(int i=0;i<NV;i++){
        while(*p==' '||*p=='\n'||*p=='\r'||*p=='\t') ++p;
        if(*p=='v') ++p;
        double x=strtod(p,&p), y=strtod(p,&p), z=strtod(p,&p);
        V[i].p=Vec3(x,y,z); V[i].sum=V[i].p; V[i].sc=V[i].p; V[i].sr=0.0; V[i].cnt=1; V[i].ver=1; V[i].active=true; V[i].q.clear();
    }
    for(int i=0;i<NF;i++){
        while(*p==' '||*p=='\n'||*p=='\r'||*p=='\t') ++p;
        if(*p=='f') ++p;
        int a=(int)strtol(p,&p,10)-1, b=(int)strtol(p,&p,10)-1, c=(int)strtol(p,&p,10)-1;
        F[i].a=a; F[i].b=b; F[i].c=c; F[i].active=true;
        V[a].inc.push_back(i); V[b].inc.push_back(i); V[c].inc.push_back(i);
    }
    activeVerts=NV;
}
static bool recomputeFace(int fid){
    Face &f=F[fid];
    Vec3 p0=V[f.a].p, p1=V[f.b].p, p2=V[f.c].p;
    Vec3 cr=crossv(p1-p0,p2-p0); double l=normv(cr);
    if(l<=1e-300){ f.n=Vec3(0,0,0); return false; }
    f.n=cr/l; return true;
}

static double complexityP95 = 0.0;
static double edgeP90 = 0.0;
static void init_quadrics_and_tol(){
    double xmin=1e100,ymin=1e100,zmin=1e100,xmax=-1e100,ymax=-1e100,zmax=-1e100;
    for(int i=0;i<NV;i++){
        const Vec3&p=V[i].p;
        xmin=min(xmin,p.x); ymin=min(ymin,p.y); zmin=min(zmin,p.z);
        xmax=max(xmax,p.x); ymax=max(ymax,p.y); zmax=max(zmax,p.z);
        V[i].q.clear();
    }
    double dx=xmax-xmin, dy=ymax-ymin, dz=zmax-zmin;
    double diag=sqrt(dx*dx+dy*dy+dz*dz);
    hausTol=0.05*diag;
    hausTolSafe=hausTol*(1.0-2e-10);

    vector<Vec3> nsum(NV, Vec3(0,0,0));
    vector<double> asum(NV, 0.0), area(NF, 0.0);
    vector<float> elens;
    elens.reserve((size_t)min<long long>(3LL*NF, 7000000LL));

    for(int i=0;i<NF;i++){
        Face &f=F[i];
        Vec3 p0=V[f.a].p, p1=V[f.b].p, p2=V[f.c].p;
        Vec3 cr=crossv(p1-p0,p2-p0); double l=normv(cr); area[i]=l;
        if((long long)elens.size() < 7000000LL){
            elens.push_back((float)normv(p1-p0));
            elens.push_back((float)normv(p2-p1));
            elens.push_back((float)normv(p0-p2));
        }
        if(l<=0){ f.n=Vec3(0,0,0); continue; }
        Vec3 n=cr/l; f.n=n;
        Vec3 wn=n*l;
        nsum[f.a]+=wn; nsum[f.b]+=wn; nsum[f.c]+=wn;
        asum[f.a]+=l; asum[f.b]+=l; asum[f.c]+=l;
    }

    if(!elens.empty()){
        size_t k=(size_t)(0.90*(elens.size()-1));
        nth_element(elens.begin(), elens.begin()+k, elens.end());
        edgeP90=elens[k];
    }

    vector<double> var(NV,0.0), vars;
    vars.reserve(NV);
    for(int i=0;i<NV;i++) if(asum[i]>0){
        double val=1.0 - normv(nsum[i])/asum[i];
        if(val<0) val=0;
        if(val>1) val=1;
        var[i]=val; vars.push_back(val);
    }
    if(!vars.empty()){
        size_t k=(size_t)(0.95*(vars.size()-1));
        nth_element(vars.begin(), vars.begin()+k, vars.end());
        complexityP95=vars[k];
    }

    const double alpha = 700.0;
    for(int i=0;i<NF;i++){
        Face &f=F[i]; double l=area[i]; if(l<=0) continue;
        Vec3 p0=V[f.a].p; Vec3 n=f.n; double d=-dotv(n,p0);
        double sal=(var[f.a]+var[f.b]+var[f.c])/3.0;
        double mult=1.0 + alpha*min(0.02, sal);
        double w=l*mult;
        V[f.a].q.addPlane(n.x,n.y,n.z,d,w);
        V[f.b].q.addPlane(n.x,n.y,n.z,d,w);
        V[f.c].q.addPlane(n.x,n.y,n.z,d,w);
    }
    vmark.assign(NV,0);
    fmark.assign(NF,0);
}
static bool solve_optimal(const Quadric&q, Vec3& out){
    double a00=q.q[0], a01=q.q[1], a02=q.q[2];
    double a11=q.q[4], a12=q.q[5], a22=q.q[7];
    double b0=-q.q[3], b1=-q.q[6], b2=-q.q[8];
    double det = a00*(a11*a22-a12*a12) - a01*(a01*a22-a12*a02) + a02*(a01*a12-a11*a02);
    double scale = fabs(a00)+fabs(a11)+fabs(a22)+fabs(a01)+fabs(a02)+fabs(a12);
    if(fabs(det) < 1e-12*max(1.0,scale*scale*scale)) return false;
    double c00 = (a11*a22-a12*a12);
    double c01 = -(a01*a22-a02*a12);
    double c02 = (a01*a12-a02*a11);
    double c10 = c01;
    double c11 = (a00*a22-a02*a02);
    double c12 = -(a00*a12-a01*a02);
    double c20 = c02;
    double c21 = c12;
    double c22 = (a00*a11-a01*a01);
    out.x=(c00*b0+c01*b1+c02*b2)/det;
    out.y=(c10*b0+c11*b1+c12*b2)/det;
    out.z=(c20*b0+c21*b1+c22*b2)/det;
    return isfinite(out.x)&&isfinite(out.y)&&isfinite(out.z);
}
static void sphere_union(const Vec3&c1,double r1,const Vec3&c2,double r2,Vec3&co,double&ro){
    Vec3 dvec=c2-c1; double d=normv(dvec);
    if(r1>=d+r2){ co=c1; ro=r1; return; }
    if(r2>=d+r1){ co=c2; ro=r2; return; }
    if(d<1e-300){ co=c1; ro=max(r1,r2); return; }
    ro=(d+r1+r2)*0.5;
    double t=(ro-r1)/d;
    co=c1+dvec*t;
}
static inline bool sphere_covers_point(const Vec3&p,const Vec3&c,double r){
    return normv(p-c) + r <= hausTolSafe;
}

struct Plan{ int keep, rem; Vec3 p; Vec3 sc; double sr; float cost; };
static bool best_plan(int u,int v,Plan&pl){
    if(u==v||!V[u].active||!V[v].active) return false;
    Vec3 sc; double sr; sphere_union(V[u].sc,V[u].sr,V[v].sc,V[v].sr,sc,sr);
    if(sr > hausTolSafe) return false;
    Quadric q=V[u].q+V[v].q;
    Vec3 cand[8]; int nc=0;
    Vec3 opt;
    if(solve_optimal(q,opt)) cand[nc++]=opt;
    cand[nc++]=(V[u].sum+V[v].sum)/(double)(V[u].cnt+V[v].cnt);
    cand[nc++]=(V[u].p+V[v].p)*0.5;
    cand[nc++]=V[u].p;
    cand[nc++]=V[v].p;
    cand[nc++]=sc;
    double best=1e300; Vec3 bp;
    for(int i=0;i<nc;i++){
        Vec3 p=cand[i];
        if(!isfinite(p.x)||!isfinite(p.y)||!isfinite(p.z)) continue;
        if(!sphere_covers_point(p,sc,sr)) continue;
        double e=q.eval(p);
        if(e<0 && e>-1e-18) e=0;
        double len2=norm2(V[u].p-V[v].p);
        e += 1e-12*len2;
        if(e<best){ best=e; bp=p; }
    }
    if(best>=1e299) return false;
    if(V[u].inc.size() >= V[v].inc.size()){ pl.keep=u; pl.rem=v; }
    else { pl.keep=v; pl.rem=u; }
    pl.p=bp; pl.sc=sc; pl.sr=sr; pl.cost=(float)min(best, (double)FLT_MAX/4); return true;
}

struct Cand{ float cost; int u,v; int vu,vv; };
struct CandCmp{ bool operator()(const Cand&a,const Cand&b) const { return a.cost > b.cost; } };
static vector<Cand> heapv;

static inline bool make_cand(int u,int v,Cand &c){
    if(u==v||u<0||v<0||u>=NV||v>=NV) return false;
    if(!V[u].active||!V[v].active) return false;
    Plan pl; if(!best_plan(u,v,pl)) return false;
    c = Cand{pl.cost,u,v,V[u].ver,V[v].ver};
    return true;
}
static void push_cand(int u,int v){
    Cand c; if(!make_cand(u,v,c)) return;
    heapv.push_back(c); push_heap(heapv.begin(),heapv.end(),CandCmp());
}
static void add_cand_raw(int u,int v){
    Cand c; if(!make_cand(u,v,c)) return;
    heapv.push_back(c);
}

static inline bool face_contains(const Face&f,int x){ return f.a==x||f.b==x||f.c==x; }
static void clean_incident(int v){
    if(!V[v].active) return;
    ++fstamp; if(fstamp==INT_MAX){ fill(fmark.begin(),fmark.end(),0); fstamp=1; }
    vector<int> nv; nv.reserve(V[v].inc.size());
    for(int fid: V[v].inc){
        if(fid<0||fid>=NF||!F[fid].active) continue;
        if(!face_contains(F[fid],v)) continue;
        if(fmark[fid]==fstamp) continue;
        fmark[fid]=fstamp; nv.push_back(fid);
    }
    V[v].inc.swap(nv);
}
static void collect_neighbors(int v, vector<int>&nb){
    nb.clear();
    ++vstamp; if(vstamp==INT_MAX){ fill(vmark.begin(),vmark.end(),0); vstamp=1; }
    if(V[v].inc.size()>512) clean_incident(v);
    for(int fid: V[v].inc){
        if(!F[fid].active) continue;
        Face &f=F[fid];
        if(f.a==v){
            if(vmark[f.b]!=vstamp){vmark[f.b]=vstamp; nb.push_back(f.b);}
            if(vmark[f.c]!=vstamp){vmark[f.c]=vstamp; nb.push_back(f.c);}
        }else if(f.b==v){
            if(vmark[f.a]!=vstamp){vmark[f.a]=vstamp; nb.push_back(f.a);}
            if(vmark[f.c]!=vstamp){vmark[f.c]=vstamp; nb.push_back(f.c);}
        }else if(f.c==v){
            if(vmark[f.a]!=vstamp){vmark[f.a]=vstamp; nb.push_back(f.a);}
            if(vmark[f.b]!=vstamp){vmark[f.b]=vstamp; nb.push_back(f.b);}
        }
    }
}
static bool validate_plan(const Plan&pl, vector<int>&unionFaces){
    int u=pl.keep, v=pl.rem;
    if(!V[u].active||!V[v].active||u==v) return false;
    if(V[u].inc.size()>512) clean_incident(u);
    if(V[v].inc.size()>512) clean_incident(v);

    int edgeFaces=0;
    ++fstamp; if(fstamp==INT_MAX){ fill(fmark.begin(),fmark.end(),0); fstamp=1; }
    unionFaces.clear();
    auto addFace=[&](int fid){
        if(fid<0||fid>=NF||!F[fid].active) return;
        Face &f=F[fid];
        bool cu=face_contains(f,u), cv=face_contains(f,v);
        if(!cu && !cv) return;
        if(fmark[fid]!=fstamp){ fmark[fid]=fstamp; unionFaces.push_back(fid); if(cu && cv) edgeFaces++; }
    };
    for(int fid: V[u].inc) addFace(fid);
    for(int fid: V[v].inc) addFace(fid);
    if(edgeFaces!=2) return false;

    vector<int> nu,nv;
    collect_neighbors(u,nu);
    collect_neighbors(v,nv);
    ++vstamp; if(vstamp==INT_MAX){ fill(vmark.begin(),vmark.end(),0); vstamp=1; }
    int stamp=vstamp;
    for(int w:nu) vmark[w]=stamp;
    int common=0;
    for(int w:nv) if(vmark[w]==stamp) common++;
    if(common!=2) return false;

    vector<unsigned long long> pairs; pairs.reserve(unionFaces.size());
    const double minDot = 1e-8;
    for(int fid: unionFaces){
        Face &f=F[fid];
        bool cu=face_contains(f,u), cv=face_contains(f,v);
        if(cu && cv) continue;
        int ids[3]={f.a,f.b,f.c};
        int other[2], k=0;
        Vec3 pp[3];
        for(int i=0;i<3;i++){
            int id=ids[i];
            if(id==u||id==v) pp[i]=pl.p;
            else { pp[i]=V[id].p; if(k<2) other[k++]=id; }
        }
        if(k!=2) return false;
        int x=other[0], y=other[1]; if(x>y) swap(x,y);
        pairs.push_back((unsigned long long)(unsigned int)x<<32 | (unsigned int)y);
        Vec3 cr=crossv(pp[1]-pp[0], pp[2]-pp[0]);
        double l2=norm2(cr);
        if(l2 <= 1e-30) return false;
        double l=sqrt(l2);
        double d=dotv(cr,f.n);
        if(d <= minDot*l) return false;
    }
    sort(pairs.begin(),pairs.end());
    for(size_t i=1;i<pairs.size();i++) if(pairs[i]==pairs[i-1]) return false;
    return true;
}
static void do_collapse(const Plan&pl, vector<int>&unionFaces){
    int u=pl.keep, v=pl.rem;
    for(int fid: unionFaces){
        Face &f=F[fid];
        if(!f.active) continue;
        bool cu=face_contains(f,u), cv=face_contains(f,v);
        if(cu && cv){ f.active=false; continue; }
        if(cv){ if(f.a==v) f.a=u; if(f.b==v) f.b=u; if(f.c==v) f.c=u; }
    }
    V[u].p=pl.p; V[u].sc=pl.sc; V[u].sr=pl.sr; V[u].sum += V[v].sum; V[u].cnt += V[v].cnt; V[u].q += V[v].q;
    V[v].active=false;
    V[u].ver++; V[v].ver++;
    V[u].inc.insert(V[u].inc.end(), V[v].inc.begin(), V[v].inc.end());
    if(V[u].inc.size()>1024) clean_incident(u);
    for(int fid: unionFaces){ if(F[fid].active) recomputeFace(fid); }
    activeVerts--;
}
static void push_edges_around(int u){
    if(!V[u].active) return;
    if(V[u].inc.size()>512) clean_incident(u);
    ++vstamp; if(vstamp==INT_MAX){ fill(vmark.begin(),vmark.end(),0); vstamp=1; }
    int stamp=vstamp;
    vector<int> neigh; neigh.reserve(32);
    for(int fid: V[u].inc){
        if(!F[fid].active) continue;
        Face &f=F[fid];
        int ids[3]={f.a,f.b,f.c};
        if(ids[0]!=u&&ids[1]!=u&&ids[2]!=u) continue;
        for(int i=0;i<3;i++) if(ids[i]!=u && vmark[ids[i]]!=stamp){
            vmark[ids[i]]=stamp; neigh.push_back(ids[i]);
        }
    }
    for(int w:neigh) push_cand(u,w);
}
static int choose_target(){
    if(NV<=100) return NV;
    double r_len = 0.0;
    if(edgeP90 > 0.0) r_len = (edgeP90 / 0.100) * (edgeP90 / 0.100);
    double r_complex;
    if(complexityP95 < 1.0e-4) r_complex = 0.025;
    else if(complexityP95 < 5.0e-4) r_complex = 0.035;
    else if(complexityP95 < 1.5e-3) r_complex = 0.055;
    else if(complexityP95 < 3.5e-3) r_complex = 0.070;
    else r_complex = 0.081;
    double ratio = max(r_len, r_complex);
    if(NV > 100000 && complexityP95 < 5.0e-4) ratio = min(ratio, 0.030);
    if(ratio > 0.55) ratio = 0.55;
    if(ratio < 0.015) ratio = 0.015;
    int target = max(4, (int)ceil(NV*ratio));
    return min(target,NV);
}
static void simplify(){
    init_quadrics_and_tol();
    int target=choose_target();
    if(target>=NV) return;
    heapv.reserve((size_t)min<long long>((long long)NF*3LL + 1000, 8000000LL));
    for(int i=0;i<NF;i++) if(F[i].active){
        add_cand_raw(F[i].a,F[i].b);
        add_cand_raw(F[i].b,F[i].c);
        add_cand_raw(F[i].c,F[i].a);
    }
    make_heap(heapv.begin(),heapv.end(),CandCmp());
    vector<int> unionFaces; unionFaces.reserve(128);
    long long pops=0;
    while(activeVerts>target && !heapv.empty()){
        pop_heap(heapv.begin(),heapv.end(),CandCmp()); Cand c=heapv.back(); heapv.pop_back();
        pops++;
        if((pops & 8191)==0){
            double t=chrono::duration<double>(chrono::steady_clock::now()-startTime).count();
            if(t>19.0) break;
            if(heapv.size()>20000000){
                vector<Cand> nh; nh.reserve((size_t)min<long long>((long long)activeVerts*8LL, 12000000LL));
                heapv.swap(nh);
                for(int fid=0;fid<NF;fid++) if(F[fid].active){
                    push_cand(F[fid].a,F[fid].b); push_cand(F[fid].b,F[fid].c); push_cand(F[fid].c,F[fid].a);
                }
            }
        }
        int u=c.u, v=c.v;
        if(u<0||v<0||u>=NV||v>=NV||!V[u].active||!V[v].active||u==v) continue;
        if(c.vu!=V[u].ver || c.vv!=V[v].ver){ push_cand(u,v); continue; }
        Plan pl; if(!best_plan(u,v,pl)) continue;
        if(pl.cost > c.cost * 1.25f + 1e-18f){ push_cand(u,v); continue; }
        if(!validate_plan(pl,unionFaces)) continue;
        do_collapse(pl,unionFaces);
        push_edges_around(pl.keep);
    }
}
static void save_obj(){
    vector<int> mapv(NV,-1); vector<char> used(NV,0);
    int outF=0;
    for(int i=0;i<NF;i++) if(F[i].active){
        Face &f=F[i];
        if(f.a==f.b||f.b==f.c||f.c==f.a) continue;
        Vec3 cr=crossv(V[f.b].p-V[f.a].p,V[f.c].p-V[f.a].p);
        if(norm2(cr)<=1e-30) continue;
        used[f.a]=used[f.b]=used[f.c]=1; outF++;
    }
    int outV=0; for(int i=0;i<NV;i++) if(used[i]) mapv[i]=outV++;
    string out; out.reserve((size_t)outV*48 + (size_t)outF*28 + 64);
    char line[128];
    out.append(line, snprintf(line,sizeof(line), "%d %d\n", outV, outF));
    for(int i=0;i<NV;i++) if(used[i]){
        Vec3 p=V[i].p;
        out.append(line, snprintf(line,sizeof(line), "v %.10g %.10g %.10g\n", p.x,p.y,p.z));
    }
    for(int i=0;i<NF;i++) if(F[i].active){
        Face &f=F[i];
        if(f.a==f.b||f.b==f.c||f.c==f.a) continue;
        Vec3 cr=crossv(V[f.b].p-V[f.a].p,V[f.c].p-V[f.a].p);
        if(norm2(cr)<=1e-30) continue;
        out.append(line, snprintf(line,sizeof(line), "f %d %d %d\n", mapv[f.a]+1,mapv[f.b]+1,mapv[f.c]+1));
    }
    fwrite(out.data(),1,out.size(),stdout);
}
int main(){
    startTime=chrono::steady_clock::now();
    load_obj();
    simplify();
    save_obj();
    return 0;
}
