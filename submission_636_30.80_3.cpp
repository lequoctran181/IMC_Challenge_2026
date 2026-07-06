#include <bits/stdc++.h>
using namespace std;

struct Vec3 {
    double x,y,z;
    Vec3():x(0),y(0),z(0){}
    Vec3(double X,double Y,double Z):x(X),y(Y),z(Z){}
    Vec3 operator+(const Vec3& o) const { return Vec3(x+o.x,y+o.y,z+o.z); }
    Vec3 operator-(const Vec3& o) const { return Vec3(x-o.x,y-o.y,z-o.z); }
    Vec3 operator*(double s) const { return Vec3(x*s,y*s,z*s); }
    Vec3 operator/(double s) const { return Vec3(x/s,y/s,z/s); }
    Vec3& operator+=(const Vec3& o){ x+=o.x; y+=o.y; z+=o.z; return *this; }
};
static inline double dot3(const Vec3&a,const Vec3&b){return a.x*b.x+a.y*b.y+a.z*b.z;}
static inline Vec3 cross3(const Vec3&a,const Vec3&b){return Vec3(a.y*b.z-a.z*b.y,a.z*b.x-a.x*b.z,a.x*b.y-a.y*b.x);}
static inline double norm2(const Vec3&a){return dot3(a,a);}
static inline double normv(const Vec3&a){return sqrt(norm2(a));}

struct Quadric {
    double q[10];
    Quadric(){ memset(q,0,sizeof(q)); }
    void add(const Quadric& o){ for(int i=0;i<10;i++) q[i]+=o.q[i]; }
    void addPlane(const Vec3& n, double d, double w){
        double a=n.x,b=n.y,c=n.z;
        q[0]+=w*a*a; q[1]+=w*a*b; q[2]+=w*a*c; q[3]+=w*a*d;
        q[4]+=w*b*b; q[5]+=w*b*c; q[6]+=w*b*d;
        q[7]+=w*c*c; q[8]+=w*c*d; q[9]+=w*d*d;
    }
    double eval(const Vec3&p) const {
        double x=p.x,y=p.y,z=p.z;
        return q[0]*x*x + 2*q[1]*x*y + 2*q[2]*x*z + 2*q[3]*x
             + q[4]*y*y + 2*q[5]*y*z + 2*q[6]*y
             + q[7]*z*z + 2*q[8]*z + q[9];
    }
};

struct Face {
    int v[3];
    unsigned char active;
    Vec3 n;
    double area;
};

struct Vertex {
    Vec3 p;
    Quadric q;
    double rad;
    Vec3 mn, mx;
    int ver;
    unsigned char active;
    vector<int> inc;
};

struct Candidate { double cost; Vec3 p; double rad; };

static vector<Vertex> Vv;
static vector<Face> Ff;
static int N0, F0;
static int activeV, activeF;
static double diagLen, hausLimit, radLimit;
static vector<int> markv;
static int markStamp = 1;

static inline bool face_has(const Face& f, int a){ return f.v[0]==a || f.v[1]==a || f.v[2]==a; }
static inline bool face_has2(const Face& f, int a, int b){ return face_has(f,a) && face_has(f,b); }
static inline int face_opp(const Face& f, int a, int b){ return f.v[0]+f.v[1]+f.v[2]-a-b; }

static bool recompute_face(int fid){
    Face &f = Ff[fid];
    Vec3 a=Vv[f.v[0]].p, b=Vv[f.v[1]].p, c=Vv[f.v[2]].p;
    Vec3 cr = cross3(b-a,c-a);
    double l = normv(cr);
    if(l <= 1e-24 || !isfinite(l)) { f.area = 0; f.n=Vec3(0,0,0); return false; }
    f.area = 0.5*l;
    f.n = cr / l;
    return true;
}

static vector<char> slurp_stdin(){
    vector<char> buf; buf.reserve(1<<27);
    char tmp[1<<16]; size_t n;
    while((n=fread(tmp,1,sizeof(tmp),stdin))>0) buf.insert(buf.end(),tmp,tmp+n);
    buf.push_back('\0'); return buf;
}

static inline void skip_ws(char*&p){ while(*p==' '||*p=='\n'||*p=='\r'||*p=='\t') ++p; }

static void load_mesh(){
    vector<char> buf=slurp_stdin(); char* p=buf.data();
    long nv=strtol(p,&p,10), nf=strtol(p,&p,10); N0=(int)nv; F0=(int)nf;
    Vv.resize(N0); Ff.resize(F0);
    double mnx=1e100,mny=1e100,mnz=1e100,mxx=-1e100,mxy=-1e100,mxz=-1e100;
    for(int i=0;i<N0;i++){
        skip_ws(p); if(*p=='v') ++p;
        double x=strtod(p,&p), y=strtod(p,&p), z=strtod(p,&p);
        Vv[i].p=Vec3(x,y,z); Vv[i].rad=0; Vv[i].mn=Vv[i].mx=Vv[i].p; Vv[i].ver=0; Vv[i].active=1;
        mnx=min(mnx,x); mny=min(mny,y); mnz=min(mnz,z); mxx=max(mxx,x); mxy=max(mxy,y); mxz=max(mxz,z);
    }
    for(int i=0;i<F0;i++){
        skip_ws(p); if(*p=='f') ++p;
        int a=(int)strtol(p,&p,10)-1, b=(int)strtol(p,&p,10)-1, c=(int)strtol(p,&p,10)-1;
        Ff[i].v[0]=a; Ff[i].v[1]=b; Ff[i].v[2]=c; Ff[i].active=1;
        Vv[a].inc.push_back(i); Vv[b].inc.push_back(i); Vv[c].inc.push_back(i);
    }
    diagLen=sqrt((mxx-mnx)*(mxx-mnx)+(mxy-mny)*(mxy-mny)+(mxz-mnz)*(mxz-mnz));
    hausLimit=0.05*diagLen;
    radLimit=hausLimit*0.985 - 2e-12;
    if(radLimit < 0) radLimit = hausLimit;
    activeV=N0; activeF=F0;

    for(int i=0;i<F0;i++) recompute_face(i);

    vector<Vec3> vns(N0, Vec3(0,0,0));
    for(int i=0;i<F0;i++){
        Face &f=Ff[i];
        Vec3 an = f.n * max(f.area, 1e-12);
        vns[f.v[0]] += an; vns[f.v[1]] += an; vns[f.v[2]] += an;
    }
    for(int i=0;i<N0;i++){
        double l=normv(vns[i]);
        if(l>0) vns[i]=vns[i]/l;
    }
    for(int i=0;i<F0;i++){
        Face &f=Ff[i];
        Vec3 a=Vv[f.v[0]].p;
        double d=-dot3(f.n,a);
        double curv = 0.0;
        for(int k=0;k<3;k++){
            double dd=fabs(dot3(f.n, vns[f.v[k]]));
            curv += max(0.0, 1.0-dd);
        }
        curv /= 3.0;
        double w=max(f.area, 1e-12) * (1.0 + 2.0*curv + 4.0*curv*curv);
        Quadric q; q.addPlane(f.n,d,w);
        Vv[f.v[0]].q.add(q); Vv[f.v[1]].q.add(q); Vv[f.v[2]].q.add(q);
    }

    markv.assign(N0,0);
}

static void clean_inc(int v){
    if(!Vv[v].active) { Vv[v].inc.clear(); return; }
    vector<int>& L=Vv[v].inc;
    int w=0;
    for(int fid: L){
        if(fid>=0 && fid<F0 && Ff[fid].active && face_has(Ff[fid],v)) L[w++]=fid;
    }
    L.resize(w);
}

static bool solve_optimal(const Quadric&q, Vec3& out){
    double a00=q.q[0], a01=q.q[1], a02=q.q[2];
    double a11=q.q[4], a12=q.q[5], a22=q.q[7];
    double b0=-q.q[3], b1=-q.q[6], b2=-q.q[8];
    double det = a00*(a11*a22-a12*a12) - a01*(a01*a22-a12*a02) + a02*(a01*a12-a11*a02);
    double scale = fabs(a00)+fabs(a11)+fabs(a22)+fabs(a01)+fabs(a02)+fabs(a12)+1e-30;
    if(fabs(det) < 1e-12*scale*scale*scale) return false;
    double inv00=(a11*a22-a12*a12)/det;
    double inv01=(a02*a12-a01*a22)/det;
    double inv02=(a01*a12-a02*a11)/det;
    double inv11=(a00*a22-a02*a02)/det;
    double inv12=(a01*a02-a00*a12)/det;
    double inv22=(a00*a11-a01*a01)/det;
    out.x=inv00*b0+inv01*b1+inv02*b2;
    out.y=inv01*b0+inv11*b1+inv12*b2;
    out.z=inv02*b0+inv12*b1+inv22*b2;
    return isfinite(out.x)&&isfinite(out.y)&&isfinite(out.z);
}

static inline double merged_radius(int a,int b,const Vec3&p){
    double mnx=min(Vv[a].mn.x,Vv[b].mn.x), mny=min(Vv[a].mn.y,Vv[b].mn.y), mnz=min(Vv[a].mn.z,Vv[b].mn.z);
    double mxx=max(Vv[a].mx.x,Vv[b].mx.x), mxy=max(Vv[a].mx.y,Vv[b].mx.y), mxz=max(Vv[a].mx.z,Vv[b].mx.z);
    double best=0.0;
    for(int ix=0; ix<2; ++ix) for(int iy=0; iy<2; ++iy) for(int iz=0; iz<2; ++iz){
        Vec3 c(ix?mxx:mnx, iy?mxy:mny, iz?mxz:mnz);
        best=max(best, normv(c-p));
    }
    return best;
}

static bool basic_candidate(int a,int b,Candidate& best){
    if(a==b || !Vv[a].active || !Vv[b].active) return false;
    Quadric q=Vv[a].q; q.add(Vv[b].q);
    Vec3 pa=Vv[a].p, pb=Vv[b].p;
    Vec3 opts[8]; int cnt=0;
    Vec3 opt;
    if(solve_optimal(q,opt)) opts[cnt++]=opt;
    opts[cnt++]=(pa+pb)*0.5;
    opts[cnt++]=pa; opts[cnt++]=pb;
    opts[cnt++]=pa*0.75+pb*0.25;
    opts[cnt++]=pa*0.25+pb*0.75;

    double bestCost=1e300; Vec3 bestP; double bestR=0; bool ok=false;
    for(int i=0;i<cnt;i++){
        Vec3 p=opts[i];
        double r=merged_radius(a,b,p);
        if(r > radLimit) continue;
        double c=q.eval(p);
        if(!isfinite(c)) continue;
        double edge2=norm2(pa-pb);
        if(i==0 && edge2>0){
            Vec3 mid=(pa+pb)*0.5;
            double off=norm2(p-mid)/(edge2+1e-30);
            if(off>16.0) c += off*max(q.eval(mid),1e-18);
        }
        c += 1e-13 * r*r;
        if(c<bestCost){ bestCost=c; bestP=p; bestR=r; ok=true; }
    }
    if(!ok) return false;
    best.cost=bestCost; best.p=bestP; best.rad=bestR; return true;
}

static void gather_neighbors_mark(int v, int except){
    if(++markStamp == INT_MAX){ fill(markv.begin(), markv.end(), 0); markStamp=1; }
    for(int fid: Vv[v].inc){
        Face &f=Ff[fid]; if(!f.active || !face_has(f,v)) continue;
        for(int k=0;k<3;k++){
            int u=f.v[k];
            if(u!=v && u!=except && Vv[u].active) markv[u]=markStamp;
        }
    }
}

static bool detailed_candidate(int a,int b,Candidate& best, int edgeFaces[2], int opps[2]){
    if(!basic_candidate(a,b,best)) return false;
    clean_inc(a); clean_inc(b);
    int cnt=0;
    vector<int> &La=Vv[a].inc, &Lb=Vv[b].inc;
    const vector<int> *L = (La.size() <= Lb.size()? &La : &Lb);
    for(int fid: *L){
        Face &f=Ff[fid];
        if(f.active && face_has2(f,a,b)){
            if(cnt<2){ edgeFaces[cnt]=fid; opps[cnt]=face_opp(f,a,b); }
            cnt++;
        }
    }
    if(cnt!=2) return false;
    if(opps[0]==opps[1] || opps[0]<0 || opps[1]<0 || !Vv[opps[0]].active || !Vv[opps[1]].active) return false;

    gather_neighbors_mark(a,b);
    int st=markStamp;
    for(int fid: Vv[b].inc){
        Face &f=Ff[fid]; if(!f.active || !face_has(f,b)) continue;
        for(int k=0;k<3;k++){
            int u=f.v[k];
            if(u==b || u==a || !Vv[u].active) continue;
            if(markv[u]==st && u!=opps[0] && u!=opps[1]) return false;
        }
    }

    vector<pair<int,int>> pairs;
    pairs.reserve(Vv[a].inc.size());
    for(int fid: Vv[a].inc){
        if(fid==edgeFaces[0] || fid==edgeFaces[1]) continue;
        Face &f=Ff[fid]; if(!f.active || !face_has(f,a)) continue;
        int x=-1,y=-1;
        for(int k=0;k<3;k++) if(f.v[k]!=a) { if(x<0) x=f.v[k]; else y=f.v[k]; }
        if(x>y) swap(x,y);
        pairs.emplace_back(x,y);
    }
    sort(pairs.begin(), pairs.end());
    for(int fid: Vv[b].inc){
        if(fid==edgeFaces[0] || fid==edgeFaces[1]) continue;
        Face &f=Ff[fid]; if(!f.active || !face_has(f,b)) continue;
        int x=-1,y=-1;
        for(int k=0;k<3;k++) if(f.v[k]!=b) { if(x<0) x=f.v[k]; else y=f.v[k]; }
        if(x>y) swap(x,y);
        if(binary_search(pairs.begin(), pairs.end(), make_pair(x,y))) return false;
    }

    Quadric q=Vv[a].q; q.add(Vv[b].q);
    Vec3 pa=Vv[a].p, pb=Vv[b].p;
    Vec3 opts[8]; int oc=0; Vec3 opt;
    if(solve_optimal(q,opt)) opts[oc++]=opt;
    opts[oc++]=(pa+pb)*0.5; opts[oc++]=pa; opts[oc++]=pb;
    opts[oc++]=pa*0.75+pb*0.25; opts[oc++]=pa*0.25+pb*0.75;

    double bestCost=1e300; Vec3 bestP; double bestR=0; bool ok=false;
    for(int oi=0; oi<oc; ++oi){
        Vec3 p=opts[oi]; double r=merged_radius(a,b,p); if(r>radLimit) continue;
        double cost=q.eval(p); if(!isfinite(cost)) continue;
        double penalty=0.0; bool geom=true;
        for(int pass=0; pass<2 && geom; ++pass){
            int vv = pass? b:a;
            for(int fid: Vv[vv].inc){
                if(fid==edgeFaces[0] || fid==edgeFaces[1]) continue;
                Face &f=Ff[fid]; if(!f.active || !face_has(f,vv)) continue;
                Vec3 p0=Vv[f.v[0]].p, p1=Vv[f.v[1]].p, p2=Vv[f.v[2]].p;
                if(f.v[0]==a || f.v[0]==b) p0=p;
                if(f.v[1]==a || f.v[1]==b) p1=p;
                if(f.v[2]==a || f.v[2]==b) p2=p;
                Vec3 cr=cross3(p1-p0,p2-p0);
                double l2=norm2(cr);
                if(l2 <= 1e-28 || !isfinite(l2)) { geom=false; break; }
                double l=sqrt(l2); Vec3 nn=cr/l;
                double d=dot3(nn,f.n);
                if(d < -0.05) { geom=false; break; }
                double nd = max(0.0, 1.0-d);
                penalty += nd*nd * max(f.area,1e-12);
            }
        }
        if(!geom) continue;
        double edge2=norm2(pa-pb);
        cost += penalty * 0.01;
        cost += 1e-13*r*r + 1e-16*edge2;
        if(cost<bestCost){bestCost=cost; bestP=p; bestR=r; ok=true;}
    }
    if(!ok) return false;
    best.cost=bestCost; best.p=bestP; best.rad=bestR; return true;
}

struct HeapEdge {
    double cost;
    int a,b;
    int va,vb;
    bool operator<(const HeapEdge& o) const { return cost > o.cost; }
};
static priority_queue<HeapEdge> heapq;

static inline void push_edge_basic(int a,int b){
    if(a==b || a<0 || b<0 || !Vv[a].active || !Vv[b].active) return;
    if(a>b) swap(a,b);
    Candidate c; if(!basic_candidate(a,b,c)) return;
    heapq.push({c.cost,a,b,Vv[a].ver,Vv[b].ver});
}

static void push_incident_edges(int v){
    clean_inc(v);
    for(int fid: Vv[v].inc){
        Face &f=Ff[fid]; if(!f.active) continue;
        for(int k=0;k<3;k++){ int u=f.v[k]; if(u!=v) push_edge_basic(v,u); }
    }
}

static bool collapse_edge(int a,int b,const Candidate& cand,const int edgeFaces[2]){
    int keep=a, gone=b;
    if(Vv[b].inc.size() > Vv[a].inc.size()) { keep=b; gone=a; }

    Vec3 newMn(min(Vv[keep].mn.x,Vv[gone].mn.x), min(Vv[keep].mn.y,Vv[gone].mn.y), min(Vv[keep].mn.z,Vv[gone].mn.z));
    Vec3 newMx(max(Vv[keep].mx.x,Vv[gone].mx.x), max(Vv[keep].mx.y,Vv[gone].mx.y), max(Vv[keep].mx.z,Vv[gone].mx.z));

    for(int i=0;i<2;i++) if(Ff[edgeFaces[i]].active){ Ff[edgeFaces[i]].active=0; activeF--; }

    for(int fid: Vv[gone].inc){
        Face &f=Ff[fid];
        if(!f.active || !face_has(f,gone)) continue;
        for(int k=0;k<3;k++) if(f.v[k]==gone) f.v[k]=keep;
        if(f.v[0]==f.v[1] || f.v[1]==f.v[2] || f.v[2]==f.v[0]){
            f.active=0; activeF--; continue;
        }
        if(!recompute_face(fid)) { f.active=0; activeF--; }
    }

    Vv[keep].p=cand.p;
    Vv[keep].rad=cand.rad;
    Vv[keep].mn=newMn; Vv[keep].mx=newMx;
    Vv[keep].q.add(Vv[gone].q);
    Vv[keep].ver++;
    Vv[gone].active=0; Vv[gone].ver++;
    activeV--;

    vector<int> tmp;
    tmp.reserve(Vv[keep].inc.size()+Vv[gone].inc.size());
    tmp.insert(tmp.end(), Vv[keep].inc.begin(), Vv[keep].inc.end());
    tmp.insert(tmp.end(), Vv[gone].inc.begin(), Vv[gone].inc.end());
    Vv[keep].inc.swap(tmp);
    vector<int>().swap(Vv[gone].inc);

    clean_inc(keep);
    push_incident_edges(keep);
    return true;
}

static int choose_target(){
    if(N0 <= 12) return 4;
    double r;
    if(N0 < 8000) r = 0.078;
    else if(N0 < 30000) r = 0.077;
    else if(N0 < 80000) r = 0.076;
    else if(N0 < 500000) r = 0.074;
    else r = 0.073;
    int t=(int)ceil(N0*r);
    if(t<4) t=4;
    return t;
}

static void build_initial_heap(){
    for(int i=0;i<F0;i++) if(Ff[i].active){
        push_edge_basic(Ff[i].v[0],Ff[i].v[1]);
        push_edge_basic(Ff[i].v[1],Ff[i].v[2]);
        push_edge_basic(Ff[i].v[2],Ff[i].v[0]);
    }
}

static void simplify(){
    if(N0 <= 4) return;
    int target=choose_target();
    build_initial_heap();
    long long pops=0, maxPops = (long long)F0*30 + 1000000LL;
    while(activeV > target && !heapq.empty() && pops < maxPops){
        HeapEdge e=heapq.top(); heapq.pop(); ++pops;
        int a=e.a,b=e.b;
        if(a==b || a<0 || b<0 || a>=N0 || b>=N0) continue;
        if(!Vv[a].active || !Vv[b].active) continue;
        if(Vv[a].ver!=e.va || Vv[b].ver!=e.vb){ push_edge_basic(a,b); continue; }
        int edgeFaces[2], opps[2]; Candidate cand;
        if(!detailed_candidate(a,b,cand,edgeFaces,opps)) continue;
        collapse_edge(a,b,cand,edgeFaces);
    }
}

static void save_mesh(){
    vector<int> remap(N0,-1);
    int nv=0,nf=0;
    for(int i=0;i<N0;i++) if(Vv[i].active) remap[i]=nv++;
    for(int i=0;i<F0;i++) if(Ff[i].active){
        int a=Ff[i].v[0],b=Ff[i].v[1],c=Ff[i].v[2];
        if(a!=b && b!=c && c!=a && remap[a]>=0 && remap[b]>=0 && remap[c]>=0) nf++;
    }

    string out;
    out.reserve((size_t)nv*48 + (size_t)nf*28 + 64);
    char line[128];
    out.append(line, snprintf(line,sizeof(line), "%d %d\n", nv, nf));
    for(int i=0;i<N0;i++) if(Vv[i].active){
        Vec3 p=Vv[i].p;
        out.append(line, snprintf(line,sizeof(line), "v %.10g %.10g %.10g\n", p.x,p.y,p.z));
    }
    for(int i=0;i<F0;i++) if(Ff[i].active){
        int a=Ff[i].v[0],b=Ff[i].v[1],c=Ff[i].v[2];
        if(a!=b && b!=c && c!=a && remap[a]>=0 && remap[b]>=0 && remap[c]>=0){
            out.append(line, snprintf(line,sizeof(line), "f %d %d %d\n", remap[a]+1, remap[b]+1, remap[c]+1));
        }
    }
    fwrite(out.data(),1,out.size(),stdout);
}

int main(){
    load_mesh();
    simplify();
    save_mesh();
    return 0;
}
