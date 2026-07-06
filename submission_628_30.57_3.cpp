#pragma GCC optimize("O3,unroll-loops")
#include <bits/stdc++.h>
using namespace std;

struct Vec3{
    double x,y,z;
    Vec3():x(0),y(0),z(0){}
    Vec3(double X,double Y,double Z):x(X),y(Y),z(Z){}
    Vec3 operator+(const Vec3& o)const{return Vec3(x+o.x,y+o.y,z+o.z);} 
    Vec3 operator-(const Vec3& o)const{return Vec3(x-o.x,y-o.y,z-o.z);} 
    Vec3 operator*(double s)const{return Vec3(x*s,y*s,z*s);} 
    Vec3 operator/(double s)const{return Vec3(x/s,y/s,z/s);} 
    Vec3& operator+=(const Vec3& o){x+=o.x;y+=o.y;z+=o.z;return *this;}
};
static inline double dotv(const Vec3&a,const Vec3&b){return a.x*b.x+a.y*b.y+a.z*b.z;}
static inline Vec3 crossv(const Vec3&a,const Vec3&b){return Vec3(a.y*b.z-a.z*b.y,a.z*b.x-a.x*b.z,a.x*b.y-a.y*b.x);} 
static inline double norm2(const Vec3&a){return dotv(a,a);} 
static inline double normv(const Vec3&a){return sqrt(norm2(a));}

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
        return q[0]*x*x + 2*q[1]*x*y + 2*q[2]*x*z + 2*q[3]*x
             + q[4]*y*y + 2*q[5]*y*z + 2*q[6]*y
             + q[7]*z*z + 2*q[8]*z + q[9];
    }
};

struct Face{
    int v[3];
    Vec3 n;
    float area;
    bool alive;
};
struct Vertex{
    Vec3 p;
    Quadric q;
    double radius=0.0;
    float imp=0.0f;
    vector<int> faces;
    vector<int> neigh;
    uint32_t ver=0;
    bool alive=true;
};
struct EdgeOcc{
    uint64_t key;
    int f;
    bool operator<(const EdgeOcc& o)const{return key<o.key;}
};
struct HeapItem{
    double cost;
    int a,b;
    uint32_t va,vb;
};
struct HeapCmp{
    bool operator()(const HeapItem&x,const HeapItem&y)const{return x.cost>y.cost;}
};

static int N,FN;
static vector<Vertex> Vv;
static vector<Face> Ff;
static vector<Vec3> OrigV;
static vector<array<int,3>> OrigF;
static double hausTol=0.05, checkTol=0.05;
static int aliveV=0, aliveF=0;
static chrono::steady_clock::time_point startTime;
static const double INF=1e100;

static vector<char> slurp_stdin(){
    vector<char> buf; buf.reserve(1<<27);
    char chunk[1<<16]; size_t n;
    while((n=fread(chunk,1,sizeof(chunk),stdin))>0) buf.insert(buf.end(),chunk,chunk+n);
    buf.push_back('\0'); return buf;
}
static inline void skipws(char*&p){ while(*p==' '||*p=='\n'||*p=='\r'||*p=='\t') ++p; }
static void load(){
    vector<char> buf=slurp_stdin(); char* p=buf.data();
    N=(int)strtol(p,&p,10); FN=(int)strtol(p,&p,10);
    Vv.resize(N); Ff.resize(FN); OrigV.resize(N); OrigF.resize(FN);
    Vec3 mn(1e9,1e9,1e9), mx(-1e9,-1e9,-1e9);
    for(int i=0;i<N;i++){
        skipws(p); if(*p=='v') ++p;
        double x=strtod(p,&p), y=strtod(p,&p), z=strtod(p,&p);
        Vv[i].p=Vec3(x,y,z); OrigV[i]=Vv[i].p;
        mn.x=min(mn.x,x); mn.y=min(mn.y,y); mn.z=min(mn.z,z);
        mx.x=max(mx.x,x); mx.y=max(mx.y,y); mx.z=max(mx.z,z);
    }
    for(int i=0;i<FN;i++){
        skipws(p); if(*p=='f') ++p;
        int a=(int)strtol(p,&p,10)-1, b=(int)strtol(p,&p,10)-1, c=(int)strtol(p,&p,10)-1;
        Ff[i].v[0]=a; Ff[i].v[1]=b; Ff[i].v[2]=c; Ff[i].alive=true;
        OrigF[i]={a,b,c};
    }
    double diag = sqrt((mx.x-mn.x)*(mx.x-mn.x)+(mx.y-mn.y)*(mx.y-mn.y)+(mx.z-mn.z)*(mx.z-mn.z));
    hausTol = 0.0490 * diag;
    checkTol = 0.04975 * diag;
    if(hausTol<=0) hausTol=1e-9;
    if(checkTol<=0) checkTol=hausTol;
    aliveV=N; aliveF=FN;
}
static inline bool face_has_vertex(const Face& f,int v){return f.v[0]==v||f.v[1]==v||f.v[2]==v;}
static inline int other_in_face(const Face&f,int a,int b){
    for(int i=0;i<3;i++){int x=f.v[i]; if(x!=a&&x!=b) return x;} return -1;
}
static inline bool recompute_face(int fi){
    Face& f=Ff[fi];
    Vec3 a=Vv[f.v[0]].p, b=Vv[f.v[1]].p, c=Vv[f.v[2]].p;
    Vec3 cr=crossv(b-a,c-a); double l=normv(cr);
    if(l<=1e-20) return false;
    f.n=cr/l; f.area=(float)(0.5*l); return true;
}
static void init_geometry(){
    vector<int> cnt(N,0);
    for(int i=0;i<FN;i++){
        recompute_face(i);
        Face& f=Ff[i];
        double w=max(1e-14,(double)f.area);
        Vec3 p=Vv[f.v[0]].p; double d=-dotv(f.n,p);
        for(int j=0;j<3;j++) Vv[f.v[j]].q.addPlane(f.n.x,f.n.y,f.n.z,d,w);
        cnt[f.v[0]]++; cnt[f.v[1]]++; cnt[f.v[2]]++;
    }
    for(int i=0;i<N;i++) Vv[i].faces.reserve(cnt[i]+4);
    for(int i=0;i<FN;i++) for(int j=0;j<3;j++) Vv[Ff[i].v[j]].faces.push_back(i);
}
static inline uint64_t edge_key(int a,int b){ if(a>b) swap(a,b); return (uint64_t)(uint32_t)a<<32 | (uint32_t)b; }
static inline void key_ab(uint64_t key,int&a,int&b){ a=(int)(key>>32); b=(int)(key&0xffffffffu); }
static void build_edges_and_importance(vector<uint64_t>& uniqueEdges){
    vector<EdgeOcc> eo; eo.reserve((size_t)FN*3);
    for(int i=0;i<FN;i++){
        int a=Ff[i].v[0], b=Ff[i].v[1], c=Ff[i].v[2];
        eo.push_back({edge_key(a,b),i}); eo.push_back({edge_key(b,c),i}); eo.push_back({edge_key(c,a),i});
    }
    sort(eo.begin(),eo.end());
    uniqueEdges.reserve(eo.size()/2+10);
    vector<int> ncnt(N,0);
    for(size_t i=0;i<eo.size();){
        size_t j=i+1; while(j<eo.size() && eo[j].key==eo[i].key) ++j;
        int a,b; key_ab(eo[i].key,a,b);
        uniqueEdges.push_back(eo[i].key); ncnt[a]++; ncnt[b]++;
        if(j-i>=2){
            Vec3 n1=Ff[eo[i].f].n, n2=Ff[eo[i+1].f].n;
            double d=max(-1.0,min(1.0,dotv(n1,n2)));
            double s=max(0.0,1.0-d);
            if(n1.x*n2.x<0) s+=0.35;
            if(n1.y*n2.y<0) s+=0.35;
            if(n1.z*n2.z<0) s+=0.35;
            if(s>0){
                float val=(float)min(4.0,s);
                if(val>Vv[a].imp) Vv[a].imp=val;
                if(val>Vv[b].imp) Vv[b].imp=val;
            }
        }
        i=j;
    }
    for(int i=0;i<N;i++) Vv[i].neigh.reserve(ncnt[i]+4);
    for(uint64_t k:uniqueEdges){int a,b; key_ab(k,a,b); Vv[a].neigh.push_back(b); Vv[b].neigh.push_back(a);}
    vector<float> old(N);
    for(int i=0;i<N;i++) old[i]=Vv[i].imp;
    for(int i=0;i<N;i++){
        float m=old[i]*0.85f;
        for(int nb:Vv[i].neigh) m=max(m, old[nb]*0.45f);
        Vv[i].imp=min(3.0f,m);
    }
}

static bool solve3x3(const double A[9], const double b[3], Vec3& x){
    double det = A[0]*(A[4]*A[8]-A[5]*A[7]) - A[1]*(A[3]*A[8]-A[5]*A[6]) + A[2]*(A[3]*A[7]-A[4]*A[6]);
    if(fabs(det)<1e-12) return false;
    double inv=1.0/det;
    double dx = b[0]*(A[4]*A[8]-A[5]*A[7]) - A[1]*(b[1]*A[8]-A[5]*b[2]) + A[2]*(b[1]*A[7]-A[4]*b[2]);
    double dy = A[0]*(b[1]*A[8]-A[5]*b[2]) - b[0]*(A[3]*A[8]-A[5]*A[6]) + A[2]*(A[3]*b[2]-b[1]*A[6]);
    double dz = A[0]*(A[4]*b[2]-b[1]*A[7]) - A[1]*(A[3]*b[2]-b[1]*A[6]) + b[0]*(A[3]*A[7]-A[4]*A[6]);
    x=Vec3(dx*inv,dy*inv,dz*inv); return isfinite(x.x)&&isfinite(x.y)&&isfinite(x.z);
}
struct MergeInfo{ double cost; Vec3 p; double r; bool ok; };
static MergeInfo compute_merge(int a,int b){
    MergeInfo best; best.cost=INF; best.ok=false; best.r=INF;
    if(a==b||!Vv[a].alive||!Vv[b].alive) return best;
    Quadric q=Vv[a].q; q.add(Vv[b].q);
    Vec3 pa=Vv[a].p, pb=Vv[b].p;
    Vec3 cand[7]; int nc=0;
    double A[9]={q.q[0],q.q[1],q.q[2], q.q[1],q.q[4],q.q[5], q.q[2],q.q[5],q.q[7]};
    double rhs[3]={-q.q[3],-q.q[6],-q.q[8]};
    Vec3 opt;
    if(solve3x3(A,rhs,opt)) cand[nc++]=opt;
    cand[nc++]=(pa+pb)*0.5;
    cand[nc++]=pa; cand[nc++]=pb;
    cand[nc++]=pa*0.75+pb*0.25;
    cand[nc++]=pa*0.25+pb*0.75;
    Vec3 dir=pb-pa;
    double l2=norm2(dir);
    if(l2>1e-20){
        double t=dotv(opt-pa,dir)/l2; t=max(0.0,min(1.0,t));
        cand[nc++]=pa+dir*t;
    }
    double invH=1.0/hausTol;
    float im=max(Vv[a].imp,Vv[b].imp);
    double impMul=1.0 + 3.5*(double)im;
    double edgeLen=sqrt(norm2(pa-pb));
    for(int i=0;i<nc;i++){
        Vec3 p=cand[i];
        if(!isfinite(p.x)||!isfinite(p.y)||!isfinite(p.z)) continue;
        double ra=Vv[a].radius + sqrt(norm2(p-pa));
        double rb=Vv[b].radius + sqrt(norm2(p-pb));
        double r=max(ra,rb);
        if(r>hausTol*2.35) continue;
        double rad=r*invH;
        double e=max(0.0,q.eval(p));
        double c=(e + 1e-11*edgeLen*edgeLen + 2e-8*rad*rad*rad*rad) * impMul;
        if(c<best.cost){best.cost=c; best.p=p; best.r=r; best.ok=true;}
    }
    return best;
}
static inline void erase_value(vector<int>&v,int x){
    for(size_t i=0;i<v.size();){ if(v[i]==x){v[i]=v.back(); v.pop_back();} else ++i; }
}
static inline void add_unique(vector<int>&v,int x){
    for(int y:v) if(y==x) return; v.push_back(x);
}
static void clean_faces(int v){
    auto& L=Vv[v].faces;
    size_t w=0;
    for(size_t i=0;i<L.size();i++){
        int f=L[i];
        if(f>=0 && Ff[f].alive && face_has_vertex(Ff[f],v)) L[w++]=f;
    }
    L.resize(w);
    if(L.size()>1){ sort(L.begin(),L.end()); L.erase(unique(L.begin(),L.end()),L.end()); }
}
static void clean_neigh_basic(int v){
    auto& L=Vv[v].neigh;
    size_t w=0;
    for(int x:L) if(x!=v && x>=0 && x<N && Vv[x].alive) L[w++]=x;
    L.resize(w);
    if(L.size()>1){ sort(L.begin(),L.end()); L.erase(unique(L.begin(),L.end()),L.end()); }
}
static void rebuild_neigh(int v){
    clean_faces(v);
    vector<int> nb; nb.reserve(Vv[v].faces.size()+4);
    for(int fi:Vv[v].faces){
        Face& f=Ff[fi];
        for(int k=0;k<3;k++){int x=f.v[k]; if(x!=v && Vv[x].alive) nb.push_back(x);} 
    }
    sort(nb.begin(),nb.end()); nb.erase(unique(nb.begin(),nb.end()),nb.end());
    Vv[v].neigh.swap(nb);
}
static bool neighbor_exists(int a,int b){
    if(Vv[a].neigh.size()>Vv[b].neigh.size()) swap(a,b);
    for(int x:Vv[a].neigh) if(x==b) return true;
    return false;
}
static bool validate_collapse(int keep,int drop,const Vec3& pnew,double rnew, vector<int>& edgeFaces, vector<int>& affected){
    if(keep==drop || !Vv[keep].alive || !Vv[drop].alive) return false;
    if(rnew>hausTol*2.35) return false;
    rebuild_neigh(keep); rebuild_neigh(drop);
    if(!neighbor_exists(keep,drop)) return false;
    clean_faces(keep); clean_faces(drop);
    edgeFaces.clear(); affected.clear();
    int opp[2]={-1,-1};
    const vector<int>& scan = (Vv[keep].faces.size()<Vv[drop].faces.size()?Vv[keep].faces:Vv[drop].faces);
    for(int fi:scan){
        Face& f=Ff[fi];
        if(f.alive && face_has_vertex(f,keep) && face_has_vertex(f,drop)){
            if((int)edgeFaces.size()<2) opp[edgeFaces.size()]=other_in_face(f,keep,drop);
            edgeFaces.push_back(fi);
        }
    }
    if(edgeFaces.size()!=2 || opp[0]<0 || opp[1]<0 || opp[0]==opp[1]) return false;
    static vector<int> mark;
    static int tag=1;
    if((int)mark.size()<N) mark.assign(N,0);
    if(++tag==INT_MAX){ fill(mark.begin(),mark.end(),0); tag=1; }
    for(int x:Vv[keep].neigh) if(Vv[x].alive) mark[x]=tag;
    int common=0; bool has0=false, has1=false;
    for(int x:Vv[drop].neigh) if(x>=0 && x<N && Vv[x].alive && mark[x]==tag){
        common++; if(x==opp[0]) has0=true; if(x==opp[1]) has1=true;
    }
    if(common!=2 || !has0 || !has1) return false;
    for(int fi:Vv[keep].faces) if(Ff[fi].alive && !(fi==edgeFaces[0]||fi==edgeFaces[1])) affected.push_back(fi);
    for(int fi:Vv[drop].faces) if(Ff[fi].alive && !(fi==edgeFaces[0]||fi==edgeFaces[1])) affected.push_back(fi);
    sort(affected.begin(),affected.end()); affected.erase(unique(affected.begin(),affected.end()),affected.end());
    const double epsArea=1e-20;
    for(int fi:affected){
        Face& f=Ff[fi];
        Vec3 p[3];
        for(int k=0;k<3;k++){
            int id=f.v[k];
            if(id==keep || id==drop) p[k]=pnew; else p[k]=Vv[id].p;
        }
        Vec3 cr=crossv(p[1]-p[0],p[2]-p[0]);
        double l=normv(cr);
        if(l<=epsArea) return false;
        Vec3 nn=cr/l;
        double d=dotv(nn,f.n);
        if(d < -0.02) return false;
    }
    return true;
}
static void push_edge(vector<HeapItem>&heap,int a,int b){
    if(a==b||!Vv[a].alive||!Vv[b].alive) return;
    if(a>b) swap(a,b);
    MergeInfo mi=compute_merge(a,b);
    if(!mi.ok) return;
    heap.push_back({mi.cost,a,b,Vv[a].ver,Vv[b].ver});
    push_heap(heap.begin(),heap.end(),HeapCmp());
}
static bool perform_collapse(int a,int b,const MergeInfo&mi, vector<HeapItem>&heap){
    int keep=a, drop=b;
    if(Vv[keep].faces.size()<Vv[drop].faces.size()) swap(keep,drop);
    if(Vv[drop].imp>Vv[keep].imp*1.35f) swap(keep,drop);
    vector<int> edgeFaces, affected;
    if(!validate_collapse(keep,drop,mi.p,mi.r,edgeFaces,affected)){
        swap(keep,drop);
        if(!validate_collapse(keep,drop,mi.p,mi.r,edgeFaces,affected)) return false;
    }
    vector<int> oldNeighbors=Vv[keep].neigh;
    oldNeighbors.insert(oldNeighbors.end(),Vv[drop].neigh.begin(),Vv[drop].neigh.end());
    sort(oldNeighbors.begin(),oldNeighbors.end()); oldNeighbors.erase(unique(oldNeighbors.begin(),oldNeighbors.end()),oldNeighbors.end());
    for(int fi:edgeFaces){ if(Ff[fi].alive){ Ff[fi].alive=false; aliveF--; } }
    for(int fi:Vv[drop].faces){
        Face& f=Ff[fi];
        if(!f.alive) continue;
        for(int k=0;k<3;k++) if(f.v[k]==drop) f.v[k]=keep;
        Vv[keep].faces.push_back(fi);
    }
    Vv[keep].p=mi.p;
    Vv[keep].radius=mi.r;
    Vv[keep].q.add(Vv[drop].q);
    Vv[keep].imp=max(Vv[keep].imp,Vv[drop].imp*0.92f);
    Vv[drop].alive=false;
    Vv[drop].faces.clear();
    Vv[drop].neigh.clear();
    Vv[keep].ver++;
    Vv[drop].ver++;
    aliveV--;
    clean_faces(keep);
    for(int fi:Vv[keep].faces) if(Ff[fi].alive) recompute_face(fi);
    rebuild_neigh(keep);
    for(int nb:oldNeighbors){
        if(nb==keep||nb==drop||nb<0||nb>=N||!Vv[nb].alive) continue;
        erase_value(Vv[nb].neigh,drop);
        bool shares=false;
        for(int fi:Vv[nb].faces){
            if(Ff[fi].alive && face_has_vertex(Ff[fi],keep) && face_has_vertex(Ff[fi],nb)){ shares=true; break; }
        }
        if(shares) add_unique(Vv[nb].neigh,keep); else erase_value(Vv[nb].neigh,keep);
        clean_neigh_basic(nb);
    }
    for(int nb:Vv[keep].neigh) push_edge(heap,keep,nb);
    return true;
}
static int target_count(){
    if(N<=20) return max(4,N-1);
    double r;
    if(N<=7000) r=0.090;
    else if(N<=30000) r=0.078;
    else if(N<=80000) r=0.072;
    else if(N<=500000) r=0.066;
    else r=0.055;
    double avgImp=0;
    for(int i=0;i<N;i++) avgImp += min(2.5f,Vv[i].imp);
    avgImp/=max(1,N);
    double adj = min(0.030, max(0.0, (avgImp-0.18)*0.010));
    r += adj;
    int t=(int)ceil(N*r);
    t=max(t,4);
    t=min(t,N);
    return t;
}
static void simplify(){
    startTime=chrono::steady_clock::now();
    if(N<4) return;
    init_geometry();
    vector<uint64_t> uniqueEdges;
    build_edges_and_importance(uniqueEdges);
    int target=target_count();
    vector<HeapItem> heap; heap.reserve(uniqueEdges.size()*2+1024);
    for(uint64_t k:uniqueEdges){int a,b; key_ab(k,a,b); MergeInfo mi=compute_merge(a,b); if(mi.ok) heap.push_back({mi.cost,min(a,b),max(a,b),Vv[min(a,b)].ver,Vv[max(a,b)].ver});}
    make_heap(heap.begin(),heap.end(),HeapCmp());
    uniqueEdges.clear();
    long long pops=0;
    const double TIME_LIMIT = (N>300000?18.8:19.6);
    while(aliveV>target && !heap.empty()){
        pop_heap(heap.begin(),heap.end(),HeapCmp());
        HeapItem it=heap.back(); heap.pop_back(); pops++;
        int a=it.a,b=it.b;
        if(a<0||b<0||a>=N||b>=N||!Vv[a].alive||!Vv[b].alive) continue;
        if(Vv[a].ver!=it.va || Vv[b].ver!=it.vb) continue;
        if(!neighbor_exists(a,b)) continue;
        MergeInfo mi=compute_merge(a,b);
        if(!mi.ok) continue;
        if(mi.cost>it.cost*1.25 + 1e-18){ push_edge(heap,a,b); continue; }
        perform_collapse(a,b,mi,heap);
        if((pops&4095)==0){
            double elapsed=chrono::duration<double>(chrono::steady_clock::now()-startTime).count();
            if(elapsed>TIME_LIMIT) break;
            if(heap.size()>12000000ULL && aliveV < target*13/10) break;
        }
    }
}
static bool validate_final(vector<int>&mapOld, vector<Vec3>&outV, vector<array<int,3>>&outF){
    mapOld.assign(N,-1); outV.clear(); outF.clear();
    vector<array<int,3>> tmp; tmp.reserve(aliveF);
    vector<unsigned char> used(N,0);
    for(int i=0;i<FN;i++) if(Ff[i].alive){
        int a=Ff[i].v[0],b=Ff[i].v[1],c=Ff[i].v[2];
        if(a<0||b<0||c<0||a>=N||b>=N||c>=N||!Vv[a].alive||!Vv[b].alive||!Vv[c].alive||a==b||b==c||c==a) return false;
        Vec3 cr=crossv(Vv[b].p-Vv[a].p,Vv[c].p-Vv[a].p);
        if(norm2(cr)<=1e-24) return false;
        tmp.push_back({a,b,c}); used[a]=used[b]=used[c]=1;
    }
    if(tmp.empty()) return false;
    outV.reserve(aliveV);
    for(int i=0;i<N;i++) if(used[i]){ mapOld[i]=(int)outV.size(); outV.push_back(Vv[i].p); }
    outF.reserve(tmp.size());
    for(auto &t:tmp) outF.push_back({mapOld[t[0]],mapOld[t[1]],mapOld[t[2]]});
    vector<uint64_t> edges; edges.reserve(outF.size()*3);
    for(auto &f:outF){
        edges.push_back(edge_key(f[0],f[1])); edges.push_back(edge_key(f[1],f[2])); edges.push_back(edge_key(f[2],f[0]));
    }
    sort(edges.begin(),edges.end());
    for(size_t i=0;i<edges.size();){ size_t j=i+1; while(j<edges.size()&&edges[j]==edges[i]) ++j; if(j-i!=2) return false; i=j; }
    return true;
}
struct GridIndex{
    double cell;
    const vector<Vec3>* pts;
    unordered_map<unsigned long long, vector<int>> buckets;
    static unsigned long long pack(int x,int y,int z){
        const int B=(1<<20);
        unsigned long long ux=(unsigned long long)((x+B)&0x1fffff);
        unsigned long long uy=(unsigned long long)((y+B)&0x1fffff);
        unsigned long long uz=(unsigned long long)((z+B)&0x1fffff);
        return (ux<<42) ^ (uy<<21) ^ uz;
    }
    inline int ci(double v) const { return (int)floor(v/cell); }
    void build(const vector<Vec3>&p,double c){
        pts=&p; cell=max(c,1e-12); buckets.clear(); buckets.reserve(p.size()*2+10);
        for(int i=0;i<(int)p.size();++i) addPoint(i,p[i]);
    }
    void addPoint(int id,const Vec3&p){
        buckets[pack(ci(p.x),ci(p.y),ci(p.z))].push_back(id);
    }
    bool nearestWithin(const Vec3&p,double maxDist2,int&best,double&bestD2) const{
        best=-1; bestD2=maxDist2;
        int ix=ci(p.x), iy=ci(p.y), iz=ci(p.z);
        for(int dx=-1;dx<=1;dx++) for(int dy=-1;dy<=1;dy++) for(int dz=-1;dz<=1;dz++){
            auto it=buckets.find(pack(ix+dx,iy+dy,iz+dz));
            if(it==buckets.end()) continue;
            for(int id:it->second){
                const Vec3&q=(*pts)[id]; double d2=norm2(p-q);
                if(d2<=bestD2){bestD2=d2; best=id;}
            }
        }
        return best>=0;
    }
    int nearestExpanded(const Vec3&p,int maxRing=8) const{
        int best=-1; double bestD2=INF;
        int ix=ci(p.x), iy=ci(p.y), iz=ci(p.z);
        for(int r=0;r<=maxRing && best<0;r++){
            for(int dx=-r;dx<=r;dx++) for(int dy=-r;dy<=r;dy++) for(int dz=-r;dz<=r;dz++){
                if(max({abs(dx),abs(dy),abs(dz)})!=r) continue;
                auto it=buckets.find(pack(ix+dx,iy+dy,iz+dz));
                if(it==buckets.end()) continue;
                for(int id:it->second){
                    const Vec3&q=(*pts)[id]; double d2=norm2(p-q);
                    if(d2<bestD2){bestD2=d2; best=id;}
                }
            }
        }
        if(best>=0) return best;
        const vector<Vec3>& P=*pts;
        for(int i=0;i<(int)P.size();++i){ double d2=norm2(p-P[i]); if(d2<bestD2){bestD2=d2; best=i;} }
        return best;
    }
};
struct RFace{int a,b,c; bool alive;};
static inline bool tri_ok(const Vec3&a,const Vec3&b,const Vec3&c){ return norm2(crossv(b-a,c-a))>1e-24; }
static bool final_edge_valid(const vector<Vec3>&verts,const vector<RFace>&faces){
    vector<unsigned long long> ed; ed.reserve(faces.size()*3);
    int liveF=0;
    vector<unsigned char> used(verts.size(),0);
    for(const auto&f:faces) if(f.alive){
        if(f.a<0||f.b<0||f.c<0||f.a>= (int)verts.size()||f.b>= (int)verts.size()||f.c>= (int)verts.size()) return false;
        if(f.a==f.b||f.b==f.c||f.c==f.a) return false;
        if(!tri_ok(verts[f.a],verts[f.b],verts[f.c])) return false;
        used[f.a]=used[f.b]=used[f.c]=1;
        ed.push_back(edge_key(f.a,f.b)); ed.push_back(edge_key(f.b,f.c)); ed.push_back(edge_key(f.c,f.a)); liveF++;
    }
    if(verts.empty()||liveF==0) return false;
    for(unsigned char u:used) if(!u) return false;
    sort(ed.begin(),ed.end());
    for(size_t i=0;i<ed.size();){ size_t j=i+1; while(j<ed.size()&&ed[j]==ed[i]) ++j; if(j-i!=2) return false; i=j; }
    return true;
}
static bool split_face_for_point(vector<Vec3>&verts, vector<RFace>&faces, vector<vector<int>>&inc, int nearV, const Vec3&p, int&newId){
    auto try_face = [&](int fi)->bool{
        if(fi<0||fi>=(int)faces.size()||!faces[fi].alive) return false;
        RFace f=faces[fi];
        if(!tri_ok(verts[f.a],verts[f.b],p) || !tri_ok(verts[f.b],verts[f.c],p) || !tri_ok(verts[f.c],verts[f.a],p)) return false;
        newId=(int)verts.size(); verts.push_back(p); inc.emplace_back();
        faces[fi].alive=false;
        int f1=(int)faces.size(); faces.push_back({f.a,f.b,newId,true});
        int f2=(int)faces.size(); faces.push_back({f.b,f.c,newId,true});
        int f3=(int)faces.size(); faces.push_back({f.c,f.a,newId,true});
        inc[f.a].push_back(f1); inc[f.a].push_back(f3);
        inc[f.b].push_back(f1); inc[f.b].push_back(f2);
        inc[f.c].push_back(f2); inc[f.c].push_back(f3);
        inc[newId].push_back(f1); inc[newId].push_back(f2); inc[newId].push_back(f3);
        return true;
    };
    if(nearV>=0 && nearV<(int)inc.size()){
        auto &L=inc[nearV];
        for(int fi:L) if(try_face(fi)) return true;
        size_t w=0; for(int fi:L) if(fi>=0&&fi<(int)faces.size()&&faces[fi].alive&&(faces[fi].a==nearV||faces[fi].b==nearV||faces[fi].c==nearV)) L[w++]=fi; L.resize(w);
        for(int fi:L) if(try_face(fi)) return true;
    }
    int tries=0;
    for(int fi=0; fi<(int)faces.size() && tries<500; ++fi) if(faces[fi].alive){ tries++; if(try_face(fi)) return true; }
    return false;
}
static bool repair_hausdorff(vector<Vec3>&outV, vector<array<int,3>>&outF){
    if(outV.empty()||outF.empty()) return false;
    GridIndex gridOrig; gridOrig.build(OrigV, checkTol);
    double tol2=checkTol*checkTol;
    for(Vec3 &p:outV){
        int bi; double bd;
        if(!gridOrig.nearestWithin(p,tol2,bi,bd)){
            int ni=gridOrig.nearestExpanded(p,12);
            if(ni<0) return false;
            p=OrigV[ni];
        }
    }
    vector<RFace> faces; faces.reserve(outF.size()+1024);
    vector<vector<int>> inc(outV.size());
    for(auto &f:outF){
        int id=(int)faces.size(); faces.push_back({f[0],f[1],f[2],true});
        inc[f[0]].push_back(id); inc[f[1]].push_back(id); inc[f[2]].push_back(id);
    }
    GridIndex gridOut; gridOut.build(outV, checkTol);
    int inserted=0;
    const int maxInsert=max(0, N - (int)outV.size());
    for(int i=0;i<N;i++){
        int bi; double bd;
        if(gridOut.nearestWithin(OrigV[i],tol2,bi,bd)) continue;
        if(inserted>=maxInsert) return false;
        int nearV=gridOut.nearestExpanded(OrigV[i],10);
        int newId=-1;
        if(!split_face_for_point(outV,faces,inc,nearV,OrigV[i],newId)) return false;
        gridOut.pts=&outV; gridOut.addPoint(newId,outV[newId]);
        inserted++;
        if((int)outV.size()>N) return false;
    }
    if(!final_edge_valid(outV,faces)) return false;
    outF.clear(); outF.reserve(faces.size());
    for(auto &f:faces) if(f.alive) outF.push_back({f.a,f.b,f.c});
    return true;
}
static void save_mesh(const vector<Vec3>&outV,const vector<array<int,3>>&outF){
    string out;
    out.reserve(outV.size()*42 + outF.size()*26 + 64);
    char line[128];
    out.append(line, snprintf(line,sizeof(line),"%zu %zu\n",outV.size(),outF.size()));
    for(const Vec3&p:outV) out.append(line, snprintf(line,sizeof(line),"v %.10g %.10g %.10g\n",p.x,p.y,p.z));
    for(auto &f:outF) out.append(line, snprintf(line,sizeof(line),"f %d %d %d\n",f[0]+1,f[1]+1,f[2]+1));
    fwrite(out.data(),1,out.size(),stdout);
}
static void save_original(){
    string out; out.reserve(OrigV.size()*42+OrigF.size()*26+64); char line[128];
    out.append(line, snprintf(line,sizeof(line),"%d %d\n",N,FN));
    for(const Vec3&p:OrigV) out.append(line, snprintf(line,sizeof(line),"v %.10g %.10g %.10g\n",p.x,p.y,p.z));
    for(auto &f:OrigF) out.append(line, snprintf(line,sizeof(line),"f %d %d %d\n",f[0]+1,f[1]+1,f[2]+1));
    fwrite(out.data(),1,out.size(),stdout);
}
int main(){
    load();
    simplify();
    vector<int> mapOld; vector<Vec3> outV; vector<array<int,3>> outF;
    if(validate_final(mapOld,outV,outF) && repair_hausdorff(outV,outF)) save_mesh(outV,outF);
    else save_original();
    return 0;
}
