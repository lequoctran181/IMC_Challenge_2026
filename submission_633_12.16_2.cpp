#include <bits/stdc++.h>
using namespace std;

struct Vec3{
    double x,y,z;
    Vec3():x(0),y(0),z(0){}
    Vec3(double X,double Y,double Z):x(X),y(Y),z(Z){}
    Vec3 operator+(const Vec3& o) const { return Vec3(x+o.x,y+o.y,z+o.z); }
    Vec3 operator-(const Vec3& o) const { return Vec3(x-o.x,y-o.y,z-o.z); }
    Vec3 operator*(double s) const { return Vec3(x*s,y*s,z*s); }
    Vec3 operator/(double s) const { return Vec3(x/s,y/s,z/s); }
    Vec3& operator+=(const Vec3& o){ x+=o.x; y+=o.y; z+=o.z; return *this; }
};
static inline double dotv(const Vec3&a,const Vec3&b){return a.x*b.x+a.y*b.y+a.z*b.z;}
static inline Vec3 crossv(const Vec3&a,const Vec3&b){return Vec3(a.y*b.z-a.z*b.y,a.z*b.x-a.x*b.z,a.x*b.y-a.y*b.x);} 
static inline double norm2(const Vec3&a){return dotv(a,a);} 
static inline double normv(const Vec3&a){return sqrt(norm2(a));}

struct Quadric{
    double q[10];
    Quadric(){ memset(q,0,sizeof(q)); }
    inline void clear(){ memset(q,0,sizeof(q)); }
    inline void add(const Quadric& o){ for(int i=0;i<10;i++) q[i]+=o.q[i]; }
    inline void addPlane(double a,double b,double c,double d,double w){
        q[0]+=w*a*a; q[1]+=w*a*b; q[2]+=w*a*c; q[3]+=w*a*d;
        q[4]+=w*b*b; q[5]+=w*b*c; q[6]+=w*b*d;
        q[7]+=w*c*c; q[8]+=w*c*d; q[9]+=w*d*d;
    }
    inline double eval(const Vec3&p) const{
        double x=p.x,y=p.y,z=p.z;
        return q[0]*x*x + 2*q[1]*x*y + 2*q[2]*x*z + 2*q[3]*x
             + q[4]*y*y + 2*q[5]*y*z + 2*q[6]*y
             + q[7]*z*z + 2*q[8]*z + q[9];
    }
};

struct Face{
    int v[3];
    float nx,ny,nz;
    unsigned char active;
};

static int NV,NF;
static vector<Vec3> P;
static vector<Face> faces;
static vector<Quadric> Q;
static vector<double> radBound;
static vector<int> clusterSize;
static vector<unsigned char> vActive;
static vector<int> versionV;
static vector<int> deg0, startInc, poolInc;
static vector<int> extraHead, extraFace, extraNext;
static vector<unsigned long long> initialEdges;
static vector<int> markV, markF;
static int stampV=1, stampF=1;
static double tolH=0.0, tolSafe=0.0, diagAABB=0.0;
static int activeVertices=0;
static const double INF_COST = 1e100;

static vector<char> slurp_stdin(){
    vector<char> buf; buf.reserve(1<<27); char tmp[1<<16]; size_t n;
    while((n=fread(tmp,1,sizeof(tmp),stdin))>0) buf.insert(buf.end(), tmp, tmp+n);
    buf.push_back('\0'); return buf;
}
static inline void skipws(char*&p){ while(*p==' '||*p=='\n'||*p=='\r'||*p=='\t') ++p; }

static inline bool faceHas(const Face& f,int v){ return f.v[0]==v || f.v[1]==v || f.v[2]==v; }
static inline int faceCorner(const Face& f,int v){ if(f.v[0]==v) return 0; if(f.v[1]==v) return 1; if(f.v[2]==v) return 2; return -1; }

static bool recomputeFace(int fid){
    Face &f=faces[fid];
    Vec3 a=P[f.v[0]], b=P[f.v[1]], c=P[f.v[2]];
    Vec3 cr=crossv(b-a,c-a);
    double n=normv(cr);
    if(!(n>1e-18)) return false;
    f.nx=(float)(cr.x/n); f.ny=(float)(cr.y/n); f.nz=(float)(cr.z/n);
    return true;
}

static void load_input(){
    vector<char> buf=slurp_stdin(); char* p=buf.data();
    NV=(int)strtol(p,&p,10); NF=(int)strtol(p,&p,10);
    P.resize(NV); faces.resize(NF); deg0.assign(NV,0);
    Vec3 mn(1e100,1e100,1e100), mx(-1e100,-1e100,-1e100);
    for(int i=0;i<NV;i++){
        skipws(p); if(*p=='v') ++p;
        double x=strtod(p,&p), y=strtod(p,&p), z=strtod(p,&p);
        P[i]=Vec3(x,y,z);
        mn.x=min(mn.x,x); mn.y=min(mn.y,y); mn.z=min(mn.z,z);
        mx.x=max(mx.x,x); mx.y=max(mx.y,y); mx.z=max(mx.z,z);
    }
    initialEdges.clear(); initialEdges.reserve((size_t)3*NF);
    for(int i=0;i<NF;i++){
        skipws(p); if(*p=='f') ++p;
        int a=(int)strtol(p,&p,10)-1;
        int b=(int)strtol(p,&p,10)-1;
        int c=(int)strtol(p,&p,10)-1;
        faces[i].v[0]=a; faces[i].v[1]=b; faces[i].v[2]=c; faces[i].active=1;
        deg0[a]++; deg0[b]++; deg0[c]++;
        int e0a=min(a,b), e0b=max(a,b);
        int e1a=min(b,c), e1b=max(b,c);
        int e2a=min(c,a), e2b=max(c,a);
        initialEdges.push_back((unsigned long long)(unsigned int)e0a<<32 | (unsigned int)e0b);
        initialEdges.push_back((unsigned long long)(unsigned int)e1a<<32 | (unsigned int)e1b);
        initialEdges.push_back((unsigned long long)(unsigned int)e2a<<32 | (unsigned int)e2b);
    }
    diagAABB=normv(mx-mn);
    tolH=0.05*diagAABB;
    tolSafe=tolH*0.9995 - 1e-12;
    if(tolSafe<=0) tolSafe=tolH*0.999;
    activeVertices=NV;
    startInc.assign(NV+1,0);
    for(int i=0;i<NV;i++) startInc[i+1]=startInc[i]+deg0[i];
    poolInc.assign(startInc[NV],0);
    vector<int> cur=startInc;
    for(int i=0;i<NF;i++){
        for(int k=0;k<3;k++){ int v=faces[i].v[k]; poolInc[cur[v]++]=i; }
    }
    Q.assign(NV, Quadric());
    for(int i=0;i<NF;i++){
        Vec3 a=P[faces[i].v[0]], b=P[faces[i].v[1]], c=P[faces[i].v[2]];
        Vec3 cr=crossv(b-a,c-a); double len=normv(cr);
        if(len<=0) { faces[i].nx=faces[i].ny=faces[i].nz=0; continue; }
        double nx=cr.x/len, ny=cr.y/len, nz=cr.z/len;
        faces[i].nx=(float)nx; faces[i].ny=(float)ny; faces[i].nz=(float)nz;
        double d=-(nx*a.x+ny*a.y+nz*a.z);
        double area=0.5*len;
        double w=max(area, 1e-16);
        Quadric fq; fq.addPlane(nx,ny,nz,d,w);
        Q[faces[i].v[0]].add(fq); Q[faces[i].v[1]].add(fq); Q[faces[i].v[2]].add(fq);
    }
    radBound.assign(NV,0.0); clusterSize.assign(NV,1); vActive.assign(NV,1); versionV.assign(NV,0);
    extraHead.assign(NV,-1); extraFace.reserve((size_t)min(8000000, max(1000, 6*NV))); extraNext.reserve(extraFace.capacity());
    markV.assign(NV,0); markF.assign(NF,0);
    sort(initialEdges.begin(), initialEdges.end());
    initialEdges.erase(unique(initialEdges.begin(), initialEdges.end()), initialEdges.end());
}

static inline void appendIncident(int v,int fid){
    int idx=(int)extraFace.size(); extraFace.push_back(fid); extraNext.push_back(extraHead[v]); extraHead[v]=idx;
}

static inline void bumpStampV(){ if(++stampV==INT_MAX){ fill(markV.begin(),markV.end(),0); stampV=1; } }
static inline void bumpStampF(){ if(++stampF==INT_MAX){ fill(markF.begin(),markF.end(),0); stampF=1; } }

static void collectData(int v, vector<int>& flist, vector<int>& nlist){
    flist.clear(); nlist.clear();
    bumpStampF(); bumpStampV();
    auto consider = [&](int fid){
        if(fid<0 || fid>=NF) return;
        if(markF[fid]==stampF) return;
        Face &f=faces[fid];
        if(!f.active) return;
        if(!faceHas(f,v)) return;
        markF[fid]=stampF; flist.push_back(fid);
        for(int k=0;k<3;k++){
            int u=f.v[k];
            if(u!=v && vActive[u] && markV[u]!=stampV){ markV[u]=stampV; nlist.push_back(u); }
        }
    };
    for(int i=startInc[v]; i<startInc[v+1]; ++i) consider(poolInc[i]);
    for(int e=extraHead[v]; e!=-1; e=extraNext[e]) consider(extraFace[e]);
}

static bool solve3x3(const Quadric& q, Vec3& sol){
    double a00=q.q[0], a01=q.q[1], a02=q.q[2];
    double a10=q.q[1], a11=q.q[4], a12=q.q[5];
    double a20=q.q[2], a21=q.q[5], a22=q.q[7];
    double b0=-q.q[3], b1=-q.q[6], b2=-q.q[8];
    double det = a00*(a11*a22-a12*a21) - a01*(a10*a22-a12*a20) + a02*(a10*a21-a11*a20);
    if(fabs(det) < 1e-14) return false;
    double id=1.0/det;
    double dx = b0*(a11*a22-a12*a21) - a01*(b1*a22-a12*b2) + a02*(b1*a21-a11*b2);
    double dy = a00*(b1*a22-a12*b2) - b0*(a10*a22-a12*a20) + a02*(a10*b2-b1*a20);
    double dz = a00*(a11*b2-b1*a21) - a01*(a10*b2-b1*a20) + b0*(a10*a21-a11*a20);
    sol=Vec3(dx*id,dy*id,dz*id);
    return isfinite(sol.x)&&isfinite(sol.y)&&isfinite(sol.z);
}

struct PosChoice{ Vec3 p; double cost; double rad; };

static void addChoice(vector<PosChoice>& choices, const Quadric& q, int a, int b, const Vec3& p){
    if(!isfinite(p.x)||!isfinite(p.y)||!isfinite(p.z)) return;
    double ra = radBound[a] + normv(P[a]-p);
    double rb = radBound[b] + normv(P[b]-p);
    double r=max(ra,rb);
    if(r > tolSafe) return;
    double c=q.eval(p);
    if(c<0 && c>-1e-12) c=0;
    if(!isfinite(c)) return;
    c += 1e-30*(r/(tolH+1e-30)) + 1e-32*norm2(P[a]-P[b]);
    choices.push_back({p,c,r});
}

static bool candidatePosition(int a,int b, PosChoice &best){
    if(!vActive[a]||!vActive[b]||a==b) return false;
    Quadric q=Q[a]; q.add(Q[b]);
    vector<PosChoice> choices; choices.reserve(7);
    addChoice(choices,q,a,b,P[a]);
    addChoice(choices,q,a,b,P[b]);
    addChoice(choices,q,a,b,(P[a]+P[b])*0.5);
    int ca=clusterSize[a], cb=clusterSize[b];
    addChoice(choices,q,a,b,(P[a]*(double)ca + P[b]*(double)cb)/(double)(ca+cb));
    Vec3 d=P[b]-P[a];
    double f0=q.eval(P[a]); double f1=q.eval(P[b]); double fm=q.eval((P[a]+P[b])*0.5);
    double m=2.0*(f1+f0-2.0*fm); double l=f1-f0-m;
    if(fabs(m)>1e-30){ double t=-l/(2.0*m); if(t>0.0 && t<1.0) addChoice(choices,q,a,b,P[a]+d*t); }
    Vec3 opt;
    if(solve3x3(q,opt)) addChoice(choices,q,a,b,opt);
    if(choices.empty()) return false;
    int bi=0; for(int i=1;i<(int)choices.size();++i) if(choices[i].cost < choices[bi].cost) bi=i;
    best=choices[bi]; return true;
}

static double estimateCost(int a,int b){
    PosChoice pc; if(!candidatePosition(a,b,pc)) return INF_COST;
    return pc.cost;
}

struct Cand{
    float cost;
    int from,to;
    int vf,vt;
};
struct CandGreater{
    bool operator()(const Cand& a,const Cand& b) const { return a.cost > b.cost; }
};

static inline Cand makeCandidate(int u,int v){
    Cand c; c.cost=INFINITY; c.from=u; c.to=v; c.vf=c.vt=-1;
    if(u==v || !vActive[u] || !vActive[v]) return c;
    int from=u, to=v;
    if(clusterSize[from] > clusterSize[to] || (clusterSize[from]==clusterSize[to] && deg0[from] > deg0[to])) swap(from,to);
    double cost=estimateCost(from,to);
    if(cost>=INF_COST/2) return c;
    c.cost=(float)cost; c.from=from; c.to=to; c.vf=versionV[from]; c.vt=versionV[to];
    return c;
}

static inline bool sameUnorderedFace(const Face& f, int a,int b,int c){
    int cnt=0;
    cnt += (f.v[0]==a || f.v[1]==a || f.v[2]==a);
    cnt += (f.v[0]==b || f.v[1]==b || f.v[2]==b);
    cnt += (f.v[0]==c || f.v[1]==c || f.v[2]==c);
    return cnt==3;
}

struct CollapseInfo{
    PosChoice pc;
    vector<int> fromFaces, toFaces, nFrom, nTo;
};

static bool faceNewNormalAfterCollapse(int fid,int from,int to,const Vec3& newp, Vec3& nn){
    Face &f=faces[fid];
    Vec3 p[3];
    for(int k=0;k<3;k++){
        int id=f.v[k];
        if(id==from || id==to) p[k]=newp; else p[k]=P[id];
    }
    if((f.v[0]==from||f.v[0]==to) && (f.v[1]==from||f.v[1]==to)) return false;
    if((f.v[1]==from||f.v[1]==to) && (f.v[2]==from||f.v[2]==to)) return false;
    if((f.v[2]==from||f.v[2]==to) && (f.v[0]==from||f.v[0]==to)) return false;
    Vec3 cr=crossv(p[1]-p[0],p[2]-p[0]);
    double len=normv(cr);
    if(!(len>1e-16)) return false;
    nn=cr/len;
    return true;
}

static bool validChoiceGeometry(int from,int to,const PosChoice& pc, CollapseInfo& info){
    for(int fid: info.fromFaces){
        Face &f=faces[fid];
        if(faceHas(f,to)) continue;
        int a=-1,b=-1;
        for(int k=0;k<3;k++) if(f.v[k]!=from){ if(a<0) a=f.v[k]; else b=f.v[k]; }
        if(a<0||b<0||a==b||a==to||b==to) return false;
        for(int gid: info.toFaces){
            if(gid==fid) continue;
            Face &g=faces[gid]; if(!g.active) continue;
            if(faceHas(g,from) && faceHas(g,to)) continue;
            if(sameUnorderedFace(g,to,a,b)) return false;
        }
    }
    const double minDot = 0.0;
    bumpStampF();
    for(int pass=0; pass<2; ++pass){
        const vector<int>& list = (pass==0?info.fromFaces:info.toFaces);
        for(int fid: list){
            if(markF[fid]==stampF) continue;
            markF[fid]=stampF;
            Face &f=faces[fid]; if(!f.active) continue;
            if(faceHas(f,from) && faceHas(f,to)) continue;
            Vec3 nn;
            if(!faceNewNormalAfterCollapse(fid,from,to,pc.p,nn)) return false;
            double od = nn.x*(double)f.nx + nn.y*(double)f.ny + nn.z*(double)f.nz;
            if(od < minDot) return false;
        }
    }
    return true;
}

static bool checkCollapse(int from,int to, CollapseInfo& info){
    if(from==to || !vActive[from] || !vActive[to]) return false;
    collectData(from, info.fromFaces, info.nFrom);
    collectData(to, info.toFaces, info.nTo);
    int edgeCnt=0;
    for(int fid: info.fromFaces){ Face &f=faces[fid]; if(faceHas(f,to)) edgeCnt++; }
    if(edgeCnt!=2) return false;
    bumpStampV();
    for(int x: info.nTo) if(x!=from) markV[x]=stampV;
    int common=0;
    for(int x: info.nFrom) if(x!=to && markV[x]==stampV) common++;
    if(common!=2) return false;
    Quadric q=Q[from]; q.add(Q[to]);
    vector<PosChoice> choices; choices.reserve(7);
    addChoice(choices,q,from,to,P[from]);
    addChoice(choices,q,from,to,P[to]);
    addChoice(choices,q,from,to,(P[from]+P[to])*0.5);
    int ca=clusterSize[from], cb=clusterSize[to];
    addChoice(choices,q,from,to,(P[from]*(double)ca + P[to]*(double)cb)/(double)(ca+cb));
    Vec3 d=P[to]-P[from];
    double f0=q.eval(P[from]); double f1=q.eval(P[to]); double fm=q.eval((P[from]+P[to])*0.5);
    double m=2.0*(f1+f0-2.0*fm); double l=f1-f0-m;
    if(fabs(m)>1e-30){ double t=-l/(2.0*m); if(t>0.0 && t<1.0) addChoice(choices,q,from,to,P[from]+d*t); }
    Vec3 opt; if(solve3x3(q,opt)) addChoice(choices,q,from,to,opt);
    if(choices.empty()) return false;
    sort(choices.begin(), choices.end(), [](const PosChoice&a,const PosChoice&b){return a.cost<b.cost;});
    for(const auto &pc: choices){
        if(validChoiceGeometry(from,to,pc,info)){ info.pc=pc; return true; }
    }
    return false;
}

static void commitCollapse(int from,int to, CollapseInfo& info){
    Vec3 newp=info.pc.p;
    for(int fid: info.fromFaces){
        Face &f=faces[fid]; if(!f.active) continue;
        if(faceHas(f,to)) { f.active=0; continue; }
        for(int k=0;k<3;k++) if(f.v[k]==from) f.v[k]=to;
        appendIncident(to,fid);
    }
    P[to]=newp;
    radBound[to]=info.pc.rad;
    Q[to].add(Q[from]);
    clusterSize[to]+=clusterSize[from];
    vActive[from]=0;
    activeVertices--;
    versionV[from]++;
    versionV[to]++;
    bumpStampF();
    for(int pass=0; pass<2; ++pass){
        const vector<int>& list=(pass==0?info.fromFaces:info.toFaces);
        for(int fid: list){
            if(markF[fid]==stampF) continue;
            markF[fid]=stampF;
            if(faces[fid].active) recomputeFace(fid);
        }
    }
}

static int chooseTarget(int n){
    if(n<=4) return n;
    if(n<=20) return max(4,n-1);
    double r;
    if(n<=1000) r=0.45;
    else if(n<=7000) r=0.27;
    else if(n<=30000) r=0.06;
    else if(n<=45000) r=0.055;
    else if(n<=70000) r=0.045;
    else if(n<=500000) r=0.025;
    else r=0.020;
    int t=(int)ceil(n*r);
    if(n>1000) t=max(t, 300);
    if(n>7000) t=max(t, 900);
    if(n>30000) t=max(t, 2200);
    if(n>70000) t=max(t, 8000);
    if(n>500000) t=max(t, 25000);
    t=max(4,min(n,t));
    return t;
}

static void simplify(){
    int target=chooseTarget(NV);
    if(activeVertices<=target) return;
    vector<Cand> heap; heap.reserve(initialEdges.size()+1024);
    for(unsigned long long key: initialEdges){
        int u=(int)(key>>32), v=(int)(key & 0xffffffffu);
        Cand c=makeCandidate(u,v);
        if(isfinite(c.cost)) heap.push_back(c);
    }
    vector<unsigned long long>().swap(initialEdges);
    make_heap(heap.begin(), heap.end(), CandGreater());
    vector<int> ftmp, ntmp;
    CollapseInfo info;
    info.fromFaces.reserve(64); info.toFaces.reserve(64); info.nFrom.reserve(64); info.nTo.reserve(64);
    long long pops=0, maxPops=(long long)max(1000000LL, 80LL*NF + 200LL*NV);
    const float planarContinueCost = 1e-18f;
    while(!heap.empty()){
        pop_heap(heap.begin(), heap.end(), CandGreater());
        Cand c=heap.back(); heap.pop_back();
        if(activeVertices<=target && (!isfinite(c.cost) || c.cost>planarContinueCost)) break;
        ++pops; if(pops>maxPops) break;
        if(!isfinite(c.cost)) continue;
        int from=c.from, to=c.to;
        if(from<0||from>=NV||to<0||to>=NV) continue;
        if(!vActive[from] || !vActive[to]) continue;
        if(c.vf!=versionV[from] || c.vt!=versionV[to]) continue;
        if(!checkCollapse(from,to,info)) continue;
        commitCollapse(from,to,info);
        collectData(to, ftmp, ntmp);
        for(int nb: ntmp){
            if(nb==to || !vActive[nb]) continue;
            Cand nc=makeCandidate(to,nb);
            if(isfinite(nc.cost)){ heap.push_back(nc); push_heap(heap.begin(), heap.end(), CandGreater()); }
        }
        if(heap.size() > (size_t)max(5000000, 18*activeVertices + 100000)){
            vector<Cand> nh; nh.reserve(heap.size()/2);
            for(const Cand &cc: heap){
                if(cc.from>=0&&cc.to>=0&&cc.from<NV&&cc.to<NV&&vActive[cc.from]&&vActive[cc.to]
                   && cc.vf==versionV[cc.from]&&cc.vt==versionV[cc.to] && isfinite(cc.cost)) nh.push_back(cc);
            }
            heap.swap(nh); make_heap(heap.begin(), heap.end(), CandGreater());
        }
    }
}

struct FastOut{
    string s;
    FastOut(){ s.reserve(1<<20); }
    ~FastOut(){ flush(); }
    inline void flush(){ if(!s.empty()){ fwrite(s.data(),1,s.size(),stdout); s.clear(); } }
    inline void append(const char* buf,int n){ if(s.size()+n > (1<<20)) flush(); s.append(buf,n); }
};

static void save_output(){
    vector<int> id(NV,-1);
    int vc=0, fc=0;
    for(int i=0;i<NV;i++) if(vActive[i]) id[i]=++vc;
    for(int i=0;i<NF;i++) if(faces[i].active) fc++;
    FastOut out; char line[128];
    out.append(line, snprintf(line,sizeof(line), "%d %d\n", vc, fc));
    for(int i=0;i<NV;i++) if(id[i]>0){
        out.append(line, snprintf(line,sizeof(line), "v %.10g %.10g %.10g\n", P[i].x,P[i].y,P[i].z));
    }
    for(int i=0;i<NF;i++) if(faces[i].active){
        Face &f=faces[i];
        int a=id[f.v[0]], b=id[f.v[1]], c=id[f.v[2]];
        out.append(line, snprintf(line,sizeof(line), "f %d %d %d\n", a,b,c));
    }
    out.flush();
}

int main(){
    load_input();
    simplify();
    save_output();
    return 0;
}
