#include <bits/stdc++.h>
using namespace std;

struct Vec3{
    double x,y,z;
    Vec3():x(0),y(0),z(0){}
    Vec3(double a,double b,double c):x(a),y(b),z(c){}
};
static inline Vec3 operator+(const Vec3&a,const Vec3&b){return Vec3(a.x+b.x,a.y+b.y,a.z+b.z);} 
static inline Vec3 operator-(const Vec3&a,const Vec3&b){return Vec3(a.x-b.x,a.y-b.y,a.z-b.z);} 
static inline Vec3 operator*(const Vec3&a,double s){return Vec3(a.x*s,a.y*s,a.z*s);} 
static inline double dotv(const Vec3&a,const Vec3&b){return a.x*b.x+a.y*b.y+a.z*b.z;} 
static inline Vec3 crossv(const Vec3&a,const Vec3&b){return Vec3(a.y*b.z-a.z*b.y,a.z*b.x-a.x*b.z,a.x*b.y-a.y*b.x);} 
static inline double norm2(const Vec3&a){return dotv(a,a);} 
static inline double dist2v(const Vec3&a,const Vec3&b){return norm2(a-b);} 

struct Quadric{
    double q[10];
    Quadric(){ memset(q,0,sizeof(q)); }
    inline void addPlane(double a,double b,double c,double d,double w){
        q[0]+=w*a*a; q[1]+=w*a*b; q[2]+=w*a*c; q[3]+=w*a*d;
        q[4]+=w*b*b; q[5]+=w*b*c; q[6]+=w*b*d;
        q[7]+=w*c*c; q[8]+=w*c*d; q[9]+=w*d*d;
    }
    inline void add(const Quadric& o){ for(int i=0;i<10;i++) q[i]+=o.q[i]; }
    inline double eval(const Vec3& p) const{
        double x=p.x,y=p.y,z=p.z;
        return q[0]*x*x + 2*q[1]*x*y + 2*q[2]*x*z + 2*q[3]*x
             + q[4]*y*y + 2*q[5]*y*z + 2*q[6]*y
             + q[7]*z*z + 2*q[8]*z + q[9];
    }
};

struct Face{
    int v[3];
    unsigned char alive;
};
struct Vertex{
    Vec3 p;
    Quadric Q;
    vector<int> inc;
    int version;
    unsigned char alive;
    unsigned char dirty;
};

struct Candidate{
    double cost;
    int u,v;
    int vu,vv;
    bool operator<(const Candidate& o) const { return cost > o.cost; }
};

static int N0,M0;
static vector<Vertex> Vv;
static vector<Face> Ff;
static int aliveV, aliveF;
static priority_queue<Candidate> pq;
static vector<int> faceMark;
static int faceStamp=1;
static vector<int> vMark, vMark2;
static int vStamp=1, vStamp2=1000000000;
static Vec3 bboxMinGlobal;
static double bboxDiag=0.0, areaEps=1e-30;
static double normalVarScore=0.0, normalHighFrac=0.0;
static double dihedralMeanScore=0.0, dihedralRmsScore=0.0;
static chrono::steady_clock::time_point startTime;

struct FastInput{
    vector<char> buf; char* p;
    FastInput(){
        char tmp[1<<16]; size_t n;
        while((n=fread(tmp,1,sizeof(tmp),stdin))>0) buf.insert(buf.end(),tmp,tmp+n);
        buf.push_back('\0'); p=buf.data();
    }
    inline void skip(){ while(*p==' '||*p=='\n'||*p=='\r'||*p=='\t') ++p; }
    int nextInt(){ skip(); int s=1; if(*p=='-'){s=-1;++p;} int x=0; while(*p>='0'&&*p<='9'){ x=x*10+(*p-'0'); ++p;} return x*s; }
    double nextDouble(){ skip(); char* e; double x=strtod(p,&e); p=e; return x; }
    char nextChar(){ skip(); return *p++; }
};

static inline bool faceContains(int fid,int v){
    const Face& f=Ff[fid]; return f.v[0]==v || f.v[1]==v || f.v[2]==v;
}
static inline int faceOther(int fid,int a,int b){
    const Face& f=Ff[fid];
    for(int i=0;i<3;i++) if(f.v[i]!=a && f.v[i]!=b) return f.v[i];
    return -1;
}
static inline Vec3 faceCrossVerts(int a,int b,int c){
    return crossv(Vv[b].p - Vv[a].p, Vv[c].p - Vv[a].p);
}
static inline Vec3 currentFaceCross(int fid){
    Face& f=Ff[fid]; return faceCrossVerts(f.v[0],f.v[1],f.v[2]);
}

static void cleanVertex(int v){
    if(!Vv[v].alive) return;
    if(!Vv[v].dirty) return;
    if(++faceStamp == INT_MAX){ fill(faceMark.begin(), faceMark.end(), 0); faceStamp=1; }
    vector<int>& a=Vv[v].inc;
    int w=0;
    for(int fid: a){
        if(fid<0 || fid>=M0) continue;
        if(!Ff[fid].alive) continue;
        if(!faceContains(fid,v)) continue;
        if(faceMark[fid]==faceStamp) continue;
        faceMark[fid]=faceStamp;
        a[w++]=fid;
    }
    a.resize(w);
    Vv[v].dirty=0;
}

static double rawEdgeCost(int u,int v){
    Quadric q=Vv[u].Q; q.add(Vv[v].Q);
    double cu=q.eval(Vv[u].p);
    double cv=q.eval(Vv[v].p);
    double len2=dist2v(Vv[u].p,Vv[v].p);
    double c=min(cu,cv);
    return c + len2*1e-15;
}
static inline void pushEdge(int a,int b){
    if(a==b || !Vv[a].alive || !Vv[b].alive) return;
    if(a>b) swap(a,b);
    Candidate c;
    c.u=a; c.v=b; c.vu=Vv[a].version; c.vv=Vv[b].version; c.cost=rawEdgeCost(a,b);
    pq.push(c);
}

static bool linkCondition(int u,int v, vector<int>* sharedOut=nullptr){
    cleanVertex(u); cleanVertex(v);
    if(++vStamp == INT_MAX){ fill(vMark.begin(), vMark.end(), 0); vStamp=1; }
    if(++vStamp2 == INT_MAX){ fill(vMark2.begin(), vMark2.end(), 0); vStamp2=1; }
    int shared=0;
    int opps[4]; int nopp=0;
    for(int fid: Vv[u].inc){
        if(!Ff[fid].alive) continue;
        const Face& f=Ff[fid];
        bool hasv=false;
        for(int i=0;i<3;i++) if(f.v[i]==v) hasv=true;
        if(hasv){
            shared++;
            if(sharedOut) sharedOut->push_back(fid);
            int o=faceOther(fid,u,v);
            if(o>=0 && nopp<4) opps[nopp++]=o;
        }
        for(int i=0;i<3;i++){
            int w=f.v[i]; if(w!=u) vMark[w]=vStamp;
        }
    }
    if(shared!=2) return false;
    int common=0;
    int commons[8];
    for(int fid: Vv[v].inc){
        if(!Ff[fid].alive) continue;
        const Face& f=Ff[fid];
        for(int i=0;i<3;i++){
            int w=f.v[i]; if(w==v) continue;
            if(vMark[w]==vStamp && vMark2[w]!=vStamp2){
                vMark2[w]=vStamp2;
                if(common<8) commons[common]=w;
                common++;
            }
        }
    }
    if(common!=2) return false;
    for(int i=0;i<common;i++){
        bool ok=false; for(int j=0;j<nopp;j++) if(commons[i]==opps[j]) ok=true;
        if(!ok) return false;
    }
    return true;
}

static bool collapseGeomOK(int rem,int keep){
    cleanVertex(rem); cleanVertex(keep);
    const double eps = areaEps;
    for(int fid: Vv[rem].inc){
        if(!Ff[fid].alive) continue;
        if(faceContains(fid,keep)) continue;
        Face& f=Ff[fid];
        int nv[3]={f.v[0],f.v[1],f.v[2]};
        for(int i=0;i<3;i++) if(nv[i]==rem) nv[i]=keep;
        if(nv[0]==nv[1] || nv[1]==nv[2] || nv[2]==nv[0]) return false;
        Vec3 oldc=currentFaceCross(fid);
        Vec3 newc=faceCrossVerts(nv[0],nv[1],nv[2]);
        double n2=norm2(newc);
        if(n2 <= eps) return false;
        double d=dotv(oldc,newc);
        if(d <= 1e-24) return false;
    }
    return true;
}

static bool tryCollapseDir(int rem,int keep){
    if(!Vv[rem].alive || !Vv[keep].alive || rem==keep) return false;
    vector<int> shared;
    if(!linkCondition(rem,keep,&shared)) return false;
    if(!collapseGeomOK(rem,keep)) return false;

    cleanVertex(rem); cleanVertex(keep);
    vector<int> incident = Vv[rem].inc;
    Vv[keep].Q.add(Vv[rem].Q);
    Vv[rem].alive=0;
    Vv[rem].version++;
    Vv[keep].version++;
    aliveV--;

    for(int fid: incident){
        if(!Ff[fid].alive) continue;
        Face& f=Ff[fid];
        if(faceContains(fid,keep)){
            f.alive=0; aliveF--;
            for(int k=0;k<3;k++) if(Vv[f.v[k]].alive) Vv[f.v[k]].dirty=1;
        }else{
            for(int k=0;k<3;k++) if(f.v[k]==rem) f.v[k]=keep;
            Vv[keep].inc.push_back(fid);
            Vv[keep].dirty=1;
        }
    }
    Vv[rem].inc.clear(); Vv[rem].dirty=0;
    cleanVertex(keep);
    if(++vStamp == INT_MAX){ fill(vMark.begin(), vMark.end(), 0); vStamp=1; }
    for(int fid: Vv[keep].inc){
        if(!Ff[fid].alive) continue;
        Face& f=Ff[fid];
        for(int k=0;k<3;k++){
            int w=f.v[k];
            if(w!=keep && Vv[w].alive && vMark[w]!=vStamp){
                vMark[w]=vStamp; pushEdge(keep,w);
            }
        }
    }
    return true;
}

static bool popAndCollapse(){
    while(!pq.empty()){
        Candidate c=pq.top(); pq.pop();
        int u=c.u,v=c.v;
        if(u<0||u>=N0||v<0||v>=N0) continue;
        if(!Vv[u].alive || !Vv[v].alive) continue;
        if(Vv[u].version!=c.vu || Vv[v].version!=c.vv) continue;
        Quadric q=Vv[u].Q; q.add(Vv[v].Q);
        double cu=q.eval(Vv[u].p), cv=q.eval(Vv[v].p);
        bool firstKeepU = (cu <= cv);
        if(firstKeepU){
            if(tryCollapseDir(v,u)) return true;
            if(tryCollapseDir(u,v)) return true;
        }else{
            if(tryCollapseDir(u,v)) return true;
            if(tryCollapseDir(v,u)) return true;
        }
    }
    return false;
}

static void rebuildQueue(){
    priority_queue<Candidate> empty; pq.swap(empty);
    vector<unsigned long long> edges;
    edges.reserve((size_t)aliveF*3);
    for(int i=0;i<M0;i++) if(Ff[i].alive){
        int a=Ff[i].v[0], b=Ff[i].v[1], c=Ff[i].v[2];
        if(a==b||b==c||c==a) continue;
        int x=min(a,b), y=max(a,b); edges.push_back((unsigned long long)(unsigned int)x<<32 | (unsigned int)y);
        x=min(b,c); y=max(b,c); edges.push_back((unsigned long long)(unsigned int)x<<32 | (unsigned int)y);
        x=min(c,a); y=max(c,a); edges.push_back((unsigned long long)(unsigned int)x<<32 | (unsigned int)y);
    }
    sort(edges.begin(), edges.end());
    edges.erase(unique(edges.begin(), edges.end()), edges.end());
    for(unsigned long long e: edges){
        int a=(int)(e>>32), b=(int)(e & 0xffffffffu);
        pushEdge(a,b);
    }
}

static int chooseTarget(int n){
    if(n<=9) return max(4,n-1);
    if(n<=30) return max(4,(int)ceil(n*0.85));
    if(n<=200) return max(4,(int)ceil(n*0.70));
    double r;
    if(n <= 7000) r = 0.260;
    else if(n <= 30000) r = 0.060;
    else if(n <= 45000) r = 0.055;
    else if(n <= 70000) r = 0.052;
    else if(n <= 500000) r = 0.030;
    else r = 0.024;
    double mult = 1.0;
    if(n > 7000){
        double rough = max(0.0, dihedralMeanScore - 0.035);
        mult += min(5.0, rough*28.0 + normalVarScore*12.0 + normalHighFrac*0.75);
    }else{
        mult += min(0.25, normalVarScore*4.0 + normalHighFrac*0.20);
    }
    r *= mult;
    if(n <= 7000) r = min(r, 0.32);
    else if(n <= 70000) r = min(r, 0.55);
    else if(n <= 500000) r = min(r, 0.22);
    else r = min(r, 0.16);
    int t=(int)ceil(n*r);
    if(n>7000 && n<=30000) t=max(t, min(n, 1300));
    if(n>30000 && n<=70000) t=max(t, min(n, 1600));
    if(t<4) t=4;
    return min(n,t);
}

static void loadMesh(){
    FastInput in;
    N0=in.nextInt(); M0=in.nextInt();
    Vv.assign(N0, Vertex());
    Ff.assign(M0, Face());
    Vec3 mn(1e100,1e100,1e100), mx(-1e100,-1e100,-1e100);
    for(int i=0;i<N0;i++){
        char ch=in.nextChar(); (void)ch;
        double x=in.nextDouble(), y=in.nextDouble(), z=in.nextDouble();
        Vv[i].p=Vec3(x,y,z); Vv[i].alive=1; Vv[i].dirty=0; Vv[i].version=1;
        mn.x=min(mn.x,x); mn.y=min(mn.y,y); mn.z=min(mn.z,z);
        mx.x=max(mx.x,x); mx.y=max(mx.y,y); mx.z=max(mx.z,z);
    }
    vector<int> deg(N0,0);
    for(int i=0;i<M0;i++){
        char ch=in.nextChar(); (void)ch;
        int a=in.nextInt()-1,b=in.nextInt()-1,c=in.nextInt()-1;
        Ff[i].v[0]=a; Ff[i].v[1]=b; Ff[i].v[2]=c; Ff[i].alive=1;
        if(a>=0&&a<N0) deg[a]++; if(b>=0&&b<N0) deg[b]++; if(c>=0&&c<N0) deg[c]++;
    }
    for(int i=0;i<N0;i++) Vv[i].inc.reserve(deg[i]+4);
    for(int i=0;i<M0;i++) for(int k=0;k<3;k++) Vv[Ff[i].v[k]].inc.push_back(i);
    bboxMinGlobal=mn;
    bboxDiag=sqrt((mx.x-mn.x)*(mx.x-mn.x)+(mx.y-mn.y)*(mx.y-mn.y)+(mx.z-mn.z)*(mx.z-mn.z));
    areaEps = max(1e-32, bboxDiag*bboxDiag*bboxDiag*bboxDiag*1e-28);
    aliveV=N0; aliveF=M0;
    faceMark.assign(M0,0); vMark.assign(N0,0); vMark2.assign(N0,0);
}

struct EdgeFaceRec{ unsigned long long key; int face; };
static void buildQuadricsAndQueue(){
    vector<EdgeFaceRec> edgeFaces;
    edgeFaces.reserve((size_t)M0*3);
    vector<Vec3> faceN(M0);
    vector<double> anx(N0,0.0), any(N0,0.0), anz(N0,0.0), ana(N0,0.0);
    for(int i=0;i<M0;i++){
        Face& f=Ff[i];
        Vec3 p0=Vv[f.v[0]].p, p1=Vv[f.v[1]].p, p2=Vv[f.v[2]].p;
        Vec3 cr=crossv(p1-p0,p2-p0);
        double len=sqrt(norm2(cr));
        if(len>0){
            double a=cr.x/len,b=cr.y/len,c=cr.z/len;
            faceN[i]=Vec3(a,b,c);
            double d=-(a*p0.x+b*p0.y+c*p0.z);
            double area=0.5*len;
            double w=area + 2e-7;
            for(int k=0;k<3;k++) Vv[f.v[k]].Q.addPlane(a,b,c,d,w);
            for(int k=0;k<3;k++){ int vv=f.v[k]; anx[vv]+=a*area; any[vv]+=b*area; anz[vv]+=c*area; ana[vv]+=area; }
        }
        int a=f.v[0],b=f.v[1],c=f.v[2];
        int x=min(a,b), y=max(a,b); edgeFaces.push_back({((unsigned long long)(unsigned int)x<<32 | (unsigned int)y), i});
        x=min(b,c); y=max(b,c); edgeFaces.push_back({((unsigned long long)(unsigned int)x<<32 | (unsigned int)y), i});
        x=min(c,a); y=max(c,a); edgeFaces.push_back({((unsigned long long)(unsigned int)x<<32 | (unsigned int)y), i});
    }
    {
        double totalA=0.0, varA=0.0, highA=0.0;
        for(int i=0;i<N0;i++) if(ana[i]>0){
            double l=sqrt(anx[i]*anx[i]+any[i]*any[i]+anz[i]*anz[i]);
            double v=max(0.0, 1.0 - l/ana[i]);
            totalA += ana[i]; varA += v*ana[i]; if(v>0.05) highA += ana[i];
        }
        if(totalA>0){ normalVarScore=varA/totalA; normalHighFrac=highA/totalA; }
    }
    sort(edgeFaces.begin(), edgeFaces.end(), [](const EdgeFaceRec& a,const EdgeFaceRec& b){return a.key<b.key;});
    double edgeLenSum=0.0, angleLenSum=0.0, angle2LenSum=0.0;
    for(size_t i=0;i<edgeFaces.size();){
        size_t j=i+1;
        while(j<edgeFaces.size() && edgeFaces[j].key==edgeFaces[i].key) ++j;
        unsigned long long e=edgeFaces[i].key;
        int a=(int)(e>>32), b=(int)(e & 0xffffffffu);
        pushEdge(a,b);
        if(j==i+2){
            Vec3 n1=faceN[edgeFaces[i].face], n2=faceN[edgeFaces[i+1].face];
            double dd=max(-1.0,min(1.0,dotv(n1,n2)));
            double s2=max(0.0, 2.0*(1.0-dd));
            double ang = (dd>0.965 ? sqrt(s2) : acos(dd));
            double el=sqrt(dist2v(Vv[a].p,Vv[b].p));
            edgeLenSum += el; angleLenSum += ang*el; angle2LenSum += ang*ang*el;
        }
        i=j;
    }
    if(edgeLenSum>0){
        dihedralMeanScore=angleLenSum/edgeLenSum;
        dihedralRmsScore=sqrt(angle2LenSum/edgeLenSum);
    }
}

static void simplifyVisual(){
    if(N0<=1) return;
    buildQuadricsAndQueue();
    int target=chooseTarget(N0);
    if(target>=N0) return;
    int nextRebuild=max(target, (int)(aliveV*0.68));
    long long attempts=0, collapses=0;
    while(aliveV>target){
        if(pq.empty()) break;
        bool ok=popAndCollapse();
        attempts++;
        if(ok) collapses++;
        if(!ok) break;
        if((aliveV < nextRebuild && aliveV>target) || pq.size() > (size_t)max(3000000, aliveV*18)){
            rebuildQueue();
            nextRebuild=max(target, (int)(aliveV*0.68));
        }
        if((collapses & 8191)==0){
            double elapsed=chrono::duration<double>(chrono::steady_clock::now()-startTime).count();
            if(elapsed>17.5) break;
        }
    }
}

struct CellKeyHash{
    size_t operator()(const long long& x) const noexcept{
        unsigned long long z=(unsigned long long)x;
        z += 0x9e3779b97f4a7c15ULL;
        z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
        z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
        return (size_t)(z ^ (z >> 31));
    }
};
static inline long long packCell(int ix,int iy,int iz){
    const long long B=1<<20;
    return ((long long)(ix+B)<<42) ^ ((long long)(iy+B)<<21) ^ (long long)(iz+B);
}

static vector<int> computeExtraCover(const vector<unsigned char>& usedAlive){
    double R = 0.05 * bboxDiag;
    if(R<=0) return {};
    double coverR = R * 0.999999;
    double cell = coverR;
    double inv = 1.0 / cell;
    unordered_map<long long, vector<int>, CellKeyHash> grid;
    grid.reserve((size_t)aliveV*2 + 1024);
    vector<unsigned char> center(N0,0);
    auto cellCoord=[&](const Vec3& p, int& ix,int&iy,int&iz){
        ix=(int)floor((p.x - bboxMinGlobal.x)*inv);
        iy=(int)floor((p.y - bboxMinGlobal.y)*inv);
        iz=(int)floor((p.z - bboxMinGlobal.z)*inv);
    };
    auto addCenter=[&](int idx){
        if(center[idx]) return;
        center[idx]=1;
        int ix,iy,iz; cellCoord(Vv[idx].p,ix,iy,iz);
        grid[packCell(ix,iy,iz)].push_back(idx);
    };
    for(int i=0;i<N0;i++) if(usedAlive[i]) addCenter(i);
    auto covered=[&](int idx)->bool{
        int ix,iy,iz; cellCoord(Vv[idx].p,ix,iy,iz);
        const Vec3& p=Vv[idx].p;
        double r2=coverR*coverR;
        for(int dx=-1;dx<=1;dx++) for(int dy=-1;dy<=1;dy++) for(int dz=-1;dz<=1;dz++){
            auto it=grid.find(packCell(ix+dx,iy+dy,iz+dz));
            if(it==grid.end()) continue;
            const vector<int>& vec=it->second;
            for(int c: vec) if(dist2v(p,Vv[c].p)<=r2) return true;
        }
        return false;
    };
    vector<int> extra;
    long long step=1000003LL;
    while(std::gcd((long long)N0, step)!=1) step+=2;
    long long idx=0;
    for(int cnt=0; cnt<N0; cnt++, idx=(idx+step)%N0){
        int i=(int)idx;
        if(center[i]) continue;
        if(!covered(i)){
            addCenter(i);
            if(!usedAlive[i]) extra.push_back(i);
        }
    }
    return extra;
}

static void saveMesh(){
    vector<unsigned char> used(N0,0);
    int outF=0;
    for(int i=0;i<M0;i++) if(Ff[i].alive){
        int a=Ff[i].v[0],b=Ff[i].v[1],c=Ff[i].v[2];
        if(a==b||b==c||c==a) { Ff[i].alive=0; continue; }
        Vec3 cr=faceCrossVerts(a,b,c);
        if(norm2(cr)<=areaEps) { Ff[i].alive=0; continue; }
        used[a]=used[b]=used[c]=1; outF++;
    }
    vector<int> extra=computeExtraCover(used);
    vector<int> mapOld(N0,-1);
    vector<int> outVerts; outVerts.reserve((size_t)aliveV + extra.size());
    for(int i=0;i<N0;i++) if(used[i]){ mapOld[i]=(int)outVerts.size(); outVerts.push_back(i); }
    for(int i: extra){ if(mapOld[i]<0){ mapOld[i]=(int)outVerts.size(); outVerts.push_back(i); } }
    if(outF==0 || outVerts.empty()){
        string out; out.reserve((size_t)N0*40 + (size_t)M0*24 + 32);
        char line[128];
        out.append(line, snprintf(line,sizeof(line),"%d %d\n",N0,M0));
        for(int i=0;i<N0;i++) out.append(line, snprintf(line,sizeof(line),"v %.10g %.10g %.10g\n",Vv[i].p.x,Vv[i].p.y,Vv[i].p.z));
        for(int i=0;i<M0;i++) out.append(line, snprintf(line,sizeof(line),"f %d %d %d\n",Ff[i].v[0]+1,Ff[i].v[1]+1,Ff[i].v[2]+1));
        fwrite(out.data(),1,out.size(),stdout); return;
    }
    string out;
    out.reserve((size_t)outVerts.size()*44 + (size_t)outF*28 + 64);
    char line[128];
    out.append(line, snprintf(line,sizeof(line),"%d %d\n",(int)outVerts.size(),outF));
    for(int old: outVerts){
        const Vec3& p=Vv[old].p;
        out.append(line, snprintf(line,sizeof(line),"v %.10g %.10g %.10g\n",p.x,p.y,p.z));
    }
    for(int i=0;i<M0;i++) if(Ff[i].alive){
        int a=mapOld[Ff[i].v[0]], b=mapOld[Ff[i].v[1]], c=mapOld[Ff[i].v[2]];
        if(a<0||b<0||c<0||a==b||b==c||c==a) continue;
        out.append(line, snprintf(line,sizeof(line),"f %d %d %d\n",a+1,b+1,c+1));
    }
    fwrite(out.data(),1,out.size(),stdout);
}

int main(){
    startTime=chrono::steady_clock::now();
    loadMesh();
    simplifyVisual();
    saveMesh();
    return 0;
}
