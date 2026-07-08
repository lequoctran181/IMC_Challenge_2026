// IMC Challenge 2026 - simplifygeometry
// C++17 single-file candidate
//
// Hybrid high-score attempt.  The main mechanism is still an adaptive r-net
// over the input vertices, but with two changes aimed at escaping the 81.x
// plateau: (1) local radii may be larger than the global radius on smooth
// low-curvature surface regions, while boundary/sharp vertices keep smaller
// radii; (2) smooth output representatives are mildly re-centered in their
// local tangent plane, which improves visual geometry without spending extra
// vertices.  Faces are obtained by contracting original triangles, then repaired
// to remove degeneracies and >2-valence edges.
//
// Compile: g++ -O2 -std=gnu++17 simplifygeometry.cpp -o simplifygeometry

#include <bits/stdc++.h>
using namespace std;

static constexpr double EPS = 1e-12;
static constexpr double AREA_EPS = 1e-22;

struct Vec3 {
    double x=0, y=0, z=0;
    Vec3() = default;
    Vec3(double X,double Y,double Z):x(X),y(Y),z(Z){}
    Vec3 operator+(const Vec3& o) const { return {x+o.x,y+o.y,z+o.z}; }
    Vec3 operator-(const Vec3& o) const { return {x-o.x,y-o.y,z-o.z}; }
    Vec3 operator*(double s) const { return {x*s,y*s,z*s}; }
    Vec3 operator/(double s) const { return {x/s,y/s,z/s}; }
    Vec3& operator+=(const Vec3& o){ x+=o.x; y+=o.y; z+=o.z; return *this; }
};
static inline double dotv(const Vec3& a,const Vec3& b){ return a.x*b.x+a.y*b.y+a.z*b.z; }
static inline Vec3 crossv(const Vec3& a,const Vec3& b){ return {a.y*b.z-a.z*b.y, a.z*b.x-a.x*b.z, a.x*b.y-a.y*b.x}; }
static inline double norm2(const Vec3& a){ return dotv(a,a); }
static inline double normv(const Vec3& a){ return sqrt(max(0.0,norm2(a))); }
static inline double dist2(const Vec3& a,const Vec3& b){ return norm2(a-b); }
static inline Vec3 normalize(const Vec3& a){ double n=normv(a); return n>EPS ? a/n : Vec3{0,0,0}; }

struct Face { int a=0,b=0,c=0; Face()=default; Face(int A,int B,int C):a(A),b(B),c(C){} };
enum class Format { RAW, OFF, OBJ };
struct Mesh { vector<Vec3> V; vector<Face> F; Format fmt=Format::RAW; int rawBase=0; };

static inline uint64_t splitmix64(uint64_t x){
    x += 0x9e3779b97f4a7c15ULL;
    x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
    x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
    return x ^ (x >> 31);
}

struct TripleKey {
    int a,b,c;
    bool operator==(const TripleKey& o) const noexcept { return a==o.a && b==o.b && c==o.c; }
};
struct TripleHash {
    size_t operator()(const TripleKey& t) const noexcept {
        uint64_t x = (uint64_t)(uint32_t)t.a * 11995408973635179863ULL
                 ^ (uint64_t)(uint32_t)t.b * 10150724397891781847ULL
                 ^ (uint64_t)(uint32_t)t.c *  7809847782465536327ULL;
        return (size_t)splitmix64(x);
    }
};
static inline TripleKey sortedTri(int a,int b,int c){
    if(a>b) swap(a,b); if(b>c) swap(b,c); if(a>b) swap(a,b); return {a,b,c};
}

struct PairKey { int a,b; bool operator==(const PairKey& o) const noexcept { return a==o.a && b==o.b; } };
struct PairHash {
    size_t operator()(const PairKey& p) const noexcept {
        return (size_t)splitmix64(((uint64_t)(uint32_t)p.a<<32) ^ (uint32_t)p.b);
    }
};

struct Bounds { Vec3 mn,mx; double diag=0; };
static Bounds computeBounds(const vector<Vec3>& V){
    Bounds b;
    if(V.empty()) return b;
    b.mn=b.mx=V[0];
    for(const Vec3& p: V){
        b.mn.x=min(b.mn.x,p.x); b.mn.y=min(b.mn.y,p.y); b.mn.z=min(b.mn.z,p.z);
        b.mx.x=max(b.mx.x,p.x); b.mx.y=max(b.mx.y,p.y); b.mx.z=max(b.mx.z,p.z);
    }
    b.diag = normv(b.mx-b.mn);
    return b;
}
static inline bool validFace(const Face& f,int n){
    return f.a>=0 && f.a<n && f.b>=0 && f.b<n && f.c>=0 && f.c<n && f.a!=f.b && f.b!=f.c && f.c!=f.a;
}

// ---------- Parsing ----------
struct Scanner {
    const string& s; size_t p=0;
    explicit Scanner(const string& S):s(S){}
    bool next(string& out){
        out.clear();
        while(p<s.size()){
            unsigned char ch=(unsigned char)s[p];
            if(isspace(ch)){ ++p; continue; }
            if(s[p]=='#'){ while(p<s.size() && s[p]!='\n') ++p; continue; }
            break;
        }
        if(p>=s.size()) return false;
        size_t st=p;
        while(p<s.size() && !isspace((unsigned char)s[p]) && s[p]!='#') ++p;
        out.assign(s.begin()+st, s.begin()+p);
        return true;
    }
    bool nextLL(long long& v){ string t; if(!next(t)) return false; try{ size_t pos=0; v=stoll(t,&pos); return pos==t.size(); }catch(...){ return false; } }
    bool nextInt(int& v){ long long x; if(!nextLL(x) || x<INT_MIN || x>INT_MAX) return false; v=(int)x; return true; }
    bool nextDouble(double& v){ string t; if(!next(t)) return false; try{ size_t pos=0; v=stod(t,&pos); return pos==t.size(); }catch(...){ return false; } }
};
static inline string ltrim_copy(string s){ size_t p=0; while(p<s.size() && isspace((unsigned char)s[p])) ++p; return s.substr(p); }

static bool parseOBJ(const string& data, Mesh& mesh){
    Mesh out; out.fmt=Format::OBJ; out.rawBase=1;
    istringstream in(data); string line;
    while(getline(in,line)){
        string t=ltrim_copy(line);
        if(t.empty() || t[0]=='#') continue;
        if(t.size()>=2 && t[0]=='v' && isspace((unsigned char)t[1])){
            istringstream ss(t.substr(1)); double x,y,z;
            if(ss>>x>>y>>z) out.V.emplace_back(x,y,z);
        }else if(t.size()>=2 && t[0]=='f' && isspace((unsigned char)t[1])){
            istringstream ss(t.substr(1)); string tok; vector<int> ids;
            while(ss>>tok){
                size_t slash=tok.find('/'); string fst = slash==string::npos ? tok : tok.substr(0,slash);
                if(fst.empty()) continue;
                long long val=0; try{ val=stoll(fst); }catch(...){ continue; }
                int id = val<0 ? (int)out.V.size() + (int)val : (int)val - 1;
                if(id>=0 && id<(int)out.V.size()) ids.push_back(id);
            }
            if(ids.size()>=3){
                for(size_t i=1;i+1<ids.size();++i) out.F.emplace_back(ids[0],ids[i],ids[i+1]);
            }
        }
    }
    if(out.V.empty()) return false;
    mesh=std::move(out); return true;
}

static bool parseMesh(const string& data, Mesh& mesh){
    Scanner sc(data); string first;
    if(!sc.next(first)) return false;
    if(first=="OFF"){
        long long nLL=0,mLL=0,eLL=0;
        if(!sc.nextLL(nLL) || !sc.nextLL(mLL) || !sc.nextLL(eLL)) return false;
        if(nLL<=0 || mLL<0 || nLL>INT_MAX || mLL>INT_MAX) return false;
        int n=(int)nLL, m=(int)mLL;
        Mesh out; out.fmt=Format::OFF; out.rawBase=0; out.V.reserve(n); out.F.reserve(m);
        for(int i=0;i<n;++i){ double x,y,z; if(!sc.nextDouble(x)||!sc.nextDouble(y)||!sc.nextDouble(z)) return false; out.V.emplace_back(x,y,z); }
        for(int i=0;i<m;++i){
            int k; if(!sc.nextInt(k)) return false; vector<int> ids; ids.reserve(max(0,k));
            for(int j=0;j<k;++j){ int x; if(!sc.nextInt(x)) return false; ids.push_back(x); }
            if((int)ids.size()>=3){
                for(int j=1;j+1<(int)ids.size();++j){
                    if(ids[0]>=0 && ids[0]<n && ids[j]>=0 && ids[j]<n && ids[j+1]>=0 && ids[j+1]<n)
                        out.F.emplace_back(ids[0],ids[j],ids[j+1]);
                }
            }
        }
        mesh=std::move(out); return true;
    }
    if(first=="v" || first=="f" || (!first.empty() && first[0]=='#')) return parseOBJ(data,mesh);

    // Raw numeric: n m, vertices, triangles.  Preserve 0/1 base.
    long long nLL=0,mLL=0;
    try{ size_t pos=0; nLL=stoll(first,&pos); if(pos!=first.size()) return parseOBJ(data,mesh); }
    catch(...){ return parseOBJ(data,mesh); }
    if(!sc.nextLL(mLL)) return false;
    if(nLL<=0 || mLL<0 || nLL>INT_MAX || mLL>INT_MAX) return false;
    int n=(int)nLL, m=(int)mLL;
    Mesh out; out.fmt=Format::RAW; out.V.reserve(n); out.F.reserve(m);
    for(int i=0;i<n;++i){ double x,y,z; if(!sc.nextDouble(x)||!sc.nextDouble(y)||!sc.nextDouble(z)) return false; out.V.emplace_back(x,y,z); }
    vector<array<long long,3>> rf; rf.reserve(m); long long mn=LLONG_MAX;
    for(int i=0;i<m;++i){ long long a,b,c; if(!sc.nextLL(a)||!sc.nextLL(b)||!sc.nextLL(c)) return false; rf.push_back({a,b,c}); mn=min(mn,min(a,min(b,c))); }
    int base = (mn==0 ? 0 : 1); out.rawBase=base;
    for(auto t: rf){ long long a=t[0]-base,b=t[1]-base,c=t[2]-base; if(a>=0&&a<n&&b>=0&&b<n&&c>=0&&c<n) out.F.emplace_back((int)a,(int)b,(int)c); }
    mesh=std::move(out); return true;
}

static void cleanFacesKeepVertices(Mesh& mesh){
    int n=(int)mesh.V.size();
    vector<Face> nf; nf.reserve(mesh.F.size());
    unordered_set<TripleKey,TripleHash> seen; seen.reserve(mesh.F.size()*2+1);
    for(const Face& f: mesh.F){
        if(!validFace(f,n)) continue;
        Vec3 cr=crossv(mesh.V[f.b]-mesh.V[f.a], mesh.V[f.c]-mesh.V[f.a]);
        if(norm2(cr)<=AREA_EPS) continue;
        TripleKey k=sortedTri(f.a,f.b,f.c);
        if(seen.insert(k).second) nf.push_back(f);
    }
    mesh.F.swap(nf);
}

// ---------- Saliency and target choice ----------
struct SaliencyData {
    vector<int> order;
    vector<float> localScale;
    vector<unsigned char> boundary;
    vector<float> curvature;
    vector<Vec3> normal;
    double complexity=0.0;
    double boundaryFrac=0.0;
    double avgCurv=0.0;
    double smoothFrac=1.0;
};

static double sampleMedianEdge(const Mesh& mesh){
    vector<double> lens; lens.reserve(min<size_t>(mesh.F.size()*3, 240000));
    int step=max(1, (int)(mesh.F.size()/80000));
    int n=(int)mesh.V.size();
    for(int i=0;i<(int)mesh.F.size();i+=step){
        const Face& f=mesh.F[i]; if(!validFace(f,n)) continue;
        double a=normv(mesh.V[f.a]-mesh.V[f.b]);
        double b=normv(mesh.V[f.b]-mesh.V[f.c]);
        double c=normv(mesh.V[f.c]-mesh.V[f.a]);
        if(a>EPS) lens.push_back(a); if(b>EPS) lens.push_back(b); if(c>EPS) lens.push_back(c);
    }
    if(lens.empty()){
        Bounds bd=computeBounds(mesh.V);
        return max(1e-9, bd.diag / max(1.0, sqrt((double)max<size_t>(1,mesh.V.size()))));
    }
    nth_element(lens.begin(), lens.begin()+lens.size()/2, lens.end());
    return max(1e-12, lens[lens.size()/2]);
}

static SaliencyData computeSaliency(const Mesh& mesh){
    int n=(int)mesh.V.size();
    SaliencyData sd;
    sd.order.resize(n); iota(sd.order.begin(), sd.order.end(), 0);
    sd.localScale.assign(n, 1.0f);
    sd.boundary.assign(n, 0);
    sd.curvature.assign(n, 0.0f);
    if(n==0) return sd;

    vector<Vec3> normal(n, {0,0,0});
    for(const Face& f: mesh.F){
        if(!validFace(f,n)) continue;
        Vec3 cr=crossv(mesh.V[f.b]-mesh.V[f.a], mesh.V[f.c]-mesh.V[f.a]);
        if(norm2(cr)<=AREA_EPS) continue;
        normal[f.a]+=cr; normal[f.b]+=cr; normal[f.c]+=cr;
    }
    for(int i=0;i<n;++i) normal[i]=normalize(normal[i]);
    sd.normal = normal;

    // Boundary/non-manifold detection via sorted edge list.  Sorting is faster and
    // much more memory-predictable than a giant unordered_map for million meshes.
    vector<uint64_t> edges; edges.reserve(mesh.F.size()*3);
    auto packEdge=[&](int a,int b)->uint64_t{ if(a>b) swap(a,b); return ((uint64_t)(uint32_t)a<<32) | (uint32_t)b; };
    for(const Face& f: mesh.F){
        if(!validFace(f,n)) continue;
        edges.push_back(packEdge(f.a,f.b)); edges.push_back(packEdge(f.b,f.c)); edges.push_back(packEdge(f.c,f.a));
    }
    sort(edges.begin(), edges.end());
    for(size_t i=0;i<edges.size();){
        size_t j=i+1; while(j<edges.size() && edges[j]==edges[i]) ++j;
        if(j-i!=2){
            int a=(int)(edges[i]>>32), b=(int)(edges[i]&0xffffffffu);
            if(a>=0&&a<n) sd.boundary[a]=1; if(b>=0&&b<n) sd.boundary[b]=1;
        }
        i=j;
    }

    // Normal-variation saliency on the one-ring.  fabs(dot) avoids false creases
    // on inputs whose face orientation is inconsistent.
    for(const Face& f: mesh.F){
        if(!validFace(f,n)) continue;
        int v[3]={f.a,f.b,f.c};
        for(int e=0;e<3;++e){
            int a=v[e], b=v[(e+1)%3];
            double da=norm2(normal[a]), db=norm2(normal[b]);
            if(da>0.5 && db>0.5){
                double d=fabs(dotv(normal[a],normal[b]));
                d=min(1.0,max(0.0,d));
                float c=(float)(1.0-d);
                if(c>sd.curvature[a]) sd.curvature[a]=c;
                if(c>sd.curvature[b]) sd.curvature[b]=c;
            }
        }
    }

    int important=0, boundaryCnt=0, smoothCnt=0;
    double curvSum=0.0;
    vector<int> bucket(n,0);
    vector<uint64_t> rnd(n,0);
    Bounds bd=computeBounds(mesh.V);
    double sx=max(1e-30, bd.mx.x-bd.mn.x), sy=max(1e-30, bd.mx.y-bd.mn.y), sz=max(1e-30, bd.mx.z-bd.mn.z);
    for(int i=0;i<n;++i){
        float c=sd.curvature[i];
        bool bnd=sd.boundary[i]!=0;
        curvSum += c;
        if(bnd) ++boundaryCnt;
        if(!bnd && c<0.012f) ++smoothCnt;
        if(bnd || c>0.08f) ++important;

        // localScale is the radius multiplier used when deciding whether this
        // vertex is already covered.  Smooth regions deliberately use larger
        // values; important features use smaller values.
        float sc=1.18f;
        if(bnd) sc=0.30f;
        else if(c>0.55f) sc=0.36f;
        else if(c>0.25f) sc=0.46f;
        else if(c>0.10f) sc=0.58f;
        else if(c>0.035f) sc=0.76f;
        else if(c<0.004f) sc=1.48f;
        else if(c<0.012f) sc=1.34f;
        else sc=1.10f;
        sd.localScale[i]=sc;

        int qb=min(9000, (int)llround((double)c*6000.0));
        bucket[i]=(bnd?22000:0)+qb;
        // Coordinate-mixed randomization: improves spatial distribution in smooth regions.
        uint64_t qx=(uint64_t)llround((mesh.V[i].x-bd.mn.x)/sx*1048575.0);
        uint64_t qy=(uint64_t)llround((mesh.V[i].y-bd.mn.y)/sy*1048575.0);
        uint64_t qz=(uint64_t)llround((mesh.V[i].z-bd.mn.z)/sz*1048575.0);
        rnd[i]=splitmix64(((uint64_t)i+1)*0x9e3779b97f4a7c15ULL ^ (qx<<42) ^ (qy<<21) ^ qz);
    }
    sd.complexity = (double)important / max(1,n);
    sd.boundaryFrac = (double)boundaryCnt / max(1,n);
    sd.avgCurv = curvSum / max(1,n);
    sd.smoothFrac = (double)smoothCnt / max(1,n);
    sort(sd.order.begin(), sd.order.end(), [&](int a,int b){
        if(bucket[a]!=bucket[b]) return bucket[a]>bucket[b];
        return rnd[a]<rnd[b];
    });
    return sd;
}

static int chooseTargetVertices(int n, int m, const SaliencyData& sd){
    if(n<=4) return n;
    double feature=min(1.0, max(0.0, sd.complexity*5.5));
    double smooth=min(1.0, max(0.0, sd.smoothFrac));
    double ratio;
    if(n < 1000) ratio = 0.54 + 0.05*feature;
    else if(n < 5000) ratio = 0.300 + 0.050*feature;
    else if(n < 20000) ratio = 0.155 + 0.055*feature;
    else if(n < 80000) ratio = 0.088 + 0.050*feature;
    else if(n < 250000) ratio = 0.060 + 0.045*feature;
    else if(n < 700000) ratio = 0.049 + 0.041*feature;
    else ratio = 0.044 + 0.038*feature;

    // Smooth scan-like meshes are where the vertex-only clarification gives the
    // most leverage; preserve more on boundary/feature-heavy or sparse meshes.
    if(smooth > 0.82 && sd.boundaryFrac < 0.05) ratio *= 0.78;
    if(smooth > 0.94 && sd.avgCurv < 0.002) ratio *= 0.70;
    if(sd.boundaryFrac > 0.12) ratio *= 1.12;

    double faceRatio = (double)m / max(1,n);
    if(faceRatio > 2.6 && n > 50000) ratio *= 0.90;
    if(faceRatio < 1.25 && n > 5000) ratio *= 1.16;

    // Safety rails: avoid going into the region that usually creates empty or
    // visually meaningless contracted meshes, while still being significantly
    // more aggressive than the prior 81.x plateau.
    double minRatio = (n<20000 ? 0.035 : (n<250000 ? 0.018 : 0.012));
    ratio=max(ratio, minRatio);
    ratio=min(ratio, 0.62);

    int target=(int)llround(n*ratio);
    int floorKeep = (n<20000 ? 64 : 420);
    target=max(target, min(n, floorKeep));
    target=min(target, n);
    return target;
}

// ---------- Adaptive r-net ----------
struct CellKey { long long x,y,z; bool operator==(const CellKey& o) const noexcept { return x==o.x&&y==o.y&&z==o.z; } };
struct CellHash {
    size_t operator()(const CellKey& k) const noexcept {
        uint64_t x=(uint64_t)k.x, y=(uint64_t)k.y, z=(uint64_t)k.z;
        return (size_t)splitmix64(x*11995408973635179863ULL ^ y*10150724397891781847ULL ^ z*7809847782465536327ULL);
    }
};

struct RNetResult {
    vector<int> selectedOrig; // output index -> original vertex index
    vector<int> repOf;        // original vertex -> output index
    double maxAssignedDist=0.0;
};

class RNetBuilder {
public:
    const vector<Vec3>& V;
    const vector<int>& order;
    const vector<float>& localScale;
    Bounds bd;
    double maxScale=1.0;
    explicit RNetBuilder(const vector<Vec3>& VV, const vector<int>& ord, const vector<float>& sc)
        : V(VV), order(ord), localScale(sc), bd(computeBounds(VV)) {
        for(float x: localScale) if(x>maxScale) maxScale=x;
        if(maxScale < 1.0) maxScale = 1.0;
    }

    CellKey keyOf(const Vec3& p, double cell) const {
        return {(long long)floor((p.x-bd.mn.x)/cell),
                (long long)floor((p.y-bd.mn.y)/cell),
                (long long)floor((p.z-bd.mn.z)/cell)};
    }

    int countOnly(double R, int earlyLimit = INT_MAX) const {
        int n=(int)V.size();
        if(R <= EPS || bd.diag <= EPS) return n;
        double cell=R*maxScale;
        unordered_map<CellKey,int,CellHash> head;
        head.reserve((size_t)min(n, max(1024, earlyLimit))*2+1);
        vector<int> selected; selected.reserve(min(n, max(1024, earlyLimit)));
        vector<int> nxt; nxt.reserve(min(n, max(1024, earlyLimit)));
        for(int id: order){
            double rr = R * (double)localScale[id];
            double rr2 = rr*rr;
            CellKey ck=keyOf(V[id],cell);
            bool covered=false;
            // Because cell >= every local radius, checking adjacent cells is enough.
            for(long long dx=-1; dx<=1 && !covered; ++dx)
                for(long long dy=-1; dy<=1 && !covered; ++dy)
                    for(long long dz=-1; dz<=1 && !covered; ++dz){
                        auto it=head.find({ck.x+dx, ck.y+dy, ck.z+dz});
                        if(it==head.end()) continue;
                        for(int s=it->second; s!=-1; s=nxt[s]){
                            if(dist2(V[id], V[selected[s]]) <= rr2){ covered=true; break; }
                        }
                    }
            if(!covered){
                int sid=(int)selected.size(); selected.push_back(id);
                auto it=head.find(ck);
                int old = (it==head.end() ? -1 : it->second);
                nxt.push_back(old);
                if(it==head.end()) head.emplace(ck,sid); else it->second=sid;
                if((int)selected.size() > earlyLimit) return (int)selected.size();
            }
        }
        return (int)selected.size();
    }

    RNetResult build(double R) const {
        int n=(int)V.size();
        RNetResult res; res.repOf.assign(n,-1);
        if(R <= EPS || bd.diag <= EPS){
            res.selectedOrig.resize(n); iota(res.selectedOrig.begin(), res.selectedOrig.end(), 0);
            iota(res.repOf.begin(), res.repOf.end(), 0);
            return res;
        }
        double cell=R*maxScale;
        unordered_map<CellKey,int,CellHash> head;
        head.reserve((size_t)n/4+1024);
        vector<int> nxt; nxt.reserve(n/8+1024);

        for(int id: order){
            double rr=R*(double)localScale[id]; double rr2=rr*rr;
            CellKey ck=keyOf(V[id],cell);
            bool covered=false;
            for(long long dx=-1; dx<=1 && !covered; ++dx)
                for(long long dy=-1; dy<=1 && !covered; ++dy)
                    for(long long dz=-1; dz<=1 && !covered; ++dz){
                        auto it=head.find({ck.x+dx, ck.y+dy, ck.z+dz});
                        if(it==head.end()) continue;
                        for(int s=it->second; s!=-1; s=nxt[s]){
                            if(dist2(V[id], V[res.selectedOrig[s]]) <= rr2){ covered=true; break; }
                        }
                    }
            if(!covered){
                int sid=(int)res.selectedOrig.size(); res.selectedOrig.push_back(id);
                auto it=head.find(ck);
                int old=(it==head.end()?-1:it->second);
                nxt.push_back(old);
                if(it==head.end()) head.emplace(ck,sid); else it->second=sid;
            }
        }

        for(int i=0;i<n;++i){
            CellKey ck=keyOf(V[i],cell);
            int best=-1; double bestD=numeric_limits<double>::infinity();
            for(long long dx=-1; dx<=1; ++dx)
                for(long long dy=-1; dy<=1; ++dy)
                    for(long long dz=-1; dz<=1; ++dz){
                        auto it=head.find({ck.x+dx, ck.y+dy, ck.z+dz});
                        if(it==head.end()) continue;
                        for(int s=it->second; s!=-1; s=nxt[s]){
                            double d=dist2(V[i], V[res.selectedOrig[s]]);
                            if(d<bestD){ bestD=d; best=s; }
                        }
                    }
            if(best<0){
                // Extremely rare unless the radius is numerically tiny.
                for(int s=0;s<(int)res.selectedOrig.size();++s){
                    double d=dist2(V[i], V[res.selectedOrig[s]]);
                    if(d<bestD){ bestD=d; best=s; }
                }
            }
            res.repOf[i]=best;
            if(bestD>res.maxAssignedDist) res.maxAssignedDist=bestD;
        }
        res.maxAssignedDist=sqrt(max(0.0,res.maxAssignedDist));
        return res;
    }
};

static double chooseRadiusForTarget(const Mesh& mesh, const SaliencyData& sd, int target){
    int n=(int)mesh.V.size();
    if(target>=n) return 0.0;
    RNetBuilder rb(mesh.V, sd.order, sd.localScale);
    double med=sampleMedianEdge(mesh);
    Bounds bd=computeBounds(mesh.V);
    double low=0.0, high=max(med*1.5, bd.diag / max(16.0, sqrt((double)n)));
    if(high<=EPS) high=max(1e-9, bd.diag*1e-6+1e-9);

    vector<pair<double,int>> trials;
    auto eval=[&](double R)->int{
        int c=rb.countOnly(R, max(n, target*3));
        trials.push_back({R,c});
        return c;
    };

    int ch=eval(high);
    int guard=0;
    while(ch>target && high < max(1e-9, bd.diag*4.0) && guard++<30){
        low=high; high*=1.45; ch=eval(high);
    }
    if(ch>target){
        // Could not reach target without absurd radius; return the largest tried.
        return high;
    }
    for(int it=0; it<11; ++it){
        double mid=(low+high)*0.5;
        int c=eval(mid);
        if(c>target) low=mid; else high=mid;
    }

    // Pick a radius whose count is closest to target, but avoid falling far below it.
    double bestR=high; long long bestScore=LLONG_MAX; int bestC=INT_MAX;
    for(auto [R,c]: trials){
        long long diff=llabs((long long)c-target);
        long long penalty = diff*1000LL + (c<target ? (long long)(target-c)*35LL : 0LL);
        if(c < target*0.82) penalty += (long long)(target*0.82 - c)*5000LL;
        if(penalty<bestScore || (penalty==bestScore && c<bestC)){
            bestScore=penalty; bestR=R; bestC=c;
        }
    }
    return bestR;
}

// ---------- Output mesh construction ----------
static inline bool addLimitedNeighbor(vector<vector<int>>& adj, int a, int b){
    if(a==b || a<0 || b<0 || a>=(int)adj.size() || b>=(int)adj.size()) return false;
    auto& v=adj[a];
    for(int x: v) if(x==b) return false;
    if(v.size()<20) v.push_back(b);
    return true;
}

static void dedupFacesKeepVertices(Mesh& mesh){
    int n=(int)mesh.V.size();
    vector<Face> nf; nf.reserve(mesh.F.size());
    unordered_set<TripleKey,TripleHash> seen; seen.reserve(mesh.F.size()*2+1);
    for(const Face& f: mesh.F){
        if(!validFace(f,n)) continue;
        Vec3 cr=crossv(mesh.V[f.b]-mesh.V[f.a], mesh.V[f.c]-mesh.V[f.a]);
        if(norm2(cr)<=AREA_EPS) continue;
        TripleKey k=sortedTri(f.a,f.b,f.c);
        if(seen.insert(k).second) nf.push_back(f);
    }
    mesh.F.swap(nf);
}

static void repairNonManifoldEdgesKeepVertices(Mesh& mesh){
    int n=(int)mesh.V.size();
    for(int iter=0; iter<4; ++iter){
        dedupFacesKeepVertices(mesh);
        vector<double> area(mesh.F.size(),0.0);
        unordered_map<PairKey, vector<int>, PairHash> ef;
        ef.reserve(mesh.F.size()*4+1);
        for(int i=0;i<(int)mesh.F.size();++i){
            const Face& f=mesh.F[i]; if(!validFace(f,n)) continue;
            area[i]=normv(crossv(mesh.V[f.b]-mesh.V[f.a], mesh.V[f.c]-mesh.V[f.a]));
            int v[3]={f.a,f.b,f.c};
            for(int e=0;e<3;++e){ int a=v[e], b=v[(e+1)%3]; if(a>b) swap(a,b); ef[{a,b}].push_back(i); }
        }
        vector<unsigned char> kill(mesh.F.size(),0); bool changed=false;
        for(auto& kv: ef){
            auto& lst=kv.second; if(lst.size()<=2) continue;
            changed=true;
            sort(lst.begin(), lst.end(), [&](int i,int j){ return area[i]>area[j]; });
            for(size_t k=2;k<lst.size();++k) kill[lst[k]]=1;
        }
        if(!changed) break;
        vector<Face> nf; nf.reserve(mesh.F.size());
        for(int i=0;i<(int)mesh.F.size();++i) if(!kill[i]) nf.push_back(mesh.F[i]);
        mesh.F.swap(nf);
    }
    dedupFacesKeepVertices(mesh);
}

static void addSupportFaces(Mesh& out, const vector<vector<int>>& adj){
    int n=(int)out.V.size();
    if(n<3) return;
    vector<unsigned char> used(n,0);
    unordered_set<TripleKey,TripleHash> seen; seen.reserve(out.F.size()*2 + n + 1);
    unordered_map<PairKey,int,PairHash> edgeUse; edgeUse.reserve(out.F.size()*4 + n + 1);
    auto ekey=[](int a,int b){ if(a>b) swap(a,b); return PairKey{a,b}; };
    for(const Face& f: out.F){
        if(validFace(f,n)){
            used[f.a]=used[f.b]=used[f.c]=1; seen.insert(sortedTri(f.a,f.b,f.c));
            edgeUse[ekey(f.a,f.b)]++; edgeUse[ekey(f.b,f.c)]++; edgeUse[ekey(f.c,f.a)]++;
        }
    }
    auto tryAdd=[&](int a,int b,int c)->bool{
        if(a==b||b==c||c==a||a<0||b<0||c<0||a>=n||b>=n||c>=n) return false;
        Vec3 cr=crossv(out.V[b]-out.V[a], out.V[c]-out.V[a]);
        if(norm2(cr)<=AREA_EPS) return false;
        TripleKey k=sortedTri(a,b,c);
        if(seen.find(k)!=seen.end()) return false;
        PairKey e1=ekey(a,b), e2=ekey(b,c), e3=ekey(c,a);
        if(edgeUse[e1]>=2 || edgeUse[e2]>=2 || edgeUse[e3]>=2) return false;
        seen.insert(k); edgeUse[e1]++; edgeUse[e2]++; edgeUse[e3]++;
        out.F.emplace_back(a,b,c); used[a]=used[b]=used[c]=1; return true;
    };

    // Local supports from contracted one-ring adjacency.
    for(int i=0;i<n;++i){
        if(used[i]) continue;
        const auto& nb = (i<(int)adj.size()?adj[i]:vector<int>());
        bool ok=false;
        for(size_t a=0;a<nb.size() && !ok;++a) for(size_t b=a+1;b<nb.size() && !ok;++b){
            ok=tryAdd(i, nb[a], nb[b]);
        }
    }

    // Global emergency supports for vertices that remain isolated. This avoids
    // repeatedly using the same anchor-anchor edge, keeping edge valence <= 2.
    vector<int> anchors;
    anchors.reserve(12);
    auto addAnchor=[&](int id){ if(id>=0 && id<n && find(anchors.begin(),anchors.end(),id)==anchors.end()) anchors.push_back(id); };
    int minx=0,maxx=0,miny=0,maxy=0,minz=0,maxz=0;
    for(int i=1;i<n;++i){
        if(out.V[i].x<out.V[minx].x) minx=i; if(out.V[i].x>out.V[maxx].x) maxx=i;
        if(out.V[i].y<out.V[miny].y) miny=i; if(out.V[i].y>out.V[maxy].y) maxy=i;
        if(out.V[i].z<out.V[minz].z) minz=i; if(out.V[i].z>out.V[maxz].z) maxz=i;
    }
    addAnchor(minx); addAnchor(maxx); addAnchor(miny); addAnchor(maxy); addAnchor(minz); addAnchor(maxz);
    for(int i=0;i<n;++i){
        if(used[i]) continue;
        bool ok=false;
        uint64_t h=splitmix64((uint64_t)i+0x123456789abcdefULL);
        for(size_t off=0; off<anchors.size() && !ok; ++off){
            size_t a=(h+off)%anchors.size();
            for(size_t off2=1; off2<anchors.size() && !ok; ++off2){
                size_t b=(a+off2)%anchors.size();
                if(anchors[a]!=i && anchors[b]!=i) ok=tryAdd(i,anchors[a],anchors[b]);
            }
        }
    }
}

static Mesh buildContractedMesh(const Mesh& in, const RNetResult& rr, const SaliencyData& sd, bool relocateSmooth){
    Mesh out; out.fmt=in.fmt; out.rawBase=in.rawBase;
    int k=(int)rr.selectedOrig.size();
    out.V.reserve(k);

    vector<Vec3> sum(k, {0,0,0});
    vector<int> cnt(k,0);
    vector<unsigned char> hasBoundary(k,0);
    vector<float> maxCurv(k,0.0f);
    for(int i=0;i<(int)in.V.size();++i){
        int r = (i<(int)rr.repOf.size() ? rr.repOf[i] : -1);
        if(r<0 || r>=k) continue;
        sum[r] += in.V[i];
        cnt[r]++;
        if(i<(int)sd.boundary.size() && sd.boundary[i]) hasBoundary[r]=1;
        if(i<(int)sd.curvature.size() && sd.curvature[i] > maxCurv[r]) maxCurv[r]=sd.curvature[i];
    }

    for(int oi=0; oi<k; ++oi){
        int orig=rr.selectedOrig[oi];
        Vec3 p=in.V[orig];
        if(relocateSmooth && cnt[oi] >= 3 && !hasBoundary[oi] && maxCurv[oi] < 0.060f){
            Vec3 cen=sum[oi] / (double)max(1,cnt[oi]);
            Vec3 d=cen-p;
            if(orig>=0 && orig<(int)sd.normal.size() && norm2(sd.normal[orig])>0.5){
                Vec3 nrm=sd.normal[orig];
                d = d - nrm * (0.90 * dotv(d,nrm)); // tangent-plane biased recentering
            }
            double alpha = (maxCurv[oi] < 0.006f ? 0.82 : 0.55);
            p = p + d*alpha;
        }
        out.V.push_back(p);
    }

    out.F.reserve(min<size_t>(in.F.size(), (size_t)k*8+1024));
    vector<vector<int>> adj(k);
    vector<unsigned char> used(k,0);
    unordered_set<TripleKey,TripleHash> seen; seen.reserve(min<size_t>(in.F.size()*2+1, (size_t)k*20+1024));
    int nIn=(int)in.V.size();
    for(const Face& f: in.F){
        if(!validFace(f,nIn)) continue;
        int a=rr.repOf[f.a], b=rr.repOf[f.b], c=rr.repOf[f.c];
        if(a<0||b<0||c<0) continue;
        addLimitedNeighbor(adj,a,b); addLimitedNeighbor(adj,b,a);
        addLimitedNeighbor(adj,b,c); addLimitedNeighbor(adj,c,b);
        addLimitedNeighbor(adj,c,a); addLimitedNeighbor(adj,a,c);
        if(a==b || b==c || c==a) continue;
        Vec3 cr=crossv(out.V[b]-out.V[a], out.V[c]-out.V[a]);
        if(norm2(cr)<=AREA_EPS) continue;
        TripleKey tk=sortedTri(a,b,c);
        if(seen.insert(tk).second){ out.F.emplace_back(a,b,c); used[a]=used[b]=used[c]=1; }
    }
    if(out.F.empty() && k>=3){
        int a=0,b=0; double best=0;
        int step=max(1,k/2000);
        for(int i=0;i<k;i+=step) for(int j=i+step;j<k;j+=step){ double d=dist2(out.V[i],out.V[j]); if(d>best){best=d; a=i; b=j;} }
        int c=-1; best=0;
        for(int i=0;i<k;i+=step){ double ar=norm2(crossv(out.V[b]-out.V[a], out.V[i]-out.V[a])); if(i!=a&&i!=b&&ar>best){best=ar;c=i;} }
        if(c<0){ for(int i=0;i<k;++i) if(i!=a&&i!=b){c=i;break;} }
        if(c>=0) out.F.emplace_back(a,b,c);
    }

    repairNonManifoldEdgesKeepVertices(out);
    addSupportFaces(out, adj);
    dedupFacesKeepVertices(out);
    return out;
}

static Mesh simplifyByAdaptiveRNet(Mesh mesh){
    cleanFacesKeepVertices(mesh);
    int n=(int)mesh.V.size();
    if(n<=4 || mesh.F.empty()) return mesh;
    SaliencyData sd=computeSaliency(mesh);
    int target=chooseTargetVertices(n, (int)mesh.F.size(), sd);
    target=max(3, min(target,n));
    if(target>=n) return mesh;

    double R=chooseRadiusForTarget(mesh, sd, target);
    RNetBuilder rb(mesh.V, sd.order, sd.localScale);
    RNetResult rr=rb.build(R);

    // If the binary search landed unexpectedly far from target because of pathological
    // coordinates, retry once with a small radius correction.
    if((int)rr.selectedOrig.size() > target*13/10 && R>EPS){
        rr=rb.build(R*1.10);
    } else if((int)rr.selectedOrig.size() < target*72/100 && R>EPS){
        rr=rb.build(R*0.90);
    }

    bool relocateSmooth = (n >= 2500);
    Mesh out=buildContractedMesh(mesh, rr, sd, relocateSmooth);
    // Structural fallback: if aggressive smooth redistribution produced too few
    // faces, rebuild without relocating representatives.
    if((out.F.empty() || out.F.size() < max<size_t>(1, out.V.size()/8)) && rr.selectedOrig.size() >= 3){
        Mesh plain = buildContractedMesh(mesh, rr, sd, false);
        if(!plain.F.empty() && plain.F.size() > out.F.size()) out = std::move(plain);
    }

    // Never return an empty or larger mesh.  The vertex set is intentionally not
    // compacted: all selected representatives participate in vertex-Hausdorff cover.
    if(out.V.size()<3 || out.F.empty() || out.V.size()>=mesh.V.size()) return mesh;
    return out;
}

static void writeMesh(const Mesh& mesh, ostream& out){
    out.setf(ios::fixed); out << setprecision(10);
    if(mesh.fmt==Format::OBJ){
        for(const Vec3& p: mesh.V) out << "v " << p.x << ' ' << p.y << ' ' << p.z << '\n';
        for(const Face& f: mesh.F) out << "f " << f.a+1 << ' ' << f.b+1 << ' ' << f.c+1 << '\n';
    }else if(mesh.fmt==Format::OFF){
        out << "OFF\n" << mesh.V.size() << ' ' << mesh.F.size() << " 0\n";
        for(const Vec3& p: mesh.V) out << p.x << ' ' << p.y << ' ' << p.z << '\n';
        for(const Face& f: mesh.F) out << "3 " << f.a << ' ' << f.b << ' ' << f.c << '\n';
    }else{
        int base=mesh.rawBase;
        out << mesh.V.size() << ' ' << mesh.F.size() << '\n';
        for(const Vec3& p: mesh.V) out << p.x << ' ' << p.y << ' ' << p.z << '\n';
        for(const Face& f: mesh.F) out << f.a+base << ' ' << f.b+base << ' ' << f.c+base << '\n';
    }
}

int main(){
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    string data((istreambuf_iterator<char>(cin)), istreambuf_iterator<char>());
    Mesh mesh;
    if(!parseMesh(data, mesh)){
        cout << data;
        return 0;
    }
    Mesh out=simplifyByAdaptiveRNet(std::move(mesh));
    writeMesh(out, cout);
    return 0;
}
