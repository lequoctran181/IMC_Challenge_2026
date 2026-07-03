#include <bits/stdc++.h>
using namespace std;

struct Vec3 {
    double x, y, z;
};
static inline Vec3 make_vec(double x=0,double y=0,double z=0){return {x,y,z};}
static inline Vec3 operator+(const Vec3&a,const Vec3&b){return {a.x+b.x,a.y+b.y,a.z+b.z};}
static inline Vec3 operator-(const Vec3&a,const Vec3&b){return {a.x-b.x,a.y-b.y,a.z-b.z};}
static inline Vec3 operator*(const Vec3&a,double s){return {a.x*s,a.y*s,a.z*s};}
static inline Vec3 operator/(const Vec3&a,double s){return {a.x/s,a.y/s,a.z/s};}
static inline double dot3(const Vec3&a,const Vec3&b){return a.x*b.x+a.y*b.y+a.z*b.z;}
static inline Vec3 cross3(const Vec3&a,const Vec3&b){return {a.y*b.z-a.z*b.y, a.z*b.x-a.x*b.z, a.x*b.y-a.y*b.x};}
static inline double norm2(const Vec3&a){return dot3(a,a);} 
static inline double norm3(const Vec3&a){return sqrt(norm2(a));}
static inline Vec3 normalized(const Vec3&a){ double n=norm3(a); return n>0? a*(1.0/n):Vec3{0,0,0}; }

struct Face { int v[3]; };

struct FastInput {
    vector<char> buf;
    char* p;
    FastInput(){
        buf.reserve(1<<27);
        char tmp[1<<16];
        size_t n;
        while((n=fread(tmp,1,sizeof(tmp),stdin))>0) buf.insert(buf.end(),tmp,tmp+n);
        buf.push_back('\0'); p=buf.data();
    }
    inline void skip(){ while(*p==' '||*p=='\n'||*p=='\r'||*p=='\t') ++p; }
    inline long nextLong(){ skip(); return strtol(p,&p,10); }
    inline double nextDouble(){ skip(); return strtod(p,&p); }
    inline char nextChar(){ skip(); return *p++; }
};

struct Quadric {
    // Symmetric 4x4 matrix, upper triangle: xx xy xz xw yy yz yw zz zw ww
    double a00=0,a01=0,a02=0,a03=0,a11=0,a12=0,a13=0,a22=0,a23=0,a33=0;
    inline void clear(){ a00=a01=a02=a03=a11=a12=a13=a22=a23=a33=0; }
    inline void addPlane(double x,double y,double z,double w,double weight=1.0){
        x*=weight; y*=weight; z*=weight; w*=weight;
        // Use sqrt(weight) outside? Here weight is already multiplicative for plane vector.
        a00 += x*x; a01 += x*y; a02 += x*z; a03 += x*w;
        a11 += y*y; a12 += y*z; a13 += y*w;
        a22 += z*z; a23 += z*w;
        a33 += w*w;
    }
    inline void addWeightedPlane(double x,double y,double z,double w,double weight){
        if(weight<=0) return;
        double s=sqrt(weight);
        addPlane(x*s,y*s,z*s,w*s,1.0);
    }
    inline Quadric& operator+=(const Quadric&o){
        a00+=o.a00; a01+=o.a01; a02+=o.a02; a03+=o.a03;
        a11+=o.a11; a12+=o.a12; a13+=o.a13; a22+=o.a22; a23+=o.a23; a33+=o.a33;
        return *this;
    }
    inline double eval(const Vec3&p) const {
        double x=p.x,y=p.y,z=p.z;
        double v = 0;
        v += a00*x*x + 2*a01*x*y + 2*a02*x*z + 2*a03*x;
        v += a11*y*y + 2*a12*y*z + 2*a13*y;
        v += a22*z*z + 2*a23*z;
        v += a33;
        return v;
    }
};

static int N=0,M=0;
static vector<Vec3> P, origP;
static vector<Face> F, origF;
static vector<unsigned char> aliveV, aliveF;
static vector<vector<int>> inc;
static vector<Quadric> Q;
static vector<double> radius_cluster;
static vector<int> cluster_size;
static int alive_vertices=0, alive_faces=0;
static double mesh_diag=1.0, hausdorff_R=0.05;
static chrono::steady_clock::time_point t0;
static vector<int> markA, markB;
static int markTag=1, seenTag=1;

static inline double elapsed(){ return chrono::duration<double>(chrono::steady_clock::now()-t0).count(); }
static inline unsigned long long edge_key(int a,int b){ if(a>b) swap(a,b); return (unsigned long long)(unsigned int)a<<32 | (unsigned int)b; }
static inline bool face_has(const Face&f,int v){ return f.v[0]==v || f.v[1]==v || f.v[2]==v; }
static inline bool face_has_id(int fid,int v){ const Face&f=F[fid]; return f.v[0]==v||f.v[1]==v||f.v[2]==v; }
static inline int third_vertex(const Face&f,int a,int b){ for(int i=0;i<3;i++){int x=f.v[i]; if(x!=a && x!=b) return x;} return -1; }
static inline array<int,3> sorted_tri(int a,int b,int c){ array<int,3>s={a,b,c}; sort(s.begin(),s.end()); return s; }
static inline bool same_unordered(const Face&f,int a,int b,int c){ return sorted_tri(f.v[0],f.v[1],f.v[2])==sorted_tri(a,b,c); }

static inline Vec3 face_cross_current(const Face&f){ return cross3(P[f.v[1]]-P[f.v[0]], P[f.v[2]]-P[f.v[0]]); }
static inline Vec3 face_cross_with_pos(const Face&f,int replaced,int keep,const Vec3&newpos){
    Vec3 a = (f.v[0]==replaced || f.v[0]==keep) ? newpos : P[f.v[0]];
    Vec3 b = (f.v[1]==replaced || f.v[1]==keep) ? newpos : P[f.v[1]];
    Vec3 c = (f.v[2]==replaced || f.v[2]==keep) ? newpos : P[f.v[2]];
    return cross3(b-a,c-a);
}

static void load_mesh(){
    FastInput in;
    N=(int)in.nextLong(); M=(int)in.nextLong();
    P.resize(N); origP.resize(N);
    aliveV.assign(N,1); aliveF.assign(M,1);
    Q.assign(N, Quadric()); radius_cluster.assign(N,0.0); cluster_size.assign(N,1);
    Vec3 mn{1e100,1e100,1e100}, mx{-1e100,-1e100,-1e100};
    for(int i=0;i<N;i++){
        (void)in.nextChar();
        P[i].x=in.nextDouble(); P[i].y=in.nextDouble(); P[i].z=in.nextDouble();
        origP[i]=P[i];
        mn.x=min(mn.x,P[i].x); mn.y=min(mn.y,P[i].y); mn.z=min(mn.z,P[i].z);
        mx.x=max(mx.x,P[i].x); mx.y=max(mx.y,P[i].y); mx.z=max(mx.z,P[i].z);
    }
    mesh_diag = norm3(mx-mn); if(!(mesh_diag>0)) mesh_diag=1.0;
    hausdorff_R = 0.05 * mesh_diag * 0.995; // tiny safety margin for printing/rounding.
    F.resize(M); origF.resize(M);
    vector<int> deg(N,0);
    for(int i=0;i<M;i++){
        (void)in.nextChar();
        int a=(int)in.nextLong()-1,b=(int)in.nextLong()-1,c=(int)in.nextLong()-1;
        F[i].v[0]=a; F[i].v[1]=b; F[i].v[2]=c; origF[i]=F[i];
        deg[a]++;deg[b]++;deg[c]++;
    }
    inc.resize(N);
    for(int i=0;i<N;i++) inc[i].reserve(deg[i]+8);
    for(int i=0;i<M;i++){ inc[F[i].v[0]].push_back(i); inc[F[i].v[1]].push_back(i); inc[F[i].v[2]].push_back(i); }
    alive_vertices=N; alive_faces=M;
    markA.assign(N,0); markB.assign(N,0);
}

struct EdgeFace { unsigned long long key; int fid; };
static void build_quadrics_and_initial_edges(vector<unsigned long long>& unique_edges){
    vector<EdgeFace> ef;
    ef.reserve((size_t)M*3);
    // Face plane quadrics. Area weighting is intentionally mild: it keeps smooth surfaces smooth,
    // but still gives tiny high-curvature triangles enough influence through the unweighted floor.
    for(int i=0;i<M;i++){
        const Face&f=F[i];
        Vec3 cr=face_cross_current(f);
        double len=norm3(cr);
        if(!(len>0)) continue;
        Vec3 n=cr*(1.0/len);
        double d=-dot3(n,P[f.v[0]]);
        double area=0.5*len;
        double w=0.35 + sqrt(max(0.0, area));
        for(int k=0;k<3;k++) Q[f.v[k]].addWeightedPlane(n.x,n.y,n.z,d,w);
        ef.push_back({edge_key(f.v[0],f.v[1]),i});
        ef.push_back({edge_key(f.v[1],f.v[2]),i});
        ef.push_back({edge_key(f.v[2],f.v[0]),i});
    }
    sort(ef.begin(), ef.end(), [](const EdgeFace&a,const EdgeFace&b){return a.key<b.key;});
    unique_edges.clear(); unique_edges.reserve(ef.size()/2);
    const double sharp_cos = cos(28.0 * acos(-1.0) / 180.0);
    for(size_t i=0;i<ef.size();){
        size_t j=i+1; while(j<ef.size() && ef[j].key==ef[i].key) ++j;
        unique_edges.push_back(ef[i].key);
        if(j-i==2){
            int f0=ef[i].fid, f1=ef[i+1].fid;
            Vec3 c0=face_cross_current(F[f0]), c1=face_cross_current(F[f1]);
            double l0=norm3(c0), l1=norm3(c1);
            if(l0>0 && l1>0){
                Vec3 n0=c0*(1.0/l0), n1=c1*(1.0/l1);
                double dc=dot3(n0,n1); dc=max(-1.0,min(1.0,dc));
                if(dc < sharp_cos){
                    int a=(int)(ef[i].key>>32), b=(int)(ef[i].key & 0xffffffffu);
                    Vec3 e=normalized(P[b]-P[a]);
                    Vec3 bp0=normalized(cross3(e,n0));
                    Vec3 bp1=normalized(cross3(e,n1));
                    double w = (dc < cos(60.0*acos(-1.0)/180.0)) ? 90.0 : 35.0;
                    if(norm2(bp0)>0){ double d=-dot3(bp0,P[a]); Q[a].addWeightedPlane(bp0.x,bp0.y,bp0.z,d,w); Q[b].addWeightedPlane(bp0.x,bp0.y,bp0.z,d,w); }
                    if(norm2(bp1)>0){ double d=-dot3(bp1,P[a]); Q[a].addWeightedPlane(bp1.x,bp1.y,bp1.z,d,w); Q[b].addWeightedPlane(bp1.x,bp1.y,bp1.z,d,w); }
                }
            }
        }
        i=j;
    }
}

static bool collect_shared(int u,int v,int shared[2],int opp[2]){
    int cnt=0;
    const vector<int>& lu = inc[u].size() <= inc[v].size()? inc[u] : inc[v];
    for(int fid: lu){
        if(!aliveF[fid]) continue;
        if(!face_has_id(fid,u) || !face_has_id(fid,v)) continue;
        if(cnt>=2) return false;
        shared[cnt]=fid; opp[cnt]=third_vertex(F[fid],u,v); cnt++;
    }
    return cnt==2 && opp[0]>=0 && opp[1]>=0 && opp[0]!=opp[1];
}

static bool link_ok(int u,int v,const int opp[2]){
    if(++markTag > 2000000000){ fill(markA.begin(),markA.end(),0); markTag=1; }
    if(++seenTag > 2000000000){ fill(markB.begin(),markB.end(),0); seenTag=1; }
    for(int fid: inc[u]) if(aliveF[fid] && face_has_id(fid,u)){
        const Face&f=F[fid];
        for(int k=0;k<3;k++){ int x=f.v[k]; if(x!=u && x!=v) markA[x]=markTag; }
    }
    int inter=0;
    for(int fid: inc[v]) if(aliveF[fid] && face_has_id(fid,v)){
        const Face&f=F[fid];
        for(int k=0;k<3;k++){
            int x=f.v[k]; if(x==u||x==v) continue;
            if(markA[x]!=markTag) continue;
            if(x!=opp[0] && x!=opp[1]) return false;
            if(markB[x]!=seenTag){ markB[x]=seenTag; inter++; }
        }
    }
    return inter==2 && markB[opp[0]]==seenTag && markB[opp[1]]==seenTag;
}

static bool solve_optimal(const Quadric&q, Vec3&out){
    // Solve A x = -b for the 3x3 leading block.
    double a00=q.a00, a01=q.a01, a02=q.a02;
    double a10=q.a01, a11=q.a11, a12=q.a12;
    double a20=q.a02, a21=q.a12, a22=q.a22;
    double b0=-q.a03,b1=-q.a13,b2=-q.a23;
    double det = a00*(a11*a22-a12*a21) - a01*(a10*a22-a12*a20) + a02*(a10*a21-a11*a20);
    double scale = fabs(a00)+fabs(a11)+fabs(a22)+fabs(a01)+fabs(a02)+fabs(a12)+1e-30;
    if(fabs(det) < 1e-12*scale*scale*scale) return false;
    double dx = b0*(a11*a22-a12*a21) - a01*(b1*a22-a12*b2) + a02*(b1*a21-a11*b2);
    double dy = a00*(b1*a22-a12*b2) - b0*(a10*a22-a12*a20) + a02*(a10*b2-b1*a20);
    double dz = a00*(a11*b2-b1*a21) - a01*(a10*b2-b1*a20) + b0*(a10*a21-a11*a20);
    out={dx/det,dy/det,dz/det};
    return isfinite(out.x)&&isfinite(out.y)&&isfinite(out.z);
}

struct CandidateEval { bool ok=false; double cost=1e300; Vec3 pos{}; double newrad=0; };

static inline double cluster_rad_for_pos(int u,int v,const Vec3&p){
    return max(radius_cluster[u] + norm3(p-P[u]), radius_cluster[v] + norm3(p-P[v]));
}

static bool would_duplicate(int keep,int rem,int fid,int a,int b,int c,int s0,int s1){
    int probe=a;
    if(inc[b].size()<inc[probe].size()) probe=b;
    if(inc[c].size()<inc[probe].size()) probe=c;
    if(inc[keep].size()<inc[probe].size()) probe=keep;
    for(int ofid: inc[probe]){
        if(!aliveF[ofid] || ofid==fid || ofid==s0 || ofid==s1) continue;
        if(face_has_id(ofid,rem)) continue;
        if(same_unordered(F[ofid],a,b,c)) return true;
    }
    return false;
}

static CandidateEval quick_candidate(int u,int v){
    CandidateEval best;
    if(!aliveV[u]||!aliveV[v]||u==v) return best;
    Quadric q=Q[u]; q += Q[v];
    vector<Vec3> cand; cand.reserve(7);
    Vec3 opt;
    if(solve_optimal(q,opt)) cand.push_back(opt);
    cand.push_back(P[u]); cand.push_back(P[v]);
    cand.push_back((P[u]+P[v])*0.5);
    double su=max(1,cluster_size[u]), sv=max(1,cluster_size[v]);
    cand.push_back((P[u]*su + P[v]*sv)/(su+sv));
    if(cand.size()>0){
        Vec3 e=P[v]-P[u]; double l2=norm2(e);
        if(l2>1e-30){
            Vec3 base = cand[0];
            double t=dot3(base-P[u],e)/l2; t=max(0.0,min(1.0,t));
            cand.push_back(P[u]+e*t);
        }
    }
    double elen2=norm2(P[u]-P[v]);
    for(const Vec3&p: cand){
        if(!isfinite(p.x)||!isfinite(p.y)||!isfinite(p.z)) continue;
        double nr=cluster_rad_for_pos(u,v,p);
        if(nr > hausdorff_R) continue;
        double c=q.eval(p); if(c<0 && c>-1e-12) c=0;
        if(!isfinite(c)) continue;
        double rc = hausdorff_R>0 ? nr/hausdorff_R : 0.0;
        c += 1e-10*elen2 + 1e-5*rc*rc;
        if(c<best.cost){ best.ok=true; best.cost=c; best.pos=p; best.newrad=nr; }
    }
    return best;
}

static CandidateEval full_evaluate_directed(int keep,int rem,const int shared[2], bool allow_relocate){
    CandidateEval base = quick_candidate(keep,rem);
    if(!base.ok) return base;
    // Directed collapse is represented by the same candidate set; no need to keep endpoint only.
    Quadric q=Q[keep]; q += Q[rem];
    vector<Vec3> cand; cand.reserve(10);
    Vec3 opt;
    if(solve_optimal(q,opt)) cand.push_back(opt);
    cand.push_back(P[keep]); cand.push_back(P[rem]);
    cand.push_back((P[keep]+P[rem])*0.5);
    double sk=max(1,cluster_size[keep]), sr=max(1,cluster_size[rem]);
    cand.push_back((P[keep]*sk + P[rem]*sr)/(sk+sr));
    if(!allow_relocate){
        // Still keep midpoint/weighted candidates when they are very close: this improves planar cleanup
        // without materially increasing Hausdorff risk.
    }
    if(!cand.empty()){
        Vec3 e=P[rem]-P[keep]; double l2=norm2(e);
        if(l2>1e-30){
            Vec3 basep=cand[0];
            double t=dot3(basep-P[keep],e)/l2; t=max(0.0,min(1.0,t));
            cand.push_back(P[keep]+e*t);
        }
    }
    CandidateEval best;
    const double area_eps = max(1e-30, mesh_diag*mesh_diag*1e-24);
    const double min_dot = -0.05; // forbid real flips, allow gradual smoothing.
    for(const Vec3&pos: cand){
        if(!isfinite(pos.x)||!isfinite(pos.y)||!isfinite(pos.z)) continue;
        double nr=cluster_rad_for_pos(keep,rem,pos);
        if(nr > hausdorff_R) continue;
        bool ok=true;
        double normal_penalty=0.0;
        int changed=0;
        // Check all faces incident to keep or rem that will survive.
        static vector<int> affected;
        affected.clear();
        affected.reserve(inc[keep].size()+inc[rem].size());
        for(int fid: inc[keep]) if(aliveF[fid] && fid!=shared[0] && fid!=shared[1] && face_has_id(fid,keep)) affected.push_back(fid);
        for(int fid: inc[rem]) if(aliveF[fid] && fid!=shared[0] && fid!=shared[1] && face_has_id(fid,rem)) affected.push_back(fid);
        sort(affected.begin(),affected.end()); affected.erase(unique(affected.begin(),affected.end()),affected.end());
        for(int fid: affected){
            Face f=F[fid];
            int nv[3]={f.v[0],f.v[1],f.v[2]};
            for(int k=0;k<3;k++) if(nv[k]==rem) nv[k]=keep;
            if(nv[0]==nv[1]||nv[0]==nv[2]||nv[1]==nv[2]) { ok=false; break; }
            Vec3 oldc=face_cross_current(F[fid]);
            double oldl=norm3(oldc);
            Vec3 pa = (nv[0]==keep?pos:P[nv[0]]);
            Vec3 pb = (nv[1]==keep?pos:P[nv[1]]);
            Vec3 pc = (nv[2]==keep?pos:P[nv[2]]);
            Vec3 newc=cross3(pb-pa,pc-pa);
            double newl=norm3(newc);
            if(!(oldl>area_eps) || !(newl>area_eps)) { ok=false; break; }
            double nd=dot3(oldc,newc)/(oldl*newl); nd=max(-1.0,min(1.0,nd));
            if(nd < min_dot) { ok=false; break; }
            if(newl < oldl*1e-8) { ok=false; break; }
            if(face_has_id(fid,rem)){
                if(would_duplicate(keep,rem,fid,nv[0],nv[1],nv[2],shared[0],shared[1])) { ok=false; break; }
            }
            normal_penalty += max(0.0, 1.0-nd);
            changed++;
        }
        if(!ok || changed==0) continue;
        double c=q.eval(pos); if(c<0 && c>-1e-12) c=0;
        if(!isfinite(c)) continue;
        double rc=hausdorff_R>0?nr/hausdorff_R:0.0;
        double el=norm2(P[keep]-P[rem]);
        c += 0.012*normal_penalty + 2e-6*rc*rc + 1e-12*el + 1e-7*changed;
        if(c<best.cost){ best.ok=true; best.cost=c; best.pos=pos; best.newrad=nr; }
    }
    return best;
}

static CandidateEval full_evaluate_edge(int u,int v,int&keep,int&rem,int shared_out[2]){
    CandidateEval none;
    keep=rem=-1; shared_out[0]=shared_out[1]=-1;
    if(u==v||!aliveV[u]||!aliveV[v]) return none;
    int shared[2]={-1,-1}, opp[2]={-1,-1};
    if(!collect_shared(u,v,shared,opp)) return none;
    if(!link_ok(u,v,opp)) return none;
    CandidateEval uv=full_evaluate_directed(u,v,shared,true);
    CandidateEval vu=full_evaluate_directed(v,u,shared,true);
    if(!uv.ok && !vu.ok) return none;
    if(vu.ok && (!uv.ok || vu.cost<uv.cost)){ keep=v; rem=u; shared_out[0]=shared[0]; shared_out[1]=shared[1]; return vu; }
    keep=u; rem=v; shared_out[0]=shared[0]; shared_out[1]=shared[1]; return uv;
}

struct PQEdge { double cost; int u,v; };
struct PQCmp { bool operator()(const PQEdge&a,const PQEdge&b) const { return a.cost > b.cost; } };
static priority_queue<PQEdge, vector<PQEdge>, PQCmp> pq;
static bool pq_push_enabled = true;

static void push_edge_quick(int u,int v){
    if(u==v||u<0||v<0||u>=N||v>=N||!aliveV[u]||!aliveV[v]) return;
    CandidateEval e=quick_candidate(u,v);
    if(e.ok && isfinite(e.cost)) pq.push({e.cost,u,v});
}

static void rebuild_incident_keep(int keep,int rem){
    vector<int> merged;
    merged.reserve(inc[keep].size()+inc[rem].size()+16);
    for(int fid: inc[keep]) if(aliveF[fid] && face_has_id(fid,keep)) merged.push_back(fid);
    for(int fid: inc[rem]) if(aliveF[fid] && face_has_id(fid,keep)) merged.push_back(fid);
    sort(merged.begin(), merged.end()); merged.erase(unique(merged.begin(),merged.end()),merged.end());
    inc[keep].swap(merged);
    vector<int>().swap(inc[rem]);
}

static void apply_collapse(int keep,int rem,const int shared[2],const Vec3&newpos,double newrad){
    for(int i=0;i<2;i++) if(shared[i]>=0 && aliveF[shared[i]]){ aliveF[shared[i]]=0; alive_faces--; }
    for(int fid: inc[rem]){
        if(!aliveF[fid] || !face_has_id(fid,rem)) continue;
        Face&f=F[fid];
        for(int k=0;k<3;k++) if(f.v[k]==rem) f.v[k]=keep;
    }
    aliveV[rem]=0; alive_vertices--;
    P[keep]=newpos;
    radius_cluster[keep]=newrad;
    cluster_size[keep]+=cluster_size[rem];
    Q[keep]+=Q[rem];
    rebuild_incident_keep(keep,rem);
    // Queue all edges touching the one-ring of keep when the heap mode is active.
    if(pq_push_enabled){
        for(int fid: inc[keep]) if(aliveF[fid]){
            const Face&f=F[fid];
            push_edge_quick(f.v[0],f.v[1]);
            push_edge_quick(f.v[1],f.v[2]);
            push_edge_quick(f.v[2],f.v[0]);
        }
    }
}


struct LenEdge { double len2; unsigned long long key; };
static int sorted_edge_pass(double time_limit, int target_vertices){
    vector<unsigned long long> keys;
    keys.reserve((size_t)alive_faces*3);
    for(int fid=0; fid<(int)F.size(); ++fid){
        if(!aliveF[fid]) continue;
        const Face&f=F[fid];
        if(!aliveV[f.v[0]]||!aliveV[f.v[1]]||!aliveV[f.v[2]]) continue;
        keys.push_back(edge_key(f.v[0],f.v[1]));
        keys.push_back(edge_key(f.v[1],f.v[2]));
        keys.push_back(edge_key(f.v[2],f.v[0]));
    }
    sort(keys.begin(), keys.end());
    keys.erase(unique(keys.begin(), keys.end()), keys.end());
    vector<LenEdge> edges;
    edges.reserve(keys.size());
    for(unsigned long long key: keys){
        int a=(int)(key>>32), b=(int)(key & 0xffffffffu);
        if(a>=0&&a<N&&b>=0&&b<N&&aliveV[a]&&aliveV[b]) edges.push_back({norm2(P[a]-P[b]),key});
    }
    sort(edges.begin(), edges.end(), [](const LenEdge&a,const LenEdge&b){ return a.len2 < b.len2; });
    int collapsed=0, checked=0;
    for(const LenEdge&e: edges){
        if(alive_vertices<=target_vertices) break;
        if((++checked & 2047)==0 && elapsed()>time_limit) break;
        int a=(int)(e.key>>32), b=(int)(e.key & 0xffffffffu);
        if(a<0||a>=N||b<0||b>=N||!aliveV[a]||!aliveV[b]||a==b) continue;
        int keep,rem,shared[2];
        CandidateEval cur=full_evaluate_edge(a,b,keep,rem,shared);
        if(!cur.ok) continue;
        apply_collapse(keep,rem,shared,cur.pos,cur.newrad);
        ++collapsed;
    }
    return collapsed;
}
static void fast_decimate_to(int target_vertices, double time_limit){
    target_vertices=max(4,target_vertices);
    pq_push_enabled=false;
    int pass=0;
    while(alive_vertices>target_vertices && elapsed()<time_limit){
        int before=alive_vertices;
        int c=sorted_edge_pass(time_limit,target_vertices);
        ++pass;
        if(c==0 || alive_vertices>=before) break;
        if(pass>=12) break;
    }
    pq_push_enabled=true;
}

static void decimate_to(int target_vertices, double time_limit){
    target_vertices=max(4,target_vertices);
    int fail_streak=0;
    while(alive_vertices>target_vertices && !pq.empty()){
        if((fail_streak & 2047)==0 && elapsed()>time_limit) break;
        PQEdge e=pq.top(); pq.pop();
        if(e.u<0||e.v<0||e.u>=N||e.v>=N||!aliveV[e.u]||!aliveV[e.v]||e.u==e.v){ fail_streak++; continue; }
        // Re-evaluate. If the edge is still too expensive relative to the stale key, reinsert once.
        int keep,rem,shared[2];
        CandidateEval cur=full_evaluate_edge(e.u,e.v,keep,rem,shared);
        if(!cur.ok){ fail_streak++; continue; }
        if(cur.cost > e.cost*1.25 + 1e-10){
            pq.push({cur.cost,e.u,e.v});
            fail_streak++;
            if(fail_streak>2000000 && alive_vertices <= target_vertices*1.15) break;
            continue;
        }
        apply_collapse(keep,rem,shared,cur.pos,cur.newrad);
        fail_streak=0;
    }
}

struct RenderMap {
    int res=0;
    vector<float> depth, nx, ny, nz;
    vector<unsigned char> mask;
};

static inline void project_view(const Vec3&p,int view,int res,double&u,double&v,double&z){
    const double D=2.5;
    const double f=800.0*((double)res/1024.0);
    const double c=0.5*(double)res;
    double sx=0,sy=0;
    if(view==0){ sx=p.y; sy=p.z; z=D-p.x; }
    else if(view==1){ sx=-p.y; sy=p.z; z=D+p.x; }
    else if(view==2){ sx=-p.x; sy=p.z; z=D-p.y; }
    else if(view==3){ sx=p.x; sy=p.z; z=D+p.y; }
    else if(view==4){ sx=p.x; sy=p.y; z=D-p.z; }
    else { sx=-p.x; sy=p.y; z=D+p.z; }
    u=f*sx/z+c; v=f*sy/z+c;
}

static void raster_tri(RenderMap&rm,const Vec3&a,const Vec3&b,const Vec3&c,const Vec3&n,int view){
    int res=rm.res;
    double x0,y0,z0,x1,y1,z1,x2,y2,z2;
    project_view(a,view,res,x0,y0,z0); project_view(b,view,res,x1,y1,z1); project_view(c,view,res,x2,y2,z2);
    if(z0<=0||z1<=0||z2<=0) return;
    int xmin=max(0,(int)floor(min(x0,min(x1,x2))-0.5));
    int xmax=min(res-1,(int)ceil(max(x0,max(x1,x2))+0.5));
    int ymin=max(0,(int)floor(min(y0,min(y1,y2))-0.5));
    int ymax=min(res-1,(int)ceil(max(y0,max(y1,y2))+0.5));
    if(xmin>xmax||ymin>ymax) return;
    double den=(y1-y2)*(x0-x2)+(x2-x1)*(y0-y2);
    if(fabs(den)<1e-14) return;
    for(int yy=ymin; yy<=ymax; ++yy){
        double py=yy+0.5;
        for(int xx=xmin; xx<=xmax; ++xx){
            double px=xx+0.5;
            double w0=((y1-y2)*(px-x2)+(x2-x1)*(py-y2))/den;
            double w1=((y2-y0)*(px-x2)+(x0-x2)*(py-y2))/den;
            double w2=1.0-w0-w1;
            if(w0<-1e-9||w1<-1e-9||w2<-1e-9) continue;
            double zp=1.0/(w0/z0+w1/z1+w2/z2);
            int idx=yy*res+xx;
            if(zp < rm.depth[idx]){
                rm.depth[idx]=(float)zp; rm.nx[idx]=(float)n.x; rm.ny[idx]=(float)n.y; rm.nz[idx]=(float)n.z; rm.mask[idx]=1;
            }
        }
    }
}

static RenderMap make_empty_map(int res){
    RenderMap m; m.res=res; int S=res*res;
    m.depth.assign(S,255.0f); m.nx.assign(S,0.0f); m.ny.assign(S,0.0f); m.nz.assign(S,0.0f); m.mask.assign(S,0);
    return m;
}

static RenderMap render_original_view(int view,int res){
    RenderMap rm=make_empty_map(res);
    for(int i=0;i<M;i++){
        const Face&f=origF[i];
        Vec3 cr=cross3(origP[f.v[1]]-origP[f.v[0]], origP[f.v[2]]-origP[f.v[0]]);
        double len=norm3(cr); if(!(len>0)) continue;
        raster_tri(rm,origP[f.v[0]],origP[f.v[1]],origP[f.v[2]],cr*(1.0/len),view);
    }
    return rm;
}
static RenderMap render_current_view(int view,int res){
    RenderMap rm=make_empty_map(res);
    for(int i=0;i<(int)F.size();i++) if(aliveF[i]){
        const Face&f=F[i];
        if(!aliveV[f.v[0]]||!aliveV[f.v[1]]||!aliveV[f.v[2]]) continue;
        Vec3 cr=face_cross_current(f);
        double len=norm3(cr); if(!(len>0)) continue;
        raster_tri(rm,P[f.v[0]],P[f.v[1]],P[f.v[2]],cr*(1.0/len),view);
    }
    return rm;
}

static vector<RenderMap> orig_cache;
static int orig_cache_res=0;
static const vector<float>* channel_ptr(const RenderMap&m,int ch){
    if(ch==0) return &m.depth;
    if(ch==1) return &m.nx;
    if(ch==2) return &m.ny;
    return &m.nz;
}
static inline float channel_value(float raw,int ch){
    if(ch==0) return raw; // depth
    return (raw+1.0f)*127.5f;
}

static double ssim_one_channel(const RenderMap&a,const RenderMap&b,const vector<unsigned char>&fg,int ch){
    int res=a.res, W=res+1;
    const vector<float>& A=*channel_ptr(a,ch);
    const vector<float>& B=*channel_ptr(b,ch);
    vector<double> sx((size_t)W*W,0), sy((size_t)W*W,0), sxx((size_t)W*W,0), syy((size_t)W*W,0), sxy((size_t)W*W,0);
    for(int y=0;y<res;y++){
        double rowsx=0, rowsy=0, rowsxx=0, rowsyy=0, rowsxy=0;
        size_t base=(size_t)(y+1)*W, prev=(size_t)y*W;
        for(int x=0;x<res;x++){
            double vx=channel_value(A[y*res+x],ch), vy=channel_value(B[y*res+x],ch);
            rowsx+=vx; rowsy+=vy; rowsxx+=vx*vx; rowsyy+=vy*vy; rowsxy+=vx*vy;
            size_t id=base+x+1, up=prev+x+1;
            sx[id]=sx[up]+rowsx; sy[id]=sy[up]+rowsy; sxx[id]=sxx[up]+rowsxx; syy[id]=syy[up]+rowsyy; sxy[id]=sxy[up]+rowsxy;
        }
    }
    auto rect=[&](const vector<double>&S,int x0,int y0,int x1,int y1)->double{
        // inclusive bounds
        size_t A0=(size_t)y0*W+x0, A1=(size_t)y0*W+(x1+1), A2=(size_t)(y1+1)*W+x0, A3=(size_t)(y1+1)*W+(x1+1);
        return S[A3]-S[A1]-S[A2]+S[A0];
    };
    const double C1=(0.01*255.0)*(0.01*255.0), C2=(0.03*255.0)*(0.03*255.0);
    const int rad=5;
    double total=0; int cnt=0;
    for(int y=0;y<res;y++) for(int x=0;x<res;x++){
        int idx=y*res+x; if(!fg[idx]) continue;
        int x0=max(0,x-rad), x1=min(res-1,x+rad), y0=max(0,y-rad), y1=min(res-1,y+rad);
        double n=(double)((x1-x0+1)*(y1-y0+1));
        double X=rect(sx,x0,y0,x1,y1), Y=rect(sy,x0,y0,x1,y1);
        double XX=rect(sxx,x0,y0,x1,y1), YY=rect(syy,x0,y0,x1,y1), XY=rect(sxy,x0,y0,x1,y1);
        double ux=X/n, uy=Y/n;
        double varx=max(0.0,XX/n-ux*ux), vary=max(0.0,YY/n-uy*uy), cov=XY/n-ux*uy;
        double num=(2*ux*uy+C1)*(2*cov+C2);
        double den=(ux*ux+uy*uy+C1)*(varx+vary+C2);
        total += (den!=0? num/den : 1.0);
        cnt++;
    }
    return cnt? total/(double)cnt : 1.0;
}

static double visual_proxy_score(int res){
    if(orig_cache_res!=res){
        orig_cache.clear(); orig_cache.reserve(6);
        for(int v=0;v<6;v++) orig_cache.push_back(render_original_view(v,res));
        orig_cache_res=res;
    }
    double total=0;
    for(int v=0;v<6;v++){
        RenderMap cur=render_current_view(v,res);
        vector<unsigned char> fg(res*res,0);
        for(int i=0;i<res*res;i++) fg[i]=(orig_cache[v].mask[i]||cur.mask[i])?1:0;
        double depth=ssim_one_channel(orig_cache[v],cur,fg,0);
        double norm=(ssim_one_channel(orig_cache[v],cur,fg,1)+ssim_one_channel(orig_cache[v],cur,fg,2)+ssim_one_channel(orig_cache[v],cur,fg,3))/3.0;
        total += 0.5*depth + 0.5*norm;
    }
    return total/6.0;
}

static bool build_compacted(vector<Vec3>&outV, vector<Face>&outF){
    vector<int> id(N,-1); outV.clear(); outF.clear();
    outV.reserve(alive_vertices); outF.reserve(alive_faces);
    double area_eps=max(1e-30,mesh_diag*mesh_diag*1e-24);
    for(int fid=0; fid<(int)F.size(); ++fid) if(aliveF[fid]){
        Face f=F[fid];
        if(f.v[0]<0||f.v[1]<0||f.v[2]<0||f.v[0]>=N||f.v[1]>=N||f.v[2]>=N) return false;
        if(!aliveV[f.v[0]]||!aliveV[f.v[1]]||!aliveV[f.v[2]]) return false;
        if(f.v[0]==f.v[1]||f.v[0]==f.v[2]||f.v[1]==f.v[2]) return false;
        Vec3 cr=face_cross_current(f); if(!(norm2(cr)>area_eps)) return false;
        Face g;
        for(int k=0;k<3;k++){
            int old=f.v[k];
            if(id[old]<0){ id[old]=(int)outV.size(); outV.push_back(P[old]); }
            g.v[k]=id[old];
        }
        outF.push_back(g);
    }
    if(outV.size()<4 || outF.empty() || outV.size()>(size_t)N) return false;
    vector<unsigned long long> edges; edges.reserve(outF.size()*3);
    for(const Face&f: outF){ edges.push_back(edge_key(f.v[0],f.v[1])); edges.push_back(edge_key(f.v[1],f.v[2])); edges.push_back(edge_key(f.v[2],f.v[0])); }
    sort(edges.begin(), edges.end());
    for(size_t i=0;i<edges.size();){ size_t j=i+1; while(j<edges.size()&&edges[j]==edges[i]) ++j; if(j-i!=2) return false; i=j; }
    return true;
}

static double sample_smoothness(){
    // Approximate fraction of sampled edges with dihedral below 25 degrees.
    vector<EdgeFace> ef; ef.reserve(min((size_t)M*3,(size_t)300000));
    int stride=max(1,M/100000);
    for(int i=0;i<M;i+=stride){ const Face&f=origF[i]; ef.push_back({edge_key(f.v[0],f.v[1]),i}); ef.push_back({edge_key(f.v[1],f.v[2]),i}); ef.push_back({edge_key(f.v[2],f.v[0]),i}); }
    sort(ef.begin(),ef.end(),[](const EdgeFace&a,const EdgeFace&b){return a.key<b.key;});
    int total=0,smooth=0; double cos25=cos(25.0*acos(-1.0)/180.0);
    for(size_t i=0;i<ef.size();){ size_t j=i+1; while(j<ef.size()&&ef[j].key==ef[i].key) ++j; if(j-i>=2){
        int f0=ef[i].fid,f1=ef[i+1].fid; Vec3 c0=cross3(origP[origF[f0].v[1]]-origP[origF[f0].v[0]],origP[origF[f0].v[2]]-origP[origF[f0].v[0]]); Vec3 c1=cross3(origP[origF[f1].v[1]]-origP[origF[f1].v[0]],origP[origF[f1].v[2]]-origP[origF[f1].v[0]]); double l0=norm3(c0),l1=norm3(c1); if(l0>0&&l1>0){ double d=dot3(c0,c1)/(l0*l1); total++; if(d>cos25) smooth++; }} i=j; }
    return total? (double)smooth/total : 0.5;
}

static vector<Vec3> bestV;
static vector<Face> bestF;
static bool haveBest=false;

static void save_current_as_best(){
    vector<Vec3> cv; vector<Face> cf;
    if(build_compacted(cv,cf)){
        bestV.swap(cv); bestF.swap(cf); haveBest=true;
    }
}

static void simplify_mesh(){
    t0=chrono::steady_clock::now();
    vector<unsigned long long> initial_edges;
    build_quadrics_and_initial_edges(initial_edges);
    const bool large_fast_mode = (N > 60000);
    if(!large_fast_mode){
        for(unsigned long long key: initial_edges){ int a=(int)(key>>32), b=(int)(key & 0xffffffffu); push_edge_quick(a,b); }
    }
    double smooth=sample_smoothness();

    int proxy_res;
    if(N <= 30000) proxy_res=512;
    else if(N <= 100000) proxy_res=384;
    else if(N <= 300000) proxy_res=192;
    else proxy_res=128;

    // The stages are intentionally aggressive. The image proxy below decides how far to keep.
    vector<double> ratios;
    if(N <= 20) ratios = {0.89,0.80,0.70};
    else if(N < 3000) ratios = {0.55,0.40,0.30,0.22,0.16,0.12,0.09};
    else if(N < 12000) ratios = {0.48,0.34,0.25,0.18,0.135,0.105,0.082};
    else if(N < 80000) ratios = {0.42,0.30,0.22,0.16,0.12,0.092,0.072};
    else if(N < 250000) ratios = {0.38,0.27,0.20,0.15,0.115,0.090,0.072};
    else ratios = {0.34,0.245,0.18,0.135,0.105,0.083,0.066};
    if(smooth > 0.94){ ratios.push_back(0.052); ratios.push_back(0.040); }

    double base_threshold;
    if(proxy_res>=384) base_threshold = 0.908;
    else if(proxy_res>=256) base_threshold = 0.922;
    else base_threshold = 0.942;
    if(smooth < 0.70) base_threshold += 0.018;
    else if(smooth > 0.94) base_threshold -= 0.010;

    int last_alive=N;
    for(size_t si=0; si<ratios.size(); ++si){
        if(elapsed()>16.2) break;
        int target=max(4,(int)ceil(N*ratios[si]));
        if(target>=alive_vertices) continue;
        if(large_fast_mode) fast_decimate_to(target, 16.6);
        else decimate_to(target, 16.6);
        if(alive_vertices >= last_alive) break;
        last_alive=alive_vertices;
        if(elapsed()>18.2){
            // If we are running out of time, keep a locally valid stage. It is usually safe because
            // every accepted collapse obeyed manifold, flip, and Hausdorff guards.
            if(!haveBest) save_current_as_best();
            break;
        }
        double score=visual_proxy_score(proxy_res);
        double threshold=base_threshold;
        if(ratios[si] < 0.075) threshold += 0.008;
        if(ratios[si] < 0.055) threshold += 0.015;
        if(score >= threshold){
            save_current_as_best();
        }else{
            break;
        }
    }
    if(!haveBest){
        // At least try a very conservative compacting pass if proxy never ran/accepted.
        if(alive_vertices < N) save_current_as_best();
    }
}

static void write_mesh_vectors(const vector<Vec3>&Vv,const vector<Face>&Ff){
    string out; out.reserve(1<<20); char line[160];
    auto append=[&](const char* s,int len){ if(out.size()+len>(1<<20)){ fwrite(out.data(),1,out.size(),stdout); out.clear(); } out.append(s,s+len); };
    int len=snprintf(line,sizeof(line),"%d %d\n",(int)Vv.size(),(int)Ff.size()); append(line,len);
    bool highprec = (int)Vv.size()*2 <= N;
    for(const Vec3&p: Vv){
        if(highprec) len=snprintf(line,sizeof(line),"v %.15g %.15g %.15g\n",p.x,p.y,p.z);
        else len=snprintf(line,sizeof(line),"v %.10g %.10g %.10g\n",p.x,p.y,p.z);
        append(line,len);
    }
    for(const Face&f:Ff){ len=snprintf(line,sizeof(line),"f %d %d %d\n",f.v[0]+1,f.v[1]+1,f.v[2]+1); append(line,len); }
    if(!out.empty()) fwrite(out.data(),1,out.size(),stdout);
}
static void write_original(){ write_mesh_vectors(origP,origF); }

int main(){
    load_mesh();
    simplify_mesh();
    if(haveBest) write_mesh_vectors(bestV,bestF);
    else write_original();
    return 0;
}
