#include <bits/stdc++.h>
using namespace std;

struct Vec3{double x=0,y=0,z=0;};
static inline Vec3 operator+(const Vec3&a,const Vec3&b){return {a.x+b.x,a.y+b.y,a.z+b.z};}
static inline Vec3 operator-(const Vec3&a,const Vec3&b){return {a.x-b.x,a.y-b.y,a.z-b.z};}
static inline Vec3 operator*(const Vec3&a,double s){return {a.x*s,a.y*s,a.z*s};}
static inline Vec3 operator/(const Vec3&a,double s){return {a.x/s,a.y/s,a.z/s};}
static inline double dot3(const Vec3&a,const Vec3&b){return a.x*b.x+a.y*b.y+a.z*b.z;}
static inline Vec3 cross3(const Vec3&a,const Vec3&b){return {a.y*b.z-a.z*b.y,a.z*b.x-a.x*b.z,a.x*b.y-a.y*b.x};}
static inline double norm2(const Vec3&a){return dot3(a,a);} 
static inline double norm3(const Vec3&a){return sqrt(max(0.0,norm2(a)));}
static inline double dist2(const Vec3&a,const Vec3&b){return norm2(a-b);} 
static inline bool finite3(const Vec3&p){return isfinite(p.x)&&isfinite(p.y)&&isfinite(p.z);} 
static inline Vec3 normalized(Vec3 a){double n=norm3(a); if(n<=1e-300)return {0,0,0}; return a/n;}

struct Face{int v[3]{};};
struct Candidate{string name; vector<Vec3> V; vector<Face> F; double proxy=-1; bool fromQem=false;};

struct FastScanner{
    vector<char> buf; char* p=nullptr;
    FastScanner(){
        char tmp[1<<16]; size_t n;
        while((n=fread(tmp,1,sizeof(tmp),stdin))>0) buf.insert(buf.end(),tmp,tmp+n);
        buf.push_back('\0'); p=buf.data();
    }
    inline void skip(){ while(*p==' '||*p=='\n'||*p=='\r'||*p=='\t') ++p; }
    string next(){ skip(); string s; while(*p && *p!=' '&&*p!='\n'&&*p!='\r'&&*p!='\t') s.push_back(*p++); return s; }
};
static inline bool isAlphaTok(const string&s){return !s.empty() && (isalpha((unsigned char)s[0]) || s[0]=='#');}
static inline int parseIndex(const string&s){ const char* c=s.c_str(); char* e=nullptr; long v=strtol(c,&e,10); return (int)v; }
static inline double parseDouble(const string&s){ return strtod(s.c_str(),nullptr); }
static inline uint64_t edgeKey(int a,int b){ if(a>b) swap(a,b); return (uint64_t)(uint32_t)a<<32 | (uint32_t)b; }
static inline array<int,3> sortTri(int a,int b,int c){ array<int,3> t{a,b,c}; sort(t.begin(),t.end()); return t; }

static int N0=0,M0=0;
static vector<Vec3> Orig;
static vector<Face> OrigF;
static Vec3 BBmin,BBmax,BBcen;
static double diagLen=1.0, epsLen=0.04975, eps2=0.0, areaEps2=1e-30;
static bool outPrefixed=true;
static chrono::steady_clock::time_point T0;
static inline double elapsed(){ return chrono::duration<double>(chrono::steady_clock::now()-T0).count(); }

static void compute_bbox(){
    if(Orig.empty()){BBmin={0,0,0};BBmax={1,1,1};BBcen={.5,.5,.5};diagLen=1;return;}
    BBmin=BBmax=Orig[0];
    for(const auto&p:Orig){
        BBmin.x=min(BBmin.x,p.x); BBmin.y=min(BBmin.y,p.y); BBmin.z=min(BBmin.z,p.z);
        BBmax.x=max(BBmax.x,p.x); BBmax.y=max(BBmax.y,p.y); BBmax.z=max(BBmax.z,p.z);
    }
    BBcen=(BBmin+BBmax)*0.5;
    diagLen=max(1e-30,norm3(BBmax-BBmin));
    epsLen=diagLen*0.04965;
    eps2=epsLen*epsLen;
    areaEps2=max(1e-32,diagLen*diagLen*diagLen*diagLen*1e-30);
}

static bool read_input(){
    FastScanner in;
    string sN=in.next(); if(sN.empty()) return false;
    N0=parseIndex(sN); M0=parseIndex(in.next());
    if(N0<=0||M0<=0) return false;
    Orig.resize(N0); OrigF.resize(M0);
    outPrefixed=false;
    for(int i=0;i<N0;i++){
        string t=in.next();
        if(isAlphaTok(t)){ outPrefixed=true; t=in.next(); }
        double x=parseDouble(t), y=parseDouble(in.next()), z=parseDouble(in.next());
        Orig[i]={x,y,z};
    }
    for(int i=0;i<M0;i++){
        string t=in.next();
        if(isAlphaTok(t)){ outPrefixed=true; t=in.next(); }
        int a=parseIndex(t), b=parseIndex(in.next()), c=parseIndex(in.next());
        --a;--b;--c;
        OrigF[i].v[0]=a; OrigF[i].v[1]=b; OrigF[i].v[2]=c;
    }
    compute_bbox();
    return true;
}

struct KDTree{
    const vector<Vec3>* pts=nullptr;
    vector<int> id;
    static inline double coord(const Vec3&p,int ax){return ax==0?p.x:(ax==1?p.y:p.z);}    
    void init(const vector<Vec3>&p){pts=&p; id.resize(p.size()); iota(id.begin(),id.end(),0); if(!id.empty()) build(0,(int)id.size(),0);}    
    void build(int l,int r,int d){ if(r-l<=1)return; int m=(l+r)>>1, ax=d%3; nth_element(id.begin()+l,id.begin()+m,id.begin()+r,[&](int a,int b){return coord((*pts)[a],ax)<coord((*pts)[b],ax);}); build(l,m,d+1); build(m+1,r,d+1); }
    void nearest_rec(const Vec3&q,int l,int r,int d,double&best)const{
        if(l>=r)return; int m=(l+r)>>1, ax=d%3; const Vec3&p=(*pts)[id[m]]; double dd=dist2(q,p); if(dd<best)best=dd;
        double df=coord(q,ax)-coord(p,ax); int l1=l,r1=m,l2=m+1,r2=r; if(df>0){swap(l1,l2);swap(r1,r2);} nearest_rec(q,l1,r1,d+1,best); if(df*df<best) nearest_rec(q,l2,r2,d+1,best);
    }
    double nearest2(const Vec3&q)const{ if(id.empty()) return 1e300; double b=1e300; nearest_rec(q,0,(int)id.size(),0,b); return b; }
};

static bool valid_closed_mesh(const vector<Vec3>&V,const vector<Face>&F,bool requireClosed=true){
    if(V.empty()||F.empty()||V.size()>(size_t)max(1,N0)) return false;
    unordered_map<uint64_t,int> ec; ec.reserve(F.size()*4+16);
    for(const Face&f:F){
        int a=f.v[0],b=f.v[1],c=f.v[2];
        if(a<0||b<0||c<0||a>=(int)V.size()||b>=(int)V.size()||c>=(int)V.size()) return false;
        if(a==b||a==c||b==c) return false;
        if(!finite3(V[a])||!finite3(V[b])||!finite3(V[c])) return false;
        if(norm2(cross3(V[b]-V[a],V[c]-V[a]))<=areaEps2) return false;
        ec[edgeKey(a,b)]++; ec[edgeKey(b,c)]++; ec[edgeKey(c,a)]++;
    }
    if(requireClosed){ for(auto &kv:ec) if(kv.second!=2) return false; }
    return true;
}

static bool hausdorff_ok(const vector<Vec3>&V,const vector<Face>&F,double mult=1.0){
    if(!valid_closed_mesh(V,F,true)) return false;
    double lim=eps2*mult*mult;
    static KDTree kdOrig; static bool built=false; if(!built){kdOrig.init(Orig); built=true;}
    KDTree kdCand; kdCand.init(V);
    for(const Vec3&p:Orig) if(kdCand.nearest2(p)>lim) return false;
    for(const Vec3&p:V) if(kdOrig.nearest2(p)>lim) return false;
    return true;
}

static double original_orientation_sign(const Vec3&center){
    double s=0; int stride=max(1,M0/200000), cnt=0;
    for(int i=0;i<M0;i+=stride){
        const Face&f=OrigF[i]; if(f.v[0]<0||f.v[1]<0||f.v[2]<0||f.v[0]>=N0||f.v[1]>=N0||f.v[2]>=N0) continue;
        Vec3 a=Orig[f.v[0]],b=Orig[f.v[1]],c=Orig[f.v[2]];
        Vec3 cr=cross3(b-a,c-a), ctr=(a+b+c)/3.0;
        s+=dot3(cr,ctr-center); cnt++;
    }
    return s>=0?1.0:-1.0;
}
static void orient_face(vector<Vec3>&V,Face&f,const Vec3&center,double sign){
    Vec3 a=V[f.v[0]],b=V[f.v[1]],c=V[f.v[2]], ctr=(a+b+c)/3.0;
    if(dot3(cross3(b-a,c-a),ctr-center)*sign<0) swap(f.v[1],f.v[2]);
}

struct Quadric{
    double q[10];
    Quadric(){ memset(q,0,sizeof(q)); }
    void addPlane(double a,double b,double c,double d,double w=1.0){
        q[0]+=w*a*a; q[1]+=w*a*b; q[2]+=w*a*c; q[3]+=w*a*d; q[4]+=w*b*b;
        q[5]+=w*b*c; q[6]+=w*b*d; q[7]+=w*c*c; q[8]+=w*c*d; q[9]+=w*d*d;
    }
    void add(const Quadric&o){for(int i=0;i<10;i++)q[i]+=o.q[i];}
    double eval(const Vec3&p)const{
        double x=p.x,y=p.y,z=p.z;
        return q[0]*x*x+2*q[1]*x*y+2*q[2]*x*z+2*q[3]*x+q[4]*y*y+2*q[5]*y*z+2*q[6]*y+q[7]*z*z+2*q[8]*z+q[9];
    }
};
struct PQNode{double cost; int a,b,va,vb; bool operator<(const PQNode&o)const{return cost>o.cost;}};
struct EdgeRef{uint64_t key; int f; bool operator<(const EdgeRef&o)const{return key<o.key||(key==o.key&&f<o.f);} };

struct QemSolver{
    int N=0,M=0, activeV=0, activeF=0;
    vector<Vec3> P;
    vector<Face> F;
    vector<unsigned char> alive, factive;
    vector<vector<int>> inc;
    vector<Quadric> Q;
    vector<int> head,tail,nxt,csz,ver,mark;
    priority_queue<PQNode> pq;
    double minCos=0.50;
    double tlimit=12.0;
    vector<int> checkpoints;
    vector<Candidate> out;
    QemSolver(double mc,double tl):minCos(mc),tlimit(tl){}
    inline bool hasv(const Face&f,int v)const{return f.v[0]==v||f.v[1]==v||f.v[2]==v;}
    inline int third(const Face&f,int a,int b)const{for(int k=0;k<3;k++){int x=f.v[k]; if(x!=a&&x!=b)return x;} return -1;}
    bool solve3(const Quadric&q,Vec3&outp){
        double a00=q.q[0],a01=q.q[1],a02=q.q[2],a11=q.q[4],a12=q.q[5],a22=q.q[7];
        double b0=-q.q[3],b1=-q.q[6],b2=-q.q[8];
        double det=a00*(a11*a22-a12*a12)-a01*(a01*a22-a12*a02)+a02*(a01*a12-a11*a02);
        double sc=max(1.0,fabs(a00)+fabs(a11)+fabs(a22)+fabs(a01)+fabs(a02)+fabs(a12));
        if(fabs(det)<1e-16*sc*sc*sc) return false;
        double dx=b0*(a11*a22-a12*a12)-a01*(b1*a22-a12*b2)+a02*(b1*a12-a11*b2);
        double dy=a00*(b1*a22-a12*b2)-b0*(a01*a22-a12*a02)+a02*(a01*b2-b1*a02);
        double dz=a00*(a11*b2-b1*a12)-a01*(a01*b2-b1*a02)+b0*(a01*a12-a11*a02);
        outp={dx/det,dy/det,dz/det}; return finite3(outp);
    }
    bool collectEdge(int u,int v,int by[2],int opp[2]){
        int cnt=0; const vector<int>&A=inc[u].size()<inc[v].size()?inc[u]:inc[v];
        for(int fid:A){ if(!factive[fid]) continue; const Face&f=F[fid]; if(!hasv(f,u)||!hasv(f,v)) continue; if(cnt>=2) return false; by[cnt]=fid; opp[cnt]=third(f,u,v); cnt++; }
        return cnt==2 && opp[0]>=0 && opp[1]>=0 && opp[0]!=opp[1];
    }
    bool linkOK(int u,int v,const int opp[2]){
        static int tok=1; if((int)mark.size()<N)mark.assign(N,0); if(++tok>2000000000){fill(mark.begin(),mark.end(),0);tok=1;} int tu=tok; if(++tok>2000000000){fill(mark.begin(),mark.end(),0);tok=1;} int tc=tok;
        for(int fid:inc[u]){ if(!factive[fid])continue; const Face&f=F[fid]; if(!hasv(f,u))continue; for(int k=0;k<3;k++){int x=f.v[k]; if(x!=u&&x!=v)mark[x]=tu;}}
        int cc=0,g0=0,g1=0;
        for(int fid:inc[v]){ if(!factive[fid])continue; const Face&f=F[fid]; if(!hasv(f,v))continue; for(int k=0;k<3;k++){int x=f.v[k]; if(x==u||x==v)continue; if(mark[x]==tu){mark[x]=tc; cc++; if(x==opp[0])g0=1; if(x==opp[1])g1=1; if(cc>2)return false;}}}
        return cc==2&&g0&&g1;
    }
    bool clusterOK(int v,const Vec3&p){ for(int m=head[v];m!=-1;m=nxt[m]) if(dist2(Orig[m],p)>eps2) return false; return true; }
    bool duplicateTri(int keep,int rem,int fid,int a,int b,int c,int r0,int r1){
        int ix=a; if(inc[b].size()<inc[ix].size())ix=b; if(inc[c].size()<inc[ix].size())ix=c; auto st=sortTri(a,b,c);
        for(int of:inc[ix]){ if(!factive[of]||of==fid||of==r0||of==r1)continue; if(hasv(F[of],rem))continue; if(sortTri(F[of].v[0],F[of].v[1],F[of].v[2])==st) return true; }
        return false;
    }
    bool checkMove(int keep,int rem,const int by[2],const Vec3&pos,double&score){
        if(!clusterOK(keep,pos)||!clusterOK(rem,pos)) return false;
        vector<int> ids; ids.reserve(inc[keep].size()+inc[rem].size());
        for(int f:inc[keep]) if(factive[f]) ids.push_back(f);
        for(int f:inc[rem]) if(factive[f]) ids.push_back(f);
        sort(ids.begin(),ids.end()); ids.erase(unique(ids.begin(),ids.end()),ids.end());
        double worst=0, qcost=0; int changed=0;
        Quadric q=Q[keep]; q.add(Q[rem]); qcost=q.eval(pos);
        for(int fid:ids){
            if(fid==by[0]||fid==by[1]) continue;
            Face f=F[fid]; bool hk=false,hr=false; for(int k=0;k<3;k++){ if(f.v[k]==keep)hk=true; if(f.v[k]==rem)hr=true; }
            if(!hk&&!hr) continue;
            int nv[3]={f.v[0],f.v[1],f.v[2]}; for(int k=0;k<3;k++) if(nv[k]==rem) nv[k]=keep;
            if(nv[0]==nv[1]||nv[0]==nv[2]||nv[1]==nv[2]) return false;
            Vec3 oldp[3]={P[f.v[0]],P[f.v[1]],P[f.v[2]]};
            Vec3 newp[3]={P[nv[0]],P[nv[1]],P[nv[2]]}; for(int k=0;k<3;k++) if(nv[k]==keep) newp[k]=pos;
            Vec3 on=cross3(oldp[1]-oldp[0],oldp[2]-oldp[0]); Vec3 nn=cross3(newp[1]-newp[0],newp[2]-newp[0]);
            double ol=norm2(on), nl=norm2(nn); if(ol<=areaEps2||nl<=areaEps2) return false;
            double d=dot3(on,nn)/sqrt(ol*nl); if(d<minCos) return false; worst=max(worst,1.0-d);
            if(hr && duplicateTri(keep,rem,fid,nv[0],nv[1],nv[2],by[0],by[1])) return false;
            changed++;
        }
        if(changed==0) return false;
        score=qcost + 1e-5*dist2(P[keep],P[rem])/(diagLen*diagLen+1e-30) + 0.02*worst;
        return true;
    }
    bool tryCollapse(int a,int b){
        if(a==b||!alive[a]||!alive[b])return false; if(dist2(P[a],P[b])>4.05*eps2) return false;
        int by[2]={-1,-1}, opp[2]={-1,-1}; if(!collectEdge(a,b,by,opp))return false; if(!linkOK(a,b,opp))return false;
        struct Choice{bool ok=false; int keep=-1,rem=-1; Vec3 pos; double cost=1e100;}; Choice best;
        auto considerDir=[&](int keep,int rem){
            Quadric q=Q[keep]; q.add(Q[rem]); Vec3 opt; Vec3 cand[7]; int cnt=0; if(solve3(q,opt))cand[cnt++]=opt; cand[cnt++]=P[keep]; cand[cnt++]=P[rem]; cand[cnt++]=(P[keep]+P[rem])*0.5; cand[cnt++]=P[keep]*0.75+P[rem]*0.25; cand[cnt++]=P[keep]*0.25+P[rem]*0.75;
            for(int i=0;i<cnt;i++){ double sc; if(checkMove(keep,rem,by,cand[i],sc)&&sc<best.cost){best.ok=true; best.keep=keep; best.rem=rem; best.pos=cand[i]; best.cost=sc;} }
        };
        considerDir(a,b); considerDir(b,a); if(!best.ok)return false; collapse(best.keep,best.rem,by,best.pos); return true;
    }
    double cheapCost(int a,int b){ Quadric q=Q[a]; q.add(Q[b]); Vec3 o; double r=1e100; if(solve3(q,o)) r=min(r,q.eval(o)); r=min(r,q.eval((P[a]+P[b])*0.5)); r=min(r,q.eval(P[a])); r=min(r,q.eval(P[b])); return r+1e-6*dist2(P[a],P[b])/(diagLen*diagLen+1e-30); }
    void pushEdge(int a,int b){ if(a==b||a<0||b<0||a>=N||b>=N||!alive[a]||!alive[b])return; if(dist2(P[a],P[b])>4.05*eps2)return; pq.push({cheapCost(a,b),a,b,ver[a],ver[b]}); }
    void compact(int v){ if(v<0||v>=N)return; auto &x=inc[v]; if(x.size()<128)return; vector<int> y; y.reserve(x.size()); for(int f:x) if(factive[f]&&hasv(F[f],v)) y.push_back(f); if(y.size()*3+32<x.size()) x.swap(y); }
    void mergeMembers(int keep,int rem){ if(head[rem]==-1)return; nxt[tail[keep]]=head[rem]; tail[keep]=tail[rem]; csz[keep]+=csz[rem]; head[rem]=tail[rem]=-1; csz[rem]=0; }
    void collapse(int keep,int rem,const int by[2],const Vec3&pos){
        Q[keep].add(Q[rem]); P[keep]=pos;
        for(int i=0;i<2;i++) if(by[i]>=0&&factive[by[i]]){factive[by[i]]=0; activeF--;}
        for(int fid:inc[rem]){
            if(!factive[fid])continue; Face&f=F[fid]; if(!hasv(f,rem))continue; if(hasv(f,keep)){factive[fid]=0; activeF--; continue;} for(int k=0;k<3;k++) if(f.v[k]==rem) f.v[k]=keep; inc[keep].push_back(fid);
        }
        alive[rem]=0; activeV--; ver[keep]++; ver[rem]++; mergeMembers(keep,rem); compact(rem); compact(keep);
        for(int fid:inc[keep]) if(factive[fid]){Face&f=F[fid]; if(!hasv(f,keep))continue; pushEdge(f.v[0],f.v[1]); pushEdge(f.v[1],f.v[2]); pushEdge(f.v[2],f.v[0]);}
    }
    Candidate exportMesh(const string&name){
        Candidate c; c.name=name; c.fromQem=true; vector<int> remap(N,-1); c.V.reserve(activeV); for(int i=0;i<N;i++) if(alive[i]){remap[i]=(int)c.V.size(); c.V.push_back(P[i]);}
        for(int i=0;i<M;i++) if(factive[i]){int a=F[i].v[0],b=F[i].v[1],d=F[i].v[2]; if(a==b||a==d||b==d)continue; if(a<0||b<0||d<0||a>=N||b>=N||d>=N)continue; if(remap[a]<0||remap[b]<0||remap[d]<0)continue; Face f; f.v[0]=remap[a]; f.v[1]=remap[b]; f.v[2]=remap[d]; c.F.push_back(f);} return c;
    }
    void init(){
        N=N0; M=M0; P=Orig; F=OrigF; activeV=N; activeF=M; alive.assign(N,1); factive.assign(M,1); inc.assign(N,{}); vector<int> deg(N,0);
        for(auto &f:F){ for(int k=0;k<3;k++) if(f.v[k]>=0&&f.v[k]<N) deg[f.v[k]]++; }
        for(int i=0;i<N;i++) inc[i].reserve(deg[i]+8);
        for(int i=0;i<M;i++) for(int k=0;k<3;k++) if(F[i].v[k]>=0&&F[i].v[k]<N) inc[F[i].v[k]].push_back(i);
        Q.assign(N,Quadric()); head.resize(N); tail.resize(N); nxt.assign(N,-1); csz.assign(N,1); ver.assign(N,0); mark.assign(N,0); for(int i=0;i<N;i++)head[i]=tail[i]=i;
        vector<EdgeRef> er; er.reserve((size_t)M*3);
        for(int i=0;i<M;i++){
            Face&f=F[i]; if(f.v[0]<0||f.v[1]<0||f.v[2]<0||f.v[0]>=N||f.v[1]>=N||f.v[2]>=N) continue;
            Vec3 n=cross3(P[f.v[1]]-P[f.v[0]],P[f.v[2]]-P[f.v[0]]); double l=norm3(n); if(l>0){n=n/l; double d=-dot3(n,P[f.v[0]]); Q[f.v[0]].addPlane(n.x,n.y,n.z,d); Q[f.v[1]].addPlane(n.x,n.y,n.z,d); Q[f.v[2]].addPlane(n.x,n.y,n.z,d);}            
            er.push_back({edgeKey(f.v[0],f.v[1]),i}); er.push_back({edgeKey(f.v[1],f.v[2]),i}); er.push_back({edgeKey(f.v[2],f.v[0]),i});
        }
        sort(er.begin(),er.end());
        for(size_t i=0;i<er.size();){size_t j=i+1; while(j<er.size()&&er[j].key==er[i].key)j++; int a=(int)(er[i].key>>32), b=(int)(er[i].key&0xffffffffu); pushEdge(a,b); i=j;}
    }
    vector<Candidate> run(){
        init(); if(N<20||M<20)return {};
        sort(checkpoints.rbegin(),checkpoints.rend()); int idx=0; int minTarget=checkpoints.empty()?max(4,N/10):checkpoints.back(); long long pops=0;
        while(activeV>minTarget&&!pq.empty()&&elapsed()<tlimit){
            PQNode n=pq.top(); pq.pop(); if(n.a<0||n.b<0||n.a>=N||n.b>=N)continue; if(!alive[n.a]||!alive[n.b])continue; if(n.va!=ver[n.a]||n.vb!=ver[n.b]){pushEdge(n.a,n.b); continue;} tryCollapse(n.a,n.b); pops++; if((pops&8191)==0&&elapsed()>tlimit)break;
            while(idx<(int)checkpoints.size()&&activeV<=checkpoints[idx]){ Candidate c=exportMesh("qem"+to_string(checkpoints[idx])); if(!c.V.empty()&&!c.F.empty()) out.push_back(move(c)); idx++; }
        }
        Candidate c=exportMesh("qem_final"); if(!c.V.empty()&&!c.F.empty()) out.push_back(move(c));
        return out;
    }
};

struct Image{int r=0; vector<float>d,nx,ny,nz; vector<unsigned char>occ;};
struct Proxy{
    int res; Vec3 cen; double scale; vector<Image> ref;
    Proxy(int r):res(r){
        cen=BBcen; double ex=BBmax.x-BBmin.x,ey=BBmax.y-BBmin.y,ez=BBmax.z-BBmin.z; scale=max(ex,max(ey,ez))*0.5; if(!(scale>1e-300))scale=1;
        for(int v=0;v<6;v++) ref.push_back(render(Orig,OrigF,v));
    }
    inline Vec3 normp(const Vec3&p)const{return (p-cen)/scale;}
    void project(const Vec3&p,int view,double&x,double&y,double&z)const{
        Vec3 q=normp(p); double sx,sy,sd;
        if(view==0){sx=q.y;sy=q.z;sd=1.08-q.x;} else if(view==1){sx=-q.y;sy=q.z;sd=1.08+q.x;} else if(view==2){sx=-q.x;sy=q.z;sd=1.08-q.y;} else if(view==3){sx=q.x;sy=q.z;sd=1.08+q.y;} else if(view==4){sx=q.x;sy=q.y;sd=1.08-q.z;} else {sx=-q.x;sy=q.y;sd=1.08+q.z;}
        const double margin=1.08; x=(sx/margin*0.5+0.5)*(res-1); y=(sy/margin*0.5+0.5)*(res-1); z=max(0.0,min(2.16,sd));
    }
    void tri(Image&im,const Vec3&A,const Vec3&B,const Vec3&C,const Vec3&n,int view)const{
        double x0,y0,z0,x1,y1,z1,x2,y2,z2; project(A,view,x0,y0,z0); project(B,view,x1,y1,z1); project(C,view,x2,y2,z2);
        int xmin=max(0,(int)floor(min(x0,min(x1,x2)))); int xmax=min(res-1,(int)ceil(max(x0,max(x1,x2))));
        int ymin=max(0,(int)floor(min(y0,min(y1,y2)))); int ymax=min(res-1,(int)ceil(max(y0,max(y1,y2)))); if(xmin>xmax||ymin>ymax)return;
        double den=(y1-y2)*(x0-x2)+(x2-x1)*(y0-y2); if(fabs(den)<1e-14)return;
        for(int yy=ymin;yy<=ymax;yy++){double py=yy+0.5; for(int xx=xmin;xx<=xmax;xx++){double px=xx+0.5; double w0=((y1-y2)*(px-x2)+(x2-x1)*(py-y2))/den; double w1=((y2-y0)*(px-x2)+(x0-x2)*(py-y2))/den; double w2=1-w0-w1; if(w0<-1e-9||w1<-1e-9||w2<-1e-9)continue; double zz=w0*z0+w1*z1+w2*z2; int id=yy*res+xx; if(zz<im.d[id]){im.d[id]=(float)zz; im.nx[id]=(float)n.x; im.ny[id]=(float)n.y; im.nz[id]=(float)n.z; im.occ[id]=1;}}}
    }
    Image render(const vector<Vec3>&V,const vector<Face>&F,int view)const{
        Image im; im.r=res; int n=res*res; im.d.assign(n,2.16f); im.nx.assign(n,0); im.ny.assign(n,0); im.nz.assign(n,0); im.occ.assign(n,0);
        for(const Face&f:F){ if(f.v[0]<0||f.v[1]<0||f.v[2]<0||f.v[0]>=(int)V.size()||f.v[1]>=(int)V.size()||f.v[2]>=(int)V.size())continue; Vec3 a=V[f.v[0]],b=V[f.v[1]],c=V[f.v[2]]; Vec3 cr=cross3(b-a,c-a); double l=norm3(cr); if(l<=1e-300)continue; tri(im,a,b,c,cr/l,view); }
        return im;
    }
    static vector<double> integral(const vector<double>&a,int r){ vector<double>I((r+1)*(r+1),0); for(int y=0;y<r;y++){double rows=0; for(int x=0;x<r;x++){rows+=a[y*r+x]; I[(y+1)*(r+1)+x+1]=I[y*(r+1)+x+1]+rows;}} return I; }
    static inline double rect(const vector<double>&I,int r,int x0,int y0,int x1,int y1){int W=r+1; return I[y1*W+x1]-I[y0*W+x1]-I[y1*W+x0]+I[y0*W+x0];}
    double ssim(const vector<double>&A,const vector<double>&B,const vector<unsigned char>&mask)const{
        int r=res,n=r*r; vector<double>A2(n),B2(n),AB(n); for(int i=0;i<n;i++){A2[i]=A[i]*A[i]; B2[i]=B[i]*B[i]; AB[i]=A[i]*B[i];}
        vector<double>IA=integral(A,r),IB=integral(B,r),IA2=integral(A2,r),IB2=integral(B2,r),IAB=integral(AB,r);
        const double C1=6.5025,C2=58.5225; double tot=0; int cnt=0; int rad=5;
        for(int y=0;y<r;y++)for(int x=0;x<r;x++){int id=y*r+x; if(!mask[id])continue; int x0=max(0,x-rad),x1=min(r,x+rad+1),y0=max(0,y-rad),y1=min(r,y+rad+1); double inv=1.0/((x1-x0)*(y1-y0)); double ma=rect(IA,r,x0,y0,x1,y1)*inv, mb=rect(IB,r,x0,y0,x1,y1)*inv; double va=max(0.0,rect(IA2,r,x0,y0,x1,y1)*inv-ma*ma), vb=max(0.0,rect(IB2,r,x0,y0,x1,y1)*inv-mb*mb), cov=rect(IAB,r,x0,y0,x1,y1)*inv-ma*mb; double val=((2*ma*mb+C1)*(2*cov+C2))/((ma*ma+mb*mb+C1)*(va+vb+C2)); if(isfinite(val)){tot+=val;cnt++;}}
        return cnt?tot/cnt:1.0;
    }
    double score(const vector<Vec3>&V,const vector<Face>&F){
        double total=0;
        for(int view=0;view<6;view++){
            Image b=render(V,F,view); const Image&a=ref[view]; int n=res*res; vector<unsigned char>mask(n); int uni=0,inter=0;
            vector<double>da(n),db(n),na(n),nb(n);
            for(int i=0;i<n;i++){ mask[i]=a.occ[i]||b.occ[i]; if(mask[i])uni++; if(a.occ[i]&&b.occ[i])inter++; da[i]=a.occ[i]?a.d[i]*(255.0/2.16):255.0; db[i]=b.occ[i]?b.d[i]*(255.0/2.16):255.0; }
            double depth=ssim(da,db,mask), normal=0;
            for(int ch=0;ch<3;ch++){ for(int i=0;i<n;i++){ float av= ch==0?a.nx[i]:(ch==1?a.ny[i]:a.nz[i]); float bv= ch==0?b.nx[i]:(ch==1?b.ny[i]:b.nz[i]); na[i]=a.occ[i]?(av+1.0)*127.5:127.5; nb[i]=b.occ[i]?(bv+1.0)*127.5:127.5; } normal+=ssim(na,nb,mask); }
            normal/=3.0; double occ=uni?((double)inter/(double)uni):1.0; double viewscore=.45*normal+.45*depth+.10*occ; if(occ<.92)viewscore-=.18*(.92-occ); total+=viewscore;
        }
        return total/6.0;
    }
};

static bool add_candidate(vector<Candidate>&C,Candidate&&c,const string&name,bool qem=false,bool doHaus=true){
    c.name=name; c.fromQem=qem; if(c.V.empty()||c.F.empty())return false; if(c.V.size()>=Orig.size())return false; if(c.V.size()>Orig.size())return false; if(!valid_closed_mesh(c.V,c.F,true))return false; if(doHaus&&!hausdorff_ok(c.V,c.F,1.0))return false; C.push_back(move(c)); return true;
}

static Candidate makeSphere(const Vec3&ctr,double r,int lat,int lon,const string&nm){
    Candidate c; c.name=nm; if(lat<4||lon<8)return c; int vc=2+(lat-1)*lon; if(vc>N0)return c; c.V.reserve(vc); c.F.reserve(2*lat*lon); double pi=acos(-1.0); c.V.push_back({ctr.x,ctr.y,ctr.z+r}); for(int i=1;i<lat;i++){double th=pi*i/lat, st=sin(th), ct=cos(th); for(int j=0;j<lon;j++){double ph=2*pi*j/lon; c.V.push_back({ctr.x+r*st*cos(ph),ctr.y+r*st*sin(ph),ctr.z+r*ct});}} c.V.push_back({ctr.x,ctr.y,ctr.z-r}); int bot=(int)c.V.size()-1; auto id=[&](int ring,int j){return 1+(ring-1)*lon+(j%lon+lon)%lon;}; double s=original_orientation_sign(ctr); auto add=[&](int a,int b,int d){Face f;f.v[0]=a;f.v[1]=b;f.v[2]=d; orient_face(c.V,f,ctr,s); c.F.push_back(f);}; for(int j=0;j<lon;j++)add(0,id(1,j+1),id(1,j)); for(int i=1;i<lat-1;i++)for(int j=0;j<lon;j++){int a=id(i,j),b=id(i,j+1),d=id(i+1,j),e=id(i+1,j+1); add(a,b,d); add(b,e,d);} for(int j=0;j<lon;j++)add(bot,id(lat-1,j),id(lat-1,j+1)); return c;
}

static void eigenSym(double A[3][3],Vec3 out[3]){
    double v[3][3]={{1,0,0},{0,1,0},{0,0,1}};
    for(int it=0;it<40;it++){int p=0,q=1; double best=fabs(A[0][1]); if(fabs(A[0][2])>best){p=0;q=2;best=fabs(A[0][2]);} if(fabs(A[1][2])>best){p=1;q=2;best=fabs(A[1][2]);} if(best<1e-18)break; double app=A[p][p],aqq=A[q][q],apq=A[p][q]; double tau=(aqq-app)/(2*apq); double t=(tau>=0?1:-1)/(fabs(tau)+sqrt(1+tau*tau)); double co=1/sqrt(1+t*t), si=t*co; for(int k=0;k<3;k++) if(k!=p&&k!=q){double akp=A[k][p],akq=A[k][q]; A[k][p]=A[p][k]=co*akp-si*akq; A[k][q]=A[q][k]=si*akp+co*akq;} A[p][p]=co*co*app-2*si*co*apq+si*si*aqq; A[q][q]=si*si*app+2*si*co*apq+co*co*aqq; A[p][q]=A[q][p]=0; for(int k=0;k<3;k++){double vp=v[k][p],vq=v[k][q]; v[k][p]=co*vp-si*vq; v[k][q]=si*vp+co*vq;}}
    int ord[3]={0,1,2}; sort(ord,ord+3,[&](int a,int b){return A[a][a]>A[b][b];}); for(int j=0;j<3;j++){int c=ord[j]; out[j]=normalized({v[0][c],v[1][c],v[2][c]});} if(dot3(cross3(out[0],out[1]),out[2])<0)out[2]=out[2]*-1;
}
struct EllFit{bool ok=false; Vec3 c,axis[3]; double r[3]; double rms=1e9,maxe=1e9;};
static EllFit fitEllipsoid(bool pca){
    EllFit fit; fit.axis[0]={1,0,0};fit.axis[1]={0,1,0};fit.axis[2]={0,0,1};
    Vec3 mean{}; for(auto&p:Orig)mean=mean+p; mean=mean/(double)max(1,N0);
    if(pca){double C[3][3]={}; int stride=max(1,N0/250000),cnt=0; for(int i=0;i<N0;i+=stride){Vec3 q=Orig[i]-mean; double x[3]={q.x,q.y,q.z}; for(int a=0;a<3;a++)for(int b=0;b<3;b++)C[a][b]+=x[a]*x[b];cnt++;} for(int a=0;a<3;a++)for(int b=0;b<3;b++)C[a][b]/=max(1,cnt); eigenSym(C,fit.axis);}    
    double lo[3]={1e100,1e100,1e100},hi[3]={-1e100,-1e100,-1e100}; for(const auto&p:Orig){for(int k=0;k<3;k++){double t=dot3(p,fit.axis[k]); lo[k]=min(lo[k],t); hi[k]=max(hi[k],t);}}
    fit.c={0,0,0}; for(int k=0;k<3;k++){double mid=.5*(lo[k]+hi[k]); fit.c=fit.c+fit.axis[k]*mid; fit.r[k]=.5*(hi[k]-lo[k]); if(fit.r[k]<=1e-12)return fit;}
    int stride=max(1,N0/250000),cnt=0; double ss=0,mx=0; for(int i=0;i<N0;i+=stride){Vec3 q=Orig[i]-fit.c; double r2=0; for(int k=0;k<3;k++){double u=dot3(q,fit.axis[k])/fit.r[k]; r2+=u*u;} double e=fabs(sqrt(max(0.0,r2))-1); ss+=e*e; mx=max(mx,e); cnt++;}
    fit.rms=sqrt(ss/max(1,cnt)); fit.maxe=mx; fit.ok=fit.rms<0.035&&fit.maxe<0.16; return fit;
}
static Candidate makeEllipsoid(const EllFit&fit,int lat,int lon,const string&nm){
    Candidate c; c.name=nm; if(!fit.ok||lat<4||lon<8)return c; int vc=2+(lat-1)*lon; if(vc>N0)return c; double pi=acos(-1.0); auto mp=[&](double x,double y,double z){return fit.c+fit.axis[0]*(fit.r[0]*x)+fit.axis[1]*(fit.r[1]*y)+fit.axis[2]*(fit.r[2]*z);}; c.V.push_back(mp(0,0,1)); for(int i=1;i<lat;i++){double th=pi*i/lat,st=sin(th),ct=cos(th); for(int j=0;j<lon;j++){double ph=2*pi*j/lon; c.V.push_back(mp(st*cos(ph),st*sin(ph),ct));}} c.V.push_back(mp(0,0,-1)); int bot=(int)c.V.size()-1; auto id=[&](int r,int j){return 1+(r-1)*lon+(j%lon+lon)%lon;}; double s=original_orientation_sign(fit.c); auto add=[&](int a,int b,int d){Face f;f.v[0]=a;f.v[1]=b;f.v[2]=d; orient_face(c.V,f,fit.c,s); c.F.push_back(f);}; for(int j=0;j<lon;j++)add(0,id(1,j+1),id(1,j)); for(int i=1;i<lat-1;i++)for(int j=0;j<lon;j++){int a=id(i,j),b=id(i,j+1),d=id(i+1,j),e=id(i+1,j+1);add(a,b,d);add(b,e,d);} for(int j=0;j<lon;j++)add(bot,id(lat-1,j),id(lat-1,j+1)); return c;
}

static void basisForAxis(int ax,Vec3&w,Vec3&u,Vec3&v){ if(ax==0){w={1,0,0};u={0,1,0};v={0,0,1};} else if(ax==1){w={0,1,0};u={1,0,0};v={0,0,1};} else {w={0,0,1};u={1,0,0};v={0,1,0};} }
struct RevFit{bool ok=false; int ax=2; double t0,t1,cu,cv,r0,r1,res;};
static RevFit fitRevolve(int ax){
    RevFit fit; fit.ax=ax; Vec3 w,u,v; basisForAxis(ax,w,u,v); double t0=1e100,t1=-1e100,lu=1e100,hu=-1e100,lv=1e100,hv=-1e100; for(auto&p:Orig){double t=dot3(p,w),x=dot3(p,u),y=dot3(p,v); t0=min(t0,t);t1=max(t1,t);lu=min(lu,x);hu=max(hu,x);lv=min(lv,y);hv=max(hv,y);} if(t1<=t0)return fit; fit.t0=t0;fit.t1=t1;fit.cu=.5*(lu+hu);fit.cv=.5*(lv+hv); double S=0,St=0,Stt=0,Sr=0,Str=0,maxr=0; int cnt=0,stride=max(1,N0/250000); for(int i=0;i<N0;i+=stride){double t=dot3(Orig[i],w),x=dot3(Orig[i],u)-fit.cu,y=dot3(Orig[i],v)-fit.cv,r=sqrt(x*x+y*y); if(r<1e-9)continue; S++;St+=t;Stt+=t*t;Sr+=r;Str+=t*r;maxr=max(maxr,r);cnt++;} if(cnt<100||maxr<=1e-12)return fit; double den=S*Stt-St*St; if(fabs(den)<1e-18)return fit; double a=(S*Str-St*Sr)/den,b=(Sr-a*St)/S; fit.r0=max(0.0,a*t0+b); fit.r1=max(0.0,a*t1+b); if(max(fit.r0,fit.r1)<maxr*.25)return fit; double ss=0,mx=0;cnt=0; for(int i=0;i<N0;i+=stride){double t=dot3(Orig[i],w),x=dot3(Orig[i],u)-fit.cu,y=dot3(Orig[i],v)-fit.cv,r=sqrt(x*x+y*y),pr=max(0.0,a*t+b); double e=fabs(r-pr); ss+=e*e; mx=max(mx,e);cnt++;} double rms=sqrt(ss/max(1,cnt)); fit.res=rms/maxr; fit.ok=fit.res<0.035 && mx<maxr*.13; return fit;
}
static Candidate makeRevolve(const RevFit&fit,int sides,const string&nm){
    Candidate c; c.name=nm; if(!fit.ok||sides<8)return c; Vec3 w,u,v; basisForAxis(fit.ax,w,u,v); Vec3 center=w*((fit.t0+fit.t1)*.5)+u*fit.cu+v*fit.cv; double s=original_orientation_sign(center); double pi=acos(-1.0); bool cone0=fit.r0<max(fit.r0,fit.r1)*.04, cone1=fit.r1<max(fit.r0,fit.r1)*.04; auto point=[&](double t,double r,int j){double th=2*pi*j/sides; return w*t+u*(fit.cu+r*cos(th))+v*(fit.cv+r*sin(th));}; auto add=[&](int a,int b,int d){Face f;f.v[0]=a;f.v[1]=b;f.v[2]=d; orient_face(c.V,f,center,s); c.F.push_back(f);};
    if(!cone0&&!cone1){c.V.push_back(w*fit.t0+u*fit.cu+v*fit.cv); c.V.push_back(w*fit.t1+u*fit.cu+v*fit.cv); int r0=(int)c.V.size(); for(int j=0;j<sides;j++)c.V.push_back(point(fit.t0,fit.r0,j)); int r1=(int)c.V.size(); for(int j=0;j<sides;j++)c.V.push_back(point(fit.t1,fit.r1,j)); if(c.V.size()>(size_t)N0){c.V.clear();c.F.clear();return c;} for(int j=0;j<sides;j++){int k=(j+1)%sides; add(r0+j,r0+k,r1+j); add(r0+k,r1+k,r1+j); add(0,r0+j,r0+k); add(1,r1+k,r1+j);} }
    else{ bool bot=cone0; double at=bot?fit.t0:fit.t1, bt=bot?fit.t1:fit.t0, br=bot?fit.r1:fit.r0; c.V.push_back(w*at+u*fit.cu+v*fit.cv); c.V.push_back(w*bt+u*fit.cu+v*fit.cv); int r=(int)c.V.size(); for(int j=0;j<sides;j++)c.V.push_back(point(bt,br,j)); if(c.V.size()>(size_t)N0){c.V.clear();c.F.clear();return c;} for(int j=0;j<sides;j++){int k=(j+1)%sides; add(0,r+j,r+k); add(1,r+k,r+j);} }
    return c;
}

struct TorFit{bool ok=false; int ax=2; double ct,cu,cv,R,r,res;};
static TorFit fitTorus(int ax){
    TorFit fit; fit.ax=ax; Vec3 w,u,v; basisForAxis(ax,w,u,v); double t0=1e100,t1=-1e100,lu=1e100,hu=-1e100,lv=1e100,hv=-1e100; for(auto&p:Orig){double t=dot3(p,w),x=dot3(p,u),y=dot3(p,v); t0=min(t0,t);t1=max(t1,t);lu=min(lu,x);hu=max(hu,x);lv=min(lv,y);hv=max(hv,y);} fit.ct=.5*(t0+t1); fit.cu=.5*(lu+hu); fit.cv=.5*(lv+hv); int stride=max(1,N0/250000),cnt=0; double minrho=1e100,maxrho=0; for(int i=0;i<N0;i+=stride){double x=dot3(Orig[i],u)-fit.cu,y=dot3(Orig[i],v)-fit.cv,rho=sqrt(x*x+y*y); minrho=min(minrho,rho);maxrho=max(maxrho,rho);cnt++;} if(maxrho<=1e-12||minrho<maxrho*.20)return fit; fit.R=.5*(minrho+maxrho); double sum=0;cnt=0; for(int i=0;i<N0;i+=stride){double t=dot3(Orig[i],w)-fit.ct,x=dot3(Orig[i],u)-fit.cu,y=dot3(Orig[i],v)-fit.cv,rho=sqrt(x*x+y*y); sum+=sqrt((rho-fit.R)*(rho-fit.R)+t*t);cnt++;} fit.r=sum/max(1,cnt); if(fit.r<=1e-12||fit.R<fit.r*1.4)return fit; double ss=0,mx=0;cnt=0; for(int i=0;i<N0;i+=stride){double t=dot3(Orig[i],w)-fit.ct,x=dot3(Orig[i],u)-fit.cu,y=dot3(Orig[i],v)-fit.cv,rho=sqrt(x*x+y*y); double e=fabs(sqrt((rho-fit.R)*(rho-fit.R)+t*t)-fit.r); ss+=e*e;mx=max(mx,e);cnt++;} fit.res=sqrt(ss/max(1,cnt))/fit.r; fit.ok=fit.res<0.045&&mx<fit.r*.20; return fit;
}
static Candidate makeTorus(const TorFit&fit,int A,int B,const string&nm){
    Candidate c; c.name=nm; if(!fit.ok||A<12||B<6||A*B>N0)return c; Vec3 w,u,v; basisForAxis(fit.ax,w,u,v); Vec3 center=w*fit.ct+u*fit.cu+v*fit.cv; double s=original_orientation_sign(center),pi=acos(-1.0); c.V.reserve(A*B); for(int i=0;i<A;i++){double th=2*pi*i/A,ct=cos(th),st=sin(th); for(int j=0;j<B;j++){double ph=2*pi*j/B,cp=cos(ph),sp=sin(ph),rad=fit.R+fit.r*cp; c.V.push_back(w*(fit.ct+fit.r*sp)+u*(fit.cu+rad*ct)+v*(fit.cv+rad*st));}} auto id=[&](int i,int j){return ((i%A+A)%A)*B+((j%B+B)%B);}; auto add=[&](int a,int b,int d){Face f;f.v[0]=a;f.v[1]=b;f.v[2]=d; orient_face(c.V,f,center,s); c.F.push_back(f);}; for(int i=0;i<A;i++)for(int j=0;j<B;j++){int a=id(i,j),b=id(i+1,j),d=id(i+1,j+1),e=id(i,j+1); add(a,b,d); add(a,d,e);} return c;
}

static bool faceMatches(const Face&f,int a,int b,int c){return sortTri(f.v[0],f.v[1],f.v[2])==sortTri(a,b,c);} 
static bool detectSphereGrid(int&R,int&V){
    if(N0<30||M0!=2*(N0-2))return false;
    for(int v=8;v<=2048;v++){ if((N0-2)%v)continue; int r=(N0-2)/v; if(r<2)continue; int step=max(1,v/64); bool ok=true; for(int j=0;j<v&&ok;j+=step){int a=1+j,b=1+(j+1)%v; if(!faceMatches(OrigF[j],0,b,a))ok=false; int off=v+2*(r-1)*v+j; int c=1+(r-1)*v+j,d=1+(r-1)*v+(j+1)%v; if(ok&&off<M0&&!faceMatches(OrigF[off],N0-1,c,d))ok=false;} if(ok){R=r;V=v;return true;} }
    return false;
}
static Candidate makeSphereGrid(int R,int V,int R2,int V2,const string&nm){
    Candidate c; c.name=nm; if(R2<2||V2<8||2+R2*V2>=N0)return c; c.V.push_back(Orig[0]); for(int i=0;i<R2;i++){int oi=1+(int)((long long)i*(R-1)/max(1,R2-1)); for(int j=0;j<V2;j++){int oj=(int)((long long)j*V/V2); c.V.push_back(Orig[1+(oi-1)*V+oj]);}} int bot=(int)c.V.size(); c.V.push_back(Orig[N0-1]); Vec3 center=BBcen; double s=original_orientation_sign(center); auto id=[&](int r,int j){return 1+r*V2+(j%V2+V2)%V2;}; auto add=[&](int a,int b,int d){Face f;f.v[0]=a;f.v[1]=b;f.v[2]=d; orient_face(c.V,f,center,s); c.F.push_back(f);}; for(int j=0;j<V2;j++)add(0,id(0,j+1),id(0,j)); for(int r=0;r<R2-1;r++)for(int j=0;j<V2;j++){int a=id(r,j),b=id(r,j+1),d=id(r+1,j),e=id(r+1,j+1); add(a,b,d); add(b,e,d);} for(int j=0;j<V2;j++)add(bot,id(R2-1,j),id(R2-1,j+1)); return c;
}

static void generate_specials(vector<Candidate>&C){
    // Sphere and ellipsoid families. Hausdorff validation is exact vertex-wise, so detectors can be generous.
    Vec3 ctr=BBcen; int stride=max(1,N0/250000),cnt=0; double sr=0,sr2=0; for(int i=0;i<N0;i+=stride){double r=norm3(Orig[i]-ctr); sr+=r; sr2+=r*r; cnt++;} double rad=cnt?sr/cnt:1; if(rad>1e-12){ double ss=0,mx=0; for(int i=0;i<N0;i+=stride){double e=fabs(norm3(Orig[i]-ctr)-rad); ss+=e*e; mx=max(mx,e);} double rel=sqrt(ss/max(1,cnt))/rad; if(rel<0.08){ int trials[][2]={{10,20},{12,24},{16,32},{20,40},{24,48},{28,56}}; for(auto&t:trials){ if(2+(t[0]-1)*t[1]>=N0)continue; Candidate c=makeSphere(ctr,rad,t[0],t[1],"sphere"); add_candidate(C,move(c),"sphere",false,true); } } }
    EllFit e0=fitEllipsoid(false), e1=fitEllipsoid(true); int etr[][2]={{12,24},{16,32},{20,40},{24,48},{28,56}}; for(auto &fit:{e0,e1}) if(fit.ok){ for(auto&t:etr){ if(2+(t[0]-1)*t[1]>=N0)continue; Candidate c=makeEllipsoid(fit,t[0],t[1],"ellipsoid"); add_candidate(C,move(c),"ellipsoid",false,true); } }
    for(int ax=0;ax<3;ax++){ RevFit rf=fitRevolve(ax); if(rf.ok){int sidesA[]={16,24,32,40,48,64}; for(int s:sidesA){Candidate c=makeRevolve(rf,s,"revolve"); add_candidate(C,move(c),"revolve",false,true);}} TorFit tf=fitTorus(ax); if(tf.ok){int tr[][2]={{32,10},{48,12},{64,16},{80,20},{96,24}}; for(auto&t:tr){Candidate c=makeTorus(tf,t[0],t[1],"torus"); add_candidate(C,move(c),"torus",false,true);}} }
    int R=0,V=0; if(detectSphereGrid(R,V)){ int rr[]={max(2,R/8),max(3,R/6),max(4,R/4),max(5,R/3)}; int vv[]={max(8,V/8),max(12,V/6),max(16,V/4),max(20,V/3)}; for(int i=0;i<4;i++){Candidate c=makeSphereGrid(R,V,rr[i],vv[i],"sphere_grid"); add_candidate(C,move(c),"sphere_grid",false,true);} }
}

static void print_mesh(const vector<Vec3>&V,const vector<Face>&F){
    static char obuf[1<<20]; setvbuf(stdout,obuf,_IOFBF,sizeof(obuf));
    printf("%d %d\n",(int)V.size(),(int)F.size());
    if(outPrefixed){ for(auto&p:V) printf("v %.15g %.15g %.15g\n",p.x,p.y,p.z); for(auto&f:F) printf("f %d %d %d\n",f.v[0]+1,f.v[1]+1,f.v[2]+1); }
    else{ for(auto&p:V) printf("%.15g %.15g %.15g\n",p.x,p.y,p.z); for(auto&f:F) printf("%d %d %d\n",f.v[0]+1,f.v[1]+1,f.v[2]+1); }
}

int main(){
    T0=chrono::steady_clock::now();
    if(!read_input()) return 0;
    Candidate original; original.name="original"; original.V=Orig; original.F=OrigF;
    vector<Candidate> C;
    generate_specials(C);
    if(elapsed()<12.5){
        vector<int> cps; int n=N0; double rs[]={0.20,0.16,0.125,0.10,0.08}; for(double r:rs){int t=max(4,(int)ceil(n*r)); if(t<n)cps.push_back(t);} sort(cps.begin(),cps.end()); cps.erase(unique(cps.begin(),cps.end()),cps.end());
        QemSolver qs(N0<5000?0.38:0.48,17.0); qs.checkpoints=cps; vector<Candidate> q=qs.run();
        for(auto &c:q) add_candidate(C,move(c),"qem",true,true);
    }
    if(C.empty()){ print_mesh(original.V,original.F); return 0; }
    int res = N0<30000 ? 176 : (N0<100000 ? 144 : 128);
    if(elapsed()>14.5) res=112;
    Proxy proxy(res);
    int best=-1; double bestRank=-1e100; double bestProxy=-1;
    for(int i=0;i<(int)C.size();i++){
        if(elapsed()>18.5) break;
        C[i].proxy=proxy.score(C[i].V,C[i].F);
        double ratio=(double)C[i].V.size()/max(1,N0);
        double rank;
        if(C[i].proxy>=0.908) rank=1000.0-ratio*100.0 + (C[i].proxy-0.908)*20.0; // among safe candidates, prefer compression.
        else if(C[i].proxy>=0.892) rank=100.0 + C[i].proxy*20.0 - ratio*10.0; // proxy can be pessimistic on silhouettes.
        else rank=C[i].proxy*10.0-ratio;
        if(rank>bestRank){bestRank=rank;best=i;bestProxy=C[i].proxy;}
    }
    if(best>=0 && (bestProxy>=0.892 || C[best].fromQem)) print_mesh(C[best].V,C[best].F);
    else print_mesh(original.V,original.F);
    return 0;
}

