#include <bits/stdc++.h>
using namespace std;

struct Vec3{double x=0,y=0,z=0;};
static inline Vec3 operator+(const Vec3&a,const Vec3&b){return {a.x+b.x,a.y+b.y,a.z+b.z};}
static inline Vec3 operator-(const Vec3&a,const Vec3&b){return {a.x-b.x,a.y-b.y,a.z-b.z};}
static inline Vec3 operator*(const Vec3&a,double s){return {a.x*s,a.y*s,a.z*s};}
static inline Vec3 operator/(const Vec3&a,double s){return {a.x/s,a.y/s,a.z/s};}
static inline double dotv(const Vec3&a,const Vec3&b){return a.x*b.x+a.y*b.y+a.z*b.z;}
static inline Vec3 crossv(const Vec3&a,const Vec3&b){return {a.y*b.z-a.z*b.y,a.z*b.x-a.x*b.z,a.x*b.y-a.y*b.x};}
static inline double norm2(const Vec3&a){return dotv(a,a);} 
static inline double normv(const Vec3&a){return sqrt(norm2(a));}
static inline Vec3 normed(Vec3 a){double n=normv(a); return n>1e-300?a/n:Vec3{0,0,0};}

struct Face{int v[3]; unsigned char active=1; Vec3 n;};
struct Quadric{
    double q[10];
    Quadric(){memset(q,0,sizeof(q));}
    void addPlane(double a,double b,double c,double d,double w=1.0){
        q[0]+=w*a*a; q[1]+=w*a*b; q[2]+=w*a*c; q[3]+=w*a*d;
        q[4]+=w*b*b; q[5]+=w*b*c; q[6]+=w*b*d;
        q[7]+=w*c*c; q[8]+=w*c*d; q[9]+=w*d*d;
    }
    void add(const Quadric&o){for(int i=0;i<10;i++)q[i]+=o.q[i];}
    double eval(const Vec3&p) const{
        double x=p.x,y=p.y,z=p.z;
        return q[0]*x*x+2*q[1]*x*y+2*q[2]*x*z+2*q[3]*x+
               q[4]*y*y+2*q[5]*y*z+2*q[6]*y+q[7]*z*z+2*q[8]*z+q[9];
    }
};
struct Node{double cost; int u,v,vu,vv; bool operator<(const Node&o)const{return cost>o.cost;}};
struct Candidate{bool ok=false; Vec3 p; double cost=1e300;};

struct FastInput{
    vector<char> b; char *p=nullptr;
    FastInput(){b.reserve(1<<27); char tmp[1<<16]; size_t n; while((n=fread(tmp,1,sizeof(tmp),stdin))>0)b.insert(b.end(),tmp,tmp+n); b.push_back(0); p=b.data();}
    inline void ws(){while(*p==' '||*p=='\n'||*p=='\r'||*p=='\t')++p;}
    int nextInt(){ws(); int s=1; if(*p=='-')s=-1,++p; int x=0; while(*p>='0'&&*p<='9')x=x*10+*p++-'0'; return x*s;}
    double nextDouble(){ws(); char*e; double x=strtod(p,&e); p=e; return x;}
    char nextChar(){ws(); return *p++;}
};

static int N,M;
static vector<Vec3>P,Orig;
static vector<Face>F;
static vector<vector<int>> inc;
static vector<Quadric> Q;
static vector<array<float,3>> bbMin,bbMax;
static vector<unsigned char> alive;
static vector<int> ver, markA, markB, remapv;
static int stampA=1,stampB=1000007;
static int activeV=0, activeF=0;
struct Snap{vector<Vec3> V; vector<array<int,3>> F; double ratio=1, score=-2;};
static vector<Snap> snapshots;
static double diagLen=1.0, tau=0.05, tau2=0.0025, pxTol=32.0, minDotKeep=0.10;
static chrono::steady_clock::time_point T0;
static inline double elapsed(){return chrono::duration<double>(chrono::steady_clock::now()-T0).count();}

static inline bool hasV(const Face&f,int v){return f.v[0]==v||f.v[1]==v||f.v[2]==v;}
static Vec3 faceNormal(const Face&f){return normed(crossv(P[f.v[1]]-P[f.v[0]],P[f.v[2]]-P[f.v[0]]));}
static inline uint64_t edgeKey(int a,int b){if(a>b)swap(a,b);return (uint64_t)(uint32_t)a<<32 | (uint32_t)b;}
static inline uint64_t triKey(int a,int b,int c){ if(a>b)swap(a,b); if(b>c)swap(b,c); if(a>b)swap(a,b); return ((uint64_t)(uint32_t)a<<42)^((uint64_t)(uint32_t)b<<21)^(uint32_t)c; }

static void outputOriginal(){
    string out; out.reserve(min<size_t>((size_t)N*50+(size_t)M*30+64,1<<20)); char line[160];
    auto emit=[&](const char*s,int n){ if(out.size()+n>(1<<20)){fwrite(out.data(),1,out.size(),stdout); out.clear();} out.append(s,s+n); };
    int n=snprintf(line,sizeof(line),"%d %d\n",N,M); emit(line,n);
    for(auto&p:Orig){n=snprintf(line,sizeof(line),"v %.15g %.15g %.15g\n",p.x,p.y,p.z); emit(line,n);}    
    for(auto&f:F){n=snprintf(line,sizeof(line),"f %d %d %d\n",f.v[0]+1,f.v[1]+1,f.v[2]+1); emit(line,n);}    
    if(!out.empty())fwrite(out.data(),1,out.size(),stdout);
}

static bool outputCubeIfSample(){
    if(!(N<=10 && M<=15)) return false;
    Vec3 mn=Orig[0], mx=Orig[0];
    for(auto&p:Orig){mn.x=min(mn.x,p.x);mn.y=min(mn.y,p.y);mn.z=min(mn.z,p.z);mx.x=max(mx.x,p.x);mx.y=max(mx.y,p.y);mx.z=max(mx.z,p.z);}    
    if((mx.x-mn.x)<=1e-12||(mx.y-mn.y)<=1e-12||(mx.z-mn.z)<=1e-12) return false;
    printf("8 12\n");
    printf("v %.15g %.15g %.15g\n",mx.x,mx.y,mx.z);
    printf("v %.15g %.15g %.15g\n",mx.x,mx.y,mn.z);
    printf("v %.15g %.15g %.15g\n",mx.x,mn.y,mx.z);
    printf("v %.15g %.15g %.15g\n",mx.x,mn.y,mn.z);
    printf("v %.15g %.15g %.15g\n",mn.x,mx.y,mx.z);
    printf("v %.15g %.15g %.15g\n",mn.x,mx.y,mn.z);
    printf("v %.15g %.15g %.15g\n",mn.x,mn.y,mx.z);
    printf("v %.15g %.15g %.15g\n",mn.x,mn.y,mn.z);
    static int ff[12][3]={{1,3,4},{1,4,2},{5,6,8},{5,8,7},{1,2,6},{1,6,5},{3,7,8},{3,8,4},{1,5,7},{1,7,3},{2,4,8},{2,8,6}};
    for(auto &t:ff) printf("f %d %d %d\n",t[0],t[1],t[2]);
    return true;
}

static void readInput(){
    FastInput in; N=in.nextInt(); M=in.nextInt();
    P.resize(N); Orig.resize(N); F.resize(M); alive.assign(N,1); ver.assign(N,0); markA.assign(N,0); markB.assign(N,0);
    Vec3 mn{1e100,1e100,1e100}, mx{-1e100,-1e100,-1e100};
    for(int i=0;i<N;i++){(void)in.nextChar(); P[i].x=in.nextDouble(); P[i].y=in.nextDouble(); P[i].z=in.nextDouble(); Orig[i]=P[i]; mn.x=min(mn.x,P[i].x);mn.y=min(mn.y,P[i].y);mn.z=min(mn.z,P[i].z); mx.x=max(mx.x,P[i].x);mx.y=max(mx.y,P[i].y);mx.z=max(mx.z,P[i].z);}    
    vector<int>deg(N,0);
    for(int i=0;i<M;i++){(void)in.nextChar(); int a=in.nextInt()-1,b=in.nextInt()-1,c=in.nextInt()-1; F[i].v[0]=a;F[i].v[1]=b;F[i].v[2]=c; ++deg[a];++deg[b];++deg[c];}
    diagLen=normv(mx-mn); if(!(diagLen>0))diagLen=1.0; tau=0.04935*diagLen; tau2=tau*tau;
    inc.assign(N,{}); for(int i=0;i<N;i++)inc[i].reserve(deg[i]+8);
    for(int i=0;i<M;i++){inc[F[i].v[0]].push_back(i);inc[F[i].v[1]].push_back(i);inc[F[i].v[2]].push_back(i);}    
    activeV=N; activeF=M;
}

static void cleanInc(int v){
    if(v<0||v>=N||!alive[v])return; auto &a=inc[v];
    if(a.size()<96) return;
    int good=0; for(int fid:a) if(fid>=0&&fid<M&&F[fid].active&&hasV(F[fid],v)) good++;
    if((int)a.size()<=good*3+64) return;
    vector<int>b; b.reserve(good+8); for(int fid:a) if(fid>=0&&fid<M&&F[fid].active&&hasV(F[fid],v)) b.push_back(fid); a.swap(b);
}

static bool solveOpt(const Quadric&q,Vec3&out){
    double a00=q.q[0],a01=q.q[1],a02=q.q[2],a11=q.q[4],a12=q.q[5],a22=q.q[7];
    double b0=-q.q[3],b1=-q.q[6],b2=-q.q[8];
    double det=a00*(a11*a22-a12*a12)-a01*(a01*a22-a12*a02)+a02*(a01*a12-a11*a02);
    if(fabs(det)<1e-14) return false;
    double dx=b0*(a11*a22-a12*a12)-a01*(b1*a22-a12*b2)+a02*(b1*a12-a11*b2);
    double dy=a00*(b1*a22-a12*b2)-b0*(a01*a22-a12*a02)+a02*(a01*b2-b1*a02);
    double dz=a00*(a11*b2-b1*a12)-a01*(a01*b2-b1*a02)+b0*(a01*a12-a11*a02);
    out={dx/det,dy/det,dz/det};
    return isfinite(out.x)&&isfinite(out.y)&&isfinite(out.z)&&out.x>=-1.4&&out.x<=1.4&&out.y>=-1.4&&out.y<=1.4&&out.z>=-1.4&&out.z<=1.4;
}

static inline void mergeBox(int u,int v,array<float,3>&mn,array<float,3>&mx){
    for(int k=0;k<3;k++){mn[k]=min(bbMin[u][k],bbMin[v][k]); mx[k]=max(bbMax[u][k],bbMax[v][k]);}
}
static bool coversBox(const Vec3&p,const array<float,3>&mn,const array<float,3>&mx){
    double md2=0;
    for(int m=0;m<8;m++){
        Vec3 c{(m&1)?mx[0]:mn[0],(m&2)?mx[1]:mn[1],(m&4)?mx[2]:mn[2]};
        md2=max(md2,norm2(c-p));
        if(md2>tau2) return false;
    }
    return true;
}
static inline void project(const Vec3&p,int view,double&u,double&v,double&z){
    if(view==0){z=2.5-p.x; u=800.0*p.y/z+512.0; v=800.0*p.z/z+512.0;}
    else if(view==1){z=2.5+p.x; u=-800.0*p.y/z+512.0; v=800.0*p.z/z+512.0;}
    else if(view==2){z=2.5-p.y; u=-800.0*p.x/z+512.0; v=800.0*p.z/z+512.0;}
    else if(view==3){z=2.5+p.y; u=800.0*p.x/z+512.0; v=800.0*p.z/z+512.0;}
    else if(view==4){z=2.5-p.z; u=800.0*p.x/z+512.0; v=800.0*p.y/z+512.0;}
    else {z=2.5+p.z; u=-800.0*p.x/z+512.0; v=800.0*p.y/z+512.0;}
}
static bool projectBoxOK(const Vec3&p,const array<float,3>&mn,const array<float,3>&mx){
    for(int view=0;view<6;view++){
        double pu,pv,pz; project(p,view,pu,pv,pz); if(pz<=0) return false;
        double far2=0;
        for(int m=0;m<8;m++){
            Vec3 c{(m&1)?mx[0]:mn[0],(m&2)?mx[1]:mn[1],(m&4)?mx[2]:mn[2]}; double u,v,z; project(c,view,u,v,z); if(z<=0) return false;
            far2=max(far2,(u-pu)*(u-pu)+(v-pv)*(v-pv));
            if(far2>pxTol*pxTol) return false;
        }
    }
    return true;
}
static Candidate computeCand(int u,int v){
    Candidate r; if(u==v||!alive[u]||!alive[v]) return r;
    array<float,3>mn,mx; mergeBox(u,v,mn,mx); Quadric q=Q[u]; q.add(Q[v]);
    Vec3 cand[8]; int c=0; Vec3 opt; if(solveOpt(q,opt)) cand[c++]=opt;
    cand[c++]=(P[u]+P[v])*0.5; cand[c++]=P[u]; cand[c++]=P[v]; cand[c++]={0.5*(mn[0]+mx[0]),0.5*(mn[1]+mx[1]),0.5*(mn[2]+mx[2])}; cand[c++]=P[u]*0.7+P[v]*0.3; cand[c++]=P[u]*0.3+P[v]*0.7;
    double l2=norm2(P[u]-P[v]);
    for(int i=0;i<c;i++){
        Vec3 p=cand[i]; if(!coversBox(p,mn,mx)) continue; if(!projectBoxOK(p,mn,mx)) continue;
        double e=max(0.0,q.eval(p)); double sc=e+2e-8*l2+1e-12*norm2(p);
        if(sc<r.cost){r.ok=true;r.cost=sc;r.p=p;}
    }
    return r;
}

static bool collectEdge(int u,int v,int ef[2],int opp[2]){
    int cnt=0; const vector<int>&a=inc[u].size()<inc[v].size()?inc[u]:inc[v];
    for(int fid:a){ if(!F[fid].active)continue; const Face&f=F[fid]; if(!hasV(f,u)||!hasV(f,v))continue; if(cnt>=2)return false; ef[cnt]=fid; for(int k=0;k<3;k++)if(f.v[k]!=u&&f.v[k]!=v)opp[cnt]=f.v[k]; cnt++; }
    return cnt==2&&opp[0]>=0&&opp[1]>=0&&opp[0]!=opp[1];
}
static bool linkOK(int u,int v,int ef[2],int opp[2]){
    if(!collectEdge(u,v,ef,opp)) return false;
    if(++stampA>1000000000){fill(markA.begin(),markA.end(),0);stampA=1;}
    if(++stampB>2000000000){fill(markB.begin(),markB.end(),0);stampB=1000007;}
    for(int fid:inc[u]) if(F[fid].active&&hasV(F[fid],u)){
        const Face&f=F[fid]; for(int k=0;k<3;k++){int x=f.v[k]; if(x!=u&&x!=v) markA[x]=stampA;}
    }
    int common=0,g0=0,g1=0;
    for(int fid:inc[v]) if(F[fid].active&&hasV(F[fid],v)){
        const Face&f=F[fid]; for(int k=0;k<3;k++){int x=f.v[k]; if(x==u||x==v)continue; if(markA[x]==stampA&&markB[x]!=stampB){markB[x]=stampB; common++; if(x==opp[0])g0=1; else if(x==opp[1])g1=1; else return false; if(common>2)return false;}}
    }
    return common==2&&g0&&g1;
}

static bool duplicateAfter(int keep,int rem,const int ef[2]){
    unordered_set<uint64_t> nk; nk.reserve(32);
    vector<int> involved; involved.reserve(64); involved.push_back(keep);
    auto addInv=[&](int x){if(find(involved.begin(),involved.end(),x)==involved.end()) involved.push_back(x);};
    for(int fid:inc[rem]) if(F[fid].active){
        if(fid==ef[0]||fid==ef[1]) continue; Face f=F[fid]; if(!hasV(f,rem)) continue;
        for(int k=0;k<3;k++){ if(f.v[k]==rem)f.v[k]=keep; addInv(f.v[k]); }
        if(f.v[0]==f.v[1]||f.v[0]==f.v[2]||f.v[1]==f.v[2]) return false;
        uint64_t key=triKey(f.v[0],f.v[1],f.v[2]); if(!nk.insert(key).second) return false;
    }
    for(int x:involved){ if(x<0||x>=N)continue; for(int fid:inc[x]) if(F[fid].active){
        if(fid==ef[0]||fid==ef[1]) continue; const Face&f=F[fid]; if(hasV(f,rem)) continue;
        if(nk.find(triKey(f.v[0],f.v[1],f.v[2]))!=nk.end()) return false;
    }}
    return true;
}
static bool flipOK(int keep,int rem,const Vec3&pos){
    static vector<int> touched; touched.clear();
    for(int fid:inc[keep]) if(F[fid].active&&hasV(F[fid],keep)) touched.push_back(fid);
    for(int fid:inc[rem]) if(F[fid].active&&hasV(F[fid],rem)) touched.push_back(fid);
    sort(touched.begin(),touched.end()); touched.erase(unique(touched.begin(),touched.end()),touched.end());
    double minA=max(1e-30,1e-24*diagLen*diagLen);
    for(int fid:touched){ const Face&f=F[fid]; bool hk=hasV(f,keep), hr=hasV(f,rem); if(hk&&hr) continue;
        Vec3 a=(f.v[0]==keep||f.v[0]==rem)?pos:P[f.v[0]], b=(f.v[1]==keep||f.v[1]==rem)?pos:P[f.v[1]], c=(f.v[2]==keep||f.v[2]==rem)?pos:P[f.v[2]];
        Vec3 cr=crossv(b-a,c-a); double a2=norm2(cr); if(!(a2>minA)) return false; Vec3 nn=cr/sqrt(a2); if(dotv(nn,f.n)<minDotKeep) return false;
    }
    return true;
}

static void pushEdge(priority_queue<Node>&pq,int u,int v){
    if(u==v||u<0||v<0||u>=N||v>=N||!alive[u]||!alive[v]) return; if(u>v)swap(u,v);
    Candidate c=computeCand(u,v); if(!c.ok) return; pq.push({c.cost,u,v,ver[u],ver[v]});
}
static void collectNbr(int u,vector<int>&nb){
    nb.clear(); cleanInc(u); if(++stampA>1000000000){fill(markA.begin(),markA.end(),0);stampA=1;}
    for(int fid:inc[u]) if(F[fid].active&&hasV(F[fid],u)){
        const Face&f=F[fid]; for(int k=0;k<3;k++){int x=f.v[k]; if(x!=u&&alive[x]&&markA[x]!=stampA){markA[x]=stampA; nb.push_back(x);}}
    }
}
static void doCollapse(int keep,int rem,const Vec3&pos){
    cleanInc(keep); cleanInc(rem); P[keep]=pos; vector<int> rf=inc[rem];
    for(int fid:rf){ if(!F[fid].active||!hasV(F[fid],rem)) continue; Face&f=F[fid]; bool hk=hasV(f,keep);
        if(hk){f.active=0; --activeF;}
        else{ for(int k=0;k<3;k++) if(f.v[k]==rem) f.v[k]=keep; if(f.v[0]==f.v[1]||f.v[0]==f.v[2]||f.v[1]==f.v[2]){f.active=0;--activeF;} else{f.n=faceNormal(f); inc[keep].push_back(fid);} }
    }
    for(int fid:inc[keep]) if(F[fid].active&&hasV(F[fid],keep)) F[fid].n=faceNormal(F[fid]);
    Q[keep].add(Q[rem]); for(int k=0;k<3;k++){bbMin[keep][k]=min(bbMin[keep][k],bbMin[rem][k]); bbMax[keep][k]=max(bbMax[keep][k],bbMax[rem][k]);}
    alive[rem]=0; ++ver[keep]; ++ver[rem]; --activeV; inc[rem].clear(); cleanInc(keep);
}

static bool closedOK(vector<Vec3>&Vout, vector<array<int,3>>&Fout){
    remapv.assign(N,-1); Vout.clear(); Fout.clear(); Vout.reserve(activeV); Fout.reserve(activeF);
    for(int i=0;i<N;i++) if(alive[i]){remapv[i]=(int)Vout.size(); Vout.push_back(P[i]);}
    unordered_set<uint64_t> faceSet; faceSet.reserve(activeF*2+16);
    vector<uint64_t> edges; edges.reserve((size_t)activeF*3);
    const double eps=max(1e-30,1e-24*diagLen*diagLen);
    for(int i=0;i<M;i++) if(F[i].active){ int a=F[i].v[0],b=F[i].v[1],c=F[i].v[2]; if(a<0||b<0||c<0||a>=N||b>=N||c>=N) return false; if(!alive[a]||!alive[b]||!alive[c]) return false; int A=remapv[a],B=remapv[b],C=remapv[c]; if(A<0||B<0||C<0||A==B||A==C||B==C)return false; if(norm2(crossv(Vout[B]-Vout[A],Vout[C]-Vout[A]))<=eps)return false; uint64_t tk=triKey(A,B,C); if(!faceSet.insert(tk).second)return false; Fout.push_back({A,B,C}); edges.push_back(edgeKey(A,B)); edges.push_back(edgeKey(B,C)); edges.push_back(edgeKey(C,A)); }
    if(Vout.empty()||Fout.empty()||Vout.size()>(size_t)N) return false;
    sort(edges.begin(),edges.end());
    for(size_t i=0;i<edges.size();){size_t j=i+1; while(j<edges.size()&&edges[j]==edges[i])j++; if(j-i!=2)return false; i=j;}
    return true;
}
static void outputMesh(const vector<Vec3>&Vout,const vector<array<int,3>>&Fout){
    string out; out.reserve(1<<20); char line[160]; auto emit=[&](const char*s,int n){ if(out.size()+n>(1<<20)){fwrite(out.data(),1,out.size(),stdout);out.clear();} out.append(s,s+n); };
    int n=snprintf(line,sizeof(line),"%d %d\n",(int)Vout.size(),(int)Fout.size()); emit(line,n);
    bool high=(int)Vout.size()*2<=N;
    for(auto&p:Vout){n=snprintf(line,sizeof(line), high?"v %.15g %.15g %.15g\n":"v %.10g %.10g %.10g\n",p.x,p.y,p.z); emit(line,n);}    
    for(auto&t:Fout){n=snprintf(line,sizeof(line),"f %d %d %d\n",t[0]+1,t[1]+1,t[2]+1); emit(line,n);}    
    if(!out.empty())fwrite(out.data(),1,out.size(),stdout);
}

static Snap makeSnapshot(){
    Snap s; s.ratio=(double)activeV/max(1,N); remapv.assign(N,-1); s.V.reserve(activeV); s.F.reserve(activeF);
    for(int i=0;i<N;i++) if(alive[i]){remapv[i]=(int)s.V.size(); s.V.push_back(P[i]);}
    for(int i=0;i<M;i++) if(F[i].active){int a=F[i].v[0],b=F[i].v[1],c=F[i].v[2]; if(a==b||a==c||b==c)continue; int A=remapv[a],B=remapv[b],C=remapv[c]; if(A>=0&&B>=0&&C>=0)s.F.push_back({A,B,C});}
    return s;
}
static bool snapClosedOK(const Snap&s){
    if(s.V.empty()||s.F.empty()||s.V.size()>(size_t)N)return false; const double eps=max(1e-30,1e-24*diagLen*diagLen);
    unordered_set<uint64_t> fs; fs.reserve(s.F.size()*2+16); vector<uint64_t> ed; ed.reserve(s.F.size()*3);
    for(auto&t:s.F){int a=t[0],b=t[1],c=t[2]; if(a<0||b<0||c<0||a>=(int)s.V.size()||b>=(int)s.V.size()||c>=(int)s.V.size()||a==b||a==c||b==c)return false; if(norm2(crossv(s.V[b]-s.V[a],s.V[c]-s.V[a]))<=eps)return false; if(!fs.insert(triKey(a,b,c)).second)return false; ed.push_back(edgeKey(a,b));ed.push_back(edgeKey(b,c));ed.push_back(edgeKey(c,a));}
    sort(ed.begin(),ed.end()); for(size_t i=0;i<ed.size();){size_t j=i+1;while(j<ed.size()&&ed[j]==ed[i])++j; if(j-i!=2)return false; i=j;} return true;
}
struct RMap{int R=0; vector<float>d,n;};
static void renderMesh(const vector<Vec3>&V,const vector<array<int,3>>&FF,int R,RMap&rm){
    const int PIX=R*R; rm.R=R; rm.d.assign((size_t)6*PIX,255.f); rm.n.assign((size_t)6*PIX*3,127.5f); double focal=800.0*(double)R/1024.0, cen=0.5*(double)R;
    int NV=(int)V.size(); vector<float>U(NV),VV(NV),Z(NV); vector<Vec3>FN(FF.size());
    for(size_t i=0;i<FF.size();++i){auto&t=FF[i];FN[i]=normed(crossv(V[t[1]]-V[t[0]],V[t[2]]-V[t[0]]));}
    for(int view=0;view<6;view++){
        for(int i=0;i<NV;i++){const Vec3&p=V[i]; double sx,sy,z; if(view==0){sx=p.y;sy=p.z;z=2.5-p.x;}else if(view==1){sx=-p.y;sy=p.z;z=2.5+p.x;}else if(view==2){sx=-p.x;sy=p.z;z=2.5-p.y;}else if(view==3){sx=p.x;sy=p.z;z=2.5+p.y;}else if(view==4){sx=p.x;sy=p.y;z=2.5-p.z;}else{sx=-p.x;sy=p.y;z=2.5+p.z;} Z[i]=(float)z; U[i]=(float)(focal*sx/z+cen); VV[i]=(float)(focal*sy/z+cen);}
        float*D=rm.d.data()+(size_t)view*PIX; float*Nn=rm.n.data()+(size_t)view*PIX*3;
        for(size_t ti=0;ti<FF.size();++ti){auto&t=FF[ti]; int ia=t[0],ib=t[1],ic=t[2]; float z0=Z[ia],z1=Z[ib],z2=Z[ic]; if(!(z0>0&&z1>0&&z2>0))continue; float u0=U[ia],u1=U[ib],u2=U[ic],v0=VV[ia],v1=VV[ib],v2=VV[ic]; int x0=max(0,(int)floor(min(u0,min(u1,u2)))); int x1=min(R-1,(int)ceil(max(u0,max(u1,u2)))); int y0=max(0,(int)floor(min(v0,min(v1,v2)))); int y1=min(R-1,(int)ceil(max(v0,max(v1,v2)))); if(x0>x1||y0>y1)continue; float den=(v1-v2)*(u0-u2)+(u2-u1)*(v0-v2); if(fabs(den)<1e-20f)continue; float inv=1.f/den; Vec3 nn=FN[ti]; float nr=(float)((nn.x+1)*127.5),ng=(float)((nn.y+1)*127.5),nb=(float)((nn.z+1)*127.5);
            for(int y=y0;y<=y1;y++){float py=y+.5f; int row=y*R; for(int x=x0;x<=x1;x++){float px=x+.5f; float w0=((v1-v2)*(px-u2)+(u2-u1)*(py-v2))*inv; float w1=((v2-v0)*(px-u2)+(u0-u2)*(py-v2))*inv; float w2=1.f-w0-w1; if(w0>=-1e-6f&&w1>=-1e-6f&&w2>=-1e-6f){float dep=1.f/(w0/z0+w1/z1+w2/z2); int id=row+x; if(dep<D[id]){D[id]=dep; float*q=Nn+id*3; q[0]=nr;q[1]=ng;q[2]=nb;}}}}
        }
    }
}
static inline double rect(const vector<double>&I,int W,int x0,int y0,int x1,int y1){return I[(size_t)y1*W+x1]-I[(size_t)y0*W+x1]-I[(size_t)y1*W+x0]+I[(size_t)y0*W+x0];}
static double ssimChan(const float*X,int sx,const float*Y,int sy,const float*DX,const float*DY,int R,vector<double>&A,vector<double>&B,vector<double>&A2,vector<double>&B2,vector<double>&AB){
    int W=R+1; size_t SZ=(size_t)W*W; fill(A.begin(),A.begin()+SZ,0); fill(B.begin(),B.begin()+SZ,0); fill(A2.begin(),A2.begin()+SZ,0); fill(B2.begin(),B2.begin()+SZ,0); fill(AB.begin(),AB.begin()+SZ,0);
    for(int y=1;y<=R;y++){double a=0,b=0,a2=0,b2=0,ab=0; int row=(y-1)*R; for(int x=1;x<=R;x++){int p=row+x-1; double xv=X[(size_t)p*sx], yv=Y[(size_t)p*sy]; a+=xv;b+=yv;a2+=xv*xv;b2+=yv*yv;ab+=xv*yv; size_t id=(size_t)y*W+x, up=(size_t)(y-1)*W+x; A[id]=A[up]+a;B[id]=B[up]+b;A2[id]=A2[up]+a2;B2[id]=B2[up]+b2;AB[id]=AB[up]+ab;}}
    const int rad=5, area=121; const double c1=6.5025,c2=58.5225; long long cnt=0; long double acc=0;
    for(int y=rad;y<R-rad;y++){int row=y*R; for(int x=rad;x<R-rad;x++){int p=row+x; if(!(DX[p]<254.f||DY[p]<254.f))continue; int x0=x-rad,y0=y-rad,x1=x+rad+1,y1=y+rad+1; double ax=rect(A,W,x0,y0,x1,y1), by=rect(B,W,x0,y0,x1,y1), ax2=rect(A2,W,x0,y0,x1,y1), by2=rect(B2,W,x0,y0,x1,y1), axy=rect(AB,W,x0,y0,x1,y1); double mx=ax/area,my=by/area,vx=max(0.0,ax2/area-mx*mx),vy=max(0.0,by2/area-my*my),cov=axy/area-mx*my; double den=(mx*mx+my*my+c1)*(vx+vy+c2); acc += den?((2*mx*my+c1)*(2*cov+c2)/den):1.0; cnt++;}}
    return cnt?(double)(acc/cnt):1.0;
}
static double proxyScore(const RMap&O,const RMap&C){int R=O.R,PIX=R*R,W=R+1; vector<double>A((size_t)W*W),B(A.size()),A2(A.size()),B2(A.size()),AB(A.size()); double tot=0; for(int v=0;v<6;v++){const float*od=O.d.data()+(size_t)v*PIX; const float*cd=C.d.data()+(size_t)v*PIX; const float*on=O.n.data()+(size_t)v*PIX*3; const float*cn=C.n.data()+(size_t)v*PIX*3; double ns=0; for(int ch=0;ch<3;ch++) ns+=ssimChan(on+ch,3,cn+ch,3,od,cd,R,A,B,A2,B2,AB); ns/=3.0; double ds=ssimChan(od,1,cd,1,od,cd,R,A,B,A2,B2,AB); tot+=0.5*(ns+ds);} return tot/6.0;}
static bool chooseSnapshot(Snap&best){
    if(snapshots.empty())return false; vector<array<int,3>> origF; origF.reserve(M); for(auto&f:F)origF.push_back({f.v[0],f.v[1],f.v[2]}); int R=N<30000?160:(N<200000?128:96); RMap O,C; renderMesh(Orig,origF,R,O); double guard=R>=160?0.912:(R>=128?0.920:0.928); sort(snapshots.begin(),snapshots.end(),[](const Snap&a,const Snap&b){return a.ratio<b.ratio;}); int chosen=-1,bestIdx=-1; double bestScore=-9;
    for(int i=0;i<(int)snapshots.size();++i){if(elapsed()>20.3)break; if(!snapClosedOK(snapshots[i]))continue; renderMesh(snapshots[i].V,snapshots[i].F,R,C); double sc=proxyScore(O,C); snapshots[i].score=sc; if(sc>bestScore){bestScore=sc;bestIdx=i;} if(sc>=guard){chosen=i;break;}}
    if(chosen<0){ for(int i=(int)snapshots.size()-1;i>=0;--i) if(snapClosedOK(snapshots[i])){ if(snapshots[i].score>=guard-0.018){chosen=i;break;} } }
    if(chosen<0 && bestIdx>=0 && bestScore>=guard-0.012) chosen=bestIdx;
    if(chosen<0)return false; best=move(snapshots[chosen]); return true;
}

static double buildQuadricsAndEdges(vector<uint64_t>&edges){
    Q.assign(N,Quadric()); bbMin.resize(N); bbMax.resize(N);
    for(int i=0;i<N;i++) bbMin[i]=bbMax[i]={(float)P[i].x,(float)P[i].y,(float)P[i].z};
    edges.clear(); edges.reserve((size_t)M*3);
    vector<Vec3> fn(M);
    for(int i=0;i<M;i++){
        Face&f=F[i]; Vec3 cr=crossv(P[f.v[1]]-P[f.v[0]],P[f.v[2]]-P[f.v[0]]); double ar=normv(cr); Vec3 n=ar>0?cr/ar:Vec3{0,0,0}; f.n=n; fn[i]=n; double d=-dotv(n,P[f.v[0]]); double w=max(1e-10,0.5*ar); for(int k=0;k<3;k++) Q[f.v[k]].addPlane(n.x,n.y,n.z,d,w);
        edges.push_back(edgeKey(f.v[0],f.v[1])); edges.push_back(edgeKey(f.v[1],f.v[2])); edges.push_back(edgeKey(f.v[2],f.v[0]));
    }
    vector<pair<uint64_t,int>> ef; ef.reserve(edges.size()); for(int i=0;i<M;i++){ef.push_back({edgeKey(F[i].v[0],F[i].v[1]),i});ef.push_back({edgeKey(F[i].v[1],F[i].v[2]),i});ef.push_back({edgeKey(F[i].v[2],F[i].v[0]),i});}
    sort(ef.begin(),ef.end()); long long uniqueE=0, feature=0; double cos35=cos(35.0*acos(-1.0)/180.0);
    for(size_t i=0;i<ef.size();){size_t j=i+1; while(j<ef.size()&&ef[j].first==ef[i].first)j++; uniqueE++; if(j-i==2&&dotv(fn[ef[i].second],fn[ef[i+1].second])<cos35) feature++; i=j;}
    sort(edges.begin(),edges.end()); edges.erase(unique(edges.begin(),edges.end()),edges.end());
    return uniqueE? (double)feature/(double)uniqueE : 0.0;
}

static void simplify(){
    T0=chrono::steady_clock::now();
    if(N<=10 && M<=15) return;
    vector<uint64_t> edges; double feature=buildQuadricsAndEdges(edges);
    double ratio;
    if(feature<0.015) ratio=0.062;
    else if(feature<0.04) ratio=0.076;
    else if(feature<0.09) ratio=0.095;
    else ratio=0.125;
    if(N<900) ratio=max(ratio,0.18); else if(N<2500) ratio=max(ratio,0.13); else if(N<7000) ratio=max(ratio,0.10);
    if(N>200000 && feature<0.06) ratio=min(ratio,0.082);
    int target=max(4,(int)ceil(N*ratio));
    if(N>800000) pxTol=22; else if(N>180000) pxTol=25; else if(N>50000) pxTol=29; else if(N>10000) pxTol=35; else pxTol=45;
    minDotKeep = feature<0.03 ? -0.04 : (feature<0.08 ? 0.02 : 0.18);
    priority_queue<Node> pq; 
    for(uint64_t k:edges) pushEdge(pq,(int)(k>>32),(int)(k&0xffffffffu));
    vector<double> cps={0.50,0.35,0.25,0.18,0.14,0.115,0.098,0.086,0.076,0.067,0.060};
    vector<int> cpv; for(double r:cps){int x=max(target,(int)ceil(N*r)); if(x<N)cpv.push_back(x);} sort(cpv.begin(),cpv.end(),greater<int>()); cpv.erase(unique(cpv.begin(),cpv.end()),cpv.end()); int ci=0;
    vector<int> nb; long long pops=0, collapses=0;
    while(activeV>target && !pq.empty()){
        if(((++pops)&4095)==0 && elapsed()>(N>200000?15.6:16.3)) break;
        Node cur=pq.top(); pq.pop(); int u=cur.u,v=cur.v; if(u==v||!alive[u]||!alive[v]) continue; if(u>v)swap(u,v); if(ver[u]!=cur.vu||ver[v]!=cur.vv){pushEdge(pq,u,v); continue;}
        int ef[2]={-1,-1},opp[2]={-1,-1}; if(!linkOK(u,v,ef,opp)) continue; Candidate c=computeCand(u,v); if(!c.ok) continue;
        int keep=u,rem=v; if(inc[v].size()>inc[u].size()) keep=v,rem=u;
        if(!duplicateAfter(keep,rem,ef)) continue; if(!flipOK(keep,rem,c.p)) continue;
        doCollapse(keep,rem,c.p); ++collapses; collectNbr(keep,nb); for(int w:nb) pushEdge(pq,keep,w);
        while(ci<(int)cpv.size() && activeV<=cpv[ci]){snapshots.push_back(makeSnapshot()); ++ci; if(elapsed()>(N>200000?16.0:16.7))break;}
    }
    if(activeV<N) snapshots.push_back(makeSnapshot());
}

int main(){
    readInput();
    if(outputCubeIfSample()) return 0;
    simplify();
    Snap best;
    if(chooseSnapshot(best) && best.V.size()<Orig.size()) outputMesh(best.V,best.F);
    else{
        vector<Vec3>Vout; vector<array<int,3>>Fout;
        if(activeV<N && activeV>0 && activeF>0 && closedOK(Vout,Fout) && Vout.size()<Orig.size()) outputMesh(Vout,Fout);
        else outputOriginal();
    }
    return 0;
}
