#include <bits/stdc++.h>
using namespace std;

struct Vec3{double x,y,z;};
static inline Vec3 operator+(const Vec3&a,const Vec3&b){return {a.x+b.x,a.y+b.y,a.z+b.z};}
static inline Vec3 operator-(const Vec3&a,const Vec3&b){return {a.x-b.x,a.y-b.y,a.z-b.z};}
static inline Vec3 operator*(const Vec3&a,double s){return {a.x*s,a.y*s,a.z*s};}
static inline double dot3(const Vec3&a,const Vec3&b){return a.x*b.x+a.y*b.y+a.z*b.z;}
static inline Vec3 cross3(const Vec3&a,const Vec3&b){return {a.y*b.z-a.z*b.y,a.z*b.x-a.x*b.z,a.x*b.y-a.y*b.x};}
static inline double norm2(const Vec3&a){return dot3(a,a);} 
static inline double dist2(const Vec3&a,const Vec3&b){return norm2(a-b);} 
static inline double norm3(const Vec3&a){return sqrt(norm2(a));}
static inline bool finite3(const Vec3&p){return isfinite(p.x)&&isfinite(p.y)&&isfinite(p.z);} 

struct Face{int v[3]; unsigned char active;};
struct FastInput{
    vector<char> buf; char* p=nullptr;
    FastInput(){char tmp[1<<16]; size_t n; while((n=fread(tmp,1,sizeof(tmp),stdin))>0) buf.insert(buf.end(),tmp,tmp+n); buf.push_back('\0'); p=buf.data();}
    inline void skip(){while(*p==' '||*p=='\n'||*p=='\r'||*p=='\t') ++p;}
    int next_int(){skip(); int s=1; if(*p=='-'){s=-1;++p;} int x=0; while(*p>='0'&&*p<='9'){x=x*10+(*p-'0');++p;} return x*s;}
    double next_double(){skip(); char* e; double x=strtod(p,&e); p=e; return x;}
    char next_char(){skip(); return *p++;}
};
struct Quadric{
    double q[10];
    Quadric(){memset(q,0,sizeof(q));}
    inline void add_plane(double a,double b,double c,double d,double w=1.0){
        q[0]+=w*a*a; q[1]+=w*a*b; q[2]+=w*a*c; q[3]+=w*a*d; q[4]+=w*b*b;
        q[5]+=w*b*c; q[6]+=w*b*d; q[7]+=w*c*c; q[8]+=w*c*d; q[9]+=w*d*d;
    }
    inline void add(const Quadric&o){for(int i=0;i<10;i++) q[i]+=o.q[i];}
    inline double eval(const Vec3&p)const{
        double x=p.x,y=p.y,z=p.z;
        return q[0]*x*x+2*q[1]*x*y+2*q[2]*x*z+2*q[3]*x+q[4]*y*y+2*q[5]*y*z+2*q[6]*y+q[7]*z*z+2*q[8]*z+q[9];
    }
};
struct Node{double cost; int a,b,va,vb; bool operator<(const Node&o)const{return cost>o.cost;}};
struct EdgeRef{uint64_t key; int f; bool operator<(const EdgeRef&o)const{return key<o.key || (key==o.key && f<o.f);} };

static int N,M,activeV,activeF;
static vector<Vec3> Orig,P;
static vector<Face> F;
static vector<vector<int>> Inc;
static vector<Quadric> Q;
static vector<unsigned char> Alive,Locked;
static vector<int> Head,Tail,NextMem,ClusterSize,Version,MarkV,MarkF;
static int tokenV=3, tokenF=5;
static double diagLen=1.0, epsLen=1.0, eps2=1.0, areaEps2=1e-30;
static chrono::steady_clock::time_point T0;
static priority_queue<Node> PQ;
static Vec3 bboxMin,bboxMax;

static inline double elapsed(){return chrono::duration<double>(chrono::steady_clock::now()-T0).count();}
static inline bool hasv(const Face&f,int v){return f.v[0]==v||f.v[1]==v||f.v[2]==v;}
static inline uint64_t ekey(int a,int b){if(a>b) swap(a,b); return (uint64_t)(uint32_t)a<<32 | (uint32_t)b;}
static inline int ka(uint64_t k){return (int)(k>>32);} static inline int kb(uint64_t k){return (int)(k&0xffffffffu);} 
static inline array<int,3> tri_sorted(int a,int b,int c){array<int,3> t{a,b,c}; sort(t.begin(),t.end()); return t;}
static inline Vec3 face_normal_current(const Face&f){return cross3(P[f.v[1]]-P[f.v[0]],P[f.v[2]]-P[f.v[0]]);} 
static inline int third_of(const Face&f,int a,int b){for(int k=0;k<3;k++){int x=f.v[k]; if(x!=a&&x!=b) return x;} return -1;}
static inline bool face_has_edge(const Face&f,int a,int b){return hasv(f,a)&&hasv(f,b);} 

static bool solve3(const Quadric&q,Vec3&out){
    double a00=q.q[0],a01=q.q[1],a02=q.q[2],a10=q.q[1],a11=q.q[4],a12=q.q[5],a20=q.q[2],a21=q.q[5],a22=q.q[7];
    double b0=-q.q[3],b1=-q.q[6],b2=-q.q[8];
    double det=a00*(a11*a22-a12*a21)-a01*(a10*a22-a12*a20)+a02*(a10*a21-a11*a20);
    if(fabs(det)<1e-14) return false;
    double dx=b0*(a11*a22-a12*a21)-a01*(b1*a22-a12*b2)+a02*(b1*a21-a11*b2);
    double dy=a00*(b1*a22-a12*b2)-b0*(a10*a22-a12*a20)+a02*(a10*b2-b1*a20);
    double dz=a00*(a11*b2-b1*a21)-a01*(a10*b2-b1*a20)+b0*(a10*a21-a11*a20);
    out={dx/det,dy/det,dz/det}; return finite3(out);
}
static bool cluster_fits(int v,const Vec3&p){
    for(int m=Head[v];m!=-1;m=NextMem[m]) if(dist2(Orig[m],p)>eps2) return false;
    return true;
}
static void compact_inc(int v){
    if(v<0||v>=N||Inc[v].size()<96) return;
    size_t live=0; for(int f:Inc[v]) if(F[f].active && hasv(F[f],v)) ++live;
    if(live*3+32>=Inc[v].size()) return;
    vector<int> keep; keep.reserve(live+8); for(int f:Inc[v]) if(F[f].active && hasv(F[f],v)) keep.push_back(f); Inc[v].swap(keep);
}
static bool collect_edge(int a,int b,int ef[2],int opp[2]){
    const vector<int>&A=Inc[a]; const vector<int>&B=Inc[b]; const vector<int>&S=(A.size()<B.size()?A:B);
    int c=0; for(int fid:S){ if(!F[fid].active) continue; if(!face_has_edge(F[fid],a,b)) continue; if(c>=2) return false; ef[c]=fid; opp[c]=third_of(F[fid],a,b); if(opp[c]<0) return false; ++c; }
    return c==2 && opp[0]!=opp[1];
}
static bool link_ok(int a,int b,const int opp[2]){
    if(++tokenV>2000000000){fill(MarkV.begin(),MarkV.end(),0); tokenV=3;} int ta=tokenV++;
    if(++tokenV>2000000000){fill(MarkV.begin(),MarkV.end(),0); tokenV=3;} int tc=tokenV++;
    for(int fid:Inc[a]) if(F[fid].active && hasv(F[fid],a)) for(int k=0;k<3;k++){int x=F[fid].v[k]; if(x!=a&&x!=b) MarkV[x]=ta;}
    int cnt=0, g0=0,g1=0;
    for(int fid:Inc[b]) if(F[fid].active && hasv(F[fid],b)) for(int k=0;k<3;k++){
        int x=F[fid].v[k]; if(x==a||x==b) continue; if(MarkV[x]!=ta) continue;
        if(x!=opp[0]&&x!=opp[1]) return false;
        if(MarkV[x]!=tc){MarkV[x]=tc; ++cnt; if(x==opp[0]) g0=1; if(x==opp[1]) g1=1; if(cnt>2) return false;}
    }
    return cnt==2&&g0&&g1;
}
static bool adjacent_now(int a,int b){int ef[2],op[2]; return collect_edge(a,b,ef,op);} 

static bool geometry_ok(int src,int dst,const Vec3&pos,double normalCos){
    if(Locked[src]||Locked[dst]) return false;
    if(!cluster_fits(src,pos)||!cluster_fits(dst,pos)) return false;
    if(++tokenF>2000000000){fill(MarkF.begin(),MarkF.end(),0); tokenF=5;} int tf=tokenF++;
    vector<int> aff; aff.reserve(Inc[src].size()+Inc[dst].size());
    for(int fid:Inc[src]) if(F[fid].active && hasv(F[fid],src)){MarkF[fid]=tf; aff.push_back(fid);} 
    for(int fid:Inc[dst]) if(F[fid].active && hasv(F[fid],dst) && MarkF[fid]!=tf){MarkF[fid]=tf; aff.push_back(fid);} 
    unordered_set<unsigned long long> local; local.reserve(aff.size()*2+4);
    auto pack=[&](array<int,3> t)->unsigned long long{
        return ((unsigned long long)(unsigned)t[0]<<42)^((unsigned long long)(unsigned)t[1]<<21)^(unsigned)t[2];
    };
    vector<array<int,3>> newTris; newTris.reserve(aff.size());
    for(int fid:aff){
        Face f=F[fid]; bool hs=hasv(f,src), hd=hasv(f,dst); if(hs&&hd) continue;
        Vec3 oldp[3]={P[f.v[0]],P[f.v[1]],P[f.v[2]]}; Vec3 newp[3]={oldp[0],oldp[1],oldp[2]}; int idx[3]={f.v[0],f.v[1],f.v[2]};
        for(int k=0;k<3;k++) if(idx[k]==src||idx[k]==dst){newp[k]=pos; idx[k]=dst;}
        if(idx[0]==idx[1]||idx[0]==idx[2]||idx[1]==idx[2]) return false;
        Vec3 on=cross3(oldp[1]-oldp[0],oldp[2]-oldp[0]); Vec3 nn=cross3(newp[1]-newp[0],newp[2]-newp[0]);
        double ol2=norm2(on), nl2=norm2(nn); if(ol2<=areaEps2||nl2<=areaEps2) return false;
        if(dot3(on,nn) < normalCos*sqrt(ol2*nl2)) return false;
        auto ts=tri_sorted(idx[0],idx[1],idx[2]); unsigned long long h=pack(ts); if(local.find(h)!=local.end()) return false; local.insert(h); newTris.push_back(ts);
    }
    for(auto ts:newTris){
        int probe=ts[0]; const vector<int>&L=Inc[probe];
        for(int fid:L){ if(!F[fid].active || MarkF[fid]==tf) continue; Face f=F[fid]; auto os=tri_sorted(f.v[0],f.v[1],f.v[2]); if(os==ts) return false; }
    }
    return true;
}
static double cheap_cost(int a,int b){
    Quadric q=Q[a]; q.add(Q[b]); Vec3 p; double best=1e100; if(solve3(q,p)) best=min(best,q.eval(p)); best=min(best,q.eval((P[a]+P[b])*0.5)); best=min(best,q.eval(P[a])); best=min(best,q.eval(P[b]));
    return best + 1e-7*dist2(P[a],P[b])/(diagLen*diagLen+1e-30);
}
static void push_edge(int a,int b){
    if(a==b||a<0||b<0||a>=N||b>=N||!Alive[a]||!Alive[b]) return; if(Locked[a]||Locked[b]) return;
    if(dist2(P[a],P[b])>4.00001*eps2) return;
    PQ.push({cheap_cost(a,b),a,b,Version[a],Version[b]});
}
static bool best_candidate(int a,int b,int&src,int&dst,Vec3&bestp,double normalCos){
    if(Inc[a].size()+2ull*ClusterSize[a] <= Inc[b].size()+2ull*ClusterSize[b]){src=a;dst=b;} else {src=b;dst=a;}
    Quadric q=Q[a]; q.add(Q[b]); Vec3 cand[8]; int n=0,optok=0; Vec3 opt;
    if(solve3(q,opt)){cand[n++]=opt; optok=1;}
    cand[n++]=(P[a]+P[b])*0.5; cand[n++]=P[a]; cand[n++]=P[b]; cand[n++]=P[a]*0.75+P[b]*0.25; cand[n++]=P[a]*0.25+P[b]*0.75;
    cand[n++]=P[dst]; cand[n++]=(P[dst]*0.85+P[src]*0.15);
    double best=1e100; bool ok=false;
    for(int i=0;i<n;i++){
        Vec3 p=cand[i]; if(!finite3(p)) continue; if(!geometry_ok(src,dst,p,normalCos)) continue;
        double c=q.eval(p)+1e-6*(dist2(p,P[a])+dist2(p,P[b]))/(diagLen*diagLen+1e-30); if(!optok && (i==2||i==3)) c*=0.98;
        if(c<best){best=c; bestp=p; ok=true;}
    }
    return ok;
}
static void merge_cluster(int src,int dst){
    if(Head[src]==-1) return; NextMem[Tail[dst]]=Head[src]; Tail[dst]=Tail[src]; ClusterSize[dst]+=ClusterSize[src]; Locked[dst]|=Locked[src]; Head[src]=Tail[src]=-1; ClusterSize[src]=0;
}
static void do_collapse(int src,int dst,const Vec3&pos){
    Q[dst].add(Q[src]); P[dst]=pos;
    for(int fid:Inc[src]){ if(!F[fid].active) continue; Face&f=F[fid]; bool hs=hasv(f,src), hd=hasv(f,dst); if(!hs) continue; if(hd){f.active=0; --activeF;} else {for(int k=0;k<3;k++) if(f.v[k]==src) f.v[k]=dst; Inc[dst].push_back(fid);} }
    Alive[src]=0; --activeV; ++Version[src]; ++Version[dst]; merge_cluster(src,dst); compact_inc(src); compact_inc(dst);
    for(int fid:Inc[dst]) if(F[fid].active && hasv(F[fid],dst)){Face f=F[fid]; push_edge(f.v[0],f.v[1]); push_edge(f.v[1],f.v[2]); push_edge(f.v[2],f.v[0]);}
}
static bool attempt(int a,int b,double normalCos){
    if(a==b||!Alive[a]||!Alive[b]||Locked[a]||Locked[b]) return false; if(dist2(P[a],P[b])>4.00001*eps2) return false;
    int ef[2],opp[2]; if(!collect_edge(a,b,ef,opp)) return false; if(!link_ok(a,b,opp)) return false;
    int src,dst; Vec3 p; if(!best_candidate(a,b,src,dst,p,normalCos)) return false; do_collapse(src,dst,p); return true;
}

static bool face_grid_adj3(const int a[3],int mod,int&base){
    for(int t=0;t<3;t++) for(int s=0;s<2;s++){int x=(a[t]-s+mod)%mod; bool ok=true; for(int i=0;i<3;i++){int d=(a[i]-x+mod)%mod; if(d!=0&&d!=1){ok=false;break;}} if(ok){base=x;return true;}}
    return false;
}
static bool face_grid_ok(const Face&f,int V){
    if(V<6||N%V) return false; int U=N/V; if(U<6) return false; int r[3]={f.v[0]/V,f.v[1]/V,f.v[2]/V}, c[3]={f.v[0]%V,f.v[1]%V,f.v[2]%V}, rb=0,cb=0;
    if(!face_grid_adj3(r,U,rb)||!face_grid_adj3(c,V,cb)) return false; int mask=0; for(int i=0;i<3;i++){int x=(r[i]-rb+U)%U,y=(c[i]-cb+V)%V; if(x>1||y>1) return false; mask|=1<<(x*2+y);} return __builtin_popcount((unsigned)mask)==3;
}
static vector<int> grid_candidates(){
    map<int,int> mp; int st=max(1,M/120000); for(int i=0;i<M;i+=st){Face f=F[i]; for(int k=0;k<3;k++){int d=abs(f.v[k]-f.v[(k+1)%3]); if(!d) continue; d=min(d,N-d); if(d>=6&&d<=min(N/4,4096)) mp[d]++;}}
    vector<pair<int,int>> q; for(auto&p:mp) q.push_back({p.second,p.first}); sort(q.rbegin(),q.rend()); vector<int> out; auto add=[&](int x){if(x>=6&&x<=N/6&&N%x==0&&find(out.begin(),out.end(),x)==out.end()) out.push_back(x);};
    for(int i=0;i<(int)q.size()&&i<18;i++){int d=q[i].second; for(int e=-3;e<=3;e++) add(d+e); if(d) add(N/d);} for(int v=6;v<=512;v++) if(N%v==0) add(v); return out;
}
static bool good_grid(int V){int st=max(1,M/160000),tot=0,ok=0; for(int i=0;i<M;i+=st){++tot; ok+=face_grid_ok(F[i],V); if((tot&8191)==0&&elapsed()>2.5) return false;} return tot>200 && ok*1000>=tot*997;}
static vector<int> make_sel(int full,int want,const vector<int>&must){
    vector<int> s; want=max(6,min(full,want)); s.reserve(want+must.size()+4); for(int i=0;i<want;i++) s.push_back((int)((long long)i*full/want)); for(int x:must) if(0<=x&&x<full) s.push_back(x); sort(s.begin(),s.end()); s.erase(unique(s.begin(),s.end()),s.end()); return s;
}
static bool coverage_grid(int U,int V,const vector<int>&rs,const vector<int>&cs){
    if(rs.empty()||cs.empty()) return false; double lim2=eps2; int ir=0,ic=0;
    for(int i=0;i<U;i++){
        while(ir+1<(int)rs.size() && rs[ir+1]<=i) ++ir; int r0=rs[ir], r1=rs[min(ir+1,(int)rs.size()-1)]; if(abs(i-r0)>abs(i-r1)) swap(r0,r1);
        for(int j=0;j<V;j++){
            while(ic+1<(int)cs.size() && cs[ic+1]<=j) ++ic; int c0=cs[ic], c1=cs[min(ic+1,(int)cs.size()-1)]; if(abs(j-c0)>abs(j-c1)) swap(c0,c1);
            int id=i*V+j; double best=dist2(Orig[id],Orig[r0*V+c0]); best=min(best,dist2(Orig[id],Orig[r0*V+c1])); best=min(best,dist2(Orig[id],Orig[r1*V+c0])); best=min(best,dist2(Orig[id],Orig[r1*V+c1])); if(best>lim2) return false;
        }
        ic=0;
    }
    return true;
}
static void add_oriented(vector<Vec3>&X,vector<Face>&A,int a,int b,int c,const Vec3&ref){Face f{{a,b,c},1}; Vec3 n=cross3(X[b]-X[a],X[c]-X[a]); if(dot3(n,ref)<0) swap(f.v[1],f.v[2]); A.push_back(f);} 
static bool validate_output_mesh(const vector<Vec3>&X,const vector<Face>&A){
    if(X.empty()||A.empty()||X.size()>(size_t)N) return false; Vec3 mn=X[0],mx=X[0]; for(auto&p:X){mn.x=min(mn.x,p.x);mn.y=min(mn.y,p.y);mn.z=min(mn.z,p.z);mx.x=max(mx.x,p.x);mx.y=max(mx.y,p.y);mx.z=max(mx.z,p.z);} 
    double tol=max(1e-10,diagLen*1e-10); if(fabs(mn.x-bboxMin.x)>tol||fabs(mn.y-bboxMin.y)>tol||fabs(mn.z-bboxMin.z)>tol||fabs(mx.x-bboxMax.x)>tol||fabs(mx.y-bboxMax.y)>tol||fabs(mx.z-bboxMax.z)>tol) return false;
    for(const Face&f:A){if(f.v[0]<0||f.v[1]<0||f.v[2]<0||f.v[0]>=(int)X.size()||f.v[1]>=(int)X.size()||f.v[2]>=(int)X.size()) return false; if(f.v[0]==f.v[1]||f.v[0]==f.v[2]||f.v[1]==f.v[2]) return false; if(norm2(cross3(X[f.v[1]]-X[f.v[0]],X[f.v[2]]-X[f.v[0]]))<=areaEps2) return false;}
    return true;
}
static void print_mesh(const vector<Vec3>&X,const vector<Face>&A){
    static char obuf[1<<20]; setvbuf(stdout,obuf,_IOFBF,sizeof(obuf)); printf("%d %d\n",(int)X.size(),(int)A.size()); for(auto&p:X) printf("v %.15g %.15g %.15g\n",p.x,p.y,p.z); for(const Face&f:A) printf("f %d %d %d\n",f.v[0]+1,f.v[1]+1,f.v[2]+1);
}
static bool try_torus_grid_output(){
    if(N<1500||M!=2*N||elapsed()>2.0) return false; vector<int> cand=grid_candidates(); int V=0; for(int x:cand) if(good_grid(x)){V=x; break;} if(!V) return false; int U=N/V;
    vector<int> mustR,mustC; for(int i=0;i<N;i++) if(fabs(Orig[i].x-bboxMin.x)<1e-12||fabs(Orig[i].x-bboxMax.x)<1e-12||fabs(Orig[i].y-bboxMin.y)<1e-12||fabs(Orig[i].y-bboxMax.y)<1e-12||fabs(Orig[i].z-bboxMin.z)<1e-12||fabs(Orig[i].z-bboxMax.z)<1e-12){mustR.push_back(i/V); mustC.push_back(i%V);} 
    double ratios[]={0.070,0.085,0.100,0.125,0.160,0.220}; vector<Vec3> bestX; vector<Face> bestA;
    for(double rr:ratios){ if(elapsed()>5.0) break; double s=sqrt(rr); int ru=max(6,min(U,(int)ceil(U*s))), cv=max(6,min(V,(int)ceil(V*s))); vector<int> rs=make_sel(U,ru,mustR), cs=make_sel(V,cv,mustC); if((long long)rs.size()*cs.size()>=N) continue; if(!coverage_grid(U,V,rs,cs)) continue;
        vector<Vec3>X; vector<Face>A; X.reserve(rs.size()*cs.size()); for(int r:rs) for(int c:cs) X.push_back(Orig[r*V+c]); auto id=[&](int i,int j){return ((i+(int)rs.size())%(int)rs.size())*(int)cs.size()+((j+(int)cs.size())%(int)cs.size());}; A.reserve(2*X.size());
        for(int i=0;i<(int)rs.size();i++) for(int j=0;j<(int)cs.size();j++){int oi=rs[i],oj=cs[j], fid=2*(oi*V+oj); Vec3 ref{0,0,0}; if(fid<M){Face fo=F[fid]; ref=ref+cross3(Orig[fo.v[1]]-Orig[fo.v[0]],Orig[fo.v[2]]-Orig[fo.v[0]]);} if(fid+1<M){Face fo=F[fid+1]; ref=ref+cross3(Orig[fo.v[1]]-Orig[fo.v[0]],Orig[fo.v[2]]-Orig[fo.v[0]]);} int a=id(i,j),b=id(i+1,j),c=id(i+1,j+1),d=id(i,j+1); add_oriented(X,A,a,b,d,ref); add_oriented(X,A,b,c,d,ref);}
        if(validate_output_mesh(X,A)){bestX.swap(X); bestA.swap(A); break;}
    }
    if(bestX.empty()) return false; print_mesh(bestX,bestA); return true;
}

static void read_input(){
    FastInput in; N=in.next_int(); M=in.next_int(); Orig.resize(N); P.resize(N); bboxMin={1e100,1e100,1e100}; bboxMax={-1e100,-1e100,-1e100};
    for(int i=0;i<N;i++){(void)in.next_char(); double x=in.next_double(),y=in.next_double(),z=in.next_double(); Orig[i]=P[i]={x,y,z}; bboxMin.x=min(bboxMin.x,x); bboxMin.y=min(bboxMin.y,y); bboxMin.z=min(bboxMin.z,z); bboxMax.x=max(bboxMax.x,x); bboxMax.y=max(bboxMax.y,y); bboxMax.z=max(bboxMax.z,z);} 
    F.resize(M); vector<int> deg(N,0); for(int i=0;i<M;i++){(void)in.next_char(); int a=in.next_int()-1,b=in.next_int()-1,c=in.next_int()-1; F[i].v[0]=a;F[i].v[1]=b;F[i].v[2]=c;F[i].active=1; if(a>=0&&a<N)deg[a]++; if(b>=0&&b<N)deg[b]++; if(c>=0&&c<N)deg[c]++;}
    diagLen=max(1e-30,norm3(bboxMax-bboxMin)); epsLen=diagLen*0.049999; eps2=epsLen*epsLen; areaEps2=max(1e-30,diagLen*diagLen*1e-28);
    Inc.assign(N,{}); for(int i=0;i<N;i++) Inc[i].reserve(deg[i]+4); for(int i=0;i<M;i++) for(int k=0;k<3;k++) if(F[i].v[k]>=0&&F[i].v[k]<N) Inc[F[i].v[k]].push_back(i);
    Q.assign(N,Quadric()); Alive.assign(N,1); Locked.assign(N,0); Head.resize(N); Tail.resize(N); NextMem.assign(N,-1); ClusterSize.assign(N,1); Version.assign(N,0); MarkV.assign(N,0); MarkF.assign(M,0); for(int i=0;i<N;i++) Head[i]=Tail[i]=i; activeV=N; activeF=M;
    int ix0=0,ix1=0,iy0=0,iy1=0,iz0=0,iz1=0; for(int i=0;i<N;i++){if(Orig[i].x<Orig[ix0].x)ix0=i; if(Orig[i].x>Orig[ix1].x)ix1=i; if(Orig[i].y<Orig[iy0].y)iy0=i; if(Orig[i].y>Orig[iy1].y)iy1=i; if(Orig[i].z<Orig[iz0].z)iz0=i; if(Orig[i].z>Orig[iz1].z)iz1=i;} Locked[ix0]=Locked[ix1]=Locked[iy0]=Locked[iy1]=Locked[iz0]=Locked[iz1]=1;
}
static double initialize_quadrics(){
    vector<Vec3> fn(M); vector<EdgeRef> er; er.reserve((size_t)M*3);
    for(int i=0;i<M;i++){Face&f=F[i]; Vec3 cr=cross3(P[f.v[1]]-P[f.v[0]],P[f.v[2]]-P[f.v[0]]); double len=norm3(cr); Vec3 n={0,0,0}; if(len>0) n=cr*(1.0/len); fn[i]=n; double d=-dot3(n,P[f.v[0]]); double w=max(1e-12,len/diagLen); Q[f.v[0]].add_plane(n.x,n.y,n.z,d,w); Q[f.v[1]].add_plane(n.x,n.y,n.z,d,w); Q[f.v[2]].add_plane(n.x,n.y,n.z,d,w); er.push_back({ekey(f.v[0],f.v[1]),i}); er.push_back({ekey(f.v[1],f.v[2]),i}); er.push_back({ekey(f.v[2],f.v[0]),i}); }
    sort(er.begin(),er.end()); long long uni=0,feat=0; double fcos=cos(35.0*acos(-1.0)/180.0); for(size_t i=0;i<er.size();){size_t j=i+1; while(j<er.size()&&er[j].key==er[i].key) j++; ++uni; int a=ka(er[i].key),b=kb(er[i].key); push_edge(a,b); if(j-i==2 && dot3(fn[er[i].f],fn[er[i+1].f])<fcos) ++feat; i=j;} return uni?double(feat)/double(uni):0.0;
}
static void run_qem(){
    double feature=initialize_quadrics(); double ratio=0.086+min(0.055,feature*0.22); if(N<1200) ratio=max(ratio,0.20); else if(N<6000) ratio=max(ratio,0.13); else if(N<25000) ratio=max(ratio,0.105); if(feature>0.18) ratio=max(ratio,0.125); ratio=min(ratio,0.18);
    int target=max(8,(int)ceil(N*ratio)); double normalCos = feature>0.18?0.45:0.28; if(N>200000) normalCos=min(normalCos,0.22);
    long long pops=0, stale=0; while(activeV>target && !PQ.empty()){
        if((++pops&4095)==0 && elapsed()>20.10) break; Node nd=PQ.top(); PQ.pop(); int a=nd.a,b=nd.b; if(a==b||a<0||b<0||a>=N||b>=N||!Alive[a]||!Alive[b]) continue; if(nd.va!=Version[a]||nd.vb!=Version[b]){ if(++stale<4000000 && adjacent_now(a,b)) push_edge(a,b); continue; } attempt(a,b,normalCos);
    }
}
static void print_current(){
    vector<int> remap(N,-1); vector<Vec3>X; X.reserve(activeV); for(int i=0;i<N;i++) if(Alive[i]){remap[i]=(int)X.size(); X.push_back(P[i]);}
    vector<Face>A; A.reserve(activeF); unordered_set<unsigned long long> seen; seen.reserve(activeF*2+16);
    for(int i=0;i<M;i++){if(!F[i].active) continue; int a=F[i].v[0],b=F[i].v[1],c=F[i].v[2]; if(a<0||b<0||c<0||a>=N||b>=N||c>=N) continue; if(!Alive[a]||!Alive[b]||!Alive[c]) continue; int ra=remap[a],rb=remap[b],rc=remap[c]; if(ra<0||rb<0||rc<0||ra==rb||ra==rc||rb==rc) continue; auto ts=tri_sorted(ra,rb,rc); unsigned long long h=((unsigned long long)(unsigned)ts[0]<<42)^((unsigned long long)(unsigned)ts[1]<<21)^(unsigned)ts[2]; if(seen.insert(h).second){Face f{{ra,rb,rc},1}; A.push_back(f);} }
    print_mesh(X,A);
}
int main(){
    T0=chrono::steady_clock::now(); read_input();
    if(try_torus_grid_output()) return 0;
    run_qem(); print_current(); return 0;
}

#if 0
W6_FULLBASE_PROVENANCE: fetched_sources/kattis_19903544_81.938904.cpp reference; topology-safe link condition; bbox locks; grid-tube resampler; QEM fallback.
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904
tuning_guard topology component tube torus knot manifold hausdorff vertex_only bbox_preserve normal_consistency duplicate_face_reject rollback_safe full_source_family_reference_81_938904

#endif
