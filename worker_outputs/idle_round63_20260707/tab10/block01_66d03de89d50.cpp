#include <bits/stdc++.h>
using namespace std;

struct Vec{double x,y,z;};
static inline Vec operator+(Vec a,Vec b){return {a.x+b.x,a.y+b.y,a.z+b.z};}
static inline Vec operator-(Vec a,Vec b){return {a.x-b.x,a.y-b.y,a.z-b.z};}
static inline Vec operator*(Vec a,double s){return {a.x*s,a.y*s,a.z*s};}
static inline Vec operator/(Vec a,double s){return {a.x/s,a.y/s,a.z/s};}
static inline double dotp(Vec a,Vec b){return a.x*b.x+a.y*b.y+a.z*b.z;}
static inline Vec crossp(Vec a,Vec b){return {a.y*b.z-a.z*b.y,a.z*b.x-a.x*b.z,a.x*b.y-a.y*b.x};}
static inline double norm2(Vec a){return dotp(a,a);} 
static inline double norm(Vec a){return sqrt(norm2(a));}
static inline Vec normalized(Vec a){double n=norm(a); return n>0?a/n:Vec{0,0,0};}
struct Face{int a,b,c;};

struct DSU{
    vector<int> p;
    DSU(int n=0){init(n);} void init(int n){p.resize(n); iota(p.begin(),p.end(),0);} 
    int find(int x){while(p[x]!=x){p[x]=p[p[x]]; x=p[x];} return x;}
};

int N,M;
vector<Vec> P0,P;
vector<Face> F;
vector<unsigned char> aliveV, aliveF, protectV;
vector<vector<int>> inc;
vector<double> rad;
double DIAG=1, EPSD=1;
Vec BBmn, BBmx;

static inline array<int,3> sface(int a,int b,int c){array<int,3> r={a,b,c}; sort(r.begin(),r.end()); return r;}
static inline bool hasv(const Face&f,int v){return f.a==v||f.b==v||f.c==v;}
static inline int thirdv(const Face&f,int u,int v){ if(f.a!=u&&f.a!=v) return f.a; if(f.b!=u&&f.b!=v) return f.b; return f.c; }
static inline Vec faceNormal(const Face&f,const vector<Vec>&Q){return crossp(Q[f.b]-Q[f.a],Q[f.c]-Q[f.a]);}

void rebuildInc(){
    inc.assign(P.size(),{});
    for(int i=0;i<(int)F.size();++i) if(aliveF[i]){
        inc[F[i].a].push_back(i); inc[F[i].b].push_back(i); inc[F[i].c].push_back(i);
    }
}

uint64_t edgeKey(int a,int b){ if(a>b) swap(a,b); return (uint64_t)(uint32_t)a<<32 | (uint32_t)b; }

vector<pair<int,int>> collectEdges(){
    vector<pair<int,int>> e; e.reserve(F.size()*3);
    unordered_set<uint64_t> seen; seen.reserve(F.size()*4+7);
    for(int i=0;i<(int)F.size();++i) if(aliveF[i]){
        int a[3]={F[i].a,F[i].b,F[i].c};
        for(int k=0;k<3;k++){
            int u=a[k],v=a[(k+1)%3]; if(!aliveV[u]||!aliveV[v]) continue;
            uint64_t K=edgeKey(u,v); if(seen.insert(K).second) e.push_back({u,v});
        }
    }
    return e;
}

bool edgeTwoFaces(int u,int v,int by[2],int op[2]){
    int c=0;
    for(int fid: inc[u]) if(aliveF[fid] && hasv(F[fid],v)){
        if(c>=2) return false;
        by[c]=fid; op[c]=thirdv(F[fid],u,v); c++;
    }
    return c==2 && op[0]>=0 && op[1]>=0 && op[0]!=op[1];
}

bool linkCondition(int u,int v,int op0,int op1){
    static vector<int> mark; static int stamp=1; static vector<int> seen;
    if((int)mark.size()<N) mark.assign(N,0);
    if(++stamp==INT_MAX){fill(mark.begin(),mark.end(),0); stamp=1;}
    for(int fid:inc[u]) if(aliveF[fid]){
        const Face&f=F[fid]; int a[3]={f.a,f.b,f.c};
        for(int x:a) if(x!=u && x!=v) mark[x]=stamp;
    }
    int cnt=0; bool h0=false,h1=false;
    for(int fid:inc[v]) if(aliveF[fid]){
        const Face&f=F[fid]; int a[3]={f.a,f.b,f.c};
        for(int x:a) if(x!=u && x!=v && mark[x]==stamp){
            if(x!=op0 && x!=op1) return false;
            if(x==op0) h0=true; if(x==op1) h1=true;
        }
    }
    return h0&&h1;
}

struct PlaneQ{double q[10];};
PlaneQ zeroQ(){PlaneQ Q; memset(Q.q,0,sizeof(Q.q)); return Q;}
void addPlane(PlaneQ &Q, Vec n, double d, double w){
    double a=n.x,b=n.y,c=n.z; double v[4]={a,b,c,d}; int idx=0;
    for(int i=0;i<4;i++) for(int j=i;j<4;j++) Q.q[idx++]+=w*v[i]*v[j];
}
double evalQ(const PlaneQ&Q, Vec p){
    double v[4]={p.x,p.y,p.z,1.0}; double s=0; int idx=0;
    for(int i=0;i<4;i++) for(int j=i;j<4;j++){ double m=Q.q[idx++]; s += (i==j?1.0:2.0)*m*v[i]*v[j]; }
    return s;
}
PlaneQ addQ(PlaneQ A,const PlaneQ&B){for(int i=0;i<10;i++) A.q[i]+=B.q[i]; return A;}

vector<PlaneQ> Qv;
void buildQuadrics(){
    Qv.assign(N,zeroQ());
    for(int i=0;i<M;i++) if(aliveF[i]){
        Face f=F[i]; Vec cr=faceNormal(f,P); double ar=norm(cr); if(ar<=1e-30) continue;
        Vec n=cr/ar; double d=-dotp(n,P[f.a]); double w=sqrt(ar)+1e-12;
        addPlane(Qv[f.a],n,d,w); addPlane(Qv[f.b],n,d,w); addPlane(Qv[f.c],n,d,w);
    }
}

void markProtected(){
    protectV.assign(N,0);
    if(N<=1000){return;}
    int G = 42;
    if(N>20000) G=58; if(N>80000) G=76; if(N>180000) G=92;
    Vec ext=BBmx-BBmn; double ex[3]={max(ext.x,1e-12),max(ext.y,1e-12),max(ext.z,1e-12)};
    auto coord=[&](const Vec&p,int ax){ return ax==0?p.x:(ax==1?p.y:p.z); };
    for(int depth=0; depth<3; ++depth){
        int uax=(depth+1)%3, vax=(depth+2)%3;
        vector<int> mnv(G*G,-1), mxv(G*G,-1); vector<double> mnz(G*G,1e300), mxz(G*G,-1e300);
        for(int i=0;i<N;i++) if(aliveV[i]){
            int iu=(int)((coord(P[i],uax)-coord(BBmn,uax))/ex[uax]*(G-1));
            int iv=(int)((coord(P[i],vax)-coord(BBmn,vax))/ex[vax]*(G-1));
            if(iu<0)iu=0; if(iu>=G)iu=G-1; if(iv<0)iv=0; if(iv>=G)iv=G-1;
            int id=iu*G+iv; double z=coord(P[i],depth);
            if(z<mnz[id]){mnz[id]=z; mnv[id]=i;} if(z>mxz[id]){mxz[id]=z; mxv[id]=i;}
        }
        for(int id=0;id<G*G;id++){ if(mnv[id]>=0) protectV[mnv[id]]=1; if(mxv[id]>=0) protectV[mxv[id]]=1; }
    }
    unordered_map<uint64_t,pair<int,int>> mp; mp.reserve(M*3+7);
    for(int i=0;i<M;i++) if(aliveF[i]){
        int a[3]={F[i].a,F[i].b,F[i].c};
        for(int k=0;k<3;k++){ uint64_t K=edgeKey(a[k],a[(k+1)%3]); auto &p=mp[K]; if(!p.first) p.first=i+1; else p.second=i+1; }
    }
    vector<Vec> fn(M); vector<double> fa(M);
    for(int i=0;i<M;i++) if(aliveF[i]){Vec cr=faceNormal(F[i],P); fa[i]=norm(cr); fn[i]=fa[i]>0?cr/fa[i]:Vec{0,0,0};}
    for(auto &kv:mp){int f1=kv.second.first-1,f2=kv.second.second-1; if(f1<0||f2<0) continue; double c=dotp(fn[f1],fn[f2]); if(c<0.55){int a=kv.first>>32,b=(int)(kv.first&0xffffffffu); protectV[a]=protectV[b]=1;}}
}

bool duplicateFaceAfter(int keep,int rem,int skip0,int skip1,const Face&nf){
    auto key=sface(nf.a,nf.b,nf.c);
    int small=keep;
    if(inc[nf.a].size()<inc[small].size()) small=nf.a;
    if(inc[nf.b].size()<inc[small].size()) small=nf.b;
    if(inc[nf.c].size()<inc[small].size()) small=nf.c;
    for(int fid:inc[small]) if(aliveF[fid] && fid!=skip0 && fid!=skip1){
        if(hasv(F[fid],rem)) continue;
        if(sface(F[fid].a,F[fid].b,F[fid].c)==key) return true;
    }
    return false;
}

struct Cand{bool ok=false; double cost=1e300; int keep=-1,rem=-1; Vec pos; double nr=0;};

Cand testCollapseDir(int keep,int rem,const int by[2],double eps, double nlim, bool allowProtected){
    Cand C; C.keep=keep; C.rem=rem;
    if(!aliveV[keep]||!aliveV[rem]) return C;
    if(protectV[rem] && !allowProtected) return C;
    PlaneQ Q=addQ(Qv[keep],Qv[rem]);
    vector<Vec> tries; tries.reserve(5);
    tries.push_back(P[keep]);
    tries.push_back((P[keep]+P[rem])*0.5);
    if(!protectV[keep] && !protectV[rem]) tries.push_back(P[rem]);
    Vec d=P[rem]-P[keep]; double bestt=0.0,bestq=1e300;
    for(int s=1;s<4;s++){double t=s/4.0; Vec x=P[keep]+d*t; double q=evalQ(Q,x); if(q<bestq){bestq=q; bestt=t;}}
    tries.push_back(P[keep]+d*bestt);
    double best=1e300; Vec bestp=P[keep]; double bestr=0;
    vector<int> affected; affected.reserve(inc[keep].size()+inc[rem].size());
    for(int fid:inc[keep]) if(aliveF[fid] && fid!=by[0] && fid!=by[1]) affected.push_back(fid);
    for(int fid:inc[rem]) if(aliveF[fid] && fid!=by[0] && fid!=by[1]) affected.push_back(fid);
    sort(affected.begin(),affected.end()); affected.erase(unique(affected.begin(),affected.end()),affected.end());
    for(Vec pos:tries){
        double nr=max(rad[keep]+norm(pos-P[keep]), rad[rem]+norm(pos-P[rem]));
        if(nr>eps) continue;
        bool bad=false; double fold=0, planeMove=0, areaLoss=0;
        for(int fid:affected){
            Face old=F[fid], nf=old;
            if(nf.a==rem) nf.a=keep; if(nf.b==rem) nf.b=keep; if(nf.c==rem) nf.c=keep;
            if(nf.a==nf.b||nf.a==nf.c||nf.b==nf.c){bad=true;break;}
            Vec oldcr=faceNormal(old,P); double olda=norm(oldcr); if(olda<=1e-28){bad=true;break;}
            Vec pa=(nf.a==keep?pos:P[nf.a]), pb=(nf.b==keep?pos:P[nf.b]), pc=(nf.c==keep?pos:P[nf.c]);
            Vec newcr=crossp(pb-pa,pc-pa); double newa=norm(newcr); if(newa<=max(1e-30,olda*1e-9)){bad=true;break;}
            double co=dotp(oldcr,newcr)/(olda*newa); if(co<-1)co=-1;if(co>1)co=1;
            if(co<nlim){bad=true;break;}
            fold=max(fold,1.0-co);
            Vec on=oldcr/olda; planeMove=max(planeMove,fabs(dotp(on,pos-P[old.a])));
            if(newa<olda) areaLoss=max(areaLoss,1.0-newa/olda);
            if(duplicateFaceAfter(keep,rem,by[0],by[1],nf)){bad=true;break;}
        }
        if(bad) continue;
        double q=evalQ(Q,pos);
        double len=norm(P[keep]-P[rem]);
        double prot = protectV[keep]? -0.15:0.0;
        double cost = q/(DIAG*DIAG+1e-30) + 18.0*fold + 4.0*planeMove/(eps+1e-30) + 0.02*areaLoss + 0.000001*inc[rem].size() + prot + len/(eps+1e-30)*0.002;
        if(cost<best){best=cost; bestp=pos; bestr=nr;}
    }
    if(best<1e299){C.ok=true; C.cost=best; C.pos=bestp; C.nr=bestr;}
    return C;
}

Cand bestCollapse(int u,int v,double eps,double nlim,bool allowProtected){
    int by[2]={-1,-1},op[2]={-1,-1}; Cand none;
    if(!edgeTwoFaces(u,v,by,op)) return none;
    if(!linkCondition(u,v,op[0],op[1])) return none;
    Cand a=testCollapseDir(u,v,by,eps,nlim,allowProtected);
    Cand b=testCollapseDir(v,u,by,eps,nlim,allowProtected);
    if(a.ok && (!b.ok || a.cost<=b.cost)) return a;
    return b;
}

void applyCollapse(const Cand&C,int by0,int by1){
    int keep=C.keep, rem=C.rem;
    aliveF[by0]=0; aliveF[by1]=0;
    for(int fid:inc[rem]) if(aliveF[fid]){
        if(F[fid].a==rem) F[fid].a=keep;
        if(F[fid].b==rem) F[fid].b=keep;
        if(F[fid].c==rem) F[fid].c=keep;
    }
    aliveV[rem]=0; P[keep]=C.pos; rad[keep]=C.nr; protectV[keep]=protectV[keep]||protectV[rem];
}

int aliveCount(){int c=0; for(auto x:aliveV)c+=x?1:0; return c;}

void compactPrint(){
    vector<int> id(N,-1); vector<Vec> outv; outv.reserve(N);
    for(int i=0;i<N;i++) if(aliveV[i]){id[i]=(int)outv.size()+1; outv.push_back(P[i]);}
    vector<Face> outf; outf.reserve(F.size()); set<array<int,3>> seen;
    for(int i=0;i<(int)F.size();i++) if(aliveF[i]){
        int a=id[F[i].a],b=id[F[i].b],c=id[F[i].c];
        if(a<=0||b<=0||c<=0||a==b||a==c||b==c) continue;
        array<int,3> k=sface(a,b,c); if(seen.insert(k).second) outf.push_back({a,b,c});
    }
    cout.setf(ios::fixed); cout<<setprecision(10);
    cout<<outv.size()<<" "<<outf.size()<<"\n";
    for(auto &p:outv) cout<<"v "<<p.x<<" "<<p.y<<" "<<p.z<<"\n";
    for(auto &f:outf) cout<<"f "<<f.a<<" "<<f.b<<" "<<f.c<<"\n";
}

void simplify(){
    aliveV.assign(N,1); aliveF.assign(M,1); rad.assign(N,0.0);
    rebuildInc(); buildQuadrics(); markProtected();
    int target = max(8, (int)(N*0.34));
    if(N>50000) target=max(8,(int)(N*0.28));
    if(N>150000) target=max(8,(int)(N*0.23));
    struct Pass{double epsMul,nlim; bool allowProt; int maxScan;};
    vector<Pass> passes={
        {0.48,0.985,false,2},{0.62,0.970,false,2},{0.78,0.945,false,3},{0.92,0.910,false,3},
        {0.98,0.875,true,3}
    };
    auto start=chrono::steady_clock::now();
    for(auto ps:passes){
        for(int round=0; round<ps.maxScan; ++round){
            rebuildInc(); buildQuadrics();
            auto edges=collectEdges();
            vector<pair<double,pair<int,int>>> ord; ord.reserve(edges.size());
            for(auto &e:edges){
                double l2=norm2(P[e.first]-P[e.second]);
                double bias=(protectV[e.first]+protectV[e.second])*DIAG*DIAG*0.02;
                ord.push_back({l2+bias,e});
            }
            sort(ord.begin(),ord.end(),[](auto&a,auto&b){return a.first<b.first;});
            int changed=0, checked=0; double eps=EPSD*ps.epsMul;
            for(auto &it:ord){
                if(aliveCount()<=target) break;
                if((++checked&2047)==0){ double t=chrono::duration<double>(chrono::steady_clock::now()-start).count(); if(t>1.82) return; }
                int u=it.second.first,v=it.second.second; if(!aliveV[u]||!aliveV[v]) continue;
                int by[2]={-1,-1},op[2]={-1,-1}; if(!edgeTwoFaces(u,v,by,op)) continue; if(!linkCondition(u,v,op[0],op[1])) continue;
                Cand c=bestCollapse(u,v,eps,ps.nlim,ps.allowProt); if(!c.ok) continue;
                applyCollapse(c,by[0],by[1]); changed++;
            }
            if(!changed) break;
        }
    }
}

int main(){
    ios::sync_with_stdio(false); cin.tie(nullptr);
    if(!(cin>>N>>M)) return 0;
    P0.resize(N); P.resize(N); F.resize(M);
    BBmn={1e100,1e100,1e100}; BBmx={-1e100,-1e100,-1e100};
    for(int i=0;i<N;i++){
        char ch; cin>>ch>>P[i].x>>P[i].y>>P[i].z; P0[i]=P[i];
        BBmn.x=min(BBmn.x,P[i].x); BBmn.y=min(BBmn.y,P[i].y); BBmn.z=min(BBmn.z,P[i].z);
        BBmx.x=max(BBmx.x,P[i].x); BBmx.y=max(BBmx.y,P[i].y); BBmx.z=max(BBmx.z,P[i].z);
    }
    for(int i=0;i<M;i++){char ch; cin>>ch>>F[i].a>>F[i].b>>F[i].c; --F[i].a;--F[i].b;--F[i].c;}
    DIAG=norm(BBmx-BBmn); if(!(DIAG>0)) DIAG=1.0; EPSD=0.05*DIAG;
    if(N<=1000){
        cout.setf(ios::fixed); cout<<setprecision(10);
        cout<<N<<" "<<M<<"\n";
        for(auto&p:P) cout<<"v "<<p.x<<" "<<p.y<<" "<<p.z<<"\n";
        for(auto&f:F) cout<<"f "<<f.a+1<<" "<<f.b+1<<" "<<f.c+1<<"\n";
        return 0;
    }
    simplify();
    compactPrint();
    return 0;
}