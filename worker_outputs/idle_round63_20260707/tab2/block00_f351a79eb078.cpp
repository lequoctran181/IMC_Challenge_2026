#include <bits/stdc++.h>
using namespace std;

struct Vec3{double x,y,z;};
static inline Vec3 operator+(const Vec3&a,const Vec3&b){return {a.x+b.x,a.y+b.y,a.z+b.z};}
static inline Vec3 operator-(const Vec3&a,const Vec3&b){return {a.x-b.x,a.y-b.y,a.z-b.z};}
static inline Vec3 operator*(const Vec3&a,double s){return {a.x*s,a.y*s,a.z*s};}
static inline Vec3 operator/(const Vec3&a,double s){return {a.x/s,a.y/s,a.z/s};}
static inline double dot3(const Vec3&a,const Vec3&b){return a.x*b.x+a.y*b.y+a.z*b.z;}
static inline Vec3 cross3(const Vec3&a,const Vec3&b){return {a.y*b.z-a.z*b.y,a.z*b.x-a.x*b.z,a.x*b.y-a.y*b.x};}
static inline double norm2(const Vec3&a){return dot3(a,a);} 
static inline double norm3(const Vec3&a){return sqrt(max(0.0,norm2(a)));}
struct Face{int a,b,c;};
struct Quadric{double q[10]; Quadric(){memset(q,0,sizeof(q));}};
static inline void qadd(Quadric&Q,const Quadric&R){for(int i=0;i<10;i++)Q.q[i]+=R.q[i];}
static inline Quadric planeQ(const Vec3&p,const Vec3&n){double a=n.x,b=n.y,c=n.z,d=-dot3(n,p); Quadric Q; double v[4]={a,b,c,d}; int k=0; for(int i=0;i<4;i++)for(int j=i;j<4;j++)Q.q[k++]=v[i]*v[j]; return Q;}
static inline double qeval(const Quadric&Q,const Vec3&p){double x=p.x,y=p.y,z=p.z;return Q.q[0]*x*x+2*Q.q[1]*x*y+2*Q.q[2]*x*z+2*Q.q[3]*x+Q.q[4]*y*y+2*Q.q[5]*y*z+2*Q.q[6]*y+Q.q[7]*z*z+2*Q.q[8]*z+Q.q[9];}

int N,M; vector<Vec3>P,P0; vector<Face>F,F0; vector<unsigned char>aliveV,aliveF,guardV; vector<double>rad; double diagL=1,epsH=1;
static inline unsigned long long ekey(int a,int b){if(a>b)swap(a,b);return (unsigned long long)(unsigned int)a<<32 | (unsigned int)b;}
static inline bool hasv(const Face&f,int v){return f.a==v||f.b==v||f.c==v;}
static inline int other(const Face&f,int u,int v){if(f.a!=u&&f.a!=v)return f.a;if(f.b!=u&&f.b!=v)return f.b;return f.c;}
static inline Vec3 fnormRaw(const Face&f,const vector<Vec3>&V){return cross3(V[f.b]-V[f.a],V[f.c]-V[f.a]);}
static inline double farea2(const Face&f,const vector<Vec3>&V){return norm3(fnormRaw(f,V));}

struct Build{
    vector<vector<int>> vf;
    vector<pair<int,int>> edges;
    unordered_map<unsigned long long,array<int,2>> ef;
    bool ok=true;
};

static Build buildTopo(){
    Build B; B.vf.assign(N,{}); B.ef.reserve((size_t)M*3+32);
    for(int i=0;i<M;i++) if(aliveF[i]){
        Face f=F[i];
        if(f.a<0||f.b<0||f.c<0||!aliveV[f.a]||!aliveV[f.b]||!aliveV[f.c]||f.a==f.b||f.a==f.c||f.b==f.c||farea2(f,P)<=diagL*diagL*1e-20){aliveF[i]=0;continue;}
        B.vf[f.a].push_back(i); B.vf[f.b].push_back(i); B.vf[f.c].push_back(i);
        int v[3]={f.a,f.b,f.c};
        for(int e=0;e<3;e++){
            unsigned long long k=ekey(v[e],v[(e+1)%3]);
            auto it=B.ef.find(k);
            if(it==B.ef.end())B.ef[k]={i,-1};
            else if(it->second[1]<0)it->second[1]=i;
            else B.ok=false;
        }
    }
    for(auto &kv:B.ef){ if(kv.second[0]<0||kv.second[1]<0)B.ok=false; int a=(int)(kv.first>>32),b=(int)(kv.first&0xffffffffu); B.edges.push_back({a,b}); }
    return B;
}

static vector<Quadric> computeQuadrics(){
    vector<Quadric> Q(N);
    for(int i=0;i<M;i++) if(aliveF[i]){
        Face f=F[i]; Vec3 cr=fnormRaw(f,P); double l=norm3(cr); if(l<=0)continue; Vec3 n=cr/l;
        Quadric q=planeQ(P[f.a],n); qadd(Q[f.a],q); qadd(Q[f.b],q); qadd(Q[f.c],q);
    }
    return Q;
}

static bool solve3(double A[3][3],double b[3],Vec3&x){
    double m[3][4]; for(int i=0;i<3;i++){for(int j=0;j<3;j++)m[i][j]=A[i][j];m[i][3]=b[i];}
    for(int c=0;c<3;c++){
        int p=c; for(int r=c+1;r<3;r++)if(fabs(m[r][c])>fabs(m[p][c]))p=r;
        if(fabs(m[p][c])<1e-12)return false; if(p!=c)for(int j=c;j<4;j++)swap(m[p][j],m[c][j]);
        double div=m[c][c]; for(int j=c;j<4;j++)m[c][j]/=div;
        for(int r=0;r<3;r++)if(r!=c){double t=m[r][c]; for(int j=c;j<4;j++)m[r][j]-=t*m[c][j];}
    }
    x={m[0][3],m[1][3],m[2][3]}; return isfinite(x.x)&&isfinite(x.y)&&isfinite(x.z);
}
static Vec3 qemPoint(const Quadric&A,const Vec3&u,const Vec3&v,bool&ok){
    double M3[3][3]={{A.q[0],A.q[1],A.q[2]},{A.q[1],A.q[4],A.q[5]},{A.q[2],A.q[5],A.q[7]}};
    double b[3]={-A.q[3],-A.q[6],-A.q[8]}; Vec3 x; ok=solve3(M3,b,x); if(!ok)return (u+v)*0.5; return x;
}

static void markGuards(){
    if(N<=8){fill(guardV.begin(),guardV.end(),1);return;}
    Vec3 mn{1e100,1e100,1e100},mx{-1e100,-1e100,-1e100};
    for(auto&p:P){mn.x=min(mn.x,p.x);mn.y=min(mn.y,p.y);mn.z=min(mn.z,p.z);mx.x=max(mx.x,p.x);mx.y=max(mx.y,p.y);mx.z=max(mx.z,p.z);}    
    vector<Vec3> vn(N,{0,0,0});
    for(auto &f:F){Vec3 cr=fnormRaw(f,P);double l=norm3(cr);if(l>0)cr=cr/l;vn[f.a]=vn[f.a]+cr;vn[f.b]=vn[f.b]+cr;vn[f.c]=vn[f.c]+cr;}
    int S=N>60000?88:(N>25000?72:(N>9000?56:42));
    struct C{double front=-1e300,back=1e300,sil=-1;int fi=-1,bi=-1,si=-1;};
    auto one=[&](int ax,int sgn){
        vector<C> grid(S*S); double rx=max(1e-30,mx.x-mn.x),ry=max(1e-30,mx.y-mn.y),rz=max(1e-30,mx.z-mn.z);
        for(int i=0;i<N;i++){
            double u,v,d; if(ax==0){u=(P[i].y-mn.y)/ry;v=(P[i].z-mn.z)/rz;d=sgn*P[i].x;}else if(ax==1){u=(P[i].x-mn.x)/rx;v=(P[i].z-mn.z)/rz;d=sgn*P[i].y;}else{u=(P[i].x-mn.x)/rx;v=(P[i].y-mn.y)/ry;d=sgn*P[i].z;}
            int iu=min(S-1,max(0,(int)(u*S))),iv=min(S-1,max(0,(int)(v*S))); C&c=grid[iu*S+iv];
            if(d>c.front){c.front=d;c.fi=i;} if(d<c.back){c.back=d;c.bi=i;}
            Vec3 n=vn[i]; double nl=norm3(n); double comp=0; if(nl>0){n=n/nl; comp=(ax==0?fabs(n.x):(ax==1?fabs(n.y):fabs(n.z)));}
            double sil=1-comp; if(sil>c.sil){c.sil=sil;c.si=i;}
        }
        for(auto&c:grid){if(c.fi>=0)guardV[c.fi]=1;if(c.bi>=0)guardV[c.bi]=1;if(c.si>=0&&c.sil>0.33)guardV[c.si]=1;}
    };
    for(int ax=0;ax<3;ax++){one(ax,1);one(ax,-1);}    
    int K=max(8,(int)(sqrt((double)N)*0.65));
    for(int ax=0;ax<3;ax++)for(int sg=-1;sg<=1;sg+=2){
        vector<pair<double,int>> a; a.reserve(N); for(int i=0;i<N;i++)a.push_back({(ax==0?P[i].x:(ax==1?P[i].y:P[i].z))*sg,i});
        nth_element(a.begin(),a.begin()+min(K,N-1),a.end()); for(int i=0;i<min(K,N);i++)guardV[a[i].second]=1;
    }
}

static bool linkOK(int u,int v,const Build&B,int &f1,int &f2,int&o1,int&o2){
    auto it=B.ef.find(ekey(u,v)); if(it==B.ef.end())return false; f1=it->second[0]; f2=it->second[1]; if(f1<0||f2<0||f1==f2)return false;
    o1=other(F[f1],u,v); o2=other(F[f2],u,v); if(o1<0||o2<0||o1==o2)return false;
    unordered_set<int> nu; nu.reserve(B.vf[u].size()*2+8);
    for(int fid:B.vf[u])if(aliveF[fid]){int a[3]={F[fid].a,F[fid].b,F[fid].c};for(int x:a)if(x!=u&&x!=v)nu.insert(x);}    
    int cnt=0; bool h1=false,h2=false; unordered_set<int> seen;
    for(int fid:B.vf[v])if(aliveF[fid]){int a[3]={F[fid].a,F[fid].b,F[fid].c};for(int x:a)if(x!=u&&x!=v&&nu.count(x)&&!seen.count(x)){seen.insert(x);cnt++;if(x==o1)h1=true;if(x==o2)h2=true;}}
    return cnt==2&&h1&&h2;
}

static bool duplicateAfter(int keep,int rem,int f1,int f2,const Build&B){
    vector<int> ids=B.vf[keep]; ids.insert(ids.end(),B.vf[rem].begin(),B.vf[rem].end()); sort(ids.begin(),ids.end()); ids.erase(unique(ids.begin(),ids.end()),ids.end());
    unordered_set<unsigned long long> st; st.reserve(ids.size()*2+16);
    for(int fid:ids) if(aliveF[fid]&&fid!=f1&&fid!=f2){
        int a=F[fid].a,b=F[fid].b,c=F[fid].c; if(a==rem)a=keep;if(b==rem)b=keep;if(c==rem)c=keep; if(a==b||a==c||b==c)continue;
        int x[3]={a,b,c}; sort(x,x+3); unsigned long long h=((unsigned long long)(unsigned int)x[0]<<42)^((unsigned long long)(unsigned int)x[1]<<21)^(unsigned int)x[2];
        if(st.count(h))return true; st.insert(h);
    }
    return false;
}

struct Cand{bool ok=false; double cost=1e100; int keep=-1,rem=-1,f1=-1,f2=-1; Vec3 pos;};
static Cand evalDir(int keep,int rem,int f1,int f2,const Build&B,const vector<Quadric>&Q,const vector<Vec3>&opts,double nmin,double pmax){
    Cand best; best.keep=keep; best.rem=rem; best.f1=f1; best.f2=f2; if(guardV[rem]&&!guardV[keep])return best; if(duplicateAfter(keep,rem,f1,f2,B))return best;
    vector<int> ids=B.vf[keep]; ids.insert(ids.end(),B.vf[rem].begin(),B.vf[rem].end()); sort(ids.begin(),ids.end()); ids.erase(unique(ids.begin(),ids.end()),ids.end());
    Quadric QQ=Q[keep]; qadd(QQ,Q[rem]);
    for(Vec3 pos:opts){
        double du=norm3(pos-P[keep]),dv=norm3(pos-P[rem]); double nr=max(rad[keep]+du,rad[rem]+dv); if(nr>epsH*0.985)continue;
        if(min(norm3(pos-P0[keep]),norm3(pos-P0[rem]))>epsH*0.985)continue;
        if(guardV[rem]&&du+dv>epsH*0.55)continue;
        bool ok=true; double worstN=0,worstP=0,areaLoss=0;
        for(int fid:ids) if(aliveF[fid]&&fid!=f1&&fid!=f2){
            Face old=F[fid], neu=old; if(neu.a==rem)neu.a=keep;if(neu.b==rem)neu.b=keep;if(neu.c==rem)neu.c=keep; if(neu.a==neu.b||neu.a==neu.c||neu.b==neu.c){ok=false;break;}
            Vec3 poa=P[old.a],pob=P[old.b],poc=P[old.c]; Vec3 no=cross3(pob-poa,poc-poa); double alo=norm3(no); if(!(alo>diagL*diagL*1e-20)){ok=false;break;}
            Vec3 pa=(neu.a==keep?pos:P[neu.a]),pb=(neu.b==keep?pos:P[neu.b]),pc=(neu.c==keep?pos:P[neu.c]); Vec3 nn=cross3(pb-pa,pc-pa); double al=norm3(nn); if(!(al>diagL*diagL*1e-20)){ok=false;break;}
            double nd=dot3(no,nn)/(alo*al); if(nd<nmin){ok=false;break;} worstN=max(worstN,1-nd);
            Vec3 un=no/alo; double pd=fabs(dot3(un,pos-P[rem])); if(pd>pmax){ok=false;break;} worstP=max(worstP,pd/pmax); areaLoss=max(areaLoss,max(0.0,1.0-al/alo));
            if(al<alo*1e-6){ok=false;break;}
        }
        if(!ok)continue;
        double cost=qeval(QQ,pos)+0.025*nr/epsH+20*worstN+0.4*worstP+0.08*areaLoss+(guardV[keep]?0.08:0)+(guardV[rem]?0.8:0);
        if(cost<best.cost){best.ok=true;best.cost=cost;best.pos=pos;}
    }
    return best;
}

static Cand evalEdge(int u,int v,const Build&B,const vector<Quadric>&Q,double nmin,double pmax){
    Cand none; if(!aliveV[u]||!aliveV[v])return none; int f1,f2,o1,o2; if(!linkOK(u,v,B,f1,f2,o1,o2))return none;
    Quadric QQ=Q[u]; qadd(QQ,Q[v]); bool sok=false; Vec3 qp=qemPoint(QQ,P[u],P[v],sok); Vec3 mid=(P[u]+P[v])*0.5;
    vector<Vec3> opts={P[u],P[v],mid}; if(sok && norm3(qp-mid)<diagL*0.2)opts.push_back(qp);
    Cand a=evalDir(u,v,f1,f2,B,Q,opts,nmin,pmax), b=evalDir(v,u,f1,f2,B,Q,opts,nmin,pmax);
    if(a.ok&&(!b.ok||a.cost<=b.cost))return a; return b;
}

static void applyCollapse(const Cand&c){
    P[c.keep]=c.pos; rad[c.keep]=max(rad[c.keep]+norm3(c.pos-P[c.keep]), rad[c.rem]+norm3(c.pos-P[c.rem]));
    // previous line used updated P; correct conservative recompute by old values is handled below via max with edge length budget surrogate
    rad[c.keep]=min(epsH*0.999, max(rad[c.keep], rad[c.rem]+norm3(c.pos-P[c.rem])));
    guardV[c.keep]=guardV[c.keep]||guardV[c.rem];
    aliveF[c.f1]=0; aliveF[c.f2]=0;
    for(int fid=0;fid<M;fid++) if(aliveF[fid]){
        if(F[fid].a==c.rem)F[fid].a=c.keep; if(F[fid].b==c.rem)F[fid].b=c.keep; if(F[fid].c==c.rem)F[fid].c=c.keep;
        if(F[fid].a==F[fid].b||F[fid].a==F[fid].c||F[fid].b==F[fid].c)aliveF[fid]=0;
    }
    aliveV[c.rem]=0;
}

static bool validateClosed(){
    unordered_map<unsigned long long,int> cnt; cnt.reserve((size_t)M*3+16); int fn=0;
    for(int i=0;i<M;i++) if(aliveF[i]){
        Face f=F[i]; if(!aliveV[f.a]||!aliveV[f.b]||!aliveV[f.c]||f.a==f.b||f.a==f.c||f.b==f.c||farea2(f,P)<=diagL*diagL*1e-20)return false; fn++;
        int v[3]={f.a,f.b,f.c}; for(int e=0;e<3;e++)cnt[ekey(v[e],v[(e+1)%3])]++;
    }
    if(fn==0)return false; for(auto &kv:cnt)if(kv.second!=2)return false; return true;
}

static void simplify(){
    if(N<=8)return; vector<unsigned char> origAV=aliveV,origAF=aliveF; vector<Vec3> origP=P; vector<Face> origF=F;
    markGuards();
    auto t0=chrono::steady_clock::now(); int rounds=0;
    double nmins[]={0.9997,0.9985,0.996,0.992,0.985,0.975};
    double lens[]={0.018,0.026,0.035,0.045,0.055,0.065};
    while(rounds<6){
        Build B=buildTopo(); if(!B.ok)break; vector<Quadric> Q=computeQuadrics();
        vector<pair<double,pair<int,int>>> order; order.reserve(B.edges.size());
        for(auto&e:B.edges){double l=norm3(P[e.first]-P[e.second]); if(l<=diagL*lens[rounds]) order.push_back({l,e});}
        sort(order.begin(),order.end()); int changed=0, checked=0;
        for(auto &it:order){
            if(++checked>280000)break; int u=it.second.first,v=it.second.second; if(!aliveV[u]||!aliveV[v])continue;
            Cand c=evalEdge(u,v,B,Q,nmins[rounds],diagL*(0.004+0.0025*rounds));
            if(!c.ok)continue; Vec3 oldKeep=P[c.keep], oldRem=P[c.rem]; double oldrk=rad[c.keep], oldrr=rad[c.rem];
            applyCollapse(c); rad[c.keep]=max(oldrk+norm3(c.pos-oldKeep),oldrr+norm3(c.pos-oldRem)); P[c.keep]=c.pos;
            changed++;
            if((changed&127)==0){ B=buildTopo(); if(!B.ok)goto done; Q=computeQuadrics(); }
            if(chrono::duration<double>(chrono::steady_clock::now()-t0).count()>1.9)goto done;
        }
        if(changed<16)rounds++; else if(rounds<3)rounds++; else rounds++;
        if(chrono::duration<double>(chrono::steady_clock::now()-t0).count()>1.9)break;
    }
    done:
    if(!validateClosed()){aliveV=origAV;aliveF=origAF;P=origP;F=origF;}
}

static void output(){
    vector<int> id(N,-1); int vc=0,fc=0; for(int i=0;i<N;i++)if(aliveV[i])id[i]=++vc;
    for(int i=0;i<M;i++)if(aliveF[i]&&aliveV[F[i].a]&&aliveV[F[i].b]&&aliveV[F[i].c]&&F[i].a!=F[i].b&&F[i].a!=F[i].c&&F[i].b!=F[i].c&&farea2(F[i],P)>diagL*diagL*1e-20)fc++;
    cout<<vc<<' '<<fc<<'\n'; cout.setf(ios::fixed); cout<<setprecision(10);
    for(int i=0;i<N;i++)if(aliveV[i])cout<<"v "<<P[i].x<<' '<<P[i].y<<' '<<P[i].z<<'\n';
    for(int i=0;i<M;i++)if(aliveF[i]&&aliveV[F[i].a]&&aliveV[F[i].b]&&aliveV[F[i].c]&&F[i].a!=F[i].b&&F[i].a!=F[i].c&&F[i].b!=F[i].c&&farea2(F[i],P)>diagL*diagL*1e-20)cout<<"f "<<id[F[i].a]<<' '<<id[F[i].b]<<' '<<id[F[i].c]<<'\n';
}

int main(){
    ios::sync_with_stdio(false); cin.tie(nullptr);
    if(!(cin>>N>>M))return 0; P.resize(N);P0.resize(N);F.resize(M);F0.resize(M); Vec3 mn{1e100,1e100,1e100},mx{-1e100,-1e100,-1e100}; char c;
    for(int i=0;i<N;i++){cin>>c>>P[i].x>>P[i].y>>P[i].z;P0[i]=P[i];mn.x=min(mn.x,P[i].x);mn.y=min(mn.y,P[i].y);mn.z=min(mn.z,P[i].z);mx.x=max(mx.x,P[i].x);mx.y=max(mx.y,P[i].y);mx.z=max(mx.z,P[i].z);}    
    for(int i=0;i<M;i++){cin>>c>>F[i].a>>F[i].b>>F[i].c;--F[i].a;--F[i].b;--F[i].c;F0[i]=F[i];}
    diagL=norm3(mx-mn); if(!(diagL>0))diagL=1; epsH=diagL*0.0492; aliveV.assign(N,1);aliveF.assign(M,1);guardV.assign(N,0);rad.assign(N,0.0);
    simplify(); output(); return 0;
}