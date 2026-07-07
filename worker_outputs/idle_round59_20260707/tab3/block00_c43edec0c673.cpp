#include <bits/stdc++.h>
using namespace std;

struct Vec3{double x,y,z;};
static inline Vec3 operator+(Vec3 a,Vec3 b){return {a.x+b.x,a.y+b.y,a.z+b.z};}
static inline Vec3 operator-(Vec3 a,Vec3 b){return {a.x-b.x,a.y-b.y,a.z-b.z};}
static inline Vec3 operator*(Vec3 a,double s){return {a.x*s,a.y*s,a.z*s};}
static inline Vec3 operator/(Vec3 a,double s){return {a.x/s,a.y/s,a.z/s};}
static inline double dotv(Vec3 a,Vec3 b){return a.x*b.x+a.y*b.y+a.z*b.z;}
static inline Vec3 crossv(Vec3 a,Vec3 b){return {a.y*b.z-a.z*b.y,a.z*b.x-a.x*b.z,a.x*b.y-a.y*b.x};}
static inline double n2(Vec3 a){return dotv(a,a);} static inline double norm(Vec3 a){return sqrt(max(0.0,n2(a)));}
static inline Vec3 unit(Vec3 a){double l=norm(a); return l>1e-300?a/l:Vec3{0,0,0};}
struct Face{int a,b,c; bool on=true;};
struct Q{double q[10]; Q(){memset(q,0,sizeof(q));}
 void addPlane(double a,double b,double c,double d,double w=1){q[0]+=w*a*a;q[1]+=w*a*b;q[2]+=w*a*c;q[3]+=w*a*d;q[4]+=w*b*b;q[5]+=w*b*c;q[6]+=w*b*d;q[7]+=w*c*c;q[8]+=w*c*d;q[9]+=w*d*d;}
 void add(const Q&o){for(int i=0;i<10;i++)q[i]+=o.q[i];}
 double eval(Vec3 p)const{double x=p.x,y=p.y,z=p.z;return q[0]*x*x+2*q[1]*x*y+2*q[2]*x*z+2*q[3]*x+q[4]*y*y+2*q[5]*y*z+2*q[6]*y+q[7]*z*z+2*q[8]*z+q[9];}
};
static inline uint64_t keye(int a,int b){if(a>b)swap(a,b);return (uint64_t)(uint32_t)a<<32 | (uint32_t)b;}
static inline array<int,3> keyt(int a,int b,int c){array<int,3> r{a,b,c}; sort(r.begin(),r.end()); return r;}
struct Cand{double score; int u,v,vu,vv; bool operator<(Cand const&o)const{return score>o.score;}};

int N,M; vector<Vec3>P0,P; vector<Face>F,F0; vector<vector<int>> inc; vector<char> alive; vector<int> ver, markv; int stamp=1; vector<Q> quad; vector<array<double,3>> bminv,bmaxv; vector<int> csz; double diagLen=1, tau=1, tau2=1; chrono::steady_clock::time_point T0;

static inline bool time_ok(){return chrono::duration<double>(chrono::steady_clock::now()-T0).count()<18.8;}
static inline bool hasv(Face const&f,int v){return f.a==v||f.b==v||f.c==v;}
static inline Vec3 fnorm(Face const&f){return unit(crossv(P[f.b]-P[f.a],P[f.c]-P[f.a]));}
static inline Vec3 fnorm0(Face const&f){return unit(crossv(P0[f.b]-P0[f.a],P0[f.c]-P0[f.a]));}
static bool solveQ(const Q&q,Vec3&out){
 double a00=q.q[0],a01=q.q[1],a02=q.q[2],a11=q.q[4],a12=q.q[5],a22=q.q[7];
 double b0=-q.q[3],b1=-q.q[6],b2=-q.q[8];
 double det=a00*(a11*a22-a12*a12)-a01*(a01*a22-a12*a02)+a02*(a01*a12-a11*a02);
 if(fabs(det)<1e-15)return false;
 double dx=b0*(a11*a22-a12*a12)-a01*(b1*a22-a12*b2)+a02*(b1*a12-a11*b2);
 double dy=a00*(b1*a22-a12*b2)-b0*(a01*a22-a12*a02)+a02*(a01*b2-b1*a02);
 double dz=a00*(a11*b2-b1*a12)-a01*(a01*b2-b1*a02)+b0*(a01*a12-a11*a02);
 out={dx/det,dy/det,dz/det}; return isfinite(out.x)&&isfinite(out.y)&&isfinite(out.z);
}
static bool coverBox(Vec3 p,array<double,3> mn,array<double,3> mx,double factor=1.0){
 double lim=tau2*factor*factor;
 for(int m=0;m<8;m++){double x=(m&1)?mx[0]:mn[0],y=(m&2)?mx[1]:mn[1],z=(m&4)?mx[2]:mn[2]; if((p.x-x)*(p.x-x)+(p.y-y)*(p.y-y)+(p.z-z)*(p.z-z)>lim)return false;}
 return true;
}
static void project(Vec3 p,int view,double&u,double&v,double&z){
 const double D=2.5, foc=800.0; double sx,sy; int ax=view/2, s=(view&1)?-1:1;
 if(ax==0){sx=p.y; sy=p.z; z=D-s*p.x; if(s<0)sx=-sx;}
 else if(ax==1){sx=p.x; sy=p.z; z=D-s*p.y; if(s>0)sx=-sx;}
 else {sx=p.x; sy=p.y; z=D-s*p.z; if(s<0)sx=-sx;}
 u=foc*sx/z+512.0; v=foc*sy/z+512.0;
}
static bool viewBox(Vec3 p,array<double,3> mn,array<double,3> mx,double pix){
 for(int view=0;view<6;view++){
   double pu,pv,pz; project(p,view,pu,pv,pz); if(pz<=0)return false;
   double far2=0;
   for(int m=0;m<8;m++){Vec3 c{(m&1)?mx[0]:mn[0],(m&2)?mx[1]:mn[1],(m&4)?mx[2]:mn[2]}; double u,v,z; project(c,view,u,v,z); if(z<=0)return false; far2=max(far2,(u-pu)*(u-pu)+(v-pv)*(v-pv));}
   if(far2>pix*pix)return false;
 }
 return true;
}
static void cleanInc(int v){ if(!alive[v])return; vector<int> r; r.reserve(inc[v].size()); for(int id:inc[v]) if(id>=0&&id<(int)F.size()&&F[id].on&&hasv(F[id],v)) r.push_back(id); inc[v].swap(r); }
static bool edgeFaces(int a,int b,int ids[2],int opp[2]){
 int c=0; if(inc[a].size()>inc[b].size())swap(a,b);
 for(int id:inc[a]) if(F[id].on && hasv(F[id],a)&&hasv(F[id],b)){
   if(c>=2)return false; ids[c]=id; Face &f=F[id]; int o=f.a^f.b^f.c^a^b; opp[c]=o; c++;
 }
 return c==2 && opp[0]!=opp[1];
}
static bool linkOK(int a,int b){
 int ids[2],opp[2]; if(!edgeFaces(a,b,ids,opp))return false;
 if(++stamp>1000000000){fill(markv.begin(),markv.end(),0); stamp=1;}
 int s1=stamp++, s2=stamp++;
 for(int id:inc[a]) if(F[id].on&&hasv(F[id],a)){Face &f=F[id]; int vs[3]={f.a,f.b,f.c}; for(int x:vs) if(x!=a&&x!=b) markv[x]=s1;}
 int cnt=0, g0=0,g1=0;
 for(int id:inc[b]) if(F[id].on&&hasv(F[id],b)){Face &f=F[id]; int vs[3]={f.a,f.b,f.c}; for(int x:vs) if(x!=a&&x!=b&&markv[x]==s1){markv[x]=s2; cnt++; if(x==opp[0])g0=1; if(x==opp[1])g1=1; if(cnt>2)return false;}}
 return cnt==2&&g0&&g1;
}
static bool duplicateAfter(int keep,int rem,int idskip0,int idskip1){
 set<array<int,3>> seen;
 for(int i=0;i<(int)F.size();i++) if(F[i].on&&i!=idskip0&&i!=idskip1){
   int a=F[i].a,b=F[i].b,c=F[i].c; if(a==rem)a=keep;if(b==rem)b=keep;if(c==rem)c=keep; if(a==b||b==c||a==c)return true; auto k=keyt(a,b,c); if(!seen.insert(k).second)return true;
 }
 return false;
}
static bool normalOK(int keep,int rem,Vec3 np,double coslim){
 static vector<int> ids; ids.clear(); ids.insert(ids.end(),inc[keep].begin(),inc[keep].end()); ids.insert(ids.end(),inc[rem].begin(),inc[rem].end()); sort(ids.begin(),ids.end()); ids.erase(unique(ids.begin(),ids.end()),ids.end());
 for(int id:ids) if(F[id].on){Face f=F[id]; bool hk=hasv(f,keep), hr=hasv(f,rem); if(hk&&hr)continue; Vec3 old=fnorm(f); int vs[3]={f.a,f.b,f.c}; for(int&i:vs)if(i==rem)i=keep; Vec3 a=(vs[0]==keep)?np:P[vs[0]],b=(vs[1]==keep)?np:P[vs[1]],c=(vs[2]==keep)?np:P[vs[2]]; Vec3 cr=crossv(b-a,c-a); if(n2(cr)<1e-28*diagLen*diagLen*diagLen*diagLen)return false; Vec3 nw=unit(cr); if(dotv(old,nw)<coslim)return false; }
 return true;
}
static bool candPos(int u,int v,Vec3 &best,double &score){
 array<double,3> mn,mx; for(int k=0;k<3;k++){mn[k]=min(bminv[u][k],bminv[v][k]); mx[k]=max(bmaxv[u][k],bmaxv[v][k]);}
 Q q=quad[u]; q.add(quad[v]); Vec3 opt; vector<Vec3> cs; if(solveQ(q,opt))cs.push_back(opt); cs.push_back(P[u]); cs.push_back(P[v]); cs.push_back((P[u]+P[v])*0.5); cs.push_back({(mn[0]+mx[0])*0.5,(mn[1]+mx[1])*0.5,(mn[2]+mx[2])*0.5}); cs.push_back(P[u]*0.75+P[v]*0.25); cs.push_back(P[u]*0.25+P[v]*0.75);
 double pix = (N<3000?36:(N<30000?28:(N<120000?22:17)));
 bool ok=false; score=1e300;
 for(Vec3 p:cs){ if(!coverBox(p,mn,mx,0.998))continue; if(!viewBox(p,mn,mx,pix))continue; double s=q.eval(p)+1e-9*n2(P[u]-P[v])+1e-7*(csz[u]+csz[v]); if(s<score){score=s;best=p;ok=true;}}
 return ok;
}
static void pushEdge(priority_queue<Cand>&pq,int a,int b){ if(a==b||!alive[a]||!alive[b])return; Vec3 p; double s; if(candPos(a,b,p,s)) pq.push({s,a,b,ver[a],ver[b]}); }
static bool collapse(priority_queue<Cand>&pq,int a,int b){
 if(a==b||!alive[a]||!alive[b])return false; int ids[2],opp[2]; if(!edgeFaces(a,b,ids,opp))return false; if(!linkOK(a,b))return false; Vec3 np; double sc; if(!candPos(a,b,np,sc))return false;
 int keep=a, rem=b; if(inc[b].size()+csz[b]*2 > inc[a].size()+csz[a]*2){keep=b;rem=a;}
 if(duplicateAfter(keep,rem,ids[0],ids[1]))return false; if(!normalOK(keep,rem,np,cos((N<8000?42:32)*acos(-1.0)/180.0)))return false;
 for(int id:ids) if(F[id].on){F[id].on=false;}
 vector<int> rlist=inc[rem]; P[keep]=np; quad[keep].add(quad[rem]); csz[keep]+=csz[rem];
 for(int k=0;k<3;k++){bminv[keep][k]=min(bminv[keep][k],bminv[rem][k]); bmaxv[keep][k]=max(bmaxv[keep][k],bmaxv[rem][k]);}
 for(int id:rlist) if(F[id].on&&hasv(F[id],rem)){ if(F[id].a==rem)F[id].a=keep; if(F[id].b==rem)F[id].b=keep; if(F[id].c==rem)F[id].c=keep; inc[keep].push_back(id); }
 alive[rem]=0; ver[keep]++; ver[rem]++; cleanInc(keep); inc[rem].clear();
 static vector<int> nb; nb.clear(); if(++stamp>1000000000){fill(markv.begin(),markv.end(),0);stamp=1;} int st=stamp++;
 for(int id:inc[keep]) if(F[id].on){int vs[3]={F[id].a,F[id].b,F[id].c}; for(int x:vs)if(x!=keep&&alive[x]&&markv[x]!=st){markv[x]=st; nb.push_back(x);}}
 for(int x:nb)pushEdge(pq,keep,x);
 return true;
}
static int aliveCount(){int c=0; for(char a:alive)c+=a; return c;}
static bool validate(vector<int>&id, vector<Face>&of){
 id.assign(N,-1); int n=0; for(int i=0;i<N;i++)if(alive[i])id[i]=n++;
 of.clear(); set<array<int,3>> ts; unordered_map<uint64_t,int> ec; ec.reserve(F.size()*3+9);
 for(auto &f:F) if(f.on){ if(!alive[f.a]||!alive[f.b]||!alive[f.c])return false; int a=id[f.a],b=id[f.b],c=id[f.c]; if(a<0||b<0||c<0||a==b||b==c||a==c)return false; Vec3 cr=crossv(P[f.b]-P[f.a],P[f.c]-P[f.a]); if(n2(cr)<1e-30*diagLen*diagLen*diagLen*diagLen)return false; auto kt=keyt(a,b,c); if(!ts.insert(kt).second)return false; of.push_back({a,b,c,true}); ec[keye(a,b)]++;ec[keye(b,c)]++;ec[keye(c,a)]++; }
 if(n<4||of.size()<4)return false; for(auto &kv:ec)if(kv.second!=2)return false; return true;
}
static bool tryPrimitiveSphere(){
 if(N<900)return false; Vec3 mn=P0[0],mx=P0[0],ctr{0,0,0}; for(auto p:P0){mn.x=min(mn.x,p.x);mn.y=min(mn.y,p.y);mn.z=min(mn.z,p.z);mx.x=max(mx.x,p.x);mx.y=max(mx.y,p.y);mx.z=max(mx.z,p.z);ctr=ctr+p;} ctr=ctr/(double)N; double ex=mx.x-mn.x,ey=mx.y-mn.y,ez=mx.z-mn.z; if(min({ex,ey,ez})<max({ex,ey,ez})*0.965)return false; double r=0; for(auto p:P0)r+=norm(p-ctr); r/=N; double mxerr=0,rm=0; for(auto p:P0){double e=fabs(norm(p-ctr)-r); mxerr=max(mxerr,e); rm+=e*e;} rm=sqrt(rm/N); if(mxerr>tau*0.35||rm>tau*0.08)return false;
 int lat=(N<3000?16:(N<20000?20:24)), lon=lat*2; vector<Vec3> nv; vector<Face> nf; nv.push_back({ctr.x,ctr.y,ctr.z+r}); for(int i=1;i<lat;i++){double th=acos(-1.0)*i/lat; for(int j=0;j<lon;j++){double ph=2*acos(-1.0)*j/lon; nv.push_back({ctr.x+r*sin(th)*cos(ph),ctr.y+r*sin(th)*sin(ph),ctr.z+r*cos(th)});}} int bot=nv.size(); nv.push_back({ctr.x,ctr.y,ctr.z-r}); auto vid=[&](int ring,int j){return 1+(ring-1)*lon+(j%lon+lon)%lon;}; auto add=[&](int a,int b,int c){nf.push_back({a,b,c,true});}; for(int j=0;j<lon;j++)add(0,vid(1,j),vid(1,j+1)); for(int i=1;i<lat-1;i++)for(int j=0;j<lon;j++){add(vid(i,j),vid(i+1,j),vid(i,j+1));add(vid(i,j+1),vid(i+1,j),vid(i+1,j+1));} for(int j=0;j<lon;j++)add(bot,vid(lat-1,j+1),vid(lat-1,j)); if((int)nv.size()>=N)return false; P=nv; F=nf; N=P.size(); alive.assign(N,1); return true;
}
int main(){
 ios::sync_with_stdio(false); cin.tie(nullptr); T0=chrono::steady_clock::now(); if(!(cin>>N>>M))return 0; P.resize(N);P0.resize(N); string s; for(int i=0;i<N;i++){cin>>s>>P[i].x>>P[i].y>>P[i].z; P0[i]=P[i];} F.resize(M);F0.resize(M); for(int i=0;i<M;i++){cin>>s>>F[i].a>>F[i].b>>F[i].c; --F[i].a;--F[i].b;--F[i].c; F0[i]=F[i];}
 if(N==8&&M==12){cout<<"8 12\n"; cout.setf(ios::fixed); cout<<setprecision(9); for(auto p:P0)cout<<"v "<<p.x<<" "<<p.y<<" "<<p.z<<"\n"; for(auto f:F0)cout<<"f "<<f.a+1<<" "<<f.b+1<<" "<<f.c+1<<"\n"; return 0;}
 Vec3 mn=P0[0],mx=P0[0]; for(auto p:P0){mn.x=min(mn.x,p.x);mn.y=min(mn.y,p.y);mn.z=min(mn.z,p.z);mx.x=max(mx.x,p.x);mx.y=max(mx.y,p.y);mx.z=max(mx.z,p.z);} diagLen=norm(mx-mn); if(diagLen<=0)diagLen=1; tau=.05*diagLen; tau2=tau*tau;
 vector<Vec3> saveP=P; vector<Face> saveF=F; int saveN=N;
 if(tryPrimitiveSphere()){
   cout<<N<<" "<<F.size()<<"\n"; cout.setf(ios::fixed); cout<<setprecision(10); for(auto p:P)cout<<"v "<<p.x<<" "<<p.y<<" "<<p.z<<"\n"; for(auto f:F)cout<<"f "<<f.a+1<<" "<<f.b+1<<" "<<f.c+1<<"\n"; return 0;
 }
 N=saveN; P=saveP; F=saveF;
 alive.assign(N,1); ver.assign(N,0); markv.assign(N,0); csz.assign(N,1); bminv.resize(N); bmaxv.resize(N); for(int i=0;i<N;i++){bminv[i]={P[i].x,P[i].y,P[i].z}; bmaxv[i]=bminv[i];}
 inc.assign(N,{}); vector<int> deg(N); for(auto &f:F){deg[f.a]++;deg[f.b]++;deg[f.c]++;} for(int i=0;i<N;i++)inc[i].reserve(deg[i]+8); for(int i=0;i<M;i++){inc[F[i].a].push_back(i);inc[F[i].b].push_back(i);inc[F[i].c].push_back(i);} quad.assign(N,Q()); for(auto &f:F){Vec3 cr=crossv(P[f.b]-P[f.a],P[f.c]-P[f.a]); double ar=norm(cr); if(ar<=1e-300)continue; Vec3 n=cr/ar; double d=-dotv(n,P[f.a]); double w=max(1e-12,ar); quad[f.a].addPlane(n.x,n.y,n.z,d,w);quad[f.b].addPlane(n.x,n.y,n.z,d,w);quad[f.c].addPlane(n.x,n.y,n.z,d,w);}
 priority_queue<Cand> pq; unordered_set<uint64_t> es; es.reserve(M*3+9); for(auto &f:F){es.insert(keye(f.a,f.b));es.insert(keye(f.b,f.c));es.insert(keye(f.c,f.a));} for(auto k:es)pushEdge(pq,(int)(k>>32),(int)(uint32_t)k);
 int target=max(4,(int)ceil(N*(N<1000?0.22:(N<5000?0.16:(N<30000?0.115:0.085))))); long long pops=0;
 while(!pq.empty()&&aliveCount()>target&&time_ok()){Cand c=pq.top();pq.pop(); if(!alive[c.u]||!alive[c.v]||ver[c.u]!=c.vu||ver[c.v]!=c.vv)continue; collapse(pq,c.u,c.v); if((++pops&4095)==0){for(int i=0;i<N;i++)if(alive[i]&&inc[i].size()>256)cleanInc(i);} }
 vector<int> id; vector<Face> of; if(!validate(id,of)){cout<<saveN<<" "<<M<<"\n"; cout.setf(ios::fixed); cout<<setprecision(10); for(auto p:P0)cout<<"v "<<p.x<<" "<<p.y<<" "<<p.z<<"\n"; for(auto f:F0)cout<<"f "<<f.a+1<<" "<<f.b+1<<" "<<f.c+1<<"\n"; return 0;}
 cout<<aliveCount()<<" "<<of.size()<<"\n"; cout.setf(ios::fixed); cout<<setprecision(10); for(int i=0;i<N;i++)if(alive[i])cout<<"v "<<P[i].x<<" "<<P[i].y<<" "<<P[i].z<<"\n"; for(auto f:of)cout<<"f "<<f.a+1<<" "<<f.b+1<<" "<<f.c+1<<"\n";
 return 0;
}