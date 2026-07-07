#include <bits/stdc++.h>
using namespace std;

struct Vec3{double x,y,z;};
static inline Vec3 operator+(Vec3 a,Vec3 b){return {a.x+b.x,a.y+b.y,a.z+b.z};}
static inline Vec3 operator-(Vec3 a,Vec3 b){return {a.x-b.x,a.y-b.y,a.z-b.z};}
static inline Vec3 operator*(Vec3 a,double s){return {a.x*s,a.y*s,a.z*s};}
static inline Vec3 operator/(Vec3 a,double s){return {a.x/s,a.y/s,a.z/s};}
static inline double dotv(Vec3 a,Vec3 b){return a.x*b.x+a.y*b.y+a.z*b.z;}
static inline Vec3 crossv(Vec3 a,Vec3 b){return {a.y*b.z-a.z*b.y,a.z*b.x-a.x*b.z,a.x*b.y-a.y*b.x};}
static inline double n2(Vec3 a){return dotv(a,a);} 
static inline double normv(Vec3 a){return sqrt(max(0.0,n2(a)));}
static inline Vec3 unitv(Vec3 a){double l=normv(a); return l>1e-300?a/l:Vec3{0,0,0};}
struct Face{int a,b,c; bool on=true;};
struct Quadric{double q[10]; Quadric(){memset(q,0,sizeof(q));}
  void addPlane(double a,double b,double c,double d,double w=1.0){q[0]+=w*a*a;q[1]+=w*a*b;q[2]+=w*a*c;q[3]+=w*a*d;q[4]+=w*b*b;q[5]+=w*b*c;q[6]+=w*b*d;q[7]+=w*c*c;q[8]+=w*c*d;q[9]+=w*d*d;}
  void add(const Quadric&o){for(int i=0;i<10;i++)q[i]+=o.q[i];}
  double eval(Vec3 p)const{double x=p.x,y=p.y,z=p.z;return q[0]*x*x+2*q[1]*x*y+2*q[2]*x*z+2*q[3]*x+q[4]*y*y+2*q[5]*y*z+2*q[6]*y+q[7]*z*z+2*q[8]*z+q[9];}
};
static inline uint64_t ekey(int a,int b){ if(a>b) swap(a,b); return (uint64_t)(uint32_t)a<<32 | (uint32_t)b; }
static inline array<int,3> tkey(int a,int b,int c){ array<int,3> r{a,b,c}; sort(r.begin(),r.end()); return r; }
struct HeapEdge{double cost; int u,v,vu,vv; bool operator<(HeapEdge const&o)const{return cost>o.cost;}};

int N,M; vector<Vec3> P0,P; vector<Face> F0,F; vector<vector<int>> inc; vector<char> alive; vector<int> ver, markv, clusterSize; int stamp=1; vector<Quadric> Qv; vector<array<double,3>> mnBox,mxBox; double diagLen=1.0,tau=1.0,tau2=1.0; chrono::steady_clock::time_point startTime;
static inline bool timeLeft(){return chrono::duration<double>(chrono::steady_clock::now()-startTime).count()<18.7;}
static inline bool hasv(Face const&f,int v){return f.a==v||f.b==v||f.c==v;}
static inline Vec3 faceNormal(const Face&f,const vector<Vec3>&V){return unitv(crossv(V[f.b]-V[f.a],V[f.c]-V[f.a]));}
static void compactInc(int v){ if(v<0||v>=N||!alive[v])return; vector<int> r; r.reserve(inc[v].size()); for(int id:inc[v]) if(id>=0&&id<(int)F.size()&&F[id].on&&hasv(F[id],v)) r.push_back(id); inc[v].swap(r); }
static bool solveOptimal(const Quadric&q,Vec3&out){
  double a00=q.q[0],a01=q.q[1],a02=q.q[2],a11=q.q[4],a12=q.q[5],a22=q.q[7]; double b0=-q.q[3],b1=-q.q[6],b2=-q.q[8];
  double det=a00*(a11*a22-a12*a12)-a01*(a01*a22-a12*a02)+a02*(a01*a12-a11*a02); if(fabs(det)<1e-15) return false;
  double dx=b0*(a11*a22-a12*a12)-a01*(b1*a22-a12*b2)+a02*(b1*a12-a11*b2);
  double dy=a00*(b1*a22-a12*b2)-b0*(a01*a22-a12*a02)+a02*(a01*b2-b1*a02);
  double dz=a00*(a11*b2-b1*a12)-a01*(a01*b2-b1*a02)+b0*(a01*a12-a11*a02);
  out={dx/det,dy/det,dz/det}; return isfinite(out.x)&&isfinite(out.y)&&isfinite(out.z);
}
static bool coverBox(Vec3 p,const array<double,3>&mn,const array<double,3>&mx,double mul){
  double lim=tau2*mul*mul; for(int m=0;m<8;m++){double x=(m&1)?mx[0]:mn[0],y=(m&2)?mx[1]:mn[1],z=(m&4)?mx[2]:mn[2]; double d=(p.x-x)*(p.x-x)+(p.y-y)*(p.y-y)+(p.z-z)*(p.z-z); if(d>lim) return false;} return true;
}
static void project1024(Vec3 p,int view,double&u,double&v,double&z){
  const double D=2.5,Fo=800.0; int ax=view/2,sg=(view&1)?-1:1; double sx,sy;
  if(ax==0){sx=p.y;sy=p.z;z=D-sg*p.x;if(sg<0)sx=-sx;} else if(ax==1){sx=p.x;sy=p.z;z=D-sg*p.y;if(sg>0)sx=-sx;} else {sx=p.x;sy=p.y;z=D-sg*p.z;if(sg<0)sx=-sx;} u=Fo*sx/z+512.0; v=Fo*sy/z+512.0;
}
static bool viewBoxOK(Vec3 p,const array<double,3>&mn,const array<double,3>&mx,double pix){
  for(int view=0;view<6;view++){double pu,pv,pz; project1024(p,view,pu,pv,pz); if(!(pz>0)) return false; double far2=0; for(int m=0;m<8;m++){Vec3 c{(m&1)?mx[0]:mn[0],(m&2)?mx[1]:mn[1],(m&4)?mx[2]:mn[2]}; double u,v,z; project1024(c,view,u,v,z); if(!(z>0)) return false; far2=max(far2,(u-pu)*(u-pu)+(v-pv)*(v-pv));} if(far2>pix*pix) return false;} return true;
}
static bool edgeFaces(int a,int b,int ids[2],int opp[2]){
  int cnt=0; if(inc[a].size()>inc[b].size()) swap(a,b); for(int id:inc[a]) if(F[id].on&&hasv(F[id],a)&&hasv(F[id],b)){ if(cnt>=2) return false; ids[cnt]=id; Face &f=F[id]; int o=-1; if(f.a!=a&&f.a!=b)o=f.a; if(f.b!=a&&f.b!=b)o=f.b; if(f.c!=a&&f.c!=b)o=f.c; opp[cnt]=o; cnt++; } return cnt==2&&opp[0]>=0&&opp[1]>=0&&opp[0]!=opp[1];
}
static bool linkCondition(int a,int b){
  int ids[2],opp[2]; if(!edgeFaces(a,b,ids,opp)) return false; if(++stamp>1000000000){fill(markv.begin(),markv.end(),0); stamp=1;} int s1=stamp++,s2=stamp++; for(int id:inc[a]) if(F[id].on&&hasv(F[id],a)){int vs[3]={F[id].a,F[id].b,F[id].c}; for(int x:vs) if(x!=a&&x!=b) markv[x]=s1;} int c=0,g0=0,g1=0; for(int id:inc[b]) if(F[id].on&&hasv(F[id],b)){int vs[3]={F[id].a,F[id].b,F[id].c}; for(int x:vs) if(x!=a&&x!=b&&markv[x]==s1){markv[x]=s2; c++; if(x==opp[0])g0=1; if(x==opp[1])g1=1; if(c>2)return false;}} return c==2&&g0&&g1;
}
static bool duplicateOrDegenerateAfter(int keep,int rem,int d0,int d1){
  set<array<int,3>> S; for(int i=0;i<(int)F.size();i++) if(F[i].on&&i!=d0&&i!=d1){int a=F[i].a,b=F[i].b,c=F[i].c; if(a==rem)a=keep;if(b==rem)b=keep;if(c==rem)c=keep; if(a==b||b==c||a==c) return true; auto k=tkey(a,b,c); if(!S.insert(k).second) return true;} return false;
}
static bool normalGuard(int keep,int rem,Vec3 np,double cosLimit){
  static vector<int> ids; ids.clear(); ids.insert(ids.end(),inc[keep].begin(),inc[keep].end()); ids.insert(ids.end(),inc[rem].begin(),inc[rem].end()); sort(ids.begin(),ids.end()); ids.erase(unique(ids.begin(),ids.end()),ids.end()); double minA=max(1e-32,1e-28*diagLen*diagLen*diagLen*diagLen);
  for(int id:ids) if(F[id].on){Face f=F[id]; bool hk=hasv(f,keep),hr=hasv(f,rem); if(hk&&hr) continue; Vec3 oldN=faceNormal(f,P); int vs[3]={f.a,f.b,f.c}; for(int&i:vs) if(i==rem) i=keep; Vec3 a=(vs[0]==keep)?np:P[vs[0]],b=(vs[1]==keep)?np:P[vs[1]],c=(vs[2]==keep)?np:P[vs[2]]; Vec3 cr=crossv(b-a,c-a); if(n2(cr)<=minA) return false; Vec3 nw=unitv(cr); if(dotv(oldN,nw)<cosLimit) return false;} return true;
}
static bool candidatePos(int u,int v,Vec3&best,double&score){
  array<double,3> mn,mx; for(int k=0;k<3;k++){mn[k]=min(mnBox[u][k],mnBox[v][k]); mx[k]=max(mxBox[u][k],mxBox[v][k]);}
  Quadric q=Qv[u]; q.add(Qv[v]); vector<Vec3> cand; Vec3 opt; if(solveOptimal(q,opt)) cand.push_back(opt); cand.push_back(P[u]); cand.push_back(P[v]); cand.push_back((P[u]+P[v])*0.5); cand.push_back({(mn[0]+mx[0])*0.5,(mn[1]+mx[1])*0.5,(mn[2]+mx[2])*0.5}); cand.push_back(P[u]*0.7+P[v]*0.3); cand.push_back(P[u]*0.3+P[v]*0.7);
  double pix=N<3000?38:(N<30000?30:(N<120000?24:18)); bool ok=false; score=1e300; for(Vec3 p:cand){ if(!coverBox(p,mn,mx,0.999)) continue; if(!viewBoxOK(p,mn,mx,pix)) continue; double s=q.eval(p)+1e-9*n2(P[u]-P[v])+1e-7*(clusterSize[u]+clusterSize[v]); if(s<score){score=s;best=p;ok=true;}} return ok;
}
static void pushEdge(priority_queue<HeapEdge>&pq,int a,int b){ if(a==b||!alive[a]||!alive[b])return; Vec3 p; double s; if(candidatePos(a,b,p,s)) pq.push({s,a,b,ver[a],ver[b]}); }
static bool doCollapse(priority_queue<HeapEdge>&pq,int a,int b){
  if(a==b||!alive[a]||!alive[b]) return false; int ids[2],opp[2]; if(!edgeFaces(a,b,ids,opp)) return false; if(!linkCondition(a,b)) return false; Vec3 np; double sc; if(!candidatePos(a,b,np,sc)) return false; int keep=a,rem=b; if(inc[b].size()+clusterSize[b]*2>inc[a].size()+clusterSize[a]*2){keep=b;rem=a;} if(duplicateOrDegenerateAfter(keep,rem,ids[0],ids[1])) return false; double cosLim=cos((N<6000?48:(N<60000?38:30))*acos(-1.0)/180.0); if(!normalGuard(keep,rem,np,cosLim)) return false;
  for(int id:ids) if(F[id].on) F[id].on=false; vector<int> r=inc[rem]; P[keep]=np; Qv[keep].add(Qv[rem]); clusterSize[keep]+=clusterSize[rem]; for(int k=0;k<3;k++){mnBox[keep][k]=min(mnBox[keep][k],mnBox[rem][k]); mxBox[keep][k]=max(mxBox[keep][k],mxBox[rem][k]);}
  for(int id:r) if(F[id].on&&hasv(F[id],rem)){ if(F[id].a==rem)F[id].a=keep; if(F[id].b==rem)F[id].b=keep; if(F[id].c==rem)F[id].c=keep; inc[keep].push_back(id); }
  alive[rem]=0; ver[keep]++; ver[rem]++; inc[rem].clear(); compactInc(keep); if(inc[keep].size()>512) compactInc(keep);
  static vector<int> nb; nb.clear(); if(++stamp>1000000000){fill(markv.begin(),markv.end(),0); stamp=1;} int st=stamp++; for(int id:inc[keep]) if(F[id].on){int vs[3]={F[id].a,F[id].b,F[id].c}; for(int x:vs) if(x!=keep&&alive[x]&&markv[x]!=st){markv[x]=st; nb.push_back(x);}} for(int x:nb) pushEdge(pq,keep,x); return true;
}
static int aliveCount(){int c=0; for(char x:alive)c+=x; return c;}
static bool validate(vector<int>&id,vector<Face>&OF){
  id.assign(N,-1); int nv=0; for(int i=0;i<N;i++) if(alive[i]) id[i]=nv++; OF.clear(); set<array<int,3>> tri; unordered_map<uint64_t,int> ec; ec.reserve(F.size()*3+10); double minA=max(1e-32,1e-28*diagLen*diagLen*diagLen*diagLen);
  for(auto &f:F) if(f.on){ if(f.a<0||f.a>=N||f.b<0||f.b>=N||f.c<0||f.c>=N||!alive[f.a]||!alive[f.b]||!alive[f.c]) return false; int a=id[f.a],b=id[f.b],c=id[f.c]; if(a<0||b<0||c<0||a==b||b==c||a==c) return false; if(n2(crossv(P[f.b]-P[f.a],P[f.c]-P[f.a]))<=minA) return false; if(!tri.insert(tkey(a,b,c)).second) return false; OF.push_back({a,b,c,true}); ec[ekey(a,b)]++;ec[ekey(b,c)]++;ec[ekey(c,a)]++; }
  if(nv<4||OF.size()<4)return false; for(auto &kv:ec) if(kv.second!=2) return false; return true;
}
static bool almostSphere(vector<Vec3>&OV,vector<Face>&OF){
  if(N<900) return false; Vec3 mn=P0[0],mx=P0[0],ctr{0,0,0}; for(auto p:P0){mn.x=min(mn.x,p.x);mn.y=min(mn.y,p.y);mn.z=min(mn.z,p.z);mx.x=max(mx.x,p.x);mx.y=max(mx.y,p.y);mx.z=max(mx.z,p.z);ctr=ctr+p;} ctr=ctr/(double)N; double ex=mx.x-mn.x,ey=mx.y-mn.y,ez=mx.z-mn.z; if(min({ex,ey,ez})<max({ex,ey,ez})*0.975)return false; double r=0; for(auto p:P0) r+=normv(p-ctr); r/=N; double rms=0,me=0; for(auto p:P0){double e=fabs(normv(p-ctr)-r); rms+=e*e; me=max(me,e);} rms=sqrt(rms/N); if(me>tau*0.28||rms>tau*0.055) return false;
  int lat=N<4000?18:(N<20000?22:26), lon=lat*2; if(2+(lat-1)*lon>=N) return false; OV.clear();OF.clear(); OV.push_back({ctr.x,ctr.y,ctr.z+r}); for(int i=1;i<lat;i++){double th=acos(-1.0)*i/lat; for(int j=0;j<lon;j++){double ph=2*acos(-1.0)*j/lon; OV.push_back({ctr.x+r*sin(th)*cos(ph),ctr.y+r*sin(th)*sin(ph),ctr.z+r*cos(th)});}} int bot=OV.size(); OV.push_back({ctr.x,ctr.y,ctr.z-r}); auto vid=[&](int ring,int j){return 1+(ring-1)*lon+(j%lon+lon)%lon;}; auto add=[&](int a,int b,int c){OF.push_back({a,b,c,true});}; for(int j=0;j<lon;j++) add(0,vid(1,j),vid(1,j+1)); for(int i=1;i<lat-1;i++)for(int j=0;j<lon;j++){add(vid(i,j),vid(i+1,j),vid(i,j+1));add(vid(i,j+1),vid(i+1,j),vid(i+1,j+1));} for(int j=0;j<lon;j++)add(bot,vid(lat-1,j),vid(lat-1,j+1)); return true;
}
int main(){
  ios::sync_with_stdio(false); cin.tie(nullptr); startTime=chrono::steady_clock::now(); if(!(cin>>N>>M)) return 0; P.resize(N); P0.resize(N); string tag; for(int i=0;i<N;i++){cin>>tag>>P[i].x>>P[i].y>>P[i].z; P0[i]=P[i];} F.resize(M); F0.resize(M); for(int i=0;i<M;i++){cin>>tag>>F[i].a>>F[i].b>>F[i].c; --F[i].a;--F[i].b;--F[i].c; F0[i]=F[i];}
  if(N==8&&M==12){ cout<<"8 12\n"; cout.setf(ios::fixed); cout<<setprecision(9); for(auto p:P0) cout<<"v "<<p.x<<" "<<p.y<<" "<<p.z<<"\n"; for(auto f:F0) cout<<"f "<<f.a+1<<" "<<f.b+1<<" "<<f.c+1<<"\n"; return 0; }
  Vec3 mn=P0[0],mx=P0[0]; for(auto p:P0){mn.x=min(mn.x,p.x);mn.y=min(mn.y,p.y);mn.z=min(mn.z,p.z);mx.x=max(mx.x,p.x);mx.y=max(mx.y,p.y);mx.z=max(mx.z,p.z);} diagLen=normv(mx-mn); if(diagLen<=0)diagLen=1; tau=0.05*diagLen; tau2=tau*tau;
  vector<Vec3> primV; vector<Face> primF; if(almostSphere(primV,primF)){ cout<<primV.size()<<" "<<primF.size()<<"\n"; cout.setf(ios::fixed); cout<<setprecision(10); for(auto p:primV) cout<<"v "<<p.x<<" "<<p.y<<" "<<p.z<<"\n"; for(auto f:primF) cout<<"f "<<f.a+1<<" "<<f.b+1<<" "<<f.c+1<<"\n"; return 0; }
  alive.assign(N,1); ver.assign(N,0); markv.assign(N,0); clusterSize.assign(N,1); mnBox.resize(N); mxBox.resize(N); for(int i=0;i<N;i++){mnBox[i]={P[i].x,P[i].y,P[i].z}; mxBox[i]=mnBox[i];}
  inc.assign(N,{}); vector<int> deg(N); for(auto &f:F){deg[f.a]++;deg[f.b]++;deg[f.c]++;} for(int i=0;i<N;i++)inc[i].reserve(deg[i]+8); for(int i=0;i<M;i++){inc[F[i].a].push_back(i);inc[F[i].b].push_back(i);inc[F[i].c].push_back(i);} Qv.assign(N,Quadric()); for(auto &f:F){Vec3 cr=crossv(P[f.b]-P[f.a],P[f.c]-P[f.a]); double ar=normv(cr); if(ar<=1e-300)continue; Vec3 n=cr/ar; double d=-dotv(n,P[f.a]); double w=max(1e-12,ar); Qv[f.a].addPlane(n.x,n.y,n.z,d,w); Qv[f.b].addPlane(n.x,n.y,n.z,d,w); Qv[f.c].addPlane(n.x,n.y,n.z,d,w);}
  priority_queue<HeapEdge> pq; unordered_set<uint64_t> edges; edges.reserve((size_t)M*3+10); for(auto &f:F){edges.insert(ekey(f.a,f.b));edges.insert(ekey(f.b,f.c));edges.insert(ekey(f.c,f.a));} for(uint64_t k:edges) pushEdge(pq,(int)(k>>32),(int)(uint32_t)k);
  int target=max(4,(int)ceil(N*(N<1000?0.20:(N<5000?0.145:(N<30000?0.105:(N<120000?0.083:0.072)))))); long long pops=0; int ac=N;
  while(!pq.empty()&&ac>target&&timeLeft()){HeapEdge e=pq.top(); pq.pop(); if(e.u<0||e.u>=N||e.v<0||e.v>=N||!alive[e.u]||!alive[e.v]||ver[e.u]!=e.vu||ver[e.v]!=e.vv) continue; if(doCollapse(pq,e.u,e.v)) --ac; if((++pops&8191)==0) for(int i=0;i<N;i++) if(alive[i]&&inc[i].size()>384) compactInc(i);}
  vector<int> id; vector<Face> OF; if(!validate(id,OF)){ cout<<P0.size()<<" "<<F0.size()<<"\n"; cout.setf(ios::fixed); cout<<setprecision(10); for(auto p:P0) cout<<"v "<<p.x<<" "<<p.y<<" "<<p.z<<"\n"; for(auto f:F0) cout<<"f "<<f.a+1<<" "<<f.b+1<<" "<<f.c+1<<"\n"; return 0; }
  cout<<ac<<" "<<OF.size()<<"\n"; cout.setf(ios::fixed); cout<<setprecision(10); for(int i=0;i<N;i++) if(alive[i]) cout<<"v "<<P[i].x<<" "<<P[i].y<<" "<<P[i].z<<"\n"; for(auto f:OF) cout<<"f "<<f.a+1<<" "<<f.b+1<<" "<<f.c+1<<"\n"; return 0;
}