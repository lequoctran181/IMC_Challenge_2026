#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <functional>
#include <map>
#include <numeric>
#include <queue>
#include <set>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>
using namespace std;

struct Vec3{double x,y,z;};
static inline Vec3 operator+(const Vec3&a,const Vec3&b){return {a.x+b.x,a.y+b.y,a.z+b.z};}
static inline Vec3 operator-(const Vec3&a,const Vec3&b){return {a.x-b.x,a.y-b.y,a.z-b.z};}
static inline Vec3 operator*(const Vec3&a,double s){return {a.x*s,a.y*s,a.z*s};}
static inline double dotv(const Vec3&a,const Vec3&b){return a.x*b.x+a.y*b.y+a.z*b.z;}
static inline Vec3 crossv(const Vec3&a,const Vec3&b){return {a.y*b.z-a.z*b.y,a.z*b.x-a.x*b.z,a.x*b.y-a.y*b.x};}
static inline double n2(const Vec3&a){return dotv(a,a);} static inline double nrm(const Vec3&a){return sqrt(n2(a));}
struct Face{int v[3];};

struct FastInput{
    vector<char> b; char* p;
    FastInput(){b.reserve(1<<27); char c[1<<16]; size_t n; while((n=fread(c,1,sizeof(c),stdin))>0)b.insert(b.end(),c,c+n); b.push_back(0); p=b.data();}
    inline void ws(){while(*p==' '||*p=='\n'||*p=='\r'||*p=='\t')++p;}
    long nextLong(){ws(); return strtol(p,&p,10);} double nextDouble(){ws(); return strtod(p,&p);} char nextChar(){ws(); return *p++;}
};

static int N,M,liveFaces; static double diagLen=1,EPS=1;
static vector<Vec3>P,OrigP; static vector<Face>F,OrigF;
static vector<unsigned char>Vlive,Flive; static vector<double>errv; static vector<vector<int>>Adj;
static chrono::steady_clock::time_point T0;
static inline double elapsed(){return chrono::duration<double>(chrono::steady_clock::now()-T0).count();}
static inline bool hasv(int fid,int x){const Face&f=F[fid];return f.v[0]==x||f.v[1]==x||f.v[2]==x;}
static inline int otherInFace(int fid,int a,int b){const Face&f=F[fid];for(int i=0;i<3;i++){int x=f.v[i];if(x!=a&&x!=b)return x;}return -1;}
static inline uint64_t edgeKey(int a,int b){if(a>b)swap(a,b);return (uint64_t)(uint32_t)a<<32 | (uint32_t)b;}
static inline uint64_t faceKey(int a,int b,int c){if(a>b)swap(a,b);if(b>c)swap(b,c);if(a>b)swap(a,b);return ((uint64_t)a<<42)^((uint64_t)b<<21)^(uint64_t)c;}
static Vec3 faceNormalRaw(const Face&f){return crossv(P[f.v[1]]-P[f.v[0]],P[f.v[2]]-P[f.v[0]]);} 

static void load(){FastInput in;N=(int)in.nextLong();M=(int)in.nextLong();P.resize(N);OrigP.resize(N);Vlive.assign(N,1);errv.assign(N,0);Vec3 mn{1e100,1e100,1e100},mx{-1e100,-1e100,-1e100};for(int i=0;i<N;i++){in.nextChar();P[i].x=in.nextDouble();P[i].y=in.nextDouble();P[i].z=in.nextDouble();OrigP[i]=P[i];mn.x=min(mn.x,P[i].x);mn.y=min(mn.y,P[i].y);mn.z=min(mn.z,P[i].z);mx.x=max(mx.x,P[i].x);mx.y=max(mx.y,P[i].y);mx.z=max(mx.z,P[i].z);}diagLen=nrm(mx-mn);if(!(diagLen>0))diagLen=1;EPS=.049*diagLen;F.resize(M);OrigF.resize(M);Flive.assign(M,1);vector<int>deg(N);for(int i=0;i<M;i++){in.nextChar();F[i].v[0]=(int)in.nextLong()-1;F[i].v[1]=(int)in.nextLong()-1;F[i].v[2]=(int)in.nextLong()-1;OrigF[i]=F[i];deg[F[i].v[0]]++;deg[F[i].v[1]]++;deg[F[i].v[2]]++;}Adj.resize(N);for(int i=0;i<N;i++)Adj[i].reserve(deg[i]);for(int i=0;i<M;i++){Adj[F[i].v[0]].push_back(i);Adj[F[i].v[1]].push_back(i);Adj[F[i].v[2]].push_back(i);}liveFaces=M;}

static void restoreOriginal(){P=OrigP;F=OrigF;Vlive.assign(N,1);Flive.assign(M,1);errv.assign(N,0);vector<int>deg(N);for(const auto&f:F){deg[f.v[0]]++;deg[f.v[1]]++;deg[f.v[2]]++;}Adj.assign(N,{});for(int i=0;i<N;i++)Adj[i].reserve(deg[i]);for(int i=0;i<M;i++){Adj[F[i].v[0]].push_back(i);Adj[F[i].v[1]].push_back(i);Adj[F[i].v[2]].push_back(i);}liveFaces=M;}

static bool edgeFaces(int u,int v,int by[2],int op[2]){int cnt=0;const vector<int>&A=(Adj[u].size()<Adj[v].size()?Adj[u]:Adj[v]);for(int fid:A){if(!Flive[fid])continue;if(!hasv(fid,u)||!hasv(fid,v))continue;if(cnt>=2)return false;by[cnt]=fid;op[cnt]=otherInFace(fid,u,v);cnt++;}return cnt==2&&op[0]>=0&&op[1]>=0&&op[0]!=op[1];}

static vector<int>markv,mark2; static int stamp1=1,stamp2=1;
static bool linkOK(int u,int v,const int op[2]){if(markv.empty()){markv.assign(N,0);mark2.assign(N,0);} if(++stamp1==INT_MAX){fill(markv.begin(),markv.end(),0);stamp1=1;} if(++stamp2==INT_MAX){fill(mark2.begin(),mark2.end(),0);stamp2=1;}for(int fid:Adj[u])if(Flive[fid]&&hasv(fid,u)){const Face&f=F[fid];for(int i=0;i<3;i++){int x=f.v[i];if(x!=u&&x!=v)markv[x]=stamp1;}}int cnt=0;for(int fid:Adj[v])if(Flive[fid]&&hasv(fid,v)){const Face&f=F[fid];for(int i=0;i<3;i++){int x=f.v[i];if(x==u||x==v||markv[x]!=stamp1)continue;if(x!=op[0]&&x!=op[1])return false;if(mark2[x]!=stamp2){mark2[x]=stamp2;cnt++;}}}return cnt==2&&mark2[op[0]]==stamp2&&mark2[op[1]]==stamp2;}

static bool sameFace(const Face&f,int a,int b,int c){return faceKey(f.v[0],f.v[1],f.v[2])==faceKey(a,b,c);} 
static bool wouldDuplicate(int rem,const vector<int>&skip,const vector<Face>&nf){vector<uint64_t>keys;keys.reserve(nf.size());for(const auto&g:nf){uint64_t k=faceKey(g.v[0],g.v[1],g.v[2]);if(find(keys.begin(),keys.end(),k)!=keys.end())return true;keys.push_back(k);}for(const auto&g:nf){int a=g.v[0],b=g.v[1],c=g.v[2];int ix=a;if(Adj[b].size()<Adj[ix].size())ix=b;if(Adj[c].size()<Adj[ix].size())ix=c;for(int fid:Adj[ix]){if(!Flive[fid])continue;if(find(skip.begin(),skip.end(),fid)!=skip.end())continue;if(hasv(fid,rem))continue;if(sameFace(F[fid],a,b,c))return true;}}return false;}

struct Eval{bool ok=false;double cost=1e100,newErr=0;};
struct Params{double rad,plane,cosang,minArea;};
static Eval evalCollapse(int keep,int rem,const int by[2],const Params&pa){Eval r;if(!Vlive[keep]||!Vlive[rem])return r;double d=nrm(P[keep]-P[rem]);r.newErr=max(errv[keep],errv[rem]+d);if(r.newErr>pa.rad||r.newErr>EPS)return r;vector<int>skip={by[0],by[1]};vector<Face>nf;nf.reserve(Adj[rem].size());double worstPlane=0,worstCos=0;int changed=0;for(int fid:Adj[rem]){if(!Flive[fid]||!hasv(fid,rem))continue;if(fid==by[0]||fid==by[1])continue;if(hasv(fid,keep))return r;Face old=F[fid],g=old;for(int i=0;i<3;i++)if(g.v[i]==rem)g.v[i]=keep;if(g.v[0]==g.v[1]||g.v[0]==g.v[2]||g.v[1]==g.v[2])return r;Vec3 co=faceNormalRaw(old), cn=faceNormalRaw(g);double lo=nrm(co),ln=nrm(cn);if(!(lo>pa.minArea&&ln>pa.minArea))return r;double cs=dotv(co,cn)/(lo*ln);if(cs<pa.cosang)return r;Vec3 un=co*(1.0/lo);double pl=fabs(dotv(un,P[keep]-P[old.v[0]]));if(pl>pa.plane)return r;worstPlane=max(worstPlane,pl);worstCos=max(worstCos,1-cs);nf.push_back(g);changed++;}
if(changed==0)return r;if(wouldDuplicate(rem,skip,nf))return r;r.ok=true;r.cost=r.newErr/(pa.rad+1e-30)+50*worstPlane/(pa.plane+1e-30)+500*worstCos+1e-6*Adj[rem].size();return r;}

static void cleanAdj(int v){vector<int>z;z.reserve(Adj[v].size());for(int fid:Adj[v])if(fid>=0&&fid<(int)F.size()&&Flive[fid]&&hasv(fid,v))z.push_back(fid);sort(z.begin(),z.end());z.erase(unique(z.begin(),z.end()),z.end());Adj[v].swap(z);} 
static void mergeAdj(int keep,int rem){vector<int>z;z.reserve(Adj[keep].size()+Adj[rem].size());for(int fid:Adj[keep])if(fid>=0&&fid<(int)F.size()&&Flive[fid]&&hasv(fid,keep))z.push_back(fid);for(int fid:Adj[rem])if(fid>=0&&fid<(int)F.size()&&Flive[fid]&&hasv(fid,keep))z.push_back(fid);sort(z.begin(),z.end());z.erase(unique(z.begin(),z.end()),z.end());Adj[keep].swap(z);Adj[rem].clear();}
static bool collapseEdge(int a,int b,const Params&pa){if(a==b||!Vlive[a]||!Vlive[b])return false;int by[2]={-1,-1},op[2]={-1,-1};if(!edgeFaces(a,b,by,op))return false;if(!linkOK(a,b,op))return false;Eval ab=evalCollapse(a,b,by,pa),ba=evalCollapse(b,a,by,pa);if(!ab.ok&&!ba.ok)return false;int keep=a,rem=b;Eval e=ab;if(ba.ok&&(!ab.ok||ba.cost<ab.cost)){keep=b;rem=a;e=ba;}for(int i=0;i<2;i++)if(Flive[by[i]]){Flive[by[i]]=0;liveFaces--;}for(int fid:Adj[rem]){if(!Flive[fid]||!hasv(fid,rem))continue;Face&f=F[fid];for(int i=0;i<3;i++)if(f.v[i]==rem)f.v[i]=keep;}Vlive[rem]=0;errv[keep]=e.newErr;mergeAdj(keep,rem);return true;}

static int countLiveVertices(){int c=0;for(unsigned char x:Vlive)c+=x?1:0;return c;}
static bool passCollapse(const Params&pa,bool stopNoChange){int changed=0;for(int fid=0;fid<(int)F.size();fid++){if((fid&4095)==0&&elapsed()>19.2)break;if(!Flive[fid])continue;Face f=F[fid];if(!Vlive[f.v[0]]||!Vlive[f.v[1]]||!Vlive[f.v[2]])continue;array<array<double,3>,3> e{{{n2(P[f.v[0]]-P[f.v[1]]),(double)f.v[0],(double)f.v[1]},{n2(P[f.v[1]]-P[f.v[2]]),(double)f.v[1],(double)f.v[2]},{n2(P[f.v[2]]-P[f.v[0]]),(double)f.v[2],(double)f.v[0]}}};sort(e.begin(),e.end(),[](auto&x,auto&y){return x[0]<y[0];});for(auto &q:e){if(collapseEdge((int)q[1],(int)q[2],pa)){changed++;break;}}}
return changed>0||!stopNoChange;}

struct FullState{vector<Vec3>P;vector<Face>F;vector<unsigned char>Vlive,Flive;vector<double>errv;vector<vector<int>>Adj;int liveFaces;};
static FullState saveState(){return {P,F,Vlive,Flive,errv,Adj,liveFaces};}
static void loadState(const FullState&s){P=s.P;F=s.F;Vlive=s.Vlive;Flive=s.Flive;errv=s.errv;Adj=s.Adj;liveFaces=s.liveFaces;}

static int edgeCountOutside(int a,int b,const vector<int>&skip){const vector<int>&A=(Adj[a].size()<Adj[b].size()?Adj[a]:Adj[b]);int c=0;for(int fid:A){if(!Flive[fid])continue;if(find(skip.begin(),skip.end(),fid)!=skip.end())continue;if(hasv(fid,a)&&hasv(fid,b))c++;if(c>2)break;}return c;}
static bool faceExistsOutside(const Face&g,const vector<int>&skip){int a=g.v[0],b=g.v[1],c=g.v[2];int ix=a;if(Adj[b].size()<Adj[ix].size())ix=b;if(Adj[c].size()<Adj[ix].size())ix=c;for(int fid:Adj[ix])if(Flive[fid]&&find(skip.begin(),skip.end(),fid)==skip.end()&&sameFace(F[fid],a,b,c))return true;return false;}
struct Patch{int v=-1,best=-1;vector<int>inc,ring;vector<Face>nf;double score=0,newErr=0;};
static bool inspectPatch(int v,Patch&pc){pc=Patch();pc.v=v;if(!Vlive[v])return false;for(int fid:Adj[v])if(Flive[fid]&&hasv(fid,v))pc.inc.push_back(fid);int k=(int)pc.inc.size();if(k<3||k>7)return false;vector<int>nb;vector<pair<int,int>>be;Vec3 nsum{0,0,0};double minA=max(1e-30,1e-24*diagLen*diagLen);for(int fid:pc.inc){const Face&f=F[fid];int o[2],t=0;for(int i=0;i<3;i++)if(f.v[i]!=v)o[t++]=f.v[i];if(t!=2||o[0]==o[1]||!Vlive[o[0]]||!Vlive[o[1]])return false;for(int j=0;j<2;j++)if(find(nb.begin(),nb.end(),o[j])==nb.end())nb.push_back(o[j]);int a=o[0],b=o[1];if(a>b)swap(a,b);if(find(be.begin(),be.end(),make_pair(a,b))!=be.end())return false;be.push_back({a,b});Vec3 cr=faceNormalRaw(f);double l=nrm(cr);if(!(l>0))return false;nsum=nsum+cr*(1.0/l);}if((int)nb.size()!=k)return false;double nl=nrm(nsum);if(!(nl>1e-12))return false;Vec3 normal=nsum*(1.0/nl);double oldCos=cos(12.0*acos(-1)/180.0);for(int fid:pc.inc){Vec3 cr=faceNormalRaw(F[fid]);double l=nrm(cr);if(!(l>0)||dotv(cr*(1.0/l),normal)<oldCos)return false;}vector<vector<int>>g(k);auto idx=[&](int x){for(int i=0;i<k;i++)if(nb[i]==x)return i;return -1;};for(auto e:be){int a=idx(e.first),b=idx(e.second);if(a<0||b<0||a==b)return false;g[a].push_back(b);g[b].push_back(a);}for(int i=0;i<k;i++)if(g[i].size()!=2)return false;vector<int>ord;int cur=0,pre=-1;for(int s=0;s<k;s++){ord.push_back(cur);int nx=(g[cur][0]==pre?g[cur][1]:g[cur][0]);pre=cur;cur=nx;}if(cur!=0)return false;for(int id:ord)pc.ring.push_back(nb[id]);double bestD=1e100;for(int u:pc.ring){double d=nrm(P[v]-P[u]);if(errv[v]+d<bestD){bestD=errv[v]+d;pc.best=u;}}if(pc.best<0||bestD>EPS*.995)return false;pc.newErr=max(errv[pc.best],bestD);auto isBoundary=[&](int a,int b){if(a>b)swap(a,b);return find(be.begin(),be.end(),make_pair(a,b))!=be.end();};double newCos=cos(20.0*acos(-1)/180.0);vector<uint64_t>newKeys;for(int i=1;i+1<k;i++){Face h{{pc.ring[0],pc.ring[i],pc.ring[i+1]}};Vec3 cr=faceNormalRaw(h);double l=nrm(cr);if(!(l>minA))return false;if(dotv(cr,normal)<0){swap(h.v[1],h.v[2]);cr=faceNormalRaw(h);l=nrm(cr);}if(!(l>minA)||dotv(cr*(1.0/l),normal)<newCos)return false;if(faceExistsOutside(h,pc.inc))return false;uint64_t fk=faceKey(h.v[0],h.v[1],h.v[2]);if(find(newKeys.begin(),newKeys.end(),fk)!=newKeys.end())return false;newKeys.push_back(fk);int ee[3][2]={{h.v[0],h.v[1]},{h.v[1],h.v[2]},{h.v[2],h.v[0]}};for(int q=0;q<3;q++){int a=ee[q][0],b=ee[q][1];int cnt=edgeCountOutside(a,b,pc.inc);if(isBoundary(a,b)){if(cnt!=1)return false;}else if(cnt!=0)return false;}pc.nf.push_back(h);}pc.score=pc.newErr/EPS+0.001*k;return true;}
static bool applyPatch(const Patch&pc){if(!Vlive[pc.v])return false;Patch p;if(!inspectPatch(pc.v,p))return false;if(p.ring!=pc.ring||p.best!=pc.best)return false;for(int fid:p.inc)if(Flive[fid]){Flive[fid]=0;liveFaces--;}Vlive[p.v]=0;Adj[p.v].clear();errv[p.best]=p.newErr;for(const Face&h:p.nf){int id=(int)F.size();F.push_back(h);Flive.push_back(1);liveFaces++;for(int i=0;i<3;i++)Adj[h.v[i]].push_back(id);}return true;}

static bool validateMesh(){vector<uint64_t>edges,facesK;edges.reserve((size_t)liveFaces*3);facesK.reserve(liveFaces);double minA=max(1e-30,1e-24*diagLen*diagLen);int lf=0;for(int i=0;i<(int)F.size();i++)if(Flive[i]){lf++;Face f=F[i];for(int j=0;j<3;j++)if(f.v[j]<0||f.v[j]>=N||!Vlive[f.v[j]])return false;if(f.v[0]==f.v[1]||f.v[0]==f.v[2]||f.v[1]==f.v[2])return false;if(n2(faceNormalRaw(f))<=minA)return false;facesK.push_back(faceKey(f.v[0],f.v[1],f.v[2]));edges.push_back(edgeKey(f.v[0],f.v[1]));edges.push_back(edgeKey(f.v[1],f.v[2]));edges.push_back(edgeKey(f.v[2],f.v[0]));}if(lf==0)return false;sort(facesK.begin(),facesK.end());if(adjacent_find(facesK.begin(),facesK.end())!=facesK.end())return false;sort(edges.begin(),edges.end());for(size_t i=0;i<edges.size();){size_t j=i+1;while(j<edges.size()&&edges[j]==edges[i])j++;if(j-i!=2)return false;i=j;}return countLiveVertices()<=N;}

static bool patchBandPass(){if(!(N>=49843&&N<50625)||elapsed()>18.0)return false;FullState S=saveState();int before=countLiveVertices();vector<Patch>C;C.reserve(4000);for(int v=0;v<N;v++){if((v&511)==0&&elapsed()>18.55)break;Patch p;if(inspectPatch(v,p))C.push_back(p);}if(C.empty()){loadState(S);return false;}sort(C.begin(),C.end(),[](const Patch&a,const Patch&b){return a.score<b.score;});int cap=max(60,min(1800,before*35/1000)),rm=0;for(const Patch&p:C){if(rm>=cap||elapsed()>19.05)break;if(applyPatch(p))rm++;}int after=countLiveVertices();if(rm<5||after>=before||!validateMesh()){loadState(S);return false;}return true;}

static void simplify(){T0=chrono::steady_clock::now();double d=diagLen;vector<Params> ps;auto add=[&](double r,double pl,double deg){ps.push_back({min(r*d,EPS*.995),pl*d,cos(deg*acos(-1)/180.0),max(1e-30,1e-24*d*d)});};add(.004,.001,1.0);add(.008,.002,1.8);add(.012,.0032,2.8);if(N<1000){add(.018,.0045,4.0);add(.024,.0060,5.4);}else if(N<60000){add(.021,.0058,4.8);add(.027,.0075,6.2);add(.030,.0085,6.8);}else if(N<180000){add(.024,.0068,5.6);add(.031,.0090,7.2);}else{add(.026,.0078,6.2);add(.034,.0105,8.0);}for(size_t i=0;i<ps.size()&&elapsed()<18.9;i++){bool ch=passCollapse(ps[i],true);if(!ch&&i+1>3)break;}patchBandPass();if(!validateMesh())restoreOriginal();}

static void output(){vector<int>id(N,-1);int vn=0,fn=0;for(int i=0;i<N;i++)if(Vlive[i])id[i]=vn++;for(int i=0;i<(int)F.size();i++)if(Flive[i])fn++;string out;out.reserve((size_t)vn*38+(size_t)fn*24+64);char buf[128];out.append(buf,snprintf(buf,sizeof(buf),"%d %d\n",vn,fn));for(int i=0;i<N;i++)if(Vlive[i])out.append(buf,snprintf(buf,sizeof(buf),"v %.10g %.10g %.10g\n",P[i].x,P[i].y,P[i].z));for(int i=0;i<(int)F.size();i++)if(Flive[i])out.append(buf,snprintf(buf,sizeof(buf),"f %d %d %d\n",id[F[i].v[0]]+1,id[F[i].v[1]]+1,id[F[i].v[2]]+1));fwrite(out.data(),1,out.size(),stdout);} 
int main(){load();simplify();output();return 0;}
