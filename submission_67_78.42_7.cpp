#include <algorithm>
#include <array>
#include <chrono>
#include <utility>
#include <queue>
#include <cstdint>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <climits>
#include <string>
#include <vector>
using namespace std;
struct Vec3{
double x=0.0,y=0.0,z=0.0;
}
;
static inline Vec3 operator-(const Vec3&a,const Vec3&b){
return{
a.x-b.x,a.y-b.y,a.z-b.z}
;
}
static inline Vec3 operator+(const Vec3&a,const Vec3&b){
return{
a.x+b.x,a.y+b.y,a.z+b.z}
;
}
static inline Vec3 operator*(const Vec3&a,double s){
return{
a.x*s,a.y*s,a.z*s}
;
}
static inline double dot3(const Vec3&a,const Vec3&b){
return a.x*b.x+a.y*b.y+a.z*b.z;
}
static inline Vec3 cross3(const Vec3&a,const Vec3&b){
return{
a.y*b.z-a.z*b.y,a.z*b.x-a.x*b.z,a.x*b.y-a.y*b.x,}
;
}
static inline double norm2(const Vec3&a){
return dot3(a,a);
}
static inline double norm3(const Vec3&a){
return sqrt(norm2(a));
}
struct Face{
int v[3]{
}
;
}
;
struct FastInput{
vector<char>buf;
char*p=nullptr;
FastInput(){
buf.reserve(1<<27);
char chunk[1<<16];
size_t n=0;
while((n=fread(chunk,1,sizeof(chunk),stdin))>0){
buf.insert(buf.end(),chunk,chunk+n);
}
buf.push_back('\0');
p=buf.data();
}
inline void skip_ws(){
while(*p==' '||*p=='\n'||*p=='\r'||*p=='\t')++p;
}
long next_long(){
skip_ws();
return strtol(p,&p,10);
}
double next_double(){
skip_ws();
return strtod(p,&p);
}
char next_char(){
skip_ws();
return*p++;
}
}
;
static int N=0,M=0;
static vector<Vec3>P;
static vector<Vec3>originalP;
static vector<Face>faces;
static vector<Face>AR;
static vector<unsigned char>BU;
static vector<unsigned char>BR;
static vector<double>BF;
static vector<vector<int>>Y;
static int BE=0;
static vector<int>EI;
static vector<int>BZ;
static int IM=1;
static int GP=1;
static double CL=1.0;
static chrono::steady_clock::time_point HJ;
static bool AS=false;
static double BL=0.0;
static double AL=0.0;
static double Z=1.0;
static double AH=1.0;
static double BG=1.0;
static inline bool AC(int fid,int v){
const Face&f=faces[fid];
return f.v[0]==v||f.v[1]==v||f.v[2]==v;
}
static inline int HK(int fid,int a,int b){
const Face&f=faces[fid];
for(int i=0;
i<3;
++i){
int x=f.v[i];
if(x!=a&&x!=b)return x;
}
return-1;
}
static inline array<int,3>CW(int a,int b,int c){
array<int,3>s{
a,b,c}
;
sort(s.begin(),s.end());
return s;
}
static inline bool IA(int fid,int a,int b,int c){
const Face&f=faces[fid];
return CW(f.v[0],f.v[1],f.v[2])==CW(a,b,c);
}
static inline Vec3 CD(int a,int b,int c){
return cross3(P[b]-P[a],P[c]-P[a]);
}
static inline Vec3 FK(int id,int moved_id,const Vec3&moved_pos){
return id==moved_id?moved_pos:P[id];
}
static inline Vec3 HL(int a,int b,int c,int moved_id,const Vec3&moved_pos){
Vec3 pa=FK(a,moved_id,moved_pos);
Vec3 pb=FK(b,moved_id,moved_pos);
Vec3 pc=FK(c,moved_id,moved_pos);
return cross3(pb-pa,pc-pa);
}
static inline Vec3 II(int fid){
const Face&f=faces[fid];
return CD(f.v[0],f.v[1],f.v[2]);
}
static inline Vec3 GQ(int fid){
Vec3 cr=II(fid);
double len=norm3(cr);
if(!(len>0.0))return{
0.0,0.0,0.0}
;
return cr*(1.0/len);
}
static inline unsigned long long ED(int a,int b){
if(a>b)swap(a,b);
return(unsigned long long)(unsigned int)a<<32|(unsigned int)b;
}
static void JC(){
FastInput in;
N=(int)in.next_long();
M=(int)in.next_long();
P.resize(N);
originalP.resize(N);
BU.assign(N,1);
BF.assign(N,0.0);
Vec3 mn{
1e100,1e100,1e100}
;
Vec3 mx{
-1e100,-1e100,-1e100}
;
for(int i=0;
i<N;
++i){
(void)in.next_char();
P[i].x=in.next_double();
P[i].y=in.next_double();
P[i].z=in.next_double();
originalP[i]=P[i];
mn.x=min(mn.x,P[i].x);
mn.y=min(mn.y,P[i].y);
mn.z=min(mn.z,P[i].z);
mx.x=max(mx.x,P[i].x);
mx.y=max(mx.y,P[i].y);
mx.z=max(mx.z,P[i].z);
}
CL=norm3(mx-mn);
if(!(CL>0.0))CL=1.0;
faces.resize(M);
AR.resize(M);
BR.assign(M,1);
vector<int>deg(N,0);
for(int i=0;
i<M;
++i){
(void)in.next_char();
int a=(int)in.next_long()-1;
int b=(int)in.next_long()-1;
int c=(int)in.next_long()-1;
faces[i].v[0]=a;
faces[i].v[1]=b;
faces[i].v[2]=c;
AR[i]=faces[i];
++deg[a];
++deg[b];
++deg[c];
}
Y.resize(N);
for(int i=0;
i<N;
++i)Y[i].reserve(deg[i]);
for(int i=0;
i<M;
++i){
Y[faces[i].v[0]].push_back(i);
Y[faces[i].v[1]].push_back(i);
Y[faces[i].v[2]].push_back(i);
}
BE=M;
EI.assign(N,0);
BZ.assign(N,0);
}
static bool FE(int u,int v,int BY[2],int opp[2]){
int cnt=0;
for(int fid:Y[u]){
if(!BR[fid])continue;
if(!AC(fid,u)||!AC(fid,v))continue;
if(cnt>=2)return false;
BY[cnt]=fid;
opp[cnt]=HK(fid,u,v);
++cnt;
}
if(cnt!=2)return false;
if(opp[0]<0||opp[1]<0||opp[0]==opp[1])return false;
return true;
}
static bool ET(int u,int v,const int opp[2]){
if(++IM>2000000000){
fill(EI.begin(),EI.end(),0);
IM=1;
}
if(++GP>2000000000){
fill(BZ.begin(),BZ.end(),0);
GP=1;
}
for(int fid:Y[u]){
if(!BR[fid]||!AC(fid,u))continue;
const Face&f=faces[fid];
for(int i=0;
i<3;
++i){
int x=f.v[i];
if(x!=u&&x!=v)EI[x]=IM;
}
}
int IG=0;
for(int fid:Y[v]){
if(!BR[fid]||!AC(fid,v))continue;
const Face&f=faces[fid];
for(int i=0;
i<3;
++i){
int x=f.v[i];
if(x==u||x==v)continue;
if(EI[x]!=IM)continue;
if(x!=opp[0]&&x!=opp[1])return false;
if(BZ[x]!=GP){
BZ[x]=GP;
++IG;
}
}
}
return IG==2&&BZ[opp[0]]==GP&&BZ[opp[1]]==GP;
}
struct AE{
double AI=0.0;
double BB=0.0;
double BQ=1.0;
double W=0.0;
double AT=1.0;
double AJ=0.0;
bool AA=false;
}
;
struct JI{
bool ok=false;
bool IU=false;
double cost=1e100;
double CS=0.0;
Vec3 EM{
}
;
}
;
struct AP{
vector<Vec3>P;
vector<Face>faces;
vector<unsigned char>BU;
vector<unsigned char>BR;
vector<double>BF;
vector<vector<int>>Y;
int BE=0;
}
;
static AP AD(){
AP s;
s.P=P;
s.faces=faces;
s.BU=BU;
s.BR=BR;
s.BF=BF;
s.Y=Y;
s.BE=BE;
return s;
}
static void restore_state(const AP&s){
P=s.P;
faces=s.faces;
BU=s.BU;
BR=s.BR;
BF=s.BF;
Y=s.Y;
BE=s.BE;
}
static bool FF(int keep,int rem,int fid,int a,int b,int c,int JT,int JU){
int IX=keep;
if((int)Y[a].size()<(int)Y[IX].size())IX=a;
if((int)Y[b].size()<(int)Y[IX].size())IX=b;
if((int)Y[c].size()<(int)Y[IX].size())IX=c;
for(int JE:Y[IX]){
if(!BR[JE]||JE==fid||JE==JT||JE==JU)continue;
if(AC(JE,rem))continue;
if(IA(JE,a,b,c))return true;
}
return false;
}
static JI FT(int keep,int rem,const int BY[2],const AE&params){
JI res;
res.EM=P[keep];
const double JL=norm3(P[keep]-P[rem]);
res.CS=max(BF[keep],BF[rem]+JL);
double DL=0.0;
double CU=0.0;
double BI=0.0;
int HD=0;
bool IU=true;
for(int fid:Y[rem]){
if(!BR[fid]||!AC(fid,rem))continue;
if(fid==BY[0]||fid==BY[1])continue;
if(AC(fid,keep))return res;
Face f=faces[fid];
int BK[3]={
f.v[0],f.v[1],f.v[2]}
;
for(int i=0;
i<3;
++i){
if(BK[i]==rem)BK[i]=keep;
}
if(BK[0]==BK[1]||BK[0]==BK[2]||BK[1]==BK[2]){
return res;
}
Vec3 FC=CD(f.v[0],f.v[1],f.v[2]);
Vec3 GU=CD(BK[0],BK[1],BK[2]);
double DG=norm3(FC);
double EF=norm3(GU);
if(!(DG>params.AJ)||!(EF>params.AJ))return res;
if(EF<DG*1e-7)return res;
double ndot=dot3(FC,GU)/(DG*EF);
if(ndot>1.0)ndot=1.0;
if(ndot<-1.0)ndot=-1.0;
if(ndot<params.BQ)return res;
Vec3 IW=FC*(1.0/DG);
double FD=fabs(dot3(IW,P[keep]-P[f.v[0]]));
if(FD>params.BB&&res.CS>params.AI)return res;
if(FF(keep,rem,fid,BK[0],BK[1],BK[2],BY[0],BY[1])){
return res;
}
DL=max(DL,FD);
CU=max(CU,1.0-ndot);
BI=max(BI,max(0.0,1.0-EF/DG));
if(FD>params.W||ndot<params.AT)IU=false;
++HD;
}
if(HD==0)return res;
if(!IU){
if(res.CS>params.AI)return res;
if(DL>params.BB)return res;
}
res.ok=true;
res.IU=IU;
double HU=params.AI>0.0?res.CS/params.AI:0.0;
double IL=params.BB>0.0?DL/params.BB:0.0;
res.cost=(IU?-1000.0:0.0)+0.80*HU+0.80*IL+250.0*CU+0.02*BI+0.0005*HD;
return res;
}
static JI FP(int keep,int rem,const int BY[2],const AE&params,const Vec3&EM){
JI res;
res.EM=EM;
const double JM=norm3(EM-P[keep]);
const double JP=norm3(EM-P[rem]);
res.CS=max(BF[keep]+JM,BF[rem]+JP);
if(res.CS>params.AI)return res;
static vector<int>DS;
DS.clear();
DS.reserve(Y[keep].size()+Y[rem].size());
for(int fid:Y[keep]){
if(BR[fid]&&fid!=BY[0]&&fid!=BY[1]&&AC(fid,keep)){
DS.push_back(fid);
}
}
for(int fid:Y[rem]){
if(BR[fid]&&fid!=BY[0]&&fid!=BY[1]&&AC(fid,rem)){
DS.push_back(fid);
}
}
sort(DS.begin(),DS.end());
DS.erase(unique(DS.begin(),DS.end()),DS.end());
if(DS.empty())return res;
double DL=0.0;
double CU=0.0;
double BI=0.0;
int HD=0;
for(int fid:DS){
const bool JH=AC(fid,rem);
Face f=faces[fid];
int BK[3]={
f.v[0],f.v[1],f.v[2]}
;
for(int i=0;
i<3;
++i){
if(BK[i]==rem)BK[i]=keep;
}
if(BK[0]==BK[1]||BK[0]==BK[2]||BK[1]==BK[2]){
return res;
}
Vec3 FC=CD(f.v[0],f.v[1],f.v[2]);
Vec3 GU=HL(BK[0],BK[1],BK[2],keep,EM);
double DG=norm3(FC);
double EF=norm3(GU);
if(!(DG>params.AJ)||!(EF>params.AJ))return res;
if(EF<DG*1e-7)return res;
double ndot=dot3(FC,GU)/(DG*EF);
if(ndot>1.0)ndot=1.0;
if(ndot<-1.0)ndot=-1.0;
if(ndot<params.BQ)return res;
Vec3 IW=FC*(1.0/DG);
double FD=fabs(dot3(IW,EM-P[f.v[0]]));
if(FD>params.BB)return res;
if(JH&&FF(keep,rem,fid,BK[0],BK[1],BK[2],BY[0],BY[1])){
return res;
}
DL=max(DL,FD);
CU=max(CU,1.0-ndot);
BI=max(BI,max(0.0,1.0-EF/DG));
++HD;
}
if(HD==0)return res;
res.ok=true;
res.IU=false;
double HU=params.AI>0.0?res.CS/params.AI:0.0;
double IL=params.BB>0.0?DL/params.BB:0.0;
res.cost=0.35+0.95*HU+0.95*IL+300.0*CU+0.03*BI+0.0007*HD;
return res;
}
static void HM(int keep,int rem){
vector<int>merged;
merged.reserve(Y[keep].size()+Y[rem].size());
for(int fid:Y[keep]){
if(BR[fid]&&AC(fid,keep))merged.push_back(fid);
}
for(int fid:Y[rem]){
if(BR[fid]&&AC(fid,keep))merged.push_back(fid);
}
sort(merged.begin(),merged.end());
merged.erase(unique(merged.begin(),merged.end()),merged.end());
Y[keep].swap(merged);
vector<int>().swap(Y[rem]);
}
static void HN(int keep,int rem,const int BY[2],double CS,const Vec3&EM){
for(int i=0;
i<2;
++i){
if(BR[BY[i]]){
BR[BY[i]]=0;
--BE;
}
}
for(int fid:Y[rem]){
if(!BR[fid]||!AC(fid,rem))continue;
Face&f=faces[fid];
for(int i=0;
i<3;
++i){
if(f.v[i]==rem)f.v[i]=keep;
}
}
BU[rem]=0;
P[keep]=EM;
BF[keep]=CS;
HM(keep,rem);
}
static bool GD(int u,int v,const AE&params){
if(u==v||!BU[u]||!BU[v])return false;
int BY[2]={
-1,-1}
;
int opp[2]={
-1,-1}
;
if(!FE(u,v,BY,opp))return false;
if(!ET(u,v,opp))return false;
JI uv=FT(u,v,BY,params);
JI vu=FT(v,u,BY,params);
if(params.AA){
Vec3 mid=(P[u]+P[v])*0.5;
JI um=FP(u,v,BY,params,mid);
JI vm=FP(v,u,BY,params,mid);
if(um.ok&&(!uv.ok||um.cost<uv.cost))uv=um;
if(vm.ok&&(!vu.ok||vm.cost<vu.cost))vu=vm;
}
if(!uv.ok&&!vu.ok)return false;
if(vu.ok&&(!uv.ok||vu.cost<uv.cost)){
HN(v,u,BY,vu.CS,vu.EM);
}
else{
HN(u,v,BY,uv.CS,uv.EM);
}
return true;
}
static inline double elapsed_seconds(){
return chrono::duration<double>(chrono::steady_clock::now()-HJ).count();
}
static bool EZ(){
AS=false;
if(N<400||M<300)return false;
const int target_faces=40000;
const int stride=max(1,M/target_faces);
const int AW=120000;
const double smooth_cos=cos(10.0*acos(-1.0)/180.0);
const double coarse_cos=cos(30.0*acos(-1.0)/180.0);
const double sharp_cos=cos(22.0*acos(-1.0)/180.0);
const double very_sharp_cos=cos(45.0*acos(-1.0)/180.0);
int sampled=0;
int smooth=0;
int coarse=0;
int sharp=0;
int very_sharp=0;
int bad=0;
for(int fid=0;
fid<M&&sampled<AW;
fid+=stride){
const Face&f=faces[fid];
const int e[3][2]={
{
f.v[0],f.v[1]}
,{
f.v[1],f.v[2]}
,{
f.v[2],f.v[0]}
}
;
for(int k=0;
k<3&&sampled<AW;
++k){
int BY[2]={
-1,-1}
;
int opp[2]={
-1,-1}
;
if(!FE(e[k][0],e[k][1],BY,opp)){
++bad;
continue;
}
Vec3 n0=GQ(BY[0]);
Vec3 n1=GQ(BY[1]);
double ndot=dot3(n0,n1);
if(ndot>1.0)ndot=1.0;
if(ndot<-1.0)ndot=-1.0;
++sampled;
if(ndot>smooth_cos)++smooth;
if(ndot>coarse_cos)++coarse;
if(ndot<sharp_cos)++sharp;
if(ndot<very_sharp_cos)++very_sharp;
}
}
const int min_samples=N<3000?250:1000;
if(sampled<min_samples)return false;
const double smooth_ratio=(double)smooth/(double)sampled;
const double coarse_ratio=(double)coarse/(double)sampled;
const double sharp_ratio=(double)sharp/(double)sampled;
const double very_sharp_ratio=(double)very_sharp/(double)sampled;
const double bad_ratio=(double)bad/(double)(sampled+bad);
AS=true;
BL=smooth_ratio;
AL=coarse_ratio;
Z=sharp_ratio;
AH=very_sharp_ratio;
BG=bad_ratio;
const bool very_smooth=smooth_ratio>=0.985&&sharp_ratio<=0.010&&bad_ratio<=0.010;
const bool small_coarse_curved=N<3000&&coarse_ratio>=0.995&&sharp_ratio<=0.080&&very_sharp_ratio<=0.005&&bad_ratio<=0.010;
return very_smooth||small_coarse_curved;
}
struct CR{
vector<double>JJ;
vector<Vec3>JF;
vector<unsigned char>JN;
}
;
static inline void GX(const Vec3&p,int view,int res,double&u,double&v,double&z){
constexpr double D=2.5;
const double f=800.0*((double)res/1024.0);
const double c=0.5*(double)res;
double sx=0.0,sy=0.0;
if(view==0){
sx=p.y;
sy=p.z;
z=D-p.x;
}
else if(view==1){
sx=-p.y;
sy=p.z;
z=D+p.x;
}
else if(view==2){
sx=-p.x;
sy=p.z;
z=D-p.y;
}
else if(view==3){
sx=p.x;
sy=p.z;
z=D+p.y;
}
else if(view==4){
sx=p.x;
sy=p.y;
z=D-p.z;
}
else{
sx=-p.x;
sy=p.y;
z=D+p.z;
}
u=f*sx/z+c;
v=f*sy/z+c;
}
static void FU(CR&rm,int res,const Vec3&a,const Vec3&b,const Vec3&c,const Vec3&unit_normal,int view){
double x0,y0,z0,x1,y1,z1,x2,y2,z2;
GX(a,view,res,x0,y0,z0);
GX(b,view,res,x1,y1,z1);
GX(c,view,res,x2,y2,z2);
if(z0<=0.0||z1<=0.0||z2<=0.0)return;
int xmin=max(0,(int)floor(min(x0,min(x1,x2))-0.5));
int xmax=min(res-1,(int)ceil(max(x0,max(x1,x2))+0.5));
int ymin=max(0,(int)floor(min(y0,min(y1,y2))-0.5));
int ymax=min(res-1,(int)ceil(max(y0,max(y1,y2))+0.5));
if(xmin>xmax||ymin>ymax)return;
const double den=(y1-y2)*(x0-x2)+(x2-x1)*(y0-y2);
if(fabs(den)<1e-12)return;
for(int yy=ymin;
yy<=ymax;
++yy){
const double py=(double)yy+0.5;
for(int xx=xmin;
xx<=xmax;
++xx){
const double px=(double)xx+0.5;
const double w0=((y1-y2)*(px-x2)+(x2-x1)*(py-y2))/den;
const double w1=((y2-y0)*(px-x2)+(x0-x2)*(py-y2))/den;
const double w2=1.0-w0-w1;
if(w0<-1e-9||w1<-1e-9||w2<-1e-9)continue;
const double zp=1.0/(w0/z0+w1/z1+w2/z2);
const int idx=yy*res+xx;
if(zp<rm.JJ[idx]){
rm.JJ[idx]=zp;
rm.JF[idx]=unit_normal;
rm.JN[idx]=1;
}
}
}
}
static CR HO(int view,int res){
CR rm;
rm.JJ.assign(res*res,255.0);
rm.JF.assign(res*res,Vec3{
0.0,0.0,0.0}
);
rm.JN.assign(res*res,0);
for(int fid=0;
fid<M;
++fid){
const Face&f=AR[fid];
const Vec3&a=originalP[f.v[0]];
const Vec3&b=originalP[f.v[1]];
const Vec3&c=originalP[f.v[2]];
Vec3 cr=cross3(b-a,c-a);
double len=norm3(cr);
if(!(len>0.0))continue;
FU(rm,res,a,b,c,cr*(1.0/len),view);
}
return rm;
}
static CR IB(int view,int res){
CR rm;
rm.JJ.assign(res*res,255.0);
rm.JF.assign(res*res,Vec3{
0.0,0.0,0.0}
);
rm.JN.assign(res*res,0);
for(int fid=0;
fid<(int)faces.size();
++fid){
if(!BR[fid])continue;
const Face&f=faces[fid];
Vec3 cr=II(fid);
double len=norm3(cr);
if(!(len>0.0))continue;
FU(rm,res,P[f.v[0]],P[f.v[1]],P[f.v[2]],cr*(1.0/len),view);
}
return rm;
}
static inline double HP(const Vec3&n,int ch){
if(ch==0)return(n.x+1.0)*127.5;
if(ch==1)return(n.y+1.0)*127.5;
return(n.z+1.0)*127.5;
}
template<class Getter>static double FV(const CR&a,const CR&b,const vector<unsigned char>&fg,int res,int win,Getter get_value){
const int rad=win/2;
const double c1=(0.01*255.0)*(0.01*255.0);
const double c2=(0.03*255.0)*(0.03*255.0);
double total=0.0;
int count=0;
for(int y=0;
y<res;
++y){
for(int x=0;
x<res;
++x){
const int BS=y*res+x;
if(!fg[BS])continue;
double sx=0.0,sy=0.0,sxx=0.0,syy=0.0,sxy=0.0;
int n=0;
for(int dy=-rad;
dy<=rad;
++dy){
int yy=y+dy;
if(yy<0)yy=0;
if(yy>=res)yy=res-1;
for(int dx=-rad;
dx<=rad;
++dx){
int xx=x+dx;
if(xx<0)xx=0;
if(xx>=res)xx=res-1;
const int idx=yy*res+xx;
const double vx=get_value(a,idx);
const double vy=get_value(b,idx);
sx+=vx;
sy+=vy;
sxx+=vx*vx;
syy+=vy*vy;
sxy+=vx*vy;
++n;
}
}
const double inv=1.0/(double)n;
const double ux=sx*inv;
const double uy=sy*inv;
const double varx=max(0.0,sxx*inv-ux*ux);
const double vary=max(0.0,syy*inv-uy*uy);
const double cov=sxy*inv-ux*uy;
const double num=(2.0*ux*uy+c1)*(2.0*cov+c2);
const double den=(ux*ux+uy*uy+c1)*(varx+vary+c2);
total+=num/den;
++count;
}
}
return count>0?total/(double)count:1.0;
}
static double visual_proxy_score(int res){
constexpr int WIN=11;
double total=0.0;
for(int view=0;
view<6;
++view){
CR orig=HO(view,res);
CR simp=IB(view,res);
vector<unsigned char>fg(res*res,0);
for(int i=0;
i<res*res;
++i)fg[i]=(orig.JN[i]||simp.JN[i])?1:0;
double normal_score=0.0;
for(int ch=0;
ch<3;
++ch){
normal_score+=FV(orig,simp,fg,res,WIN,[ch](const CR&m,int idx){
return HP(m.JF[idx],ch);
}
);
}
normal_score/=3.0;
double depth_score=FV(orig,simp,fg,res,WIN,[](const CR&m,int idx){
return m.JJ[idx];
}
);
total+=0.5*normal_score+0.5*depth_score;
}
return total/6.0;
}
struct CE{
bool ok=false;
Vec3 BS{
}
;
Vec3 BP[3]{
}
;
double AK[3]{
}
;
double GF=1e100;
double EU=1e100;
double FA=1e100;
}
;
static inline Vec3 JS(const Vec3&a,double s){
return{
a.x/s,a.y/s,a.z/s}
;
}
static inline double GE(const Vec3&p,const Vec3&a){
return dot3(p,a);
}
static double BC(const Vec3&BS){
const int AW=100000;
const int stride=max(1,M/AW);
double total=0.0;
int sampled=0;
for(int fid=0;
fid<M;
fid+=stride){
const Face&f=AR[fid];
const Vec3&a=originalP[f.v[0]];
const Vec3&b=originalP[f.v[1]];
const Vec3&c=originalP[f.v[2]];
Vec3 cr=cross3(b-a,c-a);
Vec3 ctr=(a+b+c)*(1.0/3.0);
const double s=dot3(cr,ctr-BS);
if(fabs(s)>1e-18){
total+=s;
++sampled;
}
}
if(sampled==0)return 1.0;
return total>=0.0?1.0:-1.0;
}
static void BD(vector<Vec3>&verts,Face&f,const Vec3&BS,double CC){
const Vec3&a=verts[f.v[0]];
const Vec3&b=verts[f.v[1]];
const Vec3&c=verts[f.v[2]];
Vec3 cr=cross3(b-a,c-a);
Vec3 ctr=(a+b+c)*(1.0/3.0);
if(dot3(cr,ctr-BS)*CC<0.0)swap(f.v[1],f.v[2]);
}
static bool AF(const vector<Vec3>&X,const vector<Face>&AY){
if(X.empty()||AY.empty())return false;
if((int)X.size()>N)return false;
const double IY=max(1e-30,1e-24*CL*CL);
for(const Face&f:AY){
for(int k=0;
k<3;
++k){
if(f.v[k]<0||f.v[k]>=(int)X.size())return false;
}
if(f.v[0]==f.v[1]||f.v[0]==f.v[2]||f.v[1]==f.v[2])return false;
Vec3 cr=cross3(X[f.v[1]]-X[f.v[0]],X[f.v[2]]-X[f.v[0]]);
if(!(norm2(cr)>IY))return false;
}
P.assign(N,Vec3{
}
);
for(int i=0;
i<(int)X.size();
++i)P[i]=X[i];
BU.assign(N,0);
BF.assign(N,0.0);
for(int i=0;
i<(int)X.size();
++i)BU[i]=1;
faces=AY;
BR.assign(faces.size(),1);
BE=(int)faces.size();
vector<int>deg(N,0);
for(const Face&f:faces){
++deg[f.v[0]];
++deg[f.v[1]];
++deg[f.v[2]];
}
Y.assign(N,{
}
);
for(int i=0;
i<N;
++i)Y[i].reserve(deg[i]);
for(int i=0;
i<(int)faces.size();
++i){
Y[faces[i].v[0]].push_back(i);
Y[faces[i].v[1]].push_back(i);
Y[faces[i].v[2]].push_back(i);
}
return true;
}
static int count_output_vertices_estimate();
static bool IC(){
if(N<1000||elapsed_seconds()>16.4)return false;
const int safe=count_output_vertices_estimate();
if(safe<=8||safe>20)return false;
Vec3 mn=originalP[0],mx=originalP[0];
for(const Vec3&p:originalP){
mn.x=min(mn.x,p.x);
mn.y=min(mn.y,p.y);
mn.z=min(mn.z,p.z);
mx.x=max(mx.x,p.x);
mx.y=max(mx.y,p.y);
mx.z=max(mx.z,p.z);
}
double ex=mx.x-mn.x,ey=mx.y-mn.y,ez=mx.z-mn.z;
if(ex<=CL*1e-9||ey<=CL*1e-9||ez<=CL*1e-9)return false;
vector<Vec3>vs={
{
mn.x,mn.y,mn.z}
,{
mx.x,mn.y,mn.z}
,{
mx.x,mx.y,mn.z}
,{
mn.x,mx.y,mn.z}
,{
mn.x,mn.y,mx.z}
,{
mx.x,mn.y,mx.z}
,{
mx.x,mx.y,mx.z}
,{
mn.x,mx.y,mx.z}
}
;
vector<Face>fs(12);
int t[12][3]={
{
0,1,2}
,{
0,2,3}
,{
4,6,5}
,{
4,7,6}
,{
0,4,5}
,{
0,5,1}
,{
1,5,6}
,{
1,6,2}
,{
2,6,7}
,{
2,7,3}
,{
3,7,4}
,{
3,4,0}
}
;
for(int i=0;
i<12;
++i){
fs[i].v[0]=t[i][0];
fs[i].v[1]=t[i][1];
fs[i].v[2]=t[i][2];
}
AP st=AD();
bool keep=false;
if(AF(vs,fs)&&elapsed_seconds()<17.0)keep=visual_proxy_score(512)>=0.965;
if(!keep){
restore_state(st);
return false;
}
return true;
}
static bool EQ(const CE&fit,int lat,int lon,vector<Vec3>&verts,vector<Face>&AB){
if(!fit.ok||lat<4||lon<8)return false;
const int vertex_count=2+(lat-1)*lon;
if(vertex_count>N)return false;
verts.clear();
AB.clear();
verts.reserve(vertex_count);
AB.reserve(2*lon*(lat-1));
const double CC=BC(fit.BS);
auto make_point=[&](double x,double y,double z){
Vec3 p=fit.BS;
p=p+fit.BP[0]*(fit.AK[0]*x);
p=p+fit.BP[1]*(fit.AK[1]*y);
p=p+fit.BP[2]*(fit.AK[2]*z);
return p;
}
;
verts.push_back(make_point(0.0,0.0,1.0));
verts.push_back(make_point(0.0,0.0,-1.0));
auto ring_id=[&](int ring,int j){
j=(j%lon+lon)%lon;
return 2+(ring-1)*lon+j;
}
;
for(int ring=1;
ring<=lat-1;
++ring){
const double theta=acos(-1.0)*(double)ring/(double)lat;
const double st=sin(theta);
const double ct=cos(theta);
for(int j=0;
j<lon;
++j){
const double phi=2.0*acos(-1.0)*(double)j/(double)lon;
verts.push_back(make_point(st*cos(phi),st*sin(phi),ct));
}
}
auto add_face=[&](int a,int b,int c){
Face f;
f.v[0]=a;
f.v[1]=b;
f.v[2]=c;
BD(verts,f,fit.BS,CC);
AB.push_back(f);
}
;
for(int j=0;
j<lon;
++j){
add_face(0,ring_id(1,j),ring_id(1,j+1));
}
for(int ring=1;
ring<=lat-2;
++ring){
for(int j=0;
j<lon;
++j){
const int a=ring_id(ring,j);
const int b=ring_id(ring+1,j);
const int c=ring_id(ring+1,j+1);
const int d=ring_id(ring,j+1);
add_face(a,b,c);
add_face(a,c,d);
}
}
for(int j=0;
j<lon;
++j){
add_face(1,ring_id(lat-1,j+1),ring_id(lat-1,j));
}
return true;
}
static CE EJ(const Vec3 basis[3]){
CE fit;
for(int k=0;
k<3;
++k)fit.BP[k]=basis[k];
double lo[3]={
1e100,1e100,1e100}
;
double hi[3]={
-1e100,-1e100,-1e100}
;
for(const Vec3&p:originalP){
for(int k=0;
k<3;
++k){
const double t=GE(p,fit.BP[k]);
lo[k]=min(lo[k],t);
hi[k]=max(hi[k],t);
}
}
fit.BS=Vec3{
}
;
for(int k=0;
k<3;
++k){
const double mid=0.5*(lo[k]+hi[k]);
fit.BS=fit.BS+fit.BP[k]*mid;
fit.AK[k]=0.5*(hi[k]-lo[k]);
if(!(fit.AK[k]>1e-12))return fit;
}
const int AW=200000;
const int stride=max(1,N/AW);
int sampled=0;
double sum_abs=0.0;
double sum_sq=0.0;
double max_abs=0.0;
for(int i=0;
i<N;
i+=stride){
const Vec3 q=originalP[i]-fit.BS;
double r2=0.0;
for(int k=0;
k<3;
++k){
const double u=GE(q,fit.BP[k])/fit.AK[k];
r2+=u*u;
}
const double err=fabs(sqrt(max(0.0,r2))-1.0);
sum_abs+=err;
sum_sq+=err*err;
max_abs=max(max_abs,err);
++sampled;
}
if(sampled==0)return fit;
fit.GF=sum_abs/(double)sampled;
fit.EU=sqrt(sum_sq/(double)sampled);
fit.FA=max_abs;
const double max_allowed=(N<5000?0.010:0.014);
const double rms_allowed=(N<5000?0.0035:0.0045);
fit.ok=fit.FA<=max_allowed&&fit.EU<=rms_allowed;
return fit;
}
static void ER(double a[3][3],Vec3 out_axis[3]){
double v[3][3]={
{
1.0,0.0,0.0}
,{
0.0,1.0,0.0}
,{
0.0,0.0,1.0}
}
;
for(int it=0;
it<32;
++it){
int p=0,q=1;
double best=fabs(a[0][1]);
if(fabs(a[0][2])>best){
p=0;
q=2;
best=fabs(a[0][2]);
}
if(fabs(a[1][2])>best){
p=1;
q=2;
best=fabs(a[1][2]);
}
if(best<1e-18)break;
const double app=a[p][p];
const double aqq=a[q][q];
const double apq=a[p][q];
const double tau=(aqq-app)/(2.0*apq);
const double t=(tau>=0.0?1.0:-1.0)/(fabs(tau)+sqrt(1.0+tau*tau));
const double c=1.0/sqrt(1.0+t*t);
const double s=t*c;
for(int k=0;
k<3;
++k){
if(k==p||k==q)continue;
const double akp=a[k][p];
const double akq=a[k][q];
a[k][p]=a[p][k]=c*akp-s*akq;
a[k][q]=a[q][k]=s*akp+c*akq;
}
a[p][p]=c*c*app-2.0*s*c*apq+s*s*aqq;
a[q][q]=s*s*app+2.0*s*c*apq+c*c*aqq;
a[p][q]=a[q][p]=0.0;
for(int k=0;
k<3;
++k){
const double vkp=v[k][p];
const double vkq=v[k][q];
v[k][p]=c*vkp-s*vkq;
v[k][q]=s*vkp+c*vkq;
}
}
int ord[3]={
0,1,2}
;
sort(ord,ord+3,[&](int lhs,int rhs){
return a[lhs][lhs]>a[rhs][rhs];
}
);
for(int j=0;
j<3;
++j){
const int col=ord[j];
Vec3 e{
v[0][col],v[1][col],v[2][col]}
;
const double len=norm3(e);
out_axis[j]=len>0.0?JS(e,len):Vec3{
}
;
}
if(dot3(cross3(out_axis[0],out_axis[1]),out_axis[2])<0.0){
out_axis[2]=out_axis[2]*-1.0;
}
}
static CE GA(){
Vec3 identity[3]={
{
1.0,0.0,0.0}
,{
0.0,1.0,0.0}
,{
0.0,0.0,1.0}
}
;
CE best=EJ(identity);
Vec3 mean{
}
;
for(const Vec3&p:originalP)mean=mean+p;
mean=mean*(1.0/max(1,N));
double cov[3][3]={
}
;
for(const Vec3&p:originalP){
const Vec3 q=p-mean;
const double x[3]={
q.x,q.y,q.z}
;
for(int i=0;
i<3;
++i){
for(int j=0;
j<3;
++j)cov[i][j]+=x[i]*x[j];
}
}
for(int i=0;
i<3;
++i){
for(int j=0;
j<3;
++j)cov[i][j]/=(double)max(1,N);
}
Vec3 pca_axis[3];
ER(cov,pca_axis);
if(norm2(pca_axis[0])>0.0&&norm2(pca_axis[1])>0.0&&norm2(pca_axis[2])>0.0){
CE pca=EJ(pca_axis);
if(pca.ok&&(!best.ok||pca.EU<best.EU))best=pca;
}
return best;
}
static int count_output_vertices_estimate();
static CE GY(){
CE fit;
if(N<900||originalP.empty())return fit;
Vec3 mn=originalP[0];
Vec3 mx=originalP[0];
for(const Vec3&p:originalP){
mn.x=min(mn.x,p.x);
mn.y=min(mn.y,p.y);
mn.z=min(mn.z,p.z);
mx.x=max(mx.x,p.x);
mx.y=max(mx.y,p.y);
mx.z=max(mx.z,p.z);
}
fit.BS=(mn+mx)*0.5;
const double ex=mx.x-mn.x;
const double ey=mx.y-mn.y;
const double ez=mx.z-mn.z;
const double max_extent=max(ex,max(ey,ez));
const double min_extent=min(ex,min(ey,ez));
if(!(max_extent>1e-12)||min_extent<max_extent*0.975)return fit;
const int AW=240000;
const int stride=max(1,N/AW);
double sum_r=0.0;
int sampled=0;
for(int i=0;
i<N;
i+=stride){
sum_r+=norm3(originalP[i]-fit.BS);
++sampled;
}
if(sampled==0)return fit;
const double mean_r=sum_r/(double)sampled;
if(!(mean_r>1e-12))return fit;
double sum_sq=0.0;
double max_abs=0.0;
for(int i=0;
i<N;
i+=stride){
const double dev=fabs(norm3(originalP[i]-fit.BS)-mean_r);
sum_sq+=dev*dev;
max_abs=max(max_abs,dev);
}
const double rel_rms=sqrt(sum_sq/(double)sampled)/mean_r;
const double rel_max=max_abs/mean_r;
const double rms_limit=(N<5000?0.0060:0.0045);
const double max_limit=(N<5000?0.0260:0.0180);
if(rel_rms>rms_limit||rel_max>max_limit)return fit;
fit.ok=true;
fit.BP[0]={
1.0,0.0,0.0}
;
fit.BP[1]={
0.0,1.0,0.0}
;
fit.BP[2]={
0.0,0.0,1.0}
;
fit.AK[0]=fit.AK[1]=fit.AK[2]=mean_r;
fit.GF=rel_rms;
fit.EU=rel_rms;
fit.FA=rel_max;
return fit;
}
static bool GI(){
if(elapsed_seconds()>13.5)return false;
const CE fit=GY();
if(!fit.ok)return false;
int lat=20;
int lon=40;
int AG=256;
double proxy_threshold=0.930;
if(N<1500){
lat=20;
lon=40;
proxy_threshold=0.950;
}
else if(N<3000){
lat=20;
lon=40;
proxy_threshold=0.915;
}
else if(N<8000){
lat=20;
lon=40;
proxy_threshold=0.905;
}
else if(N<20000){
lat=24;
lon=48;
proxy_threshold=0.930;
}
else if(N<40000){
lat=28;
lon=56;
AG=192;
proxy_threshold=0.945;
}
else if(N<100000){
lat=28;
lon=56;
AG=192;
proxy_threshold=0.945;
}
else{
lat=28;
lon=56;
AG=128;
proxy_threshold=0.950;
}
vector<Vec3>X;
vector<Face>AY;
if(!EQ(fit,lat,lon,X,AY))return false;
const int safe_vertices=count_output_vertices_estimate();
if(safe_vertices<=0)return false;
if((int)X.size()>=safe_vertices)return false;
if((int)X.size()*100>safe_vertices*94)return false;
AP AM=AD();
bool keep=false;
if(AF(X,AY)&&elapsed_seconds()<16.5){
const double proxy=visual_proxy_score(AG);
keep=proxy>=proxy_threshold;
}
if(!keep)restore_state(AM);
return keep;
}
static inline bool DM(int fid,int a,int b,int c){
const Face&f=AR[fid];
return CW(f.v[0],f.v[1],f.v[2])==CW(a,b,c);
}
static bool GJ(int&R,int&V){
if(N<300||M!=2*(N-2))return false;
for(int v=8;
v<=1024;
++v){
if((N-2)%v)continue;
int r=(N-2)/v;
if(r<3)continue;
int step=max(1,v/64);
bool ok=true;
for(int j=0;
j<v&&ok;
j+=step){
int a=1+j,b=1+(j+1)%v;
if(!DM(j,0,b,a))ok=false;
int off=v+2*(r-1)*v+j;
int c=1+(r-1)*v+j,d=1+(r-1)*v+(j+1)%v;
if(ok&&!DM(off,N-1,c,d))ok=false;
}
int span=max(1,(r-1)*v/192);
for(int q=0;
q<(r-1)*v&&ok;
q+=span){
int rr=q/v,j=q-rr*v;
int a=1+rr*v+j,b=1+rr*v+(j+1)%v,c=1+(rr+1)*v+j,d=1+(rr+1)*v+(j+1)%v;
int f=v+2*(rr*v+j);
if(!DM(f,a,b,c)||!DM(f+1,b,d,c))ok=false;
}
if(ok){
R=r;
V=v;
return true;
}
}
return false;
}
static bool GR(int R,int V,int R2,int V2,vector<Vec3>&verts,vector<Face>&fs){
if(R2<3||V2<8||2+R2*V2>=N)return false;
verts.clear();
fs.clear();
verts.reserve(2+R2*V2);
fs.reserve(2*R2*V2);
verts.push_back(originalP[0]);
Vec3 mn=originalP[0],mx=originalP[0];
for(const Vec3&p:originalP){
mn.x=min(mn.x,p.x);
mn.y=min(mn.y,p.y);
mn.z=min(mn.z,p.z);
mx.x=max(mx.x,p.x);
mx.y=max(mx.y,p.y);
mx.z=max(mx.z,p.z);
}
Vec3 ctr=(mn+mx)*0.5;
double s=BC(ctr);
for(int i=0;
i<R2;
++i){
int oi=1+(int)((long long)i*(R-1)/max(1,R2-1));
for(int j=0;
j<V2;
++j){
int oj=(int)((long long)j*V/V2);
verts.push_back(originalP[1+(oi-1)*V+oj]);
}
}
int bot=(int)verts.size();
verts.push_back(originalP[N-1]);
auto rid=[&](int r,int j){
return 1+(r-1)*V2+(j%V2+V2)%V2;
}
;
auto add=[&](int a,int b,int c){
Face f;
f.v[0]=a;
f.v[1]=b;
f.v[2]=c;
BD(verts,f,ctr,s);
fs.push_back(f);
}
;
for(int j=0;
j<V2;
++j)add(0,rid(1,j+1),rid(1,j));
for(int r=1;
r<R2;
++r)
for(int j=0;
j<V2;
++j){
int a=rid(r,j),b=rid(r,j+1),c=rid(r+1,j),d=rid(r+1,j+1);
add(a,b,c);
add(b,d,c);
}
for(int j=0;
j<V2;
++j)add(bot,rid(R2,j),rid(R2,j+1));
return true;
}
static bool FL(){
if(elapsed_seconds()>13.2)return false;
int R=0,V=0;
if(!GJ(R,V))return false;
Vec3 mn=originalP[0],mx=originalP[0];
for(const Vec3&p:originalP){
mn.x=min(mn.x,p.x);
mn.y=min(mn.y,p.y);
mn.z=min(mn.z,p.z);
mx.x=max(mx.x,p.x);
mx.y=max(mx.y,p.y);
mx.z=max(mx.z,p.z);
}
double ex=mx.x-mn.x,ey=mx.y-mn.y,ez=mx.z-mn.z,hi=max(ex,max(ey,ez)),lo=min(ex,min(ey,ez));
if(!(hi>1e-12)||lo<0.70*hi)return false;
int safe=count_output_vertices_estimate();
if(safe<=0)return false;
struct T{
int r,c,res;
double th;
}
;
T ts[4]={
{
max(5,R/8),max(12,V/8),512,0.935}
,{
max(5,R/6),max(12,V/6),512,0.925}
,{
max(6,R/5),max(16,V/5),512,0.920}
,{
max(8,R/4),max(20,V/4),384,0.940}
}
;
AP st=AD();
vector<Vec3>vs;
vector<Face>fs;
for(const T&t:ts){
if(elapsed_seconds()>16.4)break;
int vc=2+t.r*t.c;
if(vc>=safe||vc>N)continue;
restore_state(st);
if(!GR(R,V,t.r,t.c,vs,fs))continue;
bool keep=false;
if(AF(vs,fs)&&elapsed_seconds()<16.8)keep=visual_proxy_score(t.res)>=t.th;
if(keep)return true;
}
restore_state(st);
return false;
}
struct AxisRevolveFit{
bool ok=false;
int BP=2;
double center_u=0.0;
double center_v=0.0;
double t0=0.0;
double t1=0.0;
double r0=0.0;
double r1=0.0;
double residual=1e100;
}
;
static inline void BM(const Vec3&p,int BP,double&t,double&u,double&v){
if(BP==0){
t=p.x;
u=p.y;
v=p.z;
}
else if(BP==1){
t=p.y;
u=p.x;
v=p.z;
}
else{
t=p.z;
u=p.x;
v=p.y;
}
}
static inline Vec3 BA(int BP,double t,double u,double v){
if(BP==0)return{
t,u,v}
;
if(BP==1)return{
u,t,v}
;
return{
u,v,t}
;
}
struct CT{
bool ok=false;
Vec3 w{
}
,u{
}
,v{
}
;
double cu=0.0,cv=0.0,t0=0.0,t1=0.0;
double rx0=0.0,rx1=0.0,ry0=0.0,ry1=0.0;
double residual=1e100;
}
;
static inline Vec3 IN(Vec3 a){
double l=norm3(a);
return l>1e-300?a*(1.0/l):Vec3{
}
;
}
static inline void CZ(const Vec3&p,const Vec3&w,const Vec3&u,const Vec3&v,double&t,double&x,double&y){
t=dot3(p,w);
x=dot3(p,u);
y=dot3(p,v);
}
static inline Vec3 DA(const Vec3&w,const Vec3&u,const Vec3&v,double t,double x,double y){
return w*t+u*x+v*y;
}
static bool IO(Vec3 BP[3]){
if(N<600||originalP.empty())return false;
Vec3 mean{
}
;
for(const Vec3&p:originalP)mean=mean+p;
mean=mean*(1.0/max(1,N));
double cov[3][3]={
}
;
const int stride=max(1,N/240000);
int sampled=0;
for(int i=0;
i<N;
i+=stride){
Vec3 q=originalP[i]-mean;
double x[3]={
q.x,q.y,q.z}
;
for(int a=0;
a<3;
++a)
for(int b=0;
b<3;
++b)cov[a][b]+=x[a]*x[b];
++sampled;
}
if(sampled<50)return false;
for(int a=0;
a<3;
++a)
for(int b=0;
b<3;
++b)cov[a][b]/=(double)sampled;
ER(cov,BP);
return norm2(BP[0])>0.5&&norm2(BP[1])>0.5&&norm2(BP[2])>0.5;
}
static void ID(Vec3 w,Vec3 hint,Vec3&u,Vec3&v){
w=IN(w);
u=hint-w*dot3(hint,w);
if(norm2(u)<1e-20){
Vec3 h=fabs(w.x)<0.75?Vec3{
1,0,0}
:(fabs(w.y)<0.75?Vec3{
0,1,0}
:Vec3{
0,0,1}
);
u=h-w*dot3(h,w);
}
u=IN(u);
v=IN(cross3(w,u));
}
static CT GB(Vec3 w,Vec3 hint){
CT fit;
fit.w=IN(w);
ID(fit.w,hint,fit.u,fit.v);
if(N<1000||norm2(fit.w)<0.5||norm2(fit.u)<0.5||norm2(fit.v)<0.5)return fit;
double min_t=1e100,max_t=-1e100,min_u=1e100,max_u=-1e100,min_v=1e100,max_v=-1e100;
for(const Vec3&p:originalP){
double t,x,y;
CZ(p,fit.w,fit.u,fit.v,t,x,y);
min_t=min(min_t,t);
max_t=max(max_t,t);
min_u=min(min_u,x);
max_u=max(max_u,x);
min_v=min(min_v,y);
max_v=max(max_v,y);
}
if(!(max_t>min_t))return fit;
fit.t0=min_t;
fit.t1=max_t;
fit.cu=0.5*(min_u+max_u);
fit.cv=0.5*(min_v+max_v);
const double len=max_t-min_t;
const int stride=max(1,N/240000);
double max_r=0.0;
for(int i=0;
i<N;
i+=stride){
double t,x,y;
CZ(originalP[i],fit.w,fit.u,fit.v,t,x,y);
x-=fit.cu;
y-=fit.cv;
max_r=max(max_r,sqrt(x*x+y*y));
}
if(!(max_r>1e-12)||len<max_r*0.35)return fit;
const int B=18;
double sx[B],sy[B];
int bc[B];
for(int i=0;
i<B;
++i)sx[i]=sy[i]=0.0,bc[i]=0;
const double eps=max_r*0.035;
int CI=0,EX=0;
for(int i=0;
i<N;
i+=stride){
double t,x,y;
CZ(originalP[i],fit.w,fit.u,fit.v,t,x,y);
x-=fit.cu;
y-=fit.cv;
double r=sqrt(x*x+y*y);
if(r<=eps){
double near_end=min(fabs(t-min_t),fabs(t-max_t));
if(near_end>len*0.035)return fit;
++EX;
continue;
}
int b=(int)((t-min_t)/len*B);
if(b<0)b=0;
else if(b>=B)b=B-1;
sx[b]=max(sx[b],fabs(x));
sy[b]=max(sy[b],fabs(y));
++bc[b];
++CI;
}
if(CI<200||EX>CI/8+8)return fit;
double S=0,St=0,Stt=0,Sx=0,Stx=0,Sy=0,Sty=0;
int bins=0;
for(int b=0;
b<B;
++b)
if(bc[b]>0&&sx[b]>eps&&sy[b]>eps){
double s=((double)b+0.5)/(double)B;
S+=1.0;
St+=s;
Stt+=s*s;
Sx+=sx[b];
Stx+=s*sx[b];
Sy+=sy[b];
Sty+=s*sy[b];
++bins;
}
if(bins<5)return fit;
double den=S*Stt-St*St;
if(fabs(den)<1e-18)return fit;
double bx=(S*Stx-St*Sx)/den,ax=(Sx-bx*St)/S;
double by=(S*Sty-St*Sy)/den,ay=(Sy-by*St)/S;
fit.rx0=ax;
fit.rx1=ax+bx;
fit.ry0=ay;
fit.ry1=ay+by;
double mr=max(max(fit.rx0,fit.rx1),max(fit.ry0,fit.ry1));
double nr=min(min(fit.rx0,fit.rx1),min(fit.ry0,fit.ry1));
if(!(mr>1e-12)||nr<mr*0.10)return fit;
double sum_sq=0.0,max_abs=0.0;
int checked=0;
for(int i=0;
i<N;
i+=stride){
double t,x,y;
CZ(originalP[i],fit.w,fit.u,fit.v,t,x,y);
x-=fit.cu;
y-=fit.cv;
double r=sqrt(x*x+y*y);
if(r<=eps)continue;
double s=(t-min_t)/len;
double rx=fit.rx0+(fit.rx1-fit.rx0)*s,ry=fit.ry0+(fit.ry1-fit.ry0)*s;
if(rx<=mr*0.08||ry<=mr*0.08)return fit;
double e=fabs(sqrt((x*x)/(rx*rx)+(y*y)/(ry*ry))-1.0);
sum_sq+=e*e;
max_abs=max(max_abs,e);
++checked;
}
if(checked<200)return fit;
fit.residual=sqrt(sum_sq/(double)checked);
if(fit.residual>0.022||max_abs>0.085)return fit;
fit.ok=true;
return fit;
}
static bool FQ(const CT&fit,int sides,vector<Vec3>&verts,vector<Face>&AB){
if(!fit.ok||sides<12)return false;
verts.clear();
AB.clear();
verts.reserve(2+2*sides);
AB.reserve(4*sides);
auto make=[&](double t,double rx,double ry,int j){
double th=2.0*acos(-1.0)*(double)j/(double)sides;
return DA(fit.w,fit.u,fit.v,t,fit.cu+rx*cos(th),fit.cv+ry*sin(th));
}
;
Vec3 BS=DA(fit.w,fit.u,fit.v,0.5*(fit.t0+fit.t1),fit.cu,fit.cv);
double sign=BC(BS);
verts.push_back(DA(fit.w,fit.u,fit.v,fit.t0,fit.cu,fit.cv));
verts.push_back(DA(fit.w,fit.u,fit.v,fit.t1,fit.cu,fit.cv));
const int lo=(int)verts.size();
for(int j=0;
j<sides;
++j)verts.push_back(make(fit.t0,fit.rx0,fit.ry0,j));
const int hi=(int)verts.size();
for(int j=0;
j<sides;
++j)verts.push_back(make(fit.t1,fit.rx1,fit.ry1,j));
if((int)verts.size()>N)return false;
auto add_face=[&](int a,int b,int c){
Face f;
f.v[0]=a;
f.v[1]=b;
f.v[2]=c;
BD(verts,f,BS,sign);
AB.push_back(f);
}
;
for(int j=0;
j<sides;
++j){
int k=(j+1)%sides,l0=lo+j,l1=lo+k,h0=hi+j,h1=hi+k;
add_face(l0,l1,h1);
add_face(l0,h1,h0);
add_face(0,l1,l0);
add_face(1,h0,h1);
}
return true;
}
static bool FR(){
if(elapsed_seconds()>13.7||N<1000||N>120000)return false;
const int safe_vertices=count_output_vertices_estimate();
if(safe_vertices<=0)return false;
Vec3 ax[3];
if(!IO(ax))return false;
CT best;
for(int k=0;
k<3;
++k){
CT fit=GB(ax[k],ax[(k+1)%3]);
if(fit.ok&&(!best.ok||fit.residual<best.residual))best=fit;
}
if(!best.ok)return false;
const int AQ[]={
24,32,40,48,56,64,80,96}
;
AP AM=AD();
vector<Vec3>X;
vector<Face>AY;
for(int sides:AQ){
if(elapsed_seconds()>16.2)break;
if(2+2*sides>=safe_vertices)continue;
if(!FQ(best,sides,X,AY))continue;
restore_state(AM);
bool keep=false;
if(AF(X,AY)&&elapsed_seconds()<16.7){
int AG=(N<30000?512:256);
double threshold=(N<30000?0.920:0.935);
keep=visual_proxy_score(AG)>=threshold;
}
if(keep)return true;
}
restore_state(AM);
return false;
}
static inline bool BV(const Face&f,int a,int b,int c){
array<int,3>x{
f.v[0],f.v[1],f.v[2]}
,y{
a,b,c}
;
sort(x.begin(),x.end());
sort(y.begin(),y.end());
return x==y;
}
static bool GC(int&U,int&V){
if(N<300||M!=2*N)return false;
for(int vv=6;
vv<=192;
++vv){
if(N%vv)continue;
int uu=N/vv;
if(uu<vv||uu<8)continue;
int step=max(1,N/384),checked=0;
bool ok=true;
for(int q=0;
q<N&&checked<384;
q+=step,++checked){
int i=q/vv,j=q-i*vv,ni=(i+1==uu?0:i+1),nj=(j+1==vv?0:j+1);
int a=i*vv+j,b=i*vv+nj,c=ni*vv+j,d=ni*vv+nj;
int f0=2*q,f1=f0+1;
if(f1>=M){
ok=false;
break;
}
bool m0=BV(AR[f0],a,b,c)||BV(AR[f0],b,c,d)||BV(AR[f0],a,c,d)||BV(AR[f0],a,b,d);
bool m1=BV(AR[f1],a,b,c)||BV(AR[f1],b,c,d)||BV(AR[f1],a,c,d)||BV(AR[f1],a,b,d);
if(!(m0&&m1)){
ok=false;
break;
}
}
if(ok){
U=uu;
V=vv;
return true;
}
}
return false;
}
static bool GK(int U,int V,int U2,int V2,bool flip,vector<Vec3>&verts,vector<Face>&AB){
if(U2<6||V2<6||U2*V2>=N)return false;
verts.clear();
AB.clear();
verts.reserve(U2*V2);
AB.reserve(U2*V2*2);
for(int i=0;
i<U2;
++i){
int oi=(int)((long long)i*U/U2);
for(int j=0;
j<V2;
++j){
int oj=(int)((long long)j*V/V2);
verts.push_back(originalP[oi*V+oj]);
}
}
auto id=[&](int i,int j){
if(i==U2)i=0;
if(j==V2)j=0;
return i*V2+j;
}
;
auto add=[&](int a,int b,int c){
Face f;
f.v[0]=a;
f.v[1]=b;
f.v[2]=c;
if(flip)swap(f.v[1],f.v[2]);
AB.push_back(f);
}
;
for(int i=0;
i<U2;
++i)
for(int j=0;
j<V2;
++j){
int a=id(i,j),b=id(i,j+1),c=id(i+1,j),d=id(i+1,j+1);
add(a,b,c);
add(b,d,c);
}
return true;
}
static bool FG(){
if(elapsed_seconds()>13.5)return false;
int U=0,V=0;
if(!GC(U,V))return false;
int safe_vertices=count_output_vertices_estimate();
if(safe_vertices<=0)return false;
struct Trial{
int u,v;
}
;
Trial AQ[8]={
{
max(8,U/6),max(6,V/3)}
,{
max(8,U/5),max(6,V/2)}
,{
max(8,U/4),max(6,(V*7+19)/20)}
,{
max(10,U/3),max(8,V/3)}
,{
max(12,U/3),max(10,V/2)}
,{
max(14,U/2),max(10,V/2)}
,{
max(14,U/2),V}
,{
U,V}
}
;
AP AM=AD();
vector<Vec3>verts;
vector<Face>fs;
for(const Trial&tr:AQ){
if(elapsed_seconds()>16.5)break;
int vc=tr.u*tr.v;
if(vc<=0||vc>=safe_vertices||vc>N)continue;
for(int flip=0;
flip<2;
++flip){
restore_state(AM);
if(!GK(U,V,tr.u,tr.v,flip!=0,verts,fs))continue;
bool keep=false;
if(AF(verts,fs)&&elapsed_seconds()<16.8){
int res=(N<30000?512:256);
double th=(N<1500?0.925:(N<30000?0.925:0.935));
keep=visual_proxy_score(res)>=th;
}
if(keep)return true;
}
}
restore_state(AM);
return false;
}
struct EW{
bool ok=false;
int BP=2;
double center_t=0.0;
double center_u=0.0;
double center_v=0.0;
double major=0.0;
double minor=0.0;
double rms_rel=1e100;
double max_rel=1e100;
}
;
static EW IE(int BP){
EW fit;
fit.BP=BP;
if(N<600||originalP.empty())return fit;
double min_t=1e100,max_t=-1e100;
double min_u=1e100,max_u=-1e100;
double min_v=1e100,max_v=-1e100;
for(const Vec3&p:originalP){
double t,u,v;
BM(p,BP,t,u,v);
min_t=min(min_t,t);
max_t=max(max_t,t);
min_u=min(min_u,u);
max_u=max(max_u,u);
min_v=min(min_v,v);
max_v=max(max_v,v);
}
if(!(max_t>min_t)||!(max_u>min_u)||!(max_v>min_v))return fit;
fit.center_t=0.5*(min_t+max_t);
fit.center_u=0.5*(min_u+max_u);
fit.center_v=0.5*(min_v+max_v);
double min_r=1e100;
double max_r=0.0;
for(const Vec3&p:originalP){
double t,u,v;
BM(p,BP,t,u,v);
const double du=u-fit.center_u;
const double dv=v-fit.center_v;
const double r=sqrt(du*du+dv*dv);
min_r=min(min_r,r);
max_r=max(max_r,r);
}
if(!(max_r>min_r)||!(min_r>1e-10))return fit;
const double major=0.5*(max_r+min_r);
const double minor_r=0.5*(max_r-min_r);
const double minor_t=0.5*(max_t-min_t);
const double minor=0.5*(minor_r+minor_t);
if(!(major>1e-10)||!(minor>1e-10))return fit;
if(major<minor*1.35)return fit;
if(fabs(minor_r-minor_t)>minor*0.22)return fit;
const int AW=240000;
const int stride=max(1,N/AW);
double sum_sq=0.0;
double max_abs=0.0;
int sampled=0;
for(int i=0;
i<N;
i+=stride){
double t,u,v;
BM(originalP[i],BP,t,u,v);
const double du=u-fit.center_u;
const double dv=v-fit.center_v;
const double rho=sqrt(du*du+dv*dv);
const double tube=sqrt((rho-major)*(rho-major)+(t-fit.center_t)*(t-fit.center_t));
const double err=fabs(tube-minor);
sum_sq+=err*err;
max_abs=max(max_abs,err);
++sampled;
}
if(sampled<200)return fit;
fit.major=major;
fit.minor=minor;
fit.rms_rel=sqrt(sum_sq/(double)sampled)/minor;
fit.max_rel=max_abs/minor;
const double rms_limit=(N<3000?0.018:0.012);
const double max_limit=(N<3000?0.080:0.055);
fit.ok=fit.rms_rel<=rms_limit&&fit.max_rel<=max_limit;
return fit;
}
static bool HG(const EW&fit,int CK,int BW,vector<Vec3>&verts,vector<Face>&AB){
if(!fit.ok||CK<12||BW<6)return false;
const int vertex_count=CK*BW;
if(vertex_count>N)return false;
verts.clear();
AB.clear();
verts.reserve(vertex_count);
AB.reserve(vertex_count*2);
const double pi=acos(-1.0);
auto torus_normal_at=[&](const Vec3&p){
double t,u,v;
BM(p,fit.BP,t,u,v);
const double du=u-fit.center_u;
const double dv=v-fit.center_v;
const double rho=sqrt(du*du+dv*dv);
if(!(rho>1e-12)||!(fit.minor>1e-12))return Vec3{
0.0,0.0,0.0}
;
const double cu=du/rho;
const double sv=dv/rho;
double cp=(rho-fit.major)/fit.minor;
double sp=(t-fit.center_t)/fit.minor;
const double len=sqrt(cp*cp+sp*sp);
if(len>1e-12){
cp/=len;
sp/=len;
}
return BA(fit.BP,sp,cp*cu,cp*sv);
}
;
double orient_sum=0.0;
int orient_count=0;
const int AW=100000;
const int stride=max(1,M/AW);
for(int fid=0;
fid<M;
fid+=stride){
const Face&f=AR[fid];
const Vec3&a=originalP[f.v[0]];
const Vec3&b=originalP[f.v[1]];
const Vec3&c=originalP[f.v[2]];
Vec3 cr=cross3(b-a,c-a);
const double clen=norm3(cr);
if(!(clen>0.0))continue;
const Vec3 ctr=(a+b+c)*(1.0/3.0);
Vec3 pred=torus_normal_at(ctr);
const double plen=norm3(pred);
if(!(plen>0.0))continue;
orient_sum+=dot3(cr*(1.0/clen),pred*(1.0/plen));
++orient_count;
}
const bool flip_orientation=orient_count>0&&orient_sum<0.0;
for(int i=0;
i<CK;
++i){
const double theta=2.0*pi*(double)i/(double)CK;
const double ct=cos(theta);
const double st=sin(theta);
for(int j=0;
j<BW;
++j){
const double phi=2.0*pi*(double)j/(double)BW;
const double cp=cos(phi);
const double sp=sin(phi);
const double radial=fit.major+fit.minor*cp;
verts.push_back(BA(fit.BP,fit.center_t+fit.minor*sp,fit.center_u+radial*ct,fit.center_v+radial*st));
}
}
auto id=[&](int i,int j){
i=(i%CK+CK)%CK;
j=(j%BW+BW)%BW;
return i*BW+j;
}
;
auto add_face=[&](int a,int b,int c){
Face f;
f.v[0]=a;
f.v[1]=b;
f.v[2]=c;
if(flip_orientation)swap(f.v[1],f.v[2]);
AB.push_back(f);
}
;
for(int i=0;
i<CK;
++i){
for(int j=0;
j<BW;
++j){
const int a=id(i,j);
const int b=id(i+1,j);
const int c=id(i+1,j+1);
const int d=id(i,j+1);
add_face(a,b,c);
add_face(a,c,d);
}
}
return true;
}
static bool HH(){
if(elapsed_seconds()>13.8||N<600||N>30000)return false;
const int safe_vertices=count_output_vertices_estimate();
if(safe_vertices<=0)return false;
EW best;
for(int BP=0;
BP<3;
++BP){
const EW fit=IE(BP);
if(fit.ok&&(!best.ok||fit.rms_rel<best.rms_rel))best=fit;
}
if(!best.ok)return false;
struct TorusTrial{
int CK;
int BW;
int AG;
double threshold;
double keep_ratio;
}
;
vector<TorusTrial>AQ;
if(N<900){
AQ.push_back({
36,10,256,0.965,0.98}
);
}
else if(N<1500){
AQ.push_back({
48,14,512,0.930,0.98}
);
AQ.push_back({
56,16,256,0.945,0.98}
);
}
else if(N<5000){
AQ.push_back({
56,16,512,0.925,0.98}
);
AQ.push_back({
64,16,512,0.940,0.98}
);
AQ.push_back({
72,18,256,0.965,0.98}
);
AQ.push_back({
80,20,256,0.970,0.95}
);
}
else if(N<15000){
AQ.push_back({
80,20,256,0.950,0.90}
);
AQ.push_back({
96,24,256,0.960,0.98}
);
}
else{
AQ.push_back({
104,24,256,0.955,0.90}
);
AQ.push_back({
112,28,256,0.965,0.98}
);
}
AP AM=AD();
vector<Vec3>X;
vector<Face>AY;
for(const TorusTrial&AZ:AQ){
if(elapsed_seconds()>16.0)break;
const int AV=AZ.CK*AZ.BW;
if(AV<=0||AV>=safe_vertices)continue;
if((double)AV>(double)safe_vertices*AZ.keep_ratio)continue;
if(!HG(best,AZ.CK,AZ.BW,X,AY))continue;
restore_state(AM);
bool keep=false;
if(AF(X,AY)&&elapsed_seconds()<16.5){
const double proxy=visual_proxy_score(AZ.AG);
keep=proxy>=AZ.threshold;
}
if(keep)return true;
}
restore_state(AM);
return false;
}
static AxisRevolveFit GZ(int BP){
AxisRevolveFit fit;
fit.BP=BP;
if(N<1000||originalP.empty())return fit;
double min_t=1e100,max_t=-1e100;
double min_u=1e100,max_u=-1e100;
double min_v=1e100,max_v=-1e100;
for(const Vec3&p:originalP){
double t,u,v;
BM(p,BP,t,u,v);
min_t=min(min_t,t);
max_t=max(max_t,t);
min_u=min(min_u,u);
max_u=max(max_u,u);
min_v=min(min_v,v);
max_v=max(max_v,v);
}
if(!(max_t>min_t))return fit;
fit.center_u=0.5*(min_u+max_u);
fit.center_v=0.5*(min_v+max_v);
fit.t0=min_t;
fit.t1=max_t;
double max_r=0.0;
for(const Vec3&p:originalP){
double t,u,v;
BM(p,BP,t,u,v);
const double du=u-fit.center_u;
const double dv=v-fit.center_v;
max_r=max(max_r,sqrt(du*du+dv*dv));
}
if(!(max_r>1e-10))return fit;
const int AW=240000;
const int stride=max(1,N/AW);
const double radial_eps=max_r*0.055;
double s_t=0.0,s_r=0.0,s_tt=0.0,s_tr=0.0;
int CI=0;
int EX=0;
for(int i=0;
i<N;
i+=stride){
double t,u,v;
BM(originalP[i],BP,t,u,v);
const double du=u-fit.center_u;
const double dv=v-fit.center_v;
const double r=sqrt(du*du+dv*dv);
if(r<=radial_eps){
const double near_end=min(fabs(t-min_t),fabs(t-max_t));
if(near_end>(max_t-min_t)*0.030)return fit;
++EX;
continue;
}
s_t+=t;
s_r+=r;
s_tt+=t*t;
s_tr+=t*r;
++CI;
}
if(CI<200)return fit;
const double den=(double)CI*s_tt-s_t*s_t;
if(fabs(den)<1e-18)return fit;
const double a=((double)CI*s_tr-s_t*s_r)/den;
const double b=(s_r-a*s_t)/(double)CI;
double r0=a*min_t+b;
double r1=a*max_t+b;
if(r0<-0.035*max_r||r1<-0.035*max_r)return fit;
if(fabs(r0)<0.035*max_r)r0=0.0;
if(fabs(r1)<0.035*max_r)r1=0.0;
if(max(r0,r1)<0.25*max_r)return fit;
double sum_sq=0.0;
double max_abs=0.0;
int checked=0;
for(int i=0;
i<N;
i+=stride){
double t,u,v;
BM(originalP[i],BP,t,u,v);
const double du=u-fit.center_u;
const double dv=v-fit.center_v;
const double r=sqrt(du*du+dv*dv);
if(r<=radial_eps)continue;
const double pred=max(0.0,a*t+b);
const double err=fabs(r-pred);
sum_sq+=err*err;
max_abs=max(max_abs,err);
++checked;
}
if(checked<200)return fit;
const double rms=sqrt(sum_sq/(double)checked);
if(rms>max_r*0.0060||max_abs>max_r*0.030)return fit;
if(EX>CI/3)return fit;
fit.r0=r0;
fit.r1=r1;
fit.residual=rms/max_r;
fit.ok=true;
return fit;
}
static bool GL(const AxisRevolveFit&fit,int sides,vector<Vec3>&verts,vector<Face>&AB){
if(!fit.ok||sides<12)return false;
const double eps=max(fit.r0,fit.r1)*1e-6;
const bool cone0=fit.r0<=eps;
const bool cone1=fit.r1<=eps;
if(cone0&&cone1)return false;
const Vec3 BS=BA(fit.BP,0.5*(fit.t0+fit.t1),fit.center_u,fit.center_v);
const double sign=BC(BS);
const double pi=acos(-1.0);
verts.clear();
AB.clear();
auto make=[&](double t,double r,int j){
const double th=2.0*pi*(double)j/(double)sides;
return BA(fit.BP,t,fit.center_u+r*cos(th),fit.center_v+r*sin(th));
}
;
auto add_face=[&](int a,int b,int c){
Face f;
f.v[0]=a;
f.v[1]=b;
f.v[2]=c;
BD(verts,f,BS,sign);
AB.push_back(f);
}
;
if(!cone0&&!cone1){
const int bottom_center=0;
const int top_center=1;
verts.push_back(BA(fit.BP,fit.t0,fit.center_u,fit.center_v));
verts.push_back(BA(fit.BP,fit.t1,fit.center_u,fit.center_v));
const int bottom_ring=(int)verts.size();
for(int j=0;
j<sides;
++j)verts.push_back(make(fit.t0,fit.r0,j));
const int top_ring=(int)verts.size();
for(int j=0;
j<sides;
++j)verts.push_back(make(fit.t1,fit.r1,j));
if((int)verts.size()>N)return false;
for(int j=0;
j<sides;
++j){
const int j2=(j+1)%sides;
const int b0=bottom_ring+j;
const int b1=bottom_ring+j2;
const int t0=top_ring+j;
const int t1=top_ring+j2;
add_face(b0,b1,t0);
add_face(b1,t1,t0);
add_face(bottom_center,b0,b1);
add_face(top_center,t1,t0);
}
}
else{
const bool apex_at_bottom=cone0;
const double apex_t=apex_at_bottom?fit.t0:fit.t1;
const double base_t=apex_at_bottom?fit.t1:fit.t0;
const double base_r=apex_at_bottom?fit.r1:fit.r0;
const int apex=0;
const int base_center=1;
verts.push_back(BA(fit.BP,apex_t,fit.center_u,fit.center_v));
verts.push_back(BA(fit.BP,base_t,fit.center_u,fit.center_v));
const int ring=(int)verts.size();
for(int j=0;
j<sides;
++j)verts.push_back(make(base_t,base_r,j));
if((int)verts.size()>N)return false;
for(int j=0;
j<sides;
++j){
const int j2=(j+1)%sides;
add_face(apex,ring+j,ring+j2);
add_face(base_center,ring+j2,ring+j);
}
}
return true;
}
struct DU{
bool ok=false;
int BP=2;
double center_t=0.0;
double center_u=0.0;
double center_v=0.0;
double AK=0.0;
double BH=0.0;
double rms_rel=1e100;
double max_rel=1e100;
}
;
static DU HI(int BP){
DU fit;
fit.BP=BP;
if(N<1500||originalP.empty())return fit;
double min_t=1e100,max_t=-1e100;
double min_u=1e100,max_u=-1e100;
double min_v=1e100,max_v=-1e100;
for(const Vec3&p:originalP){
double t,u,v;
BM(p,BP,t,u,v);
min_t=min(min_t,t);
max_t=max(max_t,t);
min_u=min(min_u,u);
max_u=max(max_u,u);
min_v=min(min_v,v);
max_v=max(max_v,v);
}
const double eu=max_u-min_u;
const double ev=max_v-min_v;
const double et=max_t-min_t;
if(!(eu>1e-12)||!(ev>1e-12)||!(et>1e-12))return fit;
if(min(eu,ev)<max(eu,ev)*0.985)return fit;
fit.center_t=0.5*(min_t+max_t);
fit.center_u=0.5*(min_u+max_u);
fit.center_v=0.5*(min_v+max_v);
fit.AK=0.25*(eu+ev);
fit.BH=0.5*et-fit.AK;
if(!(fit.AK>1e-12)||fit.BH<=fit.AK*0.18)return fit;
const int AW=240000;
const int stride=max(1,N/AW);
const double rel_limit=(N<20000?0.018:0.014);
const double max_limit=(N<20000?0.060:0.045);
int sampled=0;
int side_samples=0;
int cap_samples=0;
int failed=0;
double sum_sq=0.0;
double max_abs_rel=0.0;
for(int i=0;
i<N;
i+=stride){
double t,u,v;
BM(originalP[i],BP,t,u,v);
const double dt=t-fit.center_t;
const double du=u-fit.center_u;
const double dv=v-fit.center_v;
const double rho=sqrt(du*du+dv*dv);
const double at=fabs(dt);
double rel=1e100;
if(at<=fit.BH){
rel=fabs(rho-fit.AK)/fit.AK;
if(at<fit.BH*0.70)++side_samples;
}
else{
const double q=sqrt(rho*rho+(at-fit.BH)*(at-fit.BH));
rel=fabs(q-fit.AK)/fit.AK;
if(at>fit.BH+fit.AK*0.20)++cap_samples;
}
if(at>fit.BH+fit.AK*1.015)rel=1e100;
if(!(rel<=max_limit)){
++failed;
if(failed>max(2,sampled/60+2))return fit;
}
if(rel<1e90){
sum_sq+=rel*rel;
max_abs_rel=max(max_abs_rel,rel);
}
++sampled;
}
if(sampled<200)return fit;
if(side_samples<max(20,sampled/25))return fit;
if(cap_samples<max(20,sampled/25))return fit;
fit.rms_rel=sqrt(sum_sq/(double)sampled);
fit.max_rel=max_abs_rel;
fit.ok=fit.rms_rel<=rel_limit&&fit.max_rel<=max_limit;
return fit;
}
static int FH(int sides,int CN,int CJ){
if(sides<8||CN<3||CJ<1)return 0;
const int DV=2*CN+CJ-1;
return 2+DV*sides;
}
static bool GS(const DU&fit,int sides,int CN,int CJ,vector<Vec3>&verts,vector<Face>&AB){
if(!fit.ok)return false;
const int DV=2*CN+CJ-1;
const int vertex_count=FH(sides,CN,CJ);
if(vertex_count<=0||vertex_count>N)return false;
verts.clear();
AB.clear();
verts.reserve(vertex_count);
AB.reserve(2*DV*sides);
const Vec3 BS=BA(fit.BP,fit.center_t,fit.center_u,fit.center_v);
const double CC=BC(BS);
const double pi=acos(-1.0);
auto make=[&](double rr,double w,int j){
const double th=2.0*pi*(double)j/(double)sides;
return BA(fit.BP,fit.center_t+w,fit.center_u+rr*cos(th),fit.center_v+rr*sin(th));
}
;
verts.push_back(BA(fit.BP,fit.center_t+fit.BH+fit.AK,fit.center_u,fit.center_v));
vector<int>ring_start;
ring_start.reserve(DV);
auto add_ring=[&](double rr,double w){
ring_start.push_back((int)verts.size());
for(int j=0;
j<sides;
++j)verts.push_back(make(rr,w,j));
}
;
for(int i=1;
i<=CN;
++i){
const double phi=0.5*pi*(double)i/(double)CN;
add_ring(fit.AK*sin(phi),fit.BH+fit.AK*cos(phi));
}
for(int i=1;
i<=CJ;
++i){
const double t=(double)i/(double)CJ;
add_ring(fit.AK,fit.BH*(1.0-2.0*t));
}
for(int i=1;
i<CN;
++i){
const double phi=0.5*pi+0.5*pi*(double)i/(double)CN;
add_ring(fit.AK*sin(phi),-fit.BH+fit.AK*cos(phi));
}
const int bottom=(int)verts.size();
verts.push_back(BA(fit.BP,fit.center_t-fit.BH-fit.AK,fit.center_u,fit.center_v));
if((int)verts.size()!=vertex_count||(int)ring_start.size()!=DV)return false;
auto ring=[&](int r,int j){
return ring_start[r]+((j%sides+sides)%sides);
}
;
auto add_face=[&](int a,int b,int c){
Face f;
f.v[0]=a;
f.v[1]=b;
f.v[2]=c;
BD(verts,f,BS,CC);
AB.push_back(f);
}
;
for(int j=0;
j<sides;
++j)add_face(0,ring(0,j+1),ring(0,j));
for(int r=0;
r+1<DV;
++r){
for(int j=0;
j<sides;
++j){
const int a=ring(r,j);
const int b=ring(r,j+1);
const int c=ring(r+1,j);
const int d=ring(r+1,j+1);
add_face(a,b,c);
add_face(b,d,c);
}
}
for(int j=0;
j<sides;
++j)add_face(bottom,ring(DV-1,j),ring(DV-1,j+1));
return true;
}
static bool GT(){
if(elapsed_seconds()>13.8||N<1500)return false;
const int safe_vertices=count_output_vertices_estimate();
if(safe_vertices<=0)return false;
DU best;
const int axes[3]={
2,0,1}
;
for(int i=0;
i<3;
++i){
const DU fit=HI(axes[i]);
if(fit.ok&&(!best.ok||fit.rms_rel<best.rms_rel))best=fit;
}
if(!best.ok)return false;
struct CapsuleTrial{
int sides;
int CN;
int CJ;
int AG;
double threshold;
double keep_ratio;
}
;
vector<CapsuleTrial>AQ;
if(N<20000){
AQ.push_back({
24,8,1,512,0.925,0.92}
);
AQ.push_back({
32,10,1,512,0.930,0.97}
);
}
else if(N<100000){
AQ.push_back({
28,8,1,256,0.925,0.92}
);
AQ.push_back({
36,10,1,256,0.932,0.97}
);
}
else{
AQ.push_back({
32,8,1,192,0.928,0.92}
);
AQ.push_back({
40,10,1,192,0.935,0.97}
);
}
AP AM=AD();
vector<Vec3>X;
vector<Face>AY;
for(const CapsuleTrial&AZ:AQ){
if(elapsed_seconds()>16.0)break;
const int AV=FH(AZ.sides,AZ.CN,AZ.CJ);
if(AV<=0||AV>=safe_vertices)continue;
if((double)AV>(double)safe_vertices*AZ.keep_ratio)continue;
if(!GS(best,AZ.sides,AZ.CN,AZ.CJ,X,AY))continue;
restore_state(AM);
bool keep=false;
if(AF(X,AY)&&elapsed_seconds()<16.5){
const double proxy=visual_proxy_score(AZ.AG);
keep=proxy>=AZ.threshold;
}
if(keep)return true;
}
restore_state(AM);
return false;
}
static bool EO(){
if(elapsed_seconds()>13.8)return false;
const int safe_vertices=count_output_vertices_estimate();
if(safe_vertices<=0)return false;
for(int BP=0;
BP<3&&elapsed_seconds()<13.8;
++BP){
const AxisRevolveFit fit=GZ(BP);
if(!fit.ok)continue;
const bool cone_like=min(fit.r0,fit.r1)<=max(fit.r0,fit.r1)*0.08;
const int sides=cone_like?32:32;
vector<Vec3>X;
vector<Face>AY;
if(!GL(fit,sides,X,AY))continue;
if((int)X.size()>=safe_vertices)continue;
if((int)X.size()*100>safe_vertices*90)continue;
AP AM=AD();
bool keep=false;
if(AF(X,AY)&&elapsed_seconds()<16.5){
const int AG=(N<20000?512:256);
const double proxy=visual_proxy_score(AG);
keep=proxy>=0.945;
}
if(keep)return true;
restore_state(AM);
}
return false;
}
static bool FS(){
if(N<900||elapsed_seconds()>13.0)return false;
const CE fit=GA();
if(!fit.ok)return false;
int lat=20;
int lon=40;
int AG=512;
double proxy_threshold=0.930;
if(N<1500){
lat=20;
lon=40;
AG=256;
proxy_threshold=0.950;
}
else if(N<3000){
lat=20;
lon=40;
AG=256;
proxy_threshold=0.905;
}
else if(N<5000){
lat=20;
lon=40;
AG=256;
proxy_threshold=0.904;
}
else if(N<30000){
lat=20;
lon=40;
proxy_threshold=0.930;
}
else if(N<100000){
lat=22;
lon=44;
AG=256;
proxy_threshold=0.925;
}
else{
lat=22;
lon=44;
AG=256;
proxy_threshold=0.925;
}
vector<Vec3>X;
vector<Face>AY;
if(!EQ(fit,lat,lon,X,AY))return false;
const int safe_vertices=count_output_vertices_estimate();
if(safe_vertices<=0)return false;
if((int)X.size()>=safe_vertices)return false;
if((int)X.size()*100>safe_vertices*92)return false;
AP AM=AD();
bool keep=false;
if(AF(X,AY)&&elapsed_seconds()<16.5){
const double proxy=visual_proxy_score(AG);
keep=proxy>=proxy_threshold;
}
if(!keep)restore_state(AM);
return keep;
}
static int count_output_vertices_estimate(){
vector<unsigned char>used(N,0);
int cnt=0;
for(int fid=0;
fid<(int)faces.size();
++fid){
if(!BR[fid])continue;
const Face&f=faces[fid];
for(int k=0;
k<3;
++k){
int v=f.v[k];
if(v>=0&&v<N&&BU[v]&&!used[v]){
used[v]=1;
++cnt;
}
}
}
return cnt;
}
struct StrictSphereFit{
bool ok=false;
Vec3 BS{
}
;
double AK=1.0;
}
;
static StrictSphereFit HA(){
StrictSphereFit s;
if(N<900||originalP.empty())return s;
Vec3 mn=originalP[0],mx=originalP[0];
for(const Vec3&p:originalP){
mn.x=min(mn.x,p.x);
mn.y=min(mn.y,p.y);
mn.z=min(mn.z,p.z);
mx.x=max(mx.x,p.x);
mx.y=max(mx.y,p.y);
mx.z=max(mx.z,p.z);
}
const double ex=mx.x-mn.x;
const double ey=mx.y-mn.y;
const double ez=mx.z-mn.z;
const double max_extent=max(ex,max(ey,ez));
const double min_extent=min(ex,min(ey,ez));
if(!(max_extent>1e-12)||min_extent<max_extent*0.982)return s;
s.BS=(mn+mx)*0.5;
const int AW=200000;
const int stride=max(1,N/AW);
double sum_r=0.0;
double sum_r2=0.0;
int sampled=0;
for(int i=0;
i<N;
i+=stride){
const double r=norm3(originalP[i]-s.BS);
sum_r+=r;
sum_r2+=r*r;
++sampled;
}
if(sampled==0)return s;
const double inv=1.0/(double)sampled;
const double mean_r=sum_r*inv;
if(!(mean_r>1e-12))return s;
const double std_r=sqrt(max(0.0,sum_r2*inv-mean_r*mean_r));
double max_abs=0.0;
for(int i=0;
i<N;
i+=stride){
max_abs=max(max_abs,fabs(norm3(originalP[i]-s.BS)-mean_r));
}
const double rel_std=std_r/mean_r;
const double rel_max=max_abs/mean_r;
const double std_limit=N<5000?1.2e-4:1.8e-4;
const double max_limit=N<5000?8.5e-4:1.2e-3;
if(rel_std>std_limit||rel_max>max_limit)return s;
s.AK=mean_r;
s.ok=true;
return s;
}
static bool HQ(const StrictSphereFit&sphere,int lat,int lon,vector<Vec3>&verts,vector<Face>&AB){
const int BN=2+(lat-1)*lon;
if(!sphere.ok||lat<4||lon<8||BN>N)return false;
verts.clear();
AB.clear();
verts.reserve(BN);
AB.reserve(2*lat*lon);
const double pi=acos(-1.0);
const Vec3 c=sphere.BS;
const double r=sphere.AK;
const double CC=BC(c);
verts.push_back({
c.x,c.y,c.z+r}
);
for(int i=1;
i<lat;
++i){
const double phi=pi*(double)i/(double)lat;
const double z=cos(phi);
const double rr=sin(phi);
for(int j=0;
j<lon;
++j){
const double th=2.0*pi*(double)j/(double)lon;
verts.push_back({
c.x+r*rr*cos(th),c.y+r*rr*sin(th),c.z+r*z}
);
}
}
verts.push_back({
c.x,c.y,c.z-r}
);
auto vid=[&](int ring,int j){
return 1+(ring-1)*lon+((j%lon+lon)%lon);
}
;
auto add_face=[&](int a,int b,int cc){
Face f;
f.v[0]=a;
f.v[1]=b;
f.v[2]=cc;
BD(verts,f,c,CC);
AB.push_back(f);
}
;
const int bottom=BN-1;
for(int j=0;
j<lon;
++j)add_face(0,vid(1,j+1),vid(1,j));
for(int i=1;
i<lat-1;
++i){
for(int j=0;
j<lon;
++j){
const int a=vid(i,j);
const int b=vid(i,j+1);
const int cc=vid(i+1,j);
const int d=vid(i+1,j+1);
add_face(a,b,cc);
add_face(b,d,cc);
}
}
for(int j=0;
j<lon;
++j)add_face(bottom,vid(lat-1,j),vid(lat-1,j+1));
return true;
}
static bool GM(){
if(N<900||N>=50000||elapsed_seconds()>13.0)return false;
const StrictSphereFit sphere=HA();
if(!sphere.ok)return false;
const int safe_vertices=count_output_vertices_estimate();
if(safe_vertices<=0)return false;
struct SphereTrial{
int lat;
int lon;
int AG;
double threshold;
double keep_ratio;
}
;
vector<SphereTrial>AQ;
if(N<1500){
AQ.push_back({
12,24,1024,0.922,0.94}
);
AQ.push_back({
14,28,768,0.928,0.96}
);
AQ.push_back({
16,32,512,0.934,0.98}
);
}
else if(N<5000){
AQ.push_back({
16,32,512,0.925,0.92}
);
AQ.push_back({
18,36,512,0.935,0.96}
);
}
else if(N<20000){
AQ.push_back({
18,36,512,0.900,0.92}
);
AQ.push_back({
20,40,512,0.935,0.98}
);
}
else if(N<50000){
AQ.push_back({
22,44,256,0.920,0.92}
);
AQ.push_back({
24,48,256,0.935,0.96}
);
AQ.push_back({
28,56,256,0.945,0.98}
);
}
AP AM=AD();
vector<Vec3>X;
vector<Face>AY;
for(const SphereTrial&AZ:AQ){
if(elapsed_seconds()>16.0)break;
const int AV=2+(AZ.lat-1)*AZ.lon;
if(AV<=0||AV>=safe_vertices)continue;
if((double)AV>(double)safe_vertices*AZ.keep_ratio)continue;
if(!HQ(sphere,AZ.lat,AZ.lon,X,AY))continue;
restore_state(AM);
bool keep=false;
if(AF(X,AY)&&elapsed_seconds()<16.5){
keep=visual_proxy_score(AZ.AG)>=AZ.threshold;
}
if(keep)return true;
}
restore_state(AM);
return false;
}
namespace DC{
struct Vec3{
double x=0.0,y=0.0,z=0.0;
}
;
static inline Vec3 operator+(const Vec3&a,const Vec3&b){
return{
a.x+b.x,a.y+b.y,a.z+b.z}
;
}
static inline Vec3 operator-(const Vec3&a,const Vec3&b){
return{
a.x-b.x,a.y-b.y,a.z-b.z}
;
}
static inline Vec3 operator*(const Vec3&a,double k){
return{
a.x*k,a.y*k,a.z*k}
;
}
static inline double dot3(const Vec3&a,const Vec3&b){
return a.x*b.x+a.y*b.y+a.z*b.z;
}
static inline Vec3 cross3(const Vec3&a,const Vec3&b){
return{
a.y*b.z-a.z*b.y,a.z*b.x-a.x*b.z,a.x*b.y-a.y*b.x}
;
}
static inline double norm2(const Vec3&a){
return dot3(a,a);
}
static inline double dist2(const Vec3&a,const Vec3&b){
return norm2(a-b);
}
static inline bool normalize(Vec3&a){
double n2=norm2(a);
if(n2<=1e-30)return false;
double inv=1.0/sqrt(n2);
a=a*inv;
return true;
}
struct Face{
int v[3];
unsigned char active=1;
}
;
struct Quadric{
double q[10];
Quadric(){
memset(q,0,sizeof(q));
}
void add_plane(double a,double b,double c,double d,double w=1.0){
q[0]+=w*a*a;
q[1]+=w*a*b;
q[2]+=w*a*c;
q[3]+=w*a*d;
q[4]+=w*b*b;
q[5]+=w*b*c;
q[6]+=w*b*d;
q[7]+=w*c*c;
q[8]+=w*c*d;
q[9]+=w*d*d;
}
void add(const Quadric&o){
for(int i=0;
i<10;
++i)q[i]+=o.q[i];
}
double eval(const Vec3&p)const{
double x=p.x,y=p.y,z=p.z;
return q[0]*x*x+2.0*q[1]*x*y+2.0*q[2]*x*z+2.0*q[3]*x+q[4]*y*y+2.0*q[5]*y*z+2.0*q[6]*y+q[7]*z*z+2.0*q[8]*z+q[9];
}
}
;
struct FastInput{
vector<char>buf;
char*p=nullptr;
void read_all(){
char tmp[1<<16];
size_t n;
while((n=fread(tmp,1,sizeof(tmp),stdin))>0){
buf.insert(buf.end(),tmp,tmp+n);
}
buf.push_back('\0');
p=buf.data();
}
inline void skip_ws(){
while(*p==' '||*p=='\n'||*p=='\r'||*p=='\t')++p;
}
int next_int(){
skip_ws();
int sign=1;
if(*p=='-'){
sign=-1;
++p;
}
int x=0;
while(*p>='0'&&*p<='9'){
x=x*10+(*p-'0');
++p;
}
return x*sign;
}
double next_double(){
skip_ws();
char*e=nullptr;
double x=strtod(p,&e);
p=e;
return x;
}
char next_char(){
skip_ws();
return*p++;
}
}
;
struct EdgeFace{
uint64_t key;
int face;
bool operator<(const EdgeFace&o)const{
if(key!=o.key)return key<o.key;
return face<o.face;
}
}
;
struct Node{
double cost;
int u,v;
int vu,vv;
bool operator<(const Node&o)const{
return cost>o.cost;
}
}
;
static int N,M;
static vector<Vec3>Orig;
static vector<Vec3>P;
static vector<Face>F;
static vector<vector<int>>Y;
static vector<Quadric>Q;
static vector<unsigned char>active_vertex;
static vector<int>DH,DI,next_member,cluster_size;
static vector<int>DW;
static vector<int>BZ;
static int mark_token=7;
static int active_vertices=0;
static int BE=0;
static double hausdorff_limit2=0.0;
static double diagonal_len=0.0;
static int AV=0;
static double target_ratio_override=-1.0;
static priority_queue<Node>pq;
static chrono::steady_clock::time_point start_time;
static double TIME_LIMIT_SECONDS=19.25;
static const double MIN_NORMAL_COS=0.55;
static inline uint64_t edge_key(int a,int b){
if(a>b)swap(a,b);
return(uint64_t)(uint32_t)a<<32|(uint32_t)b;
}
static inline int key_a(uint64_t key){
return(int)(key>>32);
}
static inline int key_b(uint64_t key){
return(int)(key&0xffffffffu);
}
static inline bool time_left(){
auto now=chrono::steady_clock::now();
double elapsed=chrono::duration<double>(now-start_time).count();
return elapsed<TIME_LIMIT_SECONDS;
}
static inline bool AC(const Face&f,int v){
return f.v[0]==v||f.v[1]==v||f.v[2]==v;
}
static inline bool face_has_edge(const Face&f,int a,int b){
return AC(f,a)&&AC(f,b);
}
static inline int third_vertex(const Face&f,int a,int b){
for(int k=0;
k<3;
++k){
int x=f.v[k];
if(x!=a&&x!=b)return x;
}
return-1;
}
static inline Vec3 face_normal_from_indices(int a,int b,int c){
return cross3(P[b]-P[a],P[c]-P[a]);
}
static bool are_adjacent(int a,int b){
const vector<int>&ia=Y[a];
const vector<int>&ib=Y[b];
const vector<int>&small=(ia.size()<ib.size())?ia:ib;
for(int fid:small){
if(!F[fid].active)continue;
if(face_has_edge(F[fid],a,b))return true;
}
return false;
}
static bool collect_edge_opposites(int a,int b,int opp[2]){
int cnt=0;
const vector<int>&ia=Y[a];
const vector<int>&ib=Y[b];
const vector<int>&small=(ia.size()<ib.size())?ia:ib;
for(int fid:small){
if(!F[fid].active)continue;
const Face&f=F[fid];
if(!face_has_edge(f,a,b))continue;
if(cnt>=2)return false;
int t=third_vertex(f,a,b);
if(t<0)return false;
opp[cnt++]=t;
}
if(cnt!=2)return false;
return opp[0]!=opp[1];
}
static bool ET(int a,int b){
int opp[2];
if(!collect_edge_opposites(a,b,opp))return false;
if(mark_token>2000000000){
fill(BZ.begin(),BZ.end(),0);
mark_token=7;
}
int token_u=mark_token++;
int token_common=mark_token++;
for(int fid:Y[a]){
if(!F[fid].active)continue;
const Face&f=F[fid];
if(!AC(f,a))continue;
for(int k=0;
k<3;
++k){
int x=f.v[k];
if(x!=a&&x!=b)BZ[x]=token_u;
}
}
int common_count=0;
int got0=0,got1=0;
for(int fid:Y[b]){
if(!F[fid].active)continue;
const Face&f=F[fid];
if(!AC(f,b))continue;
for(int k=0;
k<3;
++k){
int x=f.v[k];
if(x==a||x==b)continue;
if(BZ[x]==token_u){
BZ[x]=token_common;
++common_count;
if(x==opp[0])got0=1;
if(x==opp[1])got1=1;
if(common_count>2)return false;
}
}
}
return common_count==2&&got0&&got1;
}
static bool cluster_can_move_to(int v,const Vec3&to){
for(int m=DH[v];
m!=-1;
m=next_member[m]){
if(dist2(Orig[m],to)>hausdorff_limit2)return false;
}
return true;
}
static bool normals_ok_after_collapse(int a,int b,const Vec3&to){
auto scan=[&](int src)->bool{
for(int fid:Y[src]){
if(!F[fid].active)continue;
const Face&f=F[fid];
bool has_a=AC(f,a);
bool has_b=AC(f,b);
if(!has_a&&!has_b)continue;
if(has_a&&has_b)continue;
Vec3 old_p[3]={
P[f.v[0]],P[f.v[1]],P[f.v[2]]}
;
Vec3 new_p[3]={
old_p[0],old_p[1],old_p[2]}
;
for(int k=0;
k<3;
++k){
if(f.v[k]==a||f.v[k]==b)new_p[k]=to;
}
Vec3 old_n=cross3(old_p[1]-old_p[0],old_p[2]-old_p[0]);
Vec3 new_n=cross3(new_p[1]-new_p[0],new_p[2]-new_p[0]);
double old_len2=norm2(old_n);
double new_len2=norm2(new_n);
if(old_len2<=1e-30||new_len2<=1e-30)return false;
double d=dot3(old_n,new_n);
double limit=MIN_NORMAL_COS*sqrt(old_len2*new_len2);
if(d<=limit)return false;
}
return true;
}
;
return scan(a)&&scan(b);
}
static bool solve_optimal_position(const Quadric&q,Vec3&out){
double a00=q.q[0],a01=q.q[1],a02=q.q[2];
double a10=q.q[1],a11=q.q[4],a12=q.q[5];
double a20=q.q[2],a21=q.q[5],a22=q.q[7];
double b0=-q.q[3],b1=-q.q[6],b2=-q.q[8];
double det=a00*(a11*a22-a12*a21)-a01*(a10*a22-a12*a20)+a02*(a10*a21-a11*a20);
if(fabs(det)<1e-14)return false;
double dx=b0*(a11*a22-a12*a21)-a01*(b1*a22-a12*b2)+a02*(b1*a21-a11*b2);
double dy=a00*(b1*a22-a12*b2)-b0*(a10*a22-a12*a20)+a02*(a10*b2-b1*a20);
double dz=a00*(a11*b2-b1*a21)-a01*(a10*b2-b1*a20)+b0*(a10*a21-a11*a20);
out={
dx/det,dy/det,dz/det}
;
return isfinite(out.x)&&isfinite(out.y)&&isfinite(out.z);
}
static bool candidate_ok(int a,int b,const Vec3&pos){
if(!cluster_can_move_to(a,pos))return false;
if(!cluster_can_move_to(b,pos))return false;
return normals_ok_after_collapse(a,b,pos);
}
static bool best_collapse_position(int a,int b,Vec3&best_pos,double&best_cost){
Quadric q=Q[a];
q.add(Q[b]);
Vec3 opt;
Vec3 cand[6];
int cnt=0;
if(solve_optimal_position(q,opt))cand[cnt++]=opt;
cand[cnt++]=(P[a]+P[b])*0.5;
cand[cnt++]=P[a];
cand[cnt++]=P[b];
cand[cnt++]=P[a]*0.75+P[b]*0.25;
cand[cnt++]=P[a]*0.25+P[b]*0.75;
best_cost=1e100;
bool ok=false;
for(int i=0;
i<cnt;
++i){
const Vec3&pos=cand[i];
if(!candidate_ok(a,b,pos))continue;
double c=q.eval(pos)+0.0003*(dist2(pos,P[a])+dist2(pos,P[b]));
if(c<best_cost){
best_cost=c;
best_pos=pos;
ok=true;
}
}
return ok;
}
static double cheap_edge_cost(int a,int b){
Quadric q=Q[a];
q.add(Q[b]);
Vec3 opt;
double best=1e100;
if(solve_optimal_position(q,opt))best=min(best,q.eval(opt));
best=min(best,q.eval((P[a]+P[b])*0.5));
best=min(best,q.eval(P[a]));
best=min(best,q.eval(P[b]));
return best+0.0003*dist2(P[a],P[b]);
}
static void push_edge(int a,int b){
if(a==b)return;
if(!active_vertex[a]||!active_vertex[b])return;
double d2=dist2(P[a],P[b]);
if(d2>4.00001*hausdorff_limit2)return;
double c=cheap_edge_cost(a,b);
pq.push({
c,a,b,DW[a],DW[b]}
);
}
static void compact_incident(int v){
vector<int>&ids=Y[v];
if(ids.size()<128)return;
size_t alive=0;
for(int fid:ids){
if(F[fid].active&&AC(F[fid],v))++alive;
}
if(alive*3+32>=ids.size())return;
vector<int>keep;
keep.reserve(alive+8);
for(int fid:ids){
if(F[fid].active&&AC(F[fid],v))keep.push_back(fid);
}
ids.swap(keep);
}
static void merge_members(int src,int dst){
if(DH[src]==-1)return;
next_member[DI[dst]]=DH[src];
DI[dst]=DI[src];
cluster_size[dst]+=cluster_size[src];
DH[src]=DI[src]=-1;
cluster_size[src]=0;
}
static void do_collapse(int src,int dst,const Vec3&pos){
Q[dst].add(Q[src]);
P[dst]=pos;
for(int fid:Y[src]){
if(!F[fid].active)continue;
Face&f=F[fid];
bool has_src=false,has_dst=false;
for(int k=0;
k<3;
++k){
if(f.v[k]==src)has_src=true;
if(f.v[k]==dst)has_dst=true;
}
if(!has_src)continue;
if(has_dst){
f.active=0;
--BE;
}
else{
for(int k=0;
k<3;
++k){
if(f.v[k]==src)f.v[k]=dst;
}
Y[dst].push_back(fid);
}
}
active_vertex[src]=0;
--active_vertices;
++DW[src];
++DW[dst];
merge_members(src,dst);
compact_incident(src);
compact_incident(dst);
for(int fid:Y[dst]){
if(!F[fid].active)continue;
const Face&f=F[fid];
if(!AC(f,dst))continue;
push_edge(f.v[0],f.v[1]);
push_edge(f.v[1],f.v[2]);
push_edge(f.v[2],f.v[0]);
}
}
static bool attempt_collapse(int a,int b){
if(a==b)return false;
if(!active_vertex[a]||!active_vertex[b])return false;
if(dist2(P[a],P[b])>4.00001*hausdorff_limit2)return false;
if(!ET(a,b))return false;
Vec3 pos;
double cost=0.0;
if(!best_collapse_position(a,b,pos,cost))return false;
int src,dst;
size_t wa=Y[a].size()+(size_t)cluster_size[a]*2;
size_t wb=Y[b].size()+(size_t)cluster_size[b]*2;
if(wa<=wb){
src=a;
dst=b;
}
else{
src=b;
dst=a;
}
do_collapse(src,dst,pos);
return true;
}
static void JC(){
FastInput in;
in.read_all();
N=in.next_int();
M=in.next_int();
P.resize(N);
Orig.resize(N);
F.resize(M);
Vec3 mn{
1e100,1e100,1e100}
;
Vec3 mx{
-1e100,-1e100,-1e100}
;
for(int i=0;
i<N;
++i){
(void)in.next_char();
P[i].x=in.next_double();
P[i].y=in.next_double();
P[i].z=in.next_double();
Orig[i]=P[i];
mn.x=min(mn.x,P[i].x);
mn.y=min(mn.y,P[i].y);
mn.z=min(mn.z,P[i].z);
mx.x=max(mx.x,P[i].x);
mx.y=max(mx.y,P[i].y);
mx.z=max(mx.z,P[i].z);
}
vector<int>deg(N,0);
for(int i=0;
i<M;
++i){
(void)in.next_char();
int a=in.next_int()-1;
int b=in.next_int()-1;
int c=in.next_int()-1;
F[i].v[0]=a;
F[i].v[1]=b;
F[i].v[2]=c;
++deg[a];
++deg[b];
++deg[c];
}
Vec3 d=mx-mn;
diagonal_len=sqrt(norm2(d));
double hausdorff_limit=0.05*diagonal_len*0.999999;
hausdorff_limit2=hausdorff_limit*hausdorff_limit;
Y.resize(N);
for(int i=0;
i<N;
++i)Y[i].reserve(deg[i]+8);
for(int i=0;
i<M;
++i){
Y[F[i].v[0]].push_back(i);
Y[F[i].v[1]].push_back(i);
Y[F[i].v[2]].push_back(i);
}
Q.assign(N,Quadric());
active_vertex.assign(N,1);
DH.resize(N);
DI.resize(N);
next_member.assign(N,-1);
cluster_size.assign(N,1);
DW.assign(N,0);
BZ.assign(N,0);
for(int i=0;
i<N;
++i){
DH[i]=DI[i]=i;
}
active_vertices=N;
BE=M;
}
static double initialize_quadrics_and_edges(){
vector<Vec3>face_normals(M);
vector<EdgeFace>CP;
CP.reserve((size_t)M*3);
for(int i=0;
i<M;
++i){
Face&f=F[i];
Vec3 n=face_normal_from_indices(f.v[0],f.v[1],f.v[2]);
if(!normalize(n))n={
0.0,0.0,0.0}
;
face_normals[i]=n;
double dd=-dot3(n,P[f.v[0]]);
Q[f.v[0]].add_plane(n.x,n.y,n.z,dd);
Q[f.v[1]].add_plane(n.x,n.y,n.z,dd);
Q[f.v[2]].add_plane(n.x,n.y,n.z,dd);
CP.push_back({
edge_key(f.v[0],f.v[1]),i}
);
CP.push_back({
edge_key(f.v[1],f.v[2]),i}
);
CP.push_back({
edge_key(f.v[2],f.v[0]),i}
);
}
sort(CP.begin(),CP.end());
long long unique_edges=0;
long long feature_edges=0;
const double feature_cos=cos(35.0*acos(-1.0)/180.0);
for(size_t i=0;
i<CP.size();
){
size_t j=i+1;
while(j<CP.size()&&CP[j].key==CP[i].key)++j;
++unique_edges;
if(j-i==2){
double d=dot3(face_normals[CP[i].face],face_normals[CP[i+1].face]);
if(d<feature_cos)++feature_edges;
}
int a=key_a(CP[i].key);
int b=key_b(CP[i].key);
push_edge(a,b);
i=j;
}
vector<EdgeFace>().swap(CP);
vector<Vec3>().swap(face_normals);
if(unique_edges==0)return 0.0;
return(double)feature_edges/(double)unique_edges;
}
static void choose_target(double feature_ratio){
double ratio=0.089+0.035*min(0.22,feature_ratio);
if(N<=8000)ratio=max(ratio,0.095);
if(N<=1000)ratio=max(ratio,0.160);
if(N<=50)ratio=0.01;
if(target_ratio_override>0.0){
ratio=max(0.025,min(0.500,target_ratio_override));
}
else{
ratio=max(0.086,min(0.115,ratio));
}
AV=max(4,(int)ceil((double)N*ratio));
}
static void GN(){
double feature_ratio=initialize_quadrics_and_edges();
choose_target(feature_ratio);
long long pops=0;
while(active_vertices>AV&&!pq.empty()){
if((++pops&4095)==0&&!time_left())break;
Node cur=pq.top();
pq.pop();
int a=cur.u,b=cur.v;
if(a==b)continue;
if(!active_vertex[a]||!active_vertex[b])continue;
if(cur.vu!=DW[a]||cur.vv!=DW[b]){
if(are_adjacent(a,b))push_edge(a,b);
continue;
}
if(!attempt_collapse(a,b))continue;
}
}
static void JD(){
vector<int>remap(N,-1);
int BN=0;
for(int i=0;
i<N;
++i){
if(active_vertex[i])remap[i]=BN++;
}
int AB=0;
for(int i=0;
i<M;
++i){
if(!F[i].active)continue;
int a=F[i].v[0],b=F[i].v[1],c=F[i].v[2];
if(a==b||b==c||c==a)continue;
if(remap[a]<0||remap[b]<0||remap[c]<0)continue;
++AB;
}
static char outbuf[1<<20];
setvbuf(stdout,outbuf,_IOFBF,sizeof(outbuf));
printf("%d %d\n",BN,AB);
for(int i=0;
i<N;
++i){
if(!active_vertex[i])continue;
printf("v %.15g %.15g %.15g\n",P[i].x,P[i].y,P[i].z);
}
for(int i=0;
i<M;
++i){
if(!F[i].active)continue;
int a=F[i].v[0],b=F[i].v[1],c=F[i].v[2];
if(a==b||b==c||c==a)continue;
if(remap[a]<0||remap[b]<0||remap[c]<0)continue;
printf("f %d %d %d\n",remap[a]+1,remap[b]+1,remap[c]+1);
}
}
static bool EV(const vector<::Vec3>&DQ,const vector<::Face>&source_faces,double parent_diag,vector<::Vec3>&BN,vector<::Face>&AB,double seconds_limit,double target_ratio=-1.0){
start_time=chrono::steady_clock::now();
TIME_LIMIT_SECONDS=max(0.25,seconds_limit);
target_ratio_override=target_ratio;
N=(int)DQ.size();
M=(int)source_faces.size();
if(N<4||M<4)return false;
Orig.assign(N,Vec3{
}
);
P.assign(N,Vec3{
}
);
F.assign(M,Face{
}
);
Vec3 mn{
1e100,1e100,1e100}
;
Vec3 mx{
-1e100,-1e100,-1e100}
;
for(int i=0;
i<N;
++i){
P[i]={
DQ[i].x,DQ[i].y,DQ[i].z}
;
Orig[i]=P[i];
mn.x=min(mn.x,P[i].x);
mn.y=min(mn.y,P[i].y);
mn.z=min(mn.z,P[i].z);
mx.x=max(mx.x,P[i].x);
mx.y=max(mx.y,P[i].y);
mx.z=max(mx.z,P[i].z);
}
vector<int>deg(N,0);
for(int i=0;
i<M;
++i){
F[i].v[0]=source_faces[i].v[0];
F[i].v[1]=source_faces[i].v[1];
F[i].v[2]=source_faces[i].v[2];
F[i].active=1;
if(F[i].v[0]<0||F[i].v[0]>=N||F[i].v[1]<0||F[i].v[1]>=N||F[i].v[2]<0||F[i].v[2]>=N)return false;
++deg[F[i].v[0]];
++deg[F[i].v[1]];
++deg[F[i].v[2]];
}
Vec3 d=mx-mn;
diagonal_len=parent_diag>0.0?parent_diag:sqrt(norm2(d));
double hausdorff_limit=0.05*diagonal_len*0.999999;
hausdorff_limit2=hausdorff_limit*hausdorff_limit;
Y.assign(N,{
}
);
for(int i=0;
i<N;
++i)Y[i].reserve(deg[i]+8);
for(int i=0;
i<M;
++i){
Y[F[i].v[0]].push_back(i);
Y[F[i].v[1]].push_back(i);
Y[F[i].v[2]].push_back(i);
}
Q.assign(N,Quadric());
active_vertex.assign(N,1);
DH.resize(N);
DI.resize(N);
next_member.assign(N,-1);
cluster_size.assign(N,1);
DW.assign(N,0);
BZ.assign(N,0);
mark_token=7;
for(int i=0;
i<N;
++i)DH[i]=DI[i]=i;
active_vertices=N;
BE=M;
pq=priority_queue<Node>();
GN();
if(active_vertices<=0||BE<=0)return false;
vector<int>remap(N,-1);
BN.clear();
AB.clear();
BN.reserve(active_vertices);
for(int i=0;
i<N;
++i){
if(!active_vertex[i])continue;
remap[i]=(int)BN.size();
BN.push_back(::Vec3{
P[i].x,P[i].y,P[i].z}
);
}
for(int i=0;
i<M;
++i){
if(!F[i].active)continue;
int a=F[i].v[0],b=F[i].v[1],c=F[i].v[2];
if(a==b||b==c||c==a)continue;
if(a<0||b<0||c<0||a>=N||b>=N||c>=N)continue;
if(remap[a]<0||remap[b]<0||remap[c]<0)continue;
::Face f;
f.v[0]=remap[a];
f.v[1]=remap[b];
f.v[2]=remap[c];
AB.push_back(f);
}
return!BN.empty()&&!AB.empty()&&BN.size()<DQ.size();
}
}
namespace QX{
using namespace std;
struct FastScanner {

    static const size_t BUFSIZE = 1 << 20;

    int idx = 0, size = 0;

    char buf[BUFSIZE];

    inline char getch() {

        if (idx >= size) {

            size = (int)fread(buf, 1, BUFSIZE, stdin);

            idx = 0;

            if (!size) return 0;

        }

        return buf[idx++];

    }

    template<class T> bool readInt(T &out) {

        char c = getch();

        while (c && (c==' ' || c=='\n' || c=='\r' || c=='\t')) c = getch();

        if (!c) return false;

        T sign = 1, x = 0;

        if (c=='-') sign = -1, c = getch();

        while (c>='0' && c<='9') {
 x = x*10 + (c-'0');
 c = getch();
 }

        out = x * sign;

        return true;

    }

    inline bool readPrefix(char want) {

        char c = getch();

        while (c && (c==' ' || c=='\n' || c=='\r' || c=='\t')) c = getch();

        return c == want;

    }

    bool readDouble(double &out) {

        char c = getch();

        while (c && (c==' ' || c=='\n' || c=='\r' || c=='\t')) c = getch();

        if (!c) return false;

        int sign = 1;

        if (c=='-') sign = -1, c = getch();

        double x = 0.0;

        while (c>='0' && c<='9') {
 x = x*10.0 + (double)(c-'0');
 c = getch();
 }

        if (c=='.') {

            double base = 0.1;
 c = getch();

            while (c>='0' && c<='9') {
 x += base * (double)(c-'0');
 base *= 0.1;
 c = getch();
 }

        }

        if (c=='e' || c=='E') {

            int esign=1, expv=0;
 c = getch();

            if (c=='-') esign=-1, c=getch();
 else if (c=='+') c=getch();

            while (c>='0' && c<='9') {
 expv = expv*10 + (c-'0');
 c=getch();
 }

            x *= pow(10.0, esign * expv);

        }

        out = sign * x;

        return true;

    }

}
;


struct Vec3 {

    double x, y, z;

    Vec3():x(0),y(0),z(0){
}

    Vec3(double X,double Y,double Z):x(X),y(Y),z(Z){
}

    inline Vec3 operator+(const Vec3& o) const {
 return Vec3(x+o.x,y+o.y,z+o.z);
 }

    inline Vec3 operator-(const Vec3& o) const {
 return Vec3(x-o.x,y-o.y,z-o.z);
 }

    inline Vec3 operator*(double s) const {
 return Vec3(x*s,y*s,z*s);
 }

    inline Vec3 operator/(double s) const {
 return Vec3(x/s,y/s,z/s);
 }

}
;

struct Tri {
 int a,b,c;
 }
;

static inline double dotv(const Vec3& a,const Vec3& b){
 return a.x*b.x+a.y*b.y+a.z*b.z;
 }

static inline Vec3 crossv(const Vec3& a,const Vec3& b){
 return Vec3(a.y*b.z-a.z*b.y, a.z*b.x-a.x*b.z, a.x*b.y-a.y*b.x);
 }

static inline double norm2(const Vec3& a){
 return dotv(a,a);
 }

static inline double normv(const Vec3& a){
 return sqrt(norm2(a));
 }


struct Quadric {

    
    double q[10];

    Quadric(){
 memset(q,0,sizeof(q));
 }

    inline void add(const Quadric& o){
 for(int i=0;
i<10;
i++) q[i]+=o.q[i];
 }

    static Quadric plane(double a,double b,double c,double d,double w) {

        Quadric Q;

        Q.q[0]=w*a*a;
 Q.q[1]=w*a*b;
 Q.q[2]=w*a*c;
 Q.q[3]=w*a*d;

        Q.q[4]=w*b*b;
 Q.q[5]=w*b*c;
 Q.q[6]=w*b*d;

        Q.q[7]=w*c*c;
 Q.q[8]=w*c*d;
 Q.q[9]=w*d*d;

        return Q;

    }

    inline double eval(const Vec3& p) const {

        double x=p.x,y=p.y,z=p.z;

        return q[0]*x*x + 2*q[1]*x*y + 2*q[2]*x*z + 2*q[3]*x
             + q[4]*y*y + 2*q[5]*y*z + 2*q[6]*y
             + q[7]*z*z + 2*q[8]*z + q[9];

    }

}
;

static inline Quadric qsum(const Quadric& a,const Quadric& b){
 Quadric r=a;
 r.add(b);
 return r;
 }


struct Face {

    int v[3];

    unsigned char active;

    Vec3 n;

}
;

struct Candidate {

    double score;

    int u, v, vu, vv;

    bool operator<(Candidate const& o) const {
 return score > o.score;
 }

}
;

struct EvalCand {
 bool ok;
 double score;
 double rawError;
 Vec3 pos;
 }
;

struct Snapshot {
 vector<Vec3> V;
 vector<Tri> F;
 int nv=0;
 double ratio=1.0;
 double ssim=-1.0;
 }
;

struct RenderMaps {
 int R=0;
 vector<float> depth;
 vector<float> norm;
 }
;


static int N, M;

static vector<Vec3> P;

static vector<Face> F;

static vector<vector<int>> inc;

static vector<Quadric> Qv;

static vector<array<float,3>> bbMin, bbMax;

static vector<unsigned char> alive;

static vector<int> ver, markA, markB, tmpId;

static int stampA=1, stampB=1;

static double hausTau=0.0, hausTau2=0.0;

static int activeVertices=0, activeFaces=0;

static double visualPixelTol=18.0, visualPixelTol2=324.0;

static vector<Snapshot> snapshots;

static Snapshot finalSnap;

static bool haveFinalSnap=false;


static inline bool face_contains(int fid,int u){
 const Face& f=F[fid];
 return f.v[0]==u || f.v[1]==u || f.v[2]==u;
 }

static inline Vec3 computeFaceNormal(const Face& f){

    Vec3 cr = crossv(P[f.v[1]]-P[f.v[0]], P[f.v[2]]-P[f.v[0]]);

    double l=normv(cr);
 if(l<=1e-300) return Vec3(0,0,0);
 return cr/l;

}


static void cleanupIncident(int u){

    if(!alive[u]) return;

    vector<int>& lst=inc[u];
 int wr=0;

    for(int fid: lst) if(fid>=0 && fid<M && F[fid].active && face_contains(fid,u)) lst[wr++]=fid;

    lst.resize(wr);

}

static void maybeCleanupIncident(int u){

    if(!alive[u]) return;

    if(inc[u].size()>96){

        int good=0;
 for(int fid: inc[u]) if(F[fid].active && face_contains(fid,u)) ++good;

        if((int)inc[u].size()>good*3+64) cleanupIncident(u);

    }

}

static inline void mergedBBox(int u,int v,array<float,3>& mn,array<float,3>& mx){

    mn[0]=min(bbMin[u][0],bbMin[v][0]);
 mn[1]=min(bbMin[u][1],bbMin[v][1]);
 mn[2]=min(bbMin[u][2],bbMin[v][2]);

    mx[0]=max(bbMax[u][0],bbMax[v][0]);
 mx[1]=max(bbMax[u][1],bbMax[v][1]);
 mx[2]=max(bbMax[u][2],bbMax[v][2]);

}

static inline bool coversBBox(const Vec3& p,const array<float,3>& mn,const array<float,3>& mx){

    double maxd2=0.0;

    for(int mask=0;
mask<8;
mask++){

        double x=(mask&1)?mx[0]:mn[0], y=(mask&2)?mx[1]:mn[1], z=(mask&4)?mx[2]:mn[2];

        double dx=p.x-x,dy=p.y-y,dz=p.z-z;
 double d2=dx*dx+dy*dy+dz*dz;

        if(d2>maxd2) maxd2=d2;

    }

    return maxd2 <= hausTau2;

}

static inline void project1024(const Vec3& p,int view,double& u,double& v,double& dep){

    int ax=view/2;
 int sg = (view&1)?-1:1;
 
    double sx, sy;

    if(ax==0){
 sx=p.y;
 sy=p.z;
 dep=2.5 - sg*p.x;
 if(sg<0) sx=-sx;
 }

    else if(ax==1){
 sx=p.x;
 sy=p.z;
 dep=2.5 - sg*p.y;
 if(sg>0) sx=-sx;
 }

    else {
 sx=p.x;
 sy=p.y;
 dep=2.5 - sg*p.z;
 if(sg<0) sx=-sx;
 }

    u = 800.0*sx/dep + 512.0;
 v = 800.0*sy/dep + 512.0;

}

static inline bool visualBBoxOK(const Vec3& p,const array<float,3>& mn,const array<float,3>& mx){

    
    double pu[6], pv[6], pd;

    for(int view=0;
view<6;
view++) project1024(p,view,pu[view],pv[view],pd);

    double maxd2=0.0;

    for(int mask=0;
mask<8;
mask++){

        Vec3 c((mask&1)?mx[0]:mn[0], (mask&2)?mx[1]:mn[1], (mask&4)?mx[2]:mn[2]);

        for(int view=0;
view<6;
view++){

            double u,v,d;
 project1024(c,view,u,v,d);

            double du=u-pu[view], dv=v-pv[view];
 double d2=du*du+dv*dv;

            if(d2>maxd2) maxd2=d2;

            if(maxd2>visualPixelTol2) return false;

        }

    }

    return true;

}


static bool solveOptimal(const Quadric& q,Vec3& out){

    double a00=q.q[0], a01=q.q[1], a02=q.q[2];

    double a11=q.q[4], a12=q.q[5], a22=q.q[7];

    double b0=q.q[3], b1=q.q[6], b2=q.q[8];

    double det=a00*(a11*a22-a12*a12)-a01*(a01*a22-a12*a02)+a02*(a01*a12-a11*a02);

    if(fabs(det)<1e-14) return false;

    double inv00=(a11*a22-a12*a12)/det;

    double inv01=(a02*a12-a01*a22)/det;

    double inv02=(a01*a12-a02*a11)/det;

    double inv11=(a00*a22-a02*a02)/det;

    double inv12=(a01*a02-a00*a12)/det;

    double inv22=(a00*a11-a01*a01)/det;

    out.x=-(inv00*b0+inv01*b1+inv02*b2);

    out.y=-(inv01*b0+inv11*b1+inv12*b2);

    out.z=-(inv02*b0+inv12*b1+inv22*b2);

    if(!isfinite(out.x)||!isfinite(out.y)||!isfinite(out.z)) return false;

    return out.x>=-2.0&&out.x<=2.0&&out.y>=-2.0&&out.y<=2.0&&out.z>=-2.0&&out.z<=2.0;

}

static EvalCand computeCandidate(int u,int v){

    EvalCand ec{
false,1e300,1e300,Vec3()}
;

    if(u==v || !alive[u] || !alive[v]) return ec;

    if(u>v) swap(u,v);

    array<float,3> mn,mx;
 mergedBBox(u,v,mn,mx);

    Quadric q=qsum(Qv[u],Qv[v]);

    Vec3 cand[7];
 int cc=0;
 Vec3 opt;

    if(solveOptimal(q,opt)) cand[cc++]=opt;

    cand[cc++]=P[u];
 cand[cc++]=P[v];

    cand[cc++]=(P[u]+P[v])*0.5;

    cand[cc++]=Vec3((mn[0]+mx[0])*0.5,(mn[1]+mx[1])*0.5,(mn[2]+mx[2])*0.5);

    
    cand[cc++]=P[u]*0.67 + P[v]*0.33;

    cand[cc++]=P[u]*0.33 + P[v]*0.67;

    double len2=norm2(P[u]-P[v]);

    for(int i=0;
i<cc;
i++){

        Vec3 p=cand[i];

        if(!coversBBox(p,mn,mx)) continue;

        double e=q.eval(p);
 if(e<0 && e>-1e-12) e=0;

        
        
        
        
        double sc=e + 2e-10*len2 + 5e-15*norm2(p);

        if(sc<ec.score){
 ec.ok=true;
 ec.score=sc;
 ec.rawError=e;
 ec.pos=p;
 }

    }

    return ec;

}

static bool linkConditionOK(int u,int v){

    if(u==v || !alive[u] || !alive[v]) return false;

    maybeCleanupIncident(u);
 maybeCleanupIncident(v);

    if(++stampA==INT_MAX){
 fill(markA.begin(),markA.end(),0);
 stampA=1;
 }

    if(++stampB==INT_MAX){
 fill(markB.begin(),markB.end(),0);
 stampB=1;
 }

    int edgeFaces=0;

    for(int fid: inc[u]){

        Face& f=F[fid];
 if(!f.active || !face_contains(fid,u)) continue;

        bool hasv=false;

        for(int k=0;
k<3;
k++){
 int w=f.v[k];
 if(w==v) hasv=true;
 if(w!=u) markA[w]=stampA;
 }

        if(hasv) ++edgeFaces;

    }

    if(edgeFaces!=2) return false;

    int common=0;

    for(int fid: inc[v]){

        Face& f=F[fid];
 if(!f.active || !face_contains(fid,v)) continue;

        for(int k=0;
k<3;
k++){

            int w=f.v[k];
 if(w==u || w==v) continue;

            if(markA[w]==stampA && markB[w]!=stampB){
 markB[w]=stampB;
 if(++common>2) return false;
 }

        }

    }

    return common==2;

}

static bool flipOK(int keep,int rem,const Vec3& newPos,double minDot){

    static vector<int> touched;
 touched.clear();

    for(int fid: inc[keep]) if(F[fid].active && face_contains(fid,keep)) touched.push_back(fid);

    for(int fid: inc[rem]) if(F[fid].active && face_contains(fid,rem)) touched.push_back(fid);

    sort(touched.begin(),touched.end());
 touched.erase(unique(touched.begin(),touched.end()),touched.end());

    const double minArea2=1e-28;

    for(int fid: touched){

        Face& f=F[fid];
 if(!f.active) continue;

        bool hasK=false,hasR=false;
 for(int k=0;
k<3;
k++){
 if(f.v[k]==keep) hasK=true;
 if(f.v[k]==rem) hasR=true;
 }

        if(hasK && hasR) continue;

        int a=f.v[0],b=f.v[1],c=f.v[2];

        Vec3 pa=(a==keep||a==rem)?newPos:P[a];

        Vec3 pb=(b==keep||b==rem)?newPos:P[b];

        Vec3 pc=(c==keep||c==rem)?newPos:P[c];

        Vec3 cr=crossv(pb-pa,pc-pa);
 double a2=norm2(cr);

        if(!(a2>minArea2) || !isfinite(a2)) return false;

        Vec3 nn=cr/sqrt(a2);

        if(dotv(nn,f.n)<minDot) return false;

    }

    return true;

}

static void pushEdge(priority_queue<Candidate>& pq,int u,int v){

    if(u==v || !alive[u] || !alive[v]) return;
 if(u>v) swap(u,v);

    EvalCand ec=computeCandidate(u,v);
 if(!ec.ok) return;

    pq.push({
ec.score,u,v,ver[u],ver[v]}
);

}

static void collectNeighbors(int u,vector<int>& out){

    out.clear();
 if(!alive[u]) return;
 maybeCleanupIncident(u);

    if(++stampA==INT_MAX){
 fill(markA.begin(),markA.end(),0);
 stampA=1;
 }

    for(int fid: inc[u]){

        Face& f=F[fid];
 if(!f.active || !face_contains(fid,u)) continue;

        for(int k=0;
k<3;
k++){
 int w=f.v[k];
 if(w!=u && alive[w] && markA[w]!=stampA){
 markA[w]=stampA;
 out.push_back(w);
 }
 }

    }

}

static void doCollapse(int keep,int rem,const Vec3& newPos,priority_queue<Candidate>& pq){

    cleanupIncident(keep);
 cleanupIncident(rem);

    vector<int> remFaces=inc[rem];

    for(int fid: remFaces){

        Face& f=F[fid];
 if(!f.active || !face_contains(fid,rem)) continue;

        bool hasKeep=false;
 for(int k=0;
k<3;
k++) if(f.v[k]==keep) hasKeep=true;

        if(hasKeep){
 f.active=0;
 --activeFaces;
 }

        else{

            for(int k=0;
k<3;
k++) if(f.v[k]==rem) f.v[k]=keep;

            if(f.v[0]==f.v[1] || f.v[1]==f.v[2] || f.v[2]==f.v[0]){
 f.active=0;
 --activeFaces;
 }

            else{
 P[keep]=newPos;
 f.n=computeFaceNormal(f);
 inc[keep].push_back(fid);
 }

        }

    }

    P[keep]=newPos;

    for(int fid: inc[keep]){
 Face& f=F[fid];
 if(f.active && face_contains(fid,keep)) f.n=computeFaceNormal(f);
 }

    Qv[keep].add(Qv[rem]);

    for(int d=0;
d<3;
d++){
 bbMin[keep][d]=min(bbMin[keep][d],bbMin[rem][d]);
 bbMax[keep][d]=max(bbMax[keep][d],bbMax[rem][d]);
 }

    alive[rem]=0;
 ++ver[keep];
 ++ver[rem];
 --activeVertices;
 inc[rem].clear();
 cleanupIncident(keep);

    static vector<int> neigh;
 collectNeighbors(keep,neigh);
 for(int w: neigh) pushEdge(pq,keep,w);

}


static Snapshot makeSnapshot(double ratio){

    Snapshot s;
 s.ratio=ratio;
 s.nv=activeVertices;

    tmpId.assign(N,-1);
 s.V.reserve(activeVertices);

    for(int i=0;
i<N;
i++) if(alive[i]){
 tmpId[i]=(int)s.V.size();
 s.V.push_back(P[i]);
 }

    s.F.reserve(activeFaces);

    for(int i=0;
i<M;
i++) if(F[i].active){

        int a=F[i].v[0],b=F[i].v[1],c=F[i].v[2];

        if(a!=b && b!=c && c!=a && tmpId[a]>=0 && tmpId[b]>=0 && tmpId[c]>=0) s.F.push_back({
tmpId[a],tmpId[b],tmpId[c]}
);

    }

    s.nv=(int)s.V.size();

    return s;

}



static int evalR=512, evalPix=512*512;

static constexpr int VIEWS=6;


static void renderMesh(const vector<Vec3>& verts,const vector<Tri>& tris,RenderMaps& maps,int R){

    int PIX=R*R;
 double focal=800.0*(double)R/1024.0, center=0.5*(double)R;

    maps.R=R;
 maps.depth.assign((size_t)VIEWS*PIX,255.0f);
 maps.norm.assign((size_t)VIEWS*PIX*3,127.5f);

    int nv=(int)verts.size();
 int nt=(int)tris.size();

    vector<float> U(nv), VV(nv), Z(nv);

    vector<Vec3> fn(nt);

    for(int i=0;
i<nt;
i++){

        const Tri& t=tris[i];
 Vec3 cr=crossv(verts[t.b]-verts[t.a],verts[t.c]-verts[t.a]);
 double l=normv(cr);

        fn[i]=(l>1e-300)?cr/l:Vec3(0,0,0);

    }

    for(int view=0;
view<VIEWS;
view++){

        int ax=view/2;
 int sg=(view&1)?-1:1;

        for(int i=0;
i<nv;
i++){

            const Vec3& p=verts[i];
 double sx,sy,dep;

            if(ax==0){
 sx=p.y;
 sy=p.z;
 dep=2.5 - sg*p.x;
 if(sg<0) sx=-sx;
 }

            else if(ax==1){
 sx=p.x;
 sy=p.z;
 dep=2.5 - sg*p.y;
 if(sg>0) sx=-sx;
 }

            else {
 sx=p.x;
 sy=p.y;
 dep=2.5 - sg*p.z;
 if(sg<0) sx=-sx;
 }

            Z[i]=(float)dep;
 U[i]=(float)(focal*sx/dep + center);
 VV[i]=(float)(focal*sy/dep + center);

        }

        float* zbuf = maps.depth.data() + (size_t)view*PIX;

        float* nbuf = maps.norm.data() + (size_t)view*PIX*3;

        for(int ti=0;
ti<nt;
ti++){

            const Tri& t=tris[ti];
 int ia=t.a,ib=t.b,ic=t.c;

            float u0=U[ia],u1=U[ib],u2=U[ic];
 float v0=VV[ia],v1=VV[ib],v2=VV[ic];

            float z0=Z[ia],z1=Z[ib],z2=Z[ic];

            if(!(z0>0 && z1>0 && z2>0)) continue;

            float minx=min(u0,min(u1,u2)), maxx=max(u0,max(u1,u2));

            float miny=min(v0,min(v1,v2)), maxy=max(v0,max(v1,v2));

            int x0i=max(0,(int)floor(minx));
 int x1i=min(R-1,(int)ceil(maxx));

            int y0i=max(0,(int)floor(miny));
 int y1i=min(R-1,(int)ceil(maxy));

            if(x0i>x1i || y0i>y1i) continue;

            float den=(v1-v2)*(u0-u2)+(u2-u1)*(v0-v2);

            if(fabs(den)<1e-20f) continue;

            float invDen=1.0f/den;

            Vec3 n=fn[ti];

            float nr=(float)((n.x+1.0)*127.5), ng=(float)((n.y+1.0)*127.5), nb=(float)((n.z+1.0)*127.5);

            for(int y=y0i;
y<=y1i;
y++){

                float py=(float)y+0.5f;

                int row=y*R;

                for(int x=x0i;
x<=x1i;
x++){

                    float px=(float)x+0.5f;

                    float w0=((v1-v2)*(px-u2)+(u2-u1)*(py-v2))*invDen;

                    float w1=((v2-v0)*(px-u2)+(u0-u2)*(py-v2))*invDen;

                    float w2=1.0f-w0-w1;

                    if(w0>=-1e-6f && w1>=-1e-6f && w2>=-1e-6f){

                        float dep=1.0f/(w0/z0+w1/z1+w2/z2);

                        int idx=row+x;

                        if(dep<zbuf[idx]){

                            zbuf[idx]=dep;

                            float* q=nbuf+idx*3;
 q[0]=nr;
 q[1]=ng;
 q[2]=nb;

                        }

                    }

                }

            }

        }

    }

}


static inline double rectSum(const vector<double>& I,int stride,int x0,int y0,int x1,int y1){

    return I[(size_t)y1*stride+x1]-I[(size_t)y0*stride+x1]-I[(size_t)y1*stride+x0]+I[(size_t)y0*stride+x0];

}

static double ssimChannel(const float* X,int strideX,const float* Y,int strideY,const float* dX,const float* dY,int R,
                          vector<double>& IX,vector<double>& IY,vector<double>& IX2,vector<double>& IY2,vector<double>& IXY){

    int W=R+1;
 size_t SZ=(size_t)W*W;

    fill(IX.begin(),IX.begin()+SZ,0.0);
 fill(IY.begin(),IY.begin()+SZ,0.0);
 fill(IX2.begin(),IX2.begin()+SZ,0.0);
 fill(IY2.begin(),IY2.begin()+SZ,0.0);
 fill(IXY.begin(),IXY.begin()+SZ,0.0);

    for(int y=1;
y<=R;
y++){

        double sx=0,sy=0,sx2=0,sy2=0,sxy=0;
 int row=(y-1)*R;

        for(int x=1;
x<=R;
x++){

            int p=row+(x-1);

            double xv=X[(size_t)p*strideX], yv=Y[(size_t)p*strideY];

            sx+=xv;
 sy+=yv;
 sx2+=xv*xv;
 sy2+=yv*yv;
 sxy+=xv*yv;

            size_t id=(size_t)y*W+x, up=(size_t)(y-1)*W+x;

            IX[id]=IX[up]+sx;
 IY[id]=IY[up]+sy;
 IX2[id]=IX2[up]+sx2;
 IY2[id]=IY2[up]+sy2;
 IXY[id]=IXY[up]+sxy;

        }

    }

    const int rad=5, area=121;
 const double c1=6.5025, c2=58.5225;

    long long cnt=0;
 long double acc=0.0L;

    for(int y=rad;
y<R-rad;
y++){

        int row=y*R;

        for(int x=rad;
x<R-rad;
x++){

            int p=row+x;

            if(!(dX[p]<254.0f || dY[p]<254.0f)) continue;

            int x0=x-rad, y0=y-rad, x1=x+rad+1, y1=y+rad+1;

            double sx=rectSum(IX,W,x0,y0,x1,y1), sy=rectSum(IY,W,x0,y0,x1,y1);

            double sx2=rectSum(IX2,W,x0,y0,x1,y1), sy2=rectSum(IY2,W,x0,y0,x1,y1), sxy=rectSum(IXY,W,x0,y0,x1,y1);

            double mux=sx/area, muy=sy/area;

            double vx=sx2/area-mux*mux, vy=sy2/area-muy*muy, cov=sxy/area-mux*muy;

            if(vx<0 && vx>-1e-6) vx=0;
 if(vy<0 && vy>-1e-6) vy=0;

            double denom=(mux*mux+muy*muy+c1)*(vx+vy+c2);

            double val=(denom!=0.0)?((2*mux*muy+c1)*(2*cov+c2)/denom):1.0;

            acc += val;
 ++cnt;

        }

    }

    return cnt? (double)(acc/cnt) : 1.0;

}

static double renderSSIM(const RenderMaps& orig,const RenderMaps& cand){

    int R=orig.R, PIX=R*R, W=R+1;
 size_t SZ=(size_t)W*W;

    vector<double> IX(SZ),IY(SZ),IX2(SZ),IY2(SZ),IXY(SZ);

    double total=0.0;

    for(int view=0;
view<VIEWS;
view++){

        const float* od=orig.depth.data()+(size_t)view*PIX;

        const float* cd=cand.depth.data()+(size_t)view*PIX;

        double ns=0.0;

        for(int ch=0;
ch<3;
ch++){

            const float* on=orig.norm.data()+(size_t)view*PIX*3+ch;

            const float* cn=cand.norm.data()+(size_t)view*PIX*3+ch;

            ns += ssimChannel(on,3,cn,3,od,cd,R,IX,IY,IX2,IY2,IXY);

        }

        ns/=3.0;

        double ds = ssimChannel(od,1,cd,1,od,cd,R,IX,IY,IX2,IY2,IXY);

        total += 0.5*(ns+ds);

    }

    return total/6.0;

}


static void simplifyCore(){

    if(N<=4 || M==0){

        finalSnap.V=P;

        finalSnap.F.clear();

        for(auto& f:F) finalSnap.F.push_back({
f.v[0],f.v[1],f.v[2]}
);

        haveFinalSnap=true;

        return;

    }

    double xmin=P[0].x,xmax=P[0].x,ymin=P[0].y,ymax=P[0].y,zmin=P[0].z,zmax=P[0].z;

    for(int i=1;
i<N;
i++){
 xmin=min(xmin,P[i].x);
 xmax=max(xmax,P[i].x);
 ymin=min(ymin,P[i].y);
 ymax=max(ymax,P[i].y);
 zmin=min(zmin,P[i].z);
 zmax=max(zmax,P[i].z);
 }

    double diag=sqrt((xmax-xmin)*(xmax-xmin)+(ymax-ymin)*(ymax-ymin)+(zmax-zmin)*(zmax-zmin));

    hausTau=0.05*diag*0.992;
 hausTau2=hausTau*hausTau;

    if(N<8000) visualPixelTol=28.0;
 else if(N<60000) visualPixelTol=23.0;
 else if(N<200000) visualPixelTol=19.0;
 else visualPixelTol=16.0;

    visualPixelTol2=visualPixelTol*visualPixelTol;


    inc.assign(N,{
}
);
 for(int i=0;
i<N;
i++) inc[i].reserve(8);

    activeFaces=M;
 activeVertices=N;

    for(int i=0;
i<M;
i++){

        F[i].active=1;
 F[i].n=computeFaceNormal(F[i]);

        inc[F[i].v[0]].push_back(i);
 inc[F[i].v[1]].push_back(i);
 inc[F[i].v[2]].push_back(i);

    }

    Qv.assign(N,Quadric());

    for(int i=0;
i<M;
i++){

        Face& f=F[i];
 Vec3 p0=P[f.v[0]],p1=P[f.v[1]],p2=P[f.v[2]];

        Vec3 cr=crossv(p1-p0,p2-p0);
 double area2=normv(cr);
 if(area2<=1e-300) continue;

        Vec3 n=cr/area2;
 double d=-dotv(n,p0);

        
        double w=max(1e-8,0.5*area2) + 2e-7;

        Quadric q=Quadric::plane(n.x,n.y,n.z,d,w);

        Qv[f.v[0]].add(q);
 Qv[f.v[1]].add(q);
 Qv[f.v[2]].add(q);

    }

    bbMin.resize(N);
 bbMax.resize(N);

    for(int i=0;
i<N;
i++){
 bbMin[i]={
(float)P[i].x,(float)P[i].y,(float)P[i].z}
;
 bbMax[i]=bbMin[i];
 }

    alive.assign(N,1);
 ver.assign(N,0);
 markA.assign(N,0);
 markB.assign(N,0);


    vector<uint64_t> edges;
 edges.reserve((size_t)M*3);

    auto addEdge=[&](int a,int b){
 if(a>b) swap(a,b);
 edges.push_back(((uint64_t)(uint32_t)a<<32)|(uint32_t)b);
 }
;

    for(int i=0;
i<M;
i++){
 addEdge(F[i].v[0],F[i].v[1]);
 addEdge(F[i].v[1],F[i].v[2]);
 addEdge(F[i].v[2],F[i].v[0]);
 }

    sort(edges.begin(),edges.end());
 edges.erase(unique(edges.begin(),edges.end()),edges.end());

    vector<Candidate> heapVec;
 heapVec.reserve(edges.size());

    for(uint64_t key: edges){
 int u=(int)(key>>32), v=(int)(key&0xffffffffu);
 EvalCand ec=computeCandidate(u,v);
 if(ec.ok) heapVec.push_back({
ec.score,u,v,0,0}
);
 }

    edges.clear();
 edges.shrink_to_fit();

    priority_queue<Candidate> pq(less<Candidate>(),move(heapVec));


    vector<double> ratios={
0.50,0.35,0.25,0.22,0.20,0.19,0.18,0.17,0.16,0.15,0.14,0.13,0.12,0.115,0.110,0.105,0.100,0.095,0.090,0.085,0.080,0.075,0.070,0.065,0.060,0.055}
;

    vector<int> cps;

    int last=-1;

    for(double r: ratios){
 int c=max(4,(int)ceil(N*r));
 if(c<N && c!=last){
 cps.push_back(c);
 last=c;
 }
 }

    sort(cps.begin(),cps.end(),greater<int>());

    int cpIdx=0;
 int minTarget=cps.empty()?max(4,N/10):cps.back();


    while(!pq.empty() && activeVertices>minTarget){

        Candidate c=pq.top();
 pq.pop();
 int u=c.u,v=c.v;

        if(u==v || !alive[u] || !alive[v]) continue;
 if(u>v) swap(u,v);

        if(ver[u]!=c.vu || ver[v]!=c.vv) continue;

        EvalCand ec=computeCandidate(u,v);
 if(!ec.ok) continue;

        if(!linkConditionOK(u,v)) continue;

        int keep=u, rem=v;
 if(inc[v].size()>inc[u].size()){
 keep=v;
 rem=u;
 }

        
        if(!flipOK(keep,rem,ec.pos,-0.02)) continue;

        doCollapse(keep,rem,ec.pos,pq);

        while(cpIdx<(int)cps.size() && activeVertices<=cps[cpIdx]){

            snapshots.push_back(makeSnapshot((double)activeVertices/(double)N));

            ++cpIdx;

        }

    }

    
    if(snapshots.empty() || snapshots.back().nv!=activeVertices) snapshots.push_back(makeSnapshot((double)activeVertices/(double)N));

}


static void chooseByInternalSSIM(const RenderMaps& origMaps){

    if(snapshots.empty()){
 finalSnap=makeSnapshot((double)activeVertices/(double)N);
 haveFinalSnap=true;
 return;
 }

    int R=origMaps.R;

    double guard = (R>=1000 ? 0.904 : (R>=700 ? 0.910 : 0.916));

    RenderMaps cand;

    int bestIdx=-1;
 double bestScore=-1.0;

    auto evalOne = [&](int idx)->double{

        if(idx<0 || idx>=(int)snapshots.size()) return -1.0;

        if(snapshots[idx].ssim >= -0.5) return snapshots[idx].ssim;

        renderMesh(snapshots[idx].V,snapshots[idx].F,cand,R);

        double sc=renderSSIM(origMaps,cand);

        snapshots[idx].ssim=sc;

        if(sc>bestScore){
 bestScore=sc;
 bestIdx=idx;
 }

        return sc;

    }
;


    
    
    vector<double> desired={
0.055,0.065,0.075,0.085,0.095,0.105,0.115,0.130,0.160,0.200,0.250,0.350,0.500}
;

    vector<int> coarse;

    vector<unsigned char> used(snapshots.size(),0);

    for(double d: desired){

        int bi=-1;
 double bd=1e100;

        for(int i=0;
i<(int)snapshots.size();
i++){

            double diff=fabs(snapshots[i].ratio-d);

            if(diff<bd){
 bd=diff;
 bi=i;
 }

        }

        if(bi>=0 && !used[bi]){
 used[bi]=1;
 coarse.push_back(bi);
 }

    }

    sort(coarse.begin(),coarse.end(),[&](int a,int b){
 return snapshots[a].ratio < snapshots[b].ratio;
 }
);


    int prevFail=-1;

    for(int idx: coarse){

        double sc=evalOne(idx);

        if(sc>=guard){

            int chosen=idx;

            if(prevFail!=-1){

                
                
                for(int j=prevFail-1;
j>=idx;
j--){

                    double sj=evalOne(j);

                    if(sj>=guard){
 chosen=j;
 break;
 }

                }

            }

            finalSnap=std::move(snapshots[chosen]);
 haveFinalSnap=true;
 return;

        }

        prevFail=idx;

    }

    if(bestIdx<0){

        for(int i=0;
i<(int)snapshots.size();
i++) evalOne(i);

    }

    if(bestIdx>=0){
 finalSnap=std::move(snapshots[bestIdx]);
 haveFinalSnap=true;
 return;
 }

    finalSnap=makeSnapshot((double)activeVertices/(double)N);
 haveFinalSnap=true;

}



static bool build(const vector<::Vec3>&SV,const vector<::Face>&SF,vector<::Vec3>&OV,vector<::Face>&OF){

    N=(int)SV.size();
M=(int)SF.size();
if(N<=4||M==0)return false;

    P.assign(N,Vec3());
for(int i=0;
i<N;
i++)P[i]=Vec3(SV[i].x,SV[i].y,SV[i].z);

    F.assign(M,Face());
for(int i=0;
i<M;
i++){
F[i].v[0]=SF[i].v[0];
F[i].v[1]=SF[i].v[1];
F[i].v[2]=SF[i].v[2];
F[i].active=1;
}

    snapshots.clear();
finalSnap=Snapshot();
haveFinalSnap=false;
stampA=1;
stampB=1;
tmpId.clear();

    evalR=1024;
evalPix=evalR*evalR;
vector<Tri>origTris;
origTris.reserve(M);
for(int i=0;
i<M;
i++)origTris.push_back({
F[i].v[0],F[i].v[1],F[i].v[2]}
);

    RenderMaps origMaps;
renderMesh(P,origTris,origMaps,evalR);
origTris.clear();
origTris.shrink_to_fit();

    simplifyCore();
inc.clear();
inc.shrink_to_fit();
Qv.clear();
Qv.shrink_to_fit();
bbMin.clear();
bbMin.shrink_to_fit();
bbMax.clear();
bbMax.shrink_to_fit();
alive.clear();
alive.shrink_to_fit();
ver.clear();
ver.shrink_to_fit();
markA.clear();
markA.shrink_to_fit();
markB.clear();
markB.shrink_to_fit();
tmpId.clear();
tmpId.shrink_to_fit();

    chooseByInternalSSIM(origMaps);
if(!haveFinalSnap||finalSnap.V.empty()||finalSnap.F.empty()||finalSnap.V.size()>=SV.size())return false;

    OV.clear();
OF.clear();
OV.reserve(finalSnap.V.size());
OF.reserve(finalSnap.F.size());
for(auto&p:finalSnap.V)OV.push_back(::Vec3{
p.x,p.y,p.z}
);
for(auto&t:finalSnap.F){
::Face f;
f.v[0]=t.a;
f.v[1]=t.b;
f.v[2]=t.c;
OF.push_back(f);
}
return true;

}
}

static void GN(){
HJ=chrono::steady_clock::now();
const double d=CL;
const double AJ=max(1e-24,1e-18*d*d);
const double W=max(1e-11,1e-9*d);
const double CM=1.0-1e-10;
const bool EG=EZ();
const bool medium_detail_profile=AS&&N>=8000&&N<50000&&BL>=0.900&&AL>=0.990&&Z>=0.003&&Z<=0.050&&AH<=0.010&&BG<=0.005;
const bool EA=AS&&N>=50000&&AL>=0.900&&AH<=0.050&&Z>=0.020&&Z<=0.160&&BG<=0.005;
const bool BJ=AS&&N>=400&&N<8000&&AL>=0.995&&AH<=0.006&&BG<=0.005&&(BL<=0.760||Z<=0.006||Z>=0.050);
const bool small_detail_smooth_profile=AS&&N>=1000&&N<8000&&!BJ&&BL>=0.800&&AL>=0.980&&Z>=0.005&&Z<=0.060&&AH<=0.012&&BG<=0.005;
const bool small_proxy_chase_profile=BJ||small_detail_smooth_profile;
const bool DP=medium_detail_profile||EA;
const bool BX=small_detail_smooth_profile&&N>=1500&&N<8000&&Z>=0.010&&AH<=0.010;
const bool tiny_bumpy_chase_profile=AS&&N>=900&&N<1500&&!BJ&&BG<=0.006&&AL>=0.820&&Z>=0.010&&Z<=0.160&&AH<=0.070;
const bool ES=AS&&N>=60000&&N<=260000&&BG<=0.006&&BL<=0.820&&AL>=0.860&&AL<=0.980&&Z>=0.070&&Z<=0.180&&AH>=0.015&&AH<=0.070;
const bool attempt_aggressive_profile=(N>=200000)||EG||DP||BX;
const bool FY=EG&&N>=3000&&N<=50000;
AP CG;
if(attempt_aggressive_profile||small_proxy_chase_profile||tiny_bumpy_chase_profile||ES)CG=AD();
const bool small_visual_restart_profile=N>=300&&N<3000;
AP small_visual_restart_state;
if(small_visual_restart_profile)small_visual_restart_state=AD();
double EE=0.020*d;
double FJ=0.0058*d;
double DR=4.4;
if(N>=200000){
EE=0.027*d;
FJ=0.0082*d;
DR=6.1;
}
else if(N>=30000){
EE=0.024*d;
FJ=0.0072*d;
DR=5.3;
}
else if(N<=1000){
EE=0.014*d;
FJ=0.0036*d;
DR=3.0;
}
vector<AE>passes;
auto add_pass=[&](double AK,double plane,double AU,bool AA=false){
AE p;
p.AI=AK;
p.BB=plane;
p.BQ=cos(AU*acos(-1.0)/180.0);
p.W=W;
p.AT=CM;
p.AJ=AJ;
p.AA=AA;
passes.push_back(p);
}
;
add_pass(0.004*d,0.0010*d,1.0);
add_pass(0.008*d,0.0020*d,1.8);
add_pass(0.012*d,0.0032*d,2.8);
add_pass(EE,FJ,DR);
add_pass(EE,FJ,DR);
constexpr double WORK_SECONDS=16.5;
auto AX=[&](const vector<AE>&HZ,bool allow_early_break)->bool{
for(const AE&params:HZ){
int HY=0;
for(int fid=0;
fid<(int)faces.size();
++fid){
if((fid&4095)==0&&elapsed_seconds()>WORK_SECONDS)return false;
if(!BR[fid])continue;
Face f=faces[fid];
if(!BU[f.v[0]]||!BU[f.v[1]]||!BU[f.v[2]])continue;
struct EdgeCandidate{
double JQ;
int a;
int b;
}
;
EdgeCandidate e[3]={
{
norm2(P[f.v[0]]-P[f.v[1]]),f.v[0],f.v[1]}
,{
norm2(P[f.v[1]]-P[f.v[2]]),f.v[1],f.v[2]}
,{
norm2(P[f.v[2]]-P[f.v[0]]),f.v[2],f.v[0]}
,}
;
sort(e,e+3,[](const EdgeCandidate&lhs,const EdgeCandidate&rhs){
return lhs.JQ<rhs.JQ;
}
);
for(const EdgeCandidate&cand:e){
if(GD(cand.a,cand.b,params)){
++HY;
break;
}
}
}
if(allow_early_break&&HY==0&&params.AI>=EE*0.999)break;
}
return true;
}
;
if(!AX(passes,true))return;
if(small_proxy_chase_profile&&elapsed_seconds()<12.4){
AP AM=AD();
const int safe_vertices=count_output_vertices_estimate();
restore_state(CG);
vector<AE>chase_passes;
auto DD=[&](double AK,double plane,double AU,bool AA=false){
AE p;
p.AI=AK;
p.BB=plane;
p.BQ=cos(AU*acos(-1.0)/180.0);
p.W=W;
p.AT=CM;
p.AJ=AJ;
p.AA=AA;
chase_passes.push_back(p);
}
;
double CX=0.026*d;
double DJ=0.0078*d;
double DK=5.8;
double AO=92.0;
double proxy_threshold=0.970;
int chase_proxy_res=512;
if(BJ){
if(N<=1000&&Z>=0.050){
CX=0.050*d;
DJ=0.0150*d;
DK=10.0;
AO=96.0;
proxy_threshold=0.970;
}
else if(N<3000&&Z<=0.006){
CX=0.034*d;
DJ=0.0102*d;
DK=7.7;
AO=92.0;
proxy_threshold=0.955;
}
else{
CX=0.030*d;
DJ=0.0090*d;
DK=6.8;
AO=92.0;
proxy_threshold=0.965;
}
}
else if(small_detail_smooth_profile){
if(N<3000){
CX=0.034*d;
DJ=0.0102*d;
DK=7.7;
AO=96.0;
proxy_threshold=0.952;
chase_proxy_res=(N<1500?768:512);
}
else{
CX=0.024*d;
DJ=0.0072*d;
DK=5.3;
AO=92.0;
proxy_threshold=0.965;
}
}
DD(0.004*d,0.0010*d,1.0);
DD(0.008*d,0.0020*d,1.8);
DD(0.012*d,0.0032*d,2.8);
DD(CX,DJ,DK);
DD(CX,DJ,DK);
if(BJ&&N>=1500){
DD(CX,DJ,DK);
}
bool keep_chase=false;
if(AX(chase_passes,true)&&elapsed_seconds()<16.2){
const int chase_vertices=count_output_vertices_estimate();
if(chase_vertices>0&&safe_vertices>0&&(double)chase_vertices*100.0<=(double)safe_vertices*AO){
const double proxy=visual_proxy_score(chase_proxy_res);
keep_chase=proxy>=proxy_threshold;
}
}
if(!keep_chase)restore_state(AM);
}
if(N>=300&&N<=6000&&elapsed_seconds()<10.8){
struct SmallVisualTrial{
double AK;
double plane;
double AU;
double proxy_threshold;
double CB;
int GG;
bool AA;
}
;
vector<SmallVisualTrial>BT;
const bool CV=AS&&AL>=0.960&&AH<=0.035&&BG<=0.010;
const bool CQ=AS&&AL>=0.990&&AH<=0.012&&BG<=0.006;
if(N<1000){
BT.push_back({
0.026*d,0.0078*d,5.8,0.975,99.0,1,false}
);
BT.push_back({
0.038*d,0.0114*d,8.4,CV?0.950:0.965,99.0,2,CV}
);
BT.push_back({
0.050*d,0.0150*d,10.8,CQ?0.935:0.955,98.0,2,CQ}
);
}
else if(N<3000){
BT.push_back({
0.028*d,0.0084*d,6.2,0.970,99.0,1,false}
);
BT.push_back({
0.038*d,0.0114*d,8.4,CV?0.945:0.960,98.8,2,CV}
);
BT.push_back({
0.049*d,0.0147*d,10.8,CQ?0.930:0.950,98.0,2,CQ}
);
}
else{
BT.push_back({
0.030*d,0.0090*d,6.8,0.965,99.0,1,false}
);
BT.push_back({
0.041*d,0.0123*d,9.2,CV?0.940:0.955,98.5,1,CV}
);
BT.push_back({
0.050*d,0.0150*d,11.0,CQ?0.925:0.945,97.8,1,CQ}
);
}
for(const SmallVisualTrial&AZ:BT){
if(elapsed_seconds()>14.2)break;
AP trial_safe_state=AD();
const int before_vertices=count_output_vertices_estimate();
if(before_vertices<=0)break;
vector<AE>trial_passes;
AE p;
p.AI=AZ.AK;
p.BB=AZ.plane;
p.BQ=cos(AZ.AU*acos(-1.0)/180.0);
p.W=W;
p.AT=CM;
p.AJ=AJ;
p.AA=AZ.AA;
for(int i=0;
i<AZ.GG;
++i)trial_passes.push_back(p);
bool keep_visual_trial=false;
if(AX(trial_passes,true)&&elapsed_seconds()<16.2){
const int after_vertices=count_output_vertices_estimate();
if(after_vertices>0&&(double)after_vertices*100.0<=(double)before_vertices*AZ.CB){
const double proxy=visual_proxy_score(512);
keep_visual_trial=proxy>=AZ.proxy_threshold;
}
}
if(!keep_visual_trial){
restore_state(trial_safe_state);
if(AZ.proxy_threshold<=0.940)break;
}
}
}
if(BX&&elapsed_seconds()<12.2){
AP AM=AD();
const int safe_vertices=count_output_vertices_estimate();
struct DetailChaseTrial{
double small_radius;
double small_plane;
double small_angle;
double big_radius;
double big_plane;
double big_angle;
int GG;
int AG;
double proxy_threshold;
double CB;
}
;
vector<DetailChaseTrial>detail_trials;
detail_trials.push_back({
0.030*d,0.0090*d,6.8,0.032*d,0.0096*d,7.2,3,1024,0.900,98.0}
);
detail_trials.push_back({
0.024*d,0.0072*d,5.8,0.026*d,0.0078*d,6.2,2,512,0.905,92.0}
);
bool accepted_detail_chase=false;
for(const DetailChaseTrial&AZ:detail_trials){
if(elapsed_seconds()>13.4)break;
restore_state(CG);
vector<AE>detail_chase_passes;
auto add_detail_chase_pass=[&](double AK,double plane,double AU){
AE p;
p.AI=AK;
p.BB=plane;
p.BQ=cos(AU*acos(-1.0)/180.0);
p.W=W;
p.AT=CM;
p.AJ=AJ;
detail_chase_passes.push_back(p);
}
;
double detail_chase_radius=AZ.small_radius;
double detail_chase_plane=AZ.small_plane;
double detail_chase_angle=AZ.small_angle;
if(N>=3000){
detail_chase_radius=AZ.big_radius;
detail_chase_plane=AZ.big_plane;
detail_chase_angle=AZ.big_angle;
}
add_detail_chase_pass(0.004*d,0.0010*d,1.0);
add_detail_chase_pass(0.008*d,0.0020*d,1.8);
add_detail_chase_pass(0.012*d,0.0032*d,2.8);
for(int i=0;
i<AZ.GG;
++i){
add_detail_chase_pass(detail_chase_radius,detail_chase_plane,detail_chase_angle);
}
bool keep_detail_chase=false;
if(safe_vertices>0&&AX(detail_chase_passes,true)&&elapsed_seconds()<16.1){
const int BO=count_output_vertices_estimate();
if(BO>0&&(double)BO*100.0<=(double)safe_vertices*AZ.CB){
keep_detail_chase=visual_proxy_score(AZ.AG)>=AZ.proxy_threshold;
}
}
if(keep_detail_chase){
accepted_detail_chase=true;
break;
}
restore_state(AM);
}
if(!accepted_detail_chase)restore_state(AM);
}
if(small_visual_restart_profile&&elapsed_seconds()<12.0){
AP CF=AD();
const int AN=count_output_vertices_estimate();
restore_state(small_visual_restart_state);
vector<AE>restart_passes;
auto DE=[&](double AK,double plane,double AU,bool AA=false){
AE p;
p.AI=AK;
p.BB=plane;
p.BQ=cos(AU*acos(-1.0)/180.0);
p.W=W;
p.AT=CM;
p.AJ=AJ;
p.AA=AA;
restart_passes.push_back(p);
}
;
double restart_radius=0.033*d;
double restart_plane=0.0099*d;
double restart_angle=7.5;
double restart_proxy_threshold=0.960;
double restart_required_percent=98.0;
if(N<700){
restart_radius=0.033*d;
restart_plane=0.0099*d;
restart_angle=7.5;
restart_proxy_threshold=0.965;
restart_required_percent=98.0;
}
else if(N<1500){
restart_radius=0.034*d;
restart_plane=0.0102*d;
restart_angle=7.8;
restart_proxy_threshold=0.942;
restart_required_percent=97.0;
}
else{
restart_radius=0.036*d;
restart_plane=0.0108*d;
restart_angle=8.2;
restart_proxy_threshold=0.940;
restart_required_percent=97.0;
}
DE(0.004*d,0.0010*d,1.0);
DE(0.008*d,0.0020*d,1.8);
DE(0.012*d,0.0032*d,2.8);
DE(restart_radius,restart_plane,restart_angle);
DE(restart_radius,restart_plane,restart_angle);
bool keep_restart=false;
if(AN>0&&AX(restart_passes,true)&&elapsed_seconds()<16.2){
const int restart_vertices=count_output_vertices_estimate();
if(restart_vertices>0&&restart_vertices<AN&&(double)restart_vertices*100.0<=(double)AN*restart_required_percent){
const double proxy=visual_proxy_score(512);
keep_restart=proxy>=restart_proxy_threshold;
}
}
if(!keep_restart)restore_state(CF);
}
if(tiny_bumpy_chase_profile&&elapsed_seconds()<12.3){
AP CF=AD();
const int AN=count_output_vertices_estimate();
restore_state(CG);
vector<AE>tiny_passes;
auto DT=[&](double AK,double plane,double AU){
AE p;
p.AI=AK;
p.BB=plane;
p.BQ=cos(AU*acos(-1.0)/180.0);
p.W=W;
p.AT=CM;
p.AJ=AJ;
tiny_passes.push_back(p);
}
;
DT(0.004*d,0.0010*d,1.0);
DT(0.008*d,0.0020*d,1.8);
DT(0.012*d,0.0032*d,2.8);
DT(0.032*d,0.0096*d,7.8);
DT(0.032*d,0.0096*d,7.8);
DT(0.032*d,0.0096*d,7.8);
bool keep_tiny=false;
if(AN>0&&AX(tiny_passes,true)&&elapsed_seconds()<16.2){
const int tiny_vertices=count_output_vertices_estimate();
if(tiny_vertices>0&&tiny_vertices<AN&&(double)tiny_vertices*100.0<=(double)AN*96.0){
keep_tiny=visual_proxy_score(768)>=0.935;
}
}
if(!keep_tiny)restore_state(CF);
}
if(attempt_aggressive_profile&&elapsed_seconds()<11.5){
AP AM=AD();
const int safe_vertices=count_output_vertices_estimate();
restore_state(CG);
vector<AE>GH;
auto CO=[&](double AK,double plane,double AU,bool AA=false){
AE p;
p.AI=AK;
p.BB=plane;
p.BQ=cos(AU*acos(-1.0)/180.0);
p.W=W;
p.AT=CM;
p.AJ=AJ;
p.AA=AA;
GH.push_back(p);
}
;
double aggressive_radius=0.0265*d;
double aggressive_plane=0.0079*d;
double aggressive_angle=5.9;
if(N>=200000){
aggressive_radius=0.033*d;
aggressive_plane=0.0100*d;
aggressive_angle=7.4;
}
else if(N<10000){
aggressive_radius=0.024*d;
aggressive_plane=0.0072*d;
aggressive_angle=5.3;
}
else if(N<30000){
aggressive_radius=0.0265*d;
aggressive_plane=0.0079*d;
aggressive_angle=5.9;
}
CO(0.004*d,0.0010*d,1.0);
CO(0.008*d,0.0020*d,1.8);
CO(0.012*d,0.0032*d,2.8);
CO(aggressive_radius,aggressive_plane,aggressive_angle);
CO(aggressive_radius,aggressive_plane,aggressive_angle);
bool FO=false;
if(AX(GH,true)&&elapsed_seconds()<14.0){
if(BX&&!EG&&!DP){
const int BO=count_output_vertices_estimate();
if(BO>0&&safe_vertices>0&&(double)BO*100.0<=(double)safe_vertices*92.0&&elapsed_seconds()<16.2){
const double proxy=visual_proxy_score(512);
FO=proxy>=0.915;
}
}
else{
const double proxy=visual_proxy_score(64);
FO=proxy>=0.955;
}
}
if(!FO)restore_state(AM);
}
if(N>=200000&&elapsed_seconds()<12.5){
AP CF=AD();
const int AN=count_output_vertices_estimate();
restore_state(CG);
vector<AE>IQ;
auto EB=[&](double AK,double plane,double AU){
AE p;
p.AI=AK;
p.BB=plane;
p.BQ=cos(AU*acos(-1.0)/180.0);
p.W=W;
p.AT=CM;
p.AJ=AJ;
IQ.push_back(p);
}
;
EB(0.004*d,0.0010*d,1.0);
EB(0.008*d,0.0020*d,1.8);
EB(0.012*d,0.0032*d,2.8);
EB(0.045*d,0.0135*d,10.2);
EB(0.045*d,0.0135*d,10.2);
bool JB=false;
if(AX(IQ,true)&&elapsed_seconds()<15.8){
const int ultra_vertices=count_output_vertices_estimate();
if(ultra_vertices>0&&ultra_vertices*100<=AN*96){
const double proxy=visual_proxy_score(64);
JB=proxy>=0.955;
}
}
if(!JB)restore_state(CF);
}
if(DP&&elapsed_seconds()<12.5){
AP CF=AD();
const int AN=count_output_vertices_estimate();
restore_state(CG);
vector<AE>IH;
auto EN=[&](double AK,double plane,double AU){
AE p;
p.AI=AK;
p.BB=plane;
p.BQ=cos(AU*acos(-1.0)/180.0);
p.W=W;
p.AT=CM;
p.AJ=AJ;
IH.push_back(p);
}
;
double GO=0.0310*d;
double HB=0.00930*d;
double HC=6.95;
double detail_threshold=0.955;
int detail_proxy_res=128;
int final_repeats=3;
if(EA){
if(N>=200000){
GO=0.050*d;
HB=0.0150*d;
HC=11.0;
}
else{
GO=0.047*d;
HB=0.0141*d;
HC=10.5;
}
detail_threshold=0.940;
final_repeats=4;
}
EN(0.004*d,0.0010*d,1.0);
EN(0.008*d,0.0020*d,1.8);
EN(0.012*d,0.0032*d,2.8);
for(int i=0;
i<final_repeats;
++i){
EN(GO,HB,HC);
}
bool IS=false;
if(AX(IH,true)&&elapsed_seconds()<16.2){
const int detail_vertices=count_output_vertices_estimate();
if(detail_vertices>0&&AN>0&&detail_vertices*100<=AN*96){
const double proxy=visual_proxy_score(detail_proxy_res);
IS=proxy>=detail_threshold;
}
}
if(!IS)restore_state(CF);
}
if(EA&&N>=30000&&N<100000&&elapsed_seconds()<12.8){
AP CF=AD();
const int AN=count_output_vertices_estimate();
restore_state(CG);
vector<AE>IR;
auto EY=[&](double AK,double plane,double AU){
AE p;
p.AI=AK;
p.BB=plane;
p.BQ=cos(AU*acos(-1.0)/180.0);
p.W=W;
p.AT=CM;
p.AJ=AJ;
IR.push_back(p);
}
;
EY(0.004*d,0.0010*d,1.0);
EY(0.008*d,0.0020*d,1.8);
EY(0.012*d,0.0032*d,2.8);
for(int i=0;
i<4;
++i){
EY(0.055*d,0.0165*d,12.2);
}
bool JA=false;
if(AX(IR,true)&&elapsed_seconds()<16.3){
const int boost_vertices=count_output_vertices_estimate();
if(boost_vertices>0&&AN>0&&boost_vertices*100<=AN*95){
const double proxy=visual_proxy_score(128);
JA=proxy>=0.960;
}
}
if(!JA)restore_state(CF);
}
if(ES&&elapsed_seconds()<13.2){
AP boost_safe_state=AD();
const int boost_start_vertices=count_output_vertices_estimate();
vector<AE>GW;
auto FZ=[&](double AK,double plane,double AU,bool AA=false){
AE p;
p.AI=AK;
p.BB=plane;
p.BQ=cos(AU*acos(-1.0)/180.0);
p.W=W;
p.AT=CM;
p.AJ=AJ;
p.AA=AA;
GW.push_back(p);
}
;
double DO=0.036*d;
double DX=0.0108*d;
double DY=8.0;
double AO=88.0;
double proxy_threshold=0.955;
int AG=384;
bool relocate_on_first=false;
if(N<8000){
DO=0.048*d;
DX=0.0144*d;
DY=10.5;
AO=96.0;
proxy_threshold=0.935;
AG=512;
relocate_on_first=BL>=0.900&&Z<=0.035;
}
else if(N<30000){
DO=0.052*d;
DX=0.0156*d;
DY=11.5;
AO=95.0;
proxy_threshold=0.930;
AG=384;
relocate_on_first=BL>=0.930&&Z<=0.040;
}
else if(N<80000){
DO=0.050*d;
DX=0.0150*d;
DY=11.1;
AO=95.0;
proxy_threshold=0.910;
AG=128;
}
else if(N<140000){
DO=0.050*d;
DX=0.0150*d;
DY=11.1;
AO=95.0;
proxy_threshold=0.910;
AG=128;
}
else{
DO=0.056*d;
DX=0.0168*d;
DY=12.4;
AO=95.0;
proxy_threshold=0.920;
AG=96;
}
FZ(DO,DX,DY,relocate_on_first);
FZ(DO,DX,DY,false);
bool HV=false;
if(boost_start_vertices>0&&AX(GW,true)&&elapsed_seconds()<16.6){
const int boosted_vertices=count_output_vertices_estimate();
if(boosted_vertices>0&&(double)boosted_vertices*100.0<=(double)boost_start_vertices*AO){
HV=visual_proxy_score(AG)>=proxy_threshold;
}
}
if(!HV)restore_state(boost_safe_state);
}
if(FY&&elapsed_seconds()<12.5){
AP endpoint_state=AD();
const int endpoint_vertices=count_output_vertices_estimate();
restore_state(CG);
vector<AE>HF;
auto CY=[&](double AK,double plane,double AU,bool AA=false){
AE p;
p.AI=AK;
p.BB=plane;
p.BQ=cos(AU*acos(-1.0)/180.0);
p.W=W;
p.AT=CM;
p.AJ=AJ;
p.AA=AA;
HF.push_back(p);
}
;
double FN=0.024*d;
double FW=0.0072*d;
double FX=5.3;
if(N>=10000){
FN=0.0265*d;
FW=0.0079*d;
FX=5.9;
}
CY(0.004*d,0.0010*d,1.0);
CY(0.008*d,0.0020*d,1.8);
CY(0.012*d,0.0032*d,2.8);
CY(FN,FW,FX,true);
CY(FN,FW,FX,true);
bool HW=false;
if(AX(HF,true)&&elapsed_seconds()<15.5){
const int relocated_vertices=count_output_vertices_estimate();
if(relocated_vertices>0&&relocated_vertices*100<=endpoint_vertices*94){
const double proxy=visual_proxy_score(64);
HW=proxy>=0.970;
}
}
if(!HW)restore_state(endpoint_state);
}
if(elapsed_seconds()<13.0){
(void)GM();
}
if(elapsed_seconds()<13.0){
(void)FS();
}
if(elapsed_seconds()<13.8){
(void)GI();
}
if(elapsed_seconds()<13.8){
(void)FL();
}
if(elapsed_seconds()<13.8){
(void)HH();
}
if(elapsed_seconds()<13.8){
(void)EO();
}
if(elapsed_seconds()<13.8){
(void)GT();
}
if(AS&&N>=1000&&N<60000&&BG<=0.010&&AL>=0.760&&AH<=0.090&&elapsed_seconds()<12.6){
AP generic_safe_state=AD();
const int generic_start_vertices=count_output_vertices_estimate();
vector<AE>HX;
auto DF=[&](double AK,double plane,double AU,bool AA=false){
AE p;
p.AI=AK;
p.BB=plane;
p.BQ=cos(AU*acos(-1.0)/180.0);
p.W=W;
p.AT=CM;
p.AJ=AJ;
p.AA=AA;
HX.push_back(p);
}
;
double DZ=0.026*d;
double EK=0.0078*d;
double EL=6.0;
const double CB=96.0;
double proxy_threshold=0.925;
int AG=1024;
if(N>=6000){
DZ=0.034*d;
EK=0.0102*d;
EL=7.8;
proxy_threshold=0.925;
AG=1024;
}
if(N>=20000){
DZ=0.040*d;
EK=0.0120*d;
EL=9.0;
proxy_threshold=0.930;
AG=512;
}
if(N>=30000&&N<60000){
DZ=0.043*d;
EK=0.0129*d;
EL=9.6;
proxy_threshold=0.925;
AG=512;
}
if(Z>0.120||AH>0.050){
DZ*=0.82;
EK*=0.82;
EL*=0.86;
proxy_threshold+=0.020;
}
DF(0.004*d,0.0010*d,1.0,false);
DF(0.008*d,0.0020*d,1.8,false);
DF(0.012*d,0.0032*d,2.8,false);
DF(DZ,EK,EL,true);
DF(DZ,EK,EL,true);
bool IP=false;
if(generic_start_vertices>0&&AX(HX,true)&&elapsed_seconds()<16.4){
const int generic_vertices=count_output_vertices_estimate();
if(generic_vertices>0&&generic_vertices<generic_start_vertices&&(double)generic_vertices*100.0<=(double)generic_start_vertices*CB){
IP=visual_proxy_score(AG)>=proxy_threshold;
}
}
if(!IP)restore_state(generic_safe_state);
}
if(N>=900&&N<=60000&&AS&&BG<=0.012&&AL>=0.760&&AH<=0.120&&elapsed_seconds()<9.8){
AP QS=AD();
const int QN=count_output_vertices_estimate();
vector<Vec3>QV;
vector<Face>QF;
bool QK=false;
if(QN>0&&QX::build(originalP,AR,QV,QF)&&(int)QV.size()<QN&&((int)QV.size()*100<=QN*88)){
if(AF(QV,QF)&&elapsed_seconds()<16.7){
double QP=visual_proxy_score(512);
QK=QP>=0.935;
}
}
if(!QK)restore_state(QS);
}
if(N>=1000&&N<=120000&&elapsed_seconds()<13.2){
AP FM=AD();
const int DN=count_output_vertices_estimate();
vector<Vec3>CH;
vector<Face>EP;
const double heap_budget=(N<30000?4.2:(N<60000?5.6:5.2));
const bool W1=AS&&BG<=0.004&&AL>=0.995&&AH<=0.006;
bool EH=false;
if(DN>0&&N>=30000){
const double IK[]={
0.020,0.025,0.030,0.040,0.055,0.070,0.085,0.105,-1.0}
;
for(double ratio:IK){
if(ratio<0.025&&ratio>0.0&&!W1)continue;
if(N>=60000&&(ratio==0.085||ratio==0.105))continue;
if(elapsed_seconds()>=17.0)break;
CH.clear();
EP.clear();
restore_state(FM);
bool built=false;
if(ratio<0.0){
built=DC::EV(originalP,AR,CL,CH,EP,heap_budget);
}
else{
built=DC::EV(originalP,AR,CL,CH,EP,heap_budget,ratio);
}
if(!built||(int)CH.size()>=DN||elapsed_seconds()>=17.2)continue;
if(!AF(CH,EP)||elapsed_seconds()>=17.4)continue;
const double proxy=visual_proxy_score(1024);
const double proxy_threshold=(ratio<0.0?0.902:(ratio<0.025?0.940:((N<60000&&ratio>=0.085)?0.925:0.920)));
if(proxy>=proxy_threshold){
EH=true;
break;
}
}
}
else if(DN>0){
const double JG[]={
0.020,0.025,0.030,0.040,0.055,-1.0,0.080,0.115,0.180,0.350,0.500}
;
for(double ratio:JG){
if(ratio<0.025&&ratio>0.0&&!W1)continue;
if(elapsed_seconds()>=17.0)break;
CH.clear();
EP.clear();
restore_state(FM);
if(!DC::EV(originalP,AR,CL,CH,EP,heap_budget,ratio)){
continue;
}
if((int)CH.size()>=DN||elapsed_seconds()>=17.2)continue;
if(!AF(CH,EP)||elapsed_seconds()>=17.4)continue;
const double proxy=visual_proxy_score(1024);
const double proxy_threshold=(ratio<0.0?0.904:(ratio<0.025?0.940:0.920));
if(proxy>=proxy_threshold){
EH=true;
break;
}
}
}
if(!EH)restore_state(FM);
}
if(N>=30000&&N<60000&&AS&&BG<=0.012&&AL>=0.760&&AH<=0.120&&elapsed_seconds()<13.6){
AP W4H=AD();
const int W4HN=count_output_vertices_estimate();
vector<Vec3>W4HV;
vector<Face>W4HF;
const double W4HR[]={
0.135,0.160,0.200}
;
bool W4HK=false;
for(double ratio:W4HR){
if(elapsed_seconds()>=16.8)break;
W4HV.clear();
W4HF.clear();
restore_state(W4H);
if(!DC::EV(originalP,AR,CL,W4HV,W4HF,4.8,ratio))continue;
if(W4HN<=0||(int)W4HV.size()>=W4HN||elapsed_seconds()>=17.0)continue;
if((int)W4HV.size()*100>W4HN*94)continue;
if(!AF(W4HV,W4HF)||elapsed_seconds()>=17.2)continue;
const double W4HP=visual_proxy_score(1024);
const double W4HT=(ratio<=0.160?0.945:0.955);
if(W4HP>=W4HT){
W4HK=true;
break;
}
}
if(!W4HK)restore_state(W4H);
}
if(N>=300&&N<=70000&&elapsed_seconds()<12.9){
AP FI=AD();
const int DB=count_output_vertices_estimate();
auto HR=[&](const AE&params,int GG)->bool{
for(int rep=0;
rep<GG;
++rep){
if(elapsed_seconds()>16.2)return false;
vector<unsigned long long>IT;
IT.reserve((size_t)BE*3);
for(int fid=0;
fid<(int)faces.size();
++fid){
if((fid&8191)==0&&elapsed_seconds()>16.2)return false;
if(!BR[fid])continue;
const Face&f=faces[fid];
if(!BU[f.v[0]]||!BU[f.v[1]]||!BU[f.v[2]])continue;
IT.push_back(ED(f.v[0],f.v[1]));
IT.push_back(ED(f.v[1],f.v[2]));
IT.push_back(ED(f.v[2],f.v[0]));
}
sort(IT.begin(),IT.end());
IT.erase(unique(IT.begin(),IT.end()),IT.end());
struct HE{
double JQ;
unsigned long long key;
}
;
vector<HE>CP;
CP.reserve(IT.size());
for(unsigned long long key:IT){
int a=(int)(key>>32);
int b=(int)(key&0xffffffffu);
if(a>=0&&a<N&&b>=0&&b<N&&BU[a]&&BU[b]){
CP.push_back({
norm2(P[a]-P[b]),key}
);
}
}
sort(CP.begin(),CP.end(),[](const HE&lhs,const HE&rhs){
return lhs.JQ<rhs.JQ;
}
);
int HY=0;
for(const HE&e:CP){
if((HY&511)==0&&elapsed_seconds()>16.2)return false;
int a=(int)(e.key>>32);
int b=(int)(e.key&0xffffffffu);
if(a>=0&&a<N&&b>=0&&b<N&&BU[a]&&BU[b]){
if(GD(a,b,params))++HY;
}
}
if(HY==0)break;
}
return true;
}
;
if(DB>0){
struct HT{
double AK;
double plane;
double AU;
int GG;
int AG;
double proxy_threshold;
double CB;
bool JO;
bool JV;
}
;
vector<HT>AQ;
const bool FB=AS&&BG<=0.012&&AL>=0.760&&AH<=0.110;
if(N<1000){
AQ.push_back({
0.036*d,0.0140*d,10.0,2,768,0.928,99.0,true,false}
);
AQ.push_back({
0.048*d,0.0270*d,20.0,4,768,FB?0.906:0.925,98.0,true,true}
);
}
else if(N<3000){
AQ.push_back({
0.040*d,0.0180*d,13.5,3,1024,0.916,99.0,true,false}
);
AQ.push_back({
0.049*d,0.0340*d,28.0,6,1024,FB?0.903:0.922,98.0,true,true}
);
}
else if(N<8000){
AQ.push_back({
0.044*d,0.0220*d,17.0,4,1024,0.913,99.0,true,false}
);
AQ.push_back({
0.049*d,0.0400*d,34.0,7,1024,FB?0.904:0.924,98.0,true,true}
);
}
else if(N<25000){
AQ.push_back({
0.047*d,0.0250*d,19.0,4,1024,0.914,99.0,true,false}
);
AQ.push_back({
0.049*d,0.0410*d,34.0,6,1024,FB?0.902:0.930,98.0,true,true}
);
}
else if(N<60000){
AQ.push_back({
0.047*d,0.0230*d,17.5,3,512,0.925,99.0,true,false}
);
AQ.push_back({
0.049*d,0.0350*d,28.0,4,512,FB?0.914:0.940,98.0,true,true}
);
AQ.push_back({
0.050*d,0.0380*d,30.0,5,512,0.932,97.0,true,true}
);
}
else if(N<200000){
AQ.push_back({
0.047*d,0.0220*d,16.0,3,256,0.938,99.0,true,false}
);
AQ.push_back({
0.049*d,0.0320*d,24.0,4,256,FB?0.930:0.948,98.0,true,N<120000}
);
}
else{
AQ.push_back({
0.047*d,0.0220*d,16.0,3,128,0.950,99.0,true,false}
);
AQ.push_back({
0.049*d,0.0300*d,22.0,4,128,FB?0.940:0.955,98.0,true,false}
);
}
AP IV=FI;
int IF=DB;
for(const HT&AZ:AQ){
if(elapsed_seconds()>14.5)break;
restore_state(FI);
AE p;
p.AI=AZ.AK;
p.BB=AZ.plane;
p.BQ=cos(AZ.AU*acos(-1.0)/180.0);
p.W=W;
p.AT=CM;
p.AJ=AJ;
p.AA=AZ.JO;
bool ran=false;
if(AZ.JV){
vector<AE>JR;
AE w=p;
w.AI=min(p.AI,0.020*d);
w.BB=min(p.BB,0.0060*d);
w.BQ=cos(min(AZ.AU,4.8)*acos(-1.0)/180.0);
w.AA=false;
JR.push_back(w);
ran=AX(JR,true)&&HR(p,AZ.GG);
}
else{
vector<AE>HZ;
for(int i=0;
i<AZ.GG;
++i)HZ.push_back(p);
ran=AX(HZ,true);
}
bool IZ=false;
if(ran&&elapsed_seconds()<16.5){
const int BO=count_output_vertices_estimate();
if(BO>0&&BO<IF&&(double)BO*100.0<=(double)DB*AZ.CB){
const double proxy=visual_proxy_score(AZ.AG);
IZ=proxy>=AZ.proxy_threshold;
}
if(IZ){
IV=AD();
IF=BO;
}
}
}
restore_state(IV);
}
}
if(elapsed_seconds()<16.6){
(void)FR();
}
if(elapsed_seconds()<16.6){
(void)FG();
}
if(elapsed_seconds()<16.8){
(void)IC();
}
}
static void CA(string&out,const char*line,int len){
if(out.size()+(size_t)len>(1<<20)){
fwrite(out.data(),1,out.size(),stdout);
out.clear();
}
out.append(line,line+len);
}
static void IJ(){
string out;
out.reserve(1<<20);
char line[128];
int len=snprintf(line,sizeof(line),"%d %d\n",N,M);
CA(out,line,len);
for(int i=0;
i<N;
++i){
len=snprintf(line,sizeof(line),"v %.10g %.10g %.10g\n",originalP[i].x,originalP[i].y,originalP[i].z);
CA(out,line,len);
}
for(int fid=0;
fid<M;
++fid){
const Face&f=AR[fid];
len=snprintf(line,sizeof(line),"f %d %d %d\n",f.v[0]+1,f.v[1]+1,f.v[2]+1);
CA(out,line,len);
}
if(!out.empty())fwrite(out.data(),1,out.size(),stdout);
}
static bool HS(vector<int>&EC,vector<Face>&AB){
vector<int>id(N,-1);
EC.clear();
AB.clear();
EC.reserve(N);
AB.reserve(BE);
const double IY=max(1e-30,1e-24*CL*CL);
for(int fid=0;
fid<(int)faces.size();
++fid){
if(!BR[fid])continue;
const Face&f=faces[fid];
for(int k=0;
k<3;
++k){
if(f.v[k]<0||f.v[k]>=N||!BU[f.v[k]])return false;
}
if(f.v[0]==f.v[1]||f.v[0]==f.v[2]||f.v[1]==f.v[2])return false;
Vec3 cr=II(fid);
if(!(norm2(cr)>IY))return false;
Face JK;
for(int k=0;
k<3;
++k){
int old=f.v[k];
if(id[old]<0){
id[old]=(int)EC.size();
EC.push_back(old);
}
JK.v[k]=id[old];
}
AB.push_back(JK);
}
if(EC.empty()||AB.empty()||(int)EC.size()>N)return false;
vector<unsigned long long>CP;
CP.reserve(AB.size()*3);
for(const Face&f:AB){
CP.push_back(ED(f.v[0],f.v[1]));
CP.push_back(ED(f.v[1],f.v[2]));
CP.push_back(ED(f.v[2],f.v[0]));
}
sort(CP.begin(),CP.end());
for(size_t i=0;
i<CP.size();
){
size_t j=i+1;
while(j<CP.size()&&CP[j]==CP[i])++j;
if(j-i!=2)return false;
i=j;
}
return true;
}
static void JD(){
vector<int>EC;
vector<Face>AB;
if(!HS(EC,AB)){
IJ();
return;
}
string out;
out.reserve(1<<20);
char line[128];
int len=snprintf(line,sizeof(line),"%d %d\n",(int)EC.size(),(int)AB.size());
CA(out,line,len);
const bool GV=(int)EC.size()*2<=N;
for(int old:EC){
if(GV){
len=snprintf(line,sizeof(line),"v %.15g %.15g %.15g\n",P[old].x,P[old].y,P[old].z);
}
else{
len=snprintf(line,sizeof(line),"v %.10g %.10g %.10g\n",P[old].x,P[old].y,P[old].z);
}
CA(out,line,len);
}
for(const Face&f:AB){
len=snprintf(line,sizeof(line),"f %d %d %d\n",f.v[0]+1,f.v[1]+1,f.v[2]+1);
CA(out,line,len);
}
if(!out.empty())fwrite(out.data(),1,out.size(),stdout);
}
int main(){
JC();
GN();
JD();
return 0;
}
