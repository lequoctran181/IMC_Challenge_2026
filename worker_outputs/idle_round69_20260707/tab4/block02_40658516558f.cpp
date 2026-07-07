#include <bits/stdc++.h>
using namespace std;

struct V{
double x,y,z;
}
;

static inline V operator+(const V&a,const V&b){
return{
a.x+b.x,a.y+b.y,a.z+b.z}
;
}

static inline V operator-(const V&a,const V&b){
return{
a.x-b.x,a.y-b.y,a.z-b.z}
;
}

static inline V operator*(const V&a,double s){
return{
a.x*s,a.y*s,a.z*s}
;
}

static inline V operator/(const V&a,double s){
return{
a.x/s,a.y/s,a.z/s}
;
}

static inline double dotp(const V&a,const V&b){
return a.x*b.x+a.y*b.y+a.z*b.z;
}

static inline V crossp(const V&a,const V&b){
return{
a.y*b.z-a.z*b.y,a.z*b.x-a.x*b.z,a.x*b.y-a.y*b.x}
;
}

static inline double n2(const V&a){
return dotp(a,a);
}
static inline double norm(const V&a){
return sqrt(max(0.0,n2(a)));
}

static inline V unit(V a){
double l=norm(a);
return l>1e-300?a/l:V{
1,0,0}
;
}

static inline double clampd(double x,double a,double b){
return x<a?a:(x>b?b:x);
}
 
struct F{
int a,b,c;
}
;
struct Mesh{
vector<V>v;
vector<F>f;
}
;

int N,M;
Mesh Orig;
V bminv,bmaxv,ctr;
double diagv=1,tau=1,tau2=1;
chrono::steady_clock::time_point T0;

static inline double tim(){
return chrono::duration<double>(chrono::steady_clock::now()-T0).count();
}

struct In{
vector<char>b;
char*p;
In(){
char t[1<<16];
size_t n;
b.reserve(1<<26);
while((n=fread(t,1,sizeof(t),stdin))>0)b.insert(b.end(),t,t+n);
b.push_back(0);
p=b.data();
}
inline void skip(){
while(*p==' '||*p=='\n'||*p=='\r'||*p=='\t')++p;
}
int ni(){
skip();
int s=1;
if(*p=='-')s=-1,++p;
int x=0;
while(*p>='0'&&*p<='9')x=x*10+*p++-'0';
return x*s;
}
double nd(){
skip();
char*e;
double x=strtod(p,&e);
p=e;
return x;
}
char nc(){
skip();
return *p?*p++:0;
}
}
;

static inline uint64_t ek(int a,int b){
if(a>b)swap(a,b);
return(uint64_t)(uint32_t)a<<32|(uint32_t)b;
}

static inline array<int,3> fk(int a,int b,int c){
array<int,3>x{
a,b,c}
;
sort(x.begin(),x.end());
return x;
}

static inline bool samef(const F&f,int a,int b,int c){
return fk(f.a,f.b,f.c)==fk(a,b,c);
}
 
static V fcross(const Mesh&m,const F&f){
return crossp(m.v[f.b]-m.v[f.a],m.v[f.c]-m.v[f.a]);
}

static double vol6(const Mesh&m){
long double s=0;
for(auto &f:m.f)s+=dotp(m.v[f.a],crossp(m.v[f.b],m.v[f.c]));
return(double)s;
}

static void orientLike(Mesh&m){
if(vol6(Orig)*vol6(m)<0)for(auto &f:m.f)swap(f.b,f.c);
}
 
static bool read(){
In in;
N=in.ni();
M=in.ni();
if(N<=0||M<=0)return false;
Orig.v.resize(N);
Orig.f.resize(M);
bminv={
1e100,1e100,1e100}
;
bmaxv={
-1e100,-1e100,-1e100}
;
for(int i=0;
i<N;
i++){
char c=in.nc();
if(c!='v'&&c!='V')--in.p;
V p{
in.nd(),in.nd(),in.nd()}
;
Orig.v[i]=p;
bminv.x=min(bminv.x,p.x);
bminv.y=min(bminv.y,p.y);
bminv.z=min(bminv.z,p.z);
bmaxv.x=max(bmaxv.x,p.x);
bmaxv.y=max(bmaxv.y,p.y);
bmaxv.z=max(bmaxv.z,p.z);
}
for(int i=0;
i<M;
i++){
char c=in.nc();
if(c!='f'&&c!='F')--in.p;
Orig.f[i]={
in.ni()-1,in.ni()-1,in.ni()-1}
;
}
diagv=norm(bmaxv-bminv);
if(!(diagv>0))diagv=1;
ctr=(bminv+bmaxv)*0.5;
tau=0.05*diagv*0.996;
tau2=tau*tau;
return true;
}

static bool valid(const Mesh&m){
int n=m.v.size();
if(n<=0||m.f.empty()||n>N)return false;
double ae=max(1e-30,1e-24*diagv*diagv);
vector<unsigned char>use(n);
vector<uint64_t>ed;
ed.reserve(m.f.size()*3);
vector<array<int,3>>fs;
fs.reserve(m.f.size());
for(auto &f:m.f){
if(f.a<0||f.b<0||f.c<0||f.a>=n||f.b>=n||f.c>=n||f.a==f.b||f.a==f.c||f.b==f.c)return false;
if(n2(fcross(m,f))<=ae)return false;
use[f.a]=use[f.b]=use[f.c]=1;
ed.push_back(ek(f.a,f.b));
ed.push_back(ek(f.b,f.c));
ed.push_back(ek(f.c,f.a));
fs.push_back(fk(f.a,f.b,f.c));
}
for(int i=0;
i<n;
i++)if(!use[i])return false;
sort(fs.begin(),fs.end());
if(adjacent_find(fs.begin(),fs.end())!=fs.end())return false;
sort(ed.begin(),ed.end());
for(size_t i=0;
i<ed.size();
){
size_t j=i+1;
while(j<ed.size()&&ed[j]==ed[i])j++;
if(j-i!=2)return false;
i=j;
}
return true;
}

static void printM(const Mesh&m){
printf("%d %d\n",(int)m.v.size(),(int)m.f.size());
for(auto&p:m.v)printf("v %.15g %.15g %.15g\n",p.x,p.y,p.z);
for(auto&f:m.f)printf("f %d %d %d\n",f.a+1,f.b+1,f.c+1);
}
 

struct Render{
int R;
vector<float>d,n;
}
;

static inline void proj(const V&p,int view,int R,double&x,double&y,double&z){
double D=2.5,f=800.0*R/1024.0,c=0.5*R,sx,sy;
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
x=f*sx/z+c;
y=f*sy/z+c;
}

static void render(const Mesh&m,int R,Render&rm){
int P=R*R;
rm.R=R;
rm.d.assign(6*P,255.0f);
rm.n.assign(6*P*3,127.5f);
vector<float>X(m.v.size()),Y(m.v.size()),Z(m.v.size());
for(int view=0;
view<6;
view++){
for(size_t i=0;
i<m.v.size();
i++){
double x,y,z;
proj(m.v[i],view,R,x,y,z);
X[i]=x;
Y[i]=y;
Z[i]=z;
}
float*db=rm.d.data()+view*P;
float*nb=rm.n.data()+view*P*3;
for(auto&f:m.f){
int ia=f.a,ib=f.b,ic=f.c;
float x0=X[ia],x1=X[ib],x2=X[ic],y0=Y[ia],y1=Y[ib],y2=Y[ic],z0=Z[ia],z1=Z[ib],z2=Z[ic];
if(!(z0>0&&z1>0&&z2>0))continue;
int xmin=max(0,(int)floor(min(x0,min(x1,x2))-0.5)),xmax=min(R-1,(int)ceil(max(x0,max(x1,x2))+0.5));
int ymin=max(0,(int)floor(min(y0,min(y1,y2))-0.5)),ymax=min(R-1,(int)ceil(max(y0,max(y1,y2))+0.5));
if(xmin>xmax||ymin>ymax)continue;
float den=(y1-y2)*(x0-x2)+(x2-x1)*(y0-y2);
if(fabs(den)<1e-20)continue;
V nn=unit(fcross(m,f));
float nr=(nn.x+1)*127.5,ng=(nn.y+1)*127.5,nbv=(nn.z+1)*127.5,iv=1.0f/den;
for(int yy=ymin;
yy<=ymax;
yy++){
float py=yy+0.5f;
int row=yy*R;
for(int xx=xmin;
xx<=xmax;
xx++){
float px=xx+0.5f;
float w0=((y1-y2)*(px-x2)+(x2-x1)*(py-y2))*iv;
float w1=((y2-y0)*(px-x2)+(x0-x2)*(py-y2))*iv;
float w2=1-w0-w1;
if(w0>=-1e-6f&&w1>=-1e-6f&&w2>=-1e-6f){
float dep=1.0f/(w0/z0+w1/z1+w2/z2);
int id=row+xx;
if(dep<db[id]){
db[id]=dep;
float*q=nb+id*3;
q[0]=nr;
q[1]=ng;
q[2]=nbv;
}
}
}
}
}
}
}

static inline double rsum(const vector<double>&I,int W,int x0,int y0,int x1,int y1){
return I[(size_t)y1*W+x1]-I[(size_t)y0*W+x1]-I[(size_t)y1*W+x0]+I[(size_t)y0*W+x0];
}

static double ssim(const float*A,int sa,const float*B,int sb,const float*DA,const float*DB,int R,vector<double>&IA,vector<double>&IB,vector<double>&IA2,vector<double>&IB2,vector<double>&IAB){
int W=R+1;
size_t SZ=(size_t)W*W;
fill(IA.begin(),IA.begin()+SZ,0);
fill(IB.begin(),IB.begin()+SZ,0);
fill(IA2.begin(),IA2.begin()+SZ,0);
fill(IB2.begin(),IB2.begin()+SZ,0);
fill(IAB.begin(),IAB.begin()+SZ,0);
for(int y=1;
y<=R;
y++){
double a=0,b=0,a2=0,b2=0,ab=0;
int row=(y-1)*R;
for(int x=1;
x<=R;
x++){
int p=row+x-1;
double av=A[(size_t)p*sa],bv=B[(size_t)p*sb];
a+=av;
b+=bv;
a2+=av*av;
b2+=bv*bv;
ab+=av*bv;
size_t id=(size_t)y*W+x,up=(size_t)(y-1)*W+x;
IA[id]=IA[up]+a;
IB[id]=IB[up]+b;
IA2[id]=IA2[up]+a2;
IB2[id]=IB2[up]+b2;
IAB[id]=IAB[up]+ab;
}
}
int rad=5,area=121;
double c1=6.5025,c2=58.5225;
long long cnt=0;
long double acc=0;
for(int y=rad;
y<R-rad;
y++)for(int x=rad;
x<R-rad;
x++){
int p=y*R+x;
if(!(DA[p]<254.0f||DB[p]<254.0f))continue;
int x0=x-rad,y0=y-rad,x1=x+rad+1,y1=y+rad+1;
double sa1=rsum(IA,W,x0,y0,x1,y1),sb1=rsum(IB,W,x0,y0,x1,y1),sa2=rsum(IA2,W,x0,y0,x1,y1),sb2=rsum(IB2,W,x0,y0,x1,y1),sab=rsum(IAB,W,x0,y0,x1,y1);
double ma=sa1/area,mb=sb1/area,va=max(0.0,sa2/area-ma*ma),vb=max(0.0,sb2/area-mb*mb),cv=sab/area-ma*mb;
double den=(ma*ma+mb*mb+c1)*(va+vb+c2);
acc+=den?((2*ma*mb+c1)*(2*cv+c2)/den):1.0;
cnt++;
}
return cnt?(double)(acc/cnt):1.0;
}

static double scoreR(const Render&A,const Render&B){
int R=A.R,P=R*R,W=R+1;
vector<double>I1((size_t)W*W),I2(I1.size()),I3(I1.size()),I4(I1.size()),I5(I1.size());
double tot=0;
for(int v=0;
v<6;
v++){
const float*ad=A.d.data()+v*P,*bd=B.d.data()+v*P;
double ns=0;
for(int c=0;
c<3;
c++)ns+=ssim(A.n.data()+v*P*3+c,3,B.n.data()+v*P*3+c,3,ad,bd,R,I1,I2,I3,I4,I5);
ns/=3;
double ds=ssim(ad,1,bd,1,ad,bd,R,I1,I2,I3,I4,I5);
tot+=0.5*(ns+ds);
}
return tot/6;
}

static int resFor(){
if(N<2500)return 512;
if(N<30000)return 384;
if(N<120000)return 224;
if(N<350000)return 144;
return 96;
}

static double guardFor(){
if(N<2500)return .930;
if(N<30000)return .925;
if(N<120000)return .920;
if(N<350000)return .915;
return .910;
}


struct GH{
double cell,inv,r2;
unordered_map<long long,vector<int>>mp;
const vector<V>*P;
static uint64_t mix(uint64_t x){
x+=0x9e3779b97f4a7c15ULL;
x=(x^(x>>30))*0xbf58476d1ce4e5b9ULL;
x=(x^(x>>27))*0x94d049bb133111ebULL;
return x^(x>>31);
}
static long long key(int a,int b,int c){
uint64_t h=mix((uint64_t)(int64_t)a);
h^=mix((uint64_t)(int64_t)b)*3;
h^=mix((uint64_t)(int64_t)c)*7;
return(long long)h;
}
void addp(int i,const V&p){
int x=floor(p.x*inv),y=floor(p.y*inv),z=floor(p.z*inv);
mp[key(x,y,z)].push_back(i);
}
void build(const vector<V>&v,double r){
P=&v;
cell=max(r,1e-30);
inv=1.0/cell;
r2=r*r;
mp.clear();
mp.reserve(v.size()*2+7);
for(int i=0;
i<(int)v.size();
i++)addp(i,v[i]);
}
bool near(const V&q)const{
int x=floor(q.x*inv),y=floor(q.y*inv),z=floor(q.z*inv);
for(int dx=-1;
dx<=1;
dx++)for(int dy=-1;
dy<=1;
dy++)for(int dz=-1;
dz<=1;
dz++){
auto it=mp.find(key(x+dx,y+dy,z+dz));
if(it==mp.end())continue;
for(int id:it->second)if(n2((*P)[id]-q)<=r2)return true;
}
return false;
}
}
;

static vector<V> vertexNormals(){
vector<V> n(N,{
0,0,0}
);
for(auto&f:Orig.f){
V cr=fcross(Orig,f);
n[f.a]=n[f.a]+cr;
n[f.b]=n[f.b]+cr;
n[f.c]=n[f.c]+cr;
}
double s=vol6(Orig)>=0?1.0:-1.0;
for(int i=0;
i<N;
i++){
V u=unit(n[i]*s);
if(norm(n[i])<1e-20){
V r=Orig.v[i]-ctr;
u=unit(r);
}
n[i]=u;
}
return n;
}

static void addTet(Mesh&m,const V&base,const V&outn){
V n=unit(outn);
V a=fabs(n.x)<0.8?V{
1,0,0}
:V{
0,1,0}
;
V t1=unit(crossp(n,a));
V t2=unit(crossp(n,t1));
double e=max(1e-10*diagv,1e-5*tau);
int s=m.v.size();
m.v.push_back(base);
m.v.push_back(base+t1*e);
m.v.push_back(base+t2*e);
m.v.push_back(base-n*e);
m.f.push_back({
s,s+2,s+1}
);
m.f.push_back({
s,s+1,s+3}
);
m.f.push_back({
s,s+3,s+2}
);
m.f.push_back({
s+1,s+2,s+3}
);
}

static bool addWitnesses(Mesh&cand,const vector<V>&vn,int*addedOut=nullptr){
double shift=0.030*diagv;
double eps=max(1e-12,0.05*diagv*0.997);
GH gh;
gh.build(cand.v,eps);
vector<V> centers,norms;
centers.reserve(4096);
norms.reserve(4096);
vector<V> all=cand.v;
all.reserve(cand.v.size()+N/10);
for(int i=0;
i<N;
i++){
if(gh.near(Orig.v[i]))continue;
V c=Orig.v[i]-vn[i]*shift;
centers.push_back(c);
norms.push_back(vn[i]);
all.push_back(c);
gh.addp((int)all.size()-1,c);
if((int)cand.v.size()+4*(int)centers.size()>N)return false;
if(tim()>19.0)return false;
}
for(size_t i=0;
i<centers.size();
i++)addTet(cand,centers[i],norms[i]);
if(addedOut)*addedOut=centers.size();
if((int)cand.v.size()>N)return false;
return valid(cand);
}

static bool outputToOriginalOK(const Mesh&m){
GH gh;
gh.build(Orig.v,tau);
for(auto&p:m.v)if(!gh.near(p))return false;
return true;
}


static bool detectLat(int&R,int&S){
if(N<10||M!=2*(N-2))return false;
for(int s=8;
s<=4096;
s++){
if((N-2)%s)continue;
int r=(N-2)/s;
if(r<3)continue;
bool ok=1;
int step=max(1,s/80);
for(int j=0;
j<s&&ok;
j+=step){
int a=1+j,b=1+(j+1)%s;
if(!samef(Orig.f[j],0,b,a))ok=0;
}
int span=max(1,(r-1)*s/256);
for(int q=0;
q<(r-1)*s&&ok;
q+=span){
int rr=q/s,j=q%s;
int a=1+rr*s+j,b=1+rr*s+(j+1)%s,c=1+(rr+1)*s+j,d=1+(rr+1)*s+(j+1)%s,f=s+2*q;
if(f+1>=M||!samef(Orig.f[f],a,b,c)||!samef(Orig.f[f+1],b,d,c))ok=0;
}
int bot=s+2*(r-1)*s;
for(int j=0;
j<s&&ok;
j+=step){
int c=1+(r-1)*s+j,d=1+(r-1)*s+(j+1)%s;
if(bot+j>=M||!samef(Orig.f[bot+j],N-1,c,d))ok=0;
}
if(ok){
R=r;
S=s;
return true;
}
}
return false;
}

static Mesh buildLat(int R,int S,int r2,int s2){
Mesh m;
if(r2<2||s2<8)return m;
m.v.push_back(Orig.v[0]);
for(int i=0;
i<r2;
i++){
int oi=1+(int)llround((double)i*(R-1)/max(1,r2-1));
oi=max(1,min(R,oi));
for(int j=0;
j<s2;
j++){
int oj=(int)((long long)j*S/s2);
m.v.push_back(Orig.v[1+(oi-1)*S+oj]);
}
}
int bot=m.v.size();
m.v.push_back(Orig.v[N-1]);
auto id=[&](int i,int j){
return 1+i*s2+((j%s2+s2)%s2);
}
;
auto ad=[&](int a,int b,int c){
m.f.push_back({
a,b,c}
);
}
;
for(int j=0;
j<s2;
j++)ad(0,id(0,j+1),id(0,j));
for(int i=0;
i+1<r2;
i++)for(int j=0;
j<s2;
j++){
int a=id(i,j),b=id(i,j+1),c=id(i+1,j),d=id(i+1,j+1);
ad(a,b,c);
ad(b,d,c);
}
for(int j=0;
j<s2;
j++)ad(bot,id(r2-1,j),id(r2-1,j+1));
orientLike(m);
return m;
}

static bool detectTor(int&U,int&S){
if(N<36||M!=2*N)return false;
for(int s=6;
s<=2048;
s++){
if(N%s)continue;
int u=N/s;
if(u<6)continue;
bool ok=1;
int step=max(1,N/700),cnt=0;
for(int q=0;
q<N&&cnt<700&&ok;
q+=step,cnt++){
int i=q/s,j=q%s,ni=(i+1==u?0:i+1),nj=(j+1==s?0:j+1);
int a=i*s+j,b=i*s+nj,c=ni*s+j,d=ni*s+nj,f=2*q;
if(f+1>=M){
ok=0;
break;
}
bool m0=samef(Orig.f[f],a,b,c)||samef(Orig.f[f],b,d,c)||samef(Orig.f[f],a,c,d)||samef(Orig.f[f],a,b,d);
bool m1=samef(Orig.f[f+1],a,b,c)||samef(Orig.f[f+1],b,d,c)||samef(Orig.f[f+1],a,c,d)||samef(Orig.f[f+1],a,b,d);
if(!(m0&&m1))ok=0;
}
if(ok){
U=u;
S=s;
return true;
}
}
return false;
}

static Mesh buildTor(int U,int S,int u2,int s2){
Mesh m;
if(u2<4||s2<4)return m;
for(int i=0;
i<u2;
i++){
int oi=(int)((long long)i*U/u2);
for(int j=0;
j<s2;
j++){
int oj=(int)((long long)j*S/s2);
m.v.push_back(Orig.v[oi*S+oj]);
}
}
auto id=[&](int i,int j){
return ((i%u2+u2)%u2)*s2+((j%s2+s2)%s2);
}
;
for(int i=0;
i<u2;
i++)for(int j=0;
j<s2;
j++){
int a=id(i,j),b=id(i,j+1),c=id(i+1,j),d=id(i+1,j+1);
m.f.push_back({
a,b,c}
);
m.f.push_back({
b,d,c}
);
}
orientLike(m);
return m;
}


struct Ell{
bool ok=0;
V c,ax[3];
double r[3],rms=1e9,mx=1e9;
}
;

static void jac(double A[3][3],V ax[3]){
double v[3][3]={
{
1,0,0}
,{
0,1,0}
,{
0,0,1}
}
;
for(int it=0;
it<50;
it++){
int p=0,q=1;
double best=fabs(A[0][1]);
if(fabs(A[0][2])>best)p=0,q=2,best=fabs(A[0][2]);
if(fabs(A[1][2])>best)p=1,q=2,best=fabs(A[1][2]);
if(best<1e-18)break;
double app=A[p][p],aqq=A[q][q],apq=A[p][q],tau=(aqq-app)/(2*apq);
double t=(tau>=0?1.0:-1.0)/(fabs(tau)+sqrt(1+tau*tau)),c=1/sqrt(1+t*t),s=t*c;
for(int k=0;
k<3;
k++)if(k!=p&&k!=q){
double akp=A[k][p],akq=A[k][q];
A[k][p]=A[p][k]=c*akp-s*akq;
A[k][q]=A[q][k]=s*akp+c*akq;
}
A[p][p]=c*c*app-2*s*c*apq+s*s*aqq;
A[q][q]=s*s*app+2*s*c*apq+c*c*aqq;
A[p][q]=A[q][p]=0;
for(int k=0;
k<3;
k++){
double vp=v[k][p],vq=v[k][q];
v[k][p]=c*vp-s*vq;
v[k][q]=s*vp+c*vq;
}
}
int o[3]={
0,1,2}
;
sort(o,o+3,[&](int a,int b){
return A[a][a]>A[b][b];
}
);
for(int i=0;
i<3;
i++){
int col=o[i];
ax[i]=unit({
v[0][col],v[1][col],v[2][col]}
);
}
if(dotp(crossp(ax[0],ax[1]),ax[2])<0)ax[2]=ax[2]*-1;
}

static Ell fitEll(){
Ell e;
if(N<200)return e;
V mean{
0,0,0}
;
int st=max(1,N/240000),cnt=0;
for(int i=0;
i<N;
i+=st)mean=mean+Orig.v[i],cnt++;
mean=mean/max(1,cnt);
double C[3][3]{
}
;
for(int i=0;
i<N;
i+=st){
V q=Orig.v[i]-mean;
double x[3]={
q.x,q.y,q.z}
;
for(int a=0;
a<3;
a++)for(int b=0;
b<3;
b++)C[a][b]+=x[a]*x[b];
}
for(int a=0;
a<3;
a++)for(int b=0;
b<3;
b++)C[a][b]/=max(1,cnt);
jac(C,e.ax);
double lo[3]={
1e100,1e100,1e100}
,hi[3]={
-1e100,-1e100,-1e100}
;
for(auto&p:Orig.v)for(int k=0;
k<3;
k++){
double t=dotp(p,e.ax[k]);
lo[k]=min(lo[k],t);
hi[k]=max(hi[k],t);
}
e.c={
0,0,0}
;
for(int k=0;
k<3;
k++){
double mid=(lo[k]+hi[k])*0.5;
e.r[k]=(hi[k]-lo[k])*0.5;
if(e.r[k]<=1e-12)return e;
e.c=e.c+e.ax[k]*mid;
}
double mr=max(e.r[0],max(e.r[1],e.r[2])),nr=min(e.r[0],min(e.r[1],e.r[2]));
if(nr<mr*.12)return e;
double ss=0,mx=0;
cnt=0;
for(int i=0;
i<N;
i+=st){
V q=Orig.v[i]-e.c;
double s=0;
for(int k=0;
k<3;
k++){
double u=dotp(q,e.ax[k])/e.r[k];
s+=u*u;
}
double er=fabs(sqrt(max(0.0,s))-1);
ss+=er*er;
mx=max(mx,er);
cnt++;
}
e.rms=sqrt(ss/max(1,cnt));
e.mx=mx;
e.ok=e.rms<(N<5000?.018:.012)&&e.mx<(N<5000?.080:.055);
return e;
}

static Mesh buildEll(const Ell&e,int lat,int lon){
Mesh m;
if(!e.ok||lat<4||lon<8)return m;
auto pt=[&](double x,double y,double z){
return e.c+e.ax[0]*(e.r[0]*x)+e.ax[1]*(e.r[1]*y)+e.ax[2]*(e.r[2]*z);
}
;
m.v.push_back(pt(0,0,1));
for(int i=1;
i<lat;
i++){
double th=M_PI*i/lat,st=sin(th),ct=cos(th);
for(int j=0;
j<lon;
j++){
double ph=2*M_PI*j/lon;
m.v.push_back(pt(st*cos(ph),st*sin(ph),ct));
}
}
int bot=m.v.size();
m.v.push_back(pt(0,0,-1));
auto id=[&](int r,int j){
return 1+(r-1)*lon+((j%lon+lon)%lon);
}
;
for(int j=0;
j<lon;
j++)m.f.push_back({
0,id(1,j+1),id(1,j)}
);
for(int r=1;
r<lat-1;
r++)for(int j=0;
j<lon;
j++){
int a=id(r,j),b=id(r,j+1),c=id(r+1,j),d=id(r+1,j+1);
m.f.push_back({
a,b,c}
);
m.f.push_back({
b,d,c}
);
}
for(int j=0;
j<lon;
j++)m.f.push_back({
bot,id(lat-1,j),id(lat-1,j+1)}
);
orientLike(m);
return m;
}

static Mesh buildBox(int S){
Mesh m;
map<tuple<int,int,int>,int>mp;
V ext=bmaxv-bminv;
auto get=[&](int i,int j,int k){
auto key=make_tuple(i,j,k);
auto it=mp.find(key);
if(it!=mp.end())return it->second;
V p{
bminv.x+ext.x*i/S,bminv.y+ext.y*j/S,bminv.z+ext.z*k/S}
;
int id=m.v.size();
mp[key]=id;
m.v.push_back(p);
return id;
}
;
auto add=[&](int a,int b,int c){
m.f.push_back({
a,b,c}
);
}
;
for(int i=0;
i<S;
i++)for(int j=0;
j<S;
j++){
int a=get(i,j,0),b=get(i+1,j,0),c=get(i+1,j+1,0),d=get(i,j+1,0);
add(a,c,b);
add(a,d,c);
a=get(i,j,S);
b=get(i+1,j,S);
c=get(i+1,j+1,S);
d=get(i,j+1,S);
add(a,b,c);
add(a,c,d);
}
for(int j=0;
j<S;
j++)for(int k=0;
k<S;
k++){
int a=get(0,j,k),b=get(0,j+1,k),c=get(0,j+1,k+1),d=get(0,j,k+1);
add(a,c,b);
add(a,d,c);
a=get(S,j,k);
b=get(S,j+1,k);
c=get(S,j+1,k+1);
d=get(S,j,k+1);
add(a,b,c);
add(a,c,d);
}
for(int i=0;
i<S;
i++)for(int k=0;
k<S;
k++){
int a=get(i,0,k),b=get(i+1,0,k),c=get(i+1,0,k+1),d=get(i,0,k+1);
add(a,b,c);
add(a,c,d);
a=get(i,S,k);
b=get(i+1,S,k);
c=get(i+1,S,k+1);
d=get(i,S,k+1);
add(a,c,b);
add(a,d,c);
}
orientLike(m);
return m;
}


struct Simp{
vector<V>P;
vector<F>Fv;
vector<vector<int>>inc;
vector<unsigned char>av,af;
int actv,actf;
Simp(){
P=Orig.v;
Fv=Orig.f;
av.assign(N,1);
af.assign(M,1);
actv=N;
actf=M;
rebuild();
}
bool hs(int id,int v){
auto&f=Fv[id];
return f.a==v||f.b==v||f.c==v;
}
void rebuild(){
inc.assign(N,{
}
);
vector<int>d(N);
for(int i=0;
i<M;
i++)if(af[i])d[Fv[i].a]++,d[Fv[i].b]++,d[Fv[i].c]++;
for(int i=0;
i<N;
i++)if(av[i])inc[i].reserve(d[i]+4);
for(int i=0;
i<M;
i++)if(af[i])inc[Fv[i].a].push_back(i),inc[Fv[i].b].push_back(i),inc[Fv[i].c].push_back(i);
}
bool ef(int u,int v,int E[2],int O[2]){
int c=0;
for(int id:inc[u])if(af[id]&&hs(id,v)){
if(c>=2)return false;
E[c]=id;
auto&f=Fv[id];
int a[3]={
f.a,f.b,f.c}
;
O[c]=-1;
for(int x:a)if(x!=u&&x!=v)O[c]=x;
c++;
}
return c==2&&O[0]>=0&&O[1]>=0&&O[0]!=O[1];
}
bool link(int u,int v,int o0,int o1){
static vector<int>mk;
static int st=1;
if((int)mk.size()!=N)mk.assign(N,0);
if(++st==INT_MAX)fill(mk.begin(),mk.end(),0),st=1;
for(int id:inc[u])if(af[id]){
int a[3]={
Fv[id].a,Fv[id].b,Fv[id].c}
;
for(int x:a)if(x!=u&&x!=v)mk[x]=st;
}
int g0=0,g1=0,c=0;
for(int id:inc[v])if(af[id]){
int a[3]={
Fv[id].a,Fv[id].b,Fv[id].c}
;
for(int x:a)if(x!=u&&x!=v&&mk[x]==st){
if(x==o0)g0=1;
else if(x==o1)g1=1;
else return false;
c++;
mk[x]=0;
}
}
return g0&&g1&&c==2;
}
bool dup(int k,int r,int s0,int s1,int a,int b,int c){
auto key=fk(a,b,c);
vector<int>ids=inc[k];
ids.insert(ids.end(),inc[r].begin(),inc[r].end());
sort(ids.begin(),ids.end());
ids.erase(unique(ids.begin(),ids.end()),ids.end());
for(int id:ids)if(af[id]&&id!=s0&&id!=s1){
F g=Fv[id];
if(g.a==r)g.a=k;
if(g.b==r)g.b=k;
if(g.c==r)g.c=k;
if(g.a==g.b||g.a==g.c||g.b==g.c)continue;
if(fk(g.a,g.b,g.c)==key)return true;
}
return false;
}
bool ev(int k,int r,int s0,int s1,double md,double ar){
vector<int>ids=inc[r];
for(int id:ids)if(af[id]&&id!=s0&&id!=s1){
F old=Fv[id],nw=old;
if(nw.a==r)nw.a=k;
if(nw.b==r)nw.b=k;
if(nw.c==r)nw.c=k;
if(nw.a==nw.b||nw.a==nw.c||nw.b==nw.c)return false;
V co=crossp(P[old.b]-P[old.a],P[old.c]-P[old.a]),cn=crossp(P[nw.b]-P[nw.a],P[nw.c]-P[nw.a]);
double ao=norm(co),an=norm(cn);
if(ao<1e-20*diagv*diagv||an<1e-20*diagv*diagv)return false;
if(dotp(co,cn)<md*ao*an)return false;
if(an<ao*ar)return false;
if(dup(k,r,s0,s1,nw.a,nw.b,nw.c))return false;
}
return true;
}
bool col(int u,int v,double md,double ar){
if(!av[u]||!av[v]||u==v)return false;
int E[2],O[2];
if(!ef(u,v,E,O)||!link(u,v,O[0],O[1]))return false;
bool uv=ev(u,v,E[0],E[1],md,ar),vu=ev(v,u,E[0],E[1],md,ar);
if(!uv&&!vu)return false;
int k=u,r=v;
if(vu&&!uv)k=v,r=u;
if(af[E[0]])af[E[0]]=0,actf--;
if(af[E[1]])af[E[1]]=0,actf--;
for(int id:inc[r])if(af[id]&&hs(id,r)){
if(Fv[id].a==r)Fv[id].a=k;
if(Fv[id].b==r)Fv[id].b=k;
if(Fv[id].c==r)Fv[id].c=k;
if(Fv[id].a==Fv[id].b||Fv[id].a==Fv[id].c||Fv[id].b==Fv[id].c)af[id]=0,actf--;
else inc[k].push_back(id);
}
av[r]=0;
actv--;
vector<int>().swap(inc[r]);
return true;
}
vector<pair<double,uint64_t>> edges(double L){
unordered_set<uint64_t>seen;
seen.reserve(actf*3+1);
vector<pair<double,uint64_t>>e;
double L2=L*L;
for(int i=0;
i<M;
i++)if(af[i]){
int a[3]={
Fv[i].a,Fv[i].b,Fv[i].c}
;
for(int k=0;
k<3;
k++){
int u=a[k],v=a[(k+1)%3];
if(!av[u]||!av[v])continue;
uint64_t key=ek(u,v);
if(seen.insert(key).second){
double l2=n2(P[u]-P[v]);
if(l2<=L2)e.push_back({
l2,key}
);
}
}
}
sort(e.begin(),e.end());
return e;
}
Mesh snap(){
Mesh m;
vector<int>id(N,-1);
for(int i=0;
i<N;
i++)if(av[i]){
id[i]=m.v.size();
m.v.push_back(P[i]);
}
for(int i=0;
i<M;
i++)if(af[i]){
int a=Fv[i].a,b=Fv[i].b,c=Fv[i].c;
if(a!=b&&a!=c&&b!=c&&id[a]>=0&&id[b]>=0&&id[c]>=0)m.f.push_back({
id[a],id[b],id[c]}
);
}
return m;
}
vector<Mesh> run(){
vector<Mesh>sn;
struct Pp{
double L,md,ar;
}
;
vector<Pp>ph={
{
.010,.999,.35}
,{
.018,.996,.20}
,{
.030,.985,.10}
,{
.045,.955,.04}
,{
.065,.900,.015}
,{
.085,.780,.005}
}
;
for(auto p:ph){
if(tim()>13.2)break;
auto es=edges(p.L*diagv);
int done=0;
for(auto &x:es){
if((done&511)==0&&tim()>13.8)break;
int u=x.second>>32,v=x.second&0xffffffffu;
if(col(u,v,p.md,p.ar))done++;
}
rebuild();
Mesh s=snap();
if(valid(s))sn.push_back(move(s));
if(done<8)break;
}
return sn;
}
}
;


static bool tryCandidate(Mesh raw,const vector<V>&vn,const Render&origR,int R,double guard,Mesh&best,int&bestN){
if(raw.v.empty()||raw.f.empty()||raw.v.size()>=Orig.v.size())return false;
orientLike(raw);
if(!valid(raw))return false;
Mesh cand=raw;
int w=0;
if(!addWitnesses(cand,vn,&w))return false;
if(!outputToOriginalOK(cand))return false;
if((int)cand.v.size()>=bestN)return false;
Render rr;
render(cand,R,rr);
double sc=scoreR(origR,rr);
double need=guard;
if(w>0)need+=0.008;
if(sc>=need){
best=move(cand);
bestN=best.v.size();
return true;
}
return false;
}

int main(){
T0=chrono::steady_clock::now();
if(!read())return 0;
if(N<=20){
printM(Orig);
return 0;
}
vector<V>vn=vertexNormals();
int R=resFor();
double guard=guardFor();
Render origR;
render(Orig,R,origR);
Mesh best=Orig;
int bestN=N;
bool ok=false;

 {
Mesh b=buildBox(1);
if(tim()<16.0)ok|=tryCandidate(b,vn,origR,R,guard+0.015,best,bestN);
}
 
 {
Ell e=fitEll();
if(e.ok){
vector<pair<int,int>>tr={
{
8,16}
,{
10,20}
,{
12,24}
,{
16,32}
,{
20,40}
,{
24,48}
}
;
for(auto [la,lo]:tr)if(tim()<16.0)ok|=tryCandidate(buildEll(e,la,lo),vn,origR,R,guard,best,bestN);
}
}

 {
int r,s;
if(detectLat(r,s)){
vector<pair<int,int>>tr={
{
max(2,r/12),max(8,s/12)}
,{
max(2,r/8),max(8,s/8)}
,{
max(3,r/6),max(10,s/6)}
,{
max(4,r/4),max(12,s/4)}
,{
max(5,r/3),max(16,s/3)}
}
;
for(auto [a,b]:tr)if(tim()<16.2)ok|=tryCandidate(buildLat(r,s,a,b),vn,origR,R,guard-0.004,best,bestN);
}
}

 {
int u,s;
if(detectTor(u,s)){
vector<pair<int,int>>tr={
{
max(6,u/10),max(6,s/2)}
,{
max(6,u/8),max(6,s*2/3)}
,{
max(6,u/6),max(6,s*3/4)}
,{
max(6,u/4),max(6,s)}
}
;
for(auto [a,b]:tr)if(tim()<16.2)ok|=tryCandidate(buildTor(u,s,a,b),vn,origR,R,guard-0.004,best,bestN);
}
}

 if(tim()<15.0){
Simp sim;
auto snaps=sim.run();
sort(snaps.begin(),snaps.end(),[](const Mesh&a,const Mesh&b){
return a.v.size()<b.v.size();
}
);
int lim=0;
for(auto &s:snaps){
if(tim()>18.4||lim++>8)break;
ok|=tryCandidate(s,vn,origR,R,guard,best,bestN);
}
}

 if(!ok||!valid(best)||!outputToOriginalOK(best)||best.v.size()>Orig.v.size())best=Orig;
printM(best);
return 0;
}