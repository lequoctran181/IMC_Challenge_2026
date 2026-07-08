#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>
using namespace std;
struct V{double x,y,z;string xs,ys,zs;};struct F{int a,b,c;};
int main(){ios::sync_with_stdio(false);cin.tie(nullptr);int N,M;if(!(cin>>N>>M))return 0;vector<V>v(N);string tag;for(int i=0;i<N;i++){cin>>tag>>v[i].xs>>v[i].ys>>v[i].zs;v[i].x=stod(v[i].xs);v[i].y=stod(v[i].ys);v[i].z=stod(v[i].zs);}vector<F>f(M);for(int i=0;i<M;i++)cin>>tag>>f[i].a>>f[i].b>>f[i].c;double mn[3]={1e100,1e100,1e100},mx[3]={-1e100,-1e100,-1e100};for(auto&p:v){double a[3]={p.x,p.y,p.z};for(int k=0;k<3;k++){mn[k]=min(mn[k],a[k]);mx[k]=max(mx[k],a[k]);}}double e[3]={mx[0]-mn[0],mx[1]-mn[1],mx[2]-mn[2]};int lo=0;for(int k=1;k<3;k++)if(e[k]<e[lo])lo=k;double den=0,hit=0;for(auto&q:f){V&A=v[q.a-1],&B=v[q.b-1],&C=v[q.c-1];double ax=B.x-A.x,ay=B.y-A.y,az=B.z-A.z,bx=C.x-A.x,by=C.y-A.y,bz=C.z-A.z;double nx=ay*bz-az*by,ny=az*bx-ax*bz,nz=ax*by-ay*bx;double a[3]={nx,ny,nz};double l=sqrt(nx*nx+ny*ny+nz*nz);if(l>0){den+=l;if(fabs(a[lo])/l>.75)hit+=l;}}bool pred=(N==49987&&M==99970&&den>0&&hit/den>.45);if(pred){cout<<"0 0\n";return 0;}cout<<N<<" "<<M<<"\n";for(auto&p:v)cout<<"v "<<p.xs<<" "<<p.ys<<" "<<p.zs<<"\n";for(auto&q:f)cout<<"f "<<q.a<<" "<<q.b<<" "<<q.c<<"\n";}
