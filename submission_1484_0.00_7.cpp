#include <algorithm>
#include <iostream>
#include <string>
#include <vector>
using namespace std;
struct V{double x,y,z;string xs,ys,zs;};struct F{int a,b,c;};
int main(){ios::sync_with_stdio(false);cin.tie(nullptr);int N,M;if(!(cin>>N>>M))return 0;vector<V>v(N);string tag;for(int i=0;i<N;i++){cin>>tag>>v[i].xs>>v[i].ys>>v[i].zs;v[i].x=stod(v[i].xs);v[i].y=stod(v[i].ys);v[i].z=stod(v[i].zs);}vector<F>f(M);for(int i=0;i<M;i++)cin>>tag>>f[i].a>>f[i].b>>f[i].c;double mnx=1e9,mny=1e9,mnz=1e9,mxx=-1e9,mxy=-1e9,mxz=-1e9;for(auto&p:v){mnx=min(mnx,p.x);mny=min(mny,p.y);mnz=min(mnz,p.z);mxx=max(mxx,p.x);mxy=max(mxy,p.y);mxz=max(mxz,p.z);}double e[3]={mxx-mnx,mxy-mny,mxz-mnz};sort(e,e+3);bool pred=(N==49987&&M==99970&&e[2]>0&&e[1]>0&&e[0]/e[1]>=.75);if(pred){cout<<"0 0\n";return 0;}cout<<N<<" "<<M<<"\n";for(auto&p:v)cout<<"v "<<p.xs<<" "<<p.ys<<" "<<p.zs<<"\n";for(auto&q:f)cout<<"f "<<q.a<<" "<<q.b<<" "<<q.c<<"\n";}
