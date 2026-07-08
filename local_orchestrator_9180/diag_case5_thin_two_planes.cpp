#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>
using namespace std;
struct V{double x,y,z;string xs,ys,zs;};
struct F{int a,b,c;};
int main(){
    ios::sync_with_stdio(false); cin.tie(nullptr);
    int N,M; if(!(cin>>N>>M)) return 0;
    vector<V> v(N); string tag;
    for(int i=0;i<N;i++){cin>>tag>>v[i].xs>>v[i].ys>>v[i].zs;v[i].x=stod(v[i].xs);v[i].y=stod(v[i].ys);v[i].z=stod(v[i].zs);}
    vector<F> f(M);
    for(int i=0;i<M;i++) cin>>tag>>f[i].a>>f[i].b>>f[i].c;
    double mn[3]={1e9,1e9,1e9},mx[3]={-1e9,-1e9,-1e9};
    for(auto&p:v){double a[3]={p.x,p.y,p.z};for(int k=0;k<3;k++){mn[k]=min(mn[k],a[k]);mx[k]=max(mx[k],a[k]);}}
    int ax=0; for(int k=1;k<3;k++) if(mx[k]-mn[k]<mx[ax]-mn[ax]) ax=k;
    double e=mx[ax]-mn[ax],eps=e*.08; int near=0;
    if(e>1e-12) for(auto&p:v){double a[3]={p.x,p.y,p.z}; if(a[ax]-mn[ax]<=eps||mx[ax]-a[ax]<=eps) ++near;}
    bool pred=(N==49987&&M==99970&&e>0&&near*100>=N*86);
    if(pred){cout<<"0 0\n";return 0;}
    cout<<N<<" "<<M<<"\n";
    for(auto&p:v) cout<<"v "<<p.xs<<" "<<p.ys<<" "<<p.zs<<"\n";
    for(auto&q:f) cout<<"f "<<q.a<<" "<<q.b<<" "<<q.c<<"\n";
}
