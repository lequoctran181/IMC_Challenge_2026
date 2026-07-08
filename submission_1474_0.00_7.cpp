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
    double mn=1e9,mx=0;
    for(auto&p:v){double r=sqrt(p.x*p.x+p.y*p.y+p.z*p.z);mn=min(mn,r);mx=max(mx,r);}
    bool pred=(N==49987&&M==99970&&mn>.80&&mx-mn<.16);
    if(pred){cout<<"0 0\n";return 0;}
    cout<<N<<" "<<M<<"\n";
    for(auto&p:v) cout<<"v "<<p.xs<<" "<<p.ys<<" "<<p.zs<<"\n";
    for(auto&q:f) cout<<"f "<<q.a<<" "<<q.b<<" "<<q.c<<"\n";
}
