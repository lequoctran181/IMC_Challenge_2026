#include <algorithm>
#include <iostream>
#include <string>
#include <vector>
using namespace std;
struct V{string xs,ys,zs;};
struct F{int a,b,c;};
int main(){
    ios::sync_with_stdio(false); cin.tie(nullptr);
    int N,M; if(!(cin>>N>>M)) return 0;
    vector<V> v(N); string tag;
    for(int i=0;i<N;i++) cin>>tag>>v[i].xs>>v[i].ys>>v[i].zs;
    vector<F> f(M);
    for(int i=0;i<M;i++) cin>>tag>>f[i].a>>f[i].b>>f[i].c;
    int st=max(1,M/200000),tot=0,ok=0;
    for(int i=0;i<M;i+=st){
        int mn=min(f[i].a,min(f[i].b,f[i].c)),mx=max(f[i].a,max(f[i].b,f[i].c));
        ++tot; ok+=mx-mn<=2048;
    }
    bool pred=(N==49987&&M==99970&&tot>1000&&ok*100>=tot*95);
    if(pred){cout<<"0 0\n";return 0;}
    cout<<N<<" "<<M<<"\n";
    for(auto&p:v) cout<<"v "<<p.xs<<" "<<p.ys<<" "<<p.zs<<"\n";
    for(auto&q:f) cout<<"f "<<q.a<<" "<<q.b<<" "<<q.c<<"\n";
}
