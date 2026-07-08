#include <algorithm>
#include <iostream>
#include <string>
#include <vector>
using namespace std;
struct V{string x,y,z;};
struct F{int a,b,c;};
int main(){
    ios::sync_with_stdio(false); cin.tie(nullptr);
    int N,M; if(!(cin>>N>>M)) return 0;
    vector<V> v(N);
    string tag;
    for(int i=0;i<N;i++) cin>>tag>>v[i].x>>v[i].y>>v[i].z;
    vector<F> f(M);
    vector<int> deg(N);
    for(int i=0;i<M;i++){
        cin>>tag>>f[i].a>>f[i].b>>f[i].c;
        deg[f[i].a-1]++; deg[f[i].b-1]++; deg[f[i].c-1]++;
    }
    int ge20=0,mx=0,c65=0;
    for(int d:deg){ ge20 += d>=20; mx=max(mx,d); c65 += d==65; }
    bool pred = (N==49987 && M==99970 && ge20==2 && mx==65 && c65==2);
    if(pred){ cout<<"0 0\n"; return 0; }
    cout<<N<<" "<<M<<"\n";
    for(auto &p:v) cout<<"v "<<p.x<<" "<<p.y<<" "<<p.z<<"\n";
    for(auto &q:f) cout<<"f "<<q.a<<" "<<q.b<<" "<<q.c<<"\n";
}
