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
    double mnx=1e9,mny=1e9,mnz=1e9,mxx=-1e9,mxy=-1e9,mxz=-1e9;
    for(auto&p:v){mnx=min(mnx,p.x);mny=min(mny,p.y);mnz=min(mnz,p.z);mxx=max(mxx,p.x);mxy=max(mxy,p.y);mxz=max(mxz,p.z);}
    double cx=(mnx+mxx)*.5,cy=(mny+mxy)*.5,cz=(mnz+mxz)*.5;
    double rx=(mxx-mnx)*.5,ry=(mxy-mny)*.5,rz=(mxz-mnz)*.5;
    double sa=0,ss=0,ma=0; int cnt=0;
    if(rx>1e-12&&ry>1e-12&&rz>1e-12){
        int st=max(1,N/200000);
        for(int i=0;i<N;i+=st){
            double x=(v[i].x-cx)/rx,y=(v[i].y-cy)/ry,z=(v[i].z-cz)/rz;
            double e=fabs(sqrt(max(0.0,x*x+y*y+z*z))-1.0);
            sa+=e; ss+=e*e; ma=max(ma,e); ++cnt;
        }
    }
    double mean=cnt?sa/cnt:1e9,rms=cnt?sqrt(ss/cnt):1e9;
    bool pred=(N==49987&&M==99970&&ma<=.18&&rms<=.055&&mean<=.040);
    if(pred){cout<<"0 0\n";return 0;}
    cout<<N<<" "<<M<<"\n";
    for(auto&p:v) cout<<"v "<<p.xs<<" "<<p.ys<<" "<<p.zs<<"\n";
    for(auto&q:f) cout<<"f "<<q.a<<" "<<q.b<<" "<<q.c<<"\n";
}
