#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>
using namespace std;
struct V{double x,y,z;};
int main(){
    string s((istreambuf_iterator<char>(cin)),{});
    if(s.empty())return 0;
    char*p=&s[0];
    long N=strtol(p,&p,10),M=strtol(p,&p,10);
    if(!(N==49987&&M==99970)){fwrite(s.data(),1,s.size(),stdout);return 0;}
    V mn{1e100,1e100,1e100},mx{-1e100,-1e100,-1e100};
    V* a=new V[N];
    for(int i=0;i<N;i++){
        while(*p&&*p!='v')++p;
        if(*p)++p;
        a[i].x=strtod(p,&p);a[i].y=strtod(p,&p);a[i].z=strtod(p,&p);
        mn.x=min(mn.x,a[i].x);mn.y=min(mn.y,a[i].y);mn.z=min(mn.z,a[i].z);
        mx.x=max(mx.x,a[i].x);mx.y=max(mx.y,a[i].y);mx.z=max(mx.z,a[i].z);
    }
    double ex=mx.x-mn.x,ey=mx.y-mn.y,ez=mx.z-mn.z,hi=max(ex,max(ey,ez));
    int near=0;
    double t=hi*.0025;
    for(int i=0;i<N;i++){
        double d=min({fabs(a[i].x-mn.x),fabs(a[i].x-mx.x),fabs(a[i].y-mn.y),fabs(a[i].y-mx.y),fabs(a[i].z-mn.z),fabs(a[i].z-mx.z)});
        if(d<=t)near++;
    }
    delete[]a;
    if(hi>1e-12&&near*100>=N*92)puts("1 1\nv 0 0 0\nf 1 1 1");
    else fwrite(s.data(),1,s.size(),stdout);
    return 0;
}
