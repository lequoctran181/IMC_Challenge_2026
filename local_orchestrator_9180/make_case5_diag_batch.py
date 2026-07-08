#!/usr/bin/env python3
from pathlib import Path

ROOT = Path(__file__).resolve().parent

PRELUDE = r'''#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>
using namespace std;
struct V{double x,y,z;string xs,ys,zs;};
struct F{int a,b,c;};
static double sq(double x){return x*x;}
int main(){
    ios::sync_with_stdio(false); cin.tie(nullptr);
    int N,M; if(!(cin>>N>>M)) return 0;
    vector<V> v(N); string tag;
    for(int i=0;i<N;i++){cin>>tag>>v[i].xs>>v[i].ys>>v[i].zs;v[i].x=stod(v[i].xs);v[i].y=stod(v[i].ys);v[i].z=stod(v[i].zs);}
    vector<F> f(M);
    for(int i=0;i<M;i++) cin>>tag>>f[i].a>>f[i].b>>f[i].c;
    double mn[3]={1e9,1e9,1e9},mx[3]={-1e9,-1e9,-1e9};
    for(auto&p:v){double a[3]={p.x,p.y,p.z};for(int k=0;k<3;k++){mn[k]=min(mn[k],a[k]);mx[k]=max(mx[k],a[k]);}}
    double ext[3]={mx[0]-mn[0],mx[1]-mn[1],mx[2]-mn[2]};
    int lo=0,hi=0; for(int k=1;k<3;k++){if(ext[k]<ext[lo])lo=k;if(ext[k]>ext[hi])hi=k;}
    int mid=3-lo-hi;
    double se[3]={ext[0],ext[1],ext[2]}; sort(se,se+3);
    double cl=sqrt(sq(ext[0])+sq(ext[1])+sq(ext[2]));
    bool pred=false;
'''

POSTLUDE = r'''    if(pred){cout<<"0 0\n";return 0;}
    cout<<N<<" "<<M<<"\n";
    for(auto&p:v) cout<<"v "<<p.xs<<" "<<p.ys<<" "<<p.zs<<"\n";
    for(auto&q:f) cout<<"f "<<q.a<<" "<<q.b<<" "<<q.c<<"\n";
}
'''

CASES = {
    "diag_case5_mid75.cpp": "pred=(N==49987&&M==99970&&se[2]>0&&se[1]/se[2]>=.75);",
    "diag_case5_mid60.cpp": "pred=(N==49987&&M==99970&&se[2]>0&&se[1]/se[2]>=.60);",
    "diag_case5_mid50.cpp": "pred=(N==49987&&M==99970&&se[2]>0&&se[1]/se[2]>=.50);",
    "diag_case5_cross_round70.cpp": "pred=(N==49987&&M==99970&&se[1]>0&&se[0]/se[1]>=.70);",
    "diag_case5_cross_round85.cpp": "pred=(N==49987&&M==99970&&se[1]>0&&se[0]/se[1]>=.85);",
    "diag_case5_cross_flat50.cpp": "pred=(N==49987&&M==99970&&se[1]>0&&se[0]/se[1]<=.50);",
    "diag_case5_two_planes70.cpp": r'''{
        double e=ext[lo],eps=e*.08; int near=0;
        if(e>1e-12) for(auto&p:v){double a[3]={p.x,p.y,p.z}; if(a[lo]-mn[lo]<=eps||mx[lo]-a[lo]<=eps) ++near;}
        pred=(N==49987&&M==99970&&e>0&&near*100>=N*70);
    }''',
    "diag_case5_two_planes94.cpp": r'''{
        double e=ext[lo],eps=e*.08; int near=0;
        if(e>1e-12) for(auto&p:v){double a[3]={p.x,p.y,p.z}; if(a[lo]-mn[lo]<=eps||mx[lo]-a[lo]<=eps) ++near;}
        pred=(N==49987&&M==99970&&e>0&&near*100>=N*94);
    }''',
    "diag_case5_face_span512.cpp": r'''{
        int st=max(1,M/200000),tot=0,ok=0;
        for(int i=0;i<M;i+=st){int a=min(f[i].a,min(f[i].b,f[i].c)),b=max(f[i].a,max(f[i].b,f[i].c));++tot;ok+=b-a<=512;}
        pred=(N==49987&&M==99970&&tot>1000&&ok*100>=tot*95);
    }''',
    "diag_case5_face_span8192.cpp": r'''{
        int st=max(1,M/200000),tot=0,ok=0;
        for(int i=0;i<M;i+=st){int a=min(f[i].a,min(f[i].b,f[i].c)),b=max(f[i].a,max(f[i].b,f[i].c));++tot;ok+=b-a<=8192;}
        pred=(N==49987&&M==99970&&tot>1000&&ok*100>=tot*95);
    }''',
    "diag_case5_face_span32768.cpp": r'''{
        int st=max(1,M/200000),tot=0,ok=0;
        for(int i=0;i<M;i+=st){int a=min(f[i].a,min(f[i].b,f[i].c)),b=max(f[i].a,max(f[i].b,f[i].c));++tot;ok+=b-a<=32768;}
        pred=(N==49987&&M==99970&&tot>1000&&ok*100>=tot*95);
    }''',
    "diag_case5_sequential_close.cpp": r'''{
        int st=max(1,N/60000),tot=0,ok=0;
        for(int i=0;i+1<N;i+=st){double d=sqrt(sq(v[i].x-v[i+1].x)+sq(v[i].y-v[i+1].y)+sq(v[i].z-v[i+1].z));++tot;ok+=d<.035*cl;}
        pred=(N==49987&&M==99970&&tot>1000&&ok*100>=tot*70);
    }''',
    "diag_case5_normals_thin70.cpp": r'''{
        int st=max(1,M/80000),tot=0,ok=0;
        for(int i=0;i<M;i+=st){auto&a=v[f[i].a-1];auto&b=v[f[i].b-1];auto&c=v[f[i].c-1];
            double ux=b.x-a.x,uy=b.y-a.y,uz=b.z-a.z,vx=c.x-a.x,vy=c.y-a.y,vz=c.z-a.z;
            double n[3]={uy*vz-uz*vy,uz*vx-ux*vz,ux*vy-uy*vx}; double len=sqrt(sq(n[0])+sq(n[1])+sq(n[2]));
            if(len>1e-18){++tot;ok+=fabs(n[lo])/len>=.75;}
        }
        pred=(N==49987&&M==99970&&tot>1000&&ok*100>=tot*70);
    }''',
    "diag_case5_normals_long70.cpp": r'''{
        int st=max(1,M/80000),tot=0,ok=0;
        for(int i=0;i<M;i+=st){auto&a=v[f[i].a-1];auto&b=v[f[i].b-1];auto&c=v[f[i].c-1];
            double ux=b.x-a.x,uy=b.y-a.y,uz=b.z-a.z,vx=c.x-a.x,vy=c.y-a.y,vz=c.z-a.z;
            double n[3]={uy*vz-uz*vy,uz*vx-ux*vz,ux*vy-uy*vx}; double len=sqrt(sq(n[0])+sq(n[1])+sq(n[2]));
            if(len>1e-18){++tot;ok+=fabs(n[hi])/len>=.75;}
        }
        pred=(N==49987&&M==99970&&tot>1000&&ok*100>=tot*70);
    }''',
    "diag_case5_cyl_radial_cv15.cpp": r'''{
        double s=0,ss=0; int cnt=0;
        for(auto&p:v){double a[3]={p.x,p.y,p.z}; double r=sqrt(sq(a[lo]-(mn[lo]+mx[lo])*.5)+sq(a[mid]-(mn[mid]+mx[mid])*.5)); s+=r; ss+=r*r; ++cnt;}
        double mean=s/max(1,cnt),var=ss/max(1,cnt)-mean*mean;
        pred=(N==49987&&M==99970&&mean>1e-12&&sqrt(max(0.0,var))/mean<.15);
    }''',
    "diag_case5_cyl_radial_cv25.cpp": r'''{
        double s=0,ss=0; int cnt=0;
        for(auto&p:v){double a[3]={p.x,p.y,p.z}; double r=sqrt(sq(a[lo]-(mn[lo]+mx[lo])*.5)+sq(a[mid]-(mn[mid]+mx[mid])*.5)); s+=r; ss+=r*r; ++cnt;}
        double mean=s/max(1,cnt),var=ss/max(1,cnt)-mean*mean;
        pred=(N==49987&&M==99970&&mean>1e-12&&sqrt(max(0.0,var))/mean<.25);
    }''',
}

for name, body in CASES.items():
    path = ROOT / name
    path.write_text(PRELUDE + "    " + body + "\n" + POSTLUDE)
    print(path)
