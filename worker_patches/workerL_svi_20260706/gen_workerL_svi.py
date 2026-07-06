#!/usr/bin/env python3
import re,sys,pathlib
if len(sys.argv)!=3:
    print("usage: python3 gen_workerL_svi.py fetched_sources/19901232.cpp workerL_svi.cpp",file=sys.stderr);sys.exit(2)
p=pathlib.Path(sys.argv[1]);s=p.read_text()
ins=r'''namespace SVI{struct C{double z;int i;};static Face mf(int a,int b,int c){Face f;f.v[0]=a;f.v[1]=b;f.v[2]=c;return f;}static void add(vector<Vec3>&X,vector<Face>&F,Vec3 c,Vec3 n,double s){double l=norm3(n);if(!(l>1e-12))return;n=n*(1./l);Vec3 a=fabs(n.x)<.7?Vec3{1,0,0}:Vec3{0,1,0},u=cross3(n,a);l=norm3(u);if(!(l>1e-12))return;u=u*(s/l);Vec3 v=cross3(n,u);l=norm3(v);if(!(l>1e-12))return;v=v*(s/l);Vec3 w=n*(s*.28);int o=X.size();X.push_back(c+w);X.push_back(c+u);X.push_back(c+v);X.push_back(c-u);X.push_back(c-v);X.push_back(c-w);int T[8][3]={{0,1,2},{0,2,3},{0,3,4},{0,4,1},{5,2,1},{5,3,2},{5,4,3},{5,1,4}};for(auto&t:T)F.push_back(mf(o+t[0],o+t[1],o+t[2]));}static bool run(){if(N<30000||es()>14.3)return 0;AP S=AD();int base=cove();if(base<200)return 0;vector<Vec3>NN(N,{0,0,0});for(int i=0;i<M;i++){Face f=AR[i];Vec3 cr=cross3(originalP[f.v[1]]-originalP[f.v[0]],originalP[f.v[2]]-originalP[f.v[0]]);NN[f.v[0]]=NN[f.v[0]]+cr;NN[f.v[1]]=NN[f.v[1]]+cr;NN[f.v[2]]=NN[f.v[2]]+cr;if((i&262143)==0&&es()>15.1){rs(S);return 0;}}int GG[3]={20,26,32};for(int gi=0;gi<3&&es()<16.4;gi++){int G=GG[gi];vector<C>B(6*G*G,{1e100,-1});for(int i=0;i<N;i++){Vec3 p=originalP[i];for(int w=0;w<6;w++){double sx,sy,z;if(w==0){sx=p.y;sy=p.z;z=2.5-p.x;}else if(w==1){sx=-p.y;sy=p.z;z=2.5+p.x;}else if(w==2){sx=-p.x;sy=p.z;z=2.5-p.y;}else if(w==3){sx=p.x;sy=p.z;z=2.5+p.y;}else if(w==4){sx=p.x;sy=p.y;z=2.5-p.z;}else{sx=-p.x;sy=p.y;z=2.5+p.z;}if(z<=0)continue;int x=(int)((800.*sx/z+512.)*G/1024.),y=(int)((800.*sy/z+512.)*G/1024.);if(x<0||x>=G||y<0||y>=G)continue;C&b=B[(w*G+y)*G+x];if(z<b.z){b.z=z;b.i=i;}}if((i&262143)==0&&es()>16.0)break;}vector<int>mk(N,0),id;for(auto&b:B)if(b.i>=0&&!mk[b.i]){mk[b.i]=1;id.push_back(b.i);}if((int)id.size()*6>=base*97/100)continue;vector<Vec3>X;vector<Face>F;X.reserve(id.size()*6);F.reserve(id.size()*8);double ss=CL/(G*2.15);for(int q:id){Vec3 n=NN[q];double l=norm3(n);if(!(l>1e-12))n=originalP[q];add(X,F,originalP[q],n,ss);}if((int)X.size()>=base||X.empty())continue;if(AF(X,F)&&W5::strong_validator()&&cove()<base&&es()<18.4&&vps(256)>.944)return 1;rs(S);}return 0;}}'''
a='static void CA(string&out,const char*line,int len){'
if s.count(a)!=1:
    print("fail-closed: insertion anchor missing",file=sys.stderr);sys.exit(1)
s=s.replace(a,ins+a)
old='MIDEC::run();WK::run();'
new='MIDEC::run();if(SVI::run()){JD();return 0;}WK::run();'
if s.count(old)!=1:
    print("fail-closed: main anchor missing",file=sys.stderr);sys.exit(1)
s=s.replace(old,new)
macro_after='#define svs safe_vertices\n'
macros='#define S0 size\n#define B0 begin\n#define E0 end\n#define VC vector\n#define _S static\n#define _D double\n#define _I int\ntypedef unsigned long long ULL;\n'
if macro_after not in s:
    print("fail-closed: macro anchor missing",file=sys.stderr);sys.exit(1)
s=s.replace(macro_after,macro_after+macros,1)
s=s.replace('.push_back(','.pb(').replace('.size()','.S0()').replace('.begin()','.B0()').replace('.end()','.E0()')
s=s.replace('vector<','VC<').replace('unsigned long long','ULL')
s=re.sub(r'\bstatic\b','_S',s)
s=re.sub(r'\bdouble\b','_D',s)
s=re.sub(r'\bint\b','_I',s)
# Do not allow macro definitions to self-corrupt.
s=s.replace('#define _S _S','#define _S static').replace('#define _D _D','#define _D double').replace('#define _I _I','#define _I int')
s=s.replace('typedef ULL ULL;','typedef unsigned long long ULL;')
if 'SVI::run()' not in s or 'int main' in s:
    pass
n=len(s.encode())
if n>131072:
    print("fail-closed: emitted source too large",n,file=sys.stderr);sys.exit(1)
pathlib.Path(sys.argv[2]).write_text(s)
print(n)
