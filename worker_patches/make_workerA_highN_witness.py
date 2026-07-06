#!/usr/bin/env python3
import sys,re,collections,string

LIM=131072
if len(sys.argv)<3:
    print("usage: python3 make_workerA_highN_witness.py fetched_sources/19901232.cpp workerA_highN_witness.cpp",file=sys.stderr)
    sys.exit(2)

src=open(sys.argv[1],"r",encoding="utf-8").read()
orig=len(src.encode())

HN=r'''namespace HNWX{struct T{double a,b,u,p,c;int g,r;};static int V(vector<Vec3>&X,vector<Face>&A){vector<int>m(N,-1);X.clear();A.clear();X.reserve(cove());A.reserve(BE);for(int i=0;i<(int)faces.size();++i)if(BR[i]){Face f=faces[i],g;bool o=1;for(int k=0;k<3;++k){int v=f.v[k];if(v<0||v>=N||!BU[v])o=0;else{if(m[v]<0){m[v]=X.size();X.pb(P[v]);}g.v[k]=m[v];}}if(o&&g.v[0]!=g.v[1]&&g.v[1]!=g.v[2]&&g.v[0]!=g.v[2])A.pb(g);}return X.size();}static int C(double A,double B,double U,int G){AE p;p.AI=A*CL;p.BB=B*CL;p.BQ=cos(U*acos(-1.)/180.);p.W=max(1e-11,1e-9*CL);p.AT=1.-1e-10;p.AJ=max(1e-30,1e-24*CL*CL);p.AA=false;int R=0;for(int q=0;q<G&&es()<17.4;++q)for(int i=0;i<(int)faces.size()&&es()<17.4;++i){if(!BR[i])continue;Face f=faces[i];int a=f.v[0],b=f.v[1],c=f.v[2];if(!BU[a]||!BU[b]||!BU[c])continue;double x=norm2(P[a]-P[b]),y=norm2(P[b]-P[c]),z=norm2(P[c]-P[a]);if(y<x&&y<z){if(GD(b,c,p))++R;}else if(z<x&&z<y){if(GD(c,a,p))++R;}else if(GD(a,b,p))++R;}return R;}static void E(vector<Vec3>&X,vector<Face>&A,const Vec3&p){double e=CL*8e-5;int s=X.size();X.pb(p);X.pb({p.x+e,p.y,p.z});X.pb({p.x,p.y+e,p.z});X.pb({p.x,p.y,p.z+e});int q[4][3]={{0,2,1},{0,1,3},{0,3,2},{1,2,3}};for(int i=0;i<4;++i){Face f;f.v[0]=s+q[i][0];f.v[1]=s+q[i][1];f.v[2]=s+q[i][2];A.pb(f);}}static bool W(vector<Vec3>&X,vector<Face>&A,int S,double Cc){Vec3 mn=originalP[0],mx=originalP[0];for(const Vec3&p:originalP){mn.x=min(mn.x,p.x);mn.y=min(mn.y,p.y);mn.z=min(mn.z,p.z);mx.x=max(mx.x,p.x);mx.y=max(mx.y,p.y);mx.z=max(mx.z,p.z);}double t=.049*CL,t2=t*t;int nx=(int)((mx.x-mn.x)/t)+3,ny=(int)((mx.y-mn.y)/t)+3,nz=(int)((mx.z-mn.z)/t)+3;vector<vector<int>>G(nx*ny*nz);auto I=[&](const Vec3&p,int&x,int&y,int&z){x=(int)((p.x-mn.x)/t);y=(int)((p.y-mn.y)/t);z=(int)((p.z-mn.z)/t);if(x<0)x=0;if(y<0)y=0;if(z<0)z=0;if(x>=nx)x=nx-1;if(y>=ny)y=ny-1;if(z>=nz)z=nz-1;};auto K=[&](int x,int y,int z){return(x*ny+y)*nz+z;};auto A0=[&](int j){int x,y,z;I(X[j],x,y,z);G[K(x,y,z)].pb(j);};for(int i=0;i<(int)X.size();++i)A0(i);auto N0=[&](const Vec3&p){int x,y,z;I(p,x,y,z);for(int a=max(0,x-1);a<=min(nx-1,x+1);++a)for(int b=max(0,y-1);b<=min(ny-1,y+1);++b)for(int c=max(0,z-1);c<=min(nz-1,z+1);++c)for(int j:G[K(a,b,c)])if(norm2(X[j]-p)<=t2)return true;return false;};for(int i=0;i<N&&es()<18.1;++i)if(!N0(originalP[i])){int s=X.size();E(X,A,originalP[i]);A0(s);if((double)X.size()*100.>S*Cc||X.size()>(size_t)N)return false;}return X.size()<(size_t)S&&(double)X.size()*100.<=S*Cc;}static bool Q(const AP&S,int SN,const T&t){rs(S);C(t.a,t.b,t.u,t.g);vector<Vec3>X;vector<Face>A;int v=V(X,A);if(v<=0||v>=SN){rs(S);return false;}if(!W(X,A,SN,t.c)){rs(S);return false;}bool ok=false;if(AF(X,A)&&W5::strong_validator()&&es()<18.5)ok=cove()<SN&&vps(t.r)>=t.p;if(!ok)rs(S);return ok;}static bool run(){if(N<200000||N>1100000||!AS||BG>.02||AL<.6||AH>.22||es()>16.9)return false;int S=cove();if(S<2000)return false;AP B=AD();T t[]={{.18,.04,30.,.935,75.,1,64},{.28,.055,42.,.925,60.,1,64},{.12,.032,24.,.945,85.,2,96}};for(auto&x:t)if(Q(B,S,x))return true;rs(B);return false;}}'''

old_main='int main(){JC();GN();if(!W2G::run())W2C::run();W5::post_patch_pass();VIMP::run();MIDEC::run();WK::run();B16::R(39000,60000,220,-7,192,.96,18.05);if(0&&"B16P515A"){}B16::R(39000,60000,76,-10,192,.96,18.35);for(int i=2;i--&&N>47500&N<6e4&&es()<18.6;)WK::run();if(N<50625&&es()<18.9)WK::run();JD();}'
new_main='int main(){JC();GN();if(!HNWX::run()){if(!W2G::run())W2C::run();W5::post_patch_pass();VIMP::run();MIDEC::run();WK::run();for(int i=2;i--&&N>47500&N<6e4&&es()<18.6;)WK::run();if(N<50625&&es()<18.9)WK::run();}JD();}'

src,n=re.subn(r'namespace B16\{.*?\}\}(?=int main\(\))','',src,1,flags=re.S)
if n!=1:
    raise SystemExit("B16 namespace anchor not found exactly once; aborting")

if old_main not in src:
    raise SystemExit("current-best main anchor not found; aborting")
src=src.replace(old_main,new_main,1)

anchor='static void CA(string&out'
if anchor not in src:
    raise SystemExit("CA anchor not found; aborting")
src=src.replace(anchor,HN+anchor,1)

kw=set('alignas alignof and and_eq asm auto bitand bitor bool break case catch char char16_t char32_t class compl const constexpr const_cast continue decltype default delete do double dynamic_cast else enum explicit export extern false float for friend goto if inline int long mutable namespace new noexcept not not_eq nullptr operator or or_eq private protected public register reinterpret_cast return short signed sizeof static static_assert struct switch template this thread_local throw true try typedef typeid typename union unsigned using virtual void volatile wchar_t while xor xor_eq'.split())
std=set('algorithm array chrono utility cstdint queue cmath cstdio cstdlib cstring string vector std size_t uint64_t uint32_t int64_t int32_t FILE stdin stdout stderr sort min max swap fill unique adjacent_find lower_bound find abs ceil floor sqrt cos sin acos fabs hypot pow isfinite strtod strtol fread fwrite printf snprintf setvbuf memset reserve resize assign push_back insert erase end begin data clear empty size front back top pop push shrink_to_fit priority_queue less greater move pair make_pair steady_clock time_point duration count now remove_cv remove_reference type _IOFBF'.split())
protect=kw|std|{'main','include','define','first','second','append'}

tok=re.compile(r'"(?:\\.|[^"\\])*"|\'(?:\\.|[^\'\\])*\'|(?<![A-Za-z0-9_])[A-Za-z_]\w*',re.S)
ids=[]
for m in tok.finditer(src):
    s=m.group(0)
    if s[0] not in '"\'':
        ids.append(s)
cnt=collections.Counter(ids)
used=set(ids)|protect

alphabet=string.ascii_letters
tail=string.ascii_letters+string.digits
def names():
    for a in alphabet:
        yield a
    for a in alphabet:
        for b in tail:
            yield a+b
    for a in alphabet:
        for b in tail:
            for c in tail:
                yield a+b+c

def good_target(x):
    return x not in used and x not in protect and not re.match(r'^[0-9]',x)

cands=[x for x,c in cnt.items() if len(x)>=3 and x not in protect and not x.startswith('__')]
cands.sort(key=lambda x:cnt[x]*(len(x)-1),reverse=True)

mp={}
ng=names()
for x in cands:
    while True:
        y=next(ng)
        if good_target(y):
            break
    if len(y)<len(x):
        mp[x]=y
        used.add(y)

def repl(m):
    s=m.group(0)
    if s[0] in '"\'':
        return s
    return mp.get(s,s)

out=tok.sub(repl,src)
final=len(out.encode())

if final>LIM:
    raise SystemExit(f"generated source over limit: {final}>{LIM}; minifier saved {orig-final} bytes before limit check")

open(sys.argv[2],"w",encoding="utf-8",newline="").write(out)
print(f"orig={orig} final={final} net_delta={final-orig} saved_vs_orig={orig-final} limit={LIM}")
print("wrote",sys.argv[2])
