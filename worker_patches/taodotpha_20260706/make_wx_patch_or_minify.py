#!/usr/bin/env python3
import sys,re
L=131072
s=open(sys.argv[1],"r",encoding="utf-8").read() if len(sys.argv)>1 else sys.stdin.read()
def be(t,p):
    b=t.find("{",p)
    if b<0:return-1
    d=0;i=b;n=len(t)
    while i<n:
        c=t[i]
        if c in"'\"":
            q=c;i+=1
            while i<n:
                if t[i]=="\\":i+=2;continue
                if t[i]==q:i+=1;break
                i+=1
            continue
        if c=="/"and i+1<n and t[i+1]=="/":
            j=t.find("\n",i+2);i=n if j<0 else j+1;continue
        if c=="/"and i+1<n and t[i+1]=="*":
            j=t.find("*/",i+2);i=n if j<0 else j+2;continue
        if c=="{":d+=1
        elif c=="}":
            d-=1
            if d==0:return i+1
        i+=1
    return-1
def ns(t,x):
    p=t.find("namespace "+x+"{")
    if p<0:return None
    q=be(t,p)
    return None if q<0 else(p,q)
def wx(usew5):
    pd='' if usew5 else r'''static double pd(const Vec3&p,const Vec3&a,const Vec3&b,const Vec3&c){Vec3 ab=b-a,ac=c-a,ap=p-a;double d1=dot3(ab,ap),d2=dot3(ac,ap);if(d1<=0&&d2<=0)return norm2(ap);Vec3 bp=p-b;double d3=dot3(ab,bp),d4=dot3(ac,bp);if(d3>=0&&d4<=d3)return norm2(bp);double vc=d1*d4-d3*d2;if(vc<=0&&d1>=0&&d3<=0){double v=d1/(d1-d3);return norm2(p-(a+ab*v));}Vec3 cp=p-c;double d5=dot3(ab,cp),d6=dot3(ac,cp);if(d6>=0&&d5<=d6)return norm2(cp);double vb=d5*d2-d1*d6;if(vb<=0&&d2>=0&&d6<=0){double w=d2/(d2-d6);return norm2(p-(a+ac*w));}double va=d3*d6-d5*d4;if(va<=0&&(d4-d3)>=0&&(d5-d6)>=0){double w=(d4-d3)/((d4-d3)+(d5-d6));return norm2(p-(b+(c-b)*w));}Vec3 n=cross3(ab,ac);double nn=norm2(n);if(!(nn>0))return 1e100;double x=dot3(p-a,n);return x*x/nn;}'''
    dist='W5::point_tri_d2' if usew5 else 'pd'
    return r'''namespace WX{struct C{int v;vector<int>i,r;vector<Face>f;double s;};static bool h(const vector<int>&v,int x){return find(v.begin(),v.end(),x)!=v.end();}static bool hk(const vector<unsigned long long>&v,unsigned long long x){return binary_search(v.begin(),v.end(),x);}'''+pd+r'''static int ec(int a,int b,const vector<int>&sk){int c=0;for(int fid:Y[a])if(BR[fid]&&!h(sk,fid)&&AC(fid,a)&&AC(fid,b)&&++c>2)return c;return c;}static bool fx(const Face&g,const vector<int>&sk){int a=g.v[0],b=g.v[1],c=g.v[2],x=a;if(Y[b].size()<Y[x].size())x=b;if(Y[c].size()<Y[x].size())x=c;for(int fid:Y[x])if(BR[fid]&&!h(sk,fid)&&IA(fid,a,b,c))return true;return false;}static int id(const vector<int>&v,int x){for(int i=0;i<(int)v.size();++i)if(v[i]==x)return i;return-1;}static bool fan(C&o,int rt,const Vec3&n,double pl,vector<Face>&of,double&bs){int k=o.r.size();vector<Face>nf;nf.reserve(k-2);double ca=cos(48.*acos(-1.)/180.),mi=max(1e-30,1e-24*CL*CL);for(int t=1;t+1<k;++t){Face g;g.v[0]=o.r[rt];g.v[1]=o.r[(rt+t)%k];g.v[2]=o.r[(rt+t+1)%k];if(g.v[0]==g.v[1]||g.v[0]==g.v[2]||g.v[1]==g.v[2])return false;Vec3 cr=cross3(P[g.v[1]]-P[g.v[0]],P[g.v[2]]-P[g.v[0]]);if(dot3(cr,n)<0){swap(g.v[1],g.v[2]);cr=cr*(-1.);}double a2=norm2(cr),l=sqrt(max(0.,a2));if(!(a2>mi)||dot3(cr*(1./l),n)<ca||fx(g,o.i))return false;nf.pb(g);}vector<unsigned long long>bd,ed;bd.reserve(k);ed.reserve(3*nf.size());for(int t=0;t<k;++t)bd.pb(ED(o.r[t],o.r[(t+1)%k]));sort(bd.begin(),bd.end());for(const Face&g:nf){ed.pb(ED(g.v[0],g.v[1]));ed.pb(ED(g.v[1],g.v[2]));ed.pb(ED(g.v[2],g.v[0]));}sort(ed.begin(),ed.end());int bc=0;for(size_t a=0;a<ed.size();){size_t b=a+1;while(b<ed.size()&&ed[b]==ed[a])++b;int cn=b-a,x=ed[a]>>32,y=ed[a]&0xffffffffu,out=ec(x,y,o.i);bool is=hk(bd,ed[a]);if(is){if(cn!=1||out!=1)return false;++bc;}else if(cn!=2||out!=0)return false;a=b;}if(bc!=k)return false;double d2=1e100;for(const Face&g:nf)d2=min(d2,'''+dist+r'''(P[o.v],P[g.v[0]],P[g.v[1]],P[g.v[2]]));if(!(d2<1e90))return false;double d=sqrt(max(0.,d2)),lim=(N<30000?.0485:.0475)*CL;if(BF[o.v]+d>lim)return false;double sc=(BF[o.v]+d)/CL+.12*pl/CL+.002*k;if(sc<bs){bs=sc;of=nf;return true;}return false;}static bool ck(int v,C&o){o=C();o.v=v;o.s=1e100;if(v<0||v>=N||!BU[v])return false;for(int fid:Y[v])if(BR[fid]&&AC(fid,v))o.i.pb(fid);int k=o.i.size();if(k<5||k>7)return false;vector<int>nb;vector<pair<int,int>>ee;Vec3 ns{0,0,0};for(int fid:o.i){const Face&f=faces[fid];int q[2],t=0;for(int j=0;j<3;++j)if(f.v[j]!=v)q[t++]=f.v[j];if(t!=2||q[0]==q[1]||!BU[q[0]]||!BU[q[1]])return false;if(!h(nb,q[0]))nb.pb(q[0]);if(!h(nb,q[1]))nb.pb(q[1]);ee.pb({q[0],q[1]});Vec3 cr=II(fid);double l=norm3(cr);if(!(l>0))return false;ns=ns+cr*(1./l);}if((int)nb.size()!=k)return false;double nl=norm3(ns);if(!(nl>1e-12))return false;Vec3 n=ns*(1./nl);for(int fid:o.i)if(dot3(GQ(fid),n)<.58)return false;vector<vector<int>>ad(k);for(auto e:ee){int a=id(nb,e.first),b=id(nb,e.second);if(a<0||b<0||a==b||h(ad[a],e.second)||h(ad[b],e.first))return false;ad[a].pb(e.second);ad[b].pb(e.first);}for(int a=0;a<k;++a)if(ad[a].size()!=2)return false;int st=nb[0],pr=-1,cu=st;for(int step=0;step<k;++step){o.r.pb(cu);int x=id(nb,cu),nx=ad[x][0]==pr?ad[x][1]:ad[x][0];pr=cu;cu=nx;}if(cu!=st)return false;vector<int>zz=o.r;sort(zz.begin(),zz.end());if(adjacent_find(zz.begin(),zz.end())!=zz.end())return false;Vec3 ct=P[v];for(int x:o.r)ct=ct+P[x];ct=ct*(1./(k+1));double pl=fabs(dot3(P[v]-ct,n)),dm=0;for(int x:o.r){pl=max(pl,fabs(dot3(P[x]-ct,n)));dm=max(dm,norm3(P[x]-P[v]));}if(pl>(N<30000?.020:.016)*CL||dm>(N<30000?.150:.120)*CL)return false;vector<Face>bf;double bs=1e100;for(int r=0;r<k;++r)fan(o,r,n,pl,bf,bs);if(bf.empty())return false;o.f=bf;o.s=bs;return true;}static bool ap(const C&c){if(!BU[c.v])return false;C d;if(!ck(c.v,d))return false;for(int fid:d.i)if(BR[fid]){BR[fid]=0;--BE;}BU[d.v]=0;Y[d.v].clear();for(const Face&g:d.f){int z=faces.size();faces.pb(g);BR.pb(1);++BE;Y[g.v[0]].pb(z);Y[g.v[1]].pb(z);Y[g.v[2]].pb(z);}return true;}static bool run(){static int once=0;if(once++||N<8000||N>47500||es()>17.65)return false;if(AS&&(BG>.024||AL<.48||AH>.32))return false;AP S=AD();int SN=cove();if(SN<=0)return false;vector<C>cs;cs.reserve(min(N,6000));for(int v=0;v<N&&es()<18.10;++v)if(BU[v]){C c;if(ck(v,c))cs.pb(c);}if(cs.empty()){rs(S);return false;}sort(cs.begin(),cs.end(),[](const C&a,const C&b){return a.s!=b.s?a.s<b.s:a.v<b.v;});int cap=max(12,min(900,SN*76/10000+9)),rm=0;for(const C&c:cs){if(rm>=cap||es()>18.42)break;if(ap(c))++rm;}int AN=cove();double dr=AN>0?(double)(SN-AN)/SN:0;bool keep=false;if(rm>0&&AN>0&&AN<SN&&dr>=.00045&&dr<=.045&&W5::strong_validator()&&es()<18.80){double px=vps(192);keep=px>=(N>18000&&N<30000?.956:.962);}if(!keep){rs(S);return false;}return true;}}namespace WK{static bool run(){bool a=WKO::run();bool b=WX::run();return a||b;}}'''
def addwx(t):
    if "namespace WX{" in t:return t
    sp=ns(t,"WK")
    if not sp:return None
    a,b=sp
    old=t[a:b]
    if not old.startswith("namespace WK{"):return None
    return t[:a]+old.replace("namespace WK{","namespace WKO{",1)+wx("point_tri_d2" in t)+t[b:]
def min1(t):
    o=[];i=0;n=len(t);bol=True
    while i<n:
        c=t[i]
        if bol and c=="#":
            j=t.find("\n",i)
            if j<0:j=n
            line=t[i:j].replace("#include <","#include<")
            o.append(line)
            if j<n:o.append("\n")
            i=j+1;bol=True;continue
        bol=False
        if c in"'\"":
            q=c;o.append(c);i+=1
            while i<n:
                o.append(t[i])
                if t[i]=="\\"and i+1<n:
                    i+=1;o.append(t[i])
                elif t[i]==q:
                    i+=1;break
                i+=1
            continue
        if c=="/"and i+1<n and t[i+1]=="/":
            j=t.find("\n",i+2)
            if j<0:break
            i=j+1;bol=True;continue
        if c=="/"and i+1<n and t[i+1]=="*":
            j=t.find("*/",i+2);i=n if j<0 else j+2;continue
        if c.isspace():
            j=i+1
            while j<n and t[j].isspace():j+=1
            p=o[-1]if o else"";nx=t[j]if j<n else""
            if(p.isalnum()or p=="_"or p==".")and(nx.isalnum()or nx=="_"or nx=="."):o.append(" ")
            if"\n"in t[i:j]:bol=True
            i=j;continue
        o.append(c);i+=1
    r="".join(o)
    r=re.sub(r'(?<![\w.])0\.(\d+)',r'.\1',r)
    r=re.sub(r'(?<![\w.])([1-9][0-9]*)\.0(?![\w.])',r'\1.',r)
    return r
kw=set('alignas alignof asm auto bool break case catch char char8_t char16_t char32_t class const consteval constexpr constinit const_cast continue co_await co_return co_yield decltype default delete do double dynamic_cast else enum explicit export extern false float for friend goto if inline int long mutable namespace new noexcept nullptr operator private protected public register reinterpret_cast requires return short signed sizeof static static_assert static_cast struct switch template this thread_local throw true try typedef typeid typename union unsigned using virtual void volatile wchar_t while'.split())
lib=set('std vector array pair string priority_queue queue map sort stable_sort min max sqrt fabs abs floor ceil cos acos sin atan2 hypot isfinite strtod strtol fread fwrite printf snprintf setvbuf memset memcpy memcmp malloc free lower_bound upper_bound adjacent_find binary_search find fill unique reverse swap move less greater begin end size resize reserve assign clear empty data insert erase push_back pop_back front back top pop push first second chrono steady_clock now time_point duration count uint64_t uint32_t uint16_t uint8_t int64_t int32_t size_t INT_MAX LLONG_MAX UINT_MAX _IOFBF stdin stdout stderr FILE EOF NULL true false nullptr pb cove es pth rs vps svs main include define ifdef ifndef endif elif else'.split())
def masks(t):
    m=[False]*len(t);i=0;n=len(t);bol=True
    while i<n:
        if bol and i<n and t[i]=="#":
            j=t.find("\n",i);j=n if j<0 else j
            for k in range(i,j):m[k]=True
            i=j+1;bol=True;continue
        bol=False
        if i<n and t[i] in"'\"":
            q=t[i];st=i;i+=1
            while i<n:
                if t[i]=="\\":
                    i+=2;continue
                if t[i]==q:
                    i+=1;break
                i+=1
            for k in range(st,min(i,n)):m[k]=True
            continue
        if i<n and t[i]=="\n":bol=True
        i+=1
    return m
def rename(t):
    mm=masks(t);tok=[];prot=set();cnt={}
    for z in re.finditer(r'[A-Za-z_][A-Za-z0-9_]*',t):
        a,b=z.span();x=z.group()
        if any(mm[a:b]):prot.add(x);continue
        cnt[x]=cnt.get(x,0)+1
        j=a-1
        while j>=0 and t[j].isspace():j-=1
        if j>=0 and t[j]=="." or j>=1 and t[j-1:j+1]=="->" or j>=1 and t[j-1:j+1]=="::":prot.add(x)
    prot|=kw|lib
    allid=set(cnt)|prot
    abc="abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ"
    names=[]
    for a in abc:
        for b in abc+"0123456789":
            y=a+b
            if y not in allid and y not in kw and y not in lib:names.append(y)
    for a in abc:
        for b in abc+"0123456789":
            for c in abc+"0123456789":
                y=a+b+c
                if y not in allid and y not in kw and y not in lib:names.append(y)
    cand=[]
    for x,c in cnt.items():
        if x in prot or len(x)<4:continue
        cand.append(((len(x)-2)*c,x))
    cand.sort(reverse=True)
    mp={}
    ni=0
    for sav,x in cand:
        if ni>=len(names):break
        if sav<=2:continue
        mp[x]=names[ni];ni+=1
    if not mp:return t
    out=[];i=0;n=len(t)
    for z in re.finditer(r'[A-Za-z_][A-Za-z0-9_]*',t):
        a,b=z.span();out.append(t[i:a])
        x=z.group()
        out.append(x if any(mm[a:b]) else mp.get(x,x))
        i=b
    out.append(t[i:])
    return "".join(out)
def mainpatch(t):
    olds=[
"int main(){JC();GN();if(!W2G::run())W2C::run();W5::post_patch_pass();VIMP::run();MIDEC::run();WK::run();for(int i=2;i--&&N>47500&N<6e4&&es()<18.6;)WK::run();if(N<50625&&es()<18.9)WK::run();JD();}",
"int main(){JC();GN();if(!W2G::run())W2C::run();W5::post_patch_pass();VIMP::run();MIDEC::run();WKO::run();for(int i=2;i--&&N>47500&N<6e4&&es()<18.6;)WKO::run();if(N<50625&&es()<18.9)WKO::run();JD();}"
]
    new="int main(){JC();GN();if(!W2G::run())W2C::run();W5::post_patch_pass();VIMP::run();MIDEC::run();WK::run();for(int i=2;i--&&N>47500&N<6e4&&es()<18.6;)WK::run();if(N<50625&&es()<18.9)WK::run();JD();}"
    for o in olds:
        if o in t:return t.replace(o,new,1)
    return t
def fit(t):
    xs=[]
    a=mainpatch(t);xs.append(a)
    b=addwx(a)
    if b:xs.insert(0,b)
    for x in xs:
        y=min1(x)
        if len(y.encode())<=L:return y
        z=rename(y)
        z=min1(z)
        if len(z.encode())<=L:return z
    y=min1(mainpatch(t))
    if len(y.encode())<=L:return y
    z=min1(rename(y))
    if len(z.encode())<=L:return z
    return r'''#include<bits/stdc++.h>
using namespace std;int main(){ios::sync_with_stdio(0);cin.tie(0);int N,M;if(!(cin>>N>>M))return 0;vector<string>v(N),f(M);string c;double x,y,z;for(int i=0;i<N;i++){cin>>c>>x>>y>>z;ostringstream o;o<<setprecision(17)<<"v "<<x<<" "<<y<<" "<<z<<"\n";v[i]=o.str();}int a,b,d;for(int i=0;i<M;i++){cin>>c>>a>>b>>d;ostringstream o;o<<"f "<<a<<" "<<b<<" "<<d<<"\n";f[i]=o.str();}cout<<N<<" "<<M<<"\n";for(auto&s:v)cout<<s;for(auto&s:f)cout<<s;}'''
o=fit(s)
if len(o.encode())>L:o=min1(o)
sys.stdout.write(o)