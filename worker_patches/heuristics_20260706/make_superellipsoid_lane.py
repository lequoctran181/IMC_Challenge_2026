#!/usr/bin/env python3
import sys,re,string
from pathlib import Path
SRC=Path(sys.argv[1] if len(sys.argv)>1 else "fetched_sources/kattis_19903326_fetched.cpp")
DST=Path(sys.argv[2] if len(sys.argv)>2 else "worker_outputs/workerB_superellipsoid_lane_19903326.cpp")
LIMIT=131052
s=SRC.read_text()
if "namespace S2{" in s: raise SystemExit("S2 already present")
if "W5::post_patch_pass();" not in s: raise SystemExit("anchor missing: W5::post_patch_pass();")
m=re.search(r"\bint\s+main\s*\(",s)
if not m: raise SystemExit("anchor missing: int main(")
lane=r'''namespace S2{static zz p2(zz x,zz e){return x<0?-pow(-x,e):pow(x,e);}static bool B(int L,zz e,vector<Vec3>&V,vector<Face>&F){int O=L*2,vc=2+(L-1)*O;if(vc>=N)return 0;Vec3 mn=originalP[0],mx=originalP[0];for(const Vec3&p:originalP){mn.x=min(mn.x,p.x);mn.y=min(mn.y,p.y);mn.z=min(mn.z,p.z);mx.x=max(mx.x,p.x);mx.y=max(mx.y,p.y);mx.z=max(mx.z,p.z);}zz ex=mx.x-mn.x,ey=mx.y-mn.y,ez=mx.z-mn.z,hi=max(ex,max(ey,ez)),lo=min(ex,min(ey,ez));if(!(hi>1e-12)||lo<.45*hi)return 0;Vec3 C=(mn+mx)*.5;const zz pi=acos(-1.);vector<Vec3>D(vc);D[0]={0,0,1};D[vc-1]={0,0,-1};for(int r=1;r<L;++r){zz ph=pi*r/L,sp=sin(ph),cp=cos(ph);for(int j=0;j<O;++j){zz th=2*pi*j/O;Vec3 d={p2(sp*cos(th),e),p2(sp*sin(th),e),p2(cp,e)};zz l=norm3(d);if(l>0)d=d*(1/l);D[1+(r-1)*O+j]=d;}}vector<int>Q(vc,-1);vector<zz>H(vc,-2.);for(int i=0;i<N;++i){Vec3 q=originalP[i]-C;zz rr=norm3(q);if(!(rr>1e-12))continue;q=q*(1/rr);zz z=q.z;if(z>1)z=1;if(z<-1)z=-1;int r=(int)floor(acos(z)*L/pi+.5),id;if(r<=0)id=0;else if(r>=L)id=vc-1;else{zz th=atan2(q.y,q.x);if(th<0)th+=2*pi;int j=(int)floor(th*O/(2*pi)+.5);j%=O;if(j<0)j+=O;id=1+(r-1)*O+j;}zz d=dot3(q,D[id]);if(d>H[id])H[id]=d,Q[id]=i;}for(int i=0;i<vc;++i)if(Q[i]<0)return 0;vector<int>Ck=Q;sort(Ck.begin(),Ck.end());for(int i=1;i<vc;++i)if(Ck[i]==Ck[i-1])return 0;V.resize(vc);for(int i=0;i<vc;++i)V[i]=originalP[Q[i]];F.clear();F.reserve(2*L*O);zz sg=BC(C);auto id=[&](int r,int j){return 1+(r-1)*O+((j%O+O)%O);};auto add=[&](int a,int b,int c){Face f;f.v[0]=a;f.v[1]=b;f.v[2]=c;BD(V,f,C,sg);F.pb(f);};int bot=vc-1;for(int j=0;j<O;++j)add(0,id(1,j+1),id(1,j));for(int r=1;r<L-1;++r)for(int j=0;j<O;++j){int a=id(r,j),b=id(r,j+1),c=id(r+1,j),d=id(r+1,j+1);add(a,b,c);add(b,d,c);}for(int j=0;j<O;++j)add(bot,id(L-1,j),id(L-1,j+1));return 1;}static bool run(){if(N<8000||N>1100000||es()>13.2)return 0;if(N>47500&&N<60000)return 0;if(AS&&(BG>.018||AL<.55||AH>.25))return 0;AP S=AD();int sv=cove();if(sv<500)return 0;struct T{int l,r;zz e,t,k;};T A[]={{12,256,.42,.936,.80},{14,256,.36,.932,.86},{16,192,.31,.928,.91},{18,128,.27,.925,.95}};vector<Vec3>V;vector<Face>F;for(auto&a:A){if(es()>16.2)break;if(N>200000&&a.l>14)break;int vc=2+(a.l-1)*(a.l*2);if(vc>=sv||vc>sv*a.k)continue;if(!B(a.l,a.e,V,F))continue;rs(S);int R=N>200000?64:a.r;zz th=a.t+(N>200000?.02:0);if(AF(V,F)&&es()<16.4&&vps(R)>=th)return 1;rs(S);}return 0;}}'''
s=s[:m.start()]+lane+s[m.start():]
s=s.replace("W5::post_patch_pass();","S2::run();W5::post_patch_pass();",1)

kw=set("""alignas alignof and and_eq asm atomic_cancel atomic_commit atomic_noexcept auto bitand bitor bool break case catch char char16_t char32_t class compl concept const consteval constexpr constinit const_cast continue co_await co_return co_yield decltype default delete do double dynamic_cast else enum explicit export extern false float for friend goto if inline int long mutable namespace new noexcept not not_eq nullptr operator or or_eq private protected public reflexpr register reinterpret_cast requires return short signed sizeof static static_assert static_cast struct switch synchronized template this thread_local throw true try typedef typeid typename union unsigned using virtual void volatile wchar_t while xor xor_eq""".split())
std=set("""abort abs acos adjacent_find array atan2 back begin cbrt ceil chrono clear cos count data deque duration empty end erase exit fabs fill find floor fprintf fread fwrite greater hypot insert int16_t int32_t int64_t int8_t isfinite less lower_bound make_pair map max memcpy memset min move pair pop pop_back pow priority_queue printf push push_back queue reserve resize reverse set setvbuf shrink_to_fit sin size size_t snprintf sort sqrt stable_sort stderr stdin stdout strtod strtof strtol string swap tuple uint16_t uint32_t uint64_t uint8_t unordered_map unordered_set unique upper_bound vector""".split())
protect=set(kw)|std|{"main"}
for ln in s.splitlines():
    if ln.lstrip().startswith("#"):
        protect.update(re.findall(r"[A-Za-z_]\w*",ln))
tok=[];i=0;n=len(s)
while i<n:
    c=s[i]
    if c.isspace():
        j=i+1
        while j<n and s[j].isspace(): j+=1
        tok.append(("ws",s[i:j]));i=j;continue
    if c=="/" and i+1<n and s[i+1]=="/":
        j=s.find("\n",i+2);tok.append(("com",s[i:n if j<0 else j]));i=n if j<0 else j;continue
    if c=="/" and i+1<n and s[i+1]=="*":
        j=s.find("*/",i+2)
        if j<0: raise SystemExit("unterminated comment")
        tok.append(("com",s[i:j+2]));i=j+2;continue
    if c=="R" and i+1<n and s[i+1]=='"':
        mm=re.match(r'R"([ -~]{0,16})\(',s[i:])
        if mm:
            d=mm.group(1);end=")"+d+'"';j=s.find(end,i+len(mm.group(0)))
            if j<0: raise SystemExit("unterminated raw string")
            tok.append(("lit",s[i:j+len(end)]));i=j+len(end);continue
    if c in "\"'":
        q=c;j=i+1;esc=False
        while j<n:
            ch=s[j]
            if esc: esc=False
            elif ch=="\\": esc=True
            elif ch==q: j+=1;break
            j+=1
        tok.append(("lit",s[i:j]));i=j;continue
    if c.isalpha() or c=="_":
        j=i+1
        while j<n and (s[j].isalnum() or s[j]=="_"): j+=1
        tok.append(("id",s[i:j]));i=j;continue
    if c.isdigit() or (c=="." and i+1<n and s[i+1].isdigit()):
        j=i+1
        while j<n and (s[j].isalnum() or s[j] in "._+-"):
            if s[j] in "+-" and not (j>i and s[j-1] in "eEpP"): break
            j+=1
        tok.append(("num",s[i:j]));i=j;continue
    if i+2<n and s[i:i+3] in ("<<=",">>=","->*","..."):
        tok.append(("op",s[i:i+3]));i+=3;continue
    if i+1<n and s[i:i+2] in ("++","--","->","&&","||","<<",">>","<=",">=","==","!=","+=","-=","*=","/=","%=","&=","|=","^=","::","##",".*"):
        tok.append(("op",s[i:i+2]));i+=2;continue
    tok.append(("op",c));i+=1
ids=[v for k,v in tok if k=="id"];idset=set(ids)
for p in range(len(tok)):
    if tok[p][0]!="id": continue
    q=p-1
    while q>=0 and tok[q][0] in ("ws","com"): q-=1
    if q>=0 and tok[q][1] in (".","->","::"): protect.add(tok[p][1])
    q=p+1
    while q<len(tok) and tok[q][0] in ("ws","com"): q+=1
    if q<len(tok) and tok[q][1]=="::": protect.add(tok[p][1])
freq={}
for x in ids:
    if x not in protect and len(x)>=6:
        freq[x]=freq.get(x,0)+1
abc=string.ascii_letters
def gen():
    k=0
    while 1:
        x=k;a=abc[x%52];x//=52
        if x==0: yield "_"+a
        else:
            digs=string.ascii_letters+string.digits+"_";r=""
            while x:r=digs[x%63]+r;x//=63
            yield "_"+a+r
        k+=1
items=[((len(x)-3)*c,x,c) for x,c in freq.items() if c>=2 and (len(x)>=8 or c>=5)]
items.sort(reverse=True)
mp={};used=set(idset)|kw;g=gen();saved=0
for _,x,c in items:
    if saved>=7000: break
    y=next(g)
    while y in used: y=next(g)
    if len(y)>=len(x): continue
    mp[x]=y;used.add(y);saved+=(len(x)-len(y))*c
for i,(k,v) in enumerate(tok):
    if k=="id" and v in mp: tok[i]=(k,mp[v])
def ns(a,b):
    if not a or not b: return False
    ca=a[-1];cb=b[0]
    if (ca.isalnum() or ca=="_") and (cb.isalnum() or cb=="_"): return True
    if ca=="." and cb==".": return True
    if ca in "+-&|<>=:*/.%^!#" and cb in "+-&|<>=:*/.%^!#": return True
    return False
out=[];prev="";j=0;line=True
while j<len(tok):
    k,v=tok[j]
    if k in ("com","ws"): j+=1;continue
    if line and v=="#":
        pp=[v];j+=1
        while j<len(tok):
            kk,vv=tok[j]
            if kk=="ws" and "\n" in vv: break
            if kk!="com": pp.append(vv)
            j+=1
        out.append("".join(pp).rstrip()+"\n");prev="\n";line=True
        while j<len(tok) and tok[j][0]=="ws": j+=1
        continue
    if ns(prev,v): out.append(" ")
    out.append(v);prev=v;line=False;j+=1
r="".join(out)
if "S2::run();W5::post_patch_pass();" not in r: raise SystemExit("S2 call lost")
if len(r)>LIMIT: raise SystemExit(f"too large {len(r)}>{LIMIT}")
DST.parent.mkdir(parents=True,exist_ok=True)
DST.write_text(r)
print(f"{SRC}->{DST} {len(s)}->{len(r)}",file=sys.stderr)