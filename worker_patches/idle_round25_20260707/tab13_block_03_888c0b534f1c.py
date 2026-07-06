#!/usr/bin/env python3
import argparse, hashlib, os, re, string, subprocess, sys
from pathlib import Path

LIMIT=131072
DEFAULTS=[
 "fetched_sources/kattis_19903544_81.938904.cpp",
 "worker_outputs/retarget_19903544_20260706/q35_on_19903544.cpp",
 "fetched_sources/kattis_19903326_fetched.cpp",
 "fetched_sources/kattis_19903326.cpp",
 "submission_608_81.93_7.cpp",
 "submission_597_81.93_7.cpp",
 "submission_585_81.93_7.cpp",
 "submission_563_81.93_7.cpp",
 "submission_543_81.93_7.cpp",
]
SAMPLE="""9 14
v 0.5 0.5 0.5
v 0.5 0.5 -0.5
v 0.5 -0.5 0.5
v 0.5 -0.5 -0.5
v -0.5 0.5 0.5
v -0.5 0.5 -0.5
v -0.5 -0.5 0.5
v -0.5 -0.5 -0.5
v 0.5 0.49 0.49
f 1 3 9
f 1 9 2
f 9 3 4
f 9 4 2
f 5 6 8
f 5 8 7
f 1 2 6
f 1 6 5
f 3 7 8
f 3 8 4
f 1 5 7
f 1 7 3
f 2 4 8
f 2 8 6
"""
REQ=["static AP AD()","static void rs(const AP&s)","W5::strong_validator","visual_proxy_score","count_output_vertices_estimate","static vector<Vec3>originalP","static vector<Face>AR","static bool GD(int u,int v,const AE&params)","struct AE{","AF(","JD();"]
LANE=r'''namespace DIV{static Face f3(int a,int b,int c){Face f;f.v[0]=a;f.v[1]=b;f.v[2]=c;return f;}static double d2(const Vec3&a,const Vec3&b){double x=a.x-b.x,y=a.y-b.y,z=a.z-b.z;return x*x+y*y+z*z;}static Vec3 nr(Vec3 a){double l=sqrt(norm2(a));return l>1e-12?a*(1./l):Vec3{0,0,0};}static bool ar(const Vec3&a,const Vec3&b,const Vec3&c){return norm2(cross3(b-a,c-a))>max(1e-30,CL*CL*CL*CL*1e-24);}struct GH{Vec3 mn,mx,ce;double r,r2,h;int nx,ny,nz;vector<vector<int>>b;int cl(int x,int n){return x<0?0:x>=n?n-1:x;}int ix(double x){return cl((int)((x-mn.x)/h),nx);}int iy(double y){return cl((int)((y-mn.y)/h),ny);}int iz(double z){return cl((int)((z-mn.z)/h),nz);}int ky(int x,int y,int z){return(x*ny+y)*nz+z;}bool in(double R){if(!N)return 0;r=R;r2=R*R;mn=mx=originalP[0];for(const Vec3&p:originalP){mn.x=min(mn.x,p.x);mn.y=min(mn.y,p.y);mn.z=min(mn.z,p.z);mx.x=max(mx.x,p.x);mx.y=max(mx.y,p.y);mx.z=max(mx.z,p.z);}ce=(mn+mx)*.5;double sx=max(1e-12,mx.x-mn.x),sy=max(1e-12,mx.y-mn.y),sz=max(1e-12,mx.z-mn.z);h=max(R,max(max(sx,sy),sz)/110.);for(int it=0;it<7;it++){nx=max(1,(int)(sx/h)+3);ny=max(1,(int)(sy/h)+3);nz=max(1,(int)(sz/h)+3);if(1LL*nx*ny*nz<=1450000)break;h*=1.29;}if(1LL*nx*ny*nz>1900000)return 0;b.assign((size_t)nx*ny*nz,{});for(int i=0;i<N;i++)b[ky(ix(originalP[i].x),iy(originalP[i].y),iz(originalP[i].z))].push_back(i);return 1;}void mark(const Vec3&p,vector<unsigned char>&m,int&cc){int X=ix(p.x),Y=iy(p.y),Z=iz(p.z);for(int x=X-1;x<=X+1;x++)if(x>=0&&x<nx)for(int y=Y-1;y<=Y+1;y++)if(y>=0&&y<ny)for(int z=Z-1;z<=Z+1;z++)if(z>=0&&z<nz)for(int q:b[ky(x,y,z)])if(!m[q]&&d2(originalP[q],p)<=r2){m[q]=1;++cc;}}int best(int s,const vector<unsigned char>&m){int X=ix(originalP[s].x),Y=iy(originalP[s].y),Z=iz(originalP[s].z),bi=s,bc=-1;for(int x=X-1;x<=X+1;x++)if(x>=0&&x<nx)for(int y=Y-1;y<=Y+1;y++)if(y>=0&&y<ny)for(int z=Z-1;z<=Z+1;z++)if(z>=0&&z<nz)for(int cand:b[ky(x,y,z)]){int c=0,XX=ix(originalP[cand].x),YY=iy(originalP[cand].y),ZZ=iz(originalP[cand].z);for(int a=XX-1;a<=XX+1;a++)if(a>=0&&a<nx)for(int bb=YY-1;bb<=YY+1;bb++)if(bb>=0&&bb<ny)for(int cc=ZZ-1;cc<=ZZ+1;cc++)if(cc>=0&&cc<nz)for(int q:b[ky(a,bb,cc)])if(!m[q]&&d2(originalP[q],originalP[cand])<=r2)++c;if(c>bc){bc=c;bi=cand;}}return bi;}long long key(const Vec3&p){return ky(ix(p.x),iy(p.y),iz(p.z));}};static vector<Vec3> normals(){vector<Vec3>v(N,{0,0,0});for(const Face&f:AR){Vec3 c=cross3(originalP[f.v[1]]-originalP[f.v[0]],originalP[f.v[2]]-originalP[f.v[0]]);v[f.v[0]]=v[f.v[0]]+c;v[f.v[1]]=v[f.v[1]]+c;v[f.v[2]]=v[f.v[2]]+c;}return v;}static Vec3 bury(int i,const GH&g,const vector<Vec3>&vn,double s){Vec3 n=nr(vn[i]);Vec3 p=originalP[i];if(norm2(n)<1e-20)n=nr(p-g.ce);if(dot3(n,p-g.ce)<0)n=n*(-1.);Vec3 q=p-n*s;if(d2(q,p)>.0485*.0485*CL*CL)q=p-n*(.0485*CL);return q;}static bool bp(vector<Vec3>&X,vector<Face>&Q,const vector<Vec3>&C,int l,int r){int n=r-l;if(n<5){double e=max(1e-8,CL*3e-6);for(int i=l;i<r;i++){Vec3 p=C[i];int s=X.size();X.push_back(p);X.push_back(Vec3{p.x+e,p.y,p.z});X.push_back(Vec3{p.x,p.y+e,p.z});X.push_back(Vec3{p.x,p.y,p.z+e});Q.push_back(f3(s,s+2,s+1));Q.push_back(f3(s,s+1,s+3));Q.push_back(f3(s,s+3,s+2));Q.push_back(f3(s+1,s+2,s+3));}return 1;}int a=l,b=l+1;double bd=-1;for(int i=l;i<r;i++)for(int j=i+1;j<r;j++){double d=d2(C[i],C[j]);if(d>bd)bd=d,a=i,b=j;}Vec3 ax=nr(C[b]-C[a]),t=fabs(ax.x)<.7?Vec3{1,0,0}:Vec3{0,1,0},e1=nr(cross3(ax,t)),e2=cross3(ax,e1),ce{0,0,0};if(norm2(e1)<1e-20){for(int i=l;i<r;i++)bp(X,Q,C,i,i+1);return 1;}for(int i=l;i<r;i++)ce=ce+C[i];ce=ce*(1./n);vector<int>R;for(int i=l;i<r;i++)if(i!=a&&i!=b)R.push_back(i);sort(R.begin(),R.end(),[&](int u,int v){Vec3 U=C[u]-ce,V=C[v]-ce;return atan2(dot3(U,e2),dot3(U,e1))<atan2(dot3(V,e2),dot3(V,e1));});for(int i=0;i<(int)R.size();i++){int u=R[i],v=R[(i+1)%R.size()];if(!ar(C[a],C[u],C[v])||!ar(C[b],C[v],C[u])){for(int j=l;j<r;j++)bp(X,Q,C,j,j+1);return 1;}}int s=X.size();for(int i=l;i<r;i++)X.push_back(C[i]);for(int i=0;i<(int)R.size();i++){int u=R[i]-l,v=R[(i+1)%R.size()]-l;Q.push_back(f3(s+a-l,s+u,s+v));Q.push_back(f3(s+b-l,s+v,s+u));}return 1;}static bool pull(vector<Vec3>&X,vector<Face>&Q){vector<int>id(N,-1);X.clear();Q.clear();X.reserve(cove());Q.reserve(BE);for(int i=0;i<(int)faces.size();i++)if(BR[i]){Face f=faces[i],g;for(int k=0;k<3;k++){int v=f.v[k];if(v<0||v>=N||!BU[v])return 0;if(id[v]<0){id[v]=X.size();X.push_back(P[v]);}g.v[k]=id[v];}if(g.v[0]!=g.v[1]&&g.v[0]!=g.v[2]&&g.v[1]!=g.v[2])Q.push_back(g);}return!X.empty()&&!Q.empty();}struct O{long long k;Vec3 p;};static bool cover(vector<Vec3>&X,vector<Face>&Q,double sh,int cap,double lim){GH g;if(!g.in(.0490*CL))return 0;vector<Vec3>vn=normals();vector<unsigned char>m(N,0);int cc=0;for(const Vec3&p:X){g.mark(p,m,cc);if(es()>lim)return 0;}vector<O>A;A.reserve(max(16,N/700));for(int i=0;i<N&&cc<N;i++){if((i&2047)==0&&es()>lim)return 0;if(!m[i]){int q=g.best(i,m);Vec3 p=bury(q,g,vn,sh);A.push_back({g.key(p),p});g.mark(p,m,cc);if((int)X.size()+(int)A.size()>cap)return 0;}}if(cc<N)return 0;sort(A.begin(),A.end(),[](const O&a,const O&b){return a.k<b.k;});vector<Vec3>C;C.reserve(A.size());for(auto&o:A)C.push_back(o.p);int ch=N>180000?42:34;for(int i=0;i<(int)C.size();i+=ch){if(es()>lim)return 0;bp(X,Q,C,i,min((int)C.size(),i+ch));if((int)X.size()>cap)return 0;}return(int)X.size()<=cap;}static int sweep(const AE&p,int tgt,double lim){int r=0,S=faces.size();for(int pass=0;pass<5&&cove()>tgt&&es()<lim;pass++){int off=(pass*9973)%max(1,S);for(int ii=0;ii<S&&cove()>tgt&&es()<lim;ii++){int i=(ii+off)%S;if(!BR[i])continue;Face f=faces[i];if(GD(f.v[0],f.v[1],p))++r;if(cove()<=tgt||es()>lim)break;if(GD(f.v[1],f.v[2],p))++r;if(cove()<=tgt||es()>lim)break;if(GD(f.v[2],f.v[0],p))++r;}}return r;}struct Pm{double ai,bb,bq,w,at,need;int pct,res;};static bool put(vector<Vec3>&X,vector<Face>&Q,int best,int res,double need){AP S=AD();bool ok=0;if((int)X.size()<best&&AF(X,Q)&&W5::strong_validator()&&cove()<best){double v=vps(res);if(v>=need&&(v>.970||es()>20.55||vps(min(768,res*2))>=need-.006))ok=1;}if(ok)return 1;rs(S);for(auto&f:Q)swap(f.v[1],f.v[2]);if((int)X.size()<best&&AF(X,Q)&&W5::strong_validator()&&cove()<best){double v=vps(res);if(v>=need&&(v>.970||es()>20.55||vps(min(768,res*2))>=need-.006))ok=1;}if(!ok)rs(S);return ok;}static bool run(){if(N<1500||es()>20.25)return 0;int base=cove();if(base<96||base>=N)return 0;AP Start=AD(),Best;int BV=base;bool have=0;Pm T[]={{.075,.020,.965,.012,.985,.936,78,512},{.115,.034,.920,.020,.965,.925,61,384},{.170,.060,.840,.036,.900,.914,45,256},{.250,.095,.660,.058,.820,.906,31,192}};for(Pm&t:T){if(es()>20.45)break;rs(Start);AE p;p.AI=t.ai*CL;p.BB=t.bb*CL;p.BQ=t.bq;p.W=t.w*CL;p.AT=t.at;p.AJ=max(1e-30,CL*CL*1e-24);p.AA=false;int tgt=max(12,base*t.pct/100);sweep(p,tgt,20.10);vector<Vec3>X;vector<Face>Q;if(!pull(X,Q)||X.size()>=BV)continue;double H[5]={.018,.026,.034,.041,.046};for(double h:H){if(es()>20.55)break;vector<Vec3>YV=X;vector<Face>YF=Q;if(!cover(YV,YF,h*CL,min(N-1,BV-1),20.62))continue;if((int)YV.size()>=BV)continue;if(put(YV,YF,BV,t.res,t.need)){Best=AD();BV=cove();have=1;}}}if(have){rs(Best);return 1;}rs(Start);return 0;}}'''
OPS3=("<<=",">>=","->*","...")
OPS2=("++","--","->","&&","||","<<",">>","<=",">=","==","!=","+=","-=","*=","/=","%=","&=","|=","^=","::","##",".*")
OPC=set("+-*&|<>=:*/.%^!#")
def die(s): raise SystemExit("FAIL_CLOSED: "+s)
def choose(p):
    if p:
        q=Path(p)
        if not q.exists(): die("source not found: "+p)
        return q
    for x in DEFAULTS:
        q=Path(x)
        if q.exists(): return q
    for pat in ("*19903544*.cpp","*81.93*_7*.cpp","*19903326*.cpp"):
        for q in Path(".").rglob(pat):
            try: s=q.read_text(errors="ignore")
            except Exception: continue
            if "W5::strong_validator" in s and "visual_proxy_score" in s and "static AP AD()" in s: return q
    die("no current-best source found; pass --src fetched_sources/kattis_19903544_81.938904.cpp")
def find_main(src):
    m=re.search(r"\bint\s+main\s*\(\s*\)\s*\{",src)
    if not m: die("main not found")
    o=src.find("{",m.start());d=1;i=o+1;q=None;esc=False
    while i<len(src):
        c=src[i]
        if q:
            if esc: esc=False
            elif c=="\\": esc=True
            elif c==q: q=None
        else:
            if c in "'\"": q=c
            elif c=="/" and i+1<len(src) and src[i+1]=="/":
                j=src.find("\n",i+2);i=len(src) if j<0 else j
            elif c=="/" and i+1<len(src) and src[i+1]=="*":
                j=src.find("*/",i+2);i=len(src) if j<0 else j+1
            elif c=="{": d+=1
            elif c=="}":
                d-=1
                if d==0:return m.start(),o,i
        i+=1
    die("main brace scan failed")
def patch(src):
    for r in REQ:
        if r not in src: die("missing current-best anchor: "+r)
    if "namespace DIV{" in src or "DIV::run()" in src: die("already patched")
    ms,o,e=find_main(src);body=src[o+1:e];j=body.rfind("JD();")
    if j<0: die("JD hook not found")
    return src[:ms]+LANE+"int main(){"+body[:j]+"DIV::run();"+body[j:]+"}"+src[e+1:]
def lex(s):
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
            if j<0: die("unterminated comment")
            tok.append(("com",s[i:j+2]));i=j+2;continue
        if c=="R" and i+1<n and s[i+1]=='"':
            m=re.match(r'R"([ -~]{0,16})\(',s[i:])
            if m:
                end=")"+m.group(1)+'"';j=s.find(end,i+len(m.group(0)))
                if j<0: die("unterminated raw string")
                tok.append(("lit",s[i:j+len(end)]));i=j+len(end);continue
        if c in "'\"":
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
        if i+2<n and s[i:i+3] in OPS3: tok.append(("op",s[i:i+3]));i+=3;continue
        if i+1<n and s[i:i+2] in OPS2: tok.append(("op",s[i:i+2]));i+=2;continue
        tok.append(("op",c));i+=1
    return tok
def need(a,b):
    if not a or not b:return False
    x=a[-1];y=b[0]
    return ((x.isalnum() or x=="_") and (y.isalnum() or y=="_")) or (x=="." and y==".") or (x in OPC and y in OPC)
def asm(tok):
    out=[];prev="";line=True;i=0
    while i<len(tok):
        k,v=tok[i]
        if k in ("ws","com"):i+=1;continue
        if line and v=="#":
            pp=[v];i+=1
            while i<len(tok):
                kk,vv=tok[i]
                if kk=="ws" and "\n" in vv: break
                if kk!="com": pp.append(vv)
                i+=1
            out.append("".join(pp).rstrip()+"\n");prev="\n";line=True
            while i<len(tok) and tok[i][0]=="ws": i+=1
            continue
        if need(prev,v): out.append(" ")
        out.append(v);prev=v;line=False;i+=1
    return "".join(out)
def minify(src):
    src=re.sub(r"^(?:#include\s*<[^>\n]+>\s*\n)+","#include<bits/stdc++.h>\n",src,1)
    src=src.replace('if(0&&"B16P515A"){}',"").replace("nullptr","0").replace("0.0",".0").replace("1.0","1.")
    tok=lex(src)
    kw=set("alignas alignof and and_eq asm auto bitand bitor bool break case catch char char16_t char32_t class compl const constexpr const_cast continue decltype default delete do double dynamic_cast else enum explicit export extern false float for friend goto if inline int long mutable namespace new noexcept not not_eq nullptr operator or or_eq private protected public register reinterpret_cast return short signed sizeof static static_assert static_cast struct switch template this thread_local throw true try typedef typeid typename union unsigned using virtual void volatile wchar_t while xor xor_eq".split())
    std=set("abort abs acos adjacent_find array atan2 begin cbrt ceil chrono clear cos count data deque empty end erase exit fabs fill find floor fprintf fread fwrite greater hypot insert int16_t int32_t int64_t int8_t isfinite less lower_bound make_pair map max memcpy memset min move pair pop pop_back pow priority_queue printf push push_back queue reserve resize reverse set setvbuf sin size size_t snprintf sort sqrt stable_sort stderr stdin stdout strtod strtol string swap tuple uint16_t uint32_t uint64_t uint8_t unordered_map unordered_set unique upper_bound vector first second".split())
    protect=kw|std|{"main"}
    for ln in src.splitlines():
        if ln.lstrip().startswith("#"): protect.update(re.findall(r"[A-Za-z_]\w*",ln))
    ids=[v for k,v in tok if k=="id"];idset=set(ids)
    for p,(k,v) in enumerate(tok):
        if k!="id": continue
        q=p-1
        while q>=0 and tok[q][0] in ("ws","com"): q-=1
        if q>=0 and tok[q][1] in (".","->","::"): protect.add(v)
        q=p+1
        while q<len(tok) and tok[q][0] in ("ws","com"): q+=1
        if q<len(tok) and tok[q][1]=="::": protect.add(v)
    freq={}
    for x in ids:
        if x not in protect and len(x)>=6: freq[x]=freq.get(x,0)+1
    abc=string.ascii_letters
    def gen():
        k=0;digs=string.ascii_letters+string.digits+"_"
        while True:
            x=k;a=abc[x%52];x//=52;r=""
            while x:r=digs[x%63]+r;x//=63
            k+=1;yield "_"+a+r
    items=sorted([((len(x)-3)*c,x,c) for x,c in freq.items() if c>=2 and (len(x)>=8 or c>=5)],reverse=True)
    mp={};used=set(idset)|kw;g=gen();saved=0
    for _,x,c in items:
        if saved>47000: break
        y=next(g)
        while y in used: y=next(g)
        if len(y)<len(x): mp[x]=y;used.add(y);saved+=(len(x)-len(y))*c
    for i,(k,v) in enumerate(tok):
        if k=="id" and v in mp: tok[i]=(k,mp[v])
    out=asm(tok)
    if len(out.encode())<=LIMIT:return out
    return macro(out)
def macro(src):
    tok=lex(src);idset={v for k,v in tok if k=="id"}
    words=["double","static","return","vector","const","bool","void","false","true","unsigned","int"]
    names=[];i=0
    while len(names)<len(words):
        nm="Q"+string.ascii_letters[i%52];i+=1
        if nm not in idset:names.append(nm)
    mp=dict(zip(words,names))
    for i,(k,v) in enumerate(tok):
        if k=="id" and v in mp: tok[i]=(k,mp[v])
    body=asm(tok);inc="#include<bits/stdc++.h>\n"
    if body.startswith(inc): body=body[len(inc):]
    return inc+"".join("#define %s %s\n"%(mp[w],w) for w in words)+body
def sample(exe):
    r=subprocess.run([exe if os.sep in exe else "./"+exe],input=SAMPLE,text=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE,timeout=20)
    if r.returncode:
        sys.stderr.write(r.stderr);die("sample run failed")
    first=r.stdout.splitlines()[0].strip() if r.stdout.splitlines() else ""
    if first!="8 12":
        print(r.stdout[:500],file=sys.stderr);die("sample first line "+repr(first))
def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--src")
    ap.add_argument("--out",default="submit_div_cover.cpp")
    ap.add_argument("--exe",default="submit_div_cover")
    ap.add_argument("--cxx",default=os.environ.get("CXX","g++"))
    ap.add_argument("--no-compile",action="store_true")
    args=ap.parse_args()
    p=choose(args.src)
    src=p.read_text()
    if len(src.encode())>LIMIT: die("base source over 131072 bytes")
    code=minify(patch(src))
    for t in ("namespace DIV","DIV::run()","W5::strong_validator()","visual_proxy_score"):
        if t not in code: die("lost token "+t)
    if len(code.encode())>LIMIT: die("generated source too large: %d"%len(code.encode()))
    Path(args.out).write_text(code)
    print("source=%s"%p)
    print("output=%s bytes=%d sha256=%s"%(args.out,len(code.encode()),hashlib.sha256(code.encode()).hexdigest()))
    if not args.no_compile:
        cmd=[args.cxx,"-std=c++17","-O2","-pipe",args.out,"-o",args.exe]
        print("compile="+" ".join(cmd))
        r=subprocess.run(cmd,text=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
        if r.stdout:print(r.stdout,end="")
        if r.stderr:print(r.stderr,end="",file=sys.stderr)
        if r.returncode:die("compile failed")
        sample(args.exe)
        print("sample_first_line=8 12")
if __name__=="__main__":main()
