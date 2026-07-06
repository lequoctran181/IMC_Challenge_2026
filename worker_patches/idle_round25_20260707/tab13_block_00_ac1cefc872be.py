#!/usr/bin/env python3
import argparse, hashlib, os, re, shutil, string, subprocess, sys, tempfile
from pathlib import Path

LIMIT = 131072
DEFAULTS = [
    "fetched_sources/kattis_19903326_fetched.cpp",
    "fetched_sources/kattis_19903326.cpp",
    "fetched_sources/kattis_19903153_81.93_7_worker_breakthrough_BOX3_from_543.cpp",
    "fetched_sources/kattis_19902839_81.93_7.cpp",
    "fetched_sources/kattis_19902206_81.93_7.cpp",
    "fetched_sources/kattis_19902388_81.93_7.cpp",
    "fetched_sources/19901232.cpp",
    "fetched_sources/19901322.cpp",
    "submission_608_81.93_7.cpp",
    "submission_597_81.93_7.cpp",
    "submission_585_81.93_7.cpp",
    "submission_580_81.93_7.cpp",
    "submission_563_81.93_7.cpp",
    "submission_543_81.93_7.cpp",
]
REQ = [
    ("snapshot", ["static AP AD()", "AP AD()"]),
    ("restore", ["static void rs(const AP&s)", "static void rs(const AP &s)", "void rs(const AP&s)"]),
    ("validator", ["W5::strong_validator"]),
    ("proxy", ["vps(", "visual_proxy_score"]),
    ("count", ["cove(", "count_output_vertices_estimate"]),
    ("original", ["originalP"]),
    ("install", ["AF("]),
    ("collapse", ["GD("]),
    ("faces", ["faces"]),
    ("face alive", ["BR"]),
    ("vertex alive", ["BU"]),
    ("current vertices", ["P["]),
    ("read", ["GN();"]),
    ("write", ["JD();"]),
]

LANE = r'''namespace R6C{static Face ff(int a,int b,int c){Face f;f.v[0]=a;f.v[1]=b;f.v[2]=c;return f;}static bool ar(const Vec3&a,const Vec3&b,const Vec3&c){return norm2(cross3(b-a,c-a))>max(1e-30,1e-24*CL*CL);}static bool gt(const Vec3&a,const Vec3&b,const Vec3&c,const Vec3&d){return ar(a,b,c)&&ar(a,b,d)&&ar(a,c,d)&&ar(b,c,d);}static void tiny(vector<Vec3>&X,vector<Face>&Q,const Vec3&p){double e=max(1e-8,CL*7e-6);int s=(int)X.size();X.push_back(p);X.push_back(Vec3{p.x+e,p.y,p.z});X.push_back(Vec3{p.x,p.y+e,p.z});X.push_back(Vec3{p.x,p.y,p.z+e});Q.push_back(ff(s,s+2,s+1));Q.push_back(ff(s,s+1,s+3));Q.push_back(ff(s,s+3,s+2));Q.push_back(ff(s+1,s+2,s+3));}static void tet(vector<Vec3>&X,vector<Face>&Q,Vec3 a,Vec3 b,Vec3 c,Vec3 d){double e=max(1e-9,CL*1e-7);b.x+=e;c.y+=e;d.z+=e;if(!gt(a,b,c,d)){tiny(X,Q,a);tiny(X,Q,b);tiny(X,Q,c);tiny(X,Q,d);return;}int s=(int)X.size();X.push_back(a);X.push_back(b);X.push_back(c);X.push_back(d);Q.push_back(ff(s,s+2,s+1));Q.push_back(ff(s,s+1,s+3));Q.push_back(ff(s,s+3,s+2));Q.push_back(ff(s+1,s+2,s+3));}struct G{Vec3 mn,mx;double r,r2,c;int nx,ny,nz;vector<vector<int>>b;int cl(int x,int n){return x<0?0:(x>=n?n-1:x);}int ix(double x){return cl((int)floor((x-mn.x)/c),nx);}int iy(double y){return cl((int)floor((y-mn.y)/c),ny);}int iz(double z){return cl((int)floor((z-mn.z)/c),nz);}int key(int x,int y,int z){return(z*ny+y)*nx+x;}bool init(double R){if(N<=0||originalP.empty())return 0;r=R;r2=R*R;mn=mx=originalP[0];for(const Vec3&p:originalP){mn.x=min(mn.x,p.x);mn.y=min(mn.y,p.y);mn.z=min(mn.z,p.z);mx.x=max(mx.x,p.x);mx.y=max(mx.y,p.y);mx.z=max(mx.z,p.z);}double sx=max(1e-12,mx.x-mn.x),sy=max(1e-12,mx.y-mn.y),sz=max(1e-12,mx.z-mn.z),sp=max(sx,max(sy,sz));c=max(R,sp/128.);for(int it=0;it<8;++it){nx=max(1,(int)(sx/c)+3);ny=max(1,(int)(sy/c)+3);nz=max(1,(int)(sz/c)+3);if(1LL*nx*ny*nz<=1900000)break;c*=1.32;}if(1LL*nx*ny*nz>2600000)return 0;b.assign((size_t)nx*ny*nz,{});for(int i=0;i<N;++i)b[key(ix(originalP[i].x),iy(originalP[i].y),iz(originalP[i].z))].push_back(i);return 1;}void mark(const Vec3&p,vector<unsigned char>&C,int&cc){int X=ix(p.x),Y=iy(p.y),Z=iz(p.z);for(int z=Z-1;z<=Z+1;++z)if(z>=0&&z<nz)for(int y=Y-1;y<=Y+1;++y)if(y>=0&&y<ny)for(int x=X-1;x<=X+1;++x)if(x>=0&&x<nx)for(int q:b[key(x,y,z)])if(!C[q]&&norm2(originalP[q]-p)<=r2){C[q]=1;++cc;}}};static bool pull(vector<Vec3>&X,vector<Face>&Q){vector<int>id(N,-1);X.clear();Q.clear();X.reserve(max(16,cove()));Q.reserve(faces.size());for(int i=0;i<(int)faces.size();++i)if(BR[i]){Face f=faces[i],g;bool ok=1;for(int k=0;k<3;++k){int v=f.v[k];if(v<0||v>=N||!BU[v]){ok=0;break;}if(id[v]<0){id[v]=(int)X.size();X.push_back(P[v]);}g.v[k]=id[v];}if(ok&&g.v[0]!=g.v[1]&&g.v[0]!=g.v[2]&&g.v[1]!=g.v[2])Q.push_back(g);}return !X.empty()&&!Q.empty();}static unsigned long long sp(unsigned x){unsigned long long v=x&2097151u;v=(v|v<<32)&0x1f00000000ffffULL;v=(v|v<<16)&0x1f0000ff0000ffULL;v=(v|v<<8)&0x100f00f00f00f00fULL;v=(v|v<<4)&0x10c30c30c30c30c3ULL;v=(v|v<<2)&0x1249249249249249ULL;return v;}static unsigned long long mk(const Vec3&p,const Vec3&a,const Vec3&b){auto cv=[&](double x,double l,double h){double t=(x-l)/max(1e-30,h-l);if(t<0)t=0;if(t>1)t=1;return(unsigned)(t*2097151.);};return sp(cv(p.x,a.x,b.x))|(sp(cv(p.y,a.y,b.y))<<1)|(sp(cv(p.z,a.z,b.z))<<2);}static Vec3 inw(const Vec3&p,const Vec3&c,double s){Vec3 d=c-p;double l=sqrt(norm2(d));if(l>1e-12)return p+d*(s/l);return p;}static double md4(const vector<Vec3>&C,int i){double m=0;for(int a=0;a<4;++a)for(int b=a+1;b<4;++b)m=max(m,sqrt(norm2(C[i+a]-C[i+b])));return m;}static bool cover(vector<Vec3>&X,vector<Face>&Q,int cap,double R,double sh,double gm,double lim){if((int)X.size()>=cap||cap<4)return 0;G g;if(!g.init(R))return 0;Vec3 cen{.5*(g.mn.x+g.mx.x),.5*(g.mn.y+g.mx.y),.5*(g.mn.z+g.mx.z)};vector<unsigned char>Cv(N,0);int cc=0;for(int i=0;i<(int)X.size();++i){if((i&4095)==0&&es()>lim)return 0;g.mark(X[i],Cv,cc);}vector<Vec3>C;C.reserve(min(N,cap));for(int i=0;i<N&&cc<N;++i){if((i&4095)==0&&es()>lim)return 0;if(!Cv[i]){Vec3 q=inw(originalP[i],cen,sh);C.push_back(q);g.mark(q,Cv,cc);if((int)X.size()+(int)C.size()>cap)return 0;}}if(cc<N)return 0;struct O{unsigned long long k;Vec3 p;};vector<O>Ov;Ov.reserve(C.size());for(const Vec3&p:C)Ov.push_back({mk(p,g.mn,g.mx),p});sort(Ov.begin(),Ov.end(),[](const O&a,const O&b){return a.k<b.k;});for(int i=0;i<(int)C.size();++i)C[i]=Ov[i].p;for(int i=0;i<(int)C.size();){if((i&255)==0&&es()>lim)return 0;if(i+3<(int)C.size()&&md4(C,i)<=gm&&gt(C[i],C[i+1],C[i+2],C[i+3])){tet(X,Q,C[i],C[i+1],C[i+2],C[i+3]);i+=4;}else{tiny(X,Q,C[i]);++i;}if((int)X.size()>cap)return 0;}return(int)X.size()<=cap;}static bool tr(const AP&B,int S,int rat,double R,double sh,double gm,int res,double need,double lim){rs(B);if(S<100||es()>lim-.35)return 0;vector<Vec3>X;vector<Face>Q;if(!pull(X,Q)||X.empty()||Q.empty()){rs(B);return 0;}if((int)X.size()>=S){rs(B);return 0;}int cap=min(N-1,max((int)X.size()+4,(int)((long long)S*rat/100)));if(cap>=S)cap=S-1;if(!cover(X,Q,cap,R*CL,sh*CL,gm*CL,lim-.18)){rs(B);return 0;}bool ok=0;if(AF(X,Q)&&W5::strong_validator()&&cove()<S&&es()<lim+.08){double p=vps(res);ok=p>=need;if(ok&&res<384&&cove()*100<S*62&&es()<20.55)ok=vps(384)>=need-.010;if(ok&&cove()*100<S*42&&es()<20.85)ok=vps(512)>=need-.018;}if(!ok)rs(B);return ok;}static bool run(){if(N<2000||es()>20.40)return 0;int S=cove();if(S<128||S>N)return 0;AP B=AD();int r0=N>180000?96:(N>70000?128:192);if(es()>19.55){rs(B);return 0;}double b=vps(r0);double n0=max(.904,b-.004),n1=max(.900,b-.010);bool ok=0;if(N<60000)ok=tr(B,S,78,.04920,.010,.092,256,n0,20.15)||tr(B,S,66,.04905,.014,.128,256,n1,20.55);else if(N<180000)ok=tr(B,S,74,.04910,.012,.105,192,n0,20.20)||tr(B,S,60,.04890,.017,.145,256,n1,20.62);else ok=tr(B,S,76,.04915,.012,.115,128,n0,20.18)||tr(B,S,58,.04885,.018,.160,192,n1,20.60);if(!ok)rs(B);return ok;}}'''

SAMPLE = b"""9 14
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

KW = set("""alignas alignof and and_eq asm atomic_cancel atomic_commit atomic_noexcept auto bitand bitor bool break case catch char char16_t char32_t class compl concept const consteval constexpr constinit const_cast continue co_await co_return co_yield decltype default delete do double dynamic_cast else enum explicit export extern false float for friend goto if inline int long mutable namespace new noexcept not not_eq nullptr operator or or_eq private protected public reflexpr register reinterpret_cast requires return short signed sizeof static static_assert static_cast struct switch synchronized template this thread_local throw true try typedef typeid typename union unsigned using virtual void volatile wchar_t while xor xor_eq""".split())
STD = set("""abort abs acos adjacent_find array atan2 begin cbrt ceil chrono clear cos count data deque duration empty end erase exit fabs fill find floor fprintf fread fwrite greater hypot insert int16_t int32_t int64_t int8_t isfinite less lower_bound make_pair map max memcpy memset min move pair pop pop_back pow priority_queue printf push push_back queue reserve resize reverse set setvbuf shrink_to_fit sin size size_t snprintf sort sqrt stable_sort stderr stdin stdout strtod strtof strtol string swap tuple uint16_t uint32_t uint64_t uint8_t unordered_map unordered_set unique upper_bound vector puts perror getenv system""".split())
ID = re.compile(r"[A-Za-z_][A-Za-z0-9_]*")

def die(s):
    raise SystemExit("FAIL_CLOSED: " + s)

def choose_src(arg):
    if arg:
        p = Path(arg)
        if not p.exists():
            die("explicit source not found: " + arg)
        return p
    for s in DEFAULTS:
        p = Path(s)
        if p.exists():
            return p
    hits = []
    for pat in ("*19903326*.cpp", "*81.93*_7*.cpp", "submission_543_81.93_7.cpp", "19901232.cpp", "19901322.cpp"):
        hits.extend(Path(".").rglob(pat))
    seen = []
    for p in hits:
        if p not in seen:
            seen.append(p)
    for p in seen:
        try:
            x = p.read_text(encoding="utf-8", errors="ignore")
        except Exception:
            continue
        if "int main" in x and "W5::strong_validator" in x and ("vps(" in x or "visual_proxy_score" in x) and len(x.encode()) <= LIMIT:
            return p
    die("no current-best source found; pass fetched_sources/kattis_19903326_fetched.cpp explicitly")

def require(src):
    for label, alts in REQ:
        if not any(a in src for a in alts):
            die("missing current-best anchor: " + label)
    if "namespace R6C{" in src or "R6C::run(" in src:
        die("R6C already installed")

def match_brace(s, open_pos):
    d = 0
    i = open_pos
    q = None
    esc = False
    while i < len(s):
        c = s[i]
        if q:
            if esc:
                esc = False
            elif c == "\\":
                esc = True
            elif c == q:
                q = None
        else:
            if c in ("'", '"'):
                q = c
            elif c == "/" and i + 1 < len(s) and s[i + 1] == "/":
                j = s.find("\n", i + 2)
                i = len(s) if j < 0 else j
            elif c == "/" and i + 1 < len(s) and s[i + 1] == "*":
                j = s.find("*/", i + 2)
                i = len(s) if j < 0 else j + 1
            elif c == "{":
                d += 1
            elif c == "}":
                d -= 1
                if d == 0:
                    return i + 1
        i += 1
    return -1

def main_span(src):
    m = re.search(r"\bint\s+main\s*\(\s*\)\s*\{", src)
    if not m:
        die("main() not found")
    o = src.find("{", m.start())
    e = match_brace(src, o)
    if e < 0:
        die("main brace did not close")
    return m.start(), e

def patch(src):
    require(src)
    a, b = main_span(src)
    main = src[a:b]
    g = main.find("GN();")
    j = main.rfind("JD();")
    if g < 0 or j < 0 or j < g:
        die("main lacks GN/JD order")
    main2 = main[:j] + "R6C::run();" + main[j:]
    return src[:a] + LANE + main2 + src[b:]

def tokenise(src):
    tok = []
    i = 0
    n = len(src)
    bol = True
    pp = False
    while i < n:
        c = src[i]
        if bol:
            j = i
            while j < n and src[j] in " \t":
                j += 1
            pp = j < n and src[j] == "#"
            bol = False
        if c == "\n":
            tok.append(("ws", c, pp)); bol = True; pp = False; i += 1; continue
        if c.isspace():
            j = i + 1
            while j < n and src[j].isspace() and src[j] != "\n":
                j += 1
            tok.append(("ws", src[i:j], pp)); i = j; continue
        if pp:
            j = i + 1
            while j < n and src[j] != "\n":
                j += 1
            tok.append(("pp", src[i:j], True)); i = j; continue
        if c in ("'", '"'):
            q = c
            st = i
            i += 1
            esc = False
            while i < n:
                ch = src[i]
                if esc:
                    esc = False
                elif ch == "\\":
                    esc = True
                elif ch == q:
                    i += 1
                    break
                i += 1
            tok.append(("lit", src[st:i], False)); continue
        if c == "/" and i + 1 < n and src[i + 1] == "/":
            j = src.find("\n", i)
            tok.append(("com", src[i:n if j < 0 else j], False))
            i = n if j < 0 else j
            continue
        if c == "/" and i + 1 < n and src[i + 1] == "*":
            j = src.find("*/", i + 2)
            j = n if j < 0 else j + 2
            tok.append(("com", src[i:j], False)); i = j; continue
        m = ID.match(src, i)
        if m:
            tok.append(("id", m.group(0), False)); i = m.end(); continue
        if c.isdigit() or (c == "." and i + 1 < n and src[i + 1].isdigit()):
            j = i + 1
            while j < n and (src[j].isalnum() or src[j] in "._+-"):
                if src[j] in "+-" and not (j > i and src[j - 1] in "eEpP"):
                    break
                j += 1
            tok.append(("num", src[i:j], False)); i = j; continue
        if i + 2 < n and src[i:i+3] in ("<<=", ">>=", "->*", "..."):
            tok.append(("op", src[i:i+3], False)); i += 3; continue
        if i + 1 < n and src[i:i+2] in ("++","--","->","&&","||","<<",">>","<=",">=","==","!=","+=","-=","*=","/=","%=","&=","|=","^=","::","##",".*"):
            tok.append(("op", src[i:i+2], False)); i += 2; continue
        tok.append(("op", c, False)); i += 1
    return tok

def names(src):
    return set(m.group(0) for m in ID.finditer(src))

def minify(src):
    used = names(src) | KW | STD
    pool = [f"Q{i}" for i in range(170)] + [f"Z{i}" for i in range(170)] + [f"Y{i}" for i in range(170)]
    vals = ["double","static","const","inline","vector","unsigned","namespace","struct","template","typename","return false","return true","continue","return","int","bool","void","class"]
    mac = []
    for val in vals:
        name = next((x for x in pool if x not in used), None)
        if name is None:
            die("macro pool exhausted")
        used.add(name)
        mac.append((name, val))
    pos = src.find("using namespace std;")
    if pos < 0:
        die("missing using namespace std")
    src = src[:pos] + "".join(f"#define {a} {b}\n" for a,b in mac) + src[pos:]
    tok = tokenise(src)
    protect = set(KW) | set(STD) | {"main"} | {a for a,_ in mac}
    for k,v,pp in tok:
        if k == "pp":
            protect.update(m.group(0) for m in ID.finditer(v))
    for i,(k,v,pp) in enumerate(tok):
        if k != "id":
            continue
        p = i - 1
        while p >= 0 and tok[p][0] in ("ws","com"):
            p -= 1
        q = i + 1
        while q < len(tok) and tok[q][0] in ("ws","com"):
            q += 1
        if (p >= 0 and tok[p][1] in (".","->","::")) or (q < len(tok) and tok[q][1] == "::"):
            protect.add(v)
    freq = {}
    for k,v,pp in tok:
        if k == "id" and v not in protect and len(v) >= 6:
            freq[v] = freq.get(v, 0) + 1
    alpha = string.ascii_letters
    tail = string.ascii_letters + string.digits + "_"
    def gen():
        k = 0
        while True:
            x = k
            s = "_" + alpha[x % 52]
            x //= 52
            while x:
                s += tail[x % 63]
                x //= 63
            k += 1
            yield s
    g = gen()
    rmap = {}
    used2 = names(src) | KW | {a for a,_ in mac}
    items = sorted([((len(x)-3)*c, x, c) for x,c in freq.items() if c >= 2 and (len(x) >= 8 or c >= 5)], reverse=True)
    saved = 0
    for _, x, c in items:
        if saved > 26000:
            break
        y = next(g)
        while y in used2:
            y = next(g)
        if len(y) < len(x):
            rmap[x] = y
            used2.add(y)
            saved += (len(x) - len(y)) * c
    mkw = {v:a for a,v in mac if " " not in v}
    ret_false = next(a for a,v in mac if v == "return false")
    ret_true = next(a for a,v in mac if v == "return true")
    def out_id(x):
        x = rmap.get(x, x)
        return mkw.get(x, x)
    def need_space(a, b):
        if not a or not b:
            return False
        ca, cb = a[-1], b[0]
        if (ca.isalnum() or ca == "_") and (cb.isalnum() or cb == "_"):
            return True
        if ca == "." and cb == ".":
            return True
        if ca in "+-&|<>=:*/.%^!#" and cb in "+-&|<>=:*/.%^!#":
            return True
        return False
    out = []
    prev = ""
    line = True
    i = 0
    while i < len(tok):
        k, v, pp = tok[i]
        if k in ("ws","com"):
            if k == "ws" and "\n" in v:
                line = True
                prev = "\n"
            i += 1
            continue
        if k == "pp":
            if not line:
                out.append("\n")
            out.append(v.rstrip() + "\n")
            prev = "\n"
            line = True
            i += 1
            continue
        if k == "id" and v == "return":
            j = i + 1
            while j < len(tok) and tok[j][0] in ("ws","com"):
                j += 1
            if j < len(tok) and tok[j][0] == "id" and tok[j][1] in ("false","true"):
                cur = ret_false if tok[j][1] == "false" else ret_true
                if need_space(prev, cur):
                    out.append(" ")
                out.append(cur)
                prev = cur
                line = False
                i = j + 1
                continue
        cur = out_id(v) if k == "id" else v
        if need_space(prev, cur):
            out.append(" ")
        out.append(cur)
        prev = cur
        line = False
        i += 1
    return "".join(out)

def build(srcp, outp):
    raw = srcp.read_text(encoding="utf-8", errors="strict")
    if len(raw.encode()) > LIMIT:
        die(f"input source already exceeds limit: {len(raw.encode())}>{LIMIT}")
    patched = patch(raw).replace('if(0&&"B16P515A"){}', '')
    out = minify(patched)
    if len(out.encode()) > LIMIT:
        out = minify(patched.replace("if(N<50625&&es()<18.9)WK::run();", ""))
    if len(out.encode()) > LIMIT:
        die(f"generated source too large: {len(out.encode())}>{LIMIT}")
    if "R6C::run()" not in out:
        die("R6C hook missing after minify")
    outp.write_text(out, encoding="utf-8", newline="")
    return out

def compile_and_gate(outp, cxx, static, sample_gate):
    if not shutil.which(cxx):
        die("compiler not found: " + cxx)
    exe = outp.with_suffix("")
    cmd = [cxx, "-std=c++17", "-O2", "-pipe"]
    if static:
        cmd += ["-static", "-s"]
    cmd += [str(outp), "-o", str(exe)]
    print("compile=" + " ".join(cmd))
    subprocess.run(cmd, check=True)
    print("compile_ok=" + str(exe))
    if sample_gate:
        r = subprocess.run([str(exe)], input=SAMPLE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=20)
        if r.returncode != 0:
            sys.stderr.write(r.stderr.decode("utf-8", "replace"))
            die("sample run failed")
        first = r.stdout.splitlines()[0].decode("ascii", "replace") if r.stdout else ""
        print("sample_first_line=" + first)
        if first.strip() != "8 12":
            die("sample first line mismatch; expected 8 12")

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("src", nargs="?", default=None)
    ap.add_argument("-o", "--out", default="r6c_groupcover_submit.cpp")
    ap.add_argument("--no-compile", action="store_true")
    ap.add_argument("--no-static", action="store_true")
    ap.add_argument("--no-sample-gate", action="store_true")
    ap.add_argument("--cxx", default=os.environ.get("CXX", "g++"))
    args = ap.parse_args()
    srcp = choose_src(args.src)
    outp = Path(args.out)
    out = build(srcp, outp)
    print("base=" + str(srcp))
    print("output=" + str(outp))
    print("output_bytes=" + str(len(out.encode())))
    print("source_limit=" + str(LIMIT))
    print("sha256=" + hashlib.sha256(out.encode()).hexdigest())
    if not args.no_compile:
        compile_and_gate(outp, args.cxx, not args.no_static, not args.no_sample_gate)

if __name__ == "__main__":
    main()
