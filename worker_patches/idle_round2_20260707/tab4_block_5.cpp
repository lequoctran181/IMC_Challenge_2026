#!/usr/bin/env python3
# gen_twf_heightfield.py
#
# Usage:
#   python3 gen_twf_heightfield.py \
#     --src fetched_sources/kattis_19903544_81.938904.cpp \
#     --out submission_twf_heightfield.cpp \
#     --exe submission_twf_heightfield.check
#
# The generator writes exactly one submit-ready C++17 source, compiles it,
# and runs the embedded sample gate expecting first line "8 12".

import argparse
import hashlib
import os
import re
import string
import subprocess
import sys
from pathlib import Path

LIMIT = 131072

DEFAULTS = [
    "fetched_sources/kattis_19903544_81.938904.cpp",
    "fetched_sources/kattis_19903326_fetched.cpp",
    "fetched_sources/kattis_19903326.cpp",
    "fetched_sources/kattis_19903153_81.93_7_worker_breakthrough_BOX3_from_543.cpp",
    "fetched_sources/kattis_19902206_81.93_7.cpp",
    "fetched_sources/kattis_19902388_81.93_7.cpp",
    "submission_608_81.93_7.cpp",
    "submission_597_81.93_7.cpp",
    "submission_585_81.93_7.cpp",
    "submission_563_81.93_7.cpp",
    "submission_543_81.93_7.cpp",
    "fetched_sources/19901232.cpp",
]

REQ = [
    "W5::post_patch_pass();",
    "static AP AD()",
    "static void rs(const AP&s)",
    "visual_proxy_score",
    "count_output_vertices_estimate",
    "namespace W5",
    "static vector<Vec3>originalP",
    "static vector<Face>AR",
    "static double CL=1.",
]

SAMPLE = """9 14
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

LANE = r'''namespace TWF{static double gp(const Vec3&p,int k){return k==0?p.x:(k==1?p.y:p.z);}static void sp(Vec3&p,int k,double v){if(k==0)p.x=v;else if(k==1)p.y=v;else p.z=v;}static Face mf(int a,int b,int c){Face f;f.v[0]=a;f.v[1]=b;f.v[2]=c;return f;}static bool af(vector<Vec3>&X,vector<Face>&F,int a,int b,int c,const Vec3&o){if(a<0||b<0||c<0||a==b||a==c||b==c)return 0;Vec3 cr=cross3(X[b]-X[a],X[c]-X[a]);if(norm2(cr)<1e-28*max(1.,CL*CL))return 0;Face f=mf(a,b,c);Vec3 m=(X[a]+X[b]+X[c])*(1./3.);if(dot3(cr,m-o)<0)swap(f.v[1],f.v[2]);F.pb(f);return 1;}static bool ck(vector<Vec3>&X,vector<Face>&F,int base,int lim){if(X.empty()||F.empty()||(int)X.size()>=base||(int)X.size()>lim||es()>18.45)return 0;AP S=AD();bool ok=0;if(AF(X,F)&&W5::strong_validator()&&cove()<base&&cove()<lim&&es()<18.75){int after=cove();double q=N>180000?.9915:.9955;if(after*100<base*36)q-=.003;if(after*100<max(1,N)*3)q-=.002;double p=vps(512);if(p>=q&&(p>.998||es()>18.60||vps(768)>=q-.002))ok=1;}if(ok)return 1;rs(S);for(auto&f:F)swap(f.v[1],f.v[2]);if(AF(X,F)&&W5::strong_validator()&&cove()<base&&cove()<lim&&es()<18.75){int after=cove();double q=N>180000?.9915:.9955;if(after*100<base*36)q-=.003;if(after*100<max(1,N)*3)q-=.002;double p=vps(512);if(p>=q&&(p>.998||es()>18.60||vps(768)>=q-.002))ok=1;}if(!ok)rs(S);return ok;}static bool build(int h,int U,int V,int base,int lim){if(U<5||V<5||2*U*V>=base||2*U*V>lim||es()>17.70)return 0;int a=(h+1)%3,b=(h+2)%3,C=U*V;Vec3 mn=originalP[0],mx=originalP[0];for(const Vec3&p:originalP){mn.x=min(mn.x,p.x);mn.y=min(mn.y,p.y);mn.z=min(mn.z,p.z);mx.x=max(mx.x,p.x);mx.y=max(mx.y,p.y);mx.z=max(mx.z,p.z);}double A=gp(mx,a)-gp(mn,a),B=gp(mx,b)-gp(mn,b),H=gp(mx,h)-gp(mn,h);if(!(A>.05*CL&&B>.05*CL&&H>.012*CL))return 0;vector<double>lo(C,1e100),hi(C,-1e100);vector<unsigned char>oc(C,0);int occ=0;for(int i=0;i<N;i++){const Vec3&p=originalP[i];int x=(int)((gp(p,a)-gp(mn,a))/A*(U-1)+.5),y=(int)((gp(p,b)-gp(mn,b))/B*(V-1)+.5);if(x<0)x=0;if(x>=U)x=U-1;if(y<0)y=0;if(y>=V)y=V-1;int id=y*U+x;if(!oc[id]){oc[id]=1;occ++;}double z=gp(p,h);if(z<lo[id])lo[id]=z;if(z>hi[id])hi[id]=z;if((i&262143)==0&&es()>17.85)return 0;}int need=N>180000?C*72/100:C*83/100;if(occ<need)return 0;for(int y=0;y<V;y++)for(int x=0;x<U;x++){int id=y*U+x;if(oc[id])continue;int bi=-1,bd=999999;for(int r=1;r<=7&&bi<0;r++)for(int yy=y-r;yy<=y+r;yy++)if(yy>=0&&yy<V)for(int xx=x-r;xx<=x+r;xx++)if(xx>=0&&xx<U){int k=yy*U+xx;if(oc[k]){int d=(xx-x)*(xx-x)+(yy-y)*(yy-y);if(d<bd){bd=d;bi=k;}}}if(bi<0)return 0;lo[id]=lo[bi];hi[id]=hi[bi];}double av=0,mnB=1e100;int thin=0;for(int y=0;y<V;y++)for(int x=0;x<U;x++){int id=y*U+x;double t=hi[id]-lo[id];av+=t;if(t<.004*CL)thin++;if(x==0||y==0||x==U-1||y==V-1)mnB=min(mnB,t);}av/=max(1,C);if(av<.020*CL||thin*100>C*(N>180000?6:10)||mnB<.003*CL)return 0;vector<Vec3>X(2*C);auto gid=[&](int y,int x){return y*U+x;};for(int y=0;y<V;y++)for(int x=0;x<U;x++){int id=gid(y,x);double pa=gp(mn,a)+A*x/(U-1),pb=gp(mn,b)+B*y/(V-1);Vec3 p;p.x=p.y=p.z=0;sp(p,a,pa);sp(p,b,pb);sp(p,h,hi[id]);X[id]=p;sp(p,h,lo[id]);X[C+id]=p;}vector<Face>F;F.reserve(4*(U-1)*(V-1)+4*(U+V));Vec3 o=(mn+mx)*.5;for(int y=0;y+1<V;y++)for(int x=0;x+1<U;x++){int p=gid(y,x),q=gid(y,x+1),r=gid(y+1,x+1),s=gid(y+1,x);if(!af(X,F,p,q,r,o)||!af(X,F,p,r,s,o))return 0;int bp=C+p,bq=C+q,br=C+r,bs=C+s;if(!af(X,F,bp,br,bq,o)||!af(X,F,bp,bs,br,o))return 0;}for(int y=0;y+1<V;y++){int p=gid(y,0),q=gid(y+1,0),bp=C+p,bq=C+q;if(!af(X,F,p,q,bq,o)||!af(X,F,p,bq,bp,o))return 0;p=gid(y,U-1);q=gid(y+1,U-1);bp=C+p;bq=C+q;if(!af(X,F,p,bq,q,o)||!af(X,F,p,bp,bq,o))return 0;}for(int x=0;x+1<U;x++){int p=gid(0,x),q=gid(0,x+1),bp=C+p,bq=C+q;if(!af(X,F,p,bq,q,o)||!af(X,F,p,bp,bq,o))return 0;p=gid(V-1,x);q=gid(V-1,x+1);bp=C+p;bq=C+q;if(!af(X,F,p,q,bq,o)||!af(X,F,p,bq,bp,o))return 0;}return ck(X,F,base,lim);}static bool run(){if(N<2500||N>1100000||es()>15.65)return 0;if(AS&&(BG>.075||AL<.34||AH>.55))return 0;int base=cove();if(base<700||base>=N)return 0;int lim1=N<60000?max(64,N*9/100):max(128,N*45/1000);int lim=min(base*52/100,lim1);if(lim<96)return 0;double asp[3];Vec3 mn=originalP[0],mx=originalP[0];for(const Vec3&p:originalP){mn.x=min(mn.x,p.x);mn.y=min(mn.y,p.y);mn.z=min(mn.z,p.z);mx.x=max(mx.x,p.x);mx.y=max(mx.y,p.y);mx.z=max(mx.z,p.z);}for(int h=0;h<3&&es()<17.95;h++){int a=(h+1)%3,b=(h+2)%3;double A=gp(mx,a)-gp(mn,a),B=gp(mx,b)-gp(mn,b);asp[h]=A/max(1e-12,B);}int scales[5]={70,85,100,120,145};for(int z=0;z<5&&es()<18.05;z++)for(int h=0;h<3&&es()<18.05;h++){int total=lim*scales[z]/100;if(total<96)continue;int cells=max(25,total/2);double ar=sqrt(max(.20,min(5.,asp[h])));int U=max(5,(int)(sqrt((double)cells)*ar+.5));int V=max(5,cells/max(1,U));if(U*V*2>lim){double s=sqrt((double)lim/(2.*U*V));U=max(5,(int)(U*s));V=max(5,(int)(V*s));}if(U*V*2>=base||U*V*2>lim)continue;if(build(h,U,V,base,lim))return 1;}return 0;}}'''

OPS3 = ("<<=", ">>=", "->*", "...")
OPS2 = ("++", "--", "->", "&&", "||", "<<", ">>", "<=", ">=", "==", "!=", "+=", "-=", "*=", "/=", "%=", "&=", "|=", "^=", "::", "##", ".*")
OPC = set("+-*&|<>=:*/.%^!#")

def die(msg):
    print("TWF generator abort: " + msg, file=sys.stderr)
    sys.exit(2)

def choose_src(arg):
    if arg:
        p = Path(arg)
        if not p.exists():
            die("explicit --src not found: " + arg)
        return p
    for s in DEFAULTS:
        p = Path(s)
        if p.exists():
            return p
    hits = []
    for pat in ("*19903544*.cpp", "*19903326*.cpp", "*81.93*_7*.cpp", "submission_563_81.93_7.cpp", "submission_543_81.93_7.cpp", "19901232.cpp"):
        hits += list(Path(".").rglob(pat))[:25]
    seen = []
    for p in hits:
        if p not in seen:
            seen.append(p)
    for p in seen:
        try:
            x = p.read_text(errors="ignore")
        except Exception:
            continue
        if all(r in x for r in REQ[:6]) and len(x.encode()) <= LIMIT:
            return p
    die("no current-best source found; pass --src fetched_sources/kattis_19903544_81.938904.cpp")

def find_main(src):
    m = re.search(r"\bint\s+main\s*\(\s*\)\s*\{", src)
    if not m:
        die("main() not found")
    o = src.find("{", m.start())
    i = o + 1
    depth = 1
    quote = None
    esc = False
    n = len(src)
    while i < n:
        c = src[i]
        if quote:
            if esc:
                esc = False
            elif c == "\\":
                esc = True
            elif c == quote:
                quote = None
        else:
            if c in "\"'":
                quote = c
            elif c == "/" and i + 1 < n and src[i + 1] == "/":
                j = src.find("\n", i + 2)
                i = n if j < 0 else j
            elif c == "/" and i + 1 < n and src[i + 1] == "*":
                j = src.find("*/", i + 2)
                i = n if j < 0 else j + 1
            elif c == "{":
                depth += 1
            elif c == "}":
                depth -= 1
                if depth == 0:
                    return m.start(), o, i
        i += 1
    die("main brace did not close")

def patch_main(src):
    ms, o, e = find_main(src)
    body = src[o + 1:e]
    if "TWF::run()" in body or "namespace TWF{" in src:
        die("source already has TWF")
    if "JC();GN();" not in body:
        die("main is not current-best family: missing JC();GN();")
    anchor = "W5::post_patch_pass();"
    p = body.find(anchor)
    jd = body.rfind("JD();")
    if p < 0 or jd < 0 or p > jd:
        die("main missing W5/JD patch anchors")
    p += len(anchor)
    nb = body[:p] + "if(!TWF::run()){" + body[p:jd] + "}" + body[jd:]
    return src[:ms] + LANE + "int main(){" + nb + "}" + src[e + 1:]

def lex(s):
    tok = []
    i = 0
    n = len(s)
    while i < n:
        c = s[i]
        if c.isspace():
            j = i + 1
            while j < n and s[j].isspace():
                j += 1
            tok.append(("ws", s[i:j]))
            i = j
            continue
        if c == "/" and i + 1 < n and s[i + 1] == "/":
            j = s.find("\n", i + 2)
            tok.append(("com", s[i:n if j < 0 else j]))
            i = n if j < 0 else j
            continue
        if c == "/" and i + 1 < n and s[i + 1] == "*":
            j = s.find("*/", i + 2)
            if j < 0:
                die("unterminated comment")
            tok.append(("com", s[i:j + 2]))
            i = j + 2
            continue
        if c == "R" and i + 1 < n and s[i + 1] == '"':
            mm = re.match(r'R"([ -~]{0,16})\(', s[i:])
            if mm:
                d = mm.group(1)
                end = ")" + d + '"'
                j = s.find(end, i + len(mm.group(0)))
                if j < 0:
                    die("unterminated raw string")
                tok.append(("lit", s[i:j + len(end)]))
                i = j + len(end)
                continue
        if c in "\"'":
            q = c
            j = i + 1
            esc = False
            while j < n:
                ch = s[j]
                if esc:
                    esc = False
                elif ch == "\\":
                    esc = True
                elif ch == q:
                    j += 1
                    break
                j += 1
            tok.append(("lit", s[i:j]))
            i = j
            continue
        if c.isalpha() or c == "_":
            j = i + 1
            while j < n and (s[j].isalnum() or s[j] == "_"):
                j += 1
            tok.append(("id", s[i:j]))
            i = j
            continue
        if c.isdigit() or (c == "." and i + 1 < n and s[i + 1].isdigit()):
            j = i + 1
            while j < n and (s[j].isalnum() or s[j] in "._+-"):
                if s[j] in "+-" and not (j > i and s[j - 1] in "eEpP"):
                    break
                j += 1
            tok.append(("num", s[i:j]))
            i = j
            continue
        if i + 2 < n and s[i:i + 3] in OPS3:
            tok.append(("op", s[i:i + 3]))
            i += 3
            continue
        if i + 1 < n and s[i:i + 2] in OPS2:
            tok.append(("op", s[i:i + 2]))
            i += 2
            continue
        tok.append(("op", c))
        i += 1
    return tok

def need_space(a, b):
    if not a or not b:
        return False
    x = a[-1]
    y = b[0]
    if (x.isalnum() or x == "_") and (y.isalnum() or y == "_"):
        return True
    if x == "." and y == ".":
        return True
    if x in OPC and y in OPC:
        return True
    return False

def assemble(tok):
    out = []
    prev = ""
    j = 0
    line = True
    while j < len(tok):
        k, v = tok[j]
        if k in ("com", "ws"):
            j += 1
            continue
        if line and v == "#":
            pp = [v]
            j += 1
            while j < len(tok):
                kk, vv = tok[j]
                if kk == "ws" and "\n" in vv:
                    break
                if kk != "com":
                    pp.append(vv)
                j += 1
            out.append("".join(pp).rstrip() + "\n")
            prev = "\n"
            line = True
            while j < len(tok) and tok[j][0] == "ws":
                j += 1
            continue
        if need_space(prev, v):
            out.append(" ")
        out.append(v)
        prev = v
        line = False
        j += 1
    return "".join(out)

def minify(src):
    src = re.sub(r"^(?:#include\s*<[^>\n]+>\s*\n)+", "#include<bits/stdc++.h>\n", src, count=1)
    src = re.sub(r"\bnullptr\b", "0", src)
    src = re.sub(r"\binline\s+", "", src)
    src = src.replace('if(0&&"B16P515A"){}', "")
    src = src.replace("0.0", ".0").replace("1.0", "1.")
    tok = lex(src)
    kw = set(
        "alignas alignof and and_eq asm auto bitand bitor bool break case catch char char16_t char32_t class compl const constexpr "
        "const_cast continue decltype default delete do double dynamic_cast else enum explicit export extern false float for friend goto if "
        "inline int long mutable namespace new noexcept not not_eq nullptr operator or or_eq private protected public register reinterpret_cast "
        "return short signed sizeof static static_assert static_cast struct switch template this thread_local throw true try typedef typeid "
        "typename union unsigned using virtual void volatile wchar_t while xor xor_eq"
        .split()
    )
    std = set(
        "abort abs acos adjacent_find array atan2 back begin cbrt ceil chrono clear cos count data deque empty end erase exit fabs fill "
        "find floor fprintf fread fwrite greater hypot insert int16_t int32_t int64_t int8_t isfinite less lower_bound make_pair map max "
        "memcpy memset min move pair pop pop_back pow priority_queue printf push push_back queue reserve resize reverse set setvbuf sin size "
        "size_t snprintf sort sqrt stable_sort stderr stdin stdout strtod strtol string swap tuple uint16_t uint32_t uint64_t uint8_t "
        "unordered_map unordered_set unique upper_bound vector __builtin_popcount __builtin_popcountll remove_reference remove_const "
        "remove_cv is_same enable_if conditional decay numeric_limits"
        .split()
    )
    protect = set(kw) | std | {"main"}
    for ln in src.splitlines():
        if ln.lstrip().startswith("#"):
            protect.update(re.findall(r"[A-Za-z_]\w*", ln))
    ids = [v for k, v in tok if k == "id"]
    idset = set(ids)
    for p, (k, v) in enumerate(tok):
        if k != "id":
            continue
        q = p - 1
        while q >= 0 and tok[q][0] in ("ws", "com"):
            q -= 1
        if q >= 0 and tok[q][1] in (".", "->", "::"):
            protect.add(v)
        q = p + 1
        while q < len(tok) and tok[q][0] in ("ws", "com"):
            q += 1
        if q < len(tok) and tok[q][1] == "::":
            protect.add(v)
    freq = {}
    for x in ids:
        if x not in protect and len(x) >= 6:
            freq[x] = freq.get(x, 0) + 1
    abc = string.ascii_letters
    def gen():
        k = 0
        while True:
            x = k
            a = abc[x % 52]
            x //= 52
            if x == 0:
                yield "_" + a
            else:
                digs = string.ascii_letters + string.digits + "_"
                r = ""
                while x:
                    r = digs[x % 63] + r
                    x //= 63
                yield "_" + a + r
            k += 1
    items = [((len(x) - 3) * c, x, c) for x, c in freq.items() if c >= 2 and (len(x) >= 8 or c >= 5)]
    items.sort(reverse=True)
    mp = {}
    used = set(idset) | kw
    g = gen()
    saved = 0
    for _, x, c in items:
        if saved >= 32000:
            break
        y = next(g)
        while y in used:
            y = next(g)
        if len(y) >= len(x):
            continue
        mp[x] = y
        used.add(y)
        saved += (len(x) - len(y)) * c
    for i, (k, v) in enumerate(tok):
        if k == "id" and v in mp:
            tok[i] = (k, mp[v])
    out = assemble(tok)
    if len(out.encode()) <= LIMIT:
        return out
    return macro_pack(out)

def macro_pack(src):
    tok = lex(src)
    idset = {v for k, v in tok if k == "id"}
    words = ["double", "static", "return", "vector", "const", "bool", "void", "false", "true", "unsigned"]
    names = []
    k = 0
    while len(names) < len(words):
        nm = "Q" + string.ascii_letters[k % 52]
        k += 1
        if nm not in idset:
            names.append(nm)
    mp = dict(zip(words, names))
    for i, (kind, val) in enumerate(tok):
        if kind == "id" and val in mp:
            tok[i] = (kind, mp[val])
    body = assemble(tok)
    inc = "#include<bits/stdc++.h>\n"
    if body.startswith(inc):
        body = body[len(inc):]
    defs = "".join("#define %s %s\n" % (mp[w], w) for w in words)
    return inc + defs + body

def validate_generated(code):
    for x in ("namespace TWF", "TWF::run()", "W5::strong_validator()", "visual_proxy_score", "if(!TWF::run())"):
        if x not in code:
            die("post-generation validation lost token: " + x)
    if len(code.encode()) > LIMIT:
        die("generated source too large: %d > %d" % (len(code.encode()), LIMIT))

def run_sample(exe):
    cmd = exe if os.sep in exe else "./" + exe
    r = subprocess.run([cmd], input=SAMPLE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, timeout=10)
    if r.returncode != 0:
        if r.stderr:
            print(r.stderr, file=sys.stderr, end="")
        die("sample run failed with exit " + str(r.returncode))
    first = r.stdout.splitlines()[0].strip() if r.stdout.splitlines() else ""
    if first != "8 12":
        print(r.stdout[:500], file=sys.stderr)
        die("sample gate expected first line 8 12, got " + repr(first))

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--src")
    ap.add_argument("--out", default="submission_twf_heightfield.cpp")
    ap.add_argument("--exe", default="submission_twf_heightfield.check")
    ap.add_argument("--cxx", default=os.environ.get("CXX", "g++"))
    ap.add_argument("--no-compile", action="store_true")
    ap.add_argument("--no-sample", action="store_true")
    args = ap.parse_args()

    src_path = choose_src(args.src)
    src = src_path.read_text()
    sz = len(src.encode())
    if sz > LIMIT:
        die("source already over Kattis limit: %d" % sz)
    for r in REQ:
        if r not in src:
            die("missing current-best anchor: " + r)

    patched = patch_main(src)
    out = minify(patched)
    validate_generated(out)

    op = Path(args.out)
    op.write_text(out)
    print("source=" + str(src_path))
    print("input_bytes=" + str(sz))
    print("output=" + str(op))
    print("output_bytes=" + str(len(out.encode())))
    print("sha256=" + hashlib.sha256(out.encode()).hexdigest())

    if not args.no_compile:
        cmd = [args.cxx, "-std=c++17", "-O2", "-pipe", "-static", "-s", str(op), "-o", args.exe]
        print("compile_command=" + " ".join(cmd))
        r = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        if r.stdout:
            print(r.stdout, end="")
        if r.stderr:
            print(r.stderr, end="", file=sys.stderr)
        if r.returncode != 0:
            die("compile failed with exit " + str(r.returncode))
        print("compile_ok=" + args.exe)
        if not args.no_sample:
            run_sample(args.exe)
            print("sample_gate=ok first_line_8_12")

if __name__ == "__main__":
    main()