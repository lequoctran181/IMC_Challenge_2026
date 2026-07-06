#!/usr/bin/env python3
import argparse, hashlib, os, re, shutil, subprocess, sys
from pathlib import Path

LIMIT = 131072
PREFERRED = [
    ("fetched_sources/kattis_19903544_81.938904.cpp", "0248732b691263bfe6dc01760bd8b45b1fe05c56"),
    ("submission_580_81.93_7.cpp", "ebeaa582b2d06dc271c0b2171ce639b6d3ab23e9"),
    ("submission_563_81.93_7.cpp", "5d1031432ced9352b06517bd541461a1a8fc9cc5"),
    ("fetched_sources/19901232.cpp", "7b4c5386a74cef7d3d77b0e4053285409baa74f1"),
    ("submission_543_81.93_7.cpp", "7b4c5386a74cef7d3d77b0e4053285409baa74f1"),
]

LANE = r'''namespace SRF{struct T{double h,s,n;int p,r;};static Vec3 U(Vec3 a){double l=norm3(a);return l>1e-300?a*(1./l):Vec3{1,0,0};}static Vec3 ctr(){Vec3 mn=originalP[0],mx=mn;for(auto&p:originalP){mn.x=min(mn.x,p.x);mn.y=min(mn.y,p.y);mn.z=min(mn.z,p.z);mx.x=max(mx.x,p.x);mx.y=max(mx.y,p.y);mx.z=max(mx.z,p.z);}return(mn+mx)*.5;}static bool add(vector<Vec3>&X,vector<Face>&F,const Vec3&c,Vec3 n,double s,double t,const Vec3&o){n=U(n);if(dot3(n,c-o)<0)n=n*(-1);Vec3 a=fabs(n.x)<.72?Vec3{1,0,0}:(fabs(n.y)<.72?Vec3{0,1,0}:Vec3{0,0,1});Vec3 u=U(cross3(n,a)),v=U(cross3(n,u));int q=X.size();if(q+6>N)return 0;X.pb(c+n*t);X.pb(c-n*t);X.pb(c+u*s);X.pb(c+v*s);X.pb(c-u*s);X.pb(c-v*s);int z[8][3]={{0,2,3},{0,3,4},{0,4,5},{0,5,2},{1,3,2},{1,4,3},{1,5,4},{1,2,5}};for(int i=0;i<8;i++){Face f;f.v[0]=q+z[i][0];f.v[1]=q+z[i][1];f.v[2]=q+z[i][2];F.pb(f);}return 1;}static bool build(const T&t,int base,vector<Vec3>&X,vector<Face>&F){if(N<25000)return 0;double h=t.h*CL,h2=h*h;if(!(h>0))return 0;Vec3 mn=originalP[0],mx=mn;for(auto&p:originalP){mn.x=min(mn.x,p.x);mn.y=min(mn.y,p.y);mn.z=min(mn.z,p.z);mx.x=max(mx.x,p.x);mx.y=max(mx.y,p.y);mx.z=max(mx.z,p.z);}int nx=max(1,(int)((mx.x-mn.x)/h)+3),ny=max(1,(int)((mx.y-mn.y)/h)+3),nz=max(1,(int)((mx.z-mn.z)/h)+3);long long cc=1LL*nx*ny*nz;if(cc>5000000LL)return 0;auto ix=[&](double x,int n,double m){int q=(int)((x-m)/h)+1;return q<0?0:(q>=n?n-1:q);};auto key=[&](int x,int y,int z){return(z*ny+y)*nx+x;};vector<int>head((size_t)cc,-1),nxt,sel;sel.reserve(base/6+16);nxt.reserve(base/6+16);int lim=max(16,base*t.p/600);for(int i=0;i<N;i++){int X0=ix(originalP[i].x,nx,mn.x),Y0=ix(originalP[i].y,ny,mn.y),Z0=ix(originalP[i].z,nz,mn.z);bool ok=0;for(int z=max(0,Z0-1);z<=min(nz-1,Z0+1)&&!ok;z++)for(int y=max(0,Y0-1);y<=min(ny-1,Y0+1)&&!ok;y++)for(int x=max(0,X0-1);x<=min(nx-1,X0+1)&&!ok;x++)for(int e=head[key(x,y,z)];e!=-1;e=nxt[e])if(norm2(originalP[sel[e]]-originalP[i])<=h2){ok=1;break;}if(!ok){if((int)sel.size()>=lim)return 0;nxt.pb(head[key(X0,Y0,Z0)]);head[key(X0,Y0,Z0)]=(int)sel.size();sel.pb(i);}if((i&262143)==0&&es()>17.85)return 0;}if(sel.empty()||6*(int)sel.size()>=base)return 0;vector<Vec3>Nn(N,{0,0,0});for(const Face&f:AR){Vec3 cr=cross3(originalP[f.v[1]]-originalP[f.v[0]],originalP[f.v[2]]-originalP[f.v[0]]);Nn[f.v[0]]=Nn[f.v[0]]+cr;Nn[f.v[1]]=Nn[f.v[1]]+cr;Nn[f.v[2]]=Nn[f.v[2]]+cr;if((f.v[0]&262143)==0&&es()>18.35)return 0;}X.clear();F.clear();X.reserve(sel.size()*6);F.reserve(sel.size()*8);Vec3 o=ctr();double s=t.s*CL,tt=max(1e-7*CL,s*.035);if(s*s*2+tt*tt>.00245*CL*CL)return 0;for(int id:sel){Vec3 n=Nn[id];if(norm2(n)<1e-300)n=originalP[id]-o;if(!add(X,F,originalP[id],n,s,tt,o))return 0;if(((int)X.size()&65535)==0&&es()>18.75)return 0;}return !X.empty()&&!F.empty()&&(int)X.size()<base;}static bool keep(vector<Vec3>&X,vector<Face>&F,AP&B,int base,const T&t){if(X.empty()||F.empty()||(int)X.size()*100>=base*t.p||es()>19.10)return 0;rs(B);bool ok=0;if(AF(X,F)&&W5::strong_validator()&&cove()*100<base*t.p&&es()<19.55){double p=vps(t.r);ok=p>=t.n;if(ok&&t.r<512&&es()<18.95)ok=vps(512)>=t.n-.006;if(ok&&t.r<768&&N<180000&&es()<18.50)ok=vps(768)>=t.n-.009;}if(!ok)rs(B);return ok;}static bool try1(const T&t,AP&B,int base){vector<Vec3>X;vector<Face>F;if(!build(t,base,X,F))return 0;return keep(X,F,B,base,t);}static bool run(){if(N<25000||es()>16.90)return 0;if(AS&&(BG>.04||AL<.45||AH>.40))return 0;AP B=AD();int base=cove();if(base<1500)return 0;T a[4];int n=0;if(N>=500000){a[n++]={.044,.030,.992,68,384};a[n++]={.036,.026,.991,80,512};a[n++]={.030,.022,.990,90,512};}else if(N>=160000){a[n++]={.040,.028,.993,70,512};a[n++]={.033,.024,.992,82,512};a[n++]={.027,.020,.991,92,512};}else if(N>=60000){a[n++]={.035,.024,.994,76,512};a[n++]={.029,.021,.993,88,512};}else{a[n++]={.030,.021,.995,82,512};a[n++]={.025,.018,.994,92,512};}for(int i=0;i<n&&es()<19.05;i++){if(try1(a[i],B,base))return 1;rs(B);}rs(B);return 0;}}'''

MAC = [
    ("O0", "double"), ("O1", "static"), ("O2", "const"), ("O3", "inline"),
    ("O4", "vector"), ("O5", "unsigned"), ("O6", "namespace"), ("O7", "struct"),
    ("O8", "template"), ("O9", "typename"), ("P0", "return false"),
    ("P1", "return true"), ("P2", "continue"), ("P3", "return"),
    ("P4", "int"), ("P5", "bool"), ("P6", "void"), ("P7", "class"),
]
ID = re.compile(r"[A-Za-z_][A-Za-z0-9_]*")
OPS3 = ("<<=", ">>=", "->*", "...")
OPS2 = ("++", "--", "->", "&&", "||", "<<", ">>", "<=", ">=", "==", "!=", "+=", "-=", "*=", "/=", "%=", "&=", "|=", "^=", "::", "##", ".*")
OPC = set("+-*&|<>=:*/.%^!#")

def die(msg):
    raise SystemExit("srf_exact_best_generator_abort: " + msg)

def bsha(data: bytes) -> str:
    return hashlib.sha1(b"blob " + str(len(data)).encode() + b"\0" + data).hexdigest()

def scan(s: str):
    t = []
    i = 0
    n = len(s)
    bol = True
    while i < n:
        c = s[i]
        if c in " \t\r\n":
            if c == "\n":
                bol = True
            i += 1
            continue
        if bol and c == "#":
            j = s.find("\n", i)
            if j < 0:
                j = n
            t.append(("pp", s[i:j] + "\n"))
            i = j + 1
            bol = True
            continue
        bol = False
        if c == "/" and i + 1 < n and s[i + 1] == "/":
            j = s.find("\n", i + 2)
            if j < 0:
                break
            i = j + 1
            bol = True
            continue
        if c == "/" and i + 1 < n and s[i + 1] == "*":
            j = s.find("*/", i + 2)
            if j < 0:
                die("unterminated comment")
            i = j + 2
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
            t.append(("lit", s[i:j]))
            i = j
            continue
        m = ID.match(s, i)
        if m:
            t.append(("id", m.group(0)))
            i = m.end()
            continue
        if c.isdigit() or (c == "." and i + 1 < n and s[i + 1].isdigit()):
            j = i + 1
            while j < n and (s[j].isalnum() or s[j] in "._+-"):
                if s[j] in "+-" and not (j > i and s[j - 1] in "eEpP"):
                    break
                j += 1
            t.append(("num", s[i:j]))
            i = j
            continue
        z = s[i:i + 3]
        if z in OPS3:
            t.append(("op", z))
            i += 3
            continue
        z = s[i:i + 2]
        if z in OPS2:
            t.append(("op", z))
            i += 2
            continue
        t.append(("op", c))
        i += 1
    return t

def ids(s: str):
    return [v for k, v in scan(s) if k == "id"]

def need_space(a: str, b: str) -> bool:
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

def transform(s: str) -> str:
    mp = {b: a for a, b in MAC if " " not in b}
    tok = scan(s)
    out = []
    prev = ""
    i = 0
    while i < len(tok):
        k, v = tok[i]
        if k == "pp":
            if out and out[-1] != "\n":
                out.append("\n")
            out.append(v)
            prev = "\n"
            i += 1
            continue
        if k == "id" and v == "return" and i + 1 < len(tok) and tok[i + 1][0] == "id" and tok[i + 1][1] in ("false", "true"):
            v = "P0" if tok[i + 1][1] == "false" else "P1"
            i += 1
        elif k == "id":
            v = mp.get(v, v)
        if need_space(prev, v):
            out.append(" ")
        out.append(v)
        prev = v
        i += 1
    return "".join(out)

def find_main(s: str):
    m = re.search(r"\bint\s+main\s*\(\s*\)\s*\{", s)
    if not m:
        die("missing main anchor")
    start = m.start()
    open_brace = s.find("{", m.start(), m.end() + 8)
    if open_brace < 0:
        die("main opening brace missing")
    d = 0
    q = None
    esc = False
    for i in range(open_brace, len(s)):
        c = s[i]
        if q:
            if esc:
                esc = False
            elif c == "\\":
                esc = True
            elif c == q:
                q = None
        else:
            if c in "\"'":
                q = c
            elif c == "{":
                d += 1
            elif c == "}":
                d -= 1
                if d == 0:
                    return start, i + 1
    die("unmatched main")

def choose_src(arg: str):
    if arg:
        p = Path(arg)
        return p, dict(PREFERRED).get(str(p).replace(os.sep, "/"))
    for path, sh in PREFERRED:
        p = Path(path)
        if p.exists():
            return p, sh
    die("no current-best base found; run from repo root or pass fetched_sources/kattis_19903544_81.938904.cpp explicitly")

def insert_lane(s: str) -> str:
    req = [
        "static AP AD()", "static void rs(const AP&s)", "static bool AF(",
        "W5::strong_validator", "static double vps(int res)", "static int cove()",
        "static void CA(string&out", "namespace W5{",
    ]
    for x in req:
        if x not in s:
            die("missing functional anchor: " + x)
    for x in ["namespace SRF{", "SRF::run()", "namespace VXH{", "VXH::run()", "namespace K3{", "K3::run()", "namespace C5R{", "C5R::run()"]:
        if x in s:
            die("forbidden stale route anchor present: " + x)
    ca = re.search(r"static\s+void\s+CA\s*\(\s*string\s*&\s*out\s*,\s*const\s+char\s*\*\s*line\s*,\s*int\s+len\s*\)\s*\{", s)
    if not ca:
        die("CA insertion anchor not found")
    s = s[:ca.start()] + LANE + s[ca.start():]
    a, b = find_main(s)
    m = s[a:b]
    j = m.rfind("JD();")
    if j < 0:
        die("JD fallback anchor not found in main")
    ins = "if(SRF::run()){JD();return 0;}"
    m = m[:j] + ins + m[j:]
    return s[:a] + m + s[b:]

def build_output(src: str) -> str:
    patched = insert_lane(src)
    col = [x for x, _ in MAC if x in set(ids(patched))]
    if col:
        die("macro collision: " + ",".join(col))
    u = patched.find("using namespace std;")
    if u < 0:
        die("missing using namespace std")
    patched = patched[:u] + "".join("#define %s %s\n" % (a, b) for a, b in MAC) + patched[u:]
    out = transform(patched)
    if out.count("main(){") != 1:
        die("post-transform main anchor lost")
    if "SRF::run()" not in out:
        die("post-transform SRF call lost")
    if "VXH::run()" in out or "K3::run()" in out or "C5R::run()" in out:
        die("forbidden old route present after transform")
    n = len(out.encode("utf-8"))
    if n > LIMIT:
        die("generated source over limit: %d>%d" % (n, LIMIT))
    return out

def compile_and_sample(outp: Path, sample: Path, static: bool):
    if not shutil.which("g++"):
        die("g++ not found")
    exe = outp.with_suffix("")
    cmd = ["g++", "-std=c++17", "-O2", "-pipe"]
    if static:
        cmd += ["-static", "-s"]
    cmd += [str(outp), "-o", str(exe)]
    subprocess.run(cmd, check=True)
    print("compile=" + " ".join(cmd))
    if sample.exists():
        with sample.open("rb") as f:
            ans = subprocess.check_output([str(exe)], stdin=f, timeout=12)
        first = ans.splitlines()[0].decode("ascii", "replace") if ans else ""
        print("sample_first_line=" + first)
        if first.strip() != "8 12":
            die("sample first line mismatch; expected 8 12")
    else:
        print("sample_first_line=skipped_no_" + str(sample))

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("src", nargs="?", default=None)
    ap.add_argument("-o", "--out", default="srf_exact_best_breakthrough.cpp")
    ap.add_argument("--allow-sha-mismatch", action="store_true")
    ap.add_argument("--no-compile", action="store_true")
    ap.add_argument("--no-static", action="store_true")
    ap.add_argument("--sample", default="sample.in")
    args = ap.parse_args()

    inp, expected = choose_src(args.src)
    raw = inp.read_bytes()
    sh = bsha(raw)
    if expected and sh != expected and not args.allow_sha_mismatch:
        die("base blob sha mismatch got %s expected %s for %s" % (sh, expected, inp))
    src = raw.decode("utf-8")
    out = build_output(src)
    outp = Path(args.out)
    outp.write_text(out, encoding="utf-8", newline="")
    print("base=" + str(inp))
    print("base_blob_sha=" + sh)
    print("output=" + str(outp))
    print("output_bytes=%d" % len(out.encode("utf-8")))
    print("source_limit=%d" % LIMIT)
    print("delta_vs_base=%+d" % (len(out.encode("utf-8")) - len(raw)))
    if not args.no_compile:
        compile_and_sample(outp, Path(args.sample), static=not args.no_static)

if __name__ == "__main__":
    main()
