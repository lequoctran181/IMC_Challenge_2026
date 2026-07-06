#!/usr/bin/env python3
import hashlib
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

LIMIT = 131072

LANE = r'''namespace PCL2{
struct E{double l;unsigned long long k;};
static unsigned long long K(int a,int b){if(a>b)swap(a,b);return(unsigned long long)(unsigned int)a<<32|(unsigned int)b;}
static void ce(vector<E>&ed,double ml,double tl){
    ed.clear();
    vector<unsigned long long>ks;
    ks.reserve((size_t)max(1,BE)*3);
    for(int i=0;i<(int)faces.size();++i){
        if((i&8191)==0&&es()>tl)return;
        if(!BR[i])continue;
        Face f=faces[i];
        if(f.v[0]<0||f.v[1]<0||f.v[2]<0||f.v[0]>=N||f.v[1]>=N||f.v[2]>=N)continue;
        if(!BU[f.v[0]]||!BU[f.v[1]]||!BU[f.v[2]])continue;
        ks.push_back(K(f.v[0],f.v[1]));
        ks.push_back(K(f.v[1],f.v[2]));
        ks.push_back(K(f.v[2],f.v[0]));
    }
    sort(ks.begin(),ks.end());
    ks.erase(unique(ks.begin(),ks.end()),ks.end());
    const double ml2=ml*ml;
    ed.reserve(ks.size());
    for(unsigned long long k:ks){
        int a=(int)(k>>32),b=(int)(k&0xffffffffu);
        if(a<0||b<0||a>=N||b>=N||!BU[a]||!BU[b])continue;
        double l=norm2(P[a]-P[b]);
        if(l<=ml2)ed.push_back({l,k});
    }
    sort(ed.begin(),ed.end(),[](const E&a,const E&b){return a.l<b.l;});
}
static bool ec(int tgt,double ai,double bb,double deg,double tl){
    if(tgt<4)return false;
    AE p;
    p.AI=ai*CL;
    p.BB=bb*CL;
    p.BQ=cos(deg*acos(-1.)/180.);
    p.W=max(1e-12,1e-10*CL);
    p.AT=-.45;
    p.AJ=max(1e-30,1e-24*CL*CL);
    p.AA=false;
    int cur=cove(),start=cur,rounds=0;
    vector<E>ed;
    while(cur>tgt&&rounds<4&&es()<tl){
        ce(ed,min(.09,ai*1.85)*CL,tl-.15);
        if(ed.empty())break;
        int ok=0,seen=0;
        for(const E&e:ed){
            if(cur<=tgt||es()>tl)break;
            int a=(int)(e.k>>32),b=(int)(e.k&0xffffffffu);
            if(a<0||b<0||a>=N||b>=N||!BU[a]||!BU[b])continue;
            ++seen;
            if(GD(a,b,p)){--cur;++ok;}
            if((seen&4095)==0&&es()>tl)break;
        }
        if(ok==0)break;
        ++rounds;
    }
    return cove()<start;
}
static bool ok(int base,double th,int res){
    int now=cove();
    if(now<=0||now>=base)return false;
    if(es()>20.05)return false;
    if(!W5::strong_validator())return false;
    double q=vps(res);
    if(q<th)return false;
    if(res<512&&N<140000&&es()<20.45){
        double q2=vps(512);
        if(q2<th-.010)return false;
    }
    return cove()<base;
}
static bool tr(const AP&B,AP&best,int&bv,int pct,double ai,double bb,double deg,double th,int res,double tl){
    rs(B);
    int base=bv;
    int tgt=max(4,base*pct/100);
    if(tgt>=base)return false;
    if(!ec(tgt,ai,bb,deg,tl))return false;
    int cv=cove();
    if(cv>=base)return false;
    if(!ok(base,th,res))return false;
    cv=cove();
    if(cv<bv){
        bv=cv;
        best=AD();
        return true;
    }
    return false;
}
static bool run(){
    if(N<800||es()>19.15)return false;
    int base=cove();
    if(base<40)return false;
    AP B=AD(),best=B;
    int bv=base;
    if(N<30000){
        tr(B,best,bv,82,.048,.050,58.,.930,512,19.65);
        tr(B,best,bv,72,.049,.052,64.,.922,512,19.95);
        tr(B,best,bv,62,.049,.056,72.,.916,512,20.20);
    }else if(N<70000){
        tr(B,best,bv,84,.047,.050,58.,.932,512,19.62);
        tr(B,best,bv,74,.049,.053,66.,.923,512,19.95);
        tr(B,best,bv,63,.049,.058,76.,.916,512,20.20);
        tr(B,best,bv,54,.049,.062,82.,.912,512,20.38);
    }else if(N<180000){
        tr(B,best,bv,82,.047,.050,58.,.935,256,19.62);
        tr(B,best,bv,70,.049,.055,70.,.928,256,19.95);
        tr(B,best,bv,58,.049,.060,80.,.923,256,20.25);
    }else{
        tr(B,best,bv,80,.047,.050,58.,.940,256,19.60);
        tr(B,best,bv,68,.049,.055,70.,.934,256,19.95);
        tr(B,best,bv,56,.049,.061,82.,.928,256,20.25);
    }
    if(bv<base){
        rs(best);
        return true;
    }
    rs(B);
    return false;
}
}'''

SAMPLE_INPUT = """9 14
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
STD = set("""abort abs acos adjacent_find array atan2 begin cbrt ceil chrono clear cos count data deque duration empty end erase exit fabs fill find floor fprintf fread fwrite greater hypot insert int16_t int32_t int64_t int8_t isfinite less lower_bound make_pair map max memcpy memset min move pair pop pop_back pow priority_queue printf push push_back queue reserve resize reverse set setvbuf shrink_to_fit sin size size_t snprintf sort sqrt stable_sort stderr stdin stdout strtod strtof strtol string swap tuple uint16_t uint32_t uint64_t uint8_t unordered_map unordered_set unique upper_bound vector""".split())
ID = re.compile(r"[A-Za-z_][A-Za-z0-9_]*")
OPS3 = ("<<=", ">>=", "->*", "...")
OPS2 = ("++", "--", "->", "&&", "||", "<<", ">>", "<=", ">=", "==", "!=", "+=", "-=", "*=", "/=", "%=", "&=", "|=", "^=", "::", "##", ".*")
OPC = set("+-*&|<>=:*/.%^!#")


def die(msg: str) -> None:
    raise SystemExit("FAIL_CLOSED: " + msg)


def match_brace(s: str, o: int) -> int:
    d = 0
    i = o
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
            if c in "\"'":
                q = c
            elif c == "{":
                d += 1
            elif c == "}":
                d -= 1
                if d == 0:
                    return i + 1
        i += 1
    return -1


def find_main(s: str):
    p = s.find("int main(){")
    if p < 0:
        m = re.search(r"\bint\s+main\s*\(\s*\)\s*\{", s)
        if not m:
            die("missing int main")
        p = m.start()
    o = s.find("{", p)
    e = match_brace(s, o)
    if e < 0:
        die("unmatched main")
    return p, e


def scan(src: str):
    t = []
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
            t.append(("ws", c, pp))
            bol = True
            pp = False
            i += 1
            continue
        if c.isspace():
            j = i + 1
            while j < n and src[j].isspace() and src[j] != "\n":
                j += 1
            t.append(("ws", src[i:j], pp))
            i = j
            continue
        if pp:
            j = i + 1
            while j < n and src[j] != "\n":
                j += 1
            t.append(("pp", src[i:j], True))
            i = j
            continue
        if c in "\"'":
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
            t.append(("lit", src[st:i], False))
            continue
        if c == "/" and i + 1 < n and src[i + 1] == "/":
            j = src.find("\n", i)
            if j < 0:
                j = n
            t.append(("com", src[i:j], False))
            i = j
            continue
        if c == "/" and i + 1 < n and src[i + 1] == "*":
            j = src.find("*/", i + 2)
            if j < 0:
                die("unterminated comment")
            t.append(("com", src[i:j + 2], False))
            i = j + 2
            continue
        m = ID.match(src, i)
        if m:
            t.append(("id", m.group(0), False))
            i = m.end()
            continue
        if c.isdigit() or (c == "." and i + 1 < n and src[i + 1].isdigit()):
            j = i + 1
            while j < n and (src[j].isalnum() or src[j] in "._+-"):
                if src[j] in "+-" and not (j > i and src[j - 1] in "eEpP"):
                    break
                j += 1
            t.append(("num", src[i:j], False))
            i = j
            continue
        if i + 2 < n and src[i:i + 3] in OPS3:
            t.append(("op", src[i:i + 3], False))
            i += 3
            continue
        if i + 1 < n and src[i:i + 2] in OPS2:
            t.append(("op", src[i:i + 2], False))
            i += 2
            continue
        t.append(("op", c, False))
        i += 1
    return t


def token_names(src: str):
    return set(v for k, v, _ in scan(src) if k == "id")


def need_space(a: str, b: str) -> bool:
    if not a or not b:
        return False
    x, y = a[-1], b[0]
    if (x.isalnum() or x == "_") and (y.isalnum() or y == "_"):
        return True
    if x == "." and y == ".":
        return True
    if x in OPC and y in OPC:
        return True
    return False


def choose_macros(src: str):
    vals = [
        "double", "static", "const", "inline", "vector", "unsigned",
        "namespace", "struct", "template", "typename", "return false",
        "return true", "continue", "return", "int", "bool", "void", "class",
    ]
    used = token_names(src) | KW | STD
    names = []
    for p in ("P", "Q", "Z", "Y", "X", "W"):
        for i in range(100):
            names.append(f"{p}{i}")
    out = []
    for v in vals:
        name = None
        for cand in names:
            if cand not in used:
                name = cand
                break
        if name is None:
            die("macro pool exhausted")
        used.add(name)
        out.append((name, v))
    return out


def minify_with_macros(src: str) -> str:
    macros = choose_macros(src)
    pos = src.find("using namespace std;")
    if pos < 0:
        die("missing using namespace std")
    block = "".join(f"#define {a} {b}\n" for a, b in macros)
    src = src[:pos] + block + src[pos:]
    tok = scan(src)

    one_word = {b: a for a, b in macros if " " not in b}
    ret_false = next(a for a, b in macros if b == "return false")
    ret_true = next(a for a, b in macros if b == "return true")

    out = []
    prev = ""
    i = 0
    line = True
    while i < len(tok):
        k, v, pp = tok[i]
        if k in ("ws", "com"):
            if k == "ws" and "\n" in v:
                line = True
                prev = "\n"
            i += 1
            continue
        if k == "pp":
            if not line:
                out.append("\n")
            out.append(v.rstrip() + "\n")
            line = True
            prev = "\n"
            i += 1
            continue
        if k == "id" and v == "return":
            j = i + 1
            while j < len(tok) and tok[j][0] in ("ws", "com"):
                j += 1
            if j < len(tok) and tok[j][0] == "id" and tok[j][1] in ("false", "true"):
                cur = ret_false if tok[j][1] == "false" else ret_true
                if need_space(prev, cur):
                    out.append(" ")
                out.append(cur)
                prev = cur
                line = False
                i = j + 1
                continue
        cur = one_word.get(v, v) if k == "id" else v
        if need_space(prev, cur):
            out.append(" ")
        out.append(cur)
        prev = cur
        line = False
        i += 1
    return "".join(out)


def patch_source(src: str) -> str:
    if "namespace PCL2{" in src or "PCL2::run()" in src:
        die("PCL2 already installed")
    required = [
        "static AP AD()",
        "static void rs(",
        "static bool GD(int u,int v,const AE&params)",
        "struct AE",
        "W5::strong_validator",
        "vps(",
        "cove(",
        "originalP",
        "int main",
    ]
    for r in required:
        if r not in src:
            die("missing required anchor " + r)
    a, b = find_main(src)
    main_src = src[a:b]
    if "JD();" not in main_src:
        die("main has no JD() output anchor")
    main_new = main_src.replace("JD();", "PCL2::run();JD();")
    patched = src[:a] + LANE + main_new + src[b:]
    patched = patched.replace('if(0&&"B16P515A"){}', "")
    if "PCL2::run();JD();" not in patched:
        die("patch insertion failed")
    return patched


def build_from(inp: Path) -> str:
    src = inp.read_text(encoding="utf-8")
    patched = patch_source(src)
    out = minify_with_macros(patched)
    if len(out.encode()) > LIMIT:
        patched2 = patched.replace("if(N<50625&&es()<18.9)WK::run();", "")
        out = minify_with_macros(patched2)
    if len(out.encode()) > LIMIT:
        patched3 = patched.replace("if(N<50625&&es()<18.90)WK::run();", "")
        out = minify_with_macros(patched3)
    n = len(out.encode())
    if n > LIMIT:
        die(f"output too large {n}>{LIMIT}")
    if "PCL2::run()" not in out:
        die("PCL2 call lost after minify")
    return out


def compile_gate(cpp: Path, exe: Path):
    gpp = shutil.which("g++")
    if not gpp:
        die("g++ not found for compile gate")
    cmd = [gpp, "-std=c++17", "-O2", "-pipe", "-static", "-s", str(cpp), "-o", str(exe)]
    subprocess.run(cmd, check=True, timeout=120)
    return cmd


def sample_gate(exe: Path):
    proc = subprocess.run(
        [str(exe)],
        input=SAMPLE_INPUT.encode(),
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        timeout=20,
        check=True,
    )
    out = proc.stdout.decode(errors="replace").strip().splitlines()
    if not out:
        die("sample gate produced empty output")
    first = out[0].split()
    if len(first) != 2:
        die("sample gate first line malformed: " + out[0])
    try:
        v = int(first[0])
        f = int(first[1])
    except ValueError:
        die("sample gate first line is not integer pair: " + out[0])
    if not (1 <= v <= 9 and f > 0):
        die(f"sample gate invalid header V={v} F={f}")
    if len(out) < 1 + v + f:
        die("sample gate output has too few lines")
    return v, f


def main():
    outp = Path(sys.argv[2]) if len(sys.argv) > 2 else Path("pcl2_postfinal_submit.cpp")
    if len(sys.argv) > 1:
        tries = [Path(sys.argv[1])]
    else:
        candidates = [
            "fetched_sources/kattis_19903326_fetched.cpp",
            "submission_563_81.93_7.cpp",
            "fetched_sources/kattis_19902839_81.93_7.cpp",
            "fetched_sources/19901322.cpp",
            "fetched_sources/19901232.cpp",
            "submission_543_81.93_7.cpp",
        ]
        tries = [Path(x) for x in candidates if Path(x).exists()]
        if not tries:
            die("no base source found; pass current-best input.cpp explicitly")
    errors = []
    final = None
    used = None
    for p in tries:
        try:
            final = build_from(p)
            used = p
            break
        except SystemExit as e:
            errors.append(f"{p}: {e}")
            if len(sys.argv) > 1:
                raise
    if final is None:
        die("all candidate bases failed: " + " | ".join(errors))
    outp.write_text(final, encoding="utf-8")
    exe = outp.with_suffix("")
    cmd = compile_gate(outp, exe)
    sv, sf = sample_gate(exe)
    print(f"input={used}")
    print(f"output={outp}")
    print(f"bytes={len(final.encode())} limit={LIMIT}")
    print(f"sha256={hashlib.sha256(final.encode()).hexdigest()}")
    print("compile=" + " ".join(cmd))
    print(f"sample_ok header={sv} {sf}")
    print("route=PCL2 post-final guarded edge-collapse tournament; runs only before output, accepts only if cove() is below the already-computed current-best fallback and AF/strong_validator/vps gates pass.")


if __name__ == "__main__":
    main()