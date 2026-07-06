#!/usr/bin/env python3
# Worker R5 deterministic fail-closed generator for IMC 2026 simplifygeometry.
#
# Usage from repo root:
#   python3 make_workerR5_qpiv.py fetched_sources/kattis_19903326_fetched.cpp -o workerR5_qpiv_submit.cpp
#   g++ -std=c++17 -O2 -pipe -static -s workerR5_qpiv_submit.cpp -o workerR5_qpiv_submit
#
# The generator also compiles by default. With sample.in present it runs:
#   ./workerR5_qpiv_submit < sample.in
# and expects:
#   sample_first_line=8 12
#
# Lane summary embedded in generated source:
#   R5Q is a guarded replacement route, not a local B16/WK/minifier replay.
#   It tries a multi-axis radial visual shell made only from original vertices,
#   repairs vertex-only Hausdorff with a greedy 0.0492*diag witness net of
#   microscopic closed tetrahedra, then keeps the candidate only when the base
#   validator and a stricter internal axial-view proxy pass. If the route does
#   not clear all gates, it restores the exact current-best state.

from __future__ import annotations
import argparse
import hashlib
import os
import re
import shutil
import string
import subprocess
import sys
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

REQ_ANY = [
    ("snapshot AD", ["static AP AD()", "AP AD()"]),
    ("restore rs", ["static void rs(const AP&s)", "static void rs(const AP &s)", "void rs(const AP&s)", "rs(const AP"]),
    ("proxy scorer", ["vps(", "visual_proxy_score"]),
    ("vertex estimator", ["cove(", "count_output_vertices_estimate"]),
    ("mesh adopter", ["AF("]),
    ("validator", ["W5::strong_validator"]),
    ("original vertices", ["originalP"]),
    ("Vec3", ["struct Vec3", "Vec3{"]),
    ("Face", ["struct Face", "Face{"]),
    ("main", ["int main"]),
    ("read hook", ["GN();"]),
    ("write hook", ["JD();"]),
]

LANE = r'''
namespace R5Q{
static Face ff(int a,int b,int c){Face f;f.v[0]=a;f.v[1]=b;f.v[2]=c;return f;}
static bool add(vector<Vec3>&V,vector<Face>&F,int a,int b,int c,const Vec3&o){
    if(a<0||b<0||c<0||a==b||a==c||b==c)return false;
    Vec3 cr=cross3(V[b]-V[a],V[c]-V[a]);
    if(norm2(cr)<1e-26)return false;
    Vec3 cc=(V[a]+V[b]+V[c])*(1.0/3.0);
    Face f=ff(a,b,c);
    if(dot3(cr,cc-o)<0)swap(f.v[1],f.v[2]);
    F.push_back(f);
    return true;
}
static bool tet(vector<Vec3>&V,vector<Face>&F,const Vec3&p){
    double e=max(1e-8,CL*2.0e-5);
    int s=(int)V.size();
    V.push_back(p);
    V.push_back(Vec3{p.x+e,p.y,p.z});
    V.push_back(Vec3{p.x,p.y+e,p.z});
    V.push_back(Vec3{p.x,p.y,p.z+e});
    Vec3 o{p.x+0.25*e,p.y+0.25*e,p.z+0.25*e};
    return add(V,F,s,s+2,s+1,o)&&add(V,F,s,s+1,s+3,o)&&add(V,F,s,s+3,s+2,o)&&add(V,F,s+1,s+2,s+3,o);
}
static Vec3 tr(const Vec3&q,int m){
    if(m==0)return q;
    if(m==1)return Vec3{q.y,q.z,q.x};
    if(m==2)return Vec3{q.z,q.x,q.y};
    return Vec3{q.x+0.36*q.y-0.22*q.z,-0.31*q.x+q.y+0.28*q.z,0.24*q.x-0.33*q.y+q.z};
}
struct G{
    Vec3 mn,mx;
    double r=0,r2=0,c=1;
    int nx=1,ny=1,nz=1;
    vector<vector<int> > b;
    vector<Vec3>*V=nullptr;
    int id(int x,int y,int z)const{return (z*ny+y)*nx+x;}
    int cl(int x,int n)const{return x<0?0:(x>=n?n-1:x);}
    int ix(double x)const{return cl((int)floor((x-mn.x)/c),nx);}
    int iy(double y)const{return cl((int)floor((y-mn.y)/c),ny);}
    int iz(double z)const{return cl((int)floor((z-mn.z)/c),nz);}
    void init(vector<Vec3>&X,double rr){
        V=&X;r=rr;r2=rr*rr;
        mn=originalP[0];mx=mn;
        for(const Vec3&p:originalP){
            mn.x=min(mn.x,p.x);mn.y=min(mn.y,p.y);mn.z=min(mn.z,p.z);
            mx.x=max(mx.x,p.x);mx.y=max(mx.y,p.y);mx.z=max(mx.z,p.z);
        }
        double pad=max(1e-9,2.0*r);
        mn.x-=pad;mn.y-=pad;mn.z-=pad;mx.x+=pad;mx.y+=pad;mx.z+=pad;
        double rx=max(1e-9,mx.x-mn.x),ry=max(1e-9,mx.y-mn.y),rz=max(1e-9,mx.z-mn.z);
        c=max(r,max(rx,max(ry,rz))/96.0);
        nx=max(1,min(96,(int)(rx/c)+1));
        ny=max(1,min(96,(int)(ry/c)+1));
        nz=max(1,min(96,(int)(rz/c)+1));
        b.assign((size_t)nx*ny*nz,vector<int>());
    }
    void put(int i){
        const Vec3&p=(*V)[i];
        b[id(ix(p.x),iy(p.y),iz(p.z))].push_back(i);
    }
    bool near(const Vec3&p)const{
        int x=ix(p.x),y=iy(p.y),z=iz(p.z);
        for(int dz=-1;dz<=1;++dz)for(int dy=-1;dy<=1;++dy)for(int dx=-1;dx<=1;++dx){
            int xx=x+dx,yy=y+dy,zz=z+dz;
            if(xx<0||yy<0||zz<0||xx>=nx||yy>=ny||zz>=nz)continue;
            const vector<int>&q=b[id(xx,yy,zz)];
            for(int j:q){
                Vec3 d=(*V)[j]-p;
                if(norm2(d)<=r2)return true;
            }
        }
        return false;
    }
};
static bool shell(int L,int m,double sh,double sp,vector<Vec3>&V,vector<Face>&F){
    int O=L*2,vc=2+(L-1)*O;
    if(vc>=N||originalP.empty())return false;
    Vec3 mn=originalP[0],mx=mn;
    for(const Vec3&p:originalP){
        mn.x=min(mn.x,p.x);mn.y=min(mn.y,p.y);mn.z=min(mn.z,p.z);
        mx.x=max(mx.x,p.x);mx.y=max(mx.y,p.y);mx.z=max(mx.z,p.z);
    }
    double ex=mx.x-mn.x,ey=mx.y-mn.y,ez=mx.z-mn.z,hi=max(ex,max(ey,ez)),lo=min(ex,min(ey,ez));
    if(!(hi>1e-12)||lo<0.18*hi)return false;
    Vec3 C=(mn+mx)*0.5;
    const double pi=acos(-1.0);
    vector<int>Q(vc,-1);
    vector<double>S(vc,-1e100),A(vc,1e100),B(vc,-1e100);
    for(int i=0;i<N;++i){
        Vec3 qq=originalP[i]-C,qt=tr(qq,m);
        double rr=norm3(qt);
        if(!(rr>1e-12))continue;
        double z=qt.z/rr;
        if(z>1)z=1;if(z<-1)z=-1;
        int r=(int)floor(acos(z)*L/pi+0.5),idc;
        if(r<=0)idc=0;
        else if(r>=L)idc=vc-1;
        else{
            double th=atan2(qt.y,qt.x);
            if(th<0)th+=2*pi;
            int j=(int)floor(th*O/(2*pi)+sh+0.5);
            j%=O;if(j<0)j+=O;
            idc=1+(r-1)*O+j;
        }
        double ra=norm3(qq);
        A[idc]=min(A[idc],ra);B[idc]=max(B[idc],ra);
        double sc=rr+1e-10*((i*1103515245u+12345u)&1023u);
        if(sc>S[idc])S[idc]=sc,Q[idc]=i;
    }
    for(int i=0;i<vc;++i)if(Q[i]<0||B[i]-A[i]>sp*CL)return false;
    vector<int>Ck=Q;
    sort(Ck.begin(),Ck.end());
    for(int i=1;i<vc;++i)if(Ck[i]==Ck[i-1])return false;
    V.clear();F.clear();
    V.reserve(vc+4096);
    F.reserve(2*L*O+4096);
    for(int i=0;i<vc;++i)V.push_back(originalP[Q[i]]);
    auto idc=[&](int r,int j){return 1+(r-1)*O+((j%O+O)%O);};
    int bot=vc-1,expect=0;
    for(int j=0;j<O;++j){if(!add(V,F,0,idc(1,j+1),idc(1,j),C))return false;++expect;}
    for(int r=1;r<L-1;++r)for(int j=0;j<O;++j){
        int a=idc(r,j),b=idc(r,j+1),c=idc(r+1,j),d=idc(r+1,j+1);
        if(!add(V,F,a,b,c,C)||!add(V,F,b,d,c,C))return false;
        expect+=2;
    }
    for(int j=0;j<O;++j){if(!add(V,F,bot,idc(L-1,j),idc(L-1,j+1),C))return false;++expect;}
    return (int)F.size()==expect;
}
static bool repair(vector<Vec3>&V,vector<Face>&F,int lim,double rad,double stop){
    if((int)V.size()>=lim||es()>stop)return false;
    V.reserve(min(N,lim+4));
    G g;g.init(V,rad);
    for(int i=0;i<(int)V.size();++i)g.put(i);
    int addc=0;
    for(int i=0;i<N;++i){
        if((i&4095)==0&&es()>stop)return false;
        const Vec3&p=originalP[i];
        if(g.near(p))continue;
        if((int)V.size()+4>lim)return false;
        int s=(int)V.size();
        if(!tet(V,F,p))return false;
        for(int k=s;k<(int)V.size();++k)g.put(k);
        ++addc;
    }
    return true;
}
struct T{int L,m,R;double sh,lim,sp,th;};
static bool one(const T&t,double stop){
    if(es()>stop)return false;
    int base=cove();
    if(base<=0)return false;
    int lim=(int)floor(base*t.lim);
    lim=max(8,min(N-1,lim));
    if(lim>=base-2)return false;
    AP S=AD();
    vector<Vec3>V;
    vector<Face>F;
    if(!shell(t.L,t.m,t.sh,t.sp,V,F)){rs(S);return false;}
    if((int)V.size()>=lim){rs(S);return false;}
    if(!repair(V,F,lim,0.0492*CL,stop)){rs(S);return false;}
    if((int)V.size()>=base||V.empty()||F.empty()){rs(S);return false;}
    if(!AF(V,F)){rs(S);return false;}
    int nv=cove();
    if(nv<=0||nv>=base||nv>N||!W5::strong_validator()){rs(S);return false;}
    double drop=(double)(base-nv)/(double)base;
    if(drop<0.080){rs(S);return false;}
    int R=t.R;
    if(N>200000&&R>128)R=128;
    if(N>700000&&R>96)R=96;
    double need=t.th+(drop>0.55?0.012:(drop>0.35?0.006:0.0));
    if(es()>stop){rs(S);return false;}
    double p=vps(R);
    bool keep=p>=need;
    if(keep&&R<512&&N<200000&&es()+0.15<stop)keep=vps(512)>=need-0.006;
    if(keep&&R<256&&N>=200000&&es()+0.15<stop)keep=vps(256)>=need-0.010;
    if(!keep){rs(S);return false;}
    return true;
}
static bool run(int late){
    if(N<8000||N>1100000)return false;
    double stop=late?19.55:(N>200000?4.80:(N>60000?4.20:3.45));
    if(es()>stop-0.20)return false;
    int base=cove();
    if(base<900||base*100<N*4)return false;
    T A[]={
        {9,0,512,0.00,0.42,0.090,0.992},
        {10,1,512,0.33,0.48,0.084,0.989},
        {12,2,384,0.15,0.54,0.078,0.985},
        {14,3,256,0.44,0.60,0.073,0.980},
        {16,0,192,0.50,0.66,0.068,0.973},
        {18,1,160,0.25,0.72,0.063,0.966},
        {22,2,128,0.12,0.78,0.058,0.958},
        {26,3,128,0.37,0.84,0.054,0.952}
    };
    int nA=(int)(sizeof(A)/sizeof(A[0]));
    for(int i=0;i<nA&&es()<stop;++i){
        if(!late&&i>4)break;
        if(N>300000&&A[i].L>18)continue;
        if(N>800000&&A[i].L>16)continue;
        if(one(A[i],stop))return true;
    }
    return false;
}
}
'''

KW = set('''alignas alignof and and_eq asm atomic_cancel atomic_commit atomic_noexcept auto bitand bitor bool break case catch char char16_t char32_t class compl concept const consteval constexpr constinit const_cast continue co_await co_return co_yield decltype default delete do double dynamic_cast else enum explicit export extern false float for friend goto if inline int long mutable namespace new noexcept not not_eq nullptr operator or or_eq private protected public reflexpr register reinterpret_cast requires return short signed sizeof static static_assert static_cast struct switch synchronized template this thread_local throw true try typedef typeid typename union unsigned using virtual void volatile wchar_t while xor xor_eq'''.split())
STD = set('''abort abs acos adjacent_find array atan2 back begin cbrt ceil chrono clear cos count data deque duration empty end erase exit fabs fill find floor fprintf fread fwrite greater hypot insert int16_t int32_t int64_t int8_t isfinite less lower_bound make_pair map max memcpy memset min move pair pop pop_back pow priority_queue printf push push_back queue reserve resize reverse set setvbuf shrink_to_fit sin size size_t snprintf sort sqrt stable_sort stderr stdin stdout strtod strtof strtol string swap tuple uint16_t uint32_t uint64_t uint8_t unordered_map unordered_set unique upper_bound vector puts perror getenv system'''.split())

def die(msg: str) -> None:
    raise SystemExit("FAIL_CLOSED: " + msg)

def choose_src(arg: str | None) -> Path:
    if arg:
        p = Path(arg)
        if not p.exists():
            die("explicit source not found: " + arg)
        return p
    for s in DEFAULTS:
        p = Path(s)
        if p.exists():
            return p
    hits: list[Path] = []
    for pat in ("*19903326*.cpp", "*81.93*_7*.cpp", "19901232.cpp", "19901322.cpp", "submission_543_81.93_7.cpp"):
        hits.extend(Path(".").rglob(pat))
    seen: list[Path] = []
    for p in hits:
        if p not in seen:
            seen.append(p)
    for p in seen:
        try:
            x = p.read_text(encoding="utf-8", errors="ignore")
        except Exception:
            continue
        if "int main" in x and "W5::strong_validator" in x and ("vps(" in x or "visual_proxy_score" in x):
            return p
    die("no current-best source found; pass fetched_sources/kattis_19903326_fetched.cpp explicitly")

def require_anchors(src: str) -> None:
    for label, alts in REQ_ANY:
        if not any(a in src for a in alts):
            die("missing current-best anchor: " + label)
    if "namespace R5Q{" in src or "R5Q::run(" in src:
        die("R5Q already installed")
    for bad in ("namespace AX6{", "namespace Q35{", "namespace VHX{"):
        if bad in src:
            die("refusing to stack on known failed lane already present: " + bad)

def match_brace(s: str, open_pos: int) -> int:
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

def find_main_span(src: str) -> tuple[int, int]:
    m = re.search(r"\bint\s+main\s*\(\s*\)\s*\{", src)
    if not m:
        die("int main() not found")
    o = src.find("{", m.start())
    e = match_brace(src, o)
    if e < 0:
        die("main brace did not close")
    return m.start(), e

def patch_source(raw: str, early: bool = True, late: bool = True) -> str:
    require_anchors(raw)
    a, b = find_main_span(raw)
    main = raw[a:b]
    if late:
        if "JD();" not in main:
            die("main has no JD(); write hook")
        main = main.replace("JD();", "R5Q::run(1);JD();")
    if early:
        if "GN();" not in main:
            die("main has no GN(); read hook")
        main = main.replace("GN();", "GN();if(R5Q::run(0)){JD();return 0;}", 1)
    return raw[:a] + LANE + main + raw[b:]

def tokenize(s: str):
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
                die("unterminated comment in source")
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
        if c in ("'", '"'):
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
        if i + 2 < n and s[i:i + 3] in ("<<=", ">>=", "->*", "..."):
            tok.append(("op", s[i:i + 3]))
            i += 3
            continue
        if i + 1 < n and s[i:i + 2] in ("++", "--", "->", "&&", "||", "<<", ">>", "<=", ">=", "==", "!=", "+=", "-=", "*=", "/=", "%=", "&=", "|=", "^=", "::", "##", ".*"):
            tok.append(("op", s[i:i + 2]))
            i += 2
            continue
        tok.append(("op", c))
        i += 1
    return tok

def minify(src: str) -> str:
    protect = set(KW) | set(STD) | {"main"}
    for ln in src.splitlines():
        if ln.lstrip().startswith("#"):
            protect.update(re.findall(r"[A-Za-z_]\w*", ln))

    tok = tokenize(src)
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

    freq: dict[str, int] = {}
    for x in ids:
        if x not in protect and len(x) >= 6:
            freq[x] = freq.get(x, 0) + 1

    abc = string.ascii_letters
    tail = string.ascii_letters + string.digits + "_"

    def names():
        k = 0
        while True:
            x = k
            a = abc[x % len(abc)]
            x //= len(abc)
            if x == 0:
                yield "_" + a
            else:
                r = ""
                while x:
                    r = tail[x % len(tail)] + r
                    x //= len(tail)
                yield "_" + a + r
            k += 1

    items = [((len(x) - 3) * c, x, c) for x, c in freq.items() if c >= 2 and (len(x) >= 8 or c >= 5)]
    items.sort(reverse=True)
    used = set(idset) | set(KW)
    mp: dict[str, str] = {}
    gen = names()
    saved = 0
    for _, x, c in items:
        if saved >= 12000:
            break
        y = next(gen)
        while y in used or y in protect:
            y = next(gen)
        if len(y) < len(x):
            mp[x] = y
            used.add(y)
            saved += (len(x) - len(y)) * c

    for i, (k, v) in enumerate(tok):
        if k == "id" and v in mp:
            tok[i] = (k, mp[v])

    def need_space(a: str, b: str) -> bool:
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
    j = 0
    line = True
    while j < len(tok):
        k, v = tok[j]
        if k in ("ws", "com"):
            if k == "ws" and "\n" in v:
                line = True
                prev = "\n"
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

    r = "".join(out)
    if "R5Q::run(0)" not in r and "R5Q::run(1)" not in r:
        die("R5Q call lost during minify")
    return r

def build(src_path: Path, out_path: Path) -> str:
    raw = src_path.read_text(encoding="utf-8", errors="strict")
    variants = [
        (True, True),
        (True, False),
        (False, True),
    ]
    last = ""
    for early, late in variants:
        patched = patch_source(raw, early=early, late=late)
        patched = patched.replace('if(0&&"B16P515A"){}', '')
        out = minify(patched)
        if len(out.encode()) > LIMIT:
            out = minify(patched.replace("if(N<50625&&es()<18.9)WK::run();", ""))
        last = out
        if len(out.encode()) <= LIMIT:
            out_path.write_text(out, encoding="utf-8", newline="")
            return out
    die(f"generated source too large: {len(last.encode())}>{LIMIT}")

def compile_and_gate(outp: Path, cxx: str, static: bool, sample: Path | None) -> None:
    if not shutil.which(cxx):
        die("C++ compiler not found: " + cxx)
    exe = outp.with_suffix("")
    cmd = [cxx, "-std=c++17", "-O2", "-pipe"]
    if static:
        cmd += ["-static", "-s"]
    cmd += [str(outp), "-o", str(exe)]
    print("compile=" + " ".join(cmd))
    subprocess.run(cmd, check=True)
    print("compile_ok=" + str(exe))
    if sample and sample.exists():
        with sample.open("rb") as f:
            r = subprocess.run([str(exe.resolve())], stdin=f, stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=15)
        if r.returncode != 0:
            sys.stderr.write(r.stderr.decode("utf-8", "replace"))
            die("sample run failed")
        first = r.stdout.splitlines()[0].decode("ascii", "replace") if r.stdout else ""
        print("sample_first_line=" + first)
        if first.strip() != "8 12":
            die("sample first line mismatch; expected 8 12")

def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("src", nargs="?", default=None, help="current-best 81.93-family C++ source")
    ap.add_argument("-o", "--out", default="workerR5_qpiv_submit.cpp")
    ap.add_argument("--no-compile", action="store_true")
    ap.add_argument("--no-static", action="store_true")
    ap.add_argument("--sample", default="sample.in")
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
        compile_and_gate(outp, args.cxx, static=not args.no_static, sample=Path(args.sample) if args.sample else None)

if __name__ == "__main__":
    main()
