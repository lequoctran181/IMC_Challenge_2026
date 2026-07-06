#!/usr/bin/env python3
from pathlib import Path
import re, sys, hashlib

SRC = ["submission_543_81.93_7.cpp", "fetched_sources/19901232.cpp"]
OUT = "worker_breakthrough_BOX3_from_543.cpp"
LIMIT = 131072

MAIN_RE = re.compile(
    r'int main\(\)\{JC\(\);GN\(\);'
    r'if\(!W2G::run\(\)\)W2C::run\(\);'
    r'W5::post_patch_pass\(\);VIMP::run\(\);MIDEC::run\(\);WK::run\(\);'
    r'B16::R\(39000,60000,220,-7,192,\.96,18\.05\);'
    r'if\(0&&"B16P515A"\)\{\}'
    r'B16::R\(39000,60000,76,-10,192,\.96(?:1)?,18\.35\);'
    r'for\(int i=2;i--&&N>47500&N<6e4&&es\(\)<18\.6;\)WK::run\(\);'
    r'if\(N<50625&&es\(\)<18\.9\)WK::run\(\);JD\(\);\}'
)

BOX3 = r'''namespace B3{bool run(){if(N<800||es()>17.25)return 0;Vec3 mn=originalP[0],mx=mn;for(const Vec3&p:originalP){mn.x=min(mn.x,p.x);mn.y=min(mn.y,p.y);mn.z=min(mn.z,p.z);mx.x=max(mx.x,p.x);mx.y=max(mx.y,p.y);mx.z=max(mx.z,p.z);}if(min(mx.x-mn.x,min(mx.y-mn.y,mx.z-mn.z))<=.03*CL)return 0;double L=.0035*CL,sm=0;int h[6]={},bad=0,tot=0,st=max(1,N/160000);for(int i=0;i<N;i+=st){const Vec3&p=originalP[i];double d[6]={fabs(p.x-mn.x),fabs(p.x-mx.x),fabs(p.y-mn.y),fabs(p.y-mx.y),fabs(p.z-mn.z),fabs(p.z-mx.z)};int b=0;for(int j=1;j<6;++j)if(d[j]<d[b])b=j;if(d[b]>L)++bad;++h[b];sm+=d[b];++tot;}if(tot<80||bad*1000>tot*3||sm>L*.18*tot)return 0;for(int i=0;i<6;++i)if(h[i]<max(2,tot/3000))return 0;vector<Vec3>V={{mn.x,mn.y,mn.z},{mx.x,mn.y,mn.z},{mx.x,mx.y,mn.z},{mn.x,mx.y,mn.z},{mn.x,mn.y,mx.z},{mx.x,mn.y,mx.z},{mx.x,mx.y,mx.z},{mn.x,mx.y,mx.z}};vector<Face>F;int t[12][3]={{0,1,2},{0,2,3},{4,6,5},{4,7,6},{0,4,5},{0,5,1},{1,5,6},{1,6,2},{2,6,7},{2,7,3},{3,7,4},{3,4,0}};Vec3 C=(mn+mx)*.5;double sg=BC(C);for(int i=0;i<12;++i){Face f;f.v[0]=t[i][0];f.v[1]=t[i][1];f.v[2]=t[i][2];BD(V,f,C,sg);F.pb(f);}AP S=AD();bool ok=cove()>8&&AF(V,F)&&es()<18.05&&vps(N>140000?256:512)>=(N>140000?.970:.982);if(ok)return 1;rs(S);return 0;}}
'''

def die(msg):
    print("FAIL:", msg, file=sys.stderr)
    sys.exit(2)

src_name = None
src = None
for s in SRC:
    p = Path(s)
    if p.exists():
        src_name = s
        src = p.read_text()
        break
if src is None:
    die("missing current-best source: submission_543_81.93_7.cpp or fetched_sources/19901232.cpp")

in_bytes = len(src.encode())
if not (128000 <= in_bytes <= LIMIT):
    die(f"unexpected source size {in_bytes}; refusing stale/non-current source")

required = [
    "W5::post_patch_pass();VIMP::run();MIDEC::run();WK::run();",
    "B16::R(39000,60000,220,-7,192,.96,18.05);",
    'if(0&&"B16P515A"){}',
    "for(int i=2;i--&&N>47500&N<6e4&&es()<18.6;)WK::run();",
    "if(N<50625&&es()<18.9)WK::run();JD();",
    "namespace B16",
    "namespace WK",
    "static bool HS(vector<int>&EC,vector<Face>&AB)",
]
for x in required:
    if x not in src:
        die("missing required 543 anchor: " + x)

if not re.search(r"B16::R\(39000,60000,76,-10,192,\.96(?:1)?,18\.35\);", src):
    die("missing 543 second B16 pass 76,-10")
if "-13" in src or "B01" in src or "B02" in src:
    die("forbidden stale/crash token present")
if "BOX3::run()" in src:
    die("already patched")

m = MAIN_RE.search(src)
if not m or len(MAIN_RE.findall(src)) != 1:
    die("exact main anchor mismatch")

old_main = m.group(0)
prefix = "int main(){JC();GN();"
suffix = "JD();}"
body = old_main[len(prefix):-len(suffix)]
new_main = prefix + "if(!B3::run()){" + body + "}JD();}"

patched = src[:m.start()] + BOX3 + new_main + src[m.end():]

# behavior-preserving byte recovery. Keep true/false literals intact: replacing them
# can break C++ lambda return-type deduction when a lambda also returns integer 0.
patched = re.sub(r"\bnullptr\b", "0", patched)
patched = re.sub(r"\binline\s+", "", patched)
patched = patched.replace('if(0&&"B16P515A"){}', "")
patched = patched.replace("0.0", ".0").replace("1.0", "1.")

out_bytes = len(patched.encode())
if "B3::run()" not in patched:
    die("patch insertion failed")
if out_bytes > LIMIT:
    die(f"output too large: {out_bytes} > {LIMIT}")
for x in ["39000,60000,76,-10,192,.96", "N>47500&N<6e4", "N<50625", "B3"]:
    if x not in patched:
        die("post-patch validation lost token: " + x)

Path(OUT).write_text(patched)
print(f"source={src_name}")
print(f"output={OUT}")
print(f"input_bytes={in_bytes}")
print(f"output_bytes={out_bytes}")
print(f"freed_or_spent={in_bytes-out_bytes}")
print(f"sha256={hashlib.sha256(patched.encode()).hexdigest()}")
print(f"compile_command=g++ -std=c++17 -O2 -pipe -static -s {OUT} -o worker_breakthrough_BOX3_from_543")
print("score_signal=non-box hidden set should tie 81.934570/7; bbox/cuboid structural family should collapse to 8v/12f and move materially upward; validator/proxy rollback preserves 7/7 on false positives")
