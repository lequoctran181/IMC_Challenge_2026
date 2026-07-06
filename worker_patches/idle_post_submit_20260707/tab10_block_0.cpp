#!/usr/bin/env python3
# R3 deterministic generator.
# Usage:
#   python3 gen_r3.py fetched_sources/19901232.cpp r3.cpp
#   g++ -std=gnu++17 -O2 -pipe -static -s r3.cpp -o r3
# Sample gate expectation:
#   generator fails closed on missing anchors/oversize; generated source should compile cleanly.
#   On ordinary samples, fallback preserves the base accepted output; on hidden, R3 accepts only AF+strong_validator+vps-gated remeshes.

from pathlib import Path
import sys, re, hashlib

LIMIT = 131072
BASE_CANDIDATES = [
    "submission_563_81.93_7.cpp",
    "submission_543_81.93_7.cpp",
    "submission_576_81.93_7.cpp",
    "submission_580_81.93_7.cpp",
    "submission_585_81.93_7.cpp",
    "submission_597_81.93_7.cpp",
    "submission_608_81.93_7.cpp",
    "fetched_sources/19901232.cpp",
    "fetched_sources/19901322.cpp",
    "fetched_sources/19901342.cpp",
    "fetched_sources/kattis_19902206_81.93_7.cpp",
    "fetched_sources/kattis_19902388_81.93_7.cpp",
    "fetched_sources/kattis_19902774_81.93_7_external_tie.cpp",
]

LANE = r'''
namespace R3{struct QK{long long k;int v;};static Face mf(int a,int b,int c){Face f;f.v[0]=a;f.v[1]=b;f.v[2]=c;return f;}static Vec3 mnv,mxv,cen;static double dx,dy,dz,sc;static inline long long pack(int a,int b,int c,int s){return(((long long)(a+2048)&4095)<<52)|(((long long)(b+2048)&4095)<<40)|(((long long)(c+2048)&4095)<<28)|s;}static inline int side(const Vec3&n){double ax=fabs(n.x),ay=fabs(n.y),az=fabs(n.z);if(ax>=ay&&ax>=az)return n.x>=0?0:1;if(ay>=az)return n.y>=0?2:3;return n.z>=0?4:5;}static bool ck(vector<Vec3>&X,vector<Face>&F,double q){if(X.size()<4||F.size()<4)return 0;AP S=AD();if(AF(X,F)&&W5::strong_validator()&&vps(512)>=q)return 1;rs(S);for(auto&f:F)swap(f.v[1],f.v[2]);if(AF(X,F)&&W5::strong_validator()&&vps(512)>=q)return 1;rs(S);return 0;}static void bbox(){mnv={1e100,1e100,1e100};mxv={-1e100,-1e100,-1e100};for(auto&p:originalP){mnv.x=min(mnv.x,p.x);mnv.y=min(mnv.y,p.y);mnv.z=min(mnv.z,p.z);mxv.x=max(mxv.x,p.x);mxv.y=max(mxv.y,p.y);mxv.z=max(mxv.z,p.z);}cen=(mnv+mxv)*.5;dx=mxv.x-mnv.x;dy=mxv.y-mnv.y;dz=mxv.z-mnv.z;sc=max(dx,max(dy,dz));if(!(sc>0))sc=1;}static bool grid(int g,int split,double q,int cap){if(es()>17.8||N<1200)return 0;bbox();vector<Vec3>vn(N,{0,0,0});for(auto&f:AR){Vec3 n=cross3(originalP[f.v[1]]-originalP[f.v[0]],originalP[f.v[2]]-originalP[f.v[0]]);for(int j=0;j<3;j++)vn[f.v[j]]=vn[f.v[j]]+n;}vector<QK>A;A.reserve(N);double inv=(g-1)/sc;for(int i=0;i<N;i++){Vec3 p=originalP[i]-mnv;int ix=max(0,min(g-1,(int)(p.x*inv+.5))),iy=max(0,min(g-1,(int)(p.y*inv+.5))),iz=max(0,min(g-1,(int)(p.z*inv+.5)));int s=split?side(vn[i]):0;A.push_back({pack(ix,iy,iz,s),i});}sort(A.begin(),A.end(),[](const QK&a,const QK&b){return a.k<b.k;});vector<int>mp(N,-1),rep;vector<Vec3>X;for(size_t i=0;i<A.size();){size_t j=i+1;while(j<A.size()&&A[j].k==A[i].k)j++;Vec3 c{0,0,0};for(size_t t=i;t<j;t++)c=c+originalP[A[t].v];c=c*(1.0/(j-i));int best=A[i].v;double bd=1e100;for(size_t t=i;t<j;t++){double d=norm2(originalP[A[t].v]-c);if(d<bd){bd=d;best=A[t].v;}}int id=(int)X.size();X.push_back(originalP[best]);for(size_t t=i;t<j;t++)mp[A[t].v]=id;i=j;}if((int)X.size()>cap||X.size()<4)return 0;vector<Face>F;F.reserve(M);vector<unsigned long long>seen;seen.reserve(M);auto triKey=[&](int a,int b,int c){int x[3]={a,b,c};sort(x,x+3);return((unsigned long long)x[0]<<42)^((unsigned long long)x[1]<<21)^(unsigned long long)x[2];};for(auto&of:AR){int a=mp[of.v[0]],b=mp[of.v[1]],c=mp[of.v[2]];if(a==b||b==c||a==c)continue;unsigned long long k=triKey(a,b,c);seen.push_back(k);F.push_back(mf(a,b,c));}if(F.size()<X.size()*2/3)return 0;vector<int>ord(F.size());for(int i=0;i<(int)ord.size();++i)ord[i]=i;sort(ord.begin(),ord.end(),[&](int a,int b){return seen[a]<seen[b];});vector<Face>G;G.reserve(F.size());for(int ii=0;ii<(int)ord.size();){int jj=ii+1;while(jj<(int)ord.size()&&seen[ord[jj]]==seen[ord[ii]])jj++;if(jj==ii+1)G.push_back(F[ord[ii]]);ii=jj;}if(G.size()<X.size())return 0;return ck(X,G,q);}static bool shell(int u,int v,double q){if(es()>17.8||N<5000)return 0;bbox();vector<Vec3>X;vector<Face>F;auto add=[&](Vec3 p){X.push_back(p);return(int)X.size()-1;};auto P0=[&](double a,double b,int ax,int sg){Vec3 p=cen;double A=a*2-1,B=b*2-1;if(ax==0){p.x=sg>0?mxv.x:mnv.x;p.y=cen.y+A*dy*.5;p.z=cen.z+B*dz*.5;}if(ax==1){p.y=sg>0?mxv.y:mnv.y;p.x=cen.x+A*dx*.5;p.z=cen.z+B*dz*.5;}if(ax==2){p.z=sg>0?mxv.z:mnv.z;p.x=cen.x+A*dx*.5;p.y=cen.y+B*dy*.5;}return p;};int base[6];vector<vector<int>>id(6,vector<int>((u+1)*(v+1)));for(int s=0;s<6;s++){int ax=s/2,sg=s%2?1:-1;for(int i=0;i<=u;i++)for(int j=0;j<=v;j++)id[s][i*(v+1)+j]=add(P0((double)i/u,(double)j/v,ax,sg));}auto at=[&](int s,int i,int j){return id[s][i*(v+1)+j];};for(int s=0;s<6;s++){for(int i=0;i<u;i++)for(int j=0;j<v;j++){int a=at(s,i,j),b=at(s,i+1,j),c=at(s,i+1,j+1),d=at(s,i,j+1);if(s%2)F.push_back(mf(a,b,c)),F.push_back(mf(a,c,d));else F.push_back(mf(a,c,b)),F.push_back(mf(a,d,c));}}return ck(X,F,q);}static bool run(){if(es()>16.5)return 0;int n=N;int c8=max(64,n*8/100),c10=max(96,n*10/100),c12=max(128,n*12/100),c15=max(160,n*15/100);if(n>=18000){if(grid(18,1,.955,c8))return 1;if(grid(20,1,.960,c10))return 1;if(grid(22,1,.965,c12))return 1;if(grid(24,1,.970,c15))return 1;}if(n>=5000&&n<18000){if(grid(16,1,.955,c10))return 1;if(grid(18,1,.965,c12))return 1;if(grid(20,1,.970,c15))return 1;}if(n>=30000){if(shell(18,18,.975))return 1;if(shell(24,16,.980))return 1;}return 0;}}
'''

MACROS = [
    ("O0","double"),("O1","static"),("O2","const"),("O3","inline"),
    ("O4","vector"),("O5","unsigned"),("O6","namespace"),("O7","struct"),
    ("O8","template"),("O9","typename"),("P0","return false"),
    ("P1","return true"),("P2","continue"),("P3","return"),
    ("P4","int"),("P5","bool"),("P6","void"),("P7","class")
]
ID = re.compile(r"[A-Za-z_][A-Za-z0-9_]*")

def die(s):
    raise SystemExit("R3_FAIL_CLOSED: " + s)

def toks(src):
    out=[];i=0;n=len(src);bol=True;pp=False
    while i<n:
        c=src[i]
        if bol:
            j=i
            while j<n and src[j] in " \t": j+=1
            pp=j<n and src[j]=="#"; bol=False
        if c=="\n": bol=True; pp=False; i+=1; continue
        if pp: i+=1; continue
        if c in "\"'":
            q=c; i+=1
            while i<n:
                if src[i]=="\\": i+=2; continue
                if src[i]==q: i+=1; break
                i+=1
            continue
        if c=="/" and i+1<n and src[i+1]=="/":
            j=src.find("\n",i); i=n if j<0 else j; continue
        if c=="/" and i+1<n and src[i+1]=="*":
            j=src.find("*/",i+2); i=n if j<0 else j+2; continue
        m=ID.match(src,i)
        if m: out.append(m.group(0)); i=m.end(); continue
        i+=1
    return out

def macroize(src):
    if "#define O0 double" in src:
        return src
    used=set(toks(src))
    coll=[a for a,_ in MACROS if a in used]
    if coll:
        return src
    u=src.find("using namespace std;")
    if u<0:
        u=src.find("using O6 std;")
        if u<0: return src
    defs="".join("#define %s %s\n"%(a,b) for a,b in MACROS)
    src=src[:u]+defs+src[u:]
    mp={b:a for a,b in MACROS if " " not in b}
    out=[];i=0;n=len(src);bol=True;pp=False
    while i<n:
        c=src[i]
        if bol:
            j=i
            while j<n and src[j] in " \t": j+=1
            pp=j<n and src[j]=="#"; bol=False
        if c=="\n": out.append(c); bol=True; pp=False; i+=1; continue
        if pp: out.append(c); i+=1; continue
        if c in "\"'":
            q=c; st=i; i+=1
            while i<n:
                if src[i]=="\\": i+=2; continue
                if src[i]==q: i+=1; break
                i+=1
            out.append(src[st:i]); continue
        if c=="/" and i+1<n and src[i+1]=="/":
            j=src.find("\n",i); j=n if j<0 else j; out.append(src[i:j]); i=j; continue
        if c=="/" and i+1<n and src[i+1]=="*":
            j=src.find("*/",i+2); j=n if j<0 else j+2; out.append(src[i:j]); i=j; continue
        m=ID.match(src,i)
        if m:
            w=m.group(0)
            if w=="return":
                k=m.end()
                while k<n and src[k] in " \t\r\n": k+=1
                m2=ID.match(src,k)
                if m2 and m2.group(0)=="false":
                    out.append("P0"); i=m2.end(); continue
                if m2 and m2.group(0)=="true":
                    out.append("P1"); i=m2.end(); continue
            out.append(mp.get(w,w)); i=m.end(); continue
        out.append(c); i+=1
    return "".join(out)

def patch(src):
    need = ["static AP AD()","static void rs(","AF(","W5::strong_validator","vps(","int main(){JC();"]
    for t in need:
        if t not in src:
            die("missing required anchor " + repr(t))
    if "namespace R3{" in src:
        die("R3 already present")
    src = src.replace("int main(){JC();", LANE + "int main(){JC();", 1)
    if "GN();if(!W2G::run())" in src:
        src = src.replace("GN();if(!W2G::run())", "GN();if(!R3::run())if(!W2G::run())", 1)
    elif "GN();W2G::run();" in src:
        src = src.replace("GN();W2G::run();", "GN();if(!R3::run())W2G::run();", 1)
    elif "GN();" in src:
        src = src.replace("GN();", "GN();R3::run();", 1)
    else:
        die("missing GN call anchor")
    return src

def main():
    inp = Path(sys.argv[1]) if len(sys.argv) > 1 else next((Path(p) for p in BASE_CANDIDATES if Path(p).exists()), None)
    if inp is None or not inp.exists():
        die("missing base source; pass current accepted 81.93-family cpp as argv[1]")
    outp = Path(sys.argv[2]) if len(sys.argv) > 2 else Path("r3.cpp")
    raw = inp.read_text()
    base_len = len(raw.encode())
    s = patch(raw)
    if len(s.encode()) > LIMIT:
        s = macroize(s)
    size = len(s.encode())
    if size > LIMIT:
        die("generated source too large: %d > %d; base=%d delta=%d" % (size, LIMIT, base_len, size-base_len))
    if "R3::run()" not in s or "namespace R3{" not in s:
        die("patch lost R3")
    outp.write_text(s)
    print("wrote=%s bytes=%d base=%d delta=%d sha256=%s" % (outp, size, base_len, size-base_len, hashlib.sha256(s.encode()).hexdigest()))
    print("compile: g++ -std=gnu++17 -O2 -pipe -static -s %s -o r3" % outp)
    print("signal: fallback-preserving; R3 only commits AF+W5::strong_validator+vps512-gated grid-normal remesh/shell candidates")

if __name__ == "__main__":
    main()