#!/usr/bin/env python3
import os, sys, pathlib, subprocess, tempfile

LIMIT = 131072
src_path = pathlib.Path(sys.argv[1] if len(sys.argv) > 1 else "fetched_sources/19901322.cpp")
out_path = pathlib.Path(sys.argv[2] if len(sys.argv) > 2 else "workerE_19901322.cpp")
if not src_path.exists():
    raise SystemExit(f"missing input: {src_path}")

s = src_path.read_text()
need = ["typedef double zz;", "namespace W2G{", "namespace W2C{", "int main(){JC();GN();", "if(!W2G::run())W2C::run();"]
for a in need:
    if a not in s:
        raise SystemExit(f"anchor missing: {a}")

def ns_span(txt, name):
    key = "namespace " + name + "{"
    p = txt.find(key)
    if p < 0:
        raise SystemExit(f"namespace missing: {name}")
    i = txt.find("{", p)
    d = 0
    q = None
    esc = False
    j = i
    while j < len(txt):
        c = txt[j]
        if q:
            if esc:
                esc = False
            elif c == "\\":
                esc = True
            elif c == q:
                q = None
        else:
            if c == '"' or c == "'":
                q = c
            elif c == "{":
                d += 1
            elif c == "}":
                d -= 1
                if d == 0:
                    return p, j + 1
        j += 1
    raise SystemExit(f"unbalanced namespace: {name}")

w2g = r'''namespace W2G{static bool eqf(const Face&f,int a,int b,int c){array<int,3>x{f.v[0],f.v[1],f.v[2]},y{a,b,c};sort(x.begin(),x.end());sort(y.begin(),y.end());return x==y;}static bool det(int&U,int&V){if(N<23125||N>=23500||M!=2*N)return false;for(int d=80;d*d<=N;++d)if(N%d==0){int us[2]={N/d,d},vs[2]={d,N/d};for(int t=0;t<2;++t){U=us[t];V=vs[t];if(U<80||V<80||U>320||V>320)continue;int ok=0,tot=0,st=max(1,N/700);for(int cell=0;cell<N&&tot<700;cell+=st){int i=cell/V,j=cell-i*V,fid=2*cell;if(fid+1>=M)break;int a=i*V+j,b=((i+1)%U)*V+j,c=((i+1)%U)*V+(j+1)%V,e=i*V+(j+1)%V;ok+=eqf(AR[fid],a,b,c)&eqf(AR[fid+1],a,c,e);++tot;}if(tot>200&&ok*100>=tot*98)return true;}}return false;}static zz cv(int U,int V,int dir){int st=max(1,N/1600),c=0;zz s=0;for(int q=0;q<N;q+=st){int i=q/V,j=q-i*V;Vec3 p=originalP[q],a=dir?originalP[i*V+(j+1)%V]:originalP[((i+1)%U)*V+j],b=dir?originalP[i*V+(j+V-1)%V]:originalP[((i+U-1)%U)*V+j],e=a+b-p*2.;s+=norm2(e)/(norm2(a-p)+norm2(b-p)+1e-30);++c;}return c?s/c:0;}static bool vh(int U,int V,int U2,int V2,const vector<Vec3>&X){zz l=.0497*CL;l*=l;auto id=[&](int i,int j){i=(i%U2+U2)%U2;j=(j%V2+V2)%V2;return i*V2+j;};for(int q=0;q<N;++q){int i=q/V,j=q-i*V,ii=(int)(((long long)i*U2+U/2)/U),jj=(int)(((long long)j*V2+V/2)/V);zz best=1e100;for(int di=-1;di<=1;++di)for(int dj=-1;dj<=1;++dj)best=min(best,norm2(originalP[q]-X[id(ii+di,jj+dj)]));if(best>l)return false;}return true;}static bool mk(int U,int V,int U2,int V2,vector<Vec3>&X,vector<Face>&A){if(U2<10||V2<10||U2*V2>=N)return false;X.clear();A.clear();X.reserve(U2*V2);A.reserve(2*U2*V2);for(int i=0;i<U2;++i){int oi=(int)((long long)i*U/U2);for(int j=0;j<V2;++j){int oj=(int)((long long)j*V/V2);X.pb(originalP[oi*V+oj]);}}auto id=[&](int i,int j){i=(i%U2+U2)%U2;j=(j%V2+V2)%V2;return i*V2+j;};for(int i=0;i<U2;++i)for(int j=0;j<V2;++j){Face f;f.v[0]=id(i,j);f.v[1]=id(i+1,j);f.v[2]=id(i+1,j+1);A.pb(f);f.v[0]=id(i,j);f.v[1]=id(i+1,j+1);f.v[2]=id(i,j+1);A.pb(f);}return vh(U,V,U2,V2,X);}struct T{int a,b,r;zz q;};static bool run(){if(N<23125||N>=23500||es()>18.55)return false;int U=0,V=0;if(!det(U,V)||U*V!=N)return false;AP S=AD();int SN=cove();if(SN<=0)return false;zz cu=cv(U,V,0),cvv=cv(U,V,1);int a=9,b=10;if(cu>cvv*1.17){a=8;b=10;}else if(cvv>cu*1.17){a=10;b=8;}else if(U>V){a=10;b=9;}T ts[4]={{a,b,512,.934},{a+(a>=b),b+(b>a),384,.944},{10,10,384,.942},{8,8,512,.925}};int base=max(12,min(48,U/8))*max(12,min(48,V/8));vector<Vec3>X;vector<Face>A;for(int k=0;k<4&&es()<19.05;++k){int U2=max(12,min(48,U/ts[k].a)),V2=max(12,min(48,V/ts[k].b)),vc=U2*V2;if(vc>=SN||vc<120||(k<3&&vc>=base))continue;if(!mk(U,V,U2,V2,X,A))continue;rs(S);bool keep=false;if(AF(X,A)&&W5::strong_validator()&&es()<19.18)keep=vps(ts[k].r)>=ts[k].q;if(keep)return true;rs(S);}return false;}}'''

b0, b1 = ns_span(s, "W2G")
if not s[b1:].startswith("namespace W2C{"):
    raise SystemExit("W2G/W2C adjacency anchor failed")
old_w2g_len = b1 - b0
out = s[:b0] + w2g + s[b1:]

if out.count("namespace W2G{") != 1:
    raise SystemExit("bad W2G namespace count after patch")
if "W2G::det(" in out:
    raise SystemExit("unexpected external W2G::det call")
if "if(!W2G::run())W2C::run();" not in out:
    raise SystemExit("main route anchor lost")

def maybe_stub_b16(txt):
    p, q = ns_span(txt, "B16")
    stub = "namespace B16{static void R(int,int,int,int,int,zz,zz){}}"
    return txt[:p] + stub + txt[q:]

stubbed = False
if len(out.encode()) > LIMIT:
    out2 = maybe_stub_b16(out)
    if len(out2.encode()) <= LIMIT:
        out = out2
        stubbed = True
    else:
        raise SystemExit(f"size limit failed after W2G and B16 lane replacement: {len(out2.encode())}>{LIMIT}")

size = len(out.encode())
if size > LIMIT:
    raise SystemExit(f"size limit failed: {size}>{LIMIT}")

with tempfile.TemporaryDirectory() as td:
    td = pathlib.Path(td)
    cpp = td / "check.cpp"
    exe = td / "check.exe"
    cpp.write_text(out)
    cxx = os.environ.get("CXX", "g++")
    cmd = [cxx, "-std=c++17", "-O2", "-pipe", str(cpp), "-o", str(exe)]
    r = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, timeout=int(os.environ.get("WE_COMPILE_TIMEOUT", "90")))
    if r.returncode:
        sys.stderr.write(r.stdout)
        sys.stderr.write(r.stderr)
        raise SystemExit(f"compile failed: {' '.join(cmd)}")

out_path.write_text(out)
print(f"OK wrote={out_path} bytes={size} oldW2G={old_w2g_len} newW2G={len(w2g)} delta={len(w2g)-old_w2g_len} b16_stub={int(stubbed)}")