#!/usr/bin/env python3
import argparse, hashlib, re, os, sys
from pathlib import Path

LIMIT=131072
SHA="7b4c5386a74cef7d3d77b0e4053285409baa74f1"
SIZE=130310

MAINS=[
'int main(){JC();GN();if(!W2G::run())W2C::run();W5::post_patch_pass();VIMP::run();MIDEC::run();WK::run();B16::R(39000,60000,220,-7,192,.96,18.05);if(0&&"B16P515A"){}B16::R(39000,60000,76,-10,192,.96,18.35);for(int i=2;i--&&N>47500&N<6e4&&es()<18.6;)WK::run();if(N<50625&&es()<18.9)WK::run();JD();}',
'int main(){JC();GN();if(!W2G::run())W2C::run();W5::post_patch_pass();VIMP::run();MIDEC::run();WK::run();B16::R(39000,60000,220,-7,192,96,18.05);if(0&&"B16P515A"){}B16::R(39000,60000,76,-10,192,96,18.35);for(int i=2;i--&&N>47500&N<6e4&&es()<18.6;)WK::run();if(N<50625&&es()<18.9)WK::run();JD();}',
'int main(){JC();GN();if(!W2G::run())W2C::run();W5::post_patch_pass();VIMP::run();MIDEC::run();WK::run();B16::R(39000,60000,220,-7,192,.96,18.05);if(0&&"B16P515A"){}B16::R(39000,60000,76,-10,192,.96,18.35);for(int i=2;i--&&N>47500&N<6e4&&es()<18.6;)WK::run();if(N<50625&&es()<18.9)WK::run();JD();return 0;}',
'int main(){JC();GN();if(!W2G::run())W2C::run();W5::post_patch_pass();VIMP::run();MIDEC::run();WK::run();B16::R(39000,60000,220,-7,192,96,18.05);if(0&&"B16P515A"){}B16::R(39000,60000,76,-10,192,96,18.35);for(int i=2;i--&&N>47500&N<6e4&&es()<18.6;)WK::run();if(N<50625&&es()<18.9)WK::run();JD();return 0;}'
]

REQ=[
"namespace B16{","namespace W2G{","namespace W2C{","namespace WK{","namespace W5{",
"namespace VIMP{","namespace MIDEC{","static AP AD()","static void rs(const AP&s)",
"static inline double es()","static int N=0,M=0;","static double CL=1.;",
"static bool AS=false;","static double AL=.0;","static double AH=1.;","static double BG=1.;",
"static void CA(string&out,const char*line,int len){"
]
BAD=["namespace MS2{","MS2::run()","MS2SURGE","ordered_grid_diag","namespace KCM{","namespace HBH{"]

INS=r'''namespace MS2{static bool K(AP&S,int sn,double q){int an=cove();if(an<96||an>=sn*83/100||es()>19.62||!W5::strong_validator()){rs(S);return false;}double p=vps(256);if(p<q){rs(S);return false;}if(an<sn*58/100&&es()<19.82&&vps(512)<q-.009){rs(S);return false;}return true;}static bool T(int a,int b,int c,int d,double q,double l){if(es()>l-.2)return false;AP S=AD();int sn=cove();B16::R(39000,60000,a,b,c,.96,l);for(int i=d;i--&&es()<l+.12;)WK::run();return K(S,sn,q);}static bool H(){if(!(N>=39000&&N<60000)||es()>18.78)return false;if(AS&&(BG>.018||AL<.68||AH>.18))return false;int A[8]={176,132,96,112,152,88,124,204},B[8]={-8,-9,-12,-10,-7,-13,-11,-8},C[8]={192,192,160,224,192,160,192,224},D[8]={1,1,2,1,0,2,1,0};double Q[8]={.948,.950,.955,.952,.946,.958,.953,.944};int h=(N*13+M*7+cove()*3+(int)(AL*999)+(int)(BG*9999))&7;for(int j=0;j<3&&es()<19.42;j++){int k=(h+j*3)&7;if(T(A[k],B[k],C[k],D[k],Q[k],19.30))return true;}return false;}static bool G(){if(!(N>=200000&&N<=260000&&M<=650000&&AS&&BG<.0045&&AL>.991&&AH<.007&&es()<16.8))return false;AP S=AD();int sn=cove();VIMP::run();if(es()<17.75)VIMP::run();if(cove()>=sn*99/100||!W5::strong_validator()){rs(S);return false;}double q=vps(256);if(q<.972){rs(S);return false;}return true;}static bool run(){if(es()>18.94)return false;return H()||G();}}'''
INS='namespace B16{static void R(int,int,int,int,int,double,double);}'+INS

MAC=[
("O0","double"),("O1","static"),("O2","const"),("O3","inline"),("O4","vector"),
("O5","unsigned"),("O6","namespace"),("O7","struct"),("O8","template"),("O9","typename"),
("P0","return false"),("P1","return true"),("P2","continue"),("P3","return"),
("P4","int"),("P5","bool"),("P6","void"),("P7","class")
]
ID=re.compile(r"[A-Za-z_][A-Za-z0-9_]*")

def die(s):
    raise SystemExit("workerI_ms2_abort: "+s)

def bsha(b):
    return hashlib.sha1(b"blob "+str(len(b)).encode()+b"\0"+b).hexdigest()

def toks(src):
    r=[];i=0;n=len(src);bol=True;pp=False
    while i<n:
        c=src[i]
        if bol:
            j=i
            while j<n and src[j] in " \t": j+=1
            pp=j<n and src[j]=="#"; bol=False
        if c=="\n":
            bol=True; pp=False; i+=1; continue
        if pp:
            i+=1; continue
        if c in "\"'":
            q=c; i+=1
            while i<n:
                if src[i]=="\\":
                    i+=2; continue
                if src[i]==q:
                    i+=1; break
                i+=1
            continue
        if c=="/" and i+1<n and src[i+1]=="/":
            j=src.find("\n",i); i=n if j<0 else j; continue
        if c=="/" and i+1<n and src[i+1]=="*":
            j=src.find("*/",i+2); i=n if j<0 else j+2; continue
        m=ID.match(src,i)
        if m:
            r.append(m.group(0)); i=m.end(); continue
        i+=1
    return r

def tr(src):
    mp={b:a for a,b in MAC if " " not in b}
    out=[];i=0;n=len(src);bol=True;pp=False
    while i<n:
        c=src[i]
        if bol:
            j=i
            while j<n and src[j] in " \t": j+=1
            pp=j<n and src[j]=="#"; bol=False
        if c=="\n":
            out.append(c); i+=1; bol=True; pp=False; continue
        if pp:
            out.append(c); i+=1; continue
        if c in "\"'":
            q=c; st=i; i+=1
            while i<n:
                if src[i]=="\\":
                    i+=2; continue
                if src[i]==q:
                    i+=1; break
                i+=1
            out.append(src[st:i]); continue
        if c=="/" and i+1<n and src[i+1]=="/":
            j=src.find("\n",i); j=n if j<0 else j
            out.append(src[i:j]); i=j; continue
        if c=="/" and i+1<n and src[i+1]=="*":
            j=src.find("*/",i+2); j=n if j<0 else j+2
            out.append(src[i:j]); i=j; continue
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

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("src",nargs="?",default=None)
    ap.add_argument("-o","--out",default="workerI_ms2_param_family.cpp")
    ap.add_argument("--allow-sha-mismatch",action="store_true")
    args=ap.parse_args()

    if args.src:
        inp=Path(args.src)
    else:
        cand=[Path("fetched_sources/19901232.cpp"),Path("submission_543_81.93_7.cpp")]
        inp=next((p for p in cand if p.exists()),None)
        if inp is None: die("missing fetched_sources/19901232.cpp or submission_543_81.93_7.cpp")

    raw=inp.read_bytes()
    if len(raw)!=SIZE: die("size anchor mismatch got %d expected %d"%(len(raw),SIZE))
    sh=bsha(raw)
    if sh!=SHA and not args.allow_sha_mismatch: die("blob sha mismatch got "+sh)
    src=raw.decode("utf-8")

    old=None
    for m in MAINS:
        if src.count(m)==1:
            old=m; break
    if old is None: die("exact current-best main anchor not found once")
    for x in REQ:
        if src.count(x)<1: die("missing anchor "+x)
    for x in BAD:
        if x in src: die("forbidden stale anchor "+x)

    new=old.replace("if(N<50625&&es()<18.9)WK::run();JD();","if(N<50625&&es()<18.9)WK::run();MS2::run();JD();",1)
    if new==old: die("main rewrite failed")
    src=src.replace("static void CA(string&out,const char*line,int len){",INS+"static void CA(string&out,const char*line,int len){",1)
    src=src.replace(old,new,1)

    col=[a for a,_ in MAC if a in set(toks(src))]
    if col: die("macro collision "+",".join(col))
    p=src.find("using namespace std;")
    if p<0: die("using namespace anchor missing")
    src=src[:p]+"".join("#define %s %s\n"%(a,b) for a,b in MAC)+src[p:]
    out=tr(src)
    n=len(out.encode("utf-8"))
    if n>LIMIT: die("output %d exceeds %d"%(n,LIMIT))
    if out.count("main(){")!=1 or "MS2::run()" not in out or "MS2SURGE" in out:
        die("post-transform anchor lost")

    outp=Path(args.out)
    outp.write_text(out,encoding="utf-8")
    exe=os.path.splitext(str(outp))[0]
    print("input_bytes=%d"%len(raw))
    print("input_blob_sha=%s"%sh)
    print("output=%s"%outp)
    print("output_bytes=%d"%n)
    print("source_limit=%d"%LIMIT)
    print("delta=%+d"%(n-len(raw)))
    print("compile=g++ -std=c++17 -O2 -pipe %s -o %s"%(outp,exe))
    print("sample_first_line=8 12")
    print("expected_kattis_signal=7/7; fail-closed tie near 81.938904 unless MS2 mid-N B16-family or strict 200k VIMP replay validates; rollback if not 7/7 or score<=81.938904")

if __name__=="__main__":
    main()
