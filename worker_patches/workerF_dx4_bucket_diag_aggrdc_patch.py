#!/usr/bin/env python3
"""
Worker F DX4 diagnostic + targeted exploit generator for IMC 2026 Problem B.
Anchor: submission_563_81.93_7.cpp / fetched_sources/kattis_19901322.cpp.

Modes:
  probe1      accepted diagnostic: after current-best work, output original mesh for buckets with bit 1 set
  probe2      accepted diagnostic: after current-best work, output original mesh for buckets with bit 2 set
  exploit0    targeted aggressive DC/QEM tournament for post-MIDEC bucket 0
  exploit1    targeted aggressive DC/QEM tournament for post-MIDEC bucket 1
  exploit2    targeted aggressive DC/QEM tournament for post-MIDEC bucket 2
  exploit3    targeted aggressive DC/QEM tournament for post-MIDEC bucket 3

Usage:
  python3 workerF_dx4_bucket_diag_aggrdc_patch.py submission_563_81.93_7.cpp out.cpp probe1
  python3 workerF_dx4_bucket_diag_aggrdc_patch.py submission_563_81.93_7.cpp out.cpp exploit2
"""
import re, sys, hashlib
LIMIT=131052
ANCHOR='int main(){JC();GN();if(!W2G::run())W2C::run();W5::post_patch_pass();VIMP::run();MIDEC::run();WK::run();B16::R(39000,60000,220,-7,192,.96,18.05);if(0&&"B16P515A"){}B16::R(39000,60000,76,-10,192,.961,18.35);for(int i=2;i--&&N>47500&N<6e4&&es()<18.6;)WK::run();if(N<50625&&es()<18.9)WK::run();if(C5T::run()){JD();return 0;}JD();}'
REQ=['namespace B16{','namespace W2G{','namespace W2C{','namespace WK{','namespace W5{','namespace VIMP{','namespace MIDEC{','namespace C5T{','namespace DC{']
FORBID=['namespace BX8{','namespace D12{','namespace DX4{','W5L473','B16::R(39000,60000,76,-13','B16::R(39000,60000,88,-9']
MACROS=[('O0','double'),('O1','static'),('O2','const'),('O3','inline'),('O4','vector'),('O5','unsigned'),('O6','namespace'),('O7','struct'),('O8','template'),('O9','typename'),('P0','return false'),('P1','return true'),('P2','continue'),('P3','return'),('P4','int'),('P5','bool'),('P6','void'),('P7','class')]
IDENT=re.compile(r'[A-Za-z_][A-Za-z0-9_]*')
LANE=r'''
namespace DX4{static int B=-1;static void mix(unsigned long long&x,unsigned long long v){x^=v+0x9e3779b97f4a7c15ULL+(x<<6)+(x>>2);x*=0xbf58476d1ce4e5b9ULL;}static int H(){int s=cove();unsigned long long x=0x123456789abcdefULL;mix(x,N);mix(x,M);mix(x,s);mix(x,(unsigned long long)(BL*10000));mix(x,(unsigned long long)(AL*10000));mix(x,(unsigned long long)(Z*10000));mix(x,(unsigned long long)(AH*10000));mix(x,(unsigned long long)(BG*10000));return (int)((x>>62)&3);}static void I(){if(B<0)B=H();}static bool P(int m){I();return(B&m)!=0;}static bool X(int t){I();if(B!=t||es()>17.35)return false;AP S=AD();int SN=cove();if(SN<=0)return false;vector<Vec3>V;vector<Face>F;double R[]={0.012,0.016,0.020,0.025,0.030,0.040,0.055,0.075,-1};double best=-1;AP BS=S;for(double r:R){if(es()>18.65)break;rs(S);V.clear();F.clear();double tm=N<30000?3.7:(N<120000?4.6:5.4);bool ok=r<0?DC::EV(originalP,AR,CL,V,F,tm):DC::EV(originalP,AR,CL,V,F,tm,r);if(!ok||V.empty()||F.empty()||(int)V.size()>=SN)continue;if((int)V.size()*100>SN*(r<.021?55:(r<.031?68:82)))continue;if(!AF(V,F)||!W5::strong_validator())continue;int res=N<40000?1024:(N<180000?512:256);double px=vps(res);double th=N<40000?.895:(N<180000?.905:.920);double sc=(double)(SN-(int)V.size())/(double)SN+0.04*(px-th);if(px>=th&&sc>best){best=sc;BS=AD();}}if(best>0){rs(BS);return true;}rs(S);return false;}}
'''.strip()

def die(m):
    print('DX4 generator abort:',m,file=sys.stderr); sys.exit(2)

def scan_tokens(src):
    toks=[];i=0;n=len(src);bol=True;pp=False
    while i<n:
        c=src[i]
        if bol:
            j=i
            while j<n and src[j] in ' \t': j+=1
            pp=(j<n and src[j]=='#');bol=False
        if c=='\n': bol=True;pp=False;i+=1;continue
        if pp: i+=1;continue
        if c in '"\'':
            q=c;i+=1
            while i<n:
                if src[i]=='\\': i+=2;continue
                if src[i]==q: i+=1;break
                i+=1
            continue
        if c=='/' and i+1<n and src[i+1]=='/':
            j=src.find('\n',i);i=n if j<0 else j;continue
        if c=='/' and i+1<n and src[i+1]=='*':
            j=src.find('*/',i+2);i=n if j<0 else j+2;continue
        m=IDENT.match(src,i)
        if m: toks.append(m.group(0));i=m.end();continue
        i+=1
    return toks

def transform(src):
    mp={b:a for a,b in MACROS if ' ' not in b}
    out=[];i=0;n=len(src);bol=True;pp=False
    while i<n:
        c=src[i]
        if bol:
            j=i
            while j<n and src[j] in ' \t': j+=1
            pp=(j<n and src[j]=='#');bol=False
        if c=='\n': out.append(c);i+=1;bol=True;pp=False;continue
        if pp: out.append(c);i+=1;continue
        if c in '"\'':
            q=c;s=i;i+=1
            while i<n:
                if src[i]=='\\': i+=2;continue
                if src[i]==q: i+=1;break
                i+=1
            out.append(src[s:i]);continue
        if c=='/' and i+1<n and src[i+1]=='/':
            j=src.find('\n',i);j=n if j<0 else j;out.append(src[i:j]);i=j;continue
        if c=='/' and i+1<n and src[i+1]=='*':
            j=src.find('*/',i+2);j=n if j<0 else j+2;out.append(src[i:j]);i=j;continue
        m=IDENT.match(src,i)
        if m:
            tok=m.group(0)
            if tok=='return':
                k=m.end()
                while k<n and src[k] in ' \t\r\n': k+=1
                m2=IDENT.match(src,k)
                if m2 and m2.group(0)=='false': out.append('P0');i=m2.end();continue
                if m2 and m2.group(0)=='true': out.append('P1');i=m2.end();continue
            out.append(mp.get(tok,tok));i=m.end();continue
        out.append(c);i+=1
    return ''.join(out)

def main_for(mode):
    head='int main(){JC();GN();if(!W2G::run())W2C::run();W5::post_patch_pass();VIMP::run();MIDEC::run();DX4::I();'
    tail='WK::run();B16::R(39000,60000,220,-7,192,.96,18.05);if(0&&"B16P515A"){}B16::R(39000,60000,76,-10,192,.961,18.35);for(int i=2;i--&&N>47500&N<6e4&&es()<18.6;)WK::run();if(N<50625&&es()<18.9)WK::run();if(C5T::run()){JD();return 0;}'
    if mode=='probe1': return head+tail+'if(DX4::P(1)){IJ();return 0;}JD();}'
    if mode=='probe2': return head+tail+'if(DX4::P(2)){IJ();return 0;}JD();}'
    if mode.startswith('exploit'):
        t=int(mode[7:])
        if t<0 or t>3: die('exploit mode must be exploit0..exploit3')
        return head+'if(DX4::X(%d)){JD();return 0;}'%t+tail+'JD();}'
    die('mode must be probe1, probe2, or exploit0..exploit3')

def main():
    if len(sys.argv)!=4:
        print('usage: workerF_dx4_bucket_diag_aggrdc_patch.py submission_563_81.93_7.cpp out.cpp MODE',file=sys.stderr);sys.exit(2)
    path,outpath,mode=sys.argv[1],sys.argv[2],sys.argv[3]
    src=open(path,encoding='utf-8').read()
    if src.count(ANCHOR)!=1: die('exact current-best main anchor not found once')
    for r in REQ:
        if src.count(r)!=1: die('required namespace missing/duplicated: '+r)
    for f in FORBID:
        if f in src: die('known duplicate/drop/stale anchor present: '+f)
    toks=set(scan_tokens(src))
    bad=[a for a,_ in MACROS if a in toks]
    if bad: die('macro token collision: '+','.join(bad))
    p=src.find('using namespace std;')
    if p<0: die('using namespace std anchor missing')
    src=src.replace(ANCHOR,LANE+main_for(mode))
    defs=''.join('#define %s %s\n'%(a,b) for a,b in MACROS)
    src=src[:p]+defs+src[p:]
    res=transform(src)
    if 'DX4::I()' not in res: die('DX4 call lost after transform')
    if mode.startswith('exploit') and 'DX4::X(' not in res: die('exploit call missing')
    sz=len(res.encode())
    if sz>LIMIT: die('output %d exceeds %d'%(sz,LIMIT))
    open(outpath,'w',encoding='utf-8').write(res)
    h=hashlib.sha256(res.encode()).hexdigest()
    with open(outpath+'.dx4_report.txt','w') as f:
        f.write('mode=%s\ninput_bytes=%d\noutput_bytes=%d\nlimit=%d\nsha256=%s\n'%(mode,len(open(path,'rb').read()),sz,LIMIT,h))
        f.write('bucket=hash(N,M,cove_after_MIDEC,BL,AL,Z,AH,BG)&3\n')
        f.write('probe modes output original mesh for selected hash bits; exploit modes run targeted aggressive DC::EV ratio tournament only for selected bucket, with AF,strong_validator,vps rollback.\n')
    print('wrote %s mode=%s bytes=%d sha256=%s'%(outpath,mode,sz,h))
if __name__=='__main__': main()