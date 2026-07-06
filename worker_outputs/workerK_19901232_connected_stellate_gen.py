#!/usr/bin/env python3
import re,sys,hashlib
from pathlib import Path
LIMIT=131052
MAIN='int main(){JC();GN();if(!W2G::run())W2C::run();W5::post_patch_pass();VIMP::run();MIDEC::run();WK::run();B16::R(39000,60000,220,-7,192,.96,18.05);if(0&&"B16P515A"){}B16::R(39000,60000,76,-10,192,.96,18.35);for(int i=2;i--&&N>47500&N<6e4&&es()<18.6;)WK::run();if(N<50625&&es()<18.9)WK::run();JD();}'
NEWM='int main(){JC();GN();if(!W2G::run())W2C::run();W5::post_patch_pass();VIMP::run();MIDEC::run();WK::run();B16::R(39000,60000,220,-7,192,.96,18.05);if(0&&"B16P515A"){}B16::R(39000,60000,76,-10,192,.96,18.35);for(int i=2;i--&&N>47500&N<6e4&&es()<18.6;)WK::run();KSV::run();if(N<50625&&es()<18.9)WK::run();JD();}'
REQ=['namespace B16{','namespace W2G{','namespace W2C{','namespace WK{','namespace W5{','namespace VIMP{','namespace MIDEC{']
FORBID=['namespace KSV{','namespace BX8{','namespace UVH{','WD232','DX4','WITNESS_TET','B16::R(39000,60000,88,-9','B16::R(39000,60000,96,-11','B16::R(39000,60000,120,-11','B16::R(39000,60000,88,-12','B16::R(39000,60000,99,-11','B16::R(39000,60000,76,-13']
INS=r'''namespace KSV{struct E{double l;unsigned long long k;};static Face F(int a,int b,int c){Face f;f.v[0]=a;f.v[1]=b;f.v[2]=c;return f;}static void A(const Face&f){int id=(int)faces.size();faces.pb(f);BR.pb(1);++BE;Y[f.v[0]].pb(id);Y[f.v[1]].pb(id);Y[f.v[2]].pb(id);}static bool S(int tgt,double lim){AE p;p.AI=.115*CL;p.BB=.070*CL;p.BQ=cos(64.*acos(-1.)/180.);p.W=max(1e-11,1e-9*CL);p.AT=.15;p.AJ=max(1e-30,1e-24*CL*CL);p.AA=false;vector<unsigned long long>ks;ks.reserve((size_t)BE*3);for(int i=0;i<(int)faces.size();++i){if((i&8191)==0&&es()>lim-.55)return false;if(!BR[i])continue;Face f=faces[i];if(!BU[f.v[0]]||!BU[f.v[1]]||!BU[f.v[2]])continue;ks.pb(ED(f.v[0],f.v[1]));ks.pb(ED(f.v[1],f.v[2]));ks.pb(ED(f.v[2],f.v[0]));}sort(ks.begin(),ks.end());ks.erase(unique(ks.begin(),ks.end()),ks.end());vector<E>ed;ed.reserve(ks.size());for(auto k:ks){int a=(int)(k>>32),b=(int)(k&0xffffffffu);if(a>=0&&b>=0&&a<N&&b<N&&BU[a]&&BU[b])ed.pb({norm2(P[a]-P[b]),k});}sort(ed.begin(),ed.end(),[](const E&a,const E&b){return a.l<b.l;});int r=0;for(auto&e:ed){if(cove()<=tgt||es()>lim)break;int a=(int)(e.k>>32),b=(int)(e.k&0xffffffffu);if(a>=0&&b>=0&&a<N&&b<N&&BU[a]&&BU[b]&&GD(a,b,p))++r;}return r>0;}struct G{Vec3 mn,mx;double r,r2,c;int nx,ny,nz;vector<vector<int>>B;int cl(int x,int n){return x<0?0:(x>=n?n-1:x);}int ix(double x){return cl((int)((x-mn.x)/c),nx);}int iy(double y){return cl((int)((y-mn.y)/c),ny);}int iz(double z){return cl((int)((z-mn.z)/c),nz);}int key(int x,int y,int z){return(z*ny+y)*nx+x;}void init(double R){r=R;r2=R*R;mn=mx=originalP[0];for(auto&p:originalP){mn.x=min(mn.x,p.x);mn.y=min(mn.y,p.y);mn.z=min(mn.z,p.z);mx.x=max(mx.x,p.x);mx.y=max(mx.y,p.y);mx.z=max(mx.z,p.z);}c=max(R,1e-9);nx=max(1,(int)((mx.x-mn.x)/c)+1);ny=max(1,(int)((mx.y-mn.y)/c)+1);nz=max(1,(int)((mx.z-mn.z)/c)+1);B.assign((size_t)nx*ny*nz,{});for(int i=0;i<N;i++)B[key(ix(originalP[i].x),iy(originalP[i].y),iz(originalP[i].z))].pb(i);}void mark(const Vec3&p,vector<unsigned char>&C,int&cc){int X=ix(p.x),Y0=iy(p.y),Z0=iz(p.z);for(int z=Z0-1;z<=Z0+1;z++)if(z>=0&&z<nz)for(int y=Y0-1;y<=Y0+1;y++)if(y>=0&&y<ny)for(int x=X-1;x<=X+1;x++)if(x>=0&&x<nx)for(int q:B[key(x,y,z)])if(!C[q]&&norm2(originalP[q]-p)<=r2){C[q]=1;++cc;}}};static int nearv(const Vec3&p,const vector<int>&V){int b=-1;double d=1e100;for(int v:V){double e=norm2(P[v]-p);if(e<d){d=e;b=v;}}return b;}static bool st(int src,int id,const Vec3&p,const vector<int>&base){int v=nearv(p,base),fid=-1;if(v<0)return false;for(int f:Y[v])if(BR[f]&&AC(f,v)){fid=f;break;}if(fid<0)return false;Face o=faces[fid];P[id]=p;BU[id]=1;BF[id]=0;Y[id].clear();BR[fid]=0;--BE;A(F(o.v[0],o.v[1],id));A(F(o.v[1],o.v[2],id));A(F(o.v[2],o.v[0],id));return true;}static bool C(int cap){G g;g.init(.0492*CL);vector<unsigned char>Cv(N,0);int cc=0;vector<int>base,fr;for(int i=0;i<N;i++)if(BU[i]){base.pb(i);g.mark(P[i],Cv,cc);}else fr.pb(i);int fp=0;for(int i=0;i<N&&cc<N;i++)if(!Cv[i]){if(cove()+1>cap||fp>=(int)fr.size()||es()>19.55)return false;if(!st(i,fr[fp++],originalP[i],base))return false;g.mark(originalP[i],Cv,cc);}return cc==N;}static bool run(){if(!(N>=50625&&N<60000)||es()>18.72)return false;if(AS&&(BG>.014||AL<.70||AH>.16))return false;AP R=AD();int sn=cove();if(sn<1000)return false;int tgt=max(96,sn*46/100);if(!S(tgt,19.18)||cove()>sn*86/100){rs(R);return false;}int cap=max(128,min(sn*78/100,N-1));if(!C(cap)){rs(R);return false;}int an=cove();if(!(an>0&&an<sn*80/100&&W5::strong_validator()&&es()<19.85)){rs(R);return false;}double q=vps(256);if(q<(an<sn*62/100?.928:.946)){rs(R);return false;}return true;}}'''
MAC=[('O0','double'),('O1','static'),('O2','const'),('O3','inline'),('O4','vector'),('O5','unsigned'),('O6','namespace'),('O7','struct'),('O8','template'),('O9','typename'),('P0','return false'),('P1','return true'),('P2','continue'),('P3','return'),('P4','int'),('P5','bool'),('P6','void'),('P7','class')]
ID=re.compile(r'[A-Za-z_][A-Za-z0-9_]*')
def die(m): raise SystemExit('workerK stellate abort: '+m)
def toks(src):
    r=[];i=0;n=len(src);bol=True;pp=False
    while i<n:
        c=src[i]
        if bol:
            j=i
            while j<n and src[j] in ' \t':j+=1
            pp=j<n and src[j]=='#';bol=False
        if c=='\n':bol=True;pp=False;i+=1;continue
        if pp:i+=1;continue
        if c in '"\'':
            q=c;i+=1
            while i<n:
                if src[i]=='\\':i+=2;continue
                if src[i]==q:i+=1;break
                i+=1
            continue
        if c=='/' and i+1<n and src[i+1]=='/':
            j=src.find('\n',i);i=n if j<0 else j;continue
        if c=='/' and i+1<n and src[i+1]=='*':
            j=src.find('*/',i+2);i=n if j<0 else j+2;continue
        m=ID.match(src,i)
        if m:r.append(m.group(0));i=m.end();continue
        i+=1
    return r
def tr(src):
    mp={b:a for a,b in MAC if ' ' not in b};out=[];i=0;n=len(src);bol=True;pp=False
    while i<n:
        c=src[i]
        if bol:
            j=i
            while j<n and src[j] in ' \t':j+=1
            pp=j<n and src[j]=='#';bol=False
        if c=='\n':out.append(c);i+=1;bol=True;pp=False;continue
        if pp:out.append(c);i+=1;continue
        if c in '"\'':
            q=c;st=i;i+=1
            while i<n:
                if src[i]=='\\':i+=2;continue
                if src[i]==q:i+=1;break
                i+=1
            out.append(src[st:i]);continue
        if c=='/' and i+1<n and src[i+1]=='/':
            j=src.find('\n',i);j=n if j<0 else j;out.append(src[i:j]);i=j;continue
        if c=='/' and i+1<n and src[i+1]=='*':
            j=src.find('*/',i+2);j=n if j<0 else j+2;out.append(src[i:j]);i=j;continue
        m=ID.match(src,i)
        if m:
            w=m.group(0)
            if w=='return':
                k=m.end()
                while k<n and src[k] in ' \t\r\n':k+=1
                m2=ID.match(src,k)
                if m2 and m2.group(0)=='false':out.append('P0');i=m2.end();continue
                if m2 and m2.group(0)=='true':out.append('P1');i=m2.end();continue
            out.append(mp.get(w,w));i=m.end();continue
        out.append(c);i+=1
    return ''.join(out)
def main():
    if len(sys.argv)>1: inp=Path(sys.argv[1])
    else:
        cand=[Path('fetched_sources/19901232.cpp'),Path('submission_543_81.93_7.cpp')]
        inp=next((p for p in cand if p.exists()),None)
        if inp is None: die('missing fetched_sources/19901232.cpp or submission_543_81.93_7.cpp')
    outp=Path(sys.argv[2]) if len(sys.argv)>2 else Path('workerK_19901232_connected_stellate.cpp')
    src=inp.read_text(encoding='utf-8')
    if src.count(MAIN)!=1: die('exact 19901232 main anchor not found once')
    for x in REQ:
        if src.count(x)!=1: die('required namespace missing/duplicated: '+x)
    for x in FORBID:
        if x in src: die('forbidden stale/drop anchor present: '+x)
    if src.count('static void CA(string&out,const char*line,int len){')!=1: die('CA anchor mismatch')
    src=src.replace('static void CA(string&out,const char*line,int len){',INS+'static void CA(string&out,const char*line,int len){',1).replace(MAIN,NEWM,1)
    col=[a for a,_ in MAC if a in set(toks(src))]
    if col: die('macro token collision: '+','.join(col))
    pos=src.find('using namespace std;')
    if pos<0: die('using anchor missing')
    src=src[:pos]+''.join(f'#define {a} {b}\n' for a,b in MAC)+src[pos:]
    out=tr(src)
    n=len(out.encode())
    if n>LIMIT: die(f'output {n} bytes exceeds {LIMIT}')
    if out.count('main(){')!=1 or 'KSV::run()' not in out: die('post-transform anchor lost')
    outp.write_text(out,encoding='utf-8')
    print(f'wrote {outp} bytes={n} sha256={hashlib.sha256(out.encode()).hexdigest()}')
if __name__=='__main__': main()