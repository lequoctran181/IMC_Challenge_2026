#!/usr/bin/env python3
import re,sys,hashlib
from pathlib import Path
from collections import Counter
LIMIT=131072
LANE=r'''namespace HS{struct E{double l;unsigned long long k;};static Vec3 U(Vec3 a){double n=norm3(a);return n>1e-30?a*(1./n):Vec3{0,0,1};}static Face F(int a,int b,int c){Face f;f.v[0]=a;f.v[1]=b;f.v[2]=c;return f;}static int AF(Face f){int id=(int)faces.size();faces.pb(f);BR.pb(1);++BE;M=(int)faces.size();Y[f.v[0]].pb(id);Y[f.v[1]].pb(id);Y[f.v[2]].pb(id);return id;}static bool LF(int i){if(i<0||i>=(int)faces.size()||!BR[i])return false;Face f=faces[i];return f.v[0]>=0&&f.v[1]>=0&&f.v[2]>=0&&f.v[0]<N&&f.v[1]<N&&f.v[2]<N&&BU[f.v[0]]&&BU[f.v[1]]&&BU[f.v[2]]&&f.v[0]!=f.v[1]&&f.v[0]!=f.v[2]&&f.v[1]!=f.v[2];}static Vec3 CEN(const Face&f){return(P[f.v[0]]+P[f.v[1]]+P[f.v[2]])*(1./3.);}static double A2(const Face&f){return norm2(cross3(P[f.v[1]]-P[f.v[0]],P[f.v[2]]-P[f.v[0]]));}static bool EC(int tgt,double lim,double ai,double bb,double ang){AE p;p.AI=ai*CL;p.BB=bb*CL;p.BQ=cos(ang*acos(-1.)/180.);p.W=max(1e-12,1e-10*CL);p.AT=.015;p.AJ=max(1e-30,1e-25*CL*CL);p.AA=false;int st=cove(),pass=0;while(cove()>tgt&&pass++<2&&es()<lim){vector<unsigned long long>ks;ks.reserve((size_t)BE*3);for(int i=0;i<(int)faces.size();++i){if((i&8191)==0&&es()>lim-.35)return cove()<st;if(!LF(i))continue;Face f=faces[i];ks.pb(ED(f.v[0],f.v[1]));ks.pb(ED(f.v[1],f.v[2]));ks.pb(ED(f.v[2],f.v[0]));}sort(ks.begin(),ks.end());ks.erase(unique(ks.begin(),ks.end()),ks.end());vector<E>ed;ed.reserve(ks.size());for(auto k:ks){int a=(int)(k>>32),b=(int)(k&0xffffffffu);if(a>=0&&b>=0&&a<N&&b<N&&BU[a]&&BU[b])ed.pb({norm2(P[a]-P[b]),k});}sort(ed.begin(),ed.end(),[](const E&a,const E&b){return a.l<b.l;});int ok=0;for(auto&e:ed){if(cove()<=tgt||es()>lim)break;int a=(int)(e.k>>32),b=(int)(e.k&0xffffffffu);if(a>=0&&b>=0&&a<N&&b<N&&BU[a]&&BU[b]&&GD(a,b,p))++ok;}if(!ok)break;}return cove()<st;}struct G{Vec3 mn,mx;double r,r2,c;int nx,ny,nz;vector<vector<int>>B;int cl(int x,int n){return x<0?0:(x>=n?n-1:x);}int ix(double x){return cl((int)((x-mn.x)/c),nx);}int iy(double y){return cl((int)((y-mn.y)/c),ny);}int iz(double z){return cl((int)((z-mn.z)/c),nz);}int key(int x,int y,int z){return(z*ny+y)*nx+x;}void init(double R){r=R;r2=R*R;mn=mx=originalP[0];for(auto&p:originalP){mn.x=min(mn.x,p.x);mn.y=min(mn.y,p.y);mn.z=min(mn.z,p.z);mx.x=max(mx.x,p.x);mx.y=max(mx.y,p.y);mx.z=max(mx.z,p.z);}c=max(R,1e-9);nx=max(1,(int)((mx.x-mn.x)/c)+1);ny=max(1,(int)((mx.y-mn.y)/c)+1);nz=max(1,(int)((mx.z-mn.z)/c)+1);if((long long)nx*ny*nz>900000){nx=ny=nz=1;c=max({mx.x-mn.x,mx.y-mn.y,mx.z-mn.z})+1.;}B.assign((size_t)nx*ny*nz,{});for(int i=0;i<N;i++)B[key(ix(originalP[i].x),iy(originalP[i].y),iz(originalP[i].z))].pb(i);}void mark(const Vec3&p,vector<unsigned char>&C,int&cc){int X=ix(p.x),Y0=iy(p.y),Z0=iz(p.z);for(int z=Z0-1;z<=Z0+1;z++)if(z>=0&&z<nz)for(int y=Y0-1;y<=Y0+1;y++)if(y>=0&&y<ny)for(int x=X-1;x<=X+1;x++)if(x>=0&&x<nx)for(int q:B[key(x,y,z)])if(!C[q]&&norm2(originalP[q]-p)<=r2){C[q]=1;++cc;}}};static vector<Vec3> VN(){vector<Vec3>n(N);Vec3 cen{};for(auto&p:originalP)cen=cen+p;cen=cen*(1./max(1,N));double s=0;for(auto&f:AR){if(f.v[0]<0||f.v[1]<0||f.v[2]<0||f.v[0]>=N||f.v[1]>=N||f.v[2]>=N)continue;Vec3 cr=cross3(originalP[f.v[1]]-originalP[f.v[0]],originalP[f.v[2]]-originalP[f.v[0]]);Vec3 ce=(originalP[f.v[0]]+originalP[f.v[1]]+originalP[f.v[2]])*(1./3.);s+=dot3(cr,ce-cen);n[f.v[0]]=n[f.v[0]]+cr;n[f.v[1]]=n[f.v[1]]+cr;n[f.v[2]]=n[f.v[2]]+cr;}double sg=s>=0?1.:-1.;for(int i=0;i<N;i++){Vec3 o=U(n[i]*sg);Vec3 r=U(originalP[i]-cen);if(norm2(n[i])<1e-24||dot3(o,r)<-.25)o=r;n[i]=o;}return n;}static bool COV(vector<int>&id,vector<Vec3>&pos,int cap,double H,double off,double lim){G g;g.init(H);vector<unsigned char>C(N,0);int cc=0;for(int i=0;i<N;i++)if(BU[i]){if((i&4095)==0&&es()>lim)return false;g.mark(P[i],C,cc);}vector<int>fr;fr.reserve(N-cove());for(int i=0;i<N;i++)if(!BU[i])fr.pb(i);vector<Vec3>nr=VN();int fp=0,base=cove(),step=max(1,N/190000);for(int o=0;o<step&&cc<N;o++)for(int i=o;i<N&&cc<N;i+=step)if(!C[i]){if(base+(int)id.size()+1>cap||fp>=(int)fr.size()||es()>lim)return false;Vec3 p=originalP[i]-nr[i]*off;id.pb(fr[fp++]);pos.pb(p);g.mark(p,C,cc);}return cc==N;}static int HOLE(const vector<Vec3>&pos){if(pos.empty())return-1;Vec3 q=pos[0];int h=-1;double bs=1e100,cl2=CL*CL,cl4=cl2*cl2+1e-300;for(int i=0;i<(int)faces.size();++i){if(!LF(i))continue;double a=A2(faces[i]);if(!(a>1e-26*cl4))continue;Vec3 c=CEN(faces[i]);double sc=a/cl4+.015*norm2(c-q)/(cl2+1e-300);if(sc<bs){bs=sc;h=i;}}return h;}static void ORD(vector<int>&ord,const vector<Vec3>&pos,Vec3 cur){int W=pos.size();ord.clear();ord.reserve(W);vector<unsigned char>u(W,0);for(int k=0;k<W;k++){int b=-1;double bd=1e100;for(int i=0;i<W;i++)if(!u[i]){double d=norm2(pos[i]-cur);if(d<bd){bd=d;b=i;}}if(b<0)break;u[b]=1;ord.pb(b);cur=pos[b];}}static bool ATT(const vector<int>&id,const vector<Vec3>&pos,double lim){if(id.empty())return true;int hf=HOLE(pos);if(hf<0)return false;Face cell=faces[hf];BR[hf]=0;--BE;vector<int>ord;ORD(ord,pos,CEN(cell));int cellFace=-1;for(int t=0;t<(int)ord.size();++t){if(es()>lim)return false;if(cellFace>=0){if(!LF(cellFace))return false;cell=faces[cellFace];BR[cellFace]=0;--BE;}int w=id[ord[t]];P[w]=pos[ord[t]];BU[w]=1;BF[w]=0;Y[w].clear();int f0=AF(F(cell.v[0],cell.v[1],w));int f1=AF(F(cell.v[1],cell.v[2],w));int f2=AF(F(cell.v[2],cell.v[0],w));if(t+1<(int)ord.size()){Vec3 nx=pos[ord[t+1]];int fs[3]={f0,f1,f2};double bd=1e100;cellFace=f0;for(int j=0;j<3;j++){double d=norm2(CEN(faces[fs[j]])-nx);if(d<bd){bd=d;cellFace=fs[j];}}}}return true;}static bool TRY(const AP&S,int tr,int capr,double off,double ai,double bb,double ang,double qmin,double lim){rs(S);int sn=cove();if(!EC(max(80,sn*tr/100),lim-.75,ai,bb,ang))return false;int mid=cove();if(mid>=sn*92/100||mid<24)return false;vector<int>id;vector<Vec3>pos;int cap=max(96,min(sn*capr/100,N-1));if(!COV(id,pos,cap,.0492*CL,off*CL,lim-.42))return false;if((int)id.size()<8||mid+(int)id.size()>=sn*95/100)return false;if(!ATT(id,pos,lim-.12))return false;int an=cove();if(!(an>0&&an<sn*94/100&&es()<lim))return false;if(!W5::strong_validator())return false;if(es()>20.15)return false;double q=vps(512);if(q<qmin)return false;return true;}static bool run(){if(!(N>=18000&&N<90000)||es()>18.35)return false;if(AS&&(BG>.028||AL<.50||AH>.32))return false;int sn=cove();if(sn<1100||sn>N*9/10)return false;AP S=AD(),B;bool ok=false;int bc=sn;struct Pm{int tr,cr;double of,ai,bb,ag,q;};Pm ps[4]={{34,72,.021,.30,.19,78,.942},{26,66,.018,.40,.24,86,.930},{19,58,.015,.54,.33,94,.918},{12,50,.012,.78,.45,104,.904}};for(auto&p:ps){if(es()>18.95)break;if(TRY(S,p.tr,p.cr,p.of,p.ai,p.bb,p.ag,p.q,19.75)){int c=cove();if(c<bc){bc=c;B=AD();ok=true;}}}if(ok){rs(B);return true;}rs(S);return false;}}'''
ID=re.compile(r'[A-Za-z_][A-Za-z0-9_]*')
def die(s): raise SystemExit('hs_gen: '+s)
def bsha(b): return hashlib.sha1(b'blob '+str(len(b)).encode()+b'\0'+b).hexdigest()
def mb(s,o):
    d=0;i=o;q=None;esc=False
    while i<len(s):
        c=s[i]
        if q:
            if esc: esc=False
            elif c=='\\': esc=True
            elif c==q: q=None
        else:
            if c in '"\'': q=c
            elif c=='{': d+=1
            elif c=='}':
                d-=1
                if d==0: return i+1
        i+=1
    return -1
def fm(s):
    m=re.search(r'(^|[^A-Za-z0-9_])(?:signed\s+)?int\s+main\s*\([^)]*\)\s*\{',s)
    if not m: die('main not found')
    a=m.start(0)+(1 if m.group(1) else 0);o=s.find('{',a);e=mb(s,o)
    if e<0: die('unmatched main')
    return a,e
def nl(s,i):
    n=len(s);j=i
    if s[j]=='0' and j+1<n and s[j+1] in 'xX':
        j+=2
        while j<n and (s[j].isdigit() or s[j].lower() in 'abcdef'): j+=1
        if j<n and s[j] in 'pP':
            j+=1
            if j<n and s[j] in '+-': j+=1
            while j<n and s[j].isdigit(): j+=1
    else:
        seen=False
        while j<n:
            c=s[j]
            if c.isdigit() or c=='.': seen=True;j+=1;continue
            if c in 'eEpP' and seen:
                j+=1
                if j<n and s[j] in '+-': j+=1
                continue
            break
    while j<n and s[j].isalpha(): j+=1
    return j
def scan(s,rep=None,cnt=None,used=None,defs=None):
    out=[];i=0;n=len(s);bol=True;pp=False
    while i<n:
        c=s[i]
        if bol:
            j=i
            while j<n and s[j] in ' \t': j+=1
            pp=j<n and s[j]=='#'; bol=False
        if c=='\n':
            out.append(c);i+=1;bol=True;pp=False;continue
        if (c.isdigit() or (c=='.' and i+1<n and s[i+1].isdigit())) and not pp:
            j=nl(s,i);out.append(s[i:j]);i=j;continue
        if c in '"\'':
            q=c;st=i;i+=1
            while i<n:
                if s[i]=='\\': i+=2;continue
                if s[i]==q: i+=1;break
                i+=1
            out.append(s[st:i]);continue
        if not pp and c=='/' and i+1<n and s[i+1]=='/':
            j=s.find('\n',i);i=n if j<0 else j;continue
        if not pp and c=='/' and i+1<n and s[i+1]=='*':
            j=s.find('*/',i+2);i=n if j<0 else j+2;continue
        m=ID.match(s,i)
        if m:
            w=m.group(0)
            if pp and defs is not None and len(out)>=1 and ''.join(out[-1:]).endswith('#define '): defs.add(w)
            if used is not None: used.add(w)
            if cnt is not None and not pp: cnt[w]+=1
            out.append(rep.get(w,w) if rep and not pp else w)
            i=m.end();continue
        out.append(c);i+=1
    return ''.join(out)
def names(used):
    a='ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz'
    k=0
    while 1:
        for c in a:
            z=c+str(k)
            if z not in used: yield z
        k+=1
def adddefs(s,defs):
    if not defs: return s
    block=''.join('#define %s %s\n'%(a,b) for a,b in defs)
    pos=0
    while True:
        m=re.match(r'[ \t]*#[^\n]*(?:\n|$)',s[pos:])
        if not m: break
        pos+=m.end()
    return s[:pos]+block+s[pos:]
def compact(s):
    c=Counter();u=set();d=set()
    scan(s,cnt=c,used=u,defs=d)
    u|=d
    ban={'main'}
    cand=[]
    for w,n in c.items():
        if w in ban or w.startswith('__') or len(w)<3: continue
        cand.append(((len(w)-2)*n,w,n))
    cand.sort(reverse=True)
    g=names(u);rep={};defs=[]
    for _,w,n in cand:
        nm=next(g);gain=(len(w)-len(nm))*n-(9+len(nm)+len(w))
        if gain<=6: continue
        rep[w]=nm;defs.append((nm,w))
        if len(defs)>=560: break
    return adddefs(scan(s,rep=rep),defs)
def patch(s):
    req=['AD()','rs(','GD(','ED(','W5::strong_validator','vps','cove','originalP','faces','BR','BU','Y','AR','CL']
    miss=[x for x in req if x not in s]
    if miss: die('missing anchors '+','.join(miss))
    if 'namespace HS{' in s: die('already patched')
    a,e=fm(s);m=s[a:e];p=m.rfind('JD();')
    if p<0:
        r=m.rfind('return')
        p=r if r>=0 else len(m)-1
    return s[:a]+LANE+m[:p]+'HS::run();'+m[p:]+s[e:]
def pick():
    if len(sys.argv)>1 and sys.argv[1]!='-':
        p=Path(sys.argv[1])
        if not p.exists(): die('missing '+str(p))
        return p,p.read_text()
    for x in ('fetched_sources/kattis_19903544_81.938904.cpp','submission_576_81.93_7.cpp','submission_563_81.93_7.cpp','fetched_sources/19901232.cpp'):
        p=Path(x)
        if p.exists(): return p,p.read_text()
    die('pass current-best source path')
def main():
    p,s=pick()
    outp=Path(sys.argv[2]) if len(sys.argv)>2 else Path('hs.cpp')
    raw=patch(s)
    out=compact(raw)
    if len(out.encode())>LIMIT:
        die('patched source too large raw=%d compact=%d limit=%d'%(len(raw.encode()),len(out.encode()),LIMIT))
    outp.write_text(out)
    print('base=%s'%p)
    print('bytes=%d limit=%d base_blob_sha=%s'%(len(out.encode()),LIMIT,bsha(s.encode())))
    print('compile=g++ -std=c++17 -O2 -pipe -static -s %s -o hs'%outp)
    print('sample=./hs < sample.in | head -1   # expect 8 12')
if __name__=='__main__': main()
