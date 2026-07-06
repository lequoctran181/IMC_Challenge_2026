#!/usr/bin/env python3
# QVC quotient-remesh generator.
# Usage:
#   python3 gen_qvc.py current_best.cpp submit.cpp
#   g++ -std=c++17 -O2 -pipe -static -s submit.cpp -o submit
# Sample gate expectation: official sample first line remains `8 12`; if QVC fails runtime validator/proxy, it restores the current-best state before JD().

import re,sys,hashlib
from pathlib import Path

LIMIT=131072
ID=re.compile(r'[A-Za-z_][A-Za-z0-9_]*')

LANE=r'''namespace QVC{struct SS{vector<Vec3>P;vector<Face>F;vector<unsigned char>U,R;vector<double>B;vector<vector<int>>Y;int E,M;};static SS sv(){SS s;s.P=P;s.F=faces;s.U=BU;s.R=BR;s.B=BF;s.Y=Y;s.E=BE;s.M=M;return s;}static void ld(const SS&s){P=s.P;faces=s.F;BU=s.U;BR=s.R;BF=s.B;Y=s.Y;BE=s.E;M=s.M;}struct G{Vec3 mn,mx;double h;int nx,ny,nz;vector<vector<int>>b;int cl(int x,int n){return x<0?0:(x>=n?n-1:x);}int ix(double x){return cl((int)((x-mn.x)/h),nx);}int iy(double y){return cl((int)((y-mn.y)/h),ny);}int iz(double z){return cl((int)((z-mn.z)/h),nz);}int ky(int x,int y,int z){return(z*ny+y)*nx+x;}void init(double H){h=max(H,1e-12);mn=mx=originalP[0];for(auto&p:originalP){mn.x=min(mn.x,p.x);mn.y=min(mn.y,p.y);mn.z=min(mn.z,p.z);mx.x=max(mx.x,p.x);mx.y=max(mx.y,p.y);mx.z=max(mx.z,p.z);}nx=max(1,(int)((mx.x-mn.x)/h)+1);ny=max(1,(int)((mx.y-mn.y)/h)+1);nz=max(1,(int)((mx.z-mn.z)/h)+1);if((long long)nx*ny*nz>900000){nx=ny=nz=1;}b.assign((size_t)nx*ny*nz,{});}void add(int i){Vec3 p=originalP[i];b[ky(ix(p.x),iy(p.y),iz(p.z))].push_back(i);}};static unsigned long long ed(int a,int b){if(a>b)swap(a,b);return((unsigned long long)(unsigned int)a<<32)|(unsigned int)b;}struct T{int a,b,c;Face f;};static bool man(int V,const vector<Face>&fs,const vector<int>&sel){if(V<=3||fs.size()<4)return false;double e=max(1e-30,1e-24*CL*CL);vector<unsigned long long>x;x.reserve(fs.size()*3);for(auto&f:fs){if(f.v[0]<0||f.v[1]<0||f.v[2]<0||f.v[0]>=V||f.v[1]>=V||f.v[2]>=V)return false;if(f.v[0]==f.v[1]||f.v[0]==f.v[2]||f.v[1]==f.v[2])return false;Vec3 A=originalP[sel[f.v[0]]],B=originalP[sel[f.v[1]]],C=originalP[sel[f.v[2]]];if(norm2(cross3(B-A,C-A))<=e)return false;x.push_back(ed(f.v[0],f.v[1]));x.push_back(ed(f.v[1],f.v[2]));x.push_back(ed(f.v[2],f.v[0]));}sort(x.begin(),x.end());for(size_t i=0;i<x.size();){size_t j=i+1;while(j<x.size()&&x[j]==x[i])++j;if(j-i!=2)return false;i=j;}return true;}static bool mk(double R,double qmin,double lim){double r=R*CL,r2=r*r;G og;og.init(r);for(int i=0;i<N;i++)og.add(i);vector<unsigned char>cov(N,0);vector<int>sel;sel.reserve(N/4+16);int cc=0;for(int i=0;i<N&&es()<lim-.45;i++)if(!cov[i]){sel.push_back(i);Vec3 p=originalP[i];int X=og.ix(p.x),Y0=og.iy(p.y),Z0=og.iz(p.z);for(int z=Z0-1;z<=Z0+1;z++)if(z>=0&&z<og.nz)for(int y=Y0-1;y<=Y0+1;y++)if(y>=0&&y<og.ny)for(int x=X-1;x<=X+1;x++)if(x>=0&&x<og.nx)for(int q:og.b[og.ky(x,y,z)])if(!cov[q]&&norm2(originalP[q]-p)<=r2){cov[q]=1;++cc;}}if(cc<N||sel.size()<12||sel.size()>N*88/100)return false;G sg;sg.init(r);vector<int>rev(N,-1);for(int i=0;i<(int)sel.size();i++){rev[sel[i]]=i;sg.add(sel[i]);}vector<int>mp(N,-1);for(int i=0;i<N;i++){Vec3 p=originalP[i];int X=sg.ix(p.x),Y0=sg.iy(p.y),Z0=sg.iz(p.z),bi=-1;double bd=1e100;for(int z=Z0-1;z<=Z0+1;z++)if(z>=0&&z<sg.nz)for(int y=Y0-1;y<=Y0+1;y++)if(y>=0&&y<sg.ny)for(int x=X-1;x<=X+1;x++)if(x>=0&&x<sg.nx)for(int s:sg.b[sg.ky(x,y,z)]){double d=norm2(originalP[i]-originalP[s]);if(d<bd){bd=d;bi=rev[s];}}if(bi<0||bd>r2*1.000001)return false;mp[i]=bi;}vector<T>tmp;tmp.reserve(AR.size());double ae=max(1e-30,1e-24*CL*CL);for(auto&o:AR){int a=mp[o.v[0]],b=mp[o.v[1]],c=mp[o.v[2]];if(a<0||b<0||c<0||a==b||a==c||b==c)continue;Vec3 A=originalP[sel[a]],B=originalP[sel[b]],C=originalP[sel[c]];if(norm2(cross3(B-A,C-A))<=ae)continue;int x=a,y=b,z=c;if(x>y)swap(x,y);if(y>z)swap(y,z);if(x>y)swap(x,y);tmp.push_back({x,y,z,{{a,b,c}}});}sort(tmp.begin(),tmp.end(),[](const T&u,const T&v){return u.a!=v.a?u.a<v.a:u.b!=v.b?u.b<v.b:u.c<v.c;});vector<Face>nf;nf.reserve(tmp.size());for(size_t i=0;i<tmp.size();){size_t j=i+1;while(j<tmp.size()&&tmp[j].a==tmp[i].a&&tmp[j].b==tmp[i].b&&tmp[j].c==tmp[i].c)++j;nf.push_back(tmp[i].f);i=j;}int V=sel.size();if(nf.size()<V*2/3||nf.size()>V*8||!man(V,nf,sel))return false;P.assign(N,Vec3{});BU.assign(N,0);BF.assign(N,0);Y.assign(N,vector<int>());for(int i=0;i<V;i++){P[i]=originalP[sel[i]];BU[i]=1;}faces=nf;BR.assign(faces.size(),1);BE=M=(int)faces.size();for(int i=0;i<M;i++){Y[faces[i].v[0]].push_back(i);Y[faces[i].v[1]].push_back(i);Y[faces[i].v[2]].push_back(i);}if(!W5::strong_validator())return false;double q=vps(256);return q>=qmin;}static bool run(){if(N<3500||es()>18.55)return false;SS S=sv(),B;bool ok=false;int sc=cove(),bc=N+1;double rr[5]={.049,.046,.043,.040,.037};double qq[5]={.925,.935,.945,.952,.960};for(int i=0;i<5&&es()<19.55;i++){ld(S);if(mk(rr[i],qq[i],19.80)){int c=cove();if(c<bc){bc=c;B=sv();ok=true;}}}if(ok&&bc<sc*93/100){ld(B);return true;}ld(S);return false;}}'''

def die(m): raise SystemExit('qvc_gen_abort: '+m)

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
    p=s.find('int main(){')
    if p<0: return None
    o=s.find('{',p); e=mb(s,o)
    if e<0: return None
    return p,e

def toks(s):
    r=[];i=0;n=len(s);bol=True;pp=False
    while i<n:
        c=s[i]
        if bol:
            j=i
            while j<n and s[j] in ' \t': j+=1
            pp=j<n and s[j]=='#';bol=False
        if c=='\n': bol=True;pp=False;i+=1;continue
        if pp: i+=1;continue
        if c in '"\'':
            q=c;i+=1
            while i<n:
                if s[i]=='\\': i+=2;continue
                if s[i]==q: i+=1;break
                i+=1
            continue
        if c=='/' and i+1<n and s[i+1]=='/':
            j=s.find('\n',i);i=n if j<0 else j;continue
        if c=='/' and i+1<n and s[i+1]=='*':
            j=s.find('*/',i+2);i=n if j<0 else j+2;continue
        m=ID.match(s,i)
        if m: r.append(m.group(0));i=m.end();continue
        i+=1
    return set(r)

def tr(src,mac):
    mp={b:a for a,b in mac if ' ' not in b};out=[];i=0;n=len(src);bol=True;pp=False
    rf=[a for a,b in mac if b=='return false'][0]
    rt=[a for a,b in mac if b=='return true'][0]
    while i<n:
        c=src[i]
        if bol:
            j=i
            while j<n and src[j] in ' \t': j+=1
            pp=j<n and src[j]=='#';bol=False
        if c=='\n': out.append(c);i+=1;bol=True;pp=False;continue
        if pp: out.append(c);i+=1;continue
        if c in '"\'':
            q=c;st=i;i+=1
            while i<n:
                if src[i]=='\\': i+=2;continue
                if src[i]==q: i+=1;break
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
                while k<n and src[k] in ' \t\r\n': k+=1
                m2=ID.match(src,k)
                if m2 and m2.group(0)=='false': out.append(rf);i=m2.end();continue
                if m2 and m2.group(0)=='true': out.append(rt);i=m2.end();continue
            out.append(mp.get(w,w));i=m.end();continue
        out.append(c);i+=1
    return ''.join(out)

def patch(src):
    need=['originalP','faces','AR','BU','BR','BF','Y','BE','CL','JD();','namespace W5','strong_validator','vps']
    if any(x not in src for x in need) or 'namespace QVC{' in src: return src,False
    mm=fm(src)
    if not mm: return src,False
    a,b=mm; m=src[a:b]; j=m.rfind('JD();')
    if j<0: return src,False
    nm=m[:j]+'QVC::run();'+m[j:]
    return src[:a]+LANE+nm+src[b:],True

def main():
    if len(sys.argv)<2: die('usage: python3 gen_qvc.py current_best.cpp [submit.cpp]')
    inp=Path(sys.argv[1])
    outp=Path(sys.argv[2]) if len(sys.argv)>2 else Path('qvc_submit.cpp')
    src=inp.read_text(encoding='utf-8')
    pat,changed=patch(src)
    variants=[
        [(f'O{i}',w) for i,w in enumerate(['double','static','const','inline','vector','unsigned','namespace','struct','template','typename'])]+[(f'P{i}',w) for i,w in enumerate(['return false','return true','continue','return','int','bool','void','class'])],
        [(f'Q{i}',w) for i,w in enumerate(['double','static','const','inline','vector','unsigned','namespace','struct','template','typename'])]+[(f'R{i}',w) for i,w in enumerate(['return false','return true','continue','return','int','bool','void','class'])],
        [(f'Z{i}',w) for i,w in enumerate(['double','static','const','inline','vector','unsigned','namespace','struct','template','typename'])]+[(f'Y{i}',w) for i,w in enumerate(['return false','return true','continue','return','int','bool','void','class'])]
    ]
    best=pat
    for mac in variants:
        names=[a for a,b in mac]
        if any(n in toks(pat) for n in names): continue
        u=pat.find('using namespace std;')
        if u<0: continue
        block=''.join('#define %s %s\n'%(a,b) for a,b in mac)
        cand=tr(pat[:u]+block+pat[u:],mac)
        if len(cand.encode())<len(best.encode()): best=cand
    if len(best.encode())>LIMIT:
        best=src
        changed=False
    outp.write_text(best,encoding='utf-8')
    print('wrote=%s'%outp)
    print('patched=%d'%int(changed and 'QVC::run()' in best))
    print('bytes=%d'%len(best.encode()))
    print('sha256=%s'%hashlib.sha256(best.encode()).hexdigest())

if __name__=='__main__':
    main()