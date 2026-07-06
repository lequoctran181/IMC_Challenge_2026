#!/usr/bin/env python3
import re,sys,hashlib
from pathlib import Path
LIMIT=131052
MAIN='int main(){JC();GN();if(!W2G::run())W2C::run();W5::post_patch_pass();VIMP::run();MIDEC::run();WK::run();B16::R(39000,60000,220,-7,192,.96,18.05);if(0&&"B16P515A"){}B16::R(39000,60000,76,-10,192,.96,18.35);for(int i=2;i--&&N>47500&N<6e4&&es()<18.6;)WK::run();if(N<50625&&es()<18.9)WK::run();JD();}'
NEWM='int main(){JC();GN();if(!W2G::run())W2C::run();W5::post_patch_pass();VIMP::run();MIDEC::run();KCM::run();WK::run();B16::R(39000,60000,220,-7,192,.96,18.05);if(0&&"B16P515A"){}B16::R(39000,60000,76,-10,192,.96,18.35);for(int i=2;i--&&N>47500&N<6e4&&es()<18.6;)WK::run();if(N<50625&&es()<18.9)WK::run();JD();}'
REQ=['namespace B16{','namespace W2G{','namespace W2C{','namespace WK{','namespace W5{','namespace VIMP{','namespace MIDEC{']
FORBID=['namespace KCM{','namespace KREV{','namespace KSV{','namespace BX8{','namespace UVH{','WD232','DX4','WITNESS_TET','connected_stellate','HBH']
INS=r'''namespace KCM{static Face F(int a,int b,int c){Face f;f.v[0]=a;f.v[1]=b;f.v[2]=c;return f;}static void Af(vector<Vec3>&V,vector<Face>&Q,int a,int b,int c,const Vec3&o){Face f=F(a,b,c);Vec3 cr=cross3(V[b]-V[a],V[c]-V[a]),ct=(V[a]+V[b]+V[c])*(1./3.);if(dot3(cr,ct-o)<0)swap(f.v[1],f.v[2]);Q.pb(f);}static Vec3 cen(){Vec3 mn=originalP[0],mx=mn;for(auto&p:originalP){mn.x=min(mn.x,p.x);mn.y=min(mn.y,p.y);mn.z=min(mn.z,p.z);mx.x=max(mx.x,p.x);mx.y=max(mx.y,p.y);mx.z=max(mx.z,p.z);}return(mn+mx)*.5;}static int fc(const Vec3&q,double&u,double&v){double ax=fabs(q.x),ay=fabs(q.y),az=fabs(q.z),m=max(ax,max(ay,az));if(m<=0){u=v=0;return 0;}if(ax>=ay&&ax>=az){u=q.y/ax;v=q.z/ax;return q.x>=0?0:1;}if(ay>=az){u=q.x/ay;v=q.z/ay;return q.y>=0?2:3;}u=q.x/az;v=q.y/az;return q.z>=0?4:5;}static bool H(const vector<Vec3>&V){double r=.049*CL,r2=r*r,c=max(r,1e-9);Vec3 mn=originalP[0],mx=mn;for(auto&p:originalP){mn.x=min(mn.x,p.x);mn.y=min(mn.y,p.y);mn.z=min(mn.z,p.z);mx.x=max(mx.x,p.x);mx.y=max(mx.y,p.y);mx.z=max(mx.z,p.z);}int nx=max(1,(int)((mx.x-mn.x)/c)+1),ny=max(1,(int)((mx.y-mn.y)/c)+1),nz=max(1,(int)((mx.z-mn.z)/c)+1);auto ix=[&](double x,int n,double m){int q=(int)((x-m)/c);return q<0?0:(q>=n?n-1:q);};auto key=[&](int x,int y,int z){return(z*ny+y)*nx+x;};vector<vector<int>>B((size_t)nx*ny*nz);for(int i=0;i<(int)V.size();i++)B[key(ix(V[i].x,nx,mn.x),ix(V[i].y,ny,mn.y),ix(V[i].z,nz,mn.z))].pb(i);auto nearV=[&](const Vec3&p){int X=ix(p.x,nx,mn.x),Y=ix(p.y,ny,mn.y),Z=ix(p.z,nz,mn.z);for(int z=Z-1;z<=Z+1;z++)if(z>=0&&z<nz)for(int y=Y-1;y<=Y+1;y++)if(y>=0&&y<ny)for(int x=X-1;x<=X+1;x++)if(x>=0&&x<nx)for(int i:B[key(x,y,z)])if(norm2(V[i]-p)<=r2)return true;return false;};int st=max(1,N/260000);for(int i=0;i<N;i+=st)if(!nearV(originalP[i]))return false;B.assign((size_t)nx*ny*nz,{});for(int i=0;i<N;i++)B[key(ix(originalP[i].x,nx,mn.x),ix(originalP[i].y,ny,mn.y),ix(originalP[i].z,nz,mn.z))].pb(i);auto nearO=[&](const Vec3&p){int X=ix(p.x,nx,mn.x),Y=ix(p.y,ny,mn.y),Z=ix(p.z,nz,mn.z);for(int z=Z-1;z<=Z+1;z++)if(z>=0&&z<nz)for(int y=Y-1;y<=Y+1;y++)if(y>=0&&y<ny)for(int x=X-1;x<=X+1;x++)if(x>=0&&x<nx)for(int i:B[key(x,y,z)])if(norm2(originalP[i]-p)<=r2)return true;return false;};for(auto&p:V)if(!nearO(p))return false;return true;}static bool det(const Vec3&o){if(N<140000||AS&&BG>.012)return false;int R=22,C[6]={0};vector<int>seen(6*R*R,0);int st=max(1,N/180000);for(int i=0;i<N;i+=st){double u,v;int f=fc(originalP[i]-o,u,v);int x=min(R-1,max(0,(int)((u+1)*.5*(R-1)+.5))),y=min(R-1,max(0,(int)((v+1)*.5*(R-1)+.5))),id=(f*R+y)*R+x;if(!seen[id])seen[id]=1,++C[f];}for(int f=0;f<6;f++)if(C[f]<R*R*55/100)return false;return true;}static bool build(int R,vector<Vec3>&V,vector<Face>&Q){Vec3 o=cen();if(!det(o))return false;vector<int>B[6];for(int f=0;f<6;f++)B[f].assign(R*R,-1);vector<double>S[6];for(int f=0;f<6;f++)S[f].assign(R*R,-1e100);for(int i=0;i<N;i++){Vec3 q=originalP[i]-o;double u,v;int f=fc(q,u,v);int x=min(R-1,max(0,(int)((u+1)*.5*(R-1)+.5))),y=min(R-1,max(0,(int)((v+1)*.5*(R-1)+.5))),id=y*R+x;double s=norm2(q);if(s>S[f][id])S[f][id]=s,B[f][id]=i;}vector<int>mp(R*R*R,-1);auto pt=[&](int f,int i,int j){int X,Y,Z;if(f==0)X=R-1,Y=i,Z=j;else if(f==1)X=0,Y=i,Z=j;else if(f==2)X=i,Y=R-1,Z=j;else if(f==3)X=i,Y=0,Z=j;else if(f==4)X=i,Y=j,Z=R-1;else X=i,Y=j,Z=0;int kk=(Z*R+Y)*R+X;if(mp[kk]>=0)return mp[kk];double best=1e100;int bi=-1,rad=0;for(;rad<5&&bi<0;rad++)for(int y=j-rad;y<=j+rad;y++)if(y>=0&&y<R)for(int x=i-rad;x<=i+rad;x++)if(x>=0&&x<R){int id=B[f][y*R+x];if(id>=0){double e=(x-i)*(x-i)+(y-j)*(y-j)-1e-9*S[f][y*R+x];if(e<best){best=e;bi=id;}}}if(bi<0)bi=N/2;int id=(int)V.size();V.pb(originalP[bi]);mp[kk]=id;return id;};V.clear();Q.clear();for(int f=0;f<6;f++)for(int i=0;i+1<R;i++)for(int j=0;j+1<R;j++){int a=pt(f,i,j),b=pt(f,i+1,j),c=pt(f,i+1,j+1),d=pt(f,i,j+1);if(a==b||a==c||a==d||b==c||b==d||c==d)return false;Af(V,Q,a,b,c,o);Af(V,Q,a,c,d,o);}return V.size()<N&&H(V);}static bool T(int R,double need){vector<Vec3>V;vector<Face>Q;if(!build(R,V,Q))return false;AP S=AD();bool ok=false;if(AF(V,Q)&&W5::strong_validator()&&es()<19.6){double p=vps(256);ok=p>=need&&cove()*100<N*(N>400000?10:16);if(ok&&N<500000&&es()<19.85)ok=vps(512)>=need-.008;}if(!ok)rs(S);return ok;}static bool run(){if(N<140000||es()>18.25)return false;AP S=AD();if(T(N>400000?30:34,N>400000?.918:.930))return true;rs(S);if(es()<18.8&&T(N>400000?38:42,N>400000?.932:.944))return true;rs(S);return false;}}'''
MAC=[('O0','double'),('O1','static'),('O2','const'),('O3','inline'),('O4','vector'),('O5','unsigned'),('O6','namespace'),('O7','struct'),('O8','template'),('O9','typename'),('P0','return false'),('P1','return true'),('P2','continue'),('P3','return'),('P4','int'),('P5','bool'),('P6','void'),('P7','class')]
ID=re.compile(r'[A-Za-z_][A-Za-z0-9_]*')
def die(m): raise SystemExit('workerK cubemap-shell abort: '+m)
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
    outp=Path(sys.argv[2]) if len(sys.argv)>2 else Path('workerK_19901232_cubemap_shell.cpp')
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
    if out.count('main(){')!=1 or 'KCM::run()' not in out: die('post-transform anchor lost')
    outp.write_text(out,encoding='utf-8')
    print(f'wrote {outp} bytes={n} sha256={hashlib.sha256(out.encode()).hexdigest()}')
if __name__=='__main__': main()