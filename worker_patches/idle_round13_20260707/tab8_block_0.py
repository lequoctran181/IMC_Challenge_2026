#!/usr/bin/env python3
# R1b generator.
# Usage:   python3 r1b_gen.py current_best_81_938904.cpp r1b.cpp
# Compile: g++ -std=c++17 -O2 -pipe -static -s r1b.cpp -o r1b
# Sample gate: sample cube remains a no-op; first output line should be exactly: 8 12
import re,sys
from pathlib import Path
from collections import Counter
LIMIT=131072
FULL=r'''namespace R1{static Face mf(int a,int b,int c){Face f;f.v[0]=a;f.v[1]=b;f.v[2]=c;return f;}static unsigned long long ek(int a,int b){if(a>b)swap(a,b);return(unsigned long long)(unsigned int)a<<32|(unsigned int)b;}static bool ar(int a,int b,int c){return norm2(cross3(P[b]-P[a],P[c]-P[a]))>1e-28*CL*CL*CL*CL;}static void adf(Face f){int id=(int)faces.size();faces.push_back(f);BR.push_back(1);++BE;M=(int)faces.size();Y[f.v[0]].push_back(id);Y[f.v[1]].push_back(id);Y[f.v[2]].push_back(id);}static void of(int a,int b,int c,Vec3 r){Face f=mf(a,b,c);if(dot3(cross3(P[b]-P[a],P[c]-P[a]),r)<0)swap(f.v[1],f.v[2]);adf(f);}struct H{double R,R2,C;Vec3 mn;vector<pair<long long,int>>T;vector<int>E;int ix(double x){return(int)floor((x-mn.x)/C);}int iy(double y){return(int)floor((y-mn.y)/C);}int iz(double z){return(int)floor((z-mn.z)/C);}long long ky(int x,int y,int z){return((long long)(x+1048576)<<42)^((long long)(y+1048576)<<21)^(long long)(z+1048576);}void init(double r){R=r;R2=r*r;C=max(r,1e-12*CL);mn={1e100,1e100,1e100};for(auto&p:originalP){mn.x=min(mn.x,p.x);mn.y=min(mn.y,p.y);mn.z=min(mn.z,p.z);}T.clear();E.clear();T.reserve(N);for(int i=0;i<N;i++)if(BU[i])T.push_back({ky(ix(P[i].x),iy(P[i].y),iz(P[i].z)),i});sort(T.begin(),T.end());}bool cl(int i,const Vec3&p){return BU[i]&&norm2(P[i]-p)<=R2;}bool hit(const Vec3&p){int X=ix(p.x),Y0=iy(p.y),Z=iz(p.z);for(int z=Z-1;z<=Z+1;z++)for(int y=Y0-1;y<=Y0+1;y++)for(int x=X-1;x<=X+1;x++){long long k=ky(x,y,z);auto it=lower_bound(T.begin(),T.end(),make_pair(k,-1));for(;it!=T.end()&&it->first==k;++it)if(cl(it->second,p))return 1;}for(int i:E)if(cl(i,p))return 1;return 0;}void add(int i){E.push_back(i);if(E.size()>255){for(int q:E)T.push_back({ky(ix(P[q].x),iy(P[q].y),iz(P[q].z)),q});sort(T.begin(),T.end());E.clear();}}};static int near(const Vec3&p,const vector<int>&B){int r=-1,st=max(1,(int)B.size()/5000);double d=1e100;for(int i=0;i<(int)B.size();i+=st){int v=B[i];if(!BU[v])continue;double e=norm2(P[v]-p);if(e<d){d=e;r=v;}}return r;}static bool plant(int id,const Vec3&p,const vector<int>&B){int v=near(p,B);if(v<0)return 0;P[id]=p;for(int ff:Y[v])if(BR[ff]&&AC(ff,v)){Face o=faces[ff];Vec3 r=cross3(P[o.v[1]]-P[o.v[0]],P[o.v[2]]-P[o.v[0]]);if(norm2(r)<=1e-28*CL*CL*CL*CL)continue;if(!ar(o.v[0],o.v[1],id)||!ar(o.v[1],o.v[2],id)||!ar(o.v[2],o.v[0],id))continue;BR[ff]=0;--BE;BU[id]=1;BF[id]=0;Y[id].clear();of(o.v[0],o.v[1],id,r);of(o.v[1],o.v[2],id,r);of(o.v[2],o.v[0],id,r);return 1;}return 0;}static bool cover(int cap,double rr,double lim){H h;h.init(rr*CL);vector<int>B,F;B.reserve(cove());F.reserve(N-cove());for(int i=0;i<N;i++)if(BU[i])B.push_back(i);else F.push_back(i);int fp=0,cn=(int)B.size();for(int i=0;i<N;i++){if((i&4095)==0&&es()>lim)return 0;if(h.hit(originalP[i]))continue;while(fp<(int)F.size()&&BU[F[fp]])++fp;if(cn+1>cap||fp>=(int)F.size())return 0;int id=F[fp++];if(!plant(id,originalP[i],B))return 0;B.push_back(id);h.add(id);++cn;}return 1;}static bool dec(int tgt,double h,double lim){AE p;p.AI=h*CL;p.BB=(.38*h+.004)*CL;p.BQ=cos((36.+520.*h)*acos(-1.)/180.);p.W=max(1e-12,1e-10*CL);p.AT=.10;p.AJ=max(1e-30,1e-24*CL*CL);p.AA=1;for(int it=0;it<2&&cove()>tgt&&es()<lim;it++){vector<unsigned long long>K;K.reserve((size_t)BE*3);for(int i=0;i<(int)faces.size();i++){if((i&8191)==0&&es()>lim)return 1;if(!BR[i])continue;Face f=faces[i];if(!BU[f.v[0]]||!BU[f.v[1]]||!BU[f.v[2]])continue;K.push_back(ek(f.v[0],f.v[1]));K.push_back(ek(f.v[1],f.v[2]));K.push_back(ek(f.v[2],f.v[0]));}sort(K.begin(),K.end());K.erase(unique(K.begin(),K.end()),K.end());vector<pair<double,unsigned long long>>E;E.reserve(K.size());for(auto k:K){int a=(int)(k>>32),b=(int)(k&0xffffffffu);if(a>=0&&b>=0&&a<N&&b<N&&BU[a]&&BU[b])E.push_back({norm2(P[a]-P[b]),k});}sort(E.begin(),E.end());for(auto&e:E){if(cove()<=tgt||es()>lim)break;int a=(int)(e.second>>32),b=(int)(e.second&0xffffffffu);if(a>=0&&b>=0&&a<N&&b<N&&BU[a]&&BU[b])GD(a,b,p);}}return 1;}static bool tr(double tm,double h,double r,double cm,double q,double lim){AP S=AD();int b=cove();if(!dec(max(64,(int)(b*tm)),h,lim-.65)){rs(S);return 0;}if(!cover(max(80,(int)(b*cm)),r,lim-.30)){rs(S);return 0;}int a=cove();if(a<b&&es()<lim&&W5::strong_validator()&&vps(512)>=q)return 1;rs(S);return 0;}static bool run(){if(es()>18.15||N<12000||N>125000)return 0;int b=cove();if(b<600||b>N*95/100)return 0;if(tr(.33,.115,.062,.58,.928,19.78))return 1;if(tr(.42,.095,.055,.67,.938,19.70))return 1;if(tr(.52,.075,.048,.78,.950,19.60))return 1;return 0;}}'''
LITE=r'''namespace R1{static unsigned long long ek(int a,int b){if(a>b)swap(a,b);return(unsigned long long)(unsigned int)a<<32|(unsigned int)b;}static bool run(){if(es()>18.6||N<12000||N>125000)return 0;AP S=AD();int b=cove();AE p;p.AI=.105*CL;p.BB=.050*CL;p.BQ=cos(52.*acos(-1.)/180.);p.W=max(1e-12,1e-10*CL);p.AT=.15;p.AJ=max(1e-30,1e-24*CL*CL);p.AA=1;vector<unsigned long long>K;K.reserve((size_t)BE*3);for(int i=0;i<(int)faces.size();i++){if((i&8191)==0&&es()>19.25)break;if(!BR[i])continue;Face f=faces[i];if(!BU[f.v[0]]||!BU[f.v[1]]||!BU[f.v[2]])continue;K.push_back(ek(f.v[0],f.v[1]));K.push_back(ek(f.v[1],f.v[2]));K.push_back(ek(f.v[2],f.v[0]));}sort(K.begin(),K.end());K.erase(unique(K.begin(),K.end()),K.end());vector<pair<double,unsigned long long>>E;E.reserve(K.size());for(auto k:K){int a=(int)(k>>32),c=(int)(k&0xffffffffu);if(a>=0&&c>=0&&a<N&&c<N&&BU[a]&&BU[c])E.push_back({norm2(P[a]-P[c]),k});}sort(E.begin(),E.end());int t=max(80,b*82/100);for(auto&e:E){if(cove()<=t||es()>19.45)break;int a=(int)(e.second>>32),c=(int)(e.second&0xffffffffu);if(a>=0&&c>=0&&a<N&&c<N&&BU[a]&&BU[c])GD(a,c,p);}if(cove()<b&&W5::strong_validator()&&vps(512)>=.970)return 1;rs(S);return 0;}}'''
NOOP=r'''namespace R1{static bool run(){return 0;}}'''
ID=re.compile(r'[A-Za-z_][A-Za-z0-9_]*')
def fail(m): raise SystemExit('r1b_gen: '+m)
def find_main(s):
    m=re.search(r'(^|[^A-Za-z0-9_])(?:signed\s+)?int\s+main\s*\([^)]*\)\s*\{',s)
    if not m: return None
    o=s.find('{',m.start())
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
                if d==0: return m.start(0)+(1 if m.group(1) else 0),i+1,o
        i+=1
    return None
def inspos(body):
    p=body.rfind('JD();')
    if p>=0: return p
    calls=list(re.finditer(r'\b[A-Za-z_][A-Za-z0-9_:]*\s*\(\s*\)\s*;',body))
    calls=[c for c in calls if 'R1::run' not in c.group(0)]
    if calls: return calls[-1].start()
    r=body.rfind('return')
    return r if r>=0 else len(body)-1
def has_lane_api(s):
    need=['AD()','rs(','strong_validator','vps','cove','originalP','BU','BR','faces','Y']
    return all(x in s for x in need) and re.search(r'\bGD\b',s) is not None
def inject(s,lane):
    fm=find_main(s)
    if not fm or not has_lane_api(s): return s,0
    a,b,o=fm
    main=s[a:b]
    rel=inspos(main)
    return s[:a]+lane+main[:rel]+'R1::run();'+main[rel:]+s[b:],1
def scan(s,replace=None,counts=None,alltok=None):
    out=[];i=0;n=len(s);bol=True;pp=False
    while i<n:
        c=s[i]
        if bol:
            j=i
            while j<n and s[j] in ' \t': j+=1
            pp=j<n and s[j]=='#'; bol=False
        if c=='\n':
            out.append(c); i+=1; bol=True; pp=False; continue
        if c in '"\'':
            q=c; st=i; i+=1
            while i<n:
                if s[i]=='\\': i+=2; continue
                if s[i]==q: i+=1; break
                i+=1
            out.append(s[st:i]); continue
        if not pp and c=='/' and i+1<n and s[i+1]=='/':
            j=s.find('\n',i); i=n if j<0 else j; continue
        if not pp and c=='/' and i+1<n and s[i+1]=='*':
            j=s.find('*/',i+2); i=n if j<0 else j+2; continue
        m=ID.match(s,i)
        if m:
            w=m.group(0)
            if alltok is not None: alltok.add(w)
            if not pp and counts is not None: counts[w]+=1
            out.append(replace.get(w,w) if (replace and not pp) else w)
            i=m.end(); continue
        out.append(c); i+=1
    return ''.join(out)
def names(used):
    first='QZXCVBNMLKJHGFDSAPOIUYTREWqzxcvbnmlkjhgfdsa'
    more='0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ'
    for a in first:
        for b in more:
            z=a+b
            if z not in used: yield z
    k=0
    while 1:
        z='Qq'+str(k)
        if z not in used: yield z
        k+=1
def add_macros(s,defs):
    if not defs: return s
    block=''.join('#define %s %s\n'%(a,b) for a,b in defs)
    pos=0
    for m in re.finditer(r'^[ \t]*#[^\n]*(?:\n|$)',s,re.M):
        if m.start()==pos: pos=m.end()
        elif m.start()<=pos: pos=max(pos,m.end())
    if pos==0:
        m=re.search(r'\busing\b|\btypedef\b|\bstruct\b|\bnamespace\b|\bstatic\b',s)
        pos=m.start() if m else 0
    return s[:pos]+block+s[pos:]
def compact(s):
    counts=Counter(); used=set()
    scan(s,counts=counts,alltok=used)
    gen=names(used)
    pairs=[]
    for w,c in counts.items():
        if c<3 or len(w)<4: continue
        if w.startswith('__'): continue
        pairs.append((c*(len(w)-2),w,c))
    pairs.sort(reverse=True)
    repl={}; defs=[]
    for _,w,c in pairs:
        nm=next(gen)
        gain=(len(w)-len(nm))*c-(9+len(nm)+len(w))
        if gain<=4: continue
        repl[w]=nm; defs.append((nm,w))
        if len(defs)>=420: break
    t=scan(s,replace=repl)
    return add_macros(t,defs)
def build(base):
    variants=[('full',FULL),('lite',LITE),('noop',NOOP)]
    logs=[]
    for name,lane in variants:
        s,ok=inject(base,lane)
        if not ok and name!='noop':
            logs.append('%s disabled: missing structural/API anchor'%name); continue
        raw=len(s.encode())
        c=compact(s)
        n=len(c.encode())
        logs.append('%s raw=%d compact=%d'%(name,raw,n))
        if n<=LIMIT: return c,name,logs
    c=compact(base)
    if len(c.encode())<=LIMIT: return c,'base_only',logs
    if len(base.encode())<=LIMIT: return base,'base_raw',logs
    fail('cannot fit any output; '+'; '.join(logs))
def pick():
    if len(sys.argv)>1 and sys.argv[1]!='-':
        p=Path(sys.argv[1])
        if not p.exists(): fail('missing input '+str(p))
        return p,p.read_text()
    cands=['submission_576_81.93_7.cpp','submission_601_81.93_7.cpp','submission_615_81.93_7.cpp','submission_612_81.93_7.cpp','submission_608_81.93_7.cpp','submission_563_81.93_7.cpp','fetched_sources/19901232.cpp','fetched_sources/kattis_19902774_81.93_7_external_tie.cpp','worker_outputs/tab4_x92_20260706/x92_breakthrough.cpp']
    for x in cands:
        p=Path(x)
        if p.exists(): return p,p.read_text()
    fail('pass the current-best source path as argv[1]')
def main():
    p,base=pick()
    outp=Path(sys.argv[2]) if len(sys.argv)>2 else Path('r1b.cpp')
    out,name,logs=build(base)
    outp.write_text(out)
    print('base=%s'%p)
    print('lane=%s'%name)
    print('bytes=%d limit=%d'%(len(out.encode()),LIMIT))
    for z in logs: print(z)
    print('compile=g++ -std=c++17 -O2 -pipe -static -s %s -o r1b'%outp)
    print('sample_gate=first line 8 12; hidden fallback unchanged unless R1 validator/vps commits')
if __name__=='__main__': main()