#!/usr/bin/env python3
# compile generated: python3 r1_gen.py [base.cpp] r1.cpp && g++ -std=c++17 -O2 -pipe -static -s r1.cpp -o r1
# sample gate expectation: sample cube stays a no-op with first line 8 12; hidden fallback is selected 81.93-family base, R1 only commits after strong_validator()+vps guard.
import sys
from pathlib import Path
LIMIT=131072
LANE=r'''namespace R1{static Face mf(int a,int b,int c){Face f;f.v[0]=a;f.v[1]=b;f.v[2]=c;return f;}static unsigned long long ek(int a,int b){if(a>b)swap(a,b);return(unsigned long long)(unsigned int)a<<32|(unsigned int)b;}static void af(Face f){int id=(int)faces.size();faces.push_back(f);BR.push_back(1);++BE;M=(int)faces.size();Y[f.v[0]].push_back(id);Y[f.v[1]].push_back(id);Y[f.v[2]].push_back(id);}static bool ao(int a,int b,int c){return norm2(cross3(P[b]-P[a],P[c]-P[a]))>1e-28*CL*CL*CL*CL;}static void afo(int a,int b,int c,Vec3 r){Face f=mf(a,b,c);if(dot3(cross3(P[b]-P[a],P[c]-P[a]),r)<0)swap(f.v[1],f.v[2]);af(f);}struct G{double R,R2,C;Vec3 mn;vector<pair<long long,int>>H;vector<int>E;int ix(double x){return(int)floor((x-mn.x)/C);}long long ky(int x,int y,int z){return((long long)(x+1048576)<<42)^((long long)(y+1048576)<<21)^(long long)(z+1048576);}void init(double r){R=r;R2=r*r;C=max(r,1e-12*CL);mn={1e100,1e100,1e100};for(auto&p:originalP){mn.x=min(mn.x,p.x);mn.y=min(mn.y,p.y);mn.z=min(mn.z,p.z);}H.clear();E.clear();H.reserve(N);for(int i=0;i<N;i++)if(BU[i])H.push_back({ky(ix(P[i].x),ix(P[i].y),ix(P[i].z)),i});sort(H.begin(),H.end());}bool close(int id,const Vec3&p){return BU[id]&&norm2(P[id]-p)<=R2;}bool hit(const Vec3&p){int X=ix(p.x),Y0=ix(p.y),Z=ix(p.z);for(int z=Z-1;z<=Z+1;z++)for(int y=Y0-1;y<=Y0+1;y++)for(int x=X-1;x<=X+1;x++){auto k=ky(x,y,z);auto it=lower_bound(H.begin(),H.end(),make_pair(k,-1));for(;it!=H.end()&&it->first==k;++it)if(close(it->second,p))return 1;}for(int id:E)if(close(id,p))return 1;return 0;}void add(int id){E.push_back(id);if(E.size()>256){for(int q:E)H.push_back({ky(ix(P[q].x),ix(P[q].y),ix(P[q].z)),q});sort(H.begin(),H.end());E.clear();}}};static int nv(const Vec3&p,const vector<int>&B){int r=-1,st=max(1,(int)B.size()/7000);double d=1e100;for(int i=0;i<(int)B.size();i+=st){int v=B[i];if(!BU[v])continue;double e=norm2(P[v]-p);if(e<d){d=e;r=v;}}return r;}static bool pl(int id,const Vec3&p,const vector<int>&B){int v=nv(p,B);if(v<0)return 0;P[id]=p;for(int ff:Y[v])if(BR[ff]&&AC(ff,v)){Face o=faces[ff];Vec3 r=cross3(P[o.v[1]]-P[o.v[0]],P[o.v[2]]-P[o.v[0]]);if(norm2(r)<=1e-28*CL*CL*CL*CL)continue;if(!ao(o.v[0],o.v[1],id)||!ao(o.v[1],o.v[2],id)||!ao(o.v[2],o.v[0],id))continue;BR[ff]=0;--BE;BU[id]=1;BF[id]=0;Y[id].clear();afo(o.v[0],o.v[1],id,r);afo(o.v[1],o.v[2],id,r);afo(o.v[2],o.v[0],id,r);return 1;}return 0;}static bool cv(int cap,double rr,double lim){G g;g.init(rr*CL);vector<int>B,F;B.reserve(cove());F.reserve(N-cove());for(int i=0;i<N;i++)if(BU[i])B.push_back(i);else F.push_back(i);if(B.empty())return 0;int fp=0;for(int i=0;i<N;i++){if((i&4095)==0&&es()>lim)return 0;if(g.hit(originalP[i]))continue;while(fp<(int)F.size()&&BU[F[fp]])++fp;if(cove()+1>cap||fp>=(int)F.size())return 0;int id=F[fp++];if(!pl(id,originalP[i],B))return 0;B.push_back(id);g.add(id);}return 1;}static bool rd(int tgt,double h,double lim){AE p;p.AI=h*CL;p.BB=(.42*h+.004)*CL;p.BQ=cos((38.+520.*h)*acos(-1.)/180.);p.W=max(1e-12,1e-10*CL);p.AT=.12;p.AJ=max(1e-30,1e-24*CL*CL);p.AA=1;int s=cove();for(int it=0;it<3&&cove()>tgt&&es()<lim;it++){vector<unsigned long long>ks;ks.reserve((size_t)BE*3);for(int i=0;i<(int)faces.size();i++){if((i&8191)==0&&es()>lim)return cove()<s;if(!BR[i])continue;Face f=faces[i];if(!BU[f.v[0]]||!BU[f.v[1]]||!BU[f.v[2]])continue;ks.push_back(ek(f.v[0],f.v[1]));ks.push_back(ek(f.v[1],f.v[2]));ks.push_back(ek(f.v[2],f.v[0]));}sort(ks.begin(),ks.end());ks.erase(unique(ks.begin(),ks.end()),ks.end());vector<pair<double,unsigned long long>>E;E.reserve(ks.size());for(auto k:ks){int a=(int)(k>>32),b=(int)(k&0xffffffffu);if(a>=0&&b>=0&&a<N&&b<N&&BU[a]&&BU[b])E.push_back({norm2(P[a]-P[b]),k});}sort(E.begin(),E.end());for(auto&e:E){if(cove()<=tgt||es()>lim)break;int a=(int)(e.second>>32),b=(int)(e.second&0xffffffffu);if(a>=0&&b>=0&&a<N&&b<N&&BU[a]&&BU[b])GD(a,b,p);}}return cove()<s;}static bool tr(double tm,double h,double r,double cm,double q,double lim){AP S=AD();int b=cove(),t=max(64,(int)(b*tm));if(!rd(t,h,lim-.60)){rs(S);return 0;}int cap=max(80,(int)(b*cm));if(!cv(cap,r,lim-.25)){rs(S);return 0;}int a=cove();if(!(a<b&&a<=cap&&es()<lim&&W5::strong_validator()&&vps(512)>=q)){rs(S);return 0;}return 1;}static bool run(){if(es()>18.35||N<12000||N>125000)return 0;int b=cove();if(b<600||b>N*95/100)return 0;if(tr(.34,.110,.060,.58,.928,19.78))return 1;if(tr(.42,.095,.055,.66,.936,19.72))return 1;if(tr(.50,.080,.050,.74,.944,19.66))return 1;if(tr(.58,.065,.044,.82,.952,19.58))return 1;return 0;}}'''
CANDS=['submission_576_81.93_7.cpp','submission_601_81.93_7.cpp','submission_615_81.93_7.cpp','submission_612_81.93_7.cpp','submission_608_81.93_7.cpp','fetched_sources/kattis_19902774_81.93_7_external_tie.cpp','fetched_sources/kattis_19902761_81.93_7_external_tie.cpp','submission_563_81.93_7.cpp','fetched_sources/19901232.cpp','worker_outputs/tab4_x92_20260706/x92_breakthrough.cpp']
def die(s): raise SystemExit('R1 generator abort: '+s)
def match(s,o):
    d=0;q=None;esc=False;i=o
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
def find_main(s):
    p=s.find('int main(){')
    if p<0: die('missing compact int main anchor')
    o=s.find('{',p); e=match(s,o)
    if e<0: die('unmatched main')
    return p,e
def patch(s):
    for x in ['originalP','GD(','W5::strong_validator','vps(','cove','AD()','rs(','JD();']:
        if x not in s: die('missing required token '+x)
    if 'namespace R1{' in s: die('R1 already present')
    a,b=find_main(s); m=s[a:b]; j=m.rfind('JD();')
    if j<0: die('missing JD call in main')
    return s[:a]+LANE+m[:j]+'R1::run();'+m[j:]+s[b:]
def choose(explicit):
    tried=[]
    paths=[explicit] if explicit else [p for p in CANDS if Path(p).exists()]
    if not paths: die('no base source found; pass current 81.93-family source as argv[1]')
    for p in paths:
        try:
            s=Path(p).read_text()
            out=patch(s)
            n=len(out.encode())
            tried.append((p,n))
            if n<=LIMIT: return p,out,n
        except SystemExit as e:
            tried.append((p,str(e)))
            if explicit: raise
    die('all patched outputs exceed/fail limit: '+repr(tried))
def main():
    explicit=Path(sys.argv[1]) if len(sys.argv)>1 and sys.argv[1]!='-' else None
    outp=Path(sys.argv[2]) if len(sys.argv)>2 else Path('r1_connected_cover.cpp')
    base,out,n=choose(explicit)
    outp.write_text(out)
    print('base=%s'%base)
    print('wrote=%s bytes=%d limit=%d'%(outp,n,LIMIT))
    print('compile=g++ -std=c++17 -O2 -pipe -static -s %s -o r1_connected_cover'%outp)
    print('sample_gate=first line should remain 8 12 on cube/no-op; hidden signal: no-op tie if R1 rejects, otherwise guarded lower vertex count')
if __name__=='__main__': main()