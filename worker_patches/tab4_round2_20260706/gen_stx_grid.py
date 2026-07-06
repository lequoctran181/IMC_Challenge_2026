#!/usr/bin/env python3
# Usage:
#   python3 gen_stx_grid.py --src fetched_sources/kattis_19903326_fetched.cpp --out submission_stx_grid.cpp --exe submission_stx_grid.check
# The generator writes one C++17 file, compiles it, and runs the embedded sample gate.
import argparse, hashlib, os, re, string, subprocess, sys
from pathlib import Path

LIMIT = 131072
DEFAULTS = [
    'fetched_sources/kattis_19903326_fetched.cpp',
    'fetched_sources/kattis_19903326.cpp',
    'fetched_sources/kattis_19903153_81.93_7_worker_breakthrough_BOX3_from_543.cpp',
    'fetched_sources/kattis_19902206_81.93_7.cpp',
    'fetched_sources/kattis_19902388_81.93_7.cpp',
    'submission_608_81.93_7.cpp', 'submission_597_81.93_7.cpp',
    'submission_585_81.93_7.cpp', 'submission_563_81.93_7.cpp',
    'submission_543_81.93_7.cpp', 'fetched_sources/19901232.cpp'
]
REQ = [
    'W5::post_patch_pass();', 'static AP AD()', 'static void rs(const AP&s)',
    'visual_proxy_score', 'count_output_vertices_estimate', 'namespace W5',
    'static vector<Vec3>originalP', 'static vector<Face>AR', 'static double CL=1.'
]
SAMPLE = '''9 14
v 0.5 0.5 0.5
v 0.5 0.5 -0.5
v 0.5 -0.5 0.5
v 0.5 -0.5 -0.5
v -0.5 0.5 0.5
v -0.5 0.5 -0.5
v -0.5 -0.5 0.5
v -0.5 -0.5 -0.5
v 0.5 0.49 0.49
f 1 3 9
f 1 9 2
f 9 3 4
f 9 4 2
f 5 6 8
f 5 8 7
f 1 2 6
f 1 6 5
f 3 7 8
f 3 8 4
f 1 5 7
f 1 7 3
f 2 4 8
f 2 8 6
'''

LANE = r'''namespace STX{static Face mf(int a,int b,int c){Face f;f.v[0]=a;f.v[1]=b;f.v[2]=c;return f;}static bool adj(const int a[3],int m,int&b){for(int t=0;t<3;t++)for(int s=0;s<2;s++){int x=(a[t]-s+m)%m;bool ok=1;for(int i=0;i<3;i++){int d=(a[i]-x+m)%m;if(d!=0&&d!=1){ok=0;break;}}if(ok){b=x;return 1;}}return 0;}static bool adl(const int a[3],int m,int&b){for(int t=0;t<3;t++)for(int s=0;s<2;s++){int x=a[t]-s;if(x<0||x+1>=m)continue;bool ok=1;for(int i=0;i<3;i++){int d=a[i]-x;if(d!=0&&d!=1){ok=0;break;}}if(ok){b=x;return 1;}}return 0;}static bool con(int x,const vector<int>&v){return find(v.begin(),v.end(),x)!=v.end();}static void af(vector<Face>&F,const vector<Vec3>&X,Face f,const Vec3&r){Vec3 cr=cross3(X[f.v[1]]-X[f.v[0]],X[f.v[2]]-X[f.v[0]]);if(norm2(cr)<1e-28)return;if(dot3(cr,r)<0)swap(f.v[1],f.v[2]);F.pb(f);}static vector<Vec3> vn(){vector<Vec3>v(N,{0,0,0});for(const Face&f:AR){Vec3 cr=cross3(originalP[f.v[1]]-originalP[f.v[0]],originalP[f.v[2]]-originalP[f.v[0]]);v[f.v[0]]=v[f.v[0]]+cr;v[f.v[1]]=v[f.v[1]]+cr;v[f.v[2]]=v[f.v[2]]+cr;}return v;}static bool put(vector<Vec3>&X,vector<Face>&F,int base){if(X.empty()||F.empty()||(int)X.size()>=base||es()>18.72)return 0;AP B=AD();bool ok=0;if(AF(X,F)&&W5::strong_validator()&&cove()<base&&es()<18.82){int after=cove();double dr=(double)(base-after)/max(1,base),need=N>180000?.984:.993;if(dr>.45)need-=.006;else if(dr>.28)need-=.003;double p=vps(512);if(p>=need&&(p>.997||es()>18.58||vps(768)>=need-.003))ok=1;}if(!ok)rs(B);return ok;}static vector<int> cand(int n0,bool sph){map<int,int>mp;int st=max(1,M/120000);for(int i=0;i<M;i+=st){int a[3]={AR[i].v[0],AR[i].v[1],AR[i].v[2]};for(int k=0;k<3;k++){int d=abs(a[k]-a[(k+1)%3]);if(!d)continue;d=min(d,n0-d);if(d>=6&&d<=n0/3)mp[d]++;}}vector<pair<int,int>>q;for(auto&p:mp)q.pb({p.second,p.first});sort(q.rbegin(),q.rend());vector<int>r;auto add=[&](int s){if(s>=8&&s<=n0/3&&(sph?((N-2)%s==0):(N%s==0))&&!con(s,r))r.pb(s);};for(int i=0;i<(int)q.size()&&i<16;i++){int d=q[i].second;for(int e=-3;e<=3;e++)add(d+e);if(d)add(n0/d);}return r;}static bool ft(const Face&f,int S){if(S<8||N%S)return 0;int U=N/S;if(U<8)return 0;int r[3]={f.v[0]/S,f.v[1]/S,f.v[2]/S},c[3]={f.v[0]%S,f.v[1]%S,f.v[2]%S},ra=0,ca=0;if(!adj(r,U,ra)||!adj(c,S,ca))return 0;int m=0;for(int i=0;i<3;i++){int x=(r[i]-ra+U)%U,y=(c[i]-ca+S)%S;if(x>1||y>1)return 0;m|=1<<(x*2+y);}return __builtin_popcount((unsigned)m)==3;}static bool gt(int S){int st=max(1,M/180000),tot=0,ok=0;for(int i=0;i<M;i+=st){tot++;ok+=ft(AR[i],S);if((tot&8191)==0&&es()>17.25)return 0;}return tot>400&&ok*10000>=tot*9997;}static void mt(int S,int U2,int S2,const vector<Vec3>&n,vector<Vec3>&X,vector<Face>&F){int U=N/S;X.clear();F.clear();X.reserve(U2*S2);F.reserve(2*U2*S2);vector<int>src(U2*S2);for(int i=0;i<U2;i++){int oi=(long long)i*U/U2;for(int j=0;j<S2;j++){int oj=(long long)j*S/S2,id=oi*S+oj;src[i*S2+j]=id;X.pb(originalP[id]);}}auto id=[&](int i,int j){return((i+U2)%U2)*S2+((j+S2)%S2);};for(int i=0;i<U2;i++)for(int j=0;j<S2;j++){int a=id(i,j),b=id(i+1,j),c=id(i+1,j+1),d=id(i,j+1);af(F,X,mf(a,b,d),n[src[a]]+n[src[b]]+n[src[d]]);af(F,X,mf(b,c,d),n[src[b]]+n[src[c]]+n[src[d]]);}}static bool tor(int base){if(M!=2*N)return 0;vector<int>ss=cand(N,0);int S=0;for(int s:ss){if(gt(s)){S=s;break;}}if(!S)return 0;int U=N/S;vector<Vec3>n=vn();int tg1[]={4096,6144,8192,10240,12288,14336,16384,20480,24576,0};int tg2[]={12288,16384,24576,32768,49152,65536,0};int*tg=N>160000?tg2:tg1;for(int z=0;tg[z]&&es()<18.28;z++){int t=tg[z];double ar=sqrt((double)U/S);int U2=max(10,min(U,(int)(sqrt((double)t)*ar+.5)));int S2=max(10,min(S,t/max(1,U2)));int nv=U2*S2;if(nv>=base||nv>=N)continue;vector<Vec3>X;vector<Face>F;mt(S,U2,S2,n,X,F);if(put(X,F,base))return 1;}return 0;}static bool fs(const Face&f,int S){if(S<8||(N-2)%S)return 0;int R=(N-2)/S;if(R<3)return 0;int pole=0,pv=-1,o[3],oc=0;for(int i=0;i<3;i++){int v=f.v[i];if(v==0||v==N-1){pole++;pv=v;}else o[oc++]=v-1;}if(pole>1)return 0;if(pole==1){if(oc!=2)return 0;int r0=o[0]/S,r1=o[1]/S,c[3]={o[0]%S,o[1]%S,o[1]%S},b;if(pv==0&&r0==0&&r1==0&&adj(c,S,b))return 1;if(pv==N-1&&r0==R-1&&r1==R-1&&adj(c,S,b))return 1;return 0;}int rr[3],cc[3];for(int i=0;i<3;i++){int v=f.v[i]-1;rr[i]=v/S;cc[i]=v%S;}int rb=0,cb=0;if(!adl(rr,R,rb)||!adj(cc,S,cb))return 0;int m=0;for(int i=0;i<3;i++){int x=rr[i]-rb,y=(cc[i]-cb+S)%S;if(x>1||y>1)return 0;m|=1<<(x*2+y);}return __builtin_popcount((unsigned)m)==3;}static bool gs(int S){if(M!=2*(N-2))return 0;int st=max(1,M/180000),tot=0,ok=0;for(int i=0;i<M;i+=st){tot++;ok+=fs(AR[i],S);if((tot&8191)==0&&es()>17.25)return 0;}return tot>400&&ok*10000>=tot*9996;}static void ms(int S,int R2,int S2,const vector<Vec3>&n,vector<Vec3>&X,vector<Face>&F){int R=(N-2)/S;X.clear();F.clear();X.reserve(2+R2*S2);F.reserve(2*S2*R2);vector<int>src(2+R2*S2);src[0]=0;X.pb(originalP[0]);for(int i=0;i<R2;i++){int oi=R2==1?R/2:(long long)i*(R-1)/(R2-1);for(int j=0;j<S2;j++){int oj=(long long)j*S/S2,id=1+oi*S+oj;src[1+i*S2+j]=id;X.pb(originalP[id]);}}int bot=X.size();src.pb(N-1);X.pb(originalP[N-1]);auto id=[&](int i,int j){return 1+i*S2+((j%S2+S2)%S2);};for(int j=0;j<S2;j++)af(F,X,mf(0,id(0,j+1),id(0,j)),n[0]+n[src[id(0,j+1)]]+n[src[id(0,j)]]);for(int i=0;i+1<R2;i++)for(int j=0;j<S2;j++){int a=id(i,j),b=id(i,j+1),c=id(i+1,j),d=id(i+1,j+1);af(F,X,mf(a,b,c),n[src[a]]+n[src[b]]+n[src[c]]);af(F,X,mf(b,d,c),n[src[b]]+n[src[d]]+n[src[c]]);}for(int j=0;j<S2;j++)af(F,X,mf(bot,id(R2-1,j),id(R2-1,j+1)),n[N-1]+n[src[id(R2-1,j)]]+n[src[id(R2-1,j+1)]]);}static bool sph(int base){if(M!=2*(N-2))return 0;vector<int>ss=cand(N-2,1);int S=0;for(int s:ss){if(gs(s)){S=s;break;}}if(!S)return 0;int R=(N-2)/S;vector<Vec3>n=vn();int tg1[]={4098,6146,8194,10242,12290,14338,16386,20482,0};int tg2[]={12290,16386,24578,32770,49154,65538,0};int*tg=N>160000?tg2:tg1;for(int z=0;tg[z]&&es()<18.28;z++){int t=tg[z]-2;double ar=sqrt((double)R/S);int R2=max(3,min(R,(int)(sqrt((double)t)*ar+.5)));int S2=max(8,min(S,t/max(1,R2)));int nv=2+R2*S2;if(nv>=base||nv>=N)continue;vector<Vec3>X;vector<Face>F;ms(S,R2,S2,n,X,F);if(put(X,F,base))return 1;}return 0;}static bool run(){if(N<1500||N>1100000||es()>16.95)return 0;if(N>180000&&es()>13.65)return 0;if(AS&&(BG>.075||AL<.32||AH>.58))return 0;int base=cove();if(base<500||base>=N)return 0;AP B=AD();bool ok=0;if(M==2*N)ok=tor(base);if(!ok&&M==2*(N-2))ok=sph(base);if(!ok)rs(B);return ok;}}'''

def die(msg):
    print('STX generator abort: ' + msg, file=sys.stderr)
    sys.exit(2)

def choose_src(arg):
    if arg:
        p = Path(arg)
        if not p.exists():
            die('explicit --src not found: ' + arg)
        return p
    for s in DEFAULTS:
        p = Path(s)
        if p.exists():
            return p
    hits = []
    for pat in ('*19903326*.cpp', '*81.93*_7*.cpp', 'submission_563_81.93_7.cpp', 'submission_543_81.93_7.cpp', '19901232.cpp'):
        hits += list(Path('.').rglob(pat))[:25]
    seen = []
    for p in hits:
        if p not in seen:
            seen.append(p)
    for p in seen:
        try:
            x = p.read_text(errors='ignore')
        except Exception:
            continue
        if all(r in x for r in REQ[:6]) and len(x.encode()) <= LIMIT:
            return p
    die('no current-best source found; pass --src fetched_sources/kattis_19903326_fetched.cpp')

def find_main(src):
    m = re.search(r'\bint\s+main\s*\(\s*\)\s*\{', src)
    if not m:
        die('main() not found')
    o = src.find('{', m.start())
    i = o + 1
    d = 1
    q = None
    esc = False
    n = len(src)
    while i < n:
        c = src[i]
        if q:
            if esc:
                esc = False
            elif c == '\\':
                esc = True
            elif c == q:
                q = None
        else:
            if c in '"\'':
                q = c
            elif c == '/' and i + 1 < n and src[i+1] == '/':
                j = src.find('\n', i + 2)
                i = n if j < 0 else j
            elif c == '/' and i + 1 < n and src[i+1] == '*':
                j = src.find('*/', i + 2)
                i = n if j < 0 else j + 1
            elif c == '{':
                d += 1
            elif c == '}':
                d -= 1
                if d == 0:
                    return m.start(), o, i
        i += 1
    die('main brace did not close')

def patch_main(src):
    ms, o, e = find_main(src)
    body = src[o+1:e]
    if 'STX::run()' in body or 'namespace STX{' in src:
        die('source already has STX')
    if 'JC();GN();' not in body:
        die('main is not current-best family: missing JC();GN();')
    jd = body.rfind('JD();')
    if jd < 0:
        die('main missing final JD();')
    anchors = [
        'B16::R(39000,60000,220,-7,192,.96,18.05);',
        'B16::R(39000,60000,220,-7,192,96,18.05);'
    ]
    pos = -1
    for a in anchors:
        pos = body.find(a)
        if pos >= 0 and pos < jd:
            break
    if pos >= 0:
        nb = body[:pos] + 'if(!STX::run()){' + body[pos:jd] + '}' + body[jd:]
    else:
        pp = body.find('W5::post_patch_pass();')
        if pp >= 0 and pp < jd:
            pp += len('W5::post_patch_pass();')
            nb = body[:pp] + 'STX::run();' + body[pp:]
        else:
            nb = body[:jd] + 'STX::run();' + body[jd:]
    return src[:ms] + LANE + 'int main(){' + nb + '}' + src[e+1:]

OPS3 = ('<<=', '>>=', '->*', '...')
OPS2 = ('++', '--', '->', '&&', '||', '<<', '>>', '<=', '>=', '==', '!=', '+=', '-=', '*=', '/=', '%=', '&=', '|=', '^=', '::', '##', '.*')
OPC = set('+-*&|<>=:*/.%^!#')
ID_RE = re.compile(r'[A-Za-z_][A-Za-z0-9_]*')

def lex(s):
    tok = []
    i = 0
    n = len(s)
    while i < n:
        c = s[i]
        if c.isspace():
            j = i + 1
            while j < n and s[j].isspace():
                j += 1
            tok.append(('ws', s[i:j]))
            i = j
            continue
        if c == '/' and i + 1 < n and s[i+1] == '/':
            j = s.find('\n', i + 2)
            tok.append(('com', s[i:n if j < 0 else j]))
            i = n if j < 0 else j
            continue
        if c == '/' and i + 1 < n and s[i+1] == '*':
            j = s.find('*/', i + 2)
            if j < 0:
                die('unterminated comment')
            tok.append(('com', s[i:j+2]))
            i = j + 2
            continue
        if c == 'R' and i + 1 < n and s[i+1] == '"':
            mm = re.match(r'R"([ -~]{0,16})\(', s[i:])
            if mm:
                d = mm.group(1)
                end = ')' + d + '"'
                j = s.find(end, i + len(mm.group(0)))
                if j < 0:
                    die('unterminated raw string')
                tok.append(('lit', s[i:j+len(end)]))
                i = j + len(end)
                continue
        if c in '"\'':
            q = c
            j = i + 1
            esc = False
            while j < n:
                ch = s[j]
                if esc:
                    esc = False
                elif ch == '\\':
                    esc = True
                elif ch == q:
                    j += 1
                    break
                j += 1
            tok.append(('lit', s[i:j]))
            i = j
            continue
        if c.isalpha() or c == '_':
            j = i + 1
            while j < n and (s[j].isalnum() or s[j] == '_'):
                j += 1
            tok.append(('id', s[i:j]))
            i = j
            continue
        if c.isdigit() or (c == '.' and i + 1 < n and s[i+1].isdigit()):
            j = i + 1
            while j < n and (s[j].isalnum() or s[j] in '._+-'):
                if s[j] in '+-' and not (j > i and s[j-1] in 'eEpP'):
                    break
                j += 1
            tok.append(('num', s[i:j]))
            i = j
            continue
        if i + 2 < n and s[i:i+3] in OPS3:
            tok.append(('op', s[i:i+3]))
            i += 3
            continue
        if i + 1 < n and s[i:i+2] in OPS2:
            tok.append(('op', s[i:i+2]))
            i += 2
            continue
        tok.append(('op', c))
        i += 1
    return tok

def need(a, b):
    if not a or not b:
        return False
    x = a[-1]
    y = b[0]
    if (x.isalnum() or x == '_') and (y.isalnum() or y == '_'):
        return True
    if x == '.' and y == '.':
        return True
    if x in OPC and y in OPC:
        return True
    return False

def assemble(tok):
    out = []
    prev = ''
    j = 0
    line = True
    while j < len(tok):
        k, v = tok[j]
        if k in ('com', 'ws'):
            j += 1
            continue
        if line and v == '#':
            pp = [v]
            j += 1
            while j < len(tok):
                kk, vv = tok[j]
                if kk == 'ws' and '\n' in vv:
                    break
                if kk != 'com':
                    pp.append(vv)
                j += 1
            out.append(''.join(pp).rstrip() + '\n')
            prev = '\n'
            line = True
            while j < len(tok) and tok[j][0] == 'ws':
                j += 1
            continue
        if need(prev, v):
            out.append(' ')
        out.append(v)
        prev = v
        line = False
        j += 1
    return ''.join(out)

def minify(src):
    src = re.sub(r'^(?:#include\s*<[^>\n]+>\s*\n)+', '#include<bits/stdc++.h>\n', src, count=1)
    src = re.sub(r'\bnullptr\b', '0', src)
    src = re.sub(r'\binline\s+', '', src)
    src = src.replace('if(0&&"B16P515A"){}', '')
    src = src.replace('0.0', '.0').replace('1.0', '1.')
    tok = lex(src)
    kw = set('alignas alignof and and_eq asm auto bitand bitor bool break case catch char char16_t char32_t class compl const constexpr const_cast continue decltype default delete do double dynamic_cast else enum explicit export extern false float for friend goto if inline int long mutable namespace new noexcept not not_eq nullptr operator or or_eq private protected public register reinterpret_cast return short signed sizeof static static_assert static_cast struct switch template this thread_local throw true try typedef typeid typename union unsigned using virtual void volatile wchar_t while xor xor_eq'.split())
    std = set('abort abs acos array atan2 back begin cbrt ceil chrono clear cos count data deque empty end erase exit fabs fill find floor fprintf fread fwrite greater hypot insert int16_t int32_t int64_t int8_t isfinite less lower_bound make_pair map max memcpy memset min move pair pop pop_back pow priority_queue printf push push_back queue reserve resize reverse set setvbuf sin size size_t snprintf sort sqrt stable_sort stderr stdin stdout strtod strtol string swap tuple uint16_t uint32_t uint64_t uint8_t unordered_map unordered_set unique upper_bound vector'.split())
    protect = set(kw) | std | {'main'}
    for ln in src.splitlines():
        if ln.lstrip().startswith('#'):
            protect.update(re.findall(r'[A-Za-z_]\w*', ln))
    ids = [v for k, v in tok if k == 'id']
    idset = set(ids)
    for p, (k, v) in enumerate(tok):
        if k != 'id':
            continue
        q = p - 1
        while q >= 0 and tok[q][0] in ('ws', 'com'):
            q -= 1
        if q >= 0 and tok[q][1] in ('.', '->', '::'):
            protect.add(v)
        q = p + 1
        while q < len(tok) and tok[q][0] in ('ws', 'com'):
            q += 1
        if q < len(tok) and tok[q][1] == '::':
            protect.add(v)
    freq = {}
    for x in ids:
        if x not in protect and len(x) >= 6:
            freq[x] = freq.get(x, 0) + 1
    abc = string.ascii_letters
    def gen():
        k = 0
        while True:
            x = k
            a = abc[x % 52]
            x //= 52
            if x == 0:
                yield '_' + a
            else:
                digs = string.ascii_letters + string.digits + '_'
                r = ''
                while x:
                    r = digs[x % 63] + r
                    x //= 63
                yield '_' + a + r
            k += 1
    items = [((len(x)-3)*c, x, c) for x, c in freq.items() if c >= 2 and (len(x) >= 8 or c >= 5)]
    items.sort(reverse=True)
    mp = {}
    used = set(idset) | kw
    g = gen()
    saved = 0
    for _, x, c in items:
        if saved >= 26000:
            break
        y = next(g)
        while y in used:
            y = next(g)
        if len(y) >= len(x):
            continue
        mp[x] = y
        used.add(y)
        saved += (len(x) - len(y)) * c
    for i, (k, v) in enumerate(tok):
        if k == 'id' and v in mp:
            tok[i] = (k, mp[v])
    out = assemble(tok)
    if len(out.encode()) <= LIMIT:
        return out
    return macro_pack(out)

def macro_pack(src):
    tok = lex(src)
    idset = {v for k, v in tok if k == 'id'}
    words = ['double', 'static', 'return', 'vector', 'const', 'int', 'bool', 'void', 'false', 'true', 'unsigned']
    names = []
    k = 0
    while len(names) < len(words):
        nm = 'Q' + string.ascii_letters[k % 52]
        k += 1
        if nm not in idset:
            names.append(nm)
    mp = dict(zip(words, names))
    for i, (kind, val) in enumerate(tok):
        if kind == 'id' and val in mp:
            tok[i] = (kind, mp[val])
    body = assemble(tok)
    inc = '#include<bits/stdc++.h>\n'
    if body.startswith(inc):
        body = body[len(inc):]
    defs = ''.join('#define %s %s\n' % (mp[w], w) for w in words)
    return inc + defs + body

def validate_generated(code):
    for x in ('namespace STX', 'STX::run()', 'W5::strong_validator()', 'visual_proxy_score'):
        if x not in code:
            die('post-generation validation lost token: ' + x)
    if len(code.encode()) > LIMIT:
        die('generated source too large: %d > %d' % (len(code.encode()), LIMIT))

def run_sample(exe):
    cmd = exe if os.sep in exe else './' + exe
    r = subprocess.run([cmd], input=SAMPLE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, timeout=10)
    if r.returncode != 0:
        if r.stderr:
            print(r.stderr, file=sys.stderr, end='')
        die('sample run failed with exit ' + str(r.returncode))
    first = r.stdout.splitlines()[0].strip() if r.stdout.splitlines() else ''
    if first != '8 12':
        print(r.stdout[:500], file=sys.stderr)
        die('sample gate expected first line 8 12, got ' + repr(first))

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--src')
    ap.add_argument('--out', default='submission_stx_grid.cpp')
    ap.add_argument('--exe', default='submission_stx_grid.check')
    ap.add_argument('--cxx', default=os.environ.get('CXX', 'g++'))
    ap.add_argument('--no-compile', action='store_true')
    ap.add_argument('--no-sample', action='store_true')
    args = ap.parse_args()
    src_path = choose_src(args.src)
    src = src_path.read_text()
    sz = len(src.encode())
    if sz > LIMIT:
        die('source already over Kattis limit: %d' % sz)
    for r in REQ:
        if r not in src:
            die('missing current-best anchor: ' + r)
    patched = patch_main(src)
    out = minify(patched)
    validate_generated(out)
    op = Path(args.out)
    op.write_text(out)
    print('source=' + str(src_path))
    print('input_bytes=' + str(sz))
    print('output=' + str(op))
    print('output_bytes=' + str(len(out.encode())))
    print('sha256=' + hashlib.sha256(out.encode()).hexdigest())
    if not args.no_compile:
        cmd = [args.cxx, '-std=c++17', '-O2', '-pipe', '-static', '-s', str(op), '-o', args.exe]
        print('compile_command=' + ' '.join(cmd))
        r = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        if r.stdout:
            print(r.stdout, end='')
        if r.stderr:
            print(r.stderr, end='', file=sys.stderr)
        if r.returncode != 0:
            die('compile failed with exit ' + str(r.returncode))
        print('compile_ok=' + args.exe)
        if not args.no_sample:
            run_sample(args.exe)
            print('sample_gate=ok first_line_8_12')

if __name__ == '__main__':
    main()