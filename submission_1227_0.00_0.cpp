#include <bits/stdc++.h>
using namespace std;

struct P{double x,y,z;};
static inline P operator+(const P&a,const P&b){return {a.x+b.x,a.y+b.y,a.z+b.z};}
static inline P operator-(const P&a,const P&b){return {a.x-b.x,a.y-b.y,a.z-b.z};}
static inline P operator*(const P&a,double s){return {a.x*s,a.y*s,a.z*s};}
static inline double dotp(const P&a,const P&b){return a.x*b.x+a.y*b.y+a.z*b.z;}
static inline P crossp(const P&a,const P&b){return {a.y*b.z-a.z*b.y,a.z*b.x-a.x*b.z,a.x*b.y-a.y*b.x};}
static inline double norm2(const P&a){return dotp(a,a);} 
static inline double dist2p(const P&a,const P&b){return norm2(a-b);} 

struct Face{int a,b,c; unsigned char on=1;};

static inline unsigned long long ekey(int a,int b){ if(a>b) swap(a,b); return (unsigned long long)(unsigned int)a<<32 | (unsigned int)b; }

struct Mesh{ vector<P> p; vector<Face> f; };

bool valid_edges(const Mesh& m){
    unordered_map<unsigned long long,int> cnt; cnt.reserve(m.f.size()*3+10);
    double scale=1.0;
    for(auto &q:m.p) scale=max(scale, fabs(q.x)+fabs(q.y)+fabs(q.z));
    double tiny=1e-28*scale*scale;
    for(auto &t:m.f){
        if(t.a<0||t.b<0||t.c<0||t.a>= (int)m.p.size()||t.b>= (int)m.p.size()||t.c>= (int)m.p.size()) return false;
        if(t.a==t.b||t.b==t.c||t.c==t.a) return false;
        if(norm2(crossp(m.p[t.b]-m.p[t.a],m.p[t.c]-m.p[t.a]))<=tiny) return false;
        cnt[ekey(t.a,t.b)]++; cnt[ekey(t.b,t.c)]++; cnt[ekey(t.c,t.a)]++;
    }
    for(auto &kv:cnt) if(kv.second!=2) return false;
    return !m.p.empty() && !m.f.empty();
}

struct I3{int x,y,z;};
static inline unsigned long long pack3(int x,int y,int z){return ((unsigned long long)(unsigned int)x<<42)^((unsigned long long)(unsigned int)y<<21)^(unsigned int)z;}
static inline I3 unpack3(unsigned long long k){return {(int)(k>>42),(int)((k>>21)&((1u<<21)-1)),(int)(k&((1u<<21)-1))};}

struct VoxelBuilder{
    const vector<P>& pts; P mn,mx; double L,D,eps;
    VoxelBuilder(const vector<P>&p,P mn_,P mx_,double D_,double eps_):pts(p),mn(mn_),mx(mx_),D(D_),eps(eps_){L=max({mx.x-mn.x,mx.y-mn.y,mx.z-mn.z}); if(L<=0)L=1;}
    bool build(int R,int budget,Mesh& out){
        double h=L/R;
        if(h*sqrt(3.0)>eps*0.985) return false;
        P org={mn.x-0.501*h,mn.y-0.501*h,mn.z-0.501*h};
        unordered_set<unsigned long long> occ; occ.reserve(pts.size()*2/3+1000);
        for(const P& q:pts){
            int ix=(int)floor((q.x-org.x)/h), iy=(int)floor((q.y-org.y)/h), iz=(int)floor((q.z-org.z)/h);
            if(ix<0)ix=0; if(iy<0)iy=0; if(iz<0)iz=0;
            occ.insert(pack3(ix,iy,iz));
        }
        if(occ.empty()) return false;
        // Very fragmented candidates are unlikely to be useful and too expensive.
        if((long long)occ.size()*8 > (long long)max(budget,1)*12) return false;
        unordered_set<unsigned long long> vis; vis.reserve(occ.size()*2+10);
        vector<unsigned long long> cells, stack; cells.reserve(4096); stack.reserve(4096);
        vector<P> vp; vector<Face> ff; vp.reserve(min((size_t)budget*2, occ.size()*8)); ff.reserve(occ.size()*12);
        const int dx[6]={-1,1,0,0,0,0},dy[6]={0,0,-1,1,0,0},dz[6]={0,0,0,0,-1,1};
        auto add_comp=[&](const vector<unsigned long long>& comp)->bool{
            unordered_map<unsigned long long,int> id; id.reserve(comp.size()*5+16);
            auto vid=[&](int x,int y,int z)->int{
                unsigned long long k=pack3(x,y,z); auto it=id.find(k); if(it!=id.end()) return it->second;
                int r=(int)vp.size(); id.emplace(k,r); vp.push_back({org.x+h*x,org.y+h*y,org.z+h*z}); return r;
            };
            unordered_set<unsigned long long> incomp; incomp.reserve(comp.size()*2+1); for(auto k:comp) incomp.insert(k);
            for(auto ck:comp){
                I3 c=unpack3(ck);
                for(int s=0;s<6;s++){
                    unsigned long long nk=pack3(c.x+dx[s],c.y+dy[s],c.z+dz[s]);
                    if(occ.find(nk)!=occ.end()) continue;
                    int A,B,C,Dv;
                    int x=c.x,y=c.y,z=c.z;
                    if(s==0){A=vid(x,y,z);B=vid(x,y,z+1);C=vid(x,y+1,z+1);Dv=vid(x,y+1,z);} 
                    else if(s==1){A=vid(x+1,y,z);B=vid(x+1,y+1,z);C=vid(x+1,y+1,z+1);Dv=vid(x+1,y,z+1);} 
                    else if(s==2){A=vid(x,y,z);B=vid(x+1,y,z);C=vid(x+1,y,z+1);Dv=vid(x,y,z+1);} 
                    else if(s==3){A=vid(x,y+1,z);B=vid(x,y+1,z+1);C=vid(x+1,y+1,z+1);Dv=vid(x+1,y+1,z);} 
                    else if(s==4){A=vid(x,y,z);B=vid(x,y+1,z);C=vid(x+1,y+1,z);Dv=vid(x+1,y,z);} 
                    else {A=vid(x,y,z+1);B=vid(x+1,y,z+1);C=vid(x+1,y+1,z+1);Dv=vid(x,y+1,z+1);} 
                    ff.push_back({A,B,C,1}); ff.push_back({A,C,Dv,1});
                }
            }
            return (int)vp.size() <= budget*3/2+100;
        };
        for(auto k:occ){
            if(vis.find(k)!=vis.end()) continue;
            cells.clear(); stack.clear(); stack.push_back(k); vis.insert(k);
            while(!stack.empty()){
                auto u=stack.back(); stack.pop_back(); cells.push_back(u); I3 c=unpack3(u);
                for(int s=0;s<6;s++){ auto nk=pack3(c.x+dx[s],c.y+dy[s],c.z+dz[s]); if(occ.find(nk)!=occ.end() && vis.insert(nk).second) stack.push_back(nk); }
            }
            if(!add_comp(cells)) return false;
            if((int)vp.size()>budget*3/2+100) return false;
        }
        Mesh m; m.p.swap(vp); m.f.swap(ff);
        if((int)m.p.size()>budget) return false;
        if(!valid_edges(m)) return false;
        out=move(m); return true;
    }
};

struct Quadric{
    double q[10];
    Quadric(){memset(q,0,sizeof q);} 
    void addPlane(double a,double b,double c,double d){
        q[0]+=a*a; q[1]+=a*b; q[2]+=a*c; q[3]+=a*d; q[4]+=b*b; q[5]+=b*c; q[6]+=b*d; q[7]+=c*c; q[8]+=c*d; q[9]+=d*d;
    }
    Quadric& operator+=(const Quadric&o){for(int i=0;i<10;i++)q[i]+=o.q[i]; return *this;}
    double eval(const P&p)const{
        double x=p.x,y=p.y,z=p.z;
        return q[0]*x*x+2*q[1]*x*y+2*q[2]*x*z+2*q[3]*x+q[4]*y*y+2*q[5]*y*z+2*q[6]*y+q[7]*z*z+2*q[8]*z+q[9];
    }
};

struct Item{double c; int a,b; bool operator<(Item const&o)const{return c>o.c;}};

struct Decimator{
    vector<P> p; vector<Face> f; vector<vector<int>> inc; vector<unsigned char> alive; vector<float> rad; vector<Quadric> Q;
    double eps,D,tinyArea; int aliveCnt; chrono::steady_clock::time_point st;
    Decimator(const Mesh&m,double eps_,double D_):p(m.p),f(m.f),eps(eps_),D(D_){
        int n=p.size(); inc.assign(n,{}); alive.assign(n,1); rad.assign(n,0); Q.assign(n,Quadric()); aliveCnt=n; tinyArea=1e-24*D*D*D*D;
        for(int i=0;i<(int)f.size();i++){
            inc[f[i].a].push_back(i); inc[f[i].b].push_back(i); inc[f[i].c].push_back(i);
            P nrm=crossp(p[f[i].b]-p[f[i].a],p[f[i].c]-p[f[i].a]); double l=sqrt(norm2(nrm)); if(l>0){nrm=nrm*(1.0/l); double d=-dotp(nrm,p[f[i].a]); Q[f[i].a].addPlane(nrm.x,nrm.y,nrm.z,d); Q[f[i].b].addPlane(nrm.x,nrm.y,nrm.z,d); Q[f[i].c].addPlane(nrm.x,nrm.y,nrm.z,d);}        
        }
    }
    inline bool fhas(int id,int v)const{return f[id].on&&(f[id].a==v||f[id].b==v||f[id].c==v);} 
    int edgeFaces(int a,int b,int ef[2],int op[2]){
        int c=0;
        for(int id:inc[a]) if(f[id].on && (f[id].a==b||f[id].b==b||f[id].c==b)){
            if(c<2){ef[c]=id; int x=f[id].a,y=f[id].b,z=f[id].c; op[c]=(x!=a&&x!=b)?x:((y!=a&&y!=b)?y:z);} c++;
            if(c>2) return c;
        }
        return c;
    }
    void neigh(int v, vector<int>&out){
        out.clear();
        for(int id:inc[v]) if(f[id].on){
            int a=f[id].a,b=f[id].b,c=f[id].c;
            if(a==v){ if(alive[b])out.push_back(b); if(alive[c])out.push_back(c); }
            else if(b==v){ if(alive[a])out.push_back(a); if(alive[c])out.push_back(c); }
            else if(c==v){ if(alive[a])out.push_back(a); if(alive[b])out.push_back(b); }
        }
        sort(out.begin(),out.end()); out.erase(unique(out.begin(),out.end()),out.end());
    }
    bool linkOK(int a,int b){
        int ef[2],op[2]; if(edgeFaces(a,b,ef,op)!=2) return false;
        vector<int>A,B,C; neigh(a,A); neigh(b,B); C.reserve(min(A.size(),B.size()));
        set_intersection(A.begin(),A.end(),B.begin(),B.end(),back_inserter(C));
        sort(op,op+2); C.erase(unique(C.begin(),C.end()),C.end());
        return C.size()==2 && C[0]==op[0] && C[1]==op[1];
    }
    bool bestDir(int a,int b,int&keep,int&rem,double&cost){
        if(!alive[a]||!alive[b]) return false;
        double d=sqrt(dist2p(p[a],p[b])); Quadric S=Q[a]; S+=Q[b];
        bool oka=max((double)rad[a],(double)rad[b]+d)<=eps, okb=max((double)rad[b],(double)rad[a]+d)<=eps;
        if(!oka&&!okb) return false;
        double ca=oka?S.eval(p[a])+1e-12*d*d:1e100, cb=okb?S.eval(p[b])+1e-12*d*d:1e100;
        if(ca<=cb){keep=a;rem=b;cost=ca;} else {keep=b;rem=a;cost=cb;}
        return true;
    }
    bool flipOK(int keep,int rem){
        for(int id:inc[rem]) if(f[id].on && !fhas(id,keep)){
            int a=f[id].a,b=f[id].b,c=f[id].c; P oldn=crossp(p[b]-p[a],p[c]-p[a]);
            if(a==rem)a=keep; else if(b==rem)b=keep; else if(c==rem)c=keep;
            if(a==b||b==c||c==a) return false;
            P newn=crossp(p[b]-p[a],p[c]-p[a]); double nn=norm2(newn), oo=norm2(oldn);
            if(nn<=tinyArea) return false;
            if(oo>0 && dotp(oldn,newn)<-1e-12*sqrt(oo*nn)) return false;
        }
        return true;
    }
    bool can(int keep,int rem){
        if(!alive[keep]||!alive[rem]) return false;
        double d=sqrt(dist2p(p[keep],p[rem])); if(max((double)rad[keep],(double)rad[rem]+d)>eps) return false;
        if(inc[keep].size()>512) cleanInc(keep); if(inc[rem].size()>512) cleanInc(rem);
        if(inc[keep].size()+inc[rem].size()>220) return false;
        if(!linkOK(keep,rem)) return false;
        return flipOK(keep,rem);
    }
    void cleanInc(int v){
        auto &x=inc[v];
        vector<int> y; y.reserve(x.size());
        for(int id:x) if(f[id].on && (f[id].a==v||f[id].b==v||f[id].c==v)) y.push_back(id);
        sort(y.begin(),y.end()); y.erase(unique(y.begin(),y.end()),y.end()); x.swap(y);
    }
    void collapse(int keep,int rem, priority_queue<Item>&pq){
        double d=sqrt(dist2p(p[keep],p[rem])); rad[keep]=(float)max((double)rad[keep],(double)rad[rem]+d); Q[keep]+=Q[rem]; alive[rem]=0; aliveCnt--;
        for(int id:inc[rem]) if(f[id].on){
            bool hk=(f[id].a==keep||f[id].b==keep||f[id].c==keep);
            if(f[id].a==rem) f[id].a=keep; if(f[id].b==rem) f[id].b=keep; if(f[id].c==rem) f[id].c=keep;
            if(hk||f[id].a==f[id].b||f[id].b==f[id].c||f[id].c==f[id].a){f[id].on=0;}
            else inc[keep].push_back(id);
        }
        inc[rem].clear(); cleanInc(keep);
        vector<int> nb; neigh(keep,nb);
        for(int w:nb){int k,r; double c; if(bestDir(keep,w,k,r,c)) pq.push({c,keep,w});}
    }
    Mesh run(int target){
        st=chrono::steady_clock::now();
        vector<unsigned long long> es; es.reserve(f.size()*3);
        for(auto &t:f){es.push_back(ekey(t.a,t.b)); es.push_back(ekey(t.b,t.c)); es.push_back(ekey(t.c,t.a));}
        sort(es.begin(),es.end()); es.erase(unique(es.begin(),es.end()),es.end());
        priority_queue<Item> pq;
        for(auto k:es){int a=(int)(k>>32),b=(int)(k&0xffffffffu),ke,re; double c; if(bestDir(a,b,ke,re,c)) pq.push({c,a,b});}
        long long pops=0;
        while(aliveCnt>target && !pq.empty()){
            Item it=pq.top(); pq.pop(); pops++;
            if((pops&4095)==0){ double sec=chrono::duration<double>(chrono::steady_clock::now()-st).count(); if(sec>19.2) break; }
            int a=it.a,b=it.b; if(a<0||b<0||a>= (int)p.size()||b>= (int)p.size()||!alive[a]||!alive[b]) continue;
            int keep,rem; double c; if(!bestDir(a,b,keep,rem,c)) continue;
            if(!can(keep,rem)) continue;
            collapse(keep,rem,pq);
        }
        vector<int> mp(p.size(),-1); vector<P> np; vector<Face> nf;
        for(auto &t:f) if(t.on){
            int vs[3]={t.a,t.b,t.c}; bool ok=true; for(int j=0;j<3;j++) if(!alive[vs[j]]) ok=false;
            if(!ok||vs[0]==vs[1]||vs[1]==vs[2]||vs[2]==vs[0]) continue;
            int ids[3]; for(int j=0;j<3;j++){ if(mp[vs[j]]<0){mp[vs[j]]=np.size(); np.push_back(p[vs[j]]);} ids[j]=mp[vs[j]]; }
            nf.push_back({ids[0],ids[1],ids[2],1});
        }
        return {move(np),move(nf)};
    }
};

void outputMesh(const Mesh&m){
    cout.setf(ios::fixed); cout<<setprecision(9);
    cout<<m.p.size()<<' '<<m.f.size()<<'\n';
    for(auto &p:m.p) cout<<p.x<<' '<<p.y<<' '<<p.z<<'\n';
    for(auto &t:m.f) cout<<t.a+1<<' '<<t.b+1<<' '<<t.c+1<<'\n';
}

int main(){
    ios::sync_with_stdio(false); cin.tie(nullptr);
    int V,F; if(!(cin>>V>>F)) return 0;
    Mesh in; in.p.resize(V); in.f.resize(F);
    P mn={1e100,1e100,1e100}, mx={-1e100,-1e100,-1e100};
    for(int i=0;i<V;i++){cin>>in.p[i].x>>in.p[i].y>>in.p[i].z; mn.x=min(mn.x,in.p[i].x); mn.y=min(mn.y,in.p[i].y); mn.z=min(mn.z,in.p[i].z); mx.x=max(mx.x,in.p[i].x); mx.y=max(mx.y,in.p[i].y); mx.z=max(mx.z,in.p[i].z);}    
    for(int i=0;i<F;i++){cin>>in.f[i].a>>in.f[i].b>>in.f[i].c; --in.f[i].a;--in.f[i].b;--in.f[i].c; in.f[i].on=1;}
    double D=sqrt((mx.x-mn.x)*(mx.x-mn.x)+(mx.y-mn.y)*(mx.y-mn.y)+(mx.z-mn.z)*(mx.z-mn.z)); if(D<=0)D=1;
    double eps=0.0492*D;
    int dream=max(4,(int)floor(V*0.0820));
    int soft=max(4,(int)floor(V*0.0950));
    Mesh best;
    // Breakout path: a vertex-set-safe voxel shell. It is selected only when it can be fine enough and still near top-score size.
    if(V>=90000){
        VoxelBuilder vb(in.p,mn,mx,D,eps);
        int Lc=(int)ceil(max({mx.x-mn.x,mx.y-mn.y,mx.z-mn.z})/(eps/sqrt(3.0)*0.985));
        int base=max(Lc,(int)floor(sqrt(max(1,soft)/5.2)));
        vector<int> cand;
        for(double mul: {1.35,1.22,1.10,1.00,0.92,0.84,0.76}) cand.push_back(max(Lc,(int)round(base*mul)));
        cand.push_back(Lc); sort(cand.rbegin(),cand.rend()); cand.erase(unique(cand.begin(),cand.end()),cand.end());
        for(int R:cand){ Mesh m; if(vb.build(R,soft,m)){ if((int)m.p.size()>=max(100,(int)(V*0.055))) { best=move(m); break; } } }
        if(!best.p.empty() && (int)best.p.size()<=soft && (int)best.p.size()<V){ outputMesh(best); return 0; }
    }
    // Conservative topology-preserving fallback with endpoint-only QEM collapses and explicit vertex-Hausdorff radius bound.
    int target=dream;
    Decimator dec(in,eps,D);
    Mesh m=dec.run(target);
    if((int)m.p.size()<V && valid_edges(m)){ outputMesh(m); return 0; }
    outputMesh(in); return 0;
}
