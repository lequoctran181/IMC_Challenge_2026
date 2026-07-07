#include <bits/stdc++.h>
using namespace std;

struct Vec3{double x,y,z;};
static inline Vec3 operator+(Vec3 a,Vec3 b){return {a.x+b.x,a.y+b.y,a.z+b.z};}
static inline Vec3 operator-(Vec3 a,Vec3 b){return {a.x-b.x,a.y-b.y,a.z-b.z};}
static inline Vec3 operator*(Vec3 a,double s){return {a.x*s,a.y*s,a.z*s};}
static inline Vec3 operator/(Vec3 a,double s){return {a.x/s,a.y/s,a.z/s};}
static inline Vec3& operator+=(Vec3& a,Vec3 b){a.x+=b.x;a.y+=b.y;a.z+=b.z;return a;}
static inline double dot3(Vec3 a,Vec3 b){return a.x*b.x+a.y*b.y+a.z*b.z;}
static inline Vec3 cross3(Vec3 a,Vec3 b){return {a.y*b.z-a.z*b.y,a.z*b.x-a.x*b.z,a.x*b.y-a.y*b.x};}
static inline double norm2(Vec3 a){return dot3(a,a);} static inline double norm3(Vec3 a){return sqrt(max(0.0,norm2(a)));}
static inline Vec3 unit(Vec3 a){double l=norm3(a);return l>1e-300?a/l:Vec3{0,0,0};}

struct FaceIn{int a,b,c;};
struct Face{int v[3]; unsigned char alive; Vec3 n;};
struct MVert{Vec3 p; double rad=0; int ver=0; unsigned char alive=1; vector<int> inc;};
struct ECand{double cost; int a,b,va,vb,keep; bool operator<(ECand const&o)const{return cost>o.cost;}};

static int N0=0,F0=0; static vector<Vec3> P0; static vector<FaceIn> Fin; static Vec3 Bmn,Bmx; static double diagLen=1,hausTol=1,hausTol2=1,areaFloor=1e-30; static clock_t T0;

static inline double elapsed(){return double(clock()-T0)/CLOCKS_PER_SEC;}
static inline unsigned long long ekey(int a,int b){if(a>b)swap(a,b);return (unsigned long long)(unsigned int)a<<32 | (unsigned int)b;}
static inline double triArea2(const vector<Vec3>&V,int a,int b,int c){return norm2(cross3(V[b]-V[a],V[c]-V[a]));}
static inline Vec3 faceN(Vec3 a,Vec3 b,Vec3 c){return unit(cross3(b-a,c-a));}

static vector<char> slurp(){vector<char> r; char buf[1<<16]; size_t n; while((n=fread(buf,1,sizeof(buf),stdin))>0)r.insert(r.end(),buf,buf+n); r.push_back(0); return r;}
static inline void sws(char*&p){while(*p==' '||*p=='\n'||*p=='\r'||*p=='\t')++p;}
static bool readInput(){
    auto buf=slurp(); char*p=buf.data(); N0=(int)strtol(p,&p,10); F0=(int)strtol(p,&p,10); if(N0<=0||F0<=0)return false;
    P0.resize(N0); Fin.reserve(F0); Bmn={1e100,1e100,1e100}; Bmx={-1e100,-1e100,-1e100};
    for(int i=0;i<N0;i++){sws(p); if(*p=='v'||*p=='V')++p; double x=strtod(p,&p),y=strtod(p,&p),z=strtod(p,&p); P0[i]={x,y,z}; Bmn.x=min(Bmn.x,x);Bmn.y=min(Bmn.y,y);Bmn.z=min(Bmn.z,z);Bmx.x=max(Bmx.x,x);Bmx.y=max(Bmx.y,y);Bmx.z=max(Bmx.z,z);} 
    for(int i=0;i<F0;i++){sws(p); if(*p=='f'||*p=='F')++p; int a=(int)strtol(p,&p,10)-1,b=(int)strtol(p,&p,10)-1,c=(int)strtol(p,&p,10)-1; if(a<0||b<0||c<0||a>=N0||b>=N0||c>=N0||a==b||b==c||a==c)return false; Fin.push_back({a,b,c});}
    diagLen=norm3(Bmx-Bmn); if(!(diagLen>0))diagLen=1; hausTol=0.05*diagLen*0.997; hausTol2=hausTol*hausTol; areaFloor=max(1e-32,diagLen*diagLen*1e-30); return true;
}
static void outputIdentity(){printf("%d %d\n",N0,F0); for(auto&p:P0)printf("v %.17g %.17g %.17g\n",p.x,p.y,p.z); for(auto&f:Fin)printf("f %d %d %d\n",f.a+1,f.b+1,f.c+1);} 

struct PGrid{
    Vec3 mn,mx; double cell=1,r2=1; int nx=1,ny=1,nz=1; vector<vector<int>> bin; const vector<Vec3>*Q=nullptr;
    int cl(int x,int n)const{return x<0?0:(x>=n?n-1:x);} int ix(double x)const{return cl((int)((x-mn.x)/cell),nx);} int iy(double y)const{return cl((int)((y-mn.y)/cell),ny);} int iz(double z)const{return cl((int)((z-mn.z)/cell),nz);} int key(int x,int y,int z)const{return (z*ny+y)*nx+x;}
    void init(const vector<Vec3>&A,double R){Q=&A;r2=R*R;mn=mx=A[0];for(auto&p:A){mn.x=min(mn.x,p.x);mn.y=min(mn.y,p.y);mn.z=min(mn.z,p.z);mx.x=max(mx.x,p.x);mx.y=max(mx.y,p.y);mx.z=max(mx.z,p.z);} cell=max(R,1e-12); nx=max(1,(int)((mx.x-mn.x)/cell)+1); ny=max(1,(int)((mx.y-mn.y)/cell)+1); nz=max(1,(int)((mx.z-mn.z)/cell)+1); if((long long)nx*ny*nz>1700000){nx=ny=nz=1;cell=max({mx.x-mn.x,mx.y-mn.y,mx.z-mn.z})+1;} bin.assign((size_t)nx*ny*nz,{}); for(int i=0;i<(int)A.size();i++)bin[key(ix(A[i].x),iy(A[i].y),iz(A[i].z))].push_back(i);} 
    bool near(Vec3 p)const{int X=ix(p.x),Y=iy(p.y),Z=iz(p.z); for(int z=Z-1;z<=Z+1;z++)if(z>=0&&z<nz)for(int y=Y-1;y<=Y+1;y++)if(y>=0&&y<ny)for(int x=X-1;x<=X+1;x++)if(x>=0&&x<nx)for(int id:bin[key(x,y,z)])if(norm2((*Q)[id]-p)<=r2)return true; return false;}
};

static bool validateAndPrint(const vector<Vec3>&V,const vector<array<int,3>>&Fout){
    if(V.empty()||Fout.empty()||(int)V.size()>N0)return false;
    PGrid gin,gout; gin.init(P0,hausTol*1.0005); gout.init(V,hausTol*1.0005);
    for(auto&p:V)if(!gin.near(p))return false; for(auto&p:P0)if(!gout.near(p))return false;
    unordered_map<unsigned long long,int> ec; ec.reserve(Fout.size()*4+10); vector<unsigned char> used(V.size(),0);
    for(auto&t:Fout){int a=t[0],b=t[1],c=t[2]; if(a<0||b<0||c<0||a>=(int)V.size()||b>=(int)V.size()||c>=(int)V.size()||a==b||b==c||a==c)return false; if(triArea2(V,a,b,c)<=areaFloor)return false; used[a]=used[b]=used[c]=1; ec[ekey(a,b)]++; ec[ekey(b,c)]++; ec[ekey(c,a)]++;}
    for(auto u:used)if(!u)return false; for(auto&kv:ec)if(kv.second!=2)return false;
    printf("%d %d\n",(int)V.size(),(int)Fout.size()); for(auto&p:V)printf("v %.17g %.17g %.17g\n",p.x,p.y,p.z); for(auto&t:Fout)printf("f %d %d %d\n",t[0]+1,t[1]+1,t[2]+1); return true;
}
static void addOriented(vector<Vec3>&V,vector<array<int,3>>&F,int a,int b,int c,Vec3 ref){Vec3 cr=cross3(V[b]-V[a],V[c]-V[a]); if(dot3(cr,ref)<0)swap(b,c); F.push_back({a,b,c});}

static bool tryEllipsoid(){
    Vec3 C=(Bmn+Bmx)*0.5; Vec3 R=(Bmx-Bmn)*0.5; if(R.x<=0||R.y<=0||R.z<=0)return false; double rmin=min({R.x,R.y,R.z}),rmax=max({R.x,R.y,R.z}); if(rmin<diagLen*1e-6)return false;
    double rms=0,mx=0; int cnt=0; for(int i=0;i<N0;i+=max(1,N0/120000)){Vec3 d=P0[i]-C; double s=sqrt((d.x*d.x)/(R.x*R.x)+(d.y*d.y)/(R.y*R.y)+(d.z*d.z)/(R.z*R.z)); if(s<=1e-12)return false; Vec3 q=C+d/s; double e=norm3(P0[i]-q); rms+=e*e; mx=max(mx,e); cnt++;} rms=sqrt(rms/max(1,cnt)); if(!(rms<hausTol*0.105 && mx<hausTol*0.45))return false;
    double step=hausTol*0.78; int nlon=max(12,(int)ceil(2*M_PI*rmax/step)); int nlat=max(7,(int)ceil(M_PI*rmax/step)); nlon=min(nlon,96); nlat=min(nlat,48); long long nv=2ll+(long long)(nlat-1)*nlon; if(nv>=N0*82ll/100||nv>N0)return false;
    vector<Vec3> V; vector<array<int,3>> F; V.reserve((size_t)nv); int north=0,south=1; V.push_back({C.x,C.y,C.z+R.z}); V.push_back({C.x,C.y,C.z-R.z}); vector<vector<int>> id(nlat-1, vector<int>(nlon));
    for(int i=1;i<nlat;i++){double th=M_PI*i/nlat; double st=sin(th),ct=cos(th); for(int j=0;j<nlon;j++){double ph=2*M_PI*j/nlon; id[i-1][j]=V.size(); V.push_back({C.x+R.x*st*cos(ph),C.y+R.y*st*sin(ph),C.z+R.z*ct});}}
    auto en=[&](Vec3 p){Vec3 d=p-C; return unit({d.x/(R.x*R.x),d.y/(R.y*R.y),d.z/(R.z*R.z)});};
    for(int j=0;j<nlon;j++){int a=id[0][j],b=id[0][(j+1)%nlon]; addOriented(V,F,north,a,b,en((V[north]+V[a]+V[b])/3.0));}
    for(int i=0;i<nlat-2;i++)for(int j=0;j<nlon;j++){int a=id[i][j],b=id[i][(j+1)%nlon],c=id[i+1][j],d=id[i+1][(j+1)%nlon]; addOriented(V,F,a,c,d,en((V[a]+V[c]+V[d])/3.0)); addOriented(V,F,a,d,b,en((V[a]+V[d]+V[b])/3.0));}
    for(int j=0;j<nlon;j++){int a=id[nlat-2][j],b=id[nlat-2][(j+1)%nlon]; addOriented(V,F,south,b,a,en((V[south]+V[a]+V[b])/3.0));}
    return validateAndPrint(V,F);
}

static Vec3 permSet(int ax, double u,double v,double w, Vec3 C){ if(ax==2)return {C.x+u,C.y+v,C.z+w}; if(ax==1)return {C.x+u,C.y+w,C.z+v}; return {C.x+w,C.y+u,C.z+v}; }
static bool tryTorus(){
    Vec3 C=(Bmn+Bmx)*0.5; double E[3]={(Bmx.x-Bmn.x)*0.5,(Bmx.y-Bmn.y)*0.5,(Bmx.z-Bmn.z)*0.5};
    for(int ax=0;ax<3;ax++){int a=(ax+1)%3,b=(ax+2)%3; double ea=E[a],eb=E[b],eh=E[ax]; if(ea<=0||eb<=0||eh<=0)continue; if(fabs(ea-eb)>0.10*max(ea,eb))continue; double r=eh,R=(ea+eb)*0.5-r; if(!(R>r*1.08&&r>diagLen*1e-5))continue; double rms=0,mx=0;int cnt=0; for(int i=0;i<N0;i+=max(1,N0/120000)){Vec3 d=P0[i]-C; double coord[3]={d.x,d.y,d.z}; double rho=sqrt(coord[a]*coord[a]+coord[b]*coord[b]); double q=sqrt((rho-R)*(rho-R)+coord[ax]*coord[ax]); double e=fabs(q-r); rms+=e*e; mx=max(mx,e); cnt++;} rms=sqrt(rms/max(1,cnt)); if(!(rms<hausTol*0.095&&mx<hausTol*0.42))continue;
        int nu=max(12,(int)ceil(2*M_PI*(R+r)/(hausTol*0.78))); int nv=max(8,(int)ceil(2*M_PI*r/(hausTol*0.78))); nu=min(nu,128); nv=min(nv,64); if((long long)nu*nv>=N0*82ll/100)continue; vector<Vec3> V; vector<array<int,3>> F; V.reserve(nu*nv); for(int i=0;i<nu;i++){double u=2*M_PI*i/nu,cu=cos(u),su=sin(u); for(int j=0;j<nv;j++){double vv=2*M_PI*j/nv,cv=cos(vv),sv=sin(vv); double x=(R+r*cv)*cu,y=(R+r*cv)*su,z=r*sv; V.push_back(permSet(ax,x,y,z,C));}}
        for(int i=0;i<nu;i++)for(int j=0;j<nv;j++){int i2=(i+1)%nu,j2=(j+1)%nv; int p00=i*nv+j,p10=i2*nv+j,p01=i*nv+j2,p11=i2*nv+j2; Vec3 ref=unit(V[p00]-C); addOriented(V,F,p00,p10,p11,ref); addOriented(V,F,p00,p11,p01,ref);} if(validateAndPrint(V,F))return true;
    } return false;
}

static bool tryCylinder(){
    Vec3 C=(Bmn+Bmx)*0.5; double E[3]={(Bmx.x-Bmn.x)*0.5,(Bmx.y-Bmn.y)*0.5,(Bmx.z-Bmn.z)*0.5};
    for(int ax=0;ax<3;ax++){int a=(ax+1)%3,b=(ax+2)%3; double h=E[ax],ra=E[a],rb=E[b]; if(h<=0||ra<=0||rb<=0||fabs(ra-rb)>0.08*max(ra,rb))continue; double r=(ra+rb)*0.5; if(h<r*0.7)continue; double rms=0,mx=0;int cnt=0; for(int i=0;i<N0;i+=max(1,N0/120000)){Vec3 d=P0[i]-C; double q[3]={d.x,d.y,d.z}; double rho=sqrt(q[a]*q[a]+q[b]*q[b]); double az=fabs(q[ax]); double e; if(az<=h&&rho<=r)e=min(fabs(rho-r),fabs(az-h)); else if(az<=h)e=fabs(rho-r); else if(rho<=r)e=fabs(az-h); else e=sqrt((rho-r)*(rho-r)+(az-h)*(az-h)); rms+=e*e; mx=max(mx,e);cnt++;} rms=sqrt(rms/max(1,cnt)); if(!(rms<hausTol*0.085&&mx<hausTol*0.36))continue;
        int seg=max(12,(int)ceil(2*M_PI*r/(hausTol*0.78))); int hz=max(1,(int)ceil(2*h/(hausTol*0.9))); int rr=max(1,(int)ceil(r/(hausTol*0.9))); seg=min(seg,160); hz=min(hz,120); rr=min(rr,80); long long nv=(long long)(hz+1)*seg + 2ll*max(0,rr-1)*seg + 2; if(nv>=N0*82ll/100)continue; vector<Vec3> V; vector<array<int,3>> F; vector<vector<int>> side(hz+1, vector<int>(seg)); auto pos=[&](double x,double y,double z){return permSet(ax,x,y,z,C);};
        for(int iz=0;iz<=hz;iz++){double z=-h+2*h*iz/hz; for(int s=0;s<seg;s++){double ph=2*M_PI*s/seg; side[iz][s]=V.size(); V.push_back(pos(r*cos(ph),r*sin(ph),z));}}
        for(int iz=0;iz<hz;iz++)for(int s=0;s<seg;s++){int s2=(s+1)%seg; Vec3 ref=unit(V[side[iz][s]]-pos(0,0,-h+2*h*iz/hz)); addOriented(V,F,side[iz][s],side[iz+1][s],side[iz+1][s2],ref); addOriented(V,F,side[iz][s],side[iz+1][s2],side[iz][s2],ref);} 
        for(int cap=0;cap<2;cap++){double z=cap?h:-h; Vec3 ref=permSet(ax,0,0,cap?1:-1,{0,0,0}); vector<vector<int>> ring(rr+1); ring[rr]=side[cap?hz:0]; int center=V.size(); V.push_back(pos(0,0,z)); ring[0]={center}; for(int k=1;k<rr;k++){double rad=r*k/rr; ring[k].resize(seg); for(int s=0;s<seg;s++){double ph=2*M_PI*s/seg; ring[k][s]=V.size(); V.push_back(pos(rad*cos(ph),rad*sin(ph),z));}} if(rr==1){for(int s=0;s<seg;s++){int s2=(s+1)%seg; addOriented(V,F,center,ring[rr][s],ring[rr][s2],ref);}} else {for(int s=0;s<seg;s++){int s2=(s+1)%seg; addOriented(V,F,center,ring[1][s],ring[1][s2],ref);} for(int k=1;k<rr;k++)for(int s=0;s<seg;s++){int s2=(s+1)%seg; addOriented(V,F,ring[k][s],ring[k+1][s],ring[k+1][s2],ref); addOriented(V,F,ring[k][s],ring[k+1][s2],ring[k][s2],ref);}}}
        if(validateAndPrint(V,F))return true;
    } return false;
}

static vector<MVert> MV; static vector<Face> MF; static int aliveV,aliveF; static vector<int> markV; static int mtok=1;
static bool mfHas(const Face&f,int v){return f.v[0]==v||f.v[1]==v||f.v[2]==v;} static int mfOther(const Face&f,int a,int b){for(int k=0;k<3;k++){int x=f.v[k];if(x!=a&&x!=b)return x;}return -1;}
static void mClean(int u){if(!MV[u].alive){vector<int>().swap(MV[u].inc);return;}vector<int>r;for(int id:MV[u].inc)if(id>=0&&id<(int)MF.size()&&MF[id].alive&&mfHas(MF[id],u))r.push_back(id);sort(r.begin(),r.end());r.erase(unique(r.begin(),r.end()),r.end());MV[u].inc.swap(r);} 
static bool mEdgeFaces(int u,int v,int ef[2],int op[2]){mClean(u);int c=0;for(int id:MV[u].inc){Face&f=MF[id];if(f.alive&&mfHas(f,v)){if(c>=2)return false;ef[c]=id;op[c]=mfOther(f,u,v);c++;}}return c==2&&op[0]>=0&&op[1]>=0&&op[0]!=op[1];}
static bool mLink(int u,int v,int o0,int o1){if(++mtok==INT_MAX){fill(markV.begin(),markV.end(),0);mtok=1;}for(int id:MV[u].inc){Face&f=MF[id];if(!f.alive)continue;for(int k=0;k<3;k++){int w=f.v[k];if(w!=u&&w!=v&&MV[w].alive)markV[w]=mtok;}}vector<int> com;for(int id:MV[v].inc){Face&f=MF[id];if(!f.alive)continue;for(int k=0;k<3;k++){int w=f.v[k];if(w==u||w==v||!MV[w].alive)continue;if(markV[w]==mtok){if(w!=o0&&w!=o1)return false;if(find(com.begin(),com.end(),w)==com.end())com.push_back(w);}}}return com.size()==2;}
static bool mNormalOK(int keep,int rem,int ef0,int ef1){double minDot=N0>50000?-0.05:(N0>15000?0.04:0.15); auto chk=[&](int id){if(id==ef0||id==ef1)return true;Face&f=MF[id];if(!f.alive)return true;int a=f.v[0],b=f.v[1],c=f.v[2];if(a==rem)a=keep;if(b==rem)b=keep;if(c==rem)c=keep;if(a==b||b==c||a==c)return false;double ar=triArea2(P0,a,b,c);if(ar<=areaFloor)return false;Vec3 n=faceN(P0[a],P0[b],P0[c]);return dot3(n,f.n)>minDot;};for(int id:MV[keep].inc)if(!chk(id))return false;for(int id:MV[rem].inc)if(!chk(id))return false;return true;}
static bool mMakeCand(int a,int b,ECand&c){if(a==b||!MV[a].alive||!MV[b].alive)return false;int ef[2],op[2];if(!mEdgeFaces(a,b,ef,op)||!mLink(a,b,op[0],op[1]))return false;double dab=norm3(MV[a].p-MV[b].p);int keep=-1,rem=-1;double ra=max(MV[a].rad,MV[b].rad+dab), rb=max(MV[b].rad,MV[a].rad+dab);if(ra<=hausTol*0.985){keep=a;rem=b;} if(rb<ra&&rb<=hausTol*0.985){keep=b;rem=a;} if(keep<0)return false;if(!mNormalOK(keep,rem,ef[0],ef[1]))return false;c={dab*dab+MV[rem].rad*0.01,a,b,MV[a].ver,MV[b].ver,keep};return true;}
static void mPushAround(int u,priority_queue<ECand>&pq){if(!MV[u].alive)return;mClean(u);if(++mtok==INT_MAX){fill(markV.begin(),markV.end(),0);mtok=1;}for(int id:MV[u].inc){Face&f=MF[id];if(!f.alive)continue;for(int k=0;k<3;k++){int w=f.v[k];if(w==u||!MV[w].alive||markV[w]==mtok)continue;markV[w]=mtok;ECand c;if(mMakeCand(min(u,w),max(u,w),c))pq.push(c);}}}
static bool mCollapse(ECand old,priority_queue<ECand>&pq){ECand c;if(!mMakeCand(old.a,old.b,c))return false;int keep=c.keep,rem=(keep==c.a?c.b:c.a);int ef[2],op[2];if(!mEdgeFaces(keep,rem,ef,op)||!mLink(keep,rem,op[0],op[1])||!mNormalOK(keep,rem,ef[0],ef[1]))return false;double d=norm3(MV[keep].p-MV[rem].p);MV[keep].rad=max(MV[keep].rad,MV[rem].rad+d);for(int i=0;i<2;i++)if(MF[ef[i]].alive){MF[ef[i]].alive=0;aliveF--;}for(int id:MV[rem].inc){Face&f=MF[id];if(!f.alive)continue;for(int k=0;k<3;k++)if(f.v[k]==rem)f.v[k]=keep;}MV[rem].alive=0;MV[rem].ver++;MV[keep].ver++;aliveV--;MV[keep].inc.insert(MV[keep].inc.end(),MV[rem].inc.begin(),MV[rem].inc.end());vector<int>().swap(MV[rem].inc);mClean(keep);for(int id:MV[keep].inc){Face&f=MF[id];if(f.alive)f.n=faceN(P0[f.v[0]],P0[f.v[1]],P0[f.v[2]]);}mPushAround(keep,pq);vector<int> around=MV[keep].inc;for(int id:around){Face&f=MF[id];if(f.alive)for(int k=0;k<3;k++)if(f.v[k]!=keep)mPushAround(f.v[k],pq);}return true;}
static bool tryCollapseFallback(){MV.assign(N0,{});MF.clear();MF.reserve(F0);for(int i=0;i<N0;i++){MV[i].p=P0[i];MV[i].alive=1;MV[i].rad=0;}for(auto&t:Fin){Face f;f.v[0]=t.a;f.v[1]=t.b;f.v[2]=t.c;f.alive=1;f.n=faceN(P0[t.a],P0[t.b],P0[t.c]);int id=MF.size();MF.push_back(f);MV[t.a].inc.push_back(id);MV[t.b].inc.push_back(id);MV[t.c].inc.push_back(id);}aliveV=N0;aliveF=F0;markV.assign(N0,0);priority_queue<ECand>pq;vector<unsigned long long>edges;edges.reserve((size_t)F0*3);for(auto&f:MF){edges.push_back(ekey(f.v[0],f.v[1]));edges.push_back(ekey(f.v[1],f.v[2]));edges.push_back(ekey(f.v[2],f.v[0]));}sort(edges.begin(),edges.end());edges.erase(unique(edges.begin(),edges.end()),edges.end());for(auto k:edges){int a=k>>32,b=k&0xffffffffu;ECand c;if(mMakeCand(a,b,c))pq.push(c);}double ratio=N0>100000?0.08:(N0>40000?0.12:(N0>12000?0.18:0.28));int target=max(12,(int)(N0*ratio));double lim=17.5;int pops=0;while(!pq.empty()&&aliveV>target&&elapsed()<lim){ECand c=pq.top();pq.pop();if(!MV[c.a].alive||!MV[c.b].alive||MV[c.a].ver!=c.va||MV[c.b].ver!=c.vb)continue;mCollapse(c,pq);if((++pops&4095)==0&&elapsed()>lim)break;}vector<int>id(N0,-1),rev;vector<Vec3> Vout;for(int i=0;i<N0;i++)if(MV[i].alive){id[i]=Vout.size();rev.push_back(i);Vout.push_back(P0[i]);}vector<array<int,3>>Fout;for(auto&f:MF)if(f.alive){int a=id[f.v[0]],b=id[f.v[1]],c=id[f.v[2]];if(a<0||b<0||c<0||a==b||b==c||a==c)return false;Fout.push_back({a,b,c});}if((int)Vout.size()>=N0)return false;return validateAndPrint(Vout,Fout);}

int main(){T0=clock(); if(!readInput())return 0; if(N0<=8||F0<=12){outputIdentity();return 0;} if(tryTorus())return 0; if(tryCylinder())return 0; if(tryEllipsoid())return 0; if(tryCollapseFallback())return 0; outputIdentity(); return 0;}