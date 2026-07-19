// SPDX-License-Identifier: MIT
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include "strict_mesh_io.hpp"
using namespace std;

using V3=meshio::V3;
using F=meshio::F;
using Mesh=meshio::Mesh;
static inline V3 operator-(V3 a,V3 b){return{a.x-b.x,a.y-b.y,a.z-b.z};}
static inline V3 crossp(V3 a,V3 b){return{a.y*b.z-a.z*b.y,a.z*b.x-a.x*b.z,a.x*b.y-a.y*b.x};}
static inline double dotp(V3 a,V3 b){return a.x*b.x+a.y*b.y+a.z*b.z;}
static inline double norm(V3 a){return sqrt(dotp(a,a));}

struct R{vector<double>d;vector<V3>n;vector<unsigned char>fg;};

static inline void proj(const V3&p,int view,int res,double&u,double&v,double&z){
    constexpr double D=2.5;
    const double f=800.0*((double)res/1024.0), c=0.5*res;
    double sx,sy;
    if(view==0){sx=p.y;sy=p.z;z=D-p.x;}
    else if(view==1){sx=-p.y;sy=p.z;z=D+p.x;}
    else if(view==2){sx=-p.x;sy=p.z;z=D-p.y;}
    else if(view==3){sx=p.x;sy=p.z;z=D+p.y;}
    else if(view==4){sx=p.x;sy=p.y;z=D-p.z;}
    else{sx=-p.x;sy=p.y;z=D+p.z;}
    u=f*sx/z+c;v=f*sy/z+c;
}

static void tri(R&rm,int res,V3 a,V3 b,V3 c,V3 un,int view){
    double x0,y0,z0,x1,y1,z1,x2,y2,z2;
    proj(a,view,res,x0,y0,z0);proj(b,view,res,x1,y1,z1);proj(c,view,res,x2,y2,z2);
    if(z0<=0||z1<=0||z2<=0)return;
    int xmin=max(0,(int)floor(min({x0,x1,x2})-.5)),xmax=min(res-1,(int)ceil(max({x0,x1,x2})+.5));
    int ymin=max(0,(int)floor(min({y0,y1,y2})-.5)),ymax=min(res-1,(int)ceil(max({y0,y1,y2})+.5));
    if(xmin>xmax||ymin>ymax)return;
    double den=(y1-y2)*(x0-x2)+(x2-x1)*(y0-y2);
    if(fabs(den)<1e-12)return;
    for(int yy=ymin;yy<=ymax;yy++){
        double py=yy+.5;
        for(int xx=xmin;xx<=xmax;xx++){
            double px=xx+.5;
            double w0=((y1-y2)*(px-x2)+(x2-x1)*(py-y2))/den;
            double w1=((y2-y0)*(px-x2)+(x0-x2)*(py-y2))/den;
            double w2=1-w0-w1;
            if(w0<-1e-9||w1<-1e-9||w2<-1e-9)continue;
            double zp=1/(w0/z0+w1/z1+w2/z2);
            int id=yy*res+xx;
            if(zp<rm.d[id]){rm.d[id]=zp;rm.n[id]=un;rm.fg[id]=1;}
        }
    }
}

static R render(const Mesh&m,int view,int res){
    R r;r.d.assign(res*res,255);r.n.assign(res*res,{0,0,0});r.fg.assign(res*res,0);
    for(auto f:m.f){
        V3 a=m.p[f.a],b=m.p[f.b],c=m.p[f.c];
        V3 cr=crossp(b-a,c-a);
        double l=norm(cr);
        if(!(l>0)){fprintf(stderr,"internal error: validated face became degenerate\n");exit(4);}
        tri(r,res,a,b,c,{cr.x/l,cr.y/l,cr.z/l},view);
    }
    return r;
}

static inline double normal_val(const V3&n,int ch){
    if(ch==0)return(n.x+1)*127.5;
    if(ch==1)return(n.y+1)*127.5;
    return(n.z+1)*127.5;
}

template<class G>
static double ssim(const R&a,const R&b,const vector<unsigned char>&fg,int res,G get){
    const int rad=5;
    const double c1=(.01*255)*(.01*255),c2=(.03*255)*(.03*255);
    const int stride=res+1;
    const size_t sz=(size_t)stride*stride;
    vector<double> sx(sz),sy(sz),sxx(sz),syy(sz),sxy(sz);
    for(int y=1;y<=res;y++){
        double ax=0,ay=0,axx=0,ayy=0,axy=0;
        for(int x=1;x<=res;x++){
            int p=(y-1)*res+x-1;
            double vx=get(a,p),vy=get(b,p);
            ax+=vx;ay+=vy;axx+=vx*vx;ayy+=vy*vy;axy+=vx*vy;
            size_t id=(size_t)y*stride+x,up=id-stride;
            sx[id]=sx[up]+ax;sy[id]=sy[up]+ay;
            sxx[id]=sxx[up]+axx;syy[id]=syy[up]+ayy;sxy[id]=sxy[up]+axy;
        }
    }
    auto rect=[&](const vector<double>&v,int x0,int y0,int x1,int y1){
        return v[(size_t)y1*stride+x1]-v[(size_t)y0*stride+x1]
             -v[(size_t)y1*stride+x0]+v[(size_t)y0*stride+x0];
    };
    double tot=0;int cnt=0;
    for(int y=rad;y<res-rad;y++)for(int x=rad;x<res-rad;x++){
        int id=y*res+x;
        if(!fg[id])continue;
        int x0=x-rad,y0=y-rad,x1=x+rad+1,y1=y+rad+1;
        constexpr double inv=1.0/121.0;
        double tx=rect(sx,x0,y0,x1,y1),ty=rect(sy,x0,y0,x1,y1);
        double ux=tx*inv,uy=ty*inv;
        double vx=max(0.0,rect(sxx,x0,y0,x1,y1)*inv-ux*ux);
        double vy=max(0.0,rect(syy,x0,y0,x1,y1)*inv-uy*uy);
        double cv=rect(sxy,x0,y0,x1,y1)*inv-ux*uy;
        tot+=((2*ux*uy+c1)*(2*cv+c2))/((ux*ux+uy*uy+c1)*(vx+vy+c2));
        cnt++;
    }
    return cnt?tot/cnt:1.0;
}

int main(int argc,char**argv){
    if(argc<3||argc>4){fprintf(stderr,"usage: %s original.obj simplified.obj [resolution=1024]\n",argv[0]);return 1;}
    try{
    int res=argc==4?meshio::parse_resolution_or_throw(argv[3]):1024;
    if(res!=1024) fprintf(stderr,"warning: non-1024 resolution is screening-only, not release evidence\n");
    Mesh a=meshio::read_mesh_or_throw(argv[1]),b=meshio::read_mesh_or_throw(argv[2]);
    double total=0;
    for(int v=0;v<6;v++){
        R ra=render(a,v,res),rb=render(b,v,res);
        vector<unsigned char>fg(res*res);
        for(int i=0;i<res*res;i++)fg[i]=ra.fg[i]||rb.fg[i];
        double ns=0;
        for(int ch=0;ch<3;ch++)ns+=ssim(ra,rb,fg,res,[ch](const R&r,int i){return normal_val(r.n[i],ch);});
        ns/=3;
        double ds=ssim(ra,rb,fg,res,[](const R&r,int i){return r.d[i];});
        total+=.5*ns+.5*ds;
    }
    printf("%.12f\n",total/6);
    }catch(const exception& e){fprintf(stderr,"error: %s\\n",e.what());return 2;}
}
