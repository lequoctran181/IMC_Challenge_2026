// SPDX-License-Identifier: MIT
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>
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

static unsigned char byte_value(double value){
    return (unsigned char)llround(max(0.0,min(255.0,value)));
}

static void write_ppm(const string& path,int res,const vector<unsigned char>&rgb){
    FILE* fp=fopen(path.c_str(),"wb");
    if(!fp) throw runtime_error("cannot write diagnostic image: "+path);
    fprintf(fp,"P6\n%d %d\n255\n",res,res);
    if(fwrite(rgb.data(),1,rgb.size(),fp)!=rgb.size()){
        fclose(fp);throw runtime_error("short write for diagnostic image: "+path);
    }
    fclose(fp);
}

static void dump_diagnostics(const R&a,const R&b,int res,const string&prefix){
    const size_t pixels=(size_t)res*res;
    vector<unsigned char>nr(pixels*3),nc(pixels*3),ne(pixels*3);
    vector<unsigned char>dr(pixels*3),dc(pixels*3),de(pixels*3);
    double depth_lo=1e100,depth_hi=-1e100,error_hi=0;
    for(size_t i=0;i<pixels;++i){
        if(a.fg[i]){depth_lo=min(depth_lo,a.d[i]);depth_hi=max(depth_hi,a.d[i]);}
        if(b.fg[i]){depth_lo=min(depth_lo,b.d[i]);depth_hi=max(depth_hi,b.d[i]);}
        if(a.fg[i]&&b.fg[i]) error_hi=max(error_hi,fabs(a.d[i]-b.d[i]));
    }
    if(!(depth_hi>depth_lo)){depth_lo=0;depth_hi=1;}
    if(!(error_hi>0)) error_hi=1;
    for(size_t i=0;i<pixels;++i){
        auto set_normal=[&](vector<unsigned char>&out,const R&r){
            if(!r.fg[i]){out[3*i]=out[3*i+1]=out[3*i+2]=245;return;}
            out[3*i]=byte_value((r.n[i].x+1)*127.5);
            out[3*i+1]=byte_value((r.n[i].y+1)*127.5);
            out[3*i+2]=byte_value((r.n[i].z+1)*127.5);
        };
        auto set_depth=[&](vector<unsigned char>&out,const R&r){
            if(!r.fg[i]){out[3*i]=out[3*i+1]=out[3*i+2]=245;return;}
            double value=255*(r.d[i]-depth_lo)/(depth_hi-depth_lo);
            out[3*i]=out[3*i+1]=out[3*i+2]=byte_value(value);
        };
        set_normal(nr,a);set_normal(nc,b);set_depth(dr,a);set_depth(dc,b);
        if(a.fg[i]&&b.fg[i]){
            double cosine=max(-1.0,min(1.0,dotp(a.n[i],b.n[i])));
            double t=acos(cosine)/3.14159265358979323846;
            ne[3*i]=byte_value(255*min(1.0,3*t));
            ne[3*i+1]=byte_value(255*(1-min(1.0,3*t)));
            ne[3*i+2]=40;
            double d=min(1.0,fabs(a.d[i]-b.d[i])/error_hi);
            de[3*i]=byte_value(255*d);de[3*i+1]=byte_value(210*(1-d));de[3*i+2]=35;
        }else if(a.fg[i]||b.fg[i]){
            ne[3*i]=255;ne[3*i+1]=0;ne[3*i+2]=255;
            de[3*i]=255;de[3*i+1]=0;de[3*i+2]=255;
        }else{
            ne[3*i]=ne[3*i+1]=ne[3*i+2]=245;
            de[3*i]=de[3*i+1]=de[3*i+2]=245;
        }
    }
    write_ppm(prefix+"_normal_reference.ppm",res,nr);
    write_ppm(prefix+"_normal_candidate.ppm",res,nc);
    write_ppm(prefix+"_normal_error.ppm",res,ne);
    write_ppm(prefix+"_depth_reference.ppm",res,dr);
    write_ppm(prefix+"_depth_candidate.ppm",res,dc);
    write_ppm(prefix+"_depth_error.ppm",res,de);
}

// The hidden Armadillo case is normalized by the contestant before its
// direct-coordinate payload is emitted.  Keep the standard evaluator
// byte-for-byte equivalent by default, but allow large local candidate
// inventories to be audited in that exact normalized frame without first
// materializing a transformed OBJ for every candidate.
static void normalize_like_arm_input(Mesh& reference,Mesh& candidate){
    V3 lo=reference.p.front(),hi=lo;
    for(const V3& p:reference.p){
        lo.x=min(lo.x,p.x);lo.y=min(lo.y,p.y);lo.z=min(lo.z,p.z);
        hi.x=max(hi.x,p.x);hi.y=max(hi.y,p.y);hi.z=max(hi.z,p.z);
    }
    V3 center{.5*(lo.x+hi.x),.5*(lo.y+hi.y),.5*(lo.z+hi.z)};
    double scale=0;
    for(const V3& p:reference.p)scale=max(scale,norm(p-center));
    auto apply=[&](Mesh& mesh){for(V3& p:mesh.p){
        p.x=round((p.x-center.x)/scale*1e9)/1e9;
        p.y=round((p.y-center.y)/scale*1e9)/1e9;
        p.z=round((p.z-center.z)/scale*1e9)/1e9;
    }};
    apply(reference);apply(candidate);
}

static inline void proj(const V3&p,int view,int res,double&u,double&v,double&z){
    constexpr double D=2.5;
    static const double base_f=[](){const char* s=getenv("FOCAL");return s?atof(s):800.0;}();
    const double f=base_f*((double)res/1024.0), c=0.5*res;
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
    if(getenv("GAUSSIAN")){
        // Diagnostic for the common 11x11, sigma=1.5 SSIM convention.  The
        // published checker currently specifies a uniform window, so this is
        // deliberately opt-in and leaves the default evaluator unchanged.
        double kernel[11],sum=0;
        for(int i=-rad;i<=rad;i++){kernel[i+rad]=exp(-.5*i*i/(1.5*1.5));sum+=kernel[i+rad];}
        for(double&w:kernel)w/=sum;
        const size_t area=(size_t)res*res;
        vector<double> raw[5],horizontal[5];
        for(int k=0;k<5;k++){raw[k].resize(area);horizontal[k].resize(area);}
        for(size_t p=0;p<area;p++){
            double x=get(a,p),y=get(b,p);
            raw[0][p]=x;raw[1][p]=y;raw[2][p]=x*x;raw[3][p]=y*y;raw[4][p]=x*y;
        }
        for(int y=0;y<res;y++)for(int x=rad;x<res-rad;x++){
            size_t p=(size_t)y*res+x;
            for(int k=0;k<5;k++)for(int d=-rad;d<=rad;d++)horizontal[k][p]+=kernel[d+rad]*raw[k][p+d];
        }
        double total=0;int count=0;
        for(int y=rad;y<res-rad;y++)for(int x=rad;x<res-rad;x++){
            size_t p=(size_t)y*res+x;if(!fg[p])continue;
            double s[5]={};
            for(int k=0;k<5;k++)for(int d=-rad;d<=rad;d++)s[k]+=kernel[d+rad]*horizontal[k][p+(size_t)d*res];
            double vx=max(0.0,s[2]-s[0]*s[0]),vy=max(0.0,s[3]-s[1]*s[1]),cv=s[4]-s[0]*s[1];
            total+=((2*s[0]*s[1]+c1)*(2*cv+c2))/((s[0]*s[0]+s[1]*s[1]+c1)*(vx+vy+c2));
            count++;
        }
        return count?total/count:1.0;
    }
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
    if(argc<3||argc>8){fprintf(stderr,"usage: %s original.obj simplified.obj [resolution=1024] [--dump-prefix path --dump-view 0..5]\n",argv[0]);return 1;}
    try{
    int res=argc>=4?meshio::parse_resolution_or_throw(argv[3]):1024;
    string dump_prefix;int dump_view=-1;
    for(int i=4;i<argc;){
        string option=argv[i++];
        if(i>=argc) throw runtime_error("missing value for "+option);
        if(option=="--dump-prefix") dump_prefix=argv[i++];
        else if(option=="--dump-view"){
            char* end=nullptr;long value=strtol(argv[i++],&end,10);
            if(*end!='\0'||value<0||value>5) throw runtime_error("dump view must be in 0..5");
            dump_view=(int)value;
        }
        else throw runtime_error("unknown option: "+option);
    }
    if((dump_prefix.empty())!=(dump_view<0)||dump_view>5) throw runtime_error("dump requires --dump-prefix and --dump-view 0..5");
    if(res!=1024) fprintf(stderr,"warning: non-1024 resolution is screening-only, not release evidence\n");
    Mesh a=meshio::read_mesh_or_throw(argv[1]),b=meshio::read_mesh_or_throw(argv[2]);
    if(getenv("NORMALIZE_ARM_RAW"))normalize_like_arm_input(a,b);
    double total=0, normal_total=0, depth_total=0;
    for(int v=0;v<6;v++){
        R ra=render(a,v,res),rb=render(b,v,res);
        vector<unsigned char>fg(res*res);
        for(int i=0;i<res*res;i++)fg[i]=ra.fg[i]||rb.fg[i];
        double ns=0;
        for(int ch=0;ch<3;ch++)ns+=ssim(ra,rb,fg,res,[ch](const R&r,int i){return normal_val(r.n[i],ch);});
        ns/=3;
        double ds=ssim(ra,rb,fg,res,[](const R&r,int i){return r.d[i];});
        if(v==dump_view) dump_diagnostics(ra,rb,res,dump_prefix);
        fprintf(stderr,"view=%d normal=%.12f depth=%.12f combined=%.12f\n",v,ns,ds,.5*ns+.5*ds);
        normal_total+=ns;
        depth_total+=ds;
        total+=.5*ns+.5*ds;
    }
    printf("normal=%.12f depth=%.12f combined=%.12f\n",normal_total/6,depth_total/6,total/6);
    }catch(const exception& e){fprintf(stderr,"error: %s\\n",e.what());return 2;}
}
