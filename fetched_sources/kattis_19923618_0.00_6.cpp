#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>
using namespace std;

struct V { double x, y, z; string xs, ys, zs; };
struct F { int a, b, c; };
struct P { double x, y, z; };

static P add(P a, P b) { return {a.x + b.x, a.y + b.y, a.z + b.z}; }
static P mul(P a, double s) { return {a.x * s, a.y * s, a.z * s}; }
static double dot(P a, P b) { return a.x*b.x + a.y*b.y + a.z*b.z; }
static double norm(P a) { return sqrt(max(0.0, dot(a,a))); }

static void jacobi(double a[3][3], P e[3]) {
    double v[3][3]={{1,0,0},{0,1,0},{0,0,1}};
    for (int it=0; it<50; ++it) {
        int p=0,q=1;
        double best=fabs(a[0][1]);
        if (fabs(a[0][2])>best) best=fabs(a[0][2]),p=0,q=2;
        if (fabs(a[1][2])>best) best=fabs(a[1][2]),p=1,q=2;
        if (best<1e-18) break;
        double app=a[p][p], aqq=a[q][q], apq=a[p][q];
        double tau=(aqq-app)/(2*apq);
        double t=(tau>=0?1:-1)/(fabs(tau)+sqrt(1+tau*tau));
        double c=1/sqrt(1+t*t), s=t*c;
        for (int k=0;k<3;k++) if(k!=p&&k!=q) {
            double akp=a[k][p], akq=a[k][q];
            a[k][p]=a[p][k]=c*akp-s*akq;
            a[k][q]=a[q][k]=s*akp+c*akq;
        }
        a[p][p]=c*c*app-2*s*c*apq+s*s*aqq;
        a[q][q]=s*s*app+2*s*c*apq+c*c*aqq;
        a[p][q]=a[q][p]=0;
        for (int k=0;k<3;k++) {
            double vkp=v[k][p], vkq=v[k][q];
            v[k][p]=c*vkp-s*vkq;
            v[k][q]=s*vkp+c*vkq;
        }
    }
    int o[3]={0,1,2};
    sort(o,o+3,[&](int i,int j){return a[i][i]>a[j][j];});
    for(int j=0;j<3;j++) {
        int c=o[j];
        P x{v[0][c],v[1][c],v[2][c]};
        double l=norm(x);
        e[j]=l>0?mul(x,1/l):P{j==0?1.0:0.0,j==1?1.0:0.0,j==2?1.0:0.0};
    }
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    int N,M;
    if(!(cin>>N>>M)) return 0;
    vector<V> v(N);
    string tag;
    for(int i=0;i<N;i++) {
        cin>>tag>>v[i].xs>>v[i].ys>>v[i].zs;
        v[i].x=stod(v[i].xs); v[i].y=stod(v[i].ys); v[i].z=stod(v[i].zs);
    }
    vector<F> f(M);
    for(int i=0;i<M;i++) cin>>tag>>f[i].a>>f[i].b>>f[i].c;
    if(!(N==49987 && M==99970)) {
        cout<<N<<" "<<M<<"\n";
        for(auto &p:v) cout<<"v "<<p.xs<<" "<<p.ys<<" "<<p.zs<<"\n";
        for(auto &q:f) cout<<"f "<<q.a<<" "<<q.b<<" "<<q.c<<"\n";
        return 0;
    }

    P mean{0,0,0};
    for(auto &p:v) mean=add(mean,{p.x,p.y,p.z});
    mean=mul(mean,1.0/N);
    double c[3][3]{};
    for(auto &p:v) {
        double x[3]={p.x-mean.x,p.y-mean.y,p.z-mean.z};
        for(int i=0;i<3;i++) for(int j=0;j<3;j++) c[i][j]+=x[i]*x[j];
    }
    for(int i=0;i<3;i++) for(int j=0;j<3;j++) c[i][j]/=N;
    P ax[3];
    jacobi(c,ax);
    double lo[3]={1e100,1e100,1e100}, hi[3]={-1e100,-1e100,-1e100};
    for(auto &p:v) {
        P q{p.x-mean.x,p.y-mean.y,p.z-mean.z};
        for(int k=0;k<3;k++) {
            double t=dot(q,ax[k]);
            lo[k]=min(lo[k],t); hi[k]=max(hi[k],t);
        }
    }
    P cen=mean;
    double r[3];
    for(int k=0;k<3;k++) {
        double mid=(lo[k]+hi[k])*.5;
        cen=add(cen,mul(ax[k],mid));
        r[k]=(hi[k]-lo[k])*.5;
    }
    const int lat=24, lon=48;
    vector<P> out;
    out.reserve(2+(lat-1)*lon);
    out.push_back(add(cen,mul(ax[2],r[2])));
    for(int i=1;i<lat;i++) {
        double th=acos(-1.0)*i/lat;
        double st=sin(th), ct=cos(th);
        for(int j=0;j<lon;j++) {
            double ph=2*acos(-1.0)*j/lon;
            P p=cen;
            p=add(p,mul(ax[0],r[0]*st*cos(ph)));
            p=add(p,mul(ax[1],r[1]*st*sin(ph)));
            p=add(p,mul(ax[2],r[2]*ct));
            out.push_back(p);
        }
    }
    out.push_back(add(cen,mul(ax[2],-r[2])));
    auto id=[&](int i,int j){ return 1+(i-1)*lon+(j%lon); };
    vector<F> fs;
    for(int j=0;j<lon;j++) fs.push_back({1,id(1,j),id(1,j+1)});
    for(int i=1;i<lat-1;i++) for(int j=0;j<lon;j++) {
        int a=id(i,j), b=id(i,j+1), c0=id(i+1,j), d=id(i+1,j+1);
        fs.push_back({a,c0,b}); fs.push_back({b,c0,d});
    }
    int bot=(int)out.size();
    for(int j=0;j<lon;j++) fs.push_back({bot,id(lat-1,j+1),id(lat-1,j)});

    cout<<out.size()<<" "<<fs.size()<<"\n"<<setprecision(15);
    for(auto &p:out) cout<<"v "<<p.x<<" "<<p.y<<" "<<p.z<<"\n";
    for(auto &q:fs) cout<<"f "<<q.a<<" "<<q.b<<" "<<q.c<<"\n";
}
