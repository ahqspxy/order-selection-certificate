// Candidate verifier for scale-normalized source-mass boxes.
// A universal result additionally needs complete coverage and a formal audit.
#include <algorithm>
#include <array>
#include <vector>
#include <map>
#include <cstdio>
#include <cmath>
#include "../higher_ratio/higher_interval.hpp"
using std::array;using std::vector;
static const int ends[5]={0,1,5,8,16};
static I alpha[4][25],beta0[4][25],zcoef[4][25];
static I GAMMA,C;
struct Row {I r,e,D;array<I,4> cp,capacity,prefix;};
struct Cell {Row l,u;I dt,I0,I1;int col,completed;};
struct Config {vector<Cell> cells;int n,prefix;I terminal;};
static std::map<int,Config> cache;
static Row row(int k,int den,int col){
 Row a;a.r=rational_interval(k,den);a.e=expI(a.r);a.D=C*a.e+I(1)+a.r;
 for(int g=0;g<4;g++){
  a.cp[g]=I(alpha[g][col].lo)*a.e+I(beta0[g][col].lo);
  I begin=rational_interval(2*ends[g],5),end=rational_interval(2*ends[g+1],5);
  a.capacity[g]=(5*k>=2*ends[g]*den)?I(1):expI(a.r-begin);
  a.prefix[g]=expI(a.r-end)-I(1);
 }
 return a;
}
static Config &config(int den){
 auto &q=cache[den];if(!q.cells.empty())return q;
 q.n=10*den;q.prefix=32*den/5;q.cells.reserve(q.n);
 const int sub=40;
 auto integrand=[&](int k){I r=rational_interval(k,sub*den),D=C*expI(r)+I(1)+r;return r/(D*D);};
 for(int k=0;k<q.n;k++){
  Cell a;a.col=5*k/(2*den);a.l=row(k,den,a.col);a.u=row(k+1,den,a.col);a.dt=rational_interval(1,den);
  a.completed=0;while(a.completed<4&&5*k>=2*ends[a.completed+1]*den)a.completed++;
  if(k<q.prefix){
   I sum=I(.5)*(integrand(sub*k)+integrand(sub*(k+1)));
   for(int j=1;j<sub;j++)sum=sum+integrand(sub*k+j);
   I dx=rational_interval(1,sub*den),tr=dx*sum;
   I M=(I(4)+I(8)*I(a.u.r.hi))/(I(a.l.D.lo)*I(a.l.D.lo));
   I err=M*a.dt*dx*dx/I(12);a.I0=I(std::max(0.,(tr-err).lo),(tr+err).hi);
   a.I1=((a.l.r+I(1))/a.l.D-(a.u.r+I(1))/a.u.D)/C;
  }
  q.cells.push_back(a);
 }
 q.terminal=expI(I(-20))*I(265)/(I(16)*C*C);
 return q;
}
struct Mass {
 double lo[4],hi[4],H,E,F,eps;
 Mass(const double *a,const double *b,double epsilon){
  I total(0);eps=epsilon;
  for(int g=0;g<4;g++){
   I len=rational_interval(2*(ends[g+1]-ends[g]),5);
   lo[g]=std::max(0.,(I(a[g])*len).lo);hi[g]=std::max(0.,(I(b[g])*len).hi);total=total+I(hi[g]);
  }
  H=(I(eps)*total).hi;E=expI(I(H)).hi;F=(I(1)+I(H)).hi;
 }
};
static double dlower(const Row &r,const Mass &m){
 I ans(0);for(int g=0;g<4;g++)ans=ans+I(std::min((r.cp[g]*I(m.lo[g])).lo,(r.cp[g]*I(m.hi[g])).lo));
 return ans.lo;
}
static double globalR(const Row &r,const Mass &m){
 double eta=0;for(int g=0;g<4;g++)eta=std::max(eta,-r.cp[g].lo);
 return ((I(1)+I(eta))*I(m.F)).hi;
}
static bool budget_ok(const Row &r,const Mass &m,int completed,double dl){
 I v(dl);for(int g=0;g<completed;g++)v=v+I(m.lo[g])*r.prefix[g];return v.lo>=0;
}
static double negativeR(const Row &r,const Mass &m,int completed,double dl){
 I sat(0);double sig=1.;
 for(int j=completed-1;j>=0;j--){
  I cap(0);for(int g=0;g<=j;g++)cap=cap+I(m.lo[g]);
  if(cap.lo>0){
   I root=I(1)+(I(-dl)-sat)/I(cap.lo);
   sig=std::max(sig,root.hi);
  }
  sat=sat+I(m.lo[j])*r.prefix[j];
 }
 return (I(sig)*I(m.F)).hi;
}
static double positiveR(const Row &r,const Mass &m,double dl){
 I sum(0),mass(0);double ratio=expI(r.r-rational_interval(32,5)).hi;
 I target=I(dl)/I(m.E);
 for(int g=0;g<4;g++){
  mass=mass+I(m.hi[g]);sum=sum+I(m.hi[g])*r.capacity[g];
  if(mass.lo<=0)continue;
  I root=(sum-target)/I(mass.lo);
  ratio=std::max(ratio,root.hi);
 }
 return (I(ratio)*I(m.E)).hi;
}
static double chord(const Cell &a,double left,double right){
 I gap=a.u.e-a.l.e;
 I W=(a.I1-a.l.e*a.I0)/gap;
 W=I(std::max(0.,W.lo),std::min(a.I0.hi,W.hi));
 return (I(left)*a.I0+(I(right)-I(left))*W).hi;
}
static double prefix_cell(const Cell &a,const Mass &m){
 double dl=dlower(a.l,m),dr=dlower(a.u,m);
 double best=chord(a,globalR(a.l,m),globalR(a.u,m));
 if(dl>=0&&dr>=0)best=std::min(best,chord(a,positiveR(a.l,m,dl),positiveR(a.u,m,dr)));
 else if(budget_ok(a.l,m,a.completed,dl)&&budget_ok(a.u,m,a.completed,dr))
  best=std::min(best,chord(a,negativeR(a.l,m,a.completed,dl),negativeR(a.u,m,a.completed,dr)));
 return best;
}
static double sqrt_upper(double x){
 audit_sqrt++;audit_require(std::isfinite(x)&&x>=0,"sqrt finite nonnegative input");
 if(x<=0)return 0;
 audit_require(x>=std::numeric_limits<double>::min(),"sqrt excludes subnormal positive input");
 audit_sqrt_min=std::min(audit_sqrt_min,x);audit_sqrt_max=std::max(audit_sqrt_max,x);
 double q=up(std::sqrt(x));unsigned long long steps=0;
 while((I(q)*I(q)).lo<x){q=up(q);steps++;audit_require(steps<=1024,"sqrt correction bounded by 1024 steps");}
 audit_sqrt_steps_max=std::max(audit_sqrt_steps_max,steps);return q;
}
static double quadratic(double A,double B,double z,double upper,double r,double DL,double DU,double V){
 double split=std::min(upper,std::max(0.,-z));I DD=I(DL)*I(DL);
 double ans=(I(r)*(I(A)+I(B)*I(split))/DD).hi;
 if(upper>split){
  I KV=DD*I(V)/I(DU);
  I bb=I(B)*(I(r)-I(z))-I(A)-KV;
  I cc=I(A)*(I(r)-I(z))-KV*I(z);
  // Upper quadratic -B*h^2 + bb.hi*h + cc.hi on h>=0.
  double q;
  if(bb.hi<(I(2)*I(B)*I(split)).lo){I h(split);q=(-I(B)*h*h+I(bb.hi)*h+I(cc.hi)).hi;}
  else if(bb.hi>(I(2)*I(B)*I(upper)).hi){I h(upper);q=(-I(B)*h*h+I(bb.hi)*h+I(cc.hi)).hi;}
  else q=(I(cc.hi)+I(bb.hi)*I(bb.hi)/(I(4)*I(B))).hi;
  ans=std::max(ans,(I(q)/DD).hi);
 }
 return ans;
}
static double tail(const Config &q,const Mass &m){
 double V=q.terminal.hi;
 for(int j=q.n-1;j>=q.prefix;j--){
  const Cell &a=q.cells[j];double dl=std::min(dlower(a.l,m),dlower(a.u,m));
  // d_lower is concave in exp(r), hence its minimum is at an endpoint.
  double z=0;for(int g=0;g<4;g++)z=dn(z+std::min((zcoef[g][a.col]*I(m.lo[g])).lo,(zcoef[g][a.col]*I(m.hi[g])).lo));
  double dd=(I(m.eps)*I(std::min(0.,dl))).lo,zz=(I(m.eps)*I(std::min(0.,z))).lo;
  double sigma=std::max(globalR(a.l,m),globalR(a.u,m));sigma=(I(sigma)/I(m.F)).hi;
  // A single budget comparison at the left endpoint is valid throughout the cell.
  if(budget_ok(a.l,m,a.completed,dl))sigma=std::min(sigma,(I(negativeR(a.l,m,a.completed,dl))/I(m.F)).hi);
  double hmax=(I(m.H)+a.u.r-rational_interval(32,5)).hi;
  double H=quadratic((I(1)-I(dd)).hi,1.,zz,hmax,a.u.r.hi,a.l.D.lo,a.u.D.hi,V);
  if(dd<0){double cap=std::min(hmax,sqrt_upper((-I(2)*I(dd)).hi));
   H=std::max(H,quadratic(sigma,sigma,zz,cap,a.u.r.hi,a.l.D.lo,a.u.D.hi,V));}
  I step=a.dt*I(H);
  if(H<0)step=step/(I(1)+a.dt*a.u.r/I(a.l.D.lo));
  V=(I(V)+step).hi;
 }
 return V;
}
extern "C" int initialize(const char *path){
 strict_mode=true;if(std::fesetround(FE_TONEAREST)!=0)return 1;
 FILE *f=std::fopen(path,"r");if(!f)return 2;
 for(int i=0;i<16;i++)for(int j=0;j<25;j++)if(fscanf(f,"%lf %lf",&wi[i][j].lo,&wi[i][j].hi)!=2)return 3;
 for(int i=0;i<16;i++)for(int j=0;j<26;j++)if(fscanf(f,"%lf %lf",&ci[i][j].lo,&ci[i][j].hi)!=2)return 4;
 fclose(f);
 for(int g=0;g<4;g++)for(int i=ends[g];i<ends[g+1];i++){
  for(int j=0;j<25;j++)if(wi[i][j].lo!=wi[ends[g]][j].lo||wi[i][j].hi!=wi[ends[g]][j].hi)return 5;
  for(int j=0;j<26;j++)if(ci[i][j].lo!=ci[ends[g]][j].lo||ci[i][j].hi!=ci[ends[g]][j].hi)return 6;
 }
 GAMMA=I(74500001)/I(100000000);C=(I(1)-GAMMA)/GAMMA;
 for(int g=0;g<4;g++)for(int j=0;j<25;j++){
  int i=ends[g];alpha[g][j]=(ci[i][j]+wi[i][j])*expI(-rational_interval(2*j,5));beta0[g][j]=-wi[i][j];zcoef[g][j]=wi[i][j];
 }
 cache.clear();return 0;
}
extern "C" double evaluate_normalized(const double *lo,const double *hi,int den,double epsilon,double *diagnostic){
 if(den<=0||den%5||epsilon<=0)return 1e30;
 const Config &q=config(den);Mass m(lo,hi,epsilon);I value(0);
 for(int j=0;j<q.prefix;j++)value=value+I(prefix_cell(q.cells[j],m));
 double t=tail(q,m);
 if(diagnostic){diagnostic[0]=value.hi;diagnostic[1]=t;diagnostic[2]=q.prefix;}
 return (GAMMA+value+I(t)).hi;
}
