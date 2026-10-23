// Strict four-block specialization. A box bound becomes a global theorem
// only with complete coverage, mathematical audit and independent replay.
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <functional>
#include <map>
#include <tuple>
#include "../higher_ratio/higher_interval.hpp"
#include <vector>
using std::vector; using std::array;
namespace source_tail {
#include "normalized_finite.cpp"
}
// The certified exponential is deterministic for a fixed interval input.
// Clear on every evaluation so memory use cannot grow across mass boxes.
static std::map<std::pair<double,double>,I> exponential_memo;
static I memo_expI(I a){
 if(!strict_mode)return expI(a);
 auto key=std::make_pair(a.lo,a.hi);auto it=exponential_memo.find(key);
 if(it!=exponential_memo.end())return it->second;
 I result=expI(a);exponential_memo.emplace(key,result);return result;
}
#define expI memo_expI

static double w[16][25],cc[16][26];
static int kernel_generation=0;
static double gain_beta=0;
static constexpr int NG=GROUP_COUNT;
static const int gs[NG+1]={0,1,5,8,16};
static const double ts[NG+1]={0,0.4,2,3.2,6.4};
static int phase_at(double r){int p=0;while(p<NG&&r>=ts[p+1])p++;return p;}
static double clip(double x,double a,double b){return std::max(a,std::min(b,x));}
struct Knapsack {
 int n; array<double,16> cost,cap; double negcap;
 void init(const double *c,const double *q,int a,int b){
  n=b-a;array<int,16> order;for(int i=0;i<n;i++)order[i]=a+i;
  std::sort(order.begin(),order.begin()+n,[&](int x,int y){return c[x]<c[y];});
  negcap=0;for(int i=0;i<n;i++){cost[i]=c[order[i]];cap[i]=q[order[i]];if(cost[i]<0)negcap+=cap[i];}
 }
 double exact(double mass)const{
  double v=0;for(int i=0;i<n&&mass>0;i++){double q=std::min(mass,cap[i]);v+=q*cost[i];mass-=q;}return v;
 }
 double lower(double lo,double hi)const{
  if(!strict_mode)return exact(clip(negcap,lo,hi));
  double rl=lo,ru=hi,v=0;
  for(int i=0;i<n;i++){
   double amount=cost[i]<0?std::min(ru,cap[i]):std::min(rl,cap[i]);
   if(cost[i]!=0&&amount>0)v=dn(v+dn(cost[i]*amount));
   rl=rl==0?0:std::max(0.,dn(rl-cap[i]));ru=ru==0?0:std::max(0.,up(ru-cap[i]));
  }return v;
 }
};
struct Row {
 double r,D,scale,blockc[NG]; int phase;
 Knapsack dc[NG][2],zc[NG][2];
};
static void coefficients(double r,int j,double *c,double *z){
 for(int i=0;i<16;i++){
  if(j>=25){c[i]=z[i]=0;continue;}
  z[i]=w[i][j];c[i]=cc[i][j]+std::expm1(r-.4*j)*(cc[i][j]+w[i][j]);
 }
}
static Row pointrow(double r,int j,double gamma){
 Row out;out.r=r;out.D=(1-gamma)/gamma*std::exp(r)+1+r;
 out.phase=phase_at(r);
 double c[16],z[16],cp[16],cf[16];coefficients(r,j,c,z);
 double eta=0;for(int i=0;i<16;i++){eta=std::max(eta,-c[i]);cp[i]=clip(r-.4*i,0,.4);cf[i]=.4-cp[i];}
 for(int g=0;g<NG;g++)out.blockc[g]=c[gs[g]];
 out.scale=std::min(std::exp(r),1+eta);
 for(int g=0;g<NG;g++){
  out.dc[g][0].init(c,cp,gs[g],gs[g+1]);out.dc[g][1].init(c,cf,gs[g],gs[g+1]);
  out.zc[g][0].init(z,cp,gs[g],gs[g+1]);out.zc[g][1].init(z,cf,gs[g],gs[g+1]);
 }return out;
}
struct PrefixCost {
 Knapsack kn;I base,sumlo;
 double lower(double massL,double massU)const{
  double lo=std::max(0.,(I(massL)-sumlo).lo);
  double hi=std::max(lo,(I(massU)-sumlo).hi);
  return (base+I(kn.lower(lo,hi))).lo;
 }
};
struct Cell {
 double rL,rU,DL,DU,scale,td,tz,nz;int phase;
 PrefixCost prefixcost[2][3];
 I co[NG],zo[NG],co_mid[NG],dev[NG],dcenter,zcenter,dcenter_mid,ddeviation,zdeviation,quadratic_remainder;
 I rmid,Dmid,drange;
 Knapsack dc[NG][2],zc[NG][2],totald[NG],totalz[NG],negz[NG];
};
static Cell intervalrow(double rL,double rU,int j,double gamma){
 Cell out;I L(rL),R(rU);
 if(strict_mode){L=I(dn(rL),up(rL));R=I(dn(rU),up(rU));if(rL==0)L=I(0);}
 out.rL=L.lo;out.rU=R.hi;
 I ga=gamma_interval(gamma),c=(I(1)-ga)/ga;
 out.DL=(c*expI(L)+I(1)+L).lo;out.DU=(c*expI(R)+I(1)+R).hi;
 double rm=(rL+rU)/2;out.phase=phase_at(rm);
 out.rmid=I(rm);out.Dmid=c*expI(out.rmid)+I(1)+out.rmid;out.drange=I(out.rL,out.rU)-out.rmid;
 double cL[16],cU[16],z[16],dummy[16],cp[16],cf[16];
 if(!strict_mode){coefficients(rL,j,cL,z);coefficients(rU,j,cU,dummy);}
 else {
  I left=rational_interval(2*j,5),el=expI(L-left)-I(1),er=expI(R-left)-I(1);
  for(int i=0;i<16;i++){
   if(j>=25){cL[i]=cU[i]=z[i]=0;continue;}
   cL[i]=(ci[i][j]+el*(ci[i][j]+wi[i][j])).lo;
   cU[i]=(ci[i][j]+er*(ci[i][j]+wi[i][j])).lo;z[i]=wi[i][j].lo;
  }
 }
 double eta=0;for(int i=0;i<16;i++){
  cL[i]=std::min(cL[i],cU[i]);eta=std::max(eta,-cL[i]);
  cp[i]=clip((R-rational_interval(2*i,5)).hi,0,rational_interval(2,5).hi);
  cf[i]=clip((rational_interval(2*(i+1),5)-L).hi,0,rational_interval(2,5).hi);
 }
 out.scale=std::min(expI(R).hi,up(1+eta));
 for(int g=0;g<NG;g++){
  out.dc[g][0].init(cL,cp,gs[g],gs[g+1]);out.dc[g][1].init(cL,cf,gs[g],gs[g+1]);
  out.zc[g][0].init(z,cp,gs[g],gs[g+1]);out.zc[g][1].init(z,cf,gs[g],gs[g+1]);
 }
 double cap[16],minus[16];for(int i=0;i<16;i++){cap[i]=rational_interval(2,5).hi;minus[i]=j>=25?0:-(strict_mode?wi[i][j].hi:z[i]);}
 for(int g=0;g<NG;g++){
  out.totald[g].init(cL,cap,gs[g],gs[g+1]);out.totalz[g].init(z,cap,gs[g],gs[g+1]);out.negz[g].init(minus,cap,gs[g],gs[g+1]);
 }
 for(int g=0;g<NG;g++){
  if(j>=25){out.co[g]=out.zo[g]=I(0);continue;}
  int i=gs[g];I left=rational_interval(2*j,5);
  I a=ci[i][j]+(expI(L-left)-I(1))*(ci[i][j]+wi[i][j]);
  I b=ci[i][j]+(expI(R-left)-I(1))*(ci[i][j]+wi[i][j]);
  out.co[g]=I(std::min(a.lo,b.lo),std::max(a.hi,b.hi));out.zo[g]=wi[i][j];
  out.co_mid[g]=ci[i][j]+(expI(out.rmid-left)-I(1))*(ci[i][j]+wi[i][j]);
 }
 return out;
}
struct Mass {
 double lo[NG],hi[NG],mid[NG],HL[NG+1],HU[NG+1],HM[NG+1];
 Mass(const double *l,const double *u){HL[0]=HU[0]=HM[0]=0;for(int g=0;g<NG;g++){
  I length=rational_interval(2*(gs[g+1]-gs[g]),5);
  lo[g]=(I(l[g])*length).lo;hi[g]=(I(u[g])*length).hi;mid[g]=(lo[g]+hi[g])/2;
  lo[g]=std::max(0.,lo[g]);hi[g]=std::max(0.,hi[g]);
  HL[g+1]=HL[g]==0&&lo[g]==0?0:dn(HL[g]+lo[g]);HU[g+1]=HU[g]==0&&hi[g]==0?0:up(HU[g]+hi[g]);HM[g+1]=HM[g]+mid[g];
  HL[g+1]=std::max(0.,HL[g+1]);}}
 double low(double r)const{int p=phase_at(r);return p==NG?HL[NG]:std::max(HL[p],(I(HL[p+1])-rational_interval(2*gs[p+1],5)+I(r)).lo);}
 double high(double r)const{int p=phase_at(r);return p==NG?(I(HU[NG])+I(r)-rational_interval(32,5)).hi:std::min(HU[p+1],(I(HU[p])+I(r)-rational_interval(2*gs[p],5)).hi);}
 double anchor(double r)const{int p=phase_at(r);return p==NG?HM[NG]:std::max(HM[p],HM[p+1]-ts[p+1]+r);}
};
static double surplus(double h){return h<.01?h*h*(.5+h*(1./6+h*(1./24+h*(1./120+h/720)))):std::expm1(h)-h;}

static double gain_c=(1-.74501)/.74501,negative_z=0.;
static I gain_ci;
static const Mass *gain_mass=nullptr;
static I beta_interval(double rL,double rU,double hL,double len,double DU){
 I lo=I(rL)-I(len);double b=.2545;
 if(lo.lo>=1.7){I dd=gain_ci*expI(I(lo.lo))+I(1)+I(lo.lo);b=std::max(b,(I(1)-I(2)*I(lo.lo)/dd).lo);}
 else if(rU<=1.5){I dd=gain_ci*expI(I(rU))+I(1)+I(rU);b=std::max(b,(I(1)-I(2)*I(rU)/dd).lo);}
 I mass=I(std::max((I(hL)-I(len)/I(2)).lo,gain_mass->low(std::max(0.,lo.lo))))-I(negative_z);
 return I(b)+I(std::max(0.,mass.lo))/I(DU);
}
static double positive_series(double a,double l,double beta){
 double sum=0,term=l;
 for(int j=0;j<=6;j++){sum+=term*(a/(j+1.)-l/(j+2.));term*=beta*l/(j+1.);}
 return sum;
}
static I positive_series_I(double a,double l,I beta){
 I sum(0),term(l);
 for(int j=0;j<=6;j++){sum=sum+term*(I(a)/I(j+1)-I(l)/I(j+2));term=term*beta*I(l)/I(j+1);}
 return sum;
}
static double gain(double r,double h,double d,double zeta,double D){
 double delta=surplus(h)+d,hb=h+zeta,av=std::max(0.,r-hb);
 if(delta<=0||hb<=0||av<=0)return 0;
 double len=std::min(std::log1p(delta/hb),av),b=.2545,s=r-len;
 if(s>=1.7)b=std::max(b,1-2*s/(gain_c*std::exp(s)+1+s));
 else if(r<=1.5)b=std::max(b,1-2*r/D);
 b+=std::max(0.,std::max(h-len/2,gain_mass->low(std::max(0.,s)))-negative_z)/D;
 return positive_series(av,len,b)/(D*D);
}
static void zeta_levels(const double *M,double *out,bool certified=false){
 double cap[16],cost[16];for(int i=0;i<16;i++)cap[i]=certified?rational_interval(2,5).hi:.4;
 double running=0;for(int j=0;j<25;j++){
  for(int i=0;i<16;i++)cost[i]=-(certified?wi[i][j].hi:w[i][j]);
  double val=0;for(int g=0;g<NG;g++){Knapsack kn;kn.init(cost,cap,gs[g],gs[g+1]);val=(M==gain_mass->hi)?up(val-kn.lower(gain_mass->lo[g],gain_mass->hi[g])):val-kn.exact(M[g]);}
  out[j]=running=std::max(val,running);
 }out[25]=out[24];
}
static double zeta(double r,int j,const double *levels){if(j>=25)return levels[25];double f=clip((r-.4*j)/.4,0,1);return levels[j]*(1-f)+levels[j+1]*f;}
struct Potential {int first; vector<double> p; double at(int j)const{return p.at(j-first);} };
static Potential potential(const Row &row,const Mass &m,double zz,int hden,int tden){
 double dh=1./hden,dt=1./tden,b=m.anchor(row.r);
 int first=std::max(0,int(std::floor(std::min(m.low(std::max(0.,row.r-dt)),b)*hden))-2);
 int last=int(std::ceil(std::max(m.high(row.r+dt),b)*hden))+2;
 Potential out;out.first=first;out.p.resize(last-first+1);double accum=0,old=0;
 double d=0;for(int g=0;g<NG;g++)d+=row.blockc[g]*m.mid[g];
 for(int idx=first;idx<=last;idx++){
  double h=double(idx)/hden;
  double val=gain(row.r,h,d,zz,row.D);
  if(idx>first)accum+=dh*.5*(old+val);
  out.p[idx-first]=accum;old=val;
 }
 double bf=b*hden;int bi=int(std::floor(bf));double f=bf-bi;
 double offset=out.at(bi)*(1-f)+out.at(bi+1)*f;
 for(double &x:out.p)x-=offset;for(size_t q=1;q<out.p.size();q++)out.p[q]=std::max(out.p[q],out.p[q-1]);return out;
}
static void mass_ranges(int g,const Cell &c,const Mass &m,double hl,double hu,double &pl,double &pu,double &fl,double &fu){
 if(g<c.phase){pl=m.lo[g];pu=m.hi[g];fl=fu=0;return;}
 if(g>c.phase){pl=pu=0;fl=m.lo[g];fu=m.hi[g];return;}
 pl=clip(dn(hl-m.HU[g]),0,m.lo[g]);pu=clip(up(hu-m.HL[g]),0,m.hi[g]);
 fl=std::max(0.,std::max(dn(m.lo[g]-pu),dn(m.HL[g+1]-hu)));
 fu=std::max(0.,std::min(up(m.hi[g]-pl),up(m.HU[g+1]-hl)));
}
static int exp_hden=0;
static vector<I> *exp_hgrid=nullptr;
static void prepare_hgrid(int hden,int end){
 static std::map<std::pair<int,int>,vector<I>> cache;
 auto key=std::make_pair(hden,end);auto &v=cache[key];
 if(v.empty()){v.reserve((end+2)*hden+1);for(int i=0;i<=(end+2)*hden;i++)v.push_back(expI(rational_interval(i,hden)));}
 exp_hden=hden;exp_hgrid=&v;
}
static double surplus_lower(double h){
 if(!strict_mode)return surplus(h);
 int j=std::max(0,int(std::floor(h*exp_hden)));
 I x=rational_interval(j,exp_hden);
 if(j>0&&x.hi>h){j--;x=rational_interval(j,exp_hden);}
 I ex=exp_hgrid->at(j);double d=std::max(0.,(I(h)-x).lo);
 return ((ex-I(1)-x)+I(std::max(0.,(ex-I(1)).lo))*I(d)).lo;
}
static double residual(const Cell &c,const Mass &m,double hl,double hu,double W,double Tlow,double PU,double PhL,double PrL,double zz){
 double dl=c.td,nz=c.nz,zl=c.tz;
 double gapL=std::max(0.,(I(c.rU)-rational_interval(2*gs[c.phase],5)).hi);
 double gapR=std::max(0.,(rational_interval(2*gs[c.phase+1],5)-I(c.rL)).hi);
 for(int side=0;side<2;side++){
  double low=side?hl:std::max(0.,(I(hl)-I(gapL)).lo);
  double high=side?(I(hu)+I(gapR)).hi:hu;
  dl=std::max(dl,c.prefixcost[side][0].lower(low,high));
  zl=std::max(zl,c.prefixcost[side][1].lower(low,high));
  nz=std::max(nz,c.prefixcost[side][2].lower(low,high));
 }
 double kl=clip(dn(hl+zl),0,c.rU);
 double ku=clip(up(hu-nz),kl,c.rU);
 double U=(I(1)+I(hu)-I(dl)).hi;
 double E=surplus_lower(hl);
 if(dn(E+dl)<0){
  double basic=U;
  double rho=(I(c.scale)*(I(1)+I(hu))).hi;
  // If tau < r, kernel comparison gives tau >= rho_time, hence
  // exp(h(tau))-1 >= exp(h(r))/scale - 1.
  // The area deficit -Delta therefore bounds the negative inverse slack.
  double denominator=((I(E)+I(1)+I(hl))/I(c.scale)-I(1)).lo;
  if(denominator>0){
   double deficit=std::max(0.,-(I(E)+I(dl)).lo);
   double continuous=(I(basic)+I(deficit)/I(denominator)).hi;
   rho=std::min(rho,continuous);
  }

  // If tau < r, B0(r)-B(r) >= e^(r-tau)(exp(h(tau))-1)-h(tau).
  // Thus R <= (h(tau)-d)/(1-exp(-h(tau))). This is decreasing in
  // h(tau) on the negative-slack branch. The mass cones give a lower bound.
  double shift=std::max(0.,(I(c.scale)-I(1)).hi);
  double timeLow=std::max(0.,(I(c.rL)-I(shift)).lo);
  double massLow=std::max(0.,std::max(m.low(timeLow),(I(hl)-I(shift)).lo));
  if(massLow>0){
   double zz=(I(massLow)+I(std::max(0.,surplus_lower(massLow)))).lo;
   if(zz>0){
    double cap=((I(massLow)-I(dl))*(I(1)+I(zz))/I(zz)).hi;
    rho=std::min(rho,cap);
   }
  }
  U=std::max(basic,rho);

 }
 U=std::max(U,0.);
 // Monotonicity gives G_L from r_L, h_U+zeta_U, and Delta_L.
 double dlt=dn(E+dl),hb=up(hu+zz),av=std::max(0.,dn(c.rL-hb)),GL=0;
 if(dlt>0&&hb>0&&av>0){double el=std::min(log1p_lower(dn(dlt/hb)),av);GL=(positive_series_I(av,el,beta_interval(c.rL,c.rU,hl,el,c.DU))/(I(c.DU)*I(c.DU))).lo;GL=std::max(0.,GL);}
 double V=(I(Tlow)+I(W)-I(PU)).lo,den=V>=0?c.DU:c.DL;
 I rp=I(c.rL,c.rU)*(I(U)-I(1));
 auto f=[&](double k){double N=(rp-I(k)*I(U)).hi;I den2=N>=0?I(c.DL)*I(c.DL):I(c.DU)*I(c.DU);return (I(N)/den2-I(k)*I(V)/I(den)).hi;};
 double H=std::max(f(kl),f(ku));
 if(dlt>=0){
  // Keep the same mass vector in the two linear corrections. The only
  // discarded mass correlation is a quadratic, explicitly enclosed term.
  I h((hl+hu)/2),dh=I(hl,hu)-h,r(c.rL,c.rU),D(c.DL,c.DU),D2=D*D;
  I d0=c.dcenter,z0=c.zcenter,U0=I(1)+h-d0,k0=h+z0;
  I N0=r*(h-d0)-k0*U0;
  I base=N0/D2-k0*I(V)/D;
  // At fixed state and center masses, d0'=d0+z0 and D'=D-r.
  // The mean-value bound keeps time dependence in this common base term.
  I Dprime=D-r;
  I derivative=(h-d0+(k0-r)*(d0+z0)-I(2)*N0*Dprime/D+k0*I(V)*Dprime)/D2;
  I dm=c.dcenter_mid,Dm=c.Dmid;
  I midpoint=(c.rmid*(h-dm)-k0*(I(1)+h-dm))/(Dm*Dm)-k0*I(V)/Dm;
  double timed=(midpoint+derivative*c.drange).hi;
  base=I(std::min(base.hi,timed));
  for(int g=0;g<NG;g++){
   I linear=(-(r-k0)*c.co[g]-U0*c.zo[g])/D2-c.zo[g]*I(V)/D;
   base=base+I((linear*c.dev[g]).hi);
  }
  I hlinear=(r-I(1)-I(2)*h+d0-z0)/D2-I(V)/D;
  base=base+I((hlinear*dh).hi);
  // The remaining pure-state quadratic is -dh^2 and can be discarded.
  base=base+(c.quadratic_remainder+dh*(c.ddeviation-c.zdeviation))/D2;
  H=std::min(H,base.hi);
 }
 return (I(H)-I(std::min(GL,PhL))-I(PrL)).hi;
}
struct Config {vector<Row> rows;vector<Cell> cells;vector<double> T;double T0upper;};
static Config &configuration(double gamma,int tden,int end){
 static std::map<std::tuple<int,int,bool,double>,Config> cache;static int saved_generation=-1;
if(saved_generation!=kernel_generation){cache.clear();saved_generation=kernel_generation;}
 auto key=std::make_tuple(tden,end,strict_mode,gamma);Config &out=cache[key];
 if(!out.T.empty())return out;
 int n=end*tden;double dt=1./tden;
 out.rows.reserve(n+1);out.cells.reserve(n);
 for(int j=0;j<=n;j++){int col=j/(tden*2/5);out.rows.push_back(pointrow(double(j)/tden,col,gamma));
  if(j<n)out.cells.push_back(intervalrow(double(j)/tden,double(j+1)/tden,col,gamma));}
 out.T.assign(n+1,0.);double upper=0;
 I ga=gamma_interval(gamma),c=(I(1)-ga)/ga,dx=rational_interval(1,20*tden);
 auto fi=[&](int k){I r=rational_interval(k,20*tden),D=c*expI(r)+I(1)+r;return r/(D*D);};
 for(int j=n-1;j>=0;j--){
  I sum=I(.5)*(fi(20*j)+fi(20*(j+1)));
  for(int k=1;k<20;k++)sum=sum+fi(20*j+k);
  I integral=dx*sum;
  const Cell &cell=out.cells[j];
  I M=(I(4)+I(8)*I(cell.rU))/(I(cell.DL)*I(cell.DL));
  I error=M*rational_interval(1,tden)*dx*dx/I(12);
  if(!strict_mode)error=I(0);
  out.T[j]=std::max(0.,(I(out.T[j+1])+integral-error).lo);
  upper=(I(upper)+integral+error).hi;
 }
 out.T0upper=upper;return out;
}
extern "C" void set_gain_beta(double b){gain_beta=b;}
extern "C" int init_kernel(const char *path){kernel_generation++;FILE *f=std::fopen(path,"r");if(!f)return 1;
 for(int i=0;i<16;i++)for(int j=0;j<25;j++)if(fscanf(f,"%lf",&w[i][j])!=1)return 2;
 for(int i=0;i<16;i++)for(int j=0;j<26;j++)if(fscanf(f,"%lf",&cc[i][j])!=1)return 3;
 fclose(f);return 0;}
extern "C" int init_strict(const char *path){
 static_assert(std::numeric_limits<double>::is_iec559&&std::numeric_limits<double>::digits==53);
 if(std::fesetround(FE_TONEAREST)!=0)return 4;
 FILE *f=std::fopen(path,"r");if(!f)return 1;
 for(int i=0;i<16;i++)for(int j=0;j<25;j++)if(fscanf(f,"%lf %lf",&wi[i][j].lo,&wi[i][j].hi)!=2)return 2;
 for(int i=0;i<16;i++)for(int j=0;j<26;j++)if(fscanf(f,"%lf %lf",&ci[i][j].lo,&ci[i][j].hi)!=2)return 3;
 fclose(f);
 for(int g=0;g<NG;g++)for(int i=gs[g];i<gs[g+1];i++){
  for(int j=0;j<25;j++)if(wi[i][j].lo!=wi[gs[g]][j].lo||wi[i][j].hi!=wi[gs[g]][j].hi)return 5;
  for(int j=0;j<26;j++)if(ci[i][j].lo!=ci[gs[g]][j].lo||ci[i][j].hi!=ci[gs[g]][j].hi)return 6;
 }
 strict_mode=true;return source_tail::initialize(path);
}
extern "C" double evaluate(const double *lo,const double *hi,double gamma,int tden,int hden,int end,double *diagnostic){

 if(gamma!=.74500001)return 1e30;
 exponential_memo.clear();
 I gg=gamma_interval(gamma);gain_ci=(I(1)-gg)/gg;gain_c=(1-gamma)/gamma;
 Mass m(lo,hi);gain_mass=&m;negative_z=0.;
 double caps[16],cost[16];for(int i=0;i<16;i++)caps[i]=rational_interval(2,5).hi;
 for(int j=0;j<25;j++){
  for(int i=0;i<16;i++)cost[i]=strict_mode?wi[i][j].lo:w[i][j];
  I val(0);for(int g=0;g<NG;g++){Knapsack kn;kn.init(cost,caps,gs[g],gs[g+1]);val=val+I(kn.lower(m.lo[g],m.hi[g]));}
  negative_z=std::max(negative_z,-val.lo);
 }
double zm[26],zu[26];zeta_levels(m.mid,zm);zeta_levels(m.hi,zu,strict_mode);
 int n=32*tden/5;double dt=1./tden;
 Config &cfg=configuration(gamma,tden,end);
 const vector<Row> &rows=cfg.rows;vector<Cell> &cells=cfg.cells;const vector<double> &T=cfg.T;
 for(Cell &cell:cells){
  cell.td=cell.tz=cell.nz=0;
  for(int g=0;g<NG;g++){
   cell.td=dn(cell.td+cell.totald[g].lower(m.lo[g],m.hi[g]));
   cell.tz=dn(cell.tz+cell.totalz[g].lower(m.lo[g],m.hi[g]));
   cell.nz=dn(cell.nz+cell.negz[g].lower(m.lo[g],m.hi[g]));
  }
  cell.dcenter=cell.zcenter=cell.dcenter_mid=I(0);I dd(0),dz(0);
  for(int g=0;g<NG;g++){
   cell.dev[g]=I(m.lo[g],m.hi[g])-I(m.mid[g]);
   cell.dcenter=cell.dcenter+cell.co[g]*I(m.mid[g]);cell.zcenter=cell.zcenter+cell.zo[g]*I(m.mid[g]);
   cell.dcenter_mid=cell.dcenter_mid+cell.co_mid[g]*I(m.mid[g]);
   dd=dd+cell.co[g]*cell.dev[g];dz=dz+cell.zo[g]*cell.dev[g];
  }
  cell.ddeviation=dd;cell.zdeviation=dz;cell.quadratic_remainder=dd*dz;
  if(cell.phase<NG)for(int side=0;side<2;side++)for(int kind=0;kind<3;kind++){
   int cutoff=cell.phase+side;PrefixCost &pc=cell.prefixcost[side][kind];
   pc.base=I(0);pc.sumlo=I(0);double costs[16],capacity[16];
   for(int g=0;g<NG;g++){
    double cost=kind==0?cell.totald[g].cost[0]:(kind==1?cell.totalz[g].cost[0]:cell.negz[g].cost[0]);
    if(g<cutoff){
     costs[g]=cost;capacity[g]=std::max(0.,(I(m.hi[g])-I(m.lo[g])).hi);
     pc.base=pc.base+I(cost)*I(m.lo[g]);pc.sumlo=pc.sumlo+I(m.lo[g]);
    }else pc.base=pc.base+I(std::min((I(cost)*I(m.lo[g])).lo,(I(cost)*I(m.hi[g])).lo));
   }
   pc.kn.init(costs,capacity,0,cutoff);
  }
 }

 if(strict_mode)prepare_hgrid(hden,end);
 Potential next=potential(rows[n],m,zeta(6.4,16,zm),hden,tden);
 I ga=gamma_interval(gamma),ccert=(I(1)-ga)/ga;
 source_tail::FiniteMass tailmass(lo,hi,1.,1.);
 double tailupper=source_tail::finite_tail(source_tail::config(tden),tailmass);
 double W=(I(tailupper)-I(T[n])).hi;
 double terminalP=0;for(double p:next.p)terminalP=std::max(terminalP,p);W=up(W+terminalP);
 double maxham=-1e100;long long boxes=0;int previous_best=0;
 for(int ix=n-1;ix>=0;ix--){
  const Cell &cell=cells[ix];int col=ix/(tden*2/5);
  Potential prev=potential(rows[ix],m,zeta(cell.rL,col,zm),hden,tden);
  double zU=zeta(cell.rU,col,zu);if(strict_mode)zU=zu[std::min(col,25)];
  int a=std::max(0,int(std::floor(dn(m.low(cell.rL)*hden))));
  int b=std::max(a+1,int(std::ceil(up(m.high(cell.rU)*hden))));double Hmax=-1e100;
  vector<double> ph(4*(b-a)+4),pr(ph.size());
  std::function<void(int,int,int)> build=[&](int node,int l,int u){
   if(u-l==1){
    double p00=prev.at(l),p01=prev.at(u),p10=next.at(l),p11=next.at(u);
    ph[node]=dn(std::min(dn(p01-p00),dn(p11-p10))*hden);
    pr[node]=dn(std::min(dn(p10-p00),dn(p11-p01))*tden);
   }else{int mid=(l+u)/2;build(2*node,l,mid);build(2*node+1,mid,u);ph[node]=std::min(ph[2*node],ph[2*node+1]);pr[node]=std::min(pr[2*node],pr[2*node+1]);}
  };build(1,a,b);
  auto bound=[&](int node,int l,int u){
   double hl=std::max(rational_interval(l,hden).lo,m.low(cell.rL));
   double hu=std::min(rational_interval(u,hden).hi,m.high(cell.rU));
   double PU=std::max(prev.at(u),next.at(u));boxes++;
   return residual(cell,m,hl,hu,W,T[ix+1],PU,ph[node],pr[node],zU);
  };
  auto at_leaf=[&](int j){
   int node=1,l=a,u=b;while(u-l>1){int mid=(l+u)/2;if(j<mid){node*=2;u=mid;}else{node=2*node+1;l=mid;}}
   double val=bound(node,l,u);if(val>Hmax){Hmax=val;previous_best=j;}
  };
  at_leaf(int(clip(previous_best,a,b-1)));at_leaf((a+b-1)/2);at_leaf(a);at_leaf(b-1);
  const double tolerance=2e-7;
  std::function<void(int,int,int)> search=[&](int node,int l,int u){
   double val=bound(node,l,u);if(val<=up(Hmax+tolerance))return;
   if(u-l==1){Hmax=val;previous_best=l;return;}
   int mid=(l+u)/2;search(2*node,l,mid);search(2*node+1,mid,u);
  };search(1,a,b);Hmax=up(Hmax+tolerance);
  maxham=std::max(maxham,Hmax);
  I idt=rational_interval(1,tden);I increment=idt*I(Hmax);
  if(Hmax<0)increment=increment/(I(1)+idt*(I(cell.rU)/I(cell.DL)));
  W=up(W+increment.hi);next=std::move(prev);
 }
 if(diagnostic){diagnostic[0]=double(boxes);diagnostic[1]=maxham;}
 return (gamma_interval(gamma)+I(cfg.T0upper)+I(W)).hi;
}
