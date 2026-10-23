// A finite-scale refinement of the independently certified boundary bound.
// Research verifier: completion and formal replay still required.
#include "normalized_cert.cpp"
struct FiniteMass:Mass {
 double scaleL,HL[5],HU[5],EW[4],EH[5],surplus[5];
 FiniteMass(const double*a,const double*b,double sl,double su):Mass(a,b,su){
  scaleL=sl;HL[0]=HU[0]=0;EH[0]=1;surplus[0]=0;
  I cl(0),cu(0);
  for(int g=0;g<4;g++){
   I el=expI(-I(sl)*cl);
   if(sl>0)EW[g]=(el*(I(1)-expI(-I(sl)*I(hi[g])))/I(sl)).hi;
   else EW[g]=hi[g];
   EW[g]=std::max(0.,EW[g]);
   cl=cl+I(lo[g]);cu=cu+I(hi[g]);
   HL[g+1]=std::max(0.,(I(sl)*cl).lo);HU[g+1]=(I(su)*cu).hi;
   EH[g+1]=expI(I(HU[g+1])).hi;
   I h(HL[g+1]);surplus[g+1]=std::max(0.,(h*h/I(2)+h*h*h/I(6)).lo);
  }
 }
};
static double positive_finite(const Row&r,const FiniteMass&m,int phase,double dl){
 I target=I(dl)+I(m.surplus[phase])/I(m.eps);
 I u[4];for(int g=0;g<4;g++)u[g]=g<=phase?I(m.EH[phase+1]):r.capacity[g]*I(m.EH[g]);
 double ans=(expI(r.r-rational_interval(32,5))*I(m.EH[4])).hi;
 // Max over all subsets avoids any assumption on the order of enlarged knots.
 for(int mask=1;mask<16;mask++){
  I mass(0),sum(0);for(int g=0;g<4;g++)if(mask&(1<<g)){mass=mass+I(m.EW[g]);sum=sum+I(m.EW[g])*u[g];}
  if(mass.lo>0)ans=std::max(ans,((sum-target)/I(mass.lo)).hi);
 }
 return ans;
}
static double finite_prefix_cell(const Cell&a,const FiniteMass&m){
 int phase=a.completed;
 double dl=dlower(a.l,m),dr=dlower(a.u,m);
 double factor=((I(1)+I(m.HU[phase+1]))/I(m.F)).hi;
 double best=chord(a,(I(globalR(a.l,m))*I(factor)).hi,(I(globalR(a.u,m))*I(factor)).hi);
 I gain=I(m.surplus[phase])/I(m.eps);
 if((I(dl)+gain).lo>=0&&(I(dr)+gain).lo>=0)
  best=std::min(best,chord(a,positive_finite(a.l,m,phase,dl),positive_finite(a.u,m,phase,dr)));
 if(budget_ok(a.l,m,phase,dl)&&budget_ok(a.u,m,phase,dr))
  best=std::min(best,chord(a,(I(negativeR(a.l,m,phase,dl))*I(factor)).hi,(I(negativeR(a.u,m,phase,dr))*I(factor)).hi));
 return best;
}
static double finite_quadratic(double A,double B,double z,double lower,double upper,double r,double DL,double DU,double V){
 if(upper<lower)return -1e100;
 double split=std::min(upper,std::max(lower,-z));I DD=I(DL)*I(DL);
 double ans=-1e100;
 if(lower<=-z)ans=(I(r)*(I(A)+I(B)*I(split))/DD).hi;
 if(upper>=std::max(lower,-z)){
  I KV=DD*I(V)/I(DU),bb=I(B)*(I(r)-I(z))-I(A)-KV,cc=I(A)*(I(r)-I(z))-KV*I(z);
  double q;
  if(bb.hi<(I(2)*I(B)*I(split)).lo){I h(split);q=(-I(B)*h*h+I(bb.hi)*h+I(cc.hi)).hi;}
  else if(bb.hi>(I(2)*I(B)*I(upper)).hi){I h(upper);q=(-I(B)*h*h+I(bb.hi)*h+I(cc.hi)).hi;}
  else q=(I(cc.hi)+I(bb.hi)*I(bb.hi)/(I(4)*I(B))).hi;
  ans=std::max(ans,(I(q)/DD).hi);
 }
 return ans;
}
static double finite_tail(const Config&q,const FiniteMass&m){
 double V=q.terminal.hi;
 for(int j=q.n-1;j>=q.prefix;j--){
  const Cell&a=q.cells[j];double dl=std::min(dlower(a.l,m),dlower(a.u,m));
  I z(0);for(int g=0;g<4;g++)z=z+I(std::min((zcoef[g][a.col]*I(m.lo[g])).lo,(zcoef[g][a.col]*I(m.hi[g])).lo));
  double dd=(I(dl>=0?m.scaleL:m.eps)*I(dl)).lo;
  double zz=(I(z.lo>=0?m.scaleL:m.eps)*I(z.lo)).lo;
  double sigma=std::max(globalR(a.l,m),globalR(a.u,m));sigma=(I(sigma)/I(m.F)).hi;
  if(budget_ok(a.l,m,a.completed,dl))sigma=std::min(sigma,(I(negativeR(a.l,m,a.completed,dl))/I(m.F)).hi);
  double hmax=(I(m.HU[4])+a.u.r-rational_interval(32,5)).hi;
  double H=finite_quadratic((I(1)-I(dd)).hi,1.,zz,m.HL[4],hmax,a.u.r.hi,a.l.D.lo,a.u.D.hi,V);
  if(dd<0){double cap=std::min(hmax,sqrt_upper((-I(2)*I(dd)).hi));
   if(cap>=m.HL[4])H=std::max(H,finite_quadratic(sigma,sigma,zz,m.HL[4],cap,a.u.r.hi,a.l.D.lo,a.u.D.hi,V));}
  I step=a.dt*I(H);if(H<0)step=step/(I(1)+a.dt*a.u.r/I(a.l.D.lo));V=(I(V)+step).hi;
 }
 return V;
}
extern "C" double evaluate_finite(const double *lo,const double *hi,int den,double scaleL,double scaleU,double*diagnostic){
 if(den<=0||den%5||scaleL<0||scaleU<=0||scaleL>scaleU||scaleU>1)return 1e30;
 const Config&q=config(den);FiniteMass m(lo,hi,scaleL,scaleU);I value(0);double survival=1.;
 for(int j=0;j<q.prefix;j++){
  const Cell&a=q.cells[j];I kb(0);
  for(int g=0;g<4;g++){
   I co=zcoef[g][a.col]+I(g<a.completed?1:0);
   kb=kb+I(std::min((co*I(m.lo[g])).lo,(co*I(m.hi[g])).lo));
  }
  double kl=std::max(0.,(I(scaleL)*I(kb.lo)).lo),factor=1.;
  if(kl>0)factor=std::max(0.,(I(1)-I(kl)/I(a.u.r.hi)).hi);
  value=value+I(survival)*I(factor)*I(finite_prefix_cell(a,m));
  survival=(I(survival)/(I(1)+a.dt*I(kl)/I(a.u.D.hi))).hi;
 }
 double t=finite_tail(q,m);
 if(diagnostic){diagnostic[0]=value.hi;diagnostic[1]=t;diagnostic[2]=survival;}
 return (GAMMA+value+I(survival)*I(t)).hi;
}
