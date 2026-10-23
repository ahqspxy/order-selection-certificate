// Scalar continuation with state-dependent normalized inverse bounds.
#include "normalized_finite.cpp"
static std::map<std::pair<double,double>,I> refined_exp_cache;
static I refined_exp(I x){
 auto key=std::make_pair(x.lo,x.hi);auto it=refined_exp_cache.find(key);
 if(it!=refined_exp_cache.end())return it->second;
 I y=expI(x);refined_exp_cache.emplace(key,y);return y;
}
#define expI refined_exp
struct RefinedMass:FiniteMass {
 double normDen[5],absDen[5];
 RefinedMass(const double*a,const double*b,double sl,double su):FiniteMass(a,b,sl,su){
  normDen[0]=absDen[0]=0;I total(0);
  for(int g=1;g<=4;g++){
   total=total+I(hi[g-1]);
   normDen[g]=sl>0?((I(1)-expI(-I(sl)*I(total.hi)))/I(sl)).hi:total.hi;
   absDen[g]=std::min(HU[g],(I(1)-expI(-I(HU[g]))).hi);
  }
 }
};
static double anchor_cap(const Row&r,const RefinedMass&m,int g){
 I e=g<4?r.capacity[g]:r.prefix[3]+I(1);
 return (e*I(m.EH[g])).hi;
}
struct Pair {double a,b;};
static double positive_state(const Row&r,const FiniteMass&m,int phase,double dl,double eh,double gain){
 I target=I(dl)+I(gain);I u[4];
 for(int g=0;g<4;g++)u[g]=g<=phase?I(eh):r.capacity[g]*I(m.EH[g]);
 double ans=(expI(r.r-rational_interval(32,5))*I(m.EH[4])).hi;
 for(int mask=1;mask<16;mask++){
  I mass(0),sum(0);for(int g=0;g<4;g++)if(mask&(1<<g)){mass=mass+I(m.EW[g]);sum=sum+I(m.EW[g])*u[g];}
  if(mass.lo>0)ans=std::max(ans,((sum-target)/I(mass.lo)).hi);
 }
 return ans;
}
static double absolute_lower(double dl,const FiniteMass&m){
 return (I(dl>=0?m.scaleL:m.eps)*I(dl)).lo;
}
static double mass_cone_low(double t,const FiniteMass&m){
 double ans=0;
 for(int g=0;g<=4;g++){
  I gap=rational_interval(2*ends[g],5)-I(t);
  ans=std::max(ans,(I(m.HL[g])-I(std::max(0.,gap.hi))).lo);
 }
 return ans;
}
static Pair state_R(const Cell&a,const RefinedMass&m,double hl,double hu){
 int phase=a.completed;double dl=dlower(a.l,m),dr=dlower(a.u,m);
 double factor=((I(1)+I(hu))/I(m.F)).hi;
 Pair best{(I(globalR(a.l,m))*I(factor)).hi,(I(globalR(a.u,m))*I(factor)).hi};
 double score=chord(a,best.a,best.b);
 auto use=[&](Pair p){double v=chord(a,p.a,p.b);if(v<score){score=v;best=p;}};
 if(budget_ok(a.l,m,phase,dl)&&budget_ok(a.u,m,phase,dr))
  use({(I(negativeR(a.l,m,phase,dl))*I(factor)).hi,(I(negativeR(a.u,m,phase,dr))*I(factor)).hi});
 I h(hl);double gain=std::max(0.,((h*h/I(2)+h*h*h/I(6))/I(m.eps)).lo);
 if((I(dl)+I(gain)).lo>=0&&(I(dr)+I(gain)).lo>=0){
  double ddl=(I(dl>=0?m.scaleL:m.eps)*I(dl)).lo,ddr=(I(dr>=0?m.scaleL:m.eps)*I(dr)).lo;
  use({(I(1)+I(hu)-I(ddl)).hi,(I(1)+I(hu)-I(ddr)).hi});
  double eh=expI(I(hu)).hi;
  use({positive_state(a.l,m,phase,dl,eh,gain),positive_state(a.u,m,phase,dr,eh,gain)});
  // Truncate the rank envelope at each future source boundary. Each
  // before/after maximum is convex in exp(r) on this state/time cell.
  double surplus=std::max(0.,(h*h/I(2)+h*h*h/I(6)).lo);
  double deltaL=(I(surplus)+I(absolute_lower(dl,m))).lo;
  double deltaR=(I(surplus)+I(absolute_lower(dr,m))).lo;
  if(deltaL>=0&&deltaR>=0)for(int g=phase+1;g<=4;g++){
   I anchor=rational_interval(2*ends[g],5);
   if(a.u.r.lo>anchor.hi||m.HU[g]<=0)continue;
   double denom=m.absDen[g];
   if(denom<=0)continue;
   double left=std::max((I(eh)-I(deltaL)/I(denom)).hi,anchor_cap(a.l,m,g));
   double right=std::max((I(eh)-I(deltaR)/I(denom)).hi,anchor_cap(a.u,m,g));
   use({left,right});
  }
  // The normalized form remains valid down to scale zero: the slope
  // denominator divided by s is at most (1-exp(-s_L*C_g^+))/s_L.
  for(int g=phase+1;g<=4;g++)if(m.normDen[g]>0){
   double left=std::max((I(eh)-(I(dl)+I(gain))/I(m.normDen[g])).hi,anchor_cap(a.l,m,g));
   double right=std::max((I(eh)-(I(dr)+I(gain))/I(m.normDen[g])).hi,anchor_cap(a.u,m,g));
   use({left,right});
  }

 }
 // On the negative inverse branch, a positive lower bound for h(tau)
 // gives (h(tau)-d)/(1-exp(-h(tau))). Its dependence on d is affine.
 // Keep the positive-branch basic bound by taking their maximum.
 double dmin=std::min(dl,dr);
 double sigma=(I(std::max(globalR(a.l,m),globalR(a.u,m)))/I(m.F)).hi;
 if(budget_ok(a.l,m,phase,dmin))sigma=std::min(sigma,(I(negativeR(a.l,m,phase,dmin))/I(m.F)).hi);
 double shift=std::max(0.,(I(sigma)-I(1)).hi);
 double t=std::max(0.,(a.l.r-I(shift)).lo);
 double hlTau=std::max(mass_cone_low(t,m),std::max(0.,(I(hl)-I(shift)).lo));
 if(hlTau>0){
  I x(hlTau);double den=std::max(0.,(x+x*x/I(2)+x*x*x/I(6)).lo);
  if(den>0){
   double al=absolute_lower(dl,m),ar=absolute_lower(dr,m);
   double lcap=((x-I(al))*(I(1)+I(1)/I(den))).hi;
   double rcap=((x-I(ar))*(I(1)+I(1)/I(den))).hi;
   use({std::max((I(1)+I(hu)-I(al)).hi,lcap),std::max((I(1)+I(hu)-I(ar)).hi,rcap)});
  }
 }
 return best;
}
extern "C" double evaluate_refined(const double*lo,const double*hi,int den,int hbins,double scaleL,double scaleU,double*diagnostic){
 if(den<=0||den%5||hbins<=0||scaleL<0||scaleU<=0||scaleL>scaleU||scaleU>1)return 1e30;
 refined_exp_cache.clear();
 const Config&q=config(den);RefinedMass m(lo,hi,scaleL,scaleU);double V=finite_tail(q,m);double terminal=V;
 for(int j=q.prefix-1;j>=0;j--){
  const Cell&a=q.cells[j];int p=a.completed;
  I begin=rational_interval(2*ends[p],5),end=rational_interval(2*ends[p+1],5);
  double hlow=std::max(m.HL[p],(I(m.HL[p+1])+a.l.r-end).lo);
  double hhigh=std::min(m.HU[p+1],(I(m.HU[p])+a.u.r-begin).hi);
  hlow=std::max(0.,hlow);hhigh=std::max(hlow,hhigh);
  I dz(0);for(int g=0;g<4;g++)dz=dz+I(std::min((zcoef[g][a.col]*I(m.lo[g])).lo,(zcoef[g][a.col]*I(m.hi[g])).lo));
  double zl=(I(dz.lo>=0?scaleL:scaleU)*I(dz.lo)).lo;
  double Ua=-1e100,Ub=-1e100;
  for(int b=0;b<hbins;b++){
   I l=I(hlow)+(I(hhigh)-I(hlow))*rational_interval(b,hbins);
   I u=I(hlow)+(I(hhigh)-I(hlow))*rational_interval(b+1,hbins);
   double hl=std::max(0.,l.lo),hu=u.hi;
   Pair R=state_R(a,m,hl,hu);
   double kl=std::max(0.,(I(hl)+I(zl)).lo);
   if(kl>a.u.r.hi)continue; // this bin cannot contain an admissible state
   I factor=I(1)-I(kl)/I(a.u.r.hi);
   I credit=I(kl)*I(V)*I(a.l.D.lo)/I(a.u.r.hi);
   Ua=std::max(Ua,(factor*I(R.a)-credit).hi);Ub=std::max(Ub,(factor*I(R.b)-credit).hi);
  }
  double value=chord(a,Ua,Ub);
  I delta=a.dt*a.u.r/I(a.l.D.lo);
  I error=I(std::min(1.,delta.hi))*I(std::max(0.,-std::min(Ua,Ub)))*a.I0;
  V=(I(V)+I(value)+error).hi;
 }
 if(diagnostic){diagnostic[0]=terminal;diagnostic[1]=q.prefix;diagnostic[2]=hbins;}
 return (GAMMA+I(V)).hi;
}
