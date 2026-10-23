#pragma once
#include <cmath>
#include <algorithm>
#include <limits>
#include <cfenv>
#include <cstdint>
#include <cstring>

#include <cstdio>
#include <cstdlib>
static unsigned long long audit_intervals=0,audit_divisions=0,audit_exp=0,audit_sqrt=0,audit_sqrt_steps_max=0;
static double audit_div_min=INFINITY,audit_exp_max=0,audit_sqrt_min=INFINITY,audit_sqrt_max=0;
static void audit_require(bool ok,const char*what){if(!ok){std::fprintf(stderr,"CERT_GUARD_FAILURE: %s\n",what);std::fflush(stderr);std::abort();}}
static void audit_interval(double l,double u){audit_intervals++;audit_require(std::isfinite(l)&&std::isfinite(u)&&l<=u,"finite ordered interval");}
extern "C" void audit_reset(){audit_intervals=audit_divisions=audit_exp=audit_sqrt=audit_sqrt_steps_max=0;audit_div_min=audit_sqrt_min=INFINITY;audit_exp_max=audit_sqrt_max=0;}
extern "C" void audit_get_counts(unsigned long long*out){out[0]=audit_intervals;out[1]=audit_divisions;out[2]=audit_exp;out[3]=audit_sqrt;out[4]=audit_sqrt_steps_max;}
extern "C" void audit_get_ranges(double*out){out[0]=audit_div_min;out[1]=audit_exp_max;out[2]=audit_sqrt_min;out[3]=audit_sqrt_max;}

static bool strict_mode=false;
static inline double adjacent(double x,bool upward){
 std::uint64_t bits;std::memcpy(&bits,&x,sizeof(bits));
 const std::uint64_t mag=bits&UINT64_C(0x7fffffffffffffff);
 if(mag!=0&&(mag&UINT64_C(0x7ff0000000000000))!=UINT64_C(0x7ff0000000000000)){
  bool negative=(bits>>63)!=0;
  if(upward!=negative)++bits;else --bits;
  std::memcpy(&x,&bits,sizeof(x));return x;
 }
 return std::nextafter(x,upward?INFINITY:-INFINITY);
}
static double dn(double x){return strict_mode?adjacent(x,false):x;}
static double up(double x){return strict_mode?adjacent(x,true):x;}
struct I{double lo,hi;I(double x=0):lo(x),hi(x){audit_interval(x,x);}I(double l,double u):lo(l),hi(u){audit_interval(l,u);}};
static I operator+(I a,I b){if(a.lo==0&&a.hi==0)return b;if(b.lo==0&&b.hi==0)return a;return {dn(a.lo+b.lo),up(a.hi+b.hi)};}
static I operator-(I a){return {-a.hi,-a.lo};}
static I operator-(I a,I b){return a+-b;}
static I operator*(I a,I b){if((a.lo==0&&a.hi==0)||(b.lo==0&&b.hi==0))return I(0);double p[]={a.lo*b.lo,a.lo*b.hi,a.hi*b.lo,a.hi*b.hi};return {dn(*std::min_element(p,p+4)),up(*std::max_element(p,p+4))};}
static I operator/(I a,I b){audit_divisions++;audit_require(b.lo>0 && std::isfinite(b.hi),"strictly positive denominator");audit_div_min=std::min(audit_div_min,b.lo);if(a.lo==0&&a.hi==0)return I(0);return a*I(dn(1/b.hi),up(1/b.lo));}
static I exprange_positive(double l,double u){audit_exp++;audit_require(0<=l&&l<=u&&u<98,"positive exponential Taylor range < 98");audit_exp_max=std::max(audit_exp_max,u);
 double tl=1,tu=1,sl=1,su=1;
 for(int n=1;n<=96;n++){tl=dn(dn(tl*l)/n);tu=up(up(tu*u)/n);sl=dn(sl+tl);su=up(su+tu);}
 double tail=up(up(tu*u)/97);tail=up(tail/dn(1-up(u/98)));
 return {sl,up(su+tail)};
}
static I expI(I x){
 if(!strict_mode)return {std::exp(x.lo),std::exp(x.hi)};
 if(x.lo==0&&x.hi==0)return I(1);
 if(x.lo>=0)return exprange_positive(x.lo,x.hi);
 if(x.hi<=0){I e=exprange_positive(-x.hi,-x.lo);return I(1)/e;}
 return {expI(I(x.lo)).lo,expI(I(x.hi)).hi};
}
static double log1p_lower(double x){
 if(!strict_mode)return std::log1p(x);
 double t=dn(1+x);if(t<=1)return 0;
 int exponent;double y=std::frexp(t,&exponent)*2;exponent--;
 auto series=[](double z){double term=z,sum=z,z2=dn(z*z);for(int n=1;n<20;n++){term=dn(term*z2);sum=dn(sum+dn(term/(2*n+1)));}return dn(2*sum);};
 double z=dn(dn(y-1)/up(y+1));double ln2=series(dn(1./3));
 return std::max(0.,dn(dn(exponent*ln2)+series(std::max(0.,z))));
}
static I gamma_interval(double gamma){if(!strict_mode)return I(gamma);if(gamma==.74301)return I(74301)/I(100000);if(gamma==.7431)return I(7431)/I(10000);if(gamma==.74401)return I(74401)/I(100000);if(gamma==.74501)return I(74501)/I(100000);if(gamma==.744001)return I(744001)/I(1000000);if(gamma==.745001)return I(745001)/I(1000000);return I(dn(gamma),up(gamma));}
static I rational_interval(int n,int d){return I(n)/I(d);}
static I wi[16][25],ci[16][26];
