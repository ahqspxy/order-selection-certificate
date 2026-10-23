#!/usr/bin/env python3
"""Exact-rational cross-checks of the certificate's arithmetic primitives.

This supports the mathematical arithmetic proof; it is not a substitute for
the full regional replay. Only the Python standard library and C++17 are used.
"""
import argparse
import ctypes as C
from fractions import Fraction as F
import hashlib
import json
import math
from pathlib import Path
import random
import subprocess

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--package',type=Path,required=True)
    parser.add_argument('--out',type=Path,required=True)
    args=parser.parse_args()
    root,out=args.package.resolve(),args.out.resolve()
    out.mkdir(parents=True,exist_ok=True)
    src=out/'arithmetic_oracle_api.cpp'
    src.write_text('''#include "normalized_cert.cpp"
extern "C" void oracle_basic(int op,double a,double b,double c,double d,double*out){
 strict_mode=true;I x(a,b),y(c,d),z;
 if(op==0)z=x+y;else if(op==1)z=x-y;else if(op==2)z=x*y;else z=x/y;
 out[0]=z.lo;out[1]=z.hi;
}
extern "C" void oracle_exp(double x,double*out){strict_mode=true;I y=expI(I(x));out[0]=y.lo;out[1]=y.hi;}
extern "C" double oracle_log(double x){strict_mode=true;return log1p_lower(x);}
extern "C" double oracle_sqrt(double x){strict_mode=true;return sqrt_upper(x);}
''')
    binary=out/'arithmetic_oracle.so'
    command=['g++','-std=c++17','-O3','-fno-fast-math','-ffp-contract=off','-fPIC','-shared',
        str(src),'-I',str(root/'source/guarded/robust_745'),'-o',str(binary)]
    subprocess.run(command,check=True,capture_output=True)
    lib=C.CDLL(str(binary)); ptr=C.POINTER(C.c_double)
    lib.oracle_basic.argtypes=[C.c_int,*([C.c_double]*4),ptr]
    lib.oracle_exp.argtypes=[C.c_double,ptr]
    for name in ['oracle_log','oracle_sqrt']:
        f=getattr(lib,name);f.argtypes=[C.c_double];f.restype=C.c_double
    rng=random.Random(74500001)
    values=[0.,-0.,1.,-1.,math.nextafter(1.,0.),math.nextafter(1.,2.),
            2.**-1022,2.**-1074,2.**-500,2.**500,1e-21,1e-9,0.74500001]
    values += [math.ldexp(rng.uniform(-1,1),rng.randint(-500,500)) for _ in range(400)]
    count=0
    for _ in range(4000):
        a,b=sorted(rng.sample(values,2));c,d=sorted(rng.sample(values,2))
        for op in range(4):
            if op==3 and (c<=0 or c<2.**-500):continue
            A,B,Cc,D=map(F,(a,b,c,d))
            exact=([A+Cc,B+D] if op==0 else [A-D,B-Cc] if op==1 else
                   [x*y for x in (A,B) for y in (Cc,D)] if op==2 else
                   [x/y for x in (A,B) for y in (Cc,D)])
            # Ignore cases outside the finite domain required by the guards.
            if max(map(abs,exact))>F(float.fromhex('0x1.fffffffffffffp+1022')):continue
            result=(C.c_double*2)();lib.oracle_basic(op,a,b,c,d,result)
            assert F(result[0])<=min(exact)<=max(exact)<=F(result[1])
            count+=1
    def exp_positive(x):
        term=total=F(1)
        for n in range(1,201):term*=x/n;total+=term
        rest=term*x/201/(1-x/202)
        return total,total+rest
    exp_inputs=[0.,2.**-52,1e-12,1e-6,.2545,.4,1.,6.4,10.,20.]
    exp_inputs += [rng.uniform(0,20) for _ in range(70)]
    for x in exp_inputs:
        lo,hi=exp_positive(F(x))
        for arg,L,U in [(x,lo,hi),(-x,1/hi,1/lo)]:
            result=(C.c_double*2)();lib.oracle_exp(arg,result)
            assert F(result[0])<=L<=U<=F(result[1]),arg
    def log_series(z):
        term=total=z
        for n in range(1,100):term*=z*z;total+=term/(2*n+1)
        return 2*total
    ln2lower=log_series(F(1,3))
    log_inputs=[0.,2.**-1074,2.**-1022,2.**-52,1e-12,.001,.4,1.,20.,1e12]
    log_inputs += [2.**rng.uniform(-40,40) for _ in range(80)]
    for x in log_inputs:
        v=1+F(x);exponent=0
        while v>=2:v/=2;exponent+=1
        lower=exponent*ln2lower+log_series((v-1)/(v+1))
        value=lib.oracle_log(x)
        assert 0<=F(value)<=lower,x
    sqrt_inputs=[0.,2.**-1022,2.**-500,1e-30,1e-16,.4,1.,20.,2.**500]
    sqrt_inputs += [2.**rng.uniform(-1000,500) for _ in range(100)]
    for x in sqrt_inputs:
        q=lib.oracle_sqrt(x)
        assert q>=0 and F(q)*F(q)>=F(x),x
    report={
        'status':'PASS_EXACT_RATIONAL_ARITHMETIC_CROSS_CHECKS',
        'basic_interval_cases':count,'exponential_cases':2*len(exp_inputs),
        'logarithm_lower_bound_cases':len(log_inputs),'square_root_upper_bound_cases':len(sqrt_inputs),
        'reference_arithmetic':'Python Fraction, exact 200-term exponential bounds and 100-term positive logarithm lower bounds',
        'source_header_sha256':hashlib.sha256((root/'source/guarded/higher_ratio/higher_interval.hpp').read_bytes()).hexdigest(),
        'compiler_command':command,
        'scope':'Primitive cross-checks plus the manuscript arithmetic argument; not a full regional replay.',
    }
    (out/'arithmetic_check.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))

if __name__=='__main__':main()
