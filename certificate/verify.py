#!/usr/bin/env python3
"""Check the complete certificate, or rebuild and recompute every numerical region.

python3 verify.py                         # exact geometry, kernel and stored bounds
python3 verify.py --recompute --out ../run # full interval replay, no stored bounds used

Only Python's standard library and a C++17 compiler are required. Run without -O.
The replay output directory must be outside the distributed package.
"""
import sys
sys.dont_write_bytecode = True
import argparse, collections, concurrent.futures as cf, ctypes as C
import itertools, json, lzma, math, multiprocessing as mp, os, platform
from pathlib import Path
import subprocess, time
from fractions import Fraction as F
import interval_engine as E
from kernel_check import verify_kernel
from helpers.partition_audit import partition
from helpers.geometry_audit import geometry, validate_mapping

ORDER=('boundary','atlas','global')
FILES={k:k+'_cover.json.xz' for k in ORDER}
EXPECTED_COUNTS={'boundary':32529,'atlas':377458,'global':222958}
NUMERIC_COUNTS={'boundary':32529,'atlas':377458,'global':222958}
GLOBAL_METHODS={'refined_dp','fasttight','jointtight','jointtime'}
EXPECTED_GLOBAL_METHODS={'refined_dp':175903,'fasttight':4515,'jointtight':8765,'jointtime':33775}

def read_json(path):
    with lzma.open(path,'rt') if path.suffix=='.xz' else path.open() as f:
        return json.load(f)

def identities(root):
    return {'inputs':{k:E.sha(root/'data'/FILES[k]) for k in ORDER},
            'sources':{str(p.relative_to(root)):E.sha(p) for p in sorted((root/'source').rglob('*')) if p.is_file()},
            'engine':E.sha(root/'interval_engine.py')}

def check_manifest(root):
    manifest=read_json(root/'SHA256SUMS.json')
    actual={str(p.relative_to(root)) for p in root.rglob('*') if p.is_file() and p!=root/'SHA256SUMS.json'}
    E.require(actual==set(manifest),'Manifest file set mismatch')
    for name,digest in manifest.items():
        path=Path(name)
        E.require(not path.is_absolute() and '..' not in path.parts and not (root/path).is_symlink(),'Unsafe manifest path')
        E.require(E.sha(root/path)==digest,'Hash mismatch: '+name)
    return len(manifest)

def read_covers(root):
    result={}
    for k in ORDER:
        doc=read_json(root/'data'/FILES[k])
        E.require(doc['schema']=='compact_order_selection_'+k+'_v1','Wrong descriptor schema')
        E.require(F(doc['gamma'])==E.GAMMA,'Wrong target')
        result[k]=doc['rows']
        E.require(len(result[k])==EXPECTED_COUNTS[k],'Wrong number of '+k+' rows')
    return result

def valid_resolution(res,scalar=False):
    if scalar:
        E.require(type(res) is int and res>0 and res%5==0,'Invalid scalar resolution')
    else:
        E.require(type(res) is list and len(res)==2 and all(type(x) is int and x>0 for x in res) and res[0]%5==0,'Invalid paired resolution')

def valid_box(lo,hi,pivot=None):
    E.require(len(lo)==len(hi)==4,'Wrong parameter dimension')
    E.require(all(type(x) in (int,float) and math.isfinite(x) for x in lo+hi),'Nonfinite parameter')
    E.require(all(0<=F(a)<=F(b)<=1 for a,b in zip(lo,hi)),'Invalid parameter box')
    if pivot is not None:
        E.require(type(pivot) is int and 0<=pivot<4 and lo[pivot]==hi[pivot]==1,'Invalid normalized pivot')

def check_geometry(covers):
    result={}
    b=covers['boundary']
    for row in b:
        E.require(set(row)=={'lo','hi','pivot','method','resolution'},'Unexpected boundary fields')
        valid_box(row['lo'],row['hi'],row['pivot'])
        E.require(row['method'] in {'normalized_finite','normalized_dp'},'Unknown boundary method')
        valid_resolution(row['resolution'],row['method']=='normalized_finite')
    result['boundary']=[]
    for pivot in range(4):
        axes=[g for g in range(4) if g!=pivot]
        rows=[r for r in b if r['pivot']==pivot]
        result['boundary'].append({'pivot':pivot,**partition([([r['lo'][g] for g in axes],[r['hi'][g] for g in axes]) for r in rows],3)})
    groups=collections.defaultdict(list)
    for row in covers['atlas']:
        E.require(set(row)=={'alo','ahi','pivot','exponent','param_lo','param_hi','scale','resolution'},'Unexpected atlas fields')
        validate_mapping(row)
        groups[row['pivot'],row['exponent']].append(row)
    E.require(set(groups)=={(p,e) for p in range(4) for e in range(-14,-10)},'Incomplete atlas roots')
    result['atlas']=[geometry(groups[p,e],p,e) for p in range(4) for e in range(-14,-10)]
    g=covers['global'];methods=collections.Counter()
    for row in g:
        valid_box(row['lo'],row['hi']);methods[row['method']]+=1
        E.require(set(row)=={'lo','hi','method','resolution'},'Unexpected global fields')
        E.require(row['method'] in GLOBAL_METHODS,'Unknown global method')
        E.require(max(map(F,row['lo']))>=F(1,1024),'Global box enters the interior of the small-scale cube')
        valid_resolution(row['resolution'])
    E.require(dict(methods)==EXPECTED_GLOBAL_METHODS,'Global method counts mismatch')
    result['global']=partition([([0]*4,[2.**-10]*4)]+[(r['lo'],r['hi']) for r in g],4)
    result['global_methods']=dict(methods)
    return result

def call_spec(kind,row):
    if kind=='atlas':
        return ('refined_dp',row['alo'],row['ahi'],row['resolution'],*row['scale'])
    if kind=='boundary':
        return (row['method'],row['lo'],row['hi'],row['resolution'],0.,2.**-14)
    return (row['method'],row['lo'],row['hi'],row['resolution'],1.,1.)

def assemble(covers,values):
    E.require(set(values)==set(ORDER),'Incomplete result kinds')
    maxima={};counts={}
    for kind in ORDER:
        E.require(len(values[kind])==len(covers[kind]),'Result length mismatch')
        best=F(0);best_index=None;count=0
        for i,(row,encoded) in enumerate(zip(covers[kind],values[kind])):
            E.require(type(encoded) is str,'Missing numerical endpoint')
            x=float.fromhex(encoded)
            E.require(math.isfinite(x) and x.hex()==encoded,'Noncanonical finite endpoint')
            v=F(x);E.require(0<=v<E.UNIFORM,'Endpoint does not strictly pass the uniform threshold')
            if v>best:best,best_index=v,i
            count+=1
        E.require(count==NUMERIC_COUNTS[kind],'Wrong numerical count')
        maxima[kind]={'index':best_index,'upper_exact':str(best)};counts[kind]=count
    zero=E.analytic_zero_bound()
    small=max(zero,F(maxima['boundary']['upper_exact']),F(maxima['atlas']['upper_exact']))
    final=max(small,F(maxima['global']['upper_exact']))
    E.require(final<E.UNIFORM<1,'Final strict inequality failed')
    return {'numerical_counts':counts,'numerical_calls':sum(counts.values()),'maxima':maxima,
            'analytic_zero_upper_exact':str(zero),'small_scale_upper_exact':str(small),
            'maximum_upper_exact':str(final),'uniform_upper_exact':str(E.UNIFORM),
            'margin_to_one_exact':str(1-final),'implicit_small_scale_cubes':1}

def initialize_worker(root,out,binary_hashes):
    global ENGINES,ENV,ROOT,OUT
    ROOT,OUT=Path(root),Path(out);ENGINES={}
    for method,digest in binary_hashes.items():
        E.require(E.sha(OUT/(method+'.so'))==digest,'Worker binary mismatch')
    ENV=C.CDLL(str(OUT/'certification_environment.so'))
    ENV.certification_environment.argtypes=[];ENV.certification_environment.restype=C.c_int
    E.require(ENV.certification_environment()==0,'Unsupported worker environment')

def run_batch(batch):
    records=[]
    for kind,index,row in batch:
        spec=call_spec(kind,row);method=spec[0]
        if method not in ENGINES:
            ENGINES[method]=E.Engine(method,OUT/(method+'.so'),ROOT/'source/guarded/robust_745/bulk_r1_exact',ENV)
        ans=ENGINES[method].call(*spec[1:])
        E.require(F(ans['upper_exact'])<E.UNIFORM,'Non-strict numerical endpoint')
        records.append([kind,index,E.canonical_hash(row),ans['upper_hex'],ans['seconds']])
    return {'records':records,'pid':os.getpid(),'engines':{k:{'identity':v.identity,'guards':v.guard_statistics()} for k,v in ENGINES.items()}}

def write_results(path,identity,values):
    with lzma.open(path,'wt',preset=6) as f:
        json.dump({'schema':'compact_order_selection_results_v1','identity':identity,'values':values},f,separators=(',',':'))
        f.write('\n')

def replay(root,covers,args):
    out=args.out.resolve();E.require(out!=root and root not in out.parents,'Replay output must be outside the package')
    out.mkdir(parents=True,exist_ok=True)
    E.require(args.resume or not any(out.iterdir()),'Output is not empty; use --resume')
    identity=identities(root)
    run_identity={**identity,'runner':E.sha(__file__),'compiler':subprocess.check_output([args.cxx,'--version'],text=True),'flags':E.FLAGS}
    report_path=out/'replay_report.json';journal=out/'results.jsonl'
    values={k:[None]*len(covers[k]) for k in ORDER};seen=set()
    if args.resume:
        previous=read_json(report_path)
        E.require(previous['run_identity']==run_identity,'Resume inputs or implementation changed')
        with journal.open() as f:
            for line in f:
                kind,i,digest,encoded,seconds=json.loads(line)
                E.require(kind in covers and type(i) is int and 0<=i<len(covers[kind]),'Invalid resumed index')
                E.require((kind,i) not in seen,'Duplicate resumed row')
                E.require(digest==E.canonical_hash(covers[kind][i]),'Resumed input mismatch')
                value=float.fromhex(encoded)
                E.require(math.isfinite(value) and value.hex()==encoded and 0<=F(value)<E.UNIFORM,'Invalid resumed endpoint')
                values[kind][i]=encoded;seen.add((kind,i))
    report={'status':'RUNNING','run_identity':run_identity,'started_utc':E.now(),'workers':args.workers,'resumed_calls':len(seen),'builds':[],'worker_checks':{},'completed_calls':len(seen)}
    E.save(report_path,report);include=root/'source/guarded/robust_745'
    report['builds'].append(E.build(args.cxx,include/'certification_environment.cpp',out/'certification_environment.so',include,False,out))
    for method in E.SOURCES:
        src=out/(method+'_api.cpp');src.write_text(E.wrapper(method))
        report['builds'].append(E.build(args.cxx,src,out/(method+'.so'),include,method in E.POTENTIAL,out))
    hashes={m:E.sha(out/(m+'.so')) for m in [*E.SOURCES,'certification_environment']}
    jobs=((k,i,r) for k in ORDER for i,r in enumerate(covers[k]) if (k,i) not in seen)
    def batches():
        while batch:=list(itertools.islice(jobs,64)):yield batch
    started=time.monotonic();last=started;batch_iterator=iter(batches())
    try:
        with cf.ProcessPoolExecutor(args.workers,mp_context=mp.get_context('spawn'),initializer=initialize_worker,initargs=(str(root),str(out),hashes)) as pool,journal.open('a') as stream:
            pending={pool.submit(run_batch,b) for b in itertools.islice(batch_iterator,2*args.workers)}
            while pending:
                done,pending=cf.wait(pending,return_when=cf.FIRST_COMPLETED)
                for future in done:
                    batch=future.result();report['worker_checks'][str(batch['pid'])]=batch['engines']
                    for record in batch['records']:
                        kind,i,digest,encoded,seconds=record
                        E.require((kind,i) not in seen and digest==E.canonical_hash(covers[kind][i]),'Worker identity mismatch')
                        values[kind][i]=encoded;seen.add((kind,i));stream.write(json.dumps(record,separators=(',',':'))+'\n')
                    try:pending.add(pool.submit(run_batch,next(batch_iterator)))
                    except StopIteration:pass
                if time.monotonic()-last>=30:
                    stream.flush();report.update(completed_calls=len(seen),elapsed_seconds=time.monotonic()-started)
                    E.save(report_path,report);print('RECOMPUTED',len(seen),'/',sum(NUMERIC_COUNTS.values()),flush=True);last=time.monotonic()
        summary=assemble(covers,values)
        E.require(identities(root)==identity,'Inputs changed during replay')
        write_results(out/'results.json.xz',identity,values)
        report.update(status='PASS_COMPLETE_NUMERICAL_REPLAY',completed_calls=len(seen),elapsed_seconds=time.monotonic()-started,finished_utc=E.now(),summary=summary,results_sha256=E.sha(out/'results.json.xz'),journal_sha256=E.sha(journal))
        E.save(report_path,report);return report
    except BaseException as exc:
        report.update(status='INCOMPLETE',completed_calls=len(seen),error=repr(exc),elapsed_seconds=time.monotonic()-started)
        E.save(report_path,report);raise

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--recompute',action='store_true');parser.add_argument('--out',type=Path)
    parser.add_argument('--resume',action='store_true');parser.add_argument('--workers',type=int,default=6)
    parser.add_argument('--cxx',default='g++');parser.add_argument('--report',type=Path)
    args=parser.parse_args();E.require(__debug__,'Run without Python -O')
    E.require(args.workers>0,'Invalid worker count');E.require(not args.resume or args.recompute,'Resume requires recompute')
    root=Path(__file__).resolve().parent;started=time.monotonic()
    file_count=check_manifest(root);covers=read_covers(root)
    kernel=verify_kernel(root);geometry_result=check_geometry(covers)
    if args.recompute:
        E.require(args.out is not None,'Full replay requires --out')
        report=replay(root,covers,args)
        report['geometry']=geometry_result;report['kernel']=kernel
        E.save(args.out.resolve()/'replay_report.json',report)
    else:
        data=read_json(root/'results.json.xz')
        E.require(data['schema']=='compact_order_selection_results_v1' and data['identity']==identities(root),'Stored result/input identity mismatch')
        report={'status':'PASS_COMPLETE_STORED_CERTIFICATE','actual_integrations_recomputed':0,'files_checked':file_count,'kernel':kernel,'geometry':geometry_result,'summary':assemble(covers,data['values']),'elapsed_seconds':time.monotonic()-started}
    if args.report:
        path=args.report.resolve();E.require(path!=root and root not in path.parents,'Report must be outside package');E.save(path,report)
    print(json.dumps({k:report[k] for k in ['status','summary']},indent=2))

if __name__=='__main__':main()
