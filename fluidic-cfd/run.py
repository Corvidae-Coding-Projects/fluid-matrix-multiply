"""Reproduce the CFD runs. Requires Python 3, NumPy, and g++ with OpenMP."""
from pathlib import Path
import argparse,os,subprocess

ROOT=Path(__file__).resolve().parent
p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--case',default='all',choices=['all','e0','e1','e2','e3','mix','signed','zero','duct','custom'])
p.add_argument('--scale',type=int,default=1,choices=[1,2,3,4])
p.add_argument('--threads',type=int,default=min(8,os.cpu_count() or 1))
p.add_argument('--pressure',type=float,default=1.,help='Pressure scale in Pa; normalized inputs are multiplied by this')
p.add_argument('--steps',type=int,default=None)
a=p.parse_args();os.chdir(ROOT);(ROOT/'results').mkdir(exist_ok=True)
if a.case=='all' and a.pressure<=0:p.error('Matrix calibration requires a positive pressure scale')
subprocess.run(['g++','-O3','-fopenmp','-std=c++17','solver.cpp','-o','solver-rebuilt'],check=True)
env={**os.environ,'OMP_NUM_THREADS':str(a.threads)}
def run(case,scale,tol='1e-6'):
    prefix=f'results/{case}-s{scale}'
    command=['./solver-rebuilt',str(scale),case,prefix,str(a.steps or 40000*scale*scale),tol,str(a.pressure)]
    subprocess.run(command,env=env,check=True)
if a.case=='all':
    for c in ['e0','e1','e2','e3','mix','signed','zero']:run(c,1)
    run('mix',2)
    for s in [1,2]:run('duct',s,'1e-8')
    subprocess.run(['python3','analyze.py'],check=True)
    subprocess.run(['python3','build-viewer.py'],check=True)
else:
    run(a.case,a.scale,'1e-8' if a.case=='duct' else '1e-6')
