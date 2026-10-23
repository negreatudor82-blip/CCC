#!/usr/bin/env python3
"""Optional MILP refinement of existing Level 3 solutions.

The standalone C++ solver does not require this script. This optional offline
helper requires NumPy and SciPy. It improves complete solutions from exact
prefix inventories, retains each incumbent unless exact replay finds a cheaper
valid replacement, and atomically writes complete output paragraphs.

The MILP is a bounded aggregate planning heuristic, not a proof of schedule
feasibility or global optimality. Startup constraints use exactly conserved
chemical elements. Every proposed schedule is separately replayed; reported
purchase costs always come from that replay, never the relaxation objective.

Example:
  python3 optimizer/refine_with_milp.py inputs/in_level-3_4-enormous.txt \
      solutions/out_level-3_4-enormous.txt /tmp/refined.txt \
      --seconds-per-solve 5 --rounds 2
"""
from __future__ import annotations
import os
# Keep this optional refinement from taking all cores when multiple inputs run.
os.environ.setdefault('OMP_NUM_THREADS','1')
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
import argparse,re,time,json,warnings,tempfile,math
from pathlib import Path
from collections import defaultdict,Counter
import numpy as np
from scipy.optimize import milp,Bounds,LinearConstraint
from scipy.sparse import csc_matrix
warnings.filterwarnings('ignore',message='Unrecognized options detected')

def parse(path):
 lines=iter(s.strip() for s in Path(path).read_text().splitlines() if s.strip())
 prices={};reactions=[];names=set()
 for _ in range(int(next(lines))):
  p,n=next(lines).split(' -> ');prices[n]=int(p);names.add(n)
 def side(text):
  result={}
  for term in text.split(','):
   q,n=term.split();q=int(q)
   if q<=0:raise ValueError('Reaction quantities must be positive')
   result[n]=result.get(n,0)+q;names.add(n)
  return result
 for _ in range(int(next(lines))):
  text=next(lines);a,b=text.split(' -> ');reactions.append((side(a),side(b),text))
 targets=[]
 for _ in range(int(next(lines))):
  k,n=next(lines).split();targets.append((int(k),n));names.add(n)
 if list(lines):raise ValueError('Unexpected trailing input')
 if any(p<0 for p in prices.values()):raise ValueError('Negative prices are unsupported')
 return prices,reactions,targets,sorted(names)

class Optimizer:
 def __init__(self,path,base):
  self.p,self.rx,self.cases,self.names=parse(str(path));self.idx={n:i for i,n in enumerate(self.names)}
  self.bn=list(self.p);self.nr=len(self.rx);self.nb=len(self.bn);self.nm=len(self.names)
  self.rid={s:i for i,(_,_,s) in enumerate(self.rx)}
  self.best=[[self.rid[s] for s in block.splitlines()] for block in Path(base).read_text().strip().split('\n\n')]
  if len(self.best)!=len(self.cases):raise ValueError('Incumbent target paragraph count does not match input')
  self.costs=[self.replay(s,self.cases[i])[0] for i,s in enumerate(self.best)]
  if any(c is None for c in self.costs):raise ValueError('Incumbent contains an invalid reaction sequence')
  atoms={}
  for n in self.names:
   d=defaultdict(int);parts=re.findall(r'([A-Z][a-z]*)([0-9]*)',n)
   if ''.join(a+b for a,b in parts)!=n:raise ValueError(n)
   for e,q in parts:d[e]+=int(q or 1)
   atoms[n]=d
  self.pools=[]
  for e in set(e for a in atoms.values() for e in a)-{'C','H','O','N','S','Cl'}:
   if not all(sum(q*atoms[n][e] for n,q in a.items())==sum(q*atoms[n][e] for n,q in b.items()) for a,b,_ in self.rx):continue
   needs=[sum(q*atoms[n][e] for n,q in a.items()) for a,_,_ in self.rx]
   for q in sorted(set(needs)-{0}):self.pools.append((e,q,[r for r,v in enumerate(needs) if v>=q]))
  self.atoms=atoms
 def replay(self,seq,case=None,start=None):
  inv=defaultdict(int,start or {});cost=0
  for r in seq:
   for n,q in self.rx[r][0].items():
    short=max(0,q-inv[n])
    if short and n not in self.p:return None,inv
    cost+=short*self.p.get(n,0);inv[n]+=short-q
   for n,q in self.rx[r][1].items():inv[n]+=q
  if case and inv[case[1]]<case[0]:return None,inv
  return cost,inv
 def solve(self,initial,target,k,support,limit,full=False,cuts=None):
  # Buy variables are integer; reactions are integer. Activation stock cuts
  # are generic conserved-element necessities, not schedule feasibility claims.
  active=[]
  for e,q,rs in self.pools:
   existing=sum(v*self.atoms[n][e] for n,v in initial.items())
   if existing<q:active.append((e,q,rs,existing))
  nv=self.nr+self.nb+len(active)
  mat=np.zeros((self.nm,nv));rhs=np.zeros(self.nm)
  for n,v in initial.items():rhs[self.idx[n]]=-v
  rhs[self.idx[target]]+=k
  for r,(a,b,_) in enumerate(self.rx):
   for n,q in a.items():mat[self.idx[n],r]-=q
   for n,q in b.items():mat[self.idx[n],r]+=q
  for j,n in enumerate(self.bn):mat[self.idx[n],self.nr+j]=1
  extra=[];lb=[]
  for ai,(e,q,rs,existing) in enumerate(active):
   pos=self.nr+self.nb+ai
   for r in rs:
    if not full and not support.get(r):continue
    row=np.zeros(nv);row[r]=-1;row[pos]=5000;extra.append(row);lb.append(0)
   row=np.zeros(nv);row[pos]=-q
   for j,n in enumerate(self.bn):row[self.nr+j]=self.atoms[n][e]
   extra.append(row);lb.append(-existing)
  for reactions,molecules in cuts or []:
   row=np.zeros(nv)
   for r in reactions:row[r]=1
   for n in molecules:row[self.nr+self.bn.index(n)]=1
   extra.append(row);lb.append(1)
  A=csc_matrix(np.vstack([mat,extra]) if extra else mat)
  obj=np.r_[np.full(self.nr,1e-5),[self.p[n] for n in self.bn],np.zeros(len(active))]
  upper=np.r_[[5000 if full else max(0,2*support.get(r,0)+2) if support.get(r) else 0 for r in range(self.nr)],np.full(self.nb,np.inf),np.ones(len(active))]
  res=milp(obj,integrality=np.ones(nv),bounds=Bounds(0,upper),constraints=LinearConstraint(A,np.r_[rhs,lb],np.inf),options={'time_limit':limit,'mip_rel_gap':0,'threads':1})
  if res.x is None:return None
  return np.rint(res.x[:self.nr]).astype(int).tolist(),dict(zip(self.bn,np.rint(res.x[self.nr:self.nr+self.nb]).astype(int))),res.fun,getattr(res,'mip_dual_bound',None)
 def schedule(self,counts,budget,initial,ci,base,mode=0):
  # Use conceptual purchases to find an order, then separately replay the
  # actual reaction-only output with just-in-time purchases.
  k,target=self.cases[ci];stock=defaultdict(int,initial)
  for n,q in budget.items():stock[n]+=q
  left=list(counts);seq=[];source={n:-1 for n in self.p};val={n:float(v) for n,v in self.p.items()}
  for r in base:
   a,b,_=self.rx[r]
   if all(n in source for n in a):
    cc=sum(q*val[n] for n,q in a.items())
    for n,q in b.items():
     if n not in source:source[n]=r;val[n]=cc/q
  positions=defaultdict(list)
  for pos,r in reversed(list(enumerate(base))):positions[r].append(pos)
  def fire(r):
   for n,q in self.rx[r][0].items():
    if stock[n]<q:return False
   for n,q in self.rx[r][0].items():stock[n]-=q
   for n,q in self.rx[r][1].items():stock[n]+=q
   seq.append(r)
   if left[r]>0:left[r]-=1
   if positions[r]:positions[r].pop()
   return True
  def ensure(n,q,stack):
   if stock[n]>=q:return True
   if n in self.p:stock[n]=q;return True
   if n in stack or source.get(n,-1)<0:return False
   r=source[n];stack=stack|{n}
   for _ in range(10000):
    if stock[n]>=q:return True
    for _gather in range(100):
     for t,v in self.rx[r][0].items():
      if not ensure(t,v,stack):return False
     if all(stock[t]>=v for t,v in self.rx[r][0].items()):break
    else:return False
    if not fire(r):return False
   return False
  for _ in range(10000):
   if stock[target]>=k:break
   todo=[r for r,v in enumerate(left) if v>0]
   if not todo:
    if not ensure(target,k,set()):return None
    break
   ready=[]
   for r in todo:
    a,b,_=self.rx[r]
    if all(stock[n]>=q for n,q in a.items()):
     order=positions[r][-1] if positions[r] else 999999
     helpful=sum(q*val.get(n,0) for n,q in b.items())
     ready.append(((order if mode==0 else -helpful if mode==1 else r),r))
   if ready:fire(min(ready)[1]);continue
   if mode==0:
    r=min(todo,key=lambda r:positions[r][-1] if positions[r] else 999999)
    choices=[((q-stock[n])*val.get(n,1e30),n,q) for n,q in self.rx[r][0].items() if stock[n]<q]
   else:choices=[((q-stock[n])*val.get(n,1e30),n,q) for r in todo for n,q in self.rx[r][0].items() if stock[n]<q]
   if not choices:return None
   _,n,q=min(choices)
   if not ensure(n,q,set()):return None
  else:return None
  c,iv=self.replay(seq,self.cases[ci],initial)
  return (c,seq) if c is not None else None
 def write(self,out):
  # Serialize fully and replace atomically; interrupted refinement retains the
  # last validated output rather than leaving a partial file.
  output=Path(out);output.parent.mkdir(parents=True,exist_ok=True)
  text='\n\n'.join('\n'.join(self.rx[r][2] for r in seq) for seq in self.best)+'\n'
  with tempfile.NamedTemporaryFile('w',encoding='utf-8',dir=output.parent,delete=False) as handle:
   temp=Path(handle.name);handle.write(text)
  try:os.replace(temp,output)
  finally:
   if temp.exists():temp.unlink()

 def optimize(self,out,limit=1,only=None,rounds=1,fractions=None,restricted_only=False):
  self.write(out)
  for iteration in range(rounds):
   round_before=sum(self.costs)
   for ci in range(len(self.cases)):
    if only is not None and ci not in only:continue
    k,target=self.cases[ci]
    before=self.costs[ci]
    for fraction in fractions or [0.8,0.6,0.4,0.2,0.0]:
     base=list(self.best[ci]);cut=int(len(base)*fraction);prefix=base[:cut]
     paid,stock=self.replay(prefix)
     if stock[target]>=k:continue
     for full in ([False] if restricted_only else [False,True]):
      result=self.solve(stock,target,k,Counter(base[cut:]),limit,full)
      if result is None:continue
      counts,budget,obj,bound=result
      for mode in [0,1,2]:
       sched=self.schedule(counts,budget,stock,ci,base,mode)
       if sched is None:continue
       cost,tail=sched;total=paid+cost
       if total<self.costs[ci]:
        candidate=prefix+tail;actual,_=self.replay(candidate,self.cases[ci]);assert actual==total
        self.best[ci]=candidate;self.costs[ci]=total;self.write(out)
        print('IMPROVE',ci+1,'from',before,'to',total,'cut',fraction,'full',full,'obj',round(obj,3),'bound',bound,flush=True)
    print('CASE',ci+1,self.costs[ci],'total',sum(self.costs),flush=True)
   if sum(self.costs)==round_before:break
  return self.costs

def main():
 parser=argparse.ArgumentParser(description=__doc__)
 parser.add_argument('input',type=Path)
 parser.add_argument('incumbent',type=Path)
 parser.add_argument('output',type=Path)
 parser.add_argument('--seconds-per-solve',type=float,default=5)
 parser.add_argument('--rounds',type=int,default=1)
 parser.add_argument('--restricted-only',action='store_true',help='Keep reactions already used by the incumbent; faster local refinement')
 parser.add_argument('--fractions',default='0.8,0.6,0.4,0.2,0',help='Comma-separated retained prefix fractions')
 parser.add_argument('--cases',help='Optional comma-separated, 1-based case numbers')
 args=parser.parse_args()
 if args.input.resolve()==args.output.resolve():parser.error('Input and output paths must differ')
 if not math.isfinite(args.seconds_per_solve) or args.seconds_per_solve<=0:parser.error('Time budget must be positive')
 if args.rounds<1:parser.error('At least one round is required')
 optimizer=Optimizer(args.input,args.incumbent)
 selected=set(int(x)-1 for x in args.cases.split(',')) if args.cases else None
 if selected is not None and any(i<0 or i>=len(optimizer.cases) for i in selected):parser.error('Case number out of range')
 before=sum(optimizer.costs)
 fractions=[float(x) for x in args.fractions.split(',')]
 if not fractions or any(not math.isfinite(f) or f<0 or f>=1 for f in fractions):parser.error('Prefix fractions must be in [0,1)')
 optimizer.optimize(args.output,args.seconds_per_solve,selected,args.rounds,fractions,args.restricted_only)
 print('Verified total:',before,'->',sum(optimizer.costs))

if __name__=='__main__':main()
