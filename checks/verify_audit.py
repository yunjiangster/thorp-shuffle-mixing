#!/usr/bin/env python3
"""Independent finite checks for the Thorp V3 audit / V4 improvement.
Exact combinatorial assertions use integers/Fraction. Tests involving logarithms
are floating-point regressions, explicitly separated in the output.
Requires Python 3.10+, NumPy. No network access; deterministic random seed.
"""
from __future__ import annotations
import collections, functools, itertools, json, math, pathlib, time
from fractions import Fraction as Q
import numpy as np
ROOT=pathlib.Path(__file__).resolve().parent
L=math.log(2.0)

@functools.lru_cache(None)
def weights(j:int)->tuple[Q,...]:
    if j==0:return (Q(1),)
    a=1<<(j.bit_length()-1);i=j-a
    return (Q(1,2*a),)*a+tuple(x/2 for x in weights(i))

def networks(d:int):
    n=1<<d
    edges=[[(x,x^(1<<i)) for x in range(n) if not(x>>i&1)] for i in range(d)]
    for bits in range(1<<(d*n//2)):
        pos=list(range(n)); adjacency=[0]*n;offset=0
        for layer in edges:
            inv=[0]*n
            for a,x in enumerate(pos):inv[x]=a
            for x,y in layer:
                a,b=inv[x],inv[y]
                adjacency[a]|=1<<b;adjacency[b]|=1<<a
            image=list(range(n))
            for j,(x,y) in enumerate(layer):
                if bits>>(offset+j)&1:image[x],image[y]=image[y],image[x]
            pos=[image[x] for x in pos];offset+=n//2
        yield tuple(pos),tuple(adjacency)

@functools.lru_cache(None)
def candidates(permutation:tuple[int,...], j:int):
    """Recursive C_j and symmetric Q_j recovered from a realized endpoint map.
    The tests independently verify dependence ONLY on the exposed higher endpoints.
    Q_j is represented by immutable tuples of rational transition rows.
    """
    n=len(permutation)
    if n==1:return (0,),((Q(1),),)
    a=n//2
    upper=tuple(permutation[a+i]%a for i in range(a))
    lower=tuple(permutation[i]%a for i in range(a))
    if j<a:
        bits={upper[i]:1-permutation[a+i]//a for i in range(a)}
        C,M=candidates(lower,j)
        return tuple(x+a*bits[x] for x in C),M
    i=j-a;C0,M0=candidates(upper,i)
    C=tuple(x+a*b for x in C0 for b in range(2));index={x:r for r,x in enumerate(C)}
    rows=[]
    for x in C0:
        xi=C0.index(x)
        for b in range(2):
            row=[Q(0)]*len(C)
            row[index[x+a*(1-b)]]+=Q(a-i,a+i)
            if i:
                for yi,y in enumerate(C0):
                    for c in range(2):row[index[y+a*c]]+=Q(i,a+i)*M0[xi][yi]
            rows.append(tuple(row))
    return C,tuple(rows)

def exact_geometry():
    report={};total_iso=total_cand=total_na=total_route=0
    for d in (1,2,3):
        n=1<<d;net=list(networks(d));M=len(net)
        assert len(set(p for p,g in net))==M, 'One-sweep injectivity failed'
        for perm,graph in net:
            assert all(row.bit_count()==d for row in graph)
            assert sum(row.bit_count() for row in graph)==n*d
            for B in range(1,1<<n):
                v=B.bit_count();e=sum((graph[x]&B).bit_count() for x in range(n) if B>>x&1)//2
                assert 1<<(2*e)<=v**v
                total_iso+=1
        occ=np.zeros((1<<n,1<<n),dtype=np.int64)
        for perm,graph in net:
            images=[0]*(1<<n)
            for mask in range(1,1<<n):
                low=mask&-mask;images[mask]=images[mask^low]|(1<<perm[low.bit_length()-1])
            occ[np.arange(1<<n),images]+=1
        inc=occ.copy()
        for i in range(n):
            for B in range(1<<n):
                if not(B>>i&1):inc[:,B]+=inc[:,B|(1<<i)]
        for A in range(1<<n):
            k=A.bit_count()
            for B in range(1<<n):
                v=B.bit_count()
                assert int(inc[A,B])*n**v<=M*k**v
                total_na+=1
        for j in range(n):
            grouped={};kweights=[Q(0)]*(j+1);W=[Q(0)]*(j+1)
            for perm,graph in net:
                C,Qmat=candidates(perm,j);H=perm[j+1:]
                index={x:r for r,x in enumerate(C)}
                canonical=(tuple(sorted(C)),tuple((x,y,Qmat[index[x]][index[y]]) for x in sorted(C) for y in sorted(C)))
                if H not in grouped:grouped[H]=[canonical,collections.Counter()]
                assert grouped[H][0]==canonical, ('nonmeasurable',d,j,H)
                grouped[H][1][perm[j]]+=1
                assert len(C)==1<<j.bit_count()
                assert perm[j] in C and not(set(C)&set(H))
                assert all(sum(row)==1 for row in Qmat)
                assert all(Qmat[r][s]==Qmat[s][r] for r in range(len(C)) for s in range(len(C)))
                if j:
                    assert all(Qmat[r][r]==0 for r in range(len(C)))
                    inv={x:i for i,x in enumerate(perm)}
                    for r,b in enumerate(C):
                        w=Qmat[index[perm[j]]][r]
                        if w:assert inv[b]<j;kweights[inv[b]]+=w
                for k in range(j+1):
                    if perm[k] in C:W[k]+=Q(1,len(C))
                total_cand+=1
            for canonical,counts in grouped.values():
                assert set(counts)==set(canonical[0]) and len(set(counts.values()))==1
            assert tuple(w/M for w in W)==weights(j)
            if j:assert all(w/M==Q(1,j) for w in kweights[:j]) and kweights[j]==0
        for h in range(1,min(3,n)+1):
            omega=list(itertools.permutations(range(n),h))
            for x in omega:
                counts=collections.Counter(tuple(perm[a] for a in x) for perm,g in net)
                for y in omega:
                    feasible=True;collisions=0
                    for i in range(d+1):
                        loc=[(y[a]&((1<<i)-1))|(x[a]&~((1<<i)-1)) for a in range(h)]
                        if len(set(loc))!=h:feasible=False;break
                        if i<d:
                            collisions+=sum((loc[a]^(1<<i))==loc[b] for a in range(h) for b in range(a+1,h))
                    if feasible:
                        exp=d*n//2-d*h+collisions
                        assert exp>=0 and counts[y]==1<<exp
                    else:assert counts[y]==0
                    total_route+=1
        report[str(d)]={'cards':n,'networks':M,'distinct_sweep_permutations':M}
    report.update(induced_subgraph_checks=total_iso,candidate_network_label_checks=total_cand,
                  inclusion_inequalities=total_na,ordered_routing_pairs=total_route)
    return report

@functools.lru_cache(None)
def partitions(n:int,cap:int|None=None):
    if n==0:return ((),)
    if cap is None:cap=n
    return tuple((k,)+rest for k in range(min(n,cap),1-1,-1) for rest in partitions(n-k,k))
@functools.lru_cache(None)
def dim(shape:tuple[int,...])->int:
    prod=1
    for r,row in enumerate(shape):
        for c in range(row):prod*=row-c+sum(s>c for s in shape[r+1:])
    return math.factorial(sum(shape))//prod
@functools.lru_cache(None)
def mult(shape:tuple[int,...],remaining:int)->int:
    if sum(shape)==remaining:return int(shape==(remaining,) or not shape and remaining==0)
    ans=0
    for r in range(len(shape)):
        if r+1==len(shape) or shape[r]>shape[r+1]:
            a=list(shape);a[r]-=1;a=tuple(x for x in a if x)
            if remaining==0 or a and a[0]>=remaining:ans+=mult(a,remaining)
    return ans

def exact_representation():
    total=identities=0
    for n in range(2,25):
        for k in range(1,n//2+1):
            for mu in partitions(k):
                shape=(n-k,)+mu
                assert dim(shape)<=math.comb(n,k)*dim(mu)
                for h in range(k,n-k+1):
                    assert mult(shape,n-h)==math.comb(h,k)*dim(mu)
                    total+=1
    for n in range(2,13):
        for h in range(n+1):
            assert sum(dim(lam)*mult(lam,n-h) for lam in partitions(n))==math.factorial(n)//math.factorial(n-h)
            identities+=1
    return {'exact_multiplicities':total,'ordered_module_dimension_identities':identities}

def kl(p,q):
    mask=p>0
    if np.any(q[mask]<=0):return float('inf')
    return float(np.dot(p[mask],np.log(p[mask]/q[mask])))

def scalar_and_entropy():
    rng=np.random.default_rng(20261002);num=0;min_margin=math.inf
    for ell in (2,3,4,8,16,64,256):
        u=np.full(ell,1/ell)
        laws=[u]
        for mass in (1e-9,1e-5,.001,.01,.1,.5,.9,1-1e-10,1.):
            p=(1-mass)*u.copy();p[0]+=mass;laws.append(p)
        for conc in (.001,.03,.3,1.,10.,1000.):
            for _ in range(25):
                x=rng.gamma(conc,1,size=ell)
                if x.sum():laws.append(x/x.sum())
        for p in laws:
            V=max(0.,1-np.sqrt(p).sum()**2/ell)
            for q in (p,u,laws[int(rng.integers(len(laws)))]):
                C=kl(p,q)
                if not math.isfinite(C):continue
                E=kl(p,u)
                for beta in (.2,.5,math.log2(1.5)-.001):
                    R=ell**beta*float(np.sum(q**(1+beta)))
                    for delta in (.001,.1,.25,.75):
                        c=1+1/beta;a=1/(beta*(1-delta));b=math.log(8/delta)+2+c
                        rhs=c*C+(a*math.log(max(R,1)/V)+b)*V if V else c*C
                        margin=rhs-E;assert margin>=-2e-11
                        min_margin=min(min_margin,margin);num+=1
    wm=0
    for beta in (.1,.3,.5,math.log2(1.5)-.001):
        Cbeta=1/(3-2**(1+beta))
        for d in range(1,9):
            n=1<<d
            total=sum((j+1)**beta*sum(float(w)**(1+beta) for w in weights(j)) for j in range(n))
            assert total<=Cbeta*n*(1+1e-12);wm+=1
    perms=list(itertools.permutations(range(4)));index={p:i for i,p in enumerate(perms)}
    invs=[]
    for p in perms:
        inv=[0]*4
        for a,x in enumerate(p):inv[x]=a
        invs.append(tuple(inv))
    transition=np.zeros((24,24))
    for i,p in enumerate(perms):
        for x,g in networks(2):transition[i,index[tuple(x[p[a]] for a in range(4))]]+=1/16
    inputs=[np.full(24,1/24)]
    for i in range(24):
        for eps in (1e-8,.001,.1,.9,1.):
            v=np.full(24,(1-eps)/24);v[i]+=eps;inputs.append(v)
    for c in (.003,.1,1.,100.):
        for _ in range(50):
            v=rng.gamma(c,1,size=24)
            if v.sum():inputs.append(v/v.sum())
    losschecks=0;min_loss=math.inf
    for law in inputs:
        E=kl(law,np.full(24,1/24));out=law@transition;Delta=E-kl(out,np.full(24,1/24))
        S=Csum=chain=0.
        for j in range(1,4):
            buckets=collections.defaultdict(list)
            for i,inv in enumerate(invs):
                if law[i]>0:buckets[inv[j+1:]].append(i)
            for ids in buckets.values():
                mass=sum(law[i] for i in ids);p=np.zeros((j+1,4))
                for i in ids:
                    for k in range(j+1):p[k,invs[i][k]]+=law[i]/mass
                pj=p[j];q=np.array([float(w) for w in weights(j)])@p
                V=max(0.,1-np.sqrt(pj).sum()**2/(j+1))
                S+=mass*V;Csum+=mass*kl(pj,q)
                mask=pj>0;chain+=mass*float(np.sum(pj[mask]*np.log((j+1)*pj[mask])))
        assert abs(chain-E)<2e-11
        assert Delta+2e-11>=L*S and Delta+2e-11>=Csum
        if E>1e-10:
            assert E<=Delta*(100+12*max(0,math.log(4/E)))+1e-10
        min_loss=min(min_loss,Delta-L*S,Delta-Csum);losschecks+=1
    return {'scalar_tests':num,'minimum_scalar_slack':min_margin,'endpoint_moment_tests':wm,
            'S4_input_laws':losschecks,'minimum_two_loss_slack':min_loss,
            'arithmetic':'IEEE double; tolerance 2e-11, not exact certificates'}

def improved_clock(d:int,s:float,eta:float):
    J=math.ceil((s+math.log(1/eta))/L);Z=L*J
    beta=math.log2(1.5)-1/(Z+4);delta=1/(Z+4)
    c=1+1/beta;a=1/(beta*(1-delta));b=math.log(8/delta)+2+c
    Cbeta=1/(3-2**(beta+1));D=a/L
    F=c+b/L+D*(max(0,math.log(Cbeta*L))+math.log(100+12*Z)+1/math.e)
    R=math.ceil(100*max(0,math.log(d*L)))+sum(math.ceil(L*(F+D*j*L)) for j in range(1,J+1))
    return R,{'J':J,'D':D,'F':F,'beta':beta,'delta':delta}

def parameters():
    A=1/(2*math.log(1.5));z=(17+math.sqrt(97))/8
    C=(2+1/(z-2))*A**(1/3)*(z*L)**(2/3)
    report=[]
    for d in (10**4,10**6,10**8,10**10,10**12):
        for eps in (.25,.01):
            eta=eps/8;v=(z*L*d/A)**(1/3);s=v+4*math.log(d+2)+100
            logH=d*L-s;assert logH>10
            ma=2+math.ceil(z*L*(d+1)/(2*v))
            mb=2+math.ceil(z*L*(d+1)/(2*(z-2)*v))
            R,info=improved_clock(d,s,eta);t=max(R,2*ma);u=2*mb
            logaN=math.log(8*L*d)+L*(d+1)/(2*(ma-1))
            logbN=math.log(8*L*d)+L*(d+1)/(2*(mb-1))
            logdelta=math.log(16)+3+logaN+.5*logbN-.5*s
            assert logaN>=math.log(16) and logbN>=math.log(16)
            assert logaN-s<math.log(1/16) and logbN-s<math.log(1/16)
            assert logdelta<0 and s>10
            logH=d*L
            assert logH-s>=math.log(2) and logH-s+math.log(eta)>=0
            assert logH-s>=math.log(1/(math.e*eta)+math.log(1/(1-eta)))
            assert s>math.log(10)
            low=(1-eta)**-2*math.exp(logdelta)/(-math.expm1(logdelta))
            high_upper=math.exp(4.5-s/2)
            finite_TV_upper=2*eta+.5*math.sqrt(low+high_upper)
            assert finite_TV_upper<eps
            report.append({'d':d,'epsilon':eps,'s':s,'log_delta_upper':logdelta,
                           'finite_TV_upper_using_relaxed_high_term':finite_TV_upper,
                           'R':R,'regular_block_sweeps':t,'actual_last_block_sweeps':u,
                           'scaled_layers':(2*t+u)/d**(2/3),**info})
    return {'leading_constant':C,'equal_three_block_constant':3*A**(1/3)*(3*L)**(2/3),
            'previous_constant':18*(L/2)**(2/3),'entropy_clock_coefficient':A,'z':z,
            'tests':report,'arithmetic':'floating point log-domain checks, not simulation'}

if __name__=='__main__':
    started=time.time();report={}
    for name,fn in [('exact_geometry',exact_geometry),('exact_representation',exact_representation),
                    ('floating_scalar_and_entropy',scalar_and_entropy),('floating_parameters',parameters)]:
        print('Running',name,flush=True);report[name]=fn();print(name,'PASS',flush=True)
    report['elapsed_seconds']=time.time()-started;report['status']='PASS'
    (ROOT/'audit_results.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))
