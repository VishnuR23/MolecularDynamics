#!/usr/bin/env python3
"""Independent Ewald/SPC-E reference, used to cross-check the C++ engine.

Deliberately slow and naive: O(N^2) real space, an explicit triple loop over
k-vectors, no optimisation. It exists to be obviously correct rather than
fast, so that agreement with src/moldyn/ is evidence rather than coincidence
-- two implementations sharing a clever trick can share its bug.

Reproduces all 28 published values (4 configurations x 6 energy components,
plus totals) from the NIST Standard Reference Simulation Website. Run it:

    python3 tools/spce_reference.py

One convention is not stated on the NIST page and is easy to get wrong:
intramolecular distances need the minimum image. Some molecules in the
sample configurations are stored split across the periodic boundary -- the
oxygen on one side, its hydrogens wrapped to the other -- so raw separations
reach 19.9 A where the bond is 1.0 A. Using them leaves Eintra 4.8% low
while every other term still looks correct.

No dependencies beyond the standard library.
"""
import math
e=1.602176565e-19; eps0=8.854187817e-12; kB=1.3806488e-23
COUL=e*e/(4*math.pi*eps0)*1e10/kB
q=0.42380; SIG=3.16555789; EPS=78.19743111; RC=10.0

def load(p):
    lines=open(p).read().splitlines()
    L=[float(v) for v in lines[0].split()]; M=int(lines[1])
    at=[(float(f[1]),float(f[2]),float(f[3]),f[4]) for f in (l.split() for l in lines[2:]) if len(f)>=5]
    return L,M,at

def energies(path):
    L,M,at=load(path); V=L[0]*L[1]*L[2]; N=len(at)
    def mi(d,k): return d-L[k]*round(d/L[k])
    chg=[(-2*q if a[3]=='O' else q) for a in at]
    pos=[(a[0],a[1],a[2]) for a in at]
    alpha=5.6/min(L)
    O=[i for i in range(N) if at[i][3]=='O']
    Ed=0.0
    for ii in range(len(O)):
        xi,yi,zi=pos[O[ii]]
        for jj in range(ii+1,len(O)):
            xj,yj,zj=pos[O[jj]]
            dx=mi(xi-xj,0); dy=mi(yi-yj,1); dz=mi(zi-zj,2)
            r2=dx*dx+dy*dy+dz*dz
            if r2>=RC*RC: continue
            s6=(SIG*SIG/r2)**3
            Ed+=4*EPS*(s6*s6-s6)
    rho=M/V; s3=(SIG/RC)**3
    Elrc=(8/3)*math.pi*M*rho*EPS*SIG**3*(s3**3/3-s3)
    Er=0.0
    for i in range(N):
        mi_=i//3; xi,yi,zi=pos[i]; ci=chg[i]
        for j in range(i+1,N):
            if j//3==mi_: continue
            xj,yj,zj=pos[j]
            dx=mi(xi-xj,0); dy=mi(yi-yj,1); dz=mi(zi-zj,2)
            r2=dx*dx+dy*dy+dz*dz
            if r2>=RC*RC: continue
            r=math.sqrt(r2)
            Er+=COUL*ci*chg[j]*math.erfc(alpha*r)/r
    kmax=5; Ef=0.0
    for nx in range(-kmax,kmax+1):
        for ny in range(-kmax,kmax+1):
            for nz in range(-kmax,kmax+1):
                n2=nx*nx+ny*ny+nz*nz
                if n2==0 or n2>=27: continue
                kx,ky,kz=2*math.pi*nx/L[0],2*math.pi*ny/L[1],2*math.pi*nz/L[2]
                k2=kx*kx+ky*ky+kz*kz
                re=im=0.0
                for i in range(N):
                    ph=kx*pos[i][0]+ky*pos[i][1]+kz*pos[i][2]
                    re+=chg[i]*math.cos(ph); im+=chg[i]*math.sin(ph)
                Ef+=COUL*(2*math.pi/V)*math.exp(-k2/(4*alpha*alpha))/k2*(re*re+im*im)
    Es=-COUL*alpha/math.sqrt(math.pi)*sum(c*c for c in chg)
    Ei=0.0
    for m in range(M):
        i=3*m; cs=[-2*q,q,q]
        for a,b in ((0,1),(0,2),(1,2)):
            dx=mi(pos[i+a][0]-pos[i+b][0],0); dy=mi(pos[i+a][1]-pos[i+b][1],1); dz=mi(pos[i+a][2]-pos[i+b][2],2)
            r=math.sqrt(dx*dx+dy*dy+dz*dz)
            Ei-=COUL*cs[a]*cs[b]*math.erf(alpha*r)/r
    return Ed,Elrc,Er,Ef,Es,Ei,Ed+Elrc+Er+Ef+Es+Ei

REF={1:(9.95387e4,-8.23715e2,-5.58889e5,6.27009e3,-2.84469e6,2.80999e6,-4.88604e5),
     2:(1.93712e5,-3.29486e3,-1.19295e6,6.03495e3,-5.68938e6,5.61998e6,-1.06590e6),
     3:(3.54344e5,-7.41343e3,-1.96297e6,5.24461e3,-8.53407e6,8.42998e6,-1.71488e6),
     4:(4.48593e5,-1.37286e4,-3.57226e6,7.58785e3,-1.42235e7,1.41483e7,-3.20501e6)}
names=["Edisp","ELRC","Ereal","Efourier","Eself","Eintra","Etotal"]
print(f"{'cfg':>3} {'term':>9} {'ours':>15} {'NIST':>15} {'rel dev':>10}")
worst=0.0
for c in (1,2,3,4):
    got=energies(f'data/nist/spce/spce_sample_config_periodic{c}.txt')
    for n,g,r in zip(names,got,REF[c]):
        dev=abs(g-r)/abs(r); worst=max(worst,dev)
        print(f"{c:>3} {n:>9} {g:15.6g} {r:15.6g} {dev:10.2e}")
print(f"\nworst relative deviation across all 28 values: {worst:.2e}")
print("(NIST prints 6 significant figures, so ~1e-6 is the rounding floor)")
