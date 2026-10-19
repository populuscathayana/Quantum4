"""Selected numeric fchk fields and independently projected vibrational Hessian."""
from pathlib import Path
import re
import numpy as np

BOHR_A=0.529177210903
HARTREE_KCAL=627.5094740631

def fchk_field(text,label):
    lines=text.splitlines()
    pattern=re.compile(r'^'+re.escape(label)+r'\s+([IR])\s+(?:N=\s*(\d+)|(.+))$')
    for i,line in enumerate(lines):
        m=pattern.match(line)
        if not m:continue
        conv=int if m[1]=='I' else lambda x:float(x.replace('D','E'))
        if m[3]:return conv(m[3].strip())
        n=int(m[2]);values=[]
        for following in lines[i+1:]:
            values.extend(conv(x) for x in following.split())
            if len(values)>=n:break
        assert len(values)==n,(label,n,len(values))
        return np.array(values)
    raise KeyError(label)

def stationary_data(path):
    text=Path(path).read_text()
    get=lambda name:fchk_field(text,name)
    z=get('Atomic numbers');n=len(z);mass=get('Real atomic weights')
    xyz=get('Current cartesian coordinates').reshape(n,3)
    grad=get('Cartesian Gradient').reshape(n,3)
    lower=get('Cartesian Force Constants');h=np.zeros((3*n,3*n))
    h[np.tril_indices(3*n)]=lower;h=h+h.T-np.diag(h.diagonal())
    centered=xyz-np.average(xyz,axis=0,weights=mass)
    rigid=[]
    for axis in np.eye(3):
        rigid.append((np.broadcast_to(axis,(n,3))*np.sqrt(mass)[:,None]).ravel())
        rigid.append((np.cross(axis,centered)*np.sqrt(mass)[:,None]).ravel())
    u,s,_=np.linalg.svd(np.array(rigid).T,full_matrices=True)
    rank=int(np.sum(s>1e-8));assert rank==6
    q=u[:,rank:];sqrt_m=np.repeat(np.sqrt(mass),3)
    hmw=h/np.outer(sqrt_m,sqrt_m)
    eigenvalues,eigenvectors=np.linalg.eigh(q.T@hmw@q)
    # sqrt(E_h / (a_0^2 u)) / (2*pi*c), with c in cm/s.
    factor=np.sqrt(4.3597447222071e-18/(5.29177210903e-11**2*1.66053906660e-27))/(2*np.pi*29979245800)
    frequencies=np.sign(eigenvalues)*np.sqrt(np.abs(eigenvalues))*factor
    modes=(q@eigenvectors)/sqrt_m[:,None]
    return {'Z':z,'mass_amu':mass,'xyz_A':xyz*BOHR_A,'gradient_Eh_bohr':grad,
       'energy_hartree':get('Total Energy'),'S2':get('S**2'),
       'S2_after':get('S**2 after annihilation'),'electrons':get('Number of electrons'),
       'alpha':get('Number of alpha electrons'),'beta':get('Number of beta electrons'),
       'frequencies_cm_1':frequencies,'cartesian_modes':modes,'rigid_rank':rank,
       'fchk_frequencies_cm_1':get('Vib-E2')[:3*n-6]}
