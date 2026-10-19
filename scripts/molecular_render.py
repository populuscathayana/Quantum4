"""Source-backed static figures and actual IRC/Hessian motion animations."""
import argparse
import hashlib
import io
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from PIL import Image
from gaussian_data import stationary_data,HARTREE_KCAL

ROOT=Path(__file__).resolve().parents[1]
FIG=ROOT/'experiment-report/figures'
ANIM=ROOT/'analysis/visuals'
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,'axes.spines.top':False,'axes.spines.right':False,
                     'svg.fonttype':'none','savefig.facecolor':'white'})
BLUE='#235b8e';GOLD='#bf8b29';GRAY='#30363d';PINK='#bc537b'

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def align(x,ref):
    x=np.asarray(x);ref=np.asarray(ref)
    a=x-x.mean(axis=0);b=ref-ref.mean(axis=0)
    u,_,vt=np.linalg.svd(a.T@b)
    d=np.diag([1,1,np.linalg.det(u@vt)]);rot=u@d@vt
    assert np.linalg.det(rot)>.999999
    return a@rot

def molecule(ax,xyz,z,title,limits,mode=None):
    carbon=[i for i,v in enumerate(z) if v==6]
    # Fixed chemical skeleton; C1-C2 opening is shown as a dashed contact.
    edges=[(0,2),(1,2),(1,8),(0,3),(0,4),(1,5),(2,6),(2,7),(8,9),(8,10)]
    for i,j in edges:
        p=xyz[[i,j]];ax.plot(*p.T,color='#646c76',lw=3,solid_capstyle='round')
    p=xyz[[0,1]];ax.plot(*p.T,color=GOLD,lw=2,ls='--')
    for i in np.argsort(xyz[:,2]):
        c=z[i]==6;ax.scatter(*xyz[i],s=195 if c else 90,c=GRAY if c else '#eef1f4',edgecolors='#737c85',linewidths=.7,depthshade=True)
    carbon_center=xyz[carbon].mean(axis=0)
    for i in carbon:
        p=xyz[i].copy();offset=p-carbon_center;offset[2]=0
        p+=.32*offset/max(np.linalg.norm(offset),.1);p[2]+=.12
        ax.text(*p,f'C{i+1}',ha='center',va='center',fontsize=11,color=BLUE,fontweight='bold')
    if mode is not None:
        for i in carbon:
            ax.quiver(*xyz[i],*mode[i],length=.75,color=PINK,linewidth=2,arrow_length_ratio=.25)
    ax.set(xlim=limits[0],ylim=limits[1],zlim=limits[2]);ax.set_box_aspect([1,1,1],zoom=1.6)
    ax.view_init(elev=75,azim=-90);ax.set_proj_type('ortho');ax.set_axis_off();ax.set_title(title,fontsize=12,pad=3)

def frame(x,z,title,limits):
    fig=plt.figure(figsize=(6.4,5),dpi=100);ax=fig.add_subplot(projection='3d')
    molecule(ax,x,z,title,limits)
    d12=np.linalg.norm(x[0]-x[1]);d29=np.linalg.norm(x[1]-x[8])
    fig.text(.5,.035,f'C1-C2 = {d12:.3f} Å     C2-C9 = {d29:.3f} Å',ha='center',color=GRAY)
    fig.subplots_adjust(left=0,right=1,bottom=.07,top=.91)
    buf=io.BytesIO();fig.savefig(buf,format='png');plt.close(fig);buf.seek(0)
    return Image.open(buf).convert('RGB')
