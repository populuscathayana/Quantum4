"""Regenerate reviewed scientific figures into a fresh output directory."""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from verify_public_data import ROOT,load
from gaussian_data import stationary_data
from molecular_render import molecule,frame,align

BLUE='#235b8e';GOLD='#bb852d';GRAY='#68727c';PINK='#aa4c73'

def render(dest,animations=False):
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,'axes.spines.top':False,'axes.spines.right':False})
    def save(fig,name):
        fig.savefig(dest/name,dpi=220);plt.close(fig)
    scan=load('analysis/report-review-20261019/pm7-complete-audit.json')
    rows=scan['joined_points'];x=[r['C1_C2_A'] for r in rows];y=[r['relative_kcal_mol'] for r in rows]
    fig,ax=plt.subplots(figsize=(8,4.8),layout='constrained')
    ax.plot(x,y,'o-',color=BLUE);k=scan['opening_peak_index'];ax.scatter(x[k],y[k],color=GOLD,s=80,zorder=3)
    ax.annotate(f'Opening-region maximum\n{x[k]:.2f} Å; {y[k]:.2f} kcal/mol',xy=(x[k],y[k]),xytext=(1.58,13),arrowprops={'arrowstyle':'-','color':GRAY})
    ax.axvline(2.56,color=GRAY,ls=':',lw=1);ax.text(2.75,13,'Repaired tail',ha='center',color=GRAY)
    ax.set(xlabel='Constrained C1–C2 distance (Å)',ylabel='PM7 relative electronic energy (kcal/mol)',ylim=(min(y)-1,max(y)+2))
    ax.grid(axis='y',alpha=.15);save(fig,'pm7-scan-complete.png')
    spin=load('analysis/irc-20261019/spin-audit.json')['rows'];q=[r['q_sqrt_amu_bohr'] for r in spin]
    fig,axes=plt.subplots(2,1,figsize=(8,7),sharex=True,layout='constrained')
    axes[0].plot(q,[r['relative_E_kcal_mol'] for r in spin],color=BLUE)
    axes[0].set(ylabel='Electronic energy (kcal/mol)',title='IRC: optimized ring radical R = 0')
    for i,c,ls in [(0,BLUE,'-'),(1,PINK,'--'),(2,'#72833f','-.'),(8,GOLD,':')]:
        axes[1].plot(q,[r['Mulliken_populations'][i]['spin'] for r in spin],color=c,ls=ls,label=f'C{i+1}')
    axes[1].set(xlabel='Mass-weighted IRC coordinate (√amu bohr)',ylabel='Mulliken spin population');axes[1].legend(ncol=4,frameon=False)
    for ax in axes:ax.grid(axis='y',alpha=.15);ax.axvline(0,color=GRAY,ls=':',lw=.8)
    save(fig,'irc-energy-spin.png')
    fig,axes=plt.subplots(2,1,figsize=(8,6.5),sharex=True,layout='constrained')
    for b,c,ls in [('C1-C2',BLUE,'-'),('C2-C9',GOLD,'--'),('C1-C3',GRAY,':'),('C2-C3',PINK,'-.')]:
        axes[0].plot(q,[r['distances_A'][b] for r in spin],color=c,ls=ls,label=b)
    axes[0].set_ylabel('Bond distance (Å)');axes[0].legend(ncol=2,frameon=False)
    axes[1].plot(q,[r['S2'] for r in spin],color=BLUE);axes[1].axhline(.75,color=GRAY,ls='--',label='Pure doublet: 0.75')
    axes[1].set(ylabel='⟨S²⟩',xlabel='Mass-weighted IRC coordinate (√amu bohr)');axes[1].legend(frameon=False)
    for ax in axes:ax.grid(axis='y',alpha=.15);ax.axvline(0,color=GRAY,ls=':',lw=.8)
    save(fig,'irc-bond-spinquality.png')
    states=load('analysis/thermochemistry-20261019/audit.json')['states']
    def ladder(ax,values,title,color,ylim):
        for i,v in enumerate(values):
            ax.plot([i-.21,i+.21],[v,v],color=color,lw=2.5)
            ax.text(i,v+.6,f'{v:+.2f}' if i else '0.00',ha='center')
            if i<2:ax.plot([i+.21,i+.79],[v,values[i+1]],color=color,ls='--',lw=1)
        ax.set(title=title,xticks=[0,1,2],xticklabels=['R','TS','P'],ylim=ylim,xlim=(-.5,2.5))
        ax.axhline(0,color=GRAY,lw=.7);ax.grid(axis='y',alpha=.15)
    fig,axes=plt.subplots(1,2,figsize=(10,4.8),sharey=True,layout='constrained')
    for ax,key,title,col in zip(axes,['E','G_corrected'],['Electronic energy','Gibbs free energy (298.15 K, 1 atm)'],[BLUE,PINK]):
        ladder(ax,[states[s]['relative_kcal_mol'][key] for s in ['R','TS','P']],title,col,(-5,11))
    axes[0].set_ylabel('Relative energy (kcal/mol)');save(fig,'energy-gibbs-ladders.png')
    methods=load('analysis/report-review-20261019/method-comparison-data.json')['rows']
    fig,axes=plt.subplots(2,2,figsize=(8.4,6.6),sharey=True,layout='constrained')
    for ax,r,col in zip(axes.ravel(),methods,[BLUE,GRAY,GOLD,PINK]):
        ladder(ax,[0,r['forward'],r['reaction']],r['method'],col,(-5,15.5))
    for ax in axes[:,0]:ax.set_ylabel('Relative electronic energy (kcal/mol)')
    save(fig,'electronic-method-comparison.png')
    fig,axes=plt.subplots(1,2,figsize=(10,4.6),layout='constrained')
    for ax,s,title in zip(axes,['TS','P'],['Forward activation','Reaction R → P']):
        d=states[s]['relative_kcal_mol'];E=d['E'];EZ=d['E_ZPE'];H=d['H'];G=d['G_corrected']
        changes=[E,EZ-E,H-EZ,G-H,G];current=0
        assert abs(sum(changes[:4])-G)<1e-10
        for i,v in enumerate(changes):
            start=0 if i in [0,4] else current;end=start+v
            ax.bar(i,v,bottom=start,color=GRAY if i in [0,4] else BLUE if v<0 else GOLD,width=.62)
            ax.text(i,max(start,end)+.24,f'{v:+.2f}',ha='center',fontsize=10)
            if i<3:ax.plot([i+.31,i+.69],[end,end],color=GRAY,ls=':',lw=.8)
            current=end
        ax.set(title=title,xticks=range(5),xticklabels=['ΔE','ΔZPE','ΔHₜₕ','−TΔS','ΔG'],ylim=(-4.3,1.2) if s=='P' else (-.5,10.2))
        ax.axhline(0,color=GRAY,lw=.8);ax.grid(axis='y',alpha=.15)
    axes[0].set_ylabel('Energy contribution (kcal/mol)');save(fig,'thermal-decomposition.png')
    ts=stationary_data(ROOT/'data/calculations/ts-opt/calc.fchk');ref=ts['xyz_A'];z=ts['Z'];centered=ref-ref.mean(axis=0)
    mode=ts['cartesian_modes'][:,0].reshape(11,3);mode/=np.max(np.linalg.norm(mode,axis=1))
    branches=[load(f'analysis/irc-20261019/{b}.json') for b in ['reverse','forward']]
    path=list(reversed(branches[0]['points'][1:]))+branches[1]['points']
    poses=[align(p['xyz_A'],ref) for p in path];assert len(poses)==56
    extent=max(np.ptp(np.concatenate(poses+[centered+.18*mode,centered-.18*mode]),axis=0))/2+.65
    limits=[(-extent,extent)]*3
    fig=plt.figure(figsize=(6.4,5));ax=fig.add_subplot(projection='3d')
    molecule(ax,centered,z,f'Opening mode: {ts["frequencies_cm_1"][0]:.2f} cm⁻¹',limits,mode)
    fig.text(.5,.04,'Mode sign and displayed amplitude are arbitrary',ha='center',fontsize=10)
    save(fig,'ts-opening-mode.png')
    fig=plt.figure(figsize=(12,4.5))
    for j,(idx,title) in enumerate([(0,'Ring-side IRC endpoint'),(17,'Transition state'),(55,'Open-chain IRC endpoint')]):
        molecule(fig.add_subplot(1,3,j+1,projection='3d'),poses[idx],z,title,limits)
    fig.subplots_adjust(left=0,right=1,wspace=0,bottom=.02,top=.9);save(fig,'irc-key-geometries.png')
    if animations:
        fs=[frame(x,z,f'IRC q = {p["coordinate_sqrt_amu_bohr"]:+.3f} √amu bohr',limits) for x,p in zip(poses,path)]
        fs[0].save(dest/'irc-opening.gif',save_all=True,append_images=fs[1:],duration=120,loop=0)
        fs=[frame(centered+.18*np.sin(t)*mode,z,f'Opening mode {ts["frequencies_cm_1"][0]:.2f} cm⁻¹',limits) for t in np.linspace(0,2*np.pi,32,endpoint=False)]
        fs[0].save(dest/'ts-opening-mode.gif',save_all=True,append_images=fs[1:],duration=80,loop=0)
