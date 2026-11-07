import numpy as np, pandas as pd, json, itertools
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager as fm
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Rectangle
from scipy.cluster import hierarchy as sch
from scipy.stats import binomtest

import os
for f in ['Carlito-Regular','Carlito-Bold','Carlito-Italic','Carlito-BoldItalic']:
    p=f'/usr/share/fonts/truetype/crosextra/{f}.ttf'
    if os.path.exists(p): fm.fontManager.addfont(p)
mpl.rcParams.update({'font.family':['Carlito','DejaVu Sans'],'font.size':8,'axes.labelsize':8.5,'axes.titlesize':9,
    'xtick.labelsize':7.5,'ytick.labelsize':7.5,'legend.fontsize':7.2,'axes.linewidth':0.7,
    'xtick.major.width':0.7,'ytick.major.width':0.7,'xtick.major.size':3,'ytick.major.size':3,
    'axes.spines.top':False,'axes.spines.right':False,'savefig.dpi':300,'figure.dpi':150,
    'pdf.fonttype':42,'ps.fonttype':42})

val  = pd.read_pickle('log2.pkl'); flag = pd.read_pickle('flag.pkl'); clip = pd.read_pickle('clip.pkl')
qc = pd.read_pickle('qc.pkl'); gs = pd.read_pickle('gs.pkl'); D = pd.read_pickle('delta.pkl'); Dc = pd.read_pickle('delta_cens.pkl')
sim = pd.read_pickle('simgrid.pkl'); res = json.load(open('res.json'))
conds = list(val.columns); muts = conds[1:]
CC = {'WT':'#222222','R175H':'#0072B2','L194F':'#56B4E9','R273H':'#D55E00','R280K':'#E69F00','E285K':'#CC79A7'}
stars = ['DLG2AS (family)','lincRNA-p21','lincRNA-VLDLR','SNHG4']
TC = {'DLG2AS (family)':'#009E73','lincRNA-p21':'#0072B2','lincRNA-VLDLR':'#D55E00','SNHG4':'#CC79A7'}
downs = ['DISC2 (family)','H19 antisense','mascRNA']
UPC, DNC = '#B2182B', '#2166AC'
W = 7.2
def pl(ax,t,dx=-0.10,dy=1.04):
    ax.text(dx,dy,t,transform=ax.transAxes,fontsize=11,fontweight='bold',va='bottom',ha='left')
def save(fig,name):
    fig.savefig(f'../figures/{name}.png',dpi=300,bbox_inches='tight',facecolor='white')
    fig.savefig(f'../figures/{name}.pdf',bbox_inches='tight',facecolor='white'); plt.close(fig)

# ================= FIG 1 workflow + QC =================
fig = plt.figure(figsize=(W,4.4))
gs0 = fig.add_gridspec(2,2,height_ratios=[0.62,1.35],hspace=0.62,wspace=0.34)
ax = fig.add_subplot(gs0[0,:]); ax.set_xlim(0,100); ax.set_ylim(0,20); ax.axis('off'); pl(ax,'A',-0.01,1.0)
boxes = [("HEK293T cells","SV40 large-T antigen;\ntransient transfection",'#EAF2F8'),
         ("7 constructs","vector, WT p53, R175H,\nL194F, R273H, R280K,\nE285K",'#EAF2F8'),
         ("lncRNA qPCR array","SBI LncProfiler\nRA900A-1; 90 human\nlncRNAs",'#FDF2E3'),
         ("ΔΔCt vs 18S rRNA","control = 1.00; one value\nper target per condition\n(no replicates)",'#FDF2E3'),
         ("Deposited output","33 targets with fold\nchange (33/90 = 36.7%)",'#EAF6F0')]
x=0.5; w=18.2; gap=2.1
for i,(t,s_,c) in enumerate(boxes):
    ax.add_patch(FancyBboxPatch((x,1.2),w,17.6,boxstyle='round,pad=0.15,rounding_size=0.9',fc=c,ec='#555',lw=0.7))
    ax.text(x+w/2,15.6,t,ha='center',va='center',fontweight='bold',fontsize=7.6)
    ax.text(x+w/2,8.0,s_,ha='center',va='center',fontsize=6.4,linespacing=1.3)
    if i<len(boxes)-1:
        ax.add_patch(FancyArrowPatch((x+w+0.15,10),(x+w+gap-0.15,10),arrowstyle='-|>',mutation_scale=7,lw=0.9,color='#444'))
    x+=w+gap
# panel B digitization QC
axb = fig.add_subplot(gs0[1,0]); pl(axb,'B',-0.2)
lim=[0.8,4.2]
axb.plot(lim,lim,color='#999',lw=0.8,ls='--',zorder=1)
for t,g in qc.groupby('target'):
    axb.scatter(g.exact,g.digitized,s=22,color=TC[t],edgecolor='white',lw=0.4,label=t,zorder=3)
axb.set_xlim(lim); axb.set_ylim(lim)
axb.set_xlabel('Exact value in report (log$_2$ FC)'); axb.set_ylabel('Recovered from heatmap (log$_2$ FC)')
axb.text(0.04,0.95,f'n = {len(qc)} uncensored cells\nMAE = {qc.err.abs().mean():.3f} log$_2$\nmax |error| = {qc.err.abs().max():.3f} log$_2$\nr = {np.corrcoef(qc.exact,qc.digitized)[0,1]:.4f}',transform=axb.transAxes,va='top',fontsize=7)
axb.legend(loc='lower right',frameon=False,handletextpad=0.2,fontsize=6.6)
# panel C counts
axc = fig.add_subplot(gs0[1,1]); pl(axc,'C',-0.2)
rep = np.array([16,12,12,17,22,27]); recov=(clip>0).sum().values
xx=np.arange(6); bw=0.38
axc.bar(xx-bw/2,rep,bw,color='#BBBBBB',label='Reported in source report',edgecolor='#555',lw=0.5)
axc.bar(xx+bw/2,recov,bw,color='#0072B2',label='Recovered from reconstructed matrix',edgecolor='#555',lw=0.5)
for i in range(6): axc.text(xx[i]+bw/2,recov[i]+0.5,str(recov[i]),ha='center',fontsize=6.8)
axc.set_xticks(xx); axc.set_xticklabels(conds); axc.set_ylabel('Targets above 1.0 (of 33)'); axc.set_ylim(0,34)
axc.legend(frameon=False,loc='upper left',fontsize=6.8)
save(fig,'Fig1_workflow_QC')

# ================= FIG 2 heatmap =================
Z = sch.linkage(clip.values,'average')
order = sch.dendrogram(Z,no_plot=True)['leaves']
rnames = list(clip.index[order]); H = val.loc[rnames]; Fm = flag.loc[rnames]
fig = plt.figure(figsize=(W,8.6))
g = fig.add_gridspec(1,2,width_ratios=[0.9,6],wspace=0.02,left=0.02,right=0.80,top=0.935,bottom=0.115)
axd = fig.add_subplot(g[0]); axh = fig.add_subplot(g[1])
dd = sch.dendrogram(Z,orientation='left',no_plot=False,ax=axd,color_threshold=0,above_threshold_color='#666',link_color_func=lambda k:'#666')
axd.invert_yaxis(); axd.axis('off')
im = axh.imshow(H.clip(-5,5).values,aspect='auto',cmap='RdBu_r',vmin=-5,vmax=5,interpolation='nearest')
axh.set_xticks(range(6)); axh.set_xticklabels(conds,fontsize=8.5,fontweight='bold'); axh.xaxis.tick_top()
axh.set_yticks(range(len(rnames))); axh.set_yticklabels(rnames,fontsize=7.8); axh.yaxis.tick_right()
for lab in axh.get_yticklabels():
    t=lab.get_text()
    if t in stars: lab.set_color(UPC); lab.set_fontweight('bold')
    elif t in downs: lab.set_color(DNC); lab.set_fontweight('bold')
for s_ in axh.spines.values(): s_.set_visible(False)
axh.tick_params(length=0)
for i,r in enumerate(rnames):
    for j,c in enumerate(conds):
        if Fm.loc[r,c]=='censored' or (Fm.loc[r,c]=='exact' and abs(H.loc[r,c])>=5):
            axh.plot(j,i,'o',ms=2.4,mfc='white',mec='black',mew=0.35)
axh.set_xticks(np.arange(-.5,6,1),minor=True); axh.set_yticks(np.arange(-.5,len(rnames),1),minor=True)
axh.grid(which='minor',color='white',lw=0.6); axh.tick_params(which='minor',length=0)
axcb = fig.add_axes([0.30,0.088,0.36,0.012])
cb = fig.colorbar(im,cax=axcb,orientation='horizontal'); cb.set_label('log$_2$ fold change vs vector control',fontsize=8); cb.outline.set_linewidth(0.5)
cb.set_ticks([-5,-4,-2,0,2,4,5]); cb.set_ticklabels(['≤ −5','−4','−2','0','2','4','≥ 5'])
fig.text(0.02,0.020,'○ = censored at colour-scale limit (|log$_2$ FC| ≥ 5, i.e. ≥ 32-fold); true magnitude unknown from the heatmap alone.\nRed labels: ≥ 2-fold up in all six conditions. Blue labels: ≥ 2-fold down in all six conditions. Rows ordered by average-linkage clustering.',fontsize=7,va='center')
save(fig,'Fig2_heatmap_33targets')

# ================= FIG 3 four lncRNAs =================
fig,(a,b) = plt.subplots(1,2,figsize=(W,3.2),gridspec_kw={'wspace':0.42,'width_ratios':[1.05,1]})
pl(a,'A',-0.14); pl(b,'B',-0.16)
for t in stars:
    a.plot(range(6),val.loc[t].values,'-o',color=TC[t],lw=1.3,ms=4.2,label=t.replace(' (family)',''),mec='white',mew=0.5)
a.axhline(1,color='#888',lw=0.7,ls=':'); a.text(5.25,1.08,'2×',fontsize=6.8,color='#666',ha='right')
a.set_xticks(range(6)); a.set_xticklabels(conds); a.set_ylabel('log$_2$ fold change vs vector control'); a.set_ylim(0,10.3); a.set_xlim(-0.3,5.3)
a2=a.twinx(); a2.set_ylim(a.get_ylim()); tk=[0,2,4,6,8,10]; a2.set_yticks(tk); a2.set_yticklabels([f'{2**t:g}×' for t in tk]); a2.spines['right'].set_visible(True); a2.spines['right'].set_linewidth(0.7)
a.legend(frameon=False,loc='upper right',ncol=1,fontsize=7,handlelength=1.4)
a.annotate('594×',(0,val.loc['SNHG4','WT']),xytext=(0.35,9.75),fontsize=7,arrowprops=dict(arrowstyle='-',lw=0.5,color='#666'))
bw=0.16
for k,m in enumerate(muts):
    b.bar(np.arange(4)+(k-2)*bw,[D.loc[t,m] for t in stars],bw,color=CC[m],label=m,edgecolor='white',lw=0.3)
b.axhline(0,color='black',lw=0.7)
b.set_xticks(range(4)); b.set_xticklabels([t.replace(' (family)','').replace('lincRNA-','lincRNA-\n') for t in stars],fontsize=7.3)
b.set_ylabel('Mutant minus WT p53 (log$_2$ FC)'); b.legend(frameon=False,ncol=3,fontsize=6.8,loc='upper left',columnspacing=0.8,handlelength=1)
b.set_ylim(-4.2,4.9)
b.text(0.99,0.985,'Control-independent contrast:\nthe shared vector value cancels',transform=b.transAxes,ha='right',va='top',fontsize=6.6,color='#555',style='italic')
save(fig,'Fig3_four_lncRNAs')

# ================= FIG 4 robustness =================
fig,axs = plt.subplots(1,3,figsize=(W,2.75),gridspec_kw={'wspace':0.62,'width_ratios':[1.1,1.25,0.9]})
a,b,c = axs; pl(a,'A',-0.2); pl(b,'B',-0.2); pl(c,'C',-0.3)
Ts=[1.0,1.25,1.5,2,3,4]; med=clip.median(axis=0)
def counts(m):
    u=[];d=[]
    for T in Ts:
        l=np.log2(T)
        u.append(int(((m>l).all(axis=1) if T==1 else (m>=l-1e-9).all(axis=1)).sum()))
        d.append(int(((m<-l).all(axis=1) if T==1 else (m<=-l+1e-9).all(axis=1)).sum()))
    return u,d
u0,d0 = counts(val); u1,d1 = counts(val.sub(med,axis=1))
xs=np.arange(len(Ts))
a.plot(xs,u0,'-o',color=UPC,lw=1.4,ms=4,label='Up, all six')
a.plot(xs,d0,'-s',color=DNC,lw=1.4,ms=4,label='Down, all six')
a.plot(xs,u1,'--o',color=UPC,lw=1,ms=3,mfc='white',label='Up, centred*')
a.plot(xs,d1,'--s',color=DNC,lw=1,ms=3,mfc='white',label='Down, centred*')
a.set_xticks(xs); a.set_xticklabels(['>1','≥1.25','≥1.5','≥2','≥3','≥4']); a.set_xlabel('Fold threshold in every condition')
a.set_ylabel('Targets meeting criterion'); a.set_ylim(0,8.2); a.legend(frameon=False,fontsize=6.6,loc='upper right',handlelength=1.8)
a.text(0.0,-0.30,'*after subtracting each condition\'s median log$_2$ FC',transform=a.transAxes,fontsize=6.3,color='#555',style='italic')
a.axvline(3,color='#999',lw=0.6,ls=':')
yy=np.arange(6)
b.errorbar(gs.frac,yy,xerr=[gs.frac-gs.lo,gs.hi-gs.frac],fmt='none',ecolor='#555',lw=0.9,capsize=2)
b.scatter(gs.frac,yy,c=[CC[x] for x in gs.cond],s=26,zorder=3,edgecolor='white',lw=0.4)
b.axvline(0.5,color='#888',lw=0.7,ls='--'); b.set_yticks(yy); b.set_yticklabels(gs.cond); b.invert_yaxis(); b.set_xlim(0.1,1.62)
b.set_xticks([0.25,0.5,0.75,1.0]); b.spines['bottom'].set_bounds(0.1,1.0)
b.set_xlabel('Fraction of 33 targets above 1.0',x=0.35)
b.text(1.12,-0.85,'up',fontsize=6.4,ha='center',color='#555'); b.text(1.45,-0.85,'sign-test p',fontsize=6.4,ha='center',color='#555')
for i_,r in gs.iterrows():
    b.text(1.12,i_,f'{int(r.up)}/33',va='center',ha='center',fontsize=6.8)
    b.text(1.45,i_,('<0.001' if r.p<0.001 else f'{r.p:.2f}'),va='center',ha='center',fontsize=6.8,fontweight='bold' if r.p<0.05 else 'normal')
c.bar(range(6),gs.median_log2,color=[CC[x] for x in gs.cond],edgecolor='white',lw=0.4)
c.axhline(0,color='black',lw=0.7); c.set_xticks(range(6)); c.set_xticklabels(conds,rotation=45,ha='right')
c.set_ylabel('Median log$_2$ FC across 33 targets')
c2=c.twinx(); c2.set_ylim(c.get_ylim()); tk=[-1,0,1,2,3]; c2.set_yticks(tk); c2.set_yticklabels([f'{2**t:g}×' for t in tk]); c2.spines['right'].set_visible(True); c2.spines['right'].set_linewidth(0.7)
save(fig,'Fig4_robustness_global_shift')

# ================= FIG 5 shared-control =================
fig,axs = plt.subplots(1,3,figsize=(W,3.0),gridspec_kw={'wspace':0.6,'width_ratios':[1.1,1,1]})
a,b,c = axs; pl(a,'A',-0.2); pl(b,'B',-0.24); pl(c,'C',-0.24)
rm = clip.median(axis=1).sort_values()
cols=[UPC if t in stars else (DNC if t in downs else '#9a9a9a') for t in rm.index]
a.bar(range(33),rm.values,color=cols,width=0.8)
a.axhspan(-1,1,color='#EEEEEE',zorder=0); a.axhline(0,color='black',lw=0.6)
a.set_xticks([]); a.set_xlabel('33 targets, ranked'); a.set_ylabel('Median log$_2$ FC across six conditions')
a.set_ylim(-5.6,5.6)
a.text(31.5,-0.8,'±2-fold',fontsize=6.4,color='#777',ha='right',va='top')
a.text(1,4.9,'3 down in all six',color=DNC,fontsize=7,fontweight='bold',va='top'); a.text(26,5.45,'4 up in all six',color=UPC,fontsize=7,fontweight='bold',ha='center',va='top')
a.text(0.5,-0.13,'Comparable numbers of consistently up (4) and\nconsistently down (3) targets fit a shared control error',transform=a.transAxes,ha='center',va='top',fontsize=6.4,style='italic',color='#555')
for st,ls in zip([0.3,0.6,1.0],['-','--',':']):
    s=sim[sim.s_t==st]; b.plot(s.s_c,(s.E_up+s.E_down)/2,ls,color='#444',lw=1.3,label=f'σ$_t$ = {st}')
b.axhline(4,color=UPC,lw=1); b.axhline(3,color=DNC,lw=1)
b.text(3.0,4.2,'observed up (4)',color=UPC,fontsize=6.6,ha='right'); b.text(3.0,2.5,'observed down (3)',color=DNC,fontsize=6.6,ha='right')
b.set_xlabel('σ$_c$: error in control value (Ct cycles)'); b.set_ylabel('Expected targets meeting criterion\n(≥ 2-fold in all six), no true effect')
b.legend(frameon=False,fontsize=6.6,loc='upper left',title='per-well noise',title_fontsize=6.6); b.set_xlim(0,3.05); b.set_ylim(0,11)
ct=np.linspace(15,36,200)
for lod,ls,col in [(35,'-','#0072B2'),(37,'--','#0072B2'),(40,':','#0072B2')]:
    c.plot(ct,2.0**(lod-ct),ls,color=col,lw=1.3,label=f'LOD Ct = {lod}')
ex=[('SNHG4 WT',594.01,UPC),('Emx2os E285K',1007.95,'#7B3294'),('Zeb2NAT R273H',153.06,'#E08214')]
for n,f,col in ex:
    c.axhline(f,color=col,lw=0.6,ls='--',alpha=0.8); x35=35-np.log2(f)
    c.plot([x35],[f],'o',ms=4.2,color=col,zorder=4)
yt=6.5
c.text(15.4,yt*1.9,'Sample Ct needed (control at LOD 35)',fontsize=6.3,color='#333',fontweight='bold',va='bottom')
for k,(n,f,col) in enumerate(ex):
    c.text(15.4,yt/(1.75**k),f'{n}: Ct ≤ {35-np.log2(f):.1f}',fontsize=6.4,color=col,va='bottom')
c.set_yscale('log'); c.set_ylim(1,1e5); c.set_xlim(15,36)
c.set_xlabel('Target Ct in p53-expressing sample'); c.set_ylabel('Maximum measurable fold change')
c.legend(frameon=False,fontsize=6.4,loc='upper right',bbox_to_anchor=(1.0,1.0))
save(fig,'Fig5_shared_control_and_detection_limit')

# ================= FIG 6 control-free mutant vs WT =================
Dm = D.copy(); order6 = Dm.mean(axis=1).sort_values(ascending=False).index
Dm = Dm.loc[order6]; Cm = Dc.loc[order6]
cand_up = [t for t in Dm.index if (Dm.loc[t]>=2).sum()>=3]; cand_dn=[t for t in Dm.index if (Dm.loc[t]<=-2).sum()>=3]
json.dump({'cand_up':cand_up,'cand_dn':cand_dn},open('cands.json','w'))
fig = plt.figure(figsize=(5.8,8.4))
g = fig.add_gridspec(1,1,left=0.02,right=0.78,top=0.935,bottom=0.15)
ax=fig.add_subplot(g[0])
cmap = plt.get_cmap('PuOr_r').copy(); cmap.set_bad('#D9D9D9')
im=ax.imshow(np.ma.masked_invalid(Dm.clip(-6,6).values),aspect='auto',cmap=cmap,vmin=-6,vmax=6,interpolation='nearest')
ax.set_xticks(range(5)); ax.set_xticklabels(muts,fontsize=8.5,fontweight='bold'); ax.xaxis.tick_top()
ax.set_yticks(range(33)); ax.set_yticklabels(Dm.index,fontsize=7.6); ax.yaxis.tick_right(); ax.tick_params(length=0)
for s_ in ax.spines.values(): s_.set_visible(False)
for lab in ax.get_yticklabels():
    t=lab.get_text()
    if t in cand_up: lab.set_color('#8C510A'); lab.set_fontweight('bold')
    if t in cand_dn: lab.set_color('#542788'); lab.set_fontweight('bold')
for i,r in enumerate(Dm.index):
    for j,m in enumerate(muts):
        if np.isnan(Dm.loc[r,m]): ax.text(j,i,'n.e.',ha='center',va='center',fontsize=5.8,color='#555')
        elif Cm.loc[r,m] or abs(Dm.loc[r,m])>6: ax.plot(j,i,'o',ms=2.3,mfc='white',mec='black',mew=0.35)
ax.set_xticks(np.arange(-.5,5,1),minor=True); ax.set_yticks(np.arange(-.5,33,1),minor=True); ax.grid(which='minor',color='white',lw=0.6); ax.tick_params(which='minor',length=0)
axcb = fig.add_axes([0.20,0.108,0.40,0.011])
cb=fig.colorbar(im,cax=axcb,orientation='horizontal'); cb.set_label('log$_2$ (mutant FC / WT FC)',fontsize=8); cb.outline.set_linewidth(0.5)
fig.text(0.02,0.030,'Every condition is divided by the same vector value, so mutant/WT ratios are free of that control.\n○ = one term is censored, so the cell is a bound rather than an estimate (clipped at ±6). n.e. = not estimable (both terms censored).\nBrown labels: ≥ 4-fold above WT in ≥ 3 of 5 mutants. Purple labels: ≥ 4-fold below WT in ≥ 3 of 5.',fontsize=6.7,va='center')
save(fig,'Fig6_mutant_vs_WT_control_free')

# ================= FIG 7 condition similarity =================
raw = clip.corr(method='spearman'); resid = clip.sub(clip.median(axis=1),axis=0).corr(method='spearman')
ordc = ['WT','R175H','L194F','E285K','R273H','R280K']
fig,(a,b) = plt.subplots(1,2,figsize=(W,3.5),gridspec_kw={'wspace':0.38})
pl(a,'A',-0.16); pl(b,'B',-0.16)
for ax,m,ttl,vm,cm in [(a,raw,'Raw fold changes\n(shared-control effect included)',(0,1),'viridis'),(b,resid,'After removing each target\'s median\n(condition-specific signal only)',(-0.6,0.6),'RdBu_r')]:
    mm = m.loc[ordc,ordc]
    im=ax.imshow(mm.values,cmap=cm,vmin=vm[0],vmax=vm[1])
    ax.set_xticks(range(6)); ax.set_xticklabels(ordc,rotation=45,ha='right'); ax.set_yticks(range(6)); ax.set_yticklabels(ordc)
    for i in range(6):
        for j in range(6):
            v=mm.values[i,j]; dark = (v>0.55 if cm=='viridis' else abs(v)>0.4)
            ax.text(j,i,f'{v:.2f}',ha='center',va='center',fontsize=6.8,color='white' if (cm=='viridis' and v<0.6) or (cm!='viridis' and abs(v)>0.4) else 'black')
    ax.set_title(ttl,fontsize=8,pad=5); ax.tick_params(length=0)
    for s in ax.spines.values(): s.set_visible(False)
    cb=fig.colorbar(im,ax=ax,fraction=0.046,pad=0.03); cb.set_label('Spearman ρ',fontsize=7.5); cb.outline.set_linewidth(0.4)
b.add_patch(Rectangle((3.5,3.5),2,2,fill=False,ec='black',lw=1.4))
b.text(0.5,-0.30,'Box: {R273H, R280K} is the highest-contrast 2-vs-4 split of the six conditions\n(rank 1 of 15 possible splits; exact permutation p = 1/15 = 0.067, exploratory).\nCentring imposes a baseline of about −0.2, so values near −0.2 mean “unrelated”.',transform=b.transAxes,ha='center',va='top',fontsize=6.5,style='italic',color='#444')
save(fig,'Fig7_condition_similarity')
print('done'); print(cand_up); print(cand_dn)
