import numpy as np, pandas as pd, itertools, json
from scipy import stats
val  = pd.read_pickle('log2.pkl')      # exact where known, else digitized (saturating at +/-5)
flag = pd.read_pickle('flag.pkl')
clip = pd.read_pickle('clip.pkl')
conds = list(val.columns)
out = {}

# ---- global-median re-centering sensitivity ----
med = clip.median(axis=0)
cen = clip.sub(med, axis=1)       # global-median normalised
stars = ['DLG2AS (family)','lincRNA-p21','lincRNA-VLDLR','SNHG4']
vs = val.loc[stars].sub(med,axis=1)  # exact values, re-centred
out['median_shift_log2'] = med.round(3).to_dict()
out['stars_min_after_centering_FC'] = float(2**vs.min().min())
out['stars_min_cell'] = [stars[int(np.argmin(vs.min(axis=1)))], conds[int(np.argmin(vs.min(axis=0)))]]
print('stars min FC after centering', 2**vs.min().min())
print((2**vs).round(2))
# up in all six after centering at >=2x (use clip for non-exact)
cv = val.copy(); cv = cv.sub(med,axis=1)
for T in [1.0,2.0]:
    l=np.log2(T)
    up=(cv>l).all(axis=1) if T==1 else (cv>=l).all(axis=1)
    dn=(cv<-l).all(axis=1) if T==1 else (cv<=-l).all(axis=1)
    print('after centering T=%s up-all-six:'%T, list(cv.index[up]), ' down-all-six:', list(cv.index[dn]))
# up counts after centering
print('up counts after centering', (cv>0).sum().values)

# ---- control-free mutant-vs-WT (use exact values; censored remain bounds) ----
muts = conds[1:]
D = val[muts].sub(val['WT'],axis=0)
bothc = flag[muts].eq('censored') & flag[['WT']].eq('censored').values & (np.sign(val[muts]).values==np.sign(val[['WT']]).values)
D = D.mask(bothc)
D.to_pickle('delta.pkl')
censD = (flag[muts].eq('censored')) | (flag[['WT']].eq('censored').values)
censD.to_pickle('delta_cens.pkl')
print(D.loc[stars].round(2))

# ---- block-structure test on residual correlations ----
res = clip.sub(clip.median(axis=1),axis=0)
rho = res.corr(method='spearman')
def contrast(B):
    A=[c for c in conds if c not in B]
    within=[rho.loc[a,b] for a,b in itertools.combinations(A,2)]+[rho.loc[a,b] for a,b in itertools.combinations(B,2)]
    between=[rho.loc[a,b] for a in A for b in B]
    return np.mean(within)-np.mean(between)
splits=[(B,contrast(list(B))) for B in itertools.combinations(conds,2)]
splits.sort(key=lambda x:-x[1])
for B,s in splits[:5]: print(B, round(s,3))
obs = contrast(['R273H','R280K'])
rank = 1+sum(s>obs+1e-12 for B,s in splits)
print('observed {R273H,R280K} rank', rank,'/15 contrast',round(obs,3))
out['block_rank']=rank; out['block_contrast']=float(obs)
# ---- within-target spread -> anchor for treated-well noise ----
unc = ~(flag.eq('censored')).any(axis=1)
rsd = clip.loc[unc].std(axis=1,ddof=1)
print('uncensored rows',int(unc.sum()),'median within-target SD (log2)',rsd.median().round(2))
out['median_row_sd']=float(rsd.median())

# ---- simulation: no true effect, shared control error ----
rng=np.random.default_rng(7)
def sim(sc, st, nsim=4000, T=1.0, ntar=33, k=6):
    nu=nd=0
    d=rng.normal(0,sc,(nsim,ntar,1)); e=rng.normal(0,st,(nsim,ntar,k))
    x=-d+e
    up=(x>=T).all(axis=2).sum(axis=1); dn=(x<=-T).all(axis=2).sum(axis=1)
    return up.mean(), dn.mean(), (up>=4).mean(), ((up>=4)&(dn>=3)).mean()
grid=[]
for st in [0.3,0.6,1.0]:
    for sc in [0.25,0.5,1,1.5,2,2.5,3]:
        u,dn,p4,p43=sim(sc,st); grid.append((st,sc,u,dn,p4,p43))
g=pd.DataFrame(grid,columns=['s_t','s_c','E_up','E_down','P_up>=4','P_up>=4&down>=3'])
print(g.round(3).to_string())
g.to_pickle('simgrid.pkl')
json.dump(out,open('res.json','w'),indent=1,default=lambda o: o.item() if hasattr(o,'item') else str(o))

# ---- write processed tables ----
rows=[]
for t in val.index:
    for c in val.columns:
        rows.append(dict(target=t,condition=c,log2FC=round(val.loc[t,c],3),fold_change=round(2**val.loc[t,c],3),source=flag.loc[t,c],note='bound (|log2FC|>=5)' if flag.loc[t,c]=='censored' else ''))
pd.DataFrame(rows).to_csv('../../data/processed/Supplementary_Table_S1_reconstructed_matrix.csv',index=False)
D.round(3).to_csv('../../data/processed/mutant_minus_WT_log2.csv')
