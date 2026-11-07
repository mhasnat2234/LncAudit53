import numpy as np, pandas as pd, json
from scipy import stats
rows = json.load(open('labels.json'))['rows']
conds = ['WT','R175H','L194F','R273H','R280K','E285K']
M = np.load('log2_digitized.npy')[:,1:]          # drop control column
df = pd.DataFrame(M, index=rows, columns=conds)

# ---- exact overrides from the report tables ----
exact = {'DLG2AS (family)':[3.64,7.99,7.66,30.80,49.19,41.37],
         'lincRNA-p21':[10.34,3.01,3.77,4.52,6.01,14.38],
         'lincRNA-VLDLR':[8.68,3.80,2.40,6.07,7.86,46.19],
         'SNHG4':[594.01,224.10,139.25,62.06,59.76,95.62]}
flag = pd.DataFrame('digitized', index=rows, columns=conds)
cens = (df.abs() >= 4.9)
flag[cens] = 'censored'
dig = df.copy()
for k,v in exact.items():
    df.loc[k] = np.log2(v); flag.loc[k] = 'exact'
df.loc['Emx2os','E285K'] = np.log2(1007.95); flag.loc['Emx2os','E285K']='exact'
df.loc['Zeb2NAT','R273H'] = np.log2(153.06); flag.loc['Zeb2NAT','R273H']='exact'
# clipped matrix for analyses (censored values stay at +/-5 = lower bound on |effect|)
clip = df.clip(-5,5)
# QC: digitized vs exact on uncensored exact cells
qc=[]
for k,v in exact.items():
    for c,x in zip(conds,np.log2(v)):
        if abs(dig.loc[k,c])<4.9: qc.append((k,c,x,dig.loc[k,c]))
qc=pd.DataFrame(qc,columns=['target','cond','exact','digitized'])
qc['err']=qc.digitized-qc.exact
print('QC n=%d MAE=%.3f maxabs=%.3f bias=%+.3f r=%.4f'%(len(qc),qc.err.abs().mean(),qc.err.abs().max(),qc.err.mean(),np.corrcoef(qc.exact,qc.digitized)[0,1]))
print('control col digitised mean (true 0):', np.load('log2_digitized.npy')[:,0].mean())
print('censored cells:', int(cens.sum().sum()), 'of', cens.size)

# ---- 1. up counts ----
up = (clip>0).sum(); print('up',up.values, 'down', (clip<0).sum().values)
# sign tests (global shift diagnostic)
from scipy.stats import binomtest
res=[]
for c in conds:
    k=int((clip[c]>0).sum()); n=33
    bt=binomtest(k,n,0.5); ci=bt.proportion_ci(confidence_level=0.95,method='exact')
    res.append((c,k,k/n,ci.low,ci.high,bt.pvalue,clip[c].median()))
gs=pd.DataFrame(res,columns=['cond','up','frac','lo','hi','p','median_log2'])
print(gs.round(4))
# ---- 2. threshold sensitivity ----
ths=[1.0,1.25,1.5,2,3,4]
tab=[]
for T in ths:
    l=np.log2(T)
    upall=(clip>=l).all(axis=1) if T>1 else (clip>0).all(axis=1)
    dnall=(clip<=-l).all(axis=1) if T>1 else (clip<0).all(axis=1)
    tab.append((T,int(upall.sum()),int(dnall.sum()),list(clip.index[upall]),list(clip.index[dnall])))
for t in tab: print(t)
# ---- 3. shared control diagnostics ----
rowmed=clip.median(axis=1); rowsd=clip.std(axis=1,ddof=1)
tot=((clip-clip.values.mean())**2).values.sum()
rowmean=clip.mean(axis=1)
between=((rowmean-clip.values.mean())**2*6).sum()
print('variance share of target (row) effect: %.3f'%(between/tot))
# ---- 4. control-free: mutant minus WT ----
D = clip[['R175H','L194F','R273H','R280K','E285K']].sub(clip['WT'],axis=0)
print(D.round(2).sort_values('E285K'))
# ---- 5. condition similarity ----
rho_raw=clip.corr(method='spearman')
cen=clip.sub(clip.median(axis=1),axis=0)
rho_res=cen.corr(method='spearman')
print(rho_raw.round(2)); print(rho_res.round(2))
df.to_pickle('log2.pkl'); flag.to_pickle('flag.pkl'); clip.to_pickle('clip.pkl'); qc.to_pickle('qc.pkl'); gs.to_pickle('gs.pkl')
