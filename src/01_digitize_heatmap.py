import numpy as np, json
from PIL import Image
im = np.array(Image.open('../../data/raw/heatmap_original_from_report.png').convert('RGB')).astype(float)
# LUT from colorbar
y0,y1 = 68,1210
cx = 1574
ys = np.arange(y0,y1+1)
lut_col = im[ys, cx-10:cx+10, :].mean(axis=1)
lut_val = 5 - (ys - y0)/(y1-y0)*10
# grid
x0,x1 = 260,1468
n_r, n_c = 33, 7
cw = (x1-x0+1)/n_c; ch = (y1-y0+1)/n_r
rows = ["21A","Air","Alpha 280","ANRIL","anti-NOS2A","BC200","DISC2 (family)","DLG2AS (family)","E2F4 antisense","Emx2os","Evf1 and EVF2","GAS5-family","H19 antisense","HOTAIR","HOTTIP","HULC","L1PA16","lincRNA-p21","lincRNA-RoR","lincRNA-SFMBT2","lincRNA-VLDLR","LUST","Malat1","mascRNA","MEG9","MER11C","PCGEM1","PR antisense transcripts","snaR","SNHG4","Tmevpg1","Zeb2NAT","Zfhx2as"]
cols = ["pcDNA6A-Flag","Wtp53","R175H","L194F","R273H","R280K","E285K"]
M = np.zeros((n_r,n_c)); D = np.zeros((n_r,n_c))
for i in range(n_r):
    for j in range(n_c):
        cxp = int(x0 + (j+0.5)*cw); cyp = int(y0 + (i+0.5)*ch)
        patch = im[cyp-6:cyp+7, cxp-30:cxp+31, :].reshape(-1,3)
        c = np.median(patch, axis=0)
        d = np.linalg.norm(lut_col - c, axis=1)
        k = d.argmin()
        M[i,j] = lut_val[k]; D[i,j] = d[k]
print('max colour distance', D.max())
np.save('log2_digitized.npy', M)
json.dump({'rows':rows,'cols':cols}, open('labels.json','w'))
np.set_printoptions(precision=2, suppress=True, linewidth=200)
for r,row in zip(rows,M): print(f'{r:26s}', row)
# validation
up = (M[:,1:]>0).sum(axis=0); print('up counts (log2>0):', up, ' expected [16 12 12 17 22 27]')
known = {'DLG2AS (family)':[3.64,7.99,7.66,30.80,49.19,41.37],'lincRNA-p21':[10.34,3.01,3.77,4.52,6.01,14.38],'lincRNA-VLDLR':[8.68,3.80,2.40,6.07,7.86,46.19],'SNHG4':[594.01,224.10,139.25,62.06,59.76,95.62]}
for k,v in known.items():
    i = rows.index(k); print(k, 'digitized log2', M[i,1:], ' true log2', np.round(np.log2(v),2))
