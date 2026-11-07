from reportlab.lib.pagesizes import A4
from reportlab.lib.units import inch
from reportlab.lib import colors
from reportlab.lib.styles import ParagraphStyle
from reportlab.lib.enums import TA_JUSTIFY, TA_LEFT
from reportlab.platypus import *
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.pdfbase.pdfmetrics import registerFontFamily
from PIL import Image as PILImage
import os
d='/usr/share/fonts/truetype/crosextra/'
if os.path.exists(d+'Carlito-Regular.ttf'):
    fl=[('C','Carlito-Regular'),('C-B','Carlito-Bold'),('C-I','Carlito-Italic'),('C-BI','Carlito-BoldItalic')]
else:
    d='/usr/share/fonts/truetype/dejavu/'
    fl=[('C','DejaVuSans'),('C-B','DejaVuSans-Bold'),('C-I','DejaVuSans-Oblique'),('C-BI','DejaVuSans-BoldOblique')]
for n,f in fl:
    pdfmetrics.registerFont(TTFont(n,d+f+'.ttf'))
registerFontFamily('C',normal='C',bold='C-B',italic='C-I',boldItalic='C-BI')
NAVY=colors.HexColor('#1F3A5F'); GREY=colors.HexColor('#555555'); LIGHT=colors.HexColor('#EEF2F7')
body=ParagraphStyle('b',fontName='C',fontSize=10,leading=13.6,alignment=TA_JUSTIFY,spaceAfter=6)
small=ParagraphStyle('s',parent=body,fontSize=8.6,leading=11,alignment=TA_LEFT,textColor=GREY,spaceAfter=10)
h1=ParagraphStyle('h1',fontName='C-B',fontSize=14,leading=17,textColor=NAVY,spaceBefore=12,spaceAfter=6)
h2=ParagraphStyle('h2',fontName='C-B',fontSize=11,leading=14,textColor=NAVY,spaceBefore=8,spaceAfter=4)
title=ParagraphStyle('t',fontName='C-B',fontSize=21,leading=25,textColor=NAVY,spaceAfter=6)
sub=ParagraphStyle('sub',fontName='C',fontSize=11,leading=14,textColor=GREY,spaceAfter=10)
cell=ParagraphStyle('c',fontName='C',fontSize=8.2,leading=10.2)
cellb=ParagraphStyle('cb',parent=cell,fontName='C-B',textColor=colors.white)
ref=ParagraphStyle('r',fontName='C',fontSize=8.4,leading=10.6,leftIndent=14,firstLineIndent=-14,spaceAfter=2.5)
bul=ParagraphStyle('bul',parent=body,leftIndent=14,bulletIndent=2,spaceAfter=3)
FW=A4[0]-1.3*inch
def P(t,s=body): return Paragraph(t,s)
def fig(name,cap,maxh=8.2*inch):
    p=f'../figures/{name}.png'; w,h=PILImage.open(p).size
    sc=min(FW/w,maxh/h)
    return KeepTogether([Image(p,w*sc,h*sc),Spacer(1,3),P(cap,small)])
def tbl(data,widths,hdr=True):
    rows=[[P(c,cellb if (hdr and i==0) else cell) if isinstance(c,str) else c for c in r] for i,r in enumerate(data)]
    t=Table(rows,colWidths=widths,repeatRows=1)
    t.setStyle(TableStyle([('BACKGROUND',(0,0),(-1,0),NAVY),('VALIGN',(0,0),(-1,-1),'TOP'),
      ('ROWBACKGROUNDS',(0,1),(-1,-1),[colors.white,LIGHT]),('LINEBELOW',(0,0),(-1,-1),0.25,colors.HexColor('#C9D1DC')),
      ('TOPPADDING',(0,0),(-1,-1),3),('BOTTOMPADDING',(0,0),(-1,-1),3),('LEFTPADDING',(0,0),(-1,-1),4),('RIGHTPADDING',(0,0),(-1,-1),4)]))
    return t
def footer(c,doc):
    c.saveState(); c.setFont('C',8); c.setFillColor(GREY)
    c.drawString(0.65*inch,0.4*inch,'GSE159185: p53-dependent lncRNA qPCR array, secondary analysis'); c.drawRightString(A4[0]-0.65*inch,0.4*inch,f'{doc.page}'); c.restoreState()
S=[]
S+= [Spacer(1,20),P('Mutant and wild-type p53 and lncRNA expression in HEK293T cells',title),
 P('Secondary analysis of public qPCR array GSE159185: literature-gap review, robustness tests and control-independent re-analysis',sub),
 HRFlowable(width='100%',thickness=1.2,color=NAVY,spaceAfter=8)]
S+= [P('Abstract',h1),
P('GSE159185 holds a 90-assay lncRNA qPCR array (SBI LncProfiler) run on HEK293T cells transiently expressing wild-type p53 or one of five hotspot mutants (R175H, L194F, R273H, R280K, E285K), each compared with an empty-vector control [1,2]. The deposited data give one fold-change value per target per condition, with no replicates, for 33 of the 90 targets. The source study reported four lncRNAs (DLG2-AS1, VLDLR-AS1, lincRNA-p21, SNHG4) raised by every p53 construct [1]. We merged the two earlier report versions, reviewed the primary paper and the literature on each highlighted lncRNA, and identified gaps that can be addressed with the data in hand. Because the raw files could not be downloaded, we rebuilt the 33 × 6 fold-change matrix from the report tables and the original heatmap (mean error 0.03 log<sub>2</sub> units against 14 known values; the reported up/down counts were reproduced exactly).'),
P('Five findings follow. (i) The “four lncRNAs” result depends on an unstated threshold: with any fold above 1.0, six targets are up in all conditions, and with a 2-fold cut-off there are four. (ii) Three targets (DISC2, H19 antisense, mascRNA) are at least 5-fold down in all six conditions; this mirror group is not discussed. (iii) E285K shifts the median of all 33 targets up by 2.9 log<sub>2</sub> units (about 7-fold; 27 of 33 targets up, sign-test p = 0.0003), which points to a normalisation or input offset rather than target-specific biology. (iv) Because every condition shares one control value, agreement across conditions is not independent evidence; a control error of about 1.2–2 Ct cycles alone could produce equal numbers of consistently up and down targets. (v) Contrasts of each mutant against wild-type p53 are free of the control and show that SNHG4 induction is strongest with wild-type p53, whereas DLG2-AS1, Emx2os, Zeb2NAT, HOTTIP and others rise more with specific mutants. None of this replaces replicated validation; Section 6 gives a design for it.')]
S+= [P('1. Background and literature',h1),
P('<b>Source study.</b> He et al. screened p53-responsive lncRNAs in HEK293T cells to understand how ERα redirects mutant p53 between the <i>lincRNA-p21</i> and <i>DDB2</i> promoters in breast cancer [1]. The array was the SBI LncProfiler RA900A-1 (90 human lncRNAs). The paper states that p53 constructs affected about one-third (36%) of the lncRNAs and that four were raised by all of them. Only lincRNA-p21 was followed up: confirmed with two primer sets, dependent on endogenous R280K in MDA-MB-231 cells, and traced to a G-quadruplex rather than a canonical p53 response element in its promoter [1]. We note that 33/90 = 36.7%, close to the quoted 36%. That is consistent with the percentage describing targets that yielded fold-change values (a detection filter) rather than a count of regulated genes, but we could not confirm how the paper derived it.'),
P('<b>Mutant p53 and lncRNAs.</b> Reviews describe mutant p53 gain-of-function acting through new transcription-factor partnerships and non-canonical DNA binding, with lncRNAs as downstream effectors [3]. Examples include R273H-regulated lncRNAs in colorectal stem-like cells [4] and LINC01605 in breast and ovarian lines [5]. The source paper itself notes that lincRNA-p21 level is not tied to <i>TP53</i> status in colorectal or lung tumours [1].'),
P('<b>The highlighted lncRNAs.</b> SNHG4 is oncogenic or prognostic in many tumours, mostly through miRNA-sponge mechanisms [7,8]; one lung study links it to KDM3A, p53-K372 methylation and p21 [9]; a 2021 review found no breast-cancer study [8]. DLG2-AS1 was found on the same SBI 90-lncRNA panel to be lower in lung adenocarcinoma, without cis-regulation of DLG2 [6]. VLDLR-AS1 (linc-VLDLR) has been linked to extracellular-vesicle responses to chemotherapy [16] and is lower in colorectal tumours [15]. Emx2os is tumour-suppressive in lung, prostate and cervical settings [12–14] but raised in ovarian cancer. Zeb2NAT retains the <i>ZEB2</i> 5′-UTR intron and raises ZEB2 protein [10]; its knockdown reduces EMT features in MCF7 and MDA-MB-231 cells [11]. We found no study connecting any of these except lincRNA-p21 to mutant p53.'),
P('<b>Methodological context.</b> HEK293T cells carry SV40 large T antigen, which binds p53 and blocks its transactivation [17,18]; this complicates reading wild-type and mutant p53 effects. Single-reference ΔΔCt quantification [20], reporting standards [19], multi-gene reference normalisation [21] and global-median normalisation for large panels [22] are the relevant standards for the issues raised below.')]
S+= [P('2. Data and reconstruction',h1),
P('Only the two report PDFs were available; GEO could not be accessed from the analysis environment. The tables give exact fold changes for four lncRNAs (24 values) and the three largest values; the remaining cells come from the original heatmap, whose colour scale was read back through its own colour bar. Against the 14 uncensored known cells the recovery error was 0.034 log<sub>2</sub> units on average (maximum 0.073; r = 0.9996), and the number of targets above 1.0 matched the report in all six conditions (Figure 1). The heatmap saturates at |log<sub>2</sub> FC| = 5; 52 cells reached that limit, 12 of which are resolved exactly by the report tables, and the other 40 are treated as bounds. Values near zero are less certain (about ±0.1 log<sub>2</sub>). The full matrix, with a source flag for every cell, is Supplementary Table S1.'),
fig('Fig1_workflow_QC','<b>Figure 1. Experiment and data reconstruction.</b> (A) Design of GSE159185. (B) Heatmap-derived values against exact values from the report tables. (C) Targets above 1.0 per condition, reported versus recovered.',5.2*inch),
fig('Fig2_heatmap_33targets','<b>Figure 2. Fold changes of all 33 targets.</b> Reconstructed from the deposited values; every cell is relative to the empty-vector control (log<sub>2</sub> scale). Circles mark cells at the colour-scale limit.')]
S+= [P('3. Results',h1),
P('3.1 The four-lncRNA result and its mirror image',h2),
P('The four lncRNAs are 2.4–594-fold above control in all six conditions, SNHG4 the largest (594-fold in wild-type p53) and lincRNA-p21 the smallest (3.0-fold at lowest; Figure 3A). The claim that exactly four targets behave this way holds only at a 2-fold threshold. Requiring any increase (above 1.0) in every condition also admits lincRNA-SFMBT2 and MEG9, which sit at roughly 1.4–1.5-fold in their weakest conditions (Figure 4A). Three targets move the other way: DISC2, H19 antisense and mascRNA are at least 5-fold lower than control in all six conditions, a group as consistent as the four up-regulated ones and not mentioned in the source.'),
fig('Fig3_four_lncRNAs','<b>Figure 3. The four consistently induced lncRNAs.</b> (A) Fold change against vector control, exact values from the report (right axis shows fold change). (B) Each mutant relative to wild-type p53; this contrast does not depend on the control value. SNHG4 and, except for E285K, lincRNA-p21 are lower with mutants than with wild-type p53; DLG2-AS1 is higher, most with R273H, R280K and E285K.',3.4*inch),
P('3.2 A global offset in E285K and R280K',h2),
P('Across 33 targets, an unbiased experiment would be expected to show roughly half up and half down. Wild-type p53, R273H, R175H and L194F do (12–17 of 33 up). E285K does not: 27 of 33 are up (95% CI 65–93%; exact sign test p = 0.0003) and the median target is raised 7.2-fold (log<sub>2</sub> +2.85); R280K shows a smaller shift (22 of 33 up, median 2.7-fold, p = 0.08) (Figure 4B,C). A shift of this kind, common to unrelated targets, is what a change in the 18S reference or in RNA input would produce. The targets are not independent and the set is limited to detected targets, so this is a diagnostic, not a significance claim. Subtracting each condition’s median (global-median normalisation [22]) leaves the four lncRNAs at 2.0-fold or more everywhere, though lincRNA-p21 in E285K falls to exactly 1.99-fold and lincRNA-SFMBT2 enters the all-six set at the 2-fold threshold. The up/down result for E285K, including the report’s remark that it is the “most up” condition, should not be read as a biological feature without that correction.'),
fig('Fig4_robustness_global_shift','<b>Figure 4. Robustness of the headline result.</b> (A) Targets meeting an all-six-conditions criterion at increasing thresholds, before and after median-centring. (B) Fraction of targets above 1.0 per condition with exact 95% intervals and two-sided sign-test p. (C) Median log<sub>2</sub> fold change per condition.',3.0*inch),
P('3.3 One control value stands behind every fold change',h2),
P('Each condition is divided by the same empty-vector measurement, so an error in that single value pushes all six conditions the same way. Agreement across conditions, the basis of the “four stars”, therefore cannot separate a true shared p53 effect from a control artefact. We simulated 33 targets with no true effect, a control error of SD σ<sub>c</sub> Ct cycles, and independent per-well noise σ<sub>t</sub>. The observed pattern, four targets up and three down by at least 2-fold in every condition, is reproduced when σ<sub>c</sub> is about 1.2–2 cycles (σ<sub>t</sub> = 0.3–1.0; Figure 5B). Errors of that size are plausible for single wells with high Ct values, so the simulation shows that the result is not distinguishable from noise, not that it is noise. The largest values face a second limit: a 594-fold change (SNHG4, wild-type p53) needs a 9.2-cycle difference, so with the control at a detection limit of Ct 35 the p53 sample would need a Ct of 25.8 or lower; Emx2os E285K (1,008-fold) needs 25.0 or lower and Zeb2NAT R273H (153-fold) 27.7 or lower (Figure 5C). Raw Ct values would settle this; the deposited fold changes cannot. The detection limit of 35–40 is an assumption.'),
fig('Fig5_shared_control_and_detection_limit','<b>Figure 5. Shared-control and detection-limit diagnostics.</b> (A) Targets ranked by median fold change across conditions. (B) Simulated number of targets passing the criterion with no real effect, against control error; horizontal lines mark the observed counts. (C) Maximum measurable fold change against sample Ct for three assumed detection limits, with the three largest observed values.',3.2*inch),
P('3.4 Mutant-versus-wild-type contrasts that do not depend on the control',h2),
P('The ratio of a mutant’s fold change to wild-type p53’s removes the control value (assuming equal reference behaviour). On that basis seven targets are at least 4-fold higher with mutants than with wild-type p53 in at least three of five mutants: GAS5-family, Emx2os, Zeb2NAT, Air, LUST, DLG2-AS1 and HOTTIP; three are lower: PCGEM1, SNHG4 and snaR (Figure 6). Several of these rest on censored cells and are bounds. The pattern implies that “shared induction” and “mutant-specific gain” are separate groups, and that only the first was highlighted. SNHG4, the largest shared response, is attenuated by every mutant (by 2.7–10-fold), a loss-of-function-like profile.'),
fig('Fig6_mutant_vs_WT_control_free','<b>Figure 6. Mutant relative to wild-type p53, log<sub>2</sub> scale.</b> Rows sorted by mean contrast. Brown and purple labels mark the candidate groups in the text.',7.4*inch),
P('3.5 Do the mutants form groups?',h2),
P('The source report describes a “variant-dependent pattern”. Raw correlations between conditions are all positive (Spearman ρ 0.42–0.81), but that is expected from the shared control (Figure 7A). After removing each target’s median, the condition-specific signal remains weak: wild-type p53, R175H and L194F are mildly positively related (ρ 0.11–0.25), R273H and R280K to each other (ρ 0.10), and both are negatively related to wild-type p53 (−0.35 and −0.59); centring alone sets a baseline of about −0.2 (Figure 7B). Of the 15 ways to place two of six conditions in a separate block, {R273H, R280K} has the strongest contrast (exact permutation p = 1/15 = 0.067). The sources reviewed identify R175H as a conformational and R273H as a DNA-contact mutant; the class of the other three was not established here, so the split is a hypothesis to test, not a finding.'),
fig('Fig7_condition_similarity','<b>Figure 7. Similarity between conditions.</b> Spearman correlation of the 33 fold changes, raw (A) and after removing each target’s median (B).',3.6*inch)]
S+= [P('4. Literature gaps and what this analysis addressed',h1)]
gaps=[['Gap identified','Evidence','Addressed here?'],
['Highlighted lncRNAs other than lincRNA-p21 were never followed up; none other has a mutant p53 link in the literature found [1,7–16]','SNHG4, DLG2-AS1, VLDLR-AS1 studied only as tumour markers or ceRNAs','Partly: mutant-versus-wild-type contrasts (Fig. 3B, 6); validation still needed'],
['“Exactly four” result has no stated threshold','Six targets pass at any increase; four at 2-fold','Yes (Fig. 4A)'],
['Down-regulated group ignored','DISC2, H19 antisense, mascRNA ≥ 5-fold down in all six','Yes (Fig. 2, 4A, 5A)'],
['Global shifts in E285K / R280K not checked','27/33 up, median +2.85 log<sub>2</sub>; p = 0.0003','Diagnosed; correction applied as a sensitivity test (Fig. 4)'],
['Single shared control and no replicates','One value per cell; common denominator','Quantified by simulation (Fig. 5B); cannot be removed without replicates'],
['Extreme fold changes vs detection limit','594-fold and 1,008-fold imply 9–10 cycle differences','Bounded (Fig. 5C); needs raw Ct'],
['Variant-dependent claim untested','Raw correlations inflated by the shared control','Tested on residuals (Fig. 7); exploratory'],
['HEK293T with SV40 large T antigen; no p53 level data per construct','Large T antigen blocks p53 transactivation [17,18]','Not addressable computationally; in roadmap'],
['Which targets form the “36%”','33/90 = 36.7%','Hypothesis only; raw file needed']]
S+= [tbl(gaps,[2.6*inch,2.4*inch,FW-5.0*inch]),Spacer(1,8)]
S+= [P('5. Limitations',h1),
P('The matrix is rebuilt from a report, not from the raw GEO files, and 40 cells are bounds. All tests treat the 33 detected targets as the population, ignore the 57 assays without values, and cannot support significance claims for any single lncRNA. The simulation assumes normally distributed Ct errors, and the detection limits are assumptions. Mutant structural classes beyond R175H and R273H were not established. Methodological references [19–22] are standard works cited from general knowledge and were not re-checked in this session.')]
S+= [P('6. Validation roadmap',h1),
P('A minimal confirmatory design, reusing cell lines from the source study [1]:',body),
P('<b>1. Replicated re-run.</b> At least three independent transfections per construct, including an empty-vector control in triplicate, with raw Ct deposited. Normalise to two or more reference genes (or the global median) and report p53 protein per construct.',bul),
P('<b>2. A p53-null background without large T antigen</b> (for example H1299 or Saos-2) alongside HEK293T, to separate p53 effects from large-T effects.',bul),
P('<b>3. Targets.</b> The four shared lncRNAs (SNHG4 first), the mutant-biased set (Emx2os, Zeb2NAT, DLG2-AS1, HOTTIP, LUST) and the down group (DISC2, H19 antisense, mascRNA) with primers verified outside the array.',bul),
P('<b>4. Endogenous check.</b> Knock down mutant p53 in MDA-MB-231 (R280K), T-47D (L194F) and SK-BR-3 (R175H) and test whether the same lncRNAs fall, as was done for lincRNA-p21 [1].',bul),
P('<b>5. Mechanism.</b> For the shared set, test promoter G-quadruplex binding by ChIP, as shown for lincRNA-p21 [1], rather than assuming a p53 response element.',bul)]
S+= [P('References',h1)]
refs=['He YH, Yeh MH, Chen HF, et al. ERα determines the chemo-resistant function of mutant p53 involving the switch between lincRNA-p21 and DDB2 expressions. Mol Ther Nucleic Acids 2021;25:536–553. doi:10.1016/j.omtn.2021.07.022',
'NCBI GEO accession GSE159185: real-time qPCR of wild-type and mutant p53 transiently overexpressed in HEK-293T cells.',
'Di Agostino S. The impact of mutant p53 in the non-coding RNA world. Biomolecules 2020;10:472.',
'Zhao Y, Li Y, Sheng J, et al. P53-R273H mutation enhances colorectal cancer stemness through regulating specific lncRNAs. J Exp Clin Cancer Res 2019;38:379.',
'LINC01605 is a novel target of mutant p53 in breast and ovarian cancer cell lines. Int J Mol Sci 2023;24:13736.',
'Arenas AM, Cuadros M, Andrades A, et al. LncRNA DLG2-AS1 as a novel biomarker in lung adenocarcinoma. Cancers 2020;12:2080.',
'Wu J, Liu T, Sun L, Zhang S, Dong G. SNHG4 promotes renal cell carcinoma tumorigenesis and invasion by sponging miR-204-5p and upregulating RUNX2. Cancer Cell Int 2020. doi:10.1186/s12935-020-01606-z',
'Long noncoding RNA SNHG4: a novel target in human diseases. Cancer Cell Int 2021 (PMC8557547).',
'lncRNA SNHG4 in non-small cell lung cancer: miR-let-7e and KDM3A/p21 axis (PMC7995486).',
'Beltran M, Puig I, Peña C, et al. A natural antisense transcript regulates Zeb2/Sip1 gene expression during Snail1-induced EMT. Genes Dev 2008;22:756–769.',
'Eroğlu Güneş C, Güçlü E, Vural H, Kurar E. Knockdown of lncRNA ZEB2NAT suppresses EMT, metastasis and proliferation in breast cancer cells. Gene 2021;805:145904.',
'Ma L, Zhang L, Li L, Zhao L. The function of lncRNA EMX2OS/miR-653-5p in lung adenocarcinoma. Open Med 2023. doi:10.1515/med-2023-0686',
'LncRNA EMX2OS, regulated by TCF12, interacts with FUS in prostate cancer cells. OncoTargets Ther 2020;13 (PMC7398891).',
'Cheng Q, Hu T, Zhang W, et al. Down-regulation of lncRNA EMX2OS in cervical squamous cell carcinoma. Folia Biol 2024.',
'Moqadami A, Ahmadi A, Khalaj-Kondori M. lncRNA VLDLR-AS1 gene expression in colorectal cancer. 2023. doi:10.5812/jjcmb-141522',
'Takahashi K, et al. Involvement of extracellular vesicle long noncoding RNA (linc-VLDLR) in tumor cell responses to chemotherapy. Mol Cancer Res 2014 (cited via LNCipedia).',
'Bargonetti J, Reynisdottir I, Friedman PN, Prives C. Site-specific binding of wild-type p53 to cellular DNA is inhibited by SV40 T antigen and mutant p53. Genes Dev 1992;6:1886–1898.',
'T-antigen of SV40 blocks p53 transactivation but not p53 specific binding to DNA. Int J Oncol 1994;5:945–953.',
'Bustin SA, et al. The MIQE guidelines. Clin Chem 2009;55:611–622.',
'Livak KJ, Schmittgen TD. Analysis of relative gene expression data using real-time quantitative PCR and the 2<super>−ΔΔCt</super> method. Methods 2001;25:402–408.',
'Vandesompele J, et al. Accurate normalization of real-time quantitative RT-PCR data by geometric averaging of multiple internal control genes. Genome Biol 2002;3:research0034.',
'Mestdagh P, et al. A novel and universal method for microRNA RT-qPCR data normalization. Genome Biol 2009;10:R64.']
for i,r in enumerate(refs,1): S.append(P(f'[{i}] {r}',ref))
doc=SimpleDocTemplate('../../report/GSE159185_Final_Report.pdf',pagesize=A4,leftMargin=0.65*inch,rightMargin=0.65*inch,topMargin=0.6*inch,bottomMargin=0.65*inch,title='GSE159185 secondary analysis',author='')
doc.build(S,onFirstPage=footer,onLaterPages=footer)
