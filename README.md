# LncAudit53

> **Description:** A reproducible audit of a public p53 lncRNA qPCR array (GEO GSE159185): literature-gap review, robustness tests and control-independent re-analysis of wild-type and mutant p53 effects in HEK293T cells.

**Topics:** `p53` `lncRNA` `qPCR` `GEO` `GSE159185` `reproducible-research` `bioinformatics` `secondary-analysis`

**Authors:** Muhammad Hasnat, Hammad Anwar

Literature-gap review, robustness testing and control-independent re-analysis of the public qPCR array dataset
[GSE159185](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE159185) (wild-type and five hotspot mutant p53 constructs, HEK293T cells,
SBI LncProfiler 90-lncRNA array). This is a **computational secondary analysis of deposited data; no new wet-lab experiments were performed.**

**Main deliverable:** [`report/GSE159185_Final_Report.pdf`](report/GSE159185_Final_Report.pdf)

## Key findings

| # | Finding | Figure |
|---|---------|--------|
| 1 | The "four lncRNAs up in every condition" result depends on an unstated threshold: six targets pass at any increase, four at 2-fold. | Fig. 4A |
| 2 | Three targets (DISC2, H19 antisense, mascRNA) are at least 5-fold down in all six conditions; this mirror group is not discussed in the source. | Fig. 2, 5A |
| 3 | E285K shifts the median of all 33 targets by +2.85 log2 (about 7-fold; 27/33 up, sign test p = 0.0003), suggesting a normalisation or input offset. | Fig. 4B,C |
| 4 | Every condition shares one control value, so cross-condition agreement is not independent evidence. A control error of about 1.2-2 Ct cycles reproduces the observed 4-up / 3-down pattern in simulation. | Fig. 5B |
| 5 | Mutant-vs-wild-type contrasts (control-independent) show SNHG4 is strongest with wild-type p53, while Emx2os, Zeb2NAT, DLG2-AS1, HOTTIP and others rise more with specific mutants. | Fig. 3B, 6 |

## Repository layout

```
data/raw/         Original report PDFs and the heatmap image extracted from them (inputs)
data/processed/   Reconstructed 33 x 6 matrix with per-cell source flags; mutant-minus-WT table
src/              Numbered scripts 01-05 (digitise, analyse, simulate, plot, build report)
results/figures/  All figures, PNG (300 dpi) and vector PDF
report/           Final PDF report
docs/             Data provenance, methods, and anticipated questions
run_all.sh        Re-creates everything from data/raw/
```

## Reproduce

```bash
pip install -r requirements.txt
./run_all.sh
```

Tested with Python 3.12. Fonts: Carlito if installed, otherwise DejaVu Sans.

## Data provenance (please read)

The raw GEO files (`GSE159185_non-normalized_data.txt.gz`, `GSE159185_fold-change_data.txt.gz`) were **not** available to this analysis.
The 33 x 6 fold-change matrix was **reconstructed** from (a) exact values printed in the report tables
(4 lncRNAs x 6 conditions plus the three largest values) and (b) the original heatmap, read back through its own colour bar.
Every cell in `data/processed/Supplementary_Table_S1_reconstructed_matrix.csv` carries a `source` flag
(`exact`, `digitized`, or `censored` = bound at |log2 FC| >= 5). Validation: mean absolute error 0.034 log2 units on 14 known cells
(max 0.073, r = 0.9996), and the number of targets above 1.0 matched the report in all six conditions.
Full details: [`docs/DATA_PROVENANCE.md`](docs/DATA_PROVENANCE.md). Anyone can replace the reconstructed matrix with the raw GEO
values to confirm or correct the results.

## Limitations

No replicates (one value per target per condition), so no significance claims are made for individual lncRNAs.
Forty cells are bounds. Detection limits used in the simulation (Ct 35-40) are assumptions. See report Section 5.

## Citation

See [`CITATION.cff`](CITATION.cff). Please also cite the source study: He YH et al., *Mol Ther Nucleic Acids* 2021;25:536-553,
doi:10.1016/j.omtn.2021.07.022, and the GEO accession GSE159185.

## Licence

Code: MIT. Report text, figures and processed tables: CC BY 4.0. See [`LICENSE`](LICENSE).
The GEO dataset and the source study remain the property of their original authors.

## Acknowledgements

Analysis and drafting were assisted by an AI tool (Claude, Anthropic). The authors are responsible for the content.
# LncAudit53
