# Data provenance

## Inputs (data/raw/)
| File | What it is |
|------|-----------|
| `source_report_v1.pdf` | First version of the descriptive report (no heatmap on its page 4). |
| `source_report_v2_light.pdf` | Second version, identical text, includes the original heatmap. |
| `heatmap_original_from_report.png` | Heatmap extracted from v2 with `pdfimages -png -f 4 -l 4` (1800 x 1366 px). |

## What was NOT available
`GSE159185_non-normalized_data.txt.gz` and `GSE159185_fold-change_data.txt.gz` could not be downloaded in the analysis environment
(NCBI blocked automated access). All numbers therefore trace to the report tables and heatmap.

## How the matrix was reconstructed (src/01, src/02)
1. Heatmap frame and colour bar located by pixel search (colour bar spans log2 FC -5 to +5).
2. For each of the 33 x 7 cells, the median colour of a central patch was matched to the nearest colour-bar entry.
3. Exact values from report tables override the heatmap for: DLG2AS (family), lincRNA-p21, lincRNA-VLDLR, SNHG4 (all six conditions);
   Emx2os in E285K (1,007.95); Zeb2NAT in R273H (153.06).
4. Cells at the colour-scale limit (|log2 FC| >= 4.9) are flagged `censored` unless exactly known; they are used as bounds only.

## Validation
- 14 uncensored cells with known exact values: mean absolute error 0.034 log2, maximum 0.073, r = 0.9996.
- Targets above 1.0 per condition: recovered 16, 12, 12, 17, 22, 27, identical to the report.
- Control column read as 0.07 (true value 0); set to exactly 0 in analysis. Precision near zero is about +/-0.1 log2.
- 52 cells reached the colour limit; 12 are resolved exactly by the report tables, 40 remain bounds.

## To confirm with raw data
Replace `data/processed/Supplementary_Table_S1_reconstructed_matrix.csv` values with log2 of the GEO fold-change file and re-run
`src/03` to `src/05`. Raw Ct values would also allow checking the detection-limit concern (report Fig. 5C).
