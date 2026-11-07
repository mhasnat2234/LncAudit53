# Anticipated questions

**Why is the data reconstructed instead of the raw GEO files?** The raw files could not be downloaded in the analysis environment. The
reconstruction is validated (see DATA_PROVENANCE.md) and every cell is flagged by source, so reviewers can see which values are exact.

**Are any claims of statistical significance made for individual lncRNAs?** No. The data have one value per target per condition.
The only p-values are a diagnostic sign test across targets (global shift) and an exploratory exact permutation over condition splits.

**Is this a new experiment?** No. It is a computational secondary analysis of public deposited data.

**Are the Ct detection limits real measurements?** No, they are assumptions used to bound plausibility; raw Ct values are needed.

**Who wrote the paper whose data are analysed?** He et al., Mol Ther Nucleic Acids 2021. We do not claim authorship of that study or dataset.

**Was AI used?** Yes, analysis and drafting were assisted by Claude (Anthropic); the authors are responsible for the content.
