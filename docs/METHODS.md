# Methods summary

- **Global-shift diagnostic:** per condition, share of 33 targets above 1.0 with exact (Clopper-Pearson) 95% CI and two-sided
  exact binomial sign test against 0.5; median log2 FC. Diagnostic only: targets are not independent and are limited to detected targets.
- **Median-centring:** each condition's median log2 FC subtracted (global-median normalisation) and the all-six criteria re-evaluated.
- **Threshold sensitivity:** number of targets with FC at or above T in every condition, T from >1 to >=4, up and down.
- **Shared-control simulation:** 33 targets, no true effect; log2 FC = -delta + epsilon, delta ~ N(0, sigma_c) shared by a target's six
  conditions, epsilon ~ N(0, sigma_t) per condition; 4,000 simulations per grid point (seed 7).
- **Detection-limit envelope:** maximum measurable fold change = 2^(LOD - sample Ct) for assumed LOD Ct 35, 37, 40.
- **Control-free contrast:** log2(mutant FC / WT FC) per target; cells where both terms are censored in the same direction are not estimable.
  Candidate groups: at least 4-fold above (or below) WT in at least 3 of 5 mutants.
- **Condition similarity:** Spearman correlation of 33 fold changes, raw and after subtracting each target's median; all 15 two-versus-four
  splits of the conditions ranked by within-minus-between block correlation (exact permutation).
