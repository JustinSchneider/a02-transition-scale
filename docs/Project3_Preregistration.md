# Pre-Registration: The a₀/2 Transition Scale

## Empirical Validation Across Independent Rotation Curve Datasets

**Author:** Justin Schneider  
**Date:** 09 APR 2026
**Repository:** https://github.com/JustinSchneider/a02-transition-scale
**Status:** Pre-registration — analysis plan locked prior to acquisition of THINGS data

---

## Purpose

This document records my pre-specified analysis plan for Paper 3 of the Rational Taper series. It exists to establish a clean line between hypothesis generation (Paper 2) and hypothesis testing (Paper 3).

In Paper 2 (Schneider 2026b), I found that the median total centripetal acceleration at the RT taper radius Rt across 98 spatially resolved SPARC galaxies is 6.51 × 10⁻¹¹ m/s², within 8.5% of a₀/2 = 6.00 × 10⁻¹¹ m/s². This result was not predicted in advance — it was found by examining the data. Paper 3 is the confirmatory test. The analysis decisions recorded here were made before I downloaded any THINGS rotation curve data.

---

## Pre-Specified Hypothesis

The median total centripetal acceleration at the RT taper radius, after controlling for BTFR covariance, aligns with a₀/2 in datasets independent of SPARC.

Specifically: in the THINGS dataset, after removing the primary baryonic mass trend using the externally calibrated slope from Paper 2, the residual g(Rt) distribution will center within ±10% of a₀/2 at the median.

---

## Locked Analysis Definitions

These definitions cannot change after THINGS data is acquired. Any departure from them must be reported as a deviation in the paper.

**1. Definition of g(Rt)**  
g(Rt) = V²_model(Rt) / Rt, where V_model is the fitted RT model velocity evaluated at the taper radius, not the interpolated observed velocity. Rationale: V_model is a smoothed, noise-reduced estimator consistent with how Rt is defined. The observed-velocity version g_obs(Rt) = V²_obs,interp(Rt) / Rt will be computed and reported as a supplementary robustness check.

**2. Definition of "spatially resolved"**  
A galaxy's taper radius is spatially resolved if Rt < R_max, where R_max is the maximum observed galactocentric radius in the rotation curve data. Galaxies with Rt ≥ R_max are excluded from the g(Rt) analysis.

**3. BTFR slope (fixed — not refitted to Paper 3 data)**  
α = 0.238, taken from the empirical scaling log g(Rt) ∝ 0.238 × log M_bar found in Paper 2 (Schneider 2026b, Section 6). This slope is treated as an external calibration and will not be re-estimated from Paper 3 data. Re-estimating would bias residuals toward zero by construction.

**4. BTFR intercept anchor**  
The intercept is fixed by requiring the Paper 2 SPARC resolved sample (N=98) to center on a₀/2:

β = log₁₀(a₀/2) − α × log₁₀(M̃_bar)

where M̃_bar = 5.80 × 10⁹ M☉ (median baryonic mass of the Paper 2 resolved sample, N=98, computed from SPARC 3.6 µm luminosities with Y_disk=0.5, Y_bulge=0.7, and 1.33×MHI for helium), and a₀/2 = 6.00 × 10⁻¹¹ m/s².

Numerical value: β = −12.55

The BTFR reference line for residual computation is therefore:
log₁₀ g_trend = 0.238 × log₁₀(M_bar) − 12.55

**5. Primary statistical test**  
Wilcoxon signed-rank test applied to the residuals Δlog g = log₁₀(g(Rt)) − log₁₀(g_trend), using the THINGS dataset as the out-of-sample test. The null hypothesis is that residuals center on zero (BTFR fully accounts for g(Rt), no additional signal). A two-sided p < 0.05 is the threshold for rejecting the null.

**6. Alignment threshold**  
The a₀/2 alignment is considered replicated if the median g(Rt) in the confirmatory dataset is within ±10% of a₀/2 after BTFR residual removal, consistent with the ±10% corridor used in Paper 2.

**7. Dataset roles**

| Dataset              | Role                   | Rationale                                                                                           |
| -------------------- | ---------------------- | --------------------------------------------------------------------------------------------------- |
| SPARC (N=175)        | Baseline anchor        | Internal consistency check; reproduces Paper 2 result before proceeding                             |
| THINGS (N=34)        | Primary confirmatory   | Independent HI kinematics; partial SPARC overlap enables cross-pipeline comparison                  |
| LITTLE THINGS (N=41) | Secondary confirmatory | Dwarf-regime stress test; extends mass baseline below SPARC                                         |
| PROBES (~2700)       | Exploratory only       | Optical photometry introduces ≥0.2 dex M/L uncertainty; not cited as confirmatory evidence for a₀/2 |

**8. Constrained model BIC test**  
I will implement a one-parameter constrained RT model in which Rt is determined at each fit evaluation by numerically solving g(Rt) = a₀/2. The BIC of this constrained model (k=1) will be compared to the free RT model (k=2). The primary BIC test uses only galaxies satisfying n ≥ 20 radial points AND Rt/R_max < 0.5 (taper clearly resolved). The BIC competitiveness threshold is |ΔBIC| < 2 at the median of this subsample.

---

## Pre-Analysis Power Assessment

Prior to downloading THINGS data, I ran a power simulation (10,000 iterations, seed=42) using the Wilcoxon signed-rank test at p < 0.05 (two-sided). Two noise models were evaluated:

- **Conservative (raw scatter):** 0.409 dex — half the 16th–84th percentile range of log₁₀(g(Rt)/a₀/2) in the Paper 2 resolved sample
- **Realistic (BTFR-residual scatter):** 0.277 dex — half the 16th–84th percentile range after removing the fixed-slope BTFR trend

The realistic model is more appropriate for the pre-specified test (which operates on BTFR residuals), but the conservative model is reported as a lower bound on power.

**Conservative noise model (0.409 dex):**

| N   | Type I rate | Power (0.1 dex) | Power (0.2 dex) | Decision |
| --- | ----------- | --------------- | --------------- | -------- |
| 15  | 0.049       | 0.138           | 0.401           | PILOT    |
| 20  | 0.048       | 0.169           | 0.533           | PILOT    |
| 25  | 0.046       | 0.202           | 0.626           | PRIMARY  |
| 34  | 0.048       | 0.268           | 0.773           | PRIMARY  |

**Realistic noise model (0.277 dex):**

| N   | Type I rate | Power (0.1 dex) | Power (0.2 dex) | Decision |
| --- | ----------- | --------------- | --------------- | -------- |
| 15  | 0.049       | 0.242           | 0.708           | PRIMARY  |
| 20  | 0.046       | 0.317           | 0.842           | PRIMARY  |
| 25  | 0.048       | 0.392           | 0.925           | PRIMARY  |
| 34  | 0.045       | 0.513           | 0.981           | PRIMARY  |

**Interpretation:** Under the conservative noise model, THINGS requires N ≥ 25 resolved galaxies to reach 60% power. Under the realistic model, even N = 15 exceeds the threshold. Given that Paper 2 resolved 70.3% of SPARC galaxies, THINGS (N=34) is expected to yield ~20–24 resolved galaxies. This places the realistic-model power at 84–93% and the conservative-model power at 53–63%.

**Decision:** THINGS is designated as the **primary** confirmatory dataset, with the caveat that if the resolved subsample falls below N = 20, results should be interpreted with the conservative power estimate in mind. If the resolved count falls below N = 15, THINGS will be redesignated as pilot data and LITTLE THINGS combined with PROBES will serve as the primary statistical test.

---

## What I Will Not Do

- I will not test other fractions of a₀ (a₀/3, a₀/4, etc.) as formal hypotheses. If the empirical median differs from a₀/2, I will report what it is.
- I will not refit the BTFR slope to Paper 3 data.
- I will not adjust the alignment threshold after seeing results.
- I will not cite PROBES results as confirmatory evidence for the a₀/2 alignment.
- I will not add post-hoc exclusion criteria after examining THINGS g(Rt) values.

---

## Possible Outcomes

All three outcomes are publishable. I am not pre-committed to a preferred result.

**Strong replication:** THINGS residuals center on a₀/2 after BTFR removal; constrained model is BIC-competitive; THINGS and SPARC agree galaxy-by-galaxy for overlap objects. Interpretation: a₀/2 is a genuine kinematic transition scale.

**Partial replication:** Alignment survives BTFR control in SPARC but weakens in THINGS; constrained model loses to free model. Interpretation: population-level central tendency, not a universal per-galaxy constraint.

**BTFR artifact:** Residuals center on zero after BTFR removal; alignment disappears in THINGS. Interpretation: a₀/2 is a demographic artifact of SPARC's sample composition. This is a valid and important null result.

---

## Deviations

Any deviation from this pre-registered plan will be explicitly disclosed in the paper, including the reason for the deviation and its potential effect on the results.
