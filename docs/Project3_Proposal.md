# Research Proposal: The a₀/2 Transition Scale

## Empirical Validation Across Independent Rotation Curve Datasets

_Building on Schneider (2026a, 2026b) — Paper 3 of the Rational Taper Series_

---

## Central Question

Does the median total acceleration at the RT taper radius Rₜ align with a₀/2 in independent datasets, and does this alignment survive controls for BTFR covariance — or is it a SPARC sample artifact?

---

## 1. Background and Motivation

The Rational Taper model contains no reference to a₀. The taper radius Rₜ is derived purely from kinematic fits, motivated by asymptotic stability. Yet Paper 2 reports that the median total centripetal acceleration at Rₜ across 98 spatially resolved SPARC galaxies is:

**g(Rₜ) = 6.51 × 10⁻¹¹ m s⁻² ≈ a₀/2 (offset: 8.5%)**

This is unexpected for three reasons. First, a₀ appears nowhere in the RT model's construction or fitting. Second, the alignment holds across three decades in baryonic mass. Third, the factor of two is not arbitrary in MOND phenomenology — a₀/2πG appears independently as a characteristic surface density in galactic scaling relations.

Paper 2 explicitly defers interpretation and flags two unresolved concerns this proposal addresses:

**Concern 1 — BTFR Covariance.** The iso-luminosity tracks at Rₜ follow g ∝ 1/Rₜ with slope 0.514, consistent with BTFR expectations. Is the a₀/2 alignment genuine, or does it follow trivially from the sample's mass demographics?

**Concern 2 — Sample Specificity.** SPARC is a single catalog with specific selection effects and photometric methodology. The alignment must survive in independent datasets to be credible.

---

## 2. Research Questions

1. Does the a₀/2 alignment survive in independent datasets (THINGS, LITTLE THINGS, PROBES) with different selection functions and baryonic decomposition methods?
2. Is the alignment a BTFR artifact? After removing the primary mass-velocity trend, does the residual g(Rₜ) distribution still center on a₀/2?
3. Does a constrained RT model fixing Rₜ by g(Rₜ) = a₀/2 perform comparably to the free-Rₜ model by BIC? If so, the alignment is physically load-bearing, not statistical coincidence.
4. How does g(Rₜ) vary across morphology, surface brightness, and environment? Is the alignment universal?

---

## 3. Datasets

| Dataset       | N     | Role                     | Baryonic Quality                 | SPARC Overlap |
| ------------- | ----- | ------------------------ | -------------------------------- | ------------- |
| SPARC         | 175   | Primary anchor           | Gold standard (3.6μm Spitzer)    | —             |
| THINGS        | 34    | Independent replication  | High (independent HI reductions) | Partial       |
| LITTLE THINGS | 41    | Dwarf-regime stress test | Moderate (UV+optical SFH M/L)    | Minimal       |
| PROBES        | ~2700 | Statistical breadth      | Variable (optical photometry)    | Partial       |

**SPARC** anchors the analysis. All QC procedures from Papers 1 and 2 carry forward. The baseline g(Rₜ) result must be reproduced as an internal consistency check before proceeding.

**THINGS** is the critical replication test. Partial overlap with SPARC allows direct galaxy-by-galaxy comparison of fitted Rₜ and g(Rₜ) between independent pipelines. Agreement strengthens the replication claim; systematic offsets indicate pipeline artifacts.

**LITTLE THINGS** stress-tests the dwarf regime, where MOND phenomenology is strongest and the RT model's low-mass behavior is least constrained. If a₀/2 holds there, it is not a massive-galaxy artifact.

**PROBES** provides statistical power for morphological and environmental stratification, but noisier baryonic decompositions mean its results should be treated as exploratory and cross-checked against SPARC+THINGS.

**Key risk:** SPARC's photometric quality is what makes its baryonic decompositions trustworthy. Weakened alignment in THINGS or PROBES must be evaluated against their increased baryonic mass uncertainty before concluding the signal has failed.

---

## 4. Methods

### 4.1 BTFR Covariance Test

This is the central methodological challenge: separating a genuine physical signal from a demographic artifact.

1. Fit the free-Rₜ RT model to all galaxies in each dataset using the existing pipeline. Compute g(Rₜ) = V²ₜₒₜ/Rₜ for all spatially resolved galaxies.
2. Fit the primary BTFR power-law trend: log g = α log M_bar + β. Compute residuals Δlog g = log g(Rₜ) − fitted trend.
3. Test whether the residual distribution centers on zero (BTFR artifact) or on a nonzero offset corresponding to a₀/2 (physical alignment). Primary statistical test: signed-rank test against zero.
4. Stratify by baryonic mass quartile to determine whether the alignment holds across the mass spectrum or is driven by a particular regime.

### 4.2 Constrained Model BIC Test

If the a₀/2 alignment is physical, enforcing it should not cost much fit quality. The constrained RT model fixes Rₜ by requiring:

**g(Rₜ) ≡ a₀/2 → Rₜ = 2V²ₜₒₜ / a₀**

This reduces the RT model from two free parameters (ω, Rₜ) to one (ω). BIC comparison between free and constrained models is the critical test:

- If constrained model is BIC-competitive: the a₀/2 condition is not costing fit quality — strong evidence the alignment is physically meaningful.
- If constrained model loses substantially: fitted Rₜ values are not well-described by g(Rₜ) = a₀/2 individually, and the alignment is a population-level central tendency only.

Both outcomes are scientifically informative.

### 4.3 Cross-Dataset Comparison

For galaxies appearing in both SPARC and THINGS, compare fitted Rₜ and g(Rₜ) directly. The primary metric is the median g(Rₜ) and its 16th–84th percentile range. Secondary: fraction of galaxies with g(Rₜ) within a₀/2 ± 10%, and Spearman correlation of residuals with structural properties.

### 4.4 Exploratory: Morphological and Environmental Dependence

Using PROBES for statistical power, investigate whether g(Rₜ) varies systematically with morphological type, central surface brightness (LSB vs. HSB, following the Paper 2 Σ₀ analysis), and isolation metrics where available. These are explicitly hypothesis-generating, not confirmatory.

---

## 5. Possible Outcomes

All three outcomes are publishable. The paper's contribution is closing the open question raised in Paper 2, not confirming a preferred result.

**Strong replication.** Residuals center on a₀/2 after BTFR removal; constrained model is BIC-competitive; THINGS agrees with SPARC galaxy-by-galaxy. _Interpretation:_ a₀/2 is a genuine kinematic transition scale with physical content; strong quantitative constraint for theoretical models.

**Partial replication.** Alignment survives BTFR control in SPARC but weakens in THINGS; constrained model loses to free model. _Interpretation:_ alignment is a population-level central tendency, not a universal per-galaxy constraint; physically suggestive but not decisive.

**BTFR artifact.** Residuals center on zero after BTFR removal; alignment disappears in THINGS. _Interpretation:_ a₀/2 is a demographic artifact of SPARC's sample composition. Still publishable as a rigorous null result that definitively closes the question.

---

## 6. Project Phases

| Phase | Title                | Key Tasks                                                                         | Deliverable                            |
| ----- | -------------------- | --------------------------------------------------------------------------------- | -------------------------------------- |
| 1     | SPARC Baseline       | Re-run existing pipeline; reproduce Paper 2 g(Rₜ) result                          | Confirmed baseline distribution        |
| 2     | BTFR Covariance Test | Fit BTFR trend; compute residuals; signed-rank test; mass-quartile stratification | Physical vs. artifact determination    |
| 3     | Constrained Model    | Implement g(Rₜ) = a₀/2 constraint; refit SPARC; BIC comparison                    | Per-galaxy and catalog-level BIC table |
| 4     | THINGS Replication   | Acquire data; build baryonic pipeline; galaxy-by-galaxy comparison                | Independent g(Rₜ) in 34 galaxies       |
| 5     | LITTLE THINGS        | Acquire data; fit RT model; assess dwarf-regime behavior                          | Dwarf-regime stress test result        |
| 6     | PROBES Breadth       | Acquire and quality-filter; morphology/environment stratification                 | Exploratory trend analysis             |
| 7     | Synthesis            | Integrate results; assess outcome scenario; draft manuscript                      | Submitted Paper 3                      |

---

## 7. Scope Boundaries

**In scope:**

- Empirical characterization of g(Rₜ) distribution across datasets
- BIC test of constrained vs. free RT model
- Explicit BTFR covariance control and residual analysis
- Cross-dataset galaxy-by-galaxy pipeline comparison (THINGS)
- Exploratory morphological and environmental stratification

**Out of scope:**

- Proposing a physical mechanism for the a₀/2 alignment
- Claiming the RT model explains dark matter or falsifies MOND
- Extension to pressure-supported systems (ellipticals, clusters)
- Redshift evolution (requires data not yet available)
- Building a new theoretical framework

---

## 8. Anticipated Contribution

Regardless of outcome, this paper will:

- Close the open question explicitly flagged in Paper 2
- Demonstrate RT model generalizability beyond SPARC
- Provide the most rigorous characterization of the RT transition scale to date, constituting a quantitative constraint any physical theory must satisfy
- Connect RT phenomenology to the MOND acceleration scale literature in a way that is either strongly supportive (strong replication) or precisely delimiting (artifact result)
