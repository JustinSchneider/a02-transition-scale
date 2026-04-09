# Methodology: The a0/2 Transition Scale

**Author:** Justin Schneider (independent researcher)
**Last Updated:** April 2026
**Status:** Pre-analysis gates complete; THINGS data acquired (de Blok et al. 2008, 19 galaxies)

---

## 1. Scientific Background

In Paper 2 (Schneider 2026b), the Rational Taper model was benchmarked against NFW and MOND across 175 SPARC galaxies. A post-hoc finding emerged: the median total centripetal acceleration at the RT taper radius across 98 spatially resolved galaxies is 6.51 x 10^-11 m/s^2, within 8.5% of a0/2 = 6.00 x 10^-11 m/s^2. This alignment was not predicted in advance.

Paper 3 is the confirmatory test. The central question: **is the a0/2 alignment a genuine kinematic transition scale, or a SPARC sample artifact?**

The analysis plan was pre-registered before acquiring any THINGS data (commit dc24e10, 09 APR 2026). All definitions, thresholds, and statistical tests were locked prior to seeing new data. See `docs/Project3_Preregistration.md` for the full pre-registered plan.

---

## 2. Datasets

| Dataset                            | N     | Role                   | Rationale                                            |
| ---------------------------------- | ----- | ---------------------- | ---------------------------------------------------- |
| SPARC (Lelli et al. 2016)          | 175   | Baseline anchor        | Internal consistency; reproduces Paper 2 result      |
| THINGS (de Blok et al. 2008)       | 19    | Primary confirmatory   | Independent HI kinematics; 13-galaxy SPARC overlap   |
| LITTLE THINGS (Hunter et al. 2012) | 41    | Secondary confirmatory | Dwarf-regime stress test                             |
| PROBES (~2700)                     | ~2700 | Exploratory only       | Optical M/L uncertainty >= 0.2 dex; not confirmatory |

### SPARC Data Files

- **Rotation curves:** `data/raw/MassModels_Lelli2016c.mrt` — per-galaxy radial profiles of V_obs, V_err, V_gas, V_disk, V_bulge
- **Galaxy metadata:** `data/raw/SPARC_Lelli2016c.mrt` — distance, inclination, luminosity, HI mass, quality flag
- **Bulge luminosities:** `data/raw/Bulges.mrt` — L_bulge at 3.6 um for 175 galaxies (32 nonzero)
- **Paper 2 fit results:** `data/galaxy_dynamics.db` (SQLite), `data/supplemental/Tournament_Results.csv`

### THINGS/SPARC Overlap

The THINGS survey parent sample comprises 34 galaxies (Walter et al. 2008), but rotation curves were published for only 19 of them (de Blok et al. 2008). Of these 19, 13 appear in SPARC and 6 are THINGS-only. Full cross-reference in `data/things_sparc_overlap.csv`. For overlap galaxies, SPARC baryonic profiles are used with THINGS rotation curves (isolates the kinematic data while holding the baryonic baseline constant). Non-overlap galaxies require independent baryonic decomposition from S4G + THINGS HI maps.

---

## 3. Baryonic Velocity Computation

All models use an identical baryonic mass decomposition (Lelli et al. 2016, Eq. 2):

```
V_bary(R) = sqrt( |V_gas| * V_gas + Y_d * |V_disk| * V_disk + Y_b * |V_bulge| * V_bulge )
```

The |V| \* V sign-preserving convention handles negative SPARC velocity contributions (e.g., central gas depressions) without producing imaginary numbers.

**Fixed mass-to-light ratios (same for all models):**

| Parameter | Symbol  | Value |
| --------- | ------- | ----- |
| Disk M/L  | Y_disk  | 0.5   |
| Bulge M/L | Y_bulge | 0.7   |

SPARC provides V_disk and V_bulge at Y = 1; the mass-to-light ratios enter as direct multipliers on the V^2 terms.

**Implementation:** `src/physics.py :: compute_v_bary()`

---

## 4. The Rational Taper Model (Free, k = 2)

**Free parameters:** omega (km/s/kpc), R_t (kpc)

```
V_model(R) = V_bary(R) + omega * R / (1 + R / R_t)
```

- At small R: correction -> omega \* R (linear)
- At large R: correction -> omega \* R_t = V_sat (constant; flat rotation curve recovery)
- **Additive coupling** (not quadrature), consistent with Papers 1 and 2

**Fitting details:**

- Multi-start optimizer: 4 initial conditions [(5,5), (10,2), (2,15), (20,R_max)]; retain lowest chi^2
- Bounds: omega in [0, 200] km/s/kpc, R_t in [0.1, 5 * R_max] kpc
- `scipy.optimize.curve_fit` with `absolute_sigma=True` (Levenberg-Marquardt)
- Derived quantity: V_sat = omega \* R_t (not a fit parameter)

**Implementation:** `src/physics.py :: fit_rational_taper()`

### Transition Diagnostics

At R = R_t, the taper correction equals omega \* R_t / 2 (exactly half the saturation velocity). The key diagnostic quantity is the total centripetal acceleration at the taper radius:

```
g(R_t) = V_model(R_t)^2 / R_t
```

where V_model(R_t) = V_bary(R_t) + omega * R_t / 2. This uses the *model\* velocity, not the observed velocity — the model is a smoothed, noise-reduced estimator consistent with how R_t is defined. The observed-velocity version g_obs(R_t) = V_obs,interp(R_t)^2 / R_t is computed and reported as a supplementary robustness check.

A galaxy's taper radius is **spatially resolved** if R_t < R_max (fitted taper inside the data range). Galaxies with R_t >= R_max are excluded from the g(R_t) analysis.

**Implementation:** `src/physics.py :: compute_transition_diagnostics()`

**Baryonic interpolation** at arbitrary radii uses the signed-square convention: interpolate V^2_bary (with sign), then recover velocity as sign(x) \* sqrt(|x|). This avoids imaginary values when profiles cross zero.

**Implementation:** `src/physics.py :: interpolate_v_bary()`

---

## 5. The Constrained RT Model (k = 1)

**Free parameter:** omega only. R_t is derived from the constraint g(R_t) = a0/2.

### The Constraint Equation

Imposing g(R_t) = a0/2 at R = R_t:

```
[V_bary(R_t) + omega * R_t / 2]^2 / R_t = a0/2
```

This is a transcendental equation in **both** omega and R_t simultaneously, because V_bary(R_t) depends on R_t through the rotation curve profile. You cannot fix R_t from Paper 2's median and fit only omega — that fixes the wrong thing.

### Implementation Strategy

The constrained model is a 1D optimization over omega with an inner numerical root-find over R_t at each evaluation:

1. For each candidate omega, solve the constraint equation for R_t using `scipy.optimize.brentq`
2. Scan the full radial range on a fine grid (200 points) for all sign changes in the constraint function
3. Use `brentq` on each bracketed interval
4. Take the smallest positive root (physically: innermost transition)
5. The outer optimizer (`scipy.optimize.minimize_scalar`, bounded) finds the omega that minimizes chi^2

**Multiple roots:** The constraint equation may have more than one sign change, especially in galaxies with non-monotonic V_bary. The smallest positive root is used; galaxies with multiple roots are flagged.

**No-solution galaxies:** When no sign change exists in the data range for any omega, the constraint cannot be satisfied. These galaxies are tracked as a result, not discarded.

**BIC:** k = 1 for this model (only omega is free; R_t is derived).

**Implementation:** `src/physics.py :: find_constrained_rt()`, `src/physics.py :: fit_constrained_rt()`

### Gate 6 Validation Results (SPARC Baseline)

The constrained model was run on all 123 resolved SPARC galaxies (including 25 that were resolved but excluded from the N=98 g(R_t) analysis due to interpolation failures):

| Metric                                  | Value          |
| --------------------------------------- | -------------- |
| No valid R_t solution                   | 92 / 123 (75%) |
| Converged                               | 31 / 123 (25%) |
| Constrained wins on BIC                 | 1 / 31         |
| Free model wins on BIC                  | 30 / 31        |
| Median delta_BIC (constrained - free)   | +96.88         |
| High-quality subsample median delta_BIC | +162.01 (N=5)  |

**Interpretation:** The constraint g(R_t) = a0/2 is too restrictive for individual galaxies. The population median aligns with a0/2, but no per-galaxy physical law operates. This is the "population-level central tendency, not a per-galaxy law" outcome anticipated in the risk register. Both outcomes are publishable; this result sharpens the Paper 3 narrative toward testing whether the _population_ alignment replicates.

---

## 6. BTFR Covariance Test

The a0/2 alignment could be an artifact of baryonic mass scaling (the baryonic Tully-Fisher relation). To control for this:

**Fixed-slope residuals:**

```
log10(g_trend) = 0.238 * log10(M_bar) - 12.55
```

- Slope alpha = 0.238: from Paper 2's empirical scaling, fixed as external calibration (NOT refitted to Paper 3 data)
- Intercept beta = -12.55: anchored so the Paper 2 SPARC resolved sample (N=98, median M_bar = 5.80 x 10^9 Msun) centers on a0/2

**Residuals:**

```
delta_log_g = log10(g(R_t)) - log10(g_trend)
```

**Primary test:** Wilcoxon signed-rank test on THINGS residuals (out-of-sample). If residuals center on zero: BTFR fully explains g(R_t) in the new dataset (artifact). If residuals center on a nonzero offset corresponding to a0/2: alignment is genuine.

**Implementation:** Numerical anchors verified in `scripts/power_assessment.py`. Statistical test to be implemented in Phase 2 analysis.

---

## 7. Model Comparison: BIC

```
BIC = chi^2 + k * ln(n)
```

where chi^2 is the total weighted chi-squared, k is the number of free parameters, and n is the number of data points.

| Model          | k   | BIC penalty at n=15 |
| -------------- | --- | ------------------- |
| Constrained RT | 1   | 2.71                |
| Free RT        | 2   | 5.42                |

**Competitiveness threshold:** |delta_BIC| < 2 at the median of the high-quality resolved subsample (n >= 20, R_t/R_max < 0.5). From Gate 6 validation, this threshold is **not met** on SPARC — the constrained model loses decisively.

**Implementation:** `src/physics.py :: compute_bic()`

---

## 8. Database Schema

The SQLite database (`data/galaxy_dynamics.db`) uses three tables inherited from Paper 2:

- **`galaxies`** (175 rows) — Metadata: distance, inclination, luminosity, disk scale length, quality flag
- **`radial_profiles`** (3391 rows) — Per-radius velocity components from SPARC
- **`model_fits`** (700 + Paper 3 additions) — One row per (galaxy, model) fit

The `model_fits` table uses generic `param1`/`param2` columns whose semantics depend on `model_name`:

| model_name       | param1            | param2             | k   |
| ---------------- | ----------------- | ------------------ | --- |
| `rational_taper` | omega (km/s/kpc)  | R_t (kpc)          | 2   |
| `constrained_rt` | omega (km/s/kpc)  | R_t (kpc, derived) | 1   |
| `nfw`            | c (concentration) | V_200 (km/s)       | 2   |
| `mond_fixed`     | a_0 (fixed)       | --                 | 0   |
| `mond_free`      | a_0 (best-fit)    | --                 | 1   |

**Implementation:** `src/database.py`

---

## 9. Fitting Conventions

- **Optimizer:** `scipy.optimize.curve_fit` with `absolute_sigma=True` (free RT, Levenberg-Marquardt); `scipy.optimize.minimize_scalar` with bounded method (constrained RT)
- **Error handling:** Zero/negative errors replaced with minimum nonzero error (floor: 1.0 km/s). **Implementation:** `src/physics.py :: _safe_errors()`
- **Multi-start:** Free RT uses 4 initial conditions; retain lowest chi^2
- **Convergence:** `curve_fit` RuntimeError -> converged=False; `minimize_scalar` with no valid R_t -> converged=False
- **Constrained root-finding:** `scipy.optimize.brentq` on 200-point grid scan for sign changes; smallest positive root preferred

---

## 10. Software Architecture

```
Raw Data (MRT, CSV, SQLite)
    |
    v
src/ingest.py  -- load profiles and metadata from DB and SPARC files
    |
    v
src/fit.py  -- orchestration: Phase 1 baseline, Gate 6 validation, CLI
    |  calls
    v
src/physics.py  -- all model equations, fitting, diagnostics
    |  stores via
    v
src/database.py  -- SQLAlchemy ORM, query helpers
    |
    v
notebooks/  -- analysis, figures (query DB, call src/ functions)
```

**Core rule:** No physics logic in notebooks. All model equations live in `src/physics.py`. Notebooks call `src` functions and never implement models inline.

### CLI Usage

```bash
python -m src.fit --baseline       # Phase 1: verify Paper 2 g(Rt) result
python -m src.fit --gate6          # Gate 6: constrained model on SPARC
python -m src.fit --gate6 --force  # Rerun, deleting existing constrained fits
```

### Standalone Scripts

- `scripts/power_assessment.py` — Gate 2-3: numerical anchor verification and THINGS power simulation

---

## 11. Analysis Phases

| Phase        | Description                                                                        | Status                           |
| ------------ | ---------------------------------------------------------------------------------- | -------------------------------- |
| Pre-analysis | Gates 1-6: lock plan, verify anchors, power simulation, overlap, constrained model | **COMPLETE**                     |
| 1            | SPARC baseline: reproduce Paper 2 g(R_t) result                                    | **COMPLETE** (N=98, 8.5% offset) |
| 2            | BTFR covariance test: fixed-slope residuals + signed-rank test                     | Not started                      |
| 3            | Constrained model BIC comparison on SPARC                                          | **COMPLETE** (not competitive)   |
| 4            | THINGS replication: galaxy-by-galaxy comparison with SPARC                         | Not started                      |
| 5            | LITTLE THINGS: dwarf-regime stress test                                            | Not started                      |
| 6            | PROBES: exploratory morphological/environmental stratification                     | Not started                      |
| 7            | Synthesis and manuscript                                                           | Not started                      |

---

## 12. Pre-Analysis Power Assessment

A power simulation (10,000 iterations, seed=42) was run before THINGS data acquisition using the Wilcoxon signed-rank test at p < 0.05 (two-sided). Two noise models were evaluated:

- **Conservative (raw scatter):** 0.409 dex (half 16th-84th percentile range of log10(g/a0_half))
- **Realistic (BTFR-residual scatter):** 0.277 dex (after fixed-slope trend removal)

At the expected THINGS resolved sample size (N ~ 20-24), the conservative model gives 53-63% power and the realistic model gives 84-93% power to detect a 0.2 dex departure from a0/2. THINGS is designated as the primary confirmatory dataset, with fallback to pilot if resolved N < 15.

**Post-acquisition update:** De Blok et al. (2008) published rotation curves for only 19 of the 34 THINGS galaxies, yielding an expected resolved N ~ 13 (assuming ~70% resolution rate). This falls below the pre-registered N = 15 threshold. We retain THINGS as the primary confirmatory dataset rather than triggering the fallback, because (a) extracting missing rotation curves would violate the independent-pipeline premise, and (b) elevating PROBES would violate Risk 9 (exploratory-only due to >= 0.2 dex M/L uncertainty). The power loss is accepted and will be explicitly documented in the manuscript. See `docs/internal/deviations_log.md` for the full deviation record.

**Implementation:** `scripts/power_assessment.py`

---

## 13. Unit Conversion

Accelerations appear in two unit systems throughout the codebase:

| Unit             | Symbol  | Value of a0  | Context                                            |
| ---------------- | ------- | ------------ | -------------------------------------------------- |
| km^2 s^-2 kpc^-1 | A0_MOND | 3703.0       | Internal computation (natural for rotation curves) |
| m s^-2           | --      | 1.2 x 10^-10 | Reporting and comparison to literature             |

Conversion: multiply km^2/s^2/kpc by 1e6 / 3.0857e19 to get m/s^2.

**Constants defined in:** `src/physics.py` (A0_MOND, A0_HALF, ACCEL_TO_MKS)

---

## 14. References

1. de Blok, W. J. G. et al. (2008). AJ, 136, 2648. "High-Resolution Rotation Curves and Galaxy Mass Models from THINGS."
2. Hunter, D. A. et al. (2012). AJ, 144, 134. "LITTLE THINGS."
3. Kass, R. E. & Raftery, A. E. (1995). JASA, 90, 773. "Bayes Factors."
4. Lelli, F., McGaugh, S. S., & Schombert, J. M. (2016). AJ, 152, 157. "SPARC: Mass Models for 175 Disk Galaxies."
5. Schneider, J. (2026a). Ap&SS (submitted). "A Rational Taper Model for Galaxy Rotation Curves."
6. Schneider, J. (2026b). (submitted). "Rational Taper Validation: A Four-Model Comparison Across 175 SPARC Galaxies."
7. Schwarz, G. (1978). Annals of Statistics, 6, 461. "Estimating the Dimension of a Model."
8. Walter, F. et al. (2008). AJ, 136, 2563. "THINGS: The HI Nearby Galaxy Survey."
