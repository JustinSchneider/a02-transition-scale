# Pre-Registration Deviations Log

This file records all deviations from the pre-registered analysis plan
(`docs/Project3_Preregistration.md`, commit dc24e10, 09 APR 2026).

---

## Deviation 1: Primary Dataset Fallback Override

**Date:** 09 APR 2026
**Pre-registered rule:** If the THINGS resolved sample N < 15, demote THINGS to pilot
status and consider elevating PROBES.

### Description

The expected resolved sample size for THINGS has dropped to N ~ 13. De Blok et al.
(2008, AJ 136, 2648) published high-resolution rotation curves for only 19 of the
34 THINGS survey galaxies (Walter et al. 2008, AJ 136, 2563). Applying the ~70%
spatial resolution rate observed in SPARC yields an expected N_resolved ~ 13, which
falls below the pre-registered N = 15 threshold.

### Action taken

We are consciously overriding the fallback plan. THINGS is retained as the primary
confirmatory dataset despite the reduced sample size.

### Justification

Two options were considered:

2. **Trigger the fallback: demote THINGS to pilot, elevate PROBES to primary.**
   Rejected: PROBES uses optical spectroscopy with M/L uncertainty >= 0.2 dex (Risk 9
   in the risk register). Elevating it to a confirmatory role would undermine the
   analysis — any a0/2 alignment could be attributed to M/L systematics rather than
   genuine kinematics.

3. **Accept the power loss on a clean, independent dataset.**
   Adopted: This is the most scientifically rigorous path. The power reduction will be explicitly reported in the manuscript. If the THINGS result is inconclusive due to sample size, this is itself an honest result that future surveys (e.g., WALLABY) can address with larger samples.

### Impact on analysis

- Power at N=13 (conservative scatter 0.409 dex): ~45% for 0.2 dex effect
- Power at N=13 (realistic BTFR-residual scatter 0.277 dex): ~70% for 0.2 dex effect
- The Wilcoxon signed-rank test remains appropriate at N=13
- LITTLE THINGS (Block 4) provides an independent secondary test unaffected by this deviation


---

## Deviation 2: NGC 4826 Excluded from THINGS Sample

**Date:** 09 APR 2026
**Pre-registered rule:** Include all THINGS galaxies with i > 30 deg that pass
morphological QC.

### Description

NGC 4826 (M64, the "Evil Eye Galaxy") is excluded from the THINGS analysis sample
due to its counter-rotating gas disk. The de Blok et al. (2008) rotation curve shows
negative V_rot for the inner ~50 arcsec, an 80-arcsec data gap, then positive V_rot
from 130 arcsec outward. The RT model assumes a single kinematic component and
cannot meaningfully fit a counter-rotating system.

Additionally, NGC 4826 lacks a de Blok mass model file, so no baryonic decomposition
is available without independent derivation.

### Justification

The exclusion is based on pre-existing kinematic complexity (documented in the
literature) and data availability, not on any knowledge of the RT fit outcome.
This decision was committed in NB05 before NB06 (baryonic decomposition) was written.

### Impact

THINGS non-overlap sample reduced from 6 to 5 galaxies. Combined with Deviation 3
(NGC 3627), the non-overlap sample is 4 galaxies.

---

## Deviation 3: NGC 3627 Excluded from THINGS Sample

**Date:** 09 APR 2026
**Pre-registered rule:** Include all THINGS galaxies with i > 30 deg that pass
morphological QC.

### Description

NGC 3627 (M66) is excluded from the THINGS analysis sample due to its disturbed
morphology as a member of the Leo Triplet (interacting with NGC 3623 and NGC 3628).
The tidal interaction compromises the assumption of an equilibrium rotation curve
that the RT model requires.

Additionally, NGC 3627 lacks a de Blok mass model file, so baryonic decomposition
would require independent derivation from S4G photometry for a tidally perturbed
system � introducing substantial systematic uncertainty.

### Justification

The exclusion is based on documented morphological disturbance and data availability,
not on any knowledge of the RT fit outcome. This decision was committed in NB05
before NB06 was written.

### Impact

THINGS non-overlap sample reduced from 5 to 4 galaxies (combined with Deviation 2).
Total THINGS sample: 13 overlap + 4 non-overlap = 17 galaxies.

---

## Deviation 4: Initial Guess Bounds Fix in fit_rational_taper()

**Date:** 09 APR 2026
**Pre-registered rule:** N/A (code implementation detail, not analysis plan)

### Description

The multi-start optimizer in `src/physics.py::fit_rational_taper()` uses four
hard-coded initial guesses: `(5, 5), (10, 2), (2, 15), (20, R_max)`. For galaxies
with small R_max (specifically NGC 2976 with R_max = 2.25 kpc), the third starting
point `(omega=2, Rt=15)` exceeds the upper bound `Rt_upper = 5 * R_max = 11.24`.
This caused `scipy.optimize.curve_fit` to raise a `ValueError` that was not caught
(only `RuntimeError` was caught).

### Action taken

Two changes applied to `fit_rational_taper()`:
1. Initial guesses are now clamped to within the parameter bounds before being
   passed to `curve_fit`.
2. The exception handler now catches both `RuntimeError` and `ValueError`.

### Impact

The fix was applied before THINGS fitting (NB07) but after the SPARC baseline fits
(NB01, produced in Paper 2). The SPARC fits in the database were produced by the
original code, which never encountered this bug because all SPARC galaxies have
R_max large enough that the initial guesses fell within bounds. The fix is
therefore backwards-compatible: no SPARC results are affected. The THINGS fits use
the corrected code, but the correction only changes which starting points are
attempted — the optimizer converges to the same minimum if any starting point
reaches it.

---

## Deviation 5: Two Resolved THINGS Galaxies Have Non-Computable g(Rt)

**Date:** 09 APR 2026
**Pre-registered rule:** Compute g(Rt) for all spatially resolved galaxies
(Rt < R_max).

### Description

Two THINGS galaxies are formally spatially resolved (Rt < R_max) but have fitted
Rt values below their innermost radial data point (R_min), making baryonic velocity
interpolation at Rt impossible:

- **NGC 3031 (M81):** Rt = 1.53 kpc, R_min = 2.43 kpc. The large inner data gap
  is a known property of the de Blok (2008) rotation curve for this galaxy (the
  inner ~2.4 kpc is excluded due to bar contamination in the velocity field).
- **NGC 4736 (M94):** Rt = 0.22 kpc, R_min = 0.45 kpc. The fit converged with
  omega at the upper bound (200 km/s/kpc), indicating a degenerate solution where
  the taper correction is concentrated at very small radii. This is a non-physical
  fit artifact rather than a meaningful kinematic measurement.

### Action taken

Both galaxies are excluded from the g(Rt) distribution analysis. They are retained
in the database and reported in NB07 as "resolved with NaN g(Rt)" for transparency.

### Impact

The usable THINGS g(Rt) sample is reduced from 10 resolved galaxies to 8 with
valid g(Rt) values. Combined with Deviation 1 (expected N ~ 13 resolved), the
effective sample is smaller than anticipated. This further reduces statistical
power for the primary Wilcoxon test in NB09.
