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
   analysis â€” any a0/2 alignment could be attributed to M/L systematics rather than
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
system — introducing substantial systematic uncertainty.

### Justification

The exclusion is based on documented morphological disturbance and data availability,
not on any knowledge of the RT fit outcome. This decision was committed in NB05
before NB06 was written.

### Impact

THINGS non-overlap sample reduced from 5 to 4 galaxies (combined with Deviation 2).
Total THINGS sample: 13 overlap + 4 non-overlap = 17 galaxies.
