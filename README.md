# The a₀/2 Transition Scale: Empirical Validation Across Independent Rotation Curve Datasets

**Paper 3 of the Rational Taper Series**

Justin Schneider (independent researcher)

---

## Summary

Paper 2 of the Rational Taper (RT) series reported a post-hoc finding: the median total centripetal acceleration at the RT taper radius across 98 spatially resolved SPARC galaxies is 6.51 × 10⁻¹¹ m s⁻², within 8.5% of the MOND acceleration scale a₀/2. This alignment was not predicted by the model.

This repository contains the analysis for Paper 3, which tests whether that alignment survives in independent datasets and controls for baryonic Tully-Fisher relation (BTFR) covariance. The analysis plan was [pre-registered](docs/Project3_Preregistration.md) before acquiring any new data.

**Result:** The a₀/2 alignment is a demographic artifact of the SPARC sample's mass distribution. The BTFR trend crosses a₀/2 at log₁₀(M_bar) ≈ 9.78, which coincides with SPARC's median baryonic mass. Independent datasets at different mass scales show systematic deviations: THINGS (+45.6%) and LITTLE THINGS (−72.4%) both fail the pre-registered ±10% alignment threshold.

## The Rational Taper Model

```
V_model(R) = V_bary(R) + (ω · R) / (1 + R/Rₜ)
```

- **ω** (omega): kinematic correction rate [km/s/kpc]
- **Rₜ**: taper (transition) radius [kpc]
- **g(Rₜ) = V²_total(Rₜ) / Rₜ**: total centripetal acceleration at the taper radius

The RT model is introduced in Paper 1 (Schneider 2026a) and benchmarked in Paper 2 (Schneider 2026b).

## Datasets

| Dataset | N | Role | Reference |
|---|---|---|---|
| SPARC | 175 | Baseline anchor | Lelli et al. (2016) |
| THINGS | 19 | Primary confirmatory | de Blok et al. (2008) |
| LITTLE THINGS | 41 | Secondary confirmatory | Oh et al. (2015) |

## Repository Structure

```
data/raw/              SPARC, THINGS, and LITTLE THINGS public data files
data/supplemental/     Paper 2 output files (Tournament_Results.csv, etc.)
data/galaxy_dynamics.db   SQLite database with all profiles and model fits
docs/                  Pre-registration, proposal, and deviations log
notebooks/             Jupyter notebooks NB01–NB18 (analysis pipeline)
scripts/               Standalone scripts (power assessment, data fetching)
src/                   Analysis pipeline code
tests/                 Unit and integration tests
```

## Analysis Pipeline

The analysis proceeds in six blocks of Jupyter notebooks:

| Block | Notebooks | Purpose |
|---|---|---|
| 1 | NB01–NB03 | SPARC baseline: reproduce Paper 2 result, BTFR covariance test, constrained model |
| 2 | NB04–NB06 | THINGS ingestion, overlap pairing, baryonic decomposition |
| 3 | NB07–NB10 | THINGS RT fitting, cross-validation, primary confirmatory tests |
| 4 | NB11–NB13 | LITTLE THINGS dwarf-regime stress test |
| 5 | NB14–NB16 | PROBES (skipped; null result was decisive — notebooks not created) |
| 6 | NB17–NB18 | Cross-dataset synthesis and publication figures |

## Source Modules

| Module | Purpose |
|---|---|
| `src/physics.py` | RT model equations, fitting, transition diagnostics |
| `src/database.py` | SQLAlchemy ORM and queries |
| `src/ingest.py` | Data loading from files and database |
| `src/btfr.py` | BTFR covariance analysis |
| `src/fit.py` | Fitting orchestration and CLI |
| `src/utils.py` | Logging and path utilities |

## Requirements

- Python 3.9+
- numpy, scipy, pandas, matplotlib, sqlalchemy

Install dependencies:

```bash
pip install -r requirements.txt
```

## Reproducing the Analysis

Run notebooks NB01 through NB18 in order. Each notebook is self-contained and documents its inputs, methods, and outputs.

The fitting CLI can also be used directly:

```bash
python -m src.fit --baseline    # Phase 1: SPARC baseline
python -m src.fit --gate6       # Gate 6: constrained model validation
```

Run tests:

```bash
pytest tests/
```

## Key Results

- **SPARC baseline** (N=98 resolved): median g(Rₜ) = 6.51 × 10⁻¹¹ m s⁻², 8.3% offset from a₀/2
- **THINGS** (N=8 valid): median g(Rₜ) = 8.73 × 10⁻¹¹ m s⁻², +45.6% from a₀/2
- **LITTLE THINGS** (N=14 valid): median g(Rₜ) = 1.66 × 10⁻¹¹ m s⁻², −72.4% from a₀/2
- **Pooled** (N=22): median BTFR residual = −0.082 dex, Wilcoxon p = 0.166
- **BTFR crossing**: g(Rₜ) trend crosses a₀/2 at log₁₀(M_bar) = 9.78, coinciding with SPARC's median mass
- **Conclusion**: The a₀/2 alignment is mass-dependent and driven by SPARC demographics, not a universal kinematic scale

## Pre-registration

The full analysis plan, including all definitions, thresholds, and statistical tests, was locked before acquiring THINGS data. See [docs/Project3_Preregistration.md](docs/Project3_Preregistration.md). Deviations from the pre-registered plan are documented in [docs/deviations_log.md](docs/deviations_log.md).

## Related Work

- **Paper 1** — Introduces the Rational Taper model (Schneider 2026a, under review at *Astrophysics and Space Science*)
- **Paper 2** — Benchmarks RT against NFW and MOND across 175 SPARC galaxies (Schneider 2026b)

## License

This work is dedicated to the public domain under [CC0 1.0 Universal](LICENSE). Please cite appropriately if you use this work.
