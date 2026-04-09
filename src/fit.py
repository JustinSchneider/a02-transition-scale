"""Fitting orchestration for Paper 3.

Provides:
  - Phase 1 baseline reproduction (verify Paper 2 g(Rt) result from existing fits)
  - Gate 6 validation (constrained RT model on SPARC resolved sample)
  - CLI interface for running these steps

Usage:
    python -m src.fit --baseline       # Phase 1: verify Paper 2 g(Rt)
    python -m src.fit --gate6          # Gate 6: constrained model validation
    python -m src.fit --gate6 --force  # Rerun, deleting existing constrained fits
"""

import argparse
import sys

import numpy as np
import pandas as pd

from src.database import (
    get_engine, get_session, query_profiles_as_dataframe,
    query_fits_as_dataframe, insert_model_fit, delete_fits,
)
from src.physics import (
    A0_HALF, ACCEL_TO_MKS, compute_v_bary,
    compute_transition_diagnostics, fit_constrained_rt,
)
from src.utils import setup_logger

logger = setup_logger(__name__)


# ---------------------------------------------------------------------------
# Phase 1: Baseline reproduction
# ---------------------------------------------------------------------------

def run_phase1_baseline(db_path: str = None) -> pd.DataFrame:
    """Reproduce Paper 2's g(Rt) result from existing free RT fits.

    Reads the existing rational_taper fits from the database, computes g(Rt)
    for each resolved galaxy, and verifies N=98 and median ~ 6.51e-11 m/s^2.

    Returns DataFrame with one row per resolved galaxy.
    """
    engine = get_engine(db_path)
    session = get_session(engine)

    # Get all converged RT fits
    fits_df = query_fits_as_dataframe(session, model_name="rational_taper")
    fits_df = fits_df[fits_df["converged"] == True]

    results = []

    for _, fit_row in fits_df.iterrows():
        gid = fit_row["galaxy_id"]
        omega = fit_row["param1"]
        Rt = fit_row["param2"]

        if not (np.isfinite(omega) and np.isfinite(Rt) and Rt > 0):
            continue

        # Load profile
        prof = query_profiles_as_dataframe(session, gid)
        if prof.empty:
            continue

        radius = prof["radius_kpc"].values
        r_max = radius.max()

        # Spatially resolved check
        if Rt >= r_max:
            continue

        v_bary = prof["v_baryon_total"].values
        diag = compute_transition_diagnostics(radius, v_bary, omega, Rt)

        if not np.isfinite(diag["g_obs"]):
            continue

        g_mks = diag["g_obs"] * ACCEL_TO_MKS
        results.append({
            "galaxy_id": gid,
            "omega": omega,
            "Rt": Rt,
            "Rt_over_Rmax": Rt / r_max,
            "n_points": fit_row["n_points"],
            "g_Rt_kpc": diag["g_obs"],
            "g_Rt_mks": g_mks,
            "bic_free": fit_row["bic"],
        })

    session.close()
    df = pd.DataFrame(results)

    # Report
    a0_half_mks = A0_HALF * ACCEL_TO_MKS
    median_g = df["g_Rt_mks"].median()

    print("\n" + "=" * 60)
    print("  Phase 1: SPARC Baseline Reproduction")
    print("=" * 60)
    print(f"  Resolved galaxies (Rt < R_max): {len(df)}")
    print(f"  Median g(Rt): {median_g:.4e} m/s^2")
    print(f"  a0/2:         {a0_half_mks:.4e} m/s^2")
    print(f"  Ratio:        {median_g / a0_half_mks:.4f}")
    print(f"  Offset:       {(median_g / a0_half_mks - 1) * 100:+.1f}%")

    if len(df) == 98 and abs(median_g / a0_half_mks - 1.083) < 0.01:
        print("\n  BASELINE REPRODUCED: matches Paper 2 (N=98, ~8.3% offset)")
    else:
        print(f"\n  WARNING: expected N=98, got N={len(df)}")

    return df


# ---------------------------------------------------------------------------
# Gate 6: Constrained model validation
# ---------------------------------------------------------------------------

def run_gate6_validation(db_path: str = None, force: bool = False) -> pd.DataFrame:
    """Run the constrained RT model on the Paper 2 SPARC resolved sample.

    For each of the ~98 resolved galaxies:
      1. Load radial profile
      2. Fit constrained RT (g(Rt) = a0/2)
      3. Store result in database
      4. Compare BIC against free RT

    Returns DataFrame with per-galaxy results.
    """
    engine = get_engine(db_path)
    session = get_session(engine)

    # First get the resolved sample from the baseline
    fits_df = query_fits_as_dataframe(session, model_name="rational_taper")
    fits_df = fits_df[fits_df["converged"] == True]

    if force:
        n_deleted = delete_fits(session, model_name="constrained_rt")
        if n_deleted > 0:
            logger.info("Deleted %d existing constrained_rt fits", n_deleted)

    results = []
    n_no_solution = 0
    n_fitted = 0
    n_skipped = 0

    for _, fit_row in fits_df.iterrows():
        gid = fit_row["galaxy_id"]
        omega_free = fit_row["param1"]
        Rt_free = fit_row["param2"]

        if not (np.isfinite(omega_free) and np.isfinite(Rt_free) and Rt_free > 0):
            continue

        prof = query_profiles_as_dataframe(session, gid)
        if prof.empty:
            continue

        radius = prof["radius_kpc"].values
        r_max = radius.max()

        if Rt_free >= r_max:
            continue

        # Check if already fitted
        if not force:
            existing = query_fits_as_dataframe(session, galaxy_id=gid,
                                               model_name="constrained_rt")
            if not existing.empty:
                n_skipped += 1
                continue

        v_obs = prof["v_obs"].values
        v_err = prof["v_err"].values
        v_bary = prof["v_baryon_total"].values

        # Handle NaN v_err
        if np.any(pd.isna(v_err)):
            v_err = np.where(pd.isna(v_err), 1.0, v_err)

        # Fit constrained model
        result = fit_constrained_rt(
            radius, v_obs, v_err, v_bary,
            galaxy_id=gid,
        )

        # Store in database
        insert_model_fit(session, result.to_dict())

        if not result.converged:
            n_no_solution += 1

        n_fitted += 1

        row = {
            "galaxy_id": gid,
            "converged": result.converged,
            "omega_constrained": result.param1,
            "Rt_constrained": result.param2,
            "bic_constrained": result.bic,
            "bic_free": fit_row["bic"],
            "n_points": result.n_points,
            "Rt_free": Rt_free,
            "Rt_over_Rmax": Rt_free / r_max,
        }

        if result.converged and np.isfinite(result.bic) and np.isfinite(fit_row["bic"]):
            row["delta_bic"] = result.bic - fit_row["bic"]  # positive = free wins
        else:
            row["delta_bic"] = float("nan")

        results.append(row)

    session.close()
    df = pd.DataFrame(results)

    # Report
    print("\n" + "=" * 60)
    print("  Gate 6: Constrained Model Validation (g(Rt) = a0/2)")
    print("=" * 60)
    print(f"  Galaxies fitted:     {n_fitted}")
    print(f"  Skipped (existing):  {n_skipped}")
    print(f"  No valid Rt:         {n_no_solution}")
    print(f"  Converged:           {df['converged'].sum()}")

    if len(df) > 0 and df["converged"].any():
        converged = df[df["converged"]]
        print(f"\n  --- BIC Comparison (constrained - free) ---")
        print(f"  Median dBIC (all converged):  {converged['delta_bic'].median():+.2f}")
        print(f"  Mean dBIC:                    {converged['delta_bic'].mean():+.2f}")

        # High-quality subsample: n >= 20 AND Rt/R_max < 0.5
        hq = converged[(converged["n_points"] >= 20) & (converged["Rt_over_Rmax"] < 0.5)]
        if len(hq) > 0:
            print(f"\n  --- High-Quality Subsample (n>=20, Rt/R_max<0.5) ---")
            print(f"  N galaxies:          {len(hq)}")
            print(f"  Median dBIC:         {hq['delta_bic'].median():+.2f}")
            bic_competitive = abs(hq["delta_bic"].median()) < 2
            print(f"  |median dBIC| < 2:   {'YES' if bic_competitive else 'NO'}"
                  f" ({'competitive' if bic_competitive else 'not competitive'})")

        # Summary counts
        wins = (converged["delta_bic"] < 0).sum()
        losses = (converged["delta_bic"] > 0).sum()
        ties = (converged["delta_bic"] == 0).sum()
        print(f"\n  Constrained wins: {wins}  |  Free wins: {losses}  |  Ties: {ties}")

    return df


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description="Paper 3 fitting pipeline")
    parser.add_argument("--baseline", action="store_true",
                        help="Phase 1: reproduce Paper 2 g(Rt) baseline")
    parser.add_argument("--gate6", action="store_true",
                        help="Gate 6: run constrained RT on SPARC resolved sample")
    parser.add_argument("--force", action="store_true",
                        help="Delete and refit existing results")
    parser.add_argument("--db", type=str, default=None,
                        help="Path to galaxy_dynamics.db")

    args = parser.parse_args()

    if not (args.baseline or args.gate6):
        parser.print_help()
        sys.exit(1)

    if args.baseline:
        run_phase1_baseline(args.db)

    if args.gate6:
        run_gate6_validation(args.db, force=args.force)


if __name__ == "__main__":
    main()
