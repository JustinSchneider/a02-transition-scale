"""
Paper 3 Pre-Analysis Gate 2-3: Numerical Anchor Verification & Power Assessment

Computes the two Gate 2 anchors (median M_bar and g(Rt) scatter) from the Paper 2
resolved sample, then runs the Gate 3 power simulation for the THINGS signed-rank test.

Output is intended to be recorded in docs/Project3_Preregistration.md before any
THINGS data is acquired.

Usage:
    python scripts/power_assessment.py
"""

import numpy as np
from scipy.interpolate import interp1d
from scipy.stats import wilcoxon
import sqlite3
import os

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
Y_DISK = 0.5       # disk mass-to-light ratio (3.6 um)
Y_BULGE = 0.7      # bulge mass-to-light ratio (3.6 um)
HE_FACTOR = 1.33   # helium correction for HI gas mass
A0 = 1.2e-10       # MOND acceleration scale [m/s^2]
A0_HALF = A0 / 2   # alignment target [m/s^2]
KPC_TO_M = 3.0857e19  # 1 kpc in metres
ALPHA_BTFR = 0.238 # fixed BTFR slope from Paper 2

# Simulation parameters
N_SIM = 10_000
RNG_SEED = 42
SAMPLE_SIZES = [15, 20, 25, 34]
OFFSETS_DEX = [0.0, 0.1, 0.2]

# Paths (relative to project root)
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DB_PATH = os.path.join(BASE_DIR, "data", "galaxy_dynamics.db")
SPARC_PATH = os.path.join(BASE_DIR, "data", "raw", "SPARC_Lelli2016c.mrt")
BULGE_PATH = os.path.join(BASE_DIR, "data", "raw", "Bulges.mrt")


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------
def parse_sparc_table(path):
    """Parse SPARC MRT table for L[3.6] and MHI per galaxy."""
    sparc = {}
    header_done = False
    sep_count = 0
    with open(path) as f:
        for line in f:
            if line.startswith("---") or line.startswith("==="):
                sep_count += 1
                if sep_count >= 4:
                    header_done = True
                continue
            if not header_done:
                continue
            parts = line.split()
            if len(parts) < 17:
                continue
            try:
                name = parts[0]
                L36 = float(parts[7])   # 10^9 Lsun
                MHI = float(parts[13])  # 10^9 Msun
                sparc[name] = {"L36": L36, "MHI": MHI}
            except (ValueError, IndexError):
                continue
    return sparc


def parse_bulge_table(path):
    """Parse SPARC bulge luminosity table."""
    bulges = {}
    with open(path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) >= 2:
                try:
                    bulges[parts[0]] = float(parts[1])  # 10^9 Lsun
                except ValueError:
                    continue
    return bulges


def compute_mbar(sparc_entry, bulge_lum_1e9):
    """Total baryonic mass [Msun] from SPARC luminosity and gas mass."""
    L_total = sparc_entry["L36"] * 1e9   # Lsun
    L_bul = bulge_lum_1e9 * 1e9          # Lsun
    L_disk = L_total - L_bul
    MHI = sparc_entry["MHI"] * 1e9       # Msun
    return Y_DISK * L_disk + Y_BULGE * L_bul + HE_FACTOR * MHI


# ---------------------------------------------------------------------------
# Gate 2: Numerical anchors from Paper 2 resolved sample
# ---------------------------------------------------------------------------
def compute_anchors():
    """Compute median M_bar, beta, and scatter from the Paper 2 resolved sample."""
    sparc = parse_sparc_table(SPARC_PATH)
    bulges = parse_bulge_table(BULGE_PATH)

    conn = sqlite3.connect(DB_PATH)
    c = conn.cursor()

    # Resolved RT galaxies: Rt < R_max, converged
    c.execute("""
        SELECT m.galaxy_id, m.param1 AS omega, m.param2 AS Rt,
               MAX(r.radius_kpc) AS r_max
        FROM model_fits m
        JOIN radial_profiles r ON m.galaxy_id = r.galaxy_id
        WHERE m.model_name = 'rational_taper' AND m.converged = 1
        GROUP BY m.galaxy_id
        HAVING m.param2 < MAX(r.radius_kpc) AND m.param2 > 0
    """)
    resolved = c.fetchall()

    g_values = []
    m_bar_values = []

    for gid, omega, Rt, r_max in resolved:
        if gid not in sparc:
            continue

        c.execute(
            "SELECT radius_kpc, v_baryon_total "
            "FROM radial_profiles WHERE galaxy_id = ? ORDER BY radius_kpc",
            (gid,),
        )
        profile = c.fetchall()
        radii = np.array([p[0] for p in profile])
        v_bary = np.array([p[1] for p in profile])

        if len(radii) < 2 or Rt < radii.min() or Rt > radii.max():
            continue

        # g(Rt) = V_model(Rt)^2 / Rt
        v_bary_at_Rt = interp1d(radii, v_bary, kind="linear")(Rt)
        v_model_at_Rt = v_bary_at_Rt + omega * Rt / 2.0
        g_Rt_kpc = v_model_at_Rt ** 2 / Rt          # km^2/s^2/kpc
        g_Rt_mks = g_Rt_kpc * 1e6 / KPC_TO_M        # m/s^2

        M_bar = compute_mbar(sparc[gid], bulges.get(gid, 0.0))

        g_values.append(g_Rt_mks)
        m_bar_values.append(M_bar)

    conn.close()

    g_arr = np.array(g_values)
    m_bar_arr = np.array(m_bar_values)

    # Raw scatter: half the 16th-84th percentile range of log10(g / a0_half)
    log_resid_raw = np.log10(g_arr) - np.log10(A0_HALF)
    p16, p84 = np.percentile(log_resid_raw, [16, 84])
    raw_scatter = (p84 - p16) / 2.0

    # BTFR residual scatter
    M_bar_median = np.median(m_bar_arr)
    beta = np.log10(A0_HALF) - ALPHA_BTFR * np.log10(M_bar_median)
    log_g_trend = ALPHA_BTFR * np.log10(m_bar_arr) + beta
    btfr_resid = np.log10(g_arr) - log_g_trend
    bp16, bp84 = np.percentile(btfr_resid, [16, 84])
    btfr_scatter = (bp84 - bp16) / 2.0

    return {
        "n_resolved": len(g_arr),
        "median_g_Rt": np.median(g_arr),
        "median_m_bar": M_bar_median,
        "beta": beta,
        "raw_scatter": raw_scatter,
        "btfr_scatter": btfr_scatter,
    }


# ---------------------------------------------------------------------------
# Gate 3: Power simulation
# ---------------------------------------------------------------------------
def simulate_power(true_log_median_offset, n_galaxies, log_scatter, rng, n_sim=N_SIM):
    """
    Simulate power of a Wilcoxon signed-rank test.

    Parameters
    ----------
    true_log_median_offset : float
        True median offset in dex (0.0 = null hypothesis).
    n_galaxies : int
        Sample size.
    log_scatter : float
        Standard deviation of log-normal noise in dex.
    rng : np.random.Generator
    n_sim : int
        Number of Monte Carlo iterations.

    Returns
    -------
    float
        Rejection rate at p < 0.05 (two-sided).
    """
    rejections = 0
    for _ in range(n_sim):
        sample = rng.normal(true_log_median_offset, log_scatter, n_galaxies)
        _, p = wilcoxon(sample)
        if p < 0.05:
            rejections += 1
    return rejections / n_sim


def run_power_simulation(scatter_label, log_scatter, rng):
    """Run the full power grid for one scatter model."""
    print(f"\n{'=' * 65}")
    print(f"  Noise model: {scatter_label} = {log_scatter:.3f} dex")
    print(f"  {N_SIM:,} simulations per cell, seed={RNG_SEED}")
    print(f"{'=' * 65}")

    # Header
    header = f"{'N':>4}"
    for offset in OFFSETS_DEX:
        if offset == 0.0:
            header += f"  {'Type I (0.0)':>14}"
        else:
            header += f"  {'Power('+str(offset)+')':>14}"
    print(header)
    print("-" * len(header))

    results = {}
    for n in SAMPLE_SIZES:
        row = f"{n:>4}"
        results[n] = {}
        for offset in OFFSETS_DEX:
            power = simulate_power(offset, n, log_scatter, rng)
            results[n][offset] = power
            row += f"  {power:>14.3f}"
        print(row)

    return results


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    print("=" * 65)
    print("  Paper 3 Pre-Analysis: Gate 2-3")
    print("  Numerical Anchor Verification & THINGS Power Assessment")
    print("=" * 65)

    # --- Gate 2 ---
    anchors = compute_anchors()

    print(f"\n--- Gate 2: Numerical Anchors (N = {anchors['n_resolved']} resolved) ---\n")
    print(f"  Median g(Rt):          {anchors['median_g_Rt']:.4e} m/s^2")
    print(f"  a0/2:                  {A0_HALF:.4e} m/s^2")
    print(f"  Ratio (median / a0_2): {anchors['median_g_Rt'] / A0_HALF:.4f}")
    print()
    print(f"  Median M_bar:          {anchors['median_m_bar']:.4e} Msun")
    print(f"  log10(Median M_bar):   {np.log10(anchors['median_m_bar']):.4f}")
    print()
    print(f"  BTFR intercept (beta): {anchors['beta']:.4f}")
    print(f"    = log10(a0/2) - {ALPHA_BTFR} * log10(M_bar_median)")
    print(f"    = {np.log10(A0_HALF):.4f} - {ALPHA_BTFR} * "
          f"{np.log10(anchors['median_m_bar']):.4f}")
    print()
    print(f"  Raw g(Rt) scatter [half(p84-p16)]:  {anchors['raw_scatter']:.4f} dex")
    print(f"  BTFR residual scatter [half(p84-p16)]: {anchors['btfr_scatter']:.4f} dex")

    # --- Gate 3 ---
    print("\n\n--- Gate 3: THINGS Power Simulation ---")

    rng = np.random.default_rng(RNG_SEED)

    results_raw = run_power_simulation(
        "Raw scatter", anchors["raw_scatter"], rng
    )
    results_btfr = run_power_simulation(
        "BTFR-residual scatter", anchors["btfr_scatter"], rng
    )

    # --- Decision ---
    print("\n\n--- Decision ---\n")
    print("  Threshold: power >= 60% to detect 0.2 dex departure")
    print("  Expected THINGS resolved N ~ 15-25 (34 total, ~60-70% resolved)\n")

    for label, results in [("Raw scatter", results_raw),
                           ("BTFR-residual scatter", results_btfr)]:
        print(f"  [{label}]")
        for n in SAMPLE_SIZES:
            power_02 = results[n][0.2]
            status = "PRIMARY" if power_02 >= 0.60 else "PILOT"
            print(f"    N={n:>2}: Power(0.2 dex) = {power_02:.1%}  -->  "
                  f"THINGS = {status}")
        print()


if __name__ == "__main__":
    main()
