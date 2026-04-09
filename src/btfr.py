"""Baryonic Tully-Fisher Relation (BTFR) covariance test pipeline.

Tests whether the g(Rt) ~ a0/2 alignment is an artifact of the mass-dependent
BTFR trend.  The fixed-slope approach uses Paper 2's calibration (alpha=0.238,
beta=-12.55) so that no parameters are fit to the test sample.

Reused by NB02 (SPARC), NB09 (THINGS), and NB13 (LITTLE THINGS).
"""

import numpy as np
from scipy import stats

from src.ingest import load_sparc_metadata, load_bulge_luminosities, compute_mbar
from src.utils import setup_logger

logger = setup_logger(__name__)

# Locked pre-registration values (Paper 2 calibration)
BTFR_ALPHA = 0.238
BTFR_BETA = -12.55


def compute_btfr_residuals(
    g_Rt_mks: np.ndarray,
    m_bar_msun: np.ndarray,
    alpha: float = BTFR_ALPHA,
    beta: float = BTFR_BETA,
) -> np.ndarray:
    """Compute fixed-slope BTFR residuals in log space.

    residual = log10(g_Rt) - (alpha * log10(M_bar) + beta)

    Args:
        g_Rt_mks: Centripetal acceleration at Rt in m/s^2.
        m_bar_msun: Total baryonic mass in Msun.
        alpha: Fixed BTFR slope (default 0.238, locked from Paper 2).
        beta: Fixed BTFR intercept (default -12.55).

    Returns:
        Array of log10 residuals (only finite entries retained).

    Raises:
        ValueError: If no valid (finite, positive) entries remain.
    """
    g = np.asarray(g_Rt_mks, dtype=np.float64)
    m = np.asarray(m_bar_msun, dtype=np.float64)

    valid = (g > 0) & (m > 0) & np.isfinite(g) & np.isfinite(m)
    if not np.any(valid):
        raise ValueError("No valid (positive, finite) g_Rt or M_bar entries")

    log_g = np.log10(g[valid])
    log_m = np.log10(m[valid])
    log_g_trend = alpha * log_m + beta

    return log_g - log_g_trend


def run_btfr_covariance_test(
    g_Rt_mks: np.ndarray,
    m_bar_msun: np.ndarray,
    alpha: float = BTFR_ALPHA,
    beta: float = BTFR_BETA,
) -> dict:
    """Run the full BTFR covariance test on a resolved galaxy sample.

    Computes fixed-slope residuals, measures scatter reduction, and runs
    a Wilcoxon signed-rank test for median offset from zero.

    Scatter metric: half the 16th-84th percentile range (equivalent to
    1-sigma for a Gaussian distribution).

    Returns:
        Dict with keys: residuals, scatter_raw, scatter_residual,
        wilcoxon_stat, wilcoxon_pvalue, n_galaxies, median_residual.
    """
    g = np.asarray(g_Rt_mks, dtype=np.float64)
    m = np.asarray(m_bar_msun, dtype=np.float64)

    valid = (g > 0) & (m > 0) & np.isfinite(g) & np.isfinite(m)
    g_valid = g[valid]

    residuals = compute_btfr_residuals(g_valid, m[valid], alpha, beta)

    # Raw scatter: half(p84 - p16) of log10(g_Rt)
    log_g = np.log10(g_valid)
    p16_raw, p84_raw = np.percentile(log_g, [16, 84])
    scatter_raw = (p84_raw - p16_raw) / 2.0

    # Residual scatter: same metric after trend removal
    p16_res, p84_res = np.percentile(residuals, [16, 84])
    scatter_residual = (p84_res - p16_res) / 2.0

    # Wilcoxon signed-rank test: are residuals centered on zero?
    nonzero = residuals[residuals != 0]
    if len(nonzero) >= 10:
        stat, pvalue = stats.wilcoxon(nonzero, alternative="two-sided")
    else:
        logger.warning("Too few nonzero residuals (%d) for Wilcoxon test", len(nonzero))
        stat, pvalue = float("nan"), float("nan")

    return {
        "residuals": residuals,
        "scatter_raw": float(scatter_raw),
        "scatter_residual": float(scatter_residual),
        "wilcoxon_stat": float(stat),
        "wilcoxon_pvalue": float(pvalue),
        "n_galaxies": int(np.sum(valid)),
        "median_residual": float(np.median(residuals)),
    }


def compute_mbar_for_sample(
    galaxy_ids: list[str],
    sparc_mrt_path: str = None,
    bulge_path: str = None,
) -> dict[str, float]:
    """Compute M_bar for a list of galaxies from SPARC metadata.

    Loads SPARC metadata and bulge luminosities once, then calls
    compute_mbar() for each galaxy.

    Returns:
        Dict mapping galaxy_id -> M_bar in Msun.
        Galaxies not found in SPARC metadata are omitted with a warning.
    """
    sparc = load_sparc_metadata(sparc_mrt_path)
    bulges = load_bulge_luminosities(bulge_path)

    result = {}
    n_missing = 0

    for gid in galaxy_ids:
        if gid not in sparc:
            n_missing += 1
            continue
        meta = sparc[gid]
        l_bulge = bulges.get(gid, 0.0)
        result[gid] = compute_mbar(meta["L36"], meta["MHI"], l_bulge)

    if n_missing > 0:
        logger.warning(
            "%d of %d galaxies not found in SPARC metadata",
            n_missing, len(galaxy_ids),
        )

    return result
