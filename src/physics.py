"""Core physics equations for rotation curve decomposition and model fitting.

Paper 3 models:
  - Rational Taper (free):       2 parameters (omega, R_t)
  - Rational Taper (constrained): 1 parameter  (omega; R_t derived from g(Rt)=a0/2)

Baryonic velocity convention from Lelli et al. (2016):
  V_bary = sqrt(|V_gas|*V_gas + Upsilon_d*|V_disk|*V_disk + Upsilon_b*|V_bulge|*V_bulge)

Key references:
  - Lelli, McGaugh, & Schombert (2016) AJ 152, 157  -- SPARC data and Eq. 2
  - Schneider (2026a) Ap&SS (submitted)              -- RT model introduction
  - Schneider (2026b) (submitted)                    -- RT validation, g(Rt) ~ a0/2
"""

from dataclasses import dataclass, field
from typing import Optional

import numpy as np
from scipy.optimize import brentq, minimize_scalar, curve_fit

from src.utils import setup_logger

logger = setup_logger(__name__)

# ---------------------------------------------------------------------------
# Physical constants
# ---------------------------------------------------------------------------

# MOND acceleration scale: a_0 = 1.2e-10 m/s^2 in rotation curve units
# a_0 [km^2/s^2/kpc] = 1.2e-10 * 3.0857e19 / 1e6 = 3703
A0_MOND = 3703.0  # km^2 s^-2 kpc^-1

# a_0/2 -- the alignment target from Paper 2
A0_HALF = A0_MOND / 2.0  # 1851.5 km^2 s^-2 kpc^-1

# Conversion factor: km^2/s^2/kpc -> m/s^2
KPC_TO_M = 3.0857e19  # 1 kpc in metres
ACCEL_TO_MKS = 1e6 / KPC_TO_M  # multiply kpc-units by this to get m/s^2


# ---------------------------------------------------------------------------
# Baryonic velocity computation
# ---------------------------------------------------------------------------


def compute_v_bary(
    v_gas: np.ndarray,
    v_disk: np.ndarray,
    v_bulge: np.ndarray,
    upsilon_disk: float = 0.5,
    upsilon_bulge: float = 0.7,
    gas_scale: float = 1.0,
) -> np.ndarray:
    """Compute total baryonic velocity using Lelli et al. (2016) Eq. 2 convention.

    V_bary(R) = sqrt(|V_gas|*V_gas + Y_d * |V_disk|*V_disk + Y_b * |V_bulge|*V_bulge)

    The |V|*V pattern preserves sign: if V is negative (e.g., gas central
    depression), its squared contribution is negative, reducing V_bary.

    SPARC provides V_disk and V_bulge at Upsilon=1, so Y enters as a direct
    multiplier on the V^2 term.

    Returns:
        Total baryonic velocity (km/s). Always non-negative.
    """
    v_gas = np.asarray(v_gas, dtype=np.float64)
    v_disk = np.asarray(v_disk, dtype=np.float64)
    v_bulge = np.asarray(v_bulge, dtype=np.float64)

    v2_gas = gas_scale * np.abs(v_gas) * v_gas
    v2_disk = upsilon_disk * np.abs(v_disk) * v_disk
    v2_bulge = upsilon_bulge * np.abs(v_bulge) * v_bulge

    v2_total = v2_gas + v2_disk + v2_bulge
    return np.sqrt(np.maximum(0.0, v2_total))


# ---------------------------------------------------------------------------
# BIC
# ---------------------------------------------------------------------------


def compute_bic(n_points: int, k_params: int, chi_squared: float) -> float:
    """Compute the Bayesian Information Criterion: BIC = chi^2 + k * ln(n)."""
    return chi_squared + k_params * np.log(n_points)


# ---------------------------------------------------------------------------
# Shared error-handling helper
# ---------------------------------------------------------------------------


def _safe_errors(v_err: np.ndarray, galaxy_id: str = "unknown") -> np.ndarray:
    """Replace zero/negative errors with the minimum nonzero error (floor 1.0 km/s)."""
    positive_errs = v_err[v_err > 0]
    min_err = float(np.min(positive_errs)) if len(positive_errs) > 0 else 1.0
    v_err_safe = np.where(v_err > 0, v_err, min_err)
    if np.any(v_err <= 0):
        logger.warning(
            "%s: replaced %d zero/negative errors with %.2f km/s",
            galaxy_id, int(np.sum(v_err <= 0)), min_err,
        )
    return v_err_safe


# ---------------------------------------------------------------------------
# Result dataclass
# ---------------------------------------------------------------------------


@dataclass
class ModelFitResult:
    """Container for results from any single-model fit."""

    galaxy_id: str
    model_name: str
    n_params: int
    chi_squared: float
    reduced_chi_squared: float
    bic: float
    residuals_rmse: float
    n_points: int
    converged: bool
    flag_v_obs_lt_v_bary: bool
    method_version: str
    upsilon_disk: float
    upsilon_bulge: float
    v_bary: np.ndarray = field(repr=False)
    v_model: np.ndarray = field(repr=False)
    residuals: np.ndarray = field(repr=False)
    param1: Optional[float] = None
    param1_err: Optional[float] = None
    param2: Optional[float] = None
    param2_err: Optional[float] = None

    def to_dict(self) -> dict:
        """Convert to dict suitable for database insertion (excludes arrays)."""
        return {
            "galaxy_id": self.galaxy_id,
            "model_name": self.model_name,
            "n_params": self.n_params,
            "chi_squared": self.chi_squared,
            "reduced_chi_squared": self.reduced_chi_squared,
            "bic": self.bic,
            "residuals_rmse": self.residuals_rmse,
            "n_points": self.n_points,
            "converged": self.converged,
            "flag_v_obs_lt_v_bary": self.flag_v_obs_lt_v_bary,
            "method_version": self.method_version,
            "upsilon_disk": self.upsilon_disk,
            "upsilon_bulge": self.upsilon_bulge,
            "param1": self.param1,
            "param1_err": self.param1_err,
            "param2": self.param2,
            "param2_err": self.param2_err,
        }


# ---------------------------------------------------------------------------
# Baryonic interpolation (signed-square convention)
# ---------------------------------------------------------------------------


def interpolate_v_bary(
    radius_kpc: np.ndarray,
    v_baryon_total: np.ndarray,
    r_query: float,
) -> float:
    """Safely interpolate baryonic velocity at an arbitrary radius.

    Uses signed-square convention: interpolates V^2_bary (with sign) then
    recovers velocity as sign(x) * sqrt(|x|). Avoids imaginary values when
    the profile crosses zero (e.g., inner gas depressions).

    Returns NaN if r_query is outside the profile range.
    """
    if r_query < radius_kpc[0] or r_query > radius_kpc[-1]:
        return float("nan")
    v2_signed = np.sign(v_baryon_total) * v_baryon_total ** 2
    v2_at_r = float(np.interp(r_query, radius_kpc, v2_signed))
    return float(np.sign(v2_at_r) * np.sqrt(abs(v2_at_r)))


# ---------------------------------------------------------------------------
# Free Rational Taper model (k=2)
# ---------------------------------------------------------------------------


def rt_model_velocity(
    radius: np.ndarray,
    v_bary: np.ndarray,
    omega: float,
    r_t: float,
) -> np.ndarray:
    """Compute RT model velocity: V_model = V_bary + omega*R / (1 + R/Rt)."""
    return v_bary + omega * radius / (1.0 + radius / r_t)


def fit_rational_taper(
    radius: np.ndarray,
    v_obs: np.ndarray,
    v_err: np.ndarray,
    v_bary: np.ndarray,
    galaxy_id: str = "unknown",
    method_version: str = "v1_rational_taper",
    upsilon_disk: float = 0.5,
    upsilon_bulge: float = 0.7,
    omega_bounds: tuple = (0.0, 200.0),
    rt_bounds: tuple = (0.1, None),
) -> ModelFitResult:
    """Fit the free Rational Taper model (2 parameters: omega, Rt).

    Multi-start optimizer with 4 initial conditions; retains lowest chi^2.
    Returns ModelFitResult with param1=omega, param2=Rt.
    """
    radius = np.asarray(radius, dtype=np.float64)
    v_obs = np.asarray(v_obs, dtype=np.float64)
    v_err = np.asarray(v_err, dtype=np.float64)
    v_bary = np.asarray(v_bary, dtype=np.float64)

    n_points = len(radius)
    v_err_safe = _safe_errors(v_err, galaxy_id)
    flag_v_obs_lt_v_bary = bool(np.any(v_obs < v_bary))

    r_max = float(np.max(radius))
    rt_upper = rt_bounds[1] if rt_bounds[1] is not None else 5.0 * r_max

    def _model(r, omega, r_t):
        return v_bary + omega * r / (1.0 + r / r_t)

    initial_guesses = [(5.0, 5.0), (10.0, 2.0), (2.0, 15.0), (20.0, r_max)]

    # Clamp initial guesses to bounds
    lo = [omega_bounds[0], rt_bounds[0]]
    hi = [omega_bounds[1], rt_upper]
    clamped_guesses = []
    for p0 in initial_guesses:
        p0_clamped = [
            max(lo[i] + 1e-6, min(p0[i], hi[i] - 1e-6))
            for i in range(2)
        ]
        clamped_guesses.append(tuple(p0_clamped))

    best_chi2 = np.inf
    best_popt = None
    best_pcov = None
    converged = False

    for p0 in clamped_guesses:
        try:
            popt, pcov = curve_fit(
                _model, radius, v_obs, p0=list(p0), sigma=v_err_safe,
                absolute_sigma=True,
                bounds=([omega_bounds[0], rt_bounds[0]], [omega_bounds[1], rt_upper]),
                maxfev=5000,
            )
            v_trial = _model(radius, *popt)
            res = v_obs - v_trial
            chi2_trial = float(np.sum((res / v_err_safe) ** 2))
            if chi2_trial < best_chi2:
                best_chi2 = chi2_trial
                best_popt = popt
                best_pcov = pcov
                converged = True
        except (RuntimeError, ValueError):
            pass

    if not converged:
        logger.warning("%s [RT free]: no convergence", galaxy_id)
        nan = float("nan")
        return ModelFitResult(
            galaxy_id=galaxy_id, model_name="rational_taper", n_params=2,
            chi_squared=nan, reduced_chi_squared=nan, bic=nan,
            residuals_rmse=nan, n_points=n_points, converged=False,
            flag_v_obs_lt_v_bary=flag_v_obs_lt_v_bary,
            method_version=method_version, upsilon_disk=upsilon_disk,
            upsilon_bulge=upsilon_bulge,
            v_bary=v_bary, v_model=np.full_like(radius, nan),
            residuals=np.full_like(radius, nan),
            param1=nan, param1_err=nan, param2=nan, param2_err=nan,
        )

    omega_best = float(best_popt[0])
    rt_best = float(best_popt[1])
    perr = np.sqrt(np.diag(best_pcov))
    v_model = _model(radius, omega_best, rt_best)
    residuals = v_obs - v_model
    chi2 = best_chi2
    dof = max(n_points - 2, 1)

    bic = compute_bic(n_points, 2, chi2)

    return ModelFitResult(
        galaxy_id=galaxy_id, model_name="rational_taper", n_params=2,
        chi_squared=chi2, reduced_chi_squared=chi2 / dof, bic=bic,
        residuals_rmse=float(np.sqrt(np.mean(residuals ** 2))),
        n_points=n_points, converged=True,
        flag_v_obs_lt_v_bary=flag_v_obs_lt_v_bary,
        method_version=method_version, upsilon_disk=upsilon_disk,
        upsilon_bulge=upsilon_bulge,
        v_bary=v_bary, v_model=v_model, residuals=residuals,
        param1=omega_best, param1_err=float(perr[0]),
        param2=rt_best, param2_err=float(perr[1]),
    )


# ---------------------------------------------------------------------------
# Transition diagnostics
# ---------------------------------------------------------------------------


def compute_transition_diagnostics(
    radius_kpc: np.ndarray,
    v_bary_profile: np.ndarray,
    omega: float,
    R_t: float,
) -> dict:
    """Compute physical quantities at the RT transition radius Rt.

    At R = Rt the taper correction = omega*Rt/2 (half of V_sat).

    Returns dict with keys:
        v_bary_rt, v_corr_rt, v_total_rt, g_obs, g_bary,
        eta_additive, eta_quadrature.
    All values NaN if Rt <= 0 or interpolation fails.
    """
    _nan = {
        "v_bary_rt": float("nan"), "v_corr_rt": float("nan"),
        "v_total_rt": float("nan"), "g_obs": float("nan"),
        "g_bary": float("nan"), "eta_additive": float("nan"),
        "eta_quadrature": float("nan"),
    }

    if not (np.isfinite(omega) and np.isfinite(R_t) and R_t > 0):
        return _nan

    v_bary_rt = interpolate_v_bary(radius_kpc, v_bary_profile, R_t)
    if not np.isfinite(v_bary_rt):
        return _nan

    v_corr_rt = omega * R_t / 2.0
    v_total_rt = v_bary_rt + v_corr_rt

    g_obs = v_total_rt ** 2 / R_t
    g_bary = v_bary_rt ** 2 / R_t if v_bary_rt != 0 else float("nan")

    if not np.isfinite(g_bary) or g_bary == 0:
        return {**_nan, "v_bary_rt": v_bary_rt, "v_corr_rt": v_corr_rt,
                "v_total_rt": v_total_rt, "g_obs": g_obs}

    eta_additive = g_obs / g_bary
    eta_quadrature = (v_bary_rt ** 2 + v_corr_rt ** 2) / R_t / g_bary

    return {
        "v_bary_rt": v_bary_rt, "v_corr_rt": v_corr_rt,
        "v_total_rt": v_total_rt, "g_obs": g_obs, "g_bary": g_bary,
        "eta_additive": eta_additive, "eta_quadrature": eta_quadrature,
    }


# ---------------------------------------------------------------------------
# Constrained RT model (k=1): g(Rt) = a0/2
# ---------------------------------------------------------------------------


def find_constrained_rt(
    omega: float,
    radius_kpc: np.ndarray,
    v_bary_profile: np.ndarray,
    a0_half: float = A0_HALF,
    n_scan: int = 200,
) -> tuple[float, int]:
    """For a given omega, find Rt such that g(Rt) = a0/2.

    The constraint equation is:
        [V_bary(Rt) + omega*Rt/2]^2 / Rt = a0/2

    This is transcendental in Rt because V_bary(Rt) depends on Rt through
    the rotation curve profile.

    Strategy:
        1. Scan the constraint function on a fine grid across the data range
        2. Identify all sign changes (potential roots)
        3. Use brentq on each bracketed interval
        4. Return the smallest positive root (physically: innermost transition)

    Args:
        omega: Kinematic correction rate (km/s/kpc).
        radius_kpc: Radial profile positions (kpc), sorted ascending.
        v_bary_profile: Baryonic velocity at each radius (km/s).
        a0_half: Target acceleration in km^2/s^2/kpc. Default A0_HALF.
        n_scan: Number of grid points for sign-change scan.

    Returns:
        (Rt, n_roots): Rt in kpc (or np.nan if no solution), and the
        number of roots found in the data range.
    """
    r_min = float(radius_kpc[0])
    r_max = float(radius_kpc[-1])

    # Small inset to avoid edge interpolation issues
    r_lo = r_min * 1.01 if r_min > 0 else r_min + 0.01
    r_hi = r_max * 0.99

    if r_lo >= r_hi:
        return float("nan"), 0

    def constraint_eq(Rt):
        v_bary_at_Rt = interpolate_v_bary(radius_kpc, v_bary_profile, Rt)
        if not np.isfinite(v_bary_at_Rt):
            return float("nan")
        v_total = v_bary_at_Rt + omega * Rt / 2.0
        return v_total ** 2 / Rt - a0_half

    # Scan for sign changes
    r_grid = np.linspace(r_lo, r_hi, n_scan)
    f_grid = np.array([constraint_eq(r) for r in r_grid])

    # Find all sign changes
    roots = []
    for i in range(len(f_grid) - 1):
        if np.isfinite(f_grid[i]) and np.isfinite(f_grid[i + 1]):
            if f_grid[i] * f_grid[i + 1] < 0:
                try:
                    root = brentq(constraint_eq, r_grid[i], r_grid[i + 1])
                    roots.append(root)
                except (ValueError, RuntimeError):
                    pass

    n_roots = len(roots)
    if n_roots == 0:
        return float("nan"), 0

    # Return smallest positive root (innermost transition)
    return min(roots), n_roots


def fit_constrained_rt(
    radius: np.ndarray,
    v_obs: np.ndarray,
    v_err: np.ndarray,
    v_bary: np.ndarray,
    galaxy_id: str = "unknown",
    method_version: str = "v1_constrained",
    upsilon_disk: float = 0.5,
    upsilon_bulge: float = 0.7,
    omega_bounds: tuple = (0.01, 200.0),
    a0_half: float = A0_HALF,
) -> ModelFitResult:
    """Fit the constrained RT model: g(Rt) = a0/2 (1 free parameter: omega).

    For each candidate omega, Rt is determined by numerically solving the
    constraint equation [V_bary(Rt) + omega*Rt/2]^2/Rt = a0/2. The model
    velocity is then V_model = V_bary + omega*R/(1+R/Rt), and chi^2 is
    computed against v_obs. The optimizer finds the omega that minimizes chi^2.

    BIC uses k=1 (only omega is free; Rt is derived from the constraint).

    Args:
        radius: Radial positions (kpc).
        v_obs: Observed velocities (km/s).
        v_err: Velocity errors (km/s).
        v_bary: Pre-computed baryonic velocity (km/s).
        galaxy_id: Galaxy identifier.
        omega_bounds: (lower, upper) bounds for omega.
        a0_half: Target acceleration for the constraint (km^2/s^2/kpc).

    Returns:
        ModelFitResult with model_name="constrained_rt", n_params=1,
        param1=omega, param2=Rt (derived).
    """
    radius = np.asarray(radius, dtype=np.float64)
    v_obs = np.asarray(v_obs, dtype=np.float64)
    v_err = np.asarray(v_err, dtype=np.float64)
    v_bary = np.asarray(v_bary, dtype=np.float64)

    n_points = len(radius)
    v_err_safe = _safe_errors(v_err, galaxy_id)
    flag_v_obs_lt_v_bary = bool(np.any(v_obs < v_bary))

    def chi2_for_omega(omega):
        Rt, _ = find_constrained_rt(omega, radius, v_bary, a0_half)
        if not np.isfinite(Rt) or Rt <= 0:
            return 1e30  # No valid Rt -> reject this omega
        v_model = v_bary + omega * radius / (1.0 + radius / Rt)
        res = (v_obs - v_model) / v_err_safe
        return float(np.sum(res ** 2))

    # Optimize over omega
    result = minimize_scalar(
        chi2_for_omega,
        bounds=omega_bounds,
        method="bounded",
        options={"xatol": 1e-6, "maxiter": 500},
    )

    omega_best = result.x
    chi2_best = result.fun

    # Check if the best solution actually has a valid Rt
    Rt_best, n_roots = find_constrained_rt(omega_best, radius, v_bary, a0_half)

    if not np.isfinite(Rt_best) or chi2_best >= 1e29:
        logger.warning("%s [RT constrained]: no valid Rt solution found", galaxy_id)
        nan = float("nan")
        return ModelFitResult(
            galaxy_id=galaxy_id, model_name="constrained_rt", n_params=1,
            chi_squared=nan, reduced_chi_squared=nan, bic=nan,
            residuals_rmse=nan, n_points=n_points, converged=False,
            flag_v_obs_lt_v_bary=flag_v_obs_lt_v_bary,
            method_version=method_version, upsilon_disk=upsilon_disk,
            upsilon_bulge=upsilon_bulge,
            v_bary=v_bary, v_model=np.full_like(radius, nan),
            residuals=np.full_like(radius, nan),
            param1=nan, param1_err=nan, param2=nan, param2_err=nan,
        )

    v_model = v_bary + omega_best * radius / (1.0 + radius / Rt_best)
    residuals = v_obs - v_model
    chi2 = chi2_best
    dof = max(n_points - 1, 1)
    reduced_chi2 = chi2 / dof
    rmse = float(np.sqrt(np.mean(residuals ** 2)))
    bic = compute_bic(n_points, 1, chi2)

    if n_roots > 1:
        logger.info(
            "%s [RT constrained]: %d roots found; using smallest Rt=%.3f kpc",
            galaxy_id, n_roots, Rt_best,
        )

    logger.info(
        "%s [RT constrained]: omega=%.4f  Rt=%.4f  chi2_r=%.2f  RMSE=%.2f  roots=%d",
        galaxy_id, omega_best, Rt_best, reduced_chi2, rmse, n_roots,
    )

    return ModelFitResult(
        galaxy_id=galaxy_id, model_name="constrained_rt", n_params=1,
        chi_squared=chi2, reduced_chi_squared=reduced_chi2, bic=bic,
        residuals_rmse=rmse, n_points=n_points, converged=True,
        flag_v_obs_lt_v_bary=flag_v_obs_lt_v_bary,
        method_version=method_version, upsilon_disk=upsilon_disk,
        upsilon_bulge=upsilon_bulge,
        v_bary=v_bary, v_model=v_model, residuals=residuals,
        param1=omega_best, param1_err=None,
        param2=Rt_best, param2_err=None,
    )
