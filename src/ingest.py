"""Data loading utilities for Paper 3.

Paper 3 reads from the existing galaxy_dynamics.db (populated in Paper 2)
and SPARC raw data files. Also provides loaders for THINGS (de Blok et al. 2008)
and LITTLE THINGS (Oh et al. 2015) external datasets.
"""

from pathlib import Path

import numpy as np
import pandas as pd

from src.database import get_engine, get_session, query_profiles_as_dataframe
from src.utils import get_project_root, setup_logger

logger = setup_logger(__name__)


# ---------------------------------------------------------------------------
# Database loaders
# ---------------------------------------------------------------------------

def load_galaxy_profiles(galaxy_id: str, db_path: str = None) -> pd.DataFrame:
    """Load radial profiles for a single galaxy from the database."""
    engine = get_engine(db_path)
    session = get_session(engine)
    try:
        return query_profiles_as_dataframe(session, galaxy_id)
    finally:
        session.close()


def load_all_galaxy_ids(db_path: str = None) -> list[str]:
    """Return a sorted list of all galaxy IDs in the database."""
    engine = get_engine(db_path)
    with engine.connect() as conn:
        result = conn.execute(
            __import__("sqlalchemy").text("SELECT galaxy_id FROM galaxies ORDER BY galaxy_id")
        )
        return [row[0] for row in result]


# ---------------------------------------------------------------------------
# SPARC metadata loaders (for M_bar computation)
# ---------------------------------------------------------------------------

def load_sparc_metadata(sparc_mrt_path: str = None) -> dict:
    """Parse SPARC MRT table for L[3.6] and MHI per galaxy.

    Returns dict mapping galaxy_id -> {'L36': float, 'MHI': float}
    where L36 is in 10^9 Lsun and MHI is in 10^9 Msun.
    """
    if sparc_mrt_path is None:
        sparc_mrt_path = str(get_project_root() / "data" / "raw" / "SPARC_Lelli2016c.mrt")

    sparc = {}
    header_done = False
    sep_count = 0

    with open(sparc_mrt_path) as f:
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
                sparc[parts[0]] = {"L36": float(parts[7]), "MHI": float(parts[13])}
            except (ValueError, IndexError):
                continue

    logger.info("Loaded %d galaxies from SPARC metadata", len(sparc))
    return sparc


def load_bulge_luminosities(bulge_path: str = None) -> dict:
    """Parse SPARC bulge luminosity table.

    Returns dict mapping galaxy_id -> L_bulge in 10^9 Lsun.
    """
    if bulge_path is None:
        bulge_path = str(get_project_root() / "data" / "raw" / "Bulges.mrt")

    bulges = {}
    with open(bulge_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) >= 2:
                try:
                    bulges[parts[0]] = float(parts[1])
                except ValueError:
                    continue

    logger.info("Loaded %d bulge luminosities", len(bulges))
    return bulges


def compute_mbar(
    l36_1e9: float,
    mhi_1e9: float,
    l_bulge_1e9: float = 0.0,
    upsilon_disk: float = 0.5,
    upsilon_bulge: float = 0.7,
    helium_factor: float = 1.33,
) -> float:
    """Compute total baryonic mass in Msun.

    M_bar = Y_disk * (L_total - L_bulge) + Y_bulge * L_bulge + 1.33 * MHI
    """
    l_disk = (l36_1e9 - l_bulge_1e9) * 1e9  # Lsun
    l_bul = l_bulge_1e9 * 1e9
    m_gas = mhi_1e9 * 1e9 * helium_factor
    return upsilon_disk * l_disk + upsilon_bulge * l_bul + m_gas


# ---------------------------------------------------------------------------
# THINGS data loaders
# ---------------------------------------------------------------------------

# De Blok et al. (2008) distances in Mpc for all 19 THINGS rotation-curve
# galaxies. Used for arcsec -> kpc conversion. For overlap galaxies, SPARC
# database distances are preferred for consistency with baryonic profiles.
THINGS_DISTANCES_MPC = {
    "DDO154": 4.04, "IC2574": 4.02, "NGC925": 9.16, "NGC2366": 3.44,
    "NGC2403": 3.22, "NGC2841": 14.1, "NGC2903": 8.9, "NGC2976": 3.56,
    "NGC3031": 3.63, "NGC3198": 13.8, "NGC3521": 10.7, "NGC3621": 6.64,
    "NGC3627": 10.05, "NGC4736": 4.66, "NGC4826": 7.50, "NGC5055": 10.1,
    "NGC6946": 5.9, "NGC7331": 14.7, "NGC7793": 3.91,
}


def get_things_distances() -> dict:
    """Return de Blok et al. (2008) distances for all 19 THINGS galaxies."""
    return dict(THINGS_DISTANCES_MPC)


def load_things_overlap_csv(csv_path: str = None) -> tuple[list[dict], list[dict]]:
    """Parse the THINGS/SPARC cross-reference CSV.

    Returns (overlap, nonoverlap) where each entry is a dict with keys:
        things_name  — name with spaces (e.g. "NGC 2403")
        sparc_name   — name without spaces (e.g. "NGC2403"), or None
        file_name    — name without spaces for file lookups
    """
    if csv_path is None:
        csv_path = str(get_project_root() / "data" / "things_sparc_overlap.csv")

    df = pd.read_csv(csv_path)
    overlap, nonoverlap = [], []
    for _, row in df.iterrows():
        things_name = row["THINGS_Name"]
        sparc_name = row["SPARC_Name"] if row["In_SPARC"] else None
        file_name = things_name.replace(" ", "")
        entry = {"things_name": things_name, "sparc_name": sparc_name,
                 "file_name": file_name}
        if row["In_SPARC"]:
            overlap.append(entry)
        else:
            nonoverlap.append(entry)

    logger.info("THINGS overlap: %d, non-overlap: %d", len(overlap), len(nonoverlap))
    return overlap, nonoverlap


def load_things_rotation_curve(curve_path: str, distance_mpc: float) -> pd.DataFrame:
    """Parse a de Blok et al. (2008) .curve.02 rotation curve file.

    Uses columns 1 (radius arcsec), 2 (V_rot km/s), 7 (combined error km/s),
    9 (inclination deg) per the THINGS README.

    Parameters
    ----------
    curve_path : path to .curve.02 file
    distance_mpc : galaxy distance in Mpc (for arcsec -> kpc conversion)

    Returns
    -------
    DataFrame with columns: radius_kpc, v_obs, v_err, inclination
    """
    rows = []
    with open(curve_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 10:
                continue
            r_arcsec = float(parts[0])
            v_rot = float(parts[1])
            v_err = float(parts[6])
            incl = float(parts[8])
            r_kpc = r_arcsec * (distance_mpc * 1000.0) / 206265.0
            rows.append({
                "radius_kpc": r_kpc,
                "v_obs": v_rot,
                "v_err": v_err,
                "inclination": incl,
            })

    df = pd.DataFrame(rows)
    logger.info("Loaded THINGS rotation curve: %s (%d points, R_max=%.2f kpc)",
                Path(curve_path).stem, len(df),
                df["radius_kpc"].max() if len(df) > 0 else 0.0)
    return df


def load_things_mass_model(model_path: str) -> pd.DataFrame:
    """Parse a de Blok et al. (2008) .ISO.fix.REV.dat mass model file.

    Reads ROTMAS output: radius_kpc, v_gas, v_disk, v_bulge, v_obs, v_err.
    Skips the !# header block.

    Returns
    -------
    DataFrame with columns: radius_kpc, v_gas, v_disk, v_bulge, v_obs, v_err
    """
    rows = []
    with open(model_path, "rb") as f:
        for raw_line in f:
            line = raw_line.decode("utf-8", errors="replace").strip()
            if not line or line.startswith("!"):
                continue
            parts = line.split()
            if len(parts) < 6:
                continue
            try:
                rows.append({
                    "radius_kpc": float(parts[0]),
                    "v_gas": float(parts[1]),
                    "v_disk": float(parts[2]),
                    "v_bulge": float(parts[3]),
                    "v_obs": float(parts[4]),
                    "v_err": float(parts[5]),
                })
            except ValueError:
                continue

    df = pd.DataFrame(rows)
    logger.info("Loaded THINGS mass model: %s (%d points)",
                Path(model_path).stem, len(df))
    return df


# ---------------------------------------------------------------------------
# LITTLE THINGS data loaders (Oh et al. 2015)
# ---------------------------------------------------------------------------

def load_little_things_galaxies(csv_path: str = None) -> pd.DataFrame:
    """Load LITTLE THINGS galaxy metadata from Oh et al. (2015) galaxies.csv.

    Replaces sentinel values (1e+20) for missing stellar masses with NaN.

    Returns DataFrame with all columns from the VizieR galaxies table.
    Key columns: Name, Dist (Mpc), i (deg), Mgas (1e7 Msun),
    MstarK (1e7 Msun), MstarSED (1e7 Msun), Rmax (kpc).
    """
    if csv_path is None:
        csv_path = str(
            get_project_root() / "data" / "raw" / "LITTLE_THINGS" / "galaxies.csv"
        )
    df = pd.read_csv(csv_path)
    # Replace sentinel 1e+20 with NaN for missing stellar masses
    for col in ["MstarK", "MstarSED"]:
        if col in df.columns:
            df.loc[df[col] > 1e10, col] = np.nan
    logger.info("Loaded %d LITTLE THINGS galaxies", len(df))
    return df


def load_little_things_rotcurves(
    dmbar_path: str = None,
    dm_path: str = None,
) -> tuple[dict[str, pd.DataFrame], dict[str, pd.DataFrame]]:
    """Load and de-normalize LITTLE THINGS rotation curves.

    Reads rotdmbar.csv (total V_obs) and rotdm.csv (DM-only, baryonic
    subtracted). Filters to Type='Data' rows and converts from normalized
    to physical units using per-galaxy R0.3 (kpc) and V0.3 (km/s).

    Returns:
        (dmbar_dict, dm_dict) where each maps galaxy Name -> DataFrame
        with columns: radius_kpc, v_kms, v_err_kms
    """
    lt_dir = get_project_root() / "data" / "raw" / "LITTLE_THINGS"
    if dmbar_path is None:
        dmbar_path = str(lt_dir / "rotdmbar.csv")
    if dm_path is None:
        dm_path = str(lt_dir / "rotdm.csv")

    def _load_and_denorm(path: str) -> dict[str, pd.DataFrame]:
        raw = pd.read_csv(path)
        data = raw[raw["Type"] == "Data"].copy()
        result = {}
        for name, grp in data.groupby("Name"):
            r03 = grp["R0.3"].iloc[0]
            v03 = grp["V0.3"].iloc[0]
            result[name] = pd.DataFrame({
                "radius_kpc": grp["R"].values * r03,
                "v_kms": grp["V"].values * v03,
                "v_err_kms": grp["e_V"].values * v03,
            })
        return result

    dmbar = _load_and_denorm(dmbar_path)
    dm = _load_and_denorm(dm_path)
    logger.info("Loaded LITTLE THINGS rotation curves: %d dmbar, %d dm",
                len(dmbar), len(dm))
    return dmbar, dm


def derive_little_things_vbary(
    dmbar_df: pd.DataFrame,
    dm_df: pd.DataFrame,
) -> pd.DataFrame:
    """Derive baryonic velocity from total and DM-only rotation curves.

    V_bary = sqrt(max(0, V_total^2 - V_dm^2)) at each radius.

    The DM curve is interpolated onto the total curve's radial grid.

    Returns DataFrame with columns:
        radius_kpc, v_obs, v_err, v_baryon_total, n_clamped
    where n_clamped is the number of points where V_dm > V_total (clamped to 0).
    """
    r_obs = dmbar_df["radius_kpc"].values
    v_obs = dmbar_df["v_kms"].values
    v_err = dmbar_df["v_err_kms"].values

    r_dm = dm_df["radius_kpc"].values
    v_dm_raw = dm_df["v_kms"].values

    # Interpolate DM velocity onto the observed radial grid
    v_dm_interp = np.interp(r_obs, r_dm, v_dm_raw)

    # V_bary^2 = V_total^2 - V_dm^2
    v_bary_sq = v_obs ** 2 - v_dm_interp ** 2
    n_clamped = int(np.sum(v_bary_sq < 0))
    v_bary = np.sqrt(np.maximum(0.0, v_bary_sq))

    return pd.DataFrame({
        "radius_kpc": r_obs,
        "v_obs": v_obs,
        "v_err": v_err,
        "v_baryon_total": v_bary,
    }), n_clamped


def compute_mbar_little_things(
    mgas_1e7: float,
    mstar_1e7: float,
    helium_factor: float = 1.33,
) -> float:
    """Compute baryonic mass from LITTLE THINGS table values.

    Args:
        mgas_1e7: HI gas mass in units of 10^7 Msun.
        mstar_1e7: Stellar mass in units of 10^7 Msun (NaN if missing).
        helium_factor: Helium correction (default 1.33).

    Returns:
        Total baryonic mass in Msun.
    """
    m_gas = mgas_1e7 * 1e7
    m_star = mstar_1e7 * 1e7 if np.isfinite(mstar_1e7) else 0.0
    return helium_factor * m_gas + m_star
