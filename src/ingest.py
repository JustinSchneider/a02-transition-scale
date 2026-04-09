"""Data loading utilities for Paper 3.

Paper 3 reads from the existing galaxy_dynamics.db (populated in Paper 2)
and SPARC raw data files. No new data ingestion pipeline is needed for the
SPARC baseline; this module provides helpers for loading profiles and metadata.
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
