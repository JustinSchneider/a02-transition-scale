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
