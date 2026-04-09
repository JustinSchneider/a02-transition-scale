"""Tests for data loading utilities in src/ingest.py."""

import numpy as np
import pytest

from src.ingest import (
    load_sparc_metadata, load_bulge_luminosities, compute_mbar,
    load_things_overlap_csv, load_things_rotation_curve,
    load_things_mass_model, get_things_distances,
)
from src.utils import get_project_root


# ---------------------------------------------------------------------------
# Test load_sparc_metadata
# ---------------------------------------------------------------------------

class TestLoadSparcMetadata:
    def test_returns_175_galaxies(self):
        """SPARC catalog has exactly 175 galaxies."""
        sparc = load_sparc_metadata()
        assert len(sparc) == 175

    def test_entry_structure(self):
        """Each entry has L36 and MHI as positive floats."""
        sparc = load_sparc_metadata()
        for gid, meta in sparc.items():
            assert "L36" in meta, f"{gid} missing L36"
            assert "MHI" in meta, f"{gid} missing MHI"
            assert meta["L36"] > 0, f"{gid} has non-positive L36"
            assert meta["MHI"] > 0, f"{gid} has non-positive MHI"

    def test_known_galaxy_present(self):
        """NGC2403 should be in the SPARC catalog."""
        sparc = load_sparc_metadata()
        assert "NGC2403" in sparc


# ---------------------------------------------------------------------------
# Test load_bulge_luminosities
# ---------------------------------------------------------------------------

class TestLoadBulgeLuminosities:
    def test_returns_175_entries(self):
        """Bulge file has one entry per SPARC galaxy."""
        bulges = load_bulge_luminosities()
        assert len(bulges) == 175

    def test_values_non_negative(self):
        """All bulge luminosities are non-negative."""
        bulges = load_bulge_luminosities()
        for gid, l_bulge in bulges.items():
            assert l_bulge >= 0, f"{gid} has negative L_bulge"

    def test_some_nonzero(self):
        """At least some galaxies have nonzero bulge luminosities."""
        bulges = load_bulge_luminosities()
        nonzero = sum(1 for v in bulges.values() if v > 0)
        assert nonzero > 0


# ---------------------------------------------------------------------------
# Test compute_mbar
# ---------------------------------------------------------------------------

class TestComputeMbar:
    def test_disk_only(self):
        """Disk-only galaxy: M_bar = Y_d * L + 1.33 * MHI."""
        # L36 = 1.0 (10^9 Lsun), MHI = 0.5 (10^9 Msun), L_bulge = 0
        m = compute_mbar(1.0, 0.5, 0.0)
        expected = 0.5 * 1e9 + 1.33 * 0.5e9  # 5e8 + 6.65e8 = 1.165e9
        np.testing.assert_allclose(m, expected, rtol=1e-10)

    def test_with_bulge(self):
        """Bulge contribution adds (Y_b - Y_d) * L_bulge relative to disk-only."""
        m_no_bulge = compute_mbar(2.0, 1.0, 0.0)
        m_with_bulge = compute_mbar(2.0, 1.0, 0.5)
        # Difference: (0.7 - 0.5) * 0.5e9 = 0.1e9
        delta = m_with_bulge - m_no_bulge
        np.testing.assert_allclose(delta, 0.2 * 0.5e9, rtol=1e-10)

    def test_zero_mass(self):
        """Zero luminosity and gas gives zero M_bar."""
        m = compute_mbar(0.0, 0.0, 0.0)
        assert m == 0.0


# ---------------------------------------------------------------------------
# Test THINGS overlap CSV
# ---------------------------------------------------------------------------

class TestLoadThingsOverlapCSV:
    def test_overlap_count(self):
        """13 THINGS galaxies overlap with SPARC."""
        overlap, nonoverlap = load_things_overlap_csv()
        assert len(overlap) == 13

    def test_nonoverlap_count(self):
        """6 THINGS galaxies are not in SPARC."""
        overlap, nonoverlap = load_things_overlap_csv()
        assert len(nonoverlap) == 6

    def test_ngc2403_in_overlap(self):
        """NGC2403 should be in the overlap list with matching SPARC name."""
        overlap, _ = load_things_overlap_csv()
        sparc_names = [e["sparc_name"] for e in overlap]
        assert "NGC2403" in sparc_names

    def test_entry_has_file_name(self):
        """Every entry has a file_name field (no spaces)."""
        overlap, nonoverlap = load_things_overlap_csv()
        for entry in overlap + nonoverlap:
            assert "file_name" in entry
            assert " " not in entry["file_name"]


# ---------------------------------------------------------------------------
# Test THINGS rotation curve parser
# ---------------------------------------------------------------------------

class TestLoadThingsRotationCurve:
    @pytest.fixture
    def ngc2403_rc(self):
        root = get_project_root()
        path = root / "data" / "raw" / "THINGS" / "Curves" / "NGC2403.curve.02"
        return load_things_rotation_curve(str(path), distance_mpc=3.16)

    def test_row_count(self, ngc2403_rc):
        """NGC2403 has 288 radial points."""
        assert len(ngc2403_rc) == 288

    def test_radius_in_kpc(self, ngc2403_rc):
        """First radius should be ~0.061 kpc (4 arcsec at 3.16 Mpc)."""
        expected = 4.0 * (3.16 * 1000.0) / 206265.0
        np.testing.assert_allclose(ngc2403_rc["radius_kpc"].iloc[0], expected, rtol=1e-6)

    def test_velocities_positive(self, ngc2403_rc):
        """NGC2403 V_obs values should all be positive."""
        assert (ngc2403_rc["v_obs"] > 0).all()

    def test_errors_positive(self, ngc2403_rc):
        """Combined errors should all be positive."""
        assert (ngc2403_rc["v_err"] > 0).all()

    def test_columns(self, ngc2403_rc):
        """DataFrame has expected columns."""
        assert set(ngc2403_rc.columns) == {"radius_kpc", "v_obs", "v_err", "inclination"}

    def test_rmax_reasonable(self, ngc2403_rc):
        """R_max should be in the range 10-25 kpc for NGC2403."""
        rmax = ngc2403_rc["radius_kpc"].max()
        assert 10 < rmax < 25


# ---------------------------------------------------------------------------
# Test THINGS mass model parser
# ---------------------------------------------------------------------------

class TestLoadThingsMassModel:
    @pytest.fixture
    def ngc2403_mm(self):
        root = get_project_root()
        path = root / "data" / "raw" / "THINGS" / "MassModels" / "NGC2403.ISO.fix.REV.dat"
        return load_things_mass_model(str(path))

    def test_row_count(self, ngc2403_mm):
        """NGC2403 mass model has 288 points (matches rotation curve)."""
        assert len(ngc2403_mm) == 288

    def test_first_radius(self, ngc2403_mm):
        """First radius matches file value exactly."""
        np.testing.assert_allclose(ngc2403_mm["radius_kpc"].iloc[0], 0.06252, rtol=1e-4)

    def test_columns(self, ngc2403_mm):
        """DataFrame has expected columns."""
        expected = {"radius_kpc", "v_gas", "v_disk", "v_bulge", "v_obs", "v_err"}
        assert set(ngc2403_mm.columns) == expected

    def test_velocities_reasonable(self, ngc2403_mm):
        """V_obs values are in physically reasonable range (0-300 km/s)."""
        assert ngc2403_mm["v_obs"].max() < 300
        assert ngc2403_mm["v_obs"].min() > -10  # allow small negatives at center


# ---------------------------------------------------------------------------
# Test THINGS distance lookup
# ---------------------------------------------------------------------------

class TestGetThingsDistances:
    def test_count(self):
        """19 THINGS galaxies in the distance table."""
        dists = get_things_distances()
        assert len(dists) == 19

    def test_all_positive(self):
        """All distances are positive."""
        dists = get_things_distances()
        assert all(d > 0 for d in dists.values())

    def test_ngc2403_distance(self):
        """NGC2403 distance matches de Blok (2008)."""
        dists = get_things_distances()
        np.testing.assert_allclose(dists["NGC2403"], 3.22, rtol=0.01)
