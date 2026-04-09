"""Tests for data loading utilities in src/ingest.py."""

import numpy as np
import pytest

from src.ingest import load_sparc_metadata, load_bulge_luminosities, compute_mbar


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
