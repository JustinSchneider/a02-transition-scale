"""Tests for the BTFR covariance test pipeline."""

import numpy as np
import pytest

from src.btfr import (
    compute_btfr_residuals,
    run_btfr_covariance_test,
    compute_mbar_for_sample,
    BTFR_ALPHA,
    BTFR_BETA,
)


# ---------------------------------------------------------------------------
# Test compute_btfr_residuals
# ---------------------------------------------------------------------------

class TestComputeBtfrResiduals:
    def test_anchor_point(self):
        """At the Paper 2 anchor (g=6.51e-11, M=5.80e9), residual ~ 0."""
        g = np.array([6.51e-11])
        m = np.array([5.80e9])
        residuals = compute_btfr_residuals(g, m)
        # beta was derived so that g_trend at M=5.80e9 equals a0/2 = 6.00e-11
        # so residual = log10(6.51e-11) - log10(6.00e-11) = log10(6.51/6.00) ~ 0.035
        assert abs(residuals[0]) < 0.1

    def test_known_value(self):
        """Verify residual matches manual computation."""
        g = np.array([1e-10])
        m = np.array([1e10])
        log_g = np.log10(1e-10)  # -10
        log_g_trend = BTFR_ALPHA * np.log10(1e10) + BTFR_BETA  # 0.238*10 - 12.55 = -10.17
        expected = log_g - log_g_trend
        residuals = compute_btfr_residuals(g, m)
        np.testing.assert_allclose(residuals[0], expected, rtol=1e-10)

    def test_filters_invalid(self):
        """Zero and negative values are filtered out."""
        g = np.array([1e-10, 0.0, -1e-10, 1e-11])
        m = np.array([1e10, 1e9, 1e9, 1e9])
        residuals = compute_btfr_residuals(g, m)
        assert len(residuals) == 2  # only first and last valid

    def test_all_invalid_raises(self):
        """All-invalid input raises ValueError."""
        g = np.array([0.0, -1.0])
        m = np.array([1e10, 1e10])
        with pytest.raises(ValueError):
            compute_btfr_residuals(g, m)

    def test_multiple_values(self):
        """Multiple valid entries produce correct-length output."""
        rng = np.random.default_rng(42)
        g = rng.uniform(1e-11, 1e-9, 50)
        m = rng.uniform(1e8, 1e11, 50)
        residuals = compute_btfr_residuals(g, m)
        assert len(residuals) == 50


# ---------------------------------------------------------------------------
# Test run_btfr_covariance_test
# ---------------------------------------------------------------------------

class TestRunBtfrCovarianceTest:
    def test_result_structure(self):
        """Returned dict has all required keys."""
        rng = np.random.default_rng(42)
        g = rng.uniform(1e-11, 1e-9, 50)
        m = rng.uniform(1e8, 1e11, 50)
        result = run_btfr_covariance_test(g, m)
        expected_keys = {
            "residuals", "scatter_raw", "scatter_residual",
            "wilcoxon_stat", "wilcoxon_pvalue", "n_galaxies", "median_residual",
        }
        assert set(result.keys()) == expected_keys
        assert result["n_galaxies"] == 50

    def test_scatter_reduction(self):
        """Residual scatter should be less than raw scatter when a trend exists."""
        rng = np.random.default_rng(99)
        log_m = rng.uniform(8, 11, 100)
        # Embed a BTFR-like trend in g
        log_g = BTFR_ALPHA * log_m + BTFR_BETA + rng.normal(0, 0.3, 100)
        g = 10**log_g
        m = 10**log_m
        result = run_btfr_covariance_test(g, m)
        assert result["scatter_residual"] < result["scatter_raw"]

    def test_centered_residuals(self):
        """Residuals drawn from a symmetric distribution should not reject H0."""
        rng = np.random.default_rng(7)
        log_m = rng.uniform(8, 11, 80)
        # Generate g exactly on the trend line + symmetric noise
        log_g = BTFR_ALPHA * log_m + BTFR_BETA + rng.normal(0, 0.2, 80)
        g = 10**log_g
        m = 10**log_m
        result = run_btfr_covariance_test(g, m)
        assert result["wilcoxon_pvalue"] > 0.01  # not strongly significant


# ---------------------------------------------------------------------------
# Test compute_mbar_for_sample
# ---------------------------------------------------------------------------

class TestComputeMbarForSample:
    def test_known_galaxy(self):
        """A known SPARC galaxy returns a positive M_bar."""
        # NGC2403 is in every SPARC release
        result = compute_mbar_for_sample(["NGC2403"])
        assert "NGC2403" in result
        assert result["NGC2403"] > 0
        assert np.isfinite(result["NGC2403"])

    def test_missing_galaxy(self):
        """A non-existent galaxy is omitted, not an error."""
        result = compute_mbar_for_sample(["FAKE_GALAXY_XYZ"])
        assert "FAKE_GALAXY_XYZ" not in result
        assert len(result) == 0

    def test_mixed(self):
        """Mix of real and fake galaxies returns only the real ones."""
        result = compute_mbar_for_sample(["NGC2403", "FAKE_GALAXY_XYZ"])
        assert "NGC2403" in result
        assert "FAKE_GALAXY_XYZ" not in result
        assert len(result) == 1
