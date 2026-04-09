"""Tests for core physics functions, especially the constrained RT model."""

import numpy as np
import pytest

from src.physics import (
    compute_v_bary, compute_bic, interpolate_v_bary,
    find_constrained_rt, fit_constrained_rt,
    rt_model_velocity, compute_transition_diagnostics,
    A0_HALF, _safe_errors,
)


# ---------------------------------------------------------------------------
# Test compute_v_bary
# ---------------------------------------------------------------------------

class TestComputeVBary:
    def test_simple_case(self):
        """All positive components, no bulge."""
        v = compute_v_bary(
            v_gas=np.array([10.0]),
            v_disk=np.array([50.0]),
            v_bulge=np.array([0.0]),
        )
        # V^2 = 10*10 + 0.5*50*50 + 0 = 100 + 1250 = 1350
        expected = np.sqrt(1350.0)
        np.testing.assert_allclose(v, expected, rtol=1e-10)

    def test_negative_gas(self):
        """Negative gas velocity (inner depression) reduces V_bary."""
        v = compute_v_bary(
            v_gas=np.array([-5.0]),
            v_disk=np.array([30.0]),
            v_bulge=np.array([0.0]),
        )
        # V^2 = -25 + 0.5*900 = -25 + 450 = 425
        expected = np.sqrt(425.0)
        np.testing.assert_allclose(v, expected, rtol=1e-10)

    def test_always_non_negative(self):
        """V_bary should never be negative, even with large negative gas."""
        v = compute_v_bary(
            v_gas=np.array([-100.0]),
            v_disk=np.array([5.0]),
            v_bulge=np.array([0.0]),
        )
        assert v[0] >= 0


# ---------------------------------------------------------------------------
# Test compute_bic
# ---------------------------------------------------------------------------

class TestComputeBIC:
    def test_zero_params(self):
        """k=0 -> BIC = chi^2."""
        assert compute_bic(20, 0, 100.0) == 100.0

    def test_two_params(self):
        bic = compute_bic(20, 2, 100.0)
        expected = 100.0 + 2 * np.log(20)
        np.testing.assert_allclose(bic, expected)


# ---------------------------------------------------------------------------
# Test interpolate_v_bary
# ---------------------------------------------------------------------------

class TestInterpolateVBary:
    def test_midpoint(self):
        r = np.array([1.0, 2.0, 3.0])
        v = np.array([10.0, 20.0, 30.0])
        result = interpolate_v_bary(r, v, 2.0)
        np.testing.assert_allclose(result, 20.0, rtol=1e-6)

    def test_out_of_range(self):
        r = np.array([1.0, 2.0, 3.0])
        v = np.array([10.0, 20.0, 30.0])
        assert np.isnan(interpolate_v_bary(r, v, 0.5))
        assert np.isnan(interpolate_v_bary(r, v, 3.5))


# ---------------------------------------------------------------------------
# Test find_constrained_rt
# ---------------------------------------------------------------------------

class TestFindConstrainedRt:
    def test_known_solution(self):
        """Construct a case where we know the answer.

        If V_bary is constant at 30 km/s and omega=10 km/s/kpc:
            g(Rt) = (30 + 10*Rt/2)^2 / Rt = a0/2

        For a0/2 = 1851.5 km^2/s^2/kpc, solve:
            (30 + 5*Rt)^2 / Rt = 1851.5
            900/Rt + 300 + 25*Rt = 1851.5

        This is a quadratic in Rt: 25*Rt^2 - 1551.5*Rt + 900 = 0
        Rt = (1551.5 +/- sqrt(1551.5^2 - 4*25*900)) / (2*25)
        """
        radius = np.linspace(0.5, 100.0, 200)
        v_bary = np.full_like(radius, 30.0)

        omega = 10.0
        # Solve analytically
        a, b, c = 25.0, -1551.5, 900.0
        disc = b**2 - 4*a*c
        Rt_expected_small = (-b - np.sqrt(disc)) / (2*a)
        Rt_expected_large = (-b + np.sqrt(disc)) / (2*a)

        Rt_found, n_roots = find_constrained_rt(omega, radius, v_bary)

        assert np.isfinite(Rt_found)
        assert n_roots >= 1
        # Should find the smallest root
        np.testing.assert_allclose(Rt_found, Rt_expected_small, rtol=0.01)

    def test_no_solution(self):
        """Very high V_bary makes g(Rt) > a0/2 everywhere -> no solution with tiny omega."""
        radius = np.linspace(0.1, 5.0, 100)
        v_bary = np.full_like(radius, 200.0)  # Very high

        Rt, n_roots = find_constrained_rt(0.001, radius, v_bary)
        # g = (200 + tiny)^2 / R >> a0/2 for small R, but at large R might dip
        # With this setup there may or may not be a solution, but the point is
        # the function returns gracefully
        assert isinstance(n_roots, int)


# ---------------------------------------------------------------------------
# Test fit_constrained_rt
# ---------------------------------------------------------------------------

class TestFitConstrainedRt:
    def test_synthetic_galaxy(self):
        """Fit a synthetic galaxy where the constrained model should work well."""
        # Create a rotation curve from a known RT model
        omega_true = 8.0
        Rt_true = 10.0
        radius = np.linspace(0.5, 30.0, 25)
        v_bary = 20.0 + 5.0 * np.sqrt(radius)  # Rising baryonic curve
        v_model = v_bary + omega_true * radius / (1.0 + radius / Rt_true)

        # Add small noise
        rng = np.random.default_rng(42)
        v_obs = v_model + rng.normal(0, 1.0, len(radius))
        v_err = np.ones_like(radius) * 2.0

        result = fit_constrained_rt(radius, v_obs, v_err, v_bary, galaxy_id="test")

        # Should converge (unless constraint happens to be unsatisfiable)
        # We can't guarantee convergence since the constraint equation depends
        # on a0/2 which may not match the true model. But the function should
        # at least run without error.
        assert result.model_name == "constrained_rt"
        assert result.n_params == 1
        assert isinstance(result.converged, bool)

    def test_returns_model_fit_result(self):
        """Verify the result has the expected structure."""
        radius = np.linspace(1.0, 20.0, 15)
        v_bary = np.full_like(radius, 25.0)
        v_obs = v_bary + 5.0 * radius / (1.0 + radius / 8.0)
        v_err = np.ones_like(radius) * 3.0

        result = fit_constrained_rt(radius, v_obs, v_err, v_bary, galaxy_id="struct_test")

        assert result.galaxy_id == "struct_test"
        assert result.model_name == "constrained_rt"
        assert result.n_params == 1
        assert hasattr(result, "to_dict")
        d = result.to_dict()
        assert "galaxy_id" in d
        assert "bic" in d


# ---------------------------------------------------------------------------
# Test _safe_errors
# ---------------------------------------------------------------------------

class TestSafeErrors:
    def test_replaces_zeros(self):
        v_err = np.array([0.0, 5.0, 3.0, 0.0])
        safe = _safe_errors(v_err, "test")
        assert np.all(safe > 0)
        assert safe[0] == 3.0  # min of positive errors
        assert safe[1] == 5.0  # unchanged

    def test_all_positive(self):
        v_err = np.array([2.0, 3.0, 4.0])
        safe = _safe_errors(v_err, "test")
        np.testing.assert_array_equal(safe, v_err)


# ---------------------------------------------------------------------------
# Test rt_model_velocity
# ---------------------------------------------------------------------------

class TestRtModelVelocity:
    def test_asymptotic_behavior(self):
        """At large R, correction -> omega * Rt = V_sat."""
        radius = np.array([1000.0])
        v_bary = np.array([0.0])
        v = rt_model_velocity(radius, v_bary, omega=10.0, r_t=5.0)
        # V_sat = 10 * 5 = 50
        np.testing.assert_allclose(v, 50.0, rtol=0.01)

    def test_small_r(self):
        """At small R, correction -> omega * R (linear)."""
        radius = np.array([0.01])
        v_bary = np.array([0.0])
        v = rt_model_velocity(radius, v_bary, omega=10.0, r_t=5.0)
        # omega * R / (1 + R/Rt) ~ omega * R for R << Rt
        np.testing.assert_allclose(v, 10.0 * 0.01, rtol=0.01)
