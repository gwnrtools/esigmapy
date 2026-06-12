import pytest

from esigmapy.legacy import FitMOmegaIMRAttachmentNonSpinning as Fit

_BASE = 1.0 / 6**1.5


@pytest.fixture(autouse=True)
def reset_called_once():
    """Reset the class-level flag before each test so logging fires on first call."""
    Fit.called_once = False
    yield
    Fit.called_once = False


# ---------------------------------------------------------------------------
# fit_quadratic_poly
# ---------------------------------------------------------------------------


class TestFitQuadraticPoly:
    def test_zero_coeffs_returns_base(self):
        assert Fit.fit_quadratic_poly(0.25, [0.0, 0.0]) == pytest.approx(_BASE)

    def test_known_value(self):
        a1, a2, eta = 1.0, 2.0, 0.25
        expected = _BASE * (1.0 + a1 * eta + a2 * eta**2)
        assert Fit.fit_quadratic_poly(eta, [a1, a2]) == pytest.approx(
            expected, rel=1e-14
        )

    def test_wrong_coeff_count_raises(self):
        with pytest.raises(AssertionError):
            Fit.fit_quadratic_poly(0.25, [1.0])

    def test_positive(self):
        assert Fit.fit_quadratic_poly(0.25, [0.0, 0.0]) > 0


# ---------------------------------------------------------------------------
# fit_cubic_poly
# ---------------------------------------------------------------------------


class TestFitCubicPoly:
    def test_zero_coeffs_returns_base(self):
        assert Fit.fit_cubic_poly(0.25, [0.0, 0.0, 0.0]) == pytest.approx(_BASE)

    def test_known_value(self):
        a1, a2, a3, eta = 2.0, -1.0, 0.5, 0.2
        expected = _BASE * (1.0 + a1 * eta + a2 * eta**2 + a3 * eta**3)
        assert Fit.fit_cubic_poly(eta, [a1, a2, a3]) == pytest.approx(
            expected, rel=1e-14
        )

    def test_wrong_coeff_count_raises(self):
        with pytest.raises(AssertionError):
            Fit.fit_cubic_poly(0.25, [1.0, 2.0])


# ---------------------------------------------------------------------------
# fit_ratio_poly_44
# ---------------------------------------------------------------------------


class TestFitRatioPoly44:
    _Z6 = [0.0] * 6

    def test_zero_coeffs_returns_base(self):
        assert Fit.fit_ratio_poly_44(0.25, self._Z6) == pytest.approx(_BASE)

    def test_equal_num_denom_returns_base(self):
        # num == denom → ratio == 1 → result == _BASE regardless of coeffs
        cs = [1.0, 2.0, 3.0, 1.0, 2.0, 3.0]
        assert Fit.fit_ratio_poly_44(0.25, cs) == pytest.approx(_BASE, rel=1e-14)

    def test_known_value(self):
        a1, a2, a3, b1, b2, b3, eta = 1.0, 0.0, 0.0, 0.5, 0.0, 0.0, 0.2
        expected = _BASE * (1.0 + a1 * eta) / (1.0 + b1 * eta)
        assert Fit.fit_ratio_poly_44(eta, [a1, a2, a3, b1, b2, b3]) == pytest.approx(
            expected, rel=1e-14
        )

    def test_wrong_coeff_count_raises(self):
        with pytest.raises(AssertionError):
            Fit.fit_ratio_poly_44(0.25, [1.0] * 5)


# ---------------------------------------------------------------------------
# fit_ratio_sqrt_poly_44
# ---------------------------------------------------------------------------


class TestFitRatioSqrtPoly44:
    _Z6 = [0.0] * 6

    def test_zero_coeffs_returns_base(self):
        assert Fit.fit_ratio_sqrt_poly_44(0.25, self._Z6) == pytest.approx(_BASE)

    def test_uses_sqrt_eta_not_eta(self):
        # With only a1 nonzero: sqrt variant scales by eta^0.5, linear by eta
        coeffs = [1.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        eta = 0.16  # sqrt(0.16) = 0.4, so results differ
        assert Fit.fit_ratio_sqrt_poly_44(eta, coeffs) != pytest.approx(
            Fit.fit_ratio_poly_44(eta, coeffs)
        )

    def test_known_value(self):
        a1, eta = 1.0, 0.16
        s_eta = eta**0.5
        expected = _BASE * (1.0 + a1 * s_eta)
        assert Fit.fit_ratio_sqrt_poly_44(eta, [a1, 0, 0, 0, 0, 0]) == pytest.approx(
            expected, rel=1e-14
        )

    def test_wrong_coeff_count_raises(self):
        with pytest.raises(AssertionError):
            Fit.fit_ratio_sqrt_poly_44(0.25, [1.0] * 4)


# ---------------------------------------------------------------------------
# fit_ratio_sqrt_hyb1_poly_44 — documents copy-paste bug
# ---------------------------------------------------------------------------


class TestFitRatioSqrtHyb1Poly44:
    """fit_ratio_sqrt_hyb1_poly_44 computes s_eta but never uses it;
    the polynomial uses eta instead, making it identical to fit_ratio_poly_44."""

    _Z6 = [0.0] * 6

    def test_zero_coeffs_returns_base(self):
        assert Fit.fit_ratio_sqrt_hyb1_poly_44(0.25, self._Z6) == pytest.approx(_BASE)

    def test_identical_to_fit_ratio_poly_44(self):
        # Bug: s_eta is assigned but eta is used in the body
        coeffs = [1.0, 0.5, 0.2, 0.3, 0.1, 0.05]
        eta = 0.2
        assert Fit.fit_ratio_sqrt_hyb1_poly_44(eta, coeffs) == pytest.approx(
            Fit.fit_ratio_poly_44(eta, coeffs), rel=1e-14
        )

    def test_wrong_coeff_count_raises(self):
        with pytest.raises(AssertionError):
            Fit.fit_ratio_sqrt_hyb1_poly_44(0.25, [1.0] * 7)


# ---------------------------------------------------------------------------
# fit_ratio_poly_43
# ---------------------------------------------------------------------------


class TestFitRatioPoly43:
    _Z5 = [0.0] * 5

    def test_zero_coeffs_returns_base(self):
        assert Fit.fit_ratio_poly_43(0.25, self._Z5) == pytest.approx(_BASE)

    def test_known_value(self):
        a1, a2, a3, b1, b2, eta = 2.0, 0.0, 0.0, 1.0, 0.0, 0.1
        expected = _BASE * (1.0 + a1 * eta) / (1.0 + b1 * eta)
        assert Fit.fit_ratio_poly_43(eta, [a1, a2, a3, b1, b2]) == pytest.approx(
            expected, rel=1e-14
        )

    def test_wrong_coeff_count_raises(self):
        with pytest.raises(AssertionError):
            Fit.fit_ratio_poly_43(0.25, [1.0] * 6)


# ---------------------------------------------------------------------------
# fit_ratio_sqrt_poly_43
# ---------------------------------------------------------------------------


class TestFitRatioSqrtPoly43:
    _Z5 = [0.0] * 5

    def test_zero_coeffs_returns_base(self):
        assert Fit.fit_ratio_sqrt_poly_43(0.25, self._Z5) == pytest.approx(_BASE)

    def test_uses_sqrt_eta_not_eta(self):
        coeffs = [1.0, 0.0, 0.0, 0.0, 0.0]
        eta = 0.16
        assert Fit.fit_ratio_sqrt_poly_43(eta, coeffs) != pytest.approx(
            Fit.fit_ratio_poly_43(eta, coeffs)
        )

    def test_known_value(self):
        a1, eta = 1.0, 0.16
        s_eta = eta**0.5
        expected = _BASE * (1.0 + a1 * s_eta)
        assert Fit.fit_ratio_sqrt_poly_43(eta, [a1, 0, 0, 0, 0]) == pytest.approx(
            expected, rel=1e-14
        )

    def test_wrong_coeff_count_raises(self):
        with pytest.raises(AssertionError):
            Fit.fit_ratio_sqrt_poly_43(0.25, [1.0] * 6)


# ---------------------------------------------------------------------------
# fit_ratio_sqrt_hyb1_poly_43
# ---------------------------------------------------------------------------


class TestFitRatioSqrtHyb1Poly43:
    _Z5 = [0.0] * 5

    def test_zero_coeffs_returns_base(self):
        assert Fit.fit_ratio_sqrt_hyb1_poly_43(0.25, self._Z5) == pytest.approx(_BASE)

    def test_known_value(self):
        # numerator uses eta * s_eta = eta^(3/2)
        a1, a2, a3, b1, b2, eta = 2.0, 0.0, 0.0, 0.0, 0.0, 0.25
        s_eta = eta**0.5
        expected = _BASE * (1.0 + a1 * eta * s_eta)
        assert Fit.fit_ratio_sqrt_hyb1_poly_43(
            eta, [a1, a2, a3, b1, b2]
        ) == pytest.approx(expected, rel=1e-14)

    def test_distinct_from_fit_ratio_poly_43(self):
        # numerator exponents differ (3/2, 5/2, 7/2 vs 1, 2, 3) → different results
        coeffs = [1.0, 0.0, 0.0, 0.0, 0.0]
        eta = 0.16
        assert Fit.fit_ratio_sqrt_hyb1_poly_43(eta, coeffs) != pytest.approx(
            Fit.fit_ratio_poly_43(eta, coeffs)
        )

    def test_wrong_coeff_count_raises(self):
        with pytest.raises(AssertionError):
            Fit.fit_ratio_sqrt_hyb1_poly_43(0.25, [1.0] * 4)


# ---------------------------------------------------------------------------
# fit_ratio_poly_34
# ---------------------------------------------------------------------------


class TestFitRatioPoly34:
    _Z5 = [0.0] * 5

    def test_zero_coeffs_returns_base(self):
        assert Fit.fit_ratio_poly_34(0.25, self._Z5) == pytest.approx(_BASE)

    def test_known_value(self):
        a1, a2, b1, b2, b3, eta = 1.0, 0.5, 2.0, 0.0, 0.0, 0.2
        expected = _BASE * (1.0 + a1 * eta + a2 * eta**2) / (1.0 + b1 * eta)
        assert Fit.fit_ratio_poly_34(eta, [a1, a2, b1, b2, b3]) == pytest.approx(
            expected, rel=1e-14
        )

    def test_wrong_coeff_count_raises(self):
        with pytest.raises(AssertionError):
            Fit.fit_ratio_poly_34(0.25, [1.0] * 6)


# ---------------------------------------------------------------------------
# Logging patch smoke test
# ---------------------------------------------------------------------------


class TestLoggingPatch:
    def test_first_call_does_not_raise(self):
        # called_once is False (reset by fixture) → logging.info fires
        # conftest.py patches logging into legacy's namespace → no NameError
        Fit.fit_quadratic_poly(0.25, [0.0, 0.0])
