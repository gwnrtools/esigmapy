import numpy as np
import pytest
from unittest.mock import MagicMock, patch

from esigmapy.mr_generator import (
    check_available_mr_approximants,
    _generate_lalsim_modes,
    get_mr_modes,
    LALSIM_APPROXIMANTS,
    PYSEOBNR_APPROXIMANTS,
    SUPPORTED_MR_APPROXIMANTS,
)


class TestCheckAvailableMrApproximants:
    def test_lalsim_approximants_valid(self):
        for approx in LALSIM_APPROXIMANTS:
            check_available_mr_approximants(approx)  # must not raise

    def test_pyseobnr_approximants_valid(self):
        for approx in PYSEOBNR_APPROXIMANTS:
            check_available_mr_approximants(approx)  # must not raise

    def test_unsupported_approximant_raises(self):
        with pytest.raises(IOError, match="cannot generate"):
            check_available_mr_approximants("IMRPhenomD")

    def test_empty_string_raises(self):
        with pytest.raises(IOError):
            check_available_mr_approximants("")

    def test_case_sensitive(self):
        with pytest.raises(IOError):
            check_available_mr_approximants("nrsur7dq4")


class TestApproximantConstants:
    def test_supported_is_union_of_lalsim_and_pyseobnr(self):
        assert set(SUPPORTED_MR_APPROXIMANTS) == set(LALSIM_APPROXIMANTS) | set(
            PYSEOBNR_APPROXIMANTS
        )

    def test_lalsim_and_pyseobnr_disjoint(self):
        assert not set(LALSIM_APPROXIMANTS) & set(PYSEOBNR_APPROXIMANTS)

    def test_lalsim_contains_nrsur7dq4(self):
        assert "NRSur7dq4" in LALSIM_APPROXIMANTS

    def test_lalsim_contains_seobnrv4phm(self):
        assert "SEOBNRv4PHM" in LALSIM_APPROXIMANTS

    def test_pyseobnr_contains_seobnrv5hm(self):
        assert "SEOBNRv5HM" in PYSEOBNR_APPROXIMANTS

    def test_pyseobnr_contains_seobnrv5phm(self):
        assert "SEOBNRv5PHM" in PYSEOBNR_APPROXIMANTS


# ---------------------------------------------------------------------------
# _generate_lalsim_modes — Tier 2 (LAL mocked)
# ---------------------------------------------------------------------------


def _mode_node(l, m, data, nxt=None):
    """Build one node of the SphHarmTimeSeries linked list that LAL returns."""
    node = MagicMock()
    node.l = l
    node.m = m
    node.mode.data.data = np.asarray(data, dtype=complex)
    node.next = nxt
    return node


class TestGenerateLalsimModesWithMock:
    """Tests for _generate_lalsim_modes: mode filtering and data extraction."""

    _N = 64

    def test_bad_approximant_raises_ioerror(self):
        # Real lalsimulation does not have this attribute → AttributeError → IOError
        with pytest.raises(IOError, match="not available"):
            _generate_lalsim_modes(
                10,
                10,
                20,
                1.0 / 4096,
                0,
                0,
                0,
                1,
                20,
                modes_to_use=[(2, 2)],
                approximant="NotARealApproximant99999",
            )

    def test_returns_only_requested_modes(self):
        # Linked list has (2,2) and (3,3); only (2,2) is requested
        node_33 = _mode_node(3, 3, np.ones(self._N, dtype=complex))
        node_22 = _mode_node(2, 2, np.ones(self._N, dtype=complex), nxt=node_33)
        with patch("esigmapy.mr_generator.ls") as mock_ls:
            mock_ls.NRSur7dq4 = 900
            mock_ls.SimInspiralChooseTDModes.return_value = node_22
            result = _generate_lalsim_modes(
                10,
                10,
                20,
                1.0 / 4096,
                0,
                0,
                0,
                1,
                20,
                modes_to_use=[(2, 2)],
                approximant="NRSur7dq4",
            )
        assert set(result.keys()) == {(2, 2)}

    def test_returns_all_requested_modes(self):
        # Linked list has (2,2) and (3,3); both requested
        node_33 = _mode_node(3, 3, np.ones(self._N, dtype=complex))
        node_22 = _mode_node(2, 2, np.ones(self._N, dtype=complex), nxt=node_33)
        with patch("esigmapy.mr_generator.ls") as mock_ls:
            mock_ls.NRSur7dq4 = 900
            mock_ls.SimInspiralChooseTDModes.return_value = node_22
            result = _generate_lalsim_modes(
                10,
                10,
                20,
                1.0 / 4096,
                0,
                0,
                0,
                1,
                20,
                modes_to_use=[(2, 2), (3, 3)],
                approximant="NRSur7dq4",
            )
        assert set(result.keys()) == {(2, 2), (3, 3)}

    def test_mode_data_is_numpy_array(self):
        node = _mode_node(2, 2, np.ones(self._N, dtype=complex))
        node.next = None
        with patch("esigmapy.mr_generator.ls") as mock_ls:
            mock_ls.NRSur7dq4 = 900
            mock_ls.SimInspiralChooseTDModes.return_value = node
            result = _generate_lalsim_modes(
                10,
                10,
                20,
                1.0 / 4096,
                0,
                0,
                0,
                1,
                20,
                modes_to_use=[(2, 2)],
                approximant="NRSur7dq4",
            )
        assert isinstance(result[(2, 2)], np.ndarray)

    def test_unrequested_mode_in_list_is_ignored(self):
        # Linked list has (4,4) only; request is (2,2) — result should be empty
        node = _mode_node(4, 4, np.ones(self._N, dtype=complex))
        node.next = None
        with patch("esigmapy.mr_generator.ls") as mock_ls:
            mock_ls.NRSur7dq4 = 900
            mock_ls.SimInspiralChooseTDModes.return_value = node
            result = _generate_lalsim_modes(
                10,
                10,
                20,
                1.0 / 4096,
                0,
                0,
                0,
                1,
                20,
                modes_to_use=[(2, 2)],
                approximant="NRSur7dq4",
            )
        assert (2, 2) not in result


# ---------------------------------------------------------------------------
# get_mr_modes — Tier 2 (internal helpers mocked)
# ---------------------------------------------------------------------------


class TestGetMrModesWithMock:
    """Tests for the Python routing and defaulting logic in get_mr_modes."""

    _FAKE = {(2, 2): np.ones(64, dtype=complex)}

    def test_f_ref_defaults_to_f_lower(self):
        with patch("esigmapy.mr_generator._generate_lalsim_modes") as mock_gen:
            mock_gen.return_value = self._FAKE
            get_mr_modes(
                10, 10, 20.0, 1.0 / 4096, modes_to_use=[(2, 2)], approximant="NRSur7dq4"
            )
        assert mock_gen.call_args.kwargs["f_ref"] == 20.0

    def test_coa_phase_defaults_to_zero(self):
        with patch("esigmapy.mr_generator._generate_lalsim_modes") as mock_gen:
            mock_gen.return_value = self._FAKE
            get_mr_modes(
                10,
                10,
                20.0,
                1.0 / 4096,
                coa_phase=None,
                modes_to_use=[(2, 2)],
                approximant="NRSur7dq4",
            )
        assert mock_gen.call_args.kwargs["coa_phase"] == 0.0

    def test_lalsim_approximant_routes_to_lalsim(self):
        with patch(
            "esigmapy.mr_generator._generate_lalsim_modes"
        ) as mock_lalsim, patch(
            "esigmapy.mr_generator._generate_pyseobnr_modes"
        ) as mock_pyseobnr:
            mock_lalsim.return_value = self._FAKE
            get_mr_modes(
                10, 10, 20.0, 1.0 / 4096, modes_to_use=[(2, 2)], approximant="NRSur7dq4"
            )
        mock_lalsim.assert_called_once()
        mock_pyseobnr.assert_not_called()

    def test_pyseobnr_approximant_routes_to_pyseobnr(self):
        with patch(
            "esigmapy.mr_generator._generate_pyseobnr_modes"
        ) as mock_pyseobnr, patch(
            "esigmapy.mr_generator._generate_lalsim_modes"
        ) as mock_lalsim:
            mock_pyseobnr.return_value = self._FAKE
            get_mr_modes(
                10,
                10,
                20.0,
                1.0 / 4096,
                modes_to_use=[(2, 2)],
                approximant="SEOBNRv5HM",
            )
        mock_pyseobnr.assert_called_once()
        mock_lalsim.assert_not_called()
