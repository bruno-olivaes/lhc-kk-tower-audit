"""Unit tests for the LHC Kaluza-Klein Tower Bayesian & Multiplicity Audit."""

import math
from pathlib import Path

from src.kk_audit_engine import (
    ATLAS_PRD_108,
    CMS_EXO_21_010_TABLE2,
    CMS_PAS_EXO_25_004,
    audit_bic_and_five_gates,
    audit_cross_experiment_replication,
    audit_look_elsewhere,
    audit_spectral_hierarchy,
    run_full_audit,
)


def test_spectral_hierarchy_bessel_and_discrepancy() -> None:
    """Verify Bessel J1 prediction alignment and 14.82x hierarchy tension."""
    spec = audit_spectral_hierarchy()
    assert len(spec["j1_zeros"]) == 4
    assert math.isclose(spec["j1_zeros"][0], 3.8317, abs_tol=1e-3)
    # Graviton predictions anchored at 380 GeV match claimed T700, T1000, T1322 within ~1%
    for rel in spec["graviton_rel_error_pct"]:
        assert abs(rel) < 1.0
    # Naive single-throat RS ratio vs observed 380/3600 exhibits >14x discrepancy
    assert math.isclose(spec["hierarchy_discrepancy_factor"], 14.82, abs_tol=0.1)


def test_look_elsewhere_trial_factors() -> None:
    """Verify exact Gaussian survival probabilities and CMS trial-factor dilution."""
    lee = audit_look_elsewhere()
    assert len(lee["cms_table2_audit"]) == len(CMS_EXO_21_010_TABLE2)
    assert math.isclose(lee["peak_8600_trials_factor"], 1139.4, abs_tol=2.0)
    assert math.isclose(lee["peak_3600_trials_factor"], 289.1, abs_tol=2.0)
    assert lee["peak_8600_width_adjusted_global_z"] < 1.6


def test_cross_experiment_poisson_replication() -> None:
    """Verify exact luminosity split and Poisson null-observation probability."""
    rep = audit_cross_experiment_replication()
    expected_total = 138.0 + ATLAS_PRD_108["lumi_fb_inv"] + CMS_PAS_EXO_25_004["lumi_fb_inv"]
    assert math.isclose(rep["lumi_total_fb"], expected_total, rel_tol=1e-6)
    assert math.isclose(rep["cms_run2_fraction"], 138.0 / 368.0, abs_tol=1e-4)
    assert math.isclose(rep["binomial_all_in_cms2_if_2_events"], (138.0 / 368.0) ** 2, abs_tol=1e-4)

    narrow_86 = rep["scenarios_8600GeV"][0]
    assert math.isclose(narrow_86["expected_followup_events_230fb"], 2.3 * 230.0 / 138.0, abs_tol=1e-2)
    assert math.isclose(narrow_86["p_zero_followup"], math.exp(-2.3 * 230.0 / 138.0), abs_tol=1e-4)


def test_bic_and_five_gates_reject_signal() -> None:
    """Verify Gate 2-5 failures and decisive negative Delta-BIC."""
    bic = audit_bic_and_five_gates()
    assert math.isclose(bic["two_delta_lnL_local"], 30.42, abs_tol=1e-2)
    assert math.isclose(bic["two_delta_lnL_global"], 7.40, abs_tol=1e-2)
    assert bic["delta_bic_local_k9"] < -20.0
    assert bic["delta_bic_global_k9"] < -40.0
    assert bic["all_gates_passed"] is False
    assert bic["five_gates"]["gate_1_resolvability"]["pass"] is True
    assert bic["five_gates"]["gate_2_amplitude_and_replication"]["pass"] is False
    assert bic["five_gates"]["gate_3_parameter_precision_and_hierarchy"]["pass"] is False
    assert bic["five_gates"]["gate_4_multiplicity_lee"]["pass"] is False
    assert bic["five_gates"]["gate_5_bayesian_bic"]["pass"] is False


def test_full_audit_artifacts_written(tmp_path: Path) -> None:
    """Verify end-to-end persistence of JSON metrics and figures."""
    summary = run_full_audit(tmp_path)
    assert (tmp_path / "results" / "kk_audit_summary.json").exists()
    assert (tmp_path / "figures" / "fig1_kk_tower_bayesian_audit.pdf").exists()
    assert (tmp_path / "figures" / "fig1_kk_tower_bayesian_audit.png").exists()
    assert len(summary["figures"]) == 2
