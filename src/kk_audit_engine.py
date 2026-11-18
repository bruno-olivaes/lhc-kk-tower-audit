"""
Bayesian and Multiplicity Audit of LHC Kaluza-Klein Multi-Resonance Claims.

Evaluates Le Yaouanc & Richard (arXiv:2610.31458v3; arXiv:2606.29922v3) against
primary CMS (arXiv:2206.09997v3 / JHEP 07 (2026) 155 Addendum; CMS-PAS-EXO-25-004)
and ATLAS (arXiv:2307.14944v2 / Phys. Rev. D 108 (2023) 112005) paired-dijet results
using The Matrioska Protocol's five integrity gates.
"""

from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.special import jn_zeros
from scipy.stats import norm, poisson

# ---------------------------------------------------------------------------
# Primary Verbatim Experimental & Phenomenological Constants
# ---------------------------------------------------------------------------

CMS_EXO_21_010_TABLE2: list[dict[str, float]] = [
    # Source: CMS Collaboration, arXiv:2206.09997v3 Table 2 (p. 35, 138 fb^-1 @ 13 TeV)
    {"mY_TeV": 8.6, "mX_TeV": 2.15, "alpha": 0.25, "width_pct": 1.5, "sigma_BA_pb": 1.68e-5, "events": 2.3, "z_local": 3.9, "z_global": 1.6},
    {"mY_TeV": 8.6, "mX_TeV": 2.15, "alpha": 0.25, "width_pct": 5.0, "sigma_BA_pb": 2.04e-5, "events": 2.8, "z_local": 3.8, "z_global": 1.5},
    {"mY_TeV": 8.6, "mX_TeV": 2.15, "alpha": 0.25, "width_pct": 10.0, "sigma_BA_pb": 2.22e-5, "events": 3.1, "z_local": 3.6, "z_global": 1.4},
    {"mY_TeV": 3.6, "mX_TeV": 1.00, "alpha": 0.29, "width_pct": 1.5, "sigma_BA_pb": 0.83e-3, "events": 114.0, "z_local": 3.0, "z_global": 0.5},
    {"mY_TeV": 3.6, "mX_TeV": 1.00, "alpha": 0.29, "width_pct": 5.0, "sigma_BA_pb": 1.18e-3, "events": 162.0, "z_local": 3.6, "z_global": 1.6},
    {"mY_TeV": 3.6, "mX_TeV": 1.00, "alpha": 0.29, "width_pct": 10.0, "sigma_BA_pb": 1.53e-3, "events": 210.0, "z_local": 3.9, "z_global": 2.2},
]

CMS_EXPERIMENT_META: dict[str, Any] = {
    "lumi_fb_inv": 138.0,
    "sqrt_s_TeV": 13.0,
    "alpha_bins": 13,
    "bg_only_combined_p": 0.18,
    "gold_plated_8TeV_events": 2,
    "third_high_mass_event_m4j_TeV": 5.8,
    "nonresonant_mX_TeV": 0.95,
    "nonresonant_z_local": 3.6,
    "nonresonant_z_global": 2.5,
}

ATLAS_PRD_108: dict[str, Any] = {
    # Source: ATLAS Collaboration, Phys. Rev. D 108 (2023) 112005 (arXiv:2307.14944v2)
    "lumi_fb_inv": 140.0,
    "sqrt_s_TeV": 13.0,
    "bg_fit_p_min": 0.74,
    "bg_fit_p_max": 1.00,
    "max_m4j_observed_TeV": 6.6,
    "max_m4j_m2j_avg_TeV": 2.2,
    "gold_plated_8TeV_events": 0,
    "max_bumphunter_z_global_m4j": 0.53,
    "max_bumphunter_m4j_GeV": 3200.0,
    "max_bumphunter_z_global_m2j": 1.98,
    "max_bumphunter_m2j_GeV": 800.0,
}

CMS_PAS_EXO_25_004: dict[str, Any] = {
    # Partial Run 3 luminosity (~90 fb^-1 at 13.6 TeV) as cited in arXiv:2610.31458v3
    "lumi_fb_inv": 90.0,
    "sqrt_s_TeV": 13.6,
    "gold_plated_8TeV_events": 0,
}

LE_YAOUANC_RICHARD_2026: dict[str, Any] = {
    "claimed_graviton_GeV": [380.0, 700.0, 1000.0, 1322.0],
    "claimed_gluon_GeV": [3600.0, 8600.0],
    "claimed_quark_GeV": [1040.0, 2150.0],
    "rs_vector_roots": [2.45, 5.57, 8.70, 11.84],  # Quoted in arXiv:2610.31458v3 Sec. II.2
    "posthoc_top_suppression_factor": (5.0, 10.0),
    "posthoc_extra_branes": 2,
}


def audit_spectral_hierarchy() -> dict[str, Any]:
    """Compare claimed KK graviton and gluon spectra against Bessel J_1 and RS roots."""
    j1_zeros = jn_zeros(1, 4)
    g1_mass = LE_YAOUANC_RICHARD_2026["claimed_graviton_GeV"][0]
    predicted_graviton = [float(g1_mass * z / j1_zeros[0]) for z in j1_zeros]
    claimed_graviton = LE_YAOUANC_RICHARD_2026["claimed_graviton_GeV"]
    graviton_rel_err_pct = [
        100.0 * (obs - pred) / pred
        for obs, pred in zip(claimed_graviton, predicted_graviton)
    ]

    rs_vec = LE_YAOUANC_RICHARD_2026["rs_vector_roots"]
    gluon_obs = LE_YAOUANC_RICHARD_2026["claimed_gluon_GeV"]
    gluon_ratio_obs = gluon_obs[1] / gluon_obs[0]
    gluon_ratio_rs = rs_vec[1] / rs_vec[0]

    # Single-throat RS ratio of lightest KK graviton (x_G = 3.8317) to KK gauge/gluon (2.45)
    naive_grav_to_gluon = float(j1_zeros[0] / rs_vec[0])
    obs_grav_to_gluon = claimed_graviton[0] / gluon_obs[0]
    hierarchy_discrepancy = naive_grav_to_gluon / obs_grav_to_gluon

    return {
        "j1_zeros": [float(x) for x in j1_zeros],
        "claimed_graviton_GeV": claimed_graviton,
        "bessel_j1_predicted_GeV": [round(x, 2) for x in predicted_graviton],
        "graviton_rel_error_pct": [round(x, 3) for x in graviton_rel_err_pct],
        "gluon_ratio_observed": round(gluon_ratio_obs, 4),
        "gluon_ratio_rs_predicted": round(gluon_ratio_rs, 4),
        "naive_rs_graviton_to_gluon_ratio": round(naive_grav_to_gluon, 4),
        "observed_graviton_to_gluon_ratio": round(obs_grav_to_gluon, 4),
        "hierarchy_discrepancy_factor": round(hierarchy_discrepancy, 2),
    }


def audit_look_elsewhere() -> dict[str, Any]:
    """Quantify Look-Elsewhere Effect trial factors across CMS Table 2 hypotheses."""
    rows: list[dict[str, float]] = []
    for entry in CMS_EXO_21_010_TABLE2:
        p_loc = float(norm.sf(entry["z_local"]))
        p_glob = float(norm.sf(entry["z_global"]))
        n_trials_ratio = p_glob / p_loc
        n_trials_poisson = -math.log(max(1e-12, 1.0 - p_glob)) / p_loc
        rows.append({
            **entry,
            "p_local": p_loc,
            "p_global": p_glob,
            "n_trials_ratio": round(n_trials_ratio, 1),
            "n_trials_poisson": round(n_trials_poisson, 1),
        })

    best_86 = rows[0]  # 8.6 TeV, 1.5% width: 3.9s -> 1.6s
    best_36 = rows[5]  # 3.6 TeV, 10% width: 3.9s -> 2.2s

    # Conservative Bonferroni adjustment across selecting 2 distinct (alpha, width) peaks post-hoc
    # Note: CMS global p-value already scans M_Y and alpha within fixed width; choosing the best width
    # across {1.5%, 5%, 10%} adds an effective width-scan factor conservatively bounded by 3.
    p_comb_bonferroni = min(1.0, best_86["p_global"] * best_36["p_global"])  # naive if independent
    z_86_width_adj = float(norm.isf(min(0.9999, 1.0 - (1.0 - best_86["p_global"]) ** 3)))
    z_36_width_adj = float(norm.isf(min(0.9999, 1.0 - (1.0 - best_36["p_global"]) ** 3)))

    return {
        "cms_table2_audit": rows,
        "peak_8600_trials_factor": best_86["n_trials_ratio"],
        "peak_3600_trials_factor": best_36["n_trials_ratio"],
        "peak_8600_width_adjusted_global_z": round(z_86_width_adj, 3),
        "peak_3600_width_adjusted_global_z": round(z_36_width_adj, 3),
        "naive_independent_p_product": p_comb_bonferroni,
    }


def audit_cross_experiment_replication() -> dict[str, Any]:
    """
    Evaluate Poisson & binomial replication tension of the 8.6 TeV excess across
    CMS Run 2 (138 fb^-1), ATLAS Run 2 (140 fb^-1), and CMS Run 3 partial (90 fb^-1).
    """
    l_cms2 = CMS_EXPERIMENT_META["lumi_fb_inv"]
    l_atlas2 = ATLAS_PRD_108["lumi_fb_inv"]
    l_cms3 = CMS_PAS_EXO_25_004["lumi_fb_inv"]
    l_follow = l_atlas2 + l_cms3
    l_total = l_cms2 + l_follow

    k_cms2 = CMS_EXPERIMENT_META["gold_plated_8TeV_events"]
    k_follow = ATLAS_PRD_108["gold_plated_8TeV_events"] + CMS_PAS_EXO_25_004["gold_plated_8TeV_events"]

    # Conditional binomial probability that all k observed high-mass events land in CMS Run 2
    p_lumi_cms2 = l_cms2 / l_total
    prob_all_in_cms2_n2 = p_lumi_cms2 ** k_cms2
    prob_all_in_cms2_n3 = p_lumi_cms2 ** (k_cms2 + 1)

    scenarios: list[dict[str, float]] = []
    for row in CMS_EXO_21_010_TABLE2[:3]:  # Three width hypotheses at 8.6 TeV
        n_cms2_fit = row["events"]
        mu_follow = n_cms2_fit * (l_follow / l_cms2)
        p_zero_follow = float(poisson.pmf(k_follow, mu_follow))
        dilution_ratio = l_cms2 / l_total
        combined_sigma_ba = row["sigma_BA_pb"] * dilution_ratio
        scenarios.append({
            "width_pct": row["width_pct"],
            "cms_run2_best_fit_events": n_cms2_fit,
            "expected_followup_events_230fb": round(mu_follow, 3),
            "p_zero_followup": round(p_zero_follow, 5),
            "cms_sigma_BA_pb": row["sigma_BA_pb"],
            "lumi_diluted_combined_sigma_BA_pb": combined_sigma_ba,
        })

    return {
        "lumi_cms_run2_fb": l_cms2,
        "lumi_atlas_run2_fb": l_atlas2,
        "lumi_cms_run3_partial_fb": l_cms3,
        "lumi_followup_fb": l_follow,
        "lumi_total_fb": l_total,
        "cms_run2_fraction": round(p_lumi_cms2, 4),
        "binomial_all_in_cms2_if_2_events": round(prob_all_in_cms2_n2, 4),
        "binomial_all_in_cms2_if_3_events": round(prob_all_in_cms2_n3, 4),
        "scenarios_8600GeV": scenarios,
    }


def _compute_delta_bic(two_delta_ln_l: float, k_params: int, n_bins: int) -> float:
    """Compute Delta BIC = BIC_null - BIC_signal = 2*Delta_lnL - k*ln(N)."""
    return two_delta_ln_l - k_params * math.log(n_bins)


def audit_bic_and_five_gates() -> dict[str, Any]:
    """
    Apply The Matrioska Protocol's 5 integrity gates and BIC penalty for post-hoc
    multi-brane parameter inflation against the Standard Model null.
    """
    # Conservative effective bin count across 13 alpha bins x 30 active tetrajet resolution bins
    n_bins = 390
    ln_n = math.log(n_bins)

    # Best local and LEE-adjusted global test statistics from the two CMS Table 2 peaks
    two_dlnl_local = 3.9**2 + 3.9**2   # 30.42
    two_dlnl_global = 1.6**2 + 2.2**2  # 7.40

    # Parameter accounting for Le Yaouanc & Richard multi-brane post-hoc KK fit:
    # Base 9: m_G1, m_g1, extra-brane hierarchy, c_t suppression, gg/gamma suppression, alpha_1, alpha_2, Gamma_1, Gamma_2
    # Extended 11: + c_eR, mixing angle sin(theta')
    k_base = 9
    k_ext = 11

    delta_bic_local_k9 = _compute_delta_bic(two_dlnl_local, k_base, n_bins)
    delta_bic_local_k11 = _compute_delta_bic(two_dlnl_local, k_ext, n_bins)
    delta_bic_global_k9 = _compute_delta_bic(two_dlnl_global, k_base, n_bins)
    delta_bic_global_k11 = _compute_delta_bic(two_dlnl_global, k_ext, n_bins)

    # Parameter sensitivity sweep k in [2..11], N in [100..600]
    k_grid = list(range(2, 12))
    n_grid = [100, 200, 390, 600]
    bic_grid = {
        str(n): [round(_compute_delta_bic(two_dlnl_local, k, n), 2) for k in k_grid]
        for n in n_grid
    }

    lee = audit_look_elsewhere()
    rep = audit_cross_experiment_replication()

    gates = {
        "gate_1_resolvability": {
            "pass": True,
            "summary": "Jet quartet masses (3.6 TeV, 8.6 TeV; 1.0 TeV, 2.15 TeV) lie cleanly above jet-resolution thresholds.",
        },
        "gate_2_amplitude_and_replication": {
            "pass": False,
            "summary": (
                f"Zero gold-plated 8.6 TeV events in 230 fb^-1 of ATLAS Run 2 + CMS Run 3 follow-up "
                f"(P(0 | mu={rep['scenarios_8600GeV'][0]['expected_followup_events_230fb']})="
                f"{rep['scenarios_8600GeV'][0]['p_zero_followup']:.4f}; combined rate diluted by {rep['cms_run2_fraction']:.3f}x)."
            ),
        },
        "gate_3_parameter_precision_and_hierarchy": {
            "pass": False,
            "summary": (
                "Single-throat RS graviton-to-gluon ratio fails by 14.82x (0.1056 obs vs 1.5640 pred), "
                "requiring two post-hoc branes, 5-10x top suppression, and independent alpha/width tuning."
            ),
        },
        "gate_4_multiplicity_lee": {
            "pass": False,
            "summary": (
                f"CMS Look-Elsewhere drops 8.6 TeV from 3.9s to 1.6s (p_global=0.0548, trials~{lee['peak_8600_trials_factor']:.0f}) "
                f"and 3.6 TeV from 3.9s to 2.2s (p_global=0.0139); width-selection drops 8.6 TeV further to {lee['peak_8600_width_adjusted_global_z']}s."
            ),
        },
        "gate_5_bayesian_bic": {
            "pass": False,
            "summary": (
                f"Delta-BIC_local = {delta_bic_local_k9:.2f} (k=9, N=390) and Delta-BIC_global = {delta_bic_global_k9:.2f}; "
                "both severely violate the +10.0 threshold."
            ),
        },
    }

    return {
        "n_bins": n_bins,
        "ln_n": round(ln_n, 4),
        "two_delta_lnL_local": round(two_dlnl_local, 2),
        "two_delta_lnL_global": round(two_dlnl_global, 2),
        "delta_bic_local_k9": round(delta_bic_local_k9, 2),
        "delta_bic_local_k11": round(delta_bic_local_k11, 2),
        "delta_bic_global_k9": round(delta_bic_global_k9, 2),
        "delta_bic_global_k11": round(delta_bic_global_k11, 2),
        "bic_local_sweep": {"k_grid": k_grid, "by_n_bins": bic_grid},
        "five_gates": gates,
        "all_gates_passed": all(g["pass"] for g in gates.values()),
    }


def generate_audit_figure(out_dir: str | Path) -> list[str]:
    """Generate 4-panel publication-grade summary figure in PDF and PNG formats."""
    out_path = Path(out_dir)
    out_path.mkdir(parents=True, exist_ok=True)

    spec = audit_spectral_hierarchy()
    rep = audit_cross_experiment_replication()
    bic = audit_bic_and_five_gates()

    plt.rcParams.update({
        "font.size": 9.5,
        "axes.titlesize": 10.5,
        "axes.titleweight": "bold",
        "axes.labelsize": 9.5,
        "legend.fontsize": 8.0,
    })

    fig, axes = plt.subplots(2, 2, figsize=(11.5, 8.6), constrained_layout=True)

    # Panel (a): Graviton Bessel J1 scaling vs Single-Throat Gluon Contradiction
    ax = axes[0, 0]
    modes = np.arange(1, 5)
    ax.plot(modes, spec["bessel_j1_predicted_GeV"], "o--", color="#2563eb", label="Bessel $J_1$ pred. (anchored at 380 GeV)")
    ax.plot(modes, spec["claimed_graviton_GeV"], "s", markersize=7, color="#0f172a", label="Claimed KK gravitons (380-1322 GeV)")
    naive_g1_pred = 380.0 / spec["naive_rs_graviton_to_gluon_ratio"]
    ax.axhline(3600.0, color="#dc2626", linestyle="-", linewidth=1.6, label="Claimed $g_{KK1}$ (3600 GeV)")
    ax.axhline(naive_g1_pred, color="#f59e0b", linestyle=":", linewidth=1.8, label=f"Naive 1-throat RS $g_{{KK1}}$ ({naive_g1_pred:.0f} GeV)")
    ax.annotate("", xy=(2.5, 3600.0), xytext=(2.5, naive_g1_pred),
                arrowprops=dict(arrowstyle="<->", color="#dc2626", lw=1.3, alpha=0.8))
    ax.text(2.58, 2250, "14.82x hierarchy gap\n(needs 2 extra branes)",
            fontsize=8.2, color="#991b1b", bbox=dict(boxstyle="round,pad=0.3", fc="#fee2e2", ec="#dc2626"))
    ax.set_xticks(modes)
    ax.set_xlabel("Kaluza-Klein Radial Mode $n$")
    ax.set_ylabel("Mass Scale [GeV]")
    ax.set_title("(a) Bessel $J_1$ Graviton Fit vs. Gluon Hierarchy Gap")
    ax.set_ylim(0, 4100)
    ax.legend(loc="upper left", frameon=True)
    ax.grid(alpha=0.25)

    # Panel (b): CMS Table 2 Local vs Global Significance
    ax = axes[0, 1]
    labels = ["8.6 TeV\n(1.5%)", "8.6 TeV\n(5%)", "8.6 TeV\n(10%)", "3.6 TeV\n(1.5%)", "3.6 TeV\n(5%)", "3.6 TeV\n(10%)"]
    x = np.arange(len(labels))
    w = 0.35
    z_loc = [r["z_local"] for r in CMS_EXO_21_010_TABLE2]
    z_glob = [r["z_global"] for r in CMS_EXO_21_010_TABLE2]
    ax.bar(x - w / 2, z_loc, width=w, color="#94a3b8", label="CMS Local Significance ($Z_{\\mathrm{local}}$)")
    ax.bar(x + w / 2, z_glob, width=w, color="#2563eb", label="CMS Global Significance ($Z_{\\mathrm{global}}$, LEE)")
    ax.axhline(3.0, color="#f59e0b", linestyle="--", linewidth=1.2, label="3$\\sigma$ Evidence Threshold")
    ax.axhline(5.0, color="#dc2626", linestyle=":", linewidth=1.2, label="5$\\sigma$ Discovery Threshold")
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=8)
    ax.set_ylabel("Significance [$\\sigma$]")
    ax.set_title("(b) Look-Elsewhere Attenuation Across Widths (CMS Table 2)")
    ax.set_ylim(0, 5.5)
    ax.legend(loc="upper right", frameon=True)
    ax.grid(axis="y", alpha=0.25)

    # Panel (c): Cross-Experiment Poisson Zero-Observation Tension
    ax = axes[1, 0]
    mu_vals = np.linspace(0.5, 6.0, 200)
    p_zero_curve = np.exp(-mu_vals)
    ax.plot(mu_vals, p_zero_curve * 100.0, color="#0f172a", linewidth=2.0, label="$P(N_{\\mathrm{follow}}=0 \\mid \\mu) = e^{-\\mu}$")
    colors = ["#2563eb", "#7c3aed", "#dc2626"]
    for sc, col in zip(rep["scenarios_8600GeV"], colors):
        mu = sc["expected_followup_events_230fb"]
        p0 = sc["p_zero_followup"] * 100.0
        ax.scatter([mu], [p0], color=col, s=55, zorder=5,
                   label=f"$\\Gamma/M={sc['width_pct']:.1f}\\%$: $\\mu={mu:.2f}$, $P(0)={p0:.2f}\\%$")
    ax.axhline(5.0, color="#dc2626", linestyle="--", linewidth=1.1, label="$p = 5\\%$ threshold")
    ax.set_xlabel("Expected 8.6 TeV Events in 230 $\\mathrm{fb}^{-1}$ Follow-up (ATLAS + CMS Run 3)")
    ax.set_ylabel("Null Observation Probability $P(N=0)$ [%]")
    ax.set_title("(c) Cross-Experiment Replication Tension ($k_{\\mathrm{follow}} = 0$)")
    ax.set_ylim(0, 65)
    ax.legend(loc="upper right", frameon=True)
    ax.grid(alpha=0.25)

    # Panel (d): Delta BIC Phase Diagram
    ax = axes[1, 1]
    k_grid = np.array(bic["bic_local_sweep"]["k_grid"])
    for n_str, col in [("100", "#10b981"), ("390", "#2563eb"), ("600", "#7c3aed")]:
        ax.plot(k_grid, bic["bic_local_sweep"]["by_n_bins"][n_str], "o-", color=col,
                label=f"Local $2\\Delta\\ln\\mathcal{{L}}=30.42$ ($N={n_str}$ bins)")
    k_eff = np.array(k_grid)
    ax.plot(k_eff, [bic["two_delta_lnL_global"] - k * math.log(390) for k in k_eff], "s--", color="#dc2626",
            label="LEE-adjusted $2\\Delta\\ln\\mathcal{L}=7.40$ ($N=390$ bins)")
    ax.axhline(10.0, color="#16a34a", linestyle="-", linewidth=1.5, label="Matrioska Gate 5 Target ($\\Delta\\mathrm{BIC} > +10$)")
    ax.axhline(0.0, color="#64748b", linestyle=":", linewidth=1.0)
    ax.scatter([9], [bic["delta_bic_local_k9"]], color="#0f172a", s=70, zorder=6)
    ax.annotate(f"$k=9$: $\\Delta\\mathrm{{BIC}}={bic['delta_bic_local_k9']:.1f}$",
                xy=(9, bic["delta_bic_local_k9"]), xytext=(7.2, 14.5),
                arrowprops=dict(arrowstyle="->", color="#0f172a", lw=1.2), fontsize=8.5, fontweight="bold",
                bbox=dict(boxstyle="round,pad=0.25", fc="#ffffff", ec="#94a3b8", alpha=0.95))
    ax.set_xlabel("Post-Hoc Tuned Parameter Count $k$")
    ax.set_ylabel("$\\Delta\\mathrm{BIC} = \\mathrm{BIC}_{\\mathrm{null}} - \\mathrm{BIC}_{\\mathrm{KK}}$")
    ax.set_title("(d) Bayesian Over-Parametrization Penalty (Gate 5)")
    ax.set_ylim(-65, 25)
    ax.legend(loc="lower left", fontsize=7.5, frameon=True)
    ax.grid(alpha=0.25)

    saved: list[str] = []
    for ext in ("pdf", "png"):
        target = out_path / f"fig1_kk_tower_bayesian_audit.{ext}"
        fig.savefig(target, dpi=220)
        saved.append(str(target))
    plt.close(fig)
    return saved


def run_full_audit(base_dir: str | Path | None = None) -> dict[str, Any]:
    """Run all audit modules and persist JSON + figure artifacts when base_dir is supplied."""
    summary = {
        "spectral_hierarchy": audit_spectral_hierarchy(),
        "look_elsewhere": audit_look_elsewhere(),
        "cross_experiment_replication": audit_cross_experiment_replication(),
        "bic_and_five_gates": audit_bic_and_five_gates(),
    }
    if base_dir is not None:
        root = Path(base_dir)
        results_dir = root / "results"
        results_dir.mkdir(parents=True, exist_ok=True)
        with open(results_dir / "kk_audit_summary.json", "w", encoding="utf-8") as f:
            json.dump(summary, f, indent=2)
        summary["figures"] = generate_audit_figure(root / "figures")
    return summary


if __name__ == "__main__":
    repo_root = Path(__file__).resolve().parents[1]
    res = run_full_audit(repo_root)
    print(json.dumps(res["bic_and_five_gates"]["five_gates"], indent=2))
    print("Delta BIC (local, k=9):", res["bic_and_five_gates"]["delta_bic_local_k9"])
    print("Delta BIC (global, k=9):", res["bic_and_five_gates"]["delta_bic_global_k9"])
    print("All gates passed:", res["bic_and_five_gates"]["all_gates_passed"])
