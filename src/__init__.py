"""Bayesian & Multiplicity Audit Engine for LHC Kaluza-Klein Tower Claims."""

from .kk_audit_engine import (
    ATLAS_PRD_108,
    CMS_EXO_21_010_TABLE2,
    CMS_PAS_EXO_25_004,
    LE_YAOUANC_RICHARD_2026,
    audit_bic_and_five_gates,
    audit_cross_experiment_replication,
    audit_look_elsewhere,
    audit_spectral_hierarchy,
    generate_audit_figure,
    run_full_audit,
)

__all__ = [
    "ATLAS_PRD_108",
    "CMS_EXO_21_010_TABLE2",
    "CMS_PAS_EXO_25_004",
    "LE_YAOUANC_RICHARD_2026",
    "audit_spectral_hierarchy",
    "audit_look_elsewhere",
    "audit_cross_experiment_replication",
    "audit_bic_and_five_gates",
    "generate_audit_figure",
    "run_full_audit",
]
