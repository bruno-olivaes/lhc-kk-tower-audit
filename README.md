# LHC Kaluza–Klein Tower Bayesian & Multiplicity Audit (Paper 8)

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![ORCID](https://img.shields.io/badge/ORCID-0009--0003--5351--446X-green.svg)](https://orcid.org/0009-0003-5351-446X)
[![Portfolio](https://img.shields.io/badge/Portfolio-The%20Matrioska%20Protocol-blueviolet)](https://bruno-olivaes.github.io)

**Companion Quantitative Audit to *The Matrioska Protocol* (`arXiv:2610.31458v3` vs. Primary CERN CMS & ATLAS Results)**  
**Author:** Bruno Luçardo Olivaes (`brunoolivaes@gmail.com`)

---

## Overview

This repository audits recent phenomenological claims of a multi-sector Randall–Sundrum Kaluza–Klein (KK) spectrum at the LHC—spanning KK graviton candidates ($380\text{--}1322\text{ GeV}$) and heavy KK gluons $g_{\mathrm{KK}}(3.6\text{ TeV}) \to q_{\mathrm{KK}}(1.0\text{ TeV})$ and $g_{\mathrm{KK}}(8.6\text{ TeV}) \to q_{\mathrm{KK}}(2.15\text{ TeV})$ ([`arXiv:2610.31458v3`](https://arxiv.org/abs/2610.31458), [`arXiv:2606.29922v3`](https://arxiv.org/abs/2606.29922))—against verbatim experimental figures from:

- **CMS Collaboration:** [`arXiv:2206.09997v3`](https://arxiv.org/abs/2206.09997) (*JHEP 07 (2023) 161* + Addendum *JHEP 07 (2026) 155*, $138\text{ fb}^{-1}$, Table 2)
- **ATLAS Collaboration:** [`arXiv:2307.14944v2`](https://arxiv.org/abs/2307.14944) (*Phys. Rev. D 108 (2023) 112005*, $140\text{ fb}^{-1}$)
- **CMS Run 3 Partial:** `CMS-PAS-EXO-25-004` ($\sim 90\text{ fb}^{-1}$ @ $13.6\text{ TeV}$, cited in `arXiv:2610.31458v3`)

---

## Key Audited Results (5-Gate Protocol)

| Integrity Gate | Verdict | Audited Metric |
|:---|:---:|:---|
| **Gate 1 — Resolvability** | PASS | $3.6\text{ TeV}$ and $8.6\text{ TeV}$ masses reside above detector jet-resolution thresholds. |
| **Gate 2 — Cross-Experiment Replication** | **FAIL** | $0$ gold-plated $8.6\text{ TeV}$ events across $230\text{ fb}^{-1}$ follow-up ($P(0 \mid \mu = 3.833) = 0.0216$; $0.375\times$ dilution). |
| **Gate 3 — Hierarchy & Economy** | **FAIL** | Naive 1-throat RS graviton/gluon ratio disagrees by **$14.82\times$** ($0.1056$ obs vs. $1.5640$ pred), driving post-hoc brane & BR tunings. |
| **Gate 4 — Global Multiplicity (LEE)** | **FAIL** | CMS $Z_{\mathrm{local}} = 3.9\sigma \to Z_{\mathrm{global}} = 1.6\sigma$ ($8.6\text{ TeV}$, $N_{\mathrm{trials}} \approx 1{,}139$) and $2.2\sigma$ ($3.6\text{ TeV}$, $N_{\mathrm{trials}} \approx 289$). |
| **Gate 5 — Bayesian Information ($\Delta\mathrm{BIC}$)** | **FAIL** | $\Delta\mathrm{BIC}_{\mathrm{local}}(k=9) = \mathbf{-23.28}$, $\Delta\mathrm{BIC}_{\mathrm{global}}(k=9) = \mathbf{-46.30}$ (target: $> +10.0$). |

---

## Reproducibility

```bash
python src/kk_audit_engine.py
python -m pytest -v
```

Generated artifacts:
- `results/kk_audit_summary.json`
- `figures/fig1_kk_tower_bayesian_audit.pdf` / `.png`
- `manuscript/Paper8_LHC_KK_Tower_Audit.pdf`
