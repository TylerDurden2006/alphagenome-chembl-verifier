# Deterministic Cross-Verification of Non-Coding Regulatory Variants and Chemogenomic Target Tractability

**Author:** Rishi Mehra  
**Affiliation:** Department of Computer Science and Engineering, I.K. Gujral Punjab Technical University, Kapurthala, Punjab, India  
**Repository:** https://github.com/TylerDurden2006/alphagenome-chembl-verifier  
**Date:** November 16, 2026  

---

## Abstract

Prioritizing causal non-coding regulatory variants and mapping them to druggable protein targets remains bottlenecked by fragmented, manual reconciliation across genomic, population-genetic, and chemogenomic databases. We present `alphagenome-chembl-verifier`, a lightweight, deterministic Python verification pipeline that enforces reproducible evidence gates across regulatory variant impact scores (AlphaGenome), population allele frequency rarity constraints (gnomAD), reviewed protein target resolution (UniProtKB), and sub-micromolar ligand bioactivity (ChEMBL). By separating deterministic verification rules from upstream hypothesis generation, the toolkit reduces false-positive target nominations in dry-lab multi-agent discovery workflows.

## 1. Introduction

Recent advances in sequence-to-function models such as AlphaGenome enable high-throughput prediction of chromatin accessibility, splicing disruption, and transcriptional regulation for non-coding variants. However, translating a high predicted regulatory score into an actionable therapeutic hypothesis requires verifying that: (1) the variant is sufficiently rare in human populations (gnomAD), (2) the affected gene resolves to a reviewed human protein entry (UniProtKB), and (3) the target exhibits tractable chemical matter (ChEMBL).

## 2. Methodology

Given a normalized `VariantEvidence` record $(v, g, s_{\text{impact}}, f_{\text{AF}})$, `verify_target_evidence` applies four deterministic checks:
1. **Input Bounds Validation:** Verifies $s_{\text{impact}} \in [0, 1]$ and $f_{\text{AF}} \in [0, 1]$.
2. **Regulatory & Rarity Gate:** Requires $s_{\text{impact}} \ge 0.70$ and $f_{\text{AF}} \le 0.01$.
3. **Target Resolution:** Requires a reviewed human UniProtKB primary accession via `rest.uniprot.org`.
4. **Chemogenomic Tractability:** Queries ChEMBL bioactivity records (`standard_units=nM`) and counts ligands with $\text{IC}_{50}, K_i, \text{ or } K_d \le 1000\text{ nM}$. Targets passing all gates with $\ge 3$ potent ligands are classified as `TIER_1_HIGH_CONFIDENCE`.

## 3. Implementation & Reproducibility

The package is implemented in pure Python 3.10+ using the standard library (`urllib.request`, `dataclasses`, `unittest`) with zero external runtime dependencies, ensuring long-term reproducibility and straightforward integration into agentic research environments such as Google Antigravity.

## 4. Availability

Source code, automated tests, and `CITATION.cff` metadata are openly available under the MIT License at `https://github.com/TylerDurden2006/alphagenome-chembl-verifier`.
