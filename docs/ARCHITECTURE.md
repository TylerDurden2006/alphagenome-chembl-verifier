# Architecture: Multi-Agent Hypothesis Generation + Deterministic Verification

`alphagenome-chembl-verifier` is designed as the **deterministic verification layer** inside agentic scientific discovery workflows (e.g., Google Antigravity `/run-workflow` and `/teamwork` pipelines powered by Gemini 4 Carbon).

## System Design Principle

Large-horizon agentic reasoning excels at synthesizing cross-disciplinary hypotheses from literature (`arXiv`, `OpenAlex`) and sequence-to-function models (`AlphaGenome`). However, multi-step synthesis can introduce:
1. Unresolved or non-canonical gene symbols,
2. Common benign variants misclassified as causal drivers, and
3. Targets lacking chemical matter or clinical tractability.

By decoupling **probabilistic hypothesis generation** (agentic layer) from **deterministic evidence gating** (`verify_target_evidence`), every nominated target is backed by auditable, reproducible API provenance:

```
[Agentic Layer: Google Antigravity + Gemini 4 Carbon]
   ├─ Query AlphaGenome Atlas (non-coding chromatin/splicing impact)
   ├─ Query gnomAD (population allele frequency constraint)
   └─ Propose candidate (variant_id, gene_symbol, impact_score, gnomad_af)
            │
            ▼
[Deterministic Layer: alphagenome-chembl-verifier]
   ├─ 1. Ensembl REST (rest.ensembl.org)        -> Canonical ENSG ID
   ├─ 2. OpenTargets GraphQL (platform)         -> Clinical / Small-Molecule Tractability
   ├─ 3. UniProtKB REST (rest.uniprot.org)      -> Reviewed Human Accession
   └─ 4. ChEMBL REST (ebi.ac.uk/chembl)         -> Sub-micromolar Bioactivity (<= 1000 nM)
            │
            ▼
[Audited Output: TIER_1_HIGH_CONFIDENCE | TIER_2_VALIDATED | UNVERIFIED]
```
