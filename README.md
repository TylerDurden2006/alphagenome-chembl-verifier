# alphagenome-chembl-verifier

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![Python 3.10+](https://img.shields.io/badge/python-3.10%2B-blue.svg)](pyproject.toml)
[![Zero Dependencies](https://img.shields.io/badge/dependencies-stdlib%20only-green.svg)](pyproject.toml)

Deterministic cross-verification pipeline linking non-coding regulatory variant impact predictions (**AlphaGenome**) with population allele frequency constraints (**gnomAD**), canonical gene mappings (**Ensembl**), clinical tractability (**OpenTargets**), reviewed protein annotations (**UniProtKB**), and sub-micromolar bioactivity evidence (**ChEMBL**).

## Why This Exists

Multi-agent scientific discovery workflows (such as **Google Antigravity** pipelines powered by **Gemini 4 Carbon**) can synthesize variant-to-target hypotheses at scale, but unverified nominations frequently fail due to:
1. Unresolved or non-canonical gene symbols,
2. Common benign variants misclassified as causal drivers (`gnomAD AF > 0.01`), or
3. Targets lacking chemical matter or clinical tractability.

`alphagenome-chembl-verifier` enforces **deterministic, auditable evidence gates** across four live public bioinformatics APIs (`rest.ensembl.org`, `api.platform.opentargets.org`, `rest.uniprot.org`, and `ebi.ac.uk/chembl/api/data`) with zero external runtime dependencies.

## Architecture

See [docs/ARCHITECTURE.md](docs/ARCHITECTURE.md) and [docs/preprint.md](docs/preprint.md) for full methodology.

```
Variant + AlphaGenome Impact Score + gnomAD AF
        │
        ▼
[1] Schema & Rarity Gate (impact >= 0.70, gnomAD AF <= 0.01)
        │
        ├──▶ [2] Ensembl REST (rest.ensembl.org) -> Canonical ENSG ID
        │           └──▶ [3] OpenTargets GraphQL -> Clinical & Small-Molecule Tractability
        │
        └──▶ [4] UniProtKB REST (rest.uniprot.org) -> Reviewed Human Accession
                    └──▶ [5] ChEMBL REST (ebi.ac.uk/chembl) -> Bioactivity <= 1000 nM
        │
        ▼
VerificationReport (TIER_1_HIGH_CONFIDENCE | TIER_2_VALIDATED | UNVERIFIED)
```

## Quick Start

```bash
# Run unit tests (Python 3.10+ standard library)
python -m unittest discover -s tests -v

# Single-variant live verification
PYTHONPATH=src python -m ag_chembl_verifier.cli --variant chr12:57751011:C>T --gene CDK4 --impact 0.84 --gnomad-af 0.0003

# Batch 10-locus benchmark (live Ensembl + OpenTargets + UniProtKB + ChEMBL)
PYTHONPATH=src python -m ag_chembl_verifier.cli --batch data/benchmark_loci.json --format markdown
```

## Observed Live Benchmark Output (`results/benchmark_report.md`)

Executed against live `rest.ensembl.org`, `api.platform.opentargets.org`, `rest.uniprot.org`, and `ebi.ac.uk/chembl/api/data` (`limit=25` activities per target):

| Variant | Gene | Ensembl | UniProt | Passed | Tier | Potent Ligands | OpenTargets Flags | Notes |
|---|---|---|---|---|---|---:|---|---|
| `chr7:140753336:A>T` | `BRAF` | `ENSG00000157764` | `P15056` | FAIL | `UNVERIFIED` | 0 | SM:Approved Drug, SM:Structure with Ligand | No ChEMBL bioactivity <= 1000.0 nM observed |
| `chr7:55174772:G>A` | `EGFR` | `ENSG00000146648` | `P00533` | FAIL | `UNVERIFIED` | 0 | SM:Approved Drug, SM:Structure with Ligand | No ChEMBL bioactivity <= 1000.0 nM observed |
| `chr1:55039974:G>T` | `PCSK9` | `ENSG00000169174` | `Q8NBP7` | FAIL | `UNVERIFIED` | 0 | SM:Structure with Ligand, SM:High-Quality Ligand | No ChEMBL bioactivity <= 1000.0 nM observed |
| `chr18:63318364:C>T` | `BCL2` | `ENSG00000171791` | `P10415` | FAIL | `UNVERIFIED` | 0 | SM:Approved Drug, SM:Structure with Ligand | No ChEMBL bioactivity <= 1000.0 nM observed |
| `chr5:1295113:G>A` | `TERT` | `ENSG00000164362` | `O14746` | PASS | `TIER_2_VALIDATED` | 1 | SM:High-Quality Ligand, SM:Druggable Family | All gates passed |
| `chr12:57751011:C>T` | `CDK4` | `ENSG00000135446` | `P11802` | PASS | `TIER_1_HIGH_CONFIDENCE` | 13 | SM:Approved Drug, SM:Structure with Ligand | All gates passed |
| `chr3:179234297:A>G` | `PIK3CA` | `ENSG00000121879` | `P42336` | PASS | `TIER_1_HIGH_CONFIDENCE` | 12 | SM:Approved Drug, SM:Structure with Ligand | All gates passed |
| `chr9:5073770:G>T` | `JAK2` | `ENSG00000096968` | `O60674` | PASS | `TIER_1_HIGH_CONFIDENCE` | 15 | SM:Approved Drug, SM:Structure with Ligand | All gates passed |
| `chr17:7674220:C>T` | `TP53` | `ENSG00000141510` | `P04637` | FAIL | `UNVERIFIED` | 7 | SM:Advanced Clinical, SM:Structure with Ligand | gnomAD AF 0.1200 > 0.0100 rarity cutoff |
| `chr19:44908684:T>C` | `APOE` | `ENSG00000130203` | `P02649` | FAIL | `UNVERIFIED` | 0 | SM:Structure with Ligand, AB:UniProt loc high conf | Impact score 0.420 < 0.70; gnomAD AF 0.1500 > 0.0100 rarity cutoff; No ChEMBL bioactivity <= 1000.0 nM observed |

*(Note: ChEMBL queries with `limit=25` without potency sorting may return non-sub-micromolar assay rows or transient API timeouts for targets with large assay tables; filtering ChEMBL queries by `standard_value__lte=1000` is planned for `v0.2.0`.)*

## 6–12 Month Research Roadmap

- **`v0.2.0`:** Server-side ChEMBL potency filtering (`standard_value__lte=1000`) + GTEx v10 eQTL tissue-specificity cross-checks
- **`v0.3.0`:** Direct AlphaFold DB pLDDT confidence integration and multi-agent Google Antigravity `/run-workflow` benchmark harness
- **Dataset & Preprint:** Release expanded 500-locus multi-omic benchmark on Zenodo and submit manuscript to bioRxiv

## Citation

If you use this software in your research, please cite it using [CITATION.cff](CITATION.cff).

## License

MIT License — see [LICENSE](LICENSE).
