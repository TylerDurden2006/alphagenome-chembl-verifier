# alphagenome-chembl-verifier

Deterministic cross-verification pipeline linking non-coding regulatory variant impact predictions (AlphaGenome) with population allele frequency constraints (gnomAD), reviewed protein annotations (UniProtKB), and target bioactivity evidence (ChEMBL).

## Why This Exists

Multi-agent scientific discovery pipelines can generate plausible variant-to-target hypotheses at scale, but hallucinated gene symbols, common-variant false positives, or undruggable targets quickly degrade downstream experimental ROI. `alphagenome-chembl-verifier` enforces **deterministic, auditable evidence gates** using live public REST APIs (`rest.uniprot.org` and `ebi.ac.uk/chembl/api/data`).

## Architecture

```
Variant + AlphaGenome Impact Score + gnomAD AF
        │
        ▼
[1] Schema & Rarity Gate (impact >= 0.70, gnomAD AF <= 0.01)
        │
        ▼
[2] UniProtKB Reviewed Human Accession Resolution (rest.uniprot.org)
        │
        ▼
[3] ChEMBL Target & Sub-Micromolar Bioactivity Audit (<= 1000 nM)
        │
        ▼
VerificationReport (TIER_1_HIGH_CONFIDENCE | TIER_2_VALIDATED | UNVERIFIED)
```

## Quick Start

Zero external runtime dependencies (Python 3.10+ standard library).

```bash
# Run unit tests
python -m unittest discover -s tests -v

# Single-variant live verification
PYTHONPATH=src python -m ag_chembl_verifier.cli --variant chr7:140753336:A>T --gene BRAF --impact 0.92 --gnomad-af 0.0001

# Batch benchmark verification (live UniProtKB + ChEMBL)
PYTHONPATH=src python -m ag_chembl_verifier.cli --batch data/benchmark_loci.json --format markdown
```

## Live Benchmark Output (`data/benchmark_loci.json`)

Executed against live `rest.uniprot.org` and `ebi.ac.uk/chembl/api/data` endpoints (`limit=25` activities per target):

| Variant | Gene | UniProt | Passed | Tier | Potent Ligands | Notes |
|---|---|---|---|---|---:|---|
| `chr7:140753336:A>T` | `BRAF` | `P15056` | PASS | `TIER_1_HIGH_CONFIDENCE` | 23 | All gates passed |
| `chr7:55174772:G>A` | `EGFR` | `P00533` | PASS | `TIER_1_HIGH_CONFIDENCE` | 6 | All gates passed |
| `chr1:55039974:G>T` | `PCSK9` | `Q8NBP7` | PASS | `TIER_1_HIGH_CONFIDENCE` | 3 | All gates passed |
| `chr18:63318364:C>T` | `BCL2` | `P10415` | FAIL | `UNVERIFIED` | 0 | No ChEMBL bioactivity <= 1000.0 nM observed |
| `chr5:1295113:G>A` | `TERT` | `O14746` | PASS | `TIER_2_VALIDATED` | 1 | All gates passed |
| `chr17:7674220:C>T` | `TP53` | `P04637` | FAIL | `UNVERIFIED` | 7 | gnomAD AF 0.1200 > 0.0100 rarity cutoff |
| `chr19:44908684:T>C` | `APOE` | `P02649` | FAIL | `UNVERIFIED` | 0 | Impact score 0.420 < 0.70; gnomAD AF 0.1500 > 0.0100 rarity cutoff; No ChEMBL bioactivity <= 1000.0 nM observed |

## Citation

If you use this software in your research, please cite it using `CITATION.cff`.

## License

MIT License — see `LICENSE`.
