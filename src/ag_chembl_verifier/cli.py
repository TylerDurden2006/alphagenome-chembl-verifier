"""Command-line interface for alphagenome-chembl-verifier."""

from __future__ import annotations

import argparse
import json
from dataclasses import asdict

from .batch import (
    format_markdown_table,
    load_benchmark_loci,
    reports_to_json,
    run_batch_verification,
)
from .clients import fetch_chembl_target_activities, fetch_uniprot_accession
from .verifier import VariantEvidence, verify_target_evidence


def build_parser() -> argparse.ArgumentParser:
    """Create the CLI argument parser."""
    parser = argparse.ArgumentParser(
        prog="ag-verify",
        description="Cross-verify regulatory variants against UniProt and ChEMBL.",
    )
    parser.add_argument("--batch", help="Path to benchmark JSON file for batch verification")
    parser.add_argument(
        "--format",
        choices=["json", "markdown"],
        default="json",
        help="Output format (default: json)",
    )
    parser.add_argument("--variant", help="Variant ID (e.g. chr7:140753336:A>T)")
    parser.add_argument("--gene", help="Target HGNC gene symbol (e.g. BRAF)")
    parser.add_argument(
        "--impact",
        type=float,
        help="Normalized AlphaGenome regulatory impact score in [0.0, 1.0]",
    )
    parser.add_argument(
        "--gnomad-af",
        type=float,
        default=None,
        help="Optional gnomAD allele frequency in [0.0, 1.0]",
    )
    parser.add_argument(
        "--offline",
        action="store_true",
        help="Skip live UniProt/ChEMBL HTTPS calls",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    """Run CLI verification and print JSON or Markdown report."""
    parser = build_parser()
    args = parser.parse_args(argv)

    if args.batch:
        loci = load_benchmark_loci(args.batch)
        reports = run_batch_verification(loci, offline=args.offline)
        print(format_markdown_table(reports) if args.format == "markdown" else reports_to_json(reports))
        return 0

    if not (args.variant and args.gene and args.impact is not None):
        parser.error("Either --batch or (--variant, --gene, --impact) must be provided.")

    uniprot_acc = None if args.offline else fetch_uniprot_accession(args.gene)
    activities = (
        []
        if (args.offline or not uniprot_acc)
        else fetch_chembl_target_activities(uniprot_acc)
    )
    evidence = VariantEvidence(
        variant_id=args.variant,
        gene_symbol=args.gene,
        alphagenome_impact_score=args.impact,
        gnomad_af=args.gnomad_af,
        uniprot_accession=uniprot_acc,
    )
    report = verify_target_evidence(evidence, activities)
    print(json.dumps(asdict(report), indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
