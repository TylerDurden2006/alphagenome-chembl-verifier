"""Batch verification runner and Markdown/JSON report generator."""

from __future__ import annotations

import json
from dataclasses import asdict
from pathlib import Path
from typing import Callable

from .clients import fetch_chembl_target_activities, fetch_uniprot_accession
from .verifier import BioactivityRecord, VariantEvidence, VerificationReport, verify_target_evidence

UniprotFetcher = Callable[[str], str | None]
ActivityFetcher = Callable[[str], list[BioactivityRecord]]


def load_benchmark_loci(path: str | Path) -> list[dict]:
    """Load a benchmark loci JSON file."""
    payload = json.loads(Path(path).read_text(encoding="utf-8"))
    return list(payload.get("loci", []))


def run_batch_verification(
    loci: list[dict],
    *,
    offline: bool = False,
    uniprot_fetcher: UniprotFetcher = fetch_uniprot_accession,
    activity_fetcher: ActivityFetcher = fetch_chembl_target_activities,
) -> list[VerificationReport]:
    """Verify a list of variant-gene dictionaries and return VerificationReports."""
    reports: list[VerificationReport] = []
    for entry in loci:
        gene = str(entry["gene_symbol"])
        uniprot_acc = entry.get("uniprot_accession")
        activities: list[BioactivityRecord] = []

        if not offline:
            if not uniprot_acc:
                uniprot_acc = uniprot_fetcher(gene)
            if uniprot_acc:
                activities = activity_fetcher(uniprot_acc)

        evidence = VariantEvidence(
            variant_id=str(entry["variant_id"]),
            gene_symbol=gene,
            alphagenome_impact_score=float(entry["alphagenome_impact_score"]),
            gnomad_af=entry.get("gnomad_af"),
            uniprot_accession=uniprot_acc,
        )
        reports.append(verify_target_evidence(evidence, activities))
    return reports


def format_markdown_table(reports: list[VerificationReport]) -> str:
    """Format verification reports as a GitHub-flavored Markdown table."""
    lines = [
        "| Variant | Gene | UniProt | Passed | Tier | Potent Ligands | Notes |",
        "|---|---|---|---|---|---:|---|",
    ]
    for r in reports:
        notes = "; ".join(r.reasons) if r.reasons else "All gates passed"
        lines.append(
            f"| `{r.variant_id}` | `{r.gene_symbol}` | `{r.uniprot_accession or '-'}` | "
            f"{'PASS' if r.passed else 'FAIL'} | `{r.tier}` | {r.potent_ligands} | {notes} |"
        )
    return "\n".join(lines)


def reports_to_json(reports: list[VerificationReport]) -> str:
    """Serialize verification reports to formatted JSON."""
    return json.dumps([asdict(r) for r in reports], indent=2)
