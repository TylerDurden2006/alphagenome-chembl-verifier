"""Lightweight standard-library REST clients for UniProtKB and ChEMBL."""

from __future__ import annotations

import json
import urllib.parse
import urllib.request
from typing import Any

from .verifier import BioactivityRecord

DEFAULT_TIMEOUT_SEC = 15


def _get_json(url: str, timeout: int = DEFAULT_TIMEOUT_SEC) -> dict[str, Any]:
    """Fetch JSON from an HTTPS endpoint using Python's standard library."""
    req = urllib.request.Request(
        url,
        headers={
            "Accept": "application/json",
            "User-Agent": "alphagenome-chembl-verifier/0.1.0",
        },
    )
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        return json.loads(resp.read().decode("utf-8"))


def fetch_uniprot_accession(gene_symbol: str, organism_id: int = 9606) -> str | None:
    """Resolve a human gene symbol to its primary reviewed UniProtKB accession."""
    query = urllib.parse.quote(
        f"gene_exact:{gene_symbol} AND organism_id:{organism_id} AND reviewed:true"
    )
    url = (
        f"https://rest.uniprot.org/uniprotkb/search"
        f"?query={query}&fields=accession&format=json&size=1"
    )
    data = _get_json(url)
    results = data.get("results", [])
    if not results:
        return None
    return results[0].get("primaryAccession")


def fetch_chembl_target_activities(
    uniprot_accession: str,
    limit: int = 25,
) -> list[BioactivityRecord]:
    """Fetch standardized IC50/Ki/Kd bioactivity records for a UniProt accession."""
    target_url = (
        "https://www.ebi.ac.uk/chembl/api/data/target.json"
        f"?target_components__accession={urllib.parse.quote(uniprot_accession)}&limit=1"
    )
    target_data = _get_json(target_url)
    targets = target_data.get("targets", [])
    if not targets:
        return []

    target_id = targets[0]["target_chembl_id"]
    act_url = (
        "https://www.ebi.ac.uk/chembl/api/data/activity.json"
        f"?target_chembl_id={urllib.parse.quote(target_id)}"
        f"&standard_units=nM&limit={int(limit)}"
    )
    act_data = _get_json(act_url)
    records: list[BioactivityRecord] = []
    for row in act_data.get("activities", []):
        raw_val = row.get("standard_value")
        if raw_val is None:
            continue
        try:
            val_nm = float(raw_val)
        except (TypeError, ValueError):
            continue
        records.append(
            BioactivityRecord(
                chembl_target_id=target_id,
                molecule_chembl_id=str(row.get("molecule_chembl_id", "")),
                standard_type=str(row.get("standard_type", "UNKNOWN")),
                standard_value_nm=val_nm,
                assay_type=str(row.get("assay_type", "B")),
            )
        )
    return records
