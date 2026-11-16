"""Lightweight standard-library REST/GraphQL clients for Ensembl, UniProtKB, ChEMBL, and OpenTargets."""

from __future__ import annotations

import json
import time
import urllib.error
import urllib.parse
import urllib.request
from typing import Any

from .verifier import BioactivityRecord

DEFAULT_TIMEOUT_SEC = 30


def _get_json(url: str, timeout: int = DEFAULT_TIMEOUT_SEC, retries: int = 2) -> Any:
    """Fetch JSON from an HTTPS GET endpoint with bounded retries."""
    req = urllib.request.Request(
        url,
        headers={
            "Accept": "application/json",
            "User-Agent": "alphagenome-chembl-verifier/0.1.0",
        },
    )
    for attempt in range(retries + 1):
        try:
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                return json.loads(resp.read().decode("utf-8"))
        except (urllib.error.URLError, TimeoutError):
            if attempt == retries:
                return None
            time.sleep(0.5 * (attempt + 1))
    return None


def _post_json(
    url: str,
    payload: dict[str, Any],
    timeout: int = DEFAULT_TIMEOUT_SEC,
    retries: int = 2,
) -> Any:
    """Send a JSON POST request with bounded retries."""
    req = urllib.request.Request(
        url,
        data=json.dumps(payload).encode("utf-8"),
        headers={
            "Content-Type": "application/json",
            "Accept": "application/json",
            "User-Agent": "alphagenome-chembl-verifier/0.1.0",
        },
    )
    for attempt in range(retries + 1):
        try:
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                return json.loads(resp.read().decode("utf-8"))
        except (urllib.error.URLError, TimeoutError):
            if attempt == retries:
                return None
            time.sleep(0.5 * (attempt + 1))
    return None


def fetch_ensembl_gene_id(gene_symbol: str, species: str = "homo_sapiens") -> str | None:
    """Resolve an HGNC gene symbol to its primary Ensembl Gene ID (ENSG...)."""
    url = (
        f"https://rest.ensembl.org/xrefs/symbol/{urllib.parse.quote(species)}/"
        f"{urllib.parse.quote(gene_symbol)}?content-type=application/json"
    )
    data = _get_json(url)
    if not data:
        return None
    for item in data:
        if item.get("type") == "gene" and str(item.get("id", "")).startswith("ENSG"):
            return str(item["id"])
    return None


def fetch_opentargets_tractability(ensembl_gene_id: str) -> tuple[str, ...]:
    """Fetch positive clinical/small-molecule tractability flags from OpenTargets GraphQL."""
    query = {
        "query": (
            '{ target(ensemblId: "'
            + ensembl_gene_id
            + '") { tractability { label modality value } } }'
        )
    }
    data = _post_json("https://api.platform.opentargets.org/api/v4/graphql", query)
    if not data:
        return ()
    target = (data.get("data") or {}).get("target") or {}
    labels: list[str] = []
    for entry in target.get("tractability") or []:
        if entry.get("value") is True:
            labels.append(f"{entry.get('modality')}:{entry.get('label')}")
    return tuple(labels[:5])


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
    if not data:
        return None
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
    if not target_data:
        return []
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
    if not act_data:
        return []
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
