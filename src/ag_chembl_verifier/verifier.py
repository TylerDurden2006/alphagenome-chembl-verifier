"""Deterministic cross-verification rules for regulatory variants and chemogenomic targets."""

from __future__ import annotations

from dataclasses import dataclass, field


@dataclass(frozen=True)
class VariantEvidence:
    """Normalized regulatory variant evidence record."""

    variant_id: str
    gene_symbol: str
    alphagenome_impact_score: float
    gnomad_af: float | None = None
    uniprot_accession: str | None = None
    ensembl_gene_id: str | None = None
    tractability_labels: tuple[str, ...] = field(default_factory=tuple)

    def validate(self) -> list[str]:
        """Return validation errors for malformed input fields."""
        errors: list[str] = []
        if not self.variant_id.strip():
            errors.append("variant_id must be non-empty")
        if not self.gene_symbol.strip():
            errors.append("gene_symbol must be non-empty")
        if not (0.0 <= self.alphagenome_impact_score <= 1.0):
            errors.append("alphagenome_impact_score must be in [0.0, 1.0]")
        if self.gnomad_af is not None and not (0.0 <= self.gnomad_af <= 1.0):
            errors.append("gnomad_af must be in [0.0, 1.0]")
        return errors


@dataclass(frozen=True)
class BioactivityRecord:
    """Normalized ChEMBL bioactivity measurement for a target."""

    chembl_target_id: str
    molecule_chembl_id: str
    standard_type: str
    standard_value_nm: float
    assay_type: str = "B"


@dataclass(frozen=True)
class VerificationReport:
    """Audit output summarizing variant-to-target tractability checks."""

    variant_id: str
    gene_symbol: str
    ensembl_gene_id: str | None
    uniprot_accession: str | None
    passed: bool
    tier: str
    potent_ligands: int
    tractability_labels: tuple[str, ...] = field(default_factory=tuple)
    reasons: tuple[str, ...] = field(default_factory=tuple)


def verify_target_evidence(
    variant: VariantEvidence,
    activities: list[BioactivityRecord],
    *,
    impact_threshold: float = 0.70,
    max_gnomad_af: float = 0.01,
    potency_nm_cutoff: float = 1000.0,
) -> VerificationReport:
    """Verify whether a regulatory variant maps to a tractably druggable target.

    Applies deterministic filters across variant impact score, population allele
    frequency (gnomAD), UniProt resolution, and ChEMBL sub-micromolar bioactivity.
    """
    errors = variant.validate()
    if errors:
        return VerificationReport(
            variant_id=variant.variant_id,
            gene_symbol=variant.gene_symbol,
            ensembl_gene_id=variant.ensembl_gene_id,
            uniprot_accession=variant.uniprot_accession,
            passed=False,
            tier="INVALID_INPUT",
            potent_ligands=0,
            tractability_labels=variant.tractability_labels,
            reasons=tuple(errors),
        )

    reasons: list[str] = []

    if variant.alphagenome_impact_score < impact_threshold:
        reasons.append(
            f"Impact score {variant.alphagenome_impact_score:.3f} < {impact_threshold:.2f}"
        )

    if variant.gnomad_af is not None and variant.gnomad_af > max_gnomad_af:
        reasons.append(
            f"gnomAD AF {variant.gnomad_af:.4f} > {max_gnomad_af:.4f} rarity cutoff"
        )

    if not variant.uniprot_accession:
        reasons.append("Missing reviewed UniProtKB accession mapping")

    potent = [
        act
        for act in activities
        if act.standard_value_nm > 0 and act.standard_value_nm <= potency_nm_cutoff
    ]
    if not potent:
        reasons.append(
            f"No ChEMBL bioactivity <= {potency_nm_cutoff:.1f} nM observed"
        )

    passed = len(reasons) == 0
    if passed and len(potent) >= 3:
        tier = "TIER_1_HIGH_CONFIDENCE"
    elif passed:
        tier = "TIER_2_VALIDATED"
    else:
        tier = "UNVERIFIED"

    return VerificationReport(
        variant_id=variant.variant_id,
        gene_symbol=variant.gene_symbol,
        ensembl_gene_id=variant.ensembl_gene_id,
        uniprot_accession=variant.uniprot_accession,
        passed=passed,
        tier=tier,
        potent_ligands=len(potent),
        tractability_labels=variant.tractability_labels,
        reasons=tuple(reasons),
    )
