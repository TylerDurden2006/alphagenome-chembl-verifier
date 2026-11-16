"""Unit tests for deterministic variant-to-target verification."""

import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from ag_chembl_verifier.verifier import (  # noqa: E402
    BioactivityRecord,
    VariantEvidence,
    verify_target_evidence,
)


class VerifierTests(unittest.TestCase):
    """Test suite for verify_target_evidence."""

    def test_high_confidence_pass(self) -> None:
        variant = VariantEvidence(
            variant_id="chr7:140753336:A>T",
            gene_symbol="BRAF",
            alphagenome_impact_score=0.92,
            gnomad_af=0.0001,
            uniprot_accession="P15056",
        )
        activities = [
            BioactivityRecord("CHEMBL5145", "CHEMBL1229517", "IC50", 31.0),
            BioactivityRecord("CHEMBL5145", "CHEMBL1336", "IC50", 45.0),
            BioactivityRecord("CHEMBL5145", "CHEMBL502835", "Kd", 12.5),
        ]
        report = verify_target_evidence(variant, activities)
        self.assertTrue(report.passed)
        self.assertEqual(report.tier, "TIER_1_HIGH_CONFIDENCE")
        self.assertEqual(report.potent_ligands, 3)

    def test_rejects_common_variant(self) -> None:
        variant = VariantEvidence(
            variant_id="chr1:12345:C>G",
            gene_symbol="TP53",
            alphagenome_impact_score=0.85,
            gnomad_af=0.12,
            uniprot_accession="P04637",
        )
        activities = [BioactivityRecord("CHEMBL4096", "CHEMBL1", "IC50", 50.0)]
        report = verify_target_evidence(variant, activities)
        self.assertFalse(report.passed)
        self.assertEqual(report.tier, "UNVERIFIED")
        self.assertTrue(any("gnomAD AF" in r for r in report.reasons))

    def test_rejects_out_of_range_score(self) -> None:
        variant = VariantEvidence(
            variant_id="chr1:12345:C>G",
            gene_symbol="EGFR",
            alphagenome_impact_score=1.5,
            uniprot_accession="P00533",
        )
        report = verify_target_evidence(variant, [])
        self.assertFalse(report.passed)
        self.assertEqual(report.tier, "INVALID_INPUT")


if __name__ == "__main__":
    unittest.main()
