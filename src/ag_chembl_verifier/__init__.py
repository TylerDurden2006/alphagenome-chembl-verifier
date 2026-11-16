"""alphagenome-chembl-verifier: Cross-verification for regulatory variants and chemogenomic targets."""

from .verifier import (
    BioactivityRecord,
    VariantEvidence,
    VerificationReport,
    verify_target_evidence,
)

__version__ = "0.1.0"
__all__ = [
    "BioactivityRecord",
    "VariantEvidence",
    "VerificationReport",
    "verify_target_evidence",
]
