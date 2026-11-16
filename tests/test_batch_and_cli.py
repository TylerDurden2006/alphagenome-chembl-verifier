"""Unit tests for batch verification, benchmark loading, and CLI."""

import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from ag_chembl_verifier.batch import (  # noqa: E402
    format_markdown_table,
    load_benchmark_loci,
    reports_to_json,
    run_batch_verification,
)
from ag_chembl_verifier.cli import main as cli_main  # noqa: E402
from ag_chembl_verifier.verifier import BioactivityRecord  # noqa: E402


class BatchAndCliTests(unittest.TestCase):
    """Test suite for batch pipeline and CLI."""

    def test_benchmark_loci_file_loads_and_runs_with_mocked_apis(self) -> None:
        bench_path = Path(__file__).resolve().parents[1] / "data" / "benchmark_loci.json"
        loci = load_benchmark_loci(bench_path)
        self.assertGreaterEqual(len(loci), 10)

        def fake_uniprot(gene: str) -> str:
            return f"P_{gene}"

        def fake_ensembl(gene: str) -> str:
            return f"ENSG_{gene}"

        def fake_tractability(ensg: str) -> tuple[str, ...]:
            return ("SM:Clinical Precedence",)

        def fake_chembl(acc: str) -> list[BioactivityRecord]:
            return [
                BioactivityRecord("CHEMBL1", "MOL1", "IC50", 25.0),
                BioactivityRecord("CHEMBL1", "MOL2", "IC50", 80.0),
                BioactivityRecord("CHEMBL1", "MOL3", "Kd", 150.0),
            ]

        reports = run_batch_verification(
            loci,
            offline=False,
            uniprot_fetcher=fake_uniprot,
            activity_fetcher=fake_chembl,
            ensembl_fetcher=fake_ensembl,
            tractability_fetcher=fake_tractability,
        )
        self.assertEqual(len(reports), len(loci))
        for entry, report in zip(loci, reports):
            self.assertEqual(report.passed, entry["expected_pass"], msg=entry["gene_symbol"])
            self.assertTrue(report.ensembl_gene_id.startswith("ENSG_"))

        md = format_markdown_table(reports)
        self.assertIn("TIER_1_HIGH_CONFIDENCE", md)
        self.assertIn("UNVERIFIED", md)
        self.assertIn("SM:Clinical Precedence", md)

        js_out = reports_to_json(reports)
        self.assertIn("ensembl_gene_id", js_out)

    def test_cli_offline_single_variant(self) -> None:
        code = cli_main(
            [
                "--variant",
                "chr7:140753336:A>T",
                "--gene",
                "BRAF",
                "--impact",
                "0.92",
                "--gnomad-af",
                "0.0001",
                "--offline",
            ]
        )
        self.assertEqual(code, 0)

    def test_cli_offline_batch_runs(self) -> None:
        bench_path = Path(__file__).resolve().parents[1] / "data" / "benchmark_loci.json"
        code = cli_main(["--batch", str(bench_path), "--offline", "--format", "markdown"])
        self.assertEqual(code, 0)


if __name__ == "__main__":
    unittest.main()
