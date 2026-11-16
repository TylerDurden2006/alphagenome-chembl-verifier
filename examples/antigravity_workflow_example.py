"""Example integration showing how an agentic Google Antigravity workflow
invokes deterministic cross-verification on candidate regulatory variants.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from ag_chembl_verifier.batch import format_markdown_table, load_benchmark_loci, run_batch_verification


def main() -> None:
    bench_path = Path(__file__).resolve().parents[1] / "data" / "benchmark_loci.json"
    loci = load_benchmark_loci(bench_path)
    reports = run_batch_verification(loci[:3], offline=False)
    print(format_markdown_table(reports))


if __name__ == "__main__":
    main()
