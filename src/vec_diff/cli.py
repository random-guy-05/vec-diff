from __future__ import annotations

import argparse
import sys
from pathlib import Path

from rich.console import Console
from rich.table import Table

from . import __version__
from .compare import compare_files
from .report import markdown, write_json

EQUIVALENT = {"BYTE_IDENTICAL", "SCORER_CONTENT_IDENTICAL", "ROW_ORDER_ONLY", "RIGID_FRAME_ONLY"}


def parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="vec-diff",
        description="Scorer-aware comparison of two Virtual Embryo Challenge .h5ad predictions.",
    )
    p.add_argument("a", type=Path)
    p.add_argument("b", type=Path)
    p.add_argument("--task", required=True, choices=["T1", "T2", "T3"])
    p.add_argument("--top-genes", type=int, default=15)
    p.add_argument(
        "--spatial-sample",
        type=int,
        default=512,
        help="Max cells per cloud for pairwise-distance shape sketch.",
    )
    p.add_argument(
        "--json",
        type=Path,
        dest="json_path",
        help="Write full machine-readable report.",
    )
    p.add_argument("--markdown", type=Path, help="Write Markdown report.")
    p.add_argument(
        "--fail-if-equivalent",
        action="store_true",
        help="Exit 3 if files are scorer-equivalent under a recognized invariance.",
    )
    p.add_argument(
        "--fail-if-material",
        action="store_true",
        help="Exit 4 if scorer-relevant content materially differs.",
    )
    p.add_argument("--version", action="version", version=f"vec-diff {__version__}")
    return p


def _print(report: dict) -> None:
    c = Console()
    c.print(f"[bold]{report['classification']}[/bold]")
    c.print(report["interpretation"])

    expression = report["expression"]
    table = Table(title="Expression")
    table.add_column("measure")
    table.add_column("value")
    for key in [
        "expression_multiset_overlap",
        "pseudobulk_pearson",
        "pseudobulk_mean_abs_delta",
        "rowwise_mean_abs_diff",
        "rowwise_rmse",
    ]:
        value = expression.get(key)
        table.add_row(key, "—" if value is None else f"{value:.6g}")
    c.print(table)

    if report.get("spatial"):
        spatial = report["spatial"]
        spatial_table = Table(title="Spatial")
        spatial_table.add_column("measure")
        spatial_table.add_column("value")
        for key in [
            "paired_multiset_overlap",
            "log_scale_b_over_a",
            "kabsch_relative_rmsd",
            "pairwise_distance_wasserstein_relative",
        ]:
            value = spatial.get(key)
            spatial_table.add_row(key, "—" if value is None else f"{value:.6g}")
        c.print(spatial_table)


def main(argv: list[str] | None = None) -> int:
    args = parser().parse_args(argv)
    for path in [args.a, args.b]:
        if not path.exists():
            raise SystemExit(f"file not found: {path}")

    report = compare_files(
        args.a,
        args.b,
        args.task,
        args.top_genes,
        args.spatial_sample,
    )
    _print(report)

    if args.json_path:
        args.json_path.parent.mkdir(parents=True, exist_ok=True)
        write_json(report, args.json_path)
    if args.markdown:
        args.markdown.parent.mkdir(parents=True, exist_ok=True)
        args.markdown.write_text(markdown(report), encoding="utf-8")
    if args.fail_if_equivalent and report["classification"] in EQUIVALENT:
        return 3
    if args.fail_if_material and report["classification"] == "MATERIAL":
        return 4
    return 0


if __name__ == "__main__":
    sys.exit(main())
