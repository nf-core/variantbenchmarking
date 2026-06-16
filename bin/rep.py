#!/usr/bin/env python3
# coding=utf-8
#
# Create simple reports from GA4GH benchmarking results
#
# Usage:
#
#   For usage instructions run with option --help
#
# Original author:
#   Peter Krusche <pkrusche@illumina.com>
#
# Python 3 port + CDN asset loader for nf-core/variantbenchmarking:
#   - All `print` statements updated to print()
#   - `urllib2` replaced with `urllib.request`
#   - `gzip.GzipFile` open mode corrected for text output
#   - FileSystemLoader replaced with CdnFallbackLoader so that the vendored
#     lib/d3 asset is no longer required in bin/.
#     D3 is fetched from CDN at render time and inlined into the HTML output.
#     Bijou CSS is vendored directly into report.css.
#     All other GA4GH assets (helpers.js, roc.js, report.css, and all
#     *.jinja2.html templates) are loaded from the filesystem via TEMPLATEDIR,
#     which is resolved relative to this script — keeping the only required
#     layout as:
#
#       bin/
#         rep.py                        <- this file
#         ga4gh_report_assets/
#           report/
#             metrics.py               <- GA4GH Python helper (unchanged)
#           html/
#             report.jinja2.html
#             roc_plot.jinja2.html
#             metrics_table.jinja2.html
#             helpers.js
#             roc.js
#             report.css               <- includes vendored bijou CSS
#
#     No lib/bijou/ or lib/d3/ directory is needed anywhere.

import sys
import os
import json
import gzip
import copy
import argparse
import urllib.request

import jinja2

# ---------------------------------------------------------------------------
# Path setup
# ---------------------------------------------------------------------------
SCRIPTDIR   = os.path.abspath(os.path.dirname(__file__))
ASSETDIR    = os.path.join(SCRIPTDIR, "ga4gh_report_assets")
TEMPLATEDIR = os.path.join(ASSETDIR, "html")
LIBDIR      = ASSETDIR  # report/metrics.py lives at ga4gh_report_assets/report/metrics.py

sys.path.insert(0, LIBDIR)

import report.metrics

# ---------------------------------------------------------------------------
# CDN asset map
#
# Keys match the {% include "..." %} paths used in report.jinja2.html.
# Values are the CDN URLs to fetch and inline instead.
#
# D3 v3.5.x is pinned — roc.js uses the v3 API (d3.svg.axis, d3.behavior.zoom, etc.).
# ---------------------------------------------------------------------------

CDN_ASSETS: dict[str, str] = {
    "lib/d3/d3.min.js": (
        "https://cdnjs.cloudflare.com/ajax/libs/d3/3.5.17/d3.min.js"
    ),
}


class CdnFallbackLoader(jinja2.BaseLoader):
    """Jinja2 loader that serves CDN assets inline and falls back to disk.

    Any template name listed in CDN_ASSETS is fetched from the corresponding
    URL at first use and cached in memory for the lifetime of this loader.
    All other template names are delegated to a standard FileSystemLoader
    pointing at TEMPLATEDIR.

    This means the Jinja2 template files themselves (*.jinja2.html, helpers.js,
    roc.js, report.css) still live on disk in src/html/, while the heavy
    third-party libs are pulled from CDN without needing to be vendored.
    """

    def __init__(self, template_dir: str) -> None:
        self._fs = jinja2.FileSystemLoader(searchpath=template_dir)
        self._cache: dict[str, str] = {}

    def get_source(
        self, environment: jinja2.Environment, template: str
    ) -> tuple[str, str | None, object]:
        if template in CDN_ASSETS:
            if template not in self._cache:
                url = CDN_ASSETS[template]
                print(f"[rep.py] Fetching {template} from CDN: {url}", file=sys.stderr)
                try:
                    with urllib.request.urlopen(url, timeout=30) as resp:
                        self._cache[template] = resp.read().decode("utf-8")
                except Exception as exc:
                    raise jinja2.TemplateNotFound(template) from exc
            # Return (source, filename, uptodate).
            # filename=None and uptodate=lambda:True disables mtime checks.
            return self._cache[template], None, lambda: True

        # Delegate everything else (templates, helpers.js, roc.js, report.css)
        # to the filesystem loader.
        return self._fs.get_source(environment, template)

    def list_templates(self) -> list[str]:
        return self._fs.list_templates()


# ---------------------------------------------------------------------------
# Metrics extraction (unchanged logic, Python 3 style)
# ---------------------------------------------------------------------------

def extract_metrics(metrics: list) -> dict:
    """Extract metrics and return data dicts suitable for the Jinja2 template.

    Separates summary values (tables) from ROC datapoints (curves).

    :param metrics: list of metric objects as returned by report.metrics
    :return: dict with keys snp_table, snp_roc, indel_table, indel_roc, qq_fields
    """
    data_subset = [
        d.to_dict()
        for d in metrics
        if d.type in ("SNP", "INDEL")
        and d.filter in ("ALL", "PASS")
        and d.genotype == "*"
    ]

    data_subset_snp   = [d for d in data_subset if d["type"] == "SNP"   and d["subtype"] == "*"]
    data_subset_indel = [d for d in data_subset if d["type"] == "INDEL"]

    data_subset_snp.sort(key=lambda x: x["subset"])

    def _indel_type_order(t: str) -> int:
        ordering = {
            "D1_5":    1, "I1_5":    2, "C1_5":    3,
            "D6_15":   4, "I6_15":   5, "C6_15":   6,
            "D16_PLUS":7, "I16_PLUS":8, "C16_PLUS":9,
        }
        return ordering.get(t, 1000)

    data_subset_indel.sort(key=lambda x: (x["subset"], _indel_type_order(x["subtype"])))

    data_subset_snp_roc   = [copy.copy(d) for d in data_subset_snp   if d["subtype"] == "*" and d["subset"] == "*"]
    data_subset_indel_roc = [copy.copy(d) for d in data_subset_indel if d["subtype"] == "*" and d["subset"] == "*"]

    # Remove ROC arrays from table data (they only go into the curve plots)
    for d in data_subset_snp:
        del d["roc"]
    for d in data_subset_indel:
        del d["roc"]

    qq_fields = list({x.qq_field for x in metrics})

    return {
        "snp_table":   json.dumps(json.dumps(data_subset_snp)),
        "snp_roc":     json.dumps(json.dumps(data_subset_snp_roc)),
        "indel_table": json.dumps(json.dumps(data_subset_indel)),
        "indel_roc":   json.dumps(json.dumps(data_subset_indel_roc)),
        "qq_fields":   qq_fields,
    }


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(description="Create a variant calling report.")

    parser.add_argument(
        "input",
        nargs="*",
        help=(
            "Input file in GA4GH metrics CSV format. "
            "To label multiple results use: "
            "rep.py gatk-3_vcfeval-giab:gatk3.roc.all.csv.gz -o test.html "
            "(label format: <method>_<comparison>:<file>)"
        ),
    )
    parser.add_argument(
        "-o", "--output",
        required=True,
        help="Output HTML file name (e.g. 'report' or 'report.html').",
    )
    parser.add_argument(
        "--roc-max-datapoints",
        dest="roc_datapoints",
        type=int,
        default=1000,
        help="Maximum number of data points per ROC curve.",
    )
    parser.add_argument(
        "--roc-resolution",
        dest="roc_diff",
        type=float,
        default=0.005,
        help="Minimum precision/recall difference between ROC datapoints.",
    )
    parser.add_argument(
        "--min-recall",
        dest="min_recall",
        type=float,
        default=0.2,
        help="Minimum recall for ROC curves.",
    )
    parser.add_argument(
        "--min-precision",
        dest="min_precision",
        type=float,
        default=0.0,
        help="Minimum precision for ROC curves.",
    )

    args = parser.parse_args()

    # Determine output path / handle gzip output
    output = args.output
    if output.endswith(".gz"):
        output_handle = gzip.open(output, "wt", encoding="utf-8")
    else:
        if not output.endswith(".html"):
            output += ".html"
        output_handle = open(output, "w", encoding="utf-8")  # noqa: WPS515

    # 1. Read input files
    metrics = []
    for i in args.input:
        parts = i.split(":")
        method_label  = "default"
        cmethod_label = "default"

        if len(parts) <= 1:
            rfiles = [parts[0]]
        else:
            rfiles = parts[1:]
            labels = parts[0].split("_")
            if len(labels) > 0:
                method_label = labels[0]
            if len(labels) > 1:
                cmethod_label = labels[1]

        print(
            f"[rep.py] Reading {rfiles} as method={method_label!r} / comparison={cmethod_label!r}",
            file=sys.stderr,
        )

        row_metrics = report.metrics.read_qfy_csv(
            rfiles,
            method=method_label,
            cmethod=cmethod_label,
            roc_metrics=["METRIC.Precision", "METRIC.Recall"],
            roc_diff=args.roc_diff,
            max_data_points=args.roc_datapoints,
            minmax={
                "METRIC.Precision": {"min": args.min_precision},
                "METRIC.Recall":    {"min": args.min_recall},
            },
        )
        metrics.extend(row_metrics)

    if not metrics:
        raise RuntimeError("No inputs specified or no metrics could be read.")

    # 2. Subset data
    template_vars = extract_metrics(metrics)

    # 3. Render template via CdnFallbackLoader
    loader = CdnFallbackLoader(TEMPLATEDIR)
    env    = jinja2.Environment(loader=loader)
    template = env.get_template("report.jinja2.html")

    with output_handle:
        template.stream(**template_vars).dump(output_handle)

    print(f"[rep.py] Report written to: {output if not output.endswith('.gz') else args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
