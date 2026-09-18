#!/usr/bin/env python3

# Copyright 2024 - GHGA
# Author: Kuebra Narci - @kubranarci

import argparse
import sys
import pandas as pd
import plotly.express as px

SEABORN_CB_PALETTE = [
    "#0173b2", "#d55e00", "#029e73", "#de8f05", "#cc78bc",
    "#ca9161", "#fbafe4", "#949494", "#ece133", "#56b4e9"
]

TERM_MAPPING = {
    "GT": "Genotype Matching",
    "BASEPAIR": "Sequence Matching",
    "Snv": "SNVs",
    "SNP": "SNVs",
    "INDEL": "Indels",
    "Insertion": "Insertions",
    "Deletion": "Deletions",
    "JointIndel": "Joint Indels",
    "PASS": "PASS Filter Only"
}

def apply_global_styling(fig, x_title, y_title, height=None):
    """Applies bold text, distinct subplot borders, and removes internal titles."""
    layout_dict = dict(
        template="plotly_white",
        font=dict(size=18, color="black"),
        legend=dict(
            title_text="<b>Tool</b>",
            font=dict(size=18),
            title_font=dict(size=20)
        ),
        margin=dict(l=80, r=40, t=50, b=80),
        autosize=True
    )
    if height:
        layout_dict["height"] = height

    fig.update_layout(**layout_dict)

    fig.update_xaxes(
        title_text=f"<b>{x_title}</b>" if x_title else None,
        tickfont=dict(size=18),
        showline=True, linewidth=1.5, linecolor='black', mirror=True
    )
    fig.update_yaxes(
        title_text=f"<b>{y_title}</b>" if y_title else None,
        tickfont=dict(size=18),
        showline=True, linewidth=1.5, linecolor='black', mirror=True
    )

    fig.update_xaxes(tickformat="<b>%{text}</b>")

    return fig

def generate_plots(df, benchmark, clean_prefix, title_suffix, show_labels, color_map):
    title_tp = "Variant Comparison Metrics" if title_suffix == "Overall" else f"Variant Comparison Metrics - {title_suffix}"
    title_f1 = "F1 Score" if title_suffix == "Overall" else f"F1 Score - {title_suffix}"
    title_pr = "Precision vs Recall" if title_suffix == "Overall" else f"Precision vs Recall - {title_suffix}"

    id_vars = ["Tool"]
    melted = df.melt(id_vars=id_vars)

    tp_vars = ["TP_comp", "FP", "FN"]
    tp_data = melted[melted["variable"].isin(tp_vars)].copy()
    tp_data["value"] = pd.to_numeric(tp_data["value"])
    tp_data["variable"] = pd.Categorical(tp_data["variable"], categories=tp_vars, ordered=True)
    tp_data = tp_data.sort_values(["variable", "Tool"])

    metric_data = melted[melted["variable"] == "F1"].copy()
    metric_data["value"] = pd.to_numeric(metric_data["value"])

    plot_dict = {
        "tp": None, "tp_title": title_tp,
        "f1": None, "f1_title": title_f1,
        "pr": None, "pr_title": title_pr
    }

    # 1. Visualize TP_comp, FP, FN
    if not tp_data.empty:
        fig_tp = px.line(
            tp_data, x="Tool", y="value", color="Tool", facet_col="variable",
            markers=True, color_discrete_map=color_map,
            facet_col_spacing=0.08
        )
        fig_tp = apply_global_styling(fig_tp, "Tool", None, height=500)
        fig_tp.update_yaxes(matches=None, automargin=True, showticklabels=True, tickformat="~s")

        facet_annotations = [a for a in fig_tp.layout.annotations if a.text and "=" in a.text]
        yaxes = list(fig_tp.select_yaxes())

        for yaxis, ann in zip(yaxes, facet_annotations):
            var_name = ann.text.split("=")[-1].replace("TP_comp", "TP")
            yaxis.title.text = f"<b>{var_name} count</b>"
            ann.text = ""

        if show_labels:
            def format_label(val):
                if val >= 1000:
                    return f"<b>{val/1000:g}k</b>"
                return f"<b>{val}</b>"

            fig_tp.update_traces(
                textposition="top center",
                text=[format_label(v) for v in tp_data["value"]],
                mode="lines+markers+text",
                textfont=dict(size=16)
            )
        else:
            fig_tp.update_traces(marker=dict(size=10), line=dict(width=4))

        fig_tp.update_layout(hovermode="x unified")
        fig_tp.update_xaxes(tickvals=tp_data["Tool"].unique(), ticktext=[f"<b>{t}</b>" for t in tp_data["Tool"].unique()])

        plot_dict["tp"] = fig_tp

    # 2. Visualize F1
    if not metric_data.empty:
        fig_f1 = px.scatter(
            metric_data, x="Tool", y="value", color="Tool",
            color_discrete_map=color_map
        )
        fig_f1 = apply_global_styling(fig_f1, "Tool", "F1 Score", height=450)
        fig_f1.update_yaxes(range=[0, 1.05])

        if show_labels:
            fig_f1.update_traces(
                textposition="top center",
                text=["<b>" + str(round(v, 3)) + "</b>" for v in metric_data["value"]],
                mode="markers+text",
                textfont=dict(size=16)
            )
        else:
            fig_f1.update_traces(marker=dict(size=12))

        fig_f1.update_layout(hovermode="x unified")
        fig_f1.update_xaxes(tickvals=metric_data["Tool"].unique(), ticktext=[f"<b>{t}</b>" for t in metric_data["Tool"].unique()])
        plot_dict["f1"] = fig_f1

    # 3. Visualize Precision vs Recall
    if "Recall" in df.columns and "Precision" in df.columns:
        fig_pr = px.scatter(
            df, x="Recall", y="Precision", color="Tool",
            color_discrete_map=color_map
        )
        fig_pr = apply_global_styling(fig_pr, "Recall", "Precision", height=450)
        fig_pr.update_xaxes(range=[-0.02, 1.02])
        fig_pr.update_yaxes(range=[-0.02, 1.02])
        fig_pr.update_traces(marker=dict(size=12))
        fig_pr.update_layout(hovermode="closest")
        plot_dict["pr"] = fig_pr

    return plot_dict

def main():
    parser = argparse.ArgumentParser(description="Generate interactive metrics plots")
    parser.add_argument("summary_file", help="Input summary CSV file")
    parser.add_argument("benchmark", help="Benchmark name")
    parser.add_argument("--output", default="metrics_report.html", help="Output HTML file")
    parser.add_argument("--labels", action="store_true", help="Show point labels on plots")
    args = parser.parse_args()

    try:
        df = pd.read_csv(args.summary_file)
    except Exception as e:
        print(f"Error reading {args.summary_file}: {e}")
        sys.exit(1)

    if args.benchmark == "rtgtools" and "Threshold" in df.columns:
        mask = df["Threshold"].astype(str).str.strip().str.lower() == "none"
        if mask.any():
            df = df[mask]

    unique_tools = sorted(df["Tool"].dropna().unique())
    global_color_map = {tool: SEABORN_CB_PALETTE[i % len(SEABORN_CB_PALETTE)] for i, tool in enumerate(unique_tools)}

    possible_groupings = ["Type", "Filter", "StatsType", "Comparison", "Region"]
    split_cols = [col for col in possible_groupings if col in df.columns]

    all_groups = []

    if split_cols:
        grouped = df.groupby(split_cols)
        for name, group in grouped:
            parts = list(name) if isinstance(name, tuple) else [name]
            parts_str = [str(p) for p in parts]
            valid_parts = [p for p in parts_str if p.lower() not in ["all", "none"]]

            if valid_parts:
                clean_prefix = "_".join(valid_parts)
                pretty_parts = [TERM_MAPPING.get(p, p) for p in valid_parts]
                title_suffix = " | ".join(pretty_parts)
            else:
                clean_prefix = "overall"
                title_suffix = "Overall"

            all_groups.append(generate_plots(group, args.benchmark, clean_prefix, title_suffix, args.labels, global_color_map))
    else:
        all_groups.append(generate_plots(df, args.benchmark, "overall", "Overall", args.labels, global_color_map))

    if not all_groups:
        print("No valid data found to generate plots.")
        sys.exit(0)

    html_blocks = []
    plotly_js_included = False
    plot_config = {'responsive': True}

    for group in all_groups:
        group_html = '<div class="group-container">'

        if group["tp"]:
            tp_html = group["tp"].to_html(full_html=False, include_plotlyjs='cdn' if not plotly_js_included else False, config=plot_config)
            plotly_js_included = True
            group_html += f'<div class="plot-full"><div class="html-title">{group["tp_title"]}</div>{tp_html}</div>'

        if group["f1"] or group["pr"]:
            group_html += '<div class="flex-row">'
            if group["f1"]:
                f1_html = group["f1"].to_html(full_html=False, include_plotlyjs='cdn' if not plotly_js_included else False, config=plot_config)
                plotly_js_included = True
                group_html += f'<div class="plot-half"><div class="html-title">{group["f1_title"]}</div>{f1_html}</div>'
            if group["pr"]:
                pr_html = group["pr"].to_html(full_html=False, include_plotlyjs='cdn' if not plotly_js_included else False, config=plot_config)
                plotly_js_included = True
                group_html += f'<div class="plot-half"><div class="html-title">{group["pr_title"]}</div>{pr_html}</div>'
            group_html += '</div>'

        group_html += '</div>'
        html_blocks.append(group_html)

    html_template = f"""
    <!DOCTYPE html>
    <html>
    <head>
        <title>Metrics Report</title>
        <style>
            body {{ font-family: -apple-system, sans-serif; background-color: #f8f9fa; margin: 20px; }}
            .group-container {{ margin-bottom: 50px; border-bottom: 2px solid #dee2e6; padding-bottom: 20px; }}
            .plot-full {{ width: 100%; margin-bottom: 20px; background: white; padding: 20px; border-radius: 8px; box-shadow: 0 2px 4px rgba(0,0,0,0.05); box-sizing: border-box; }}
            .flex-row {{ display: flex; gap: 20px; flex-wrap: wrap; width: 100%; box-sizing: border-box; }}
            .plot-half {{ flex: 1; min-width: 480px; background: white; padding: 20px; border-radius: 8px; box-shadow: 0 2px 4px rgba(0,0,0,0.05); box-sizing: border-box; }}
            .plot-half .plotly-graph-div {{ width: 100% !important; }}
            .plot-full .plotly-graph-div {{ width: 100% !important; }}
            .html-title {{ font-size: 22px; font-weight: bold; color: #2c3e50; margin-bottom: 15px; text-align: center; word-wrap: break-word; }}
        </style>
    </head>
    <body>
        {''.join(html_blocks)}
    </body>
    </html>
    """

    with open(args.output, "w") as f:
        f.write(html_template)

if __name__ == "__main__":
    main()
