#!/usr/bin/env python3
# Usage:
#   python plot_job_durations_over_time.py [job_durations.csv] [out.png]
#   Defaults: csv='job_durations.csv', out='job_durations_over_time.png'

import sys, os
import math
import pandas as pd
import matplotlib.pyplot as plt

def load_data(csv_path: str) -> pd.DataFrame:
    df = pd.read_csv(csv_path)
    if "timestamp_iso" not in df.columns or "duration_seconds" not in df.columns:
        raise SystemExit("CSV must include 'timestamp_iso' and 'duration_seconds'.")
    df["timestamp"] = pd.to_datetime(df["timestamp_iso"], errors="coerce")
    df["duration"] = pd.to_numeric(df["duration_seconds"], errors="coerce")
    df = df.dropna(subset=["timestamp", "duration"])
    # Optional: keep only completed/skipped for cleaner trends
    if "status" in df.columns:
        df = df[df["status"].isin(["completed", "skipped", ""])]
    return df

def plot_scatter_by_jobtype(df: pd.DataFrame, out_png: str) -> None:
    plt.figure(figsize=(12, 6), dpi=150)
    job_types = sorted(df["job_type"].dropna().unique())
    colors = plt.cm.get_cmap("tab10", len(job_types))

    for i, jt in enumerate(job_types):
        sub = df[df["job_type"] == jt].sort_values("timestamp")
        if sub.empty:
            continue
        plt.scatter(
            sub["timestamp"], sub["duration"],
            s=12, alpha=0.35, color=colors(i), label=jt, edgecolors="none"
        )
        # Rolling median (by observation count) to show trend
        if len(sub) >= 10:
            sub = sub.copy()
            sub["roll_med"] = sub["duration"].rolling(window=max(10, len(sub)//20), min_periods=5).median()
            plt.plot(sub["timestamp"], sub["roll_med"], color=colors(i), linewidth=1.8, alpha=0.9)

    plt.title("Job durations over time")
    plt.xlabel("Time")
    plt.ylabel("Duration (s)")
    plt.legend(loc="upper left", ncols=2, fontsize=8)
    plt.tight_layout()
    plt.savefig(out_png)
    plt.close()

def plot_facet_by_jobtype(df: pd.DataFrame, out_png: str) -> None:
    try:
        import seaborn as sns
    except Exception:
        return  # seaborn not available; skip facet
    job_types = sorted(df["job_type"].dropna().unique())
    if not job_types:
        return

    # Limit facets to avoid huge grids
    max_facets = min(len(job_types), 9)
    keep = set(job_types[:max_facets])
    dff = df[df["job_type"].isin(keep)].copy()

    g = sns.FacetGrid(
        dff.sort_values("timestamp"),
        col="job_type", col_wrap=3, sharex=False, sharey=False, height=3.0
    )
    g.map_dataframe(sns.scatterplot, x="timestamp", y="duration", s=10, alpha=0.35, edgecolor=None)
    g.set_axis_labels("Time", "Duration (s)")
    g.set_titles("{col_name}")
    plt.tight_layout()
    plt.savefig(out_png)
    plt.close()

def main():
    csv_path = sys.argv[1] if len(sys.argv) > 1 else "job_durations.csv"
    out_png = sys.argv[2] if len(sys.argv) > 2 else "job_durations_over_time.png"
    facet_png = os.path.splitext(out_png)[0] + "_facet.png"

    if not os.path.exists(csv_path):
        raise SystemExit(f"File not found: {csv_path}")

    df = load_data(csv_path)
    if df.empty:
        raise SystemExit("No valid rows to plot.")

    plot_scatter_by_jobtype(df, out_png)
    plot_facet_by_jobtype(df, facet_png)
    print(f"Wrote: {out_png}")
    if os.path.exists(facet_png):
        print(f"Wrote: {facet_png}")

if __name__ == "__main__":
    main()
