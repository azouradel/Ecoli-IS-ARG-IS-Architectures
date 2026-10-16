
#!/usr/bin/env python3
"""General-purpose visualization of genomic feature counts."""

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import pandas as pd


def visualize(input_file, output_file, top_n):
    """Plot the most frequent genomic feature types."""
    df = pd.read_csv(input_file, sep="\t")

    required = {"Feature_type", "Feature_count"}
    missing = required - set(df.columns)

    if missing:
        raise ValueError(
            f"Missing required columns: {sorted(missing)}"
        )

    df["Feature_count"] = pd.to_numeric(
        df["Feature_count"], errors="raise"
    )

    if df["Feature_count"].isna().any():
        raise ValueError("Feature_count contains missing values.")

    if (df["Feature_count"] < 0).any():
        raise ValueError("Feature_count must be non-negative.")

    summary = (
        df.groupby("Feature_type")["Feature_count"]
        .sum()
        .sort_values(ascending=False)
        .head(top_n)
        .sort_values()
    )

    if summary.empty:
        raise ValueError("No features available for plotting.")

    fig, ax = plt.subplots(
        figsize=(9, max(4, len(summary) * 0.45))
    )

    summary.plot.barh(ax=ax)

    ax.set_xlabel("Number of annotations")
    ax.set_ylabel("Genomic feature type")
    ax.set_title("Distribution of genomic features")

    fig.tight_layout()

    output_path = Path(output_file)
    output_path.parent.mkdir(
        parents=True, exist_ok=True
    )

    fig.savefig(output_path, dpi=300)
    plt.close(fig)

    print(f"Figure saved to: {output_path}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Visualize genomic feature counts."
    )

    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument(
        "--top", type=int, default=20
    )

    args = parser.parse_args()

    if args.top < 1:
        parser.error("--top must be at least 1")

    visualize(args.input, args.output, args.top)
