
#!/usr/bin/env python3
"""General-purpose analysis of genomic feature annotations."""

import argparse
from pathlib import Path

import pandas as pd


def summarize_features(input_file, output_file):
    """Summarize annotated features by genome and feature type."""
    df = pd.read_csv(input_file, sep="\t")

    required = {"Genome_ID", "Feature_type"}
    missing = required - set(df.columns)

    if missing:
        raise ValueError(f"Missing required columns: {sorted(missing)}")

    summary = (
        df.groupby(["Genome_ID", "Feature_type"])
        .size()
        .reset_index(name="Feature_count")
    )

    Path(output_file).parent.mkdir(parents=True, exist_ok=True)
    summary.to_csv(output_file, sep="\t", index=False)
    print(f"Summary saved to: {output_file}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Summarize genomic feature annotations."
    )
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    summarize_features(args.input, args.output)
