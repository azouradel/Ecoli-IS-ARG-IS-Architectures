
#!/usr/bin/env python3
"""General statistical comparison of genomic feature lengths."""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf
from scipy.stats import mannwhitneyu


def analyze(input_file, output_file):
    df = pd.read_csv(input_file, sep="\t")

    required = {"Genome_ID", "Length_bp", "Location"}
    missing = required - set(df.columns)

    if missing:
        raise ValueError(f"Missing columns: {sorted(missing)}")

    df = df.copy()
    df["Length_bp"] = pd.to_numeric(
        df["Length_bp"], errors="coerce"
    )

    df = df.dropna(
        subset=["Genome_ID", "Length_bp", "Location"]
    )
    df = df[df["Length_bp"] >= 0].copy()

    allowed = {"Plasmid", "Chromosomal"}
    if not set(df["Location"]).issubset(allowed):
        raise ValueError(
            "Location must be Plasmid or Chromosomal."
        )

    df["Plasmid"] = (
        df["Location"] == "Plasmid"
    ).astype(int)

    if df["Plasmid"].nunique() != 2:
        raise ValueError("Both location groups are required.")

    df["log_Length"] = np.log1p(df["Length_bp"])

    plasmid = df.loc[
        df["Plasmid"] == 1, "Length_bp"
    ]
    chromosomal = df.loc[
        df["Plasmid"] == 0, "Length_bp"
    ]

    u_stat, mw_p = mannwhitneyu(
        plasmid, chromosomal, alternative="two-sided"
    )

    model = smf.gee(
        "log_Length ~ Plasmid",
        groups="Genome_ID",
        data=df,
        family=sm.families.Gaussian(),
        cov_struct=sm.cov_struct.Exchangeable(),
    )
    fit = model.fit()

    beta = float(fit.params["Plasmid"])
    se = float(fit.bse["Plasmid"])

    result = pd.DataFrame([{
        "Plasmid_N": len(plasmid),
        "Chromosomal_N": len(chromosomal),
        "Plasmid_median_bp": plasmid.median(),
        "Chromosomal_median_bp": chromosomal.median(),
        "Mann_Whitney_U": u_stat,
        "Mann_Whitney_P": mw_p,
        "GEE_beta": beta,
        "GEE_ratio": np.exp(beta),
        "CI95_lower": np.exp(beta - 1.96 * se),
        "CI95_upper": np.exp(beta + 1.96 * se),
        "GEE_P": float(fit.pvalues["Plasmid"]),
    }])

    Path(output_file).parent.mkdir(
        parents=True, exist_ok=True
    )
    result.to_csv(output_file, sep="\t", index=False)
    print(result.to_string(index=False))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Compare genomic feature lengths."
    )
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    analyze(args.input, args.output)
