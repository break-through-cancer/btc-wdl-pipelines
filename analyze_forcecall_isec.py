#!/usr/bin/env python3

import argparse
import numpy as np
import pandas as pd
from cyvcf2 import VCF
from matplotlib_venn import venn2
import matplotlib.pyplot as plt


def get_af(v, sample_idx):
    # FORMAT/AF for selected sample
    try:
        af = v.format("AF")
        if af is not None:
            return float(af[sample_idx][0])
    except Exception:
        pass

    # FORMAT/AD for selected sample
    try:
        ad = v.format("AD")
        if ad is not None and len(ad[sample_idx]) >= 2:
            ref_depth = ad[sample_idx][0]
            alt_depth = ad[sample_idx][1]
            total = ref_depth + alt_depth
            if total > 0:
                return alt_depth / total
    except Exception:
        pass

    return np.nan


def get_dp(v, sample_idx):
    # FORMAT/DP for selected sample
    try:
        dp = v.format("DP")
        if dp is not None:
            return float(dp[sample_idx][0])
    except Exception:
        pass

    return np.nan


def variant_type(v):
    ref = v.REF
    alts = v.ALT

    if len(alts) != 1:
        return "MULTIALLELIC"

    alt = alts[0]

    if len(ref) == 1 and len(alt) == 1:
        return "SNV"
    elif len(ref) != len(alt):
        return "INDEL"
    else:
        return "MNV_or_complex"


def print_header_info(path, label):
    vcf = VCF(path)

    print(f"\n{'='*60}")
    print(f"{label} HEADER INFO")
    print(f"{'='*60}")

    print("\nSamples:")
    for s in vcf.samples:
        print(f"  - {s}")

    info_fields = []
    format_fields = []
    contigs = []

    for h in vcf.header_iter():
        try:
            if h["HeaderType"] == "INFO":
                info_fields.append(h["ID"])
            elif h["HeaderType"] == "FORMAT":
                format_fields.append(h["ID"])
            elif h["HeaderType"] == "CONTIG":
                contigs.append(h["ID"])
        except Exception:
            pass

    print("\nINFO fields:")
    print(", ".join(sorted(info_fields)))

    print("\nFORMAT fields:")
    print(", ".join(sorted(format_fields)))

    print("\nContigs:")
    print(contigs[:10])
    print(f"Number of contigs: {len(contigs)}")

    if any(c.startswith("chr") for c in contigs):
        print("Detected chr-prefixed contigs")
    else:
        print("Detected non-chr-prefixed contigs")


def parse_vcf(path, label, tumor_sample):
    rows = []
    vcf = VCF(path)

    if tumor_sample not in vcf.samples:
        raise ValueError(
            f"Tumor sample '{tumor_sample}' not found in {label}. "
            f"Available samples: {vcf.samples}"
        )

    sample_idx = vcf.samples.index(tumor_sample)

    print(f"\nUsing tumor sample for {label}: {tumor_sample}")
    print(f"Sample index: {sample_idx}")

    for v in vcf:
        alt = ",".join(v.ALT)
        key = f"{v.CHROM}:{v.POS}:{v.REF}:{alt}"

        rows.append({
            "variant_id": key,
            "CHROM": v.CHROM,
            "POS": v.POS,
            "REF": v.REF,
            "ALT": alt,
            f"QUAL_{label}": v.QUAL,
            f"FILTER_{label}": v.FILTER or "PASS",
            f"AF_{label}": get_af(v, sample_idx),
            f"DP_{label}": get_dp(v, sample_idx),
            "variant_type": variant_type(v),
        })

    return pd.DataFrame(rows)


def print_basic_summary(df, label):
    print(f"\n{'='*60}")
    print(f"{label} BASIC SUMMARY")
    print(f"{'='*60}")

    print(f"\nTotal variants: {len(df):,}")

    print("\nVariant types:")
    print(df["variant_type"].value_counts())

    filter_col = f"FILTER_{label}"

    print("\nFilter status:")
    print(df[filter_col].value_counts())

    print("\nTop chromosomes:")
    print(df["CHROM"].value_counts().head(10))

    multiallelic = df["ALT"].str.contains(",").sum()
    print(f"\nMultiallelic variants: {multiallelic:,}")

    pass_count = (df[filter_col] == "PASS").sum()
    print(f"PASS variants: {pass_count:,}")
    print(f"PASS rate: {100 * pass_count / len(df):.2f}%")

    af_col = f"AF_{label}"
    dp_col = f"DP_{label}"

    print("\nAF summary:")
    print(df[af_col].describe(percentiles=[0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99]))

    print("\nDP summary:")
    print(df[dp_col].describe(percentiles=[0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99]))

    print("\nAF by variant type:")
    print(
        df.groupby("variant_type")[af_col]
        .describe(percentiles=[0.5, 0.9, 0.99])
        [["count", "mean", "50%", "90%", "99%", "max"]]
    )

    print("\nDP by variant type:")
    print(
        df.groupby("variant_type")[dp_col]
        .describe(percentiles=[0.5, 0.9, 0.99])
        [["count", "mean", "50%", "90%", "99%", "max"]]
    )

def print_overlap_summary(df1, df2, name1, name2):
    print(f"\n{'='*60}")
    print("VARIANT OVERLAP SUMMARY")
    print(f"{'='*60}")

    set1 = set(df1["variant_id"])
    set2 = set(df2["variant_id"])

    shared = set1 & set2
    only1 = set1 - set2
    only2 = set2 - set1

    n1 = len(set1)
    n2 = len(set2)
    n_shared = len(shared)

    print(f"\n{name1} total variants: {n1:,}")
    print(f"{name2} total variants: {n2:,}")
    print(f"Shared variants: {n_shared:,}")

    print(f"\nUnique to {name1}: {len(only1):,}")
    print(f"Unique to {name2}: {len(only2):,}")

    print(f"\nPercent shared:")
    print(f"{name1}: {100 * n_shared / n1:.2f}%")
    print(f"{name2}: {100 * n_shared / n2:.2f}%")

    jaccard = n_shared / len(set1 | set2)

    print(f"\nJaccard similarity: {jaccard:.4f}")

        # -------------------------
    # Venn diagram
    # -------------------------

    plt.figure(figsize=(8, 8))

    venn2(
        subsets=(
            len(only1),   # unique to set1
            len(only2),   # unique to set2
            n_shared      # shared
        ),
        set_labels=(name1, name2)
    )

    plt.title("Variant Overlap")

    plt.savefig(
        f"{name1}_vs_{name2}_venn.png",
        dpi=300,
        bbox_inches="tight"
    )

    print(f"\nSaved venn diagram:")
    print(f"{name1}_vs_{name2}_venn.png")

    # ---- Variant type overlap ----

    merged = df1.merge(
        df2,
        on="variant_id",
        how="outer",
        indicator=True,
        suffixes=(f"_{name1}", f"_{name2}")
    )

    def overlap_label(x):
        if x == "both":
            return "shared"
        elif x == "left_only":
            return f"{name1}_only"
        else:
            return f"{name2}_only"

    merged["overlap_status"] = merged["_merge"].apply(overlap_label)

    # use variant type from whichever side exists
    merged["variant_type_final"] = (
        merged.get("variant_type_" + name1)
        .combine_first(merged.get("variant_type_" + name2))
    )

    print("\nOverlap by variant type:")
    print(
        merged.groupby(["overlap_status", "variant_type_final"])
        .size()
        .reset_index(name="count")
    )

def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("--vcf1", required=True)
    parser.add_argument("--vcf2", required=True)

    parser.add_argument("--name1", default="file1")
    parser.add_argument("--name2", default="file2")

    parser.add_argument("--tumor_sample1", required=True)
    parser.add_argument("--tumor_sample2", required=True)

    args = parser.parse_args()

    print_header_info(args.vcf1, args.name1)
    print_header_info(args.vcf2, args.name2)

    print("\nReading VCFs...")

    df1 = parse_vcf(args.vcf1, args.name1, args.tumor_sample1)
    df2 = parse_vcf(args.vcf2, args.name2, args.tumor_sample2)

    print_basic_summary(df1, args.name1)
    print_basic_summary(df2, args.name2)
    print_overlap_summary(df1, df2, args.name1, args.name2)


if __name__ == "__main__":
    main()