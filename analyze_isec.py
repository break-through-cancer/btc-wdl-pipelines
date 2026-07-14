#!/usr/bin/env python3

import argparse
import numpy as np
import pandas as pd
from cyvcf2 import VCF
from matplotlib_venn import venn2
import matplotlib.pyplot as plt
import os
from scipy.stats import pearsonr

#Dont average AF? 
# add all the refs, add all the alts, and then average 


def get_af(v, sample_idx):
    try:
        af = v.format("AF")
        if af is not None:
            return float(af[sample_idx][0])
    except Exception:
        pass

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


def choose_non_pbmc_sample(vcf, label):

    print(f"\n=== SAMPLE SELECTION FOR {label} ===")

    print("\nAll samples in VCF:")
    for i, s in enumerate(vcf.samples):
        print(f"  [{i}] {s}")

    pbmc_samples = [
        s for s in vcf.samples
        if "PBMC" in s or s.split(":")[-1].endswith("PBMC")
    ]

    tumor_samples = [
        s for s in vcf.samples
        if "PBMC" not in s and not s.split(":")[-1].endswith("PBMC")
    ]

    print("\nDetected PBMC/normal samples:")
    for s in pbmc_samples:
        print(f"  - {s}")

    print("\nDetected non-PBMC/tumor samples:")
    for s in tumor_samples:
        print(f"  - {s}")

    if len(tumor_samples) == 0:
        raise ValueError(
            f"No non-PBMC sample found for {label}. "
            f"Available samples: {vcf.samples}"
        )

    chosen = tumor_samples[0]
    sample_idx = vcf.samples.index(chosen)

    print(f"\nUsing sample for AF/DP analysis:")
    print(f"  Sample name : {chosen}")
    print(f"  Sample index: {sample_idx}")

    return sample_idx


def parse_vcf(path, label):
    rows = []
    vcf = VCF(path)

    sample_idx = choose_non_pbmc_sample(vcf, label)

    for i, v in enumerate(vcf, start=1):

        if i % 100000 == 0:
            print(f"{label}: parsed {i:,} variants")

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

    plt.figure(figsize=(8, 8))

    venn2(
        subsets=(
            len(only1),
            len(only2),
            n_shared
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

def save_af_dp_distributions(df, label, outdir):
    af_col = f"AF_{label}"
    dp_col = f"DP_{label}"

    os.makedirs(outdir, exist_ok=True)

    af = df[af_col].dropna()

    plt.figure(figsize=(8,6))
    plt.hist(af, bins=100)
    plt.xlabel("Allele Fraction")
    plt.ylabel("Variant Count")
    plt.title(f"{label} AF Distribution")
    plt.savefig(
        os.path.join(outdir, f"{label}_af_distribution.png"),
        dpi=300,
        bbox_inches="tight"
    )
    plt.close()

    dp = df[dp_col].dropna()

    plt.figure(figsize=(8,6))
    plt.hist(np.log10(dp + 1), bins=100)
    plt.xlabel("log10(DP + 1)")
    plt.ylabel("Variant Count")
    plt.title(f"{label} DP Distribution")
    plt.savefig(
        os.path.join(outdir, f"{label}_dp_distribution.png"),
        dpi=300,
        bbox_inches="tight"
    )
    plt.close()

def save_variant_type_counts(df, label, outdir):

    counts = df["variant_type"].value_counts()

    plt.figure(figsize=(8,6))

    counts.plot(kind="bar")

    plt.ylabel("Count")
    plt.title(f"{label} Variant Types")

    plt.savefig(
        os.path.join(outdir, f"{label}_variant_types.png"),
        dpi=300,
        bbox_inches="tight"
    )

    plt.close()

def save_af_by_variant_type(df, label, outdir):

    af_col = f"AF_{label}"

    data = []

    labels = []

    for vt in sorted(df["variant_type"].unique()):

        vals = (
            df.loc[df["variant_type"] == vt, af_col]
            .dropna()
        )

        if len(vals):
            data.append(vals)
            labels.append(vt)

    plt.figure(figsize=(10,6))

    plt.boxplot(
        data,
        tick_labels=labels,
        showfliers=False
    )

    plt.ylabel("AF")
    plt.title(f"{label} AF by Variant Type")

    plt.savefig(
        os.path.join(outdir, f"{label}_af_by_type.png"),
        dpi=300,
        bbox_inches="tight"
    )

    plt.close()

def save_dp_by_variant_type(df, label, outdir):

    dp_col = f"DP_{label}"

    data = []
    labels = []

    for vt in sorted(df["variant_type"].unique()):

        vals = (
            np.log10(
                df.loc[df["variant_type"] == vt, dp_col]
                .dropna()
                + 1
            )
        )

        if len(vals):
            data.append(vals)
            labels.append(vt)

    plt.figure(figsize=(10,6))

    plt.boxplot(
        data,
        tick_labels=labels,
        showfliers=False
    )

    plt.ylabel("log10(DP+1)")
    plt.title(f"{label} DP by Variant Type")

    plt.savefig(
        os.path.join(outdir, f"{label}_dp_by_type.png"),
        dpi=300,
        bbox_inches="tight"
    )

    plt.close()

def save_shared_variant_plots(df1, df2, name1, name2, outdir):

    merged = df1.merge(
        df2,
        on="variant_id",
        how="inner",
        suffixes=(f"_{name1}", f"_{name2}")
    )

    print(f"\nShared variants for AF comparison: {len(merged):,}")

    af1 = merged[f"AF_{name1}"]
    af2 = merged[f"AF_{name2}"]

    mask = (
        af1.notna() &
        af2.notna()
    )

    af1 = af1[mask]
    af2 = af2[mask]

    if len(af1):

        corr = pearsonr(af1, af2)[0]

        diff = np.abs(af1 - af2)

        print(f"Shared variants with AF values: {len(af1):,}")
        print(f"Pearson r: {corr:.4f}")
        print(f"Median |AF1-AF2|: {np.median(diff):.4f}")
        print(f"Mean |AF1-AF2|: {np.mean(diff):.4f}")

        # ============================================================
        # Hexbin AF concordance plot
        # ============================================================

        plt.figure(figsize=(8, 8))

        hb = plt.hexbin(
            af1,
            af2,
            gridsize=100,
            bins="log",
            mincnt=1
        )

        plt.plot(
            [0, 1],
            [0, 1],
            "r--",
            linewidth=1,
            label="y=x"
        )

        plt.colorbar(hb, label="log10(count)")

        plt.xlabel(f"{name1} AF")
        plt.ylabel(f"{name2} AF")

        plt.title(
            f"Shared Variant AF Concordance\n"
            f"Pearson r={corr:.3f}"
        )

        plt.legend()

        plt.savefig(
            os.path.join(outdir, "shared_af_hexbin.png"),
            dpi=300,
            bbox_inches="tight"
        )

        plt.close()

        # ============================================================
        # AF difference histogram
        # ============================================================

        plt.figure(figsize=(8, 6))

        plt.hist(
            diff,
            bins=100
        )

        plt.axvline(
            np.median(diff),
            linestyle="--",
            linewidth=2,
            label=f"median={np.median(diff):.3f}"
        )

        plt.xlabel("|AF1 - AF2|")
        plt.ylabel("Count")

        plt.title(
            "Shared Variant AF Differences"
        )

        plt.legend()

        plt.savefig(
            os.path.join(outdir, "shared_af_difference.png"),
            dpi=300,
            bbox_inches="tight"
        )

        plt.close()

        # ============================================================
        # AF density distributions
        # ============================================================

        plt.figure(figsize=(8, 6))

        plt.hist(
            af1,
            bins=100,
            density=True,
            alpha=0.5,
            label=name1
        )

        plt.hist(
            af2,
            bins=100,
            density=True,
            alpha=0.5,
            label=name2
        )

        plt.xlabel("Allele Fraction")
        plt.ylabel("Density")

        plt.title(
            "Shared Variant AF Distribution"
        )

        plt.legend()

        plt.savefig(
            os.path.join(outdir, "shared_af_density.png"),
            dpi=300,
            bbox_inches="tight"
        )

        plt.close()
def save_overlap_type_plot(df1, df2, name1, name2, outdir):

    merged = df1.merge(
        df2,
        on="variant_id",
        how="outer",
        indicator=True,
        suffixes=(f"_{name1}", f"_{name2}")
    )

    def classify(x):
        if x == "both":
            return "shared"
        elif x == "left_only":
            return f"{name1}_only"
        else:
            return f"{name2}_only"

    merged["status"] = merged["_merge"].apply(classify)

    merged["variant_type_final"] = (
        merged.get(f"variant_type_{name1}")
        .combine_first(
            merged.get(f"variant_type_{name2}")
        )
    )

    counts = (
        merged
        .groupby(
            ["status", "variant_type_final"]
        )
        .size()
        .unstack(fill_value=0)
    )

    counts.plot(
        kind="bar",
        figsize=(10,6)
    )

    plt.ylabel("Variant Count")

    plt.title(
        "Overlap Status by Variant Type"
    )

    plt.savefig(
        os.path.join(
            outdir,
            "overlap_variant_types.png"
        ),
        dpi=300,
        bbox_inches="tight"
    )

    plt.close()
def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("--vcf1", required=True)
    parser.add_argument("--vcf2", required=True)

    parser.add_argument("--name1", default="file1")
    parser.add_argument("--name2", default="file2")

    parser.add_argument(
        "--plot_dir",
        default="plots",
        help="Directory to save output plots"
    )

    args = parser.parse_args()
    print_header_info(args.vcf1, args.name1)
    print_header_info(args.vcf2, args.name2)

    print("\nReading VCFs...")

    df1 = parse_vcf(args.vcf1, args.name1)
    df2 = parse_vcf(args.vcf2, args.name2)

    print_basic_summary(df1, args.name1)
    print_basic_summary(df2, args.name2)

    print_overlap_summary(df1, df2, args.name1, args.name2)
    save_af_dp_distributions(
        df1,
        args.name1,
        args.plot_dir
        )

    save_af_dp_distributions(
        df2,
        args.name2,
        args.plot_dir
    )

    save_variant_type_counts(
        df1,
        args.name1,
        args.plot_dir
    )

    save_variant_type_counts(
        df2,
        args.name2,
        args.plot_dir
    )

    save_af_by_variant_type(
        df1,
        args.name1,
        args.plot_dir
    )

    save_af_by_variant_type(
        df2,
        args.name2,
        args.plot_dir
    )

    save_dp_by_variant_type(
        df1,
        args.name1,
        args.plot_dir
    )

    save_dp_by_variant_type(
        df2,
        args.name2,
        args.plot_dir
    )

    save_shared_variant_plots(
        df1,
        df2,
        args.name1,
        args.name2,
        args.plot_dir
    )

    save_overlap_type_plot(
        df1,
        df2,
        args.name1,
        args.name2,
        args.plot_dir
    )


if __name__ == "__main__":
    main()