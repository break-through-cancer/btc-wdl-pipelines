#!/usr/bin/env python3

import argparse
import gzip
import os
from pathlib import Path
from collections import defaultdict

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def open_vcf(path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def get_variant_type(ref, alt):
    if len(ref) == 1 and len(alt) == 1:
        return "SNV"
    elif len(ref) != len(alt):
        return "INDEL"
    else:
        return "MNV_or_complex"


def is_transition(ref, alt):
    return (ref, alt) in {
        ("A", "G"), ("G", "A"),
        ("C", "T"), ("T", "C")
    }


def parse_sample_metrics(info, fmt_keys, sample_vals, alt_index):
    info_dict = {}
    for item in info.split(";"):
        if "=" in item:
            k, v = item.split("=", 1)
            info_dict[k] = v
        else:
            info_dict[item] = True

    fmt = dict(zip(fmt_keys, sample_vals))

    dp = np.nan
    af = np.nan
    tlod = np.nan

    if "DP" in fmt and fmt["DP"] not in [".", ""]:
        try:
            dp = float(fmt["DP"])
        except ValueError:
            pass
    elif "DP" in info_dict:
        try:
            dp = float(str(info_dict["DP"]).split(",")[0])
        except ValueError:
            pass

    if "AF" in fmt and fmt["AF"] not in [".", ""]:
        vals = fmt["AF"].split(",")
        if alt_index < len(vals):
            try:
                af = float(vals[alt_index])
            except ValueError:
                pass

    elif "AD" in fmt and fmt["AD"] not in [".", ""]:
        vals = fmt["AD"].split(",")
        try:
            ref_count = float(vals[0])
            alt_count = float(vals[alt_index + 1])
            total = ref_count + alt_count
            if total > 0:
                af = alt_count / total
        except Exception:
            pass

    if "TLOD" in info_dict:
        vals = str(info_dict["TLOD"]).split(",")
        if alt_index < len(vals):
            try:
                tlod = float(vals[alt_index])
            except ValueError:
                pass

    return af, dp, tlod


def parse_vcf(path, sample_name=None):
    rows = []
    header_sample = None

    with open_vcf(path) as f:
        for line in f:
            if line.startswith("##"):
                continue

            if line.startswith("#CHROM"):
                parts = line.rstrip("\n").split("\t")
                if len(parts) > 9:
                    header_sample = parts[9]
                continue

            parts = line.rstrip("\n").split("\t")
            if len(parts) < 8:
                continue

            chrom, pos, vid, ref, alts, qual, filt, info = parts[:8]

            fmt_keys = []
            sample_vals = []

            if len(parts) > 9:
                fmt_keys = parts[8].split(":")
                sample_vals = parts[9].split(":")

            for alt_index, alt in enumerate(alts.split(",")):
                variant_id = f"{chrom}:{pos}:{ref}:{alt}"
                vtype = get_variant_type(ref, alt)

                af, dp, tlod = parse_sample_metrics(
                    info,
                    fmt_keys,
                    sample_vals,
                    alt_index
                )

                rows.append({
                    "variant_id": variant_id,
                    "chrom": chrom,
                    "pos": int(pos),
                    "ref": ref,
                    "alt": alt,
                    "filter": filt,
                    "is_pass": filt in ["PASS", "."],
                    "variant_type": vtype,
                    "AF": af,
                    "DP": dp,
                    "TLOD": tlod,
                    "is_transition": is_transition(ref, alt) if vtype == "SNV" else np.nan,
                })

    if sample_name is None:
        sample_name = header_sample or Path(path).name.replace(".vcf.gz", "").replace(".vcf", "")

    df = pd.DataFrame(rows)
    df["sample"] = sample_name
    return df


def summarize_sample(df):
    snvs = df[df["variant_type"] == "SNV"]

    transitions = snvs["is_transition"].sum()
    transversions = len(snvs) - transitions

    titv = np.nan
    if transversions > 0:
        titv = transitions / transversions

    return {
        "sample": df["sample"].iloc[0],
        "total_variants": len(df),
        "pass_variants": df["is_pass"].sum(),
        "snv_count": (df["variant_type"] == "SNV").sum(),
        "indel_count": (df["variant_type"] == "INDEL").sum(),
        "mnv_or_complex_count": (df["variant_type"] == "MNV_or_complex").sum(),
        "median_AF": df["AF"].median(),
        "mean_AF": df["AF"].mean(),
        "median_DP": df["DP"].median(),
        "mean_DP": df["DP"].mean(),
        "median_TLOD": df["TLOD"].median(),
        "mean_TLOD": df["TLOD"].mean(),
        "TiTv": titv,
    }


def save_hist(df, column, outpath, title, xlabel, log_y=False):
    vals = df[column].dropna()
    if len(vals) == 0:
        return

    plt.figure(figsize=(8, 6))
    plt.hist(vals, bins=80)
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel("Variant count")
    if log_y:
        plt.yscale("log")
    plt.tight_layout()
    plt.savefig(outpath, dpi=200)
    plt.close()


def save_bar(series, outpath, title, ylabel):
    plt.figure(figsize=(10, 6))
    series.plot(kind="bar")
    plt.title(title)
    plt.ylabel(ylabel)
    plt.xticks(rotation=45, ha="right")
    plt.tight_layout()
    plt.savefig(outpath, dpi=200)
    plt.close()


def save_heatmap(matrix, outpath, title):
    plt.figure(figsize=(10, 8))
    plt.imshow(matrix, aspect="auto")
    plt.colorbar(label="Jaccard index")
    plt.xticks(range(len(matrix.columns)), matrix.columns, rotation=90)
    plt.yticks(range(len(matrix.index)), matrix.index)
    plt.title(title)
    plt.tight_layout()
    plt.savefig(outpath, dpi=200)
    plt.close()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--vcf_dir", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--merged_vcf", default=None)
    parser.add_argument("--pass_only", action="store_true")
    args = parser.parse_args()

    outdir = Path(args.outdir)
    plotdir = outdir / "plots"
    outdir.mkdir(parents=True, exist_ok=True)
    plotdir.mkdir(parents=True, exist_ok=True)

    vcf_paths = sorted(
        list(Path(args.vcf_dir).glob("*.vcf")) +
        list(Path(args.vcf_dir).glob("*.vcf.gz"))
    )

    if not vcf_paths:
        raise ValueError(f"No VCFs found in {args.vcf_dir}")

    print(f"Found {len(vcf_paths)} VCFs")

    all_dfs = []

    for path in vcf_paths:
        print(f"Parsing {path}")
        df = parse_vcf(path)

        if args.pass_only:
            df = df[df["is_pass"]].copy()

        all_dfs.append(df)

        sample = df["sample"].iloc[0]

        save_hist(
            df,
            "AF",
            plotdir / f"{sample}.AF_hist.png",
            f"{sample} AF Distribution",
            "Allele fraction",
            log_y=True
        )

        save_hist(
            df,
            "DP",
            plotdir / f"{sample}.DP_hist.png",
            f"{sample} DP Distribution",
            "Depth",
            log_y=True
        )

        save_hist(
            df,
            "TLOD",
            plotdir / f"{sample}.TLOD_hist.png",
            f"{sample} TLOD Distribution",
            "TLOD",
            log_y=True
        )

    all_variants = pd.concat(all_dfs, ignore_index=True)

    all_variants.to_csv(outdir / "all_sample_variants.long.csv", index=False)

    summary = pd.DataFrame([summarize_sample(df) for df in all_dfs])
    summary.to_csv(outdir / "per_sample_summary.csv", index=False)

    print("\nPer-sample summary:")
    print(summary)

    save_bar(
        summary.set_index("sample")["total_variants"],
        plotdir / "total_variants_per_sample.png",
        "Total Variants per Sample",
        "Variant count"
    )

    save_bar(
        summary.set_index("sample")["pass_variants"],
        plotdir / "pass_variants_per_sample.png",
        "PASS Variants per Sample",
        "PASS variant count"
    )

    type_counts = (
        all_variants
        .groupby(["sample", "variant_type"])
        .size()
        .reset_index(name="count")
    )
    type_counts.to_csv(outdir / "variant_type_counts.csv", index=False)

    presence = (
        all_variants
        .assign(present=1)
        .pivot_table(
            index="variant_id",
            columns="sample",
            values="present",
            aggfunc="max",
            fill_value=0
        )
        .astype(int)
    )

    presence.to_csv(outdir / "variant_presence_matrix.csv")

    presence["support_count"] = presence.sum(axis=1)

    support_hist = (
        presence["support_count"]
        .value_counts()
        .sort_index()
        .reset_index()
    )
    support_hist.columns = ["n_samples_present", "variant_count"]
    support_hist.to_csv(outdir / "support_count_histogram.csv", index=False)

    save_bar(
        support_hist.set_index("n_samples_present")["variant_count"],
        plotdir / "support_count_histogram.png",
        "How Many Samples Each Variant Appears In",
        "Variant count"
    )

    samples = [df["sample"].iloc[0] for df in all_dfs]
    jaccard = pd.DataFrame(index=samples, columns=samples, dtype=float)
    overlap_rows = []

    sample_sets = {
        sample: set(all_variants.loc[all_variants["sample"] == sample, "variant_id"])
        for sample in samples
    }

    for s1 in samples:
        for s2 in samples:
            set1 = sample_sets[s1]
            set2 = sample_sets[s2]

            shared = len(set1 & set2)
            union = len(set1 | set2)
            only_s1 = len(set1 - set2)
            only_s2 = len(set2 - set1)

            jac = shared / union if union else np.nan
            jaccard.loc[s1, s2] = jac

            overlap_rows.append({
                "sample_1": s1,
                "sample_2": s2,
                "sample_1_total": len(set1),
                "sample_2_total": len(set2),
                "shared": shared,
                "sample_1_only": only_s1,
                "sample_2_only": only_s2,
                "union": union,
                "jaccard": jac,
                "shared_fraction_of_sample_1": shared / len(set1) if len(set1) else np.nan,
                "shared_fraction_of_sample_2": shared / len(set2) if len(set2) else np.nan,
            })

    jaccard.to_csv(outdir / "pairwise_jaccard_matrix.csv")
    pd.DataFrame(overlap_rows).to_csv(outdir / "pairwise_overlap_counts.csv", index=False)

    save_heatmap(
        jaccard,
        plotdir / "pairwise_jaccard_heatmap.png",
        "Pairwise Jaccard Similarity"
    )

    if args.merged_vcf:
        print(f"\nParsing merged VCF: {args.merged_vcf}")
        merged_df = parse_vcf(args.merged_vcf, sample_name="MERGED")

        if args.pass_only:
            merged_df = merged_df[merged_df["is_pass"]].copy()

        merged_set = set(merged_df["variant_id"])

        merged_rows = []

        for sample in samples:
            sample_set = sample_sets[sample]
            shared = sample_set & merged_set

            merged_rows.append({
                "sample": sample,
                "sample_total": len(sample_set),
                "merged_total": len(merged_set),
                "shared_with_merged": len(shared),
                "sample_only": len(sample_set - merged_set),
                "merged_only": len(merged_set - sample_set),
                "fraction_sample_retained_in_merged": len(shared) / len(sample_set) if sample_set else np.nan,
                "fraction_merged_seen_in_sample": len(shared) / len(merged_set) if merged_set else np.nan,
            })

        merged_compare = pd.DataFrame(merged_rows)
        merged_compare.to_csv(outdir / "merged_vs_each_sample_retention.csv", index=False)

        save_bar(
            merged_compare.set_index("sample")["fraction_sample_retained_in_merged"],
            plotdir / "fraction_sample_retained_in_merged.png",
            "Fraction of Each Sample's Variants Retained in Merged VCF",
            "Fraction retained"
        )

        original_support = presence[["support_count"]].copy()
        original_support["in_merged"] = original_support.index.isin(merged_set)

        retention_by_support = (
            original_support
            .groupby("support_count")["in_merged"]
            .agg(["sum", "count"])
            .reset_index()
        )
        retention_by_support["retention_rate"] = (
            retention_by_support["sum"] / retention_by_support["count"]
        )

        retention_by_support.columns = [
            "n_original_samples_present",
            "n_retained_in_merged",
            "n_original_variants",
            "retention_rate"
        ]

        retention_by_support.to_csv(
            outdir / "merged_retention_by_original_support_count.csv",
            index=False
        )

        save_bar(
            retention_by_support.set_index("n_original_samples_present")["retention_rate"],
            plotdir / "merged_retention_by_support_count.png",
            "Merged Retention Rate by Original Sample Support Count",
            "Retention rate"
        )

        merged_small = merged_df[["variant_id", "AF", "DP", "TLOD"]].rename(
            columns={
                "AF": "AF_MERGED",
                "DP": "DP_MERGED",
                "TLOD": "TLOD_MERGED"
            }
        )

        af_compare_rows = []

        for df in all_dfs:
            sample = df["sample"].iloc[0]

            sample_small = df[["variant_id", "AF", "DP", "TLOD"]].rename(
                columns={
                    "AF": f"AF_{sample}",
                    "DP": f"DP_{sample}",
                    "TLOD": f"TLOD_{sample}"
                }
            )

            merged_sample = sample_small.merge(
                merged_small,
                on="variant_id",
                how="inner"
            )

            merged_sample.to_csv(
                outdir / f"{sample}.merged_shared_metrics.csv",
                index=False
            )

            af1 = merged_sample[f"AF_{sample}"]
            af2 = merged_sample["AF_MERGED"]

            mask = af1.notna() & af2.notna()

            if mask.sum() > 1:
                corr = np.corrcoef(af1[mask], af2[mask])[0, 1]
                median_abs_diff = np.median(np.abs(af1[mask] - af2[mask]))
            else:
                corr = np.nan
                median_abs_diff = np.nan

            af_compare_rows.append({
                "sample": sample,
                "shared_variants_with_AF": int(mask.sum()),
                "AF_pearson_r_with_merged": corr,
                "median_abs_AF_diff_with_merged": median_abs_diff,
            })

            if mask.sum() > 0:
                plt.figure(figsize=(7, 7))
                plt.scatter(af1[mask], af2[mask], s=2, alpha=0.25)
                plt.xlabel(f"{sample} AF")
                plt.ylabel("Merged AF")
                plt.title(f"{sample} vs Merged AF")
                plt.tight_layout()
                plt.savefig(plotdir / f"{sample}.AF_vs_merged_scatter.png", dpi=200)
                plt.close()

                delta = af2[mask] - af1[mask]
                plt.figure(figsize=(8, 6))
                plt.hist(delta, bins=80)
                plt.xlabel("Merged AF - Sample AF")
                plt.ylabel("Variant count")
                plt.title(f"{sample} AF Shift in Merged VCF")
                plt.tight_layout()
                plt.savefig(plotdir / f"{sample}.AF_delta_merged_minus_sample.png", dpi=200)
                plt.close()

        pd.DataFrame(af_compare_rows).to_csv(
            outdir / "merged_af_comparison_summary.csv",
            index=False
        )

    print(f"\nDone. Results written to: {outdir}")


if __name__ == "__main__":
    main()