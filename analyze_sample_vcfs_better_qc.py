#!/usr/bin/env python3

import argparse
import gzip
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def open_vcf(path):
    """Open plain-text or gzipped VCF."""
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def get_variant_type(ref, alt):
    """Classify a biallelic record by REF/ALT length."""
    if len(ref) == 1 and len(alt) == 1:
        return "SNV"
    elif len(ref) != len(alt):
        return "INDEL"
    else:
        return "MNV_or_complex"


def is_transition(ref, alt):
    """Return True for A<->G or C<->T SNVs."""
    return (ref, alt) in {
        ("A", "G"), ("G", "A"),
        ("C", "T"), ("T", "C")
    }


def parse_sample_metrics(info, fmt_keys, sample_vals, alt_index):
    """
    Parse AF, DP, and TLOD for one ALT allele.

    AF is pulled from FORMAT/AF if present. If AF is missing, it is calculated
    from FORMAT/AD as alt_depth / (ref_depth + alt_depth).

    DP is pulled from FORMAT/DP first, then INFO/DP.

    TLOD is a Mutect2 INFO field. It can contain one value per ALT allele, so
    we use alt_index to select the matching TLOD for multiallelic records.
    """
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
    """Parse one VCF into one row per ALT allele."""
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
    """Return one summary row for a sample."""
    snvs = df[df["variant_type"] == "SNV"]
    transitions = snvs["is_transition"].sum()
    transversions = len(snvs) - transitions
    titv = transitions / transversions if transversions > 0 else np.nan

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


def save_hist(df, column, outpath, title, xlabel, log_y=False, transform=None):
    """Save a single-sample histogram."""
    vals = df[column].dropna()
    if transform == "log10_plus_1":
        vals = np.log10(vals + 1)
    vals = vals[np.isfinite(vals)]
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
    plt.savefig(outpath, dpi=250)
    plt.close()


def save_bar(series, outpath, title, ylabel):
    """Save a barplot from a pandas Series."""
    plt.figure(figsize=(10, 6))
    series.plot(kind="bar")
    plt.title(title)
    plt.ylabel(ylabel)
    plt.xticks(rotation=45, ha="right")
    plt.tight_layout()
    plt.savefig(outpath, dpi=250)
    plt.close()


def save_overlay_density_by_sample(all_variants, column, outpath, title, xlabel, transform=None, max_samples=20):
    """
    Overlay one normalized histogram per sample on the same axes.

    density=True means each sample's histogram is scaled so the total area under
    that sample's curve is 1. This is useful when samples have different numbers
    of variants because it compares distribution shape instead of raw counts.
    """
    plt.figure(figsize=(9, 6))

    n_plotted = 0
    for sample, sub in all_variants.groupby("sample"):
        vals = sub[column].dropna()
        if transform == "log10_plus_1":
            vals = np.log10(vals + 1)
        vals = vals[np.isfinite(vals)]
        if len(vals) == 0:
            continue

        plt.hist(
            vals,
            bins=80,
            density=True,
            histtype="step",
            linewidth=1.5,
            label=sample
        )
        n_plotted += 1
        if n_plotted >= max_samples:
            break

    if n_plotted == 0:
        plt.close()
        return

    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel("Density")
    plt.legend(fontsize=7, frameon=False)
    plt.tight_layout()
    plt.savefig(outpath, dpi=250)
    plt.close()


def save_boxplot_by_sample(all_variants, column, outpath, title, ylabel, transform=None):
    """Save a compact sample-level boxplot for AF/DP/TLOD distributions."""
    samples = []
    data = []

    for sample, sub in all_variants.groupby("sample"):
        vals = sub[column].dropna()
        if transform == "log10_plus_1":
            vals = np.log10(vals + 1)
        vals = vals[np.isfinite(vals)]
        if len(vals) > 0:
            samples.append(sample)
            data.append(vals)

    if not data:
        return

    plt.figure(figsize=(max(10, len(samples) * 0.5), 6))
    plt.boxplot(data, labels=samples, showfliers=False)
    plt.title(title)
    plt.ylabel(ylabel)
    plt.xticks(rotation=45, ha="right")
    plt.tight_layout()
    plt.savefig(outpath, dpi=250)
    plt.close()


def save_hexbin(x, y, outpath, title, xlabel, ylabel, xlog=False, ylog=False):
    """Save a hexbin plot for dense relationships like TLOD vs DP or TLOD vs AF."""
    mask = x.notna() & y.notna() & np.isfinite(x) & np.isfinite(y)
    x = x[mask]
    y = y[mask]

    if xlog:
        x = np.log10(x + 1)
        xlabel = f"log10({xlabel} + 1)"
    if ylog:
        y = np.log10(y + 1)
        ylabel = f"log10({ylabel} + 1)"

    if len(x) == 0:
        return

    plt.figure(figsize=(8, 6))
    hb = plt.hexbin(x, y, gridsize=80, bins="log", mincnt=1)
    plt.colorbar(hb, label="log10(variant count)")
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.tight_layout()
    plt.savefig(outpath, dpi=250)
    plt.close()


def save_overlay_density_by_status(status_df, column, outpath, title, xlabel, transform=None):
    """
    Overlay distributions for biologically interpretable groups:
    sample_only, shared_with_merged, and merged_only.

    This is usually more informative than making one plot per sample-vs-merged
    comparison because it asks: do retained/shared calls look higher quality than
    calls that are dropped or only appear after merging/force-calling?
    """
    plt.figure(figsize=(9, 6))
    plotted = 0

    for status in ["sample_only", "shared_with_merged", "merged_only"]:
        vals = status_df.loc[status_df["status"] == status, column].dropna()
        if transform == "log10_plus_1":
            vals = np.log10(vals + 1)
        vals = vals[np.isfinite(vals)]
        if len(vals) == 0:
            continue

        plt.hist(
            vals,
            bins=100,
            density=True,
            histtype="step",
            linewidth=2,
            label=f"{status} (n={len(vals):,})"
        )
        plotted += 1

    if plotted == 0:
        plt.close()
        return

    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel("Density")
    plt.legend(frameon=False)
    plt.tight_layout()
    plt.savefig(outpath, dpi=250)
    plt.close()


def save_retention_by_metric_bins(status_df, metric, outpath, title, xlabel, bins):
    """
    Plot the fraction of original sample calls retained/shared with the merged VCF
    across AF/DP/TLOD bins.

    This is often more useful than Jaccard because it asks whether stronger calls
    are preferentially retained.
    """
    original = status_df[status_df["status"].isin(["sample_only", "shared_with_merged"])].copy()
    original = original.dropna(subset=[metric])
    original = original[np.isfinite(original[metric])]
    if len(original) == 0:
        return

    original["bin"] = pd.cut(original[metric], bins=bins, include_lowest=True)
    rates = (
        original
        .groupby("bin", observed=True)["status"]
        .apply(lambda s: (s == "shared_with_merged").mean())
        .reset_index(name="retention_rate")
    )
    rates["bin_label"] = rates["bin"].astype(str)
    rates.to_csv(outpath.with_suffix(".csv"), index=False)

    plt.figure(figsize=(10, 6))
    plt.plot(rates["bin_label"], rates["retention_rate"], marker="o")
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel("Fraction retained/shared with merged VCF")
    plt.xticks(rotation=45, ha="right")
    plt.ylim(0, 1)
    plt.tight_layout()
    plt.savefig(outpath, dpi=250)
    plt.close()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--vcf_dir", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--merged_vcf", default=None)
    parser.add_argument("--pass_only", action="store_true")
    parser.add_argument(
        "--skip_jaccard_heatmap",
        action="store_true",
        help="Do not make the pairwise Jaccard heatmap. The CSV is still written."
    )
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

        if len(df) == 0:
            print(f"WARNING: {path} had no variants after filtering; skipping")
            continue

        all_dfs.append(df)
        sample = df["sample"].iloc[0]

        # Keep per-sample histograms, but the more useful plots are the cohort
        # overlays created after all samples are parsed.
        save_hist(df, "AF", plotdir / f"{sample}.AF_hist.png", f"{sample} AF Distribution", "Allele fraction", log_y=True)
        save_hist(df, "DP", plotdir / f"{sample}.DP_hist.png", f"{sample} DP Distribution", "Depth", log_y=True)
        save_hist(df, "TLOD", plotdir / f"{sample}.TLOD_hist.png", f"{sample} TLOD Distribution", "TLOD", log_y=True)

    if not all_dfs:
        raise ValueError("No variants found after parsing/filtering")

    all_variants = pd.concat(all_dfs, ignore_index=True)
    all_variants.to_csv(outdir / "all_sample_variants.long.csv", index=False)

    summary = pd.DataFrame([summarize_sample(df) for df in all_dfs])
    summary.to_csv(outdir / "per_sample_summary.csv", index=False)

    print("\nPer-sample summary:")
    print(summary)

    save_bar(summary.set_index("sample")["total_variants"], plotdir / "total_variants_per_sample.png", "Total Variants per Sample", "Variant count")
    save_bar(summary.set_index("sample")["pass_variants"], plotdir / "pass_variants_per_sample.png", "PASS Variants per Sample", "PASS variant count")
    save_bar(summary.set_index("sample")["median_TLOD"], plotdir / "median_TLOD_per_sample.png", "Median TLOD per Sample", "Median TLOD")

    # Cohort-level overlays: these put all samples on one graph, which is usually
    # much easier to interpret than one isolated plot per sample.
    save_overlay_density_by_sample(all_variants, "AF", plotdir / "overlay_AF_density_by_sample.png", "AF Distribution Overlay by Sample", "Allele fraction")
    save_overlay_density_by_sample(all_variants, "DP", plotdir / "overlay_log10DP_density_by_sample.png", "Depth Distribution Overlay by Sample", "log10(DP + 1)", transform="log10_plus_1")
    save_overlay_density_by_sample(all_variants, "TLOD", plotdir / "overlay_TLOD_density_by_sample.png", "TLOD Distribution Overlay by Sample", "TLOD")

    save_boxplot_by_sample(all_variants, "AF", plotdir / "AF_boxplot_by_sample.png", "AF by Sample", "Allele fraction")
    save_boxplot_by_sample(all_variants, "DP", plotdir / "log10DP_boxplot_by_sample.png", "Depth by Sample", "log10(DP + 1)", transform="log10_plus_1")
    save_boxplot_by_sample(all_variants, "TLOD", plotdir / "TLOD_boxplot_by_sample.png", "TLOD by Sample", "TLOD")

    save_hexbin(all_variants["AF"], all_variants["TLOD"], plotdir / "cohort_TLOD_vs_AF_hexbin.png", "Cohort TLOD vs AF", "AF", "TLOD")
    save_hexbin(all_variants["DP"], all_variants["TLOD"], plotdir / "cohort_TLOD_vs_log10DP_hexbin.png", "Cohort TLOD vs DP", "DP", "TLOD", xlog=True)

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
        .pivot_table(index="variant_id", columns="sample", values="present", aggfunc="max", fill_value=0)
        .astype(int)
    )
    presence.to_csv(outdir / "variant_presence_matrix.csv")
    presence["support_count"] = presence.sum(axis=1)

    support_hist = presence["support_count"].value_counts().sort_index().reset_index()
    support_hist.columns = ["n_samples_present", "variant_count"]
    support_hist.to_csv(outdir / "support_count_histogram.csv", index=False)
    save_bar(support_hist.set_index("n_samples_present")["variant_count"], plotdir / "support_count_histogram.png", "How Many Samples Each Variant Appears In", "Variant count")

    samples = [df["sample"].iloc[0] for df in all_dfs]
    sample_sets = {
        sample: set(all_variants.loc[all_variants["sample"] == sample, "variant_id"])
        for sample in samples
    }

    # Keep pairwise overlap counts as a CSV because they are useful for auditing,
    # but the heatmap is optional because Jaccard can be visually uninformative
    # when all values are uniformly low.
    jaccard = pd.DataFrame(index=samples, columns=samples, dtype=float)
    overlap_rows = []
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

    if not args.skip_jaccard_heatmap:
        plt.figure(figsize=(10, 8))
        plt.imshow(jaccard, aspect="auto")
        plt.colorbar(label="Jaccard index")
        plt.xticks(range(len(jaccard.columns)), jaccard.columns, rotation=90)
        plt.yticks(range(len(jaccard.index)), jaccard.index)
        plt.title("Pairwise Jaccard Similarity")
        plt.tight_layout()
        plt.savefig(plotdir / "pairwise_jaccard_heatmap.png", dpi=250)
        plt.close()

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
        save_bar(merged_compare.set_index("sample")["fraction_sample_retained_in_merged"], plotdir / "fraction_sample_retained_in_merged.png", "Fraction of Each Sample's Variants Retained in Merged VCF", "Fraction retained")

        original_support = presence[["support_count"]].copy()
        original_support["in_merged"] = original_support.index.isin(merged_set)
        retention_by_support = original_support.groupby("support_count")["in_merged"].agg(["sum", "count"]).reset_index()
        retention_by_support["retention_rate"] = retention_by_support["sum"] / retention_by_support["count"]
        retention_by_support.columns = ["n_original_samples_present", "n_retained_in_merged", "n_original_variants", "retention_rate"]
        retention_by_support.to_csv(outdir / "merged_retention_by_original_support_count.csv", index=False)
        save_bar(retention_by_support.set_index("n_original_samples_present")["retention_rate"], plotdir / "merged_retention_by_support_count.png", "Merged Retention Rate by Original Sample Support Count", "Retention rate")

        # Build one pooled status table instead of making only separate sample-vs-merged plots.
        # This creates interpretable groups across the whole cohort:
        #   sample_only = variants seen in original sample VCFs but not in merged
        #   shared_with_merged = variants seen in original sample VCFs and merged
        #   merged_only = variants seen in merged but not in any original sample VCF
        original_ids = set(presence.index)
        original_status = all_variants.copy()
        original_status["status"] = np.where(
            original_status["variant_id"].isin(merged_set),
            "shared_with_merged",
            "sample_only"
        )

        merged_only = merged_df[~merged_df["variant_id"].isin(original_ids)].copy()
        merged_only["status"] = "merged_only"

        status_df = pd.concat([original_status, merged_only], ignore_index=True, sort=False)
        status_df.to_csv(outdir / "variant_metrics_by_overlap_status.csv", index=False)

        status_summary = (
            status_df
            .groupby("status")
            .agg(
                n_variants=("variant_id", "count"),
                median_AF=("AF", "median"),
                median_DP=("DP", "median"),
                median_TLOD=("TLOD", "median"),
                mean_AF=("AF", "mean"),
                mean_DP=("DP", "mean"),
                mean_TLOD=("TLOD", "mean")
            )
            .reset_index()
        )
        status_summary.to_csv(outdir / "overlap_status_metric_summary.csv", index=False)

        save_overlay_density_by_status(status_df, "AF", plotdir / "overlay_AF_density_by_overlap_status.png", "AF Distribution by Overlap Status", "Allele fraction")
        save_overlay_density_by_status(status_df, "DP", plotdir / "overlay_log10DP_density_by_overlap_status.png", "Depth Distribution by Overlap Status", "log10(DP + 1)", transform="log10_plus_1")
        save_overlay_density_by_status(status_df, "TLOD", plotdir / "overlay_TLOD_density_by_overlap_status.png", "TLOD Distribution by Overlap Status", "TLOD")

        # These plots answer: are variants with stronger evidence more likely to
        # be retained in the merged VCF?
        save_retention_by_metric_bins(
            status_df,
            "AF",
            plotdir / "retention_rate_by_AF_bin.png",
            "Merged Retention Rate by AF Bin",
            "AF bin",
            bins=[0, 0.01, 0.03, 0.05, 0.10, 0.25, 0.50, 1.0]
        )
        save_retention_by_metric_bins(
            status_df,
            "TLOD",
            plotdir / "retention_rate_by_TLOD_bin.png",
            "Merged Retention Rate by TLOD Bin",
            "TLOD bin",
            bins=[-np.inf, 3, 6.3, 10, 20, 50, 100, np.inf]
        )
        save_retention_by_metric_bins(
            status_df,
            "DP",
            plotdir / "retention_rate_by_DP_bin.png",
            "Merged Retention Rate by DP Bin",
            "DP bin",
            bins=[0, 10, 20, 30, 50, 100, 200, 500, np.inf]
        )

    print(f"\nDone. Results written to: {outdir}")


if __name__ == "__main__":
    main()
