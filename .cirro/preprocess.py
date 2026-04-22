import json
import pandas as pd
from cirro.helpers.preprocess_dataset import PreprocessDataset


def extract_bams(ds):
    df = ds.files.copy()
    df["file"] = df["file"].astype(str)

    df = df[df["file"].str.endswith(".bam") | df["file"].str.endswith(".bam.bai")]

    bam_map = {}

    for sample, group in df.groupby("sample"):
        bam = ""
        bai = ""

        for f in group["file"]:
            if f.endswith(".bam") and not f.endswith(".bam.bai"):
                bam = f
            elif f.endswith(".bam.bai"):
                bai = f

        if bam:
            bam_map[str(sample)] = {"bam": bam, "bai": bai}

    if not bam_map:
        raise ValueError("No BAMs found in dataset")

    return bam_map


def main():
    ds = PreprocessDataset.from_running()

    print("=== ds.files preview ===")
    print(ds.files.head(20).to_string(index=False))

    bam_map = extract_bams(ds)

    normal_bam = None
    normal_bai = None
    tumor_runs = []

    for sample, files in bam_map.items():
        sample_name = str(sample)

        if "PBMC" in sample_name.upper():
            normal_bam = files["bam"]
            normal_bai = files["bai"]
        else:
            tumor_runs.append({
                "output_prefix": sample_name,
                "tumor_reads": files["bam"],
                "tumor_reads_index": files["bai"],
                "tumor_sample_name": sample_name,
            })

    if not tumor_runs:
        raise ValueError("No tumor BAM found")

    mutect_runs = []
    for run in tumor_runs:
        mutect_runs.append({
            "output_prefix": run["output_prefix"],
            "tumor_reads": run["tumor_reads"],
            "tumor_reads_index": run["tumor_reads_index"],
            "normal_reads": normal_bam,
            "normal_reads_index": normal_bai,
            "tumor_sample_name": run["tumor_sample_name"],
        })

    ds.add_param("mutect_runs", mutect_runs)

    print("\nFinal parameters:")
    print(json.dumps(ds.params, indent=2, default=str))


if __name__ == "__main__":
    main()