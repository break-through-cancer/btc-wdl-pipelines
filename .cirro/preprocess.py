import json
import pandas as pd
from cirro.helpers.preprocess_dataset import PreprocessDataset


def extract_bams(ds):
    df = ds.files.copy()
    df["file"] = df["file"].astype(str)

    print("\n=== ALL FILES IN DATASET ===")
    print(ds.files.to_string(index=False))

    df = df[df["file"].str.endswith(".bam") | df["file"].str.endswith(".bam.bai")]

    print(f"\n=== BAM/BAI FILES ONLY ({len(df)} rows) ===")
    print(df.to_string(index=False))

    bam_map = {}
    for sample, group in df.groupby("sample"):
        bam = ""
        bai = ""
        for f in group["file"]:
            if f.endswith(".bam") and not f.endswith(".bam.bai"):
                bam = f
            elif f.endswith(".bam.bai"):
                bai = f

        print(f"\n  Sample '{sample}':")
        print(f"    BAM: {bam or '*** MISSING ***'}")
        print(f"    BAI: {bai or '*** MISSING ***'}")

        if bam:
            bam_map[str(sample)] = {"bam": bam, "bai": bai}
        else:
            print(f"    !! Skipping '{sample}' — no BAM found")

    if not bam_map:
        raise ValueError("No BAMs found in dataset")

    print(f"\n=== BAM MAP SUMMARY ({len(bam_map)} samples) ===")
    for s, f in bam_map.items():
        print(f"  {s}: {f}")

    return bam_map


def build_mutect_runs(samples: pd.DataFrame, bam_map: dict) -> list[dict]:
    runs = []

    samples = samples.copy()
    samples["status"] = samples["status"].astype(int)

    print("\n=== SAMPLESHEET USED FOR PAIRING ===")
    print(samples.to_string(index=False))

    print(f"\n=== PAIRING TUMOR + NORMAL (grouped by patient) ===")

    for patient, group in samples.groupby("patient"):
        tumor_rows  = group[group["status"] == 1]
        normal_rows = group[group["status"] == 0]

        print(f"\n  Patient '{patient}':")
        print(f"    Tumor rows:  {tumor_rows['sample'].tolist()}")
        print(f"    Normal rows: {normal_rows['sample'].tolist()}")

        normal_bam = normal_bai = normal_sample_name = None
        if not normal_rows.empty:
            n_sample = str(normal_rows.iloc[0]["sample"])
            if n_sample in bam_map:
                normal_bam         = bam_map[n_sample]["bam"]
                normal_bai         = bam_map[n_sample]["bai"]
                normal_sample_name = n_sample
                print(f"    Normal BAM:  {normal_bam}")
            else:
                print(f"    !! Normal sample '{n_sample}' not found in bam_map — tumor-only mode")
        else:
            print(f"    No normal rows — tumor-only mode")

        for _, trow in tumor_rows.iterrows():
            t_sample = str(trow["sample"])
            if t_sample not in bam_map:
                print(f"    !! Tumor sample '{t_sample}' not in bam_map — skipping")
                continue

            output_prefix = f"{patient}__{t_sample}"
            print(f"    -> Building run: {output_prefix}")

            run = {
                "output_prefix":     output_prefix,
                "tumor_reads":       bam_map[t_sample]["bam"],
                "tumor_reads_index": bam_map[t_sample]["bai"],
                "tumor_sample_name": t_sample,
            }

            if normal_bam:
                run["normal_reads"]        = normal_bam
                run["normal_reads_index"]  = normal_bai
                run["normal_sample_name"]  = normal_sample_name
                print(f"       Tumor:  {run['tumor_reads']}")
                print(f"       Normal: {run['normal_reads']}")
            else:
                print(f"       Tumor:  {run['tumor_reads']}")
                print(f"       Normal: (none — tumor-only)")

            runs.append(run)

    if not runs:
        raise ValueError("No valid tumor/normal pairs could be built.")
    return runs


def main():
    ds = PreprocessDataset.from_running()

    bam_map = extract_bams(ds)

    if "status" not in ds.samplesheet.columns:
        print("\n!! No 'status' column in samplesheet — falling back to PBMC name heuristic")
        rows = []
        for sample in bam_map:
            status = 0 if "PBMC" in sample.upper() else 1
            rows.append({"sample": sample, "patient": sample, "status": status})
            print(f"  '{sample}' -> status={status}")
        samplesheet = pd.DataFrame(rows)
    else:
        samplesheet = ds.samplesheet

    runs = build_mutect_runs(samplesheet, bam_map)

    print(f"\n=== FINAL RUNS ({len(runs)} total) ===")
    for i, r in enumerate(runs):
        print(f"\n  Run {i+1}: {r['output_prefix']}")
        print(f"    tumor_reads:       {r['tumor_reads']}")
        print(f"    tumor_reads_index: {r['tumor_reads_index']}")
        print(f"    tumor_sample_name: {r['tumor_sample_name']}")
        print(f"    normal_reads:      {r.get('normal_reads', '(none)')}")
        print(f"    normal_reads_index:{r.get('normal_reads_index', '(none)')}")

    # Sanity checks
    print("\n=== SANITY CHECKS ===")
    prefixes = [r["output_prefix"] for r in runs]
    if len(prefixes) != len(set(prefixes)):
        print("!! DUPLICATE output_prefixes detected — this will cause collisions!")
    else:
        print(f"  output_prefixes: all {len(prefixes)} are unique")

    missing_bai = [r["output_prefix"] for r in runs if not r.get("tumor_reads_index")]
    if missing_bai:
        print(f"  !! Runs missing tumor BAI: {missing_bai}")
    else:
        print(f"  tumor BAIs: all present")

    ds.add_param("mutect_runs", runs)

    print("\n=== FINAL PARAMS EMITTED ===")
    print(json.dumps(ds.params, indent=2, default=str))


if __name__ == "__main__":
    main()
# import json
# from cirro.helpers.preprocess_dataset import PreprocessDataset
# from cirro.api.models.s3_path import S3Path

# # def process_bams(ds):
# #     for base_name, group in ds.files.groupby("sample"):
       
# #         for f in group["file"]:
# #             if f.endswith(".bam") and not f.endswith(".bam.bai"):
# #                 if "PBMC" in base_name:
# #                     normal_bam = f
# #                 else:
# #                     tumor_bam = f
# #             elif f.endswith(".bam.bai"):
# #                 if "PBMC" in base_name:
# #                     normal_bai = f
# #                 else:
# #                     tumor_bai = f

# #         if tumor_bam and tumor_bai:
# #             print(tumor_bai, tumor_bam)
# #         if normal_bam and normal_bai:
# #             print(normal_bai, normal_bam)

# def main():
#     ds = PreprocessDataset.from_running()
#     print(ds.files)
#     # process_bams(ds)

# if __name__ == "__main__":
#     main()

