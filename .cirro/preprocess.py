import json
import pandas as pd
from cirro.helpers.preprocess_dataset import PreprocessDataset
import pysam

def extract_bams(ds):
    df = ds.files.copy()
    df["file"] = df["file"].astype(str)

    # Keep only BAM + BAI
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

    bam_path = ds.files[ds.files['file'].str.endswith('.bam')]['file'].iloc[0]

    with pysam.AlignmentFile(bam_path, 'rb', check_sq=False) as bam:
        rg_samples = list({rg['SM'] for rg in bam.header.to_dict().get('RG', [])})

    assert len(rg_samples) == 1, f"Expected 1 SM tag, got: {rg_samples}"
    tumor_sample = rg_samples[0]

    ds.add_param('tumor_sample', tumor_sample)


    bam_map = extract_bams(ds)

    tumor_bam = None
    tumor_bai = None
    normal_bam = None
    normal_bai = None

    # Simple rule: PBMC = normal, everything else = tumor
    for sample, files in bam_map.items():
        if "PBMC" in sample.upper(): # switch on the sample type column "Status"
            normal_bam = files["bam"]
            normal_bai = files["bai"]
        else:
            tumor_bam = files["bam"]
            tumor_bai = files["bai"]
            tumor_sample_name = str(sample)

    if not tumor_bam:
        raise ValueError("No tumor BAM found")

    # Always set tumor
    ds.add_param("tumor_reads", tumor_bam)
    ds.add_param("tumor_reads_index", tumor_bai)
    ds.add_param("tumor_sample_name", tumor_sample_name)

    # Only set normal if present
    if normal_bam:
        ds.add_param("normal_reads", normal_bam)
        ds.add_param("normal_reads_index", normal_bai)
        print("Matched normal detected.")
    else:
        print("No normal detected. Tumor-only mode.")

    print("\nFinal parameters:")
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

