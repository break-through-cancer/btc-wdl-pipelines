import cmd
import json
import pandas as pd
from cirro.helpers.preprocess_dataset import PreprocessDataset
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

    # cmd = f"samtools view -H {bam_path}"
    # header = subprocess.check_output(cmd, shell=True, text=True)

    # samples = set()
    # for line in header.splitlines():
    #     if line.startswith("@RG"):
    #         for field in line.split("\t"):
    #             if field.startswith("SM:"):
    #                 samples.add(field.replace("SM:", ""))

    # if len(samples) != 1:
    #     raise ValueError(f"Expected 1 SM tag, got: {samples}")

    # tumor_sample = list(samples)[0]

    # ds.add_param('tumor_sample', tumor_sample)


    bam_map = extract_bams(ds)

    tumor_bam = None
    tumor_bai = None
    normal_bam = None
    normal_bai = None
    tumor_sample_name = None

    for sample, files in bam_map.items():
        if "PBMC" in sample.upper():
            normal_bam = files["bam"]
            normal_bai = files["bai"]
        else:
            tumor_bam = files["bam"]
            tumor_bai = files["bai"]
            tumor_sample_name = str(sample)

    if not tumor_bam:
        raise ValueError("No tumor BAM found")

    mutect_runs = [
        {
            "output_prefix": tumor_sample_name,
            "tumor_reads": tumor_bam,
            "tumor_reads_index": tumor_bai,
            "normal_reads": normal_bam,
            "normal_reads_index": normal_bai,
            "tumor_sample_name": tumor_sample_name,
        }
    ]

    ds.add_param("mutect_runs", mutect_runs)

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

