import json
from cirro.helpers.preprocess_dataset import PreprocessDataset
from cirro.api.models.s3_path import S3Path

def process_bams(ds):
    for base_name, group in ds.files.groupby("sample"):
       
        for f in group["file"]:
            if f.endswith(".bam") and not f.endswith(".bam.bai"):
                if "PBMC" in base_name:
                    normal_bam = f
                else:
                    tumor_bam = f
            elif f.endswith(".bam.bai"):
                if "PBMC" in base_name:
                    normal_bai = f
                else:
                    tumor_bai = f

        if tumor_bam and tumor_bai:
            print(tumor_bai, tumor_bam)
        if normal_bam and normal_bai:
            print(normal_bai, normal_bam)

def main():
    ds = PreprocessDataset.from_running()
    print(ds.files)
    print(ds.samples)
    print(ds.params)
    process_bams(ds)

if __name__ == "__main__":
    main()

