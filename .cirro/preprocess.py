#!/usr/bin/env python3
import json
from cirro.helpers.preprocess_dataset import PreprocessDataset

PROCESS_INPUT_FILE = "process-input.json"
WORKFLOW_PREFIX = "Mutect2"

def main():
    ds = PreprocessDataset.from_running()

    # Load process-input.json
    with open(PROCESS_INPUT_FILE) as f:
        process_inputs = json.load(f)

    all_inputs = []

    df = ds.files
    for sample_name, group in df.groupby("sample"):
        tumor_bam = None
        tumor_index = None
        normal_bam = None
        normal_index = None

        for f in group["file"]:
            if "tumor" in f.lower():
                if f.endswith(".bam"):
                    tumor_bam = f
                elif f.endswith(".bai"):
                    tumor_index = f
            elif "normal" in f.lower():
                if f.endswith(".bam"):
                    normal_bam = f
                elif f.endswith(".bai"):
                    normal_index = f

        if not tumor_bam or not tumor_index:
            raise ValueError(f"Missing tumor BAM or index for sample {sample_name}")

        # combine workflow-level inputs with sample-specific BAM paths
        input_dict = process_inputs.copy()
        input_dict.update({
            f"{WORKFLOW_PREFIX}.tumor_reads": tumor_bam,
            f"{WORKFLOW_PREFIX}.tumor_reads_index": tumor_index,
        })
        if normal_bam and normal_index:
            input_dict[f"{WORKFLOW_PREFIX}.normal_reads"] = normal_bam
            input_dict[f"{WORKFLOW_PREFIX}.normal_reads_index"] = normal_index

        all_inputs.append(input_dict)

    # Write out a JSON for each run
    for i, input_dict in enumerate(all_inputs):
        with open(f"inputs.{i}.json", "w") as f:
            json.dump(input_dict, f, indent=4)

    print(f"Generated {len(all_inputs)} input JSON(s) for {WORKFLOW_PREFIX}.")

if __name__ == "__main__":
    main()
