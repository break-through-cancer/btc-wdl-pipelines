#!/usr/bin/env python3
import json
from cirro.helpers.preprocess_dataset import PreprocessDataset

PROCESS_INPUT_FILE = "process-input.json"
WORKFLOW_PREFIX = "Mutect2"

def generate_mutect2_inputs(ds: PreprocessDataset, process_inputs: dict):
    """
    ds: PreprocessDataset object from Cirro
    process_inputs: dictionary loaded from process-input.json
    """
    all_inputs = []

    # user-provided tumor BAM
    user_tumor_bam = process_inputs.get("tumor_reads")
    user_tumor_index = process_inputs.get("tumor_reads_index")

    # Group files by sample
    df = ds.files
    for sample_name, group in df.groupby("sample"):
        dataset_tumor_bam = None
        dataset_tumor_index = None
        dataset_normal_bam = None
        dataset_normal_index = None

        # Scan all files in this sample
        for f in group["file"]:
            f_path = getattr(f, "path", str(f))  # get string path
            fname = f_path.lower()
            if dataset_tumor_bam is None and fname.endswith(".bam") and "tumor" in fname:
                dataset_tumor_bam = f_path
            elif dataset_tumor_index is None and fname.endswith(".bai") and "tumor" in fname:
                dataset_tumor_index = f_path
            elif dataset_normal_bam is None and fname.endswith(".bam") and "normal" in fname:
                dataset_normal_bam = f_path
            elif dataset_normal_index is None and fname.endswith(".bai") and "normal" in fname:
                dataset_normal_index = f_path

        # Prefer user-provided tumor BAMs
        final_tumor_bam = user_tumor_bam or dataset_tumor_bam
        final_tumor_index = user_tumor_index or dataset_tumor_index
        final_normal_bam = dataset_normal_bam
        final_normal_index = dataset_normal_index

        if not final_tumor_bam or not final_tumor_index:
            raise ValueError(f"Missing tumor BAM or index for sample {sample_name}")

        # combine workflow-level inputs with sample-specific BAM paths
        input_dict = process_inputs.copy()
        input_dict[f"{WORKFLOW_PREFIX}.tumor_reads"] = final_tumor_bam
        input_dict[f"{WORKFLOW_PREFIX}.tumor_reads_index"] = final_tumor_index
        if final_normal_bam and final_normal_index:
            input_dict[f"{WORKFLOW_PREFIX}.normal_reads"] = final_normal_bam
            input_dict[f"{WORKFLOW_PREFIX}.normal_reads_index"] = final_normal_index

        all_inputs.append(input_dict)

    return all_inputs


def main():
    # Load running dataset via Cirro
    ds = PreprocessDataset.from_running()

    # Load process-input.json
    with open(PROCESS_INPUT_FILE) as f:
        process_inputs = json.load(f)

    # Generate all input dictionaries
    all_inputs = generate_mutect2_inputs(ds, process_inputs)

    # Write JSON files for each sample/run
    for i, input_dict in enumerate(all_inputs):
        with open(f"inputs.{i}.json", "w") as f:
            json.dump(input_dict, f, indent=4)

    print(f"Generated {len(all_inputs)} input JSON(s) for {WORKFLOW_PREFIX}.")


if __name__ == "__main__":
    main()