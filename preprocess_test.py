#!/usr/bin/env python3
import json
from typing import List

WORKFLOW_PREFIX = "Mutect2"
PROCESS_INPUT_FILE = "process-input.json"

def generate_mutect2_inputs(files: List, process_inputs: dict):
    """
    files: list of file objects or paths. Each file should have a 'path' attribute or be a string.
    process_inputs: dictionary loaded from process-input.json
    """
    all_inputs = []

    # user-provided tumor BAM
    user_tumor_bam = process_inputs.get(WORKFLOW_PREFIX + ".tumor_reads")
    user_tumor_index = process_inputs.get(WORKFLOW_PREFIX + ".tumor_reads_index")

    dataset_tumor_bam = None
    dataset_tumor_index = None
    dataset_normal_bam = None
    dataset_normal_index = None

    # scan dataset for BAM/BAI files
    for f in files:
        f_path = getattr(f, "path", str(f))
        fname = f_path.lower()
        if dataset_tumor_bam is None and fname.endswith(".bam") and "tumor" in fname:
            dataset_tumor_bam = f_path
        elif dataset_tumor_index is None and fname.endswith(".bai") and "tumor" in fname:
            dataset_tumor_index = f_path
        elif dataset_normal_bam is None and fname.endswith(".bam") and "normal" in fname:
            dataset_normal_bam = f_path
        elif dataset_normal_index is None and fname.endswith(".bai") and "normal" in fname:
            dataset_normal_index = f_path

    # final BAMs: prefer user input if available
    final_tumor_bam = user_tumor_bam or dataset_tumor_bam
    final_tumor_index = user_tumor_index or dataset_tumor_index
    final_normal_bam = dataset_normal_bam
    final_normal_index = dataset_normal_index

    if final_tumor_bam is None or final_tumor_index is None:
        raise ValueError("Missing tumor BAM or index (neither user input nor dataset has it).")

    input_dict = process_inputs.copy()
    input_dict[f"{WORKFLOW_PREFIX}.tumor_reads"] = final_tumor_bam
    input_dict[f"{WORKFLOW_PREFIX}.tumor_reads_index"] = final_tumor_index

    if final_normal_bam and final_normal_index:
        input_dict[f"{WORKFLOW_PREFIX}.normal_reads"] = final_normal_bam
        input_dict[f"{WORKFLOW_PREFIX}.normal_reads_index"] = final_normal_index
        print('got to normals')

    all_inputs.append(input_dict)
    return all_inputs


# ------------------------------
# Example usage in a notebook
# ------------------------------
if __name__ == "__main__":
    # Load process-input.json
    with open(".cirro/process-input.json") as f:
        process_inputs = json.load(f)

    # Assume `files` is a list of objects returned from your dataset.list_files()
    # Example: files = dataset.list_files()
    all_inputs = generate_mutect2_inputs(files, process_inputs)

    # Write out JSONs
    for i, input_dict in enumerate(all_inputs):
        with open(f"inputs.{i}.json", "w") as f:
            json.dump(input_dict, f, indent=4)

    print(f"Generated {len(all_inputs)} input JSON(s) for {WORKFLOW_PREFIX}.")
