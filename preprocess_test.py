import json
import pandas as pd

WORKFLOW_PREFIX = "Mutect2"

import json

WORKFLOW_PREFIX = "Mutect2"

def generate_mutect2_inputs(files, process_inputs):
    """
    files: list of file-like objects (e.g., from portal.list_files())
    process_inputs: dict loaded from process-input.json
    """
    all_inputs = []

    # User-provided tumor BAM (can be None)
    user_tumor_bam = process_inputs.get(WORKFLOW_PREFIX+"tumor_reads")

    # Scan dataset files for a BAM if user didn't provide one
    dataset_tumor_bam = None
    dataset_tumor_index = None
    for f in files:
        f_path = getattr(f, "path", str(f))  # get path string
        if dataset_tumor_bam is None and f_path.endswith(".bam"):
            dataset_tumor_bam = f_path
        elif dataset_tumor_index is None and f_path.endswith(".bai"):
            dataset_tumor_index = f_path

    final_tumor_bam = user_tumor_bam or dataset_tumor_bam
    final_tumor_index = dataset_tumor_index  # still use dataset index

    input_dict = process_inputs.copy()
    input_dict[f"{WORKFLOW_PREFIX}.tumor_reads"] = final_tumor_bam
    input_dict[f"{WORKFLOW_PREFIX}.tumor_reads_index"] = final_tumor_index

    all_inputs.append(input_dict)
    return all_inputs
