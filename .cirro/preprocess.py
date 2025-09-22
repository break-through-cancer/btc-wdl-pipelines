#!/usr/bin/env python3
import json
from cirro.helpers.preprocess_dataset import PreprocessDataset

PROCESS_INPUT_FILE = "process-input.json"
WORKFLOW_PREFIX = "Mutect2"

def yield_single_inputs(ds: PreprocessDataset):
    df = ds.files  # ds.files is a DataFrame with columns ["sample", "file"]
    for base_name, group in df.groupby("sample"):
        if "PBMC" in base_name:
                continue

        sample_to_analyze = None
        sample_to_analyze_index = None
        normal_bam = None
        normal_bai = None

        for f in group["file"]:
            if f.endswith(".bam") and not f.endswith(".bam.bai"):
                if "PBMC" in base_name:
                    normal_bam = f
                else:
                    sample_to_analyze = f
            elif f.endswith(".bam.bai"):
                if "PBMC" in base_name:
                    normal_bai = f
                else:
                    sample_to_analyze_index = f

        if sample_to_analyze and sample_to_analyze_index:
            yield {
                f"{WORKFLOW_PREFIX}.input_bam": sample_to_analyze,
                f"{WORKFLOW_PREFIX}.input_bam_index": sample_to_analyze_index
            }
        
        if normal_bam and normal_bai:
            yield {
                f"{WORKFLOW_PREFIX}.normal_bam": normal_bam,
                f"{WORKFLOW_PREFIX}.normal_bam_index": normal_bai
            }


def setup_inputs(ds: PreprocessDataset):
    # Make a combined set of inputs with each of the BAM files
    all_inputs = [
        {
            **single_input,
            **{
                kw: val
                for kw, val in ds.params.items()
                if kw.startswith(WORKFLOW_PREFIX)
            }
        }
        for single_input in yield_single_inputs(ds)
    ]

    # Raise an error if no inputs are found
    assert len(all_inputs) > 0, "No inputs found -- stopping execution"

    # complete set of inputs
    write_json("inputs.json", all_inputs)

    #individual
    for i, input in enumerate(all_inputs):
        write_json(f"inputs.{i}.json", input)

def write_json(fp, obj, indent=4) -> None:

    with open(fp, "wt") as handle:
        json.dump(obj, handle, indent=indent)


def main():
    """Primary entrypoint for the script"""

    # Get information on the analysis launched by the user
    ds = PreprocessDataset.from_running()

    # # Set up the options.json file
    # setup_options(ds)

    setup_inputs(ds)

if __name__ == "__main__":
    main()