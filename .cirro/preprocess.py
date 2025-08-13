#!/usr/bin/env python3

import json
from cirro.helpers.preprocess_dataset import PreprocessDataset
from cirro.api.models.s3_path import S3Path

WORKFLOW_PREFIX = "Funcotator"

def setup_options(ds: PreprocessDataset):

    # Set up the scriptBucketName, which is needed by the workflow
    # to stage analysis scripts
    ds.add_param(
        "scriptBucketName",
        S3Path(ds.params['final_workflow_outputs_dir']).bucket
    )

    # Isolate the options arguments for the workflow
    # Define a new dictionary which contains all of the items
    # from `ds.params` which do not start with the workflow
    # prefix
    options = {
        kw: val
        for kw, val in ds.params.items()
        if not kw.startswith(WORKFLOW_PREFIX)
    }

    # Write out to the options.json file
    write_json("options.json", options)

def yield_single_inputs(ds: PreprocessDataset):
    df = ds.files  # assuming ds.files is a DataFrame with columns ["sample", "file"]

    # Group by sample (participant_id)
    for base_name, group in df.groupby("sample"):
        # base_name should be in the format {tumor}_vs_{normal}
        ids = base_name.split('_vs_')
        case_id,control_id = ids[0],ids[1]

        variant_vcf_to_funcotate = None
        variant_vcf_to_funcotate_index = None

        for f in group["file"]:
            if f.endswith(".mutect2.filtered.vcf.gz"):
                variant_vcf_to_funcotate = f
            elif f.endswith(".mutect2.filtered.vcf.gz.tbi"):
                variant_vcf_to_funcotate_index = f

        if variant_vcf_to_funcotate and variant_vcf_to_funcotate_index:
            yield {
                f"{WORKFLOW_PREFIX}.Funcotate.case_id": case_id,
                f"{WORKFLOW_PREFIX}.Funcotate.control_id": control_id,
                f"{WORKFLOW_PREFIX}.output_file_base_name": base_name,
                f"{WORKFLOW_PREFIX}.variant_vcf_to_funcotate": variant_vcf_to_funcotate,
                f"{WORKFLOW_PREFIX}.variant_vcf_to_funcotate_index": variant_vcf_to_funcotate_index,
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

    # Write out the complete set of inputs
    write_json("inputs.json", all_inputs)

    # Write out each individual file pair
    for i, input in enumerate(all_inputs):
        write_json(f"inputs.{i}.json", input)

def write_json(fp, obj, indent=4) -> None:

    with open(fp, "wt") as handle:
        json.dump(obj, handle, indent=indent)


def main():
    """Primary entrypoint for the script"""

    # Get information on the analysis launched by the user
    ds = PreprocessDataset.from_running()

    # Set up the options.json file
    setup_options(ds)

    setup_inputs(ds)

if __name__ == "__main__":
    main()