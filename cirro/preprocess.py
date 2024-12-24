#!/usr/bin/env python3

# Import the packages needed to run the code below
import json
from cirro.helpers.preprocess_dataset import PreprocessDataset
from cirro.api.models.s3_path import S3Path


def setup_options_inputs(ds: PreprocessDataset):

    # Set up the scriptBucketName, which is needed by the workflow
    # to stage analysis scripts
    ds.add_param(
        "scriptBucketName",
        S3Path(ds.params['final_workflow_outputs_dir']).bucket
    )

    # Isolate the options arguments for the workflow
    # Define a new dictionary which contains all of the items
    # from `ds.params` which do not start with the workflow
    # prefix "HapCNA"
    options = {
        kw: val
        for kw, val in ds.params.items()
        if not kw.startswith(("HapCNA", "CheckSamplesUnique"))
    }

    inputs_1 = {
        kw: val
        for kw, val in ds.params.items()
        if kw.startswith("HapCNA")
    }

    # inputs = {
    #     kw: val
    #     for kw, val in ds.params.items()
    #     if kw.startswith("CheckSamplesUnique")
    # }

    # Write out to the options.json file
    write_json("options.json", options)
    # write_json("inputs.json", inputs)
    write_json("inputs.1.json", inputs_1)


def yield_single_inputs(ds: PreprocessDataset) -> dict:
    """
    This function is used to identify each of the BAM/BAI pairs
    which have been provided as part of the input dataset by the user.
    The output of this function will be the a yielded series of
    struct objects with the `input_bam` and `input_bam_index` attributes.
    """

    # The ds.files object is a DataFrame with columns for `sample` and `file`
    # For example:
    # | sample  | file            |
    # |---------|-----------------|
    # | sampleA | sampleA.bam     |
    # | sampleA | sampleA.bam.bai |
    # | sampleB | sampleB.bam     |
    # | sampleB | sampleB.bam.bai |

    # Iterate over each batch of files with a distinct sample
    for sample, files in ds.files.groupby("sample"):

        # Set up an empty struct
        dat = dict(input_bam=None, input_bam_index=None)

        # Iterate over each file
        for file in files['file'].values:
            for suffix, kw in [(".bam", "input_bam"), (".bai", "input_bam_index")]:
                if file.endswith(suffix):
                    if dat[kw] is not None:
                        raise ValueError(f"Multiple '{suffix}' files found for sample {sample}")
                    dat[kw] = file

        ds.logger.info(f"Sample: {sample}")
        if dat["input_bam"] is None:
            ds.logger.info("No BAM file found, skipping")
            continue
        if dat["input_bam_index"] is None:
            ds.logger.info("No BAM Index file found, skipping")
            continue

        ds.logger.info(f"BAM: {dat['input_bam']}")
        ds.logger.info(f"BAM Index: {dat['input_bam_index']}")

        # Add the workflow prefix
        yield {
            f"HapCNA.{kw}": val
            for kw, val in dat.items()
        }


def setup_inputs(ds: PreprocessDataset):

    # Make a combined set of inputs with each of the BAM files
    all_inputs = [
        {
            **single_input,
            **{
                kw: val
                for kw, val in ds.params.items()
                if kw.startswith("HapCNA")
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
    setup_options_inputs(ds)
    # Set up the inputs files
    # setup_inputs(ds)


if __name__ == "__main__":
    main()