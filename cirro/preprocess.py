#!/usr/bin/env python3

# Import the packages needed to run the code below
import json
from cirro.helpers.preprocess_dataset import PreprocessDataset
from cirro.api.models.s3_path import S3Path

WORKFLOW_PREFIX = "HapCNA"

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

def setup_inputs_single(ds: PreprocessDataset):

    # Isolate the options arguments for the workflow
    # Define a new dictionary which contains all of the items
    # from `ds.params` which do not start with the workflow
    # prefix "HapCNA"

    inputs = {
        kw: val
        for kw, val in ds.params.items()
        if kw.startswith("HapCNA")
    }

    # Write out to the options.json file
    write_json("inputs.1.json", inputs)

def yield_single_inputs(ds: PreprocessDataset) -> dict:
    """
    This function processes the samplesheet and identifies each of the BAM/BAI pairs
    for tumor and normal samples, and assigns the corresponding files to the appropriate
    fields based on the input dataset.
    """
    
    for patient, group in ds.samplesheet.groupby("patient"):
        # Separate normal and tumor samples
        normal_samples = group[group["status"] == 0]
        tumor_samples = group[group["status"] == 1]

        # Skip patients with no normal or no tumor samples
        if normal_samples.empty or tumor_samples.empty:
            ds.logger.info(f"Patient {patient} has insufficient samples (Normal: {len(normal_samples)}, Tumor: {len(tumor_samples)}); skipping.")
            continue

        # Process each normal sample
        for _, normal_sample in normal_samples.iterrows():
            normal_sample_id = normal_sample["sample"]
            normal_sex = normal_sample["sex"]

            # Find the corresponding BAM/BAI files for the normal sample
            normal_files = ds.files[ds.files["sample"] == normal_sample_id]
            normal_dat = {
                "normal_sample_id": normal_sample_id,
                "normal_input_bam": None,
                "normal_input_bam_index": None,
            }

            for file in normal_files["file"].values:
                for suffix, key in [(".bam", "normal_input_bam"), (".bai", "normal_input_bam_index")]:
                    if file.endswith(suffix):
                        if normal_dat[key] is not None:
                            raise ValueError(f"Multiple '{suffix}' files found for normal sample {normal_sample_id}")
                        normal_dat[key] = file

            # Ensure both BAM and BAM index are found for the normal sample
            if normal_dat["normal_input_bam"] is None or normal_dat["normal_input_bam_index"] is None:
                ds.logger.info(f"Normal sample {normal_sample_id} is missing BAM or BAM index files; skipping.")
                continue

            # Process each tumor sample for the patient
            for _, tumor_sample in tumor_samples.iterrows():
                tumor_sample_id = tumor_sample["sample"]
                tumor_sex = tumor_sample["sex"]

                # Find the corresponding BAM/BAI files for the tumor sample
                tumor_files = ds.files[ds.files["sample"] == tumor_sample_id]
                tumor_dat = {
                    "tumor_sample_id": tumor_sample_id,
                    "tumor_input_bam": None,
                    "tumor_input_bam_index": None,
                }
                for file in tumor_files["file"].values:
                    for suffix, key in [(".bam", "tumor_input_bam"), (".bai", "tumor_input_bam_index")]:
                        if file.endswith(suffix):
                            if tumor_dat[key] is not None:
                                raise ValueError(f"Multiple '{suffix}' files found for tumor sample {tumor_sample_id}")
                            tumor_dat[key] = file

                # Ensure both BAM and BAM index are found for the tumor sample
                if tumor_dat["tumor_input_bam"] is None or tumor_dat["tumor_input_bam_index"] is None:
                    ds.logger.info(f"Tumor sample {tumor_sample_id} is missing BAM or BAM index files; skipping.")
                    continue

                # Validate that the sex matches between the normal and tumor samples
                if normal_sex != tumor_sex:
                    ds.logger.info(f"Sex mismatch for patient {patient}: Normal ({normal_sex}) vs Tumor ({tumor_sex}); skipping.")
                    continue

                if normal_sex.strip().upper() in ['F', 'XX']:
                    sample_sex = 'female'
                elif normal_sex.strip().upper() in ['M', 'XY']:
                    sample_sex = 'male'
                else:
                    raise ValueError(f"sample_sex value '{normal_sex}' for patient {patient} is invalid.")

                # Create the participant_id
                participant_id = f"{normal_sample_id}_{tumor_sample_id}"

                # Yield the data as a dictionary with workflow_prefix
                yield {
                    f"{WORKFLOW_PREFIX}.{key}": value
                    for key, value in {
                        "participant_id": participant_id,
                        "sample_sex": sample_sex,
                        **normal_dat,
                        **tumor_dat,
                    }.items()
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

    # Set up the inputs files
    input_type = ds.params.get("input_type")
    ds.remove_param("input_type")

    if input_type == "single":
        setup_inputs_single(ds)
    else:
        form_params = [
            "HapCNA.participant_id",
            "HapCNA.normal_input_bam",
            "HapCNA.normal_input_bam_index",
            "HapCNA.normal_sample_id",
            "HapCNA.tumor_input_bam",
            "HapCNA.tumor_input_bam_index",
            "HapCNA.tumor_sample_id",
            "HapCNA.sample_sex"
        ]
        
        for param in form_params:
            ds.remove_param(param, force=True)

        setup_inputs(ds)


if __name__ == "__main__":
    main()