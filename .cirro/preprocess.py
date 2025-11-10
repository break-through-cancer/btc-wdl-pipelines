import json
from cirro.helpers.preprocess_dataset import PreprocessDataset
from cirro.api.models.s3_path import S3Path

PROCESS_INPUT_FILE = "process-input.json"
WORKFLOW_PREFIX = "Mutect2"

def write_json(fp, obj, indent=4):
    with open(fp, "wt") as handle:
        json.dump(obj, handle, indent=indent)

def collapse_arrays(obj):
    if isinstance(obj, list):
        return [collapse_arrays(o) for o in obj]
    elif isinstance(obj, dict):
        new = {}
        for k, v in obj.items():
            if isinstance(v, list) and len(v) == 1 and not k.endswith("tumor_reads") and not k.endswith("normal_reads"):
                new[k] = v[0]
            else:
                new[k] = v
        return new
    else:
        return obj

def yield_single_inputs(ds: PreprocessDataset):
    df = ds.files
    for base_name, group in df.groupby("sample"):
        if "PBMC" in base_name:
            continue
        tumor_bam = None
        tumor_bai = None
        normal_bam = None
        normal_bai = None
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
            yield {
                f"{WORKFLOW_PREFIX}.tumor_reads": tumor_bam,
                f"{WORKFLOW_PREFIX}.tumor_reads_index": tumor_bai,
            }
        if normal_bam and normal_bai:
            yield {
                f"{WORKFLOW_PREFIX}.normal_reads": normal_bam,
                f"{WORKFLOW_PREFIX}.normal_reads_index": normal_bai,
            }

def yield_joint_inputs(ds: PreprocessDataset):
    df = ds.files
    for base_name, group in df.groupby("sample"):
        if "PBMC" in base_name:
            continue
        tumor_bams = []
        tumor_bais = []
        normal_bams = []
        normal_bais = []

        for f in group["file"]:
            if f.endswith(".bam") and not f.endswith(".bam.bai"):
                if "PBMC" in f or "normal" in f.lower():
                    normal_bams.append(f)
                else:
                    tumor_bams.append(f)
            elif f.endswith(".bam.bai"):
                if "PBMC" in f or "normal" in f.lower():
                    normal_bais.append(f)
                else:
                    tumor_bais.append(f)

        inputs = {}
        if tumor_bams and tumor_bais:
            inputs[f"{WORKFLOW_PREFIX}.tumor_reads"] = tumor_bams
            inputs[f"{WORKFLOW_PREFIX}.tumor_reads_index"] = tumor_bais
        if normal_bams and normal_bais:
            inputs[f"{WORKFLOW_PREFIX}.normal_reads"] = normal_bams
            inputs[f"{WORKFLOW_PREFIX}.normal_reads_index"] = normal_bais

        if inputs:
            yield inputs

def setup_inputs(ds: PreprocessDataset):
    # Load the process-input.json file if it exists
    try:
        with open(PROCESS_INPUT_FILE) as f:
            process_inputs = json.load(f)
    except FileNotFoundError:
        process_inputs = {}

    all_inputs = []

    # Decide which generator to use based on joint_calling
    if getattr(ds.params, "joint_calling", False):
        input_generator = yield_single_inputs(ds)
        print("Using single inputs")
    else:
        input_generator = yield_joint_inputs(ds)
        print("Using joint inputs")

    # Iterate over the chosen generator
    for input_dict in input_generator:
        # Merge everything together:
        # 1. fields from process-input.json
        # 2. workflow-specific parameters from ds.params
        # 3. tumor/normal files
        combined = {
            **process_inputs,
            **{k: v for k, v in ds.params.items() if k.startswith(WORKFLOW_PREFIX)},
            **input_dict
        }
        all_inputs.append(combined)

    assert all_inputs, "No inputs found -- stopping execution"
    all_inputs = collapse_arrays(all_inputs)

    # Write individual JSON files
    for i, inp in enumerate(all_inputs):
        write_json(f"inputs.{i}.json", inp)

    print("Inputs written:", all_inputs)
    return all_inputs

def setup_options(ds: PreprocessDataset):
    if "out_dir" in ds.params:
        ds.add_param("scriptBucketName", S3Path(ds.params['out_dir']).bucket)
    else:
        ds.add_param("scriptBucketName", "mock-bucket")

    options = {k: v for k, v in ds.params.items() if not k.startswith(WORKFLOW_PREFIX)}
    options = collapse_arrays(options)
    write_json("options.json", options)
    print("Options written:", options)
    return options

def main():
    ds = PreprocessDataset.from_running()
    setup_options(ds)
    setup_inputs(ds)

if __name__ == "__main__":
    main()
