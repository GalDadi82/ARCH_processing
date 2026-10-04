"""Link each job's per-sample folders into All_Samples/{2_mapped,5_vcfs}.

Port of the "Create All_Samples/5_vcfs and 2_mapped folders" notebook cell. Writes a
manifest of what was linked so Snakemake has a concrete output to depend on.
"""

import os
import re
import sys

import pandas as pd
import yaml

analysis_dir = snakemake.params.analysis_dir
output_base_dir = snakemake.params.output_base_dir
manifest_path = snakemake.output.manifest

SUB_DIRS = ["2_mapped", "5_vcfs"]
JOB_RE = re.compile(r"^job\d+$")

with open(snakemake.log[0], "w") as log:

    def emit(msg):
        print(msg)
        log.write(msg + "\n")

    for sub_dir in SUB_DIRS:
        os.makedirs(os.path.join(output_base_dir, sub_dir), exist_ok=True)

    records = []
    job_dirs = sorted(d for d in os.listdir(analysis_dir) if JOB_RE.match(d))
    if not job_dirs:
        sys.exit(f"No job* directories found under {analysis_dir}")

    for job_num in job_dirs:
        sample_yaml_path = os.path.join(analysis_dir, job_num, "Sample_file_path_info.yaml")
        if not os.path.exists(sample_yaml_path):
            emit(f"Skipping {job_num}: no Sample_file_path_info.yaml")
            continue

        with open(sample_yaml_path) as f:
            sample_names = yaml.safe_load(f).keys()

        for sample_name in sample_names:
            for sub_dir in SUB_DIRS:
                src_path = os.path.join(analysis_dir, job_num, "paired", sub_dir, sample_name)
                dst_path = os.path.join(output_base_dir, sub_dir, sample_name)

                if not os.path.islink(dst_path):
                    emit(f"Linking {sample_name} from {src_path} to {dst_path}")
                    os.symlink(src_path, dst_path)

                records.append(
                    {
                        "Sample": sample_name,
                        "job": job_num,
                        "sub_dir": sub_dir,
                        "source": src_path,
                        "link": dst_path,
                        "source_exists": os.path.exists(src_path),
                    }
                )

    df = pd.DataFrame(records)
    df.to_csv(manifest_path, sep="\t", index=False)

    n_samples = df["Sample"].nunique() if not df.empty else 0
    emit(f"Linked {n_samples} samples from {len(job_dirs)} job directories")

    missing = df[~df["source_exists"]] if not df.empty else df
    if not missing.empty:
        emit(f"WARNING: {len(missing)} link target(s) do not exist, e.g.:")
        for src in missing["source"].head(5):
            emit(f"  {src}")
