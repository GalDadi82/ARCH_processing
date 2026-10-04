"""Report per-sample ANNOVAR multianno file sizes.

Port of the notebook's file-size QC cell -- a quick way to spot samples whose ANNOVAR run
produced nothing.
"""

import os

import pandas as pd

with open(snakemake.log[0], "w") as log:

    def emit(msg):
        print(msg)
        log.write(msg + "\n")

    records = []
    for file_path in snakemake.input.txt:
        # .../5_vcfs/{sample}/varscan_annot_dedup_filt.hg38_multianno.txt
        sample_name = os.path.basename(os.path.dirname(file_path))
        exists = os.path.exists(file_path)
        records.append(
            {
                "Sample_Name": sample_name,
                "File_Name": os.path.basename(file_path) if exists else None,
                "File_Path": file_path if exists else None,
                "File_Size_MB": round(os.path.getsize(file_path) / (1024 * 1024), 2)
                if exists
                else None,
            }
        )

    df = pd.DataFrame(records).sort_values("Sample_Name")
    df.to_excel(snakemake.output.report, index=False)

    emit(f"Report saved to: {snakemake.output.report}")
    empty = df[df["File_Size_MB"].fillna(0) == 0]
    if not empty.empty:
        emit(f"WARNING: {len(empty)} sample(s) have an empty or missing multianno table:")
        for sample in empty["Sample_Name"]:
            emit(f"  {sample}")
    emit(df.to_string(index=False))
