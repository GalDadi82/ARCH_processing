"""Apply the Blood-2023-paper filtering to the merged mutation table.

Thin wrapper around utils/filter_mutations.py (the Python port of
Mutation_Curation/R/FilterMutations.R, verified equivalent to the R version).
"""

import sys

sys.path.insert(0, snakemake.params.utils_dir)

from filter_mutations import filter_mutations  # noqa: E402

with open(snakemake.log[0], "w") as log:
    result = filter_mutations(
        input_file=snakemake.input.merged,
        output_dir=snakemake.params.output_dir,
        run_name=snakemake.params.run_name,
        sample_column=snakemake.params.sample_column,
        dups_input_file=snakemake.params.dups_input_file,
    )

    for label, key in (("listed", "listed_df"), ("not listed", "rest_df")):
        msg = f"{label}: {len(result[key])} rows"
        print(msg)
        log.write(msg + "\n")

    # filter_mutations builds its own file names from run_name; make sure they are the
    # ones the rule declared, so a config change cannot silently desync them.
    for produced, declared in (
        (result["listed_fn"], snakemake.output.listed),
        (result["rest_fn"], snakemake.output.not_listed),
    ):
        if produced != declared:
            sys.exit(f"filter_mutations wrote {produced}, expected {declared}")
