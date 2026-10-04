# ARCH_processing

Different utils and pipelines for processing ARCH sequencing run, initially adapted from Shlush lab code.

## 1. Generating job commands (apptainer/docker)

The script *utils/gen_mip_processing_cmds.py* generates the folder for running the containerized MIP jobs under the sequencing run folder.

Its default mode is actually for genotype runs, it only needs the path to the sequencing folder:

`python ./gen_mip_processing_cmds.py /home/projects/shlush/shared/runs_analysis/NovaseqR999`

This will create a folder named *analysis_geno* under *NovaseqR999* (override with `--output-dir`) with jobs subfolders and the 2 job files for the first and second step of the pipeline.

**Note:** the generated commands now run via **apptainer** by default (`apptainer exec`, binding the VEP data directory into the container at `/vep`). Pass `--run-mode docker` to fall back to the previous docker-based commands (an LSF `-app docker-cpu` job flag, not a literal `docker run`) — docker mode only works with `--scheduler lsf` (not supported on `sge`). Apptainer `.sif` images are expected under `/home/projects/shlush/shared/singularity_aka_apptainer_containers` by default; override with `--use-container-img`.

The script full usage:

```
usage: gen_mip_processing_cmds.py [-h] [--analysis-type {genotype,arch}]
                                  [--fastq-dir FASTQ_DIR]
                                  [--output-dir OUTPUT_DIR]
                                  [--read1-str READ1_STR]
                                  [--read2-str READ2_STR]
                                  [--wildcard-str-in-sample-name WILDCARD_STR_IN_SAMPLE_NAME]
                                  [--n-samples-per-job N_SAMPLES_PER_JOB]
                                  [--job-queue JOB_QUEUE]
                                  [--cmd1-mem CMD1_MEM] [--cmd2-mem CMD2_MEM]
                                  [--n-cores1 N_CORES1] [--n-cores2 N_CORES2]
                                  [--use-container-img USE_CONTAINER_IMG]
                                  [--run-mode {docker,apptainer}]
                                  [--scheduler {lsf,sge}]
                                  [--genome-ref-fasta-path GENOME_REF_FASTA_PATH]
                                  [--vep-path VEP_PATH]
                                  run_base_dir

Generate MIP genotyping processing commands

positional arguments:
  run_base_dir          Base directory for the sequencing run

options:
  -h, --help            show this help message and exit
  --analysis-type {genotype,arch}
                        Analysis type: genotype or arch
  --fastq-dir FASTQ_DIR
                        Directory containing FASTQ files
  --output-dir OUTPUT_DIR
                        Directory to store output (default: analysis_geno / analysis_arch)
  --read1-str READ1_STR
                        Read 1 identifier string (default: _R1)
  --read2-str READ2_STR
                        Read 2 identifier string (default: _R2)
  --wildcard-str-in-sample-name WILDCARD_STR_IN_SAMPLE_NAME
                        Wildcard string in sample name to filter by (default: Geno)
  --n-samples-per-job N_SAMPLES_PER_JOB
                        Number of samples per job (default: 5)
  --job-queue JOB_QUEUE
                        Job queue name (default: 'gsla-cpu' for lsf, 'all.q' for sge)
  --cmd1-mem CMD1_MEM   Memory for command 1 (MB), default: 7000
  --cmd2-mem CMD2_MEM   Memory for command 2 (MB), default: 30000
  --n-cores1 N_CORES1   Number of cores for command 1 (default: 32)
  --n-cores2 N_CORES2   Number of cores for command 2 (default: 5)
  --use-container-img USE_CONTAINER_IMG
                        Override default container image
  --run-mode {docker,apptainer}
                        Container run mode: docker or apptainer (default: apptainer)
  --scheduler {lsf,sge}
                        Cluster job scheduler: lsf or sge (default: lsf)
  --genome-ref-fasta-path GENOME_REF_FASTA_PATH
                        Path to genome reference FASTA (default: /home/projects/shlush/shared/Homo_sapiens.GRCh38.dna_sm.primary_assembly.genomeFile/hg38.fa)
  --vep-path VEP_PATH   Path to VEP data (default: /home/projects/shlush/shared/homo_sapiens_merged_vep_104_GRCh38/)
```

So a basic *arch* run looks like this:
`python ./gen_mip_processing_cmds.py --analysis-type arch --wildcard-str-in-sample-name ARCH /home/projects/shlush/shared/runs_analysis/NovaseqR999`

(this generates apptainer-based job commands by default; add `--run-mode docker` for the previous docker-based behavior)

## 2. Full ARCH processing pipeline (notebook)

All ARCH processing steps now live in a single notebook, *notebooks/arch_processing.ipynb* — merging sample VCFs, annovar annotation, Blood-2023-paper filtering, multi-sample/COSMIC recurrence filtering, CALR mutation detection, and the mutation email alert. (Previously split across *arch_processing_step1.ipynb* and *arch_processing_step2_multi_sample_muts.ipynb* — both removed.) You should **update** the paths at the top of the notebook, and execute cells sequentially.

The notebook's stages, in order:
- Set up per-run paths and the `5_vcfs`/`2_mapped` sample-linking folders; generate QC/statistics via *utils/mip_run_statistics.py*.
- Annotate each sample's varscan VCF with annovar (an interactive cell, or an LSF-per-sample version using *run_annovar_one_lsf.sh*), then convert + merge all samples into *merged_varscan_annotations.txt*.
- Apply the Blood-2023-paper filtering (Python port of *Mutation_Curation/R/FilterMutations.R*, via *utils/filter_mutations.py*), producing *FilterMutations_{seq_run}_arch_listed.tsv*.
- Remove recurrent mutations without COSMIC/ClinVar support, apply the final VAF filter (`> 0.05`) and drop known artefact/germline calls — the master output is *FilterMutations_{seq_run}_arch_listed_FilteredRecurrent_VAF_0.05.tsv*.
- Detect CALR mutations and collect results (see section 3 below).
- Cross-reference against Denver-cohort normals.
- Run the mutation email alert (see section 4 below).

Notes:
- Annovar script and database have been updated — this changes the output format (column names, new columns).
- You can run annovar in the notebook (annotation cell) but it's best to submit a job per sample (the LSF version cell).
- Mutations that passed all filters up to the Blood-2023 filtering are in the *..._filtered_yet_not_listed.tsv* file; the fully filtered set is in *..._listed.tsv*.
- We still don't have the Denver-dataset poisson model to fully filter false mutations, so manually inspect the output — be suspicious of mutations shared across many samples with a relatively uniform VAF.

## 2b. Same pipeline as a Snakemake workflow

*workflow/* runs everything section 2 describes without the notebook: one command, resumable,
one cluster job per sample instead of the notebook's bsub-submission cells. The notebook is
unchanged and still works — use whichever you prefer.

Setup (snakemake is not in any of the existing conda envs):

```
conda env create -f workflow/envs/arch.yaml
conda activate arch_snakemake
```

Run a whole sequencing run:

```
snakemake -s workflow/Snakefile --profile workflow/profiles/lsf --config seq_run=NovaseqR132
```

All the parameters that used to live in the notebook's first cell are now in
*workflow/config/config.yaml* (run name, panel, QC cutoffs, annovar databases, VAF
thresholds, the excluded-mutation list, CALR settings, Denver paths, alert config).
Edit that file, or override single values with `--config key=value`.

Useful flags and targets:
- `-n` — dry run, prints what would be done. Always worth running first.
- `snakemake ... qc` — stop after the QC report and sample selection. Use this when sample
  names need the manual `EXPERIMENT_PANEL_SAMPLE_INDEX` fix described above: run `qc`, fix
  *{panel}_{seq_run}_statistics.csv*, then run the full pipeline (`qc.rebuild` stays `false`,
  so your edit is not overwritten).
- `snakemake ... mutations` — everything except the alert email.
- Other targets: `annotate`, `calr`, `denver`.
- Without `--profile`, it runs locally: `snakemake -s workflow/Snakefile --cores 8`.

Notes:
- The sample list is discovered at runtime by a Snakemake *checkpoint*, so a single
  invocation runs QC, picks the samples and then fans out per-sample ANNOVAR and CALR jobs.
- Outputs land in the same `All_Samples/5_vcfs` folder, under the same file names the
  notebook produced. Per-rule logs go to `5_vcfs/logs/`.
- One extra intermediate appears that the notebook kept only in memory:
  *FilterMutations_{seq_run}_arch_listed.tsv* (the recurrence step's output).
- The alert email is part of the default run, as in the notebook. `snakemake ... mutations`
  or `alert.enabled: false` stops before it. Delete *mutation_alert_{seq_run}.sent* to
  re-send.
- Re-running only redoes what is out of date. To force a step, delete its output (or use
  `--forcerun <rule>`).

## 3. CALR mutation processing

Two scripts detect the two known CALR ARCH-relevant mutation types. Both queue one LSF job per sample (submit mode) and are also run automatically as part of *arch_processing.ipynb* (section 2 above).

- *utils/CALR_type1_from_fastq.py* — detects the CALR Type-1 52bp deletion (del52) directly from FASTQ, by matching junction/breakpoint motifs built from a seed WT sequence.
  Submit: `python CALR_type1_from_fastq.py --submit --fastq_dir /path/to/fastqs` (filters samples via `--sample_pattern`, default `ARCH`).
  Per-sample output: *calr_results/{sample}_CALR.tsv* (`WT_count`, `DEL52_count`, `VAF_DEL52`, ...).
- *utils/CALR_type2_from_bam.py* — detects the CALR Type-2 5bp insertion (ins5) from BAM alignments, scanning a fixed chr19 pileup window for a CIGAR insertion matching a known motif, and reporting the most-supported breakpoint position.
  Submit: `python CALR_type2_from_bam.py --bam_dir /path/to/samples` (same `--sample_pattern` filtering).
  Per-sample output: *calr_results/{sample}_CALR_type2.tsv* (`CHR`, `POS`, `POS_frac`, `INS5_count`, `VAF_INS5`, ...).

Merge either script's per-sample outputs with e.g.:
`awk '(NR==1)||(FNR>1)' calr_results/*_CALR.tsv > CALR_type1_fromFASTQ_merged.tsv`

(the notebook instead collects per-sample files directly into *CALR_type1_results_{seq_run}.tsv* / *CALR_type2_results_{seq_run}.tsv*).

## 4. Mutation email alert

*utils/mutation_alert.py*, configured via *utils/mutation_alert_config.yaml*, sends an email whenever:
- a mutation matches one of the gene/variant filters listed under `target_mutations` in the config (e.g. a specific `Gene.refGene`, or a specific `{CHR, POS, REF, ALT}`/`Selected.AA.change` combination), or
- a CALR Type-1/Type-2 result exceeds `min_CALR_del52_VAF` / `min_CALR_ins5_VAF` (default 0.05 each).

Edit `target_mutations` in *utils/mutation_alert_config.yaml* to add/remove genes or variants to watch for; edit `recipients`/`smtp` to change who gets alerted and how mail is sent. The email body summarizes each triggered section and attaches full CSVs (target mutations, CALR type-1, CALR type-2).

Run standalone with:
`python mutation_alert.py --mutations-file /path/to/FilterMutations_..._VAF_0.05.tsv --config mutation_alert_config.yaml --run-name <seq_run> --calr-type1-file CALR_type1_results_<seq_run>.tsv --calr-type2-file CALR_type2_results_<seq_run>.tsv`

(the notebook runs this automatically as its final cell, once the master mutation file and both CALR result files exist.)
