# Copyright © 2026 CLISEQ LTD. All rights reserved.
# Generate folders and commands for processing MIP genotyping sequencing data
import os
import sys
import argparse

def _read_samples_file(samples_file, read1_str):
    """Read the requested sample names from samples_file (one per line, '#' comments ignored)."""
    requested = []
    with open(samples_file) as f:
        for line in f:
            # tolerate csv/tsv exports: take the first field of the line
            name = line.split(',')[0].split('\t')[0].strip()
            if name and not name.startswith('#'):
                # tolerate a whole R1 file name instead of a sample name
                requested.append(os.path.basename(name).split(read1_str)[0])
    if not requested:
        raise ValueError(f"No sample names found in --samples-file {samples_file}")
    return requested


def _match_sample(requested_name, available_samples):
    """Match a requested sample name against the sample names found in the FASTQ dir.

    Tries an exact match first, then a unique prefix match, then a unique substring match,
    so that a sample sheet holding 'S123' still matches 'S123_S7_L001'.
    """
    if requested_name in available_samples:
        return [requested_name]
    prefix_matches = [s for s in available_samples if s.startswith(requested_name)]
    if prefix_matches:
        return prefix_matches
    return [s for s in available_samples if requested_name in s]


def collect_samples(fastq_path, read1_str, read2_str, wildcard_str_in_sample_name, samples_file=None):
    """
    Build the ordered list of (sample_name, filename1, filename2) tuples to process.

    Without samples_file, all FASTQ pairs whose R1 file contains wildcard_str_in_sample_name
    are taken. With samples_file, only the listed samples are taken, in the order of the file,
    and every listed sample that has no R1 (or no matching R2) file is reported before aborting.
    """
    # Map every R1 FASTQ found in the directory to its sample name
    fastqs_by_sample = {}
    for filename1 in sorted(os.listdir(fastq_path)):
        if read1_str not in filename1 or 'Undetermined' in filename1:
            continue
        if not (filename1.endswith('.fastq') or filename1.endswith('.fastq.gz')):
            continue
        sample_name = filename1.split(read1_str)[0]
        if sample_name in fastqs_by_sample:
            print(f"Warning: several R1 files map to sample {sample_name}, "
                  f"using {fastqs_by_sample[sample_name]} and ignoring {filename1}")
            continue
        fastqs_by_sample[sample_name] = filename1

    missing = []
    if samples_file is None:
        selected = [(sample_name, filename1) for sample_name, filename1 in fastqs_by_sample.items()
                    if wildcard_str_in_sample_name in filename1]
    else:
        selected = []
        already_selected = set()
        for requested_name in _read_samples_file(samples_file, read1_str):
            matches = _match_sample(requested_name, fastqs_by_sample)
            if not matches:
                missing.append(f"{requested_name}: no R1 file matching it in {fastq_path}")
            elif len(matches) > 1:
                missing.append(f"{requested_name}: ambiguous, matches {', '.join(sorted(matches))}")
            elif matches[0] in already_selected:
                print(f"Warning: sample {matches[0]} requested more than once, keeping it once")
            else:
                already_selected.add(matches[0])
                selected.append((matches[0], fastqs_by_sample[matches[0]]))

    # Resolve the R2 mate of every selected sample
    samples = []
    for sample_name, filename1 in selected:
        filename2 = filename1.replace(read1_str, read2_str)
        if not os.path.exists(os.path.join(fastq_path, filename2)):
            missing.append(f"{sample_name}: R1 file {filename1} has no R2 mate {filename2}")
            continue
        samples.append((sample_name, filename1, filename2))

    if missing:
        print(f'Error: {len(missing)} sample(s) cannot be processed:', file=sys.stderr)
        for problem in missing:
            print(f'  - {problem}', file=sys.stderr)
        raise ValueError(f"{len(missing)} sample(s) missing from {fastq_path}, see the list above")

    if not samples:
        raise ValueError(f"No samples to process were found in {fastq_path}")

    return samples


def generate_genotype_cmds(run_base_dir,
                           analysis_type = 'genotype',
                           fastq_dir = 'fastq', 
                           output_dir = 'analysis_geno',
                           read1_str = '_R1',
                           read2_str = '_R2',
                           n_samples_per_job = 5,
                           wildcard_str_in_sample_name = 'Geno',
                           samples_file = None,
                           downsample_to = None,
                           job_queue = 'gsla-cpu',
                           cmd1_mem = 7000,
                           cmd2_mem = 30000,
                           n_cores1 = 32,
                           n_cores2 = 5,
                           apptainer_containers_base_dir = '/home/projects/shlush/shared/singularity_aka_apptainer_containers',
                           container_img = None, #'ops:5000/genotyping_v3.0.0:v1.16.18',
                           run_mode = 'apptainer',
                           scheduler = 'lsf',
                           genome_ref_fasta_path = '/home/projects/shlush/shared/Homo_sapiens.GRCh38.dna_sm.primary_assembly.genomeFile/hg38.fa',
                           vep_path = '/home/projects/shlush/shared/homo_sapiens_merged_vep_104_GRCh38/'
                           ):
    """
    Generate directories and command lines for MIP genotyping data processing.

    Parameters:
    run_base_dir (str): Base directory for the sequencing run.
    analysis_type (str): Analysis type: genotype or arch.
    fastq_dir (str): Directory containing FASTQ files.
    output_dir (str): Directory to store output.
    read1_str (str): Read 1 identifier string.
    read2_str (str): Read 2 identifier string.
    wildcard_str_in_sample_name (str): Wildcard string in sample name to filter by.
    samples_file (str): Path to a file listing the samples to run (one sample name per line).
                        Overrides wildcard_str_in_sample_name; aborts if an R1/R2 file of a
                        listed sample is missing from the FASTQ directory.
    downsample_to (float): Downsample every sample to this many million reads (None = no downsampling).
    n_samples_per_job (int): Number of samples to process per job.
    job_queue (str): Job queue name.
    cmd1_mem (int): Memory for command 1 (MB).
    cmd2_mem (int): Memory for command 2 (MB).
    n_cores1 (int): Number of cores for command 1.
    n_cores2 (int): Number of cores for command 2.
    container_img (str): SIF exported docker image.
    run_mode (str): Container run mode, either 'docker' or 'apptainer'.
    scheduler (str): Cluster job scheduler, either 'lsf' or 'sge'.
    genome_ref_fasta_path (str): Path to genome reference FASTA.
    vep_path (str): Path to VEP data.

    Returns:
    Creates folders per job and the commands for cluster submission (currently breaking down to 2 steps)
    """
    if scheduler not in ('lsf', 'sge'):
        raise ValueError(f"Invalid scheduler: {scheduler}, must be 'lsf' or 'sge'")
    if run_mode == 'docker' and scheduler == 'sge':
        raise ValueError("run_mode='docker' is not supported with scheduler='sge' "
                          "(SGE has no equivalent of LSF's -app docker-cpu integration); use run_mode='apptainer'")

    if samples_file is not None and not os.path.exists(samples_file):
        raise ValueError(f"--samples-file {samples_file} does not exist")

    # Create output directory if it doesn't exist
    output_path = os.path.join(run_base_dir, output_dir)
    if os.path.exists(output_path):
        print(f'Error: output folder {output_path} already exists, delete it manually if you want this command to overwrite it, Aborting...')
        return

    # Resolve the samples to run before creating anything, so a bad --samples-file leaves no leftovers
    fastq_path = os.path.join(run_base_dir, fastq_dir)
    samples = collect_samples(fastq_path, read1_str, read2_str, wildcard_str_in_sample_name, samples_file)
    if samples_file is not None:
        print(f'Selected {len(samples)} sample(s) from {samples_file} (ignoring --wildcard-str-in-sample-name)')

    os.makedirs(output_path, exist_ok=True)

    # bsub -q gsla-cpu -R "rusage[mem=30000] span[hosts=1]" -n 5 -J NovaseqR85_genotype_10 -app docker-cpu -env LSB_CONTAINER_IMAGE=ops:5000/genotyping_v3.0.0:v1.16.18 analyze --cores 5 --samples-path /home/projects/shlush/shared/runs_analysis/NovaseqR85/analysis_genotype/job10/Sample_file_path_info.yaml        --output-folder /home/projects/shlush/shared/runs_analysis/NovaseqR85/analysis_genotype/job10        --analysis NovaseqR85_genotype_10 --reference /home/projects/shlush/shared/Homo_sapiens.GRCh38.dna_sm.primary_assembly.genomeFile/hg38.fa                                                                                        --variant-calling on -v
    # bsub -q medium   -R "rusage[mem=30000] span[hosts=1]" -n 5 -J MiseqR4_1              -app docker-cpu -env LSB_CONTAINER_IMAGE=ops:5000/genotyping_v3.0.0:v1.16.18 analyze --cores 5 --samples-path /home/projects/shlush/USER/shared_shlush/runs_analysis/MiseqR4/analysis_geno/job1/Sample_file_path_info.yaml --output-folder /home/projects/shlush/USER/shared_shlush/runs_analysis/MiseqR4/analysis_geno/job1 --analysis MiseqR4_1              --reference /home/projects/shlush/shared/Homo_sapiens.GRCh38.dna_sm.primary_assembly.genomeFile/hg38.fa --config vep='{dir: /home/projects/shlush/shared/homo_sapiens_merged_vep_104_GRCh38/}' -v --variant-calling on

    # List to hold command lines
    cmds1 = []
    cmds2 = []
    sample_bams = []

    # Container invocation and the --config entries are the same for every job, so build them once.
    # The analyze command accepts a single --config flag, so vep and downsample share it.
    config_entries = []
    if run_mode == 'apptainer':
        container_path = os.path.join(apptainer_containers_base_dir, container_img)
        container_prefix = f"apptainer exec --bind {vep_path}:/vep {container_path} "
        # no vep config entry: default vep dir is /vep, we bind the vep_path to /vep in the container, as it's problematic to pass it as a parameter to the containerized analyze command, so we just bind it to /vep and let the analyze command use the default /vep path
    elif run_mode == 'docker':
        container_prefix = f"-app docker-cpu -env LSB_CONTAINER_IMAGE={container_img} "
        config_entries.append(f"vep='{{dir: {vep_path}}}'")
    else:
        raise ValueError(f"Invalid run_mode: {run_mode}, must be 'docker' or 'apptainer'")

    if downsample_to is not None:
        # X is given in millions of reads and passed through as such, e.g. --downsample-to 1 -> downsample=1
        downsample_val = int(downsample_to) if float(downsample_to).is_integer() else downsample_to
        config_entries.append(f"downsample={downsample_val}")
        print(f'Downsampling every sample to {downsample_val}M reads')

    config_flag = f"--config {' '.join(config_entries)} " if config_entries else ""

    # Iterate over the selected samples
    n_samples = 0
    for sample_name, filename1, filename2 in samples:
        job_num = n_samples // n_samples_per_job + 1
        print(f'Processing file: {filename1} into job {job_num}')

        # New job - create job directory and commands
        if n_samples % n_samples_per_job == 0:
            # Create job directory
            job_dir = os.path.join(output_path, f'job{job_num}')
            os.makedirs(job_dir, exist_ok=True)

            job_id = f"{analysis_type}_{os.path.basename(run_base_dir)}_{job_num}"
            if scheduler == 'lsf':
                submit_bin = "bsub"
                job_name_flag = f"-J {job_id}"
                mem1_flag = f'-R "rusage[mem={cmd1_mem}] span[hosts=1]"'
                mem2_flag = f'-R "rusage[mem={cmd2_mem}] span[hosts=1]"'
                cores1_flag = f"-n {n_cores1}"
                cores2_flag = f"-n {n_cores2}"
            else:  # sge
                submit_bin = "qsub -V -b y"
                job_name_flag = f"-N {job_id}"
                mem1_flag = f"-l mem_free={cmd1_mem}mb"
                mem2_flag = f"-l mem_free={cmd2_mem}mb"
                cores1_flag = f"-pe smp {n_cores1}"
                cores2_flag = f"-pe smp {n_cores2}"

            cmd1 = (f"{submit_bin} "
                   f"-q {job_queue} "
                   f"{mem1_flag} "
                   f"{cores1_flag} "
                   f"{job_name_flag} "
                   f"-o {os.path.join(job_dir, f'{os.path.basename(run_base_dir)}_{job_num}_J1.out')} "
                   f"-e {os.path.join(job_dir, f'{os.path.basename(run_base_dir)}_{job_num}_J1.err')} "
                   f"{container_prefix}"
                   "analyze "
                   f"--cores {n_cores1} "
                   f"--samples-path {os.path.join(job_dir, 'Sample_file_path_info.yaml')} "
                   f"--output-folder {job_dir} "
                   f"--analysis {job_id} "
                   f"--reference {genome_ref_fasta_path} "
                   f"{config_flag}"
                   "-v")
            cmds1.append(cmd1)

            cmd2 = cmd1.replace(mem1_flag, mem2_flag) \
                       .replace(cores1_flag, cores2_flag) \
                       .replace(f'--cores {n_cores1}', f'--cores {n_cores2}') \
                       .replace('_J1.', '_J2.') \
                       .replace(' -v', ' --variant-calling on -v')
            cmds2.append(cmd2)

        # Append sample info to job's sample file
        sample_file_path = os.path.join(job_dir, 'Sample_file_path_info.yaml')
        with open(sample_file_path, 'a') as sample_file:
            sample_file.write(f"{sample_name}:\n")
            sample_file.write(f"  R1: {os.path.join(fastq_path, filename1)}\n")
            sample_file.write(f"  R2: {os.path.join(fastq_path, filename2)}\n")

        sample_bams.append(os.path.join(job_dir, 'paired', '2_mapped', sample_name, f"{sample_name}.filtered.sorted.rg.realigned.rn.bam"))

        n_samples += 1

    # Write commands to files
    # Commands for step 1
    cmds1_file = os.path.join(output_path, f'{analysis_type}_cmds_step1.sh')
    with open(cmds1_file, 'w') as f:
        for cmd in cmds1:
            f.write('sleep 1\n')
            f.write(cmd + '\n')
    
    # Commands for step 2
    cmds2_file = os.path.join(output_path, f'{analysis_type}_cmds_step2.sh')
    with open(cmds2_file, 'w') as f:
        for cmd in cmds2:
            f.write('sleep 1\n')
            f.write(cmd + '\n')

    # Write sample bams to file
    with open(os.path.join(output_path, "sample_bams.txt"), 'w') as f:
        for sample_bam in sample_bams:
            f.write(sample_bam + '\n')
    
    # Define your input and output filenames
    samples_file = os.path.join(output_path, "sample_bams.txt")
    array_ext = 'bsub' if scheduler == 'lsf' else 'qsub'
    array_file = os.path.join(output_path, f"submit_lib_size_estimate_array.{array_ext}")

    # 2. Define the array-job scheduler directives and job/task-index env vars
    if scheduler == 'lsf':
        array_header = f"""#BSUB -J "picard_arr[1-{n_samples}]"
    #BSUB -q short
    #BSUB -n 2
    #BSUB -R "span[ptile=1]"
    #BSUB -M 16384
    #BSUB -R "rusage[mem=16384]"
    #BSUB -o {output_path}/logs/picard_%J_%I.out
    #BSUB -e {output_path}/logs/picard_%J_%I.err"""
        task_var, jobid_var = "LSB_JOBINDEX", "LSB_JOBID"
    else:  # sge; reuses job_queue since there's no separate SGE mapping for the LSF-hardcoded 'short' queue
        array_header = f"""#$ -N picard_arr
    #$ -t 1-{n_samples}
    #$ -q {job_queue}
    #$ -pe smp 2
    #$ -l mem_free=16384mb
    #$ -o {output_path}/logs/picard_$JOB_ID_$TASK_ID.out
    #$ -e {output_path}/logs/picard_$JOB_ID_$TASK_ID.err
    #$ -cwd
    #$ -V"""
        task_var, jobid_var = "SGE_TASK_ID", "JOB_ID"

    # 3. Define the array script using an f-string
    # Note: Double curly braces {{ }} are used to safely escape bash variables like ${LSB_JOBINDEX} in Python f-strings.
    array_content = f"""#!/bin/bash
    {array_header}

    # module load java/1.8
    module load picard

    # Ensure output directories exist
    mkdir -p {output_path}/logs {output_path}/results

    # Grab the BAM path corresponding to this job's unique array index
    SAMPLE_BAM=$(sed -n "${{{task_var}}}p" {samples_file})

    # Extract a clean base name so each job writes a unique output file
    BASE_NAME=$(basename "$SAMPLE_BAM" .bam)
    OUTPUT_METRICS="{output_path}/results/${{BASE_NAME}}_complexity_metrics.txt"

    # Create a unique scratch directory for this specific sub-job
    export PICARD_TMP="/tmp/picard_${{{jobid_var}}}_${{{task_var}}}"
    mkdir -p "$PICARD_TMP"

    echo "Starting Array Index: ${task_var}"
    echo "Processing File: $SAMPLE_BAM"
    echo "Writing Output To: $OUTPUT_METRICS"

    # Run Picard safely
    java -Xmx14g -XX:ParallelGCThreads=2 -jar $EBROOTPICARD/picard.jar EstimateLibraryComplexity \\
        -I "$SAMPLE_BAM" \\
        -O "$OUTPUT_METRICS" \\
        -TMP_DIR "$PICARD_TMP" \\
        -VALIDATION_STRINGENCY SILENT \\
        -MAX_GROUP_RATIO 10000 \\
        -MIN_IDENTICAL_BASES 10

    # Clean up the temporary scratch space
    rm -rf "$PICARD_TMP"

    echo "Index ${task_var} finished successfully."
    """

    # 4. Write the script to a file
    with open(array_file, "w") as f:
        f.write(array_content)

    print(f"Successfully generated '{array_file}'.")
    if scheduler == 'lsf':
        print(f"You can now submit it to the cluster by running: bsub < {array_file}")
    else:
        print(f"You can now submit it to the cluster by running: qsub {array_file}")
    return

def main():
    
    parser = argparse.ArgumentParser(description='Generate MIP genotyping processing commands')
    parser.add_argument('run_base_dir', help='Base directory for the sequencing run')
    parser.add_argument('--analysis-type', default='genotype', choices=['genotype', 'arch'], help='Analysis type: genotype or arch')
    parser.add_argument('--fastq-dir', default='fastq', help='Directory containing FASTQ files')
    parser.add_argument('--output-dir', default='analysis_geno', help='Directory to store output')
    parser.add_argument('--read1-str', default='_R1', help='Read 1 identifier string')
    parser.add_argument('--read2-str', default='_R2', help='Read 2 identifier string')
    parser.add_argument('--wildcard-str-in-sample-name', default='Geno', help='Wildcard string in sample name to filter by')
    parser.add_argument('--samples-file', default=None, help='File listing the samples to run (one sample name per line). Overrides --wildcard-str-in-sample-name, aborts if the R1/R2 file of a listed sample is missing from the FASTQ dir')
    parser.add_argument('--downsample-to', type=float, default=None, help='Downsample every sample to this many million reads (adds downsample=X to --config, default: no downsampling)')
    parser.add_argument('--n-samples-per-job', type=int, default=5, help='Number of samples per job')
    parser.add_argument('--job-queue', default=None, help="Job queue name (default: 'gsla-cpu' for lsf, 'all.q' for sge)")
    parser.add_argument('--cmd1-mem', type=int, default=7000, help='Memory for command 1 (MB)')
    parser.add_argument('--cmd2-mem', type=int, default=30000, help='Memory for command 2 (MB)')
    parser.add_argument('--n-cores1', type=int, default=32, help='Number of cores for command 1')
    parser.add_argument('--n-cores2', type=int, default=5, help='Number of cores for command 2')
    # parser.add_argument('--container-img', default='ops:5000/genotyping_v3.0.0:v1.16.18', help='Docker container image')
    parser.add_argument('--use-container-img', default='', help='Override default container image')
    parser.add_argument('--run-mode', default='apptainer', choices=['docker', 'apptainer'], help='Container run mode: docker or apptainer')
    parser.add_argument('--scheduler', default='lsf', choices=['lsf', 'sge'], help='Cluster job scheduler: lsf or sge')
    parser.add_argument('--genome-ref-fasta-path', default='/home/projects/shlush/shared/Homo_sapiens.GRCh38.dna_sm.primary_assembly.genomeFile/hg38.fa', help='Path to genome reference FASTA')
    parser.add_argument('--vep-path', default='/home/projects/shlush/shared/homo_sapiens_merged_vep_104_GRCh38/', help='Path to VEP data')
    
    args = parser.parse_args()
    
    if args.analysis_type == 'genotype':
        # output_dir = 'analysis_geno'
        if args.use_container_img == '':
            if args.run_mode == 'docker':
                container_img = 'ops:5000/genotyping_v3.0.0:v1.16.18'
            else:
                container_img = 'genotyping_v3.0.0_v1.16.18.sif'
        else:
            container_img = args.use_container_img
        # container_img = 'ops:5000/genotyping_v3.0.0:v1.17.13'
    elif args.analysis_type == 'arch':
        # output_dir = 'analysis_arch'
        # container_img = 'ops:5000/arch_v5.0.3c:v1.13.6'
        # container_img = 'ops:5000/archv7.1.1_isp500:v1.16.16'
        if args.use_container_img == '':
            if args.run_mode == 'docker':
                container_img = 'ops:5000/myeloid_melignancies_isp510:v1.17.13'
            else:
                container_img = 'myeloid_melignancies_isp510_v1.17.13.sif'
        else:
            container_img = args.use_container_img

    else:
        raise ValueError(f'Invalid analysis type: {args.analysis_type}')

    if args.job_queue is not None:
        job_queue = args.job_queue
    else:
        job_queue = 'gsla-cpu' if args.scheduler == 'lsf' else 'all.q'

    try:
        generate_genotype_cmds(
            args.run_base_dir,
            analysis_type=args.analysis_type,
            fastq_dir=args.fastq_dir,
            output_dir=args.output_dir,
            read1_str=args.read1_str,
            read2_str=args.read2_str,
            wildcard_str_in_sample_name=args.wildcard_str_in_sample_name,
            samples_file=args.samples_file,
            downsample_to=args.downsample_to,
            n_samples_per_job=args.n_samples_per_job,
            job_queue=job_queue,
            cmd1_mem=args.cmd1_mem,
            cmd2_mem=args.cmd2_mem,
            n_cores1=args.n_cores1,
            n_cores2=args.n_cores2,
            container_img=container_img,
            run_mode=args.run_mode,
            scheduler=args.scheduler,
            genome_ref_fasta_path=args.genome_ref_fasta_path,
            vep_path=args.vep_path
        )
    except ValueError as err:
        # invalid arguments / missing samples are user errors, report them without a traceback
        print(f'Error: {err}, Aborting...', file=sys.stderr)
        sys.exit(1)


if __name__ == '__main__':
    main()


# bsub -q medium -R "rusage[mem=7000] span[hosts=1]" -n 32 -J NovaseqR84_1       -app docker-cpu -env LSB_CONTAINER_IMAGE=ops:5000/genotyping_v3.0.0:v1.16.18 analyze --cores 32 --samples-path /home/projects/shlush/shared/runs_analysis/NovaseqR84/analysis_geno_Nat/job1/Sample_file_path_info.yaml --output-folder /home/projects/shlush/shared/runs_analysis/NovaseqR84/analysis_geno_Nat/job1 --analysis NovaseqR84_1 --reference /home/projects/shlush/shared/Homo_sapiens.GRCh38.dna_sm.primary_assembly.genomeFile/hg38.fa --config vep="{dir: /home/projects/shlush/shared/homo_sapiens_merged_vep_104_GRCh38/}" -v
# bsub -q medium -R "rusage[mem=7000] span[hosts=1]" -n 32 -J NEX45_Sequentify_1 -app docker-cpu -env LSB_CONTAINER_IMAGE=ops:5000/genotyping_v3.0.0:v1.16.18 analyze --cores 32 --samples-path /home/projects/shlush/USER/shared_shlush/runs_analysis/NEX45_Sequentify/analysis_geno/job1/Sample_file_path_info.yaml --output-folder /home/projects/shlush/USER/shared_shlush/runs_analysis/NEX45_Sequentify/analysis_geno/job1 --analysis NEX45_Sequentify_1 --reference /home/projects/shlush/shared/Homo_sapiens.GRCh38.dna_sm.primary_assembly.genomeFile/hg38.fa --config vep='{dir: /home/projects/shlush/shared/homo_sapiens_merged_vep_104_GRCh38/}' -v
# arch
# bsub -q medium -R   "rusage[mem=5000] span[hosts=1]"  -n 32 -J Nili5              -app docker-cpu -env LSB_CONTAINER_IMAGE=ops:5000/arch_v5.0.3c:v1.13.6       analyze --cores 32 --samples-path /home/labs/shlush/shared/runs_analysis/NovaseqR40/analysis_ARCH/job5/Sample_file_path_info.yaml         --output-folder /home/labs/shlush/shared/runs_analysis/NovaseqR40/analysis_ARCH/job5/        --analysis  Nili5             --reference /home/labs/shlush/shared/Homo_sapiens.GRCh38.dna_sm.primary_assembly.genomeFile/hg38.fa     --vep  /home/labs/shlush/shared/homo_sapiens_merged_vep_104_GRCh38/  -v
# bsub -q gsla-cpu -R "rusage[mem=7000] span[hosts=1]" -n 32 -J NovaseqR85_arch_1 -app docker-cpu   -env LSB_CONTAINER_IMAGE=ops:5000/archv7.1.1_isp500:v1.16.16 analyze --cores 32 --samples-path /home/projects/shlush/shared/runs_analysis/NovaseqR85/analysis_arch/job1/Sample_file_path_info.yaml     --output-folder /home/projects/shlush/shared/runs_analysis/NovaseqR85/analysis_arch/job1     --analysis NovaseqR85_arch_1  --reference /home/projects/shlush/shared/Homo_sapiens.GRCh38.dna_sm.primary_assembly.genomeFile/hg38.fa --config "vep={dir: /home/projects/shlush/shared/homo_sapiens_merged_vep_104_GRCh38/}" -v