#!/usr/bin/env Rscript

# Initialisation ----

# Load libraries
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(logr))
suppressPackageStartupMessages(library(glue))

source("mut_utils.R")

# Get user arguments

parser = ArgumentParser(description = 'Filter mutations from multiple samples (annovar output) by VAF, mutation type (e.g. stop gain/loss, slicing), coordinate within specific genes and a curated list of mutations (hotspots).')

input_parser = parser$add_argument_group('Input parameters', '')
input_parser$add_argument('-input-file', type = 'character', nargs = 1, required = TRUE, 
                          help = 'Input mutation files (tab-delimited, merged from multiple samples, expecting annovar output columns).')
#input_parser$add_argument('-input-file', type = 'character', nargs = 1, required = FALSE, default = "/Users/eyal-l01/wexac_shared/runs_analysis/NovaseqR98/analysis_arch/All_Samples/5_vcfs/merged_varscan_annotations.txt",
#                          help = 'Input mutation files (tab-delimited, merged from multiple samples, expecting annovar output columns).')

input_parser$add_argument('-gene-column', type = 'character', nargs = 1, required = FALSE, default = 'Gene.refGene',
                          help = 'Gene column name. Default: Gene.refGene')
input_parser$add_argument('-chr-column', type = 'character', nargs = 1, required = FALSE, default = 'CHR',
                          help = 'Chromosome column name. Default: CHR')
input_parser$add_argument('-pos-column', type = 'character', nargs = 1, required = FALSE, default = 'POS',
                          help = 'Coordinate column name. Default: POS')
input_parser$add_argument('-sample-column', type = 'character', nargs = 1, required = FALSE, default = 'SAMPLE_ID',
                          help = 'Sample column name. Default: SAMPLE_ID')
#input_parser$add_argument('-sample-column', type = 'character', nargs = 1, required = FALSE, default = 'Sample_Name',
#                          help = 'Sample column name. Default: SAMPLE_ID')

input_parser$add_argument('-listed-mutations-file-name', type = 'character', nargs = 1, required = FALSE, default = "../data/Vlasschaert_Blood_2023_data.xlsx",
                          help = "XLSX file with listed mutaions ('Include' sheet) and those to exclude ('Exclude'). Default: ../data/Vlasschaert_Blood_2023_data.xlsx, Mandatory parameter.")

input_parser = parser$add_argument_group('Output parameters', '')
input_parser$add_argument('-output-dir', type = 'character', nargs = 1, required = TRUE, help = 'Output directory, multiple output files will be generated in it.')
#input_parser$add_argument('-output-dir', type = 'character', nargs = 1, required = FALSE, help = 'Output directory, multiple output files will be generated in it.', default = "/Users/eyal-l01/wexac_shared/runs_analysis/NovaseqR98/analysis_arch/All_Samples/5_vcfs")
input_parser$add_argument('-run-name', type = 'character', nargs = 1, required = TRUE, help = 'Name of run (e.g. NovaseqR84_ARCH), will be part of output files.')
#input_parser$add_argument('-run-name', type = 'character', nargs = 1, required = FALSE, help = 'Name of run (e.g. NovaseqR84_ARCH), will be part of output files.', default='NovaseqR98_arch')
input_parser$add_argument('-overwrite', type = 'logical', nargs = 1, required = FALSE, default = TRUE, help = 'Overwrite output files if exist. Default: TRUE')

filter_parser = parser$add_argument_group('Filtering parameters', '')
filter_parser$add_argument('-min-VAF', type = 'double', nargs = 1, required = FALSE, default = 0.005,
                           help = 'Min VAF to consider (per sample). Default: 0.005')
filter_parser$add_argument('-min-depth', type = 'integer', nargs = 1, required = FALSE, default = 20,
                           help = 'Min depth to consider (per sample). Default: 20')
filter_parser$add_argument('-max-gnomad41-AF', type = 'double', nargs = 1, required = FALSE, default = 0.001,
                           help = 'Max gnomad41 AF to consider (per sample). Default: 0.001')
filter_parser$add_argument('-remove-intronic', type = 'logical', nargs = 1, required = FALSE, default = TRUE,
                           help = 'Remove intronic mutations. Default: TRUE, will keep only exonic and splicing in Func.refGene')
#filter_parser$add_argument('-get-splicing', type = 'logical', nargs = 1, required = FALSE, default = TRUE,
#                           help = 'Get splicing mutations. Default: TRUE')
#filter_parser$add_argument('-get-stop-mutations', type = 'logical', nargs = 1, required = FALSE, default = TRUE,
#                           help = 'Get stop gain/loss mutations. Default: TRUE')
filter_parser$add_argument('-remove-synonymous', type = 'logical', nargs = 1, required = FALSE, default = TRUE,
                           help = 'Remove synonymous mutations. Default: TRUE')
filter_parser$add_argument('-include-non-syn-by-loci', type = 'logical', nargs = 1, required = FALSE, default = TRUE,
                           help = 'Include any AA change if it is in a listed AA position. Default: TRUE')

filter_parser$add_argument('-dups-input-file', type = 'logical', nargs = 1, required = FALSE, default = FALSE,
                           help = 'Both dup samples mutations in the same row, default: TRUE')
filter_parser$add_argument("-sample-specific-columns", type = 'character', nargs = 1, required = FALSE, default = 'VAF;Depth;Sample_Name',
                           help = 'Columns with sample specific values (semicolon delimited), relevant if -dups-input-file is FALSE. Default: VAF;Depth;Sample_Name')
filter_parser$add_argument('-min-supporting-samples', type = 'integer', nargs = 1, required = FALSE, default = 1,
                           help = 'Min number of samples supporting a mutation to be considered. Default: 1')
filter_parser$add_argument('-min-avg-VAF', type = 'double', nargs = 1, required = FALSE, default = -1,
                           help = 'Min average VAF across samples. Default: -1, meaning will adjust min VAF: 2 / (min-depth + depth)')

# TODO: Currently not used, update after Denver initial readout ----
filter_parser$add_argument('-over-represented-aa-changes-threshold', type = 'integer', nargs = 1, required = FALSE, default = 1e6, 
                           help = 'Remove over-represented AA changes that appear in more (>) samples than this threshold. Default: 5')
filter_parser$add_argument('-over-represented-splicing-mutations-threshold', type = 'integer', nargs = 1, required = FALSE, default = 1e6,
                           help = 'Remove over-represented splicing mutations that appear in more (>) samples than this threshold. Default: 5')

# str(commandArgs(trailingOnly = TRUE))
# for (i in seq_along(commandArgs(trailingOnly = TRUE))) {
#   if (commandArgs(trailingOnly = TRUE)[i] == "-h" | commandArgs(trailingOnly = TRUE)[i] == "--help") {
#     parser$print_help()
#     quit(status = 0)
#   }
# } 

args <- parser$parse_args()
# str(args)

# Open output
if (!dir.exists(args$output_dir)) {
  dir.create(args$output_dir, recursive = TRUE)
} 
output_fn = file.path(args$output_dir, glue("FilterMutations_{args$run_name}_listed.tsv"))
rest_fn = file.path(args$output_dir, glue("FilterMutations_{args$run_name}_filtered_yet_not_listed.tsv"))
log_fn = file.path(args$output_dir, glue("FilterMutations_{args$run_name}.log"))

if (!args$overwrite) {
  for (out_fn in c(output_fn , rest_fn, log_fn)) {
    if (file.exists(out_fn)) {
      stop(glue("Output file {out_fn} exists and overwrite is set to FALSE!"))
    }
  }
}

# Log function wrapper
log_wrap <- function(message, log_type = 'info') {
  if (log_type == 'info') {
    log_print(message)
  } else if (log_type == 'error') {
    log_print(message)
    stop(message)
  }
}
log_line = paste(rep("=", 80), collapse = "")

log_open(log_fn, logdir = FALSE)
log_wrap(glue("Starting FilterMutations for run {args$run_name}\n"))
log_wrap("User arguments:")
for (arg_name in names(args)) {
  log_wrap(glue("{arg_name}: {args[[arg_name]]}"))
}
log_wrap(log_line)

# Check data frame after filtering
check_df <- function(df) {
  if (nrow(df) == 0) {
    log_wrap("No mutations left after filtering!", log_type = 'error')
  }
}

# Read input file
df <- fread(args$input_file, data.table = FALSE)
required_columns <- c(args$gene_column, args$sample_column, args$chr_column, args$pos_column, 'REF', 'ALT', 'VAF', 'gnomad41_genome_AF', 'Func.refGene', 'AAChange.refGene')
if (args$dups_input_file) {
  required_columns <- c(required_columns, 'VAF_DUP', 'DEPTH_DUP')
} else {
  sample_specific_cols <- unlist(strsplit(args$sample_specific_columns, split = ";"))
  required_columns <- c(required_columns, sample_specific_cols)
}
for (col_names in required_columns) {
  if (!col_names %in% colnames(df)) {
    log_wrap(glue("Column {col_names} not found in input file!"), log_type = 'error')
  }
}
log_wrap(glue("Input file {args$input_file} read: {nrow(df)} rows, {ncol(df)} columns"))
check_df(df)

# Filter mutations ----
df_f <- df

# Filter by VAF
vaf_cols = if(args$dups_input_file) { c('VAF', 'VAF_DUP') } else {'VAF' }
for (c_vaf_col in vaf_cols) {
  if (!is.numeric(df_f[, c_vaf_col])) {
    log_wrap(glue("WARNING: {c_vaf_col} column is not numeric, checking if it contains percentage strings..."))
    if (all(grepl("%$", df_f[, c_vaf_col]))) {
      df_f[, c_vaf_col] <- as.numeric(gsub("%", "", df_f[, c_vaf_col])) / 100
      log_wrap(glue("Converted {c_vaf_col} percentage strings to numeric values."))
    } else {
      log_wrap(glue("ERROR: {c_vaf_col} column is not numeric and does not contain percentage strings"), log_type = 'error')
    }
  }
  df_f <- filter(df_f, !!sym(c_vaf_col) >= args$min_VAF)
  log_wrap(glue("After {c_vaf_col} >= {args$min_VAF} filtering: {nrow(df_f)} rows"))
  check_df(df_f)
}

# Filter by depth ----
depth_cols = if(args$dups_input_file) { c('DEPTH', 'DEPTH_DUP') } else {'Depth' }
for (c_depth_col in depth_cols) {
  df_f <- filter(df_f, !!sym(c_depth_col) >= args$min_depth)
  log_wrap(glue("After {c_depth_col} >= {args$min_depth} filtering: {nrow(df_f)} rows"))
  check_df(df_f)
}

# Filter intronic events (if required by user)
if (args$remove_intronic) {
  df_f <- filter(df_f, Func.refGene %in% c('exonic', 'splicing'))
  log_wrap(glue("After removing intronic mutations, left with exonic/splicing: {nrow(df_f)} rows"))
  check_df(df_f)
}

# Filter synonymous events (if required by user)
if (args$remove_synonymous) {
  df_f <- filter(df_f, ExonicFunc.refGene != "synonymous SNV")
  log_wrap(glue("After removing synonymous mutations: {nrow(df_f)} rows"))
  check_df(df_f)
}

# Annotate mutations
browser()
df_f <- annotate_mutations_df(input_df = df_f, 
                              listed_mutations_fn = args$listed_mutations_file_name,
                              mane_ifn = "../data/MANE.GRCh38.v1.4.summary.txt.gz",
                              include_non_syn_by_loci = args$include_non_syn_by_loci)
log_wrap("Annotated mutations with listed/excluded mutations and MANE information")

# Get splicing events (if required by user)
# if (args$get_splicing) {
#   splicing_df <- filter(df_f, grepl("splicing", Func.refGene, ignore.case = TRUE))
#   log_wrap(glue("Splicing mutations found: {nrow(splicing_df)} rows"))
# }

# Get stop gain/loss events (if required by user)
# if (args$get_stop_mutations) {
#   stop_df <- filter(df_f, grepl("stopgain|stoploss", ExonicFunc.refGene, ignore.case = TRUE))
#   log_wrap(glue("Stop gain/loss mutations found: {nrow(stop_df)} rows"))
# } 

# Filter excluded mutations
n_excluded <- sum(df_f$Has.excluded)
df_f <- filter(df_f, !Has.excluded)
log_wrap(glue("After removing {n_excluded} mutations in the exclusion list left with: {nrow(df_f)} rows"))
check_df(df_f)

# fwrite(df_f, gsub(".tsv", "_debug.tsv", output_fn), sep = "\t", na = "NA", quote = FALSE)

# Get mean VAF across samples
if (args$dups_input_file) {
  df_f <- df_f %>%
    rowwise() %>%
    mutate(avg_VAF = mean(c_across(c(VAF, VAF_DUP)), na.rm = TRUE),
           mean_depth = mean(c_across(c(DEPTH, DEPTH_DUP))))
  log_wrap("Calculated average VAF and DEPTH across duplicate samples")
  stopifnot(args$min_supporting_samples <= 2)
} else {
  # Assuming single sample per donor
  df_f <- mutate(df_f, avg_VAF = VAF, mean_depth = Depth)
  # sample_specific_columns <- unlist(strsplit(args$sample_specific_columns, split = ";"))
  # stopifnot(length(sample_specific_columns) > 0 && 
  #             all(sample_specific_columns %in% colnames(df_f)) && all(c('VAF', 'Depth') %in% sample_specific_columns))
  # df_f <- group_by(df_f, across(-all_of(sample_specific_columns))) %>%
  #   summarise(avg_VAF = mean(VAF, na.rm = TRUE),
  #             n_samples = n(), 
  #             mean_depth = mean(Depth, na.rm=T)) %>%
  #   ungroup() %>% 
  #   filter(n_samples >= args$min_supporting_samples)
  log_wrap("Calculated average VAF and Depth across duplicate samples")
}
# Filter listed mutations with high enough mean_VAF
f_min_avg_VAF = df_f$avg_VAF >= args$min_avg_VAF
if (args$min_avg_VAF == -1) {
  f_min_avg_VAF = df_f$avg_VAF >= (args$min_depth / 10) / (args$min_depth + df_f$mean_depth)
} 

df_ff <- df_f[df_f$Has.listed & f_min_avg_VAF, ]
log_wrap(glue("Listed mutations found: {nrow(df_ff)} rows"))
df_ff_suspected_common <- mutate(df_ff, numeric_gnomad41_AF = as.numeric(gnomad41_genome_AF)) %>%
  filter(!is.na(numeric_gnomad41_AF) & numeric_gnomad41_AF > args$max_gnomad41_AF)
if (nrow(df_ff_suspected_common) > 0) {
  sus_info = apply(df_ff_suspected_common[, c('Gene.refGene', 'Func.refGene', 'ExonicFunc.refGene', 'Selected.AA.change', 'avg_VAF')], 1, paste0, collapse=';')
  sus_joined = paste0(sus_info, collapse=',')
  log_wrap(glue("Found {nrow(df_ff)} listed mutations with population AF above {args.max_gnomad41_AF} cutoff: {sus_joined}"))
}

df_f_rest <- df_f[!df_f$Has.listed | !f_min_avg_VAF, ]

# Filter rest by gnomad41 AF
df_f_rest <- mutate(df_f_rest, numeric_gnomad41_AF = as.numeric(gnomad41_genome_AF)) %>%
  filter(is.na(numeric_gnomad41_AF) | numeric_gnomad41_AF <= args$max_gnomad41_AF)
log_wrap(glue("Non-listed mutations after gnomad41_AF <= {args$max_gnomad41_AF} filtering: {nrow(df_f_rest)} rows"))
check_df(df_f)

log_wrap(glue("Non-listed mutations yet survived filters so far found: {nrow(df_f_rest)} rows"))

# Write output 
fwrite(df_ff, output_fn, sep = "\t", na = "NA", quote = FALSE)
fwrite(df_f_rest, rest_fn, sep = "\t", na = "NA", quote = FALSE)

