# Description: Utility functions for mutation annotation and filtering
library(tidyverse)
library(data.table)
library(glue)
library(readxl)

# Annotation of transcripts matching mutations. Expecting annovar output (with the AAChange.refGene, Func.refGene and ExonicFunc.refGene columns)
#
# Mutations are tested by transcript (each might have several) against a curated list of listed/excluded mutations,
# currently based on Vlasschaert et al., Blood 2023 paper, but expected to be updated over time.
#
# Listed mutation table should have these columns:
#   - Gene
#   - Accession: NM_* transcript
#   - Mutation_types: For (Exonic)Func.refGene filtering, multiple values sep by ';'. Possible values: frameshift;nonsense;splicing;missense;nonframeshift_deletion;nonframeshift_insertion
#   - Restrict_to_exons: Further filtering by exons, optional, multiple exons sep by ';'
#   -	Restrict_to_AA_range: Further filtering by protein AA region, optional, a START-END pair per range, multiple regions sep by ';' 
#   - Restrict_to_DNA_range: similar to AA range restriction but by DNA gene coordinates
#   - Restrict_missense_to_AA_change: Further filtering of missense mutations by specific AA changes, optional, multiple changes sep by ';'
#
# Excluded mutations table should have columns Gene and Variant (same format as AAChange.refGene)
#
# Selected transcript+mutation if they appear in the list, if not, selecting the one with the top MANE status (if MANE file provided by user), otherwise the last transcript in the 'aa_change_col'.

# Annotation adds these columns to the input data frame: 
# 1. Has.protein.coding: Has a valid transcript (NM_*)
# 2. Top.MANE.status: Highest MANE status among transcripts (MANE plus clinical, MANE select, None), if MANE file supplied by user
# 3. Has.listed.: any of the transcripts matches a transcript+mutation in the 'included' list
# 4. Has.excluded: any of the transcripts matches a transcript+mutation in the 'excluded' list
# 5. Selected.transcript: see above
# 5. Selected.AA.change: see above
# 6. Last.AA.Change: AA change of the last transcript in 'aa_change_col' (to check for backward compatibility/sanity if it is always the MANE select transcript)
#
# Output: Beside the 3 columns above, a new column named 'aa_change_col'.Annotated with information on transcripts (comma delimited 
# multiple entries, same order as transcripts in the input 'aa_change_col' column): 
#    (Gene : transcript : exon : cDNA change str : protein change str : AA change : AA position : MANE status : Listed mutation : Excluded mutation)
annotate_mutations_df_orig <- function(input_df, 
                                  listed_mutations_fn = "../data/Vlasschaert_Blood_2023_data.xlsx",
                                  mane_ifn = NULL, 
                                  mane_status_priorities = c('MANE Plus Clinical', 'MANE Select', 'None'),
                                  ignore_transcript_id_version = TRUE,
                                  mut_info_col = "AAChange.refGene",
                                  non_exonic_mut_info_col = 'GeneDetail.refGene',
                                  merged_mut_info_col = paste0(mut_info_col, '.Annotated'),
                                  func_col = "Func.refGene",
                                  exonic_func_col = "ExonicFunc.refGene",
                                  mut_info_sep = ",",
                                  overwrite = T) {
  
  # Check if input_df already contains annotation columns
  if (any(colnames(input_df) %in% c('Has.protein.coding', 'Top.MANE.status', 'Has.listed', 'Has.excluded', 'Selected.transcript', 'Selected.AA.change', 'Last.AA.Change', merged_mut_info_col)) & !overwrite) {
    stop("Input data frame already contains some of the output annotation columns, set overwrite = TRUE to overwrite them.")
  }
  
  # Check if mut_info_col exists in input_df
  if (!(mut_info_col %in% colnames(input_df))) {
    stop(glue("Column {mut_info_col} not found in input data frame."))
  }
  if (!(non_exonic_mut_info_col %in% colnames(input_df))) {
    stop(glue("Column {non_exonic_mut_info_col} not found in input data frame."))
  }
  
  # Get AA change of last occurrence of "p." in mut_info_col
  input_df <- mutate(input_df, Last.AA.Change := gsub(".*:p.", "", !!sym(mut_info_col)))
  
  # Store original columns for consolidated output data frame                   
  orig_columns <- setdiff(colnames(input_df), mut_info_col)
                     
  # Separate mut_info_col column into multiple rows using ',' as a separator
  df_exonic <- filter(input_df, !!sym(mut_info_col) != '.') %>%
    tidyr::separate_longer_delim(cols = all_of(mut_info_col), delim = mut_info_sep) %>%
    separate(all_of(mut_info_col), into = paste0("sub__", c("gene", "transcript", "exon", "cDNA_change", "protein_change")), sep = ":", remove=F) %>%
    mutate(sub__protein_change = gsub("^p.", "", sub__protein_change),
           sub__aa_position = as.integer(gsub("^[A-Za-z]+", "", gsub("[A-Za-z]+$", "", sub__protein_change))),
           sub__cDNA_change = gsub("^c.", "", sub__cDNA_change), 
           sub__cDNA_start = as.integer(gsub("_[0-9]+", "", gsub("[A-Za-z]", "", sub__cDNA_change))),
           sub__cDNA_end = as.integer(gsub("[0-9]+_", "", gsub("[A-Za-z]", "", sub__cDNA_change))),
           Has.protein.coding = grepl("NM_", sub__transcript),
           mutation_func = ifelse(Func.refGene == 'splicing', 'splicing', ExonicFunc.refGene)
           ) %>%
    tidyr::separate_rows(mutation_func, sep = ";") %>%
    mutate(row_num = row_number())
  
  df_non_exonic <- filter(input_df, func_col == 'splicing') %>%
    tidyr::separate_longer_delim(cols = all_of(non_exonic_mut_info_col), delim = mut_info_sep) %>%
    separate(all_of(non_exonic_mut_info_col), into = paste0("sub__", c("transcript", "exon", "cDNA_change")), sep = ":", remove=F) %>%
    mutate(sub__cDNA_change = gsub("^c.", "", sub__cDNA_change), 
           sub__cDNA_start = as.integer(gsub("_[0-9]+", "", gsub("[A-Za-z]", "", sub__cDNA_change))),
           sub__cDNA_end = as.integer(gsub("[0-9]+_", "", gsub("[A-Za-z]", "", sub__cDNA_change))),
           mutation_func = 'splicing',
           sub__gene = 'None', sub__protein_change = 'None', sub__aa_position = -1, Has.protein.coding = F
    ) %>%
    tidyr::separate_rows(mutation_func, sep = ";") %>%
    mutate(row_num = nrow(df_exonic) + row_number())
  df <- rbind(df_exonic, df_non_exonic)
  
  df$transcript_for_join <- if(ignore_transcript_id_version) {
      gsub("\\.[0-9]+$", "", df$sub__transcript) } else
      { df$sub__transcript }

  # Add MANE select information (if mane_ifn provided)
  if (file.exists(mane_ifn)) {
    mane = fread(mane_ifn)
    mane$transcript_for_join <- if(ignore_transcript_id_version) {
      gsub("\\.[0-9]+$", "", mane$RefSeq_nuc) } else
      { mane$RefSeq_nuc }
    
    # Join with MANE data
    df <- left_join(df, dplyr::select(mane, transcript_for_join, MANE_status)) %>%
      mutate(MANE_status = ifelse(is.na(MANE_status), 'None', MANE_status),
             MANE_status = factor(MANE_status, levels = mane_status_priorities, ordered = TRUE))
  } else {
    df$MANE_status <- NA
  }
  
  # Parse listed mutations information
  listed_mutations <- read_excel(path = listed_mutations_fn, sheet = 'Include') %>%
    filter(Include == 'Yes') %>%
    dplyr::select(Gene, Accession, Mutation_types, Restrict_to_exons, Restrict_to_AA_range, Restrict_to_DNA_range, Restrict_missense_to_AA_change) %>%
    tidyr::separate_rows(Mutation_types, sep = ";")
  
  stopifnot(all(unique(listed_mutations$Mutation_types) %in% 
                  c('frameshift', 'nonsense', 'splicing', 'missense', 'nonframeshift deletion', 'nonframeshift insertion')))
  
  # Inflate table by different mutation types and restricting values
  mut_dict = c('frameshift' = 'frameshift insertion;frameshift deletion', 
               'nonsense' = 'stopgain;stoploss', 
               'splicing' = 'splicing', 
               'missense' = 'nonsynonymous SNV', 
               'nonframeshift_deletion' = 'nonframeshift_deletion', 
               'nonframeshift_insertion' = 'nonframeshift_insertion')
  
  listed_mutations <- mutate(listed_mutations, mutation_func = mut_dict[Mutation_types]) %>%
    tidyr::separate_rows(mutation_func, sep = ";") %>%
    tidyr::separate_rows(Restrict_to_exons, sep = ";") %>%
    tidyr::separate_rows(Restrict_to_AA_range, sep = ";") %>%
    tidyr::separate_rows(Restrict_to_DNA_range, sep = ";") %>%
    tidyr::separate_rows(Restrict_missense_to_AA_change, sep = ";") %>%
    mutate(AA_range_start = gsub("-.*", "", Restrict_to_AA_range),
           AA_range_end   = gsub(".*-", "", Restrict_to_AA_range),
           DNA_range_start = gsub("-.*", "", Restrict_to_DNA_range),
           DNA_range_end   = gsub(".*-", "", Restrict_to_DNA_range))
  
  # Join with listed mutations data and filter
  message("I'm sane!!!")
  browser()
  listed_df <- inner_join(df, listed_mutations, 
              by = c(transcript_for_join = 'Accession',
                     mutation_func = 'mutation_func'),
              relationship = "many-to-many") %>%
    filter(is.na(Restrict_to_exons) | sub__exon == paste0('exon', Restrict_to_exons)) %>%
    filter(is.na(Restrict_to_AA_range) | (sub__aa_position >= AA_range_start & sub__aa_position <= AA_range_end)) %>%
    filter(is.na(Restrict_to_DNA_range) | (sub__cDNA_start >= DNA_range_start & sub__cDNA_end <= DNA_range_end) ) %>%
    filter(is.na(Restrict_missense_to_AA_change) | (mutation_func != 'nonsynonymous SNV') | (mutation_func == 'nonsynonymous SNV' & sub__protein_change == Restrict_missense_to_AA_change))
  
  browser()
  
  rest_df <- anti_join(df, listed_df, by = c("row_num")) %>% 
    mutate()
  
  df <- bind_rows(
    mutate(dplyr::select(listed_df, colnames(df), row_num), is_listed = 'Listed'),
    mutate(dplyr::select(rest_df,   colnames(df), row_num), is_listed = 'Not listed')
  ) %>% 
    arrange(row_num)
    
  
  # Mark excluded mutations
  excluded_mutations <- read_excel(path = listed_mutations_fn, sheet = 'Exclude') %>%
    dplyr::select('Gene', 'Variant')
  colnames(excluded_mutations) = c('sub__gene', 'gene_info_col')
  
  df <- df %>%
    left_join(mutate(excluded_mutations, is_excluded = T)) %>%
    mutate(is_excluded = ifelse(is.na(is_excluded), 'Not excluded', 'Excluded'))
  
  # Create consolidated data frame, making sure to keep original order of rows
  output_df <- mutate(df, merged__sub = glue("{sub__gene}:{sub__transcript}:{sub__exon}:{sub__cDNA_change}:{sub__protein_change}:{sub__aa_position}:{MANE_status}:{is_listed}:{is_excluded}")) %>%
    group_by(across(all_of(orig_columns))) %>%
    summarise(Has.protein.coding = any(Has.protein.coding),
              Top.MANE.status = mane_status_priorities[min(as.numeric(MANE_status))],
              Has.listed = any(is_listed == 'Listed'),
              Has.excluded = any(is_excluded == 'Excluded'),
              Selected.transcript = ifelse(Has.listed, 
                                           sub__transcript[which(is_listed == 'Listed')[1]], 
                                           ifelse(Top.MANE.status == 'None', 
                                                  sub__transcript[n()], 
                                                  sub__transcript[which(MANE_status == Top.MANE.status)[1]])),
              Selected.AA.change = ifelse(Has.listed, 
                                          sub__protein_change[which(is_listed == 'Listed')[1]], 
                                          ifelse(Top.MANE.status == 'None', 
                                                 sub__protein_change[n()], 
                                                 sub__protein_change[which(MANE_status == Top.MANE.status)[1]])),
              !!(merged_mut_info_col) := paste0(merged__sub, collapse = mut_info_sep),
              !!(mut_info_col) := paste0(!!sym(mut_info_col), collapse = mut_info_sep)) %>%
    ungroup %>%
    right_join(dplyr::select(input_df, all_of(orig_columns)), by = orig_columns)

  return(output_df)
}

# Annotation adds these columns to the input data frame: 
# 1. Has.protein.coding: Has a valid transcript (NM_*)
# 2. Top.MANE.status: Highest MANE status among transcripts (MANE plus clinical, MANE select, None), if MANE file supplied by user
# 3. Has.listed.: any of the transcripts matches a transcript+mutation in the 'included' list
# 4. Has.excluded: any of the transcripts matches a transcript+mutation in the 'excluded' list
# 5. Selected.transcript: see above
# 5. Selected.AA.change: see above
# 6. Last.AA.Change: AA change of the last transcript in 'aa_change_col' (to check for backward compatibility/sanity if it is always the MANE select transcript)
#
# Output: Beside the 3 columns above, a new column named 'aa_change_col'.Annotated with information on transcripts (comma delimited 
# multiple entries, same order as transcripts in the input 'aa_change_col' column): 
#    (Gene : transcript : exon : cDNA change str : protein change str : AA change : AA position : MANE status : Listed mutation : Excluded mutation)
annotate_mutations_df <- function(input_df, 
                                  listed_mutations_fn = "../data/Vlasschaert_Blood_2023_data.xlsx",
                                  mane_ifn = NULL, 
                                  mane_status_priorities = c('MANE Plus Clinical', 'MANE Select', 'None'),
                                  ignore_transcript_id_version = TRUE,
                                  mut_info_col = "AAChange.refGene",
                                  non_exonic_mut_info_col = 'GeneDetail.refGene',
                                  merged_mut_info_col = paste0(mut_info_col, '.Annotated'),
                                  func_col = "Func.refGene",
                                  exonic_func_col = "ExonicFunc.refGene",
                                  mut_info_sep = ",",
                                  non_exonic_mut_info_sep = ';',
                                  include_non_syn_by_loci = T,
                                  overwrite = T) {
  
  # Check if input_df already contains annotation columns
  if (any(colnames(input_df) %in% c('Has.protein.coding', 'Top.MANE.status', 'Has.listed', 'Has.excluded', 'Selected.transcript', 'Selected.AA.change', 'Last.AA.Change', merged_mut_info_col)) & !overwrite) {
    stop("Input data frame already contains some of the output annotation columns, set overwrite = TRUE to overwrite them.")
  }
  
  # Check if mut_info_col exists in input_df
  if (!(mut_info_col %in% colnames(input_df))) {
    stop(glue("Column {mut_info_col} not found in input data frame."))
  }
  if (!(non_exonic_mut_info_col %in% colnames(input_df))) {
    stop(glue("Column {non_exonic_mut_info_col} not found in input data frame."))
  }
  
  # Get AA change of last occurrence of "p." in mut_info_col
  input_df <- mutate(input_df, Last.AA.Change := gsub(".*:p.", "", !!sym(mut_info_col)))
  
  # Store original columns for consolidated output data frame                   
  orig_columns <- setdiff(colnames(input_df), mut_info_col)
  
  # Separate rows to mutliple transcripts (for exonic and splicing) column into multiple rows using ',' as a separator
  df <- mutate(input_df, non_exonic_mut_info_col_fixed = gsub(non_exonic_mut_info_sep, mut_info_sep, !!sym(non_exonic_mut_info_col))) %>%
    mutate(transcript_info_col = case_when(
    !!sym(func_col) == 'splicing' ~ non_exonic_mut_info_col_fixed,
    !!sym(func_col) == 'exonic' ~ !!sym(mut_info_col),
    TRUE ~ '.'
  ))
  
  df <- df %>%
    tidyr::separate_longer_delim(cols = transcript_info_col, delim = mut_info_sep) %>%
    mutate(transcript_info_col = if_else(!!sym(func_col) == 'splicing', paste('None', transcript_info_col, 'p.A-1A', sep=":"), transcript_info_col)) %>%
    separate(transcript_info_col, into = paste0("sub__", c("gene", "transcript", "exon", "cDNA_change", "protein_change")), sep = ":", remove=F) %>%
    mutate(sub__protein_change = gsub("^p.", "", sub__protein_change),
           # Capture the digit run immediately after the leading AA letter(s), e.g. "W288Cfs*12" -> "288".
           # (Stripping a trailing letter run instead would miss this: frameshift suffixes like "fs*12"
           # end in digits, not letters, so the AA position would never be extracted.)
           sub__aa_position = as.integer(stringr::str_match(sub__protein_change, "^[A-Za-z]+(\\d+)")[, 2]),
           sub__cDNA_change = gsub("^c.", "", sub__cDNA_change),
           sub__cDNA_start = as.integer(gsub("_[0-9]+", "", gsub("[A-Za-z]", "", sub__cDNA_change))),
           sub__cDNA_end = as.integer(gsub("[0-9]+_", "", gsub("[A-Za-z]", "", sub__cDNA_change))),
           Has.protein.coding = grepl("NM_", sub__transcript),
           mutation_func = ifelse(!!sym(func_col) == 'splicing', 'splicing', !!sym(exonic_func_col))
    ) %>%
    tidyr::separate_rows(mutation_func, sep = ";") %>%
    mutate(row_num = row_number())

  df$transcript_for_join <- if(ignore_transcript_id_version) {
    gsub("\\.[0-9]+$", "", df$sub__transcript) } else
    { df$sub__transcript }
  
  # Add MANE select information (if mane_ifn provided)
  if (file.exists(mane_ifn)) {
    mane = fread(mane_ifn)
    mane$transcript_for_join <- if(ignore_transcript_id_version) {
      gsub("\\.[0-9]+$", "", mane$RefSeq_nuc) } else
      { mane$RefSeq_nuc }
    
    # Join with MANE data
    df <- left_join(df, dplyr::select(mane, transcript_for_join, MANE_status)) %>%
      mutate(MANE_status = ifelse(is.na(MANE_status), 'None', MANE_status),
             MANE_status = factor(MANE_status, levels = mane_status_priorities, ordered = TRUE))
  } else {
    df$MANE_status <- NA
  }
  
  # Parse listed mutations information
  listed_mutations <- read_excel(path = listed_mutations_fn, sheet = 'Include') %>%
    filter(Include == 'Yes') %>%
    dplyr::select(Gene, Accession, Mutation_types, Restrict_to_exons, Restrict_to_AA_range, Restrict_to_DNA_range, Restrict_missense_to_AA_change) %>%
    tidyr::separate_rows(Mutation_types, sep = ";")
  
  stopifnot(all(unique(listed_mutations$Mutation_types) %in% 
                  c('frameshift', 'nonsense', 'splicing', 'missense', 'nonframeshift deletion', 'nonframeshift insertion')))
  
  # Inflate table by different mutation types and restricting values
  mut_dict = c('frameshift' = 'frameshift insertion;frameshift deletion', 
               'nonsense' = 'stopgain;stoploss', 
               'splicing' = 'splicing', 
               'missense' = 'nonsynonymous SNV', 
               'nonframeshift_deletion' = 'nonframeshift_deletion', 
               'nonframeshift_insertion' = 'nonframeshift_insertion')
  
  listed_mutations <- mutate(listed_mutations, mutation_func = mut_dict[Mutation_types]) %>%
    tidyr::separate_rows(mutation_func, sep = ";") %>%
    tidyr::separate_rows(Restrict_to_exons, sep = ";") %>%
    tidyr::separate_rows(Restrict_to_AA_range, sep = ";") %>%
    tidyr::separate_rows(Restrict_to_DNA_range, sep = ";") %>%
    tidyr::separate_rows(Restrict_missense_to_AA_change, sep = ";") %>%
    mutate(AA_range_start = as.numeric(gsub("-.*", "", Restrict_to_AA_range)),
           AA_range_end   = as.numeric(gsub(".*-", "", Restrict_to_AA_range)),
           DNA_range_start = as.numeric(gsub("-.*", "", Restrict_to_DNA_range)),
           DNA_range_end   = as.numeric(gsub(".*-", "", Restrict_to_DNA_range)),
           AA_change_start = as.numeric(gsub("[A-Z]", "", Restrict_missense_to_AA_change)))
  
  # Join with listed mutations data and filter
  listed_df <- inner_join(df, listed_mutations, 
                          by = c(transcript_for_join = 'Accession',
                                 mutation_func = 'mutation_func'),
                          relationship = "many-to-many") %>%
    filter(is.na(Restrict_to_exons) | sub__exon == paste0('exon', Restrict_to_exons)) %>%
    filter(is.na(Restrict_to_AA_range) | (sub__aa_position >= AA_range_start & sub__aa_position <= AA_range_end)) %>%
    filter(is.na(Restrict_to_DNA_range) | (sub__cDNA_start >= DNA_range_start & sub__cDNA_end <= DNA_range_end) ) %>%
    filter(is.na(Restrict_missense_to_AA_change) | 
             (mutation_func != 'nonsynonymous SNV') | 
             (mutation_func == 'nonsynonymous SNV' & (sub__protein_change == Restrict_missense_to_AA_change)))
  
  listed_df <- inner_join(df, listed_mutations, 
                          by = c(transcript_for_join = 'Accession',
                                 mutation_func = 'mutation_func'),
                          relationship = "many-to-many") %>%
    filter(is.na(Restrict_to_exons) | sub__exon == paste0('exon', Restrict_to_exons)) %>%
    filter(is.na(Restrict_to_AA_range) | (sub__aa_position >= AA_range_start & sub__aa_position <= AA_range_end)) %>%
    filter(is.na(Restrict_to_DNA_range) | (sub__cDNA_start >= DNA_range_start & sub__cDNA_end <= DNA_range_end) ) %>%
    filter(is.na(Restrict_missense_to_AA_change) | 
             (mutation_func != 'nonsynonymous SNV') | 
             (mutation_func == 'nonsynonymous SNV' & ((sub__protein_change == Restrict_missense_to_AA_change) | (include_non_syn_by_loci & sub__aa_position == AA_change_start))))
  
  
  rest_df <- anti_join(df, listed_df, by = c("row_num")) 

  df <- bind_rows(
    mutate(dplyr::select(listed_df, colnames(df), row_num), is_listed = 'Listed'),
    mutate(dplyr::select(rest_df,   colnames(df), row_num), is_listed = 'Not listed')
  ) %>% 
    arrange(row_num)
  
  # Mark excluded mutations (for exonic AA changes drop the protein part which may have different formats)
  df <- df %>%
    mutate(AAChange.refGene.for_join = if_else(mutation_func == 'splicing', '.', transcript_info_col), 
           GeneDetail.refGene.for_join = if_else(mutation_func == 'splicing', glue("{sub__transcript}:{sub__exon}:c.{sub__cDNA_change}"), '.')) %>%
    mutate(AAChange.refGene.for_join = gsub(":p\\..+", "", AAChange.refGene.for_join))
    
   
  excluded_mutations <- read_excel(path = listed_mutations_fn, sheet = 'Exclude') %>%
    dplyr::select('Gene.refGene', 'AAChange.refGene', 'GeneDetail.refGene') %>%
    mutate(AAChange.refGene = gsub(":p\\..+", "", AAChange.refGene))
  
  colnames(excluded_mutations) = c('Gene.refGene', 'AAChange.refGene.for_join', 'GeneDetail.refGene.for_join')
  
  excluded_mutations <- excluded_mutations %>%
    tidyr::separate_rows('AAChange.refGene.for_join', sep = mut_info_sep) %>%
    tidyr::separate_rows('GeneDetail.refGene.for_join', sep = non_exonic_mut_info_sep)
  
  df <- df %>%
    left_join(mutate(excluded_mutations, is_excluded = T)) %>%
    mutate(is_excluded = ifelse(is.na(is_excluded), 'Not excluded', 'Excluded'))
  
  
  # Create consolidated data frame, making sure to keep original order of rows
  output_df <- mutate(df, merged__sub = glue("{sub__gene}:{sub__transcript}:{sub__exon}:{sub__cDNA_change}:{sub__protein_change}:{sub__aa_position}:{MANE_status}:{is_listed}:{is_excluded}")) %>%
    group_by(across(all_of(orig_columns))) %>%
    summarise(Has.protein.coding = any(Has.protein.coding),
              Top.MANE.status = mane_status_priorities[min(as.numeric(MANE_status))],
              Has.listed = any(is_listed == 'Listed'),
              Has.excluded = any(is_excluded == 'Excluded'),
              Selected.transcript = ifelse(Has.listed, 
                                           sub__transcript[which(is_listed == 'Listed')[1]], 
                                           ifelse(Top.MANE.status == 'None', 
                                                  sub__transcript[n()], 
                                                  sub__transcript[which(MANE_status == Top.MANE.status)[1]])),
              Selected.AA.change = ifelse(Has.listed, 
                                          sub__protein_change[which(is_listed == 'Listed')[1]], 
                                          ifelse(Top.MANE.status == 'None', 
                                                 sub__protein_change[n()], 
                                                 sub__protein_change[which(MANE_status == Top.MANE.status)[1]])),
              !!(merged_mut_info_col) := paste0(merged__sub, collapse = mut_info_sep),
              !!(mut_info_col) := paste0(!!sym(mut_info_col), collapse = mut_info_sep)) %>%
    ungroup %>%
    right_join(dplyr::select(input_df, all_of(orig_columns)), by = orig_columns)
  
  return(output_df)
}
# Function to check if all target mutations are covered by the panel
check_panel_covers_targets <- function()
{  
  # Placeholder for future implementation
  
}

# old -----

# Annotation of transcripts matching mutations. Expecting transcript matches in the 'aa_change_col' (e.g. AAChange.refSeq from annovar) 
# with the following structure: gene, transcript, exon, cDNA_change, protein_change.
#
# Annotation adds these columns to the input data frame: 
# 1. Has.protein.coding: Has a valid transcript (NM_*)
# 2. Top.MANE.status: Highest MANE status among transcripts (MANE plus clinical, MANE select, None), if MANE file supplied by user
# 3. Selected.AA.change: AA change of the transcript with the highest MANE status (if multiple transcripts have the same MANE status, the first one is selected) or the last AA change in the 'aa_change_col' if no MANE transcript exist
# 4. Has.hotspot: Any of the transcripts matches a hotspot AA change (if hotspots supplied by user)
# 5. In.hotspot.region: mutation falls within any of the hotspot gene regions (if supplied by user)
# 6. Last.AA.Change: AA change of the last transcript in 'aa_change_col' (to check for backward compatibility/sanity if it is always the MANE select transcript)
#
# Output: Beside the 3 columns above, a new column named 'aa_change_col'.Annotated with information on transcripts (comma delimited 
# multiple entries, same order as transcripts in the input 'aa_change_col' column): 
#    (Gene : transcript : exon : cDNA change str : protein change str : AA change : AA position : MANE status : Is hotspot : Matched cDNA hotspot region)
annotate_mutations_df_based_on_LS_funcs <- function(input_df, 
                                                    hotspots_df = NULL, 
                                                    hotspots_gene_col = 'Gene',
                                                    hotspots_aa_change_col = 'AA_change',
                                                    gene_hotspot_regions_df = NULL,
                                                    regions_gene_col = 'Gene',
                                                    regions_start_col = 'Region_Start',
                                                    regions_end_col = 'Region_End',
                                                    mane_ifn = NULL, 
                                                    mane_status_priorities = c('MANE Plus Clinical', 'MANE Select', 'None'),
                                                    ignore_transcript_id_version = TRUE,
                                                    mut_info_col = "AAChange.refGene", 
                                                    merged_mut_info_col = paste0(mut_info_col, '.Annotated'),
                                                    mut_info_sep = ",",
                                                    overwrite = T) {
  # Check if input_df already contains annotation columns
  if (any(colnames(input_df) %in% c('Has.protein.coding', 'Top.MANE.status', 'Has.hotspot', 'Last.AA.Change', merged_mut_info_col)) & !overwrite) {
    stop("Input data frame already contains annotation columns, set overwrite = TRUE to overwrite them.")
  }
  
  # Check if mut_info_col exists in input_df
  if (!(mut_info_col %in% colnames(input_df))) {
    stop(glue("Column {mut_info_col} not found in input data frame."))
  }
  
  # Get AA change of last occurancee of "p." in mut_info_col
  input_df <- mutate(input_df, Last.AA.Change := gsub(".*:p.", "", !!sym(mut_info_col)))
  
  # Store original columns for consolidated output data frame                   
  orig_columns <- setdiff(colnames(input_df), mut_info_col)
  
  # Separate mut_info_col column into multiple rows using ',' as a separator
  df <- tidyr::separate_longer_delim(input_df, cols = all_of(mut_info_col), delim = mut_info_sep) %>%
    separate(all_of(mut_info_col), into = paste0("sub__", c("gene", "transcript", "exon", "cDNA_change", "protein_change")), sep = ":", remove=F) %>%
    mutate(sub__protein_change = gsub("^p.", "", sub__protein_change),
           sub__aa_position = as.integer(gsub("^[A-Za-z]+", "", gsub("[A-Za-z]+$", "", sub__protein_change))),
           sub__cDNA_change = gsub("^c.", "", sub__cDNA_change), 
           sub__cDNA_start = as.integer(gsub("_[0-9]+", "", gsub("[A-Za-z]", "", sub__cDNA_change))),
           sub__cDNA_end = as.integer(gsub("[0-9]+_", "", gsub("[A-Za-z]", "", sub__cDNA_change))),
           Has.protein.coding = grepl("NM_", sub__transcript))
  
  df$transcript_for_join <- if(ignore_transcript_id_version) {
    gsub("\\.[0-9]+$", "", df$sub__transcript) } else
    { df$sub__transcript }
  
  # Add MANE select information (if mane_ifn provided)
  if (file.exists(mane_ifn)) {
    mane = fread(mane_ifn)
    mane$transcript_for_join <- if(ignore_transcript_id_version) {
      gsub("\\.[0-9]+$", "", mane$RefSeq_nuc) } else
      { mane$RefSeq_nuc }
    
    
    # Join with MANE data
    df <- left_join(df, dplyr::select(mane, transcript_for_join, MANE_status)) %>%
      mutate(MANE_status = ifelse(is.na(MANE_status), 'None', MANE_status),
             MANE_status = factor(MANE_status, levels = mane_status_priorities, ordered = TRUE))
  } else {
    df$MANE_status <- NA
  }
  
  # Add hotspot information (if provided)
  if (!is.null(hotspots_df)) {
    hotspots_df <- dplyr::select(hotspots_df, all_of(c(hotspots_gene_col, hotspots_aa_change_col))) %>%
      distinct() %>%
      mutate(is_hotspot = T) 
    df <- left_join(df, hotspots_df, by=c(sub__gene = hotspots_gene_col, sub__protein_change = hotspots_aa_change_col)) %>%
      mutate(is_hotspot = ifelse(is.na(is_hotspot), 'Not Hotspot', 'Hotspot')) 
  } else {
    df$is_hotspot <- 'NA'
  }
  
  
  # Add mutations in specific gene regions (if provided)
  if (!is.null(gene_hotspot_regions_df)) {
    matched_df <- inner_join(df, gene_hotspot_regions_df, by=c(sub__gene = regions_gene_col)) %>%
      filter((!is.na(sub__cDNA_start)  & sub__cDNA_start >= !!sym(regions_start_col) & sub__cDNA_start <= !!sym(regions_end_col) ) |
               (!is.na(sub__cDNA_end)  & sub__cDNA_end >= !!sym(regions_start_col)  & sub__cDNA_end <= !!sym(regions_end_col))) %>%
      mutate(matched_cDNA_region = paste(!!sym(regions_start_col) , !!sym(regions_end_col), sep='-')) %>%
      dplyr::select(-!!sym(regions_start_col) , -!!sym(regions_end_col))
    
    df <- left_join(df, matched_df) %>%
      mutate(matched_cDNA_region = ifelse(is.na(matched_cDNA_region), 'None', matched_cDNA_region))
  } else {
    df$matched_cDNA_region = 'NA'
  }
  
  # Create consolidated data frame, making sure to keep original order of rows
  output_df <- mutate(df, merged__sub = glue("{sub__gene}:{sub__transcript}:{sub__exon}:{sub__cDNA_change}:{sub__protein_change}:{sub__aa_position}:{MANE_status}:{is_hotspot}:{matched_cDNA_region}")) %>%
    group_by(across(all_of(orig_columns))) %>%
    summarise(Has.protein.coding = any(!is.na(sub__transcript)),
              Top.MANE.status = mane_status_priorities[min(as.numeric(MANE_status))],
              Selected.AA.change = ifelse(Top.MANE.status == 'None', Last.AA.Change, sub__protein_change[which(MANE_status == Top.MANE.status)[1]]),
              Has.hotspot = any(is_hotspot == 'Hotspot'),
              In.hotspot.region = any(!(matched_cDNA_region %in% c('NA', 'None'))),
              !!(merged_mut_info_col) := paste0(merged__sub, collapse = mut_info_sep),
              !!(mut_info_col) := paste0(!!sym(mut_info_col), collapse = mut_info_sep),
              !!(non_exonic_mut_info_col) := paste0(!!sym(non_exonic_mut_info_col), collapse = mut_info_sep)) %>%
    ungroup %>%
    right_join(dplyr::select(input_df, all_of(orig_columns)), by = orig_columns)
  
  return(output_df)
}

# Function: extract amino acid position from the last occurrence of "p." in AAChange.refGene  
extract_aa_position <- function(aa_change) {  
  if (is.na(aa_change) || aa_change == "") {  
    return(NA)  
  }  
  
  # Find the last occurrence of "p."  
  last_p_pos <- stringr::str_locate_all(aa_change, "p\\.")[[1]]  
  if (nrow(last_p_pos) == 0) {  
    return(NA)  
  }  
  
  last_p_start <- last_p_pos[nrow(last_p_pos), 1]  
  aa_change_part <- substr(aa_change, last_p_start, nchar(aa_change))  
  
  # Extract position using regex  
  position_match <- stringr::str_match(aa_change_part, "p\\.[A-Za-z]*(\\d+)[A-Za-z_\\*]*")  
  if (is.na(position_match[1, 2])) {  
    return(NA)  
  }  
  
  return(as.numeric(position_match[1, 2]))  
}  

# Extract the specific AA change from AAChange.refGene  
extract_specific_aa_change <- function(aa_change) {  
  if (is.na(aa_change) || aa_change == "") {  
    return(NA)  
  }  
  
  # Find the last occurrence of "p."  
  last_p_pos <- stringr::str_locate_all(aa_change, "p\\.")[[1]]  
  if (nrow(last_p_pos) == 0) {  
    return(NA)  
  }  
  
  last_p_pos <- last_p_pos[nrow(last_p_pos), "start"]  
  
  # Extract everything after the last "p."  
  aa_change_str <- stringr::str_sub(aa_change, last_p_pos + 2)  
  
  # Extract the specific AA change (e.g., R882H)  
  specific_change <- stringr::str_extract(aa_change_str, "[A-Z]\\d+[A-Z\\*]")  
  
  return(specific_change)  
}  

matches_target_string <- function(aa_change, targets) {  
  if (is.na(aa_change) || length(targets) == 0) {  
    return(FALSE)  
  }  
  
  for (target in targets) {  
    if (grepl(target, aa_change, fixed = TRUE)) {  
      return(TRUE)  
    }  
  }  
  
  return(FALSE)  
}  

