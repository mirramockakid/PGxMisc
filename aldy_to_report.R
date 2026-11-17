#!/usr/bin/env Rscript
# Collapse Aldy *.out files into one star-allele table (star_alleles.csv).
# Usage: Rscript aldy_to_report.R --path ./results

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
})

report_columns <- c("Sample", "Gene", "SolutionID", "Major", "Minor")

# Return the value following --path (supports "--path dir" and "--path=dir").
parse_path_arg <- function(args) {
  inline <- grep("^--path=", args, value = TRUE)
  if (length(inline) == 1) {
    return(sub("^--path=", "", inline))
  }
  flag_index <- which(args == "--path")
  if (length(flag_index) != 1 || flag_index == length(args)) {
    stop("Usage: Rscript aldy_to_report.R --path <directory of .out files>", call. = FALSE)
  }
  args[flag_index + 1]
}

# Read one Aldy output file and return its distinct star-allele rows.
read_aldy_file <- function(file_path) {
  # Lines starting with "#" are the "#Sample..." header and the per-diplotype
  # "#Solution N: ..." comments. Data rows end with a variable number of
  # trailing fields (Code/Status are sometimes empty), which makes readr warn
  # about ragged rows, so trim each row to the first five fields before parsing.
  data_lines <- grep("^#", readLines(file_path), value = TRUE, invert = TRUE)
  first_five_fields <- sub("^(([^\t]*\t){4}[^\t]*).*$", "\\1", data_lines)

  read_tsv(
    I(first_five_fields),
    col_names = report_columns,
    col_types = cols(.default = col_character()),
    progress = FALSE
  ) |>
    distinct()
}

# Split a diplotype column into single alleles and return per-gene count and frequency.
allele_frequencies <- function(star_alleles, diplotype_column, delimiter, allele_label) {
  # Every solution of an ambiguous sample is counted, so such samples add
  # 2 alleles per solution. The denominator is the number of alleles counted
  # per gene (not 2 x samples), so frequencies within a gene sum to 1.
  # Homozygous diplotypes (e.g. *1/*1) correctly contribute two copies.
  star_alleles |>
    select(Gene, allele = all_of(diplotype_column)) |>
    separate_longer_delim(allele, delim = delimiter) |>
    count(Gene, allele, name = "Count") |>
    group_by(Gene) |>
    mutate(Frequency = round(Count / sum(Count), 4)) |>
    ungroup() |>
    arrange(Gene, desc(Count), allele) |>
    rename({{ allele_label }} := allele)
}

main <- function() {
  input_dir <- parse_path_arg(commandArgs(trailingOnly = TRUE))
  if (!dir.exists(input_dir)) {
    stop("Directory not found: ", input_dir, call. = FALSE)
  }

  aldy_files <- list.files(input_dir, pattern = "\\.out$", full.names = TRUE)
  if (length(aldy_files) == 0) {
    stop("No .out files found in: ", input_dir, call. = FALSE)
  }

  # A sample can carry several equally likely solutions; all are kept so the
  # ambiguity remains visible in the report.
  star_alleles <- lapply(aldy_files, read_aldy_file) |>
    bind_rows() |>
    distinct() |>
    arrange(Sample, Gene, as.integer(SolutionID))

  # Ambiguity is per sample *and* gene: a sample genotyped for several genes
  # should only be flagged for the gene(s) with competing solutions.
  ambiguous <- star_alleles |>
    group_by(Sample, Gene) |>
    filter(n_distinct(SolutionID) > 1) |>
    ungroup()

  allele_freq <- allele_frequencies(star_alleles, "Major", "/", "Allele")
  # Minor sub-alleles are ";"-separated and not necessarily in the same order
  # as the Major alleles; order doesn't matter for counting.
  suballele_freq <- allele_frequencies(star_alleles, "Minor", ";", "Suballele")

  write_csv(star_alleles, "star_alleles.csv")
  write_csv(allele_freq, "star_allele_freq.csv")
  write_csv(suballele_freq, "suballele_freq.csv")
  # Always written (header-only when nothing is ambiguous) so downstream
  # steps can rely on the file existing.
  write_csv(ambiguous, "ambiguous.csv")
  message(
    "Wrote ", nrow(star_alleles), " rows from ", length(aldy_files),
    " files to star_alleles.csv\n",
    "Wrote ", nrow(ambiguous), " rows (",
    nrow(distinct(ambiguous, Sample, Gene)),
    " sample/gene pairs with >1 solution) to ambiguous.csv\n",
    "Wrote ", nrow(allele_freq), " alleles to star_allele_freq.csv\n",
    "Wrote ", nrow(suballele_freq), " sub-alleles to suballele_freq.csv"
  )
}

main()
