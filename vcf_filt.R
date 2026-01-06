#!/usr/bin/env Rscript

# Load required libraries
suppressPackageStartupMessages({
  library(optparse)
  library(tidyverse)
  library(vcfR)
})

# Define command-line options
option_list <- list(
  make_option(c("-v", "--vcf_file"), 
              type = "character", 
              default = NULL,
              help = "Input VCF file path [required]", 
              metavar = "FILE"),
  make_option(c("-q", "--qual_filter"), 
              type = "numeric", 
              default = 99,
              help = "Minimum quality score threshold [default = %default]", 
              metavar = "NUMBER"),
  make_option(c("-m", "--median_filt"), 
              type = "numeric", 
              default = 0.5,
              help = "Fraction of median depth for filtering [default = %default]", 
              metavar = "NUMBER"),
  make_option(c("-e", "--exclude_chrom"), 
              type = "character", 
              default = "S288C_chrmt",
              help = "Chromosome to exclude [default = %default]", 
              metavar = "STRING"),
  make_option(c("-o", "--out_file_vcf"), 
              type = "character", 
              default = NULL,
              help = "Output VCF file name [default = auto-generated from input]", 
              metavar = "FILE"),
  make_option(c("-c", "--out_file_csv"), 
              type = "character", 
              default = NULL,
              help = "Output CSV file name with CHROM and POS [default = auto-generated from input]", 
              metavar = "FILE"),
  make_option(c("--csv_only"), 
              action = "store_true",
              default = FALSE,
              help = "Only output CSV file (skip VCF output) [default = %default]"),
  make_option(c("--vcf_only"), 
              action = "store_true",
              default = FALSE,
              help = "Only output VCF file (skip CSV output) [default = %default]")
)

# Parse arguments
opt_parser <- OptionParser(
  option_list = option_list,
  description = "VCF SNP filtering of SK1 seq data"
)
opt <- parse_args(opt_parser)

# Check required arguments
if (is.null(opt$vcf_file)) {
  print_help(opt_parser)
  stop("Error: --vcf_file is required", call. = FALSE)
}

# Check if input file exists
if (!file.exists(opt$vcf_file)) {
  stop(paste("Error: Input file", opt$vcf_file, "does not exist"), call. = FALSE)
}

# Check for conflicting options
if (opt$csv_only && opt$vcf_only) {
  stop("Error: Cannot specify both --csv_only and --vcf_only", call. = FALSE)
}

# Generate output filenames if not provided
if (is.null(opt$out_file_vcf)) {
  opt$out_file_vcf <- paste0(str_extract(basename(opt$vcf_file), "^[\\w,-,_]+"), "_filtered.vcf.gz")
}

if (is.null(opt$out_file_csv)) {
  opt$out_file_csv <- paste0(str_extract(basename(opt$vcf_file), "^[\\w,-,_]+"), "_filtered.csv")
}

# Print parameters
cat("VCF SNP Filtering\n")
cat("=================\n")
cat("Input VCF file:", opt$vcf_file, "\n")
cat("Quality filter:", opt$qual_filter, "\n")
cat("Median depth fraction:", opt$median_filt, "\n")
cat("Exclude chromosome:", opt$exclude_chrom, "\n")
if (!opt$csv_only) {
  cat("Output VCF file:", opt$out_file_vcf, "\n")
}
if (!opt$vcf_only) {
  cat("Output CSV file:", opt$out_file_csv, "\n")
}
cat("\n")

# Load SNPs from seq data
cat("Loading VCF file...\n")
SK1 <- read.vcfR(opt$vcf_file)
cat("Loaded", nrow(SK1@fix), "variants\n\n")

# Get median depth
cat("Calculating median depth...\n")
depths <- extract.info(SK1, element = "DP", as.numeric = TRUE)
sk1_median <- median(depths, na.rm = TRUE)
cat("Median depth:", sk1_median, "\n")
cat("Depth threshold:", sk1_median * opt$median_filt, "\n\n")

# SNP filters
cat("Applying filters...\n")

# Filter for only biallelics (single ALT)
SK1_filt <- SK1[is.biallelic(SK1)]
cat("After biallelic filter:", nrow(SK1_filt@fix), "variants\n")

# Filter for SNPs only (remove indels)
SK1_filt <- SK1_filt[!is.indel(SK1_filt)]
cat("After indel filter:", nrow(SK1_filt@fix), "variants\n")

# Filter by quality
SK1_filt <- SK1_filt[as.numeric(SK1_filt@fix[, "QUAL"]) > opt$qual_filter]
cat("After quality filter (>", opt$qual_filter, "):", nrow(SK1_filt@fix), "variants\n")

# Filter by depth
depths_filt <- extract.info(SK1_filt, element = "DP", as.numeric = TRUE)
SK1_filt <- SK1_filt[depths_filt >= sk1_median * opt$median_filt]
cat("After depth filter:", nrow(SK1_filt@fix), "variants\n")

# Filter for AC == 2
ac2 <- extract.info(SK1_filt, element = "AC", as.numeric = TRUE)
SK1_filt <- SK1_filt[ac2 == 2]
cat("After allele count filter (AC == 2):", nrow(SK1_filt@fix), "variants\n")

# Remove excluded chromosome
SK1_filt <- SK1_filt[SK1_filt@fix[, "CHROM"] != opt$exclude_chrom]
cat("After excluding", opt$exclude_chrom, ":", nrow(SK1_filt@fix), "variants\n\n")

# Write output files
if (!opt$csv_only) {
  cat("Writing filtered VCF file to:", opt$out_file_vcf, "\n")
  write.vcf(SK1_filt, opt$out_file_vcf)
}

if (!opt$vcf_only) {
  cat("Writing CSV file with CHROM and POS to:", opt$out_file_csv, "\n")
  # Extract CHROM and POS from the filtered VCF
  snp_positions <- tibble(
    `#CHROM` = SK1_filt@fix[, "CHROM"],
    POS = as.numeric(SK1_filt@fix[, "POS"])
  )
  write_csv(snp_positions, opt$out_file_csv)
}

cat("\nFiltering complete!\n")