#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse)
})

# CLI options
option_list <- list(
  make_option(
    c("--single_cell_r_object"),
    type    = "character",
    default = NULL,
    help    = "Path to the single cell data R Seurat object with normalized data already ran RDS file [required]",
    metavar = "FILE"
  ),
  make_option(
    c("--all_quants"),
    type    = "character",
    default = NULL,
    help    = "Path to the all_quants TSV/TSV.gz file [required]",
    metavar = "FILE"
  ),
  make_option(
    c("--seed"),
    type    = "integer",
    default = 42L,
    help    = "Random seed for pseudo-donor assignment [default: %default]",
    metavar = "INT"
  ),
  make_option(
    c("--output"),
    type    = "character",
    default = NULL,
    help    = "Output file path for music_props_df (TSV) [required]",
    metavar = "FILE"
  )
)

opt_parser <- OptionParser(
  usage       = "%prog [options]",
  option_list = option_list,
  description = "MuSiC cell-type deconvolution for Pf bulk RNA-seq data."
)
opt <- parse_args(opt_parser)

# Validate required args
missing_args <- c()
if (is.null(opt$single_cell_r_object))  missing_args <- c(missing_args, "--single_cell_r_object")
if (is.null(opt$all_quants)) missing_args <- c(missing_args, "--all_quants")
if (is.null(opt$output))     missing_args <- c(missing_args, "--output")

if (length(missing_args) > 0) {
  print_help(opt_parser)
  stop(paste("Missing required argument(s):", paste(missing_args, collapse = ", ")), call. = FALSE)
}

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(MuSiC)
  library(Biobase)
  library(SingleCellExperiment)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(validate)
})

# Helper function
process_qunat_sf_for_3d7 <- function(quant_fnp, accetable_genes) {
  quant_dat <- readr::read_tsv(quant_fnp, show_col_types = FALSE)

  quant_dat_3d7 <- quant_dat %>%
    filter(grepl("^PF3D7", Name)) %>%
    mutate(in_sc = Name %in% accetable_genes) %>%
    mutate(new_name = ifelse(in_sc, Name, gsub("\\..*", "", Name))) %>%
    mutate(new_name_in_sc = new_name %in% accetable_genes)

  quant_dat_3d7_in_sc <- quant_dat_3d7 %>%
    filter(new_name_in_sc)

  return(quant_dat_3d7_in_sc)
}

# Validate all_quants schema
message("Reading and validating all_quants file...")

# Peek at the raw file first to check required columns exist before processing
raw_quants <- readr::read_tsv(opt$all_quants, show_col_types = FALSE)

quants_rules <- validator(
  has_Name    = is.character(Name),
  has_TPM     = is.double(TPM),
  has_sample  = is.character(sample),
  Name_not_na = !is.na(Name),
  TPM_not_na  = !is.na(TPM),
  TPM_non_neg = TPM >= 0
)

quants_check <- confront(raw_quants, quants_rules)
quants_summary <- summary(quants_check)

if (any(!quants_summary$passes & quants_summary$items > 0)) {
  message("Validation results for all_quants:")
  print(quants_summary)
  failing <- quants_summary[!quants_summary$passes & quants_summary$items > 0, ]
  stop(paste(
    "all_quants file failed validation for rule(s):",
    paste(failing$name, collapse = ", ")
  ), call. = FALSE)
}
message("all_quants validation passed.")

# Load SC reference
message("Loading nf54 SCT object...")
single_cell_r_object <- readRDS(opt$single_cell_r_object)

sc_counts <- GetAssayData(single_cell_r_object, assay = "RNA", layer = "counts")
rownames(sc_counts) <- gsub("PF3D7-", "PF3D7_", rownames(sc_counts))
single_cell_gene_names <- rownames(sc_counts)

cell_labels <- single_cell_r_object$STAGE_HR2

# Process bulk quants
message("Processing bulk quants...")
all_quants <- process_qunat_sf_for_3d7(opt$all_quants, single_cell_gene_names)

all_quants_sp <- all_quants %>%
  ungroup() %>%
  select(new_name, sample, TPM) %>%
  spread(sample, TPM, fill = 0)

common_genes <- all_quants_sp$new_name

# Build SingleCellExperiment with pseudo-donors
set.seed(opt$seed)
n_pseudodonors <- 5
pseudo_donors <- paste0("donor_", sample(
  seq_len(n_pseudodonors),
  ncol(sc_counts),
  replace = TRUE
))

sc_sce <- SingleCellExperiment(
  assays  = list(counts = as.matrix(sc_counts[common_genes, ])),
  colData = DataFrame(
    cellType    = cell_labels,
    pseudoDonor = pseudo_donors
  )
)

# Build bulk ExpressionSet
bulk_mat <- as.matrix(all_quants_sp[, 2:ncol(all_quants_sp)])
rownames(bulk_mat) <- all_quants_sp$new_name
bulk_eset <- ExpressionSet(assayData = bulk_mat)

# Run MuSiC
message("Running MuSiC deconvolution...")
music_result <- music_prop(
  bulk.mtx = bulk_mat,
  sc.sce   = sc_sce,
  clusters = "cellType",
  samples  = "pseudoDonor",
  verbose  = TRUE
)

# Extract and write results
music_props <- music_result$Est.prop.weighted

music_props_df <- as.data.frame(music_props) %>%
  tibble::rownames_to_column("sample") %>%
  tidyr::pivot_longer(-sample, names_to = "stage", values_to = "proportion")

message(paste("Writing output to:", opt$output))
readr::write_tsv(music_props_df, opt$output)
message("Done.")


