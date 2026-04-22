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
    c("--output_course"),
    type    = "character",
    default = NULL,
    help    = "Output file path for InstaPrism_props_df for the course label (TSV) [required]",
    metavar = "FILE"
  ),
  make_option(
    c("--output_fine"),
    type    = "character",
    default = NULL,
    help    = "Output file path for InstaPrism_props_df for the fine label (TSV) [required]",
    metavar = "FILE"
  ),
  make_option(
    c("--coarse_cell_label"),
    type    = "character",
    default = NULL,
    help    = "coarse cell label from the single cell data [required]"
  ),
  make_option(
    c("--fine_cell_label"),
    type    = "character",
    default = NULL,
    help    = "fine cell label from the single cell data [required]"
  ),
  make_option(
    c("--top_diff_gene_amount"),
    type    = "integer",
    default = 200L,
    help    = "the number of top genes per stage to use when doing top diff [default: %default]",
    metavar = "INT"
  ),
  make_option(
    c("--use_top_diff"),
    action = "store_true",
    default = FALSE,
    help    = "whether to use the top diff genes instead of all the genes [default: %default]",
  )
)

opt_parser <- OptionParser(
  usage       = "%prog [options]",
  option_list = option_list,
  description = "InstaPrism cell-type deconvolution for Pf bulk RNA-seq data."
)
opt <- parse_args(opt_parser)

# Validate required args
missing_args <- c()
if (is.null(opt$single_cell_r_object))  missing_args <- c(missing_args, "--single_cell_r_object")
if (is.null(opt$all_quants)) missing_args <- c(missing_args, "--all_quants")
if (is.null(opt$output_course))     missing_args <- c(missing_args, "--output_course")
if (is.null(opt$output_fine))     missing_args <- c(missing_args, "--output_fine")
if (is.null(opt$coarse_cell_label))     missing_args <- c(missing_args, "--coarse_cell_label")
if (is.null(opt$fine_cell_label))     missing_args <- c(missing_args, "--fine_cell_label")

if (length(missing_args) > 0) {
  print_help(opt_parser)
  stop(paste("Missing required argument(s):", paste(missing_args, collapse = ", ")), call. = FALSE)
}

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(InstaPrism)
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

if(opt$use_top_diff){
  message("get top diff genes")
  # Set identity to fine stage labels
  Idents(single_cell_r_object) <- opt$fine_cell_label

  # FindAllMarkers finds top markers for EVERY stage vs all others
  # this is the key function for finding stage-specific genes
  all_stage_markers <- FindAllMarkers(
    single_cell_r_object,
    test.use        = "wilcox",
    min.pct         = 0.25,
    logfc.threshold = 0.5,
    assay           = "SCT",   # counts need to be SCT set
    recorrect_umi   = FALSE,
    only.pos        = TRUE     # only upregulated markers, these are most useful for deconvolution
  )

  # Clean up and rank
  all_stage_markers_tidy <- all_stage_markers %>%
    as_tibble() %>%
    filter(p_val_adj < 0.05) %>%
    arrange(cluster, desc(avg_log2FC))

  # How many markers per stage?
  all_stage_markers_tidy_count = all_stage_markers_tidy %>% dplyr::count(cluster)

  # Top markers per stage
  all_stage_markers_tidy_top_per_stage = all_stage_markers_tidy %>%
    group_by(cluster) %>%
    slice_max(avg_log2FC, n = opt$top_diff_gene_amount)

  sc_counts = sc_counts[sort(unique(all_stage_markers_tidy_top_per_stage$gene)),]
}


rownames(sc_counts) <- gsub("PF3D7-", "PF3D7_", rownames(sc_counts))
single_cell_gene_names <- rownames(sc_counts)

cell_type_labels <- single_cell_r_object@meta.data |> pull(opt$coarse_cell_label)
cell_state_labels <- single_cell_r_object@meta.data |> pull(opt$fine_cell_label)

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

# refPrepare expects a gene x cells matrix with cell type labels
message("Preparing for InstaPrism deconvolution...")

ip_ref <- refPrepare(
  sc_Expr           = as.matrix(sc_counts[common_genes, ]),  # genes x cells
  cell.type.labels  = cell_type_labels,    # coarse
  cell.state.labels = cell_state_labels    # fine
)

# prepare bulk - needs to be genes x samples
bulk_for_ip <- as.matrix(all_quants_sp[, 2:ncol(all_quants_sp)])
rownames(bulk_for_ip) = all_quants_sp$new_name

# run InstaPrism
message("Running InstaPrism deconvolution...")
ip_result <- InstaPrism(
  bulk_Expr  = bulk_for_ip,
  refPhi_cs  = ip_ref
)

# extract cell type fractions
# Cell type level (coarse)
ip_props_type <- t(ip_result@Post.ini.ct@theta)
ip_props_type_df <- as.data.frame(ip_props_type) %>%
  tibble::rownames_to_column("sample") %>%
  tidyr::pivot_longer(-sample, names_to = "stage", values_to = "proportion")

# Cell state level (fine)
ip_props_state <- t(ip_result@Post.ini.cs@theta)
ip_props_state_df = as.data.frame(ip_props_state) %>%
      tibble::rownames_to_column("sample") %>%
      tidyr::pivot_longer(-sample, names_to = "stage", values_to = "proportion")

message(paste("Writing output to:", opt$output_fine, "and", opt$output_course))
readr::write_tsv(ip_props_state_df, opt$output_fine)
readr::write_tsv(ip_props_type_df, opt$output_course)
message("Done.")

