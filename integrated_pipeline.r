# run_complete_pipeline.R

#' Complete Integrated Pipeline: VCF Validation → VEP/CIViC Annotation → KOIOS → VRS → Phenopackets
#'
#' This script integrates the preprocessing step with the existing KOIOS-VRS pipeline.
#' Usage: Rscript inst/scripts/run_complete_pipeline.R

# --- SETUP ---

# Load required libraries
if (!requireNamespace("httr", quietly = TRUE)) stop("Package 'httr' required")
if (!requireNamespace("jsonlite", quietly = TRUE)) stop("Package 'jsonlite' required")
if (!requireNamespace("vcfR", quietly = TRUE)) stop("Package 'vcfR' required")

# Source both pipeline scripts
source("R/vcf_preprocessor.R")
source("R/koios_pipeline.R")

# --- CONFIGURATION ---

# Input files (example: melanoma sample)
raw_vcf <- "inst/extdata/melanoma_sample.vcf"
cleaned_vcf <- "inst/extdata/melanoma_sample_cleaned.vcf"
output_base <- "melanoma_final"

# Pipeline parameters
AF_THRESHOLD <- 0.01  # 1% allele frequency filter
GENOME_ASSEMBLY <- "GRCh38"

# --- EXECUTION ---

cat("╔════════════════════════════════════════════════════════════╗\n")
cat("║  COMPLETE CANCER VARIANT ANNOTATION & REPORTING PIPELINE  ║\n")
cat("╚════════════════════════════════════════════════════════════╝\n\n")

# PHASE 1: VCF Preprocessing and Re-annotation
cat("━━━ PHASE 1: VCF Validation and Re-annotation ━━━\n\n")

tryCatch({
  preprocess_and_annotate_vcf(
    vcf_path = raw_vcf,
    output_vcf_path = cleaned_vcf,
    genome_assembly = GENOME_ASSEMBLY,
    vep_timeout = 120,
    civic_timeout = 60
  )
}, error = function(e) {
  cat("ERROR in preprocessing phase:\n")
  cat(e$message, "\n")
  stop("Pipeline aborted.")
})

cat("\n━━━ PHASE 2: KOIOS Enrichment and VRS Integration ━━━\n\n")

# PHASE 2: KOIOS-VRS Pipeline
tryCatch({
  results <- koios_vrs_pipeline(
    vcf_path = cleaned_vcf,
    output_base_name = output_base,
    af_threshold = AF_THRESHOLD
  )
  
  cat("\n╔════════════════════════════════════════════════════════════╗\n")
  cat("║                  PIPELINE COMPLETED                        ║\n")
  cat("╚════════════════════════════════════════════════════════════╝\n\n")
  
  cat("📄 Output Files Generated:\n")
  cat("  • Cleaned VCF:      ", cleaned_vcf, "\n")
  cat("  • Phenopackets:     ", paste0(output_base, "_Phenopackets.json"), "\n")
  cat("  • VUS Report:       ", paste0(output_base, "_VUS.csv"), "\n")
  cat("  • Actionable Report:", paste0(output_base, "_Note.csv"), "\n\n")
  
  cat("📊 Results Summary:\n")
  cat("  • VUS variants:     ", nrow(results$vus_df), "\n")
  cat("  • Actionable vars:  ", nrow(results$note_df), "\n")
  
}, error = function(e) {
  cat("ERROR in KOIOS-VRS phase:\n")
  cat(e$message, "\n")
  stop("Pipeline aborted.")
})

cat("\n✓ All phases completed successfully.\n")
