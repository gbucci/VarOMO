# Integration Guide: VCF Preprocessing with VEP/CIViC

## Overview

This guide describes the new **preprocessing module** added to the `koiosVrsCancerPipe` package. The module ensures that input VCFs are properly validated, cleaned, and re-annotated using authoritative sources (Ensembl VEP and CIViC) before entering the KOIOS-VRS pipeline.

## Why Preprocessing?

The original pipeline assumed VCFs were pre-annotated with ClinVar/gnomAD data. This created issues when:
- VCFs had inconsistent or missing annotation fields
- Annotation sources were outdated or non-standard
- Clinical significance classifications were ambiguous

The new preprocessing step **standardizes** all inputs by:
1. **Validating** VCF format and genome assembly (hg38/GRCh38)
2. **Stripping** all previous annotations to avoid conflicts
3. **Re-annotating** with fresh data from VEP (gnomAD frequencies) and CIViC (clinical significance)
4. **Ensuring compatibility** with downstream KOIOS filtering requirements

## Architecture Changes

### New File Structure

```
koios_vrs_pipeline/
├── R/
│   ├── koios_pipeline.R           # Original KOIOS-VRS pipeline (unchanged)
│   └── vcf_preprocessor.R         # NEW: VEP/CIViC annotation module
└── inst/scripts/
    ├── run_example.R              # Original single-phase runner
    └── run_complete_pipeline.R    # NEW: Integrated two-phase runner
```

### Pipeline Flow

```
Raw VCF (any state)
    ↓
[PHASE 1: vcf_preprocessor.R]
    ├─ Validate VCF format
    ├─ Strip old annotations
    ├─ VEP REST API (gnomAD_AF)
    └─ CIViC API (CLNSIG)
    ↓
Cleaned & Annotated VCF
    ↓
[PHASE 2: koios_pipeline.R]
    ├─ KOIOS enrichment (HGVSg, OMOP)
    ├─ VRS ID retrieval
    ├─ Frequency filtering (AF < 1%)
    └─ Phenopackets v2 generation
    ↓
Final Outputs (JSON, CSV)
```

## API Integration Details

### Ensembl VEP REST API

**Purpose:** Obtain population allele frequencies (gnomAD) and consequence predictions

**Endpoint:** `https://rest.ensembl.org/vep/human/region/`

**Key Parameters:**
- `canonical=1` - Returns only canonical transcripts
- `vcf_string=1` - Includes VCF-formatted response
- `variant_class=1` - Adds variant classification

**Data Extracted:**
- `gnomAD_AF` - Maximum allele frequency from gnomADg or gnomADe
- `most_severe_consequence` - VEP consequence term (e.g., missense_variant)

**Rate Limiting:** 0.3 seconds between requests (conservative)

### CIViC API (Clinical Interpretation of Variants in Cancer)

**Purpose:** Determine clinical significance for cancer-relevant variants

**Endpoint:** `https://civicdb.org/api/graphql`

**Query Method:** GraphQL query by genomic coordinates

**Significance Mapping:**
```
CIViC Term              → VCF CLNSIG
─────────────────────────────────────
Pathogenic/Oncogenic    → Pathogenic
Likely Pathogenic       → Likely_Pathogenic
Benign/Neutral          → Benign
Likely Benign           → Likely_Benign
Unknown/No match        → VUS
```

**Rate Limiting:** 0.5 seconds between requests (CIViC more restrictive)

## Compatibility with KOIOS Filters

The preprocessing module is **specifically designed** to produce VCF fields that match KOIOS pipeline expectations:

### Field Alignment

| VCF INFO Field | Source | KOIOS Pipeline Usage |
|----------------|--------|----------------------|
| `gnomAD_AF` | VEP REST API | → `koios_gnomAD_AF` → AF filtering (< 1%) |
| `CLNSIG` | CIViC API | → `concept_id` mapping → VUS vs. Actionable |
| `VEP_Consequence` | VEP REST API | Informational (not filtered) |

### Critical Filter: Allele Frequency

**Original Code (koios_pipeline.R:70-76):**
```r
freq_col <- "koios_gnomAD_AF"
if (freq_col %in% names(vcf.df)) {
    vcf.df[[freq_col]] <- suppressWarnings(as.numeric(vcf.df[[freq_col]]))
    exclude_freq <- !is.na(vcf.df[[freq_col]]) & vcf.df[[freq_col]] > af_threshold
    variants_for_vrs <- vcf.df[!exclude_freq, ]
}
```

**Preprocessing Ensures:**
- All variants have a `gnomAD_AF` INFO field (or `.` if unavailable)
- KOIOS maps this to `koios_gnomAD_AF` during `processClinGen()`
- Variants with AF > 1% are excluded **before VRS queries** (saves API calls)

### Clinical Significance Mapping

**Original Code (koios_pipeline.R:66-67):**
```r
vcf.df$concept_id[is.na(vcf.df$concept_id)] <- default_vus_id
vcf.df$concept_name[is.na(vcf.df$concept_name)] <- default_vus_name
```

**Preprocessing Ensures:**
- `CLNSIG` is populated for all variants (at minimum "VUS")
- KOIOS maps this to OMOP `concept_id`:
  - Pathogenic → specific OMOP ID
  - VUS → `1028197L` (default_vus_id)
- Output CSVs correctly split VUS from Actionable variants

## Usage

### Option 1: Integrated Pipeline (Recommended)

```r
# Run complete two-phase pipeline
source("inst/scripts/run_complete_pipeline.R")
```

This automatically:
1. Validates and re-annotates the raw VCF
2. Runs KOIOS-VRS pipeline on cleaned VCF
3. Generates all outputs

### Option 2: Manual Two-Phase Execution

```r
# Phase 1: Preprocessing only
source("R/vcf_preprocessor.R")

preprocess_and_annotate_vcf(
  vcf_path = "path/to/raw_sample.vcf",
  output_vcf_path = "path/to/cleaned_sample.vcf",
  genome_assembly = "GRCh38"
)

# Phase 2: Original KOIOS-VRS pipeline
source("R/koios_pipeline.R")

koios_vrs_pipeline(
  vcf_path = "path/to/cleaned_sample.vcf",
  output_base_name = "sample_output",
  af_threshold = 0.01
)
```

### Option 3: Preprocessing Only (Testing)

```r
source("R/vcf_preprocessor.R")

preprocess_and_annotate_vcf(
  vcf_path = "inst/extdata/melanoma_sample.vcf",
  output_vcf_path = "melanoma_cleaned.vcf"
)

# Inspect cleaned VCF before proceeding
```

## Testing the Integration

### Test with Provided Example

```bash
# From project root
Rscript inst/scripts/run_complete_pipeline.R
```

**Expected Output:**
```
━━━ PHASE 1: VCF Validation and Re-annotation ━━━
Step 1: Validating VCF format...
  ✓ VCF format validated
  ✓ Found 10 variant records

Step 2: Removing previous annotations...
  ✓ Stripped 10 records to baseline format

Step 3: Preparing variants for VEP annotation...
  ✓ Parsed 10 variants

Step 4: Annotating with Ensembl VEP REST API...
  Progress: 10 / 10
  ✓ VEP annotation complete

Step 5: Annotating with CIViC API...
  Progress: 10 / 10
  ✓ CIViC annotation complete

Step 6: Building annotated VCF...
  ✓ Annotated VCF written to: inst/extdata/melanoma_sample_cleaned.vcf

━━━ PHASE 2: KOIOS Enrichment and VRS Integration ━━━
[KOIOS processing output...]

✓ Pipeline completed successfully
```

### Validation Checks

1. **Cleaned VCF has standardized headers:**
   ```bash
   grep "^##INFO=<ID=gnomAD_AF" melanoma_sample_cleaned.vcf
   grep "^##INFO=<ID=CLNSIG" melanoma_sample_cleaned.vcf
   ```

2. **All variants have annotations:**
   ```bash
   grep -v "^#" melanoma_sample_cleaned.vcf | cut -f8 | grep "gnomAD_AF"
   ```

3. **Output CSVs match expectations:**
   - `_VUS.csv` contains only variants with `concept_id == 1028197L`
   - `_Note.csv` contains Pathogenic/Likely_Pathogenic variants
   - Both have `vrs_id` columns populated (where available)

## API Requirements & Rate Limits

### Ensembl VEP
- **No authentication required** (public REST API)
- **Rate limit:** ~15 requests/second (we use 0.3s intervals = safe)
- **Timeout:** 120 seconds per request (configurable)
- **Documentation:** https://rest.ensembl.org/documentation/info/vep_region_get