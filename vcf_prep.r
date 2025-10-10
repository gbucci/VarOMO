# vcf_preprocessor.R

#' VCF Preprocessing and Re-annotation Pipeline
#'
#' Validates VCF format, removes previous annotations, and re-annotates using
#' Ensembl VEP REST API and CIViC API. Ensures compatibility with downstream
#' KOIOS pipeline requirements (gnomAD_AF, CLNSIG fields).
#'
#' @param vcf_path Path to input VCF file (must be hg38/GRCh38)
#' @param output_vcf_path Path for cleaned and re-annotated VCF output
#' @param genome_assembly Genome assembly version (default: "GRCh38")
#' @param vep_timeout VEP API timeout in seconds (default: 120)
#' @param civic_timeout CIViC API timeout in seconds (default: 60)
#' @return Invisible TRUE on success, stops on error
#' @export
preprocess_and_annotate_vcf <- function(vcf_path, 
                                        output_vcf_path,
                                        genome_assembly = "GRCh38",
                                        vep_timeout = 120,
                                        civic_timeout = 60) {
  
  # --- 1. VALIDATE INPUT VCF ---
  
  cat("Step 1: Validating VCF format...\n")
  
  if (!file.exists(vcf_path)) {
    stop("VCF file not found: ", vcf_path)
  }
  
  # Read VCF
  vcf_lines <- readLines(vcf_path)
  
  # Check basic VCF structure
  if (!any(grepl("^##fileformat=VCF", vcf_lines))) {
    stop("Invalid VCF: Missing ##fileformat header")
  }
  
  header_idx <- which(grepl("^#CHROM", vcf_lines))
  if (length(header_idx) != 1) {
    stop("Invalid VCF: Missing or multiple #CHROM header lines")
  }
  
  # Validate genome assembly
  ref_lines <- vcf_lines[grepl("^##reference=", vcf_lines)]
  if (length(ref_lines) > 0) {
    if (!any(grepl("hg38|GRCh38", ref_lines, ignore.case = TRUE))) {
      warning("VCF reference is not explicitly hg38/GRCh38. Proceeding with caution.")
    }
  } else {
    warning("No ##reference header found. Assuming hg38/GRCh38.")
  }
  
  # Extract header and data
  header_line <- vcf_lines[header_idx]
  data_start <- header_idx + 1
  
  if (data_start > length(vcf_lines)) {
    stop("VCF contains no variant records")
  }
  
  cat("  ✓ VCF format validated\n")
  cat("  ✓ Found", length(vcf_lines) - data_start + 1, "variant records\n")
  
  # --- 2. CLEAN VCF (Remove Previous Annotations) ---
  
  cat("\nStep 2: Removing previous annotations...\n")
  
  # Keep only essential headers
  essential_headers <- c(
    "##fileformat=VCFv4.2",
    paste0("##fileDate=", format(Sys.Date(), "%Y%m%d")),
    paste0("##reference=", genome_assembly),
    '##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">',
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
    '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic Depths">',
    '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">'
  )
  
  # Parse data lines and strip INFO field
  data_lines <- vcf_lines[data_start:length(vcf_lines)]
  cleaned_data <- lapply(data_lines, function(line) {
    fields <- strsplit(line, "\t")[[1]]
    if (length(fields) < 8) return(NULL)
    
    # Keep CHROM, POS, ID, REF, ALT, QUAL, FILTER; reset INFO to "."
    fields[8] <- "."
    paste(fields, collapse = "\t")
  })
  
  cleaned_data <- Filter(Negate(is.null), cleaned_data)
  
  cat("  ✓ Stripped", length(vcf_lines) - header_idx - 1, "records to baseline format\n")
  
  # --- 3. PARSE VARIANTS FOR ANNOTATION ---
  
  cat("\nStep 3: Preparing variants for VEP annotation...\n")
  
  variants <- lapply(cleaned_data, function(line) {
    fields <- strsplit(line, "\t")[[1]]
    list(
      chrom = sub("^chr", "", fields[1]),  # VEP expects without 'chr'
      pos = fields[2],
      id = fields[3],
      ref = fields[4],
      alt = fields[5],
      original_line = line
    )
  })
  
  cat("  ✓ Parsed", length(variants), "variants\n")
  
  # --- 4. VEP ANNOTATION (REST API) ---
  
  cat("\nStep 4: Annotating with Ensembl VEP REST API...\n")
  cat("  (This may take several minutes for large VCFs)\n")
  
  vep_base_url <- "https://rest.ensembl.org"
  
  annotate_variant_vep <- function(var, assembly = "GRCh38", timeout = 120) {
    # Construct VEP input format: "chrom start end allele/allele strand"
    pos_num <- as.numeric(var$pos)
    allele_string <- paste(var$ref, var$alt, sep = "/")
    
    # VEP REST endpoint
    endpoint <- paste0("/vep/human/region/", var$chrom, ":", pos_num, "-", pos_num, ":", allele_string, "/1")
    url <- paste0(vep_base_url, endpoint)
    
    response <- tryCatch({
      httr::GET(
        url,
        httr::timeout(timeout),
        httr::add_headers(
          "Content-Type" = "application/json"
        ),
        query = list(
          canonical = 1,
          vcf_string = 1,
          variant_class = 1
        )
      )
    }, error = function(e) {
      warning("VEP API error for ", var$chrom, ":", var$pos, " - ", e$message)
      return(NULL)
    })
    
    if (is.null(response) || httr::http_error(response)) {
      return(list(gnomad_af = NA, consequence = NA))
    }
    
    content <- tryCatch({
      jsonlite::fromJSON(httr::content(response, as = "text", encoding = "UTF-8"))
    }, error = function(e) {
      return(NULL)
    })
    
    if (is.null(content) || length(content) == 0) {
      return(list(gnomad_af = NA, consequence = NA))
    }
    
    # Extract gnomAD AF (using gnomADg or gnomADe)
    gnomad_af <- NA
    if (!is.null(content$colocated_variants)) {
      for (cv in content$colocated_variants) {
        if (!is.null(cv$frequencies)) {
          if (!is.null(cv$frequencies$gnomADg)) {
            af_values <- unlist(cv$frequencies$gnomADg)
            gnomad_af <- max(af_values, na.rm = TRUE)
            break
          } else if (!is.null(cv$frequencies$gnomADe)) {
            af_values <- unlist(cv$frequencies$gnomADe)
            gnomad_af <- max(af_values, na.rm = TRUE)
            break
          }
        }
      }
    }
    
    # Extract most severe consequence
    consequence <- NA
    if (!is.null(content$most_severe_consequence)) {
      consequence <- content$most_severe_consequence
    }
    
    list(gnomad_af = gnomad_af, consequence = consequence)
  }
  
  vep_annotations <- vector("list", length(variants))
  for (i in seq_along(variants)) {
    if (i %% 10 == 0) {
      cat("  Progress:", i, "/", length(variants), "\n")
    }
    vep_annotations[[i]] <- annotate_variant_vep(variants[[i]], assembly = genome_assembly, timeout = vep_timeout)
    Sys.sleep(0.3)  # Rate limiting
  }
  
  cat("  ✓ VEP annotation complete\n")
  
  # --- 5. CIViC ANNOTATION (Clinical Significance) ---
  
  cat("\nStep 5: Annotating with CIViC API for clinical significance...\n")
  
  civic_base_url <- "https://civicdb.org/api"
  
  annotate_variant_civic <- function(var, timeout = 60) {
    # CIViC search by coordinate
    chrom_clean <- sub("^chr", "", var$chrom)
    
    endpoint <- "/graphql"
    url <- paste0(civic_base_url, endpoint)
    
    # GraphQL query for variant at position
    query <- sprintf('
    {
      variants(
        entrezSymbol: "*"
        startPosition: %s
        endPosition: %s
        referenceBuild: GRCH38
      ) {
        edges {
          node {
            clinicalSignificance
            variantTypes {
              name
            }
          }
        }
      }
    }', var$pos, var$pos)
    
    response <- tryCatch({
      httr::POST(
        url,
        httr::timeout(timeout),
        httr::add_headers("Content-Type" = "application/json"),
        body = jsonlite::toJSON(list(query = query), auto_unbox = TRUE)
      )
    }, error = function(e) {
      return(NULL)
    })
    
    if (is.null(response) || httr::http_error(response)) {
      return(NA)
    }
    
    content <- tryCatch({
      jsonlite::fromJSON(httr::content(response, as = "text", encoding = "UTF-8"))
    }, error = function(e) {
      return(NULL)
    })
    
    # Extract clinical significance
    if (!is.null(content$data$variants$edges) && length(content$data$variants$edges) > 0) {
      sig <- content$data$variants$edges[[1]]$node$clinicalSignificance
      if (!is.null(sig) && length(sig) > 0) {
        return(sig[1])
      }
    }
    
    return(NA)
  }
  
  civic_annotations <- vector("character", length(variants))
  for (i in seq_along(variants)) {
    if (i %% 10 == 0) {
      cat("  Progress:", i, "/", length(variants), "\n")
    }
    civic_annotations[i] <- annotate_variant_civic(variants[[i]], timeout = civic_timeout)
    Sys.sleep(0.5)  # Rate limiting (CIViC is more restrictive)
  }
  
  # Map CIViC significance to ClinVar-like terms
  map_civic_to_clnsig <- function(civic_sig) {
    if (is.na(civic_sig)) return("VUS")
    
    civic_lower <- tolower(civic_sig)
    if (grepl("pathogenic|oncogenic|sensitiv", civic_lower)) return("Pathogenic")
    if (grepl("likely.*pathogenic", civic_lower)) return("Likely_Pathogenic")
    if (grepl("benign|neutral", civic_lower)) return("Benign")
    if (grepl("likely.*benign", civic_lower)) return("Likely_Benign")
    
    return("VUS")
  }
  
  clnsig_values <- sapply(civic_annotations, map_civic_to_clnsig)
  
  cat("  ✓ CIViC annotation complete\n")
  
  # --- 6. CONSTRUCT ANNOTATED VCF ---
  
  cat("\nStep 6: Building annotated VCF...\n")
  
  # Add annotation headers
  annotation_headers <- c(
    '##INFO=<ID=gnomAD_AF,Number=A,Type=Float,Description="gnomAD Allele Frequency from VEP">',
    '##INFO=<ID=CLNSIG,Number=.,Type=String,Description="Clinical Significance from CIViC">',
    '##INFO=<ID=VEP_Consequence,Number=.,Type=String,Description="Most severe consequence from VEP">'
  )
  
  # Add contig headers (extract unique chroms)
  unique_chroms <- unique(sapply(variants, function(v) paste0("chr", v$chrom)))
  contig_headers <- sprintf('##contig=<ID=%s>', unique_chroms)
  
  all_headers <- c(essential_headers, annotation_headers, contig_headers, header_line)
  
  # Build annotated data lines
  annotated_lines <- vector("character", length(variants))
  for (i in seq_along(variants)) {
    fields <- strsplit(variants[[i]]$original_line, "\t")[[1]]
    
    # Build INFO field
    info_parts <- c()
    
    # Add DP if present in original FORMAT fields
    if (length(fields) >= 10) {
      format_fields <- strsplit(fields[9], ":")[[1]]
      sample_values <- strsplit(fields[10], ":")[[1]]
      dp_idx <- which(format_fields == "DP")
      if (length(dp_idx) > 0 && length(sample_values) >= dp_idx) {
        info_parts <- c(info_parts, paste0("DP=", sample_values[dp_idx]))
      }
    }
    
    # Add gnomAD_AF
    if (!is.na(vep_annotations[[i]]$gnomad_af) && !is.infinite(vep_annotations[[i]]$gnomad_af)) {
      info_parts <- c(info_parts, sprintf("gnomAD_AF=%.6f", vep_annotations[[i]]$gnomad_af))
    } else {
      info_parts <- c(info_parts, "gnomAD_AF=.")
    }
    
    # Add CLNSIG
    info_parts <- c(info_parts, paste0("CLNSIG=", clnsig_values[i]))
    
    # Add VEP consequence
    if (!is.na(vep_annotations[[i]]$consequence)) {
      info_parts <- c(info_parts, paste0("VEP_Consequence=", vep_annotations[[i]]$consequence))
    }
    
    fields[8] <- paste(info_parts, collapse = ";")
    
    # Restore chr prefix
    fields[1] <- paste0("chr", variants[[i]]$chrom)
    
    annotated_lines[i] <- paste(fields, collapse = "\t")
  }
  
  # --- 7. WRITE OUTPUT VCF ---
  
  output_content <- c(all_headers, annotated_lines)
  writeLines(output_content, output_vcf_path)
  
  cat("  ✓ Annotated VCF written to:", output_vcf_path, "\n")
  
  # --- 8. VALIDATION SUMMARY ---
  
  cat("\n=== Annotation Summary ===\n")
  cat("Total variants:", length(variants), "\n")
  cat("With gnomAD AF:", sum(!is.na(sapply(vep_annotations, function(x) x$gnomad_af))), "\n")
  cat("Clinical Significance breakdown:\n")
  print(table(clnsig_values))
  cat("\n✓ Preprocessing complete. Output ready for KOIOS pipeline.\n")
  
  invisible(TRUE)
}
