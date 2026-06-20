# ============================================================================
# RSTUDIO GENOME ASSEMBLY PIPELINE FOR SRR25083113
# ============================================================================
# Interactive pipeline designed for RStudio on macOS ARM (Apple Silicon).
# QUAST and Prokka are not natively supported on macOS ARM, so this pipeline
# uses SeqKit (quality assessment) and a modular Prodigal + Barrnap + ARAGORN
# workflow (structural annotation) in their place, matching the Weblem-6
# report methodology.
#
# Just source this file or run sections interactively!
# ============================================================================

cat("\n")
cat("================================================================\n")
cat("   GENOME ASSEMBLY PIPELINE - RSTUDIO VERSION                  \n")
cat("   Run this directly in RStudio!                               \n")
cat("================================================================\n\n")

# ============================================================================
# RSTUDIO SETUP INSTRUCTIONS
# ============================================================================
cat("This script will:\n")
cat("  1) Download SRR25083113 from SRA (fastq-dump)\n")
cat("  2) Assemble the genome with SPAdes\n")
cat("  3) Evaluate assembly with SeqKit\n")
cat("  4) Annotate with Prodigal (CDS) + Barrnap (rRNA) + ARAGORN (tRNA)\n")
cat("  5) Summarize results in a text report\n")
cat("     (For the publication figures, run 02_generate_figures.R\n")
cat("      after this script completes.)\n\n")

cat("Make sure you have these tools installed (e.g. via conda):\n")
cat("  spades.py, seqkit, prodigal, barrnap, aragorn, fastq-dump\n\n")

cat("If you use conda, for example:\n")
cat("  conda install -c bioconda spades seqkit prodigal barrnap aragorn sra-tools\n\n")

cat("After sourcing this file, run:\n")
cat("  run_pipeline_interactive()\n")
cat("to execute the full pipeline.\n\n")

# ============================================================================
# AUTOMATIC PACKAGE INSTALLATION
# ============================================================================

cat("===================================================================\n")
cat("STEP 1: Installing/Loading Required Packages\n")
cat("===================================================================\n\n")

install_if_needed <- function(packages, bioc = FALSE) {
  for (pkg in packages) {
    if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
      cat(sprintf("Installing %s...\n", pkg))
      if (bioc) {
        if (!require("BiocManager", quietly = TRUE)) {
          install.packages("BiocManager", repos = "https://cran.r-project.org")
        }
        BiocManager::install(pkg, update = FALSE, ask = FALSE, quietly = TRUE)
      } else {
        install.packages(pkg, repos = "https://cran.r-project.org", quiet = TRUE)
      }
      library(pkg, character.only = TRUE)
    }
  }
}

cat("Installing CRAN packages (if needed)...\n")
cran_pkgs <- c("tidyverse", "ggplot2", "cowplot")
install_if_needed(cran_pkgs, bioc = FALSE)

cat("\nInstalling Bioconductor packages (required for figure generation)...\n")
bioc_pkgs <- c("Biostrings")
tryCatch({
  install_if_needed(bioc_pkgs, bioc = TRUE)
  BIOC_AVAILABLE <- TRUE
}, error = function(e) {
  cat("Warning: Biostrings not available - figure generation script will fail\n")
  BIOC_AVAILABLE <<- FALSE
})

cat("\nPackages loaded.\n\n")

# ============================================================================
# CHECK EXTERNAL TOOLS (INTERACTIVE)
# ============================================================================

cat("===================================================================\n")
cat("STEP 2: Checking External Tools\n")
cat("===================================================================\n\n")

check_tool <- function(tool) {
  result <- suppressWarnings(system2("which", tool, stdout = TRUE, stderr = FALSE))
  exists <- length(result) > 0 && !grepl("not found", result[1], ignore.case = TRUE)

  if (exists) {
    cat(sprintf("[OK] %s found\n", tool))
    return(TRUE)
  } else {
    cat(sprintf("[MISSING] %s not found\n", tool))
    return(FALSE)
  }
}

tools_status <- list(
  spades   = check_tool("spades.py"),
  seqkit   = check_tool("seqkit"),
  prodigal = check_tool("prodigal"),
  barrnap  = check_tool("barrnap"),
  aragorn  = check_tool("aragorn"),
  sra      = check_tool("fastq-dump")
)

all_tools_available <- all(unlist(tools_status))

cat("\n")

if (!all_tools_available) {
  cat("SOME TOOLS ARE MISSING\n\n")
  cat("Install them via conda:\n")
  cat("  conda install -c bioconda spades seqkit prodigal barrnap aragorn sra-tools\n\n")
}

cat("-------------------------------------------------------------------\n\n")

# ============================================================================
# CONFIGURATION
# ============================================================================

cat("===================================================================\n")
cat("STEP 3: Setting Configuration\n")
cat("===================================================================\n\n")

CONFIG <- list(
  sra_accession = "SRR25083113",
  base_dir      = "SRR25083113_analysis",
  threads       = 4
)

cat(sprintf("Dataset:       %s\n", CONFIG$sra_accession))
cat(sprintf("Output:        %s/\n", CONFIG$base_dir))
cat(sprintf("Threads:       %d\n\n", CONFIG$threads))

dirs <- c(
  CONFIG$base_dir,
  file.path(CONFIG$base_dir, "raw_data"),
  file.path(CONFIG$base_dir, "assembly"),
  file.path(CONFIG$base_dir, "quality"),
  file.path(CONFIG$base_dir, "annotation"),
  file.path(CONFIG$base_dir, "reports"),
  file.path(CONFIG$base_dir, "plots")
)

for (d in dirs) {
  if (!dir.exists(d)) {
    dir.create(d, recursive = TRUE, showWarnings = FALSE)
  }
}

# ============================================================================
# PIPELINE FUNCTIONS
# ============================================================================

#' Download SRA data (SRR25083113)
download_data <- function() {
  cat("===================================================================\n")
  cat("STEP 4: Downloading Data\n")
  cat("===================================================================\n\n")

  raw_dir <- file.path(CONFIG$base_dir, "raw_data")
  r1 <- file.path(raw_dir, paste0(CONFIG$sra_accession, "_1.fastq.gz"))
  r2 <- file.path(raw_dir, paste0(CONFIG$sra_accession, "_2.fastq.gz"))

  if (file.exists(r1) && file.exists(r2)) {
    cat("Data already downloaded.\n\n")
    return(list(r1 = r1, r2 = r2))
  }

  cat("Downloading from NCBI SRA...\n\n")

  cmd <- sprintf("fastq-dump --split-files --gzip --outdir %s %s",
                 raw_dir, CONFIG$sra_accession)

  result <- system(cmd)

  if (result == 0 && file.exists(r1) && file.exists(r2)) {
    cat("\nDownload complete.\n\n")
    return(list(r1 = r1, r2 = r2))
  } else {
    stop("Download failed. Check internet connection and SRA Toolkit.")
  }
}

#' Run SPAdes assembly
#' NOTE: careful mode is intentionally disabled (default SPAdes mode) to
#' match the Weblem-6 report methodology: "Parameters: default (careful
#' mode disabled)".
run_assembly <- function(fastq_files) {
  cat("===================================================================\n")
  cat("STEP 5: Genome Assembly (SPAdes, default mode)\n")
  cat("===================================================================\n\n")

  assembly_dir <- file.path(CONFIG$base_dir, "assembly")
  contigs <- file.path(assembly_dir, "contigs.fasta")

  if (file.exists(contigs)) {
    cat("Assembly already exists.\n\n")
    return(contigs)
  }

  cat("Running SPAdes (default parameters)...\n\n")

  cmd <- sprintf("spades.py -1 %s -2 %s -o %s -t %d",
                 fastq_files$r1, fastq_files$r2, assembly_dir, CONFIG$threads)

  system(cmd)

  if (file.exists(contigs)) {
    cat("\nAssembly complete.\n\n")
    return(contigs)
  } else {
    stop("Assembly failed (contigs.fasta not found).")
  }
}

#' Run SeqKit quality assessment (QUAST unavailable on macOS ARM)
run_quality <- function(contigs) {
  cat("===================================================================\n")
  cat("STEP 6: Quality Assessment (SeqKit)\n")
  cat("===================================================================\n\n")

  quality_dir <- file.path(CONFIG$base_dir, "quality")
  report_file <- file.path(quality_dir, "report.tsv")

  cmd <- sprintf("seqkit stats -a -T %s > %s", contigs, report_file)
  system(cmd)

  report <- read.table(report_file,
                        header = TRUE,
                        sep = "",        # split on any whitespace
                        check.names = TRUE)

  # check.names = TRUE mangles "GC(%)" into "GC..."
  metrics <- list(
    contigs = report$num_seqs[1],
    length  = report$sum_len[1],
    largest = report$max_len[1],
    n50     = report$N50[1],
    gc      = report$GC...[1]
  )

  cat("\n=== QUALITY METRICS ===\n")
  cat(sprintf("Contigs:        %d\n", metrics$contigs))
  cat(sprintf("Total length:   %.2f Mb\n", metrics$length / 1e6))
  cat(sprintf("Largest contig: %d bp\n", metrics$largest))
  cat(sprintf("N50:            %.1f kb\n", metrics$n50 / 1000))
  cat(sprintf("GC content:     %.2f%%\n\n", metrics$gc))

  cat("Quality assessment complete.\n\n")

  return(metrics)
}

#' Run structural annotation: Prodigal (CDS) + Barrnap (rRNA) + ARAGORN (tRNA)
#' (Prokka is incompatible with macOS M1 / ARM)
run_annotation <- function(contigs) {
  cat("===================================================================\n")
  cat("STEP 7: Genome Annotation (Prodigal + Barrnap + ARAGORN)\n")
  cat("===================================================================\n\n")

  anno_dir <- file.path(CONFIG$base_dir, "annotation")
  proteins_faa <- file.path(anno_dir, "proteins.faa")
  genes_ffn    <- file.path(anno_dir, "genes.ffn")
  rrna_gff     <- file.path(anno_dir, "rRNAs.gff")
  trna_txt     <- file.path(anno_dir, "tRNAs.txt")

  cat("Running Prodigal (CDS prediction)...\n")
  system(sprintf("prodigal -i %s -a %s -d %s", contigs, proteins_faa, genes_ffn))

  cat("Running Barrnap (rRNA detection)...\n")
  system(sprintf("barrnap --kingdom bac %s > %s", contigs, rrna_gff))

  cat("Running ARAGORN (tRNA detection)...\n\n")
  system(sprintf("aragorn -t %s -o %s", contigs, trna_txt))

  if (file.exists(proteins_faa) && file.exists(rrna_gff) && file.exists(trna_txt)) {
    cds_count  <- length(Biostrings::readAAStringSet(proteins_faa))

    rrna_tbl   <- read.table(rrna_gff, sep = "\t", header = FALSE,
                              comment.char = "#", stringsAsFactors = FALSE)
    rrna_count <- nrow(rrna_tbl)

    trna_lines <- readLines(trna_txt)
    trna_count <- sum(grepl("tRNA", trna_lines, ignore.case = TRUE))

    stats <- list(cds = cds_count, rRNA = rrna_count, tRNA = trna_count)

    cat("\n=== ANNOTATION RESULTS ===\n")
    cat(sprintf("CDS:   %d\n", stats$cds))
    cat(sprintf("rRNA:  %d\n", stats$rRNA))
    cat(sprintf("tRNA:  %d\n\n", stats$tRNA))

    return(stats)
  } else {
    cat("Warning: one or more annotation output files not found.\n")
    return(NULL)
  }
}

#' Generate final text report
generate_report <- function(metrics, annotation) {
  cat("===================================================================\n")
  cat("STEP 8: Generating Report\n")
  cat("===================================================================\n\n")

  report_file <- file.path(CONFIG$base_dir, "reports", "ANALYSIS_REPORT.txt")

  sink(report_file)
  cat("===================================================================\n")
  cat("       GENOME ASSEMBLY REPORT\n")
  cat("===================================================================\n\n")
  cat(sprintf("Dataset: %s\n", CONFIG$sra_accession))
  cat(sprintf("Date: %s\n\n", Sys.Date()))

  cat("ASSEMBLY QUALITY:\n")
  cat(sprintf("  Contigs:        %d\n", metrics$contigs))
  cat(sprintf("  Total length:   %.2f Mb\n", metrics$length / 1e6))
  cat(sprintf("  Largest contig: %d bp\n", metrics$largest))
  cat(sprintf("  N50:            %.1f kb\n", metrics$n50 / 1000))
  cat(sprintf("  GC content:     %.2f%%\n\n", metrics$gc))

  if (!is.null(annotation)) {
    cat("ANNOTATION:\n")
    cat(sprintf("  CDS:   %d\n", annotation$cds))
    cat(sprintf("  rRNA:  %d\n", annotation$rRNA))
    cat(sprintf("  tRNA:  %d\n\n", annotation$tRNA))
  }

  cat("===================================================================\n")
  sink()

  cat("Report saved to:", report_file, "\n\n")
}

# ============================================================================
# INTERACTIVE FUNCTIONS FOR RSTUDIO
# ============================================================================

#' Run complete pipeline interactively
#' @export
run_pipeline_interactive <- function() {
  cat("\nStarting Interactive Pipeline\n\n")

  tryCatch({
    fastq_files <- download_data()
    contigs     <- run_assembly(fastq_files)
    metrics     <- run_quality(contigs)
    annotation  <- run_annotation(contigs)
    generate_report(metrics, annotation)

    cat("\n")
    cat("================================================================\n")
    cat("              ANALYSIS COMPLETE                                \n")
    cat("================================================================\n\n")

    cat("Results:\n")
    cat(sprintf("  Contigs:       %d\n", metrics$contigs))
    cat(sprintf("  N50:           %.1f kb\n", metrics$n50 / 1000))
    cat(sprintf("  Genome Size:   %.2f Mb\n", metrics$length / 1e6))
    if (!is.null(annotation)) {
      cat(sprintf("  Genes (CDS):   %d\n", annotation$cds))
    }
    cat(sprintf("\n  Output folder: %s/\n\n", CONFIG$base_dir))

    cat("Next step: run scripts/02_generate_figures.R from inside\n")
    cat(sprintf("  %s/\n", CONFIG$base_dir))
    cat("to generate the publication figures (Figures 1-5).\n\n")

  }, error = function(e) {
    cat("\nError:", e$message, "\n")
    cat("\nTroubleshooting:\n")
    cat("1. Make sure all tools are installed\n")
    cat("2. Check your internet connection\n")
    cat("3. Confirm SRA Toolkit (fastq-dump) is in PATH\n\n")
  })
}

# ============================================================================
# DISPLAY INSTRUCTIONS
# ============================================================================

cat("\n")
cat("================================================================\n")
cat("              READY TO RUN IN RSTUDIO                          \n")
cat("================================================================\n\n")

cat("WHAT TO DO NEXT:\n")
cat("-------------------------------------------------------------------\n")
cat("\nIn RStudio Console, type:\n\n")
cat("   run_pipeline_interactive()\n\n")
cat("Or click 'Source' to run automatically when sourcing this file.\n\n")

# Auto-run if sourced non-interactively (e.g., Rscript)
if (!interactive()) {
  cat("Auto-running pipeline...\n\n")
  run_pipeline_interactive()
}
