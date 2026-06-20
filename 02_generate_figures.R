# ============================================================================
# FIGURE GENERATION: SRR25083113 GENOME ASSEMBLY & ANNOTATION
# ============================================================================
# Run this AFTER scripts/01_genome_assembly_pipeline.R has completed.
#
# Working directory must be the analysis output folder, i.e.:
#   setwd("<path-to>/SRR25083113_analysis")
# which contains: assembly/contigs.fasta, annotation/proteins.faa,
# annotation/rRNAs.gff, annotation/tRNAs.txt
#
# Install once if needed:
#   install.packages(c("tidyverse", "cowplot"))
#   BiocManager::install("Biostrings")
# ============================================================================

library(tidyverse)
library(cowplot)
library(Biostrings)

plots_dir <- "plots"
if (!dir.exists(plots_dir)) dir.create(plots_dir, recursive = TRUE)

# ============================================================================
# FIGURE 1: Assembly Quality Summary
# Panel A: contig length distribution | Panel B: GC content vs contig length
# ============================================================================

contigs <- Biostrings::readDNAStringSet("assembly/contigs.fasta")

contig_df <- tibble(
  contig = names(contigs),
  length = width(contigs),
  gc     = letterFrequency(contigs, c("G", "C"), as.prob = TRUE) %>% rowSums() * 100
)

# Sanity check against report metrics (5,420,912 bp / 325,905 bp largest)
sum(contig_df$length)
max(contig_df$length)
median(contig_df$length)

p_len_hist <- ggplot(contig_df, aes(x = length)) +
  geom_histogram(bins = 50) +
  scale_x_log10() +
  labs(
    title = "Contig Length Distribution",
    x = "Contig length (bp, log10 scale)",
    y = "Count"
  ) +
  theme_minimal(base_size = 14)

p_gc_vs_len <- ggplot(contig_df, aes(x = length, y = gc)) +
  geom_point(alpha = 0.5, size = 1) +
  scale_x_log10() +
  labs(
    title = "GC Content vs Contig Length",
    x = "Contig length (bp, log10 scale)",
    y = "GC (%)"
  ) +
  theme_minimal(base_size = 14)

fig1_assembly <- plot_grid(p_len_hist, p_gc_vs_len, ncol = 2, labels = c("A", "B"))
ggsave(file.path(plots_dir, "fig1_assembly_qc.png"), fig1_assembly,
       width = 10, height = 5, dpi = 300)

# ============================================================================
# FIGURE 2: Cumulative Assembly Curve
# ============================================================================

contig_sorted <- contig_df %>%
  arrange(desc(length)) %>%
  mutate(
    index    = row_number(),
    cum_len  = cumsum(length),
    cum_frac = cum_len / sum(length)
  )

p_cum <- ggplot(contig_sorted, aes(x = index, y = cum_frac)) +
  geom_step() +
  labs(
    title = "Cumulative Assembly Curve",
    x = "Contig rank (largest to smallest)",
    y = "Cumulative genome fraction"
  ) +
  theme_minimal(base_size = 14)

ggsave(file.path(plots_dir, "fig2_cumulative_assembly.png"), p_cum,
       width = 6, height = 4, dpi = 300)

# ============================================================================
# Shared annotation parsing (used by Figures 3 and 4)
# ============================================================================

proteins <- Biostrings::readAAStringSet("annotation/proteins.faa")
cds_count <- length(proteins)

rrna_gff <- read.table(
  "annotation/rRNAs.gff",
  sep = "\t", header = FALSE, comment.char = "#", stringsAsFactors = FALSE
)
colnames(rrna_gff) <- c("seqid", "source", "type", "start", "end",
                         "score", "strand", "phase", "attributes")
rrna_count <- nrow(rrna_gff)

trna_lines <- readLines("annotation/tRNAs.txt")
trna_count <- sum(grepl("tRNA", trna_lines, ignore.case = TRUE))

# Sanity check against report metrics (CDS=5513, tRNA=103, rRNA=9)
cds_count; rrna_count; trna_count

# ============================================================================
# FIGURE 3: Genome Feature Composition
# Panel A: counts | Panel B: proportions
# ============================================================================

feature_df <- tibble(
  Feature = c("CDS", "tRNA", "rRNA"),
  Count   = c(cds_count, trna_count, rrna_count)
)

p_feat_bar <- ggplot(feature_df, aes(x = Feature, y = Count, fill = Feature)) +
  geom_col(show.legend = FALSE) +
  labs(title = "Genome Feature Composition", x = NULL, y = "Count") +
  theme_minimal(base_size = 14)

p_feat_pie <- ggplot(feature_df, aes(x = "", y = Count, fill = Feature)) +
  geom_col(color = "white") +
  coord_polar(theta = "y") +
  labs(title = "Genome Feature Composition") +
  theme_void(base_size = 14)

fig3_features <- plot_grid(p_feat_bar, p_feat_pie, ncol = 2, labels = c("A", "B"))
ggsave(file.path(plots_dir, "fig3_gene_features.png"), fig3_features,
       width = 10, height = 5, dpi = 300)

# ============================================================================
# FIGURE 4: rRNA Feature Distribution (16S / 23S / 5S) - clean/high-contrast
# ============================================================================

rrna_gff$rRNA_type <- rrna_gff$attributes %>%
  str_extract("(16S|23S|5S)_rRNA") %>%
  str_replace("_rRNA", "")

rrna_counts <- rrna_gff %>%
  count(rRNA_type) %>%
  arrange(desc(n))

print(rrna_counts)

p_rrna <- ggplot(rrna_counts, aes(x = reorder(rRNA_type, n), y = n, fill = rRNA_type)) +
  geom_col(width = 0.65, color = "white") +
  coord_flip() +
  scale_fill_manual(values = c("16S" = "#2E86AB",
                                "23S" = "#F18F01",
                                "5S"  = "#9BC53D")) +
  labs(
    title = "Distribution of rRNA Features",
    subtitle = "Counts of 16S, 23S, and 5S rRNA genes detected by Barrnap",
    x = "rRNA Type",
    y = "Count"
  ) +
  theme_minimal(base_size = 15) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "none",
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  )

ggsave(file.path(plots_dir, "fig4_rrna_types_clean.png"), p_rrna,
       width = 7, height = 5, dpi = 300)

# ============================================================================
# FIGURE 5: tRNA Gene Distribution Across Top 20 Contigs - high-contrast
# ============================================================================

current_contig <- NA_character_
trna_tbl <- tibble()

for (ln in trna_lines) {
  if (str_starts(ln, ">")) {
    current_contig <- str_remove(ln, "^>")
  } else if (str_detect(ln, "tRNA")) {
    aa        <- str_match(ln, "tRNA-([A-Za-z]+)")[, 2]
    anticodon <- str_match(ln, "\\(([a-z]{3})\\)")[, 2]
    coords    <- str_match(ln, "\\[([0-9]+),([0-9]+)\\]")

    trna_tbl <- trna_tbl %>% add_row(
      contig    = current_contig,
      aa        = aa,
      anticodon = anticodon,
      start     = as.integer(coords[, 2]),
      end       = as.integer(coords[, 3])
    )
  }
}

trna_tbl$length <- abs(trna_tbl$end - trna_tbl$start) + 1

top20 <- trna_tbl %>%
  count(contig) %>%
  slice_max(n, n = 20)

aa_palette <- c(
  "Ala" = "#0072B2", "Arg" = "#D55E00", "Asn" = "#009E73", "Asp" = "#CC79A7",
  "Cys" = "#F0E442", "Gln" = "#56B4E9", "Glu" = "#E69F00", "Gly" = "#000000",
  "His" = "#009E73", "Ile" = "#0072B2", "Leu" = "#D55E00", "Lys" = "#F0E442",
  "Met" = "#56B4E9", "Phe" = "#CC79A7", "Pro" = "#000000", "Ser" = "#009E73",
  "Thr" = "#E69F00", "Trp" = "#0072B2", "Tyr" = "#D55E00", "Val" = "#56B4E9"
)

p_trna <- ggplot(
  trna_tbl %>% filter(contig %in% top20$contig),
  aes(x = start, y = reorder(contig, start), color = aa)
) +
  geom_segment(
    aes(xend = end, yend = reorder(contig, start)),
    linewidth = 2.5, alpha = 0.95
  ) +
  scale_color_manual(values = aa_palette) +
  labs(
    title = "tRNA Gene Distribution Across Top 20 Contigs",
    subtitle = "High-contrast color palette for clear amino acid interpretation",
    x = "Genomic Position (bp)",
    y = "Contig",
    color = "Amino Acid"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 16),
    legend.title = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank()
  )

ggsave(file.path(plots_dir, "fig5_trna_per_contig_highcontrast.png"), p_trna,
       width = 12, height = 7, dpi = 300)

cat("\nAll figures saved to:", plots_dir, "\n")
cat("  fig1_assembly_qc.png\n")
cat("  fig2_cumulative_assembly.png\n")
cat("  fig3_gene_features.png\n")
cat("  fig4_rrna_types_clean.png\n")
cat("  fig5_trna_per_contig_highcontrast.png\n")
