#!usr/bin/env Rscript


# Import libraries
library(Biostrings)
library(tidyverse)
library(cowplot)
library(patchwork)
library(fastqcr)
library(RColorBrewer)
library(ggpubr)
library(ggbeeswarm)

# Import auxiliary scripts
source("../scripts/aux_plot-theme.R")
source("../scripts/aux_functions.R")

# Scales and palettes
scale_fill_batch  <- purrr::partial(scale_fill_brewer, name = "Sample Group",
                                      palette = "Set1")
scale_color_batch <- purrr::partial(scale_color_brewer, name = "Sample Group",
                                      palette = "Set1")
scale_shape_batch <- purrr::partial(scale_shape_discrete, name = "Sample Group")
scale_fill_ribo <- purrr::partial(scale_fill_brewer, name="Ribodepletion",
                                palette="Dark2")
scale_color_ribo <- purrr::partial(scale_color_brewer, name="Ribodepletion",
                                palette="Dark2")
scale_shape_ribo <- purrr::partial(scale_shape_discrete, name="Ribodepletion")
scale_color_delivery <- purrr::partial(scale_color_brewer, palette="Dark2", name="Delivery")
scale_fill_delivery <- purrr::partial(scale_fill_brewer, palette="Dark2", name="Delivery")

#==============================================================================
# 1. Set input directories
#==============================================================================

# Set input dir
data_dir <- "../data"

# Structured subsubdirectories
log_dir <- file.path(data_dir, "logging")
input_dir <- file.path(data_dir, "input")
results_dir <- file.path(data_dir, "results")

#==============================================================================
# 2. Prepare metadata and groups
#==============================================================================

# Get metadata paths
metadata_path <- file.path(data_dir, "bio_sample_table.tsv")

# Import metadata
metadata <- read_tsv(metadata_path, show_col_types = FALSE) %>%
  mutate(sample = `Sample Name`) %>%
  select(-`Sample Name`)

#==============================================================================
# 3. Import and clean QC data
#==============================================================================

# Data input paths
basic_stats_paths <- file.path(results_dir, "qc_basic_stats.tsv")
quality_base_stats_paths <- file.path(results_dir, "qc_quality_base_stats.tsv")

# Import QC data
stages <- c("raw_concat", "cleaned")

basic_stats <- read_tsv(basic_stats_paths, show_col_types = FALSE) %>%
  inner_join(metadata, by="sample") %>%
  mutate(stage = factor(stage, levels = stages),
         sample = fct_inorder(sample))

quality_base_stats <- read_tsv(quality_base_stats_paths, show_col_types = FALSE) %>%
  inner_join(metadata, by="sample") %>%
  mutate(stage = factor(stage, levels = stages),
         sample = fct_inorder(sample),
         date = as.Date(str_extract(sample, "\\d{4}-\\d{2}-\\d{2}"), format="%Y-%m-%d"),
         sequencer = ifelse(date >= as.Date("2024-02-25"), "NovaSeq X", "NovaSeq 6000"),
         read_pair = str_split(file, "_") %>% sapply(dplyr::last))

basic_stats_raw <- basic_stats %>% filter(stage == "raw_concat")

# Count reads
read_counts_total <- basic_stats_raw %>%
  summarize(
    rmin = min(n_read_pairs),
    rmax = max(n_read_pairs),
    rmean = mean(n_read_pairs),
    rtot = sum(n_read_pairs),
    btot = sum(n_bases_approx),
    .groups = "drop"
  )

print(read_counts_total)
basic_stats_raw_metrics <- basic_stats_raw %>%
  mutate(n_read_pairs = replace_na(n_read_pairs, 0),
         n_bases_approx = replace_na(n_bases_approx, 0),
         date = as.Date(str_extract(sample, "\\d{4}-\\d{2}-\\d{2}"), format="%Y-%m-%d"),
         sequencer = ifelse(date >= as.Date("2024-02-25"), "NovaSeq X", "NovaSeq 6000"),
         date_display = format(date, "%Y-%m-%d"),
         date_display = factor(date_display, levels = unique(date_display[order(date)]))) %>%
  select(date_display,
         sequencer,
         `# Read pairs` = n_read_pairs,
         `Total base pairs` = n_bases_approx) %>%
  # First gather the numeric columns only
  pivot_longer(cols = c(`# Read pairs`, `Total base pairs`),
               names_to = "metric",
               values_to = "value") %>%
  mutate(metric = fct_inorder(metric))

# Set up plot using the Dark2 palette
g_basic <- ggplot(basic_stats_raw_metrics,
                 aes(x = date_display, y = value, fill = sequencer)) +
  geom_col(width = 0.3) +
  scale_x_discrete(name = "Date") +
  scale_y_continuous(expand = c(0,0),
                    limits = function(x) c(0, max(x) * 1.1)) +  # Set upper limit to 10% above max value
  scale_fill_brewer(palette = "Dark2", name = "Sequencing\nMachine") +  # Using Dark2 palette for sequencer colors
  facet_grid(metric ~ ., scales = "free", space = "free_x", switch = "y") +
  theme_tilt +
  theme(
    axis.title.y = element_blank(),
    strip.text.y = element_text(face = "plain"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
g_basic

ggsave("fig_1.png", g_basic, width = 10, height = 5)

#==============================================================================
# 4. Figure 2: Quality metrics
#==============================================================================

g_qual <- ggplot(mapping = aes(color = sequencer, linetype = read_pair,
                              group = interaction(sample, read_pair))) +
  scale_linetype_discrete(name = "Read Pair") +
  scale_color_brewer(palette = "Dark2", name = "Sequencing\nMachine") +
  guides(color = guide_legend(nrow = 1, byrow = TRUE),
         linetype = guide_legend(nrow = 1, byrow = TRUE)) +
  theme_base


quality_base_stats_raw <- quality_base_stats %>% filter(stage == "raw_concat")

g_quality_base <- g_qual +
  geom_hline(yintercept=25, linetype="dashed", color="red") +
  geom_hline(yintercept=30, linetype="dashed", color="red") +
  geom_line(aes(x=position, y=mean_phred_score), data=quality_base_stats_raw) +
  scale_y_continuous(name="Mean Phred score", expand=c(0,0), limits=c(10,45)) +
  scale_x_continuous(name="Position", limits=c(0,NA),
                     breaks=seq(0,140,20), expand=c(0,0))

g_quality_base

ggsave("fig_2.png", g_quality_base, width = 8, height = 3)
