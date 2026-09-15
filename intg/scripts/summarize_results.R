suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

cmd_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", cmd_args, value = TRUE)
if (length(file_arg) != 1) stop("Could not determine script path")
script_path <- normalizePath(sub("^--file=", "", file_arg))
intg_dir <- dirname(dirname(script_path))
setwd(intg_dir)

fin <- file.path("output", "gdi_intg_long.csv")
if (!file.exists(fin)) stop("Missing ", fin, ". Run the gdi analysis first.")

raw <- read_csv(fin, show_col_types = FALSE)
chromosomes <- paste0("ch", 1:9)

results <- raw |>
  mutate(
    chromosome = factor(chromosome, levels = chromosomes),
    series = case_when(
      comparison == "intg_N_vs_CA" & focal_population == "intg_nGL_Nslope" ~ "N relative to CA",
      comparison == "intg_N_vs_CA" & focal_population == "intg_CAswGL" ~ "CA relative to N",
      comparison == "CA_vs_hook" & focal_population == "intg_CAswGL" ~ "CA relative to hook",
      comparison == "CA_vs_hook" & focal_population == "hook" ~ "hook relative to CA",
      TRUE ~ NA_character_
    ),
    series = factor(
      series,
      levels = c("N relative to CA", "CA relative to N", "CA relative to hook", "hook relative to CA")
    )
  )

wide <- results |>
  mutate(column = case_when(
    comparison == "intg_N_vs_CA" & focal_population == "intg_nGL_Nslope" ~ "gdi_intg_N_vs_CA__intg_nGL_Nslope",
    comparison == "intg_N_vs_CA" & focal_population == "intg_CAswGL" ~ "gdi_intg_N_vs_CA__intg_CAswGL",
    comparison == "CA_vs_hook" & focal_population == "intg_CAswGL" ~ "gdi_CA_vs_hook__intg_CAswGL",
    comparison == "CA_vs_hook" & focal_population == "hook" ~ "gdi_CA_vs_hook__hook",
    TRUE ~ NA_character_
  )) |>
  select(chromosome, column, gdi) |>
  pivot_wider(names_from = column, values_from = gdi) |>
  arrange(chromosome) |>
  mutate(chromosome = as.character(chromosome))

write_csv(wide, file.path("output", "gdi_intg.csv"))

summary_table <- results |>
  group_by(series) |>
  summarise(
    mean_gdi = mean(gdi),
    median_gdi = median(gdi),
    min_gdi = min(gdi),
    max_gdi = max(gdi),
    .groups = "drop"
  )
write_csv(summary_table, file.path("output", "gdi_intg_summary.csv"))

p <- ggplot(results, aes(x = series, y = gdi)) +
  geom_violin(trim = FALSE, alpha = 0.25) +
  geom_jitter(width = 0.08, height = 0, size = 2) +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 3, fill = "white") +
  geom_hline(yintercept = 0.2, linetype = "dashed") +
  geom_hline(yintercept = 0.7, linetype = "dashed") +
  coord_cartesian(ylim = c(0, 1)) +
  labs(
    x = NULL,
    y = "gdi",
    title = "Genealogical divergence across chromosomes",
    subtitle = "Points are chromosome-specific estimates; diamonds show means"
  ) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 25, hjust = 1),
    panel.grid.minor = element_blank()
  )

dir.create(file.path("output", "figures"), recursive = TRUE, showWarnings = FALSE)
ggsave(file.path("output", "figures", "gdi_intg_violin.pdf"), p, width = 8, height = 5)
ggsave(file.path("output", "figures", "gdi_intg_violin.png"), p, width = 8, height = 5, dpi = 300)

print(wide)
print(summary_table)
message("Wrote corrected gdi_intg.csv, gdi_intg_summary.csv, and violin plots in output/figures/")
