#!/usr/bin/env Rscript
## =====================================================================
## make_boxplots.R  --  build the Figure-2 replacement panels
## ---------------------------------------------------------------------
## Reads results/grid_<scenario>_B<BB>.csv (from run_grid_array.R) and writes,
## in figures/, ONE self-contained figure per scenario laid out exactly as:
##
##        |            Coverage prob.            |          Avg. length          |
##   B=1  | COPP tCOPP-W EIF(I) EIF(M) EIF(P)    | COPP tCOPP-W EIF(I) EIF(M) EIF(P)
##   B=20 | COPP tCOPP-W EIF(I) EIF(M) EIF(P)    | COPP tCOPP-W EIF(I) EIF(M) EIF(P)
##
##   -> 5 CP boxes + 5 AL boxes per B row, two B rows  = 20 boxes / scenario
##   -> 3 scenarios                                    = 60 boxes total
##
## Files:
##   fig_scenario_a.(png|pdf), fig_scenario_b.(png|pdf), fig_scenario_c.(png|pdf)
##       -> drop these straight into the (a)/(b)/(c) \includegraphics slots
##   fig_MASTER.(png|pdf)          -> all three scenarios stacked (one image)
##   results_combined.csv, grid_summary_medians.csv
##
## Boxes are coloured by source-mode so the "internal source only and multiple
## sources" caption is visible at a glance:
##   COPP / tCOPP-W / tCOPP-EIF(I) = Internal only
##   tCOPP-EIF(M) / tCOPP-EIF(P)   = Multiple sources
## =====================================================================

suppressMessages({ library(ggplot2); library(dplyr); library(tidyr) })

in_dir  <- "results"
fig_dir <- "figures"
if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)

method_levels <- c("COPP", "tCOPP-W", "tCOPP-EIF(I)", "tCOPP-EIF(M)", "tCOPP-EIF(P)")
metric_levels <- c("Coverage probability", "Average length")
scen_labels   <- c(a = "(a) CCOD holds for all sources",
                   b = "(b) CCOD violated for external sources",
                   c = "(c) CCOD violated (complex)")
mode_fill     <- c("Internal only" = "#8FB8DE", "Multiple sources" = "#E8A87C")
alpha_nominal <- 0.9

## ---- optional layout packages (graceful fallback if absent) ---------
have_patch <- requireNamespace("patchwork", quietly = TRUE)
have_grid  <- requireNamespace("gridExtra", quietly = TRUE)

## ---- load available blocks -----------------------------------------
files <- list.files(in_dir, pattern = "^grid_[abc]_B[0-9]+\\.csv$", full.names = TRUE)
if (length(files) == 0) stop("No grid_*.csv in '", in_dir, "'. Run run_grid_array.R first.")
dat <- bind_rows(lapply(files, read.csv, stringsAsFactors = FALSE))
dat <- dat[dat$method %in% method_levels, ]
dat$mode_lab <- ifelse(dat$mode == "internal", "Internal only", "Multiple sources")
write.csv(dat, file.path(fig_dir, "results_combined.csv"), row.names = FALSE)
message(sprintf("Loaded %d rows across %d block(s).", nrow(dat), length(files)))

## long form (CP + AL stacked into one 'value' column)
long <- dat %>%
  pivot_longer(c(CP, AL), names_to = "metric_raw", values_to = "value") %>%
  mutate(metric  = factor(ifelse(metric_raw == "CP", "Coverage probability", "Average length"),
                          levels = metric_levels),
         method  = factor(method, levels = method_levels),
         mode_lab= factor(mode_lab, levels = names(mode_fill)),
         Blab    = factor(paste0("B = ", B), levels = paste0("B = ", sort(unique(B)))))

base_theme <- theme_bw(base_size = 11) +
  theme(axis.text.x  = element_text(angle = 45, hjust = 1, face = "bold", size = 8),
        plot.title   = element_text(hjust = 0.5, size = 12),
        strip.text   = element_text(face = "bold"),
        legend.position = "top", legend.title = element_blank())

## one CP-or-AL panel for a given scenario & B
panel <- function(dsB, metric_name, show_y_title = TRUE) {
  d <- dsB %>% filter(metric == metric_name)
  g <- ggplot(d, aes(method, value, fill = mode_lab)) +
    geom_boxplot(outlier.size = 0.5) +
    scale_fill_manual(values = mode_fill, drop = FALSE) +
    labs(x = NULL, y = if (show_y_title) metric_name else NULL) +
    base_theme
  if (metric_name == "Coverage probability")
    g <- g + geom_hline(yintercept = alpha_nominal, linewidth = 0.6)
  g
}

## assemble the 2x2 (B rows x metric cols) figure for one scenario
scenario_figure <- function(sc) {
  ds <- long %>% filter(scenario == sc)
  Bs <- sort(unique(ds$B))
  ttl <- scen_labels[sc]

  if (have_patch) {
    library(patchwork)
    panels <- list()
    for (b in Bs) {
      dB <- ds %>% filter(B == b)
      cp <- panel(dB, "Coverage probability") + ggtitle(paste0("B = ", b))
      al <- panel(dB, "Average length")        + ggtitle(paste0("B = ", b))
      panels <- c(panels, list(cp, al))
    }
    wrap_plots(panels, ncol = 2, byrow = TRUE) +
      plot_layout(guides = "collect") +
      plot_annotation(title = ttl,
                      theme = theme(plot.title = element_text(hjust = 0.5, size = 13),
                                    legend.position = "top"))
  } else {
    ## fallback: single ggplot, facet_grid(B ~ metric) with free scales
    ggplot(ds, aes(method, value, fill = mode_lab)) +
      geom_boxplot(outlier.size = 0.5) +
      geom_hline(data = data.frame(metric = factor("Coverage probability", levels = metric_levels),
                                   y = alpha_nominal),
                 aes(yintercept = y), linewidth = 0.6, inherit.aes = FALSE) +
      scale_fill_manual(values = mode_fill, drop = FALSE) +
      facet_grid(Blab ~ metric, scales = "free_y") +
      labs(x = NULL, y = NULL, title = ttl) +
      base_theme
  }
}

## ---- write the three scenario figures ------------------------------
for (sc in intersect(c("a","b","c"), unique(long$scenario))) {
  fig <- scenario_figure(sc)
  for (ext in c("png","pdf"))
    ggsave(file.path(fig_dir, paste0("fig_scenario_", sc, ".", ext)),
           fig, width = 9, height = 6.5, dpi = 220)
  message("wrote fig_scenario_", sc)
}

## ---- combined master (all three scenarios stacked) -----------------
if (have_patch) {
  library(patchwork)
  figs <- lapply(intersect(c("a","b","c"), unique(long$scenario)), scenario_figure)
  master <- Reduce(`/`, figs)
  for (ext in c("png","pdf"))
    ggsave(file.path(fig_dir, paste0("fig_MASTER.", ext)), master,
           width = 9, height = 6.5 * length(figs), dpi = 200, limitsize = FALSE)
  message("wrote fig_MASTER")
} else if (have_grid) {
  library(gridExtra)
  figs <- lapply(intersect(c("a","b","c"), unique(long$scenario)), scenario_figure)
  g <- do.call(gridExtra::arrangeGrob, c(figs, ncol = 1))
  for (ext in c("png","pdf"))
    ggsave(file.path(fig_dir, paste0("fig_MASTER.", ext)), g,
           width = 9, height = 6.5 * length(figs), dpi = 200, limitsize = FALSE)
  message("wrote fig_MASTER (gridExtra)")
} else {
  message("patchwork/gridExtra not installed -> skipped combined master; per-scenario figures written.")
}

## ---- median summary table ------------------------------------------
summ <- dat %>% group_by(scenario, B, method, mode) %>%
  summarise(CP_med = round(median(CP, na.rm = TRUE), 3),
            AL_med = round(median(AL, na.rm = TRUE), 3),
            n = sum(!is.na(CP)), .groups = "drop")
summ$method <- factor(summ$method, levels = method_levels)
summ <- summ[order(summ$scenario, summ$B, summ$method), ]
write.csv(summ, file.path(fig_dir, "grid_summary_medians.csv"), row.names = FALSE)
print(as.data.frame(summ))
message("\nDone. Figures + summaries in '", fig_dir, "/'.")
