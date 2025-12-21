# make_plots_wager.R
# Plots performance of EP / DR / R / T learners from CSVs written by EPlearnerCATE_wager.R

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(stringr)
})
 BASE_PATH <- Sys.getenv("BASE_PATH", unset = NA_character_)
if (is.na(BASE_PATH) || !nzchar(BASE_PATH)) {
  this_file <- tryCatch(normalizePath(sys.frame(1)$ofile, winslash = "/"), error = function(e) NA_character_)
  if (is.na(this_file)) stop("Set BASE_PATH env var, e.g. Sys.setenv(BASE_PATH='/path/to/project_root')")
  BASE_PATH <- normalizePath(dirname(this_file), winslash = "/")
} else {
  BASE_PATH <- normalizePath(path.expand(BASE_PATH), winslash = "/")
}


 RESULTS_DIR <- file.path(BASE_PATH, "experiment_results")  # <- new convention
 FIGURES_DIR <- file.path(BASE_PATH, "experiment_results")      # keep existing relative structure
 PLOTS_DIR   <- file.path(FIGURES_DIR, "plots")

 dir.create(PLOTS_DIR,   recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# Helper: read & combine all CSVs
# ------------------------------------------------------------
# Expects files named like:
# simsEP_n=2000_d=5_sigma=1_base_learner=xg_sim_type=A_nsims=100.csv
all_csvs <- list.files(RESULTS_DIR, pattern = "^simsEP_.*\\.csv$", full.names = TRUE)
if (length(all_csvs) == 0L) {
  stop(
    "No result CSVs found in: ", RESULTS_DIR,
    "\nDid you run the sims and write the CSV outputs?"
  )
}

read_one <- function(f) {
  dt <- tryCatch(fread(f), error = function(e) NULL)
  if (is.null(dt) || nrow(dt) == 0L) return(NULL)

  # file-derived fallback fields if missing
  info <- str_match(
    basename(f),
    "n=([0-9]+)_d=([0-9]+)_sigma=([0-9\\.]+)_base_learner=([a-z]+)_sim_type=([A-D])_nsims=([0-9]+)"
  )
  if (!"n"        %in% names(dt)) dt[, n := as.numeric(info[2])]
  if (!"d"        %in% names(dt)) dt[, d := as.numeric(info[3])]
  if (!"sigma"    %in% names(dt)) dt[, sigma := as.numeric(info[4])]
  if (!"lrnr"     %in% names(dt)) dt[, lrnr := info[5]]      # base_learner
  if (!"sim_type" %in% names(dt)) dt[, sim_type := info[6]]
  dt[]
}

res_list <- lapply(all_csvs, read_one)
res_list <- res_list[!vapply(res_list, is.null, logical(1))]
results  <- rbindlist(res_list, use.names = TRUE, fill = TRUE)
if (nrow(results) == 0L) stop("Parsed 0 rows from results CSVs. Check file names/format.")

# ------------------------------------------------------------
# Remove iterations where DR-Learner explodes (then drop across methods)
# ------------------------------------------------------------
abs_cut <- 1e3   # absolute cutoff for dr_risk
rel_cut <- 50    # relative cutoff vs cell median

by_keys <- c("sim_type", "lrnr", "n")

dr_ref <- results[, .(dr_median = median(dr_risk, na.rm = TRUE)), by = by_keys]
res2 <- merge(results, dr_ref, by = by_keys, all.x = TRUE)

res2[, dr_bad := (!is.na(dr_risk)) & (
  dr_risk > abs_cut |
    (is.finite(dr_median) & dr_median > 0 & dr_risk > rel_cut * dr_median)
)]

bad_iters <- unique(res2[dr_bad == TRUE, .(sim_type, lrnr, n, iter, dr_risk)])
if (nrow(bad_iters)) {
  fwrite(bad_iters, file.path(PLOTS_DIR, "dropped_iters_dr_excessive.csv"))
  message(sprintf(
    "Dropping %d iters where DR risk exceeded thresholds. Logged to dropped_iters_dr_excessive.csv",
    nrow(bad_iters)
  ))
}

results <- res2[dr_bad == FALSE, !(c("dr_median", "dr_bad")), with = FALSE]

# ------------------------------------------------------------
# Tidy to long format (EP/DR/R/T + CF when present)
# ------------------------------------------------------------
meas <- intersect(
  c("ep_risk", "dr_risk", "r_risk", "t_risk", "causal_forest_risk"),
  names(results)
)

long <- melt(
  results,
  id.vars = c("iter", "sim_type", "n", "d", "sigma", "lrnr"),
  measure.vars = meas,
  variable.name = "method",
  value.name = "mse"
)

method_map <- c(
  ep_risk            = "EP-Learner (ours)",
  dr_risk            = "DR-Learner",
  r_risk             = "R-Learner",
  t_risk             = "T-Learner",
  causal_forest_risk = "Causal Forest"
)
long[, method := factor(
  method_map[method],
  levels = c("EP-Learner (ours)", "DR-Learner", "R-Learner", "T-Learner", "Causal Forest")
)]

# keep CF only for rf; drop NA rows
long <- long[!(method == "Causal Forest" & lrnr != "rf") & !is.na(mse) & !is.na(method)]
long[, method := droplevels(method)]

blr_map <- c(xg = "XGBoost", rf = "Random Forest", gam = "GAM")
long[, base_learner := factor(blr_map[lrnr], levels = c("XGBoost", "Random Forest", "GAM"))]

long[, `:=`(n = as.numeric(n), sigma = as.numeric(sigma))]
long[, sim_type := factor(sim_type, levels = c("A", "B", "C", "D"))]

# ------------------------------------------------------------
# Aggregate across iterations
# ------------------------------------------------------------
summ <- long[
  , .(
    mean_mse = mean(mse, na.rm = TRUE),
    se_mse   = sd(mse,  na.rm = TRUE) / sqrt(sum(!is.na(mse))),
    n_iter   = sum(!is.na(mse))
  ),
  by = .(sim_type, base_learner, n, method)
][n_iter > 0]

summ[, `:=`(
  lo = pmax(mean_mse - 1.96 * se_mse, .Machine$double.eps),
  hi = mean_mse + 1.96 * se_mse
)]
summ <- summ[!is.na(base_learner)]

# ------------------------------------------------------------
# Tables mirroring the plot
# ------------------------------------------------------------
setkeyv(summ, c("sim_type", "base_learner", "n", "method"))

tbl_long <- summ[, .(sim_type, base_learner, n, method, mean_mse, se_mse, n_iter)]
fwrite(tbl_long, file.path(PLOTS_DIR, "performance_table_long.csv"))

tbl_wide <- dcast(tbl_long, sim_type + base_learner + n ~ method, value.var = "mean_mse")
setcolorder(tbl_wide, c("sim_type", "base_learner", "n",
                        "EP-Learner (ours)", "DR-Learner", "R-Learner", "T-Learner"))
fwrite(tbl_wide, file.path(PLOTS_DIR, "performance_table_wide.csv"))

fmt_num <- function(x) formatC(x, format = "e", digits = 2)
tbl_pretty_src <- copy(summ)[
  , .(cell = sprintf("%s (%s) [%d]", fmt_num(mean_mse), fmt_num(se_mse), n_iter)),
  by = .(sim_type, base_learner, n, method)
]
tbl_wide_pretty <- dcast(tbl_pretty_src, sim_type + base_learner + n ~ method, value.var = "cell")
setcolorder(tbl_wide_pretty, c("sim_type", "base_learner", "n",
                               "EP-Learner (ours)", "DR-Learner", "R-Learner", "T-Learner"))
fwrite(tbl_wide_pretty, file.path(PLOTS_DIR, "performance_table_wide_pretty.csv"))

# ------------------------------------------------------------
# Plot helpers
# ------------------------------------------------------------
base_theme <- theme_bw() +
  theme(
    axis.text  = element_text(size = 12),
    axis.title = element_text(size = 13, face = "bold"),
    strip.text = element_text(size = 12, face = "bold"),
    legend.title = element_text(size = 12),
    legend.text  = element_text(size = 11)
  )

make_plot <- function(dt, title = NULL) {
  ggplot(dt, aes(x = n, y = mean_mse, color = method, linetype = method, group = method)) +
    geom_line() +
    geom_point(size = 1.8) +
    scale_y_log10() +
    scale_x_log10(breaks = sort(unique(dt$n))) +
    labs(
      x = "Sample Size (n)", y = "Mean Squared Error (MSE)",
      color = "Method", linetype = "Method", title = title
    ) +
    base_theme +
    theme(
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.box = "horizontal"
    )
}

# ------------------------------------------------------------
# 1) Full grid: sim_type (rows) × base_learner (cols)
# ------------------------------------------------------------
p_all <- make_plot(summ) +
  facet_grid(sim_type ~ base_learner, scales = "free_y") +
  theme(panel.spacing.x = unit(1.2, "cm"))

ggsave(file.path(PLOTS_DIR, "performance_all_facets.pdf"), p_all, width = 9.5, height = 7.5)
