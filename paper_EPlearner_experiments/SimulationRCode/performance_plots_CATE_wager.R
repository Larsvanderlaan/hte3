# make_plots_wager.R
# Plots performance of EP / DR / R / T learners from CSVs written by EPlearnerCATE_wager.R

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(stringr)
})

# ---------- Paths ----------
results_dir <- "~/repos/hte3/paper_EPlearner_experiments/results"
plots_dir   <- file.path(results_dir, "plots_wager")
if (!dir.exists(plots_dir)) dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)

# ---------- Helper: read & combine all CSVs ----------
# Expects files named like:
# simsEP_n=2000_d=5_sigma=1_base_learner=xg_sim_type=A_nsims=100.csv
all_csvs <- list.files(results_dir, pattern = "^simsEP_.*\\.csv$", full.names = TRUE)
if (length(all_csvs) == 0L) {
  stop("No result CSVs found in: ", results_dir,
       "\nDid you run the sims and write the CSV outputs?")
}

read_one <- function(f) {
  dt <- tryCatch(fread(f), error = function(e) NULL)
  if (is.null(dt) || nrow(dt) == 0L) return(NULL)
  # add file-derived fields as fallback if missing
  info <- str_match(basename(f),
                    "n=([0-9]+)_d=([0-9]+)_sigma=([0-9\\.]+)_base_learner=([a-z]+)_sim_type=([A-D])_nsims=([0-9]+)")
  if (!"n"     %in% names(dt)) dt[, n := as.numeric(info[2])]
  if (!"d"     %in% names(dt)) dt[, d := as.numeric(info[3])]
  if (!"sigma" %in% names(dt)) dt[, sigma := as.numeric(info[4])]
  if (!"lrnr"  %in% names(dt)) dt[, lrnr := info[6]]  # base_learner
  if (!"sim_type" %in% names(dt)) dt[, sim_type := info[7]]
  dt[]
}

res_list <- lapply(all_csvs, read_one)
res_list <- res_list[!vapply(res_list, is.null, logical(1))]
results  <- rbindlist(res_list, use.names = TRUE, fill = TRUE)

if (nrow(results) == 0L) stop("Parsed 0 rows from results CSVs. Check file names/format.")

# ---------- Tidy to long format (EP/DR/R/T) ----------

# ---------- Remove iterations where DR-Learner explodes ----------
# Tunables:
abs_cut <- 1e3     # absolute MSE cutoff for DR (e.g., > 1e3)
rel_cut <- 50      # relative cutoff vs. within-cell median DR (e.g., > 50× median)

# Define the "cell" = (sim_type, base_learner, n)
by_keys <- c("sim_type", "lrnr", "n")

# Median DR per cell (robust baseline)
dr_ref <- results[
  , .(dr_median = median(dr_risk, na.rm = TRUE)),
  by = by_keys
]

res2 <- merge(results, dr_ref, by = by_keys, all.x = TRUE)

# Flag bad DR rows by absolute OR relative rule
res2[, dr_bad := (!is.na(dr_risk)) & (
  dr_risk > abs_cut |
    (is.finite(dr_median) & dr_median > 0 & dr_risk > rel_cut * dr_median)
)]

# List which (iter) got dropped per cell
bad_iters <- unique(res2[dr_bad == TRUE, .(sim_type, lrnr, n, iter, dr_risk)])
if (nrow(bad_iters)) {
  if (!dir.exists(plots_dir)) dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)
  fwrite(bad_iters, file.path(plots_dir, "dropped_iters_dr_excessive.csv"))
  message(sprintf("Dropping %d iters where DR risk exceeded thresholds. Logged to dropped_iters_dr_excessive.csv",
                  nrow(bad_iters)))
}

# Remove those iterations from ALL methods in that cell
results <- res2[dr_bad == FALSE, !(c("dr_median","dr_bad")), with = FALSE]



# ---------- Tidy to long format (EP/DR/R/T + CF when present) ----------
meas <- intersect(
  c("ep_risk", "dr_risk", "r_risk", "t_risk", "causal_forest_risk"),
  names(results)
)

long <- melt(
  results,
  id.vars = c("iter", "sim_type", "n", "d", "sigma", "lrnr"),
  measure.vars = meas,
  variable.name = "method",
  value.name   = "mse"
)

# Map method codes to nice labels (ONCE)
method_map <- c(
  ep_risk            = "EP-Learner (ours)",
  dr_risk            = "DR-Learner",
  r_risk             = "R-Learner",
  t_risk             = "T-Learner",
  causal_forest_risk = "Causal Forest"
)
long[, method := factor(method_map[method],
                        levels = c("EP-Learner (ours)","DR-Learner","R-Learner","T-Learner","Causal Forest"))]

# Keep CF only for rf; drop NA metrics/labels and unused levels
long <- long[!(method == "Causal Forest" & lrnr != "rf") & !is.na(mse) & !is.na(method)]
long[, method := droplevels(method)]

# Map base_learner codes to nicer facet labels
blr_map <- c(xg = "XGBoost", rf = "Random Forest", gam = "GAM")
long[, base_learner := factor(blr_map[lrnr], levels = c("XGBoost","Random Forest","GAM"))]

# Ensure types
long[, `:=`(n = as.numeric(n), sigma = as.numeric(sigma))]
long[, sim_type := factor(sim_type, levels = c("A","B","C","D"))]

# ---------- Aggregate across iterations ----------
# Compute mean and standard error across iter for each cell
summ <- long[
  , .(
    mean_mse = mean(mse, na.rm = TRUE),
    se_mse   = sd(mse,  na.rm = TRUE) / sqrt(sum(!is.na(mse))),
    n_iter   = sum(!is.na(mse))
  ),
  by = .(sim_type, base_learner, n, method)
][n_iter > 0]   # <— drop all-NA cells (e.g., CF when not rf)


# Optional: 95% CI (on linear scale; we’ll plot on log y)
summ[, `:=`(
  lo = pmax(mean_mse - 1.96 * se_mse, .Machine$double.eps),
  hi = mean_mse + 1.96 * se_mse
)]
summ <- summ[!is.na(summ$base_learner),]



# ---------- Tables mirroring the plot ----------
# Order rows like the facets (A–D by rows; XGBoost/RF/GAM by cols), and then n
setkeyv(summ, c("sim_type", "base_learner", "n", "method"))

# 1) Long (tidy) table
tbl_long <- copy(summ)[
  , .(sim_type, base_learner, n, method, mean_mse, se_mse, n_iter)
]
fwrite(tbl_long, file.path(plots_dir, "performance_table_long.csv"))

# 2) Wide table (means only)
tbl_wide <- dcast(
  tbl_long,
  sim_type + base_learner + n ~ method,
  value.var = "mean_mse"
)
# Optional ordering of columns
setcolorder(tbl_wide, c("sim_type","base_learner","n",
                        "EP-Learner (ours)","DR-Learner","R-Learner","T-Learner"))
fwrite(tbl_wide, file.path(plots_dir, "performance_table_wide.csv"))

# 3) Wide, pretty cells: "mean (SE) [n]"
fmt_num <- function(x) formatC(x, format = "e", digits = 2)  # scientific, 2 sig figs
tbl_pretty_src <- copy(summ)[
  , .(
    cell = sprintf("%s (%s) [%d]", fmt_num(mean_mse), fmt_num(se_mse), n_iter)
  ),
  by = .(sim_type, base_learner, n, method)
]
tbl_wide_pretty <- dcast(
  tbl_pretty_src,
  sim_type + base_learner + n ~ method,
  value.var = "cell"
)
setcolorder(tbl_wide_pretty, c("sim_type","base_learner","n",
                               "EP-Learner (ours)","DR-Learner","R-Learner","T-Learner"))
fwrite(tbl_wide_pretty, file.path(plots_dir, "performance_table_wide_pretty.csv"))



# ---------- Plotting helpers ----------
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
    labs(x = "Sample Size (n)", y = "Mean Squared Error (MSE)", color = "Method",
         linetype = "Method", fill = "Method", title = title) +
    base_theme +
    theme(
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.box = "horizontal"
    )
}


# ---------- 1) Full grid: sim_type (rows) × base_learner (cols) ----------
p_all <- make_plot(summ) + facet_grid(sim_type ~ base_learner, scales = "free_y") + theme(
  panel.spacing.x = unit(1.2, "cm")  # increase horizontal spacing between facets
)

ggsave(file.path(plots_dir, "performance_all_facets.pdf"), p_all, width = 9.5, height = 7.5)

# ---------- 2) Per base learner: facet by sim_type ----------
for (bl in levels(summ$base_learner)) {
  dt_bl <- summ[base_learner == bl]
  if (nrow(dt_bl) == 0L) next
  p_bl <- make_plot(dt_bl, title = bl) + facet_wrap(~ sim_type, nrow = 1, scales = "free_y")
  fn <- sprintf("performance_by_simtype_%s.pdf", gsub(" ", "_", tolower(bl)))
  ggsave(file.path(plots_dir, fn), p_bl, width = 10, height = 3.4)
}

# ---------- 3) Per sim_type: facet by base_learner ----------
for (st in levels(summ$sim_type)) {
  dt_st <- summ[sim_type == st]
  if (nrow(dt_st) == 0L) next
  p_st <- make_plot(dt_st, title = paste("Setup", st)) + facet_wrap(~ base_learner, nrow = 1, scales = "free_y")
  fn <- sprintf("performance_by_baselearner_setup_%s.pdf", st)
  ggsave(file.path(plots_dir, fn), p_st, width = 10, height = 3.4)
}

# ---------- 4) Small multiplies: one PDF per base_learner × sim_type ----------
for (bl in levels(summ$base_learner)) {
  for (st in levels(summ$sim_type)) {
    dt_bls <- summ[base_learner == bl & sim_type == st]
    if (nrow(dt_bls) == 0L) next
    p_bls <- make_plot(dt_bls, title = paste(bl, "- Setup", st))
    fn <- sprintf("performance_%s_setup_%s.pdf",
                  gsub(" ", "_", tolower(bl)), st)
    ggsave(file.path(plots_dir, fn), p_bls, width = 4.5, height = 3.6)
  }
}

# ---------- 5) Table of winners (best mean MSE) ----------
winners <- summ[, .SD[which.min(mean_mse)], by = .(sim_type, base_learner)]
fwrite(winners, file.path(plots_dir, "winners_by_cell.csv"))

message("Done. Plots written to: ", plots_dir)
