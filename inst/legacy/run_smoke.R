legacy_ref <- Sys.getenv("HTE3_LEGACY_REF", unset = "legacy-paper-repro")
repo_root <- normalizePath(".", winslash = "/")
paper_root <- file.path(repo_root, "paper_EPlearner_experiments")
temp_lib <- file.path(tempdir(), "hte3-legacy-lib")
dir.create(temp_lib, recursive = TRUE, showWarnings = FALSE)

old_libpaths <- .libPaths()
on.exit(.libPaths(old_libpaths), add = TRUE)
.libPaths(c(temp_lib, old_libpaths))

if (!requireNamespace("remotes", quietly = TRUE)) {
  stop("The 'remotes' package is required to run the legacy smoke test.")
}

remotes::install_gitlocal(
  repo = repo_root,
  ref = legacy_ref,
  upgrade = "never",
  dependencies = TRUE
)

run_legacy_script <- function(script_path, bindings) {
  env <- list2env(bindings, parent = baseenv())
  sys.source(script_path, envir = env)
}

assert_metric_file <- function(path, required_columns) {
  if (!file.exists(path)) {
    stop(sprintf("Expected smoke-test output was not created: %s", path))
  }

  results <- data.table::fread(path)
  missing_columns <- setdiff(required_columns, names(results))
  if (length(missing_columns) > 0L) {
    stop(sprintf("Smoke-test output is missing columns: %s", paste(missing_columns, collapse = ", ")))
  }

  metric_columns <- intersect(required_columns, names(results))
  metric_values <- unlist(results[, ..metric_columns], use.names = FALSE)
  metric_values <- metric_values[is.finite(metric_values)]

  if (length(metric_values) == 0L) {
    stop(sprintf("Smoke-test output contained no finite metric values: %s", path))
  }

  if (any(metric_values < 0 | metric_values > 100, na.rm = TRUE)) {
    stop(sprintf("Smoke-test metrics fell outside the expected range in %s", path))
  }
}

Sys.setenv(BASE_PATH = file.path(paper_root, "code_for_main_experiments_section_7"))

run_legacy_script(
  file.path(paper_root, "code_for_main_experiments_section_7", "experiment_CATE.R"),
  list(
    n = 50,
    nsims = 2,
    hard = FALSE,
    pos = TRUE,
    sim_type = "CATElow",
    base_learner = "gam"
  )
)

assert_metric_file(
  file.path(paper_root, "experiment_results", "simsEP_n=50_pos=TRUE_hard=FALSE_base_learner=gam_sim_type=CATElow_nsims=2.csv"),
  c("ep_risk", "dr_risk", "r_risk", "t_risk")
)

run_legacy_script(
  file.path(paper_root, "code_for_main_experiments_section_7", "experiment_logCRR.R"),
  list(
    n = 50,
    nsims = 2,
    hard = FALSE,
    pos = TRUE,
    base_learner = "gam"
  )
)

assert_metric_file(
  file.path(paper_root, "experiment_results", "simsEP_n=50_pos=TRUE_hard=FALSE_base_learner=gam_sim_type=crr_nsims=2.csv"),
  c("ep_risk", "ipw_risk", "t_risk")
)
