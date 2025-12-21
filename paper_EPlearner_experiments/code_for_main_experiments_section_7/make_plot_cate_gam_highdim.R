library(ggplot2)
library(data.table)
library(stringr)

# ----------------------------
# Reproducible base path + dirs
# ----------------------------
BASE_PATH <- Sys.getenv("BASE_PATH", unset = NA_character_)
if (is.na(BASE_PATH) || !nzchar(BASE_PATH)) {
  this_file <- tryCatch(normalizePath(sys.frame(1)$ofile, winslash = "/"), error = function(e) NA_character_)
  if (is.na(this_file)) stop("Set BASE_PATH env var, e.g. Sys.setenv(BASE_PATH='/path/to/project_root')")
  BASE_PATH <- normalizePath(dirname(this_file), winslash = "/")
} else {
  BASE_PATH <- normalizePath(path.expand(BASE_PATH), winslash = "/")
}

RESULTS_DIR   <- file.path(BASE_PATH, "experiment_results")              # <- your new convention
FIGURES_DIR   <- file.path(BASE_PATH, "mainSimResults")                  # keep your existing relative structure
PLOTS_DIR     <- file.path(FIGURES_DIR, "plots")

dir.create(RESULTS_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(FIGURES_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(PLOTS_DIR,   showWarnings = FALSE, recursive = TRUE)

# ----------------------------
# Inputs
# ----------------------------
ns <- sort(c(500, 1000, 2500, 5000))
hard_list <- c(TRUE, FALSE)
pos_list  <- c(TRUE, FALSE)
use_oracle_sieve <- FALSE

# Helper: load .RData safely from RESULTS_DIR
load_simresults <- function(n, hard, pos) {
  # Original: mainSimResults2/mainSimResults2/simsCATE{hard}{pos}n{n}_gam_highdim
  # Make it explicit and reproducible:
  rel_path <- file.path(
    "mainSimResults2", "mainSimResults2",
    paste0("simsCATE", hard, pos, "n", n, "_gam_highdim", ".RData")
  )
  f <- file.path(RESULTS_DIR, rel_path)

  if (!file.exists(f)) stop(sprintf("Missing file: %s", f))
  load(f)
  if (!exists("simresults", inherits = FALSE)) stop(sprintf("File loaded but 'simresults' not found: %s", f))
  simresults <- get("simresults", inherits = FALSE)
  simresults
}

for (pos in pos_list) {
  for (hard in hard_list) {
    try({

      sims_list <- lapply(ns, function(n) {
        try({

          simresults <- load_simresults(n = n, hard = hard, pos = pos)
          simresults <- simresults[sapply(simresults, is.list)]

          onestepbenchoracle <- rowMeans(do.call(cbind, lapply(simresults, `[[`, "CATEonestepbenchoracle")))
          onestepbench       <- rowMeans(do.call(cbind, lapply(simresults, `[[`, "CATEonestepbench")))

          substrisks   <- rowMeans(do.call(cbind, lapply(simresults, `[[`, "risk_subst")))
          cvsubstrisks <- rowMeans(do.call(cbind, lapply(simresults, `[[`, "risk_subst_cv")))

          lrnr_names <- names(simresults[[1]]$CATEonestepbench)
          lrnr_names <- unlist(lapply(lrnr_names, function(name) {
            paste0(name, c("_no_sieve.plugin", paste0("_fourier_basis_", 1:4, "_plugin")))
          }))
          iter <- rep(seq_along(simresults), each = length(lrnr_names))

          cvrisksDRoracle <- unlist(lapply(simresults, function(item) item$sieve$cvrisksDRoracle))
          cvrisksDR <- unlist(lapply(seq_along(simresults), function(index) {
            as.vector(simresults[[index]]$sieve$cvrisksDR)
          }))
          risks_oracle <- unlist(lapply(simresults, function(item) item$sieve$risks_oracle))

          dt <- data.table(iter, lrnr_full = lrnr_names, cvrisksDR, cvrisksDRoracle, risks_oracle)

          dt$degree <- as.numeric(stringr::str_match(dt$lrnr_full, "fourier_basis_([0-9]+)")[, 2])
          dt$degree[grep("no_sieve", dt$lrnr_full)] <- 0

          dt$lrnr[grep("gam3", dt$lrnr_full)] <- "gam3"
          dt$lrnr[grep("gam4", dt$lrnr_full)] <- "gam4"
          dt$lrnr[grep("gam5", dt$lrnr_full)] <- "gam5"
          dt$lrnr[grep("glm",  dt$lrnr_full)] <- "glm"
          dt$lrnr[grep("earth", dt$lrnr_full)] <- "earth"
          dt$lrnr[grep("rpart", dt$lrnr_full)] <- "rpart"
          dt$lrnr[grep("ranger_500_TRUE_none_1_7",  dt$lrnr_full)] <- "ranger_7"
          dt$lrnr[grep("ranger_500_TRUE_none_1_13", dt$lrnr_full)] <- "ranger_13"
          dt$lrnr[grep("ranger_500_TRUE_none_1_10", dt$lrnr_full)] <- "ranger_10"

          dt$lrnr <- gsub("Lrnr_", "", dt$lrnr)
          dt$lrnr <- gsub("_fourier_basis.+", "", dt$lrnr)
          dt$lrnr <- gsub(".fourier_basis.+", "", dt$lrnr)
          dt$lrnr <- gsub("_no_sieve.+", "", dt$lrnr)
          dt$lrnr <- gsub(".no_sieve.+", "", dt$lrnr)
          dt$lrnr <- gsub("_fourier.basis.+", "", dt$lrnr)
          dt$lrnr <- gsub(".fourier.basis.+", "", dt$lrnr)

          dt$type[!is.na(as.numeric(dt$degree))] <- "Sieve-Plugin"
          dt$type[ is.na(as.numeric(dt$degree))] <- dt$degree[is.na(as.numeric(dt$degree))]

          dt <- dt[dt$degree > 0]

          if (!use_oracle_sieve) {
            dt[!is.na(as.numeric(dt$degree)),
               risks_oracle := risks_oracle[which.min(cvrisksDR)],
               by = c("lrnr", "type", "iter")]
          } else {
            dt[!is.na(as.numeric(dt$degree)),
               risks_oracle := min(risks_oracle),
               by = c("iter", "lrnr", "type")]
          }
          dt[, risks_best := mean(risks_oracle), by = c("lrnr", "type")]

          # Add non-sieve baselines
          dt <- rbind(dt, data.table(risks_oracle = onestepbench,       lrnr_full = names(onestepbench), degree = "One-step"),       fill = TRUE)
          dt <- rbind(dt, data.table(risks_oracle = onestepbenchoracle, lrnr_full = names(onestepbench), degree = "Oracle one-step"), fill = TRUE)
          dt <- rbind(dt, data.table(risks_oracle = cvsubstrisks,       lrnr_full = names(onestepbench), degree = "Substitution-CV"), fill = TRUE)
          dt <- rbind(dt, data.table(risks_oracle = substrisks,         lrnr_full = names(onestepbench), degree = "Substitution"),    fill = TRUE)

          # Clean labels for baselines
          dt$lrnr <- dt$lrnr_full
          dt$lrnr[grep("gam3", dt$lrnr_full)] <- "gam3"
          dt$lrnr[grep("gam4", dt$lrnr_full)] <- "gam4"
          dt$lrnr[grep("gam5", dt$lrnr_full)] <- "gam5"
          dt$lrnr[grep("glm",  dt$lrnr_full)] <- "glm"
          dt$lrnr[grep("earth", dt$lrnr_full)] <- "earth"
          dt$lrnr[grep("rpart", dt$lrnr_full)] <- "rpart"
          dt$lrnr[grep("ranger_500_TRUE_none_1_7",  dt$lrnr_full)] <- "ranger_7"
          dt$lrnr[grep("ranger_500_TRUE_none_1_13", dt$lrnr_full)] <- "ranger_13"
          dt$lrnr[grep("ranger_500_TRUE_none_1_10", dt$lrnr_full)] <- "ranger_10"

          dt$lrnr <- gsub("Lrnr_", "", dt$lrnr)
          dt$lrnr <- gsub("_fourier_basis.+", "", dt$lrnr)
          dt$lrnr <- gsub(".fourier_basis.+", "", dt$lrnr)
          dt$lrnr <- gsub("_no_sieve.+", "", dt$lrnr)
          dt$lrnr <- gsub(".no_sieve.+", "", dt$lrnr)
          dt$lrnr <- gsub("_fourier.basis.+", "", dt$lrnr)
          dt$lrnr <- gsub(".fourier.basis.+", "", dt$lrnr)

          dt$type[!is.na(as.numeric(dt$degree))] <- "Sieve-Plugin"
          dt$type[ is.na(as.numeric(dt$degree))] <- dt$degree[is.na(as.numeric(dt$degree))]

          dt2 <- unique(dt[, .(lrnr, risks_best, type)])
          dt2$n <- n
          return(dt2)

        })
      })

      sims_list <- sims_list[sapply(sims_list, is.data.frame)]
      dt <- rbindlist(sims_list)

      dt$lrnr[grep("glm", dt$lrnr)] <- "GLM"
      dt$lrnr[grep("earth", dt$lrnr)] <- "MARS (earth)"
      s <- stringr::str_match(dt$lrnr[grep("gam", dt$lrnr)], "s(.+)_")[, 2]
      dt$lrnr[grep("gam", dt$lrnr)] <- paste0("GAM (s=", s, ")")
      dt$lrnr[grep("GAM", dt$lrnr)][is.na(s)] <- "GAM (s=cv)"

      dt_tmp <- copy(dt)
      dt_tmp <- dt_tmp[dt_tmp$type != "Substitution"]
      dt_tmp[type == "Sieve-Plugin",     type := "EP-Learner (*)"]
      dt_tmp[type == "One-step",         type := "DR-Learner"]
      dt_tmp[type == "Oracle one-step",  type := "Oracle DR-Learner"]
      dt_tmp[type == "Substitution-CV",  type := "T-Learner (CV)"]

      # ---- Plot 1 (saved under BASE_PATH/mainSimResults/...) ----
      plt <- ggplot(dt_tmp, aes(x = n, y = risks_best, group = type, color = type, linetype = type)) +
        geom_line(size = 0.5) +
        facet_wrap(~lrnr, scales = "free") +
        theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 1)) +
        scale_y_log10(limits = c(min(1e-1, min(dt_tmp$risks_best)), max(dt_tmp$risks_best))) +
        scale_x_log10(breaks = c(500, 1000, 2500, 5000, 10000)) +
        xlab("Sample Size (n)") + ylab("Mean-Squared-Error (MSE)") +
        theme_bw() +
        labs(color = "Method", group = "Method", linetype = "Method") +
        theme(axis.text = element_text(size = 12),
              legend.text = element_text(size = 12),
              legend.title = element_text(size = 12),
              axis.title = element_text(size = 14, face = "bold"),
              legend.justification = c(0.9, 0),
              legend.position = c(0.9, 0))

      ggsave(
        filename = file.path(FIGURES_DIR, paste0("performancePlot_CATE_GAM_highDim_pos=", pos, "hard=", hard, ".pdf")),
        plot = plt, width = 8, height = 7
      )

      # ---- Plot 2 (main paper subset) ----
      dt_tmp2 <- copy(dt)
      dt_tmp2$lrnr[grep("glm", dt_tmp2$lrnr)] <- "GLM"
      dt_tmp2$lrnr[grep("earth", dt_tmp2$lrnr)] <- "MARS (earth)"
      s <- stringr::str_match(dt_tmp2$lrnr[grep("gam", dt_tmp2$lrnr)], "s(.+)_")[, 2]
      dt_tmp2$lrnr[grep("gam", dt_tmp2$lrnr)] <- paste0("GAM (s=", s, ")")
      dt_tmp2$lrnr[grep("GAM", dt_tmp2$lrnr)][is.na(s)] <- "GAM (s=cv)"

      keep <- dt_tmp2$lrnr %in% c("MARS (earth)", "GAM (s=1)", "GAM (s=3)", "GAM (s=cv)")
      dt_tmp2 <- dt_tmp2[keep, ]
      dt_tmp2 <- dt_tmp2[dt_tmp2$type != "Substitution"]
      dt_tmp2[type == "Sieve-Plugin",     type := "EP-Learner (*)"]
      dt_tmp2[type == "One-step",         type := "DR-Learner"]
      dt_tmp2[type == "Oracle one-step",  type := "Oracle DR-Learner"]
      dt_tmp2[type == "Substitution-CV",  type := "T-Learner (CV)"]

      plt2 <- ggplot(dt_tmp2, aes(x = n, y = risks_best, group = type, color = type, linetype = type)) +
        geom_line(size = 0.5) +
        facet_wrap(~lrnr, scales = "free") +
        theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 1)) +
        scale_y_log10(limits = c(min(1e-1, min(dt_tmp2$risks_best)), max(dt_tmp2$risks_best))) +
        scale_x_log10(breaks = c(500, 1000, 2500, 5000, 10000)) +
        xlab("Sample Size (n)") + ylab("Mean-Squared-Error (MSE)") +
        theme_bw() +
        labs(color = "Method", group = "Method", linetype = "Method") +
        theme(axis.text = element_text(size = 12),
              legend.text = element_text(size = 12),
              legend.title = element_text(size = 12),
              axis.title = element_text(size = 14, face = "bold"),
              legend.justification = c(0.9, 0),
              legend.position = c(0.95, 0.23))

      ggsave(
        filename = file.path(FIGURES_DIR, paste0("performancePlot_CATE_GAM_highDim_pos=", pos, "hard=", hard, "_mainpaper.pdf")),
        plot = plt2, width = 8, height = 7
      )

      # ---- Plot 3: per-learner panels saved to mainSimResults/plots ----
      dt_tmp3 <- as.data.frame(dt_tmp2)
      dt_tmp3 <- dt_tmp3[dt_tmp3$type != "Oracle DR-Learner", ]

      for (lrnr in unique(dt_tmp3$lrnr)) {
        plt3 <- ggplot(dt_tmp3[dt_tmp3$lrnr %in% lrnr, ],
                       aes(x = n, y = risks_best, group = type, color = type, linetype = type)) +
          geom_line(size = 0.75) +
          facet_wrap(~lrnr, scales = "free") +
          theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 1)) +
          scale_y_log10(limits = c(min(1e-1, min(dt_tmp3$risks_best)), max(dt_tmp3$risks_best))) +
          scale_x_log10(breaks = c(500, 1000, 2500, 5000, 10000)) +
          xlab("Sample Size (n)") + ylab("Mean-Squared-Error (MSE)") +
          theme_bw() +
          labs(color = "Method", group = "Method", linetype = "Method") +
          theme(axis.text = element_text(size = 14),
                strip.text.x = element_text(size = 20),
                legend.text = element_text(size = 12),
                legend.title = element_text(size = 12),
                axis.title = element_text(size = 12, face = "bold"),
                legend.position = "none")

        ggsave(
          filename = file.path(PLOTS_DIR, paste0("performancePlot_CATE_GAM_highDim_pos=", pos, "hard=", hard, "_", lrnr, ".pdf")),
          plot = plt3, width = 4, height = 4
        )
      }

    })
  }
}
