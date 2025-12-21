# ------------------------------------------------------------
# Paths / reproducibility (single source of truth)
# ------------------------------------------------------------
BASE_PATH <- Sys.getenv("BASE_PATH", unset = NA_character_)
if (is.na(BASE_PATH) || !nzchar(BASE_PATH)) {
  this_file <- tryCatch(normalizePath(sys.frame(1)$ofile, winslash = "/"),
                        error = function(e) NA_character_)
  if (is.na(this_file)) {
    stop("Set BASE_PATH env var, e.g. Sys.setenv(BASE_PATH='/path/to/project_root')")
  }
  BASE_PATH <- normalizePath(dirname(this_file), winslash = "/")
} else {
  BASE_PATH <- normalizePath(path.expand(BASE_PATH), winslash = "/")
}


RESULTS_DIR <- file.path(BASE_PATH, "experiment_results")  # <- new convention
FIGURES_DIR <- file.path(BASE_PATH, "experiment_results")      # keep existing relative structure
PLOTS_DIR   <- file.path(FIGURES_DIR, "plots")

dir.create(PLOTS_DIR,   recursive = TRUE, showWarnings = FALSE)

# Canonical helpers (so paths are never scattered)
sim_file <- function(hard, pos, n) {
  file.path(RESULTS_DIR, "mainSimResults2", "mainSimResults2",
            paste0("simsCATE", hard, pos, "n", n, "_xgboost_highDim"))
}

out_file <- function(stem, pos, hard) {
  file.path(FIGURES_DIR, paste0(stem, "pos=", pos, "hard=", hard, ".pdf"))
}

out_plot_file <- function(stem, pos, hard, lrnr) {
  file.path(PLOTS_DIR, paste0(stem, "pos=", pos, "hard=", hard, "_", lrnr, ".pdf"))
}

# ------------------------------------------------------------
# Your script continues below (unchanged except for load/ggsave paths)
# ------------------------------------------------------------
library(ggplot2)
library(data.table)

ns <- sort(c(500, 1000, 2500, 5000))
hard_list <- c(TRUE, FALSE)
pos_list  <- c(TRUE, FALSE)
use_oracle_sieve <- FALSE

for (pos in pos_list) {
  for (hard in hard_list) {
    try({

      sims_list <- lapply(ns, function(n) {
        try({
          load(sim_file(hard, pos, n))
          simresults <- get("simresults")
          simresults <- simresults[sapply(simresults, is.list)]

          onestepbenchoracle <- rowMeans(do.call(cbind, lapply(simresults, `[[`, "CATEonestepbenchoracle")))
          onestepbench       <- rowMeans(do.call(cbind, lapply(simresults, `[[`, "CATEonestepbench")))
          causalforestrisks  <- rowMeans(do.call(cbind, lapply(simresults, `[[`, "risk_cf")))
          substrisks         <- rowMeans(do.call(cbind, lapply(simresults, `[[`, "risk_subst")))
          cvsubstrisks       <- rowMeans(do.call(cbind, lapply(simresults, `[[`, "risk_subst_cv")))

          lrnr_names <- names(simresults[[1]]$CATEonestepbench)
          lrnr_names <- unlist(lapply(lrnr_names, function(name) {
            paste0(name, c("_no_sieve.plugin", paste0("_fourier_basis_", 1:4, "_plugin")))
          }))

          iter <- rep(1:length(simresults), each = length(lrnr_names))

          cvrisksDRoracle <- unlist(lapply(simresults, function(item) item$sieve$cvrisksDRoracle))
          cvrisksDR <- unlist(lapply(seq_along(simresults), function(index) {
            item <- simresults[[index]]
            as.vector(item$sieve$cvrisksDR)
          }))
          risks_oracle <- unlist(lapply(simresults, function(item) item$sieve$risks_oracle))

          dt <- data.table(iter, lrnr_full = lrnr_names, cvrisksDR, cvrisksDRoracle, risks_oracle)
          dt$degree <- as.numeric(stringr::str_match(dt$lrnr_full, "fourier_basis_([0-9]+)")[, 2])
          dt$degree[grep("no_sieve", dt$lrnr_full)] <- 0

          dt$lrnr[grep("gam3", dt$lrnr_full)] <- "gam3"
          dt$lrnr[grep("gam4", dt$lrnr_full)] <- "gam4"
          dt$lrnr[grep("gam5", dt$lrnr_full)] <- "gam5"
          dt$lrnr[grep("glm",  dt$lrnr_full)] <- "glm"
          dt$lrnr[grep("earth",dt$lrnr_full)] <- "earth"
          dt$lrnr[grep("rpart",dt$lrnr_full)] <- "rpart"
          dt$lrnr[grep("ranger_500_TRUE_none_1_7",  dt$lrnr_full)] <- "ranger_7"
          dt$lrnr[grep("ranger_500_TRUE_none_1_13", dt$lrnr_full)] <- "ranger_13"
          dt$lrnr[grep("ranger_500_TRUE_none_1_10", dt$lrnr_full)] <- "ranger_10"

          dt$lrnr <- gsub("Lrnr_", "", dt$lrnr)
          dt$lrnr <- gsub("_fourier_basis.+", "", dt$lrnr)

          dt$type[!is.na(as.numeric(dt$degree))] <- "Sieve-Plugin"
          dt$type[ is.na(as.numeric(dt$degree))] <- dt$degree[is.na(as.numeric(dt$degree))]

          dt <- dt[dt$degree > 0]

          tmp <- data.table(risks_oracle = onestepbench, lrnr_full = names(onestepbench), degree = "One-step")
          dt <- rbind(dt, tmp, fill = TRUE)
          tmp <- data.table(risks_oracle = onestepbenchoracle, lrnr_full = names(onestepbench), degree = "Oracle one-step")
          dt <- rbind(dt, tmp, fill = TRUE)
          tmp <- data.table(risks_oracle = cvsubstrisks, lrnr_full = names(onestepbench), degree = "Substitution-CV")
          dt <- rbind(dt, tmp, fill = TRUE)
          tmp <- data.table(risks_oracle = causalforestrisks, lrnr_full = names(onestepbench),
                            lrnr = "Causal-Forest", degree = "Causal-Forest")
          dt <- rbind(dt, tmp, fill = TRUE)
          tmp <- data.table(risks_oracle = substrisks, lrnr_full = names(onestepbench), degree = "Substitution")
          dt <- rbind(dt, tmp, fill = TRUE)

          dt$lrnr <- dt$lrnr_full
          dt$lrnr[grep("gam3", dt$lrnr_full)] <- "gam3"
          dt$lrnr[grep("gam4", dt$lrnr_full)] <- "gam4"
          dt$lrnr[grep("gam5", dt$lrnr_full)] <- "gam5"
          dt$lrnr[grep("glm",  dt$lrnr_full)] <- "glm"
          dt$lrnr[grep("earth",dt$lrnr_full)] <- "earth"
          dt$lrnr[grep("rpart",dt$lrnr_full)] <- "rpart"
          dt$lrnr[grep("ranger_500_TRUE_none_1_7",  dt$lrnr_full)] <- "ranger_7"
          dt$lrnr[grep("ranger_500_TRUE_none_1_13", dt$lrnr_full)] <- "ranger_13"
          dt$lrnr[grep("ranger_500_TRUE_none_1_10", dt$lrnr_full)] <- "ranger_10"

          dt$lrnr <- gsub("Lrnr_", "", dt$lrnr)
          dt$lrnr <- gsub("_fourier_basis.+", "", dt$lrnr)

          dt$type[!is.na(as.numeric(dt$degree))] <- "Sieve-Plugin"
          dt$type[ is.na(as.numeric(dt$degree))] <- dt$degree[is.na(as.numeric(dt$degree))]

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
          dt2 <- unique(dt[, .(lrnr, risks_best, type)])
          dt2$n <- n
          dt2
        })
      })

      sims_list <- sims_list[sapply(sims_list, is.data.frame)]
      dt <- rbindlist(sims_list)

      # --- XGBoost facet plot
      dt_tmp <- dt[grep("xgboost", dt$lrnr), ]
      max_depth <- stringr::str_match(dt_tmp$lrnr, "([0-9]+)_0$")[, 2]
      max_depth[is.na(max_depth)] <- "cv"
      dt_tmp$lrnr <- paste0("xgboost (", "max_depth=", max_depth, ")")
      dt_tmp <- dt_tmp[max_depth %in% c("1", "2", "3", "5", "cv"), ]
      dt_tmp <- rbind(dt_tmp, dt[grep("Causal-Forest", dt$lrnr), ])
      dt_tmp <- dt_tmp[!(dt_tmp$type == "Substitution"), ]
      dt_tmp[(dt_tmp$type == "Sieve-Plugin"), "type"] <- "EP-Learner (*)"
      dt_tmp[(dt_tmp$type == "One-step"), "type"] <- "DR-Learner"
      dt_tmp[(dt_tmp$type == "Oracle one-step"), "type"] <- "Oracle DR-Learner"
      dt_tmp[(dt_tmp$type == "Substitution-CV"), "type"] <- "T-Learner (CV)"

      plt <- ggplot(dt_tmp, aes(x = n, y = risks_best, group = type, color = type, linetype = type)) +
        geom_line() +
        facet_wrap(~lrnr, scales = "free") +
        theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 1)) +
        scale_y_log10(limits = range(dt_tmp$risks_best)) +
        scale_x_log10(breaks = c(500, 1000, 2500, 5000, 10000)) +
        xlab("Sample Size (n)") + ylab("Mean-Squared-Error (MSE)") +
        theme_bw() + labs(color = "Method", group = "Method", linetype = "Method") +
        theme(axis.text = element_text(size = 12),
              legend.text = element_text(size = 12),
              legend.title = element_text(size = 12),
              axis.title = element_text(size = 14, face = "bold")) +
        theme(legend.justification = c(0.9, 0), legend.position = c(0.9, 0))

      ggsave(out_file("performancePlot_CATE_highDim_xgboost", pos, hard), plot = plt, width = 8, height = 7)

      # --- Ranger facet plot
      dt_tmp <- dt[!(dt$lrnr %in% c("glm", "earth", "gam3", "gam4", "gam5")), ]
      dt_tmp <- dt_tmp[grep("rf", dt$lrnr), ]
      max_depth <- stringr::str_match(dt_tmp$lrnr, "([0-9]+)_xg$")[, 2]
      max_depth[is.na(max_depth)] <- "cv"
      dt_tmp$lrnr <- paste0("ranger (", "max_depth=", max_depth, ")")
      dt_tmp <- dt_tmp[max_depth %in% c("3", "5", "7", "9", "cv"), ]
      dt_tmp$lrnr <- factor(dt_tmp$lrnr, levels = paste0("ranger (", "max_depth=", c("3", "5", "7", "9", "cv"), ")"))
      dt_tmp <- rbind(dt_tmp, dt[grep("causalforest", dt$lrnr), ])
      dt_tmp <- dt_tmp[!(dt_tmp$type == "Substitution"), ]
      dt_tmp[(dt_tmp$type == "Sieve-Plugin"), "type"] <- "EP-Learner (*)"
      dt_tmp[(dt_tmp$type == "One-step"), "type"] <- "DR-Learner"
      dt_tmp[(dt_tmp$type == "Oracle one-step"), "type"] <- "Oracle DR-Learner"
      dt_tmp[(dt_tmp$type == "Substitution-CV"), "type"] <- "T-Learner (CV)"

      plt <- ggplot(dt_tmp, aes(x = n, y = risks_best, group = type, color = type, linetype = type)) +
        geom_line() +
        facet_wrap(~lrnr, scales = "free") +
        theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 1)) +
        scale_y_log10(limits = range(dt_tmp$risks_best)) +
        scale_x_log10(breaks = c(500, 1000, 2500, 5000, 10000)) +
        xlab("Sample Size (n)") + ylab("Mean-Squared-Error (MSE)") +
        theme_bw() + labs(color = "Method", group = "Method", linetype = "Method") +
        theme(axis.text = element_text(size = 12),
              legend.text = element_text(size = 12),
              legend.title = element_text(size = 12),
              axis.title = element_text(size = 14, face = "bold")) +
        theme(legend.justification = c(0.9, 0), legend.position = c(0.9, 0))

      ggsave(out_file("performancePlot_CATE_highDim_ranger", pos, hard), plot = plt, width = 8, height = 7)

      # --- Main-paper combined plot (keep your original naming)
      dt_tmp <- dt[!(dt$lrnr %in% c("glm", "earth", "gam3", "gam4", "gam5")), ]
      dt_tmp <- dt_tmp[grep("rf", dt$lrnr), ]
      max_depth <- stringr::str_match(dt_tmp$lrnr, "([0-9]+)_xg$")[, 2]
      max_depth[is.na(max_depth)] <- "cv"
      dt_tmp$lrnr <- paste0("ranger (", "max_depth=", max_depth, ")")
      dt_tmp <- dt_tmp[max_depth %in% c("3", "7", "cv"), ]
      dt_tmp$lrnr <- factor(dt_tmp$lrnr, levels = paste0("ranger (", "max_depth=", c("3", "5", "7", "9", "cv"), ")"))
      dt_tmp <- rbind(dt_tmp, dt[grep("causalforest", dt$lrnr), ])
      dt_tmp <- dt_tmp[!(dt_tmp$type == "Substitution"), ]
      dt_tmp[(dt_tmp$type == "Sieve-Plugin"), "type"] <- "EP-Learner (*)"
      dt_tmp[(dt_tmp$type == "One-step"), "type"] <- "DR-Learner"
      dt_tmp[(dt_tmp$type == "Oracle one-step"), "type"] <- "Oracle DR-Learner"
      dt_tmp[(dt_tmp$type == "Substitution-CV"), "type"] <- "T-Learner (CV)"
      dt_rf <- dt_tmp

      dt_tmp <- dt[grep("xgboost", dt$lrnr), ]
      max_depth <- stringr::str_match(dt_tmp$lrnr, "([0-9]+)_0$")[, 2]
      max_depth[is.na(max_depth)] <- "cv"
      dt_tmp$lrnr <- paste0("xgboost (", "max_depth=", max_depth, ")")
      dt_tmp <- dt_tmp[max_depth %in% c("1", "5", "cv"), ]
      dt_tmp <- rbind(dt_tmp, dt[grep("Causal-Forest", dt$lrnr), ])
      dt_tmp <- dt_tmp[!(dt_tmp$type == "Substitution"), ]
      dt_tmp[(dt_tmp$type == "Sieve-Plugin"), "type"] <- "EP-Learner (*)"
      dt_tmp[(dt_tmp$type == "One-step"), "type"] <- "DR-Learner"
      dt_tmp[(dt_tmp$type == "Oracle one-step"), "type"] <- "Oracle DR-Learner"
      dt_tmp[(dt_tmp$type == "Substitution-CV"), "type"] <- "T-Learner (CV)"
      dt_xg <- dt_tmp

      dt_comb <- rbind(dt_xg, dt_rf)

      plt <- ggplot(dt_comb, aes(x = n, y = risks_best, group = type, color = type, linetype = type)) +
        geom_line() +
        facet_wrap(~lrnr, scales = "free") +
        theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 1)) +
        scale_y_log10(limits = range(dt_comb$risks_best)) +
        scale_x_log10(breaks = c(500, 1000, 2500, 5000, 10000)) +
        xlab("Sample Size (n)") + ylab("Mean-Squared-Error (MSE)") +
        theme_bw() + labs(color = "Method", group = "Method", linetype = "Method") +
        theme(axis.text = element_text(size = 12),
              legend.text = element_text(size = 12),
              legend.title = element_text(size = 12),
              axis.title = element_text(size = 14, face = "bold")) +
        theme(legend.justification = c(0.9, 0), legend.position = c(0.95, 0.2))

      ggsave(
        file.path(FIGURES_DIR, paste0("performancePlot_CATE_highDim_tree", "pos=", pos, "hard=", hard, "mainpaper.pdf")),
        plot = plt, width = 8, height = 7
      )

      # --- Per-learner small plots
      for (lrnr in unique(dt$lrnr)) {
        dt_df <- as.data.frame(dt)
        dt_df <- dt_df[dt_df$type != "Oracle DR-Learner", ]

        plt <- ggplot(dt_df[dt_df$lrnr %in% lrnr, ],
                      aes(x = n, y = risks_best, group = type, color = type, linetype = type)) +
          geom_line(size = 0.75) +
          theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 1)) +
          scale_y_log10(limits = c(min(1e-1, min(dt_df$risks_best)), max(dt_df$risks_best))) +
          scale_x_log10(breaks = c(500, 1000, 2500, 5000, 10000)) +
          facet_wrap(~lrnr, scales = "free") +
          xlab("Sample Size (n)") + ylab("Mean-Squared-Error (MSE)") +
          theme_bw() + labs(color = "Method", group = "Method", linetype = "Method") +
          theme(axis.text = element_text(size = 14),
                strip.text.x = element_text(size = 20),
                legend.text = element_text(size = 12),
                legend.title = element_text(size = 12),
                axis.title = element_text(size = 12, face = "bold"))

        labels   <- c("Causal-Forest", "DR-Learner", "EP-Learner (*)", "T-Learner (CV)")
        colors   <- c("#619CFF", "#00BA38", "#F8766D", "#E76BF3")
        linetype <- c("longdash", "dashed", "solid", "dotted")
        names(colors) <- labels
        names(linetype) <- labels

        plt <- plt +
          scale_colour_manual(values = colors) +
          scale_linetype_manual(values = linetype) +
          theme(legend.position = "none")

        ggsave(out_plot_file("performancePlot_CATE_highDim_tree", pos, hard, lrnr),
               plot = plt, width = 4, height = 4)
      }

    })
  }
}
