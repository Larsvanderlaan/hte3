# ------------------------------------------------------------
# Paths / reproducibility (single source of truth)
# ------------------------------------------------------------
BASE_PATH <- Sys.getenv("BASE_PATH", unset = NA_character_)
if (is.na(BASE_PATH) || !nzchar(BASE_PATH)) {
  this_file <- tryCatch(normalizePath(sys.frame(1)$ofile, winslash = "/"),
                        error = function(e) NA_character_)
  if (is.na(this_file)) stop("Set BASE_PATH env var, e.g. Sys.setenv(BASE_PATH='/path/to/project_root')")
  BASE_PATH <- normalizePath(dirname(this_file), winslash = "/")
} else {
  BASE_PATH <- normalizePath(path.expand(BASE_PATH), winslash = "/")
}

RESULTS_DIR <- file.path(BASE_PATH, "experiment_results")  # where the sims live
FIGURES_DIR <- file.path(BASE_PATH, "experiment_results")      # where you save PDFs (keep existing convention)
PLOTS_DIR   <- file.path(FIGURES_DIR, "plots")


dir.create(PLOTS_DIR,   recursive = TRUE, showWarnings = FALSE)

# Canonical helpers
sim_file <- function(hard, pos, n) {
  file.path(RESULTS_DIR, "mainSimResults", "mainSimResults",
            paste0("simsLRR", hard, pos, "n", n, "_xgboost"))
}
fig_file <- function(stem, pos, hard) {
  file.path(FIGURES_DIR, paste0(stem, "pos=", pos, "hard=", hard, ".pdf"))
}
plot_file <- function(stem, pos, hard, lrnr) {
  file.path(PLOTS_DIR, paste0(stem, "pos=", pos, "hard=", hard, "_", lrnr, ".pdf"))
}

# ------------------------------------------------------------
library(ggplot2)
library(data.table)

ns <- sort(c(5000, 500, 1000, 2500))
hard_list <- c(FALSE)
pos_list  <- c(TRUE)
use_oracle_sieve <- FALSE

for (pos in pos_list) {
  for (hard in hard_list) {

    ({
      sims_list <- lapply(ns, function(n) {

        # load + compute dt2 (exactly your logic, just path replaced)
        try({ load(sim_file(hard, pos, n)) })

        simresults <- get("simresults")

        keep <- sapply(simresults, function(item) {
          try({
            force(item$risk_IPW)
            TRUE
          })
          FALSE
        })
        simresults <- simresults[keep]

        ipwrisks   <- rowMeans(do.call(cbind, lapply(simresults, `[[`, "risk_IPW")))
        substrisks <- rowMeans(do.call(cbind, lapply(simresults, `[[`, "risk_subst")))

        lrnr_names <- names(simresults[[1]]$risk_IPW)
        lrnr_names[grep("CV", lrnr_names)] <- c("Lrnr_xgboost_cv", "Lrnr_rf_cv")
        lrnr_names_no_sieve <- lrnr_names

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

        dt$type[!is.na(as.numeric(dt$degree))] <- "sieve"
        dt$type[ is.na(as.numeric(dt$degree))] <- dt$degree[is.na(as.numeric(dt$degree))]

        dt <- dt[dt$degree > 0]

        # add baselines
        tmp <- data.table(risks_oracle = ipwrisks, lrnr_full = lrnr_names_no_sieve,
                          lrnr = lrnr_names_no_sieve, degree = "ipw")
        dt <- rbind(dt, tmp, fill = TRUE)

        tmp <- data.table(risks_oracle = substrisks, lrnr_full = lrnr_names_no_sieve, degree = "subst")
        dt <- rbind(dt, tmp, fill = TRUE)

        # re-derive lrnr label after rbind
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

        dt$type[!is.na(as.numeric(dt$degree))] <- "sieve"
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

      dt <- rbindlist(sims_list)

      options(repr.plot.width = 20, repr.plot.height = 10)

      dt[(dt$type == "sieve"), "type"] <- "EP-learner (*)"
      dt[(dt$type == "ipw"),   "type"] <- "IPW-learner"
      dt[(dt$type == "subst"), "type"] <- "T-learner"
      dt <- dt[(dt$type != "xgboost-Ensemble"), ]

      # -----------------------------
      # XGBoost facet plot
      # -----------------------------
      dt_tmp <- dt
      dt_tmp <- dt_tmp[grep("xgboost", dt$lrnr), ]
      max_depth <- stringr::str_match(dt_tmp$lrnr, "xgboost_10_1_([0-9]+)")[, 2]
      max_depth[is.na(max_depth)] <- "cv"
      dt_tmp$lrnr <- paste0("xgboost (", "max_depth=", max_depth, ")")
      dt_tmp <- dt_tmp[max_depth %in% c("1", "2", "4", "6", "cv"), ]
      dt_tmp <- rbind(dt_tmp, dt[grep("causalforest", dt$lrnr), ])

      plt <- ggplot(dt_tmp, aes(x = n, y = risks_best, group = type, color = type, linetype = type)) +
        geom_line() +
        facet_wrap(~lrnr, scales = "free") +
        theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 1)) +
        ylab("MSE") +
        scale_y_log10(limits = range(dt_tmp$risks_best)) +
        scale_x_log10(breaks = c(500, 1000, 2500, 5000, 10000)) +
        xlab("Sample Size (n)") + ylab("Mean-Squared-Error (MSE)") +
        theme_bw() + labs(color = "Method", group = "Method", linetype = "Method")

      ggsave(fig_file("performancePlot_xgboost_LRR", pos, hard), plot = plt, width = 8, height = 7)

      # -----------------------------
      # RF facet plot
      # -----------------------------
      dt_tmp <- dt[!(dt$lrnr %in% c("glm", "earth", "gam3", "gam4", "gam5")), ]
      dt_tmp <- dt_tmp[grep("rf", dt$lrnr), ]
      max_depth <- stringr::str_match(dt_tmp$lrnr, "([0-9]+)_xg$")[, 2]
      max_depth[is.na(max_depth)] <- "cv"
      dt_tmp$lrnr <- paste0("Random forest (", "max_depth=", max_depth, ")")
      dt_tmp <- dt_tmp[max_depth %in% c("5", "7", "9", "cv"), ]

      plt <- ggplot(dt_tmp, aes(x = n, y = risks_best, group = type, color = type, linetype = type)) +
        geom_line() +
        facet_wrap(~lrnr, scales = "free") +
        theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 1)) +
        ylab("MSE") +
        scale_y_log10(limits = range(dt_tmp$risks_best)) +
        scale_x_log10(breaks = c(500, 1000, 2500, 5000, 10000)) +
        xlab("Sample Size (n)") + ylab("Mean-Squared-Error (MSE)") +
        theme_bw() + labs(color = "Method", group = "Method", linetype = "Method")

      ggsave(fig_file("performancePlot_ranger_LRR", pos, hard), plot = plt, width = 8, height = 7)

      dt_rf <- dt_tmp

      # -----------------------------
      # Small per-learner plots (same outputs, but path fixed)
      # -----------------------------
      dt2 <- rbindlist(sims_list)
      dt2[(dt2$type == "sieve"), "type"] <- "EP-learner (*)"
      dt2[(dt2$type == "ipw"),   "type"] <- "IPW-learner"
      dt2[(dt2$type == "subst"), "type"] <- "T-learner"

      dt_tmp <- dt2[grep("xgboost", dt2$lrnr), ]
      max_depth <- stringr::str_match(dt_tmp$lrnr, "([0-9]+)_0")[, 2]
      max_depth[is.na(max_depth)] <- "cv"
      dt_tmp$lrnr <- paste0("xgboost (", "max_depth=", max_depth, ")")
      dt_tmp <- dt_tmp[max_depth %in% c("1", "5", "cv"), ]
      dt_xg <- dt_tmp

      dt_small_all <- rbind(dt_xg, dt_rf)

      for (lrnr in unique(dt_small_all$lrnr)) {
        dt_small <- as.data.frame(dt_small_all)
        dt_small <- dt_small[dt_small$type != "Oracle DR-Learner", ]

        plt <- ggplot(dt_small[dt_small$lrnr %in% lrnr, ],
                      aes(x = n, y = risks_best, group = type, color = type, linetype = type)) +
          geom_line(size = 0.75) +
          theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 1)) +
          ylab("MSE") +
          scale_y_log10(limits = c(min(1e-1, min(dt_small$risks_best)), max(dt_small$risks_best))) +
          scale_x_log10(breaks = c(500, 1000, 2500, 5000, 10000)) +
          facet_wrap(~lrnr, scales = "free") +
          xlab("Sample Size (n)") + ylab("Mean-Squared-Error (MSE)") +
          theme_bw() + labs(color = "Method", group = "Method", linetype = "Method") +
          theme(axis.text = element_text(size = 14),
                strip.text.x = element_text(size = 16),
                legend.text = element_text(size = 12),
                legend.title = element_text(size = 12),
                axis.title = element_text(size = 12, face = "bold"))

        labels <- c("IPW-learner", "EP-learner (*)", "T-learner")
        colors <- c("#619CFF", "#00BA38", "#F8766D", "#E76BF3")[-1]
        linetypes <- c("longdash", "dashed", "solid", "dotted")[-1]
        names(colors) <- labels
        names(linetypes) <- labels

        plt <- plt +
          scale_colour_manual(values = colors) +
          scale_linetype_manual(values = linetypes) +
          theme(legend.key.height = unit(0.5, "cm"),
                legend.key.width  = unit(1, "cm"),
                legend.position = "none")

        ggsave(plot_file("performancePlot_LRR_tree_", pos, hard, lrnr),
               plot = plt, width = 4, height = 4)
      }

    })
  }
}
