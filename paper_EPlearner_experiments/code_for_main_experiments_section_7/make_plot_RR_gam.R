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

RESULTS_DIR <- file.path(BASE_PATH, "experiment_results")  # where sims live
FIGURES_DIR <- file.path(BASE_PATH, "experiment_results")      # keep your existing output folder name here
PLOTS_DIR   <- file.path(FIGURES_DIR, "plots")


dir.create(PLOTS_DIR,   recursive = TRUE, showWarnings = FALSE)

# Canonical helpers (so there are no hard-coded strings later)
sim_file <- function(hard, pos, n) {
  file.path(RESULTS_DIR, "mainSimResults", "mainSimResults",
            paste0("simsLRR", hard, pos, "n", n, "_gam"))
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

ns <- sort(c(500, 1000, 2500, 5000))
ns_label <- ns
names(ns_label) <- ns

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

          cvsubstrisks <- rowMeans(do.call(cbind, lapply(simresults, `[[`, "risk_subst_cv")))
          ipwrisks     <- rowMeans(do.call(cbind, lapply(simresults, `[[`, "risk_IPW")))
          substrisks   <- rowMeans(do.call(cbind, lapply(simresults, `[[`, "risk_subst")))

          if (is.na(substrisks[length(substrisks)])) {
            substrisks[length(substrisks)] <- substrisks[1]
          }

          lrnr_names <- names(simresults[[1]]$risk_IPW)
          lrnr_names[grep("CV", lrnr_names)] <- "Lrnr_gam_scv_x"
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
          dt$lrnr <- gsub(".fourier_basis.+", "", dt$lrnr)
          dt$lrnr <- gsub("_no_sieve.+", "", dt$lrnr)
          dt$lrnr <- gsub(".no_sieve.+", "", dt$lrnr)
          dt$lrnr <- gsub("_fourier.basis.+", "", dt$lrnr)
          dt$lrnr <- gsub(".fourier.basis.+", "", dt$lrnr)

          dt$type[!is.na(as.numeric(dt$degree))] <- "sieve"
          dt$type[ is.na(as.numeric(dt$degree))] <- dt$degree[is.na(as.numeric(dt$degree))]

          dt <- dt[dt$degree > 0]

          # ---- add baselines (no sieve)
          tmp <- data.table(risks_oracle = ipwrisks,     lrnr_full = lrnr_names_no_sieve, degree = "IPW")
          dt <- rbind(dt, tmp, fill = TRUE)

          tmp <- data.table(risks_oracle = cvsubstrisks, lrnr_full = lrnr_names_no_sieve, degree = "xgboost-Ensemble")
          dt <- rbind(dt, tmp, fill = TRUE)

          tmp <- data.table(risks_oracle = substrisks,   lrnr_full = lrnr_names_no_sieve, degree = "Substitution")
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
          dt$lrnr <- gsub(".fourier_basis.+", "", dt$lrnr)
          dt$lrnr <- gsub("_no_sieve.+", "", dt$lrnr)
          dt$lrnr <- gsub(".no_sieve.+", "", dt$lrnr)
          dt$lrnr <- gsub("_fourier.basis.+", "", dt$lrnr)
          dt$lrnr <- gsub(".fourier.basis.+", "", dt$lrnr)

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
      })

      sims_list <- sims_list[sapply(sims_list, is.data.frame)]
      dt <- rbindlist(sims_list)

      # ---- plot prep (same as your logic)
      dt <- dt[-grep("glm", dt$lrnr)]
      s <- stringr::str_match(dt$lrnr, "s(.+)_")[, 2]
      s[is.na(s)] <- "cv"
      dt$lrnr <- paste0("GAM (s=", s, ")")
      dt <- dt[s %in% c("1", "3", "5", "cv"), ]

      dt_tmp <- dt
      dt_tmp <- dt_tmp[(dt_tmp$type != "xgboost-Ensemble"), ]
      dt_tmp[(dt_tmp$type == "sieve"),        "type"] <- "EP-learner (*)"
      dt_tmp[(dt_tmp$type == "IPW"),          "type"] <- "IPW-learner"
      dt_tmp[(dt_tmp$type == "xgboost-Ensemble"), "type"] <- "T-learner (CV)"
      dt_tmp[(dt_tmp$type == "Substitution"), "type"] <- "T-learner"

      plt <- ggplot(dt_tmp, aes(x = n, y = risks_best, group = type, color = type, linetype = type)) +
        geom_line() +
        facet_wrap(~lrnr, scales = "free") +
        theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 1)) +
        ylab("MSE") +
        scale_y_log10(limits = c(min(1e-1, min(dt_tmp$risks_best)), max(dt_tmp$risks_best))) +
        scale_x_log10(breaks = c(500, 1000, 2500, 5000, 10000)) +
        xlab("Sample Size (n)") + ylab("Mean-Squared-Error (MSE)") +
        theme_bw() + labs(color = "Method", group = "Method", linetype = "Method")

      ggsave(fig_file("performancePlot_GAM_LRR_", pos, hard), plot = plt, width = 8, height = 7)

      for (lrnr in unique(dt_tmp$lrnr)) {
        dt_small <- as.data.frame(dt_tmp)
        dt_small <- dt_small[dt_small$type != "Oracle DR-Learner", ]

        plt <- ggplot(dt_small[dt_small$lrnr %in% lrnr, ],
                      aes(x = n, y = risks_best, group = type, color = type, linetype = type)) +
          geom_line(size = 0.75) +
          theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 1)) +
          ylab("MSE") +
          scale_y_log10(limits = c(min(1e-1, min(dt_small$risks_best)), max(dt_small$risks_best))) +
          scale_x_log10(breaks = c(500, 1000, 2500, 5000, 5000, 10000)) +
          facet_wrap(~lrnr, scales = "free") +
          xlab("Sample Size (n)") + ylab("Mean-Squared-Error (MSE)") +
          theme_bw() + labs(color = "Method", group = "Method", linetype = "Method") +
          theme(axis.text = element_text(size = 14),
                strip.text.x = element_text(size = 20),
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
          theme(legend.key.size = unit(1.5, "cm"),
                legend.title = element_blank(),
                legend.direction = "horizontal",
                legend.position = "none")

        ggsave(plot_file("performancePlot_LRR_GAM_", pos, hard, lrnr),
               plot = plt, width = 4, height = 4)
      }

    })
  }
}
