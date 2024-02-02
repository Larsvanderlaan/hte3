library(ggplot2)
library(data.table)
ns <- c(    500, 1000, 2500, 5000)
ns <- sort(ns)
hard_list <-   c(T,F)
pos_list <-  c(T,F)
use_oracle_sieve <-F


for(pos in pos_list){
  for(hard in hard_list) {
    try({
    sims_list <- lapply(ns, function(n) {
      try({load(paste0("mainSimResults2/mainSimResults2/simsCATE", hard, pos,  "n", n, "_xgboost"))
      simresults <- get(paste0("simresults"))
      simresults <- simresults[sapply(simresults, is.list)]

      onestepbenchoracle <- rowMeans(do.call(cbind, lapply(simresults, `[[`, "CATEonestepbenchoracle")))
      onestepbench  <- rowMeans(do.call(cbind, lapply(simresults, `[[`, "CATEonestepbench")))

      causalforestrisks <- rowMeans(do.call(cbind, lapply(simresults, `[[`, "risk_cf")))
      substrisks  <- rowMeans(do.call(cbind, lapply(simresults, `[[`, "risk_subst")))
      cvsubstrisks  <- rowMeans(do.call(cbind, lapply(simresults, `[[`, "risk_subst_cv")))

      uniq_degrees <- c("0_0", "1_0", "2_0", "3_0", "1_1")


      lrnr_names <- names(simresults[[1]]$CATEonestepbench) #simresults[[1]]$sieve[[1]]
      lrnr_names <- unlist(lapply(lrnr_names, function(name) {
        paste0(name , c( "_no_sieve.plugin", paste0("_fourier_basis_", uniq_degrees, "_plugin")))
      }))
      iter <- rep(1:length(simresults), each = length(lrnr_names))
      cvrisksDRoracle <- unlist( lapply(simresults, function(item) {
        item$sieve$cvrisksDRoracle
      }))

      cvrisksDR <- unlist(lapply(seq_along(simresults), function(index) {
        item <- simresults[[index]]
        as.vector(item$sieve$cvrisksDR)

      }))


      risks_oracle <- unlist( lapply(simresults, function(item) {
        item$sieve$risks_oracle
      }))

      dt <- data.table(iter, lrnr_full = lrnr_names, cvrisksDR, cvrisksDRoracle, risks_oracle)
      dt$degree <- as.numeric(stringr::str_match(dt$lrnr_full, "fourier_basis_([0-9]+)")[,2])
      dt$degree[grep("no_sieve", dt$lrnr_full)] <- 0

      dt$lrnr[ grep("gam3", dt$lrnr_full)] <- "gam3"
      dt$lrnr[ grep("gam4", dt$lrnr_full)] <- "gam4"
      dt$lrnr[ grep("gam5", dt$lrnr_full)] <- "gam5"
      dt$lrnr[ grep("glm", dt$lrnr_full)] <- "glm"
      dt$lrnr[ grep("earth", dt$lrnr_full)] <- "earth"
      dt$lrnr[ grep("earth", dt$lrnr_full)] <- "earth"
      dt$lrnr[ grep("rpart", dt$lrnr_full)] <- "rpart"
      dt$lrnr[ grep("ranger_500_TRUE_none_1_7", dt$lrnr_full)] <- "ranger_7"
      dt$lrnr[ grep("ranger_500_TRUE_none_1_13", dt$lrnr_full)] <- "ranger_13"
      dt$lrnr[ grep("ranger_500_TRUE_none_1_10", dt$lrnr_full)] <- "ranger_10"

      ## dt$lrnr[ grep("xgboost_20_1_7", dt$lrnr_full)] <- "xgboost_7"
      #  dt$lrnr[ grep("xgboost_20_1_5", dt$lrnr_full)] <- "xgboost_5"
      # dt$lrnr[ grep("xgboost_20_1_3", dt$lrnr_full)] <- "xgboost_3"
      dt$lrnr <- gsub("Lrnr_", "", dt$lrnr)
      dt$lrnr <- gsub("_fourier_basis.+", "", dt$lrnr)


      dt$type[!is.na(as.numeric(dt$degree))] <- "Sieve-Plugin"
      dt$type[is.na(as.numeric(dt$degree))] <- dt$degree[is.na(as.numeric(dt$degree))]

      dt <- dt[dt$degree > 0]

      tmp <- dt[, cv_sieve_risk := risks_oracle[which.min(cvrisksDR)], by = c("lrnr", "iter")]
      tmp <- tmp[, oracle_sieve_risk := risks_oracle[which.min(risks_oracle)], by = c("lrnr", "iter")]
      tmp <- tmp[!duplicated(paste0(degree, lrnr, iter )),]


       tmp <- dt[, cv_sieve_risk := which.min(cvrisksDR), by = c("lrnr", "iter")]
       tmp <- tmp[, oracle_sieve_risk := which.min(risks_oracle), by = c("lrnr", "iter")]
       tmp <- tmp[!duplicated(paste0(degree, lrnr, iter )),]



      ###### LATER


      #dt <- data.table(lrnr_full = lrnr_names, cvrisksDRoracle,cvrisksDR,  risks_oracle)
     # dt$degree <- as.numeric(stringr::str_match(dt$lrnr_full, "fourier_basis_([0-9]+)")[,2])
      #dt$degree[grep("no_sieve", dt$lrnr_full)] <- 0
       tmp <- data.table(risks_oracle = onestepbench, lrnr_full = names(onestepbench), degree = "One-step")
       dt <- rbind(dt, tmp, fill = T)
       tmp <- data.table(risks_oracle = onestepbenchoracle, lrnr_full = names(onestepbench), degree = "Oracle one-step")
       dt <- rbind(dt, tmp, fill = T)
       tmp <- data.table(risks_oracle = cvsubstrisks, lrnr_full =names(onestepbench), degree = "Substitution-CV")
       dt <- rbind(dt, tmp, fill = T)
       tmp <- data.table(risks_oracle = causalforestrisks, lrnr_full = names(onestepbench), lrnr = "Causal-Forest", degree = "Causal-Forest")
       dt <- rbind(dt, tmp, fill = T)
       tmp <- data.table(risks_oracle = substrisks, lrnr_full = names(onestepbench), degree = "Substitution")
       dt <- rbind(dt, tmp, fill = T)
      dt$lrnr <- dt$lrnr_full
      dt$lrnr[ grep("gam3", dt$lrnr_full)] <- "gam3"
      dt$lrnr[ grep("gam4", dt$lrnr_full)] <- "gam4"
      dt$lrnr[ grep("gam5", dt$lrnr_full)] <- "gam5"
      dt$lrnr[ grep("glm", dt$lrnr_full)] <- "glm"
      dt$lrnr[ grep("earth", dt$lrnr_full)] <- "earth"
      dt$lrnr[ grep("earth", dt$lrnr_full)] <- "earth"
      dt$lrnr[ grep("rpart", dt$lrnr_full)] <- "rpart"
      dt$lrnr[ grep("ranger_500_TRUE_none_1_7", dt$lrnr_full)] <- "ranger_7"
      dt$lrnr[ grep("ranger_500_TRUE_none_1_13", dt$lrnr_full)] <- "ranger_13"
      dt$lrnr[ grep("ranger_500_TRUE_none_1_10", dt$lrnr_full)] <- "ranger_10"

      ## dt$lrnr[ grep("xgboost_20_1_7", dt$lrnr_full)] <- "xgboost_7"
      #  dt$lrnr[ grep("xgboost_20_1_5", dt$lrnr_full)] <- "xgboost_5"
      # dt$lrnr[ grep("xgboost_20_1_3", dt$lrnr_full)] <- "xgboost_3"

      dt$lrnr <- gsub("Lrnr_", "", dt$lrnr)
      dt$lrnr <- gsub("_fourier_basis.+", "", dt$lrnr)

      dt$type[!is.na(as.numeric(dt$degree))] <- "Sieve-Plugin"
      dt$type[is.na(as.numeric(dt$degree))] <- dt$degree[is.na(as.numeric(dt$degree))]



      if(!use_oracle_sieve){
        dt[!is.na(as.numeric(dt$degree)), risks_oracle := risks_oracle[which.min(cvrisksDR)], by = c("lrnr", "type", "iter")]

      } else {
        dt[!is.na(as.numeric(dt$degree)), risks_oracle := min(risks_oracle), by = c("iter", "lrnr", "type")]
      }
      dt[, risks_best := mean(risks_oracle), by = c("lrnr", "type")]

      #dt[is.na(as.numeric(dt$degree)), risks_best := risks_oracle, by = c("lrnr", "type")]

      dt2 <- dt[,c("lrnr", "risks_best", "type"), with = F]
      dt2 <- unique(dt2)
      dt2$n <- n
      return(dt2)
    })
    })
    sims_list <- sims_list[sapply(sims_list, is.data.frame)]

    dt <- rbindlist(sims_list)


    dt_tmp<-dt
    dt_tmp <- dt_tmp[ grep("xgboost", dt$lrnr), ]
    max_depth <- stringr::str_match(dt_tmp$lrnr, "([0-9]+)_0$")[,2]
    max_depth[is.na(max_depth)] <- "cv"
    dt_tmp$lrnr <- paste0("xgboost (", "max_depth=", max_depth,")")
    dt_tmp <- dt_tmp[max_depth %in% c("1", "2", "3", "5", "cv"),]
    dt_tmp <- rbind(dt_tmp, dt[grep("Causal-Forest", dt$lrnr),])
    dt_tmp <- dt_tmp[!(dt_tmp$type == "Substitution"),]
    dt_tmp[(dt_tmp$type == "Sieve-Plugin"),"type"] <- "EP-Learner (*)"
    dt_tmp[(dt_tmp$type == "One-step"),"type"] <- "DR-Learner"
    dt_tmp[(dt_tmp$type == "Oracle one-step"),"type"] <- "Oracle DR-Learner"
    dt_tmp[(dt_tmp$type == "Substitution-CV"),"type"] <- "T-Learner (CV)"

    plt <- ggplot(dt_tmp, aes(x = n, y = risks_best, group = type, color = type, linetype = type)) + geom_line() +
      facet_wrap(~lrnr, scales = "free") + theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1)) + ylab("MSE") + scale_y_log10(limits = range(dt_tmp$risks_best))  +  scale_x_log10(breaks = c(500, 1000, 2500, 5000, 10000))
    plt <- plt + xlab("Sample Size (n)") + ylab("Mean-Squared-Error (MSE)") + theme_bw() + labs(color = "Method", group = "Method", linetype = "Method")
    plt <- plt +  theme_bw() + theme(axis.text=element_text(size=12),
                                     legend.text=element_text(size=12),
                                     legend.title=element_text(size=12),
                                     axis.title=element_text(size=14,face="bold"))
    plt <- plt + theme(legend.justification = c(0.9, 0), legend.position = c(0.9, 0))

    ggsave(paste0("mainSimResults/performancePlot_CATE_xgboost", "pos=",pos, "hard=",hard, ".pdf"), width = 8, height = 7)


    dt_tmp<-dt[!(dt$lrnr %in% c("glm", "earth", "gam3", "gam4", "gam5")),]
    dt_tmp <- dt_tmp[grep("rf", dt$lrnr), ]
     max_depth <- stringr::str_match(dt_tmp$lrnr, "([0-9]+)_xg$")[,2]
    max_depth[is.na(max_depth)] <- "cv"
    dt_tmp$lrnr <- paste0("ranger (", "max_depth=", max_depth,")")
    dt_tmp <- dt_tmp[max_depth %in% c("3", "5", "7", "9", "cv"),]
    dt_tmp$lrnr <- factor(dt_tmp$lrnr, levels = paste0("ranger (", "max_depth=", c("3", "5", "7", "9", "cv"),")"))

    dt_tmp <- rbind(dt_tmp, dt[grep("causalforest", dt$lrnr),])

    dt_tmp <- dt_tmp[!(dt_tmp$type == "Substitution"),]
    dt_tmp[(dt_tmp$type == "Sieve-Plugin"),"type"] <- "EP-Learner (*)"
    dt_tmp[(dt_tmp$type == "One-step"),"type"] <- "DR-Learner"
    dt_tmp[(dt_tmp$type == "Oracle one-step"),"type"] <- "Oracle DR-Learner"
    dt_tmp[(dt_tmp$type == "Substitution-CV"),"type"] <- "T-Learner (CV)"

    plt <- ggplot(dt_tmp, aes(x = n, y = risks_best, group = type, color = type, linetype = type)) + geom_line() +
      facet_wrap(~lrnr, scales = "free") + theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1)) + ylab("MSE") + scale_y_log10(limits = range(dt_tmp$risks_best))  +  scale_x_log10(breaks = c(500, 1000, 2500, 5000, 10000))
    plt <- plt + xlab("Sample Size (n)") + ylab("Mean-Squared-Error (MSE)") + theme_bw()+ labs(color = "Method", group = "Method", linetype = "Method")
    plt <- plt +  theme_bw() + theme(axis.text=element_text(size=12),
                                     legend.text=element_text(size=12),
                                     legend.title=element_text(size=12),
                                     axis.title=element_text(size=14,face="bold"))
    plt <- plt + theme(legend.justification = c(0.9, 0), legend.position = c(0.9, 0))

      ggsave(paste0("mainSimResults/performancePlot_CATE_ranger", "pos=",pos, "hard=",hard, ".pdf"), width = 8, height = 7)



      dt_tmp<-dt[!(dt$lrnr %in% c("glm", "earth", "gam3", "gam4", "gam5")),]
      dt_tmp <- dt_tmp[grep("rf", dt$lrnr), ]
      max_depth <- stringr::str_match(dt_tmp$lrnr, "([0-9]+)_xg$")[,2]
      max_depth[is.na(max_depth)] <- "cv"
      dt_tmp$lrnr <- paste0("ranger (", "max_depth=", max_depth,")")
      dt_tmp <- dt_tmp[max_depth %in% c( "3", "7" ,"cv"),]
      dt_tmp$lrnr <- factor(dt_tmp$lrnr, levels = paste0("ranger (", "max_depth=", c("3", "5", "7", "9", "cv"),")"))

      dt_tmp <- rbind(dt_tmp, dt[grep("causalforest", dt$lrnr),])

      dt_tmp <- dt_tmp[!(dt_tmp$type == "Substitution"),]
      dt_tmp[(dt_tmp$type == "Sieve-Plugin"),"type"] <- "EP-Learner (*)"
      dt_tmp[(dt_tmp$type == "One-step"),"type"] <- "DR-Learner"
      dt_tmp[(dt_tmp$type == "Oracle one-step"),"type"] <- "Oracle DR-Learner"
      dt_tmp[(dt_tmp$type == "Substitution-CV"),"type"] <- "T-Learner (CV)"
      dt_rf <- dt_tmp
      print(dt_rf)
      dt <- rbindlist(sims_list)


      dt_tmp<-dt
      dt_tmp <- dt_tmp[ grep("xgboost", dt$lrnr), ]
      max_depth <- stringr::str_match(dt_tmp$lrnr, "([0-9]+)_0$")[,2]
      max_depth[is.na(max_depth)] <- "cv"
      dt_tmp$lrnr <- paste0("xgboost (", "max_depth=", max_depth,")")
      dt_tmp <- dt_tmp[max_depth %in% c("1", "5", "cv"),]
      dt_tmp <- rbind(dt_tmp, dt[grep("Causal-Forest", dt$lrnr),])
      dt_tmp <- dt_tmp[!(dt_tmp$type == "Substitution"),]
      dt_tmp[(dt_tmp$type == "Sieve-Plugin"),"type"] <- "EP-Learner (*)"
      dt_tmp[(dt_tmp$type == "One-step"),"type"] <- "DR-Learner"
      dt_tmp[(dt_tmp$type == "Oracle one-step"),"type"] <- "Oracle DR-Learner"
      dt_tmp[(dt_tmp$type == "Substitution-CV"),"type"] <- "T-Learner (CV)"
      dt_xg <- dt_tmp

      dt <- rbind( dt_xg,dt_rf)
      plt <- ggplot(dt, aes(x = n, y = risks_best, group = type, color = type, linetype = type)) + geom_line() +
        facet_wrap(~lrnr, scales = "free") + theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1)) + ylab("MSE") + scale_y_log10(limits = range(dt$risks_best))  +  scale_x_log10(breaks = c(500, 1000, 2500, 5000, 10000))
      plt <- plt + xlab("Sample Size (n)") + ylab("Mean-Squared-Error (MSE)") + theme_bw()+ labs(color = "Method", group = "Method", linetype = "Method")
      plt <- plt +  theme_bw() + theme(axis.text=element_text(size=12),
                                       legend.text=element_text(size=12),
                                       legend.title=element_text(size=12),
                                       axis.title=element_text(size=14,face="bold"))
      plt <- plt + theme(legend.justification = c(0.9, 0), legend.position = c(0.95, 0.2))

      ggsave(paste0("mainSimResults/performancePlot_CATE_tree", "pos=",pos, "hard=",hard, "mainpaper.pdf"), width = 8, height = 7)



      for(lrnr in unique(dt$lrnr)) {
        print(lrnr)
        dt <- as.data.frame(dt)
        dt <- dt[dt$type != "Oracle DR-Learner",]
        dt_tmp <- dt
        plt <- ggplot(dt[dt$lrnr %in% lrnr,], aes(x = n, y = risks_best, group = type, color = type, linetype = type)) + geom_line(size = 0.75)  + theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1)) + ylab("MSE") + scale_y_log10(limits = c(min(1e-1, min(dt_tmp$risks_best)), max(dt_tmp$risks_best)))  +  scale_x_log10(breaks = c(500, 1000, 2500, 5000, 10000)) +
          facet_wrap(~lrnr, scales = "free")
        plt <- plt + xlab("Sample Size (n)") + ylab("Mean-Squared-Error (MSE)") + theme_bw() + labs(color = "Method", group = "Method", linetype = "Method")
        plt <- plt +  theme_bw() + theme(axis.text=element_text(size=14),
                                         strip.text.x = element_text(size = 20),
                                         legend.text=element_text(size=12),
                                         legend.title=element_text(size=12),
                                         axis.title=element_text(size=12,face="bold"))
       # plt <- plt + theme(legend.justification = c(0.6, 0), legend.position = c(0.75, 0.6))
        #plt <- plt +  scale_colour_manual(labels = c("Causal-Forest", "DR-Learner", "EP-Learner (*)", "T-Learner (CV)" ), values =   c("#619CFF", "#00BA38", "#F8766D", "#E76BF3"))
        #plt <- plt + scale_linetype_manual(labels = c("Causal-Forest", "DR-Learner", "EP-Learner (*)", "T-Learner (CV)" ), values = c("longdash" ,"dashed" , "solid", "dotted"))
        labels <- c("Causal-Forest", "DR-Learner", "EP-Learner (*)", "T-Learner (CV)" )
        colors <- c("#619CFF", "#00BA38", "#F8766D", "#E76BF3")
        linetypes <- c("longdash" ,"dashed" , "solid", "dotted")
        names(colors) <- labels
        names(linetypes) <- labels
        plt <- plt +  scale_colour_manual(  values =   colors)
        plt <- plt + scale_linetype_manual( values = linetypes)
        plt <- plt + theme(legend.key.size = unit(2, 'cm'), legend.title = element_blank(), legend.direction  = "horizontal") #legend.position = "horizontal")
        ggsave(paste0("mainSimResults/plots/performancePlot_CATE_tree_", "pos=",pos, "hard=",hard, "_", lrnr,".pdf"), width = 10, height = 10)

        plt <- plt +
          theme(legend.key.height= unit(0.5, 'cm'),
                legend.key.width= unit(1, 'cm'))  + theme(legend.position = "none")
        ggsave(paste0("mainSimResults/plots/performancePlot_CATE_tree_", "pos=",pos, "hard=",hard, "_", lrnr,".pdf"), width = 4, height = 4)
      }

})
  }
}
