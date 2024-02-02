library(ggplot2)
library(data.table)
ns <- c(  500, 1000, 2500 , 5000)
ns <- sort(ns)
hard_list <-   c(T,F)
pos_list <-  c(T,F)
use_oracle_sieve <- F
for(pos in pos_list){
  for(hard in hard_list) {
    try({
    sims_list <- lapply(ns, function(n) {
      try({load(paste0("mainSimResults2/mainSimResults2/simsCATE", hard, pos,  "n", n, "_gam"))
      simresults <- get(paste0("simresults" ))
      simresults <- simresults[sapply(simresults, is.list)]

      onestepbenchoracle <- rowMeans(do.call(cbind, lapply(simresults, `[[`, "CATEonestepbenchoracle")))
      onestepbench  <- rowMeans(do.call(cbind, lapply(simresults, `[[`, "CATEonestepbench")))

       substrisks  <- rowMeans(do.call(cbind, lapply(simresults, `[[`, "risk_subst")))
       cvsubstrisks  <- rowMeans(do.call(cbind, lapply(simresults, `[[`, "risk_subst_cv")))



       degree <- (stringr::str_match(simresults[[1]]$sieve$sieve_names, "fourier[_.]basis_([0-9]_[0-9])")[,2])
       if(any(is.na(degree))){
       degree[is.na(degree)] <- (stringr::str_match(simresults[[1]]$sieve$sieve_names[is.na(degree)], "fourier[._]basis_([0-9])")[,2])
}
     #degree[grep("no_sieve", lrnrs_full)] <- "0"
      uniq_degrees <- sort(unique(degree))


      lrnr_names <- names(simresults[[1]]$CATEonestepbench) #simresults[[1]]$sieve[[1]]
      lrnr_names <- unlist(lapply(lrnr_names, function(name) {
        paste0(name , paste0("_fourier_basis_", uniq_degrees, "_plugin")) #c( "_no_sieve.plugin", paste0("_fourier_basis_", uniq_degrees, "_plugin")))
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
      dt$degree <- (stringr::str_match(dt$lrnr_full, "fourier[._]basis_([0-9]_[0-9])")[,2])
      if(any(is.na(dt$degree))){
      dt$degree[is.na(dt$degree)] <- (stringr::str_match(dt$lrnr_full[is.na(dt$degree)], "fourier[._]basis_([0-9])")[,2])
}
      #dt$degree[grep("no_sieve", dt$lrnr_full)] <- 0

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
      dt$lrnr <- gsub(".fourier_basis.+", "", dt$lrnr)
      dt$lrnr <- gsub("_no_sieve.+", "", dt$lrnr)
      dt$lrnr <- gsub(".no_sieve.+", "", dt$lrnr)
      dt$lrnr <- gsub("_fourier.basis.+", "", dt$lrnr)
      dt$lrnr <- gsub(".fourier.basis.+", "", dt$lrnr)


      dt$type  <- "Sieve-Plugin" #[!is.na(as.numeric(dt$degree))]
      #dt$type[is.na(as.numeric(dt$degree))] <- dt$degree[is.na(as.numeric(dt$degree))]

      dt <- dt[dt$degree != "0"]

     # tmp <- dt[, cv_sieve_risk := risks_oracle[which.min(cvrisksDR)], by = c("lrnr", "iter")]
    #  tmp <- tmp[, oracle_sieve_risk := risks_oracle[which.min(risks_oracle)], by = c("lrnr", "iter")]
    #  tmp <- tmp[!duplicated(paste0(degree, lrnr, iter )),]


     #  tmp <- dt[, cv_sieve_risk := which.min(cvrisksDR), by = c("lrnr", "iter")]
    #   tmp <- tmp[, oracle_sieve_risk := which.min(risks_oracle), by = c("lrnr", "iter")]
    #   tmp <- tmp[!duplicated(paste0(degree, lrnr, iter )),]
    # dt <- tmp


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
       tmp <- data.table(risks_oracle = substrisks, lrnr_full = names(onestepbench), degree = "Substitution")
       dt <- rbind(dt, tmp, fill = T)
       dt$type[is.na(dt$type)]  <- dt$degree[is.na(dt$type)]

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
      dt$lrnr <- gsub(".fourier_basis.+", "", dt$lrnr)
      dt$lrnr <- gsub("_no_sieve.+", "", dt$lrnr)
      dt$lrnr <- gsub(".no_sieve.+", "", dt$lrnr)
      dt$lrnr <- gsub("_fourier.basis.+", "", dt$lrnr)
      dt$lrnr <- gsub(".fourier.basis.+", "", dt$lrnr)

    #  dt$type[!is.na(as.numeric(dt$degree))] <- "Sieve-Plugin"
     # dt$type[is.na(as.numeric(dt$degree))] <- dt$degree[is.na(as.numeric(dt$degree))]

      print(unique(dt$lrnr_full))

      if(!use_oracle_sieve){
        dt[dt$type=="Sieve-Plugin", risks_oracle := risks_oracle[which.min(cvrisksDR)], by = c("lrnr", "type", "iter")]

      } else {
        #dt[!is.na(as.numeric(dt$degree)), risks_oracle := min(risks_oracle), by = c("iter", "lrnr", "type")]
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
    dt$lrnr[grep("glm", dt$lrnr)] <- "GLM"
    dt$lrnr[grep("earth", dt$lrnr)] <- "MARS (earth)"
    s <- stringr::str_match(dt$lrnr[grep("gam",dt$lrnr)], "s(.+)_")[,2]
    dt$lrnr[grep("gam",dt$lrnr)] <- paste0("GAM (s=", s, ")")
    dt$lrnr[grep("GAM",dt$lrnr)][is.na(s)] <- paste0("GAM (s=cv)")

    dt_tmp<-dt



    dt_tmp <- dt_tmp[!(dt_tmp$type == "Substitution"),]
    dt_tmp[(dt_tmp$type == "Sieve-Plugin"),"type"] <- "EP-Learner (*)"
    dt_tmp[(dt_tmp$type == "One-step"),"type"] <- "DR-Learner"
    dt_tmp[(dt_tmp$type == "Oracle one-step"),"type"] <- "Oracle DR-Learner"
    dt_tmp[(dt_tmp$type == "Substitution-CV"),"type"] <- "T-Learner (CV)"

    plt <- ggplot(dt_tmp, aes(x = n, y = risks_best, group = type, color = type, linetype = type)) + geom_line(size = 0.5) +
      facet_wrap(~lrnr, scales = "free") + theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1)) + ylab("MSE") + scale_y_log10(limits = c(min(1e-1, min(dt_tmp$risks_best)), max(dt_tmp$risks_best)))  +  scale_x_log10(breaks = c(500, 1000, 2500, 5000, 10000))
    plt <- plt + xlab("Sample Size (n)") + ylab("Mean-Squared-Error (MSE)") + theme_bw() + labs(color = "Method", group = "Method", linetype = "Method")
    plt <- plt +  theme_bw() + theme(axis.text=element_text(size=12),
                                     legend.text=element_text(size=12),
                                     legend.title=element_text(size=12),
                                     axis.title=element_text(size=14,face="bold"))
    plt <- plt + theme(legend.justification = c(0.9, 0), legend.position = c(0.9, 0))

    ggsave(paste0("mainSimResults/performancePlot_CATE_GAM_", "pos=",pos, "hard=",hard, ".pdf"), width = 8, height = 7)

    dt <- rbindlist(sims_list)
    dt$lrnr[grep("glm", dt$lrnr)] <- "GLM"
    dt$lrnr[grep("earth", dt$lrnr)] <- "MARS (earth)"
    s <- stringr::str_match(dt$lrnr[grep("gam",dt$lrnr)], "s(.+)_")[,2]
    dt$lrnr[grep("gam",dt$lrnr)] <- paste0("GAM (s=", s, ")")
    dt$lrnr[grep("GAM",dt$lrnr)][is.na(s)] <- paste0("GAM (s=cv)")
    dt_tmp <- dt
    keep <- dt_tmp$lrnr %in% c("MARS (earth)", "GAM (s=1)", "GAM (s=3)", "GAM (s=cv)")
    dt_tmp <- dt_tmp[keep,]
    dt_tmp <- dt_tmp[!(dt_tmp$type == "Substitution"),]
    dt_tmp[(dt_tmp$type == "Sieve-Plugin"),"type"] <- "EP-Learner (*)"
    dt_tmp[(dt_tmp$type == "One-step"),"type"] <- "DR-Learner"
    dt_tmp[(dt_tmp$type == "Oracle one-step"),"type"] <- "Oracle DR-Learner"
    dt_tmp[(dt_tmp$type == "Substitution-CV"),"type"] <- "T-Learner (CV)"


    plt <- ggplot(dt_tmp, aes(x = n, y = risks_best, group = type, color = type, linetype = type)) + geom_line(size = 0.75) +
      facet_wrap(~lrnr, scales = "free") + theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1)) + ylab("MSE") + scale_y_log10(limits = c(min(1e-1, min(dt_tmp$risks_best)), max(dt_tmp$risks_best)))  +  scale_x_log10(breaks = c(500, 1000, 2500, 5000, 10000))
    plt <- plt + xlab("Sample Size (n)") + ylab("Mean-Squared-Error (MSE)") + theme_bw() + labs(color = "Method", group = "Method", linetype = "Method")
    plt <- plt +  theme_bw() + theme(axis.text=element_text(size=12),
                                     legend.text=element_text(size=12),
                                     legend.title=element_text(size=12),
                                     axis.title=element_text(size=14,face="bold"))
    plt <- plt + theme(legend.justification = c(0.9, 0), legend.position = c(0.95, 0.23))
plt
    ggsave(paste0("mainSimResults/performancePlot_CATE_GAM_", "pos=",pos, "hard=",hard, "_mainpaper.pdf"), width = 8, height = 7)



    for(lrnr in unique(dt_tmp$lrnr)) {
      dt_tmp <- as.data.frame(dt_tmp)

      dt_tmp <- dt_tmp[dt_tmp$type != "Oracle DR-Learner",]
      plt <- ggplot(dt_tmp[dt_tmp$lrnr %in% lrnr,], aes(x = n, y = risks_best, group = type, color = type, linetype = type)) + geom_line(size = 0.75)  + theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1)) + ylab("MSE") + scale_y_log10(limits = c(min(1e-1, min(dt_tmp$risks_best)), max(dt_tmp$risks_best)))  +  scale_x_log10(breaks = c(500, 1000, 2500, 5000, 10000)) +
        facet_wrap(~lrnr, scales = "free")
      plt <- plt + xlab("Sample Size (n)") + ylab("Mean-Squared-Error (MSE)") + theme_bw() + labs(color = "Method", group = "Method", linetype = "Method")
      plt <- plt +  theme_bw() + theme(axis.text=element_text(size=14),
                                       strip.text.x = element_text(size = 20),
                                       legend.text=element_text(size=12),
                                       legend.title=element_text(size=12),
                                       axis.title=element_text(size=12,face="bold"))
      plt <- plt + theme(legend.justification = c(0.1, 0), legend.position = c(0.1, 0.1))
      labels <- c("Causal-Forest", "DR-Learner", "EP-Learner (*)", "T-Learner (CV)" )[-1]
      colors <- c("#619CFF", "#00BA38", "#F8766D", "#E76BF3")[-1]
      linetypes <- c("longdash" ,"dashed" , "solid", "dotted")[-1]
      names(colors) <- labels
      names(linetypes) <- labels
      plt <- plt +  scale_colour_manual(  values =   colors)
      plt <- plt + scale_linetype_manual(  values = linetypes)
      plt <- plt +
        theme(legend.key.height= unit(0.5, 'cm'),
              legend.key.width= unit(1, 'cm')) + theme(legend.position = "none")
      ggsave(paste0("mainSimResults/plots/performancePlot_CATE_GAM_", "pos=",pos, "hard=",hard, "_", lrnr,".pdf"), width = 4, height = 4)

    }

    })
  }
}
