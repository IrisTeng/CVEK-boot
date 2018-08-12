setwd('/home/wd45/KMR0812/output/Temp')

library(ggplot2)
library(dplyr)
library(reshape2)

if (!dir.exists("./plots/")) {
  dir.create("./plots/")
} 

lambda_mean <- NULL
sigma2_expmean <- NULL
# tau_expmean <- NULL
# dif_mean <- NULL
# V0_inv_dif_logmean <- NULL
# K12_tr_mean <- NULL
A_tr_mean <- NULL

for (d in 1:5) {
  
  lambda_sum <- NULL
  sigma2_expsum <- NULL
  # tau_expsum <- NULL
  # dif_sum <- NULL
  # V0_inv_dif_logsum <- NULL
  # K12_tr_sum <- NULL
  A_tr_sum <- NULL
  for (mode in c("loocv", "AICc", "GCVc", "gmpml")) {
    for (strategy in c("average", "exp", "erm")) {
      for (int_effect in seq(0, 1, .1)) {
        
        res_key <-
          paste0("mode", mode, 
                 "_strategy", strategy, 
                 "_int", int_effect,
                 "_K", d,
                 ".RData")
        
        load(res_key)

        for (i in 1:200) {
          
          temp <- res[[i]]
          
          lambda_sum <- rbind(lambda_sum, c(as.numeric(temp$lam), 
                                            int_effect, mode, strategy))
          
          sigma2_expsum <- rbind(sigma2_expsum, c(exp(as.numeric(temp$sigma2_hat)), 
                                                 int_effect, mode, strategy))
          
          # tau_expsum <- rbind(tau_expsum, c(exp(as.numeric(temp$tau_hat)), 
          #                                   int_effect, mode, strategy))
          # 
          # dif_sum <- rbind(dif_sum, c(as.numeric(temp$dif_sum), 
          #                             int_effect, mode, strategy))
          # 
          # V0_inv_dif_temp <- ifelse(is.finite(as.numeric(temp$V0_inv_dif_sum)), 
          #                           log(as.numeric(temp$V0_inv_dif_sum)), 999)
          # 
          # V0_inv_dif_logsum <- rbind(V0_inv_dif_logsum, c(V0_inv_dif_temp, 
          #                                                 int_effect, mode, strategy))
          # 
          # K12_tr_sum <- rbind(K12_tr_sum, c(as.numeric(temp$K12_tr), 
          #                                   int_effect, mode, strategy))
          
          A_tr_sum <- rbind(A_tr_sum, c(as.numeric(temp$A_tr), 
                                        int_effect, mode, strategy))
        }
        
        lambda_mean <- rbind(lambda_mean, c(mean(as.numeric(lambda_sum[, 1])), 
                                            int_effect, mode, strategy, d))
        
        sigma2_expmean <- rbind(sigma2_expmean, c(mean(as.numeric(sigma2_expsum[, 1])), 
                                            int_effect, mode, strategy, d))
        
        # tau_expmean <- rbind(tau_expmean, c(mean(as.numeric(tau_expsum[, 1])), 
        #                                     int_effect, mode, strategy, d))
        # 
        # dif_mean <- rbind(dif_mean, c(mean(as.numeric(dif_sum[, 1])), 
        #                               int_effect, mode, strategy, d))
        # 
        # V0_inv_dif_meantemp <- mean(as.numeric(V0_inv_dif_logsum[
        #   which(as.numeric(V0_inv_dif_logsum[, 1]) != 999), 1]))
        # 
        # V0_inv_dif_logmean <- rbind(V0_inv_dif_logmean, c(V0_inv_dif_meantemp, 
        #                                           int_effect, mode, strategy, d))
        # 
        # K12_tr_mean <- rbind(K12_tr_mean, c(mean(as.numeric(K12_tr_sum[, 1])), 
        #                               int_effect, mode, strategy, d))
        
        A_tr_mean <- rbind(A_tr_mean, c(mean(as.numeric(A_tr_sum[, 1])), 
                                            int_effect, mode, strategy, d))
      }
    }
  }
  
  stat_ss <- list()
  stat_ss[[1]] <- lambda_sum
  stat_ss[[2]] <- sigma2_expsum
  # stat_ss[[3]] <- tau_expsum
  # stat_ss[[4]] <- dif_sum
  # stat_ss[[5]] <- V0_inv_dif_logsum
  # stat_ss[[6]] <- K12_tr_sum
  stat_ss[[3]] <- A_tr_sum
  
  for (i in 1:3) {
    
    file <- stat_ss[[i]]
    colnames(file) <- c("value", "int_effect", "mode", "strategy")
    file <- as.data.frame(file)
    file <- mutate(file, ind = ifelse(as.numeric(as.character(int_effect)) > .5, 1, 0))
    file$value <- as.numeric(as.character(file$value))
    
    file$mode <- as.factor(file$mode)
    file$ind <- as.factor(file$ind)
    file$strategy <- as.factor(file$strategy)
    
    fig <- ggplot(data = file, aes(x = value, linetype = ind, color = int_effect))+ 
      geom_density(alpha = .1) + facet_grid(strategy ~ mode)+
      labs(x = 'obsv_stat', y = 'density') + 
      theme_set(theme_bw()) + theme(panel.grid = element_blank()) + 
      theme(axis.title.x = element_text(size = 15), axis.title.y = element_text(size = 15), 
            plot.title = element_text(size = 16), axis.text.x=element_text(size = 12), 
            axis.text.y=element_text(size = 12), panel.grid =element_blank()) + 
      theme(legend.title = element_text(size = 12, face = "bold")) +
      theme(legend.text = element_text(size = 12))
    
    if (i == 2) {
      if (d %in% 1:2) {
        fig <- fig + xlim(1, 1.02)
      } else if (d == 4) {
        fig <- fig + xlim(1, 1.01)
      } else {
        fig <- fig + xlim(1, 1.001)
      }
    } else if (i == 1) {
        fig <- fig + xlim(0, .1)
    } else if (i == 3) {
      if (d %in% c(1, 2, 4)) {
        fig <- fig + xlim(18.84, 19)
      }
    }
    
    ggsave(paste0("./plots/library_", d, "_density_", i, ".pdf"), 
           fig, width = 14, height = 9)
    
  }
}

# stat_mean <- list()
# stat_mean[[1]] <- lambda_mean
# stat_mean[[2]] <- sigma2_expmean
# stat_mean[[3]] <- tau_expmean
# stat_mean[[4]] <- dif_mean
# stat_mean[[5]] <- V0_inv_dif_logmean
# stat_mean[[6]] <- K12_tr_mean
# stat_mean[[7]] <- A_tr_mean
# 
# for (i in 1:7) {
#   
#   file_mean <- stat_mean[[i]]
#   file_mean <- as.data.frame(file_mean)
#   colnames(file_mean) <- c("mean", "int_effect", "mode", "strategy", "d")
#   
#   file_mean$mean <- as.numeric(as.character(file_mean$mean))
#   file_mean$int_effect <- as.numeric(as.character(file_mean$int_effect))
#   file_mean$mode <- as.factor(file_mean$mode)
#   file_mean$d <- as.factor(file_mean$d)
#   file_mean$strategy <- as.factor(file_mean$strategy)
#   
#   fig <- ggplot(data = file_mean, aes(x = int_effect, y = mean, color = d))+ 
#     facet_grid(strategy ~ mode) + geom_point(size = 0.6) + geom_line() +
#     labs(x = 'int_effect', y = 'mean') + 
#     theme_set(theme_bw()) + theme(panel.grid = element_blank()) + 
#     theme(axis.title.x = element_text(size = 15), axis.title.y = element_text(size = 15), 
#           plot.title = element_text(size = 16), axis.text.x=element_text(size = 12), 
#           axis.text.y=element_text(size = 12), panel.grid =element_blank()) + 
#     theme(legend.title = element_text(size = 12, face = "bold")) +
#     theme(legend.text = element_text(size = 12))
#   ggsave(paste0("./plots/mean_", i, ".pdf"),
#          fig, width = 14, height = 9)
# }
