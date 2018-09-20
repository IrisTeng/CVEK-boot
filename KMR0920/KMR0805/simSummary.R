setwd('/home/wd45/KMR0806/Temp')

library(ggplot2)
library(e1071)
library(dplyr)
library(reshape2)

if (!dir.exists("./plots/")) {
  dir.create("./plots/")
} 


null_summary <- NULL
tau_expmean <- NULL
dif_mean <- NULL
V0_inv_logmean <- NULL
K12_mean <- NULL
for (d in 1:5) {
  
  test_sum <- NULL
  weight_sum <- NULL
  tau_expsum <- NULL
  dif_sum <- NULL
  V0_inv_logsum <- NULL
  K12_sum <- NULL
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
        null_stat <- NULL
        weight_mat <- NULL
        for (i in 1:200) {
          
          temp <- res[[i]]
          test_sum <- rbind(test_sum, c(as.numeric(temp$score_chi), 
                                        int_effect, mode, strategy))
          
          null_stat <- c(null_stat, temp$bs_test)
          weight_mat <- rbind(weight_mat, temp$u_weight)
          
          tau_expsum <- rbind(tau_expsum, c(exp(as.numeric(temp$tau_ret)), 
                                            int_effect, mode, strategy))
          
          dif_sum <- rbind(dif_sum, c(as.numeric(temp$dif_det), 
                                      int_effect, mode, strategy))
          
          V0_inv_temp <- ifelse(is.finite(as.numeric(temp$V0_inv_det)), 
                                as.numeric(temp$V0_inv_det), 999)
          
          V0_inv_logsum <- rbind(V0_inv_logsum, c(V0_inv_temp, 
                                                  int_effect, mode, strategy))
          
          K12_sum <- rbind(K12_sum, c(as.numeric(temp$K12_det), 
                                      int_effect, mode, strategy))
        }
        
        null_summary <- rbind(null_summary, c(mean(null_stat), 
                                              var(null_stat), skewness(null_stat), 
                                              int_effect, mode, strategy, d))
        
        weight_sum <- rbind(weight_sum, c(colMeans(weight_mat), 
                                          int_effect, mode, strategy))
        
        tau_expmean <- rbind(tau_expmean, c(mean(as.numeric(tau_expsum[, 1])), 
                                            int_effect, mode, strategy, d))
        dif_mean <- rbind(dif_mean, c(mean(as.numeric(dif_sum[, 1])), 
                                      int_effect, mode, strategy, d))
        
        V0_inv_meantemp <- mean(as.numeric(V0_inv_logsum[
          which(as.numeric(V0_inv_logsum[, 1]) != 999), 1]))
        V0_inv_logmean <- rbind(V0_inv_logmean, c(V0_inv_meantemp, 
                                                  int_effect, mode, strategy, d))
        K12_mean <- rbind(K12_mean, c(mean(as.numeric(K12_sum[, 1])), 
                                      int_effect, mode, strategy, d))
      }
    }
  }
  
  stat_ss <- list()
  stat_ss[[1]] <- tau_expsum
  stat_ss[[2]] <- dif_sum
  stat_ss[[3]] <- V0_inv_logsum
  stat_ss[[4]] <- K12_sum
  stat_ss[[5]] <- test_sum
  
  for (i in 1:5) {
    
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
    
    if (i == 5) {
      fig <- fig + xlim(0, 4e4)
    } else if (i == 1) {
      if (d %in% c(1, 2)) {
        fig <- fig + xlim(0, log(1.3))
      } else if (d == 4) {
        fig <- fig + xlim(0, log(1.13))
      } else {
        fig <- fig + xlim(0, log(1.04))
      }
    }
    
    ggsave(paste0("./plots/library_", d, "_density_", i, ".pdf"), 
           fig, width = 14, height = 9)
    
  }
  
  weight_sum <- as.data.frame(weight_sum)
  colnames(weight_sum)[c((ncol(weight_sum) - 2):ncol(weight_sum))] <-
    c("int_effect", "mode", "strategy")
  
  weight_sum <- weight_sum %>% mutate_at(c(1:(ncol(weight_sum) - 2)),
                                         funs(as.character)) %>% mutate_at(c(1:(ncol(weight_sum) - 2)),
                                                                           funs(as.numeric))
  
  weight_melt <- melt(weight_sum, id.vars = c('int_effect', 'mode', 'strategy'),
                      variable.name = 'ind', value.name = 'weight')
  
  weight_melt$mode <- as.factor(weight_melt$mode)
  weight_melt$ind <- as.factor(weight_melt$ind)
  weight_melt$strategy <- as.factor(weight_melt$strategy)
  
  fig <- ggplot(data = weight_melt, aes(x = int_effect, y = weight, color = ind))+
    facet_grid(strategy ~ mode)+ geom_line()+ geom_point(size = 0.6) + 
    labs(x = 'int_effect', y = 'weight', title = paste0("library_", d)) +
    theme_set(theme_bw()) + theme(panel.grid = element_blank()) +
    theme(axis.title.x = element_text(size = 15), axis.title.y = element_text(size = 15),
          plot.title = element_text(size = 16), axis.text.x=element_text(size = 12),
          axis.text.y=element_text(size = 12), panel.grid =element_blank()) +
    theme(legend.title = element_text(size = 12, face = "bold")) +
    theme(legend.text = element_text(size = 12))
  ggsave(paste0("./plots/library_", d, "_weight.pdf"),
         fig, width = 14, height = 9)
  
}

# geom_errorbar(width = .1, aes(y = weight, ymax=weight, ymin=weight, color = ind)) 
stat_mean <- list()
stat_mean[[1]] <- tau_expmean
stat_mean[[2]] <- dif_mean
stat_mean[[3]] <- V0_inv_logmean
stat_mean[[4]] <- K12_mean

for (i in 1:4) {
  
  file_mean <- stat_mean[[i]]
  file_mean <- as.data.frame(file_mean)
  colnames(file_mean) <- c("mean", "int_effect", "mode", "strategy", "d")
  
  file_mean$mean <- as.numeric(as.character(file_mean$mean))
  file_mean$int_effect <- as.numeric(as.character(file_mean$int_effect))
  file_mean$mode <- as.factor(file_mean$mode)
  file_mean$d <- as.factor(file_mean$d)
  file_mean$strategy <- as.factor(file_mean$strategy)
  
  fig <- ggplot(data = file_mean, aes(x = int_effect, y = mean, color = d))+ 
    facet_grid(strategy ~ mode) + geom_point(size = 0.6) + geom_line() +
    labs(x = 'int_effect', y = 'mean') + 
    theme_set(theme_bw()) + theme(panel.grid = element_blank()) + 
    theme(axis.title.x = element_text(size = 15), axis.title.y = element_text(size = 15), 
          plot.title = element_text(size = 16), axis.text.x=element_text(size = 12), 
          axis.text.y=element_text(size = 12), panel.grid =element_blank()) + 
    theme(legend.title = element_text(size = 12, face = "bold")) +
    theme(legend.text = element_text(size = 12))
  ggsave(paste0("./plots/mean_", i, ".pdf"),
         fig, width = 14, height = 9)
}


colnames(null_summary) <- c("mean", "var", "skew", 
                            "int_effect", "mode", "strategy", "d")
null_summary <- as.data.frame(null_summary)
null_summary$mean <- as.numeric(as.character(null_summary$mean))
null_summary$var <- as.numeric(as.character(null_summary$var))
null_summary$skew <- as.numeric(as.character(null_summary$skew))
null_summary$int_effect <- as.numeric(as.character(null_summary$int_effect))
null_summary$mode <- as.factor(null_summary$mode)
null_summary$d <- as.factor(null_summary$d)
null_summary$strategy <- as.factor(null_summary$strategy)

## mean
fig <- ggplot(data = null_summary, aes(x = int_effect, y = mean, color = d))+ 
  geom_line() + facet_grid(strategy ~ mode)+ geom_point(size = 0.6) + 
  labs(x = 'int_effect', y = 'mean') + 
  theme_set(theme_bw()) + theme(panel.grid = element_blank()) + 
  theme(axis.title.x = element_text(size = 15), axis.title.y = element_text(size = 15), 
        plot.title = element_text(size = 16), axis.text.x=element_text(size = 12), 
        axis.text.y=element_text(size = 12), panel.grid =element_blank()) + 
  theme(legend.title = element_text(size = 12, face = "bold")) +
  theme(legend.text = element_text(size = 12))
ggsave("./plots/null_dist_mean.pdf", 
       fig, width = 14, height = 9)

## var
fig <- ggplot(data = null_summary, aes(x = int_effect, y = var, color = d))+ 
  geom_line() + facet_grid(strategy ~ mode)+ geom_point(size = 0.6) + 
  labs(x = 'int_effect', y = 'var') + 
  theme_set(theme_bw()) + theme(panel.grid = element_blank()) + 
  theme(axis.title.x = element_text(size = 15), axis.title.y = element_text(size = 15), 
        plot.title = element_text(size = 16), axis.text.x=element_text(size = 12), 
        axis.text.y=element_text(size = 12), panel.grid =element_blank()) + 
  theme(legend.title = element_text(size = 12, face = "bold")) +
  theme(legend.text = element_text(size = 12))
ggsave("./plots/null_dist_var.pdf", 
       fig, width = 14, height = 9)

## skew
fig <- ggplot(data = null_summary, aes(x = int_effect, y = skew, color = d))+ 
  geom_line() + facet_grid(strategy ~ mode)+ geom_point(size = 0.6) + 
  labs(x = 'int_effect', y = 'skew') + 
  theme_set(theme_bw()) + theme(panel.grid = element_blank()) + 
  theme(axis.title.x = element_text(size = 15), axis.title.y = element_text(size = 15), 
        plot.title = element_text(size = 16), axis.text.x=element_text(size = 12), 
        axis.text.y=element_text(size = 12), panel.grid =element_blank()) + 
  theme(legend.title = element_text(size = 12, face = "bold")) +
  theme(legend.text = element_text(size = 12))
ggsave("./plots/null_dist_skew.pdf", 
       fig, width = 14, height = 9)
