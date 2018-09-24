setwd('/home/xz172/KMR0922_3/output/Temp')

library(ggplot2)
library(dplyr)
library(reshape2)

if (!dir.exists("./plots/")) {
  dir.create("./plots/")
} 

sum_stat <- NULL
for (method in c("polynomial", "rbf")) {
  for (l in c(.5, 1, 1.5)) {
    for (d in 1) {
      weight_sum <- NULL
      for (mode in c("loocv", "AICc", "GCVc", "gmpml")) {
        for (strategy in c("avg")) {
          for (int_effect in c(seq(0, .29, .05), seq(.3, 1, .1))) {
            
            res_key <-
              paste0("true", method, 
                     "_p", 2 * l, 
                     "_l", l,
                     "_mode", mode, 
                     "_strategy", strategy, 
                     "_int", int_effect,
                     "_K", d,
                     ".RData")
            
            load(res_key)
            weight_mat <- NULL
            
            pvalue_temp <- NULL
            pvalue2_temp <- NULL
            lam_temp <- NULL
            K_temp <- NULL
            A_temp <- NULL
            V0_invlog_temp <- NULL
            train_temp <- NULL
            test_temp <- NULL
            trtst_temp <- NULL
            for (i in 1:200) {
              weight_mat <- rbind(weight_mat, c(as.numeric(res[[i]]$u_weight)))
              
              pvalue_temp <- c(pvalue_temp, res[[i]]$pvalue)
              pvalue2_temp <- c(pvalue2_temp, res[[i]]$pvalue2)
              lam_temp <- c(lam_temp, res[[i]]$lam)
              K_temp <- c(K_temp, res[[i]]$K_tr)
              A_temp <- c(A_temp, res[[i]]$A_tr)
              V0_invlog_temp <- c(V0_invlog_temp, log(as.numeric(res[[i]]$V0_inv_tr)))
              train_temp <- c(train_temp, res[[i]]$train_RMSE)
              test_temp <- c(test_temp, res[[i]]$test_RMSE)
              trtst_temp <- c(trtst_temp, res[[i]]$trtst_ratio)
            }
            
            weight_sum <- rbind(weight_sum, c(colMeans(weight_mat), 
                                              int_effect, mode, strategy))
            
            sim_power <- sum(pvalue_temp < .05) / 200
            sim2_power <- sum(pvalue2_temp < .05) / 200
            lam_mean <- round(median(lam_temp), 4)
            K_mean <- round(mean(K_temp), 4)
            A_mean <- round(mean(A_temp), 4)
            V0_invlog_mean <- round(mean(V0_invlog_temp), 4)
            train_mean <- round(mean(train_temp), 4)
            test_mean <- round(mean(test_temp), 4)
            trtst_mean <- round(mean(trtst_temp), 4)
            
            sum_stat <- rbind(sum_stat, c(as.character(method), 
                                          2 * l, l, int_effect, 
                                          as.character(mode), 
                                          as.character(strategy), 
                                          d, sim_power, sim2_power, 
                                          lam_mean, K_mean, 
                                          A_mean, V0_invlog_mean, 
                                          train_mean, test_mean, 
                                          trtst_mean))
            
          }
        }
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
      weight_melt$int_effect <- as.numeric(as.character(weight_melt$int_effect))
      
      fig <- ggplot(data = weight_melt, aes(x = int_effect, y = weight, color = ind))+
        facet_grid(strategy ~ mode)+ geom_line()+ geom_point(size = 0.6) + 
        labs(x = 'int_effect', y = 'weight', title = paste0("library_", d)) +
        theme_set(theme_bw()) + theme(panel.grid = element_blank()) +
        theme(axis.title.x = element_text(size = 15), axis.title.y = element_text(size = 15),
              plot.title = element_text(size = 16), axis.text.x=element_text(size = 12),
              axis.text.y=element_text(size = 12), panel.grid =element_blank()) +
        theme(legend.title = element_text(size = 12, face = "bold")) +
        theme(legend.text = element_text(size = 12))
      ggsave(paste0("./plots/library_", d, "_true", as.character(method), "_p", 2 * l, 
                    "_l", l, "_weight.pdf"),
             fig, width = 14, height = 9)
      
    }
  }
}


for (method in c("polynomial", "rbf")) {
  for (l in c(.5, 1, 1.5)) {
    for (d in 2:5) {
      weight_sum <- NULL
      for (mode in c("loocv", "AICc", "GCVc", "gmpml")) {
        for (strategy in c("avg", "exp", "erm")) {
          for (int_effect in c(seq(0, .29, .05), seq(.3, 1, .1))) {
            
            res_key <-
              paste0("true", method, 
                     "_p", 2 * l, 
                     "_l", l,
                     "_mode", mode, 
                     "_strategy", strategy, 
                     "_int", int_effect,
                     "_K", d,
                     ".RData")
            
            load(res_key)
            weight_mat <- NULL
            
            pvalue_temp <- NULL
            pvalue2_temp <- NULL
            lam_temp <- NULL
            K_temp <- NULL
            A_temp <- NULL
            V0_invlog_temp <- NULL
            train_temp <- NULL
            test_temp <- NULL
            trtst_temp <- NULL
            for (i in 1:200) {
              weight_mat <- rbind(weight_mat, c(as.numeric(res[[i]]$u_weight)))
              
              pvalue_temp <- c(pvalue_temp, res[[i]]$pvalue)
              pvalue2_temp <- c(pvalue2_temp, res[[i]]$pvalue2)
              lam_temp <- c(lam_temp, res[[i]]$lam)
              K_temp <- c(K_temp, res[[i]]$K_tr)
              A_temp <- c(A_temp, res[[i]]$A_tr)
              V0_invlog_temp <- c(V0_invlog_temp, log(as.numeric(res[[i]]$V0_inv_tr)))
              train_temp <- c(train_temp, res[[i]]$train_RMSE)
              test_temp <- c(test_temp, res[[i]]$test_RMSE)
              trtst_temp <- c(trtst_temp, res[[i]]$trtst_ratio)
            }
            
            weight_sum <- rbind(weight_sum, c(colMeans(weight_mat), 
                                              int_effect, mode, strategy))
            
            sim_power <- sum(pvalue_temp < .05) / 200
            sim2_power <- sum(pvalue2_temp < .05) / 200
            lam_mean <- round(median(lam_temp), 4)
            K_mean <- round(mean(K_temp), 4)
            A_mean <- round(mean(A_temp), 4)
            V0_invlog_mean <- round(mean(V0_invlog_temp), 4)
            train_mean <- round(mean(train_temp), 4)
            test_mean <- round(mean(test_temp), 4)
            trtst_mean <- round(mean(trtst_temp), 4)
            
            sum_stat <- rbind(sum_stat, c(as.character(method), 
                                          2 * l, l, int_effect, 
                                          as.character(mode), 
                                          as.character(strategy), 
                                          d, sim_power, sim2_power,
                                          lam_mean, K_mean, 
                                          A_mean, V0_invlog_mean, 
                                          train_mean, test_mean, 
                                          trtst_mean))

          }
        }
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
      weight_melt$int_effect <- as.numeric(as.character(weight_melt$int_effect))
      
      fig <- ggplot(data = weight_melt, aes(x = int_effect, y = weight, color = ind))+
        facet_grid(strategy ~ mode)+ geom_line()+ geom_point(size = 0.6) + 
        labs(x = 'int_effect', y = 'weight', title = paste0("library_", d)) +
        theme_set(theme_bw()) + theme(panel.grid = element_blank()) +
        theme(axis.title.x = element_text(size = 15), axis.title.y = element_text(size = 15),
              plot.title = element_text(size = 16), axis.text.x=element_text(size = 12),
              axis.text.y=element_text(size = 12), panel.grid =element_blank()) +
        theme(legend.title = element_text(size = 12, face = "bold")) +
        theme(legend.text = element_text(size = 12))
      ggsave(paste0("./plots/library_", d, "_true", as.character(method), "_p", 2 * l, 
                    "_l", l, "_weight.pdf"),
             fig, width = 14, height = 9)
      
    }
  }
}

for (method in "matern") {
  for (p in 1:2) {
    for (l in c(.5, 1, 1.5)) {
      for (d in 1) {
        weight_sum <- NULL
        for (mode in c("loocv", "AICc", "GCVc", "gmpml")) {
          for (strategy in c("avg")) {
            for (int_effect in c(seq(0, .29, .05), seq(.3, 1, .1))) {
              
              res_key <-
                paste0("true", method, 
                       "_p", p, 
                       "_l", l,
                       "_mode", mode, 
                       "_strategy", strategy, 
                       "_int", int_effect,
                       "_K", d,
                       ".RData")
              
              load(res_key)
              weight_mat <- NULL
              
              pvalue_temp <- NULL
              pvalue2_temp <- NULL
              lam_temp <- NULL
              K_temp <- NULL
              A_temp <- NULL
              V0_invlog_temp <- NULL
              train_temp <- NULL
              test_temp <- NULL
              trtst_temp <- NULL
              for (i in 1:200) {
                weight_mat <- rbind(weight_mat, c(as.numeric(res[[i]]$u_weight)))
                
                pvalue_temp <- c(pvalue_temp, res[[i]]$pvalue)
                pvalue2_temp <- c(pvalue2_temp, res[[i]]$pvalue2)
                lam_temp <- c(lam_temp, res[[i]]$lam)
                K_temp <- c(K_temp, res[[i]]$K_tr)
                A_temp <- c(A_temp, res[[i]]$A_tr)
                V0_invlog_temp <- c(V0_invlog_temp, log(as.numeric(res[[i]]$V0_inv_tr)))
                train_temp <- c(train_temp, res[[i]]$train_RMSE)
                test_temp <- c(test_temp, res[[i]]$test_RMSE)
                trtst_temp <- c(trtst_temp, res[[i]]$trtst_ratio)
              }
              
              weight_sum <- rbind(weight_sum, c(colMeans(weight_mat), 
                                                int_effect, mode, strategy))
              
              sim_power <- sum(pvalue_temp < .05) / 200
              sim2_power <- sum(pvalue2_temp < .05) / 200
              lam_mean <- round(median(lam_temp), 4)
              K_mean <- round(mean(K_temp), 4)
              A_mean <- round(mean(A_temp), 4)
              V0_invlog_mean <- round(mean(V0_invlog_temp), 4)
              train_mean <- round(mean(train_temp), 4)
              test_mean <- round(mean(test_temp), 4)
              trtst_mean <- round(mean(trtst_temp), 4)
              
              sum_stat <- rbind(sum_stat, c(as.character(method), 
                                            p, l, int_effect, 
                                            as.character(mode), 
                                            as.character(strategy), 
                                            d, sim_power, sim2_power, 
                                            lam_mean, K_mean, 
                                            A_mean, V0_invlog_mean, 
                                            train_mean, test_mean, 
                                            trtst_mean))
            }
          }
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
        weight_melt$int_effect <- as.numeric(as.character(weight_melt$int_effect))
        
        fig <- ggplot(data = weight_melt, aes(x = int_effect, y = weight, color = ind))+
          facet_grid(strategy ~ mode)+ geom_line()+ geom_point(size = 0.6) + 
          labs(x = 'int_effect', y = 'weight', title = paste0("library_", d)) +
          theme_set(theme_bw()) + theme(panel.grid = element_blank()) +
          theme(axis.title.x = element_text(size = 15), axis.title.y = element_text(size = 15),
                plot.title = element_text(size = 16), axis.text.x=element_text(size = 12),
                axis.text.y=element_text(size = 12), panel.grid =element_blank()) +
          theme(legend.title = element_text(size = 12, face = "bold")) +
          theme(legend.text = element_text(size = 12))
        ggsave(paste0("./plots/library_", d, "_true", as.character(method), "_p", p, 
                      "_l", l, "_weight.pdf"),
               fig, width = 14, height = 9)
        
      }
    }
  }
}


for (method in "matern") {
  for (p in 1:2) {
    for (l in c(.5, 1, 1.5)) {
      for (d in 2:5) {
        weight_sum <- NULL
        for (mode in c("loocv", "AICc", "GCVc", "gmpml")) {
          for (strategy in c("avg", "exp", "erm")) {
            for (int_effect in c(seq(0, .29, .05), seq(.3, 1, .1))) {
              
              res_key <-
                paste0("true", method, 
                       "_p", p, 
                       "_l", l,
                       "_mode", mode, 
                       "_strategy", strategy, 
                       "_int", int_effect,
                       "_K", d,
                       ".RData")
              
              load(res_key)
              weight_mat <- NULL
              
              pvalue_temp <- NULL
              pvalue2_temp <- NULL
              lam_temp <- NULL
              K_temp <- NULL
              A_temp <- NULL
              V0_invlog_temp <- NULL
              train_temp <- NULL
              test_temp <- NULL
              trtst_temp <- NULL
              for (i in 1:200) {
                weight_mat <- rbind(weight_mat, c(as.numeric(res[[i]]$u_weight)))
                
                pvalue_temp <- c(pvalue_temp, res[[i]]$pvalue)
                pvalue2_temp <- c(pvalue2_temp, res[[i]]$pvalue2)
                lam_temp <- c(lam_temp, res[[i]]$lam)
                K_temp <- c(K_temp, res[[i]]$K_tr)
                A_temp <- c(A_temp, res[[i]]$A_tr)
                V0_invlog_temp <- c(V0_invlog_temp, log(as.numeric(res[[i]]$V0_inv_tr)))
                train_temp <- c(train_temp, res[[i]]$train_RMSE)
                test_temp <- c(test_temp, res[[i]]$test_RMSE)
                trtst_temp <- c(trtst_temp, res[[i]]$trtst_ratio)
              }
              
              weight_sum <- rbind(weight_sum, c(colMeans(weight_mat), 
                                                int_effect, mode, strategy))
              
              sim_power <- sum(pvalue_temp < .05) / 200
              sim2_power <- sum(pvalue2_temp < .05) / 200
              lam_mean <- round(median(lam_temp), 4)
              K_mean <- round(mean(K_temp), 4)
              A_mean <- round(mean(A_temp), 4)
              V0_invlog_mean <- round(mean(V0_invlog_temp), 4)
              train_mean <- round(mean(train_temp), 4)
              test_mean <- round(mean(test_temp), 4)
              trtst_mean <- round(mean(trtst_temp), 4)
              
              sum_stat <- rbind(sum_stat, c(as.character(method), 
                                            p, l, int_effect, 
                                            as.character(mode), 
                                            as.character(strategy), 
                                            d, sim_power, sim2_power, 
                                            lam_mean, K_mean, 
                                            A_mean, V0_invlog_mean, 
                                            train_mean, test_mean, 
                                            trtst_mean))
            }
          }
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
        weight_melt$int_effect <- as.numeric(as.character(weight_melt$int_effect))
        
        fig <- ggplot(data = weight_melt, aes(x = int_effect, y = weight, color = ind))+
          facet_grid(strategy ~ mode)+ geom_line()+ geom_point(size = 0.6) + 
          labs(x = 'int_effect', y = 'weight', title = paste0("library_", d)) +
          theme_set(theme_bw()) + theme(panel.grid = element_blank()) +
          theme(axis.title.x = element_text(size = 15), axis.title.y = element_text(size = 15),
                plot.title = element_text(size = 16), axis.text.x=element_text(size = 12),
                axis.text.y=element_text(size = 12), panel.grid =element_blank()) +
          theme(legend.title = element_text(size = 12, face = "bold")) +
          theme(legend.text = element_text(size = 12))
        ggsave(paste0("./plots/library_", d, "_true", as.character(method), "_p", p, 
                      "_l", l, "_weight.pdf"),
               fig, width = 14, height = 9)
        
      }
    }
  }
}

write.csv(sum_stat, file = "./plots/sim_result.csv", row.names = F, quote = F)
