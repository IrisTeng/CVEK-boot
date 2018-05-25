# plot
setwd('/Users/iristeng/Desktop/KMR/KMR0511/output')
library(ggplot2)
file <- read.table("simulation_power_2.txt")
colnames(file) <- c('int.effect', 'method', 'mode', 'strategy', 'n', 'power')
file$p <- rep(rep(1:3, rep(32, 3)), 3)
file$para <- c(rep(paste0("p=", 1:3), rep(32, 3)), 
               rep(paste0("l=", c(0.5, 1, 1.5)), rep(32, 3)), 
               rep(paste0("nu=", c(0.5, 1.5, 2.5), ", l=", c(0.5, 1, 1.5)), 
                   rep(32, 3)))

file$test <- ifelse(file$mode == 'loocv' & file$strategy == 'loocv', 1, 
                    ifelse(file$mode == 'loocv' & file$strategy == 'average', 2, 
                           ifelse(file$mode == 'AICc' & file$strategy == 'loocv', 3, 
                                  ifelse(file$mode == 'AICc' & file$strategy == 'average', 4, 
                                         ifelse(file$mode == 'GCVc' & file$strategy == 'loocv', 5, 
                                                ifelse(file$mode == 'GCVc' & file$strategy == 'average', 6, 
                                                       ifelse(file$mode == 'gmpml' & file$strategy == 'loocv', 7, 8)))))))

file$test <- factor(file$test, levels = 1:8,
                          labels = c('loocv, erm','loocv, average',
                                     'AICc, erm','AICc, average',
                                     'GCVc, erm','GCVc, average',
                                     'gmpml, erm','gmpml, average'))
file$p <- as.factor(file$p)
file$method <- as.factor(file$method)
ggplot(data = file, aes(x = int.effect, y = power, shape = test, color = test))+
  geom_point(size = 0.8) + geom_line() +
  scale_shape_manual(values = seq(0, 8)) + 
  geom_hline(yintercept = 0.05) + facet_grid(p ~ method)+
  labs(x = 'interaction', y = 'probability') +
  theme_set(theme_bw()) + theme(panel.grid = element_blank()) + 
  theme(axis.title.x = element_text(size = 15), axis.title.y = element_text(size = 15), 
        plot.title = element_text(size = 16), axis.text.x=element_text(size = 12), 
        axis.text.y=element_text(size = 12), panel.grid =element_blank()) + 
  theme(legend.title = element_text(size = 12, face = "bold")) +
  theme(legend.text = element_text(size = 12))

## 11:9