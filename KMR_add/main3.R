library(mvtnorm)
library(MASS)
library(snowfall)
library(psych)

source('./source/KernelGenerate.R')
source('./source/LooCV.R')
source('./source/gpr.R')
source('./source/util.R')

# n <- 200
# formula <- Y ~ X1 + X2
# label.names <- list(X1 = c("x1", "x2"), X2 = c("x3", "x4"))
# rawData <- OriginalData2(size = n, label.names,
#                          method = "linear", int.effect = 0)
# res <- gpr(formula, label.names, rawData, "linear",
#            l = 3, lambda = exp(seq(-5, 5, 1)), null.hypo = TRUE)
# yhat <- res$intercept + res$K1 %*% res$alpha[1:n] + res$K2 %*% res$alpha[-(1:n)]
# 
# plot(rawData$Y, yhat)
# abline(a=0, b=1)
# # 
# linear.pvalue <-
#   KernelBoot(formula, label.names, data = rawData,
#              method = "linear", l=3,
#              lambda = exp(seq(-5, 5, 1)), B = 100)
# 
# gaussian.pvalue <-
#   KernelBoot(formula, label.names, data = rawData,
#              method = "gaussian", l=3,
#              lambda = exp(seq(-5, 5, 1)), B = 100)

verify <- 
  function(i){
    rawData <- OriginalData2(size = n, label.names, 
                             method = method, int.effect = int.effect)
    linear.pvalue <- KernelBoot(formula, label.names, data = rawData,
                                method = "linear", l=l,
                                lambda = exp(seq(-5, 5, 1)), B = B)
    gaussian.pvalue <- KernelBoot(formula, label.names, data = rawData, 
                                  method = "gaussian", l=l, 
                                  lambda = exp(seq(-5, 5, 1)), B = B)
    return(c(linear.pvalue,gaussian.pvalue))
  }

# produce plot similar to the CVEK paper. 
# (x-axis: interaction effect size, y-axis: ratio of p-value smaller than 0.05). 
# vary the interaction effect size between seq(0, 0.1, 0.02), 
# vary the effect size multiplication factor within  (1, 3, 5, 10), 
# vary the noise level within (0.05, 0.1, 0.25, 0.5, 1). 
# Keep bootstrap sample size at 200. Keep b11 = 2 (lower bound of effect size).
# number of repetition M fix at 1000.
# Repeat such experiment for linear and ridge regression.

sfInit(parallel = T, cpus = 20)
sfLibrary(mvtnorm)
sfLibrary(MASS)
sfLibrary(psych)
sfSource('./source/KernelGenerate.R')
sfSource('./source/LooCV.R')
sfSource('./source/gpr.R')
sfSource('./source/util.R')
formula <- Y ~ X1 + X2
label.names <- list(X1 = c("x1", "x2"), X2 = c("x3", "x4"))
result <- NULL
for(n in seq(50, 250, 50)){
  for (M in c(200)){
    for (method in c("linear", "gaussian")){
      for(int.effect in seq(0, 1, 0.1)){
        for(l in c(0.6, 1, 3)){
          for (B in c(100)){
            sfExport("formula", "label.names", "int.effect",
                     "n", "B", "M", "l", "method", "KernelBoot", "OriginalData2",
                     "GenericFormula", "gpr", "MultiScore2") 
            system.time(res <- sfSapply(1 : M, verify))
            write.table(t(res), file = "simulation_power_pvalue.txt",
                        row.names = F, col.names = F, append = T)            
            res2 <- apply(res, 1, function(x){sum(x < 0.05) / M})
            result <- rbind(result, c(int.effect, l, method, n, res2))
            cat(c(int.effect, l, method, n, res2),
                file = "simulation_power.txt", append = T, "\n")
            cat("Finished:interaction effect size=", int.effect, 
                "n=", n, "l=", l, "\n")
          }
        }
      }
    }
  }
}
write.csv(result, file = "simulation_power.csv", row.names = F, quote = F)

# plot
# file <- read.csv(file="/Users/iristeng/Desktop/simulation_power.csv")
# #colnames(file) <- c('beta_int','b11','b12','b21','b22','n','B','prob')
# file$effect_ind <- ifelse(file$b11==2,1,
#                           ifelse(file$b11==6,2,ifelse(file$b11==10,3,4)))
# library(ggplot2)
# file$effect_ind <- factor(file$effect_ind, levels = 1:4,
#                           labels = c('effect size1','effect size2',
#                                      'effect size3','effect size4'))
# file$e <- factor(file$e, labels = c('e=0.05','e=0.10','e=0.25','e=0.50','e=1.00'))
# file$n <- as.factor(file$n)
# p <- ggplot(data = file, aes(x=beta_int, y=linear, shape=n, color=n))+
#   geom_point(size=0.8)+geom_line()+
#   geom_hline(yintercept=0.05)+facet_grid(effect_ind~e)+
#   labs(x='beta interaction', y='probability')+
#   theme_set(theme_bw())+theme(panel.grid=element_blank())
# p
# 
# q <- ggplot(data = file, aes(x=beta_int, y=ridge, shape=n, color=n))+
#   geom_point(size=0.8)+geom_line()+
#   geom_hline(yintercept=0.05)+facet_grid(effect_ind~e)+
#   labs(x='beta interaction', y='probability')+
#   theme_set(theme_bw())+theme(panel.grid=element_blank())
# q