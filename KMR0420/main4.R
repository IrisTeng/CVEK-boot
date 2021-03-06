setwd('/home/wd45/KMR0415')
library(mvtnorm)
library(MASS)
library(snowfall)
library(psych)
library(limSolve)
source('./source/KernelGenerate.R')
source('./source/Tuning.R')
source('./source/Ensemble.R')
source('./source/util3.R')


verify <- 
  function(temp){
    rawData <- OriginalData2(size = n, label.names, 
                             method = method, int.effect = int.effect)
    dd <- KernelBoot(formula, label.names, Kernlist, mode, strategy, rawData)
    
    pvalue <- dd[[1]]
    # bs.test <- dd[[2]]
    # original.test <- dd[[3]]
    return(pvalue)
  }



sfInit(parallel = T, cpus = 20)
sfLibrary(mvtnorm)
sfLibrary(MASS)
sfLibrary(psych)
sfLibrary(limSolve)
sfSource('./source/KernelGenerate.R')
sfSource('./source/Tuning.R')
sfSource('./source/Ensemble.R')
sfSource('./source/util3.R')
formula <- Y ~ X1 + X2
label.names <- list(X1 = c("x1", "x2"), X2 = c("x3", "x4"))

Kernlist <- NULL
# intercept kernel
Kernlist <- c(Kernlist, KernelGenerate('intercept'))
# linear kernel
Kernlist <- c(Kernlist, KernelGenerate('linear'))
# polynomial kernel, p = 2
Kernlist <- c(Kernlist, KernelGenerate('polynomial', p = 2))
# polynomial kernel, p = 3
Kernlist <- c(Kernlist, KernelGenerate('polynomial', p = 3))
# rbf kernel, l = .6
Kernlist <- c(Kernlist, KernelGenerate('rbf', l = .6))
# rbf kernel, l = 1
Kernlist <- c(Kernlist, KernelGenerate('rbf', l = 1))
# rbf kernel, l = 2
Kernlist <- c(Kernlist, KernelGenerate('rbf', l = 2))
# rbf kernel, l = 3
Kernlist <- c(Kernlist, KernelGenerate('rbf', l = 3))

result <- NULL
for(n in c(50, 100)){
  for (M in c(200)){
    for (method in c("linear", "rbf", "polynomial")){
      for (int.effect in c(0, 0.1, 0.2, 0.3)){
        for (B in c(100)){
          for(mode in c("loocv", "AICc", "GCVc", "gmpml")){
            for(strategy in c('loocv', 'average')){
              sfExport("formula", "label.names", "int.effect", "method", 
                       "mode", "strategy", "n", "B", "M", "Kernlist") 
              system.time(res <- sfSapply(1:M, verify))
              # write.table(t(res), file = "simulation_dd.txt",
              #             row.names = F, col.names = F, append = T)  
              res2 <- sum(res < 0.05) / M
              # res2 <- sapply(res, function(x){sum(x < 0.05) / M})
              result <- rbind(result, c(n, method, int.effect, mode, strategy, res2))
              cat(c(int.effect, method, mode, strategy, n, res2),
                  file = "simulation_power.txt", append = T, "\n")
              cat("Finished:interaction effect size=", int.effect, 
                  "n=", n, "\n")
            }
          }
        }
      }
    }
  }
}


write.csv(result, file = "simulation_power.csv", row.names = F, quote = F)


