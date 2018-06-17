setwd('/home/wd45/KMR0617')
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
    rawData <- OriginalData2(size = n, label.names, l = l, p = p,
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

# 1st library: True Kernel
result <- NULL
for(n in c(100)){
  for (M in c(200)){
    for (method in c("matern")){
      for(p1 in 2)
        for (int.effect in seq(0, 1, .1)){
          for (B in c(100)){
            for(mode in c("loocv")){
              for(strategy in c('exp', 'average')){
                if(method == "matern"){
                  p <- p1 - 1
                  for(l in c(0.5, 1, 1.5)){
                    Kernlist <- list(KernelGenerate(method = method, l = l, p = p))
                    sfExport("formula", "label.names", "int.effect", "method", 
                             "mode", "strategy", "n", "B", "M", "l", "p", 
                             "Kernlist") 
                    system.time(res <- sfSapply(1:M, verify))
                    res2 <- sum(res < 0.05) / M
                    # res2 <- sapply(res, function(x){sum(x < 0.05) / M})
                    result <- rbind(result, c(n, method, p, l, int.effect,
                                              mode, strategy, res2))
                    cat(c(method, p, l, int.effect, mode, strategy, n, res2),
                        file = "simulation_power3.txt", append = T, "\n")
                    cat("Finished:interaction effect size=", int.effect, 
                        "n=", n, "\n")
                  }
                }
                else{
                  l <- p1 / 2
                  p <- p1
                  
                  Kernlist <- list(KernelGenerate(method = method, l = l, p = p))
                  sfExport("formula", "label.names", "int.effect", "method", 
                           "mode", "strategy", "n", "B", "M", "l", "p", 
                           "Kernlist") 
                  system.time(res <- sfSapply(1:M, verify))
                  res2 <- sum(res < 0.05) / M
                  #res2 <- sum(res < 0.05) / M
                  # res2 <- sapply(res, function(x){sum(x < 0.05) / M})
                  result <- rbind(result, c(n, method, p, l, int.effect,
                                            mode, strategy, res2))
                  cat(c(method, p, l, int.effect, mode, strategy, n, res2),
                      file = "simulation_power3.txt", append = T, "\n")
                  cat("Finished:interaction effect size=", int.effect, 
                      "n=", n, "\n")
                }
                
              }
            }
          }
        }
    }
  }
}

write.csv(result, file = "simulation_power3.csv", row.names = F, quote = F)