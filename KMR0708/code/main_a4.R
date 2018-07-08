setwd('/home/wd45/KMR0708')
library(mvtnorm)
library(MASS)
library(snowfall)
library(psych)
library(limSolve)
source('./source/kernelGenerate.R')
source('./source/tuning.R')
source('./source/kernelTest.R')
source('./source/ensemble.R')
source('./source/util.R')


verify <- 
  function(temp){
    rawData <- dataGenerate(size = n, label_names, method = method, int_effect = int_effect, l = l, p = p)
    pvalue <- kernelTest(formula, label_names, Kernlist, rawData, mode, strategy, test = "asym")
    return(pvalue)
  }



sfInit(parallel = T, cpus = 20)
sfLibrary(mvtnorm)
sfLibrary(MASS)
sfLibrary(psych)
sfLibrary(limSolve)
sfSource('./source/kernelGenerate.R')
sfSource('./source/tuning.R')
sfSource('./source/ensemble.R')
sfSource('./source/kernelTest.R')
sfSource('./source/util.R')
formula <- Y ~ X1 + X2
label_names <- list(X1 = c("x1", "x2"), X2 = c("x3", "x4"))

# 4th library: 3 polynomials with different degree
Kernlist <- NULL
Kernlist <- c(Kernlist, kernelGenerate('polynomial', p = 1))
Kernlist <- c(Kernlist, kernelGenerate('polynomial', p = 2))
Kernlist <- c(Kernlist, kernelGenerate('polynomial', p = 3))
Kernlist <- c(Kernlist, kernelGenerate('rbf', l = .6))
Kernlist <- c(Kernlist, kernelGenerate('rbf', l = 1))
Kernlist <- c(Kernlist, kernelGenerate('rbf', l = 2))

result <- NULL
for(n in c(100)){
  for (M in c(200)){
    for (method in c("polynomial", "rbf", "matern")){
      for(p1 in 1:3)
        for (int_effect in seq(0, .3, .1)){
          for (B in c(100)){
            for(mode in c("loocv", "AICc", "GCVc", "gmpml")){
              for(strategy in c('exp')){
                if(method == "matern"){
                  p <- p1 - 1
                  for(l in c(0.5, 1, 1.5)){
                    # Kernlist <- list(KernelGenerate(method = method, l = l, p = p))
                    sfExport("formula", "label_names", "int_effect", "method", 
                             "mode", "strategy", "n", "B", "M", "l", "p", 
                             "Kernlist") 
                    system.time(res <- sfSapply(1:M, verify))
                    res2 <- sum(res < 0.05) / M
                    # res2 <- sapply(res, function(x){sum(x < 0.05) / M})
                    result <- rbind(result, c(n, method, p, l, int_effect,
                                              mode, strategy, res2))
                    cat(c(method, p, l, int_effect, mode, strategy, n, res2),
                        file = "simulation_power_a4.txt", append = T, "\n")
                    cat("Finished:interaction effect size=", int_effect, 
                        "n=", n, "\n")
                  }
                }
                else{
                  l <- p1 / 2
                  p <- p1
                  
                  # Kernlist <- list(KernelGenerate(method = method, l = l, p = p))
                  sfExport("formula", "label_names", "int_effect", "method", 
                           "mode", "strategy", "n", "B", "M", "l", "p", 
                           "Kernlist") 
                  system.time(res <- sfSapply(1:M, verify))
                  res2 <- sum(res < 0.05) / M
                  #res2 <- sum(res < 0.05) / M
                  # res2 <- sapply(res, function(x){sum(x < 0.05) / M})
                  result <- rbind(result, c(n, method, p, l, int_effect,
                                            mode, strategy, res2))
                  cat(c(method, p, l, int_effect, mode, strategy, n, res2),
                      file = "simulation_power_a4.txt", append = T, "\n")
                  cat("Finished:interaction effect size=", int_effect, 
                      "n=", n, "\n")
                }
                
              }
            }
          }
        }
    }
  }
}

write.csv(result, file = "simulation_power_a4.csv", row.names = F, quote = F)
