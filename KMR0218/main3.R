library(mvtnorm)
library(MASS)
library(snowfall)
library(psych)

source('./source/KernelGenerate.R')
source('./source/LooCV.R')
source('./source/gpr.R')
source('./source/util.R')



verify <- 
  function(temptemp){
    rawData <- data
    rawData[, 1] <- rawData[, 1]  + rnorm(n, 0, 0.01)
    res <- KernelBoot(formula, label.names, data = rawData, 
                      method = method, l=l, 
                      lambda = exp(seq(-5, 5, 1)), B = B)
    
    #gaussian.pvalue <- res[[1]]
    result.test <- res[[2]]
    bs.test <- result.test[1, ]
    perturb.test <- result.test[2, ]				  			
    
    original.test <-  res[[3]]
    
    return(c(original.test,bs.test,perturb.test))
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
for(n in c(50,100,200,300,400,500)){
  for (M in c(100)){
    for (method in c("quadratic")){
      for(int.effect in 0){
        for(l in 1){
          for (B in c(100)){
            l <- 1
            X1 <- rmvnorm(n,
                          mean = rep(0, length(label.names[[1]])),
                          sigma = diag(length(label.names[[1]])))
            X2 <- rmvnorm(n,
                          mean = rep(0, length(label.names[[2]])),
                          sigma = diag(length(label.names[[2]])))
            
            Kern <- KernelGenerate(method, l = 1)
            w1 <- rnorm(n)
            w2 <- w1
            K1 <- Kern(X1, X1)
            K2 <- Kern(X2, X2)
            h0 <- K1 %*% w1 + K2 %*% w2
            h0 <- h0 / sqrt(sum(h0 ^ 2))
            bias <- rnorm(1)
            Y <- h0 + bias + rnorm(n, 0, 0.01)
            data<- as.data.frame(cbind(Y, X1, X2))
            colnames(data) <- c("Y", label.names[[1]], label.names[[2]])

            sfExport("formula", "label.names", "int.effect",
                     "n", "B", "M", "l", "method", "KernelBoot", "OriginalData2", 
                     "GenericFormula", "gpr", "MultiScore2", "data") 
            system.time(res <- sfSapply(1 : M, verify))
            write.table(t(res), file = "quadratic_stat.txt",
                        row.names = F, col.names = F, append = T)            
            #res2 <- apply(res, 1, function(x){sum(x < 0.05) / M})
            #result <- rbind(result, c(int.effect, l, method, n, res2))
            #cat(c(int.effect, l, method, n, res2),
            #    file = "simulation_power.txt", append = T, "\n")
            cat("Finished:interaction effect size=", int.effect, 
                "n=", n, "l=", l, "\n")
          }
        }
      }
    }
  }
}
#write.csv(result, file = "simulation_power.csv", row.names = F, quote = F)
