library(mvtnorm)
library(MASS)
library(snowfall)
library(psych)

source('./source/KernelGenerate.R')
source('./source/LooCV.R')
source('./source/gpr.R')
source('./source/util3.R')



verify <- 
  function(temp){
    if (!((simu.method == "primal") && 
          (generate.method == "gaussian"))){
      # generate data
      rawData <- OriginalData2(size = n, label.names, simu.method = simu.method,
                               method = generate.method, int.effect = int.effect)
      
      if (lambda.method == "No.Re"){
        # no regularization
        linear.res <- LinearBoost(formula, label.names, data = rawData,
                                  method = "linear", B = B)
        
        quadratic.res <- LinearBoost(formula, label.names, data = rawData, 
                                     method = "quadratic", B = B)	
        
        return(c(linear.pvalue = linear.res[[1]], 
                 quadratic.pvalue = quadratic.res[[1]],
                 linear.true.stat = linear.res[[3]],
                 quadratic.true.stat = quadratic.res[[3]],
                 linear.boost.stat = linear.res[[2]],
                 quadratic.boost.stat = quadratic.res[[2]])
        )
      }
      else if (lambda.method == "Loocv.Re"){
        # Loocv regularization
        linear.res <- KernelBoot(formula, label.names, data = rawData,
                                 method = "linear", l=l, B = B)
        
        gaussian.res <- KernelBoot(formula, label.names, data = rawData, 
                                   method = "gaussian", l=l, B = B)
        
        quadratic.res <- KernelBoot(formula, label.names, data = rawData, 
                                    method = "quadratic", l=l, B = B)	
        
        return(c(linear.pvalue = linear.res[[1]], 
                 gaussian.pvalue = gaussian.res[[1]], 
                 quadratic.pvalue = quadratic.res[[1]],
                 linear.true.stat = linear.res[[3]],
                 gaussian.true.stat = gaussian.res[[3]],
                 quadratic.true.stat = quadratic.res[[3]],
                 linear.boost.stat = linear.res[[2]],
                 gaussian.boost.stat = gaussian.res[[2]],
                 quadratic.boost.stat = quadratic.res[[2]])
        )
      }
      else
        stop("lambda method must be No or Loocv Regularization!")
    }
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
for(simu.method in c("dual", "primal")) {
  for(lambda.method in c("No.Re", "Loocv.Re")) {
    for(n in c(200)) {
      for (M in c(200)) {
        for (generate.method in c("linear", "gaussian", "quadratic")) {
          for (int.effect in seq(0, 0.5, 0.05)) {
            for (l in 1) {
              for (B in c(200)) {
                sfExport("formula", "label.names", "int.effect", "simu.method", 
                         "lambda.method","n", "B", "M", "l", "generate.method", 
                         "KernelBoot", "OriginalData2",
                         "GenericFormula", "gpr", "MultiScore2",
                         "LinearBoost", "LinearTestStat") 
                system.time(res <- sfSapply(1:M, verify))
                file.name <- paste(paste(simu.method, lambda.method, 
                                         n, generate.method, int.effect, sep = "_"), 
                                   ".txt", sep = "")
                write.table(t(res), file = file.name,
                            row.names = F, col.names = F)
                if (lambda.method == "No.Re"){
                  res2 <- apply(res[c(1, 2), ], 1, 
                                function(x){sum(x < 0.05) / M
                                })
                  cat(c(simu.method, lambda.method, 
                        n, generate.method, int.effect, res2),
                      file = "simulation_power.txt", append = T, "\n")
                }
                else{
                  res2 <- apply(res[c(1, 2, 3), ], 1, 
                                function(x){sum(x < 0.05) / M
                                })
                  cat(c(simu.method, lambda.method, 
                        n, generate.method, int.effect, res2),
                      file = "simulation_power.txt", append = T, "\n")
                }
                
                cat("Finished:interaction effect size=", int.effect, 
                    "n=", n, "l=", l, "\n")
              }
            }
          }
        }
      }
    }
  }
}

#write.csv(result, file = "simulation_power.csv", row.names = F, quote = F)

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