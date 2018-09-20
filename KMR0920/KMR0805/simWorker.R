source('./source/defineModel.R')
source('./source/ensemble.R')
source('./source/estimation.R')
source('./source/kernelGenerate.R')
source('./source/testing.R')
source('./source/tuning.R')
source('./source/util.R')

simWorker <- function(M = 200, n = 100, B = 100, 
                      method = "rbf", int_effect = 0, 
                      Sigma = 0, l = 1, p = 2, eps = .01, 
                      mode = "loocv", strategy = "erm", 
                      beta = 1, test = "boot", d = 1, 
                      lambda = exp(seq(-5, 5))){
  K1 <- NULL
  K1 <- c(K1, kernelGenerate(method = as.character(method), l = l, p = p))
  
  K2 <- NULL
  K2 <- c(K2, kernelGenerate('polynomial', p = 1))
  K2 <- c(K2, kernelGenerate('polynomial', p = 2))
  K2 <- c(K2, kernelGenerate('polynomial', p = 3))
  
  K3 <- NULL
  K3 <- c(K3, kernelGenerate('rbf', l = .6))
  K3 <- c(K3, kernelGenerate('rbf', l = 1))
  K3 <- c(K3, kernelGenerate('rbf', l = 2))
  
  K4 <- NULL
  K4 <- c(K4, kernelGenerate('polynomial', p = 1))
  K4 <- c(K4, kernelGenerate('polynomial', p = 2))
  K4 <- c(K4, kernelGenerate('polynomial', p = 3))
  K4 <- c(K4, kernelGenerate('rbf', l = .6))
  K4 <- c(K4, kernelGenerate('rbf', l = 1))
  K4 <- c(K4, kernelGenerate('rbf', l = 2))
  
  K5 <- NULL
  K5 <- c(K5, kernelGenerate('matern', p = 0))
  K5 <- c(K5, kernelGenerate('matern', p = 1))
  K5 <- c(K5, kernelGenerate('matern', p = 2))
  K5 <- c(K5, kernelGenerate('rbf', l = .6))
  K5 <- c(K5, kernelGenerate('rbf', l = 1))
  K5 <- c(K5, kernelGenerate('rbf', l = 2))
  
  Kerns <- list()
  Kerns[[1]] <- K1
  Kerns[[2]] <- K2
  Kerns[[3]] <- K3
  Kerns[[4]] <- K4
  Kerns[[5]] <- K5
  
  kern <- Kerns[[d]]
  
  res <- vector("list", length = M)
  filename <-
    paste0("mode", mode, 
           "_strategy", strategy, 
           "_int", int_effect,
           "_K", d,
           ".RData")
  header <- 
    paste0("mode=", mode, 
           ", strategy=", strategy, 
           ", int=", int_effect, 
           ", K=", d, ":")
  
  for (i in 1:M){
    res[[i]] <- 
      tryCatch(
        data.sim <-
          dataGenerate(n, label_names, method = method,
                       int_effect = int_effect, l = l, 
                       p = p, eps = eps),
        error = function(e) e
      )
    
    if ("error" %in% class(res[[i]])){
      print(paste0(header, "Iteration ", i, "/", M,
                   " | ", res[[i]]$message , " >_<"))
    } else {
      fit <- defineModel(formula, label_names, data.sim, kern)
      
      
      #### 2. Fitting ####
      res[[i]] <- tryCatch(
        res_tst <- 
          testing(formula_int, label_names, 
                  fit$Y, fit$X1, fit$X2, 
                  kern,  mode, strategy, 
                  beta, test, lambda, B),
        error = function(e) e
      )
      
      if(all(class(res[[i]]) == "list")){
        print(paste0(header, "Iteration ", i, "/", M, 
                     ", pval = ", 
                     paste(sapply(res[[i]]$pvalue, 
                                  function(p) round(p, 3)), 
                           collapse = "/")
        )
        )
        
      } else {
        print(paste0(header, "Iteration ", i, "/", M,
                     " | ", res[[i]]$message , " >_<"))
      }
    }
    
    flush.console()
    
    if (!dir.exists("./output/")){
      dir.create("./output/")
      dir.create("./output/Temp/")
    } else if (!dir.exists("./output/Temp/")){
      dir.create("./output/Temp/")
    }
  }
  save(res, file = paste0("./output/Temp/", filename))
  
  res2 <- NULL
  for (i in 1:M) {
    res2 <- c(res2, res[[i]]$pvalue)
  }
  sim_power <- sum(res2 < .05) / M
  cat(c(as.character(method), p, l, int_effect, as.character(mode), 
        as.character(strategy), n, as.character(beta), d, sim_power),
      file = "sim_power_boot.txt", append = T, "\n")
  
  return(filename)
}

