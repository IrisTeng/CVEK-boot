
simWorker <- function(M = 200, n = 100, B = 200,
                      method = "rbf", int_effect = 0, 
                      Sigma = 0, l = 1, p = 2, eps = .01, 
                      mode = "loocv", strategy = "erm", 
                      beta = 1, test = "boot", d = 1,
                      lambda = exp(seq(-10, 5, .5))) {
  
  K1 <- data.frame(method = as.character(method),
                   Sigma = 0, l = l, p = p)
  K1$method <- as.character(K1$method)
  
  K2 <- data.frame(method = c("polynomial",
                              "polynomial", "polynomial"), 
                   Sigma = rep(0, 3), l = c(.5, 1, 1.5), p = 1:3)
  K2$method <- as.character(K2$method)
  
  K3 <- data.frame(method = c("rbf", "rbf", "rbf"),
                   Sigma = rep(0, 3), l = c(.6, 1, 2), p = 1:3)
  K3$method <- as.character(K3$method)
  
  K4 <- data.frame(method = c("polynomial",
                              "polynomial", "polynomial", 
                              "rbf", "rbf", "rbf"), 
                   Sigma = rep(0, 6), l = c(1, 1, 1, .6, 1, 2), p = 1:6)
  K4$method <- as.character(K4$method)
  
  K5 <- data.frame(method = c("matern",
                              "matern", "matern", 
                              "rbf", "rbf", "rbf"), 
                   Sigma = rep(0, 6), l = c(1, 1, 1, .6, 1, 2),
                   p = c(0, 1, 2, 2, 2, 2))
  K5$method <- as.character(K5$method)
  
  Kerns <- list()
  Kerns[[1]] <- K1
  Kerns[[2]] <- K2
  Kerns[[3]] <- K3
  Kerns[[4]] <- K4
  Kerns[[5]] <- K5
  
  kern <- Kerns[[d]]
  
  res <- vector("list", length = M)
  filename <-
    paste0("true", method, 
           "_p", p,
           "_l", l,
           "_mode", mode, 
           "_strategy", strategy, 
           "_int", int_effect,
           "_K", d,
           ".RData")
  header <- 
    paste0("true=", method, 
           "_p", p, 
           "_l", l,
           "mode=", mode, 
           ", strategy=", strategy, 
           ", int=", int_effect, 
           ", K=", d, ":")
  
  for (i in 1:M){
    res[[i]] <- 
      tryCatch(
        data.sim <-
          generate_data(n, label_names, method = as.character(method),
                        int_effect = int_effect, l = l, 
                        p = p, eps = eps),
        error = function(e) e
      )
    
    if ("error" %in% class(res[[i]])){
      print(paste0(header, "Iteration ", i, "/", M,
                   " | ", res[[i]]$message , " >=<"))
    } else {
      fit <- define_model(formula, label_names, data.sim, kern)
      
      #### 2. Fitting ####
      res[[i]] <- tryCatch(
        res_tst <- 
          testing(formula_int, label_names, 
                  fit$Y, fit$X1, fit$X2, 
                  fit$kern_list,  as.character(mode), as.character(strategy), 
                  beta, as.character(test), lambda, B),
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

