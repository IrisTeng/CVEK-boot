
simWorker <- function(M = 200, n = 100, B = 100, 
                      method = "rbf", int_effect = 0, 
                      Sigma = 0, l = 1, p = 2, eps = .01, 
                      mode = "loocv", strategy = "erm", 
                      beta = 1, test = "boot", d = 1, 
                      lambda = exp(seq(-5, 5))) {
  
  K1 <- data.frame(method = c("intercept", as.character(method)), 
                   Sigma = rep(0, 2), l = c(1, l), p = c(2, p))
  K1$method <- as.character(K1$method)
  
  K2 <- data.frame(method = c("intercept", "polynomial", 
                              "polynomial", "polynomial"), 
                   Sigma = rep(0, 4), l = c(1, .5, 1, 1.5), p = 0:3)
  K2$method <- as.character(K2$method)
  
  K3 <- data.frame(method = c("intercept", "rbf", "rbf", "rbf"), 
                   Sigma = rep(0, 4), l = c(1, .6, 1, 2), p = 0:3)
  K3$method <- as.character(K3$method)
  
  K4 <- data.frame(method = c("intercept", "polynomial", 
                              "polynomial", "polynomial", 
                              "rbf", "rbf", "rbf"), 
                   Sigma = rep(0, 7), l = c(1, 1, 1, 1, .6, 1, 2), p = 0:6)
  K4$method <- as.character(K4$method)
  
  K5 <- data.frame(method = c("intercept", "matern", 
                              "matern", "matern", 
                              "rbf", "rbf", "rbf"), 
                   Sigma = rep(0, 7), l = c(1, 1, 1, 1, .6, 1, 2), 
                   p = c(1, 0, 1, 2, 2, 2, 2))
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
        data_all <-
          generate_data(200, label_names, method = as.character(method),
                        int_effect = int_effect, l = l, 
                        p = p, eps = eps),
        error = function(e) e
      )
    
    if ("error" %in% class(res[[i]])){
      print(paste0(header, "Iteration ", i, "/", M,
                   " | ", res[[i]]$message , " >=<"))
    } else {
      data.sim <- data_all[1:100, ]
      data_test <- data_all[101:200, ]
      fit <- define_model(formula, label_names, data.sim, kern)
      fit_test <- define_model(formula, label_names, data_test, kern)
      
      #### 2. Fitting ####
      res[[i]] <- tryCatch(
        res_tst <- 
          testing(formula_int, label_names, 
                  fit$Y, fit$X1, fit$X2, 
                  fit$kern_list,  as.character(mode), as.character(strategy), 
                  beta, as.character(test), lambda, B, data_test, fit_test),
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
  
  pvalue_temp <- NULL
  lam_temp <- NULL
  K_temp <- NULL
  A_temp <- NULL
  V0_invlog_temp <- NULL
  trtst_temp <- NULL
  for (i in 1:M) {
    pvalue_temp <- c(pvalue_temp, res[[i]]$pvalue)
    lam_temp <- c(lam_temp, res[[i]]$lam)
    K_temp <- c(K_temp, res[[i]]$K_tr)
    A_temp <- c(A_temp, res[[i]]$A_tr)
    V0_invlog_temp <- c(V0_invlog_temp, log(res[[i]]$V0_inv_tr))
    trtst_temp <- c(trtst_temp, res[[i]]$trtst_ratio)
  }
  
  sim_power <- sum(pvalue_temp < .05) / M
  lam_mean <- round(mean(lam_temp), 4)
  K_mean <- round(mean(K_temp), 4)
  A_mean <- round(mean(A_temp), 4)
  V0_invlog_mean <- round(mean(V0_invlog_temp), 4)
  trtst_mean <- round(mean(trtst_temp), 4)
  cat(c(as.character(method), p, l, int_effect, as.character(mode), 
        as.character(strategy), n, as.character(beta), d, sim_power, 
        lam_mean, K_mean, A_mean, V0_invlog_mean, trtst_mean),
      file = "sim_power_boot.txt", append = T, "\n")
  
  return(filename)
}

