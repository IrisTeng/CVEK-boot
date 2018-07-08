kernelTest <- 
  function(formula, label_names, Kernlist, 
           data = NULL, mode = "loocv", strategy = "erm", beta = 1, 
           test = "boot", lambda = exp(seq(-5, 5)), B = 100){
    # Compute score tests comparing a fitted model and a more general alternative model.
    #
    # Args:
    #   formula:  A symbolic description of the model to be fitted.
    #   label.names: A character string indicating all the interior variables
    #                included in each predictor.
    #   data: A dataframe to be fitted.
    #   method: A character string indicating which Kernel is to be computed.
    #   l: A numeric number indicating the hyperparameter of a specific kernel.
    #   lambda: A numeric string specifying the range of noise to be chosen.
    #           The lower limit of lambda must be above 0.
    #   B: The number of bootstrap replicates.
    #
    # Returns:
    #   bs.pvalue: A number indicating the test P-value.
    
    # prepare data
    Y <- data[, as.character(attr(terms(formula), "variables"))[2]]
    
    # extract the information from the given formula to construct a true formula
    generic_formula0 <- 
      genericFormula(formula, label_names)$generic_formula
    
    # extract design matrix
    X <- model.matrix(generic_formula0, data)[, -1]
    n <- nrow(X)
    Xm <- colMeans(X)
    p <- ncol(X)
    X <- X - rep(Xm, rep(n, p))
    Xscale <- drop(rep(1 / n, n) %*% X ^ 2) ^ .5
    X <- X / rep(Xscale, rep(n, p))
    
    X1 <- X[, c(1:length(label_names[[1]]))]
    X2 <- X[, c((length(label_names[[1]]) + 1):
                  (length(label_names[[1]]) + length(label_names[[2]])))]
    data0 <- as.data.frame(cbind(Y, X))
    
    result <- ensemble(formula, label_names, Kernlist, data0, 
                       mode, strategy, beta, lambda)
    lam <- result[[1]]
    beta0 <- result[[2]]
    alpha0 <- result[[3]]
    K_gpr <- result[[4]]
    u_weight <- result[[5]]
    sigma2_hat <- noiseEstimate(Y, lam, beta0, alpha0, K_gpr)
    meanY <- K_gpr %*% alpha0 + beta0
    
    if (test == "boot"){
      # conduct bootstrap
      bs_test <- sapply(1:B, function(k){
        Ystar <- meanY + rnorm(n, sd = sqrt(sigma2_hat)) 
        dat <- cbind(Ystar, X)
        colnames(dat)[1] <- "Y"
        dat <- as.data.frame(dat)
        multiScore(Y ~ X1 + X2 + X1 * X2, label_names, Kernlist, dat, 
                   mode, strategy, beta, lambda)
      })
      
      # assemble test statistic
      original_test <- 
        multiScore(Y ~ X1 + X2 + X1 * X2, label_names, Kernlist, data0, 
                   mode, strategy, beta, lambda)
      
      pvalue <- sum(as.numeric(original_test) <= bs_test) / B
    }
    else if (test == "asym"){
      score_chi <- 
        multiScore(Y ~ X1 + X2 + X1 * X2, label_names, Kernlist, data0, 
                   mode, strategy, beta, lambda)
      X12 <- NULL
      for (i in 1:length(label_names[[1]])){
        X12 <- cbind(X12, X1[, i] * X2)
      }
      
      tau_hat <- sigma2_hat / lam
      K0 <- K_gpr
      K12 <- X12 %*% t(X12)
      V0_inv <- ginv(tau_hat * K0 + sigma2_hat * diag(n))
      one <- rep(1, n)
      P0_mat <- V0_inv - V0_inv %*%
        one %*% ginv(t(one) %*% V0_inv %*% one) %*% t(one) %*% V0_inv
      
      drV0_tau <- K0
      drV0_sigma2 <- diag(n)
      drV0_del <- tau_hat * K12
      
      I0 <- infoMat(P0_mat, 
                    mat_del = drV0_del, mat_sigma2 = drV0_sigma2, 
                    mat_tau = drV0_tau)
      
      #Effective Info for delta
      tot_dim <- ncol(I0)
      I_deldel <-  
        I0[1, 1] - 
        I0[1, 2:tot_dim] %*% ginv(I0[2:tot_dim, 2:tot_dim]) %*% I0[2:tot_dim, 1] 
      
      #
      md <- tau_hat * tr(K12 %*% P0_mat) / 2
      
      m_chi <- I_deldel / (2 * md)
      d_chi <- md / m_chi
      
      pvalue <- 1 - pchisq(score_chi / m_chi, d_chi)
    }
    else
      stop("test must be boot or asym!")
    
    return(pvalue)
  }
