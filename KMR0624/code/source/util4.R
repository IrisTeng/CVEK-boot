GenericFormula <- 
  function(formula, label.names){
    # transform the format of predictors from vectors to single variables.
    #
    # Args:
    #   formula:  A symbolic description of the model to be fitted.
    #   label.names: A character string indicating all the interior variables
    #                included in each predictor.
    #
    # Returns:
    #   generic.formula: The specific formula of the model to be fitted.
    #   length.main.effect: A numeric number indicating the number 
    #                       of main effect predictors.
    
    formula.factors <- attr(terms(formula), "factors")
    generic.formula <- Y ~ 1
    length.main.effect <- 0
    for (i in 1:dim(formula.factors)[2]){
      terms.names <- 
        rownames(formula.factors)[which(formula.factors[, i] == 1)]
      if (length(terms.names) == 1){
        generic.formula <- 
          update.formula(generic.formula,
                         as.formula(paste(as.character(
                           attr(terms(formula), "variables"))[2],
                           paste(label.names[[terms.names]], collapse=" + "), sep=" ~ .+")))
        length.main.effect <- length.main.effect + length(label.names[[terms.names]])
      } else {
        interaction.formula <- 
          paste("(", paste(label.names[[terms.names[1]]], collapse=" + "), ")", sep="")
        for (j in 2:length(terms.names)){
          interaction.formula <- 
            paste(interaction.formula, "*(",
                  paste(label.names[[terms.names[j]]], collapse=" + "), ")", sep=" ")
        }
        generic.formula <- 
          update.formula(generic.formula,
                         as.formula(paste(as.character(
                           attr(terms(formula), "variables"))[2], interaction.formula, sep=" ~ .+")
                         )
          )
      }
    }
    return(list(generic.formula =  generic.formula,
                length.main.effect = length.main.effect))
  }


OriginalData2 <- 
  function(size, label.names, method = NULL, int.effect = 0,
           l = l, p = p, eps = 0.01){
    # Generate original data.
    #
    # Args:
    #   size: A numeric number specifying the number of observations.
    #   label.names: A character string indicating all the interior variables
    #                included in each predictor.
    #   method: A character string indicating which Kernel is to be computed.
    #   int.effect: A numeric number specifying the size of interaction.
    #   eps: A numeric number indicating the size of noise.
    #
    # Returns:
    #   data: A dataframe including response and predictors.
    
    # X~N(0,I)
      X1 <- rmvnorm(n = size,
                    mean = rep(0, length(label.names[[1]])),
                    sigma = diag(length(label.names[[1]])))
      X2 <- rmvnorm(n = size,
                    mean = rep(0, length(label.names[[2]])),
                    sigma = diag(length(label.names[[2]])))
      
      # X.int <- NULL
      # for (i in 1:length(label.names[[1]])){
      #   X.int <- cbind(X.int, X1[, i] * X2)
      # }
      Kern <- KernelGenerate(method, l = l, p = p)
      w1 <- rnorm(size)
      # w2 <- rnorm(size)
      w2 <- w1
      w12 <- rnorm(size)
      K1 <- Kern(X1, X1)
      K2 <- Kern(X2, X2)
      K1 <- K1 / tr(K1)
      K2 <- K2 / tr(K2)
      h0 <- K1 %*% w1 + K2 %*% w2
      h0 <- h0 / sqrt(sum(h0 ^ 2))

      h1.prime <- (K1 * K2) %*% w12
      Ks <- svd(K1 + K2)
      if(length(Ks$d / sum(Ks$d) > 0.001) > 0){
        len <- length(Ks$d[Ks$d / sum(Ks$d) > 0.001])
        U0 <- Ks$u[, 1:len]
        h1.prime.hat <- fitted(lm(h1.prime ~ U0))
        h1 <- h1.prime - h1.prime.hat
        if(all(h1 == 0)){
          warning("interaction term colinear with main-effect space!")
          h1 <- h1.prime
          h1 <- h1 / sqrt(sum(h1 ^ 2))
        }
        else
          h1 <- h1 / sqrt(sum(h1 ^ 2))
      }
      else{
        warning("largest eigen value smaller than 0.001!")
        h1 <- h1.prime
        h1 <- h1 / sqrt(sum(h1 ^ 2))
      }
      # h1 <- (K1 * K2) %*% w12
      
      Y <- h0 + int.effect * h1 + rnorm(1) + rnorm(size, 0, eps)
      ## run again
      # Y <- (Y - mean(Y)) / range(Y)
      data <- as.data.frame(cbind(Y, X1, X2))
      colnames(data) <- c("Y", label.names[[1]], label.names[[2]])
      return(data)
  }


NoiseEstimate <- function(y, lambda.hat, beta.hat, alpha.hat, K.hat){
  # An implementation of Gaussian processes for estimating noise.
  #
  # Args:
  #   y: A vector of response from original data.
  #   lambda.hat: Lambda calculated from the original data.
  #   beta.hat: Intercept calculated from the original data.
  #   alpha.hat: A vector of parameters calculated from the original data.
  #   K.hat: Kernel matrix calculated from the original data.

  #
  # Returns:
  #   sigma2.hat: A numeric number indicating the estimated noise.
  
  n <- nrow(K.hat)
  V <- lambda.hat * diag(n) + K.hat
  one <- rep(1, n)
  P.x <- one %*% ginv(t(one) %*% ginv(V) %*% one) %*% t(one) %*% ginv(V)
  P.k <- K.hat %*% ginv(V) %*% (diag(n) - P.x)
  A <- P.x + P.k
  sigma2.hat <- sum((y - beta.hat - K.hat %*% alpha.hat) ^ 2) / (n - tr(A))
  return(sigma2.hat)
 }


MultiScore2 <- 
  function(formula, label.names, 
           data = NULL, l = 1, lambda = exp(seq(-5, 5, 1))){
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
    #
    # Returns:
    #   test.stat: A numeric number indicating the score test statistics.
    
    # prepare data and estimate under null
    y <- data[, as.character(attr(terms(formula), "variables"))[2]]
    re <- GenericFormula(formula, label.names)
    generic.formula0 <- re[[1]]
    len <- re[[2]]
    X <- model.matrix(generic.formula0, data)[, -1]
    n <- nrow(X)
    X1 <- X[, c(1:length(label.names[[1]]))]
    X2 <- X[, c((length(label.names[[1]]) + 1):len)]
    X12 <- X[, c((len + 1):dim(X)[2])]
    
    # X12 <- NULL
    # for (i in 1:length(label.names[[1]])){
    #   X12 <- cbind(X12, X1[, i] * X2)
    # }
    
    # result <- gpr(Y ~ X1 + X2, label.names, data, method, l, lambda)
    result <- Ensemble(Y ~ X1 + X2, label.names, Kernlist, beta = 1,
                       data, lambda, mode = mode, strategy = strategy)
    lam <- result[[1]]
    beta0 <- result[[2]]
    alpha0 <- result[[3]]
    K.gpr <- result[[4]]
    sigma2.hat <- NoiseEstimate(y, lam, beta0, alpha0, K.gpr)
    tau.hat <- sigma2.hat / lam
      
    # compute score statistic
    # Kern <- KernelGenerate(method, l)
    # K1 <- Kern(X1, X1)
    # K2 <- Kern(X2, X2)
    # K0 <- K1 + K2
    K0 <- K.gpr
    K12 <- X12 %*% t(X12)
    
    V0_inv <- ginv(tau.hat * K0 + sigma2.hat * diag(n))
    
    test.stat <- tau.hat * t(y - beta0) %*% V0_inv %*%
      K12 %*% V0_inv %*% (y - beta0) / 2
    
    return(test.stat)
  }


InfoMat <- 
  function(P0.mat,
           mat.del = NULL, mat.sigma2 = NULL, 
           mat.tau = NULL)
  {
    I0 <- matrix(NA, 3, 3)
    
    I0[1,1] <- tr( P0.mat %*% mat.del %*% P0.mat %*% mat.del )/2  ##del.del
    #
    I0[1,2] <- tr( P0.mat %*% mat.del %*% P0.mat %*% mat.sigma2 )/2  ##del.sig
    I0[2,1] <- I0[1,2]
    #
    I0[1,3] <- tr( P0.mat %*% mat.del %*% P0.mat %*% mat.tau )/2  ##del.lam
    I0[3,1] <- I0[1,3]
    #
    #
    I0[2,2] <- tr( P0.mat %*% mat.sigma2 %*% P0.mat %*% mat.sigma2 )/2  ##sig.sig
    #
    I0[2,3] <-  tr( P0.mat %*% mat.sigma2 %*% P0.mat %*% mat.tau )/2  ##sig.lam
    I0[3,2] <- I0[2,3]
    #
    I0[3,3] <-  tr( P0.mat %*% mat.tau %*% P0.mat %*% mat.tau )/2  ##lam.lam
    
    return(I0)
  }


KernelTest <- 
  function(formula, label.names, Kernlist, mode, strategy, test = "boot",
           data = NULL, l = 1, lambda = exp(seq(-5, 5, 1)), B = 100){
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
    generic.formula0 <- 
      GenericFormula(formula,label.names)$generic.formula
    
    # extract design matrix
    X <- model.matrix(generic.formula0, data)[, -1]
    n <- nrow(X)
    Xm <- colMeans(X)
    p <- ncol(X)
    X <- X - rep(Xm, rep(n, p))
    Xscale <- drop(rep(1 / n, n) %*% X ^ 2) ^ 0.5
    X <- X / rep(Xscale, rep(n, p))
    
    X1 <- X[, c(1:length(label.names[[1]]))]
    X2 <- X[, c((length(label.names[[1]]) + 1):
                  (length(label.names[[1]]) + length(label.names[[2]])))]
    data0 <- as.data.frame(cbind(Y, X))
    
    # result <- gpr(formula, label.names, data0, method, l, lambda)
    result <- Ensemble(formula, label.names, Kernlist, data0, beta = 1,
                       lambda, mode = mode, strategy = strategy)
    lam <- result[[1]]
    beta0 <- result[[2]]
    alpha0 <- result[[3]]
    K.gpr <- result[[4]]
	  u.weight <- result[[5]]
    sigma2.hat <- NoiseEstimate(Y, lam, beta0, alpha0, K.gpr)
    mean.Y <- K.gpr %*% alpha0 + beta0
    
    if (test == "boot"){
      # conduct bootstrap
      bs.test <- sapply(1:B, function(k){
        Y.star <- mean.Y + rnorm(n, sd = sqrt(sigma2.hat)) 
        dat <- cbind(Y.star, X)
        colnames(dat)[1] <- "Y"
        dat <- as.data.frame(dat)
        MultiScore2(Y ~ X1 + X2 + X1 * X2, label.names, 
                    data = dat, l, lambda)
      })
      
      # assemble test statistic
      original.test <- 
        MultiScore2(Y ~ X1 + X2 + X1 * X2, label.names, 
                    data = data0, l, lambda)
      
      bs.pvalue <- sum(as.numeric(original.test) <= bs.test) / B
      
      result <- list(pvalue = bs.pvalue, bs.test = bs.test, 
                     original.test = original.test, u.weight = u.weight)
    }
    else if (test == "asym"){
      score.chi <- 
        MultiScore2(Y ~ X1 + X2 + X1 * X2, label.names, 
                    data = data0, l, lambda)
      X12 <- NULL
      for (i in 1:length(label.names[[1]])){
        X12 <- cbind(X12, X1[, i] * X2)
      }
      
      tau.hat <- sigma2.hat / lam
      K0 <- K.gpr
      K12 <- X12 %*% t(X12)
      V0_inv <- ginv(tau.hat * K0 + sigma2.hat * diag(n))
      one <- rep(1, n)
      P0.mat <- V0_inv - V0_inv %*%
        one %*% ginv(t(one) %*% V0_inv %*% one) %*% t(one) %*% V0_inv
      
      drV0.tau <- K0
      drV0.sigma2 <- diag(n)
      drV0.del <- tau.hat * K12
      
      I0 <- InfoMat(P0.mat, 
                    mat.del = drV0.del, mat.sigma2 = drV0.sigma2, 
                    mat.tau = drV0.tau)
      
      #Effective Info for delta
      tot.dim <- ncol(I0)
      I.deldel <-  
        I0[1,1] - 
        I0[1,2:tot.dim] %*% ginv(I0[2:tot.dim,2:tot.dim]) %*% I0[2:tot.dim, 1] 
      
      #
      md <- tau.hat * tr(K12 %*% P0.mat) / 2
      
      m.chi <- I.deldel / (2 * md)
      d.chi <- md / m.chi
      
      pval.chi <- 1 - pchisq(score.chi / m.chi, d.chi)
      
      result <- list(pvalue = pval.chi, stat = score.chi, 
                     m = m.chi, d = d.chi)
    }
    else
      stop("test must be boot or asym!")
    
    return(result)
  }

