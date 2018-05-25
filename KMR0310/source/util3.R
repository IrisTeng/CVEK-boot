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
  function(size, label.names, simu.method = "dual",
           method = NULL, int.effect = 0, eps = 0.01){
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
    if (simu.method == "primal") {
      n <- size
      p1 <- length(label.names[[1]])
      p2 <- length(label.names[[2]])
      X <- matrix(rnorm(n = n * (p1 + p2)), nrow = n)
      # produce interaction term
      X1 <- X[, 1:p1]
      X2 <- X[, (1 + p1) : (p1 + p2)]
      
      X12 <- NULL
      for (i in 1:p1){
        X12 <- cbind(X12, X1[, i] * X2)
      }

      if (method == "linear"){
        X <- cbind(1, X)
        
        P.X <- X %*% ginv(t(X) %*% X) %*% t(X)
        
        beta <- rnorm(p1 + p2 + 1)
        beta12 <- rnorm(p1 * p2)
        
        h0 <- X %*% beta
        h12 <- X12 %*% beta12
        h0 <- h0 / sqrt(sum(h0 ^ 2))
        h12 <- h12 / sqrt(sum(h12 ^ 2))
        
        eps <- rnorm(n, sd = 0.01)
        
        Y <- h0 + int.effect * (diag(n) - P.X) %*% h12 + eps
        
        data <- as.data.frame(cbind(Y, X1, X2))
        colnames(data) <- c("Y", label.names[[1]], label.names[[2]])
        return(data)
      }
      else if(method == "quadratic") {
        X <- cbind(1, X, X ^ 2)
        
        P.X <- X %*% ginv(t(X) %*% X) %*% t(X)
        beta <- rnorm(2 * (p1 + p2) + 1)
        
        X12 <- cbind(X12, X12 ^ 2)
        beta12 <- rnorm(2 * (p1 * p2))
        
        h0 <- X %*% beta
        h12 <- X12 %*% beta12
        h0 <- h0 / sqrt(sum(h0 ^ 2))
        h12 <- h12 / sqrt(sum(h12 ^ 2))
        
        eps <- rnorm(n, sd = 0.01)
        
        Y <- h0 + int.effect * (diag(n) - P.X) %*% h12 + eps
        
        data <- as.data.frame(cbind(Y, X1, X2))
        colnames(data) <- c("Y", label.names[[1]], label.names[[2]])
        return(data)
      }
      else
        stop("In primal form, generation of data must be linear or quadratic!")
    } 
    else if(simu.method == "dual") {
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
      Kern <- KernelGenerate(method, l = 1)
      w1 <- rnorm(size)
      # w2 <- rnorm(size)
      w2 <- w1
      w12 <- rnorm(size)
      K1 <- Kern(X1, X1)
      K2 <- Kern(X2, X2)
      h0 <- K1 %*% w1 + K2 %*% w2
      h0 <- h0 / sqrt(sum(h0 ^ 2))
      
      h1.prime <- (K1 * K2) %*% w12
      Ks <- svd(Kern(X1, X1) + Kern(X2, X2))
      if(length(Ks$d / sum(Ks$d) > 0.001) > 0){
        len <- length(Ks$d[Ks$d / sum(Ks$d) > 0.001])
        # len <- length(Ks$d / sum(Ks$d) > 0.001)
        U0 <- Ks$u[, 1:len]
        h1.prime.hat <- fitted(lm(h1.prime ~ U0))
        h1 <- h1.prime - h1.prime.hat
        if(all(h1 == 0)){
          stop("interaction term colinear with main-effect space!")
          h1 <- h1.prime
          h1 <- h1 / sqrt(sum(h1 ^ 2))
        }
        else
          h1 <- h1 / sqrt(sum(h1 ^ 2))
      }
      else{
        stop("largest eigen value smaller than 0.001!")
        h1 <- h1.prime
        h1 <- h1 / sqrt(sum(h1 ^ 2))
      }
      
      Y <- h0 + int.effect * h1 + rnorm(1) + rnorm(size, 0, eps)
      data <- as.data.frame(cbind(Y, X1, X2))
      colnames(data) <- c("Y", label.names[[1]], label.names[[2]])
      return(data)
    } 
    else {
      stop("simulation method is wrong!")
    }
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
           data = NULL, method = NULL, l = 1, lambda = exp(seq(-5, 5, 1))){
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
    
    result <- gpr(Y ~ X1 + X2, label.names, data, method, l, lambda, null.hypo = T)
    sigma2.n <- result[[1]]
    beta0 <- result[[2]]
    alpha0 <- result[[3]]
    K.gpr <- result[[4]]
    noise.hat <- NoiseEstimate(y, sigma2.n, beta0, alpha0, K.gpr)
    tau.hat <- noise.hat / sigma2.n
      
    # compute score statistic
    Kern <- KernelGenerate(method, l)
    K1 <- Kern(X1, X1)
    K2 <- Kern(X2, X2)
    K0 <- K1 + K2
    K12 <- X12 %*% t(X12)
    V0 <- tau.hat * K0 + noise.hat * diag(n)
    
    test.stat <- tau.hat * t(y - beta0) %*% ginv(V0) %*% K12 %*% ginv(V0) %*% (y - beta0)
    
    # Jeremiah: heuristic procedure for test statistic correction, 
    # make sure to delete after adding estimation procedure for sigma
    # sigma_true <- 0.1
    # tau_true <- sigma_true/sigma2.n
    # test.stat <- test.stat/tau_true
    
    return(test.stat)
  }



KernelBoot <- 
  function(formula, label.names, data = NULL, 
           method = NULL, l = 1, lambda = exp(seq(-5, 5, 1)), B = 100){
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
    # Ym <- mean(Y)
    p <- ncol(X)
    X <- X - rep(Xm, rep(n, p))
    # Y <- Y - Ym
    Xscale <- drop(rep(1 / n, n) %*% X ^ 2) ^ 0.5
    X <- X / rep(Xscale, rep(n, p))
    # Yscale <- drop(rep(1 / n, n) %*% Y ^ 2) ^ 0.5
    # Y <- Y / rep(Yscale, n)
    
    X1 <- X[, c(1:length(label.names[[1]]))]
    X2 <- X[, c((length(label.names[[1]]) + 1):
                  (length(label.names[[1]]) + length(label.names[[2]])))]
    data0 <- as.data.frame(cbind(Y, X))
    
    result <- gpr(formula, label.names, data0, method, l, lambda, null.hypo = T)
    sigma2.n <- result[[1]]
    beta0 <- result[[2]]
    alpha0 <- result[[3]]
    K.gpr <- result[[4]]
    noise.hat <- NoiseEstimate(Y, sigma2.n, beta0, alpha0, K.gpr)
    mean.Y <- K.gpr %*% alpha0 + beta0
    
    # conduct bootstrap
    bs.test <- sapply(1:B, function(k){
      Y.star <- mean.Y + rnorm(n, sd = sqrt(noise.hat)) 
      dat <- cbind(Y.star, X)
      colnames(dat)[1] <- "Y"
      dat <- as.data.frame(dat)
      MultiScore2(Y ~ X1 + X2 + X1 * X2, label.names, 
                  data = dat, method, l, lambda)
    })
    
    # assemble test statistic
    original.test <- 
      MultiScore2(Y ~ X1 + X2 + X1 * X2, label.names, 
                  data = data0, method, l, lambda)
				  
    bs.pvalue <- sum(as.numeric(original.test) <= bs.test) / B

    return(list(pvalue = bs.pvalue, bs.test = bs.test, original.test = original.test))
  }

LinearBoost <- 
  function(formula, label.names, data = NULL, 
           method = NULL, B = 100){
    # Compute score tests comparing a fitted model and a more general alternative model.
    #
    # Args:
    #   formula:  A symbolic description of the model to be fitted.
    #   label.names: A character string indicating all the interior variables
    #                included in each predictor.
    #   data: A dataframe to be fitted.
    #   method: A character string indicating inear or quardratic is to be computed.
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
    # Ym <- mean(Y)
    p <- ncol(X)
    X <- X - rep(Xm, rep(n, p))
    # Y <- Y - Ym
    Xscale <- drop(rep(1 / n, n) %*% X ^ 2) ^ 0.5
    X <- X / rep(Xscale, rep(n, p))
    # Yscale <- drop(rep(1 / n, n) %*% Y ^ 2) ^ 0.5
    # Y <- Y / rep(Yscale, n)
    
    X1 <- X[, c(1:length(label.names[[1]]))]
    X2 <- X[, c((length(label.names[[1]]) + 1):
                  (length(label.names[[1]]) + length(label.names[[2]])))]
    data0 <- as.data.frame(cbind(Y, X))
    
    if (method == "linear"){
      lm <- lm(Y ~ X)
      beta0.est <- coef(lm)[1]
      y.est <- predict(lm)
      eps.est <- summary(lm)$sigma
      y.boot <- y.est + rnorm(n, sd = eps.est)
      
      true.stat <- LinearTestStat(Y, X, beta0.est)
      
      boot.stat <- sapply(1:B, function(k){
        y.boot <- y.est + rnorm(n, sd = eps.est)
        LinearTestStat(y.boot, X, beta0.est)
      })
    }
    else if (method == "quadratic"){
      X <- cbind(X, X ^ 2)
      Xscale <- drop(rep(1 / n, n) %*% X ^ 2) ^ 0.5
      X <- X / rep(Xscale, rep(n, 2 * p))
      
      lm <- lm(Y ~ X)
      beta0.est <- coef(lm)[1]
      y.est <- predict(lm)
      eps.est <- summary(lm)$sigma
      y.boot <- y.est + rnorm(n, sd = eps.est)
      
      true.stat <- LinearTestStat(Y, X, beta0.est)
      
      boot.stat <- sapply(1:B, function(k){
        y.boot <- y.est + rnorm(n, sd = eps.est)
        LinearTestStat(y.boot, X, beta0.est)
      })
    }
    else
      stop ("method must be linear or quadratic!")
    
    pvalue <- sum(as.numeric(true.stat) <= boot.stat) / B
    return(list(pvalue = pvalue, boot.stat = boot.stat, true.stat = true.stat))
  }

LinearTestStat <- 
  function(Y, X, beta0.est){
    p <- ncol(X) / 2
    X1 <- X[, 1:(p)]
    X2 <- X[, (1 + p):(2 * p)]
    
    X12 <- NULL
    for (i in 1:p){
      X12 <- cbind(X12, X1[, i] * X2)
    }
    X <- cbind(1, X)
    V.inv <- ginv(X %*% t(X))
    K12 <- X12 %*% t(X12)
    y.center <- matrix(Y - beta0.est, ncol = 1)
    
    t(y.center) %*% V.inv %*% K12 %*% V.inv %*% y.center
  }

