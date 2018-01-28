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
  function(size, label.names, 
           method = NULL, int.effect = 0, eps = 0.1){
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
    
    Kern <- KernelGenerate(method, l = 1)
    w1 <- rnorm(size)
    # w2 <- rnorm(size)
    w2 <- w1
    w12 <- rnorm(size)
    Y <- Kern(X1, X1) %*% w1 + Kern(X2, X2) %*% w2 + 
      int.effect * (Kern(X1, X1) * Kern(X2, X2)) %*% w12 + 
      rnorm(1) + rnorm(size, 0, eps)
    data <- as.data.frame(cbind(Y, X1, X2))
    colnames(data) <- c("Y", label.names[[1]], label.names[[2]])
    return(data)
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
    
    result <- gpr(Y ~ X1 + X2, label.names, data, method, l, lambda, null.hypo = F)
    sigma2.n <- result[[1]]
    beta0 <- result[[2]]
    
    # compute score statistic
    Kern <- KernelGenerate(method, l)
    K1 <- Kern(X1, X1)
    K2 <- Kern(X2, X2)
    K0 <- K1 + K2
    K12 <- X12 %*% t(X12)
    V0 <- K0 + sigma2.n * diag(n)
    test.stat <- t(y - beta0) %*% ginv(V0) %*% K12 %*% ginv(V0) %*% (y - beta0)
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
    Ym <- mean(Y)
    p <- ncol(X)
    X <- X - rep(Xm, rep(n, p))
    Y <- Y - Ym
    Xscale <- drop(rep(1/n, n) %*% X^2)^0.5
    X <- X/rep(Xscale, rep(n, p))
    
    X1 <- X[, c(1:length(label.names[[1]]))]
    X2 <- X[, c((length(label.names[[1]]) + 1):
                  (length(label.names[[1]]) + length(label.names[[2]])))]
    data0 <- as.data.frame(cbind(Y, X))
    
    result <- gpr(formula, label.names, data0, method, l, lambda, null.hypo = T)
    sigma2.n <- result[[1]]
    beta0 <- result[[2]]
    alpha <- result[[3]]
    Kern <- KernelGenerate(method, l)
    
    inds <- 1:n
    
    # conduct bootstrap
    bs.test <- sapply(1:B, function(k){
      bs.ind <- sample(inds, size = n, replace = T)
      X.star <- data0[bs.ind, -1]
      X.star1 <- X.star[, c(1:length(label.names[[1]]))]
      X.star2 <- X.star[, c((length(label.names[[1]]) + 1):
                              (length(label.names[[1]]) + length(label.names[[2]])))]
      K.star <- Kern(X.star1, X1) + Kern(X.star2, X2)
      mean.Y <- K.star %*% alpha + beta0
      Y <- rnorm(n, mean = mean.Y, sd = sqrt(sigma2.n))
      dat <- cbind(Y, X.star)
      MultiScore2(Y ~ X1 + X2 + X1 * X2, label.names, 
                  data = dat, method, l, lambda)
    })
    
    # assemble test statistic
    original.test <- MultiScore2(Y ~ X1 + X2 + X1 * X2, label.names, 
                                 data = data0, method, l, lambda)
    
    bs.pvalue <- sum(as.numeric(original.test) <= bs.test) / B
    return(bs.pvalue)
  }