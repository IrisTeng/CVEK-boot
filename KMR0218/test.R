GprTest <- 
  function(formula, label.names, 
           data = NULL, method = NULL, l = 1){
    # An implementation of Gaussian processes for regression.
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
    #   sigma.n: A numeric number chosen as the best size of noise via LooCV.
    #   intercept: A numberic number indicating the bias weight.
    #   alpha: A numeric string indicating the projection vector.
    
    # prepare data
    y <- data[, as.character(attr(terms(formula), "variables"))[2]]
    
    # extract the information from the given formula to construct a true formula
    generic.formula0 <- 
      GenericFormula(formula, label.names)$generic.formula
    
    # extract design matrix
    X <- model.matrix(generic.formula0, data)[, -1]
    n <- nrow(X)
    
    # estimation
    Kern <- KernelGenerate(method, l)
    K <- Kern(X, X)
    
    # # visualizaton playground
    # x <- matrix(seq(-10, 10, 0.1), ncol = 1)
    # K <- Kern(x, x)
    # plot(x, y, type = "n", ylim = c(-5, 5))
    # for (i in 1:10){
    #   y <- K %*% rnorm(length(x))
    #   lines(x, y, col = i)
    # }
    # #
    
    
    # lambda0 <- LooCV(formula, label.names, data, method, l, lambda)
    
    K1 <- cbind(1, K)
    K2 <- cbind(0, rbind(0, K))
    
    lambda0 <- 0.6
    theta <- ginv(lambda0 * K2 + t(K1) %*% K1) %*% t(K1) %*% y
    beta0 <- theta[1]
    alpha <- theta[-1]
    # one <- rep(1, n)
    # bb <- ginv(t(one) %*% (diag(n) - K %*% ginv(K + lambda0 * diag(n))) %*% one) %*% 
    #   (t(one) %*% (diag(n) - K %*% ginv(K + lambda0 * diag(n))) %*% y)
    y.star <- K %*% alpha + beta0
    
    return(list(sigma2.n = lambda0, intercept = beta0, alpha = alpha, y.hat = y.star, y=y, K=K))
  }


#label.names <- list(X1=c("x1","x2"),X2=c("x3","x4"))
#rawData <- OriginalData2(2,5,1,3,20,label_names=label_names,beta_int=0.1,scale=10,eps=1)
K0 <- GprTest(Y~X1+X2,label.names,data=rawData4,"gaussian",l=2)
noise.hat <- NoiseEstimate(K0$y, K0$sigma2.n, K0$intercept, K0$alpha, K.hat=K0$K)

gp <- gausspr(rawData4[,-1],rawData4[,1], kernel="rbfdot", 
              kpar=list(sigma=1/8), fit=T,scaled=FALSE, var=noise.hat)
cbind(K0$y.hat,fitted(gp))

rawData <- OriginalData2(200, label.names, "linear", int.effect = 1)
rawData2 <- OriginalData2(20, label.names, "gaussian", int.effect = 0)
rawData3 <- OriginalData2(20, label.names, "linear", int.effect = 1)
rawData4 <- OriginalData2(20, label.names, "gaussian", int.effect = 1)


rawData <- OriginalData2(size = n, label.names, method = "gaussian", 
                         int.effect = 0, eps = 0.01)
res <- gpr(formula, label.names, rawData, "gaussian", 
           l = 1, lambda = exp(seq(-5, 5, 1)), null.hypo = TRUE)
yhat <- res$intercept + res$K %*% res$alpha
plot(rawData$Y, yhat)
abline(a=0, b=1)

################################
# Experiment

Kern <- KernelGenerate(method, l = 1)
w1 <- rnorm(size)
# w2 <- rnorm(size)
w2 <- w1
w12 <- rnorm(size)
h0 <- K1 %*% w1 + K2 %*% w2
h0 <- h0 / sqrt(sum(h0 ^ 2))

trial_len <- 1e3
trial_container <- vector("numeric", trial_len)
for (trial in 1:trial_len){
  X1 <- rmvnorm(n = size,
                mean = rep(0, length(label.names[[1]])),
                sigma = diag(length(label.names[[1]])))
  X2 <- rmvnorm(n = size,
                mean = rep(0, length(label.names[[2]])),
                sigma = diag(length(label.names[[2]])))
  K1 <- Kern(X1, X1)
  K2 <- Kern(X2, X2)
  h1.prime <- (K1 * K2) %*% w12
  trial_container[trial] <- 
    eigen(K1 + K2, only.values=1)$values[1]
}

pb <- txtProgressBar(0, trial_len, style = 3)
for (trial in 1:trial_len){
  setTxtProgressBar(pb, trial)
  OriginalData2(size = 50, label.names, method = "gaussian", 
                int.effect = 0, eps = 0.01)
}


# 20180220
SanityCheck <- function(size = 200, method = NULL, 
                        D = 1000, eps = .01, B = 100){
  formula <- Y ~ X1 + X2
  label.names <- list(X1 = c("x1", "x2"), X2 = c("x3", "x4"))
  l <- 1
  lambda <- exp(seq(-5, 5, 1))
  X1 <- rmvnorm(n = size,
                mean = rep(0, length(label.names[[1]])),
                sigma = diag(length(label.names[[1]])))
  X2 <- rmvnorm(n = size,
                mean = rep(0, length(label.names[[2]])),
                sigma = diag(length(label.names[[2]])))
  
  Kern <- KernelGenerate(method, l = 1)
  w1 <- rnorm(size)
  w2 <- w1
  K1 <- Kern(X1, X1)
  K2 <- Kern(X2, X2)
  h0 <- K1 %*% w1 + K2 %*% w2
  h0 <- h0 / sqrt(sum(h0 ^ 2))
  bias <- rnorm(1)
  ds <- list()
  for (d in 1:D){
    set.seed(d)
    Y <- h0 + bias + rnorm(size, 0, eps)
    data <- as.data.frame(cbind(Y, X1, X2))
    colnames(data) <- c("Y", label.names[[1]], label.names[[2]])
    ds[[d]] <- data
  }
  
  original.test <- sapply(1:D, function(d){
    MultiScore2(Y ~ X1 + X2 + X1 * X2, label.names, 
                data = ds[[d]], method, l, lambda)
  }
  )
  
  test.stat <- matrix(NA, nrow = B, ncol = D)
  for(d in 1:D){
    test.stat[, d] <- TestBoot(formula, label.names, 
                               data = ds[[d]], method = method, B = B)
  }
  
  return(list(original = original.test, bootmat <- test.stat))
}



TestBoot <- 
  function(formula, label.names, data = NULL, 
           method = NULL, l = 1, lambda = exp(seq(-5, 5, 1)), B = 100){

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
      # K.star <- (as.matrix(X.star1) %*% t(X1)) + (as.matrix(X.star2) %*% t(X2))
      mean.Y <- K.star %*% alpha0 + beta0
      Y.star <- mean.Y + rnorm(n, mean = 0, sd = sqrt(noise.hat)) 
      dat <- cbind(Y = Y.star, X.star)
      MultiScore2(Y ~ X1 + X2 + X1 * X2, label.names, 
                  data = dat, method, l, lambda)
    })

    return(bs.test)
  }