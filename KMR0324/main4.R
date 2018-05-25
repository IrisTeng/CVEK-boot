setwd('/Users/iristeng/Desktop/KMR/KMR0324/source')
library(mvtnorm)
library(MASS)
library(snowfall)
library(psych)
library(limSolve)

formula <- Y ~ X1 + X2
label.names <- list(X1 = c("x1","x2"), X2 = c("x3","x4"))

Kernlist <- NULL
# intercept kernel
Kernlist <- c(Kernlist, KernelGenerate('intercept'))
# linear kernel
Kernlist <- c(Kernlist, KernelGenerate('linear'))
# polynomial kernel, p = 2
Kernlist <- c(Kernlist, KernelGenerate('polynomial', p = 2))
# polynomial kernel, p = 3
Kernlist <- c(Kernlist, KernelGenerate('polynomial', p = 3))
# polynomial kernel, p = 4
# Kernlist <- c(Kernlist, KernelGenerate('polynomial', p = 4))
# rbf kernel, l = .3
# Kernlist <- c(Kernlist, KernelGenerate('rbf', l = .3))
# rbf kernel, l = .6
Kernlist <- c(Kernlist, KernelGenerate('rbf', l = .6))
# rbf kernel, l = 1
Kernlist <- c(Kernlist, KernelGenerate('rbf', l = 1))
# rbf kernel, l = 2
Kernlist <- c(Kernlist, KernelGenerate('rbf', l = 2))
# rbf kernel, l = 3
Kernlist <- c(Kernlist, KernelGenerate('rbf', l = 3))
# rbf kernel, l = 4
# Kernlist <- c(Kernlist, KernelGenerate('rbf', l = 4))
# rbf kernel, l = 5
# Kernlist <- c(Kernlist, KernelGenerate('rbf', l = 5))




data <- OriginalData2(200, label.names, "linear")
res <- Ensemble(formula, label.names, Kernlist, data)
dd <- KernelBoot(formula, label.names, Kernlist, data)

Ensemble2 <-
  function(formula, label.names, Kernlist,
           data = NULL, lambda = exp(seq(-5, 5, 1))){

    y <- data[, as.character(attr(terms(formula), "variables"))[2]]

    # extract the information from the given formula to construct a true formula
    re <- GenericFormula(formula, label.names)
    generic.formula0 <- re[[1]]
    len <- re[[2]]

    # extract design matrix
    X <- model.matrix(generic.formula0, data)[, -1]
    n <- nrow(X)

    # estimation
    # D = 10
    X1 <- X[, c(1:length(label.names[[1]]))]
    X2 <- X[, c((length(label.names[[1]]) + 1):len)]

    D <- length(Kernlist)
    # lambda.mat <- matrix(0, nrow = D, ncol = n + 1)
    error.mat <- matrix(0, nrow = n, ncol = D)

    for (d in seq(D)){
      Kern <- Kernlist[[d]]
      K <- Kern(X1, X1) + Kern(X2, X2)
      if (length(lambda) != 1){
        lambda0 <- LooCV(y, K, lambda)
        K1 <- cbind(1, K)
        K2 <- cbind(0, rbind(0, K))
        theta <- ginv(lambda0 * K2 + t(K1) %*% K1) %*% t(K1) %*% y
        # error.mat[, d] <- y - K1 %*% theta
        beta0 <- theta[1]
        M <- K %*% ginv(K + lambda0 * diag(n))
        error.mat[, d] <- (diag(n) - M) %*% (y - beta0) / diag(diag(n) - M)
      }
    }

    A <- error.mat
    B <- rep(0, n)
    E <- rep(1, D)
    F <- 1
    G <- diag(D)
    H <- rep(0, D)
    u.hat <- lsei(A, B, E = E, F = F, G = G, H = H)$X

    return(u.hat)
  }

## without interaction effect
w <- NULL
pb <- txtProgressBar(0, 20, style = 3)
for(method in c('linear', 'polynomial', 'rbf')){
  for(i in seq(20)){
    setTxtProgressBar(pb, i)
    data0 <- OriginalData2(100, label.names, method = method)
    w <- rbind(w, Ensemble2(formula, label.names, Kernlist, data = data0))
  }
}

w <- as.data.frame(round(w, 4))
w$generate.method <- rep(c('linear', 'polynomial', 'rbf'), rep(20, 3))
colnames(w) <- c('intercept', 'linear', 'polynomial p=2', 'polynomial p=3',
                 paste0('rbf l=', c(0.6, 1, 2, 3)), 'generate.method')
write.csv(w, file = 'weights.csv')


## with interaction effect=5
w <- NULL
pb <- txtProgressBar(0, 20, style = 3)
for(method in c('linear', 'polynomial', 'rbf')){
  for(i in seq(20)){
    setTxtProgressBar(pb, i)
    data0 <- OriginalData2(100, label.names, method = method, int.effect = 5)
    w <- rbind(w, Ensemble2(formula, label.names, Kernlist, data = data0))
  }
}

w <- as.data.frame(round(w, 4))
w$generate.method <- rep(c('linear', 'polynomial', 'rbf'), rep(20, 3))
colnames(w) <- c('intercept', 'linear', 'polynomial p=2', 'polynomial p=3',
                 paste0('rbf l=', c(0.6, 1, 2, 3)), 'generate.method')
write.csv(w, file = 'weights_int.csv')
