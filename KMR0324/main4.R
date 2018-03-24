library(mvtnorm)
library(MASS)
library(snowfall)
library(psych)
library(limSolve)

formula <- Y~X1+X2
label.names <- list(X1=c("x1","x2"),X2=c("x3","x4"))
rawData <- OriginalData2(100, label.names, "linear", int.effect = 0)

Kernlist <- NULL
# linear kernel
Kernlist <- c(Kernlist, KernelGenerate('linear'))
# polynomial kernel, p = 2
Kernlist <- c(Kernlist, KernelGenerate('polynomial', p = 2))
# polynomial kernel, p = 3
Kernlist <- c(Kernlist, KernelGenerate('polynomial', p = 3))
# polynomial kernel, p = 4
Kernlist <- c(Kernlist, KernelGenerate('polynomial', p = 4))
# rbf kernel, l = .3
Kernlist <- c(Kernlist, KernelGenerate('rbf', l = .3))
# rbf kernel, l = .6
Kernlist <- c(Kernlist, KernelGenerate('rbf', l = .6))
# rbf kernel, l = 1
Kernlist <- c(Kernlist, KernelGenerate('rbf', l = 1))
# rbf kernel, l = 2
Kernlist <- c(Kernlist, KernelGenerate('rbf', l = 2))
# rbf kernel, l = 3
Kernlist <- c(Kernlist, KernelGenerate('rbf', l = 3))
# rbf kernel, l = 4
Kernlist <- c(Kernlist, KernelGenerate('rbf', l = 4))
# rbf kernel, l = 4
Kernlist <- c(Kernlist, KernelGenerate('rbf', l = 5))
# rbf kernel, l = 4
Kernlist <- c(Kernlist, KernelGenerate('rbf', l = 6))

Ensemble <- 
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
    lambda.mat <- matrix(0, nrow = D, ncol = n + 1)
    
    for (i in 1:D){
      Kern <- Kernlist[[i]]
      K <- Kern(X1, X1) + Kern(X2, X2)
      if (length(lambda) != 1) 
        lambda.mat[i, ] <- LooCV(y, K, lambda)
    }
    
    A <- t(lambda.mat[, 2:(n + 1)])
    B <- rep(0, 100)
    E <- rep(1, D)
    F <- 1
    G <- diag(D)
    H <- rep(0, D)
    u.hat <- lsei(A, B, E = E, F = F, G = G, H = H)$X
    return(u.hat)
  }

data <- OriginalData2(100, label.names, "polynomial", int.effect = 3)
h <- Ensemble(formula, label.names, Kernlist, data)
