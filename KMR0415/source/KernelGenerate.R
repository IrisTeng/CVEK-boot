KernelGenerate <-
  function(method = NULL, Sigma = NULL, l = 1, p = 2){
    # Generate kernel function.
    #
    # Args:
    #   method: A character string indicating which Kernel is to be computed.
    #   l: A numeric number indicating the hyperparameter of a specific kernel.
    #   p: For polynomial, p is the power;
    #      for matern, v = p + 1 / 2;
    #      for rational, alpha = p.
    #   Sigma: Covariance matrix for neural network
    #
    # Returns:
    #   Kern: A matrix indicating the expected kernel function.
    
    if (method == "intercept")
      SE <- function(xp, xq, Sigma, l, p) 1
    else if (method == "linear")
      SE <- function(xp, xq, Sigma, l, p) t(xp) %*% xq
    else if (method == "polynomial")
      SE <- function(xp, xq, Sigma, l, p) (t(xp) %*% xq + 1) ^ p
    else if (method == "rbf")
      SE <- function(xp, xq, Sigma, l, p) exp(- sum((xp - xq) ^ 2) / (2 * l ^ 2))
    else if (method == "matern")
      SE <- function(xp, xq, Sigma, l, p){
        r <- sqrt(sum((xp - xq) ^ 2))
        v <- p + 1 / 2
        s <- 0
        for (i in 0:p) {
          s <- s + factorial(p + i) / (factorial(i) * factorial(p - i)) * 
            (sqrt(8 * v) * r / l) ^ (p - i)
        }
        k.v <- exp(-sqrt(2 * v) * r / l) * gamma(p + 1) / gamma(2 * p - 1) * s
        return(k.v)
      }
    else if (method == "rational")
      SE <- function(xp, xq, Sigma, l, p){
        r <- sqrt(sum((xp - xq) ^ 2))
        k.a <- (1 + r ^ 2 / (2 * p * l ^ 2)) ^ (- p)
        return(k.a)
      }
    else if (method == "nn")
      SE <- function(xp, xq, Sigma, l, p){
        xp <- c(1, xp)
        xq <- c(1, xq)
        if(is.null(Sigma)) Sigma <- diag(length(xp))
        s <- 2 * t(xp) %*% Sigma %*% xq / (sqrt((1 + 2 * t(xp) %*% Sigma %*% xp) 
                                                * (1 + 2 * t(xq) %*% Sigma %*% xq)))
        k.n <- asin(s)
        return(k.n)
      }
    else
      stop("method must be one kind of kernel!")
    
    Kern <- function(X2, X1) apply(X1, 1, function(xp){
      apply(X2, 1, function(xq){
        SE(xp, xq, Sigma, l, p)
      })
    })
    #Kern <- Kern / tr(Kern)
    return(Kern)
  }
