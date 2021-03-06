KernelGenerate <-
  function(method = NULL, l = 1){
    # Generate kernel function.
    #
    # Args:
    #   method: A character string indicating which Kernel is to be computed.
    #   l: A numeric number indicating the hyperparameter of a specific kernel.
    #
    # Returns:
    #   Kern: A matrix indicating the expected kernel function.
    
    if (method == "linear")
      SE <- function(xp, xq, l) t(xp) %*% xq
    else if (method == "quadratic")
      SE <- function(xp, xq, l) (t(xp) %*% xq + 1) ^ 2
    else if (method == "gaussian")
      SE <- function(xp, xq, l) exp(- sum((xp - xq) ^ 2) / (2 * l ^ 2))
    else
      stop("method must be linear, quadratic or gaussian")
    
    Kern <- function(X2, X1) apply(X1, 1, function(xp){
      apply(X2, 1, function(xq){
        SE(xp, xq, l = l)
      })
    })
    return(Kern)
  }
