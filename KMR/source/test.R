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
    
    return(list(sigma2.n = lambda0, intercept = beta0, alpha = alpha, y.hat = y.star))
  }


#label_names <- list(X1=c("x1","x2"),X2=c("x3","x4"))
#rawData <- originalData(2,5,1,3,20,label_names=label_names,beta_int=0.1,scale=10,eps=1)
K0 <- GprTest(Y~X1+X2,label.names,data=rawData4,"gaussian",l=2)

gp <- gausspr(rawData4[,-1],rawData4[,1], kernel="rbfdot", 
              kpar=list(sigma=1/8), fit=T,scaled=FALSE, var=0.6)
cbind(K0$y.hat,fitted(gp))

rawData <- originalData2(20, label.names, "linear", int_effect = 0)
rawData2 <- originalData2(20, label.names, "gaussian", int_effect = 0)
rawData3 <- originalData2(20, label.names, "linear", int_effect = 20)
rawData4 <- originalData2(20, label.names, "gaussian", int_effect = 20)
