genericFormula <- 
  function(formula, label_names){
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
    
    formula_factors <- attr(terms(formula), "factors")
    generic_formula <- Y ~ 1
    length_main <- 0
    for (i in 1:dim(formula_factors)[2]){
      terms_names <- 
        rownames(formula_factors)[which(formula_factors[, i] == 1)]
      if (length(terms_names) == 1){
        generic_formula <- 
          update.formula(generic_formula,
                         as.formula(paste(as.character(
                           attr(terms(formula), "variables"))[2],
                           paste(label_names[[terms_names]], collapse=" + "), sep=" ~ .+")))
        length_main <- length_main + length(label_names[[terms_names]])
      } else {
        interaction_formula <- 
          paste("(", paste(label_names[[terms_names[1]]], collapse=" + "), ")", sep="")
        for (j in 2:length(terms_names)){
          interaction_formula <- 
            paste(interaction_formula, "*(",
                  paste(label_names[[terms_names[j]]], collapse=" + "), ")", sep=" ")
        }
        generic_formula <- 
          update.formula(generic_formula,
                         as.formula(paste(as.character(
                           attr(terms(formula), "variables"))[2], interaction_formula, sep=" ~ .+")
                         )
          )
      }
    }
    return(list(generic_formula =  generic_formula, length_main = length_main))
  }


dataGenerate <- 
  function(size, label_names, method = "rbf", int_effect = 0,
           l = 1, p = 2, eps = 0.01){
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
    
    X1 <- rmvnorm(n = size,
                  mean = rep(0, length(label_names[[1]])),
                  sigma = diag(length(label_names[[1]])))
    X2 <- rmvnorm(n = size,
                  mean = rep(0, length(label_names[[2]])),
                  sigma = diag(length(label_names[[2]])))
    
    Kern <- kernelGenerate(method = method, l = l, p = p)
    w1 <- rnorm(size)
    w2 <- w1
    w12 <- rnorm(size)
    K1 <- Kern(X1, X1)
    K2 <- Kern(X2, X2)
    K1 <- K1 / tr(K1)
    K2 <- K2 / tr(K2)
    h0 <- K1 %*% w1 + K2 %*% w2
    h0 <- h0 / sqrt(sum(h0 ^ 2))
    
    h1_prime <- (K1 * K2) %*% w12
    Ks <- svd(K1 + K2)
    if(length(Ks$d / sum(Ks$d) > 0.001) > 0){
      len <- length(Ks$d[Ks$d / sum(Ks$d) > 0.001])
      U0 <- Ks$u[, 1:len]
      h1_prime_hat <- fitted(lm(h1_prime ~ U0))
      h1 <- h1_prime - h1_prime_hat
      if(all(h1 == 0)){
        warning("interaction term colinear with main-effect space!")
        h1 <- h1_prime
        h1 <- h1 / sqrt(sum(h1 ^ 2))
      }
      else
        h1 <- h1 / sqrt(sum(h1 ^ 2))
    }
    else{
      warning("largest eigen value smaller than 0.001!")
      h1 <- h1_prime
      h1 <- h1 / sqrt(sum(h1 ^ 2))
    }
    
    Y <- h0 + int_effect * h1 + rnorm(1) + rnorm(size, 0, eps)
    data <- as.data.frame(cbind(Y, X1, X2))
    colnames(data) <- c("Y", label_names[[1]], label_names[[2]])
    return(data)
  }


noiseEstimate <- function(Y, lambda_hat, beta_hat, alpha_hat, K_hat){
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
  
  n <- nrow(K_hat)
  V <- lambda_hat * diag(n) + K_hat
  one <- rep(1, n)
  Px <- one %*% ginv(t(one) %*% ginv(V) %*% one) %*% t(one) %*% ginv(V)
  Pk <- K_hat %*% ginv(V) %*% (diag(n) - Px)
  A <- Px + Pk
  sigma2_hat <- sum((Y - beta_hat - K_hat %*% alpha_hat) ^ 2) / (n - tr(A))
  return(sigma2_hat)
 }


multiScore <- 
  function(formula_int, label_names, Kernlist, data, 
           mode, strategy, beta, lambda){
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
    y <- data[, as.character(attr(terms(formula_int), "variables"))[2]]
    re <- genericFormula(formula_int, label_names)
    generic_formula0 <- re[[1]]
    len <- re[[2]]
    X <- model.matrix(generic_formula0, data)[, -1]
    n <- nrow(X)
    X1 <- X[, c(1:length(label_names[[1]]))]
    X2 <- X[, c((length(label_names[[1]]) + 1):len)]
    X12 <- X[, c((len + 1):dim(X)[2])]
    
    result <- ensemble(Y ~ X1 + X2, label_names, Kernlist, data, 
                       mode, strategy, beta, lambda)
    
    lam <- result[[1]]
    beta0 <- result[[2]]
    alpha0 <- result[[3]]
    K_gpr <- result[[4]]
    sigma2_hat <- noiseEstimate(y, lam, beta0, alpha0, K_gpr)
    tau_hat <- sigma2_hat / lam

    K0 <- K_gpr
    K12 <- X12 %*% t(X12)
    
    V0_inv <- ginv(tau_hat * K0 + sigma2_hat * diag(n))
    
    test_stat <- tau_hat * t(y - beta0) %*% V0_inv %*%
      K12 %*% V0_inv %*% (y - beta0) / 2
    
    return(test_stat)
  }


infoMat <- 
  function(P0_mat, mat_del = NULL, mat_sigma2 = NULL, mat_tau = NULL){
    I0 <- matrix(NA, 3, 3)
    
    I0[1, 1] <- tr(P0_mat %*% mat_del %*% P0_mat %*% mat_del) / 2  ##del.del
    #
    I0[1, 2] <- tr(P0_mat %*% mat_del %*% P0_mat %*% mat_sigma2) / 2  ##del.sig
    I0[2, 1] <- I0[1, 2]
    #
    I0[1, 3] <- tr(P0_mat %*% mat_del %*% P0_mat %*% mat_tau) / 2  ##del.lam
    I0[3, 1] <- I0[1, 3]
    #
    #
    I0[2, 2] <- tr(P0_mat %*% mat_sigma2 %*% P0_mat %*% mat_sigma2) / 2  ##sig.sig
    #
    I0[2, 3] <-  tr(P0_mat %*% mat_sigma2 %*% P0_mat %*% mat_tau) / 2  ##sig.lam
    I0[3, 2] <- I0[2, 3]
    #
    I0[3, 3] <-  tr(P0_mat %*% mat_tau %*% P0_mat %*% mat_tau)/ 2  ##lam.lam
    
    return(I0)
  }

