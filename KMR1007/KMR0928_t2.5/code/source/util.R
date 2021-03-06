#' From Vectors to Single Variables
#'
#' Transform format of predictors from vectors to single variables.
#'
#'
#' @param formula (formula) A symbolic description of the model to be fitted.
#' @param label_names (list) A character string indicating all the interior 
#' variables included in each predictor.
#' @return \item{generic_formula}{(formula) A symbolic description of the 
#' model written in single variables format.}
#'
#' \item{length_main}{(integer) A numeric value indicating the length of 
#' main effects.}
#' @author Wenying Deng
#' @examples
#'
#'
#' generic_formula0 <- generate_formula(formula = Y ~ X1 + X2,
#' label_names = list(X1 = c("x1", "x2"), X2 = c("x3", "x4")))
#'
#'
#' @export generate_formula
generate_formula <-
  function(formula, label_names) {
    
    formula_factors <- attr(terms(formula), "factors")
    generic_formula <- Y ~ 1
    length_main <- 0
    for (i in 1:dim(formula_factors)[2]) {
      terms_names <-
        rownames(formula_factors)[which(formula_factors[, i] == 1)]
      if (length(terms_names) == 1) {
        generic_formula <-
          update.formula(generic_formula,
                         as.formula(paste(as.character(
                           attr(terms(formula), "variables"))[2],
                           paste(label_names[[terms_names]], 
                                 collapse=" + "), sep=" ~ .+")))
        length_main <- length_main + length(label_names[[terms_names]])
      } else {
        interaction_formula <-
          paste("(", paste(label_names[[terms_names[1]]], 
                           collapse=" + "), ")", sep="")
        for (j in 2:length(terms_names)) {
          interaction_formula <-
            paste(interaction_formula, "*(",
                  paste(label_names[[terms_names[j]]], 
                        collapse=" + "), ")", sep=" ")
        }
        generic_formula <-
          update.formula(generic_formula,
                         as.formula(paste(as.character(
                           attr(terms(formula), "variables"))[2], 
                           interaction_formula, sep=" ~ .+")
                         )
          )
      }
    }
    
    list(generic_formula =  generic_formula, length_main = length_main)
  }




#' Generating Original Data
#'
#' Generate original data based on specific kernels.
#'
#' This function generates with a specific kernel. The argument int_effect
#' represents the strength of interaction relative to the main effect since all
#' sampled functions have been standardized to have unit norm.
#'
#' @param n (integer) A numeric number specifying the number of observations.
#' @param label_names (list) A character string indicating all the interior 
#' variables included in each predictor.
#' @param method (character) A character string indicating which kernel is 
#' to be computed.
#' @param int_effect (numeric) A numeric number specifying the size of interaction.
#' @param l (numeric) A numeric number indicating the hyperparameter 
#' (flexibility) of a specific kernel.
#' @param p (integer) For polynomial, p is the power; for matern, v = p + 1 / 2; 
#' for rational, alpha = p.
#' @param eps (numeric) A numeric number indicating the size of noise.
#' @return \item{data}{(dataframe, n*P) A dataframe to be fitted.}
#' @author Wenying Deng
#' @examples
#'
#'
#' data <- generate_data(n = 100, label_names =
#' list(X1 = c("x1", "x2"), X2 = c("x3", "x4")),
#' method = "rbf", int_effect = 0, l = 1, p = 2, eps = .01)
#'
#'
#' @export generate_data
generate_data <-
  function(n, label_names, 
           method = "rbf", int_effect = 0,
           l = 1, p = 2, eps = .01) {
    
    X1 <- rmvnorm(n = n,
                  mean = rep(0, length(label_names[[1]])),
                  sigma = diag(length(label_names[[1]])))
    X2 <- rmvnorm(n = n,
                  mean = rep(0, length(label_names[[2]])),
                  sigma = diag(length(label_names[[2]])))
    
    kern <- generate_kernel(method = method, l = l, p = p)
    w1 <- rnorm(n)
    w2 <- w1
    w12 <- rnorm(n)
    K1 <- kern(X1, X1)
    K2 <- kern(X2, X2)
    K1 <- K1 / tr(K1)
    K2 <- K2 / tr(K2)
    h0 <- K1 %*% w1 + K2 %*% w2
    h0 <- h0 / sqrt(sum(h0 ^ 2))
    
    h1_prime <- (K1 * K2) %*% w12
    Ks <- svd(K1 + K2)
    if (length(Ks$d / sum(Ks$d) > .001) > 0) {
      len <- length(Ks$d[Ks$d / sum(Ks$d) > .001])
      U0 <- Ks$u[, 1:len]
      h1_prime_hat <- fitted(lm(h1_prime ~ U0))
      h1 <- h1_prime - h1_prime_hat
      if (all(h1 == 0)) {
        warning("interaction term colinear with main-effect space!")
        h1 <- h1_prime
        h1 <- h1 / sqrt(sum(h1 ^ 2))
      } else {
        h1 <- h1 / sqrt(sum(h1 ^ 2))
      }
    } else {
      warning("largest eigen value smaller than 0.001!")
      h1 <- h1_prime
      h1 <- h1 / sqrt(sum(h1 ^ 2))
    }
    
    Y <- h0 + int_effect * h1 + rnorm(1) + rt(n, 2.5)
    data <- as.data.frame(cbind(Y, X1, X2))
    colnames(data) <- c("Y", label_names[[1]], label_names[[2]])
    
    data
  }





#' Estimating Noise
#'
#' An implementation of Gaussian processes for estimating noise.
#'
#'
#' @param Y (vector of length n) Reponses of the dataframe.
#' @param lambda_hat (numeric) The selected tuning parameter based on the 
#' estimated ensemble kernel matrix.
#' @param beta_hat (numeric) Estimated bias of the model.
#' @param alpha_hat (vector of length n) Estimated coefficients of the estimated 
#' ensemble kernel matrix.
#' @param K_hat (matrix, n*n) Estimated ensemble kernel matrix.
#' @return \item{sigma2_hat}{(numeric) The estimated noise of the fixed effects.}
#' 
#' \item{SSE}{(numeric) The estimated noise of the fixed effects.}
#' 
#' \item{A}{(matrix) The estimated noise of the fixed effects.}
#' @author Wenying Deng
#' @references Jeremiah Zhe Liu and Brent Coull. Robust Hypothesis Test for
#' Nonlinear Effect with Gaus- sian Processes. October 2017.
#' @examples
#'
#'
#' sigma2_hat <- estimate_noise(Y, lam, beta0, alpha0, K_gpr)
#'
#'
#' @export estimate_noise
estimate_noise <- function(Y, lambda_hat, beta_hat, alpha_hat, K_hat) {

  n <- nrow(K_hat)
  V <- lambda_hat * diag(n) + K_hat
  one <- rep(1, n)
  Px <- one %*% ginv(t(one) %*% ginv(V) %*% one) %*% t(one) %*% ginv(V)
  Pk <- K_hat %*% ginv(V) %*% (diag(n) - Px)
  A <- Px + Pk
  sigma2_hat <- sum((Y - beta_hat - K_hat %*% alpha_hat) ^ 2) / (n - tr(A) - 1)
  SSE <- sum((Y - beta_hat - K_hat %*% alpha_hat) ^ 2)
  
  list(sigma2_hat = sigma2_hat, SSE = SSE, A = A)
 }



#' Computing Score Test Statistics.
#'
#' Compute score test statistics.
#'
#' The test statistic is distributed as a scaled Chi-squared distribution.
#'
#' @param n (integer) A numeric number specifying the number of observations.
#' @param Y (vector of length n) Reponses of the dataframe.
#' @param X12 (dataframe, n*(p1\*p2)) The interaction items of first and second 
#' types of factors in the dataframe.
#' @param beta0 (numeric) Estimated bias of the model.
#' @param K_gpr (matrix, n*n) Estimated ensemble kernel matrix.
#' @param sigma2_hat (numeric) The estimated noise of the fixed effects.
#' @param tau_hat (numeric) The estimated noise of the random effects.
#' @return \item{test_stat}{(numeric) The computed test statistic.}
#' 
#' \item{V0_inv}{(numeric) The computed test statistic.}
#' @author Wenying Deng
#' @references Arnab Maity and Xihong Lin. Powerful tests for detecting a gene
#' effect in the presence of possible gene-gene interactions using garrote
#' kernel machines. December 2011.
#' @examples
#'
#'
#' compute_stat(n = 100, Y, X12, beta0, K_gpr, sigma2_hat, tau_hat)
#'
#'
#' @export compute_stat
compute_stat <-
  function(n, Y, X12, beta0, K_gpr, sigma2_hat, tau_hat) {

    K0 <- K_gpr
    K12 <- X12 %*% t(X12)
    
    V0_inv <- ginv(tau_hat * K0 + sigma2_hat * diag(n))
    
    # lam <- sigma2_hat / tau_hat
    # V0_inv <- ginv(K0 + lam * diag(n))
    
    test_stat <- tau_hat * t(Y - beta0) %*% V0_inv %*%
      K12 %*% V0_inv %*% (Y - beta0) / 2
    
    # test_stat <- t(Y - beta0) %*% V0_inv %*% 
    #   K12 %*% V0_inv %*% (Y - beta0) / 2
    
    list(test_stat = test_stat, V0_inv = V0_inv)
  }




#' Computing Information Matrices
#'
#' Compute information matrices based on block matrices.
#'
#' This function gives the information value of the interaction strength.
#'
#' @param P0_mat (matrix, n*n) Scale projection matrix under REML.
#' @param mat_del (matrix, n*n) Derivative of the scale covariance matrix of Y 
#' with respect to delta.
#' @param mat_sigma2 (matrix, n*n) Derivative of the scale covariance matrix 
#' of Y with respect to sigma2.
#' @param mat_tau (matrix, n*n) Derivative of the scale covariance matrix of Y 
#' with respect to tau.
#' @return \item{I0}{(matrix, n*n) The computed information value.}
#' @author Wenying Deng
#' @references Arnab Maity and Xihong Lin. Powerful tests for detecting a gene
#' effect in the presence of possible gene-gene interactions using garrote
#' kernel machines. December 2011.
#' @examples
#'
#'
#' I0 <- compute_info(P0_mat, mat_del = drV0_del,
#' mat_sigma2 = drV0_sigma2, mat_tau = drV0_tau)
#'
#'
#' @export compute_info
compute_info <-
  function(P0_mat, mat_del = NULL, mat_sigma2 = NULL, mat_tau = NULL) {
    
    I0 <- matrix(NA, 3, 3)

    I0[1, 1] <- tr(P0_mat %*% mat_del %*% P0_mat %*% mat_del) / 2  
    I0[1, 2] <- tr(P0_mat %*% mat_del %*% P0_mat %*% mat_sigma2) / 2
    I0[2, 1] <- I0[1, 2]
    
    I0[1, 3] <- tr(P0_mat %*% mat_del %*% P0_mat %*% mat_tau) / 2
    I0[3, 1] <- I0[1, 3]
  
    I0[2, 2] <- tr(P0_mat %*% mat_sigma2 %*% P0_mat %*% mat_sigma2) / 2
    
    I0[2, 3] <-  tr(P0_mat %*% mat_sigma2 %*% P0_mat %*% mat_tau) / 2  
    I0[3, 2] <- I0[2, 3]
    
    I0[3, 3] <-  tr(P0_mat %*% mat_tau %*% P0_mat %*% mat_tau) / 2  

    I0
  }

