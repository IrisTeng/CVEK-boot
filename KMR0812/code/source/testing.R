#' Conducting Score Tests for Interaction
#'
#' Conduct score tests comparing a fitted model and a more general alternative
#' model.
#'
#' There are two tests available here:
#'
#' \bold{Asymptotic Test}
#'
#' This is based on the classical variance component test to construct a
#' testing procedure for the hypothesis about Gaussian process function.
#'
#' \bold{Bootstrap Test}
#'
#' When it comes to small sample size, we can use bootstrap test instead, which
#' can give valid tests with moderate sample sizes and requires similar
#' computational effort to a permutation test.
#'
#' @param formula_int A symbolic description of the model with interaction.
#' @param label_names A character string indicating all the interior variables
#' included in each predictor.
#' @param Y Reponses of the dataframe.
#' @param X1 The first type of factor in the dataframe (could contains several
#' subfactors).
#' @param X2 The second type of factor in the dataframe (could contains several
#' subfactors).
#' @param Kernlist The kernel library containing several kernels given by user.
#' @param mode A character string indicating which tuning parameter criteria is
#' to be used.
#' @param strategy A character string indicating which ensemble strategy is to
#' be used.
#' @param beta A numeric value specifying the parameter when strategy = "exp".
#' @param test A character string indicating which test is to be used.
#' @param lambda A numeric string specifying the range of noise to be chosen.
#' The lower limit of lambda must be above 0.
#' @param B A numeric value indicating times of resampling when test = "boot".
#' @return \item{pvalue}{p-value of the test.}
#' @author Wenying Deng
#' @seealso method: \code{\link{kernelGenerate}}
#'
#' mode: \code{\link{tuning}}
#'
#' strategy: \code{\link{ensemble}}
#' @references Xihong Lin. Variance component testing in generalised linear
#' models with random effects. June 1997.
#'
#' Arnab Maity and Xihong Lin. Powerful tests for detecting a gene effect in
#' the presence of possible gene-gene interactions using garrote kernel
#' machines. December 2011.
#'
#' Petra Bu ̊zˇkova ́, Thomas Lumley, and Kenneth Rice. Permutation and
#' parametric bootstrap tests for gene-gene and gene-environment interactions.
#' January 2011.
#' @examples
#'
#'
#' ##testing(formula_int = Y ~ X1 * X2,
#' ##label_names = list(X1 = c("x1", "x2"), X2 = c("x3", "x4")),
#' ##Y, X1, X2, Kernlist, mode = "loocv", strategy = "erm",
#' ##beta = 1, test = "boot", lambda = exp(seq(-5, 5)), B = 100)
#'
#'
#' @export testing
testing <- function(formula_int, label_names, Y, X1, X2, Kernlist,
                    mode = "loocv", strategy = "erm", beta = 1,
                    test = "boot", lambda = exp(seq(-5, 5)), B = 100){

  re <- genericFormula(formula_int, label_names)
  generic_formula0 <- re$generic_formula
  len <- re$length_main
  data <- as.data.frame(cbind(Y, X1, X2))
  colnames(data) <- c("Y", label_names[[1]], label_names[[2]])
  X <- model.matrix(generic_formula0, data)[, -1]
  X12 <- X[, c((len + 1):dim(X)[2])]
  n <- length(Y)

  result <- estimation(Y, X1, X2, Kernlist, mode, strategy, beta, lambda)
  lam <- result[[1]]
  beta0 <- result[[2]]
  alpha0 <- result[[3]]
  K_gpr <- result[[4]]
  u_weight <- result[[5]]
  noise_hat <- noiseEstimate(Y, lam, beta0, alpha0, K_gpr)
  sigma2_hat <- noise_hat[[1]]
  A <- noise_hat[[2]]
  tau_hat <- sigma2_hat / lam

  if (test == "boot"){
    # conduct bootstrap
    meanY <- K_gpr %*% alpha0 + beta0
    
    bs_test <- sapply(1:B, function(k) {
      
      Ystar <- meanY + rnorm(n, sd = sqrt(sigma2_hat))
      scoreStat(n, Ystar, X12, beta0, sigma2_hat, tau_hat, K_gpr)$test_stat
    })

    # assemble test statistic
    original_ans <-
      scoreStat(n, Y, X12, beta0, sigma2_hat, tau_hat, K_gpr)
    
    V0_inv_mat <- original_ans$V0_inv
    K12_mat <- original_ans$K12
    dif_vec <- Y - beta0
    V0_inv_dif <- V0_inv_mat %*% dif_vec
    
    dif_sum <- sum(dif_vec ^ 2)
    V0_inv_dif_sum <- sum(V0_inv_dif ^ 2)
    K12_tr <- tr(K12_mat)
    A_tr <- tr(A)
    
    original_test <- original_ans$test_stat

    pvalue <- sum(as.numeric(original_test) <= bs_test) / B
    
    out <- list(pvalue = pvalue, u_weight = u_weight,
                score_chi = original_test, bs_test = bs_test,
                tau_hat = tau_hat, lam = lam, 
                sigma2_hat = sigma2_hat, dif_sum = dif_sum, 
                V0_inv_dif_sum = V0_inv_dif_sum, 
                K12_tr = K12_tr, A_tr = A_tr, 
                dif_vec = dif_vec, V0_inv_dif = V0_inv_dif, 
                K12_mat = K12_mat, A = A)
  }
  else if (test == "asym"){
    asym_ans <-
      scoreStat(n, Y, X12, beta0, sigma2_hat, tau_hat, K_gpr)
    
    # dif_vec <- asym_ans$dif
    # V0_inv_mat <- asym_ans$V0_inv
    # K12_mat <- asym_ans$K12
    
    score_chi <- asym_ans$test_stat

    K0 <- K_gpr
    K12 <- X12 %*% t(X12)
    V0_inv <- ginv(tau_hat * K0 + sigma2_hat * diag(n))
    one <- rep(1, n)
    P0_mat <- V0_inv - V0_inv %*%
      one %*% ginv(t(one) %*% V0_inv %*% one) %*% t(one) %*% V0_inv

    drV0_tau <- K0
    drV0_sigma2 <- diag(n)
    drV0_del <- tau_hat * K12

    I0 <- infoMat(P0_mat,
                  mat_del = drV0_del, mat_sigma2 = drV0_sigma2,
                  mat_tau = drV0_tau)

    #Effective Info for delta
    tot_dim <- ncol(I0)
    I_deldel <-
      I0[1, 1] -
      I0[1, 2:tot_dim] %*% ginv(I0[2:tot_dim, 2:tot_dim]) %*% I0[2:tot_dim, 1]

    md <- tau_hat * tr(K12 %*% P0_mat) / 2

    m_chi <- I_deldel / (2 * md)
    d_chi <- md / m_chi

    pvalue <- 1 - pchisq(score_chi / m_chi, d_chi)
    out <- list(pvalue = pvalue, u_weight = u_weight,
                score_chi = score_chi, bs_test = 0, 
                tau_hat = tau_hat, lam = lam, 
                sigma2_hat = sigma2_hat, dif_sum = 0, 
                V0_inv_dif_sum = 0, 
                K12_tr = 0, A_tr = 0, 
                dif_vec = 0, V0_inv_dif = 0, 
                K12_mat = 0, A = 0)
  }
  else
    stop("test must be boot or asym!")

  return(out)
}
