# Hello World!
# from vectors to single variables
generic_formula <- function(formula, label_names = label_names){
  formula_factors <- attr(terms(formula),"factors")
  generic_formula <- Y~1
  length_main_effect <- 0
  for(i in 1:dim(formula_factors)[2]){
    terms_names <- rownames(formula_factors)[which(formula_factors[,i]==1)]
    if(length(terms_names)==1){
      generic_formula <- update.formula(generic_formula,as.formula(paste(as.character(attr(terms(formula),"variables"))[2],
                                                                         paste(label_names[[terms_names]], collapse=" + "), sep=" ~ .+")))
      length_main_effect <- length_main_effect+length(label_names[[terms_names]])
    }
    else{
      interaction_formula <- paste("(",paste(label_names[[terms_names[1]]], collapse=" + "),")",sep="")
      for(j in 2:length(terms_names)){
        interaction_formula <- paste(interaction_formula,"*(",paste(label_names[[terms_names[j]]], collapse=" + "),")",sep=" ")
      }
      generic_formula <- update.formula(generic_formula,as.formula(paste(as.character(attr(terms(formula),"variables"))[2],interaction_formula, sep=" ~ .+")))
    }
  }
  return(list(generic_formula=generic_formula,length_main_effect=length_main_effect))
}

#generic_formula(Y~X1+X2+X3+X2*X3,list(X1=c("x1", "x2","x3"),X2=c("x4","x5"),X3=c("x6","x7")))

# Step 1: Implement Linear Regression from scratch. (using matrix multiplication)
linearReg <- function(formula, data = NULL,label_names = label_names){
  y <- data[,as.character(attr(terms(formula),"variables"))[2]]
  # extract the information from the given formula to construct a true formula
  generic_formula0=generic_formula(formula=formula,label_names = label_names)$generic_formula
  # extract design matrix
  x <- model.matrix(generic_formula0 ,data=data)
  if(abs(det(t(x)%*%x))<1e-6)
    stop("X must be a non-singluar matrix")
  beta <- solve(t(x)%*%x)%*%t(x)%*%y
  residual <- y-x%*%beta
  se <- sqrt(sum(residual^2)/(length(y)-dim(x)[2]))
  se.beta <- se*sqrt(diag(solve(t(x)%*%x)))
  t.value <- as.vector(beta/se.beta)
  p.value <- 2*pnorm(abs(t.value), lower.tail = F)
  ret <- list(beta=as.vector(beta), se.beta=se.beta, t.value=t.value, p.value=p.value, se=se)
  return(ret)
}

# linearReg(Y~X1+X2,data=dat,label_names=label_names)$beta
# linearReg(Y~X1+X2+X1*X2,data=dat,label_names=label_names)$beta

# Step 2: Derive score test for multivariate interaction term. Implement in R. 
# Check with R implementation, see if results agree.
multiScore <- function(formula, data = NULL,label_names = label_names){
  y0 <- data[,as.character(attr(terms(formula),"variables"))[2]]
  generic_formula0=generic_formula(formula=formula,label_names = label_names)$generic_formula
  len <- generic_formula(formula=formula,label_names = label_names)$length_main_effect
  x <- model.matrix(generic_formula0 ,data=data)
  x_main <- x[,c(0:len+1)]
  if(abs(det(t(x_main)%*%x_main))<1e-6)
    stop("X must be a non-singluar matrix")
  beta0 <- solve(t(x_main)%*%x_main)%*%t(x_main)%*%y0
  y=y0-x_main%*%beta0
  se <- sqrt(sum(y^2)/length(y0))
  U0 <- t(x)%*%y/se^2
  I0 <- t(x)%*%x/se^2
  score <- t(U0)%*%solve(I0)%*%U0
  k <- dim(x)[2]-(len+1)
  p.value <- pchisq(score, k, lower.tail = F)
  ret <- list(test.stat=score, df=k, p.value=p.value)
  return(ret)
}

# set.seed(306)
# dat <- as.data.frame(matrix(runif(240, -10, 10), ncol = 8))
# colnames(dat) <- c("Y", "x1", "x2","x3","x4","x5","x6","x7")
# label_names <- list(X1=c("x1","x2","x3"),X2=c("x4","x5","x6","x7"))
# 
# multiScore(Y~X1+X2+X1*X2,data=dat,label_names=label_names)
# fit <- glm(Y ~ x1 + x2 + x3 + x4 + x5 + x6 + x7, family = gaussian, data = dat)
# attach(dat)
# x99 <- cbind(x1*x4, x1*x5 , x1*x6 ,x1*x7 , x2*x4 , x2*x5 , x2*x6 , x2*x7 , x3*x4 , x3*x5 , x3*x6 ,x3*x7)
# detach(dat)
# summary(mdscore(fit, x99))

# Step 3: Implement the bootstrap testing procedure using Score test statistic.
# genarate original data, from which we conduct bootstrap resampling
originalData <- function(X1_mean, X1_sd, X2_mean, X2_sd, size, label_names = label_names, beta_int=0){
  X1 <- matrix(rnorm(size*length(label_names[[1]]), X1_mean, X1_sd), ncol = length(label_names[[1]]))
  X2 <- matrix(rnorm(size*length(label_names[[2]]), X2_mean, X2_sd), ncol = length(label_names[[2]]))
  beta1 <- runif(length(label_names[[1]]), 2, 5)
  beta2 <- runif(length(label_names[[2]]), 1, 3)
  if(beta_int!=0){
    X_int <- NULL
    for (i in 1:length(label_names[[1]])){
      X_int <- cbind(X_int, X1[,i]*X2)
    }
    X <- cbind(rep(1,size), X1, X2, X_int)
    beta <- t(cbind(2,t(beta1),t(beta2),t(rep(beta_int, length(label_names[[1]])*length(label_names[[2]])))))
  }
  else{
    X <- cbind(rep(1,size), X1, X2)
    beta <- t(cbind(2,t(beta1),t(beta2)))
  }
  Y <- X%*%beta + rnorm(size, 0, 1)
  data <- as.data.frame(cbind(Y,X1,X2))
  colnames(data) <- c("Y", label_names[[1]],label_names[[2]])
  return(data)
}
# rawData <- originalData(2,2,0,1,100,label_names = label_names)
# rawData2 <- originalData(2,2,0,1,100,label_names = label_names,beta_int = .1)

linearBoot <- function(formula, data = NULL,label_names = label_names, B=100){
  coef <- linearReg(formula,data=data,label_names=label_names)$beta
  sd <- linearReg(formula,data=data,label_names=label_names)$se
  n <- dim(data)[1]
  inds <- 1:n
  bs_test <- NULL
  for (j in 1:B){
    bs_ind <- sample(inds, size = n, replace = T)
    dat <- data[bs_ind,-1]
    dat <- cbind(rep(1,n),dat)
    Y <- as.matrix(dat)%*%coef + rnorm(n, 0, sd)
    dat <- cbind(Y,dat[,-1])
    bs_test <- c(bs_test, multiScore(Y~X1+X2+X1*X2,data=dat,label_names=label_names)$test.stat)
  }
  original_test <- multiScore(Y~X1+X2+X1*X2,data=data,label_names=label_names)$test.stat
  bs_pvalue <- sum(as.numeric(original_test)<=bs_test)/B
  return(bs_pvalue)
}


repeatBoot <- function(formula, label_names = label_names,beta=0, n, B=100, M=200){
  pvalue <- NULL
  for (i in 1:M){
    rawData <- originalData(0,1,0,1,size=n,label_names, beta_int=beta)
    pvalue <- c(pvalue, linearBoot(formula, data = rawData, label_names = label_names, B=B))
  }
  prob <- sum(pvalue<0.05)/M
  return(prob)
}

label_names <- list(X1=c("x1","x2"),X2=c("x3","x4"))

# for (MM in seq(200, 1000, 100)){
#   for (nn in seq(100, 1000, 100)){
#     prob <- repeatBoot(Y~X1+X2, label_names = label_names, n=nn, B=500, M=MM)
#   }
#   print(MM,nn,prob)
# }
# 
# for (MM in seq(200, 1000, 100)){
#   for (nn in seq(100, 1000, 100)){
#     prob <- repeatBoot(Y~X1+X2, label_names = label_names, n=nn, B=1000, M=MM)
#   }
#   print(MM,nn,prob)
# }

for (b in seq(0, 10, 0.1)){
  for (nn in seq(100, 1000, 100)){
    prob <- repeatBoot(Y~X1+X2, label_names = label_names, beta=b,n=nn, B=500, M=200)
  }
  print(b,nn,prob)
}

for (b in seq(0, 10, 0.1)){
  for (nn in seq(100, 1000, 100)){
    prob <- repeatBoot(Y~X1+X2, label_names = label_names, beta=b,n=nn, B=1000, M=200)
  }
  print(b,nn,prob)
}