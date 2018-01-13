# Hello World!
# from vectors to single variables
generic_formula <- 
  function(formula, label_names = label_names){
    formula_factors <- attr(terms(formula),"factors")
    generic_formula <- Y~1
    length_main_effect <- 0
    for(i in 1:dim(formula_factors)[2]){
      terms_names <- 
        rownames(formula_factors)[which(formula_factors[,i]==1)]
      if(length(terms_names)==1){
        generic_formula <- 
          update.formula(generic_formula,
                         as.formula(paste(as.character(
                           attr(terms(formula),"variables"))[2],
                           paste(label_names[[terms_names]], collapse=" + "), sep=" ~ .+")))
        length_main_effect <- length_main_effect+length(label_names[[terms_names]])
      }
      else{
        interaction_formula <- 
          paste("(",paste(label_names[[terms_names[1]]], collapse=" + "),")",sep="")
        for(j in 2:length(terms_names)){
          interaction_formula <- 
            paste(interaction_formula,"*(",
                  paste(label_names[[terms_names[j]]], collapse=" + "),")",sep=" ")
        }
        generic_formula <- 
          update.formula(generic_formula,
                         as.formula(paste(as.character(
                           attr(terms(formula),"variables"))[2], interaction_formula, sep=" ~ .+")
                         )
          )
      }
    }
    return(list(generic_formula=generic_formula,
                length_main_effect=length_main_effect))
  }

#generic_formula(Y~X1+X2+X3+X2*X3,list(X1=c("x1", "x2","x3"),X2=c("x4","x5"),X3=c("x6","x7")))

# Step 1: Implement Linear Regression from scratch. (using matrix multiplication)
linearReg <- 
  function(formula, data = NULL,label_names = label_names){
    # prepare data
    y <- data[,as.character(attr(terms(formula),"variables"))[2]]
    
    # extract the information from the given formula to construct a true formula
    generic_formula0 <- 
      generic_formula(formula=formula,label_names = label_names)$generic_formula
    
    # extract design matrix
    x <- model.matrix(generic_formula0 ,data=data)
    if(abs(det(t(x)%*%x))<1e-6)
      stop("X must be a non-singluar matrix")
    
    # estimation
    beta <- solve(t(x)%*%x)%*%t(x)%*%y
    residual <- y-x%*%beta
    se <- sqrt(sum(residual^2)/(length(y)-dim(x)[2]))
    se.beta <- se*sqrt(diag(solve(t(x)%*%x)))
    t.value <- as.vector(beta/se.beta)
    p.value <- 2*pt(abs(t.value), length(y)-dim(x)[2], lower.tail = F)
    
    # result
    ret <- list(beta=as.vector(beta), 
                se.beta=se.beta, 
                t.value=t.value, 
                p.value=p.value, 
                se=se)
    return(ret)
  }

# linearReg(Y~X1+X2,data=dat,label_names=label_names)$beta
# linearReg(Y~X1+X2+X1*X2,data=dat,label_names=label_names)$beta

# Step 12: Implement Ridge Regression from scratch. (using matrix multiplication)
library(psych)
lmRidge <- 
  function (formula, data=NULL, lambda = 0,label_names = label_name) {
    # prepare data
    Y <- data[,as.character(attr(terms(formula),"variables"))[2]]
    
    # extract the information from the given formula to construct a true formula
    generic_formula0 <- 
      generic_formula(formula=formula,label_names = label_names)$generic_formula
    
    # extract design matrix
    X <- model.matrix(generic_formula0, data=data)
    if(abs(det(t(X)%*%X))<1e-6)
      stop("X must be a non-singluar matrix")
    
    # standardize predictors and response
    n <- nrow(X)
    Xm <- colMeans(X[, -1])
    Ym <- mean(Y)
    p <- ncol(X) - 1
    X <- X[, -1] - rep(Xm, rep(n, p))
    Y <- Y - Ym
    Xscale <- drop(rep(1/n, n) %*% X^2)^0.5
    X <- X/rep(Xscale, rep(n, p))
    
    # vector of GCV values
    # GCV <- sapply(lambda, function(k){
    #   A <- X%*%solve((t(X)%*%X+k*diag(p)))%*%t(X)
    #   (sum(((diag(n)-A)%*%Y)^2)/(tr(diag(n)-A))^2)
    # })
    
    Loocv <- sapply(lambda, function(k){
      A <- X%*%solve((t(X)%*%X+k*diag(p)))%*%t(X)
      sum(((diag(n)-A)%*%Y/diag(diag(n)-A))^2)
    })
    
    # calculate coefficients
    # lambda0 <- lambda[which(GCV==min(GCV))]
    lambda0 <- lambda[which(Loocv==min(Loocv))]
    beta <- solve(t(X)%*%X+lambda0*diag(p))%*%t(X)%*%Y/Xscale
    inte <- Ym-Xm%*%beta
    residual <- Y-X%*%beta
    se <- sqrt(sum(residual^2)/(n-p))
    ret <- list(beta=as.vector(c(inte,beta)), 
                lambda=lambda0, 
                se=se)
    return(ret)
  }

# library(MASS)
# test1, consistent
# label_names <- list(X1=c("x1","x2"),X2=c("x3","x4"))
# rawData <- originalData(2,5,1,3,20,label_names = label_names,beta_int=0.1,scale=10,eps=1)
# ridge0 <- lmRidge(Y~X1+X2,data=rawData,lambda = seq(0, 1, .01),label_names=label_names)
# ridge_test <- lm.ridge (Y ~ x1+x2+x3+x4, lambda = seq(0, 1, .01), data=rawData)
# 
# test2, inconsistent at the beginning
# N <- 20
# x1 <- runif(n=N)
# x2 <- runif(n=N)
# x3 <- runif(n=N)
# x4 <- runif(n=N)
# ep <- rnorm(n=N)
# y <- x1 + x2 + ep
# data0 <- as.data.frame(cbind(y,x1,x2,x3,x4))
# ridge1 <- lmRidge(y~X1+X2,data=data0,lambda = seq(0, 100, 1),label_names=label_names)
# ridge2 <- lm.ridge (y ~ x1+x2+x3+x4, lambda = seq(0, 100, 1))


# Step 2: Derive score test for multivariate interaction term. Implement in R. 
# Check with R implementation, see if results agree.
multiScore <- 
  function(formula, data = NULL,label_names = label_names){
    # prepare data and estimate under null
    y0 <- data[,as.character(attr(terms(formula),"variables"))[2]]
    generic_formula0=generic_formula(formula=formula,label_names = label_names)$generic_formula
    len <- generic_formula(formula=formula,label_names = label_names)$length_main_effect
    x <- model.matrix(generic_formula0 ,data=data)
    x_main <- x[,c(1:(len+1))]
    if(abs(det(t(x_main)%*%x_main))<1e-6)
      stop("X must be a non-singluar matrix")
    beta0 <- solve(t(x_main)%*%x_main)%*%t(x_main)%*%y0
    y <- y0-x_main%*%beta0
    se <- sqrt(sum(y^2)/length(y0))
    
    # compute score statistic
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
originalData <- 
  function(b1_l=2, b1_u=5, b2_l=1, b2_u=3, size, label_names = label_names, beta_int=0, scale=1,eps=0.1){
    X1 <- matrix(rnorm(size*length(label_names[[1]]), 0, 1), ncol = length(label_names[[1]]))
    X2 <- matrix(rnorm(size*length(label_names[[2]]), 0, 1), ncol = length(label_names[[2]]))
    beta1 <- runif(length(label_names[[1]]), b1_l, b1_u)*scale
    beta2 <- runif(length(label_names[[2]]), b2_l, b2_u)*scale
    if(beta_int!=0){
      X_int <- NULL
      for (i in 1:length(label_names[[1]])){
        X_int <- cbind(X_int, X1[,i]*X2)
      }
      X <- cbind(1, X1, X2, X_int)
      beta <- t(cbind(2,t(beta1),t(beta2),t(rep(beta_int, length(label_names[[1]])*length(label_names[[2]])))))
    }
    else{
      X <- cbind(1, X1, X2)
      beta <- t(cbind(2,t(beta1),t(beta2)))
    }
    Y <- X%*%beta + rnorm(size, 0, eps)
    data <- as.data.frame(cbind(Y,X1,X2))
    colnames(data) <- c("Y", label_names[[1]],label_names[[2]])
    return(data)
  }
# rawData <- originalData(2,5,1,3,100,label_names = label_names)
# rawData2 <- originalData(2,5,1,3,100,label_names = label_names,beta_int = .1)

# conduct bootstrap for linear regression
linearBoot <- 
  function(formula, data = NULL,label_names = label_names, B=100){
    coef <- linearReg(formula,data=data,label_names=label_names)$beta
    sd <- linearReg(formula,data=data,label_names=label_names)$se
    n <- dim(data)[1]
    inds <- 1:n
    
    # conduct bootstrap
    bs_test <- sapply(1:B, function(k){
      bs_ind <- sample(inds, size = n, replace = T)
      dat <- data[bs_ind,-1]
      dat <- cbind(1,dat)
      Y <- as.matrix(dat)%*%coef + rnorm(n, 0, sd)
      dat <- cbind(Y,dat[,-1])
      multiScore(Y~X1+X2+X1*X2,data=dat,label_names=label_names)$test.stat
    })
    
    # assemble test statistic
    original_test <- multiScore(Y~X1+X2+X1*X2,data=data,label_names=label_names)$test.stat
    bs_pvalue <- sum(as.numeric(original_test)<=bs_test)/B
    return(bs_pvalue)
  }

# conduct bootstrap for ridge regression
ridgeBoot <- 
  function(formula, data = NULL,label_names = label_names,
           lambda = seq(0, 1, .01), B=100){
    coef <- lmRidge(formula,data=data,lambda = lambda,label_names=label_names)$beta
    sd <- lmRidge(formula,data=data,lambda = lambda,label_names=label_names)$se
    n <- dim(data)[1]
    inds <- 1:n
    
    # conduct bootstrap
    bs_test <- sapply(1:B, function(k){
      bs_ind <- sample(inds, size = n, replace = T)
      dat <- data[bs_ind,-1]
      dat <- cbind(1,dat)
      Y <- as.matrix(dat)%*%coef + rnorm(n, 0, sd)
      dat <- cbind(Y,dat[,-1])
      multiScore(Y~X1+X2+X1*X2,data=dat,label_names=label_names)$test.stat
    })
    
    # assemble test statistic
    original_test <- multiScore(Y~X1+X2+X1*X2,data=data,label_names=label_names)$test.stat
    bs_pvalue <- sum(as.numeric(original_test)<=bs_test)/B
    return(bs_pvalue)
  }


# Rerun p-value test in high-signal (i.e. large effect size) 
# low-noise (low variance for residual term, e.g. 0.01) scenario. 
# In order to verify result is correct.
# fix B=1000, M=1000, n=200
verify <- 
  function(i){
    # pvalue <- sapply(1:1000, function(k){
    rawData <- originalData(2,5,1,3,size=n,label_names, beta_int=beta_int,scale=scale,eps=e)
    linear_pvalue <- linearBoot(formula, data = rawData, label_names = label_names, B=B)
    ridge_pvalue <- ridgeBoot(formula, data = rawData, 
                              label_names = label_names, lambda = seq(0, 1, .01),B=B)
    # })
    # prob <- sum(pvalue<0.05)/M
    # return(prob)
    return(c(linear_pvalue,ridge_pvalue))
  }

# verify_linear <- 
#   function(i){
#     # pvalue <- sapply(1:1000, function(k){
#     rawData <- originalData(2,5,1,3,size=n,label_names, beta_int=beta_int,eps=e)
#     pvalue <- linearBoot(formula, data = rawData, 
#                        label_names = label_names, B=B)
#     # })
#     # prob <- sum(pvalue<0.05)/M
#     # return(prob)
#     return(pvalue)
#   }

# library(snowfall)
# sfInit(parallel=T,cpus=20)
# formula <- Y~X1+X2
# label_names <- list(X1=c("x1","x2"),X2=c("x3","x4"))
# result <- NULL
# for(B in c(1000)){
#   for (M in c(1000)){
#     for (n in c(200)){
#       for(beta_int in c(0)){
#         for(b11 in seq(2,10,2)){
#           for (e in seq(0.01,0.09,0.02)) {
#             b12 <- b11+3
#             b21 <- b11-1
#             b22 <- b21+2
#             sfExport("formula","label_names","n","B","M","beta_int",
#                      "b11","b12","b21","b22","e","linearBoot","ridgeBoot",
#                      "originalData","generic_formula","linearReg","lmRidge",
#                      "multiScore")
#             system.time(res <- sfSapply(1:M,verify_linear))
#             write.table(t(res),
#                         file="simulation_pvalue.txt",
#                         row.names=F,col.names=F,append=T)
#             res2 <- apply(res,1,function(x){sum(x<0.05)/M})
#             result <- rbind(result,c(B,M,n,beta_int,b11,b12,b21,b22,e,res2))
#             cat(c(B,M,n,beta_int,b11,b12,b21,b22,e,res2),
#                 file="simulation_type1.txt",append=T,"\n")
#             cat("Finished:B=",B," M=",M,"n=",n,"beta_int=",beta_int,"\n")
#           }
#         }
#       }
#     }
#   }
# }
# write.csv(result,file="simulation_type1.csv",row.names = F,quote = F)

# produce plot similar to the CVEK paper. 
# (x-axis: interaction effect size, y-axis: ratio of p-value smaller than 0.05). 
# vary the interaction effect size between seq(0, 0.1, 0.02), 
# vary the effect size multiplication factor within  (1, 3, 5, 10), 
# vary the noise level within (0.05, 0.1, 0.25, 0.5, 1). 
# Keep bootstrap sample size at 200. Keep b11 = 2 (lower bound of effect size).
# number of repetition M fix at 1000.
# Repeat such experiment for linear and ridge regression.
library(snowfall)
sfInit(parallel=T,cpus=20)
formula <- Y~X1+X2
label_names <- list(X1=c("x1","x2"),X2=c("x3","x4"))
result <- NULL
for (M in c(1000)){
  for (e in c(0.05,0.1,0.25,0.5,1)){
    for(beta_int in seq(0,0.1,0.02)){
      for(scale in c(1,3,5,10)){
        for(n in seq(200,1000,200)){
          for (B in c(200)){
            sfExport("formula","label_names","beta_int",
                     "n","B","M","scale","e","ridgeBoot","linearBoot","originalData",
                     "generic_formula","lmRidge","linearReg","multiScore") 
            system.time(res <- sfSapply(1:M,verify))
            write.table(t(res),file="simulation_power_pvalue.txt",
                        row.names=F,col.names=F,append=T)            
            res2 <- apply(res,1,function(x){sum(x<0.05)/M})
            result <- rbind(result,c(beta_int,c(2,5,1,3)*scale,e,n,res2))
            cat(c(beta_int,c(2,5,1,3)*scale,e,n,res2),
                file="simulation_power.txt",append=T,"\n")
            cat("Finished:effect size=",c(2,5,1,3)*scale,"n=",n,"e=",e,"\n")
          }
        }
      }
    }
  }
}
write.csv(result,file="simulation_power.csv",row.names = F,quote = F)

# plot
file <- read.csv(file="/Users/iristeng/Desktop/simulation_power.csv")
#colnames(file) <- c('beta_int','b11','b12','b21','b22','n','B','prob')
file$effect_ind <- ifelse(file$b11==2,1,
                          ifelse(file$b11==6,2,ifelse(file$b11==10,3,4)))
library(ggplot2)
file$effect_ind <- factor(file$effect_ind, levels = 1:4,
                          labels = c('effect size1','effect size2',
                                     'effect size3','effect size4'))
file$e <- factor(file$e, labels = c('e=0.05','e=0.10','e=0.25','e=0.50','e=1.00'))
file$n <- as.factor(file$n)
p <- ggplot(data = file, aes(x=beta_int, y=linear, shape=n, color=n))+
  geom_point(size=0.8)+geom_line()+
  geom_hline(yintercept=0.05)+facet_grid(effect_ind~e)+
  labs(x='beta interaction', y='probability')+
  theme_set(theme_bw())+theme(panel.grid=element_blank())
p

q <- ggplot(data = file, aes(x=beta_int, y=ridge, shape=n, color=n))+
  geom_point(size=0.8)+geom_line()+
  geom_hline(yintercept=0.05)+facet_grid(effect_ind~e)+
  labs(x='beta interaction', y='probability')+
  theme_set(theme_bw())+theme(panel.grid=element_blank())
q
