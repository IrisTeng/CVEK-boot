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
# dat <- as.data.frame(matrix(round(runif(240, -10, 10), digits = 2), ncol = 8))
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
