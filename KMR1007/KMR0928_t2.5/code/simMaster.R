setwd('/home/wd45/KMR0928')
library(magrittr)
library(MASS)
library(survey)
library(CompQuadForm)
library(mvtnorm)
library(dplyr)
library(psych)
library(limSolve)
# library(CVEK)
rm(list = ls())
source('simWorker.R')
source('./source/define_model.R')
source('./source/ensemble.R')
source('./source/estimation.R')
source('./source/generate_kernel.R')
source('./source/testing.R')
source('./source/tuning.R')
source('./source/util.R')

# config read in 
config_all <- read.csv("settings.txt")

# read config
print("Mr Handy: How may be of service, master?")
print("Mr Handy: Oh a 'argument', how wonderful..")
args <- commandArgs(trailingOnly = TRUE)
config_idx <- as.numeric(args)
#eval(parse(text = args))
print(sprintf("Mr Handy: You see the config index '%d'..", config_idx))
print("Mr Handy: and who gets to read all this mumble jumble? Me, that's who...")

# extract config and execute
setTable <- config_all[config_idx, ]

label_names <- list(X1 = c("x1", "x2"), X2 = c("x3", "x4"))
formula <- Y ~ X1 + X2
formula_int <- Y ~ X1 * X2

for (i in 1:nrow(setTable)){
  # extract command
  n <- setTable$n[i]
  method <- setTable$method[i]
  p <- setTable$p[i]
  l <- setTable$l[i]
  int_effect <- setTable$int_effect[i]
  mode <- setTable$mode[i]
  strategy <- setTable$strategy[i]
  beta <- setTable$beta[i]
  d <- setTable$d[i]
  
  # execute command
  taskname <-
    simWorker(M = 200, n = n, B = 200,
      method = method, int_effect = int_effect, 
      Sigma = 0, l = l, p = p, eps = .01, 
      mode = mode, strategy = strategy, 
      beta = beta, test = "boot", d = d,
      lambda = exp(seq(-10, 5, .5)))
  
  # sign out sheet
  write(
    paste(config_idx, taskname, Sys.time(), 
          collapse = "\t\t"), 
    file="sign_out.txt", append=TRUE)
  print("Mr Handy: there, 'sign_out.txt' signed. You hava a --")
}

print("Mr Handy: Here kitty-kitty-kitty....")
