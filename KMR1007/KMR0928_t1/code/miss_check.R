setwd('/home/wd45/KMR1002_c/output/Temp')

for (method in c("polynomial", "rbf")) {
  for (l in c(.5, 1, 1.5)) {
    for (d in 1) {
      for (mode in c("AIC", "BIC", "GCV", "AICc", "GCVc", "gmpml", "loocv")) {
        for (strategy in c("avg")) {
          for (int_effect in c(seq(0, .29, .05), seq(.3, 1, .1))) {
            res_key <-
              paste0("true", method, 
                     "_p", 2 * l, 
                     "_l", l,
                     "_mode", mode, 
                     "_strategy", strategy, 
                     "_int", int_effect,
                     "_K", d,
                     ".RData")
            if (!file.exists(res_key)) {
              print(res_key)
            }
          }
        }
      }
    }
  }
}

for (method in c("polynomial", "rbf")) {
  for (l in c(.5, 1, 1.5)) {
    for (d in 2:5) {
      for (mode in c("AIC", "BIC", "GCV", "AICc", "GCVc", "gmpml", "loocv")) {
        for (strategy in c("avg", "exp", "erm")) {
          for (int_effect in c(seq(0, .29, .05), seq(.3, 1, .1))) {
            res_key <-
              paste0("true", method, 
                     "_p", 2 * l, 
                     "_l", l,
                     "_mode", mode, 
                     "_strategy", strategy, 
                     "_int", int_effect,
                     "_K", d,
                     ".RData")
            if (!file.exists(res_key)) {
              print(res_key)
            }
          }
        }
      }
    }
  }
}

for (method in "matern") {
  for (p in 1:2) {
    for (l in c(.5, 1, 1.5)) {
      for (d in 1) {
        for (mode in c("AIC", "BIC", "GCV", "AICc", "GCVc", "gmpml", "loocv")) {
          for (strategy in c("avg")) {
            for (int_effect in c(seq(0, .29, .05), seq(.3, 1, .1))) {
              res_key <-
                paste0("true", method, 
                       "_p", p, 
                       "_l", l,
                       "_mode", mode, 
                       "_strategy", strategy, 
                       "_int", int_effect,
                       "_K", d,
                       ".RData")
              if (!file.exists(res_key)) {
                print(res_key)
              }
            }
          }
        }
      }
    }
  }
}

for (method in "matern") {
  for (p in 1:2) {
    for (l in c(.5, 1, 1.5)) {
      for (d in 2:5) {
        for (mode in c("AIC", "BIC", "GCV", "AICc", "GCVc", "gmpml", "loocv")) {
          for (strategy in c("avg", "exp", "erm")) {
            for (int_effect in c(seq(0, .29, .05), seq(.3, 1, .1))) {
              res_key <-
                paste0("true", method, 
                       "_p", p, 
                       "_l", l,
                       "_mode", mode, 
                       "_strategy", strategy, 
                       "_int", int_effect,
                       "_K", d,
                       ".RData")
              if (!file.exists(res_key)) {
                print(res_key)
              }
            }
          }
        }
      }
    }
  }
}
