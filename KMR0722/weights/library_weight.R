setwd("/Users/dorabeedeng/Desktop/KMR/KMR0722/weights")
library(dplyr)
boot_weight2 <- read.table("boot_weight_2.txt", header = T)
boot_weight3 <- read.table("boot_weight_3.txt", header = T)
boot_weight4 <- read.table("boot_weight_4.txt", header = T)
boot_weight5 <- read.table("boot_weight_5.txt", header = T)
asym_weight2 <- read.table("asym_weight_2.txt", header = T)
asym_weight3 <- read.table("asym_weight_3.txt", header = T)
asym_weight4 <- read.table("asym_weight_4.txt", header = T)
asym_weight5 <- read.table("asym_weight_5.txt", header = T)

Boot_Weight <- list()
Boot_Weight[[1]] <- boot_weight2 %>% group_by(p, int_effect, mode, strategy) %>% 
  summarise_at(vars(colnames(boot_weight2)[1:3]), funs(mean(., na.rm = TRUE)))

Boot_Weight[[2]] <- boot_weight3 %>% group_by(p, int_effect, mode, strategy) %>% 
  summarise_at(vars(colnames(boot_weight3)[1:3]), funs(mean(., na.rm = TRUE)))

Boot_Weight[[3]] <- boot_weight4 %>% group_by(p, int_effect, mode, strategy) %>% 
  summarise_at(vars(colnames(boot_weight4)[1:6]), funs(mean(., na.rm = TRUE)))

Boot_Weight[[4]] <- boot_weight5 %>% group_by(p, int_effect, mode, strategy) %>% 
  summarise_at(vars(colnames(boot_weight5)[1:6]), funs(mean(., na.rm = TRUE)))


Asym_Weight <- list()
Asym_Weight[[1]] <- asym_weight2 %>% group_by(p, int_effect, mode, strategy) %>% 
  summarise_at(vars(colnames(asym_weight2)[1:3]), funs(mean(., na.rm = TRUE)))

Asym_Weight[[2]] <- asym_weight3 %>% group_by(p, int_effect, mode, strategy) %>% 
  summarise_at(vars(colnames(asym_weight3)[1:3]), funs(mean(., na.rm = TRUE)))

Asym_Weight[[3]] <- asym_weight4 %>% group_by(p, int_effect, mode, strategy) %>% 
  summarise_at(vars(colnames(asym_weight4)[1:6]), funs(mean(., na.rm = TRUE)))

Asym_Weight[[4]] <- asym_weight5 %>% group_by(p, int_effect, mode, strategy) %>% 
  summarise_at(vars(colnames(asym_weight5)[1:6]), funs(mean(., na.rm = TRUE)))

View(Boot_Weight[[1]])
View(Boot_Weight[[2]])
View(Boot_Weight[[3]])
View(Boot_Weight[[4]])

View(Asym_Weight[[1]])
View(Asym_Weight[[2]])
View(Asym_Weight[[3]])
View(Asym_Weight[[4]])
