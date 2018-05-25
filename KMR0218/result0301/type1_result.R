# type 1 error
#
num <- c(50, 100, 200, 300, 400)
method <- c("linear", "gaussian", "quadratic")
B <- 100
D <- 100
j <- 1
result <- NULL
for(m in method){
  input.address <- "./input/"
  stat <- read.table(paste(input.address, m, "_stat.txt", sep = ""))
  for(i in 1:length(num)){
    original.stat <- as.numeric(stat[(D * (i - 1) + 1):(D * i), 1])
    bs.stat <- stat[(D * (i - 1) + 1):(D * i), 2:(B + 1)]
    perturb.stat <- stat[(D * (i - 1) + 1):(D * i), (B + 2):(2 * B + 1)]
    bs.pvalue <- sapply(1:D, function(b){
      sum(original.stat[b] <= bs.stat[b, ]) / B
    })
    perturb.pvalue <- sapply(1:D, function(b){
      sum(original.stat[b] <= perturb.stat[b, ]) / B
    })
    bs.power <- sum(bs.pvalue < 0.05) / D
    perturb.power <- sum(perturb.pvalue < 0.05) / D
    result <- rbind(result, c(num[i], m, bs.power, perturb.power))
  }
}
colnames(result) <- c("n", "method", "bs.type1", "perturb.type1")
output.address <- "./output/"
write.table(result, file = paste(output.address, "_type1_error.txt"), 
            row.names = F, quote = F)
