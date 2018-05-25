# plot
# bs.stat vs original.stat
num <- c(50, 100, 200, 300, 400)
method <- c("linear", "gaussian", "quadratic")
B <- 100
D <- 100
j <- 1 # jth bs.stat in D copy data
for(m in method){
  output.address <- "./output/"
  jpeg(file = paste(output.address, m, "_original.stat", ".jpeg", sep = ""), 
       width = 1200, height = 800)
  input.address <- "./input/"
  stat <- read.table(paste(input.address, m, "_stat.txt", sep = ""))
  par(mfrow = c(2, 3)) 
  for(i in 1:length(num)){
    original.stat <- as.numeric(stat[(D * (i - 1) + 1):(D * i), 1])
    bs.stat <- as.numeric(stat[j + D * (i - 1), 2:(B + 1)])
    perturb.stat <- as.numeric(stat[j + D * (i - 1), (B + 2):(2 * B + 1)])
    # hist(perturb.stat)
    plot(density(bs.stat), main = paste(m, ": sample size =", num[i]), col = 'black')
    lines(density(original.stat), col = 'blue')
  }
  dev.off()
}

# plot
# bs.stat vs perturb.stat
num <- c(50, 100, 200, 300, 400)
method <- c("linear", "gaussian", "quadratic")
B <- 100
D <- 100
j <- 1
for(m in method){
  output.address <- "./output/"
  jpeg(file = paste(output.address, m, "_perturb.stat", ".jpeg", sep = ""), 
       width = 1200, height = 800)
  input.address <- "./input/"
  stat <- read.table(paste(input.address, m, "_stat.txt", sep = ""))
  par(mfrow = c(2, 3)) 
  for(i in 1:length(num)){
    original.stat <- as.numeric(stat[(D * (i - 1) + 1):(D * i), 1])
    bs.stat <- as.numeric(stat[j + D * (i - 1), 2:(B + 1)])
    perturb.stat <- as.numeric(stat[j + D * (i - 1), (B + 2):(2 * B + 1)])
    # hist(perturb.stat)
    plot(density(bs.stat), main = paste(m, ": sample size =", num[i]), col = 'black')
    lines(density(perturb.stat), col = 'blue')
  }
  dev.off()
}
