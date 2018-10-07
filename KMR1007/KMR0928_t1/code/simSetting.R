# generate command setting table
setting_total0 <- 
  expand.grid(n = 100, 
              method = c("polynomial", "rbf"),
              p = 1:3, 
              int_effect = c(seq(0, .29, .05), seq(.3, 1, .1)),
              mode = c("AIC", "BIC", "GCV", "AICc", "GCVc", "gmpml", "loocv"),
              strategy = "avg",
              beta = "min",
              d = 1)

setting_total1 <- 
  expand.grid(n = 100, 
              method = c("polynomial", "rbf"),
              p = 1:3, 
              int_effect = c(seq(0, .29, .05), seq(.3, 1, .1)),
              mode = c("AIC", "BIC", "GCV", "AICc", "GCVc", "gmpml", "loocv"),
              strategy = c("avg", "erm", "exp"),
              beta = "min",
              d = 2:5)
setting_total1 <- rbind(setting_total0, setting_total1)
setting_total1$l <- setting_total1$p / 2
setting_total1 <- setting_total1[, c(1:3, 9, 4:8)]

setting_total3 <- 
  expand.grid(n = 100, 
              method = c("matern"),
              p = 1:2, l = c(.5, 1, 1.5), 
              int_effect = c(seq(0, .29, .05), seq(.3, 1, .1)),
              mode = c("AIC", "BIC", "GCV", "AICc", "GCVc", "gmpml", "loocv"),
              strategy = "avg",
              beta = "min",
              d = 1)
setting_total2 <- 
  expand.grid(n = 100, 
              method = c("matern"),
              p = 1:2, l = c(.5, 1, 1.5), 
              int_effect = c(seq(0, .29, .05), seq(.3, 1, .1)),
              mode = c("AIC", "BIC", "GCV", "AICc", "GCVc", "gmpml", "loocv"),
              strategy = c("avg", "erm", "exp"),
              beta = "min",
              d = 2:5)
setting_total2 <- rbind(setting_total3, setting_total2)

setting_total <- rbind(setting_total1, setting_total2)

write.table(setting_total, "settings.txt", 
            sep = ",", row.names = FALSE)

# remove previous bash script
file_handle <- "command_exec"
file_names <- 
  grep(file_handle, list.files(), value = TRUE)
file_names_bkp <- 
  paste0("./", file_names)

for (name in c(file_names, file_names_bkp)){
  file.remove(name)
}


# generate cluster bash script
n_run_per_worker <- 1
n_exec <- 
  ceiling(nrow(setting_total) / (n_run_per_worker * 1e3))
cluster_type <- c("SLURM", "LSF")[1]

for (exec_id in 1:n_exec) {
  n_rask_rest <- 
    min(nrow(setting_total)/n_run_per_worker - 
          (exec_id-1)*1e3, 1e3)
  if (cluster_type == "SLURM"){
    string <- 
      paste0(
        "#!/bin/bash
#SBATCH --array=", (exec_id-1) * 1e3 + 1, 
        "-", (exec_id-1) * 1e3 + n_rask_rest, 
        "	    		#job array list goes 1,2,3...n
#SBATCH -c 1
#SBATCH -N 1
#SBATCH -t 0-12:00			#job run 12 hour
#SBATCH -p short			#submit to 'short' queue
#SBATCH --mem=6370  		# use 4 GB memory
#SBATCH -o %A_%a.out
#SBATCH -e %A_%a.err
#SBATCH --mail-type=END      #Type of email notification
#SBATCH --mail-user=wdeng@hsph.harvard.edu
Rscript './simMaster.R' ${SLURM_ARRAY_TASK_ID}"
      )
  } else if (cluster_type == "LSF"){
    string <- 
      paste0(
        "#!/bin/bash
        #BSUB -q short					#submit to 'short' queue
        #BSUB -n 1						#each  job run on 1 core
        #BSUB -W 12:00					#job run 12 hour
        #BSUB -J jobArray[", (exec_id-1) * 1e3 + 1, 
        "-", (exec_id-1) * 1e3 + n_rask_rest, 
        "]		#job array list goes 1,2,3...n
#BSUB -o './log/out_%I.txt' 			#lsf output file
#BSUB -e './log/err_%I.txt' 			#lsf error file
#BSUB -R 'rusage[mem=4096]'		#use 4GB memory
Rscript './simMaster.R' $LSB_JOBINDEX"
      )
  }
  write(string, file = paste0("command_exec", exec_id, ".sh"))
  write(string, file = paste0("./command_exec", exec_id, ".sh"))
}

