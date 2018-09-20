# generate command setting table
setting_total <- 
  expand.grid(n = 100, 
              method = "polynomial",
              p = 3, l = 1.5, 
              int_effect = seq(0, 1, .1),
              mode = c("AICc", "GCVc", "gmpml", "loocv"),
              strategy = c("average", "exp", "erm"),
              beta = "min",
              d = 1:5)

write.table(setting_total, "settings.txt", 
            sep = ",", row.names = FALSE)



# setting_total <- read.table("settings.txt", header = T)

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
#SBATCH -e wd45_%j.err
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
