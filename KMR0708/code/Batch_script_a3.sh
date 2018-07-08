#!/bin/sh
#SBATCH --ntasks-per-node=16
#SBATCH --nodes=1
#SBATCH -p short #Partition to submit to, another option is general
#SBATCH --time=12:00:00
#SBATCH --mem-per-cpu=1024
#SBATCH --job-name=TestJob
#SBATCH --error=TestJob.%J.stdout
#SBATCH --output=TestJob.%J.stderr
#SBATCH --mail-type=END      #Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=wdeng@hsph.harvard.edu  #Email to which notifications will be sent

/home/wd45/bin/R-3.4.3 CMD BATCH /home/wd45/KMR0708/main_a3.R


