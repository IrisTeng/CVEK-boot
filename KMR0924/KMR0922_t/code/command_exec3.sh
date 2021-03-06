#!/bin/bash
#SBATCH --array=2001-3000	    		#job array list goes 1,2,3...n
#SBATCH -c 1
#SBATCH -N 1
#SBATCH -t 0-12:00			#job run 12 hour
#SBATCH -p short			#submit to 'short' queue
#SBATCH --mem=6370  		# use 4 GB memory
#SBATCH -o %A_%a.out
#SBATCH -e %A_%a.err
#SBATCH --mail-type=END      #Type of email notification
#SBATCH --mail-user=wdeng@hsph.harvard.edu
Rscript './simMaster.R' ${SLURM_ARRAY_TASK_ID}
