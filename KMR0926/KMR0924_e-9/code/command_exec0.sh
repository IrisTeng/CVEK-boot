#!/bin/bash
#SBATCH -c 1
#SBATCH -N 1
#SBATCH -t 0-12:00			#job run 12 hour
#SBATCH -p short			#submit to 'short' queue
#SBATCH --mem=917  		# use 4 GB memory
#SBATCH -o wd45_%j.out
#SBATCH -e wd45_%j.err
#SBATCH --mail-type=END      #Type of email notification
#SBATCH --mail-user=wdeng@hsph.harvard.edu
Rscript './simSummary.R'
