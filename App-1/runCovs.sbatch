#!/bin/bash
#################
#set a job name  
#SBATCH --job-name=gen_batch
#################  
#a file for job output, you can check job progress
#SBATCH --output=rm_t.out
#################
# a file for errors from the job
#SBATCH --error=rm_t.err
#################
#time you think you need; default is one hour
#SBATCH --time=48:00:00
#################
#quality of service; think of it as job priority
#SBATCH --qos=normal
#################
#use the QSU job partition
#SBATCH --partition=manishad
#################
#number of nodes you are requesting
#SBATCH --nodes=1
#################
#memory per node; default is 4000 MB
#SBATCH --mem=64000
#you could use --mem-per-cpu; they mean what we are calling cores
#################
#get emailed about job BEGIN, END, and FAIL
#SBATCH --mail-type=ALL
#################
#who to send email to; please change to your email
#SBATCH  --mail-user=ariadnag@stanford.edu
#################
#task to run per node; each node has 16 cores
# CHANGE TO 16 TO USE ALL CORES
#SBATCH --ntasks=1
#################
#SBATCH --cpus-per-task=1
#now run normal batch commands
                
ml load R
srun R -f /share/PI/manishad/power/generatesbatch.R --args 10 10000 150 './betas.csv' '' '' '03_23_2017__11_13_53'