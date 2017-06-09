################################## 
# RUN SIMULATION 
#################################

# GENERATE SBATCH TO TRIGGER SIMULATION

genSbatchCov <- function(n.simulations, n.subj, n.obs, .betas, .surv, .cens, .sim_folder)
{
  return(paste(
"#!/bin/bash
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
srun R -f /share/PI/manishad/power/generatesbatch.R --args ",n.simulations," ",n.subj," ",n.obs," '",
.betas,"' '",.surv,"' '",.cens,"' '", .sim_folder,"'",sep="")
  )
}

# GENERATE R FILE TO COPY EVERYTHING TO SHERLOCK

genRfile <- function (folder_name)
{
  return(paste(
"# Open terminal first and type:
# kinit ariadnag@stanford.edu (replace with your own id)
# (hit enter) type your password then (hit enter) and do
# ssh ariadnag@sherlock.stanford.edu
# once you've logged to sherlock do the steps below (copy in the terminal window)

cd /share/PI/manishad/power 
mkdir ",folder_name,  

"\n \n # replace the id below with your own user_id and run these lines directly from R
setwd(\"path/to/your/Dropbpx/download\") 
system( \"scp -r ./* ariadnag@sherlock.stanford.edu:/share/PI/manishad/power/",folder_name ,"\" )

# double check that all the files got copied, in the terminal copy this line
ls -l ./", folder_name, "\n

# In terminal type this line
rm ./",folder_name ,"/open_this_file_first.R \n

# Start running the simulation
sbatch /share/PI/manishad/power/",folder_name,"/runCovs.sbatch -p manishad

# If you need to check the status of the simulation you can type in the terminal while logged-in in Sherlock 
sacct

# Once the simulation is fully completed you can delete the simulation files using
cd /share/PI/manishad/power 
rm -R ./",folder_name,"

", sep=""


))
}