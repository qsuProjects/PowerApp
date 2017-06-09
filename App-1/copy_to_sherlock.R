# Open terminal first and type:
# kinit ariadnag@stanford.edu (replace with your own id)
# (hit enter) type your password then (hit enter) and do
# ssh ariadnag@sherlock.stanford.edu
# once you've logged to sherlock do the steps below (copy in the terminal window)

cd /share/PI/manishad/power 
mkdir 02_28_2017__12_35_30
 
 # replace the id below with your own user_id and run these lines directly from R
setwd("path/to/your/Dropbpx/download") 
system( "scp -r ./* ariadnag@sherlock.stanford.edu:/share/PI/manishad/power/02_28_2017__12_35_30" )

# double check that all the files got copied, in the terminal copy this line
ls -l ./02_28_2017__12_35_30



# In terminal type this line
rm ./02_28_2017__12_35_30/copy_to_sherlock.R 


# Start running the simulation
sbatch /share/PI/manishad/power/02_28_2017__12_35_30/runCovs.sbatch -p manishad

# If you need to check the status of the simulation you can type in the terminal while logged-in in Sherlock 
sacct

