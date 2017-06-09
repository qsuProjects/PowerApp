# Open terminal first and type:
# kinit ariadnag@stanford.edu (replace with your own id)
# (hit enter) type your password then (hit enter) and do
# ssh ariadnag@sherlock.stanford.edu
# once you've logged to sherlock do the steps below (copy in the terminal window)

cd /share/PI/manishad/power 
mkdir 03_23_2017__11_13_53
 
 # replace the id below with your own user_id and run these lines directly from R
setwd("path/to/your/Dropbpx/download") 
system( "scp -r ./* ariadnag@sherlock.stanford.edu:/share/PI/manishad/power/03_23_2017__11_13_53" )

# double check that all the files got copied, in the terminal copy this line
ls -l ./03_23_2017__11_13_53


# In terminal type this line
rm ./03_23_2017__11_13_53/open_this_file_first.R 


# Start running the simulation
sbatch /share/PI/manishad/power/03_23_2017__11_13_53/runCovs.sbatch -p manishad

# If you need to check the status of the simulation you can type in the terminal while logged-in in Sherlock 
sacct

# Once the simulation is fully completed you can delete the simulation files using
cd /share/PI/manishad/power 
rm -R ./03_23_2017__11_13_53

