

# load command line arguments
args = commandArgs(trailingOnly = TRUE)
print(args)

# the number of clusters to use
n_clusters = as.numeric(args[1])

# name prefix for datasets
name_prefix = args[2]

# covariates file
.covs  = args[3]

# betas file
.betas = args[4]

# survival file
.surv = args[5]

# censoring times file
.cens = args[6]

# simulation file
.sim_folder = args[7]

# simulation number
n = args[8]

if(.surv == ''){.surv <- NULL}
if(.cens == ''){.cens <- NULL}

print(.covs)
print(.betas)
print(.surv)
print(.cens)


# load required libraries
require(plyr)
require(doSNOW)
require(car)

number_cores <- 16
#set up parallelization
cl<-makeCluster(number_cores)
registerDoSNOW(cl)

number_cores
getDoParWorkers()


#give each worker an ID 
#this is optional, but occasionally is useful 
#in situations where you need workers to write data
#every worker writing to the same file will result in 
#uncaught errors and partially overwritten data
#worker ids enable each worker to write to its own set of 
#files
clusterApply(cl, seq(along=cl), function(id) WORKER.ID <<- paste0("worker_", id))

  
  l_ply( c( 1:n), .parallel=T, function(.item, .covs, .betas, .surv, .cens, name_prefix, .sim_folder, n) {
#  foreach( c( 1:1 ))  %dopar%  function(.covs, .betas, .surv, .cens, name_prefix, .sim_folder){
    #l_ply iterates through each element of its first argument, here c(1:1000), and passes each element to
    #function() as .item    
    #instead of c(1:1000), you could pass l_ply() a list where each element is a set of subject specific means
    #then each worker/iteration would expand those means into a full set of data for that subject
    #avoid workers writing to the same file
    #if each worker/iteration generates data for one subject, you will also need to write code that stitches
    #all subject data together into one data file.
    
    #everything in here will be performed by workers in parallel    
    #packages and source functions used by workers should be loaded within the this block
    
    #generate outcome 
    setwd("/share/PI/manishad/power/")
    source("simulate_outcome_functions.R", local=TRUE)
    
    setwd(paste("/share/PI/manishad/power/", .sim_folder, sep=""))
    outcome <- simulateCERCovariates(.covs, .betas, .surv, .cens, "id", n)
    setwd(paste("/share/PI/manishad/power/",.sim_folder ,"/datasets", sep = ""))
    name_prefix= paste( name_prefix, "outcome", sep="_" )
    write.csv(outcome, name_prefix)
    
  }, .covs, .betas, .surv, .cens, name_prefix, .sim_folder, n)
  






