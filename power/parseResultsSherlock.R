args = commandArgs(trailingOnly = TRUE)
print(args)

# simulation folder
.sim_folder = args[1]


if(FALSE) {
  rm(list = ls())
  source("~/.Rprofile")
}


path_to_folder <- paste("/share/PI/manishad/power/",.sim_folder ,sep="")
setwd(path_to_folder)

list.files(paste("/share/PI/manishad/power/",.sim_folder ,sep=""))
openFilesInDirectory <- function(directory,
                                 match_string,
                                 merge = FALSE,
                                 sep = ",",
                                 na.strings = "NA",
                                 header = T,
                                 fill=T,
                                 skip = 0) {
  require(plyr)
  file_array <-  paste0(directory, "/", list.files(directory)[grep(pattern=match_string, list.files(directory))])
  
  data_list <- llply(file_array, function(file_path, delim_str) {
    cat(file_path, "\n")
    to_return <- read.table(file = file_path, header = header, sep = sep, stringsAsFactors = FALSE, fill=fill, quote="\"", na.strings = na.strings, skip = skip )
    
    if(nrow(to_return) > 0 ) {
      to_return["loaded_file_name"] <- tail(strsplit(file_path, "/")[[1]],1)
      return(to_return)
    } else {
      cat(tail(strsplit(file_path, "/")[[1]],1), " was empty\n")
      warning(tail(strsplit(file_path, "/")[[1]],1), " was empty\n")
      to_return["loaded_file_name"] <- NULL
      return(to_return)
      
    } 
  }, delim_str)
  
  if(merge |  length(file_array) == 1) {
    data_list <- ldply(data_list, identity)
  }
  
  return(data_list)
}

# system("scp kikapp@sherlock:/scratch/PI/manishad/PCORI/power/modelFit/* ~/shexport/PCORI/power/modelFit/", 
       # intern = TRUE, ignore.stdout = FALSE, ignore.stderr = FALSE)

# setwd("~/shexport/PCORI/power/")
library(plyr)
# library(knitr)
library(rmarkdown)
# library(dplyr)
# library(ggplot2)
# library(reshape2)

#load model results files
all_results <- openFilesInDirectory(directory = paste0(path_to_folder, "/datasets"), sep = ",", match_string = "results", 
                     merge = TRUE, header = TRUE)

input_params <- openFilesInDirectory(directory = paste0(path_to_folder), sep = ",", match_string = "input_param", 
                                     merge = TRUE, header = TRUE)

effect_params <- openFilesInDirectory(directory = paste0(path_to_folder), sep = ",", match_string = "^parameters", 
                                     merge = TRUE, header = TRUE)


all_results$coef <- as.numeric(all_results$coef)
all_results$se <- as.numeric(all_results$se)
all_results$p <- as.numeric(all_results$p)
all_results$beta <- as.numeric(all_results$beta)

all_results$lcl <- all_results$coef - 1.96 * all_results$se
all_results$ucl <- all_results$coef + 1.96 * all_results$se
all_results$covered <- as.numeric((all_results$lcl <= all_results$beta) & 
                                    (all_results$ucl >= all_results$beta))

betas <- unique(all_results[ ,c("var", "beta")])
test_vars <- unlist(strsplit(input_params$exp.test, "-") )
non_test_vars <- unique(all_results$var[ !(all_results$var %in% test_vars) ] )

to_test <- all_results[all_results$var %in% test_vars, ]
to_test$rep_comp <- rep(1:input_params$no.sims, each = length(test_vars))
  
sig_level <- input_params$sign

#summaryize results by rep
rep_summary_gen <- ddply(to_test, .(rep_comp), function(.rep_data, .sig_level) {
  
  
  .test_vars <- sort(unique(.rep_data$var))
  .n_to_test <- length(unique(.rep_data$var))
  
  .tor <- data.frame(rep = .rep_data$rep_comp[1])
  .tor[["all significant variables were detected"]] = as.numeric(mean(.rep_data$p <= .sig_level) == 1)
  .tor[["Mean proportin of significant variables detected"]] = mean(.rep_data$p <= .sig_level)
  .tor[["number of significant variables detected"]] = sum(.rep_data$p <= .sig_level)
  
  return(.tor)
 
}, sig_level)

rep_summary_vars <- ddply(to_test, .(rep_comp), function(.rep_data, .sig_level) {
  
  
  .test_vars <- sort(unique(.rep_data$var))
  .n_to_test <- length(unique(.rep_data$var))
  
  .tor <- data.frame(rep = .rep_data$rep_comp[1])
  
  for( .var_num in 1:length(.test_vars) ) {
    ..var_name <- paste0(.test_vars[[.var_num]], " detected")
    
    .tor[[..var_name]] <- as.numeric(sum( .rep_data$p <= .sig_level & .rep_data$var == .test_vars[[.var_num]]) )
    
  }
  
  return(.tor)
  
}, sig_level)


rep_summary_at_ex <- ddply(to_test, .(rep_comp), function(.rep_data, .sig_level) {
  
  
  .test_vars <- sort(unique(.rep_data$var))
  .n_to_test <- length(unique(.rep_data$var))
  
  .tor <- data.frame(rep = .rep_data$rep_comp[1])
  
  for( .var_num in 1:.n_to_test ) {
    
    ..var_name_exactly <- paste0("exactly ", .var_num, " significant effect detected")
    ..var_name_atleast <- paste0("at least ", .var_num, " significant effect detected")
    ..var_name <- paste0(.test_vars[[.var_num]], " detected")
    
    .tor[[..var_name_exactly]] <- as.numeric(sum( .rep_data$p <= .sig_level) == .var_num)
    .tor[[..var_name_atleast]] <- as.numeric(sum( .rep_data$p <= .sig_level) >= .var_num)
    
  }
  
  return(.tor)
  
}, sig_level)


#summary(rep_summary$n_nonsig_detected)

# write.table(x = rep_summary, file = paste0(path_to_folder, "repSummary.csv"), 
            # sep = ",", col.names = TRUE, row.names = FALSE)




output_file <- 
"---
title: 'Power Simulation Preliminary Results'
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = FALSE)
```


## Description of Simulation

These results were generated using the following parameters:

* Data Sets:"


n_sims <- paste0("\n     + ", input_params$no.sims, " input data sets")
n_subj <- paste0("\n     + ", input_params$subj, " subjects each")
n_obs <- paste0("\n     + ", input_params$obs, " observations for each subject")

var_text <- "\n* The following variables were included with the indicated effect sizes:"
for(..var in c(test_vars, non_test_vars)) {
  var_text <- paste0(var_text,
                      "\n     + ", ..var, ": ", betas$beta[betas$var == ..var])  
}

models_text <- "\n* Using the following model(s)"
for(..model in unique(all_results$type)) {
  models_text <- paste0(models_text,
                      "\n     + ", ..model)  
}

gen_summary_text <- paste0("\n* Power Simulation Results:",
                   "\n     + Proportion of reps where all significant variables were detected: ", 
                   round(mean(rep_summary_gen[,3]),3),
                   "\n     + Mean proportion of significant variables detected: ", 
                   round(mean(rep_summary_gen[,4]),3),
                   "\n     + Mean number of significant variables detected: ", 
                   round(mean(rep_summary_gen[,5]),3)
)

var_summary_text <- ""
for(..var in names(rep_summary_vars)[-c(1:2)] ) {
  
  var_summary_text <- paste0(var_summary_text,
         "\n     + Proportion of reps where ", ..var, ": ", round(mean(rep_summary_vars[[..var]]),3) )
  
}

at_ex_test <- ""
for(..var in names(rep_summary_at_ex)[-c(1:2)] ) {
  
  at_ex_test <- paste0(at_ex_test,
                             "\n     + Proportion of reps where ", ..var, ": ", round(mean(rep_summary_at_ex[[..var]]),3) )
  
}

#compile output
output_text <- paste0(output_file, 
                      n_sims, 
                      n_subj, 
                      n_obs, 
                      var_text, 
                      models_text, 
                      gen_summary_text, 
                      var_summary_text, 
                      at_ex_test)


cat(output_text, file = "output.Rmd", append = FALSE)
rmarkdown::render("output.Rmd")
