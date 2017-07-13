library(shiny)
library(shinyAce)
library(mailR)
library(shinyjs)
library(rdrop2)

shinyServer(
  function(input, output, session) {
    
    
    # HANDLERS TO DOWNLOAD TEMPLATES WHEN THE LINK IS CLICKED 
     output$download_cov <- downloadHandler(
       filename = "effects.csv" ,
       content = function(file) {
         write.csv(
           read.csv("./files/betas.csv"), 
           file, na = "",
           row.names = FALSE,
           col.names = FALSE)})
     
     output$download_param <- downloadHandler(
       filename = "parameters.csv" ,
       content = function(file) {
         write.csv(
           read.csv("./files/parameters.csv"), 
           file, na = "",
           row.names = FALSE,
           col.names = FALSE)})
     
     output$download_catparam <- downloadHandler(
       filename = "categorical_parameters.csv" ,
       content = function(file) {
         write.csv(
           read.csv("./files/categorical_parameters.csv"), 
           file, na = "",
           row.names = FALSE,
           col.names = FALSE)})
     
     output$download_pcor <- downloadHandler(
       filename = "pcor.csv" ,
       content = function(file) {
         write.csv(
           read.csv("./files/pcor.csv"), 
           file, na = "",
           row.names = FALSE,
           col.names = FALSE)})
     
     output$download_wcor <- downloadHandler(
       filename = "wcor.csv" ,
       content = function(file) {
         write.csv(
           read.csv("./files/wcor.csv"), 
           file, na = "",
           row.names = FALSE,
           col.names = FALSE)})
     
     ######################################################
     

     #####################################################
     ########## CONNECT TO DROPBOX TO STORE PERSISTENT FILES ########
     
     # To set up this for the first time do the following steps: 
     # 1- library(rdrop2) 
     # 2- token <- drop_auth() , this will open an authentication window on the browser 
     # 3 - saveRDS(token, "droptoken.rds"), an authentication file will be created that can be store locally in the App-1 folder and
     # be used every time that the application needs to authenticate Dropbox.
     
     # read it back with readRDS
     token <- readRDS("droptoken.rds")
     # Then pass the token to each drop_ function
     drop_acc(dtoken = token)
     
     ###### authentication completed now let's create a function to those files #####
     saveData <- function(fileName, outputDir) {
       drop_upload(fileName, dest = outputDir)
     }
     
     check_input <- function(){
       
       validate(
         need(input$file6, 'Effects file (Box 3) is missing'),
         need(input$file7, 'Covariates file (Box 2) is missing'), 
         #need(input$file1, 'Cat parameters file (Box 2) is missing'),
         need(input$file3, 'PCOR file (Box 2) is missing'),
         need(input$file2, 'WCOR file (Box 2) is missing'), 
         need(input$text, 'An email address must be provided'),
         need(input$sig, 'Significance (Box 1) must be specified'),
         need(input$numsu, 'Number of subjects (Box 1) must be specified'), 
         need(input$numsi, 'Number of simulations (Box 1) must be specified'), 
         need(input$drugs, 'Exposures to be tested (Box 1) must be specified'), 
         need(input$numob, 'Number of observations (Box 2) must be specified')

       )   
       
     }
     
     output$error <- renderText({ NULL })

     # GENERATE SBATCH FILES #############################
     source("runSimulation.R")
     
    
      # DO STUFF after "run simulation" button has been clicked
     observeEvent(input$run, {
       
       # VALIDATE REQUIRED FIELDS
       output$error <- renderPrint({
         check_input()
       })
       if (!is.null(check_input()))
       {
          return(NULL)
        }

       # Name files with current Pacific time   
       Sys.setenv(TZ="America/Los_Angeles")
       outputDir <- format(Sys.time(), "%m_%d_%Y__%H_%M_%S")
       print(outputDir)
        
       ##### SAVE FILES SUBMITTED BY THE USER ################
       
       withProgress(message = 'Saving simulation', value = 1, {
         
        # Read betas
         betas <- input$file6
         write.csv(read.csv(betas$datapath),"betas.csv", na = "",
                   row.names = FALSE,
                   col.names = FALSE) 
         saveData("betas.csv",paste('Sim Data/',outputDir, sep = ""))
        
        
        # Read normal covariates
        cov <- input$file7
        # re-organize data
        param <- read.csv(cov$datapath)
        param <- rbind(param[param$type=="static.binary",], param[param$type!="static.binary",])
        
        # Add needed rows and columns 
        param$name <- as.character(param$name)
        param$type <- as.character(param$type)
        param <- rbind(c("id", "id", NA, NA, NA, NA, NA),param)  # adding id row 
        
        # adding XX.var columns 
        param$across.SD <- as.numeric(param$across.SD)
        param$across.var <- param$across.SD**2
        param$within.sd <- as.numeric(param$within.sd)
        param$within.var <- param$within.sd**2
        
        # add categorical variables if any
        cat_cov <- input$file1
        if (!is.null(cat_cov))
        {
          cat_param <- read.csv(cat_cov$datapath)
          levels <- as.factor(cat_param$level)
          for(i in levels(levels)) #add cat names
          {
            param <- rbind(param, c(i, "cat.static", NA, NA, NA, NA, NA, NA, NA))
          }
          
          write.csv(read.csv(cat_cov$datapath),"categorical_parameters.csv", na = "",
                  row.names = FALSE,
                  col.names = FALSE) 
          saveData("categorical_parameters.csv",paste('Sim Data/',outputDir, sep = ""))
        }
        
        write.csv(param,"parameters.csv", na = "",
                  row.names = FALSE,
                  col.names = FALSE) 
        saveData("parameters.csv",paste('Sim Data/',outputDir, sep = ""))
        
        # Read pcor
        pcor <- input$file3
        write.csv(read.csv(pcor$datapath),"pcor.csv", na = "",
                  row.names = FALSE,
                  col.names = FALSE) 
        saveData("pcor.csv",paste('Sim Data/',outputDir, sep = ""))
        
        # Read wcor
        wcor <- input$file2
        write.csv(read.csv(wcor$datapath),"wcor.csv", na = "",
                  row.names = FALSE,
                  col.names = FALSE) 
        saveData("wcor.csv",paste('Sim Data/',outputDir, sep = ""))
        
        drop_create(paste('Sim Data/', outputDir,"/sbatch_files" ,sep = ""))
        drop_create(paste('Sim Data/', outputDir,"/datasets" ,sep = ""))
        
        
        # Read surv 
        surv <- input$file5
        if (!is.null(surv)) {
          # User uploaded a file yet
          write.csv(read.csv(surv$datapath),"survival.csv", na = "",
                    row.names = FALSE,
                    col.names = FALSE) 
          saveData("survival.csv",paste('Sim Data/',outputDir, sep = ""))
          surv <- paste("./survival.csv", sep = "")
        }
        
        # Read surv 
        cens <- input$file4
        if (!is.null(cens)) {
          # User uploaded a file yet
          write.csv(read.csv(cens$datapath),"censoring.csv", na = "",
                    row.names = FALSE,
                    col.names = FALSE) 
          saveData("censoring.csv",paste('Sim Data/',outputDir, sep = ""))
          cens <- paste("./censoring.csv", sep = "")
        }
        
        # WRITE FILE WITH SIMULATION USER INPUT 
        
        drugs <- gsub(",", "-", input$drugs) # replace commas by hyphen 
        drugs <- gsub(" ", "", drugs) # remove space
        df<-NULL; 
        df <- rbind(df,c(input$sig,input$numsu,input$numob,input$text,input$numsi,drugs))
        colnames(df)<- c('sign','subj','obs', 'email', 'no.sims', 'exp.test')
        write.csv(df, "input_params.csv", row.names=FALSE)
        saveData("input_params.csv",paste('Sim Data/',outputDir, sep = ""))
        
        path = "./"
        runfile_path = paste(path, "runCovs.sbatch", sep="")
        if (!is.na(runfile_path)) {
          outfile_lines <- paste(genSbatchCov(input$numsi,input$numsu,input$numob, paste("./betas.csv", sep = ""), surv, cens, outputDir))
          cat(paste0(outfile_lines, collapse = "\n"), file = runfile_path)
        }
        saveData("runCovs.sbatch",paste('Sim Data/',outputDir, sep = ""))
        
        path = "./"
        runfile_path = paste(path, "open_this_file_first.R", sep="")
        if (!is.na(runfile_path)) {
          outfile_lines <- paste(genRfile(outputDir))
          cat(paste0(outfile_lines, collapse = "\n"), file = runfile_path)
        }
        saveData("open_this_file_first.R",paste('Sim Data/',outputDir, sep = ""))

       })
       

       sender <- "qsu.cer.pcori@gmail.com"
       recipients <- input$text
       send.mail(from = sender,
                 to = recipients,
                 subject="CER simulation submitted",
                 body = "We received your simulation parameters. Your simulation should start 
                 running shortly and you should receive a file with results as soon as it completes",
                 smtp = list(host.name = "smtp.gmail.com", port = 465, 
                             user.name="qsu.cer.pcori@gmail.com", passwd="@@pcfaicerS", ssl=TRUE),
                 authenticate = TRUE,
                 send = TRUE) 
       
       sender <- "qsu.cer.pcori@gmail.com"
       recipients <- "qsu.cer.pcori@gmail.com"
       send.mail(from = sender,
                 to = recipients,
                 subject="New simulation",
                 body = "A simulation was just submitted. Please check Dropbox",
                 smtp = list(host.name = "smtp.gmail.com", port = 465, 
                             user.name="qsu.cer.pcori@gmail.com", passwd="@@pcfaicerS", ssl=TRUE),
                 authenticate = TRUE,
                 send = TRUE) 
       
       text("main", "<div style=\"height:900px;\"><img src = \"logo_all.png\" /><div id=\"header\">
             <h2>\"Power Calculator for Associations in Comparative Effectiveness Research Studies\"</h2>
            </div><div class=\"centered\"><table class=\"thanks_table\"><tr><td><b>Thank you!</b></td></tr><tr><td>You have successfully submitted your simulation. <br>You should receive a confirmation email
            shortly!</td></tr></table></div></div>")
             
     })

  }
)