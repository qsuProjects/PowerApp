library(shiny)
library(shinyjs)

shinyUI(fluidPage(id="page",
  
  includeCSS("style.css"),
  useShinyjs(),  # Set up shinyjs
  
  mainPanel(id="main",
  fluidRow(
  
  div(
  div(img(src = "logo.png"), style="float:left;"),  
  div(span(id="pcori_img", "Developed with the support of"), img(src="pcori_logo.jpeg", width="35%", height="35%", style="float:right;"), style="width:360px; position:relative; float:right; ")
  ,style ="margin:0; padding:0;")
  ),
  
  fluidRow(
  div(id="header",
  h2(titlePanel("Power Calculator for Associations in Comparative Effectiveness Research Studies"))
 ),
  helpText("The following is a simulation-based power calculator for CER studies where covariates can be correlated and either time-varying or time-invariant, 
           and the outcome is right-censored and represents time to an event. The calculator can evaluate power when testing an association between exposure(s) 
           and outcome for various prevalence(s) of the exposure(s). Hover on the", img(src = "question_mark.png", width="15px", height="15px"), 
            "button for a description of the item and in this", a("help file", href="#",style="font-weight:bold; text-decoration: underline;"), 
            "for a detailed explanation of all the user inputs. An example is provided for the specific values currently displayed along with the template files provided. The help 
             file also includes the output produced for the example.", style="color:#777777; margin-left:25px; margin-right:15px;"),
            
            #href="helpfile.html", target="_blank"
  hr()
  ),
      fluidRow(
        column(3, 
               div(id="col2", p(img(src="one.png", id="img_no", width="25px", height="25px"), h4(a("General"))),
              
               div(id="num2", p(h5("Significance Level", style="float:left;"))), 
               div(id="num", numericInput("sig", label = h5(""), value = 0.05, min=0, max=1, step = 0.01), style="margin-top:-20px;"),
               
               div(id="num2", p(h5("Number of Subjects", style="float:left;"),div(class="tt", span(class="tt_text","Expected number of subjects in the study.")))), 
               div(id="num", numericInput("numsu", label = h5(""), value = 10000, min=2), style="margin-top:-20px;"), 
               
               div(id="num2",p(h5("List of exposures to be tested", style="float:left;"), div(class="tt", style="margin-top:-25px; margin-left:210px;",
                                                                                                                                      span(class="tt_text",style="width:270px;", "Indicate which variables (separated by commas). 
                                                                                                                                           The variable names in this list should be the same as the names of 
                                                                                                                                           the variables in the uploaded file under “List of effect sizes of 
                                                                                                                                           variables included in the final Cox model” in Box 3 where the associations are defined. "
                                                                                                                                           )))),
               
               div(id="num", textInput("drugs", label="", value= "d1,d2"), style="margin-top:-18px;"),
                   
               
               
               div(id="num2", p(h5("Number of simulations", style="float:left;"))), 
               div(id="num", numericInput("numsi", label = h5(""), value = 10, min=1), style="margin-top:-20px;"))
              ),
        column(3, 
               div(id="col1", p(img(src="two.png", id="img_no", width="25px", height="25px"), h4(a("Exposure and Covariates"))),
               
                div(id="num2", p(h5("Maximum number of Observations per Subject", style="float:left;"),div(class="tt", style="margin-top:-25px; margin-left:90px;", span(class="tt_text", 
                                                                                                                                                                         style="width:270px;",
                                                                                                                                                                         "The number of observations generated for each subject. 
                                                                                                                                                                         This is the maximum number of records for a given subject.")))), 
                div(id="num", numericInput("numob", label = h5(""), value = 150, min = 5), style="margin-top:-20px;"),     
                   
               div(id="num",p(h5("Specification of Continuous and Binary covariates", style="float:left;"), div(class="tt", style="margin-top:-25px; margin-left:125px;", 
                                                                                        span(class="tt_text", style="width:270px;",
                                                                                            "Use the template provided to generate your own file defining the attributes 
                                                                                            of the exposures and covariates needed in the simulation. The algorithm requires 
                                                                                            means and standard deviations for continuous variables and proportions for binary 
                                                                                            variables. A description of the various data types that can be generated is provided 
                                                                                            in the help file.")))), 
               div(id="num",
                downloadLink('download_param', 'Download continuous/binary covariates template', class="download_temp"),  
                fileInput("file7", label=""), style="margin-top:-10px;"),
               
               div(id="num",p(h5("Specification of Categorical Covariates", style="float:left;"), div(class="tt", style="margin-top:-23px; margin-left:85px;",
                                                                                     span(class="tt_text", style="width:270px;",
                                                                                          "The algorithm allows for one categorical variable that is defined 
                                                                                          using a multinomial model. Provide a file with parameter effects that 
                                                                                          describe the relationship between each category and all other covariates 
                                                                                          and exposures. See template for an example. "))), style="margin-top:-20px;"),

               div(id="num",
               downloadLink('download_catparam', 'Download categorical covariates template', class="download_temp"),  
               fileInput("file1", label=""), style="margin-top:-10px;"),
               
               div(id="num2", p(h5("Population Correlation", style="float:left;"),div(class="tt", 
                                                                                      span(class="tt_text", style="width:270px;",
                                                                                           "Provide correlations between all continuous and binary variables defined in the 
                                                                                           uploaded file under Specification of continuous and binary covariates. See template 
                                                                                           for example."))), style="margin-top:-20px;"),
               div(id="num2",
                   downloadLink('download_pcor', 'Download pcor template', class="download_temp"),  
                   fileInput("file3", label = ""), style="margin-top:-10px;"),
               
               div(id="num2",p(h5("Within Subject Correlation", style="float:left;"),div(class="tt", 
                                                                                       span(class="tt_text", style="width:270px;",
                                                                                            "Within-Subject correlation: Provide correlations between all continuous and binary 
                                                                                            variables used for the generation of the time-varying observations for each subject. 
                                                                                            See template for example."))), style="margin-top:-20px;"),
               div(id="num2",
               downloadLink('download_wcor', 'Download wcor template', class="download_temp"),  
               fileInput("file2", label = ""), style="margin-top:-10px;"))
               ),
        
        column(3, 
               div(id="col2",p(img(src="three.png", id="img_no", width="25px", height="25px"), h4(a("Time to Event Outcome"))), 
                
                div(id="num",p(h5("List of effect sizes of variables included in final Cox model", style="float:left;"), div(class="tt", style="margin-top:-23px; margin-left:190px;",
                                                                                                                                          span(class="tt_text",style="width:270px;",
                                                                                                                                               "Provide a file with two rows as in the template. 
                                                                                                                                               One row with variable names and a second row with 
                                                                                                                                               effect sizes. Names should be the same  
                                                                                                                                               as in the uploaded file under Specification of Continuous 
                                                                                                                                               and Binary Covariates.")))),
                   
                div(id="num",
                       downloadLink('download_cov', 'Download effect size template', class="download_temp"),  
                       fileInput("file6", label=""), style="margin-top:-10px;"),
                   
                   
                div(id="num2", p(h5("Data Set with Survival Times", style="float:left;"),div(class="tt", 
                                                                            span(class="tt_text", style="width:270px;",
                                                                                 "Provide a data set with the number of columns equal to the number 
                                                                                 of simulations being performed (defined in Box 1) where each column 
                                                                                 contains N survival times and N equals the expected number of subjects 
                                                                                 in the study (defined in Box 1)."))), style="margin-top:-20px;"),
                div(id="num",
                  fileInput("file5", label = ""), style="margin-top:-20px;"),
                
               
               div(id="num2", p(h5("Data Set with Censoring Times", style="float:left;"),div(class="tt", 
                                                                        span(class="tt_text", style="width:270px;",
                                                                             "Provide a data set with the number of columns equal to the number of simulations
                                                                             being performed (defined in Box 1) where each column contains N censoring times 
                                                                             and N equals the expected number of subjects in the study (defined in Box 1).")))),
               div(id="num",
                   fileInput("file4", label = ""), style="margin-top:-20px;")
               )
               ),
        
        column(3, 
            div(id="col2", p(img(src="four.png", id="img_no", width="25px", height="25px"), h4(a("Email"))), 
            div(id="num", helpText("After inputting all required information, enter your email address and click the button below. Once the 
                     simulation is complete, the results will be  emailed to you.")), 
            div(id="num", textInput("text", label=h5("Email Address"), value = 'name@stanford.edu')),
            div(id="num", actionButton("run", "Run Simulation", icon("paper-plane"), 
                                       style="color: #fff; background-color: #880000; border-color: #000000")),
            div(id = "num", textOutput("error"), style="color:#DC143C; margin-top: 20px; text-align: justify; text-justify: inter-word;")
            ))
        
    ),
    
    fluidRow(
    #textOutput("text1"),
    br(), br(), br(),
    p(id="copyright_text", "Copyright Stanford University - 2016 - 2017. All rights reserved.")
    ))
  )
  )

