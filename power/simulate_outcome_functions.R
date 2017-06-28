###########################################
#simulateCERCovariates takes covariate data,
#uses those data to generate survival times,
#merges those survival times with covariate data


simulateCERCovariates <- function(.covs, .betas, .surv_times = NULL, .cens_times = NULL, .id_var = "id",
                                  .rep = NULL) {
  
  .rep <- as.numeric(.rep)
  
  require(PermAlgo)
  require(dplyr)
  
  if (is.null(.rep) ) {.rep <- 1}
  ################################################
  # load/check covariates
  ################################################
  
  if( is.character(.covs) ) {
    if( !file.exists(.covs) ) {
      stop( paste0(".covs is a character object and assumed to be a path, ", 
                   "but references a file that does not exist\n",
                   ".covs is: ", .covs, "\n") )
    } else {
      .covs_path <- .covs   
      cat( "Attempting to load covariates from ", .covs_path, "\n" )
      .covs <- read.csv(file = .covs_path, stringsAsFactors = FALSE)
      cat( "Covariates loaded.\n" )
    }
  } else if( !is.data.frame( .covs ) ) {
    stop( paste0(".covs is not resolveable. It should be one of the following: \n",
                 "- a path to a .csv file containing named columns of numeric values\n",
                 "- a data.frame\n") )
  }
  
  ################################################
  # load/check betas
  ################################################
  
  if( is.character(.betas) ) {
    if( !file.exists(.betas) ) {
      stop( paste0(".betas is a character object and assumed to be a path, ", 
                   "but references a file that does not exist\n",
                   ".betas is: ", .betas, "\n") )
    } else {
      .betas_path <- .betas
      .betas <- read.csv(file = .betas_path)
      cat( "Betas loaded from ", .covs_path, "\n" )
    }
  } else if( !is.data.frame( .betas_path ) ) {
    stop( paste0(".betas_path is not resolveable. It should be one of the following: \n",
                 "- a path to a .csv file containing named columns of numeric values with one row\n",
                 "- a data.frame\n") )
  }
  
  if( nrow(.betas) != 1 ) {
    stop( paste0(".betas has ", nrow(.betas),  " rows, but should only have one\n")) 
  }
  
  # make sure all variables in .betas are in .covs
  if( length( setdiff( names(.betas), names(.covs) ) ) != 0 ) {
    stop( paste0(".betas contains variables which are not present in .covs\n",
                 "These variables are:\n",
                 paste0(setdiff( names(.betas), names(.covs)), collapse = "\n")))
  }               
  # check for non-numeric variables
  .covs_classes <- colClasses(.covs)
  .covs_classes <- .covs_classes[.covs_classes$col %in% names(.betas), ]
  .not_numeric <- .covs_classes$col[ !grepl("integer|numeric", .covs_classes$class)]
  
  if( length(.not_numeric) != 0 ) {
    stop( paste0(".betas refers to variables in .covs which are not numeric or integer\n",
                 "These variables are:\n",
                 paste0(.not_numeric, collapse = "\n") ) )
  }
  
  # determine how many subjects are in .covs
  if( is.null( .covs[[.id_var]] ) ) { 
    stop( paste0(".id_var ", .id_var, " is NULL\n",
                 ".id_var must uniquely identify subjects in .covs\n" ) )
  } else {
    .n_subjects <- length( unique( .covs[[.id_var]] ) )
  }
  
  # determine how many observations per subject
  .n_obs <- unique( table( .covs[[.id_var]] ) )
  
  if( length(.n_obs) != 1 ) { 
    stop( paste0("Non-uniform follow-up times in .covs\n",
                 "Each subject must have the same number of observations.\n",
                 "Subjects with the following number of observations are present: ", 
                 paste0(.n_obs, collapse = ", "), "\n" ), "\n" )
  }
  
  ################################################
  # validate and format survival times 
  ################################################
  if( is.character(.surv_times) ) {
    if( !file.exists(.surv_times) ) {
      stop( paste0(".surv_time is a character object and assumed to be a path, ", 
                   "but references a file that does not exist\n",
                   ".surv_time is: ", .surv_times) )
    } else {
      .surv_times_vector <- read.csv(file = .surv_times)[ , .rep]
      cat( "Survival times loaded from ", .surv_times, "\n" )
    }
  } else if( is.numeric(.surv_times) ) {
    
    .surv_times_vector <- .surv_times
    
  } else if( is.null(.surv_times) ) {
    cat( paste0( "No survival times specified via .surv_time.\nSurvival times will be drawn from ",
                 "U(0, ", .n_obs, ")\nthen rounded to the nearest integer\n" ) )
    .surv_times_vector <- as.integer( ceiling( runif(.n_subjects, 0, .n_obs) ) )
  } else {
    stop( paste0(".surv_time is not resolveable. It should be one of the following: \n",
                 "- a path to a .csv file containing a single column of numeric values\n",
                 "- a numeric vector\n",
                 "- NULL, which will result in a U(0, ", .n_obs, ") distributions of survival times\n") ) 
  }
  if( !is.integer(.surv_times_vector) ) {
    cat( ".surv_times does not contain integers. Values will be rounded to the nearest integer\n" )
    .surv_times_vector <- as.integer( ceiling(.surv_times_vector) )
  }
  if( min(.surv_times_vector) < 0 ) {
    cat( ".surv_times contains negative values, which will be set to zero\n" )
    .surv_times_vector[.surv_times_vector < 0] <- 0
  }
  if( length(.surv_times_vector) != .n_subjects ) {
    stop( paste0("Data have ", .n_subjects, " subjects, but only ", length(.surv_times_vector),  
                 " survival times were supplied.\n",
                 "Number of survival times should be equal to the number of subjects.\n") )
  }
  ################################################
  # validate and format censoring times 
  ################################################
  if( is.character(.cens_times) ) {
    if( !file.exists(.cens_times) ) {
      stop( paste0(".cens_times is a character object and assumed to be a path, ", 
                   "but references a file that does not exist\n",
                   ".cens_times is: ", .surv_times, "\n") )
    } else {
      .cens_times_vector <- read.csv(file = .cens_times)[ , .rep]
      cat( "Censoring times loaded from ", .cens_times, "\n" )
    }
  } else if( is.numeric(.cens_times) ) {
    
    .cens_times_vector <- .cens_times
    
  } else if( is.null(.cens_times) ) {
    cat( paste0( "No censoring times specified via .cens_times.\nCensoring times will be drawn from ",
                 "U(0, ", .n_obs, ")\n then rounded up to the nearest integer\n" ) )
    .cens_times_vector <- as.integer( ceiling( runif(.n_subjects, 0, .n_obs)) )
  } else if (.cens_times == "none") {
    .cens_times_vector <- rep(.n_obs + 1, .n_subjects)
  } else {
    stop( paste0(".cens_times is not resolveable. It should be one of the following: \n",
                 "- a path to a .csv file containing a single column of numeric values\n",
                 "- a numeric vector\n",
                 "- \"none\", which will produce uncensored data\n",
                 "- NULL, which will result in a U(0, ", .n_obs, ") distribution of censoring times\n") ) 
  }
  if( !is.integer(.cens_times_vector) ) {
    cat( ".cens_times does not contain integers. Values will be rounded up to the nearest integer\n" )
    .cens_times_vector <- ceiling(.cens_times_vector)
  }
  if( min(.cens_times_vector) < 0 ) {
    cat( ".cens_times contains negative values, which will be set to zero\n" )
    .cens_times_vector[.cens_times_vector < 0] <- 0
  }
  if( length(.cens_times_vector) != .n_subjects ) {
    stop( paste0("Data have ", .n_subjects, " subjects, but only ", length(.cens_times_vector),  
                 " censoring times were supplied.\n",
                 "Number of censoring times should be equal to the number of subjects.\n") )
  }
  #generate t0 variables for merging
  .covs$t0 <- rep(0:(.n_obs - 1), .n_subjects)
  
  #generate survival times
  .surv_data <- permalgorithm(numSubjects = .n_subjects, maxTime = .n_obs, 
                              Xmat = as.matrix(.covs[ , names(.betas)]), XmatNames = names(.betas),
                              beta = t(.betas), censorRandom = .cens_times_vector,
                              eventRandom = .surv_times_vector )
  
  names(.surv_data)[c(1,2,4,5)] <- c(.id_var, "d", "t0", "t")
  
  .data_to_merge <- .covs[ ,c(.id_var, "t0", setdiff( names(.covs), names(.surv_data) ) ) ]
  
  .surv_data <- merge(x = .surv_data,
                      y = .data_to_merge,
                      by = c(.id_var, "t0") )
  
  #sort by id and time
  .surv_data <- .surv_data[order( .surv_data[[.id_var]], .surv_data$t0), ]
  
  #order columns
  .surv_data <- .surv_data[ , c(.id_var, "d", "t0", "t", "Fup", 
                                setdiff( names(.surv_data), c(.id_var, "d", "t0", "t", "Fup") ) ) ]
  return(.surv_data)
}

################################################
# colClasses returns the classes of each column 
# in a data frame
################################################

colClasses <- function(.df) {  
  
  .cols <- rep(NA, length(.df[1,]) )
  .class <- rep(NA, length(.df[1,]) )
  
  for (.col in 1:length(.df[1,]) ) {
    .cols[.col] <- names(.df)[.col]
    .class[.col] <- class(.df[[.col]])
  }
  return(data.frame(col = .cols, class = .class, stringsAsFactors = F) )
}
