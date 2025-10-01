# mean.n function ----
#similar to mean.n in SPSS
#can generate scale means that are a little more lenient with NAs
#for example if you have a 15-item scale and you have a participant that completed only 14
#of those items, rowMeans() command would yield an NA for that observation
#the "n" in this function should be the number of NAs you will allow
#df is a dataframe of the items in the scale
#Example: mean.n(df,3) will generate a vector that is the mean of each row in the data frame as long as
#there are no more than 3 NAs

mean.n<-function(df,n){
  means <- apply(as.matrix(df), 1, mean, na.rm = TRUE)
  notvalid <- apply(as.matrix(df), 1, function(df) sum(is.na(df)))
  ifelse(notvalid <= n, means, NA)
}

#  Search Function ----
#found this function online
#example:varlist(data,pattern="hsns",exclude="5|oj")
#above example identifies all variable names that include 'hsns' but
#excludes those that have a 5 or "oj". Pipe key "|" is the "or" operator
#other operators (i.e,. "&") can also be used.
varlist <- function (pattern = NULL, df = NULL, exclude = NULL,
                     type = c("numeric", "factor", "character", "double", "logical", "integer"),
                     ignore.case = TRUE) {

  # --- Input Validation ---
  if (is.null(pattern)) {
    stop("Argument 'pattern' must be provided (regular expression string).")
  }
  if (!is.character(pattern) || length(pattern) != 1) {
    stop("Argument 'pattern' must be a single character string (regex).")
  }
  if (!is.null(exclude) && (!is.character(exclude) || length(exclude) != 1) ) {
    stop("Argument 'exclude' must be NULL or a single character string (regex).")
  }
  if (!is.logical(ignore.case) || length(ignore.case) != 1) {
    stop("Argument 'ignore.case' must be TRUE or FALSE.")
  }
  # Allow type to be NULL or missing, meaning default to all types
  if (!missing(type) && !is.null(type) && !is.character(type)) {
    stop("Argument 'type' must be a character vector of valid types or NULL.")
  }


  # --- Dataframe Handling ---
  # If df argument is NULL, try to get 'df' from the calling environment
  if (is.null(df)) {
    if (exists("df", envir = parent.frame())) {
      df <- get("df", envir = parent.frame())
      # Optional: Inform the user that the default 'df' is being used
      # message("Using dataframe 'df' from the calling environment.")
    } else {
      stop("Argument 'df' was not provided and no object named 'df' was found in the calling environment.")
    }
  }

  # Check if the result is actually a dataframe
  if (!is.data.frame(df)) {
    stop("The object assigned to 'df' is not a data frame.")
  }

  # --- Variable Type Filtering ---
  # Get all variable names first
  all_vars <- names(df)
  if (length(all_vars) == 0) {
    return(character(0)) # Handle empty dataframe
  }

  # Filter by type
  vars_filtered_by_type <- character(0)
  valid_types_provided <- FALSE

  # Define the base R type checking functions
  type_checks <- list(
    numeric = is.numeric,   # Includes integer and double
    double = is.double,
    integer = is.integer,
    factor = is.factor,
    character = is.character,
    logical = is.logical
  )

  # Determine which types to actually check based on user input
  # Use intersect to only consider valid types specified in the function definition
  types_to_check <- intersect(tolower(type), names(type_checks))

  if (length(types_to_check) > 0) {
    valid_types_provided <- TRUE
    # Apply the checks for the selected types
    # Note: Using `is.numeric` will catch both integer and double if 'numeric' is in types_to_check
    selected_vars_list <- list()
    for (t in types_to_check) {
      # Apply the corresponding check function to each column
      matches <- sapply(df, type_checks[[t]])
      selected_vars_list[[t]] <- all_vars[matches]
    }
    # Combine results from different type checks and remove duplicates
    vars_filtered_by_type <- unique(unlist(selected_vars_list, use.names = FALSE))

  } else {
    # If 'type' argument was missing, NULL, or empty after validation
    # Default behavior: consider all variables regardless of type
    vars_filtered_by_type <- all_vars
    # Or if type was specified but none were valid after intersect():
    if (!missing(type) && !is.null(type)) {
      warning("None of the specified 'type' values were valid. Considering all variable types. Valid types are: ",
              paste(names(type_checks), collapse=", "))
    }
  }

  # If after type filtering (if any types were specified) we have no variables, return early
  if (valid_types_provided && length(vars_filtered_by_type) == 0) {
    # message("No variables found matching the specified types.") # Optional message
    return(character(0))
  }


  # --- Pattern Matching and Exclusion ---
  # Filter by include pattern
  if (!is.null(pattern)) {
    vars_matched <- vars_filtered_by_type[grepl(pattern, vars_filtered_by_type, ignore.case = ignore.case)]
  } else {
    # Should not happen due to initial check, but for robustness:
    vars_matched <- vars_filtered_by_type
  }


  # Filter by exclude pattern
  if (!is.null(exclude)) {
    # Ensure exclude pattern doesn't remove *all* matches if exclude is too broad
    vars_to_exclude <- grepl(exclude, vars_matched, ignore.case = ignore.case)
    final_list <- vars_matched[!vars_to_exclude]
  } else {
    final_list <- vars_matched
  }

  return(final_list)
}

#More examples of how varlist function works, df = dataframe object:
#can also look up data on the grepl function, all of that should work
## All variable starting with cred:
# varlist(df,pattern="^cred")

### All variables labeled neo followed by any number (e.g., 'neo1','neo2', etc.):
# varlist(df,pattern="neo[1-9]")

## All numeric variables:
# varlist(df,type="numeric")

## All factor variable except variable gb and variables starting with c:
# varlist(df,type="factor",exclude=c("gb|^c"))

## Can use in conjunction with sapply:
# sapply(df[,varlist(df,type="numeric",pattern="credit")], summary)


#####cor function####
#the skeleton of this function I got from a different student here. It gives flagged correlations.
#I modified a little so if you find it doesn't work for one purpose or another it is probably my fault.
#for the most part you can ignore this function because it just feeds into my r_table() function that is below.

cor_table<- function(x,vars=NULL,with=NULL,flag=TRUE,strict=FALSE,round = 2, type = 'pearson'){
  # require(Hmisc,warn.conflicts=TRUE)
  # require(polycor, warn.conflicts = T)
  # I added type = 'mixed' option this to allow for calculating the correlation
  # matrix with a mix of continuous, ordinal, and dichotomous data. Note that
  # the p-value calculation is going to be close but not spot-on because it
  # looks like pearson p-values are typically calculated using a t-test with
  # n-2 df while (according to chatgpt) polychoric correlation p-values are
  # typically based on the z distribution - the type = 'mixed' option currently
  # calculates all p-values using the z-distribution, which will be very close
  # to the t assuming reasonable sample size.
  if(type == 'mixed'){
    data <- data.frame(x)
    Rmat <- polycor::hetcor(data = data, use = 'pairwise.complete.obs')
    zmat <- Rmat$correlations/Rmat$std.errors
    Rmat$r <- Rmat$correlations
    Rmat$P <- p_of_z(zmat)
    RmatP <- Rmat$P
    ComRmat <- round(Rmat$correlations, 5)
    originalRmat<-Rmat$correlations
    ComRPmat <- round(RmatP,3)
    originalRPmat <- RmatP
    ComRnmat<- Rmat$n
  }else{
    data<-as.matrix(x)
    Rmat <- Hmisc::rcorr(data, type = type)
    RmatP <- Rmat$P
    ComRmat <- round(Rmat$r, 5)
    originalRmat<-Rmat$r
    ComRPmat <- round(Rmat$P,3)
    originalRPmat <- Rmat$P # added this line to avoid flagging issues due to rounding errors
    ComRnmat<- Rmat$n
  }

  if(round == 2){
    numformat <- function(val) { sub("^(-?)0.", "\\1.", sprintf("%.2f", val)) } #remove leading zero function
  }

  if(round == 3){
    numformat <- function(val) { sub("^(-?)0.", "\\1.", sprintf("%.3f", val)) } #remove leading zero function
  }

  if(round == 4){
    numformat <- function(val) { sub("^(-?)0.", "\\1.", sprintf("%.4f", val)) } #remove leading zero function
  }

  if(round == 1){
    numformat <- function(val) { sub("^(-?)0.", "\\1.", sprintf("%.f", val)) } #remove leading zero function
  }

  ComRmat<-matrix(numformat(ComRmat),ncol=ncol(ComRmat),nrow=nrow(ComRmat),dimnames = list(rownames(ComRmat),colnames(ComRmat)))

  #this for statement changes the NAs in the p-value table to 1's
  #if else statements are in the form of: ifelse(test, yes, no)
  for(i in 1:nrow(ComRmat)) {
    for (g in 1:ncol(ComRmat)){
      ifelse(is.na(ComRPmat[i,g]), ComRPmat[i,g] <- 1, ComRPmat[i,g] <- ComRPmat[i,g])
    }
  }

  ComRmatFin <- matrix(nrow=nrow(Rmat$r), ncol=ncol(Rmat$r)) #makes an empty matrix of correct size

  if(flag==TRUE){
    if(strict==FALSE){
      for(i in 1:nrow(ComRmat)) {
        for (g in 1:ncol(ComRmat)){
          # in the following statements replaced ComRPmat with originalRPmat
          ifelse(originalRPmat[i,g] <= .01, ComRmatFin[i,g] <- paste(ComRmat[i,g], '**',sep=""),
                 ifelse(originalRPmat[i,g] <= .05, ComRmatFin[i,g] <- paste(ComRmat[i,g], '*',sep=""),
                        ifelse(originalRPmat[i,g] <= .1, ComRmatFin[i,g] <- paste(ComRmat[i,g], "'",sep=""),
                               ifelse(i==g,ComRmatFin[i,g]<-paste('n=',ComRnmat[i,g],sep=""),
                                      ComRmatFin[i,g] <- ComRmat[i,g]))))
        }
      }
    }
    if(strict==TRUE){
      for(i in 1:nrow(ComRmat)) {
        for (g in 1:ncol(ComRmat)){
          # in the following statements replaced ComRPmat with originalRPmat
          ifelse(originalRPmat[i,g] <= .01, ComRmatFin[i,g] <- paste(ComRmat[i,g], '*',sep=""),
                 ifelse(i==g,ComRmatFin[i,g]<-paste('n=',ComRnmat[i,g],sep=""),
                        ComRmatFin[i,g] <- ComRmat[i,g]))
        }
      }
    }
  }

  if(flag==FALSE){

    for(i in 1:nrow(ComRmat)) {
      for (g in 1:ncol(ComRmat)){
        ComRmatFin[i,g] <- ComRmat[i,g]
      }
    }
  }

  ComRmatFin <- as.data.frame(ComRmatFin,stringsAsFactors=FALSE)

  ComRPmat<-as.data.frame(ComRPmat)
  ComRnmat<-as.data.frame(ComRnmat)
  originalRmat<-as.data.frame(originalRmat)

  names(ComRmat)<-colnames(ComRmat)
  row.names(ComRmat)<-row.names(ComRmat)

  names(originalRmat)<-colnames(ComRmat)
  row.names(originalRmat)<-row.names(ComRmat)

  names(ComRmatFin) <- colnames(ComRmat)
  row.names(ComRmatFin) <- row.names(ComRmat)

  names(ComRPmat) <- colnames(ComRmat)
  row.names(ComRPmat) <- row.names(ComRmat)

  names(ComRnmat) <- colnames(ComRmat)
  row.names(ComRnmat) <- row.names(ComRmat)

  if(!is.null(vars)&is.null(with)){
    ComRmatFin<-t(ComRmatFin[vars,colnames(ComRmat)[!colnames(ComRmat)%in%vars],drop=FALSE])
    ComRPmat<-t(ComRPmat[vars,colnames(ComRmat[!colnames(ComRmat)%in%vars]),drop=FALSE])
    ComRnmat<-t(ComRnmat[vars,colnames(ComRmat[!colnames(ComRmat)%in%vars]),drop=FALSE])
    originalRmat<-t(originalRmat[vars,colnames(ComRmat[!colnames(ComRmat)%in%vars]),drop=FALSE])
    ComRmat<-t(ComRmat[vars,colnames(ComRmat[!colnames(ComRmat)%in%vars]),drop=FALSE])
  }
  if(is.null(vars)&!is.null(with)){
    ComRmatFin<-ComRmatFin[colnames(ComRmat)[!colnames(ComRmat) %in% with],with,drop=FALSE]
    ComRPmat<-ComRPmat[colnames(ComRmat)[!colnames(ComRmat) %in% with],with,drop=FALSE]
    ComRnmat<-ComRnmat[colnames(ComRmat)[!colnames(ComRmat) %in% with],with,drop=FALSE]
    originalRmat<-originalRmat[colnames(ComRmat)[!colnames(ComRmat) %in% with],with,drop=FALSE]
    ComRmat<-ComRmat[colnames(ComRmat)[!colnames(ComRmat) %in% with],with,drop=FALSE]

  }
  if(!is.null(vars)&!is.null(with)){
    ComRmatFin<-ComRmatFin[vars,with,drop=FALSE]
    ComRPmat<-ComRPmat[vars,with,drop=FALSE]
    ComRnmat<-ComRnmat[vars,with,drop=FALSE]
    originalRmat<-originalRmat[vars,with,drop=FALSE]
    ComRmat<-ComRmat[vars,with,drop=FALSE]

  }
  list(r_table=ComRmatFin,cors=originalRmat,p_table=ComRPmat,n_table=ComRnmat)

}


#### Package Coding for r_table function ####
#this is the function I use for correlations, it does what the cor_table function does but it is coded like a
#package so it can present a little cleaner output, although I rarely use it for that purpose.

# x has to be a dataframe or matrix (something that can be coerced to matrix)
# Example for use: r_table(test) - includes flags by default
# the diagonal gives the N for each variable

#the "with" command I use for an abbreviated table
#it works the same way as 'with' in the syntax of SPSS if you have ever used that
#basically, if you just want correlate the neo domains with all of your outcomes,
#r_table(df,with=c("n","e","o","a","c")) would give you a table with 5 columns showing the correlations
#of the NEO domains with your outcomes, but not any of the correlations of the neo domains with each other.
# Example: r_table(df,with=c("n","e","o","a","c")) - or - r_table(df, with = names(df[,1:5])) -
# or - r_table(df, with = varlist(df,pattern="hsns"))

# can't really remember how I was intending to use the "vars" command,
# just tried it and it wasn't really what I was expecting
# so I might have fucked that part of the code up

# flag you can use false if you want to remove the flags e.g. r_table(test,with=c("a","b"),flag=FALSE)
# can also just call the original cor_table output:
# r_table(df, with = c('var1','var2','var3'))$cors
# that is actually the way I typically use this function and the way I recommend using it if you are
# going to just be pasting into excel as I typtically am.


r_table<- function(x, #data frame that can be coerced to a matrix
                   vars=NULL, #select the variables to be used
                   with=NULL, #the variables that will appear across the top
                   flag=TRUE, #flag significance
                   strict=FALSE, #If TRUE, it will use .01 cutoff for significance
                   round = 2, #how many decimal places to round to when printing
                   type = 'pearson',
                   ...)
  UseMethod("r_table")

r_table.default<-function(x, #data frame that can be coerced to a matrix
                          vars=NULL, #select the variables to be used
                          with=NULL, #the variables that will appear across the top
                          flag=TRUE, #flag significance
                          strict=FALSE, #If TRUE, it will use .01 cutoff for significance
                          round = 2, #how many decimal places to round to
                          type = 'pearson',
                          ...)
{
  r_table<-cor_table(x,vars=vars,with=with,flag=flag,strict=strict,round=round, type = type)
  class(r_table)<-"r_table"
  r_table
}

print.r_table<-function(x,...){         #removed p-values table from print.r_table
  cat("Correlation Table:\n\n")
  print(x$r_table)
  cat("\nIf strict = FALSE: .01 **; .05 *\nIf strict = TRUE: .01 *\n")
}

summary.r_table<-function(x,...){
  cat("Correlation Table:\n\n")
  print(x$r_table)
  cat("\nIf strict = FALSE: .01 **; .05 *\nIf strict = TRUE: .01 *\n")
}

#### Double Entry ICC ####
# Need to fix to allow input of dataframe with 2 columns
double.icc <- function(v1, v2){
  cor(x = c(v1,v2), y = c(v2,v1))
}

#### Number Format ####
#don't use this one much anymore, found it online, it can be used to remove leading zeros when printing tables
#but it turns them into character vectors so it can lead to problems, I think I have only ever used it within
#the context of the r_table function.
numformat <- function(val) { sub("^(-?)0.", "\\1.", sprintf("%.2f", val)) } #remove leading zero

#### NAs per row ####
#Identify the number of NAs in a row
numNAs<-function(x){
  sum(is.na(x))
}

# Comparing columns, including NAs #####
compareNA <- function(v1,v2) {
  # This function returns TRUE wherever elements are the same, including NA's,
  # and false everywhere else.
  same <- (v1 == v2)  |  (is.na(v1) & is.na(v2))
  same[is.na(same)] <- FALSE
  return(same)
}

#### IRT function ####

irt.func<-function(dataframe, irt.model = NULL, factors = 1, seed = 123, ...){ #dataframe should not include ID variable
  if(apply(dataframe, 1, numNAs) %>% equals(length(dataframe[1,])) %>% any){
    print("Warning, some rows removed for being empty")
    missing.rows.list <- remove.missing.rows(dataframe)
  }
  dataframe <- data.frame(missing.rows.list$df)
  irt<-mirt(dataframe,factors,itemtype = irt.model, technical = list(removeEmptyRows = TRUE)) #using irt
  scores<-fscores(irt,method='EAP', full.scores=TRUE,scores.only=TRUE) #EAP estimation method for the scores
  saved.scores <- fscores(irt, method = 'EAP', full.scores = TRUE, full.scores.SE = TRUE)
  saved.scores <- data.frame(saved.scores)
  if(any(is.na(data.frame(dataframe)))){
    set.seed(seed)
    fulldataframe<-imputeMissing(irt,scores) #just imputing  the data one time
    firt<-mirt(fulldataframe,1,itemtype = irt.model, technical = list(removeEmptyRows = TRUE)) #save imputed dataset
    m2<-M2(firt)#save M2 statistics
    coefs<-coef(firt, simplify=TRUE) #save item parameters
    items<-itemfit(firt,simplify=TRUE) #save item fit statistics
    output<-list("model" = firt, "m2" = m2, "scores" = saved.scores, 'coefs' = coefs, 'itemfit' = items)
  } else{
    m2<-M2(irt)#save M2 statistics
    coefs<-coef(irt, simplify=TRUE) #save item parameters
    items<-itemfit(irt,simplify=TRUE) #save item fit statistics
    output<-list("model" = irt, "m2" = m2, "scores" = saved.scores, 'coefs' = coefs, 'itemfit' = items)
  }
  return(output)
}

##### Copy dataframe to Excel through clipboard #######
write.excel <- function(x,row.names=FALSE,col.names=TRUE,...) {
  write.table(x,"clipboard-10000",sep="\t",na = '', row.names=row.names,col.names=ifelse(row.names == T & col.names == T, NA, col.names),...)
}

#if dataframe is df:
#write.excel(df)
#then ctrl + v into excel spreadsheet

##### Copy data from Excel to R through clipboard #######
read.excel <- function(header=TRUE,...) {
  read.table("clipboard",sep="\t",header=header,...)
}

#after copying to clipboard
# df <- read.excel()

##### Function to automatically load packages i use regularly ######
#.First <- function(x){
#  library(plyr)
#  library(dplyr)
#  library(dtplyr)
#  library(tibble)
#  library(Hmisc)
#  library(psych)
#  library(mirt)
#  library(magrittr)
#}

.First <- function(x){
  library(MASS)
  library(tidyverse)
  library(magrittr)
  library(Hmisc)
  library(psych)
  library(mirt)
  library(jtools)
}

#### copy dataframe to excel on mac ####

write.excel.mac <- function(x,row.names=TRUE,col.names=TRUE,...) {
  clip <- pipe("pbcopy", "w")
  write.table(x,file = clip,sep="\t",row.names=row.names,col.names=col.names,...)
  close(clip)
}

# p-value from z-score ----
p_of_z <- function(x, two_tail = T){
  if(two_tail == T){
    2*(1-pnorm(abs(x)))
  }else{
    1-pnorm(abs(x))
  }
}

##### calculate size of correlation necessary to reach particular level of significance ######

r.crit <- function(tails = 2, p = .01, n){
  t.crit <- qt(1-(p/tails),n-2)
  r.crit <- sqrt(((t.crit)^2)/((n-2)+((t.crit)^2)))
  return(r.crit)
}

#### calculate critical t value #####
#I never use this function, but I found it online and saved it because it is an example of how to implement user imputs in a function

#critical.t <- function(){
#  cat("\n","\bEnter Alpha Level","\n")
#  alpha<-scan(n=1,what = double(0),quiet=T)
#  cat("\n","\b1 Tailed or 2 Tailed:\nEnter either 1 or 2","\n")
#  tt <- scan(n=1,what = double(0),quiet=T)
#  cat("\n","\bEnter Number of Observations","\n")
#  n <- scan(n=1,what = double(0),quiet=T)
#  cat("\n\nCritical Value =",qt(1-(alpha/tt), n-2), "\n")
#}


###### remove entirely missing rows from a dataframe ######

remove.missing.rows <- function(df){
  all.missing.rows <- which(apply(df,1,numNAs) == length(df))
  cat('removed ', length(all.missing.rows), ' observations: ', all.missing.rows, "\n")
  if(length(all.missing.rows) == 0){
    new.df <- df
    cat('returned original dataframe', '\n')
  }
  else{
    new.df <- df[-c(all.missing.rows),]
    cat('returned new dataframe', '\n')
  }
  list(df = new.df, old.df = df, missing.rows = all.missing.rows)
}

#### Paste IRT results into Modfit program - needs to be modified ####
# I think we need to divide by the constant

paste.modfit <- function(data, model.results){
  cat('\n','Recoding data','\n')
  min <- min(data, na.rm = T)
  if(min != 0){
    data <- data - min
    cat('\n','Response values must start at 0','\n')
    cat('\n','Subtracted ',min,' from data','\n')
    cat('\n','Current minimum','\n')
    print(apply(data, 2, min))
    cat('\n','Current maximum','\n')
    print(apply(data, 2, max))
  }
  if(numNAs(data)>0){
    all.missing.rows <- which(apply(data,1,numNAs) == length(data))
    if(length(all.missing.rows) > 0){
      data <- remove.missing.rows(data)
      data <- data$df
    }
    cat('\n','Recoded ',numNAs(data), ' missing data points to 9','\n')
    data[is.na(data)] <- 9
  }
  cat('\n','Number of items: ', length(data),'\n')
  cat('\n','Number of persons: ', nrow(data),'\n')
  coefs <- coef(model.results, simplify=TRUE, IRTpars = T)
  coefs$items[,'a'] <- coefs$items[,'a']/1.702
  write.excel(coefs$items, row.names = F, col.names = F)
  cat("\n","\bPaste item parameters","\n")
  cat("\n","\bCoef a divided by constant 1.702","\n")
  cat("\n","\bEnter 1 when complete","\n")
  continue <- scan(n=1, what = numeric(0), quiet = T)
  if(continue == 1){
    write.excel(data, row.names = F, col.names = F)
    cat("\n","\bPaste item responses","\n")
  }
}


###### write .dat file for GGUM #######
# This function isn't really necessary anymore now that GGUM is in the R IRT package
#write.ggum <- function(dataframe, append=FALSE, quote=FALSE, sep="", na="-9", rownames = FALSE, colnames = FALSE, rowCol = NULL, justify = 'right', formatInfo = TRUE, quoteInfo=TRUE, width = 8, eol="\n", qmethod=c("escape", "double"),  scientific=TRUE, ...) {
#  require(gdata)
#  new.dataframe <- remove.missing.rows(dataframe)
#  dataframe <- new.dataframe$df
#  write.fwf(x = dataframe, file = file.choose(), append = append, quote = quote, sep = sep, na = na, rownames = rownames, colnames = colnames, rowCol = rowCol, justify = justify, formatInfo = formatInfo, quoteInfo=quoteInfo, width = width, eol=eol, qmethod=qmethod,  scientific=scientific, ...)
#}

write.ggum <- function(dataframe, id.var = 'participant',append=FALSE, quote=FALSE, sep="", na="-9", rownames = FALSE, colnames = FALSE, rowCol = NULL, justify = 'right', formatInfo = TRUE, quoteInfo=TRUE, width = 8, eol="\n", qmethod=c("escape", "double"),  scientific=TRUE, ...) {
  require(gdata)
  scale.only <- dataframe[,!names(dataframe) %in% id.var]
  new.dataframe <- remove.missing.rows(scale.only)
  if(length(new.dataframe$missing.rows) > 0){
    data <- data.frame('participant' = dataframe[-c(new.dataframe$missing.rows),id.var], new.dataframe$df)
  }
  if(length(new.dataframe$missing.rows) == 0){
    data <- dataframe
  }
  write.fwf(x = data, file = file.choose(), append = append, quote = quote, sep = sep, na = na, rownames = rownames, colnames = colnames, rowCol = rowCol, justify = justify, formatInfo = formatInfo, quoteInfo=quoteInfo, width = width, eol=eol, qmethod=qmethod,  scientific=scientific, ...)
}


###### read a .csv downloaded from qualtrics ###########
#found the base of this function online

read.qualtrics.csv <- function(filename, stringsAsFactors = FALSE, debrief = FALSE, ...) {
  n <- read.csv(filename, nrows = 1, stringsAsFactors = FALSE)
  dat <- read.csv(filename, header = FALSE, skip = 2, stringsAsFactors = stringsAsFactors, ...)
  names(dat) <- names(n)
  names(dat)[1:10] <- n[1,1:10]
  for(i in seq_along(dat)) {
    attr(dat[,i], "question") <- n[1,i]
  }
  dat <- select(dat, -c(ResponseID,
                        ResponseSet,
                        Name,
                        ExternalDataReference,
                        EmailAddress,
                        IPAddress,
                        Status,
                        Finished,
                        LocationLatitude,
                        LocationLongitude,
                        LocationAccuracy,
                        X))
  if(debrief == TRUE){
    dat <- dat[-c(which(dat$debrief == 2)),]
  }
  dat$StartDate <- strptime(dat$StartDate,format = '%Y-%m-%d %H:%M:%S')
  dat$EndDate <- strptime(dat$EndDate,format = '%Y-%m-%d %H:%M:%S')
  time <- difftime(dat$EndDate,dat$StartDate, units = 'mins')
  participant <- 1:nrow(dat)
  dat <- data.frame('participant' = participant, 'start' = dat$StartDate, 'finish' = dat$EndDate, 'time.minutes' = time, select(dat, -c(StartDate,EndDate)))
  dat
}


### Compute CFI from fa() output #####
#retreived from https://gist.github.com/tonosan/cb7581f3459ae7c4217a
#formula seems to match that provided by David A. Kenny on his website
fa.CFI<-function(x){
  nombre<-paste(x,"CFI",sep = ".")
  nombre<-
    ((x$null.chisq-x$null.dof)-(x$STATISTIC-x$dof))/(x$null.chisq-x$null.dof)
  return(nombre)
}

### test of skew, I don't think this is right, see Crawley 2012 R book in Zotero and http://www.real-statistics.com/tests-normality-and-symmetry/analysis-skewness-kurtosis/ ####
#test.skew <- function(x){
#  x <- x[!is.na(x)]
#  m3 <- sum((x-mean(x, na.rm = T))^3)/length(x)
#  s3 <- sqrt(var(x))^3
#  skew <- m3/s3
#  skew.se <- sqrt(6/length(x))
#  t.test <- skew/skew.se
#  p <- 1-pt(t.test,length(x)-2)
#  return(list('skew' = skew, 't.test' = t.test, 'p' = p))
#}

#### Test correlation differences ######
# functional supply loop for testing correlation differences
# requires output from rcorr() function - labeled in below syntax as "cors"
# rcorr() yields list of 3 dataframes, cors$r, cors$n, cors$p
# this loop assumes you are comparing cor differences across the first four variables in the correlation matrix

#sapply(rownames(cors$r)[5:length(rownames(cors$r))], function(third.var){
#  lapply(1:3, function(x){
#    lapply(1:(4-x), function(y){
#      out <- cocor.dep.groups.overlap(cors$r[third.var,x], cors$r[third.var,(x+y)], cors$r[x,(x+y)], min(cors$n[third.var,(x+y)], cors$n[third.var,(x+y)], cors$n[x,(x+y)]), test = 'meng1992')@meng1992$p.value
#      names(out) = paste0(x,'v',(x+y))
#      return(out)
#    }) %>% do.call(c, .)
#  }) %>% do.call(c, .)
#}) %>% t %>% round(3)

##### comparing nomological nets by correlating correlation profiles ####
profile_correlation<-function(x){
  require(Hmisc)
  a<-matrix(as.numeric(as.matrix(x)),nrow(x),ncol(x),dimnames = list(rownames(x),colnames(x)))
  mat<-matrix(NA,ncol(x),ncol(x),dimnames = list(colnames(x),colnames(x)))
  for(i in 1:ncol(x)){
    for (e in 1:ncol(x)){
      mat[i,e]<-rcorr(a[,i],a[,e])$r[1,2]
    }
  }
  return(mat)
}
#Example syntax:
# cors<-r_table(agreevarsdf,with = varlist(agreevarsdf,pattern = 'F[1-8]'), flag = FALSE)$cors
# profile_r<-profile_correlation(cors) # this is for the supplementary material

###### rIP function - not ready for use ########
#https://github.com/MAHDLab/rIP
#This is an R code project for detecting likely responsese from server farms on MTurk surveys.
#
#Takes as its input an array of IPs and the user's X-Key, passes these to iphub.info, and returns a dataframe with the ip (used for merging),
#country code, country name, asn, isp, block, and hostname.
#
#Especially important in this is the variable "block", which gives a score indicating whether the IP address is likely from a server
#farm and should be excluded from the data. It is codes 0 if the IP is residential/unclassified (i.e. safe IP),
#1 if the IP is non-residential IP (hostping provider, proxy, etc. -- should likely be excluded), and 2 for non-residential and residential IPs
#(more stringent, may flag innocent respondents).
#
#The recommendation from iphub.info is to block or exclude those who score block = 1.
#
#Credit to @tylerburleigh for pointing out the utility of iphub.info. His method for incorporating this information into Qualtrics
#surveys can be found here: https://twitter.com/tylerburleigh/status/1042528912511848448?s=19.

# @export

#getIPinfo <- function(ips, key) {
#  options(stringsAsFactors = FALSE)
#  url <- "http://v2.api.iphub.info/ip/"
#  pb <- txtProgressBar(min = 0, max = length(ips), style = 3)
#  ipDF <- c()
#  for (i in 1:length(ips)) {
#    ipInfo <- httr::GET(paste0(url, ips[i]), httr::add_headers(`X-Key` = key))
#    infoVector <- unlist(httr::content(ipInfo))
#    ipDF <- rbind(ipDF, infoVector)
#    setTxtProgressBar(pb, i)
#  }
#  close(pb)
#  ipDF <- data.frame(ipDF)
#  rownames(ipDF) <- NULL
#  return(ipDF)
#}

### Example code for recoding variable names using find and replace #####
# https://dplyr.tidyverse.org/reference/rename.html - more examples of below renaming function
#df <- rename_with(df, ~ gsub('_','',.x, fixed = T), #rename using a function, replace "_" with ""
#                  starts_with('srp'), ends_with('_r')) #select the variables to apply the function to.

### Example code for renaming variables inside a pipe chain ####
# [dataframe/matrix] %>% 'colnames<-'(character.vector.of.names)

#### trim leading and trailing spaces from character vectors ######
trim <- function (x) gsub("^\\s+|\\s+$", "", x)

#### Paste IRT results into Modfit program - needs to be modified ####
# Example of function that takes input in console
# I think we need to divide by the constant

paste.modfit <- function(data, model.results){
  cat('\n','Recoding data','\n')
  min <- min(data, na.rm = T)
  if(min != 0){
    data <- data - min
    cat('\n','Response values must start at 0','\n')
    cat('\n','Subtracted ',min,' from data','\n')
    cat('\n','Current minimum','\n')
    print(apply(data, 2, min))
    cat('\n','Current maximum','\n')
    print(apply(data, 2, max))
  }
  if(numNAs(data)>0){
    all.missing.rows <- which(apply(data,1,numNAs) == length(data))
    if(length(all.missing.rows) > 0){
      data <- remove.missing.rows(data)
      data <- data$df
    }
    cat('\n','Recoded ',numNAs(data), ' missing data points to 9','\n')
    data[is.na(data)] <- 9
  }
  cat('\n','Number of items: ', length(data),'\n')
  cat('\n','Number of persons: ', nrow(data),'\n')
  coefs <- coef(model.results, simplify=TRUE, IRTpars = T)
  coefs$items[,'a'] <- coefs$items[,'a']/1.702
  write.excel(coefs$items, row.names = F, col.names = F)
  cat("\n","\bPaste item parameters","\n")
  cat("\n","\bCoef a divided by constant 1.702","\n")
  cat("\n","\bEnter 1 when complete","\n")
  continue <- scan(n=1, what = numeric(0), quiet = T)
  if(continue == 1){
    write.excel(data, row.names = F, col.names = F)
    cat("\n","\bPaste item responses","\n")
  }
}

#### Duplicate Item Removal #####

# sample data
#library(readr)
#library(psych)
#library(dplyr)
#full.df <- read_csv("C:/Users/micha/Google Drive/A UGA Stuff/Research Stuff/Structure of Narcissism/clean nar structure dataset no worker id - 5-28-18.csv")
#full.df <- select(full.df, nar1_1:nar8_32)

# Function for identify correlations greater than cutoff to be used in make.cors_df
# Returns correlation matrix with only variables with at least one correlation greater than the cutoff
# Upper diag of the correlation matrix is NA
# Recently modified to be able to take correlation matrix as input as well
make.dat_cors <- function(full.df, cut, cor, dat.is.cor = F){
  require(magrittr)
  require(qgraph)
  if(dat.is.cor == F){
    if(cor == 'cor'){
      dat_cors <- cor(full.df, use = 'pairwise')
    }
    if(cor == 'tet'){
      dat_cors <- tetrachoric(clinical.df)$rho
    }
    if(cor == 'auto'){
      dat_cors <- cor_auto(full.df, missing = 'pairwise')
    }
  }
  if(dat.is.cor == T){
    dat_cors <- as.matrix(full.df)
  }
  dat_cors[upper.tri(dat_cors, diag = TRUE)] <- NA
  dat_cors <- data.frame(dat_cors)
  dat_cors[,
           sapply(dat_cors,function(x){which(x>= cut)}, simplify = TRUE) %>%
             sapply(length) %>%
             equals(0) %>%
             which() %>%
             names %>% c()
  ] <- NULL
  return(dat_cors)
}

# Uses output from make.dat_cors function
# Generates list of length ncol(dat_cors) to be used in make.cors_df
# Each element matches the variable with those that it overlaps with
overlap.var.list <- function(dat_cors, cut){
  #var.list <- sapply(dat_cors,function(x){which(x>= cut)}, simplify = TRUE)
  var.list <- lapply(dat_cors,function(x){which(x>= cut)})
  return(var.list)
}

#overlap.named.var.list <- function(dat_cors, var.list){
#  vars.with.overlap<-list()
#  for(i in 1:length(var.list)){
#    vars.with.overlap[[i]]<-row.names(dat_cors[var.list[[i]],])
#  }
#  names(vars.with.overlap)<-names(var.list)
#  return(vars.with.overlap)
#}

# Uses list generated from overlap.var.list and labels elements within each list
# To be used in make.cors_df
overlap.named.var.list <- function(dat_cors, var.list){
  vars.with.overlap<-list()
  if(length(var.list) == 1){
    vars.with.overlap[[1]]<-row.names(dat_cors)[var.list[[1]]]
  } else{
    for(i in 1:length(var.list)){
      vars.with.overlap[[i]]<-row.names(dat_cors[var.list[[i]],])
    }
  }
  names(vars.with.overlap)<-names(var.list)
  return(vars.with.overlap)
}

# Uses make.dat_cors, overlap.var.list, overlap.named.var.list
# Input is dataframe, cutoff, correlation type, and (optionally) dataframe with two columns matching variable name (column 1) with item content (column 2)
# Generates dataframe with one row for each correlation with magnitude greater than cutoff
make.cors_df <- function(full.df, cut, item.content = NULL, cor = 'cor', dat.is.cor = F){
  dat_cors <- make.dat_cors(full.df, cut, cor, dat.is.cor)
  var.list <- overlap.var.list(dat_cors, cut)
  if(length(var.list) == 0){
    cat('No correlations greater than or equal to ', cut)
  }
  if(length(var.list) > 0){
    named.var.list <- overlap.named.var.list(dat_cors, var.list)
    overlapping.vars<-data.frame('var1' = rep(names(named.var.list),c(sapply(named.var.list,length))),'var2' = unlist(named.var.list))
    char.overlapping.vars <- data.frame('var1' = as.character(overlapping.vars[[1]]), 'var2' = as.character(overlapping.vars[[2]]), stringsAsFactors = FALSE)
    overlapping.vars$var1_num <- rep(which(names(full.df) %in% names(var.list)),c(sapply(var.list,length)))
    overlapping.vars$var2_num <- unlist(var.list)
    overlapping.vars$r <- NA
    for(i in 1:length(overlapping.vars[[1]])){
      overlapping.vars$r[i]<-dat_cors[char.overlapping.vars[i,2],char.overlapping.vars[i,1]]
    }
    if(!is.null(item.content)){
      item.content$var1 <- item.content[,1]
      item.content$var2 <- item.content[,1]
      overlapping.vars <- left_join(overlapping.vars,data.frame('content_v1' = item.content[,2], 'var1' = item.content$var1), by = 'var1')
      overlapping.vars <- left_join(overlapping.vars,data.frame('content_v2' = item.content[,2], 'var2' = item.content$var2), by = 'var2')
    }
    cors_df <- overlapping.vars
    return(cors_df)
  }
}

# To be used in remove.duplicate.items
# Uses output from make.cors_df
# Excludes given item from dataframe
pull.item <- function(cors_df, item){
  #if(exists('removed.items') == TRUE){
  #  removed.items <<- c(removed.items, item)
  #}
  #if(get0('removed.items', ifnotfound = FALSE) == FALSE){
  #  removed.items <<- c(item)
  #}
  #items <- c(item)
  cors_df <- filter(cors_df, var1 != item)
  cors_df <- filter(cors_df, var2 != item)
  return(list('cors_df' = cors_df, 'item' = item))
}

# To be used in remove.duplicate.items
# Calculates total number of correlations each item has in make.cors_df output
calculate.count <- function(cors_df){
  temp.1.count <-table(cors_df$var1) %>% data.frame(stringsAsFactors = FALSE)
  temp.2.count <- table(cors_df$var2) %>% data.frame(stringsAsFactors = FALSE)
  names(temp.1.count) <- c('var','freq1')
  names(temp.2.count) <- c('var','freq2')
  temp.cor.count <- full_join(temp.1.count, temp.2.count, by = 'var')
  temp.cor.count$var <- as.character(temp.cor.count$var)
  temp.cor.count <- tibble(temp.cor.count)
  temp.cor.count$freq1[is.na(temp.cor.count$freq1)] <- 0
  temp.cor.count$freq2[is.na(temp.cor.count$freq2)] <- 0
  temp.cor.count <- mutate(temp.cor.count, total = freq1 + freq2)
  return(arrange(temp.cor.count, desc(total)))
}

# Uses pull.item and calculate.count
# Selects items for removal iteratively until there are no more correlations in the make.cors_df table
# Items with the greatest number of correlations are prioritized, random selection is used when ties occur
# Provides vector with names of variables to be removed
remove.duplicate.items <- function(cors_df, seed = 123){
  temp.cors_df <- cors_df
  removed.items <- c()
  set.seed(seed)
  while(length(temp.cors_df[[1]]) > 0){
    calculate.count(temp.cors_df) %>% top_n(1, total) %>% print()
    item <- calculate.count(temp.cors_df) %>% top_n(1, total) %>% sample_n(1) %>% .[[1]]
    cat('\n','remove ', item,'\n')
    temp.cors_df <- pull.item(temp.cors_df, item)$cors_df
    removed.items <- c(removed.items, item)
  }
  return(removed.items)
}

### Bass-ackward Syntax ####
factor.analyses <- function(df, max.factors = 9, rotate = 'promax', fm = 'pa', alpha = .05, ...){ #can pass other inputs (i.e., rotation, estimation method) to fa function
  require(psych)
  lapply(1:max.factors,function(x, df, rotate = rotate, fm = fm, alpha = alpha, ...){
    fa(r = df, nfactors = x, rotate = rotate, fm = fm, alpha = alpha, ...)
  }, df, rotate = rotate, fm = fm, alpha = alpha, ...)
}

extract.structures <- function(fa.list, order.names = F){
  structure.list <- lapply(fa.list, function(x){unclass(x$Structure)})
  if(order.names == T){
    order.list <- lapply(structure.list, function(x){
      data.frame(x) %>% names() %>% gsub('PA','',., fixed = T) %>% as.numeric() %>% order()
    })
    structure.list <- lapply(1:length(structure.list), function(x){
      structure.list[[x]] <- structure.list[[x]][,c(order.list[[x]])]
    })
  }
  structure.dat <- do.call(data.frame, structure.list)
  suffix <- lapply(1:length(fa.list), function(x){
    1:x
  }) %>% do.call(c, .)
  prefix <- lapply(1:length(fa.list),function(x){
    rep(x,x)
  }) %>% do.call(c,.)
  f.names <- paste0('F',prefix,'.',suffix)
  names(structure.dat) <- f.names
  return(structure.dat)
}

extract.loadings <- function(fa.list, order.names = F){
  loadings.list <- lapply(fa.list, function(x){unclass(x$loadings)})
  if(order.names == T){
    order.list <- lapply(loadings.list, function(x){
      data.frame(x) %>% names() %>% gsub('PA','',., fixed = T) %>% as.numeric() %>% order()
    })
    loadings.list <- lapply(1:length(loadings.list), function(x){
      loadings.list[[x]] <- loadings.list[[x]][,c(order.list[[x]])]
    })
  }
  loadings.dat <- do.call(data.frame, loadings.list)
  suffix <- lapply(1:length(fa.list), function(x){
    1:x
  }) %>% do.call(c, .)
  prefix <- lapply(1:length(fa.list),function(x){
    rep(x,x)
  }) %>% do.call(c,.)
  f.names <- paste0('F',prefix,'.',suffix)
  names(loadings.dat) <- f.names
  return(loadings.dat)
}

extract.vaccounted <- function(fa.list, order.names = F){
  vaccounted.list <- lapply(fa.list, function(x){unclass(x$Vaccounted)})
  vaccounted.list[[1]] <- rbind(vaccounted.list[[1]], 'Cumulative Var' = vaccounted.list[[1]][2,], 'Proportion Explained' = 1, 'Cumulative Proportion' = 1)
  if(order.names == T){
    order.list <- lapply(vaccounted.list, function(x){
      data.frame(x) %>% names() %>% gsub('PA','',., fixed = T) %>% as.numeric() %>% order()
    })
    vaccounted.list <- lapply(1:length(vaccounted.list), function(x){
      vaccounted.list[[x]] <- vaccounted.list[[x]][,c(order.list[[x]])]
    })
  }
  vaccounted.dat <- do.call(data.frame, vaccounted.list)
  suffix <- lapply(1:length(fa.list), function(x){
    1:x
  }) %>% do.call(c, .)
  prefix <- lapply(1:length(fa.list),function(x){
    rep(x,x)
  }) %>% do.call(c,.)
  f.names <- paste0('F',prefix,'.',suffix)
  names(vaccounted.dat) <- f.names

  vec <- lapply(fa.list, function(x){x$Vaccounted['Proportion Var',] %>% sum()}) %>%
    do.call(rbind, .)

  v.list <- list('Vac' = data.frame('Vaccounted' = vec), 'full' = vaccounted.dat)
  return(v.list)
}

#extract.vaccounted <- function(fa.list){
#  vec <- lapply(fa.list, function(x){x$Vaccounted['Proportion Var',] %>% sum()}) %>%
#    do.call(rbind, .)
#  return(data.frame('Vaccounted' = vec))
#}

extract.scores <- function(fa.list, order.names = F, method = 'Thurstone', dat = NULL){
  if(method == 'Thurstone'){
    scores.list <- lapply(fa.list, function(x){x$scores})
  }
  if(method == 'tenBerge'){
    if(is.null(dat)){
      print('Need to include raw data to calculate tenBerge scores')
    }
    scores.list <- lapply(fa.list, function(x){
      factor.scores(dat, x, method = 'tenBerge')$scores
    })
  }
  if(order.names == T){
    order.list <- lapply(scores.list, function(x){
      data.frame(x) %>% names() %>% gsub('PA','',., fixed = T) %>% as.numeric() %>% order()
    })
    scores.list <- lapply(1:length(scores.list), function(x){
      scores.list[[x]] <- scores.list[[x]][,c(order.list[[x]])]
    })
  }
  scores.dat <- do.call(data.frame, scores.list)
  suffix <- lapply(1:length(fa.list), function(x){
    1:x
  }) %>% do.call(c, .)
  prefix <- lapply(1:length(fa.list),function(x){
    rep(x,x)
  }) %>% do.call(c,.)
  f.names <- paste0('F',prefix,'.',suffix)
  names(scores.dat) <- f.names
  return(scores.dat)
}

# This is the original function designed to work with dataframes that were entirely numeric
# avg.list <- function(list){
#   Reduce("+", list)/length(list)
# }

# This was developed with the help of generative AI, to be used as part of a different project
# (needed a way to summarize across a list of lavaan outputs that include both numeric and character vectors)
# I think this should work just as well with lists that are all numeric, but it has not been checked for that purpose
# avg.list <- function(list){
#
#   # Extract column names
#   cols <- colnames(list[[1]])
#
#   # Split between character and numeric columns
#   char_cols <- cols[sapply(list[[1]], is.character)]
#   num_cols <- setdiff(cols, char_cols)
#
#   # Initialize output
#   out <- list[[1]][char_cols]
#
#   # Average numeric columns
#   for(c in num_cols){
#     out[[c]] <- Reduce("+", lapply(list, "[[", c))/length(list)
#   }
#
#   as.data.frame(out)
# }

# Realized the above requires there to be named columns, what follows is a integration of the two functions:
avg.list <- function(list){

  # Extract column names
  cols <- colnames(list[[1]])

  if(is.null(cols)){
    Reduce("+", list)/length(list)
  } else {
    # Split between character and numeric columns
    char_cols <- cols[sapply(list[[1]], is.character)]
    num_cols <- setdiff(cols, char_cols)

    # Initialize output
    out <- list[[1]][char_cols]

    # Average numeric columns
    for(c in num_cols){
      out[[c]] <- Reduce("+", lapply(list, "[[", c))/length(list)
    }

    as.data.frame(out)
  }
}

#Can pull specific parts of EFA output using sapply
#sapply(fa.list, function(x){x$RMSEA})

fit.table <- function(fa.list, vss.result = NULL, df = NULL){
  fit <- data.frame('Factors' = 1:length(fa.list))
  if(is.null(vss.result) & is.null(df)){
    print('Either vss object or raw data matrix must be provided for MAP analysis')
    map <- NA
  }
  if(!is.null(vss.result)){
    if(!is.null(df)){
      print('raw dataframe ignored as vss.results were provided')
    }
    map <- vss.result$map[1:length(fa.list)]
  }
  if(is.null(vss.result) & !is.null(df)){
    vss <- vss(df, n = length(fa.list),
               rotate = fa.list[[2]]$rotation,
               fm = fa.list[[2]]$fm,
               plot = F)
    map <- vss$map[1:length(fa.list)]
  }
  fit$VarAccounted <- sapply(fa.list, function(x){x$Vaccounted['Proportion Var',] %>% sum()})
  fit$Delta.Var.Account <- c(NA,
                             sapply(2:length(fit$VarAccounted), function(x){
                               fit$VarAccounted[x] - fit$VarAccounted[x-1]
                             }))
  rmsea <- t(sapply(fa.list, function(x){x$RMSEA}))
  confidence <- rmsea[1,'confidence']
  cat('RMSEA confidence = ', confidence, '\n')
  fit$RMSEA <- rmsea[,'RMSEA']
  fit$lower.CI <- rmsea[,'lower']
  fit$upper.CI <- rmsea[,'upper']
  fit$CI.overlap <- c(sapply(1:(length(fa.list)-1), function(x){
    ifelse((fit$lower.CI[x]-fit$upper.CI[x+1]) > 0, 0, 1)
  }),NA)
  fit$MAP <- map
  fit$BIC <- sapply(fa.list, function(x){x$BIC})
  return(fit)
}

### Would love to add HTML tabling function #####
# See below for possible resources
# https://www.htmlwidgets.org/showcase_datatables.html
# https://rstudio.github.io/DT/

### Percent of Missing data #####
# identify percent of missing data for each variable
p.missing <- function(df, sort = F){
  p <- unlist(lapply(df, function(x) sum(is.na(x))))/nrow(df)
  if(sort == T){
    p.missing <- sort(p, decreasing = T)
  } else{
    p.missing <- p
  }
  return(p.missing)
}

### Parallel Analysis with 95th %ile of Simulated Data ######
# Need to convert to a function
#set.seed(123)
#parallel.analysis <- fa.parallel(p.df, fm = 'pa', fa = 'both', n.iter = 1000, plot = F) #14 components
#
#data.frame('Factors' = (1:20),
#           'Eigenvalues' = parallel.analysis$pc.values[1:20],
#           'Sim_avg' = parallel.analysis$pc.sim[1:20],
#           'Sim_95' = (parallel.analysis$values %>%
#                         data.frame %>%
#                         select(CSim1:CSim20) %>%
#                         sapply(function(x){
#                           quantile(x, .95)
#                         })),
#           'one' = rep(1,20))  %>% write.excel(row.names = F, col.names = F)

#### Save and install list of R packages ######
# don't need to use this if packages are stored in a user library (as they are at VA) - just copy and paste folders
## Save list of R packages
#tmp = installed.packages()
#installedpackages = as.vector(tmp[is.na(tmp[,"Priority"]), 1])
#save(installedpackages, file="~/installed_packages.rda") # installs under this pc, documents folder

## Install R packages from list
#load("~/installed_packages.rda")
#
#for (count in 1:length(installedpackages)) {
#  install.packages(installedpackages[count])
#}


### Example code for programatically making new variables with standard name #####
#Line of code selects variables and makes new recoded variables with the suffix “binary”
#Tractsbin <- tracts %>%  mutate_at(vars(num_range("NSI_",1:4), num_range("NSI_",6:10),num_range("NSI_",12:22)), .funs = funs(binary = car::recode(., "0:1 = 0; 2:4 = 1"))) %>%
#  mutate_at(vars(PCL1,PCL11,PCL17), .funs = funs(binary = car::recode(., "1:2 = 0; 3:5 = 1")))

#### Calculate Mode from vector #####
## This if vector is multimodal, this function will just return the first value
#Mode <- function(x, na.rm = FALSE) {
#  if(na.rm){
#    x = x[!is.na(x)]
#  }
#
#  ux <- unique(x)
#  return(ux[which.max(tabulate(match(x, ux)))])
#}

# This is a modification of the above that will return the mean of the modes in the event that multiple values appear an equal number of times
Mode <- function(x, na.rm = FALSE) {
  if(na.rm){
    x = x[!is.na(x)]
  }
  ux <- unique(x)
  tabs <- tabulate(match(x, ux))
  max <- max(tabs)
  getmode <- ux[which(tabs == max)]
  ifelse(length(getmode) > 1,
         return(mean(getmode)),
         return(getmode))
}

#### Force a specific number of decimals #####
specify_decimal <- function(x, k) trimws(format(round(x, k), nsmall=k))

#### Number of parameters estimated ####

num.params <- function(observed.vars, num.factors){
  p <- observed.vars
  m <- num.factors
  a <- p*m + m*(m+1)/2 + p - m^2
  b <- p*(p+1)/2
  print('Number of parameters must be less than or equal to Vcov')
  return(c('Num.Params' = a, 'Num.Vcov' = b))
}

#### Syntax for taking vector output (e.g., result of call to names()) and convert to list separated by commas
# Example:
# names() %>% dput()

#### Function for identifying the edge strength of edges identified in NCT ouptut #######
# Needs to be generalized to take names
#edge.value <- function(nct){
#  sig <- filter(nct$einv.pvals, nct$einv.pvals$'p-value' < .05)
#  nw1 <- nct$nw1
#  nw2 <- nct$nw2
#  colnames(nw1) <- c("bor_ang","bor_aff","bor_emp","bor_idd","bor_par","bor_aband","bor_sib","bor_imp","bor_rel" )
#  rownames(nw1) <- c("bor_ang","bor_aff","bor_emp","bor_idd","bor_par","bor_aband","bor_sib","bor_imp","bor_rel" )
#  colnames(nw2) <- c("bor_ang","bor_aff","bor_emp","bor_idd","bor_par","bor_aband","bor_sib","bor_imp","bor_rel" )
#  rownames(nw2) <- c("bor_ang","bor_aff","bor_emp","bor_idd","bor_par","bor_aband","bor_sib","bor_imp","bor_rel")
#  result <- sapply(1:nrow(sig), function(x){
#    data.frame('NW1' = nw1[sig[x,"Var1"],sig[x,"Var2"]],
#               'NW2' = nw2[sig[x,"Var1"],sig[x,"Var2"]])
#  }) %>% t()
#  cbind(sig, result) %>% return()
#}

#Calculate edge differences ####
# Allows you to say this edge is greater than or less than X% of edges
calc.edge.dif <- function(x, alpha = .05, statistics = 'edge',
                          order = 'sample',
                          decreasing = T,
                          node.labels = NULL,
                          onlyNonZero = T,
                          bonferroni = F){
  cent <- x$bootTable %>% filter(type %in% statistics) %>% dplyr::select(name,id,value,type)

  if (onlyNonZero == T){
    include <-     unique(x$sampleTable$id[x$sampleTable$type %in% statistics & x$sampleTable$value != 0])
    cent <- cent %>% filter(id %in% include)
  } else {
    include <-     unique(x$sampleTable$id[x$sampleTable$type %in% statistics])
    cent <- cent %>% filter(id %in% include)
  }

  if (bonferroni == T){
    nInclude <- length(include)
    alpha <- alpha / (nInclude*(nInclude-1)/2)
    if (verbose) message(paste0("Significance level (alpha) set to: ",format(signif(alpha,2),scientific = FALSE)))
  }

  #if (verbose){
  #  exp <- expAlpha(alpha,length(x$boots))
  #  if (verbose) message(paste0("Expected significance level given number of bootstrap samples is approximately: ",format(signif(exp,2),scientific = FALSE)))
  #}

  fullTable <- expand.grid(name = unique(cent$name),id1=unique(cent$id),id2=unique(cent$id),type = unique(cent$type),
                           stringsAsFactors = FALSE)

  Quantiles <- fullTable %>%
    left_join(dplyr::select(cent,name,id1=id,value1=value,type),by=c("name","id1","type")) %>%
    left_join(dplyr::select(cent,name,id2=id,value2=value,type),by=c("name","id2","type"))  %>%
    group_by(id1,id2,type) %>%
    summarize(lower = quantile(value2-value1,alpha/2),upper = quantile(value2-value1,1-alpha/2)) %>%
    mutate(contain0 = 0 >= lower & 0 <= upper)

  #bootmean:
  bootMeans <- x$bootTable %>% filter(type %in% statistics) %>% rename(id1=id) %>%
    group_by(id1,type) %>% summarize(mean = mean(value,na.rm=TRUE))

  sample <- x$sampleTable %>% filter(type %in% statistics) %>% dplyr::select(id1=id,value,type) %>%
    left_join(bootMeans,by=c("id1","type"))

  # Now for every node: minimal node equal to....
  DF <-  sample %>% group_by(type) %>%
    mutate(rank = order(order(value,mean))) %>% arrange(rank)

  DF2 <- DF %>% filter(type == statistics[[1]])

  if (onlyNonZero == T){
    include <-     x$sampleTable$id[x$sampleTable$type %in% statistics & x$sampleTable$value != 0]
    DF2 <- DF2 %>% filter(id1 %in% include)
    DF <- DF %>% filter(id1 %in% include)
  }

  if (order == "sample"){
    levels <- DF2$id1[order(DF2$value,decreasing = !decreasing)]
  } else if (order == "mean"){
    levels <- DF2$id1[order(DF2$rank, decreasing = !decreasing)]
  } else  if (order == "id"){
    levels <- gtools::mixedsort(unique(include)) #changed this from sample$id1 - wasn't previously limiting to only non-zero edges
  }

  Quantiles$id1 <- factor(Quantiles$id1,levels=levels)
  Quantiles$id2 <- factor(Quantiles$id2,levels=levels)
  Quantiles$fill <- ifelse(Quantiles$id1 == Quantiles$id2, "same",
                           ifelse(Quantiles$contain0,"nonsig",
                                  ifelse(Quantiles$lower < 0 & Quantiles$upper < 0, 'greater',
                                         ifelse(Quantiles$lower > 0 & Quantiles$upper > 0, 'less', NA))))
  percent.dif <- group_by(Quantiles, id1) %>%
    summarize(greater = (sum(fill == 'greater'))/(length(levels(Quantiles$id1))-1),
              less = (sum(fill == 'less'))/(length(levels(Quantiles$id1))-1))
  return(percent.dif)
}

# Improve accuracy of data import ----
# Imporve read_csv by increasing the number of rows used for guessing the type of variable
# vasterling.df <- read_csv(here("data/vasterling_dat.csv"), guess_max = 3078)

# Foundation for counting the number of edges each node has with a magnitude greater than some cutoff ####
#pcl.btw.net.no.t$graph %>% apply(.,c(1,2),function(x){ifelse(abs(x) > .05, return(1), return(0))}) %>% apply(.,2,sum)

#### Attempts at generating figure from Bass-ackward results #####
# diagram function (https://cran.r-project.org/web/packages/diagram/vignettes/diagram.pdf) should also be considered
# I believe I made a partially functional function using that for the Narcissism project

# https://github.com/rich-iannone/DiagrammeR - use this as online vignette doesn't appear updated.
# example DiagrammeR syntax follows:

#ndf <-
#create_node_df(
#  n = 4,
#  label = c("a", "b", "c", "d"),
#  type  = "lower",
#  style = "filled",
#  color = "aqua",
#  shape = c("circle", "circle",
#            "rectangle", "rectangle"),
#  data = c(3.5, 2.6, 9.4, 2.7)
#)
#edf <-
#  create_edge_df(
#    label = c('1','.9','.85'),
#    from = c(1, 2, 3),
#    to   = c(4, 3, 1),
#    rel  = "leading_to"
#  )
#
#get_edge_df(graph)
#
#graph <-
#  create_graph(
#    nodes_df = ndf,
#    edges_df = edf
#  ) %>%
#  set_node_attrs(
#    node_attr = "fontname",
#    values = "Helvetica"
#  ) %>%
#  set_edge_attrs(
#    edge_attr = "color",
#    values = "blue"
#  ) %>%
#  set_edge_attrs(
#    edge_attr = "arrowsize",
#    values = 2
#  )
#
#graph %>% render_graph()

# Consider using bassAckward function syntax as partial template
# View(bassAckward.diagram)

#### FUNCTION IS NOT FINISHED AND WON'T WORK
# ___c. Attempt 3 #####
# Bass-ackward list is clps.bass
#clps.bass[[3]]$Vaccounted
#extract.vaccounted(clps.bass)
#
#
#scores.cor <- extract.scores(clps.bass) %>% cor(use = 'pairwise')
#scores.cor
#
#from <- c(rep(1,2), rep(2,3), rep(3,3))
#to <- c(rep(c(2:3),1), rep(c(4:6),2))
#edges <- tibble(from, to)
#label <- c()
#for(i in 1:8){
#  label <- c(label,
#             scores.cor[edges$to[i], edges$from[i]])
#}
#
#grViz('
#  digraph bass{
## node defintions
#
#node[fontname = Times, shape = rectangle];
##n0 [label = "Total Variance"; fixedsize = T; width = "@@12"] #
#n1 [label = "@@1"] #; fixedsize = T; width = "@@13"
#n2 [label = "@@2"] #; fixedsize = T; width = "@@14"
#n3 [label = "@@3"] #; fixedsize = T; width = "@@15"
#n4 [label = "@@4"] #; fixedsize = T; width = "@@16"
#n5 [label = "@@5"] #; fixedsize = T; width = "@@17"
#n6 [label = "@@6"] #; fixedsize = T; width = "@@18"
#
#subgraph cluster1
#   {
#       style = invis;
#       n1;
#   }
#
#subgraph cluster2
#   {
#       style = invis;
#       n2; n3;
#   }
#
#subgraph cluster23
#   {
#       style = invis;
#       n4; n5; n6;
#   }
#
#n1 -> n2 [label = "@@7"]
#n1 -> n3 [label = "@@8"]
#n2 -> n4 [label = "@@9"]
#n2 -> n6 [label = "@@10"]
#n3 -> n5 [label = "@@11"]
##n3 -> n6 [label = "@@12"]
#  }
#
#[1]: names(extract.loadings(clps.bass)[1])
#[2]: names(extract.loadings(clps.bass)[2])
#[3]: names(extract.loadings(clps.bass)[3])
#[4]: names(extract.loadings(clps.bass)[4])
#[5]: names(extract.loadings(clps.bass)[5])
#[6]: names(extract.loadings(clps.bass)[6])
#[7]: round(label[1],2)
#[8]: round(label[2],2)
#[9]: round(label[3],2)
#[10]: round(label[5],2)
#[11]: round(label[7],2)
#')
#
## ___d. Attempt 4 ######
#library(Gmisc)
#library(glue)
#library(htmlTable)
#library(grid)
#library(magrittr)
#
#clps.bass[[3]]$Vaccounted
#vacc <- extract.vaccounted(clps.bass)[[2]]
#
#
#scores.cor <- extract.scores(clps.bass) %>% cor(use = 'pairwise')
#scores.cor
#
#from <- c(rep(1,2), rep(2,3), rep(3,3))
#to <- c(rep(c(2:3),1), rep(c(4:6),2))
#edges <- tibble(from, to)
#label <- c()
#for(i in 1:8){
#  label <- c(label,
#             scores.cor[edges$to[i], edges$from[i]])
#}
#
#grid.newpage()
#
#vert <- spreadVertical(
#  'one' = boxGrob('one'),
#  'two' = boxGrob('two'),
#  'three' = boxGrob('three'),
#  .from = .025,
#  .to = .90
#)
#
#bot <- spreadHorizontal(
#  '1' = boxGrob("F3.1", width = vacc$F3.1[2], height = .05, y = coords(vert$one)$y),
#  '3' = boxGrob("F3.3", width = vacc$F3.3[2], height = .05, y = coords(vert$one)$y),
#  '2' = boxGrob("F3.2", width = vacc$F3.2[2], height = .05, y = coords(vert$one)$y),
#  .from = .025,
#  .to = .975
#)
#
#boxGrob("Total Variance", width = 1, height = .05, y = .975)
#bot
#
#f2.1 <- boxGrob("F2.1", width = vacc$F2.1[2], height = .05, y = coords(vert$two)$y, x = (coords(bot$`1`)$x+coords(bot$`3`)$x)/2)
#f2.2 <- boxGrob("F2.2", width = vacc$F2.2[2], height = .05, y = coords(vert$two)$y, x = coords(bot$`2`)$x)
#f2.1
#f2.2
#
#boxGrob("F1.1", width = vacc$F1.1[2], height = .05, y = coords(vert$three)$y, x = (coords(f2.1)$x + coords(f2.2)$x)/2)
#
# can't figure out how to label a connection

##### Figuring out ggplot defaults #####
#temp <- data.frame(x = 1:2, y = rep(20:1, each = 2), grp = factor(rep(1:20, each = 2)))
#temp
#
## plot
#p <- ggplot(data = temp, aes(x = x, y = y, linetype = grp, colour = grp)) +
#  geom_line() +
#  geom_text(aes(x = 0.95, label = grp)) +
#  theme_classic() +
#  theme(axis.title = element_blank(),
#        axis.text = element_blank(),
#        axis.ticks = element_blank(),
#        axis.line = element_blank(),
#        legend.position = "none")
#p
#
## Just two lines, see below for pulling data
#temp <- data.frame(x = 1:2, y = rep(20:19, each = 2), grp = factor(rep(1:2, each = 2)))
#temp
#
## plot
#p <- ggplot(data = temp, aes(x = x, y = y, linetype = grp, colour = grp)) +
#  geom_line() +
#  geom_text(aes(x = 0.95, label = grp)) +
#  theme_classic() +
#  theme(axis.title = element_blank(),
#        axis.text = element_blank(),
#        axis.ticks = element_blank(),
#        axis.line = element_blank(),
#        legend.position = "none")
#p
#
#g <- ggplot_build(p)
#g$data[[1]]

# Should make function to easily plot relationship between two variables ----
# ggplot(out, aes(x=var.1, y=var.2)) + geom_point() + coord_fixed() + geom_smooth(method='lm')

# SEM Plot figure ----
# This syntax resulted in a pretty decent CFA plot; might be able to build something from that
# the five_factor_cfa object is a lavaan model object

# semPaths(five_factor_cfa, 'model','estimates', 'lisrel',
#          'tree2',
#          intercepts = F, nCharNodes = 0,
#          sizeMan = 3, sizeMan2 = 3,
#          sizeLat = 8, sizeLat2 = 5,
#          fixedStyle = c('black', 2),
#          freeStyle = 'black',
#          #layoutSplit = T,
#          #levels = c(1,1.1), #can't figure out how to make levels closer using this function
#          nDigits = 2,
#          label.cex = 2,
#          #filetype = 'pdf',
#          #width = 10, height = 6,
#          mar = c(3,3,3,3),
#          normalize = T)

# Calculate from summary ----
# The following is example syntax I used to calculate McFadden's r-square from the summary results of logistic regression
# The with() function seems useful here, I need to understand that better
# glm(SCID1_PTSD_Currentdx_T2 ~  scale(PCL5tot_T2) + scale(PHQMDDsev_T2), data = fa.df, family = 'binomial') %>% summary %>% with(1-deviance/null.deviance)

# Detect warnings ----
# Detect Warnings in lavaan models
#' Catch errors and warnings and store them for subsequent evaluation
#'
#' Factory modified from a version written by Martin Morgan on Stack Overflow (see below).
#' Factory generates a function which is appropriately wrapped by error handlers.
#' If there are no errors and no warnings, the result is provided.
#' If there are warnings but no errors, the result is provided with a warn attribute set.
#' If there are errors, the result is a list with the elements of warn and err.
#' This is a nice way to recover from a problems that may have occurred during loop evaluation or during cluster usage.
#' Check the references for additional related functions.
#' I have not included the other factory functions included in the original Stack Overflow
#' answer because they did not play well with the return item as an S4 object.
#' @param fun The function to be turned into a factory
#' @return The result of the function given to turn into a factory.  If this function was in error "An error as occurred" as a character element.  factory-error and factory-warning attributes may also be set as appropriate.
#' @references
#' \url{http://stackoverflow.com/questions/4948361/how-do-i-save-warnings-and-errors-as-output-from-a-function}
#' @author Martin Morgan; Modified by Russell S. Pierce; Modified by Michael Crowe
#' @examples
#' f.log <- factory(log)
#' f.log("a")
#' f.as.numeric <- factory(as.numeric)
#' f.as.numeric(c("a","b",1))
#' @export
factory <- function (fun) {
  errorOccurred <- FALSE
  library(data.table)
  function(...) {
    warn <- err <- NULL
    res <- withCallingHandlers(tryCatch(fun(...), error = function(e) {
      err <<- conditionMessage(e)
      errorOccurred <<- TRUE
      NULL
    }), warning = function(w) {
      warn <<- append(warn, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
    # The following three lines were part of the original syntax that I downloaded, but including that
    # appears to replace the results of the function with the character string, which makes it difficult
    # to check to see what the error actually was
    # if (errorOccurred) {
    #   res <- "An error occurred in the factory function"
    # }

    if (is.character(warn)) {
      data.table::setattr(res,"factory-warning",warn)
    } else {
      data.table::setattr(res,"factory-warning",NULL)
    }

    if (is.character(err)) {
      data.table::setattr(res,"factory-error",err)
    } else {
      data.table::setattr(res, "factory-error", NULL)
    }
    return(res)
  }
}

.has <- function(x, what) {
  !is.null(attr(x,what))
}

hasWarning <- function(x) {.has(x, "factory-warning")}
hasError <- function(x) {.has(x, "factory-error")}
isClean <- function(x) {!(hasError(x) | hasWarning(x))}

# # attr([object],"factory-warning") # to see warning
# # Need to use this function as a wrapper, so if running CFA models, you need to use:
# f.cfa <- factory(cfa)

# # example:
# fit1s <- map(balanced_dfs[1:100], function(x){
#   f.cfa(mod, data = x, missing = 'pairwise', estimator = 'WLSMV', group = "white_black")
# })
# beep()
#
# fit1s_no_warnings <- map(fit1s, ~attr(.x, "factory-warning")) %>%
#   sapply(is.null)

# Compare observed cor matrix to model implied cor matrix ----
# from Thomas: You might appreciate this: I made a function to make a correlation table based on an observed vs model-implied correlation matrix: the upper diagonal shows the observed correlations, and the lower diagonal shows the model-implied correlations, and the color of each cell is proportional to the difference for the respective correlation between observed and model-implied. the only arguments you need are 1) the fitted lavaan object from cfa() or whatever and 2) the observed correlation matrix created with cor()
create_comparison_table <- function(fitted_lavaan_object, observed_cor_matrix) {

  implied_corrmat <- lavaan::lavInspect(fitted_lavaan_object, 'cov.ov', add.class=FALSE) |> cov2cor()

  sortcol_impcormat <- implied_corrmat[, order(colnames(implied_corrmat))]
  sortrow_impcormat <- sortcol_impcormat[order(rownames(sortcol_impcormat)), ]

  model_implied_cor <- sortrow_impcormat
  observed_cor <- observed_cor_matrix

  combined_cor <- observed_cor
  combined_cor[lower.tri(combined_cor)] <- model_implied_cor[lower.tri(model_implied_cor)]

  difference_matrix <- abs(model_implied_cor - observed_cor)

  double_df_with_differences <- combined_cor |> tibble::as_tibble(rownames = "var") |>
    dplyr::bind_cols(
      difference_matrix |> tibble::as_tibble() |> dplyr::rename_with(~ paste0("diff_", .x))
    )

  gt_table <- gt::gt(double_df_with_differences) %>%
    gt::fmt_number(columns = tidyselect::everything(), decimals = 2) %>%
    gt::data_color(
      columns = tidyselect::starts_with('diff_'),
      fn = scales::col_numeric(
        palette = "Blues",
        domain = c(min(difference_matrix), max(difference_matrix))
      ),
      target_columns = !tidyselect::starts_with('diff_') & !tidyselect::starts_with('var')
    ) |>
    gt::cols_hide(columns = tidyselect::starts_with('diff_')) |>
    gt::tab_header(
      title = gt::md("Observed (upper) vs Model-implied (lower) Correlations"),
      subtitle = "Color corresponds to the magnitude of the difference between observed and model-implied matrices"
    ) |>
    gt::cols_align('center')

  return(gt_table)
}

# Install package from github ----
# reminder of how to install package directly from github
# remotes::install_github("yrosseel/lavaan")

# Superscripts for correlation differences ----
# Note that tests value is dataframe of correlation difference test p-values
#          1v2          1v3          1v4          2v3          2v4          3v4
# 2.871070e-01 4.112781e-02 5.599676e-11 3.602023e-01 7.344505e-09 8.336107e-06

# # Syntax I have used previously to generate the "tests" dataframe for this function follows:
# cors <- select(df, c('part','obs', 'decl', 'DRRICombattot_T1', criteria)) %>%
#   as.matrix %>%
#   rcorr(., type = 'spearman')
#
# ## Test correlation differences
# cor.dif <- sapply(rownames(cors$r)[5:length(rownames(cors$r))], function(third.var){
#   lapply(1:3, function(x){
#     lapply(1:(4-x), function(y){
#       out <- cocor.dep.groups.overlap(cors$r[third.var,x], cors$r[third.var,(x+y)], cors$r[x,(x+y)], min(cors$n[third.var,(x+y)], cors$n[third.var,(x+y)], cors$n[x,(x+y)]), test = 'meng1992')@meng1992$p.value
#       names(out) = paste0(x,'v',(x+y))
#       return(out)
#     }) %>% do.call(c, .)
#   }) %>% do.call(c, .)
# }) %>% t

# Can apply to dataframe using apply function as in:
# apply(cor.dif, 1, labs, nvar = 4, pvalue = .01) %>% t

labs <- function(tests, nvar = 4, pvalue = .01){
  "%next%" <- function(x, y) x[!x %in% y][1] #--  x without y
  labels <- vector('character',nvar)
  names(labels) <- paste0('v',1:nvar)
  mat <- matrix(nrow = nvar, ncol = nvar)
  mat[lower.tri(mat, diag = F)] <- tests
  mat[upper.tri(mat, diag = F)] <- t(mat)[upper.tri(t(mat),diag = F)]
  used.letters <- c()
  if(any(mat <= pvalue, na.rm = T)){
    for(i in 1:nvar){
      if(any(mat[,i] <= pvalue, na.rm = T)){
        labels[[i]] <- paste0(labels[[i]],unique(labels[which(mat[,i] > pvalue)]), collapse = '') # assign labels for vars that are the same
        if(str_length(labels[which(mat[,i] < pvalue)] %>% str_split('') %>% unlist() %>% paste0(collapse = '|')) > 0){
          labels[[i]] <- str_replace(labels[[i]],
                                     labels[which(mat[,i] < pvalue)] %>% str_split('') %>% unlist() %>% paste0(collapse = '|'), # any of the labels it can't use
                                     '')
        }
        if(str_length(labels[[i]]) == 0){
          used.letters <- c(used.letters,letters %next% used.letters)
          labels[[i]] <- tail(used.letters, 1)
        }
      }
    }
    for(i in nvar:1){
      if(any(mat[,i] > pvalue, na.rm = T)){
        for(j in which(mat[,i] > pvalue)){
          if(!any((labels[[i]] %>% str_split('') %>% unlist()) %in% (labels[[j]] %>% str_split('') %>% unlist())) & labels[[j]] != ""){
            labels[[j]] <- paste0(labels[[j]], labels[[i]], collapse = '')
          }
        }
      }
    }
    return(labels)
  }else{
    return(labels)
  }
}

# # What follows is he results of a gpt-4 prompt in which I asked it to provide an improved version of the above function
# # I tried running it using data from the atrocities project and it resulted in an error - didn't spend any time trying to correct the function
# assignLabelsBasedOnPValue <- function(tests, nvar = 4, pvalue = .01){
#   if(!is.numeric(tests)) stop("'tests' must be a numeric vector or matrix")
#   if(length(tests) != nvar*(nvar-1)/2) stop("'tests' must be of length nvar*(nvar-1)/2")
#
#   nextElement <- function(x, y) x[!x %in% y][1] # x without y
#
#   labels <- paste0('v',1:nvar)
#   names(labels) <- labels
#
#   mat <- matrix(nrow = nvar, ncol = nvar)
#   mat[lower.tri(mat, diag = FALSE)] <- tests
#   mat[upper.tri(mat, diag = FALSE)] <- t(mat)[upper.tri(t(mat),diag = FALSE)]
#
#   used.letters <- c()
#
#   if(any(mat <= pvalue, na.rm = TRUE)){
#     for(i in 1:nvar){
#       if(any(mat[,i] <= pvalue, na.rm = TRUE)){
#         labels[i] <- paste(labels[i], unique(labels[mat[,i] > pvalue]), collapse = '')
#         forbidden_labels <- labels[mat[,i] < pvalue]
#         if(nchar(forbidden_labels) > 0){
#           labels[i] <- gsub(paste(strsplit(forbidden_labels, "")[[1]], collapse = '|'), '', labels[i])
#         }
#         if(nchar(labels[i]) == 0){
#           used.letters <- c(used.letters, nextElement(letters, used.letters))
#           labels[i] <- tail(used.letters, 1)
#         }
#       }
#     }
#     for(i in nvar:1){
#       if(any(mat[,i] > pvalue, na.rm = TRUE)){
#         for(j in which(mat[,i] > pvalue)){
#           if(!any(strsplit(labels[i], "")[[1]] %in% strsplit(labels[j], "")[[1]]) & labels[j] != ""){
#             labels[j] <- paste(labels[j], labels[i], collapse = '')
#           }
#         }
#       }
#     }
#   }
#   return(labels)
# }

# Calculate H index of factor reliability ----
# Rodriguez, A., Reise, S. P., & Haviland, M. G. (2016). Evaluating bifactor models: Calculating and interpreting statistical indices. Psychological Methods, 21, 137–150. https://doi.org/10.1037/met0000045
# l is vector of fully standardized factor loadings

h_index <- function(l){
  ratio <- (l^2)/(1-(l^2))
  sum_ratio <- sum(ratio)
  denom <- 1 + (1/sum_ratio)
  h <- 1/denom
  h
}

# function for copying and pasting file paths ----
# By default it saves the path to clipboard, setting value to false makes it return the path in string format
get_path <- function(to_clipboard = T){
  path_from_clipboard <- readClipboard()
  # Remove leading and trailing quotes
  trimmed_path <- gsub('^"|"$', '', path_from_clipboard)
  # Replace backslashes with forward slashes
  new_path <- gsub("\\\\", "/", trimmed_path, fixed = F)
  if(to_clipboard == T){
    writeClipboard(new_path, format = 13)
  }
  if(to_clipboard == F)
    return(new_path)
}

# Use
# 1) find file in windows explorer, shift+right_click, select copy as path
# get_path() # running this will save the correct file path in the clipboard

# Calculate RMSEAd ----
# based on Savalei et al. 2024

# Define the compare_nested_models_RMSEA function
# compare_nested_models_RMSEA <- function(models, model_names = names(models)) {

# Input can be a list of models or individual lavaan models
# e.g.,
# compare_nested_models_RMSEA(fit0, fit1)
# or
# model_list <- list(fit0, fit1, fit2)
# compare_nested_models_RMSEA(model_list)
compare_nested_models_RMSEA <- function(..., model_names = NULL, CI = .90) {
  require(lavaan)

  # Capture all input arguments
  models_input <- list(...)

  # Check if the first argument is a list of models
  if (length(models_input) == 1 && is.list(models_input[[1]]) &&
      all(sapply(models_input[[1]], function(x) inherits(x, "lavaan")))) {
    models <- models_input[[1]]
  } else {
    # Assume that individual models are passed directly
    models <- models_input
  }

  # Check if at least two models are provided
  if (length(models) < 2) {
    stop("Please provide at least two nested lavaan models for comparison.")
  }

  # Assign default model names if not provided
  if (is.null(model_names)) {
    model_names <- paste0("Model_", 1:length(models))
  } else {
    # Check if the number of names matches the number of models
    if (length(model_names) != length(models)) {
      stop("Length of model_names must match the number of models provided.")
    }
  }

  # Ensure that all elements in the models list are lavaan objects
  if (!all(sapply(models, function(x) inherits(x, "lavaan")))) {
    stop("All elements in the models list must be lavaan model objects.")
  }

  # Extract number of groups from each model
  n_groups <- sapply(models, function(x) lavaan::inspect(x, "ngroups"))

  # Check if all models have the same number of groups
  if(length(unique(n_groups)) != 1){
    stop("All models must have the same number of groups.")
  }

  G <- n_groups[1]  # Number of groups

  # Internal RMSEA.CI function
  RMSEA.CI <- function(T, df, N, G, CI){

    lower.lambda <- function(lambda) {
      (pchisq(T, df=df, ncp=lambda) - (1 - (1-CI)/2))
    }
    upper.lambda <- function(lambda) {
      (pchisq(T, df=df, ncp=lambda) - ((1-CI)/2))
    }

    # RMSEA CI Lower Bound
    lambda.l <- try(uniroot(f=lower.lambda, lower=0, upper=T)$root, silent=TRUE)
    if(inherits(lambda.l, "try-error")) {
      lambda.l <- NA
      RMSEA.CI.l <- NA  # Could consider setting lower CI to 0 if if uniroot fails
    } else {
      if(lambda.l < 0){
        RMSEA.CI.l <- 0
      } else {
        RMSEA.CI.l <- sqrt(lambda.l * G / ((N - 1) * df))
      }
    }

    # RMSEA CI Upper Bound
    N.RMSEA <- max(N, T * 4)
    lambda.u <- try(uniroot(f=upper.lambda, lower=0, upper=N.RMSEA)$root, silent=TRUE)
    if(inherits(lambda.u, "try-error")) {
      lambda.u <- NA
      RMSEA.CI.u <- NA
    } else {
      if(lambda.u < 0){
        RMSEA.CI.u <- 0
      } else {
        RMSEA.CI.u <- sqrt(lambda.u * G / ((N - 1) * df))
      }
    }

    RMSEA.CI <- c(RMSEA.CI.l, RMSEA.CI.u)
    names(RMSEA.CI) <- c("Lower_CI", "Upper_CI")
    return(RMSEA.CI)
  }

  # Internal pvals function
  pvals <- function(T, df, N, G){
    RMSEA0 <- c(0, 0.01, 0.05, 0.06, 0.08, 0.1)
    eps0 <- df * RMSEA0^2 / G
    nonc <- eps0 * (N - G)
    pvals <- pchisq(T, df=df, ncp=nonc)
    names(pvals) <- c("RMSEA>0", "RMSEA>.01", "RMSEA>.05", "RMSEA>.06", "RMSEA>.08", "RMSEA>.10")
    return(pvals)
  }

  # Extract fit indices, sample size, and number of parameters from each model
  fit_indices <- data.frame(
    Model = model_names,
    ChiSq = sapply(models, function(x) fitMeasures(x, "chisq")),
    df = sapply(models, function(x) fitMeasures(x, "df")),
    RMSEA = sapply(models, function(x) fitMeasures(x, "rmsea")),
    stringsAsFactors = FALSE
  )

  # Extract sample size (assuming all models have the same sample size)
  # changed "nobs" to "ntotal" as "nobs" didn't work when there were multiple groups
  # my read of the Savalei et al. 2024 manuscript is that all of the calculations using sample size use total sample rather than group sample
  N_values <- sapply(models, function(x) inspect(x, "ntotal"))
  if(length(unique(N_values)) != 1){
    stop("All models must have the same sample size (ntotal).")
  }
  N <- N_values[1]

  # Initialize vectors to store comparison results
  comparison <- paste(model_names[-length(model_names)], "vs", model_names[-1])
  Delta_ChiSq <- vector("numeric", length = length(models) - 1)
  Delta_df <- vector("numeric", length = length(models) - 1)
  Delta_RMSEA <- vector("numeric", length = length(models) - 1)
  RMSEA_d <- vector("numeric", length = length(models) - 1)
  RMSEA_d_Lower_CI <- vector("numeric", length = length(models) - 1)
  RMSEA_d_Upper_CI <- vector("numeric", length = length(models) - 1)
  pvals_list <- vector("list", length = length(models) - 1)

  # Loop through consecutive model pairs
  for(i in 2:length(models)) {
    # Chi-square and df differences
    Delta_ChiSq[i - 1] <- fit_indices$ChiSq[i] - fit_indices$ChiSq[i - 1]
    Delta_df[i - 1] <- fit_indices$df[i] - fit_indices$df[i - 1]

    # G is number of groups
    # Already extracted before the loop

    # Compute RMSEAd
    numerator <- Delta_ChiSq[i - 1] - Delta_df[i - 1]
    denominator <- Delta_df[i - 1] * (N - G) / G

    if(denominator <= 0) {
      warning(paste("Non-positive denominator for RMSEAd computation between",
                    model_names[i - 1], "and", model_names[i], ". Setting RMSEAd to NA."))
      rmsea_d_value <- NA
    } else {
      if(numerator <= 0){
        rmsea_d_value <- 0
      } else {
        rmsea_d_value <- sqrt(numerator / denominator)
      }
    }
    RMSEA_d[i - 1] <- rmsea_d_value

    # Compute RMSEA CI for RMSEAd
    if(!is.na(rmsea_d_value)) {
      ci <- RMSEA.CI(T = Delta_ChiSq[i - 1], df = Delta_df[i - 1], N = N, G = G, CI = CI)
      RMSEA_d_Lower_CI[i - 1] <- ci["Lower_CI"]
      RMSEA_d_Upper_CI[i - 1] <- ci["Upper_CI"]
    } else {
      RMSEA_d_Lower_CI[i - 1] <- NA
      RMSEA_d_Upper_CI[i - 1] <- NA
    }

    # Equivalence tests based on RMSEAd
    if(!is.na(rmsea_d_value)) {
      pvals_comp <- pvals(T = Delta_ChiSq[i - 1], df = Delta_df[i - 1], N = N, G = G)
      pvals_list[[i - 1]] <- round(pvals_comp, 3)
    } else {
      pvals_list[[i - 1]] <- rep(NA, 6)
      names(pvals_list[[i - 1]]) <- c("RMSEA>0", "RMSEA>.01", "RMSEA>.05",
                                      "RMSEA>.06", "RMSEA>.08", "RMSEA>.10")
    }

    # Compute Delta RMSEA
    Delta_RMSEA[i - 1] <- fit_indices$RMSEA[i] - fit_indices$RMSEA[i - 1]
  }

  # Compile p-values into separate columns
  pvals_df <- do.call(rbind, pvals_list)
  rownames(pvals_df) <- NULL

  # Create the final comparison data frame
  comparison_df <- data.frame(
    Comparison = comparison,
    Delta_ChiSq = Delta_ChiSq,
    Delta_df = Delta_df,
    Delta_RMSEA = Delta_RMSEA,
    RMSEA_d = RMSEA_d,
    RMSEA_d_Lower_CI = RMSEA_d_Lower_CI,
    RMSEA_d_Upper_CI = RMSEA_d_Upper_CI,
    #p_RMSEA_gt0 = pvals_df[, "RMSEA>0"],
    #p_RMSEA_gt0.01 = pvals_df[, "RMSEA>.01"],
    p_RMSEA_gt0.05 = pvals_df[, "RMSEA>.05"],
    p_RMSEA_gt0.06 = pvals_df[, "RMSEA>.06"],
    p_RMSEA_gt0.08 = pvals_df[, "RMSEA>.08"],
    p_RMSEA_gt0.10 = pvals_df[, "RMSEA>.10"],
    stringsAsFactors = FALSE,
    row.names = NULL
  )

  return(comparison_df)
}

# Plot spline with confidence intervals ----
# function from Thomas Crowe (personal communication)
fit_and_plot_spline_interaction <- function(DATASET,
                                            OUTCOME, NUMERICVAR, CATEGORICAL_INTERACTION_VAR = NULL,
                                            NKNOTS=5, print_simplenames=TRUE, overwrite_model_names=FALSE) {
  # Capture expressions and convert to strings
  NUMERICVAR_expr <- rlang::enquo(NUMERICVAR)
  OUTCOME_expr <- rlang::enquo(OUTCOME)

  # Convert to character strings for use in functions that require character input
  NUMERICVAR_str <- rlang::as_name(NUMERICVAR_expr)
  OUTCOME_str <- rlang::as_name(OUTCOME_expr)

  # Check if CATEGORICAL_INTERACTION_VAR is NULL
  has_interaction <- !rlang::quo_is_null(rlang::enquo(CATEGORICAL_INTERACTION_VAR))

  if (has_interaction) {
    CATEGORICAL_INTERACTION_VAR_expr <- rlang::enquo(CATEGORICAL_INTERACTION_VAR)
    CATEGORICAL_INTERACTION_VAR_str <- rlang::as_name(CATEGORICAL_INTERACTION_VAR_expr)
  }

  # Check inputs
  if (!is.data.frame(DATASET)) {
    stop("DATASET must be a data frame")
  }

  if (!NUMERICVAR_str %in% names(DATASET)) {
    stop(paste("Variable", NUMERICVAR_str, "not found in the dataset"))
  }

  if (has_interaction && !CATEGORICAL_INTERACTION_VAR_str %in% names(DATASET)) {
    stop(paste("Variable", CATEGORICAL_INTERACTION_VAR_str, "not found in the dataset"))
  }

  if (!OUTCOME_str %in% names(DATASET)) {
    stop(paste("Variable", OUTCOME_str, "not found in the dataset"))
  }

  # Extract the numeric variable as a vector
  numvar_vector <- DATASET[[NUMERICVAR_str]]

  if (!is.numeric(numvar_vector)) {
    stop(paste(NUMERICVAR_str, "must be a numeric variable"))
  }

  # Embedded getquantiles function
  getquantiles <- function(numericvar, nknots) {
    if (!is.numeric(numericvar)) {
      stop("Input must be a numeric vector")
    }
    if (!(nknots %in% 3:7)) {
      stop("nknots must be between 3 and 7")
    }
    quantiles <- switch(as.character(nknots),
                        "3" = c(0.10, 0.50, 0.90),
                        "4" = c(0.05, 0.35, 0.65, 0.95),
                        "5" = c(0.05, 0.275, 0.50, 0.725, 0.95),
                        "6" = c(0.05, 0.23, 0.41, 0.59, 0.77, 0.95),
                        "7" = c(.025, .1833, .3417, .5, .6583, .8167, .975))
    knot_quantiles <- quantile(numericvar, probs = quantiles, na.rm = TRUE)
    result <- list(
      inner = knot_quantiles[-c(1, length(knot_quantiles))],
      outer = c(min(knot_quantiles), max(knot_quantiles)),
      all = knot_quantiles |> as.double(),
      quants = quantiles
    )
    return(result)
  }

  # Get quantiles for spline with nknots
  quant_result <- getquantiles(numvar_vector, nknots = NKNOTS)

  # Check outcome type and fit appropriate model
  outcome_var <- DATASET[[OUTCOME_str]]

  # Create formula based on whether we have an interaction term
  if (has_interaction) {
    formula_str <- paste(OUTCOME_str, "~", CATEGORICAL_INTERACTION_VAR_str, "* splines::ns(",
                         NUMERICVAR_str, ", knots = c(",
                         paste(quant_result$inner, collapse = ", "),
                         "), Boundary.knots = c(",
                         paste(quant_result$outer, collapse = ", "), "))")
  } else {
    formula_str <- paste(OUTCOME_str, "~ splines::ns(",
                         NUMERICVAR_str, ", knots = c(",
                         paste(quant_result$inner, collapse = ", "),
                         "), Boundary.knots = c(",
                         paste(quant_result$outer, collapse = ", "), "))")
  }

  actual_formula <- as.formula(formula_str)

  if (is.numeric(outcome_var) || is.double(outcome_var)) {
    # Fit linear model for numeric outcome
    call_expr <- substitute(
      lm(formula = FORMULA, data = DATANAME),
      list(FORMULA = actual_formula, DATANAME = substitute(DATASET))
    )
    mymodel <- eval(call_expr)

  } else if ((is.factor(outcome_var) || is.character(outcome_var)) &&
             length(unique(na.omit(outcome_var))) == 2) {
    # Fit binomial GLM for binary outcome
    call_expr <- substitute(
      glm(formula = FORMULA, family = binomial(link = "logit"), data = DATANAME),
      list(FORMULA = actual_formula, DATANAME = substitute(DATASET))
    )
    mymodel <- eval(call_expr)

  } else if ((is.factor(outcome_var) || is.character(outcome_var)) &&
             length(unique(na.omit(outcome_var))) > 2) {
    stop("OUTCOME has more than two levels. Only binary or numeric outcomes are supported.")
  } else {
    stop("OUTCOME must be numeric or a binary factor/character variable")
  }

  # Generate plots using marginaleffects
  if (has_interaction) {
    # With interaction: plot both variables
    predsplot <- marginaleffects::plot_predictions(
      mymodel,
      condition = c(NUMERICVAR_str, CATEGORICAL_INTERACTION_VAR_str)
    )

    comparisonplot <- marginaleffects::plot_comparisons(
      mymodel,
      variables = CATEGORICAL_INTERACTION_VAR_str,
      condition = NUMERICVAR_str) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dotted", color = "#FF9999", size = 1) +
      ggplot2::annotate("text", x = -Inf, y = 0, label = "No difference",
                        hjust = -0.1, vjust = -0.5, color = "#FF9999", size = 3.5)

    # Combine plots using patchwork
    combined_plot <- patchwork:::`|.ggplot`(predsplot, comparisonplot) +
      patchwork::plot_layout(ncol = 2)
  } else {
    # Without interaction: just plot the numeric variable
    predsplot <- marginaleffects::plot_predictions(
      mymodel,
      condition = NUMERICVAR_str
    )

    # No comparison plot or combined plot when there's no interaction
    comparisonplot <- NULL
    combined_plot <- predsplot
  }

  # Create simplified names for coefficients
  coef_names <- names(coef(mymodel))
  spline_pattern <- paste0("splines::ns\\(", NUMERICVAR_str,
                           ", knots = c\\([^)]+\\), Boundary\\.knots = c\\([^)]+\\)\\)")
  simplified_names <- gsub(spline_pattern, paste0("ns(", NUMERICVAR_str, ")"), coef_names)

  if (overwrite_model_names) {
    # Actually modify the model's coefficient names
    names(mymodel$coefficients) <- simplified_names

    # For lm/glm objects, also update terms and other components
    if (inherits(mymodel, "lm") || inherits(mymodel, "glm")) {
      # Update model matrix column names if they exist
      if (!is.null(mymodel$model) && is.data.frame(mymodel$model)) {
        spline_cols <- grep(spline_pattern, colnames(mymodel$model))
        if (length(spline_cols) > 0) {
          simplified_model_cols <- gsub(spline_pattern, paste0("ns(", NUMERICVAR_str, ")"),
                                        colnames(mymodel$model)[spline_cols])
          colnames(mymodel$model)[spline_cols] <- simplified_model_cols
        }
      }

      # Update terms
      term_labels <- attr(mymodel$terms, "term.labels")
      simplified_terms <- gsub(spline_pattern, paste0("ns(", NUMERICVAR_str, ")"), term_labels)
      attr(mymodel$terms, "term.labels") <- simplified_terms

      # Update contrasts if present
      if (!is.null(mymodel$contrasts)) {
        names(mymodel$contrasts) <- gsub(spline_pattern,
                                         paste0("ns(", NUMERICVAR_str, ")"),
                                         names(mymodel$contrasts))
      }

      # Update xlevels if present
      if (!is.null(mymodel$xlevels)) {
        names(mymodel$xlevels) <- gsub(spline_pattern,
                                       paste0("ns(", NUMERICVAR_str, ")"),
                                       names(mymodel$xlevels))
      }
    }

    # Print the modified model summary
    print(summary(mymodel))
  } else if (print_simplenames) {
    # Create a temporary clean model for printing only
    clean_model <- mymodel
    names(clean_model$coefficients) <- simplified_names

    # For glm objects
    if (inherits(clean_model, "lm") || inherits(clean_model, "glm")) {
      term_labels <- attr(clean_model$terms, "term.labels")
      simplified_terms <- gsub(spline_pattern, paste0("ns(", NUMERICVAR_str, ")"), term_labels)
      attr(clean_model$terms, "term.labels") <- simplified_terms
    }

    # Print the clean model summary
    print(summary(clean_model))
  } else {
    # Print the original model summary
    print(summary(mymodel))
  }

  # Print the plot(s)
  if (has_interaction) {
    print(combined_plot)
  } else {
    print(predsplot)
  }

  # Create result list
  result <- list(
    mymodel = mymodel,
    predsplot = predsplot,
    has_interaction = has_interaction
  )

  # Add interaction-specific elements only if they exist
  if (has_interaction) {
    result$comparisonplot <- comparisonplot
    result$combined_plot <- combined_plot
  }

  # Set the class to allow proper printing
  class(result) <- c(class(result), "spline_interaction_model")

  # Return the result object invisibly
  return(invisible(result))
}

# Example usage
df591 |> fit_and_plot_spline_interaction(droppedout, NUMERICVAR = age, CATEGORICAL_INTERACTION_VAR = treatment_group)
# result should be a model summary (with splines) and two corresponding plots patched together into one,
# or just one plot if you don't specify a categorical interaction variable
# if you want to be more conservative, you can set NKNOTS=4 or NKNOTS=3 to get fewer bends

# Attempt at parallel processing wrapper ----
# Function doesn't currently work correctly
#' Run a Function in Parallel with Progress Bar and Optional Beep
#'
#' This function provides a convenient wrapper for parallelized application of any function
#' using \code{future_lapply} (for single input) or \code{future_mapply} (for multiple inputs),
#' while displaying a progress bar. It supports flexible core allocation and optionally plays a sound when complete.
#'
#' @param FUN The function to apply. Should take as its first argument the elements of the input (as in \code{lapply}) or the arguments to be mapped over (as in \code{mapply}).
#' @param ... Arguments to be passed as inputs to \code{FUN}, typically one or more lists or vectors to iterate over.
#' @param steps Number of steps for the progress bar. If \code{NULL} (default), will be automatically set to the length of the first input.
#' @param future_mode \code{"lapply"} (default) for a single mapped input, or \code{"mapply"} for mapping over multiple inputs simultaneously.
#' @param show_beep Logical, whether to play a sound upon completion. Default is \code{TRUE}. Requires the \code{beepr} package.
#' @param workers Number of parallel workers. Default is \code{parallelly::availableCores() - 1}.
#' @param progress_format A string template for the progress bar, passed to \code{progressr::handler_progress}.
#' @param .seed Seed for parallel random number generation. Default is \code{NULL}.
#'
#' @return A list of the results from \code{FUN} applied in parallel to the input(s).
#'
#' @examples
#' # lapply example: squaring numbers in parallel
#' x <- 1:10
#' squares <- with_progress_parallel(FUN = function(x) x^2, x)
#'
#' # mapply example: adding two vectors in parallel
#' a <- 1:5
#' b <- 6:10
#' sums <- with_progress_parallel(FUN = function(x, y) x + y, a, b, future_mode = "mapply")
#'
#' @import future.apply
#' @import future
#' @importFrom progressr handlers handler_progress with_progress progressor
#' @importFrom parallelly availableCores
#' @importFrom beepr beep
#' @export
with_progress_parallel <- function(FUN, ...,
                                   steps = NULL,
                                   future_mode = c("lapply", "mapply"),
                                   show_beep = TRUE,
                                   workers = parallelly::availableCores() - 1,
                                   progress_format = "[:bar] :percent :eta :message",
                                   .seed = NULL) {
  progressr::handlers(progressr::handler_progress(format = progress_format))
  future::plan(multisession, workers = workers)

  old_opt <- getOption("progressr.enable")
  options(progressr.enable = TRUE)
  on.exit({
    options(progressr.enable = old_opt)
    plan(sequential)
  }, add = TRUE)

  args <- list(...)
  future_mode <- match.arg(future_mode)
  if (is.null(steps)) steps <- length(args[[1]])
  result <- NULL

  progressr::with_progress({
    p <- progressr::progressor(steps = steps)
    if (future_mode == "lapply") {
      result <- future.apply::future_lapply(
        args[[1]],
        function(x, ...) {
          p()
          FUN(x, ...)
        },
        ...,
        future.seed = .seed
      )
    } else if (future_mode == "mapply") {
      result <- future.apply::future_mapply(
        function(...) {
          p()
          FUN(...)
        },
        ...,
        future.seed = .seed,
        SIMPLIFY = FALSE
      )
    }
  })
  future::plan(sequential)
  options(old_o)
  if (show_beep && requireNamespace("beepr", quietly = TRUE)) {
    beepr::beep()
  }
  return(result)
}


# Plot to compare variables over time ----
# Credit to Thomas Crow who used Gemini 2.5 pro to assist
# Note that there is an odd parsing issue where you need to highlight the whole function or run it from the final curly brace in order for it to be interpeted as a complete command

# You can add an arbitrary number of variables AND, if you want, add a “grouping variable” (e.g., treatment group)
# which will separate each plotted variable into (at least) two lines, one per level of the grouping variable.
# It can also divide up time in different ways (using the “cut_” series of functions from ggplot2).

#' Plot Multiple Variables Over Time, Optionally Grouped (Updated)
#'
#' Plots one or more variables from a data frame over a specified time variable,
#' optionally grouping the lines for each variable by the levels of a discrete
#' grouping variable. By default, plots against each unique value of the time
#' variable. Provides an option to bin the time variable if it is numeric.
#' Handles both continuous (plotting means with SE) and binary
#' (plotting proportions) outcomes, requiring all plotted variables (`...`) to be
#' of the same type. Offers options for data standardization.
#'
#' @param data A data frame or tibble containing the data.
#' @param time The bare (unquoted) name of the column representing the time variable.
#'        Can be numeric, Date, POSIXt, or factor (ordered or unordered).
#' @param ... Bare (unquoted) names of the columns representing the variables to plot.
#'        At least one variable must be provided. All variables must be of the same
#'        type (either all binary or all continuous/numeric).
#' @param grouping_var Optional. The bare (unquoted) name of a discrete column
#'        to group the lines by within each variable plotted in `...`.
#'        Must have 5 or fewer unique non-NA values. If more than 3 variables are
#'        plotted in `...`, `grouping_var` must be binary (exactly 2 unique non-NA values).
#'        Defaults to `NULL` (no grouping). Rows with NA in `grouping_var` are excluded.
#' @param bin_time Logical. If `TRUE`, bins a numeric `time` variable using the
#'        `cut_method` and `cut_into_ngroups`. If `FALSE` (Default), plots against
#'        each unique value of the `time` variable. Ignored if `time` is not numeric.
#' @param cut_method Character string. Method for binning numeric time variable if
#'        `bin_time = TRUE`. Accepts "interval" (Default), "number", "width".
#'        "Width" behaves poorly.
#' @param cut_into_ngroups Integer. Target number of bins for time variable cutting
#'        if `bin_time = TRUE`. Defaults to 5. For `cut_method = "width"`, this
#'        value specifies the width rather than the number of groups.
#' @param time_cut_threshold Integer. If `bin_time = FALSE` and the number of unique
#'        numeric time points exceeds this threshold, a warning about potential
#'        plot clutter is issued. Defaults to 20. Ignored if `bin_time = TRUE`.
#' @param standardize Logical. If TRUE, standardizes continuous variables before plotting.
#'        Standardization (overall, first_value, first_bin) is calculated based on the
#'        data for each variable *ignoring* the `grouping_var` to provide a consistent
#'        baseline across groups. See `standardize_ref`. Ignored for binary variables.
#'        Defaults to FALSE.
#' @param standardize_ref Character string specifying the reference for standardization
#'        if `standardize = TRUE`. Accepts: "overall" (Default), "first_value", "first_bin".
#'        Note: "first_bin" requires `bin_time = TRUE` and a numeric time variable.
#'        Fallback logic to "overall" applies if reference method fails.
#' @param warn_diff_scale_ratio Numeric. If plotting continuous variables without standardization,
#'        a warning is issued if the ratio of the largest range to the smallest range among
#'        the plotted variables exceeds this value. Set to `Inf` to disable. Defaults to 10.
#'
#' @return A ggplot object.
#'
#' @import dplyr
#' @import ggplot2
#' @import rlang
#' @import tidyr
#' @import forcats
#' @import scales
#' @import readr
#'
#' @examples
#' \dontrun{
#' # --- Setup Example Data ---
#' set.seed(123)
#' time_points_few <- 1:10 # Fewer time points for default plotting
#' time_points_many <- 1:30 # More time points to trigger warning/binning need
#' n_obs_per_time <- 20
#' data_ex_few <- data.frame(
#'   week = rep(time_points_few, each = n_obs_per_time),
#'   subject_id = 1:(length(time_points_few) * n_obs_per_time),
#'   treatment = sample(c("Placebo", "DrugX"),
#'                      length(time_points_few) * n_obs_per_time, replace = TRUE),
#'   site = sample(paste0("Site ", LETTERS[1:3]),
#'                  length(time_points_few) * n_obs_per_time, replace = TRUE)
#' ) %>%
#'   mutate(
#'     biomarker1 = 10 + 0.5 * week + rnorm(n(), 0, 2) + ifelse(treatment == "DrugX", 5, 0),
#'     biomarker2 = 100 - 1 * week + rnorm(n(), 0, 10) + ifelse(treatment == "DrugX", -20, 0),
#'     responder_f = factor(rnorm(n()) > (0 - 0.1*week + ifelse(treatment=="DrugX", -0.5, 0)),
#'                          levels=c(FALSE,TRUE),labels=c("No","Yes"))
#'   )
#' data_ex_many <- data.frame(
#'   day = rep(time_points_many, each = n_obs_per_time),
#'   subject_id = 1:(length(time_points_many) * n_obs_per_time),
#'   treatment = sample(c("Placebo", "DrugX"),
#'                      length(time_points_many) * n_obs_per_time, replace = TRUE)
#' ) %>%
#' mutate(
#'     biomarker1 = 10 + 0.5 * day + rnorm(n(), 0, 2) + ifelse(treatment == "DrugX", 5, 0),
#'     biomarker2 = 100 - 1 * day + rnorm(n(), 0, 10) + ifelse(treatment == "DrugX", -20, 0)
#' )
#'
#' # --- Examples ---
#'
#' # 1. Default: Plot unique time points (week has 10 unique values)
#' plot_vars_over_time(data_ex_few, week, biomarker1, biomarker2, grouping_var = treatment)
#'
#' # 2. Default with many time points: > time_cut_threshold (30 > 20), issues warning
#' plot_vars_over_time(data_ex_many, day, biomarker1, biomarker2, grouping_var = treatment)
#'
#' # 3. Explicit binning: Use bin_time = TRUE for many time points
#' plot_vars_over_time(data_ex_many, day, biomarker1, biomarker2, grouping_var = treatment,
#'                             bin_time = TRUE, cut_into_ngroups = 6)
#'
#' # 4. Binary variables, default plotting (unique time points)
#' plot_vars_over_time(data_ex_few, week, responder_f, grouping_var = treatment)
#'
#' # 5. Binary variables, explicit binning
#' plot_vars_over_time(data_ex_few, week, responder_f, grouping_var = treatment,
#'                             bin_time = TRUE, cut_into_ngroups = 3)
#'
#' # 6. Standardization with default unique time points
#' plot_vars_over_time(data_ex_few, week, biomarker1, biomarker2, grouping_var = treatment,
#'                         standardize = TRUE, standardize_ref = "first_value")
#'
#' # 7. Standardization with binned time
#' plot_vars_over_time(data_ex_many, day, biomarker1, biomarker2, grouping_var = treatment,
#'                        bin_time = TRUE, standardize = TRUE, standardize_ref = "first_bin")
#' # Try first_bin without bin_time = TRUE (should warn and fallback)
#' plot_vars_over_time(data_ex_many, day, biomarker1, biomarker2, grouping_var = treatment,
#'                            standardize = TRUE, standardize_ref = "first_bin")
#'
#' }
plot_vars_over_time <- function(data,
                                time,
                                ..., # Variables to plot
                                grouping_var = NULL, # Optional grouping variable
                                bin_time = FALSE, # NEW: control binning explicitly
                                cut_method = "interval",
                                cut_into_ngroups = 5,
                                time_cut_threshold = 20, # Threshold for WARNING when bin_time=FALSE
                                standardize = FALSE,
                                standardize_ref = "overall", # Accepts "overall", "first_value", "first_bin"
                                warn_diff_scale_ratio = 10
) {

  ## Load required libraries ----
  # Assuming loaded via @import or library() calls elsewhere

  ## Capture variable names using NSE ----
  time_quo <- enquo(time)
  vars_quos <- enquos(...)
  grouping_var_quo <- enquo(grouping_var) # Capture grouping var

  time_name <- quo_name(time_quo)
  vars_names <- map_chr(vars_quos, quo_name)
  grouping_var_name <- if (!quo_is_null(grouping_var_quo)) quo_name(grouping_var_quo) else NULL

  if (length(vars_names) == 0) {
    stop("You must provide at least one variable name in ... to plot.", call. = FALSE)
  }

  ## Input checks ----
  if (!is.data.frame(data)) {
    stop("Argument 'data' must be a data frame or tibble.", call. = FALSE)
  }
  required_vars <- c(time_name, vars_names)
  if (!is.null(grouping_var_name)) {
    required_vars <- c(required_vars, grouping_var_name)
  }
  missing_vars <- setdiff(required_vars, names(data))
  if (length(missing_vars) > 0) {
    stop("The following variables are not found in the data: ",
         paste(missing_vars, collapse = ", "), call. = FALSE)
  }
  # Check new bin_time parameter
  if (!is.logical(bin_time) || length(bin_time) != 1) {
    stop("Argument 'bin_time' must be TRUE or FALSE.", call. = FALSE)
  }
  # Other checks remain mostly the same
  if (!is.logical(standardize) || length(standardize) != 1) {
    stop("Argument 'standardize' must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.numeric(warn_diff_scale_ratio) || length(warn_diff_scale_ratio) != 1 || warn_diff_scale_ratio <= 0) {
    stop("Argument 'warn_diff_scale_ratio' must be a single positive number.", call. = FALSE)
  }
  if (!is.character(cut_method) || length(cut_method) != 1 || !cut_method %in% c("interval", "number", "width")) {
    stop("Argument 'cut_method' must be one of 'interval', 'number', or 'width'.", call. = FALSE)
  }
  if (!is.numeric(cut_into_ngroups) || length(cut_into_ngroups) != 1 || cut_into_ngroups <= 0 || floor(cut_into_ngroups) != cut_into_ngroups) {
    stop("Argument 'cut_into_ngroups' must be a single positive integer.", call. = FALSE)
  }
  if (!is.numeric(time_cut_threshold) || length(time_cut_threshold) != 1 || time_cut_threshold <= 0 || floor(time_cut_threshold) != time_cut_threshold) {
    stop("Argument 'time_cut_threshold' must be a single positive integer.", call. = FALSE)
  }
  valid_std_refs <- c("overall", "first_value", "first_bin")
  if (!standardize_ref %in% valid_std_refs) {
    stop("Argument 'standardize_ref' must be one of: ",
         paste(valid_std_refs, collapse=", "), ".", call. = FALSE)
  }
  if(standardize_ref %in% c("first_value", "first_bin") && standardize == FALSE) {
    # This is just informational, not an error
    # message("Note: 'standardize_ref = \"", standardize_ref, "\"' is ignored because 'standardize = FALSE'.")
  }
  # Check dependency for first_bin standardization
  if(standardize && standardize_ref == "first_bin" && !bin_time) {
    message("Note: 'standardize_ref = \"first_bin\"' requires 'bin_time = TRUE'. Standardization may fallback to 'overall' if binning doesn't occur or fails.")
  }


  # Grouping variable specific checks (remain same)
  n_grouping_levels <- NA
  if (!is.null(grouping_var_name)) {
    grouping_col_vector <- data[[grouping_var_name]]
    n_grouping_levels <- n_distinct(grouping_col_vector, na.rm = TRUE)

    if (n_grouping_levels > 5) {
      stop("Argument 'grouping_var' ('", grouping_var_name,
           "') must be discrete with 5 or fewer unique non-NA values. Found: ",
           n_grouping_levels, ".", call. = FALSE)
    }
    if(n_grouping_levels < 1) {
      stop("Argument 'grouping_var' ('", grouping_var_name,
           "') appears to have no non-NA values after accounting for NAs in time.", call. = FALSE)
    }
    if (length(vars_names) > 3 && n_grouping_levels > 2) {
      stop("When plotting more than 3 variables (", length(vars_names),
           " requested), 'grouping_var' must be binary (exactly 2 unique non-NA levels) ",
           "to avoid excessive plot clutter. '", grouping_var_name,
           "' has ", n_grouping_levels, " levels.", call. = FALSE)
    }
    message("Grouping plot by variable: '", grouping_var_name, "' (", n_grouping_levels, " levels).")
  }


  ## Select necessary columns and initial processing ----
  select_quos <- c(
    time_col = time_quo,
    vars_quos,
    if (!is.null(grouping_var_name)) rlang::list2(!!grouping_var_name := grouping_var_quo) else NULL
  )

  plot_data_wide <- data %>%
    ungroup() %>%
    select(!!!select_quos) %>%
    filter(!is.na(time_col))

  # Filter out NAs in grouping variable if it exists
  if (!is.null(grouping_var_name)){
    initial_rows <- nrow(plot_data_wide)
    plot_data_wide <- plot_data_wide %>%
      filter(!is.na(.data[[grouping_var_name]]))
    removed_rows <- initial_rows - nrow(plot_data_wide)
    if(removed_rows > 0){
      message("Removed ", removed_rows, " rows with NA in the grouping variable '", grouping_var_name, "'.")
    }
    # Re-check grouping levels after filtering NAs in time_col (edge case)
    if (nrow(plot_data_wide) > 0) {
      grouping_col_vector <- plot_data_wide[[grouping_var_name]]
      n_grouping_levels_after_time_filter <- n_distinct(grouping_col_vector, na.rm = TRUE)
      if(n_grouping_levels_after_time_filter < 1) {
        stop("Argument 'grouping_var' ('", grouping_var_name,
             "') has no non-NA values remaining after filtering NA in time variable '", time_name, "'.", call. = FALSE)
      }
      # Convert grouping var to factor to ensure consistent levels/order in plot
      plot_data_wide <- plot_data_wide %>%
        mutate(!!grouping_var_name := factor(.data[[grouping_var_name]]))
    } else {
      stop("No valid data remaining after filtering NAs in time variable '", time_name, "'.", call. = FALSE)
    }
  }


  if(nrow(plot_data_wide) == 0) {
    stop("No valid data remaining after filtering NAs in time ",
         if(!is.null(grouping_var_name)) paste0("and '", grouping_var_name, "'."),
         call. = FALSE)
  }

  ## Handle Time Variable Processing ----
  binning_info_string <- ""
  time_was_binned <- FALSE # Flag to track if binning actually occurred
  plot_time_label <- time_name # Default label

  time_col_vector <- plot_data_wide$time_col # Get the vector once
  time_is_numeric <- is.numeric(time_col_vector)
  n_distinct_time <- n_distinct(time_col_vector) # Use the vector

  # Determine if binning *should* be applied based on parameter and type
  apply_binning_condition <- bin_time && time_is_numeric

  if (apply_binning_condition) {
    ### Apply binning if requested and time is numeric ----
    message("Time variable '", time_name, "' is numeric and bin_time=TRUE.",
            "\nBinning into ", cut_into_ngroups, " groups using cut_", cut_method, "().")

    binning_info_string <- paste0("\nTime variable '", time_name, "' binned into ",
                                  cut_into_ngroups, " groups using cut_", cut_method, "().")
    time_was_binned <- TRUE # Set the flag

    # The cut_* functions handle cases with few unique values gracefully
    plot_data_wide <- plot_data_wide %>%
      mutate(
        time_plot = switch(cut_method,
                           "interval" = cut_interval(time_col, n = cut_into_ngroups, dig.lab = 4),
                           "number"   = cut_number(time_col, n = cut_into_ngroups, dig.lab = 4),
                           "width"    = cut_width(time_col, width = cut_into_ngroups, dig.lab = 4)
        )
      )
    plot_time_label <- paste0(time_name, " (Binned)")

  } else {
    ### Plot unique time points (default or time not numeric) ----
    time_was_binned <- FALSE # Ensure flag is false

    # Issue warning if plotting many unique numeric points without binning
    if (!bin_time && time_is_numeric && n_distinct_time > time_cut_threshold) {
      warning("Time variable '", time_name, "' has ", n_distinct_time,
              " unique numeric values (> threshold of ", time_cut_threshold,"). Plotting all points.",
              "\n  Consider setting `bin_time = TRUE` if the plot is cluttered or slow.", call. = FALSE)
    } else if (bin_time && !time_is_numeric) {
      message("Note: `bin_time = TRUE` ignored because time variable '", time_name, "' is not numeric. Plotting unique values/levels.")
    }

    # Handle factor time variables (reordering if unordered and possible)
    time_is_factor <- is.factor(time_col_vector)
    time_is_ordered <- is.ordered(time_col_vector)

    if (time_is_factor && !time_is_ordered) {
      message("Time variable '", time_name, "' is an unordered factor. Attempting to determine logical order...")
      original_levels <- levels(time_col_vector)
      numeric_levels <- suppressWarnings(as.numeric(original_levels))

      # Try different ordering methods
      ordered_levels <- NULL
      if (!any(is.na(numeric_levels))) {
        message("... Ordering factor levels numerically.")
        ordered_levels <- original_levels[order(numeric_levels)]
      } else {
        num_extract <- suppressWarnings(readr::parse_number(original_levels))
        if (!any(is.na(num_extract)) && length(unique(num_extract)) == length(original_levels)) {
          message("... Ordering factor levels based on extracted numbers.")
          ordered_levels <- original_levels[order(num_extract)]
        } else {
          # Try Date parsing as last resort
          date_levels <- suppressWarnings(as.Date(original_levels)) # Try common formats
          if (!any(is.na(date_levels))) {
            message("... Ordering factor levels as Dates.")
            ordered_levels <- original_levels[order(date_levels)]
          } else {
            message("... Could not determine a clear numeric/date order. Using existing factor level order (or alphabetical if none specified).")
          }
        }
      }
      # Apply ordering if found
      if (!is.null(ordered_levels)) {
        plot_data_wide <- plot_data_wide %>%
          mutate(time_col = factor(time_col, levels = ordered_levels, ordered = TRUE))
      }
      plot_data_wide <- plot_data_wide %>% mutate(time_plot = time_col)
    } else {
      # Assign original time column to time_plot for numeric/date/posix/ordered factors/character
      plot_data_wide <- plot_data_wide %>% mutate(time_plot = time_col)
    }
    # plot_time_label remains time_name (set as default above)
  }

  # Determine if the resulting time_plot axis is discrete AFTER potential processing
  # Character time cols will fall through to here and be treated as discrete.
  is_discrete_time_axis <- is.factor(plot_data_wide$time_plot) || is.character(plot_data_wide$time_plot)
  # Determine if axis text might need angling
  needs_angling <- is_discrete_time_axis || time_was_binned || (!bin_time && time_is_numeric && n_distinct_time > time_cut_threshold)


  ## Determine the first time point VALUE for potential standardization ----
  # (Logic remains identical, uses plot_data_wide$time_col - the original time)
  first_time_val <- NA
  if (standardize && standardize_ref == "first_value") {
    time_col_for_first_val <- plot_data_wide$time_col # Use original time
    if (is.factor(time_col_for_first_val)) {
      # Use the potentially re-ordered levels if it's a factor
      first_time_val <- levels(time_col_for_first_val)[1]
    } else if (is.numeric(time_col_for_first_val) || inherits(time_col_for_first_val, "Date") || inherits(time_col_for_first_val, "POSIXt")) {
      first_time_val <- min(time_col_for_first_val, na.rm = TRUE)
    } else {
      warning("Cannot determine first time point for 'first_value' standardization (time type not supported). Defaulting to 'overall'.", call. = FALSE)
      standardize_ref <- "overall"
    }
    if(is.na(first_time_val) || (is.numeric(first_time_val) && is.infinite(first_time_val)) ){
      warning("Could not determine a valid first time point value for 'first_value' standardization. Defaulting to 'overall'.", call. = FALSE)
      standardize_ref <- "overall"
    }
  }


  ## Pivot data to long format ----
  pivot_cols <- vars_names
  plot_data_long <- plot_data_wide %>%
    pivot_longer(
      cols = all_of(pivot_cols),
      names_to = "variable",
      values_to = "value"
    )

  ## Determine Variable Types ----
  # (Logic remains identical)
  var_types <- plot_data_long %>%
    # Ensure NAs in value don't break type detection completely
    filter(!is.na(value)) %>%
    group_by(variable) %>%
    # Use first non-NA value to guess base type robustly + N distinct
    summarise(
      n_distinct_val = n_distinct(value, na.rm = TRUE),
      first_val = first(value), # Get an example value
      .groups = "drop"
    ) %>%
    # Add base type check more directly
    mutate(
      is_num = is.numeric(first_val),
      is_logi = is.logical(first_val),
      # Add check if it's factor/char explicitly, useful if numbers were stored as char
      is_fac_char = is.factor(first_val) || is.character(first_val)
    )

  # Refine categorization: Ensure numeric/continuous vars are actually numeric/logical
  all_binary <- all(var_types$n_distinct_val == 2)
  all_continuous <- all(var_types$n_distinct_val > 2 & (var_types$is_num | var_types$is_logi))
  vars_too_few_levels <- var_types$variable[var_types$n_distinct_val < 2]

  if (length(vars_too_few_levels) > 0) {
    warning("The following variables have fewer than 2 distinct non-NA values: ",
            paste(vars_too_few_levels, collapse = ", "), ". Check data. These may cause issues or empty plots.", call. = FALSE)
  }

  # Stricter check: If supposed to be binary, ensure not numeric like 1, 5
  if (all_binary) {
    non_logical_binary <- var_types %>% filter(n_distinct_val == 2 & !is_logi & !is_fac_char)
    if(nrow(non_logical_binary) > 0){
      message("Note: Variables treated as binary include: ", paste(non_logical_binary$variable, collapse=", "),". Ensure these represent two categories.")
    }
  }
  # Stricter check: If supposed to be continuous, ensure not factor/char
  if (all_continuous) {
    non_numeric_continuous <- var_types %>% filter(n_distinct_val > 2 & is_fac_char)
    if(nrow(non_numeric_continuous) > 0){
      stop("Plotting failed. Variables intended as continuous appear to be factors or characters: ",
           paste(non_numeric_continuous$variable, collapse=", "), ". Please convert them to numeric.", call. = FALSE)
    }
  }


  ## Choose Plot Type & Process ----

  plot_title <- if(is.null(grouping_var_name)) "Comparison Over Time" else paste0("Comparison Over Time by '", grouping_var_name, "'")
  gg <- NULL # Initialize plot object

  if (all_binary) {
    ### Plotting Proportions for Binary Variables ----
    # (Logic remains largely same, uses 'time_plot' correctly)
    message(paste0("All variables appear binary. Plotting proportions",
                   if (!is.null(grouping_var_name)) paste0(" grouped by '", grouping_var_name, "'.") else "."))
    if(standardize) message("Note: 'standardize = TRUE' is ignored for binary variables.")

    # Determine reference level (ensure robust for various binary types)
    get_ref_level <- function(v) {
      unique_vals <- unique(na.omit(v))
      if (length(unique_vals) != 2) return(NA)

      # Attempt to order logically: FALSE < TRUE, 0 < 1, "No" < "Yes" etc.
      if(is.logical(unique_vals)) return(TRUE)
      if(is.numeric(unique_vals)) return(max(unique_vals))

      # Try common factor level pairings (case-insensitive)
      lvls_lower <- tolower(as.character(unique_vals))
      if(all(c("no", "yes") %in% lvls_lower)) return(unique_vals[which(lvls_lower=="yes")[1]])
      if(all(c("false", "true") %in% lvls_lower)) return(unique_vals[which(lvls_lower=="true")[1]])
      if(all(c("0", "1") %in% lvls_lower)) return(unique_vals[which(lvls_lower=="1")[1]])

      # Fallback: sort and take second element (may not be semantically meaningful)
      # Convert to character for consistent sorting
      return(sort(as.character(unique_vals))[2])
    }

    ref_levels_list <- plot_data_long %>%
      group_by(variable) %>%
      summarise(ref_level = get_ref_level(value), .groups="drop") %>%
      filter(!is.na(ref_level)) # Filter NAs if get_ref_level failed

    if(nrow(ref_levels_list) == 0){
      stop("Could not determine a valid binary reference level for any variable.", call.=FALSE)
    } else if(nrow(ref_levels_list) != length(vars_names)){
      failed_vars <- setdiff(vars_names, ref_levels_list$variable)
      warning("Could not determine binary reference level for variable(s): ", paste(failed_vars, collapse=", "), ". These variable(s) will be excluded.", call. = FALSE)
      # Proceed only with variables where ref level was found
      plot_data_long <- plot_data_long %>% filter(variable %in% ref_levels_list$variable)
      if(nrow(plot_data_long) == 0) stop("No data remaining after excluding variables without determinable reference level.", call.=FALSE)
    }

    # Define grouping columns for summary
    summary_group_cols <- c("time_plot", "variable", "ref_level", if (!is.null(grouping_var_name)) grouping_var_name)

    # Summarize data
    summary_data <- plot_data_long %>%
      inner_join(ref_levels_list, by = "variable") %>%
      mutate(value_char = as.character(value), ref_level_char = as.character(ref_level)) %>% # Compare as chars for robustness
      group_by(across(all_of(summary_group_cols))) %>%
      summarise(
        proportion = mean(value_char == ref_level_char, na.rm = TRUE),
        n = sum(!is.na(value)),
        .groups = "drop"
      ) %>%
      mutate(
        proportion = if_else(n == 0 | is.nan(proportion), NA_real_, proportion),
        variable_label = paste0(variable, " (", ref_level, ")")
      )

    subtitle_prop <- paste0("Showing proportion of the reference level (typically second unique level) for each variable.",
                            binning_info_string)

    # Base plot
    gg <- ggplot(summary_data, aes(x = time_plot, y = proportion, color = variable_label))

    # Add grouping aesthetics
    if (!is.null(grouping_var_name)) {
      gg <- gg + aes(linetype = .data[[grouping_var_name]],
                     group = interaction(variable, .data[[grouping_var_name]]))
    } else {
      gg <- gg + aes(group = variable)
    }

    gg <- gg +
      geom_line(linewidth = 1, na.rm = TRUE) +
      geom_point(size = 2.5, na.rm = TRUE) +
      scale_y_continuous(limits = c(0, 1), labels = scales::percent_format()) +
      labs(
        title = plot_title,
        subtitle = subtitle_prop,
        x = plot_time_label,
        y = "Proportion",
        color = "Variable (Reference Level)",
        linetype = if(!is.null(grouping_var_name)) grouping_var_name else NULL
      ) +
      theme_minimal(base_size = 12) +
      theme(axis.text.x = element_text(angle = if(needs_angling) 45 else 0,
                                       hjust = if(needs_angling) 1 else 0.5),
            legend.position = "bottom",
            legend.box="vertical")


  } else if (all_continuous) {
    ### Plotting Means for Continuous Variables ----
    # (Logic remains largely same, uses 'time_plot' correctly)
    message(paste0("All variables appear continuous/numeric. Plotting means with SE bars",
                   if (!is.null(grouping_var_name)) paste0(" grouped by '", grouping_var_name, "'.") else "."))

    # Check scale differences (as before)
    var_ranges <- plot_data_long %>%
      group_by(variable) %>%
      summarise(range_val = diff(range(as.numeric(value), na.rm = TRUE)), .groups = "drop") %>%
      filter(is.finite(range_val) & range_val > 1e-9) # Use 1e-9 to avoid division by almost zero

    if(nrow(var_ranges) > 1 && !is.infinite(warn_diff_scale_ratio)) {
      min_range <- min(var_ranges$range_val)
      max_range <- max(var_ranges$range_val)
      if (min_range > 1e-9) { # Avoid division by zero if range is tiny
        scale_ratio <- max_range / min_range
        if (scale_ratio > warn_diff_scale_ratio) {
          warning(sprintf("Ranges differ significantly (max/min ratio = %.1f > %.1f). Consider `standardize = TRUE`.", scale_ratio, warn_diff_scale_ratio), call. = FALSE)
        }
      }
    }

    y_label <- "Mean Value"
    plot_subtitle <- "Points represent mean, error bars represent +/- 1 Standard Error (SE)."
    standardization_details <- ""

      #### Apply standardization if requested ----
    if (standardize) {
      original_standardize_ref <- standardize_ref
      standardization_applied_variables <- character(0)
      apply_standardization <- FALSE

      # Prepare data for standardization calc (needs original time for 'first_value', potentially binned for 'first_bin')
      data_for_std_calc <- plot_data_long %>%
        mutate(original_time_col_for_std = time_col,
               plot_time_for_std = time_plot) # Keep both time versions

      actual_standardize_method <- standardize_ref
      fallback_triggered <- FALSE

      #### Reference-based standardization ----
      if (actual_standardize_method %in% c("first_bin", "first_value")) {
        ref_time_description <- ""
        ref_filter_expr <- NULL

        if (actual_standardize_method == "first_bin") {
          # Use the time_was_binned flag here
          if (!time_was_binned) {
            warning("'standardize_ref = \"first_bin\"' requires successful time binning (use `bin_time = TRUE` with numeric time). Falling back to 'overall'.", call. = FALSE)
            actual_standardize_method <- "overall"; fallback_triggered <- TRUE
          } else {
            # Ensure time_plot is factor (result of cut_*)
            if (!is.factor(data_for_std_calc$plot_time_for_std)) {
              warning("Cannot apply 'first_bin' standardization because the binned time variable is not a factor. Falling back to 'overall'.", call.=FALSE)
              actual_standardize_method <- "overall"; fallback_triggered <- TRUE
            } else {
              first_bin_level <- levels(data_for_std_calc$plot_time_for_std)[1]
              if(is.na(first_bin_level)) {
                warning("Could not determine first bin level for 'first_bin' standardization. Falling back to 'overall'.", call.=FALSE)
                actual_standardize_method <- "overall"; fallback_triggered <- TRUE
              } else {
                ref_time_description <- paste("first bin '", format(first_bin_level), "'", sep="")
                ref_filter_expr <- rlang::expr(plot_time_for_std == !!first_bin_level)
              }
            }
          }
        } else { # first_value
          if(is.na(first_time_val) || (is.numeric(first_time_val) && is.infinite(first_time_val))){
            warning("Invalid or missing original first time value for 'first_value' std. Falling back to 'overall'.", call. = FALSE)
            actual_standardize_method <- "overall"; fallback_triggered <- TRUE
          } else {
            ref_time_description <- paste("first time value '", format(first_time_val), "'", sep="")
            # Compare carefully: allow for floating point inaccuracy if numeric, exact match if factor/date
            if(is.numeric(first_time_val) && is.numeric(data_for_std_calc$original_time_col_for_std)){
              # Use tolerance for numeric comparison
              time_tolerance <- .Machine$double.eps ^ 0.5
              ref_filter_expr <- rlang::expr(abs(original_time_col_for_std - !!first_time_val) < !!time_tolerance)
            } else {
              # Expect exact match for Date, POSIXt, factor, character
              ref_filter_expr <- rlang::expr(original_time_col_for_std == !!first_time_val)
            }
          }
        }

        if (!fallback_triggered && !is.null(ref_filter_expr)) {
          ref_stats <- data_for_std_calc %>%
            filter(!!ref_filter_expr) %>%
            group_by(variable) %>%
            summarise(ref_mean = mean(as.numeric(value), na.rm = TRUE),
                      ref_sd = sd(as.numeric(value), na.rm = TRUE),
                      ref_n = sum(!is.na(value) & is.finite(as.numeric(value))), # Count only valid numbers for ref calc
                      .groups = "drop") %>%
            mutate(is_valid_ref = ref_n > 1 & !is.na(ref_sd) & ref_sd > 1e-9 & is.finite(ref_mean)) # sd > 0 essential

          valid_ref_vars <- ref_stats %>% filter(is_valid_ref)
          invalid_ref_vars <- setdiff(vars_names, valid_ref_vars$variable)

          if (length(invalid_ref_vars) > 0) {
            warning("Cannot standardize relative to ", ref_time_description, " for: ", paste(invalid_ref_vars, collapse=", "), ". Insufficient valid data, zero variance, or non-finite mean at reference point. Falling back to 'overall' for these.", call. = FALSE)
            # Selectively apply? Or fallback all? Let's fallback all for simplicity here. Could be refined to partial.
            actual_standardize_method <- "overall"; fallback_triggered <- TRUE
            # Alternative: Standardize only valid_ref_vars? More complex merging later.
          } else if (nrow(valid_ref_vars) == 0) {
            warning("Cannot standardize relative to ", ref_time_description, " (no variables had valid reference data). Falling back to 'overall'.", call. = FALSE)
            actual_standardize_method <- "overall"; fallback_triggered <- TRUE
          } else {
            message("Applying standardization relative to ", ref_time_description, ".")
            plot_data_long <- plot_data_long %>%
              left_join(select(ref_stats, variable, ref_mean, ref_sd), by = "variable") %>%
              mutate(value = (as.numeric(value) - ref_mean) / ref_sd) %>%
              select(-ref_mean, -ref_sd)
            y_label <- paste0("Value (Std. ref: ", ref_time_description, ")")
            standardization_details <- paste0("Std: ", ref_time_description) # Short form for subtitle
            apply_standardization <- TRUE; standardization_applied_variables <- vars_names
          }
        }
      } # End reference-based block


      #### Overall standardization (if primary or fallback) ----
      if (actual_standardize_method == "overall") {
        if (original_standardize_ref == "overall" && !fallback_triggered) {
          message("Standardizing variables overall (Z-scores based on all time points).")
        } else if (fallback_triggered) {
          message("Applying 'overall' standardization due to fallback.")
        }

        overall_stats <- data_for_std_calc %>%
          group_by(variable) %>%
          summarise(overall_mean = mean(as.numeric(value), na.rm = TRUE),
                    overall_sd = sd(as.numeric(value), na.rm = TRUE),
                    n_total = sum(!is.na(value) & is.finite(as.numeric(value))),
                    .groups = "drop") %>%
          mutate(can_standardize = n_total > 1 & !is.na(overall_sd) & overall_sd > 1e-9 & is.finite(overall_mean))

        can_std_vars <- overall_stats %>% filter(can_standardize)
        cannot_std_vars <- setdiff(vars_names, can_std_vars$variable)

        if (length(cannot_std_vars) > 0) {
          warning("Cannot apply overall standardization for: ", paste(cannot_std_vars, collapse=", "), ". Insufficient data or zero variance. These variables remain unstandardized.", call. = FALSE)
          if (nrow(can_std_vars) > 0) {
            # Standardize the ones we can
            plot_data_long <- plot_data_long %>%
              left_join(select(can_std_vars, variable, overall_mean, overall_sd), by = "variable") %>%
              # Only apply std if overall_sd is available (from the join)
              mutate(value = ifelse(!is.na(overall_sd), (as.numeric(value) - overall_mean) / overall_sd, as.numeric(value))) %>%
              select(-overall_mean, -overall_sd)
            y_label <- "Mean Value (Partial Overall Std.)"; standardization_details <- "Std: Overall, partial"
            apply_standardization <- TRUE; standardization_applied_variables <- can_std_vars$variable
          } else { # None could be standardized
            warning("Overall standardization skipped entirely for all variables.", call.=FALSE)
            y_label <- "Mean Value"; standardization_details <- "Std: Failed"
            apply_standardization <- FALSE # Ensure flag is correct
          }
        } else { # Standardize all
          plot_data_long <- plot_data_long %>%
            left_join(select(overall_stats, variable, overall_mean, overall_sd), by = "variable") %>%
            mutate(value = (as.numeric(value) - overall_mean) / overall_sd) %>%
            select(-overall_mean, -overall_sd)
          y_label <- "Standardized Mean Value (Z-score)"; standardization_details <- "Std: Overall"
          apply_standardization <- TRUE; standardization_applied_variables <- vars_names
        }
      } # End overall standardisation block

      # Update subtitle based on final outcome
      if(apply_standardization) plot_subtitle <- paste0(plot_subtitle, " (", standardization_details, ")")
      else if(standardize) plot_subtitle <- paste0(plot_subtitle, " (Std. Failed)") # If requested but failed
    } # End if(standardize)

    # Summarize data (means, SE) - ensure value is numeric after potential std
    summary_group_cols <- c("time_plot", "variable", if (!is.null(grouping_var_name)) grouping_var_name)
    summary_data <- plot_data_long %>%
      mutate(value = as.numeric(value)) %>% # Ensure numeric before summary
      filter(is.finite(value)) %>% # Remove Inf/-Inf from summary calcs
      group_by(across(all_of(summary_group_cols))) %>%
      summarise(
        mean = mean(value, na.rm = TRUE),
        sd = sd(value, na.rm = TRUE),
        n = sum(!is.na(value)), # Count non-NA finite values
        .groups = "drop"
      ) %>%
      mutate(
        se = sd / sqrt(n),
        # Handle cases with n=0 or n=1 where mean/se are NA/NaN/Inf
        se = if_else(n > 1 & !is.na(se) & is.finite(se), se, NA_real_),
        mean = if_else(n > 0 & !is.na(mean) & is.finite(mean), mean, NA_real_)
      ) %>%
      filter(!is.na(mean)) # Remove points where mean could not be computed

    if(nrow(summary_data) == 0) {
      stop("No valid summary data could be calculated (e.g., all values NA, Inf, or insufficient N). Cannot generate plot.", call. = FALSE)
    }

    # Create the plot
    gg <- ggplot(summary_data, aes(x = time_plot, y = mean, color = variable))

    # Add grouping aesthetics
    if (!is.null(grouping_var_name)) {
      gg <- gg + aes(linetype = .data[[grouping_var_name]],
                     shape = .data[[grouping_var_name]], # Add shape mapping
                     group = interaction(variable, .data[[grouping_var_name]]))
    } else {
      gg <- gg + aes(group = variable)
    }

    # Calculate error bar Width based on whether axis is discrete/binned
    error_bar_width <- if(is_discrete_time_axis || time_was_binned) 0.1 else 0

    gg <- gg +
      geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
                    width = error_bar_width,
                    linewidth = 0.7,
                    na.rm = TRUE) + # Error bars first, behind line/point
      geom_line(linewidth = 1, na.rm = TRUE) +
      geom_point(size = 2.5, na.rm = TRUE) + # Apply shape here, works even if shape=NULL
      labs(
        title = plot_title,
        subtitle = paste0(plot_subtitle, binning_info_string),
        x = plot_time_label,
        y = y_label,
        color = "Variable",
        linetype = if(!is.null(grouping_var_name)) grouping_var_name else NULL,
        shape = if(!is.null(grouping_var_name)) grouping_var_name else NULL
      ) +
      theme_minimal(base_size = 12) +
      theme(axis.text.x = element_text(angle = if(needs_angling) 45 else 0,
                                       hjust = if(needs_angling) 1 else 0.5),
            legend.position = "bottom",
            legend.box = "vertical")

  } else {
    ### Mixed or invalid types ----
    type_summary <- var_types %>%
      mutate(type = case_when(
        n_distinct_val == 2 ~ "Binary",
        n_distinct_val > 2 & (is_num | is_logi) ~ "Continuous",
        is_fac_char ~ "Factor/Character", # Catch factors/chars explicitly
        TRUE ~ "Other/Insufficient Non-NA Levels"
      )) %>% select(variable, type) %>% tidyr::unite("summary", variable, type, sep=": ") %>% pull(summary) %>% paste(collapse=", ")

    stop("Plotting failed. All variables provided via '...' must be either binary (e.g., logical, factor/char with 2 levels) or numeric/continuous.\n",
         "Detected variable types: ", type_summary, call. = FALSE)
  }

  # Return the final plot
  if(is.null(gg)) {
    stop("Plotting failed unexpectedly after type checking.", call. = FALSE)
  }

  return(gg)
}
