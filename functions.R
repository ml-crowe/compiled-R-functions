####mean.n function####
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

#### Search Function ####
#found this function online
#example:varlist(data,pattern="hsns",exclude="5|oj")
#above example identifies all variable names that include 'hsns' but
#excludes those that have a 5 or "oj". Pipe key "|" is the "or" operator
#other operators (i.e,. "&") can also be used.
varlist <- function (df=NULL, pattern=NULL, exclude=NULL, type=c("numeric","factor","character", "double", "logical", "integer"), ignore.case=TRUE) {
  vars <- character(0)
  if (any(type %in% c("numeric", "double"))) {
    vars <- c(vars,names(df)[sapply(df,is.double)]) #had is.numeric and separate is.double, but that was resulting in duplicates
  }
  if (any(type %in% "factor")) {
    vars <- c(vars,names(df)[sapply(df,is.factor)])
  }
  if (any(type %in% "character")) {
    vars <- c(vars,names(df)[sapply(df,is.character)])
  }
  if (any(type %in% "logical")) {
    vars <- c(vars,names(df)[sapply(df,is.logical)])
  }
  if (any(type %in% "integer")) {
    vars <- c(vars,names(df)[sapply(df,is.integer)])
  }
  if(!is.null(exclude)){
    if(ignore.case==TRUE){
      list<-vars[(!vars %in% vars[grepl(vars,pattern=exclude,ignore.case=TRUE)]) & grepl(vars,pattern=pattern, ignore.case=TRUE)]
    }
    else{
      list<-vars[(!vars %in% vars[grepl(vars,pattern=exclude,ignore.case=FALSE)]) & grepl(vars,pattern=pattern, ignore.case=FALSE)]
    }
  }
  if(is.null(exclude)){
    if(ignore.case==TRUE){
      list<-vars[(!vars %in% exclude) & grepl(vars,pattern=pattern, ignore.case=TRUE)]
    }
    else{
      list<-vars[(!vars %in% exclude) & grepl(vars,pattern=pattern, ignore.case=FALSE)]
    }
  }
  list
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

search <- function(x = ''){
  varlist(df, pattern = x)
}

#####cor function####
#the skeleton of this function I got from a different student here. It gives flagged correlations.
#I modified a little so if you find it doesn't work for one purpose or another it is probably my fault.
#for the most part you can ignore this function because it just feeds into my r_table() function that is below.

cor_table<- function(x,vars=NULL,with=NULL,flag=TRUE,strict=FALSE,round = 2, type = 'pearson'){
  require(Hmisc,warn.conflicts=TRUE)
  data<-as.matrix(x)
  Rmat <- rcorr(data, type = type)
  RmatP <- Rmat$P
  ComRmat <- round(Rmat$r, 5)
  originalRmat<-Rmat$r
  ComRPmat <- round(Rmat$P,3)
  originalRPmat <- Rmat$P # added this line to avoid flagging issues due to rounding errors
  ComRnmat<- Rmat$n
  
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
  write.table(x,"clipboard-10000",sep="\t",na = '', row.names=row.names,col.names=col.names,...)
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

#### Syntax for taking vector output (e.g., result of call to names()) and convert to list separated by columns
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
#' If there are errors, the result retutrns is a list with the elements of warn and err.
#' This is a nice way to recover from a problems that may have occurred during loop evaluation or during cluster usage.
#' Check the references for additional related functions.
#' I have not included the other factory functions included in the original Stack Overflow answer because they did not play well with the return item as an S4 object.
#' @export
#' @param fun The function to be turned into a factory
#' @return The result of the function given to turn into a factory.  If this function was in error "An error as occurred" as a character element.  factory-error and factory-warning attributes may also be set as appropriate.
#' @references
#' \url{http://stackoverflow.com/questions/4948361/how-do-i-save-warnings-and-errors-as-output-from-a-function}
#' @author Martin Morgan; Modified by Russell S. Pierce
#' @examples 
#' f.log <- factory(log)
#' f.log("a")
#' f.as.numeric <- factory(as.numeric)
#' f.as.numeric(c("a","b",1))
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
