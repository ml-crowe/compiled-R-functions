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
varlist <- function (df=NULL,type=c("numeric","factor","character", "double", "logical", "integer"), pattern=NULL, exclude=NULL, ignore.case=TRUE) {
  vars <- character(0)
  if (any(type %in% "numeric")) {
    vars <- c(vars,names(df)[sapply(df,is.numeric)])
  }
  if (any(type %in% "factor")) {
    vars <- c(vars,names(df)[sapply(df,is.factor)])
  }
  if (any(type %in% "character")) {
    vars <- c(vars,names(df)[sapply(df,is.character)])
  }
  if (any(type %in% "double")) {
    vars <- c(vars,names(df)[sapply(df,is.double)])
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

cor_table<- function(x,vars=NULL,with=NULL,flag=TRUE,strict=FALSE,round = 2){
  require(Hmisc,warn.conflicts=TRUE)
  data<-as.matrix(x)
  Rmat <- rcorr(data)
  RmatP <- Rmat$P
  ComRmat <- round(Rmat$r, 5)
  originalRmat<-Rmat$r
  ComRPmat <- round(Rmat$P,3)
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
          ifelse(ComRPmat[i,g] <= .01, ComRmatFin[i,g] <- paste(ComRmat[i,g], '**',sep=""),
                 ifelse(ComRPmat[i,g] <= .05, ComRmatFin[i,g] <- paste(ComRmat[i,g], '*',sep=""),
                        ifelse(ComRPmat[i,g] <= .1, ComRmatFin[i,g] <- paste(ComRmat[i,g], "'",sep=""),
                               ifelse(i==g,ComRmatFin[i,g]<-paste('n=',ComRnmat[i,g],sep=""),
                                      ComRmatFin[i,g] <- ComRmat[i,g]))))
        }
      }
    }
    if(strict==TRUE){
      for(i in 1:nrow(ComRmat)) {
        for (g in 1:ncol(ComRmat)){
          ifelse(ComRPmat[i,g] <= .01, ComRmatFin[i,g] <- paste(ComRmat[i,g], '*',sep=""),
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
                   round = 2,...) #how many decimal places to round to when printing
  UseMethod("r_table")

r_table.default<-function(x, #data frame that can be coerced to a matrix
                          vars=NULL, #select the variables to be used
                          with=NULL, #the variables that will appear across the top
                          flag=TRUE, #flag significance
                          strict=FALSE, #If TRUE, it will use .01 cutoff for significance
                          round = 2,...) #how many decimal places to round to
{
  r_table<-cor_table(x,vars=vars,with=with,flag=flag,strict=strict,round=round)
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

#### Potentially valuable coding scheme ######

##### attempt at programmatically comparing nomological nets across factors #####
#temp <- list()
#multiple.cor.comp <- function(outcome, dataframe, factors){
#  for(i in 1:(factors-1)){
#    for(g in 2:factors){
#      if(i != g){
#        if(i < g){
#          result <- as.formula(paste("~F", factors, '.', i, '+',outcome,"|F", factors, '.', g, '+',outcome, sep = '')) %>%
#            cocor(dataframe, alternative = 't', test = 'meng1992', alpha = .01)
#          temp[[i]] <- result@meng1992$p.value
#          #print(paste(i, 'v', g)) #was able to get this function to print the comparisons, but that is as close as I could get
#        }
#      }
#    }
#  }
#  return(temp)
#}

#example attempt
#multiple.cor.comp('n',outcomesdf, 3)
#only provided two of the three comparisons because of limited temp[[i]] assignment

#possible that the the function could be simplified with following:
#cor.comp <- function(outcome, dataframe, factors, var1, var2){
#  y=as.name(outcome)
#  comparison <- as.formula(paste('~F',factors,'.', var1, '+', y, '|F', factors, '.', var2, '+', y, sep = '')) %>%
#    cocor(dataframe, alternative = 't',test = 'meng1992',alpha = .01)
#  return(comparison@meng1992$p.value)
#}


#also don't know if attempt would work in a supply loop as it was intented to be used - see other example below:

#f3comps<-data.frame(
#  'AEvN' = sapply(names(outcomesdf[1:6]),function(x){
#    y=as.name(x)
#    round(cocor(as.formula(paste("~F3.1 +",y,"| F3.2 +",y)),outcomesdf,alternative = 't',test = 'meng1992',alpha = .01)@meng1992$p.value,4)
#  }),
#
#  'AEvA' = sapply(names(outcomesdf[1:6]),function(x){
#    y=as.name(x)
#    round(cocor(as.formula(paste("~F3.1 +",y,"| F3.3 +",y)),outcomesdf,alternative = 't',test = 'meng1992',alpha = .01)@meng1992$p.value,4)
#  }),
#
#  'AvN' = sapply(names(outcomesdf[1:6]),function(x){
#    y=as.name(x)
#    round(cocor(as.formula(paste("~F3.2 +",y,"| F3.3 +",y)),outcomesdf,alternative = 't',test = 'meng1992',alpha = .01)@meng1992$p.value,4)
#  })
#)

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
make.dat_cors <- function(full.df, cut, cor){
  require(magrittr)
  require(qgraph)
    if(cor == 'cor'){
    dat_cors <- cor(full.df, use = 'pairwise')
  }
  if(cor == 'tet'){
    dat_cors <- tetrachoric(clinical.df)$rho
  }
  if(cor == 'auto'){
    dat_cors <- cor_auto(full.df, missing = 'pairwise')
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
make.cors_df <- function(full.df, cut, item.content = NULL, cor = 'cor'){
  dat_cors <- make.dat_cors(full.df, cut, cor)
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

extract.structures <- function(fa.list){
  structure.list <- lapply(fa.list, function(x){unclass(x$Structure)})
  order.list <- lapply(structure.list, function(x){
    data.frame(x) %>% names() %>% gsub('PA','',., fixed = T) %>% as.numeric() %>% order()
  })
  structure.list <- lapply(1:length(structure.list), function(x){
    structure.list[[x]] <- structure.list[[x]][,c(order.list[[x]])]
  })
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

extract.loadings <- function(fa.list){
  loadings.list <- lapply(fa.list, function(x){unclass(x$loadings)})
  order.list <- lapply(loadings.list, function(x){
    data.frame(x) %>% names() %>% gsub('PA','',., fixed = T) %>% as.numeric() %>% order()
  })
  loadings.list <- lapply(1:length(loadings.list), function(x){
    loadings.list[[x]] <- loadings.list[[x]][,c(order.list[[x]])]
  })
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
                    
extract.vaccounted <- function(fa.list){
  vaccounted.list <- lapply(fa.list, function(x){unclass(x$Vaccounted)})
  vaccounted.list[[1]] <- rbind(vaccounted.list[[1]], 'Cumulative Var' = vaccounted.list[[1]][2,], 'Proportion Explained' = 1, 'Cumulative Proportion' = 1)
  order.list <- lapply(vaccounted.list, function(x){
    data.frame(x) %>% names() %>% gsub('PA','',., fixed = T) %>% as.numeric() %>% order()
  })
  vaccounted.list <- lapply(1:length(vaccounted.list), function(x){
    vaccounted.list[[x]] <- vaccounted.list[[x]][,c(order.list[[x]])]
  })
  
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

extract.scores <- function(fa.list){
  scores.list <- lapply(fa.list, function(x){x$scores})
  order.list <- lapply(scores.list, function(x){
    data.frame(x) %>% names() %>% gsub('PA','',., fixed = T) %>% as.numeric() %>% order()
  })
  scores.list <- lapply(1:length(scores.list), function(x){
    scores.list[[x]] <- scores.list[[x]][,c(order.list[[x]])]
  })
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

avg.list <- function(list){
  Reduce("+", list)/length(list)
}

#Can pull specific parts of EFA output using sapply
#sapply(fa.list, function(x){x$RMSEA})

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
# don't need to use this if packages are stored in a user library
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
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
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

<<<<<<< HEAD
#### Attempts at generating figure from Bass-ackward results  #####
# --- NOTE FUNCTION IS NOT FINISHED AND WON'T WORK
=======
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

# Foundation for counting the number of edges each node has with a magnitude greater than some cutoff ####
#pcl.btw.net.no.t$graph %>% apply(.,c(1,2),function(x){ifelse(abs(x) > .05, return(1), return(0))}) %>% apply(.,2,sum)

#### Attempts at generating figure from Bass-ackward results #####
>>>>>>> 3b254b8a00f46ecfd562240209b732261d4dc993
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

