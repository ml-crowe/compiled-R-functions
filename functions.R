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
varlist <- function (df=NULL,type=c("numeric","factor","character"), pattern=NULL, exclude=NULL, ignore.case=TRUE) {
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
#quick function that I load into all of my projects were identifying the number of NAs in a row
numNAs<-function(x){
  sum(is.na(x))
}

#### IRT function ####
#Not regularly used

#irt.func<-function(dataframe, irt.model = NULL, factors = 1, seed = 123, ...){ #dataframe should not include ID #variable
#  if(apply(dataframe, 1, numNAs) %>% equals(length(dataframe[1,])) %>% any){
#    print("Warning, some rows removed for being empty")
#    missing.rows.list <- remove.missing.rows(dataframe)
#  }
#  dataframe <- data.frame(missing.rows.list$df)
#  irt<-mirt(dataframe,factors,itemtype = irt.model, technical = list(removeEmptyRows = TRUE)) #using irt
#  scores<-fscores(irt,method='EAP', full.scores=TRUE,scores.only=TRUE) #EAP estimation method for the scores
#  saved.scores <- fscores(irt, method = 'EAP', full.scores = TRUE, full.scores.SE = TRUE)
#  saved.scores <- data.frame(saved.scores)
#  if(any(is.na(data.frame(dataframe)))){
#    set.seed(seed)
#    fulldataframe<-imputeMissing(irt,scores) #just imputing  the data one time
#    firt<-mirt(fulldataframe,1,itemtype = irt.model, technical = list(removeEmptyRows = TRUE)) #save imputed dataset
#    m2<-M2(firt)#save M2 statistics
#    coefs<-coef(firt, simplify=TRUE) #save item parameters
#    items<-itemfit(firt,simplify=TRUE) #save item fit statistics
#    output<-list("model" = firt, "m2" = m2, "scores" = saved.scores, 'coefs' = coefs, 'itemfit' = items)
#  } else{
#    m2<-M2(irt)#save M2 statistics
#    coefs<-coef(irt, simplify=TRUE) #save item parameters
#    items<-itemfit(irt,simplify=TRUE) #save item fit statistics
#    output<-list("model" = irt, "m2" = m2, "scores" = saved.scores, 'coefs' = coefs, 'itemfit' = items)
#  }
#  return(output)
#}

##### Copy dataframe to Excel through clipboard #######
write.excel <- function(x,row.names=FALSE,col.names=TRUE,...) {
  write.table(x,"clipboard-10000",sep="\t",row.names=row.names,col.names=col.names,...)
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

###test of skew, I don't think this is right, see Crawley 2012 R book in Zotero and http://www.real-statistics.com/tests-normality-and-symmetry/analysis-skewness-kurtosis/ ####
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


#also don't know if attempt would work in a supply loop was it was intented to be used - see other example below:

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
