#' Creation of summary-level riskset tables from individual-level data
#'
#' This function creates summary-level riskset tables from individual-level data at each data-contributing site, under both the conventional inverse probability weights and the stabilized weights. A logistic regression model is used to estimate the propensity scores with an option to do weight truncation. This is a function for data-contributing sites.
#'@param data An individual-level dataset in the form of R data frame.
#'@param indA A column name indicating the exposure/treatment variable.
#'@param indX A vector of column names indicating the baseline covariates to be included in the propensity score model.
#'@param indStatus A column name indicating the non-censoring status (1 if observed and 0 if censored).
#'@param indTime A column name indicating the outcome variable, i.e., min(true event time, censoring time).
#'@param truncP A value between 0 and 1 indicating the percentile used for weight truncation. The default is 1, corresponding to no weight truncation.
#'@param shareEventTime If the data-contributing site would like to share a column of event times, then set \code{shareEventTime="yes"}. If not, then set \code{shareEventTime="no"}. By default,  \code{shareEventTime="no"}.
#'@return A list of two summary-level riskset tables to be shared with the analysis center, under both the conventional inverse probability weights (i.e., $ipwTable) and the stabilized weights (i.e., $stabTable). If \code{shareEventTime="no"}, then a riskset table has 8 columns: total weights of exposed cases or events, total weights of cases or events, total weights of exposed individuals, total weights of unexposed individuals, total squared weights of exposed cases or events, total squared weights of unexposed cases or events, total squared weights of exposed individuals, total squared weights of unexposed individuals. If \code{shareEventTime="yes"}, then a riskset table further has a column of event times.
#'@import stats
#'
#'@examples
#'#load an example dataset site1.RData in the package
#'#site1 contains individual-level data of data-contributing site 1
#'data(site1)
#'#data-contributing site 1 creates its two summary-level riskset tables
#'#with logistic propensity score model A~X1+X2+X3+X4+X5
#'#no weight truncation
#'#agree to share a column of event times
#'rsTb1=createRisksetTable(data=site1,indA="A",indX=c("X1","X2","X3","X4","X5"),
#'                         indStatus="status",indTime="time",truncP=1,shareEventTime="yes")
#'#print the first six rows of riskset table using conventional weights
#'head(rsTb1$ipwTable)
#'#print the first six rows of riskset table using stabilized weights
#'head(rsTb1$stabTable)
#'
#'@export
createRisksetTable<-function(data,indA,indX,indStatus,indTime,truncP=1,shareEventTime="no"){

  dat=data

  dat$A=dat[,indA]
  dat$time=dat[,indTime]
  dat$status=dat[,indStatus]

  nX=length(indX)+1
  covX0=dat[,indX]
  A=dat$A
  psmd=glm(A~.,family="binomial",data=as.data.frame(cbind(A,covX0)))
  psfit=predict(psmd, type = "response")

  dat$wt=dat$A/psfit+(1-dat$A)/(1-psfit)
  prevA=mean(dat$A)
  dat$swt=prevA*dat$A/psfit+(1-prevA)*(1-dat$A)/(1-psfit)

  wtQ=quantile(dat$wt,truncP)
  dat$wt[which(dat$wt>=wtQ)]=wtQ

  swtQ=quantile(dat$swt,truncP)
  dat$swt[which(dat$swt>=swtQ)]=swtQ


  eventPa=subset(dat,dat$status==1)
  eventimes=unique(eventPa$time)

  RScol3=rep(0,length(eventimes))
  RScol4=RScol3
  RScol1=RScol3
  RScol2=RScol3

  RScol5=RScol3
  RScol6=RScol3
  RScol7=RScol3
  RScol8=RScol3

  for (i in 1:length(eventimes)){
    idc=which(dat$time==eventimes[i]&dat$status==1)
    RD1=dat[idc,]

    RScol1[i]=sum(RD1$wt*RD1$A)

    RScol2[i]=sum(RD1$wt)

    RScol5[i]=sum((RD1$wt)^2*RD1$A)

    RScol6[i]=sum((RD1$wt)^2*(1-RD1$A))

    RSind=which(dat$time>=eventimes[i])
    RS1=dat[RSind,]

    RScol3[i]=sum(RS1$wt*RS1$A)
    RScol4[i]=sum(RS1$wt*(1-RS1$A))

    RScol7[i]=sum((RS1$wt)^2*RS1$A)
    RScol8[i]=sum((RS1$wt)^2*(1-RS1$A))
  }

  if(shareEventTime=="no"){
    risksetFull1Share=data.frame(sumEC=RScol1,sumC=RScol2,sumE=RScol3,sumUnE=RScol4,
                                 sumSquareEC=RScol5,sumSquareUnEC=RScol6,sumSquareE=RScol7,sumSquareUnE=RScol8)

  } else if(shareEventTime=="yes"){
    risksetFull1Share=data.frame(sumEC=RScol1,sumC=RScol2,sumE=RScol3,sumUnE=RScol4,
                                 sumSquareEC=RScol5,sumSquareUnEC=RScol6,sumSquareE=RScol7,sumSquareUnE=RScol8,
                                 eventTime=eventimes)
  } else {
    stop('shareEventTime must be "yes" or "no"')
  }



  eventPa=subset(dat,dat$status==1)
  eventimes=unique(eventPa$time)


  RScol3=rep(0,length(eventimes))
  RScol4=RScol3
  RScol1=RScol3
  RScol2=RScol3

  RScol5=RScol3
  RScol6=RScol3
  RScol7=RScol3
  RScol8=RScol3



  for (i in 1:length(eventimes)){
    idc=which(dat$time==eventimes[i]&dat$status==1)
    RD1=dat[idc,]

    RScol1[i]=sum(RD1$swt*RD1$A)

    RScol2[i]=sum(RD1$swt)

    RScol5[i]=sum((RD1$swt)^2*RD1$A)

    RScol6[i]=sum((RD1$swt)^2*(1-RD1$A))

    RSind=which(dat$time>=eventimes[i])
    RS1=dat[RSind,]

    RScol3[i]=sum(RS1$swt*RS1$A)
    RScol4[i]=sum(RS1$swt*(1-RS1$A))

    RScol7[i]=sum((RS1$swt)^2*RS1$A)
    RScol8[i]=sum((RS1$swt)^2*(1-RS1$A))
  }


  if(shareEventTime=="no"){
    risksetFulls1Share=data.frame(sumEC=RScol1,sumC=RScol2,sumE=RScol3,sumUnE=RScol4,
                                  sumSquareEC=RScol5,sumSquareUnEC=RScol6,sumSquareE=RScol7,sumSquareUnE=RScol8)

  } else{
    risksetFulls1Share=data.frame(sumEC=RScol1,sumC=RScol2,sumE=RScol3,sumUnE=RScol4,
                                  sumSquareEC=RScol5,sumSquareUnEC=RScol6,sumSquareE=RScol7,sumSquareUnE=RScol8,
                                  eventTime=eventimes)
  }


  res=list(ipwTable=risksetFull1Share,stabTable=risksetFulls1Share)
  res

}
