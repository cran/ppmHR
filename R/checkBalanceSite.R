#' Checking of site-specific covariate balance before and after weighting
#'
#' This function checks for site-specific covariate balance before and after weighting using the mean of baseline covariates in the original unweighted and the inverse probability weighted samples. A logistic regression model is used to estimate the propensity scores with an option to do weight truncation. This is a function for data-contributing sites.
#'@param data An individual-level dataset in the form of R data frame.
#'@param indA A column name indicating the exposure/treatment variable.
#'@param indX A vector of column names indicating the baseline covariates to be included in the propensity score model.
#'@param truncP A value between 0 and 1 indicating the percentile used for weight truncation. The default is 1, corresponding to no weight truncation.
#'@return A data frame of balance checking results. Each covariate has a row of six numbers: the mean of this covariate for exposed individuals in the original sample, for unexposed individuals in the original sample, for exposed individuals in the weighted sample with conventional weights, for unexposed individuals in the weighted sample with conventional weights, for exposed individuals in the weighted sample with stabilized weights, and for unexposed individuals in the weighted sample with stabilized weights.
#'@import stats
#'
#'@examples
#'#load an example dataset site1.RData in the package
#'#site1 contains individual-level data of data-contributing site 1
#'data(site1)
#'#data-contributing site 1 checks covariate balance before and after weighting,
#'#with logistic propensity score model A~X1+X2+X3+X4+X5
#'#no weight truncation
#'checkBalanceSite(data=site1,indA="A",indX=c("X1","X2","X3","X4","X5"),truncP=1)
#'#with truncation: set weights larger than the 90% quantile of original weights to the 90% quantile
#'checkBalanceSite(data=site1,indA="A",indX=c("X1","X2","X3","X4","X5"),truncP=0.9)
#'
#'@export
checkBalanceSite<-function(data,indA,indX,truncP=1){

  dat=data
  N=nrow(dat)
  dat$A=dat[,indA]
  trtIndex=which(dat$A==1)

  covX0=dat[,indX]
  A=dat$A
  psmd=glm(A~.,family="binomial",data=as.data.frame(cbind(A,covX0)))
  psfit=predict(psmd, type = "response")

  wt=dat$A/psfit+(1-dat$A)/(1-psfit)
  prevA=mean(dat$A)
  swt=prevA*dat$A/psfit+(1-prevA)*(1-dat$A)/(1-psfit)

  wtQ=quantile(wt,truncP)
  wt[which(wt>=wtQ)]=wtQ

  swtQ=quantile(swt,truncP)
  swt[which(swt>=swtQ)]=swtQ

  covasTrt=covX0[trtIndex,]
  MeanTrt=apply(covasTrt,2,mean)
  covasCon=covX0[-trtIndex,]
  MeanCon=apply(covasCon,2,mean)

  wtcovas=diag(wt)%*%as.matrix(covX0)
  wtcovasTrt=wtcovas[trtIndex,]
  wtMeanTrt=apply(wtcovasTrt,2,sum)/sum(wt[trtIndex])
  wtcovasCon=wtcovas[-trtIndex,]
  wtMeanCon=apply(wtcovasCon,2,sum)/sum(wt[-trtIndex])

  swtcovas=diag(swt)%*%as.matrix(covX0)
  swtcovasTrt=swtcovas[trtIndex,]
  swtMeanTrt=apply(swtcovasTrt,2,sum)/sum(swt[trtIndex])
  swtcovasCon=swtcovas[-trtIndex,]
  swtMeanCon=apply(swtcovasCon,2,sum)/sum(swt[-trtIndex])

  balanceTable=data.frame(MeanTrt,MeanCon,wtMeanTrt,wtMeanCon,swtMeanTrt,swtMeanCon)
  rownames(balanceTable)=indX
  colnames(balanceTable)=c("exposed","unexposed",
                           "exposed.ipw","unexposed.ipw",
                           "exposed.stab","unexposed.stab")

  balanceTable
}





