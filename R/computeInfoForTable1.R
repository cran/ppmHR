#' Computation of summary-level information needed for the creation of a "Table 1" of baseline patient characteristics
#'
#' Computation of summary-level information needed for the creation of a "Table 1" that shows the baseline patient characteristics in the original unweighted and the inverse probability weighted samples. A logistic regression model is used to estimate propensity scores with an option to do weight truncation. This is a function for data-contributing sites.
#'@param data An individual-level dataset in the form of R data frame.
#'@param indA A column name indicating the exposure/treatment variable.
#'@param indX A vector of column names indicating the baseline covariates to be included in the propensity score model.
#'@param truncP A value between 0 and 1 indicating the percentile used for weight truncation. The default is 1, corresponding to no weight truncation.
#'@return A data frame of summary-level information needed for the creation of a "Table 1" of baseline patient characteristics. Each covariate has a row of 19 numbers: covariate type ("yes" if it is binary, "no" if it is continuous/count), number of exposed individuals (same across rows), covariate mean of exposed individuals, covariate mean square of exposed individuals, total conventional weights of exposed individuals, weighted covariate mean of exposed individuals with conventional weights, weighted covariate mean square of exposed individuals with conventional weights, total stabilized weights of exposed individuals, weighted covariate mean of exposed individuals with stabilized weights, weighted covariate mean square of exposed individuals with stabilized weights, number of unexposed individuals (same across rows), covariate mean of unexposed individuals, covariate mean square of unexposed individuals, total conventional weights of unexposed individuals, weighted covariate mean of unexposed individuals with conventional weights, weighted covariate mean square of unexposed individuals with conventional weights, total stabilized weights of unexposed individuals, weighted covariate mean of unexposed individuals with stabilized weights, and weighted covariate mean square of unexposed individuals with stabilized weights.
#'@import stats
#'
#'@examples
#'#load an example dataset site1.RData in the package
#'#site1 contains individual-level data of data-contributing site 1
#'data(site1)
#'#site 1 computes summary-level information needed for creating the "Table 1"
#'#with logistic propensity score model A~X1+X2+X3+X4+X5
#'#no weight truncation
#'computeInfoForTable1(data=site1,indA="A",indX=c("X1","X2","X3","X4","X5"),truncP=1)
#'#with truncation: set weights larger than the 90% quantile of original weights to the 90% quantile
#'computeInfoForTable1(data=site1,indA="A",indX=c("X1","X2","X3","X4","X5"),truncP=0.9)
#'
#'@export
computeInfoForTable1<-function(data,indA,indX,truncP=1){

  dat=data
  N=nrow(dat)
  dat$A=dat[,indA]
  trtIndex=which(dat$A==1)

  Ntrt=length(trtIndex)
  Ncon=N-Ntrt

  covX0=dat[,indX]

  isBinary<-function(vec){
    numFalse=sum(1-as.integer(vec==vec^2))
    ifelse(numFalse==0,"yes","no")
  }

  binary=as.vector(apply(covX0,2,isBinary))

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


  wtrtSum=sum(wt[trtIndex])
  swtrtSum=sum(swt[trtIndex])
  wconSum=sum(wt[-trtIndex])
  swconSum=sum(swt[-trtIndex])

  covasTrt=covX0[trtIndex,]
  MeanTrt=apply(covasTrt,2,mean)
  MeanSqTrt=apply(covasTrt^2,2,mean)

  covasCon=covX0[-trtIndex,]
  MeanCon=apply(covasCon,2,mean)
  MeanSqCon=apply(covasCon^2,2,mean)

  wtcovas=diag(wt)%*%as.matrix(covX0)
  wtcovasTrt=wtcovas[trtIndex,]
  wtMeanTrt=apply(wtcovasTrt,2,sum)/sum(wt[trtIndex])
  wtcovasCon=wtcovas[-trtIndex,]
  wtMeanCon=apply(wtcovasCon,2,sum)/sum(wt[-trtIndex])

  wtcovas=diag(wt)%*%((as.matrix(covX0))^2)
  wtcovasTrt=wtcovas[trtIndex,]
  wtMeanSqTrt=apply(wtcovasTrt,2,sum)/sum(wt[trtIndex])
  wtcovasCon=wtcovas[-trtIndex,]
  wtMeanSqCon=apply(wtcovasCon,2,sum)/sum(wt[-trtIndex])

  swtcovas=diag(swt)%*%as.matrix(covX0)
  swtcovasTrt=swtcovas[trtIndex,]
  swtMeanTrt=apply(swtcovasTrt,2,sum)/sum(swt[trtIndex])
  swtcovasCon=swtcovas[-trtIndex,]
  swtMeanCon=apply(swtcovasCon,2,sum)/sum(swt[-trtIndex])

  swtcovas=diag(swt)%*%((as.matrix(covX0))^2)
  swtcovasTrt=swtcovas[trtIndex,]
  swtMeanSqTrt=apply(swtcovasTrt,2,sum)/sum(swt[trtIndex])
  swtcovasCon=swtcovas[-trtIndex,]
  swtMeanSqCon=apply(swtcovasCon,2,sum)/sum(swt[-trtIndex])

  infoTable=data.frame(binary,Ntrt,MeanTrt,MeanSqTrt,wtrtSum,wtMeanTrt,wtMeanSqTrt,
                          swtrtSum,swtMeanTrt,swtMeanSqTrt,
                          Ncon,MeanCon,MeanSqCon,wconSum,wtMeanCon,wtMeanSqCon,
                          swconSum,swtMeanCon,swtMeanSqCon)
  rownames(infoTable)=indX
  colnames(infoTable)=c("binary","nE","XmeanE","XmeanSqE",
                        "sumWtE.ipw","XmeanE.ipw","XmeanSqE.ipw",
                        "sumWtE.stab","XmeanE.stab","XmeanSqE.stab",
                        "nUnE","XmeanUnE","XmeanSqUnE",
                        "sumWtUnE.ipw","XmeanUnE.ipw","XmeanSqUnE.ipw",
                        "sumWtUnE.stab","XmeanUnE.stab","XmeanSqUnE.stab")

  infoTable
}






