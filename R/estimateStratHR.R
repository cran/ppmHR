#' Privacy-protecting estimation of overall hazard ratios using inverse probability weighted Cox models stratified on data-contributing sites
#'
#' This function estimates the overall hazard ratios using inverse probability weighted Cox models stratified on data-contributing sites, under both the conventional inverse probability weights and the stabilized weights. The robust sandwich variance estimation method is used for estimating the variance of the log hazard ratio estimates. The Breslow method is used to handle tied events. This is a function for the analysis center.
#'@param tableList A list of summary-level riskset tables shared by data-contributing sites. Each data-contributing site provides a list of two riskset tables; the first using the conventional inverse probability weights and the second using the stabilized weights. Data-contributing sites can obtain their summary-level riskset tables using the function \code{\link[ppmHR]{createRisksetTable}} in this package.
#'@param initialHR An initial value for hazard ratio when solving the stratified weighted partial likelihood score equation. The default is 1.
#'@param endpoint A value of the end of follow-up used to conduct sensitivity analysis. Observed events in the original data that occur after this value will be censored. The default is Inf, which means that we use the original data without conducting sensitivity analysis. For riskset tables that do not provide event times, endpoint should be left as the default Inf.
#'@param confidence A confidence level between 0 and 1. The default is 0.95 corresponding to a 95 per cent confidence interval.
#'@return A matrix of inference results from inverse probability weighted Cox models stratified on data-contributing sites. The first and the second rows report the log hazard ratio estimate and associated robust sandwich standard error, hazard ratio estimate and associated normality-based confidence interval, under the conventional inverse probability weights and the stabilized weights, respectively.
#'@import nleqslv
#'@import stats
#'
#'@examples
#'#load example datasets in the package
#'#site1-3 contain individual-level data of data-contributing sites 1-3
#'data(site1)
#'data(site2)
#'data(site3)
#'#data-contributing sites 1-3 create summary-level riskset tables
#'#with logistic propensity score model A~X1+X2+X3+X4+X5
#'#no weight truncation
#'#agree to share event times
#'rsTb1=createRisksetTable(data=site1,indA="A",indX=c("X1","X2","X3","X4","X5"),
#'                         indStatus="status",indTime="time",truncP=1,shareEventTime="yes")
#'rsTb2=createRisksetTable(data=site2,indA="A",indX=c("X1","X2","X3","X4","X5"),
#'                         indStatus="status",indTime="time",truncP=1,shareEventTime="yes")
#'rsTb3=createRisksetTable(data=site3,indA="A",indX=c("X1","X2","X3","X4","X5"),
#'                         indStatus="status",indTime="time",truncP=1,shareEventTime="yes")
#'#analysis center estimates hazard ratio in a stratified IPW Cox model
#'#using summary-level riskset tables rsTb1-3 shared by data-contributing sites
#'estimateStratHR(list(rsTb1,rsTb2,rsTb3),initialHR=1,endpoint=Inf,confidence=0.95)
#'#sensitivity analysis at endpoint 20
#'estimateStratHR(list(rsTb1,rsTb2,rsTb3),initialHR=1,endpoint=20,confidence=0.95)
#'
#'@export
estimateStratHR<-function(tableList,initialHR=1,endpoint=Inf,confidence=0.95){

  nDP=length(tableList)

  for (k in 1:nDP){
    if(dim(tableList[[k]]$ipwTable)[2]==9){
      timeK=tableList[[k]]$ipwTable$eventTime
      timeInd=which(timeK<=endpoint)
      tableList[[k]]$ipwTable=tableList[[k]]$ipwTable[timeInd,]
      tableList[[k]]$stabTable= tableList[[k]]$stabTable[timeInd,]
    }else if(dim(tableList[[k]]$ipwTable)[2]==8&endpoint!=Inf){
      stop('endpoint should be left as default when some riskset tables do not share event times')
    }
    else{
    }
  }


  for (k in 1:nDP){
    risksetFullShare=tableList[[k]]$ipwTable
    risksetFullsShare=tableList[[k]]$stabTable

    orderDP=order(risksetFullShare$sumSquareE,risksetFullShare$sumSquareUnE,decreasing=TRUE)
    tableList[[k]]$ipwTable=risksetFullShare[orderDP,]

    orderDP=order(risksetFullsShare$sumSquareE,risksetFullsShare$sumSquareUnE,decreasing=TRUE)
    tableList[[k]]$stabTable=risksetFullsShare[orderDP,]

  }



  coxScore<-function(HR){

    y=rep(0, nDP)

    for (k in 1:nDP){
      risksetFull=tableList[[k]]$ipwTable
      y[k]=sum(risksetFull$sumEC-risksetFull$sumC*(risksetFull$sumE*HR)/
                 (risksetFull$sumE*HR+risksetFull$sumUnE))

    }

    sum(y)
  }

  HRest=nleqslv(initialHR, coxScore)$x


  coxScores<-function(HR){
    y=rep(0, nDP)

    for (k in 1:nDP){
      risksetFull=tableList[[k]]$stabTable
      y[k]=sum(risksetFull$sumEC-risksetFull$sumC*(risksetFull$sumE*HR)/
                 (risksetFull$sumE*HR+risksetFull$sumUnE))

    }

    sum(y)
  }

  HRests=nleqslv(initialHR, coxScores)$x


  AAvec=rep(0, nDP)
  q1DPvec=rep(0, nDP)
  q2DPvec=rep(0, nDP)
  q3DPvec=rep(0, nDP)
  q4DPvec=rep(0, nDP)
  q5DPvec=rep(0, nDP)
  q6DPvec=rep(0, nDP)

  for (k in 1:nDP){

    risksetFull=tableList[[k]]$ipwTable
    S0DP=risksetFull$sumE*HRest+risksetFull$sumUnE
    S1DP=risksetFull$sumE*HRest
    AAvec[k]=sum(risksetFull$sumC*(S1DP/S0DP-(S1DP/S0DP)^2))

    q1DPvec[k]=sum((1-S1DP/S0DP)^2*risksetFull$sumSquareEC+(S1DP/S0DP)^2*risksetFull$sumSquareUnEC)


    sumSquEdiffDP=risksetFull$sumSquareE-c(risksetFull$sumSquareE[-1],0)
    sumSquUnEdiffDP=risksetFull$sumSquareUnE-c(risksetFull$sumSquareUnE[-1],0)
    cumsum1EDP=cumsum(risksetFull$sumC/S0DP)
    cumsum2EDP=cumsum(risksetFull$sumC*S1DP/(S0DP^2))

    q2DPvec[k]=(exp(2*log(HRest)))*sum(sumSquEdiffDP*cumsum1EDP^2)
    q3DPvec[k]=(exp(2*log(HRest)))*sum(sumSquEdiffDP*cumsum2EDP^2)+sum(sumSquUnEdiffDP*cumsum2EDP^2)
    q4DPvec[k]=HRest*sum(risksetFull$sumSquareEC*(1-S1DP/S0DP)*cumsum1EDP)
    q5DPvec[k]=HRest*sum(risksetFull$sumSquareEC*(1-S1DP/S0DP)*cumsum2EDP)+sum(risksetFull$sumSquareUnEC*(0-S1DP/S0DP)*cumsum2EDP)
    q6DPvec[k]=(exp(2*log(HRest)))*sum(sumSquEdiffDP*cumsum1EDP*cumsum2EDP)

  }

  AA=sum(AAvec)

  q=sum(q1DPvec)+sum(q2DPvec)+sum(q3DPvec)-2*sum(q4DPvec)+2*sum(q5DPvec)-2*sum(q6DPvec)


  robuVar=q/(AA^2)

  robuStdErr=robuVar^0.5


  AAvec=rep(0, nDP)
  q1DPvec=rep(0, nDP)
  q2DPvec=rep(0, nDP)
  q3DPvec=rep(0, nDP)
  q4DPvec=rep(0, nDP)
  q5DPvec=rep(0, nDP)
  q6DPvec=rep(0, nDP)

  for (k in 1:nDP){
    risksetFulls=tableList[[k]]$stabTable

    S0DP=risksetFulls$sumE*HRests+risksetFulls$sumUnE
    S1DP=risksetFulls$sumE*HRests
    AAvec[k]=sum(risksetFulls$sumC*(S1DP/S0DP-(S1DP/S0DP)^2))

    q1DPvec[k]=sum((1-S1DP/S0DP)^2*risksetFulls$sumSquareEC+(S1DP/S0DP)^2*risksetFulls$sumSquareUnEC)

    sumSquEdiffDP=risksetFulls$sumSquareE-c(risksetFulls$sumSquareE[-1],0)
    sumSquUnEdiffDP=risksetFulls$sumSquareUnE-c(risksetFulls$sumSquareUnE[-1],0)
    cumsum1EDP=cumsum(risksetFulls$sumC/S0DP)
    cumsum2EDP=cumsum(risksetFulls$sumC*S1DP/(S0DP^2))

    q2DPvec[k]=(exp(2*log(HRests)))*sum(sumSquEdiffDP*cumsum1EDP^2)
    q3DPvec[k]=(exp(2*log(HRests)))*sum(sumSquEdiffDP*cumsum2EDP^2)+sum(sumSquUnEdiffDP*cumsum2EDP^2)
    q4DPvec[k]=HRests*sum(risksetFulls$sumSquareEC*(1-S1DP/S0DP)*cumsum1EDP)
    q5DPvec[k]=HRests*sum(risksetFulls$sumSquareEC*(1-S1DP/S0DP)*cumsum2EDP)+sum(risksetFulls$sumSquareUnEC*(0-S1DP/S0DP)*cumsum2EDP)
    q6DPvec[k]=(exp(2*log(HRests)))*sum(sumSquEdiffDP*cumsum1EDP*cumsum2EDP)
  }

  AA=sum(AAvec)
  q=sum(q1DPvec)+sum(q2DPvec)+sum(q3DPvec)-2*sum(q4DPvec)+2*sum(q5DPvec)-2*sum(q6DPvec)

  robuVar=q/(AA^2)
  robuStdErrS=robuVar^0.5


  lowProp=log(HRest)-qnorm(1-(1-confidence)/2)*robuStdErr
  upProp=log(HRest)+qnorm(1-(1-confidence)/2)*robuStdErr

  lowPropS=log(HRests)-qnorm(1-(1-confidence)/2)*robuStdErrS
  upPropS=log(HRests)+qnorm(1-(1-confidence)/2)*robuStdErrS

  mtd=c("IPW","stabilized")
  est=c(log(HRest),log(HRests))
  hrest=exp(est)
  se=c(robuStdErr,robuStdErrS)
  low=exp(c(lowProp,lowPropS))
  up=exp(c(upProp,upPropS))
  output=cbind(est,se,hrest,low,up)
  colnames(output)=c("log HR Estimate","Standard Error","HR Estimate",
                     paste("HR ", confidence*100,"% CI", "-low", sep =""),
                     paste("HR ", confidence*100,"% CI", "-up", sep =""))
  rownames(output)=mtd
  output

}




