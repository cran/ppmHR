#' Privacy-protecting estimation of site-specific hazard ratios using inverse probability weighted Cox models
#'
#' This function allows privacy-protecting estimation of the site-specific hazard ratios using inverse probability weighted Cox models, under both the conventional inverse probability weights and the stabilized weights. The robust sandwich variance estimation method is used for estimating the variance of the log hazard ratio estimates. The Breslow method is used to handle tied events. This is a function for the analysis center.
#'@param tableList A list of summary-level riskset tables shared by data-contributing sites. Each data-contributing site provides a list of two riskset tables; the first using the conventional inverse probability weights and the second using the stabilized weights. Data-contributing sites can obtain their summary-level riskset tables using the function \code{\link[ppmHR]{createRisksetTable}} in this package.
#'@param initialHR An initial value for hazard ratio when solving the site-specific weighted partial likelihood score equation. The default is 1.
#'@param endpoint A value of the end of follow-up used to conduct sensitivity analysis. Observed events in the original data that occur after this value will be censored. The default is Inf, which means that we use the original data without conducting sensitivity analysis. For riskset tables that do not provide event times, endpoint should be left as the default Inf.
#'@param confidence A confidence level between 0 and 1. The default is 0.95 corresponding to a 95 per cent confidence interval.
#'@return A matrix of inference results from inverse probability weighted Cox models for each data-contributing site. Each site has two rows of results, where the first and the second rows report the log hazard ratio estimate and associated robust sandwich standard error, hazard ratio estimate and associated normality-based confidence interval, under the conventional inverse probability weights and the stabilized weights, respectively.
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
#'#do not share event times
#'rsTb1=createRisksetTable(data=site1,indA="A",indX=c("X1","X2","X3","X4","X5"),
#'                         indStatus="status",indTime="time",truncP=1,shareEventTime="no")
#'rsTb2=createRisksetTable(data=site2,indA="A",indX=c("X1","X2","X3","X4","X5"),
#'                         indStatus="status",indTime="time",truncP=1,shareEventTime="no")
#'rsTb3=createRisksetTable(data=site3,indA="A",indX=c("X1","X2","X3","X4","X5"),
#'                         indStatus="status",indTime="time",truncP=1,shareEventTime="no")
#'#analysis center estimates site-specific hazard ratios in IPW Cox models
#'#for all data-contributing sites
#'#using summary-level riskset tables rsTb1-3 shared by data-contributing sites
#'estimateSiteHRs(list(rsTb1,rsTb2,rsTb3),initialHR=1,endpoint=Inf,confidence=0.95)
#'
#'@export
estimateSiteHRs<-function(tableList,initialHR=1,endpoint=Inf,confidence=0.95){

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


    risksetFullShare=tableList[[k]]$ipwTable
    risksetFullsShare=tableList[[k]]$stabTable

    orderDP=order(risksetFullShare$sumSquareE,risksetFullShare$sumSquareUnE,decreasing=TRUE)
    tableList[[k]]$ipwTable=risksetFullShare[orderDP,]

    orderDP=order(risksetFullsShare$sumSquareE,risksetFullsShare$sumSquareUnE,decreasing=TRUE)
    tableList[[k]]$stabTable=risksetFullsShare[orderDP,]

  }

  HRest=rep(0, nDP)
  HRests=rep(0, nDP)
  robuStdErr=rep(0, nDP)
  robuStdErrS=rep(0, nDP)

  output=NULL
  siteTag=NULL

  for (k in 1:nDP){

    coxScore<-function(HR){

      risksetFull=tableList[[k]]$ipwTable
      y=sum(risksetFull$sumEC-risksetFull$sumC*(risksetFull$sumE*HR)/
              (risksetFull$sumE*HR+risksetFull$sumUnE))


      y
    }

    HRest[k]=nleqslv(initialHR, coxScore)$x

    coxScores<-function(HR){

      risksetFull=tableList[[k]]$stabTable
      y=sum(risksetFull$sumEC-risksetFull$sumC*(risksetFull$sumE*HR)/
              (risksetFull$sumE*HR+risksetFull$sumUnE))


      y
    }

    HRests[k]=nleqslv(initialHR, coxScores)$x


    risksetFull=tableList[[k]]$ipwTable
    S0DP=risksetFull$sumE*HRest[k]+risksetFull$sumUnE
    S1DP=risksetFull$sumE*HRest[k]
    AA=sum(risksetFull$sumC*(S1DP/S0DP-(S1DP/S0DP)^2))

    q1DP=sum((1-S1DP/S0DP)^2*risksetFull$sumSquareEC+(S1DP/S0DP)^2*risksetFull$sumSquareUnEC)


    sumSquEdiffDP=risksetFull$sumSquareE-c(risksetFull$sumSquareE[-1],0)
    sumSquUnEdiffDP=risksetFull$sumSquareUnE-c(risksetFull$sumSquareUnE[-1],0)
    cumsum1EDP=cumsum(risksetFull$sumC/S0DP)
    cumsum2EDP=cumsum(risksetFull$sumC*S1DP/(S0DP^2))

    q2DP=(exp(2*log(HRest[k])))*sum(sumSquEdiffDP*cumsum1EDP^2)
    q3DP=(exp(2*log(HRest[k])))*sum(sumSquEdiffDP*cumsum2EDP^2)+sum(sumSquUnEdiffDP*cumsum2EDP^2)
    q4DP=HRest[k]*sum(risksetFull$sumSquareEC*(1-S1DP/S0DP)*cumsum1EDP)
    q5DP=HRest[k]*sum(risksetFull$sumSquareEC*(1-S1DP/S0DP)*cumsum2EDP)+sum(risksetFull$sumSquareUnEC*(0-S1DP/S0DP)*cumsum2EDP)
    q6DP=(exp(2*log(HRest[k])))*sum(sumSquEdiffDP*cumsum1EDP*cumsum2EDP)

    q=q1DP+q2DP+q3DP-2*q4DP+2*q5DP-2*q6DP

    robuVar=q/(AA^2)

    robuStdErr[k]=robuVar^0.5


    risksetFulls=tableList[[k]]$stabTable

    S0DP=risksetFulls$sumE*HRests[k]+risksetFulls$sumUnE
    S1DP=risksetFulls$sumE*HRests[k]
    AA=sum(risksetFulls$sumC*(S1DP/S0DP-(S1DP/S0DP)^2))

    q1DP=sum((1-S1DP/S0DP)^2*risksetFulls$sumSquareEC+(S1DP/S0DP)^2*risksetFulls$sumSquareUnEC)

    sumSquEdiffDP=risksetFulls$sumSquareE-c(risksetFulls$sumSquareE[-1],0)
    sumSquUnEdiffDP=risksetFulls$sumSquareUnE-c(risksetFulls$sumSquareUnE[-1],0)
    cumsum1EDP=cumsum(risksetFulls$sumC/S0DP)
    cumsum2EDP=cumsum(risksetFulls$sumC*S1DP/(S0DP^2))

    q2DP=(exp(2*log(HRests[k])))*sum(sumSquEdiffDP*cumsum1EDP^2)
    q3DP=(exp(2*log(HRests[k])))*sum(sumSquEdiffDP*cumsum2EDP^2)+sum(sumSquUnEdiffDP*cumsum2EDP^2)
    q4DP=HRests[k]*sum(risksetFulls$sumSquareEC*(1-S1DP/S0DP)*cumsum1EDP)
    q5DP=HRests[k]*sum(risksetFulls$sumSquareEC*(1-S1DP/S0DP)*cumsum2EDP)+sum(risksetFulls$sumSquareUnEC*(0-S1DP/S0DP)*cumsum2EDP)
    q6DP=(exp(2*log(HRests[k])))*sum(sumSquEdiffDP*cumsum1EDP*cumsum2EDP)

    q=q1DP+q2DP+q3DP-2*q4DP+2*q5DP-2*q6DP

    robuVar=q/(AA^2)
    robuStdErrS[k]=robuVar^0.5


    lowProp=log(HRest[k])-qnorm(1-(1-confidence)/2)*robuStdErr[k]
    upProp=log(HRest[k])+qnorm(1-(1-confidence)/2)*robuStdErr[k]

    lowPropS=log(HRests[k])-qnorm(1-(1-confidence)/2)*robuStdErrS[k]
    upPropS=log(HRests[k])+qnorm(1-(1-confidence)/2)*robuStdErrS[k]


    est=c(log(HRest[k]),log(HRests[k]))
    hrest=exp(est)
    se=c(robuStdErr[k],robuStdErrS[k])
    low=exp(c(lowProp,lowPropS))
    up=exp(c(upProp,upPropS))
    outputAdd=cbind(est,se,hrest,low,up)

    rownames(outputAdd)=c(paste("Site ",k,"-IPW",sep = ""),paste("Site ",k,"-stabilized",sep = ""))

    output=rbind(output,outputAdd)

    colnames(output)=c("log HR Estimate","Standard Error","HR Estimate",
                       paste("HR ", confidence*100,"% CI", "-low", sep =""),
                       paste("HR ", confidence*100,"% CI", "-up", sep =""))

  }

  output

}


