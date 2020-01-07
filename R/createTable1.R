#' Creation of an overall "Table 1" of baseline patient characteristics before and after weighting
#'
#' Creation of an overall "Table 1" of baseline patient characteristics in the original unweighted and the inverse probability weighted samples, using only summary-level information shared by data-contributing sites. This is a function for the analysis center.
#'@param XsummaryList A list of summary-level information tables shared by data-contributing sites. Each site provides one table, which can be obtained by using the function \code{\link[ppmHR]{computeInfoForTable1}} in this package.
#'@param digits An integer indicating the number of decimal places for the generated "Table 1". The default is 2.
#'@return A data frame of the "Table 1". Each row represents a covariate. The seven columns represent the covariate type ("yes" if it is binary, "no" if it is continuous or count), the exposed group in the original sample, the unexposed group in the original sample, the exposed group in the weighted sample with conventional weights, the unexposed group in the weighted sample with conventional weights, the exposed group in the weighted sample with stabilized weights, and the unexposed group in the weighted sample with stabilized weights. For cells of binary covariates, values are count(percentage). For cells of continuous or count covariates, values are mean(standard deviation).
#'@import stats
#'
#'@examples
#'#load example datasets in the package
#'#site1-3 contain individual-level data of data-contributing sites 1-3
#'data(site1)
#'data(site2)
#'data(site3)
#'#sites 1-3 compute summary-level information needed for creating the "Table 1"
#'#with logistic propensity score model A~X1+X2+X3+X4+X5
#'#no weight truncation
#'Xsummary1=computeInfoForTable1(data=site1,indA="A",indX=c("X1","X2","X3","X4","X5"),truncP=1)
#'Xsummary2=computeInfoForTable1(data=site2,indA="A",indX=c("X1","X2","X3","X4","X5"),truncP=1)
#'Xsummary3=computeInfoForTable1(data=site3,indA="A",indX=c("X1","X2","X3","X4","X5"),truncP=1)
#'#analysis center creates the "Table 1"
#'#using summary-level information Xsummary1-3 shared by data-contributing sites
#'#display the table with 3 decimal places
#'createTable1(list(Xsummary1,Xsummary2,Xsummary3),digits=3)
#'#analysis center can also generate site-specific "Table 1"
#'#for example, for site 1
#'createTable1(list(Xsummary1),digits=3)
#'
#'@export
createTable1<-function(XsummaryList,digits=2){

  nDP=length(XsummaryList)

  for (k in 1:nDP){
    colnames(XsummaryList[[k]])=c("binary","nE","XmeanE","XmeanSqE",
                          "sumWtE.ipw","XmeanE.ipw","XmeanSqE.ipw",
                          "sumWtE.stab","XmeanE.stab","XmeanSqE.stab",
                          "nUnE","XmeanUnE","XmeanSqUnE",
                          "sumWtUnE.ipw","XmeanUnE.ipw","XmeanSqUnE.ipw",
                          "sumWtUnE.stab","XmeanUnE.stab","XmeanSqUnE.stab")
  }

  XsummaryListCopy=XsummaryList

  for (k in 1:nDP){
    Xsummary=XsummaryList[[k]]

    Xsummary$XE=(Xsummary$XmeanE)*(Xsummary$nE)
    Xsummary$ipwXE=(Xsummary$XmeanE.ipw)*(Xsummary$sumWtE.ipw)
    Xsummary$stabXE=(Xsummary$XmeanE.stab)*(Xsummary$sumWtE.stab)

    Xsummary$X2E=(Xsummary$XmeanSqE)*(Xsummary$nE)
    Xsummary$ipwX2E=(Xsummary$XmeanSqE.ipw)*(Xsummary$sumWtE.ipw)
    Xsummary$stabX2E=(Xsummary$XmeanSqE.stab)*(Xsummary$sumWtE.stab)

    Xsummary$XUnE=(Xsummary$XmeanUnE)*(Xsummary$nUnE)
    Xsummary$ipwXUnE=(Xsummary$XmeanUnE.ipw)*(Xsummary$sumWtUnE.ipw)
    Xsummary$stabXUnE=(Xsummary$XmeanUnE.stab)*(Xsummary$sumWtUnE.stab)

    Xsummary$X2UnE=(Xsummary$XmeanSqUnE)*(Xsummary$nUnE)
    Xsummary$ipwX2UnE=(Xsummary$XmeanSqUnE.ipw)*(Xsummary$sumWtUnE.ipw)
    Xsummary$stabX2UnE=(Xsummary$XmeanSqUnE.stab)*(Xsummary$sumWtUnE.stab)

    XsummaryList[[k]]=subset(Xsummary, select=c("nE","sumWtE.ipw","sumWtE.stab",
                                                "nUnE","sumWtUnE.ipw","sumWtUnE.stab",
                                                "XE", "ipwXE","stabXE","X2E","ipwX2E","stabX2E",
                                                "XUnE","ipwXUnE","stabXUnE","X2UnE","ipwX2UnE","stabX2UnE"))


  }


  sumList=Reduce("+", XsummaryList)

  sumList$binary=XsummaryListCopy[[1]]$binary

  nCov=nrow(XsummaryList[[1]])

  output=as.data.frame(matrix(0,nCov,7))
  rownames(output)=rownames(XsummaryList[[1]])
  output[,1]=XsummaryListCopy[[1]]$binary
  colnames(output)=c("binary","exposed","unexposed",
                     "exposed.ipw","unexposed.ipw",
                     "exposed.stab","unexposed.stab")

  for (i in 1:nCov){

    if (sumList$binary[i]=="yes") {

      count=sumList$XE[i]
      prev=(sumList$XE[i])/(sumList$nE[i])*100
      count=round(count,digits=digits)
      prev=round(prev,digits=digits)
      output$exposed[i]=paste(count, "(",prev,"%)", sep ="")

      count=sumList$XUnE[i]
      prev=(sumList$XUnE[i])/(sumList$nUnE[i])*100
      count=round(count,digits=digits)
      prev=round(prev,digits=digits)
      output$unexposed[i]=paste(count, "(",prev,"%)", sep ="")

      count=sumList$ipwXE[i]
      prev=(sumList$ipwXE[i])/(sumList$sumWtE.ipw[i])*100
      count=round(count,digits=digits)
      prev=round(prev,digits=digits)
      output$exposed.ipw[i]=paste(count, "(",prev,"%)", sep ="")

      count=sumList$ipwXUnE[i]
      prev=(sumList$ipwXUnE[i])/(sumList$sumWtUnE.ipw[i])*100
      count=round(count,digits=digits)
      prev=round(prev,digits=digits)
      output$unexposed.ipw[i]=paste(count, "(",prev,"%)", sep ="")

      count=sumList$stabXE[i]
      prev=(sumList$stabXE[i])/(sumList$sumWtE.stab[i])*100
      count=round(count,digits=digits)
      prev=round(prev,digits=digits)
      output$exposed.stab[i]=paste(count, "(",prev,"%)", sep ="")

      count=sumList$stabXUnE[i]
      prev=(sumList$stabXUnE[i])/(sumList$sumWtUnE.stab[i])*100
      count=round(count,digits=digits)
      prev=round(prev,digits=digits)
      output$unexposed.stab[i]=paste(count, "(",prev,"%)", sep ="")

    } else {
      avg=(sumList$XE[i])/(sumList$nE[i])
      std=((sumList$X2E[i]-(avg^2)*(sumList$nE[i]))/(sumList$nE[i]-1))^0.5
      avg=round(avg,digits=digits)
      std=round(std,digits=digits)
      output$exposed[i]=paste(avg, "(",std,")", sep ="")

      avg=(sumList$XUnE[i])/(sumList$nUnE[i])
      std=((sumList$X2UnE[i]-(avg^2)*(sumList$nUnE[i]))/(sumList$nUnE[i]-1))^0.5
      avg=round(avg,digits=digits)
      std=round(std,digits=digits)
      output$unexposed[i]=paste(avg, "(",std,")", sep ="")

      avg=(sumList$ipwXE[i])/(sumList$sumWtE.ipw[i])
      std=((sumList$ipwX2E[i]-(avg^2)*(sumList$sumWtE.ipw[i]))/(sumList$sumWtE.ipw[i]-1))^0.5
      avg=round(avg,digits=digits)
      std=round(std,digits=digits)
      output$exposed.ipw[i]=paste(avg, "(",std,")", sep ="")

      avg=(sumList$ipwXUnE[i])/(sumList$sumWtUnE.ipw[i])
      std=((sumList$ipwX2UnE[i]-(avg^2)*(sumList$sumWtUnE.ipw[i]))/(sumList$sumWtUnE.ipw[i]-1))^0.5
      avg=round(avg,digits=digits)
      std=round(std,digits=digits)
      output$unexposed.ipw[i]=paste(avg, "(",std,")", sep ="")

      avg=(sumList$stabXE[i])/(sumList$sumWtE.stab[i])
      std=((sumList$stabX2E[i]-(avg^2)*(sumList$sumWtE.stab[i]))/(sumList$sumWtE.stab[i]-1))^0.5
      avg=round(avg,digits=digits)
      std=round(std,digits=digits)
      output$exposed.stab[i]=paste(avg, "(",std,")", sep ="")

      avg=(sumList$stabXUnE[i])/(sumList$sumWtUnE.stab[i])
      std=((sumList$stabX2UnE[i]-(avg^2)*(sumList$sumWtUnE.stab[i]))/(sumList$sumWtUnE.stab[i]-1))^0.5
      avg=round(avg,digits=digits)
      std=round(std,digits=digits)
      output$unexposed.stab[i]=paste(avg, "(",std,")", sep ="")

    }

  }


  output

}





