\name{ppmHR-package}
\alias{ppmHR-package}
\docType{package}
\title{
Privacy-Protecting Hazard Ratio Estimation in Distributed Data Networks
}
\description{
This package allows implementation of a one-step, privacy-protecting method for estimating the overall and site-specific hazard ratios using inverse probability weighted Cox models in multi-center distributed data network studies, as proposed by Shu, Yoshida, Fireman, and Toh (2019). The method only requires sharing of summary-level riskset tables instead of individual-level data by data-contributing sites. Both the conventional inverse probability weights and the stabilized weights are implemented.
}
\details{
The \code{ppmHR} package implements a one-step, privacy-protecting method for estimating overall and site-specific hazard ratios using inverse probability weighted Cox models, under both the conventional inverse probability weights and the stabilized weights, as proposed by Shu, Yoshida, Fireman, and Toh (2019). The function \code{\link[ppmHR]{checkBalanceSite}} allows the data-contributing sites to check their site-specific covariate balance before and after weighting. The function \code{\link[ppmHR]{computeInfoForTable1}} allows the data-contributing sites to compute summary-level information that will be used to create the "Table 1" of baseline patient characteristics before and after weighting. The function \code{\link[ppmHR]{createTable1}} allows the analysis center to create the "Table 1" of baseline patient characteristics before and after weighting, using only summary-level information shared by data-contributing sites. The function \code{\link[ppmHR]{createRisksetTable}} allows the data-contributing sites to create summary-level riskset tables using their site-specific individual-level data. The function \code{\link[ppmHR]{estimateStratHR}} allows the analysis center to estimate the overall hazard ratios and robust sandwich variances using inverse probability weighted Cox models stratified on data-contributing sites, using only the summary-level riskset tables provided by the data-contributing sites. The function \code{\link[ppmHR]{estimateSiteHRs}} allows the analysis center to estimate the site-specific hazard ratios and robust sandwich variances using inverse probability weighted Cox models, using only the summary-level riskset tables provided by the data-contributing sites.
}
\author{
Di Shu and Sengwee Toh

Maintainer: Di Shu <shudi1991@gmail.com>
}
\references{
Shu D, Yoshida K, Fireman BH, Toh S (2019).  Inverse probability weighted Cox model in multi-site studies without sharing individual-level data. \emph{Statistical Methods in Medical Research} <doi:10.1177/0962280219869742>
}
\keyword{ package }

