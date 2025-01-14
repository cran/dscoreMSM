% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dscore.R
\name{dscore}
\alias{dscore}
\title{Survival Proximity Score matching for MSM}
\usage{
dscore(status, data, prob, m, n, method = "euclidean")
}
\arguments{
\item{status}{status column name in the survival data}

\item{data}{survival data}

\item{prob}{threshold probability}

\item{m}{starting column number}

\item{n}{ending column number}

\item{method}{distance metric name e.g. "euclidean","minkowski","canberra"}
}
\value{
list with newdataset updated using dscore
}
\description{
function for survival proximity score matching in multistate model with three state.
}
\examples{
\donttest{
##s
data(simulated_data)
udata<-dscore(status="status",data=simulated_data,prob=0.65,m=4,n=7)
##
}
}
\references{
Vishwakarma, G. K., Bhattacherjee, A., Rajbongshi, B. K., & Tripathy, A. (2024). Censored imputation of time to event outcome through survival proximity score method. \emph{Journal of Computational and Applied Mathematics}, 116103;

Bhattacharjee, A., Vishwakarma, G. K., Tripathy, A., & Rajbongshi, B. K. (2024). Competing risk multistate censored data modeling by propensity score matching method. \emph{Scientific Reports}, 14(1), 4368.
}
\seealso{
\link{cphGM},\link{simfdata}
}
\author{
Atanu Bhattacharjee, Bhrigu Kumar Rajbongshi and Gajendra K. Vishwakarma
}
