% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stareg.R
\name{stareg}
\alias{stareg}
\title{An empirical Bayes approach for replicability analysis across two studies}
\usage{
stareg(pa, pb, init.pi0 = TRUE)
}
\arguments{
\item{pa}{A numeric vector of p-values from study 1.}

\item{pb}{A numeric vector of p-values from study 2.}

\item{init.pi0}{A logistic value for deciding whether to initialize the prior probabilities based on the estimates of pi0's. If true, estimate the marginal pi0's in two studies using qvalue; otherwise, specify pi0_pa = pi_pb = 0.9.}
}
\value{
A list:
\item{Lfdr}{The estimated local false discovery rate for replicability null.}
\item{fdr}{The adjusted Lfdr values based on the step-up procedure for FDR control.}
\item{xi00}{An estimate of the prior probability for joint state (0, 0) in two studies.}
\item{xi01}{An estimate of the prior probability for joint state (0, 1) in two studies.}
\item{xi10}{An estimate of the prior probability for joint state (1, 0) in two studies.}
\item{xi11}{An estimate of the prior probability for joint state (1, 1) in two studies.}
\item{f1}{A non-parametric estimate for the non-null probability density function in study 1.}
\item{f2}{A non-parametric estimate for the non-null probability density function in study 2.}
}
\description{
An empirical Bayes approach for replicability analysis across two studies
}
\examples{
# Simulate p-values in two studies
m = 10000
h = sample(0:3, m, replace = TRUE, prob = c(0.9, 0.025, 0.025, 0.05))
states1 = rep(0, m); states2 = rep(0, m)
states1[which(h==2|h==3)] = 1; states2[which(h==1|h==3)] = 1
z1 = rnorm(m, states1*2, 1)
z2 = rnorm(m, states2*3, 1)
p1 = 1 - pnorm(z1); p2 = 1 - pnorm(z2)
# Run STAREG to identify replicable signals
res.stareg = stareg(p1, p2)
sig.idx = which(res.stareg$fdr <= 0.05)

}
