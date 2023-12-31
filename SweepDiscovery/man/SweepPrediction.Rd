% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SweepFeatures.R
\name{SweepPrediction}
\alias{SweepPrediction}
\title{SweepPrediction}
\usage{
SweepPrediction(Datapath, Chromosome, Start, Finish)
}
\arguments{
\item{Datapath}{Datasets}

\item{Chromosome}{Chromosome Number}

\item{Start}{Start Number}

\item{Finish}{Finish Number}
}
\value{
\itemize{
  \item Prediction: Results
}
}
\description{
SweepPrediction
}
\examples{
\donttest{
library("SweepDiscovery")
Path<-file.path(system.file("exdata/ExampleVCF.vcf.gz", package = "SweepDiscovery"))
pred<-SweepPrediction(Datapath=Path,"1",20253,1976067)
}
}
\references{
\itemize{
\item Pavlidis, P., Alachiotis, N. A survey of methods and tools to detect recent and strong positive selection. J of Biol Res-Thessaloniki 24, 7 (2017). https://doi.org/10.1186/s40709-017-0064-0
}
}
