\name{plotSim}
\alias{plotSim}
\title{
Plot the similarty between corelation values and P vals for all genes. The top 5 genes are labeled.
}
\description{
Plot the similarty between corelation values and P val;
}
\usage{
plotSim(dx,dy,clr=NULL, plot=TRUE)
}
\arguments{
    \item{dx}{a matrix of p vals}
    \item{dy}{a matrix of correlation vals}
    \item{clr}{Desired range of color}
    \item{plot}{ default plot =TRUE}
}
\value{a plot of similarity}
\author{
Ying Hu,Trinh Nguyen
}
\examples{
library(BiocParallel)
data (OntargetM)
set.seed (12345)
All.Drugs <- OntargetM$DrugMetadata[,"broad_id_trimmed"]
Sample.Drugs <- sample(All.Drugs, 5)
KO.GES <- OntargetM$avana_CRISPR
sec.prism <- OntargetM$secondary_prism
sim.out <- bplapply(Sample.Drugs,function(x) computeCor(x,sec.prism,KO.GES))
names(sim.out) <- Sample.Drugs
P.Values=vapply(sim.out, function(x) x[,1],FUN.VALUE=numeric(nrow(sim.out[[1]])))
estimate.cor.values=vapply(sim.out, function(x) x[,2],FUN.VALUE=numeric(nrow(sim.out[[1]])))
par(mar=c(4,4,5,2), xpd=TRUE, mfrow=c(3,3));
plotSim(dx=P.Values,dy=estimate.cor.values);
}
