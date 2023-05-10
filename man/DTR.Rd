\name{DTR}
\alias{DTR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot of drug target response for primary or secondary targets.
}
\description{
Generates a plot of drug target response for primary or secondary targets based on a drug of interest and a gene of interest.
}
\usage{
DTR(DN = DOI, GN = GOI, Pred = Pred, Exp = Exp, DRS = DRS, GES = GES, CutOff = 3,plot = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{DN}{
Drug of interest
}
  \item{GN}{
Gene of interest
}
  \item{Pred}{
Prediction object, an output result from prediction whether it is a primary target or secondary target
}
  \item{Exp}{
Expression matrix
}
  \item{DRS}{
Drug response matrix
}
  \item{GES}{
Gene Effect Scores
}
  \item{plot}{
Whether users want to plot, default is true
}
  \item{CutOff}{
Cutoff value for gene expression of gene of interest, high or low
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
vialbility after KO vs drug response in GOI low vs high
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
 sanjusinha7
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
data (OntargetM)
set.seed (12345)
S.Drugs <- c('K70301465','K09951645')
KO.GES <- OntargetM$avana_CRISPR
sec.prism <- OntargetM$secondary_prism
d.expr <- OntargetM$expression_20Q4
sim <- mclapply(S.Drugs,function(x) GetSim(x,DRS=sec.prism, GES=KO.GES),mc.cores = 2)
names(sim) <- S.Drugs
Meta.data <- OntargetM$DrugMetadata
DrugTargetSim <- PredTarget(sim,Meta.data)
Drug.Gene.max.sim <- PredMaxSim(Sim.GES.DRS=sim, D.M = Meta.data)
identical ( DrugTargetSim[,1],Drug.Gene.max.sim[,1] )
Pred.d <-cbind (DrugTargetSim,Drug.Gene.max.sim )
DOI = 'ibrutinib'
GOI ='BTK'
DTR ( DN=DOI,GN=GOI,Pred=Pred.d, Exp=d.expr,DRS= sec.prism,GES= KO.GES,CutOff= 2)

}
\keyword{ aplot }
\keyword{ regression }
