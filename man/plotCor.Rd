\name{plotCor}
\alias{plotCor}
\title{
Plot the correlation
}
\description{
Plot the correlation of the predicted target
}
\usage{
plotCor(DN = Drugname, GN = GOI, Pred = Pred, DRS = DRS, GES = GES, plot = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{DN}{
Drug Name
}
  \item{GN}{
Gene Name
}
  \item{Pred}{
Output object from the prediction method
}
  \item{DRS}{
Drug response score
}
  \item{GES}{
Gene Effect scores
}
  \item{plot}{
Whether to generate the plot, defaults to True
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Correlation plot
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
plotCor ( DN=DOI,GN=GOI,Pred=Pred.d, DRS= sec.prism,GES= KO.GES)
}

\keyword{ aplot }
