\name{DMB}
\alias{DMB}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{

Plot of drug target response for mutant and non-mutant groups.
}
\description{

Generates a plot of Drug target response for primary or secondary targets based on a drug of interest.
}
\usage{
DMB(DN = Drugname, GN = GOI, Pred = Pred, Mutant = Mutant, DRS = DRS, GES = GES, plot = TRUE)
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
  \item{Mutant}{
Mutant matrix
}
  \item{DRS}{
Drug response matrix
}
  \item{GES}{
Gene Effect Scores
}
  \item{plot}{
%%     ~~Describe \code{plot} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
vialbility after KO vs drug response in mutant.
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

S.Drugs <- c('K70301465','K09951645')
KO.GES <- OntargetM$avana_CRISPR
sec.prism <- OntargetM$secondary_prism
d.mt <- OntargetM$mutations_mat
sim <- mclapply(S.Drugs,function(x) GetSim(x,DRS=sec.prism, GES=KO.GES),mc.cores = 2)
names(sim) <- S.Drugs
Meta.data <- OntargetM$DrugMetadata
DrugTargetSim <- PredTarget(sim,Meta.data)
Drug.Gene.max.sim <- PredMaxSim(Sim.GES.DRS=sim, D.M = Meta.data)
identical ( DrugTargetSim[,1],Drug.Gene.max.sim[,1] )
Pred.d <-cbind (DrugTargetSim,Drug.Gene.max.sim )
DOI = 'dabrafenib'
GOI = 'BRAF'
DMB (DN=DOI,GN=GOI,Pred=Pred.d,Mutant=d.mt,DRS= sec.prism,GES= KO.GES)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
\keyword{ aplot }
\keyword{ regression }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
