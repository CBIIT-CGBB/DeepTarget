\name{PredMaxSim}
\alias{PredMaxSim}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 Predict the most similar gene to the drug response
}
\description{

Predicts the gene that has the most similarity associated with drug’s response scores from the set of all genes.

}
\usage{
PredMaxSim (Sim.GES.DRS = Sim.GES.DRS, D.M = Drug.Metadata)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{Sim.GES.DRS}{
Similarity between Drug’s response scores and Gene effect scores from Knock-out method such as CRISPR.
}
\item{D.M}{
Drug Metadata
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A table of drug information and the most correlated gene across all genes.
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
All.Drugs <- OntargetM$DrugMetadata[,"broad_id_trimmed"]
S.Drugs <- sample(All.Drugs, 5)
KO.GES <- OntargetM$avana_CRISPR
sec.prism <- OntargetM$secondary_prism
sim <- mclapply(S.Drugs,function(x) GetSim(x,DRS=sec.prism, GES=KO.GES),mc.cores = 2)
names(sim) <- S.Drugs
Meta.data <- OntargetM$DrugMetadata
Drug.Gene.max.sim <- PredMaxSim(Sim.GES.DRS=sim, D.M = Meta.data)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
\keyword{ methods }
\keyword{ manip }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
