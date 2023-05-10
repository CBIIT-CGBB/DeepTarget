\name{DoPWY}
\alias{DoPWY}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Provide a probability score for each pathway to be the MOA of a drug
}
\description{
Predicts a Primary Target at a pathwya Level. It next finds the pathways that are most enriched in the genes with high DKS scores. It does this by performing a pathway enrichment test on the ranked gene list by DKS score. The output is a matrix of pathway-level probabilities for each drug to be the primary MOA
}
\usage{
DoPWY(Sim.GES.DRS = Sim.GES.DRS, D.M = Drug.Metadata)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{Sim.GES.DRS}{
The list of result from "CorCRISPRDrugR" function.
}
\item{D.M}{
meta data from drug
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
a dataframe contain the pathway level probability to be a primary MOA
%%  ~~ If necessary, more details than the description above ~~
    \item{drug1 }{a dataframe contain the pathway level probability to be a primary MOA}
    \item{drug2 }{a dataframe contain the pathway level probability to be a primary MOA}
    \item{drugN }{a dataframe contain the pathway level probability to be a primary MOA}
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
Pwy.Enr <- DoPWY(Sim.GES.DRS=sim,D.M = Meta.data)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
\keyword{ manip }
\keyword{ distribution }
\keyword{ methods }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
