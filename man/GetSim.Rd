\name{GetSim}
\alias{GetSim}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Compute a correlation between the every gene crispr KO vs each drug response.
}
\description{

Computes correlations between the viability of cell lines after CRISPR knockout of each gene and the viability of the same cell lines after drug treatment.
}
\usage{
GetSim(drugName = DrugName, DRS = drugResponseScore, GES = GeneEffectScores)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{drugName}{
     Drug Name
}
  \item{DRS}{
      Drug’s response scores
}
  \item{GES}{
      Gene effect scores from knock-out method such as CRISPR.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A dataframe containing genes with their correlation values and P values associated with response scores from a given drug ID.
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
#saveRDS(sim,
  #      file = paste('drugVScrispr.Similarity.List.RDS', Sys.Date(), '.RDS', sep=''))
}
\keyword{ methods }
\keyword{ manip }
