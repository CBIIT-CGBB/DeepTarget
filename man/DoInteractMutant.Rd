\name{DoInteractMutant}
\alias{DoInteractMutant}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Compute interaction between the drug and KO expression
}
\description{
Computes the interaction between the drug and KO expression in terms of mutant vs non-mutant.
}
\usage{

DoInteractMutant(Predtargets = Drug.targeted.sim, Mutant = Mutant, DRS = DRS, GES = GES)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
 \item{Predtargets}{
  A table of drug information and their most targeted gene with stats of correlation}
\item{Mutant}{
 Mutant matrix}
\item{DRS}{
 Drug scores}
\item{GES}{
 Gene effect scores from KO method}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A list of drug names with their interaction values from two groups mutant and non-mutant
    \item{drug1 }{interaction with estimate and P vals}
    \item{drug2 }{interaction with estimate and P vals}
    \item{drugN }{interaction with estimate and P vals}
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
DrugTargetSim <- PredTarget(sim, D.M = Meta.data)
d.mt <- OntargetM$mutations_mat
MutantInteract <- DoInteractMutant (Predtargets=DrugTargetSim,Mutant=d.mt,DRS=sec.prism,GES=KO.GES)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
\keyword{ manip }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
