\name{DoInteractExp}
\alias{DoInteractExp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Compute the interaction between the drug and KO expression
}
\description{
Computes the interaction between the drug and KO expression in terms of lower vs higher expression
}
\usage{
DoInteractExp(Predtargets = Drug.sim, Exp = Exp, DRS = DRS, GES = GES, CutOff = 3)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
 \item{Predtargets}{
  A table of drug information and their most targeted gene with stats of correlation}

\item{Exp}{
 Expression matrix}

\item{DRS}{
 Drug scores}

\item{GES}{
 Gene effect scores from KO method}

\item{CutOff}{
 desired cut-off for low expression}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A list of drug names with their interaction values from two groups low and high expression.
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
d.expr <- OntargetM$expression_20Q4
ExpInteract <- DoInteractExp (DrugTargetSim,d.expr,sec.prism,KO.GES,CutOff = 2)
}
\keyword{ manip }
