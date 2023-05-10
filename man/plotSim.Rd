\name{plotSim}
\alias{plotSim}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot the similarty between corelation values and P values;
}
\description{
Generates a plot of the similarty between corelation values and P values;
}
\usage{
plotSim(dx = d.Pval, dy = d.corr, clr = NULL, plot = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dx}{
A matrix of p values
}
  \item{dy}{
A matrix of correlation values
}
  \item{clr}{
Desired color
}
  \item{plot}{
Whether to generate the plot, defaults to True.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Ying Hu
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
cor=sapply(sim, function(x) x[,2])
Pval=sapply(sim, function(x) x[,1])
pdf ("SimPlotKnownTargets.pdf")
par(mar=c(4,4,5,2), xpd=TRUE, mfrow=c(1,1));
color <- colorRampPalette(c("lightblue",'darkblue'));
plotSim(dx=Pval,dy=cor,clr=color);
dev.off()
}
            
\keyword{ aplot }
