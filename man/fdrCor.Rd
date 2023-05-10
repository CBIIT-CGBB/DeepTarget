\name{fdrCor}
\alias{fdrCor}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 Calculate False Discovery Rate (FDR)
}
\description{
This function is used to calculate for FDR based on the list of p values.
}
\usage{
 fdrCor(test_list)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{test_list}{
A list of p values.
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
## The function is currently defined as
function (test_list)
{
    p.adjust(test_list, method = "fdr")
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
\keyword{ utilities }
\keyword{ internal }
\keyword{ htest }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
