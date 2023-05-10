\name{line2user}
\alias{line2user}
\title{
A utility function to correctly place the colorBar.
}
\description{
This function is used in conjunction with the colorBar function to create color bar placed in top middle of a plot.
}
\usage{
line2user(line, side)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{line}{
%%     ~~Describe \code{line} here~~
}
  \item{side}{
%%     ~~Describe \code{side} here~~
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
Ying Hu <yhu@mail.nih.gov>
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
function (line, side)
{
    lh <- par("cin")[2] * par("cex") * par("lheight")
    x_off <- diff(grconvertX(0:1, "inches", "user"))
    y_off <- diff(grconvertY(0:1, "inches", "user"))
    switch(side, `1` = par("usr")[3] - line * y_off * lh, `2` = par("usr")[1] -
        line * x_off * lh, `3` = par("usr")[4] + line * y_off *
        lh, `4` = par("usr")[2] + line * x_off * lh, stop("side must be 1, 2, 3, or 4",
        call. = FALSE))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
\keyword{ utilities }
\keyword{ internal }
\keyword{ aplot }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
