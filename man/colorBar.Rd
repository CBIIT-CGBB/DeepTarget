\name{colorBar}
\alias{colorBar}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Create a color key bar.
}
\description{
Utility method to create a color key bar for a plot.
}
\usage{
colorBar(xl, yb, xr, yt, v.min, v.max, col.l, alpha = 1, color.spec = "hsv", title = "")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{xl}{
x axis, left
}
  \item{yb}{
y axis, bottom
}
  \item{xr}{
x, right
}
  \item{yt}{
y, top
}
  \item{v.min}{
mininum value
}
  \item{v.max}{
maximum value
}
  \item{col.l}{
%%     ~~Describe \code{col.l} here~~
}
  \item{alpha}{
%%     ~~Describe \code{alpha} here~~
}
  \item{color.spec}{
%%     ~~Describe \code{color.spec} here~~
}
  \item{title}{
Color bar title
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

#par(mar=c(4,4,5,2), xpd=TRUE);
#plot (NULL);
#rbPal <- colorRampPalette(c("lightblue",'darkblue'))
#colorBar(6, -0.02, 9, 0.02, col.l=rbPal(100), v.min=-1, v.max=1);

## The function is currently defined as
function (xl, yb, xr, yt, v.min, v.max, col.l, alpha = 1, color.spec = "hsv",
    title = "")
{
    nticks <- 11
    ticks <- seq(v.min, v.max, len = nticks)
    tmp.v <- seq(v.min, v.max, length.out = 100)
    lut <- col.l
    xs <- (xr - xl)/(length(lut) - 1)
    xl.old <- xl
    for (i in 1:(length(lut) - 1)) {
        rect(xl.old, yb, xl.old + xs, yt, col = lut[i], border = NA)
        xl.old <- xl.old + xs
    }
    v.med <- (v.max + v.min)/2
    xm <- (xl + xr)/2
    v.max <- round(v.max, 2)
    v.min <- round(v.min, 2)
    v.med <- round(v.med, 2)
    x.d <- xr - xl
    y.d <- yt - yb
    text(xl, yt + x.d, title, cex = 1)
    text(xr, yb - y.d, v.max, cex = 0.6)
    text(xl, yb - y.d, v.min, cex = 0.6)
    text(xm, yb - y.d, v.med, cex = 0.6)
    segments(xl, yb - y.d/5, xr, yb - y.d/5, col = "black", lwd = 0.8)
    segments(xl, yb - y.d/5, xl, yb - y.d/2, col = "black", lwd = 0.8)
    segments(xm, yb - y.d/5, xm, yb - y.d/2, col = "black", lwd = 0.8)
    segments(xr, yb - y.d/5, xr, yb - y.d/2, col = "black", lwd = 0.8)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
\keyword{ aplot }
\keyword{ color }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
