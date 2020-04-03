\name{tess2spat}
\alias{tess2spat}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Convert Voronoi tesselation tiles to a shapefile.
}
\description{
Convert Voronoi tesselation tiles to a shapefile.
}
\usage{
tess2spat(obj, idvec=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{
Voronoi tesselation object created through function \code{deldir}.
}
  \item{idvec}{
Optional vector of ids for output shapefile polygons.
}
}
\details{
}
\value{
Object of class \code{SpatialPolygons}.
}
\references{
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

library(deldir)
library(sp)

vortess <- deldir(x=runif(8), y=runif(8), plotit=FALSE, suppressMsge=TRUE)
par(mfcol=c(1,2))
plot.deldir(vortess, wlines="tess", xlim=c(0,1), ylim=c(0,1))
vortess_shape <- tess2spat(obj=vortess)
plot(vortess_shape)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ spatial }% use one of  RShowDoc("KEYWORDS")