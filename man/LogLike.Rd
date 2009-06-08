\name{LogLike}
\alias{LogLike}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Individual log-likelihood}
\description{
'LogLike' calculates the individual log-likelihood which is needed for 'Clarketest' and 'Vuongtest'.
}
\usage{
LogLike(Yin, fm.X, region, regmodel, thinning = 1, burnin = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Yin}{response vector of length n.}
  \item{fm.X}{formula for mean design.}
  \item{region}{region of each observation.}
  \item{regmodel}{the output of model fits obtained by using 'est.sc'.}
  \item{thinning}{number of MCMC iterations between two writing steps (if thinning=1 (default), all
states are calculated whereas if e.g. thinning=10 only each 10 iteration is calculated)}
  \item{burnin}{number of steps dropped from the chain to allow for a burn-in phase. Defaults to 1.}
}
%\details{
%  ~~ If necessary, more details than the description above ~~
%}
\value{
  \item{ll }{is a matrix of dimension (n x (totalit+2)), where each 
column represents one iteration and contains the likelihood 
contributions of each observation.}
}
%\references{ ~put references to the literature/web site here ~ }
%\author{ ~~who you are~~ }
%\note{ ~~further notes~~ 

% ~Make other sections like Warning with \section{Warning }{....} ~
%}
%\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }
\examples{
data(sim.Yin)
data(sim.fm.X)
data(sim.region)
data(sim.gmat)
data(sim.nmat)

poi <- est.sc(sim.Yin, ~1+sim.fm.X[,1]+sim.fm.X[,2], sim.region, 
model="Poi", sim.gmat, sim.nmat, 200)

ll.poi <- LogLike(sim.Yin, ~1+sim.fm.X[,1]+sim.fm.X[,2], sim.region, poi)

# log likelihood for the single iterations
apply(ll.poi$ll,2,sum)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% ONLY ONE keyword per line
