\name{DIC}
\alias{DIC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Deviance information criterion}
\description{
Deviance information criterion without standardizing constant.
}
\usage{
DIC(Yin, fm.X, region, regmodel, burnin = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Yin}{response vector of length n.}
  \item{fm.X}{formula for mean design.}
  \item{region}{region of each observation.}
  \item{regmodel}{the output of model fits obtained by using 'est.sc'.}
  \item{burnin}{number of steps dropped from the chain to allow for a burn-in phase. Defaults to 1.}
}
%\details{
%  ~~ If necessary, more details than the description above ~~
%}
\value{
  \item{DIC }{deviance information criterion.}
  \item{mean.deviance }{posterior mean of deviance.}
  \item{p.D }{effective number of parameters.}
}
\references{Gelman, Andrew, John B. Carlin, Hall S. Stern, and Donald B. Rubin (2003). Bayesian Data Analysis, Second Edition. Chapman \& Hall/CRC.

Spiegelhalter, David J.; Nicola G. Best, Bradley P. Carlin, and Angelika van der Linde (2002). Bayesian measures of model complexity and fit. Journal of the Royal Statistical Society, Series B. 64 (4), 583-639.}
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
model="Poi", sim.gmat, sim.nmat, 3)
DIC.poi <- DIC(sim.Yin, ~1+sim.fm.X[,1]+sim.fm.X[,2], sim.region, poi)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% ONLY ONE keyword per line
