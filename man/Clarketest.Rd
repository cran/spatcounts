\name{Clarketest}
\alias{Clarketest}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Clarke's test for non-nested model comparison}
\description{
'Clarketest' suggests the better of two (not necessarily nested) models according to Clarke's statistic for the parameters in each of the iterations.
}
\usage{
Clarketest(LogLike1, LogLike2, alpha = 0.05, p = NULL, q = NULL, 
correction = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{LogLike1, LogLike2}{the output of two model fits obtained by using 'LogLike'.}
  \item{alpha}{significance level, defaults to 0.05.}
  \item{p, q}{the number of estimated coefficients in models LogLike1 and Loglike2, respectively.}
  \item{correction}{boolean, if TRUE (default), the Schwarz correction will be used on the differences of log-likelihoods.}
}
%\details{
%  ~~ If necessary, more details than the description above ~~
%}
%\value{
%  ~Describe the value returned
%  If it is a LIST, use
%  \item{comp1 }{Description of 'comp1'}
%  \item{comp2 }{Description of 'comp2'}
%  ...
%}
\references{Clarke, Kevin A. (2003). Nonparametric Model Discrimination in International Relations. Journal
of Conflict Resolution 47(1), 72-93.

Schwarz, G. (1978). Estimating the Dimension of a Model. Annals of Statistics 6, 461-464.}

%\author{ ~~who you are~~ }
%\note{ ~~further notes~~ 

% ~Make other sections like Warning with \section{Warning }{....} ~
%}
\seealso{Vuongtest}
\examples{
data(sim.Yin)
data(sim.fm.X)
data(sim.region)
data(sim.gmat)
data(sim.nmat)

poi <- est.sc(sim.Yin, ~1+sim.fm.X[,1]+sim.fm.X[,2], sim.region, 
model="Poi", sim.gmat, sim.nmat, 3)
nb <- est.sc(sim.Yin, ~1+sim.fm.X[,1]+sim.fm.X[,2], sim.region, 
model="NB", sim.gmat, sim.nmat, 3)

DIC.poi <- DIC(sim.Yin, ~1+sim.fm.X[,1]+sim.fm.X[,2], sim.region, poi)
DIC.nb <- DIC(sim.Yin, ~1+sim.fm.X[,1]+sim.fm.X[,2], sim.region, nb)

ll.poi <- LogLike(sim.Yin, ~1+sim.fm.X[,1]+sim.fm.X[,2], sim.region, poi)
ll.nb <- LogLike(sim.Yin, ~1+sim.fm.X[,1]+sim.fm.X[,2], sim.region, nb)

Clarke.poi.nb <- Clarketest(ll.poi, ll.nb, alpha = 0.05, p = DIC.poi$p.D, 
q = DIC.nb$p.D, correction = TRUE)
}
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }%ONLY ONE keyword per line
