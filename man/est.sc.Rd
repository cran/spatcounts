\name{est.sc}
\alias{est.sc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Fitting spatial count regression models}
\description{
MCMC algorithm for the Poisson, the GP, the NB, the ZIP and the ZIGP regression
models with spatial random effects.
}
\usage{
est.sc(Yin, fm.X, region, model = "Poi", gmat, nmat, totalit, 
fm.ga = TRUE, t.i = NULL, phi0 = 1, omega0 = 0, r0 = 1, 
beta0 = NULL, gamma0 = NULL, sigma0 = 1, psi0 = 1, Tau = 10, 
alpha = 2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Yin}{response vector of length n.}
  \item{fm.X}{formula for mean design.}
  \item{region}{region of each observation as vector of length n.}
  \item{model}{the regression model. Currently, the regression models "Poi", "NB", "GP", "ZIP" and "ZIGP" are supported. Defaults to 'Poi'.}
  \item{gmat}{spatial adjacency matrix, where entry (i,j) is 1 if region i is a 
neighbor of region j and 0 else. See data(sim.gmat) for an example.}
  \item{nmat}{matrix containing the number of neighbors of each region (last 
column) and the neighbors of each region (first (maximual number of 
neighbours) columns), filled up by zeros. See data(sim.nmat) for an example.}
  \item{totalit}{number of MCMC iterations, i.e. length of the Markov chain.}
  \item{fm.ga}{should the spatial random effects be included (defaults to TRUE)?}
  \item{t.i}{exposure vector.}
  \item{phi0}{starting value for the over-dispersion parameter for GP and ZIGP model.}
  \item{omega0}{starting value for the extra proportion for ZIP and ZIGP model.}
  \item{r0}{starting value for the scale paramter for NB model.}
  \item{beta0}{starting values for the regression parameters.}
  \item{gamma0}{starting values for the spatial paramters.}
  \item{sigma0}{starting value for the spatial hyperparamter of CAR prior.}
  \item{psi0}{starting value for the spatial hyperparamter of CAR prior.}
  \item{Tau}{modifiable normal prior for the regression parameters with variance Tau$^2$.}
  \item{alpha}{modifiable prior distribution of hyperparamter psi (suggested values: 2, 1.5, 1, 0.5, 0).}
}
%\details{
%  ~~ If necessary, more details than the description above ~~
%}
\value{
  \item{acceptb }{acceptance rate for the regression parameters beta.}
  \item{acceptga }{range of the acceptance rate for the spatial parameters gamma.}
  \item{acceptphi }{acceptance rate for the GP and ZIGP model specific dispersion parameter phi.}
  \item{acceptomega }{acceptance rate for the ZIP and ZIGP model specific extra proportion omega.}
  \item{acceptr }{acceptance rate for the NB model specific scale parameter r.}
  \item{acceptpsi }{acceptance rate for the spatial hyperparameter psi.}
  \item{beta }{are the parameter estimates for the regression parameters beta.}
  \item{gamma }{are the parameter estimates for the spatial parameters gamma.}
  \item{invsigsq }{are the parameter estimates for the inverse spatial hyperparameter sigma$^2$.}
  \item{psi }{are the parameter estimates for the spatial hyperparameter psi.}
  \item{phi }{are the parameter estimates for the GP and ZIGP model specific dispersion parameter phi.}
  \item{omega }{are the parameter estimates for the ZIP and ZIGP model specific extra proportion omega.}
  \item{r }{are the parameter estimates for the NB model specific scale parameter r.}
  \item{Coefficients }{are the number of parameter estimates.}
  \item{t.i }{exposure vector.}
}
\references{Gschloessl, Susanne (2006). Hierarchical Bayesian spatial regression models with applications to non-life insurance. Dissertation, Centre of Mathematical Sciences, Munich University of Technology, Chair in Mathematical Statistics, Munich University of Technology,
Boltzmannstr. 3, D-85748 Garching near Munich.

Masterthesis: 
Schabenberger, Holger (2009). Spatial count regression models with applications 
to health insurance data. ("http://www-m4.ma.tum.de/Diplarb/").

Czado, C., Erhardt, V., Min, A., Wagner, S. (2007).
Zero-inflated generalized Poisson models with regression effects on the mean, 
dispersion and zero-inflation level applied to patent outsourcing rates.
Statistical Modelling 7 (2), 125-153. 
}
%\author{ ~~who you are~~ }
%\note{ ~~further notes~~ 

% ~Make other sections like Warning with \section{Warning }{....} ~
%}
\seealso{ R-package ZIGP for fitting GP, ZIP, ZIGP regression models using MLE. }
\examples{
data(sim.Yin)
data(sim.fm.X)
data(sim.region)
data(sim.gmat)
data(sim.nmat)
# true parameters for generating this data:
# beta.true = c(-1, 0.4, 1.5)
# gamma.true = vector of spatial effects according to the CAR model with mean 0, psi = 3 and sigma = 1
# range of gamma.true = c(-0.851, 0.8405) 

poi <- est.sc(sim.Yin, ~1+sim.fm.X[,1]+sim.fm.X[,2], sim.region,
model="Poi", sim.gmat, sim.nmat, totalit=200)

# posterior means not considering a burn-in or thinning of iterations
apply(poi$beta,1,mean)
apply(poi$gamma,1,mean)


# Compare Poisson to different model classes
nb <- est.sc(sim.Yin, ~1+sim.fm.X[,1]+sim.fm.X[,2], sim.region, model="NB", sim.gmat, sim.nmat, 200)

gp <- est.sc(sim.Yin, ~1+sim.fm.X[,1]+sim.fm.X[,2], sim.region, model="GP", sim.gmat, sim.nmat, 200)

zip <- est.sc(sim.Yin, ~1+sim.fm.X[,1]+sim.fm.X[,2], sim.region, model="ZIP", sim.gmat, sim.nmat, 200)

zigp <- est.sc(sim.Yin, ~1+sim.fm.X[,1]+sim.fm.X[,2], sim.region, model="ZIGP", sim.gmat, sim.nmat, 200)

DIC.poi <- DIC(sim.Yin, ~1+sim.fm.X[,1]+sim.fm.X[,2], sim.region, poi)
DIC.nb <- DIC(sim.Yin, ~1+sim.fm.X[,1]+sim.fm.X[,2], sim.region, nb)
DIC.gp <- DIC(sim.Yin, ~1+sim.fm.X[,1]+sim.fm.X[,2], sim.region, gp)
DIC.zip <- DIC(sim.Yin, ~1+sim.fm.X[,1]+sim.fm.X[,2], sim.region, zip)
DIC.zigp <- DIC(sim.Yin, ~1+sim.fm.X[,1]+sim.fm.X[,2], sim.region, zigp)

ll.poi <- LogLike(sim.Yin, ~1+sim.fm.X[,1]+sim.fm.X[,2], sim.region, poi)
ll.nb <- LogLike(sim.Yin, ~1+sim.fm.X[,1]+sim.fm.X[,2], sim.region, nb)
ll.gp <- LogLike(sim.Yin, ~1+sim.fm.X[,1]+sim.fm.X[,2], sim.region, gp)
ll.zip <- LogLike(sim.Yin, ~1+sim.fm.X[,1]+sim.fm.X[,2], sim.region, zip)
ll.zigp <- LogLike(sim.Yin, ~1+sim.fm.X[,1]+sim.fm.X[,2], sim.region, zigp)

Vuong.poi.nb <- Vuongtest(ll.poi, ll.nb, alpha = 0.05, p = DIC.poi$p.D, q = DIC.nb$p.D, correction = TRUE)
Vuong.poi.gp <- Vuongtest(ll.poi, ll.gp, alpha = 0.05, p = DIC.poi$p.D, q = DIC.gp$p.D, correction = TRUE)
Vuong.poi.zip <- Vuongtest(ll.poi, ll.zip, alpha = 0.05, p = DIC.poi$p.D, q = DIC.zip$p.D, correction = TRUE)
Vuong.poi.zigp <- Vuongtest(ll.poi, ll.zigp, alpha = 0.05, p = DIC.poi$p.D, q = DIC.zigp$p.D, correction = TRUE)

Clarke.poi.nb <- Clarketest(ll.poi, ll.nb, alpha = 0.05, p = DIC.poi$p.D, q = DIC.nb$p.D, correction = TRUE)
Clarke.poi.gp <- Clarketest(ll.poi, ll.gp, alpha = 0.05, p = DIC.poi$p.D, q = DIC.gp$p.D, correction = TRUE)
Clarke.poi.zip <- Clarketest(ll.poi, ll.zip, alpha = 0.05, p = DIC.poi$p.D, q = DIC.zip$p.D, correction = TRUE)
Clarke.poi.zigp <- Clarketest(ll.poi, ll.zigp, alpha = 0.05, p = DIC.poi$p.D, q = DIC.zigp$p.D, correction = TRUE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% ONLY ONE keyword per line
