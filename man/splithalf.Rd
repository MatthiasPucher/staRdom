% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/parafac_functions.R
\name{splithalf}
\alias{splithalf}
\title{Running a Split-Half analysis on a PARAFAC model}
\usage{
splithalf(
  eem_list,
  comps,
  splits = NA,
  rand = FALSE,
  normalise = TRUE,
  nstart = 20,
  cores = parallel::detectCores(logical = FALSE),
  maxit = 2500,
  ctol = 10^(-7),
  rescale = TRUE,
  strictly_converging = FALSE,
  verbose = FALSE,
  ...
)
}
\arguments{
\item{eem_list}{eemlist containing sample data}

\item{comps}{number of desired components}

\item{splits}{optional, list of 4 numerical vectors containing the sample numbers for A,B,C and D sample subsets}

\item{rand}{logical, splits are randomised}

\item{normalise}{state whether EEM data should be normalised in advance}

\item{nstart}{number of random starts}

\item{cores}{number of parallel calculations (e.g. number of physical cores in CPU)}

\item{maxit}{maximum iterations for PARAFAC algorithm}

\item{ctol}{Convergence tolerance (R^2 change)}

\item{rescale}{rescale splithalf models to Fmax, see \code{\link[staRdom]{eempf_rescaleBC}}}

\item{strictly_converging}{calculate nstart converging models and take the best. Please see \code{\link[staRdom]{eem_parafac}}.}

\item{verbose}{states whether you want additional information during calculation}

\item{...}{additional parameters that are passed on to \code{\link[multiway]{parafac}}}
}
\value{
data frame containing components of the splithalf models
}
\description{
The samples are split into four subsamples: A,B,C,D. Subsamples are then combined and compared: AB vs. CD, AC vs. BD, AD vs. BC. The results show graphs from the components of each of the 6 models.
}
\details{
Split data sets can be split suboptimal and cause low TCCs. Therefore, subsamples are recombined in 3 different ways and a TCC close to 1 in only one split combination per component is already a positive result. Check the split sets to check for sample independency.
}
\examples{
\donttest{
data(eem_list)

splithalf <- splithalf(eem_list, comps = 6, verbose = TRUE, cores = 2)
splithalf_plot(splithalf)

# Similarity of splits using SSCs
sscs <- splithalf_tcc(splithalf)
}
}
\seealso{
\code{\link[staRdom]{splithalf_plot}}, \code{\link[staRdom]{splithalf_tcc}}
}
