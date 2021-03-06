% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/parafac_functions.R
\name{eempf_excomp}
\alias{eempf_excomp}
\title{Extracting components of a PARAFAC model}
\usage{
eempf_excomp(pfmodel, comps)
}
\arguments{
\item{pfmodel}{parafac model}

\item{comps}{vector with numbers of components to extract}
}
\value{
list
}
\description{
Extracting components of a PARAFAC model
}
\examples{
data(pf_models)
pfmodel <- pf4[[1]]
comps <- eempf_excomp(pfmodel,c(1,3))
}
