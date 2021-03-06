% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sim.taxonomy.R
\name{sim.cryptic.species}
\alias{sim.cryptic.species}
\title{Simulate cryptic species on a taxonomy object}
\usage{
sim.cryptic.species(species, kappa)
}
\arguments{
\item{species}{Taxonomy object.}

\item{kappa}{Probability that speciation event is cryptic.}
}
\value{
An object of class taxonomy. Note the origin or root can not be cryptic.
}
\description{
Simulate cryptic species on a taxonomy object
}
\examples{
t = ape::rtree(10)
sp = sim.taxonomy(t, 1)
sim.cryptic.species(sp, 0.5)

}
\seealso{
\code{\link{taxonomy}}
}
