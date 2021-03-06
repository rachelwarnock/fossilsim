% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BEAST.output.R
\name{fossils.to.BEAST.start.tree}
\alias{fossils.to.BEAST.start.tree}
\title{Create a suitable starting tree for a DPPDIV style FBD analysis in BEAST2}
\usage{
fossils.to.BEAST.start.tree(tree, fossils, complete = FALSE)
}
\arguments{
\item{tree}{an object of class "phylo", representing the tree upon which the
fossil occurrences were simulated.}

\item{fossils}{an object of class "fossils" that corresponds to fossil
occurrences for the "tree" argument.}

\item{complete}{logical, if TRUE then the tree are built for the complete
tree, if FALSE then the tree is built for the crown clades only.}
}
\value{
a string representing the starting tree in newick format.
}
\description{
If complete = FALSE, only the extant taxa are used to construct the tree,
resulting in a DPPDIV style analysis in which the extant topology is fixed
and fossils can float in the tree.
}
\examples{
t = TreeSim::sim.bd.taxa(10, 1, 0.1, 0.05)[[1]]
f = sim.fossils.poisson(0.1, t, root.edge = FALSE)
fossils.to.BEAST.start.tree(t,f, complete = FALSE)
}
