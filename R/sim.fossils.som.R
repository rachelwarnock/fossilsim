#' Simulate fossils under a Poisson sampling model with rate heterogeneity.
#'
#'This model allows users to simulate fossils on a given phylogentic tree with a desired rate distribution.
#'Users can specifiy a rate distribution from which rates for all the edges in the tree (including the root edge, if available)
#'will be drawn. Thus, each edge is applied a different rate. Varying rates are also applied to two edges with a common ancestor. Fossil numbers are simulated proportional to the branch length of the edge.
#'
#' @param tree Phylo object.
#' @param rd The desired rate distribution from which the rates for individual edges would be
#' sampled.\cr (default = rlnorm(nrow(t$edge)+1,0.1,0.1)).\cr
#' The total number of rates to be sampled is equal to the number of edges (including the root
#' edge, if available) in the tree.
#' @param root.edge If TRUE include the root edge (default = TRUE).
#' @return A list with an object of class fossils and a named vector of rates applied on each edge:\cr \cr
#' fossils: sp = node labels and h = ages. \cr \cr
#' rates = rates. A named vector with the total number of rates equal to the sum of number of edges in the tree and root
#' edge (if root.edge = F or the tree does not have root edge, the length equals the number of edges). The names in this vector correspond
#' in a row-wise manner to the edge table of the input phylogentic tree (rates["1"] is the 1st row in the edge table). Thus, the first element of this vector is the rate
#' applied to the first edge of the tree (1st row of the edge table). The last element of this vector is the rate applied to the
#' root edge, if available.
#' @examples
#' # simulate tree
#' t<-TreeSim::sim.bd.taxa(4,1,1,0.10)[[1]]
#'
#' # simulate fossils
#' f<-sim.fossils.poisson.rh(t)
#'
#' #plot fossils on tree
#' plot(f$fossils, t)
#'
#' #rates vector
#' f$rates
#' @export
sim.fossils.poisson.rh <- function(tree, rd = rlnorm(nrow(tree$edge)+1,0.1,0.1), root.edge=TRUE){

  rates <- numeric() #keeps track of the rates applied to the edges of the tree. Also returned as an output object

  node.ages <- n.ages(tree)

  fossils <- data.frame(h=numeric(),sp=numeric())

  root <- length(tree$tip.label) + 1

  if(root.edge && exists("root.edge",tree) ){
    if(nrow(tree$edge)+1 == length(rd)){ #if the total number of rates to be drawn is equal to number of all edges and root edge
      lineages <- c(tree$edge[,2], root)
      rates <- c(rd) #rd is user specified rate distribution
      names(rates) <- c(1:nrow(tree$edge),"RootEdge") #give names to elements of the rates
    }
    else stop("Please Check if the number of rates to be drawn is equal to Number of Edges in the tree (Include the root edge)")
    }
  else{
    lineages <- tree$edge[,2]
    rates <- c(rd)[1:length(lineages)] #take rates only for the edges in the tree (no root edge, no extra rates)
    names(rates) <- c(1:nrow(tree$edge))
  }

  for (i in lineages){ # internal nodes + tips

    if(i == root){

      row <- which(lineages == i)
      # root age
      a <- which(names(node.ages) == root)
      lineage.end <- node.ages[[a]]

      # origin time
      b <- tree$root.edge
      lineage.start <- lineage.end + b

    } else {

      # work out the max age of the lineage (e.g. when that lineage became extant)
      # & get ancestor
      row <- which(tree$edge[,2] == i)
      ancestor <- tree$edge[,1][row]

      # get the age of the ancestor
      a <- which(names(node.ages) == ancestor)
      lineage.start <- node.ages[[a]]

      # work out the min age of the lineage (e.g. when that lineage became extinct)
      # & get the branch length

      b <- tree$edge.length[row]
      lineage.end <- lineage.start - b # branch length
    }

    # sample fossil numbers from the Poisson distribution
    rand <- rpois(1,b*rates[row]) #applying individual rates to each edge b

    if(rand > 0){
      h <- runif(rand, min = lineage.end, max = lineage.start)
      fossils <- rbind(fossils, data.frame(h = h,sp = i))
    }
  }

  fossils <- fossils(fossils, age = "continuous", speciation.mode = "symmetric")
  return(list(fossils = fossils, rates = rates)) # in this data frame h=fossil age and sp=lineage
  # EOF
}

#' Simulate fossils under a Poisson sampling model with rate heterogeneity and dependency.
#'
#'This models allows the user to simulate fossils on a given phylogentic tree with the desired rate distribution.
#'Users can also specifiy different rate distributions for root edge and all other edges.The rates applied
#'on all the edges (except the root) is obtained as a multiplicative factor of the rate applied on the
#'ancesteral of the current edge and the rate distribution applied on the current edge \cr
#'rate.current.edge = rate.ancester.edge * edged \cr
#'The rates applied on 2 edges with a common ancestor will be the same. Fossil numbers are simulated proportial to branch length of the edge.
#'
#'
#' @param tree Phylo object.
#' @param rootd The desired rate distribution from which a single rate for the root edge would be
#' sampled (default = rlnorm(1,0.1,0.1)).
#' @param edged The desired rate distribution from which rates for all the other edges in the tree would be
#' sampled (default = rlnorm(1,0.1,0.1)).
#' @param root.edge If TRUE include the root edge (default = TRUE).
#' @return A list with an object of class fossils and a named vector of rates applied on each edge:\cr \cr
#' fossils: sp = node labels and h = ages. \cr \cr
#' rates = rates. A named vector with the total number of rates is the sum of number of edges in the tree and root
#' edge. The first rate is applied on the root edge, while the rest of rates are applied on all other edges. The names in this vector
#' correspond in a row-wise manner to the edge table of the input phylogentic tree (rates["1"] is the 1st row in the edge table).
#' @examples
#' # simulate tree
#' t<-TreeSim::sim.bd.taxa(4,1,1,0.10)[[1]]
#'
#' # simulate fossils
#' f<-sim.fossils.poisson.rde(t)
#'
#' #plot fossils on tree
#' plot(f$fossils, t)
#'
#' #the rates vector
#' f$rates
#' @export

sim.fossils.poisson.rde <- function(tree, rootd = rlnorm(1,0.1,0.1), edged = rlnorm(1,0.1,0.1), root.edge = TRUE){

  if(length(rootd) > 1 | length(edged) > 1)
    stop("Please provide a rate distribution for a single value (n = 1). Example: rlnorm(1,0.1,0.1)")

  else{

  rates <- c() #keeps track of the rates applied to the edges of the tree. Also returned as an output object.
  ancestor <- c() #internal vector which keeps track of the ancestor nodes on which rates have already been applied. The vector grows as the tree is traversed from root to present time.
  fossils <- data.frame(h = numeric(), sp = numeric())

  node.ages <- n.ages(tree)

  # Check if root edge is present and apply poisson sampling on it.
  if(exists("root.edge",tree)){

    root <- length(tree$tip.label) + 1
    a <- which(names(node.ages) == root)
    lineage.end <- node.ages[[a]]

    b <- tree$root.edge #b is the branch length
    lineage.start <- lineage.end + b

    # sample fossil numbers from the Poisson distribution
    rates <- c(rates, rootd) #rootd is the user-defined rate distribution for the root edge
    names(rates) <- c("RootEdge")#gives names to the elements of the rate vector. The first element is the root edge.
    rand = rpois(1, b*rates) #apply the above rate to the root edge

    if(rand > 0){
      h <- runif(rand, min=lineage.end, max=lineage.start)
      fossils <- rbind(fossils, data.frame(h=h,sp=root))

    }
  }

  else stop("Please input a tree with a root edge.") #This works only for trees with root edge. When the input tree does not have any root edge, the function stops.

  for(i in tree$edge[,2]){

    child1 <- i #child1 is the node currently being investigated by the function.
    row <- which(tree$edge[,2] == child1)
    ancestor1 <- tree$edge[,1][row] #ancestor1 is the ancestor of the current node.

    #Check if the current ancestor is already present in the internal ancestor vector.
    if(ancestor1 %in% ancestor){
      ancestor <- c(ancestor, ancestor1) #If yes, add the current ancestor to the vector
      rates <- c(rates, rates[which(ancestor1 == ancestor) + 1]) #use the rate already applied on that ancestor
      names(rates)[row + 1] <- row
      rates <- rates[!is.na(rates)]
      next #move to the next edge.
    }

    else {

      row2 <- which(tree$edge[,1] == ancestor1)[2]
      child2  <- tree$edge[,2][row2]

      ancestor <- c(ancestor, ancestor1) #add the new ancestor to the internal vector

      # get the age of the ancestor

      a <- which(names(node.ages) == ancestor1)
      lineage.start = node.ages[[a]]

      # work out the min age of the lineage (e.g. when that lineage became extinct)
      # & get the branch length

      b = c(tree$edge.length[row], tree$edge.length[row2]) # The 2 edges that arise from the current ancestor.
      lineage.end = lineage.start - b # branch length

      rates <- c(rates, rates[row]*edged) #edged is the user defined rate distribution for all edges except the root. This rate is a multiplicative factor of the ancesteral rate.

      names(rates)[row + 1] <- row #gives names to the elements of the rate vector. The names from 2nd to the last element of the rates vector correspond to the edges of the tree.

      # sample fossil numbers for each edge from the Poisson distribution

      rand1 <- rpois(1,b[1]*rates[row + 1])
      rand2 <- rpois(1,b[2]*rates[row + 1])

      if(rand1 > 0){
        h = runif(rand1, min = lineage.end[1], max = lineage.start)
        fossils <- rbind(fossils,data.frame(h=h,sp=child1))
      }

      if(rand2 > 0){
        h = runif(rand2, min = lineage.end[2], max = lineage.start)
        fossils <- rbind(fossils,data.frame(h=h,sp=child2))
      }
    }
  }

  fossils <- fossils(fossils, age = "continuous", speciation.mode = "symmetric")
  return(list(fossils = fossils, rates = rates)) # in this data frame h=fossil age and sp=lineage
  # EOF

  }
}
