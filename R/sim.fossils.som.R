#' Simulate fossils under a Poisson sampling model with individual rates for every edge
#'
#' @param tree Phylo object.
#' @param rd The desired rate distribution from which the rates for individual edges would be
#' sampled (default = rlnorm(nrow(t$edge)+1,0.1,0.1)).
#' The total number of rates to be sampled is equal to the number of edges (include the root
#' edge) in the tree.
#' @param root.edge If TRUE include the root edge (default = TRUE).
#' @return A list with an object of class fossils and a vector of rates applied on each edge.
#' sp = node labels. h = ages. rates = rates.
#' @examples
#' # simulate tree
#' t<-TreeSim::sim.bd.taxa(4,1,1,0.10)[[1]]
#' # simulate fossils
#' f<-sim.fossils.poisson.rh(t)
#' plot(f$fossils, t)
#' f$rates
#' @export
sim.fossils.poisson.rh <- function(tree, rd = rlnorm(nrow(tree$edge)+1,0.1,0.1), root.edge=TRUE){

  rates <- numeric() #keeps track of the rates applied to the edges of the tree. Also returned as an output object

  node.ages <- n.ages(tree)

  fossils <- data.frame(h=numeric(),sp=numeric())

  root <- length(tree$tip.label) + 1

  if(root.edge && exists("root.edge",tree) ){

    lineages <- c(tree$edge[,2], root)

  } else lineages <- tree$edge[,2]

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
    rates <- c(rates, rd) #rd is the user-defined rate distribution
    rates <- rates[1:(nrow(tree$edge)+1)] #Limits the size of the vector to total number of edges in the tree including the root edge.

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

#' Simulate fossils under a Poisson sampling model rates varying for each edge and dependent on its ancestor
#'
#' @param tree Phylo object.
#' @param rd The desired rate distribution from which the rates for individual edges would be
#' sampled (default = rlnorm(1,0.1,0.1)).
#' @param root.edge If TRUE include the root edge (default = TRUE).
#' @return A list with an object of class fossils and a vector of rates applied on each edge.
#' sp = node labels. h = ages. rates = rates, The total number of rates is the sum of number of edges in the tree and root
#' edge
#' @examples
#' # simulate tree
#' t<-TreeSim::sim.bd.taxa(4,1,1,0.10)[[1]]
#' # simulate fossils
#' f<-sim.fossils.poisson.rde(t)
#' plot(f$fossils, t)
#' f$rates
#' @export

sim.fossils.poisson.rde <- function(tree, rd = rlnorm(1,0.1,0.1), root.edge = TRUE){

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
    rates <- c(rates, rd) #rd is the user-defined rate distribution
    rand = rpois(1, b*rates) #apply the above rate to the root edge

    if(rand > 0){
      h <- runif(rand, min=lineage.end, max=lineage.start)
      fossils <- rbind(fossils, data.frame(h=h,sp=root))

    }
  }

  else stop("Please input a rooted tree") #This works only for rooted trees. When the input tree does not have any roots, the function stops.

  for(i in tree$edge[,2]){

    extant1 <- i #extant 1 is the node currently being investigated by the function.
    row <- which(tree$edge[,2] == extant1)
    ancestor1 <- tree$edge[,1][row] #ancestor1 is the ancestor of the current node.

    #Check if the current ancestor is already present in the internal ancestor vector.
    if(ancestor1 %in% ancestor){
      ancestor <- c(ancestor, ancestor1) #If yes, add the current ancestor to the vector
      rates <- c(rates, rates[which(ancestor1 == ancestor)[1] + 1]) #use the rate already applied on that ancestor
      next #move to the next edge.
    }

    else {

      row2 <- which(tree$edge[,1] == ancestor1)[2]
      extant2  <- tree$edge[,2][row2]

      ancestor <- c(ancestor, ancestor1) #add the new ancestor to the internal vector

      # get the age of the ancestor

      a <- which(names(node.ages) == ancestor1)
      lineage.start = node.ages[[a]]

      # work out the min age of the lineage (e.g. when that lineage became extinct)
      # & get the branch length

      b = c(tree$edge.length[row], tree$edge.length[row2]) # The 2 edges that arise from the current ancestor.
      lineage.end = lineage.start - b # branch length

      rates <- c(rates, rates[row]*rd) #rd is the user defined rate distribution. This rate is a multiplicative factor of the ancesteral rate.

      # sample fossil numbers for each edge from the Poisson distribution

      rand1 <- rpois(1,b[1]*rates[row + 1])
      rand2 <- rpois(1,b[2]*rates[row + 1])

      if(rand1 > 0){
        h = runif(rand1, min = lineage.end[1], max = lineage.start)
        fossils <- rbind(fossils,data.frame(h=h,sp=extant1))
      }

      if(rand2 > 0){
        h = runif(rand2, min = lineage.end[2], max = lineage.start)
        fossils <- rbind(fossils,data.frame(h=h,sp=extant2))
      }
    }
  }
  fossils <- fossils(fossils, age = "continuous", speciation.mode = "symmetric")
  return(list(fossils = fossils, rates = rates)) # in this data frame h=fossil age and sp=lineage
  # EOF
}
