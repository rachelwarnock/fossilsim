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
#' f<-sim.fossils.poisson(t)
#' plot(f$fossils, t)
#' f$rates
#' @export
sim.fossils.poisson.rh <- function(tree, rd = rlnorm(nrow(tree$edge)+1,0.1,0.1), root.edge=TRUE){

  rates <- numeric()

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
    rates <- c(rates, rd)
    rates <- rates[1:(nrow(tree$edge)+1)]

    rand <- rpois(1,b*rates[row])

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

  rates <- c()
  ancestor <- c()
  fossils <- data.frame(h = numeric(), sp = numeric())

  node.ages <- n.ages(tree)

  if(exists("root.edge",tree)){

    root <- length(tree$tip.label) + 1
    a <- which(names(node.ages) == root)
    lineage.end <- node.ages[[a]]

    b <- tree$root.edge #b is the branch length
    lineage.start <- lineage.end + b

    # sample fossil numbers from the Poisson distribution
    rates <- c(rates, rd)
    rand = rpois(1, b*rates)

    #cat("The random poisson number is : ", rand, "\n")

    if(rand > 0){
      h <- runif(rand, min=lineage.end, max=lineage.start)
      #print(h)
      fossils <- rbind(fossils, data.frame(h=h,sp=root))
      #print(fossils)
    }
  }

  else stop("Please input a rooted tree")


  for(i in tree$edge[,2]){

    extant1 <- i
    row <- which(tree$edge[,2] == extant1)
    ancestor1 <- tree$edge[,1][row]
    if(ancestor1 %in% ancestor){
      #cat(ancestor1, "is already present in", ancestor, "\n")
      ancestor <- c(ancestor, ancestor1)
      rates <- c(rates, rates[which(ancestor1 == ancestor)[1] + 1])
      next
    }

    #if(row == row2){

    # rates <- c(rates, rates[row])
    #next()

    #}

    else {

      row2 <- which(tree$edge[,1] == ancestor1)[2]
      extant2  <- tree$edge[,2][row2]

      ancestor <- c(ancestor, ancestor1)

      # get the age of the ancestor

      a <- which(names(node.ages) == ancestor1)
      lineage.start = node.ages[[a]]

      # work out the min age of the lineage (e.g. when that lineage became extinct)
      # & get the branch length

      b = c(tree$edge.length[row], tree$edge.length[row2])
      lineage.end = lineage.start - b # branch length

      rates <- c(rates, rates[row]*rd)

      # sample fossil numbers from the Poisson distribution

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
