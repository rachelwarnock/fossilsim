#' Simulate fossils under a Poisson sampling model with rate heterogeneity.
#'
#'This model allows users to simulate fossils on a given phylogenetic tree with a desired rate distribution.
#'Users can specifiy a rate distribution from which rates for all the edges in the tree (including the root edge, if available)
#'will be drawn. Thus, each edge is applied a different rate. Varying rates are also applied to two edges with a common ancestor. Fossil numbers are simulated proportional to the branch length of the edge.
#'
#' @param tree Phylo object.
#' @param rates The desired rate distribution from which the rates for individual edges would be
#' sampled.\cr (default = rlnorm(nrow(t$edge)+1,0.1,0.1)).\cr
#' The total number of rates to be sampled is equal to the number of edges (including the root
#' edge, if available) in the tree.
#' @param root.edge If TRUE include the root edge if present (default = TRUE).
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
sim.fossils.poisson.rh <- function(tree, rates = rlnorm(nrow(tree$edge)+1,0.1,0.1), root.edge=TRUE){
  
  node.ages <- n.ages(tree)
  
  fossils <- data.frame(h=numeric(),sp=numeric())
  
  root <- length(tree$tip.label) + 1
  
  if(root.edge && exists("root.edge",tree) ){
    if(nrow(tree$edge)+1 == length(rd)){ #if the total number of rates to be drawn is equal to number of all edges and root edge
      lineages <- c(tree$edge[,2], root)
      names(rates) <- c(1:nrow(tree$edge),"RootEdge") #give names to elements of the rates
    }
    else stop("Please Check if the number of rates to be drawn is equal to Number of Edges in the tree (Include the root edge)")
  }
  else{
    lineages <- tree$edge[,2]
    rates <- rates[1:length(lineages)] #take rates only for the edges in the tree (no root edge, no extra rates)
    names(rates) <- c(1:nrow(tree$edge))
  }
  
  for (i in lineages){ # internal nodes + tips
    
    if(i == root){
      
      row <- length(lineages)
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
}

#' Simulate fossils under a Poisson sampling model with rate heterogeneity and dependency.
#'
#'This models allows the user to simulate fossils on a given phylogentic tree with the desired rate distribution.
#'Users can also specifiy different rate distributions for root edge and all other edges.The rates applied
#'on all the edges (except the root) is obtained as a multiplicative factor of the rate applied on the
#'ancestral of the current edge and the rate distribution applied on the current edge \cr
#'rate.current.edge = rate.ancestor.edge * ratio_d \cr
#'The rates applied on 2 edges with a common ancestor will be the same. Fossil numbers are simulated proportial to branch length of the edge.
#'
#'
#' @param tree Phylo object.
#' @param origin_d The rate distribution from which a single rate for the root edge would be
#' sampled (default = function() {rlnorm(1,0.1,0.1)}).
#' @param ratio_d The rate distribution from which ratios for non-root edges will be sampled 
#' (default = function() {rlnorm(1,0.1,0.1)}).
#' @param root.edge If TRUE include the root edge in the fossil sampling if present (default = TRUE).
#' @return A list with an object of class fossils and a named vector of rates applied on each edge:\cr \cr
#' fossils: sp = node labels and h = ages. \cr \cr
#' rates = rates. A named vector with the total number of rates is the sum of number of edges in the tree and root
#' edge. The first rate is applied on the root edge if present and root.edge=T, otherwise it just serves as basis for the other rates. 
#' The rest of the rates are applied on all other edges. The names in this vector
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

sim.fossils.poisson.rde <- function(tree, origin_d = function() {rlnorm(1,0.1,0.1)}, ratio_d = function() {rlnorm(1,0.1,0.1)}, root.edge = TRUE){
  
  if(length(ratio_d()) > 1 || length(origin_d()))
    stop("Please provide a rate distribution for a single value (n = 1). Example: rlnorm(1,0.1,0.1)")
  
  rates <- origin_d() #keeps track of the rates applied to the edges of the tree. Also returned as an output object.
  names(rates) <- "RootEdge" #gives names to the elements of the rate vector. The first element is the root edge.
  fossils <- data.frame(h = numeric(), sp = numeric())  # in this data frame h=fossil age and sp=lineage
  
  node.ages <- n.ages(tree)
  root <- length(tree$tip.label) + 1
  
  # Check if root edge is present and apply poisson sampling on it.
  if(root.edge && exists("root.edge",tree)) {
    a <- which(names(node.ages) == root)
    lineage.end <- node.ages[[a]]
    
    b <- tree$root.edge #b is the branch length
    lineage.start <- lineage.end + b
    
    # sample fossil numbers from the Poisson distribution
    rand = rpois(1, b*rates) #apply the above rate to the root edge
    
    if(rand > 0){
      h <- runif(rand, min=lineage.end, max=lineage.start)
      fossils <- rbind(fossils, data.frame(h=h,sp=root))
      
    }
  }
  
  aux = function(ancestor, row, current_rate, list_fr) {
    a <- which(names(node.ages) == ancestor)
    lineage.start = node.ages[[a]]
    
    node = tree$edge[row,2]
    lineage.end = lineage.start - tree$edge.length[row] # branch length
    
    new_rate = current_rate * ratio_d() #ratio_d is the user defined rate distribution for all edges except the root. This rate is a multiplicative factor of the ancesteral rate.
    list_fr$rates <- c(rates, new_rate) 
    names(list_fr$rates) = c(names(list_fr$rates), row) #gives names to the elements of the rate vector. The names from 2nd to the last element of the rates vector correspond to the edges of the tree.
    
    # sample fossil numbers for each edge from the Poisson distribution
    rand1 <- rpois(1,tree$edge.length[row]*new_rate)
    
    if(rand1 > 0){
      h = runif(rand1, min = lineage.end[1], max = lineage.start)
      list_fr$fossils <- rbind(list_fr$fossils,data.frame(h=h,sp=node))
    }
    
    desc_edges = which(tree$edge[,1] == node)
    for(child in desc_edges) {
      list_fr = aux(node, child, new_rate, list_fr)
    }
    return(list_fr)
  }
  
  list_fr = list(fossils = fossils, rates = rates)
  desc_edges = which(tree$edge[,1] == root)
  for(child in desc_edges) {
    list_fr = aux(root, child, rates, list_fr)
  }
  
  list_fr$fossils <- fossils(list_fr$fossils, age = "continuous", speciation.mode = "symmetric")
  return(list_fr)
}
