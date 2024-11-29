#' This function is to simulate discrete character data at the tips and internal nodes
#' @param data Trait data associated with the tips of the tree
#' @param tree Tree with branches that represent genetic distance associated with the character data.
#' @param time.tree Tree with branches that represent time associated with the character data.
#' @param br.rate clock Rate, currently can only be strict clock (a single rate)
#' @param k Number of states
#' @param trait.num The number of traits to simulate
#' @param ancestral T/F Whether or not to output information from internal nodes
#' @export


sim.morpho <- function(tree = NULL, time.tree= NULL, br.rates = NULL,
                       k = 2, trait.num = 2, ancestral = FALSE){

  # check that a tree is provided
  if (is.null(tree) && is.null(time.tree)) stop(cat("Must provide a tree object"))
  # check that there are more then two states
  if (k < 2) stop(cat("Trait data must have more than 1 state"))



  ## if provided with time tree, need to transform branches in genetic distance
  ## rates can be a single value or a vector for each branch
  if (is.null(tree) && !is.null(time.tree)){
    tree <- time.tree

    if(is.null(br.rates)){
      print("No branch rate provide, using default of 0.1 for all branches")
      tree$edge.length <- time.tree$edge.length * 0.1
    } else {
    tree$edge.length <- time.tree$edge.length * br.rates
    }
  }



  ### for a symmetric Q matrix (default?)
  Q <- symmetric.Q.matrix(k)
  trait.num <- trait.num

  states <- as.character(c(0:(k-1)))

  # reorder branches in the phylogeny
  #tree.ordered <- ape::reorder.phylo(tree, "postorder")
  tree.ordered <- reorder(tree)

  ##reorder the nodes on the time tree to match the format of the genetic distance tree
  if (!is.null(time.tree)) {
    time.tree.order <- reorder(time.tree)
  } else {time.tree.order = NULL}

  edge <- tree.ordered$edge
  num.nodes <- max(edge)
  num.tips <- ape::Ntip(tree.ordered)

  # create an empty matrix
  output <- matrix(NA, num.nodes, trait.num)

  # define the node labels
  parent <- as.integer(edge[, 1])
  child <- as.integer(edge[, 2])

  # identify the root
  root<- as.integer(parent[!match(parent, child, 0)][1])

  # simulate the root state
  root.state <- sample(states, trait.num, replace = TRUE, prob = rep(1/k, k)) #later prob would need to take into account for bf

  # add the root state to the matrix
  output[root, ] <- root.state

  # edge lengths
  bl <- tree.ordered$edge.length

  # loop through all branches
  for (i in seq_along(bl)){

    from <- parent[i]
    to <- child[i]

    # generate a matrix of probabilities
    #P <- getP(tl[i], eig, rate)[[1]]
    # probability of change along each branch
    #P  <- lapply(bl[i], function(l) ape::matexpo(Q[[1]] * l))
    P  <- ape::matexpo(Q * bl[i])

    # avoid numerical problems for larger P and small bl
    if (any(P < 0)) P[P < 0] <- 0

    for (j in 1:k) {
      # boolean, if current state is TRUE then we calc the prob of changing from that state. Useful loop is sites > 1, otherwise some redundancy
      ind <- output[from, ] == states[j]
      output[to, ind] <- sample(states, sum(ind), replace = TRUE, prob = P[, j])
    }
  }

  tip.labels <- c(tree.ordered$tip.label, as.character( (num.tips + 1):num.nodes))
  rownames(output) <- tip.labels

  # do we want information about the node
  if (ancestral) {
    output_node.seq <- output[setdiff(rownames(output), tree.ordered$tip.label), , drop = FALSE]
    output <- output[tree.ordered$tip.label,  ,drop = FALSE]

  } else {
   output <- output[tree.ordered$tip.label,  ,drop = FALSE]
   output_node.seq <- NULL
}
  # formatting for morpho object
  tip_names <- rownames(output)
  tip_sequence = list()

  #  create list of simulated traits
  for ( i in 1:length(tip_names)){
    tip_sequence[[tip_names[i]]] <- output[tip_names[i],]
  }

  if (ancestral){
    # formatting for morpho object
    node_names <- rownames(output_node.seq)
    node_sequence = list()

    #  create list of simulated traits
    for ( i in 1:length(node_names)){
      node_sequence[[node_names[i]]] <- output_node.seq[node_names[i],]
    }
  } else {
    node_sequence <- NULL
  }
  # create morpho object
  sim.output <- as.morpho(data = tip_sequence, tree = tree.ordered, model = "Mk",
                          node.seq = node_sequence, time.tree = time.tree.order )


  return(sim.output)
}



