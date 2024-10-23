sim.test <- function(tree, k = 2, trait.num = 2, ancestral = FALSE){

  ### for a symmetric Q matrix (default?)
  Q <- symmetric_Q_matrix(k)
  trait.num <- trait.num

  states <- as.character(c(0:(k-1)))

  # reorder branches in the phylogeny
  #tree.ordered <- ape::reorder.phylo(tree, "postorder")
  tree.ordered <- reorder(tree)
  edge <- tree.ordered$edge
  num.nodes <- max(edge)
  num.tips <- ape::Ntip(tree.ordered)

  # create an empty matrix
  output <- matrix(NA, num.nodes, trait.num)

  # define the node labels
  parent <- as.integer(edge[, 1])
  child <- as.integer(edge[, 2])

  # identify the root
  root <- as.integer(parent[!match(parent, child, 0)][1])

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

  # remove ancestral sequences
  if (!ancestral) output <- output[tree.ordered$tip.label, , drop = FALSE]

  # formatting for morpho object
   tip_sequences <- rownames(output)
   sequence = list()

  #  create list of simulated traits
   for ( i in 1:length(tip_sequences)){
     sequence[[tip_sequences[i]]] <- output[tip_sequences[i],]
   }

  # create morpho object
  sim.output <- as.morpho(sequence, tree.ordered, "Mk")


  return(sim.output)
}



set.seed(123)
phy <- ape::rtree(15)
plot(phy)
simulated_morpho <- sim.test(phy, k = 5, trait.num = 10)
simulated_morpho

