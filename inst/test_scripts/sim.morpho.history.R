#' Simulate characters along branches in a tree
#'
#' @description
#' This function simulates discrete character data along the branches of a phylogentic tree. It can be used with
#' either a time tree or a tree with branch lengths in genetic distance. If using a time tree
#' branch rates can be specified, either as one values for all branches or as a vector with
#' different rates per branch. If no branch rates are specified a default of 1 is applied to
#' all branches.
#' @param tree Tree with branches that represent genetic distance associated with the character data.
#' @param time.tree Tree with branches that represent time associated with the character data.
#' @param br.rates Clock rates per branch, currently can only be strict clock (a single rate)
#' @param k Number of trait states.
#' @param trait.num The number of traits to simulate
#'
#' @return An object of class morpho.
#'
#' @export

#' @examples
#'  # simulated tree
#'  phy <- ape::rtree(10)
#'
#'  # simulate characters along the branches of the tree
#'  continuous_traits <- sim.morpho.completeprocess(phy, k=4, trait.num =12)
#'


sim.morpho.history <- function(tree = NULL, time.tree= NULL, br.rates = NULL,
                                       k = 2, trait.num = 2){





  # check that a tree is provided
  if (is.null(tree) && is.null(time.tree)) stop("Must provide a tree object")
  if (k < 2) stop("Trait data must have more than 1 state")


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



  ### for a symmetric Q matrix
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

  # create an result matrix
  continuous_traits <- list()
  state_at_tips <- matrix(nrow = num.tips, ncol = trait.num)
  rownames(state_at_tips) <- tree.ordered$tip.label

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

    # get the number of transitions along the branch.
    # this will need to be updated for asymmetric transition matrices


    ## cureently calculates the probability of a single transition. I think.
    # given x amout of time what state at the end. Does not allow for multiple changes.
    # is that wrong?
    # TODO: 1. calculate the number of changes according to a poissoin process
    #       2. what are those changes?
    #       3. where along the branch are they found
    #       4. ensure that the design will allow for asymmetric transitions

    # best way to specify these?
    if (Qmat == "Equal"){}
    else if (Qmat == "F81"){}


    current_rate <- Q[(parent_state+ 1),(parent_state+1)] * -1

    # rate needs to change for an asymmetric - current rate at which there is any
    # transition
    rand = rpois(1, bl[i]*current_rate)
    h = runif(rand, min = 0, max = bl[i])
    h <- sort(h)

# now we need to calculate what those transitions are
    P  <- ape::matexpo(Q * bl[i])

    # avoid numerical problems for larger P and small bl
    if (any(P < 0)) P[P < 0] <- 0

    for (j in 1:k) {
      # boolean, if current state is TRUE then we calc the prob of changing from that state. Useful loop is sites > 1, otherwise some redundancy
      ind <- output[from, ] == states[j]
      output[to, ind] <- sample(states, sum(ind), replace = TRUE, prob = P[, j])
    }
  }








}
