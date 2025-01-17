#' Simulate characters along branches in a tree
#'
#' @description
#' This function simulates discrete character data along the branches of a phylogentic tree. It can be used with
#' either a time tree or a tree with branch lengths in genetic distance. If using a time tree
#' branch rates can be specified, either as one values for all branches or as a vector with
#' different rates per branch. If no branch rates are specified a default of 1 is applied to
#' all branches.
#' @param tree Tree with branches that represent genetic distance associated with the character data.
#' @param ACRV Allow for among character rate variation. The default here is set to NULL. Can supply arguments "gamma" or "lognormal"
#' @param variable Simulate only varying characters. The default here is set to FALSE
#' @param time.tree Tree with branches that represent time associated with the character data.
#' @param br.rates Clock rates per branch, currently can only be strict clock (a single rate)
#' @param k Number of trait states.
#' @param trait.num The number of traits to simulate
#' @param ancestral Return the states at all ancestral nodes. Default set to false
#'
#' @return An object of class morpho.
#'
#' @export

#' @examples
#'  # simulated tree
#'  phy <- ape::rtree(10)
#'
#'  # simulate characters along the branches of the tree
#'  transition_history <-  sim.morpho.history(tree = phy, k = 3, trait.num = 4, ancestral = TRUE)
#'


sim.morpho.history <- function(tree = NULL, time.tree= NULL, ACRV = NULL, br.rates = NULL,  variable = FALSE, ancestral = FALSE,
                                       k = 2, trait.num = 2, alpha.gamma = 1, ncats.gamma = 4){


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
  #trait.num <- trait.num

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

  # define the node labels
  parent <- as.integer(edge[, 1])
  child <- as.integer(edge[, 2])
  tips <- setdiff(child,parent)
  nodes <- setdiff(parent,tips)


  # create an result matrix
  continuous_traits <- list()
  state_at_tips <- matrix(nrow = num.tips, ncol = trait.num)
  rownames(state_at_tips) <- tree.ordered$tip.label
  state_at_nodes <- matrix(nrow = tree$Nnode, ncol = trait.num)
  rownames(state_at_nodes) <- nodes
  ACRV_rate <- matrix(ncol = trait.num, nrow = 1)

  # identify the root
  root<- as.integer(parent[!match(parent, child, 0)][1])

  # simulate the root state using symmetric Q-matrix
  root.state <- sample(states, trait.num, replace = TRUE, prob = rep(1/k, k)) #later prob would need to take into account for bf
  state_at_nodes[as.character(root),] <- as.numeric(root.state)

  # edge lengths
  bl <- tree.ordered$edge.length

  ### Among character rate variation
  if(!is.null(ACRV)){
    if (ACRV == "gamma"){
      gamma_rates <- get_gamma_rates(alpha.gamma, ncats.gamma)
    }
    else if(ACRV == "lognormal"){

    }
  }

   for (tr in 1:trait.num){
    # message(tr)
    repeat {
      #message("conditions not met repreat")
    ### Among character rate variation
    if(!is.null(ACRV)){
      if (ACRV == "gamma"){
        trait_rate <- gamma_rates[sample(1:ncats.gamma, 1)]
        Q_r <- Q * trait_rate
      }
      else if(ACRV == "lognormal"){

      }
    } else {
      Q_r <- Q
    }

      # container for the transitions
  transitions <- matrix(ncol = 3, nrow= 0)
  colnames(transitions) <- c("edge", "state", "hmin")
  # loop through all branches
  for (i in seq_along(bl)){

    from <- parent[i]
    to <- child[i]

      current_state <- as.numeric(unname(state_at_nodes[as.character(from), tr]))

    # best way to specify these?
   # if (Qmat == "Equal"){}
   # else if (Qmat == "F81"){}

    # calculate the rate of any transition occurring
    current_rate <- Q_r[(current_state+ 1),(current_state+1)] * -1

    # how many transitions
    rand = rpois(1, bl[i]*current_rate)

    # if there are any trnaisitons along this branch
    if (rand > 0 ){
    h = runif(rand, min = 0, max = bl[i])
    h <- sort(h)

    # now we need to calculate what those transitions are
    P  <- ape::matexpo(Q_r * bl[i])

    # avoid numerical problems for larger P and small bl
    if (any(P < 0)) P[P < 0] <- 0
    # calculated the number of changes there for we need to ensure that a new state is sampled.
     state_changes <- c()
    new_state <- current_state
    for (r in 1:rand){
      while(as.numeric(new_state == current_state) ){
        new_state <- sample(states, 1, replace = TRUE, prob = P[,(current_state+1)])
      }
      state_changes <- rbind(state_changes, new_state)
      current_state <- as.numeric(new_state)

    }

    for (q in 1:rand){
      add_t <- c(i, state_changes[q], h[q])
      transitions <- rbind(transitions, as.numeric(add_t))

    }
}

    if (to %in% tips) state_at_tips[tree$tip.label[to],tr] <- current_state
    if (to %in% nodes) state_at_nodes[as.character(to),tr] <- current_state
}

   if (!variable ){
       break
   }

    if (length(unique(state_at_tips[,tr])) > 1 & length(unique(state_at_nodes[,tr])) > 1 ){
     break
    }
    }

     if (!is.null(ACRV)) ACRV_rate[tr] <- trait_rate
      continuous_traits[[tr]] <- as.data.frame(transitions)
  }

  ## create morpho object
  # formatting for morpho object
  tip_names <- rownames(state_at_tips)
  tip_sequence = list()

  #  create list of simulated traits
  for ( i in 1:length(tip_names)){
    tip_sequence[[tip_names[i]]] <-state_at_tips[tip_names[i],]
  }


  if (ancestral){
    # formatting for morpho object
    node_names <- rownames(state_at_nodes)
    node_sequence = list()

    #  create list of simulated traits
    for ( i in 1:length(node_names)){
      node_sequence[[node_names[i]]] <- state_at_nodes[node_names[i],]
    }
  } else {
    node_sequence <- NULL
  }

 if(is.null(ACRV)){
    ACRV_rate <- NULL}

  sim.output <- as.morpho(data = tip_sequence, tree = tree.ordered, model = "Mk",
                          time.tree = time.tree, continuous_traits= continuous_traits,
                          root.states = root.state, node.seq = node_sequence,  ACRV_rate =  ACRV_rate  )

  return(sim.output )


}




