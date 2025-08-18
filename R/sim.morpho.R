#' Simulate characters along branches in a tree
#'
#' @description
#' This function simulates discrete character data along the branches of a phylogentic tree. It can be used with
#' either a time tree or a tree with branch lengths in evolutionary distance. If using a time tree
#' branch rates can be specified, either as one values for all branches or as a vector with
#' different rates per branch. If no branch rates are specified a default of 0.1 is applied to
#' all branches.
#' @param tree A phylogenetic tree (class "phylo") with branches representing genetic distance.
#' @param time.tree A phylogenetic tree (class "phylo") with branches representing time.
#' @param br.rates Clock rates per branch. Can be a single value (strict clock) or a vector of rates.
#' @param k Number of trait states (integer ≥ 2). Can be a vector if using partitions.
#' @param trait.num The total number of traits to simulate (integer > 0).
#' @param partition Vector specifying the number of traits per partition.
#' @param ACRV Among-character rate variation. Default is `NULL`. If supplied, must be `"gamma"`.
#' @param variable If `TRUE`, simulate only varying characters. Default is `FALSE`.
#' @param ancestral If `TRUE`, return the states at all ancestral nodes. Default is `FALSE`.
#' @param partition Specify the number of traits per partition
#' @param fossil Fossil object (from `FossilSim`) to simulate morphological characters.
#' @param define_gamma_rates Vector of gamma rate categories for the simulation.
#' @param alpha.gamma Shape parameter α for the gamma distribution.
#' @param ncats.gamma Number of gamma rate categories.
#' @param define_Q Q matrix for simulation. Must be square and rows must sum to zero.
#'
#' @return An object of class morpho.
#'
#' @export

#' @examples
#'# simulated tree
#' phy <- ape::rtree(10)
#'
#'# simulate characters along the branches of the tree
#'transition_history <-  sim.morpho(tree = phy,
#'                                          k = c(2,3,4),
#'                                          trait.num = 20,
#'                                          ancestral = TRUE,
#'                                          partition = c(10,5,5),
#'                                          ACRV = "gamma",
#'                                          variable = TRUE,
#'                                          ncats.gamma = 4)
#'
#'
#' # To simulate ordered characters:
#' # First define a Q-matrix. The following is for ordered characters where transitions can only occur
#' # between states 0 and 1 and 1 and 2
#'
#' ord_Q <- matrix(c(
#' -0.5, 0.5, 0.0,
#' 0.3333333, -0.6666667, 0.3333333,
#' 0.0, 0.5, -0.5
#' ), nrow = 3, byrow = TRUE)
#'
#' # This Q matrix can be then used to simualte character data.
#'
#' transition_history_2 <-  sim.morpho(tree = phy,
#'                                     k = 3,
#'                                     trait.num = 20,
#'                                     ancestral = TRUE,
#'                                     ACRV = "gamma",
#'                                     variable = TRUE,
#'                                     ncats.gamma = 4,
#'                                    define_Q = ord_Q)



sim.morpho <- function(tree = NULL, time.tree= NULL, ACRV = NULL, br.rates = NULL,  variable = FALSE, ancestral = FALSE,
                               k = 2, partition = NULL, trait.num = NULL, fossil = NULL,
                               alpha.gamma = 1, ncats.gamma = 4, define_gamma_rates = NULL, define_Q = NULL){


  if (is.null(tree) && is.null(time.tree))
    stop ("Must provide a tree object")

  if (any(k <2)) stop("Need to simulate at least 2 states")

  if (is.null(trait.num)) stop ("Specify the total number of traits to partition")

  if(!is.null(ACRV) && ACRV != "gamma")
    stop("Rate variation can only be modeled using a gamma distribution" )

  if(is.null(partition) && length(k) != 1)
    stop("Data being simulated under 1 partition, supply 1 character state for k")

  if(!is.null(partition) && length(k) != length(partition))
    stop("Need to specify a Q matrix size for each partition")

  if(!is.null(partition) &&  sum(partition) != trait.num)
    stop("The total number characters in partition and trait.num must match")

  if (!is.null(define_Q)) {
    if (!is.matrix(define_Q)) stop("`define_Q` must be a matrix.")
    if (nrow(define_Q) != ncol(define_Q)) stop("`define_Q` must be square.")
    if (any(abs(rowSums(define_Q)) > 1e-6)) {
      stop("Incorrect Q matrix specified: rows must sum to zero.")
    }
  }

  if (!is.logical(variable)) stop("`variable` must be TRUE or FALSE.")
  if (!is.logical(ancestral)) stop("`ancestral` must be TRUE or FALSE.")
  if (!is.null(fossil) && !is.data.frame(fossil)) {
    stop("`fossil` must be a FossilSim object.")
  }


  ## if provided with time tree, need to transform branches in genetic distance
  ## rates can be a single value or a vector for each branch
  if (is.null(tree) && !is.null(time.tree)){
    tree <- time.tree

    if(is.null(br.rates)){
      print("No branch rate provide, using default of 0.1 for all branches")
      br.rates <- rep(0.1, length(tree$edge.length))
    }

      if(length(br.rates) == 1){
        br.rates <- rep(br.rates, length(tree$edge.length))
      }

      tree$edge.length <- time.tree$edge.length * br.rates
  }


  tree.ordered <- reorder(tree)

  ##reorder the nodes on the time tree to match the format of the genetic distance tree
  if (!is.null(time.tree)) {
    time.tree.order <- reorder(time.tree)
  } else {time.tree.order = NULL}

  # counter for trait number
  tr.num <- 1
  edge <- tree.ordered$edge
  num.nodes <- max(edge)
  num.tips <- ape::Ntip(tree.ordered)
  rs <- c()


  # define the node labels
  parent <- as.integer(edge[, 1])
  child <- as.integer(edge[, 2])
  tips <- setdiff(child,parent)
  nodes <- setdiff(parent,tips)


  # create an result matrix
  transition_history <- list()
  state_at_tips <- matrix(nrow = num.tips, ncol = trait.num)
  rownames(state_at_tips) <- tree.ordered$tip.label
  state_at_nodes <- matrix(nrow = tree$Nnode, ncol = trait.num)
  rownames(state_at_nodes) <- nodes
  ACRV_rate <- matrix(ncol = trait.num, nrow = 1)

  # identify the root
  root<- as.integer(parent[!match(parent, child, 0)][1])


  ## current loop is a bit weird and depends on there being partitions.
  ## if traits set using one partition need to set partition <- trait.num

  if(is.null(partition)) partition <- trait.num

  ### Among character rate variation
  # eventually may need to put into the state function if partitions are not linked
  if(!is.null(ACRV)){
    if (ACRV == "gamma"){
      if (!is.null(define_gamma_rates)){
        gamma_rates <- define_gamma_rates
      } else {
        gamma_rates <- get_gamma_rates(alpha.gamma, ncats.gamma)
      }
    }
    else if(ACRV == "lognormal"){

    }
  }

  ## start the partition loop
  for (part in 1:length(partition)){
    part.trait.num <- partition[part]

    ##start loop for number of different states
    # for (stat in 1:length(k))
    ### for a symmetric Q matrix
    part_k <- k[part]

    if (is.null(define_Q)){
      Q <- symmetric.Q.matrix(part_k)
    } else {Q <- define_Q}

    states <- as.character(c(0:(part_k -1)))
    bl <- tree.ordered$edge.length


    for (tr in 1:part.trait.num){
      # message(tr)
      repeat {

        # simulate the root state using symmetric Q-matrix
        root.state <- sample(states, 1, replace = TRUE, prob = rep(1/part_k , part_k)) #later prob would need to take into account for bf
        state_at_nodes[as.character(root),tr.num] <- as.numeric(root.state)

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

          total_wait <- 0
          from <- parent[i]
          to <- child[i]

          current_state <- as.numeric(unname(state_at_nodes[as.character(from), tr.num]))

          # Time tracking
          time_remaining <- bl[i]
            while (time_remaining > 0) {
            # Get the rates for the current state
            rates <- Q_r[current_state+1, ]
            rates[current_state+1] <- 0 # No transition to the same state

             # calculate the rate of any transition occurring
            total_rate <- -Q_r[current_state+1, current_state+1]
            waiting_time <- rexp(1, rate = total_rate)

            # Transition to a new state
            time_remaining <- time_remaining - waiting_time
            if (time_remaining <= 0) {
              break
            }

            # Transition to a new state
            next_state <- sample(states, 1, prob = rates / total_rate)
            current_state <- as.numeric(next_state)
            add_t <- c(i, next_state, waiting_time + total_wait)
            transitions <- rbind(transitions, as.numeric(add_t))
            total_wait <- total_wait + waiting_time

            }

          if (to %in% tips) state_at_tips[tree$tip.label[to],tr.num] <- current_state
          if (to %in% nodes) state_at_nodes[as.character(to),tr.num] <- current_state
        }

        if (!variable ){
          break
        }

        if (length(unique(state_at_tips[,tr.num])) > 1 & length(unique(state_at_nodes[,tr.num])) > 1 ){
          break
        }
      }

      if (!is.null(ACRV)) ACRV_rate[tr.num] <- which(gamma_rates == trait_rate)
      transition_history[[tr.num]] <- as.data.frame(transitions)
      tr.num <- tr.num + 1
      rs <- rbind(rs, root.state)
    }

  }

  ### fossil object
   if(!is.null(fossil)){
    f.morpho <- fossil[fossil$hmin != 0, ]
    f.morpho$ape.branch <-NA
    f.morpho$specimen <- NA
    f.morpho$specimen <- seq(1,length(f.morpho$sp))

    state_at_fossils <- matrix(nrow = length(f.morpho$sp), ncol = trait.num)
    rownames(state_at_fossils) <- f.morpho$specimen


    for ( i in 1:length(f.morpho$edge)){
      child <-  f.morpho$edge[i]
      f.morpho$ape.branch[i] <- which(tree.ordered$edge[,2] == child)
    }


    tree.age <- max(ape::node.depth.edgelength(time.tree))

    for ( tr.num in 1:trait.num){

      ##go through each fossil for trait n
      for (spec in 1:length(f.morpho$specimen)){

        branch <- f.morpho$ape.branch[[spec]]
        parent <- tree.ordered$edge[branch,1]
        current_state <- state_at_nodes[[which(rownames(state_at_nodes) == parent),tr.num]]

        ## does this branch have any changes?

        if (!(f.morpho$ape.branch[spec]  %in% transition_history[[tr.num]][[1]])){
          state_at_fossils[spec,tr.num] <- current_state
        } else {
          time_position_fossil <-  tree.age -   f.morpho$hmin[spec] + time.tree$root.edge
          changes <- which(transition_history[[tr.num]][[1]] == f.morpho$ape.branch[spec])
          changes_along_edge <- transition_history[[tr.num]][changes,]

          ## get the age of the partent node
          node_age <- ape::node.depth.edgelength(time.tree)[time.tree$edge[branch,1]]
          # get the time of change
          changes_along_edge$time <- (changes_along_edge[,3] / br.rates[branch]) +  node_age + time.tree$root.edge


        ## remove changes after fossil occurance
        later <- which(changes_along_edge[,4] > time_position_fossil)
         if (!length(later)== 0) changes_along_edge <- changes_along_edge[-later,]


        if ( length(changes_along_edge[,1]) == 0) {
          state_at_fossils[spec,tr.num] <- current_state
        } else {
          tran <- which(changes_along_edge$hmin == max(changes_along_edge$hmin))
          state_at_fossils[spec,tr.num] <- changes_along_edge[tran,2]
        }
      }
    }

       # formatting for morpho object
    fossil_names <- paste0(f.morpho$specimen,"_", f.morpho$ape.branch)
    fossil_sequence = list()
    rownames(state_at_fossils) <- fossil_names
    #  create list of simulated traits
    for ( i in 1:length(fossil_names)){
      fossil_sequence[[fossil_names[i]]] <- state_at_fossils[fossil_names[i],]
     }
   }
  }

  ## create morpho object

  seq <- list(NA,NA,NA)
  names(seq) <- c("tips","nodes", "SA")
  tip_names <- rownames(state_at_tips)
  tip_sequence = list()

  #  create list of simulated traits
  for ( i in 1:length(tip_names)){
    tip_sequence[[tip_names[i]]] <-state_at_tips[tip_names[i],]
  }

  seq[["tips"]] <- tip_sequence


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

  seq[["nodes"]] <- node_sequence

  if(is.null(fossil)) {
    fossil_sequence <- NULL
    f.morpho <- NULL
  }

  seq[["SA"]] <- fossil_sequence

  if(is.null(ACRV)){
    ACRV_rate <- NULL}

  if(!exists("gamma_rates")){
    gamma_rates <- NULL
  }


  # define model used
  if (isTRUE(variable)) {
    var <- "V"
  } else {
    var <- NULL
  }
  model_components <- paste0("Mk", ACRV , var, "_Part:", partition, "states:", k)

  trees <- list(NA, NA, NA)
  names(trees) <- c("EvolTree", "TimeTree", "BrRates")
  trees[["EvolTree"]] <-  tree.ordered
  trees[["TimeTree"]] <- time.tree
  trees[["BrRates"]] <- br.rates

  model <- list(NA,NA,NA)
  names(model) <- c("Specified", "RateVar", "RateVarTrait")
  model[["Specified"]] <- model_components
  model[["RateVar"]] <- gamma_rates
  model[["RateVarTrait"]] <- ACRV_rate

  sim.output <- as.morpho(data = seq, trees = trees, model = model,
                           transition_history = transition_history,
                          root.states = rs,  fossil = f.morpho)

  return(sim.output )


}




