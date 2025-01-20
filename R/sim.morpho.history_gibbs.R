#' Simulate characters along branches in a tree
#'
#' @description
#' This function simulates discrete character data along the branches of a phylogentic tree. It can be used with
#' either a time tree or a tree with branch lengths in genetic distance. If using a time tree
#' branch rates can be specified, either as one values for all branches or as a vector with
#' different rates per branch. If no branch rates are specified a default of 1 is applied to
#' all branches.
#' @param tree Tree with branches that represent genetic distance associated with the character data.
#' @param ACRV Allow for among character rate variation. The default here is set to NULL. Can supply arguments "gamma" or "lognormal".
#' @param variable Simulate only varying characters. The default here is set to FALSE.
#' @param time.tree Tree with branches that represent time associated with the character data.
#' @param br.rates Clock rates per branch, currently can only be strict clock (a single rate).
#' @param k Number of trait states.
#' @param trait.num The number of traits to simulate.
#' @param ancestral Return the states at all ancestral nodes. Default set to FALSE.
#' @param define_gamma_rates You can specify the gamma rate categories you want to use for the simulation.
#'
#' @return An object of class morpho.
#'
#' @export

#' @examples
#'# simulated tree
#' phy <- ape::rtree(10)
#'
#'# simulate characters along the branches of the tree
#'transition_history <-  sim.morpho.history(tree = phy,
#'                                          k = c(2,3,4),
#'                                          trait.num = 20,
#'                                          ancestral = TRUE,
#'                                          partition = c(10,5,5),
#'                                          ACRV = "gamma",
#'                                          variable = TRUE,
#'                                          ncats.gamma = 4)



sim.morpho.history_gibbs <- function(tree = NULL, time.tree= NULL, ACRV = NULL, br.rates = NULL,  variable = FALSE, ancestral = FALSE,
                               k = 2, partition = NULL, trait.num = 2,
                               alpha.gamma = 1, ncats.gamma = 4, define_gamma_rates = NULL, define_Q = NULL){


  # check that a tree is provided
  #if (is.null(tree) && is.null(time.tree)) stop("Must provide a tree object")
  #if (k < 2) stop("Trait data must have more than 1 state")


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

  # reorder branches in the phylogeny
  #tree.ordered <- ape::reorder.phylo(tree, "postorder")
  tree.ordered <- reorder(tree)

  ##reorder the nodes on the time tree to match the format of the genetic distance tree
  if (!is.null(time.tree)) {
    time.tree.order <- reorder(time.tree)
  } else {time.tree.order = NULL}

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
  continuous_traits <- list()
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


    # edge lengths
    bl <- tree.ordered$edge.length


    for (tr in 1:part.trait.num){
      # message(tr)
      repeat {

        # simulate the root state using symmetric Q-matrix
        root.state <- sample(states, 1, replace = TRUE, prob = rep(1/part_k , part_k)) #later prob would need to take into account for bf
        state_at_nodes[as.character(root),tr.num] <- as.numeric(root.state)



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
            #time_remaining <- time_remaining - waiting_time
            next_state <- sample(states, 1, prob = rates / total_rate)
            #states <- c(states, next_state)
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
      continuous_traits[[tr.num]] <- as.data.frame(transitions)
      tr.num <- tr.num + 1
      rs <- rbind(rs, root.state)
    }

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

  if(!exists("gamma_rates")){
    gamma_rates <- NULL
  }

  # define model used
  if (isTRUE(variable)) {
    var <- "V"
  } else {
    var <- NULL
  }
  model <- paste0("Mk", ACRV , var, "_Part:", partition, "states:", k)

  sim.output <- as.morpho(data = tip_sequence, tree = tree.ordered, model = model,
                          time.tree = time.tree, continuous_traits= continuous_traits,
                          root.states = rs, node.seq = node_sequence,  ACRV_rate =  ACRV_rate,
                          gamma_rates = gamma_rates  )

  return(sim.output )


}




