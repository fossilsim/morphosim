#' Calculates statistics for a morpho object
#'
#' @description
#' This function computes three key pieces of information:
#' 1. The Consistency Index (CI) and Retention Index (RI) based on
#' the tip sequence data.
#' 2. Convergent traits, identifying traits that have evolved independently
#'  multiple times.
#' 3. Summary information about the size and structure of the tree.
#'
#' @param data Morpho object
#'
#' @export
#'
#' @examples
#' summary <- stats_morpho(data = morpho_data)

stats_morpho <- function(data){

  morpho_summary <- vector("list", 3)
  names(morpho_summary) <- c("Statistics", "Convergent_Traits", "Tree")
  ## consistency index & retention index
  phylev <- unique(unlist(unique(data$sequences$tips)))
  tree <- data$trees$EvolTree
  data_ph <-  phangorn::phyDat(data$sequences$tips, type="USER", levels=phylev, return.index = TRUE)

  morpho_summary[["Statistics"]] <- data.frame(
    CI =  phangorn::CI(tree, data_ph),
    RI = phangorn::RI(tree, data_ph)
  )


  ## convergent characters
  morpho_summary[["Convergent_Traits"]] <-convergent_evol(data = data)

  ## number SA
  if (!is.null(data$sequences$SA)){
  num_sa <- length(data$sequences$SA)
  } else {
    num_sa <- 0
  }
  ## number extant tips
  tip_depths <- node.depth.edgelength(data$trees$TimeTree)[1:length(data$trees$TimeTree$tip.label)]
  tree_height <- max(node.depth.edgelength(data$trees$TimeTree))
  extant_tips <- length(data$trees$TimeTree$tip.label[abs(tip_depths - tree_height) < 1e-8])

  ## number extinct tips
  extinct_tips <- length(data$trees$TimeTree$tip.label) - extant_tips

  morpho_summary[["Tree"]] <- data.frame(
    Extanat =  extant_tips,
    Extinct = extinct_tips,
    Sampled_Ancestors = num_sa
  )

  return(morpho_summary)

}

#' Determines the number of convergently evolved traits
#' @description
#' This function determines which traits have evoled through convergent evolution
#'
#' @param data a morpho object
#'
#'  @export
#'

convergent_evol <- function(data = NULL){

  dat <- data[["transition_history"]]
  tree <- data$trees$EvolTree
  tip_states <- as.data.frame(data$sequences$tips)
  ntax <- length(tree$tip.label)
  convergent_traits <- matrix(ncol = 3)
  colnames(convergent_traits) <- c("trait", "state", "num.transitions")

  for (convT in 1:length(dat)){
    temp <- dat[[convT]]
    trait_n <- tip_states[,convT]


    trans_trait <- matrix(nrow =  length(tree$tip.label), ncol = 2)
    rownames(trans_trait) <- tree$tip.label
    colnames(trans_trait)<- c("transition", "state")
    ## number our transitions
    temp$num <- seq(1, length(temp$edge), 1)

    for(tip in 1:length(tree$tip.label)){
      # go through the tree and get the last transistion associated with each tip
      route_n <- find_path_to_tip(tree, tree$tip.label[tip])
      ## where are the transitions? starting at the route
      for (l in seq(length(route_n[,1]), 1)){
        bran <- which(tree$edge[, 1] == route_n[l, "parent"] &
                        tree$edge[, 2] == route_n[l, "child"])

        ## is there a change on this branch?
        if (bran %in% temp$edge){
          temp_sub <- subset(temp, temp$edge == bran)
          max_tran <- max(temp_sub$hmin)

          trans_trait[tree$tip.label[tip], "transition"] <-
            temp_sub$num[temp_sub$hmin ==  max_tran]

          trans_trait[tree$tip.label[tip], "state"] <-
            temp_sub$state[temp_sub$hmin ==  max_tran]
          break
        }
      }
    }

    # make any unchanged tips equal to root (0) and root state
    # identify rows with any NA
    rows_na <- apply(trans_trait, 1, function(row) any(is.na(row)))

    # fill those rows with values
    trans_trait[rows_na, ] <- matrix(
      rep(c(0, data$root.state[convT]), sum(rows_na)),
      nrow = sum(rows_na),
      byrow = TRUE
    )
    states <- unique(trans_trait[,"state"])
    for (s in 1:length(states)){
      by_state <- trans_trait[which(trans_trait[,"state"] == states[s])]
      if (length(unique(by_state)) > 1){
        n <- as.numeric(c(convT,states[s], length(unique(by_state))))
        convergent_traits <-rbind(convergent_traits, n)
      }
    }
  }

convergent_traits <- convergent_traits[-1,]
if (length(convergent_traits) > 0){
  if (length(convergent_traits["trait"])){
    convergent_traits <- as.data.frame(t(convergent_traits))
  } else {
    convergent_traits <- as.data.frame(convergent_traits)
  }
rownames(convergent_traits) <- seq(1,length(convergent_traits$trait), 1)
  }
  return(convergent_traits)
}


#' Determines the route (nodes and branches) for a tip in a phylogenetic tree
#' @description
#' This function traverses the tree to determine the evolutionary path (branches)
#' from root to a given tip
#'
#' @param tree phylogenetic tree
#' @param tip tip label
#'
#' @export
#'
#' @examples
#'route_n <- find_path_to_tip(tree, "t2")
#'
find_path_to_tip <- function(tree, tip) {
  # Ensure the tip is valid
  if (!tip %in% tree$tip.label) {
    stop("The specified tip is not found in the tree.")
  }

  # Get the edge matrix and root node
  edge <- tree$edge
  root <- ape::Ntip(tree) + 1

  # the node number of the tip
  tip_node <- which(tree$tip.label == tip)

  # Start from the tip node and move upward toward the root
  path <- tip_node
  current_node <- tip_node

  # loop until reaching the root node
  while (current_node != root) {
    # find the parent of the current node
    parent_node <- edge[edge[, 2] == current_node, 1]
    path <- c(parent_node, path)
    current_node <- parent_node
  }
  length(path)-1
  output <- matrix(nrow =length(path)-1, ncol=2)
  colnames(output) <- c("parent", "child")

  for (i in 1:length(path)-1){
    output[i,] <- c(path[i], path[i+1])
  }

  return(output)
}
