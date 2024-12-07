#' Plot full evolutionary history
#'
#'#' @description
#' This function creates a plot showing continuous evolution of discrete traits
#' @param data A morpho object
#' @param timetree TRUE or FALSE Indicate whether you want to plot a time tree or not. default FALSE
#' @param trait.num What trait number you want to visualize
#' @param br.rates Required if you provide a time tree
#' @import ape
#' @import FossilSim
#' @export

#' @examples

#' #set.seed(123)
#' # Non-time tree
#' phy <- ape::rtree(10)
#'
#' # simulate characters along the branches of the tree
#' continuous_traits <- sim.morpho.completeprocess(phy, k=4, trait.num =12)
#'
#' plot.morpho(data = continuous_traits, trait = 1, timetree = FALSE)
#'
#' # time tree
#' tree <- TreeSim::sim.bd.taxa(8,1,1,0.5)[[1]]
#' continuous_traits <- sim.morpho.completeprocess(time.tree = tree, br.rates = 0.2, k = 3, trait.num = 5)
#' plot.morpho(data = continuous_traits, trait = 1, timetree = TRUE, br.rates = 0.2 )


plot.morpho <- function(data = NULL, trait = NULL, timetree = FALSE, br.rates = NULL, col = c("#fdfdfd", "lightgray", "lightblue", "pink", "yellow", "green", "orange"), col.timescale = "darkgrey") {
  col_cont <- col
  ## Are we using a time tree?
  if (timetree) {
    plot(data$time.tree)
  } else {
    plot(data$tree)
  }

  if (class(trait) != "numeric") {
    print("Please select a viable integer for 'trait'")
  }

  if (trait >= length(data$sequences[[1]])) {
    print("Your selected character does not exist, please choose a lower integer")
  }


  # This is for sure not the right way to do it but will change asap!!
  tree_plot_info <- get("last_plot.phylo", envir = .PlotPhyloEnv)

  edge_start_x <- tree_plot_info$xx[data$tree$edge[, 1]]
  edge_end_x <- tree_plot_info$xx[data$tree$edge[, 2]]

  yy <- tree_plot_info$yy

  edge <- data$tree$edge
  parent <- as.integer(edge[, 1])
  child <- as.integer(edge[, 2])

  # Identify the root
  root <- as.integer(parent[!match(parent, child, 0)][1])


  ## Which trait?
  df <- data$continuous_traits[trait][[1]]

  if (nrow(df) > 0) {
    # Loop over the dataframe to add points
    for (i in 1:nrow(df)) {
      branch <- as.numeric(df$edge[i])

      if (timetree) {
        actual_position <- as.numeric(df$hmin[i]) * br.rates
        position <- (actual_position / (data$tree$edge.length[branch] * br.rates))
      } else {
        actual_position <- as.numeric(df$hmin[i])
        position <- (actual_position / data$tree$edge.length[branch])
      }

      # Calculate the point's coordinates along the branch
      point_x <- edge_start_x[branch] + position * (edge_end_x[branch] - edge_start_x[branch])

      if (timetree) {
        point_y <- yy[data$time.tree[["edge"]][branch, 2]]
      } else {
        point_y <- yy[data$tree[["edge"]][branch, 2]]
      }
      paint <- as.numeric(df$state[i]) + 1
      points(point_x, point_y, pch = 22, col = "black", bg = col_cont[paint], cex = 4)
      text(point_x, point_y, labels = as.numeric(df$state[i]))
    }
    points(0.05, yy[root], pch = 22, col = "black", bg = col_cont[as.numeric(data$root.states[trait])+1], cex = 4)
    text(0.05, yy[root], label = as.numeric(data$root.states[trait]))

  } else {
    text(0.05, yy[root], label = as.numeric(data$root.states[trait]))
    message("No transitions in this state across taxa")
  }

  # Add a timescale below the plot

if (timetree) {
  tree.age <- max(ape::node.depth.edgelength(data$time.tree))

  axis_labels <- c(tree.age,0)
  axis(1, at = c(0, tree.age), labels = round(axis_labels, 2), line = 1,col = col.timescale, lwd = 3, cex.axis = 1.3, col.axis = col.timescale)
  mtext("Time before present", side = 1, line = 2.5, cex = 1.3, col = col.timescale)
}
}
